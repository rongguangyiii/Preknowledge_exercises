
#include "solver/include/flowSolver.h"
#include "mesh/include/mesh.h"
#include "utility/include/log.h"
#include "iniCondition/include/GlobalData.h"
#include "utility/include/tools.h"
#include <cmath> 
#include <array>
#include <iostream>
#include <fstream>

void FlowSlover::solve()
{
	spdlog::info("Start solving...");
	Initializeflow();
	timeAdvance();
	spdlog::info("End solving...");
}

void FlowSlover::Initializeflow()
{
	std::vector<Mesh::nodePtr>& nodesVe = mesh_->getnodelist();
	for (auto& curnode : nodesVe)
	{
		auto& Privalue = curnode->getPrimitiveValue();
		if (Node::Nodetype::leftbound == curnode->getType())
		{
			Privalue[0] = 2.80;   //rho = 1.0
			Privalue[1] = 0.5;  //u   = 10.0
			Privalue[2] = 0.0;   //v   = 0.0
			Privalue[3] = 1.0;   //p   = 1.0
			//privalue[4] = 1.0;
		}
		else if (Node::Nodetype::leftbound != curnode->getType())
		{
			Privalue[0] = 1.4;
			Privalue[1] = 0.0;
			Privalue[2] = 0.0;
			Privalue[3] = 1.0;
			//privalue[4] = 1.0;
		}
		else
		{
			spdlog::error("Initializeflow error! unexpected node type!");
		}
		curnode->toConservedform();
	}

}

void FlowSlover::timeAdvance()
{
	size_t timestep = 0;
	writeTecplotFile(timestep); //输出第0步用于检查.
	//residual_ = 40.0;  // 初始化残差
	do
	{
		if (caltime_ > 30.0 || timestep > 500000)
			break;
		++timestep;
		//1.欧拉时间推进
		solverEuler();
		//2.当前步残差计算,并更新原始变量
		computeResidual();
		std::cout << "Timestep: " << timestep << " Residual: " << residual_ << std::endl;
		//spdlog::info("Timestep: {} Residual: {}", timestep, residual_);
		//3.输出当前时间步的结果
		writeTecplotFile(timestep);
	} while (residual_ > 1e-10);
}

void FlowSlover::solverEuler()
{
	residual_ = 0.0;
	double dt = 0.01;//dt后续再算
	std::vector<Mesh::nodePtr>& nodesVe = mesh_->getnodelist();
	for (auto& it : nodesVe)
	{
		//it->toPrimitiveform(); 
		auto thetype = it->getType();
		if (Node::Nodetype::inner != thetype)
			continue;
		std::vector<double> rhs = computeRHS(it);
		std::vector<double> U = it->getConservedValue();

		const size_t calsize = 4;
		std::vector<double> Unew(calsize, 0.0);
		double jacobian = it->getcoordTrans().jacob_;
		for (size_t ivar = 0; ivar < calsize; ivar++)
		{
			Unew[ivar] = U[ivar] - jacobian * rhs[ivar] * dt;
		}
		residual_ += std::pow(rhs[1], 2);//这只是一个量的残差，
		it->setConservedValue(Unew);
	}
	caltime_ += dt;
	//residual_ = std::sqrt(residual_);
	BoundExtrapolate();
}

std::vector<double> FlowSlover::computeRHS(std::shared_ptr<Node>& it)
{
	const size_t calsize = 4;
	const std::vector<std::shared_ptr<Node>>& neighbor = it->getNeighbor();
	std::vector<Mesh::nodePtr> nodeTemplateX = { neighbor[2],it,neighbor[0] };
	std::vector<Mesh::nodePtr> nodeTemplateY = { neighbor[3],it,neighbor[1] };
	std::vector<double> fluxplusX(calsize, 0.0);
	std::vector<double> fluxminusX(calsize, 0.0);
	split(nodeTemplateX, fluxplusX, fluxminusX);
	std::vector<double> fluxplusY(calsize, 0.0);
	std::vector<double> fluxminusY(calsize, 0.0);
	split(nodeTemplateY, fluxplusY, fluxminusY);

	std::vector<double> rhs(calsize, 0.0);
	//RHS = F(i+1/2) - F(i-1/2) + G(j+1/2)) - G(j-1/2)
	for (size_t iFlux = 0; iFlux < calsize; iFlux++)
	{
		rhs[iFlux] = fluxplusX[iFlux] - fluxminusX[iFlux] + fluxplusY[iFlux] - fluxminusY[iFlux];
	}

	return rhs;
}

void FlowSlover::BoundExtrapolate()
{
	std::vector<Mesh::nodePtr>& nodesVe = mesh_->getnodelist();
	for (auto& it : nodesVe)
	{
		auto thetype = it->getType();
		if (Node::Nodetype::inner == thetype || Node::Nodetype::leftbound == thetype)
			continue;
		const std::vector<std::shared_ptr<Node>>& neighbor = it->getNeighbor();
		std::vector<double> Unew;
		if (Node::Nodetype::downbound == thetype)
		{
			Unew = neighbor[1]->getConservedValue();
		}
		else if (Node::Nodetype::rightbound == thetype)
		{
			Unew = neighbor[1]->getConservedValue();
		}
		else if (Node::Nodetype::upbound == thetype)
		{
			Unew = neighbor[2]->getConservedValue();
		}
		else
		{
			spdlog::error("unexpected Extrapolate node type!");
		}
		it->setConservedValue(Unew);

	}
}

void FlowSlover::split(const std::vector<std::shared_ptr<Node>>& nodeTemplate, std::vector<double>& fluxplus, std::vector<double>& fluxminus)
{
	/* ----------------------------------------------------------------------------------------------------------------------------
		以ξ方向为例说明以下变量命名规则：
		halfmins : i-1/2        halfplus : i+1/2
		i-1/2半点处左值: halfMinsLeft : Q^-(i-1/2)   i-1/2半点处右值: halfMinsRight : Q^+(i-1/2)
		i+1/2半点处左值: halfPlusLeft : Q^-(i+1/2)   i+1/2半点处右值: halfPlusRight : Q^+(i+1/2)
		正通量: fluxPositive : F^+		负通量: fluxNegative : F^-

		则有：
		Q(i-1/2)处左值的正，负通量：fluxPositive_halfMinsLeft : F^+(Q^-(i-1/2))		fluxNegative_halfMinsLeft : F^-(Q^-(i-1/2))
		Q(i-1/2)处右值的正，负通量：fluxPositive_halfMinsRight : F^+(Q^+(i-1/2))	fluxNegative_halfMinsRight : F^-(Q^+(i-1/2))

		同理：
		Q(i+1/2)处右值的正，负通量：fluxPositive_halfPlusLeft : F^+(Q^-(i+1/2))		fluxNegative_halfPlusLeft : F^-(Q^-(i+1/2))
		Q(i+1/2)处右值的正，负通量：fluxPositive_halfPlusRight : F^+(Q^+(i+1/2))	fluxNegative_halfPlusRight : F^-(Q^+(i+1/2))

		nodeTemplate[0] : Q_left(i-1/2) ,nodeTemplate[1] : Q_right(i-1/2)
		nodeTemplate[1] : Q_left(i+1/2) ,nodeTemplate[3] : Q_right(i+1/2)
	 ---------------------------------------------------------------------------------------------------------------------------- */

	 //step 1. F(i-1/2)
	std::vector<double> fluxPositive_halfMinsLeft(4, 0.0);
	std::vector<double> fluxNegative_halfMinsLeft(4, 0.0);
	std::vector<double> fluxPositive_halfMinsRight(4, 0.0);
	std::vector<double> fluxNegative_halfMinsRight(4, 0.0);

	std::vector<double> rhs(4, 0.0);
	steger(nodeTemplate[0], nodeTemplate[1]->getcoordTrans(), fluxPositive_halfMinsLeft, fluxNegative_halfMinsLeft);//F_n(Q_left(i-1/2)),F_p(Q_left(i-1/2))
	steger(nodeTemplate[1], nodeTemplate[1]->getcoordTrans(), fluxPositive_halfMinsRight, fluxNegative_halfMinsRight);//F_n(Q_right(i-1/2)),F_p(Q_right(i-1/2))

	for (size_t iFlux = 0; iFlux < fluxminus.size(); iFlux++)
	{
		fluxminus[iFlux] = fluxPositive_halfMinsLeft[iFlux] + fluxNegative_halfMinsRight[iFlux];
	}

	//step 2. F(i+1/2)
	std::vector<double> fluxPositive_halfPlusLeft(4, 0.0);
	std::vector<double> fluxNegative_halfPlusLeft(4, 0.0);
	std::vector<double> fluxPositive_halfPlusRight(4, 0.0);
	std::vector<double> fluxNegative_halfPlusRight(4, 0.0);
	steger(nodeTemplate[1], nodeTemplate[1]->getcoordTrans(), fluxPositive_halfPlusLeft, fluxNegative_halfPlusLeft);//F_n(Q_left(i+1/2)),F_p(Q_left(i+1/2))
	steger(nodeTemplate[2], nodeTemplate[1]->getcoordTrans(), fluxPositive_halfPlusRight, fluxNegative_halfPlusRight);//F_n(Q_right(i+1/2)),F_p(Q_right(i+1/2))

	for (size_t iFlux = 0; iFlux < fluxplus.size(); iFlux++)
	{
		fluxplus[iFlux] = fluxPositive_halfPlusLeft[iFlux] + fluxNegative_halfPlusRight[iFlux];
	}


}

void FlowSlover::steger(const std::shared_ptr<Node>& nodeptr, const Node::TransCoef& trans, std::vector<double>& fp, std::vector<double>& fn)
{
	double ro = nodeptr->rho();
	double uu = nodeptr->u();
	double vv = nodeptr->v();
	double ee = nodeptr->E();
	double pp = nodeptr->p();
	//Node::TransCoef trans = nodeptr->getcoordTrans();
	double k_x = trans.ksi_x_;//冻结系数这里还要考虑，不能直接用邻居点的。
	double k_y = trans.ksi_y_;
	double k_t = trans.ksi_t_;
	double jj = trans.jacob_;
	double gama = GlobalData::GetDouble("refGamma");
	double epls = 1.0e-10;  // 熵修正函数中的较小正数

	//double pp = (gama - 1.0) * (ee - ro * (uu * uu + vv * vv) / 2.0);  // 压力
	double cc = std::sqrt(gama * pp / ro);  // 声速
	checkNaN(cc);

	double d_k = std::sqrt(k_x * k_x + k_y * k_y);  // 分母
	double ub = k_x * uu + k_y * vv + k_t;  // 分子坐标变化
	double uc = ub / d_k;  // 变化速度
	double k1 = k_x / d_k;
	double k2 = k_y / d_k;
	double k3 = k_t / d_k;

	// LAMBDA+
	double LAMBDA11 = 0.5 * (ub + std::sqrt(ub * ub + epls * epls));  // 特征值
	double LAMBDA12 = LAMBDA11;
	double LAMBDA13 = 0.5 * ((ub + cc * d_k) + std::sqrt((ub + cc * d_k) * (ub + cc * d_k) + epls * epls));
	double LAMBDA14 = 0.5 * ((ub - cc * d_k) + std::sqrt((ub - cc * d_k) * (ub - cc * d_k) + epls * epls));

	// LAMBDA-
	double LAMBDA21 = 0.5 * (ub - std::sqrt(ub * ub + epls * epls));
	double LAMBDA22 = LAMBDA21;
	double LAMBDA23 = 0.5 * ((ub + cc * d_k) - std::sqrt((ub + cc * d_k) * (ub + cc * d_k) + epls * epls));
	double LAMBDA24 = 0.5 * ((ub - cc * d_k) - std::sqrt((ub - cc * d_k) * (ub - cc * d_k) + epls * epls));

	// 正向质量通量
	double fmass12 = (LAMBDA13 + LAMBDA14 - 2.0 * LAMBDA11) * (0.5 * ro / gama);
	double fmass11 = ro * LAMBDA11 + fmass12;

	fp[0] = fmass11;
	fp[1] = fmass11 * uu + (LAMBDA13 - LAMBDA14) * (0.5 * k1 * ro * cc / gama);
	fp[2] = fmass11 * vv + (LAMBDA13 - LAMBDA14) * (0.5 * k2 * ro * cc / gama);
	fp[3] = fmass11 * ee + fmass12 * pp + (LAMBDA13 - LAMBDA14) * (0.5 * uc * ro * cc / gama);

	// LAMBDA没有除以jj，因此fp要除以jj
	for (int i = 0; i < 4; ++i) {
		fp[i] /= jj;
	}

	// 反向质量通量
	double fmass22 = (LAMBDA23 + LAMBDA24 - 2.0 * LAMBDA21) * (0.5 * ro / gama);
	double fmass21 = ro * LAMBDA21 + fmass22;

	fn[0] = fmass21;
	fn[1] = fmass21 * uu + (LAMBDA23 - LAMBDA24) * (0.5 * k1 * ro * cc / gama);
	fn[2] = fmass21 * vv + (LAMBDA23 - LAMBDA24) * (0.5 * k2 * ro * cc / gama);
	fn[3] = fmass21 * ee + fmass22 * pp + (LAMBDA23 - LAMBDA24) * (0.5 * uc * ro * cc / gama);

	// 同样fn要除以jj
	for (int i = 0; i < 4; ++i) {
		fn[i] /= jj;
	}

	checkNaN(fp);
	checkNaN(fn);
}

void FlowSlover::computeResidual()
{

	//residual_ = 0.0;
	std::vector<Mesh::nodePtr>& nodesVe = mesh_->getnodelist();
	for (auto& it : nodesVe)
	{
		it->toPrimitiveform();//检查顺序！
	}

	//for (int i = 1; i < nx; ++i) {
	//    for (int j = 1; j < ny; ++j) {
	//        residual += std::pow(u_new[i][j] - u[i][j], 2);
	//    }
	//}
	residual_ = std::sqrt(residual_);
}

void FlowSlover::vanLeer(const std::shared_ptr<Node>& nodeptr, const Node::TransCoef& trans, std::vector<double>& fp, std::vector<double>& fn)
{
	double ro = nodeptr->rho();
	double uu = nodeptr->u();
	double vv = nodeptr->v();
	double ee = nodeptr->E();
	//Node::TransCoef trans = nodeptr->getcoordTrans();
	double k_x = trans.ksi_x_;
	double k_y = trans.ksi_y_;
	double k_t = trans.ksi_t_;
	double jj = trans.jacob_;
	double gama = GlobalData::GetDouble("refGamma");
	//double epls = 1.0e-10;  // 熵修正函数中的较小正数

	double pp = (gama - 1.0) * (ee - ro * (uu * uu + vv * vv) / 2.0);  // 压力
	double a = std::sqrt(gama * pp / ro);  // 声速

	double d_k = std::sqrt(k_x * k_x + k_y * k_y);  // 分母
	double ub = k_x * uu + k_y * vv + k_t;  // 分子（坐标变化）
	double uc = ub / d_k;  // 当代速度
	double k1 = k_x / d_k;
	double k2 = k_y / d_k;
	double k3 = k_t / d_k;

	double mach = uc / a;
	if (mach >= 1.0) {  // 局部超声速
		fp[0] = (ro * uc) * (d_k / jj);  // positive direction
		fp[1] = (ro * uu * uc + pp * k1) * (d_k / jj);
		fp[2] = (ro * vv * uc + pp * k2) * (d_k / jj);
		fp[3] = ((ee + pp) * uc - pp * k3) * (d_k / jj);
		std::fill(fn.begin(), fn.end(), 0.0); // negative direction
		//fn.fill(0.0); 
	}
	else if (mach <= -1.0) {  // 局部超声速
		//fp.fill(0.0);
		std::fill(fp.begin(), fp.end(), 0.0);
		fn[0] = (ro * uc) * (d_k / jj);  // negative direction
		fn[1] = (ro * uu * uc + pp * k1) * (d_k / jj);
		fn[2] = (ro * vv * uc + pp * k2) * (d_k / jj);
		fn[3] = ((ee + pp) * uc - pp * k3) * (d_k / jj);
	}
	else {  // 局部亚声速
		double fmassp = (ro * a * std::pow(mach + 1.0, 2)) / 4.0;
		fp[0] = fmassp * (d_k / jj);
		fp[1] = (k1 * (-uc + 2.0 * a) / gama + uu) * fmassp * (d_k / jj);
		fp[2] = (k2 * (-uc + 2.0 * a) / gama + vv) * fmassp * (d_k / jj);
		fp[3] = (((1.0 - gama) * uc * uc + 2.0 * (gama - 1.0) * uc * a + 2.0 * a * a) / (gama * gama - 1.0) +
			(uu * uu + vv * vv) / 2.0 - k3 * (-uc + 2.0 * a) / gama) * fmassp * (d_k / jj);

		double fmassn = -ro * a * std::pow(mach - 1.0, 2) / 4.0;
		fn[0] = fmassn * (d_k / jj);
		fn[1] = (k1 * (-uc - 2.0 * a) / gama + uu) * fmassn * (d_k / jj);
		fn[2] = (k2 * (-uc - 2.0 * a) / gama + vv) * fmassn * (d_k / jj);
		fn[3] = (((1.0 - gama) * uc * uc - 2.0 * (gama - 1.0) * uc * a + 2.0 * a * a) / (gama * gama - 1.0) +
			(uu * uu + vv * vv) / 2.0 - k3 * (-uc - 2.0 * a) / gama) * fmassn * (d_k / jj);
	}
	checkNaN(fp);
	checkNaN(fn);
}

void FlowSlover::checkNaN(std::vector<double> value)
{
	for (const auto& per : value)
	{
		checkNaN(per);
	}
}

void FlowSlover::checkNaN(const double value)
{
	if (std::isnan(value))
	{
		spdlog::error("the result is NaN !!!");
		throw std::runtime_error("Errorthe result is NaN.");
	}
}

void FlowSlover::writeTecplotFile(const size_t timestep) const
{
	if (timestep % 100 != 0)
		return;
	std::string basePath = "result";
	Tools::createFolder(basePath);
	std::string filename = basePath + "/solution_" + std::to_string(timestep) + ".dat";
	std::ofstream os;
	os.open(filename);
	using std::endl;
	os << "TITLE =  MESH SOLUTION DATA " << std::endl;
	os << "VARIABLES = x, y, r, u, v, p" << std::endl;

	size_t outElementnum = mesh_->getelementnum();
	size_t outNodenum = mesh_->getnodenum();
	if (mesh_->getnodenum() == 0) {
		spdlog::error("the node num is 0!");
		return;
	}
	os << "ZONE t=\" MESH \"" << std::endl;
	os << "solutiontime = " << timestep << std::endl;
	os << "N=" << outNodenum << ", E=" << outElementnum << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
	for (size_t iNode = 0; iNode < outNodenum; ++iNode)
	{
		auto& currentNode = mesh_->getnode(iNode);
		os << currentNode->getCoord().x() << "  " << currentNode->getCoord().y() << "  "
			<< currentNode->rho() << "  " << currentNode->u() << "  "
			<< currentNode->v() << "  " << currentNode->p() << "  " << endl;
	}

	for (size_t iElem = 0; iElem < outElementnum; ++iElem)
	{
		auto& currentElement = mesh_->getelement(iElem);
		vector<size_t>updateNodeIndex;
		const auto& elementNode = currentElement->getnodelist();
		for (size_t iNode = 0; iNode < elementNode.size(); ++iNode)
		{
			auto& currenNode = elementNode[iNode];
			os << currenNode->getIndex() + 1 << "  ";
		}
		if (elementNode.size() == 3)
			os << elementNode[2]->getIndex() + 1;
		os << endl;
	}
	os.close();
	return;

}