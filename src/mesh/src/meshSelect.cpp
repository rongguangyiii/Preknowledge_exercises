
#include"mesh/include/meshsSelect.h"
#include "iniCondition/include/GlobalData.h"
#include "utility/include/tools.h"

MeshsSelect::MeshsSelect(std::shared_ptr<Mesh>& mesh)
	: mesh_(mesh),
	x(GlobalData::GetInt("xNodeNum"), std::vector<double>(GlobalData::GetInt("yNodeNum"), 0.0)),
	y(GlobalData::GetInt("xNodeNum"), std::vector<double>(GlobalData::GetInt("yNodeNum"), 0.0))
{
	// 这里你可以执行其他初始化代码
}

void MeshsSelect::coordTrans2ord()
{
	double XK, XI, YK, YI, JJ;
	size_t mx = GlobalData::GetInt("xNodeNum");
	size_t my = GlobalData::GetInt("yNodeNum");
	for (size_t i = 0; i < mx; ++i) {
		for (size_t j = 0; j < my; ++j) {
			if (i == 0 && j != 0 && j != my - 1) {  // 1. 左边界 (i == 0) 非角点 (0 < j < my)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}
			else if (i == mx - 1 && j != 0 && j != my - 1) {  // 2. 右边界 (i == mx) 非角点 (0 < j < my)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}
			else if (i == 0 && j == 0) {  // 3. 左下角 (i == 0, j == 0)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
			}
			else if (i == 0 && j == my - 1) {  // 4. 左上角 (i == 0, j == my)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
			}
			else if (i == mx - 1 && j == my - 1) {  // 5. 右上角 (i == mx, j == my)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
			}
			else if (i == mx - 1 && j == 0) {  // 6. 右下角 (i == mx, j == 0)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
			}
			else if (j == 0) {  // 7. 下边界 (j == 0) 非角点 (0 < i < mx)
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
			}
			else if (j == my - 1) {  // 8. 上边界 (j == my) 非角点 (0 < i < mx)
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
			}
			else {  // 9. 内部点
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}

			JJ = XK * YI - XI * YK;
			double ksi_x = YI / JJ;
			double ksi_y = -XI / JJ;
			double eta_x = -YK / JJ;
			double eta_y = XK / JJ;
			double jacobian = 1.0 / JJ;

			const size_t perNodeindex = (j * mx) + i;
			auto& curnode = mesh_->getnode(perNodeindex);
			curnode->setCoordTrans(ksi_x, ksi_y, eta_x, eta_y, jacobian);

		}
	}
}

void MeshsSelect::genNeiborNode()
{
	size_t xNodeNum = GlobalData::GetInt("xNodeNum");
	size_t yNodeNum = GlobalData::GetInt("yNodeNum");
	size_t i, j;
	auto& nodelist = mesh_->getnodelist();
	for (size_t iNode = 0; iNode < mesh_->getnodelist().size(); ++iNode)
	{
		auto& currentNode = mesh_->getnode(iNode);
		Node::Nodetype currentNodeType = Node::Nodetype::inner;
		std::vector<std::shared_ptr<Node>>& neiborNode = currentNode->getNeighbor();

		i = iNode % xNodeNum;
		j = iNode / xNodeNum;
		std::vector<size_t> neiborNodeIndex = { iNode + 1, iNode + xNodeNum,iNode - 1 , iNode - xNodeNum };//设置邻居节点的索引，右，上，左，下
		if (j == 0)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin() + 3);
			currentNodeType = Node::Nodetype::downbound;
		}

		if (j == yNodeNum - 1)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin() + 1);
			currentNodeType = Node::Nodetype::upbound;
		}

		if (i == 0)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin() + 2);
			currentNodeType = Node::Nodetype::leftbound;
		}

		if (i == xNodeNum - 1)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin());
			currentNodeType = Node::Nodetype::rightbound;
		}

		//if ((i == xNodeNum - 1 || i == 0) && (j == 0 || j == yNodeNum - 1))
		//{
		//	currentNodeType = Node::Nodetype::corner;
		//}
		currentNode->setType(currentNodeType);

		//neiborNode.resize(neiborNodeIndex.size());
		for (const auto it:neiborNodeIndex)
			neiborNode.emplace_back(mesh_->getnode(it));
	}
}

void MeshsSelect::genmesh()
{
	gennode();
	genNeiborNode();
	coordTrans2ord();
	genelement();
}


void UniformMesh::gennode()
{
	const size_t mx = GlobalData::GetInt("xNodeNum");
	const size_t my = GlobalData::GetInt("yNodeNum");
	const double xa = GlobalData::GetDouble("xMin");
	const double ya = GlobalData::GetDouble("yMin");
	const double xb = GlobalData::GetDouble("xMax");
	const double yb = GlobalData::GetDouble("yMax");
	const double dx = (xb - xa) / (mx - 1);
	const double dy = (yb - ya) / (my - 1);
	auto& meshNodevec = mesh_->getnodelist();
	size_t count = 0;
	for (size_t j = 0; j < my; j++)
	{
		for (size_t i = 0; i < mx; i++)
		{
			const double x = xa + i * dx;
			const double y = ya + j * dy;
			this->x[i][j] = x;
			this->y[i][j] = y;
			if (i == 0)
			{
				meshNodevec.emplace_back(std::make_shared<Node>(Coordinate{ x, y }, count++, mesh_, Node::Nodetype::leftbound));
			}
			else
			{
				meshNodevec.emplace_back(std::make_shared<Node>(Coordinate{ x, y }, count++, mesh_));
			}
		}
	}

}

void UniformMesh::genelement()
{
	auto& meshElevec = mesh_->getelementlist();
	const size_t mx = GlobalData::GetInt("xNodeNum");
	const size_t my = GlobalData::GetInt("yNodeNum");
	size_t eleNum = 0;
	for (size_t j = 0; j < my - 1; j++)
	{
		//spdlog::debug("第 {} 行的单元",j);
		for (size_t i = 0; i < mx - 1; i++)
		{
			const size_t perNodeindex = (j * mx) + i;
			/*spdlog::debug("当前单元的四个节点编号：{},{},{},{}", perNodeindex, perNodeindex + 1,
				perNodeindex + mx + 1, perNodeindex + mx);*/
			auto& node_1 = mesh_->getnode(perNodeindex);
			auto& node_2 = mesh_->getnode(perNodeindex + 1);
			auto& node_3 = mesh_->getnode(perNodeindex + mx + 1);
			auto& node_4 = mesh_->getnode(perNodeindex + mx);
			//std::vector<std::shared_ptr<Node>> perNodevec = { node_1, node_2, node_3, node_4 };
			meshElevec.emplace_back(std::make_shared<Element>(eleNum++, node_1, node_2, node_3, node_4));

		}
	}
}

void SemicircularMesh::gennode()
{
	const double R0 = 1;
	const double R1 = 5;
	size_t mx = GlobalData::GetInt("xNodeNum");
	size_t my = GlobalData::GetInt("yNodeNum");
	auto& meshNodevec = mesh_->getnodelist();
	size_t count = 0;
	for (size_t j = 0; j < my; j++)
	{
		for (size_t i = 0; i < mx; i++)
		{
			double R = R1 - j * (R1 - R0) / (my - 1);
			double theta = Tools::PI * i / (mx - 1);
			double x = -R * sin(theta);
			double y = R * cos(theta);
			if (j == 0)
			{
				meshNodevec.emplace_back(std::make_shared<Node>(Coordinate{ x, y }, count++, mesh_, Node::Nodetype::leftbound));
			}
			else
			{
				meshNodevec.emplace_back(std::make_shared<Node>(Coordinate{ x, y }, count++, mesh_));
			}
		}
	}

}

void SemicircularMesh::genelement()
{
	auto& meshElevec = mesh_->getelementlist();
	const size_t mx = GlobalData::GetInt("xNodeNum");
	const size_t my = GlobalData::GetInt("yNodeNum");
	size_t eleNum = 0;
	for (size_t j = 0; j < my - 1; j++)
	{
		//spdlog::debug("第 {} 行的单元",j);
		for (size_t i = 0; i < mx - 1; i++)
		{
			const size_t perNodeindex = (j * mx) + i;
			/*spdlog::debug("当前单元的四个节点编号：{},{},{},{}", perNodeindex, perNodeindex + 1,
				perNodeindex + mx + 1, perNodeindex + mx);*/
			auto& node_1 = mesh_->getnode(perNodeindex);
			auto& node_2 = mesh_->getnode(perNodeindex + 1);
			auto& node_3 = mesh_->getnode(perNodeindex + mx + 1);
			auto& node_4 = mesh_->getnode(perNodeindex + mx);
			//std::vector<std::shared_ptr<Node>> perNodevec = { node_1, node_2, node_3, node_4 };
			meshElevec.emplace_back(std::make_shared<Element>(eleNum++, node_1, node_2, node_3, node_4));

		}
	}
}