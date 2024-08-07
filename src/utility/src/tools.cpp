
#include "utility/include/tools.h"
#include "utility/include/log.h"
#include <vector>
#include <filesystem>  
#include <fstream>

const Coordinate Tools::Cross(const Coordinate& vA, const Coordinate& vB) {
	return Coordinate(
		vA.y() * vB.z() - vA.z() * vB.y(), // i分量  
		vA.z() * vB.x() - vA.x() * vB.z(), // j分量  
		vA.x() * vB.y() - vA.y() * vB.x()  // k分量  
	);
}

/*-------------------------------------------------------------------
*  Function: isPointInPolygon_moreInner
*  Purpose: 判断点是否在多边形内部版本1：包括与边界重合或顶点重合的点都认为在内部
*  Arguments:
*    p - 目标判断点
*    polygon - 有序封闭多边形
*  Returns:
*    boll
-------------------------------------------------------------------*/
bool Tools::isPointInPolygon_moreInner(const Coordinate& p, const std::vector<Coordinate>& polygon) {
	const double EPSILON = 1e-10;
	bool inside = false;
	size_t j = polygon.size() - 1;

	for (size_t i = 0; i < polygon.size(); i++) {
		// Check if point is on vertex
		if ((std::fabs(p.x() - polygon[i].x()) < EPSILON) && (std::fabs(p.y() - polygon[i].y()) < EPSILON)) {
			return true; // Point is on a vertex
		}

		// Check if point is on edge
		double minX = std::min(polygon[i].x(), polygon[j].x());
		double maxX = std::max(polygon[i].x(), polygon[j].x());
		double minY = std::min(polygon[i].y(), polygon[j].y());
		double maxY = std::max(polygon[i].y(), polygon[j].y());

		if ((minX <= p.x() && p.x() <= maxX) && (minY <= p.y() && p.y() <= maxY)) {
			double dx = polygon[j].x() - polygon[i].x();
			double dy = polygon[j].y() - polygon[i].y();
			if (std::fabs(dy) < EPSILON) { // horizontal line
				if (std::fabs(p.y() - polygon[i].y()) < EPSILON) {
					return true; // Point is on horizontal edge
				}
			}
			else if (std::fabs(dx) < EPSILON) { // vertical line
				if (std::fabs(p.x() - polygon[i].x()) < EPSILON) {
					return true; // Point is on vertical edge
				}
			}
			else {
				double slope = dy / dx;
				double intercept = polygon[i].y() - slope * polygon[i].x();
				if (std::fabs(p.y() - (slope * p.x() + intercept)) < EPSILON) {
					return true; // Point is on non-vertical edge
				}
			}
		}

		// Check if point is inside using ray-casting algorithm
		if (((polygon[i].y() <= p.y()) && (p.y() < polygon[j].y())) ||
			((polygon[j].y() <= p.y()) && (p.y() < polygon[i].y()))) {
			double intersectionX = polygon[i].x() + (p.y() - polygon[i].y()) * (polygon[j].x() - polygon[i].x()) / (polygon[j].y() - polygon[i].y());
			if (intersectionX < p.x()) {
				inside = !inside;
			}
		}

		j = i;
	}
	return inside;
}

/*-------------------------------------------------------------------
*  Function: isPointInPolygon_moreOuter
*  Purpose: 判断点是否在多边形内部版本2：包括与边界重合或顶点重合的点都认为在外部
*  Arguments:
*    p - 目标判断点
*    polygon - 有序封闭多边形
*  Returns:
*    boll
-------------------------------------------------------------------*/
bool Tools::isPointInPolygon_moreOuter(const Coordinate& p, const std::vector<Coordinate>& polygon) {
	const double EPSILON = 1e-10;
	bool inside = false;
	size_t j = polygon.size() - 1;

	for (size_t i = 0; i < polygon.size(); i++) {
		// Check if point is on vertex
		if ((std::fabs(p.x() - polygon[i].x()) < EPSILON) && (std::fabs(p.y() - polygon[i].y()) < EPSILON)) {
			return false; // Point is on a vertex
		}

		// Check if point is on edge
		double minX = std::min(polygon[i].x(), polygon[j].x());
		double maxX = std::max(polygon[i].x(), polygon[j].x());
		double minY = std::min(polygon[i].y(), polygon[j].y());
		double maxY = std::max(polygon[i].y(), polygon[j].y());

		if ((minX <= p.x() && p.x() <= maxX) && (minY <= p.y() && p.y() <= maxY)) {
			double dx = polygon[j].x() - polygon[i].x();
			double dy = polygon[j].y() - polygon[i].y();
			if (std::fabs(dy) < EPSILON) { // horizontal line
				if (std::fabs(p.y() - polygon[i].y()) < EPSILON) {
					return false; // Point is on horizontal edge
				}
			}
			else if (std::fabs(dx) < EPSILON) { // vertical line
				if (std::fabs(p.x() - polygon[i].x()) < EPSILON) {
					return false; // Point is on vertical edge
				}
			}
			else {
				double slope = dy / dx;
				double intercept = polygon[i].y() - slope * polygon[i].x();
				if (std::fabs(p.y() - (slope * p.x() + intercept)) < EPSILON) {
					return false; // Point is on non-vertical edge
				}
			}
		}

		// Check if point is inside using ray-casting algorithm
		if (((polygon[i].y() <= p.y()) && (p.y() < polygon[j].y())) ||
			((polygon[j].y() <= p.y()) && (p.y() < polygon[i].y()))) {
			double intersectionX = polygon[i].x() + (p.y() - polygon[i].y()) * (polygon[j].x() - polygon[i].x()) / (polygon[j].y() - polygon[i].y());
			if (intersectionX < p.x()) {
				inside = !inside;
			}
		}

		j = i;
	}
	return inside;
}

//是否在两点之间
bool Tools::IsInTwoPoint(const Coordinate& P1, const Coordinate& P2, const Coordinate& judgeP)
{
	double coor_x_max = 0, coor_x_min = 0;
	if (P1.x() > P2.x())
	{
		coor_x_max = P1.x();
		coor_x_min = P2.x();
	}
	else
	{
		coor_x_max = P2.x();
		coor_x_min = P1.x();
	}
	double coor_y_max = 0, coor_y_min = 0;
	if (P1.y() > P2.y())
	{
		coor_y_max = P1.y();
		coor_y_min = P2.y();
	}
	else
	{
		coor_y_max = P2.y();
		coor_y_min = P1.y();
	}
	if ((judgeP.x() < coor_x_max && judgeP.x() > coor_x_min) || (judgeP.y() < coor_y_max && judgeP.y() > coor_y_min))
	{
		return true;
	}
	else if (abs(judgeP.x() - coor_x_max) < 1E13 || abs(judgeP.x() - coor_x_min) < 1E13 ||
		abs(judgeP.y() - coor_y_max) < 1E13 || abs(judgeP.y() - coor_y_min) < 1E13)
	{
		//坐标相等的情况
		return true;
	}
	else
	{
		return false;
	}
}

//点到直线的距离
double Tools::distanceToLine(const Coordinate& A, const Coordinate& B, const Coordinate& P)
{
	double numerator = std::abs((B.y() - A.y()) * P.x() - (B.x() - A.x()) * P.y() + B.x() * A.y() - A.x() * B.y());
	double denominator = std::sqrt((B.y() - A.y()) * (B.y() - A.y()) + (B.x() - A.x()) * (B.x() - A.x()));
	return numerator / denominator;
}

//该函数判断两个向量是否同向
bool Tools::IsSameDirection(const Coordinate& vec1, const Coordinate& vec2)
{
	double dot = vec1.x() * vec2.x() + vec1.y() * vec2.y() + vec1.z() * vec2.z();
	if (dot > 0)
		return true;
	else
		return false;
}

//单位化向量
Coordinate Tools::unitNormalize(const Coordinate& v) {
	double len = sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
	// 确定向量不是零向量
	if (len > 0) {
		return Coordinate{ v.x() / len, v.y() / len, v.z() / len };
	}
	else {
		// 如果是零向量，返回一个无效的向量或者抛出异常
		return Coordinate{ 0, 0, 0 };
	}
}



/*-------------------------------------------------------------------
*  Function: createFolder()
*  Purpose: 在指定路径创建文件夹，存在则不创建
*  Usage: createFolder("tempData")
*  Arguments:
*    basePath - 文件夹路径
*  Returns:
*    null
-------------------------------------------------------------------*/
void Tools::createFolder(const std::string& basePath)
{
	// 构造完整的文件夹路径  
	std::filesystem::path folderPath = basePath;

	// 检查文件夹是否已经存在  
	if (!std::filesystem::exists(folderPath)) {
		// 创建文件夹  
		std::filesystem::create_directories(folderPath);
		spdlog::info("Folder '{}' created successfully!", basePath);
	}
}


//在二维情况下，给定三个Coordinate，判断这三个点是否共圆
bool Tools::isCocircular(const Coordinate& p1, const Coordinate& p2, const Coordinate& p3) {
	// 使用行列式判断三个点是否共圆  
	// 计算行列式的值：如果三个点共圆，那么这个行列式的值理论上应为零
	double det =
		(p1.x() * p1.x() + p1.y() * p1.y()) * (p2.y() - p3.y()) * (p1.x() - p2.x()) +
		(p2.x() * p2.x() + p2.y() * p2.y()) * (p3.y() - p1.y()) * (p2.x() - p3.x()) +
		(p3.x() * p3.x() + p3.y() * p3.y()) * (p1.y() - p2.y()) * (p3.x() - p1.x());

	// 使用相对误差判断行列式值是否足够小，这里取一个相对比例例如1e-9作为阈值
	// 注意：实际应用中，阈值的选择可能需要根据具体情况调整
	double threshold = std::fabs(det) / (1 + std::fabs(p1.x() * p1.x() + p1.y() * p1.y()) +
		std::fabs(p2.x() * p2.x() + p2.y() * p2.y()) +
		std::fabs(p3.x() * p3.x() + p3.y() * p3.y()));
	return threshold < 1e-8;
}

//在二维情况下，给定三个不共线的Coordinate，计算这三个点共圆的圆心
Coordinate Tools::getCocircularCenter(const Coordinate& p1, const Coordinate& p2, const Coordinate& p3)
{
	// 向量定义
	Coordinate v1 = p2 - p1;
	Coordinate v2 = p3 - p2;
	Coordinate v3 = p1 - p3;

	// 计算向量的叉乘（外积）模长的两倍，作为三角形面积的2倍
	double areaTimes2 = v1.x() * v2.y() - v1.y() * v2.x();

	// 防止除以零的情况（三点共线），此时areaTimes2应为0
	if (std::fabs(areaTimes2) < 1e-12) {
		throw std::invalid_argument("The three points are collinear and cannot define a unique circle.");
	}

	// 计算每个点到另外两点连线的向量叉乘模长的和，再除以2*areaTimes2得到圆心坐标
	double x = ((v1.norm() * v1.norm() * (p2.y() - p3.y()) + v2.norm() * v2.norm() * (p3.y() - p1.y()) + v3.norm() * v3.norm() * (p1.y() - p2.y())) / (2 * areaTimes2));
	double y = ((v1.norm() * v1.norm() * (p3.x() - p2.x()) + v2.norm() * v2.norm() * (p1.x() - p3.x()) + v3.norm() * v3.norm() * (p2.x() - p1.x())) / (2 * areaTimes2));

	return Coordinate(x, y);
}

/*-------------------------------------------------------------------
*  Function: areCollinear(Coordinate p1, Coordinate p2, Coordinate p3)
*  Purpose: 检查三角形的三个点是否共线
*  Arguments:
*   Coordinate p1,p2,p3 - 三角形的三个顶点
*  Returns:
*    bool - 是否共线的布尔值
-------------------------------------------------------------------*/
bool Tools::isCollinear(const Coordinate& p1, const Coordinate& p2, const Coordinate& p3)
{
	// 使用斜率法判断三点是否共线  
	// 如果斜率相同或任意两点重合（在容差范围内），则它们共线  
	constexpr double epsilon = std::numeric_limits<double>::epsilon();
	double dx1 = p2.x() - p1.x();
	double dy1 = p2.y() - p1.y();
	double dx2 = p3.x() - p1.x();
	double dy2 = p3.y() - p1.y();

	// 使用斜率来判断是否共线  
	if (std::abs(dx1 * dy2 - dx2 * dy1) < epsilon) {
		return true;
	}
	// 检查是否接近0（即检查是否有重合的点）  
	if (std::abs(dx1) < epsilon && std::abs(dx2) < epsilon) {
		// 如果x坐标都接近相同，则检查y坐标是否也接近相同  
		return std::abs(dy1) < epsilon && std::abs(dy2) < epsilon;
	}

	return false;
}

//将角度标准化到0到π之间。
double Tools::normalizeAngle(double angleRad) 
{
	angleRad = fmod(angleRad, 2 * PI);
	if (angleRad < 0) {
		angleRad += 2 * PI;
	}
	if (angleRad > PI) {
		angleRad = 2 * PI - angleRad;
	}
	return angleRad;
}

//根据angle旋转坐标A
Coordinate Tools::rotateNormal(const Coordinate& A, const double angle) {
	double newX = A.x() * cos(angle) - A.y() * sin(angle);
	double newY = A.x() * sin(angle) + A.y() * cos(angle);
	return Coordinate(newX, newY);
}

/*-------------------------------------------------------------------------------
*  Function: genLinePoints()
*  Purpose: 在给定的两个点中间添加点
*  Usage:
*	 genLinePoints(coor1, coor2, 0.0);
*	 genLinePoints(coor1, coor2, 0.8);
*  Arguments:
*    pt1 - 起始点
*    pt2 - 结束点
*    delta - 控制长度：当长度为0.0时，即不添加点。
*  Returns:
*    null
-------------------------------------------------------------------------------*/
std::vector<Coordinate> Tools::genLinePoints(const Coordinate& pt1, const Coordinate& pt2, const double delta)
{
	const double EPSILON = 1e-12;//用于避免由于浮点数的表示误差导致的比较问题,出现在两个浮点数相等的情况，尽管这样的情况几乎不会发生。
	if (delta < EPSILON)
		return { };

	std::vector<Coordinate> tempPointVec;
	double dis = distance(pt1, pt2);
	double dx = (pt2.x() - pt1.x()) * delta / dis;
	double dy = (pt2.y() - pt1.y()) * delta / dis;
	tempPointVec.push_back(pt1);
	//当剩下的线段长度大于等于delta的1.5倍，就继续生成点。这样做是为了避免最后一个点距离pt2太近。
	while (distance(tempPointVec[tempPointVec.size() - 1], pt2) > 1.5 * delta + EPSILON)
		tempPointVec.push_back({ tempPointVec[tempPointVec.size() - 1].x() + dx,tempPointVec[tempPointVec.size() - 1].y() + dy, 0, CoorType::patch });
	tempPointVec.push_back(pt2);
	return tempPointVec;
}

//返回三个有序点构成的夹角，前提设定B点一定在A点和C点之间
const double Tools::angle_ABC(const Coordinate& ptA, const Coordinate& ptB, const Coordinate& ptC)
{
	auto V1 = ptA - ptB;
	auto V2 = ptC - ptB;
	double angleRad_V1 = atan2(V1.y(), V1.x());
	double angleRad_V2 = atan2(V2.y(), V2.x());
	double gap_V1_V2 = Tools::normalizeAngle(angleRad_V2 - angleRad_V1) * 180.0 / PI;
	return gap_V1_V2;
}

Coordinate Tools::returnPolycenter(const std::vector<Coordinate>& points)
{
	// 计算凸多边形的中心点
	Coordinate center;
	for (const auto& point : points)
		center += point;
	center /= points.size();
	return center;
}

void Tools::OutCompleteMesh2Tecplot(const std::shared_ptr<Mesh>& meshin)
{
	std::string basePath = "tempData";
	Tools::createFolder(basePath);
	std::ofstream os;
	os.open("tempData/Mesh_1.dat");
	using std::endl;
	os << "TITLE = MESH" << std::endl;
	os << "variables=x,y" << std::endl;

	size_t outElementnum = meshin->getelementnum();
	size_t outNodenum = meshin->getnodenum();
	if (meshin->getnodenum() == 0) {
		spdlog::error("the node num is 0!");
		return;
	}
	os << "ZONE t=\" MESH \"" << std::endl;
	os << "DATAPACKING = BLOCK" << std::endl;
	os << "N=" << outNodenum << ", E=" << outElementnum << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
	for (size_t iNode = 0; iNode < outNodenum; ++iNode)
	{
		auto& currentNode = meshin->getnode(iNode);
		os << currentNode->getCoord().x() << "  " << currentNode->getCoord().y() << "  " << endl;
	}

	for (size_t iElem = 0; iElem < outElementnum; ++iElem)
	{
		auto& currentElement = meshin->getelement(iElem);
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


// 定义一个函数来计算点相对于中心点的极角
double Tools::calculateAngle(const Coordinate& point, const Coordinate& center) {
	double angle = atan2(point.y() - center.y(), point.x() - center.x());
	// 将角度调整到0到2pi之间
	if (angle < 0) 
		angle += 2 * PI;
	return angle;
}

bool Tools::isVectorCollinear(const Coordinate& a, const Coordinate& b)
{
	// 归一化向量:a,b
	double epsilon = 1e-8;

	// 计算点积
	double dot_product = a.x() * b.x() + a.y() * b.y() + a.z() * b.z();

	// 检查点积是否接近1或-1
	return std::fabs(dot_product - 1.0) < epsilon || std::fabs(dot_product + 1.0) < epsilon;
}

double Tools::distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

double Tools::distance(const Coordinate& pt1, const Coordinate& pt2)
{
	return distance(pt1.x(), pt1.y(), pt1.z(), pt2.x(), pt2.y(), pt2.z());
}