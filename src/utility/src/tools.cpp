
#include "utility/include/tools.h"
#include "utility/include/log.h"
#include <vector>
#include <filesystem>  
#include <fstream>

const Coordinate Tools::Cross(const Coordinate& vA, const Coordinate& vB) {
	return Coordinate(
		vA.y() * vB.z() - vA.z() * vB.y(), // i����  
		vA.z() * vB.x() - vA.x() * vB.z(), // j����  
		vA.x() * vB.y() - vA.y() * vB.x()  // k����  
	);
}

/*-------------------------------------------------------------------
*  Function: isPointInPolygon_moreInner
*  Purpose: �жϵ��Ƿ��ڶ�����ڲ��汾1��������߽��غϻ򶥵��غϵĵ㶼��Ϊ���ڲ�
*  Arguments:
*    p - Ŀ���жϵ�
*    polygon - �����ն����
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
*  Purpose: �жϵ��Ƿ��ڶ�����ڲ��汾2��������߽��غϻ򶥵��غϵĵ㶼��Ϊ���ⲿ
*  Arguments:
*    p - Ŀ���жϵ�
*    polygon - �����ն����
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

//�Ƿ�������֮��
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
		//������ȵ����
		return true;
	}
	else
	{
		return false;
	}
}

//�㵽ֱ�ߵľ���
double Tools::distanceToLine(const Coordinate& A, const Coordinate& B, const Coordinate& P)
{
	double numerator = std::abs((B.y() - A.y()) * P.x() - (B.x() - A.x()) * P.y() + B.x() * A.y() - A.x() * B.y());
	double denominator = std::sqrt((B.y() - A.y()) * (B.y() - A.y()) + (B.x() - A.x()) * (B.x() - A.x()));
	return numerator / denominator;
}

//�ú����ж����������Ƿ�ͬ��
bool Tools::IsSameDirection(const Coordinate& vec1, const Coordinate& vec2)
{
	double dot = vec1.x() * vec2.x() + vec1.y() * vec2.y() + vec1.z() * vec2.z();
	if (dot > 0)
		return true;
	else
		return false;
}

//��λ������
Coordinate Tools::unitNormalize(const Coordinate& v) {
	double len = sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
	// ȷ����������������
	if (len > 0) {
		return Coordinate{ v.x() / len, v.y() / len, v.z() / len };
	}
	else {
		// �����������������һ����Ч�����������׳��쳣
		return Coordinate{ 0, 0, 0 };
	}
}



/*-------------------------------------------------------------------
*  Function: createFolder()
*  Purpose: ��ָ��·�������ļ��У������򲻴���
*  Usage: createFolder("tempData")
*  Arguments:
*    basePath - �ļ���·��
*  Returns:
*    null
-------------------------------------------------------------------*/
void Tools::createFolder(const std::string& basePath)
{
	// �����������ļ���·��  
	std::filesystem::path folderPath = basePath;

	// ����ļ����Ƿ��Ѿ�����  
	if (!std::filesystem::exists(folderPath)) {
		// �����ļ���  
		std::filesystem::create_directories(folderPath);
		spdlog::info("Folder '{}' created successfully!", basePath);
	}
}


//�ڶ�ά����£���������Coordinate���ж����������Ƿ�Բ
bool Tools::isCocircular(const Coordinate& p1, const Coordinate& p2, const Coordinate& p3) {
	// ʹ������ʽ�ж��������Ƿ�Բ  
	// ��������ʽ��ֵ����������㹲Բ����ô�������ʽ��ֵ������ӦΪ��
	double det =
		(p1.x() * p1.x() + p1.y() * p1.y()) * (p2.y() - p3.y()) * (p1.x() - p2.x()) +
		(p2.x() * p2.x() + p2.y() * p2.y()) * (p3.y() - p1.y()) * (p2.x() - p3.x()) +
		(p3.x() * p3.x() + p3.y() * p3.y()) * (p1.y() - p2.y()) * (p3.x() - p1.x());

	// ʹ���������ж�����ʽֵ�Ƿ��㹻С������ȡһ����Ա�������1e-9��Ϊ��ֵ
	// ע�⣺ʵ��Ӧ���У���ֵ��ѡ�������Ҫ���ݾ����������
	double threshold = std::fabs(det) / (1 + std::fabs(p1.x() * p1.x() + p1.y() * p1.y()) +
		std::fabs(p2.x() * p2.x() + p2.y() * p2.y()) +
		std::fabs(p3.x() * p3.x() + p3.y() * p3.y()));
	return threshold < 1e-8;
}

//�ڶ�ά����£��������������ߵ�Coordinate�������������㹲Բ��Բ��
Coordinate Tools::getCocircularCenter(const Coordinate& p1, const Coordinate& p2, const Coordinate& p3)
{
	// ��������
	Coordinate v1 = p2 - p1;
	Coordinate v2 = p3 - p2;
	Coordinate v3 = p1 - p3;

	// ���������Ĳ�ˣ������ģ������������Ϊ�����������2��
	double areaTimes2 = v1.x() * v2.y() - v1.y() * v2.x();

	// ��ֹ���������������㹲�ߣ�����ʱareaTimes2ӦΪ0
	if (std::fabs(areaTimes2) < 1e-12) {
		throw std::invalid_argument("The three points are collinear and cannot define a unique circle.");
	}

	// ����ÿ���㵽�����������ߵ��������ģ���ĺͣ��ٳ���2*areaTimes2�õ�Բ������
	double x = ((v1.norm() * v1.norm() * (p2.y() - p3.y()) + v2.norm() * v2.norm() * (p3.y() - p1.y()) + v3.norm() * v3.norm() * (p1.y() - p2.y())) / (2 * areaTimes2));
	double y = ((v1.norm() * v1.norm() * (p3.x() - p2.x()) + v2.norm() * v2.norm() * (p1.x() - p3.x()) + v3.norm() * v3.norm() * (p2.x() - p1.x())) / (2 * areaTimes2));

	return Coordinate(x, y);
}

/*-------------------------------------------------------------------
*  Function: areCollinear(Coordinate p1, Coordinate p2, Coordinate p3)
*  Purpose: ��������ε��������Ƿ���
*  Arguments:
*   Coordinate p1,p2,p3 - �����ε���������
*  Returns:
*    bool - �Ƿ��ߵĲ���ֵ
-------------------------------------------------------------------*/
bool Tools::isCollinear(const Coordinate& p1, const Coordinate& p2, const Coordinate& p3)
{
	// ʹ��б�ʷ��ж������Ƿ���  
	// ���б����ͬ�����������غϣ����ݲΧ�ڣ��������ǹ���  
	constexpr double epsilon = std::numeric_limits<double>::epsilon();
	double dx1 = p2.x() - p1.x();
	double dy1 = p2.y() - p1.y();
	double dx2 = p3.x() - p1.x();
	double dy2 = p3.y() - p1.y();

	// ʹ��б�����ж��Ƿ���  
	if (std::abs(dx1 * dy2 - dx2 * dy1) < epsilon) {
		return true;
	}
	// ����Ƿ�ӽ�0��������Ƿ����غϵĵ㣩  
	if (std::abs(dx1) < epsilon && std::abs(dx2) < epsilon) {
		// ���x���궼�ӽ���ͬ������y�����Ƿ�Ҳ�ӽ���ͬ  
		return std::abs(dy1) < epsilon && std::abs(dy2) < epsilon;
	}

	return false;
}

//���Ƕȱ�׼����0����֮�䡣
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

//����angle��ת����A
Coordinate Tools::rotateNormal(const Coordinate& A, const double angle) {
	double newX = A.x() * cos(angle) - A.y() * sin(angle);
	double newY = A.x() * sin(angle) + A.y() * cos(angle);
	return Coordinate(newX, newY);
}

/*-------------------------------------------------------------------------------
*  Function: genLinePoints()
*  Purpose: �ڸ������������м���ӵ�
*  Usage:
*	 genLinePoints(coor1, coor2, 0.0);
*	 genLinePoints(coor1, coor2, 0.8);
*  Arguments:
*    pt1 - ��ʼ��
*    pt2 - ������
*    delta - ���Ƴ��ȣ�������Ϊ0.0ʱ��������ӵ㡣
*  Returns:
*    null
-------------------------------------------------------------------------------*/
std::vector<Coordinate> Tools::genLinePoints(const Coordinate& pt1, const Coordinate& pt2, const double delta)
{
	const double EPSILON = 1e-12;//���ڱ������ڸ������ı�ʾ���µıȽ�����,������������������ȵ��������������������������ᷢ����
	if (delta < EPSILON)
		return { };

	std::vector<Coordinate> tempPointVec;
	double dis = distance(pt1, pt2);
	double dx = (pt2.x() - pt1.x()) * delta / dis;
	double dy = (pt2.y() - pt1.y()) * delta / dis;
	tempPointVec.push_back(pt1);
	//��ʣ�µ��߶γ��ȴ��ڵ���delta��1.5�����ͼ������ɵ㡣��������Ϊ�˱������һ�������pt2̫����
	while (distance(tempPointVec[tempPointVec.size() - 1], pt2) > 1.5 * delta + EPSILON)
		tempPointVec.push_back({ tempPointVec[tempPointVec.size() - 1].x() + dx,tempPointVec[tempPointVec.size() - 1].y() + dy, 0, CoorType::patch });
	tempPointVec.push_back(pt2);
	return tempPointVec;
}

//������������㹹�ɵļнǣ�ǰ���趨B��һ����A���C��֮��
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
	// ����͹����ε����ĵ�
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


// ����һ���������������������ĵ�ļ���
double Tools::calculateAngle(const Coordinate& point, const Coordinate& center) {
	double angle = atan2(point.y() - center.y(), point.x() - center.x());
	// ���Ƕȵ�����0��2pi֮��
	if (angle < 0) 
		angle += 2 * PI;
	return angle;
}

bool Tools::isVectorCollinear(const Coordinate& a, const Coordinate& b)
{
	// ��һ������:a,b
	double epsilon = 1e-8;

	// ������
	double dot_product = a.x() * b.x() + a.y() * b.y() + a.z() * b.z();

	// ������Ƿ�ӽ�1��-1
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