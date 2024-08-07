/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
!	@file	tools.h
!	@brief	useful tools function.
!	@author	Liu Guangying.
!   @date  2024.08.06
!   @location DaLian
\*---------------------------------------------------------------------------*/
#ifndef TOOLS_H  
#define TOOLS_H  
#include <vector>
#include <string>
#include "mesh/include/mesh.h"
#include"Coordinate/include/Coordinate.h"

namespace Tools
{
	const double PI = 3.141592653589793238462643383279502884197169399;//圆周率pi

	void createFolder(const std::string& basePath);

	const Coordinate Cross(const Coordinate& vA, const Coordinate& vB);

	bool isPointInPolygon_moreInner(const Coordinate& p, const std::vector<Coordinate>& polygon);

	bool isPointInPolygon_moreOuter(const Coordinate& p, const std::vector<Coordinate>& polygon);

	//是否在两点之间
	bool IsInTwoPoint(const Coordinate& P1, const Coordinate& P2, const Coordinate& judgeP);

	//点到直线的距离
	double distanceToLine(const Coordinate& A, const Coordinate& B, const Coordinate& P);

	bool IsSameDirection(const Coordinate& vec1, const Coordinate& vec2);

	void OutCompleteMesh2Tecplot(const std::shared_ptr<Mesh>& grid);

	bool isCollinear(const Coordinate& p1, const Coordinate& p2, const Coordinate& p3);

	bool isCocircular(const Coordinate& p1, const Coordinate& p2, const Coordinate& p3);

	Coordinate getCocircularCenter(const Coordinate& p1, const Coordinate& p2, const Coordinate& p3);

	double normalizeAngle(double angleRad);

	std::vector<Coordinate> genLinePoints(const Coordinate& pt1, const Coordinate& pt2, const double delta);

	Coordinate rotateNormal(const Coordinate& A, const double angle);

	const double angle_ABC(const Coordinate& ptA, const Coordinate& ptB, const Coordinate& ptC);

	Coordinate returnPolycenter(const std::vector<Coordinate>& points);

	// 定义一个函数来计算点相对于中心点的极角
	double calculateAngle(const Coordinate& point, const Coordinate& center);
	
	bool isVectorCollinear(const Coordinate& a, const Coordinate& b);

	Coordinate unitNormalize(const Coordinate& v);

	double distance(const Coordinate& pt1, const Coordinate& pt2);

	double distance(double x1, double y1, double z1, double x2, double y2, double z2);
}


#endif // TOOLS_H