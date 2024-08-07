/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
!	@file		Coor.h
!	@brief	Coordinate information.
!	@author	liu guanying.
\*---------------------------------------------------------------------------*/
#pragma once
#include<vector>

using std::vector;

enum class CoorType
{
	ordinary,
	patch,
};

class Coordinate
{
public:
	Coordinate(double x = 0, double y = 0, double z = 0, CoorType type = CoorType::ordinary) :x_(x), y_(y), z_(z), coorType_(type) {}
	Coordinate(const Coordinate& pt) :x_(pt.x_), y_(pt.y_), z_(pt.z_), coorType_(pt.coorType_) {}
	const double& x() const { return x_; }
	const double& y() const { return y_; }
	const double& z() const { return z_; }
	void SetCoord(const double& x = 0, const double &y = 0, const double &z = 0);
	bool operator ==(const Coordinate& pt)const;
	double norm() const;
	void setCoorType(CoorType type) { coorType_ = type; }
	CoorType getCoorType() const { return coorType_; }
	Coordinate& operator+=(const Coordinate& pt);
	Coordinate& operator-=(const Coordinate& pt);
	Coordinate& operator/=(double scalar);
	Coordinate operator-(const Coordinate& pt)const;
	Coordinate operator+(const Coordinate& pt)const;
	bool operator<(const Coordinate& other) const;
private:
	double x_, y_, z_;
	CoorType coorType_;
};

struct CoorHash {
	std::size_t operator()(const Coordinate& pt) const;
};

