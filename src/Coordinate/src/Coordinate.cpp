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
#include"Coordinate/include/Coordinate.h"
#include<algorithm>
#include<cmath>
#include<stdexcept>

void Coordinate::SetCoord(const double& x, const double& y, const double& z)
{
	x_ = x;
	y_ = y;
	z_ = z;
}

bool Coordinate::operator ==(const Coordinate& pt)const
{
	double delta = 1e-14;
	return fabs(x_ - pt.x_) < delta && fabs(y_ - pt.y_) < delta && fabs(z_ - pt.z_) < delta;
};

bool Coordinate::operator<(const Coordinate& other) const {
	if (x_ != other.x_) return x_ < other.x_;
	if (y_ != other.y_) return y_ < other.y_;
	return z_ < other.z_;
}

Coordinate& Coordinate::operator+=(const Coordinate& pt)
{
	x_ += pt.x();
	y_ += pt.y();
	z_ += pt.z();
	return *this;
}
Coordinate& Coordinate::operator-=(const Coordinate& pt)
{
	x_ -= pt.x();
	y_ -= pt.y();
	z_ -= pt.z();
	return *this;
}
Coordinate Coordinate::operator-(const Coordinate& pt) const
{
	return Coordinate(x_ - pt.x(), y_ - pt.y(), z_ - pt.z());
}
Coordinate Coordinate::operator+(const Coordinate& pt) const
{
	return Coordinate(x_ + pt.x(), y_ + pt.y(), z_ + pt.z());
}
Coordinate& Coordinate::operator/=(double scalar)
{
	if (scalar != 0)
	{
		x_ /= scalar;
		y_ /= scalar;
		z_ /= scalar;
	}
	else
	{
		// 处理除数为零的情况  
		throw std::invalid_argument("Cannot divide by zero.");
	}
	return *this;
}

double Coordinate::norm() const
{
	return std::sqrt(x_ * x_ + y_ * y_);
}

size_t CoorHash::operator()(const Coordinate& pt) const {

	// 直接使用浮点数来进行哈希运算可能会导致不同的浮点数具有相同的哈希值，或者相近的浮点数具有截然不同的哈希值。 
	// 为减少由于浮点数精度带来的影响，使用固定的精度来量化坐标值 
	const double precision = 1e-7; // 例如，使用小数点后8位的精度，由于文件数据只有小数点后7位.  
	double qx = std::round(pt.x() / precision) * precision;
	double qy = std::round(pt.y() / precision) * precision;
	double qz = std::round(pt.z() / precision) * precision;

	// 使用量化后的坐标值来计算哈希  
	auto hash1 = std::hash<double>{}(qx);
	auto hash2 = std::hash<double>{}(qy);
	auto hash3 = std::hash<double>{}(qz);

	// 结合三个哈希值以产生最终的哈希  
	return hash1 ^ (hash2 << 1) ^ (hash3 << 2);
}

