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
		// �������Ϊ������  
		throw std::invalid_argument("Cannot divide by zero.");
	}
	return *this;
}

double Coordinate::norm() const
{
	return std::sqrt(x_ * x_ + y_ * y_);
}

size_t CoorHash::operator()(const Coordinate& pt) const {

	// ֱ��ʹ�ø����������й�ϣ������ܻᵼ�²�ͬ�ĸ�����������ͬ�Ĺ�ϣֵ����������ĸ��������н�Ȼ��ͬ�Ĺ�ϣֵ�� 
	// Ϊ�������ڸ��������ȴ�����Ӱ�죬ʹ�ù̶��ľ�������������ֵ 
	const double precision = 1e-7; // ���磬ʹ��С�����8λ�ľ��ȣ������ļ�����ֻ��С�����7λ.  
	double qx = std::round(pt.x() / precision) * precision;
	double qy = std::round(pt.y() / precision) * precision;
	double qz = std::round(pt.z() / precision) * precision;

	// ʹ�������������ֵ�������ϣ  
	auto hash1 = std::hash<double>{}(qx);
	auto hash2 = std::hash<double>{}(qy);
	auto hash3 = std::hash<double>{}(qz);

	// ���������ϣֵ�Բ������յĹ�ϣ  
	return hash1 ^ (hash2 << 1) ^ (hash3 << 2);
}

