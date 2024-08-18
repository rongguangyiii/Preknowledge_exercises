/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
!	@file		node.h
!	@brief	node info
!	@author	liu guanying.
\*---------------------------------------------------------------------------*/

#pragma once
#include "Coordinate/include/Coordinate.h"
#include "mesh/include/mesh.h"
#include <cstddef>  // For size_t
#include <map>
#include <string>
#include <memory>

class Mesh;

class Node
{
public:
	enum Nodetype { unset, inner, leftbound, rightbound, upbound, downbound, corner};
	Node(const Coordinate& coor, const size_t index, const std::shared_ptr<Mesh>& ptr, Nodetype type = unset);
	Node() = delete;
	Node(const Node& other) = delete;
	Node& operator=(const Node& other);
	~Node() {}
	const Coordinate& getCoord() const { return nodecoor_; }
	size_t getIndex() const { return nodeindex_; }
	bool operator==(const Node& other) const;
	std::vector<double>& getPrimitiveValue() { return primitiveValue_; }
	std::vector<double>& getConservedValue() { return conservedValue_; }
	void toConservedform();
	void toPrimitiveform();
	void setType(Nodetype typein) { nodetype_ = typein; }
	void setCoordTrans(double inksi_x, double inksi_y, double ineta_x, double ineat_y, double inJacobian);
	Nodetype getType()const { return nodetype_; }
	std::vector<std::shared_ptr<Node>>& getNeighbor() { return neighborNode_; }
private:
	Coordinate nodecoor_;
	size_t nodeindex_;
	Nodetype nodetype_;
	std::vector<double> primitiveValue_;
	std::vector<double> conservedValue_;
	std::vector<std::shared_ptr<Node>> neighborNode_;
	std::shared_ptr<Mesh> meshptr_;
public:
	struct TransCoef
	{
		double ksi_x_;
		double ksi_y_;
		double ksi_t_;
		double eta_x_;
		double eta_y_;
		double eta_t_;
		double jacob_;
	}  coordTrans_;
	const TransCoef& getcoordTrans()const { return coordTrans_; }
	void setConservedValue(const std::vector<double>& value) { conservedValue_ = value; }
	void setPrimitiveValue(const std::vector<double>& value) { primitiveValue_ = value; }
	double& r() { return primitiveValue_[0]; }
	double& u() { return primitiveValue_[1]; }
	double& v() { return primitiveValue_[2]; }
	//double& w() { return primitiveValue_[3]; }
	double& p() { return primitiveValue_[3]; }

	double& rho() { return conservedValue_[0]; }
	double& rhou() { return conservedValue_[1]; }
	double& rhov() { return conservedValue_[2]; }
	//double& rhow() { return conservedValue_[3]; }
	double& E() { return conservedValue_[3]; }
	void checkNaN(const double value);
};
