/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
!	@file		mesh.h
!	@brief	more mesh info
!	@author	liu guanying.
\*---------------------------------------------------------------------------*/
#pragma once
#include "mesh/include/node.h"
#include"mesh/include/element.h"
#include <vector>
#include <memory>

class Node;
class Element;

class Mesh
{
public:
	Mesh();
	~Mesh() {};
	using nodePtr = std::shared_ptr<Node>;
	using elementPtr= std::shared_ptr<Element>;
	nodePtr& getnode(size_t index);
	elementPtr& getelement(size_t index);
	std::vector<nodePtr>& getnodelist();
	std::vector<elementPtr>& getelementlist();
	size_t getnodenum();
	size_t getelementnum();
	void calCoordTrans2ord();

private:
	enum meshtype { uniform, semicircular };
	std::vector<nodePtr> nodelist_;
	std::vector<elementPtr> elementlist_;


};


