/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
!	@file	element.h
!	@brief	elenent info
!	@author	liu guanying.
£¡  @date 2024.8.6
\*---------------------------------------------------------------------------*/
#pragma once
#include "mesh/include/node.h"
#include<memory>
#include<vector>

class Element
{
public:
	using nodePtr = std::shared_ptr<Node>;

	Element(std::vector<nodePtr>& nodevec, size_t index);
	Element(size_t index, nodePtr& node_1, nodePtr& node_2, nodePtr& node_3, nodePtr& node_4);
	~Element();
	std::vector<nodePtr>& getnodelist();
private:
	std::vector<nodePtr> eleNodevec_;
	size_t eleindex_;

};

