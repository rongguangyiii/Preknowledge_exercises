/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
!	@file		genmesh.h
!	@brief	mesh generation.
!	@author	liu guanying.
!	@date 2024.8.6
\*---------------------------------------------------------------------------*/
#pragma once
#include "mesh/include/mesh.h"
#include<string>
#include<memory>

class Genmesh
{
public:
	Genmesh(std::shared_ptr<Mesh>& it) :mesh_(it) {}
	~Genmesh() {}
	void gennode();
	void genelement();

private:
	std::shared_ptr<Mesh>& mesh_;

};
