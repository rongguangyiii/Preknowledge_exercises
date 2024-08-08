/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
!	@file		DUMPFILE.h
!	@brief	Output grid information to display in tecplot.
!	@author	liu guanying.
!   @date   2024.8.8
\*---------------------------------------------------------------------------*/
#pragma once
#include "mesh/include/mesh.h"
#include <fstream>
#include <memory>

class DumpFile
{
public:
	DumpFile(std::shared_ptr<Mesh>& outmesh);
	~DumpFile() {}
	void outputmesh();
private:
	std::shared_ptr<Mesh>& mesh_;
};