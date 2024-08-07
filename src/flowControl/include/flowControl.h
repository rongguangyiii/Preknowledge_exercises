/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
 @file		FLOWCONTRLO.h
 @brief	   flow Control Process.
 @author	liu guanying.
 @date     2024.8.6
\ * -------------------------------------------------------------------------- - */

#ifndef FLOWCONTROL_H
#define FLOWCONTROL_H
#include "mesh/include/mesh.h"
#include "utility/include/log.h"
#include <string>

class FlowControl
{
public:
	FlowControl(const std::string& inname);
	~FlowControl() {}

	void start();
	void preProcess();
	void postProcess();
private:
	void readconfig();
	void outputMesh();
	std::string inputlFileName_;
	std::shared_ptr<Mesh> mesh_;
};



#endif