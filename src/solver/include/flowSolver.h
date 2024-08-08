/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
 @brief	   flow Slover Process.
 @author	liu guanying.
 @date     2024.8.8
\ * -------------------------------------------------------------------------- - */
#include "mesh/include/mesh.h"

class FlowSlover
{
public:
	FlowSlover(std::shared_ptr<Mesh>& mesh) : mesh_(mesh) {}
	~FlowSlover() {}
	void solve();
	void Initializeflow();
private:
	std::shared_ptr<Mesh>& mesh_;

};
