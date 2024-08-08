
#include "solver/include/flowSolver.h"
#include "mesh/include/mesh.h"
#include <cmath> 

void FlowSlover::solve()
{
	Initializeflow();
	//to do ....
}

void FlowSlover::Initializeflow()
{
	std::vector<Mesh::nodePtr>& nodesVe = mesh_->getnodelist();
	for (auto& curnode : nodesVe)
	{
		double x = curnode->getCoord().x();
		auto& privalue = curnode->getPrimitiveValue();
		if (x < 2.0 || std::fabs(x - 2.0) < 1E-6)
			privalue["rho"] = 2.8;
		else
			privalue["rho"] = 1.4;
		privalue["u"] = 0.5;
		privalue["v"] = 0.0;
		privalue["w"] = 0.0;
		privalue["p"] = 1.0;
		curnode->toConservedform();
	}
}

