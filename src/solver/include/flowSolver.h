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
#include "iniCondition/include/GlobalData.h"

class FlowSlover
{
public:
	FlowSlover(std::shared_ptr<Mesh>& mesh) : mesh_(mesh), caltime_(0.0),residual_(1.0) {}
	~FlowSlover() {}
	void solve();
	void Initializeflow();
	void timeAdvance();
	void writeTecplotFile(const size_t timestep)const;
	std::vector<double> computeRHS(std::shared_ptr<Node>& it);
	void computeResidual();
	void BoundExtrapolate();
	void solverEuler();

private:
	void split(const std::vector<std::shared_ptr<Node>>& nodeTemplate, std::vector<double>& fluxplus, std::vector<double>& fluxminus);
	void steger(const std::shared_ptr<Node>& nodeptr, const Node::TransCoef& trans, std::vector<double>& fp, std::vector<double>& fn);
	void vanLeer(const std::shared_ptr<Node>& nodeptr, const Node::TransCoef& trans, std::vector<double>& fp, std::vector<double>& fn);
	void checkNaN(double value1);
	void checkNaN(std::vector<double> value1);
	double Dt();
	std::shared_ptr<Mesh>& mesh_;
	double residual_;
	double caltime_;

};
