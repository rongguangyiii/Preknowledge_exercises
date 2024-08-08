
#include"mesh/include/meshsSelect.h"
#include "iniCondition/include/GlobalData.h"
#include"Coordinate/include/Coordinate.h"
#include "mesh/include/node.h"
#include "utility/include/log.h"
#include "utility/include/tools.h"


void UniformMesh::gennode()
{
	const size_t mx = GlobalData::GetInt("xNodeNum");
	const size_t my = GlobalData::GetInt("yNodeNum");
	const double xa = GlobalData::GetDouble("xMin");
	const double ya = GlobalData::GetDouble("yMin");
	const double xb = GlobalData::GetDouble("xMax");
	const double yb = GlobalData::GetDouble("yMax");
	const double dx = (xb - xa) / (mx - 1);
	const double dy = (yb - ya) / (my - 1);
	auto& meshNodevec = mesh_->getnodelist();
	size_t count = 0;
	for (size_t j = 0; j < my; j++)
	{
		for (size_t i = 0; i < mx; i++)
		{
			const double x = xa + i * dx;
			const double y = ya + j * dy;
			meshNodevec.emplace_back(std::make_shared<Node>(Coordinate{ x, y }, count++));
		}
	}
}

void UniformMesh::genelement()
{
	auto& meshElevec = mesh_->getelementlist();
	const size_t mx = GlobalData::GetInt("xNodeNum");
	const size_t my = GlobalData::GetInt("yNodeNum");
	size_t eleNum = 0;
	for (size_t j = 0; j < my - 1; j++)
	{
		//spdlog::debug("第 {} 行的单元",j);
		for (size_t i = 0; i < mx - 1; i++)
		{
			const size_t perNodeindex = (j * mx) + i;
			/*spdlog::debug("当前单元的四个节点编号：{},{},{},{}", perNodeindex, perNodeindex + 1,
				perNodeindex + mx + 1, perNodeindex + mx);*/
			auto& node_1 = mesh_->getnode(perNodeindex);
			auto& node_2 = mesh_->getnode(perNodeindex + 1);
			auto& node_3 = mesh_->getnode(perNodeindex + mx + 1);
			auto& node_4 = mesh_->getnode(perNodeindex + mx);
			//std::vector<std::shared_ptr<Node>> perNodevec = { node_1, node_2, node_3, node_4 };
			meshElevec.emplace_back(std::make_shared<Element>(eleNum++, node_1, node_2, node_3, node_4));

		}
	}
}

void SemicircularMesh::gennode()
{
	const double R0 = 1;
	const double R1 = 5;
	size_t mx = GlobalData::GetInt("xNodeNum");
	size_t my = GlobalData::GetInt("yNodeNum");
	auto& meshNodevec = mesh_->getnodelist();
	size_t count = 0;
	for (size_t j = 0; j < my; j++)
	{
		for (size_t i = 0; i < mx; i++)
		{
			double R = R0 + j * (R1 - R0) / (my - 1);
			double theta = Tools::PI * i / (mx - 1);
			double x = -R * sin(theta);
			double y = -R * cos(theta);
			meshNodevec.emplace_back(std::make_shared<Node>(Coordinate{ x, y }, count++));
		}
	}

}

void SemicircularMesh::genelement()
{
	auto& meshElevec = mesh_->getelementlist();
	const size_t mx = GlobalData::GetInt("xNodeNum");
	const size_t my = GlobalData::GetInt("yNodeNum");
	size_t eleNum = 0;
	for (size_t j = 0; j < my - 1; j++)
	{
		//spdlog::debug("第 {} 行的单元",j);
		for (size_t i = 0; i < mx - 1; i++)
		{
			const size_t perNodeindex = (j * mx) + i;
			/*spdlog::debug("当前单元的四个节点编号：{},{},{},{}", perNodeindex, perNodeindex + 1,
				perNodeindex + mx + 1, perNodeindex + mx);*/
			auto& node_1 = mesh_->getnode(perNodeindex);
			auto& node_2 = mesh_->getnode(perNodeindex + 1);
			auto& node_3 = mesh_->getnode(perNodeindex + mx + 1);
			auto& node_4 = mesh_->getnode(perNodeindex + mx);
			//std::vector<std::shared_ptr<Node>> perNodevec = { node_1, node_2, node_3, node_4 };
			meshElevec.emplace_back(std::make_shared<Element>(eleNum++, node_1, node_2, node_3, node_4));

		}
	}
}