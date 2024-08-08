#include "flowControl/include/flowControl.h"
#include "iniCondition/include/GlobalData.h"
#include "utility/include/log.h"
#include "outfile/include/dumpfile.h"

FlowControl::FlowControl(const std::string& inname)
{
	mesh_ = std::make_shared<Mesh>();
	GlobalData::Init(inname);
	std::string flag = GlobalData::GetString("meshKind");
	if (flag == "Uniform")
		basemesh_ = std::make_shared<UniformMesh>(mesh_);
	else if (flag == "Semicircular")
		basemesh_ = std::make_shared<SemicircularMesh>(mesh_);
	else
		spdlog::error("the mesh key word is wrng!");
}

void FlowControl::preProcess()
{
	spdlog::info("start gen mesh.");
	basemesh_->genmesh();
	spdlog::info("currrent mesh generation is complete. ");
}

void FlowControl::postProcess()
{
	spdlog::info("start output mesh.");
	DumpFile dumpfile(mesh_);
	dumpfile.outputmesh();
	spdlog::info("Mesh dump file complete. ");
}

void FlowControl::start()
{
	//to do solver...
}
