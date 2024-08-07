#include "flowControl/include/flowControl.h"
#include "iniCondition/include/GlobalData.h"
#include "utility/include/tools.h"
#include "mesh/include/genmesh.h"
#include "utility/include/log.h"

#include <fstream>

FlowControl::FlowControl(const std::string& inname)
{
	inputlFileName_ = inname;
	mesh_ = std::make_shared<Mesh>();
	GlobalData::Init();
}

void FlowControl::preProcess()
{
	readconfig();
}

void FlowControl::postProcess()
{
	outputMesh();
}

void FlowControl::start()
{
	spdlog::info("start gen mesh.");
	Genmesh genmesh(mesh_);
	genmesh.gennode();
	genmesh.genelement();
	spdlog::info("The mesh generation is complete. ");
}

void FlowControl::readconfig()
{
	std::ifstream fin("../../config/" + inputlFileName_);
	std::string line;
	std::string dataType;
	std::string dataName;
	std::string dataValue;

	while (std::getline(fin, line))
	{
		std::string separator = " =\r\n\t#$,;\"";
		if (line.empty())
			continue;
		size_t id = line.find_first_not_of(separator);
		if (id == std::string::npos)
			continue;
		//删除前方没用的信息
		line.erase(line.begin(), line.begin() + id);
		//注释行
		if (line[0] == '!' || line[0] == '！' || line[0] == '#' || line[0] == '/')
			continue;
		id = line.find_first_of(separator);
		dataType = line.substr(0, id);
		line.erase(0, id);
		line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
		if (line.empty())
			continue;
		id = line.find_first_of("=");
		if (id == std::string::npos)
			continue;
		dataName = line.substr(0, id);
		dataValue = line.substr(id + 1);
		if (dataType == "string")
			GlobalData::Update(dataName, dataValue);
		else if (dataType == "double")
			GlobalData::Update(dataName, stod(dataValue));
		else if (dataType == "int")
			GlobalData::Update(dataName, stoi(dataValue));
		else
		{
			spdlog::warn("Unspported Data Type:{}, Name:{}, Value:{}", dataType, dataName, dataValue);
		}
	}
}

void FlowControl::outputMesh()
{
	std::string basePath = "tempData";
	Tools::createFolder(basePath);
	std::ofstream os;
	os.open("tempData/Mesh_1.dat");
	using std::endl;
	os << "TITLE = MESH" << std::endl;
	os << "variables=x,y" << std::endl;

	size_t outElementnum = mesh_->getelementnum();
	size_t outNodenum = mesh_->getnodenum();
	if (mesh_->getnodenum() == 0) {
		spdlog::error("the node num is 0!");
		return;
	}
	os << "ZONE t=\" MESH \"" << std::endl;
	os << "N=" << outNodenum << ", E=" << outElementnum << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
	for (size_t iNode = 0; iNode < outNodenum; ++iNode)
	{
		auto& currentNode = mesh_->getnode(iNode);
		os << currentNode->getCoord().x() << "  " << currentNode->getCoord().y() << "  " << endl;
	}

	for (size_t iElem = 0; iElem < outElementnum; ++iElem)
	{
		auto& currentElement = mesh_->getelement(iElem);
		vector<size_t>updateNodeIndex;
		const auto& elementNode = currentElement->getnodelist();
		for (size_t iNode = 0; iNode < elementNode.size(); ++iNode)
		{
			auto& currenNode = elementNode[iNode];
			os << currenNode->getIndex() + 1 << "  ";
		}
		if (elementNode.size() == 3)
			os << elementNode[2]->getIndex() + 1;
		os << endl;
	}
	os.close();
	return;
}