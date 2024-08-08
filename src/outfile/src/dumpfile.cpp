
#include "utility/include/log.h"
#include "outfile/include/dumpfile.h"
#include "utility/include/tools.h"
#include "iniCondition/include/GlobalData.h"


DumpFile::DumpFile(std::shared_ptr<Mesh>& outmesh) :mesh_(outmesh) 
{
	std::string basePath = "tempdata";
	Tools::createFolder(basePath);
}

void DumpFile::outputmesh()
{
	std::ofstream os;
	const std::string meshtag = GlobalData::GetString("meshKind");
	os.open("tempData/Mesh_" + meshtag + ".dat");
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