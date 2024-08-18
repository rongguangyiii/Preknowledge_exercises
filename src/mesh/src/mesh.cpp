#include "mesh/include/mesh.h"
#include "iniCondition/include/GlobalData.h"


Mesh::Mesh()
{
}

Mesh::nodePtr& Mesh::getnode(size_t index)
{
	return nodelist_.at(index);
}

std::vector<Mesh::nodePtr>& Mesh::getnodelist()
{
	return nodelist_;
}

size_t Mesh::getnodenum()
{
	return nodelist_.size();
}

std::vector<Mesh::elementPtr>& Mesh::getelementlist()
{
	return elementlist_;
}

size_t Mesh::getelementnum()
{
	return elementlist_.size();
}
Mesh::elementPtr& Mesh::getelement(size_t index)
{
	return elementlist_.at(index);
}

