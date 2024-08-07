#include"mesh/include/element.h"

Element::Element(std::vector<nodePtr>& nodevec, size_t index)
{
	eleindex_ = index;
	eleNodevec_ = nodevec;
}

Element::Element(size_t index,nodePtr& node_1, nodePtr& node_2, nodePtr& node_3, nodePtr& node_4)
{
	eleindex_ = index;
	eleNodevec_ = { node_1, node_2, node_3, node_4 };

}

Element::~Element()
{
}
std::vector<Element::nodePtr>& Element::getnodelist()
{
	return eleNodevec_;
}

