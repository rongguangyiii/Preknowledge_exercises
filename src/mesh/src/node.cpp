#include "mesh/include/node.h"

Node& Node::operator=(const Node& other) {
    if (this != &other) { 
        nodecoor_ = other.nodecoor_;
        nodeindex_ = other.nodeindex_;
    }
    return *this;
}

bool Node::operator==(const Node& other) const {
    return (nodecoor_ == other.nodecoor_) && (nodeindex_ == other.nodeindex_);
}