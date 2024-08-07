/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
!	@file		node.h
!	@brief	node info
!	@author	liu guanying.
\*---------------------------------------------------------------------------*/

#pragma once
#include "Coordinate/include/Coordinate.h"
#include <cstddef>  // For size_t

class Node
{
public:
    Node(const Coordinate& coor, const size_t index) : nodecoor_(coor), nodeindex_(index) {}
    Node(const Node& other) : nodecoor_(other.nodecoor_), nodeindex_(other.nodeindex_) {}
    Node& operator=(const Node& other);
    ~Node() {}
    const Coordinate& getCoord() const { return nodecoor_; }
    size_t getIndex() const { return nodeindex_; }
    bool operator==(const Node& other) const;
private:
    Coordinate nodecoor_;
    size_t nodeindex_;
};
