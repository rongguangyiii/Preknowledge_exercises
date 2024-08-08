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
#include <map>
#include <string>

class Node
{
public:
    Node(const Coordinate& coor, const size_t index);
    Node(const Node& other) : nodecoor_(other.nodecoor_), nodeindex_(other.nodeindex_) {}
    Node& operator=(const Node& other);
    ~Node() {}
    const Coordinate& getCoord() const { return nodecoor_; }
    size_t getIndex() const { return nodeindex_; }
    bool operator==(const Node& other) const;
    std::map<std::string, double>& getPrimitiveValue()  { return primitiveValue_; }
    std::map<std::string, double>& getConservedValue()  { return conservedValue_; }
    //double getPrimitiveValue(const std::string& name) const;
    //double getConservedValue(const std::string& name) const;
    void toConservedform();
    void toPrimitiveform();

private:
    void initValues();
    Coordinate nodecoor_;
    size_t nodeindex_;
    std::vector<double> neighbors_;
    std::map<std::string, double> primitiveValue_;
    std::map<std::string, double> conservedValue_;
};
