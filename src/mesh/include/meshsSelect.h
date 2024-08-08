/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
!	@file	element.h
!	@brief	elenent info
!	@author	liu guanying.
！  @date 2024.8.7
\*---------------------------------------------------------------------------*/
#pragma once
#include "mesh/include/mesh.h"
#include <memory>

// MeshsSelect 基类
class MeshsSelect
{
public:
    MeshsSelect(std::shared_ptr<Mesh>& mesh) : mesh_(mesh) {}
    virtual ~MeshsSelect() = default;

    virtual void gennode() = 0;
    virtual void genelement() = 0;
    //void genmesh();
protected:
    std::shared_ptr<Mesh>& mesh_;
};

// UniformMesh 派生类
class UniformMesh : public MeshsSelect
{
public:
    UniformMesh(std::shared_ptr<Mesh>& mesh) : MeshsSelect(mesh) {}
    ~UniformMesh() override = default;
    void gennode() override;
    void genelement() override;
};

// SemicircularMesh 派生类
class SemicircularMesh : public MeshsSelect
{
public:
    SemicircularMesh(std::shared_ptr<Mesh>& mesh) : MeshsSelect(mesh) {}
    ~SemicircularMesh() override = default;
    void gennode() override;
    void genelement() override;
};
