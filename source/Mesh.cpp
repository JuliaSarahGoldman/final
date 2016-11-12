/** \file Mesh.cpp */
#include "Mesh.h"

Mesh::Mesh() {};
Mesh::~Mesh() {};

void Mesh::addVertex(const Vector3& vertex) {
    m_vertexPositions.append(vertex);
};

void Mesh::addVertex(const Array<Vector3>& vertexList){ 
    m_vertexPositions.append(vertexList);   
};

void Mesh::addIndex(int index){
    m_indexArray.append(index);    
};


void Mesh::addIndex(const Array<int>& indexList){
    m_indexArray.append(indexList);
};

void Mesh::addVertex(const Vector3& vertex, int index){ 
    addVertex(vertex); 
    addIndex(index);
};

void Mesh::addVertex(const Array<Vector3>& vertexList, const Array<int>& indexList){
    addVertex(vertexList); 
    addIndex(indexList);
};

void Mesh::computeAdjacency(){
    MeshAlg::computeAdjacency(m_vertexPositions, m_indexArray, m_faceArray, m_edgeArray, m_vertexArray);
};
