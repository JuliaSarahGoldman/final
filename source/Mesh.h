/**
  \file Mesh.h
   
   Defines mesh and its operators, including edge collapse and beveling.
 */
/*
    John Freeman 
    Jose Rivas-Garcia 
    Julia Goldman 
    Matheus de Carvalho Souza
*/
#pragma once
#include <G3D/G3DAll.h>


class Mesh { 
protected:
    Array<Vector3> m_vertexPositions;
    Array<int> m_indexArray;
    Array<MeshAlg::Face> m_faceArray; 
    Array<MeshAlg::Edge> m_edgeArray; 
    Array<MeshAlg::Vertex> m_vertexArray;

public:
    void addVertex(const Vector3& vertex); 
    void addVertex(const Array<Vector3>& vertexList); 
    
    void addIndex(int index);
    void addIndex(const Array<int>& indexList);
   
    void addVertex(const Vector3& vertex, int index);
    void addVertex(const Array<Vector3>& vertexList, const Array<int>& indexList); 
   
    void computeAdjacency(); 

    void collapseEdges(const std::function<void (Array<MeshAlg::Edge>&)>& sort); 

    void bevelEdges();

    Mesh(const shared_ptr<TriTree>& triTree); 
    ~Mesh();


};

