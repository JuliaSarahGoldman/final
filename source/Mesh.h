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
    Array<Vector3int32> m_triArray;

public:
    void addVertex(const Vector3& vertex);
    void addVertex(const Array<Vector3>& vertexList);

    void addIndex(int index);
    void addIndex(const Array<int>& indexList);

    void addVertex(const Vector3& vertex, int index);
    void addVertex(const Array<Vector3>& vertexList, const Array<int>& indexList);

    void computeAdjacency(Array<MeshAlg::Face>& faceArray, Array<MeshAlg::Edge>& edgeArray, Array<MeshAlg::Vertex>& vertexArray);

    void computeNormals(const Array<MeshAlg::Face>& faceArray, const Array<MeshAlg::Edge>& edgeArray, const Array<MeshAlg::Vertex>& vertexArray,
    Array<Vector3>& vertexNormalArray, Array<Vector3>& faceNormalArray);

    void collapseEdges();

    void bevelEdges();

    void toObj(String filename);

    Mesh(const Array<Vector3>& vertexPositions, const Array<Vector3int32>& triArray);
    ~Mesh();


};

