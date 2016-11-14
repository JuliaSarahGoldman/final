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

    int edgeLength(const MeshAlg::Edge& edge);
    
    MeshAlg::Edge findMinEdge(const Array<MeshAlg::Edge>& edges);


public:
    void addVertex(const Vector3& vertex);
    void addVertex(const Array<Vector3>& vertexList);

    void addIndex(int index);
    void addIndex(const Array<int>& indexList);

    void addVertex(const Vector3& vertex, int index);
    void addVertex(const Array<Vector3>& vertexList, const Array<int>& indexList);

    void computeAdjacency(Array<MeshAlg::Face>& faceArray, Array<MeshAlg::Edge>& edgeArray, Array<MeshAlg::Vertex>& vertexArray);

    /** pre: faceArray, edgeArray and vertexArray computed by computeAdjacency
    post: vertexNormalArray and faceNormalArray filled with appropriate normal values */
    void computeNormals(const Array<MeshAlg::Face>& faceArray, const Array<MeshAlg::Edge>& edgeArray, const Array<MeshAlg::Vertex>& vertexArray,
        Array<Vector3>& vertexNormalArray, Array<Vector3>& faceNormalArray);

    void collapseEdges(int numCollapsed);

    //bump is how much we expand the panet's radius by
    void bevelEdges(float bump);

    void toObj(String filename);

    Mesh(const Array<Vector3>& vertexPositions, const Array<Vector3int32>& triArray);
    ~Mesh();


};

