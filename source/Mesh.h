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

class Mesh : public ReferenceCountedObject {
protected:
    Array<Vector3> m_vertexPositions;
    Array<int> m_indexArray;
    Array<Vector3int32> m_triArray;

    int edgeLength(const MeshAlg::Edge& edge);
    int Mesh::edgeLength(int i0, int i1);

    MeshAlg::Edge findMinEdge(const Array<MeshAlg::Edge>& edges);

    static void Mesh::merge(SmallArray<float, 6>& data, SmallArray<float, 6>& temp, int low, int middle, int high, SmallArray<int, 6>& along, SmallArray<int, 6>& temp2);
    Array<MeshAlg::Edge> minTriEdges();
    static void Mesh::mergeSortRecursive(SmallArray<float, 6>& data, SmallArray<float, 6>& temp, int low, int high, SmallArray<int, 6>& along, SmallArray<int, 6>& temp2);
    static void Mesh::mergeSort(SmallArray<float, 6>& data, SmallArray<int, 6>& along);

    int findMinEdge(int minIndex, int maxIndex);

    /** Called by collapseEdges() */
    Array<Array<int>> toCollapse(int regionSize);





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

    /** As outlined by Stanford Graphics http://graphics.stanford.edu/courses/cs468-10-fall/LectureSlides/08_Simplification.pdf
        Calls toCollapse()
        Collapses one edge for each block of regionSize indices */
    void collapseEdges(int regionSize);

    //bump is how much we expand the panet's radius by
    //Only works on completely closed and welded meshes.
    void bevelEdges(float bump);

    void toObj(String filename);

    shared_ptr<Model> toArticulatedModel(String name);

    static std::shared_ptr<Mesh> create(const Array<Vector3>& vertexPositions, const Array<Vector3int32>& triArray);
    Mesh(const Array<Vector3>& vertexPositions, const Array<Vector3int32>& triArray);
    ~Mesh();


};

