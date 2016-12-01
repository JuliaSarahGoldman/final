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
    Array<Vector3> m_vertexNormals;
    bool m_hasFakeNormals = false;

    int edgeLength(const MeshAlg::Edge& edge);
    int Mesh::edgeLength(int i0, int i1);

    /** Called by collapseEdges()
        true if not a boundary edge, if collapsing it doesn't create a manifold and doesn't flip face normals*/
    bool isCollapsable(const MeshAlg::Edge& edge, const Array<MeshAlg::Face>& faces, const Array<MeshAlg::Edge>& edges, const Array<MeshAlg::Vertex>& vertices);

    /** Returns true if the angle between the normal of the adjacent faces of elem1 is greater than that of 
        the adjacent faces of elem2 */
    //bool greaterAngle(const MeshAlg::Edge& elem1, const MeshAlg::Edge& elem2);

    MeshAlg::Edge Mesh::toCollapse(const Array<MeshAlg::Edge>& data);

    void merge(Array<MeshAlg::Edge>& data, const Array<MeshAlg::Edge>& temp, int low, int middle, int hight);
    void mergeSortRecursive(Array<MeshAlg::Edge>& data, Array<MeshAlg::Edge>& temp, int low, int hight);
    void mergeSort(Array<MeshAlg::Edge>& data);

public:
    static int edgeLength(const MeshAlg::Edge& edge, const Array<Vector3>& vertexArray);

    void computeAdjacency(Array<MeshAlg::Face>& faceArray,
        Array<MeshAlg::Edge>& edgeArray = Array<MeshAlg::Edge>(),
        Array<MeshAlg::Vertex>& vertexArray = Array<MeshAlg::Vertex>());

    /** pre: faceArray, edgeArray and vertexArray computed by computeAdjacency
    post: vertexNormalArray and faceNormalArray filled with appropriate normal values */
    void computeNormals(const Array<MeshAlg::Face>& faceArray, const Array<MeshAlg::Edge>& edgeArray, const Array<MeshAlg::Vertex>& vertexArray,
        Array<Vector3>& vertexNormalArray, Array<Vector3>& faceNormalArray);

    void computeFaceNormals(Array<Vector3>& faceNormals, bool normalize = true);

    /** As outlined by Stanford Graphics http://graphics.stanford.edu/courses/cs468-10-fall/LectureSlides/08_Simplification.pdf
        Calls isCollapsable()
        Collapses the numEdges most unimportant edges*/
    void collapseEdges(int numEdges);

    //bump is how much we expand the panet's radius by
    //Only works on completely closed and welded meshes.
    void bevelEdges(float bump);
    void bevelEdges2(float bump);

    void toObj(String filename);

    shared_ptr<Model> toArticulatedModel(String name, Color3& color) const;
    shared_ptr<Model> toArticulatedModel(String name, String anyStr, int width, int height) const;

    static std::shared_ptr<Mesh> Mesh::create(const String& filename);
    static std::shared_ptr<Mesh> create(const Array<Vector3>& vertexPositions, const Array<Vector3int32>& triArray);

    Mesh(const Array<Vector3>& vertexPositions, const Array<Vector3int32>& triArray);
    Mesh(const String& filename);

    ~Mesh();


};

