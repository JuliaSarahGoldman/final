/** \file Mesh.cpp */
#include "Mesh.h"

Mesh::Mesh(const Array<Vector3>& vertexPositions, const Array<Vector3int32>& triArray) {
    MeshBuilder builder;
    for(int i = 0; i < triArray.size(); ++i){
        Vector3int32 v = triArray[i]; 
        m_indexArray.append(v.x,v.y,v.z);
    }
    m_vertexPositions = vertexPositions; 
    m_triArray = triArray;
};

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

void Mesh::computeAdjacency(Array<MeshAlg::Face>& faceArray, Array<MeshAlg::Edge>& edgeArray, Array<MeshAlg::Vertex>& vertexArray){
    MeshAlg::computeAdjacency(m_vertexPositions, m_indexArray, faceArray, edgeArray, vertexArray);
};

/** pre: faceArray, edgeArray and vertexArray computed by computeAdjacency
    post: vertexNormalArray and faceNormalArray filled with appropriate normal values
*/
void Mesh::computeNormals(const Array<MeshAlg::Face>& faceArray, const Array<MeshAlg::Edge>& edgeArray, const Array<MeshAlg::Vertex>& vertexArray,
    Array<Vector3>& vertexNormalArray, Array<Vector3>& faceNormalArray){
     MeshAlg::computeNormals(m_vertexPositions, faceArray, vertexArray, vertexNormalArray, faceNormalArray);
};

void Mesh::toObj(String filename){
    TextOutput file(filename);
    //file.printf("g name\n");
        //loop to make vertices
    //debugPrintf(STR(%d\n), sizeof(m_vertexArray));
    for(int i = 0; i < m_vertexPositions.size(); ++i){
        file.printf(STR(v %f %f %f\n), m_vertexPositions[i].x, m_vertexPositions[i].y, m_vertexPositions[i].z);
    }

    //Loop for faces
    for(int i = 0; i < m_triArray.size(); ++i){
        file.printf(STR(f %d %d %d\n), m_triArray[i].x, m_triArray[i].y, m_triArray[i].z);
    }
    file.printf(STR(\n));
    file.commit();
}
