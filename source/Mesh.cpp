/** \file Mesh.cpp */
#include "Mesh.h"

Mesh::Mesh(const Array<Vector3>& vertexPositions, const Array<Vector3int32>& triArray) {
    MeshBuilder builder;
    for (int i = 0; i < triArray.size(); ++i) {
        Vector3int32 v = triArray[i];
        m_indexArray.append(v.x, v.y, v.z);
    }
    m_vertexPositions = vertexPositions;
    m_triArray = triArray;
};

Mesh::~Mesh() {};

void Mesh::addVertex(const Vector3& vertex) {
    m_vertexPositions.append(vertex);
};

void Mesh::addVertex(const Array<Vector3>& vertexList) {
    m_vertexPositions.append(vertexList);
};

void Mesh::addIndex(int index) {
    m_indexArray.append(index);
};

void Mesh::addIndex(const Array<int>& indexList) {
    m_indexArray.append(indexList);
};

void Mesh::addVertex(const Vector3& vertex, int index) {
    addVertex(vertex);
    addIndex(index);
};

void Mesh::addVertex(const Array<Vector3>& vertexList, const Array<int>& indexList) {
    addVertex(vertexList);
    addIndex(indexList);
};

void Mesh::computeAdjacency(Array<MeshAlg::Face>& faceArray, Array<MeshAlg::Edge>& edgeArray, Array<MeshAlg::Vertex>& vertexArray) {
    //debugPrintf(STR(%d %d \n), m_indexArray.size(), m_vertexPositions.size());
    MeshAlg::computeAdjacency(m_vertexPositions, m_indexArray, faceArray, edgeArray, vertexArray);
};

void Mesh::computeNormals(const Array<MeshAlg::Face>& faceArray, const Array<MeshAlg::Edge>& edgeArray, const Array<MeshAlg::Vertex>& vertexArray,
    Array<Vector3>& vertexNormalArray, Array<Vector3>& faceNormalArray) {
    MeshAlg::computeNormals(m_vertexPositions, faceArray, vertexArray, vertexNormalArray, faceNormalArray);
};

int Mesh::edgeLength(const MeshAlg::Edge& edge) {
    return length(m_vertexPositions[edge.vertexIndex[0]] - m_vertexPositions[edge.vertexIndex[1]]);
};


Array<int> randomIntList(int min, int max) { 
    Array<int> toReturn;
    int iters(max); 
    while(iters >= min) { 
        int r = Random::threadCommon().integer(min,max); 
        if(!toReturn.contains(r)) {
            toReturn.append(r);
            --iters;
        }    
    }
    return toReturn;
}

MeshAlg::Edge Mesh::findMinEdge(const Array<MeshAlg::Edge>& data) {
    MeshAlg::Edge minEdge(data[Random::threadCommon().integer(0,data.size()-1)]);
    Array<int> ind(randomIntList(0, data.size()-1));
  
    for (int i(0); i < ind.size(); ++i) {
        if (edgeLength(minEdge) < edgeLength(data[ind[i]])) {
            minEdge = data[ind[i]];
        }
    }
    return minEdge;
}

void Mesh::collapseEdges(int numCollapsed) {
    Array<MeshAlg::Edge> edges;
    for (int i(1); i < numCollapsed; ++i) {
        computeAdjacency(Array<MeshAlg::Face>(), edges, Array<MeshAlg::Vertex>());

        MeshAlg::Edge minEdge(findMinEdge(edges));
        int newIndex(minEdge.vertexIndex[0]);
        int toCollapse(minEdge.vertexIndex[1]);

        for (int j(0); j < m_indexArray.size(); ++j) {
            if (m_indexArray[j] == toCollapse) {
                m_indexArray[j] = newIndex;
            }
        }
    }
};

void Mesh::bevelEdges(float bump){
    //Step 1: Explode the planet
    Array<Vector3> newVertices;
    Array<int> newIndices;

    //We need normals.
    Array<MeshAlg::Face> faceArray;
    Array<MeshAlg::Edge> edgeArray;
    Array<MeshAlg::Vertex> vertexArray;
    computeAdjacency(faceArray, edgeArray, vertexArray);

    Array<Vector3> vertexNormalArray;
    Array<Vector3> faceNormalArray;
    computeNormals(faceArray, edgeArray, vertexArray, vertexNormalArray, faceNormalArray);

    //Iterate through the face, creating new vertices and a new face using those vertices for each one
    for (int i = 0; i < m_indexArray.size(); i+=3) {
       Vector3 normal( faceNormalArray[i/3]);
       Vector3 v1(m_vertexPositions[m_indexArray[i]] + bump*normal);
       Vector3 v2(m_vertexPositions[m_indexArray[i+1]] + bump*normal);
       Vector3 v3(m_vertexPositions[m_indexArray[i+2]] + bump*normal);
       newVertices.append(v1, v2, v3);
       newIndices.append(i, i+1, i+2);
    }

    //Temporary code for testing
    m_vertexPositions = newVertices;
    m_indexArray = newIndices;
}

void Mesh::toObj(String filename) {
    TextOutput file(filename);
    //file.printf("g name\n");
        //loop to make vertices
    //debugPrintf(STR(%d\n), sizeof(m_vertexArray));
    for (int i = 0; i < m_vertexPositions.size(); ++i) {
        file.printf(STR(v %f %f %f\n), m_vertexPositions[i].x, m_vertexPositions[i].y, m_vertexPositions[i].z);
    }

    //Loop for faces
    //using m_triArray
    /*
    for (int i = 0; i < m_triArray.size(); ++i) {
        file.printf(STR(f %d %d %d\n), m_triArray[i].x, m_triArray[i].y, m_triArray[i].z);
    }*/
    //using m_indexArray
    for (int i = 0; i < m_indexArray.size(); i+=3) {
        file.printf(STR(f %d %d %d\n), m_indexArray[i]+1, m_indexArray[i+1]+1, m_indexArray[i+2]+1);
    }
    file.printf(STR(\n));
    file.commit();
}

