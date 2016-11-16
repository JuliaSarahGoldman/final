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

int Mesh::edgeLength(int i0, int i1) {
    return length(m_vertexPositions[i1] - m_vertexPositions[i0]);
};


Array<int> randomIntList(int min, int max) {
    Array<int> toReturn;
    int iters(max);
    while (iters >= min) {
        int r = Random::threadCommon().integer(min, max);
        if (!toReturn.contains(r)) {
            toReturn.append(r);
            --iters;
        }
    }
    return toReturn;
}

MeshAlg::Edge Mesh::findMinEdge(const Array<MeshAlg::Edge>& data) {
    MeshAlg::Edge minEdge(data[Random::threadCommon().integer(0, data.size() - 1)]);
    Array<int> ind(randomIntList(0, data.size() - 1));

    for (int i(0); i < ind.size(); ++i) {
        if (edgeLength(minEdge) < edgeLength(data[ind[i]])) {
            minEdge = data[ind[i]];
        }
    }
    return minEdge;
}

int Mesh::findMinEdge(int minIndex, int maxIndex) {
    int min(edgeLength(m_indexArray[minIndex], m_indexArray[minIndex+1])); 
    int r = minIndex;

    for(int i(minIndex); i < maxIndex-1; ++i) { 
        int cur(edgeLength(m_indexArray[minIndex], m_indexArray[minIndex+1]));
        if(min > cur) { 
            min = cur; 
            r = i;
        }
    }

    return r;
}

Array<Array<int>> Mesh::toCollapse(int regionSize) {
    Array<Array<int>> toReturn;
    for (int i(0); i < m_indexArray.size() - 1; i += regionSize) {
        int j (findMinEdge(i, i+regionSize-1));
        toReturn.append(Array<int>(m_indexArray[j], m_indexArray[j + 1]));
    }
    return toReturn;
}

void Mesh::collapseEdges(int regionSize) {
    Array<Array<int>> collapseList(toCollapse(regionSize));
    for (int j(0); j < collapseList.size(); ++j) {
        int index0(collapseList[j][0]);
        int index1(collapseList[j][1]);

        Vector3 average = (m_vertexPositions[index1] +  m_vertexPositions[index0])/2; 
        m_vertexPositions[index0] = average; 
        m_vertexPositions[index1] = average;
    }
    Welder::weld(m_vertexPositions, Array<Vector2>(), Array<Vector3>(), m_indexArray, G3D::Welder::Settings());
};

void Mesh::bevelEdges(float bump) {
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

    //Map from vertex index to all new indices
    Array<SmallArray<int, 6>> indexMap;

    indexMap.resize(m_vertexPositions.size());

    //Map from face and old index to new index
    Array<Table<int, int>> faceIndexMap;

    faceIndexMap.resize(m_indexArray.size());

    //Iterate through the faces, creating new vertices and a new face using those vertices for each one
    for (int i = 0; i < m_indexArray.size(); i += 3) {
        Vector3 normal(faceNormalArray[i / 3]);
        Vector3 v1(m_vertexPositions[m_indexArray[i]] + bump*normal);
        Vector3 v2(m_vertexPositions[m_indexArray[i + 1]] + bump*normal);
        Vector3 v3(m_vertexPositions[m_indexArray[i + 2]] + bump*normal);
        newVertices.append(v1, v2, v3);
        newIndices.append(i, i + 1, i + 2);

        //Save vertex mapping
        indexMap[m_indexArray[i]].append(i);
        indexMap[m_indexArray[i + 1]].append(i + 1);
        indexMap[m_indexArray[i + 2]].append(i + 2);

        //Face index is i%3
        debugPrintf(STR(Mapping face %d at original vertex %d to neww vertex %d\n), i / 3, m_indexArray[i], i);
        faceIndexMap[i / 3].set(m_indexArray[i], i);
        debugPrintf(STR(Mapping face %d at original vertex %d to neww vertex %d\n), i / 3, m_indexArray[i + 1], i + 1);
        faceIndexMap[i / 3].set(m_indexArray[i + 1], i + 1);
        debugPrintf(STR(Mapping face %d at original vertex %d to neww vertex %d\n), i / 3, m_indexArray[i + 2], i + 2);
        faceIndexMap[i / 3].set(m_indexArray[i + 2], i + 2);
    }

    //Iterate through edges. For each edge, find the 4 points associated with it, via indexing. Construct 2 new triangles. 
    for (int i = 0; i < edgeArray.size(); ++i) {

        //get the faceIndexes:
        int face1 = edgeArray[i].faceIndex[0];
        int face2 = edgeArray[i].faceIndex[1];

        debugPrintf(STR(Checking face %d at original vertex %d\n), face1, edgeArray[i].vertexIndex[0]);
        //Problem- we have the vertex index, not the index index. Oh, vertexIndex is the index index.... That's hwo they're labeled
        int v1 = faceIndexMap[face1][edgeArray[i].vertexIndex[0]];
        int v2 = faceIndexMap[face1][edgeArray[i].vertexIndex[1]];
        int v3 = faceIndexMap[face2][edgeArray[i].vertexIndex[0]];
        int v4 = faceIndexMap[face2][edgeArray[i].vertexIndex[1]];

        newIndices.append(v3, v2, v1);
        newIndices.append(v4, v2, v3);

    }


    //Now iterate through the vertices
     //Assume that we're working with a topologicalically closed shape, and every vertex is in at least 3 triangles.
    for (int i = 0; i < vertexArray.size(); ++i) {
        //use indexMap[i]
        //compute radius of sphere cap
        Vector3 f1n(faceNormalArray[indexMap[i][0] / 3]);
        Vector3 f2n(faceNormalArray[indexMap[i][0] / 3]);
        float mag = f1n.magnitude()*f1n.magnitude();
        float angle = acosf(dot(f1n, f2n) / mag) / 2.0;
        float radius = sin(angle)*bump;

        //Draw polygon
        int iOff = newIndices.size();

        for (int j = 0; j < indexMap[i].size(); ++j) {
            //newVertices.append();
            newIndices.append(iOff, iOff + 1, iOff + 2);
        }
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
    for (int i = 0; i < m_indexArray.size(); i += 3) {
        file.printf(STR(f %d %d %d\n), m_indexArray[i] + 1, m_indexArray[i + 1] + 1, m_indexArray[i + 2] + 1);
    }
    file.printf(STR(\n));
    file.commit();
}

