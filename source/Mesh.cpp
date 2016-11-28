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
    computeFaceNormals(m_faceNormals);
    m_vertexNormals = Array<Vector3>();
};

Mesh::Mesh(const String& filename) {
    shared_ptr<ArticulatedModel> model = ArticulatedModel::fromFile(filename);
    const Array<ArticulatedModel::Geometry*>& geometryArray(model->geometryArray());
    const Array< ArticulatedModel::Mesh*>& meshArray(model->meshArray());
    for (int i(0); i < geometryArray[0]->cpuVertexArray.vertex.size(); ++i) {
        m_vertexPositions.append(geometryArray[0]->cpuVertexArray.vertex[i].position);
    }
    m_indexArray = meshArray[0]->cpuIndexArray;
    computeFaceNormals(m_faceNormals);
};

Mesh::~Mesh() {};

std::shared_ptr<Mesh> Mesh::create(const Array<Vector3>& vertexPositions, const Array<Vector3int32>& triArray) {
    return std::shared_ptr<Mesh>(new Mesh(vertexPositions, triArray));
}

std::shared_ptr<Mesh> Mesh::create(const String& filename) {
    return std::shared_ptr<Mesh>(new Mesh(filename));
}

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

void Mesh::computeFaceNormals(Array<Vector3>& faceNormals, bool normalize) {
    Array<MeshAlg::Face> faceArray;
    computeAdjacency(faceArray);
    MeshAlg::computeFaceNormals(m_vertexPositions, faceArray, faceNormals, normalize);
};


int Mesh::edgeLength(const MeshAlg::Edge& edge) {
    return length(m_vertexPositions[edge.vertexIndex[0]] - m_vertexPositions[edge.vertexIndex[1]]);
};

int Mesh::edgeLength(int i0, int i1) {
    return length(m_vertexPositions[i1] - m_vertexPositions[i0]);
};

int computeAngle(const Vector3& v1, const Vector3& v2) {
    return G3D::acos(v1.dot(v2) / (length(v1)*length(v2)));
}

bool Mesh::greaterAngle(const MeshAlg::Edge& elem1, const MeshAlg::Edge& elem2) {
    if (elem1.boundary() || elem2.boundary()) {
        return edgeLength(elem1) < edgeLength(elem2);
    }
    int x1(computeAngle(m_faceNormals[elem1.faceIndex[0]], m_faceNormals[elem1.faceIndex[1]]));
    int x2(computeAngle(m_faceNormals[elem2.faceIndex[0]], m_faceNormals[elem2.faceIndex[1]]));
    return x1 > x2;
}

MeshAlg::Edge Mesh::toCollapse(const Array<MeshAlg::Edge>& data) {

    computeFaceNormals(m_faceNormals);

    MeshAlg::Edge toReturn(data[Random::threadCommon().integer(0, data.size() - 1)]);
    while (toReturn.boundary()) {
        toReturn = data[Random::threadCommon().integer(0, data.size() - 1)];
    }

    Vector3 fn0(m_faceNormals[toReturn.faceIndex[0]]);
    Vector3 fn1(m_faceNormals[toReturn.faceIndex[1]]);

    int maxAngle(computeAngle(fn0, fn1));

    for (int i(0); i < data.size(); ++i) {
        MeshAlg::Edge curEdge(data[i]);
        if (!curEdge.boundary()) {
            fn0 = m_faceNormals[curEdge.faceIndex[0]];
            fn1 = m_faceNormals[curEdge.faceIndex[1]];

            int curAngle(computeAngle(fn0, fn1));
            if (maxAngle < curAngle) {
                toReturn = curEdge;
            }
            else if (maxAngle == curAngle) {
                toReturn = edgeLength(toReturn) < edgeLength(curEdge) ? toReturn : curEdge;
            }
        }
    }
    return toReturn;
}




void Mesh::collapseEdges(int numEdges) {
    Array<MeshAlg::Edge> edges;
    computeAdjacency(Array<MeshAlg::Face>(), edges);

    for (int i(0); i < edges.size(); ++i) {
        if (edges[i].boundary()) {
            edges.remove(i);
        }
    }
    mergeSort(edges);
    numEdges = min(numEdges, edges.size());
    for (int x(0); x < numEdges; ++x) {
        int index0(edges[x].vertexIndex[0]);
        int index1(edges[x].vertexIndex[1]);
        
        for (int j(0); j < m_indexArray.size(); ++j) {
            if (m_indexArray[j] == index0 || m_indexArray[j] == index1) {
                m_indexArray[j] = min(index0, index1);
            }
        }

    }

    //Array<MeshAlg::Edge> edges;
    //for (int i(0); i < numEdges; ++i) {
    //    computeAdjacency(Array<MeshAlg::Face>(), edges);

    //    if (edges.size() < 1) return;
    //    numEdges = min(numEdges, edges.size());

    //    MeshAlg::Edge collapsedEdge(toCollapse(edges));

    //    int index0(collapsedEdge.vertexIndex[0]);
    //    int index1(collapsedEdge.vertexIndex[1]);

    //    for (int j(0); j < m_indexArray.size(); ++j) {
    //        if (m_indexArray[j] == index0 || m_indexArray[j] == index1) {
    //            m_indexArray[j] = min(index0, index1);
    //        }
    //    }
    //}
};

void Mesh::merge(Array<MeshAlg::Edge>& data, const Array<MeshAlg::Edge>& temp, int low, int middle, int high) {
    int ri = low;
    int ti = low;
    int di = middle;

    // While two lists are not empty, merge smaller value
    while (ti < middle && di <= high) {
        if (greaterAngle(data[di], temp[ti])) {
            data[ri++] = data[di++]; // smaller is in high data
        }
        else {
            data[ri++] = temp[ti++]; // smaller is in temp
        }
    }

    // Possibly some values left in temp array
    while (ti < middle) {
        data[ri++] = temp[ti++];
    }
    // ...or some values left in correct place in data array
};

void Mesh::mergeSortRecursive(Array<MeshAlg::Edge>& data, Array<MeshAlg::Edge>& temp, int low, int high) {
    int n = high - low + 1;
    int middle = low + n / 2;

    if (n < 2) return;
    // move lower half of data into temporary storage
    for (int i = low; i < middle; i++) {
        temp[i] = data[i];
    }

    // Sort lower half of array 
    mergeSortRecursive(temp, data, low, middle - 1);
    // sort upper half of array 
    mergeSortRecursive(data, temp, middle, high);
    // merge halves together
    merge(data, temp, low, middle, high);
};

void Mesh::mergeSort(Array<MeshAlg::Edge>& data) {
    Array<MeshAlg::Edge> temp;
    temp.resize(data.size());
    mergeSortRecursive(data, temp, 0, data.size() - 1);
};

class AngledVertex {
public:
    float angle;
    int index;
    bool operator>(const AngledVertex& other) const {
        return angle > other.angle;
    }
    bool operator<(const AngledVertex& other) const {
        return angle < other.angle;
    }
};
/*
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
        //debugPrintf(STR(Mapping face %d at original vertex %d to neww vertex %d\n), i / 3, m_indexArray[i], i);
        faceIndexMap[i / 3].set(m_indexArray[i], i);
        //debugPrintf(STR(Mapping face %d at original vertex %d to neww vertex %d\n), i / 3, m_indexArray[i + 1], i + 1);
        faceIndexMap[i / 3].set(m_indexArray[i + 1], i + 1);
        //debugPrintf(STR(Mapping face %d at original vertex %d to neww vertex %d\n), i / 3, m_indexArray[i + 2], i + 2);
        faceIndexMap[i / 3].set(m_indexArray[i + 2], i + 2);
    }

    //Iterate through edges. For each edge, find the 4 points associated with it, via indexing. Construct 2 new triangles. 
    for (int i = 0; i < edgeArray.size(); ++i) {

        //get the faceIndexes:
        int face1 = edgeArray[i].faceIndex[0];
        int face2 = edgeArray[i].faceIndex[1];

        //debugPrintf(STR(Checking face %d at original vertex %d\n), face1, edgeArray[i].vertexIndex[0]);
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
        /*
        //compute radius of sphere cap
        Vector3 f1n(faceNormalArray[indexMap[i][0]/3]);
        Vector3 f2n(faceNormalArray[indexMap[i][0]/3]);
        float mag = f1n.magnitude()*f1n.magnitude();
        float angle = acosf(dot(f1n,f2n)/mag)/2.0;
        float radius = sin(angle)*bump;
        
        //Draw polygon
        //int iOff = newIndices.size();


       //1. project them into the plane of the normal (i.e., generate an arbitrary coordinate from from the normal as the z axis; 
       //G3D has several routines for this; and then call cframe::pointToObjectSpace and only keep the xy coordinates or whatever the permutation is)
        Vector3 vNorm(vertexNormalArray[i]);

        Vector3 x(1, 0, 0);
        if (abs(dot(x, vNorm) > .9)) {
            x = Vector3(0, 1, 0);
        }
        Vector3 y = normalize(vNorm.cross(x));
        x = y.cross(vNorm);

        Array<AngledVertex> angles;
        //2. use atan2 to compute the angle in that plane to each point
        Vector3 sum(0.0, 0, 0);
        for (int j = 0; j < indexMap[i].size(); ++j) {
            sum += newVertices[indexMap[i][j]];
        }
        Vector3 C = (1.0 / indexMap[i].size())*sum;
        for (int j = 0; j < indexMap[i].size(); ++j) {
            int myIndex = indexMap[i][j];
            Vector3 a(newVertices[myIndex] - C);
            AngledVertex av;
            av.angle = atan2(dot(y, a), dot(x, a));
            av.index = myIndex;
            angles.append(av);
        }
        //3. sort the points by angle
        angles.sort(SORT_INCREASING);

        //4. connect them all in that order! 
        int point1 = angles[0].index;
        for (int j = 1; j < indexMap[i].size() - 1; ++j) {
            newIndices.append(point1, angles[j].index, angles[j + 1].index);
        }
    }

    m_vertexPositions = newVertices;
    m_indexArray = newIndices;
}*/

//Bevel Edges Without Blowing Out the planet
void Mesh::bevelEdges2(float bump) {

    //Step 1: Explode the planet
    Array<Vector3> newVertices;
    Array<int> newIndices;

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
        Vector3 v1(m_vertexPositions[m_indexArray[i]]);
        Vector3 v2(m_vertexPositions[m_indexArray[i + 1]]);
        Vector3 v3(m_vertexPositions[m_indexArray[i + 2]]);

        //Calculate vector for v1:
        //Take midpoint of v2 and v3
        Vector3 mid((v2 + v3) / 2.0);
        //Subtract from v1 to create new vector
        Vector3 vec(v1 - mid);
        //shift v1 by new vector * bump in direction
        v1 = v1 - normalize(vec)*bump;

        mid = (v1 + v3) / 2.0;
        vec = v2 - mid;
        v2 = v2 - normalize(vec)*bump;

        mid = (v2 + v1) / 2.0;
        vec = v3 - mid;
        v3 = v3 - normalize(vec)*bump;

        newVertices.append(v1, v2, v3);
        newIndices.append(i, i + 1, i + 2);

        //Save vertex mapping
        indexMap[m_indexArray[i]].append(i);
        indexMap[m_indexArray[i + 1]].append(i + 1);
        indexMap[m_indexArray[i + 2]].append(i + 2);

        //Face index is i%3
        //debugPrintf(STR(Mapping face %d at original vertex %d to neww vertex %d\n), i / 3, m_indexArray[i], i);
        faceIndexMap[i / 3].set(m_indexArray[i], i);
        //debugPrintf(STR(Mapping face %d at original vertex %d to neww vertex %d\n), i / 3, m_indexArray[i + 1], i + 1);
        faceIndexMap[i / 3].set(m_indexArray[i + 1], i + 1);
        //debugPrintf(STR(Mapping face %d at original vertex %d to neww vertex %d\n), i / 3, m_indexArray[i + 2], i + 2);
        faceIndexMap[i / 3].set(m_indexArray[i + 2], i + 2);
    }

    //Iterate through edges. For each edge, find the 4 points associated with it, via indexing. Construct 2 new triangles. 
    for (int i = 0; i < edgeArray.size(); ++i) {

        //get the faceIndexes:
        int face1 = edgeArray[i].faceIndex[0];
        int face2 = edgeArray[i].faceIndex[1];

        //debugPrintf(STR(Checking face %d at original vertex %d\n), face1, edgeArray[i].vertexIndex[0]);
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

        //Draw polygon
        //int iOff = newIndices.size();


       //1. project them into the plane of the normal (i.e., generate an arbitrary coordinate from from the normal as the z axis; 
       //G3D has several routines for this; and then call cframe::pointToObjectSpace and only keep the xy coordinates or whatever the permutation is)
        Vector3 vNorm(vertexNormalArray[i]);

        Vector3 x(1, 0, 0);
        if (abs(dot(x, vNorm) > .9)) {
            x = Vector3(0, 1, 0);
        }
        Vector3 y = normalize(vNorm.cross(x));
        x = y.cross(vNorm);

        Array<AngledVertex> angles;
        //2. use atan2 to compute the angle in that plane to each point
        Vector3 sum(0.0, 0, 0);
        for (int j = 0; j < indexMap[i].size(); ++j) {
            sum += newVertices[indexMap[i][j]];
        }
        Vector3 C = (1.0 / indexMap[i].size())*sum;
        for (int j = 0; j < indexMap[i].size(); ++j) {
            int myIndex = indexMap[i][j];
            Vector3 a(newVertices[myIndex] - C);
            AngledVertex av;
            av.angle = atan2(dot(y, a), dot(x, a));
            av.index = myIndex;
            angles.append(av);
        }
        //3. sort the points by angle
        angles.sort(SORT_INCREASING);

        //4. connect them all in that order! 
        int point1 = angles[0].index;
        for (int j = 1; j < indexMap[i].size() - 1; ++j) {
            newIndices.append(point1, angles[j].index, angles[j + 1].index);
        }
    }

    //Fill the vertex normal array for the new vertices
    //Since they're in order, each vertex index divided by 3 should give the appropriate face index
    m_vertexNormals.resize(newVertices.size());
    for (int i = 0; i < newVertices.size(); ++i){
        m_vertexNormals[i] = faceNormalArray[i/3];
    }
    m_hasFakeNormals = true;

    m_vertexPositions = newVertices;
    m_indexArray = newIndices;
}

void Mesh::toObj(String filename) {
    TextOutput file(filename + ".obj");
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

void Mesh::merge(SmallArray<float, 6>& data, SmallArray<float, 6>& temp, int low, int middle, int high, SmallArray<int, 6>& along, SmallArray<int, 6>& temp2) {
    int ri = low;
    int ti = low;
    int di = middle;

    // While two lists are not empty, merge smaller value
    while (ti < middle && di <= high) {
        if (data[di] < temp[ti]) {
            along[ri] = along[di];
            data[ri++] = data[di++]; // smaller is in high data
        }
        else {
            along[ri] = along[ti];
            data[ri++] = temp[ti++]; // smaller is in temp
        }
    }

    // Possibly some values left in temp array
    while (ti < middle) {
        along[ri] = along[ti];
        data[ri++] = temp[ti++];
    }
    // ...or some values left in correct place in data array
}

void Mesh::mergeSortRecursive(SmallArray<float, 6>& data, SmallArray<float, 6>& temp, int low, int high, SmallArray<int, 6>& along, SmallArray<int, 6>& temp2) {
    int n = high - low + 1;
    int middle = low + n / 2;

    if (n < 2) return;
    // move lower half of data into temporary storage
    for (int i = low; i < middle; i++) {
        temp[i] = data[i];
        temp2[i] = along[i];
    }

    // Sort lower half of array 
    mergeSortRecursive(temp, data, low, middle - 1, along, temp2);
    // sort upper half of array 
    mergeSortRecursive(data, temp, middle, high, along, temp2);
    // merge halves together
    merge(data, temp, low, middle, high, along, temp2);
}

void Mesh::mergeSort(SmallArray<float, 6>& data, SmallArray<int, 6>& along) {
    SmallArray<float, 6> newArray;
    newArray.resize(data.size());
    SmallArray<int, 6> newArray2;
    newArray2.resize(data.size());
    mergeSortRecursive(data, newArray, 0, data.size() - 1, along, newArray2);
}

shared_ptr<Model> Mesh::toArticulatedModel(String name, Color3& color) {
    String anyStr = "UniversalMaterial::Specification { lambertian = Color3(" + (String)std::to_string(color.r) + ", " + (String)std::to_string(color.g) + ", " + (String)std::to_string(color.b) + "); }";
    return toArticulatedModel(name, anyStr, 1, 1);
}

shared_ptr<Model> Mesh::toArticulatedModel(String name, String anyStr, int width, int height) {
    const shared_ptr<ArticulatedModel>& model = ArticulatedModel::createEmpty(name);
    ArticulatedModel::Part*     part = model->addPart("root");
    ArticulatedModel::Geometry* geometry = model->addGeometry("geom");
    ArticulatedModel::Mesh*     mesh = model->addMesh("mesh", part, geometry);

    //Any any = Any::parse("UniversalMaterial::Specification { lambertian = Color3(" + String(color.r) + ", " + String(color.g) + ", " + String(color.b) + "); }");
    //String test = (String) std::to_string(color.r);
    //String anyStr = "UniversalMaterial::Specification { lambertian = Color3(" + (String)std::to_string(color.r) + ", " + (String)std::to_string(color.g) + ", " + (String)std::to_string(color.b) + "); }";
    Any any = Any::parse(anyStr);
    mesh->material = UniversalMaterial::create(any);
    /*        PARSE_ANY(
            UniversalMaterial::Specification {
                /*lambertian = Texture::Specification {
                    filename = "image/checker-32x32-1024x1024.png";
                    // Orange
                    encoding = Color3(1.0, 0.7, 0.15);
                };

                glossy     = Color4(Color3(0.01), 0.2);
                lambertian = Color3(0,1,.2);
            }));*/

    Array<CPUVertexArray::Vertex>& vertexArray = geometry->cpuVertexArray.vertex;
    Array<int>& indexArray = mesh->cpuIndexArray;

    for (int i = 0; i < m_vertexPositions.size(); ++i) {
        CPUVertexArray::Vertex& v = vertexArray.next();
        v.position = m_vertexPositions[i];

        //fix
        //v.texCoord0 = Vector2(0,0);

        // Set to NaN to trigger automatic vertex normal and tangent computation
        if (m_hasFakeNormals){
            v.normal = m_vertexNormals[i];
        }
        else{
            v.normal = Vector3::nan();
        }
        v.tangent = Vector4::nan();

        Vector3 vertex = m_vertexPositions[i];

        Vector3 d = (vertex - Vector3(0, 0, 0)).unit();

        float nx = width * (0.5f + atanf(d.z / d.x) / (2.0f*pif()));
        float ny = height * (0.5f - asinf(d.y) * 1 / pif());

        int ix = (int)abs((int)nx % width);
        int iy = (int)abs((int)ny % height);
        //v.texCoord0 = Vector2(nx, ny);
        //v.texCoord0 = Vector2(ix, iy);
        v.texCoord0 = Vector2((ix*1.0)/width, (iy*1.0)/height);


    }
    for (int i = 0; i < m_indexArray.size(); ++i) {
        indexArray.append(m_indexArray[i]);
    }

    ArticulatedModel::CleanGeometrySettings geometrySettings;
    geometrySettings.allowVertexMerging = false;
    model->cleanGeometry(geometrySettings);

    return model;
}

