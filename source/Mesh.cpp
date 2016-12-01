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
    //Welder::weld(m_vertexPositions, Array<Point2>(), Array<Vector3>(), m_indexArray, Welder::Settings());
};

Mesh::Mesh(const String& filename) {
    shared_ptr<ArticulatedModel> model = ArticulatedModel::fromFile(filename);
    const Array<ArticulatedModel::Geometry*>& geometryArray(model->geometryArray());
    const Array< ArticulatedModel::Mesh*>& meshArray(model->meshArray());
    for (int i(0); i < geometryArray[0]->cpuVertexArray.vertex.size(); ++i) {
        m_vertexPositions.append(geometryArray[0]->cpuVertexArray.vertex[i].position);
    }
    m_indexArray = meshArray[0]->cpuIndexArray;
    //Welder::weld(m_vertexPositions, Array<Point2>(), Array<Vector3>(), m_indexArray, Welder::Settings());
};

Mesh::~Mesh() {};

std::shared_ptr<Mesh> Mesh::create(const Array<Vector3>& vertexPositions, const Array<Vector3int32>& triArray) {
    return std::shared_ptr<Mesh>(new Mesh(vertexPositions, triArray));
}

std::shared_ptr<Mesh> Mesh::create(const String& filename) {
    return std::shared_ptr<Mesh>(new Mesh(filename));
}

void Mesh::computeAdjacency(Array<MeshAlg::Face>& faceArray, Array<MeshAlg::Edge>& edgeArray, Array<MeshAlg::Vertex>& vertexArray) {
    //debugPrintf(STR(%d %d \n), m_indexArray.size(), m_vertexPositions.size());
    //Welder::weld(m_vertexPositions, Array<Point2>(), Array<Vector3>(), m_indexArray, Welder::Settings());
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

static int edgeLength(const MeshAlg::Edge& edge, const Array<Vector3>& vertexArray) {
    return length(vertexArray[edge.vertexIndex[1]] - vertexArray[edge.vertexIndex[0]]);
};

int Mesh::edgeLength(const MeshAlg::Edge& edge) {
    return length(m_vertexPositions[edge.vertexIndex[0]] - m_vertexPositions[edge.vertexIndex[1]]);
};

int Mesh::edgeLength(int i0, int i1) {
    return length(m_vertexPositions[i1] - m_vertexPositions[i0]);
}

static bool isManifold(const MeshAlg::Edge& edge, const Array<MeshAlg::Edge>& edges, const Array<MeshAlg::Vertex>& vertices) {
    // Check how many end points of the adjacent edges to edge coincide. 
    // A manifold will be formed if there are more than two such endpoints. 
    MeshAlg::Vertex v0(vertices[edge.vertexIndex[0]]);
    SmallArray<int, 6> adjacentEdges0(v0.edgeIndex); //indices of edges adjacent to vertex 0

    MeshAlg::Vertex v1(vertices[edge.vertexIndex[1]]);
    SmallArray<int, 6> adjacentEdges1(v1.edgeIndex); //indices of edges adjacent to vertex 1

    Array<int> endPoints0;
    Array<int> endPoints1;

    int maxSize(max(adjacentEdges0.size(), adjacentEdges1.size()));
    // Fill both endPoints arrays
    for (int i(0); i < maxSize; ++i) {
        if (i < adjacentEdges0.size()) {
            int e0 = adjacentEdges0[i];
            endPoints0.append(edges[e0 >= 0 ? e0 : ~e0].vertexIndex[1]);
        }

        if (i < adjacentEdges1.size()) {
            int e1 = adjacentEdges1[i];
            endPoints1.append(edges[e1 >= 0 ? e1 : ~e1].vertexIndex[1]);
        }
    }

    int counter(0); // Count how many common end points there are
    for (int i(0); i < endPoints0.size(); ++i) {
        if (endPoints1.contains(endPoints0[i])) {
            counter += 1;
            if (counter > 2) { // If there are more than 2 common end points, collapsing it creates a manifold.
                return true;
            }
        }
    }
    return false;
}

// Finds faces that are not collapsed with the edge
static Array<MeshAlg::Face> findRemainingFaces(const MeshAlg::Edge& edge, const Array<MeshAlg::Vertex>& vertices, const Array<MeshAlg::Face>& faces) {
    SmallArray<int, 6> faceIndices0(vertices[edge.vertexIndex[0]].faceIndex);
    SmallArray<int, 6> faceIndices1(vertices[edge.vertexIndex[1]].faceIndex);

    Array<MeshAlg::Face> remainingFaces; // Faces that contain either vertex 0 or vertex 1 of edge, but not both

    int maxSize(max(faceIndices0.size(), faceIndices1.size()));
    for (int i(0); i < maxSize; ++i) { // Fill remainingFaces
        if (i < faceIndices0.size()) {
            if (!faceIndices1.contains(faceIndices0[i])) {
                remainingFaces.append(faces[faceIndices0[i]]);
            }
        }
        if (i < faceIndices1.size()) {
            if (!faceIndices0.contains(faceIndices1[i])) {
                remainingFaces.append(faces[faceIndices1[i]]);
            }
        }
    }
    return remainingFaces;
}

static Vector3 computeCurNormal(const MeshAlg::Face& face, const Array<Vector3>& vertices) {
    Vector3 v0(vertices[face.vertexIndex[0]]);
    Vector3 v1(vertices[face.vertexIndex[1]]);
    Vector3 v2(vertices[face.vertexIndex[2]]);

    Vector3 a(v0 - v2);
    Vector3 b(v1 - v2);

    return cross(a, b).direction();
}

static bool isDegenerate(int v0, int v1, int v2) {
    return (v0 == v1) || (v0 == v2) || (v1 == v2);
}

static Vector3 computeNewNormal(const MeshAlg::Face& face, const Array<Vector3>& vertices, const MeshAlg::Edge& edgeToCollapse) {
    int toReplace(edgeToCollapse.vertexIndex[0]);   // Vertex at index 0 of edge is replaced by vertex at index 1
    int replacement(edgeToCollapse.vertexIndex[1]);

    Array<int> vertexIndices; // Vertex indices after edge collapse
    for (int i(0); i < 3; ++i) {
        int index(face.vertexIndex[i] == toReplace ? replacement : face.vertexIndex[i]);
        vertexIndices.append(index);
    }

    // If face is degenerate, then it will be discarded after collapsing, so we don't want to prevent collapsing.
    if (isDegenerate(vertexIndices[0], vertexIndices[1], vertexIndices[2])) {
        return computeCurNormal(face, vertices);
    }

    Vector3 v0(vertices[vertexIndices[0]]);
    Vector3 v1(vertices[vertexIndices[1]]);
    Vector3 v2(vertices[vertexIndices[2]]);

    Vector3 a(v0 - v2);
    Vector3 b(v1 - v2);

    return cross(a, b).direction();
}

static bool isSignOpposite(const Vector3& v1, const Vector3& v2) {
    return float(v1.dot(v2)) < 0.0f;
}

static bool normalsFlipped(const MeshAlg::Edge& edge, const Array<MeshAlg::Face>& faces, const Array<MeshAlg::Vertex>& vertices, const Array<Vector3>& vertexPositions) {
    //Array<MeshAlg::Face> facesToCheck(findRemainingFaces(edge, vertices, faces));
    //for (int i(0); i < facesToCheck.size(); ++i) {
    //    MeshAlg::Face face(facesToCheck[i]);
    //    if (isSignOpposite(computeCurNormal(face, vertexPositions), computeNewNormal(face, vertexPositions, edge))) {
    //        return true;
    //    }
    //}

    SmallArray<int, 6> faces0(vertices[edge.vertexIndex[0]].faceIndex); // Faces containing vertex 0 of edge
    SmallArray<int, 6> faces1(vertices[edge.vertexIndex[1]].faceIndex); // Faces containing vertex 1 of edge

    // Check normal signs only for non-overlapping faces, since overlapping faces are collapsed
    for (int i(0); i < faces0.size(); ++i) {
        MeshAlg::Face face(faces[faces0[i]]);
        if (isSignOpposite(computeCurNormal(face, vertexPositions), computeNewNormal(face, vertexPositions, edge))) {
            return true;
        }
    }

    for (int i(0); i < faces1.size(); ++i) {
        MeshAlg::Face face(faces[faces1[i]]);
        if (isSignOpposite(computeCurNormal(face, vertexPositions), computeNewNormal(face, vertexPositions, edge))) {
            return true;
        }
    }
    return false;
}

bool Mesh::isCollapsable(const MeshAlg::Edge& edge, const Array<MeshAlg::Face>& faces, const Array<MeshAlg::Edge>& edges, const Array<MeshAlg::Vertex>& vertices) const {
    bool isInside(!edge.boundary());
    if (isInside) { // If the edge is an inside edge
        // Edge is not collapsable if by collapsing it we get a many fold or the normals flip
        return !isManifold(edge, edges, vertices) && !normalsFlipped(edge, faces, vertices, m_vertexPositions);
    }
    return isInside;
};

static float computeAngle(const Vector3& v1, const Vector3& v2) {
    return G3D::acos(v1.dot(v2) / (length(v1)*length(v2)));
}

static bool greaterAngle(const MeshAlg::Edge& elem1, const MeshAlg::Edge& elem2, const Array<Vector3>& faceNormals, const Array<Vector3>& vertices) {
    debugAssertM(!elem1.boundary() && !elem2.boundary(), "Boundary Edge encountered!");
    float x1(computeAngle(faceNormals[elem1.faceIndex[0]], faceNormals[elem1.faceIndex[1]]));
    float x2(computeAngle(faceNormals[elem2.faceIndex[0]], faceNormals[elem2.faceIndex[1]]));
    float l1(edgeLength(elem1, vertices));
    float l2(edgeLength(elem2, vertices));
    return (fabs(x1) / sqrt(l1)) > (fabs(x2) / sqrt(l2));
}


class ComparableEdge {
public:
    MeshAlg::Edge edge;
    Vector3 leftFaceNormal;
    Vector3 rightFaceNormal;
    Vector3 v0;
    Vector3 v1;
    bool operator>(const ComparableEdge& other) const {
        const float x1(computeAngle(leftFaceNormal, rightFaceNormal));
        const float x2(computeAngle(other.leftFaceNormal, other.rightFaceNormal));
        float l1(length(v1 - v0));
        float l2(length(other.v1 - other.v0));
        return (fabs(x1) / cbrt(l1)) > (fabs(x2) / cbrt(l2));
    }

    bool operator<(const ComparableEdge& other) const {
        const float x1(computeAngle(leftFaceNormal, rightFaceNormal));
        const float x2(computeAngle(other.leftFaceNormal, other.rightFaceNormal));
        float l1(length(v1 - v0));
        float l2(length(other.v1 - other.v0));
        return (fabs(x1) / cbrt(l1)) < (fabs(x2) / cbrt(l2));
    }
};

static void remapIndices(Array<Vector3>& vertexArray, Array<int>& indexArray, int i, int j) {
    vertexArray[i] = vertexArray.last();
    vertexArray.pop();
    int old(vertexArray.size());

    for (int k(0); k < indexArray.size(); ++k) {
        if (indexArray[k] == i && j != old) {
            indexArray[k] = j;
        }
        else if (indexArray[k] == old) {
            indexArray[k] = i;
        }
    }
}

static void removeDegenerateFaces(Array<int>& indexArray) {
    for (int i(0); i < indexArray.size(); i += 3) {
        if (isDegenerate(indexArray[i], indexArray[i + 1], indexArray[i + 2])) {
            indexArray[i] = indexArray[indexArray.size() - 3];
            indexArray[i + 1] = indexArray[indexArray.size() - 2];
            indexArray[i + 2] = indexArray[indexArray.size() - 1];

            indexArray.pop();
            indexArray.pop();
            indexArray.pop();

            i -= 3;
        }
    }
}
static void collapseOneEdge(const MeshAlg::Edge& edge, Array<Vector3>& vertexArray, Array<int>& indexArray) {
    int index0(edge.vertexIndex[0]);
    int index1(edge.vertexIndex[1]);

    remapIndices(vertexArray, indexArray, index0, index1);
    removeDegenerateFaces(indexArray);
}

void Mesh::collapseEdges(int numEdges) {
    for (int x(0); x < numEdges; ++x) {
        Array<MeshAlg::Edge> edges;
        Array<MeshAlg::Face> faces;
        Array<MeshAlg::Vertex> vertices;
        Array<Vector3> faceNormals;

        computeAdjacency(faces, edges, vertices);
        computeFaceNormals(faceNormals);

        MeshAlg::Edge edge(edges[0]);
        int collapsable(0);
        //Array<ComparableEdge> compEdges;
        for (int i(0); i < edges.size(); ++i) {
            if (isCollapsable(edges[i], faces, edges, vertices)) {
                //ComparableEdge ce; 
                //ce.edge = edges[i]; 
                //ce.leftFaceNormal = faceNormals[edges[i].faceIndex[0]]; 
                //ce.rightFaceNormal = faceNormals[edges[i].faceIndex[1]];
                //ce.v0 = m_vertexPositions[edges[i].vertexIndex[0]];
                //ce.v1 = m_vertexPositions[edges[i].vertexIndex[1]]; 
                //compEdges.append(ce);
                collapsable++;
                if (greaterAngle(edges[i], edge, faceNormals, m_vertexPositions)) {
                    edge = edges[i];
                }
            }
        }
        if (collapsable < 1) {
            return;
        }
        collapseOneEdge(edge, m_vertexPositions, m_indexArray);

        //compEdges.sort(SORT_DECREASING);
        //numEdges = min(numEdges, compEdges.size()); 
        //for(int i(0); i < numEdges; ++i) { 
        //    collapseOneEdge(compEdges[i].edge, m_vertexPositions, m_indexArray);
        //}
    }
}

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
    for (int i = 0; i < newVertices.size(); ++i) {
        m_vertexNormals[i] = faceNormalArray[i / 3];
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
        file.printf("f %d %d %d\n", m_indexArray[i] + 1, m_indexArray[i + 1] + 1, m_indexArray[i + 2] + 1);
    }
    file.printf("\n");
    file.commit();
}

shared_ptr<Model> Mesh::toArticulatedModel(String name, Color3& color) const {
    String anyStr = "UniversalMaterial::Specification { lambertian = Color3(" + (String)std::to_string(color.r) + ", " + (String)std::to_string(color.g) + ", " + (String)std::to_string(color.b) + "); }";
    return toArticulatedModel(name, anyStr, 1, 1);
}

shared_ptr<Model> Mesh::toArticulatedModel(String name, String anyStr, int width, int height) const {
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
        if (m_hasFakeNormals) {
            v.normal = m_vertexNormals[i];
        }
        else {
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
        v.texCoord0 = Vector2((ix*1.0) / width, (iy*1.0) / height);


    }
    for (int i = 0; i < m_indexArray.size(); ++i) {
        indexArray.append(m_indexArray[i]);
    }

    ArticulatedModel::CleanGeometrySettings geometrySettings;
    geometrySettings.allowVertexMerging = false;
    model->cleanGeometry(geometrySettings);

    return model;
}


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