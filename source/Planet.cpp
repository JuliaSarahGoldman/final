#include "Planet.h"

void Planet::writeSphere(String filename, float radius, int depths, shared_ptr<Array<Vector3>>& vertices, shared_ptr<Array<Vector3int32>>& faces) {
    
    makeIcohedron(radius, vertices, faces);

    for(int i(0); i < depths; ++i){
        subdivideIcoHedron(radius, vertices, faces);
    }

    //const G3D::Welder::Settings settings;
    Array<Vector3> normals;
    Array<Vector2> texture;
    Array<int> indices;
    Array<Vector3> verts;
    
    for (int i(0); i < vertices->length(); ++i) {
        verts.append(vertices->operator[](i));
    }

    for(int i(0); i < faces->length(); ++i){
        Vector3int32 face = faces->operator[](i);
        indices.append(face.x, face.y, face.z);
    }

    Welder::weld(verts, texture, normals, indices, G3D::Welder::Settings());

    faces = std::make_shared<Array<Vector3int32>>();
    for(int i(0); i < indices.size()-2; i += 3) {
        faces->append(Vector3int32(indices[i], indices[i+1], indices[i+2]));
    }

    Noise noise;
    for(int i(0); i < verts.size(); ++i) {
        Vector3 vertex = verts[i];
        float n = noise.sampleFloat(vertex.x, vertex.y, vertex.z, 1) * radius * 1000;
        debugPrintf("%f\n", n);
        verts[i] += normals[i] * n;
    }
    vertices = std::make_shared<Array<Vector3>>(verts);

    /*int numVert = vertices->size();
    int numFace = faces->size();

    TextOutput output("model/" + filename + ".obj");

    for (int i = 0; i < numVert; ++i) {
        Vector3 vertex = vertices->operator[](i);
        output.printf("v %f %f %f\n", vertex[0], vertex[1], vertex[2]);
    }

    for (int i = 0; i < numFace; ++i) {
        Vector3int32 face = faces->operator[](i);
        output.printf("f %d %d %d\n", face[0]+1, face[1]+1, face[2]+1);
    }
    output.commit(true);
    G3D::ArticulatedModel::clearCache();*/
}

void Planet::makeIcohedron(float radius, shared_ptr<Array<Vector3>>& vertices, shared_ptr<Array<Vector3int32>>& faces) {
    float t = (1.0f + sqrt(5.0f)) / 2.0f;

    vertices->append(radius * Vector3(t, 1, 0));
    vertices->append(radius * Vector3(-t, 1, 0));
    vertices->append(radius * Vector3(t, -1, 0));
    vertices->append(radius * Vector3(-t, -1, 0));
    vertices->append(radius * Vector3(1, 0, t));
    vertices->append(radius * Vector3(1, 0, -t));
    vertices->append(radius * Vector3(-1, 0, t));
    vertices->append(radius * Vector3(-1, 0, -t));
    vertices->append(radius * Vector3(0, t, 1));
    vertices->append(radius * Vector3(0, -t, 1));
    vertices->append(radius * Vector3(0, t, -1));
    vertices->append(radius * Vector3(0, -t, -1));

    faces->append(Vector3int32(0, 8, 4));
    faces->append(Vector3int32(1, 10, 7));
    faces->append(Vector3int32(2, 9, 11));
    faces->append(Vector3int32(7, 3, 1));
    faces->append(Vector3int32(0, 5, 10));
    faces->append(Vector3int32(3, 9, 6));
    faces->append(Vector3int32(3, 11, 9));
    faces->append(Vector3int32(8, 6, 4));
    faces->append(Vector3int32(2, 4, 9));
    faces->append(Vector3int32(3, 7, 11));
    faces->append(Vector3int32(4, 2, 0));
    faces->append(Vector3int32(9, 4, 6));
    faces->append(Vector3int32(2, 11, 5));
    faces->append(Vector3int32(0, 10, 8));
    faces->append(Vector3int32(5, 0, 2));
    faces->append(Vector3int32(10, 5, 7));
    faces->append(Vector3int32(1, 6, 8));
    faces->append(Vector3int32(1, 8, 10));
    faces->append(Vector3int32(6, 1, 3));
    faces->append(Vector3int32(11, 7, 5));

    /*for(int i(0); i < 6; ++i){
        subdivideIcoHedron(radius, vertices, faces);
    }
    //const G3D::Welder::Settings settings;
    Array<Vector3> normals;
    Array<Vector2> texture;
    Array<int> indices;
    Array<Vector3> verts;
    
    for (int i(0); i < vertices->length(); ++i) {
        verts.append(vertices->operator[](i));
    }

    for(int i(0); i < faces->length(); ++i){
        Vector3int32 face = faces->operator[](i);
        indices.append(face.x, face.y, face.z);
    }

    Welder::weld(verts, texture, normals, indices, G3D::Welder::Settings());

    faces = std::make_shared<Array<Vector3int32>>();
    for(int i(0); i < indices.size()-3; i += 3) {
        faces->append(Vector3int32(indices[i], indices[i+1], indices[i+2]));
    }

    for(int i(0); i < verts.size(); ++i) {
        verts[i] += normals[i] * Random::threadCommon().uniform(1.0, 1.2);
    }
    vertices = std::make_shared<Array<Vector3>>(verts);*/
}

void Planet::getMiddle(float radius, Vector3& v1, Vector3& v2, Vector3& newVector) {
    float t = (1.0f + sqrt(5.0f)) / 2.0f;
    newVector = (v2 + v1)/2.0f;
    newVector = newVector.unit();
    newVector *= sqrt(t*t+1)*radius;
}

void Planet::subdivideIcoHedron(float radius, shared_ptr<Array<Vector3>>& vertices, shared_ptr<Array<Vector3int32>>& faces) {
    float t = (1.0f + sqrt(5.0f)) / 2.0f;

    Vector3 newVec1;
    Vector3 newVec2;
    Vector3 newVec3;
    
    int numFaces = faces->length();
    
    shared_ptr<Array<Vector3int32>> newFaces = std::make_shared<Array<Vector3int32>>();

    for (int i(0); i < numFaces; ++i) {
    
        int numVertices = vertices->length();
        Vector3int32 face = faces->operator[](i);

        Vector3 vert1 = vertices->operator[](face.x);
        Vector3 vert2 = vertices->operator[](face.y);
        Vector3 vert3 = vertices->operator[](face.z);

        //Find the midpoints
        getMiddle(radius, vert1, vert2, newVec1);
        getMiddle(radius, vert2, vert3, newVec2);
        getMiddle(radius, vert3, vert1, newVec3);

        vertices->append(newVec1, newVec2, newVec3);

        int posVec1 = numVertices;
        int posVec2 = numVertices + 1;
        int posVec3 = numVertices + 2;

        //add the new faces
        newFaces->append(Vector3int32(posVec3, face.x, posVec1));
        newFaces->append(Vector3int32(posVec1, face.y, posVec2));
        newFaces->append(Vector3int32(posVec2, face.z, posVec3));
        newFaces->append(Vector3int32(posVec1, posVec2, posVec3));
    }
    faces = newFaces;
}

/*
    
    Table<Vector3, int> vertexPositions;
    for(int i(0); i < vertices->length(); ++i) {
        vertexPositions.set(vertices->operator[](i), i);
    }

            if(!vertexPositions.get(newVec1, posVec1)) {
            vertices->append(newVec1);
            vertexPositions.set(newVec1, numVertices);
            posVec1 = numVertices;
        }
        if(!vertexPositions.get(newVec2, posVec2)) {
            vertices->append(newVec2);
            vertexPositions.set(newVec2, numVertices + 1);
            posVec2 = numVertices + 1;
        }
        if(!vertexPositions.get(newVec3, posVec3)) {
            vertices->append(newVec3);
            vertexPositions.set(newVec3, numVertices + 2);
            posVec3 = numVertices + 2;
        }
*/