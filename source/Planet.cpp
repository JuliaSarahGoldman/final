#include "Planet.h"
#include <unordered_map>

void Planet::writeSphere(String filename, shared_ptr<Array<Vector3>>& vertices, shared_ptr<Array<Vector3int32>>& faces) {
    makeIcohedron(5.0f, vertices, faces);

    int numVert = vertices->size();
    int numFace = faces->size();

    String sphere = "OFF\n";
    sphere += (String) std::to_string(numVert) + " " + (String) std::to_string(numFace) + " 0\n\n";

    for(int i = 0; i < numVert; ++i) {
        Vector3 vertex = vertices->operator[](i);
        sphere += (String) std::to_string(vertex[0]) + " " + (String) std::to_string(vertex[1]) + " " + (String) std::to_string(vertex[2]) + "\n";
    }

    for(int i = 0; i < numFace; ++i) {
        Vector3 face = faces->operator[](i);
        sphere += "3 " + (String) std::to_string(face[0]) + " " + (String) std::to_string(face[1]) + " " + (String) std::to_string(face[2]) + "\n";
    }

    TextOutput output ("model/" + filename + ".off");
    output.writeSymbol(sphere);
    output.commit(true);
    G3D::ArticulatedModel::clearCache();
}

void Planet::makeIcohedron(float radius, shared_ptr<Array<Vector3>>& vertices, shared_ptr<Array<Vector3int32>>& faces) {
    float t = (1.0 + sqrt(5.0)) / 2.0;
    float s = sqrt(1 + square(t));

    vertices->append(radius * Vector3(t, 1, 0) / s);
    vertices->append(radius * Vector3(-t, 1, 0) / s);
    vertices->append(radius * Vector3(t, -1, 0) / s);
    vertices->append(radius * Vector3(-t, -1, 0) / s);
    vertices->append(radius * Vector3(1, 0, t) / s);
    vertices->append(radius * Vector3(1, 0, -t) / s);
    vertices->append(radius * Vector3(-1, 0, t) / s);
    vertices->append(radius * Vector3(-1, 0, -t) / s);
    vertices->append(radius * Vector3(0, t, 1) / s);
    vertices->append(radius * Vector3(0, -t, 1) / s);
    vertices->append(radius * Vector3(0, t, -1) / s);
    vertices->append(radius * Vector3(0, -t, -1) / s);

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
}

void Planet::getMiddle(float radius, Vector3& v1, Vector3& v2, Vector3& newVector){
    float t = (1+sqrt(5.0f))/2;
    newVector = (v2-v1)*0.5+v1;
    newVector.unit();
    newVector *= (sqrt(t * t + 1) * radius)/t;
}

void Planet::subdivideIcoHedron(float radius, shared_ptr<Array<Vector3>>& vertices, shared_ptr<Array<Vector3int32>>& faces){
    float t = (1+sqrt(5))/2;

    std::unordered_map<size_t, int> vertexPositions;
    for (int i(0); i < vertices->length(); ++i) {
        vertexPositions[vertices->operator[](i).hashCode()] = i;
    }

    Vector3 newVec1;
    Vector3 newVec2;
    Vector3 newVec3;

    for (int i(0); i < faces->length(); i += 3){
        
        Vector3int32 face = faces->operator[](i);

        Vector3 vert1 = vertices->operator[](face.x);
        Vector3 vert2 = vertices->operator[](face.y);
        Vector3 vert3 = vertices->operator[](face.z);

        try{
            int pos1 = vertexPositions.at(vert1.hashCode());
            int pos2 = vertexPositions.at(vert2.hashCode());
            int pos3 = vertexPositions.at(vert3.hashCode());
        } catch(...){
            printf("Undefined face in face array at position %d", i);
            return;
        }
        //Find the midpoints
        getMiddle(radius, vertices->operator[](faces->operator[](i).x), vertices->operator[](i + 1), newVec1);
        getMiddle(radius, vertices->operator[](i + 1), vertices->operator[](i + 2), newVec1);
        getMiddle(radius, vertices->operator[](i + 2), vertices->operator[](i + 3), newVec1);

        vertices->append(newVec1, newVec2, newVec3);

        //add the new faces
        Vector3 temp;
        faces->append(Vector3int32(0,0,0));
    }
}