#include "Planet.h"
#include <unordered_map>
#include <G3D/G3DAll.h>

void Planet::writeSphere(String filename, shared_ptr<Array<Vector3>>& vertices, shared_ptr<Array<Vector3>>& faces) {
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

void Planet::makeIcohedron(float radius, shared_ptr<Array<Vector3>>& vertices, shared_ptr<Array<Vector3>>& faces) {
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

    faces->append(Vector3(0, 8, 4));
    faces->append(Vector3(1, 10, 7));
    faces->append(Vector3(2, 9, 11));
    faces->append(Vector3(7, 3, 1));
    faces->append(Vector3(0, 5, 10));
    faces->append(Vector3(3, 9, 6));
    faces->append(Vector3(3, 11, 9));
    faces->append(Vector3(8, 6, 4));
    faces->append(Vector3(2, 4, 9));
    faces->append(Vector3(3, 7, 11));
    faces->append(Vector3(4, 2, 0));
    faces->append(Vector3(9, 4, 6));
    faces->append(Vector3(2, 11, 5));
    faces->append(Vector3(0, 10, 8));
    faces->append(Vector3(5, 0, 2));
    faces->append(Vector3(10, 5, 7));
    faces->append(Vector3(1, 6, 8));
    faces->append(Vector3(1, 8, 10));
    faces->append(Vector3(6, 1, 3));
    faces->append(Vector3(11, 7, 5));
}

void Planet::getMiddle(float radius, Vector3& v1, Vector3& v2, Vector3& newVector){
    float t = (1+sqrt(5.0f))/2;

    newVector = (v2-v1)*0.5+v1;

    newVector.unit();
    
    newVector *= (sqrt(t * t + 1) * radius)/t;
}

void Planet::subdivideIcoHedron(float radius, shared_ptr<Array<Vector3>>& vertices, shared_ptr<Array<Vector3>>& faces){
    float t = (1+sqrt(5))/2;
    shared_ptr<Array<Vector3>>& vertices = std::make_shared<Array<Vector3>>();

    std::unordered_map<Vector3, int> vertexPositions;
    for (int i(0); i < vertices->length; ++i) {
        vertexPositions[vertices->operator[](i)] = i;
    }

    Vector3 newVec1;
    Vector3 newVec2;
    Vector3 newVec3;

    for (int i(0); i < faces->length(); i += 3){
        
        //Find the midpoints
        getMiddle(radius, vertices->operator[](i), vertices->operator[](i + 1), newVec1);
        getMiddle(radius, vertices->operator[](i + 1), vertices->operator[](i + 2), newVec1);
        getMiddle(radius, vertices->operator[](i + 2), vertices->operator[](i + 3), newVec1);

        vertices->append(newVec1, newVec2, newVec3);

        //add the new faces
        Vector3 temp;
        faces->append(Vector3(0.0f,0.0f,0.0f));
    }
}