/** \file SimpleMesh.cpp */
#include "SimpleMesh.h"

SimpleMesh::SimpleMesh(Array<Vector3> vertices, Array<Vector3int32> triangles){
    //not how you set stuff in constructor- fix
    m_vertexArray = vertices;
    m_triArray = triangles;
}
void SimpleMesh::toObj(String filename){
    TextOutput file(filename);
    file.printf("g name\n");
        //loop to make vertices
    //debugPrintf(STR(%d\n), sizeof(m_vertexArray));
    for(int i = 0; i < m_vertexArray.size(); ++i){
        file.printf(STR(v %f %f %f\n), m_vertexArray[i].x, m_vertexArray[i].y, m_vertexArray[i].z);
    }

    //Loop for faces
    for(int i = 0; i < m_triArray.size(); ++i){
        file.printf(STR(f %d %d %d\n), m_triArray[i].x, m_triArray[i].y, m_triArray[i].z);
    }
    file.printf(STR(\n));
    file.commit();
}