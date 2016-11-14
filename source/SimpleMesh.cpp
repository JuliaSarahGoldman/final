/** \file SimpleMesh.cpp */
#include "SimpleMesh.h"

SimpleMesh::SimpleMesh(Array<Vector3> vertices, Array<Vector3int32> triangles){
    //not how you set stuff in constructor- fix
    m_vertexArray = vertices;
    m_triArray = triangles;
}

void SimpleMesh::toObj(String filename){

    int numVert = m_vertexArray.size();
    int numFace = m_triArray.size();

    TextOutput output(filename);

    for (int i(0); i < numVert; ++i) {
        Vector3 vertex = m_vertexArray[i];
        output.printf("v %f %f %f\n", vertex[0], vertex[1], vertex[2]);
    }

    for (int i(0); i < numFace; ++i) {
        Vector3int32 face = m_triArray[i];
        output.printf("f %d %d %d\n", face[0]+1, face[1]+1, face[2]+1);
    }
    output.commit(true);
    G3D::ArticulatedModel::clearCache();

    /*TextOutput file(filename);
    //file.printf("g name\n");
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
    file.commit();*/
}