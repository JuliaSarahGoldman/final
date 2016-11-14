/**
  \file SimpleMesh.h
   
   Defines mesh and its operators, including edge collapse and beveling.
 */
/*
    John Freeman 
    Jose Rivas-Garcia 
    Julia Goldman 
    Matheus de Carvalho Souza
*/
#pragma once
#include <G3D/G3DAll.h>

class SimpleMesh { 
protected:
    Array<Vector3> m_vertexArray;
    Array<Vector3int32> m_triArray;

public:
    void toObj(String filename);
    SimpleMesh(Array<Vector3> vertices, Array<Vector3int32> triangles);
    //~SimpleMesh();
};