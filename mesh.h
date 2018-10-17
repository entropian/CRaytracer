#pragma once
#include <stdio.h>
#include <stdint.h>
//#include "objloader/objloader.h"
#include "util/vec.h"
#include "util/constants.h"

struct OBJShape_s;
typedef struct OBJShape_s OBJShape;

typedef struct Mesh_s
{
    float* positions;
    float* normals;
    float* texcoords;    
    int* indices;
    float* face_normals;
    vec3* tangents;
    int num_positions, num_normals, num_texcoords, num_indices, num_face_normals;
    char mat_name[NAME_LENGTH];
    char mesh_name[NAME_LENGTH];    
}Mesh;

void Mesh_destroy(Mesh* mesh);
void Mesh_copyOBJShape(Mesh* mesh, const OBJShape* obj_shape);
