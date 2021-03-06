#include "mesh.h"
#include "util/util.h"
#define OBJ_LOADER_IMPLEMENTATION    
#include "objloader/objloader.h"

void Mesh_destroy(Mesh* mesh)
{
    if(mesh->num_positions > 0){free(mesh->positions);}
    if(mesh->num_normals > 0){free(mesh->normals);}
    if(mesh->num_texcoords > 0){free(mesh->texcoords);}
    if(mesh->num_indices > 0){free(mesh->indices);}
    if(mesh->num_face_normals > 0){free(mesh->face_normals);}

    mesh->num_positions = 0;
    mesh->num_normals = 0;
    mesh->num_texcoords = 0;
    mesh->num_indices = 0;
    mesh->num_face_normals = 0;
}

void Mesh_copyOBJShape(Mesh* mesh, const OBJShape* obj_shape)
{
    mesh->positions = obj_shape->positions;
    mesh->normals = obj_shape->normals;
    mesh->texcoords = obj_shape->texcoords;
    mesh->indices = obj_shape->indices;
    mesh->face_normals = NULL;
    mesh->num_positions = obj_shape->num_positions;
    mesh->num_normals = obj_shape->num_normals;
    mesh->num_texcoords = obj_shape->num_texcoords;
    mesh->num_indices = obj_shape->num_indices;
    mesh->num_face_normals = 0;
    // TODO: Platform
    //strcpy_s(mesh->mat_name, NAME_LENGTH, obj_shape->mat_name);    
    //strcpy_s(mesh->mesh_name, NAME_LENGTH, obj_shape->mesh_name);
    stringCopy(mesh->mat_name, NAME_LENGTH, obj_shape->mat_name);    
    stringCopy(mesh->mesh_name, NAME_LENGTH, obj_shape->mesh_name);
}

