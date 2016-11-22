#pragma once

#include "GLFW/glfw3.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "../util/vec.h"
#include "../util/constants.h"
#include "../shapes/shapes.h"
#include "../shapes/instanced.h"
#include "../lights.h"
#include "../materials.h"
#include "../objloader/dbuffer.h"
#include "../sampling.h"
#include "../texture.h"

bool getPresetColor(vec3 r, const char buffer[])
{
    if(strcmp(buffer, "RED") == 0)
    {
        vec3_copy(r, RED);
        return true;
    }else if(strcmp(buffer, "GREEN") == 0)
    {
        vec3_copy(r, GREEN);
        return true;
    }else if(strcmp(buffer, "BLUE") == 0)
    {
        vec3_copy(r, BLUE);
        return true;
    }else if(strcmp(buffer, "WHITE") == 0)
    {
        vec3_copy(r, WHITE);
        return true;
    }else if(strcmp(buffer, "BLACK") == 0)
    {
        vec3_copy(r, BLACK);
        return true;
    }else if(strcmp(buffer, "YELLOW") == 0)
    {
        vec3_copy(r, YELLOW);
        return true;
    }else if(strcmp(buffer, "CYAN") == 0)
    {
        vec3_copy(r, CYAN);
        return true;
    }else if(strcmp(buffer, "PINK") == 0)
    {
        vec3_copy(r, PINK);
        return true;
    }else if(strcmp(buffer, "GREY") == 0)
    {
        vec3_copy(r, GREY);
        return true;
    }
    return false;
}

bool parseVec3(vec3 r, FILE* fp)
{
    char buffer[128];
    if(!getNextTokenInFile(buffer, fp)){return false;}
    r[0] = (float)atof(buffer);
    if(!getNextTokenInFile(buffer, fp)){return false;}
    r[1] = (float)atof(buffer);
    if(!getNextTokenInFile(buffer, fp)){return false;}
    r[2] = (float)atof(buffer);
    return true;
}

bool parseColor(vec3 r, FILE* fp)
{
    char buffer[128];
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(!getPresetColor(r, buffer))
    {
        r[0] = (float)atof(buffer);
        if(!getNextTokenInFile(buffer, fp)){return false;}
        r[1] = (float)atof(buffer);
        if(!getNextTokenInFile(buffer, fp)){return false;}
        r[2] = (float)atof(buffer);
    }
    return true;
}

// Need to pass SceneTextures
Texture* parseTexture(Scene* scene, FILE* fp)
{
    char buffer[128];

    if(!getNextTokenInFile(buffer, fp)){return NULL;}    // get file name
    Texture tex;
    if(!loadTexture(&tex, buffer))
    {
        return NULL;
    }
    return Scene_addTexture(scene, &tex, buffer);
}

bool parseTextures(Material* mat, Scene* scene, FILE* fp)
{
    char buffer[128];
    mat->tex_flags = NO_TEXTURE;
    Texture* tex_ptr;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Get texture type
    if(strcmp(buffer, "END") == 0)
    {
        return true;
    }
    if(strcmp(buffer, "NOISE") == 0)
    {
        mat->tex_flags |= NOISE;
        return true;
        //if(!getNextTokenInFile(buffer, fp)){return false;}
    }

    if(strcmp(buffer, "DIFFUSE_MAP") == 0)
    {
        tex_ptr = parseTexture(scene, fp);
        if(tex_ptr)
        {
            mat->tex_array[0] = tex_ptr;
            mat->tex_flags |= DIFFUSE;
        }
        if(!getNextTokenInFile(buffer, fp)){return false;}    // Get texture type
    }

    if(strcmp(buffer, "NORMAL_MAP") == 0)
    {
        tex_ptr = parseTexture(scene, fp);
        if(tex_ptr)
        {
            mat->tex_array[1] = tex_ptr;
            mat->tex_flags |= NORMAL;
        }
    }
    return true;
}

bool parseMatEntry(Material* mat, char** name, Scene* scene, FILE* fp)
{
    char buffer[128];
    char type_name[128];
    if(!getNextTokenInFile(type_name, fp)){return false;}
    mat->mat_type = getMatTypeFromString(type_name);
    mat->tex_flags = NO_TEXTURE;
    if(mat->mat_type == INVALID_MAT_TYPE)
    {
        fprintf(stderr, "Invalid material type %s.\n", type_name);
        if(!getNextTokenInFile(buffer, fp)){return false;}
        return false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip NAME
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get name
    //strcpy_s(*name, NAME_LENGTH, buffer);
    stringCopy(*name, NAME_LENGTH, buffer);
    if(mat->mat_type == EMISSIVE)
    {
        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word COLOR
        if(!parseColor(mat->ce, fp)){return false;}

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word INTENSITY
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->ke = (float)atof(buffer);
        mat->shadow = false;
        mat->h_samples = NULL;
    }else
    {
        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SHADOWED
        if(!getNextTokenInFile(buffer, fp)){return false;}
        if(strcmp(buffer, "yes") == 0)
        {
            mat->shadow = true;
        }else
        {
            mat->shadow = false;
        }

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word AMB_COLOR
        if(!parseColor(mat->ca, fp)){return false;}

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word AMB_CONSTANT
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->ka = (float)atof(buffer);

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word DIFF_COLOR
        if(!parseColor(mat->cd, fp)){return false;}

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word DIFF_CONSTANT
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->kd = (float)atof(buffer);

        if(mat->mat_type == MATTE)
        {
            mat->h_samples = genHemisphereSamples(MULTIJITTERED, 1.0f);
            goto Texture;
        }

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SPEC_COLOR
        if(!parseColor(mat->cs, fp)){return false;}

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SPEC_CONSTANT
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->ks = (float)atof(buffer);

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SPEC_EXPONENT
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->exp = (float)atof(buffer);

        if(mat->mat_type == PHONG)
        {
            mat->h_samples = NULL;
            goto Texture;
        }

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word REF_CONSTANT
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->kr = (float)atof(buffer);
        mat->h_samples = genHemisphereSamples(MULTIJITTERED, mat->exp);

        // TODO: Participating media should be its own category
        if(mat->mat_type == REFLECTIVE || mat->mat_type == PARTICIPATING)
        {
            goto Texture;
        }

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word IOR_IN
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->ior_in = (float)atof(buffer);

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word IOR_OUT
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->ior_out = (float)atof(buffer);

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CF_IN
        if(!parseColor(mat->cf_in, fp)){return false;}    

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CF_OUT
        if(!parseColor(mat->cf_out, fp)){return false;}

        if(mat->mat_type == TRANSPARENT)
        {
            goto Texture;
        }
        return false;

    Texture:
        parseTextures(mat, scene, fp);
        return true;
    }
}

int parseMaterials(Scene* scene, FILE* fp)
{
    char buffer[128];
    while(getNextTokenInFile(buffer, fp) && strcmp(buffer, "END_MATERIALS") != 0)
    {
        if(strcmp(buffer, "MATERIAL") == 0)
        {
            Material mat;
            char* name = (char*)malloc(sizeof(char) * MAX_NAME_LENGTH);
            parseMatEntry(&mat, &name, scene, fp);
            Scene_addMaterial(scene, &mat, name);
        }
    }
    return Scene_getNumMaterials(scene);
}

bool parseSphereEntry(Object_t* obj, FILE* fp, Scene* scene)
{
    char buffer[128];
    Sphere* sphere_ptr = (Sphere*)malloc(sizeof(Sphere));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        sphere_ptr->shadow = true;
    }else
    {
        sphere_ptr->shadow = false;
    }
    
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word RADIUS
    if(!getNextTokenInFile(buffer, fp)){return false;}
    sphere_ptr->radius = (float)atof(buffer);
    
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CENTER
    if(!parseVec3(sphere_ptr->center, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word PHI
    if(!getNextTokenInFile(buffer, fp)){return false;}
    sphere_ptr->phi = (float)atof(buffer);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MIN_THETA
    if(!getNextTokenInFile(buffer, fp)){return false;}
    sphere_ptr->min_theta = (float)atof(buffer);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MAX_THETA
    if(!getNextTokenInFile(buffer, fp)){return false;}
    sphere_ptr->max_theta = (float)atof(buffer);                

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get material name
    sphere_ptr->mat = Scene_findMaterial(scene, buffer);

    obj->ptr = sphere_ptr;
    obj->type = SPHERE;
    return true;    
}

bool parsePlaneEntry(Object_t* obj, FILE* fp, Scene* scene)
{
    char buffer[128];
    Plane* plane_ptr = (Plane*)malloc(sizeof(Plane));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        plane_ptr->shadow = true;
    }else 
    {
        plane_ptr->shadow = false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word POINT
    if(!parseVec3(plane_ptr->point, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word NORMAL
    if(!parseVec3(plane_ptr->normal, fp)){return false;}
    
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get material name    
    plane_ptr->mat = Scene_findMaterial(scene, buffer);    

    obj->ptr = plane_ptr;
    obj->type = PLANE;
    return true;
}

bool parseRectEntry(Object_t* obj, FILE* fp, Scene* scene)
{
    char buffer[128];
    Rectangle* rect_ptr = (Rectangle*)malloc(sizeof(Rectangle));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        rect_ptr->shadow = true;
    }else
    {
        rect_ptr->shadow = false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word POINT
    if(!parseVec3(rect_ptr->point, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word WIDTH
    if(!parseVec3(rect_ptr->width, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word HEIGHT
    if(!parseVec3(rect_ptr->height, fp)){return false;}

    vec3_cross(rect_ptr->normal, rect_ptr->width, rect_ptr->height);
    vec3_normalize(rect_ptr->normal, rect_ptr->normal);    

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get material name    
    rect_ptr->mat = Scene_findMaterial(scene, buffer);        

    obj->ptr = rect_ptr;
    obj->type = RECTANGLE;
    return true;
}

bool parseTriangleEntry(Object_t* obj,  FILE* fp, Scene* scene)
{
    char buffer[128];
    Triangle* tri_ptr = (Triangle*)malloc(sizeof(Triangle));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        tri_ptr->shadow = true;
    }else
    {
        tri_ptr->shadow = false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word V0
    if(!parseVec3(tri_ptr->v0, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word V1
    if(!parseVec3(tri_ptr->v1, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word V2
    if(!parseVec3(tri_ptr->v2, fp)){return false;}

    calcTriangleNormal(tri_ptr->normal, tri_ptr->v0, tri_ptr->v1, tri_ptr->v2);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get material name    
    tri_ptr->mat = Scene_findMaterial(scene, buffer);        

    obj->ptr = tri_ptr;
    obj->type = TRIANGLE;
    return true;
}

bool parseBoxEntry(Object_t* obj,  FILE* fp, Scene* scene)
{
    char buffer[128];
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    bool shadow;
    if(strcmp(buffer, "yes") == 0)
    {
        shadow = true;
    }else
    {
        shadow = false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word LENGTH
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get x span
    float length = (float)atof(buffer);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word HEIGHT
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get y span
    float height = (float)atof(buffer);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word WIDTH
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get z span
    float width = (float)atof(buffer);

    vec3 location;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word LOCATION
    if(!parseVec3(location, fp)){return false;}    
    vec3 scale;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SCALE
    if(!parseVec3(scale, fp)){return false;}        
    vec3 orientation;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word ORIENTATION
    if(!parseVec3(orientation, fp)){return false;}    

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get material name    
    Material* mat = Scene_findMaterial(scene, buffer);

    mat4 inv_scale, rotation, inv_rotation, inv_translation, tmp, inv_transform;
    mat4_scale_inverse(inv_scale, scale);
    mat4_translate(inv_translation, -location[0], -location[1], -location[2]);
    eulerAngToMat4(rotation, orientation);
    mat4_invert_rotation(inv_rotation, rotation);
    mat4_mult(tmp, inv_scale, inv_rotation);
    mat4_mult(inv_transform, tmp, inv_translation);        

    obj->ptr = initBox(inv_transform, length, height, width, mat, shadow);
    obj->type = INSTANCED;
    return true;
}

bool parseOpenCylEntry(Object_t* obj,  FILE* fp, Scene* scene)
{
    char buffer[128];
    bool shadow;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        shadow = true;
    }else
    {
        shadow = false;
    }
    float phi;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word PHI
    if(!getNextTokenInFile(buffer, fp)){return false;}
    phi = (float)atof(buffer);
    vec3 location;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word LOCATION
    if(!parseVec3(location, fp)){return false;}    
    vec3 scale;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SCALE
    if(!parseVec3(scale, fp)){return false;}        
    vec3 orientation;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word ORIENTATION
    if(!parseVec3(orientation, fp)){return false;}            
    NormalType normal_type;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word NORMAL_TYPE
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "OPEN") == 0)
    {
        normal_type = OPEN;
    }else if(strcmp(buffer, "CONVEX") == 0)
    {
        normal_type = CONVEX;
    }else if(strcmp(buffer, "CONCAVE") == 0)
    {
        normal_type = CONCAVE;
    }else
    {
        fprintf(stderr, "Invalid normal type %s\n", buffer);
        return false;
    }

    Material* mat;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get material name
    mat = Scene_findMaterial(scene, buffer);

    mat4 inv_scale, rotation, inv_rotation, inv_translation, tmp, inv_transform;
    mat4_scale_inverse(inv_scale, scale);
    mat4_translate(inv_translation, -location[0], -location[1], -location[2]);
    eulerAngToMat4(rotation, orientation);
    mat4_invert_rotation(inv_rotation, rotation);
    mat4_mult(tmp, inv_scale, inv_rotation);
    mat4_mult(inv_transform, tmp, inv_translation);

    obj->ptr = initOpenCylinder(inv_transform, phi, mat, normal_type, shadow);
    obj->type = INSTANCED;
    return true;
}

bool parseSolidCylEntry(Object_t* obj,  FILE* fp, Scene* scene)
{
    char buffer[128];
    bool shadow;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        shadow = true;
    }else
    {
        shadow = false;
    }
    vec3 location;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word LOCATION
    if(!parseVec3(location, fp)){return false;}    
    vec3 scale;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SCALE
    if(!parseVec3(scale, fp)){return false;}        
    vec3 orientation;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word ORIENTATION
    if(!parseVec3(orientation, fp)){return false;}            
    NormalType normal_type;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word NORMAL_TYPE
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "OPEN") == 0)
    {
        normal_type = OPEN;
    }else if(strcmp(buffer, "CONVEX") == 0)
    {
        normal_type = CONVEX;
    }else if(strcmp(buffer, "CONCAVE") == 0)
    {
        normal_type = CONCAVE;
    }else
    {
        fprintf(stderr, "Invalid normal type %s\n", buffer);
        return false;
    }

    Material* mat;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get material name
    mat = Scene_findMaterial(scene, buffer);

    mat4 inv_scale, rotation, inv_rotation, inv_translation, tmp, inv_transform;
    mat4_scale_inverse(inv_scale, scale);
    mat4_translate(inv_translation, -location[0], -location[1], -location[2]);
    eulerAngToMat4(rotation, orientation);
    mat4_invert_rotation(inv_rotation, rotation);
    mat4_mult(tmp, inv_scale, inv_rotation);
    mat4_mult(inv_transform, tmp, inv_translation);

    obj->ptr = initSolidCylinder(inv_transform, 1.0f, 1.0f, (float)PI, mat, shadow);
    obj->type = INSTANCED;
    return true;
}

bool parseDiskEntry(Object_t* obj,  FILE* fp, Scene* scene)
{
    char buffer[128];
    Disk* disk_ptr = (Disk*)malloc(sizeof(Disk));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        disk_ptr->shadow = true;
    }else
    {
        disk_ptr->shadow = false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MIN
    if(!parseVec3(disk_ptr->center, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MAX
    if(!parseVec3(disk_ptr->normal, fp)){return false;}
    printVec3WithText("disk->normal", disk_ptr->normal);
    vec3_normalize(disk_ptr->normal, disk_ptr->normal);
    printVec3WithText("disk->normal normalized ", disk_ptr->normal);    

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word RADIUS
    if(!getNextTokenInFile(buffer, fp)){return false;}
    disk_ptr->radius = (float)atof(buffer);    

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get material name    
    disk_ptr->mat = Scene_findMaterial(scene, buffer);                        

    obj->ptr = disk_ptr;
    obj->type = DISK;
    return true;
}

bool parseTorusEntry(Object_t* obj,  FILE* fp, Scene* scene)
{
    char buffer[128];
    bool shadow = true;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        shadow = true;
    }else
    {
        shadow = false;
    }

    float swept_radius;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SWEPT_RADIUS
    if(!getNextTokenInFile(buffer, fp)){return false;}
    swept_radius = (float)atof(buffer);
    float tube_radius;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word TUBE_RADIUS
    if(!getNextTokenInFile(buffer, fp)){return false;}
    tube_radius = (float)atof(buffer);    
    float phi;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word PHI
    if(!getNextTokenInFile(buffer, fp)){return false;}
    phi = (float)atof(buffer);

    vec3 location;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word LOCATION
    if(!parseVec3(location, fp)){return false;}    
    vec3 scale;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SCALE
    if(!parseVec3(scale, fp)){return false;}        
    vec3 orientation;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word ORIENTATION
    if(!parseVec3(orientation, fp)){return false;}        
    
    Material* mat;
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!getNextTokenInFile(buffer, fp)){return false;}    // get material name
    printf("buffer name %s\n", buffer);
    mat = Scene_findMaterial(scene, buffer);

    mat4 inv_scale, rotation, inv_rotation, inv_translation, tmp, inv_transform;
    mat4_scale_inverse(inv_scale, scale);
    mat4_translate(inv_translation, -location[0], -location[1], -location[2]);
    eulerAngToMat4(rotation, orientation);
    mat4_invert_rotation(inv_rotation, rotation);
    mat4_mult(tmp, inv_scale, inv_rotation);
    mat4_mult(inv_transform, tmp, inv_translation);

    obj->ptr = initTorus(inv_transform, swept_radius, tube_radius, phi, mat, shadow);
    obj->type = INSTANCED;
    return true;
}



typedef struct
{
    bool shadow;
    bool smooth;        
    vec3 scaling;
    vec3 location;
    vec3 orientation; // Euler angle
    char mesh_name[NAME_LENGTH];
    char mat_name[NAME_LENGTH];
}MeshEntry;

// Return -1 if the file cannot be read
int parseMesh(MeshEntry* mesh_entry, OBJShape** shapes, char mesh_file_names[][NAME_LENGTH],
              int* num_file_names, FILE* fp)
{
    char buffer[128];
    if(!getNextTokenInFile(buffer, fp)){return -1;}    // Skip over the word FILE_NAME
    if(!getNextTokenInFile(buffer, fp)){return -1;}    // get file name
    bool file_not_read = true;
    for(int i = 0; i < *num_file_names; i++)
    {
        if(strcmp(buffer, mesh_file_names[i]) == 0)
        {
            file_not_read = false;
            break;
        }
    }
    
    int num_mesh = 0;
    if(file_not_read)
    {
        double start, end;
        start = glfwGetTime();
        num_mesh = loadOBJ(shapes, buffer);
        end = glfwGetTime();
        printf("Loading %s took %f sec.\n", buffer, end - start);
    }else
    {
        //strcpy_s(mesh_file_names[(*num_file_names)++], NAME_LENGTH, buffer);
        stringCopy(mesh_file_names[(*num_file_names)++], NAME_LENGTH, buffer);
    }

    if(num_mesh != -1)
    {
        int len = strcspn(buffer, ".");
        //strncpy_s(mesh_entry->mesh_name, NAME_LENGTH, buffer, len);
        stringNCopy(mesh_entry->mesh_name, NAME_LENGTH, buffer, len);
        mesh_entry->mesh_name[len] = '\0';
    }    

    if(!getNextTokenInFile(buffer, fp)){return -1;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return -1;}
    if(strcmp(buffer, "yes") == 0)
    {
        mesh_entry->shadow = true;
    }else
    {
        mesh_entry->shadow = false;
    }

    if(!getNextTokenInFile(buffer, fp)){return -1;}    // Skip over the word SMOOTH
    if(!getNextTokenInFile(buffer, fp)){return -1;}
    if(strcmp(buffer, "yes") == 0)
    {
        mesh_entry->smooth = true;
    }else
    {
        mesh_entry->smooth = false;
    }    

    if(!getNextTokenInFile(buffer, fp)){return -1;}    // Skip over the word SCALING
    if(!parseVec3(mesh_entry->scaling, fp)){return -1;}

    if(!getNextTokenInFile(buffer, fp)){return -1;}    // Skip over the word LOCATION
    if(!parseVec3(mesh_entry->location, fp)){return -1;}

    if(!getNextTokenInFile(buffer, fp)){return -1;}    // Skip over the word ORIENTATION
    if(!parseVec3(mesh_entry->orientation, fp)){return -1;}      

    if(!getNextTokenInFile(buffer, fp)){return -1;}    // Skip over the word MATERIAL
    if(!getNextTokenInFile(buffer, fp)){return -1;}
    //strcpy_s(mesh_entry->mat_name, NAME_LENGTH, buffer);
    stringCopy(mesh_entry->mat_name, NAME_LENGTH, buffer);
    
    return num_mesh;
}

bool parsePrimitive(Object_t* obj, FILE* fp, Scene* scene, const char* prim_name)
{
    bool parse_status = false;
    if(strcmp(prim_name, "SPHERE") == 0)
    {
        parse_status = parseSphereEntry(obj, fp, scene);
    }else if(strcmp(prim_name, "PLANE") == 0)
    {
        parse_status = parsePlaneEntry(obj, fp, scene);
    }else if(strcmp(prim_name, "RECTANGLE") == 0)
    {
        parse_status = parseRectEntry(obj, fp, scene);
    }else if(strcmp(prim_name, "TRIANGLE") == 0)
    {
        parse_status = parseTriangleEntry(obj, fp, scene);
    }else if(strcmp(prim_name, "BOX") == 0)
    {
        parse_status = parseBoxEntry(obj, fp, scene);
    }else if(strcmp(prim_name, "OPENCYLINDER") == 0)
    {
        parse_status = parseOpenCylEntry(obj, fp, scene);
    }else if(strcmp(prim_name, "DISK") == 0)
    {
        parse_status = parseDiskEntry(obj, fp, scene);
    }else if(strcmp(prim_name, "TORUS") == 0)
    {
        parse_status = parseTorusEntry(obj, fp, scene);
    }else if(strcmp(prim_name, "SOLIDCYLINDER") == 0)
    {
        parse_status = parseSolidCylEntry(obj, fp, scene);
    }
    return parse_status;
}
