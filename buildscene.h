#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include "util/vec.h"
#include "util/constants.h"
#include "shapes/shapes.h"
#include "shapes/instanced.h"
#include "gl/glcode.h"
#include "sampling.h"
#include "camera.h"
#include "lights.h"
#include "sceneobj.h"
#include "scenefile.h"
#include "shading.h"
#include "objloader/objloader.h"

// Assuming mesh->normals and mesh->face_normals are uninitialized
void computeMeshNormals(OBJShape* mesh)
{
    mesh->face_normals = (float*)malloc(sizeof(float) * mesh->num_indices);
    int num_face_normals = 0;
    mesh->normals = (float*)malloc(sizeof(float) * mesh->num_positions);
    mesh->num_normals = mesh->num_positions;
    for(int i = 0; i < mesh->num_normals; i++)
    {
        mesh->normals[i] = 0.0f;
    }

    for(int i = 0; i < mesh->num_indices; i += 3)
    {
        int i0, i1, i2;
        i0 = mesh->indices[i];
        i1 = mesh->indices[i+1];
        i2 = mesh->indices[i+2];        
        vec3 v0, v1, v2;
        int index;
        index = i0 * 3;
        vec3_assign(v0, mesh->positions[index], mesh->positions[index+1],
                    mesh->positions[index+2]);
        index = i1 * 3;
        vec3_assign(v1, mesh->positions[index], mesh->positions[index+1],
                    mesh->positions[index+2]);
        index = i2 * 3;
        vec3_assign(v2, mesh->positions[index], mesh->positions[index+1],
                    mesh->positions[index+2]);
        vec3 face_normal;
        calcTriangleNormal(face_normal, v0, v1, v2);
        mesh->face_normals[num_face_normals++] = face_normal[0];
        mesh->face_normals[num_face_normals++] = face_normal[1];
        mesh->face_normals[num_face_normals++] = face_normal[2];

    }
    mesh->num_face_normals = num_face_normals;


    for(int i = 0; i < mesh->num_indices; i += 3)
    {
        int normal_index, face_normal_index;
        face_normal_index = i;
        normal_index = mesh->indices[i] * 3;
        mesh->normals[normal_index] += mesh->face_normals[face_normal_index];
        mesh->normals[normal_index+1] += mesh->face_normals[face_normal_index+1];
        mesh->normals[normal_index+2] += mesh->face_normals[face_normal_index+2];

        normal_index = mesh->indices[i+1] * 3;
        mesh->normals[normal_index] += mesh->face_normals[face_normal_index];
        mesh->normals[normal_index+1] += mesh->face_normals[face_normal_index+1];
        mesh->normals[normal_index+2] += mesh->face_normals[face_normal_index+2];

        normal_index = mesh->indices[i+2] * 3;
        mesh->normals[normal_index] += mesh->face_normals[face_normal_index];
        mesh->normals[normal_index+1] += mesh->face_normals[face_normal_index+1];
        mesh->normals[normal_index+2] += mesh->face_normals[face_normal_index+2];                        
    }
    for(int i = 0; i < mesh->num_normals; i += 3)
    {
        vec3 normal = {mesh->normals[i], mesh->normals[i+1], mesh->normals[i+2]};
        vec3_normalize(normal, normal);
        mesh->normals[i] = normal[0];
        mesh->normals[i+1] = normal[1];
        mesh->normals[i+2] = normal[2];        
    }
}

void generateMeshTriangles(SceneObjects*so, const MeshEntry mesh_entry, const SceneMaterials *sm,
                           const SceneMeshes* s_meshes)
{
    mat3 rotation;
    eulerAngToMat3(rotation, mesh_entry.orientation);
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("%f ", rotation[i][j]);
        }
        printf("\n");
    }
    
    for(int i = 0; i < s_meshes->size; i++)
    {
        if(strcmp(mesh_entry.mesh_name, s_meshes->meshes[i].mesh_name) == 0)
        {
            OBJShape* mesh = &(s_meshes->meshes[i]);
            Material* mat = findMaterial(mesh_entry.mat_name, sm);
            for(int j = 0; j < mesh->num_indices; j += 3)
            {
                int i0, i1, i2;
                i0 = mesh->indices[j];
                i1 = mesh->indices[j + 1];
                i2 = mesh->indices[j + 2];

                vec3 v0, v1, v2;
                int index;
                index = i0 * 3;
                vec3_assign(v0, mesh->positions[index], mesh->positions[index+1],
                            mesh->positions[index+2]);
                index = i1 * 3;    
                vec3_assign(v1, mesh->positions[index], mesh->positions[index+1],
                            mesh->positions[index+2]);
                index = i2 * 3;
                vec3_assign(v2, mesh->positions[index], mesh->positions[index+1],
                            mesh->positions[index+2]);

                vec3_mult(v0, mesh_entry.scaling, v0);
                vec3_mult(v1, mesh_entry.scaling, v1);
                vec3_mult(v2, mesh_entry.scaling, v2);

                vec3 new_v0, new_v1, new_v2;
                mat3_mult_vec3(new_v0, rotation, v0);
                mat3_mult_vec3(new_v1, rotation, v1);
                mat3_mult_vec3(new_v2, rotation, v2);                
                
                vec3_add(new_v0, mesh_entry.location, new_v0);
                vec3_add(new_v1, mesh_entry.location, new_v1);
                vec3_add(new_v2, mesh_entry.location, new_v2);

                vec3 inv_scale = {1.0f/mesh_entry.scaling[0], 1.0f/mesh_entry.scaling[1], 1.0f/mesh_entry.scaling[2]};                
                if(!mesh_entry.smooth)
                {
                    FlatTriangle* mesh_tri = (FlatTriangle*)malloc(sizeof(FlatTriangle));
                    mesh_tri->i0 = i0;
                    mesh_tri->i1 = i1;
                    mesh_tri->i2 = i2;
                    vec3_copy(mesh_tri->v0, new_v0);
                    vec3_copy(mesh_tri->v1, new_v1);
                    vec3_copy(mesh_tri->v2, new_v2);                    
                    // Instead of computing an inverse transpose for transforming normals, I just inverted the scaling
                    // and ignored translation
                    vec3 face_normal = {mesh->face_normals[j], mesh->face_normals[j+1], mesh->face_normals[j+2]};
                    vec3 new_face_normal;
                    vec3_mult(face_normal, inv_scale, face_normal);
                    mat3_mult_vec3(new_face_normal, rotation, face_normal);
                    vec3_normalize(new_face_normal, new_face_normal);
                    vec3_copy(mesh_tri->normal, new_face_normal);
                    mesh_tri->shadow = mesh_entry.shadow;
                    mesh_tri->mesh_ptr = mesh;
                    mesh_tri->mat = mat;                                    
                    Object_t obj = {FLAT_TRIANGLE, mesh_tri};
                    SceneObjects_push_obj(so, obj);                    
                }else
                {
                    SmoothTriangle* mesh_tri = (SmoothTriangle*)malloc(sizeof(SmoothTriangle));
                    mesh_tri->i0 = i0;
                    mesh_tri->i1 = i1;
                    mesh_tri->i2 = i2;
                    vec3_copy(mesh_tri->v0, new_v0);
                    vec3_copy(mesh_tri->v1, new_v1);
                    vec3_copy(mesh_tri->v2, new_v2);                                        
                    vec3 normal, transformed_normal;
                    int index = i0 * 3;
                    vec3_assign(normal, mesh->normals[index],
                                mesh->normals[index+1], mesh->normals[index+2]);
                    vec3_mult(normal, inv_scale, normal);
                    mat3_mult_vec3(transformed_normal, rotation, normal);
                    vec3_normalize(transformed_normal, transformed_normal);
                    vec3_copy(mesh_tri->n0, transformed_normal);
                    index = i1 * 3;
                    vec3_assign(normal, mesh->normals[index],
                                mesh->normals[index+1], mesh->normals[index+2]);
                    vec3_mult(normal, inv_scale, normal);
                    mat3_mult_vec3(transformed_normal, rotation, normal);
                    vec3_normalize(transformed_normal, transformed_normal);                    
                    vec3_copy(mesh_tri->n1, transformed_normal);
                    index = i2 * 3;
                    vec3_assign(normal, mesh->normals[index],
                                mesh->normals[index+1], mesh->normals[index+2]);
                    vec3_mult(normal, inv_scale, normal);
                    mat3_mult_vec3(transformed_normal, rotation, normal);
                    vec3_normalize(transformed_normal, transformed_normal);                    
                    vec3_copy(mesh_tri->n2, transformed_normal);
                    mesh_tri->shadow = mesh_entry.shadow;
                    mesh_tri->mesh_ptr = mesh;
                    mesh_tri->mat = mat;                                    
                    Object_t obj = {SMOOTH_TRIANGLE, mesh_tri};
                    SceneObjects_push_obj(so, obj);                    
                }                
            }
        }
    }
}

void initSceneObjects(SceneObjects *so, SceneMaterials *sm, SceneMeshes* s_meshes,
                      const SceneLights* sl, const char* scenefile)
{
    FILE* fp;
#ifdef _MSC_VER
    fopen_s(&fp, scenefile, "r");
#else
    fp = fopen(scenefile, "r");
#endif
    if(fp)
    {
        int num_mat = parseMaterials(sm, fp);
        char buffer[128];
        while(getNextTokenInFile(buffer, fp))
        {
            if(buffer[0] == '#')
            {
                while(fgetc(fp) != '\n'){}
            }
            if(strcmp(buffer, "OBJECT") == 0)
            {
                getNextTokenInFile(buffer, fp);
                if(strcmp(buffer, "MESH") == 0)
                {
                    OBJShape* shapes;
                    int num_mesh;
                    char mesh_file_names[MAX_MESH][NAME_LENGTH];
                    int num_file_names = 0;
                    MeshEntry mesh_entry;
                    num_mesh = parseMesh(&mesh_entry, &shapes, mesh_file_names, &num_file_names, fp);
                    if(num_mesh > 0)
                    {
                        for(int i = 0; i < num_mesh; i++)
                        {
                            computeMeshNormals(&(shapes[i]));
                            SceneMeshes_push(s_meshes, shapes[i]);
                        }
                    }
                    free(shapes);
                    if(num_mesh != -1)
                    {
                        generateMeshTriangles(so, mesh_entry, sm, s_meshes);
                    }
                }else
                {
                    Object_t obj;                    
                    if(parsePrimitive(&obj, fp, sm, buffer))
                    {
                        SceneObjects_push_obj(so, obj);
                    }
                }
            }   
        }
    }
    
    //initSolidCylinder(so);
    //initInstanced(so);
    /*
    for(int i = 0; i < 50; ++i)
    {
        Sphere* sphere = (Sphere*)malloc(sizeof(Sphere));
        sphere->shadow = true;
        sphere->min_theta = 0;
        sphere->max_theta = PI;
        sphere->phi = PI;
        sphere->radius = 0.5;
        float r = (float)rand() / (float)RAND_MAX;
        float g = (float)rand() / (float)RAND_MAX;
        float b = (float)rand() / (float)RAND_MAX;
        vec3 color = {r, g, b};
        initDefaultMatteMat(&(sphere->mat), color);
        float x = (float)rand() / (float)RAND_MAX * 18.0f - 9.0f;
        float y = (float)rand() / (float)RAND_MAX * 18.0f - 10.0f;
        float z = (float)rand() / (float)RAND_MAX * -10.0f - 2.0f;
        vec3_assign(sphere->center, x, y, z);
        //so->obj_ptrs[so->num_obj] = sphere;
        //so->obj_types[so->num_obj] = SPHERE;
        //(so->num_obj)++;
        Object_t obj = {SPHERE, sphere};
        SceneObjects_push_obj(so, obj);
    }
    */
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == AREALIGHT)
        {
            AreaLight* area_light_ptr = (AreaLight*)(sl->light_ptrs[i]);
            Object_t obj = {area_light_ptr->obj_type, area_light_ptr->obj_ptr};
            SceneObjects_push_obj(so, obj);
        }
    }
}

void initAreaLights(SceneLights* sl)
{
    // Area light
    if(sl->num_lights == MAX_LIGHTS){return;}
    AreaLight* area_light_ptr = (AreaLight*)malloc(sizeof(AreaLight));
    area_light_ptr->intensity = 20.0f;
    vec3_copy(area_light_ptr->color, WHITE);
    vec3_assign(area_light_ptr->sample_point, 0.0f, 0.0f, 0.0f);

    /*
    // Rectangle
    Rectangle* rect = (Rectangle*)malloc(sizeof(Rectangle));
    rect->shadow = false;
    vec3_assign(rect->point, -2.0f, 0.0f, -6.0f);
    vec3_assign(rect->width, 2.0f, 0.0f, 0.0f);
    vec3_assign(rect->height, 0.0f, 4.0f, 0.0f);    
    vec3_copy(rect->normal, BACKWARD);
    vec3_copy(rect->mat.ce, area_light_ptr->color);
    rect->mat.ke = area_light_ptr->intensity;    
    rect->mat.mat_type = EMISSIVE;
    float width = sqrt(vec3_dot(rect->width, rect->width));
    float height = sqrt(vec3_dot(rect->height, rect->height));        
    
    Samples2D* samples = (Samples2D*)malloc(sizeof(Samples2D));
    samples->samples = NULL;
    samples->shuffled_indices = NULL;
    prepSample2DStruct(samples);
    genMultijitteredSamples(samples);

    area_light_ptr->pdf = 1.0f/(width * height);
    area_light_ptr->samples2D = samples;
    area_light_ptr->samples3D = NULL;
    area_light_ptr->obj_ptr = rect;
    area_light_ptr->obj_type = RECTANGLE;
    sl->shadow[sl->num_lights] = true;    
    sl->light_ptrs[sl->num_lights] = area_light_ptr;
    sl->light_types[sl->num_lights] = AREALIGHT;
    (sl->num_lights)++;
    */

    // Sphere
    Sphere* sphere = (Sphere*)malloc(sizeof(Sphere));
    sphere->shadow = false;
    //vec3_assign(sphere->center, 2.0f, 4.0f, 1.5f);
    vec3_assign(sphere->center, -0.9f, 5.0f, -3.1f);
    sphere->radius = 1.0f;
    sphere->min_theta = 0.0f;
    sphere->max_theta = (float)PI;
    sphere->phi = (float)PI;
    // TODO: fix the material pointer situation
    sphere->mat = (Material*)malloc(sizeof(Material));
    vec3_copy(sphere->mat->ce, area_light_ptr->color);    
    sphere->mat->ke = area_light_ptr->intensity;    
    sphere->mat->mat_type = EMISSIVE;

    Samples3D* samples = genHemisphereSamples(MULTIJITTERED, 1.0f);    

    area_light_ptr->pdf = 1.0f / (4.0f * (float)PI * sphere->radius * sphere->radius);
    area_light_ptr->samples2D = NULL;
    area_light_ptr->samples3D = samples;
    area_light_ptr->obj_ptr = sphere;
    area_light_ptr->obj_type = SPHERE;
    sl->shadow[sl->num_lights] = true;    
    sl->light_ptrs[sl->num_lights] = area_light_ptr;
    sl->light_types[sl->num_lights] = AREALIGHT;
    (sl->num_lights)++;
}

void initEnvLight(SceneLights* sl)
{
    if(sl->num_lights == MAX_LIGHTS){return;}
    EnvLight* env_light = (EnvLight*)malloc(sizeof(EnvLight));
    env_light->intensity = 1.0f;
    vec3_copy(env_light->color, WHITE);

    Samples3D* samples = genHemisphereSamples(MULTIJITTERED, 1.0f);

    env_light->samples3D = samples;
    sl->shadow[sl->num_lights] = true;
    sl->light_ptrs[sl->num_lights] = env_light;
    sl->light_types[sl->num_lights] = ENVLIGHT;
    (sl->num_lights)++;
    sl->env_light = env_light;
}

void initAmbLight(SceneLights *sl)
{
    sl->amb_light = (AmbientLight*)malloc(sizeof(AmbientLight));
    vec3_copy(sl->amb_light->color, WHITE);
    sl->amb_light->intensity = 0.0f;
    sl->amb_light->amb_occlusion = false;
}

void initBackgroundColor(SceneLights* sl)
{
    if(sl->env_light != NULL)
    {
        getIncRadiance(sl->bg_color, ENVLIGHT, sl->env_light);
    }else
    {
        vec3_copy(sl->bg_color, BLACK);
    }
}

void initSceneLights(SceneLights* sl)
{
    sl->env_light = NULL;    
    for(int i = 0; i < MAX_LIGHTS; i++)
    {
        sl->light_ptrs[i] = NULL;
    }
    sl->num_lights = 0;
    // Directional light
    if(sl->num_lights == MAX_LIGHTS){return;}
    DirLight* dir_light_ptr = (DirLight*)malloc(sizeof(DirLight));
    float intensity = 2.0f;
    vec3 direction = {10.0f, 10.0f, 10.0f};
    vec3_normalize(direction, direction);
    assignDirLight(dir_light_ptr, intensity, WHITE, direction);
    sl->shadow[sl->num_lights] = true;
    sl->light_ptrs[sl->num_lights] = dir_light_ptr;
    sl->light_types[sl->num_lights] = DIRECTIONAL;
    (sl->num_lights)++;

    
    // Point light
    /*
    if(sl->num_lights == MAX_LIGHTS){return;}
    PointLight* point_light_ptr = (PointLight*)malloc(sizeof(PointLight));
    intensity = 0.1f;
    vec3 point = {-3.0f, -5.0f, 0.0f};
    assignPointLight(point_light_ptr, intensity, WHITE, point);
    sl->shadow[sl->num_lights] = true;    
    sl->light_ptrs[sl->num_lights] = point_light_ptr;
    sl->light_ptrs[sl->num_lights] = POINTLIGHT;    
    (sl->num_lights)++;
    */


    initAreaLights(sl);
    //initEnvLight(sl);
    initAmbLight(sl);
    initBackgroundColor(sl);
}

void initScene(SceneObjects* so, SceneLights* sl, SceneMaterials* sm, SceneMeshes* s_meshes,
               const char* scenefile, const AccelType accel_type)
{
    // NOTE: initSceneObjects must be called after after initSceneLights if there are area lights
    initSceneLights(sl);    
    initSceneObjects(so, sm, s_meshes, sl, scenefile);
    mvNonGridObjToStart(so);
    printf("num_obj %d\n", so->num_obj);
    printf("non grid obj %d\n", so->num_non_grid_obj);

    so->accel = accel_type;
    if(so->accel == GRID)
    {
        UniformGrid* rg = UniformGrid_create(so->objects, &(so->num_obj), so->num_non_grid_obj, 2);
        so->accel_ptr = rg;
    }else if(so->accel == BVH)
    {
        BVHNode* tree;
        BVH_build(&tree, &(so->objects[so->num_non_grid_obj]), so->num_obj - so->num_non_grid_obj, 0);
        //int leaf_count = 0;
        //printBVH(tree, &leaf_count, 0);        
        so->accel_ptr = tree;
    }
}
