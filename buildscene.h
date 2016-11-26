#pragma once

#include "GLFW/glfw3.h"
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
#include "scene/scenedata.h"
#include "scene/scene.h"
#include "scene/scenefile.h"
#include "shading.h"
#include "objloader/objloader.h"
#include "mesh.h"
#include "bsphere.h"

//#define CORNELL_BOX

// Assuming mesh->normals and mesh->face_normals are uninitialized
void calcTriangleNormals(Mesh* mesh)
{
    mesh->face_normals = (float*)malloc(sizeof(float) * mesh->num_indices);
    int num_face_normals = 0;

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

    if(mesh->num_normals > 0)
    {
        return;
    }
    mesh->normals = (float*)malloc(sizeof(float) * mesh->num_positions);
    mesh->num_normals = mesh->num_positions;
    for(int i = 0; i < mesh->num_normals; i++)
    {
        mesh->normals[i] = 0.0f;
    }    
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

void calcTangentVec(Mesh* mesh)
{
    int vert_count = mesh->num_positions / 3;    
    vec3* tangents = (vec3*)malloc(sizeof(vec3) * vert_count);
    vec3* binormals = (vec3*)malloc(sizeof(vec3) * vert_count);
    int8_t* det = (int8_t*)malloc(sizeof(int8_t) * vert_count);
    for(int i = 0; i < mesh->num_positions / 3; i++)
    {
        vec3_copy(tangents[i], BLACK);
        vec3_copy(binormals[i], BLACK);
    }
    for(int i = 0; i < mesh->num_indices; i+=3)
    {
        int i0 = mesh->indices[i];
        int i1 = mesh->indices[i+1];
        int i2 = mesh->indices[i+2];
        
        vec3 v0 = {mesh->positions[i0*3], mesh->positions[i0*3 + 1], mesh->positions[i0*3 + 2]};
        vec3 v1 = {mesh->positions[i1*3], mesh->positions[i1*3 + 1], mesh->positions[i1*3 + 2]};
        vec3 v2 = {mesh->positions[i2*3], mesh->positions[i2*3 + 1], mesh->positions[i2*3 + 2]};
        vec2 uv0 = {mesh->texcoords[i0*2], mesh->texcoords[i0*2 + 1]};
        vec2 uv1 = {mesh->texcoords[i1*2], mesh->texcoords[i1*2 + 1]};
        vec2 uv2 = {mesh->texcoords[i2*2], mesh->texcoords[i2*2 + 1]};                
        
        vec3 q0, q1;
        vec3_sub(q0, v1, v0);
        vec3_sub(q1, v2, v0);
        float s0 = uv1[0] - uv0[0];
        float s1 = uv2[0] - uv0[0];
        float t0 = uv1[1] - uv0[1];
        float t1 = uv2[1] - uv0[1];

        vec3 t1q0;
        vec3_scale(t1q0, q0, t1);
        vec3 t0q1;
        vec3_scale(t0q1, q1, t0);
        vec3 neg_s1q0;
        vec3_scale(neg_s1q0, q0, -s1);
        vec3 s0q1;
        vec3_scale(s0q1, q1, s0);
        
        vec3 u_axis, v_axis;
        vec3_sub(u_axis, t1q0, t0q1);
        vec3_add(v_axis, neg_s1q0, s0q1);

        vec3_add(tangents[i0], tangents[i0], u_axis);
        vec3_add(binormals[i0], binormals[i0], v_axis);
        vec3_add(tangents[i1], tangents[i1], u_axis);
        vec3_add(binormals[i1], binormals[i1], v_axis);
        vec3_add(tangents[i2], tangents[i2], u_axis);
        vec3_add(binormals[i2], binormals[i2], v_axis);        
    }

    for(int i = 0; i < vert_count; i++)
    {
        vec3 normal = {mesh->normals[i*3], mesh->normals[i*3 + 1], mesh->normals[i*3 + 2]};

        // Ensure tangent is perpendicular to the normal
        vec3 displacement;
        vec3_scale(displacement, normal, vec3_dot(normal, tangents[i]));
        vec3_sub(tangents[i], tangents[i], displacement);
        vec3_normalize(tangents[i], tangents[i]);
    }
    // TODO: find a way to check if the mesh already has tangents and dets
    mesh->tangents = tangents;
    free(binormals);
}

int generateMeshTriangles(Scene* scene, const MeshEntry mesh_entry)
{
    mat3 rotation;
    eulerAngToMat3(rotation, mesh_entry.orientation);

    int num_mesh_found = 0;
    Mesh** meshes = Scene_findMeshes(&num_mesh_found, scene, mesh_entry.mesh_name);
    Material* mat = Scene_findMaterial(scene, mesh_entry.mat_name);
    mat3 inv_scale_mat, normal_mat;
    mat3_scale_inverse(inv_scale_mat, mesh_entry.scaling);
    mat3_mult(normal_mat, rotation, inv_scale_mat);
    mat3* normal_mat_ptr =  Scene_addTransform(scene, &normal_mat);    
    int triangle_count = 0;
    for(int i = 0; i < num_mesh_found; i++)
    {
        Mesh* mesh = meshes[i];
        triangle_count += mesh->num_indices / 3;
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
                //vec3_mult(face_normal, inv_scale, face_normal);
                mat3_mult_vec3(new_face_normal, normal_mat, face_normal);
                vec3_normalize(new_face_normal, new_face_normal);
                vec3_copy(mesh_tri->normal, new_face_normal);                     
                mesh_tri->shadow = mesh_entry.shadow;
                mesh_tri->mesh_ptr = mesh;
                mesh_tri->mat = mat;                                    
                Object_t obj = {FLAT_TRIANGLE, mesh_tri};
                Scene_addObject(scene, &obj);
            }else
            {
                SmoothTriangle* mesh_tri = (SmoothTriangle*)malloc(sizeof(SmoothTriangle));
                mesh_tri->normal_mat = normal_mat_ptr;

                mesh_tri->i0 = i0;
                mesh_tri->i1 = i1;
                mesh_tri->i2 = i2;
                vec3_copy(mesh_tri->v0, new_v0);
                vec3_copy(mesh_tri->v1, new_v1);
                vec3_copy(mesh_tri->v2, new_v2);

                mesh_tri->shadow = mesh_entry.shadow;
                mesh_tri->mesh_ptr = mesh;
                mesh_tri->mat = mat;                                    
                Object_t obj = {SMOOTH_TRIANGLE, mesh_tri};
                Scene_addObject(scene, &obj);                    
            }                
        }        
    }
    return triangle_count;
}

void initSceneObjects(Scene* scene, const char* scenefile)
{
    SceneLights* sl = &(scene->lights);
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == AREALIGHT)
        {
            AreaLight* area_light_ptr = (AreaLight*)(sl->light_ptrs[i]);
            Object_t obj = {area_light_ptr->obj_type, area_light_ptr->obj_ptr};
            Scene_addObject(scene, &obj);
        }
    }
    
    FILE* fp;
#ifdef _MSC_VER
    fopen_s(&fp, scenefile, "r");
#else
    fp = fopen(scenefile, "r");
#endif
    if(!fp)
    {
        return;
    }
        
    int num_mat = parseMaterials(scene, fp);
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
                OBJShape* shapes = NULL;
                int num_mesh;
                char mesh_file_names[MAX_MESH][NAME_LENGTH];
                int num_file_names = 0;
                MeshEntry mesh_entry;
                num_mesh = parseMesh(&mesh_entry, &shapes, mesh_file_names, &num_file_names, fp);
                for(int i = 0; i < num_mesh; i++)
                {
                    Mesh mesh;
                    Mesh_copyOBJShape(&mesh, &(shapes[i]));
                    calcTriangleNormals(&mesh);
                    Material* mat = Scene_findMaterial(scene, mesh_entry.mat_name);
                    if(mat->tex_flags & NORMAL)
                    {
                        calcTangentVec(&mesh);
                    }
                    //calcTangentVec(&mesh);
                    Scene_addMesh(scene, &mesh);
                }
                if(shapes){free(shapes);}
                if(num_mesh != -1)
                {
                    generateMeshTriangles(scene, mesh_entry);
                }
            }else
            {
                Object_t obj;
                if(parsePrimitive(&obj, fp, scene, buffer))
                {
                    Scene_addObject(scene, &obj);
                }
            }
        }
    }
}

void initAreaLights(SceneLights* sl)
{

    // Area light
    if(sl->num_lights == MAX_LIGHTS){return;}
    AreaLight* area_light_ptr = (AreaLight*)malloc(sizeof(AreaLight));
    // Cornell rectangle light intensity
    area_light_ptr->intensity = 55.0f;
    //area_light_ptr->intensity = 100.0f;
    //area_light_ptr->intensity = 20.0f;
    // Photon map intensity
    //area_light_ptr->intensity = 100.0f;
    vec3_assign(area_light_ptr->color, 1.0f, 0.85f, 0.5f);
    vec3_assign(area_light_ptr->sample_point, 0.0f, 0.0f, 0.0f);

    // Rectangle
    Rectangle* rect = (Rectangle*)malloc(sizeof(Rectangle));
    rect->mat = (Material*)malloc(sizeof(Material)); // NOTE: memory leak?
    rect->shadow = false;

    vec3_assign(rect->point, 213.0f, 547.0f, -227.0f);
    vec3_assign(rect->width, 130.0f, 0.0f, 0.0f);
    vec3_assign(rect->height, 0.0f, 0.0f, -105.0f);

    //vec3_assign(rect->point, 213.0f, 800.0f, -227.0f);
    /*
    // pm test 
    vec3_assign(rect->point, 213.0f, 900.0f, -227.0f);
    vec3_assign(rect->width, 300.0f, 0.0f, 0.0f);
    vec3_assign(rect->height, 0.0f, 0.0f, -250.0f);
    */
    vec3_copy(rect->normal, DOWN);
    vec3_copy(rect->mat->ce, area_light_ptr->color);
    rect->mat->ke = area_light_ptr->intensity;
    rect->mat->mat_type = EMISSIVE;
    float width = sqrt(vec3_dot(rect->width, rect->width));
    float height = sqrt(vec3_dot(rect->height, rect->height));

    area_light_ptr->flux = area_light_ptr->intensity * width * height * PI;

    Samples2D* unit_square_samples = (Samples2D*)malloc(sizeof(Samples2D));
    unit_square_samples->samples = NULL;
    genMultijitteredSamples(unit_square_samples);

    area_light_ptr->pdf = 1.0f/(width * height);
    area_light_ptr->samples2D = unit_square_samples;
    area_light_ptr->samples3D = NULL;
    area_light_ptr->obj_ptr = rect;
    area_light_ptr->obj_type = RECTANGLE;
    sl->shadow[sl->num_lights] = true;
    sl->light_ptrs[sl->num_lights] = area_light_ptr;
    sl->light_types[sl->num_lights] = AREALIGHT;
    (sl->num_lights)++;


    // Area light 2
    /*
    if(sl->num_lights == MAX_LIGHTS){return;}    
    area_light_ptr = (AreaLight*)malloc(sizeof(AreaLight));
    area_light_ptr->intensity = 55.0f;
    vec3_assign(area_light_ptr->color, 1.0f, 0.85f, 0.5f);
    vec3_assign(area_light_ptr->sample_point, 0.0f, 0.0f, 0.0f);

    // Rectangle
    rect = (Rectangle*)malloc(sizeof(Rectangle));
    rect->mat = (Material*)malloc(sizeof(Material)); // NOTE: memory leak?
    rect->shadow = false;
    vec3_assign(rect->point, 50.0f, 260.0f, -600.0f);
    vec3_assign(rect->width, 200.0f, 0.0f, 0.0f);
    vec3_assign(rect->height, 0.0f, 150.0f, 0.0f);
    vec3_copy(rect->normal, BACKWARD);
    vec3_copy(rect->mat->ce, area_light_ptr->color);
    rect->mat->ke = area_light_ptr->intensity;
    rect->mat->mat_type = EMISSIVE;
    width = sqrt(vec3_dot(rect->width, rect->width));
    height = sqrt(vec3_dot(rect->height, rect->height));
    area_light_ptr->flux = area_light_ptr->intensity * width * height * PI;

    area_light_ptr->pdf = 1.0f/(width * height);
    area_light_ptr->samples2D = unit_square_samples;
    area_light_ptr->samples3D = NULL;
    area_light_ptr->obj_ptr = rect;
    area_light_ptr->obj_type = RECTANGLE;
    sl->shadow[sl->num_lights] = true;
    sl->light_ptrs[sl->num_lights] = area_light_ptr;
    sl->light_types[sl->num_lights] = AREALIGHT;
    (sl->num_lights)++;
    */
    /*
    // Sphere
    AreaLight* sphere_light_ptr = (AreaLight*)malloc(sizeof(AreaLight));
    sphere_light_ptr->intensity = 1000000.0f;
    vec3_assign(sphere_light_ptr->color, 1.0f, 0.85f, 0.5f);
    vec3_assign(sphere_light_ptr->sample_point, 0.0f, 0.0f, 0.0f);

    Sphere* sphere = (Sphere*)malloc(sizeof(Sphere));
    sphere->shadow = false;
    //vec3_assign(sphere->center, -0.9f, 5.0f, -3.1f);
    //vec3_assign(sphere->center, -0.9f, 5.0f, -1.0f);
    // Photon map location
    vec3_assign(sphere->center, 500.0f, 225.0f, -290.0f);
    sphere->radius = 20.0f;
    sphere->min_theta = 0.0f;
    sphere->max_theta = (float)PI;
    sphere->phi = (float)PI;
    // TODO: fix the material pointer situation
    sphere->mat = (Material*)malloc(sizeof(Material));
    vec3_copy(sphere->mat->ce, sphere_light_ptr->color);
    sphere->mat->ke = sphere_light_ptr->intensity;
    sphere->mat->mat_type = EMISSIVE;

    Samples3D* h_samples = genHemisphereSamples(MULTIJITTERED, 1.0f);


    sphere_light_ptr->pdf = 1.0f / (4.0f * (float)PI * sphere->radius * sphere->radius);
    sphere_light_ptr->samples2D = NULL;
    sphere_light_ptr->samples3D = h_samples;
    sphere_light_ptr->obj_ptr = sphere;
    sphere_light_ptr->obj_type = SPHERE;
    sl->shadow[sl->num_lights] = true;
    sl->light_ptrs[sl->num_lights] = sphere_light_ptr;
    sl->light_types[sl->num_lights] = AREALIGHT;
    (sl->num_lights)++;
    */
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
        vec3 place_holder; // NOTE: just there to fill out the function parameters. Actual value doesn't matter.
        getIncRadiance(sl->bg_color, ENVLIGHT, sl->env_light, place_holder);
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
    float intensity = 0.0f;
    //vec3 direction = {10.0f, 10.0f, 10.0f};
    vec3 direction = {1.0f, 0.0f, 10.0f};
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
    point_light_ptr->dist_atten = true;
    //intensity = 1000000.0f;
    intensity = 100000.0f;
    point_light_ptr->flux = intensity * 4.0f * PI;
    //vec3 point = {500.0f, 225.0f, -290.0f};
    vec3 point = {450.0f, 225.0f, -225.0f};
    assignPointLight(point_light_ptr, intensity, WHITE, point);
    point_light_ptr->proj_map = NULL;
    sl->shadow[sl->num_lights] = true;    
    sl->light_ptrs[sl->num_lights] = point_light_ptr;
    sl->light_types[sl->num_lights] = POINTLIGHT;    
    (sl->num_lights)++;
    */
    /*
    if(sl->num_lights == MAX_LIGHTS){return;}
    point_light_ptr = (PointLight*)malloc(sizeof(PointLight));
    point_light_ptr->dist_atten = false;    
    intensity = 1000000.0f;
    vec3_assign(point, 500.0f, 225.0f, -290.0f);
    //vec3 point = {250.0f, 490.0f, -290.0f};
    //vec3 point = {250.0f, 490.0f, -490.0f};
    assignPointLight(point_light_ptr, intensity, WHITE, point);
    sl->shadow[sl->num_lights] = true;    
    sl->light_ptrs[sl->num_lights] = point_light_ptr;
    sl->light_types[sl->num_lights] = POINTLIGHT;    
    (sl->num_lights)++;    
    */

    //initAreaLights(sl);
    initEnvLight(sl);
    initAmbLight(sl);
    initBackgroundColor(sl);
}

void buildSceneAccel(Scene *scene)
{
    SceneObjects* so = &(scene->objects);    
    double start, end;    
    if(so->accel == GRID)
    {
        start = glfwGetTime();
        UniformGrid* rg = UniformGrid_create(so->objects, &(so->num_obj), so->num_non_grid_obj, 2);
        end = glfwGetTime();
        printf("Uniform grid build time: %f sec\n", end - start);        
        so->accel_ptr = rg;
    }else if(so->accel == BVH)
    {
        BVHNode* tree;
        double start, end;
        start = glfwGetTime();
        BVH_build(&tree, &(so->objects[so->num_non_grid_obj]), so->num_obj - so->num_non_grid_obj);
        end = glfwGetTime();
        printf("BVH build time: %f sec\n", end - start);
        //int leaf_count = 0;
        //printBVH(tree, &leaf_count, 0);        
        so->accel_ptr = tree;
    }else if(so->accel == BVH4)
    {
        BVHNode4* tree;
        double start, end;
        start = glfwGetTime();
        BVH4_build(&tree, &(so->objects[so->num_non_grid_obj]), so->num_obj - so->num_non_grid_obj);
        end = glfwGetTime();
        printf("BVH4 build time: %f sec\n", end - start);
        so->accel_ptr = tree;
    }
}

void initScene(Scene* scene, const char* scenefile, const AccelType accel_type)
{
    // NOTE: initSceneObjects must be called after after initSceneLights if there are area lights
    initSceneLights(&(scene->lights));    
    initSceneObjects(scene, scenefile);
    SceneObjects* so = &(scene->objects);
    so->accel = accel_type;
    mvNonGridObjToStart(so);
    printf("num_obj %d\n", so->num_obj);
    printf("non bounded obj %d\n", so->num_non_grid_obj);

    // NOTE: taken out to main so that projection map can be built before the
    // mesh triangles are scrambled.
    //buildSceneAccel(scene);
}
