

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
#include "imagefile.h"

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
        if(normal[0] == 0.0f && normal[1] == 0.0f && normal[2])
        {
            printf("degenerate calc normal\n");
        }        
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
    bool mesh_has_no_texcoord = false;
    if(mesh->num_texcoords == 0)
    {
        mesh_has_no_texcoord = true;
        mesh->texcoords = (float*)malloc(sizeof(float) * 2 *vert_count);
        mesh->num_texcoords = 2 * vert_count;
    }
    
    for(int i = 0; i < mesh->num_indices; i+=3)
    {
        int i0 = mesh->indices[i];
        int i1 = mesh->indices[i+1];
        int i2 = mesh->indices[i+2];
        
        vec3 v0 = {mesh->positions[i0*3], mesh->positions[i0*3 + 1], mesh->positions[i0*3 + 2]};
        vec3 v1 = {mesh->positions[i1*3], mesh->positions[i1*3 + 1], mesh->positions[i1*3 + 2]};
        vec3 v2 = {mesh->positions[i2*3], mesh->positions[i2*3 + 1], mesh->positions[i2*3 + 2]};
        vec2 uv0, uv1, uv2;
        if(mesh_has_no_texcoord)
        {
            vec2_assign(uv0, mesh->texcoords[i0*2], mesh->texcoords[i0*2 + 1]);
            vec2_assign(uv1, mesh->texcoords[i1*2], mesh->texcoords[i1*2 + 1]);
            vec2_assign(uv2, mesh->texcoords[i2*2], mesh->texcoords[i2*2 + 1]);
            vec2 zeroes = {0.0f, 0.0f};
            if(vec2_equal(uv0, zeroes) && vec2_equal(uv1, zeroes) && vec2_equal(uv2, zeroes))
            {
                vec2_assign(uv0, 0.0f, 0.0f);
                vec2_assign(uv1, 1.0f, 0.0f);
                vec2_assign(uv2, 1.0f, 1.0f);
                /*
                mesh->texcoords[i0*2] = uv0[0];
                mesh->texcoords[i0*2 + 1] = uv0[1];
                mesh->texcoords[i1*2] = uv1[0];
                mesh->texcoords[i1*2 + 1] = uv1[1];
                mesh->texcoords[i2*2] = uv2[0];
                mesh->texcoords[i2*2 + 1] = uv2[1];
                */
            }
        }else
        {
            /*
            vec2_assign(uv0, mesh->texcoords[i0*2], mesh->texcoords[i0*2 + 1]);
            vec2_assign(uv1, mesh->texcoords[i1*2], mesh->texcoords[i1*2 + 1]);
            vec2_assign(uv2, mesh->texcoords[i2*2], mesh->texcoords[i2*2 + 1]);
            */
            vec2_assign(uv0, 0.0f, 0.0f);
            vec2_assign(uv1, 1.0f, 0.0f);
            vec2_assign(uv2, 1.0f, 1.0f);
            /*
            mesh->texcoords[i0*2] = uv0[0];
            mesh->texcoords[i0*2 + 1] = uv0[1];
            mesh->texcoords[i1*2] = uv1[0];
            mesh->texcoords[i1*2 + 1] = uv1[1];
            mesh->texcoords[i2*2] = uv2[0];
            mesh->texcoords[i2*2 + 1] = uv2[1];
            */
        }
        
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
    Material* mat_default = Scene_findMaterial(scene, mesh_entry.mat_name);
    mat3 inv_scale_mat, normal_mat;
    mat3_scale_inverse(inv_scale_mat, mesh_entry.scaling);
    mat3_mult(normal_mat, rotation, inv_scale_mat);
    mat3* normal_mat_ptr =  Scene_addTransform(scene, &normal_mat);
    int triangle_count = 0;
    for(int i = 0; i < num_mesh_found; i++)
    {
        Mesh* mesh = meshes[i];
        Material* mat_for_mesh;
        // NOTE: ignoreing mtl materials for now and use material defined in scene file
        /*
        Material *mat_for_mesh = Scene_findMaterial(scene, mesh->mat_name);
        if(!mat_for_mesh)
        {

        }
        */
        mat_for_mesh = mat_default;
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
                mesh_tri->mesh_ptr = mesh;
                mesh_tri->mat = mat_for_mesh;
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

                mesh_tri->mesh_ptr = mesh;
                mesh_tri->mat = mat_for_mesh;
                Object_t obj = {SMOOTH_TRIANGLE, mesh_tri};
                Scene_addObject(scene, &obj);                    
            }                
        }        
    }
    return triangle_count;
}

void linkMaterialTextures(Scene *scene)
{
    SceneMaterials* sm = (SceneMaterials*)&(scene->materials);
    SceneTextures* st = (SceneTextures*)&(scene->textures);
    for(int i = 0; i < sm->size; i++)
    {
        Material* mat = &(sm->materials[i]);
        switch(mat->mat_type)
        {
        case MATTE:
        {
            Matte* matte = (Matte*)mat->data;
            matte->diffuse = findTexture(matte->diffuse_file_name, st);
            //matte->normal = findTexture(matte->normal_file_name, st);
        } break;
        case PLASTIC:
        {
            Plastic* plastic = (Plastic*)mat->data;
            plastic->diffuse = findTexture(plastic->diffuse_file_name, st);
        } break;
        default:
        {
            
        } break;
        }
    }    
}

FILE* setupFilmAndCamera(Film* film, Camera* camera, FILE* fp)
{
    char buffer[128];
    // Setup film
    while(strcmp(buffer, "WINDOW_WIDTH") != 0)
    {
        getNextTokenInFile(buffer, fp);
    }        
    getNextTokenInFile(buffer, fp);
    film->window_width = atoi(buffer);
    getNextTokenInFile(buffer, fp); // Skip WINDOW_HEIGHT
    getNextTokenInFile(buffer, fp);
    film->window_height = atoi(buffer);
    getNextTokenInFile(buffer, fp); // Skip IMAGE_WIDTH
    getNextTokenInFile(buffer, fp);
    film->frame_res_width = atoi(buffer);
    getNextTokenInFile(buffer, fp); // Skip IMAGE_HEIGHT
    getNextTokenInFile(buffer, fp);
    film->frame_res_height = atoi(buffer);
    getNextTokenInFile(buffer, fp); // Skip FOV
    getNextTokenInFile(buffer, fp);
    film->fov = degToRad(atof(buffer));
    film->num_pixels = film->frame_res_width * film->frame_res_height;
    
    
    // Setup camera
    initPinholeCameraDefault(camera);
    vec3 cam_pos = {0.0f, 0.0f, 0.0f};
    vec3 look_point = {0.0f, 0.0f, 0.0f};
    getNextTokenInFile(buffer, fp);
    while(strcmp(buffer, "CAMERA_POS") != 0)
    {
        getNextTokenInFile(buffer, fp);
    }
    parseVec3(cam_pos, fp);
    getNextTokenInFile(buffer, fp); // Skip over LOOK_POINT
    parseVec3(look_point, fp);
    vec3 up_vec = {0.0f, 1.0f, 0.0f};
    cameraLookAt(camera, cam_pos, look_point, up_vec);

    calcFilmDimension(film, camera);
    return fp;
}


void loadSceneFile(Scene* scene, const char* scenefile)
{
    FILE* fp;
    openFile(&fp, scenefile, "r");
    char buffer[128];

    fp = setupFilmAndCamera(&(scene->film), &(scene->camera), fp);
    
    /*
      Load textures as we parse materials, but defer storing texture pointers in materials
      until all texture loadings are done.
      As we parse materials, store texture name and type, and material name into a auxillary struct.
      Once all texture loadings are done, iterate through the aux struct, get corresponding texture
      pointers and store it in the material.
     */
    //DBuffer mat_tex_pairs = DBuffer_create(MatTexNamePair);
    int num_mat = parseMaterials(scene, fp);

    //char buffer[128];
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
                OBJShape *shapes = NULL;
                OBJMaterial *materials = NULL;
                int num_mesh = 0;
                int num_mat = 0;
                char mesh_file_names[MAX_MESH][NAME_LENGTH];
                int num_file_names = 0;
                MeshEntry mesh_entry;
                parseMesh(&mesh_entry, &shapes, &materials, &num_mesh, &num_mat,
                          mesh_file_names, &num_file_names, fp);
                if(materials){free(materials);}

                for(int i = 0; i < num_mesh; i++)
                {
                    Mesh mesh;
                    Mesh_copyOBJShape(&mesh, &(shapes[i]));
                    calcTriangleNormals(&mesh);
                    calcTangentVec(&mesh);
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
        }else if(strcmp(buffer, "ENV_LIGHT") == 0)
        {            
            EnvLight* env_light = (EnvLight*)malloc(sizeof(EnvLight));
            Samples3D* samples = genHemisphereSamples(MULTIJITTERED, 1.0f);
            env_light->samples3D = samples;            
            getNextTokenInFile(buffer, fp); // Skip TYPE
            getNextTokenInFile(buffer, fp);
            if(strcmp(buffer, "CONSTANT") == 0)
            {
                env_light->type = CONSTANT;
                getNextTokenInFile(buffer, fp); // Skip COLOR
                parseColor(env_light->color, fp);
                getNextTokenInFile(buffer, fp); // Skip INTENSITY
                getNextTokenInFile(buffer, fp);
                env_light->intensity = atof(buffer);                
            }else if(strcmp(buffer, "TEXTURE") == 0)
            {
                getNextTokenInFile(buffer, fp); // Skip COLOR
                getNextTokenInFile(buffer, fp);
                // check extension
                char extension[MAX_NAME_LENGTH];
                int end_index, start_index;
                for(end_index = 0; buffer[end_index] != '\0'; end_index++){}
                for(start_index = end_index; buffer[start_index] != '.'; start_index--){}
                start_index++;
                strncpy(extension, buffer+start_index, end_index - start_index);
                if(strcmp(extension, "exr") == 0)
                {
                    env_light->type = TEXTURE;
                    Texture* env_map = &(env_light->env_map);
                    env_map->is_float = true;
                    readRgba1(buffer, (float**)&(env_map->data), &(env_map->width), &(env_map->height));
                    env_map->comp = 3;
                }else
                {
                    env_light->type = CONSTANT;
                    vec3_copy(env_light->color, WHITE);
                }
                getNextTokenInFile(buffer, fp); // Skip INTENSITY
                getNextTokenInFile(buffer, fp);
                env_light->intensity = atof(buffer);

                /*
                  0.693251 0 -0.720696 0
                  0.720696 0 0.693251 0
                  0 1 0 0
                  0 0 0 1]
                 */

                // Transform
                /*
                env_light->transform[0][0] = 0.693251;
                env_light->transform[0][1] = 0.0f;
                env_light->transform[0][2] = -0.720696;
                //env_light->transform[0][3] = 0.0f;
                env_light->transform[1][0] = 0.720696;
                env_light->transform[1][1] = 0.0f;
                env_light->transform[1][2] = -0.693251;
                //env_light->transform[1][3] = 0.0f;
                env_light->transform[2][0] = 0.0f;
                env_light->transform[2][1] = 0.0f;
                env_light->transform[2][2] = 1.0f;
                */

                mat3_rotate_y(env_light->transform, -0.76);
                //mat3_identity(env_light->transform);
            }
            if(scene->lights.env_light)
            {
                freeSamples3D(scene->lights.env_light->samples3D);                
                free(scene->lights.env_light);
            }
            SceneLights* sl = &(scene->lights);
            sl->env_light = env_light;            
            if(env_light->intensity > 0.0f)
            {
                sl->light_ptrs[sl->num_lights] = env_light;
                sl->light_types[sl->num_lights] = ENVLIGHT;                
                (sl->num_lights)++;
            }else
            {
                vec3_copy(sl->env_light->color, BLACK);
            }
        }
    }
    linkMaterialTextures(scene);

    //DBuffer_destroy(&mat_tex_pairs);
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

    //initAreaLights(sl);
    initBackgroundColor(sl);
}

int initAreaLights(Scene* scene)
{
    SceneObjects* so = &(scene->objects);
    SceneLights* sl = &(scene->lights);
    for(int i = 0; i < so->num_obj; i++)
    {
        Object_t obj = so->objects[i];
        Material* mat = getObjectMatPtr(obj);
        if((obj.type == RECTANGLE || obj.type == SPHERE || obj.type == DISK)
           && mat->mat_type == EMISSIVE)
        {
            AreaLight* area_light_ptr = (AreaLight*)malloc(sizeof(AreaLight));
            Emissive* emissive = (Emissive*)mat->data;
            area_light_ptr->intensity = emissive->intensity;
            vec3_copy(area_light_ptr->color, emissive->color);
            if(obj.type == RECTANGLE)
            {
                Rectangle* rect = (Rectangle*)(obj.ptr);
                float width = sqrt(vec3_dot(rect->width, rect->width));
                float height = sqrt(vec3_dot(rect->height, rect->height));
                area_light_ptr->flux = area_light_ptr->intensity * width * height * PI;

                area_light_ptr->pdf = 1.0f/(width * height);
                area_light_ptr->obj_ptr = rect;
                area_light_ptr->obj_type = RECTANGLE;
            }else if(obj.type == SPHERE)
            {
                Sphere* sphere = (Sphere*)(obj.ptr);
                area_light_ptr->pdf = 1.0f / (4.0f * (float)PI * sphere->radius * sphere->radius);
                area_light_ptr->obj_ptr = sphere;
                area_light_ptr->obj_type = SPHERE;                
            }else if(obj.type == DISK)
            {
                Disk* disk = (Disk*)(obj.ptr);
                area_light_ptr->pdf = 1.0f / (PI * disk->radius * disk->radius);
                area_light_ptr->obj_ptr = disk;
                area_light_ptr->obj_type = DISK;
            }
            SceneLights_addLight(area_light_ptr, AREALIGHT, sl);
        }
    }
}

int calcCausticObjectsAABB(AABB *aabb, const SceneObjects *so)
{
    int caustic_obj_count = 0;
    Mesh *cur_mesh_ptr = NULL;
    AABB cur_aabb;
    vec3_assign(cur_aabb.max, K_EPSILON, K_EPSILON, K_EPSILON);
    vec3_assign(cur_aabb.min, -K_EPSILON, -K_EPSILON, -K_EPSILON);
    for(int i = so->num_non_grid_obj; i < so->num_obj; i++)
    {
        Object_t obj = so->objects[i];
        Material *mat = getObjectMatPtr(obj);
        if(!mat)
        {
            fprintf(stderr, "Null material pointer.\n");
            continue;
        }
        if(mat->mat_type == MIRROR || mat->mat_type == TRANSPARENT)
        {
            if(obj.type != FLAT_TRIANGLE && obj.type != SMOOTH_TRIANGLE)
            {                
                if(cur_mesh_ptr)
                {
                    aabb[caustic_obj_count] = cur_aabb;
                    caustic_obj_count++;
                    cur_mesh_ptr = NULL;
                    vec3_assign(cur_aabb.max, K_EPSILON, K_EPSILON, K_EPSILON);
                    vec3_assign(cur_aabb.min, -K_EPSILON, -K_EPSILON, -K_EPSILON);
                }
                if(getObjectAABB(&cur_aabb, obj))
                {
                    aabb[caustic_obj_count] = cur_aabb;
                    caustic_obj_count++;
                }else
                {
                    fprintf(stderr, "Cannot calculate bounding box.\n");
                }
            }else
            {
                if(obj.type == FLAT_TRIANGLE)
                {
                    FlatTriangle *triangle = (FlatTriangle*)(obj.ptr);
                    if(triangle->mesh_ptr != cur_mesh_ptr && cur_mesh_ptr)
                    {
                        aabb[caustic_obj_count] = cur_aabb;
                        caustic_obj_count++;
                        cur_mesh_ptr = triangle->mesh_ptr;
                        getObjectAABB(&cur_aabb, obj);
                    }else if(!cur_mesh_ptr)
                    {
                        cur_mesh_ptr = triangle->mesh_ptr;
                        getObjectAABB(&cur_aabb, obj);
                    }else
                    {
                        AABB tri_aabb;
                        getObjectAABB(&tri_aabb, obj);
                        addToAABB(&cur_aabb, &tri_aabb);
                    }
                }else if(obj.type == SMOOTH_TRIANGLE)
                {
                    SmoothTriangle *triangle = (SmoothTriangle*)(obj.ptr);
                    if(triangle->mesh_ptr != cur_mesh_ptr && cur_mesh_ptr)
                    {
                        aabb[caustic_obj_count] = cur_aabb;
                        caustic_obj_count++;
                        cur_mesh_ptr = triangle->mesh_ptr;
                        getObjectAABB(&cur_aabb, obj);
                    }else if(!cur_mesh_ptr)
                    {
                        cur_mesh_ptr = triangle->mesh_ptr;
                        getObjectAABB(&cur_aabb, obj);
                    }else
                    {
                        AABB tri_aabb;
                        getObjectAABB(&tri_aabb, obj);
                        addToAABB(&cur_aabb, &tri_aabb);
                    }
                }
            }
        }
    }
    if(cur_mesh_ptr)
    {
        aabb[caustic_obj_count] = cur_aabb;
        caustic_obj_count++;
    }
    return caustic_obj_count;
}

void buildSceneAccel(Scene *scene)
{
    SceneObjects* so = &(scene->objects);
    double start, end;
    if(so->accel == GRID)
    {
        start = glfwGetTime();
        const int multiplier = 3;
        UniformGrid* rg = UniformGrid_create(so->objects, &(so->num_obj), so->num_non_grid_obj, multiplier);
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

void destroySceneAccel(Scene *scene)
{
    SceneObjects *so = &(scene->objects);
    if(so->accel == GRID)
    {
        UniformGrid_destroy((UniformGrid*)(so->accel_ptr));
    }else if(so->accel == BVH)
    {
        BVH_destroy((BVHNode*)(so->accel_ptr));
    }else if(so->accel == BVH4)
    {
        BVH4_destroy((BVHNode4*)(so->accel_ptr));
    }
}

// Assumes that triangles from a mesh are in a contiguous chunk
void initMeshLights(Scene* scene)
{
    MeshLight* cur_mesh_light;
    Mesh* cur_light_mesh_ptr = NULL;
    SceneObjects* so = &(scene->objects);
    SceneLights* sl = &(scene->lights);
    int inAMeshLight = 0;
    for(int i = 0; i < so->num_obj; i++)
    {
        Object_t obj = so->objects[i];
        Material* mat = getObjectMatPtr(obj);
        if(!mat)
        {
            fprintf(stderr, "Null material pointer.\n");
            continue;
        }
        if(mat->mat_type == EMISSIVE && (obj.type == FLAT_TRIANGLE || obj.type == SMOOTH_TRIANGLE))
        {
            Mesh* tri_mesh_ptr = NULL;
            Emissive* emissive = (Emissive*)mat->data;
            vec3 v0, v1, v2;
            if(obj.type == FLAT_TRIANGLE)
            {
                FlatTriangle* tri_ptr = (FlatTriangle*)(obj.ptr);
                tri_mesh_ptr = tri_ptr->mesh_ptr;
                vec3_copy(v0, tri_ptr->v0);
                vec3_copy(v1, tri_ptr->v1);
                vec3_copy(v2, tri_ptr->v2);
            }else if(obj.type == SMOOTH_TRIANGLE)
            {
                SmoothTriangle* tri_ptr = (SmoothTriangle*)(obj.ptr);
                tri_mesh_ptr = tri_ptr->mesh_ptr;
                vec3_copy(v0, tri_ptr->v0);
                vec3_copy(v1, tri_ptr->v1);
                vec3_copy(v2, tri_ptr->v2);
            }
                
            if(!inAMeshLight)
            {
                // TODO clean up
                // Encountering a new mesh light when not previously in one
                cur_light_mesh_ptr = tri_mesh_ptr;
                cur_mesh_light = (MeshLight*)malloc(sizeof(MeshLight));
                //cur_mesh_light->intensity = mat->ke;
                //vec3_copy(cur_mesh_light->color, mat->ce);
                cur_mesh_light->intensity = emissive->intensity;
                vec3_copy(cur_mesh_light->color, emissive->color);
                MeshLight_init(cur_mesh_light, tri_mesh_ptr);
                cur_mesh_light->obj_type = obj.type;
                cur_light_mesh_ptr = tri_mesh_ptr;
                inAMeshLight = 1;                
            }else if(cur_light_mesh_ptr != tri_mesh_ptr)
            {
                // Encountering a new mesh light when already in another one
                // TODO what happens two mesh lights with the same mesh are next to each other?
                MeshLight_normalizeCDF(cur_mesh_light);
                SceneLights_addLight(cur_mesh_light, MESHLIGHT, sl);
                cur_light_mesh_ptr = tri_mesh_ptr;         
                cur_mesh_light = (MeshLight*)malloc(sizeof(MeshLight));
                cur_mesh_light->obj_type = obj.type;
                cur_mesh_light->intensity = emissive->intensity;
                vec3_copy(cur_mesh_light->color, emissive->color);
                MeshLight_init(cur_mesh_light, tri_mesh_ptr);
            }
            MeshLight_addTriangle(cur_mesh_light, obj);
        }else
        {
            if(inAMeshLight && cur_mesh_light)
            {
                MeshLight_normalizeCDF(cur_mesh_light);
                SceneLights_addLight(cur_mesh_light, MESHLIGHT, sl);
                inAMeshLight = 0;
                cur_mesh_light = NULL;
                cur_light_mesh_ptr = NULL;
            }    
        }        
    }
    if(inAMeshLight && cur_mesh_light)
    {
        MeshLight_normalizeCDF(cur_mesh_light);        
        SceneLights_addLight(cur_mesh_light, MESHLIGHT, sl);
        inAMeshLight = 0;
        cur_mesh_light = NULL;
    }
}

void preprocessLights(Scene* scene)    
{
    SceneLights* sl = &(scene->lights);
    const SceneObjects* so = &(scene->objects);
    if(sl->env_light)
    {
        if(so->accel == BVH)
        {
            BVHNode* node = (BVHNode*)(so->accel_ptr);
            vec3 dist;
            vec3_sub(dist, node->aabb.max, node->aabb.min);
            //sl->env_light->world_radius = vec3_length(dist) * 0.5f;
            sl->env_light->world_radius = vec3_length(dist) * 2.0f;
        }else if(so->accel == GRID)
        {
            UniformGrid* grid = (UniformGrid*)(so->accel_ptr);
            vec3 dist;
            vec3_sub(dist, grid->aabb.max, grid->aabb.min);
            //sl->env_light->world_radius = vec3_length(dist) * 0.5f;
            sl->env_light->world_radius = vec3_length(dist) * 2.0f;
        }else if(so->accel == BVH4)
        {
            BVHNode4* node = (BVHNode4*)(so->accel_ptr);
            vec3 box_min = {K_EPSILON, K_EPSILON, K_EPSILON};
            vec3 box_max = {-HUGEVALUE, -HUGEVALUE, -HUGEVALUE};
            for(int i = 0; i < 4; i++)
            {
                box_min[0] = min(node->bbox[i + 0], box_min[0]);
                box_min[1] = min(node->bbox[i + 4], box_min[1]);
                box_min[2] = min(node->bbox[i + 8], box_min[1]);
                box_max[0] = max(node->bbox[i + 12], box_max[0]);
                box_max[1] = max(node->bbox[i + 16], box_max[0]);
                box_max[2] = max(node->bbox[i + 20], box_max[0]);
            }
            vec3 dist;
            vec3_sub(dist, box_max, box_min);
            sl->env_light->world_radius = vec3_length(dist) * 2.0f;
        }
    }
    
    float total_power = 0.0f;
    for(int i = 0; i < sl->num_lights; i++)
    {
        float power = 0.0f;
        switch(sl->light_types[i])
        {
        case ENVLIGHT:
        {
            EnvLight* env_light = (EnvLight*)sl->light_ptrs[i];
            vec3 color;
            if(env_light->type == CONSTANT)
            {
                vec3_copy(color, env_light->color);
            }else if(env_light->type == TEXTURE)
            {
                // TODO
                vec3_copy(color, WHITE);
            }

            power = (color[0] + color[1] + color[2]) / 3.0f;
            power *= env_light->intensity * env_light->world_radius;
        } break;
        case AREALIGHT:
        {
            AreaLight* area_light = (AreaLight*)sl->light_ptrs[i];
            float area = 0.0f;
            if(area_light->obj_type == SPHERE)
            {
                Sphere* sphere = (Sphere*)area_light->obj_ptr;
                area = 4.0f * PI * sphere->radius * sphere->radius;
            }else if(area_light->obj_type == RECTANGLE)
            {
                Rectangle* rect = (Rectangle*)area_light->obj_ptr;
                area = vec3_length(rect->width) * vec3_length(rect->height);
            }else if(area_light->obj_type == DISK)
            {
                Disk* disk = (Disk*)area_light->obj_ptr;
                area = PI * disk->radius * disk->radius;
            }            
            power = (area_light->color[0] * area_light->color[1] * area_light->color[2]) / 3.0f
                * area_light->intensity * area;
        } break;
        }
        sl->power[i] = power;
        total_power += power;
    }
    for(int i = 0; i < sl->num_lights; i++)
        sl->power[i] /= total_power;
}

void initScene(Scene* scene, const char* scenefile, const AccelType accel_type)
{
    // NOTE: initSceneObjects must be called after after initSceneLights if there are area lights
    initSceneLights(&(scene->lights));    
    loadSceneFile(scene, scenefile);
    SceneObjects* so = &(scene->objects);
    so->accel = accel_type;
    mvNonGridObjToStart(so);
    initAreaLights(scene);
    initMeshLights(scene);
    printf("num_obj %d\n", so->num_obj);
    printf("non bounded obj %d\n", so->num_non_grid_obj);
    // NOTE: taken out to main so that projection map can be built before the
    // mesh triangles are scrambled.
    //buildSceneAccel(scene);
}
