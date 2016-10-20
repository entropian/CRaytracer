
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "util/vec.h"
#include "util/ray.h"
#include "util/mat.h"
#include "util/constants.h"
#include "gl/glcode.h"
#include "sampling.h"
#include "camera.h"
#include "lights.h"
#include "scene/scene.h"
#include "shading.h"
#include "buildscene.h"
#include "intersect.h"
#include "trace.h"
#include "config.h"
#include "texture.h"
#include "noise.h"
#include "util/simd.h"
#include "bsphere.h"

#define CACHE_ALIGN __declspec(align(16))

#define SHOW_PROGRESS 1
//#define PROG
#define SIMDSHIT

bool EXIT = false;
int MAX_DEPTH = 0;

vec3 cam_position = {0.0f, 0.0f, 0.0f};

// TODO: 
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(action == GLFW_PRESS)
    {
        switch(key)
        {
        case GLFW_KEY_A:
            cam_position[0] -= 0.01f;
            break;
        case GLFW_KEY_D:
            cam_position[0] += 0.01f;
            break;
        case GLFW_KEY_W:
            cam_position[2] -= 0.01f;
            break;
        case GLFW_KEY_S:
            cam_position[2] += 0.01f;
            break;
        case GLFW_KEY_Q:
            EXIT = true;
            break;
        }
    }
}

void addToBoundingSphere(vec3 center, float *radius, vec3 point)
{
    vec3 dist_vec;
    vec3_sub(dist_vec, point, center);
    float dist = vec3_length(dist_vec);
    if(dist > *radius)
    {
        float half_diff = (dist - *radius) / dist * 0.5f;
        vec3 center_mv;
        vec3_scale(center_mv, dist_vec, half_diff);
        vec3_add(center, center, center_mv);
        *radius += (dist - *radius) * 0.5f;
    }
}

int calcCausticBoundingSpheres(vec3 *centers, float *radii, const SceneObjects * const scene_objects)
{
    int caustic_obj_count = 0;
    Mesh *cur_mesh_ptr = NULL;
    vec3 mesh_sphere_center = {0.0f, 0.0f, 0.0f};
    float mesh_sphere_radius = 0;
    for(int i = scene_objects->num_non_grid_obj; i < scene_objects->num_obj; i++)
    {
        if(i == 8886)
        {
            printf("here\n");
        }
        Object_t obj = scene_objects->objects[i];
        Material *mat = getObjectMatPtr(obj);        
        if(!mat)
        {
            fprintf(stderr, "Null material pointer.\n");
            continue;
        }
        if(mat->mat_type == REFLECTIVE || mat->mat_type == TRANSPARENT)
        {
            if(obj.type != FLAT_TRIANGLE && obj.type != SMOOTH_TRIANGLE)
            {
                if(cur_mesh_ptr)
                {
                    vec3_copy(centers[caustic_obj_count], mesh_sphere_center);
                    radii[caustic_obj_count] = mesh_sphere_radius;                        
                    caustic_obj_count++;                    
                    cur_mesh_ptr = NULL;
                    vec3_assign(mesh_sphere_center, 0.0f, 0.0f, 0.0f);
                    mesh_sphere_radius = 0.0f;
                }
                vec3 center;
                float radius;
                if(calcBoundingSphere(center, &radius, obj))
                {
                    vec3_copy(centers[caustic_obj_count], center);
                    radii[caustic_obj_count] = radius;
                    caustic_obj_count++;
                }else
                {
                    fprintf(stderr, "Cannot calculate bounding sphere.\n");
                }
            }else
            {                
                if(obj.type == FLAT_TRIANGLE)
                {
                    FlatTriangle *triangle = (FlatTriangle*)(obj.ptr);
                    if(triangle->mesh_ptr != cur_mesh_ptr && cur_mesh_ptr)
                    {
                        vec3_copy(centers[caustic_obj_count], mesh_sphere_center);
                        radii[caustic_obj_count] = mesh_sphere_radius;                        
                        caustic_obj_count++;
                        cur_mesh_ptr = triangle->mesh_ptr;
                        calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);
                    }else if(!cur_mesh_ptr)
                    {
                        cur_mesh_ptr = triangle->mesh_ptr;                        
                        calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);                        
                    }else
                    {                        
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v0);
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v1);
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v2);                        
                    }
                }else if(obj.type == SMOOTH_TRIANGLE)
                {
                    SmoothTriangle *triangle = (SmoothTriangle*)(obj.ptr);
                    if(triangle->mesh_ptr != cur_mesh_ptr && cur_mesh_ptr)
                    {
                        vec3_copy(centers[caustic_obj_count], mesh_sphere_center);
                        radii[caustic_obj_count] = mesh_sphere_radius;                        
                        caustic_obj_count++;
                        cur_mesh_ptr = triangle->mesh_ptr;
                        calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);
                    }else if(!cur_mesh_ptr)
                    {
                        cur_mesh_ptr = triangle->mesh_ptr;
                        calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);
                    }else
                    {
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v0);
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v1);
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v2);
                    }
                }
            }
        }
    }
    if(cur_mesh_ptr)
    {
        vec3_copy(centers[caustic_obj_count], mesh_sphere_center);
        radii[caustic_obj_count] = mesh_sphere_radius;                        
        caustic_obj_count++;                    
    }
    return caustic_obj_count;
}

void projCircleToMap(int *map, const int theta_row, const int phi_column,
                     const vec3 center, const float radius, const vec3 light_pos)
{
    vec3 sphere_to_light;
    vec3_sub(sphere_to_light, center, light_pos);
    vec3 center_dir;
    float dist;
    vec3_normalize(center_dir, sphere_to_light);
    dist = vec3_length(sphere_to_light);
    vec3 perp_vec;
    vec3_cross(perp_vec, center_dir, JITTERED_UP);
    vec3_normalize(perp_vec, perp_vec);
    vec3_scale(perp_vec, perp_vec, radius);
    vec3 side_vec;
    vec3_add(side_vec, sphere_to_light, perp_vec);
    vec3_normalize(side_vec, side_vec);
    float angle = acos(vec3_dot(center_dir, side_vec));

    // spherical coordinates
    // atan2 output range is -PI, PI 
    float theta_per_cell = (float)(PI / (double)theta_row);
    float phi_per_cell = (float)(2.0 * PI / (double)phi_column);
    float theta = PI * 0.5 - asin(center_dir[1]);
    float phi = atan2(center_dir[2], center_dir[0]);
    int theta_index = theta / theta_per_cell;
    int phi_index = (phi + PI) / phi_per_cell;
    int theta_half_span = angle / theta_per_cell + 1;

    for(int i = theta_index - theta_half_span; i < theta_index + theta_half_span; i++)
    {
        int theta_index = i * phi_column;
        int phi_offset = 0;
        if(i < 0)
        {
            theta_index = (abs(i) - 1) * phi_column;
            phi_offset = phi_column / 2;
        }
        if(i >= theta_row)
        {
            theta_index = (theta_row - (i - theta_row) - 1) * phi_column;
            phi_offset = phi_column / 2;
        }
        float ang_from_equator = abs(i - theta_row / 2) * theta_per_cell;
        float cur_phi_per_cell = (float)cos(ang_from_equator + (theta_per_cell * 0.5f)) * phi_per_cell;
        cur_phi_per_cell = fabs(cur_phi_per_cell);
        int phi_half_span = angle / cur_phi_per_cell + 1;
        if(phi_half_span * 2 >= phi_column)
        {
            for(int j = 0; j < phi_column; j++)
            {
                map[theta_index + j] = 1;
            }
            continue;
        }
        for(int j = phi_index - phi_half_span; j < phi_index + phi_half_span; j++)
        {
            int phi_index = j;
            if(j < 0)
            {
                phi_index = phi_column + j;
            }
            if(j >= phi_column)
            {
                phi_index = j - phi_column;
            }
            phi_index = (phi_index + phi_offset) % phi_column;
            int index = theta_index + phi_index;
            assert(index > -1 && index < theta_row * phi_column);
            map[index] = 1;
        }
    } 
}


int main()
{
    ConfigParams params;
    parseConfigFile(&params);
    MAX_DEPTH = params.max_depth;
    
    GLFWwindow* window = initWindow(params.window_width, params.window_height);
    glfwSetKeyCallback(window, keyCallback);
    
    GlViewport viewport;
    initViewport(&viewport);
    
    unsigned char *image;
    int frame_res_width = params.image_width, frame_res_height = params.image_height;
    int num_pixels = frame_res_width * frame_res_height;
    image = (unsigned char*)calloc(num_pixels * 3, sizeof(char));    
    
    // Samples
    setNumSamplesAndSets(params.num_samples, params.num_sample_sets);    // This sets the number of samples and sets for every 
                                                                         // sample struct that follows


#ifdef PROG
    setInterleaved(true);
    unsigned char* set_buffer = (unsigned char*)malloc(sizeof(unsigned char) * num_pixels);
    for(unsigned int i = 0; i < num_pixels; i++)
    {
        set_buffer[i] = (unsigned char)(rand() % params.num_sample_sets);
    }    
#endif
    srand((unsigned int)time(NULL));    
    Samples2D unit_square_samples = getDefaultSamples2D();
    Samples2D disk_samples = getDefaultSamples2D();
    Samples3D h_samples = getDefaultSamples3D();
    genMultijitteredSamples(&unit_square_samples);
    mapSamplesToDisk(&disk_samples, &unit_square_samples);
    mapSamplesToHemisphere(&h_samples, &disk_samples, 1);

    // Scene data structures
    const int MAX_BSPHERES = 100;
    Scene scene = Scene_create();
    initScene(&scene, params.file_name, params.accel_type);
    buildSceneAccel(&scene);
    vec3 bsphere_centers[MAX_BSPHERES];
    float bsphere_radii[MAX_BSPHERES];
    
    int num_bsphere = calcCausticBoundingSpheres(bsphere_centers, bsphere_radii, &(scene.objects));

    // Camera
    Camera camera;
    initPinholeCameraDefault(&camera);
    
    vec3 position = {0.0f, 2.0f, 5.0f};
    vec3 look_point = {0.0f, 2.0f, 0.0f};
    vec3 up_vec = {0.0f, 1.0f, 0.0f};
    cameraLookAt(&camera, position, look_point, up_vec);

    // TODO: throw these in a struct for camera scaling
    float fov = 70.0f / 180.0f * PI;
    float frame_length = 2.0f * (sin(fov/2.0f) * camera.focal_pt_dist);
    float frame_height = frame_length * (float)(frame_res_height)/(float)(frame_res_width);    
    float pixel_length = frame_length/(float)(frame_res_width);


    vec3 light_pos = {0.0f, 0.0f, 0.0f};
    // sphere: phi = 2*PI theta = PI   
    const int theta_row = 26;
    const int phi_column = 52;
    int grid[theta_row * phi_column];
    for(int i = 0; i < theta_row * phi_column; i++)
    {
        grid[i] = 0;
    }
    Sphere sphere;
    vec3_assign(sphere.center, 0.0f, 5.0f, 0.01f);
    sphere.radius = 5.0f;    
    projCircleToMap(grid, theta_row, phi_column, sphere.center, sphere.radius, light_pos);

    for(int i = 0; i < theta_row; i++)
    {
        for(int j = 0; j < phi_column; j++)
        { 
            printf("%d ", grid[i*phi_column + j]);
        }
        printf("\n");
    }

//#endif 
    //displayImage(window, viewport, image, frame_res_width, frame_res_height);

    // Clean up
    freeSamples2D(&unit_square_samples);
    freeSamples2D(&disk_samples);
    freeSamples3D(&h_samples);
    Scene_destroy(&scene);
    Camera_destroy(&camera);

    double frames_per_sec = 10.0;
    double time_between_frames = 1.0 / frames_per_sec;
    double current_time, last_draw_time = 0.0;
    while(!EXIT)
    {
        current_time = glfwGetTime();
        if((current_time - last_draw_time) >= time_between_frames)
        {
            last_draw_time = current_time;
        }
        glfwPollEvents();
    }
    //free(image);    
    return 0;
}
