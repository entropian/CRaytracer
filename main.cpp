
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
//#define CORNELL_BOX
#include "buildscene.h"
#include "intersect.h"
//#include "trace.h"
#include "raymarch.h"
#include "config.h"
#include "texture.h"
#include "noise.h"
#include "imagefile.h"
#include "photonmap.h"
#include "projmap.h"

#define SHOW_PROGRESS 1

extern double g_traversal_time;

bool EXIT = false;
int MAX_DEPTH = 0;

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(action == GLFW_PRESS)
    {
        switch(key)
        {
        case GLFW_KEY_Q:
            EXIT = true;
            break;
        }
    }
}

int calcCausticAABB(AABB *aabb, const SceneObjects *so)
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
        if(mat->mat_type == REFLECTIVE || mat->mat_type == TRANSPARENT)
        {
            if(obj.type != FLAT_TRIANGLE && obj.type != SMOOTH_TRIANGLE)
            {
                if(cur_mesh_ptr)
                {
                    //vec3_copy(centers[caustic_obj_count], mesh_sphere_center);
                    //radii[caustic_obj_count] = mesh_sphere_radius;
                    aabb[caustic_obj_count] = cur_aabb;
                    caustic_obj_count++;
                    cur_mesh_ptr = NULL;
                    //vec3_assign(mesh_sphere_center, 0.0f, 0.0f, 0.0f);
                    //mesh_sphere_radius = 0.0f;
                    vec3_assign(cur_aabb.max, K_EPSILON, K_EPSILON, K_EPSILON);
                    vec3_assign(cur_aabb.min, -K_EPSILON, -K_EPSILON, -K_EPSILON);
                }
                //if(calcBoundingSphere(center, &radius, obj))
                if(getObjectAABB(&cur_aabb, obj))
                {
                    //vec3_copy(centers[caustic_obj_count], center);
                    //radii[caustic_obj_count] = radius;
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
                        //vec3_copy(centers[caustic_obj_count], mesh_sphere_center);
                        //radii[caustic_obj_count] = mesh_sphere_radius;
                        aabb[caustic_obj_count] = cur_aabb;
                        caustic_obj_count++;
                        cur_mesh_ptr = triangle->mesh_ptr;
                        //calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);
                        getObjectAABB(&cur_aabb, obj);
                    }else if(!cur_mesh_ptr)
                    {
                        cur_mesh_ptr = triangle->mesh_ptr;
                        //calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);
                        getObjectAABB(&cur_aabb, obj);
                    }else
                    {
                        //addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v0);
                        //addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v1);
                        //addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v2);
                        AABB tri_aabb;
                        getObjectAABB(&tri_aabb, obj);
                        addToAABB(&cur_aabb, &tri_aabb);
                    }
                }else if(obj.type == SMOOTH_TRIANGLE)
                {
                    SmoothTriangle *triangle = (SmoothTriangle*)(obj.ptr);
                    if(triangle->mesh_ptr != cur_mesh_ptr && cur_mesh_ptr)
                    {
                        //vec3_copy(centers[caustic_obj_count], mesh_sphere_center);
                        //radii[caustic_obj_count] = mesh_sphere_radius;
                        aabb[caustic_obj_count] = cur_aabb;
                        caustic_obj_count++;
                        cur_mesh_ptr = triangle->mesh_ptr;
                        //calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);
                        getObjectAABB(&cur_aabb, obj);
                    }else if(!cur_mesh_ptr)
                    {
                        cur_mesh_ptr = triangle->mesh_ptr;
                        //calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);
                        getObjectAABB(&cur_aabb, obj);
                    }else
                    {
                        //addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v0);
                        //addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v1);
                        //addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v2);
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
    int num_pixels = params.image_width * params.image_height;
    image = (unsigned char*)calloc(num_pixels * 3, sizeof(char));

    float* color_buffer = (float*)calloc(num_pixels * 3, sizeof(float));
    //LatticeNoise_init(CUBIC, 5, 1.0f, 2.0f);

    // Samples
    setNumSamplesAndSets(params.num_samples, params.num_sample_sets);    // This sets the number of samples and sets
                                                                         // for every sample struct that follows

    setInterleaved(true);
    unsigned char* set_buffer = (unsigned char*)malloc(sizeof(unsigned char) * num_pixels);
    for(unsigned int i = 0; i < num_pixels; i++)
    {
        set_buffer[i] = (unsigned char)(rand() % params.num_sample_sets);
    }

    srand((unsigned int)time(NULL));
    Samples2D unit_square_samples = getDefaultSamples2D();
    Samples2D disk_samples = getDefaultSamples2D();
    Samples3D h_samples = getDefaultSamples3D();
    genMultijitteredSamples(&unit_square_samples);
    mapSamplesToDisk(&disk_samples, &unit_square_samples);
    mapSamplesToHemisphere(&h_samples, &disk_samples, 1);

    // Scene data structures
    Scene scene = Scene_create();
    initScene(&scene, params.file_name, params.accel_type);
    AABB aabbs[MAX_MESH];
    // NOTE: calc aabb before building accel
    int num_aabb = calcCausticAABB(aabbs, &(scene.objects));
    buildSceneAccel(&scene);

    // Photon map
    const int num_photons = params.pm_config.num_photons;
    const int max_bounce = params.pm_config.photon_depth;
    const int num_caustic_photons = params.pm_config.num_caustic_photons;
    Photonmap photon_map, caustic_map;
    bool photon_map_status = false;
    float* caustic_buffer = NULL;
    if(params.trace_type == PHOTONMAP)
    {
        photon_map_status = true;
        Photonmap_init(&photon_map, num_photons, max_bounce);
        emitPhotons(&photon_map, &(scene.objects), &(scene.lights));
        Photonmap_balance(&photon_map);
        if(params.caustic_map)
        {
            caustic_buffer = (float*)calloc(num_pixels * 3, sizeof(float));
            if(!caustic_buffer)
            {
                fprintf(stderr, "Couldn't allocate memory for caustic buffer.\n");
            }
            Photonmap_init(&caustic_map, num_caustic_photons, max_bounce);
            emitCaustics(&caustic_map, &(scene.objects), &(scene.lights), aabbs, num_aabb);
            Photonmap_balance(&caustic_map);
        }
    }
    PhotonQueryVars query_vars;
    query_vars.nphotons = params.pm_config.num_estimate;
    query_vars.photon_radius = params.pm_config.photon_radius;
    query_vars.caustic_radius = params.pm_config.caustic_radius;
    query_vars.caustic = params.caustic_map;

    // Reminder
    if(scene.objects.accel == BVH4)
    {
        BVH4_destroy((BVHNode4*)(scene.objects.accel_ptr));
    }
    scene.objects.accel = GRID;
    buildSceneAccel(&scene);

    // Camera
    Camera camera;
    initPinholeCameraDefault(&camera);
    //initThinLensCameraDefault(&camera, DEFAULT_FOCAL_LENGTH, DEFAULT_LENS_RADIUS);
    //initThinLensCameraDefault(&camera, 8.2, 1, params.num_samples, params.num_sample_sets);

#ifdef CORNELL_BOX
    //vec3 position = {278.0f, 273.0f, 800.0f};
    vec3 position = {278.0f, 600.0f, 800.0f};
    vec3 look_point = {278.0f, 273.0f, 0.0f};
    /*
      // PM test
    vec3 position = {340.0f, 500.0f, 200.0f};
    vec3 look_point = {340.0f, 500.0f, 0.0f};
    */
#else
    //vec3 position = {0.0f, 2.0f, 5.0f};
    //vec3 position = {0.0f, 50.0f, 3.0f};
    //vec3 position = {-5.0f, 40.0f, 1.0f};
    //vec3 look_point = {0.0f, 1.0f, 0.0f};
    vec3 position = {-25.0f, 7.0f, 2.0f};
    vec3 look_point = {0.0f, 7.0f, 0.0f};
#endif
    vec3 up_vec = {0.0f, 1.0f, 0.0f};
    cameraLookAt(&camera, position, look_point, up_vec);
    Material *medium_mat = getMediumMatPtr(position, &(scene.objects));

    // TODO: the image buffer could be part of film
    Film film;
    film.frame_res_width = params.image_width;
    film.frame_res_height = params.image_height;
    film.num_pixels = num_pixels;
    film.fov = 70.0f / 180.0f * PI; // TODO
    calcFilmDimension(&film, &camera);
    moveSamples2D(&(film.samples), &unit_square_samples);

    // Set trace function
    float (*trace)(vec3, int, const Ray, TraceArgs*);
    if(medium_mat)
    {
        trace = getTraceMediumFunc(params.trace_type);
    }else
    {
        trace = getTraceFunc(params.trace_type);
    }

    const int num_caustic_samples = 4;
    calcCausticBuffer(caustic_buffer, &camera, &film, &scene, &caustic_map, &query_vars,
                      set_buffer, num_caustic_samples);

//#endif
    end_time = glfwGetTime();
    double sec = end_time - start_time;
    printf("%f seconds.\n", sec);

    displayImage(window, viewport, image, film.frame_res_width, film.frame_res_height);
    printf("Traversal time = %f\n", g_traversal_time);

    if(params.image_save)
    {
        PPM_write("output.ppm", image, num_pixels * 3, params.image_width, params.image_height);
        EXIT = true;
    }

    // Clean up
    //freeSamples2D(&unit_square_samples);
    freeSamples2D(&disk_samples);
    freeSamples3D(&h_samples);
    Scene_destroy(&scene);
    Camera_destroy(&camera);
    if(params.photon_map)
    {
        Photonmap_destroy(&photon_map);
        Photonmap_destroy(&caustic_map);
        free(caustic_buffer);
    }
    free(set_buffer);

    double frames_per_sec = 10.0;
    double time_between_frames = 1.0 / frames_per_sec;
    double current_time, last_draw_time = 0.0;
    while(!EXIT)
    {
        current_time = glfwGetTime();
        if((current_time - last_draw_time) >= time_between_frames)
        {
            displayImage(window, viewport, image, film.frame_res_width, film.frame_res_height);
            last_draw_time = current_time;
        }
        glfwPollEvents();
    }
    free(image);
    return 0;
}
