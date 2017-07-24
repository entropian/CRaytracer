
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <pthread.h>
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

typedef struct GlobalData_s
{
    int p;
    unsigned char* set_buffer;
    Film* film;
    Camera* camera;
    Scene* scene;
    Samples3D* h_samples;
    Photonmap* photon_map;
    Photonmap* caustic_map;
    PhotonQueryVars *query_vars;
    float* caustic_buffer;
    float* color_buffer;
    unsigned char* image;
    ConfigParams* params;
    float (*trace)(vec3, int, const Ray, TraceArgs*);
}GlobalData;
GlobalData g;

void* threadFunc(void* vargp)
{
    int* start_end = (int*)vargp;
    int start_index = start_end[0];
    int end_index = start_end[1];
    
    //for(int i = 0; i < num_pixels; i++)
    // p, set_buffer, film, camera, scene, h_samples, photon_map, query_vars
    // caustic_map, caustic_buffer, color_buffer, image
    for(int i = start_index; i < end_index; i++)
    {
        int sample_index = calcInterleavedSampleIndex(g.p, g.set_buffer[i]);
        vec3 color = {0.0f, 0.0f, 0.0f};
        vec2 imageplane_coord;
        calcImageCoord(imageplane_coord, g.film, sample_index, i);

        Ray ray;
        calcCameraRay(&ray, imageplane_coord, g.camera, sample_index);

        TraceArgs trace_args;
        //trace_args.medium_mat = medium_mat;
        trace_args.objects = &(g.scene->objects);
        trace_args.lights = &(g.scene->lights);
        trace_args.sample_index = sample_index;
        getSample3D(trace_args.h_sample, g.h_samples, sample_index);
        if(g.params->trace_type == PHOTONMAP)
        {
            trace_args.photon_map = g.photon_map;
            trace_args.query_vars = g.query_vars;
            if(g.params->caustic_map)
            {
                trace_args.caustic_map = g.caustic_map;
                vec3_assign(trace_args.caustic_rad,
                            g.caustic_buffer[i*3], g.caustic_buffer[i*3+1], g.caustic_buffer[i*3+2]);
            }
        }

        vec3 radiance = {0.0f, 0.0f, 0.0f};
        g.trace(radiance, g.params->max_depth, ray, &trace_args);
        vec3_add(color, color, radiance);

        // Planned optimizations: 
        // Irradiance caching?
        // SIMD triangle intersection for uniform grid?

        g.color_buffer[i*3] += color[0];
        g.color_buffer[i*3 + 1] += color[1];
        g.color_buffer[i*3 + 2] += color[2];
        vec3_assign(color, g.color_buffer[i*3], g.color_buffer[i*3 +1], g.color_buffer[i*3 + 2]);
        vec3_scale(color, color, 1/(float)(g.p+1));
        maxToOne(color, color);

        g.image[i*3] = (char)(color[0] * 255.0f);
        g.image[i*3 + 1] = (char)(color[1] * 255.0f);
        g.image[i*3 + 2] = (char)(color[2] * 255.0f);
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
    int num_aabb = calcCausticObjectsAABB(aabbs, &(scene.objects));
    buildSceneAccel(&scene);

    // Photon map
    const int num_photons = params.pm_config.num_photons;
    const int max_bounce = params.pm_config.photon_depth;
    const int num_caustic_photons = params.pm_config.num_caustic_photons;
    Photonmap photon_map, caustic_map;
    float* caustic_buffer = NULL;
    PhotonQueryVars query_vars;
    if(params.trace_type == PHOTONMAP)
    {
        Photonmap_init(&photon_map, num_photons, max_bounce);
        emitPhotons(&photon_map, &(scene.objects), &(scene.lights));
        Photonmap_balance(&photon_map);
        query_vars.nphotons = params.pm_config.num_estimate;
        query_vars.photon_radius = params.pm_config.photon_radius;
        query_vars.caustic_radius = params.pm_config.caustic_radius;
        query_vars.caustic = params.caustic_map;
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

    // Reminder
    destroySceneAccel(&scene);
    scene.objects.accel = GRID;
    buildSceneAccel(&scene);

    //Material *medium_mat = getMediumMatPtr(position, &(scene.objects));

    // TODO: the image buffer could be part of film
    Film film;
    film.frame_res_width = params.image_width;
    film.frame_res_height = params.image_height;
    film.num_pixels = num_pixels;
    film.fov = 70.0f / 180.0f * PI; // TODO
    Camera *camera = &(scene.camera);
    calcFilmDimension(&film, camera);
    moveSamples2D(&(film.samples), &unit_square_samples);

    // Set trace function
    float (*trace)(vec3, int, const Ray, TraceArgs*);
    //if(medium_mat)
    if(NULL)
    {
        trace = getTraceMediumFunc(params.trace_type);
    }else
    {
        trace = getTraceFunc(params.trace_type);
    }
    if(params.trace_type == PHOTONMAP && params.caustic_map)
    {
        const int num_caustic_samples = 4;
        calcCausticBuffer(caustic_buffer, camera, &film, &scene, &caustic_map, &query_vars,
                          set_buffer, num_caustic_samples);
    }

    g.set_buffer = set_buffer;
    g.film = &film;
    g.camera = camera;
    g.scene = &scene;
    g.h_samples = &h_samples;
    g.photon_map = &photon_map;
    g.caustic_map = &caustic_map;
    g.query_vars = &query_vars;
    g.caustic_buffer = caustic_buffer;
    g.color_buffer = color_buffer;
    g.image = image;
    g.params = &params;
    g.trace = trace;

    int num_patches = 4;
    int num_pixels_per_patch = num_pixels / num_patches;
    pthread_t threads[10];
    int patches[11];
    for(int i = 0; i < num_patches + 1; i++)
    {
        patches[i] = i * num_pixels_per_patch;
    }

    double start_time, end_time;
    start_time = glfwGetTime();
    end_time = start_time;
    int prev_percent = 0;

#define MULTITHREAD
#ifdef MULTITHREAD
    for(unsigned int p = 0; p < params.num_samples; p++)
    {
        g.p = p;
        // Thread stuff
        for(int i = 0; i < num_patches; i++)
        {
            pthread_create(&(threads[i]), NULL, threadFunc, &(patches[i]));
        }
        for(int i = 0; i < num_patches; i++)
        {
            pthread_join(threads[i], NULL);
        }
        
        if(SHOW_PROGRESS)
        {
            displayImage(window, viewport, image, film.frame_res_width, film.frame_res_height);
            int cur_percent = (int)((float)p / (float)(params.num_samples) * 100.0f);
            double new_end_time = glfwGetTime();
            double single_iteration_time = new_end_time - end_time;
            double whole_duration = new_end_time - start_time;
            end_time = new_end_time;
            if(cur_percent > prev_percent)
            {
                prev_percent = cur_percent;
                printf("%d%%\t%f sec\n", cur_percent, whole_duration);
            }else
            {
                printf("%f sec\t%f sec\n", single_iteration_time, whole_duration);
            }
        }
    }
#else
    for(unsigned int p = 0; p < params.num_samples; p++)
    {
        for(int i = 0; i < num_pixels; i++)
        {
            int sample_index = calcInterleavedSampleIndex(p, set_buffer[i]);
            vec3 color = {0.0f, 0.0f, 0.0f};
            vec2 imageplane_coord;
            calcImageCoord(imageplane_coord, &film, sample_index, i);

            Ray ray;
            calcCameraRay(&ray, imageplane_coord, camera, sample_index);

            TraceArgs trace_args;
            //trace_args.medium_mat = medium_mat;
            trace_args.objects = &(scene.objects);
            trace_args.lights = &(scene.lights);
            trace_args.sample_index = sample_index;
            getSample3D(trace_args.h_sample, &h_samples, sample_index);
            if(params.trace_type == PHOTONMAP)
            {
                trace_args.photon_map = &photon_map;
                trace_args.query_vars = &query_vars;
                if(params.caustic_map)
                {
                    trace_args.caustic_map = &caustic_map;
                    vec3_assign(trace_args.caustic_rad,
                                caustic_buffer[i*3], caustic_buffer[i*3+1], caustic_buffer[i*3+2]);
                }
            }

            vec3 radiance = {0.0f, 0.0f, 0.0f};
            trace(radiance, params.max_depth, ray, &trace_args);
            vec3_add(color, color, radiance);

            // Planned optimizations: 
            // Irradiance caching?
            // SIMD triangle intersection for uniform grid?

            color_buffer[i*3] += color[0];
            color_buffer[i*3 + 1] += color[1];
            color_buffer[i*3 + 2] += color[2];
            vec3_assign(color, color_buffer[i*3], color_buffer[i*3 +1], color_buffer[i*3 + 2]);
            vec3_scale(color, color, 1/(float)(p+1));
            maxToOne(color, color);

            image[i*3] = (char)(color[0] * 255.0f);
            image[i*3 + 1] = (char)(color[1] * 255.0f);
            image[i*3 + 2] = (char)(color[2] * 255.0f);
        }
        if(SHOW_PROGRESS)
        {
            displayImage(window, viewport, image, film.frame_res_width, film.frame_res_height);
            int cur_percent = (int)((float)p / (float)(params.num_samples) * 100.0f);
            double new_end_time = glfwGetTime();
            double single_iteration_time = new_end_time - end_time;
            double whole_duration = new_end_time - start_time;
            end_time = new_end_time;
            if(cur_percent > prev_percent)
            {
                prev_percent = cur_percent;
                printf("%d%%\t%f sec\n", cur_percent, whole_duration);
            }else
            {
                printf("%f sec\t%f sec\n", single_iteration_time, whole_duration);
            }
        }
    }
#endif
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
    if(params.trace_type == PHOTONMAP)
    {
        Photonmap_destroy(&photon_map);
        if(params.caustic_map)
        {
            Photonmap_destroy(&caustic_map);
            free(caustic_buffer);
        }
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
