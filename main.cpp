
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <pthread.h>
#include <unistd.h>
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
//#include "raymarch.h"
#include "config.h"
#include "texture.h"
#include "noise.h"
#include "imagefile.h"
//#include "photonmap.h"
#include "projmap.h"
#include "imagestate.h"
#include "reflection.h"

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

const int MAX_JOBS = 1000;
typedef struct JobQueue_s
{
    int indices[MAX_JOBS *2];
    int num_jobs;
    int head;
    pthread_mutex_t mtx;
}JobQueue;

void JobQueue_init(JobQueue* job_queue)
{
    job_queue->num_jobs = 0;
    job_queue->head = 0;
    job_queue->mtx = PTHREAD_MUTEX_INITIALIZER;
}

void JobQueue_addJob(JobQueue* job_queue, int start, int end)
{
    pthread_mutex_lock(&(job_queue->mtx));
    int queue_end = job_queue->num_jobs * 2;
    job_queue->indices[queue_end] = start;
    job_queue->indices[queue_end+1] = end;
    job_queue->num_jobs += 1;
    pthread_mutex_unlock(&(job_queue->mtx));
}

int JobQueue_getJob(JobQueue* job_queue, int* start, int* end)
{
    pthread_mutex_lock(&(job_queue->mtx));
    if(job_queue->head == job_queue->num_jobs * 2)
    {
        pthread_mutex_unlock(&(job_queue->mtx));
        return 0;
    }
    *start = job_queue->indices[job_queue->head];
    *end = job_queue->indices[job_queue->head + 1];
    job_queue->head += 2;
    pthread_mutex_unlock(&(job_queue->mtx));
    return 1;
}

typedef struct ThreadData_s
{
    int p;
    int prev_num_samples;
    unsigned char* set_buffer;
    Film* film;
    Camera* camera;
    Scene* scene;
    Samples3D* h_samples;
    //Photonmap* photon_map;
    //Photonmap* caustic_map;
    //PhotonQueryVars *query_vars;
    //float* caustic_buffer;
    float* color_buffer;
    unsigned char* image;
    ConfigParams* params;
    float (*trace)(vec3, int, const Ray, TraceArgs*);
    JobQueue* job_queue;
}ThreadData;

void* threadFunc(void* vargp)
{
    ThreadData* thread_data = (ThreadData*)vargp;
    Sampler sampler;
    Sampler_create(&sampler);
    sampler.cur_sample_index = thread_data->p;
    int start_index, end_index;
    while(JobQueue_getJob(thread_data->job_queue, &start_index, &end_index))
    {
        for(int i = start_index; i < end_index; i++)
        {
            Sampler_setPixel(&sampler, i);
            vec3 color = {0.0f, 0.0f, 0.0f};
            vec2 imageplane_coord;
            vec2 sample;
            Sampler_getSample(sample, &sampler);
            calcImageCoord(imageplane_coord, thread_data->film, sample, i);

            Ray ray;
            Sampler_getSample(sample, &sampler);
            calcCameraRay(&ray, imageplane_coord, thread_data->camera, sample);

            TraceArgs trace_args;
            trace_args.objects = &(thread_data->scene->objects);
            trace_args.lights = &(thread_data->scene->lights);
            trace_args.sampler = &sampler;
            /*
            if(thread_data->params->trace_type == PHOTONMAP)
            {
                trace_args.photon_map = thread_data->photon_map;
                trace_args.query_vars = thread_data->query_vars;
                if(thread_data->params->caustic_map)
                {
                    trace_args.caustic_map = thread_data->caustic_map;
                    vec3_assign(trace_args.caustic_rad, thread_data->caustic_buffer[i*3],
                                thread_data->caustic_buffer[i*3+1], thread_data->caustic_buffer[i*3+2]);
                }
            }
            */
            vec3 radiance = {0.0f, 0.0f, 0.0f};
            thread_data->trace(radiance, thread_data->params->max_depth, ray, &trace_args);
            vec3_add(color, color, radiance);

            thread_data->color_buffer[i*3] += color[0];
            thread_data->color_buffer[i*3 + 1] += color[1];
            thread_data->color_buffer[i*3 + 2] += color[2];
            vec3_assign(color, thread_data->color_buffer[i*3], thread_data->color_buffer[i*3 +1],
                        thread_data->color_buffer[i*3 + 2]);
            vec3_scale(color, color, 1/(float)(thread_data->p+1 + thread_data->prev_num_samples));
            maxToOne(color, color);

            thread_data->image[i*3] = (char)(color[0] * 255.0f);
            thread_data->image[i*3 + 1] = (char)(color[1] * 255.0f);
            thread_data->image[i*3 + 2] = (char)(color[2] * 255.0f);
        }
    }

}

void ppmToImageState()
{
    unsigned char* image;
    int width, height, size;
    PPM_read(&image, &size, &width, &height, "output.ppm");
    float* buffer = (float*)malloc(size * sizeof(float));
    int num_samples = 10000;
    int i;
    for(i = 0; i < size; i++)
    {
        buffer[i] = (float)(image[i]) / 255.0f * num_samples;
    }
    free(image);
    saveImageState(buffer, num_samples, width, height, "savestate.is");
    free(buffer);
}
    
int main(int argc, char** argv)
{
    bool using_image_state = false;
    char image_state_file[256];
    if(argc > 2)
    {
        printf("argv[1] %s\n", argv[1]);
        if(strcmp(argv[1], "-s") == 0)
        {
            using_image_state = true;
            stringCopy(image_state_file, 256, argv[2]);
        }
    }
    
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
    /*
    {
        float* buffer;
        int num_samples, width, height, size;
        readImageState(&buffer, &num_samples, &size, &width, &height, "savestate.is");
        if(size != num_pixels*3)
        {
            printf("wrong size\n");
        }
        for(int i = 0; i < size; i++)
        {
            image[i] = (unsigned char)(buffer[i] * 255.0f / num_samples);
        }
        while(1)
        {
            displayImage(window, viewport, image, width, height);
        }
    }
    */
    int prev_num_samples = 0;
    float* color_buffer = NULL;
    if(using_image_state)
    {
        printf("here\n");
        int width, height, size;
        readImageState(&color_buffer, &prev_num_samples, &size, &width, &height, image_state_file);
        if(width != params.image_width || height != params.image_height)
        {
            fprintf(stderr, "Wrong dimensions on image state.\n");
            free(color_buffer);
            prev_num_samples = 0;
            color_buffer = (float*)calloc(num_pixels * 3, sizeof(float));
        }
    }else
    {
        color_buffer = (float*)calloc(num_pixels * 3, sizeof(float));
    }
    //LatticeNoise_init(CUBIC, 5, 1.0f, 2.0f);

    createGlobalSampleObject(params.num_samples, params.num_sample_sets, num_pixels);
    
    // Scene data structures
    Scene scene = Scene_create();
    initScene(&scene, params.file_name, params.accel_type);
    AABB aabbs[MAX_CAUSTIC_OBJECTS];
    // NOTE: calc aabb before building accel
    // TODO fix calcCausticObjectsAABB
    int num_aabb = calcCausticObjectsAABB(aabbs, &(scene.objects));
    buildSceneAccel(&scene);

    // Photon map
    /*
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
            // TODO enable after fixing calcCausticObjectsAABB()
            emitCaustics(&caustic_map, &(scene.objects), &(scene.lights), aabbs, num_aabb);
            Photonmap_balance(&caustic_map);
        }
    }
    */

    // Reminder
    //destroySceneAccel(&scene);
    //scene.objects.accel = GRID;
    //buildSceneAccel(&scene);

    //Material *medium_mat = getMediumMatPtr(position, &(scene.objects));

    // TODO: the image buffer could be part of film
    Film film;
    film.frame_res_width = params.image_width;
    film.frame_res_height = params.image_height;
    film.num_pixels = num_pixels;
    //film.fov = 70.0f / 180.0f * PI; // TODO
    film.fov = 40.0f / 180.0f * PI; // TODO
    Camera *camera = &(scene.camera);
    calcFilmDimension(&film, camera);

    // Set trace function
    float (*trace)(vec3, int, const Ray, TraceArgs*);
    //if(medium_mat)
    if(NULL)
    {
        //trace = getTraceMediumFunc(params.trace_type);
    }else
    {
        trace = getTraceFunc(params.trace_type);
    }
    /*
    if(params.trace_type == PHOTONMAP && params.caustic_map)
    {
        const int num_caustic_samples = 4;
        calcCausticBuffer(caustic_buffer, camera, &film, &scene, &caustic_map, &query_vars, num_caustic_samples);
    }
    */
    ThreadData thread_data;
    thread_data.prev_num_samples = prev_num_samples;
    thread_data.film = &film;
    thread_data.camera = camera;
    thread_data.scene = &scene;
    //thread_data.photon_map = &photon_map;
    //thread_data.caustic_map = &caustic_map;
    //thread_data.query_vars = &query_vars;
    //thread_data.caustic_buffer = caustic_buffer;
    thread_data.color_buffer = color_buffer;
    thread_data.image = image;
    thread_data.params = &params;
    thread_data.trace = trace;

    int num_patches = 64;
    int num_threads = 4;
    int num_pixels_per_patch = num_pixels / num_patches;
    pthread_t threads[10];
    int patches[128];
    for(int i = 0; i < num_patches + 1; i++)
    {
        patches[i] = i * num_pixels_per_patch;
    }

    initBSDFMem(num_threads, params.max_depth+1);

    double start_time, end_time;
    start_time = glfwGetTime();
    end_time = start_time;
    int prev_percent = 0;

//#define MULTITHREAD
#ifdef MULTITHREAD
    JobQueue job_queue;
    JobQueue_init(&job_queue);
    thread_data.job_queue = &job_queue;
    for(unsigned int p = 0; p < params.num_samples; p++)
    {
        job_queue.num_jobs = 0;
        job_queue.head = 0;
        thread_data.p = p;
        for(int i = 0; i < num_patches; i++)
        {
            JobQueue_addJob(&job_queue, patches[i], patches[i+1]);
        }
        for(int i = 0; i < num_threads; i++)
        {
            pthread_create(&(threads[i]), NULL, threadFunc, &thread_data);
        }
        for(int i = 0; i < num_threads; i++)
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
                printf("%d%%\t%f sec\t%f sec\n", cur_percent, single_iteration_time, whole_duration);
            }else
            {
                printf("%f sec\t%f sec\n", single_iteration_time, whole_duration);
            }
        }
    }
#else
    Sampler sampler;
    Sampler_create(&sampler);
    for(unsigned int p = 0; p < params.num_samples; p++)
    {
        sampler.cur_sample_index = p;
        for(int i = 0; i < num_pixels; i++)
        {
            Sampler_setPixel(&sampler, i);
            vec3 color = {0.0f, 0.0f, 0.0f};
            vec2 sample;
            Sampler_getSample(sample, &sampler);
            vec2 imageplane_coord;
            calcImageCoord(imageplane_coord, &film, sample, i);

            Ray ray;
            Sampler_getSample(sample, &sampler);
            calcCameraRay(&ray, imageplane_coord, camera, sample);

            TraceArgs trace_args;
            trace_args.objects = &(scene.objects);
            trace_args.lights = &(scene.lights);
            trace_args.sampler = &sampler;
            /*
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
            */

            vec3 radiance = {0.0f, 0.0f, 0.0f};
            trace(radiance, params.max_depth, ray, &trace_args);
            //pathTraceOld(radiance, params.max_depth, ray, &trace_args);
            vec3_add(color, color, radiance);

            // Planned optimizations: 
            // Irradiance caching?
            // SIMD triangle intersection for uniform grid?

            color_buffer[i*3] += color[0];
            color_buffer[i*3 + 1] += color[1];
            color_buffer[i*3 + 2] += color[2];
            vec3_assign(color, color_buffer[i*3], color_buffer[i*3 +1], color_buffer[i*3 + 2]);
            vec3_scale(color, color, 1/(float)(p+1 + prev_num_samples));
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
                printf("%d%%\t%f sec\t%f sec\n", cur_percent, single_iteration_time, whole_duration);
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

    // Save image state
    saveImageState(color_buffer, params.num_samples + prev_num_samples, params.image_width,
                   params.image_height, "savestate.is");

    // Clean up
    freeBSDFMem();
    Scene_destroy(&scene);
    /*
    if(params.trace_type == PHOTONMAP)
    {
        Photonmap_destroy(&photon_map);
        if(params.caustic_map)
        {
            Photonmap_destroy(&caustic_map);
            free(caustic_buffer);
        }
    }
    */
    free(color_buffer);

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
