
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
#include "buildscene.h"
#include "intersect.h"
#include "trace.h"
#include "config.h"
#include "texture.h"
#include "noise.h"
#include "imagefile.h"
#include "projmap.h"
#include "imagestate.h"
#include "reflection.h"
#include "parallel.h"

bool SHOW_PROGRESS = true;
extern double g_traversal_time;
bool EXIT = false;
bool PAUSE = false;

double g_click_x, g_click_y;
bool g_left_clicked = false;

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(action == GLFW_PRESS)
    {
        switch(key)
        {
        case GLFW_KEY_Q:
            EXIT = true;
            break;
        case GLFW_KEY_P:
            PAUSE = !PAUSE;
            break;
        }
    }
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
    if(button == GLFW_MOUSE_BUTTON_LEFT)
    {
        if(action == GLFW_PRESS)
        {
            g_left_clicked = true;
            glfwGetCursorPos(window, &(g_click_x), &(g_click_y));
            printf("%f %f\n", g_click_x, g_click_y);
        }
    }
}

void calcProgress(const double start_time, double* end_time, int* prev_percent,
                  const int p, const ConfigParams* params)
{
    int cur_percent = (int)((float)p / (float)(params->num_samples) * 100.0f);
    double new_end_time = glfwGetTime();
    double single_iteration_time = new_end_time - *end_time;
    double whole_duration = new_end_time - start_time;
    *end_time = new_end_time;
    if(cur_percent > *prev_percent)
    {
        *prev_percent = cur_percent;
        printf("%d%%\t%f sec\t%f sec\n", cur_percent, single_iteration_time, whole_duration);
    }else
    {
        printf("%f sec\t%f sec\n", single_iteration_time, whole_duration);
    }        
}

void* threadFunc(void* vargp)
{
    ThreadData* thread_data = (ThreadData*)vargp;
    Sampler sampler;
    Sampler_create(&sampler);
    sampler.cur_sample_index = thread_data->p;
    Film* film = &(thread_data->scene->film);
    Camera* camera = &(thread_data->scene->camera);
    int start_index, end_index;
    while(JobQueue_getJob(thread_data->job_queue, &start_index, &end_index))
    {
        for(int i = start_index; i < end_index; i++)
        {
            // Calculate imageplane coordinate
            Sampler_setPixel(&sampler, i);
            vec3 color = {0.0f, 0.0f, 0.0f};
            vec2 imageplane_coord;
            vec2 sample;
            Sampler_getSample(sample, &sampler);
            calcImageCoord(imageplane_coord, film, sample, i);

            // Calculate primary ray
            Ray ray;
            Sampler_getSample(sample, &sampler);
            calcCameraRay(&ray, imageplane_coord, camera, sample);

            TraceArgs trace_args;
            trace_args.objects = &(thread_data->scene->objects);
            trace_args.lights = &(thread_data->scene->lights);
            trace_args.sampler = &sampler;

            // Tracing 
            vec3 radiance = {0.0f, 0.0f, 0.0f};
            thread_data->trace(radiance, thread_data->params->max_depth, ray, &trace_args);
            vec3_add(color, color, radiance);

            // Calculate variance
            if(thread_data->prev_num_samples + thread_data->p > 0)
            {
                vec3 current_sample_sum = {thread_data->color_buffer[i*3], thread_data->color_buffer[i*3 + 1],
                thread_data->color_buffer[i*3 + 2]};
                vec3 sample_avg;
                vec3_scale(sample_avg, current_sample_sum, 1.0 / (thread_data->prev_num_samples + thread_data->p));
                vec3 diff;
                vec3_sub(diff, color, sample_avg);
                thread_data->variance_buffer[i] = vec3_dot(diff, diff);                
            }

            // Add sample to buffer
            thread_data->color_buffer[i*3] += color[0];
            thread_data->color_buffer[i*3 + 1] += color[1];
            thread_data->color_buffer[i*3 + 2] += color[2];
            vec3_assign(color, thread_data->color_buffer[i*3], thread_data->color_buffer[i*3 +1],
                        thread_data->color_buffer[i*3 + 2]);
            vec3_scale(color, color, 1/(float)(thread_data->p+1 + thread_data->prev_num_samples));
            toneMap(color, color);

            thread_data->image[i*3] = (unsigned char)(color[0] * 255.0f);
            thread_data->image[i*3 + 1] = (unsigned char)(color[1] * 255.0f);
            thread_data->image[i*3 + 2] = (unsigned char)(color[2] * 255.0f);
        }
        glfwPollEvents();
        if(EXIT)
        {
            exit(0);
        }
    }
}

void pauseFunc(const float* buffer, const Film* film)
{
    while(PAUSE)
    {
        if(g_left_clicked)
        {
            float image_window_ratio = (float)film->frame_res_width / (float)film->window_width;
            int x = g_click_x * image_window_ratio;
            int y = g_click_y * image_window_ratio;
            int index = (y * film->frame_res_width + x) * 3;
            vec3 sample = {buffer[index], buffer[index+1], buffer[index+2]};
            printVec3WithText("", sample);
            g_left_clicked = false;
        }
        glfwPollEvents();
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

extern int g_intersect_count;
int main(int argc, char** argv)
{
    // Init GLFW
    if(glfwInit() != GL_TRUE)
    {
        fprintf(stderr, "Failed to initialize GLFW\n");
        //return NULL;
        //return -1;
    }
    
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

    // Scene data structures
    Scene scene = Scene_create();
    initScene(&scene, params.file_name, params.accel_type);
    AABB aabbs[MAX_CAUSTIC_OBJECTS];
    // TODO get rid of calcCausticObjectsAABB
    int num_aabb = calcCausticObjectsAABB(aabbs, &(scene.objects));
    buildSceneAccel(&scene);
    preprocessLights(&scene);

    GlViewport viewport;
    GLFWwindow* window = initWindow(scene.film.window_width, scene.film.window_height);    
    
    if(window)
    {
        glfwSetKeyCallback(window, keyCallback);
        glfwSetMouseButtonCallback(window, mouseButtonCallback);
        initViewport(&viewport);
    }else
    {
        SHOW_PROGRESS = false;
        params.image_save = true;
    }    
    // TODO: the image buffer could be part of film
    unsigned char* image;
    int num_pixels = scene.film.frame_res_width * scene.film.frame_res_height;
    image = (unsigned char*)calloc(num_pixels * 3, sizeof(char));


    int prev_num_samples = 0;
    float* color_buffer = NULL;
    double* variance_buffer = (double*)calloc(num_pixels, sizeof(double));
    if(using_image_state)
    {
        int width, height, size;
        readImageState(&color_buffer, &prev_num_samples, &size, &width, &height, image_state_file);
        if(width != scene.film.frame_res_width || height != scene.film.frame_res_height)
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

    srand(0);
    createGlobalSampleObject(params.num_samples, params.num_sample_sets, num_pixels);
    
    // Set trace function
    float (*trace)(vec3, int, const Ray, TraceArgs*);
    trace = getTraceFunc(params.trace_type);
    
    ThreadData thread_data;
    thread_data.prev_num_samples = prev_num_samples;
    thread_data.scene = &scene;
    thread_data.color_buffer = color_buffer;
    thread_data.variance_buffer = variance_buffer;
    thread_data.image = image;
    thread_data.params = &params;
    thread_data.trace = trace;

    int num_patches = 256;
    int num_threads = 3;
    int num_pixels_per_patch = num_pixels / num_patches;
    pthread_t threads[17];
    int patches[257];
    for(int i = 0; i < num_patches + 1; i++)
    {
        patches[i] = i * num_pixels_per_patch;
    }

    initBSDFMem(num_threads, params.max_depth+1);
    
    double start_time, end_time;
    start_time = glfwGetTime();
    end_time = start_time;
    int prev_percent = 0;

    initSampleLog(params.max_depth);
    
    JobQueue job_queue;
    JobQueue_init(&job_queue);
    thread_data.job_queue = &job_queue;
    // Rendering
    for(unsigned int p = 0; p < params.num_samples; p++)
    {
        if(PAUSE)
        {
            printf("Paused\n");
            pauseFunc(color_buffer, &(scene.film));
            printf("Unpaused\n");
        }
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

        calcProgress(start_time, &end_time, &prev_percent, p, &params);
        printf("Number of object intersection tests: %lld\n", g_intersect_count);
        if(SHOW_PROGRESS)
        {
            displayImage(window, viewport, image, scene.film.frame_res_width, scene.film.frame_res_height);
        }
    }
    end_time = glfwGetTime();
    double sec = end_time - start_time;
    printf("%f seconds.\n", sec);

    displayImage(window, viewport, image, scene.film.frame_res_width, scene.film.frame_res_height);
    printf("Traversal time = %f\n", g_traversal_time);

    // TEMP
    params.image_save = true;
    if(params.image_save)
    {
        PPM_write("output.ppm", image, num_pixels * 3, scene.film.frame_res_width, scene.film.frame_res_height);

        // Variance pass
        int indices[100];
        double top_variance[100];
        for(int i = 0; i < 10; i++){indices[i] = -1;}
        for(int i = 0; i < 10; i++){top_variance[i] = 0.0;}
        for(int i = 0; i < num_pixels; i++)
        {
            for(int j = 0; j < 100; j++)
            {
                if(variance_buffer[i] > top_variance[j])
                {
                    for(int k = 100; k > j; k--)
                    {
                        top_variance[k] = top_variance[k-1];
                        indices[k] = indices[k-1];
                    }
                    top_variance[j] = variance_buffer[i];
                    indices[j] = i;
                }
            }
        }

        for(int i = 0; i < 100; i++)
        {
            int count = 0;
            vec3 neighbor_sum = {0.0f, 0.0f, 0.0f};
            int index = indices[i] * 3;
            int current_index = index - scene.film.frame_res_width * 3;
            if(current_index > -1)
            {
                neighbor_sum[0] += color_buffer[current_index];
                neighbor_sum[1] += color_buffer[current_index];
                neighbor_sum[2] += color_buffer[current_index];
                count++;
            }

            current_index = index - 3;
            if(current_index > -1)
            {
                neighbor_sum[0] += color_buffer[current_index];
                neighbor_sum[1] += color_buffer[current_index];
                neighbor_sum[2] += color_buffer[current_index];
                count++;                
            }

            current_index = index + 3;
            if(current_index < num_pixels * 3)
            {
                neighbor_sum[0] += color_buffer[current_index];
                neighbor_sum[1] += color_buffer[current_index];
                neighbor_sum[2] += color_buffer[current_index];
                count++;                
            }

            current_index = index + scene.film.frame_res_width * 3;
            if(current_index < num_pixels * 3)
            {
                neighbor_sum[0] += color_buffer[current_index];
                neighbor_sum[1] += color_buffer[current_index];
                neighbor_sum[2] += color_buffer[current_index];
                count++;                
            }

            vec3 new_value;
            vec3_scale(new_value, neighbor_sum, 1.0f / count);
            color_buffer[index] = new_value[0];
            color_buffer[index+1] = new_value[1];
            color_buffer[index+2] = new_value[2];
        }

        for(int i = 0; i < num_pixels; i++)
        {
            int index = i * 3;
            vec3 color;
            color[0] = color_buffer[index];
            color[1] = color_buffer[index+1];
            color[2] = color_buffer[index+2];
            vec3_scale(color, color, 1.0f / (float)(thread_data.p + thread_data.prev_num_samples));
            toneMap(color, color);

            image[index] = (unsigned char)(color[0] * 255.0f);
            image[index+1] = (unsigned char)(color[1] * 255.0f);
            image[index+2] = (unsigned char)(color[2] * 255.0f);
        }

        PPM_write("output2.ppm", image, num_pixels * 3, scene.film.frame_res_width, scene.film.frame_res_height);
        
        EXIT = true;
    }

    // Save image state
    saveImageState(color_buffer, params.num_samples + prev_num_samples, scene.film.frame_res_width,
                   scene.film.frame_res_height, "savestate.is");

    // Clean up
    freeBSDFMem();
    free(color_buffer);
    double frames_per_sec = 10.0;
    double time_between_frames = 1.0 / frames_per_sec;
    double current_time, last_draw_time = 0.0;
    while(!EXIT)
    {
        current_time = glfwGetTime();
        if((current_time - last_draw_time) >= time_between_frames)
        {
            displayImage(window, viewport, image, scene.film.frame_res_width, scene.film.frame_res_height);
            last_draw_time = current_time;
        }
        glfwPollEvents();
    }
    Scene_destroy(&scene);    
    free(image);
    return 0;
}
