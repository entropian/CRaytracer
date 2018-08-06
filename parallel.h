#pragma once
#include <pthread.h>
#include "scene/scene.h"

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
    float* color_buffer;
    unsigned char* image;
    ConfigParams* params;
    float (*trace)(vec3, int, const Ray, TraceArgs*);
    JobQueue* job_queue;
}ThreadData;
