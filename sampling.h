#pragma once

#include <cstdlib>
#include <ctime>
#include "vec.h"

/*
 TODO:
 shuffled indices
*/
struct Samples_s
{
    int num_samples;
    int num_sets;
    vec2 *samples;
    int *shuffled_indices;
    unsigned long count;
    int jump;
} Samples_default = {0, 0, NULL, NULL, 0, 0};
typedef Samples_s Samples;


void genRegularSamples(vec2 *samples, const int num_samples)
{
    int samples_per_row = (int)sqrt(num_samples);
    int sample_index;
    for(int i = 0; i < samples_per_row; i++)
    {
        for(int j = 0; j < samples_per_row; j++)
        {
            // NOTE: instead of computing a sample index, could probably just increment an index variable.
            sample_index = i*samples_per_row + j;
            samples[sample_index][0] = (i + 0.5f) / samples_per_row;
            samples[sample_index][1] = (j + 0.5f) / samples_per_row;
        }
    }
}

void getNextSample(vec2 r, Samples *samples)
{
    srand(time(NULL));
    if(samples->count % samples->num_samples == 0)
    {
        samples->jump = (rand() % samples->num_sets) * samples->num_samples;
    }
    vec2_copy(r, samples->samples[samples->jump + samples->count++ % samples->num_samples]);
}

void prepSampleStruct(Samples *samples, const int num_samples, const int num_sets)
{
    if(samples->samples != NULL)
    {
        free(samples->samples);        
    }
    if(samples->shuffled_indices != NULL)
    {
        free(samples->shuffled_indices);
    }
    samples->count = samples->jump = 0;
    samples->samples = (vec2*)malloc(sizeof(vec2) * num_samples * num_sets);
    samples->num_samples = num_samples;
    samples->num_sets = num_sets;    
}

void genRegularSamples(Samples *samples, const int num_samples, const int num_sets)
{
    int samples_per_row = (int)sqrt(num_samples);    
    if(samples_per_row * samples_per_row != num_samples)
    {
        fprintf(stderr, "num_samples must be a perfect square.\n");
        return ;
    }
    prepSampleStruct(samples, num_samples, num_sets);
    int sample_index = 0;
    for(int p = 0; p < num_sets; p++)
    {
        for(int i = 0; i < samples_per_row; i++)
        {
            for(int j = 0; j < samples_per_row; j++)
            {
                samples->samples[sample_index][0] = (j + 0.5f) / samples_per_row;
                samples->samples[sample_index][1] = (i + 0.5f) / samples_per_row;
                sample_index++;
            }
        }
    }
}

void genMultijitteredSamples(Samples *samples, const int num_samples, const int num_sets)
{
    int samples_per_row = (int)sqrt(num_samples);    
    if(samples_per_row * samples_per_row != num_samples)
    {
        fprintf(stderr, "num_samples must be a perfect square.\n");
        return ;
    }
    prepSampleStruct(samples, num_samples, num_sets);    
    int sample_index = 0;
    for(int p = 0; p < num_sets; p++)
    {
        for(int i = 0; i < samples_per_row; i++)
        {
            for(int j = 0; j < samples_per_row; j++)
            {
                float x = (float)rand() / (float)RAND_MAX;
                float y = (float)rand() / (float)RAND_MAX;                
                x += (float)(i + j * samples_per_row);
                x /= (float)(num_samples);
                y += (float)(j + i * samples_per_row);
                y /= (float)(num_samples);
                
                samples->samples[sample_index][0] = x;
                samples->samples[sample_index][1] = y;
                sample_index++;
            }
        }
    }    
}

void freeSamples(Samples *samples)
{
    if(samples->samples != NULL)
    {
        free(samples->samples);
        samples->samples = NULL;
    }
    if(samples->shuffled_indices != NULL)
    {
        free(samples->shuffled_indices);
        samples->shuffled_indices = NULL;
    }
    samples->count = samples->jump = 0;
    samples->num_sets = 0;
    samples->num_samples = 0;
}
