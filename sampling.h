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

void getNextSample(vec2 r, Samples *samples)
{
    srand((unsigned int)time(NULL));
    if(samples->count % samples->num_samples == 0)
    {
        samples->jump = (rand() % samples->num_sets) * samples->num_samples;
    }
    vec2_copy(r, samples->samples[samples->jump + samples->shuffled_indices[samples->jump +
                                                                            samples->count++ % samples->num_samples]]);
}

void shuffleIndices(Samples *samples)
{
    samples->shuffled_indices = (int*)malloc(sizeof(int) * samples->num_samples * samples->num_sets);
    for(int i = 0; i < samples->num_sets * samples->num_samples; i++)
    {
        samples->shuffled_indices[i] = -1;
    }
    srand((unsigned int)time(NULL));    
    for(int i = 0; i < samples->num_sets; i++)
    {
        int offset = i * samples->num_samples;
        for(int j = 0; j < samples->num_samples; j++)
        {
            int random_index = rand() % samples->num_samples;
            while(samples->shuffled_indices[offset + random_index] != -1)
            {
                if(random_index == samples->num_samples -1)
                {
                    random_index = 0;
                }else
                {
                    random_index++;
                }
            }
            samples->shuffled_indices[offset + random_index] = j;
        }        
    }    
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
    shuffleIndices(samples);
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

// Returns j reflected around the decimal point in binary
float radicalInverse(int j)
{
    float x = 0.0f;
    float f = 0.5f;
    while(j)
    {
        x += f * (float)(!j & 1);
        j /= 2;
        f *= 0.5f;
    }
    return x;
}

void genHammersleySamples(Samples *samples, const int num_samples, const int num_sets)
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
        for(int i = 0; i < num_samples; i++)
        {
            samples->samples[i][0] = (float)i / (float)num_samples;
            samples->samples[i][1] = radicalInverse(i);
        }
    }
}

int hashFuncMult(const int i, const int table_size)
{
    double a = 0.6180339887;
    double ia = i * a;
    double fractional_part = ia - floor(ia);
    return (int)floor(table_size * fractional_part);
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
