#pragma once

#include <cstdlib>
#include <ctime>
#include "vec.h"

/*
 TODO:
 shuffled indices that isn't crap
 Not all sampling methods needs the same randomization?
 Think about whether the disk samples should be in the same struct
*/
struct Samples_s
{
    int num_samples;
    int num_sets;
    vec2 *samples;
    vec2 *disk_samples;                // for mapping samples from unit square to unit disk
    vec3 *hemisphere_samples;
    int *shuffled_indices;
    unsigned long count;
    unsigned long disk_count;
    int jump;
    int disk_jump;                // TODO: verify if disk_count and disk_jump are necessary
} Samples_default = {0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0};
typedef Samples_s Samples;

void freeSamples(Samples *samples)
{
    if(samples->samples != NULL)
    {
        free(samples->samples);
        samples->samples = NULL;
    }
    if(samples->disk_samples != NULL)
    {
        free(samples->disk_samples);
        samples->disk_samples = NULL;
    }
    if(samples->hemisphere_samples != NULL)
    {
        free(samples->hemisphere_samples);
        samples->hemisphere_samples = NULL;
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

void getNextSample(vec2 r, Samples *samples)
{
    srand((unsigned int)time(NULL));
    if(samples->count % samples->num_samples == 0)
    {
        samples->disk_jump = (rand() % samples->num_sets) * samples->num_samples;
    }
    vec2_copy(r, samples->samples[samples->disk_jump + samples->shuffled_indices[samples->disk_jump +
                                                                            samples->count++ % samples->num_samples]]);
}

void getNextDiskSample(vec2 r, Samples *samples)
{
    if(samples->disk_samples == NULL)
    {
        fprintf(stderr, "No disk samples.\n");
        return;
    }
    srand((unsigned int)time(NULL));
    if(samples->disk_count % samples->num_samples == 0)
    {
        samples->jump = (rand() % samples->num_sets) * samples->num_samples;
    }
    vec2_copy(r, samples->disk_samples[samples->jump + samples->shuffled_indices[samples->jump +
                                                                            samples->disk_count++ % samples->num_samples]]);
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
    if(samples->disk_samples != NULL)
    {
        free(samples->disk_samples);
    }
    if(samples->hemisphere_samples != NULL)
    {
        free(samples->hemisphere_samples);
    }    
    if(samples->shuffled_indices != NULL)
    {
        free(samples->shuffled_indices);
    }
    samples->count = samples->jump = samples->disk_jump = 0;
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
float radicalInverse(unsigned int j)
{
    float x = 0.0;
    float f = 0.5;
    while(j)
    {
        x += f * (float)(j & 1);
        j = j >> 1;
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
    printf("num_samples = %d\n", num_samples);
    for(int p = 0; p < num_sets; p++)
    {
        printf("new set\n");
        for(unsigned int i = 0; i < (unsigned)num_samples; i++)
        {
            samples->samples[i + p*num_samples][0] = (float)i / (float)num_samples;
            samples->samples[i + p*num_samples][1] = radicalInverse(i);
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

void mapSamplesToDisk(Samples *samples)
{
    if(samples->samples == NULL)
    {
        fprintf(stderr, "No samples to be mapped to disk.\n");
    }    
    int size = samples->num_sets * samples->num_samples;
    float r, phi;
    if(samples->disk_samples != NULL)
    {
        free(samples->disk_samples);
    }
    samples->disk_samples = (vec2*)malloc(sizeof(vec2) * size);
    float x, y;

    for(int i = 0; i < size; i++)
    {
        // map sample point from [0,1] to [-1,1]
        x = 2.0f * samples->samples[i][0] - 1.0f;
        y = 2.0f * samples->samples[i][1] - 1.0f;
        if(x > -y)
        {
            if(x > y)
            {
                r = x;
                phi = y/x;
            }else
            {
                r = y;
                phi = 2 - x/y;
            }
        }else
        {
            if(x < y)
            {
                r = -x;
                phi = 4 + y/x;
            }else
            {
                r = -y;
                if(y != 0.0f)
                {
                    phi = 6 - x/y;
                }else
                {
                    phi = 0.0f;
                }
            }
        }
        phi *= (float)PI/4.0f;
        samples->disk_samples[i][0] = r*cos(phi);
        samples->disk_samples[i][1] = r*sin(phi);        
    }
}

void mapSamplesToHemisphere(Samples *samples, const float e)
{
    if(samples->samples == NULL)
    {
        fprintf(stderr, "No samples to be mapped to hemisphere.\n");
    }    
    int size = samples->num_sets * samples->num_samples;
    if(samples->hemisphere_samples != NULL)
    {
        free(samples->hemisphere_samples);
    }
    samples->hemisphere_samples = (vec3*)malloc(sizeof(vec3) * size);

    for(int i = 0; i < size; i++)
    {
        float cos_phi = cos(2.0f * (float)PI * samples->samples[i][0]);
        float sin_phi = sin(2.0f * (float)PI * samples->samples[i][0]);
        float cos_theta = pow((1.0f - samples->samples[i][1]), 1.0f / (e + 1.0f));
        float sin_theta = sqrt(1.0f - cos_theta * cos_theta);
        float pu = sin_theta * cos_phi;
        float pv = sin_theta * sin_phi;
        float pw = cos_theta;
        samples->hemisphere_samples[i][0] = pu;
        samples->hemisphere_samples[i][1] = pv;
        samples->hemisphere_samples[i][2] = pw;        
    }
}
