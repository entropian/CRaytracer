#pragma once

#include <cstdlib>
#include <ctime>
#include <cstdio>
#include "vec.h"

struct Samples2D_s
{
    int num_samples;
    int num_sets;
    vec2 *samples;
    int *shuffled_indices;
    unsigned long count;
    int jump;
} Samples2D_default = {0, 0, NULL, NULL, 0, 0};
typedef Samples2D_s Samples2D;

struct Samples3D_s
{
    int num_samples;
    int num_sets;
    vec3 *samples;
    int *shuffled_indices;
    unsigned long count;
    int jump;
} Samples3D_default = {0, 0, NULL, NULL, 0, 0};
typedef Samples3D_s Samples3D;

void freeSamples2D(Samples2D *samples)
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
    samples->count = samples->jump = samples->num_sets = samples->num_samples = 0;
}

void freeSamples3D(Samples3D *samples)
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
    samples->count = samples->jump = samples->num_sets = samples->num_samples = 0;
}

void getNextSample2D(vec2 r, Samples2D *samples)
{
    if(samples->count % samples->num_samples == 0)
    {
        samples->jump = (rand() % samples->num_sets) * samples->num_samples;
    }
    vec2_copy(r, samples->samples[samples->jump + samples->shuffled_indices[samples->jump +
                                                                            samples->count++ % samples->num_samples]]);
}

void getNextSample3D(vec3 r, Samples3D *samples)
{
    if(samples->count % samples->num_samples == 0)
    {
        samples->jump = (rand() % samples->num_sets) * samples->num_samples;
    }
    vec3_copy(r, samples->samples[samples->jump + samples->shuffled_indices[samples->jump +
                                                                            samples->count++ % samples->num_samples]]);
}

void shuffleIndices(int **indices, const int num_sets, const int num_samples)
{
    *indices = (int*)malloc(sizeof(int) * num_samples * num_sets);
    for(int i = 0; i < num_sets; i++)
    {
        int offset = i * num_samples;
        for(int j = 0; j < num_samples; j++)
        {
            (*indices)[j + offset] = j;
        }
    }
    for(int i = 0; i < num_sets; i++)
    {
        int offset = i * num_samples;
        for(int j = 0; j < num_samples; j++)
        {
            int random_index = (rand() % num_samples) + offset;
            int tmp = (*indices)[j + offset];
            (*indices)[j + offset] = (*indices)[random_index];
            (*indices)[random_index] = tmp;
        }        
    }    
}

void prepSample2DStruct(Samples2D *samples, const int num_samples, const int num_sets)
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
    shuffleIndices(&(samples->shuffled_indices), num_sets, num_samples);
}

void prepSample3DStruct(Samples3D *samples, const int num_samples, const int num_sets)
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
    samples->samples = (vec3*)malloc(sizeof(vec3) * num_samples * num_sets);    
    samples->num_samples = num_samples;
    samples->num_sets = num_sets;
    shuffleIndices(&(samples->shuffled_indices), num_sets, num_samples);
}

void genRegularSamples(Samples2D *samples, const int num_samples, const int num_sets)
{
    int samples_per_row = (int)sqrt(num_samples);    
    if(samples_per_row * samples_per_row != num_samples)
    {
        fprintf(stderr, "num_samples must be a perfect square.\n");
        return ;
    }
    prepSample2DStruct(samples, num_samples, num_sets);
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

void genMultijitteredSamples(Samples2D *samples, const int num_samples, const int num_sets)
{
    int samples_per_row = (int)sqrt(num_samples);    
    if(samples_per_row * samples_per_row != num_samples)
    {
        fprintf(stderr, "num_samples must be a perfect square.\n");
        return ;
    }
    prepSample2DStruct(samples, num_samples, num_sets);    
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
    for(int p = 0; p < samples->num_sets; p++)
    {
        int row_offset = 0;
        for(int i = 0; i < samples->num_samples; i++)
        {
            if(i % samples_per_row == 0)
            {
                row_offset = i + p * samples->num_samples;
            }
            int target = rand() % samples_per_row + row_offset;
            float tmp = samples->samples[i + p * samples->num_samples][0];
            samples->samples[i + p * samples->num_samples][0] = samples->samples[target][0];
            samples->samples[target][0] = tmp;
        }
    }
    for(int p = 0; p < samples->num_sets; p++)
    {

        for(int i = 0; i < samples->num_samples; i++)
        {
            int target = rand() % samples_per_row * samples_per_row + (i % samples_per_row)
                + p * samples->num_samples;
            float tmp = samples->samples[i + p * samples->num_samples][1];
            samples->samples[i + p * samples->num_samples][1] = samples->samples[target][1];
            samples->samples[target][1] = tmp;
        }
    }
}

/*
// Book's implementation, sort of
void genMultijitteredSamples(Samples *samples, const int num_samples, const int num_sets)
{
	// num_samples needs to be a perfect square			
	int n = (int)sqrt((float)num_samples);
	float subcell_width = 1.0f / ((float) num_samples);
	
	// fill the samples array with dummy points to allow us to use the [ ] notation when we set the 
	// initial patterns
    prepSampleStruct(samples, num_samples, num_sets);        
		
	// distribute points in the initial patterns
	
	for (int p = 0; p < num_sets; p++) 
		for (int i = 0; i < n; i++)		
			for (int j = 0; j < n; j++) {
                float rand_float = (float)rand() / (float)(RAND_MAX/subcell_width);
                samples->samples[i * n + j + p * num_samples][0] = (i * n + j) * subcell_width + rand_float;
                rand_float = (float)rand() / (float)(RAND_MAX/subcell_width);                
				samples->samples[i * n + j + p * num_samples][1] = (j * n + i) * subcell_width + rand_float;
			}
	
	// shuffle x coordinates
	
	for (int p = 0; p < num_sets; p++) 
		for (int i = 0; i < n; i++)		
			for (int j = 0; j < n; j++) {
				//int k = rand_int(j, n - 1);
                int k = rand() % (n - j) + j;
                
				float t = samples->samples[i * n + j + p * num_samples][0];
				samples->samples[i * n + j + p * num_samples][0] = samples->samples[i * n + k + p * num_samples][0];
				samples->samples[i * n + k + p * num_samples][0] = t;
			}

	// shuffle y coordinates
	
	for (int p = 0; p < num_sets; p++)
		for (int i = 0; i < n; i++)		
			for (int j = 0; j < n; j++) {
				//int k = rand_int(j, n - 1);
                int k = rand() % (n - j) + j;
				float t = samples->samples[j * n + i + p * num_samples][1];
				samples->samples[j * n + i + p * num_samples][1] = samples->samples[k * n + i + p * num_samples][1];
				samples->samples[k * n + i + p * num_samples][1] = t;
		}
}
*/

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

void genHammersleySamples(Samples2D *samples, const int num_samples, const int num_sets)
{
    int samples_per_row = (int)sqrt(num_samples);    
    if(samples_per_row * samples_per_row != num_samples)
    {
        fprintf(stderr, "num_samples must be a perfect square.\n");
        return ;
    }
    prepSample2DStruct(samples, num_samples, num_sets);
    int sample_index = 0;
    for(int p = 0; p < num_sets; p++)
    {
        for(unsigned int i = 0; i < (unsigned)num_samples; i++)
        {
            samples->samples[i + p*num_samples][0] = (float)i / (float)num_samples;
            samples->samples[i + p*num_samples][1] = radicalInverse(i);
        }
    }
}

void mapSamplesToDisk(Samples2D *dest_samples, const Samples2D *src_samples)
{
    if(src_samples->samples == NULL)
    {
        fprintf(stderr, "No samples to be mapped to disk.\n");
        return;
    }
    if(src_samples == dest_samples)
    {
        fprintf(stderr, "Samples pointers can't be the same.\n");
        return;
    }
    prepSample2DStruct(dest_samples, src_samples->num_samples, src_samples->num_sets);

    float r, phi;    
    float x, y;
    unsigned int size = src_samples->num_samples * src_samples->num_sets;
    for(int i = 0; i < size; i++)
    {
        // map sample point from [0,1] to [-1,1]
        x = 2.0f * src_samples->samples[i][0] - 1.0f;
        y = 2.0f * src_samples->samples[i][1] - 1.0f;
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
        dest_samples->samples[i][0] = r*cos(phi);
        dest_samples->samples[i][1] = r*sin(phi);        
    }
}

void mapSamplesToHemisphere(Samples3D* dest_samples, Samples2D * src_samples, const float e)
{
    if(src_samples->samples == NULL)
    {
        fprintf(stderr, "No samples to be mapped to hemisphere.\n");
        return;
    }
    prepSample3DStruct(dest_samples, src_samples->num_samples, src_samples->num_sets);
    unsigned int size = src_samples->num_samples * src_samples->num_sets;    
    for(int i = 0; i < size; i++)
    {
        float cos_phi = cos(2.0f * (float)PI * src_samples->samples[i][0]);
        float sin_phi = sin(2.0f * (float)PI * src_samples->samples[i][0]);
        float cos_theta = pow((1.0f - abs(src_samples->samples[i][1])), 1.0f / (e + 1.0f));
        //float cos_theta = pow((1.0f - (samples->disk_samples[i][1])), 1.0f / (e + 1.0f));
        float sin_theta = sqrt(1.0f - cos_theta * cos_theta);
        float pu = sin_theta * cos_phi;
        float pv = sin_theta * sin_phi;
        float pw = cos_theta;
        vec3_assign(dest_samples->samples[i], pu, pv, pw);
    }
}

void drawSamples(unsigned char *image, Samples2D *samples,
                 const int frame_res_width, const int frame_res_height, const int num_pixels)
{
    for(int i = 0; i < num_pixels; i++)
    {
        image[i*3] = (char)255;
        image[i*3 + 1] = (char)255;
        image[i*3 + 2] = (char)255;        
    }
    
    for(int i = 0; i < samples->num_samples; i++)
    {
        vec2 sample;
        getNextSample2D(sample, samples);
        int x = (int)(sample[0] * frame_res_width);
        int y = (int)((1.0f - sample[1]) * frame_res_height);        
        int index = y * frame_res_width + x;
        image[index*3] = 0;
        image[index*3 + 1] = 0;
        image[index*3 + 2] = 0;
    }
}

void drawHemisphereSamples2D(unsigned char *image, Samples3D *samples,
                             const int frame_res_width, const int frame_res_height, const int num_pixels)
{
    for(int i = 0; i < num_pixels; i++)
    {
        image[i*3] = (char)255;
        image[i*3 + 1] = (char)255;
        image[i*3 + 2] = (char)255;        
    }
    for(int i = 0; i < samples->num_samples; i++)
    {
        vec3 sample;
        getNextSample3D(sample, samples);
        int x = (int)((sample[0] * 0.5f + 0.5f) * frame_res_width);
        int y = (int)((1.0f - (sample[1] * 0.5f + 0.5f)) * frame_res_height);        
        int index = y * frame_res_width + x;

        image[index*3] = (char)(sample[2] * 255.0f);
        image[index*3 + 1] = 0;
        image[index*3 + 2] = 0;
    }    
}
