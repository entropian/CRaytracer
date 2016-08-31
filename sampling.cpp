#include "sampling.h"
#include <assert.h>

const Samples2D Samples2D_default = {0, 0, NULL, NULL, 0, 0};
const Samples3D Samples3D_default = {0, 0, NULL, NULL, 0, 0};

/*
  All sampling structs share the same number of samples and sets
 */
int NUM_SAMPLES = 0;
int NUM_SAMPLE_SETS = 0;

void setNumSamplesAndSets(const int num_samples, const int num_sets)
{
    NUM_SAMPLES = num_samples;
    NUM_SAMPLE_SETS = num_sets;
}

Samples2D getDefaultSamples2D()
{
    return Samples2D_default;
}

Samples3D getDefaultSamples3D()
{
    return Samples3D_default;
}

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

// Used for rendering one pixel at a time
void getNextSample2D(vec2 r, Samples2D *samples)
{
    if(samples->count % samples->num_samples == 0)
    {
        // Randomizes which set of samples a given pixel uses
        samples->jump = (rand() % samples->num_sets) * samples->num_samples;
    }
    // sample->jump is offset for different sets
    vec2_copy(r, samples->samples[samples->jump + samples->shuffled_indices[samples->jump +
                                                                            samples->count++ % samples->num_samples]]);
}

// Used for rendering one pixel at a time
void getNextSample3D(vec3 r, Samples3D *samples)
{
    if(samples->count % samples->num_samples == 0)
    {
        samples->jump = (rand() % samples->num_sets) * samples->num_samples;
    }
    vec3_copy(r, samples->samples[samples->jump + samples->shuffled_indices[samples->jump +
                                                                            samples->count++ % samples->num_samples]]);
}

// Used for rendering the entire image one sample at a time
// Interlveaves samples across sample sets
void interleaveSampleSets2D(Samples2D* samples)
{
    vec2* new_samples = (vec2*)malloc(sizeof(vec2) * NUM_SAMPLES * NUM_SAMPLE_SETS);
    for(int i = 0; i < NUM_SAMPLES; i++)
    {
        int offset = i * NUM_SAMPLE_SETS;
        for(int j = 0; j < NUM_SAMPLE_SETS; j++)
        {
            vec2_copy(new_samples[j + offset], samples->samples[j * NUM_SAMPLES +
                                                                samples->shuffled_indices[j * NUM_SAMPLES + i]]);
        }
    }
    
    free(samples->samples);
    samples->samples = new_samples;
}

// Used for rendering the entire image one sample at a time
// Interlveaves samples across sample sets
void interleaveSampleSets3D(Samples3D* samples)
{
    vec3* new_samples = (vec3*)malloc(sizeof(vec3) * NUM_SAMPLES * NUM_SAMPLE_SETS);
    for(int i = 0; i < NUM_SAMPLES; i++)
    {
        int offset = i * NUM_SAMPLE_SETS;
        for(int j = 0; j < NUM_SAMPLE_SETS; j++)
        {
            vec3_copy(new_samples[j + offset], samples->samples[j * NUM_SAMPLES +
                                                                samples->shuffled_indices[j * NUM_SAMPLES + i]]);
        }
    }
    free(samples->samples);
    samples->samples = new_samples;
}


// Used for rendering the entire image one sample at a time
void getInterleavedSample2D(vec2 r, Samples2D* samples)
{
    if(samples->count % NUM_SAMPLE_SETS == 0)
    {
        for(int i = 0; i < NUM_SAMPLE_SETS; i++)
        {
            int random_index = rand() % NUM_SAMPLE_SETS;
            int tmp = samples->interleave_indices[i];
            samples->interleave_indices[i] = samples->interleave_indices[random_index];
            samples->interleave_indices[random_index] = tmp;
        }
    }
    
    int offset = (samples->count % (NUM_SAMPLES * NUM_SAMPLE_SETS)) / NUM_SAMPLE_SETS;
    vec2_copy(r, samples->samples[samples->interleave_indices[samples->count++ % NUM_SAMPLE_SETS] + offset]);    
    //vec2_copy(r, samples->samples[(rand() % NUM_SAMPLE_SETS) + offset]);
    //vec2_copy(r, samples->samples[(samples->count++) % (NUM_SAMPLES * NUM_SAMPLE_SETS)]);
}

// Used for rendering the entire image one sample at a time
void getInterleavedSample3D(vec3 r, Samples3D* samples)
{
    if(samples->count % NUM_SAMPLE_SETS == 0)
    {
        for(int i = 0; i < NUM_SAMPLE_SETS; i++)
        {
            int random_index = rand() % NUM_SAMPLE_SETS;
            int tmp = samples->interleave_indices[i];
            samples->interleave_indices[i] = samples->interleave_indices[random_index];
            samples->interleave_indices[random_index] = tmp;
        }
    }
    int offset = (samples->count++ % (NUM_SAMPLES * NUM_SAMPLE_SETS)) / NUM_SAMPLE_SETS;
    vec3_copy(r, samples->samples[samples->interleave_indices[samples->count++ % NUM_SAMPLE_SETS] + offset]);        
    //vec3_copy(r, samples->samples[(rand() % NUM_SAMPLE_SETS) + offset]);    
    //vec3_copy(r, samples->samples[(samples->count++) % (NUM_SAMPLES * NUM_SAMPLE_SETS)]);
}


void shuffleIndices(int **indices)
{
    *indices = (int*)malloc(sizeof(int) * NUM_SAMPLES * NUM_SAMPLE_SETS);
    for(int i = 0; i < NUM_SAMPLE_SETS; i++)
    {
        int offset = i * NUM_SAMPLES;
        for(int j = 0; j < NUM_SAMPLES; j++)
        {
            (*indices)[j + offset] = j;
        }
    }
    for(int i = 0; i < NUM_SAMPLE_SETS; i++)
    {
        int offset = i * NUM_SAMPLES;
        for(int j = 0; j < NUM_SAMPLES; j++)
        {
            int random_index = (rand() % NUM_SAMPLES) + offset;
            int tmp = (*indices)[j + offset];
            (*indices)[j + offset] = (*indices)[random_index];
            (*indices)[random_index] = tmp;
        }        
    }    
}

void prepSample2DStruct(Samples2D *samples)
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
    samples->samples = (vec2*)malloc(sizeof(vec2) * NUM_SAMPLES * NUM_SAMPLE_SETS);    
    samples->num_samples = NUM_SAMPLES;
    samples->num_sets = NUM_SAMPLE_SETS;
    shuffleIndices(&(samples->shuffled_indices));

    samples->interleave_indices = (int*)malloc(sizeof(int) * NUM_SAMPLE_SETS);
    for(int i = 0; i < NUM_SAMPLE_SETS; i++)
    {
        samples->interleave_indices[i] = i;
    }
}

void prepSample3DStruct(Samples3D *samples)
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
    samples->samples = (vec3*)malloc(sizeof(vec3) * NUM_SAMPLES * NUM_SAMPLE_SETS);    
    samples->num_samples = NUM_SAMPLES;
    samples->num_sets = NUM_SAMPLE_SETS;
    shuffleIndices(&(samples->shuffled_indices));

    samples->interleave_indices = (int*)malloc(sizeof(int) * NUM_SAMPLE_SETS);
    for(int i = 0; i < NUM_SAMPLE_SETS; i++)
    {
        samples->interleave_indices[i] = i;
    }    
}

void genRegularSamples(Samples2D *samples)
{
    int samples_per_row = (int)sqrt(NUM_SAMPLES);    
    if(samples_per_row * samples_per_row != NUM_SAMPLES)
    {
        fprintf(stderr, "num_samples must be a perfect square.\n");
        return ;
    }
    prepSample2DStruct(samples);
    int sample_index = 0;
    for(int p = 0; p < NUM_SAMPLE_SETS; p++)
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


void genMultijitteredSamples(Samples2D *samples)
{
    int samples_per_row = (int)sqrt(NUM_SAMPLES);    
    if(samples_per_row * samples_per_row != NUM_SAMPLES)
    {
        fprintf(stderr, "num_samples must be a perfect square.\n");
        return ;
    }
    prepSample2DStruct(samples);    
    int sample_index = 0;
    for(int p = 0; p < NUM_SAMPLE_SETS; p++)
    {
        for(int i = 0; i < samples_per_row; i++)
        {
            for(int j = 0; j < samples_per_row; j++)
            {
                float x = (float)rand() / (float)RAND_MAX;
                float y = (float)rand() / (float)RAND_MAX;                
                x += (float)(i + j * samples_per_row);
                x /= (float)(NUM_SAMPLES);
                y += (float)(j + i * samples_per_row);
                y /= (float)(NUM_SAMPLES);
                
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

void genHammersleySamples(Samples2D *samples)
{
    int samples_per_row = (int)sqrt(NUM_SAMPLES);    
    if(samples_per_row * samples_per_row != NUM_SAMPLES)
    {
        fprintf(stderr, "num_samples must be a perfect square.\n");
        return ;
    }
    prepSample2DStruct(samples);
    int sample_index = 0;
    for(int p = 0; p < NUM_SAMPLE_SETS; p++)
    {
        for(unsigned int i = 0; i < (unsigned)NUM_SAMPLES; i++)
        {
            samples->samples[i + p*NUM_SAMPLES][0] = (float)i / (float)NUM_SAMPLES;
            samples->samples[i + p*NUM_SAMPLES][1] = radicalInverse(i);
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
    prepSample2DStruct(dest_samples);

    float r, phi;    
    float x, y;
    unsigned int size = src_samples->num_samples * src_samples->num_sets;
    for(unsigned int i = 0; i < size; i++)
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
    prepSample3DStruct(dest_samples);
    unsigned int size = src_samples->num_samples * src_samples->num_sets;    
    for(unsigned int i = 0; i < size; i++)
    {
        float cos_phi = (float)cos(2.0f * (float)PI * src_samples->samples[i][0]);
        float sin_phi = (float)sin(2.0f * (float)PI * src_samples->samples[i][0]);
        float cos_theta = powf((1.0f - fabs(src_samples->samples[i][1])), 1.0f / (e + 1.0f));
        //float cos_theta = pow((1.0f - (src_samples->samples[i][1])), 1.0f / (e + 1.0f));
        float sin_theta = (float)sqrt(1.0f - cos_theta * cos_theta);
        float pu = sin_theta * cos_phi;
        float pv = sin_theta * sin_phi;
        float pw = cos_theta;
        vec3_assign(dest_samples->samples[i], pu, pv, pw);
    }
}

Samples3D* genHemisphereSamples(const SamplesType st, const float exp)
{
    Samples3D *samples = (Samples3D*)malloc(sizeof(Samples3D));    
    Samples2D unit_square_samples = Samples2D_default;
    Samples2D disk_samples = Samples2D_default;
    Samples3D hemisphere_samples = Samples3D_default;
    switch(st)
    {
    case REGULAR:
    {
    genRegularSamples(&unit_square_samples);
    }break;
    case MULTIJITTERED:
    {
        genMultijitteredSamples(&unit_square_samples);
    }break;
    case HAMMERSLEY:
    {
        genHammersleySamples(&unit_square_samples);
    }break;    
    }
    mapSamplesToDisk(&disk_samples, &unit_square_samples);
    mapSamplesToHemisphere(&hemisphere_samples, &disk_samples, exp);
    *samples = hemisphere_samples;
    freeSamples2D(&unit_square_samples);
    freeSamples2D(&disk_samples);    
    return samples;    
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
