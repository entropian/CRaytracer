#include "sampling.h"
#include <assert.h>

static bool INTERLEAVED = false;
static unsigned char* set_buffer;
static int count = 0;
static int jump = 0;

const Samples2D Samples2D_default = {0, 0, NULL};
const Samples3D Samples3D_default = {0, 0, NULL};

/*
  All sampling structs share the same number of samples and sets
 */
static int NUM_SAMPLES = 0;
static int NUM_SAMPLE_SETS = 0;
static int TOTAL_SAMPLES = 0;

void setNumSamplesAndSets(const int num_samples, const int num_sets)
{
    NUM_SAMPLES = num_samples;
    NUM_SAMPLE_SETS = num_sets;
    TOTAL_SAMPLES = NUM_SAMPLES * NUM_SAMPLE_SETS;
}

void setInterleaved(bool interleaved)
{
    INTERLEAVED = interleaved;
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
    samples->num_sets = samples->num_samples = 0;
}

void freeSamples3D(Samples3D *samples)
{
    if(samples->samples != NULL)
    {
        free(samples->samples);
        samples->samples = NULL;
    }
    samples->num_sets = samples->num_samples = 0;
}

void moveSamples2D(Samples2D *dst, Samples2D *src)
{
    dst->num_samples = src->num_samples;
    dst->num_sets = src->num_sets;
    dst->samples = src->samples;
    src->num_samples = 0;
    src->num_sets = 0;
    src->samples = NULL;
}

void copySamples2D(Samples2D* dst, Samples2D* src)
{
    dst->num_samples = src->num_samples;
    dst->num_sets = src->num_sets;
    dst->samples = (vec2*)malloc(sizeof(vec2) * NUM_SAMPLES * NUM_SAMPLE_SETS);
    for(int i = 0; i < (NUM_SAMPLES * NUM_SAMPLE_SETS); i++)
    {
        vec2_copy(dst->samples[i], src->samples[i]);
    }
}

int calcNextSampleIndex()
{
    if(count % NUM_SAMPLES == 0)
    {
        jump = (rand() % NUM_SAMPLE_SETS) * NUM_SAMPLES;
    }
    int sample_index = count++ % NUM_SAMPLES;
    return jump + sample_index;
}

int calcInterleavedSampleIndex(const int sample_index, const int set_index)
{
    return sample_index * NUM_SAMPLE_SETS + set_index;
}

void getSample2D(vec2 r, const Samples2D* samples, const int sample_index)
{
    vec2_copy(r, samples->samples[sample_index % TOTAL_SAMPLES]);
}

void getSample3D(vec3 r, const Samples3D* samples, const int sample_index)
{
    vec3_copy(r, samples->samples[sample_index % TOTAL_SAMPLES]);
}

void shuffleSamples(Samples2D* samples)
{
    for(int i = 0; i < samples->num_sets; i++)
    {
        int offset = i * samples->num_samples;
        for(int j = 0; j < samples->num_samples; j++)
        {
            int random_index = (rand() % samples->num_samples) + offset;
            vec2 tmp;
            vec2_copy(tmp, samples->samples[j + offset]);
            vec2_copy(samples->samples[j + offset], samples->samples[random_index]);
            vec2_copy(samples->samples[random_index], tmp);
        }
    }
    /*
    if(INTERLEAVED)
    {
        interleaveSampleSets2D(samples);
    }
    */
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
            vec2_copy(new_samples[j + offset], samples->samples[j * NUM_SAMPLES + i]);            
        }
    }
    
    free(samples->samples);
    samples->samples = new_samples;
}

void prepSample2DStruct(Samples2D *samples)
{
    if(samples->samples != NULL)
    {
        free(samples->samples);        
    }
    samples->samples = (vec2*)malloc(sizeof(vec2) * NUM_SAMPLES * NUM_SAMPLE_SETS);    
    samples->num_samples = NUM_SAMPLES;
    samples->num_sets = NUM_SAMPLE_SETS;
}

void prepSample3DStruct(Samples3D *samples)
{
    if(samples->samples != NULL)
    {
        free(samples->samples);        
    }
    samples->samples = (vec3*)malloc(sizeof(vec3) * NUM_SAMPLES * NUM_SAMPLE_SETS);    
    samples->num_samples = NUM_SAMPLES;
    samples->num_sets = NUM_SAMPLE_SETS;
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
    shuffleSamples(samples);
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
        int rand_num = rand();
        for(int i = 0; i < samples->num_samples; i++)
        {
            if(i % samples_per_row == 0)
            {
                row_offset = i + p * samples->num_samples;
            }
            //int target = rand() % samples_per_row + row_offset;
            int target = rand_num % samples_per_row + row_offset;
            float tmp = samples->samples[i + p * samples->num_samples][0];
            samples->samples[i + p * samples->num_samples][0] = samples->samples[target][0];
            samples->samples[target][0] = tmp;
        }
    }
    for(int p = 0; p < samples->num_sets; p++)
    {
        int rand_num = rand();
        for(int i = 0; i < samples->num_samples; i++)
        {
            //int target = rand() % samples_per_row * samples_per_row + (i % samples_per_row)
            int target = rand_num % samples_per_row * samples_per_row + (i % samples_per_row)
                + p * samples->num_samples;
            float tmp = samples->samples[i + p * samples->num_samples][1];
            samples->samples[i + p * samples->num_samples][1] = samples->samples[target][1];
            samples->samples[target][1] = tmp;
        }
    }
    shuffleSamples(samples);
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
    shuffleSamples(samples);
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
        int sample_index = calcNextSampleIndex();
        //getNextSample2D(sample, samples);
        getSample2D(sample, samples, sample_index);
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
        int sample_index = calcNextSampleIndex();        
        //getNextSample3D(sample, samples);
        getSample3D(sample, samples, sample_index);
        int x = (int)((sample[0] * 0.5f + 0.5f) * frame_res_width);
        int y = (int)((1.0f - (sample[1] * 0.5f + 0.5f)) * frame_res_height);        
        int index = y * frame_res_width + x;

        image[index*3] = (char)(sample[2] * 255.0f);
        image[index*3 + 1] = 0;
        image[index*3 + 2] = 0;
    }    
}

static Samples2D global_samples;
static int** permutation_arrays = NULL;

void createGlobalSampleObject(const int num_samples, const int num_sets)
{
    setNumSamplesAndSets(num_samples, num_sets);
    genMultijitteredSamples(&global_samples);
    permutation_arrays = (int**)malloc(sizeof(int*) * num_sets);
    for(int i = 0; i < num_sets; i++)
    {
        permutation_arrays[i] = (int*)malloc(sizeof(int) * num_sets);
        for(int j = 0; j < num_sets; j++)
        {
            permutation_arrays[i][j] = j;
        }
        for(int j = 0; j < num_sets; j++)
        {
            int random_index = rand() % num_sets;
            int tmp = permutation_arrays[i][j];
            permutation_arrays[i][j] = permutation_arrays[i][random_index];
            permutation_arrays[i][random_index] = tmp;
        }
    }
}

void destroyGlobalSampleObject()
{
    freeSamples2D(&global_samples);
    if(permutation_arrays)
    {
        for(int i = 0; i < NUM_SAMPLE_SETS; i++)
        {
            free(permutation_arrays[i]);
        }
        free(permutation_arrays);
    }
}


void Sampler_create(Sampler* sampler)
{
    sampler->cur_sample_index = 0;
    sampler->cur_dimension = 0;
    sampler->set_sequence = (int*)malloc(sizeof(int) * NUM_SAMPLE_SETS);
}

void Sampler_delete(Sampler* sampler)
{
    sampler->cur_sample_index = 0;
    sampler->cur_dimension = 0;
    free(sampler->set_sequence);
    sampler->set_sequence = NULL;
}

void Sampler_calcSetSquence(Sampler* sampler, const int a)
{
    //int index = a % NUM_SAMPLE_SETS;
    int index = rand() % NUM_SAMPLE_SETS;
    for(int i = 0; i < NUM_SAMPLE_SETS; i++)
    {
        sampler->set_sequence[i] = permutation_arrays[i][index];
    }
}

void Sampler_setPixel(Sampler* sampler, const int a)
{
    sampler->cur_dimension = 0;
    Sampler_calcSetSquence(sampler, a);
}

void Sampler_getSample(vec2 out, Sampler* sampler)
{
    int offset = sampler->set_sequence[(sampler->cur_dimension)++] * NUM_SAMPLES;
    int sample_index = offset + sampler->cur_sample_index;
    vec2_copy(out, global_samples.samples[sample_index]);
}

// Need test?
void mapSampleToDisk(vec2 out, const vec2 in)
{
    float phi = 2.0f * M_PI * in[0];
    float radius = sqrtf(in[1]);
    out[0] = cosf(phi) * radius;
    out[1] = sinf(phi) * radius;
}

void mapSampleToHemisphere(vec3 out, const vec2 in)
{
    /*
    float phi = 2.0f * M_PI * in[0];
    float radius = sqrtf(in[1]);
    // radius = sin(theta)
    // z = cos(theta)
    // z = sqrt(1 - sin(theta) * sin(theta));
    out[0] = cosf(phi) * radius;
    out[1] = sinf(phi) * radius;
    out[2] = sqrtf(1.0f - radius*radius);
    */
    mapSampleToDisk(out, in);
    out[2] = sqrtf(max(0, 1.0f - out[0]*out[0] - out[1]*out[1]));
}

void mapSampleWithCosPower(vec2 out, const vec2 in, const float exp)
{
    float cos_phi = (float)cos(2.0f * (float)PI * in[0]);
    float sin_phi = (float)sin(2.0f * (float)PI * in[0]);
    float cos_theta = powf((1.0f - fabs(in[1])), 1.0f / (exp + 1.0f));
    //float cos_theta = pow((1.0f - (src_samples->samples[i][1])), 1.0f / (e + 1.0f));
    float sin_theta = (float)sqrt(1.0f - cos_theta * cos_theta);
    float pu = sin_theta * cos_phi;
    float pv = sin_theta * sin_phi;
    float pw = cos_theta;
    vec3_assign(out, pu, pv, pw);
}
