#pragma once

#include <cstdlib>
#include <ctime>
#include <cstdio>
#include "util/vec.h"
#ifdef __APPLE__
#define abs(a) std::abs((a))
#endif
#include "util/util.h"


void setNumSamplesAndSets(const int num_samples, const int num_sets);

enum SamplesType
{
    REGULAR,
    MULTIJITTERED,
    HAMMERSLEY
};

struct Samples2D_s
{
    int num_samples;
    int num_sets;
    vec2 *samples;
};
typedef Samples2D_s Samples2D;

struct Samples3D_s
{
    int num_samples;
    int num_sets;
    vec3 *samples;
};
typedef Samples3D_s Samples3D;

Samples2D getDefaultSamples2D();
Samples3D getDefaultSamples3D();
void freeSamples2D(Samples2D *samples);
void freeSamples3D(Samples3D *samples);
void getNextSample2D(vec2 r, Samples2D *samples);
void getNextSample3D(vec3 r, Samples3D *samples);    
void genRegularSamples(Samples2D *samples);
void genMultijitteredSamples(Samples2D* samples);
void genHammersleySamples(Samples2D *samples);
void mapSamplesToDisk(Samples2D *dest_samples, const Samples2D *src_samples);
void mapSamplesToHemisphere(Samples3D* dest_samples, Samples2D * src_samples, const float e);
Samples3D* genHemisphereSamples(const SamplesType st, const float exp);
void drawSamples(unsigned char *image, Samples2D *samples,
                 const int frame_res_width, const int frame_res_height, const int num_pixels);
void drawHemisphereSamples2D(unsigned char *image, Samples3D *samples,
                             const int frame_res_width, const int frame_res_height, const int num_pixels);
void setInterleaved(bool interleaved);
void interleaveSampleSets2D(Samples2D* samples);
void interleaveSampleSets3D(Samples3D* samples);
void copySamples2D(Samples2D* dst, Samples2D* src);
void moveSamples2D(Samples2D *dst, Samples2D *src);
int calcNextSampleIndex();
int calcInterleavedSampleIndex(const int sample_index, const int set_index);
void getSample2D(vec2 r, const Samples2D* samples, const int sample_index);
void getSample3D(vec3 r, const Samples3D* samples, const int sample_index);
void createGlobalSampleObject(const int num_samples, const int num_sets, const int num_pixels);
void destroyGlobalSampleObject();


typedef struct Sampler_s
{
    int cur_sample_index;
    int cur_dimension;
    int* set_sequence;    
}Sampler;

void Sampler_create(Sampler* sampler);
void Sampler_delete(Sampler* sampler);
void Sampler_calcSetSequence(Sampler* sampler, const int a);
void Sampler_setPixel(Sampler* sampler, const int a);
void Sampler_getSample(vec2 out, Sampler* sampler);
void mapSampleToDisk(vec2 out, const vec2 in);
void mapSampleToHemisphere(vec3 out, const vec2 in);
void mapSampleWithCosPower(vec2 out, const vec2 in, const float exp);
void Sampler_getHemisphereSample(vec3 out, Sampler* sampler);
void Sampler_getCosPowerSample(vec3 out, Sampler* sampler, const float exp);
