#pragma once

#include <stdlib.h>
#include "util/util.h"
#include "util/vec.h"
#include <emmintrin.h>
#include <xmmintrin.h>

#define LINEAR 0
#define CUBIC 1

const int NOISE_TABLE_SIZE = 256;
const int NOISE_TABLE_MASK = NOISE_TABLE_SIZE - 1;
const int NOISE_SEED = 253;

float calcLinNoiseVal(const vec3);
float calcCubicNoiseVal(const vec3);
float calcCubicNoiseValSSE(const vec3);

typedef struct LatticeNoise_s
{
    unsigned char perm_table[NOISE_TABLE_SIZE];
    float value_table[NOISE_TABLE_SIZE];
    int num_octaves;
    float fs_min, fs_max;
    float gain, lacunarity;
}LatticeNoise;

inline unsigned char getPermIndex(const int a);
inline float getLatticeVal(const int ix, const int iy, const int iz);

/*
  nt:
  0 -- Linear interploation
  1 -- Cubic interpolation
 */
void LatticeNoise_init(const int noise_type, const int num_octaves,
                       const float gain, const float lacunarity);

// Calculate noise value with trilinear interpolation
float calcLinNoiseVal(const vec3 p);
float calcCubicNoiseValSSE(const vec3 p);
float calcCubicNoiseVal(const vec3 p);
float turbulenceNoise(const vec3 p);
float fBm(const vec3 p);
