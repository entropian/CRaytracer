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

float (*noiseFunc)(const vec3);

float calcLinNoiseVal(const vec3);
float calcCubicNoiseVal(const vec3);

typedef struct LatticeNoise_s
{
    unsigned char perm_table[NOISE_TABLE_SIZE];
    float value_table[NOISE_TABLE_SIZE];
    int num_octaves;
    float fs_min, fs_max;
    float gain, lacunarity;
}LatticeNoise;

LatticeNoise lattice_noise;

/*
  nt:
  0 -- Linear interploation
  1 -- Cubic interpolation
 */
void LatticeNoise_init(const int noise_type, const int num_octaves,
                       const float gain, const float lacunarity)
{       
    if(num_octaves < 1)
    {
        fprintf(stderr, "Number of octaves cannot be less than 1. Defaulting to 3.\n");
        lattice_noise.num_octaves = 3;
    }else
    {        
        lattice_noise.num_octaves = num_octaves;
    }

    if(gain < 0)
    {
        fprintf(stderr, "Noise gain must not be negative. Defaulting to 0.5\n");
        lattice_noise.gain = 0.5f;
    }else
    {
        lattice_noise.gain = gain;
    }

    if(lacunarity <= 0)
    {
        fprintf(stderr, "Noise lacunarity must be greater than 0. Defaulting to 2.0.\n");
        lattice_noise.lacunarity = 2.0f;
    }else
    {
        lattice_noise.lacunarity = lacunarity;
    }

    // Compute fractal sum bounds
    lattice_noise.fs_max = 0.0f;
    for(int i = 0; i < lattice_noise.num_octaves; i++)
    {
        lattice_noise.fs_max += powf(lattice_noise.gain, i);        
    }
    lattice_noise.fs_min = -lattice_noise.fs_max;

    int nt;
    if(noise_type < 0 || noise_type > 1)
    {
        fprintf(stderr, "Invalid noise type. Defaulting to linear.\n");
        nt = 0;
    }else
    {
        nt = noise_type;
    }
    
    switch(nt)
    {
    case LINEAR:
        noiseFunc = &calcLinNoiseVal;
        break;
    case CUBIC:
        noiseFunc = &calcCubicNoiseVal;
        break;
    }
    
    srand(NOISE_SEED);
    for(int i = 0; i < NOISE_TABLE_SIZE; i++)
    {
        lattice_noise.value_table[i] = 1.0f - 2.0f * ((float)rand() / (float)RAND_MAX);
        lattice_noise.perm_table[i] = i;
    }
    // Shuffle the elements of perm_table
    for(int i = 0; i < NOISE_TABLE_SIZE; i++)
    {
        int random_index = rand() % NOISE_TABLE_SIZE;
        unsigned char tmp = lattice_noise.perm_table[i];
        lattice_noise.perm_table[i] = lattice_noise.perm_table[random_index];
        lattice_noise.perm_table[random_index] = tmp;
    }    
}

inline unsigned char getPermIndex(const int a)
{
    return lattice_noise.perm_table[a & NOISE_TABLE_MASK];
}

inline float getLatticeVal(const int ix, const int iy, const int iz)
{
    int index = getPermIndex(ix + getPermIndex(iy + getPermIndex(iz)));
    return lattice_noise.value_table[index];
}

// Calculate noise value with trilinear interpolation
float calcLinNoiseVal(const vec3 p)
{
    int ix, iy, iz;
    float fx, fy, fz;    
    float d[2][2][2];
    float x0, x1, x2, x3, y0, y1, z0;

    ix = (int)floor(p[0]);
    fx = p[0] - ix;
    iy = (int)floor(p[1]);
    fy = p[1] - iy;
    iz = (int)floor(p[2]);
    fz = p[2] - iz;

    for(int k = 0; k <= 1; k++)
    {
        for(int j = 0; j <= 1; j++)
        {
            for(int i = 0; i <= 1; i++)
            {
                d[k][j][i] = getLatticeVal(ix + i, iy + j, iz + k);
            }
        }
    }

    x0 = lerp(fx, d[0][0][0], d[0][0][1]);
    x1 = lerp(fx, d[0][1][0], d[0][1][1]);
    x2 = lerp(fx, d[1][0][0], d[1][0][1]);
    x3 = lerp(fx, d[1][1][0], d[1][1][1]);
    y0 = lerp(fy, x0, x1);
    y1 = lerp(fy, x2, x3);
    z0 = lerp(fz, y0, y1);

    return z0;
}

float calcCubicNoiseVal(const vec3 p)
{
    int ix, iy, iz;
    float fx, fy, fz;
    float xknots[4], yknots[4], zknots[4];

    ix = (int)floor(p[0]);
    fx = p[0] - ix;
    iy = (int)floor(p[1]);
    fy = p[1] - iy;
    iz = (int)floor(p[2]);
    fz = p[2] - iz;

    for(int k = -1; k <= 2; k++)
    {
        for(int j = -1; j <= 2; j++)
        {
            for(int i = -1; i <= 2; i++)
            {
                xknots[i+1] = getLatticeVal(ix + i, iy + j, iz + k);
            }
            yknots[j+1] = fourKnotSpline(fx, xknots);
        }
        zknots[k+1] = fourKnotSpline(fy, yknots);        
    }
    return clamp(fourKnotSpline(fz, zknots), -1.0f, 1.0f);
}

float turbulenceNoise(const vec3 p)
{
    float amplitude = 1.0f;
    float frequency = 1.0f;
    float turbulence = 0.0f;

    vec3 p_freq;
    for(int i = 0; i < lattice_noise.num_octaves; i++)
    {
        vec3_scale(p_freq, p, frequency);            
        turbulence += amplitude * fabs(noiseFunc(p_freq));
        amplitude *= 0.5f;
        frequency *= 2.0f;
    }
    return turbulence /= lattice_noise.fs_max;    
}

float fBm(const vec3 p)
{
    float amplitude = 1.0f;
    float frequency = 1.0f;
    float fBm = 0.0f;

    vec3 p_freq;    
    for(int i = 0; i < lattice_noise.num_octaves; i++)
    {
        vec3_scale(p_freq, p, frequency);                    
        fBm += amplitude * noiseFunc(p_freq);
        amplitude *= lattice_noise.gain;
        frequency *= lattice_noise.lacunarity;        
    }
    return (fBm - lattice_noise.fs_min) / (lattice_noise.fs_max - lattice_noise.fs_min);        
}
