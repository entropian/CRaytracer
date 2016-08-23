#pragma once

#include <stdlib.h>
#include "util/util.h"
#include "util/vec.h"

const int NOISE_TABLE_SIZE = 256;
const int NOISE_TABLE_MASK = NOISE_TABLE_SIZE - 1;
const int NOISE_SEED = 253;


typedef struct LatticeNoise_s
{
    unsigned char perm_table[NOISE_TABLE_SIZE];
    float value_table[NOISE_TABLE_SIZE];
}LatticeNoise;

void LatticeNoise_init(LatticeNoise* ln)
{
    srand(NOISE_SEED);
    for(int i = 0; i < NOISE_TABLE_SIZE; i++)
    {
        ln->value_table[i] = 1.0f - 2.0f * ((float)rand() / (float)RAND_MAX);
        ln->perm_table[i] = i;
    }
    // Shuffle the elements of perm_table
    for(int i = 0; i < NOISE_TABLE_SIZE; i++)
    {
        int random_index = rand() % NOISE_TABLE_SIZE;
        unsigned char tmp = ln->perm_table[i];
        ln->perm_table[i] = ln->perm_table[random_index];
        ln->perm_table[random_index] = tmp;
    }    
}

inline unsigned char getPermIndex(const LatticeNoise* ln, const int a)
{
    return ln->perm_table[a & NOISE_TABLE_MASK];
}

inline float getLatticeVal(const LatticeNoise* ln, const int ix, const int iy, const int iz)
{
    int index = getPermIndex(ln, ix + getPermIndex(ln, iy + getPermIndex(ln, iz)));
    return ln->value_table[index];
}

// Calculate noise value with trilinear interpolation
float calcLinNoiseVal(const LatticeNoise* ln, const vec3 p)
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
                d[k][j][i] = getLatticeVal(ln, ix + i, iy + j, iz + k);
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

