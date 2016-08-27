#include "noise.h"
#include "util/math.h"

LatticeNoise lattice_noise;

float (*noiseFunc)(const vec3);

inline unsigned char getPermIndex(const int a)
{
    return lattice_noise.perm_table[a & NOISE_TABLE_MASK];
}

inline float getLatticeVal(const int ix, const int iy, const int iz)
{
    int index = getPermIndex(ix + getPermIndex(iy + getPermIndex(iz)));
    return lattice_noise.value_table[index];
}

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
        //noiseFunc = &calcCubicNoiseVal;
        noiseFunc = &calcCubicNoiseValSSE;
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

typedef union uSIMD_u
{
    __m128 m;
    float a[4];
}uSIMD;

float calcCubicNoiseValSSE(const vec3 p)    
{
    int ix, iy, iz;
    __m128 fx, fy;
    float fz;

    ix = (int)floor(p[0]);
    fx = _mm_set_ps1(p[0] - ix);    
    iy = (int)floor(p[1]);
    fy = _mm_set_ps1(p[1] - iy);
    iz = (int)floor(p[2]);
    fz = p[2] - iz;

    uSIMD k0, k1, k2, k3;
    __m128 out0, out1, out2, out3;

    for(int k = -1; k <= 2; k++)
    {
        for(int j = -1; j <= 2; j++)
        {
            k0.a[j+1] = getLatticeVal(ix-1, iy + j, iz + k);
            k1.a[j+1] = getLatticeVal(ix+0, iy + j, iz + k);
            k2.a[j+1] = getLatticeVal(ix+1, iy + j, iz + k);
            k3.a[j+1] = getLatticeVal(ix+2, iy + j, iz + k);            
        }
        switch(k)
        {
        case -1:
            out0 = fourKnotSplineSSE(&fx, &(k0.m), &(k1.m), &(k2.m), &(k3.m));
            break;
        case 0:
            out1 = fourKnotSplineSSE(&fx, &(k0.m), &(k1.m), &(k2.m), &(k3.m));            
            break;
        case 1:
            out2 = fourKnotSplineSSE(&fx, &(k0.m), &(k1.m), &(k2.m), &(k3.m));            
            break;
        case 2:
            out3 = fourKnotSplineSSE(&fx, &(k0.m), &(k1.m), &(k2.m), &(k3.m));            
            break;
        }
    }
    // Transpose the matrix formed by the out vectors.
    __m128 t1 = _mm_movelh_ps(out1, out0);
    __m128 t2 = _mm_movehl_ps(out0, out1);
    __m128 t3 = _mm_movelh_ps(out3, out2);
    __m128 t4 = _mm_movehl_ps(out2, out3);
    k0.m = _mm_shuffle_ps(t1, t3, _MM_SHUFFLE(0, 2, 0, 2));
    k1.m = _mm_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 3, 1, 3));
    k2.m = _mm_shuffle_ps(t2, t4, _MM_SHUFFLE(0, 2, 0, 2));
    k3.m = _mm_shuffle_ps(t2, t4, _MM_SHUFFLE(1, 3, 1, 3));                    

    uSIMD final_knots;
    final_knots.m  = fourKnotSplineSSE(&fy, &(k0.m), &(k1.m), &(k2.m), &(k3.m));
    return clamp(fourKnotSpline(fz, final_knots.a), -1.0f, 1.0f);
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

