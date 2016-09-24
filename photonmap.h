#include <stdlib.h>
#include <math.h>
#include "util/vec.h"
#include "scene/scenedata.h"
#include "util/ray.h"
#include "objloader/dbuffer.h"
#include "aabb.h"
#include "shapes/shapes.h"

extern float intersectTest(ShadeRec* sr, const SceneObjects* so, const Ray ray);

typedef struct Photon_s
{
    vec3 position;
    short plane;
    unsigned char theta, phi;
    vec3 power;
}Photon;

typedef struct Photonmap_s
{
    Photon* photons;
    int stored_photons;
    int half_stored_photons;
    int max_photons;
    int prev_scale;
    float costheta[256];
    float sintheta[256];
    float cosphi[256];
    float sinphi[256];
    AABB bbox;
}Photonmap;

void initPhotonmap(Photonmap* photon_map, const int max_photons)
{
    photon_map->stored_photons = 0;
    photon_map->prev_scale = 1;
    photon_map->max_photons = max_photons;
    photon_map->photons = (Photon*)malloc(sizeof(Photon) * max_photons);
    if(!photon_map->photons)
    {
        fprintf(stderr, "Not enough memory initializing photon map.\n");
        exit(-1);
    }

    vec3_assign(photon_map->bbox.min, 1e8f, 1e8f, 1e8f);
    vec3_assign(photon_map->bbox.max, -1e8f, -1e8f, -1e8f);
    for(int i = 0; i < 256; i++)
    {
        double angle = (double)i * (1.0/256.0) * PI;
        photon_map->costheta[i] = cos(angle);
        photon_map->sintheta[i] = sin(angle);
        photon_map->cosphi[i] = cos(2.0*angle);
        photon_map->sinphi[i] = sin(2.0*angle);
    }
}

void Photonmap_destroy(Photonmap* photon_map)
{
    if(photon_map->photons && photon_map->max_photons > 0)
    {
        free(photon_map->photons);
    }
    photon_map->max_photons = 0;
    photon_map->stored_photons = 0;
    photon_map->half_stored_photons = 0;
}

int emitPhotons(Photonmap* photon_map, const SceneObjects *so, const SceneLights *sl)
{    
    int point_light_count = 0;
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == POINTLIGHT)
        {
            point_light_count++;
        }
    }
    if(point_light_count == 0)
    {
        fprintf(stderr, "No point light in scene.\n");
        return 0;
    }
    
    int photons_per_light = photon_map->max_photons / point_light_count;
    vec3* sphere_samples = (vec3*)malloc(sizeof(vec3) * photons_per_light);
    int sample_count = 0;
    while(sample_count < photons_per_light)
    {
        float x, y, z;
        do{
        x = (float)rand() / (float)RAND_MAX * 2.0f - 1.0f;
        y = (float)rand() / (float)RAND_MAX * 2.0f - 1.0f;
        z = (float)rand() / (float)RAND_MAX * 2.0f - 1.0f;        
        }while((x*x + y*y + z*z) > 1);
        vec3_assign(sphere_samples[sample_count], x, y, z);
        vec3_normalize(sphere_samples[sample_count], sphere_samples[sample_count]);
        sample_count++;
    }

    int photon_count = 0;
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == POINTLIGHT)
        {
            PointLight* point_light = (PointLight*)(sl->light_ptrs[i]);
            vec3 light_power, photon_power;
            vec3_scale(light_power, point_light->color, point_light->intensity);
            vec3_scale(photon_power, light_power, 1.0f/(float)photons_per_light);
            for(int j = 0; j < photons_per_light; j++)
            {
                Ray ray;
                vec3_copy(ray.origin, point_light->point);
                vec3_copy(ray.direction, sphere_samples[j]);

                // NOTE: limiting photon emission to one bounce for now
                ShadeRec sr;
                float t = intersectTest(&sr, so, ray);
                if(t < TMAX)
                {
                    Photon* cur_photon = &(photon_map->photons[photon_count]);
                    vec3_copy(cur_photon->position, sr.hit_point);
                    int theta = int(acos(sr.wo[2]) * (256.0/PI));
                    if(theta > 255)
                    {
                        cur_photon->theta = 255;
                    }else
                    {
                        cur_photon->theta = (unsigned char)theta;
                    }
                    int phi = int(atan2(sr.wo[1], sr.wo[0]) * (256.0 / (2.0f * PI)));
                    if(phi > 255)
                    {
                        cur_photon->phi = 255;
                    }else if(phi < 0)
                    {
                        cur_photon->phi = (unsigned char)(phi + 256);
                    }else
                    {
                        cur_photon->phi = (unsigned char)phi;
                    }
                    vec3_copy(cur_photon->power, photon_power);
                    for(int k = 0; k < 3; k++)
                    {
                        if(cur_photon->position[k] < photon_map->bbox.min[k])
                        {
                            photon_map->bbox.min[k] = cur_photon->position[k];
                        }
                        if(cur_photon->position[k] > photon_map->bbox.max[k])
                        {
                            photon_map->bbox.max[k] = cur_photon->position[k];
                        }
                    }
                    photon_count++;
                }
            }
        }
    }
    photon_map->stored_photons = photon_count;
    return photon_count;
}
