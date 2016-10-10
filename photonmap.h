#include <stdlib.h>
#include <math.h>
#include "util/vec.h"
#include "scene/scenedata.h"
#include "util/ray.h"
#include "objloader/dbuffer.h"
#include "aabb.h"
#include "shapes/shapes.h"
#include "sampling.h"
#include "util/math.h"

extern float intersectTest(ShadeRec* sr, const SceneObjects* so, const Ray ray);
static const int MAX_NPHOTONS = 202;
vec3 max_power = {0.0f, 0.0f, 0.0f};
float max_component = 0.0f;
vec3 max_pos = {0.0f, 0.0f, 0.0f};
int max_bounces = 0;

typedef struct Photon_s
{
    vec3 pos;
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

void Photonmap_init(Photonmap* photon_map, const int max_photons)
{
    photon_map->stored_photons = 0;
    photon_map->prev_scale = 1;
    photon_map->max_photons = max_photons;
    photon_map->photons = (Photon*)malloc(sizeof(Photon) * (max_photons + 1));
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

bool calcNewRayAndPhotonPower(Ray *ray, vec3 photon_power, const ShadeRec *sr, const Photonmap *photon_map,
    const int sample_index)
{

    float rand_float = (float)rand() / (float)RAND_MAX; // Random float for Russian roulette
    vec3 sample, ref_ray_dir;
    // If random float is less than or equal to the corresponding reflectance coefficient,
    // then cast a new ray
    switch(sr->mat->mat_type)
    {                        
    case MATTE:
    {
        float reflectance_avg = (sr->mat->cd[0] + sr->mat->cd[1] + sr->mat->cd[2]) / 3.0f;                        
        if(rand_float > reflectance_avg){return false;}
        //getSample3D(sample, sr->mat->h_samples, sample_index++);
        getSample3D(sample, sr->mat->h_samples, sample_index);
        getVec3InLocalBasis(ref_ray_dir, sample, sr->normal);
        photon_power[0] *= sr->mat->cd[0] / reflectance_avg;
        photon_power[1] *= sr->mat->cd[1] / reflectance_avg;
        photon_power[2] *= sr->mat->cd[2] / reflectance_avg;                        
    } break;
    case REFLECTIVE:                    
    {
        float reflectance_avg = (sr->mat->cr[0] + sr->mat->cr[1] + sr->mat->cr[2]) / 3.0f;                 
        if(rand_float > reflectance_avg){return false;}                        
        vec3 reflect_dir;
        calcReflectRayDir(reflect_dir, sr->normal, ray->direction);
        getVec3InLocalBasis(ref_ray_dir, sample, reflect_dir);
        photon_power[0] *= sr->mat->cr[0] / reflectance_avg;
        photon_power[1] *= sr->mat->cr[1] / reflectance_avg;
        photon_power[2] *= sr->mat->cr[2] / reflectance_avg;                                                
    } break;
    case TRANSPARENT:
        // TODO
        return false;
    case EMISSIVE:
        return false;
    }
    vec3_copy(ray->origin, sr->hit_point);
    vec3_copy(ray->direction, ref_ray_dir);
    return true;
}

void getPointLightPhoton(Ray *ray, vec3 photon_power, const PointLight *point_light,
                         const vec3* sphere_samples, const int sample_index)
{
    vec3_scale(photon_power, point_light->color, point_light->intensity);    
    vec3_copy(ray->origin, point_light->point);
    vec3_copy(ray->direction, sphere_samples[sample_index]);
}

/*
void getRectLightPhoton(Ray *ray, vec3 photon_power, const AreaLight *area_light,
                        const vec3* sphere_samples, const int sphere_sample_index,
                        const h_samples, const int h_sample_index)
{
    vec3 hemipshere_sample, light_normal, point_on_light;
    vec2 unit_square_sample;
    
    getSample3D(hemisphere_sample, h_samples, h_sample_index);    
    getAreaLightNormal(light_normal, area_light, ORIGIN);
    getVec3InLocalBasis(ray.direction, hemisphere_sample, light_normal);
}
*/

void emitPhotons(Photonmap* photon_map, const SceneObjects *so, const SceneLights *sl)
{    
    int point_light_count = 0;
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == POINTLIGHT)
        {
            point_light_count++;
        }else if(sl->light_types[i] == AREALIGHT)
        {
            // TODO
        }        
    }
    if(point_light_count == 0)
    {
        fprintf(stderr, "No point light in scene.\n");
        return;
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

    const int max_bounce = 5;
    int stored_photons = 0;
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == POINTLIGHT)
        {
            int sample_index = 0, photon_count = 0, emitted = 0;
            int last_stored = stored_photons;
            PointLight* point_light = (PointLight*)(sl->light_ptrs[i]);
            vec3 photon_power = {0.0f, 0.0f, 0.0f};
            Ray ray;            
            bool reflected = false;
            int bounce_count = 0, sphere_sample_index = 0;
            // TODO: fix below. Emitted photons and stored photons don't match up due to bounces.
            // NOTE: How do bounced photons relate to photon power overall?
            while(photon_count < photons_per_light && stored_photons <= photon_map->max_photons)
            {
                // Contruct new primary ray if last ray wasn't reflected or max bounce exeeded
                if(!reflected || bounce_count > max_bounce)
                {
                    getPointLightPhoton(&ray, photon_power, point_light,
                                        sphere_samples, sphere_sample_index);
                    vec3_scale(photon_power, photon_power, 1.0f / photons_per_light);
                    sphere_sample_index = (sphere_sample_index + 1) % photons_per_light;
                    bounce_count = 0;
                    emitted++;
                }

                ShadeRec sr;
                float t = intersectTest(&sr, so, ray);
                if(t < TMAX)
                {

                    // Store photon if hit surface is a diffuse surface
                    if(sr.mat->mat_type == MATTE)
                    {
                        stored_photons++;                        
                        Photon* cur_photon = &(photon_map->photons[stored_photons]);
                        vec3_copy(cur_photon->pos, sr.hit_point);
                        int theta = int(acos(sr.wo[2]) * (256.0/PI));
                        if(theta > 255)
                        {
                            cur_photon->theta = 255;
                        }else
                        {
                            cur_photon->theta = (unsigned char)theta;
                        }
                        int phi = (int)(atan2(sr.wo[1], sr.wo[0]) * (256.0 / (2.0f * PI)));
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
                            if(cur_photon->pos[k] < photon_map->bbox.min[k])
                            {
                                photon_map->bbox.min[k] = cur_photon->pos[k];
                            }
                            if(cur_photon->pos[k] > photon_map->bbox.max[k])
                            {
                                photon_map->bbox.max[k] = cur_photon->pos[k];
                            }
                        }
                        photon_count++;
                    }
                    if(calcNewRayAndPhotonPower(&ray, photon_power, &sr, photon_map, sample_index))
                    {
                        reflected = true;
                        bounce_count++;                                            
                        sample_index++;
                    }else
                    {
                        reflected = false;
                    }
                }else // Ray hits nothing
                {
                    reflected = false;
                }
            }            
            //for(int k = last_stored; k <= stored_photons; k++)
            for(int k = 0; k <= photon_map->max_photons + 2; k++)
            {
                //vec3_scale(photon_map->photons[k].power, photon_map->photons[k].power, 1.0f/emitted);
            }
        }
    }
    photon_map->stored_photons = stored_photons;
    free(sphere_samples);
}

#define swap(ph, a, b) {Photon *ph2 = ph[a]; ph[a] = ph[b]; ph[b] = ph2;}

void medianSplit(Photon **p, const int start, const int end, const int median, const int axis)
{
    int left = start;
    int right = end;
    while(right > left)
    {
        const float v = p[right]->pos[axis];
        int i = left - 1;
        int j = right;
        for(;;)
        {
            while(p[++i]->pos[axis] < v)
                ;
            while(p[--j]->pos[axis] > v && j > left)
                ;
            if(i >= j)
                break;
            swap(p, i, j);
        }

        swap(p, i, right);
        if(i >= median)
        {
            right = i - 1;
        }
        if(i <= median)
        {
            left = i + 1;
        }
    }
}

void balanceSegment(Photonmap *photon_map,
                    Photon **pbal, Photon **porg, const int index, const int start, const int end)
{
    int median = 1;
    while((4*median) <= (end - start + 1))
    {
        median += median;
    }
    if((3*median) <= (end - start + 1))
    {
        median += median;
        median += start - 1;
    }else
    {
        median = end - median + 1;
    }

    int axis = 2;
    float x_extent = photon_map->bbox.max[0] - photon_map->bbox.min[0];
    float y_extent = photon_map->bbox.max[1] - photon_map->bbox.min[1];
    float z_extent = photon_map->bbox.max[2] - photon_map->bbox.min[2];

    if(x_extent > y_extent && x_extent > z_extent)
    {
        axis = 0;
    }else if(y_extent > z_extent)
    {
        axis = 1;
    }

    medianSplit(porg, start, end, median, axis);

    pbal[index] = porg[median];
    pbal[index]->plane = axis;

    if(median > start)
    {
        if(start < median - 1)
        {
            const float tmp = photon_map->bbox.max[axis];
            photon_map->bbox.max[axis] = pbal[index]->pos[axis];
            balanceSegment(photon_map, pbal, porg, 2*index, start, median-1);
            photon_map->bbox.max[axis] = tmp;            
        }else
        {
            pbal[2*index] = porg[start];
        }
    }

    if(median < end)
    {
        if(median+1 < end)
        {
            const float tmp = photon_map->bbox.min[axis];
            photon_map->bbox.min[axis] = pbal[index]->pos[axis];
            balanceSegment(photon_map, pbal, porg, 2*index+1, median+1, end);
            photon_map->bbox.min[axis] = tmp;                        
        }else
        {
            pbal[2*index+1] = porg[end];
        }
    }
}

void Photonmap_balance(Photonmap* photon_map)
{
    if(photon_map->stored_photons > 1)
    {
        Photon* photons = photon_map->photons;
        Photon** pa1 = (Photon**)malloc(sizeof(Photon*) * (photon_map->stored_photons+1));
        Photon** pa2 = (Photon**)malloc(sizeof(Photon*) * (photon_map->stored_photons+1));

        for(int i = 0; i <= photon_map->stored_photons; i++)
        {
            pa2[i] = &(photons[i]);
        }
        balanceSegment(photon_map, pa1, pa2, 1, 1, photon_map->stored_photons);
        free(pa2);

        int d, j = 1, foo = 1;
        Photon foo_photon = photons[j];

        for(int i = 1; i <= photon_map->stored_photons; i++)
        {
            d = pa1[j] - photons;
            pa1[j] = NULL;
            if(d != foo)
            {
                photons[j] = photons[d];
            }else
            {
                photons[j] = foo_photon;
                if(i < photon_map->stored_photons)
                {
                    for(; foo <= photon_map->stored_photons; foo++)
                    {
                        if(pa1[foo])
                            break;
                    }
                    foo_photon = photons[foo];
                    j = foo;
                }
                continue;
            }
            j = d;
        }
        free(pa1);
    }
    photon_map->half_stored_photons = photon_map->stored_photons / 2 - 1;
}

typedef struct NearestPhotons_s
{
    int max;
    int found;
    int got_heap;
    float pos[3];
    float *dist2;
    Photon **index;
}NearestPhotons;

void locatePhotons(NearestPhotons *const np, const Photonmap* photon_map, const int index) 
{
    Photon* photons = photon_map->photons;
    Photon *p = &(photons[index]);
    float dist1;

    if(index < photon_map->half_stored_photons)
    {
        dist1 = np->pos[p->plane] - p->pos[p->plane];
        if(dist1 > 0.0f) // If dist1 is positive, search right plane
        {
            locatePhotons(np, photon_map, 2*index+1);
            if(dist1*dist1 < np->dist2[0])
            {
                locatePhotons(np, photon_map, 2*index);
            }
        }else  // If dist1 is negative, search left first
        {
            locatePhotons(np, photon_map, 2*index);
            if(dist1*dist1 < np->dist2[0])
            {
                locatePhotons(np, photon_map, 2*index+1);
            }
        }
    }

    // Compute squared distance between current photon and np->pos
    dist1 = p->pos[0] - np->pos[0];
    float dist2 = dist1*dist1;
    dist1 = p->pos[1] - np->pos[1];
    dist2 += dist1*dist1;
    dist1 = p->pos[2] - np->pos[2];
    dist2 += dist1*dist1;

    if(dist2 < np->dist2[0])
    {
        // We found a photon :) Insert it in the candidate list
        if(np->found < np->max)
        {
            // heap is not full; use array
            if(np->max == 1) 
            {
                np->dist2[0] = dist2;
            }
            (np->found)++;    
            np->dist2[np->found] = dist2;
            np->index[np->found] = p;            
        }else
        {
            int j, parent;
            if(np->got_heap == 0)
            {
                // Build heap
                float dst2;
                Photon *phot;
                int half_found = np->found >> 1;
                for(int k = half_found; k >= 1; k--)
                {
                    parent = k; 
                    phot = np->index[k];
                    dst2 = np->dist2[k];
                    while(parent <= half_found)
                    {
                        j = parent + parent;
                        if(j < np->found && np->dist2[j] < np->dist2[j+1])
                        {
                            j++;
                        }
                        if(dst2 >= np->dist2[j])
                        {
                            break;
                        }
                        np->dist2[parent] = np->dist2[j];
                        np->index[parent] = np->index[j];
                        parent = j;
                    }
                    np->dist2[parent] = dst2;
                    np->index[parent] = phot;
                }
                np->got_heap = 1;
            }
            // Insert new photon into max heap
            // Delete largest element, insert new, and reorder the heap
            parent = 1;
            j = 2;
            while(j <= np->found)
            {
                if(j < np->found && np->dist2[j] < np->dist2[j+1])
                {
                    j++;
                }
                if(dist2 > np->dist2[j])
                {
                    break;
                }
                np->dist2[parent] = np->dist2[j];
                np->index[parent] = np->index[j];
                parent = j;
                j += j;
            }
            np->index[parent] = p;
            np->dist2[parent] = dist2;
            np->dist2[0] = np->dist2[1];

        }
    }
}

void photonDir(vec3 dir, const Photonmap *photon_map, const Photon *p)
{
    dir[0] = photon_map->sintheta[p->theta] * photon_map->cosphi[p->phi];
    dir[1] = photon_map->sintheta[p->theta] * photon_map->sinphi[p->phi];
    dir[2] = photon_map->costheta[p->theta];
}


static float dist2[MAX_NPHOTONS];
static Photon* index[MAX_NPHOTONS];

void irradEstimate(vec3 irrad, const Photonmap *photon_map, const vec3 pos, const vec3 normal,
                   const float max_dist, const int nphotons)
{
    if(nphotons > MAX_NPHOTONS)
    {
        fprintf(stderr, "nphotons exceeds MAX_NPHOTONS.\n");
        return;
    }
    irrad[0] = irrad[1] = irrad[2] = 0.0f;

    NearestPhotons np;
    vec3_copy(np.pos, pos);
    //float dist2[MAX_NPHOTONS];
    //Photon* index[MAX_NPHOTONS];
    np.dist2 = dist2;
    np.index = index;
    np.max = nphotons;
    np.found = 0;
    np.got_heap = 0;
    np.dist2[0] = max_dist * max_dist;

    locatePhotons(&np, photon_map, 1);

    if(np.found < 8)
    {
        return;
    }

    vec3 pdir;

    // sum irradiance from all photons
    for(int i = 1; i <= np.found; i++)
    {
        const Photon *p = np.index[i];
        // the photon_dir call and following can be omitted for speed
        // if the scene does not have any thin surfaces
        photonDir(pdir, photon_map, p);
        if(vec3_dot(pdir, normal) > 0.0f)
        {
            vec3_add(irrad, irrad, p->power);
        }
    }

    const float tmp = (1.0f / PI) / (np.dist2[0]);
    vec3_scale(irrad, irrad, tmp);
}
