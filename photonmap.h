#include <GLFW/glfw3.h>
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

#define SHOW_TIME 

float calcFresnelReflectance(const ShadeRec*);
float calcTransmitDir(vec3, const ShadeRec*);
extern float intersectTest(ShadeRec* sr, const SceneObjects* so, const Ray ray);
static const int MAX_NPHOTONS = 1502;

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
    int max_bounce;
    float costheta[256];
    float sintheta[256];
    float cosphi[256];
    float sinphi[256];
    AABB bbox;
}Photonmap;

typedef struct PhotonQueryVars_s
{
    bool caustic;
    int num_photons;
    int num_caustic_photons;
    int nphotons;
    float photon_radius;
    float caustic_radius;
}PhotonQueryVars;

void Photonmap_init(Photonmap* photon_map, const int max_photons, const int max_bounce)
{
    photon_map->stored_photons = 0;
    photon_map->prev_scale = 1;
    photon_map->max_photons = max_photons;
    photon_map->max_bounce = max_bounce;
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
        photon_map->costheta[i] = (float)cos(angle);
        photon_map->sintheta[i] = (float)sin(angle);
        photon_map->cosphi[i] = (float)cos(2.0*angle);
        photon_map->sinphi[i] = (float)sin(2.0*angle);
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

void calcSphereSample(vec3 sphere_sample)
{
    float x, y, z;
    do{
        x = (float)rand() / (float)RAND_MAX * 2.0f - 1.0f;
        y = (float)rand() / (float)RAND_MAX * 2.0f - 1.0f;
        z = (float)rand() / (float)RAND_MAX * 2.0f - 1.0f;
    }while((x*x + y*y + z*z) > 1);
    vec3_assign(sphere_sample, x, y, z);
    vec3_normalize(sphere_sample, sphere_sample);
}

bool calcNewRayAndPhotonPower(Ray *ray, vec3 photon_power, const float t, const ShadeRec *sr, const Photonmap *photon_map,
                              const int sample_index)
{

    float rand_float = (float)rand() / (float)RAND_MAX;
    vec3 sample, ref_ray_dir;
    switch(sr->mat.mat_type)
    {
    case MATTE:
    {
        float reflectance_avg = (sr->mat.cd[0] + sr->mat.cd[1] + sr->mat.cd[2]) / 3.0f;
        if(rand_float > reflectance_avg){return false;}
        getSample3D(sample, sr->mat.h_samples, sample_index);
        getVec3InLocalBasis(ref_ray_dir, sample, sr->normal);
        photon_power[0] *= sr->mat.cd[0] / reflectance_avg;
        photon_power[1] *= sr->mat.cd[1] / reflectance_avg;
        photon_power[2] *= sr->mat.cd[2] / reflectance_avg;
    } break;
    case REFLECTIVE:
    {
        float reflectance_avg = (sr->mat.cs[0] + sr->mat.cs[1] + sr->mat.cs[2]) / 3.0f;
        if(rand_float > reflectance_avg){return false;}
        getSample3D(sample, sr->mat.h_samples, sample_index);
        vec3 reflect_dir;
        calcReflectRayDir(reflect_dir, sr->normal, ray->direction);
        getVec3InLocalBasis(ref_ray_dir, sample, reflect_dir);
        photon_power[0] *= sr->mat.cs[0] / reflectance_avg;
        photon_power[1] *= sr->mat.cs[1] / reflectance_avg;
        photon_power[2] *= sr->mat.cs[2] / reflectance_avg;
    } break;
    case TRANSPARENT:
    {
        float kr = calcFresnelReflectance(sr);
        vec3 transmit_dir;
        calcTransmitDir(transmit_dir, sr);
        float ndotwo = vec3_dot(sr->normal, sr->wo);
        vec3 color_filter_ref, color_filter_trans;
        // NOTE: the pairing between color_filter_ref and color_filter_trans,
        // and cf_in and cf_out are the reverse of that in the whittedTrace and pathTrace function,
        // because the tracing direction is reversed.
        if(ndotwo > 0.0f)
        {
            vec3_pow(color_filter_ref, sr->mat.cf_in, t);
            vec3_pow(color_filter_trans, sr->mat.cf_out, t);
        }else
        {
            vec3_pow(color_filter_ref, sr->mat.cf_out, t);
            vec3_pow(color_filter_trans, sr->mat.cf_in, t);
        }
        if(rand_float <= kr)
        {
            calcReflectRayDir(ref_ray_dir, sr->normal, ray->direction);
            vec3_mult(photon_power, photon_power, color_filter_ref);
        }else
        {
            vec3_copy(ref_ray_dir, transmit_dir);
            vec3_mult(photon_power, photon_power, color_filter_trans);
        }
    } break;
    case EMISSIVE:
        return false;
    }
    vec3_copy(ray->origin, sr->hit_point);
    vec3_copy(ray->direction, ref_ray_dir);
    return true;
}

void getPointLightPhoton(Ray *ray, vec3 photon_power, const PointLight *point_light)
{
    vec3_scale(photon_power, point_light->color, point_light->flux);
    vec3_copy(ray->origin, point_light->point);
    calcSphereSample(ray->direction);
}

void getPointLightPhotonCaustic(Ray *ray, vec3 photon_power, const PointLight *point_light)
{
    vec3_scale(photon_power, point_light->color, point_light->intensity);
    vec3_copy(ray->origin, point_light->point);
    if(point_light->proj_map)
    {
        float rand_float = (float)rand() / (float)RAND_MAX;
        float threshold = rand_float * point_light->proj_coverage;
        float theta_per_cell = (float)(PI / (double)THETA_ROW);
        float phi_per_cell = (float)(2.0 * PI / (double)PHI_COLUMN);
        float base_area = theta_per_cell * phi_per_cell;
        float sum = 0.0f;
        int i, j;
        for(i = 0; i < THETA_ROW; i++)
        {
            float angle_from_equator = (float)abs(i - THETA_ROW * 0.5f) * theta_per_cell;
            float cell_area = base_area * cos(angle_from_equator);
            bool reached = false;
            for(j = 0; j < PHI_COLUMN; j++)
            {
                if(point_light->proj_map[i*PHI_COLUMN + j])
                {
                    sum += cell_area;
                    if(sum >= threshold)
                    {
                        reached = true;
                        break;
                    }
                }
            }
            if(reached)
            {
                break;
            }
        }
        float rand_float2 = (float)rand() / (float)RAND_MAX;
        float theta = ((float)i + rand_float) * theta_per_cell;
        float phi = ((float)j + rand_float2) * phi_per_cell;
        float x, y, z;
        y = (float)cos(theta);
        x = (float)cos(phi) * (float)sin(theta);
        z = (float)sin(phi) * (float)sin(theta);
        vec3_assign(ray->direction, x, y, z);
        vec3_normalize(ray->direction, ray->direction);
    }
}

void getRectLightPhoton(Ray *ray, vec3 photon_power, const AreaLight *area_light,
                        const Samples3D *h_samples, const unsigned int sample_index)
{
    vec3 h_sample, light_normal, point_on_light;

    getSample3D(h_sample, h_samples, sample_index);
    getAreaLightNormal(light_normal, area_light, ORIGIN);
    float rand_float = (float)rand() / (float)RAND_MAX;
    if(rand_float < 0.5)
    {
        vec3_negate(light_normal, light_normal);
    }
    getVec3InLocalBasis(ray->direction, h_sample, light_normal);

    vec2 unit_square_sample;
    getSample2D(unit_square_sample, area_light->samples2D, sample_index + (rand() % 100));
    vec3 displacement;
    Rectangle* rect = (Rectangle*)(area_light->obj_ptr);
    vec3_scale(displacement, rect->width, unit_square_sample[0]);
    vec3_add(point_on_light, rect->point, displacement);
    vec3_scale(displacement, rect->height, unit_square_sample[1]);
    vec3_add(point_on_light, point_on_light, displacement);
    vec3_copy(ray->origin, point_on_light);

    vec3_scale(photon_power, area_light->color, area_light->flux);
}

void getSphereLightPhoton(Ray *ray, vec3 photon_power, const AreaLight *area_light, const Samples3D *h_samples,
                          const unsigned int sample_index)
{
    Sphere *sphere = (Sphere*)(area_light->obj_ptr);
    vec3 normal, displacement;
    calcSphereSample(normal);
    vec3_scale(displacement, normal, sphere->radius);
    vec3_add(ray->origin, displacement, sphere->center);

    vec3 h_sample;
    getSample3D(h_sample, h_samples, sample_index);
    getVec3InLocalBasis(ray->direction, h_sample, normal);

    vec3_scale(photon_power, area_light->color, area_light->intensity);
}

void getAreaLightPhoton(Ray *ray, vec3 photon_power, AreaLight *area_light, Samples3D *h_samples, unsigned int *sample_index)
{
    if(area_light->obj_type == RECTANGLE)
    {
        getRectLightPhoton(ray, photon_power, area_light, h_samples, *sample_index);
        (*sample_index)++;
    }else if(area_light->obj_type == SPHERE)
    {
        getSphereLightPhoton(ray, photon_power, area_light, h_samples, *sample_index);
        (*sample_index)++;
    }else
    {
        fprintf(stderr, "Wrong light geometry type.\n");
        return;
    }
}

void getAreaLightPhotonCaustic(Ray *ray, vec3 photon_power, AreaLight *area_light, Samples3D *h_samples, unsigned int *sample_index)
{
    if(area_light->obj_type == RECTANGLE)
    {
        getRectLightPhoton(ray, photon_power, area_light, h_samples, *sample_index);
        (*sample_index)++;
    }else if(area_light->obj_type == SPHERE)
    {
        getSphereLightPhoton(ray, photon_power, area_light, h_samples, *sample_index);
        (*sample_index)++;
    }else
    {
        fprintf(stderr, "Wrong light geometry type.\n");
        return;
    }
}

void storePhoton(Photon* photon, Photonmap *photon_map, const vec3 photon_power, const ShadeRec *sr)
{
    vec3_copy(photon->pos, sr->hit_point);
    int theta = int(acos(sr->wo[2]) * (256.0/PI));
    if(theta > 255)
    {
        photon->theta = 255;
    }else
    {
        photon->theta = (unsigned char)theta;
    }
    int phi = (int)(atan2(sr->wo[1], sr->wo[0]) * (256.0 / (2.0f * PI)));
    if(phi > 255)
    {
        photon->phi = 255;
    }else if(phi < 0)
    {
        photon->phi = (unsigned char)(phi + 256);
    }else
    {
        photon->phi = (unsigned char)phi;
    }
    vec3_copy(photon->power, photon_power);
    for(int k = 0; k < 3; k++)
    {
        if(photon->pos[k] < photon_map->bbox.min[k])
        {
            photon_map->bbox.min[k] = photon->pos[k];
        }
        if(photon->pos[k] > photon_map->bbox.max[k])
        {
            photon_map->bbox.max[k] = photon->pos[k];
        }
    }
}

int countLight(const SceneLights *sl)
{
    int light_count = 0;
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == POINTLIGHT)
        {
            light_count++;
        }else if(sl->light_types[i] == AREALIGHT)
        {
            AreaLight *area_light = (AreaLight*)(sl->light_ptrs[i]);
            light_count++;
        }
    }
    return light_count;
}

void emitPhotons(Photonmap* photon_map, const SceneObjects *so, const SceneLights *sl)
{
#ifdef SHOW_TIME
    double start_time, end_time;
    start_time = glfwGetTime();
#endif
    Photon *photons = photon_map->photons;
    int light_count = countLight(sl);
    int photons_per_light = photon_map->max_photons / light_count;

    Samples3D *h_samples = genHemisphereSamples(MULTIJITTERED, 1.0f);

    const int max_bounce = photon_map->max_bounce;
    int stored_photons = 0;
    unsigned int sample_index = 0;
    for(int i = 0; i < sl->num_lights; i++)
    {
        void* light_ptr = sl->light_ptrs[i];
        LightType light_type = sl->light_types[i];
        if(light_type != POINTLIGHT && light_type != AREALIGHT)
        {
            light_ptr = NULL;
        }
        if(!light_ptr){continue;}

        int bounce_count = 0;
        int emitted = 0;
        bool reflected = false;
        Ray ray;
        vec3 photon_power = {0.0f, 0.0f, 0.0f};
        int last_stored = stored_photons + 1;
        while(stored_photons - last_stored < photons_per_light)
        {
            if(!reflected || bounce_count == max_bounce)
            {
                switch(light_type)
                {
                case POINTLIGHT:
                    getPointLightPhoton(&ray, photon_power, (PointLight*)light_ptr);
                    break;
                case AREALIGHT:
                    getAreaLightPhoton(&ray, photon_power, (AreaLight*)light_ptr, h_samples, &sample_index);
                    break;
                }
                bounce_count = 0;
                emitted++;
            }
            ShadeRec sr;
            float t = intersectTest(&sr, so, ray);
            if(t < TMAX)
            {
                if(sr.mat.tex_flags != NO_TEXTURE)
                {
                    updateShadeRecWithTexInfo(&sr);
                }
                /*
                if(sr.mat.mat_type == MATTE
                   && bounce_count != 0) // Excluding first bounce so the photon map is only for indirection illum
                */
                if(sr.mat.mat_type == MATTE) // Using global photon map for final gather
                {
                    // Store photon if surface is matte
                    stored_photons++;
                    Photon *cur_photon = &(photons[stored_photons]);
                    storePhoton(cur_photon, photon_map, photon_power, &sr);
                }
                // Reflect photon
                if(calcNewRayAndPhotonPower(&ray, photon_power, t, &sr, photon_map, sample_index))
                {
                    reflected = true;
                    bounce_count++;
                    sample_index++;
                }else
                {
                    reflected = false;
                }
            }else
            {
                reflected = false;
            }
        }
        for(int j = last_stored; j <= stored_photons; j++)
        {
            vec3_scale(photons[j].power, photons[j].power, 1.0f / (float)emitted);
        }
    }
    photon_map->stored_photons = stored_photons;
    freeSamples3D(h_samples);
#ifdef SHOW_TIME
    end_time = glfwGetTime();
    double seconds = end_time - start_time;
    printf("Photon emission %f seconds.\n", seconds);
#endif
}

void emitCaustics(Photonmap* photon_map, const SceneObjects *so, const SceneLights *sl,
                  const AABB *aabb, const unsigned int num_aabb)
{
#ifdef SHOW_TIME
    double start_time, end_time;
    start_time = glfwGetTime();
#endif
    Photon *photons = photon_map->photons;
    int light_count = countLight(sl);
    int photons_per_light = photon_map->max_photons / light_count;

    Samples3D *h_samples = genHemisphereSamples(MULTIJITTERED, 1.0f);

    const int max_bounce = photon_map->max_bounce;
    int stored_photons = 0;
    unsigned int sample_index = 0;
    for(int i = 0; i < sl->num_lights; i++)
    {
        void* light_ptr = sl->light_ptrs[i];
        LightType light_type = sl->light_types[i];
        if(light_type != POINTLIGHT && light_type != AREALIGHT)
        {
            light_ptr = NULL;
        }
        if(!light_ptr){continue;}

        int bounce_count = 0;
        int emitted = 0;
        bool reflected = false;
        Ray ray;
        vec3 photon_power = {0.0f, 0.0f, 0.0f};
        int last_stored = stored_photons + 1;
        while(stored_photons - last_stored < photons_per_light)
        {
            if(!reflected || bounce_count == max_bounce)
            {
                /*
                switch(light_type)
                {
                case POINTLIGHT:
                    getPointLightPhotonCaustic(&ray, photon_power, (PointLight*)light_ptr);
                    break;
                case AREALIGHT:
                    // TODO
                    getAreaLightPhotonCaustic(&ray, photon_power, (AreaLight*)light_ptr,
                                              h_samples, &sample_index);
                    break;
                }
                emitted++;
                bounce_count = 0;
                */
                // TODO
                bool hits_caustic_bb = false;
                while(!hits_caustic_bb)
                {
                    switch(light_type)
                    {
                    case POINTLIGHT:
                        getPointLightPhotonCaustic(&ray, photon_power, (PointLight*)light_ptr);
                        break;
                    case AREALIGHT:
                        // TODO
                        getAreaLightPhotonCaustic(&ray, photon_power, (AreaLight*)light_ptr,
                                                  h_samples, &sample_index);
                        break;
                    }
                    emitted++;
                    for(int j = 0; j < num_aabb; j++)
                    {
                        if(rayIntersectAABB(&(aabb[j]), ray) < TMAX)
                        {
                            hits_caustic_bb = true;
                            break;
                        }
                    }
                }
                bounce_count = 0;

            }
            ShadeRec sr;
            float t = intersectTest(&sr, so, ray);
            if(t < TMAX)
            {
                if(sr.mat.tex_flags != NO_TEXTURE)
                {
                    updateShadeRecWithTexInfo(&sr);
                }
                // Problem: the order is wrong for reflection and power adjustment
                if(sr.mat.mat_type == MATTE && reflected)
                {
                    // Store photon if surface is diffuse
                    stored_photons++;
                    Photon *cur_photon = &(photons[stored_photons]);
                    storePhoton(cur_photon, photon_map, photon_power, &sr);
                    reflected = false;
                }
                // Reflect photon
                if((sr.mat.mat_type == REFLECTIVE || sr.mat.mat_type == TRANSPARENT) &&
                   calcNewRayAndPhotonPower(&ray, photon_power, t, &sr, photon_map, sample_index))
                {
                    reflected = true;
                    bounce_count++;
                    sample_index++;
                }else
                {
                    reflected = false;
                }
            }else
            {
                reflected = false;
            }
        }
        for(int j = last_stored; j <= stored_photons; j++)
        {
            vec3_scale(photons[j].power, photons[j].power, 1.0f / (float)emitted);
        }
    }
    photon_map->stored_photons = stored_photons;
    freeSamples3D(h_samples);
#ifdef SHOW_TIME
    end_time = glfwGetTime();
    double seconds = end_time - start_time;
    printf("Caustic emission %f seconds.\n", seconds);
#endif
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

void irradEstimate(vec3 irrad, const Photonmap *photon_map, const vec3 pos, const vec3 normal,
                   const float max_dist, const int nphotons)
{
    if(nphotons > MAX_NPHOTONS)
    {
        fprintf(stderr, "nphotons exceeds MAX_NPHOTONS.\n");
        return;
    }
    irrad[0] = irrad[1] = irrad[2] = 0.0f;

    float dist2[MAX_NPHOTONS];
    Photon* index[MAX_NPHOTONS];
    NearestPhotons np;
    vec3_copy(np.pos, pos);
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

    const float tmp = (float)(1.0 / PI) / (np.dist2[0]);
    vec3_scale(irrad, irrad, tmp);
}

void calcPhotonmapComponent(vec3 color, const vec3 h_sample, const PhotonQueryVars *query_vars,
                                 const Photonmap *photon_map, const Photonmap *caustic_map,
                                 const SceneObjects *so, const ShadeRec *sr, const vec3 caustic_rad)
{
    vec3 pm_color = {0.0f, 0.0f, 0.0f};
    if(sr->mat.mat_type == DIFFUSE)
    {
        Ray new_ray;
        vec3_copy(new_ray.origin, sr->hit_point);
        getVec3InLocalBasis(new_ray.direction, h_sample, sr->normal);
        ShadeRec new_sr;
        float t = intersectTest(&new_sr, so, new_ray);
        vec3 f;

        if(t < TMAX)
        {
            if(new_sr.mat.tex_flags != NO_TEXTURE)
            {
                updateShadeRecWithTexInfo(&new_sr);
            }
            vec3 irrad = {0.0f, 0.0f, 0.0f};
            irradEstimate(irrad, photon_map, new_sr.hit_point, new_sr.normal,
                          query_vars->photon_radius, query_vars->nphotons);

            vec3_scale(f, new_sr.mat.cd, new_sr.mat.kd / (float)PI);
            vec3 inc_rad;
            vec3_mult(inc_rad, irrad, f);

            float ndotwi = vec3_dot(sr->normal, new_ray.direction);
            float pdf = ndotwi / (float)PI;
            // Incoming radiance from secondary point * cosine theta
            vec3_scale(f, sr->mat.cd, sr->mat.kd / (float)PI * ndotwi / pdf);
            vec3_mult(pm_color, inc_rad, f);
        }

        if(query_vars->caustic)
        {
            vec3_add(pm_color, pm_color, caustic_rad);
        }
    }
    vec3_copy(color, pm_color);
}

void calcCausticBuffer(float* caustic_buffer, Camera *camera, Film *film, 
                       const Scene *scene, const Photonmap *caustic_map, const PhotonQueryVars *query_vars,
                       const int num_caustic_samples)
{
    Sampler sampler;
    Sampler_create(&sampler);
    for(int p = 0; p < num_caustic_samples; p++)
    {
        sampler.cur_sample_index = p;
        for(int i = 0; i < film->num_pixels; i++)
        {
            Sampler_setPixel(&sampler, i);
            vec2 imageplane_coord;
            vec2 sample;
            Sampler_getSample(sample, &sampler);
            calcImageCoord(imageplane_coord, film, sample, i);

            Ray ray;
            Sampler_getSample(sample, &sampler);
            calcCameraRay(&ray, imageplane_coord, camera, sample);

            ShadeRec sr;
            float t = intersectTest(&sr, &(scene->objects), ray);
            vec3 pm_color = {0.0f, 0.0f, 0.0f};
            if(t < TMAX)
            {
                vec3 caustic_irrad, caustic_rad;
                irradEstimate(caustic_irrad, caustic_map, sr.hit_point, sr.normal,
                              query_vars->caustic_radius, query_vars->nphotons);
                vec3 f;
                vec3_scale(f, sr.mat.cd, sr.mat.kd / PI);
                vec3_mult(caustic_rad, caustic_irrad, f);
                vec3_add(pm_color, pm_color, caustic_rad);
            }
            caustic_buffer[i*3] += pm_color[0];
            caustic_buffer[i*3+1] += pm_color[1];
            caustic_buffer[i*3+2] += pm_color[2];
        }
    }
    for(int i = 0; i < film->num_pixels; i++)
    {
        caustic_buffer[i*3] /= num_caustic_samples;
        caustic_buffer[i*3+1] /= num_caustic_samples;
        caustic_buffer[i*3+2] /= num_caustic_samples;
    }
}
