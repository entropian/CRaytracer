#include "lights.h"

void getAreaLightNormal(vec3 r, const AreaLight* area_light_ptr, const vec3 hit_point)
{
    switch(area_light_ptr->obj_type)
    {
    case RECTANGLE:
    {
        vec3_copy(r, ((Rectangle*)(area_light_ptr->obj_ptr))->normal);
    } break;
    case SPHERE:
    {
        vec3 displacement;
        vec3_sub(displacement, hit_point, ((Sphere*)(area_light_ptr->obj_ptr))->center);
        vec3_normalize(r, displacement);
    } break;
    }
}

void assignDirLight(DirLight *dir_light, const float intensity, const vec3 color, const vec3 direction)
{
    dir_light->intensity = intensity;
    vec3_copy(dir_light->color, color);
    vec3_copy(dir_light->direction, direction);
}

void assignPointLight(PointLight *point_light, const float intensity, const vec3 color, const vec3 point)
{
    point_light->intensity = intensity;
    vec3_copy(point_light->color, color);
    vec3_copy(point_light->point, point);
}

void getLightDir(vec3 r, const LightType light_type, const void* light_ptr, const ShadeRec* sr, const int sample_index)
{
    switch(light_type)
    {
    case DIRECTIONAL:
    {
        vec3_copy(r, ((DirLight*)light_ptr)->direction);
    } break;
    case POINTLIGHT:
    { 
        vec3 displacement;
        vec3_sub(displacement, ((PointLight*)light_ptr)->point, sr->hit_point);
        vec3_normalize(r, displacement);        
    } break;
    case AREALIGHT:
    {
        // Side effect: calculates and stores the surface sample point 
        AreaLight* area_light_ptr = (AreaLight*)light_ptr;
        switch(area_light_ptr->obj_type)
        {
        case RECTANGLE:
        {
            vec2 sample;
            //getNextSample2D(sample, area_light_ptr->samples2D);
            getSample2D(sample, area_light_ptr->samples2D, sample_index);
            vec3 displacement;
            Rectangle* rect = (Rectangle*)(area_light_ptr->obj_ptr);
            vec3_scale(displacement, rect->width, sample[0]);
            vec3_add(area_light_ptr->sample_point, rect->point, displacement);
            vec3_scale(displacement, rect->height, sample[1]);
            vec3_add(area_light_ptr->sample_point, area_light_ptr->sample_point, displacement);            
            vec3_sub(displacement, area_light_ptr->sample_point, sr->hit_point);
            vec3_normalize(r, displacement);
        } break;
        case SPHERE:
        {
            vec3 h_sample;
            //getNextSample3D(h_sample, area_light_ptr->samples3D);
            getSample3D(h_sample, area_light_ptr->samples3D, sample_index);
            getVec3InLocalBasis(area_light_ptr->sample_point, h_sample, sr->normal);
            vec3_add(area_light_ptr->sample_point, area_light_ptr->sample_point,
                     ((Sphere*)(area_light_ptr->obj_ptr))->center);
            vec3 displacement;
            vec3_sub(displacement, area_light_ptr->sample_point, sr->hit_point);
            vec3_normalize(r, displacement);
        } break;
        }
    }  break;
    case ENVLIGHT:
    {
        // get sample
        // calculate orthonormal basis based on the hit point and hit normal
        // transform sample by orthonormal basis
        EnvLight* env_light_ptr = (EnvLight*)light_ptr;
        vec3 h_sample;
        //getNextSample3D(h_sample, env_light_ptr->samples3D);
        getSample3D(h_sample, env_light_ptr->samples3D, sample_index);
        getVec3InLocalBasis(r, h_sample, sr->normal);        
    } break;
    }
}

void getIncRadiance(vec3 r, const LightType light_type, const void* light_ptr)
{
    switch(light_type)
    {
    case DIRECTIONAL:
    {
        vec3_scale(r, ((DirLight*)light_ptr)->color, ((DirLight*)light_ptr)->intensity);
    } break;
    case POINTLIGHT:
    {
        vec3_scale(r, ((PointLight*)light_ptr)->color, ((PointLight*)light_ptr)->intensity);
    } break;
    case ENVLIGHT:
    {
        vec3_scale(r, ((EnvLight*)light_ptr)->color, ((EnvLight*)light_ptr)->intensity);
    } break;
    case AREALIGHT:
    {
        vec3_scale(r, ((AreaLight*)light_ptr)->color, ((AreaLight*)light_ptr)->intensity);
    } break;
    }
}

float calcLightDistance(const LightType light_type, const void* light_ptr, const vec3 hit_point)
{
    float t;
    switch(light_type)
    {
    case DIRECTIONAL:
    {
        t = TMAX;                                                                    
    } break;
    case POINTLIGHT:
    {
        vec3 light_to_hit_point;
        PointLight* point_light_ptr = (PointLight*)(light_ptr);
        vec3_sub(light_to_hit_point, point_light_ptr->point, hit_point);
        t = vec3_length(light_to_hit_point);
    }break;
    case AREALIGHT:
    {
        vec3 light_to_hit_point;
        AreaLight* area_light_ptr = (AreaLight*)(light_ptr);
        vec3_sub(light_to_hit_point, area_light_ptr->sample_point, hit_point);
        t = vec3_length(light_to_hit_point);
    } break;
    case ENVLIGHT:
    {
        t = TMAX;
    } break;
    }
    return t;
}
