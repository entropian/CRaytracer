#include "disk.h"

float rayIntersectDisk(ShadeRec* sr, Disk* disk, const Ray ray)
{
    vec3 displacement;
    vec3_sub(displacement, disk->center, ray.origin);
    float t = vec3_dot(displacement, disk->normal) / vec3_dot(ray.direction, disk->normal);
    if(t > K_EPSILON)
    {
        vec3 hit_point;
        getPointOnRay(hit_point, ray, t);
        vec3_sub(displacement, hit_point, disk->center);
        if(vec3_length(displacement) <= disk->radius)
        {
            vec3_negate(sr->wo, ray.direction);
            if(vec3_dot(sr->wo, disk->normal) < 0)
            {
                vec3_negate(sr->wo, disk->normal);
            }else
            {
                vec3_copy(sr->normal, disk->normal);
            }
            vec3_copy(sr->hit_point, hit_point);
            sr->mat = *(disk->mat);
            return t;
        }
    }
    return TMAX;
}

float shadowRayIntersectDisk(const Disk* disk, const Ray ray)
{
    vec3 displacement;
    vec3_sub(displacement, disk->center, ray.origin);
    float t = vec3_dot(displacement, disk->normal) / vec3_dot(ray.direction, disk->normal);
    if(t > K_EPSILON)
    {
        vec3 hit_point;
        getPointOnRay(hit_point, ray, t);        
        vec3_sub(displacement, hit_point, disk->center);
        if(vec3_length(displacement) <= disk->radius)
        {
            return t;
        }
    }
    return TMAX;
}
