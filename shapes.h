#pragma once

#include "vec.h"
#include "materials.h"
#include "constants.h"

enum ObjectType
{
    SPHERE,
    PLANE,
    RECTANGLE
};

typedef struct Sphere
{
    bool shadow;
    float radius;    
    vec3 center;
    Material mat;
} Sphere;

typedef struct Plane
{
    bool shadow;
    vec3 point;
    vec3 normal;
    Material mat;
} Plane;

typedef struct Rectangle
{
    bool shadow;
    vec3 point;                // Bottom left of the rectangle
    vec3 width, height;
    vec3 normal;
    Material mat;
} Rectangle;

typedef struct Ray
{
    vec3 origin;
    vec3 direction;
} Ray;

typedef struct ShadeRec
{
    bool hit_status;
    Material *mat;                // Points to the material member of an object
    vec3 hit_point;
    vec3 normal;
} ShadeRec;

void fillShadeRecSphere(ShadeRec *sr, Sphere *sphere, const Ray ray, const float t)
{
    sr->hit_status = true;
    vec3 displacement;
    vec3_scale(displacement, ray.direction, t);
    vec3_add(sr->hit_point, ray.origin, displacement);
    vec3_sub(displacement, sr->hit_point, sphere->center);
    vec3_normalize(sr->normal, displacement);
    sr->mat = &(sphere->mat);
}

float rayIntersectSphere(ShadeRec *sr, Sphere *sphere, const Ray ray)
{
    // The analytic solution is so much better than the shitty loop from before
    // though I probably should have used smaller increment for the loop
    float a = vec3_dot(ray.direction, ray.direction);
    vec3 origin_to_center, tmp;
    vec3_sub(origin_to_center, ray.origin, sphere->center);
    vec3_scale(tmp, origin_to_center, 2);
    float b = vec3_dot(tmp, ray.direction);
    float c = vec3_dot(origin_to_center, origin_to_center) - sphere->radius*sphere->radius;
    float disc = b*b - 4*a*c;

    if(disc < 0.0f)
    {
        return TMAX;
    }else
    {
        float e = sqrt(disc);
        float denom = 2*a;
        float t = (-b - e)/denom;
        if(t > k_epsilon)
        {
            fillShadeRecSphere(sr, sphere, ray, t);            
            return t;
        }

        t = (-b + e)/denom;
        if(t > k_epsilon)
        {
            fillShadeRecSphere(sr, sphere, ray, t);            
            return t;
        }
    }
    return TMAX;
}

float shadowRayIntersectSphere(Sphere *sphere, const Ray ray)
{
    if(!sphere->shadow)
    {
        return TMAX;
    }
    float a = vec3_dot(ray.direction, ray.direction);
    vec3 origin_to_center, tmp;
    vec3_sub(origin_to_center, ray.origin, sphere->center);
    vec3_scale(tmp, origin_to_center, 2);
    float b = vec3_dot(tmp, ray.direction);
    float c = vec3_dot(origin_to_center, origin_to_center) - sphere->radius*sphere->radius;
    float disc = b*b - 4*a*c;

    if(disc < 0.0f)
    {
        return TMAX;
    }else
    {
        float e = sqrt(disc);
        float denom = 2*a;
        float t = (-b - e)/denom;
        if(t > k_epsilon)
        {
            return t;
        }

        t = (-b + e)/denom;
        if(t > k_epsilon)
        {
            return t;
        }
    }
    return TMAX;
}

float rayIntersectPlane(ShadeRec *sr, Plane *plane, const Ray ray)
{
    vec3 displacement;
    vec3_sub(displacement, plane->point, ray.origin);
    float t = vec3_dot(displacement, plane->normal) / vec3_dot(ray.direction, plane->normal);
    if(t > k_epsilon)
    {
        vec3_copy(sr->normal, plane->normal);
        vec3_scale(displacement, ray.direction, t);
        vec3_add(sr->hit_point, ray.origin, displacement);
        sr->mat = &(plane->mat);
        return t;
    }
    return TMAX;
}

float shadowRayIntersectPlane(Plane *plane, const Ray ray)
{
    if(!plane->shadow)
    {
        return TMAX;
    }
    vec3 displacement;
    vec3_sub(displacement, plane->point, ray.origin);
    float t = vec3_dot(displacement, plane->normal) / vec3_dot(ray.direction, plane->normal);
    if(t > k_epsilon)
    {
        return t;
    }
    return TMAX;
}

float rayIntersectRect(ShadeRec *sr, Rectangle *rect, const Ray ray)
{
    vec3 displacement;
    vec3_sub(displacement, rect->point, ray.origin);
    float t = vec3_dot(displacement, rect->normal) / vec3_dot(ray.direction, rect->normal);
    if(t > k_epsilon)
    {
        vec3 point;
        vec3_scale(displacement, ray.direction, t);
        vec3_add(point, ray.origin, displacement);
        vec3_sub(displacement, point, rect->point);
        float dot_product = vec3_dot(displacement, rect->width);
        if(dot_product >= 0 && dot_product <= vec3_dot(rect->width, rect->width))
        {
            dot_product = vec3_dot(displacement, rect->height);
            if(dot_product >= 0 && dot_product <= vec3_dot(rect->height, rect->height))
            {
                vec3_copy(sr->normal, rect->normal);
                vec3_scale(displacement, ray.direction, t);
                vec3_copy(sr->hit_point, point);
                sr->mat = &(rect->mat);                
                return t;
            }
        }
    }
    return TMAX;
}

float shadowRayIntersectRect(Rectangle *rect, const Ray ray)
{
    vec3 displacement;
    vec3_sub(displacement, rect->point, ray.origin);
    float t = vec3_dot(displacement, rect->normal) / vec3_dot(ray.direction, rect->normal);
    if(t > k_epsilon)
    {
        vec3 point;
        vec3_scale(displacement, ray.direction, t);
        vec3_add(point, ray.origin, displacement);
        vec3_sub(displacement, point, rect->point);
        float dot_product = vec3_dot(displacement, rect->width);
        if(dot_product >= 0 && dot_product <= vec3_dot(rect->width, rect->width))
        {
            dot_product = vec3_dot(displacement, rect->height);
            if(dot_product >= 0 && dot_product <= vec3_dot(rect->height, rect->height))
            {
                return t;
            }
        }
    }
    return TMAX;
}
