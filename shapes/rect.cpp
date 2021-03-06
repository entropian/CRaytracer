#include "rect.h"

float rayIntersectRect(ShadeRec *sr, Rectangle *rect, const Ray ray)
{
    vec3 displacement;
    vec3_sub(displacement, rect->point, ray.origin);
    vec3 normal;
    vec3_copy(normal, rect->normal);
    if(vec3_dot(displacement, rect->normal) < 0.0f)
    {
        vec3_negate(normal, rect->normal);
    }else
    {
        //vec3_negate(normal, rect->normal);
        vec3_copy(normal, rect->normal);
    }

    float t = vec3_dot(displacement, normal) / vec3_dot(ray.direction, normal);
    if(t > K_EPSILON)
    {
        vec3 point;
        getPointOnRay(point, ray, t);
        vec3_sub(displacement, point, rect->point);
        //float dot_product = vec3_dot(displacement, rect->width);
        float u = vec3_dot(displacement, rect->width) / vec3_dot(rect->width, rect->width);
        //if(dot_product >= 0 && dot_product <= vec3_dot(rect->width, rect->width))
        if(u >= 0.0f && u <= 1.0f)
        {
            //dot_product = vec3_dot(displacement, rect->height);
            float v = vec3_dot(displacement, rect->height) / vec3_dot(rect->height, rect->height);
            //if(dot_product >= 0 && dot_product <= vec3_dot(rect->height, rect->height))
            if(v >= 0.0f && v <= 1.0f)
            {
                vec2_assign(sr->uv, u, v);
                vec3_normalize(sr->dpdu, rect->width);
                vec3_normalize(sr->dpdv, rect->height);
                vec3_copy(sr->normal, normal);
                vec3_scale(displacement, ray.direction, t);
                vec3_copy(sr->hit_point, point);
                vec3_negate(sr->wo, ray.direction);
                if(vec3_dot(sr->wo, sr->normal) < 0.0f)
                {
                    vec3_negate(sr->normal, sr->normal);
                    vec3_negate(sr->dpdu, sr->dpdu);
                }
                sr->mat = *(rect->mat);
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
    vec3 normal;
    if(vec3_dot(displacement, rect->normal) < 0.0f)
    {
        vec3_negate(normal, rect->normal);
    }else
    {
        vec3_negate(normal, rect->normal);
    }
    float t = vec3_dot(displacement, normal) / vec3_dot(ray.direction, normal);
    if(t > K_EPSILON)
    {
        vec3 point;
        getPointOnRay(point, ray, t);        
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
