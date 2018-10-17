#include "sphere.h"

void fillShadeRecSphere(ShadeRec *sr, Sphere *sphere, const vec3 hit_point, const Ray ray, const float t,
                        const float theta, const float phi)
{
    sr->hit_status = true;
    vec3_copy(sr->hit_point, hit_point);
    vec3 hit_point_to_center;
    vec3_sub(hit_point_to_center, sr->hit_point, sphere->center);
    vec3_normalize(sr->normal, hit_point_to_center);
    vec3_negate(sr->wo, ray.direction);
    sr->mat = *(sphere->mat);
    
    const double theta_max = PI;
    const double theta_min = 0.0;
    const double phi_max = 2.0 * PI;    
    float u = phi / phi_max;
    float v = (theta - theta_min) / (theta_max - theta_min);
    vec2_assign(sr->uv, u, v);

    float y_radius = sqrt(hit_point[0] * hit_point[0] + hit_point[2] * hit_point[2]);
    float inv_y_radius = 1.0f / y_radius;
    float cos_phi = hit_point[0] * inv_y_radius;
    float sin_phi = hit_point[2] * inv_y_radius;
    vec3_assign(sr->dpdu, -phi_max * (hit_point[2] - sphere->center[2]), 0.0f,
                phi_max * (hit_point[0] - sphere->center[0]));
    vec3 tmp = {hit_point[1] * cos_phi, -sphere->radius * sinf(theta), hit_point[1] * sin_phi};
    vec3_scale(sr->dpdv, tmp, theta_max * theta_max);
    vec3_normalize(sr->dpdu, sr->dpdu);
    vec3_normalize(sr->dpdv, sr->dpdv);
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
        float e = (float)sqrt(disc);
        float denom = 2*a;
        float t = (-b - e)/denom;
        if(t > K_EPSILON)
        {
            vec3 hit_point;            
            getPointOnRay(hit_point, ray, t);                        
            //float phi = (float)atan2(hit_point[0] - sphere->center[0], hit_point[2] - sphere->center[2]);
            float phi = (float)atan2(hit_point[2] - sphere->center[2], hit_point[0] - sphere->center[0]);
            //if(phi < 0.0f) phi += 2.0*PI;
            float theta = (float)acos(((hit_point[1] - sphere->center[1]) / sphere->radius));
            if(fabs(phi) <= sphere->phi && theta >= sphere->min_theta&&
               theta <= sphere->max_theta)
            {
                fillShadeRecSphere(sr, sphere, hit_point, ray, t, theta, phi);
                return t;
            }
        }

        t = (-b + e)/denom;
        if(t > K_EPSILON)
        {
            vec3 hit_point;
            getPointOnRay(hit_point, ray, t);
            float phi = (float)atan2(hit_point[0] - sphere->center[0], hit_point[2] - sphere->center[2]);
            //if(phi < 0.0f) phi += 2.0*PI;
            float theta = (float)acos(((hit_point[1] - sphere->center[1]) / sphere->radius));
            if(fabs(phi) <= sphere->phi && theta >= sphere->min_theta&&
               theta <= sphere->max_theta)
            {
                fillShadeRecSphere(sr, sphere, hit_point, ray, t, theta, phi);
                return t;
            }
        }
    }
    return TMAX;
}

float shadowRayIntersectSphere(Sphere *sphere, const Ray ray)
{
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
        float e = (float)sqrt(disc);
        float denom = 2*a;
        float t = (-b - e)/denom;
        if(t > K_EPSILON)
        {
            vec3 hit_point;
            getPointOnRay(hit_point, ray, t);
            float phi = (float)atan2(hit_point[0] - sphere->center[0], hit_point[2] - sphere->center[2]);
            float theta = (float)acos(((hit_point[1] - sphere->center[1]) / sphere->radius));
            if(fabs(phi) <= sphere->phi && theta >= sphere->min_theta&&
               theta <= sphere->max_theta)
            {
                return t;
            }
        }

        t = (-b + e)/denom;
        if(t > K_EPSILON)
        {
            vec3 hit_point;            
            getPointOnRay(hit_point, ray, t);
            float phi = (float)atan2(hit_point[0] - sphere->center[0], hit_point[2] - sphere->center[2]);
            float theta = (float)acos(((hit_point[1] - sphere->center[1]) / sphere->radius));
            if(fabs(phi) <= sphere->phi && theta >= sphere->min_theta&&
               theta <= sphere->max_theta)                        
            {
                return t;
            }            
        }
    }
    return TMAX;
}
