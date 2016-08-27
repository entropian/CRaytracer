#include "generic.h"

void fillShadeRecGenericOpenCylinder(ShadeRec *sr, GenericOpenCylinder* oc, const Ray ray, const vec3 hit_point)
{
    sr->hit_status = true;
    vec3_negate(sr->wo, ray.direction);
    vec3_copy(sr->hit_point, hit_point);
    switch(oc->normal_type)
    {
    case OPEN:
    {
        vec3_assign(sr->normal, hit_point[0]/oc->radius, 0.0f, hit_point[2]/oc->radius);
        if(vec3_dot(sr->wo, sr->normal) < 0)
        {
            vec3_negate(sr->normal, sr->normal);
        }
    } break;
    case CONVEX:
    {
        vec3_assign(sr->normal, hit_point[0]/oc->radius, 0.0f, hit_point[2]/oc->radius);        
    } break;
    case CONCAVE:
    {
        vec3_assign(sr->normal, (-hit_point[0])/oc->radius, 0.0f, (-hit_point[2])/oc->radius);
    } break;
    }
    sr->mat = oc->mat;
}

float rayIntersectGenericOpenCylinder(ShadeRec* sr, GenericOpenCylinder* oc, const Ray ray)
{
    // intersection equation: at^2 + bt + c = 0 quadratic
    float a = ray.direction[0]*ray.direction[0] + ray.direction[2]*ray.direction[2];
    float b = 2*(ray.origin[0]*ray.direction[0] + ray.origin[2]*ray.direction[2]);
    float c = ray.origin[0]*ray.origin[0] + ray.origin[2]*ray.origin[2] - oc->radius*oc->radius;
    float disc = b*b - 4*a*c;

    if(disc < 0.0f)
    {
        return TMAX;
    }else
    {
        float e = (float)sqrt(disc);
        float denom = 2*a;
        float t = (-b - e)/denom;
        vec3 hit_point;        
        if(t > K_EPSILON)
        {
            getPointOnRay(hit_point, ray, t);
            if(fabs(hit_point[1]) <= oc->half_height)
            {
                float phi = (float)atan2(hit_point[0], hit_point[2]);
                if(fabs(phi) <= oc->phi)
                {
                    fillShadeRecGenericOpenCylinder(sr, oc, ray, hit_point);
                    return t;                                        
                }
            }
        }
        t = (-b + e)/denom;
        if(t > K_EPSILON)
        {
            getPointOnRay(hit_point, ray, t);
            if(fabs(hit_point[1]) <= oc->half_height)
            {
                float phi = (float)atan2(hit_point[0], hit_point[2]);
                if(fabs(phi) <= oc->phi)
                {
                    fillShadeRecGenericOpenCylinder(sr, oc, ray, hit_point);
                    return t;                    
                }
            }
        }
    }
    return TMAX;
}

float shadowRayIntersectGenericOpenCylinder(const GenericOpenCylinder* oc, const Ray ray)
{
    float a = ray.direction[0]*ray.direction[0] + ray.direction[2]*ray.direction[2];
    float b = 2*(ray.origin[0]*ray.direction[0] + ray.origin[2]*ray.direction[2]);
    float c = ray.origin[0]*ray.origin[0] + ray.origin[2]*ray.origin[2] - oc->radius*oc->radius;
    float disc = b*b - 4*a*c;

    if(disc < 0.0f)
    {
        return TMAX;
    }else
    {
        float e = (float)sqrt(disc);
        float denom = 2*a;
        float t = (-b - e)/denom;
        vec3 hit_point;        
        if(t > K_EPSILON)
        {
            getPointOnRay(hit_point, ray, t);            
            if(fabs((float)hit_point[1]) <= oc->half_height)
            {
                float phi = (float)atan2(hit_point[0], hit_point[2]);
                if(fabs(phi) <= oc->phi)
                {
                    return t;
                }
            }
        }
        t = (-b + e)/denom;
        if(t > K_EPSILON)
        {
            getPointOnRay(hit_point, ray, t);            
            if(fabs(hit_point[1]) <= oc->half_height)
            {
                float phi = (float)atan2(hit_point[0], hit_point[2]);
                if(fabs(phi) <= oc->phi)
                {
                    return t;
                }
            }
        }
    }
    return TMAX;
}

void computeGenericTorusNormal(vec3 r, const GenericTorus* torus, const vec3 hit_point)
{
    float param_squared = torus->swept_radius*torus->swept_radius +
        torus->tube_radius*torus->tube_radius;

    float x = hit_point[0];
    float y = hit_point[1];
    float z = hit_point[2];
    float sum_squared = x*x + y*y + z*z;

    r[0] = 4.0f * x * (sum_squared - param_squared);
    r[1] = 4.0f * y * (sum_squared - param_squared + 2.0f * torus->swept_radius * torus->swept_radius);
    r[2] = 4.0f * z * (sum_squared - param_squared);
    vec3_normalize(r, r);
}

void calcAABBGenericTorus(AABB* aabb, GenericTorus* torus)
{
    //AABB* aabb = &(torus->aabb);
    aabb->min[0] = -(torus->swept_radius + torus->tube_radius);
    aabb->min[1] = -(torus->tube_radius);
    aabb->min[2] = aabb->min[0];

    aabb->max[0] = torus->swept_radius + torus->tube_radius;
    aabb->max[1] = torus->tube_radius;
    aabb->max[2] = aabb->max[0];
}


float rayIntersectGenericTorus(ShadeRec* sr, GenericTorus* torus, const Ray ray)
{
    double x1 = ray.origin[0]; double y1 = ray.origin[1]; double z1 = ray.origin[2];
    double d1 = ray.direction[0]; double d2 = ray.direction[1]; double d3 = ray.direction[2];
    double coeffs[5];
    double roots[4];    

    double sum_d_sqrd = d1*d1 + d2*d2 + d3*d3;
    double e = x1*x1 + y1*y1 + z1*z1
        - torus->swept_radius*torus->swept_radius - torus->tube_radius*torus->tube_radius;
    double f = x1*d1 + y1*d2 + z1*d3;
    double four_a_sqrd = 4.0f * torus->swept_radius * torus->swept_radius;

    coeffs[0] = e*e - four_a_sqrd*(torus->tube_radius * torus->tube_radius - y1*y1);
    coeffs[1] = 4.0f*f*e + 2.0f*four_a_sqrd*y1*d2;
    coeffs[2] = 2.0f*sum_d_sqrd*e + 4.0f*f*f + four_a_sqrd*d2*d2;
    coeffs[3] = 4.0f*sum_d_sqrd*f;
    coeffs[4] = sum_d_sqrd*sum_d_sqrd;

    int num_real_roots = solveQuartic(coeffs, roots);

    bool intersected = false;
    float t = FLT_MAX;

    if(num_real_roots == 0)
    {
        return TMAX;
    }

    for(int i = 0; i < num_real_roots; i++)
    {
        if(roots[i] > K_EPSILON)
        {
            intersected = true;
            if(roots[i] < t)
            {
                t = (float)roots[i];
            }
        }
    }

    if(!intersected)
    {
        return TMAX;
    }
    vec3 hit_point;
    getPointOnRay(hit_point, ray, t);
    float phi = (float)atan2(hit_point[0], hit_point[2]);
    vec3 horizontal_comp = {hit_point[0], 0.0f, hit_point[2]};
    float x = vec3_length(horizontal_comp) - torus->swept_radius;
    float theta = (float)atan2(hit_point[1], x);
    if(fabs(phi) <= torus->phi)
    {
        vec3_copy(sr->hit_point, hit_point);
        computeGenericTorusNormal(sr->normal, torus, sr->hit_point);
        vec3_negate(sr->wo, ray.direction);
        if(vec3_dot(sr->wo, sr->normal) < 0)
        {
            vec3_negate(sr->normal, sr->normal);
        }
        sr->mat = torus->mat;
        return (float)t;                                
    }

    return (float)t;
}

float shadowRayIntersectGenericTorus(const GenericTorus* torus, const Ray ray)
{
    float x1 = ray.origin[0]; float y1 = ray.origin[1]; float z1 = ray.origin[2];
    float d1 = ray.direction[0]; float d2 = ray.direction[1]; float d3 = ray.direction[2];
    double coeffs[5];
    double roots[4];    

    float sum_d_sqrd = d1*d1 + d2*d2 + d3*d3;
    float e = x1*x1 + y1*y1 + z1*z1
        - torus->swept_radius*torus->swept_radius - torus->tube_radius*torus->tube_radius;
    float f = x1*d1 + y1*d2 + z1*d3;
    float four_a_sqrd = 4.0f * torus->swept_radius * torus->swept_radius;    

    coeffs[0] = e*e - four_a_sqrd*(torus->tube_radius * torus->tube_radius - y1*y1);    
    coeffs[1] = 4.0f*f*e + 2.0f*four_a_sqrd*y1*d2;
    coeffs[2] = 2.0f*sum_d_sqrd*e + 4.0f*f*f + four_a_sqrd*d2*d2;
    coeffs[3] = 4.0f*sum_d_sqrd*f;
    coeffs[4] = sum_d_sqrd*sum_d_sqrd;

    int num_real_roots = solveQuartic(coeffs, roots);

    bool intersected = false;
    float t = FLT_MAX;

    if(num_real_roots == 0)
    {
        return TMAX;
    }

    for(int i = 0; i < num_real_roots; i++)
    {
        if(roots[i] > K_EPSILON)
        {
            intersected = true;
            if(roots[i] < t)
            {
                t = (float)roots[i];
            }
        }
    }
    if(!intersected)
    {
        return TMAX;
    }

    vec3 hit_point;
    getPointOnRay(hit_point, ray, t);
    float phi = (float)atan2(hit_point[0], hit_point[2]);
    float theta = (float)acos(min(hit_point[1] / torus->tube_radius, 1.0f));
    if(fabs(phi) <= torus->phi)        
    {
        return t;
    }

    return TMAX;
}
