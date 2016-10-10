#include "aabb.h"
/*
float rayIntersectAABB(const AABB* aabb, const Ray ray)
{
    float tx_min, ty_min, tz_min;
    float tx_max, ty_max, tz_max;

    float a = 1.0f/ray.direction[0];
    if(a >= 0)
    {
        tx_min = (aabb->min[0] - ray.origin[0]) * a;
        tx_max = (aabb->max[0] - ray.origin[0]) * a;
    }else
    {
        tx_min = (aabb->max[0] - ray.origin[0]) * a;
        tx_max = (aabb->min[0] - ray.origin[0]) * a;        
    }

    float b = 1.0f/ray.direction[1];
    if(b >= 0)
    {
        ty_min = (aabb->min[1] - ray.origin[1]) * b;
        ty_max = (aabb->max[1] - ray.origin[1]) * b;        
    }else
    {
        ty_min = (aabb->max[1] - ray.origin[1]) * b;
        ty_max = (aabb->min[1] - ray.origin[1]) * b;        
    }

    float c = 1.0f/ray.direction[2];
    if(c >= 0)
    {
        tz_min = (aabb->min[2] - ray.origin[2]) * c;
        tz_max = (aabb->max[2] - ray.origin[2]) * c;
    }else
    {
        tz_min = (aabb->max[2] - ray.origin[2]) * c;
        tz_max = (aabb->min[2] - ray.origin[2]) * c;
    }
    
    float t0, t1;

    if(tx_min > ty_min)
    {
        t0 = tx_min;
    }else
    {
        t0 = ty_min;
    }
    if(tz_min > t0)
    {
        t0 = tz_min;
    }

    //t0 = fmax(fmax(tx_min, ty_min), tz_min);    


    if(tx_max < ty_max)
    {
        t1 = tx_max;
    }else
    {
        t1 = ty_max;
    }
    if(tz_max < t1)
    {
        t1 = tz_max;
    }

    //t1 = fmin(fmin(tx_max, ty_max), tz_min);

    if(t0 < t1 && t1 > K_EPSILON)
    {
        if(t0 > K_EPSILON)
        {
            return t0;
        }else
        {
            return t1;
        }
    }else
    {
        return TMAX;
    }
    //return (t0 < t1 && t1 > K_EPSILON);
}
*/

float rayIntersectAABB(const AABB* aabb, const Ray ray)
{
    float tx_min, ty_min, tz_min;
    float tx_max, ty_max, tz_max;

    float a = 1.0f/ray.direction[0];
    if(a >= 0)
    {
        tx_min = (aabb->min[0] - ray.origin[0]) * a;
        tx_max = (aabb->max[0] - ray.origin[0]) * a;
    }else
    {
        tx_min = (aabb->max[0] - ray.origin[0]) * a;
        tx_max = (aabb->min[0] - ray.origin[0]) * a;        
    }

    float b = 1.0f/ray.direction[1];
    if(b >= 0)
    {
        ty_min = (aabb->min[1] - ray.origin[1]) * b;
        ty_max = (aabb->max[1] - ray.origin[1]) * b;        
    }else
    {
        ty_min = (aabb->max[1] - ray.origin[1]) * b;
        ty_max = (aabb->min[1] - ray.origin[1]) * b;        
    }

    float c = 1.0f/ray.direction[2];
    if(c >= 0)
    {
        tz_min = (aabb->min[2] - ray.origin[2]) * c;
        tz_max = (aabb->max[2] - ray.origin[2]) * c;
    }else
    {
        tz_min = (aabb->max[2] - ray.origin[2]) * c;
        tz_max = (aabb->min[2] - ray.origin[2]) * c;
    }
    
    float t0, t1;
    if(tx_min > ty_min)
    {
        t0 = tx_min;
    }else
    {
        t0 = ty_min;
    }

    if(tz_min > t0)
    {
        t0 = tz_min;
    }

    if(tx_max < ty_max)
    {
        t1 = tx_max;
    }else
    {
        t1 = ty_max;
    }

    if(tz_max < t1)
    {
        t1 = tz_max;
    }

    if(t0 < t1 && t1 > K_EPSILON)
    {
        if(t0 > K_EPSILON)
        {
            return t0;
        }else
        {
            return t1;
        }
    }else
    {
        return TMAX;
    }
    //return (t0 < t1 && t1 > K_EPSILON);
}


bool isInsideAABB(const AABB* aabb, const vec3 point)
{
    for(int i = 0; i < 3; i++)
    {
        if(point[i] < aabb->min[i] || point[i] > aabb->max[i])
        {
            return false;
        }
    }
    return true;
}

void addToAABB(AABB* r, const AABB* a)
{
    if(a->min[0] < r->min[0]){r->min[0] = a->min[0];}
    if(a->min[1] < r->min[1]){r->min[1] = a->min[1];}
    if(a->min[2] < r->min[2]){r->min[2] = a->min[2];}

    if(a->max[0] > r->max[0]){r->max[0] = a->max[0];}
    if(a->max[1] > r->max[1]){r->max[1] = a->max[1];}
    if(a->max[2] > r->max[2]){r->max[2] = a->max[2];}        
}
   