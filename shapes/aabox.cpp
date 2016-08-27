#include "aabox.h"

void getAABoxNormal(vec3 r, const int face_hit)
{
    switch(face_hit)
    {
    case 0:
    {
        vec3_assign(r, -1.0f, 0.0f, 0.0f);
    } break;
    case 1:
    {
        vec3_assign(r, 0.0f, -1.0f, 0.0f);
    } break;
    case 2:
    {
        vec3_assign(r, 0.0f, 0.0f, -1.0f);    // NOTE: z component is negated due to negative z being forward
    } break;
    case 3:
    {
        vec3_assign(r, 1.0f, 0.0f, 0.0f);
    } break;
    case 4:
    {
        vec3_assign(r, 0.0f, 1.0f, 0.0f);
    } break;
    case 5:
    {
        vec3_assign(r, 0.0f, 0.0f, 1.0f);    // NOTE: z component is negated due to negative z being forward
    } break;    
    }
}

float rayIntersectAABox(ShadeRec* sr, AABox* aabox, const Ray ray)
{    
    float tx_min, ty_min, tz_min;
    float tx_max, ty_max, tz_max;
    
    // NOTE: put the code below into a fuction
    float a = 1.0f/ray.direction[0];
    if(a >= 0)
    {
        tx_min = (aabox->min[0] - ray.origin[0]) * a;
        tx_max = (aabox->max[0] - ray.origin[0]) * a;
    }else
    {
        tx_min = (aabox->max[0] - ray.origin[0]) * a;
        tx_max = (aabox->min[0] - ray.origin[0]) * a;        
    }

    float b = 1.0f/ray.direction[1];
    if(b >= 0)
    {
        ty_min = (aabox->min[1] - ray.origin[1]) * b;
        ty_max = (aabox->max[1] - ray.origin[1]) * b;        
    }else
    {
        ty_min = (aabox->max[1] - ray.origin[1]) * b;
        ty_max = (aabox->min[1] - ray.origin[1]) * b;        
    }
    // NOTE: z component is negated due to negative z being forward
    float c = 1.0f/ray.direction[2];
    if(c >= 0)
    {
        tz_min = (aabox->min[2] - ray.origin[2]) * c;
        tz_max = (aabox->max[2] - ray.origin[2]) * c;
    }else
    {
        tz_min = (aabox->max[2] - ray.origin[2]) * c;
        tz_max = (aabox->min[2] - ray.origin[2]) * c;
    }

    float t0, t1;
    int face_in, face_out;
    if(tx_min > ty_min)
    {
        t0 = tx_min;
        face_in = (a >= 0.0f) ? 0 : 3;
    }else
    {
        t0 = ty_min;
        face_in = (b >= 0.0f) ? 1 : 4;        
    }
    if(tz_min > t0)
    {
        t0 = tz_min;
        face_in = (c >= 0.0f) ? 2 : 5;                
    }

    if(tx_max < ty_max)
    {
        t1 = tx_max;
        face_out = (a >= 0.0f) ? 3 : 0;
    }else
    {
        t1 = ty_max;
        face_out = (b >= 0.0f) ? 4 : 1;        
    }
    if(tz_max < t1)
    {
        t1 = tz_max;
        face_out = (c >= 0.0f) ? 5 : 2;                
    }

    float t_min = TMAX;
    if(t0 < t1 && t1 > K_EPSILON)
    {
        if(t0 > K_EPSILON)
        {
            t_min = t0;
            getAABoxNormal(sr->normal, face_in);
        }else
        {
            t_min = t1;
            getAABoxNormal(sr->normal, face_out);
        }
        if(vec3_dot(sr->normal, ray.direction) > 0.0f)
        {
            vec3_negate(sr->normal, sr->normal);
        }
        getPointOnRay(sr->hit_point, ray, t_min);
        sr->mat = aabox->mat;
        sr->hit_status = true;
        vec3_negate(sr->wo, ray.direction);
    }
    return t_min;
}

float shadowRayIntersectAABox(const AABox* aabox, const Ray ray)
{
    if(!aabox->shadow)
    {
        return TMAX;
    }
    float tx_min, ty_min, tz_min;
    float tx_max, ty_max, tz_max;

    float a = 1.0f/ray.direction[0];
    if(a >= 0)
    {
        tx_min = (aabox->min[0] - ray.origin[0]) * a;
        tx_max = (aabox->max[0] - ray.origin[0]) * a;
    }else
    {
        tx_min = (aabox->max[0] - ray.origin[0]) * a;
        tx_max = (aabox->min[0] - ray.origin[0]) * a;        
    }

    float b = 1.0f/ray.direction[1];
    if(b >= 0)
    {
        ty_min = (aabox->min[1] - ray.origin[1]) * b;
        ty_max = (aabox->max[1] - ray.origin[1]) * b;        
    }else
    {
        ty_min = (aabox->max[1] - ray.origin[1]) * b;
        ty_max = (aabox->min[1] - ray.origin[1]) * b;        
    }

    float c = 1.0f/ray.direction[2];
    if(c >= 0)
    {
        tz_min = (aabox->min[2] - ray.origin[2]) * c;
        tz_max = (aabox->max[2] - ray.origin[2]) * c;
    }else
    {
        tz_min = (aabox->max[2] - ray.origin[2]) * c;
        tz_max = (aabox->min[2] - ray.origin[2]) * c;
    }

    float t0, t1;
    int face_in, face_out;
    if(tx_min > ty_min)
    {
        t0 = tx_min;
        face_in = (a >= 0.0f) ? 0 : 3;
    }else
    {
        t0 = ty_min;
        face_in = (b >= 0.0f) ? 1 : 4;        
    }
    if(tz_min > t0)
    {
        t0 = tz_min;
        face_in = (c >= 0.0f) ? 2 : 5;                
    }

    if(tx_max < ty_max)
    {
        t1 = tx_max;
        face_out = (a >= 0.0f) ? 3 : 0;
    }else
    {
        t1 = ty_max;
        face_out = (b >= 0.0f) ? 4 : 1;        
    }
    if(tz_max < t1)
    {
        t1 = tz_max;
        face_out = (c >= 0.0f) ? 5 : 2;                
    }

    float t_min = TMAX;
    if(t0 < t1 && t1 > K_EPSILON)
    {
        if(t0 > K_EPSILON)
        {
            t_min = t0;
        }else
        {
            t_min = t1;
        }
    }
    return t_min;

}

