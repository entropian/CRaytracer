#include "shapes.h"

float rayIntersectObject(ShadeRec* sr, const Object_t obj, const Ray ray)
{
    float t = TMAX;
    switch(obj.type)
    {
    case SPHERE:
        t = rayIntersectSphere(sr, (Sphere*)obj.ptr, ray);
        break;
    case PLANE:
        t = rayIntersectPlane(sr, (Plane*)obj.ptr, ray);
        break;
    case RECTANGLE:
        t = rayIntersectRect(sr, (Rectangle*)obj.ptr, ray);
        break;
    case AABOX:
        t = rayIntersectAABox(sr, (AABox*)obj.ptr, ray);
        break;
    case TRIANGLE:
        t = rayIntersectTriangle(sr, (Triangle*)obj.ptr, ray);
        break;
    case GENERICOPENCYLINDER:
        t = rayIntersectGenericOpenCylinder(sr, (GenericOpenCylinder*)obj.ptr, ray);
        break;
    case DISK:
        t = rayIntersectDisk(sr, (Disk*)obj.ptr, ray);
        break;
    case GENERICTORUS:
        t = rayIntersectGenericTorus(sr, (GenericTorus*)obj.ptr, ray);
        break;
    case INSTANCED:
        t = rayIntersectInstanced(sr, (InstancedShape*)obj.ptr, ray);
        break;
    case COMPOUND:
        t = rayIntersectCompound(sr, (CompoundObject*)obj.ptr, ray);
        break;
    case FLAT_TRIANGLE:
        t = rayIntersectFlatTriangle(sr, (FlatTriangle*)obj.ptr, ray);
        break;
    case SMOOTH_TRIANGLE:
        t = rayIntersectSmoothTriangle(sr, (SmoothTriangle*)obj.ptr, ray);
        break;                
    }
    return t;
}

float shadowRayIntersectObject(const Object_t obj, const Ray ray)
{
    float t = TMAX;
    switch(obj.type)
    {
    case SPHERE:
        t = shadowRayIntersectSphere((Sphere*)obj.ptr, ray);
        break;
    case PLANE:
        t = shadowRayIntersectPlane((Plane*)obj.ptr, ray);
        break;
    case RECTANGLE:
        t = shadowRayIntersectRect((Rectangle*)obj.ptr, ray);
        break;
    case AABOX:
        t = shadowRayIntersectAABox((AABox*)obj.ptr, ray);
        break;
    case TRIANGLE:
        t = shadowRayIntersectTriangle((Triangle*)obj.ptr, ray);
        break;
    case GENERICOPENCYLINDER:
        t = shadowRayIntersectGenericOpenCylinder((GenericOpenCylinder*)obj.ptr, ray);
        break;
    case DISK:
        t = shadowRayIntersectDisk((Disk*)obj.ptr, ray);
        break;
    case GENERICTORUS:
        t = shadowRayIntersectGenericTorus((GenericTorus*)obj.ptr, ray);
        break;
    case INSTANCED:
        t = shadowRayIntersectInstanced((InstancedShape*)obj.ptr, ray);
        break;
    case COMPOUND:
        t = shadowRayIntersectCompound((CompoundObject*)obj.ptr, ray);
        break;
    case FLAT_TRIANGLE:
        t = shadowRayIntersectFlatTriangle((FlatTriangle*)obj.ptr, ray);
        break;
    case SMOOTH_TRIANGLE:
        t = shadowRayIntersectSmoothTriangle((SmoothTriangle*)obj.ptr, ray);
        break;                        
    }
    return t;
}

bool isGridObjType(const Object_t obj)
{
    switch(obj.type)
    {
    case SPHERE:
    case RECTANGLE:
    case AABOX:
    case TRIANGLE:
    case GENERICOPENCYLINDER:
    case DISK:
    case GENERICTORUS:
    case FLAT_TRIANGLE:
    case SMOOTH_TRIANGLE:        
        return true;
        break;
    case INSTANCED:
    {
        InstancedShape* is = (InstancedShape*)(obj.ptr);
        return isGridObjType(is->obj);
    } break;
    }
    return false;
}

void updateAABB(AABB *out, const mat4 transform, const AABB a)
{
    for(int i = 0; i < 3; i++)
    {
        out->min[i] = out->max[i] = transform[i][3];
        for(int j = 0; j < 3; j++)
        {
            float e = transform[j][i] * a.min[j];
            float f = transform[j][i] * a.max[j];
            if(e < f)
            {
                out->min[i] += e;
                out->max[i] += f; 
            }else
            {
                out->min[i] += f;
                out->max[i] += e; 
            }
        }
    }
}

void getObjectAABB(AABB* aabb, const Object_t obj)
{
    switch(obj.type)
    {
    case PLANE:
      break;
    case SPHERE:
    {
        Sphere* sphere = (Sphere*)obj.ptr;
        aabb->min[0] = sphere->center[0] - sphere->radius;
        aabb->min[1] = sphere->center[1] - sphere->radius;
        aabb->min[2] = sphere->center[2] - sphere->radius;

        aabb->max[0] = sphere->center[0] + sphere->radius;
        aabb->max[1] = sphere->center[1] + sphere->radius;
        aabb->max[2] = sphere->center[2] + sphere->radius;
    } break;
    case RECTANGLE:
    {
        // NOTE: shit code
        Rectangle* rect = (Rectangle*)obj.ptr;
        vec3 pts[4];
        vec3_copy(pts[0], rect->point);
        vec3_add(pts[1], pts[0], rect->width);
        vec3_add(pts[2], pts[0], rect->height);
        vec3 tmp;
        vec3_add(tmp, rect->width, rect->height);
        vec3_add(pts[3], pts[0], tmp);
            
        vec3_copy(aabb->min, pts[0]);
        vec3_copy(aabb->max, pts[0]);

        for(int i = 1; i < 4; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                if(pts[i][j] < aabb->min[j]){aabb->min[j] = pts[i][j];}
                if(pts[i][j] > aabb->max[j]){aabb->max[j] = pts[i][j];}
            }
        }
        for(int i = 0; i < 3; i++)
        {
            aabb->min[i] -= K_EPSILON;
            aabb->max[i] += K_EPSILON;            
        }
    } break;
    case AABOX:
    {
        AABox* aabox = (AABox*)obj.ptr;
        float gap = 0.00001f;
        vec3_assign(aabb->min, aabox->min[0] - gap, aabox->min[1] - gap, aabox->min[2] - gap);
        vec3_assign(aabb->max, aabox->max[0] + gap, aabox->max[1] + gap, aabox->max[2] + gap);        
    } break;
    case TRIANGLE:
    {
        Triangle* triangle = (Triangle*)obj.ptr;
        vec3_copy(aabb->min, triangle->v0);
        vec3_copy(aabb->max, triangle->v0);

        for(int i = 0; i < 3; i++)
        {
            if(triangle->v1[i] < aabb->min[i]){aabb->min[i] = triangle->v1[i];}
            if(triangle->v1[i] > aabb->max[i]){aabb->max[i] = triangle->v1[i];}
            if(triangle->v2[i] < aabb->min[i]){aabb->min[i] = triangle->v2[i];}
            if(triangle->v2[i] > aabb->max[i]){aabb->max[i] = triangle->v2[i];}
        }
        for(int i = 0; i < 3; i++)
        {
            aabb->min[i] -= K_FLAT_AABB;
            aabb->max[i] += K_FLAT_AABB;
        }
    } break;
    case FLAT_TRIANGLE:
    {
        FlatTriangle* mesh_tri = (FlatTriangle*)obj.ptr;
        vec3_copy(aabb->min, mesh_tri->v0);
        vec3_copy(aabb->max, mesh_tri->v0);

        for(int i = 0; i < 3; i++)
        {
            if(mesh_tri->v1[i] < aabb->min[i]){aabb->min[i] = mesh_tri->v1[i];}
            if(mesh_tri->v1[i] > aabb->max[i]){aabb->max[i] = mesh_tri->v1[i];}
            if(mesh_tri->v2[i] < aabb->min[i]){aabb->min[i] = mesh_tri->v2[i];}
            if(mesh_tri->v2[i] > aabb->max[i]){aabb->max[i] = mesh_tri->v2[i];}
        }
        for(int i = 0; i < 3; i++)
        {
            aabb->min[i] -= K_FLAT_AABB;
            aabb->max[i] += K_FLAT_AABB;
        }
    } break;
    case SMOOTH_TRIANGLE:
    {
        SmoothTriangle* mesh_tri = (SmoothTriangle*)obj.ptr;
        vec3_copy(aabb->min, mesh_tri->v0);
        vec3_copy(aabb->max, mesh_tri->v0);

        for(int i = 0; i < 3; i++)
        {
            if(mesh_tri->v1[i] < aabb->min[i]){aabb->min[i] = mesh_tri->v1[i];}
            if(mesh_tri->v1[i] > aabb->max[i]){aabb->max[i] = mesh_tri->v1[i];}
            if(mesh_tri->v2[i] < aabb->min[i]){aabb->min[i] = mesh_tri->v2[i];}
            if(mesh_tri->v2[i] > aabb->max[i]){aabb->max[i] = mesh_tri->v2[i];}
        }
        for(int i = 0; i < 3; i++)
        {
            aabb->min[i] -= K_FLAT_AABB;
            aabb->max[i] += K_FLAT_AABB;
        }
    } break;        
    case GENERICOPENCYLINDER:
    {
        GenericOpenCylinder* oc = (GenericOpenCylinder*)obj.ptr;
        aabb->min[0] = -(oc->radius);
        aabb->min[1] = -(oc->half_height);
        aabb->min[2] = -(oc->radius);
        aabb->max[0] = oc->radius;
        aabb->max[1] = oc->half_height;
        aabb->max[2] = oc->radius;
    } break;
    case DISK:
    {
        // TODO: improve this crap
        Disk* disk = (Disk*)obj.ptr;
        aabb->min[0] = disk->center[0] - disk->radius;
        aabb->min[1] = disk->center[1] - disk->radius;
        aabb->min[2] = disk->center[2] - disk->radius;

        aabb->max[0] = disk->center[0] + disk->radius;
        aabb->max[1] = disk->center[1] + disk->radius;
        aabb->max[2] = disk->center[2] + disk->radius;        
    } break;
    case GENERICTORUS:
    {
        GenericTorus* torus = (GenericTorus*)obj.ptr;
        //*aabb = torus->aabb;
        calcAABBGenericTorus(aabb, torus);
    } break;
    case INSTANCED:
    {
        InstancedShape* is = (InstancedShape*)obj.ptr;
        getObjectAABB(aabb, is->obj);
        AABB tmp;
        vec3_copy(tmp.min, aabb->min);
        vec3_copy(tmp.max, aabb->max);        
        mat4 transform;
        affine_inverse(transform, is->inv_transform);
        /*
        vec4 tmp1, tmp2;
        vec3 tmp_min, tmp_max;
        vec4_assign(tmp1, aabb->min[0], aabb->min[1], aabb->min[2], 1.0f);
        mat4_mult_vec4(tmp2, transform, tmp1);
        vec3_assign(aabb->min, tmp2[0], tmp2[1], tmp2[2]);
        vec4_assign(tmp1, aabb->max[0], aabb->max[1], aabb->max[2], 1.0f);
        mat4_mult_vec4(tmp2, transform, tmp1);
        vec3_assign(aabb->max, tmp2[0], tmp2[1], tmp2[2]);
        */
        updateAABB(aabb, transform, tmp);
    } break;
    case COMPOUND:
    {
        CompoundObject* co = (CompoundObject*)obj.ptr;
        *aabb = co->aabb;
    } break;
    }
}

Material* getObjectMatPtr(const Object_t obj)
{
    Material *mat = NULL;
    switch(obj.type)
    {
    case SPHERE:
        mat = ((Sphere*)(obj.ptr))->mat;
        break;
    case RECTANGLE:
        mat = ((Rectangle*)(obj.ptr))->mat;
        break;
    case AABOX:
        mat = ((AABox*)(obj.ptr))->mat;
        break;
    case TRIANGLE:
        mat = ((Triangle*)(obj.ptr))->mat;
        break;
    case GENERICOPENCYLINDER:
        mat = ((GenericOpenCylinder*)(obj.ptr))->mat;        
        break;
    case DISK:
        mat = ((Disk*)(obj.ptr))->mat;
        break;
    case GENERICTORUS:
        mat = ((GenericTorus*)(obj.ptr))->mat;
        break;
    case FLAT_TRIANGLE:
        mat = ((FlatTriangle*)(obj.ptr))->mat;
        break;
    case SMOOTH_TRIANGLE:
        mat = ((SmoothTriangle*)(obj.ptr))->mat;        
        break;
    case INSTANCED:
        mat = ((InstancedShape*)(obj.ptr))->mat;        
        break;
        // TODO: compound objects
    }
    return mat;
}


bool calcBoundingSphere(vec3 center, float *radius, const Object_t obj)
{    
    switch(obj.type)
    {
    case SPHERE:
    {
        Sphere *sphere = (Sphere*)(obj.ptr);
        *radius = sphere->radius + K_EPSILON;
        vec3_copy(center, sphere->center);
        return true;
    } break;
    case RECTANGLE:
    {
        Rectangle *rect = (Rectangle*)(obj.ptr);
        vec3 displacement, tmp;
        vec3_scale(displacement, rect->width, 0.5f);
        vec3_scale(tmp, rect->height, 0.5f);
        vec3_add(displacement, displacement, tmp);
        vec3_add(center, rect->point, displacement);
        *radius = vec3_length(displacement) + K_EPSILON;
        return true;
    } break;
    case AABOX:
        // TODO
        break;
    case GENERICOPENCYLINDER:        
    case TRIANGLE:
    case FLAT_TRIANGLE:
    case SMOOTH_TRIANGLE:
    case INSTANCED:
    case COMPOUND:        
    {
        AABB aabb;
        getObjectAABB(&aabb, obj);
        vec3 cross_vec;
        vec3_sub(cross_vec, aabb.max, aabb.min);
        vec3_scale(cross_vec, cross_vec, 0.5f);
        vec3_add(center, aabb.min, cross_vec);
        *radius = vec3_length(cross_vec) + K_EPSILON;
        return true;
    } break;
    case DISK:
    {
        Disk *disk = (Disk*)(obj.ptr);
        vec3_copy(center, disk->center);
        *radius = disk->radius + K_EPSILON;
        return true;
    } break;
    case GENERICTORUS:
    {
        GenericTorus *gt = (GenericTorus*)(obj.ptr);
        vec3_assign(center, 0.0f, 0.0f, 0.0f);
        *radius = gt->swept_radius + gt->tube_radius;
        return true;
    } break;
    }
    return false;
}
