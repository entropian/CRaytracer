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

bool getObjectAABB(AABB* aabb, const Object_t obj)
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
    default:
        return false;
    }
    return true;
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

int testAABBPlane(AABB* aabb, vec3 plane_normal, float plane_d)
{
    vec3 span;
    vec3_sub(span, aabb->max, aabb->min);
    //vec3 center = (box.max - box.min) * 0.5f + box.min;
    vec3 half_span;
    vec3_scale(half_span, span, 0.5f);
    vec3 center;
    vec3_add(center, aabb->min, half_span);
    
    //float r = e[0] * fabs(plane_normal[0]) + e[1] * fabs(plane_normal[1]) + e[2] * fabs(plane_normal[2]);    
    float r = half_span[0] * fabs(plane_normal[0]) + half_span[1] * fabs(plane_normal[1])
        + half_span[2] * fabs(plane_normal[2]);
    float s = vec3_dot(plane_normal, center) - plane_d;
    return fabs(s) <= r;
}
 
//int testTriangleBox(const Triangle& triangle, const Box& box, bool& hasZeroVector)
int testTriangleAABB(vec3 tv0, vec3 tv1, vec3 tv2, AABB* aabb, int has_zero_vector)
{
    //Vec3 center = (box.max - box.min) * 0.5f + box.min;
    vec3 span;
    vec3_sub(span, aabb->max, aabb->min);
    vec3 half_span;
    vec3_scale(half_span, span, 0.5f);
    vec3 center;
    vec3_add(center, half_span, aabb->min);
    //const float x_extent = (box.max[0] - box.min[0]) * 0.5f;
    //const float y_extent = (box.max[1] - box.min[1]) * 0.5f;
    //const float z_extent = (box.max[2] - box.min[2]) * 0.5f;
    float x_extent = half_span[0];
    float y_extent = half_span[1];
    float z_extent = half_span[2];
 
    //Vec3 v0 = triangle.v0 - center;
    //Vec3 v1 = triangle.v1 - center;
    //Vec3 v2 = triangle.v2 - center;    
    vec3 v0, v1, v2;
    vec3_sub(v0, tv0, center);
    vec3_sub(v1, tv1, center);
    vec3_sub(v2, tv2, center);
 
    vec3 u0 = {1.0f, 0.0f, 0.0f};
    vec3 u1 = {0.0f, 1.0f, 0.0f};
    vec3 u2 = {0.0f, 0.0f, 1.0f};
 
    //Vec3 edge_0 = v1 - v0, edge_1 = v2 - v1, edge_2 = v0 - v2;
    vec3 e0, e1, e2;
    vec3_sub(e0, v1, v0);
    vec3_sub(e1, v2, v1);
    vec3_sub(e2, v0, v2);
 
    vec3 a00, a01, a02;
    //Vec3 a00 = cross(u0, edge_0);
    vec3_cross(a00, u0, e0);
    has_zero_vector |= isZeroVector(a00);
    //Vec3 a01 = cross(u0, edge_1);
    vec3_cross(a01, u0, e1);
    has_zero_vector |= isZeroVector(a01);
    //Vec3 a02 = cross(u0, edge_2);
    vec3_cross(a02, u0, e2);
    has_zero_vector |= isZeroVector(a02);

    vec3 a10, a11, a12;
    //Vec3 a10 = cross(u1, edge_0);
    vec3_cross(a10, u1, e0);
    has_zero_vector |= isZeroVector(a10);
    //Vec3 a11 = cross(u1, edge_1);
    vec3_cross(a11, u1, e1);
    has_zero_vector |= isZeroVector(a11);
    //Vec3 a12 = cross(u1, edge_2);
    vec3_cross(a12, u1, e2);
    has_zero_vector |= isZeroVector(a12);

    vec3 a20, a21, a22;
    //Vec3 a20 = cross(u2, edge_0);
    vec3_cross(a20, u2, e0);
    has_zero_vector |= isZeroVector(a20);
    //Vec3 a21 = cross(u2, edge_1);
    vec3_cross(a21, u2, e1);
    has_zero_vector |= isZeroVector(a21);
    //Vec3 a22 = cross(u2, edge_2);
    vec3_cross(a22, u2, e2);
    has_zero_vector |= isZeroVector(a22);
 
    float p0 = vec3_dot(v0, a00);
    float p1 = vec3_dot(v1, a00);
    float p2 = vec3_dot(v2, a00);
    float r = x_extent * fabs(a00[0]) + y_extent * fabs(a00[1]) + z_extent * fabs(a00[2]);
    if (max(-max(p0, p2), min(p0, p2)) > r) return 0;
 
    p0 = vec3_dot(v0, a01);
    p1 = vec3_dot(v1, a01);
    p2 = vec3_dot(v2, a01);
    r = x_extent * fabs(a01[0]) + y_extent * fabs(a01[1]) + z_extent * fabs(a01[2]);
    if (max(-max3(p0, p1, p2), min3(p0, p1, p2)) > r) return false;
 
    p0 = vec3_dot(v0, a02);
    p1 = vec3_dot(v1, a02);
    p2 = vec3_dot(v2, a02);
    r = x_extent * fabs(a02[0]) + y_extent * fabs(a02[1]) + z_extent * fabs(a02[2]);
    if (max(-max3(p0, p1, p2), min3(p0, p1, p2)) > r) return false;
 
    p0 = vec3_dot(v0, a10);
    p1 = vec3_dot(v1, a10);
    p2 = vec3_dot(v2, a10);
    r = x_extent * fabs(a10[0]) + y_extent * fabs(a10[1]) + z_extent * fabs(a10[2]);
    if (max(-max3(p0, p1, p2), min3(p0, p1, p2)) > r) return false;
 
    p0 = vec3_dot(v0, a11);
    p1 = vec3_dot(v1, a11);
    p2 = vec3_dot(v2, a11);
    r = x_extent * fabs(a11[0]) + y_extent * fabs(a11[1]) + z_extent * fabs(a11[2]);
    if (max(-max3(p0, p1, p2), min3(p0, p1, p2)) > r) return false;
 
    p0 = vec3_dot(v0, a12);
    p1 = vec3_dot(v1, a12);
    p2 = vec3_dot(v2, a12);
    r = x_extent * fabs(a12[0]) + y_extent * fabs(a12[1]) + z_extent * fabs(a12[2]);
    if (max(-max3(p0, p1, p2), min3(p0, p1, p2)) > r) return false;
 
    p0 = vec3_dot(v0, a20);
    p1 = vec3_dot(v1, a20);
    p2 = vec3_dot(v2, a20);
    r = x_extent * fabs(a20[0]) + y_extent * fabs(a20[1]) + z_extent * fabs(a20[2]);
    if (max(-max3(p0, p1, p2), min3(p0, p1, p2)) > r) return false;
 
    p0 = vec3_dot(v0, a21);
    p1 = vec3_dot(v1, a21);
    p2 = vec3_dot(v2, a21);
    r = x_extent * fabs(a21[0]) + y_extent * fabs(a21[1]) + z_extent * fabs(a21[2]);
    if (max(-max3(p0, p1, p2), min3(p0, p1, p2)) > r) return false;
 
    p0 = vec3_dot(v0, a22);
    p1 = vec3_dot(v1, a22);
    p2 = vec3_dot(v2, a22);
    r = x_extent * fabs(a22[0]) + y_extent * fabs(a22[1]) + z_extent * fabs(a22[2]);
    if (max(-max3(p0, p1, p2), min3(p0, p1, p2)) > r) return false;
 
    if (max3(v0[0], v1[0], v2[0]) < -x_extent || min3(v0[0], v1[0], v2[0]) > x_extent) return false;
    if (max3(v0[1], v1[1], v2[1]) < -y_extent || min3(v0[1], v1[1], v2[1]) > y_extent) return false;
    if (max3(v0[2], v1[2], v2[2]) < -z_extent || min3(v0[2], v1[2], v2[2]) > z_extent) return false;
 
    //Vec3 plane_normal = normalize(vec3_cross(edge_0, edge_1));
    vec3 plane_normal;
    vec3_cross(plane_normal, e0, e1);
    vec3_normalize(plane_normal, plane_normal);
    float plane_d = vec3_dot(plane_normal, v0);
    return testAABBPlane(aabb, plane_normal, plane_d);
}
