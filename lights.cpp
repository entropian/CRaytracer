#include "lights.h"
#include "util/math.h"

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

void getLightDir(vec3 r, const LightType light_type, const void* light_ptr, const ShadeRec* sr, const vec2 sample)
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
            mapSampleToHemisphere(h_sample, sample);
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
        //getSample3D(h_sample, env_light_ptr->samples3D, sample_index);
        mapSampleToHemisphere(h_sample, sample);
        getVec3InLocalBasis(r, h_sample, sr->normal);        
    } break;
    }
}

void getIncRadiance(vec3 r, const LightType light_type, const void* light_ptr, const vec3 hit_point)
{
    switch(light_type)
    {
    case DIRECTIONAL:
    {
        vec3_scale(r, ((DirLight*)light_ptr)->color, ((DirLight*)light_ptr)->intensity);
    } break;
    case POINTLIGHT:
    {
        PointLight *point_light = (PointLight*)light_ptr;
        vec3_scale(r, point_light->color, point_light->intensity);
        if(point_light->dist_atten)
        {
            vec3 dist;
            vec3_sub(dist, point_light->point, hit_point);
            float length = vec3_length(dist);
            vec3_scale(r, r, 1.0f / (length*length));
        }
    } break;
    case ENVLIGHT:
    {
        vec3_scale(r, ((EnvLight*)light_ptr)->color, ((EnvLight*)light_ptr)->intensity);
    } break;
    case AREALIGHT:
    {
        AreaLight *area_light = (AreaLight*)light_ptr;
        vec3_scale(r, area_light->color, area_light->intensity);
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

void MeshLight_init(MeshLight* mesh_light, Mesh* mesh)
{
    mesh_light->num_triangles = 0;
    mesh_light->surface_area = 0.0f;
    printf("num pre allocated %d\n", mesh->num_indices / 3);
    mesh_light->cdf = (float*)malloc(mesh->num_indices / 3 * sizeof(float));
    mesh_light->triangles = (void**)malloc(mesh->num_indices / 3 * sizeof(void*));
}

int MeshLight_addTriangle(MeshLight* mesh_light, Object_t obj)
{
    vec3 v0, v1, v2;
    if(obj.type != mesh_light->obj_type)
    {
        fprintf(stderr, "Geometry type don't match.\n");
        return 0;
    }
    if(obj.type == FLAT_TRIANGLE)
    {
        FlatTriangle* tri_ptr = (FlatTriangle*)(obj.ptr);
        vec3_copy(v0, tri_ptr->v0);
        vec3_copy(v1, tri_ptr->v1);
        vec3_copy(v2, tri_ptr->v2);
    }else if(obj.type == SMOOTH_TRIANGLE)
    {
        SmoothTriangle* tri_ptr = (SmoothTriangle*)(obj.ptr);
        vec3_copy(v0, tri_ptr->v0);
        vec3_copy(v1, tri_ptr->v1);
        vec3_copy(v2, tri_ptr->v2);        
    }else
    {
        fprintf(stderr, "Wrong geometry type.\n");
        return 0;
    }
    vec3 e0, e1;
    vec3_sub(e0, v1, v0);
    vec3_sub(e1, v2, v0);
    vec3 cross_product;
    vec3_cross(cross_product, e0, e1);
    float tri_area = vec3_length(cross_product) * 0.5f;
    mesh_light->cdf[mesh_light->num_triangles] = tri_area + mesh_light->surface_area;
    mesh_light->triangles[mesh_light->num_triangles] = obj.ptr;
    mesh_light->surface_area += tri_area;
    mesh_light->num_triangles++;
    return 1;
}

void MeshLight_normalizeCDF(MeshLight* mesh_light)
{
    float norm_factor = 1.0f / mesh_light->cdf[mesh_light->num_triangles - 1];
    for(int i = 0; i < mesh_light->num_triangles; i++)
    {
        mesh_light->cdf[i] *= norm_factor;
    }
}

void MeshLight_destroy(MeshLight* mesh_light)
{
    if(mesh_light->cdf)
    {
        free(mesh_light->cdf);
    }
    if(mesh_light->triangles)
    {
        free(mesh_light->triangles);
    }
    mesh_light->num_triangles = 0;
    mesh_light->surface_area = 0.0f;
}

void MeshLight_genSample(vec3 sample, vec3 normal, MeshLight* mesh_light)
{
    float rand_float = (float)rand() / (float)RAND_MAX;

    int i = (mesh_light->num_triangles-1) / 2;
    int upper = mesh_light->num_triangles - 1;
    int lower = 0;
    while(1)
    {
        if(rand_float > mesh_light->cdf[i])
        {
            if(i+1 == mesh_light->num_triangles-1)
            {
                i++;
                break;
            }
            lower = i;
            i = (i + upper) / 2;
        }else
        {
            if(i > 0)
            {
                if(mesh_light->cdf[i-1] < rand_float)
                {
                    break;
                }
            }else
            {
                break;
            }
            upper = i;
            i = (i + lower) / 2;
        }
    }
   
    // p = (1 - sqrt(r1))v0 + (sqrt(r1)(1 - sqrt(r2))v1 + (r2*sqrt(r1))v2
    float r1 = (float)rand() / (float)RAND_MAX;
    float r2 = (float)rand() / (float)RAND_MAX;
    float sqrt_r1 = sqrtf(r1);
    vec3 a, b, c;
    vec3 v0, v1, v2;
    if(mesh_light->obj_type == FLAT_TRIANGLE)
    {
        FlatTriangle* tri_ptr = (FlatTriangle*)(mesh_light->triangles[i]);
        vec3_scale(a, tri_ptr->v0, 1.0 - sqrt_r1);
        vec3_scale(b, tri_ptr->v1, sqrt_r1 * (1.0f - r2));
        vec3_scale(c, tri_ptr->v2, r2 * sqrt_r1);
        vec3 a_plus_b;
        vec3_add(a_plus_b, a, b);
        vec3_add(sample, a_plus_b, c);
        vec3_copy(normal, tri_ptr->normal);
        vec3_copy(v0, tri_ptr->v0);
        vec3_copy(v1, tri_ptr->v1);
        vec3_copy(v2, tri_ptr->v2);
    }else if(mesh_light->obj_type == SMOOTH_TRIANGLE)
    {
        // TODO add interploated normal
        SmoothTriangle* tri_ptr = (SmoothTriangle*)(mesh_light->triangles[i]);
        vec3 e0, e1;
        vec3_sub(e0, tri_ptr->v1, tri_ptr->v0);
        vec3_sub(e1, tri_ptr->v2, tri_ptr->v0);
        vec3_cross(normal, e0, e1);
        vec3_normalize(normal, normal);

        vec3_scale(a, tri_ptr->v0, 1.0 - sqrt_r1);
        vec3_scale(b, tri_ptr->v1, sqrt_r1 * (1.0f - r2));
        vec3_scale(c, tri_ptr->v2, r2 * sqrt_r1);
        vec3 a_plus_b;
        vec3_add(a_plus_b, a, b);
        vec3_add(sample, a_plus_b, c);
        vec3_copy(v0, tri_ptr->v0);
        vec3_copy(v1, tri_ptr->v1);
        vec3_copy(v2, tri_ptr->v2);
    }
    /*
    AABB aabb;
    vec3_copy(aabb.min, v0);
    vec3_copy(aabb.max, v0);
    AABB_coverPoint(&aabb, v1);
    AABB_coverPoint(&aabb, v2);
    vec3 offset = {K_EPSILON, K_EPSILON, K_EPSILON};
    vec3_sub(aabb.min, aabb.min, offset);
    vec3_add(aabb.max, aabb.max, offset);

    if(!isInsideAABB(&aabb, sample))
    {
        printf("not inside aabb\n");
        printVec3WithText("sample", sample);
        printVec3WithText("v0", v0);
        printVec3WithText("v1", v1);    
        printVec3WithText("v2", v2);
    }
    */
}
void getEnvLightIncRadiance(vec3 r, const vec3 dir, EnvLight* env_light)
{
    if(env_light->type == CONSTANT)
    {
        vec3_scale(r, env_light->color, env_light->intensity);
    }else if(env_light->type == TEXTURE)
    {
        vec2 spherical;
        cartesianToSpherical(spherical, dir);
        vec2 uv;
        sphericalToUV(uv, spherical);
        vec3 tex_color;
        getTexColor(tex_color, &(env_light->env_map), uv);
        vec3_scale(r, tex_color, env_light->intensity);
    }
}

void EnvLight_destroy(EnvLight* env_light)
{
    if(env_light)
        env_light->intensity = 0.0f;
}

float AreaLight_sample_Li(vec3 Li, vec3 wi, float* t,
                          const ShadeRec *sr, const vec2 sample, const AreaLight* area_light)
{
    vec3 sample_point;
    if(area_light->obj_type == RECTANGLE)
    {
        vec3 displacement;
        Rectangle* rect = (Rectangle*)(area_light->obj_ptr);
        vec3_scale(displacement, rect->width, sample[0]);
        vec3_add(sample_point, rect->point, displacement);
        vec3_scale(displacement, rect->height, sample[1]);
        vec3_add(sample_point, sample_point, displacement);
    }else if(area_light->obj_type == SPHERE)
    {
        vec3 h_sample;
        mapSampleToHemisphere(h_sample, sample);
        getVec3InLocalBasis(sample_point, h_sample, sr->normal);
        vec3_add(sample_point, sample_point,
                 ((Sphere*)(area_light->obj_ptr))->center);
    }
    vec3 sample_to_hit_point;    
    vec3_sub(sample_to_hit_point, sample_point, sr->hit_point);
    vec3_normalize(wi, sample_to_hit_point);

    vec3_scale(Li, area_light->color, area_light->intensity);
    *t = vec3_length(sample_to_hit_point);
    return area_light->pdf;
}

float EnvLight_sample_Li(vec3 Li, vec3 wi, float* t, const ShadeRec* sr, const vec2 sample, const EnvLight* env_light)
{
    vec3 h_sample;
    mapSampleToHemisphere(h_sample, sample);
    getVec3InLocalBasis(wi, h_sample, sr->normal);
    vec3_scale(Li, env_light->color, env_light->intensity);
    *t = TMAX;
    return vec3_dot(wi, sr->normal) * INV_PI;
}

float Light_sample_Li(vec3 Li, vec3 wi, float* t,
                      const ShadeRec *sr, const vec2 sample, const void* light, const LightType type)
{
    // return light_pdf
    // Store incident radiance in Li
    switch(type)
    {
        case AREALIGHT:
        {
            AreaLight* area_light = (AreaLight*)light;
            AreaLight_sample_Li(Li, wi, t, sr, sample, area_light);
        } break;
        case ENVLIGHT:
        {
            EnvLight* env_light = (EnvLight*)light;
            EnvLight_sample_Li(Li, wi, t, sr, sample, env_light);
        } break;
    }
}
