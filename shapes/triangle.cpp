#include "triangle.h"
#include "../util/simd.h"

void calcTriangleNormal(vec3 r, const vec3 v0, const vec3 v1, const vec3 v2)
{
    // Assume the vertices are ordered counterclock wise
    vec3 tmp1, tmp2, tmp3;
    vec3_sub(tmp1, v1, v0);
    vec3_sub(tmp2, v2, v0);
    vec3_cross(tmp3, tmp1, tmp2);
    vec3_normalize(r, tmp3);    
}

float calcTriangleIntersect(float* beta_out, float* gamma_out,
                            const vec3 v0, const vec3 v1, const vec3 v2, const Ray ray)
{
    //printf("in regular\n");
    // o + td = v0 + v(v1-v0) + w(v2-v0)
    // v(v0-v1) + w(v0-v2) + td = v0 - o

    float a = v0[0] - v1[0], b = v0[0] - v2[0], c = ray.direction[0], d = v0[0] - ray.origin[0];      
    float e = v0[1] - v1[1], f = v0[1] - v2[1], g = ray.direction[1], h = v0[1] - ray.origin[1];
    float i = v0[2] - v1[2], j = v0[2] - v2[2], k = ray.direction[2], l = v0[2] - ray.origin[2];
    /*
    printf("a %f, b %f, c %f, d %f\n", a, b, c, d);
    printf("e %f, f %f, g %f, h %f\n", e, f, g, h);
    printf("i %f, j %f, k %f, l %f\n", i, j, k, l);
    */
    float m = f*k - g*j, n = h*k - g*l, p = f*l - h*j;
    float q = g*i - e*k, s = e*j - f*i;
    //printf("m %f, n %f, p %f, q %f, s %f\n", m, n, p, q, s);

    float inv_denom = 1.0f / (a*m + b*q + c*s);
    //printf("inv_denom %f\n", inv_denom);

    float e1 = d*m - b*n - c*p;
    float beta = e1 * inv_denom;
    //printf("e1 %f, beta %f\n", e1, beta);

    float r = e*l - h*i;
    float e2 = a*n + d*q + c*r;
    float gamma = e2 * inv_denom;
    //printf("r %f, e2 %f, gamma %f\n", r, e2, gamma);
    *beta_out = beta;
    *gamma_out = gamma;

    if(beta < 0.0f)
    {
        return TMAX;
    }    

    if(gamma < 0.0f)
    {
        return TMAX;
    }

    if(beta + gamma > 1.0f)
    {
        return TMAX;
    }

    /*
    if(beta < 0.0f || gamma < 0.0f || beta + gamma > 1.0f) 
    {
        return TMAX;
    }
    */
    
    float e3 = a*p - b*r + d*s;
    float t = e3 * inv_denom;
    //printf("e3 %f, t %f\n", e3, t);
    if(t < K_EPSILON)
    {
        return TMAX;
    }
    return t;
}


// Intersects a ray with four triangles
__m128 calcTriangleIntersect4(__m128* beta_out, __m128* gamma_out, const vec3_4* v0, const vec3_4* v1,
                                           const vec3_4* v2, const vec3_4* ray_o, const vec3_4* ray_d)
{
    __m128 a = _mm_sub_ps(v0->x, v1->x), b = _mm_sub_ps(v0->x, v2->x), c = ray_d->x, d = _mm_sub_ps(v0->x, ray_o->x);
    __m128 e = _mm_sub_ps(v0->y, v1->y), f = _mm_sub_ps(v0->y, v2->y), g = ray_d->y, h = _mm_sub_ps(v0->y, ray_o->y);
    __m128 i = _mm_sub_ps(v0->z, v1->z), j = _mm_sub_ps(v0->z, v2->z), k = ray_d->z, l = _mm_sub_ps(v0->z, ray_o->z);
    __m128* ptr = &a;
    //m = f*k - g*j, n = h*k - g*l, p = f*l - h*j;
    //q = g*i - e*k, s = e*j - f*i;
    __m128 m = _mm_sub_ps(_mm_mul_ps(f, k), _mm_mul_ps(g, j));
    __m128 n = _mm_sub_ps(_mm_mul_ps(h, k), _mm_mul_ps(g, l));
    __m128 p = _mm_sub_ps(_mm_mul_ps(f, l), _mm_mul_ps(h, j));
    __m128 q = _mm_sub_ps(_mm_mul_ps(g, i), _mm_mul_ps(e, k));
    __m128 s = _mm_sub_ps(_mm_mul_ps(e, j), _mm_mul_ps(f, i));
    
    // inv_denom = 1.0f / (a*m + b*q + c*s);
    __m128 one = _mm_set1_ps(1.0f);
    __m128 tmp1 = _mm_mul_ps(a, m);
    __m128 tmp2 = _mm_mul_ps(b, q);
    __m128 tmp3 = _mm_mul_ps(c, s);
    __m128 inv_denom = _mm_div_ps(one, _mm_add_ps(_mm_add_ps(tmp1, tmp2), tmp3));

    // e1 = d*m - b*n - c*p;
    tmp1 = _mm_mul_ps(d, m);
    tmp2 = _mm_mul_ps(b, n);
    tmp3 = _mm_mul_ps(c, p);
    __m128 e1 = _mm_sub_ps(_mm_sub_ps(tmp1, tmp2), tmp3);

    __m128 beta = _mm_mul_ps(e1, inv_denom);

    // r = e*l - h*i;      
    tmp1 = _mm_mul_ps(e, l);
    tmp2 = _mm_mul_ps(h, i);
    __m128 r = _mm_sub_ps(tmp1, tmp2);
    
    // e2 = a*n + d*q + c*r;
    tmp1 = _mm_mul_ps(a, n);
    tmp2 = _mm_mul_ps(d, q);
    tmp3 = _mm_mul_ps(c, r);
    __m128 e2 = _mm_add_ps(_mm_add_ps(tmp1, tmp2), tmp3);

    __m128 gamma = _mm_mul_ps(e2, inv_denom);
    *beta_out = beta;
    *gamma_out = gamma;

    __m128 zeroes = _mm_set1_ps(0.0f);
    __m128 beta_mask = _mm_cmplt_ps(zeroes, beta);
    __m128 gamma_mask = _mm_cmplt_ps(zeroes, gamma);
    __m128 sum_mask = _mm_cmplt_ps(_mm_add_ps(beta, gamma), one);
    __m128 mask = _mm_and_ps(_mm_and_ps(beta_mask, gamma_mask), sum_mask);
    unsigned int max_int = 0xffffffff;
    float *float_ptr = (float*)(&max_int);
    __m128 max_float = _mm_set1_ps(*float_ptr);    

    //float e3 = a*p - b*r + d*s;
    tmp1 = _mm_mul_ps(a, p);
    tmp2 = _mm_mul_ps(b, r);
    tmp3 = _mm_mul_ps(d, s);    
    __m128 e3 = _mm_add_ps(_mm_sub_ps(tmp1, tmp2), tmp3);
    
    __m128 t = _mm_and_ps(_mm_mul_ps(e3, inv_denom), mask);
    __m128 k_ep = _mm_set1_ps(K_EPSILON);
    __m128 k_mask = _mm_cmplt_ps(k_ep, t);
    mask = _mm_and_ps(mask, k_mask);
    __m128 inv_mask = _mm_andnot_ps(mask, max_float);
    __m128 t_max = _mm_set1_ps(TMAX);
    __m128 t_max_mask = _mm_and_ps(inv_mask, t_max);
    t = _mm_and_ps(t, mask);
    t = _mm_or_ps(t, t_max_mask);    
    return t;
}


float rayIntersectTriangle(ShadeRec* sr, Triangle* tri, const Ray ray)
{
    float gamma, beta; // For smooth triangles. Unused here.
    float t = calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);

    vec3_copy(sr->normal, tri->normal);
    getPointOnRay(sr->hit_point, ray, t);
    vec3_negate(sr->wo, ray.direction);
    if(vec3_dot(sr->wo, sr->normal) < 0.0f)
    {
        vec3_negate(sr->normal, sr->normal);
    }
    sr->mat = *(tri->mat);
    return t;
}

void interpTexcoord(vec2 uv_out,const float beta, const float gamma,
                    const Mesh* mesh_ptr, const int i0, const int i1, const int i2)
{
    vec2 uv0 = {mesh_ptr->texcoords[i0*2], mesh_ptr->texcoords[i0*2 + 1]};
    vec2 uv1 = {mesh_ptr->texcoords[i1*2], mesh_ptr->texcoords[i1*2 + 1]};
    vec2 uv2 = {mesh_ptr->texcoords[i2*2], mesh_ptr->texcoords[i2*2 + 1]};

    vec2 tmp;
    vec2_scale(uv_out, uv0, 1.0f - beta - gamma);
    vec2_scale(tmp, uv1, beta);
    vec2_add(uv_out, uv_out, tmp);
    vec2_scale(tmp, uv2, gamma);
    vec2_add(uv_out, uv_out, tmp);
}

float rayIntersectFlatTriangle(ShadeRec* sr, FlatTriangle* tri, const Ray ray)
{
    float gamma, beta; // For smooth triangles. Unused here.    
    float t = calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);

    //if(tri->mesh_ptr->num_texcoords > 0 && tri->mat->tex_flags != NO_TEXTURE)
    if(tri->mesh_ptr->num_texcoords > 0)
    {
        interpTexcoord(sr->uv, beta, gamma, tri->mesh_ptr, tri->i0, tri->i1, tri->i2);
    }    
    vec3_copy(sr->normal, tri->normal);
    getPointOnRay(sr->hit_point, ray, t);
    vec3_negate(sr->wo, ray.direction);    
    sr->mat = *(tri->mat);
    return t;
}

void interpTriangleVec3(vec3 out, const float beta, const float gamma,
                        const vec3 v0, const vec3 v1, const vec3 v2)
{
    vec3 tmp;
    vec3_scale(out, v0, 1.0f - beta - gamma);
    vec3_scale(tmp, v1, beta);
    vec3_add(out, out, tmp);
    vec3_scale(tmp, v2, gamma);
    vec3_add(out, out, tmp);
    vec3_normalize(out, out);
}

float rayIntersectSmoothTriangle(ShadeRec* sr, SmoothTriangle* tri, const Ray ray)
{
    float beta, gamma;
    float t = calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);
    if(t == TMAX){return t;}

    Mesh* mesh = tri->mesh_ptr;

    vec3 n0, n1, n2, tmp_normal;
    int index_0 = tri->i0 * 3;
    int index_1 = tri->i1 * 3;
    int index_2 = tri->i2 * 3;
    vec3_assign(n0, mesh->normals[index_0], mesh->normals[index_0+1], mesh->normals[index_0+2]);
    vec3_assign(n1, mesh->normals[index_1], mesh->normals[index_1+1], mesh->normals[index_1+2]);
    vec3_assign(n2, mesh->normals[index_2], mesh->normals[index_2+1], mesh->normals[index_2+2]);
    interpTriangleVec3(tmp_normal, beta, gamma, n0, n1, n2);
    mat3_mult_vec3(sr->normal, *(tri->normal_mat), tmp_normal);
    vec3_normalize(sr->normal, sr->normal);

    interpTexcoord(sr->uv, beta, gamma, tri->mesh_ptr, tri->i0, tri->i1, tri->i2);
    vec3 tangent, binormal, tex_normal, normal, tmp;
    interpTriangleVec3(tmp, beta, gamma,
                       mesh->tangents[tri->i0], mesh->tangents[tri->i1], mesh->tangents[tri->i2]);
    mat3_mult_vec3(tangent, *(tri->normal_mat), tmp);
    vec3_normalize(sr->dpdu, tangent);

    vec3_cross(binormal, sr->normal, tangent);
    vec3_normalize(sr->dpdv, binormal);

    //assert(sr->normal[0] != 0.0f || sr->normal[1] != 0.0f || sr->normal[2] != 0.0f);
    getPointOnRay(sr->hit_point, ray, t);
    vec3_negate(sr->wo, ray.direction);
    sr->mat = *(tri->mat);
    return t;
}



float getSmoothTriangleShadeRec(ShadeRec* sr, SmoothTriangle* tri, const Ray ray,
                                const float beta, const float gamma, const float t)
{
    //float t = calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);
    //if(t == TMAX){return t;}

    Mesh* mesh = tri->mesh_ptr;

    vec3 n0, n1, n2, tmp_normal;
    int index_0 = tri->i0 * 3;
    int index_1 = tri->i1 * 3;
    int index_2 = tri->i2 * 3;
    vec3_assign(n0, mesh->normals[index_0], mesh->normals[index_0+1], mesh->normals[index_0+2]);
    vec3_assign(n1, mesh->normals[index_1], mesh->normals[index_1+1], mesh->normals[index_1+2]);
    vec3_assign(n2, mesh->normals[index_2], mesh->normals[index_2+1], mesh->normals[index_2+2]);
    interpTriangleVec3(tmp_normal, beta, gamma, n0, n1, n2);
    mat3_mult_vec3(sr->normal, *(tri->normal_mat), tmp_normal);
    vec3_normalize(sr->normal, sr->normal);
    vec3 surface_normal;
    vec3_copy(surface_normal, sr->normal);
    //assert(surface_normal[0] != 0.0f || surface_normal[1] != 0.0f || surface_normal[2] != 0.0f);
    
    //if(tri->mesh_ptr->num_texcoords > 0 && tri->mat->tex_flags != NO_TEXTURE)
    if(tri->mesh_ptr->num_texcoords > 0)
    {
        interpTexcoord(sr->uv, beta, gamma, tri->mesh_ptr, tri->i0, tri->i1, tri->i2);
        //if(tri->mat->tex_flags & NORMAL)
        if(Material_hasNormalMap(tri->mat))
        {
            vec3 tangent, binormal, tex_normal, normal, tmp;
            interpTriangleVec3(tmp, beta, gamma,
                               mesh->tangents[tri->i0], mesh->tangents[tri->i1], mesh->tangents[tri->i2]);
            mat3_mult_vec3(tangent, *(tri->normal_mat), tmp);
            vec3_normalize(tangent, tangent);

            vec3_cross(binormal, sr->normal, tangent);
            vec3_normalize(binormal, binormal);

            //getMaterialNormalTexColor(tex_normal, tri->mat, sr->uv);
            Material_getNormalMapValue(tex_normal, tri->mat, sr->uv);
            orthoNormalTransform(normal, tangent, binormal, sr->normal, tex_normal);
            vec3_normalize(normal, normal);
            vec3_copy(sr->normal, normal);
        }
    }

    if(vec3_equal(sr->normal, BLACK))
    {
        vec3_copy(sr->normal, surface_normal);
    }
    //assert(sr->normal[0] != 0.0f || sr->normal[1] != 0.0f || sr->normal[2] != 0.0f);
    getPointOnRay(sr->hit_point, ray, t);
    vec3_negate(sr->wo, ray.direction);
    sr->mat = *(tri->mat);
    return t;
}

