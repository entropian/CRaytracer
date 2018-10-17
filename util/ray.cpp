#include "ray.h"
#include "util.h"

void getPointOnRay(vec3 r, const Ray ray, const float t)
{
    vec3 displacement;
    vec3_scale(displacement, ray.direction, t);
    vec3_add(r, ray.origin, displacement);
}

void calcReflectRayDir(vec3 reflect_dir, const vec3 normal, const vec3 incident_dir)
{
    float magnitude = -2.0f * vec3_dot(incident_dir, normal);
    vec3 adjusted_normal;
    vec3_scale(adjusted_normal, normal, magnitude);
    vec3_add(reflect_dir, incident_dir, adjusted_normal);
}

float calcTransmitDir(vec3 transmit_dir, const vec3 normal, const vec3 wo, const float ior_in,
                      const float ior_out)
{
    vec3 n;
    vec3_copy(n, normal);
    float cos_theta_i = vec3_dot(n, wo);
    float eta = ior_in / ior_out;
    if(cos_theta_i < 0.0f)
    {
        cos_theta_i = -cos_theta_i;
        vec3_negate(n, n);
        eta = 1.0f / eta;
    }

    float tmp = 1.0f - (1.0f - cos_theta_i * cos_theta_i) / (eta * eta);
    float cos_theta2 = (float)sqrt(tmp);
    vec3 tmp_vec3_1, tmp_vec3_2;
    vec3_scale(tmp_vec3_1, wo, -1.0f / eta);
    vec3_scale(tmp_vec3_2, n, cos_theta2 - cos_theta_i / eta);
    vec3_sub(transmit_dir, tmp_vec3_1, tmp_vec3_2);
    return eta;
}

bool refract(vec3 wi, const vec3 wo, const vec3 n, const float eta)
{
    float cos_theta_i = vec3_dot(n, wo);
    float sin_2_theta_i = max(0.0f, 1.0f - cos_theta_i * cos_theta_i);
    float sin_2_theta_t = eta * eta * sin_2_theta_i;
    if(sin_2_theta_t >= 1.0f) return false;
    float cos_theta_t = sqrtf(1.0f - sin_2_theta_t);
    
    vec3 perp_comp, neg_wo;
    vec3_negate(neg_wo, wo);
    vec3_scale(perp_comp, neg_wo, eta);

    float scale_factor = eta * cos_theta_i - cos_theta_t;
    vec3 parallel_comp;
    vec3_scale(parallel_comp, n, scale_factor);
    vec3_add(wi, perp_comp, parallel_comp);
    return true;
}

void resetRay(Ray* ray, const vec3 origin, const vec3 dir)
{
    vec3_copy(ray->origin, origin);
    vec3_copy(ray->direction, dir);
}
