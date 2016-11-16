#pragma once
#include "trace.h"
//float pathTrace(vec3, int, const Ray, TraceArgs trace_args);

traceFunc getTraceMediumFunc(const TraceType trace_type)
{
    traceFunc func;
    switch(trace_type)
    {
    case RAYCAST:
        func = &raycastMedium;
        break;
    case WHITTED:
        func = &whittedTraceMedium;
        break;
        /*
    case PATHTRACE:
        func = &pathTrace;
        break;
        */
    default:
        func = &raycastMedium;
    }
    return func;
}

// Returns the material pointer of the smallest partcipating medium that encloses the start point
Material* getMediumMatPtr(const vec3 start, const SceneObjects *so)
{
    // Use a stack to find the enclosing medium
    Material *stack[10];
    stack[0] = NULL;
    int top = -1;
    Ray ray;
    vec3_copy(ray.origin, start);
    vec3_assign(ray.direction, 0.0f, 0.0f, -1.0f);
    ShadeRec sr;
    float t = intersectTest(&sr, so, ray);
    int count = 0;
    while(t < TMAX)
    {
        if(sr.mat->mat_type == PARTICIPATING)
        {
            if(top == 9)
            {
                fprintf(stderr, "Too many enclosing media.\n");
                break;
            }else if(top = -1)
            {
                top++;
                stack[top] = sr.mat;
            }else
            {
                if(stack[top] == sr.mat)
                {
                    stack[top] = NULL;
                    top--;
                }else
                {
                    top++;
                    stack[top] = sr.mat;
                }
            }
        }
        getPointOnRay(ray.origin, ray, t);
        t = intersectTest(&sr, so, ray);
    }
    return stack[0];
}

float schlickPhaseFunc(const float cos_theta, const float k)
{
    float numerator = 1.0f - k * k;
    float a = 1.0f + k * cos_theta;
    float denom = 4.0f * (float)PI * (a*a);
    return numerator / denom;
}

void directIllumInScatter(vec3 radiance, const vec3 point, const float extinct_coeff, const float scatter_coeff,
                          const vec3 view_dir, TraceArgs trace_args)
{
    const SceneObjects *so = trace_args.objects;
    const SceneLights *sl = trace_args.lights;
    const int sample_index = trace_args.sample_index;

    ShadeRec tmp_sr;
    vec3_copy(tmp_sr.hit_point, point);
    Ray ray;
    vec3_copy(ray.origin, point);
    for(int j = 0; j < sl->num_lights; j++)
    {
        // construct a ray towards a light point
        getLightDir(ray.direction, sl->light_types[j], sl->light_ptrs[j], &tmp_sr, sample_index);
        if(shadowTest(j, sl, so, ray.direction, &tmp_sr)){continue;}
        ShadeRec light_sr;
        float t_light_exit = intersectTest(&light_sr, so, ray);
        float light_dist = calcLightDistance(sl->light_types[j], sl->light_ptrs[j], tmp_sr.hit_point);
        float medium_dist = min(t_light_exit, light_dist);
        vec3 light_rad = {0.0f, 0.0f, 0.0f};
        getIncRadiance(light_rad, sl->light_types[j], sl->light_ptrs[j], light_sr.hit_point);
        float rad_dec = powf((float)K_E, -extinct_coeff * medium_dist);
        vec3_scale(light_rad, light_rad, rad_dec);
        float phase = schlickPhaseFunc(vec3_dot(ray.direction, view_dir), 0.5f);
        vec3_scale(light_rad, light_rad, phase * scatter_coeff);
        vec3_add(radiance, radiance, light_rad);
    }
}

void raymarch(vec3 radiance, const Ray ray, TraceArgs trace_args)
{
    const SceneObjects *so = trace_args.objects;
    const SceneLights *sl = trace_args.lights;
    const int sample_index = trace_args.sample_index;
    vec3 h_sample;
    vec3_copy(h_sample, trace_args.h_sample);

    vec3 view_dir;
    vec3_negate(view_dir, ray.direction);
    const float extinct_coeff = 0.01f;
    const float phase_func = 1.0f / (float)(4.0 * PI);
    const float scatter_coeff = extinct_coeff * 0.5f;
    float t_seg = 5.0f;
    ShadeRec in_sr;
    // 1. FInd exiting t value and location
    float t_exit = intersectTest(&in_sr, so, ray);
    if(t_exit == TMAX){return;}
    if(in_sr.mat->mat_type == PARTICIPATING)
    {
        // 2. Calculate initial radiance
        Ray exit_ray = ray;
        vec3_copy(exit_ray.origin, in_sr.hit_point);
        vec3 init_rad = {0.0f, 0.0f, 0.0f};
        raycast(init_rad, 0, exit_ray, trace_args);
        // 3. Ray march back to front
        for(float i = t_exit; i > K_EPSILON; i -= t_seg)
        {
            // 4. At each segment, calculate direction illumination
            vec3 point;
            getPointOnRay(point, ray, i);
            vec3 total_light_rad = {0.0f, 0.0f, 0.0f};
            directIllumInScatter(total_light_rad, point, extinct_coeff, scatter_coeff, view_dir, trace_args);
            // Decrease initial radiance
            float rad_dec = powf((float)K_E, -extinct_coeff * t_seg);
            vec3_scale(init_rad, init_rad, rad_dec);
            vec3_add(init_rad, init_rad, total_light_rad);
        }
        vec3_copy(radiance, init_rad);
    }else
    {
        vec3 init_rad = {0.0f, 0.0f, 0.0f};
        if(in_sr.mat->mat_type == EMISSIVE)
        {
            vec3_scale(init_rad, in_sr.mat->ce, in_sr.mat->ke/1.0f);
        }else
        {
            for(int i = 0; i < sl->num_lights; i++)
            {
                vec3 light_dir;
                getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], &in_sr, sample_index);
                float ndotwi = vec3_dot(light_dir, in_sr.normal);
                if(ndotwi <= 0){continue;}
                if(shadowTest(i, sl, so, light_dir, &in_sr)){continue;}
                vec3 light_rad;
                getIncRadiance(light_rad, sl->light_types[i], sl->light_ptrs[i], in_sr.hit_point);
                float light_dist = calcLightDistance(sl->light_types[i], sl->light_ptrs[i], in_sr.hit_point);
                Ray light_ray;
                vec3_copy(light_ray.origin, in_sr.hit_point);
                vec3_copy(light_ray.direction, light_dir);
                ShadeRec light_sr;
                float light_t = intersectTest(&light_sr, so, light_ray);
                float dist_in_medium = 0.0f;
                if(light_t > light_dist)
                {
                    // Case where the light inside the medium or the light isn't physical
                    dist_in_medium = light_dist;
                }else if(light_sr.mat->mat_type == PARTICIPATING)
                {
                    dist_in_medium = light_t;
                }
                float rad_dec = powf((float)K_E, -extinct_coeff * dist_in_medium);
                vec3 dir_illum = {0.0f, 0.0f, 0.0f};
                vec3_scale(dir_illum, light_rad, rad_dec);
                vec3_add(init_rad, init_rad, dir_illum);
            }
        }
        // TODO: sometimes t_exit == T_MAX
        for(float i = t_exit; i > K_EPSILON; i -= t_seg)
        {
            // 4. At each segment, calculate direction illumination
            vec3 point;
            getPointOnRay(point, ray, i);
            ShadeRec tmp_sr;
            vec3_copy(tmp_sr.hit_point, point);
            vec3 total_light_rad = {0.0f, 0.0f, 0.0f};
            directIllumInScatter(total_light_rad, point, extinct_coeff, scatter_coeff, view_dir, trace_args);
            // Decrease initial radiance
            float rad_dec = powf((float)K_E, -extinct_coeff * t_seg);
            vec3_scale(init_rad, init_rad, rad_dec);
            vec3_add(init_rad, init_rad, total_light_rad);
        }
        vec3_copy(radiance, init_rad);
    }
}

void mediumMarch(vec3 out_rad, const vec3 entry_rad, const Ray ray, const float t_seg,
                 const float medium_dist, const float extinct_coeff,
                 const float scatter_coeff, TraceArgs trace_args)
{
    const SceneObjects *so = trace_args.objects;
    const SceneLights *sl = trace_args.lights;
    const int sample_index = trace_args.sample_index;
    vec3 h_sample;
    vec3_copy(h_sample, trace_args.h_sample);
    vec3 init_rad;
    vec3_copy(init_rad, entry_rad);
    vec3 view_dir;
    vec3_negate(view_dir, ray.direction);
    for(float i = medium_dist; i > K_EPSILON; i -= t_seg)
    {
        // 4. At each segment, calculate direction illumination
        vec3 point;
        getPointOnRay(point, ray, i);
        vec3 total_light_rad = {0.0f, 0.0f, 0.0f};
        directIllumInScatter(total_light_rad, point, extinct_coeff, scatter_coeff, view_dir, trace_args);
        // Decrease initial radiance
        float rad_dec = powf((float)K_E, -extinct_coeff * t_seg);
        vec3_scale(init_rad, init_rad, rad_dec);
        vec3_scale(total_light_rad, total_light_rad, (rad_dec * t_seg) + (1.0 - rad_dec)*0.5f*t_seg);
        vec3_add(init_rad, init_rad, total_light_rad);
    }
    vec3_copy(out_rad, init_rad);
}

void calcDirectIllumSurfaceInMedium(vec3 radiance, const ShadeRec *sr,
                                    const float extinct_coeff, const float scatter_coeff,
                                    TraceArgs trace_args)
{
    const SceneObjects *so = trace_args.objects;
    const SceneLights *sl = trace_args.lights;
    const int sample_index = trace_args.sample_index;
    vec3 h_sample;
    vec3_copy(h_sample, trace_args.h_sample);
    for(int i = 0; i < sl->num_lights; i++)
    {
        vec3 light_dir;
        getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], sr, sample_index);
        float ndotwi = vec3_dot(light_dir, sr->normal);
        if(ndotwi <= 0){continue;}
        if(shadowTest(i, sl, so, light_dir, sr)){continue;}
        vec3 light_rad;
        getIncRadiance(light_rad, sl->light_types[i], sl->light_ptrs[i], sr->hit_point);
        float light_dist = calcLightDistance(sl->light_types[i], sl->light_ptrs[i], sr->hit_point);
        Ray light_ray;
        vec3_copy(light_ray.origin, sr->hit_point);
        vec3_copy(light_ray.direction, light_dir);
        ShadeRec light_sr;
        float light_t = intersectTest(&light_sr, so, light_ray);
        float medium_dist = min(light_t, light_dist);
        float rad_dec = powf((float)K_E, -extinct_coeff * medium_dist);
        vec3 dir_illum = {0.0f, 0.0f, 0.0f};
        vec3_scale(dir_illum, light_rad, rad_dec);
        vec3_add(radiance, radiance, dir_illum);
    }
}

void fogmarch(vec3 radiance, const Ray ray, TraceArgs trace_args)
{
    const SceneObjects *so = trace_args.objects;
    const SceneLights *sl = trace_args.lights;
    const int sample_index = trace_args.sample_index;
    vec3 h_sample;
    vec3_copy(h_sample, trace_args.h_sample);
    vec3_copy(radiance, BLACK);
    vec3 view_dir;
    vec3_negate(view_dir, ray.direction);
    const float extinct_coeff = 0.002f;
    const float scatter_coeff = extinct_coeff * 0.6f;
    float t_seg = 5.0f;
    ShadeRec min_sr;
    float min_t = intersectTest(&min_sr, so, ray);

    if(min_t != TMAX)
    {
        vec3 init_rad = {0.0f, 0.0f, 0.0f};
        if(min_sr.mat->mat_type == EMISSIVE)
        {
            vec3_scale(init_rad, min_sr.mat->ce, min_sr.mat->ke/1.0f);
        }else
        {
            calcDirectIllumSurfaceInMedium(init_rad, &min_sr, extinct_coeff, scatter_coeff, trace_args);
        }

        // ray march back to camera
        mediumMarch(radiance, init_rad, ray, t_seg, min_t, extinct_coeff, scatter_coeff, trace_args);
    }else
    {
        const float max_fog_dist = 1000.0f;
        vec3 init_rad = {0.0f, 0.0f, 0.0f};
        mediumMarch(radiance, init_rad, ray, t_seg, max_fog_dist, extinct_coeff, scatter_coeff, trace_args);
    }
}

// Called when a ray enters or originates inside a medium
float raycastMedium(vec3 radiance, const int depth, const Ray ray, TraceArgs trace_args)
{
    const SceneObjects *so = trace_args.objects;
    const SceneLights *sl = trace_args.lights;
    const int sample_index = trace_args.sample_index;
    vec3 h_sample;
    vec3_copy(h_sample, trace_args.h_sample);
    vec3 view_dir;
    vec3_negate(view_dir, ray.direction);
    const float extinct_coeff = trace_args.medium_mat->extinct_coeff;
    const float scatter_coeff = trace_args.medium_mat->scatter_coeff;
    float t_seg = 5.0f;
    ShadeRec sr;
    float t = intersectTest(&sr, so, ray);
    if(t == TMAX){return t;} // Shouldn't happen if the ray is inside a medium
    vec3 init_rad = {0.0f, 0.0f, 0.0f};
    if(sr.mat->mat_type == PARTICIPATING) // Ray passes through the medium hitting nothing
    {
        TraceArgs new_trace_args = trace_args;
        new_trace_args.medium_mat = NULL;
        Ray new_ray;
        vec3_copy(new_ray.direction, ray.direction);
        getPointOnRay(new_ray.origin, ray, t);
        raycast(init_rad, depth-1, new_ray, new_trace_args);
    }else if(sr.mat->mat_type == EMISSIVE)
    {
        vec3_scale(init_rad, sr.mat->ce, sr.mat->ke/1.0f);
    }else
    {
        calcDirectIllumSurfaceInMedium(init_rad, &sr, extinct_coeff, scatter_coeff, trace_args);
    }
    mediumMarch(radiance, init_rad, ray, t_seg, t, extinct_coeff, scatter_coeff, trace_args);
    return t;
}

// Called when a ray enters or originates inside a medium
float whittedTraceMedium(vec3 radiance, const int depth, const Ray ray, TraceArgs trace_args)
{
    const SceneObjects *so = trace_args.objects;
    const SceneLights *sl = trace_args.lights;
    const int sample_index = trace_args.sample_index;
    vec3 h_sample;
    vec3_copy(h_sample, trace_args.h_sample);
    vec3 view_dir;
    vec3_negate(view_dir, ray.direction);
    const float extinct_coeff = 0.01f;
    const float scatter_coeff = extinct_coeff * 0.5f;
    float t_seg = 5.0f;
    ShadeRec sr;
    float t = intersectTest(&sr, so, ray);
    if(t == TMAX){return t;} // Shouldn't happen if the ray is inside a medium
    vec3 init_rad = {0.0f, 0.0f, 0.0f};
    if(sr.mat->mat_type == PARTICIPATING) // Ray passes through the medium hitting nothing
    {
        Ray new_ray;
        vec3_copy(new_ray.direction, ray.direction);
        getPointOnRay(new_ray.origin, ray, t);
        whittedTrace(init_rad, depth-1, new_ray, trace_args);
    }else if(sr.mat->mat_type == EMISSIVE)
    {
        vec3_scale(init_rad, sr.mat->ce, sr.mat->ke/1.0f);
    }else if(sr.mat->mat_type == REFLECTIVE)
    {
        Ray new_ray;
        vec3_copy(new_ray.origin, sr.hit_point);
        calcReflectRayDir(new_ray.direction, sr.normal, ray.direction);
        whittedTraceMedium(init_rad, depth+1, new_ray, trace_args);
    }else if(sr.mat->mat_type == TRANSPARENT)
    {
        float transmit_t = TMAX;
        float reflect_t = TMAX;
        vec3 reflected_illum = {0.0f, 0.0f, 0.0f};
        vec3 transmitted_illum = {0.0f, 0.0f, 0.0f};
        float ndotwo = vec3_dot(sr.normal, sr.wo);
        float kr = calcFresnelReflectance(&sr);
        float rand_float = (float)rand() / (float)RAND_MAX;
        if(rand_float <= kr) // Reflection
        {
            Ray new_ray;
            vec3_copy(new_ray.origin, sr.hit_point);
            calcReflectRayDir(new_ray.direction, sr.normal, ray.direction);
            reflect_t = whittedTraceMedium(init_rad, depth+1, new_ray, trace_args);
        }else // Transmission
        {
            vec3 transmit_dir = {0, 0, 0};
            float eta = calcTransmitDir(transmit_dir, &sr);
            float ndotwt = fabs(vec3_dot(sr.normal, transmit_dir));
            float kt = 1.0f - kr;

            Ray new_ray;
            vec3 btdf;
            vec3_scale(btdf, WHITE, kt / (eta*eta) / ndotwt);
            Ray transmitted_ray;
            vec3_copy(new_ray.origin, sr.hit_point);
            vec3_copy(new_ray.direction, transmit_dir);
            transmit_t = whittedTrace(transmitted_illum, depth+1, new_ray, trace_args);
            vec3_scale(transmitted_illum, transmitted_illum, ndotwt);
            vec3_mult(transmitted_illum, transmitted_illum, btdf);
        }
        vec3 color_filter_ref, color_filter_trans;
        if(ndotwo > 0.0f)
        {
            vec3_pow(color_filter_ref, sr.mat->cf_out, reflect_t);
            vec3_pow(color_filter_trans, sr.mat->cf_in, transmit_t);
        }else
        {
            vec3_pow(color_filter_ref, sr.mat->cf_in, reflect_t);
            vec3_pow(color_filter_trans, sr.mat->cf_out, transmit_t);
        }
        vec3_mult(reflected_illum, reflected_illum, color_filter_ref);
        vec3_mult(transmitted_illum, transmitted_illum, color_filter_trans);
        vec3_add(init_rad, init_rad, transmitted_illum);
        vec3_add(init_rad, init_rad, reflected_illum);
    }else
    {
        calcDirectIllumSurfaceInMedium(init_rad, &sr, extinct_coeff, scatter_coeff, trace_args);
    }
    mediumMarch(radiance, init_rad, ray, t_seg, t, extinct_coeff, scatter_coeff, trace_args);
    return t;
}

/*
  Duplicate function with different trace function calls
      Having the keep all verions updated whenever one of them is changed
      Too hacky
  Passing info in the radiance parameter or using the sign of the sample index or depth
      having one variable do two things and one of them being kinda obscure is risky
      Too hacky
  Global variable?
      globa states are bad 
      Too hacky
  So what's the good solution? Add a new variable in the paramters?
     There are already a lot of variables in the parameters
 */
