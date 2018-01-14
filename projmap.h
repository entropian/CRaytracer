#pragma once
#include "util/vec.h"
#include "scene/scenedata.h"
void addToBoundingSphere(vec3 center, float *radius, vec3 point)
{
    vec3 dist_vec;
    vec3_sub(dist_vec, point, center);
    float dist = vec3_length(dist_vec);
    if(dist > *radius)
    {
        float half_diff = (dist - *radius) / dist * 0.5f;
        vec3 center_mv;
        vec3_scale(center_mv, dist_vec, half_diff);
        vec3_add(center, center, center_mv);
        *radius += (dist - *radius) * 0.5f;
    }
}


int calcCausticBoundingSpheres(vec3 *centers, float *radii, const SceneObjects * const scene_objects)
{
    int caustic_obj_count = 0;
    Mesh *cur_mesh_ptr = NULL;
    vec3 mesh_sphere_center = {0.0f, 0.0f, 0.0f};
    float mesh_sphere_radius = 0;
    for(int i = scene_objects->num_non_grid_obj; i < scene_objects->num_obj; i++)
    {
        Object_t obj = scene_objects->objects[i];
        Material *mat = getObjectMatPtr(obj);        
        if(!mat)
        {
            fprintf(stderr, "Null material pointer.\n");
            continue;
        }
        if(mat->mat_type == MIRROR || mat->mat_type == TRANSPARENT)
        {
            if(obj.type != FLAT_TRIANGLE && obj.type != SMOOTH_TRIANGLE)
            {
                if(cur_mesh_ptr)
                {
                    vec3_copy(centers[caustic_obj_count], mesh_sphere_center);
                    radii[caustic_obj_count] = mesh_sphere_radius;                        
                    caustic_obj_count++;                    
                    cur_mesh_ptr = NULL;
                    vec3_assign(mesh_sphere_center, 0.0f, 0.0f, 0.0f);
                    mesh_sphere_radius = 0.0f;
                }
                vec3 center;
                float radius;
                if(calcBoundingSphere(center, &radius, obj))
                {
                    vec3_copy(centers[caustic_obj_count], center);
                    radii[caustic_obj_count] = radius;
                    caustic_obj_count++;
                }else
                {
                    fprintf(stderr, "Cannot calculate bounding sphere.\n");
                }
            }else
            {                
                if(obj.type == FLAT_TRIANGLE)
                {
                    FlatTriangle *triangle = (FlatTriangle*)(obj.ptr);
                    if(triangle->mesh_ptr != cur_mesh_ptr && cur_mesh_ptr)
                    {
                        vec3_copy(centers[caustic_obj_count], mesh_sphere_center);
                        radii[caustic_obj_count] = mesh_sphere_radius;                        
                        caustic_obj_count++;
                        cur_mesh_ptr = triangle->mesh_ptr;
                        calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);
                    }else if(!cur_mesh_ptr)
                    {
                        cur_mesh_ptr = triangle->mesh_ptr;                        
                        calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);                        
                    }else
                    {                        
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v0);
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v1);
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v2);                        
                    }
                }else if(obj.type == SMOOTH_TRIANGLE)
                {
                    SmoothTriangle *triangle = (SmoothTriangle*)(obj.ptr);
                    if(triangle->mesh_ptr != cur_mesh_ptr && cur_mesh_ptr)
                    {
                        vec3_copy(centers[caustic_obj_count], mesh_sphere_center);
                        radii[caustic_obj_count] = mesh_sphere_radius;                        
                        caustic_obj_count++;
                        cur_mesh_ptr = triangle->mesh_ptr;
                        calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);
                    }else if(!cur_mesh_ptr)
                    {
                        cur_mesh_ptr = triangle->mesh_ptr;
                        calcBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, obj);
                    }else
                    {
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v0);
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v1);
                        addToBoundingSphere(mesh_sphere_center, &mesh_sphere_radius, triangle->v2);
                    }
                }
            }
        }
    }
    if(cur_mesh_ptr)
    {
        vec3_copy(centers[caustic_obj_count], mesh_sphere_center);
        radii[caustic_obj_count] = mesh_sphere_radius;                        
        caustic_obj_count++;                    
    }
    return caustic_obj_count;
}

void projSphereToMap(char *map, const int theta_row, const int phi_column,
                     const vec3 center, const float radius, const vec3 light_pos)
{
    vec3 sphere_to_light;
    vec3_sub(sphere_to_light, center, light_pos);
    vec3 center_dir;
    float dist;
    vec3_normalize(center_dir, sphere_to_light);
    dist = vec3_length(sphere_to_light);
    vec3 perp_vec;
    vec3_cross(perp_vec, center_dir, JITTERED_UP);
    vec3_normalize(perp_vec, perp_vec);
    vec3_scale(perp_vec, perp_vec, radius);
    vec3 side_vec;
    vec3_add(side_vec, sphere_to_light, perp_vec);
    vec3_normalize(side_vec, side_vec);
    float angle = acos(vec3_dot(center_dir, side_vec));

    // spherical coordinates
    // atan2 output range is -PI, PI 
    float theta_per_cell = (float)(PI / (double)theta_row);
    float phi_per_cell = (float)(2.0 * PI / (double)phi_column);
    float theta = PI * 0.5 - asin(center_dir[1]);
    float phi = atan2(center_dir[2], center_dir[0]);
    int theta_index = theta / theta_per_cell;
    int phi_index = (phi + PI) / phi_per_cell;
    int theta_half_span = angle / theta_per_cell + 1;

    for(int i = theta_index - theta_half_span; i < theta_index + theta_half_span; i++)
    {
        int theta_index = i * phi_column;
        int phi_offset = 0;
        if(i < 0)
        {
            theta_index = (abs(i) - 1) * phi_column;
            phi_offset = phi_column / 2;
        }
        if(i >= theta_row)
        {
            theta_index = (theta_row - (i - theta_row) - 1) * phi_column;
            phi_offset = phi_column / 2;
        }
        float ang_from_equator = abs(i - theta_row / 2) * theta_per_cell;
        float cur_phi_per_cell = (float)cos(ang_from_equator + (theta_per_cell * 0.5f)) * phi_per_cell;
        cur_phi_per_cell = fabs(cur_phi_per_cell);
        int phi_half_span = angle / cur_phi_per_cell + 1;
        if(phi_half_span * 2 >= phi_column)
        {
            for(int j = 0; j < phi_column; j++)
            {
                map[theta_index + j] = 1;
            }
            continue;
        }
        for(int j = phi_index - phi_half_span; j < phi_index + phi_half_span; j++)
        {
            int phi_index = j;
            if(j < 0)
            {
                phi_index = phi_column + j;
            }
            if(j >= phi_column)
            {
                phi_index = j - phi_column;
            }
            phi_index = (phi_index + phi_offset) % phi_column;
            int index = theta_index + phi_index;
            assert(index > -1 && index < theta_row * phi_column);
            map[index] = 1;
        }
    } 
}

void buildProjMap(SceneLights *sl, const vec3 *centers, const float *radii, const int num_bspheres)
{
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == POINTLIGHT)
        {
            PointLight *point_light = (PointLight*)(sl->light_ptrs[i]);
            point_light->proj_map = (char*)calloc(THETA_ROW * PHI_COLUMN, sizeof(char));

            for(int j = 0; j < num_bspheres; j++)
            {
                projSphereToMap(point_light->proj_map, THETA_ROW, PHI_COLUMN, centers[j], radii[j], point_light->point);
            }
            
            float theta_per_cell = (float)(PI / (double)THETA_ROW);
            float phi_per_cell = (float)(2.0 * PI / (double)PHI_COLUMN);
            float base_area = theta_per_cell * phi_per_cell;
            float sum = 0;
            for(int i = 0; i < THETA_ROW; i++)
            {
                float angle_from_equator = (float)abs(i - THETA_ROW * 0.5f) * theta_per_cell;                
                float cell_area = base_area * cos(angle_from_equator);
                for(int j = 0; j < PHI_COLUMN; j++)
                {
                    if(point_light->proj_map[i*PHI_COLUMN + j])
                    {
                        sum += cell_area;
                    }
                }
            }
            point_light->proj_coverage = sum;
        }
    }
}
