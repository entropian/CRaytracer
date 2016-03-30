/*
  Things that can be swapped:
  camera --  perspective vs orthographic
  objects -- implicit surfaces vs explicit surfaces
  sampling -- different sampling algorithms

  Features to consider:
  being able to scale the camera and viewplane 

  TODO:
  gamma correction
  consider how variables should be grouped into objects?
 */

    /*
      Camera related things: fov, resolution, position, orientation.
      Be cool and throw them into a class?
      What design pattern is this?
     */

#include <GL/glew.h>

#if __GNUG__
#    include <GLFW/glfw3.h>
#else
#    include <GL/glfw3.h>
#    include <SOIL.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "vec.h"
#include "shaders.h"
#include "shapes.h"
#include "glcode.h"
#include "sampling.h"
#include "camera.h"
#include "lights.h"
#include "constants.h"
#include "sceneobj.h"
#include "lighting.h"

vec3 cam_position = {0.0f, 0.0f, 0.0f};

// TODO: not finished
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(action == GLFW_PRESS)
    {
        switch(key)
        {
        case GLFW_KEY_A:
            cam_position[0] -= 0.01f;
            break;
        case GLFW_KEY_D:
            cam_position[0] += 0.01f;
            break;
        case GLFW_KEY_W:
            cam_position[2] -= 0.01f;
            break;
        case GLFW_KEY_S:
            cam_position[2] += 0.01f;
            break;            
        }
    }
}

void initSpheres(SceneObjects *so)
{
    // TODO: tidy up material assignment
    Sphere* sphere = (Sphere*)malloc(sizeof(Sphere));
    
    vec3_assign(sphere->center, 0.0f, 0.0f, -4.0f);
    sphere->radius = 1.0f;
    sphere->shadow = true;        
    vec3_copy(sphere->mat.cd, PINK);
    vec3_copy(sphere->mat.ca, PINK);
    vec3_copy(sphere->mat.cs, PINK);    
    sphere->mat.kd = 0.6f;
    sphere->mat.ka = 1.0f;
    sphere->mat.ks = 0.4f;
    sphere->mat.exp = 5.0f;        
    sphere->mat.mat_type = MATTE;
    sphere->mat.shadow = true;
    so->obj_ptrs[so->num_obj] = sphere; 
    so->obj_types[so->num_obj] = SPHERE;
    (so->num_obj)++;

    sphere = (Sphere*)malloc(sizeof(Sphere));
    vec3_assign(sphere->center, 2.0f, 0.0f, -8.0f);
    sphere->radius = 1.0f;
    sphere->shadow = true;            
    vec3_copy(sphere->mat.cd, CYAN);
    vec3_copy(sphere->mat.ca, CYAN);
    vec3_copy(sphere->mat.cs, CYAN);    
    sphere->mat.kd = 0.6f;
    sphere->mat.ka = 1.0f;
    sphere->mat.ks = 0.4f;
    sphere->mat.exp = 10.0f;            
    sphere->mat.mat_type = PHONG;
    sphere->mat.shadow = true;
    so->obj_ptrs[so->num_obj] = sphere; 
    so->obj_types[so->num_obj] = SPHERE;
    (so->num_obj)++;    
}

void initPlanes(SceneObjects *so)
{
    Plane* plane = (Plane*)malloc(sizeof(Plane));
    plane->shadow  = true;    
    vec3_assign(plane->point, 0.0f, -1.0f, 0.0f);
    vec3_copy(plane->normal, UP);
    vec3_copy(plane->mat.cd, GREY);
    vec3_copy(plane->mat.ca, GREY);
    plane->mat.kd = 0.6f;
    plane->mat.ka = 1.0f;
    plane->mat.mat_type  = MATTE;
    plane->mat.shadow = true;    
    so->obj_ptrs[so->num_obj] = plane; 
    so->obj_types[so->num_obj] = PLANE;
    (so->num_obj)++;    
}

void initRectangles(SceneObjects *so)
{
    Rectangle* rect = (Rectangle*)malloc(sizeof(Rectangle));
    rect->shadow = true;
    vec3_assign(rect->point, -2.0f, -1.0f, -6.0f);
    // Current task: investigate how changing the width and height vector affects the rectangle's shape
    vec3_assign(rect->width, 2.0f, 0.0f, 0.0f);
    vec3_assign(rect->height, 0.0f, 4.0f, 0.0f);    
    vec3_copy(rect->normal, BACKWARD);
    vec3_copy(rect->mat.cd, YELLOW);
    vec3_copy(rect->mat.ca, YELLOW);
    rect->mat.kd = 0.6f;
    rect->mat.ka = 1.0f;
    rect->mat.mat_type  = MATTE;
    rect->mat.shadow = true;
    so->obj_ptrs[so->num_obj] = rect; 
    so->obj_types[so->num_obj] = RECTANGLE;
    (so->num_obj)++;        
}

void initSceneObjects(SceneObjects *so)
{
    for(int i = 0; i < MAX_OBJECTS; i++)
    {
        so->obj_ptrs[i] = NULL;
    }
    so->num_obj = 0;
    initSpheres(so);
    initPlanes(so);
    //initRectangles(so);
}

void initSceneLights(SceneLights* sl)
{
    for(int i = 0; i < MAX_LIGHTS; i++)
    {
        sl->light_ptrs[i] = NULL;
    }
    sl->num_lights = 0;    
    DirLight* dir_light = (DirLight*)malloc(sizeof(DirLight));
    float intensity = 3.0f;
    vec3 direction = {10.0f, 10.0f, 10.0f};
    vec3_normalize(direction, direction);
    assignDirLight(dir_light, intensity, WHITE, direction);
    sl->shadow[sl->num_lights] = true;
    sl->light_ptrs[sl->num_lights] = dir_light;
    (sl->num_lights)++;

    PointLight* point_light = (PointLight*)malloc(sizeof(PointLight));
    intensity = 0.5f;
    vec3 point = {-3.0f, -5.0f, 0.0f};
    assignPointLight(point_light, intensity, WHITE, point);
    sl->shadow[sl->num_lights] = true;    
    sl->light_ptrs[sl->num_lights] = point_light;
    (sl->num_lights)++;    
}

int main()
{
    int window_width = 900, window_height = 900;
    GLFWwindow* window = initWindow(window_width, window_height);
    glfwSetKeyCallback(window, keyCallback);
    
    GlViewport viewport;
    initViewport(&viewport);
    
    unsigned char *image;
    int frame_res_width = 900, frame_res_height = 900;
    //int frame_res_width = 360, frame_res_height = 360;    
    int num_pixels = frame_res_width * frame_res_height;
    image = (unsigned char*)calloc(num_pixels * 3, sizeof(char));

    vec3 bg_color = {1.0f, 1.0f, 1.0f};

    // Scene Objects
    SceneObjects scene_objects;
    initSceneObjects(&scene_objects);

    // Ambient light
    vec3 amb_color = {1.0f, 1.0f, 1.0f};
    float amb_ls = 1.0f;
    bool amb_occlusion = true;

    // Scene Lights
    SceneLights scene_lights;
    initSceneLights(&scene_lights);

    // Samples
    srand(time(NULL));
    int num_samples = 256;
    int num_sets = 83;    
    Samples samples = Samples_default;
    //genRegularSamples(&samples, num_samples, num_sets);
    genMultijitteredSamples(&samples, num_samples, num_sets);
    //genHammersleySamples(&samples, num_samples, num_sets);
    mapSamplesToDisk(&samples);
    mapSamplesToHemisphere(&samples, 1);

    // Camera
    Camera camera;
    initPinholeCameraDefault(&camera);
    //initThinLensCameraDefault(&camera);
    //camera.focal_length = 4.0f;
    //camera.lens_radius = 0.2f;


    vec3 position = {0.0f, 0.5f, 0.0f};
    vec3 look_point = {0.0f, 0.0f, -4.0f};
    vec3 up_vec = {0.0f, 1.0f, 0.0f};
    cameraLookAt(&camera, position, look_point, up_vec);

    // TODO: throw these in a struct for camera scaling
    float fov = 90.0f;
    float frame_length = 2.0f * (sin(fov/2.0f) * camera.focal_pt_dist);
    float frame_height = frame_length * (float)(frame_res_height)/(float)(frame_res_width);    
    float pixel_length = frame_length/(float)(frame_res_width);

//#if 0
    time_t startTime, endTime;
    time(&startTime);        
    unsigned sample_index = 0;    
    for(int i = 0; i < num_pixels; i++)
    {
        vec3 color = {0.0f, 0.0f, 0.0f};
        for(int p = 0; p < num_samples; p++)
        {
            vec2 sample, disk_sample, imageplane_coord;
            getNextShuffledSample(sample, &samples);
            imageplane_coord[0] = -frame_length/2 + pixel_length * ((float)(i % frame_res_width) + sample[0]);            
            imageplane_coord[1] = frame_height/2 - pixel_length * ((float)(i / frame_res_width) + sample[1]);
            
            Ray ray;
            switch(camera.camera_type)
            {
            case Pinhole:
                calcRayPinhole(&ray, imageplane_coord, &camera);
                break;
            case ThinLens:
                getNextDiskSample(disk_sample, &samples);            
                calcRayThinLens(&ray, imageplane_coord, disk_sample, &camera);
                break;
            }

            float min_t = TMAX,  tmp_t = TMAX;
            ShadeRec min_sr;
            min_t = intersectTest(&min_sr, &scene_objects, ray);

            // Shading
            if(min_t < TMAX)
            {
                vec3 radiance = {0.0f, 0.0f, 0.0f};
                // Ambient component
                // ka*ca * amb_inc_radiance
                vec3 amb_inc_radiance;
                vec3_scale(amb_inc_radiance, amb_color, amb_ls);
                vec3 reflectance;
                vec3_scale(reflectance, min_sr.mat->ca, min_sr.mat->ka);

                // Ambient Occlusion
                if(amb_occlusion && AOTest(&samples, &scene_objects, &min_sr) < TMAX)
                {
                    vec3 min_amb;
                    vec3_scale(min_amb, amb_color, 0.01f);
                    vec3_copy(radiance, min_amb);
                }else
                {
                    vec3_mult(radiance, amb_inc_radiance, reflectance);
                }
                
                /*
                for(int i = 0; i < scene_lights.num_lights; i++)
                {
                    vec3 light_dir;
                    getLightDir(light_dir, scene_lights.light_types[i], scene_lights.light_ptrs[i], min_sr.hit_point);
                    float ndotwi = vec3_dot(light_dir, min_sr.normal);
                    if(ndotwi >= 0)
                    {
                        // Shadow test
                        bool in_shadow = false;
                        if(scene_lights.shadow[i] && min_sr.mat->shadow)
                        {
                            Ray shadow_ray;
                            vec3_copy(shadow_ray.origin, min_sr.hit_point);
                            vec3_copy(shadow_ray.direction, light_dir);
                            min_t = shadowIntersectTest(&scene_objects, shadow_ray);                            
                            float t;                            
                            if(scene_lights.light_types[i] == DIRECTIONAL)
                            {
                                t = TMAX;                                
                            }else if(scene_lights.light_types[i] == POINT)
                            {
                                vec3 light_to_hit_point;
                                PointLight* point_light_ptr = (PointLight*)(scene_lights.light_ptrs[i]);
                                vec3_sub(light_to_hit_point, point_light_ptr->point, min_sr.hit_point);
                                t = vec3_length(light_to_hit_point);
                            }
                            if(min_t < t)
                            {
                                in_shadow = true;
                            }                            
                        }
                        if(!in_shadow)
                        {
                            // Diffuse component
                            // kd*cd/PI * inc_radiance_cos
                            vec3 inc_radiance_cos, tmp;
                            getIncRadiance(tmp, scene_lights.light_types[i], scene_lights.light_ptrs[i]);
                            vec3_scale(inc_radiance_cos, tmp, ndotwi);
                            diffuseShading(radiance, ndotwi, inc_radiance_cos,  &min_sr);
                            if(min_sr.mat->mat_type == PHONG)
                            {
                                // Specular component
                                // TODO: add wo to ShadeRec
                                vec3 wo;
                                vec3_negate(wo, ray.direction);
                                specularShading(radiance, wo, light_dir, inc_radiance_cos, &min_sr);
                            }
                        }
                    }
                }
                */                
                vec3_add(color, color, radiance);                
            }else
            {
                vec3_add(color, color, bg_color);
            }
        }
        vec3_scale(color, color, 1.0f/num_samples);
        // divide color by max component if max component > 1
        float max;
        max = (color[0] > color[1]) ? color[0] : color[1];
        max = (max > color[2]) ? max : color[2];
        if(max > 1.0f)
        {
            vec3_scale(color, color, 1.0f/max);
        }
        image[i*3] = (char)(color[0] * 255.0f);
        image[i*3 + 1] = (char)(color[1] * 255.0f);
        image[i*3 + 2] = (char)(color[2] * 255.0f);
    }
    time(&endTime);
    double sec = difftime(endTime, startTime);
    printf("%f seconds.\n", sec);        
//#endif
    displayImage(window, viewport, image, frame_res_width, frame_res_height);
    free(image);
    freeSamples(&samples);
    freeSceneObjects(&scene_objects);
    freeSceneLights(&scene_lights);

    while(true)
    {

    }
    return 0;
}
