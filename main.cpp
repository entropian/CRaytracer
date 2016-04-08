/*
  Things that can be swapped:
  camera --  perspective vs orthographic
  objects -- implicit surfaces vs explicit surfaces

  Features to consider:
  being able to scale the camera and viewplane 

  TODO:
  gamma correction
  consider how variables should be grouped into objects?

  Potential optimization:
  use float sample_x, sample_y instead of vec2 sample?
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
#include "shapes.h"
#include "glcode.h"
#include "sampling.h"
#include "camera.h"
#include "lights.h"
#include "constants.h"
#include "sceneobj.h"
#include "shading.h"
#include "init.h"

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


int main()
{
    int window_width = 900, window_height = 900;
    GLFWwindow* window = initWindow(window_width, window_height);
    glfwSetKeyCallback(window, keyCallback);
    
    GlViewport viewport;
    initViewport(&viewport);
    
    unsigned char *image;
    //int frame_res_width = 900, frame_res_height = 900;
    int frame_res_width = 360, frame_res_height = 360;    
    int num_pixels = frame_res_width * frame_res_height;
    image = (unsigned char*)calloc(num_pixels * 3, sizeof(char));
    
    // Samples
    srand(time(NULL));    
    const int num_samples = 64;
    const int num_sets = 83;
    Samples2D unit_square_samples = Samples2D_default;
    Samples2D disk_samples = Samples2D_default;
    Samples3D h_samples = Samples3D_default;
    genMultijitteredSamples(&unit_square_samples, num_samples, num_sets);
    //genRegularSamples(&unit_square_samples, num_samples, num_sets);
    mapSamplesToDisk(&disk_samples, &unit_square_samples);
    mapSamplesToHemisphere(&h_samples, &disk_samples, 1);
    
    // Ambient light
    // TODO: move this into scene lights
    vec3 amb_color = {1.0f, 1.0f, 1.0f};
    float amb_ls = 0.0f;
    bool amb_occlusion = true;

    // Scene 
    SceneLights scene_lights;
    SceneObjects scene_objects;
    initScene(&scene_objects, &scene_lights, num_samples, num_sets);

    // vec3 bg_color = {0.0f, 0.0f, 0.0f}; 
    vec3 bg_color;
    initBackgroundColor(bg_color, &scene_lights, BLACK);

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

    time_t startTime, endTime;
    time(&startTime);        
    unsigned sample_index = 0;    
    for(int i = 0; i < num_pixels; i++)
    {
        vec3 color = {0.0f, 0.0f, 0.0f};
        for(int p = 0; p < num_samples; p++)
        {
            vec2 sample, imageplane_coord;
            getNextSample2D(sample, &unit_square_samples);
            imageplane_coord[0] = -frame_length/2 + pixel_length * ((float)(i % frame_res_width) + sample[0]);
            imageplane_coord[1] = frame_height/2 - pixel_length * ((float)(i / frame_res_width) + sample[1]);
            
            Ray ray;
            switch(camera.camera_type)
            {
            case Pinhole:
            {
                calcRayPinhole(&ray, imageplane_coord, &camera);
            } break;
            case ThinLens:
            {
                vec2 disk_sample;
                getNextSample2D(disk_sample, &disk_samples);            
                calcRayThinLens(&ray, imageplane_coord, disk_sample, &camera);
            } break;
            }

            float min_t = TMAX;
            ShadeRec min_sr;
            min_t = intersectTest(&min_sr, &scene_objects, ray);

            // Shading
            if(min_t < TMAX)
            {                
                vec3 radiance = {0.0f, 0.0f, 0.0f};                
                if(min_sr.mat->mat_type == EMISSIVE)
                {
                    vec3_scale(radiance, min_sr.mat->ce, min_sr.mat->ke/1.0f);
                    maxToOne(radiance, radiance);
                }else
                {
                    // Ambient component
                    // ka*ca * amb_inc_radiance
                    vec3 amb_inc_radiance;
                    vec3_scale(amb_inc_radiance, amb_color, amb_ls);
                    vec3 reflectance;
                    vec3_scale(reflectance, min_sr.mat->ca, min_sr.mat->ka);

                    // Ambient Occlusion
                    if(amb_occlusion && AOTest(&h_samples, &scene_objects, &min_sr) < TMAX)
                    {
                        vec3 min_amb;
                        vec3_scale(min_amb, amb_color, 0.01f);
                        vec3_copy(radiance, min_amb);
                    }else
                    {
                        vec3_mult(radiance, amb_inc_radiance, reflectance);
                    }

                    for(int i = 0; i < scene_lights.num_lights; i++)
                    {
                        vec3 light_dir;
                        // Area light affected -- half implemented
                        getLightDir(light_dir, scene_lights.light_types[i], scene_lights.light_ptrs[i], &min_sr);
                        float ndotwi = vec3_dot(light_dir, min_sr.normal);
                        if(ndotwi > 0)
                        {
                            // Shadow test
                            bool in_shadow = false;
                            if(scene_lights.shadow[i] && min_sr.mat->shadow)
                            {
                                Ray shadow_ray;
                                vec3_copy(shadow_ray.origin, min_sr.hit_point);
                                vec3_copy(shadow_ray.direction, light_dir);
                                min_t = shadowIntersectTest(&scene_objects, shadow_ray);                            
                                // Area light affected
                                float t = calcLightDistance(scene_lights.light_types[i],
                                                            scene_lights.light_ptrs[i], min_sr.hit_point);
                                if(min_t < t)
                                {
                                    in_shadow = true;
                                }                            
                            }
                            if(!in_shadow)
                            {
                                if(scene_lights.light_types[i] == AREALIGHT)
                                {
                                    // TODO: Specular shading for area lights
                                    AreaLight* area_light_ptr = (AreaLight*)(scene_lights.light_ptrs[i]);
                                    areaLightShading(radiance, ndotwi, area_light_ptr, &min_sr);
                                }else if(scene_lights.light_types[i] == ENVLIGHT)
                                {
                                    EnvLight* env_light_ptr = (EnvLight*)(scene_lights.light_ptrs[i]);
                                    vec3 f, inc_radiance_cos, tmp;
                                    diffuseBRDF(f, &min_sr);
                                    getIncRadiance(tmp, scene_lights.light_types[i], scene_lights.light_ptrs[i]);
                                    vec3_scale(inc_radiance_cos, tmp, ndotwi);
                                    diffuseShading(radiance, inc_radiance_cos, &min_sr);
                                    float pdf = ndotwi / PI;
                                    vec3_scale(radiance, radiance, 1.0f/pdf);
                                }else
                                {
                                    // Diffuse component
                                    // kd*cd/PI * inc_radiance_cos
                                    vec3 inc_radiance_cos, tmp;
                                    getIncRadiance(tmp, scene_lights.light_types[i], scene_lights.light_ptrs[i]);
                                    vec3_scale(inc_radiance_cos, tmp, ndotwi);
                                    diffuseShading(radiance, inc_radiance_cos,  &min_sr);
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
                    }                    
                }
                vec3_add(color, color, radiance);
            }else
            {
                vec3_add(color, color, bg_color);
            }
        }
        vec3_scale(color, color, 1.0f/num_samples);
        maxToOne(color, color);
        image[i*3] = (char)(color[0] * 255.0f);
        image[i*3 + 1] = (char)(color[1] * 255.0f);
        image[i*3 + 2] = (char)(color[2] * 255.0f);
    }
    time(&endTime);
    double sec = difftime(endTime, startTime);
    printf("%f seconds.\n", sec);
    
    displayImage(window, viewport, image, frame_res_width, frame_res_height);

    // Clean up
    free(image);
    freeSamples2D(&unit_square_samples);
    freeSamples2D(&disk_samples);
    freeSamples3D(&h_samples);
    freeSceneObjects(&scene_objects);
    freeSceneLights(&scene_lights);

    while(true)
    {

    }
    return 0;
}
