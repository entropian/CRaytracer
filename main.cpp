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

int initSpheres(Sphere spheres[])
{
    // TODO: tidy up material assignment
    int sphere_count = 0;
    vec3_assign(spheres[sphere_count].center, -1.0f, 0.0f, -4.0f);
    spheres[sphere_count].radius = 1.0f;
    vec3_copy(spheres[sphere_count].mat.cd, PINK);
    vec3_copy(spheres[sphere_count].mat.ca, PINK);
    vec3_copy(spheres[sphere_count].mat.cs, PINK);    
    spheres[sphere_count].mat.kd = 0.6f;
    spheres[sphere_count].mat.ka = 1.0f;
    spheres[sphere_count].mat.ks = 0.4f;
    spheres[sphere_count].mat.exp = 5.0f;        
    spheres[sphere_count].mat.mat_type = Matte;
    sphere_count++;

    vec3_assign(spheres[sphere_count].center, 0.75f, 0.0f, -4.0f);
    spheres[sphere_count].radius = 1.0f;
    vec3_copy(spheres[sphere_count].mat.cd, CYAN);
    vec3_copy(spheres[sphere_count].mat.ca, CYAN);
    vec3_copy(spheres[sphere_count].mat.cs, CYAN);    
    spheres[sphere_count].mat.kd = 0.6f;
    spheres[sphere_count].mat.ka = 1.0f;
    spheres[sphere_count].mat.ks = 0.4f;
    spheres[sphere_count].mat.exp = 10.0f;            
    spheres[sphere_count].mat.mat_type = Matte;    
    sphere_count++;

    return sphere_count;
}

int main()
{
    int window_width = 1920, window_height = 1080;
    GLFWwindow* window = initWindow(1920, 1080);
    glfwSetKeyCallback(window, keyCallback);
    
    GlViewport viewport;
    initViewport(&viewport);
    
    unsigned char *image;
    int frame_res_width = 720, frame_res_height = 405;
    int num_pixels = frame_res_width * frame_res_height;
    image = (unsigned char*)calloc(num_pixels * 3, sizeof(char));

    vec3 bg_color = {35.0f/255.0f, 47.0f/255.0f, 47.0f/255.0f};

    // Objects
    Sphere spheres[MAX_SPHERES];
    int sphere_count = 0;
    sphere_count = initSpheres(spheres);

    // Ambient light
    vec3 amb_color = {1.0f, 1.0f, 1.0f};
    float amb_ls = 0.05f;

    // Directional light
    float intensity = 2.0f;
    vec3 direction = {10.0f, 10.0f, 10.0f};
    vec3_normalize(direction, direction);    
    DirectionalLight dir_light;
    assignDirLight(&dir_light, intensity, WHITE, direction);

    // Samples
    int num_samples = 16;
    int num_sets = 3;    
    Samples samples = Samples_default;
    //genRegularSamples(&samples, num_samples, num_sets);
    //genMultijitteredSamples(&samples, num_samples, num_sets);
    genHammersleySamples(&samples, num_samples, num_sets);
    mapSamplesToDisk(&samples);
    //mapSamplesToHemisphere(&samples, 5);

    // Camera
    Camera camera;
    initPinholeCameraDefault(&camera);

    vec3 position = {0.0f, 0.0f, 0.0f};
    vec3 look_point = {0.75f, 0.0f, -4.0f};
    vec3 up_vec = {0.0f, 1.0f, 0.0f};
    cameraLookAt(&camera, position, look_point, up_vec);

    // TODO: throw these in a struct for camera scaling
    float fov = 90.0f;
    float frame_length = 2.0f * (sin(fov/2.0f) * camera.focal_pt_dist);
    float frame_height = frame_length * (float)(frame_res_height)/(float)(frame_res_width);    
    float pixel_length = frame_length/(float)(frame_res_width);

    for(int i = 0; i < num_pixels; i++)
    {
        vec3 color = {0.0f, 0.0f, 0.0f};
        for(int p = 0; p < num_samples; p++)
        {
            vec2 sample, disk_sample, imageplane_coord;
            getNextSample(sample, &samples);
            imageplane_coord[0] = -frame_length/2 + pixel_length * ((float)(i % frame_res_width) + sample[0]);            
            imageplane_coord[1] = frame_height/2 - pixel_length * ((float)(i / frame_res_width) + sample[1]);
            
            Ray ray;
            switch(camera.camera_type)
            {
            case Pinhole:
                calcRayCam(&ray, imageplane_coord, &camera);
                break;
            case ThinLens:
                getNextDiskSample(disk_sample, &samples);            
                calcRayThinLens(&ray, imageplane_coord, disk_sample, &camera);
                break;
            }
            
            float min_t = TMAX,  tmp_t = TMAX;
            ShadeRec min_sr, tmp_sr;
            for(int i = 0; i < sphere_count; i++)
            {
                tmp_t = rayIntersectSphere(&tmp_sr, &(spheres[i]), ray);
                if(tmp_t < min_t)
                {
                    min_t = tmp_t;
                    min_sr = tmp_sr;
                }
            }
            if(min_t < TMAX)
            {
                // TODO: multiple lights, different materials
                vec3 radiance;
                // Ambient component
                // ka*ca * amb_inc_radiance
                vec3 amb_inc_radiance;
                vec3_scale(amb_inc_radiance, amb_color, amb_ls);
                vec3 reflectance;
                vec3_scale(reflectance, min_sr.mat->ca, min_sr.mat->ka);
                vec3_mult(radiance, amb_inc_radiance, reflectance);
                
                float ndotwi = vec3_dot(dir_light.direction, min_sr.normal);
                if(ndotwi >= 0)
                {
                    // TODO: write out the math expressions in the comments
                    // Diffuse component
                    // kd*cd/PI * inc_radiance_cos
                    vec3_scale(reflectance, min_sr.mat->cd, min_sr.mat->kd);
                    vec3 f;
                    vec3_scale(f, reflectance, 1.0f/(float)PI);
                    vec3 tmp, inc_radiance_cos;                                    
                    vec3_scale(tmp, dir_light.color, dir_light.intensity);
                    vec3_scale(inc_radiance_cos, tmp, ndotwi);
                    vec3_mult(tmp, inc_radiance_cos, f);
                    vec3_add(radiance, radiance, tmp);


                    // Specular
                    // ks*cs * cos(theta)^exp * inc_radiance_cos
                    // theta = angle between light and reflected view vector
                    vec3 reflect_dir;
                    // TODO: add wo to ShadeRec, optimize
                    // reflect vector: 2*dot(n, v)*n - v
                    vec3 wo;
                    vec3_negate(wo, ray.direction);
                    float ndotwo = vec3_dot(min_sr.normal, wo);
                    vec3_scale(tmp, min_sr.normal, ndotwo * 2);
                    vec3_add(reflect_dir, ray.direction, tmp);
                    vec3_normalize(reflect_dir, reflect_dir);
                    float rdotwi = vec3_dot(reflect_dir, dir_light.direction);
                    vec3_scale(tmp, min_sr.mat->cs, min_sr.mat->ks * pow(rdotwi, min_sr.mat->exp));
                    vec3_mult(tmp, inc_radiance_cos, tmp);
                    vec3_add(radiance, radiance, tmp);
                }
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
    displayImage(window, viewport, image, frame_res_width, frame_res_height);
    free(image);
    freeSamples(&samples);

    while(true)
    {

    }
    return 0;
}
