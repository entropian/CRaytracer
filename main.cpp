/*
  Things that can be swapped:
  camera --  perspective vs orthographic
  objects -- implicit surfaces vs explicit surfaces
  sampling -- different sampling algorithms

  Need a struct for recording intersect point info,
  like position, normal, brdf, etc.

  TODO:
  gamma correction
  consider how variables should be grouped into objects?
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

static const float TMAX = 1000.0f;
static const float k_epsilon = 0.000001f;
static const int MAX_SPHERES = 1000;

vec3 cam_position = {0.0f, 0.0f, 0.0f};
vec3 cam_orientation = {0.0f, 0.0f, 0.0f};


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

void generateRegularSamples(vec2 *samples, const int num_samples)
{
    int samples_per_row = (int)sqrt(num_samples);
    int sample_index;
    for(int i = 0; i < samples_per_row; i++)
    {
        for(int j = 0; j < samples_per_row; j++)
        {
            // NOTE: instead of computing a sample index, could probably just increment an index variable.
            sample_index = i*samples_per_row + j;
            samples[sample_index][0] = (i + 0.5f) / samples_per_row;
            samples[sample_index][1] = (j + 0.5f) / samples_per_row;
        }
    }
}

struct Ray
{
    vec3 origin;
    vec3 direction;
};

struct ShadeRec
{
    bool hit_status;
    vec3 hit_point;
    vec3 normal;
    vec3 color;                // NOTE: Temporary 
};

void fillShadeRecSphere(ShadeRec *sr, const Sphere sphere, const Ray ray, const float t)
{
    sr->hit_status = true;
    vec3 displacement;
    vec3_scale(displacement, ray.direction, t);
    vec3_add(sr->hit_point, ray.origin, displacement);
    vec3_sub(displacement, sr->hit_point, sphere.center);
    vec3_normalize(sr->normal, displacement);
    vec3_copy(sr->color, sphere.color);
}

float rayIntersectSphere(ShadeRec *sr, const Sphere sphere, const Ray ray)
{
    // The analytic solution is so much better than the shitty loop from before
    // though I probably should have used smaller increment for the loop
    float a = vec3_dot(ray.direction, ray.direction);
    vec3 origin_to_center, tmp;
    vec3_sub(origin_to_center, ray.origin, sphere.center);
    vec3_scale(tmp, origin_to_center, 2);
    float b = vec3_dot(tmp, ray.direction);
    float c = vec3_dot(origin_to_center, origin_to_center) - sphere.radius*sphere.radius;
    float disc = b*b - 4*a*c;

    if(disc < 0.0f)
    {
        return TMAX;
    }else
    {
        float e = sqrt(disc);
        float denom = 2*a;
        float t = (-b - e)/denom;
        if(t > k_epsilon)
        {
            fillShadeRecSphere(sr, sphere, ray, t);            
            return t;
        }

        t = (-b + e)/denom;
        if(t > k_epsilon)
        {
            fillShadeRecSphere(sr, sphere, ray, t);            
            return t;
        }
    }
    return TMAX;
}

int initSpheres(Sphere spheres[])
{
    int sphere_count = 0;
    vec3_assign(spheres[sphere_count].center, 0.0f, 0.0f, -4.0f);
    spheres[sphere_count].radius = 1.0f;
    vec3_assign(spheres[sphere_count].color, 255.0f, 0.0f, 0.0f);
    sphere_count++;    

    vec3_assign(spheres[sphere_count].center, 0.75f, 0.0f, -4.0f);
    spheres[sphere_count].radius = 1.0f;
    vec3_assign(spheres[sphere_count].color, 0.0f, 255.0f, 0.0f);
    sphere_count++;

    return sphere_count;
}

int main()
{
    // TODO: sort out window size vs resolution
    int width = 720, height = 405;        
    GLFWwindow* window = initWindow(1920, 1080);
    glfwSetKeyCallback(window, keyCallback);
    
    GlViewport viewport;
    initViewport(&viewport);
    
    unsigned char *image;
    image = (unsigned char*)calloc(width * height * 3, sizeof(char));

    /* TODO: Moving the camera
       Should be done after sphere is implemented so that we can see the result
       Two Vec3's. One for translation and one for euler angle. Done?
       First calculate ray direction, then rotate it with euler angle.
       Rotate and translate the position of the pixel to use as ray starting point.
       Implement euler-angle-to-matrix conversion.
     */

    /*
      Camera related things: fov, resolution, position, orientation.
      Be cool and throw them into a class?
      What design pattern is this?
     */
    /*
      Array of objects
      Array of lights?
     */

    int num_pixels = width * height;
    vec3 focal_point = {0.0f, 0.0f, 1.0f};
    float fov = 90.0f;
    float focal_dist = 1.0f;
    float frame_length = 2.0f * (sin(fov/2.0f) * focal_dist);
    float frame_height = frame_length * static_cast<float>(height)/static_cast<float>(width);    
    float pixel_length = frame_length/static_cast<float>(width);
    int num_samples = 4;
    int samples_per_row = (int)sqrt(num_samples);
    float sample_length = pixel_length / static_cast<float>(samples_per_row);

    vec2 *samples = (vec2*)malloc(sizeof(vec2) * num_samples);
    generateRegularSamples(samples, num_samples);

    Sphere spheres[MAX_SPHERES];
    int sphere_count = 
    sphere_count = initSpheres(spheres);

    for(int i = 0; i < num_pixels; i++)
    {
        vec3 color = {0.0f, 0.0f, 0.0f};
        for(int p = 0; p < num_samples; p++)
        {
            // (-length/2 + p_len * (n mod width) + p_len/2, height/2 - p_len(n / height) - p_len/2)
            float x = -frame_length/2 + pixel_length * ((float)(i % width) + samples[p][0]);
            float y = frame_height/2 - pixel_length * ((float)(i / width) - samples[p][1]);
            
            Ray ray;
            vec3_assign(ray.origin, x, y, 0.0f);
            vec3_add(ray.origin, ray.origin, cam_position);
            vec3_sub(ray.direction, ray.origin, focal_point);
            vec3_normalize(ray.direction, ray.direction);

            float min_t = tmp_t = TMAX;
            ShadeRec min_sr, tmp_sr;
            for(int i = 0; i < sphere_count; i++)
            {
                tmp_t = rayIntersectSphere(&tmp_sr, spheres[i], ray);
                if(tmp_t < min_t)
                {
                    min_t = tmp_t;
                    min_sr = tmp_sr;
                }
            }
            if(min_t < TMAX)
            {
                vec3_add(color, color, min_sr.color);
            }
        }
        vec3_scale(color, color, 1.0f/num_samples);
        image[i*3] = (char)(color[0]);
        image[i*3 + 1] = (char)(color[1]);
        image[i*3 + 2] = (char)(color[2]);
    }
    displayImage(window, viewport, image, width, height);

    while(true)
    {

    }
    return 0;
}
