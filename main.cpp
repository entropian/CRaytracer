/*
  Things that can be swapped:
  camera --  perspective vs orthographic
  objects -- implicit surfaces vs explicit surfaces
  sampling -- different sampling algorithms

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
#include "sampling.h"
#include "camera.h"

static const float TMAX = 1000.0f;
static const float k_epsilon = 0.000001f;
static const int MAX_SPHERES = 1000;

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

typedef struct Ray
{
    vec3 origin;
    vec3 direction;
} Ray;

typedef struct ShadeRec
{
    bool hit_status;
    vec3 hit_point;
    vec3 normal;
    vec3 color;                // NOTE: Temporary 
} ShadeRec;


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

void calcRayCam(Ray *ray, const float sample_x, const float sample_y, const Camera *camera)
{
    vec3 sample_loc_cam, tmp_vec_1, tmp_vec_2, focal_point;
    vec3_scale(tmp_vec_1, camera->x_axis, sample_x);
    vec3_scale(tmp_vec_2, camera->y_axis, sample_y);
    vec3_add(sample_loc_cam, tmp_vec_1, tmp_vec_2);
    vec3_scale(focal_point, camera->z_axis, camera->focal_dist);    // focal_point is rotated but not translated
    vec3_sub(ray->direction, sample_loc_cam, focal_point);
    vec3_normalize(ray->direction, ray->direction);
    vec3_add(ray->origin, focal_point, camera->position);
}

void transformPoint(vec3 r, const vec3 a, const vec3 x, const vec3 y, const vec3 z, const vec3 cam_position)
{
    vec3 tmp_vec, result_vec;
    result_vec[0] = vec3_dot(a, x[0], y[0], z[0]);
    result_vec[1] = vec3_dot(a, x[1], y[1], z[1]);
    result_vec[2] = vec3_dot(a, x[2], y[2], z[2]);    
    vec3_add(result_vec, result_vec, cam_position);
    vec3_copy(r, result_vec);
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
    int num_samples = 16;
    int num_sets = 3;

    // TODO: hammersley sampling doesn't seem to work horizontally
    Samples samples = Samples_default;
    //genRegularSamples(&samples, num_samples, num_sets);
    genMultijitteredSamples(&samples, num_samples, num_sets);
    //genHammersleySamples(&samples, num_samples, num_sets);
    //mapSamplesToDisk(&samples);
    //mapSamplesToHemisphere(&samples, 5);


    Camera camera;
    initCameraDefault(&camera);
    vec3 look_point = {0.75f, 0.0f, -4.0f};
    vec3 up_vec = {0.0f, 1.0f, 0.0f};
    cameraLookAt(&camera, cam_position, look_point, up_vec);

    Sphere spheres[MAX_SPHERES];
    int sphere_count = 0;
    sphere_count = initSpheres(spheres);

    for(int i = 0; i < num_pixels; i++)
    {
        vec3 color = {0.0f, 0.0f, 0.0f};
        vec2 sample;
        for(int p = 0; p < num_samples; p++)
        {
            getNextSample(sample, &samples);
            float x = -frame_length/2 + pixel_length * ((float)(i % width) + sample[0]);
            float y = frame_height/2 - pixel_length * ((float)(i / width) + sample[1]);            
            Ray ray;
            calcRayCam(&ray, x, y, &camera);
            
            float min_t = TMAX,  tmp_t = TMAX;
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
    free(image);
    freeSamples(&samples);

    while(true)
    {

    }
    return 0;
}
