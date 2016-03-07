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
#include "glcode.h"

static const float TMAX = 1000.0f;
static const float k_epsilon = 0.000001f;

Vec3 cam_position(0.0f, 0.0f, 0.0f);
Vec3 cam_orientation(0.0f, 0.0f, 0.0f);


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

void generateRegularSamples(Vec2 *samples, const int num_samples)
{
    int samples_per_row = static_cast<int>(sqrt(num_samples));
    for(int i = 0; i < samples_per_row; i++)
    {
        for(int j = 0; j < samples_per_row; j++)
        {
            samples[i*samples_per_row + j] = Vec2((i + 0.5f) / samples_per_row,
                                                  (j + 0.5f) / samples_per_row);
        }
    }
}

struct Ray
{
    Vec3 origin;
    Vec3 direction;
};

struct ShadeRec
{
    bool hit_status;
    Vec3 hit_point;
    Vec3 normal;
};

struct Sphere
{
    Vec3 center;
    Vec3 color;
    float radius;
};

float rayIntersectSphere(Sphere sphere, Ray ray)
{
    // The analytic solution is so much better than the shitty loop from before
    // though I probably should have used smaller increment for the loop
    float a = dot(ray.direction, ray.direction);
    float b = dot((ray.origin - sphere.center)*2, ray.direction);
    float c = norm2(ray.origin - sphere.center) - sphere.radius*sphere.radius;
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
            return t;
        }

        t = (-b + e)/denom;
        if(t > k_epsilon)
        {
            return t;
        }
    }
    return TMAX;
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

    /*
      Center of the frame is at (0, 0, 0)
      Assume the default focal distance is 1, and default aspect ratio is 16:9
      Assume resolution is 1280 x 720
      Assume default horizontal fov is 90 degrees
      The length of the frame is 2(sin(fov/2) * focal distance) = 2(sin(45)) = sqrt(2)
      The height of the frame is length * 9/16.
      The area of a pixel is length/1280 * height/720. length/1280 == height 720      
      The center of the top left pixel is (-length/2 + length/1280/2, height/2 - height/720/2)
      let p_len = length/1280
      The center of the nth pixel is (-length/2 + p_len * (n mod 1280) + p_len/2, height/2 - p_len(n / 720) - p_len/2)
      Focal point is (0, 0, 1).
      Ray direction for the nth pixel is:
      (-l/2 + p_len * (n mod 1280) + p_len/2, h/2 - p_len(n / 720) - p_len/2), 0) - (0, 0, 1)
      Ray origin is at the pixel on the frame.

     */

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
    float fov = 90.0f;
    float focal_dist = 1.0f;
    float frame_length = 2.0f * (sin(fov/2.0f) * focal_dist);
    float frame_height = frame_length * static_cast<float>(height)/static_cast<float>(width);    
    float pixel_length = frame_length/static_cast<float>(width);
    int num_samples = 4;
    int samples_per_row = sqrt(static_cast<float>(4));
    float sample_length = pixel_length / static_cast<float>(samples_per_row);

    Sphere sphere1;
    sphere1.center = Vec3(0.0f, 0.0f, -4.0f);
    sphere1.radius = 1.0f;
    sphere1.color = Vec3(255.0f, 0.0f, 0.0f);

    Sphere sphere2;
    sphere2.center = Vec3(0.75f, 0.0f, -4.0f);
    sphere2.radius = 1.0f;
    sphere2.color = Vec3(0.0f, 255.0f, 0.0f);

    for(int i = 0; i < num_pixels; i++)
    {
        Vec3 color;
        for(int p = 0; p < samples_per_row; p++)
        {
            for(int q = 0; q < samples_per_row; q++)
            {
                // (-length/2 + p_len * (n mod width) + p_len/2, height/2 - p_len(n / height) - p_len/2)
                float x = -frame_length/2 + pixel_length * static_cast<float>(i % width) + sample_length * p;
                float y = frame_height/2 - pixel_length * static_cast<float>(i / width) - sample_length * q;
                Ray ray;
                ray.origin = Vec3(x, y, 0.0f) + cam_position;;
                ray.direction = normalize(ray.origin - Vec3(0.0f, 0.0f, 1.0f));

                float sphere1_intersect = rayIntersectSphere(sphere1, ray);
                float sphere2_intersect = rayIntersectSphere(sphere2, ray);

                if(sphere1_intersect < sphere2_intersect)
                {
                    color += sphere1.color;
                }else if(sphere2_intersect  < sphere1_intersect)
                {
                    //image[i*3 + 2] = static_cast<unsigned char>(sphere2.color[2]);
                    color += sphere2.color;
                }
            }
        }
        color /= static_cast<float>(num_samples);
        image[i*3] = static_cast<float>(color[0]);
        image[i*3 + 1] = static_cast<float>(color[1]);
        image[i*3 + 2] = static_cast<float>(color[2]);                    
    }

    displayImage(window, viewport, image, width, height);

    while(true)
    {

    }
    return 0;
}
