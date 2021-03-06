#pragma once
#include <float.h>
#include "vec.h"

const int MAX_OBJECTS = 1000;
const int INITIAL_NUM_OBJECTS = 20;
const int INITIAL_NUM_TRIANGLES = 200;
const int MAX_LIGHTS = 200;
const int MAX_MESH = 200;
const int MAX_CAUSTIC_OBJECTS = 10000; // NOTE: will probabily cause problems down the line
const int DEFAULT_MATERIAL = 40;
const int MAX_NAME_LENGTH = 128;
const int NAME_LENGTH = 32;

const int DEFAULT_NUM_TEXTURES = 20;

const float DEFAULT_GLOSSINESS = 100.0f;

// Directions
const vec3 ORIGIN = {0.0f, 0.0f, 0.0f};
const vec3 UP = {0.0f, 1.0f, 0.0f};
const vec3 DOWN = {0.0f, -1.0f, 0.0f};
const vec3 LEFT = {-1.0f, 0.0f, 0.0f};
const vec3 RIGHT = {1.0f, 0.0f, 0.0f};
const vec3 FORWARD = {0.0f, 0.0f, -1.0f};
const vec3 BACKWARD = {0.0f, 0.0f, 1.0f};
const vec3 JITTERED_UP = {0.0072f, 1.0f, 0.0034f};

// RGB values
const vec3 RED = {1.0f, 0.0f, 0.0f};
const vec3 GREEN = {0.0f, 1.0f, 0.0f};
const vec3 BLUE = {0.0f, 0.0f, 1.0f};
const vec3 WHITE = {1.0f, 1.0f, 1.0f};
const vec3 BLACK = {0.0f, 0.0f, 0.0f};
const vec3 YELLOW = {1.0f, 1.0f, 0.0f};
const vec3 CYAN = {0.0f, 1.0f, 1.0f};
const vec3 PINK = {1.0f, 0.0f, 1.0f};
const vec3 GREY = {0.5f, 0.5f, 0.5f};
const vec3 MED_ORCHID = {0.729, 0.333, 0.827};

// Camera
const float DEFAULT_FOCAL_LENGTH = 3.0f;
const float DEFAULT_LENS_RADIUS = 0.2f;

const float K_EPSILON = 0.000007f;
//const float K_EPSILON = 0.005f;
const float K_SMALLVALUE = 0.0000001f;
const float TMAX = FLT_MAX;
const float HUGEVALUE = 1.0E10;
const float K_FLAT_AABB = 0.001f;

const int MAX_BOUNDING_SPHERES = 100;
const int THETA_ROW = 26;
const int PHI_COLUMN = 52;

const double K_E = 2.718281828459045;
