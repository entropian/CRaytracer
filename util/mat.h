#pragma once

#include "constants.h"
#include "vec.h"

typedef vec2 mat2[2];
typedef vec3 mat3[3];
typedef vec4 mat4[4];

void mat2_identity(mat2 r)
{
    r[0][0] = 1.0f;
    r[0][1] = 0.0f;
    r[1][0] = 0.0f;
    r[1][1] = 1.0f;
}

void mat3_identity(mat3 r)
{
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            r[i][j] = 0.0f;
        }
        r[i][i] = 1.0f;
    }
}

void mat4_identity(mat4 r)
{
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            r[i][j] = 0.0f;
        }
        r[i][i] = 1.0f;
    }
}

void initMat2(mat2 r, const vec2 a, const vec2 b)
{
    r[0][0] = a[0];
    r[0][1] = a[1];
    r[1][0] = b[0];
    r[1][1] = b[1];
}

void initMat3(mat3 r, const vec3 a, const vec3 b, const vec3 c)
{
    for(int i = 0; i < 3; i++)
    {
        r[0][i] = a[i];
        r[1][i] = b[i];
        r[2][i] = c[i];
    }
}

void initMat4(mat4 r, const vec3 a, const vec3 b, const vec3 c)
{
    mat4_identity(r);
    for(int i = 0; i < 3; i++)
    {
        r[0][i] = a[i];
        r[1][i] = b[i];
        r[2][i] = c[i];        
    }
}

void initMat4(mat4 r, const vec4 a, const vec4 b, const vec4 c, const vec4 d)
{
    for(int i = 0; i < 4; i++)
    {
        r[0][i] = a[i];
        r[1][i] = b[i];
        r[2][i] = c[i];
        r[3][i] = d[i];
    }
}

void mat2_copy(mat2 r, const mat2 a)
{
    r[0][0] = a[0][0];
    r[0][1] = a[0][1];
    r[1][0] = a[1][0];
    r[1][1] = a[1][1];
}

void mat3_copy(mat3 r, const mat3 a)
{
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            r[i][j] = a[i][j];
        }
    }
}

void mat4_copy(mat4 r, const mat4 a)
{
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            r[i][j] = a[i][j];
        }
    }
}

void mat2_transpose(mat2 r, const mat2 a)
{
    r[0][0] = a[0][0];
    r[0][1] = a[1][0];
    r[1][0] = a[0][1];
    r[1][1] = a[1][1];
}

void mat3_transpose(mat3 r, const mat3 a)
{
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            r[i][j] = a[j][i];
        }
    }
}

void mat4_transpose(mat4 r, const mat4 a)
{
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            r[i][j] = a[j][i];
        }
    }
}

void mat2_mult(mat2 r, const mat2 a, const mat2 b)
{
    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 2; j++)
        {
            r[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j];
        }
    }
}

void mat3_mult(mat3 r, const mat3 a, const mat2 b)
{
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            r[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j] + a[i][2]*b[2][j];
        }
    }
}

void mat4_mult(mat4 r, const mat4 a, const mat4 b)
{
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            r[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j] + a[i][2]*b[2][j] + a[i][3]*b[3][j];
        }
    }
}

void mat2_mult_vec2(vec2 r, const mat2 a, const vec2 b)
{
    r[0] = vec2_dot(a[0], b);
    r[1] = vec2_dot(a[1], b);
}

void mat3_mult_vec3(vec3 r, const mat3 a, const vec3 b)
{
    for(int i = 0; i < 3; i++)
    {
        r[i] = vec3_dot(a[i], b);
    }
}

void mat4_mult_vec4(vec4 r, const mat4 a, const vec4 b)
{
    for(int i = 0; i < 4; i++)
    {
        r[i] = vec4_dot(a[i], b);
    }
}

void mat2_rotate(mat2 r, const float theta)
{
    float sin_theta = (float)sin(theta);
    float cos_theta = (float)cos(theta);    
    r[0][0] = cos_theta;
    r[0][1] = sin_theta;
    r[1][0] = -cos_theta;
    r[1][1] = sin_theta;
}

void mat3_rotate_x(mat3 r, const float theta)
{
    float sin_theta = (float)sin(theta);
    float cos_theta = (float)cos(theta);
    mat3_identity(r);
    r[1][1] = cos_theta;
    r[1][2] = sin_theta;
    r[2][1] = -sin_theta;
    r[2][2] = cos_theta;
}

void mat4_rotate_x(mat4 r, float theta)
{
    float sin_theta = (float)sin(theta);
    float cos_theta = (float)cos(theta);
    mat4_identity(r);
    r[1][1] = cos_theta;
    r[1][2] = sin_theta;
    r[2][1] = -sin_theta;
    r[2][2] = cos_theta;
}

void mat3_rotate_y(mat3 r, const float theta)
{
    float sin_theta = (float)sin(theta);
    float cos_theta = (float)cos(theta);
    mat3_identity(r);
    r[0][0] = cos_theta;
    r[0][2] = sin_theta;
    r[2][0] = -cos_theta;
    r[2][2] = sin_theta;
}

void mat4_rotate_y(mat4 r, const float theta)
{
    float sin_theta = (float)sin(theta);
    float cos_theta = (float)cos(theta);
    mat4_identity(r);
    r[0][0] = cos_theta;
    r[0][2] = sin_theta;
    r[2][0] = -sin_theta;
    r[2][2] = cos_theta;
}

void mat3_rotate_z(mat3 r, const float theta)
{
    float sin_theta = (float)sin(theta);
    float cos_theta = (float)cos(theta);
    mat3_identity(r);
    r[0][0] = cos_theta;
    r[0][1] = sin_theta;
    r[1][0] = -sin_theta;
    r[1][1] = cos_theta;
}

void mat4_rotate_z(mat4 r, const float theta)
{
    float sin_theta = (float)sin(theta);
    float cos_theta = (float)cos(theta);
    mat4_identity(r);
    r[0][0] = cos_theta;
    r[0][1] = sin_theta;
    r[1][0] = -sin_theta;
    r[1][1] = cos_theta;
}

void mat4_rotate(mat4 r, const vec3 axis, const float theta)
{
    // axis as z axis for orthonormal basis
    vec3 x_axis, y_axis;
    vec3_cross(x_axis, axis, JITTERED_UP);
    vec3_normalize(x_axis, x_axis);
    vec3_cross(y_axis, x_axis, axis);

    mat4 pre, inv_pre, z_rotation, tmp;
    initMat4(pre, x_axis, y_axis, axis);
    mat4_transpose(inv_pre, pre);
    mat4_rotate_z(z_rotation, theta);    

    mat4_mult(tmp, pre, z_rotation);
    mat4_mult(r, tmp, inv_pre);
}

void mat3_invert_rotation(mat3 r, const mat3 a)
{
    mat3_transpose(r, a);
}

void mat4_invert_rotation(mat4 r, const mat4 a)
{
    mat4_transpose(r, a);
}

void mat2_scale(mat2 r, const float x, const float y)
{
    r[0][0] = x;
    r[1][1] = y;
    r[0][1] = r[1][0] = 0.0f;
}

void mat3_scale(mat3 r, const float x, const float y, const float z)
{
    mat3_identity(r);
    r[0][0] = x;
    r[1][1] = y;
    r[2][2] = z;
}

void mat4_scale(mat4 r, const float x, const float y, const float z)
{
    mat4_identity(r);
    r[0][0] = x;
    r[1][1] = y;
    r[2][2] = z;
}

void mat2_scale_inverse(mat2 r, const vec2 a)
{
    r[0][0] = 1.0f / a[0];
    r[0][1] = 0.0f;
    r[1][0] = 0.0f;
    r[1][1] = 1.0f / a[1];
}

void mat3_scale_inverse(mat3 r, const vec3 a)
{
    mat3_identity(r);
    r[0][0] = 1.0f / a[0];
    r[1][1] = 1.0f / a[1];
    r[2][2] = 1.0f / a[2];    
}

void mat4_scale_inverse(mat4 r, const vec3 a )
{
    mat4_identity(r);
    r[0][0] = 1.0f / a[0];
    r[1][1] = 1.0f / a[1];
    r[2][2] = 1.0f / a[2];
}

void mat2_invert_scale(mat2 r, const mat2 a)
{
    r[0][0] = 1.0f / a[0][0];
    r[1][1] = 1.0f / a[1][1];
    r[0][1] = r[1][0] = 0.0f;
}

void mat3_invert_scale(mat3 r, const mat3 a)
{
    mat3_identity(r);
    for(int i = 0; i < 3; i++)
    {
        r[i][i] = 1.0f / a[i][i];
    }
}

void mat4_invert_scale(mat4 r, const mat4 a)
{
    mat4_identity(r);
    for(int i = 0; i < 3; i++)
    {
        // TODO: add assert
        r[i][i] = 1.0f / a[i][i];
    }
}

void mat3_translate(mat3 r, const float x, const float y)
{
    mat3_identity(r);
    r[0][2] = x;
    r[1][2] = y;
}

void mat4_translate(mat4 r, const float x, const float y, const float z)
{
    mat4_identity(r);
    r[0][3] = x;
    r[1][3] = y;
    r[2][3] = z;
}

void mat3_invert_translation(mat3 r, const mat3 a)
{
    mat3_identity(r);
    r[0][2] = -a[0][2];
    r[1][2] = -a[1][2];
}

void mat4_invert_translation(mat4 r, const mat4 a)
{
    mat4_identity(r);
    r[0][3] = -a[0][3];
    r[1][3] = -a[1][3];
    r[2][3] = -a[2][3];    
}
