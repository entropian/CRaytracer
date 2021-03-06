#include "math.h"

float distanceSquared(const vec3 a, const vec3 b)
{
    vec3 displacement;
    vec3_sub(displacement, a, b);
    return vec3_dot(displacement, displacement);
}

void orthoNormalTransform(vec3 r, const vec3 u, const vec3 v, const vec3 w, const vec3 a)
{
    vec3 u_comp, v_comp, w_comp;
    vec3_scale(u_comp, u, a[0]);
    vec3_scale(v_comp, v, a[1]);
    vec3_scale(w_comp, w, a[2]);
    vec3_add(r, u_comp, v_comp);
    vec3_add(r, r, w_comp);    
}

void transposeTransform(vec3 r, const vec3 u, const vec3 v, const vec3 w, const vec3 a)
{
    //Vector3f(Dot(v, ss), Dot(v, ts), Dot(v, ns));
    vec3_assign(r, vec3_dot(a, u), vec3_dot(a, v), vec3_dot(a, w));
    /*
    r[0] = u[0]*a[0] + v[0]*a[1] + w[0]*a[2];
    r[1] = u[1]*a[0] + v[1]*a[1] + w[1]*a[2];
    r[2] = u[2]*a[0] + v[2]*a[1] + w[2]*a[2];
    */
}

void getVec3InLocalBasis(vec3 r, const vec3 a, const vec3 normal)
{
    vec3 u, v, w;
    vec3_copy(w, normal);
    vec3_cross(v, w, JITTERED_UP);
    vec3_normalize(v, v);
    vec3_cross(u, v, w);
    orthoNormalTransform(r, u, v, w, a);
}

void transformRay(Ray* dest_ray, const mat4 mat, const Ray src_ray)
{
    vec4 src_origin, src_direction;
    vec4 dest_origin, dest_direction;
    vec4FromVec3(src_origin, src_ray.origin, 1.0f);
    vec4FromVec3(src_direction, src_ray.direction, 0.0f);
    mat4_mult_vec4(dest_origin, mat, src_origin);
    mat4_mult_vec4(dest_direction, mat, src_direction);
    vec3FromVec4(dest_ray->origin, dest_origin);
    vec3FromVec4(dest_ray->direction, dest_direction);
}

void defaultInvTransform(mat4 r, const vec3 scaling, const vec3 axis, const float theta,
                         const vec3 translation)
{
    mat4 inv_translation, inv_rotation, inv_scale, tmp;
    mat4_translate(inv_translation, -translation[0], -translation[1], -translation[2]);
    mat4_scale_inverse(inv_scale, scaling);
    mat4_rotate(inv_rotation, axis, theta);
    mat4_mult(tmp, inv_scale, inv_rotation);
    mat4_mult(r, tmp, inv_translation);
}


void eulerAngToMat3(mat3 r, const vec3 euler_ang)
{
    float z_rad, y_rad, x_rad;
    z_rad = degToRad(euler_ang[2]);
    x_rad = degToRad(euler_ang[1]);
    y_rad = degToRad(euler_ang[0]);    
    
    mat3 z_rot, y_rot, x_rot, tmp;
    mat3_rotate_z(z_rot, z_rad);
    mat3_rotate_x(x_rot, x_rad);        
    mat3_rotate_y(y_rot, y_rad);    
    mat3_mult(tmp, x_rot, z_rot);
    mat3_mult(r, y_rot, tmp);
}

void eulerAngToMat4(mat4 r, const vec3 euler_ang)
{
    float z_rad, y_rad, x_rad;
    z_rad = degToRad(euler_ang[2]);
    x_rad = degToRad(euler_ang[1]);
    y_rad = degToRad(euler_ang[0]);    
    
    mat4 z_rot, y_rot, x_rot, tmp;
    mat4_rotate_z(z_rot, z_rad);
    mat4_rotate_x(x_rot, x_rad);        
    mat4_rotate_y(y_rot, y_rad);    
    mat4_mult(tmp, x_rot, z_rot);
    mat4_mult(r, y_rot, tmp);
}

__m128 fourKnotSplineSSE(__m128* x, __m128* k0, __m128* k1, __m128* k2, __m128* k3)
{

    __m128 coeff = _mm_set_ps1(-0.5f);
    __m128 tmp1 = _mm_mul_ps(*k0, coeff);
    coeff = _mm_set_ps1(1.5f);
    __m128 tmp2 = _mm_mul_ps(*k1, coeff);
    coeff = _mm_set_ps1(-1.5f);    
    __m128 tmp3 = _mm_mul_ps(*k2, coeff);
    coeff = _mm_set_ps1(0.5f);
    __m128 tmp4 = _mm_mul_ps(*k3, coeff);
    tmp1 = _mm_add_ps(tmp1, tmp2);
    tmp2 = _mm_add_ps(tmp3, tmp4);
    __m128 c3 = _mm_add_ps(tmp1, tmp2);

    coeff = _mm_set_ps1(-2.5f);
    tmp1 = _mm_mul_ps(coeff, *k1);
    coeff = _mm_set_ps1(2.0f);
    tmp2 = _mm_mul_ps(coeff, *k2);
    coeff = _mm_set_ps1(-0.5f);
    tmp3 = _mm_mul_ps(coeff, *k3);
    tmp1 = _mm_add_ps(*k0, tmp1);
    tmp2 = _mm_add_ps(tmp2, tmp3);
    __m128 c2 = _mm_add_ps(tmp1, tmp2);

    coeff = _mm_set_ps1(-1.0f);
    tmp1 = _mm_mul_ps(coeff, *k0);
    tmp1 = _mm_add_ps(tmp1, *k2);
    coeff = _mm_set_ps1(0.5f);
    __m128 c1 = _mm_mul_ps(coeff, tmp1);

    __m128 output = _mm_add_ps(_mm_mul_ps(_mm_add_ps(_mm_mul_ps(_mm_add_ps(_mm_mul_ps(c3, *x), c2), *x), c1), *x), *k1);
    return output;
}


/*
  Utility functions to find cubic and quartic roots.
  Copied from code for Ray Tracing from the Ground up.
  Author: Jochen Schwarze (schwarze@isa.de)
 */

#ifndef M_PI
#define M_PI PI
#endif

// You may have to experiment with EQN_EPS
// The original was 1e-9, but I use 1e-90  KS Dec 3, 2007

//#define     EQN_EPS     1e-9  
//#define     EQN_EPS     1e-30
//#define     EQN_EPS     1e-60
#define     EQN_EPS     1e-90

#define	IsZero(x)	((x) > -EQN_EPS && (x) < EQN_EPS)
/*
#ifndef CBRT
#define     cbrt(x)  ((x) > 0.0 ? pow((double)(x), 1.0/3.0) : \
			  		 ((x) < 0.0 ? -pow((double)-(x), 1.0/3.0) : 0.0))
#endif
*/
int solveQuadric(double c[3], double s[2])
{
    double p, q, D;

    /* normal form: x^2 + px + q = 0 */

    p = c[ 1 ] / (2 * c[ 2 ]);
    q = c[ 0 ] / c[ 2 ];

    D = p * p - q;

    if(IsZero(D))
    {
        s[ 0 ] = - p;
        return 1;
    }else if (D > 0)
    {
        double sqrt_D = (double)sqrt(D);

        s[ 0 ] =   sqrt_D - p;
        s[ 1 ] = - sqrt_D - p;
        return 2;
    }
    else /* if (D < 0) */
    {
        return 0;
    }
}

int solveCubic(double c[4], double s[3])
{
    int     i, num;
    double  sub;
    double  A, B, C;
    double  sq_A, p, q;
    double  cb_p, D;

    /* normal form: x^3 + Ax^2 + Bx + C = 0 */

    A = c[ 2 ] / c[ 3 ];
    B = c[ 1 ] / c[ 3 ];
    C = c[ 0 ] / c[ 3 ];

    /*  substitute x = y - A/3 to eliminate quadric term:
        x^3 +px + q = 0 */

    sq_A = A * A;
    p = 1.0f/3 * (- 1.0f/3 * sq_A + B);
    q = 1.0f/2 * (2.0f/27 * A * sq_A - 1.0f/3 * A * B + C);

    /* use Cardano's formula */

    cb_p = p * p * p;
    D = q * q + cb_p;

    if (IsZero(D)) {
		if (IsZero(q)) { /* one triple solution */
		    s[ 0 ] = 0;
		    num = 1;
		}
        else { /* one single and one double solution */
            double u = (double)cbrt(-q);
            s[ 0 ] = 2 * u;
            s[ 1 ] = - u;
            num = 2;
        }
    }
    else if (D < 0) { /* Casus irreducibilis: three real solutions */
		double phi = 1.0f/3 * (double)acos(-q / (double)sqrt(-cb_p));
		double t = 2.0f * (double)sqrt(-p);

		s[ 0 ] =   t * (double)cos(phi);
		s[ 1 ] = - t * (double)cos(phi + M_PI / 3);
		s[ 2 ] = - t * (double)cos(phi - M_PI / 3);
		num = 3;
    }
    else { /* one real solution */
		double sqrt_D = (double)sqrt(D);
		double u = (double)cbrt(sqrt_D - q);
		double v = - (double)cbrt(sqrt_D + q);

		s[ 0 ] = u + v;
		num = 1;
    }

    /* resubstitute */

    sub = 1.0f/3 * A;

    for (i = 0; i < num; ++i)
        s[ i ] -= sub;

    return num;
}

int solveQuartic(double c[5], double s[4])
{
    double  coeffs[4];
    double  z, u, v, sub;
    double  A, B, C, D;
    double  sq_A, p, q, r;
    int     i, num;

    /* normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0 */

    A = c[ 3 ] / c[ 4 ];
    B = c[ 2 ] / c[ 4 ];
    C = c[ 1 ] / c[ 4 ];
    D = c[ 0 ] / c[ 4 ];

    /*  substitute x = y - A/4 to eliminate cubic term:
        x^4 + px^2 + qx + r = 0 */

    sq_A = A * A;
    p = - 3.0f/8 * sq_A + B;
    q = 1.0f/8 * sq_A * A - 1.0f/2 * A * B + C;
    r = - 3.0f/256*sq_A*sq_A + 1.0f/16*sq_A*B - 1.0f/4*A*C + D;

    if (IsZero(r))
    {
		/* no absolute term: y(y^3 + py + q) = 0 */

		coeffs[ 0 ] = q;
		coeffs[ 1 ] = p;
		coeffs[ 2 ] = 0;
		coeffs[ 3 ] = 1;

		num = solveCubic(coeffs, s);

		s[ num++ ] = 0;
    }else
    {
		/* solve the resolvent cubic ... */

		coeffs[ 0 ] = 1.0f/2 * r * p - 1.0f/8 * q * q;
		coeffs[ 1 ] = - r;
		coeffs[ 2 ] = - 1.0f/2 * p;
		coeffs[ 3 ] = 1;
		
		(void) solveCubic(coeffs, s);

		/* ... and take the one real solution ... */

		z = s[ 0 ];

		/* ... to build two quadric equations */

		u = z * z - r;
		v = 2 * z - p;

		if (IsZero(u))
		    u = 0;
		else if (u > 0)
		    u = sqrt(u);
		else
		    return 0;

		if (IsZero(v))
		    v = 0;
		else if (v > 0)
		    v = sqrt(v);
		else
		    return 0;

		coeffs[ 0 ] = z - u;
		coeffs[ 1 ] = q < 0 ? -v : v;
		coeffs[ 2 ] = 1;

		num = solveQuadric(coeffs, s);

		coeffs[ 0 ]= z + u;
		coeffs[ 1 ] = q < 0 ? v : -v;
		coeffs[ 2 ] = 1;

		num += solveQuadric(coeffs, s + num);
	}

    /* resubstitute */

    sub = 1.0f/4 * A;

    for (i = 0; i < num; ++i)
		s[ i ] -= sub;

    return num;
}


