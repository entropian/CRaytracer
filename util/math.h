#pragma once

#include "vec.h"
#include "constants.h"

void orthoNormalTransform(vec3 r, const vec3 u, const vec3 v, const vec3 w, const vec3 a)
{
    vec3 u_comp, v_comp, w_comp;
    vec3_scale(u_comp, u, a[0]);
    vec3_scale(v_comp, v, a[1]);
    vec3_scale(w_comp, w, a[2]);
    vec3_add(r, u_comp, v_comp);
    vec3_add(r, r, w_comp);    
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

#ifndef CBRT
#define     cbrt(x)  ((x) > 0.0 ? pow((float)(x), 1.0/3.0) : \
			  		 ((x) < 0.0 ? -pow((float)-(x), 1.0/3.0) : 0.0))
#endif

int solveQuadric(float c[3], float s[2])
{
    float p, q, D;

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
        float sqrt_D = (float)sqrt(D);

        s[ 0 ] =   sqrt_D - p;
        s[ 1 ] = - sqrt_D - p;
        return 2;
    }
    else /* if (D < 0) */
    {
        return 0;
    }
}

int solveCubic(float c[4], float s[3])
{
    int     i, num;
    float  sub;
    float  A, B, C;
    float  sq_A, p, q;
    float  cb_p, D;

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
            float u = (float)cbrt(-q);
            s[ 0 ] = 2 * u;
            s[ 1 ] = - u;
            num = 2;
        }
    }
    else if (D < 0) { /* Casus irreducibilis: three real solutions */
		float phi = 1.0f/3 * (float)acos(-q / (float)sqrt(-cb_p));
		float t = 2.0f * (float)sqrt(-p);

		s[ 0 ] =   t * (float)cos(phi);
		s[ 1 ] = - t * (float)cos(phi + M_PI / 3);
		s[ 2 ] = - t * (float)cos(phi - M_PI / 3);
		num = 3;
    }
    else { /* one real solution */
		float sqrt_D = (float)sqrt(D);
		float u = (float)cbrt(sqrt_D - q);
		float v = - (float)cbrt(sqrt_D + q);

		s[ 0 ] = u + v;
		num = 1;
    }

    /* resubstitute */

    sub = 1.0f/3 * A;

    for (i = 0; i < num; ++i)
        s[ i ] -= sub;

    return num;
}

int solveQuartic(float c[5], float s[4])
{
    float  coeffs[4];
    float  z, u, v, sub;
    float  A, B, C, D;
    float  sq_A, p, q, r;
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

