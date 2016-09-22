#include <stdlib.h>
#include "util/vec.h"

typedef struct Photon_s
{
    vec3 position;
    short plane;
    unsigned char theta, phi;
    vec3 power;
}Photon;

