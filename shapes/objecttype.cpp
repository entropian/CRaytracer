#include "objecttype.h"

void printObjType(const ObjectType obj_type)
{
    switch(obj_type)
    {
    case SPHERE:
        printf("SPHERE\n");
        break;
    case RECTANGLE:
        printf("RECTANGLE\n");
        break;
    case AABOX:
        printf("AABOX\n");
        break;
    case TRIANGLE:
        printf("TRIANGLE\n");
        break;
    case GENERICOPENCYLINDER:
        printf("GENERICOPENCYLINDER\n");
        break;
    case DISK:
        printf("DISK\n");
        break;
    case GENERICTORUS:
        printf("TORUS\n");
        break;
    case FLAT_TRIANGLE:
        printf("FLAT_TRIANGLE\n");
        break;
    case SMOOTH_TRIANGLE:
        printf("SMOOTH_TRIANGLE\n");
        break;        
    }
}
