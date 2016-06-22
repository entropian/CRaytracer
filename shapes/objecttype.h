#pragma once

enum ObjectType
{
    SPHERE,
    PLANE,
    RECTANGLE,
    AABOX,
    TRIANGLE,
    GENERICOPENCYLINDER,
    DISK,
    GENERICTORUS,
    INSTANCED,
    COMPOUND,
    MESH_TRIANGLE
};

typedef struct
{
    ObjectType type;
    void* ptr;    
}Object_t;

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
    case MESH_TRIANGLE:
        printf("MESH_TRIANGLE\n");
        break;
    }
}
