#pragma once
#include <stdio.h>

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
    FLAT_TRIANGLE,
    SMOOTH_TRIANGLE
};

typedef struct
{
    ObjectType type;
    void* ptr;    
}Object_t;


void printObjType(const ObjectType obj_type);

