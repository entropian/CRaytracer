#pragma once

enum ObjectType
{
    SPHERE,
    PLANE,
    RECTANGLE,
    AABOX,
    TRIANGLE,
    OPENCYLINDER,
    DISK,
    TORUS,
    INSTANCED,
    COMPOUND
};

bool isGridObjType(const ObjectType obj_type)
{
    switch(obj_type)
    {
    case SPHERE:
    case RECTANGLE:
    case AABOX:
    case TRIANGLE:
    case OPENCYLINDER:
    case DISK:
    case TORUS:
        return true;
    }
    return false;
}

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
    case OPENCYLINDER:
        printf("OPENCYLINDER\n");
        break;
    case DISK:
        printf("DISK\n");
        break;
    case TORUS:
        printf("TORUS\n");
        break;
    }
}
