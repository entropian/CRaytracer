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
