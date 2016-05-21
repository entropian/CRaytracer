#pragma onece

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "util/vec.h"
#include "util/constants.h"
#include "shapes/shapes.h"
#include "shapes/instanced.h"
#include "lights.h"

void printSphere(const Sphere* sphere)
{
    printf("shadow %s\n", sphere->shadow ? "true" : "false");
    printf("min_theta = %f\n", sphere->min_theta);
    printf("max_theta = %f\n", sphere->max_theta);
    printf("phi = %f\n", sphere->phi);
    printf("radius = %f\n", sphere->radius);
    printf("center = %f, %f, %f\n", sphere->center[0], sphere->center[1], sphere->center[2]);    
}

void printTorus(const Torus* torus)
{
    printf("shadow %s\n", torus->shadow ? "true" : "false");
    printf("swept_radius = %f\n", torus->swept_radius);
    printf("tube_radius = %f\n", torus->tube_radius);
    printf("phi = %f\n", torus->phi);
    AABB aabb = torus->aabb;
    printf("aabb.min = %f, %f, %f\n", aabb.min[0], aabb.min[1], aabb.min[2]);
    printf("aabb.max = %f, %f, %f\n", aabb.max[0], aabb.max[1], aabb.max[2]);
}

bool getPresetColor(vec3 r, const char buffer[])
{
    if(strcmp(buffer, "RED") == 0)
    {
        vec3_copy(r, RED);
        return true;
    }else if(strcmp(buffer, "GREEN") == 0)
    {
        vec3_copy(r, GREEN);
        return true;
    }else if(strcmp(buffer, "BLUE") == 0)
    {
        vec3_copy(r, BLUE);
        return true;
    }else if(strcmp(buffer, "WHITE") == 0)
    {
        vec3_copy(r, WHITE);
        return true;
    }else if(strcmp(buffer, "BLACK") == 0)
    {
        vec3_copy(r, BLACK);
        return true;
    }else if(strcmp(buffer, "YELLOW") == 0)
    {
        vec3_copy(r, YELLOW);
        return true;
    }else if(strcmp(buffer, "CYAN") == 0)
    {
        vec3_copy(r, CYAN);
        return true;
    }else if(strcmp(buffer, "PINK") == 0)
    {
        vec3_copy(r, PINK);
        return true;
    }else if(strcmp(buffer, "GREY") == 0)
    {
        vec3_copy(r, GREY);
        return true;
    }
    return false;
}

bool parseVec3(vec3 r, FILE* fp)
{
    char buffer[128];
    if(!getNextTokenInFile(buffer, fp)){return false;}
    r[0] = (float)atof(buffer);
    if(!getNextTokenInFile(buffer, fp)){return false;}
    r[1] = (float)atof(buffer);
    if(!getNextTokenInFile(buffer, fp)){return false;}
    r[2] = (float)atof(buffer);
    return true;
}

bool parseColor(vec3 r, FILE* fp)
{
    char buffer[128];
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(!getPresetColor(r, buffer))
    {
        r[0] = (float)atof(buffer);
        if(!getNextTokenInFile(buffer, fp)){return false;}
        r[1] = (float)atof(buffer);
        if(!getNextTokenInFile(buffer, fp)){return false;}
        r[2] = (float)atof(buffer);
    }
    return true;
}

bool parseMatEntry(Material* mat, FILE* fp)
{
    char buffer[128];
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "PHONG") == 0)
    {
        mat->mat_type = PHONG;
        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SHADOWED
        if(!getNextTokenInFile(buffer, fp)){return false;}            
        if(strcmp(buffer, "yes") == 0)
        {
            mat->shadow = true;
        }else
        {
            mat->shadow = false;
        }
        
        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word DIFF_COLOR
        if(!parseColor(mat->cd, fp));

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SPEC_COLOR
        if(!parseColor(mat->cs, fp));        

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word AMB_COLOR
        if(!parseColor(mat->ca, fp));        

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word DIFF_CONSTANT
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->kd = (float)atof(buffer);

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SPEC_CONSTANT
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->ks = (float)atof(buffer);

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word AMB_CONSTANT
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->ka = (float)atof(buffer);
        return true;
    }else if(strcmp(buffer, "MATTE") == 0)
    {
        mat->mat_type = MATTE;
        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SHADOWED
        if(!getNextTokenInFile(buffer, fp)){return false;}            
        if(strcmp(buffer, "yes") == 0)
        {
            mat->shadow = true;
        }else
        {
            mat->shadow = false;
        }
        
        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word DIFF_COLOR
        if(!parseColor(mat->cd, fp));        

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word AMB_COLOR
        if(!parseColor(mat->ca, fp));                

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word DIFF_CONSTANT
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->kd = (float)atof(buffer);

        if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word AMB_CONSTANT
        if(!getNextTokenInFile(buffer, fp)){return false;}
        mat->ka = (float)atof(buffer);
        return true;        
        return true;
    }else
    {
        fprintf(stderr, "Invalid material type %s.\n", buffer);
        return false;
    }    
}

bool parseSphereEntry(Sphere** r, FILE* fp)
{
    char buffer[128];
    Sphere* sphere_ptr = (Sphere*)malloc(sizeof(Sphere));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        sphere_ptr->shadow = true;
    }else
    {
        sphere_ptr->shadow = false;
    }
    
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word RADIUS
    if(!getNextTokenInFile(buffer, fp)){return false;}
    sphere_ptr->radius = (float)atof(buffer);
    
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CENTER
    if(!parseVec3(sphere_ptr->center, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word PHI
    if(!getNextTokenInFile(buffer, fp)){return false;}
    sphere_ptr->phi = (float)atof(buffer);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MIN_THETA
    if(!getNextTokenInFile(buffer, fp)){return false;}
    sphere_ptr->min_theta = (float)atof(buffer);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MAX_THETA
    if(!getNextTokenInFile(buffer, fp)){return false;}
    sphere_ptr->max_theta = (float)atof(buffer);                

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!parseMatEntry(&(sphere_ptr->mat), fp)){return false;}

    *r = sphere_ptr;    
    return true;    
}

bool parsePlaneEntry(Plane** r, FILE* fp)
{
    char buffer[128];
    Plane* plane_ptr = (Plane*)malloc(sizeof(Plane));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        plane_ptr->shadow = true;
    }else 
    {
        plane_ptr->shadow = false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word POINT
    if(!parseVec3(plane_ptr->point, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word NORMAL
    if(!parseVec3(plane_ptr->normal, fp)){return false;}
    
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!parseMatEntry(&(plane_ptr->mat), fp)){return false;}

    *r = plane_ptr;
    return true;
}

bool parseRectEntry(Rectangle** r, FILE* fp)
{
    char buffer[128];
    Rectangle* rect_ptr = (Rectangle*)malloc(sizeof(Rectangle));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        rect_ptr->shadow = true;
    }else
    {
        rect_ptr->shadow = false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word POINT
    if(!parseVec3(rect_ptr->point, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word WIDTH
    if(!parseVec3(rect_ptr->width, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word HEIGHT
    if(!parseVec3(rect_ptr->height, fp)){return false;}

    vec3_cross(rect_ptr->normal, rect_ptr->width, rect_ptr->height);
    vec3_normalize(rect_ptr->normal, rect_ptr->normal);    

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!parseMatEntry(&(rect_ptr->mat), fp)){return false;}

    *r = rect_ptr;
    return true;
}

bool parseTriangleEntry(Triangle** r,  FILE* fp)
{
    char buffer[128];
    Triangle* tri_ptr = (Triangle*)malloc(sizeof(Triangle));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        tri_ptr->shadow = true;
    }else
    {
        tri_ptr->shadow = false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word V0
    if(!parseVec3(tri_ptr->v0, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word V1
    if(!parseVec3(tri_ptr->v1, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word V2
    if(!parseVec3(tri_ptr->v2, fp)){return false;}

    calcTriangleNormal(tri_ptr);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!parseMatEntry(&(tri_ptr->mat), fp)){return false;}

    *r = tri_ptr;
    return true;
}

bool parseAABoxEntry(AABox** r,  FILE* fp)
{
    char buffer[128];
    AABox* aabox_ptr = (AABox*)malloc(sizeof(AABox));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        aabox_ptr->shadow = true;
    }else
    {
        aabox_ptr->shadow = false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MIN
    if(!parseVec3(aabox_ptr->min, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MAX
    if(!parseVec3(aabox_ptr->max, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!parseMatEntry(&(aabox_ptr->mat), fp)){return false;}

    *r = aabox_ptr;
    return true;
}

bool parseOpenCylEntry(OpenCylinder** r,  FILE* fp)
{
    char buffer[128];
    OpenCylinder* cyl_ptr = (OpenCylinder*)malloc(sizeof(OpenCylinder));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        cyl_ptr->shadow = true;
    }else
    {
        cyl_ptr->shadow = false;
    }
    
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word HALF_HEIGHT
    if(!getNextTokenInFile(buffer, fp)){return false;}
    cyl_ptr->half_height = (float)atof(buffer);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word RADIUS
    if(!getNextTokenInFile(buffer, fp)){return false;}
    cyl_ptr->radius = (float)atof(buffer);    
    
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word PHI
    if(!getNextTokenInFile(buffer, fp)){return false;}
    cyl_ptr->phi = (float)atof(buffer);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word NORMAL_TYPE
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "OPEN") == 0)
    {
        cyl_ptr->normal_type = OPEN;
    }else if(strcmp(buffer, "CONVEX") == 0)
    {
        cyl_ptr->normal_type = CONVEX;
    }else if(strcmp(buffer, "CONCAVE") == 0)
    {
        cyl_ptr->normal_type = CONCAVE;
    }else
    {
        fprintf(stderr, "Invalid normal type %s\n", buffer);
        return false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!parseMatEntry(&(cyl_ptr->mat), fp)){return false;}    

    *r = cyl_ptr;
    return true;
}

bool parseDiskEntry(Disk** r,  FILE* fp)
{
    char buffer[128];
    Disk* disk_ptr = (Disk*)malloc(sizeof(Disk));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        disk_ptr->shadow = true;
    }else
    {
        disk_ptr->shadow = false;
    }

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MIN
    if(!parseVec3(disk_ptr->center, fp)){return false;}

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MAX
    if(!parseVec3(disk_ptr->normal, fp)){return false;}
    printVec3WithText("disk->normal", disk_ptr->normal);
    vec3_normalize(disk_ptr->normal, disk_ptr->normal);
    printVec3WithText("disk->normal normalized ", disk_ptr->normal);    

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word RADIUS
    if(!getNextTokenInFile(buffer, fp)){return false;}
    disk_ptr->radius = (float)atof(buffer);    

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!parseMatEntry(&(disk_ptr->mat), fp)){return false;}

    *r = disk_ptr;
    return true;
}

bool parseTorusEntry(Torus** r,  FILE* fp)
{
    char buffer[128];
    Torus* torus_ptr = (Torus*)malloc(sizeof(Torus));
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word CAST_SHADOW
    if(!getNextTokenInFile(buffer, fp)){return false;}
    if(strcmp(buffer, "yes") == 0)
    {
        torus_ptr->shadow = true;
    }else
    {
        torus_ptr->shadow = false;
    }
    
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word SWEPT_RADIUS
    if(!getNextTokenInFile(buffer, fp)){return false;}
    torus_ptr->swept_radius = (float)atof(buffer);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word TUBE_RADIUS
    if(!getNextTokenInFile(buffer, fp)){return false;}
    torus_ptr->tube_radius = (float)atof(buffer);    
    
    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word PHI
    if(!getNextTokenInFile(buffer, fp)){return false;}
    torus_ptr->phi = (float)atof(buffer);

    calcAABBTorus(torus_ptr);

    if(!getNextTokenInFile(buffer, fp)){return false;}    // Skip over the word MATERIAL
    if(!parseMatEntry(&(torus_ptr->mat), fp)){return false;}    

    //printTorus(torus_ptr);
    *r = torus_ptr;
    return true;
}
