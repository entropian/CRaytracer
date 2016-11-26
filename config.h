#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util/util.h"
#include "scene/scenedata.h"
#include "trace.h"

typedef struct PhotonmapConfig_s
{
    int num_photons;
    int num_caustic_photons;
    int photon_depth;
    int num_estimate;    
    float photon_radius;
    float caustic_radius;
}PhotonmapConfig;

typedef struct ConfigParams_s
{    
    unsigned int window_width;
    unsigned int window_height;
    unsigned int image_width;
    unsigned int image_height;
    unsigned int num_samples;
    unsigned int num_sample_sets;
    unsigned int max_depth;
    TraceType trace_type;
    AccelType accel_type;
    bool image_save;
    bool photon_map;
    PhotonmapConfig pm_config;
    char file_name[128];
}ConfigParams;

void parseConfigFile(ConfigParams* cp)
{
    FILE *fp;
#ifdef _MSC_VER
    fopen_s(&fp, "config.txt", "r");
#else
    fp = fopen("config.txt", "r");
#endif
    char buffer[128];
    if(!fp)
    {
        fprintf(stderr, "config.txt cannot be opened.\n");
    }
    
    while(getNextTokenInFile(buffer, fp))
    {
        if(strcmp(buffer, "scene_file") == 0)
        {
            getNextTokenInFile(buffer, fp);
            //strcpy_s(cp->file_name, NAME_LENGTH, buffer);
            stringCopy(cp->file_name, NAME_LENGTH, buffer);
        }else if(strcmp(buffer, "window_width") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->window_width = atoi(buffer);
        }else if(strcmp(buffer, "window_height") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->window_height = atoi(buffer);
        }else if(strcmp(buffer, "image_width") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->image_width = atoi(buffer);
        }else if(strcmp(buffer, "image_height") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->image_height = atoi(buffer);
        }else if(strcmp(buffer, "num_samples") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->num_samples = atoi(buffer);
        }else if(strcmp(buffer, "num_sample_sets") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->num_sample_sets = atoi(buffer);
        }else if(strcmp(buffer, "max_depth") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->max_depth = atoi(buffer);
        }else if(strcmp(buffer, "trace_type") == 0)
        {
            getNextTokenInFile(buffer, fp);
            if(strcmp(buffer, "RAYCAST") == 0)
            {
                cp->trace_type = RAYCAST;
            }else if(strcmp(buffer, "WHITTED") == 0)
            {
                cp->trace_type = WHITTED;
            }else if(strcmp(buffer, "PATHTRACE") == 0)
            {
                cp->trace_type = PATHTRACE;
            }
        }else if(strcmp(buffer, "accel_struct") == 0)
        {
            getNextTokenInFile(buffer, fp);
            if(strcmp(buffer, "BVH") == 0)
            {
                cp->accel_type = BVH;
            }else if(strcmp(buffer, "GRID") == 0)
            {
                cp->accel_type = GRID;
            }else if(strcmp(buffer, "BVH4") == 0)
            {
                cp->accel_type = BVH4;
            }else
            {
                cp->accel_type = NONE;
            }
        }else if(strcmp(buffer, "image_save") == 0)
        {
            cp->image_save = false;
            getNextTokenInFile(buffer, fp);
            if(strcmp(buffer, "yes") == 0)
            {
                cp->image_save = true;
            }
        }else if(strcmp(buffer, "photon_map") == 0)
        {
            cp->photon_map = false;
            getNextTokenInFile(buffer, fp);
            if(strcmp(buffer, "yes") == 0)
            {
                cp->photon_map = true;
            }
        }else if(strcmp(buffer, "num_photons") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->pm_config.num_photons = atoi(buffer);
        }else if(strcmp(buffer, "num_caustic_photons") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->pm_config.num_caustic_photons = atoi(buffer);
        }else if(strcmp(buffer, "photon_depth") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->pm_config.photon_depth = atoi(buffer);
        }else if(strcmp(buffer, "photon_radius") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->pm_config.photon_radius = atof(buffer);
        }else if(strcmp(buffer, "caustic_radius") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->pm_config.caustic_radius = atof(buffer);
        }else if(strcmp(buffer, "num_estimate") == 0)
        {
            getNextTokenInFile(buffer, fp);
            cp->pm_config.num_estimate = atoi(buffer);
        }        
    }
    fclose(fp);
}
