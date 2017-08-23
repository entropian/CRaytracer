#pragma once
#include <stdio.h>
#include "util/util.h"

int PPM_write(const char *file_name, const unsigned char *image, const int image_size, 
			   const int width, const int height)
{
	int num_pixel = width * height;	
	if (image_size / 3 != num_pixel)
	{
		fprintf(stderr, "Image size doesn't match image dimensions.\n");
		return 0;
	}

	FILE *fp;
    openFile(&fp, file_name, "wb");

	// Header
	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", width, height);
	fprintf(fp, "255\n");
	
	fwrite(image, sizeof(char), image_size, fp);
	fclose(fp);
	return 1;
}

int PPM_read(unsigned char** image, int* image_size, int* width, int* height, const char* file_name)
{
	FILE *fp;
    openFile(&fp, file_name, "r");

    char buffer[16];
    int image_width, image_height, tmp;
    fscanf(fp, "%s\n", buffer);
    fscanf(fp, "%d %d\n", &image_width, &image_height);
    fscanf(fp, "%d\n", &tmp);
    int size = image_width * image_height * 3;
    *image = (unsigned char*)malloc(size * sizeof(unsigned char));
    int result = fread(*image, sizeof(char), size, fp);
    if(result != size)
    {
        fprintf(stderr, "Read length doesn't match image size\n");
    }
    *image_size = result;
    *width = image_width;
    *height = image_height;
    fclose(fp);
    return 1;
}
