#pragma once
#include <stdio.h>
#include "util/util.h"

// EXR
#include "/usr/include/OpenEXR/ImfRgba.h"
#include "/usr/include/OpenEXR/ImfRgbaFile.h"
#include "/usr/include/OpenEXR/ImfArray.h"
#include "/usr/include/OpenEXR/ImathBox.h"

static void readRgba1(const char file_name[], float** data, int* width, int* height)
{
    Imf::RgbaInputFile file(file_name);
    //Imath_2_2::Box2i dw = file.dataWindow();
    IMATH_NAMESPACE::Box2i dw = file.dataWindow();

    *width = dw.max.x - dw.min.x + 1;
    *height = dw.max.y - dw.min.y + 1;
    Imf::Rgba* pixels;
    pixels = (Imf::Rgba*)malloc(sizeof(Imf::Rgba) * *width * *height);

    file.setFrameBuffer(&pixels[0] - dw.min.x - dw.min.y * *width, 1, *width);
    file.readPixels(dw.min.y, dw.max.y);

    int num_pixels = *width * *height;
    float* floats = (float*)malloc(sizeof(float) * num_pixels * 3);
    for(int i = 0; i < num_pixels; i++)
    {
        floats[i*3] = pixels[i].r;
        floats[i*3 + 1] = pixels[i].g;
        floats[i*3 + 2] = pixels[i].b;
    }
    *data = floats;
}

static int PPM_write(const char *file_name, const unsigned char *image, const int image_size, 
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

static int PPM_read(unsigned char** image, int* image_size, int* width, int* height, const char* file_name)
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
