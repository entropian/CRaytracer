#pragma once
#include <stdio.h>

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
#ifdef _MSC_VER
	fopen_s(&fp, file_name, "wb");
#else
    fp = fopen(file_name, "wb");
#endif
    if(!fp)
    {
        fprintf(stderr, "Failed to open file %s\n", file_name);
        return 0;
    }

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
#ifdef _MSC_VER
	fopen_s(&fp, file_name, "r");
#else
    fp = fopen(file_name, "r");
#endif
    if(!fp)
    {
        fprintf(stderr, "Failed to open file %s\n", file_name);
        return 0;
    }

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
