#pragma once
#include <stdio.h>

int PPM_write(const char *file_name, const unsigned char *image, const int image_size, 
			   const int width, const int height)
{
	int num_pixel = width * height;	
	if (image_size / 3 != num_pixel)
	{
		fprintf(stderr, "Image size doesn't match image dimensions.\n");
		return -1;
	}

	FILE *fp;
#ifdef _MSC_VER
	fopen_s(&fp, file_name, "wb");
#else
  fp = fopen(file_name, "wb");
#endif

	// Header
	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", width, height);
	fprintf(fp, "255\n");
	
	fwrite(image, sizeof(char), image_size, fp);
	fclose(fp);
	return 0;
}
