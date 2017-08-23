#include <stdio.h>
#include <string.h>


int saveImageState(float* color_buffer, int num_samples, int width, int height, const char* file_name)
{
    FILE *fp;
    fprintf(fp, "%d\n", num_samples);
    fprintf(fp, "%d %d\n", width, height);
    int size = width * height * 3;
    int result = fwrite(color_buffer, sizeof(float), size, fp);
    if(result != size)
    {
        fprintf(stderr, "Write error.\n");
    }
    fclose(fp);
    return 1;
}

int readImageState(float** color_buffer, int* num_samples, int* size, int* width, int* height, const char* file_name)
{
	FILE *fp;
    openFile(&fp, file_name, "r");
    int buffer_num_samples, buffer_size, buffer_width, buffer_height;
    fscanf(fp, "%d\n", &buffer_num_samples);
    fscanf(fp, "%d %d\n", &buffer_width, &buffer_height);
    buffer_size = buffer_width * buffer_height * 3;
    float* buffer = (float*)malloc(buffer_size * sizeof(float));
    int result = fread(buffer, sizeof(float), buffer_size, fp);
    printf("result %d\n", result);
    printf("buffer_size %d\n", buffer_size);
    if(result != buffer_size)
    {
        fprintf(stderr, "Read error\n");
    }
    *color_buffer = buffer;
    *num_samples = buffer_num_samples;
    *size = result;
    *width = buffer_width;
    *height = buffer_height;
    fclose(fp);
}
