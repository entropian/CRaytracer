#pragma once

#include <cstdio>
#include <cstdlib>

static const int DEFAULT_INITIAL_CAPACITY = 100;

typedef struct
{
    char* data;
    int max;            // Max capacity in bytes
    int size;           // Size in bytes
    int element_size;
}DBuffer;

#define DBuffer_size(a) ((a).size) / ((a).element_size)
#define DBuffer_max_elements(a) ((a).max) / ((a).element_size)

void DBuffer_push_f(DBuffer* dbuf, const char* new_element, const size_t element_size);
#define DBuffer_push(a, b) DBuffer_push_f((&(a)), (const char*)(&b), sizeof(b))

DBuffer DBuffer_create_f(const int element_size);
#define DBuffer_create(a) DBuffer_create_f(sizeof(a))

DBuffer DBuffer_create_cap_f(const int element_size, const int max_elements);
#define DBuffer_create_cap(a, b) DBuffer_create_cap_f(sizeof(a), b)

void DBuffer_erase(DBuffer* dbuf);
void DBuffer_assume(DBuffer* dbuf, char* data_ptr, const int size, const int max, const size_t element_size);
void* DBuffer_data_ptr_f(DBuffer* dbuf);
#define DBuffer_data_ptr(a) DBuffer_data_ptr_f((&(a)))
void DBuffer_destroy(DBuffer* dbuf);

#ifdef OBJ_LOADER_IMPLEMENTATION    
void DBuffer_push_f(DBuffer* dbuf, const char* new_element, const size_t element_size)
{
    if(dbuf->size + element_size > (unsigned)dbuf->max)
    {
        int new_max = dbuf->max * 2;
        char* new_data = (char*)realloc(dbuf->data, new_max);
        if(new_data)
        {
            dbuf->data = new_data;
            dbuf->max = new_max;
        }else
        {
            fprintf(stderr, "DBuffer out of memory.\n");
            return;
        }
    }
    for(size_t i = 0; i < element_size; i++)
    {
        dbuf->data[dbuf->size] = new_element[i];
        dbuf->size += 1;
    }
}


DBuffer DBuffer_create_f(const int element_size)
{
    DBuffer dbuf;
    dbuf.data = (char*)malloc(element_size * DEFAULT_INITIAL_CAPACITY);
    dbuf.max = element_size * DEFAULT_INITIAL_CAPACITY;
    dbuf.size = 0;
    dbuf.element_size = element_size;
    return dbuf;
}

DBuffer DBuffer_create_cap_f(const int element_size, const int max_elements)
{
    DBuffer dbuf;
    dbuf.data = (char*)malloc(element_size * max_elements);
    dbuf.max = element_size * max_elements;
    dbuf.size = 0;
    dbuf.element_size = element_size;
    return dbuf;    
}

void DBuffer_erase(DBuffer* dbuf)
{
    for(int i = 0; i < dbuf->size; ++i)
    {
        dbuf->data[i] = 0;
    }
    dbuf->size = 0;
}

void DBuffer_assume(DBuffer* dbuf, char* data_ptr, const int size, const int max, const size_t element_size)
{
    dbuf->data = data_ptr;
    dbuf->size = size * element_size;
    dbuf->max = max * element_size;
    dbuf->element_size = element_size;
}

void* DBuffer_data_ptr_f(DBuffer *dbuf)
{
    return (void*)(dbuf->data);
}

void DBuffer_destroy(DBuffer* dbuf)
{
    if(dbuf->data){free(dbuf->data);}
    dbuf->max = 0;
    dbuf->size = 0;
    dbuf->element_size = 0;
}
#endif
