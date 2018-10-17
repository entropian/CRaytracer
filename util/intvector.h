#pragma once

#include <cstdlib>
#include <cstdio>

static const int INITIAL_CAPACITY = 10;

typedef struct IntVector
{
    int* array;
    int size;
    int max;
} IntVector;

IntVector IntVector_create()
{
    IntVector r;
    //r.array = (int*)malloc(sizeof(int) * INITIAL_CAPACITY);
    r.array = NULL;
    //r.max = INITIAL_CAPACITY;
    r.max = 0;
    r.size = 0;
    return r;
}

void IntVector_destroy(IntVector* iv)
{
    free(iv->array);
    iv->max = 0;
    iv->size = 0;
}

void IntVector_push(IntVector* iv, const int element)
{
    if(iv->array == NULL)
    {
        iv->array = (int*)malloc(sizeof(int) * INITIAL_CAPACITY);
        iv->max = INITIAL_CAPACITY;
    }
    
    if(iv->size < iv->max)
    {
        iv->array[iv->size] = element;
        iv->size += 1;
    }else
    {
        int* p = (int*)realloc((void*)iv->array, iv->max * 2 * sizeof(int));
        if(p)
        {
            iv->array = p;
            iv->max *= 2;
            iv->array[iv->size] = element;
            iv->size += 1;
        }else
        {
            fprintf(stderr, "Can't reallocate memory for IntVector.\n");
        }
    }
}

int IntVector_get(const IntVector* iv, const int index)
{
    if(index >= iv->size)
    {
        fprintf(stderr, "Invalid index, returning -1.\n");
        return -1;
    }
    return iv->array[index];
}

void IntVector_set(IntVector* iv, const int index, const int element)
{
    if(index >= iv->max)
    {
        fprintf(stderr, "Invalid index.\n");
        return;
    }
    iv->array[index] = element;
}

void IntVector_clear(IntVector* iv, const int val)
{
    for(int i = 0; i < iv->max; i++)
    {
        iv->array[0] = val;
    }
}

