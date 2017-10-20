#include "mempool.h"
#include <stdio.h>
#include <stdlib.h>

int MemPool_init(MemPool* mem_pool, int element_size, int num_elements)
{    
    int block_size = element_size * num_elements;
    mem_pool->block = (unsigned char*)malloc(block_size);
    if(!mem_pool->block)
    {
        fprintf(stderr, "Cannot allocate memory for MemPool.\n");
        return 0;
    }
    mem_pool->block_size = block_size;
    mem_pool->element_size = element_size;
    mem_pool->num_elements = num_elements;

    int i;
    unsigned char* prev_ptr = mem_pool->block;
    unsigned char* cur_ptr = prev_ptr;
    for(i = 0; i < num_elements - 1; i++)
    {
        cur_ptr += element_size;
        unsigned char** address_ptr = (unsigned char**)prev_ptr;
        *address_ptr = cur_ptr;
        prev_ptr = cur_ptr;
    }
    unsigned char** address_ptr = (unsigned char**)cur_ptr;
    *address_ptr = NULL;

    mem_pool->head = mem_pool->block;
    mem_pool->tail = cur_ptr;

    return 1;
}

void MemPool_destroy(MemPool* mem_pool)
{
    free(mem_pool->block);
    mem_pool->block_size = 0;
    mem_pool->element_size = 0;
    mem_pool->num_elements = 0;
    mem_pool->head = NULL;
    mem_pool->tail = NULL;
}

unsigned char* MemPool_requestElement(MemPool* mem_pool)
{
    unsigned char* r = mem_pool->head;
    if(!r)
    {
        fprintf(stderr, "MemPool out of memory.\n");
        return NULL;
    }
    unsigned char** address_ptr = (unsigned char**)(mem_pool->head);
    mem_pool->head = *address_ptr;
    return r;
}

void MemPool_releaseElement_f(MemPool* mem_pool, unsigned char** element)
{
    if(*element < mem_pool->block || *element > mem_pool->block + mem_pool->block_size)
    {
        fprintf(stderr, "Invalid element pointer\n.");
        return;
    }
    unsigned char** address_ptr = (unsigned char**)(mem_pool->tail);
    *address_ptr = (unsigned char*)(*element);
    address_ptr = (unsigned char**)(*element);
    *address_ptr = NULL;
    mem_pool->tail = *element;
    *element = NULL;
    if(!mem_pool->head)
    {
        mem_pool->head = mem_pool->tail;
    }
}

