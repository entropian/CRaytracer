#pragma once

typedef struct MemPool_s
{
    unsigned char* block;
    unsigned char* head;
    unsigned char* tail;
    int block_size;
    int element_size;
    int num_elements;    
}MemPool;

int MemPool_init(MemPool* mem_pool, int element_size, int num_elements);
void MemPool_destroy(MemPool* mem_pool);
unsigned char* MemPool_requestElement(MemPool* mem_pool);

void MemPool_releaseElement_f(MemPool* mem_pool, unsigned char** element);
#define MemPool_releaseElement(a, b) MemPool_releaseElement_f((a), ((unsigned char**)(b)))
