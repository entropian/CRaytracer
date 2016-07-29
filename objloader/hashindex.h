#pragma once
#include <cmath>
#include <cassert>
#include <cstdlib>

const int DEFAULT_HASH_SIZE = 1024;
const int DEFAULT_INDEX_SIZE = 1024;
const int DEFAULT_HASH_GRANULARITY = 1024;
static int INVALID_INDEX[1] = {-1};

typedef struct HashIndex_s
{
    int *hash;
    int *index_chain;
    int hash_size;
    int index_size;
    int granularity;
    int hash_mask;
    int lookup_mask;    
}HashIndex;

inline void HashIndex_init(HashIndex* hash_index, const int initial_hash_size, const int initial_index_size)
{
    assert(isPowerOfTwo(initial_hash_size));    
    hash_index->hash_size = initial_hash_size;
    hash_index->index_size = initial_index_size;
    hash_index->hash = INVALID_INDEX;
    hash_index->index_chain = INVALID_INDEX;
    hash_index->granularity = DEFAULT_HASH_GRANULARITY;
    hash_index->hash_mask = hash_index->hash_size - 1;
    hash_index->lookup_mask = 0;
}

inline void HashIndex_free(HashIndex* hash_index)
{
    if(hash_index->hash != INVALID_INDEX)
    {
        free(hash_index->hash);
        hash_index->hash = INVALID_INDEX;
    }
    if(hash_index->index_chain != INVALID_INDEX)
    {
        free(hash_index->index_chain);
        hash_index->index_chain = INVALID_INDEX;
    }
    hash_index->lookup_mask = 0;
}

bool isPowerOfTwo(int x)
{
    return (x & (x - 1)) == 0 && x > 0;
}

void HashIndex_alloc(HashIndex* hash_index, const int new_hash_size, const int new_index_size)
{
    assert(isPowerOfTwo(new_hash_size));
    HashIndex_free(hash_index);
    hash_index->hash_size = new_hash_size;
    hash_index->hash = (int*)malloc(sizeof(int) * new_hash_size);
    memset(hash_index->hash, 0xff, new_hash_size * sizeof(hash_index->hash[0]));
    hash_index->index_size = new_index_size;
    hash_index->index_chain = (int*)malloc(sizeof(int) * new_index_size);
    memset(hash_index->index_chain, 0xff, new_index_size * sizeof(hash_index->index_chain[0]));
    hash_index->hash_mask = new_hash_size - 1;
    hash_index->lookup_mask = -1;
}


void HashIndex_resize_index(HashIndex* hash_index, const int new_index_size)
{
    if(new_index_size <= hash_index->index_size)
    {
        return;
    }
    
    int mod = new_index_size % hash_index->granularity;
    int new_size;
    if(!mod)
    {
        new_size = new_index_size;
    }else
    {
        new_size = new_index_size + hash_index->granularity - mod;
    }

    int* old_chain = hash_index->index_chain;
    hash_index->index_chain = (int*)malloc(sizeof(int) * new_size);
    int index_size = hash_index->index_size;
    memcpy(hash_index->index_chain, old_chain, index_size * sizeof(int));
    memset(hash_index->index_chain + index_size, 0xff, (new_size - index_size) * sizeof(int));
    free(old_chain);
    hash_index->index_size = new_size;
}

inline void HashIndex_add(HashIndex* hash_index, const int key, const int index)
{
    assert(index >= 0);
    int index_size = hash_index->index_size;
    if(hash_index->hash == INVALID_INDEX)
    {
        HashIndex_alloc(hash_index, hash_index->hash_size, index >= index_size ? index + 1 : index_size);        
    }else if(index >= index_size)
    {
        HashIndex_resize_index(hash_index, index + 1);
    }
    int h = key & hash_index->hash_mask;
    hash_index->index_chain[index] = hash_index->hash[h];
    hash_index->hash[h] = index;
}

inline void HashIndex_remove(HashIndex* hash_index, const int key, const int index)
{
    int k = key & hash_index->hash_mask;
    if(hash_index->hash == INVALID_INDEX)
    {
        return;
    }
    if(hash_index->hash[k] == index)
    {
        hash_index->hash[k] = hash_index->index_chain[index];
    }else
    {
        for(int i = hash_index->hash[k]; i != -1; i = hash_index->index_chain[i])
        {
            if(hash_index->index_chain[i] == index)
            {
                hash_index->index_chain[i] = hash_index->index_chain[index];
                break;
            }
        }
    }
    hash_index->index_chain[index] = -1;
}

inline int HashIndex_hash(const char* string)
{
    int i, hash = 0;
    for(i = 0; *string != '\0'; i++)
    {
        hash += (*string++) * (i + 119);
    }
    return hash;
}

inline int HashIndex_gen_key(const HashIndex* hash_index, const char *string)
{
    return HashIndex_hash(string) & hash_index->hash_mask;
}

inline int HashIndex_First(const HashIndex* hash_index, const int key)
{
    int hash_mask = hash_index->hash_mask;
    int lookup_mask = hash_index->lookup_mask;
    return hash_index->hash[key & hash_mask & lookup_mask];
}

inline int HashIndex_Next(const HashIndex* hash_index, const int index)
{
    assert(index >= 0 && index < hash_index->index_size);
    int lookup_mask = hash_index->lookup_mask;
    return hash_index->index_chain[index & lookup_mask];
}




