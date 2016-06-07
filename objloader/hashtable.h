#pragma once

#include <cstdint>
#include <cmath>
#include "dbuffer.h"

typedef struct
{
    uint64_t key;
    int index;
}HashNode;

typedef struct
{
    DBuffer table;
    int size;
    int used;
}HashTable;

#define HashTable_array(a) (HashNode*)((a).table.data);

HashTable HashTable_create(const int size)
{
    HashTable ht;
    ht.table = DBuffer_create_cap(HashNode, size);
    HashNode* table_ptr = HashTable_array(ht);
    for(int i = 0; i < size; i++)
    {
        table_ptr[i].key = 0;
        table_ptr[i].index = -1;        
    }
    ht.size = size;
    ht.used = 0;
    return ht;
}

void HashTable_destroy(HashTable* ht)
{
    DBuffer_destroy(&(ht->table));
    ht->size = 0;
    ht->used = 0;
}

static const double HASH_CONSTANT = 0.6180339887f;
const float LOAD_FACTOR = 0.7f;

int hashFunc(const uint64_t key, const int table_size)
{
    double tmp = key * HASH_CONSTANT;
    double fractional = tmp - floor(tmp);
    int hash_index = (int)floor(fractional * table_size);
    return hash_index;
}

int getTableIndex(HashTable* index_table, const uint64_t key)
{
    int hash_index = hashFunc(key, index_table->size);
    HashNode* table_ptr = HashTable_array(*index_table);
    for(int i = hash_index; table_ptr[i].index != -1 && i < index_table->size; i++)
    {
        if(table_ptr[i].key == key)
        {
            return table_ptr[i].index;
        }
    }
    return -1;
}

void rehash(HashTable* ht);

int insertTable(HashTable* index_table, const uint64_t key, const int index)
{
    float load = (float)(index_table->used + 1) / (float)index_table->size;
    if(load >= 0.7f)
    {
        rehash(index_table);
    }
    int hash_index = hashFunc(key, index_table->size);    
    HashNode* table_ptr = HashTable_array(*index_table);
    for(int i = hash_index; i < index_table->size; i++)
    {
        if(table_ptr[i].index == -1)
        {
            table_ptr[i].key = key;
            table_ptr[i].index = index;
            index_table->used += 1;
            return index;
        }
    }
    fprintf(stderr, "Cannot insert into table.\n");
    return -1;
}

void rehash(HashTable* ht)
{
    HashTable new_ht = HashTable_create(ht->size * 2);
    HashNode* table_ptr = HashTable_array(*ht);
    for(int i = 0, count = 0; count < ht->used && i < ht->size; i++)
    {
        if(table_ptr[i].index != -1)
        {
            insertTable(&new_ht, table_ptr[i].key, table_ptr[i].index);
        }
    }

    HashTable_destroy(ht);
    *ht = new_ht;
}
