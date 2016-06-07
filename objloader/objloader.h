#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "dbuffer.h"
#include "hashtable.h"

#define NAME_LENGTH 30

typedef struct
{
    float* positions;
    float* normals;
    float* texcoords;
    int* indices;
    float num_positions, num_normals, num_texcoords, num_indices;
    char mat_name[NAME_LENGTH];
    char mesh_name[NAME_LENGTH];
}OBJShape;

void OBJShape_destroy(OBJShape* obj_shape)
{
    if(obj_shape->positions){free(obj_shape->positions);}
    if(obj_shape->normals){free(obj_shape->normals);}
    if(obj_shape->texcoords){free(obj_shape->texcoords);}
    if(obj_shape->indices){free(obj_shape->indices);}

    obj_shape->num_positions = 0;
    obj_shape->num_normals = 0;
    obj_shape->num_texcoords = 0;
    obj_shape->num_indices = 0;    
}

typedef struct
{
    int v_idx, vn_idx, vt_idx;
}VertexIndex;

bool OBJGetLine(FILE* fp, int* read_result, char buffer[], const int buffer_size)
{
    bool not_end_of_file = true;
    bool reading = true;
    int index = 0;
    while(reading)
    {
        char c = fgetc(fp);
        switch(c)
        {
        case '\n':
            reading = false;
            break;
        case '\r':
            continue;
        case EOF:
            reading = false;
            not_end_of_file = false;
            break;
        case '\0':
            reading = false;
            not_end_of_file = false;
            break;
        default:
            buffer[index++] = c;            
        }
    }
    *read_result = index;
    buffer[index] = '\0';
    return not_end_of_file;
}

static inline int OBJParseInt(const char** str)
{
    *str += strspn(*str, " \t");
    int r = atoi(*str);
    *str += strcspn(*str, " \t\r");
    return r;
}

static inline float OBJParseFloat(const char** str)
{
    *str += strspn(*str, " \t");
    float r = (float)atof(*str);
    *str += strcspn(*str, " \t\r");
    return r;
}

static inline void OBJParseFloat2(float *x, float *y, const char** str)
{
    *x = OBJParseFloat(str);
    *y = OBJParseFloat(str);
}

static inline void OBJParseFloat3(float *x, float *y, float *z, const char** str)
{
    *x = OBJParseFloat(str);
    *y = OBJParseFloat(str);
    *z = OBJParseFloat(str);    
}

VertexIndex OBJParseFaceTriple(const char** str)
{
    VertexIndex vi = {0, 0, 0};

    *str += strspn(*str, " \t");
    vi.v_idx = atoi(*str);
    *str += strcspn(*str, "/ \t");

    if((*str)[0] != '/') // Only positions
    {
        return vi;
    }
    (*str)++;

    if((*str)[0] == '/')    // No texcoord -- v//vn
    {
        *str++;
        vi.vn_idx = atoi(*str);
        return vi;
    }
    
    vi.vt_idx = atoi(*str);
    *str += strcspn(*str, "/ \t");
    if((*str)[0] != '/')    // No normal -- v/vt 
    {
        return vi;
    }
    (*str)++;

    vi.vn_idx = atoi(*str);
    *str += strcspn(*str, "/ \t\r\0");
    return vi;
}

// TODO: length problem
void OBJParseString(char buffer[], const char** str)
{
    *str += strspn(*str, " \t");
    int length = strcspn(*str, " \t\r\0");
    strncpy(buffer, *str, length);
    buffer[length] = '\0';
    *str += length;
}

uint64_t getKey(const VertexIndex* vi)
{
    uint64_t key;
    key = vi->v_idx;
    key *= 1000000;
    key += vi->vt_idx;
    key *= 1000000;
    key += vi->vn_idx;
    return key;
}

int UpdateVertexCache(DBuffer* positions, DBuffer* normals, DBuffer* texcoords, HashTable* index_table,
                           const VertexIndex* vi, const DBuffer* in_positions, const DBuffer* in_normals,
                           const DBuffer* in_texcoords)
{
    int index;
    if((index = getTableIndex(index_table, getKey(vi))) != -1)
    {
        return index;
    }

    assert(vi->v_idx * 3 < DBuffer_size(*in_positions));
    float* in_pos_ptr = (float*)(in_positions->data);
    DBuffer_push(*positions, in_pos_ptr[vi->v_idx * 3]);
    DBuffer_push(*positions, in_pos_ptr[vi->v_idx * 3 + 1]);
    DBuffer_push(*positions, in_pos_ptr[vi->v_idx * 3 + 2]);

    if(vi->vt_idx >= 0)
    {
        float* in_texcoords_ptr = (float*)(in_texcoords->data);
        DBuffer_push(*texcoords, in_texcoords_ptr[vi->vt_idx * 2]);
        DBuffer_push(*texcoords, in_texcoords_ptr[vi->vt_idx * 2 + 1]);
    }

    if(vi->vn_idx >= 0)
    {
        float* in_normal_ptr = (float*)(in_normals->data);
        DBuffer_push(*normals, in_normal_ptr[vi->vn_idx * 3]);
        DBuffer_push(*normals, in_normal_ptr[vi->vn_idx * 3 + 1]);
        DBuffer_push(*normals, in_normal_ptr[vi->vn_idx * 3 + 2]);
    }

    index = DBuffer_size(*positions) / 3 - 1;
    insertTable(index_table, getKey(vi), index);
    return index;
}

bool exportGroupToShape(OBJShape* shape, const DBuffer* in_positions, const DBuffer* in_normals,
                        const DBuffer* in_texcoords, DBuffer* in_face_group)
{
    if(DBuffer_size(*in_face_group) == 0)
    {
        return false;
    }
    OBJShape new_shape;
    DBuffer positions = DBuffer_create(float);
    DBuffer normals = DBuffer_create(float);
    DBuffer texcoords = DBuffer_create(float);
    DBuffer indices = DBuffer_create(int);

    HashTable index_table = HashTable_create(1000);

    DBuffer* face_group_ptr = (DBuffer*)(in_face_group->data);
    for(int i = 0; i < DBuffer_size(*in_face_group); i++)
    {
        VertexIndex* vi_ptr = (VertexIndex*)(face_group_ptr[i].data);
        VertexIndex* v0 = &(vi_ptr[0]);
        VertexIndex* v1;
        VertexIndex* v2 = &(vi_ptr[1]);
        for(int j = 2; j < DBuffer_size(face_group_ptr[i]); j++)
        {
            v1 = v2;
            v2++;

            int i0 = UpdateVertexCache(&positions, &normals, &texcoords, &index_table,
                                       v0, in_positions, in_normals, in_texcoords);
            int i1 = UpdateVertexCache(&positions, &normals, &texcoords, &index_table,
                                       v1, in_positions, in_normals, in_texcoords);
            int i2 = UpdateVertexCache(&positions, &normals, &texcoords, &index_table,
                                       v2, in_positions, in_normals, in_texcoords);


            DBuffer_push(indices, i0);
            DBuffer_push(indices, i1);
            DBuffer_push(indices, i2);                    
        }
    }

    new_shape.positions = (float*)(positions.data);
    new_shape.normals = (float*)(normals.data);
    new_shape.texcoords = (float*)(texcoords.data);
    new_shape.indices = (int*)(indices.data);
    new_shape.num_positions = DBuffer_size(positions);
    new_shape.num_normals = DBuffer_size(normals);
    new_shape.num_texcoords = DBuffer_size(texcoords);
    new_shape.num_indices = DBuffer_size(indices);

    DBuffer_erase(in_face_group);
    HashTable_destroy(&index_table);

    *shape = new_shape;
    return true;
}

int loadOBJ(OBJShape** shapes, const char*  file_name)
{
    FILE* fp = fopen(file_name, "r");
    if(!fp)
    {
        fprintf(stderr, "Cannot open file %s\n", file_name);
        return -1;
    }
    int read_result;
    char line_buffer[1024];
    char cur_mat_name[128];
    cur_mat_name[0] = '\0';
    char mesh_name[128];
    int i;
    for(i = 0; file_name[i] != '.'; i++){}
    strncpy(mesh_name, file_name, i);
    mesh_name[i] = '\0';

    DBuffer obj_shapes = DBuffer_create(OBJShape);

    DBuffer in_positions = DBuffer_create(float);
    DBuffer in_normals = DBuffer_create(float);
    DBuffer in_texcoords = DBuffer_create(float);
    DBuffer in_face_group = DBuffer_create(DBuffer);    
    
    while(OBJGetLine(fp, &read_result, line_buffer, 1024))
    {
        const char* line_ptr = &(line_buffer[0]);
        line_ptr += strspn(line_ptr, " \t");

        // Position
        if(line_ptr[0] == 'v' && line_ptr[1] == ' ')
        {
            line_ptr += 1; // Skip "v"
            float x, y, z;
            OBJParseFloat3(&x, &y, &z, &line_ptr);
            DBuffer_push(in_positions, x);
            DBuffer_push(in_positions, y);
            DBuffer_push(in_positions, z);            
        }
        
        // Normal
        if(line_ptr[0] == 'v' && line_ptr[1] == 'n')
        {
            line_ptr += 2; // Skip "vn"
            float x, y, z;
            OBJParseFloat3(&x, &y, &z, &line_ptr);
            DBuffer_push(in_normals, x);
            DBuffer_push(in_normals, y);
            DBuffer_push(in_normals, z);                        
        }

        // Texcoord
        if(line_ptr[0] == 'v' && line_ptr[1] == 't')
        {
            line_ptr += 2; // Skip "vt"
            float u, v;
            OBJParseFloat2(&u, &v, &line_ptr);
            DBuffer_push(in_texcoords, u);
            DBuffer_push(in_texcoords, v);
        }

        // Face
        if(line_ptr[0] == 'f' && line_ptr[1] == ' ')
        {
            line_ptr += 1; // Skip "f"
            DBuffer face = DBuffer_create_cap(VertexIndex, 4);
            int count = 0;
            while(line_ptr[0] != '\0' && line_ptr < line_buffer + strlen(line_buffer))
            {
                count++;
                VertexIndex vi = OBJParseFaceTriple(&line_ptr);
                vi.v_idx -= 1;
                vi.vt_idx -= 1;
                vi.vn_idx -= 1;                
                DBuffer_push(face, vi);
            }            
            DBuffer_push(in_face_group, face);
        }

        if(line_ptr[0] == 'g')
        {
            OBJShape shape;
            if(exportGroupToShape(&shape, &in_positions, &in_normals, &in_texcoords, &in_face_group))
            {
                int i;
                strcpy(shape.mesh_name, mesh_name);
                strcpy(shape.mat_name, cur_mat_name);
                DBuffer_push(obj_shapes, shape);
            }
        }

        if(line_ptr[0] == 'u' && strncmp(line_ptr, "usemtl", 6) == 0)
        {
            OBJShape shape;
            if(exportGroupToShape(&shape, &in_positions, &in_normals, &in_texcoords, &in_face_group))
            {
                strcpy(shape.mesh_name, mesh_name);
                strcpy(shape.mat_name, cur_mat_name);
                DBuffer_push(obj_shapes, shape);
            }
            line_ptr += 6;    // Skip over "usemtl"
            OBJParseString(cur_mat_name, &line_ptr);
        }

        if(line_ptr[0] == 'o' && line_ptr[1] == ' ')
        {
            OBJShape shape;
            if(exportGroupToShape(&shape, &in_positions, &in_normals, &in_texcoords, &in_face_group))
            {
                DBuffer_push(obj_shapes, shape);
            }
        }
    }
    OBJShape shape;
    if(exportGroupToShape(&shape, &in_positions, &in_normals, &in_texcoords, &in_face_group))
    {
        DBuffer_push(obj_shapes, shape);
    }
    printf("num shapes %d\n", DBuffer_size(obj_shapes));
#if 0
    OBJShape* obj_shape = (OBJShape*)(obj_shapes.data);
    for(int i = 0; i < DBuffer_size(obj_shapes); i++)
    {
        
        printf("Shape %d\n", i);
        printf("position\n");
        for(int j = 0; j < obj_shape[i].num_positions; j++)
        {
            if(j > 0 && j % 3 == 0)
            {
                printf("\n");
            }
            printf("%f ", obj_shape[i].positions[j]);
        }
        printf("\n");        
        
        printf("texcoords\n");
        for(int j = 0; j < obj_shape[i].num_texcoords; j++)
        {
            if(j > 0 && j % 2 == 0)
            {
                printf("\n");
            }
            printf("%f ", obj_shape[i].texcoords[j]);
        }
        printf("\n");

        printf("normal\n");
        for(int j = 0; j < obj_shape[i].num_normals; j++)
        {
            if(j > 0 && j % 3 == 0)
            {
                printf("\n");
            }
            printf("%f ", obj_shape[i].normals[j]);
        }
        printf("\n");
        
        printf("indices\n");
        for(int j = 0; j < obj_shape[i].num_indices; j++)
        {
            if(j > 0 && j % 3 == 0)
            {
                printf("\n");
            }
            printf("%d ", obj_shape[i].indices[j]);
        }
        printf("\n");
    }
#endif    
    /*
    float* in_pos_ptr = (float*)(in_positions.data);
    for(int i = 0; i < DBuffer_size(in_positions); i++)
    {
        if(i % 3 == 0)
        {
            printf("\n");
        }        
        printf("%f ", in_pos_ptr[i]);
    }
    printf("\n");


    for(int i = 0; i < DBuffer_size(in_face_group); i++)
    {
        VertexIndex* vi_ptr = (VertexIndex*)(face_ptr[i].data);
        printf("face size %d\n", DBuffer_size(face_ptr[i]));
        printf("face max %d\n", DBuffer_max_elements(face_ptr[i]));
        for(int j = 0; j < DBuffer_size(face_ptr[i]); j++)
        {
            printf("%d/%d/%d ", vi_ptr[j].v_idx, vi_ptr[j].vt_idx, vi_ptr[j].vn_idx);
        }
        printf("\n");
    }
    printf("\n");
    */
    DBuffer* face_ptr = (DBuffer*)(in_face_group.data);
    DBuffer_destroy(&in_positions);
    DBuffer_destroy(&in_normals);
    DBuffer_destroy(&in_texcoords);

    for(int i = 0; i < DBuffer_size(in_face_group); i++)
    {
        DBuffer_destroy(&(face_ptr[i]));
    }
    DBuffer_destroy(&in_face_group);
    
    *shapes = (OBJShape*)(obj_shapes.data);
    return DBuffer_size(obj_shapes);
}
