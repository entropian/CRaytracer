#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "dbuffer.h"
#include "hashindex.h"

#define OBJ_NAME_LENGTH 32

static void stringCopy(char* dest, const int max_len, const char* src)
{
#ifdef _MSC_VER
    strcpy_s(dest, max_len, src);
#else
    strcpy(dest, src);
#endif
}

static void stringNCopy(char* dest, const int max_len, const char*src, const int len)
{
#ifdef _MSC_VER
    strncpy_s(dest, max_len, src, len);
#else
    strncpy(dest, src, len);
#endif
}

struct OBJShape_s
{
    float* positions;
    float* normals;
    float* texcoords;
    int* indices;
    //float* face_normals;
    int num_positions, num_normals, num_texcoords, num_indices;
    char mat_name[OBJ_NAME_LENGTH];
    char mesh_name[OBJ_NAME_LENGTH];
};
typedef struct OBJShape_s OBJShape;

typedef struct OBJMaterial_s
{
    float ambient[3];
    float diffuse[3];
    float specular[3];
    float emissive[3];
    
    int illum;
}OBJMaterial;

void OBJShape_destroy(OBJShape* obj_shape);

typedef struct VertexIndex_s
{
    int v_idx, vn_idx, vt_idx;
}VertexIndex;

static bool OBJGetLine(FILE* fp, int* read_result, char buffer[], const int buffer_size)
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

static inline int fixIndex(const int index, const int size)
{
    if(index > 0)
    {
        return index - 1;
    }
    if(index < 0)
    {
        return index + size;
    }
    return -1;
}

static VertexIndex OBJParseFaceTriple(const char** str)
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
        (*str)++;
        vi.vn_idx = atoi(*str);
        *str += strcspn(*str, "/ \t");        
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
static void OBJParseString(char buffer[], const char** str)
{
    *str += strspn(*str, " \t");
    int length = strcspn(*str, " \t\r\0");
    stringNCopy(buffer, 128, *str, length);
    //strncpy_s(buffer, 128, *str, length);
    buffer[length] = '\0';
    *str += length;
}

static inline void getVertexIndexString(char* string, const VertexIndex* vi)
{
#ifdef _MSC_VER
    sprintf_s(string, 256, "p%dn%dt%d", vi->v_idx, vi->vn_idx, vi->vt_idx);
#else
    sprintf(string, "p%dn%dt%d", vi->v_idx, vi->vn_idx, vi->vt_idx);
#endif
}

static inline bool vertexIndexComp(const VertexIndex* a, const VertexIndex* b)
{
    if(a->v_idx != b->v_idx){return false;}
    if(a->vn_idx != b->vn_idx){return false;}
    if(a->vt_idx != b->vt_idx){return false;}
    return true;
}

static int UpdateVertexCache(DBuffer* positions, DBuffer* normals, DBuffer* texcoords, HashIndex* hash_index,
                      DBuffer* vi_cache,
                      const VertexIndex* vi, const DBuffer* in_positions, const DBuffer* in_normals,
                      const DBuffer* in_texcoords)
{

    char vertex_string[256];
    getVertexIndexString(vertex_string, vi);
    
    VertexIndex* vi_cache_ptr = (VertexIndex*)(vi_cache->data);
    int key = HashIndex_gen_key(hash_index, vertex_string);
    int i;
    for(i = HashIndex_First(hash_index, key); i != -1; i = HashIndex_Next(hash_index, i))
    {
        if(vertexIndexComp(vi, &(vi_cache_ptr[i])))
        {
            return i;
        }
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

    int index = DBuffer_size(*positions) / 3 - 1;
    HashIndex_add(hash_index, key, index);
    DBuffer_push(*vi_cache, *vi);
    return index;
}

static bool exportGroupToShape(OBJShape* shape, const DBuffer* in_positions, const DBuffer* in_normals,
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

    DBuffer vi_cache = DBuffer_create(VertexIndex);
    HashIndex hash_index;
    HashIndex_init(&hash_index, DEFAULT_HASH_SIZE, DEFAULT_INDEX_SIZE);
    
    
    DBuffer* face_group_ptr = (DBuffer*)(in_face_group->data);
    float prev = 0.0f, current;
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
            int i0 = UpdateVertexCache(&positions, &normals, &texcoords, &hash_index, &vi_cache,
                                       v0, in_positions, in_normals, in_texcoords);
            int i1 = UpdateVertexCache(&positions, &normals, &texcoords, &hash_index, &vi_cache,
                                       v1, in_positions, in_normals, in_texcoords);
            int i2 = UpdateVertexCache(&positions, &normals, &texcoords, &hash_index, &vi_cache,
                                       v2, in_positions, in_normals, in_texcoords);            

            DBuffer_push(indices, i0);
            DBuffer_push(indices, i1);
            DBuffer_push(indices, i2);                    
        }
        current = (float)i / (DBuffer_size(*in_face_group));
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
    DBuffer_destroy(&vi_cache);
    HashIndex_free(&hash_index);

    *shape = new_shape;
    return true;
}

int loadOBJ(OBJShape** shapes, const char*  file_name);

#ifdef OBJ_LOADER_IMPLEMENTATION
int loadOBJ(OBJShape** shapes, const char*  file_name)
{
    FILE* fp;
#ifdef _MSC_VER
    fopen_s(&fp, file_name, "r");
#else
    fp = fopen(file_name, "r");
#endif
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
    stringNCopy(mesh_name, OBJ_NAME_LENGTH, file_name, i);
    mesh_name[i] = '\0';

    DBuffer obj_shapes = DBuffer_create(OBJShape);

    DBuffer in_positions = DBuffer_create(float);
    DBuffer in_normals = DBuffer_create(float);
    DBuffer in_texcoords = DBuffer_create(float);

    /*
                            .       .       .       .       .       .       .      
                            .       .       .       .       .       .       .      
                            .       .       .       .       .       .       .      
                            .       .       .       .       .       .       .
                            VI      VI      VI      VI      VI      VI      VI
                            VI      VI      VI      VI      VI      VI      VI
                            VI      VI      VI      VI      VI      VI      VI                            
      in_face_group.data -> DBuffer DBuffer DBuffer DBuffer DBuffer DBuffer DBuffer......
                    size
                    max
                    element_size
     */
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
                vi.v_idx = fixIndex(vi.v_idx, DBuffer_size(in_positions) / 3);
                vi.vt_idx = fixIndex(vi.vt_idx, DBuffer_size(in_texcoords) / 2);
                vi.vn_idx = fixIndex(vi.vn_idx, DBuffer_size(in_normals) / 3);
                /*
                if(vi.v_idx < 0)
                {
                    printf("here\n");
                }
                */
                // TODO: fix this
                if(vi.v_idx != -1)
                {
                    DBuffer_push(face, vi);
                }
            }            
            DBuffer_push(in_face_group, face);
        }

        if(line_ptr[0] == 'g')
        {
            OBJShape shape;
            if(exportGroupToShape(&shape, &in_positions, &in_normals, &in_texcoords, &in_face_group))
            {
                stringCopy(shape.mesh_name, OBJ_NAME_LENGTH, mesh_name);
                stringCopy(shape.mat_name, OBJ_NAME_LENGTH, cur_mat_name);
                DBuffer_push(obj_shapes, shape);
            }
        }

        if(line_ptr[0] == 'u' && strncmp(line_ptr, "usemtl", 6) == 0)
        {
            OBJShape shape;
            if(exportGroupToShape(&shape, &in_positions, &in_normals, &in_texcoords, &in_face_group))
            {
                stringCopy(shape.mesh_name, OBJ_NAME_LENGTH, mesh_name);
                stringCopy(shape.mat_name, OBJ_NAME_LENGTH, cur_mat_name);
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
        stringCopy(shape.mesh_name, OBJ_NAME_LENGTH, mesh_name);
        stringCopy(shape.mat_name, OBJ_NAME_LENGTH, cur_mat_name);        
        DBuffer_push(obj_shapes, shape);
    }
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
        For(int j = 0; j < DBuffer_size(face_ptr[i]); j++)
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

void OBJShape_destroy(OBJShape* obj_shape)
{
    if(obj_shape->num_positions > 0){free(obj_shape->positions);}
    if(obj_shape->num_normals > 0){free(obj_shape->normals);}
    if(obj_shape->num_texcoords > 0){free(obj_shape->texcoords);}
    if(obj_shape->num_indices > 0){free(obj_shape->indices);}

    obj_shape->num_positions = 0;
    obj_shape->num_normals = 0;
    obj_shape->num_texcoords = 0;
    obj_shape->num_indices = 0;
}
#endif
