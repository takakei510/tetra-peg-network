#ifndef TPN_UNION_FIND_H
#define TPN_UNION_FIND_H
#include <stddef.h>
typedef struct { size_t n, components; size_t *parent,*size; } UnionFind;
int uf_init(UnionFind *uf,size_t n);
void uf_free(UnionFind *uf);
size_t uf_find(UnionFind *uf,size_t x); /* x < n */
int uf_union(UnionFind *uf,size_t a,size_t b); /* 1 merged, 0 already joined, -1 invalid */
void uf_largest_two(UnionFind *uf,size_t *largest,size_t *second);
#endif
