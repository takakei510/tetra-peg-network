#include "union_find.h"
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
int uf_init(UnionFind *uf,size_t n) {
    if(!uf||!n||n>SIZE_MAX/sizeof(size_t)) return -1;
    size_t *p=malloc(n*sizeof *p),*s=malloc(n*sizeof *s);
    if(!p||!s) {free(p);free(s);return -1;}
    for(size_t i=0;i<n;i++) {p[i]=i;s[i]=1;}
    *uf=(UnionFind){n,n,p,s};return 0;
}
void uf_free(UnionFind *uf) {if(uf) {free(uf->parent);free(uf->size);*uf=(UnionFind){0};}}
size_t uf_find(UnionFind *uf,size_t x) {assert(uf&&x<uf->n); while(uf->parent[x]!=x) {uf->parent[x]=uf->parent[uf->parent[x]];x=uf->parent[x];} return x;}
int uf_union(UnionFind *uf,size_t a,size_t b) {if(!uf||a>=uf->n||b>=uf->n)return -1; a=uf_find(uf,a);b=uf_find(uf,b);if(a==b)return 0;if(uf->size[a]<uf->size[b]) {size_t t=a;a=b;b=t;}uf->parent[b]=a;uf->size[a]+=uf->size[b];uf->size[b]=0;uf->components--;return 1;}
void uf_largest_two(UnionFind *uf,size_t *largest,size_t *second) { *largest=*second=0; for(size_t i=0;i<uf->n;i++) if(uf->parent[i]==i) {size_t s=uf->size[i];if(s>*largest){*second=*largest;*largest=s;}else if(s>*second)*second=s;} }
