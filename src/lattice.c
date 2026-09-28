#include "lattice.h"
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
int lattice_init(Lattice *lat,int dim,int L) {
    if(!lat|| (dim!=2&&dim!=3)||L<2) return -1;
    size_t n=1; for(int i=0;i<dim;i++) {if(n>SIZE_MAX/(size_t)L) return -1; n*=L;}
    if(n>SIZE_MAX/sizeof(int)||n>(size_t)INT_MAX) return -1;
    int *owner=malloc(n*sizeof *owner); if(!owner) return -1;
    for(size_t i=0;i<n;i++) owner[i]=-1;
    *lat=(Lattice){dim,L,n,owner}; return 0;
}
void lattice_free(Lattice *lat) { if(lat) {free(lat->owner); *lat=(Lattice){0};} }
int lattice_index(const Lattice *lat,const int coord[3],size_t *out) {
    if(!lat||!coord||!out||!lat->owner) return -1;
    size_t index=0; for(int i=0;i<lat->dim;i++) {if(coord[i]<0||coord[i]>=lat->L) return -1; index=index*(size_t)lat->L+(size_t)coord[i];}
    *out=index; return 0;
}
int lattice_coord(const Lattice *lat,size_t index,int coord[3]) {
    if(!lat||!coord||!lat->owner||index>=lat->n_sites) return -1;
    coord[0]=coord[1]=coord[2]=0;
    for(int i=lat->dim-1;i>=0;i--) {coord[i]=(int)(index%(size_t)lat->L); index/=(size_t)lat->L;} return 0;
}
int lattice_neighbors(const Lattice *lat,size_t index,size_t out[6]) {
    int c[3]; if(!out||lattice_coord(lat,index,c)) return -1;
    int count=0; for(int d=0;d<lat->dim;d++) { int saved=c[d]; if(saved>0) {c[d]--; lattice_index(lat,c,&out[count++]);c[d]=saved;} if(saved<lat->L-1) {c[d]++;lattice_index(lat,c,&out[count++]);c[d]=saved;} } return count;
}
