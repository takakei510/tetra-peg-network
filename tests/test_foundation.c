#include "config.h"
#include "rng.h"
#include "lattice.h"
#include "union_find.h"
#include <assert.h>
#include <stdio.h>
int main(void) {
    Config cfg;char error[128];assert(config_load("configs/default.cfg",&cfg,error,sizeof error)==0);assert(cfg.dim==3&&cfg.L==32&&cfg.seed==12345);
    Rng a,b;rng_seed(&a,0);rng_seed(&b,0);for(int i=0;i<1000;i++)assert(rng_next(&a)==rng_next(&b));assert(rng_trial_seed(42,0)!=rng_trial_seed(42,1));for(int i=0;i<1000;i++)assert(rng_bounded(&a,7)<7);
    Lattice lat={0};assert(lattice_init(&lat,3,4)==0);assert(lat.n_sites==64);int c[3]={3,2,1},back[3];size_t index;assert(lattice_index(&lat,c,&index)==0&&index==57);assert(lattice_coord(&lat,index,back)==0&&back[0]==3&&back[1]==2&&back[2]==1);size_t neighbors[6];assert(lattice_neighbors(&lat,0,neighbors)==3);assert(lattice_neighbors(&lat,21,neighbors)==6);c[0]=4;assert(lattice_index(&lat,c,&index)!=0);assert(lat.owner[0]==-1);lattice_free(&lat);
    UnionFind uf={0};assert(uf_init(&uf,5)==0);assert(uf_union(&uf,0,1)==1);assert(uf_union(&uf,2,3)==1);assert(uf_union(&uf,3,4)==1);assert(uf_union(&uf,2,4)==0);size_t first,second;uf_largest_two(&uf,&first,&second);assert(uf.components==2&&first==3&&second==2);uf_free(&uf);
    puts("foundation tests passed");return 0;
}
