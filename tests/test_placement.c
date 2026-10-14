#include "placement.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

static void verify(const Lattice *lat, const MoleculeStore *s) {
    size_t occupied=0; unsigned char *seen=calloc(lat->n_sites,1);assert(seen);
    for(size_t m=0;m<s->count;++m) {
        size_t center=s->molecules[m].center_site;
        assert(center<lat->n_sites && !seen[center]++ && lat->owner[center]==(int)m);
        for(size_t a=0;a<ARM_COUNT;++a) for(size_t k=0;k<s->arm_length;++k) {
            size_t site=s->path_sites[molecule_path_index(s,m,a,k)];
            assert(site<lat->n_sites && !seen[site]++ && lat->owner[site]==(int)m);
        }
    }
    for(size_t i=0;i<lat->n_sites;++i) {assert((lat->owner[i]>=0)==(seen[i]!=0));occupied+=lat->owner[i]>=0;}
    assert(occupied==s->count*(1+ARM_COUNT*s->arm_length));free(seen);
}
int main(void) {
    NetworkConfig c;char error[256];
    assert(!network_config_load("configs/network/preview_phi010_L16_seed3.cfg",&c,error,sizeof error));
    for(size_t seed=0;seed<24;++seed) {
        Lattice a={0},b={0};MoleculeStore x={0},y={0};Rng r,t;PlacementStats sx,sy;
        assert(!lattice_init(&a,3,c.base.L) && !lattice_init(&b,3,c.base.L));
        assert(!molecule_store_init(&x,c.target_molecules,c.arm_length) && !molecule_store_init(&y,c.target_molecules,c.arm_length));
        rng_seed(&r,seed);rng_seed(&t,seed);
        assert(!placement_generate(&a,&x,&r,&c,10000,&sx));
        assert(!placement_generate(&b,&y,&t,&c,10000,&sy));
        assert(x.count==19 && sx.attempts==x.count+sx.center_occupied+sx.trapped);
        assert(sx.attempts==sy.attempts && sx.center_occupied==sy.center_occupied && sx.trapped==sy.trapped);
        assert(!memcmp(a.owner,b.owner,a.n_sites*sizeof *a.owner));
        assert(!memcmp(x.path_sites,y.path_sites,x.count*ARM_COUNT*x.arm_length*sizeof *x.path_sites));
        verify(&a,&x);
        molecule_store_free(&x);molecule_store_free(&y);lattice_free(&a);lattice_free(&b);
    }
    /* 上限1回なら19分子未達。完成軌跡以外の仮占有が残らない。 */
    Lattice lat={0};MoleculeStore s={0};Rng rng;PlacementStats stats;
    assert(!lattice_init(&lat,3,c.base.L) && !molecule_store_init(&s,c.target_molecules,c.arm_length));
    rng_seed(&rng,3);assert(placement_generate(&lat,&s,&rng,&c,1,&stats)==1);
    assert(stats.attempts==1 && s.count==1);verify(&lat,&s);
    molecule_store_free(&s);lattice_free(&lat);
    /* L2の角は隣接3点なので4腕は完成不能。trappedを数えて完全復元。 */
    c=(NetworkConfig){.base={.dim=3,.L=2,.seed=3},.target_molecules=1,.arm_length=1,.a_probability=1};
    assert(!lattice_init(&lat,3,2) && !molecule_store_init(&s,1,1));rng_seed(&rng,3);
    assert(placement_generate(&lat,&s,&rng,&c,7,&stats)==1);
    assert(!s.count && stats.trapped==7 && !stats.center_occupied);verify(&lat,&s);
    assert(placement_generate(&lat,&s,&rng,&c,0,&stats)==-1);
    molecule_store_free(&s);lattice_free(&lat);
    puts("placement tests passed (24 seeds, reproducibility, occupancy, partial and zero completion)");
}
