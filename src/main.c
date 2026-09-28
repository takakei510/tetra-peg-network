#include "config.h"
#include "lattice.h"
#include "rng.h"
#include <stdio.h>
#include <inttypes.h>
int main(int argc,char **argv) {
    if(argc>2) {fprintf(stderr,"usage: %s [config]\n",argv[0]);return 2;}
    Config config; char error[256];
    if(config_load(argc==2?argv[1]:"configs/default.cfg",&config,error,sizeof error)) {fprintf(stderr,"config: %s\n",error);return 1;}
    Lattice lattice={0};if(lattice_init(&lattice,config.dim,config.L)) {fprintf(stderr,"cannot allocate lattice (dim=%d, L=%d)\n",config.dim,config.L);return 1;}
    Rng rng; rng_seed(&rng,config.seed);
    printf("Foundation ready: dim=%d L=%d sites=%zu seed=%" PRIu64 "\n",config.dim,config.L,lattice.n_sites,config.seed);
    lattice_free(&lattice);return 0;
}
