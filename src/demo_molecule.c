/* ミニマム検証用。実験仕様が確定するまでは本体mainと分離する。 */
#include "config.h"
#include "molecule.h"
#include <inttypes.h>
int main(int argc,char **argv) {
    if (argc!=2) {fprintf(stderr,"usage: %s config\n",argv[0]);return 2;}
    Config c; char error[256];
    if (config_load(argv[1],&c,error,sizeof error)) {fprintf(stderr,"%s\n",error);return 1;}
    Lattice lat={0}; MoleculeStore store={0}; Rng rng;
    /* 暫定例：4腕・各5歩・A確率0.5。固定条件を標準エラーへ記録。 */
    if (c.dim!=3 || lattice_init(&lat,c.dim,c.L) || molecule_store_init(&store,1,5)) {
        lattice_free(&lat); return 1;
    }
    rng_seed(&rng,c.seed); GenerationStatus status=GEN_TRAPPED;
    size_t attempts=0;
    while (attempts<10000 && status!=GEN_SUCCESS) {
        ++attempts;
        status=molecule_try_generate(&lat,&store,&rng,
            (size_t)rng_bounded(&rng,lat.n_sites),0.5);
        if (status==GEN_INVALID_CONFIG || status==GEN_ALLOCATION_ERROR) break;
    }
    fprintf(stderr,"L=%d seed=%" PRIu64 " arm_length=5 r=0.5 boundary=open growth=free_neighbor order=sequential attempts=%zu status=%d\n",
        c.L,c.seed,attempts,status);
    int result=status!=GEN_SUCCESS || molecule_write_csv(stdout,&lat,&store);
    molecule_store_free(&store); lattice_free(&lat); return result?1:0;
}
