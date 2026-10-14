#include "placement.h"

int placement_generate(Lattice *lat, MoleculeStore *s, Rng *rng,
                       const NetworkConfig *c, size_t max_attempts, PlacementStats *stats) {
    if(!lat || !s || !rng || !c || !stats || !max_attempts ||
       network_config_validate(c,NULL,0) || !lat->owner ||
       lat->dim!=c->base.dim || lat->L!=c->base.L || s->count ||
       s->capacity<c->target_molecules || s->arm_length!=c->arm_length) return -1;
    /* 既存試料への追記ではない。占有点のある格子を誤って再利用しない。 */
    for(size_t site=0;site<lat->n_sites;++site) if(lat->owner[site]>=0) return -1;
    *stats=(PlacementStats){0};
    while(s->count<c->target_molecules && stats->attempts<max_attempts) {
        /* 全格子から中心候補を一様抽選。占有なら棄却して再試行する。 */
        size_t center=(size_t)rng_bounded(rng,lat->n_sites);
        ++stats->attempts;
        GenerationStatus result=molecule_try_generate(lat,s,rng,center,c->a_probability);
        stats->last_status=result;
        if(result==GEN_CENTER_OCCUPIED) ++stats->center_occupied;
        else if(result==GEN_TRAPPED) ++stats->trapped;
        else if(result!=GEN_SUCCESS) return -1; /* メモリ不足等は物理的棄却と分ける。 */
        /* 成功時のみstore.countが増える。失敗で分子IDの欠番を作らない。 */
    }
    return s->count==c->target_molecules ? 0 : 1;
}
