#include "rng.h"
#include <assert.h>
/* SplitMix64; specified unsigned overflow makes streams reproducible. */
uint64_t rng_next(Rng *r) { uint64_t z=(r->state+=UINT64_C(0x9e3779b97f4a7c15)); z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9); z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb); return z^(z>>31); }
void rng_seed(Rng *r,uint64_t seed) { r->state=seed; }
uint64_t rng_trial_seed(uint64_t master,uint64_t trial) { Rng r; rng_seed(&r,master+trial*UINT64_C(0x9e3779b97f4a7c15)); return rng_next(&r); }
uint64_t rng_bounded(Rng *r,uint64_t bound) { assert(bound); uint64_t threshold=(uint64_t)(-bound)%bound, v; do {v=rng_next(r);} while(v<threshold); return v%bound; }
double rng_uniform(Rng *r) { return (double)(rng_next(r)>>11)*0x1.0p-53; }
