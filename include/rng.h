#ifndef TPN_RNG_H
#define TPN_RNG_H
#include <stdint.h>
typedef struct { uint64_t state; } Rng;
void rng_seed(Rng *rng, uint64_t seed);
uint64_t rng_next(Rng *rng);
uint64_t rng_bounded(Rng *rng, uint64_t bound); /* bound > 0 */
double rng_uniform(Rng *rng); /* [0,1) */
uint64_t rng_trial_seed(uint64_t master_seed, uint64_t trial);
#endif
