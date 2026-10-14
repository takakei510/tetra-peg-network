/* 次の①：設定だけを読み込み・検査する。分子生成や結合はまだ呼ばない。 */
#include "network_config.h"
#include <inttypes.h>
#include <stdio.h>
int main(int argc, char **argv) {
    if(argc!=2) {fprintf(stderr,"usage: %s network.cfg\n",argv[0]);return 2;}
    NetworkConfig c; char error[256];
    if(network_config_load(argv[1],&c,error,sizeof error)) {fprintf(stderr,"%s\n",error);return 1;}
    size_t sites=(size_t)c.base.L*c.base.L*c.base.L;
    size_t required=c.target_molecules*(1+4*c.arm_length); /* 検査済みなので安全。 */
    printf("config_validated=1 dim=%d L=%d seed=%" PRIu64 " target_molecules=%zu "
        "arm_length=%zu a_probability=%.17g required_sites=%zu resolved_phi=%.17g "
        "generation_run=0",c.base.dim,c.base.L,c.base.seed,c.target_molecules,
        c.arm_length,c.a_probability,required,(double)required/sites);
    if(c.target_mode==TARGET_OCCUPIED_FRACTION)
        printf(" input_mode=occupied_fraction requested_phi=%.17g rounding=floor phi_site_budget=%zu\n",
            c.requested_phi,c.phi_site_budget);
    else printf(" input_mode=molecule_count requested_phi=not_set rounding=none\n");
    return 0;
}
