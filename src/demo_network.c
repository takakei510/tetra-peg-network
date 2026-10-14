/* 多分子実行入口。生成・探索・結合・UFは既存モジュールを呼び出す。 */
#include "placement.h"
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include "cluster.h"
#include <inttypes.h>

int main(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr, "usage: %s config candidates.csv bonds.csv cluster-members.csv cluster-summary.csv network-summary.csv placement-summary.csv max_attempts\n", argv[0]);
        return 2;
    }
    NetworkConfig config;
    char error[256];
    if (network_config_load(argv[1], &config, error, sizeof error)) {
        fprintf(stderr, "%s\n", error); return 1;
    }
    /* CLIの試行上限も、符号・小数・桁溢れを拒否して読む。 */
    for(const char *p=argv[8];*p;++p) if(!isdigit((unsigned char)*p)) return 2;
    char *end; errno=0; uintmax_t value=strtoumax(argv[8],&end,10);
    if(errno || end==argv[8] || *end || !value || value>SIZE_MAX) return 2;
    size_t max_attempts=(size_t)value;
    Lattice lat = {0};
    MoleculeStore store = {0};
    Rng rng;
    if (lattice_init(&lat, config.base.dim, config.base.L) ||
        molecule_store_init(&store, config.target_molecules, config.arm_length)) {
        fprintf(stderr, "cannot allocate lattice or molecule store\n");
        molecule_store_free(&store);
        lattice_free(&lat);
        return 1;
    }
    rng_seed(&rng, config.base.seed);
    PlacementStats placement={0};
    int placement_result=placement_generate(&lat,&store,&rng,&config,max_attempts,&placement);
    int fatal=placement_result<0;
    size_t occupied=0, type_a=0;
    for(size_t site=0;site<lat.n_sites;++site) occupied+=lat.owner[site]>=0;
    for(size_t m=0;m<store.count;++m) type_a+=store.molecules[m].type==TYPE_A;
    size_t required=config.target_molecules*(1+ARM_COUNT*config.arm_length);
    const char *input_mode=config.target_mode==TARGET_OCCUPIED_FRACTION ? "occupied_fraction" : "molecule_count";
    const char *stop_reason=fatal ? "generation_error" : placement_result ? "attempt_limit" : "target_reached";
    fprintf(stderr,"dim=3 L=%d seed=%" PRIu64 " arm_length=%zu a_probability=%.17g "
        "type_mode=random_AB boundary=open growth=free_neighbor order=sequential "
        "input_mode=%s target=%zu placed=%zu max_placement_attempts=%zu attempts=%zu "
        "center_occupied=%zu trapped=%zu occupied=%zu actual_phi=%.17g resolved_phi=%.17g "
        "type_A=%zu type_B=%zu completed=%d stop_reason=%s generation_run=1\n",
        config.base.L,config.base.seed,config.arm_length,config.a_probability,input_mode,
        config.target_molecules,store.count,max_attempts,placement.attempts,
        placement.center_occupied,placement.trapped,occupied,(double)occupied/lat.n_sites,
        (double)required/lat.n_sites,type_a,store.count-type_a,placement_result==0,stop_reason);
    /* 設定の目標と実績を独立したCSVに保存。未達も消さない。 */
    FILE *summary=fopen(argv[7],"w"); int placement_csv_error=summary==NULL;
    if(summary) {
        if(fputs("input_mode,requested_phi,target_molecules,resolved_phi,placed_molecules,occupied_sites,actual_phi,type_A,type_B,max_placement_attempts,attempts,center_occupied,trapped,completed,stop_reason\n",summary)==EOF) placement_csv_error=1;
        if(fprintf(summary,"%s,",input_mode)<0) placement_csv_error=1;
        if(config.target_mode==TARGET_OCCUPIED_FRACTION && fprintf(summary,"%.17g",config.requested_phi)<0) placement_csv_error=1;
        if(fprintf(summary,",%zu,%.17g,%zu,%zu,%.17g,%zu,%zu,%zu,%zu,%zu,%zu,%d,%s\n",
            config.target_molecules,(double)required/lat.n_sites,store.count,occupied,
            (double)occupied/lat.n_sites,type_a,store.count-type_a,max_attempts,
            placement.attempts,placement.center_occupied,placement.trapped,placement_result==0,stop_reason)<0) placement_csv_error=1;
        if(fclose(summary)) placement_csv_error=1;
    }
    if(config.target_mode==TARGET_OCCUPIED_FRACTION)
        fprintf(stderr,"requested_phi=%.17g phi_site_budget=%zu rounding=floor\n",config.requested_phi,config.phi_site_budget);
    else fprintf(stderr,"requested_phi=not_set rounding=none\n");
    /* 未完了の場合も完成分子のCSVを出すが、終了コード3（未達）とmetadataで区別する。 */
    int csv_error = molecule_write_csv(stdout, &lat, &store);
    /* 生成完了後に索引を構築。候補を保存した後、独立した結合処理へ渡す。 */
    EndpointIndex index={0}; CandidateStore candidates={0}; SearchStats stats={0};
    BondStore bonds={0}; int bond_error=0;
    int search_error=0;
    if(!fatal && placement_result==0) {
        search_error=endpoint_index_build(&lat,&store,&index);
        if(!search_error) search_error=endpoint_find_candidates(&lat,&store,&index,&candidates,&stats);
        if(!search_error) {
            FILE *fp=fopen(argv[2],"w");
            if(!fp) search_error=1;
            else {
                search_error=candidate_write_csv(fp,&lat,&store,&candidates);
                if(fclose(fp)) search_error=1;
            }
        }
        fprintf(stderr,"distance_metric=manhattan bond_radius=2 endpoints=%zu candidates=%zu offsets_tested=%zu site_lookups=%zu search_completed=%d\n",
            index.endpoint_count,candidates.count,stats.offsets_tested,stats.site_lookups,!search_error);
    }
    if(!fatal && placement_result==0 && !search_error) {
        /* 生成に続く同じ乱数系列で末端順と相手を抽選。多重結合は許可。 */
        bond_error=bond_select(&candidates,index.endpoint_count,&rng,1,&bonds);
        if(!bond_error) {
            FILE *fp=fopen(argv[3],"w");
            if(!fp) bond_error=1;
            else {
                bond_error=bond_write_csv(fp,&lat,&store,&bonds);
                if(fclose(fp)) bond_error=1;
            }
        }
        fprintf(stderr,"pairing=endpoint_shuffle_uniform_partner allow_multiple_bonds=1 "
            "bond_rng=continued_generation bonds=%zu bonded_endpoints=%zu "
            "bond_fraction=%.9g bonding_completed=%d\n",
            bonds.count,2*bonds.count,2.0*bonds.count/index.endpoint_count,!bond_error);
    }
    ClusterStore clusters={0}; int cluster_error=0;
    if(!fatal && placement_result==0 && !search_error && !bond_error) {
        /* 位置と結合の抽選には触れず、成立済みの結合からクラスタを集計。 */
        cluster_error=cluster_build(&bonds,store.count,&clusters);
        if(!cluster_error) {
            int (*writers[])(FILE *,const ClusterStore *)={cluster_members_write_csv,
                cluster_summary_write_csv,network_summary_write_csv};
            for(int k=0;k<3;++k) {
                FILE *fp=fopen(argv[4+k],"w");
                if(!fp) {cluster_error=1;break;}
                if(writers[k](fp,&clusters)) cluster_error=1;
                if(fclose(fp)) cluster_error=1;
                if(cluster_error) break;
            }
        }
        fprintf(stderr,"cluster_method=union_find cluster_id_rule=min_member_order "
            "clusters=%zu largest=%zu second_largest=%zu isolated_molecules=%zu "
            "largest_fraction=%.9g clustering_completed=%d\n",
            clusters.cluster_count,clusters.largest,clusters.second_largest,
            clusters.isolated_molecules,(double)clusters.largest/store.count,!cluster_error);
    }
    cluster_store_free(&clusters);
    bond_store_free(&bonds);
    candidate_store_free(&candidates); endpoint_index_free(&index);
    int result = fatal || placement_csv_error || csv_error || search_error || bond_error || cluster_error;
    /* 未達の試料は結合・UFを実行しない。目標濃度の完成試料と混ぜない。 */
    if(placement_result!=0) fprintf(stderr,"analysis_skipped=1\n");
    molecule_store_free(&store);
    lattice_free(&lat);
    return result ? 1 : placement_result ? 3 : 0;
}
