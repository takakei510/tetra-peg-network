/*
 * ③-B：2分子配置とランダム末端結合の動作確認用デモ。分子生成そのものはmolecule.cを再利用する。
 * 1個目の占有を残して2個目を生成し、全分子共通ownerで重複を防ぐ。
 * 可視化でA/Bを見分けられるよう、型はA、Bの順に固定する。
 * これは確認用条件であり、確率rによるA/B混合の統計試料ではない。
 */
#include "config.h"
#include "bond.h"
#include <inttypes.h>

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s config candidates.csv bonds.csv\n", argv[0]);
        return 2;
    }
    Config config;
    char error[256];
    if (config_load(argv[1], &config, error, sizeof error)) {
        fprintf(stderr, "%s\n", error);
        return 1;
    }
    if (config.dim != 3) {
        fprintf(stderr, "two-molecule demo requires dim=3\n");
        return 1;
    }
    Lattice lat = {0};
    MoleculeStore store = {0};
    Rng rng;
    if (lattice_init(&lat, config.dim, config.L) ||
        molecule_store_init(&store, 2, 5)) {
        fprintf(stderr, "cannot allocate lattice or molecule store\n");
        molecule_store_free(&store);
        lattice_free(&lat);
        return 1;
    }
    rng_seed(&rng, config.seed);
    size_t attempts = 0, occupied_rejections = 0, trapped_rejections = 0;
    int fatal = 0;
    /*
     * 完成分子だけがstore.countに含まれる。失敗してもcountは変わらず、
     * 同じ分子IDで次の中心を試す。上限は2分子合計で10000回。
     */
    while (store.count < 2 && attempts < 10000) {
        size_t center = (size_t)rng_bounded(&rng, lat.n_sites);
        double a_probability = store.count == 0 ? 1.0 : 0.0;
        ++attempts;
        GenerationStatus status = molecule_try_generate(
            &lat, &store, &rng, center, a_probability);
        if (status == GEN_CENTER_OCCUPIED) ++occupied_rejections;
        else if (status == GEN_TRAPPED) ++trapped_rejections;
        else if (status != GEN_SUCCESS) {
            fprintf(stderr, "generation error: status=%d\n", status);
            fatal = 1;
            break;
        }
    }
    /* 実績占有数を格子から数える。成功時は2*(1+4*5)=42点になる。 */
    size_t occupied = 0;
    for (size_t site = 0; site < lat.n_sites; ++site)
        occupied += lat.owner[site] >= 0;
    fprintf(stderr,
        "L=%d seed=%" PRIu64 " arm_length=5 type_mode=fixed_AB "
        "boundary=open growth=free_neighbor order=sequential "
        "target=2 placed=%zu attempts=%zu center_occupied=%zu trapped=%zu "
        "occupied=%zu phi=%.9g completed=%d\n",
        config.L, config.seed, store.count, attempts, occupied_rejections,
        trapped_rejections, occupied, (double)occupied / (double)lat.n_sites,
        !fatal && store.count == 2);
    /* 未完了の場合も完成分子のCSVを出すが、終了コード1とmetadataで区別する。 */
    int csv_error = molecule_write_csv(stdout, &lat, &store);
    /* 生成完了後に索引を構築。候補を保存した後、独立した結合処理へ渡す。 */
    EndpointIndex index={0}; CandidateStore candidates={0}; SearchStats stats={0};
    BondStore bonds={0}; int bond_error=0;
    int search_error=0;
    if(!fatal && store.count==2) {
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
    if(!fatal && store.count==2 && !search_error) {
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
            "bond_fraction=%.9g bonding_completed=%d uf_implemented=0\n",
            bonds.count,2*bonds.count,2.0*bonds.count/index.endpoint_count,!bond_error);
    }
    bond_store_free(&bonds);
    candidate_store_free(&candidates); endpoint_index_free(&index);
    int result = fatal || store.count != 2 || csv_error || search_error || bond_error;
    molecule_store_free(&store);
    lattice_free(&lat);
    return result ? 1 : 0;
}
