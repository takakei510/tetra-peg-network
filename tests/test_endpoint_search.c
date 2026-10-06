#include "endpoint_search.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/* 総当たりを独立した検証用の正解とする。実用探索とは異なる経路で距離を計算。 */
static void compare_brute_force(Lattice *lat, MoleculeStore *s) {
    EndpointIndex index={0}; CandidateStore candidates={0}; SearchStats stats={0};
    int *owners=malloc(lat->n_sites*sizeof(int)); assert(owners);
    memcpy(owners,lat->owner,lat->n_sites*sizeof(int));
    assert(!endpoint_index_build(lat,s,&index));
    assert(!endpoint_find_candidates(lat,s,&index,&candidates,&stats));
    size_t n=4*s->count, expected=0;
    unsigned char *seen=calloc(n*n+1,1); assert(seen);
    for(size_t k=0;k<candidates.count;++k) {
        size_t a=candidates.pairs[k].endpoint_a,b=candidates.pairs[k].endpoint_b;
        assert(a<b && b<n && !seen[a*n+b]);
        seen[a*n+b]=1;
    }
    for(size_t a=0;a<n;++a) for(size_t b=a+1;b<n;++b) {
        int ca[3],cb[3];
        size_t sa=s->path_sites[molecule_path_index(s,a/4,a%4,s->arm_length-1)];
        size_t sb=s->path_sites[molecule_path_index(s,b/4,b%4,s->arm_length-1)];
        assert(!lattice_coord(lat,sa,ca) && !lattice_coord(lat,sb,cb));
        int distance=0; for(int d=0;d<3;++d) distance+=abs(ca[d]-cb[d]);
        int possible=s->molecules[a/4].type!=s->molecules[b/4].type && distance>=1 && distance<=2;
        assert(seen[a*n+b]==possible); expected+=possible;
    }
    assert(candidates.count==expected && stats.offsets_tested==24*n);
    assert(stats.site_lookups<=stats.offsets_tested);
    assert(!memcmp(owners,lat->owner,lat->n_sites*sizeof(int)));
    FILE *fp=tmpfile(); assert(fp && !candidate_write_csv(fp,lat,s,&candidates)); fclose(fp);
    free(seen); free(owners); candidate_store_free(&candidates); endpoint_index_free(&index);
}

int main(void) {
    Lattice lat={0}; MoleculeStore s={0};
    assert(!lattice_init(&lat,3,8) && !molecule_store_init(&s,2,1));
    /* 距離判定専用の人工末端配置。分子生成の試料とは区別する。
     * (0,0,0)から距離1、距離2、距離3と、ユークリッドなら2以内の
     * (1,1,1)を置き、開境界の向こう側(7,0,0)を候補にしないことも確認。 */
    const int coordinates[8][3]={{0,0,0},{7,7,7},{7,6,7},{6,7,7},
                                {1,0,0},{0,2,0},{1,1,1},{7,0,0}};
    s.count=2; s.molecules[0].type=TYPE_A; s.molecules[1].type=TYPE_B;
    for(size_t id=0;id<8;++id) {
        size_t site; assert(!lattice_index(&lat,coordinates[id],&site));
        s.path_sites[id]=site; lat.owner[site]=(int)(id/4);
    }
    compare_brute_force(&lat,&s);
    EndpointIndex i={0}; CandidateStore c={0}; SearchStats stats={0};
    assert(!endpoint_index_build(&lat,&s,&i));
    assert(!endpoint_find_candidates(&lat,&s,&i,&c,&stats));
    assert(c.count==2); /* 0-4は距離1、0-5は距離2だけが候補。 */
    candidate_store_free(&c); endpoint_index_free(&i);
    s.molecules[1].type=TYPE_A; compare_brute_force(&lat,&s); /* 同型だけなら0候補 */
    size_t saved=s.path_sites[4]; s.path_sites[4]=s.path_sites[0];
    assert(endpoint_index_build(&lat,&s,&i) && !i.endpoint_at_site); /* 重複を拒否 */
    s.path_sites[4]=lat.n_sites;
    assert(endpoint_index_build(&lat,&s,&i) && !i.endpoint_at_site); /* 範囲外を拒否 */
    s.path_sites[4]=saved; s.count=0; compare_brute_force(&lat,&s);
    molecule_store_free(&s); lattice_free(&lat);
    /* 本物の分子軌跡でも照合。型の順序・境界・込み合いが異なる50 seed。 */
    for(uint64_t seed=0;seed<50;++seed) {
        assert(!lattice_init(&lat,3,8) && !molecule_store_init(&s,6,5));
        Rng rng; rng_seed(&rng,seed);
        size_t attempts=0;
        while(s.count<6 && attempts++<10000) {
            GenerationStatus status=molecule_try_generate(&lat,&s,&rng,
                (size_t)rng_bounded(&rng,lat.n_sites),s.count%2 ? 1 : 0);
            assert(status==GEN_SUCCESS || status==GEN_TRAPPED || status==GEN_CENTER_OCCUPIED);
        }
        assert(s.count==6);
        size_t before[120]; memcpy(before,s.path_sites,sizeof before);
        compare_brute_force(&lat,&s);
        assert(!memcmp(before,s.path_sites,sizeof before));
        molecule_store_free(&s); lattice_free(&lat);
    }
    puts("endpoint search tests passed (synthetic cases + 50 seeds vs brute force)");
    return 0;
}
