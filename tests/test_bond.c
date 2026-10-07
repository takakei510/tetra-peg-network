#include "bond.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>

/* 入力候補との照合と逆引きの検証。抽選手順を再実装するテストにはしない。 */
static void check(const CandidateStore *c, const BondStore *s, int multiple) {
    size_t used=0;
    for(size_t a=0;a<s->endpoint_count;++a) {
        size_t k=s->bond_at_endpoint[a];
        if(k==SIZE_MAX) continue;
        ++used; assert(k<s->count);
        assert(s->pairs[k].endpoint_a==a || s->pairs[k].endpoint_b==a);
    }
    assert(used==2*s->count);
    for(size_t k=0;k<s->count;++k) {
        Bond p=s->pairs[k]; int found=0;
        assert(p.endpoint_a<p.endpoint_b && p.endpoint_b<s->endpoint_count);
        assert(s->bond_at_endpoint[p.endpoint_a]==k && s->bond_at_endpoint[p.endpoint_b]==k);
        for(size_t j=0;j<c->count;++j)
            found|=c->pairs[j].endpoint_a==p.endpoint_a && c->pairs[j].endpoint_b==p.endpoint_b;
        assert(found);
        if(!multiple) for(size_t j=0;j<k;++j)
            assert(s->pairs[j].endpoint_a/4!=p.endpoint_a/4 || s->pairs[j].endpoint_b/4!=p.endpoint_b/4);
    }
    /* 許可時は未結合同士の候補が残らない（極大。ただし最大とは限らない）。 */
    if(multiple) for(size_t j=0;j<c->count;++j)
        assert(s->bond_at_endpoint[c->pairs[j].endpoint_a]!=SIZE_MAX ||
               s->bond_at_endpoint[c->pairs[j].endpoint_b]!=SIZE_MAX);
}

int main(void) {
    CandidatePair pairs[]={{0,5},{0,4},{1,4}};
    CandidateStore c={pairs,3,3};
    int saw_one=0,saw_two=0;
    for(uint64_t seed=0;seed<100;++seed) {
        Rng r1,r2; rng_seed(&r1,seed); rng_seed(&r2,seed);
        BondStore a={0},b={0};
        assert(!bond_select(&c,8,&r1,1,&a) && !bond_select(&c,8,&r2,1,&b));
        check(&c,&a,1); assert(a.count==1 || a.count==2);
        saw_one|=a.count==1; saw_two|=a.count==2;
        assert(a.count==b.count && r1.state==r2.state);
        assert(!memcmp(a.pairs,b.pairs,a.count*sizeof *a.pairs));
        assert(!memcmp(a.bond_at_endpoint,b.bond_at_endpoint,8*sizeof(size_t)));
        assert(!memcmp(pairs,(CandidatePair[]){{0,5},{0,4},{1,4}},sizeof pairs));
        bond_store_free(&a); bond_store_free(&b);
        rng_seed(&r1,seed);
        assert(!bond_select(&c,8,&r1,0,&a)); check(&c,&a,0); assert(a.count==1);
        bond_store_free(&a);
    }
    assert(saw_one && saw_two); /* 最大2本の候補でも抽選順により1本になり得る。 */
    Rng rng; rng_seed(&rng,3); BondStore b={0};
    CandidateStore empty={0};
    assert(!bond_select(&empty,8,&rng,1,&b) && !b.count); check(&empty,&b,1);
    bond_store_free(&b);
    assert(!bond_select(&empty,0,&rng,1,&b) && !b.count); bond_store_free(&b);
    CandidatePair bad[]={{0,4},{0,4}}; CandidateStore invalid={bad,2,2};
    uint64_t before=rng.state;
    assert(bond_select(&invalid,8,&rng,1,&b) && !b.pairs && !b.bond_at_endpoint && rng.state==before);
    bad[1]=(CandidatePair){0,8}; assert(bond_select(&invalid,8,&rng,1,&b));
    bad[1]=(CandidatePair){0,1}; assert(bond_select(&invalid,8,&rng,1,&b));
    assert(bond_select(&empty,7,&rng,1,&b));
    /* 6分子・50seedの実際の生成結果で、型・距離・占有不変も確認。 */
    for(uint64_t seed=0;seed<50;++seed) {
        Lattice lat={0}; MoleculeStore m={0}; EndpointIndex index={0};
        CandidateStore candidates={0}; SearchStats stats={0};
        assert(!lattice_init(&lat,3,8) && !molecule_store_init(&m,6,5)); rng_seed(&rng,seed);
        for(size_t attempt=0;m.count<6 && attempt<10000;++attempt) {
            GenerationStatus st=molecule_try_generate(&lat,&m,&rng,
                (size_t)rng_bounded(&rng,lat.n_sites),m.count%2?1:0);
            assert(st==GEN_SUCCESS || st==GEN_TRAPPED || st==GEN_CENTER_OCCUPIED);
        }
        assert(m.count==6);
        int owners[512]; size_t paths[120];
        memcpy(owners,lat.owner,sizeof owners); memcpy(paths,m.path_sites,sizeof paths);
        assert(!endpoint_index_build(&lat,&m,&index));
        assert(!endpoint_find_candidates(&lat,&m,&index,&candidates,&stats));
        for(int multiple=0;multiple<=1;++multiple) {
            assert(!bond_select(&candidates,index.endpoint_count,&rng,multiple,&b));
            check(&candidates,&b,multiple);
            for(size_t k=0;k<b.count;++k) {
                size_t a=b.pairs[k].endpoint_a,other=b.pairs[k].endpoint_b;
                assert(m.molecules[a/4].type!=m.molecules[other/4].type);
                int ca[3],cb[3],distance=0;
                assert(!lattice_coord(&lat,paths[(a/4*4+a%4)*5+4],ca));
                assert(!lattice_coord(&lat,paths[(other/4*4+other%4)*5+4],cb));
                for(int d=0;d<3;++d) distance+=abs(ca[d]-cb[d]);
                assert(distance>=1 && distance<=2);
            }
            FILE *fp=tmpfile(); assert(fp && !bond_write_csv(fp,&lat,&m,&b)); fclose(fp);
            assert(!memcmp(owners,lat.owner,sizeof owners) && !memcmp(paths,m.path_sites,sizeof paths));
            bond_store_free(&b);
        }
        candidate_store_free(&candidates); endpoint_index_free(&index);
        molecule_store_free(&m); lattice_free(&lat);
    }
    puts("bond tests passed (constraints, reproducibility, rejection, 50 generated seeds)");
}
