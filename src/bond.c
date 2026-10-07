#include "bond.h"
#include <stdlib.h>
#include <stdint.h>

void bond_store_free(BondStore *s) {
    if(s) { free(s->pairs); free(s->bond_at_endpoint); *s=(BondStore){0}; }
}

/* 同じ分子対の結合を禁止する場合だけ、最大4腕の既存相手を調べる。 */
static int eligible(const BondStore *s, size_t a, size_t b, int multiple) {
    if(s->bond_at_endpoint[b]!=SIZE_MAX) return 0;
    if(multiple) return 1;
    size_t molecule=a/ARM_COUNT;
    for(size_t arm=0;arm<ARM_COUNT;++arm) {
        size_t id=ARM_COUNT*molecule+arm;
        size_t bond=s->bond_at_endpoint[id];
        if(bond==SIZE_MAX) continue;
        Bond pair=s->pairs[bond];
        size_t partner=pair.endpoint_a==id ? pair.endpoint_b : pair.endpoint_a;
        if(partner/ARM_COUNT==b/ARM_COUNT) return 0;
    }
    return 1;
}

int bond_select(const CandidateStore *c, size_t n, Rng *rng, int multiple, BondStore *out) {
    if(!c || !rng || !out || out->pairs || out->bond_at_endpoint || out->count ||
       out->endpoint_count || n%ARM_COUNT || (multiple!=0 && multiple!=1) ||
       c->count>c->capacity || (c->count && !c->pairs) ||
       n>=SIZE_MAX/sizeof(size_t) || c->count>SIZE_MAX/2/sizeof(size_t) ||
       n/2>SIZE_MAX/sizeof(Bond)) return -1;
    BondStore s={0}; s.endpoint_count=n;
    /* 隣接リストを連続配列にまとめる。毎回全候補を走査しない。 */
    size_t *offset=calloc(n+1,sizeof *offset), *cursor=NULL, *neighbors=NULL, *order=NULL;
    if(!offset) goto fail;
    if(n) {
        cursor=malloc(n*sizeof *cursor); order=malloc(n*sizeof *order);
        s.bond_at_endpoint=malloc(n*sizeof *s.bond_at_endpoint);
        s.pairs=malloc((n/2)*sizeof *s.pairs); /* 1末端1結合なので最大n/2本。 */
        if(!cursor || !order || !s.bond_at_endpoint || !s.pairs) goto fail;
    }
    if(c->count) {
        neighbors=malloc(2*c->count*sizeof *neighbors);
        if(!neighbors) goto fail;
    }
    for(size_t k=0;k<c->count;++k) {
        size_t a=c->pairs[k].endpoint_a,b=c->pairs[k].endpoint_b;
        /* 末端ID順・分子間・24近傍の上限を確認してから添字に使う。 */
        if(a>=b || b>=n || a/ARM_COUNT==b/ARM_COUNT) goto fail;
        if(++offset[a+1]>24 || ++offset[b+1]>24) goto fail;
    }
    for(size_t a=0;a<n;++a) offset[a+1]+=offset[a];
    for(size_t a=0;a<n;++a) {
        cursor[a]=offset[a]; order[a]=a; s.bond_at_endpoint[a]=SIZE_MAX;
    }
    for(size_t k=0;k<c->count;++k) {
        size_t a=c->pairs[k].endpoint_a,b=c->pairs[k].endpoint_b;
        /* 同じ候補の重複は抽選確率を変えるため拒否する。最大24件を照合。 */
        for(size_t j=offset[a];j<cursor[a];++j) if(neighbors[j]==b) goto fail;
        neighbors[cursor[a]++]=b; neighbors[cursor[b]++]=a;
    }
    /* 必要な確保と入力検証を済ませてから乱数を消費する。 */
    for(size_t remaining=n;remaining>1;--remaining) {
        size_t j=(size_t)rng_bounded(rng,remaining);
        size_t saved=order[remaining-1]; order[remaining-1]=order[j]; order[j]=saved;
    }
    for(size_t k=0;k<n;++k) {
        size_t a=order[k];
        if(s.bond_at_endpoint[a]!=SIZE_MAX) continue;
        /* 既に相手として結ばれた末端は除く。空いている候補数を数える。 */
        size_t available=0;
        for(size_t j=offset[a];j<offset[a+1];++j)
            available+=eligible(&s,a,neighbors[j],multiple);
        if(!available) continue;
        size_t selected=(size_t)rng_bounded(rng,available), b=SIZE_MAX;
        for(size_t j=offset[a];j<offset[a+1];++j) {
            if(!eligible(&s,a,neighbors[j],multiple)) continue;
            if(selected==0) {b=neighbors[j];break;}
            --selected;
        }
        /* 両末端を同時に登録。以後はどちらから調べても結合済みになる。 */
        s.pairs[s.count]=(Bond){a<b?a:b,a<b?b:a};
        s.bond_at_endpoint[a]=s.count; s.bond_at_endpoint[b]=s.count; ++s.count;
    }
    free(offset); free(cursor); free(neighbors); free(order); *out=s;
    return 0;
fail:
    free(offset); free(cursor); free(neighbors); free(order); bond_store_free(&s);
    return -1;
}

int bond_write_csv(FILE *fp, const Lattice *lat, const MoleculeStore *m, const BondStore *s) {
    if(!fp || !lat || !m || !s || !m->arm_length || m->count>SIZE_MAX/ARM_COUNT ||
       s->endpoint_count!=ARM_COUNT*m->count || (s->count && !s->pairs)) return -1;
    if(fputs("bond_id,endpoint_a,molecule_a,arm_a,type_a,site_a,x_a,y_a,z_a,endpoint_b,molecule_b,arm_b,type_b,site_b,x_b,y_b,z_b,manhattan_distance\n",fp)==EOF) return -1;
    for(size_t k=0;k<s->count;++k) {
        size_t a=s->pairs[k].endpoint_a,b=s->pairs[k].endpoint_b;
        if(a>=b || b>=s->endpoint_count) return -1;
        size_t sa=m->path_sites[molecule_path_index(m,a/4,a%4,m->arm_length-1)];
        size_t sb=m->path_sites[molecule_path_index(m,b/4,b%4,m->arm_length-1)];
        int ca[3],cb[3];
        if(lattice_coord(lat,sa,ca) || lattice_coord(lat,sb,cb)) return -1;
        long long distance=0; /* 出力時だけ、成立したペアの距離を記録。 */
        for(int d=0;d<3;++d) distance+=llabs((long long)ca[d]-cb[d]);
        if(fprintf(fp,"%zu,%zu,%zu,%zu,%c,%zu,%d,%d,%d,%zu,%zu,%zu,%c,%zu,%d,%d,%d,%lld\n",
            k,a,a/4,a%4,m->molecules[a/4].type==TYPE_A?'A':'B',sa,ca[0],ca[1],ca[2],
            b,b/4,b%4,m->molecules[b/4].type==TYPE_A?'A':'B',sb,cb[0],cb[1],cb[2],distance)<0) return -1;
    }
    return ferror(fp)?-1:0;
}
