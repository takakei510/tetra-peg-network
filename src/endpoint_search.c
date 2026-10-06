#include "endpoint_search.h"
#include <stdlib.h>
#include <stdint.h>

static int valid_store(const Lattice *lat, const MoleculeStore *s) {
    return lat && lat->owner && lat->dim == 3 && s && s->arm_length &&
        s->molecules && s->path_sites && s->count <= s->capacity &&
        s->count <= SIZE_MAX / (ARM_COUNT * 24) &&
        s->capacity <= SIZE_MAX / ARM_COUNT / s->arm_length;
}

void endpoint_index_free(EndpointIndex *i) {
    if(i) { free(i->endpoint_at_site); *i = (EndpointIndex){0}; }
}
void candidate_store_free(CandidateStore *s) {
    if(s) { free(s->pairs); *s = (CandidateStore){0}; }
}

int endpoint_index_build(const Lattice *lat, const MoleculeStore *s, EndpointIndex *out) {
    if(!out || out->endpoint_at_site || !valid_store(lat,s) ||
       lat->n_sites > SIZE_MAX / sizeof(size_t)) return -1;
    EndpointIndex index = {0};
    index.endpoint_at_site = malloc(lat->n_sites * sizeof(size_t));
    if(!index.endpoint_at_site) return -1;
    index.n_sites = lat->n_sites;
    index.endpoint_count = ARM_COUNT * s->count;
    /* 密な配列は初期化O(L^3)、メモリO(L^3)。その後の位置照会はO(1)。 */
    for(size_t site=0; site<lat->n_sites; ++site) index.endpoint_at_site[site]=SIZE_MAX;
    for(size_t m=0; m<s->count; ++m) {
        if(s->molecules[m].type!=TYPE_A && s->molecules[m].type!=TYPE_B) goto fail;
        for(size_t arm=0; arm<ARM_COUNT; ++arm) {
            size_t site=s->path_sites[molecule_path_index(s,m,arm,s->arm_length-1)];
            if(site>=lat->n_sites || index.endpoint_at_site[site]!=SIZE_MAX ||
               lat->owner[site]<0 || (size_t)lat->owner[site]!=m) goto fail;
            index.endpoint_at_site[site]=ARM_COUNT*m+arm;
        }
    }
    *out=index;
    return 0;
fail:
    endpoint_index_free(&index);
    return -1;
}

static int append(CandidateStore *s, size_t a, size_t b) {
    if(s->count==s->capacity) {
        size_t limit=SIZE_MAX/sizeof(CandidatePair);
        size_t next=s->capacity ? s->capacity*2 : 16;
        if(s->capacity>limit/2 || next>limit) return -1;
        CandidatePair *pairs=realloc(s->pairs,next*sizeof *pairs);
        if(!pairs) return -1;
        s->pairs=pairs; s->capacity=next;
    }
    s->pairs[s->count++]=(CandidatePair){a,b};
    return 0;
}

int endpoint_find_candidates(const Lattice *lat, const MoleculeStore *s,
                             const EndpointIndex *index, CandidateStore *out,
                             SearchStats *stats) {
    if(!out || out->pairs || !stats || !valid_store(lat,s) || !index ||
       !index->endpoint_at_site || index->n_sites!=lat->n_sites ||
       index->endpoint_count!=ARM_COUNT*s->count) return -1;
    /*
     * 24点の座標差を探索の開始時に1回だけ作る。末端ごとに全末端との距離を
     * 計算しない。半径2は現段階の固定仕様。変更時は距離規則も記録する。
     */
    int offsets[24][3], n=0;
    for(int x=-2;x<=2;++x) for(int y=-2;y<=2;++y) for(int z=-2;z<=2;++z) {
        int distance=abs(x)+abs(y)+abs(z);
        if(distance>=1 && distance<=2) {
            offsets[n][0]=x; offsets[n][1]=y; offsets[n++][2]=z;
        }
    }
    CandidateStore candidates={0}; SearchStats measured={0};
    for(size_t a=0;a<index->endpoint_count;++a) {
        size_t site=s->path_sites[molecule_path_index(s,a/4,a%4,s->arm_length-1)];
        int c[3];
        if(lattice_coord(lat,site,c) || index->endpoint_at_site[site]!=a) goto fail;
        for(int k=0;k<n;++k) {
            ++measured.offsets_tested;
            int near[3]; int inside=1;
            for(int d=0;d<3;++d) {
                /* 加算前に境界を調べ、int上限付近でも符号付き加算を溢れさせない。 */
                int delta=offsets[k][d];
                if((delta<0 && c[d]<-delta) ||
                   (delta>0 && c[d]>lat->L-1-delta)) {inside=0;break;}
                near[d]=c[d]+delta;
            }
            if(!inside) continue;
            size_t neighbor;
            if(lattice_index(lat,near,&neighbor)) goto fail;
            ++measured.site_lookups;
            size_t b=index->endpoint_at_site[neighbor];
            if(b==SIZE_MAX) continue;
            if(b>=index->endpoint_count) goto fail;
            /* 両側から発見するためID順で片方だけ採用。分子内/同型は除外。 */
            if(a>=b || s->molecules[a/4].type==s->molecules[b/4].type) continue;
            if(append(&candidates,a,b)) goto fail;
        }
    }
    *out=candidates; *stats=measured;
    return 0;
fail:
    candidate_store_free(&candidates);
    return -1;
}

int candidate_write_csv(FILE *fp, const Lattice *lat, const MoleculeStore *s,
                        const CandidateStore *candidates) {
    if(!fp || !candidates || !valid_store(lat,s)) return -1;
    if(fputs("endpoint_a,molecule_a,arm_a,type_a,site_a,x_a,y_a,z_a,endpoint_b,molecule_b,arm_b,type_b,site_b,x_b,y_b,z_b,manhattan_distance\n",fp)==EOF) return -1;
    for(size_t k=0;k<candidates->count;++k) {
        size_t a=candidates->pairs[k].endpoint_a,b=candidates->pairs[k].endpoint_b;
        if(a>=4*s->count || b>=4*s->count) return -1;
        size_t sa=s->path_sites[molecule_path_index(s,a/4,a%4,s->arm_length-1)];
        size_t sb=s->path_sites[molecule_path_index(s,b/4,b%4,s->arm_length-1)];
        int ca[3],cb[3];
        if(lattice_coord(lat,sa,ca) || lattice_coord(lat,sb,cb)) return -1;
        int distance=abs(ca[0]-cb[0])+abs(ca[1]-cb[1])+abs(ca[2]-cb[2]);
        if(fprintf(fp,"%zu,%zu,%zu,%c,%zu,%d,%d,%d,%zu,%zu,%zu,%c,%zu,%d,%d,%d,%d\n",
            a,a/4,a%4,s->molecules[a/4].type==TYPE_A?'A':'B',sa,ca[0],ca[1],ca[2],
            b,b/4,b%4,s->molecules[b/4].type==TYPE_A?'A':'B',sb,cb[0],cb[1],cb[2],distance)<0) return -1;
    }
    return ferror(fp) ? -1 : 0;
}
