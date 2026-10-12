#include "cluster.h"
#include "union_find.h"
#include <stdint.h>
#include <stdlib.h>

void cluster_store_free(ClusterStore *s) {
    if(!s) return;
    free(s->molecule_cluster); free(s->cluster_sizes);
    free(s->cluster_bonds); free(s->bond_degree); *s=(ClusterStore){0};
}

int cluster_build(const BondStore *b, size_t n, ClusterStore *out) {
    if(!b || !out || out->molecule_cluster || out->cluster_sizes || out->cluster_bonds ||
       out->bond_degree || out->molecule_count || out->bond_count || out->cluster_count ||
       out->largest || out->second_largest || out->isolated_molecules ||
       n>SIZE_MAX/ARM_COUNT || n>SIZE_MAX/sizeof(size_t) ||
       b->endpoint_count!=ARM_COUNT*n || b->count>b->endpoint_count/2 ||
       (b->count && !b->pairs) || (n && !b->bond_at_endpoint)) return -1;
    ClusterStore s={0}; s.molecule_count=n; s.bond_count=b->count;
    /* 0分子は空ネットワーク。正のnを要求する既存uf_initは呼ばない。 */
    if(!n) { *out=s; return 0; }
    UnionFind uf={0}; size_t *root_cluster=NULL;
    s.molecule_cluster=malloc(n*sizeof *s.molecule_cluster);
    s.cluster_sizes=calloc(n,sizeof *s.cluster_sizes);
    s.cluster_bonds=calloc(n,sizeof *s.cluster_bonds);
    s.bond_degree=calloc(n,sizeof *s.bond_degree);
    root_cluster=malloc(n*sizeof *root_cluster);
    /* 全登録分子をサイズ1のクラスタとして開始。孤立分子を落とさない。 */
    if(!s.molecule_cluster || !s.cluster_sizes || !s.cluster_bonds || !s.bond_degree ||
       !root_cluster || uf_init(&uf,n)) goto fail;
    for(size_t m=0;m<n;++m) root_cluster[m]=SIZE_MAX;
    for(size_t k=0;k<b->count;++k) {
        size_t a=b->pairs[k].endpoint_a,other=b->pairs[k].endpoint_b;
        /* 1末端1結合の逆引きと照合。不正なIDをUFへ渡さない。 */
        if(a>=other || other>=b->endpoint_count || a/4==other/4 ||
           b->bond_at_endpoint[a]!=k || b->bond_at_endpoint[other]!=k) goto fail;
        ++s.bond_degree[a/4]; ++s.bond_degree[other/4];
        /* UFの要素は末端・格子点ではなく「分子」。同じrootなら統合しない。 */
        if(uf_union(&uf,a/4,other/4)<0) goto fail;
    }
    for(size_t a=0;a<b->endpoint_count;++a) {
        size_t k=b->bond_at_endpoint[a];
        if(k==SIZE_MAX) continue;
        if(k>=b->count || (b->pairs[k].endpoint_a!=a && b->pairs[k].endpoint_b!=a)) goto fail;
    }
    s.cluster_count=uf.components;
    uf_largest_two(&uf,&s.largest,&s.second_largest); /* 既存のサイズ集計も直接利用。 */
    size_t next=0;
    for(size_t m=0;m<n;++m) {
        size_t root=uf_find(&uf,m); /* 間接的な接続も同じ代表分子へたどる。 */
        /* rootは統合順で変わる。出力IDは最小所属分子の順に付け直す。 */
        if(root_cluster[root]==SIZE_MAX) root_cluster[root]=next++;
        size_t id=root_cluster[root]; s.molecule_cluster[m]=id; ++s.cluster_sizes[id];
        s.isolated_molecules+=s.bond_degree[m]==0;
    }
    for(size_t k=0;k<b->count;++k)
        ++s.cluster_bonds[s.molecule_cluster[b->pairs[k].endpoint_a/4]];
    free(root_cluster); uf_free(&uf); *out=s;
    return 0;
fail:
    free(root_cluster); uf_free(&uf); cluster_store_free(&s);
    return -1;
}

int cluster_members_write_csv(FILE *fp, const ClusterStore *s) {
    if(!fp || !s || (s->molecule_count && (!s->molecule_cluster || !s->cluster_sizes || !s->bond_degree))) return -1;
    if(fputs("molecule,cluster_id,cluster_size,bond_degree\n",fp)==EOF) return -1;
    for(size_t m=0;m<s->molecule_count;++m) {
        size_t id=s->molecule_cluster[m]; if(id>=s->cluster_count) return -1;
        if(fprintf(fp,"%zu,%zu,%zu,%zu\n",m,id,s->cluster_sizes[id],s->bond_degree[m])<0) return -1;
    }
    return ferror(fp)?-1:0;
}
int cluster_summary_write_csv(FILE *fp, const ClusterStore *s) {
    if(!fp || !s || (s->cluster_count && (!s->cluster_sizes || !s->cluster_bonds))) return -1;
    if(fputs("cluster_id,molecule_count,bond_count\n",fp)==EOF) return -1;
    for(size_t id=0;id<s->cluster_count;++id)
        if(fprintf(fp,"%zu,%zu,%zu\n",id,s->cluster_sizes[id],s->cluster_bonds[id])<0) return -1;
    return ferror(fp)?-1:0;
}
int network_summary_write_csv(FILE *fp, const ClusterStore *s) {
    if(!fp || !s) return -1;
    if(fputs("molecules,bonds,clusters,largest,second_largest,isolated_molecules,largest_fraction\n",fp)==EOF) return -1;
    if(fprintf(fp,"%zu,%zu,%zu,%zu,%zu,%zu,",s->molecule_count,s->bond_count,s->cluster_count,
        s->largest,s->second_largest,s->isolated_molecules)<0) return -1;
    /* N=0で最大クラスタ比は未定義なので、0ではなく空欄を出す。 */
    if(s->molecule_count && fprintf(fp,"%.17g",(double)s->largest/s->molecule_count)<0) return -1;
    return fputc('\n',fp)==EOF || ferror(fp)?-1:0;
}
