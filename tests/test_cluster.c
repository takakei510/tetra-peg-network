#include "cluster.h"
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static BondStore fixture(size_t n, const Bond *pairs, size_t count) {
    BondStore b={0}; b.endpoint_count=4*n; b.count=count;
    if(n) {
        b.bond_at_endpoint=malloc(4*n*sizeof(size_t)); assert(b.bond_at_endpoint);
        for(size_t a=0;a<4*n;++a) b.bond_at_endpoint[a]=SIZE_MAX;
    }
    if(count) {
        b.pairs=malloc(count*sizeof(Bond)); assert(b.pairs); memcpy(b.pairs,pairs,count*sizeof(Bond));
        for(size_t k=0;k<count;++k) {
            size_t a=pairs[k].endpoint_a,other=pairs[k].endpoint_b;
            assert(a<other && other<4*n && a/4!=other/4);
            assert(b.bond_at_endpoint[a]==SIZE_MAX && b.bond_at_endpoint[other]==SIZE_MAX);
            b.bond_at_endpoint[a]=k; b.bond_at_endpoint[other]=k;
        }
    }
    return b;
}
/* UFとは別のBFSで全分子の所属・サイズを照合する。直接結合だけに依存しない。 */
static void compare_bfs(const BondStore *b, size_t n) {
    ClusterStore s={0}; assert(!cluster_build(b,n,&s));
    unsigned char *adj=calloc(n*n+1,1); size_t *labels=malloc((n+1)*sizeof(size_t));
    size_t *queue=malloc((n+1)*sizeof(size_t)), *sizes=calloc(n+1,sizeof(size_t));
    size_t *degrees=calloc(n+1,sizeof(size_t)); assert(adj && labels && queue && sizes && degrees);
    for(size_t m=0;m<n;++m) labels[m]=SIZE_MAX;
    for(size_t k=0;k<b->count;++k) {
        size_t a=b->pairs[k].endpoint_a/4,other=b->pairs[k].endpoint_b/4;
        adj[a*n+other]=adj[other*n+a]=1; ++degrees[a]; ++degrees[other];
    }
    size_t components=0,isolated=0,largest=0,second=0;
    for(size_t first=0;first<n;++first) {
        isolated+=degrees[first]==0;
        if(labels[first]!=SIZE_MAX) continue;
        size_t head=0,tail=0; queue[tail++]=first; labels[first]=components;
        while(head<tail) {
            size_t m=queue[head++]; ++sizes[components];
            for(size_t next=0;next<n;++next) if(adj[m*n+next] && labels[next]==SIZE_MAX) {
                labels[next]=components; queue[tail++]=next;
            }
        }
        ++components;
    }
    size_t molecules_sum=0,bonds_sum=0;
    assert(s.cluster_count==components && s.isolated_molecules==isolated && s.bond_count==b->count);
    for(size_t id=0;id<components;++id) {
        assert(s.cluster_sizes[id]==sizes[id]); molecules_sum+=s.cluster_sizes[id]; bonds_sum+=s.cluster_bonds[id];
        if(sizes[id]>largest) {second=largest;largest=sizes[id];}
        else if(sizes[id]>second) second=sizes[id];
    }
    assert(molecules_sum==n && bonds_sum==b->count && s.largest==largest && s.second_largest==second);
    for(size_t m=0;m<n;++m) assert(s.molecule_cluster[m]==labels[m] && s.bond_degree[m]==degrees[m]);
    FILE *fp=tmpfile(); assert(fp && !cluster_members_write_csv(fp,&s)); fclose(fp);
    fp=tmpfile(); assert(fp && !cluster_summary_write_csv(fp,&s)); fclose(fp);
    fp=tmpfile(); assert(fp && !network_summary_write_csv(fp,&s));
    if(!n) {rewind(fp); char line[256];assert(fgets(line,sizeof line,fp));assert(fgets(line,sizeof line,fp));assert(!strcmp(line,"0,0,0,0,0,0,\n"));}
    fclose(fp);
    cluster_store_free(&s); free(adj); free(labels); free(queue); free(sizes); free(degrees);
}

int main(void) {
    /* 4分子の鎖＋同じ分子対の2本目＋2孤立分子。クラスタは4,1,1。 */
    const Bond pairs[]={{0,4},{5,8},{9,12},{1,6}};
    BondStore b=fixture(6,pairs,4); ClusterStore s={0},reverse={0};
    compare_bfs(&b,6); assert(!cluster_build(&b,6,&s));
    assert(s.cluster_count==3 && s.largest==4 && s.second_largest==1 && s.isolated_molecules==2);
    assert(s.cluster_bonds[0]==4 && s.bond_degree[0]==2 && s.bond_degree[1]==3);
    const Bond reversed[]={{1,6},{9,12},{5,8},{0,4}};
    BondStore other=fixture(6,reversed,4); assert(!cluster_build(&other,6,&reverse));
    assert(!memcmp(s.molecule_cluster,reverse.molecule_cluster,6*sizeof(size_t)));
    assert(!memcmp(s.cluster_sizes,reverse.cluster_sizes,3*sizeof(size_t)));
    cluster_store_free(&s); cluster_store_free(&reverse); bond_store_free(&other);
    /* 壊れた逆引き、不正末端ID、分子数不一致を拒否する。 */
    b.bond_at_endpoint[0]=SIZE_MAX; assert(cluster_build(&b,6,&s) && !s.molecule_cluster);
    b.bond_at_endpoint[0]=0; b.pairs[0].endpoint_b=24; assert(cluster_build(&b,6,&s));
    b.pairs[0].endpoint_b=4; assert(cluster_build(&b,5,&s)); bond_store_free(&b);
    const Bond cycle[]={{0,4},{5,8},{9,12},{1,13}}; b=fixture(4,cycle,4);
    compare_bfs(&b,4); bond_store_free(&b); /* 閉路も所属・サイズを増やさない。 */
    b=fixture(4,NULL,0);compare_bfs(&b,4);bond_store_free(&b);
    const Bond ties[]={{0,4},{8,12}};b=fixture(4,ties,2);compare_bfs(&b,4);bond_store_free(&b);
    b=fixture(0,NULL,0);compare_bfs(&b,0);bond_store_free(&b);
    /* 実際の6分子生成・候補探索・抽選の50seedについてもBFSと一致。 */
    for(uint64_t seed=0;seed<50;++seed) {
        Lattice lat={0}; MoleculeStore m={0}; EndpointIndex index={0};
        CandidateStore candidates={0}; SearchStats stats={0}; Rng rng;
        assert(!lattice_init(&lat,3,8) && !molecule_store_init(&m,6,5)); rng_seed(&rng,seed);
        for(size_t attempt=0;m.count<6 && attempt<10000;++attempt) {
            GenerationStatus st=molecule_try_generate(&lat,&m,&rng,
                (size_t)rng_bounded(&rng,lat.n_sites),m.count%2?1:0);
            assert(st==GEN_SUCCESS || st==GEN_TRAPPED || st==GEN_CENTER_OCCUPIED);
        }
        assert(m.count==6 && !endpoint_index_build(&lat,&m,&index));
        assert(!endpoint_find_candidates(&lat,&m,&index,&candidates,&stats));
        for(int multiple=0;multiple<=1;++multiple) {
            assert(!bond_select(&candidates,24,&rng,multiple,&b));
            size_t before[24]; Bond saved[12]; memcpy(before,b.bond_at_endpoint,sizeof before);
            if(b.count) memcpy(saved,b.pairs,b.count*sizeof(Bond));
            compare_bfs(&b,6);
            assert(!memcmp(before,b.bond_at_endpoint,sizeof before));
            assert(!b.count || !memcmp(saved,b.pairs,b.count*sizeof(Bond)));
            bond_store_free(&b);
        }
        candidate_store_free(&candidates);endpoint_index_free(&index);
        molecule_store_free(&m);lattice_free(&lat);
    }
    puts("cluster tests passed (synthetic cases + 50 generated seeds vs BFS)");
}
