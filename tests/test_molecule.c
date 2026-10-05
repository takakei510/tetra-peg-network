#include "molecule.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
static void check_shape(const Lattice *lat,const MoleculeStore *s) {
    size_t occupied=0;
    for(size_t j=0;j<lat->n_sites;j++) occupied+=lat->owner[j]>=0;
    assert(occupied==s->count*(1+4*s->arm_length));
    for(size_t m=0;m<s->count;m++) {
        assert(lat->owner[s->molecules[m].center_site]==(int)m);
        for(size_t a=0;a<4;a++) {
            size_t previous=s->molecules[m].center_site;
            for(size_t t=0;t<s->arm_length;t++) {
                size_t k=molecule_path_index(s,m,a,t),site=s->path_sites[k];
                int c[3],d[3]; lattice_coord(lat,previous,c); lattice_coord(lat,site,d);
                int distance=0; for(int j=0;j<3;j++) distance+=abs(c[j]-d[j]);
                assert(distance==1 && lat->owner[site]==(int)m);
                unsigned direction=s->path_directions[k];
                assert(direction<6 && d[direction/2]-c[direction/2]==(direction%2==0?1:-1));
                previous=site;
            }
        }
    }
}
int main(void) {
    Lattice a={0},b={0}; MoleculeStore x={0},y={0}; Rng r,q;
    assert(!lattice_init(&a,3,16) && !lattice_init(&b,3,16));
    assert(!molecule_store_init(&x,2,5) && !molecule_store_init(&y,2,5));
    rng_seed(&r,42); rng_seed(&q,42);
    int center[3]={8,8,8}; size_t site; lattice_index(&a,center,&site);
    assert(molecule_try_generate(&a,&x,&r,site,1)==GEN_SUCCESS);
    assert(molecule_try_generate(&b,&y,&q,site,1)==GEN_SUCCESS);
    assert(x.molecules[0].type==TYPE_A && r.state==q.state);
    assert(!memcmp(x.path_sites,y.path_sites,20*sizeof(size_t)));
    assert(!memcmp(x.path_directions,y.path_directions,20));
    check_shape(&a,&x);
    assert(molecule_try_generate(&a,&x,&r,site,0)==GEN_CENTER_OCCUPIED);
    assert(x.count==1);
    /* 既存の1分子を残し、空点を2点だけ用意。1歩進んだ後に必ずtrapped。 */
    for(size_t j=0;j<a.n_sites;j++) if(a.owner[j]==-1) a.owner[j]=99;
    int c0[3]={0,0,0},c1[3]={1,0,0};size_t s0,s1;
    lattice_index(&a,c0,&s0);lattice_index(&a,c1,&s1);
    a.owner[s0]=a.owner[s1]=-1;
    int *before=malloc(a.n_sites*sizeof(int));assert(before);
    memcpy(before,a.owner,a.n_sites*sizeof(int));
    size_t old_paths[20];memcpy(old_paths,x.path_sites,sizeof old_paths);
    assert(molecule_try_generate(&a,&x,&r,s0,0)==GEN_TRAPPED);
    assert(x.count==1 && !memcmp(before,a.owner,a.n_sites*sizeof(int)));
    assert(!memcmp(old_paths,x.path_sites,sizeof old_paths));free(before);
    assert(molecule_try_generate(&a,&x,&r,s0,-1)==GEN_INVALID_CONFIG);
    MoleculeStore bad={0};assert(molecule_store_init(&bad,SIZE_MAX,5));
    FILE *fp=tmpfile();assert(fp && !molecule_write_csv(fp,&b,&y));fclose(fp);
    /* 1個目の腕上を中心に選んだ場合も、既存分子を変更せず棄却する。 */
    size_t first_paths[20]; unsigned char first_directions[20];
    memcpy(first_paths,y.path_sites,sizeof first_paths);
    memcpy(first_directions,y.path_directions,sizeof first_directions);
    Molecule first=y.molecules[0];
    assert(molecule_try_generate(&b,&y,&q,y.path_sites[0],0)==GEN_CENTER_OCCUPIED);
    /* 既存A分子を残した格子へ、B分子を追加する。失敗時は中心を再抽選。 */
    size_t attempts=0;
    while(y.count<2 && attempts++<10000) {
        GenerationStatus status=molecule_try_generate(&b,&y,&q,
            (size_t)rng_bounded(&q,b.n_sites),0);
        assert(status==GEN_SUCCESS || status==GEN_CENTER_OCCUPIED || status==GEN_TRAPPED);
    }
    assert(y.count==2 && y.molecules[1].type==TYPE_B);
    assert(y.molecules[0].type==first.type && y.molecules[0].center_site==first.center_site);
    assert(!memcmp(first_paths,y.path_sites,sizeof first_paths));
    assert(!memcmp(first_directions,y.path_directions,sizeof first_directions));
    check_shape(&b,&y); /* 42点・owner・最近接・全分子間の重複なし。 */
    molecule_store_free(&x);molecule_store_free(&y);lattice_free(&a);lattice_free(&b);
    puts("molecule tests passed");return 0;
}
