/*
 * 完成分子の保存と4-arm kinetic SAW生成。
 * 旧percolation-model/src/random_walk.cのrun_one_walk()にある
 * 「近傍列挙→訪問点除外→候補から選択→候補ゼロでtrapped」を参考に再実装。
 * 旧関数の直接移植ではない。訪問判定を全分子共通ownerへ拡張し、
 * MSD集計を分離、4腕の軌跡保存と分子全体の失敗復元を追加した。
 * 先生の生成規則との一致は未確認。平衡SAWの一様抽出ではない。
 */
#include "molecule.h"
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
int molecule_store_init(MoleculeStore *s, size_t n, size_t l) {
    if (!s || !n || !l || n > (size_t)INT_MAX || l > SIZE_MAX/ARM_COUNT)
        return -1;
    size_t per = ARM_COUNT*l;
    if (n > SIZE_MAX/per || n > SIZE_MAX/sizeof(Molecule)) return -1;
    size_t points=n*per;
    if (points > SIZE_MAX/sizeof(size_t)) return -1;
    MoleculeStore tmp={0};
    tmp.molecules=malloc(n*sizeof *tmp.molecules);
    tmp.path_sites=malloc(points*sizeof *tmp.path_sites);
    tmp.path_directions=malloc(points*sizeof *tmp.path_directions);
    if (!tmp.molecules || !tmp.path_sites || !tmp.path_directions) {
        molecule_store_free(&tmp); return -1;
    }
    tmp.capacity=n; tmp.arm_length=l; *s=tmp; return 0;
}
void molecule_store_free(MoleculeStore *s) {
    if (!s) return;
    free(s->molecules); free(s->path_sites); free(s->path_directions);
    *s=(MoleculeStore){0};
}
size_t molecule_path_index(const MoleculeStore *s,size_t m,size_t a,size_t t) {
    return (m*ARM_COUNT+a)*s->arm_length+t;
}
GenerationStatus molecule_try_generate(Lattice *lat,MoleculeStore *s,
                                       Rng *rng,size_t center,double r) {
    if (!lat || !s || !rng || lat->dim!=3 || !lat->owner ||
        center>=lat->n_sites || !isfinite(r) || r<0 || r>1 ||
        !s->arm_length || s->count>=s->capacity || !s->molecules ||
        !s->path_sites || !s->path_directions) return GEN_INVALID_CONFIG;
    if (lat->owner[center]!=-1) return GEN_CENTER_OCCUPIED;
    /* 仮占有一覧を先に確保。以降、成功確定まで追加のmallocは行わない。 */
    size_t points=ARM_COUNT*s->arm_length;
    if (points==SIZE_MAX || points+1>SIZE_MAX/sizeof(size_t))
        return GEN_INVALID_CONFIG;
    size_t *touched=malloc((points+1)*sizeof *touched);
    if (!touched) return GEN_ALLOCATION_ERROR;
    size_t used=0, base=s->count*points;
    Molecule molecule={rng_uniform(rng)<r ? TYPE_A:TYPE_B,center};
    lat->owner[center]=(int)s->count; touched[used++]=center;
    for (size_t a=0;a<ARM_COUNT;a++) {
        size_t current=center;
        for (size_t t=0;t<s->arm_length;t++) {
            size_t candidates[6]; unsigned char directions[6]; size_t n=0;
            int xyz[3]; lattice_coord(lat,current,xyz);
            /* 列挙順は再現性の一部。+x,-x,+y,-y,+z,-zで固定する。 */
            for (unsigned char d=0;d<6;d++) {
                int c[3]={xyz[0],xyz[1],xyz[2]};
                c[d/2]+=(d%2==0 ? 1:-1);
                size_t next;
                if (!lattice_index(lat,c,&next) && lat->owner[next]==-1) {
                    candidates[n]=next; directions[n++]=d;
                }
            }
            if (!n) {
                /* この試行の点だけ戻す。既存分子のownerは変更しない。 */
                for (size_t j=0;j<used;j++) lat->owner[touched[j]]=-1;
                free(touched); return GEN_TRAPPED;
            }
            size_t choice=(size_t)rng_bounded(rng,n);
            current=candidates[choice];
            lat->owner[current]=(int)s->count; touched[used++]=current;
            size_t index=base+a*s->arm_length+t;
            /* 未登録スロットを作業領域として使用。count未満は触らない。 */
            s->path_sites[index]=current;
            s->path_directions[index]=directions[choice];
        }
    }
    s->molecules[s->count]=molecule;
    ++s->count; /* 最後に公開することで、不完全分子を登録しない。 */
    free(touched); return GEN_SUCCESS;
}
int molecule_write_csv(FILE *fp,const Lattice *lat,const MoleculeStore *s) {
    if (!fp || !lat || !s) return -1;
    if (fprintf(fp,"molecule,type,arm,step,site,x,y,z,direction\n")<0) return -1;
    for (size_t m=0;m<s->count;m++) {
        int c[3]; lattice_coord(lat,s->molecules[m].center_site,c);
        if (fprintf(fp,"%zu,%c,-1,0,%zu,%d,%d,%d,-1\n",m,
            s->molecules[m].type==TYPE_A?'A':'B',s->molecules[m].center_site,
            c[0],c[1],c[2])<0) return -1;
        for (size_t a=0;a<ARM_COUNT;a++) for (size_t t=0;t<s->arm_length;t++) {
            size_t k=molecule_path_index(s,m,a,t), site=s->path_sites[k];
            lattice_coord(lat,site,c);
            if (fprintf(fp,"%zu,%c,%zu,%zu,%zu,%d,%d,%d,%u\n",m,
                s->molecules[m].type==TYPE_A?'A':'B',a,t+1,site,
                c[0],c[1],c[2],(unsigned)s->path_directions[k])<0) return -1;
        }
    }
    return ferror(fp)?-1:0;
}
