#ifndef TPN_MOLECULE_H
#define TPN_MOLECULE_H
#include "lattice.h"
#include "rng.h"
#include <stdio.h>
enum { ARM_COUNT = 4 };
typedef enum { TYPE_A, TYPE_B } MoleculeType;
/* 中心は軌跡に含めない。方向は +x,-x,+y,-y,+z,-z を0..5で保存。 */
typedef struct { MoleculeType type; size_t center_site; } Molecule;
typedef struct {
    Molecule *molecules;
    size_t *path_sites;
    unsigned char *path_directions;
    size_t count, capacity, arm_length;
} MoleculeStore;
typedef enum { GEN_SUCCESS, GEN_CENTER_OCCUPIED, GEN_TRAPPED,
               GEN_INVALID_CONFIG, GEN_ALLOCATION_ERROR } GenerationStatus;
/* 空のstoreを初期化する。容量と腕長は正、サイズの積を検査する。 */
int molecule_store_init(MoleculeStore *store, size_t capacity, size_t arm_length);
void molecule_store_free(MoleculeStore *store);
/*
 * 指定中心で1回試行。3D・開境界・空近傍から一様選択・腕0..3の順。
 * 中心抽選は呼び出し側で行うため、固定中心での検証にも使用できる。
 * 成功時だけcountを増加。棄却時は既存占有と登録データを保持する。
 * A/Bは空中心確認後に抽選。乱数状態は棄却時も巻き戻さない。
 */
GenerationStatus molecule_try_generate(Lattice *lat, MoleculeStore *store,
                                      Rng *rng, size_t center, double a_probability);
/* index=(molecule*4+arm)*l+step。最後のstepが末端。 */
size_t molecule_path_index(const MoleculeStore *store, size_t molecule,
                           size_t arm, size_t step);
/* 全登録分子の中心と軌跡をCSVへ。step=0は中心、1..lは腕の到達点。 */
int molecule_write_csv(FILE *fp, const Lattice *lat, const MoleculeStore *store);
#endif
