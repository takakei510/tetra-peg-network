#ifndef TPN_BOND_H
#define TPN_BOND_H
#include "endpoint_search.h"
/* 結合の正本は末端ペア。座標・型はMoleculeStoreから取得する。 */
typedef CandidatePair Bond;
typedef struct {
    Bond *pairs;                 /* 添字がbond ID。成立順に保存。 */
    size_t count, endpoint_count;
    size_t *bond_at_endpoint;     /* endpoint ID → bond ID。SIZE_MAXは未結合。 */
} BondStore;
/* ③-Aで得た重複なしA-B候補を渡す。入力の候補・軌跡は変更しない。
 * 全末端をシャッフルし、未結合の候補から相手を一様に1個選ぶ。
 * allow_multiple_bonds=1なら同じ分子対の別腕結合も許可、0なら禁止。
 * outはゼロ初期化。エラー時はoutと乱数状態を変更しない。
 * 候補の距離・型の判定はendpoint_find_candidatesの責務。
 */
int bond_select(const CandidateStore *candidates, size_t endpoint_count,
                Rng *rng, int allow_multiple_bonds, BondStore *out);
void bond_store_free(BondStore *store);
int bond_write_csv(FILE *fp, const Lattice *lat, const MoleculeStore *molecules,
                   const BondStore *bonds);
#endif
