#ifndef TPN_ENDPOINT_SEARCH_H
#define TPN_ENDPOINT_SEARCH_H
#include "molecule.h"
/*
 * endpoint ID = 4*molecule ID + arm ID。軌跡の最後の点を索引へ登録する。
 * 元の軌跡や占有格子を変更しない。生成完了後に構築し、分子を追加・変更
 * したら必ず再構築する。各格子点には高々1末端（点重複禁止が前提）。
 */
typedef struct {
    size_t *endpoint_at_site; /* SIZE_MAX = 末端なし。中心や腕途中の点も同じ値。 */
    size_t n_sites, endpoint_count;
} EndpointIndex;
typedef struct { size_t endpoint_a, endpoint_b; } CandidatePair;
typedef struct {
    CandidatePair *pairs;
    size_t count, capacity;
} CandidateStore;
typedef struct {
    size_t offsets_tested; /* 境界外も含め、確認した相対位置数 */
    size_t site_lookups;   /* 境界内で索引を参照した回数 */
} SearchStats;
/* outはゼロ初期化。エラーならoutを変更しない。成功後はfree関数で解放。 */
int endpoint_index_build(const Lattice *lat, const MoleculeStore *store, EndpointIndex *out);
void endpoint_index_free(EndpointIndex *index);
/* 3D・開境界・マンハッタン距離1～2。全末端を巡回し、A-Bの組だけ保存。
 * a<bに正規化して同じペアを2回保存しない。aは必ずしもA型ではない。
 * indexとstoreは同じ生成完了状態のものを渡す。結合は一切確定しない。
 */
int endpoint_find_candidates(const Lattice *lat, const MoleculeStore *store,
                             const EndpointIndex *index, CandidateStore *out,
                             SearchStats *stats);
void candidate_store_free(CandidateStore *store);
int candidate_write_csv(FILE *fp, const Lattice *lat, const MoleculeStore *store,
                        const CandidateStore *candidates);
#endif
