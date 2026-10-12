#ifndef TPN_CLUSTER_H
#define TPN_CLUSTER_H
#include "bond.h"
/* 分子を頂点、成立した末端結合を辺として求める派生データ。 */
typedef struct {
    size_t molecule_count, bond_count, cluster_count;
    size_t largest, second_largest, isolated_molecules; /* 単位は分子数。 */
    size_t *molecule_cluster; /* 分子ID → 出力用クラスタID。 */
    size_t *cluster_sizes;    /* クラスタID → 所属分子数。 */
    size_t *cluster_bonds;    /* 多重結合も1本ずつ数えた、クラスタ内の結合数。 */
    size_t *bond_degree;      /* 分子ID → 結合した腕の本数。隣接分子数とは別。 */
} ClusterStore;
/* outはゼロ初期化。結合と逆引きの整合性を検査し、エラー時はout不変。
 * 孤立分子も含める。クラスタIDは最小所属分子IDの順に0から付ける。
 * BondStoreが正本。UFの内部root番号はCSVのクラスタIDに使わない。
 */
int cluster_build(const BondStore *bonds, size_t molecule_count, ClusterStore *out);
void cluster_store_free(ClusterStore *store);
int cluster_members_write_csv(FILE *fp, const ClusterStore *store);
int cluster_summary_write_csv(FILE *fp, const ClusterStore *store);
int network_summary_write_csv(FILE *fp, const ClusterStore *store);
#endif
