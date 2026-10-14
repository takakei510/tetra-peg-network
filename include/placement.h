#ifndef TPN_PLACEMENT_H
#define TPN_PLACEMENT_H
#include "network_config.h"
#include "molecule.h"
/* 試行上限は全分子の合計。占有中心の棄却と腕の行き止まりを分ける。 */
typedef struct {
    size_t attempts, center_occupied, trapped;
    GenerationStatus last_status;
} PlacementStats;
/* 空の格子・storeから開始。0=目標達成、1=試行上限で未達、-1=実行エラー。
 * 未達でも完成分子を残す。失敗途中の形はmolecule_try_generateで取り消す。
 * 乱数は呼出元所有。後続の結合に同じ系列を引き継ぐ。
 */
int placement_generate(Lattice *lat, MoleculeStore *store, Rng *rng,
                       const NetworkConfig *config, size_t max_attempts, PlacementStats *stats);
#endif
