#ifndef TPN_NETWORK_CONFIG_H
#define TPN_NETWORK_CONFIG_H
#include "config.h"
/* 多分子用。既存の3キーConfigを含め、生成条件を別に管理する。 */
typedef enum { TARGET_MOLECULE_COUNT, TARGET_OCCUPIED_FRACTION } NetworkTargetMode;
typedef struct {
    Config base;             /* dim、L、seed。 */
    size_t target_molecules; /* 解決済みの生成目標。占有率指定でもこの分子数へ統一。 */
    size_t arm_length;       /* 各腕の歩数。中心を含まない。 */
    NetworkTargetMode target_mode;
    double requested_phi;       /* 占有率指定の入力値。分子数指定なら0（未指定）。 */
    size_t phi_site_budget;     /* 入力占有率×L^3の切り捨て。実績占有点数ではない。 */
    double a_probability;   /* 空中心を選んだ後のA型抽選確率。B型は1-r。 */
} NetworkConfig;
/* dim,L,seed,arm_length,a_probabilityと、target_molecules/target_phiの片方が必須。未知・重複キーは拒否し、失敗時outは変更しない。 */
int network_config_load(const char *path, NetworkConfig *out, char *error, size_t error_size);
/* 数値範囲・既存配列のサイズ上限・必要占有点数を、確保せずに検査する。 */
int network_config_validate(const NetworkConfig *config, char *error, size_t error_size);
#endif
