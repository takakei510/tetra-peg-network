#ifndef TPN_NETWORK_CONFIG_H
#define TPN_NETWORK_CONFIG_H
#include "config.h"
/* 多分子用。既存の3キーConfigを含め、生成条件を別に管理する。 */
typedef struct {
    Config base;             /* dim、L、seed。 */
    size_t target_molecules; /* 完成させたい分子数。実現した分子数とは別。 */
    size_t arm_length;       /* 各腕の歩数。中心を含まない。 */
    double a_probability;   /* 空中心を選んだ後のA型抽選確率。B型は1-r。 */
} NetworkConfig;
/* 6キーはすべて必須。未知・重複キーは拒否し、失敗時outは変更しない。 */
int network_config_load(const char *path, NetworkConfig *out, char *error, size_t error_size);
/* 数値範囲・既存配列のサイズ上限・必要占有点数を、確保せずに検査する。 */
int network_config_validate(const NetworkConfig *config, char *error, size_t error_size);
#endif
