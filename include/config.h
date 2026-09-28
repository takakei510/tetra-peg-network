#ifndef TPN_CONFIG_H
#define TPN_CONFIG_H
#include <stdint.h>
#include <stddef.h>
typedef struct { int dim, L; uint64_t seed; } Config;
/* Returns 0 on success; leaves out untouched on failure. */
int config_load(const char *path, Config *out, char *error, size_t error_size);
#endif
