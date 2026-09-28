#ifndef TPN_LATTICE_H
#define TPN_LATTICE_H
#include <stddef.h>
typedef struct { int dim, L; size_t n_sites; int *owner; } Lattice; /* -1 empty; otherwise molecule ID */
int lattice_init(Lattice *lat,int dim,int L);
void lattice_free(Lattice *lat);
int lattice_index(const Lattice *lat,const int coord[3],size_t *out);
int lattice_coord(const Lattice *lat,size_t index,int coord[3]);
int lattice_neighbors(const Lattice *lat,size_t index,size_t out[6]); /* open boundary */
#endif
