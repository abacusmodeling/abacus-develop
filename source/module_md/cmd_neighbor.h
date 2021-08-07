#ifndef CMD_NEIGHBOR_H
#define CMD_NEIGHBOR_H

#include "../module_base/vector3.h"
#include "../module_cell/unitcell_pseudo.h"

class CMD_neighbor
{
public:

    CMD_neighbor();
    ~CMD_neighbor();

    Vector3<double> Cell_periodic(const Vector3<double> a, const Vector3<double> b);
    void Neighbor(UnitCell_pseudo &ucell_c);

    int **list;    // record the index of adjent atoms of every atom
    int *nlist;    // record the adjent num of every atom
    int dim;
};
#endif