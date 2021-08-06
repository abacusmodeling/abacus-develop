#ifndef CMD_NEIGHBOR_H
#define CMD_NEIGHBOR_H

class CMD_neighbor
{
public:

    CMD_neighbor();
    ~CMD_neighbor();

    double Lennard_Jones(UnitCell_pseudo &ucell_c, Grid_Driver &grid_neigh, Vector3<double> *force, matrix &stress);
    double LJ_energy(const double d);
    Vector3<double> LJ_force(const double d, const Vector3<double> dr);

};
#endif