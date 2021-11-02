#ifndef PWBASISK_H
#define PWBASISK_H

#include "pw_basis.h"
#include "../module_base/intarray.h"

//
//Special pw_basis class.
//It includes different k-points
//plane waves: <r|g,k>=1/sqrt(V)*exp(i(k+g)r)
//
class PW_Basis_K : public PW_Basis
{

public:
    PW_Basis_K();
    ~PW_Basis_K();
    void initparameters(
        bool gamma_only_in,
        double ecut_in,
        double gk_ecut_in,
        int nk_in, //number of k points in this pool
        ModuleBase::Vector3<double> *kvec_d, // Direct coordinates of k points
        int poolnproc_in, // Number of processors in this pool
        int poolrank_in, // Rank in this pool
        int distribution_type_in
    );
    void setupIndGk(); //set up igk



public:
    int nks;//number of k points in this pool
    ModuleBase::Vector3<double> *kvec_d; // Direct coordinates of k points
    ModuleBase::Vector3<double> *kvec_c; // Cartesian coordinates of k points
    ModuleBase::IntArray igk; //[nks, npw_max] map igk_local to ig_local
    int *npwk; //[nks] number of plane waves of different k-points
    int npwk_max; //max npwk among all nks k-points
    double gk_ecut; //Energy cut off for (g+k)^2/2

public:
    void init_k();//initialize some data for current k-points
    //After inik_k()
    int *GR_index; //[npw_max] map igk_local to (is,iz) of current k, used after inik_k()

public:
    //operator
    ModuleBase::Vector3<double> get_GPlusK_cartesian(const int ik, const int ig) const;
    
};
#endif //PlaneWave_K class
