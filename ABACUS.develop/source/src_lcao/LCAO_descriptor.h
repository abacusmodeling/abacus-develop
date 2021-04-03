//caoyu add 2021-03-29
#include "../src_global/intarray.h"
#include "../src_global/complexmatrix.h"

#ifndef LCAO_MATRIX_DESCRIPTOR_H
#define LCAO_MATRIX_DESCRIPTOR_H

class LCAO_Descriptor
{
public:
    LCAO_Descriptor();
    ~LCAO_Descriptor();

    void build_S_descriptor(const bool &calc_deri); //cal S_alpha_mu：overlap between lcao basis Phi and descriptor basis Alpha
    void cal_projective_DM();                       //cal PDM: S_alpha_mu * inv(Sloc) * DM * inv(Sloc) * S_nu_beta
    void cal_descriptor();                          //cal d: EIGENVALUE of PDM in block of I_n_l
    void print_descriptor();

private:
    double *S_mu_alpha; //overlap between lcao and descriptor basis
    double *PDM;        //projective density matrix
    double *d;          //descriptors
    int n_descriptor = 0;
    int des_per_atom = 0; //\sum_L{Nchi(L)*(2L+1)}
    IntArray *mu_index;
    void init_mu_index(void);
    void set_S_mu_alpha(const int &iw1_all, const int &iw2_all, const double &v);
    void print_projective_DM(ofstream &ofs, ComplexMatrix &des, const int &it, const int &ia, const int &l, const int &n);
};

#endif