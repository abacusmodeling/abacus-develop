//==========================================================
// Author: Xin Qu
// DATE : 2019-12-10
//==========================================================
#ifndef DFTU_RELAX_H
#define DFTU_RELAX_H

#include <vector>
#include "dftu_yukawa.h"
#include "../module_base/complexmatrix.h"

using namespace std;

//==========================================================
// CLASS :
// NAME : DFTU_RELAX
//==========================================================
class DFTU_RELAX : public DFTU_Yukawa
{

public:

    DFTU_RELAX();
    ~DFTU_RELAX();

    void force_stress();
    void cal_force_k(const int ik, const complex<double>* rho_VU);
    void cal_stress_k(const int ik, const complex<double>* rho_VU);
    void cal_force_gamma(const double* rho_VU);
    void cal_stress_gamma(const double* rho_VU);

    void fold_dSR_gamma(const int dim1, const int dim2, double* dSR_gamma);
    void fold_dSm_k(const int ik, const int dim, complex<double>* dSm_k);
    void fold_dSR_k(const int ik, const int dim1, const int dim2, complex<double>* dSR_gamma);

    void cal_VU_pot_mat_complex(const int spin, const bool newlocale, complex<double>* VU);
    void cal_VU_pot_mat_real(const int spin, const bool newlocale, double* VU);

    double get_onebody_eff_pot(
      const int T, const int iat,
	    const int L, const int N, const int spin, 
	    const int m0,
      const int m1,
      const int type, const bool newlocale);

    //forces and stress
    vector<vector<double>> force_dftu;      //force_dftu[iat][dim] 
    vector<vector<double>> stress_dftu;

    //transform between iwt index and it, ia, L, N and m index
    vector<vector<vector<vector<vector<int>>>>> iatlnmipol2iwt;   //iatlnm2iwt[iat][l][n][m][ipol]
    vector<int> iwt2it;                               //iwt2it[iwt]
    vector<int> iwt2l;                                //iwt2l[iwt]
    vector<int> iwt2n;                                //iwt2n[iwt]
    vector<int> iwt2m;                                //iwt2m[iwt]
    vector<int> iwt2ipol;                             //iwt2ipol[iwt]
    vector<int> iat2it; 

    //local occupancy matrix of the correlated subspace
    //locale: the out put local occupation number matrix of correlated electrons in the current electronic step
    //locale_save: the input local occupation number matrix of correlated electrons in the current electronic step
    vector<vector<vector<vector<matrix>>>> locale;            // locale[iat][l][n][spin](m1,m2)
    vector<vector<vector<vector<matrix>>>> locale_save;       // locale_save[iat][l][n][spin](m1,m2)
};

#endif
