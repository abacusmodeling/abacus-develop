#ifndef WF_IGK_H
#define WF_IGK_H

#include "tools.h"
#include "pw_basis.h"

class WF_igk 
{
	public:

    WF_igk();
    ~WF_igk();

    // for each k point , find the number of pws
    int setupIndGk(const PW_Basis &pwb,const int nks);
    
    int npwx;
    int npw;
    IntArray igk;

    // g2kin : [npw],kinetic energy for current k point
    double *g2kin;

    // Calculate kinetic energy
    void ekin(const int ik);

    double* get_qvec_cartesian(const int &ik);

    Vector3<double> get_1qvec_cartesian(const int ik,const int ig)const;

    complex<double>* get_sk(const int ik, const int it, const int ia)const;

	// pengfei 2016-11-23
    complex<double>* get_skq(int ik, int it, int ia, Vector3<double> q);

};
#endif
