#ifndef HAMILT_H
#define HAMILT_H

#include "tools.h"
#include "../src_lcao/ORB_control.h"
#include "hamilt_pw.h"

class Hamilt
{
	public:

    Hamilt();
    ~Hamilt();

    void init(const bool use_pw = true);

    void init_k(const int ik);

    void init_before_ions(void);

    void init_before_electrons(void);

    void clear_after_ions(void);

    void cinitcgg(const int ik, const int nstart,const int nbnd,
                  const ComplexMatrix &psi,ComplexMatrix &evc,double *en);

    void diago(const int &istep,const int &iter,const int &ik,
               const double *precondition,double &avg_iter);

    void cdiaghg(const int n,const int m,const ComplexMatrix &hc,const ComplexMatrix &sc,
                 const int ldh,double *e,ComplexMatrix &hvec);

    void h_1psi(const int dim,const complex<double> *psi,complex<double> *hpsi,complex<double> *spi);
    void s_1psi(const int dim,const complex<double> *psi,complex<double> *spsi);

    void h_1psi_real(const int dim,const double *psi,double *hpsi,double *spi);
    void s_1psi_real(const int dim,const double *psi, double *spsi);

    void h_psi( const int dim, const complex<double> *psi,complex<double> *hpsi);

    Hamilt_PW hpw;
	
	// mohan update 2021-02-10
	ORB_control orb_con;

private:

    bool test_exit_cond( const int &ntry, const int &notconv);

};
#endif
