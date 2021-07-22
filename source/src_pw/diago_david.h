//==========================================================
// AUTHOR : wangjp
// Data :2009-04
// Last Update:
//
// 09-05-10 modify SchmitOrth() diag_zhegvx() as static
// member function
//==========================================================

#ifndef Diago_David_H
#define Diago_David_H

#include "tools.h"

class Diago_David
{
public:

    Diago_David();
    ~Diago_David();

    static void SchmitOrth(
        const int& npw,
        const int n_band,
        const int m,
        const ComplexMatrix &psi,
        complex<double>* psi_m,
        complex<double>* spsi);

    static void diag_zhegvx(
        const int& n,
        const int& m,
        const ComplexMatrix &hc,
        const ComplexMatrix &sc,
        const int& ldh,
        double* e,
        ComplexMatrix &vc);

    void diag(
        ComplexMatrix &psi,
        double *en,
        const int &npw,
        const int &nband,
        const double *precondition,
        const int order,
        const double &eps,
        const int &maxiter,
        int &notconv,
        double &avg_iter);

private:

    int test_david;

    void cal_grad(
        const int& npw,
        const int& nbase,
        const int& notconv,
        ComplexMatrix &basis,
        ComplexMatrix &hp,
        ComplexMatrix &sp,
        const ComplexMatrix &vc,
        const int* unconv,
        const double* precondition,
        const double* en,
        complex<double>* hpsi,
        complex<double>* spsi,
        complex<double>* ppsi,
        complex<double>* respsi);

    void cal_elem(
        const int& npw,
        int& nbase,
        const int& notconv,
        const ComplexMatrix &basis,
        const ComplexMatrix &hp,
        const ComplexMatrix &sp,
        ComplexMatrix &hc,
        ComplexMatrix &sc);


    void refresh(
        const int& npw,
        const int& nband,
        int& nbase,
        const double* en,
        const ComplexMatrix &psi,
        ComplexMatrix &basis,
        ComplexMatrix &hp,
        ComplexMatrix &sp,
        ComplexMatrix &hc,
        ComplexMatrix &sc,
        ComplexMatrix &vc);

    void cal_err(
        const int& npw,
        const int& nband,
        const int& nbase,
        const ComplexMatrix &vc,
        const ComplexMatrix &hp,
        const ComplexMatrix &basis,
        const double* en,
        complex<double>* respsi);

};

#endif
