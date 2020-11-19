// AUTHOR : Mohan Chen
// UPDATE LAST TIME : 2010-1-9  
#ifndef OCCUPY_H
#define OCCUPY_H

#include "tools.h"

using namespace std;

class Occupy
{

    friend class Chi0_hilbert;           // pengfei 2016-11-23

public:

    Occupy();
    ~Occupy();

    static void calculate_weights(void);

    static void decision(const string &name,const string &smearing,const double &degauss);

    static const bool& gauss(void) {
        return use_gaussian_broadening;
    }
    static const bool& tetra(void) {
        return use_tetrahedron_method;
    }
    static const bool& fix(void) {
        return fixed_occupations;
    }

public:

    // gaussian_broadening
    static bool use_gaussian_broadening;
    static int gaussian_type;
    static double gaussian_parameter;

    // tetrahedron
    static bool use_tetrahedron_method;

    // fixed occupations
    static bool fixed_occupations;

    static void iweights(const int nks,const double *wk,const int nband,
                         const double &nelec, double **ekb, double &ef, matrix &wg, const int &is, const int *isk);

    static void gweights(const int nks,const double *wk,const int nband,const double &nelec,const double &degauss,
                         const int ngauss,double **ekb, double &ef, double &demet, matrix &wg, const int &is, const int *isk);

    static void tweights(const int nks,const int nspin,const int nband,const double &nelec,const int ntetra,
                         const matrix &tetra,double **ekb,double &ef,matrix &wg);

    static double wsweight(const Vector3<double> &r, Vector3<double> *rws,const int nrws);

private:

    static void efermig(double **ekb,const int nbnd,const int nks,const double &nelec,const double *wk,
                        const double &degauss,const int ngauss,double &ef, const int &is, const int *isk);

    static double sumkg(double **ekb,const int nband,const int nks,const double *wk,const double &degauss,
                        const int ngauss,const double &e, const int &is, const int *isk);

    static double wgauss(const double &x,const int n);

    static double w1gauss(const double &x,const int n);

    //============================
    // Needed in tweights
    //============================
    static void efermit(double **ekb,const int nband,
                        const int nks,const double &nelec,
                        const int nspin,const int ntetra,
                        const matrix &tetra, double &ef);

    static double sumkt(double **ekb,const int nband,const int nks,const int nspin,const int ntetra,
                        const matrix &tetra,const double &eup);

    static void piksort(const int n, double *a);

};

#endif
