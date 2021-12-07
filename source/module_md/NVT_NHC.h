#ifndef NVT_NHC_H
#define NVT_NHC_H

#include "verlet.h"

class NVT_NHC : public Verlet
{
public:
    NVT_NHC(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in);
    ~NVT_NHC();

    void setup();
    void first_half();
    void second_half();
    void outputMD();
    void write_restart();
    void restart();
    void integrate();

    const static int nc = 1;
    const static int nsy = 7;
    double w[nsy];         // scale evolution operator

    double *Q;             // MHC mass
    double *G;             // parameter in NHC
    double *NHCeta;        // NHC position
    double *NHCveta;       // NHC velocity

};

#endif