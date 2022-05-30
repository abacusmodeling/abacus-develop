#include "surchem.h"

void surchem::cal_totn(const UnitCell &cell, ModulePW::PW_Basis* rho_basis,
                       const complex<double> *Porter_g, complex<double> *N,
                       complex<double> *TOTN) {

    // vloc to N
    complex<double> *vloc_g = new complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(vloc_g, rho_basis->npw);

    // method 1
    GlobalC::UFFT.ToReciSpace(GlobalC::pot.vltot, vloc_g, rho_basis); // now n is vloc in Recispace
    for (int ig = 0; ig < rho_basis->npw; ig++) {
        if(rho_basis->ig_gge0==ig)    continue;
        const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI /
                           (cell.tpiba2 * rho_basis->gg[ig]);
        N[ig] = -vloc_g[ig] / fac;
    }

    if (GlobalV::MY_RANK == 0) {
        N[0] = Porter_g[0];
    }

    for (int ig = 0; ig < rho_basis->npw; ig++) {
        TOTN[ig] = N[ig] - Porter_g[ig];
    }
    delete[] vloc_g;
    return;
}