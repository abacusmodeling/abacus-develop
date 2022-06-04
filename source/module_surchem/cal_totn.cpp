#include "surchem.h"

void add_comp_chg(const UnitCell &cell, ModulePW::PW_Basis* rho_basis, double q, double l, double center, complex<double> *NG, int dim)
{
    // x dim
    if (dim == 0)
    {
        double L = cell.a1[0];
        q /= (cell.a2[1] * cell.a3[2] * l);
        ModuleBase::GlobalFunc::ZEROS(NG, rho_basis->npw);
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            if(ig == rho_basis->ig_gge0) continue;
            double GX = rho_basis->gcar[ig][0];
            double GY = rho_basis->gcar[ig][1];
            double GZ = rho_basis->gcar[ig][2];
            GX = GX * 2 * ModuleBase::PI;
            if (GY == 0 && GZ == 0 && GX != 0)
            {
                NG[ig] = exp(ModuleBase::NEG_IMAG_UNIT * GX * center) * complex<double>(2.0 * q * sin(GX * l / 2.0) / (L * GX), 0.0);
            }
        }
        NG[0] = complex<double>(q * l / L, 0.0);
    }
    // y dim
    else if (dim == 1)
    {
        double L = cell.a2[1];
        q /= (cell.a1[0] * cell.a3[2] * l);
        ModuleBase::GlobalFunc::ZEROS(NG, rho_basis->npw);
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            if(ig == rho_basis->ig_gge0) continue;
            double GX = rho_basis->gcar[ig][0];
            double GY = rho_basis->gcar[ig][1];
            double GZ = rho_basis->gcar[ig][2];
            GY = GY * 2 * ModuleBase::PI;
            if (GX == 0 && GZ == 0 && GY != 0)
            {
                NG[ig] = exp(ModuleBase::NEG_IMAG_UNIT * GY * center) * complex<double>(2.0 * q * sin(GY * l / 2.0) / (L * GY), 0.0);
            }
        }
        NG[0] = complex<double>(q * l / L, 0.0);
    }
    // z dim
    else if (dim == 2)
    {
        double L = cell.a3[2];
        q /= (cell.a1[0] * cell.a2[1] * l);
        ModuleBase::GlobalFunc::ZEROS(NG, rho_basis->npw);
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            if(ig == rho_basis->ig_gge0) continue;
            double GX = rho_basis->gcar[ig][0];
            double GY = rho_basis->gcar[ig][1];
            double GZ = rho_basis->gcar[ig][2];
            GZ = GZ * 2 * ModuleBase::PI;
            if (GX == 0 && GY == 0 && GZ != 0)
            {
                NG[ig] = exp(ModuleBase::NEG_IMAG_UNIT * GZ * center) * complex<double>(2.0 * q * sin(GZ * l / 2.0) / (L * GZ), 0.0);
            }
        }
        NG[0] = complex<double>(q * l / L, 0.0);
    }
}

void test_sum(const UnitCell &cell, ModulePW::PW_Basis* rho_basis, double* NR)
{
    double sum = 0.0;
    for (int i=0; i < rho_basis->nrxx; i++)
    {
        sum += NR[i];
    }
    cout << sum * cell.omega / rho_basis->nrxx << endl;
}

void surchem::cal_totn(const UnitCell &cell, ModulePW::PW_Basis* rho_basis,
                       const complex<double> *Porter_g, complex<double> *N,
                       complex<double> *TOTN) {
    /*
    // test compensating charge
    complex<double> *comp_reci = new complex<double>[rho_basis->npw];
    double *comp_real = new double[pwb.nrxx];
    add_comp_chg(cell, pwb, comp_q, comp_l, comp_center, comp_reci, comp_dim);
    GlobalC::UFFT.ToRealSpace(comp_reci, comp_real);
    for (int ir = 0; ir < pwb.nz; ir++)
    {
        cout << comp_real[ir] << endl;
    }
    cout << "sum" << endl;
    test_sum(cell, pwb, comp_real);
    delete[] comp_real;
    delete[] comp_reci;
    */

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
        // cout << "Ng[0]" << N[0] << endl;
    }

    for (int ig = 0; ig < rho_basis->npw; ig++) {
        TOTN[ig] = N[ig] - Porter_g[ig];
    }

    delete[] vloc_g;
    return;
}