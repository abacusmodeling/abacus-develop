// #include "../src_pw/diago_cg.h"
#include "../module_base/constants.h"
#include "../module_base/timer.h"
#include "../module_xc/xc_functional.h"
#include "../src_parallel/parallel_reduce.h"
#include "surchem.h"

#include <cmath>

ModuleBase::matrix surchem::v_correction(const UnitCell &cell,
                                         ModulePW::PW_Basis* rho_basis,
                                         const int &nspin,
                                         const double *const *const rho)
{
    ModuleBase::TITLE("surchem", "v_correction");
    ModuleBase::timer::tick("surchem", "v_correction");

    double *Porter = new double[rho_basis->nrxx];
    for (int i = 0; i < rho_basis->nrxx; i++)
        Porter[i] = 0.0;
    const int nspin0 = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; is++)
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
            Porter[ir] += rho[is][ir];

    complex<double> *Porter_g = new complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(Porter_g, rho_basis->npw);

    GlobalC::UFFT.ToReciSpace(Porter, Porter_g, rho_basis);

    complex<double> *N = new complex<double>[rho_basis->npw];
    complex<double> *TOTN = new complex<double>[rho_basis->npw];
    complex<double> *PS_TOTN = new complex<double>[rho_basis->npw];

    cal_totn(cell, rho_basis, Porter_g, N, TOTN);

    cal_pseudo(cell, rho_basis, Porter_g, PS_TOTN);

    ModuleBase::matrix v(nspin, rho_basis->nrxx);

    v += cal_vel(cell, rho_basis, TOTN, PS_TOTN, nspin);
    v += cal_vcav(cell, rho_basis, PS_TOTN, nspin);

    delete[] Porter;
    delete[] Porter_g;
    delete[] N;
    delete[] PS_TOTN;
    delete[] TOTN;

    ModuleBase::timer::tick("surchem", "v_correction");
    return v;
}

void surchem::add_comp_chg(const UnitCell &cell, ModulePW::PW_Basis* rho_basis, double q, double l, double center, complex<double> *NG, int dim)
{
    // x dim
    double tmp_q = 0.0;
    if (dim == 0)
    {
        double L = cell.a1[0];
        tmp_q = q / (cross(cell.a2, cell.a3).norm() * l);
        ModuleBase::GlobalFunc::ZEROS(NG, rho_basis->npw);
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            if(ig==rho_basis->ig_gge0)
                continue;
            double GX = rho_basis->gcar[ig][0];
            double GY = rho_basis->gcar[ig][1];
            double GZ = rho_basis->gcar[ig][2];
            GX = GX * 2 * ModuleBase::PI;
            if (GY == 0 && GZ == 0 && GX != 0)
            {
                NG[ig] = exp(ModuleBase::NEG_IMAG_UNIT * GX * center) * complex<double>(2.0 * tmp_q * sin(GX * l / 2.0) / (L * GX), 0.0);
            }
        }
        // NG[0] = complex<double>(tmp_q * l / L, 0.0);
    }
    // y dim
    else if (dim == 1)
    {
        double L = cell.a2[1];
        tmp_q = q / (cross(cell.a1, cell.a3).norm() * l);
        ModuleBase::GlobalFunc::ZEROS(NG, rho_basis->npw);
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            if(ig==rho_basis->ig_gge0)
                continue;
            double GX = rho_basis->gcar[ig][0];
            double GY = rho_basis->gcar[ig][1];
            double GZ = rho_basis->gcar[ig][2];
            GY = GY * 2 * ModuleBase::PI;
            if (GX == 0 && GZ == 0 && GY != 0)
            {
                NG[ig] = exp(ModuleBase::NEG_IMAG_UNIT * GY * center) * complex<double>(2.0 * tmp_q * sin(GY * l / 2.0) / (L * GY), 0.0);
            }
        }
        // NG[0] = complex<double>(tmp_q * l / L, 0.0);
    }
    // z dim
    else if (dim == 2)
    {
        double L = cell.a3[2];
        // cout << "area" << cross(cell.a1, cell.a2).norm() << endl;
        tmp_q = q / (cross(cell.a1, cell.a2).norm() * l);
        ModuleBase::GlobalFunc::ZEROS(NG, rho_basis->npw);
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            if(ig==rho_basis->ig_gge0)
                continue;
            double GX = rho_basis->gcar[ig][0];
            double GY = rho_basis->gcar[ig][1];
            double GZ = rho_basis->gcar[ig][2];
            GZ = GZ * 2 * ModuleBase::PI;
            if (GX == 0 && GY == 0 && GZ != 0)
            {
                NG[ig] = exp(ModuleBase::NEG_IMAG_UNIT * GZ * center) * complex<double>(2.0 * tmp_q * sin(GZ * l / 2.0) / (L * GZ), 0.0);
            }
        }
        // NG[0] = complex<double>(tmp_q * l / L, 0.0);
    }
}

ModuleBase::matrix surchem::v_compensating(const UnitCell &cell, ModulePW::PW_Basis *rho_basis)
{
    ModuleBase::TITLE("surchem", "v_compensating");
    ModuleBase::timer::tick("surchem", "v_compensating");

    complex<double> *comp_reci = new complex<double>[rho_basis->npw];
    complex<double> *phi_comp_G = new complex<double>[rho_basis->npw];
    double *phi_comp_R = new double[rho_basis->nrxx];

    ModuleBase::GlobalFunc::ZEROS(comp_reci, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(phi_comp_G, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(phi_comp_R, rho_basis->nrxx);
    // get comp chg in reci space
    add_comp_chg(cell, rho_basis, comp_q, comp_l, comp_center, comp_reci, comp_dim);
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if (rho_basis->gg[ig] >= 1.0e-12) // LiuXh 20180410
        {
            const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (cell.tpiba2 * rho_basis->gg[ig]);
            phi_comp_G[ig] = fac * comp_reci[ig];
        }
    }

    GlobalC::UFFT.ToRealSpace(phi_comp_G, phi_comp_R, rho_basis);

    ModuleBase::matrix v_comp(GlobalV::NSPIN, rho_basis->nrxx);
    if (GlobalV::NSPIN == 4)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
            v_comp(0, ir) = phi_comp_R[ir];
    }
    else
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
            for (int ir = 0; ir < rho_basis->nrxx; ir++)
                v_comp(is, ir) = phi_comp_R[ir];
    }

    delete[] comp_reci;
    delete[] phi_comp_G;
    delete[] phi_comp_R;

    ModuleBase::timer::tick("surchem", "v_compensating");
    return v_comp;
}

void test_print(double* tmp, ModulePW::PW_Basis *rho_basis)
{
    for (int i = 0; i < rho_basis->nz; i++)
    {
        cout << tmp[i] << endl;
    }
}

void surchem::test_V_to_N(ModuleBase::matrix &v, 
                const UnitCell &cell, 
                ModulePW::PW_Basis *rho_basis, 
                const double *const *const rho)
{
    double *phi_comp_R = new double[rho_basis->nrxx];
    complex<double> *phi_comp_G = new complex<double>[rho_basis->npw];
    complex<double> *comp_reci = new complex<double>[rho_basis->npw];
    double *N_real = new double[rho_basis->nrxx];

    ModuleBase::GlobalFunc::ZEROS(phi_comp_R, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(phi_comp_G, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(comp_reci, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(N_real, rho_basis->nrxx);

    for (int ir = 0; ir < rho_basis->nz; ir++)
    {
        cout << v(0, ir) << endl;
    }

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        phi_comp_R[ir] = v(0, ir);
    }

    GlobalC::UFFT.ToReciSpace(phi_comp_R, phi_comp_G, rho_basis);
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if (rho_basis->gg[ig] >= 1.0e-12) // LiuXh 20180410
        {
            const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (cell.tpiba2 * rho_basis->gg[ig]);
            comp_reci[ig] = phi_comp_G[ig] / fac;
        }
    }
    GlobalC::UFFT.ToRealSpace(comp_reci, N_real, rho_basis);

    complex<double> *vloc_g = new complex<double>[rho_basis->npw];
    complex<double> *ng = new complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(vloc_g, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(ng, rho_basis->npw);

    double* Porter = new double[rho_basis->nrxx];
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
        Porter[ir] = rho[0][ir];

    GlobalC::UFFT.ToReciSpace(GlobalC::pot.vltot,
                                  vloc_g, rho_basis); // now n is vloc in Recispace
    for (int ig = 0; ig < rho_basis->npw; ig++) {
        if (rho_basis->gg[ig] >= 1.0e-12) // LiuXh 20180410
        {
            const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI /
                               (cell.tpiba2 * rho_basis->gg[ig]);

            ng[ig] = -vloc_g[ig] / fac;
        }
    }
    double *nr = new double[rho_basis->nrxx];
    GlobalC::UFFT.ToRealSpace(ng, nr, rho_basis);

    double *diff = new double[rho_basis->nrxx];
    double *diff2 = new double[rho_basis->nrxx];
    for (int i = 0; i < rho_basis->nrxx; i++)
    {
        diff[i] = N_real[i] - nr[i];
        diff2[i] = N_real[i] - Porter[i];
    }

    for (int i = 0; i < rho_basis->nrxx;i++)
    {
        diff[i] -= Porter[i];
    }

    delete[] phi_comp_R;
    delete[] phi_comp_G;
    delete[] comp_reci;
    delete[] diff;
    delete[] vloc_g;
    delete[] Porter;
}
