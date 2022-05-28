#include "../module_base/timer.h"
#include "../module_xc/xc_functional.h"
#include "surchem.h"

void shape_gradn(const double *PS_TOTN_real, ModulePW::PW_Basis* rho_basis, double *eprime)
{

    double epr_c = 1.0 / sqrt(ModuleBase::TWO_PI) / GlobalV::sigma_k;
    double epr_z = 0;
    double min = 1e-10;
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        epr_z = log(max(PS_TOTN_real[ir], min) / GlobalV::nc_k) / sqrt(2) / GlobalV::sigma_k;
        eprime[ir] = epr_c * exp(-pow(epr_z, 2)) / max(PS_TOTN_real[ir], min);
    }
}

void eps_pot(const double *PS_TOTN_real, const complex<double> *phi, ModulePW::PW_Basis* rho_basis, double *d_eps, double *vwork)
{
    double *eprime = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(eprime, rho_basis->nrxx);

    shape_gradn(PS_TOTN_real, rho_basis, eprime);

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        eprime[ir] = eprime[ir] * (GlobalV::eb_k - 1);
    }

    ModuleBase::Vector3<double> *nabla_phi = new ModuleBase::Vector3<double>[rho_basis->nrxx];
    double *phisq = new double[rho_basis->nrxx];

    // nabla phi
    XC_Functional::grad_rho(phi, nabla_phi, GlobalC::rhopw);

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        phisq[ir] = pow(nabla_phi[ir].x, 2) + pow(nabla_phi[ir].y, 2) + pow(nabla_phi[ir].z, 2);
    }

    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        vwork[ir] = eprime[ir] * phisq[ir] / (8 * ModuleBase::PI);
    }

    delete[] eprime;
    delete[] nabla_phi;
    delete[] phisq;
}

ModuleBase::matrix surchem::cal_vel(const UnitCell &cell,
                                    ModulePW::PW_Basis* rho_basis,
                                    complex<double> *TOTN,
                                    complex<double> *PS_TOTN,
                                    int nspin)
{
    ModuleBase::TITLE("surchem", "cal_vel");
    ModuleBase::timer::tick("surchem", "cal_vel");

    double *TOTN_real = new double[rho_basis->nrxx];
    GlobalC::UFFT.ToRealSpace(TOTN, TOTN_real,rho_basis);

    // -4pi * TOTN(G)
    complex<double> *B = new complex<double>[rho_basis->npw];
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        B[ig] = -4.0 * ModuleBase::PI * TOTN[ig];
    }

    // Build a nrxx vector to DO FFT .
    double *PS_TOTN_real = new double[rho_basis->nrxx];
    GlobalC::UFFT.ToRealSpace(PS_TOTN, PS_TOTN_real,rho_basis);

    // build epsilon in real space (nrxx)
    double *epsilon = new double[rho_basis->nrxx];
    double *epsilon0 = new double[rho_basis->nrxx];
    cal_epsilon(rho_basis, PS_TOTN_real, epsilon, epsilon0);

    complex<double> *Sol_phi = new complex<double>[rho_basis->npw];
    complex<double> *Sol_phi0 = new complex<double>[rho_basis->npw];
    int ncgsol = 0;

    double *Vel = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(Vel, rho_basis->nrxx);

    // Calculate Sol_phi with epsilon.
    ncgsol = 0;
    minimize_cg(cell, rho_basis, epsilon, B, Sol_phi, ncgsol);

    ncgsol = 0;
    // Calculate Sol_phi0 with epsilon0.
    minimize_cg(cell, rho_basis, epsilon0, B, Sol_phi0, ncgsol);

    double *phi_tilda_R = new double[rho_basis->nrxx];
    double *phi_tilda_R0 = new double[rho_basis->nrxx];
    double *delta_phi_R = new double[rho_basis->nrxx];

    GlobalC::UFFT.ToRealSpace(Sol_phi, phi_tilda_R,rho_basis);
    GlobalC::UFFT.ToRealSpace(Sol_phi0, phi_tilda_R0,rho_basis);

    // the 1st item of Vel
    for (int i = 0; i < rho_basis->nrxx; i++)
    {
        delta_phi_R[i] = phi_tilda_R[i] - phi_tilda_R0[i];
        Vel[i] += delta_phi_R[i];
    }

    // calculate Ael
    double Ael = cal_Ael(cell, rho_basis, TOTN_real, delta_phi_R);

    // the 2nd item of Vel
    double *Vel2 = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(Vel2, rho_basis->nrxx);

    eps_pot(PS_TOTN_real, Sol_phi, rho_basis, epsilon, Vel2);

    for (int i = 0; i < rho_basis->nrxx; i++)
    {
        Vel[i] += Vel2[i];
    }

    ModuleBase::matrix v(nspin, rho_basis->nrxx);

    if (nspin == 4)
    {
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
        {
            v(0, ir) += Vel[ir];
        }
    }
    else
    {
        for (int is = 0; is < nspin; is++)
        {
            for (int ir = 0; ir < rho_basis->nrxx; ir++)
            {
                v(is, ir) += Vel[ir];
            }
        }
    }

    delete[] PS_TOTN_real;
    delete[] Sol_phi;
    delete[] Sol_phi0;
    delete[] B;
    delete[] epsilon;
    delete[] epsilon0;
    delete[] Vel;
    delete[] Vel2;
    delete[] TOTN_real;
    delete[] phi_tilda_R;
    delete[] phi_tilda_R0;
    delete[] delta_phi_R;

    ModuleBase::timer::tick("surchem", "cal_vel");
    return v;
}
