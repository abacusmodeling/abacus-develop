#include "dipole.h"
#include "../module_base/constants.h"
#include "../module_base/timer.h"
#include "../module_base/global_variable.h"
#include "../src_parallel/parallel_reduce.h"

int Dipole::dir = 2;
double Dipole::dipole_energy = 0.0;
double Dipole::max_pos = 0.5;
double Dipole::de_reg = 0.1;

Dipole::Dipole(){}

Dipole::~Dipole(){}

//=======================================================
// calculate dipole potential in surface calculations
//=======================================================
ModuleBase::matrix Dipole::v_dipole(const UnitCell &cell, 
                                    PW_Basis &pwb, 
                                    const int &nspin, 
                                    const double *const *const rho)
{
    ModuleBase::TITLE("Dipole", "v_dipole");
    ModuleBase::timer::tick("Dipole", "v_dipole");

    double h_inv;    // inverse of height
    if(dir == 0)
    {
        h_inv = sqrt(cell.G.e11 * cell.G.e11 + cell.G.e12 * cell.G.e12 + cell.G.e13 * cell.G.e13);
    }
    else if(dir == 1)
    {
        h_inv = sqrt(cell.G.e21 * cell.G.e21 + cell.G.e22 * cell.G.e22 + cell.G.e23 * cell.G.e23);
    }
    else if(dir = 2)
    {
        h_inv = sqrt(cell.G.e31 * cell.G.e31 + cell.G.e32 * cell.G.e32 + cell.G.e33 * cell.G.e33);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Dipole::ion_dipole", "dipole direction is wrong!");
    }

    double ion_dipole = cal_ion_dipole(cell, h_inv);
    double elec_dipole = cal_elec_dipole(cell, pwb, nspin, rho, h_inv);
    double tot_dipole = ion_dipole - elec_dipole;

    // dipole energy correction
    dipole_energy = 0.5 * ModuleBase::e2 * tot_dipole * tot_dipole * cell.omega / ModuleBase::FOUR_PI;

    // dipole potential
    ModuleBase::matrix v(nspin, pwb.nrxx);
    const int nspin0 = (nspin == 2) ? 2 : 1;

    for (int ir = 0; ir < pwb.nrxx; ++ir)
    {
        int i = ir / (pwb.ncy * pwb.nczp);
        int j = ir / pwb.nczp - i * pwb.ncy;
        int k = ir % pwb.nczp + pwb.nczp_start;
        double x = (double)i / pwb.ncx;
        double y = (double)j / pwb.ncy;
        double z = (double)k / pwb.ncz;
        ModuleBase::Vector3<double> pos(x, y, z);

        double saw = saw_function(max_pos, de_reg, pos[dir]);

        for (int is = 0; is < nspin0; is++)
        {
            v(is, ir) = saw;
        }
    }

    double fac = - ModuleBase::e2 * tot_dipole * cell.lat0 / h_inv;

    ModuleBase::timer::tick("Dipole", "v_dipole");
    return v * fac;
}


//=======================================================
// calculate dipole density in surface calculations
//=======================================================
// double Dipole::dipole(const UnitCell &cell, 
//                         PW_Basis &pwb, 
//                         const double *const *const rho,
//                         complex<double> *TOTN)
// {
//     double *Porter = new double[pwb.nrxx];

//     for (int ir = 0; ir < pwb.nrxx; ++ir)
//     {
//         int i = ir / (pwb.ncy * pwb.nczp);
//         int j = ir / pwb.nczp - i * pwb.ncy;
//         int k = ir % pwb.nczp + pwb.nczp_start;
//         double x = (double)i / pwb.ncx;
//         double y = (double)j / pwb.ncy;
//         double z = (double)k / pwb.ncz;
//         ModuleBase::Vector3<double> pos(x, y, z);

//         pos = cell.latvec * pos;

//         Porter[ir] = cell.lat0 * pos[dir];
//     }

//     complex<double> *Porter_g = new complex<double>[pwb.ngmc];
//     GlobalC::UFFT.ToReciSpace(Porter, Porter_g);

//     double m = 0;
//     for (int ig = 0; ig < pwb.ngmc; ig++)
//     {
//         m += (conj(TOTN[ig]) * Porter_g[ig]).real();
//     }

//     Parallel_Reduce::reduce_double_pool(m);

//     // height
//     const double h = cell.latvec.to_matrix()(dir, dir) * cell.lat0;

//     delete[] Porter, Porter_g;

//     return -m * h * ModuleBase::e2;
// }

double Dipole::cal_ion_dipole(const UnitCell &cell, const double &h_inv)
{
    double ion_dipole = 0;
    for(int it=0; it<cell.ntype; ++it)
    {
        double sum = 0;
        for(int ia=0; ia<cell.atoms[it].na; ++ia)
        {
            sum += saw_function(max_pos, de_reg, cell.atoms[it].taud[ia][dir]);
        }
        ion_dipole += sum * cell.atoms[it].zv;
    }
    ion_dipole *= cell.lat0 / h_inv * ModuleBase::FOUR_PI / cell.omega;

    // std::cout << "ion_dipole = " << ion_dipole << std::endl;

    return ion_dipole;
}

double Dipole::cal_elec_dipole(const UnitCell &cell, 
                            PW_Basis &pwb, 
                            const int &nspin, 
                            const double *const *const rho,
                            const double &h_inv)
{
    double elec_dipole = 0;
    const int nspin0 = (nspin == 2) ? 2 : 1;

    for (int ir = 0; ir < pwb.nrxx; ++ir)
    {
        int i = ir / (pwb.ncy * pwb.nczp);
        int j = ir / pwb.nczp - i * pwb.ncy;
        int k = ir % pwb.nczp + pwb.nczp_start;
        double x = (double)i / pwb.ncx;
        double y = (double)j / pwb.ncy;
        double z = (double)k / pwb.ncz;
        ModuleBase::Vector3<double> pos(x, y, z);

        double saw = saw_function(max_pos, de_reg, pos[dir]);

        for (int is = 0; is < nspin0; is++)
        {
            elec_dipole += rho[is][ir] * saw;
        }
    }

    Parallel_Reduce::reduce_double_pool(elec_dipole);
    elec_dipole *= cell.lat0 / h_inv * ModuleBase::FOUR_PI / pwb.ncxyz;

    // std::cout << "elec_dipole = " << elec_dipole << std::endl;

    return elec_dipole;
}

double Dipole::saw_function(const double &a, const double &b, const double &x)
{
    assert(x>=0);
    assert(x<=1);

    const double fac = 1 - b;

    if( x < a )
    {
        return x - a + 0.5 * fac;
    }
    else if( x > (a+b))
    {
        return x - a - 1 + 0.5 * fac;
    }
    else
    {
        return 0.5 * fac - fac * (x - a) / b;
    }
}