#include "dipole.h"
#include "../module_base/constants.h"
#include "../module_base/timer.h"
#include "../module_base/global_variable.h"
#include "../src_parallel/parallel_reduce.h"

int Dipole::dir = 2;
double Dipole::dipole_energy = 0.0;

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

    const double m = dipole_density(cell, pwb, nspin, rho);
    const double fac = m * ModuleBase::FOUR_PI;
    ModuleBase::Vector3<int> gird;
    gird.set(pwb.ncx, pwb.ncy, pwb.ncz);

    ModuleBase::matrix v(nspin, pwb.nrxx);
    double e_elec = 0;

    for (int ir = 0; ir < pwb.nrxx; ++ir)
    {
        int i = ir / (pwb.ncy * pwb.nczp);
        int j = ir / pwb.nczp - i * pwb.ncy;
        int k = ir % pwb.nczp + pwb.nczp_start;
        ModuleBase::Vector3<int> index;
        index.set(i,j,k);

        if (nspin == 4)
        {
            v(0, ir) = index[dir] / gird[dir] - 0.5;
            e_elec += v(0, ir) * rho[0][ir];
        }
        else
        {
            for (int is = 0; is < nspin; is++)
            {
                v(is, ir) = index[dir] / gird[dir] - 0.5;
                e_elec += v(is, ir) * rho[is][ir];
            }
        }
    }

    Parallel_Reduce::reduce_double_pool(e_elec);
    e_elec *= cell.omega / pwb.ncxyz;

    // which_ip: found iz belongs to which ip.
    int *which_ip = new int[pwb.ncz];
    ModuleBase::GlobalFunc::ZEROS(which_ip, pwb.ncz);

#ifdef __MPI
    // num_z: how many planes on processor 'ip'
    int *num_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
    for (int iz=0;iz<pwb.nbz;iz++)
    {
        int ip = iz % GlobalV::NPROC_IN_POOL;
        num_z[ip] += pwb.bz;
    }

    // start_z: start position of z in processor ip.
    int *start_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(start_z, GlobalV::NPROC_IN_POOL);
    for (int ip=1;ip<GlobalV::NPROC_IN_POOL;ip++)
    {
        start_z[ip] = start_z[ip-1]+num_z[ip-1];
    }

    for(int iz=0; iz<pwb.ncz; iz++)
    {
        for(int ip=0; ip<GlobalV::NPROC_IN_POOL; ip++)
        {
            if(iz>=start_z[GlobalV::NPROC_IN_POOL-1]) 
            {
                which_ip[iz] = GlobalV::NPROC_IN_POOL-1;
                break;
            }
            else if(iz>=start_z[ip] && iz<start_z[ip+1])
            {
                which_ip[iz] = ip;
                break;
            }
        }
    }

    delete[] num_z;
    delete[] start_z;
#endif

    double e_ion = 0;
    for(int it=0; it<cell.ntype; ++it)
    {
        for(int ia=0; ia<cell.atoms[it].na; ++ia)
        {
            int iz = cell.atoms[it].taud[ia][2] * pwb.ncz;    // FFT divided by z axis
            if(GlobalV::RANK_IN_POOL == which_ip[iz])
            {
                int ix = cell.atoms[it].taud[ia][0] * pwb.ncx;
                int iy = cell.atoms[it].taud[ia][1] * pwb.ncy;
                iz -= pwb.nczp_start;
                int index = ix * pwb.ncy * pwb.nczp + iy * pwb.nczp + iz;

                e_ion += cell.atoms[it].zv * v(0, index);
            }
        }
    }

    Parallel_Reduce::reduce_double_pool(e_ion);

    dipole_energy = 0.5 * fac * (e_ion - e_elec);
    v *= fac;

    delete[] which_ip;

    ModuleBase::timer::tick("Dipole", "v_dipole");
    return v;
}


//=======================================================
// calculate dipole density in surface calculations
//=======================================================
double Dipole::dipole_density(const UnitCell &cell, 
                                PW_Basis &pwb, 
                                const int &nspin, 
                                const double *const *const rho)
{
    // ion part
    double m_ion = 0;
    for(int it=0; it<cell.ntype; ++it)
    {
        double sum = 0;
        for(int ia=0; ia<cell.atoms[it].na; ++ia)
        {
            sum += cell.atoms[it].tau[ia][dir];
        }
        m_ion += sum * cell.atoms[it].zv;
    }
    m_ion *= cell.lat0 * ModuleBase::e2;   //  multiply ModuleBase::e2 to convert Ha to Ry ??

    // height and surface area
    const double h = cell.latvec.to_matrix()(dir, dir) * cell.lat0;
    const double area = cell.omega / h;

    m_ion /= area;


    // electron part
    double m_elec = 0;
    const int nspin0 = (nspin == 2) ? 2 : 1;

    
    for (int ir = 0; ir < pwb.nrxx; ++ir)
    {
        int i = ir / (pwb.ncy * pwb.nczp);
        int j = ir / pwb.nczp - i * pwb.ncy;
        int k = ir % pwb.nczp + pwb.nczp_start;
        ModuleBase::Vector3<int> index;
        index.set(i,j,k);

        for (int is = 0; is < nspin0; is++)
        {
            m_elec += rho[is][ir] * index[dir];
        }
    }

    Parallel_Reduce::reduce_double_pool(m_elec);

    ModuleBase::Vector3<int> gird;
    gird.set(pwb.ncx, pwb.ncy, pwb.ncz);
    
    m_elec *= h * h / pwb.ncxyz / gird[dir] * ModuleBase::e2;

    return m_ion - m_elec;
}

