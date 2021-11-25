#include "MSST.h"
#include "MD_func.h"

MSST::MSST(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in) : Verlet(MD_para_in, unit_in)
{
    std::cout << "MSST" << std::endl;

    mdp.Qmass = mdp.Qmass / pow(ModuleBase::ANGSTROM_AU, 4) / pow(ModuleBase::AU_to_MASS, 2);
    mdp.velocity = mdp.velocity * ModuleBase::ANGSTROM_AU * ModuleBase::AU_to_FS;
    mdp.viscosity = mdp.viscosity / ModuleBase::AU_to_MASS / ModuleBase::ANGSTROM_AU * ModuleBase::AU_to_FS;

    old_v = new ModuleBase::Vector3<double> [ucell.nat];
    dilation.set(1,1,1);
    omega.set(0,0,0);
    p0 = 0;
    e0 = 0;
    v0 = 1;
    totmass = 0;

    for(int i=0; i<ucell.nat; ++i)
    {
        totmass += allmass[i];
    }
}

MSST::~MSST()
{
    delete []old_v;
}

void MSST::setup()
{
    ModuleBase::TITLE("MSST", "setup");
    ModuleBase::timer::tick("MSST", "setup");

    int sd = mdp.direction;

    MD_func::force_virial(mdp, ucell, potential, force, virial);
    MD_func::kinetic_stress(ucell, vel, allmass, kinetic, stress);
    stress += virial;

    lag_pos = 0;
    v0 = ucell.omega;
    p0 = stress(sd, sd);
    e0 = potential + kinetic;

    if(kinetic > 0 && mdp.tscale > 0)
    {
        double fac1 = mdp.tscale * totmass * 2.0 * kinetic / mdp.Qmass;
        omega[sd] = -1.0 * sqrt(fac1);
        double fac2 = omega[sd] / v0;

        std::cout << "initial strain rate = " << fac2 << "    tscale = " << mdp.tscale << std::endl;

        for(int i=0; i<ucell.nat; ++i)
        {
            vel[i] *= sqrt(1.0 - mdp.tscale);
        }
    }

    MD_func::kinetic_stress(ucell, vel, allmass, kinetic, stress);
    stress += virial;

    ModuleBase::timer::tick("MSST", "setup");
}

void MSST::first_half()
{
    ModuleBase::TITLE("MSST", "first_half");
    ModuleBase::timer::tick("MSST", "first_half");

    const int sd = mdp.direction;
    double p_current = stress(sd, sd);
    const double dthalf = 0.5 * mdp.dt;

    energy_ = potential + kinetic;

    double p_msst = mdp.velocity * mdp.velocity * totmass * (v0 - ucell.omega) / (v0 * v0);
    double const_A = totmass * (p_current - p0 - p_msst) / mdp.Qmass;
    double const_B = totmass * mdp.viscosity / (mdp.Qmass * ucell.omega);

    // prevent the increase of volume
    if(ucell.omega > v0 && const_A > 0)
    {
        const_A = -const_A;
    }

    // avoid singularity at B = 0 with Taylor expansion
    double fac1 = const_B * dthalf;
    if(fac1 > 1e-6)
    {
        omega[sd] = (omega[sd] + const_A * (exp(fac1) - 1) / const_B) * exp(-fac1);
    }
    else
    {
        omega[sd] += (const_A - const_B * omega[sd]) * dthalf + 
            0.5 * (const_B * const_B * omega[sd] - const_A * const_B) * dthalf * dthalf;
    }

    vsum = vel_sum();

    // propagate velocity sum 1/2 step by temporarily propagating the velocities
    double fac2 = mdp.viscosity * pow(omega[sd], 2) / (vsum * ucell.omega);
    for(int i=0; i<ucell.nat; ++i)
    {
        old_v[i] = vel[i];
        ModuleBase::Vector3<double> const_C = force[i] / allmass[i];
        ModuleBase::Vector3<double> const_D;
        const_D.set(fac2/allmass[i], fac2/allmass[i], fac2/allmass[i]);
        const_D[sd] -= 2 * omega[sd] / ucell.omega;

        for(int k=0; k<3; ++k)
        {
            if( fabs(dthalf*const_D[k]) > 1e-6 )
            {
                double expd = exp(dthalf*const_D[k]);
                vel[i][k] = expd * ( const_C[k] + const_D[k] * vel[i][k] - const_C[k] / expd ) / const_D[k];
            }
            else
            {
                vel[i][k] += ( const_C[k] + const_D[k] * vel[i][k] ) * dthalf + 
                    0.5 * (const_D[k] * const_D[k] * vel[i][k] + const_C[k] * const_D[k] ) * dthalf * dthalf;
            }
        }
    }

    vsum = vel_sum();

    // reset the velocities
    for(int i; i<ucell.nat; ++i)
    {
        vel[i] = old_v[i];
    }

    // propagate velocities 1/2 step using the new velocity sum
    fac2 = mdp.viscosity * pow(omega[sd], 2) / (vsum * ucell.omega);
    for(int i=0; i<ucell.nat; ++i)
    {
        ModuleBase::Vector3<double> const_C = force[i] / allmass[i];
        ModuleBase::Vector3<double> const_D;
        const_D.set(fac2/allmass[i], fac2/allmass[i], fac2/allmass[i]);
        const_D[sd] -= 2 * omega[sd] / ucell.omega;

        for(int k=0; k<3; ++k)
        {
            if( fabs(dthalf*const_D[k]) > 1e-6 )
            {
                double expd = exp(dthalf*const_D[k]);
                vel[i][k] = expd * ( const_C[k] + const_D[k] * vel[i][k] - const_C[k] / expd ) / const_D[k];
            }
            else
            {
                vel[i][k] += ( const_C[k] + const_D[k] * vel[i][k] ) * dthalf + 
                    0.5 * (const_D[k] * const_D[k] * vel[i][k] + const_C[k] * const_D[k] ) * dthalf * dthalf;
            }
        }
    }

    // propagate volume 1/2 step
    double vol = ucell.omega + omega[sd] * dthalf;

    // rescale positions and change box size
    rescale(vol);

    // propagate atom positions 1 time step
    for(int i=0; i<ucell.nat; ++i)
    {
        pos[i] += vel[i] * mdp.dt;
    }
    ucell.update_pos_tau(pos);
    ucell.periodic_boundary_adjustment();

    // propagate volume 1/2 step
    vol = ucell.omega + omega[sd] * dthalf;

    // rescale positions and change box size
    rescale(vol);

    ModuleBase::timer::tick("MSST", "first_half");
}

void MSST::second_half(){}

void MSST::outputMD(){}

double MSST::extra_term()
{
    return 0;
}

double MSST::vel_sum()
{
    double vsum = 0;

    for(int i=0; i<ucell.nat; ++i)
    {
        vsum += vel[i].norm2();
    }

    return vsum;
}

void MSST::rescale(double volume)
{
    int sd = mdp.direction;

    dilation[sd] = volume/ucell.omega;
    ucell.latvec.e11 *= dilation[0];
    ucell.latvec.e22 *= dilation[1];
    ucell.latvec.e33 *= dilation[2];

    ucell.setup_cell_after_vc(GlobalV::ofs_running);
    MD_func::InitPos(ucell, pos);

    // rescale velocity
    for(int i=0; i<ucell.nat; ++i)
    {
        vel[i][sd] *= dilation[sd];
    }
}