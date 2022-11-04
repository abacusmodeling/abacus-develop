#include "lattice_change_basic.h"

#include "../src_parallel/parallel_common.h"
#include "../src_pw/global.h"

int Lattice_Change_Basic::dim = 0;
bool Lattice_Change_Basic::converged = true;
double Lattice_Change_Basic::largest_grad = 0.0;
int Lattice_Change_Basic::update_iter = 0;
int Lattice_Change_Basic::istep = 0;
int Lattice_Change_Basic::stress_step = 0;

double Lattice_Change_Basic::ediff = 0.0;
double Lattice_Change_Basic::etot = 0.0;
double Lattice_Change_Basic::etot_p = 0.0;

// double Lattice_Change_Basic::lattice_change_ini = 0.5; // default is 0.5
double Lattice_Change_Basic::lattice_change_ini = 0.01; // default is 0.5

int Lattice_Change_Basic::out_stru = 0;

void Lattice_Change_Basic::setup_gradient(double *lat, double *grad, ModuleBase::matrix &stress)
{
    ModuleBase::TITLE("Lattice_Change_Basic", "setup_gradient");

    if (INPUT.fixed_axes == "volume")
    {
        double stress_aver = (stress(0, 0) + stress(1, 1) + stress(2, 2)) / 3.0;
        stress(0, 0) = stress(0, 0) - stress_aver;
        stress(1, 1) = stress(1, 1) - stress_aver;
        stress(2, 2) = stress(2, 2) - stress_aver;
    }

    lat[0] = GlobalC::ucell.latvec.e11 * GlobalC::ucell.lat0;
    lat[1] = GlobalC::ucell.latvec.e12 * GlobalC::ucell.lat0;
    lat[2] = GlobalC::ucell.latvec.e13 * GlobalC::ucell.lat0;
    lat[3] = GlobalC::ucell.latvec.e21 * GlobalC::ucell.lat0;
    lat[4] = GlobalC::ucell.latvec.e22 * GlobalC::ucell.lat0;
    lat[5] = GlobalC::ucell.latvec.e23 * GlobalC::ucell.lat0;
    lat[6] = GlobalC::ucell.latvec.e31 * GlobalC::ucell.lat0;
    lat[7] = GlobalC::ucell.latvec.e32 * GlobalC::ucell.lat0;
    lat[8] = GlobalC::ucell.latvec.e33 * GlobalC::ucell.lat0;

    if (GlobalC::ucell.lc[0] == 1)
    {
        grad[0] = -(lat[0] * stress(0, 0) + lat[1] * stress(1, 0) + lat[2] * stress(2, 0));
        grad[1] = -(lat[0] * stress(0, 1) + lat[1] * stress(1, 1) + lat[2] * stress(2, 1));
        grad[2] = -(lat[0] * stress(0, 2) + lat[1] * stress(1, 2) + lat[2] * stress(2, 2));
    }
    if (GlobalC::ucell.lc[1] == 1)
    {
        grad[3] = -(lat[3] * stress(0, 0) + lat[4] * stress(1, 0) + lat[5] * stress(2, 0));
        grad[4] = -(lat[3] * stress(0, 1) + lat[4] * stress(1, 1) + lat[5] * stress(2, 1));
        grad[5] = -(lat[3] * stress(0, 2) + lat[4] * stress(1, 2) + lat[5] * stress(2, 2));
    }
    if (GlobalC::ucell.lc[2] == 1)
    {
        grad[6] = -(lat[6] * stress(0, 0) + lat[7] * stress(1, 0) + lat[8] * stress(2, 0));
        grad[7] = -(lat[6] * stress(0, 1) + lat[7] * stress(1, 1) + lat[8] * stress(2, 1));
        grad[8] = -(lat[6] * stress(0, 2) + lat[7] * stress(1, 2) + lat[8] * stress(2, 2));
    }

    // grad[0] = -stress(0,0);   grad[1] = -stress(0,1);  grad[2] = -stress(0,2);
    // grad[3] = -stress(1,0);   grad[4] = -stress(1,1);  grad[5] = -stress(1,2);
    // grad[6] = -stress(2,0);   grad[7] = -stress(2,1);  grad[8] = -stress(2,2);

    return;
}

void Lattice_Change_Basic::change_lattice(double *move, double *lat)
{
    ModuleBase::TITLE("Lattice_Change_Basic", "change_lattice");

    assert(move != NULL);
    assert(lat != NULL);

    /*
        std::cout<<" LATTICE CONSTANT  OLD:"<<std::endl;
        std::cout<<" "<<std::setprecision(12)<<GlobalC::ucell.latvec.e11<<"   "<<GlobalC::ucell.latvec.e12<<"
       "<<GlobalC::ucell.latvec.e13<<std::endl; std::cout<<" "<<std::setprecision(12)<<GlobalC::ucell.latvec.e21<<"
       "<<GlobalC::ucell.latvec.e22<<"   "<<GlobalC::ucell.latvec.e23<<std::endl; std::cout<<"
       "<<std::setprecision(12)<<GlobalC::ucell.latvec.e31<<"   "<<GlobalC::ucell.latvec.e32<<"
       "<<GlobalC::ucell.latvec.e33<<std::endl;
    */

    if (GlobalC::ucell.lc[0] != 0)
    {
        GlobalC::ucell.latvec.e11 = (move[0] + lat[0]) / GlobalC::ucell.lat0;
        GlobalC::ucell.latvec.e12 = (move[1] + lat[1]) / GlobalC::ucell.lat0;
        GlobalC::ucell.latvec.e13 = (move[2] + lat[2]) / GlobalC::ucell.lat0;
    }
    if (GlobalC::ucell.lc[1] != 0)
    {
        GlobalC::ucell.latvec.e21 = (move[3] + lat[3]) / GlobalC::ucell.lat0;
        GlobalC::ucell.latvec.e22 = (move[4] + lat[4]) / GlobalC::ucell.lat0;
        GlobalC::ucell.latvec.e23 = (move[5] + lat[5]) / GlobalC::ucell.lat0;
    }
    if (GlobalC::ucell.lc[2] != 0)
    {
        GlobalC::ucell.latvec.e31 = (move[6] + lat[6]) / GlobalC::ucell.lat0;
        GlobalC::ucell.latvec.e32 = (move[7] + lat[7]) / GlobalC::ucell.lat0;
        GlobalC::ucell.latvec.e33 = (move[8] + lat[8]) / GlobalC::ucell.lat0;
    }

    GlobalC::ucell.a1.x = GlobalC::ucell.latvec.e11;
    GlobalC::ucell.a1.y = GlobalC::ucell.latvec.e12;
    GlobalC::ucell.a1.z = GlobalC::ucell.latvec.e13;
    GlobalC::ucell.a2.x = GlobalC::ucell.latvec.e21;
    GlobalC::ucell.a2.y = GlobalC::ucell.latvec.e22;
    GlobalC::ucell.a2.z = GlobalC::ucell.latvec.e23;
    GlobalC::ucell.a3.x = GlobalC::ucell.latvec.e31;
    GlobalC::ucell.a3.y = GlobalC::ucell.latvec.e32;
    GlobalC::ucell.a3.z = GlobalC::ucell.latvec.e33;

    GlobalC::ucell.omega
        = abs(GlobalC::ucell.latvec.Det()) * GlobalC::ucell.lat0 * GlobalC::ucell.lat0 * GlobalC::ucell.lat0;

    GlobalC::ucell.GT = GlobalC::ucell.latvec.Inverse();
    GlobalC::ucell.G = GlobalC::ucell.GT.Transpose();
    GlobalC::ucell.GGT = GlobalC::ucell.G * GlobalC::ucell.GT;
    GlobalC::ucell.invGGT = GlobalC::ucell.GGT.Inverse();

#ifdef __MPI
    // distribute lattice vectors.
    Parallel_Common::bcast_double(GlobalC::ucell.latvec.e11);
    Parallel_Common::bcast_double(GlobalC::ucell.latvec.e12);
    Parallel_Common::bcast_double(GlobalC::ucell.latvec.e13);
    Parallel_Common::bcast_double(GlobalC::ucell.latvec.e21);
    Parallel_Common::bcast_double(GlobalC::ucell.latvec.e22);
    Parallel_Common::bcast_double(GlobalC::ucell.latvec.e23);
    Parallel_Common::bcast_double(GlobalC::ucell.latvec.e31);
    Parallel_Common::bcast_double(GlobalC::ucell.latvec.e32);
    Parallel_Common::bcast_double(GlobalC::ucell.latvec.e33);

    // distribute lattice vectors.
    Parallel_Common::bcast_double(GlobalC::ucell.a1.x);
    Parallel_Common::bcast_double(GlobalC::ucell.a1.y);
    Parallel_Common::bcast_double(GlobalC::ucell.a1.z);
    Parallel_Common::bcast_double(GlobalC::ucell.a2.x);
    Parallel_Common::bcast_double(GlobalC::ucell.a2.y);
    Parallel_Common::bcast_double(GlobalC::ucell.a2.z);
    Parallel_Common::bcast_double(GlobalC::ucell.a3.x);
    Parallel_Common::bcast_double(GlobalC::ucell.a3.y);
    Parallel_Common::bcast_double(GlobalC::ucell.a3.z);
#endif
    /*
            std::cout<<" LATTICE CONSTANT NEW: "<<std::endl;
            std::cout<<" "<<std::setprecision(12)<<GlobalC::ucell.latvec.e11<<"   "<<GlobalC::ucell.latvec.e12<<"
       "<<GlobalC::ucell.latvec.e13<<std::endl; std::cout<<" "<<std::setprecision(12)<<GlobalC::ucell.latvec.e21<<"
       "<<GlobalC::ucell.latvec.e22<<"   "<<GlobalC::ucell.latvec.e23<<std::endl; std::cout<<"
       "<<std::setprecision(12)<<GlobalC::ucell.latvec.e31<<"   "<<GlobalC::ucell.latvec.e32<<"
       "<<GlobalC::ucell.latvec.e33<<std::endl;
    */
        //--------------------------------------------
        // Print out the structure file.
        //--------------------------------------------
    std::stringstream ss;
    ss << GlobalV::global_out_dir << "STRU_ION";
    if (out_stru == 1)
    {
        ss << istep;
        GlobalC::ucell.print_cell_cif("STRU_NOW.cif");
    }
    ss << "_D";    
    GlobalC::ucell.print_stru_file(ss.str(), 2, 0);

    return;
}

void Lattice_Change_Basic::check_converged(ModuleBase::matrix &stress, double *grad)
{
    ModuleBase::TITLE("Lattice_Change_Basic", "check_converged");

    Lattice_Change_Basic::largest_grad = 0.0;
    double stress_ii_max = 0.0;

    if (GlobalC::ucell.lc[0] == 1 && GlobalC::ucell.lc[1] == 1 && GlobalC::ucell.lc[2] == 1)
    {
        for (int i = 0; i < 3; i++)
        {
            if (stress_ii_max < abs(stress(i, i)))
                stress_ii_max = abs(stress(i, i));
            for (int j = 0; j < 3; j++)
            {
                if (Lattice_Change_Basic::largest_grad < abs(stress(i, j)))
                {
                    Lattice_Change_Basic::largest_grad = abs(stress(i, j));
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < 9; i++)
        {
            if (Lattice_Change_Basic::largest_grad < abs(grad[i]))
            {
                Lattice_Change_Basic::largest_grad = abs(grad[i]);
            }
        }
    }

    double unit_transform = 0.0;
    unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    Lattice_Change_Basic::largest_grad = Lattice_Change_Basic::largest_grad * unit_transform;
    stress_ii_max = stress_ii_max * unit_transform;

    if (Lattice_Change_Basic::largest_grad == 0.0)
    {
        GlobalV::ofs_running << " largest stress is 0, no movement is possible." << std::endl;
        GlobalV::ofs_running << " it may converged, otherwise no movement of lattice parameters is allowed."
                             << std::endl;
        Lattice_Change_Basic::converged = true;
    }
    else if (GlobalC::ucell.lc[0] == 1 && GlobalC::ucell.lc[1] == 1 && GlobalC::ucell.lc[2] == 1)
    {
        // if(Lattice_Change_Basic::largest_grad < GlobalV::STRESS_THR)
        if (Lattice_Change_Basic::largest_grad < GlobalV::STRESS_THR && stress_ii_max < GlobalV::STRESS_THR)
        {
            GlobalV::ofs_running << "\n Lattice relaxation is converged!" << std::endl;
            GlobalV::ofs_running << "\n Largest gradient is = " << largest_grad << std::endl;
            Lattice_Change_Basic::converged = true;
            ++Lattice_Change_Basic::update_iter;
        }
        else
        {
            GlobalV::ofs_running << "\n Lattice relaxation is not converged yet (threshold is " << GlobalV::STRESS_THR
                                 << ")" << std::endl;
            Lattice_Change_Basic::converged = false;
        }
    }
    else
    {
        /*for(int i=0; i<9; i++)
        {
            std::cout<<"i= "<<i<<" "<<grad[i]<<std::endl;
        }*/
        if (Lattice_Change_Basic::largest_grad < 10 * GlobalV::STRESS_THR)
        {
            GlobalV::ofs_running << "\n Lattice relaxation is converged!" << std::endl;
            GlobalV::ofs_running << "\n Largest gradient is = " << largest_grad << std::endl;
            Lattice_Change_Basic::converged = true;
            ++Lattice_Change_Basic::update_iter;
        }
        else
        {
            GlobalV::ofs_running << "\n Lattice relaxation is not converged yet (threshold is " << GlobalV::STRESS_THR
                                 << ")" << std::endl;
            Lattice_Change_Basic::converged = false;
        }
    }

    return;
}

void Lattice_Change_Basic::terminate(void)
{
    ModuleBase::TITLE("Lattice_Change_Basic", "terminate");
    if (Lattice_Change_Basic::converged)
    {
        GlobalV::ofs_running << " end of lattice optimization" << std::endl;
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "stress_step", Lattice_Change_Basic::stress_step);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "update iteration", Lattice_Change_Basic::update_iter);
        /*
        GlobalV::ofs_running<<"Saving the approximate inverse hessian"<<std::endl;
        std::ofstream hess("hess.out");
        for(int i=0;i<dim;i++)
        {
            for(int j=0;j<dim;j++)
            {
                hess << inv_hess(i,j);
            }
        }
        hess.close();
        */
    }
    else
    {
        GlobalV::ofs_running << " the maximum number of steps has been reached." << std::endl;
        GlobalV::ofs_running << " end of lattice optimization." << std::endl;
    }

    return;
}

void Lattice_Change_Basic::setup_etot(const double &energy_in, const bool judgement)
{
    if (Lattice_Change_Basic::stress_step == 1)
    {
        // p == previous
        Lattice_Change_Basic::etot_p = energy_in;
        Lattice_Change_Basic::etot = energy_in;
        ediff = etot - etot_p;
    }
    else
    {
        if (judgement)
        {
            Lattice_Change_Basic::etot = energy_in;
            if (Lattice_Change_Basic::etot_p > etot)
            {
                ediff = etot - etot_p;
                Lattice_Change_Basic::etot_p = etot;
            }
            else
            {
                // this step will not be accepted
                ediff = 0.0;
            }
        }
        else // for bfgs
        {
            Lattice_Change_Basic::etot_p = etot;
            Lattice_Change_Basic::etot = energy_in;
            ediff = etot - etot_p;
        }
    }

    return;
}
