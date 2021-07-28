#include "../src_lcao/local_orbital_charge.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"


void Local_Orbital_Charge::read_dm(const int &is, const string &fn)
{
    TITLE("Local_Orbital_Charge","read_dm");
    timer::tick("Local_Orbital_Charge","read_dm");

    GlobalV::ofs_running << "\n processor 0 is reading density matrix from file < " << fn << " > " << endl;
    //xiaohui modify 2015-03-25
    //bool quit_mesia = false;
    bool quit_abacus = false;

    ifstream ifs;
    if(GlobalV::MY_RANK==0)
    {
        ifs.open(fn.c_str());
        if (!ifs)
        {
            //xiaohui modify 2015-03-25
            //quit_mesia = true;
            quit_abacus = true;
        }
        else
        {
            // if the number is not match,
            // quit the program or not.
            bool quit=false;

            string name;
            ifs >> name;

            // check lattice constant, unit is Angstrom
            CHECK_DOUBLE(ifs,GlobalC::ucell.lat0 * BOHR_TO_A,quit);
            CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e11,quit);
            CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e12,quit);
            CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e13,quit);
            CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e21,quit);
            CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e22,quit);
            CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e23,quit);
            CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e31,quit);
            CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e32,quit);
            CHECK_DOUBLE(ifs,GlobalC::ucell.latvec.e33,quit);

            for(int it=0; it<GlobalC::ucell.ntype; it++)
            {
                CHECK_STRING(ifs,GlobalC::ucell.atoms[it].label,quit);
            }

            for(int it=0; it<GlobalC::ucell.ntype; it++)
            {
                CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].na,quit);
            }

            string coordinate;
            ifs >> coordinate;

            for(int it=0; it<GlobalC::ucell.ntype; it++)
            {
                for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
                {
                    CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].x,quit);
                    CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].y,quit);
                    CHECK_DOUBLE(ifs,GlobalC::ucell.atoms[it].taud[ia].z,quit);
                }
            }

            CHECK_INT(ifs, GlobalV::NSPIN);
            if(GlobalV::NSPIN == 1||GlobalV::NSPIN == 4)
            {
                READ_VALUE(ifs, GlobalC::en.ef);
                GlobalV::ofs_running << " read in fermi energy = " << GlobalC::en.ef << endl;
            }
            else if(GlobalV::NSPIN == 2)
            {
                if(is==0)READ_VALUE(ifs, GlobalC::en.ef_up);
                else if(is==1)READ_VALUE(ifs, GlobalC::en.ef_dw);
            }
            else
            {
                WARNING_QUIT("read_dm","check nspin!");
            }
            CHECK_INT(ifs, GlobalV::NLOCAL);
            CHECK_INT(ifs, GlobalV::NLOCAL);
        }// If file exist, read in data.
    } // Finish reading the first part of density matrix.


#ifndef __MPI
    GlobalV::ofs_running << " Read SPIN = " << is+1 << " density matrix now." << endl;

    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        for(int i=0; i<GlobalV::NLOCAL; ++i)
        {
            for(int j=0; j<GlobalV::NLOCAL; ++j)
            {
                ifs >> DM[is][i][j];
            }
        }
    }
    else
    {
        WARNING_QUIT("Local_Orbital_Charge::read_dm","The nnrg should not be update");
        CHECK_INT(ifs,GlobalC::LNNR.nnrg);

        for(int i=0; i<GlobalC::LNNR.nnrg; ++i)
        {
            ifs >> DM_R[is][i];
        }
    }
#else

    // distribution of necessary data
    //xiaohui modify 2015-03-25
    //Parallel_Common::bcast_bool(quit_mesia);
    Parallel_Common::bcast_bool(quit_abacus);
    //xiaohui modify 2015-03-25
    //if(quit_mesia)
    if(quit_abacus)
    {
        WARNING_QUIT("Local_Orbital_Charge::read_dm","Can not find the density matrix file.");
    }


    if(GlobalV::NSPIN==1||GlobalV::NSPIN==4)
    {
        Parallel_Common::bcast_double(GlobalC::en.ef);
    }
    else if(GlobalV::NSPIN==2)
    {
        Parallel_Common::bcast_double(GlobalC::en.ef_up);
        Parallel_Common::bcast_double(GlobalC::en.ef_dw);
    }


    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        //GlobalV::ofs_running << " NLOCAL=" << GlobalV::NLOCAL << endl;
        //GlobalV::ofs_running << " lgd_now=" << lgd_now << endl;
        //GlobalV::ofs_running << " GlobalC::GridT.lgd=" << GlobalC::GridT.lgd << endl;

        double *tmp = new double[GlobalV::NLOCAL];
        for(int i=0; i<GlobalV::NLOCAL; ++i)
        {
            //GlobalV::ofs_running << " i=" << i << endl;
            ZEROS(tmp, GlobalV::NLOCAL);
            if(GlobalV::MY_RANK==0)
            {
                for(int j=0; j<GlobalV::NLOCAL; ++j)
                {
                    ifs >> tmp[j];
                }
            }
            Parallel_Common::bcast_double(tmp, GlobalV::NLOCAL);

            const int mu = GlobalC::GridT.trace_lo[i];
            if(mu >= 0)
            {   
                for(int j=0; j<GlobalV::NLOCAL; ++j)
                {
                    const int nu = GlobalC::GridT.trace_lo[j];
                    if(nu >= 0)
                    {
                        DM[is][mu][nu] = tmp[j];
                    }
                }
            }
        }// i
        delete[] tmp;
    }
    else
    {
        WARNING_QUIT("Local_Orbital_Charge::read_dm","not ready to readin DM_R");
    }
#endif
    if(GlobalV::MY_RANK==0) ifs.close();

    GlobalV::ofs_running << " Finish reading density matrix." << endl;

    timer::tick("Local_Orbital_Charge","read_dm");
    return;
}
