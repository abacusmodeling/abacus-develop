//wenfei 2022-1-5
//This file contains subroutines that deals with io,
//as well as interface with libnpy, since many arrays must be saved in numpy format

//There are 6 subroutines in this file,
//two of which prints in human-readable formats:
//1. print_descriptor : prints descriptors into descriptor.dat
//2. print_F_delta : prints F_delta into F_delta.dat

//four of which prints in numpy formats:
//3. save_npy_d : descriptor ->dm_eig.npy
//4. save_npy_gvx : gvx ->grad_vx.npy
//5. save_npy_e : energy
//6. save_npy_f : force
#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "npy.hpp"

//prints descriptor into descriptor.dat
void LCAO_Deepks::print_descriptor(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "print_descriptor");
    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"
       << "descriptor.dat";
    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }
    for (int ia = 0; ia < nat; ia++)
    {
        ofs << " atom_index " << ia + 1 << " n_descriptor " << this->des_per_atom << std::endl;
        for(int inl=0;inl<inlmax/nat;inl++)
        {
            int nm = 2*inl_l[inl]+1;
            for(int im=0;im<nm;im++)
            {
                const int ind=ia*inlmax/nat+inl;
                ofs << std::setprecision(10) << d_tensor[ind].index({im}).item().toDouble() << " ";
            }
        }
        ofs << std::endl << std::endl;
    }
    GlobalV::ofs_running << " Descriptors have been printed to " << ss.str() << std::endl;

    return;
}

//prints F_delta into F_delta.dat
void LCAO_Deepks::print_F_delta(const string& fname, const UnitCell_pseudo &ucell)
{
    ModuleBase::TITLE("LCAO_Deepks", "print_F_delta");

    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"<< fname ;

    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }

    ofs << "F_delta(Hatree/Bohr) from deepks model: " << std::endl;
    ofs << std::setw(12) << "type" << std::setw(12) << "atom" << std::setw(15) << "dF_x" << std::setw(15) << "dF_y" << std::setw(15) << "dF_z" << std::endl;

    for (int it = 0;it < ucell.ntype;++it)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;++ia)
        {
            int iat = ucell.itia2iat(it, ia);
            ofs << std::setw(12) << ucell.atoms[it].label << std::setw(12) << ia
                << std::setw(15) << this->F_delta(iat, 0) / 2 << std::setw(15) << this->F_delta(iat, 1) / 2
                << std::setw(15) << this->F_delta(iat, 2) / 2 << std::endl;
        }
    }

    ofs << "F_delta(eV/Angstrom) from deepks model: " << std::endl;
    ofs << std::setw(12) << "type" << std::setw(12) << "atom" << std::setw(15) << "dF_x" << std::setw(15) << "dF_y" << std::setw(15) << "dF_z" << std::endl;

    for (int it = 0;it < ucell.ntype;++it)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;++ia)
        {
            int iat = ucell.itia2iat(it, ia);
            ofs << std::setw(12) << ucell.atoms[it].label << std::setw(12)
                << ia << std::setw(15) << this->F_delta(iat, 0) * ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A
                << std::setw(15) << this->F_delta(iat, 1) * ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A
                << std::setw(15) << this->F_delta(iat, 2) * ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A << std::endl;
        }
    }

    GlobalV::ofs_running << " F_delta has been printed to " << ss.str() << std::endl;
    ofs.close();

    return;
}

//saves descriptor into dm_eig.npy
void LCAO_Deepks::save_npy_d(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_d");
    //save descriptor in .npy format
    vector<double> npy_des;
    for (int inl = 0;inl < inlmax;++inl)
    {
        int nm = 2*inl_l[inl] + 1;
        for(int im=0;im<nm;im++)
        {
            npy_des.push_back(this->d_tensor[inl].index({im}).item().toDouble());
        }
    }
    const long unsigned dshape[] = {(long unsigned) nat, (long unsigned) this->des_per_atom };
    if (GlobalV::MY_RANK == 0)
    {
        npy::SaveArrayAsNumpy("dm_eig.npy", false, 2, dshape, npy_des);
    }
    return;
}

//saves gvx into grad_vx.npy
void LCAO_Deepks::save_npy_gvx(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_gvx");
    //save grad_vx.npy (when  force label is in use)
    //unit: /Bohr
    const long unsigned gshape[] = {(long unsigned) nat, 3, nat, this->des_per_atom};
    vector<double> npy_gvx;
    for (int ibt = 0;ibt < nat;++ibt)
    {
        for (int i = 0;i < 3;i++)
        {
            for (int iat = 0;iat < nat;++iat)
            {
                for(int p=0;p<this->des_per_atom;++p)
                {
                    npy_gvx.push_back(this->gvx_tensor.index({ ibt, i, iat, p }).item().toDouble());
                }
            }
        }
    }
    npy::SaveArrayAsNumpy("grad_vx.npy", false, 4, gshape, npy_gvx);
    return;
}

//saves energy in numpy format
void LCAO_Deepks::save_npy_e(const double &e, const std::string &e_file)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_e");
    //save e_base
    const long unsigned eshape[] = { 1 };
    vector<double> npy_e;
    npy_e.push_back(e);
    npy::SaveArrayAsNumpy(e_file, false, 1, eshape, npy_e);
    return;
}

//saves force in numpy format
void LCAO_Deepks::save_npy_f(const ModuleBase::matrix &f, const std::string &f_file, const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_f");
    //save f_base
    //caution: unit: Rydberg/Bohr
    const long unsigned fshape[] = {(long unsigned) nat, 3 };
    vector<double> npy_f;
    for (int iat = 0;iat < nat;++iat)
    {
        for (int i = 0;i < 3;i++)
        {
            npy_f.push_back(f(iat, i));
        }
    }
    npy::SaveArrayAsNumpy(f_file, false, 2, fshape, npy_f);
    return;
}

#endif