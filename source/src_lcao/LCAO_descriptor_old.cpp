//keep some no longer used deepks subroutines here
#ifdef __DEEPKS

#include "LCAO_descriptor.h"
#include "LCAO_matrix.h"
#include "../module_base/lapack_connector.h"
#include "../module_base/intarray.h"
#include "../module_base/complexmatrix.h"
#include "global_fp.h"
#include "../src_pw/global.h"
#include "../src_io/winput.h"

#include <torch/script.h>
#include <torch/csrc/autograd/autograd.h>
#include <npy.hpp>
#include <torch/csrc/api/include/torch/linalg.h>

void LCAO_Descriptor::cal_v_delta(const ModuleBase::matrix& dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_v_delta");
    //1.  (dE/dD)<alpha_m'|psi_nv> (descriptor changes in every scf iter)
    this->cal_gedm(dm);
    
    //2. multiply overlap matrice and sum
    double* tmp_v1 = new double[(2 * lmaxd + 1) * GlobalV::NLOCAL];
    double* tmp_v2 = new double[GlobalV::NLOCAL *GlobalV::NLOCAL];

    ModuleBase::GlobalFunc::ZEROS(this->H_V_delta, GlobalV::NLOCAL * GlobalV::NLOCAL); //init before calculate
    
    for (int inl = 0;inl < inlmax;inl++)
    {
        ModuleBase::GlobalFunc::ZEROS(tmp_v1, (2 * lmaxd + 1) * GlobalV::NLOCAL);
        ModuleBase::GlobalFunc::ZEROS(tmp_v2, GlobalV::NLOCAL * GlobalV::NLOCAL);
        int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
        const char t = 'T';  //transpose
        const char nt = 'N'; //non transpose
        const double alpha = 1;
        const double beta = 0;
        double* a = this->gedm[inl];//[nm][nm]
        double* b = S_mu_alpha[inl];//[GlobalV::NLOCAL][nm]--trans->[nm][GlobalV::NLOCAL]
        double* c = tmp_v1;
        
        //2.1  (dE/dD)*<alpha_m'|psi_nv>
        dgemm_(&nt, &t, &nm, &GlobalV::NLOCAL, &nm, &alpha, a, &nm, b, &GlobalV::NLOCAL, &beta, c, &nm);

        //2.2  <psi_mu|alpha_m>*(dE/dD)*<alpha_m'|psi_nv>
        a = b; //[GlobalV::NLOCAL][nm]
        b = c;//[nm][GlobalV::NLOCAL]
        c = tmp_v2;//[GlobalV::NLOCAL][GlobalV::NLOCAL]
        dgemm_(&nt, &nt, &GlobalV::NLOCAL, &GlobalV::NLOCAL, &nm, &alpha, a, &GlobalV::NLOCAL, b, &nm, &beta, c, &GlobalV::NLOCAL);

        //3. sum of Inl
        for (int i = 0;i < GlobalV::NLOCAL * GlobalV::NLOCAL;++i)
        {
            this->H_V_delta[i] += c[i];
        }
    }
    delete[] tmp_v1;
    delete[] tmp_v2;

    GlobalV::ofs_running << " Finish calculating H_V_delta" << std::endl;
    return;
}

// compute the full projected density matrix for each atom
// save the matrix for each atom in order to minimize the usage of memory
// --mohan 2021-08-04
void LCAO_Descriptor::cal_dm_as_descriptor(const ModuleBase::matrix &dm)
{
	ModuleBase::TITLE("LCAO_Descriptor", "cal_proj_dm");

	for(int it=0; it<GlobalC::ucell.ntype; ++it)
	{
		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ++ia)
		{
			// compute S^T * dm * S to obtain descriptor
			// of each atom 
			// and then diagonalize it
		}
	}

	return;
}

void LCAO_Descriptor::cal_f_delta(const ModuleBase::matrix &dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_f_delta");
    this->F_delta.zero_out();
    //1. cal gedm
    this->cal_gedm(dm);
    
    //2. cal gdmx
    this->init_gdmx();
    this->cal_gdmx(dm);
    
    //3.multiply and sum for each atom
    //3.1 Pulay term 
    // \sum_{Inl}\sum_{mm'} <gedm, gdmx>_{mm'}
    //notice: sum of multiplied corresponding element(mm') , not matrix multiplication !
    int iat = 0;    //check if the index same as GlobalC::ucell.iw2iat or not !!
    for (int it = 0;it < GlobalC::ucell.ntype;++it)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;++ia)
        {
            for (int inl = 0;inl < this->inlmax;++inl)
            {
                int nm = 2 * inl_l[inl] + 1;
                for (int m1 = 0;m1 < nm;++m1)
                {
                    for (int m2 = 0; m2 < nm;++m2)
                    {
                        this->F_delta(iat, 0) += this->gedm[inl][m1 * nm + m2] * gdmx[iat][inl][m1 * nm + m2];
                        this->F_delta(iat, 1) += this->gedm[inl][m1 * nm + m2] * gdmy[iat][inl][m1 * nm + m2];
                        this->F_delta(iat, 2) += this->gedm[inl][m1 * nm + m2] * gdmz[iat][inl][m1 * nm + m2];
                    }
                }
            }//end inl
            ++iat;
        }
    }
    this->print_F_delta("F_delta_pulay_old.dat");
    this->F_delta.zero_out();
    iat = 0;
    for (int it = 0;it < GlobalC::ucell.ntype;++it)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;++ia)
        {
            //3.2 HF term
            double** ss = this->S_mu_alpha;
            double** dsx = this->DS_mu_alpha_x;
            double** dsy = this->DS_mu_alpha_y;
            double** dsz = this->DS_mu_alpha_z;
            for (int mu = 0;mu < GlobalV::NLOCAL;++mu)
            {
                for (int nu = 0;nu < GlobalV::NLOCAL;++nu)
                {
                    for (int l = 0;l <= GlobalC::ORB.Alpha[0].getLmax();++l)
                    {
                        for (int n = 0;n < GlobalC::ORB.Alpha[0].getNchi(l);++n)
                        {
                            for (int m1 = 0;m1 < 2 * l + 1;++m1)
                            {
                                for (int m2 = 0;m2 < 2 * l + 1;++m2)
                                {
                                    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx")
                                    {
                                        this->F_delta(iat, 0) -= 2*dm(mu, nu) * dsx[inl_index[it](ia, l, n)][m1 * GlobalV::NLOCAL + mu]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][m2 * GlobalV::NLOCAL + nu];
                                        this->F_delta(iat, 1) -= 2*dm(mu, nu) * dsy[inl_index[it](ia, l, n)][m1 * GlobalV::NLOCAL + mu]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][m2 * GlobalV::NLOCAL + nu];
                                        this->F_delta(iat, 2) -= 2*dm(mu, nu) * dsz[inl_index[it](ia, l, n)][m1 * GlobalV::NLOCAL + mu]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][m2 * GlobalV::NLOCAL + nu];
                                    }
                                    else
                                    {
                                        this->F_delta(iat, 0) -= 2*dm(mu, nu) * dsx[inl_index[it](ia, l, n)][mu* (2*l+1) + m1]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][nu* (2*l+1) + m2];
                                        this->F_delta(iat, 1) -= 2*dm(mu, nu) * dsy[inl_index[it](ia, l, n)][mu* (2*l+1) + m1]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][nu* (2*l+1) + m2];
                                        this->F_delta(iat, 2) -= 2*dm(mu, nu) * dsz[inl_index[it](ia, l, n)][mu* (2*l+1) + m1]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][nu* (2*l+1) + m2];
                                    }
                                }//end m2
                            }//end m1
                        }//end n
                    }//end l
                }//end nu
            }//end mu
            ++iat;
        }//end ia
    }//end it
    this->print_F_delta("F_delta_hf_old.dat");
    //3.3 Overlap term
    //somthing in NN, which not included in Hamiltonian
    /*
    for (int mu = 0;mu < GlobalV::NLOCAL;++mu)
    {
        const int iat = GlobalC::ucell.iwt2iat[mu];
        for (int nu = 0;nu < GlobalV::NLOCAL;++nu)
        {
            this->F_delta(iat, 0) += 2*(this->E_delta - this->e_delta_band)* dm(mu, nu) * GlobalC::LM.DSloc_x[mu * GlobalV::NLOCAL + nu];
            this->F_delta(iat, 1) += 2*(this->E_delta - this->e_delta_band) * dm(mu, nu) * GlobalC::LM.DSloc_y[mu * GlobalV::NLOCAL + nu];
            this->F_delta(iat, 2) += 2*(this->E_delta - this->e_delta_band) * dm(mu, nu) * GlobalC::LM.DSloc_z[mu * GlobalV::NLOCAL + nu];
        }
    }*/
    this->del_gdmx();
    return;
}
#endif