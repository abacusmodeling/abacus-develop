//wenfei add 2021 october
#ifdef __DEEPKS

#include "LCAO_descriptor.h"

//===============================
//DeePKS Part 2
//deals with application of correction dV to Hamiltonian and force
//using a new algorithm similar to the nonlocal pseudoopotential,
//that allows for parallelization, multi-k calculations, and is more efficient
//for some of them, the older version can be found in LCAO_descriptor_old
//===============================


void LCAO_Descriptor::deepks_pre_scf(const string& model_file)
{
    ModuleBase::TITLE("LCAO_Descriptor", "deepks_pre_scf");

	// load the DeePKS model from deep neural network
    this->load_model(model_file);
    
    //initialize the H matrix H_V_delta
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        delete[] this->H_V_delta;
        this->H_V_delta = new double[GlobalC::ParaO.nloc];
        ModuleBase::GlobalFunc::ZEROS(this->H_V_delta, GlobalC::ParaO.nloc);
    }
    else
    {
        H_V_delta_k = new std::complex<double>* [GlobalC::kv.nks];
        for(int ik=0;ik<GlobalC::kv.nks;ik++)
        {
            this->H_V_delta_k[ik] = new std::complex<double>[GlobalC::ParaO.nloc];
            ModuleBase::GlobalFunc::ZEROS(this->H_V_delta_k[ik], GlobalC::ParaO.nloc);
        }
    }

    //init gedm**
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    this->gedm = new double* [this->inlmax];
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->gedm[inl] = new double[pdm_size];
        ModuleBase::GlobalFunc::ZEROS(this->gedm[inl], pdm_size);
    }
    if (GlobalV::FORCE)
    {
        //init F_delta
        F_delta.create(GlobalC::ucell.nat, 3);
        //init DS_mu_alpha**
        this->DS_mu_alpha_x = new double* [this->inlmax];
        this->DS_mu_alpha_y = new double* [this->inlmax];
        this->DS_mu_alpha_z = new double* [this->inlmax];
        for (int inl = 0;inl < this->inlmax;inl++)
        {
            this->DS_mu_alpha_x[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];
            this->DS_mu_alpha_y[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];
            this->DS_mu_alpha_z[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];
            ModuleBase::GlobalFunc::ZEROS(DS_mu_alpha_x[inl], GlobalV::NLOCAL * (2 * this->lmaxd + 1));
            ModuleBase::GlobalFunc::ZEROS(DS_mu_alpha_y[inl], GlobalV::NLOCAL * (2 * this->lmaxd + 1));
            ModuleBase::GlobalFunc::ZEROS(DS_mu_alpha_z[inl], GlobalV::NLOCAL * (2 * this->lmaxd + 1));
        }
        //init DH_V_delta*
        if(GlobalV::GAMMA_ONLY_LOCAL)
        {
            delete[] DH_V_delta_x;
            delete[] DH_V_delta_y;
            delete[] DH_V_delta_z;
            this->DH_V_delta_x = new double[GlobalC::ParaO.nloc];
            this->DH_V_delta_y = new double [GlobalC::ParaO.nloc];
            this->DH_V_delta_z = new double[GlobalC::ParaO.nloc];
            ModuleBase::GlobalFunc::ZEROS(DH_V_delta_x, GlobalC::ParaO.nloc);
            ModuleBase::GlobalFunc::ZEROS(DH_V_delta_y, GlobalC::ParaO.nloc);
            ModuleBase::GlobalFunc::ZEROS(DH_V_delta_z, GlobalC::ParaO.nloc);
        }
        else
        {
            DH_V_delta_x_k = new std::complex<double>* [GlobalC::kv.nks];
            DH_V_delta_y_k = new std::complex<double>* [GlobalC::kv.nks];
            DH_V_delta_z_k = new std::complex<double>* [GlobalC::kv.nks];
            for(int ik=0;ik<GlobalC::kv.nks;ik++)
            {
                this->DH_V_delta_x_k[ik] = new std::complex<double>[GlobalC::ParaO.nloc];
                this->DH_V_delta_y_k[ik] = new std::complex<double>[GlobalC::ParaO.nloc];
                this->DH_V_delta_z_k[ik] = new std::complex<double>[GlobalC::ParaO.nloc];
                ModuleBase::GlobalFunc::ZEROS(this->DH_V_delta_x_k[ik], GlobalC::ParaO.nloc);
                ModuleBase::GlobalFunc::ZEROS(this->DH_V_delta_y_k[ik], GlobalC::ParaO.nloc);
                ModuleBase::GlobalFunc::ZEROS(this->DH_V_delta_z_k[ik], GlobalC::ParaO.nloc);
            }           
        }
    }
    if(!GlobalV::GAMMA_ONLY_LOCAL)
    {
        this->nlm_k.resize(GlobalC::ucell.nat);
    }
    return;
}

//this subroutine adds dV to the Kohn-Sham Hamiltonian
//for gamma_only calculations
void LCAO_Descriptor::add_v_delta(void)
{
    ModuleBase::TITLE("LCAO_DESCRIPTOR", "add_v_delta");

    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        for (int iw1 = 0;iw1 < GlobalV::NLOCAL;++iw1)
        {
            for (int iw2 = 0;iw2 < GlobalV::NLOCAL;++iw2)
            {
				if (!GlobalC::ParaO.in_this_processor(iw1,iw2))
				{
					continue;
				}
                
                int iw1_local=GlobalC::ParaO.trace_loc_row[iw1];
                int iw2_local=GlobalC::ParaO.trace_loc_col[iw2];
                int index = iw2_local * GlobalC::ParaO.nrow+ iw1_local;
                GlobalC::LM.set_HSgamma(iw1, iw2, this->H_V_delta[index], 'L');
            }
        }
    }
    else
    {
		ModuleBase::WARNING_QUIT("add_v_delta","should not be used for multi-k; use add_v_delta_k instead");
        //call set_HSk, complex Matrix
    }
	return;
}

//this subroutine adds dV to the Kohn-Sham Hamiltonian
//for multi-k calculations
void LCAO_Descriptor::add_v_delta_k(const int &ik)
{
    ModuleBase::TITLE("LCAO_DESCRIPTOR", "add_v_delta_k");
    ModuleBase::GlobalFunc::ZEROS(this->H_V_delta_k[ik], GlobalC::ParaO.nloc);

    const double Rcut_Alpha = GlobalC::ORB.Alpha[0].getRcut();

    for (int T0 = 0; T0 < GlobalC::ucell.ntype; T0++)
    {
		Atom* atom0 = &GlobalC::ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = GlobalC::ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GlobalC::GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GlobalC::GridD.getType(ad1);
                const int I1 = GlobalC::GridD.getNatom(ad1);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &GlobalC::ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

                ModuleBase::Vector3<double> dR1(GlobalC::GridD.getBox(ad1).x, GlobalC::GridD.getBox(ad1).y, GlobalC::GridD.getBox(ad1).z); 

				for (int ad2=0; ad2 < GlobalC::GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GlobalC::GridD.getType(ad2);
					const int I2 = GlobalC::GridD.getNatom(ad2);
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &GlobalC::ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    ModuleBase::Vector3<double> dR2(GlobalC::GridD.getBox(ad2).x, GlobalC::GridD.getBox(ad2).y, GlobalC::GridD.getBox(ad2).z);
					
					const double Rcut_AO2 = GlobalC::ORB.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * GlobalC::ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

                    const double arg = ( GlobalC::kv.kvec_d[ik] * (dR2-dR1) ) * ModuleBase::TWO_PI;
                    const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1;
						const int iw1_local = GlobalC::ParaO.trace_loc_row[iw1_all];
						if(iw1_local < 0)continue;
						const int iw1_0 = iw1/GlobalV::NPOL;
						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2;
							const int iw2_local = GlobalC::ParaO.trace_loc_col[iw2_all];
							if(iw2_local < 0)continue;
							const int iw2_0 = iw2/GlobalV::NPOL;

                            double nlm=0.0;
                            std::vector<double> nlm1 = this->nlm_k[iat][ad1][iw1_all][0];
                            std::vector<double> nlm2 = this->nlm_k[iat][ad2][iw2_all][0];

                            assert(nlm1.size()==nlm2.size());
                            int ib=0;
                            for (int L0 = 0; L0 <= GlobalC::ORB.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < GlobalC::ORB.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                    {
                                        for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                        {
                                            nlm += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[ib+m2];
                                        }
                                    }
                                    ib+=(2*L0+1);
                                }
                            }

                            assert(ib==nlm1.size());
                            
                            int iic;

                            if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                            {
                                iic=iw1_local+iw2_local*GlobalC::ParaO.nrow;
                            }
                            else
                            {
                                iic=iw1_local*GlobalC::ParaO.ncol+iw2_local;
                            }
                            this->H_V_delta_k[ik][iic] += nlm*kphase;
                            GlobalC::LM.Hloc2[iic] += nlm * kphase;

						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }
    }
}

//the pulay term in force
//for gamma_only calculations
void LCAO_Descriptor::cal_f_delta_pulay(const ModuleBase::matrix& dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_f_delta_pulay");
    //this->F_delta.zero_out();
    this->build_v_delta_alpha_new(1);
    //this->build_v_delta_mu(1);    //, if multi-k
    for (int i = 0;i < GlobalV::NLOCAL;++i)   //col, diff
    {
        const int iat = GlobalC::ucell.iwt2iat[i];//the atom whose force being calculated
        for (int j = 0;j < GlobalV::NLOCAL;++j)   //row
        {
            const int mu = GlobalC::ParaO.trace_loc_row[j];
            const int nu = GlobalC::ParaO.trace_loc_col[i];
            if (mu >= 0 && nu >= 0)
            {
                const int index = mu * GlobalC::ParaO.ncol + nu;
                this->F_delta(iat, 0) += 2.0 * dm(nu, mu) * this->DH_V_delta_x[index];
                this->F_delta(iat, 1) += 2.0 * dm(nu, mu) * this->DH_V_delta_y[index];
                this->F_delta(iat, 2) += 2.0 * dm(nu, mu) * this->DH_V_delta_z[index];
            }
        }
    }
    return;
}

//the pulay term in force
//for multi_k calculations
void LCAO_Descriptor::cal_f_delta_k_pulay(const std::vector<ModuleBase::ComplexMatrix>& dm/**<[in] density matrix*/)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_f_delta_k_pulay");
    this->build_v_delta_alpha_new(1);
    //this->build_v_delta_mu(1);    //, if multi-k
    //GlobalV::ofs_running << "dim of dm : " << dm[0].nr << " " << dm[0].nc << std::endl;

    for (int i = 0;i < GlobalV::NLOCAL;++i)   //col, diff
    {
        const int iat = GlobalC::ucell.iwt2iat[i];//the atom whose force being calculated
        for (int j = 0;j < GlobalV::NLOCAL;++j)   //row
        {
            const int mu = GlobalC::ParaO.trace_loc_row[j];
            const int nu = GlobalC::ParaO.trace_loc_col[i];
            if (mu >= 0 && nu >= 0)
            {
                const int index = mu * GlobalC::ParaO.ncol + nu;
                //GlobalV::ofs_running << "mu,nu,index : " << mu << " " << nu << " " << index << std::endl;
                
                std::complex<double> f[3];
                for(int ik=0;ik<GlobalC::kv.nks;ik++)
                {
                    f[0]+= 2.0 * dm[ik](nu, mu) * this->DH_V_delta_x_k[ik][index];
                    f[1]+= 2.0 * dm[ik](nu, mu) * this->DH_V_delta_y_k[ik][index];
                    f[2]+= 2.0 * dm[ik](nu, mu) * this->DH_V_delta_z_k[ik][index];
                }
                if(f[0].imag()>1e-12 || f[1].imag()>1e-12 || f[2].imag()>1e-12)
                {
                    GlobalV::ofs_running << "f_xyz : " << f[0] << " " << f[1] << " " << f[2] << std::endl;
                    ModuleBase::WARNING_QUIT("cal_f_delta_hf_k","Force should be real!");
                }
                this->F_delta(iat, 0) += f[0].real();
                this->F_delta(iat, 1) += f[1].real();
                this->F_delta(iat, 2) += f[2].real();
            }
        }
    }
    return;

}

//the hellman-feynmann term in force
//for gamma_only calculations
void LCAO_Descriptor::cal_f_delta_hf_new(const ModuleBase::matrix& dm, const bool isstress, ModuleBase::matrix& svnl_dalpha)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_f_delta_hf_new");
    this->F_delta.zero_out();

    const double Rcut_Alpha = GlobalC::ORB.Alpha[0].getRcut();
    for (int T0 = 0; T0 < GlobalC::ucell.ntype; T0++)
    {
		Atom* atom0 = &GlobalC::ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            int iat = GlobalC::ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom0->tau[I0] ,T0, I0);

            //=======================================================
            //Step1:	
            //saves <beta|psi>, where beta runs over L0,M0 on atom I0
            //and psi runs over atomic basis sets on the current core
            //=======================================================
            std::vector<std::unordered_map<int,std::vector<std::vector<double>>>> nlm_tot;

            //outermost loop : all adjacent atoms
            nlm_tot.resize(GlobalC::GridD.getAdjacentNum()+1);

            for (int ad=0; ad<GlobalC::GridD.getAdjacentNum()+1 ; ++ad)
            {
                const int T1 = GlobalC::GridD.getType(ad);
                const int I1 = GlobalC::GridD.getNatom(ad);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad);
                const Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int nw1_tot = atom1->nw*GlobalV::NPOL;

                //middle loop : atomic basis on current processor (either row or column)
                nlm_tot[ad].clear();

                const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;

                if (dist1 > Rcut_Alpha + Rcut_AO1)
                {
                    continue;
                }

                for (int iw1=0; iw1<nw1_tot; ++iw1)
                {                   
                    const int iw1_all = start1 + iw1;
                    const int iw1_local = GlobalC::ParaO.trace_loc_row[iw1_all];
                    const int iw2_local = GlobalC::ParaO.trace_loc_col[iw1_all];
                    if(iw1_local < 0 && iw2_local < 0)continue;
                    const int iw1_0 = iw1/GlobalV::NPOL;
                    std::vector<std::vector<double>> nlm;
                    //2D, but first dimension is only 1 here
                    //for force, the right hand side is the gradient
                    //and the first dimension is then 3
                    //inner loop : all projectors (L0,M0)
					GlobalC::UOT.snap_psialpha_half(
						nlm, 1, tau1, T1,
						atom1->iw2l[ iw1_0 ], // L1
						atom1->iw2m[ iw1_0 ], // m1
						atom1->iw2n[ iw1_0 ], // N1
						GlobalC::ucell.atoms[T0].tau[I0], T0, I0, this->inl_index); //R0,T0

                    nlm_tot[ad].insert({iw1_all,nlm});
                }//end iw
            }//end ad

            //=======================================================
            //Step2:	
            //calculate sum_(L0,M0) alpha<psi_i|alpha><alpha|psi_j>
            //and accumulate the value to Hloc_fixed(i,j)
            //=======================================================

            for (int ad1=0; ad1<GlobalC::GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GlobalC::GridD.getType(ad1);
                const int I1 = GlobalC::GridD.getNatom(ad1);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad1);
                const Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int nw1_tot = atom1->nw*GlobalV::NPOL;
                const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

                for (int ad2=0; ad2 < GlobalC::GridD.getAdjacentNum()+1 ; ad2++)
                {
                    const int T2 = GlobalC::GridD.getType(ad2);
                    const int I2 = GlobalC::GridD.getNatom(ad2);
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                    const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau(ad2);
                    const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                    const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    
                    const double Rcut_AO2 = GlobalC::ORB.Phi[T2].getRcut();
                    const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                    const double dist2 = (tau2-tau0).norm() * GlobalC::ucell.lat0;

                    if (dist1 > Rcut_Alpha + Rcut_AO1
                            || dist2 > Rcut_Alpha + Rcut_AO2)
                    {
                        continue;
                    }					

                    double r0[3];
                    double r1[3];
                    if(isstress)
                    {
                        r1[0] = ( tau1.x - tau0.x) ;
                        r1[1] = ( tau1.y - tau0.y) ;
                        r1[2] = ( tau1.z - tau0.z) ;
                        r0[0] = ( tau2.x - tau0.x) ;
                        r0[1] = ( tau2.y - tau0.y) ;
                        r0[2] = ( tau2.z - tau0.z) ;
                    }

                    for (int iw1=0; iw1<nw1_tot; ++iw1)
                    {
                        const int iw1_all = start1 + iw1;
                        const int iw1_local = GlobalC::ParaO.trace_loc_col[iw1_all];
                        if(iw1_local < 0)continue;

                        for (int iw2=0; iw2<nw2_tot; ++iw2)
                        {
                            const int iw2_all = start2 + iw2;
                            const int iw2_local = GlobalC::ParaO.trace_loc_row[iw2_all];
                            if(iw2_local < 0)continue;

                            double nlm[3]={0,0,0};
                            double nlm_t[3] = {0,0,0}; //for stress
                            std::vector<double> nlm1 = nlm_tot[ad1][iw1_all][0];
                            std::vector<std::vector<double>> nlm2;
                            nlm2.resize(3);

                            for(int dim=0;dim<3;dim++)
                            {
                                nlm2[dim] = nlm_tot[ad2][iw2_all][dim+1];
                            }

                            assert(nlm1.size()==nlm2[0].size());

                            int ib=0;
                            for (int L0 = 0; L0 <= GlobalC::ORB.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < GlobalC::ORB.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < nm; ++m1)
                                    {
                                        for (int m2 = 0; m2 < nm; ++m2)
                                        {
                                            for(int dim=0;dim<3;dim++)
                                            {                                            
                                                nlm[dim] += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[dim][ib+m2];
                                            }
                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());

                            // HF term is minus, only one projector for each atom force.
                            this->F_delta(iat, 0) -= 2 * dm(iw1_local, iw2_local) * nlm[0];
                            this->F_delta(iat, 1) -= 2 * dm(iw1_local, iw2_local) * nlm[1];
                            this->F_delta(iat, 2) -= 2 * dm(iw1_local, iw2_local) * nlm[2];

                            if(isstress)
                            {
                                nlm1 = nlm_tot[ad2][iw2_all][0];
                                for(int i=0;i<3;i++)
                                {
                                    nlm2[i] = nlm_tot[ad1][iw1_all][i+1];
                                }

                                assert(nlm1.size()==nlm2[0].size());                                

                                int ib=0;
                                for (int L0 = 0; L0 <= GlobalC::ORB.Alpha[0].getLmax();++L0)
                                {
                                    for (int N0 = 0;N0 < GlobalC::ORB.Alpha[0].getNchi(L0);++N0)
                                    {
                                        const int inl = this->inl_index[T0](I0, L0, N0);
                                        const int nm = 2*L0+1;
                                        for (int m1 = 0;m1 < nm; ++m1)
                                        {
                                            for (int m2 = 0; m2 < nm; ++m2)
                                            {
                                                for(int dim=0;dim<3;dim++)
                                                {                                            
                                                    nlm_t[dim] += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[dim][ib+m2];
                                                }
                                            }
                                        }
                                        ib+=nm;
                                    }
                                }
                                assert(ib==nlm1.size());

                                for(int ipol=0;ipol<3;ipol++)
                                {
                                    svnl_dalpha(0,ipol) -= dm(iw1_local, iw2_local) * (nlm[0] * r0[ipol] + nlm_t[0] * r1[ipol])* -1;
                                    svnl_dalpha(1,ipol) -= dm(iw1_local, iw2_local) * (nlm[1] * r0[ipol] + nlm_t[1] * r1[ipol])* -1;
                                    svnl_dalpha(2,ipol) -= dm(iw1_local, iw2_local) * (nlm[2] * r0[ipol] + nlm_t[2] * r1[ipol])* -1;
                                }
                            }
                        }//iw2
                    }//iw1
                }//ad2
            }//ad1
        }//end I0
    }//end T0

    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                svnl_dalpha(i,j) *=  GlobalC::ucell.lat0 / GlobalC::ucell.omega;
            }
        }
    }
}

//the hellman-feynmann term in force
//for multi-k calculations

void LCAO_Descriptor::cal_f_delta_hf_k_new(const std::vector<ModuleBase::ComplexMatrix>& dm/**<[in] density matrix*/, const bool isstress, ModuleBase::matrix& svnl_dalpha)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_f_delta_hf_k_new");
    this->F_delta.zero_out();
    ModuleBase::ComplexMatrix F_delta_k;
    F_delta_k.create(this->F_delta.nr,this->F_delta.nc);

    ModuleBase::ComplexMatrix svnl_dalpha_k;
    svnl_dalpha_k.create(3,3);

    const double Rcut_Alpha = GlobalC::ORB.Alpha[0].getRcut();
    for (int T0 = 0; T0 < GlobalC::ucell.ntype; T0++)
    {
		Atom* atom0 = &GlobalC::ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            int iat = GlobalC::ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom0->tau[I0] ,T0, I0);

            //=======================================================
            //Step1:	
            //saves <beta|psi>, where beta runs over L0,M0 on atom I0
            //and psi runs over atomic basis sets on the current core
            //=======================================================
            std::vector<std::unordered_map<int,std::vector<std::vector<double>>>> nlm_tot;

            //outermost loop : all adjacent atoms
            nlm_tot.resize(GlobalC::GridD.getAdjacentNum()+1);

            for (int ad=0; ad<GlobalC::GridD.getAdjacentNum()+1 ; ++ad)
            {
                const int T1 = GlobalC::GridD.getType(ad);
                const int I1 = GlobalC::GridD.getNatom(ad);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad);
                const Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int nw1_tot = atom1->nw*GlobalV::NPOL;

                //middle loop : atomic basis on current processor (either row or column)
                nlm_tot[ad].clear();

                const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;

                if (dist1 > Rcut_Alpha + Rcut_AO1)
                {
                    continue;
                }

                for (int iw1=0; iw1<nw1_tot; ++iw1)
                {                   
                    const int iw1_all = start1 + iw1;
                    const int iw1_local = GlobalC::ParaO.trace_loc_row[iw1_all];
                    const int iw2_local = GlobalC::ParaO.trace_loc_col[iw1_all];
                    if(iw1_local < 0 && iw2_local < 0)continue;
                    const int iw1_0 = iw1/GlobalV::NPOL;
                    std::vector<std::vector<double>> nlm;
                    //2D, but first dimension is only 1 here
                    //for force, the right hand side is the gradient
                    //and the first dimension is then 3
                    //inner loop : all projectors (L0,M0)
					GlobalC::UOT.snap_psialpha_half(
						nlm, 1, tau1, T1,
						atom1->iw2l[ iw1_0 ], // L1
						atom1->iw2m[ iw1_0 ], // m1
						atom1->iw2n[ iw1_0 ], // N1
						GlobalC::ucell.atoms[T0].tau[I0], T0, I0, this->inl_index); //R0,T0

                    nlm_tot[ad].insert({iw1_all,nlm});
                }//end iw
            }//end ad

            //=======================================================
            //Step2:	
            //calculate sum_(L0,M0) alpha<psi_i|alpha><alpha|psi_j>
            //and accumulate the value to Hloc_fixed(i,j)
            //=======================================================

            for (int ad1=0; ad1<GlobalC::GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GlobalC::GridD.getType(ad1);
                const int I1 = GlobalC::GridD.getNatom(ad1);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad1);
                const Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int nw1_tot = atom1->nw*GlobalV::NPOL;
                const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

                ModuleBase::Vector3<double> dR1(GlobalC::GridD.getBox(ad1).x, GlobalC::GridD.getBox(ad1).y, GlobalC::GridD.getBox(ad1).z);

                for (int ad2=0; ad2 < GlobalC::GridD.getAdjacentNum()+1 ; ad2++)
                {
                    const int T2 = GlobalC::GridD.getType(ad2);
                    const int I2 = GlobalC::GridD.getNatom(ad2);
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                    const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau(ad2);
                    const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                    const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    ModuleBase::Vector3<double> dR2(GlobalC::GridD.getBox(ad2).x, GlobalC::GridD.getBox(ad2).y, GlobalC::GridD.getBox(ad2).z);
                    
                    const double Rcut_AO2 = GlobalC::ORB.Phi[T2].getRcut();
                    const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                    const double dist2 = (tau2-tau0).norm() * GlobalC::ucell.lat0;

                    if (dist1 > Rcut_Alpha + Rcut_AO1
                            || dist2 > Rcut_Alpha + Rcut_AO2)
                    {
                        continue;
                    }

                    double r0[3];
                    double r1[3];
                    if(isstress)
                    {
                        r1[0] = ( tau1.x - tau0.x) ;
                        r1[1] = ( tau1.y - tau0.y) ;
                        r1[2] = ( tau1.z - tau0.z) ;
                        r0[0] = ( tau2.x - tau0.x) ;
                        r0[1] = ( tau2.y - tau0.y) ;
                        r0[2] = ( tau2.z - tau0.z) ;
                    }

                    for (int iw1=0; iw1<nw1_tot; ++iw1)
                    {
                        const int iw1_all = start1 + iw1;
                        const int iw1_local = GlobalC::ParaO.trace_loc_col[iw1_all];
                        if(iw1_local < 0)continue;

                        for (int iw2=0; iw2<nw2_tot; ++iw2)
                        {
                            const int iw2_all = start2 + iw2;
                            const int iw2_local = GlobalC::ParaO.trace_loc_row[iw2_all];
                            if(iw2_local < 0)continue;

                            double nlm[3]={0,0,0};
                            double nlm_t[3] = {0,0,0}; //for stress
                            std::vector<double> nlm1 = nlm_tot[ad1][iw1_all][0];
                            std::vector<std::vector<double>> nlm2;
                            nlm2.resize(3);

                            for(int dim=0;dim<3;dim++)
                            {
                                nlm2[dim] = nlm_tot[ad2][iw2_all][dim+1];
                            }

                            assert(nlm1.size()==nlm2[0].size());

                            int ib=0;
                            for (int L0 = 0; L0 <= GlobalC::ORB.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < GlobalC::ORB.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < nm; ++m1)
                                    {
                                        for (int m2 = 0; m2 < nm; ++m2)
                                        {
                                            for(int dim=0;dim<3;dim++)
                                            {                                            
                                                nlm[dim] += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[dim][ib+m2];
                                            }
                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());

                            // HF term is minus, only one projector for each atom force.
                            for(int ik=0;ik<GlobalC::kv.nks;ik++)
                            {
                                const double arg = ( GlobalC::kv.kvec_d[ik] * (dR2-dR1) ) * ModuleBase::TWO_PI;
                                const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );
                                F_delta_k(iat, 0) -= 2.0 * dm[ik](iw1_local, iw2_local) * nlm[0] * kphase;
                                F_delta_k(iat, 1) -= 2.0 * dm[ik](iw1_local, iw2_local) * nlm[1] * kphase;
                                F_delta_k(iat, 2) -= 2.0 * dm[ik](iw1_local, iw2_local) * nlm[2] * kphase;
                            }

                            if(isstress)
                            {
                                nlm1 = nlm_tot[ad2][iw2_all][0];
                                for(int i=0;i<3;i++)
                                {
                                    nlm2[i] = nlm_tot[ad1][iw1_all][i+1];
                                }

                                assert(nlm1.size()==nlm2[0].size());                                

                                int ib=0;
                                for (int L0 = 0; L0 <= GlobalC::ORB.Alpha[0].getLmax();++L0)
                                {
                                    for (int N0 = 0;N0 < GlobalC::ORB.Alpha[0].getNchi(L0);++N0)
                                    {
                                        const int inl = this->inl_index[T0](I0, L0, N0);
                                        const int nm = 2*L0+1;
                                        for (int m1 = 0;m1 < nm; ++m1)
                                        {
                                            for (int m2 = 0; m2 < nm; ++m2)
                                            {
                                                for(int dim=0;dim<3;dim++)
                                                {                                            
                                                    nlm_t[dim] += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[dim][ib+m2];
                                                }
                                            }
                                        }
                                        ib+=nm;
                                    }
                                }
                                assert(ib==nlm1.size());
                                for(int ik=0;ik<GlobalC::kv.nks;ik++)
                                {
                                    const double arg = ( GlobalC::kv.kvec_d[ik] * (dR2-dR1) ) * ModuleBase::TWO_PI;
                                    const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );                                    
                                    for(int ipol=0;ipol<3;ipol++)
                                    {
                                        svnl_dalpha_k(0,ipol) -= dm[ik](iw1_local, iw2_local) * kphase * (nlm[0] * r0[ipol] + nlm_t[0] * r1[ipol])* -1;
                                        svnl_dalpha_k(1,ipol) -= dm[ik](iw1_local, iw2_local) * kphase * (nlm[1] * r0[ipol] + nlm_t[1] * r1[ipol])* -1;
                                        svnl_dalpha_k(2,ipol) -= dm[ik](iw1_local, iw2_local) * kphase * (nlm[2] * r0[ipol] + nlm_t[2] * r1[ipol])* -1;
                                    }
                                }
                            }
                        }//iw2
                    }//iw1
                }//ad2
            }//ad1
        }//end I0
    }//end T0

    if(!F_delta_k.checkreal())
    {
        ModuleBase::WARNING_QUIT("cal_f_delta_hf_k","Force should be real!");
    }
    this->F_delta = F_delta_k.dble();

    if(isstress)
    {
        if(!svnl_dalpha_k.checkreal())
        {
            ModuleBase::WARNING_QUIT("cal_f_delta_hf_k","Stress should be real!");
        }

        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                svnl_dalpha(i,j) =  svnl_dalpha_k(i,j).real() * GlobalC::ucell.lat0 / GlobalC::ucell.omega;
            }
        }
    }
}

//calculating sum of correction band energies
//for gamma_only calculations
void LCAO_Descriptor::cal_e_delta_band(const std::vector<ModuleBase::matrix> &dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_e_delta_band");
    this->e_delta_band = 0;
    for (int i = 0; i < GlobalV::NLOCAL; ++i)
    {
        for (int j = 0; j < GlobalV::NLOCAL; ++j)
        {
            const int mu = GlobalC::ParaO.trace_loc_row[j];
            const int nu = GlobalC::ParaO.trace_loc_col[i];
            
            if (mu >= 0 && nu >= 0)
            {                
                const int index=nu*GlobalC::ParaO.nrow+mu;
                for (int is = 0; is < GlobalV::NSPIN; ++is)
                {
                    this->e_delta_band += dm[is](nu, mu) * this->H_V_delta[index];
                }
            }
        }
    }
    Parallel_Reduce::reduce_double_all(this->e_delta_band);
    return;
}

//calculating sum of correction band energies
//for multi_k calculations
void LCAO_Descriptor::cal_e_delta_band_k(const std::vector<ModuleBase::ComplexMatrix> &dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_e_delta_band");
    std::complex<double> e_delta_band_k=std::complex<double>(0.0,0.0);
    for (int i = 0; i < GlobalV::NLOCAL; ++i)
    {
        for (int j = 0; j < GlobalV::NLOCAL; ++j)
        {
            const int mu = GlobalC::ParaO.trace_loc_row[j];
            const int nu = GlobalC::ParaO.trace_loc_col[i];
            
            if (mu >= 0 && nu >= 0)
            {                
                int iic;
                if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                {
                    iic=mu+nu*GlobalC::ParaO.nrow;
                }
                else
                {
                    iic=mu*GlobalC::ParaO.ncol+nu;
                }
                for(int ik=0;ik<GlobalC::kv.nks;ik++)
                {
                    e_delta_band_k += dm[ik](nu, mu) * this->H_V_delta_k[ik][iic];
                }
            }
        }
    }
    Parallel_Reduce::reduce_complex_double_all(e_delta_band_k);
    if(e_delta_band_k.imag()>1e-12)
    {
        GlobalV::ofs_running << "e_delta_band_k : " << e_delta_band_k << std::endl;
        ModuleBase::WARNING_QUIT("e_delta_band_k","energy should be real!");
    }
    this->e_delta_band = e_delta_band_k.real();
    return;
}

//calculates sum_(L0,M0) alpha<psi_i|alpha><alpha|psi_j>
//and accumulate the value to H_V_delta(i,j)
//as well as the counterpart in forces, if calc_deri=1
void LCAO_Descriptor::build_v_delta_alpha_new(const bool& calc_deri)
{
    ModuleBase::TITLE("LCAO_Descriptor", "build_v_delta_alpha_new");
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        ModuleBase::GlobalFunc::ZEROS(this->H_V_delta,GlobalC::ParaO.nloc); //init before calculate
    }

    const double Rcut_Alpha = GlobalC::ORB.Alpha[0].getRcut();
    //same for all types of atoms
    int job;
    if(!calc_deri)
    {
        job=0;
    }
    else
    {
        job=1;
    }

    for (int T0 = 0; T0 < GlobalC::ucell.ntype; T0++)
    {
		Atom* atom0 = &GlobalC::ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = GlobalC::ucell.itia2iat(T0,I0);
			//=======================================================
			//Step1:	
			//saves <beta|psi>, where beta runs over L0,M0 on atom I0
			//and psi runs over atomic basis sets on the current core
			//=======================================================
			std::vector<std::unordered_map<int,std::vector<std::vector<double>>>> nlm_tot;

            //GlobalC::GridD.Find_atom( atom0->tau[I0] );
			const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom0->tau[I0] ,T0, I0);

			//outermost loop : all adjacent atoms
            if(GlobalV::GAMMA_ONLY_LOCAL)
            {
			    nlm_tot.resize(GlobalC::GridD.getAdjacentNum()+1);
            }
            else
            {
                this->nlm_k[iat].resize(GlobalC::GridD.getAdjacentNum()+1);
            }

            for (int ad=0; ad<GlobalC::GridD.getAdjacentNum()+1 ; ++ad)
            {
                const int T1 = GlobalC::GridD.getType(ad);
                const int I1 = GlobalC::GridD.getNatom(ad);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
				const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad);
				const Atom* atom1 = &GlobalC::ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;

				//middle loop : atomic basis on current processor (either row or column)
                if(GlobalV::GAMMA_ONLY_LOCAL)
				{
                    nlm_tot[ad].clear();
                }
                else
                {
                    this->nlm_k[iat][ad].clear();
                }
				const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;

				if (dist1 > Rcut_Alpha + Rcut_AO1)
				{
					continue;
				}

				for (int iw1=0; iw1<nw1_tot; ++iw1)
				{
					const int iw1_all = start1 + iw1;
					const int iw1_local = GlobalC::ParaO.trace_loc_row[iw1_all];
					const int iw2_local = GlobalC::ParaO.trace_loc_col[iw1_all];
					if(iw1_local < 0 && iw2_local < 0)continue;
					const int iw1_0 = iw1/GlobalV::NPOL;
					std::vector<std::vector<double>> nlm;
					//2D, but first dimension is only 1 here
					//for force, the right hand side is the gradient
					//and the first dimension is then 3
					//inner loop : all projectors (L0,M0)
					GlobalC::UOT.snap_psialpha_half(
						nlm, job, tau1, T1,
						atom1->iw2l[ iw1_0 ], // L1
						atom1->iw2m[ iw1_0 ], // m1
						atom1->iw2n[ iw1_0 ], // N1
						GlobalC::ucell.atoms[T0].tau[I0], T0, I0, this->inl_index); //R0,T0
                    if(GlobalV::GAMMA_ONLY_LOCAL)
                    {
					    nlm_tot[ad].insert({iw1_all,nlm});
                    }
                    else
                    {
                        this->nlm_k[iat][ad].insert({iw1_all,nlm});
                    }
				}//end iw
			}//end ad

            if(!GlobalV::GAMMA_ONLY_LOCAL && !calc_deri)
            {
                //saves <alpha(0)|psi(R)>, to be used later 
                continue;
            }
			//=======================================================
			//Step2:	
			//calculate sum_(L0,M0) alpha<psi_i|alpha><alpha|psi_j>
			//and accumulate the value to H_V_delta(i,j)
            //as well as the counterpart in forces, if required
			//=======================================================

			for (int ad1=0; ad1<GlobalC::GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GlobalC::GridD.getType(ad1);
                const int I1 = GlobalC::GridD.getNatom(ad1);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &GlobalC::ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

                ModuleBase::Vector3<double> dR1(GlobalC::GridD.getBox(ad1).x, GlobalC::GridD.getBox(ad1).y, GlobalC::GridD.getBox(ad1).z);

				for (int ad2=0; ad2 < GlobalC::GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GlobalC::GridD.getType(ad2);
					const int I2 = GlobalC::GridD.getNatom(ad2);
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &GlobalC::ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
					
					const double Rcut_AO2 = GlobalC::ORB.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * GlobalC::ucell.lat0;

                    ModuleBase::Vector3<double> dR2(GlobalC::GridD.getBox(ad2).x, GlobalC::GridD.getBox(ad2).y, GlobalC::GridD.getBox(ad2).z);

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}					

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1;
						const int iw1_local = GlobalC::ParaO.trace_loc_row[iw1_all];
						if(iw1_local < 0)continue;
						const int iw1_0 = iw1/GlobalV::NPOL;
						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2;
							const int iw2_local = GlobalC::ParaO.trace_loc_col[iw2_all];
							if(iw2_local < 0)continue;
							const int iw2_0 = iw2/GlobalV::NPOL;

                            if(calc_deri)
                            {
                                double nlm[3]={0,0,0};
                                std::vector<double> nlm1;
                                std::vector<std::vector<double>> nlm2;
                                if(GlobalV::GAMMA_ONLY_LOCAL)
                                {
                                    nlm1 = nlm_tot[ad1][iw1_all][0];
                                }
                                else
                                {
                                    nlm1 = nlm_k[iat][ad1][iw1_all][0];
                                }
                                
                                nlm2.resize(3);
                                for(int dim=0;dim<3;dim++)
                                {
                                    if(GlobalV::GAMMA_ONLY_LOCAL)
                                    {
                                        nlm2[dim] = nlm_tot[ad2][iw2_all][dim+1];
                                    }
                                    else
                                    {
                                        nlm2[dim] = nlm_k[iat][ad2][iw2_all][dim+1];
                                    }
                                }

                                assert(nlm1.size()==nlm2[0].size());
                                int ib=0;
                                for (int L0 = 0; L0 <= GlobalC::ORB.Alpha[0].getLmax();++L0)
                                {
                                    for (int N0 = 0;N0 < GlobalC::ORB.Alpha[0].getNchi(L0);++N0)
                                    {
                                        const int inl = this->inl_index[T0](I0, L0, N0);
                                        const int nm = 2*L0+1;
                                        for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                        {
                                            for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                            {
                                                for(int dim=0;dim<3;dim++)
                                                {
                                                    nlm[dim] += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[dim][ib+m2];
                                                }
                                            }
                                        }
                                        ib+=(2*L0+1);
                                    }
                                }

                                assert(ib==nlm1.size());

                                int index = iw1_local * GlobalC::ParaO.ncol+ iw2_local;
                                if(GlobalV::GAMMA_ONLY_LOCAL)
                                {
                                    this->DH_V_delta_x[index] += nlm[0];
                                    this->DH_V_delta_y[index] += nlm[1];
                                    this->DH_V_delta_z[index] += nlm[2];
                                }
                                else
                                {
                                    for(int ik=0;ik<GlobalC::kv.nks;ik++)
                                    {
                                        const double arg = ( GlobalC::kv.kvec_d[ik] * (dR2-dR1) ) * ModuleBase::TWO_PI;
                                        const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );
                                        this->DH_V_delta_x_k[ik][index] += nlm[0] * kphase;
                                        this->DH_V_delta_y_k[ik][index] += nlm[1] * kphase;
                                        this->DH_V_delta_z_k[ik][index] += nlm[2] * kphase;
                                    }
                                }

                            }
                            else
                            {
                                std::vector<double> nlm1 = nlm_tot[ad1][iw1_all][0];
                                std::vector<double> nlm2 = nlm_tot[ad2][iw2_all][0];

                                assert(nlm1.size()==nlm2.size());
							    double nlm=0.0;
							    int ib = 0;

                                for (int L0 = 0; L0 <= GlobalC::ORB.Alpha[0].getLmax();++L0)
                                {
                                    for (int N0 = 0;N0 < GlobalC::ORB.Alpha[0].getNchi(L0);++N0)
                                    {
                                        const int inl = this->inl_index[T0](I0, L0, N0);
                                        const int nm = 2*L0+1;
                                        for (int m01 = 0;m01 < nm;++m01)
                                        {
                                            for (int m02 = 0; m02 < nm; ++m02)
                                            {
                                                nlm += this->gedm[inl][m01*nm+m02]*nlm1[ib+m01]*nlm2[ib+m02];
                                            }
                                        }
                                        ib+=(2*L0+1);
                                    }
                                }
                                int index = iw2_local * GlobalC::ParaO.nrow+ iw1_local;     //for genelpa
                                this->H_V_delta[index] += nlm;

                                assert(ib==nlm1.size());
                            }
						}//iw2
					}//iw1
				}//ad2
			}//ad1
		}//end I0
	}//end T0
/*
    for(int iw1=0;iw1<GlobalC::ParaO.ncol;iw1++)
    {
        for(int iw2=0;iw2<GlobalC::ParaO.nrow;iw2++)
        {
            GlobalV::ofs_running << H_V_delta[iw1*GlobalC::ParaO.nrow+iw2] << " ";
        }
        GlobalV::ofs_running << std::endl;
    }
*/
    ModuleBase::timer::tick ("LCAO_gen_fixedH","build_Nonlocal_alpha_new");
	return;

}

#endif