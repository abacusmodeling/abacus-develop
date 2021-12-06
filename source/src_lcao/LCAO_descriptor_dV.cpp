//wenfei add 2021 october
#ifdef __DEEPKS

#include "LCAO_nnr.h"
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
        H_V_deltaR = new double[GlobalC::LNNR.nnr];
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
    }
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        this->nlm_save.resize(GlobalC::ucell.nat);
    }
    else
    {
        this->nlm_save_k.resize(GlobalC::ucell.nat);
    }
    return;
}

void LCAO_Descriptor::resize_nlm()
{
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        this->nlm_save.resize(GlobalC::ucell.nat);
    }
    else
    {
        this->nlm_save_k.resize(GlobalC::ucell.nat);
    }
    return;
}

//this subroutine adds dV to the Kohn-Sham Hamiltonian
//for gamma_only calculations
void LCAO_Descriptor::add_v_delta(void)
{
    ModuleBase::TITLE("LCAO_DESCRIPTOR", "add_v_delta");
    ModuleBase::timer::tick ("LCAO_gen_fixedH","add_v_delta");
    ModuleBase::GlobalFunc::ZEROS(this->H_V_delta,GlobalC::ParaO.nloc); //init before calculate

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
                            std::vector<double> nlm1 = this->nlm_save[iat][ad1][iw1_all][0];
                            std::vector<double> nlm2 = this->nlm_save[iat][ad2][iw2_all][0];

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
                                    ib+=nm;
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
                            this->H_V_delta[iic] += nlm;
                            GlobalC::LM.Hloc[iic] += nlm;
						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }
    }

    ModuleBase::timer::tick ("LCAO_DESCRIPTOR","add_v_delta");
	return;
}

//this subroutine calculates H_V_deltaR
//used in multi-k calculations
void LCAO_Descriptor::add_v_delta_k(void)
{
    ModuleBase::TITLE("LCAO_DESCRIPTOR", "add_v_delta_k");
    ModuleBase::timer::tick ("LCAO_DESCRIPTOR","add_v_delta_k");
    ModuleBase::GlobalFunc::ZEROS(this->H_V_deltaR, GlobalC::LNNR.nnr);

    const double Rcut_Alpha = GlobalC::ORB.Alpha[0].getRcut();

	int nnr = 0;
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
	ModuleBase::Vector3<float> dtau1_f, dtau2_f;
	double distance = 0.0;
	double rcut = 0.0;
	double rcut1, rcut2;
		
	//	Record_adj RA;
	//	RA.for_2d();

	// psi1
	for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
	{
		const Atom* atom1 = &GlobalC::ucell.atoms[T1];
		for (int I1 =0; I1< atom1->na; ++I1)
		{

			GlobalC::GridD.Find_atom(GlobalC::ucell, atom1->tau[I1] ,T1, I1);
			const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			tau1 = atom1->tau[I1];

			// psi2
			for (int ad2=0; ad2<GlobalC::GridD.getAdjacentNum()+1; ++ad2)
			{
				const int T2 = GlobalC::GridD.getType(ad2);
				const Atom* atom2 = &GlobalC::ucell.atoms[T2];
				
				const int I2 = GlobalC::GridD.getNatom(ad2);
				const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
				const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
				tau2 = GlobalC::GridD.getAdjacentTau(ad2);

				bool is_adj = false;

				const int rx2=GlobalC::GridD.getBox(ad2).x;
				const int ry2=GlobalC::GridD.getBox(ad2).y;
				const int rz2=GlobalC::GridD.getBox(ad2).z;

					
				dtau = tau2 - tau1;
				distance = dtau.norm2() * pow(GlobalC::ucell.lat0,2);
				// this rcut is in order to make nnr consistent 
				// with other matrix.
				rcut = pow(GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut(),2);
				if(distance < rcut) is_adj = true;
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GlobalC::GridD.getType(ad0);
						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int T0 = RA.info[iat1][ad0][3];
						//const int I0 = RA.info[iat1][ad0][4];
						//const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GlobalC::GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;

						const double distance1 = dtau1.norm2() * pow(GlobalC::ucell.lat0,2);
						const double distance2 = dtau2.norm2() * pow(GlobalC::ucell.lat0,2);

						rcut1 = pow(GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(),2);
						rcut2 = pow(GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(),2);

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							is_adj = true;
							break;
						}
					}
				}


				if(is_adj)
				{
					// < psi1 | all projectors | psi2 >
					// ----------------------------- enter the nnr increaing zone -------------------------
					for (int ad0=0; ad0 < GlobalC::GridD.getAdjacentNum()+1 ; ++ad0)
					{
						const int T0 = GlobalC::GridD.getType(ad0);
						const int I0 = GlobalC::GridD.getNatom(ad0);
						const int iat = GlobalC::ucell.itia2iat(T0,I0);

						// mohan add 2010-12-19
						if( GlobalC::ucell.infoNL.nproj[T0] == 0) continue;

						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
						tau0 = GlobalC::GridD.getAdjacentTau(ad0);

						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;
						const double distance1 = dtau1.norm2() * pow(GlobalC::ucell.lat0,2);
						const double distance2 = dtau2.norm2() * pow(GlobalC::ucell.lat0,2);

						// seems a bug here!! mohan 2011-06-17
						rcut1 = pow(GlobalC::ORB.Phi[T1].getRcut() + Rcut_Alpha,2);
						rcut2 = pow(GlobalC::ORB.Phi[T2].getRcut() + Rcut_Alpha,2);

						if(distance1 >= rcut1 || distance2 >= rcut2)
						{
							continue;
						}
						//const Atom* atom0 = &GlobalC::ucell.atoms[T0];
						const int rx0=GlobalC::GridD.getBox(ad0).x;
						const int ry0=GlobalC::GridD.getBox(ad0).y;
						const int rz0=GlobalC::GridD.getBox(ad0).z;
						key_tuple key_1(iat1,-rx0,-ry0,-rz0);
						key_tuple key_2(iat2,rx2-rx0,ry2-ry0,rz2-rz0);

                        int nnr_inner = 0;

						for (int iw1=0; iw1<atom1->nw*GlobalV::NPOL; iw1++)
						{
							const int iw1_all = start1 + iw1;
							const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
							if(mu < 0)continue; 

							// fix a serious bug: atom2[T2] -> atom2
							// mohan 2010-12-20
							for (int iw2=0; iw2<atom2->nw*GlobalV::NPOL; iw2++)
							{
								const int iw2_all = start2 + iw2;
								const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];						
								if(nu < 0)continue;
  
                                std::vector<double> nlm1 = this->nlm_save_k[iat][key_1][iw1_all][0];
                                std::vector<double> nlm2 = this->nlm_save_k[iat][key_2][iw2_all][0];
                                assert(nlm1.size()==nlm2.size());

                                int ib=0;
                                double nlm = 0.0;
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

                                this->H_V_deltaR[nnr+nnr_inner] += nlm;
                                nnr_inner++;
                            }//iw2
                        }//iw1
                    }//ad0

					//outer circle : accumulate nnr
					for (int j=0; j<atom1->nw*GlobalV::NPOL; j++)
					{
						const int j0 = j/GlobalV::NPOL;//added by zhengdy-soc
						const int iw1_all = start1 + j;
						const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
						if(mu < 0)continue; 

						// fix a serious bug: atom2[T2] -> atom2
						// mohan 2010-12-20
						for (int k=0; k<atom2->nw*GlobalV::NPOL; k++)
						{
							const int k0 = k/GlobalV::NPOL;
							const int iw2_all = start2 + k;
							const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];						
							if(nu < 0)continue;

							nnr++;
						}
					}
                }//is_adj
            }//ad2
        }//I1
    }//T1

    if( nnr!=GlobalC::LNNR.nnr)
    {
        ModuleBase::WARNING_QUIT("LCAO_gen_fixedH::build_Nonlocal_mu_new","nnr!=GlobalC::LNNR.nnr");
    }

    ModuleBase::timer::tick ("LCAO_DESCRIPTOR","add_v_delta_k");
    return;
}


//force for gamma only calculations
//Pulay and HF terms are calculated together
void LCAO_Descriptor::cal_f_delta_new(const ModuleBase::matrix& dm, const bool isstress, ModuleBase::matrix& svnl_dalpha)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_f_delta_new");
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
                    const int ibt = GlobalC::ucell.itia2iat(T2,I2);
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
                            std::vector<double> nlm1 = this->nlm_save[iat][ad1][iw1_all][0];
                            std::vector<std::vector<double>> nlm2;
                            nlm2.resize(3);

                            for(int dim=0;dim<3;dim++)
                            {
                                nlm2[dim] = this->nlm_save[iat][ad2][iw2_all][dim+1];
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

                            // Pulay term is plus, only one projector for each atom force.
                            this->F_delta(ibt, 0) += 2 * dm(iw1_local, iw2_local) * nlm[0];
                            this->F_delta(ibt, 1) += 2 * dm(iw1_local, iw2_local) * nlm[1];
                            this->F_delta(ibt, 2) += 2 * dm(iw1_local, iw2_local) * nlm[2];

                            if(isstress)
                            {
                                nlm1 = this->nlm_save[iat][ad2][iw2_all][0];
                                for(int i=0;i<3;i++)
                                {
                                    nlm2[i] = this->nlm_save[iat][ad1][iw1_all][i+1];
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
                                    svnl_dalpha(0,ipol) -= dm(iw1_local, iw2_local) * (nlm[0] * r0[ipol] + nlm_t[0] * r1[ipol])* -1.0;
                                    svnl_dalpha(1,ipol) -= dm(iw1_local, iw2_local) * (nlm[1] * r0[ipol] + nlm_t[1] * r1[ipol])* -1.0;
                                    svnl_dalpha(2,ipol) -= dm(iw1_local, iw2_local) * (nlm[2] * r0[ipol] + nlm_t[2] * r1[ipol])* -1.0;
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

    return;
}

//force for multi-k calculations
//Pulay and HF terms are calculated together

void LCAO_Descriptor::cal_f_delta_k(const std::vector<ModuleBase::ComplexMatrix>& dm/**<[in] density matrix*/, const bool isstress, ModuleBase::matrix& svnl_dalpha)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_f_delta_hf_k_new");
    ModuleBase::timer::tick("LCAO_Descriptor","cal_f_delta_hf_k_new");
    this->F_delta.zero_out();

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
                const int ibt1 = GlobalC::ucell.itia2iat(T1,I1);
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
                    const int ibt2 = GlobalC::ucell.itia2iat(T2,I2);
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
                            double dm_current;
                            std::complex<double> tmp = 0.0;
                            for(int ik=0;ik<GlobalC::kv.nks;ik++)
                            {
                                const double arg = ( GlobalC::kv.kvec_d[ik] * (dR2-dR1) ) * ModuleBase::TWO_PI;
                                const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );
                                tmp += dm[ik](iw1_local, iw2_local) * kphase;
                            }
                            if(tmp.imag()>1.0e-8)
                            {
                                GlobalV::ofs_running << "dm_current not real in force_hf" << tmp << "\n";
                            }
                            dm_current=tmp.real();


                            double nlm[3]={0,0,0};
                            double nlm_t[3] = {0,0,0}; //for stress
                            key_tuple key_1(ibt1,dR1.x,dR1.y,dR1.z);
                            key_tuple key_2(ibt2,dR2.x,dR2.y,dR2.z);
                            std::vector<double> nlm1 = this->nlm_save_k[iat][key_1][iw1_all][0];
                            std::vector<std::vector<double>> nlm2;
                            nlm2.resize(3);
                            for(int dim=0;dim<3;dim++)
                            {
                                nlm2[dim] = this->nlm_save_k[iat][key_2][iw2_all][dim+1];
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

                            // Pulay term is plus
                            this->F_delta(ibt2, 0) += 2.0 * dm_current * nlm[0];
                            this->F_delta(ibt2, 1) += 2.0 * dm_current * nlm[1];
                            this->F_delta(ibt2, 2) += 2.0 * dm_current * nlm[2];

                            // HF term is minus, only one projector for each atom force.
                            this->F_delta(iat, 0) -= 2.0 * dm_current * nlm[0];
                            this->F_delta(iat, 1) -= 2.0 * dm_current * nlm[1];
                            this->F_delta(iat, 2) -= 2.0 * dm_current * nlm[2];

                            if(isstress)
                            {
                                nlm1 = this->nlm_save_k[iat][key_2][iw2_all][0];
                                for(int i=0;i<3;i++)
                                {
                                    nlm2[i] = this->nlm_save_k[iat][key_1][iw1_all][i+1];
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
                                    svnl_dalpha(0,ipol) -= dm_current * (nlm[0] * r0[ipol] + nlm_t[0] * r1[ipol])* -1.0;
                                    svnl_dalpha(1,ipol) -= dm_current * (nlm[1] * r0[ipol] + nlm_t[1] * r1[ipol])* -1.0;
                                    svnl_dalpha(2,ipol) -= dm_current * (nlm[2] * r0[ipol] + nlm_t[2] * r1[ipol])* -1.0;
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
                svnl_dalpha(i,j) =  svnl_dalpha(i,j) * GlobalC::ucell.lat0 / GlobalC::ucell.omega;
            }
        }
    }

    ModuleBase::timer::tick("LCAO_Descriptor","cal_f_delta_hf_k_new");
    return;
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
        //ModuleBase::WARNING_QUIT("e_delta_band_k","energy should be real!");
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
    ModuleBase::timer::tick ("LCAO_Descriptor","build_v_delta_alpha_new");

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
            //Step 1 : 
			//saves <alpha|psi>, where alpha runs over all projectors
			//and psi runs over atomic basis sets on the current core
			//=======================================================

			const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom0->tau[I0] ,T0, I0);

			//outermost loop : all adjacent atoms
            if(GlobalV::GAMMA_ONLY_LOCAL)
            {
                this->nlm_save[iat].resize(GlobalC::GridD.getAdjacentNum()+1);
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

                std::unordered_map<int,std::vector<std::vector<double>>> nlm_cur;
                if(GlobalV::GAMMA_ONLY_LOCAL)
                {
                    this->nlm_save[iat][ad].clear();
                }
                else
                {
                    nlm_cur.clear();
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
						GlobalC::ucell.atoms[T0].tau[I0], T0, I0); //R0,T0

                    if(GlobalV::GAMMA_ONLY_LOCAL)
                    {
                        this->nlm_save[iat][ad].insert({iw1_all,nlm});
                    }
                    else
                    {
                        nlm_cur.insert({iw1_all,nlm});
                    }
				}//end iw

                if(!GlobalV::GAMMA_ONLY_LOCAL)
                {
                    const int ibt=GlobalC::ucell.itia2iat(T1, I1);
                    const int rx=GlobalC::GridD.getBox(ad).x;
                    const int ry=GlobalC::GridD.getBox(ad).y;
                    const int rz=GlobalC::GridD.getBox(ad).z;
                    key_tuple key_1(ibt,rx,ry,rz);
                    this->nlm_save_k[iat][key_1]=nlm_cur;
                }
			}//end ad
		}//end I0
	}//end T0

    ModuleBase::timer::tick ("LCAO_Descriptor","build_v_delta_alpha_new");
	return;

}

#endif