//keep some no longer used deepks subroutines here
#ifdef __DEEPKS

#include "LCAO_descriptor.h"

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

//Note: with the modified gdmx, we can not longer use this subroutine
//for calculating force!!!
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

//for multi-k, search adjacent atoms from mu
void LCAO_Descriptor::build_v_delta_mu(const bool& calc_deri)
{
    ModuleBase::TITLE("LCAO_Descriptor", "build_v_delta_mu");
    ModuleBase::GlobalFunc::ZEROS(this->H_V_delta, GlobalV::NLOCAL * GlobalV::NLOCAL); //init before calculate
    //timer::tick ("LCAO_gen_fixedH","build_Nonlocal_mu");

    // < phi1 | beta > < beta | phi2 >
	// phi1 is within the unitcell.
	// while beta is in the supercell.
	// while phi2 is in the supercell.

	int nnr = 0;
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
	double distance = 0.0;
	double distance1, distance2;
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
            //GlobalC::GridD.Find_atom( atom1->tau[I1] );
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom1->tau[I1] ,T1, I1);
			//const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            tau1 = atom1->tau[I1];

			// psi2
            for (int ad2=0; ad2<GlobalC::GridD.getAdjacentNum()+1; ++ad2)
			{
				const int T2 = GlobalC::GridD.getType(ad2);
				const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                
				const int I2 = GlobalC::GridD.getNatom(ad2);
				//const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                tau2 = GlobalC::GridD.getAdjacentTau(ad2);

				bool is_adj = false;
					
				dtau = tau2 - tau1;
				distance = dtau.norm() * GlobalC::ucell.lat0;
				// this rcut is in order to make nnr consistent 
				// with other matrix.
				rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
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

                        double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

                        rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Alpha[0].getRcut();
                        rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ORB.Alpha[0].getRcut();

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


							//(3) run over all projectors in nonlocal pseudopotential.
							for (int ad0=0; ad0 < GlobalC::GridD.getAdjacentNum()+1 ; ++ad0)
							{
								const int T0 = GlobalC::GridD.getType(ad0);
                                const int I0 = GlobalC::GridD.getNatom(ad0);
								tau0 = GlobalC::GridD.getAdjacentTau(ad0);

								dtau1 = tau0 - tau1;
								dtau2 = tau0 - tau2;
								distance1 = dtau1.norm() * GlobalC::ucell.lat0;
								distance2 = dtau2.norm() * GlobalC::ucell.lat0;

								// seems a bug here!! mohan 2011-06-17
								rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Alpha[0].getRcut();
								rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ORB.Alpha[0].getRcut();

								if(distance1 < rcut1 && distance2 < rcut2)
								{
									//const Atom* atom0 = &GlobalC::ucell.atoms[T0];
									double nlm[3]={0,0,0};
									if(!calc_deri)
									{
										GlobalC::UOT.snap_psialpha(
												nlm, 0, tau1, T1,
												atom1->iw2l[ j0 ], // L1
												atom1->iw2m[ j0 ], // m1
												atom1->iw2n[ j0 ], // N1
												tau2, T2,
												atom2->iw2l[ k0 ], // L2
												atom2->iw2m[ k0 ], // m2
												atom2->iw2n[ k0 ], // n2
												tau0, T0, I0,
                                                this->inl_index,
                                                this->gedm);


										if(GlobalV::GAMMA_ONLY_LOCAL)
										{
											// mohan add 2010-12-20
											if( nlm[0]!=0.0 )
											{
                                                //GlobalC::LM.set_HSgamma(iw1_all,iw2_all,nlm[0],'N');//N stands for nonlocal.
                                                int index = iw2_all * GlobalV::NLOCAL + iw1_all;     //for genelpa
                                                this->H_V_delta[index] += nlm[0];
                                            }
										}
										else
                                        {
                                            //for multi-k, not prepared yet 
                                        }
									}// calc_deri
									else // calculate the derivative
									{
										if(GlobalV::GAMMA_ONLY_LOCAL)
										{
                                            GlobalC::UOT.snap_psialpha(
                                                    nlm, 1, tau1, T1,
                                                    atom1->iw2l[ j0 ], // L1
                                                    atom1->iw2m[ j0 ], // m1
                                                    atom1->iw2n[ j0 ], // N1
                                                    tau2, T2,
                                                    atom2->iw2l[ k0 ], // L2
                                                    atom2->iw2m[ k0 ], // m2
                                                    atom2->iw2n[ k0 ], // n2
                                                    tau0, T0, I0,
                                                    this->inl_index,
                                                    this->gedm);

											// sum all projectors for one atom.
											//GlobalC::LM.set_force (iw1_all, iw2_all,	nlm[0], nlm[1], nlm[2], 'N');
										}
										else
										{
											// mohan change the order on 2011-06-17
											// origin: < psi1 | beta > < beta | dpsi2/dtau >
											//now: < psi1/dtau | beta > < beta | psi2 >
											GlobalC::UOT.snap_psialpha(
													nlm, 1, tau2, T2,
													atom2->iw2l[ k0 ], // L2
													atom2->iw2m[ k0 ], // m2
													atom2->iw2n[ k0 ], // n2
													tau1, T1,
													atom1->iw2l[ j0 ], // L1
													atom1->iw2m[ j0 ], // m1
													atom1->iw2n[ j0 ], // N1
                                                    tau0, T0, I0,
                                                    this->inl_index,
                                                    this->gedm);

											//GlobalC::LM.DHloc_fixedR_x[nnr] += nlm[0];
											//GlobalC::LM.DHloc_fixedR_y[nnr] += nlm[1];
											//GlobalC::LM.DHloc_fixedR_z[nnr] += nlm[2];
										}
									}//!calc_deri
								}// distance
							} // ad0
							++nnr;
						}// k
					} // j 
				}// end is_adj
			} // ad2
		} // I1
	} // T1

    //timer::tick ("LCAO_gen_fixedH","build_Nonlocal_mu");
	return;
}

//for GAMMA_ONLY, search adjacent atoms from I0
void LCAO_Descriptor::build_v_delta_alpha(const bool& calc_deri)
{
    ModuleBase::TITLE("LCAO_Descriptor", "build_v_delta_alpha");
    ModuleBase::GlobalFunc::ZEROS(this->H_V_delta, GlobalV::NLOCAL * GlobalV::NLOCAL); //init before calculate

    std::vector<double> Rcut;
	for(int it1=0; it1<GlobalC::ucell.ntype; ++it1)
        Rcut.push_back(GlobalC::ORB.Phi[it1].getRcut() + GlobalC::ORB.Alpha[0].getRcut());
    
    for (int T0 = 0;T0 < GlobalC::ucell.ntype;++T0)
    {
        for (int I0 = 0;I0 < GlobalC::ucell.atoms[T0].na;++I0)
        {
            const ModuleBase::Vector3<double> tau0 = GlobalC::ucell.atoms[T0].tau[I0];
            //Rcut in this function may need to be changed ?! (I think the range of adjacent atoms here should be rcut(phi)+rcut(alpha))
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau0, T0, I0);

            //adj atom pairs
            for (int ad1 = 0;ad1 < GlobalC::GridD.getAdjacentNum() + 1;++ad1)
            {
                const int T1 = GlobalC::GridD.getType(ad1);
                const int I1 = GlobalC::GridD.getNatom(ad1);
				//const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad1);
                const Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int nw1_tot = atom1->nw * GlobalV::NPOL;
                for (int ad2 = 0;ad2 < GlobalC::GridD.getAdjacentNum() + 1;++ad2)
                {
                    const int T2 = GlobalC::GridD.getType(ad2);
					const int I2 = GlobalC::GridD.getNatom(ad2);
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                    const int nw2_tot = atom2->nw * GlobalV::NPOL;

                    ModuleBase::Vector3<double> dtau10 = tau1 - tau0;
                    ModuleBase::Vector3<double> dtau20 = tau2 - tau0;
                    double distance10 = dtau10.norm() * GlobalC::ucell.lat0;
                    double distance20 = dtau20.norm() * GlobalC::ucell.lat0;
                    
                    if (distance10 < Rcut[T1] && distance20 < Rcut[T2])
                    {
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

								double nlm[3];
								nlm[0] = nlm[1] = nlm[2] = 0.0;

								if(!calc_deri)
								{
									GlobalC::UOT.snap_psialpha(
											nlm, 0, tau1, T1,
											atom1->iw2l[ iw1_0 ], // L1
											atom1->iw2m[ iw1_0 ], // m1
											atom1->iw2n[ iw1_0 ], // N1
											tau2, T2,
											atom2->iw2l[ iw2_0 ], // L2
											atom2->iw2m[ iw2_0 ], // m2
											atom2->iw2n[ iw2_0 ], // n2
											GlobalC::ucell.atoms[T0].tau[I0], T0, I0, 
                                            this->inl_index,
                                            this->gedm);
                                    int index=iw2_local*GlobalC::ParaO.nrow+iw1_local;
                                    this->H_V_delta[index] += nlm[0];
                                }
								else  // calculate force
								{
									GlobalC::UOT.snap_psialpha(
											nlm, 1, tau1, T1,
											atom1->iw2l[ iw1_0 ], // L1
											atom1->iw2m[ iw1_0 ], // m1
											atom1->iw2n[ iw1_0 ], // N1
											tau2, T2,
											atom2->iw2l[ iw2_0 ], // L2
											atom2->iw2m[ iw2_0 ], // m2
											atom2->iw2n[ iw2_0 ], // n2
											GlobalC::ucell.atoms[T0].tau[I0], T0, I0, 
                                            this->inl_index,
                                            this->gedm);
                                    //for Pulay Force
                                    //GlobalC::LM.set_force(iw1_all, iw2_all, nlm[0], nlm[1], nlm[2], 'N');
                                    const int ir = GlobalC::ParaO.trace_loc_row[ iw1_all ];
                                    const int ic = GlobalC::ParaO.trace_loc_col[ iw2_all ];
                                    const long index = ir * GlobalC::ParaO.ncol + ic;
                                    //this->DH_V_delta_x[index] += nlm[0];
                                    //this->DH_V_delta_y[index] += nlm[1];
                                    //this->DH_V_delta_z[index] += nlm[2];
                                    //not used anymore
                                }
							}// end iw2
						}// end iw1
					} // end distance
                }//end ad2
            }//end ad1
        }//end I0
    }//end T0

    /*
    for(int iw1=0;iw1<GlobalV::NLOCAL;iw1++)
    {
        for(int iw2=0;iw2<GlobalV::NLOCAL;iw2++)
        {
            GlobalV::ofs_running << H_V_delta[iw1*GlobalV::NLOCAL+iw2] << " ";
        }
        GlobalV::ofs_running << std::endl;
    }
    */
    return;
}

//The hellmann-feynmann term in force
void LCAO_Descriptor::cal_f_delta_hf(const ModuleBase::matrix& dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_f_delta_hf");
    this->F_delta.zero_out();
    for (int iat = 0; iat < GlobalC::ucell.nat; ++iat)
    {
        const int it = GlobalC::ucell.iat2it[iat];
        const int ia = GlobalC::ucell.iat2ia[iat];
        const ModuleBase::Vector3<double> tau0 = GlobalC::ucell.atoms[it].tau[ia];
		GlobalC::GridD.Find_atom(GlobalC::ucell, GlobalC::ucell.atoms[it].tau[ia] ,it, ia);
		const double Rcut_Alpha = GlobalC::ORB.Alpha[0].getRcut();

        //FOLLOWING ARE CONTRIBUTIONS FROM
        //VNL DUE TO PROJECTOR'S DISPLACEMENT
        for (int ad1 =0 ; ad1 < GlobalC::GridD.getAdjacentNum()+1; ad1++)
        {
            const int T1 = GlobalC::GridD.getType (ad1);
            const Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int I1 = GlobalC::GridD.getNatom (ad1);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau (ad1);
			const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

            for (int ad2 =0 ; ad2 < GlobalC::GridD.getAdjacentNum()+1; ad2++)
            {
                const int T2 = GlobalC::GridD.getType (ad2);
                const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                const int I2 = GlobalC::GridD.getNatom (ad2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau (ad2);
                const double Rcut_AO2 = GlobalC::ORB.Phi[T2].getRcut();

                const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                const double dist2 = (tau2-tau0).norm() * GlobalC::ucell.lat0;

                if (dist1 > Rcut_Alpha + Rcut_AO1
                        || dist2 > Rcut_Alpha + Rcut_AO2)
                {
                    continue;
                }

                for (int jj = 0; jj < GlobalC::ucell.atoms[T1].nw; jj++)
                {
                    const int iw1_all = start1 + jj;
                    const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
                    if(mu<0) continue;
                    for (int kk = 0; kk < GlobalC::ucell.atoms[T2].nw; kk++)
                    {
                        const int iw2_all = start2 + kk;
                        const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
                        if(nu<0) continue;
                    
                        double nlm[3] = {0,0,0};
                                
                        GlobalC::UOT.snap_psialpha(
                            nlm, 1,
                            tau1, T1,
                            atom1->iw2l[jj], // L2
                            atom1->iw2m[jj], // m2
                            atom1->iw2n[jj], // N2
                            tau2, T2,
                            atom2->iw2l[kk], // L1
                            atom2->iw2m[kk], // m1
                            atom2->iw2n[kk], // n1
                            tau0, it, ia, 
                            this->inl_index,
                            this->gedm); // mohan  add 2021-05-07

                        double nlm1[3] = {0,0,0};

                        const int index = mu * GlobalC::ParaO.ncol + nu;

                        // HF term is minus, only one projector for each atom force.

                        double sum_dm = 0.0;
                        //remaining: sum for is
                        this->F_delta(iat, 0) -= 2 * dm(mu, nu) * nlm[0];
                        this->F_delta(iat, 1) -= 2 * dm(mu, nu) * nlm[1];
                        this->F_delta(iat, 2) -= 2 * dm(mu, nu) * nlm[2];
                        //this->F_delta(iat, 0) -= 4 * dm(mu, nu) * nlm[0];   //2 for v_delta(not calculated togethor), 2 for e_delta
                        //this->F_delta(iat, 1) -= 4 * dm(mu, nu) * nlm[1];
                        //this->F_delta(iat, 2) -= 4 * dm(mu, nu) * nlm[2];
                    }//!kk
                }//!ad2
            }//!jj
        }//!ad1
    }//!iat
    return;
}

#endif