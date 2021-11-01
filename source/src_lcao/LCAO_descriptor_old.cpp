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
#endif