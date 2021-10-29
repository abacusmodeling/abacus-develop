//wenfei add 2021 october
#ifdef __DEEPKS

#include "LCAO_descriptor.h"
#include "LCAO_matrix.h"
#include "../module_base/lapack_connector.h"
#include "../module_base/intarray.h"
#include "../module_base/complexmatrix.h"
#include "global_fp.h"
#include "../src_pw/global.h"
#include "../src_io/winput.h"

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

            if(!GlobalV::GAMMA_ONLY_LOCAL)
            {
                //saves <alpha(0)|psi(R)>, to be used later 
                continue;
            }
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
                                this->DH_V_delta_x[index] += nlm[0];
                                this->DH_V_delta_y[index] += nlm[1];
                                this->DH_V_delta_z[index] += nlm[2];

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

void LCAO_Descriptor::cal_f_delta_hf_new(const ModuleBase::matrix& dm)
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

                        }//iw2
                    }//iw1
                }//ad2
            }//ad1
        }//end I0
    }//end T0
}

#endif