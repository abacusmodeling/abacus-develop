//wenfei 2022-1-5
//This file contains one subroutine for calculating the overlap
//between atomic basis and projector alpha : chi_mu|alpha>;
//it will be used in calculating pdm, gdmx, H_V_delta, F_delta

#ifdef __DEEPKS

#include "LCAO_deepks.h"

void LCAO_Deepks::build_psialpha(const bool& calc_deri,
    const UnitCell_pseudo &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver &GridD,
    const Parallel_Orbitals &ParaO,
    const ORB_gen_tables &UOT)
{
    ModuleBase::TITLE("LCAO_Deepks", "build_psialpha");
    ModuleBase::timer::tick ("LCAO_Deepks","build_psialpha");

    const double Rcut_Alpha = orb.Alpha[0].getRcut();
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

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);
			//=======================================================
            //Step 1 : 
			//saves <alpha|psi>, where alpha runs over all projectors
			//and psi runs over atomic basis sets on the current core
			//=======================================================

			const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

			//outermost loop : all adjacent atoms
            if(GlobalV::GAMMA_ONLY_LOCAL)
            {
                this->nlm_save[iat].resize(GridD.getAdjacentNum()+1);
            }

            for (int ad=0; ad<GridD.getAdjacentNum()+1 ; ++ad)
            {
                const int T1 = GridD.getType(ad);
                const int I1 = GridD.getNatom(ad);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
				const double Rcut_AO1 = orb.Phi[T1].getRcut();

                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad);
				const Atom* atom1 = &ucell.atoms[T1];
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

				const double dist1 = (tau1-tau0).norm() * ucell.lat0;

				if (dist1 > Rcut_Alpha + Rcut_AO1)
				{
					continue;
				}

				for (int iw1=0; iw1<nw1_tot; ++iw1)
				{
					const int iw1_all = start1 + iw1;
					const int iw1_local = ParaO.trace_loc_row[iw1_all];
					const int iw2_local = ParaO.trace_loc_col[iw1_all];
					if(iw1_local < 0 && iw2_local < 0)continue;
					const int iw1_0 = iw1/GlobalV::NPOL;
					std::vector<std::vector<double>> nlm;
					//2D, but first dimension is only 1 here
					//for force, the right hand side is the gradient
					//and the first dimension is then 3
					//inner loop : all projectors (L0,M0)
					UOT.snap_psialpha_half(
						nlm, job, tau1, T1,
						atom1->iw2l[ iw1_0 ], // L1
						atom1->iw2m[ iw1_0 ], // m1
						atom1->iw2n[ iw1_0 ], // N1
						ucell.atoms[T0].tau[I0], T0, I0); //R0,T0

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
                    const int ibt=ucell.itia2iat(T1, I1);
                    const int rx=GridD.getBox(ad).x;
                    const int ry=GridD.getBox(ad).y;
                    const int rz=GridD.getBox(ad).z;
                    key_tuple key_1(ibt,rx,ry,rz);
                    this->nlm_save_k[iat][key_1]=nlm_cur;
                }
			}//end ad
		}//end I0
	}//end T0

    ModuleBase::timer::tick ("LCAO_Deepks","build_psialpha");
	return;

}

#endif
