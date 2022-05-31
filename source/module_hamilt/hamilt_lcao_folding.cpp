#include "hamilt_lcao.h"

// be called in LCAO_Hamilt::calculate_Hk.
template<typename T, typename T1>
void HamiltLCAO<T, T1>::updateHk(const int &ik)
{
	ModuleBase::TITLE("HamiltLCAO","updateHk");
    ModuleBase::timer::tick("HamiltLCAO","updateHk");
    const Parallel_Orbitals* pv = this->LM->ParaV;

	int iat = 0;
	int index = 0;
	ModuleBase::Vector3<double> dtau;
	ModuleBase::Vector3<double> tau1;
	ModuleBase::Vector3<double> tau2;

	ModuleBase::Vector3<double> dtau1;
	ModuleBase::Vector3<double> dtau2;
	ModuleBase::Vector3<double> tau0;

#ifdef __DEEPKS
	if (GlobalV::deepks_scf)
    {
		ModuleBase::GlobalFunc::ZEROS(GlobalC::ld.H_V_delta_k[ik], pv->nloc);
	}
#endif

	for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
	{
		Atom* atom1 = &GlobalC::ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			//GlobalC::GridD.Find_atom(tau1);
			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
			Atom* atom1 = &GlobalC::ucell.atoms[T1];
			const int start = GlobalC::ucell.itiaiw2iwt(T1,I1,0);

			// (2) search among all adjacent atoms.
			for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
			{
				const int T2 = GlobalC::GridD.getType(ad);
				const int I2 = GlobalC::GridD.getNatom(ad);
				Atom* atom2 = &GlobalC::ucell.atoms[T2];

				tau2 = GlobalC::GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * GlobalC::ucell.lat0;
				double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

				bool adj = false;

				if(distance < rcut) 
				{
					adj = true;
				}
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GlobalC::GridD.getType(ad0); 
						const int I0 = GlobalC::GridD.getNatom(ad0); 
						//const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GlobalC::GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;

						double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
						double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

						double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
						double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							adj = true;
							break;
						}
					}
				}

				if(adj) // mohan fix bug 2011-06-26, should not be '<='
				{
					// (3) calculate the nu of atom (T2, I2)
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);
					//------------------------------------------------
					// exp(k dot dR)
					// dR is the index of box in Crystal coordinates
					//------------------------------------------------
					ModuleBase::Vector3<double> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z); 
					const double arg = ( GlobalC::kv.kvec_d[ik] * dR ) * ModuleBase::TWO_PI;
					//const double arg = ( GlobalC::kv.kvec_d[ik] * GlobalC::GridD.getBox(ad) ) * ModuleBase::TWO_PI;
					const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );

					//--------------------------------------------------
					// calculate how many matrix elements are in 
					// this processor.
					//--------------------------------------------------
					for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
					{
						// the index of orbitals in this processor
						const int iw1_all = start + ii;
						const int mu = pv->trace_loc_row[iw1_all];
						if(mu<0)continue;

						for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
						{
							int iw2_all = start2 + jj;
							const int nu = pv->trace_loc_col[iw2_all];

							if(nu<0)continue;
							//const int iic = mu*pv->ncol+nu;
                            int iic;
                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                            {
                                iic=mu+nu*pv->nrow;
                            }
                            else
                            {
                                iic=mu*pv->ncol+nu;
                            }

							//########################### EXPLAIN ###############################
							// 1. overlap matrix with k point
							// this->SlocR = < phi_0i | phi_Rj >, where 0, R are the cell index
							// while i,j are the orbital index.

							// 2. H_fixed=T+Vnl matrix element with k point (if Vna is not used).
							// H_fixed=T+Vnl+Vna matrix element with k point (if Vna is used).
							// this->Hloc_fixed = < phi_0i | H_fixed | phi_Rj>

							// 3. H(k) |psi(k)> = S(k) | psi(k)> 
							// Sloc2 is used to diagonalize for a give k point.
							// Hloc_fixed2 is used to diagonalize (eliminate index R).
							//###################################################################
							
							if(GlobalV::NSPIN!=4)
							{
								this->Sloc2[iic] += this->SlocR[index] * kphase;
								this->Hloc_fixed2[iic] += this->Hloc_fixedR[index] * kphase;
#ifdef __DEEPKS
								if(GlobalV::deepks_scf)
								{
									GlobalC::ld.H_V_delta_k[ik][iic] += GlobalC::ld.H_V_deltaR[index] * kphase;
								}
#endif
							}
							else
							{
								this->Sloc2[iic] += this->SlocR_soc[index] * kphase;
								this->Hloc_fixed2[iic] += this->Hloc_fixedR_soc[index] * kphase;
							}
							++index;

						}//end jj
					}//end ii
				}
			}// end ad
			++iat;
		}// end I1
	} // end T1

	assert(index==this->LM->ParaV->nnr);

	ModuleBase::timer::tick("HamiltLCAO","updateHk");
	return;
}