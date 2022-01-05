//caoyu add 2021-03-2
#ifdef __DEEPKS

#include "LCAO_deepks.h"

//===============================
//DeePKS Part 1
//deals with generation of descriptors as well as labels
//===============================

namespace GlobalC
{
    LCAO_Deepks ld;
}

LCAO_Deepks::LCAO_Deepks()
{
    alpha_index = new ModuleBase::IntArray[1];
    inl_index = new ModuleBase::IntArray[1];
    inl_l = new int[1];
    H_V_delta = new double[1];
}

LCAO_Deepks::~LCAO_Deepks()
{
    delete[] alpha_index;
    delete[] inl_index;
    delete[] inl_l;
    delete[] H_V_delta;

    //=======1. "out_descriptor" part==========
    //delete pdm**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] pdm[inl];
    }
    delete[] pdm;
    //=======2. "deepks_scf" part==========
    if (GlobalV::deepks_scf)
    {
        //delete gedm**
        for (int inl = 0;inl < this->inlmax;inl++)
        {
            delete[] gedm[inl];
        }
        delete[] gedm;
    }
}

#ifdef __MPI
void LCAO_Deepks::allsum_deepks(
    int inlmax, //first dimension
    int ndim, //second dimension
    double** mat) //the array being reduced 
{
    for(int inl=0;inl<inlmax;inl++)
    {
        Parallel_Reduce::reduce_double_all(mat[inl],ndim);
    }
}
#endif

void LCAO_Deepks::init(
    const LCAO_Orbitals &orb,
    const int nat,
    const int ntype,
    std::vector<int> na) 
{
    ModuleBase::TITLE("LCAO_Deepks", "init");

    GlobalV::ofs_running << " Initialize the descriptor index for DeePKS (lcao line)" << std::endl;

    const int lm = orb.get_lmax_d();
    const int nm = orb.get_nchimax_d();
    const int tot_inl_per_atom = orb.Alpha[0].getTotal_nchi();

    assert(lm >= 0);
    assert(nm >= 0);
    assert(tot_inl_per_atom >= 0);
    
    const int tot_inl = tot_inl_per_atom * nat;

    this->lmaxd = lm;
    this->nmaxd = nm;
    this->inlmax = tot_inl;
    GlobalV::ofs_running << " lmax of descriptor = " << this->lmaxd << std::endl;
    GlobalV::ofs_running << " nmax of descriptor= " << nmaxd << std::endl;
	GlobalV::ofs_running << " total basis (all atoms) for descriptor= " << std::endl;
    
    //init pdm**
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    this->pdm = new double* [this->inlmax];
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->pdm[inl] = new double[pdm_size];
        ModuleBase::GlobalFunc::ZEROS(this->pdm[inl], pdm_size);
    }

    // cal n(descriptor) per atom , related to Lmax, nchi(L) and m. (not total_nchi!)
	this->des_per_atom=0; // mohan add 2021-04-21
    for (int l = 0; l <= this->lmaxd; l++)
    {
        this->des_per_atom += orb.Alpha[0].getNchi(l) * (2 * l + 1);
    }

    this->n_descriptor = nat * this->des_per_atom;

    this->init_index(ntype, nat, na, tot_inl, orb);
    this->allocate_nlm(nat);

    return;
}

void LCAO_Deepks::init_index(const int ntype, const int nat, std::vector<int> na, const int Total_nchi, const LCAO_Orbitals &orb)
{
    delete[] this->alpha_index;
    this->alpha_index = new ModuleBase::IntArray[ntype];
    delete[] this->inl_index;
    this->inl_index = new ModuleBase::IntArray[ntype];
    delete[] this->inl_l;
    this->inl_l = new int[this->inlmax];
    ModuleBase::GlobalFunc::ZEROS(this->inl_l, this->inlmax);

    int inl = 0;
    int alpha = 0;
    for (int it = 0; it < ntype; it++)
    {
        this->alpha_index[it].create(
            na[it],
            this->lmaxd + 1, // l starts from 0
            this->nmaxd,
            2 * this->lmaxd + 1); // m ==> 2*l+1

        this->inl_index[it].create(
            na[it],
            this->lmaxd + 1,
            this->nmaxd); 

        GlobalV::ofs_running << " Type " << it + 1
                    << " number_of_atoms " << na[it] << std::endl;

        for (int ia = 0; ia < na[it]; ia++)
        {
            //alpha
            for (int l = 0; l < this->lmaxd + 1; l++)
            {
                for (int n = 0; n < orb.Alpha[0].getNchi(l); n++)
                {
                    for (int m = 0; m < 2 * l + 1; m++)
                    {
                        this->alpha_index[it](ia, l, n, m) = alpha;
                        alpha++;
                    }
                    this->inl_index[it](ia, l, n) = inl;
                    this->inl_l[inl] = l;
                    inl++;
                }
            }
        }//end ia
    }//end it
    assert(this->n_descriptor == alpha);
    assert(Total_nchi == inl);
    GlobalV::ofs_running << " descriptors_per_atom " << this->des_per_atom << std::endl;
    GlobalV::ofs_running << " total_descriptors " << this->n_descriptor << std::endl;
	return;
}

//this subroutine performs the calculation of projected density matrices
//pdm_m,m'=\sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>
void LCAO_Deepks::cal_projected_DM(const ModuleBase::matrix &dm, 
    const UnitCell_pseudo &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver &GridD,
    const Parallel_Orbitals &ParaO)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_projected_DM");
    ModuleBase::timer::tick("LCAO_Deepks","cal_projected_DM"); 
    if(dm.nr == 0 && dm.nc ==0)
    {
        return;
    }

    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    for(int inl=0;inl<inlmax;inl++)
    {
        ModuleBase::GlobalFunc::ZEROS(pdm[inl],pdm_size);
    }

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = orb.Phi[T1].getRcut(); 

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
					const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1;
						const int iw1_local = ParaO.trace_loc_col[iw1_all];
						if(iw1_local < 0)continue;
						const int iw1_0 = iw1/GlobalV::NPOL;

						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2;
							const int iw2_local = ParaO.trace_loc_row[iw2_all];
							if(iw2_local < 0)continue;
							const int iw2_0 = iw2/GlobalV::NPOL;

                            std::vector<double> nlm1 = this->nlm_save[iat][ad1][iw1_all][0];
                            std::vector<double> nlm2 = this->nlm_save[iat][ad2][iw2_all][0];
                            assert(nlm1.size()==nlm2.size());

                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                    {
                                        for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                        {
                                            int ind = m1*nm + m2;
                                            pdm[inl][ind] += dm(iw1_local, iw2_local)*nlm1[ib+m1]*nlm2[ib+m2];
                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());               
						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }//I0
    }//T0

#ifdef __MPI
    allsum_deepks(this->inlmax,pdm_size,this->pdm);
#endif
    ModuleBase::timer::tick("LCAO_Deepks","cal_projected_DM"); 
    return;
}

void LCAO_Deepks::cal_projected_DM_k(const std::vector<ModuleBase::ComplexMatrix>& dm,
    const UnitCell_pseudo &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver &GridD,
    const Parallel_Orbitals &ParaO,
    const K_Vectors &kv)
{

    ModuleBase::timer::tick("LCAO_Deepks","cal_projected_DM_k");

    if(dm[0].nr == 0 && dm[0].nc ==0)
    {
        return;
    }

    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    for(int inl=0;inl<inlmax;inl++)
    {
        ModuleBase::GlobalFunc::ZEROS(pdm[inl],pdm_size);
    }

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1,I1);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = orb.Phi[T1].getRcut();

                ModuleBase::Vector3<double> dR1(GridD.getBox(ad1).x, GridD.getBox(ad1).y, GridD.getBox(ad1).z); 

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
                    const int ibt2 = ucell.itia2iat(T2,I2);
					const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    ModuleBase::Vector3<double> dR2(GridD.getBox(ad2).x, GridD.getBox(ad2).y, GridD.getBox(ad2).z);
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1;
						const int iw1_local = ParaO.trace_loc_col[iw1_all];
						if(iw1_local < 0)continue;
						const int iw1_0 = iw1/GlobalV::NPOL;
						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2;
							const int iw2_local = ParaO.trace_loc_row[iw2_all];
							if(iw2_local < 0)continue;
							const int iw2_0 = iw2/GlobalV::NPOL;
 
                            double dm_current;
                            std::complex<double> tmp = 0.0;
                            for(int ik=0;ik<kv.nks;ik++)
                            {
                                const double arg = ( kv.kvec_d[ik] * (dR1-dR2) ) * ModuleBase::TWO_PI;
                                const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );
                                tmp += dm[ik](iw1_local,iw2_local)*kphase;
                            }
                            dm_current=tmp.real();

                            key_tuple key_1(ibt1,dR1.x,dR1.y,dR1.z);
                            key_tuple key_2(ibt2,dR2.x,dR2.y,dR2.z);
                            std::vector<double> nlm1 = this->nlm_save_k[iat][key_1][iw1_all][0];
                            std::vector<double> nlm2 = this->nlm_save_k[iat][key_2][iw2_all][0];
                            assert(nlm1.size()==nlm2.size());

                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                    {
                                        for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                        {
                                            int ind = m1*nm + m2;
                                            pdm[inl][ind] += dm_current*nlm1[ib+m1]*nlm2[ib+m2];
                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());               
						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }//I0
    }//T0

#ifdef __MPI
    allsum_deepks(this->inlmax,pdm_size,this->pdm);
#endif
    ModuleBase::timer::tick("LCAO_Deepks","cal_projected_DM_k");
    return;
    
}

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

//this subroutine calculates the gradient of projected density matrices
//gdmx_m,m = d/dX sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>
void LCAO_Deepks::cal_gdmx(const ModuleBase::matrix &dm,
    const UnitCell_pseudo &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver &GridD,
    const Parallel_Orbitals &ParaO)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_gdmx");
    //get DS_alpha_mu and S_nu_beta

    int size = (2 * lmaxd + 1) * (2 * lmaxd + 1);
    for (int iat = 0;iat < ucell.nat;iat++)
    {
        for (int inl = 0;inl < inlmax;inl++)
        {
            ModuleBase::GlobalFunc::ZEROS(gdmx[iat][inl], size);
            ModuleBase::GlobalFunc::ZEROS(gdmy[iat][inl], size);
            ModuleBase::GlobalFunc::ZEROS(gdmz[iat][inl], size);
        }
    }

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);//on which alpha is located
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1,I1); //on which chi_mu is located
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = orb.Phi[T1].getRcut(); 

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
					const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    const int ibt2 = ucell.itia2iat(T2,I2);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1;
						const int iw1_local = ParaO.trace_loc_row[iw1_all];
						if(iw1_local < 0)continue;
						const int iw1_0 = iw1/GlobalV::NPOL;
						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2;
							const int iw2_local = ParaO.trace_loc_col[iw2_all];
							if(iw2_local < 0)continue;
							const int iw2_0 = iw2/GlobalV::NPOL;
                            
                            std::vector<double> nlm1 = this->nlm_save[iat][ad1][iw1_all][0];
                            std::vector<std::vector<double>> nlm2 = this->nlm_save[iat][ad2][iw2_all];

                            assert(nlm1.size()==nlm2[0].size());

                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                    {
                                        for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                        {
                                            //ibt : on which iw2 is located
                                            //iat : on which alpha is located

                                            double fact_x = (nlm2[1][ib+m2] * nlm1[ib+m1] + nlm2[1][ib+m2] * nlm1[ib+m1]) * dm(iw2_local,iw1_local);
                                            double fact_y = (nlm2[2][ib+m2] * nlm1[ib+m1] + nlm2[2][ib+m2] * nlm1[ib+m1]) * dm(iw2_local,iw1_local);
                                            double fact_z = (nlm2[3][ib+m2] * nlm1[ib+m1] + nlm2[3][ib+m2] * nlm1[ib+m1]) * dm(iw2_local,iw1_local);

                                            //(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'> + (<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                            gdmx[iat][inl][m1*nm + m2] += fact_x;
                                            gdmy[iat][inl][m1*nm + m2] += fact_y;                                               
                                            gdmz[iat][inl][m1*nm + m2] += fact_z;

                                            //(<chi_mu|d/dX alpha_m>)<chi_nu|alpha_m'> + (<chi_nu|d/dX alpha_m'>)<chi_mu|alpha_m>
                                            // = -(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'> - (<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                            gdmx[ibt2][inl][m1*nm + m2] -= fact_x;                                               
                                            gdmy[ibt2][inl][m1*nm + m2] -= fact_y;                                               
                                            gdmz[ibt2][inl][m1*nm + m2] -= fact_z;

                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());
						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }//I0
    }//T0

#ifdef __MPI
    for(int iat=0;iat<ucell.nat;iat++)
    {
        allsum_deepks(this->inlmax,size,this->gdmx[iat]);
        allsum_deepks(this->inlmax,size,this->gdmy[iat]);
        allsum_deepks(this->inlmax,size,this->gdmz[iat]);
    }
#endif
    return;
}

void LCAO_Deepks::cal_gdmx_k(const std::vector<ModuleBase::ComplexMatrix>& dm,
    const UnitCell_pseudo &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver &GridD,
    const Parallel_Orbitals &ParaO,
    const K_Vectors &kv)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_gdmx");

    int size = (2 * lmaxd + 1) * (2 * lmaxd + 1);

    for (int iat = 0;iat < ucell.nat;iat++)
    {
        for (int inl = 0;inl < inlmax;inl++)
        {
            ModuleBase::GlobalFunc::ZEROS(gdmx[iat][inl], size);
            ModuleBase::GlobalFunc::ZEROS(gdmy[iat][inl], size);
            ModuleBase::GlobalFunc::ZEROS(gdmz[iat][inl], size);
        }
    }

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);//on which alpha is located
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1,I1); //on which chi_mu is located
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = orb.Phi[T1].getRcut();

                ModuleBase::Vector3<double> dR1(GridD.getBox(ad1).x, GridD.getBox(ad1).y, GridD.getBox(ad1).z); 

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
					const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    const int ibt2 = ucell.itia2iat(T2,I2);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    ModuleBase::Vector3<double> dR2(GridD.getBox(ad2).x, GridD.getBox(ad2).y, GridD.getBox(ad2).z);
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1;
						const int iw1_local = ParaO.trace_loc_row[iw1_all];
						if(iw1_local < 0)continue;
						const int iw1_0 = iw1/GlobalV::NPOL;
						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2;
							const int iw2_local = ParaO.trace_loc_col[iw2_all];
							if(iw2_local < 0)continue;
							const int iw2_0 = iw2/GlobalV::NPOL;

                            double dm_current;
                            std::complex<double> tmp = 0.0;
                            for(int ik=0;ik<kv.nks;ik++)
                            {
                                const double arg = - ( kv.kvec_d[ik] * (dR1-dR2) ) * ModuleBase::TWO_PI;
                                const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );
                                tmp += dm[ik](iw2_local,iw1_local)*kphase;
                            }
                            dm_current=tmp.real();
                            
                            key_tuple key_1(ibt1,dR1.x,dR1.y,dR1.z);
                            key_tuple key_2(ibt2,dR2.x,dR2.y,dR2.z);
                            std::vector<double> nlm1 = this->nlm_save_k[iat][key_1][iw1_all][0];
                            std::vector<std::vector<double>> nlm2 = this->nlm_save_k[iat][key_2][iw2_all];

                            assert(nlm1.size()==nlm2[0].size());

                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                    {
                                        for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                        {
                                            //ibt : on which iw2 is located
                                            //iat : on which alpha is located

                                            double fact_x = (nlm2[1][ib+m2] * nlm1[ib+m1] + nlm2[1][ib+m2] * nlm1[ib+m1]) * dm_current;
                                            double fact_y = (nlm2[2][ib+m2] * nlm1[ib+m1] + nlm2[2][ib+m2] * nlm1[ib+m1]) * dm_current;
                                            double fact_z = (nlm2[3][ib+m2] * nlm1[ib+m1] + nlm2[3][ib+m2] * nlm1[ib+m1]) * dm_current;

                                            //(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'> + (<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                            gdmx[iat][inl][m1*nm + m2] += fact_x;
                                            gdmy[iat][inl][m1*nm + m2] += fact_y;                                               
                                            gdmz[iat][inl][m1*nm + m2] += fact_z;

                                            //(<chi_mu|d/dX alpha_m>)<chi_nu|alpha_m'> + (<chi_nu|d/dX alpha_m'>)<chi_mu|alpha_m>
                                            // = -(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'> - (<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                            gdmx[ibt2][inl][m1*nm + m2] -= fact_x;                                               
                                            gdmy[ibt2][inl][m1*nm + m2] -= fact_y;                                               
                                            gdmz[ibt2][inl][m1*nm + m2] -= fact_z;

                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());
						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }//I0
    }//T0

#ifdef __MPI
    for(int iat=0;iat<ucell.nat;iat++)
    {
        allsum_deepks(this->inlmax,size,this->gdmx[iat]);
        allsum_deepks(this->inlmax,size,this->gdmy[iat]);
        allsum_deepks(this->inlmax,size,this->gdmz[iat]);
    }
#endif
    return;
}

void LCAO_Deepks::init_gdmx(const int nat)
{
    this->gdmx = new double** [nat];
    this->gdmy = new double** [nat];
    this->gdmz = new double** [nat];
    for (int iat = 0;iat < nat;iat++)
    {
        this->gdmx[iat] = new double* [inlmax];
        this->gdmy[iat] = new double* [inlmax];
        this->gdmz[iat] = new double* [inlmax];
        for (int inl = 0;inl < inlmax;inl++)
        {
            this->gdmx[iat][inl] = new double [(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            this->gdmy[iat][inl] = new double [(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            this->gdmz[iat][inl] = new double[(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            ModuleBase::GlobalFunc::ZEROS(gdmx[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
            ModuleBase::GlobalFunc::ZEROS(gdmy[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
            ModuleBase::GlobalFunc::ZEROS(gdmz[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
        }
    }
    return;
}


void LCAO_Deepks::del_gdmx(const int nat)
{
    for (int iat = 0;iat < nat;iat++)
    {
        for (int inl = 0;inl < inlmax;inl++)
        {
            delete[] this->gdmx[iat][inl];
            delete[] this->gdmy[iat][inl];
            delete[] this->gdmz[iat][inl];
        }
        delete[] this->gdmx[iat];
        delete[] this->gdmy[iat];
        delete[] this->gdmz[iat];
    }
    delete[] this->gdmx;
    delete[] this->gdmy;
    delete[] this->gdmz;
    return;
}

#endif
