#include "cell_unittest.h"

#ifdef __LCAO
namespace Test_Cell
{
	Grid_Driver GridD(GlobalV::test_deconstructor, GlobalV::test_grid_driver,GlobalV::test_grid);
}

test_cell_orb::test_cell_orb()
{}
test_cell_orb::~test_cell_orb()
{}
#else
test_cell::test_cell()
{}
test_cell::~test_cell()
{}
#endif

#ifdef __LCAO
void test_cell_orb::set_parameters()
#else
void test_cell::set_parameters()
#endif
{
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::global_pseudo_type= "auto";
	GlobalV::PSEUDORCUT = 15.0;
	GlobalV::global_out_dir="./";
	GlobalV::ofs_warning.open("warning.log");
	ofs_running.open("log.txt");

	ucell.latName = "test";
	ucell.ntype = ntype;
	return;
}

#ifdef __LCAO
void test_cell_orb::setup_cell()
#else
void test_cell::setup_cell()
#endif
{
	ucell.setup_cell(
#ifdef __LCAO
	ORB,
#endif
	"./",
	out,
	"STRU", 
	ofs_running);

	return;
}

#ifndef __LCAO
void test_cell::count_ntype()
{
	std::cout << "count number of atom types" << std::endl;
	std::ifstream ifs("STRU",std::ios::in);

	if (!ifs)
	{
		std::cout << "ERROR : file STRU does not exist" <<std::endl;
		exit(1);
	}

	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs,"ATOMIC_SPECIES");	
	
	ntype = 0;

	std::string x;
	ifs.rdstate();
	while ( ifs.good() )
	{
//read a line
		std::getline(ifs,x);

//trim white space
		const char* typeOfWhitespaces = " \t\n\r\f\v";
		x.erase(x.find_last_not_of(typeOfWhitespaces) + 1);
		x.erase(0,x.find_first_not_of(typeOfWhitespaces));

		if(x=="LATTICE_CONSTANT" || x=="NUMERICAL_ORBITAL" || x=="LATTICE_VECTORS" || x=="ATOMIC_POSITIONS") break;

		std::string tmpid=x.substr(0,1);
		if(!x.empty() && tmpid!="#") ntype++;
	}

	std::cout << "ntype : "<< ntype << std::endl;
	ifs.close();

	return;
}
#endif

#ifdef __LCAO
void test_cell_orb::build_ST_new(const char& dtype, const bool& calc_deri)
{
/*
    //array to store data
    double olm[3]={0.0,0.0,0.0};
	int nnr = 0; // used onlyh for k points.

    //\sum{T} e**{ikT} <\phi_{ia}|d\phi_{k\beta}(T)>
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
    for (int T1=0; T1<ucell.ntype; ++T1)
    {
		Atom* atom1 = &ucell.atoms[T1];
        for (int I1=0; I1<atom1->na; ++I1)
        {
			tau1 = atom1->tau[I1];
            //Test_Cell::GridD.Find_atom(tau1);
            Test_Cell::GridD.Find_atom(ucell, tau1, T1, I1);
            for (int ad = 0; ad < Test_Cell::GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = Test_Cell::GridD.getType(ad);
				const int I2 = Test_Cell::GridD.getNatom(ad);
				Atom* atom2 = &ucell.atoms[T2];
				tau2 = Test_Cell::GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
				if(distance < rcut)
				{
					int iw1_all = ucell.itiaiw2iwt( T1, I1, 0) ; //iw1_all = combined index (it, ia, iw)

					for(int jj=0; jj<atom1->nw*GlobalV::NPOL; ++jj)
					{
						const int jj0 = jj/GlobalV::NPOL;
						const int L1 = atom1->iw2l[jj0];
						const int N1 = atom1->iw2n[jj0];
						const int m1 = atom1->iw2m[jj0];

						int iw2_all = ucell.itiaiw2iwt( T2, I2, 0);//zhengdy-soc
						for(int kk=0; kk<atom2->nw*GlobalV::NPOL; ++kk)
						{
							const int kk0 = kk/GlobalV::NPOL;
							const int L2 = atom2->iw2l[kk0];
							const int N2 = atom2->iw2n[kk0];
							const int m2 = atom2->iw2m[kk0];

							// mohan add 2010-06-29
							// this is in fact the same as in build_Nonlocal_mu,
							// the difference is that here we use {L,N,m} for ccycle,
							// build_Nonlocal_mu use atom.nw for cycle.
							// so, here we use ParaO::in_this_processor,
							// in build_Non... use trace_loc_row
							// and trace_loc_col directly,
							if ( !GlobalC::ParaO.in_this_processor(iw1_all,iw2_all) )
							{
								++iw2_all;
								continue;
							}

							olm[0] = olm[1] = olm[2] = 0.0;

							std::complex<double> olm1[4]={ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO};
							std::complex<double> *olm2 = &olm1[0];
							if(!calc_deri)
							{
								// PLEASE use UOT as an input parameter of this subroutine
								// mohan add 2021-03-30
								GlobalC::UOT.snap_psipsi( olm, 0, dtype, tau1, 
										T1, L1, m1, N1, Test_Cell::GridD.getAdjacentTau(ad), 
										T2, L2, m2, N2, GlobalV::NSPIN,
										olm2//for soc
										);

								if(GlobalV::GAMMA_ONLY_LOCAL)
								{
									// mohan add 2010-06-29
									// set the value in Hloc and Sloc
									// according to trace_loc_row and trace_loc_col
									// the last paramete: 1 for Sloc, 2 for Hloc
									// and 3 for Hloc_fixed.
									GlobalC::LM.set_HSgamma(iw1_all, iw2_all, olm[0], dtype);
								}
								else // k point algorithm
								{
									// mohan add 2010-10
									// set the values in SlocR and Hloc_fixedR.
									// which is a 1D array.
									if(dtype=='S')
									{
										if(GlobalV::NSPIN!=4) GlobalC::LM.SlocR[nnr] = olm[0];
										else
										{//only has diagonal term here.
												int is = (jj-jj0*GlobalV::NPOL) + (kk-kk0*GlobalV::NPOL)*2;
											GlobalC::LM.SlocR_soc[nnr] = olm1[is];
										}
									}
									else if(dtype=='T')
									{
										if(GlobalV::NSPIN!=4) GlobalC::LM.Hloc_fixedR[nnr] = olm[0];// <phi|kin|d phi>
										else
										{//only has diagonal term here.
												int is = (jj-jj0*GlobalV::NPOL) + (kk-kk0*GlobalV::NPOL)*2;
											GlobalC::LM.Hloc_fixedR_soc[nnr] = olm1[is];
										}
									}
									++nnr;
								}
							}
							else // calculate the derivative
							{
								GlobalC::UOT.snap_psipsi( olm, 1, dtype, 
									tau1, T1, L1, m1, N1,
									Test_Cell::GridD.getAdjacentTau(ad), T2, L2, m2, N2, GlobalV::NSPIN
									);

								if(GlobalV::GAMMA_ONLY_LOCAL)
								{
									GlobalC::LM.set_force (iw1_all, iw2_all,	olm[0], olm[1], olm[2], dtype);
									if(GlobalV::STRESS) GlobalC::LM.set_stress (iw1_all, iw2_all, olm[0], olm[1], olm[2], dtype, dtau);
								}
								else // k point algorithm
								{
									if(dtype=='S')
									{
										GlobalC::LM.DSloc_Rx[nnr] = olm[0];
										GlobalC::LM.DSloc_Ry[nnr] = olm[1];
										GlobalC::LM.DSloc_Rz[nnr] = olm[2];
										if(GlobalV::STRESS)
										{
											GlobalC::LM.DH_r[nnr*3] = dtau.x;
											GlobalC::LM.DH_r[nnr*3 + 1] = dtau.y;
											GlobalC::LM.DH_r[nnr*3 + 2] = dtau.z;
										}
									}
									else if(dtype=='T')
									{
										// notice the 'sign'
										GlobalC::LM.DHloc_fixedR_x[nnr] = olm[0];
										GlobalC::LM.DHloc_fixedR_y[nnr] = olm[1];
										GlobalC::LM.DHloc_fixedR_z[nnr] = olm[2];
										if(GlobalV::STRESS)
										{
											GlobalC::LM.stvnl11[nnr] = olm[0] * dtau.x;
											GlobalC::LM.stvnl12[nnr] = olm[0] * dtau.y;
											GlobalC::LM.stvnl13[nnr] = olm[0] * dtau.z;
											GlobalC::LM.stvnl22[nnr] = olm[1] * dtau.y;
											GlobalC::LM.stvnl23[nnr] = olm[1] * dtau.z;
											GlobalC::LM.stvnl33[nnr] = olm[2] * dtau.z;
										}
									}
									++nnr;
								}
							}
							++iw2_all;
						}// nw2 
						++iw1_all;
					}// nw1
				}// distance
				else if(distance>=rcut && (!GlobalV::GAMMA_ONLY_LOCAL))
				{
					int start1 = ucell.itiaiw2iwt( T1, I1, 0);
					int start2 = ucell.itiaiw2iwt( T2, I2, 0);

					bool is_adj = false;
					for (int ad0=0; ad0 < Test_Cell::GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = Test_Cell::GridD.getType(ad0);
						//const int I0 = Test_Cell::GridD.getNatom(ad0);
						//const int iat0 = ucell.itia2iat(T0, I0);
						//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
						tau0 = Test_Cell::GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						double distance1 = dtau1.norm() * ucell.lat0;
						double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Beta[T0].get_rcut_max();
						dtau2 = tau0 - tau2;
						double distance2 = dtau2.norm() * ucell.lat0;
						double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ORB.Beta[T0].get_rcut_max();
						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							is_adj = true;
							break;
						}
					}//ad0


					if( is_adj )
					{
						for(int jj=0; jj<atom1->nw * GlobalV::NPOL; ++jj)
						{
							const int mu = GlobalC::ParaO.trace_loc_row[start1+jj];
							if(mu<0)continue; 
							for(int kk=0; kk<atom2->nw * GlobalV::NPOL; ++kk)
							{
								const int nu = GlobalC::ParaO.trace_loc_col[start2+kk];
								if(nu<0)continue;
								++nnr;
							}//kk
						}//jj
					}
				}//distance
			}// ad
		}// I1
	}// T1

	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
		if(nnr != GlobalC::LNNR.nnr)
		{
			std::cout << " nnr=" << nnr << " GlobalC::LNNR.nnr=" << GlobalC::LNNR.nnr << std::endl;
			GlobalV::ofs_running << " nnr=" << nnr << " GlobalC::LNNR.nnr=" << GlobalC::LNNR.nnr << std::endl;
			ModuleBase::WARNING_QUIT("LCAO_gen_fixedH::build_ST_new","nnr != GlobalC::LNNR.nnr");
		}
	}
*/
    return;
}
#endif
