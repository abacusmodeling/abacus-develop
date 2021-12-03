#include "FORCE_k.h"
#include "../src_pw/global.h"
#include <unordered_map>
#include <map>

Force_LCAO_k::Force_LCAO_k ()
{
}

Force_LCAO_k::~Force_LCAO_k ()
{
}

#include "LCAO_nnr.h"
// be called in Force_LCAO::start_force_calculation
void Force_LCAO_k::ftable_k (
		const bool isforce,
		const bool isstress,
		ModuleBase::matrix& foverlap,
		ModuleBase::matrix& ftvnl_dphi,
		ModuleBase::matrix& fvnl_dbeta,	
		ModuleBase::matrix& fvl_dphi,
		ModuleBase::matrix& soverlap,
		ModuleBase::matrix& stvnl_dphi,
		ModuleBase::matrix& svnl_dbeta,
		ModuleBase::matrix& svl_dphi
		)
{
    ModuleBase::TITLE("Force_LCAO_k", "ftable_k");
	ModuleBase::timer::tick("Force_LCAO_k","ftable_k");
	
	this->allocate_k();

	// calculate the energy density matrix
	// and the force related to overlap matrix and energy density matrix.
	this->cal_foverlap_k(isforce, isstress, foverlap, soverlap);

	// calculate the density matrix
	double** dm2d = new double*[GlobalV::NSPIN];
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		dm2d[is] = new double[GlobalC::LNNR.nnr];
		ModuleBase::GlobalFunc::ZEROS(dm2d[is], GlobalC::LNNR.nnr);
	}
	ModuleBase::Memory::record ("Force_LCAO_k", "dm2d", GlobalV::NSPIN*GlobalC::LNNR.nnr, "double");	
	bool with_energy = false;

	
	this->set_EDM_k(dm2d, with_energy);
	
	this->cal_ftvnl_dphi_k(dm2d, isforce, isstress, ftvnl_dphi, stvnl_dphi);


	// ---------------------------------------
	// doing on the real space grid.
	// ---------------------------------------
	this->cal_fvl_dphi_k(dm2d, isforce, isstress, fvl_dphi, svl_dphi);

	this->calFvnlDbeta(dm2d, isforce, isstress, fvnl_dbeta, svnl_dbeta, GlobalV::vnl_method);

	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		delete[] dm2d[is];
	}
	delete[] dm2d;

	//----------------------------------------------------------------
	// reduce the force according to 2D distribution of H & S matrix.
	//----------------------------------------------------------------
    if(isforce)
	{
        Parallel_Reduce::reduce_double_pool( foverlap.c, foverlap.nr * foverlap.nc);
        Parallel_Reduce::reduce_double_pool( ftvnl_dphi.c, ftvnl_dphi.nr * ftvnl_dphi.nc);
        Parallel_Reduce::reduce_double_pool( fvnl_dbeta.c, fvnl_dbeta.nr * fvnl_dbeta.nc);
        Parallel_Reduce::reduce_double_pool( fvl_dphi.c, fvl_dphi.nr * fvl_dphi.nc);
	}
    if(isstress)
    {
        Parallel_Reduce::reduce_double_pool( soverlap.c, soverlap.nr * soverlap.nc);
        Parallel_Reduce::reduce_double_pool( stvnl_dphi.c, stvnl_dphi.nr * stvnl_dphi.nc);
        Parallel_Reduce::reduce_double_pool( svnl_dbeta.c, svnl_dbeta.nr * svnl_dbeta.nc);
        Parallel_Reduce::reduce_double_pool( svl_dphi.c, svl_dphi.nr * svl_dphi.nc);
    }

	// test the force.
	/*
	std::cout << " overlap force" << std::endl;
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
	{
		const double fac = ModuleBase::Ry_to_eV / 0.529177;
		std::cout << std::setw(5) << iat+1 << std::setw(15) << foverlap[iat][0] *fac<< std::setw(15) << foverlap[iat][1]*fac << 
		std::setw(15) << foverlap[iat][2]*fac << std::endl;
	}
	*/

	this->finish_k();

	ModuleBase::timer::tick("Force_LCAO_k","ftable_k");
    return;
}

void Force_LCAO_k::allocate_k(void)
{
	ModuleBase::TITLE("Force_LCAO_k","allocate_k");
	ModuleBase::timer::tick("Force_LCAO_k","allocate_k");

	const int nnr = GlobalC::LNNR.nnr;
	//--------------------------------
    // (1) allocate for dSx dSy & dSz
	//--------------------------------
	GlobalC::LM.DSloc_Rx = new double [nnr];
    GlobalC::LM.DSloc_Ry = new double [nnr];
    GlobalC::LM.DSloc_Rz = new double [nnr];
    ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_Rx, nnr);
    ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_Ry, nnr);
    ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_Rz, nnr);
	ModuleBase::Memory::record("force_lo", "dS", nnr*3, "double");
    
	if(GlobalV::STRESS){
		GlobalC::LM.DH_r = new double [3* nnr];
		ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DH_r, 3 * nnr);
		GlobalC::LM.stvnl11 = new double [nnr];
		GlobalC::LM.stvnl12 = new double [nnr];
		GlobalC::LM.stvnl13 = new double [nnr];
		GlobalC::LM.stvnl22 = new double [nnr];
		GlobalC::LM.stvnl23 = new double [nnr];
		GlobalC::LM.stvnl33 = new double [nnr];
		ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.stvnl11,  nnr);
		ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.stvnl12,  nnr);
		ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.stvnl13,  nnr);
		ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.stvnl22,  nnr);
		ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.stvnl23,  nnr);
		ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.stvnl33,  nnr);
		ModuleBase::Memory::record("stress_lo", "dSR", nnr*6, "double");
	}

	//-----------------------------
	// calculate dS = <phi | dphi> 
	//-----------------------------
    // tips: build_ST_new --> GlobalC::ParaO.set_force 
	bool cal_deri = true;
	GlobalC::UHM.genH.build_ST_new ('S', cal_deri, GlobalC::ucell);

	//-----------------------------------------
	// (2) allocate for <phi | T + Vnl | dphi>
	//-----------------------------------------
    GlobalC::LM.DHloc_fixedR_x = new double [nnr];
    GlobalC::LM.DHloc_fixedR_y = new double [nnr];
    GlobalC::LM.DHloc_fixedR_z = new double [nnr];
    ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixedR_x, nnr);
    ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixedR_y, nnr);
    ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixedR_z, nnr);
	ModuleBase::Memory::record("force_lo", "dTVNL", nnr*3, "double");
    
    // calculate dT=<phi|kin|dphi> in LCAO
    // calculate T + VNL(P1) in LCAO basis
    GlobalC::UHM.genH.build_ST_new ('T', cal_deri, GlobalC::ucell);
	//test(GlobalC::LM.DHloc_fixedR_x,"GlobalC::LM.DHloc_fixedR_x T part");
   
   	// calculate dVnl=<phi|dVnl|dphi> in LCAO 
	this->NonlocalDphi(GlobalV::NSPIN, GlobalV::vnl_method, cal_deri);
	//test(GlobalC::LM.DHloc_fixedR_x,"GlobalC::LM.DHloc_fixedR_x Vnl part");

	ModuleBase::timer::tick("Force_LCAO_k","allocate_k");
	return;
}

void Force_LCAO_k::finish_k(void)
{
    delete [] GlobalC::LM.DSloc_Rx;
    delete [] GlobalC::LM.DSloc_Ry;
    delete [] GlobalC::LM.DSloc_Rz;
    delete [] GlobalC::LM.DHloc_fixedR_x;
    delete [] GlobalC::LM.DHloc_fixedR_y;
    delete [] GlobalC::LM.DHloc_fixedR_z;
	if(GlobalV::STRESS)
	{
		delete [] GlobalC::LM.DH_r;
		delete [] GlobalC::LM.stvnl11;
		delete [] GlobalC::LM.stvnl12;
		delete [] GlobalC::LM.stvnl13;
		delete [] GlobalC::LM.stvnl22;
		delete [] GlobalC::LM.stvnl23;
		delete [] GlobalC::LM.stvnl33;
	}
	return;
}

#include "record_adj.h"
#include "LCAO_nnr.h"
void Force_LCAO_k::set_EDM_k(double** dm2d, const bool with_energy)
{
	ModuleBase::TITLE("Force_LCAO_k","set_EDM_k");
	ModuleBase::timer::tick("Force_LCAO_k","set_EDM_k");

	ModuleBase::Vector3<double> tau1, dtau;

	//----------------------------------------------------------
	// RA will set the adjacent information for each atom
	// 2d means this adjacent information is for HPSEPS's kind
	// of division of H matrix.
	//----------------------------------------------------------
 	//xiaohui add "OUT_LEVEL", 2015-09-16
	if(GlobalV::OUT_LEVEL != "m") GlobalV::ofs_running << " Calculate the energy density matrix with k " << std::endl;
	Record_adj RA;
	RA.for_2d();

	//------------------------
	// circle for each atom
	//------------------------
	for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
		for(int I1=0; I1<atom1->na; ++I1)
		{
			const int iat = GlobalC::ucell.itia2iat(T1,I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1,I1,0);
			const int gstart = GlobalC::LNNR.nlocstart[iat];
			const int irr = GlobalC::LNNR.nlocdim[iat];//number of adjacet orbitals

			std::complex<double> **vvv = new std::complex<double>*[GlobalV::NSPIN];
         //xiaohui add 2014-03-17, add "if(irr > 0)", 
         //meaing only allocate memory when number of iat-th atom adjacet orbitals are not zero in this processor
         if(irr > 0)
         {
   			for(int is=0; is<GlobalV::NSPIN; is++)
			   {
				   vvv[is] = new std::complex<double>[ irr ];
				   ModuleBase::GlobalFunc::ZEROS(vvv[is], irr );
			   }
         }

			int ispin=0;
			std::complex<double> *dm;
			//----------------
			// circle for k
			//----------------
			for(int ik=0; ik<GlobalC::kv.nks; ++ik)
			{
				// set the spin direction
				if(GlobalV::NSPIN==2)
				{
					ispin = GlobalC::kv.isk[ik];
				}
				dm = vvv[ispin];
				//------------------
				// circle for bands
				//------------------
				for(int ib=0; ib<GlobalV::NBANDS; ++ib)
				{
					const double w1=GlobalC::wf.wg(ik,ib);
					if(w1>0)
					{
						//-----------------------------
						// start the adjacent cycle.
						//-----------------------------
						std::complex<double> *wfc = GlobalC::LOWF.WFC_K[ik][ib];
						int count = 0;
						for (int cb = 0; cb < RA.na_each[iat]; ++cb)
						{
							const int T2 = RA.info[iat][cb][3];
							const int I2 = RA.info[iat][cb][4];
							Atom* atom2 = &GlobalC::ucell.atoms[T2];

							//-----------------
							// exp[i * R * k]
							//-----------------
							const std::complex<double> phase = w1 * exp(ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT * (
										GlobalC::kv.kvec_d[ik].x * RA.info[iat][cb][0] +
										GlobalC::kv.kvec_d[ik].y * RA.info[iat][cb][1] +
										GlobalC::kv.kvec_d[ik].z * RA.info[iat][cb][2]
										) );

							const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);

							for(int jj=0; jj<atom1->nw; ++jj)
							{
								const int iw1_all = start1 + jj;
								
								// 2D division (HPSEPS)
								const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
								if(mu<0) continue;
								const int mug = GlobalC::GridT.trace_lo[iw1_all];

								for(int kk=0; kk<atom2->nw; ++kk)
								{
									const int iw2_all = start2 + kk;
									
									// 2D division (HPSEPS)
									const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
									if(nu<0) continue;
									const int nug = GlobalC::GridT.trace_lo[iw2_all];
	
									if(mug >= 0 && nug >= 0)
									{
										dm[count] += set_EDM_k_element(phase, with_energy, 
										wfc[mug], wfc[nug], GlobalC::wf.ekb[ik][ib]); 
									}
									else if( mug >= 0 && nug <= 0)
									{
										const int a4 = GlobalC::LOWF.trace_aug[iw2_all];
									
										//assert(a4>=0);

										dm[count] += set_EDM_k_element(phase, with_energy, wfc[mug], 
										GlobalC::LOWF.WFC_K_aug[ik][ib][a4], GlobalC::wf.ekb[ik][ib]); 
									}
									else if( mug <= 0 && nug >= 0)
									{
										const int a3 = GlobalC::LOWF.trace_aug[iw1_all]; 

										dm[count] += set_EDM_k_element(phase, with_energy, 
										GlobalC::LOWF.WFC_K_aug[ik][ib][a3], wfc[nug], GlobalC::wf.ekb[ik][ib]); 
									}
									else if( mug <=0 && nug <=0 )
									{
										const int a1 = GlobalC::LOWF.trace_aug[iw1_all];
										const int a2 = GlobalC::LOWF.trace_aug[iw2_all];

										dm[count] += set_EDM_k_element(phase, with_energy, 
										GlobalC::LOWF.WFC_K_aug[ik][ib][a1], GlobalC::LOWF.WFC_K_aug[ik][ib][a2], GlobalC::wf.ekb[ik][ib]); 
									}
									assert(count<irr);
									++ count;
								}//kk
							}//jj
						}// cb
//						GlobalV::ofs_running << " count = " << count << std::endl;
						assert(count == GlobalC::LNNR.nlocdim[iat]);
					}// w1
				}//ib
			}//ik

			//--------------------------------------
			// get the real value density matrix or
			// energy density matrix 
			//--------------------------------------
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{
				for(int iv=0; iv<irr; ++iv)
				{
					dm2d[is][gstart+iv] = vvv[is][iv].real();
				}
            //xiaohui add 2014-03-17, add "if(irr > 0)"
            //meaning delete memory when number of iat-th atom adjacet orbitals are not zero in this processor 
            if(irr > 0)
            {
				   delete[] vvv[is];
            }
			}
			delete[] vvv;

		}// I1
	}// T1

RA.delete_grid();//xiaohui add 2015-02-04
	ModuleBase::timer::tick("Force_LCAO_k","set_EDM_k");
	return;
}


std::complex<double> Force_LCAO_k::set_EDM_k_element(
	const std::complex<double> &phase,
	const bool with_energy,
	std::complex<double> &coef1, std::complex<double> &coef2,
	const double &ekb)
{
	std::complex<double> dm = std::complex<double>(0,0);
	//--------------------------------------
	// for energy density matrix
	// \sum E(i)*exp(iRk)*psi(mu)*psi(nu)
	//--------------------------------------
	if(with_energy)
	{
		dm += phase * ekb * conj(coef1) * coef2;
	}
	//--------------------------------------
	// for density matrix
	// \sum E(i)*psi(mu)*psi(nu)
	//--------------------------------------
	else
	{
		dm += phase * conj(coef1) * coef2 ;
	}
	return dm;
}


void Force_LCAO_k::cal_foverlap_k(
	const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& foverlap, 
	ModuleBase::matrix& soverlap)
{
	ModuleBase::TITLE("Force_LCAO_k","cal_foverlap_k");
	ModuleBase::timer::tick("Force_LCAO_k","cal_foverlap_k");

	//--------------------------------------------
	// (1) allocate energy density matrix (nnr)
	//--------------------------------------------
	double** edm2d = new double*[GlobalV::NSPIN];
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		edm2d[is] = new double[GlobalC::LNNR.nnr];
		ModuleBase::GlobalFunc::ZEROS(edm2d[is], GlobalC::LNNR.nnr);
	}
	bool with_energy = true;

	//--------------------------------------------	
	// calculate the energy density matrix here.
	//--------------------------------------------	
	this->set_EDM_k(edm2d, with_energy);

	//--------------------------------------------
    //summation \sum_{i,j} E(i,j)*dS(i,j)
    //BEGIN CALCULATION OF FORCE OF EACH ATOM
	//--------------------------------------------
	ModuleBase::Vector3<double> tau1, dtau, tau2;

	Record_adj RA;
	RA.for_2d();

	int irr = 0;
	int iat = 0;
    for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1=0; I1<atom1->na; ++I1)
        {
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1,I1,0);
			for (int cb = 0; cb < RA.na_each[iat]; ++cb)
			{
				const int T2 = RA.info[iat][cb][3];
				const int I2 = RA.info[iat][cb][4];
				const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

				Atom* atom2 = &GlobalC::ucell.atoms[T2];

				for(int jj=0; jj<atom1->nw; jj++)
				{
					const int iw1_all = start1 + jj; 

					// HPSEPS
					const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
					if(mu<0)continue;

					for(int kk=0; kk<atom2->nw; kk++)
					{
						const int iw2_all = start2 + kk;

						// HPSEPS
						const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
						if(nu<0)continue;
						//==============================================================
						// here we use 'minus', but in GlobalV::GAMMA_ONLY_LOCAL we use 'plus',
						// both are correct because the 'DSloc_Rx' is used in 'row' (-),
						// however, the 'DSloc_x' in GAMMA is used in 'col' (+),
						// mohan update 2011-06-16
						//==============================================================
						for(int is=0; is<GlobalV::NSPIN; ++is)
						{
							double edm2d2 = 2.0 * edm2d[is][irr];
							if(isforce)
							{
								foverlap(iat,0) -= edm2d2 * GlobalC::LM.DSloc_Rx[irr];
								foverlap(iat,1) -= edm2d2 * GlobalC::LM.DSloc_Ry[irr];
								foverlap(iat,2) -= edm2d2 * GlobalC::LM.DSloc_Rz[irr];
							}
							if(isstress)
							{
								for(int ipol = 0;ipol<3;ipol++)
								{
									soverlap(0,ipol) += edm2d[is][irr] * GlobalC::LM.DSloc_Rx[irr] * GlobalC::LM.DH_r[irr * 3 + ipol];
									soverlap(1,ipol) += edm2d[is][irr] * GlobalC::LM.DSloc_Ry[irr] * GlobalC::LM.DH_r[irr * 3 + ipol];
									soverlap(2,ipol) += edm2d[is][irr] * GlobalC::LM.DSloc_Rz[irr] * GlobalC::LM.DH_r[irr * 3 + ipol];
								}
							}
						}
						++irr;
					}// end kk
				}// end jj
			}// end cb
			++iat;
		}
	}

	//-----------------
	// test the force
	//-----------------
	/*
	std::cout << " overlap force" << std::endl;
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
	{
		const double fac = ModuleBase::Ry_to_eV / 0.529177;
		std::cout << std::setw(5) << iat+1 << std::setw(15) << foverlap[iat][0] *fac<< std::setw(15) << foverlap[iat][1]*fac << 
		std::setw(15) << foverlap[iat][2]*fac << std::endl;
	}
	*/
	if(isstress){
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				soverlap(i,j) *=  GlobalC::ucell.lat0 / GlobalC::ucell.omega;
			}
		}
	}

	if(irr!=GlobalC::LNNR.nnr)
	{
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"wrong irr",irr);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"wrong GlobalC::LNNR.nnr",GlobalC::LNNR.nnr);
		ModuleBase::WARNING_QUIT("Force_LCAO_k::cal_foverlap_k","irr!=GlobalC::LNNR.nnr");
	}
	
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		delete[] edm2d[is];
	}
	delete[] edm2d;

	RA.delete_grid();//xiaohui add 2015-02-04
	ModuleBase::timer::tick("Force_LCAO_k","cal_foverlap_k");
	return;
}

void Force_LCAO_k::cal_ftvnl_dphi_k(
	double** dm2d, 
	const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& ftvnl_dphi, 
	ModuleBase::matrix& stvnl_dphi)
{	
	ModuleBase::TITLE("Force_LCAO_k","cal_ftvnl_dphi");
	ModuleBase::timer::tick("Force_LCAO_k","cal_ftvnl_dphi");
	
	// get the adjacent atom's information.

//	GlobalV::ofs_running << " calculate the ftvnl_dphi_k force" << std::endl;
	Record_adj RA;
	RA.for_2d();

	int irr = 0;
    for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1=0; I1<atom1->na; ++I1)
        {
			const int iat = GlobalC::ucell.itia2iat(T1,I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1,I1,0);
			for (int cb = 0; cb < RA.na_each[iat]; ++cb)
			{
				const int T2 = RA.info[iat][cb][3];
				const int I2 = RA.info[iat][cb][4];
				const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);
				Atom* atom2 = &GlobalC::ucell.atoms[T2];

				for(int jj=0; jj<atom1->nw; ++jj)
				{
					const int iw1_all = start1 + jj; 
					const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
					if(mu<0)continue;
					for(int kk=0; kk<atom2->nw; ++kk)
					{
						const int iw2_all = start2 + kk;
						const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
						if(nu<0)continue;
						//==============================================================
						// here we use 'minus', but in GlobalV::GAMMA_ONLY_LOCAL we use 'plus',
						// both are correct because the 'DSloc_Rx' is used in 'row' (-),
						// however, the 'DSloc_x' is used in 'col' (+),
						// mohan update 2011-06-16
						//==============================================================
						for(int is=0; is<GlobalV::NSPIN; ++is)
						{
							double dm2d2 = 2.0 * dm2d[is][irr];
							if(isforce)
							{
								ftvnl_dphi(iat,0) += dm2d2 * GlobalC::LM.DHloc_fixedR_x[irr];
								ftvnl_dphi(iat,1) += dm2d2 * GlobalC::LM.DHloc_fixedR_y[irr];
								ftvnl_dphi(iat,2) += dm2d2 * GlobalC::LM.DHloc_fixedR_z[irr];
							}
							if(isstress)
							{
								stvnl_dphi(0,0) -= dm2d[is][irr] * GlobalC::LM.stvnl11[irr];
								stvnl_dphi(0,1) -= dm2d[is][irr] * GlobalC::LM.stvnl12[irr];
								stvnl_dphi(0,2) -= dm2d[is][irr] * GlobalC::LM.stvnl13[irr];
								stvnl_dphi(1,1) -= dm2d[is][irr] * GlobalC::LM.stvnl22[irr];
								stvnl_dphi(1,2) -= dm2d[is][irr] * GlobalC::LM.stvnl23[irr];
								stvnl_dphi(2,2) -= dm2d[is][irr] * GlobalC::LM.stvnl33[irr];
							}
						}
						++irr;
					}//end kk
				}//end jj
			}// end cb
		}
	}
	assert(irr==GlobalC::LNNR.nnr);
	
//	test(GlobalC::LM.DSloc_Rx);
//	test(dm2d[0],"dm2d");

	if(isstress){
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				if(i<j) stvnl_dphi(j,i) = stvnl_dphi(i,j);
				stvnl_dphi(i,j) *=  GlobalC::ucell.lat0 / GlobalC::ucell.omega;
			}
		}
	}

	RA.delete_grid();//xiaohui add 2015-02-04
	ModuleBase::timer::tick("Force_LCAO_k","cal_ftvnl_dphi");
	return;
}

	
void Force_LCAO_k::test(double* mmm, const std::string &name)
{
	if(GlobalV::NPROC!=1)return;
	std::cout << "test!" << std::endl;

	int irr = 0;
	int ca = 0;

	GlobalV::ofs_running << " Calculate the test in Force_LCAO_k" << std::endl;
	Record_adj RA;
	RA.for_2d();
	
	double *test;
	test = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
	ModuleBase::GlobalFunc::ZEROS(test, GlobalV::NLOCAL *GlobalV::NLOCAL);
	
	for(int T1=0; T1<GlobalC::ucell.ntype; T1++)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1=0; I1<atom1->na; I1++)
        {
			//const int iat = GlobalC::ucell.itia2iat(T1,I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1,I1,0);
			for (int cb = 0; cb < RA.na_each[ca]; cb++ )
			{
				const int T2 = RA.info[ca][cb][3];
				const int I2 = RA.info[ca][cb][4];
				Atom* atom2 = &GlobalC::ucell.atoms[T2];
				const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);

				for(int jj=0; jj<atom1->nw; jj++)
				{
					const int iw1_all = start1+jj;	
					for(int kk=0; kk<atom2->nw; kk++)
					{
						const int iw2_all = start2+kk;
						assert(irr<GlobalC::LNNR.nnr);
						//test[iw1_all*GlobalV::NLOCAL+iw2_all] += GlobalC::LM.DHloc_fixedR_x[irr];
						test[iw1_all*GlobalV::NLOCAL+iw2_all] += mmm[irr];
						++irr;
					}
				}
			}
			++ca;
		}
	}
		
	std::cout << "\n " << name << std::endl;
	std::cout << std::setprecision(4);
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			if( abs(test[i*GlobalV::NLOCAL+j]) > 1.0e-5)
			std::cout << std::setw(12) << test[i*GlobalV::NLOCAL+j];
			else
			std::cout << std::setw(12) << "0";
		}
		std::cout << std::endl;
	}
	delete[] test;	

RA.delete_grid();//xiaohui add 2015-02-04
	return;
}


// must consider three-center H matrix.
void Force_LCAO_k::cal_fvnl_dbeta_k(
	double** dm2d, 
	const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& fvnl_dbeta, 
	ModuleBase::matrix& svnl_dbeta)
{
	ModuleBase::TITLE("Force_LCAO_k","cal_fvnl_dbeta_k");
	ModuleBase::timer::tick("Force_LCAO_k","cal_fvnl_dbeta_k");
	int iir = 0;
	ModuleBase::Vector3<double> tau1;
	ModuleBase::Vector3<double> tau2;
	ModuleBase::Vector3<double> dtau;
	ModuleBase::Vector3<double> tau0;
	ModuleBase::Vector3<double>	dtau1;
	ModuleBase::Vector3<double> dtau2;

	double rcut;
	double distance;

	double rcut1;
	double rcut2;
	double distance1;
	double distance2;
	
	for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
	{
		const Atom* atom1 = &GlobalC::ucell.atoms[T1];

		for (int I1 =0; I1< atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			//GlobalC::GridD.Find_atom( tau1 );
			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1 ,T1, I1);
			//const int iat = GlobalC::ucell.itia2iat(T1, I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);

			for (int ad2=0; ad2<GlobalC::GridD.getAdjacentNum()+1 ; ++ad2)
			{
				const int T2 = GlobalC::GridD.getType(ad2);
				const Atom* atom2 = &GlobalC::ucell.atoms[T2];
				const int I2 = GlobalC::GridD.getNatom(ad2);
				//const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
				const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
				tau2 = GlobalC::GridD.getAdjacentTau(ad2);

				dtau = tau2 - tau1;
				distance = dtau.norm() * GlobalC::ucell.lat0;
				rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

				// check if this a adjacent atoms.
				bool is_adj = false;
				if(distance < rcut) is_adj = true;
				else if(distance >= rcut)
				{
					for (int ad0=0; ad0 < GlobalC::GridD.getAdjacentNum()+1 ; ++ad0)
					{
						const int T0 = GlobalC::GridD.getType(ad0);
						if( GlobalC::ucell.infoNL.nproj[T0] == 0) continue;
						const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GlobalC::GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						distance1 = dtau1.norm() * GlobalC::ucell.lat0;
						rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

						dtau2 = tau0 - tau2;
						distance2 = dtau2.norm() * GlobalC::ucell.lat0;
						rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(); 

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
					for (int j=0; j<atom1->nw; ++j)
					{
						const int iw1_all = start1 + j;
						const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
						if(mu < 0)continue;
						for (int k=0; k<atom2->nw; ++k)
						{
							const int iw2_all = start2 + k;
							const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
							if(nu < 0)continue;
							
							for (int ad0=0; ad0 < GlobalC::GridD.getAdjacentNum()+1 ; ++ad0)
							{
								const int T0 = GlobalC::GridD.getType(ad0);
								if( GlobalC::ucell.infoNL.nproj[T0] == 0) continue;
								const int I0 = GlobalC::GridD.getNatom(ad0);
								const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
								//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
								tau0 = GlobalC::GridD.getAdjacentTau(ad0);

								dtau1 = tau0 - tau1;
								distance1 = dtau1.norm() * GlobalC::ucell.lat0;
								rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

								dtau2 = tau0 - tau2;
								distance2 = dtau2.norm() * GlobalC::ucell.lat0;
								rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

								double r0[3];
								double r1[3];
								r1[0] = ( tau1.x - tau0.x) ;
								r1[1] = ( tau1.y - tau0.y) ;
								r1[2] = ( tau1.z - tau0.z) ;
								r0[0] = ( tau2.x - tau0.x) ;
								r0[1] = ( tau2.y - tau0.y) ;
								r0[2] = ( tau2.z - tau0.z) ;

								if(distance1 < rcut1 && distance2 < rcut2)
								{
									//const Atom* atom0 = &GlobalC::ucell.atoms[T0];
									double nlm[3]={0,0,0};

									GlobalC::UOT.snap_psibeta(
											GlobalC::ORB,
											GlobalC::ucell.infoNL,
											nlm, 1,
											tau2,
											T2,
											atom2->iw2l[ k ], // L2
											atom2->iw2m[ k ], // m2
											atom2->iw2n[ k ], // n2
											tau1,
											T1,
											atom1->iw2l[ j ], // L1
											atom1->iw2m[ j ], // m1
											atom1->iw2n[ j ], // N1
											tau0, T0, GlobalC::ucell.atoms[T0].dion, GlobalV::NSPIN,
											GlobalC::ucell.atoms[T0].d_so,
											GlobalC::ucell.atoms[T0].non_zero_count_soc[0], // index stands for spin
											GlobalC::ucell.atoms[T0].index1_soc[0],
											GlobalC::ucell.atoms[T0].index2_soc[0],
											GlobalC::ucell.atoms[T0].nproj_soc
											); // mohan  add 2021-05-07

									double nlm1[3]={0,0,0};
									if(isstress)
									{
										GlobalC::UOT.snap_psibeta(
											GlobalC::ORB,
											GlobalC::ucell.infoNL,
											nlm1, 1,
											tau1,
											T1,
											atom1->iw2l[ j ], // L1
											atom1->iw2m[ j ], // m1
											atom1->iw2n[ j ], // n1
											tau2,
											T2,
											atom2->iw2l[ k ], // L2
											atom2->iw2m[ k ], // m2
											atom2->iw2n[ k ], // N2
											tau0, T0, GlobalC::ucell.atoms[T0].dion, GlobalV::NSPIN,
											GlobalC::ucell.atoms[T0].d_so,
											GlobalC::ucell.atoms[T0].non_zero_count_soc[0], // index stands for spin
											GlobalC::ucell.atoms[T0].index1_soc[0],
											GlobalC::ucell.atoms[T0].index2_soc[0],
											GlobalC::ucell.atoms[T0].nproj_soc
											); // mohan  add 2021-05-07
									}
									/// only one projector for each atom force, but another projector for stress
									for(int is=0; is<GlobalV::NSPIN; ++is)
									{
										double dm2d2 = 2.0 * dm2d[is][iir];
										for(int jpol=0;jpol<3;jpol++)
										{
											if(isforce)
											{
												fvnl_dbeta(iat0, jpol) -= dm2d2 * nlm[jpol];
											}
											if(isstress) 
											{
												for(int ipol=0;ipol<3;ipol++)
												{
													svnl_dbeta(jpol, ipol) += dm2d[is][iir] * 
													(nlm[jpol] * r1[ipol] + nlm1[jpol] * r0[ipol]);
												}
											}
										}
									}

								}// distance
							}// ad0

							++iir;
						}// k
					}// j
				}// distance
			}// ad2
		}// I1
	}// T1

	assert( iir == GlobalC::LNNR.nnr );

	if(isstress)
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				svnl_dbeta(i,j) *=  GlobalC::ucell.lat0 / GlobalC::ucell.omega;
			}
		}
	}

	ModuleBase::timer::tick("Force_LCAO_k","cal_fvnl_dbeta_k");
	return;
}

typedef std::tuple<int,int,int,int> key_tuple;

// must consider three-center H matrix.
void Force_LCAO_k::cal_fvnl_dbeta_k_new(
	double** dm2d, 
	const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& fvnl_dbeta, 
	ModuleBase::matrix& svnl_dbeta)
{
	ModuleBase::TITLE("Force_LCAO_k","cal_fvnl_dbeta_k_new");
	ModuleBase::timer::tick("Force_LCAO_k","cal_fvnl_dbeta_k_new");

	
	for(int iat=0;iat<GlobalC::ucell.nat;iat++)
	{

		const int it = GlobalC::ucell.iat2it[iat];
		const int ia = GlobalC::ucell.iat2ia[iat];

		//Step 1 : generate <psi|beta>
		//type of atom; distance; atomic basis; projectors

		std::map<key_tuple,std::unordered_map<int,std::vector<std::vector<double>>>> nlm_tot;
		nlm_tot.clear();

		const double Rcut_Beta = GlobalC::ucell.infoNL.Beta[it].get_rcut_max();
		const ModuleBase::Vector3<double> tau = GlobalC::ucell.atoms[it].tau[ia];
        GlobalC::GridD.Find_atom(GlobalC::ucell, tau ,it, ia);

		for (int ad=0; ad<GlobalC::GridD.getAdjacentNum()+1 ; ++ad)
		{
			const int T1 = GlobalC::GridD.getType(ad);
			const int I1 = GlobalC::GridD.getNatom(ad);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

			const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad);
			const Atom* atom1 = &GlobalC::ucell.atoms[T1];
			const int nw1_tot = atom1->nw*GlobalV::NPOL;

			const ModuleBase::Vector3<double> dtau = tau1-tau;
			const double dist1 = dtau.norm2() * pow(GlobalC::ucell.lat0,2);
			if (dist1 > pow(Rcut_Beta + Rcut_AO1,2))
			{
				continue;
			}

			std::unordered_map<int,std::vector<std::vector<double>>> nlm_cur;
			nlm_cur.clear();

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
				GlobalC::UOT.snap_psibeta_half(
					GlobalC::ORB,
					GlobalC::ucell.infoNL,
					nlm, tau1, T1,
					atom1->iw2l[ iw1_0 ], // L1
					atom1->iw2m[ iw1_0 ], // m1
					atom1->iw2n[ iw1_0 ], // N1
					tau, it, 1); //R0,T0

				nlm_cur.insert({iw1_all,nlm});
			}//end iw
			const int iat1=GlobalC::ucell.itia2iat(T1, I1);
			const int rx1=GlobalC::GridD.getBox(ad).x;
			const int ry1=GlobalC::GridD.getBox(ad).y;
			const int rz1=GlobalC::GridD.getBox(ad).z;
			key_tuple key_1(iat1,rx1,ry1,rz1);
			nlm_tot[key_1]=nlm_cur;
		}//end ad

		//=======================================================
		//Step2:	
		//calculate sum_(L0,M0) beta<psi_i|beta><beta|psi_j>
		//and accumulate the value to Hloc_fixedR(i,j)
		//=======================================================
		int iir = 0;
		ModuleBase::Vector3<double> tau1;
		ModuleBase::Vector3<double> tau2;
		ModuleBase::Vector3<double> dtau;
		ModuleBase::Vector3<double> tau0;
		ModuleBase::Vector3<double>	dtau1;
		ModuleBase::Vector3<double> dtau2;

		double rcut;
		double distance;

		double rcut1;
		double rcut2;
		double distance1;
		double distance2;

		for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
		{
			const Atom* atom1 = &GlobalC::ucell.atoms[T1];

			for (int I1 =0; I1< atom1->na; ++I1)
			{
				tau1 = atom1->tau[I1];

				GlobalC::GridD.Find_atom(GlobalC::ucell, tau1 ,T1, I1);
				const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
				const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);

				for (int ad2=0; ad2<GlobalC::GridD.getAdjacentNum()+1 ; ++ad2)
				{
					const int T2 = GlobalC::GridD.getType(ad2);
					const Atom* atom2 = &GlobalC::ucell.atoms[T2];
					const int I2 = GlobalC::GridD.getNatom(ad2);
					const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
					tau2 = GlobalC::GridD.getAdjacentTau(ad2);

					dtau = tau2 - tau1;
					distance = dtau.norm() * GlobalC::ucell.lat0;
					rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

					// check if this a adjacent atoms.
					bool is_adj = false;
					if(distance < rcut) is_adj = true;
					else if(distance >= rcut)
					{
						for (int ad0=0; ad0 < GlobalC::GridD.getAdjacentNum()+1 ; ++ad0)
						{
							const int T0 = GlobalC::GridD.getType(ad0);
							if( GlobalC::ucell.infoNL.nproj[T0] == 0) continue;
							const int I0 = GlobalC::GridD.getNatom(ad0);
							//const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
							//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

							tau0 = GlobalC::GridD.getAdjacentTau(ad0);
							dtau1 = tau0 - tau1;
							distance1 = dtau1.norm() * GlobalC::ucell.lat0;
							rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

							dtau2 = tau0 - tau2;
							distance2 = dtau2.norm() * GlobalC::ucell.lat0;
							rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(); 

							if( distance1 < rcut1 && distance2 < rcut2 )
							{
								is_adj = true;
								break;
							}
						}
					}

					if(is_adj)
					{

						const int rx2=GlobalC::GridD.getBox(ad2).x;
						const int ry2=GlobalC::GridD.getBox(ad2).y;
						const int rz2=GlobalC::GridD.getBox(ad2).z;

						// < psi1 | all projectors | psi2 >
						// ----------------------------- enter the nnr increaing zone -------------------------
						for (int j=0; j<atom1->nw; ++j)
						{
							const int iw1_all = start1 + j;
							const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
							if(mu < 0)continue;
							for (int k=0; k<atom2->nw; ++k)
							{
								const int iw2_all = start2 + k;
								const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
								if(nu < 0)continue;
								
								for (int ad0=0; ad0 < GlobalC::GridD.getAdjacentNum()+1 ; ++ad0)
								{
									const int T0 = GlobalC::GridD.getType(ad0);
									if( GlobalC::ucell.infoNL.nproj[T0] == 0) continue;
									const int I0 = GlobalC::GridD.getNatom(ad0);
									if(T0!=it || I0!=ia) continue;
									
									const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
									//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
									tau0 = GlobalC::GridD.getAdjacentTau(ad0);

									dtau1 = tau0 - tau1;
									distance1 = dtau1.norm() * GlobalC::ucell.lat0;
									rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

									dtau2 = tau0 - tau2;
									distance2 = dtau2.norm() * GlobalC::ucell.lat0;
									rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

									double r0[3];
									double r1[3];
									r1[0] = ( tau1.x - tau0.x) ;
									r1[1] = ( tau1.y - tau0.y) ;
									r1[2] = ( tau1.z - tau0.z) ;
									r0[0] = ( tau2.x - tau0.x) ;
									r0[1] = ( tau2.y - tau0.y) ;
									r0[2] = ( tau2.z - tau0.z) ;

									if(distance1 < rcut1 && distance2 < rcut2)
									{
										const int rx0=GlobalC::GridD.getBox(ad0).x;
										const int ry0=GlobalC::GridD.getBox(ad0).y;
										const int rz0=GlobalC::GridD.getBox(ad0).z;
										key_tuple key1(iat1,-rx0,-ry0,-rz0);
										key_tuple key2(iat2,rx2-rx0,ry2-ry0,rz2-rz0);

										//const Atom* atom0 = &GlobalC::ucell.atoms[T0];
										double nlm[3]={0,0,0};
										std::vector<double> nlm_1 = nlm_tot[key2][iw2_all][0];
										std::vector<std::vector<double>> nlm_2;
										nlm_2.resize(3);
										for(int i=0;i<3;i++)
										{
											nlm_2[i] = nlm_tot[key1][iw1_all][i+1];
										}

										assert(nlm_1.size()==nlm_2[0].size());

										const int nproj = GlobalC::ucell.infoNL.nproj[T0];
										int ib = 0;
										for (int nb = 0; nb < nproj; nb++)
										{
											const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
											for(int m=0;m<2*L0+1;m++)
											{
												for(int ir=0;ir<3;ir++)
												{
													nlm[ir] += nlm_2[ir][ib]*nlm_1[ib]*GlobalC::ucell.atoms[T0].dion(nb,nb);
												}
												ib+=1;
											}
										}
										assert(ib==nlm_1.size());

										double nlm1[3]={0,0,0};
										if(isstress)
										{
											std::vector<double> nlm_1 = nlm_tot[key1][iw1_all][0];
											std::vector<std::vector<double>> nlm_2;
											nlm_2.resize(3);
											for(int i=0;i<3;i++)
											{
												nlm_2[i] = nlm_tot[key2][iw2_all][i+1];
											}

											assert(nlm_1.size()==nlm_2[0].size());

											const int nproj = GlobalC::ucell.infoNL.nproj[T0];
											int ib = 0;
											for (int nb = 0; nb < nproj; nb++)
											{
												const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
												for(int m=0;m<2*L0+1;m++)
												{
													for(int ir=0;ir<3;ir++)
													{
														nlm1[ir] += nlm_2[ir][ib]*nlm_1[ib]*GlobalC::ucell.atoms[T0].dion(nb,nb);
													}
													ib+=1;
												}
											}
											assert(ib==nlm_1.size());
										}
										/// only one projector for each atom force, but another projector for stress
										for(int is=0; is<GlobalV::NSPIN; ++is)
										{
											double dm2d2 = 2.0 * dm2d[is][iir];
											for(int jpol=0;jpol<3;jpol++)
											{
												if(isforce)
												{
													fvnl_dbeta(iat0, jpol) -= dm2d2 * nlm[jpol];
												}
												if(isstress) 
												{
													for(int ipol=0;ipol<3;ipol++)
													{
														svnl_dbeta(jpol, ipol) += dm2d[is][iir] * 
														(nlm[jpol] * r1[ipol] + nlm1[jpol] * r0[ipol]);
													}
												}
											}
										}

									}// distance
								}// ad0

								++iir;
							}// k
						}// j
					}// distance
				}// ad2
			}// I1
		}// T1

		assert( iir == GlobalC::LNNR.nnr );
	}//iat

	if(isstress)
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				svnl_dbeta(i,j) *=  GlobalC::ucell.lat0 / GlobalC::ucell.omega;
			}
		}
	}

	ModuleBase::timer::tick("Force_LCAO_k","cal_fvnl_dbeta_k_new");
	return;
}

// calculate the force due to < phi | Vlocal | dphi >
void Force_LCAO_k::cal_fvl_dphi_k(
	double** dm2d, 
	const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& fvl_dphi, 
	ModuleBase::matrix& svl_dphi)
{
	ModuleBase::TITLE("Force_LCAO_k","cal_fvl_dphi_k");
	ModuleBase::timer::tick("Force_LCAO_k","cal_fvl_dphi_k");

	if(!isforce&&!isstress) 
	{
		ModuleBase::timer::tick("Force_LCAO_k","cal_fvl_dphi_k");
		return;
	}
	assert(GlobalC::LM.DHloc_fixedR_x!=NULL);
	assert(GlobalC::LM.DHloc_fixedR_y!=NULL);
	assert(GlobalC::LM.DHloc_fixedR_z!=NULL);

	int istep = 1;

	// if Vna potential is not used.
	GlobalC::pot.init_pot(istep, GlobalC::pw.strucFac);


	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		GlobalV::CURRENT_SPIN = is;
//		ZEROS (GlobalC::LM.DHloc_fixedR_x, GlobalC::LNNR.nnr);
//		ZEROS (GlobalC::LM.DHloc_fixedR_y, GlobalC::LNNR.nnr);
//		ZEROS (GlobalC::LM.DHloc_fixedR_z, GlobalC::LNNR.nnr);
//		std::cout << " CURRENT_SPIN=" << GlobalV::CURRENT_SPIN << std::endl;

		for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
		}

		//--------------------------------
		// Grid integration here.
		//--------------------------------
		// fvl_dphi can not be set to zero here if Vna is used
		if(isstress||isforce) 
		{
			GlobalC::UHM.GK.svl_k_RealSpace(isforce, isstress, fvl_dphi,svl_dphi,GlobalC::pot.vr_eff1);
		}
		else if(isforce) 
		{
			GlobalC::UHM.GK.fvl_k_RealSpace(fvl_dphi,GlobalC::pot.vr_eff1);
		}
	}

	
	if(isstress){
		for(int ipol=0;ipol<3;ipol++){
			for(int jpol=0;jpol<3;jpol++){
				if(ipol < jpol) svl_dphi(jpol, ipol) = svl_dphi(ipol, jpol);
				svl_dphi(ipol, jpol) /= GlobalC::ucell.omega;
			}
		}
	}

	ModuleBase::timer::tick("Force_LCAO_k","cal_fvl_dphi_k");
	return;
}

void Force_LCAO_k::calFvnlDbeta(
	double** dm2d, 
	const bool &isforce, 
	const bool &isstress, 
	ModuleBase::matrix& fvnl_dbeta, 
	ModuleBase::matrix& svnl_dbeta,
	const int &vnl_method)
{
	ModuleBase::TITLE("Force_LCAO_k", "calFvnlDbeta");
	if(GlobalV::NSPIN==4 || vnl_method == 0)
	{
		this->cal_fvnl_dbeta_k(dm2d, isforce, isstress, fvnl_dbeta, svnl_dbeta);
	}
	else if(vnl_method == 1)
	{
		this->cal_fvnl_dbeta_k_new(dm2d, isforce, isstress, fvnl_dbeta, svnl_dbeta);
	}
	else 
    {
        ModuleBase::WARNING_QUIT("Force_LCAO_k","This method has not been implemented");
    }
}
