#include "istate_charge.h"
#include "../src_pw/global.h"
#include "../src_pw/tools.h"
#include "../module_base/scalapack_connector.h"

IState_Charge::IState_Charge(){}

IState_Charge::~IState_Charge(){}


void IState_Charge::begin(void)
{
	ModuleBase::TITLE("IState_Charge","begin");

	std::cout << " Perform |psi(i)|^2 for selected bands." << std::endl;

	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
		ModuleBase::WARNING_QUIT("IState_Charge::begin","Only available for GlobalV::GAMMA_ONLY_LOCAL now.");
	}

	int mode = 0;
	if(GlobalV::NBANDS_ISTATE > 0) mode = 1;
	else mode = 2;

	int fermi_band=0;
	int bands_below=0;
	int bands_above=0;

	// (2) cicle:
	// (2.1) calculate the selected density matrix
	// from wave functions.
	// (2.2) carry out the grid integration to
	// get the charge density.
	this->bands_picked = new int[GlobalV::NBANDS];
	ModuleBase::GlobalFunc::ZEROS(bands_picked, GlobalV::NBANDS);

	// (1) 
	// (1.1) allocate the space for GlobalC::LOWF.WFC_GAMMA

	// (1.2) read in LOWF_GAMMA.dat
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"GlobalC::LOWF.allocate_flag",GlobalC::LOWF.get_allocate_flag());	
	std::cout << " number of electrons = " << GlobalC::CHR.nelec << std::endl;

	// mohan update 2011-03-21
	// if ucell is odd, it's correct,
	// if ucell is even, it's also correct.
	// +1.0e-8 in case like (2.999999999+1)/2
	fermi_band = static_cast<int>( (GlobalC::CHR.nelec+1)/2 + 1.0e-8 ) ;
	std::cout << " number of occupied bands = " << fermi_band << std::endl;

	if(mode == 1)
	{
		bands_below = GlobalV::NBANDS_ISTATE;
		bands_above = GlobalV::NBANDS_ISTATE;

		std::cout << " plot band decomposed charge density below fermi surface with " 
			<< bands_below << " bands." << std::endl;

		std::cout << " plot band decomposed charge density above fermi surface with " 
			<< bands_above << " bands." << std::endl;

		for(int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			if( ib >= fermi_band - bands_below ) 
			{
				if( ib < fermi_band + bands_above)
				{
					bands_picked[ib] = 1;
				}
			}
		}
	}
	else if(mode == 2)
	{
		bool stop = false;
		std::stringstream ss;
		ss << GlobalV::global_out_dir << "istate.info";
		std::cout << " Open the file : " << ss.str() << std::endl; 
		if(GlobalV::MY_RANK==0)
		{
			std::ifstream ifs(ss.str().c_str());
			if(!ifs)
			{
				stop = true;
			}
			else
			{
				//int band_index;
				for(int ib=0; ib<GlobalV::NBANDS; ++ib)
				{
					ModuleBase::GlobalFunc::READ_VALUE(ifs, bands_picked[ib]);
				}
			}
		}

#ifdef __MPI
		Parallel_Common::bcast_bool(stop);
		Parallel_Common::bcast_int(bands_picked, GlobalV::NBANDS);
#endif
		if(stop)
		{
			GlobalV::ofs_warning << " Can't find the file : " << ss.str() << std::endl;
			ModuleBase::WARNING_QUIT("IState_Charge::begin","can't find the istate file.");
		}
	}

	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		if(bands_picked[ib])
		{
			std::cout << " Perform band decomposed charge density for band " << ib+1 << std::endl;
			// (1)
			// This has been done once in LOOP_ions.
			// but here we need to done for each band.
			//GlobalC::LOC.allocate_gamma(GridT);	
			
			// (2) calculate the density matrix for a partuclar 
			// band, whenever it is occupied or not.
			this->idmatrix(ib);

			// (3) zero out of charge density array. 
			for(int is=0; is<GlobalV::NSPIN; is++)
			{
				ModuleBase::GlobalFunc::ZEROS( GlobalC::CHR.rho[is], GlobalC::pw.nrxx );
			}
			
			// (4) calculate charge density for a particular 
			// band.
   			GlobalC::UHM.GG.cal_rho(GlobalC::LOC.DM);
			GlobalC::CHR.save_rho_before_sum_band(); //xiaohui add 2014-12-09
			std::stringstream ss;
			std::stringstream ss1;
			ss << GlobalV::global_out_dir << "BAND" << ib + 1 << "_CHG";
			// 0 means definitely output charge density.
			for(int is=0; is<GlobalV::NSPIN; is++)
			{
				ss1 << GlobalV::global_out_dir << "BAND" << ib + 1 << "_SPIN" << is << "_CHG.cube";
				bool for_plot = true;
				GlobalC::CHR.write_rho(GlobalC::CHR.rho_save[is], is, 0, ss.str(), 3, for_plot );
				GlobalC::CHR.write_rho_cube(GlobalC::CHR.rho_save[is], is, ss1.str(), 3);
			}
		}
	}

	delete[] bands_picked;
	return;
}


void IState_Charge::idmatrix(const int &ib)
{
	ModuleBase::TITLE("IState_Charge","idmatrix");
/*		
	for(int is=0; is<NSPIN; is++)
	{
		for (int i=0; i<GlobalV::NLOCAL; i++)
		{
			const int mu_local = GlobalC::GridT.trace_lo[i];
			if ( mu_local >= 0)
			{
				// set a pointer.
				//double *alpha = GlobalC::LOC.DM[is][mu_local];
				for (int j=i; j<GlobalV::NLOCAL; j++)
				{
					const int nu_local = GlobalC::GridT.trace_lo[j];
					if ( nu_local >= 0)
					{
						//---------------------------------------------------------
						// GlobalC::LOWF.WFC_GAMMA has been replaced by wfc_dm_2d.cpp 
						// we need to fix this function in near future.
						// -- mohan add 2021-02-09
						//---------------------------------------------------------
						ModuleBase::WARNING_QUIT("IState_Charge::idmatrix","need to update GlobalC::LOWF.WFC_GAMMA");
						// 2 stands for degeneracy.
						//alpha[nu_local] += 2.0 * GlobalC::LOWF.WFC_GAMMA[is][ib][mu_local] * GlobalC::LOWF.WFC_GAMMA[is][ib][nu_local];
					}
				}
			}
		}
	}
	return;
*/

		assert(GlobalC::wf.wg.nr==GlobalV::NSPIN);
		for(int is=0; is!=GlobalV::NSPIN; ++is)
		{
			std::vector<double> wg_local(GlobalC::ParaO.ncol,0.0);
			const int ib_local = GlobalC::ParaO.trace_loc_col[ib];

			int fermi_band=0;
			fermi_band = static_cast<int>( (GlobalC::CHR.nelec+1)/2 + 1.0e-8 ) ;

			if(ib_local>=0)
			{
				if(ib<fermi_band)
				{
					wg_local[ib_local] = GlobalC::wf.wg(is,ib);
				}
				else
				{
					wg_local[ib_local] = GlobalC::wf.wg(is,fermi_band-1);
				}//unoccupied bands, use occupation of homo
			}
		
			// wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
			ModuleBase::matrix wg_wfc(GlobalC::LOC.wfc_dm_2d.wfc_gamma[is]);
	
			for(int ir=0; ir!=wg_wfc.nr; ++ir)
			{
				LapackConnector::scal( wg_wfc.nc, wg_local[ir], wg_wfc.c+ir*wg_wfc.nc, 1 );
			}

			// C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
			const double one_float=1.0, zero_float=0.0;
			const int one_int=1;
			const char N_char='N', T_char='T';
			GlobalC::LOC.wfc_dm_2d.dm_gamma[is].create( wg_wfc.nr, wg_wfc.nc );

			pdgemm_(
				&N_char, &T_char,
				&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalC::wf.wg.nc,
				&one_float,
				wg_wfc.c, &one_int, &one_int, GlobalC::ParaO.desc,
				GlobalC::LOC.wfc_dm_2d.wfc_gamma[is].c, &one_int, &one_int, GlobalC::ParaO.desc,
				&zero_float,
				GlobalC::LOC.wfc_dm_2d.dm_gamma[is].c, &one_int, &one_int, GlobalC::ParaO.desc);
		}

		std::cout << " finished calc dm_2d : " << std::endl;

		GlobalC::LOC.cal_dk_gamma_from_2D_pub();
		
		std::cout << " finished convert : " << std::endl;

}

