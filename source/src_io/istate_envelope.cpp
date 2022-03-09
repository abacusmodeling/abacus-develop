#include "istate_envelope.h"
#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"

IState_Envelope::IState_Envelope()
{}

IState_Envelope::~IState_Envelope()
{}


void IState_Envelope::begin(Local_Orbital_wfc &lowf, Gint_Gamma &gg)
{
	ModuleBase::TITLE("IState_Envelope","begin");

	std::cout << " perform |psi(band, r)| for selected bands." << std::endl;

	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
		ModuleBase::WARNING_QUIT("IState_Envelope::begin","Only available for GAMMA_ONLY_LOCAL now.");
	}

	// (1) 
	// (1.1) allocate the space for GlobalC::LOWF.WFC_GAMMA

	// (1.2) read in LOWF_GAMMA.dat

	// mohan update 2011-03-21
	// if ucell is odd, it's correct,
	// if ucell is even, it's also correct.
	// +1.0e-8 in case like (2.999999999+1)/2
	int fermi_band = static_cast<int>( (GlobalC::CHR.nelec+1)/2 + 1.0e-8 ) ;
	int bands_below = GlobalV::NBANDS_ISTATE;
	int bands_above = GlobalV::NBANDS_ISTATE;

	std::cout << " number of electrons = " << GlobalC::CHR.nelec << std::endl;
	std::cout << " number of occupied bands = " << fermi_band << std::endl;
	std::cout << " plot band decomposed charge density below fermi surface with " 
	<< bands_below << " bands." << std::endl;
	
	std::cout << " plot band decomposed charge density above fermi surface with " 
	<< bands_above << " bands." << std::endl;
	
	// (2) cicle:
	
	// (2.1) calculate the selected density matrix
	// from wave functions.

	// (2.2) carry out the grid integration to
	// get the charge density.

	// (2.3) output the charge density in .cub format.
	this->bands_picked = new bool[GlobalV::NBANDS];
	ModuleBase::GlobalFunc::ZEROS(bands_picked, GlobalV::NBANDS);
	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		if( ib >= fermi_band - bands_below ) 
		{
			if( ib < fermi_band + bands_above)
			{
				bands_picked[ib] = true;
			}
		}
	}

	//allocate grid wavefunction for gamma_only
	std::vector<double**> wfc_gamma_grid(GlobalV::NSPIN);
	for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        wfc_gamma_grid[is] = new double* [GlobalV::NBANDS];
        for (int ib = 0;ib < GlobalV::NBANDS; ++ib)
            wfc_gamma_grid[is][ib] = new double[GlobalC::GridT.lgd];
    }

    const Parallel_Orbitals* pv = lowf.ParaV;
    
    //calculate maxnloc for bcasting 2d-wfc
    int nprocs, myid;
    MPI_Comm_size(pv->comm_2D, &nprocs);
    MPI_Comm_rank(pv->comm_2D, &myid);
    
    long maxnloc; // maximum number of elements in local matrix
	MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
    MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
    const int inc = 1;
    int naroc[2]; // maximum number of row or column
    double* work = new double[maxnloc]; // work/buffer matrix 
    
    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
	{
		if(bands_picked[ib])
		{
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{
				std::cout << " Perform envelope function for band " << ib+1 << std::endl;
				ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[is],GlobalC::pw.nrxx);	


				//---------------------------------------------------------
				// GlobalC::LOWF.WFC_GAMMA has been replaced by wfc_dm_2d.cpp 
				// we need to fix this function in near future.
				// -- mohan add 2021-02-09
				//---------------------------------------------------------
				//ModuleBase::WARNING_QUIT("IState_Charge::idmatrix","need to update LOWF.WFC_GAMMA");
				
				//convert 2d `wfc_gamma` to grid `wfc_gamma_grid`
				int info;
				for(int iprow=0; iprow<pv->dim0; ++iprow)
				{
					for(int ipcol=0; ipcol<pv->dim1; ++ipcol)
					{
						const int coord[2]={iprow, ipcol};
						int src_rank;
						MPI_Cart_rank(pv->comm_2D, coord, &src_rank);
						if(myid==src_rank)
						{
							BlasConnector::copy(pv->nloc_wfc, lowf.wfc_gamma[is].c, inc, work, inc);
							naroc[0]=pv->nrow;
							naroc[1]=pv->ncol_bands;
						}
						info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, pv->comm_2D);
						info=MPI_Bcast(work, maxnloc, MPI_DOUBLE, src_rank, pv->comm_2D);

                        info=lowf.q2WFC(myid, naroc, pv->nb,
                            pv->dim0, pv->dim1, iprow, ipcol, pv->loc_size,
                            work, wfc_gamma_grid[is]);
					}//loop ipcol
				}//loop iprow

				gg.cal_env( wfc_gamma_grid[is][ib], GlobalC::CHR.rho[is] );


				GlobalC::CHR.save_rho_before_sum_band(); //xiaohui add 2014-12-09
				std::stringstream ss;
				ss << GlobalV::global_out_dir << "BAND" << ib + 1 << "_ENV" << is+1 << "_CHG";
				// 0 means definitely output charge density.
				bool for_plot = true;
				GlobalC::CHR.write_rho(GlobalC::CHR.rho_save[is], is, 0, ss.str(), 3, for_plot );
			}
		}
	}

    delete[] work;
    delete[] bands_picked;
    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        for (int ib = 0;ib < GlobalV::NBANDS; ++ib)
            delete[] wfc_gamma_grid[is][ib];
        delete[] wfc_gamma_grid[is];
    }
    return;
}


