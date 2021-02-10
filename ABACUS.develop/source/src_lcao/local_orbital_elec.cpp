#include "local_orbital_elec.h"
#include "diago_lcao_matrix.h"
#include "src_pw/global.h"
#include "src_pw/symmetry_rho.h"
#include "input_update.h"
#include "src_pw/chi0_hilbert.h"
#include "evolve_lcao_matrix.h"
#include "dftu.h"
//
#include "../src_global/sltk_atom_arrange.h"
#include "src_lcao/lcao_nnr.h"
#include "../src_io/istate_charge.h"
#include "../src_io/istate_envelope.h"
#include "ELEC_scf.h"
#include "ELEC_nscf.h"
#include "ELEC_cbands_gamma.h"
#include "ELEC_cbands_k.h"
#include "ELEC_evolve.h"
//
#include "src_lcao/exx_abfs.h"
#include "src_lcao/exx_opt_orb.h"


void Local_Orbital_Elec::solve_elec_stru(const int &istep)
{
    TITLE("Local_Orbital_Elec","solve_elec_stru"); 

	// prepare HS matrices, prepare grid integral
	this->set_matrix_grid_index();
	// density matrix extrapolation and prepare S,T,VNL matrices 
	this->before_solver(istep);
	// do self-interaction calculations / nscf/ tddft, etc. 
	this->solver(istep);

	return;
}


void Local_Orbital_Elec::set_matrix_grid_index(void)
{
    TITLE("Local_Orbital_Elec","set_matrix_grid_index"); 

	// (1) Find adjacent atoms for each atom.
	atom_arrange::set_sr_NL();
	atom_arrange::search( SEARCH_RADIUS );
	//DONE(ofs_running,"SEARCH ADJACENT ATOMS");

	// (3) Periodic condition search for each grid.
	// here 0 means GridT is not used for Vna.
	// mohan add Vna 2012-06-13
	// because it must be consistent with gtf in lcao_vna.
	GridT.set_pbc_grid(
			pw.ncx, pw.ncy, pw.ncz,
			pw.bx, pw.by, pw.bz,
			pw.nbx, pw.nby, pw.nbz,
			pw.nbxx, pw.nbzp_start, pw.nbzp);

	// (2) If k point is used here, allocate HlocR after atom_arrange.
	if(!GAMMA_ONLY_LOCAL)
	{
		// For each atom, calculate the adjacent atoms in different cells
		// and allocate the space for H(R) and S(R).
		LNNR.cal_nnr();
		LM.allocate_HS_R(LNNR.nnr);

		// need to first calculae lgd.
		// using GridT.init.
		LNNR.cal_nnrg(GridT);
	}

	return;
}


void Local_Orbital_Elec::before_solver(const int &istep)
{
    TITLE("Local_Orbital_Elec","before_solver"); 

	// set the augmented orbitals index.
	// after ParaO and GridT, 
	// this information is used to calculate
	// the force.
	LOWF.set_trace_aug(GridT);

	// init density kernel and wave functions.
	LOC.allocate_dm_wfc(GridT);

	//======================================
	// do the charge extrapolation before the density matrix is regenerated.
	// mohan add 2011-04-08
	// because once atoms are moving out of this processor,
	// the density matrix will not map the new atomic configuration,
	//======================================
	// THIS IS A BUG, BECAUSE THE INDEX GridT.trace_lo
	// HAS BEEN REGENERATED, SO WE NEED TO
	// REALLOCATE DENSITY MATRIX FIRST, THEN READ IN DENSITY MATRIX,
	// AND USE DENSITY MATRIX TO DO RHO CALCULATION.-- mohan 2013-03-31
	//======================================
	if(pot.extra_pot=="dm" && istep>1)//xiaohui modify 2015-02-01
	{
		for(int is=0; is<NSPIN; is++)
		{
			ZEROS(CHR.rho[is], pw.nrxx);
			stringstream ssd;
			ssd << global_out_dir << "SPIN" << is + 1 << "_DM" ;
			// reading density matrix,
			LOC.read_dm(is, ssd.str() );
		}

		// calculate the charge density
		if(GAMMA_ONLY_LOCAL)
		{
			UHM.GG.cal_rho();
		}
		else
		{
			UHM.GK.calculate_charge();
		}

		// renormalize the charge density
		CHR.renormalize_rho();

		// initialize the potential
		pot.init_pot( istep-1 );
	}


	// (9) compute S, T, Vnl, Vna matrix.
	UHM.set_ion();

	return;
}

void Local_Orbital_Elec::solver(const int &istep)
{

    TITLE("Local_Orbital_Elec","solver"); 

	// Peize Lin add 2014-04-04, update 2019-04-26
	if(vdwd2.vdwD2)
	{
		vdwd2.energy();
	}
	// jiyy add 2019-05-18
	else if(vdwd3.vdwD3)
	{
		vdwd3.energy();
	}

	// self consistent calculations for electronic ground state
	if (CALCULATION=="scf" || CALCULATION=="md"
			|| CALCULATION=="relax" || CALCULATION=="cell-relax") //pengfei 2014-10-13
	{
		//Peize Lin add 2016-12-03
		switch(exx_lcao.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				exx_lcao.cal_exx_ions();
				break;
			case Exx_Global::Hybrid_Type::No:
			case Exx_Global::Hybrid_Type::Generate_Matrix:
				break;
			default:
				throw invalid_argument(TO_STRING(__FILE__)+TO_STRING(__LINE__));
		}

		// No exx
		if( Exx_Global::Hybrid_Type::No==exx_global.info.hybrid_type  )
		{
			ELEC_scf es;
			es.scf(istep-1);
		}
		else if( Exx_Global::Hybrid_Type::Generate_Matrix == exx_global.info.hybrid_type )
		{
			Exx_Opt_Orb exx_opt_orb;
			exx_opt_orb.generate_matrix();
		}
		else    // Peize Lin add 2016-12-03
		{
			ELEC_scf es;
			es.scf(istep-1);
			if( exx_global.info.separate_loop )
			{
				for( size_t hybrid_step=0; hybrid_step!=exx_global.info.hybrid_step; ++hybrid_step )
				{
					exx_global.info.set_xcfunc(xcf);
					exx_lcao.cal_exx_elec();
					
					ELEC_scf es;
					es.scf(istep-1);
					if(ELEC_scf::iter==1)     // exx converge
					{
						break;
					}
				}
			}
			else
			{
				exx_global.info.set_xcfunc(xcf);

				ELEC_scf es;
				es.scf(istep-1);
				
			}
		}
	}
	else if (CALCULATION=="nscf")
	{
		ELEC_nscf::nscf(UHM);
	}
	else if (CALCULATION=="istate")
	{
		IState_Charge ISC;
		ISC.begin();
	}
	else if (CALCULATION=="ienvelope")
	{
		IState_Envelope IEP;
		IEP.begin();
	}
	else
	{
		WARNING_QUIT("Local_Orbital_Ions::opt_ions","What's the CALCULATION.");
	}

	return;
}

