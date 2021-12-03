#include "ORB_control.h"
#include "ORB_gen_tables.h"
//#include "build_st_pw.h"

ORB_control::ORB_control()
{}

ORB_control::~ORB_control()
{}

void ORB_control::read_orb_first(
	std::ofstream &ofs_in, 
	LCAO_Orbitals &orb,
	const int &ntype, // mohan add 2021-04-26
	const int &lmax, // mohan add 2021-04-26 
	const double &lcao_ecut_in, // mohan add 2021-04-16
	const double &lcao_dk_in, // mohan add 2021-04-16
	const double &lcao_dr_in, // mohan add 2021-04-16
	const double &lcao_rmax_in, // mohan add 2021-04-16
	const int &out_descriptor,
	const int &out_r_matrix,
	const bool &force_flag, // mohan add 2021-05-07
	const int &my_rank // mohan add 2021-04-26
) 
{
    ModuleBase::TITLE("ORB_control","read_orb_first");
	ModuleBase::timer::tick("ORB_control","read_orb_first");
    
	/////////////////////////////////////////////////////////////////
	/// (1) FUNCTION : use 'info' to generate 'Numerical Orbital'
	///
	/// (1) RESULT : We have 'Numerical Orbital' for calculate S-table and T-table.
	/////////////////////////////////////////////////////////////////

	// mohan add 2021-04-16
	assert(ntype>0);
	assert(lmax>=0);
	assert(lcao_ecut_in>0.0);
	assert(lcao_dk_in>0.0);
	assert(lcao_dr_in>0.0);
	assert(lcao_rmax_in>0.0);

	// mohan add 2021-04-16
	orb.ecutwfc = lcao_ecut_in;
	orb.dk = lcao_dk_in;
	orb.dR = lcao_dr_in;
	orb.Rmax = lcao_rmax_in;
	
    orb.Read_Orbitals(
		ofs_in,
		ntype, 
		lmax, 
		out_descriptor, 
		out_r_matrix, 
		force_flag,
		my_rank);

	ModuleBase::timer::tick("ORB_control","read_orb_first");
	return;
}

void ORB_control::set_orb_tables(
	std::ofstream &ofs_in,
	ORB_gen_tables &OGT, 
	LCAO_Orbitals &orb,
	const double &lat0,
	const int &out_descriptor,
	const int &Lmax_exx,
	const int &nprojmax, 
	const int* nproj,
	const Numerical_Nonlocal* beta_) 
{
    ModuleBase::TITLE("ORB_control","set_orb_tables");
	ModuleBase::timer::tick("ORB_control","set_orb_tables");


#ifdef __NORMAL

#else
	if(GlobalV::CALCULATION=="test")
	{
		ModuleBase::timer::tick("ORB_control","set_orb_tables");
		return;
	}
#endif


    ///////////////////////////////////////////////////////////////////
    /// (2) FUNCTION : Generate Gaunt_Coefficients and S-table using OGT.init
	/// 	   Must have 'Numerical Orbital' infomation
	///
	/// (2) RESULT : we have tabulated S table for use.
    ///////////////////////////////////////////////////////////////////
    const int job0 = 3;
    /// job0 :
    /// 1. generate overlap table
    /// 2. generate kinetic table
    /// 3. generate overlap & kinetic table
    OGT.gen_tables(ofs_in, job0, orb, Lmax_exx, out_descriptor, nprojmax, nproj, beta_);
    // init lat0, in order to interpolated value from this table.

	assert(lat0>0.0);
    OGT.set_unit(lat0);

	ModuleBase::timer::tick("ORB_control","set_orb_tables");
    return;
}

void ORB_control::clear_after_ions(
	ORB_gen_tables &OGT, 
	LCAO_Orbitals &orb,
	const int &out_descriptor,
	const int* nproj_)
{
    ModuleBase::TITLE("ORB_control","clear_after_ions");
    OGT.MOT.Destroy_Table(orb);
    OGT.tbeta.Destroy_Table_Beta(orb.get_ntype(), orb.Phi, nproj_);
    
	//caoyu add 2021-03-18
    if (out_descriptor>0) 
	{
        OGT.talpha.Destroy_Table_Alpha(orb);
    }
    return;
}
