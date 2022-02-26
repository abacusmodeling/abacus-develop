#include "ORB_control.h"
#include "ORB_gen_tables.h"
#include "../module_base/timer.h"
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


void ORB_control::setup_2d_division(void)
{
    ModuleBase::TITLE("ORB_control","setup_2d_division");
    GlobalV::ofs_running << "\n SETUP THE DIVISION OF H/S MATRIX" << std::endl;
    
    // (1) calculate nrow, ncol, nloc.
    if (GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="hpseps" || GlobalV::KS_SOLVER=="scalpack" 
        || GlobalV::KS_SOLVER=="selinv" || GlobalV::KS_SOLVER=="scalapack_gvx")
    {
        GlobalV::ofs_running << " divide the H&S matrix using 2D block algorithms." << std::endl;
#ifdef __MPI
        // storage form of H and S matrices on each processor
        // is determined in 'divide_HS_2d' subroutine
        this->divide_HS_2d(DIAG_WORLD);
#else
        ModuleBase::WARNING_QUIT("LCAO_Matrix::init","diago method is not ready.");
#endif
	}
	else
	{
		// the full matrix
		this->ParaV.nloc = GlobalV::NLOCAL * GlobalV::NLOCAL;
	}

	// (2) set the trace, then we can calculate the nnr.
	// for 2d: calculate po.nloc first, then trace_loc_row and trace_loc_col
	// for O(N): calculate the three together.
	this->set_trace();
}


#ifdef __MPI
#include "src_pdiag/pdgseps.h"
#include "src_pdiag/pzgseps.h"
void ORB_control::readin(
    const std::string& fa,
	const std::string &fb,
	const int &nlocal_tot,
	double *eigen,
	double *eigvr)
{
    ModuleBase::TITLE("Pdiag_Double","readin");

    Parallel_Orbitals* pv = &this->ParaV;
    int coord[2];
    int dim[2];
    int period[2];
    int i,j,tmp1,tmp2;
    int k,loc_size,loc_pos;
    double time1,time2;

    MPI_Comm comm=DIAG_HPSEPS_WORLD,comm_2D,comm_col,comm_row,newcomm;

    dim[0]=(int)sqrt((double)GlobalV::DSIZE);

    while (GlobalV::DSIZE%dim[0]!=0)
    {
        dim[0]=dim[0]-1;
    }
    dim[1]=GlobalV::DSIZE/dim[0];

    // call mpi_creat_cart
    this->mpi_creat_cart(&pv->comm_2D,dim[0],dim[1]);

    // call mat_2d
    this->mat_2d(pv->comm_2D, nlocal_tot,nlocal_tot,pv->nb,pv->MatrixInfo);

    pv->loc_size=nlocal_tot/GlobalV::DSIZE;
    if (GlobalV::DRANK<nlocal_tot%GlobalV::DSIZE) loc_size=loc_size+1;

    GlobalV::ofs_running << " loc_size = " << loc_size;

    /*Distribute the matrix*/
    const long nloc = pv->MatrixInfo.col_num * pv->MatrixInfo.row_num;

    double *A = new double[nloc];
    double *B = new double[nloc];
    double *Z = new double[loc_size*nlocal_tot];
    ModuleBase::GlobalFunc::ZEROS(A, nloc);
    ModuleBase::GlobalFunc::ZEROS(B, nloc);
    ModuleBase::GlobalFunc::ZEROS(Z, loc_size * nlocal_tot);

    GlobalV::ofs_running << "\n Data distribution of H." << std::endl;
    this->data_distribution(pv->comm_2D,fa,nlocal_tot,pv->nb,A,pv->MatrixInfo);
    GlobalV::ofs_running << "\n Data distribution of S." << std::endl;
    this->data_distribution(pv->comm_2D,fb,nlocal_tot,pv->nb,B,pv->MatrixInfo);

    time1=MPI_Wtime();
    // call pdgseps
    char uplo = 'U';
    pdgseps(pv->comm_2D,nlocal_tot,pv->nb,A,B,Z,eigen,pv->MatrixInfo,uplo,pv->loc_size,loc_pos);
    time2=MPI_Wtime();
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"time1",time1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"time2",time2);

    //this->gath_eig(comm,n,eigvr,Z);

    GlobalV::ofs_running << "\n " << std::setw(6) << "Band" << std::setw(25) << "Ry" << std::setw(25) << " eV" << std::endl;
    for(int i=0; i<nlocal_tot; i++)
    {
        GlobalV::ofs_running << " " << std::setw(6) << i << std::setw(25) << eigen[i] << std::setw(25)<< eigen[i] * 13.6058 << std::endl;
    }

    delete[] A;
    delete[] B;
    delete[] Z;
}
#endif