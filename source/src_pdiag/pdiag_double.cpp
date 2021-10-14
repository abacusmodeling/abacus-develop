/*test generalized Stardand double precision symmetric eigenproblem*/
#include "pdiag_double.h"
#include "../module_base/lapack_connector.h"
#include "../src_pw/occupy.h"
#include "../src_pw/global.h"
#include "../src_lcao/local_orbital_charge.h"
#include "../src_io/wf_local.h"


#ifdef __MPI
extern "C"
{
    #include "Cblacs.h"
    #include "my_elpa.h"
	#include "../module_base/scalapack_connector.h"
}
#include "pdgseps.h"
#include "pzgseps.h"
#include "../module_base/lapack_connector.h"
#endif

#include "../src_external/src_test/test_function.h"

inline int globalIndex(int localIndex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock=localIndex/nblk;
    gIndex=(iblock*nprocs+myproc)*nblk+localIndex%nblk;
    return gIndex;
}


inline int localIndex(int globalIndex, int nblk, int nprocs, int& myproc)
{
    myproc=int((globalIndex%(nblk*nprocs))/nblk);
    return int(globalIndex/(nblk*nprocs))*nblk+globalIndex%nblk;
}

inline int cart2blacs(
	MPI_Comm comm_2D,
	int nprows,
	int npcols,
	int N,
	int nblk,
	int lld,
	int *desc)
{
#ifdef __MPI
    int my_blacs_ctxt;
    int myprow, mypcol;
    int *usermap=new int[nprows*npcols];
    int info=0;
    for(int i=0; i<nprows; ++i)
    {
        for(int j=0; j<npcols; ++j)
        {
            int pcoord[2]={i, j};
            MPI_Cart_rank(comm_2D, pcoord, &usermap[i+j*nprows]);
        }
    }
    MPI_Fint comm_2D_f = MPI_Comm_c2f(comm_2D);
    Cblacs_get(comm_2D_f, 0, &my_blacs_ctxt);
    Cblacs_gridmap(&my_blacs_ctxt, usermap, nprows, nprows, npcols);
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    delete[] usermap;
    int ISRC=0;
    descinit_(desc, &N, &N, &nblk, &nblk, &ISRC, &ISRC, &my_blacs_ctxt, &lld, &info);

    return my_blacs_ctxt;
#else
    return 0;
#endif
}

inline int set_elpahandle(elpa_t &handle, int *desc, int local_nrows, int local_ncols)
{
  int error;
  int nprows, npcols, myprow, mypcol;
  Cblacs_gridinfo(desc[1], &nprows, &npcols, &myprow, &mypcol);
  elpa_init(20210430);
  handle = elpa_allocate(&error);
  elpa_set_integer(handle, "na", desc[2], &error);
  elpa_set_integer(handle, "nev", desc[2], &error);

  elpa_set_integer(handle, "local_nrows", local_nrows, &error);

  elpa_set_integer(handle, "local_ncols", local_ncols, &error);

  elpa_set_integer(handle, "nblk", desc[4], &error);

  elpa_set_integer(handle, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_WORLD), &error);

  elpa_set_integer(handle, "process_row", myprow, &error);

  elpa_set_integer(handle, "process_col", mypcol, &error);

  elpa_set_integer(handle, "blacs_context", desc[1], &error);

  elpa_set_integer(handle, "cannon_for_generalized", 0, &error);
   /* Setup */
  elpa_setup(handle);   /* Set tunables */
  return 0;
}

inline int q2CTOT(
	int myid,
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	int loc_size,
	double* work,
	double** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
            if(myid==0) CTOT[igcol][igrow]=work[j*naroc[0]+i];
        }
    }
    return 0;
}

inline int q2ZLOC_WFC(
	int pos,
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	int loc_size,
	double* work,
	double* ZLOC,
	double** WFC)
{
    //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"start q2ZLOC_WFC");
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        int zcol=igcol-pos;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
            if(0<=zcol && zcol<loc_size)
            {
                ZLOC[igrow*loc_size+zcol]=work[j*naroc[0]+i];
            }
	        int mu_local=GlobalC::SGO.trace_lo_tot[igrow];
            if(mu_local>=0 && igrow<GlobalV::NBANDS)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
        }
    }
    //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"WFC was done in q2ZLOC_WFC");
    return 0;
}

inline int q2ZLOC_WFC_WFCAUG(
	int pos,
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	int loc_size,
	double* work,
	double* ZLOC,
	double** WFC,
	double** WFCAUG)
{

    std::stringstream ss;
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
/*        ss<<"local column "<<j<<" nb "<<nb<<" dim1 "<<dim1<<" mypcol "<<ipcol<<" global column (GlobalV::NBANDS) "<<igcol;
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,ss.str());
		ss.str("");*/
        if(igcol>=GlobalV::NBANDS) continue;
        int zcol=igcol-pos;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
/*	        ss<<"    local row "<<i<<" nb "<<nb<<" dim0 "<<dim0<<" myprow "<<iprow<<" global row (GlobalV::NLOCAL)"<<igrow;
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,ss.str());
			ss.str("");*/
            if(0<=zcol && zcol<loc_size)
            {
                ZLOC[igrow*loc_size+zcol]=work[j*naroc[0]+i];
            }
	        int mu_local=GlobalC::SGO.trace_lo_tot[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
	        int mu_aug=GlobalC::LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                WFCAUG[igcol][mu_aug]=work[j*naroc[0]+i];
            }
        }
    }

    //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"WFCAUG was done in q2ZLOC_WFC_WFCAUG");
    return 0;
}

inline int q2ZLOC_WFC_CTOT(
	int myid,
	int pos,
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	int loc_size,
	double* work,
	double* ZLOC,
	double** WFC,
	double** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        int zcol=igcol-pos;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
            if(0<=zcol && zcol<loc_size)
            {
                ZLOC[igrow*loc_size+zcol]=work[j*naroc[0]+i];
            }
	        int mu_local=GlobalC::SGO.trace_lo_tot[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
            if(myid==0) CTOT[igcol][igrow]=work[j*naroc[0]+i];
        }
    }
    return 0;
}

inline int q2ZLOC_WFC_WFCAUG_CTOT(
	int myid,
	int pos,
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	int loc_size,
	double* work,
	double* ZLOC,
	double** WFC,
	double** WFCAUG,
	double** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        int zcol=igcol-pos;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
            if(0<=zcol && zcol<loc_size)
            {
                ZLOC[igrow*loc_size+zcol]=work[j*naroc[0]+i];
            }
	        int mu_local=GlobalC::SGO.trace_lo_tot[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
	        int mu_aug=GlobalC::LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                WFCAUG[igcol][mu_aug]=work[j*naroc[0]+i];
            }
            if(myid==0) CTOT[igcol][igrow]=work[j*naroc[0]+i];
        }
    }
    return 0;
}

inline int q2WFC_complex(
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	std::complex<double>* work,
	std::complex<double>** WFC)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GlobalC::GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
        }
    }
    return 0;
}

inline int q2WFC_WFCAUG_complex(
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	std::complex<double>* work,
	std::complex<double>** WFC,
	std::complex<double>** WFCAUG)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GlobalC::GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
	        int mu_aug=GlobalC::LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                WFCAUG[igcol][mu_aug]=work[j*naroc[0]+i];
            }
        }
    }
    return 0;
}

inline int q2WFC_CTOT_complex(
	int myid,
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	std::complex<double>* work,
	std::complex<double>** WFC,
	std::complex<double>** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GlobalC::GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
            if(myid==0) CTOT[igcol][igrow]=work[j*naroc[0]+i];
        }
    }
    return 0;
}

inline int q2WFC_WFCAUG_CTOT_complex(
	int myid,
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	std::complex<double>* work,
	std::complex<double>** WFC,
	std::complex<double>** WFCAUG,
	std::complex<double>** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GlobalC::GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
	        int mu_aug=GlobalC::LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                WFCAUG[igcol][mu_aug]=work[j*naroc[0]+i];
            }
            if(myid==0) CTOT[igcol][igrow]=work[j*naroc[0]+i];
        }
    }
    return 0;
}

Pdiag_Double::Pdiag_Double()
{
	// default value of nb is 1,
	// but can change to larger value from input.
    nb = 1;
	MatrixInfo.row_set = new int[1];
	MatrixInfo.col_set = new int[1];
}

Pdiag_Double::~Pdiag_Double()
{
	delete[] MatrixInfo.row_set;
	delete[] MatrixInfo.col_set;
}

void Pdiag_Double::divide_HS_2d
(
#ifdef __MPI
	MPI_Comm DIAG_WORLD
#endif
)
{
	ModuleBase::TITLE("Pdiag_Double","divide_HS_2d");
	assert(GlobalV::NLOCAL>0);
	assert(GlobalV::DSIZE>0);

#ifdef __MPI
	DIAG_HPSEPS_WORLD=DIAG_WORLD;
#endif

	if(GlobalV::DCOLOR!=0) return; // mohan add 2012-01-13

	// get the 2D index of computer.
	this->dim0 = (int)sqrt((double)GlobalV::DSIZE); //mohan update 2012/01/13
	//while (GlobalV::NPROC_IN_POOL%dim0!=0)
	while (GlobalV::DSIZE%dim0!=0)
	{
		this->dim0 = dim0 - 1;
	}
	assert(dim0 > 0);
	this->dim1=GlobalV::DSIZE/dim0;

	if(testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"dim0",dim0);
	if(testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"dim1",dim1);

#ifdef __MPI
	// mohan add 2011-04-16
	if(GlobalV::NB2D==0)
	{
		if(GlobalV::NLOCAL>0) this->nb = 1;
		if(GlobalV::NLOCAL>500) this->nb = 32;
		if(GlobalV::NLOCAL>1000) this->nb = 64;
	}
	else if(GlobalV::NB2D>0)
	{
		this->nb = GlobalV::NB2D; // mohan add 2010-06-28
	}
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nb2d",nb);

	this->set_parameters();

	// call mpi_creat_cart
	this->mpi_creat_cart(&this->comm_2D,this->dim0,this->dim1);

	// call mat_2d
	this->mat_2d(this->comm_2D, GlobalV::NLOCAL, GlobalV::NLOCAL, this->nb, this->MatrixInfo);

	// mohan add 2010-06-29
	this->nrow = this->MatrixInfo.row_num;
	this->ncol = this->MatrixInfo.col_num;
	this->nloc = MatrixInfo.col_num * MatrixInfo.row_num;

	// init blacs context for genelpa
    if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")
    {
        blacs_ctxt=cart2blacs(comm_2D, dim0, dim1, GlobalV::NLOCAL, nb, nrow, desc);
    }
#else // single processor used.
	this->nb = GlobalV::NLOCAL;
	this->nrow = GlobalV::NLOCAL;
	this->ncol = GlobalV::NLOCAL;
	this->nloc = GlobalV::NLOCAL * GlobalV::NLOCAL;
	this->set_parameters();
	MatrixInfo.row_b = 1;
	MatrixInfo.row_num = GlobalV::NLOCAL;
	delete[] MatrixInfo.row_set;
	MatrixInfo.row_set = new int[GlobalV::NLOCAL];
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		MatrixInfo.row_set[i]=i;
	}
	MatrixInfo.row_pos=0;

	MatrixInfo.col_b = 1;
	MatrixInfo.col_num = GlobalV::NLOCAL;
	delete[] MatrixInfo.col_set;
	MatrixInfo.col_set = new int[GlobalV::NLOCAL];
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		MatrixInfo.col_set[i]=i;
	}
	MatrixInfo.col_pos=0;
#endif

	assert(nloc>0);
	if(testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MatrixInfo.row_num",MatrixInfo.row_num);
	if(testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MatrixInfo.col_num",MatrixInfo.col_num);
	if(testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nloc",nloc);
	return;
}

void Pdiag_Double::diago_double_begin(
	const int &ik, // k-point index
	double **wfc, // wave functions
	ModuleBase::matrix &wfc_2d, // wave functions in 2d
	double* h_mat, // hamiltonian matrix
	double* s_mat, // overlap matrix
	double* ekb) // eigenvalues for each k-point and band
{
	#ifdef TEST_DIAG
	{
		static int istep = 0;
		auto print_matrix_C = [&](const std::string &file_name, double*m)
		{
			std::ofstream ofs(file_name+"-C_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ic=0; ic<GlobalC::ParaO.ncol; ++ic)
			{
				for(int ir=0; ir<GlobalC::ParaO.nrow; ++ir)
				{
					const int index=ic*GlobalC::ParaO.nrow+ir;
					if(abs(m[index])>1E-10)
						ofs<<m[index]<<"\t";
					else
						ofs<<0<<"\t";
				}
				ofs<<std::endl;
			}
		};
		auto print_matrix_F = [&](const std::string &file_name, double*m)
		{
			std::ofstream ofs(file_name+"-F_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ir=0; ir<GlobalC::ParaO.nrow; ++ir)
			{
				for(int ic=0; ic<GlobalC::ParaO.ncol; ++ic)
				{
					const int index=ic*GlobalC::ParaO.nrow+ir;
					if(abs(m[index])>1E-10)
						ofs<<m[index]<<"\t";
					else
						ofs<<0<<"\t";
				}
				ofs<<std::endl;
			}
		};
		print_matrix_F("H_gamma", h_mat);
		print_matrix_F("S_gamma", s_mat);
		print_matrix_C("H_gamma", h_mat);
		print_matrix_C("S_gamma", s_mat);
		++istep;
	}
	#endif

#ifdef __MPI
	ModuleBase::TITLE("Pdiag_Double","diago_begin");
	assert(this->loc_size > 0);
	assert(GlobalV::NLOCAL > 0);

	char uplo='U';
	const int inc=1;

    int nprocs, myid;
    MPI_Status status;
    MPI_Comm_size(comm_2D, &nprocs);
    MPI_Comm_rank(comm_2D, &myid);

	// parallel diagonalize the
	// H | psi > = S | psi >
	// problem.
	int loc_pos;

	double* Stmp = GlobalC::LM.Sdiag;

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"start solver, GlobalV::KS_SOLVER",GlobalV::KS_SOLVER);
    if(GlobalV::KS_SOLVER=="hpseps")
    {
        double *eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);

        double* Z = new double[this->loc_size * GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(Z, this->loc_size * GlobalV::NLOCAL);

        ModuleBase::Memory::record("Pdiag_Double","Z",loc_size * GlobalV::NLOCAL,"double");
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pdgseps");
		LapackConnector::copy(nloc, s_mat, inc, Stmp, inc);
		pdgseps(comm_2D, GlobalV::NLOCAL, nb, h_mat, Stmp, Z, eigen, this->MatrixInfo, uplo, this->loc_size, loc_pos);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pdgseps");

        if(myid <= lastband_in_proc)
        {
            for(int i=0; i<loc_sizes[myid]; i++)
            {
                for(int n=0; n<GlobalV::NLOCAL; n++)
                {
                    Z_LOC[ik][n*loc_sizes[myid] + i] = Z[n*loc_sizes[myid] + i];
                }
            }
        }

        // the eigenvalues.
        //xiaohui modify 2014-06-15, move to the top
        LapackConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;
        //=====================================
        // gather the eigenvectors and
        // distribute them to each processor
        // Z is delete in gath_eig
        //=====================================

        //xiaohui modify 2014-06-18
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig");
        this->gath_eig(DIAG_HPSEPS_WORLD, GlobalV::NLOCAL, wfc, Z);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig");
        //delete[] Z; //LiuXh 20171109
	}// HPSEPS method
    else if(GlobalV::KS_SOLVER=="genelpa")
    {
        double *eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);

        long maxnloc; // maximum number of elements in local matrix
        MPI_Reduce(&nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, comm_2D);
        MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, comm_2D);
		wfc_2d.create(this->ncol,this->nrow);			// Fortran order

        int is_already_decomposed, elpa_error;
        static elpa_t handle;

        if(GlobalC::CHR.get_new_e_iteration())
        {
            ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_set");
            LapackConnector::copy(nloc, s_mat, inc, Stmp, inc);
            set_elpahandle(handle, desc, nrow, ncol);
            is_already_decomposed=0;
            ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_set");
        }
        else
        {
            is_already_decomposed=1;
        }
        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_solve");
        elpa_generalized_eigenvectors_d(handle, h_mat, Stmp, eigen, wfc_2d.c, is_already_decomposed, &elpa_error);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_solve");

    	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"K-S equation was solved by genelpa2");
        LapackConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;
	    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"eigenvalues were copied to ekb");

		// convert wave function to band distribution
			// and calculate the density matrix in the tranditional way
			// redistribute eigenvectors to wfc / wfc_aug

		ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig");
		int pos=0;
		for(int i=0; i<myid; ++i)
		{
			pos+=loc_sizes[i];
		}
		int naroc[2]; // maximum number of row or column
		double **ctot;

		if(this->out_lowf && myid==0)
		{
			ctot = new double*[GlobalV::NBANDS];
			for (int i=0; i<GlobalV::NBANDS; i++)
			{
				ctot[i] = new double[GlobalV::NLOCAL];
				ModuleBase::GlobalFunc::ZEROS(ctot[i], GlobalV::NLOCAL);
			}
			ModuleBase::Memory::record("Pdiag_Basic","ctot",GlobalV::NBANDS*GlobalV::NLOCAL,"double");
		}

        double *work=new double[maxnloc]; // work/buffer matrix
        int info;
		for(int iprow=0; iprow<dim0; ++iprow)
		{
			for(int ipcol=0; ipcol<dim1; ++ipcol)
			{
				const int coord[2]={iprow, ipcol};
				int src_rank;
				MPI_Cart_rank(comm_2D, coord, &src_rank);
				if(myid==src_rank)
				{
					LapackConnector::copy(nloc, wfc_2d.c, inc, work, inc);
					naroc[0]=nrow;
					naroc[1]=ncol;
				}
				info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, comm_2D);
				info=MPI_Bcast(work, maxnloc, MPI_DOUBLE, src_rank, comm_2D);

				if(out_lowf)
				{
					if(INPUT.new_dm==0)
					{
						// mohan delete Bfield option 2021-02-12
						info=q2ZLOC_WFC_WFCAUG_CTOT(myid, pos, naroc, nb,
							dim0, dim1, iprow, ipcol, this->loc_size,
							work, Z_LOC[ik], wfc, GlobalC::LOWF.WFC_GAMMA_aug[GlobalV::CURRENT_SPIN], ctot);
					}
					else
					{
						info=q2CTOT(myid, naroc, nb,
							dim0, dim1, iprow, ipcol, this->loc_size,
							work, ctot);
					}
				}//out_lowf
				else
				{
					// mohan update 2021-02-12, delete Bfield option
					info=q2ZLOC_WFC_WFCAUG(pos, naroc, nb,
						dim0, dim1, iprow, ipcol, this->loc_size,
						work, Z_LOC[ik], wfc, GlobalC::LOWF.WFC_GAMMA_aug[GlobalV::CURRENT_SPIN]);
				}
			}//loop ipcol
		}//loop iprow

		if(out_lowf && myid==0)
		{
			std::stringstream ss;
			ss << GlobalV::global_out_dir << "LOWF_GAMMA_S" << GlobalV::CURRENT_SPIN+1 << ".dat";
			// mohan add 2012-04-03, because we need the occupations for the
				// first iteration.
			WF_Local::write_lowf( ss.str(), ctot );//mohan add 2010-09-09
			for (int i=0; i<GlobalV::NBANDS; i++)
			{
				delete[] ctot[i];
			}
			delete[] ctot;
		}

		delete[] work;
		ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig");
	} // GenELPA method
	else if(GlobalV::KS_SOLVER=="lapack_gv")
	{
		wfc_2d.create(this->ncol, this->nrow, false);
		memcpy( wfc_2d.c, h_mat, sizeof(double)*this->ncol*this->nrow );
		ModuleBase::matrix s_tmp(this->ncol, this->nrow, false);
		memcpy( s_tmp.c, s_mat, sizeof(double)*this->ncol*this->nrow );
		std::vector<double> ekb_tmp(GlobalV::NLOCAL,0);

		const char jobz='V', uplo='U';
		const int itype=1;
		int lwork=-1, info=0;
		std::vector<double> work(1,0);
		dsygv_(&itype, &jobz, &uplo, &GlobalV::NLOCAL, wfc_2d.c, &GlobalV::NLOCAL,
			s_tmp.c, &GlobalV::NLOCAL, ekb_tmp.data(), work.data(), &lwork, &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

		lwork = work[0];
		work.resize(lwork);

		dsygv_(&itype, &jobz, &uplo, &GlobalV::NLOCAL, wfc_2d.c, &GlobalV::NLOCAL,
			s_tmp.c, &GlobalV::NLOCAL, ekb_tmp.data(), work.data(), &lwork, &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
		memcpy( ekb, ekb_tmp.data(), sizeof(double)*GlobalV::NBANDS );

		if(INPUT.new_dm==0)
		{
			throw std::domain_error("INPUT.new_dm must be 1. "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
	}
	else if(GlobalV::KS_SOLVER=="lapack_gvx")
	{
		ModuleBase::matrix h_tmp(this->ncol, this->nrow, false);
		memcpy( h_tmp.c, h_mat, sizeof(double)*this->ncol*this->nrow );
		ModuleBase::matrix s_tmp(this->ncol, this->nrow, false);
		memcpy( s_tmp.c, s_mat, sizeof(double)*this->ncol*this->nrow );
		wfc_2d.create(this->ncol, this->nrow, false);

		const char jobz='V', range='I', uplo='U';
		const int itype=1, il=1, iu=GlobalV::NBANDS;
		int M=0, lwork=-1, info=0;
		const double abstol=0;
		std::vector<double> work(1,0);
		std::vector<int> iwork(5*GlobalV::NLOCAL,0);
		std::vector<int> ifail(GlobalV::NLOCAL,0);

		dsygvx_(&itype, &jobz, &range, &uplo,
			&GlobalV::NLOCAL, h_tmp.c, &GlobalV::NLOCAL, s_tmp.c, &GlobalV::NLOCAL, NULL, NULL, &il, &iu, &abstol,
			&M, ekb, wfc_2d.c, &GlobalV::NLOCAL, work.data(), &lwork, iwork.data(), ifail.data(), &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

		lwork = work[0];
		work.resize(lwork);
		dsygvx_(&itype, &jobz, &range, &uplo,
			&GlobalV::NLOCAL, h_tmp.c, &GlobalV::NLOCAL, s_tmp.c, &GlobalV::NLOCAL, NULL, NULL, &il, &iu, &abstol,
			&M, ekb, wfc_2d.c, &GlobalV::NLOCAL, work.data(), &lwork, iwork.data(), ifail.data(), &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
		if(M!=GlobalV::NBANDS)
		{
			throw std::runtime_error("M="+ModuleBase::GlobalFunc::TO_STRING(M)+". GlobalV::NBANDS="+ModuleBase::GlobalFunc::TO_STRING(GlobalV::NBANDS)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

		if(INPUT.new_dm==0)
		{
			throw std::domain_error("INPUT.new_dm must be 1. "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
	}
	else if(GlobalV::KS_SOLVER=="scalapack_gvx")
	{
		ModuleBase::matrix h_tmp(this->ncol, this->nrow, false);
		memcpy( h_tmp.c, h_mat, sizeof(double)*this->ncol*this->nrow );
		ModuleBase::matrix s_tmp(this->ncol, this->nrow, false);
		memcpy( s_tmp.c, s_mat, sizeof(double)*this->ncol*this->nrow );
		wfc_2d.create(this->ncol, this->nrow, false);

		const char jobz='V', range='I', uplo='U';
		const int itype=1, il=1, iu=GlobalV::NBANDS, one=1;
		int M=0, NZ=0, lwork=-1, liwork=-1, info=0;
		const double abstol=0, orfac=-1;
		std::vector<double> work(1,0);
		std::vector<int> iwork(1,0);
		std::vector<int> ifail(GlobalV::NLOCAL,0);
		std::vector<int> iclustr(2*GlobalV::DSIZE);
		std::vector<double> gap(GlobalV::DSIZE);

		pdsygvx_(&itype, &jobz, &range, &uplo,
			&GlobalV::NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
			NULL, NULL, &il, &iu, &abstol,
			&M, &NZ, ekb, &orfac, wfc_2d.c, &one, &one, desc,
			work.data(), &lwork, iwork.data(), &liwork, ifail.data(), iclustr.data(), gap.data(), &info);

		GlobalV::ofs_running<<"lwork="<<work[0]<<"\t"<<"liwork="<<iwork[0]<<std::endl;
		lwork = work[0];
		work.resize(lwork,0);
		liwork = iwork[0];
		iwork.resize(liwork,0);

		pdsygvx_(&itype, &jobz, &range, &uplo,
			&GlobalV::NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
			NULL, NULL, &il, &iu, &abstol,
			&M, &NZ, ekb, &orfac, wfc_2d.c, &one, &one, desc,
			work.data(), &lwork, iwork.data(), &liwork, ifail.data(), iclustr.data(), gap.data(), &info);

		GlobalV::ofs_running<<"M="<<M<<"\t"<<"NZ="<<NZ<<std::endl;

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
		if(M!=GlobalV::NBANDS)
		{
			throw std::runtime_error("M="+ModuleBase::GlobalFunc::TO_STRING(M)+". GlobalV::NBANDS="+ModuleBase::GlobalFunc::TO_STRING(GlobalV::NBANDS)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
		if(M!=NZ)
		{
			throw std::runtime_error("M="+ModuleBase::GlobalFunc::TO_STRING(M)+". NZ="+ModuleBase::GlobalFunc::TO_STRING(NZ)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
		if(INPUT.new_dm==0)
		{
			throw std::domain_error("INPUT.new_dm must be 1. "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
	}
    //delete[] Stmp; //LiuXh 20171109
#endif

#ifdef TEST_DIAG
	{
		static int istep = 0;
		{
			std::ofstream ofs("ekb_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ib=0; ib<GlobalV::NBANDS; ++ib)
			{
				ofs<<ekb[ib]<<std::endl;
			}
		}
		{
			std::ofstream ofs("wfc-C_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			ofs<<wfc_2d<<std::endl;
		}
		{
			std::ofstream ofs("wfc-F_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			ofs<<transpose(wfc_2d)<<std::endl;
		}
		++istep;
	}
#endif

	return;
}


void Pdiag_Double::diago_complex_begin(
	const int &ik,
	std::complex<double> **wfc,
	ModuleBase::ComplexMatrix &wfc_2d,
	std::complex<double>* ch_mat,
	std::complex<double>* cs_mat,
	double *ekb)
{
    #ifdef TEST_DIAG
   	{
		static int istep = 0;
		auto print_matrix_C = [&](const std::string &file_name, std::complex<double>*m)
		{
			std::ofstream ofs(file_name+"-C_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ic=0; ic<GlobalC::ParaO.ncol; ++ic)
			{
				for(int ir=0; ir<GlobalC::ParaO.nrow; ++ir)
				{
					const int index=ic*GlobalC::ParaO.nrow+ir;
					if(std::norm(m[index])>1E-10)
                    {
                        if(std::imag(m[index])>1E-10)
						{
                            ofs<<m[index]<<"\t";
						}
                        else
						{
                            ofs<<std::real(m[index])<<"\t";
						}
                    }
					else
					{
						ofs<<0<<"\t";
					}
				}
				ofs<<std::endl;
			}
		};
		auto print_matrix_F = [&](const std::string &file_name, std::complex<double>*m)
		{
			std::ofstream ofs(file_name+"-F_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ir=0; ir<GlobalC::ParaO.nrow; ++ir)
			{
				for(int ic=0; ic<GlobalC::ParaO.ncol; ++ic)
				{
					const int index=ic*GlobalC::ParaO.nrow+ir;
					if(std::norm(m[index])>1E-10)
                    {
                        if(std::imag(m[index])>1E-10)
						{
                            ofs<<m[index]<<"\t";
						}
                        else
						{
                            ofs<<std::real(m[index])<<"\t";
						}
                    }
					else
					{
						ofs<<0<<"\t";
					}
				}
				ofs<<std::endl;
			}
		};
		print_matrix_F("H_gamma", ch_mat);
		print_matrix_F("S_gamma", cs_mat);
		print_matrix_C("H_gamma", ch_mat);
		print_matrix_C("S_gamma", cs_mat);
		++istep;
	}
    #endif

#ifdef __MPI
	ModuleBase::TITLE("Pdiag_Double","diago_complex_begin");

	char uplo='U';
	const int inc=1;

    int nprocs, myid;
    MPI_Status status;
    MPI_Comm_size(comm_2D, &nprocs);
    MPI_Comm_rank(comm_2D, &myid);

	// parallel diagonalize the
	// H | psi > = S | psi >
	// problem.
	int loc_pos;

	// because the output Stmp will be different from Sloc2, so we need to copy that.
	std::complex<double>* Stmp = GlobalC::LM.Sdiag2;

	if(GlobalV::KS_SOLVER=="hpseps")
	{
        double *eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);

        assert(loc_size > 0);
        std::complex<double>* Z = new std::complex<double>[this->loc_size * GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(Z, this->loc_size * GlobalV::NLOCAL);

        ModuleBase::Memory::record("Pdiag_Double","Z",loc_size * GlobalV::NLOCAL,"cdouble");
		int nbands_tmp = GlobalV::NBANDS;
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pzgseps");
		LapackConnector::copy(nloc, cs_mat, inc, Stmp, inc);
    	pzgseps(comm_2D, GlobalV::NLOCAL, nb, nbands_tmp, ch_mat, Stmp, Z, eigen, this->MatrixInfo, uplo, this->loc_size, loc_pos);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pzgseps");
        // the eigenvalues.
        LapackConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;

        // Z is delete in gath_eig
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");
        this->gath_eig_complex(DIAG_HPSEPS_WORLD, GlobalV::NLOCAL, wfc, Z, ik);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");
        //delete[] Z; //LiuXh 20180329, fix bug of 'double free()'
        //this->gath_full_eig_complex(DIAG_WORLD, GlobalV::NLOCAL, c, Z);
	} // HPSEPS method
    else if(GlobalV::KS_SOLVER=="genelpa")
    {
        double *eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);
        long maxnloc; // maximum number of elements in local matrix
        MPI_Reduce(&nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, comm_2D);
        MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, comm_2D);
        wfc_2d.create(this->ncol,this->nrow);            // Fortran order

        LapackConnector::copy(nloc, cs_mat, inc, Stmp, inc);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_set");
        static elpa_t handle;

        if(GlobalC::CHR.get_new_e_iteration())
        {
            set_elpahandle(handle, desc, nrow, ncol);
        }
        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_set");
        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_solve");
        int elpa_derror;
        elpa_generalized_eigenvectors_dc(handle, reinterpret_cast<double _Complex*>(ch_mat),
                                         reinterpret_cast<double _Complex*>(Stmp),
                                         eigen, reinterpret_cast<double _Complex*>(wfc_2d.c), 0, &elpa_derror);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_solve");

        // the eigenvalues.
        LapackConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;

        //change eigenvector matrix from block-cycle distribute matrix to column-divided distribute matrix
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");
        std::complex<double> *work=new std::complex<double>[maxnloc]; // work/buffer matrix
        int naroc[2]; // maximum number of row or column
        int info;
        for(int iprow=0; iprow<dim0; ++iprow)
        {
            for(int ipcol=0; ipcol<dim1; ++ipcol)
            {
                const int coord[2]={iprow, ipcol};
                int src_rank;
                MPI_Cart_rank(comm_2D, coord, &src_rank);
                if(myid==src_rank)
                {
                    LapackConnector::copy(nloc, wfc_2d.c, inc, work, inc);
                    naroc[0]=nrow;
                    naroc[1]=ncol;
                }
                info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, comm_2D);
                info=MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, comm_2D);

                if(this->out_lowf)
                {
                    std::complex<double> **ctot;
                    if(myid==0)
                    {
                        ctot = new std::complex<double>*[GlobalV::NBANDS];
                        for (int i=0; i<GlobalV::NBANDS; i++)
                        {
                            ctot[i] = new std::complex<double>[GlobalV::NLOCAL];
                            ModuleBase::GlobalFunc::ZEROS(ctot[i], GlobalV::NLOCAL);
                        }
                        ModuleBase::Memory::record("Pdiag_Basic","ctot",GlobalV::NBANDS*GlobalV::NLOCAL,"cdouble");
                    }
					// mohan update 2021-02-12, delete BFIELD option
					info=q2WFC_WFCAUG_CTOT_complex(myid, naroc, nb,
							dim0, dim1, iprow, ipcol,
							work, wfc, GlobalC::LOWF.WFC_K_aug[ik], ctot);
                    std::stringstream ss;
	                ss << GlobalV::global_out_dir << "LOWF_K_" << ik+1 << ".dat";
                    // mohan add 2012-04-03, because we need the occupations for the
                    // first iteration.
                    WF_Local::write_lowf_complex( ss.str(), ctot, ik );//mohan add 2010-09-09
                    if(myid==0)
                    {
                        for (int i=0; i<GlobalV::NBANDS; i++)
                        {
                            delete[] ctot[i];
                        }
                        delete[] ctot;
                    }
                }
                else
                {
					// mohan update 2021-02-12, delte BFIELD option
					info=q2WFC_WFCAUG_complex(naroc, nb,
							dim0, dim1, iprow, ipcol,
							work, wfc, GlobalC::LOWF.WFC_K_aug[ik]);
				}
            }
        }
        delete[] work;
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");
    } // GenELPA method
	else if(GlobalV::KS_SOLVER=="scalapack_gvx")
	{
		ModuleBase::ComplexMatrix h_tmp(this->ncol, this->nrow, false);
		memcpy( h_tmp.c, ch_mat, sizeof(std::complex<double>)*this->ncol*this->nrow );
		ModuleBase::ComplexMatrix s_tmp(this->ncol, this->nrow, false);
		memcpy( s_tmp.c, cs_mat, sizeof(std::complex<double>)*this->ncol*this->nrow );
		wfc_2d.create(this->ncol, this->nrow, false);

		const char jobz='V', range='I', uplo='U';
		const int itype=1, il=1, iu=GlobalV::NBANDS, one=1;
		int M=0, NZ=0, lwork=-1, lrwork=-1, liwork=-1, info=0;
		const double abstol=0, orfac=-1;
		std::vector<std::complex<double>> work(1,0);
		std::vector<double> rwork(1,0);
		std::vector<int> iwork(1,0);
		std::vector<int> ifail(GlobalV::NLOCAL,0);
		std::vector<int> iclustr(2*GlobalV::DSIZE);
		std::vector<double> gap(GlobalV::DSIZE);

		pzhegvx_(&itype, &jobz, &range, &uplo,
			&GlobalV::NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
			NULL, NULL, &il, &iu, &abstol,
			&M, &NZ, ekb, &orfac, wfc_2d.c, &one, &one, desc,
			work.data(), &lwork, rwork.data(), &lrwork,
			iwork.data(), &liwork, ifail.data(), iclustr.data(), gap.data(), &info);
		if (info)
			throw std::runtime_error("info=" + ModuleBase::GlobalFunc::TO_STRING(info) + ". " + ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line " + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		GlobalV::ofs_running<<"lwork="<<work[0]<<"\t"<<"liwork="<<iwork[0]<<std::endl;
		lwork = work[0].real();
		work.resize(lwork,0);
		lrwork = rwork[0];
		rwork.resize(lrwork,0);
		liwork = iwork[0];
		iwork.resize(liwork,0);
		pzhegvx_(&itype, &jobz, &range, &uplo,
			&GlobalV::NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
			NULL, NULL, &il, &iu, &abstol,
			&M, &NZ, ekb, &orfac, wfc_2d.c, &one, &one, desc,
			work.data(), &lwork, rwork.data(), &lrwork,
			iwork.data(), &liwork, ifail.data(), iclustr.data(), gap.data(), &info);
		GlobalV::ofs_running<<"M="<<M<<"\t"<<"NZ="<<NZ<<std::endl;

		if(info)
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		if(M!=GlobalV::NBANDS)
			throw std::runtime_error("M="+ModuleBase::GlobalFunc::TO_STRING(M)+". GlobalV::NBANDS="+ModuleBase::GlobalFunc::TO_STRING(GlobalV::NBANDS)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		if(M!=NZ)
			throw std::runtime_error("M="+ModuleBase::GlobalFunc::TO_STRING(M)+". NZ="+ModuleBase::GlobalFunc::TO_STRING(NZ)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));

//		if(INPUT.new_dm==0)
//			throw std::domain_error("INPUT.new_dm must be 1. "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		// the follow will be deleted after finish newdm
		{
			//change eigenvector matrix from block-cycle distribute matrix to column-divided distribute matrix
			ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");

			long maxnloc; // maximum number of elements in local matrix
			MPI_Reduce(&nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, comm_2D);
			MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, comm_2D);
			std::complex<double> *work=new std::complex<double>[maxnloc]; // work/buffer matrix

			int naroc[2]; // maximum number of row or column
			for(int iprow=0; iprow<dim0; ++iprow)
			{
				for(int ipcol=0; ipcol<dim1; ++ipcol)
				{
					const int coord[2]={iprow, ipcol};
					int src_rank;
					MPI_Cart_rank(comm_2D, coord, &src_rank);
					if(myid==src_rank)
					{
						LapackConnector::copy(nloc, wfc_2d.c, inc, work, inc);
						naroc[0]=nrow;
						naroc[1]=ncol;
					}
					info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, comm_2D);
					info=MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, comm_2D);

					if(this->out_lowf)
					{
						std::complex<double> **ctot;
						if(myid==0)
						{
							ctot = new std::complex<double>*[GlobalV::NBANDS];
							for (int i=0; i<GlobalV::NBANDS; i++)
							{
								ctot[i] = new std::complex<double>[GlobalV::NLOCAL];
								ModuleBase::GlobalFunc::ZEROS(ctot[i], GlobalV::NLOCAL);
							}
							ModuleBase::Memory::record("Pdiag_Basic","ctot",GlobalV::NBANDS*GlobalV::NLOCAL,"cdouble");
						}
						// mohan update 2021-02-12, delete BFIELD option
						info=q2WFC_WFCAUG_CTOT_complex(myid, naroc, nb,
								dim0, dim1, iprow, ipcol,
								work, wfc, GlobalC::LOWF.WFC_K_aug[ik], ctot);
						std::stringstream ss;
						ss << GlobalV::global_out_dir << "LOWF_K_" << ik+1 << ".dat";
						// mohan add 2012-04-03, because we need the occupations for the
						// first iteration.
						WF_Local::write_lowf_complex( ss.str(), ctot, ik );//mohan add 2010-09-09
						if(myid==0)
						{
							for (int i=0; i<GlobalV::NBANDS; i++)
							{
								delete[] ctot[i];
							}
							delete[] ctot;
						}
					}
					else
					{
						// mohan update 2021-02-12, delte BFIELD option
						info=q2WFC_WFCAUG_complex(naroc, nb,
								dim0, dim1, iprow, ipcol,
								work, wfc, GlobalC::LOWF.WFC_K_aug[ik]);
					}
				}
			}
			delete[] work;
			ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");
		}
	}

#endif
	return;
}



#ifdef __MPI
void Pdiag_Double::readin(
	const std::string &fa,
	const std::string &fb,
	const int &nlocal_tot,
	double *eigen,
	double *eigvr)
{
    ModuleBase::TITLE("Pdiag_Double","readin");

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
    this->mpi_creat_cart(&comm_2D,dim[0],dim[1]);

    // call mat_2d
    this->mat_2d(comm_2D,nlocal_tot,nlocal_tot,nb,MatrixInfo);

    loc_size=nlocal_tot/GlobalV::DSIZE;
    if (GlobalV::DRANK<nlocal_tot%GlobalV::DSIZE) loc_size=loc_size+1;

    GlobalV::ofs_running << " loc_size = " << loc_size;

    /*Distribute the matrix*/
    const long nloc = MatrixInfo.col_num * MatrixInfo.row_num;

    double *A = new double[nloc];
    double *B = new double[nloc];
    double *Z = new double[loc_size*nlocal_tot];
    ModuleBase::GlobalFunc::ZEROS(A, nloc);
    ModuleBase::GlobalFunc::ZEROS(B, nloc);
    ModuleBase::GlobalFunc::ZEROS(Z, loc_size * nlocal_tot);

    GlobalV::ofs_running << "\n Data distribution of H." << std::endl;
    this->data_distribution(comm_2D,fa,nlocal_tot,nb,A,MatrixInfo);
    GlobalV::ofs_running << "\n Data distribution of S." << std::endl;
    this->data_distribution(comm_2D,fb,nlocal_tot,nb,B,MatrixInfo);

    time1=MPI_Wtime();
    // call pdgseps
    char uplo = 'U';
    pdgseps(comm_2D,nlocal_tot,nb,A,B,Z,eigen,MatrixInfo,uplo,loc_size,loc_pos);
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
