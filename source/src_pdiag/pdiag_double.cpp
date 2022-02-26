/*test generalized Stardand double precision symmetric eigenproblem*/
#include "pdiag_double.h"
#include "../module_base/lapack_connector.h"
#include "../src_pw/occupy.h"
#include "../src_pw/global.h"
#include "../src_lcao/local_orbital_charge.h"
#include "../src_io/wf_local.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"


#ifdef __MPI
extern "C"
{
    #include "../module_base/blacs_connector.h"
    #include "my_elpa.h"
	#include "../module_base/scalapack_connector.h"
}
#include "pdgseps.h"
#include "pzgseps.h"
#include "../module_base/lapack_connector.h"
#endif

#include "../src_external/src_test/test_function.h"


#ifdef __MPI
inline int set_elpahandle(elpa_t &handle, const int *desc,const int local_nrows,const int local_ncols, const int nbands)
{
  int error;
  int nprows, npcols, myprow, mypcol;
  Cblacs_gridinfo(desc[1], &nprows, &npcols, &myprow, &mypcol);
  elpa_init(20210430);
  handle = elpa_allocate(&error);
  elpa_set_integer(handle, "na", desc[2], &error);
  elpa_set_integer(handle, "nev", nbands, &error);

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
#endif


inline bool ifElpaHandle(const bool& newIteration, const bool& ifNSCF)
{
    int doHandle = false;
	if(newIteration) doHandle = true;
	if(ifNSCF) doHandle = true;
	return doHandle;
}

int Pdiag_Double::out_hs = 0;
int Pdiag_Double::out_hsR = 0;
int Pdiag_Double::out_lowf = 0;

Pdiag_Double::Pdiag_Double() {}

Pdiag_Double::~Pdiag_Double(){}


void Pdiag_Double::diago_double_begin(
	const int &ik, // k-point index
	Local_Orbital_wfc &lowf,
	double* h_mat, // hamiltonian matrix
    double* s_mat, // overlap matrix
    double* Stmp, 	// because the output Stmp will be different from Sloc2, so we need to copy that.
    double* ekb) // eigenvalues for each k-point and band
{
    const Parallel_Orbitals* pv = this->ParaV = lowf.ParaV;
    
#ifdef TEST_DIAG
	{
		static int istep = 0;
		auto print_matrix_C = [&](const std::string &file_name, double*m)
		{
			std::ofstream ofs(file_name+"-C_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ic=0; ic<pv->ncol; ++ic)
			{
				for(int ir=0; ir<pv->nrow; ++ir)
				{
					const int index=ic*pv->nrow+ir;
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
			for(int ir=0; ir<pv->nrow; ++ir)
			{
				for(int ic=0; ic<pv->ncol; ++ic)
				{
					const int index=ic*pv->nrow+ir;
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
	assert(pv->loc_size > 0);
	assert(GlobalV::NLOCAL > 0);

	char uplo='U';
	const int inc=1;

    int nprocs, myid;
    MPI_Status status;
    MPI_Comm_size(pv->comm_2D, &nprocs);
    MPI_Comm_rank(pv->comm_2D, &myid);

	// parallel diagonalize the
	// H | psi > = S | psi >
	// problem.
	int loc_pos;

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"start solver, ks_solver",GlobalV::KS_SOLVER);
    if(GlobalV::KS_SOLVER=="hpseps")
    {
        double *eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);

        double* Z = new double[pv->loc_size * GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(Z, pv->loc_size * GlobalV::NLOCAL);

        ModuleBase::Memory::record("Pdiag_Double","Z", pv->loc_size * GlobalV::NLOCAL,"double");
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pdgseps");
        BlasConnector::copy(pv->nloc, s_mat, inc, Stmp, inc);
        int loc_size = pv->loc_size;
        pdgseps(pv->comm_2D, GlobalV::NLOCAL, pv->nb, h_mat, Stmp, Z, eigen, pv->MatrixInfo, uplo, loc_size, loc_pos);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pdgseps");

        if(myid <= pv->lastband_in_proc)
        {
            for(int i=0; i<pv->loc_sizes[myid]; i++)
            {
                for(int n=0; n<GlobalV::NLOCAL; n++)
                {
                    pv->Z_LOC[ik][n*pv->loc_sizes[myid] + i] = Z[n*pv->loc_sizes[myid] + i];
                }
            }
        }

        // the eigenvalues.
        //xiaohui modify 2014-06-15, move to the top
        BlasConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;
        //=====================================
        // gather the eigenvectors and
        // distribute them to each processor
        // Z is delete in gath_eig
        //=====================================

        //xiaohui modify 2014-06-18
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig");
        this->gath_eig(DIAG_HPSEPS_WORLD, GlobalV::NLOCAL, Z);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig");
        //delete[] Z; //LiuXh 20171109
	}// HPSEPS method
    else if(GlobalV::KS_SOLVER=="genelpa")
    {
        double *eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);

        long maxnloc; // maximum number of elements in local matrix
        MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
        MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
        lowf.wfc_gamma[ik].create(pv->ncol, pv->nrow);			// Fortran order
        
        static elpa_t handle;
        static bool has_set_elpa_handle = false;
        if(! has_set_elpa_handle)
        {
            set_elpahandle(handle, pv->desc, pv->nrow, pv->ncol, GlobalV::NBANDS);
            has_set_elpa_handle = true;
        }

        int is_already_decomposed;
        if(ifElpaHandle(GlobalC::CHR.get_new_e_iteration(), (GlobalV::CALCULATION=="nscf")))
        {
            ModuleBase::timer::tick("Diago_LCAO_Matrix","decompose_S");
            BlasConnector::copy(pv->nloc, s_mat, inc, Stmp, inc);
            is_already_decomposed=0;
            ModuleBase::timer::tick("Diago_LCAO_Matrix","decompose_S");
        }
        else
        {
            is_already_decomposed=1;
        }

        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_solve");
        int elpa_error;
        std::cout << "before elpa" << std::endl;
        elpa_generalized_eigenvectors_d(handle, h_mat, Stmp, eigen, lowf.wfc_gamma[ik].c, is_already_decomposed, &elpa_error);
        std::cout << "after elpa" << std::endl;
        ModuleBase::timer::tick("Diago_LCAO_Matrix", "elpa_solve");

    	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"K-S equation was solved by genelpa2");
        BlasConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;
	    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"eigenvalues were copied to ekb");

		// convert wave function to band distribution
			// and calculate the density matrix in the tranditional way
			// redistribute eigenvectors to wfc / wfc_aug

		ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig");
		int pos=0;
		for(int i=0; i<myid; ++i)
		{
			pos+=pv->loc_sizes[i];
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
		for(int iprow=0; iprow<pv->dim0; ++iprow)
		{
			for(int ipcol=0; ipcol<pv->dim1; ++ipcol)
			{
				const int coord[2]={iprow, ipcol};
				int src_rank;
				MPI_Cart_rank(pv->comm_2D, coord, &src_rank);
				if(myid==src_rank)
				{
					BlasConnector::copy(pv->nloc_wfc, lowf.wfc_gamma[ik].c, inc, work, inc);
					naroc[0]=pv->nrow;
					naroc[1]=pv->ncol_bands;
				}
				info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, pv->comm_2D);
				info=MPI_Bcast(work, maxnloc, MPI_DOUBLE, src_rank, pv->comm_2D);

				if(out_lowf)
				{
                    info=lowf.q2CTOT(myid, naroc, pv->nb,
                        pv->dim0, pv->dim1, iprow, ipcol, pv->loc_size,
                        work, ctot);
				}//out_lowf
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
		lowf.wfc_gamma[ik].create(pv->ncol_bands, pv->nrow, false);
		memcpy( lowf.wfc_gamma[ik].c, h_mat, sizeof(double)*pv->ncol*pv->nrow );
		ModuleBase::matrix s_tmp(pv->ncol, pv->nrow, false);
		memcpy( s_tmp.c, s_mat, sizeof(double)*pv->ncol*pv->nrow );
		std::vector<double> ekb_tmp(GlobalV::NLOCAL,0);

		const char jobz='V', uplo='U';
		const int itype=1;
		int lwork=-1, info=0;
		std::vector<double> work(1,0);
		dsygv_(&itype, &jobz, &uplo, &GlobalV::NLOCAL, lowf.wfc_gamma[ik].c, &GlobalV::NLOCAL,
			s_tmp.c, &GlobalV::NLOCAL, ekb_tmp.data(), work.data(), &lwork, &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

		lwork = work[0];
		work.resize(lwork);

		dsygv_(&itype, &jobz, &uplo, &GlobalV::NLOCAL, lowf.wfc_gamma[ik].c, &GlobalV::NLOCAL,
			s_tmp.c, &GlobalV::NLOCAL, ekb_tmp.data(), work.data(), &lwork, &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
		memcpy( ekb, ekb_tmp.data(), sizeof(double)*GlobalV::NBANDS );

	}
	else if(GlobalV::KS_SOLVER=="lapack_gvx")
	{
		ModuleBase::matrix h_tmp(pv->ncol, pv->nrow, false);
		memcpy( h_tmp.c, h_mat, sizeof(double)*pv->ncol*pv->nrow );
		ModuleBase::matrix s_tmp(pv->ncol, pv->nrow, false);
		memcpy( s_tmp.c, s_mat, sizeof(double)*pv->ncol*pv->nrow );
		lowf.wfc_gamma[ik].create(pv->ncol, pv->nrow, false);

		const char jobz='V', range='I', uplo='U';
		const int itype=1, il=1, iu=GlobalV::NBANDS;
		int M=0, lwork=-1, info=0;
		const double abstol=0;
		std::vector<double> work(1,0);
		std::vector<int> iwork(5*GlobalV::NLOCAL,0);
		std::vector<int> ifail(GlobalV::NLOCAL,0);

		dsygvx_(&itype, &jobz, &range, &uplo,
			&GlobalV::NLOCAL, h_tmp.c, &GlobalV::NLOCAL, s_tmp.c, &GlobalV::NLOCAL, NULL, NULL, &il, &iu, &abstol,
			&M, ekb, lowf.wfc_gamma[ik].c, &GlobalV::NLOCAL, work.data(), &lwork, iwork.data(), ifail.data(), &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

		lwork = work[0];
		work.resize(lwork);
		dsygvx_(&itype, &jobz, &range, &uplo,
			&GlobalV::NLOCAL, h_tmp.c, &GlobalV::NLOCAL, s_tmp.c, &GlobalV::NLOCAL, NULL, NULL, &il, &iu, &abstol,
			&M, ekb, lowf.wfc_gamma[ik].c, &GlobalV::NLOCAL, work.data(), &lwork, iwork.data(), ifail.data(), &info);

		if(info)
		{
			throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
		if(M!=GlobalV::NBANDS)
		{
			throw std::runtime_error("M="+ModuleBase::GlobalFunc::TO_STRING(M)+". GlobalV::NBANDS="+ModuleBase::GlobalFunc::TO_STRING(GlobalV::NBANDS)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

	}
	else if(GlobalV::KS_SOLVER=="scalapack_gvx")
	{
		diag_scalapack_gvx.pdsygvx_diag(pv->desc, pv->ncol, pv->nrow, h_mat, s_mat, ekb, lowf.wfc_gamma[ik]);		// Peize Lin add 2021.11.02
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
			ofs<<lowf.wfc_gamma[ik]<<std::endl;
		}
		{
			std::ofstream ofs("wfc-F_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			ofs<<transpose(lowf.wfc_gamma[ik])<<std::endl;
		}
		++istep;
	}
#endif

	return;
}

void Pdiag_Double::diago_complex_begin(
	const int &ik,
    Local_Orbital_wfc &lowf,
	std::complex<double>* ch_mat,
    std::complex<double>* cs_mat,
    std::complex<double>* Stmp, 	// because the output Stmp will be different from Sloc2, so we need to copy that.
    double* ekb)
{
    const Parallel_Orbitals* pv = this->ParaV = lowf.ParaV;
    
#ifdef TEST_DIAG
   	{
		static int istep = 0;
		auto print_matrix_C = [&](const std::string &file_name, std::complex<double>*m)
		{
			std::ofstream ofs(file_name+"-C_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
			for(int ic=0; ic<pv->ncol; ++ic)
			{
				for(int ir=0; ir<pv->nrow; ++ir)
				{
					const int index=ic*pv->nrow+ir;
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
			for(int ir=0; ir<pv->nrow; ++ir)
			{
				for(int ic=0; ic<pv->ncol; ++ic)
				{
					const int index=ic*pv->nrow+ir;
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
    MPI_Comm_size(pv->comm_2D, &nprocs);
    MPI_Comm_rank(pv->comm_2D, &myid);

	// parallel diagonalize the
	// H | psi > = S | psi >
	// problem.
	int loc_pos;

	if(GlobalV::KS_SOLVER=="hpseps")
	{
        double *eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);

        assert(pv->loc_size > 0);
        std::complex<double>* Z = new std::complex<double>[pv->loc_size * GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(Z, pv->loc_size * GlobalV::NLOCAL);

        ModuleBase::Memory::record("Pdiag_Double","Z",pv->loc_size * GlobalV::NLOCAL,"cdouble");
		int nbands_tmp = GlobalV::NBANDS;
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pzgseps");
        BlasConnector::copy(pv->nloc, cs_mat, inc, Stmp, inc);
        int loc_size = pv->loc_size;
        pzgseps(pv->comm_2D, GlobalV::NLOCAL, pv->nb, nbands_tmp, ch_mat, Stmp, Z, eigen, pv->MatrixInfo, uplo, loc_size, loc_pos);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","pzgseps");
        // the eigenvalues.
        BlasConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;

        // Z is delete in gath_eig
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");
        this->gath_eig_complex(DIAG_HPSEPS_WORLD, GlobalV::NLOCAL, lowf.wfc_k_grid[ik], Z, ik);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");
        //delete[] Z; //LiuXh 20180329, fix bug of 'double free()'
        //this->gath_full_eig_complex(DIAG_WORLD, GlobalV::NLOCAL, c, Z);
	} // HPSEPS method
    else if(GlobalV::KS_SOLVER=="genelpa")
    {
        double *eigen = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);
        long maxnloc; // maximum number of elements in local matrix
        MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
        MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
        lowf.wfc_k[ik].create(pv->ncol_bands, pv->nrow);            // Fortran order

        static elpa_t handle;
        static bool has_set_elpa_handle = false;
        if(! has_set_elpa_handle)
        {
            set_elpahandle(handle, pv->desc, pv->nrow, pv->ncol, GlobalV::NBANDS);
            has_set_elpa_handle = true;
        }

        BlasConnector::copy(pv->nloc, cs_mat, inc, Stmp, inc);

        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_solve");
        int elpa_derror;
        elpa_generalized_eigenvectors_dc(handle, reinterpret_cast<double _Complex*>(ch_mat),
                                         reinterpret_cast<double _Complex*>(Stmp),
                                         eigen, reinterpret_cast<double _Complex*>(lowf.wfc_k[ik].c), 0, &elpa_derror);
        ModuleBase::timer::tick("Diago_LCAO_Matrix","elpa_solve");

        // the eigenvalues.
        BlasConnector::copy(GlobalV::NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;

        //change eigenvector matrix from block-cycle distribute matrix to column-divided distribute matrix
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");
        std::complex<double> *work=new std::complex<double>[maxnloc]; // work/buffer matrix
        int naroc[2]; // maximum number of row or column
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
                    BlasConnector::copy(pv->nloc_wfc, lowf.wfc_k[ik].c, inc, work, inc);
                    naroc[0]=pv->nrow;
                    naroc[1]=pv->ncol_bands;
                }
                info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, pv->comm_2D);
                info=MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, pv->comm_2D);

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
					info=lowf.q2WFC_CTOT_complex(myid, naroc, pv->nb,
							pv->dim0, pv->dim1, iprow, ipcol,
							work, lowf.wfc_k_grid[ik], ctot);
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
					info=lowf.q2WFC_complex(naroc, pv->nb,
							pv->dim0, pv->dim1, iprow, ipcol,
							work, lowf.wfc_k_grid[ik]);
				}
            }
        }
        delete[] work;
        ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");
    } // GenELPA method
	else if(GlobalV::KS_SOLVER=="scalapack_gvx")
	{
		diag_scalapack_gvx.pzhegvx_diag(pv->desc, pv->ncol, pv->nrow, ch_mat, cs_mat, ekb, lowf.wfc_k[ik]);		// Peize Lin add 2021.11.02

		// the follow will be deleted after finish newdm
		{
			//change eigenvector matrix from block-cycle distribute matrix to column-divided distribute matrix
			ModuleBase::timer::tick("Diago_LCAO_Matrix","gath_eig_complex");

			int info;
			long maxnloc; // maximum number of elements in local matrix
			info=MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
			info=MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
			std::complex<double> *work=new std::complex<double>[maxnloc]; // work/buffer matrix

			int naroc[2]; // maximum number of row or column
			for(int iprow=0; iprow<pv->dim0; ++iprow)
			{
				for(int ipcol=0; ipcol<pv->dim1; ++ipcol)
				{
					const int coord[2]={iprow, ipcol};
					int src_rank;
					info=MPI_Cart_rank(pv->comm_2D, coord, &src_rank);
					if(myid==src_rank)
					{
						BlasConnector::copy(pv->nloc_wfc, lowf.wfc_k[ik].c, inc, work, inc);
						naroc[0]=pv->nrow;
						naroc[1]=pv->ncol_bands;
					}
					info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, pv->comm_2D);
					info=MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, pv->comm_2D);

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
						info=lowf.q2WFC_CTOT_complex(myid, naroc, pv->nb,
								pv->dim0, pv->dim1, iprow, ipcol,
								work, lowf.wfc_k_grid[ik], ctot);
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
						info=lowf.q2WFC_complex(naroc, pv->nb,
								pv->dim0, pv->dim1, iprow, ipcol,
								work, lowf.wfc_k_grid[ik]);
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


