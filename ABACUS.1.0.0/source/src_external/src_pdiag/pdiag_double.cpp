/*test generalized Stardand double precision symmetric eigenproblem*/
#include "pdiag_double.h"
#include "blas_interface.h"
#include "../../src_pw/occupy.h"
#include "../../src_pw/global.h"
//#include "../src_pw/global.h"
//xiaohui add 2014-06-20
#include "../../src_lcao/local_orbital_charge.h"

#ifdef __MPI
extern "C"
{   
    #include "Cblacs.h"
    #include "pblas.h"
    #include "scalapack.h"
    #include "my_elpa.h"
}
#include "GenELPA.h"
#include "pdgseps.h"
#include "pzgseps.h"
#endif

//#include "saveMatrix.hpp"

inline int cart2blacs(MPI_Comm comm_2D, int nprows, int npcols, int N, int nblk, int lld, int *desc, int &mpi_comm_rows, int &mpi_comm_cols)
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
    Cblacs_get(comm_2D, 0, &my_blacs_ctxt);
    Cblacs_gridmap(&my_blacs_ctxt, usermap, nprows, nprows, npcols);
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    info=elpa_get_communicators(comm_2D, myprow, mypcol, &mpi_comm_rows, &mpi_comm_cols);
    delete[] usermap;
    int ISRC=0;
    descinit_(desc, &N, &N, &nblk, &nblk, &ISRC, &ISRC, &my_blacs_ctxt, &lld, &info);

    return my_blacs_ctxt;
#else
    return 0;
#endif
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
	TITLE("Pdiag_Double","divide_HS_2d");
	assert(NLOCAL>0);
	assert(DSIZE>0);

#ifdef __MPI
	DIAG_HPSEPS_WORLD=DIAG_WORLD;
#endif

	if(DCOLOR!=0) return; // mohan add 2012-01-13

//	OUT(ofs_running,"NLOCAL",NLOCAL);
//	OUT(ofs_running,"NPROC",NPROC);

	// get the 2D index of computer.
	this->dim0 = (int)sqrt((double)DSIZE); //mohan update 2012/01/13
	//while (NPROC_IN_POOL%dim0!=0)
	while (DSIZE%dim0!=0)
	{
		this->dim0 = dim0 - 1;
	}
	assert(dim0 > 0);
	this->dim1=DSIZE/dim0;
//	testpb=1; //mohan test
	if(testpb)OUT(ofs_running,"dim0",dim0);
	if(testpb)OUT(ofs_running,"dim1",dim1);

#ifdef __MPI
	// mohan add 2011-04-16
	if(NB2D==0)
	{
		if(NLOCAL>0) this->nb = 1;
		if(NLOCAL>500) this->nb = 32;
		if(NLOCAL>1000) this->nb = 64;
	}
	else if(NB2D>0)
	{
		this->nb = NB2D; // mohan add 2010-06-28
	}
	OUT(ofs_running,"nb2d",nb);

	this->set_parameters();

	// call mpi_creat_cart
	this->mpi_creat_cart(&this->comm_2D,this->dim0,this->dim1);

	// call mat_2d
	this->mat_2d(this->comm_2D, NLOCAL, NLOCAL, this->nb, this->MatrixInfo);

	// mohan add 2010-06-29
	this->nrow = this->MatrixInfo.row_num;
	this->ncol = this->MatrixInfo.col_num;
	this->nloc = MatrixInfo.col_num * MatrixInfo.row_num;
    
	// init blacs context for genelpa
    if(KS_SOLVER=="genelpa")
    {
        blacs_ctxt=cart2blacs(comm_2D, dim0, dim1, NLOCAL, nb, nrow, desc, mpi_comm_rows, mpi_comm_cols);
    }
#else // single processor used.
	this->nb = NLOCAL;
	this->nrow = NLOCAL;
	this->ncol = NLOCAL;
	this->nloc = NLOCAL * NLOCAL;
	this->set_parameters();
	MatrixInfo.row_b = 1;
	MatrixInfo.row_num = NLOCAL;
	delete[] MatrixInfo.row_set;
	MatrixInfo.row_set = new int[NLOCAL];
	for(int i=0; i<NLOCAL; i++)
	{
		MatrixInfo.row_set[i]=i;
	}
	MatrixInfo.row_pos=0;

	MatrixInfo.col_b = 1;
	MatrixInfo.col_num = NLOCAL;
	delete[] MatrixInfo.col_set;
	MatrixInfo.col_set = new int[NLOCAL];
	for(int i=0; i<NLOCAL; i++)
	{
		MatrixInfo.col_set[i]=i;
	}	
	MatrixInfo.col_pos=0;
#endif

	assert(nloc>0);
	if(testpb)OUT(ofs_running,"MatrixInfo.row_num",MatrixInfo.row_num);
	if(testpb)OUT(ofs_running,"MatrixInfo.col_num",MatrixInfo.col_num);
	if(testpb)OUT(ofs_running,"nloc",nloc);
	return;
}

void Pdiag_Double::diago_double_begin(const int &ik, double **wfc, 
	double* h_mat, double* s_mat, double* ekb)
{
#ifdef __MPI
	TITLE("Pdiag_Double","diago_begin");
	assert(this->loc_size > 0);
	assert(NLOCAL > 0);

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
	double *eigen = new double[NLOCAL];
	ZEROS(eigen, NLOCAL);

	double* Z = new double[this->loc_size * NLOCAL];
	ZEROS(Z, this->loc_size * NLOCAL);
	
	Memory::record("Pdiag_Double","Z",loc_size * NLOCAL,"double");

	double* Stmp = LM.Sdiag;
	
	//saveMatrix("h_mat", ncol, nrow, h_mat);
	//saveMatrix("s_mat", ncol, nrow, s_mat);

	double start1 = MPI_Wtime();
    OUT(ofs_running,"start solver, KS_SOLVER",KS_SOLVER);
    if(KS_SOLVER=="hpseps")
    {	
		dcopy_(&nloc, s_mat, &inc, Stmp, &inc);
		pdgseps(comm_2D, NLOCAL, nb, h_mat, Stmp, Z, eigen, this->MatrixInfo, uplo, this->loc_size, loc_pos);
	}// HPSEPS method
    else if(KS_SOLVER=="genelpa")
    {
        int maxnloc; // maximum number of elements in local matrix
        MPI_Reduce(&nloc, &maxnloc, 1, MPI_INT, MPI_MAX, 0, comm_2D);
        MPI_Bcast(&maxnloc, 1, MPI_INT, 0, comm_2D);
        double *q=new double[nloc];
        double *work=new double[maxnloc]; // work/buffer matrix
        static int method;
        bool wantEigenVector=true;
        bool wantDebug=true;
        int info;
        int comm_2D_f=MPI_Comm_c2f(comm_2D);

        int THIS_REAL_ELPA_KERNEL_API=12;
        int useQR=0;						// may be changed to input parameter sometime

        if(chr.new_e_iteration)
        {
            timer::tick("Diago_LCAO_Matrix","pdDecomposeRightMatrix2",'G');
            method=0;			
        	//dcopy_(&nloc, s_mat, &inc, Stmp, &inc);
            for(int i=0; i<ncol; ++i)
                dcopy_(&nrow, &s_mat[i], &ncol, &Stmp[i*nrow], &inc);
	        //saveMatrix("Stmp", nrow, ncol, Stmp);

            info=pdDecomposeRightMatrix2(NLOCAL, nrow, ncol, desc,
                                        Stmp, eigen, q, work,
                                        comm_2D_f, mpi_comm_rows, mpi_comm_cols, method,
                                        THIS_REAL_ELPA_KERNEL_API, useQR);
            timer::tick("Diago_LCAO_Matrix","pdDecomposeRightMatrix2",'G');
        }
	    double* Htmp = new double[nloc];
        for(int i=0; i<ncol; ++i)
            dcopy_(&nrow, &h_mat[i], &ncol, &Htmp[i*nrow], &inc);
	    //saveMatrix("Htmp", nrow, ncol, Htmp);
        timer::tick("Diago_LCAO_Matrix","pdSolveEigen2",'G');
        info=pdSolveEigen2(NBANDS, NLOCAL, nrow, ncol, desc,
                          Htmp, Stmp, eigen, q, work,
                          comm_2D_f, mpi_comm_rows, mpi_comm_cols, method,
                          THIS_REAL_ELPA_KERNEL_API, useQR,
                          wantEigenVector, wantDebug);
        delete[] Htmp;
        timer::tick("Diago_LCAO_Matrix","pdSolveEigen2",'G');

        //change eigenvector matrix from block-cycle distribute matrix to column-divided distribute matrix
        int pos=0;
        for(int i=0; i<myid; ++i)
        {
            pos+=loc_sizes[i];
        }
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
                    dcopy_(&nloc, q, &inc, work, &inc);
                    naroc[0]=nrow;
                    naroc[1]=ncol;
                }
                info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, comm_2D);
                info=MPI_Bcast(work, maxnloc, MPI_DOUBLE, src_rank, comm_2D);
                for(int j=0; j<naroc[1]; ++j)
                {
                    int zcol=globalIndex(j, nb, dim1, ipcol)-pos;
                    if(0<=zcol && zcol<loc_size)
                    {
                        for(int i=0; i<naroc[0]; ++i)
                        {
                            int zrow=globalIndex(i, nb, dim0, iprow);
                            Z[zrow*loc_size+zcol]=work[j*naroc[0]+i];
                        }
                    }
                }
            }
        }
        delete[] q;
        delete[] work;
    } // GenELPA method		
	double end1 = MPI_Wtime();
	//xiaohui add 'OUT_LEVEL', 2015-09-16
	if(OUT_LEVEL != "m") OUT(ofs_running,"TIME OF DIAGO (Sec)",end1 - start1);

	int idsize;

	if(myid <= lastband_in_proc)
	{
		for(int i=0; i<loc_sizes[myid]; i++)
		{
			for(int n=0; n<NLOCAL; n++)
			{
				Z_LOC[ik][n*loc_sizes[myid] + i] = Z[n*loc_sizes[myid] + i];
			}
		}
	}
	
	// the eigenvalues.
	//xiaohui modify 2014-06-15, move to the top
	for(int ib=0; ib<NBANDS; ib++)
	{
		ekb[ib] = eigen[ib];
	}
	delete[] eigen;

//	cout << "\n " << setw(6) << "Band" << setw(25) << "Ry" << setw(25) << " eV" << endl;
//	for(int ib=0; ib<NLOCAL; ib++)
//	{
//		cout << " " << setw(6) << ib+1 << setw(25) << eigen[ib] << setw(25)<< eigen[ib] * Ry_to_eV << endl;
//	}

	
	//=====================================
	// gather the eigenvectors and 
	// distribute them to each processor
	// Z is delete in gath_eig
	//=====================================
	
	//xiaohui modify 2014-06-18
	this->gath_eig(DIAG_HPSEPS_WORLD, NLOCAL, wfc, Z);

//  not used anymore
// 	this->gath_full_eig(DIAG_WORLD, NLOCAL, c, Z);

//	ofs_running << setprecision(10);
#endif
	return;
}


void Pdiag_Double::diago_complex_begin(const int &ik, complex<double> **cc, 
	complex<double>* ch_mat, complex<double>* cs_mat, double *ekb)
{
#ifdef __MPI
	TITLE("Pdiag_Double","diago_complex_begin");

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
	double *eigen = new double[NLOCAL];
	ZEROS(eigen, NLOCAL);

	assert(loc_size > 0);
	complex<double>* Z = new complex<double>[this->loc_size * NLOCAL];
	ZEROS(Z, this->loc_size * NLOCAL);
	
	Memory::record("Pdiag_Double","Z",loc_size * NLOCAL,"cdouble");

	// because the output Stmp will be different from Sloc2, so we need to copy that. 
	//complex<double>* Stmp = new complex<double>[nloc];
	complex<double>* Stmp = LM.Sdiag2;
	
	//saveMatrix("ch_mat", ncol, nrow, ch_mat);
	//saveMatrix("cs_mat", ncol, nrow, cs_mat);
	
	if(KS_SOLVER=="hpseps")
	{
		int nbands_tmp = NBANDS;
		zcopy_(&nloc, cs_mat, &inc, Stmp, &inc);
    	pzgseps(comm_2D, NLOCAL, nb, nbands_tmp, ch_mat, Stmp, Z, eigen, this->MatrixInfo, uplo, this->loc_size, loc_pos);
	} // HPSEPS method
    else if(KS_SOLVER=="genelpa")
    {
        int maxnloc; // maximum number of elements in local matrix
        MPI_Reduce(&nloc, &maxnloc, 1, MPI_INT, MPI_MAX, 0, comm_2D);
        MPI_Bcast(&maxnloc, 1, MPI_INT, 0, comm_2D);
        complex<double> *q=new complex<double>[nloc];
        complex<double> *work=new complex<double>[maxnloc]; // work/buffer matrix
        static int method;
        bool wantEigenVector=true;
        bool wantDebug=true;
        int info;
        int comm_2D_f=MPI_Comm_c2f(comm_2D);

        int THIS_REAL_ELPA_KERNEL_API=9;

        if(chr.new_e_iteration)
        {
            timer::tick("Diago_LCAO_Matrix","pdDecomposeRightMatrix2",'G');
            method=0;			
        	//zcopy_(&nloc, cs_mat, &inc, Stmp, &inc);
            for(int i=0; i<ncol; ++i)
                zcopy_(&nrow, &cs_mat[i], &ncol, &Stmp[i*nrow], &inc);
            info=pzDecomposeRightMatrix2(NLOCAL, nrow, ncol, desc,
                                        Stmp, eigen, q, work,
                                        comm_2D_f, mpi_comm_rows, mpi_comm_cols, method,
                                        THIS_REAL_ELPA_KERNEL_API);
            timer::tick("Diago_LCAO_Matrix","pdDecomposeRightMatrix2",'G');
			//saveMatrix("U", nrow, ncol, Stmp);
        }
        timer::tick("Diago_LCAO_Matrix","pdSolveEigen2",'G');
	    complex<double>* Htmp = new complex<double>[nloc];
        //zcopy_(&nloc, ch_mat, &inc, Htmp, &inc);
        for(int i=0; i<ncol; ++i)
            zcopy_(&nrow, &ch_mat[i], &ncol, &Htmp[i*nrow], &inc);

        info=pzSolveEigen2(NBANDS, NLOCAL, nrow, ncol, desc,
                          Htmp, Stmp, eigen, q, work,
                          comm_2D_f, mpi_comm_rows, mpi_comm_cols, method,
                          THIS_REAL_ELPA_KERNEL_API,
                          wantEigenVector, wantDebug);
        delete[] Htmp;
        timer::tick("Diago_LCAO_Matrix","pdSolveEigen2",'G');

        //change eigenvector matrix from block-cycle distribute matrix to column-divided distribute matrix
        int pos=0;
        for(int i=0; i<myid; ++i)
        {
            pos+=loc_sizes[i];
        }
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
                    zcopy_(&nloc, q, &inc, work, &inc);
                    naroc[0]=nrow;
                    naroc[1]=ncol;
                }
                info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, comm_2D);
                info=MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, comm_2D);
                for(int j=0; j<naroc[1]; ++j)
                {
                    int zcol=globalIndex(j, nb, dim1, ipcol)-pos;
                    if(0<=zcol && zcol<loc_size)
                    {
                        for(int i=0; i<naroc[0]; ++i)
                        {
                            int zrow=globalIndex(i, nb, dim0, iprow);
                            Z[zrow*loc_size+zcol]=work[j*naroc[0]+i];
                        }
                    }
                }
            }
        }
        delete[] q;
        delete[] work;
    } // GenELPA method		

	
//	cout << " loc_pos=" << loc_pos << endl;

	//delete[] Stmp;


	/*
	ofs_running << "\n Z:" << endl;
	for(int i=0; i<loc_size; i++)
	{
		ofs_running << " BandWaveFunc " << setw(5) << i; 
		for(int j=0; j<NLOCAL; j++)
		{
			double a = norm(Z[i*NLOCAL+j]);
			if( a < abs(1.0e-6) )
			ofs_running << setw(15) << "0";
			else
			ofs_running << setw(15) << a;
		}
		ofs_running << endl;
	}
	*/
	
	// the eigenvalues.
	for(int ib=0; ib<NBANDS; ib++)
	{
		ekb[ib] = eigen[ib];
	}

//	cout << "\n " << setw(6) << "Band" << setw(25) << "Ry" << setw(25) << " eV" << endl;
//	for(int ib=0; ib<NBANDS; ib++)
//	{
//		cout << " " << setw(6) << ib+1 << setw(25) << eigen[ib] << setw(25)<< eigen[ib] * Ry_to_eV << endl;
//	}

	delete[] eigen;

	// Z is delete in gath_eig
	this->gath_eig_complex(DIAG_HPSEPS_WORLD, NLOCAL, cc, Z, ik);
 	//this->gath_full_eig_complex(DIAG_WORLD, NLOCAL, c, Z);

	/*
	for(int i=0; i<NBANDS; i++)
	{
		cout << " Band " << i;
		for(int j=0; j<NLOCAL; j++)
		{
			cout << c[i][j] << " " ;
		}
		cout << endl;
	}
	*/

#endif
	return;
}



#ifdef __MPI
void Pdiag_Double::readin(const string &fa, const string &fb, const int &nlocal_tot, double *eigen, double *eigvr)
{
    TITLE("Pdiag_Double","readin");

    int coord[2];
    int dim[2];
    int period[2];
    int i,j,tmp1,tmp2;
    int k,loc_size,loc_pos;
    double time1,time2;

    MPI_Comm comm=DIAG_HPSEPS_WORLD,comm_2D,comm_col,comm_row,newcomm;

    dim[0]=(int)sqrt((double)DSIZE);

    while (DSIZE%dim[0]!=0)
    {
        dim[0]=dim[0]-1;
    }
    dim[1]=DSIZE/dim[0];

    // call mpi_creat_cart
    this->mpi_creat_cart(&comm_2D,dim[0],dim[1]);

    // call mat_2d
    this->mat_2d(comm_2D,nlocal_tot,nlocal_tot,nb,MatrixInfo);

    loc_size=nlocal_tot/DSIZE;
    if (DRANK<nlocal_tot%DSIZE) loc_size=loc_size+1;

    ofs_running << " loc_size = " << loc_size;

    /*Distribute the matrix*/
    const int nloc = MatrixInfo.col_num * MatrixInfo.row_num;

    double *A = new double[nloc];
    double *B = new double[nloc];
    double *Z = new double[loc_size*nlocal_tot];
    ZEROS(A, nloc);
    ZEROS(B, nloc);
    ZEROS(Z, loc_size * nlocal_tot);

    ofs_running << "\n Data distribution of H." << endl;
    this->data_distribution(comm_2D,fa,nlocal_tot,nb,A,MatrixInfo);
    ofs_running << "\n Data distribution of S." << endl;
    this->data_distribution(comm_2D,fb,nlocal_tot,nb,B,MatrixInfo);

    time1=MPI_Wtime();
    // call pdgseps
    char uplo = 'U';
    pdgseps(comm_2D,nlocal_tot,nb,A,B,Z,eigen,MatrixInfo,uplo,loc_size,loc_pos);
    time2=MPI_Wtime();
    OUT(ofs_running,"time1",time1);
    OUT(ofs_running,"time2",time2);

    //this->gath_eig(comm,n,eigvr,Z);

    ofs_running << "\n " << setw(6) << "Band" << setw(25) << "Ry" << setw(25) << " eV" << endl;
    for(int i=0; i<nlocal_tot; i++)
    {
        ofs_running << " " << setw(6) << i << setw(25) << eigen[i] << setw(25)<< eigen[i] * 13.6058 << endl;
    }

	

    delete[] A;
    delete[] B;
    delete[] Z;
}
#endif

