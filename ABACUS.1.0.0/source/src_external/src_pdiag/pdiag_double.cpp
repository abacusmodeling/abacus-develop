/*test generalized Stardand double precision symmetric eigenproblem*/
#include "pdiag_double.h"
#include "../../src_pw/occupy.h"
#include "../../src_pw/global.h"
//#include "../src_pw/global.h"
//xiaohui add 2014-06-20
#include "../../src_lcao/local_orbital_charge.h"

#ifdef __MPI
#include "pdgseps.h"
#include "pzgseps.h"
#endif

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

	// parallel diagonalize the 
	// H | psi > = S | psi > 
	// problem.
	int loc_pos;
	double *eigen = new double[NLOCAL];
	ZEROS(eigen, NLOCAL);

	double* Z = new double[this->loc_size * NLOCAL];
	ZEROS(Z, this->loc_size * NLOCAL);
	
	Memory::record("Pdiag_Double","Z",loc_size * NLOCAL,"double");

	//double* Stmp = new double[nloc];
	double* Stmp = LM.Sdiag;
	for(int i=0; i<nloc; i++)
	{
		Stmp[i] = s_mat[i];
	}


	double start1 = MPI_Wtime();
	pdgseps(comm_2D, NLOCAL, nb, h_mat, Stmp, Z, eigen, this->MatrixInfo, uplo, this->loc_size, loc_pos);
	double end1 = MPI_Wtime();
	//xiaohui add 'OUT_LEVEL', 2015-09-16
	if(OUT_LEVEL != "m") OUT(ofs_running,"TIME OF DIAGO (Sec)",end1 - start1);

	//delete[] Stmp;

	int idsize;
	int nprocs,myid;
	MPI_Status status;
	MPI_Comm_size(DIAG_HPSEPS_WORLD,&nprocs);
	MPI_Comm_rank(DIAG_HPSEPS_WORLD,&myid);

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
	for(int i=0; i<nloc; i++)
	{
		Stmp[i] = cs_mat[i];
//		cout << " i=" << i << " Sloc2=" << cs_mat[i] << endl;
	}

	int nbands_tmp = NBANDS;
    pzgseps(comm_2D, NLOCAL, nb, nbands_tmp, ch_mat, Stmp, Z, eigen, this->MatrixInfo, uplo, this->loc_size, loc_pos);

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

