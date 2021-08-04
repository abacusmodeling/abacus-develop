#include "LCAO_evolve.h"
#include "../src_pw/global.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_io/write_HS.h"
#include"../input.h"
#include <complex>
//fuxiang add 2016-10-28

Evolve_LCAO_Matrix::Evolve_LCAO_Matrix(){}
Evolve_LCAO_Matrix::~Evolve_LCAO_Matrix(){}

void Evolve_LCAO_Matrix::evolve_complex_matrix(const int &ik, std::complex<double>** cc, std::complex<double>** cc_init)const
{
	TITLE("Evolve_LCAO_Matrix","evolve_complex_matrix");
	time_t time_start = time(NULL);
	GlobalV::ofs_running << " Start Time : " << ctime(&time_start);

	if (INPUT.tddft==1)
	{
/*
#ifdef __MPI
		this->using_ScaLAPACK_complex(ik, cc, cc_init);
#else
		this->using_LAPACK_complex(ik, cc, cc_init);
#endif
*/
		this->using_LAPACK_complex(ik, cc, cc_init);
	}
	else
	{
		WARNING_QUIT("Evolve_LCAO_Matrix::evolve_complex_matrix","only tddft==1 cando evolve");
	}

	time_t time_end = time(NULL);
	OUT_TIME("evolve(std::complex)", time_start, time_end);
	
	return;
}

void Evolve_LCAO_Matrix::using_LAPACK_complex(const int &ik, std::complex<double>** c, std::complex<double>** c_init)const
{
	                                                                                                                     TITLE("Evolve_LCAO_Matrix","using_LAPACK_complex");

//	Calculate the U operator

	ComplexMatrix Htmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	ComplexMatrix Stmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	for(int i=0; i<GlobalV::NLOCAL; i++)
        {
        	for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                	Htmp(i,j) = GlobalC::LM.Hloc2[i*GlobalV::NLOCAL+j];
                        Stmp(i,j) = GlobalC::LM.Sloc2[i*GlobalV::NLOCAL+j];
                }
        }

/*
	std::cout << " Htmp: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << Htmp(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;

	std::cout << " Stmp: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << Stmp(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;
*/
				
	int INFO;
/*	LapackConnector::zpotrf('U',GlobalV::NLOCAL,Stmp,GlobalV::NLOCAL,&INFO);
	for (int i=0;i<GlobalV::NLOCAL; i++)
	{
		for (int j=0; j<GlobalV::NLOCAL; j++)
		{
			Stmp(i,j)=Stmp(j,i);
		}
	}
	LapackConnector::zpotri('U',GlobalV::NLOCAL,Stmp,GlobalV::NLOCAL,&INFO);
*/

        int LWORK=3*GlobalV::NLOCAL-1; //tmp
        std::complex<double> * WORK = new std::complex<double>[LWORK];
        ZEROS(WORK, LWORK);
        int IPIV[GlobalV::NLOCAL];

        LapackConnector::zgetrf( GlobalV::NLOCAL, GlobalV::NLOCAL, Stmp, GlobalV::NLOCAL, IPIV, &INFO);
        LapackConnector::zgetri( GlobalV::NLOCAL, Stmp, GlobalV::NLOCAL, IPIV, WORK, LWORK, &INFO);

/*
        std::cout << " S^-1: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << Stmp(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;
*/

	ComplexMatrix S_plus_H(GlobalV::NLOCAL,GlobalV::NLOCAL);
	S_plus_H = Stmp*Htmp;

/*
        std::cout << " S^-1*H: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << S_plus_H(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;
*/

	ComplexMatrix Denominator(GlobalV::NLOCAL,GlobalV::NLOCAL);
	for (int i=0; i<GlobalV::NLOCAL; i++)
       	{
               	for (int j=0; j<GlobalV::NLOCAL; j++)
                {
                     /*   real(Denominator(i,j)) = -imag(S_plus_H(i,j));
                        imag(Denominator(i,j)) = real(S_plus_H(i,j));*/
                          Denominator(i,j) = std::complex<double>( -S_plus_H(i,j).imag(), S_plus_H(i,j).real() );
                }
        }
        
        ComplexMatrix Idmat(GlobalV::NLOCAL,GlobalV::NLOCAL);
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        if(i==j) Idmat(i,j) = std::complex<double>(1.0, 0.0);
                       	else Idmat(i,j) = std::complex<double>(0.0, 0.0);
                }
        }
        //double delta_t;
//      delta_t = 0.2;	//identity: fs;
        ComplexMatrix Numerator(GlobalV::NLOCAL,GlobalV::NLOCAL);
        Numerator = Idmat - 0.5*INPUT.mdp.dt*41.34*Denominator;
        Denominator = Idmat + 0.5*INPUT.mdp.dt*41.34*Denominator;

	int info;
        int lwork=3*GlobalV::NLOCAL-1; //tmp
        std::complex<double> * work = new std::complex<double>[lwork];
        ZEROS(work, lwork);
        int ipiv[GlobalV::NLOCAL];
        
        LapackConnector::zgetrf( GlobalV::NLOCAL, GlobalV::NLOCAL, Denominator, GlobalV::NLOCAL, ipiv, &info);
        LapackConnector::zgetri( GlobalV::NLOCAL, Denominator, GlobalV::NLOCAL, ipiv, work, lwork, &info);

        ComplexMatrix U_operator(GlobalV::NLOCAL,GlobalV::NLOCAL);
/*
        std::cout << " Numerator: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << Numerator(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;

        std::cout << " Denominator: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << Denominator(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;
*/

        U_operator = Numerator*Denominator;
/*
	std::cout << "U_operator Success!!!" <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << U_operator(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
	std::cout <<std::endl;
*/

// Calculate wave function at t+delta t
				
//	std::cout << "wave function coe at t+delta t !" << std::endl;

/*	for(int i=0; i<GlobalV::NBANDS; i++)
	{
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			std::cout << c[i][j] << "\t";
		}
		std::cout <<std::endl;
	}
	std::cout << std::endl;
*/

/*	for(int i=0; i<GlobalV::NBANDS; i++)
	{
                for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			std::cout << GlobalC::LOWF.WFC_K[ik][i][j] << "\t";
		}
		std::cout <<std::endl;
	}
	std::cout <<std::endl;
*/

	for(int i=0; i<GlobalV::NBANDS; i++)
	{
		std::complex<double> ccc[GlobalV::NLOCAL];
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{	
			ccc[j] = std::complex<double>(0.0,0.0);
			for(int k=0; k<GlobalV::NLOCAL; k++)
			{
				 ccc[j] += U_operator(j,k)*c_init[i][k];
			}
		}
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			c[i][j] = ccc[j];
			GlobalC::LOWF.WFC_K[ik][i][j] = ccc[j];
		}	
	}

/*	for(int i=0; i<GlobalV::NBANDS; i++)
	{
                for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			std::cout << GlobalC::LOWF.WFC_K[ik][i][j] << "\t";
		}
		std::cout <<std::endl;
	}
*/

/*      std::cout << " c: " <<std::endl;
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << c[i][j] <<"\t";
                }
                std::cout <<std::endl;
        }
        std::cout <<std::endl;
*/
/*
	for(int i=0; i<GlobalV::NBANDS; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << c[i][j] << "\t";
                }
                std::cout <<std::endl;
        }
        std::cout << std::endl;
*/

//	delete[] work;
//	delete[] ipiv;

	return;
}

/*
#include <cmath>
#include <mpi.h>
extern "C"
{
    #include "/export/soft/intel2015/composer_xe_2015.1.133/mkl/include/mkl_pblas.h"
    #include "/export/soft/intel2015/composer_xe_2015.1.133/mkl/include/mkl_blacs.h"
    #include "/export/soft/intel2015/composer_xe_2015.1.133/mkl/include/mkl_scalapack.h"
    #include "/export/soft/intel2015/composer_xe_2015.1.133/mkl/include/mkl_blas.h"
}
#include "../module_base/blas_connector.h"

int Evolve_LCAO_Matrix::using_ScaLAPACK_complex(const int &ik, std::complex<double>** c, std::complex<double>** c_init)const
{
	TITLE("Evolve_LCAO_Matrix", "using_scalapack_complex");
	ComplexMatrix Htmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	ComplexMatrix Stmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	for(int i=0; i<GlobalV::NLOCAL; i++)
        {
        		for(int j=0; j<GlobalV::NLOCAL; j++)
                	{
                		Htmp(i,j) = GlobalC::LM.Hloc2[i*GlobalV::NLOCAL+j];
                           	Stmp(i,j) = GlobalC::LM.Sloc2[i*GlobalV::NLOCAL+j];
                	}
        }

        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
                for(int j=0; j<GlobalV::NLOCAL; j++)
                {
                        std::cout << Stmp(i,j) <<"\t";
                }
                std::cout <<std::endl;
        }
	std::cout <<std::endl;

	std::cout << "11111" << std::endl;

	int nFull, nblk;
	int myid, nprocs, myprow, nprows, mypcol, npcols;
	int my_blacs_ctxt;
	int mpierr, info;
	int narows, nacols;
	int MPIROOT=0;
	char BLACS_LAYOUT='R', BLACS_SCOPE='A';
	int ISRCPROC=0;
	int desc[9];
	
	nFull=GlobalV::NLOCAL;
	nblk=GlobalV::NLOCAL;	

	// set mpi enviroment
	int argc;
	char **argv;
//	MPI_Init(&argc, &argv);
//    	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
//    	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	
	mpierr=MPI_Bcast(&nFull, 1, MPI_INT, MPIROOT, MPI_COMM_WORLD);
	mpierr=MPI_Bcast(&nblk, 1, MPI_INT, MPIROOT, MPI_COMM_WORLD);

	// set blacs parameters 
	for(npcols=int(sqrt(double(GlobalV::NPROC))); npcols>=2; --npcols)
    	{
        	if(GlobalV::NPROC%npcols==0) break;
    	}
    	nprows=GlobalV::NPROC/npcols;
	
	std::cout << "22222" <<std::endl;

//	Cblacs_get(MPI_COMM_WORLD, 0, &my_blacs_ctxt);
//    	Cblacs_gridinit(&my_blacs_ctxt, &BLACS_LAYOUT, nprows, npcols);
//	Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

//	blacs_get_(MPI_COMM_WORLD, 0, &my_blacs_ctxt);

//	blacs_gridinit_(&my_blacs_ctxt, &BLACS_LAYOUT, &nprows, &npcols);
//	blacs_gridinfo_(&my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);


	narows=numroc_(&nFull, &nblk, &myprow, &ISRCPROC, &nprows);
    	nacols=numroc_(&nFull, &nblk, &mypcol, &ISRCPROC, &npcols);
    	descinit_(desc, &nFull, &nFull, &nblk, &nblk, &ISRCPROC, &ISRCPROC, &my_blacs_ctxt, &narows, &info);


	//init main matrix

//	int subMatrixSize=narows*nacols;

	int subMatrixSize=nFull*nFull;
	std::cout << "333333" << std::endl;

	// start
	int allinfo;
    	char transa, transb, side, uplo, diag;
    	int isrc=1, jsrc=1;
    	double alpha, beta;
    	int matrixCols=nacols;
    	int lda;
	int ipiv;
	int lwork, iwork, liwork;


	double *S, *work, *a, *b;
	S=new double[subMatrixSize];
	work=new double[subMatrixSize];

	MPI_Barrier(MPI_COMM_WORLD);

	std::cout << subMatrixSize << std::endl;

	for(int i=0; i<subMatrixSize; ++i) S[i]= real(GlobalC::LM.Sloc2[i]);
//	for(int i=0; i<subMatrixSize; ++i) std::cout << S[i] <<std::endl;

	for(int i=0; i<subMatrixSize; ++i) c[i]= real(GlobalC::LM.Sloc2[i]);

	int inc=1;
	DCOPY(&subMatrixSize, &c, &inc, S, &inc);


	std::cout << "44444" << std::endl;

//	pdgetrf(&nFull, &nFull, S, &isrc, &jsrc, desc, &ipiv, &info);

	std::cout << "55555"	<<std::endl;

	MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	std::cout << "66666" <<std::endl;

	//calculate the inverse
//	pdgetri(&nFull, S, &isrc, &jsrc, desc, &ipiv, work, &lwork, &iwork, &liwork, &info);
	
	std::cout << "77777" << std::endl;

//    	pdgemm_(&S, &Htmp, &nFull, &nFull, &nFull,
//            	&alpha, work, &isrc, &jsrc, desc,
//                    	q,    &isrc, &jsrc, desc,
//            	&beta,  b,    &isrc, &jsrc, desc);

	

	MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	std::cout << "88888"	<< std::endl;

	MPI_Barrier(MPI_COMM_WORLD);

	std::cout << "99999" << std::endl;

//	blacs_gridexit_(&my_blacs_ctxt);

//    	MPI_Finalize();

	std::cout <<"end!!!!!"<< std::endl;
    	return 0;
}

int globalIndex(int localIndex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock=localIndex/nblk;
    gIndex=(iblock*nprocs+myproc)*nblk+localIndex%nblk;
    return gIndex;
}


int localIndex(int globalIndex, int nblk, int nprocs, int& myproc)
{
    myproc=int((globalIndex%(nblk*nprocs))/nblk);
    return int(globalIndex/(nblk*nprocs))*nblk+globalIndex%nblk;
}

*/
