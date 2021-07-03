//column-circle decomposition
#include"./MRRR/mr_interface.h"
#include "pdgseps.h"
#include "../src_parallel/parallel_reduce.h"
#include "pdsyg2st.h"
#include "pdsytrd.h"
#include "pdt2s.h"
#include "pdst2g.h"

void pdgseps(MPI_Comm comm_2D,int n,int nb,double *A, double *B,double *Z,double *eigen,
             LocalMatrix LM,char uplo, int &loc_size,int &loc_pos)
{
    TITLE("Parallel_Diago","pdgseps");
    int i,j,count=0,k=0,incx=1;
    int color,key;

    int dim[2],period[2],coord[2],NPROCS,IAM;
    double norm[n],err;
    double tol=0.00000001;

    MPI_Cart_get(comm_2D,2,dim,period,coord);
    IAM=coord[0]*dim[1]+coord[1];
    NPROCS=dim[0]*dim[1];

    double* diag = new double[n];
    double* off_diag = new double[n];
    double* D = new double[n];
    double* E = new double[n];
    double* Z1 = new double[loc_size*n];

	// (1)
	int wrong_input = 0;
	if(LM.col_num==0 || LM.row_num==0)
	{
		wrong_input = 1;
	}

	Parallel_Reduce::reduce_int_diag( wrong_input );

	//ofs_running << "\n wrong_input = " << wrong_input << endl;
	if(wrong_input > 0)
	{
		ofs_running << "\n col_num == 0 || row_num == 0" << endl;
		WARNING_QUIT("pdgseps","col_num == 0 || row_num == 0");
	}
//time_t s,e;
//s=clock();
    pdsyg2st(uplo,comm_2D,A,LM,B,LM,nb,n);
//e=clock();
//if(IAM==0)printf("g2st time is %d\n",e-s);
	// (2)
//s=clock();
    pdsytrd(comm_2D,LM,n,nb,A,diag,off_diag,norm,uplo);
//e=clock();
//if(IAM==0)printf("sytrd time is %d\n",e-s);
 char jobz='V',range='I';
 double vl=0.0,vu=0.0;
 //int ne=NBANDS;
 int il=1,iu=NBANDS,m=0,nzc=n;
 int *isuppz = new int[n*2];
 bool tryrac = false;
 int lwork = 30 * n;
 double *work = new double[lwork];
 int liwork = 30 * n;
 int *iwork = new int[liwork];
 int iinfo = 0;

double abstol=0.001;
int try1=0;
off_diag[n-1]=0.0;
//s=clock();
  pdstemr_mpi(DIAG_HPSEPS_WORLD, &jobz, &range, &n, diag, off_diag, &vl, &vu, &il, &iu,
                &m, eigen,Z, &n, &nzc, isuppz, &tryrac, work, &lwork, iwork,
                &liwork, &iinfo);
	//2014-07-02
	//2015-05-07
	delete[] isuppz;
	delete[] work;
	delete[] iwork;
  //free(isuppz);
  //free(work);
  //free(iwork);

//e=clock();
//if(IAM==0)printf("mrrr time is %d\n",e-s);
	// (5)
//s=clock();
    pdt2s(comm_2D,n,nb,A,Z,LM,loc_size,norm,uplo);
//e=clock();
//if(IAM==0)printf("t2s time is %d\n",e-s);
//s=clock();
    if (uplo=='l'||(uplo=='L'))
    {
        // (6)
        pdst2g(comm_2D,nb,n,B,Z,LM,loc_size,uplo);
    }

	//ofs_running << "\n local_size = " << loc_size << endl;
    if (uplo=='u'||(uplo=='U'))
    {
        for (i=0; i<loc_size; i++)
        {
            for (j=0; j<n; j++)
            {
                Z1[j*loc_size+i]=Z[i*n+j];
            }
        }
        // (7)
        pdst2g(comm_2D,nb,n,B,Z1,LM,loc_size,uplo);
        for (i=0; i<loc_size; i++)
        {
            for (j=0; j<n; j++)
            {
                Z[i*n+j]=Z1[i*n+j];
            }
        }
    }

    delete[] diag;
//e=clock();
//if(IAM==0)printf("st2g time is %d\n",e-s);
    delete[] off_diag;
    delete[] D;
    delete[] E;
    delete[] Z1;

	return;
}
