//column-circle decomposition
#include "pdt2s.h"
void pdt2s(MPI_Comm comm2D,int N_A,int NB,double *A,double *X,
           LocalMatrix loc_A,int loc_size,double *norm,char uplo)
/*  PSEPS routine (version 1.0) --
 *  Computer Network Information Center, CAS.
 *  July 15, 2006
 *  Purpose
 *  ========
 *  pdt2s computes the eigenvectors of stardand eigenproblem by the
 *  eigenvector of  tridiagonal eigenproblem matrix using Householder vector
 *  Arguments
 *  =========
 *  uplo    (global input) CHARACTER
 *          Specifies whether the upper or lower triangular part of the
 *          Hermitian matrix A is stored:
 *          = 'U':  Upper triangular
 *          = 'L':  Lower triangular
 * comm_2D  (input) MPI_Comm
 *          MPI 2D communicator
 *  N_A       (global input) integer
 *          The order of the matrix A n>=0
 *  A       (local input)double precision  pointer
 *          to point local memory to an array (loc_A.row_num* loc_A.col_num).
 *          if uplo = 'U', the elements above the first superdiagonal represent the unitary
 *          matrix Q as a product of elementary reflectors; if uplo = 'L',  the elements below the
 *          first subdiagonal represent the unitary matrix Q as a product of elementary reflectors.
 *  norm    (local input) double precision pointer .
 *          norm contains the tau which generated Householder transformation I - tau * v * v'
 *  X       (local input)float array
 *          local part of eigenvectors of tridiagonal symmetric matrix
 *  loc_A   (local input) struct Loc_A
 *          This struct avaible stores the information required to establish
 *          the mapping between an object element and its corresponding
 *          process  and memory location.
 */
{
	TITLE("Parallel_Diago","pdt2s");

	/*
	ofs_running << "\n A matrix elements : " << endl;
	for(int i=0; i<loc_A.row_num; i++)
	{
		for(int j=0; j<loc_A.col_num; j++)
		{
			ofs_running << A[j*loc_A.row_num+i] << endl;
		}
	}
	*/

    double  W[N_A][NB],U[N_A][NB],UC[NB][loc_size];
    double w[N_A],u[N_A],p[N_A],u_c[loc_A.col_num],u_r[loc_A.row_num];
    double g,scal,alpha,beta;

    MPI_Comm comm_col,comm_row;
    int i,j,k,incx=1,incy=1,m,n,kk,tmp,myid;
    int cur_col,cur_row;
    int iarow,iacol,dim[2],period[2],coord[2];
    int lda,ldb,ldc,pos,pos1,pos2;
    int bt,bn,bm,first_col,end_col;
    char transa,transb;

    MPI_Cart_get(comm2D,2,dim,period,coord);
    mpi_sub_col(comm2D,&comm_col);
    mpi_sub_row(comm2D,&comm_row);

    loc_A.row_pos=loc_A.row_num-1;
    loc_A.col_pos=loc_A.col_num-1;
    indxg2p(comm2D,NB,N_A-1,N_A-1,&iarow,&iacol);

    if (iarow==coord[0]) loc_A.row_pos--;
    if (iacol==coord[1]) loc_A.col_pos--;
    bn=N_A/NB;
    if (bn*NB<N_A) bn++;

    if (uplo=='U'||uplo=='u')
    {
        for (bt=bn-1; bt>=0; bt--)
        {
            for (i=0; i<N_A; i++)
                for (j=0; j<NB; j++)
                {
                    W[i][j]=0.0;
                    U[i][j]=0.0;
                }
            first_col=bt*NB;
            end_col=(first_col+NB-1>N_A-2)?(N_A-2):(first_col+NB-1);
            bm=end_col-first_col+1;
            for (k=end_col; k>=first_col; k--)
            {
                indxg2p(comm2D,NB,k,k,&iarow,&iacol);
                for (i=k+1; i<N_A; i++)
                    u[i]=p[i]=0.0;
                if (iacol==coord[1]) cur_col=loc_A.col_pos+1;
                else cur_col=loc_A.col_pos;
                if (iarow==coord[0])
                {
                    for (i=cur_col; i<loc_A.col_num; i++)
                        u_c[i]=A[loc_A.row_pos*loc_A.col_num+i];
                    MPI_Bcast(&u_c[cur_col],loc_A.col_num-cur_col,MPI_DOUBLE,iarow,comm_col);
                }
                else
                    MPI_Bcast(&u_c[cur_col],loc_A.col_num-cur_col,MPI_DOUBLE,iarow,comm_col);
                if (iarow==coord[0]) loc_A.row_pos--;
                if (iacol==coord[1]) loc_A.col_pos--;
                for (i=cur_col; i<loc_A.col_num; i++)
                {
                    tmp=loc_A.col_set[i];
                    p[tmp]=u_c[i];
                }
                MPI_Allreduce(&p[k+1],&u[k+1],N_A-k-1,MPI_DOUBLE,MPI_SUM,comm_row);
                for (i=k+1; i<N_A; i++)
                    U[i][k-first_col]=u[i];
                g=norm[k];
                scal=-g;
                if (k==end_col)
                    dscal_(&N_A,&scal,u,&incx);
                else
                {
                    transa='n';
                    m=N_A-k-1;
                    n=end_col-k;
                    beta=0.0;
                    pos=k-first_col+1;
                    alpha=1.0;
                    lda=NB;
                    dgemv_(&transa,&n,&m,&alpha,&U[k+1][pos],&lda,&u[k+1],&incx,&beta,w,&incy);
                    transa='t';
                    beta=1.0;
                    dgemv_(&transa,&n,&m,&alpha,&W[k+1][pos],&lda,w,&incx,&beta,&u[k+1],&incy);
                    dscal_(&N_A,&scal,u,&incx);
                }
                for (i=k+1; i<N_A; i++)
                    W[i][k-first_col]=u[i];
            }
            transa='T';
            transb='T';
            m=bm;
            kk=N_A-first_col-1;
            n=loc_size;
            alpha=1.0;
            beta=0.0;
            lda=N_A;
            ldb=NB;
            ldc=loc_size;
            pos1=first_col+1;
            pos2=first_col+1;
            //dgemm_(&transa,&transb,&n,&m,&kk,&alpha,&X[pos1],&lda,&W[pos2][0],&ldb,&beta,UC,&ldc);
              dgemm_(&transa,&transb,&n,&m,&kk,&alpha,&X[pos1],&lda,&W[pos2][0],&ldb,&beta,&UC[0][0],&ldc);
            transa='T';
            transb='T';
            m=loc_size;
            kk=bm;
            n=N_A-first_col-1;
            lda=NB;
            ldb=loc_size;
            ldc=N_A;
            beta=1.0;
            pos1=first_col+1;
            pos2=first_col+1;
            //dgemm_(&transa,&transb,&n,&m,&kk,&alpha,&U[pos1][0],&lda,UC,&ldb,&beta,&X[pos2],&ldc);
              dgemm_(&transa,&transb,&n,&m,&kk,&alpha,&U[pos1][0],&lda,&UC[0][0],&ldb,&beta,&X[pos2],&ldc);
        }
    }
    else
        for (bt=bn-1; bt>=0; bt--)
        {
            for (i=0; i<N_A; i++)
                for (j=0; j<NB; j++)
                {
                    W[i][j]=0.0;
                    U[i][j]=0.0;
                }
            first_col=bt*NB;
            end_col=(first_col+NB-1>N_A-2)?(N_A-2):(first_col+NB-1);
            bm=end_col-first_col+1;
            for (k=end_col; k>=first_col; k--)
            {
                indxg2p(comm2D,NB,k,k,&iarow,&iacol);
                for (i=k+1; i<N_A; i++)
                    u[i]=p[i]=0.0;
                if (iarow==coord[0]) cur_row=loc_A.row_pos+1;
                else cur_row=loc_A.row_pos;
                if (iacol==coord[1])
                {
                    for (i=cur_row; i<loc_A.row_num; i++)
                        u_r[i]=A[i*loc_A.col_num+loc_A.col_pos];
                    MPI_Bcast(&u_r[cur_row],loc_A.row_num-cur_row,MPI_DOUBLE,iacol,comm_row);
                }
                else
                    MPI_Bcast(&u_r[cur_row],loc_A.row_num-cur_row,MPI_DOUBLE,iacol,comm_row);
                if (iarow==coord[0]) loc_A.row_pos--;
                if (iacol==coord[1]) loc_A.col_pos--;
                for (i=cur_row; i<loc_A.row_num; i++)
                {
                    tmp=loc_A.row_set[i];
                    p[tmp]=u_r[i];
                }
                MPI_Allreduce(&p[k+1],&u[k+1],N_A-k-1,MPI_DOUBLE,MPI_SUM,comm_col);
                for (i=k+1; i<N_A; i++)
                    U[i][k-first_col]=u[i];
                g=norm[k];
                scal=-g;
                incx=1;
                if (k==end_col)
                    dscal_(&N_A,&scal,u,&incx);
                else
                {
                    transa='n';
                    m=N_A-k-2;
                    n=end_col-k;
                    beta=0.0;
                    pos=k-first_col+1;
                    incx=1;
                    incy=1;
                    alpha=1.0;
                    lda=NB;
                    dgemv_(&transa,&n,&m,&alpha,&U[k+2][pos],&lda,&u[k+2],&incx,&beta,w,&incy);
                    transa='t';
                    beta=1.0;
                    dgemv_(&transa,&n,&m,&alpha,&W[k+2][pos],&lda,w,&incx,&beta,&u[k+2],&incy);
                    u[k+1]=1.0;
                    dscal_(&N_A,&scal,u,&incx);
                }
                for (i=k+1; i<N_A; i++)
                    W[i][k-first_col]=u[i];
            }
            transa='T';
            transb='T';
            m=bm;
            kk=N_A-first_col-1;
            n=loc_size;
            alpha=1.0;
            beta=0.0;
            lda=N_A;
            ldb=NB;
            ldc=loc_size;
            pos1=first_col+1;
            pos2=first_col+1;
            //dgemm_(&transa,&transb,&n,&m,&kk,&alpha,&X[pos1],&lda,&W[pos2][0],&ldb,&beta,UC,&ldc);
              dgemm_(&transa,&transb,&n,&m,&kk,&alpha,&X[pos1],&lda,&W[pos2][0],&ldb,&beta,&UC[0][0],&ldc);
            transa='T';
            transb='T';
            m=loc_size;
            kk=bm;
            n=N_A-first_col-1;
            lda=NB;
            ldb=loc_size;
            ldc=N_A;
            beta=1.0;
            pos1=first_col+1;
            pos2=first_col+1;
            //dgemm_(&transa,&transb,&n,&m,&kk,&alpha,&U[pos1][0],&lda,UC,&ldb,&beta,&X[pos2],&ldc);
              dgemm_(&transa,&transb,&n,&m,&kk,&alpha,&U[pos1][0],&lda,&UC[0][0],&ldb,&beta,&X[pos2],&ldc);
        }
    MPI_Comm_free(&comm_col);
    MPI_Comm_free(&comm_row);
	return;
}
