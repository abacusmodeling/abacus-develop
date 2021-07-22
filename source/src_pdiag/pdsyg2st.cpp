#include "pdsyg2st.h"
#include "pdtrsm.h"

void pdsyg2st(char isuplo,MPI_Comm comm_2D,double *A,LocalMatrix loc_A,double *B,
              LocalMatrix loc_B,int NB,int N_A)

/*
 * PSEPS routine (version 1.0) --
 * Computer Network Information Center, CAS.
 * December 15, 2005
 *
 * Purpose
 * =======
 * pzheg2st transforms generalized  Hermitian eigenproblem into standard
 * eigenproblem
 * computes Cholesky factorization of an N-by-N real
 * symmetric positive definite matrix B and  L-1AL-T »òU -TAU -1 with a
 * 2D-block-cyclic in parallel. The routine introduces a new parallelization which
 * combines the Cholesky into the transformation from generalized to standard form.
 * Arguments
 * =========
 * NOTE: In the current code, matrix A employ 2D-blocked-cyclic mapping
 *       distribute.
 *
 *  isuplo  (global input) char
 *          = 'U' or 'u':  Upper triangles of sub( A ) and sub( B ) are stored;
 *          = 'L' or 'l':  Lower triangles of sub( A ) and sub( B ) are stored.
 *  comm_2D. (input) MPI_Comm
 *            2-D grid MPI communicator
 *  N_A  (global input) INTEGER
 *      The number of columns and rows to be operated on matrices A£¬N >= 0.
 *  NB  (input) INTEGER
 *       blocked size of 2D blocked cyclic matrix
 *  A       (local input/local output) double precision pointer,
 *          pointer into the local memory to an array of dimension
 *          (loc_A.row_num, loc_A.col_num).
 *          On entry, this array contains the local pieces of the
 *          N-by-N Hermitian distributed matrix sub( A ). If UPLO = 'U',
 *          the leading N-by-N upper triangular part of sub( A ) contains
 *          the upper triangular part of the matrix.  If UPLO = 'L', the
 *          leading N-by-N lower triangular part of sub( A ) contains
 *          the lower triangular part of the matrix.
 *  B       (local input/local output) double precision pointer,
 *          pointer into the local memory to an array of dimension
 *         (loc_B.row_num, loc_B.col_num).
 *          On entry, this array contains the local pieces of the
 *          N-by-N Hermitian distributed matrix sub( B). If UPLO = 'U',
 *          the leading N-by-N upper triangular part of sub( A ) contains
 *          the upper triangular part of the matrix.  If UPLO = 'L', the
 *          leading N-by-N lower triangular part of sub( A ) contains
 *          the lower triangular part of the matrix.
 *  loc_A and Loc_B:  (local input) struct Loc_A
 *          This struct avaible stores the information required to establish
 *          the mapping between an object element and its corresponding
 *          process and memory location.
 */
{
    double alpha,beta;
    double C_A[NB*NB],C_B[NB*NB],C_A1[NB*NB],C_B1[NB*NB];
    int coord[2],crd[2],dim[2],iarow,iacol,iarow1,iacol1,iarow2,iacol2;
    double U1[loc_A.row_num*NB],W1[NB*loc_A.col_num];
    double U_1[loc_A.row_num*NB],W_1[NB*loc_A.col_num];
    int itype=1;
    int in,bt,i,j,k;
    int root_id;
    int pos,ia,ja,ia1,ja1,ia2,ja2,ki,kj,ki1,kj1,ki2,kj2,info;
    int b_n,kk;
    int loc_i,loc_j,loc_i1,loc_j1,loc_i2,loc_j2,period[2];
    MPI_Comm comm_col,comm_row;
    int tag=0;
	int rp,sp;
    MPI_Request *request;
    MPI_Status  status;
    int lda,ldb,ldc;
    char transa,transb,diag,side,uplo;
    int bm,bn,bm1,bn1,m,n,gi,tmp,gj;
    MPI_Cart_get(comm_2D,2,dim,period,coord);
    mpi_sub_col(comm_2D,&comm_col);
    mpi_sub_row(comm_2D,&comm_row);
    b_n=N_A/NB;
    if (b_n*NB<N_A) b_n++;

    if (dim[0]*dim[1]==1)
    {
        if (isuplo=='u'||isuplo=='U')
        {
            uplo='l';
            lda=N_A;
            ldb=N_A;
            m=N_A;
            dpotf2_(&uplo,&m,B,&ldb,&info);
            dsygs2_(&itype,&uplo,&m,A,&lda,B,&ldb,&info);
        }
        else
        {
            uplo='u';
            lda=N_A;
            ldb=N_A;
            m=N_A;
            dpotf2_(&uplo,&m,B,&ldb,&info);
            dsygs2_(&itype,&uplo,&m,A,&lda,B,&ldb,&info);
        }
	//2015-05-07
	MPI_Comm_free(&comm_col);
	MPI_Comm_free(&comm_row);
        return;
    }
    if (isuplo=='u'||isuplo=='U')
        for (bt=0; bt<b_n; bt++)
        {
            MPI_Barrier(DIAG_HPSEPS_WORLD);
            ia=bt*NB;
            ja=bt*NB;
            ki=(ia+NB-1>N_A-1)?(N_A-1):(ia+NB-1);
            bm=ki-ia+1;
            indxg2p(comm_2D,NB,ia,ja,&iarow,&iacol);
            indxg2l(ia,ja,NB,dim[0],dim[1],&loc_i,&loc_j);
            for (i=0; i<NB; i++)
                for (j=0; j<NB; j++)
                    C_A1[i*NB+j]=C_B1[i*NB+j]=C_A[i*NB+j]=C_B[i*NB+j]=0.0;
            for (i=0; i<loc_A.row_num; i++)
                for (j=0; j<NB; j++)
                    U1[i*NB+j]=U_1[i*NB+j]=0.0;
            for (i=0; i<NB; i++)
                for (j=0; j<loc_A.col_num; j++)
                    W1[i*loc_A.col_num+j]=W_1[i*loc_A.col_num+j]=0.0;

            if (coord[0]==iarow)
            {
                if (coord[1]==iacol)
                {
                    uplo='l';
                    lda=loc_A.col_num;
                    ldb=loc_B.col_num;
                    m=bm;
                    pos=loc_i*loc_A.col_num+loc_j;
                    dpotf2_(&uplo,&m,&B[pos],&ldb,&info);
                    dsygs2_(&itype,&uplo,&m,&A[pos],&lda,&B[pos],&ldb,&info);
                    dlacpy_(&uplo,&m,&m,&A[pos],&lda,C_A,&NB);
                    dlacpy_(&uplo,&m,&m,&B[pos],&ldb,C_B,&NB);
                }
                MPI_Bcast(C_A,NB*NB,MPI_DOUBLE,iacol,comm_row);
                MPI_Bcast(C_B,NB*NB,MPI_DOUBLE,iacol,comm_row);
            }
            if (coord[1]==iacol)
            {
                loc_A.col_pos=loc_A.col_pos+bm;
                loc_B.col_pos=loc_B.col_pos+bm;
            }
            if (ki>=N_A-1) break;
            if ((coord[0]==iarow))
            {
                transa='T';
                transb='T';
                diag='n';
                side='r';
                uplo='l';
                m=bm;
                n=loc_A.col_num-loc_A.col_pos;
                lda=loc_A.col_num;
                ldb=loc_A.col_num;
                ldc=NB;
                pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
                alpha=1.0;

                //printf("(%d,%d), in pdsyg2st n=%d bm=%d, col_pos=%d, ki=%d\n",coord[0],coord[1],n,bm,loc_A.col_pos,ki);
                //printf("(%d,%d), in pdsyg2st n=%d \n",coord[0],coord[1],n);
				
				
                dtrsm_(&side,&uplo,&transa,&diag,&n,&m,&alpha,C_B,&ldc,&A[pos],&lda);
                dtrsm_(&side,&uplo,&transb,&diag,&n,&m,&alpha,C_B,&ldc,&B[pos],&ldb);


                alpha=-0.5;
                beta=1.0;
                dsymm_(&side,&uplo,&n,&m,&alpha,C_A,&ldc,&B[pos],&ldb,&beta,&A[pos],&lda);

                for (i=0; i<m; i++)
                    for (j=loc_A.col_pos; j<loc_A.col_num; j++)
                    {
                        W1[i*loc_A.col_num+j]=A[(i+loc_A.row_pos)*loc_A.col_num+j];
                        W_1[i*loc_B.col_num+j]=B[(i+loc_B.row_pos)*loc_B.col_num+j];
                    }
                loc_A.row_pos=loc_A.row_pos+m;
                loc_B.row_pos=loc_B.row_pos+m;
            }
            MPI_Bcast(W1,NB*loc_A.col_num,MPI_DOUBLE,iarow,comm_col);
            MPI_Bcast(W_1,NB*loc_A.col_num,MPI_DOUBLE,iarow,comm_col);
            for (i=bt+1; i<b_n; i++)
            {
                ia1=bt*NB;
                ja1=i*NB;
                ki1=(ia1+NB-1>N_A-1)?(N_A-1):(ia1+NB-1);
                kj1=(ja1+NB-1>N_A-1)?(N_A-1):(ja1+NB-1);
                bm1=ki1-ia1+1;
                bn1=kj1-ja1+1;
                indxg2p(comm_2D,NB,ia1,ja1,&iarow1,&iacol1);
                indxg2l(ia1,ja1,NB,dim[0],dim[1],&loc_i1,&loc_j1);
                crd[0]=iarow1;
                crd[1]=iacol1;
                MPI_Cart_rank(comm_2D,crd,&sp);
                /*diagonal blocked position*/
                ia2=ja1;
                ja2=ia1;
                indxg2p(comm_2D,NB,ia2,ja2,&iarow2,&iacol2);
                indxg2l(ia2,ja2,NB,dim[0],dim[1],&loc_i2,&loc_j2);
                crd[0]=iarow2;
                crd[1]=iacol2;
                MPI_Cart_rank(comm_2D,crd,&rp);
                lda=ldb=loc_A.col_num;
                ldc=NB;

                if (sp==rp)
                {
                    if (coord[0]==iarow1&&coord[1]==iacol1)
                    {
                        uplo=' ';
                        pos=loc_j1;
                        dlacpy_(&uplo,&bn1,&bm1,&W1[pos],&lda,C_A1,&ldc);
                        dlacpy_(&uplo,&bn1,&bm1,&W_1[pos],&lda,C_B1,&ldc);

                        for (k=0; k<bm1; k++)
                            for (j=0; j<bn1; j++)
                            {
                                U1[(j+loc_i2)*NB+k]=C_A1[k*NB+j];
                                U_1[(j+loc_i2)*NB+k]=C_B1[k*NB+j];
                            }
                    }
                }
                else if (coord[0]==iarow1&&coord[1]==iacol1)
                {
                    uplo=' ';
                    pos=loc_j1;
                    dlacpy_(&uplo,&bn1,&bm1,&W1[pos],&lda,C_A1,&ldc);
                    dlacpy_(&uplo,&bn1,&bm1,&W_1[pos],&ldb,C_B1,&ldc);
                    MPI_Send(C_A1,NB*NB,MPI_DOUBLE,rp,tag,DIAG_HPSEPS_WORLD);
                    MPI_Send(C_B1,NB*NB,MPI_DOUBLE,rp,tag,DIAG_HPSEPS_WORLD);
                }
                else if (coord[0]==iarow2&&coord[1]==iacol2)
                {
                    MPI_Recv(C_A1,NB*NB,MPI_DOUBLE,sp,tag,DIAG_HPSEPS_WORLD,&status);
                    MPI_Recv(C_B1,NB*NB,MPI_DOUBLE,sp,tag,DIAG_HPSEPS_WORLD,&status);
                    for (k=0; k<bm1; k++)
                        for (j=0; j<bn1; j++)
                        {
                            U1[(j+loc_i2)*NB+k]=C_A1[k*NB+j];
                            U_1[(j+loc_i2)*NB+k]=C_B1[k*NB+j];
                        }
                }
            }
            MPI_Bcast(U1,NB*loc_A.row_num,MPI_DOUBLE,iacol,comm_row);
            MPI_Bcast(U_1,NB*loc_A.row_num,MPI_DOUBLE,iacol,comm_row);

            transa='n';
            transb='n';
            m=loc_A.row_num-loc_A.row_pos;
            n=loc_A.col_num-loc_A.col_pos;
            k=bm;
            alpha=-1.0;
            beta=1.0;
            lda=loc_A.col_num;
            ldb=NB;
            ldc=loc_A.col_num;
            pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
            dgemm_(&transa,&transb,&n,&m,&k,&alpha,&W_1[loc_A.col_pos],&lda,
                   &U_1[loc_A.row_pos*NB],&ldb,&beta,&B[pos],&ldc);
            dgemm_(&transa,&transb,&n,&m,&k,&alpha,&W1[loc_A.col_pos],&lda,
                   &U_1[loc_A.row_pos*NB],&ldb,&beta,&A[pos],&ldc);
            dgemm_(&transa,&transb,&n,&m,&k,&alpha,&W_1[loc_A.col_pos],&lda,
                   &U1[loc_A.row_pos*NB],&ldb,&beta,&A[pos],&ldc);
            if ((coord[0]==iarow))
            {
                lda=loc_A.col_num;
                ldb=loc_A.col_num;
                m=bm;
                n=loc_A.col_num-loc_A.col_pos;
                ldc=NB;
                pos=(loc_A.row_pos-bm)*loc_A.col_num+loc_A.col_pos;
                alpha=-0.5;
                beta=1.0;
                side='r';
                uplo='l';
                dsymm_(&side,&uplo,&n,&m,&alpha,C_A,&ldc,&B[pos],&ldb,&beta,&A[pos],&lda);
            }
        }
    else
        for (bt=0; bt<b_n; bt++)
        {
            ia=bt*NB;
            ja=bt*NB;
            ki=(ia+NB-1>N_A-1)?(N_A-1):(ia+NB-1);
            bm=ki-ia+1;
            indxg2p(comm_2D,NB,ia,ja,&iarow,&iacol);
            indxg2l(ia,ja,NB,dim[0],dim[1],&loc_i,&loc_j);
            for (i=0; i<NB; i++)
                for (j=0; j<NB; j++)
                    C_A[i*NB+j]=C_B[i*NB+j]=0.0;
            for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                for (j=0; j<NB; j++)
                    U1[i*NB+j]=U_1[i*NB+j]=0.0;
            for (i=0; i<NB; i++)
                for (j=loc_A.col_pos; j<loc_A.col_num; j++)
                    W1[i*loc_A.col_num+j]=W_1[i*loc_A.col_num+j]=0.0;
            if (coord[1]==iacol)
            {
                if (coord[0]==iarow)
                {
                    uplo='u';
                    lda=loc_A.col_num;
                    ldb=loc_B.col_num;
                    pos=loc_i*loc_A.col_num+loc_j;
                    m=bm;
                    dpotf2_(&uplo,&m,&B[pos],&ldb,&info);
                    dsygs2_(&itype,&uplo,&m,&A[pos],&lda,&B[pos],&ldb,&info);

                    dlacpy_(&uplo,&m,&m,&A[pos],&lda,C_A,&NB);
                    dlacpy_(&uplo,&m,&m,&B[pos],&ldb,C_B,&NB);
                }
                MPI_Bcast(C_A,NB*NB,MPI_DOUBLE,iarow,comm_col);
                MPI_Bcast(C_B,NB*NB,MPI_DOUBLE,iarow,comm_col);
            }
            if (coord[0]==iarow)
            {
                loc_A.row_pos=loc_A.row_pos+bm;
                loc_B.row_pos=loc_B.row_pos+bm;
            }
            if (ki>=N_A-1) break;
            if ((coord[1]==iacol))
            {
                transa='T';
                transb='T';
                diag='n';
                side='l';
                uplo='u';
                m=loc_A.row_num-loc_A.row_pos;
                n=bm;
                lda=loc_A.col_num;
                ldb=loc_A.col_num;
                ldc=NB;
                pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
                alpha=1.0;

                //printf("(%d,%d), in pdsyg2st n2=%d \n",coord[0],coord[1],n);
                dtrsm_(&side,&uplo,&transa,&diag,&n,&m,&alpha,C_B,&ldc,&A[pos],&lda);
                dtrsm_(&side,&uplo,&transb,&diag,&n,&m,&alpha,C_B,&ldc,&B[pos],&ldb);
                alpha=-0.5;
                beta=1.0;
                dsymm_(&side,&uplo,&n,&m,&alpha,C_A,&ldc,&B[pos],&ldb,&beta,&A[pos],&lda);
                for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                {
                    for (j=0; j<bm; j++)
                    {
                        U1[i*NB+j]=A[i*loc_A.col_num+j+loc_A.col_pos];
                        U_1[i*NB+j]=B[i*loc_B.col_num+j+loc_B.col_pos];
                    }
                }
                loc_A.col_pos=loc_A.col_pos+bm;
                loc_B.col_pos=loc_B.col_pos+bm;
            }
            /*pos=loc_A.row_pos*NB;*/
            MPI_Bcast(U1,NB*loc_A.row_num,MPI_DOUBLE,iacol,comm_row);
            MPI_Bcast(U_1,NB*loc_A.row_num,MPI_DOUBLE,iacol,comm_row);
            for (i=bt+1; i<b_n; i++)
            {
                MPI_Barrier(DIAG_HPSEPS_WORLD);
                ia1=i*NB;
                ja1=bt*NB;
                ki1=(ia1+NB-1>N_A-1)?(N_A-1):(ia1+NB-1);
                kj1=(ja1+NB-1>N_A-1)?(N_A-1):(ja1+NB-1);
                bm1=ki1-ia1+1;
                bn1=kj1-ja1+1;
                indxg2p(comm_2D,NB,ia1,ja1,&iarow1,&iacol1);
                indxg2l(ia1,ja1,NB,dim[0],dim[1],&loc_i1,&loc_j1);
                crd[0]=iarow1;
                crd[1]=iacol1;
                MPI_Cart_rank(comm_2D,crd,&sp);
                /*diagonal blocked position*/
                ia2=ja1;
                ja2=ia1;
                indxg2p(comm_2D,NB,ia2,ja2,&iarow2,&iacol2);
                indxg2l(ia2,ja2,NB,dim[0],dim[1],&loc_i2,&loc_j2);
                crd[0]=iarow2;
                crd[1]=iacol2;
                MPI_Cart_rank(comm_2D,crd,&rp);
                if (sp==rp)
                {
                    if (coord[0]==iarow1&&coord[1]==iacol1)
                    {
                        uplo=' ';
                        pos=loc_i1*NB;
                        dlacpy_(&uplo,&bn1,&bm1,&U1[pos],&NB,C_A1,&NB);
                        dlacpy_(&uplo,&bn1,&bm1,&U_1[pos],&NB,C_B1,&NB);
                        for (k=0; k<bm1; k++)
                            for (j=0; j<bn1; j++)
                            {
                                W1[j*loc_A.col_num+k+loc_j2]=C_A1[k*NB+j];
                                W_1[j*loc_A.col_num+k+loc_j2]=C_B1[k*NB+j];
                            }
                    }
                }
                else if (coord[0]==iarow1&&coord[1]==iacol1)
                {
                    uplo=' ';
                    pos=loc_i1*NB;
                    dlacpy_(&uplo,&bn1,&bm1,&U1[pos],&NB,C_A1,&NB);
                    dlacpy_(&uplo,&bn1,&bm1,&U_1[pos],&NB,C_B1,&NB);
                    MPI_Send(C_A1,NB*NB,MPI_DOUBLE,rp,tag,DIAG_HPSEPS_WORLD);
                    MPI_Send(C_B1,NB*NB,MPI_DOUBLE,rp,tag,DIAG_HPSEPS_WORLD);
                }
                else if (coord[0]==iarow2&&coord[1]==iacol2)
                {
                    MPI_Recv(C_A1,NB*NB,MPI_DOUBLE,sp,tag,DIAG_HPSEPS_WORLD,&status);
                    MPI_Recv(C_B1,NB*NB,MPI_DOUBLE,sp,tag,DIAG_HPSEPS_WORLD,&status);
                    for (k=0; k<bm1; k++)
                        for (j=0; j<bn1; j++)
                        {
                            W1[j*loc_A.col_num+k+loc_j2]=C_A1[k*NB+j];
                            W_1[j*loc_A.col_num+k+loc_j2]=C_B1[k*NB+j];
                        }
                }
            }
            MPI_Bcast(W1,NB*loc_A.col_num,MPI_DOUBLE,iarow,comm_col);
            MPI_Barrier(comm_row);
            MPI_Bcast(W_1,NB*loc_A.col_num,MPI_DOUBLE,iarow,comm_col);
            MPI_Barrier(comm_row);
            transa='n';
            transb='n';
            m=loc_A.row_num-loc_A.row_pos;
            n=loc_A.col_num-loc_A.col_pos;
            k=bm;
            alpha=-1.0;
            beta=1.0;
            lda=loc_A.col_num;
            ldb=NB;
            ldc=loc_A.col_num;
            pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
            dgemm_(&transa,&transb,&n,&m,&k,&alpha,&W_1[loc_A.col_pos],&lda,
                   &U_1[loc_A.row_pos*NB],&ldb,&beta,&B[pos],&ldc);
            dgemm_(&transa,&transb,&n,&m,&k,&alpha,&W1[loc_A.col_pos],&lda,
                   &U_1[loc_A.row_pos*NB],&ldb,&beta,&A[pos],&ldc);
            dgemm_(&transa,&transb,&n,&m,&k,&alpha,&W_1[loc_A.col_pos],&lda,
                   &U1[loc_A.row_pos*NB],&ldb,&beta,&A[pos],&ldc);
            if ((coord[1]==iacol))
            {
                lda=loc_A.col_num;
                ldb=loc_A.col_num;
                ldc=NB;
                m=loc_A.row_num-loc_A.row_pos;
                n=bm;
                pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos-bm;
                alpha=-0.5;
                beta=1.0;
                side='l';
                uplo='u';
                dsymm_(&side,&uplo,&n,&m,&alpha,C_A,&ldc,&B[pos],&ldb,&beta,&A[pos],&lda);
            }
        }
    pdtrsm(isuplo,b_n,comm_2D,NB,N_A,A,B,loc_A,loc_B);
    MPI_Comm_free(&comm_col);
    MPI_Comm_free(&comm_row);
}

