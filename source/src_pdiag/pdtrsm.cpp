#include "pdtrsm.h"
void pdtrsm(char isuplo,int b_n,MPI_Comm comm_2D,int NB,int N_A,
            double *A,double *B, LocalMatrix loc_A,LocalMatrix loc_B)
{
    /*  PSEPS routine (version 1.0) --
     *  Computer Network Information Center, CAS.
     *  December 15, 2005
     *
     *  Purpose
     *  ========
     *  pdtrsm compute triangular matrix equation with multiple right hand sides in parallel
     *  Arguments
     *  =========
     *  isuplo (global input) CHARACTER
     *         Specifies whether the upper or lower triangular part of the Hermitian matrix A is stored: *         = 'U':  Upper triangular
     *         = 'L':  Lower triangular
     *  comm_2D. (input) MPI_Comm
     *         2-D grid MPI communicator
     *  N_A    (global input) INTEGER
     *         The number of columns and rows to be operated on matrices A£¬N >= 0.
     *  NB     (input) INTEGER
     *         blocked size of 2D blocked cyclic matrix
     *  A      (local input/local output) double precision pointer,
     *         pointer into the local memory to an array of dimension (loc_A.row_num, loc_A.col_num).
     *         On entry, this array contains the local pieces of the
     *         N-by-N Hermitian distributed matrix sub( A ). If UPLO = 'U',the leading N-by-N
     *         upper triangular part of sub( A ) contains the upper triangular part of the matrix.
     *         If UPLO = 'L', the leading N-by-N lower triangular part of sub( A ) contains the lower    *         triangular part of the matrix.
     *  B      (local input/local output) double precision pointer,
     *         pointer into the local memory to an array of dimension (loc_B.row_num, loc_B.col_num).
     *          On entry, this array contains the local pieces of the N-by-N Hermitian distributed
     *          matrix sub( B). If UPLO = 'U',the leading N-by-N upper triangular part of sub( A )
     *          contains the upper triangular part of the matrix.  If UPLO = 'L', the leading N-by-N
     *          lower triangular part of sub( A ) contains the lower triangular part of the matrix.
     *  loc_A and Loc_B:  (local input) struct Loc_A
     *          This struct avaible stores the information required to establish
     *          the mapping between an object element and its corresponding
     *          process and memory location.
     */
    int m,n,m1,n1,m2,n2,loc_m,loc_n,lda,ldb,ldc,loc_i,loc_j,loc_i1,loc_j1,loc_i2,loc_j2;
    int coord[2],dim[2],period[2];
    int i,j,bi,bj,k,k1,k2,ia,ja,ia1,ja1,ia2,ja2,km,kn,km1,kn1,km2,kn2;
    int cur_col,cur_row;
    double alpha,beta;

    char transa,transb,diag,side,uplo;
    int iarow,iacol,pos,iarow1,iacol1,iarow2,iacol2;
    double CU[NB][loc_A.col_num],CU_A[loc_A.row_num][NB];
    double CR[loc_A.row_num][NB],CR_A[NB][loc_A.col_num];

    MPI_Comm comm_row,comm_col;

    MPI_Cart_get(comm_2D,2,dim,period,coord);
    mpi_sub_col(comm_2D,&comm_col);
    mpi_sub_row(comm_2D,&comm_row);

    loc_A.row_pos=0;
    loc_A.col_pos=0;
    /*send B to the T  in  the diaganol*/
    if (isuplo=='u'||isuplo=='U')
    {
        /*Compute the bi-th row block ,ia specifies the row number,
          iarow specifies the row number of processors in the 2D grid*/
        ia=0;
        ja=0;
        indxg2p(comm_2D,NB,ia,ja,&iarow,&iacol);
        /*
         * Get the local row position and column position corresponding to
         * the handling part of BThe local column position of row block
         *  in submatrix X To be computed
         */
        if (coord[1]==iacol)
            loc_A.col_pos=loc_A.col_pos+NB;
        if (coord[0]==iarow)
            loc_A.row_pos=loc_A.row_pos+NB;
        /*The local column position of row block in submatrix X To be computed*/
        for (bj=1; bj<b_n; bj++)
        {
            for (i=0; i<NB; i++)
                for (j=0; j<loc_A.col_num; j++)
                    CU[i][j]=0;
            for (i=0; i<loc_A.row_num; i++)
                for (j=0; j<NB; j++)
                    CU_A[i][j]=0;
            /*Compute the row number of processors corresponding to the bj row block*/
            ia1=bj*NB;
            ja1=bj*NB;
            km1=(ia1+NB-1>N_A-1)?(N_A-1):(ia1+NB-1);
            kn1=(ja1+NB-1>N_A-1)?(N_A-1):(ja1+NB-1);
            m1=km1-ia1+1;
            n1=kn1-ja1+1;
            indxg2p(comm_2D,NB,ia1,ja1,&iarow1,&iacol1);
            indxg2l(ia1,ja1,NB,dim[0],dim[1],&loc_i1,&loc_j1);
            cur_col=loc_A.col_pos;
            /*Send bj row to bi row*/
            if (coord[0]==iarow1)
            {
                ldb=loc_B.col_num;
                uplo=' ';
                m=m1;
                n=loc_B.col_num-cur_col;
                pos=loc_A.row_pos*loc_B.col_num+loc_A.col_pos;
                dlacpy_(&uplo,&n,&m,&B[pos],&ldb,&CU[0][cur_col],&ldb);
            }
            MPI_Bcast(CU,NB*loc_A.col_num,MPI_DOUBLE,iarow1,comm_col);
            /*Compute the XB=A for the block part in iarow row*/
            if (coord[1]==iacol1)
            {
                transb='n';
                diag='n';
                side='l';
                uplo='l';
                m=loc_A.row_pos;;
                n=n1;
                lda=loc_A.col_num;
                ldc=loc_A.col_num;
                pos=cur_col;
                alpha=1;
                dtrsm_(&side,&uplo,&transb,&diag,&n,&m,&alpha,&CU[0][cur_col],&ldc,&A[pos],&lda);
                for (i=0; i<m; i++)
                    for (j=0; j<n; j++)
                        CU_A[i][j]=A[i*loc_A.col_num+j+cur_col];
                cur_col=cur_col+n1;

            }
            MPI_Bcast(CU_A,loc_A.row_num*NB,MPI_DOUBLE,iacol1,comm_row);
            /*update A in iarow row*/
            transa='n';
            transb='n';
            m=loc_A.row_pos;
            k=n1;
            n=loc_A.col_num-cur_col;
            alpha=-1.0;
            beta=1.0;
            lda=loc_A.col_num;
            ldb=loc_A.col_num;
            ldc=NB;
            pos=cur_col;/*loc_A.row_pos*loc_A.col_num+cur_col;*/
            //dgemm_(&transa,&transb,&n,&m,&k,&alpha,&CU[0][cur_col],&ldb,CU_A,&ldc,&beta,&A[pos],&lda);
            dgemm_(&transa,&transb,&n,&m,&k,&alpha,&CU[0][cur_col],&ldb,&CU_A[0][0],&ldc,&beta,&A[pos],&lda);

            if (coord[0]==iarow1)
                loc_A.row_pos=loc_A.row_pos+NB;
            if (coord[1]==iacol1)
                loc_A.col_pos=loc_A.col_pos+NB;
        }
    }
    else
    {
        /*
         * Compute the bi-th row block ,ia specifies the row number,
         * iarow specifies the row number of processors in the 2D grid
         */
        ia=0;
        ja=0;
        indxg2p(comm_2D,NB,ia,ja,&iarow,&iacol);
        if (coord[0]==iarow)
            loc_A.row_pos=loc_A.row_pos+NB;
        if (coord[1]==iacol)
            loc_A.col_pos=loc_A.col_pos+NB;
        /*The local column position of row block in submatrix X To be computed*/
        for (bj=1; bj<b_n; bj++)
        {
            for (i=0; i<loc_A.row_num; i++)
                for (j=0; j<NB; j++)
                    CR[i][j]=0.0;
            for (i=0; i<NB; i++)
                for (j=0; j<loc_A.col_num; j++)
                    CR_A[i][j]=0.0;
            /*Compute the row number of processors corresponding to the bj row block*/
            ia1=bj*NB;
            ja1=bj*NB;
            km1=(ia1+NB-1>N_A-1)?(N_A-1):(ia1+NB-1);
            kn1=(ja1+NB-1>N_A-1)?(N_A-1):(ja1+NB-1);
            m1=km1-ia1+1;
            n1=kn1-ja1+1;
            indxg2p(comm_2D,NB,ia1,ja1,&iarow1,&iacol1);
            indxg2l(ia1,ja1,NB,dim[0],dim[1],&loc_i1,&loc_j1);
            cur_row=loc_A.row_pos;
            /*Send bj row to bi row*/
            if (coord[1]==iacol1)
            {
                ldb=loc_B.col_num;
                uplo=' ';
                m=loc_B.row_num-cur_row;
                n=n1;
                ldc=NB;
                pos=loc_A.row_pos*loc_B.col_num+loc_A.col_pos;
                dlacpy_(&uplo,&n,&m,&B[pos],&ldb,&CR[loc_A.row_pos][0],&ldc);
            }
            MPI_Bcast(CR,NB*loc_A.row_num,MPI_DOUBLE,iacol1,comm_row);

            if (coord[0]==iarow1)
            {
                transb='n';
                diag='n';
                side='r';
                uplo='u';
                m=m1;
                n=loc_A.col_pos;
                lda=loc_A.col_num;
                ldc=NB;
                pos=cur_row*loc_A.col_num;
                alpha=1;
                dtrsm_(&side,&uplo,&transb,&diag,&n,&m,&alpha,&CR[cur_row][0],&ldc,&A[pos],&lda);
                for (i=0; i<m; i++)
                    for (j=0; j<n; j++)
                        CR_A[i][j]=A[(i+cur_row)*loc_A.col_num+j];
                cur_row=cur_row+m1;
            }
            MPI_Bcast(CR_A,loc_A.col_num*NB,MPI_DOUBLE,iarow1,comm_col);
            /*update A in iarow row*/
            transa='n';
            transb='n';
            m=loc_A.row_num-cur_row;
            k=m1;
            n=loc_A.col_pos;
            alpha=-1.0;
            beta=1.0;
            lda=loc_A.col_num;
            ldb=NB;
            ldc=loc_A.col_num;
            pos=cur_row*loc_A.col_num;
            //dgemm_(&transa,&transb,&n,&m,&k,&alpha,CR_A,&ldc,&CR[cur_row][0],&ldb,&beta,&A[pos],&lda);
            dgemm_(&transa,&transb,&n,&m,&k,&alpha,&CR_A[0][0],&ldc,&CR[cur_row][0],&ldb,&beta,&A[pos],&lda);
            if (coord[0]==iarow1)
                loc_A.row_pos=loc_A.row_pos+NB;
            if (coord[1]==iacol1)
                loc_A.col_pos=loc_A.col_pos+NB;
        }
    }
	//2015-05-07
	MPI_Comm_free(&comm_col);
	MPI_Comm_free(&comm_row);
}

