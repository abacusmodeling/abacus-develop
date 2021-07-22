#include "pzheg2st.h"
#include "pzhtrsm.h" 

void pzheg2st(char isuplo,MPI_Comm comm_2D,complex<double> *A,LocalMatrix loc_A,complex<double> *B,
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
 *  A       (local input/local output) double precision complex pointer,
 *          pointer into the local memory to an array of dimension 
 *          (loc_A.row_num, loc_A.col_num).
 *          On entry, this array contains the local pieces of the
 *          N-by-N Hermitian distributed matrix sub( A ). If UPLO = 'U',
 *          the leading N-by-N upper triangular part of sub( A ) contains
 *          the upper triangular part of the matrix.  If UPLO = 'L', the
 *          leading N-by-N lower triangular part of sub( A ) contains
 *          the lower triangular part of the matrix.
 *  B       (local input/local output) double precision complex pointer,
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
	TITLE("Parallel_Diag","pzheg2st");
    complex<double> alpha,beta;

    complex<double>* C_A = new complex<double>[NB*NB];
    complex<double>* C_B = new complex<double>[NB*NB];
    complex<double>* C_A1 = new complex<double>[NB*NB];
    complex<double>* C_B1 = new complex<double>[NB*NB];

    int coord[2],crd[2],dim[2],iarow,iacol,iarow1,iacol1,iarow2,iacol2;
    complex<double> U1[loc_A.row_num*NB],W1[NB*loc_A.col_num];
    complex<double> U_1[loc_A.row_num*NB],W_1[NB*loc_A.col_num];
    int itype=1;
    int in,bt,i,j,k;
    int root_id;
    int pos,ia,ja,ia1,ja1,ia2,ja2,ki,kj,ki1,kj1,ki2,kj2,info;
    int b_n,kk;
    int loc_i,loc_j,loc_i1,loc_j1,loc_i2,loc_j2,period[2];
    MPI_Comm comm_col,comm_row;
    int tag=0;
	int rp,sp;
    MPI_Status status1;
    MPI_Status status2;
    MPI_Status status3;
    MPI_Status status4;
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
            zpotf2_(&uplo,&m,B,&ldb,&info);
            zhegs2_(&itype,&uplo,&m,A,&lda,B,&ldb,&info);
        }
        else
        {
            uplo='u';
            lda=N_A;
            ldb=N_A;
            m=N_A;
            zpotf2_(&uplo,&m,B,&ldb,&info);
            zhegs2_(&itype,&uplo,&m,A,&lda,B,&ldb,&info);
        }
	//xiaohui add 2014-03-06, test
        MPI_Comm_free(&comm_col);
        MPI_Comm_free(&comm_row);

	//xiaohui add 2015-04-07
	delete[] C_A;
	delete[] C_B;
	delete[] C_A1;
	delete[] C_B1;

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
                {
                    C_A1[i*NB+j]=C_B1[i*NB+j]=C_A[i*NB+j]=C_B[i*NB+j]=complex<double> (0.0,0.0);
                   
                }
            for (i=0; i<loc_A.row_num; i++)
                for (j=0; j<NB; j++)
                    {
                    U1[i*NB+j]=U_1[i*NB+j]=complex<double> (0.0,0.0);
                    }
            for (i=0; i<NB; i++)
                for (j=0; j<loc_A.col_num; j++)
                    {
                    W1[i*loc_A.col_num+j]=W_1[i*loc_A.col_num+j]=complex<double> (0.0,0.0);

                    }
            if (coord[0]==iarow)
            {
                if (coord[1]==iacol)
                {
                    uplo='l';
                    lda=loc_A.col_num;
                    ldb=loc_B.col_num;
                    m=bm;
                    pos=loc_i*loc_A.col_num+loc_j;
                    zpotf2_(&uplo,&m,&B[pos],&ldb,&info);
                    zhegs2_(&itype,&uplo,&m,&A[pos],&lda,&B[pos],&ldb,&info);
                    zlacpy_(&uplo,&m,&m,&A[pos],&lda,C_A,&NB);
                    zlacpy_(&uplo,&m,&m,&B[pos],&ldb,C_B,&NB);
                }
                MPI_Bcast(C_A,NB*NB,MPI_DOUBLE_COMPLEX,iacol,comm_row);
                MPI_Bcast(C_B,NB*NB,MPI_DOUBLE_COMPLEX,iacol,comm_row);
            }
            if (coord[1]==iacol)
            {
                loc_A.col_pos=loc_A.col_pos+bm;
                loc_B.col_pos=loc_B.col_pos+bm;
            }
            if (ki>=N_A-1) break;
            if ((coord[0]==iarow))
            {
                transa='c';
                transb='c';
                diag='n';
                side='r';
                uplo='l';
                m=bm;
                n=loc_A.col_num-loc_A.col_pos;
                lda=loc_A.col_num;
                ldb=loc_A.col_num;
                ldc=NB;
                pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
 
                alpha=complex<double> (1.0,0.0);
                //printf("(%d,%d), in pdsyg2st n=%d bm=%d, col_pos=%d, ki=%d\n",coord[0],coord[1],n,bm,loc_A.col_pos,ki);
                //printf("(%d,%d), in pdsyg2st n=%d \n",coord[0],coord[1],n);	
				
                ztrsm_(&side,&uplo,&transa,&diag,&n,&m,&alpha,C_B,&ldc,&A[pos],&lda);
                ztrsm_(&side,&uplo,&transb,&diag,&n,&m,&alpha,C_B,&ldc,&B[pos],&ldb);


                alpha=complex<double> (-0.5,0.0);
                beta=complex<double> (1.0,0.0);

                zhemm_(&side,&uplo,&n,&m,&alpha,C_A,&ldc,&B[pos],&ldb,&beta,&A[pos],&lda);

                for (i=0; i<m; i++)
                    for (j=loc_A.col_pos; j<loc_A.col_num; j++)
                    {
                        W1[i*loc_A.col_num+j]=A[(i+loc_A.row_pos)*loc_A.col_num+j];
                        W_1[i*loc_B.col_num+j]=B[(i+loc_B.row_pos)*loc_B.col_num+j];
                    }
                loc_A.row_pos=loc_A.row_pos+m;
                loc_B.row_pos=loc_B.row_pos+m;
            }
            MPI_Bcast(W1,NB*loc_A.col_num,MPI_DOUBLE_COMPLEX,iarow,comm_col);

			//cout << "\n Bcast done." << endl;
			
            MPI_Bcast(W_1,NB*loc_A.col_num,MPI_DOUBLE_COMPLEX,iarow,comm_col);

			//cout << "\n Bcast done." << endl;
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

				//cout << "\n MPI_Cart_rank done. " << endl;
                lda=ldb=loc_A.col_num;
                ldc=NB;

                if (sp==rp)
                {
                    if (coord[0]==iarow1&&coord[1]==iacol1)
                    {
                        uplo=' ';
                        pos=loc_j1;
                        zlacpy_(&uplo,&bn1,&bm1,&W1[pos],&lda,C_A1,&ldc);
                        zlacpy_(&uplo,&bn1,&bm1,&W_1[pos],&lda,C_B1,&ldc);

                        for (k=0; k<bm1; k++)
                            for (j=0; j<bn1; j++)
                            {
                                U1[(j+loc_i2)*NB+k]=conj(C_A1[k*NB+j]);
                                U_1[(j+loc_i2)*NB+k]=conj(C_B1[k*NB+j]);
                            }
                    }
                }
                else if (coord[0]==iarow1&&coord[1]==iacol1)
                {
                    uplo=' ';
                    pos=loc_j1;
                    zlacpy_(&uplo,&bn1,&bm1,&W1[pos],&lda,C_A1,&ldc);
                    zlacpy_(&uplo,&bn1,&bm1,&W_1[pos],&ldb,C_B1,&ldc);
                    MPI_Send(C_A1,NB*NB,MPI_DOUBLE_COMPLEX,rp,tag,DIAG_HPSEPS_WORLD);
                    MPI_Send(C_B1,NB*NB,MPI_DOUBLE_COMPLEX,rp,tag,DIAG_HPSEPS_WORLD);

//					ofs_running << "\n MPI_Send done. " << " tag = " << tag << endl;
                }
                else if (coord[0]==iarow2&&coord[1]==iacol2)
                {
					//ofs_running << "\n MPI_Recv prepare for tag = " << tag << endl;
                    MPI_Recv(C_A1,NB*NB,MPI_DOUBLE_COMPLEX,sp,tag,DIAG_HPSEPS_WORLD,&status1);
                    MPI_Recv(C_B1,NB*NB,MPI_DOUBLE_COMPLEX,sp,tag,DIAG_HPSEPS_WORLD,&status2);

					//ofs_running << "\n MPI_Recv done. " << endl;
                    for (k=0; k<bm1; k++)
                        for (j=0; j<bn1; j++)
                        {
                            U1[(j+loc_i2)*NB+k]=conj(C_A1[k*NB+j]);
                            U_1[(j+loc_i2)*NB+k]=conj(C_B1[k*NB+j]);
                        }
                }
            }
            MPI_Bcast(U1,NB*loc_A.row_num,MPI_DOUBLE_COMPLEX,iacol,comm_row);
            MPI_Bcast(U_1,NB*loc_A.row_num,MPI_DOUBLE_COMPLEX,iacol,comm_row);

            transa='n';
            transb='n';
            m=loc_A.row_num-loc_A.row_pos;
            n=loc_A.col_num-loc_A.col_pos;
            k=bm;
            alpha=complex<double>(-1.0,0.0);
            beta=complex<double>(1.0,0.0);
            lda=loc_A.col_num;
            ldb=NB;
            ldc=loc_A.col_num;
            pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
            zgemm_(&transa,&transb,&n,&m,&k,&alpha,&W_1[loc_A.col_pos],&lda,
                   &U_1[loc_A.row_pos*NB],&ldb,&beta,&B[pos],&ldc);
            zgemm_(&transa,&transb,&n,&m,&k,&alpha,&W1[loc_A.col_pos],&lda,
                   &U_1[loc_A.row_pos*NB],&ldb,&beta,&A[pos],&ldc);
            zgemm_(&transa,&transb,&n,&m,&k,&alpha,&W_1[loc_A.col_pos],&lda,
                   &U1[loc_A.row_pos*NB],&ldb,&beta,&A[pos],&ldc);
            if ((coord[0]==iarow))
            {
                lda=loc_A.col_num;
                ldb=loc_A.col_num;
                m=bm;
                n=loc_A.col_num-loc_A.col_pos;
                ldc=NB;
                pos=(loc_A.row_pos-bm)*loc_A.col_num+loc_A.col_pos;
                alpha=complex<double> (-0.5,0.0);
                beta=complex<double> (1.0,0.0);
                side='r';
                uplo='l';
                zhemm_(&side,&uplo,&n,&m,&alpha,C_A,&ldc,&B[pos],&ldb,&beta,&A[pos],&lda);
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
                {
                    C_A[i*NB+j]=C_B[i*NB+j]=complex<double> (0.0,0.0);

                }
            for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                for (j=0; j<NB; j++)
                {
                    U1[i*NB+j]=U_1[i*NB+j]=complex<double>(0.0,0.0);
                }
            for (i=0; i<NB; i++)
                for (j=loc_A.col_pos; j<loc_A.col_num; j++)
                {
                    W1[i*loc_A.col_num+j]=W_1[i*loc_A.col_num+j]=complex<double> (0.0,0.0);
                 }
            if (coord[1]==iacol)
            {
                if (coord[0]==iarow)
                {
                    uplo='u';
                    lda=loc_A.col_num;
                    ldb=loc_B.col_num;
                    pos=loc_i*loc_A.col_num+loc_j;
                    m=bm;
                    zpotf2_(&uplo,&m,&B[pos],&ldb,&info);
                    zhegs2_(&itype,&uplo,&m,&A[pos],&lda,&B[pos],&ldb,&info);

                    zlacpy_(&uplo,&m,&m,&A[pos],&lda,C_A,&NB);
                    zlacpy_(&uplo,&m,&m,&B[pos],&ldb,C_B,&NB);
                }
                MPI_Bcast(C_A,NB*NB,MPI_DOUBLE_COMPLEX,iarow,comm_col);
                MPI_Bcast(C_B,NB*NB,MPI_DOUBLE_COMPLEX,iarow,comm_col);
            }
            if (coord[0]==iarow)
            {
                loc_A.row_pos=loc_A.row_pos+bm;
                loc_B.row_pos=loc_B.row_pos+bm;
            }
            if (ki>=N_A-1) break;
            if ((coord[1]==iacol))
            {
                transa='c';
                transb='c';
                diag='n';
                side='l';
                uplo='u';
                m=loc_A.row_num-loc_A.row_pos;
                n=bm;
                lda=loc_A.col_num;
                ldb=loc_A.col_num;
                ldc=NB;
                pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
                alpha=complex<double> (1.0,0.0);


                //printf("(%d,%d), in pdsyg2st n2=%d \n",coord[0],coord[1],n);
                ztrsm_(&side,&uplo,&transa,&diag,&n,&m,&alpha,C_B,&ldc,&A[pos],&lda);
                ztrsm_(&side,&uplo,&transb,&diag,&n,&m,&alpha,C_B,&ldc,&B[pos],&ldb);
                alpha=complex<double>(-0.5,0.0);

                beta=complex<double>(1.0,0.0);

                zhemm_(&side,&uplo,&n,&m,&alpha,C_A,&ldc,&B[pos],&ldb,&beta,&A[pos],&lda);
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
            MPI_Bcast(U1,NB*loc_A.row_num,MPI_DOUBLE_COMPLEX,iacol,comm_row);
            MPI_Bcast(U_1,NB*loc_A.row_num,MPI_DOUBLE_COMPLEX,iacol,comm_row);
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
                        zlacpy_(&uplo,&bn1,&bm1,&U1[pos],&NB,C_A1,&NB);
                        zlacpy_(&uplo,&bn1,&bm1,&U_1[pos],&NB,C_B1,&NB);
                        for (k=0; k<bm1; k++)
                            for (j=0; j<bn1; j++)
                            {
                                W1[j*loc_A.col_num+k+loc_j2]=conj(C_A1[k*NB+j]);
                                W_1[j*loc_A.col_num+k+loc_j2]=conj(C_B1[k*NB+j]);

                            }
                    }
                }
                else if (coord[0]==iarow1&&coord[1]==iacol1)
                {
                    uplo=' ';
                    pos=loc_i1*NB;
                    zlacpy_(&uplo,&bn1,&bm1,&U1[pos],&NB,C_A1,&NB);
                    zlacpy_(&uplo,&bn1,&bm1,&U_1[pos],&NB,C_B1,&NB);
                    MPI_Send(C_A1,NB*NB,MPI_DOUBLE_COMPLEX,rp,tag,DIAG_HPSEPS_WORLD);
                    MPI_Send(C_B1,NB*NB,MPI_DOUBLE_COMPLEX,rp,tag,DIAG_HPSEPS_WORLD);
                }
                else if (coord[0]==iarow2&&coord[1]==iacol2)
                {
                    MPI_Recv(C_A1,NB*NB,MPI_DOUBLE_COMPLEX,sp,tag,DIAG_HPSEPS_WORLD,&status3);
                    MPI_Recv(C_B1,NB*NB,MPI_DOUBLE_COMPLEX,sp,tag,DIAG_HPSEPS_WORLD,&status4);
                    for (k=0; k<bm1; k++)
                        for (j=0; j<bn1; j++)
                        {
                            W1[j*loc_A.col_num+k+loc_j2]=conj(C_A1[k*NB+j]);
                            W_1[j*loc_A.col_num+k+loc_j2]=conj(C_B1[k*NB+j]);
                        }
                }
            }
            MPI_Bcast(W1,NB*loc_A.col_num,MPI_DOUBLE_COMPLEX,iarow,comm_col);
            MPI_Barrier(comm_row);
            MPI_Bcast(W_1,NB*loc_A.col_num,MPI_DOUBLE_COMPLEX,iarow,comm_col);
            MPI_Barrier(comm_row);
            transa='n';
            transb='n';
            m=loc_A.row_num-loc_A.row_pos;
            n=loc_A.col_num-loc_A.col_pos;
            k=bm;
            alpha=complex<double> (-1.0,0.0);
            beta=complex<double> (1.0,0.0);
            lda=loc_A.col_num;
            ldb=NB;
            ldc=loc_A.col_num;
            pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
            zgemm_(&transa,&transb,&n,&m,&k,&alpha,&W_1[loc_A.col_pos],&lda,
                   &U_1[loc_A.row_pos*NB],&ldb,&beta,&B[pos],&ldc);
            zgemm_(&transa,&transb,&n,&m,&k,&alpha,&W1[loc_A.col_pos],&lda,
                   &U_1[loc_A.row_pos*NB],&ldb,&beta,&A[pos],&ldc);
            zgemm_(&transa,&transb,&n,&m,&k,&alpha,&W_1[loc_A.col_pos],&lda,
                   &U1[loc_A.row_pos*NB],&ldb,&beta,&A[pos],&ldc);
            if ((coord[1]==iacol))
            {
                lda=loc_A.col_num;
                ldb=loc_A.col_num;
                ldc=NB;
                m=loc_A.row_num-loc_A.row_pos;
                n=bm;
                pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos-bm;
                alpha=complex<double> (-0.5,0.0);
                beta=complex<double> (1.0,0.0);
                side='l';
                uplo='u';
                zhemm_(&side,&uplo,&n,&m,&alpha,C_A,&ldc,&B[pos],&ldb,&beta,&A[pos],&lda);
            }
        }
//   printf("A=%lf+i%lf,%lf+i%lf\n",A[0].real(),A[0].imag(),A[1].real(),A[1].imag());
	//xiaohui add 2014-03-05, test
	MPI_Comm_free(&comm_col);
	MPI_Comm_free(&comm_row);
	//xiaohui add 2014-01-22, test
	//cout<<"pzhtrsm begin"<<endl;

	pzhtrsm(isuplo,b_n,comm_2D,NB,N_A,A,B,loc_A,loc_B);
	//xiaohui add 2014-01-22,test
	//cout<<"pzhtrsm finish"<<endl;

	//xiaohui modify 2014-03-05
	//MPI_Comm_free(&comm_col);
	//MPI_Comm_free(&comm_row);

	delete[] C_A;
	delete[] C_B;
	delete[] C_A1;
	delete[] C_B1;

}
