#include "pdst2g.h"
void pdst2g(MPI_Comm comm_2D,int NB,int N_A,double *A,double *B,
            LocalMatrix loc_A,int loc_size,char isuplo)
{
    /*
     *  PSEPS routine (version 1.0) --
     *  Computer Network Information Center, CAS.
     *  December 15, 2005
     *  Purpose
     *  ======================
     *  pdst2g computes the eigenvector of a generalized symmetric eigenproblem.
     *  by computing Z= L-1W involves a triangular solve with multiple right
     *  hand sides
     *  Arguments
     *  ====================
     *  suplo    (global input) CHARACTER
     *           Specifies whether the upper or lower triangular part of the Hermitian matrix A
     *           is stored:
     *           = 'U':  Upper triangular
     *           = 'L':  Lower triangular
     *  comm2D   (input) MPI_Comm
     *           MPI 2D communicator
     *    n      (global input) integer
     *           The order of the matrix A n>=0
     *    B      (local input)double precision  pointer
     *           to point local memory to an array (loc_A.row_num* loc_A.col_num).
     *           if uplo = 'U', the elements above the first super tridiagonal is referred .
     *           if uplo = 'L',  the elements below the first subdiagonal is referred
     *   loc_A   (local input) struct Loc_A
     *           This struct avaible stores the information required to establish
     *           the mapping between an object element and its corresponding
     *           process and memory location.
     *     Z     (local input/output) double precision pointer,
     *           contains the computed eigenvectors associated with generalized
     *           Hermitian eigenproblem
     */
    int m,n,lda,ldb,ldc,loc_i,loc_j;
    int coord[2],dim[2],period[2];
    int i,j,bi,bj,k,kj,ia,ja,km,kn;
    int cur_col,cur_row;
    double alpha,beta;

    char transa,transb,diag,side,uplo;
    int iarow,iacol,pos1,pos2;
    int myid;
    int bn;
    double CU[NB][NB];

    MPI_Status status;
    MPI_Cart_get(comm_2D,2,dim,period,coord);
    loc_A.row_pos=0;
    loc_A.col_pos=0;
    bn=N_A/NB;
    if (bn*NB<N_A) bn++;
    //send B to the T  in  the diaganol
    if (isuplo=='u'||isuplo=='U')
    {
        for (bi=bn-1; bi>=0; bi--)
        {
            for (i=0; i<NB; i++)
                for (j=0; j<NB; j++)
                    CU[i][j]=0.0;
            ia=bi*NB;
            km=(ia+NB-1>N_A-1)?(N_A-1):(ia+NB-1);
            m=km-ia+1;
            for (kj=bn-1; kj>bi; kj--)
            {
                ja=kj*NB;
                kn=(ja+NB-1>N_A-1)?(N_A-1):(ja+NB-1);
                n=kn-ja+1;
                indxg2p(comm_2D,NB,ia,ja,&iarow,&iacol);
                indxg2l(ia,ja,NB,dim[0],dim[1],&loc_i,&loc_j);
                myid=iarow*dim[1]+iacol;
                if (coord[0]==iarow&&coord[1]==iacol)
                {
                    for (i=loc_i; i<loc_i+m; i++)
                        for (j=loc_j; j<loc_j+n; j++)
                            CU[i-loc_i][j-loc_j]=A[i*loc_A.col_num+j];
                    MPI_Bcast(CU,NB*NB,MPI_DOUBLE,myid,DIAG_HPSEPS_WORLD);
                }
                else
                    MPI_Bcast(CU,NB*NB,MPI_DOUBLE,myid,DIAG_HPSEPS_WORLD);

                transa='n';
                transb='n';
                k=n;
                n=loc_size;
                alpha=-1.0;
                beta=1.0;
                lda=loc_size;
                ldb=loc_size;
                ldc=NB;
                pos1=ia*loc_size;
                pos2=ja*loc_size;
                // mohan modify CU to CU[0] 2010-06-21
                dgemm_(&transa,&transb,&n,&m,&k,&alpha,&B[pos2],&lda,CU[0],&ldc,&beta,&B[pos1],&ldb);
            }
            indxg2p(comm_2D,NB,ia,ia,&iarow,&iacol);
            indxg2l(ia,ia,NB,dim[0],dim[1],&loc_i,&loc_j);

            myid=iarow*dim[1]+iacol;
            if ((coord[0]==iarow)&&(coord[1]==iacol))
            {
                for (i=loc_i; i<loc_i+m; i++)
                    for (j=loc_j; j<loc_j+m; j++)
                        CU[i-loc_i][j-loc_j]=A[i*loc_A.col_num+j];
                MPI_Bcast(CU,NB*NB,MPI_DOUBLE,myid,DIAG_HPSEPS_WORLD);
            }
            else
                MPI_Bcast(CU,NB*NB,MPI_DOUBLE,myid,DIAG_HPSEPS_WORLD);

            transb='n';
            diag='n';
            side='r';
            uplo='l';
            n=loc_size;
            pos2=ia*loc_size;

            lda=NB;
            ldb=loc_size;
            alpha=1;
            dtrsm_(&side,&uplo,&transb,&diag,&n,&m,&alpha,&CU[0][0],&lda,&B[pos2],&ldb);
        }
    }
    else
    {
        for (bi=bn-1; bi>=0; bi--)
        {
            for (i=0; i<NB; i++)
                for (j=0; j<NB; j++)
                    CU[i][j]=0.0;

            ja=bi*NB;
            kn=(ja+NB-1>N_A-1)?(N_A-1):(ja+NB-1);
            n=kn-ja+1;

            for (kj=bn-1; kj>bi; kj--)
            {
                ia=kj*NB;
                km=(ia+NB-1>N_A-1)?(N_A-1):(ia+NB-1);
                m=km-ia+1;
                indxg2p(comm_2D,NB,ia,ja,&iarow,&iacol);
                indxg2l(ia,ja,NB,dim[0],dim[1],&loc_i,&loc_j);
                myid=iarow*dim[1]+iacol;
                if (coord[0]==iarow&&coord[1]==iacol)
                {
                    for (i=loc_i; i<loc_i+m; i++)
                        for (j=loc_j; j<loc_j+n; j++)
                            CU[i-loc_i][j-loc_j]=A[i*loc_A.col_num+j];

                    MPI_Bcast(CU,NB*NB,MPI_DOUBLE,myid,DIAG_HPSEPS_WORLD);
                }
                else
                    MPI_Bcast(CU,NB*NB,MPI_DOUBLE,myid,DIAG_HPSEPS_WORLD);
                transa='n';
                transb='n';
                k=m;
                m=loc_size;
                alpha=-1.0;
                beta=1.0;
                lda=N_A;
                ldb=N_A;
                ldc=NB;
                pos1=ja;
                pos2=ia;
                dgemm_(&transa,&transb,&n,&m,&k,&alpha,CU[0],&ldc,&B[pos2],&lda,&beta,&B[pos1],&ldb);//mohan modify 2010-06-21
            }
            indxg2p(comm_2D,NB,ja,ja,&iarow,&iacol);
            indxg2l(ja,ja,NB,dim[0],dim[1],&loc_i,&loc_j);

            myid=iarow*dim[1]+iacol;
            if ((coord[0]==iarow)&&(coord[1]==iacol))
            {
                for (i=loc_i; i<loc_i+n; i++)
                    for (j=loc_j; j<loc_j+n; j++)
                        CU[i-loc_i][j-loc_j]=A[i*loc_A.col_num+j];

                MPI_Bcast(CU,NB*NB,MPI_DOUBLE,myid,DIAG_HPSEPS_WORLD);
            }
            else
                MPI_Bcast(CU,NB*NB,MPI_DOUBLE,myid,DIAG_HPSEPS_WORLD);


            transb='n';
            diag='n';
            side='l';
            uplo='u';
            m=loc_size;
            pos2=ja;

            lda=NB;
            ldb=N_A;
            alpha=1;
            dtrsm_(&side,&uplo,&transb,&diag,&n,&m,&alpha,&CU[0][0],&lda,&B[pos2],&ldb);
        }//end_for(bj)*/
    }//end_if()
}

