#include <complex>
#include <regex>
#include <fstream>
#include <cfloat>
#include <cmath>
#include <cstring>

#include <mpi.h>

#include "elpa_legacy.h"
#include "elpa_solver.h"

#include "my_math.hpp"
#include "utils.h"

int ELPA_Solver::eigenvector(double* A, double* EigenValue, double* EigenVector)
{

    int info=0;
    int success, allsuccess;
    double t;

    if(loglevel>0 && myid==0)
    {
        t=-1;
        timer(myid, "elpa_solve_evp_real_2stage", "1", t);
    }
    info=elpa_solve_evp_real_2stage(nFull, nev, A, lda, EigenValue, EigenVector, lda, nblk, nacols,
                                    mpi_comm_rows, mpi_comm_cols, comm_f,
                                    kernel_id, useQR);
    if(loglevel>0 && myid==0)
    {
        timer(myid, "elpa_solve_evp_real_2stage", "1", t);
    }
    MPI_Allreduce(&success, &allsuccess, 1, MPI_INT, MPI_MIN, comm);
    if(allsuccess == 1)
        info=0;
    else
        info=1;
    return info;
}

int ELPA_Solver::generalized_eigenvector(double* A, double* B, int& DecomposedState,
                                         double* EigenValue, double* EigenVector)
{
    int allinfo;
    int success, allsuccess;
    double t;

    if(loglevel>0 && myid==0)
    {
        t=-1;
        timer(myid, "decomposeRightMatrix", "1", t);
    }
    if(DecomposedState==0)
        allinfo=decomposeRightMatrix(B, EigenValue, EigenVector, DecomposedState);
    else
        allinfo=0;
    if(loglevel>0 && myid==0)
    {
        timer(myid, "decomposeRightMatrix", "1", t);
    }

    if(allinfo != 0)
        return allinfo;
    // transform A to A~
    if(loglevel>0 && myid==0)
    {
        t=-1;
        timer(myid, "transform A to A~", "2", t);
    }
    if(DecomposedState == 1 || DecomposedState == 2)
    {
        // calculate A*U^-1, put to work
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "A*U^-1", "2", t);
        }
        Cpdgemm('T', 'N', nFull, 1.0, A, B, 0.0, dwork.data(), desc);
        if(loglevel>1)
        {
            timer(myid, "A*U^-1", "2", t);
        }

        // calculate U^-T*(A*U^-1), put to a
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "U^-T*(A*U^-1)", "3", t);
        }
        Cpdgemm('T', 'N', nFull, 1.0, B, dwork.data(), 0.0, A, desc);
        if(loglevel>1)
        {
            timer(myid, "U^-T*(A*U^-1)", "3", t);
        }
    }
    else
    {
        // calculate B*A^T and put to work
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "B*A^T", "2", t);
        }
        Cpdgemm('N', 'T', nFull, 1.0, B, A, 0.0, dwork.data(), desc);
        if(loglevel>1)
        {
            timer(myid, "B*A^T", "2", t);
        }
        // calculate B*work^T = B*(B*A^T)^T and put to A -- origian A*x=v*B*x was transform to a*x'=v*x'
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "B*(B*A^T)^T", "3", t);
        }
        Cpdgemm('N', 'T', nFull, 1.0, B, dwork.data(), 0.0, A, desc);
        if(loglevel>1)
        {
            timer(myid, "B*(B*A^T)^T", "3", t);
        }
    }
    if(loglevel>0 && myid==0)
    {
        timer(myid, "transform A to A~", "2", t);
    }

    // calculate the eigenvalues and eigenvectors, put to ev and q
    if(loglevel>0 && myid==0)
    {
        t=-1;
        timer(myid, "elpa_solve_evp_real_2stage", "4", t);
    }
    success=elpa_solve_evp_real_2stage(nFull, nev, A, lda, EigenValue, EigenVector, lda, nblk, nacols,
                                        mpi_comm_rows, mpi_comm_cols, comm_f,
                                        kernel_id, useQR);
    if(loglevel>0 && myid==0)
    {
        timer(myid, "elpa_solve_evp_real_2stage", "4", t);
    }
    MPI_Allreduce(&success, &allsuccess, 1, MPI_INT, MPI_MIN, comm);
    if(allsuccess != 1)
        return allinfo=1;
    if(loglevel>0 && myid==0)
    {
        t=-1;
        timer(myid, "composeEigenVector", "5", t);
    }
    allinfo=composeEigenVector(DecomposedState, B, EigenVector);
    if(loglevel>0 && myid==0)
    {
        timer(myid, "composeEigenVector", "5", t);
    }
    return allinfo;
}

int ELPA_Solver::decomposeRightMatrix(double* B, double* EigenValue, double* EigenVector, int& DecomposedState)
{
    int info=0;
    int allinfo=0;
    int success;
    int allsuccess;
    double t;

    // first try cholesky decomposing
    if(nFull<CHOLESKY_CRITICAL_SIZE)
    {
        DecomposedState=1;
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "pdpotrf_", "1", t);
        }
        info=Cpdpotrf('U', nFull, B, desc);
        if(loglevel>1)
        {
            timer(myid, "pdpotrf_", "1", t);
        }
        MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MAX, comm);
        if(allinfo != 0) //pdpotrf fail, try elpa_cholesky_real
        {
            DecomposedState=2;
            if(loglevel>1)
            {
                t=-1;
                timer(myid, "elpa_cholesky_real", "2", t);
            }
            success=elpa_cholesky_real(nFull, B, narows, nblk, nacols, mpi_comm_rows, mpi_comm_cols, wantDebug);
            if(loglevel>1)
            {
                timer(myid, "elpa_cholesky_real", "2", t);
            }
            MPI_Allreduce(&success, &allsuccess, 1, MPI_INT, MPI_MIN, comm);
            if(allsuccess != 1)
                allinfo=1;
        }
    } else
    {
        DecomposedState=2;
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "elpa_cholesky_real", "1", t);
        }
        success=elpa_cholesky_real(nFull, B, narows, nblk, nacols, mpi_comm_rows, mpi_comm_cols, wantDebug);
        if(loglevel>1)
        {
            timer(myid, "elpa_cholesky_real", "1", t);
        }
        MPI_Allreduce(&success, &allsuccess, 1, MPI_INT, MPI_MIN, comm);
        if(allsuccess == 1)
            allinfo=0;
        else
            allinfo=1;
        if(allinfo != 0)
        {
            DecomposedState=1;
            if(loglevel>1)
            {
                t=-1;
                timer(myid, "pdpotrf_", "2", t);
            }
            info=Cpdpotrf('U', nFull, B, desc);
            if(loglevel>1)
            {
                timer(myid, "pdpotrf_", "2", t);
            }
            MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MAX, comm);
        }
    }

    if(allinfo==0) // calculate U^{-1}
    {
        // clear low triangle
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "clear low triangle", "1", t);
        }
        for(int j=0; j<nacols; ++j)
        {
            int jGlobal=globalIndex(j, nblk, npcols, mypcol);
            for(int i=0; i<narows; ++i)
            {
                int iGlobal=globalIndex(i, nblk, nprows, myprow);
                if(iGlobal>jGlobal) B[i+j*narows]=0;
            }
        }
        if(loglevel>1)
        {
            timer(myid, "clear low triangle", "1", t);
        }
        if(loglevel>2) saveMatrix("U.dat", nFull, B, desc, cblacs_ctxt);
        // calculate the inverse U^{-1}
        int ldb=narows;
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "invert U", "1", t);
        }
        info=elpa_invert_trm_real(nFull, B, ldb, nblk, nacols, mpi_comm_rows, mpi_comm_cols, wantDebug);
        if(loglevel>1)
        {
            timer(myid, "invert U", "1", t);
        }
        if(loglevel>2) saveMatrix("U_inv.dat", nFull, B, desc, cblacs_ctxt);
    }
    else {
    // if cholesky decomposing failed, try diagonalize
        DecomposedState=3;
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "calculate eigenvalue and eigenvector of B", "1", t);
        }
        success=elpa_solve_evp_real_2stage(nFull, nFull, B, lda, EigenValue, EigenVector, lda, nblk, nacols,
                                           mpi_comm_rows, mpi_comm_cols, comm_f,
                                           kernel_id, useQR);
        if(loglevel>1)
        {
            timer(myid, "calculate eigenvalue and eigenvector of B", "1", t);
        }
        MPI_Allreduce(&success, &allsuccess, 1, MPI_INT, MPI_MIN, comm);
        if(allsuccess != 1)
            allinfo=0;
        else
            allinfo=1;
        // calculate B^{-1/2}_{i,j}=\sum_k q_{i,k}*ev_k^{-1/2}*q_{j,k} and put to b, which will be b^-1/2
        // calculate q*ev and put to work
        for(int i=0; i<nacols; ++i)
        {
            int eidx=globalIndex(i, nblk, npcols, mypcol);
            //double ev_sqrt=1.0/sqrt(ev[eidx]);
            double ev_sqrt=EigenValue[eidx]>DBL_MIN?1.0/sqrt(EigenValue[eidx]):0;
            for(int j=0; j<narows; ++j)
                dwork[i*lda+j]=EigenVector[i*lda+j]*ev_sqrt;
        }

        // calculate qevq=qev*q^T, put to b, which is B^{-1/2}
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "qevq=qev*q^T", "2", t);
        }
        Cpdgemm('N',  'T', nFull, 1.0, dwork.data(), EigenVector, 0.0, B, desc);
        if(loglevel>1)
        {
            timer(myid, "qevq=qev*q^T", "2", t);
        }
    }
    return allinfo;
}

int ELPA_Solver::composeEigenVector(int DecomposedState, double* B, double* EigenVector)
{
    double t;
    if(DecomposedState==1 || DecomposedState==2)
    {
        // transform the eigenvectors to original general equation, let U^-1*q, and put to q
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "Cpdtrmm", "1", t);
        }
        Cpdtrmm('L', 'U', 'N', 'N', nFull, nev, 1.0, B, EigenVector, desc);
        if(loglevel>1)
        {
            timer(myid, "Cpdtrmm", "1", t);
        }
    } else
    {
        // transform the eigenvectors to original general equation, let b^T*q, and put to q
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "Cpdgemm", "1", t);
        }
        Cpdgemm('T', 'N', nFull, 1.0, B, dwork.data(), 0.0, EigenVector, desc);
        if(loglevel>1)
        {
            timer(myid, "Cpdgemm", "1", t);
        }
    }
    return 0;
}

// calculate remains of A*V - V*D
// V: eigenvector matrix
// D: Diaganal matrix of eigenvalue
// maxRemain: maximum absolute value of remains
// meanRemain: mean absolute value of remains
void ELPA_Solver::verify(double* A, double* EigenValue, double* EigenVector,
                         double &maxError, double &meanError)
{
    double* V=EigenVector;
    const int naloc=narows*nacols;
    double* D=new double[naloc];
    double* R=dwork.data();

    for(int i=0; i<naloc; ++i)
        D[i]=0;

    for(int i=0; i<nFull; ++i)
    {
        int localRow, localCol;
        int localProcRow, localProcCol;

        localRow=localIndex(i, nblk, nprows, localProcRow);
        if(myprow==localProcRow)
        {
            localCol=localIndex(i, nblk, npcols, localProcCol);
            if(mypcol==localProcCol)
            {
                int idx = localRow + localCol*narows;
                D[idx]=EigenValue[i];
            }
        }
    }

    // R=V*D
    Cpdsymm('R', 'U', nFull, 1.0, D, V, 0.0, R, desc);
    // R=A*V-V*D=A*V-R
    Cpdsymm('L', 'U', nFull, 1.0, A, V, -1.0, R, desc);
    // calculate the maximum and mean value of sum_i{R(:,i)*R(:,i)}
    double sumError=0;
    maxError=0;
    for(int i=1; i<=nev; ++i)
    {
        double E;
        Cpddot(nFull, E, R, 1, i, 1,
                         R, 1, i, 1, desc);
        sumError+=E;
        maxError=std::max(E, maxError);
    }
    meanError=sumError/nFull;
    // global mean and max Error
    delete[] D;
}


// calculate the computational error
// $ \ket{ \delta \psi_i } = (H - \epsilon_i S)\ket{\psi_i} $
// $ \delta_i = \braket{ \delta \psi_i | \delta \psi_i } $
//
// V: eigenvector matrix
// D: Diaganal matrix of eigenvalue
// maxError: maximum absolute value of error
// meanError: mean absolute value of error
void ELPA_Solver::verify(double* A, double* B, double* EigenValue, double* EigenVector,
                        double &maxError, double &meanError)
{
    double* V=EigenVector;
    const int naloc=narows*nacols;
    double* D=new double[naloc];
    double* R=new double[naloc];

    for(int i=0; i<naloc; ++i)
        D[i]=0;

    for(int i=0; i<nFull; ++i)
    {
        int localRow, localCol;
        int localProcRow, localProcCol;

        localRow=localIndex(i, nblk, nprows, localProcRow);
        if(myprow==localProcRow)
        {
            localCol=localIndex(i, nblk, npcols, localProcCol);
            if(mypcol==localProcCol)
            {
                int idx = localRow + localCol*narows;
                D[idx]=EigenValue[i];
            }
        }
    }

    // dwork=B*V
    Cpdsymm('L', 'U', nFull, 1.0, B, V, 0.0, dwork.data(), desc);
    // R=B*V*D=dwork*D
    Cpdsymm('R', 'U', nFull, 1.0, D, dwork.data(), 0.0, R, desc);
    // R=A*V-B*V*D=A*V-R
    Cpdsymm('L', 'U', nFull, 1.0, A, V, -1.0, R, desc);
    // calculate the maximum and mean value of sum_i{R(:,i)*R(:,i)}
    double sumError=0;
    maxError=0;
    for(int i=1; i<=nev; ++i)
    {
        double E;
        Cpddot(nFull, E, R, 1, i, 1,
                         R, 1, i, 1, desc);
        sumError+=E;
        maxError=maxError>E?maxError:E;
    }
    meanError=sumError/nFull;

    delete[] D;
    delete[] R;
}
