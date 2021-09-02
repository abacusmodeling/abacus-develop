#include "global.h"
#include "hamilt.h"
#include "diago_cg.h"
#include "diago_david.h"
#include "diago_cg_gpu.h"
#include "cufft.h"



Hamilt::Hamilt() {}
Hamilt::~Hamilt() {}


void Hamilt::diagH_pw(
    const int &istep,
    const int &iter,
    const int &ik,
    const double *precondition,
    double &avg_iter)
{
	TITLE("Hamilt","diagH_pw");
    timer::tick("Hamilt", "diagH_pw");
    double avg = 0.0;

	// set ik0 because of mem_saver.
	// if mem_saver is not used, ik0=ik, as usual.
	// but if mem_saver is used, ik0=0.
	int ik0 = ik;

	if(GlobalV::CALCULATION=="nscf" && GlobalC::wf.mem_saver==1)
	{
		if(GlobalV::BASIS_TYPE=="pw")
		{
			// generate PAOs first, then diagonalize to get
			// inital wavefunctions.
			GlobalC::wf.diago_PAO_in_pw_k2(ik, GlobalC::wf.evc[0]);
		}
#ifdef __LCAO
		else if(GlobalV::BASIS_TYPE=="lcao_in_pw")
		{
			GlobalC::wf.LCAO_in_pw_k(ik, GlobalC::wf.wanf2[0]);
		}
#endif
		ik0 = 0;
	}

    if(GlobalV::BASIS_TYPE=="lcao_in_pw")
    {
		if(GlobalV::KS_SOLVER=="lapack")
		{
			assert(GlobalV::NLOCAL >= GlobalV::NBANDS);
        	this->diagH_subspace(
				ik,
				GlobalV::NLOCAL,
				GlobalV::NBANDS,
				GlobalC::wf.wanf2[ik0],
				GlobalC::wf.evc[ik0],
				GlobalC::wf.ekb[ik]);
		}
		else
		{
			GlobalV::ofs_warning << " The diago_type " << GlobalV::KS_SOLVER
				<< " not implemented yet." << std::endl; //xiaohui add 2013-09-02
			WARNING_QUIT("Hamilt::diago","no implemt yet.");
		}
    }
    else
    {
        int ntry = 0;
        int notconv = 0;
        do
        {
	   		if(GlobalV::KS_SOLVER=="cg")
            {
				// qian change it, because it has been executed in diago_PAO_in_pw_k2
                if ( iter > 1 || istep > 1 ||  ntry > 0)
                {
                    this->diagH_subspace(
						ik,
						GlobalV::NBANDS,
						GlobalV::NBANDS,
						GlobalC::wf.evc[ik0],
						GlobalC::wf.evc[ik0],
						GlobalC::wf.ekb[ik]);

                    avg_iter += 1.0;
                }
                Diago_CG_GPU cg_gpu;
				bool reorder = true;

				CUFFT_COMPLEX *d_wf_evc;
                double *d_wf_ekb;
                int DIM_CG_GPU = GlobalC::kv.ngk[ik];
                int DIM_CG_GPU2 = GlobalC::wf.npwx * GlobalV::NPOL;
                double *d_precondition;

				if(GlobalV::NPOL==1)
				{
                    cudaMalloc((void**)&d_wf_evc, GlobalV::NBANDS * GlobalC::wf.npwx * sizeof(CUFFT_COMPLEX));
                    cudaMalloc((void**)&d_wf_ekb, GlobalV::NBANDS * sizeof(double));
                    cudaMalloc((void**)&d_precondition, DIM_CG_GPU * sizeof(double));

                    cudaMemcpy(d_wf_evc, GlobalC::wf.evc[ik0].c, GlobalV::NBANDS * GlobalC::wf.npwx * sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);
                    // cudaMemcpy(d_wf_ekb, wf.ekb[ik], NBANDS * sizeof(double), cudaMemcpyHostToDevice);
                    cudaMemcpy(d_precondition, precondition, DIM_CG_GPU * sizeof(double), cudaMemcpyHostToDevice);

                    cg_gpu.diag(d_wf_evc, d_wf_ekb, DIM_CG_GPU, GlobalC::wf.npwx,
                        GlobalV::NBANDS, d_precondition, GlobalV::ETHR,
                        GlobalV::DIAGO_CG_MAXITER, reorder, notconv, avg);

                    // to cpu
                    cudaMemcpy(GlobalC::wf.evc[ik0].c, d_wf_evc, GlobalV::NBANDS * GlobalC::wf.npwx * sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
                    cudaMemcpy(GlobalC::wf.ekb[ik], d_wf_ekb, GlobalV::NBANDS * sizeof(double), cudaMemcpyDeviceToHost);

                    cudaFree(d_wf_evc);
                    cudaFree(d_wf_ekb);
                    cudaFree(d_precondition);
				}
				else
				{
                    // to gpu
                    cudaMalloc((void**)&d_wf_evc, GlobalV::NBANDS * DIM_CG_GPU2 * sizeof(CUFFT_COMPLEX));
                    cudaMalloc((void**)&d_wf_ekb, GlobalV::NBANDS * sizeof(double));
                    cudaMalloc((void**)&d_precondition, DIM_CG_GPU2 * sizeof(double));

                    cudaMemcpy(d_wf_evc, GlobalC::wf.evc[ik0].c, GlobalV::NBANDS * DIM_CG_GPU2 * sizeof(CUFFT_COMPLEX), cudaMemcpyHostToDevice);
                    // cudaMemcpy(d_wf_ekb, GlobalC::wf.ekb[ik], GlobalV::NBANDS * sizeof(double), cudaMemcpyHostToDevice);
                    cudaMemcpy(d_precondition, precondition, DIM_CG_GPU2 * sizeof(double), cudaMemcpyHostToDevice);
                    // do things

                    cg_gpu.diag(d_wf_evc, d_wf_ekb, DIM_CG_GPU2, DIM_CG_GPU2,
                        GlobalV::NBANDS, d_precondition, GlobalV::ETHR,
                        GlobalV::DIAGO_CG_MAXITER, reorder, notconv, avg);

                    // to cpu
                    cudaMemcpy(GlobalC::wf.evc[ik0].c, d_wf_evc, GlobalV::NBANDS * DIM_CG_GPU2 * sizeof(CUFFT_COMPLEX), cudaMemcpyDeviceToHost);
                    cudaMemcpy(GlobalC::wf.ekb[ik], d_wf_ekb, GlobalV::NBANDS * sizeof(double), cudaMemcpyDeviceToHost);

                    cudaFree(d_wf_evc);
                    cudaFree(d_wf_ekb);
                    cudaFree(d_precondition);
				}
				// P.S. : nscf is the flag about reorder.
				// if diagH_subspace is done once,
				// we don't need to reorder the eigenvectors order.
				// if diagH_subspace has not been called,
				// we need to reorder the eigenvectors.
            }
	   		else if(GlobalV::KS_SOLVER=="dav")
        	{
				Diago_David david;
				if(GlobalV::NPOL==1)
				{
					david.diag(GlobalC::wf.evc[ik0], GlobalC::wf.ekb[ik], GlobalC::kv.ngk[ik],
						GlobalV::NBANDS, precondition, GlobalV::DIAGO_DAVID_NDIM,
				 		GlobalV::ETHR, GlobalV::DIAGO_CG_MAXITER, notconv, avg);
				}
				else
				{
					david.diag(GlobalC::wf.evc[ik0], GlobalC::wf.ekb[ik], GlobalC::wf.npwx*GlobalV::NPOL,
						GlobalV::NBANDS, precondition, GlobalV::DIAGO_DAVID_NDIM,
						GlobalV::ETHR, GlobalV::DIAGO_CG_MAXITER, notconv, avg);
				}
        	}
        	else
        	{
				WARNING_QUIT("calculate_bands","Check GlobalV::KS_SOLVER !");
        	}
            avg_iter += avg;
            ++ntry;
        }
        while ( this->test_exit_cond(ntry, notconv) );

        if ( notconv > max(5, GlobalV::NBANDS/4) )
        {
            std::cout << "\n notconv = " << notconv;
            std::cout << "\n Hamilt::diago', too many bands are not converged! \n";
        }
    }

	timer::tick("Hamilt","diagH_pw");
    return;
}


bool Hamilt::test_exit_cond(const int &ntry, const int &notconv)
{
    //================================================================
    // If this logical function is true, need to do diagH_subspace
	// and cg again.
    //================================================================

	bool scf = true;
	if(GlobalV::CALCULATION=="nscf") scf=false;

    // If ntry <=5, try to do it better, if ntry > 5, exit.
    const bool f1 = (ntry <= 5);

    // In non-self consistent calculation, do until totally converged.
    const bool f2 = ( (!scf && (notconv > 0)) );

    // if self consistent calculation, if not converged > 5,
    // using diagH_subspace and cg method again. ntry++
    const bool f3 = ( ( scf && (notconv > 5)) );
    return  ( f1 && ( f2 || f3 ) );
}

void Hamilt::diagH_subspace(
    const int ik,
    const int nstart,
    const int n_band,
    const ComplexMatrix &psi,
    ComplexMatrix &evc,
    double *en)
{
	if(nstart < n_band)
	{
		WARNING_QUIT("diagH_subspace","nstart < n_band!");
	}

    if(GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")
    {
        this->hpw.diagH_subspace(ik, nstart, n_band, psi, evc, en);
    }
    else
    {
		WARNING_QUIT("diagH_subspace","Check parameters: GlobalV::BASIS_TYPE. ");
    }
    return;
}



//====================================================================
// calculates eigenvalues and eigenvectors of the generalized problem
// Hv=eSv, with H hermitean matrix, S overlap matrix .
// On output both matrix are unchanged
// LAPACK version - uses both ZHEGV and ZHEGVX
//=====================================================================
void Hamilt::diagH_LAPACK(
	const int nstart,
	const int nbands,
	const ComplexMatrix &hc,
	const ComplexMatrix &sc,
	const int ldh, // nstart
	double *e,
	ComplexMatrix &hvec)
{
    TITLE("Hamilt","diagH_LAPACK");
	timer::tick("Hamilt","diagH_LAPACK");

    int lwork=0;
    //========================================
    // int ILAENV();
    // ILAENV returns optimal block size "nb"
    //========================================

    ComplexMatrix sdum(nstart, ldh);
    ComplexMatrix hdum;

    sdum = sc;

    const bool all_eigenvalues = (nstart == nbands);

    int nb = LapackConnector::ilaenv(1, "ZHETRD", "U", nstart, -1, -1, -1);
//  int nb = ILAENV(1,  "ZHETRD", "U", n, -1, -1, -1);
//  int nb = 32;

    if (nb < 1)
    {
        nb = std::max(1, nstart);
    }

	if (nb == 1 || nb >= nstart)
    {
        lwork = 2 * nstart; // mohan modify 2009-08-02
    }
    else
    {
        lwork = (nb + 1) * nstart;
    }

    std::complex<double> *work = new std::complex<double>[lwork];
	ZEROS(work, lwork);

    //=====================================================================
    // input s and (see below) h are copied so that they are not destroyed
    //=====================================================================

    int info = 0;
    int rwork_dim;
    if (all_eigenvalues)
    {
        rwork_dim = 3*nstart-2;
    }
    else
    {
        rwork_dim = 7*nstart;
    }

    double *rwork = new double[rwork_dim];
    ZEROS( rwork, rwork_dim );

    if (all_eigenvalues)
    {
        //===========================
        // calculate all eigenvalues
        //===========================
        hvec = hc;
        LapackConnector::zhegv(1, 'V', 'U', nstart, hvec , ldh, sdum, ldh, e, work , lwork , rwork, info);
    }
    else
    {
        //=====================================
        // calculate only m lowest eigenvalues
        //=====================================
        int *iwork = new int [5*nstart];
        int *ifail = new int[nstart];

        ZEROS(rwork,7*nstart);
        ZEROS(iwork,5*nstart);
        ZEROS(ifail,nstart);

        hdum.create(nstart, ldh);
        hdum = hc;

    	//=============================
    	// Number of calculated bands
    	//=============================
    	int mm = nbands;

        LapackConnector::zhegvx
        (
                1,      //INTEGER
                'V',    //CHARACTER*1
                'I',    //CHARACTER*1
                'U',    //CHARACTER*1
                nstart, //INTEGER
                hdum,   //COMPLEX*16 array
                ldh,    //INTEGER
                sdum,   //COMPLEX*16 array
                ldh,    //INTEGER
           		0.0,    //DOUBLE PRECISION
                0.0,    //DOUBLE PRECISION
                1,      //INTEGER
                nbands, //INTEGER
                0.0,    //DOUBLE PRECISION
                mm,     //INTEGER
                e,      //DOUBLE PRECISION array
                hvec,   //COMPLEX*16 array
                ldh,    //INTEGER
                work,   //DOUBLE array, dimension (MAX(1,LWORK))
                lwork,  //INTEGER
                rwork , //DOUBLE PRECISION array, dimension (7*N)
                iwork,  //INTEGER array, dimension (5*N)
                ifail,  //INTEGER array, dimension (N)
                info    //INTEGER
        );

        delete[] iwork;
        delete[] ifail;
    }
    delete[] rwork;
    delete[] work;

	timer::tick("Hamilt","diagH_LAPACK");
    return;
}
