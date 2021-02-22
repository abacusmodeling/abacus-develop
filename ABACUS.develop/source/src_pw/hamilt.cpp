#include "global.h"
#include "hamilt.h"
#include "diago_cg.h"
#include "diago_david.h"

Hamilt::Hamilt() {}
Hamilt::~Hamilt() 
{
	if(test_deconstructor)
	{
		cout << " ~Hamilt()" << endl;	
	}
}


void Hamilt::init_k(const int ik)
{
    this->hpw.init_k(ik);
}


void Hamilt::diago(
    const int &istep,
    const int &iter,
    const int &ik,
    const double *precondition,
    double &avg_iter)
{
	TITLE("Hamilt","diago");
	timer::tick("Hamilt","diago",'F');
//	timer::tick("Hamilt","diago");
    double avg = 0.0;

	// set ik0 because of mem_saver.
	// if mem_saver is not used, ik0=ik, as usual.
	// but if mem_saver is used, ik0=0.
	int ik0 = ik;

	// mohan add 2010-09-07
	// mohan update 2010-09-14, about wanf2
	if(CALCULATION=="nscf" && wf.mem_saver==1)
	{
		//if(LOCAL_BASIS==0) xiaohui modify 2013-09-02
		if(BASIS_TYPE=="pw") //xiaohui add 2013-09-02
		{
			// generate PAOs first, then diagonalize to get
			// inital wavefunctions.
			wf.diago_PAO_in_pw_k2(ik, wf.evc[0]);	
		}
		//else if(LOCAL_BASIS==4) xiaohui modify 2013-09-02
		else if(BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02. Attention! Maybe "BASIS_TYPE==lcao_in_pw" is just ok.
		{
			wf.LCAO_in_pw_k(ik, wf.wanf2[0]);
		}
		ik0 = 0;
	}

    //if (LOCAL_BASIS) xiaohui modify 2013-09-02
    if(BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02. Attention...
    {
		//if(DIAGO_TYPE=="lapack") xiaohui modify 2013-09-02
		if(KS_SOLVER=="lapack") //xiaohui add 2013-09-02
		{
			assert(NLOCAL >= NBANDS);
        	this->cinitcgg(ik, NLOCAL, NBANDS, wf.wanf2[ik0], wf.evc[ik0], wf.ekb[ik]);

			/*
			cout << " Check evc unit of ik = " << ik0 << endl;
			for(int ib=0; ib<NBANDS; ib++)
			{
				complex<double> tmpn = ZERO;
				for(int ig=0; ig<kv.ngk[ik0]; ig++)
				{
					tmpn += conj( wf.evc[ik0](ib,ig) ) * wf.evc[ik0](ib,ig);
				}
				cout << " tmpn = " << tmpn << endl;
			}
			*/

		}
		else
		{
			//ofs_warning << " The diago_type " << DIAGO_TYPE << " not implemented yet." << endl; xiaohui modify 2013-09-02
			ofs_warning << " The diago_type " << KS_SOLVER << " not implemented yet." << endl; //xiaohui add 2013-09-02
			WARNING_QUIT("Hamilt::diago","no implemt yet.");
		}
		/*
		for(int ik=0; ik<kv.nks; ik++)
		{
			for(int ib=0; ib<NBANDS; ib++)
			{
				cout << " " << wf.ekb[ik][ib];
			}
			cout << endl;
		}
		*/
    }
    else
    {
        int ntry = 0;
        int notconv = 0;
        do
        {	
	   		//if ( DIAGO_TYPE=="cg" ) xiaohui modify 2013-09-02
	   		if(KS_SOLVER=="cg") //xiaohui add 2013-09-02
            {			
                if ( iter > 0 || istep > 0 ||  ntry > 0)
                {
                    this->cinitcgg( ik ,NBANDS, NBANDS, wf.evc[ik0], wf.evc[ik0], wf.ekb[ik]);
                    avg_iter += 1.0;
                }
                Diago_CG cg;
    	
				bool reorder = true;
				if(NPOL==1) cg.diag(wf.evc[ik0], wf.ekb[ik], kv.ngk[ik],
						NBANDS, precondition, ETHR,
						DIAGO_CG_MAXITER, reorder, notconv, avg);
				else{
					cg.diag(wf.evc[ik0], wf.ekb[ik], wf.npwx*NPOL,
						NBANDS, precondition, ETHR,
						DIAGO_CG_MAXITER, reorder, notconv, avg);
				}
				// P.S. : nscf is the flag about reorder.
				// if cinitcgg is done once,
				// we don't need to reorder the eigenvectors order.
				// if cinitcgg has not been called,
				// we need to reorder the eigenvectors.
            }
            //else if ( DIAGO_TYPE=="dav" ) xiaohui modify 2013-09-02
	   		else if(KS_SOLVER=="dav") //xiaohui add 2013-09-02
        	{
				Diago_David david;
				if(NPOL==1) david.diag(wf.evc[ik0], wf.ekb[ik], kv.ngk[ik],
						NBANDS, precondition, DIAGO_DAVID_NDIM,
						ETHR, DIAGO_CG_MAXITER, notconv, avg);
				else{
					david.diag(wf.evc[ik0], wf.ekb[ik], wf.npwx*NPOL,
						NBANDS, precondition, DIAGO_DAVID_NDIM,
						ETHR, DIAGO_CG_MAXITER, notconv, avg);
				}
        	}
        	else
        	{
                //WARNING_QUIT("calculate_bands","Check DIAGO_TYPE !"); xiaohui modify 2013-09-02
				WARNING_QUIT("calculate_bands","Check KS_SOLVER !"); //xiaohui add 2013-09-02
        	}
            avg_iter += avg;
            ++ntry;
        }
        while ( this->test_exit_cond(ntry, notconv) );

        if ( notconv > max(5, NBANDS/4) )
        {
            cout << "\n notconv = " << notconv;
            cout << "\n Hamilt::diago', too many bands are not converged! \n";
        }
    }

	timer::tick("Hamilt","diago",'F');
    return;
}

bool Hamilt::test_exit_cond(const int &ntry, const int &notconv)
{
    //================================================================
    // If this logical function is true, need to do cinitcgg and cg
    // again.
    //================================================================

	bool scf = true;
	if(CALCULATION=="nscf") scf=false;

    // If ntry <=5, try to do it better, if ntry > 5, exit.
    const bool f1 = (ntry <= 5);

    // In non-self consistent calculation, do until totally converged.
    const bool f2 = ( (!scf && (notconv > 0)) );

    // if self consistent calculation, if not converged > 5,
    // using cinitcgg and cg method again. ntry++
    const bool f3 = ( ( scf && (notconv > 5)) );
    return  ( f1 && ( f2 || f3 ) );
}

void Hamilt::cinitcgg(
    const int ik,
    const int nstart,
    const int n_band,
    const ComplexMatrix &psi,
    ComplexMatrix &evc,
    double *en)
{
    // 1 : Local basis or not.
    //     yes : using cinitcgg
    //     no  : using cinitcgg, then cg method
    //
    // 2 : Linear scaling or not.
    //     yes : using
    //
    // 3 : Lapack using or cg method.
    //
    //
    // 4 : Sparse matrix or not.
    // 	   (Default setting : yes for LINEAR SCALING)
    //
    //     LINEAR SCALING: yes( consider set it or not ), no (never do it)
    //     LOCAL BASIS : yes( consider set it or not), no(never do it)
    //
    //     belong to STORE hamiltonian and overlap.
    //     Can used for LAPACK for check.
    //     Formally used for CG method.
    //
    //	   In LINEAR SCALING case,
    //	   there are some difference between
    //     S|psi and H|psi
    //
    // 5 : belong to CONSTRUCT hamiltonian and overlap.
    // 	   only for Local basis = 1 or 2.
    //     not for Local basis = 0
    //     can used for both cg or lapack
    //     can used for both sparsematrix or complexmatrix
    //
    //
    //     Using overlap from SIESTA method
    //     Using kinetic from SISTA method
    //     Using non-local from SISTA method
    //     Local part from SISTA method
    //
    //     Overlap, kinetic and non-local can only do once.
    //     Local part built up in each iteration

	if(nstart < n_band)
	{
		WARNING_QUIT("cinitcgg","nstart < n_band!");
	}

    //if (LINEAR_SCALING && LOCAL_BASIS) xiaohui modify 2013-09-02
    if(BASIS_TYPE=="lcao") //xiaohui add 2013-09-02
    {
        WARNING_QUIT("cinitcgg","Diago not doing here.");
    }
    //else if (!LINEAR_SCALING) xiaohui modify 2013-09-02
    else if(BASIS_TYPE=="pw" || BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02
    {
        this->hpw.cinitcgg(ik, nstart, n_band, psi, evc, en);
    }
    else
    {
        //WARNING_QUIT("cinitcgg","Check parameters: LINEAR_SCALING and LOCAL_BASIS. "); xiaohui modify 2013-09-02
	WARNING_QUIT("cinitcgg","Check parameters: BASIS_TYPE. "); //xiaohui add 2013-09-02
    }
    return;
}


void Hamilt::h_1psi(const int dim,const complex<double> *psi,complex<double> *hpsi,complex<double> *spsi)
{
    //if (LINEAR_SCALING) xiaohui modify 2013-09-02
    if(BASIS_TYPE=="lcao") //xiaohui add 2013-09-02. Attention! Maybe this change is not reasonable.
    {
        timer::tick("Hamilt","Sparse_SH");
        double *psi_real = new double[dim];
        double *hpsi_real = new double[dim];
        double *spsi_real = new double[dim];

        for (int i=0; i<dim; i++) psi_real[i] = psi[i].real();
        for (int i=0; i<dim; i++) hpsi_real[i] = hpsi[i].real();
        for (int i=0; i<dim; i++) spsi_real[i] = spsi[i].real();

		WARNING_QUIT("no sparse anymore","haha");
//      this->hon.UHM.Sparse_H.multiply_vector(psi_real, hpsi_real);
//      this->hon.UOM.Sparse_S.multiply_vector(psi_real, spsi_real);

        for (int i=0; i<dim; i++) spsi[i] = complex<double>( spsi_real[i], 0.0 );
        for (int i=0; i<dim; i++) hpsi[i] = complex<double>( hpsi_real[i], 0.0 );

        delete[] psi_real;
        delete[] hpsi_real;
        delete[] spsi_real;
        timer::tick("Hamilt","Sparse_SH");
    }
    else
    {
        this->hpw.h_1psi(dim, psi, hpsi, spsi);
    }
    return;
}

void Hamilt::s_1psi(const int dim, const complex<double> *psi, complex<double> *spsi)
{
    //if (LINEAR_SCALING) xiaohui modify 2013-09-02
    if(BASIS_TYPE=="lcao") //xiaohui add 2013-09-02. Attention! Maybe this change is not reasonable.
    {
        timer::tick("Hamilt","Sparse_S");
        double *psi_real = new double[dim];
        double *spsi_real = new double[dim];

        for (int i=0; i<dim; i++) psi_real[i] = psi[i].real();
        for (int i=0; i<dim; i++) spsi_real[i] = spsi[i].real();

		WARNING_QUIT("no sparse anymore","haha");
        //this->hon.UOM.Sparse_S.multiply_vector(psi_real, spsi_real);

        for (int i=0; i<dim; i++)
        {
            spsi[i] = complex<double>( spsi_real[i] , 0.0 );
        }

        delete[] psi_real;
        delete[] spsi_real;
        timer::tick("Hamilt","Sparse_S");
    }
    else
    {
        this->hpw.s_1psi(dim, psi, spsi);
    }
    return;
}

void Hamilt::h_psi( const int dim, const complex<double> *psi, complex<double> *hpsi)
{
    //if (LINEAR_SCALING) xiaohui modify 2013-09-02
    if(BASIS_TYPE=="lcao") //xiaohui add 2013-09-02. Attention! Maybe this change is not reasonable.
    {
        timer::tick("Hamilt","Sparse_H");
        double *psi_real = new double[dim];
        double *hpsi_real = new double[dim];

        for (int i=0; i<dim; i++) psi_real[i] = psi[i].real();
        for (int i=0; i<dim; i++) hpsi_real[i] = hpsi[i].real();

		WARNING_QUIT("no sparse anymore","haha");
        //this->hon.UHM.Sparse_H.multiply_vector(psi_real, hpsi_real);

        for (int i=0; i<dim; i++) hpsi[i] = complex<double>( hpsi[i].real(), 0.0 );

        delete[] psi_real;
        delete[] hpsi_real;
        timer::tick("Hamilt","Sparse_H");
    }
    else
    {
        this->hpw.h_psi( psi, hpsi);
    }
    return;
}

void Hamilt::cdiaghg(
	const int nstart,
	const int nbands,
	const ComplexMatrix &hc,
	const ComplexMatrix &sc,
	const int ldh, // nstart
	double *e,
	ComplexMatrix &hvec)
{
    TITLE("Hamilt","cdiaghg");
	timer::tick("Hamilt","cdiaghg");
    //====================================================================
    // calculates eigenvalues and eigenvectors of the generalized problem
    // Hv=eSv, with H hermitean matrix, S overlap matrix .
    // On output both matrix are unchanged
    // LAPACK version - uses both ZHEGV and ZHEGVX
    //=====================================================================
//	BLOCK_HERE("cdiaghg");

    int lwork;
    //========================================
    // int ILAENV();
    // ILAENV returns optimal block size "nb"
    //========================================

    ComplexMatrix sdum(nstart, ldh);
    ComplexMatrix hdum;
//	cout << "\n nstart = " << nstart << endl;
//	cout << "\n ldh = " << ldh  << endl;
//	cout << "\n nbands = " << nbands << endl;

    sdum = sc;

    const bool all_eigenvalues = (nstart == nbands);
//	cout << "\n all eigenvalue = " << all_eigenvalues << endl;

    //  cout<<"\n dimension = "<<n<<endl;
    int nb = LapackConnector::ilaenv(1, "ZHETRD", "U", nstart, -1, -1, -1);
//  int nb = ILAENV(1,  "ZHETRD", "U", n, -1, -1, -1);
//  int nb = 32;
//	cout<<"\n nb = "<<nb;

    if (nb < 1)
    {
        nb = std::max(1, nstart);
    }
    
	if (nb == 1 || nb >= nstart)
    {
		// lwork = 2 * nstart - 1; // wrong!!!!
        lwork = 2 * nstart; // chenmohan modify 2009-08-02
    }
    else
    {
        lwork = (nb + 1) * nstart;
    }

    complex<double> *work = new complex<double>[lwork];
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
//		cout << "\n ZHEGV";
        //===========================
        // calculate all eigenvalues
        //===========================
        hvec = hc;
        LapackConnector::zhegv(1, 'V', 'U', nstart, hvec , ldh, sdum, ldh, e, work , lwork , rwork, info);
    }
    else
    {
//		cout << "\n ZHEGVXXXXXXXXXXXXXXXXXX !!!";
        //=====================================
        // calculate only m lowest eigenvalues
        //=====================================
        double *work2 = new double[2*lwork];
        int *iwork = new int [5*nstart];
        int *ifail = new int[nstart];

        ZEROS(work2,2*lwork);
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

//		cout << "\n pw.gcar[0] = " << pw.gcar[0].x << " " << pw.gcar[0].y << " " << pw.gcar[0].z << endl; 

//		cout << "\n work[0] = " << work[0];
//		cout << "\n info = " << info; 
		
        //ZHEGVX(1, 'V', 'I', 'U', nstart, hdum, ldh, sdum, ldh,
          //     0.0, 0.0, 1, nbands, 0.0, mm, e, hvec, ldh,
            //   work2, lwork, rwork, iwork, ifail, info);
		
        delete[] work2;
        delete[] iwork;
        delete[] ifail;
    }
    delete[] rwork;
    delete[] work;

	timer::tick("Hamilt","cdiaghg");
    return;
}
