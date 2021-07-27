#include "numerical_basis.h"
#include "../src_pw/global.h"
#include "../module_symmetry/symmetry.h"
#include "winput.h"
#include "../module_base/math_ylmreal.h"

Numerical_Basis::Numerical_Basis() {}
Numerical_Basis::~Numerical_Basis() {}

//============================================================
// MEMBER FUNCTION :
// NAME : init
// DESCRIPTION : Two main functions:
// (1) start_from_file = true;
// Firstly, use check(1) to call bessel_basis.init
// to generate TableOne.
// Secondly readin C4 from file.
// Thirdly generate 3D atomic wfc in G space, put the
// results in psi.
//
// (2) If output overlap Q, start_from_file = false;
// Firstly, use check(0) to call bessel_basis,init
// to generate TableOne
// Secondly output overlap, use psi(evc) and jlq3d.
//============================================================
void Numerical_Basis::start_from_file_k( const int &ik, ComplexMatrix &psi)
{
    TITLE("Numerical_Basis","start_from_file_k");

    if (!this->init_label)
    {
        // 1 stands for : start_from_file
        this->bessel_basis.init( 1, pw.ecutwfc, ucell.ntype, ucell.lmax );
        this->mu_index = this->init_mu_index();
        this->init_label = true;
    }
    this->numerical_atomic_wfc(ik, kv.ngk[ik], psi);
}

// The function is called in run_fp.cpp.
void Numerical_Basis::output_overlap( const ComplexMatrix *psi)
{
    TITLE("Numerical_Basis","output_overlap");
    NEW_PART("Overlap Data For Spillage Minimization");
    
	//---------------------------------------------------------
	// if the numerical_basis hasn't been initialized yet,
    // then we initial here.
	//---------------------------------------------------------
    if (!this->init_label)
    {
        // 0 stands for : 'Faln' is not used.
        this->bessel_basis.init( 0, pw.ecutwfc, ucell.ntype, ucell.lmax );
        this->mu_index = this->init_mu_index();
        this->init_label = true;
    }

    const string command0 =  "test -d " + winput::spillage_outdir + " || mkdir " + winput::spillage_outdir;
    if(MY_RANK==0)
        system( command0.c_str() );

	for(int derivative_order=0; derivative_order<=1; ++derivative_order)            // Peize Lin add 2020.04.23
	{
        ofstream ofs;
        stringstream ss;
        // the parameter 'winput::spillage_outdir' is read from INPUTw.
        ss << winput::spillage_outdir << "/" << ucell.latName << "." << derivative_order << ".dat";
        if (MY_RANK==0)
        {
            ofs.open(ss.str().c_str());
        }

        // OVERLAP : < J_mu | Psi >
        std::vector<ComplexArray> overlap_Q(kv.nks);

        // OVERLAP : < J_mu | J_nu >
        std::vector<ComplexArray> overlap_Sq(kv.nks);

        OUT(ofs_running,"number of k points", kv.nks);
        OUT(ofs_running,"number of bands", NBANDS);
        OUT(ofs_running,"number of local orbitals", NLOCAL);
        OUT(ofs_running,"number of eigenvalues of Jl(x)", this->bessel_basis.get_ecut_number());

        // nks now is the reduced k-points.
        for (int ik=0; ik<kv.nks; ik++)
        {
            const int npw= kv.ngk[ik];
            ofs_running << " --------------------------------------------------------" << endl;
            ofs_running << " Print the overlap matrixs Q and S for this kpoint" << endl;
            ofs_running << setw(8) << "ik" << setw(8) << "npw" << endl;
            ofs_running << setw(8) << ik+1 << setw(8) << npw << endl;
            ofs_running << " --------------------------------------------------------" << endl;

            // search for all k-points.
            overlap_Q[ik] = this->cal_overlap_Q(ik, npw, psi[ik], derivative_order);
            DONE(ofs_running,"cal_overlap_Q");

            // (2) generate Sq matrix if necessary.
            if (winput::out_spillage == 2)
            {
                overlap_Sq[ik] = this->cal_overlap_Sq( ik, npw, derivative_order );
                DONE(ofs_running,"cal_overlap_Sq");
            }
        }

        const matrix overlap_V = this->cal_overlap_V(psi, derivative_order);		// Peize Lin add 2020.04.23

    #ifdef __MPI
        for (int ik=0; ik<kv.nks; ik++)
        {
            Parallel_Reduce::reduce_complex_double_pool( overlap_Q[ik].ptr, overlap_Q[ik].getSize() );
            Parallel_Reduce::reduce_complex_double_pool( overlap_Sq[ik].ptr, overlap_Sq[ik].getSize() );
        }
        Parallel_Reduce::reduce_double_pool(overlap_V.c, overlap_V.nr*overlap_V.nc);		// Peize Lin add 2020.04.23
    #endif

        this->output_info( ofs, bessel_basis );
        
        this->output_k( ofs );

        this->output_overlap_Q( ofs, overlap_Q );

        if (winput::out_spillage == 2)
        {
            this->output_overlap_Sq(ss.str(), ofs, overlap_Sq);
        }

        this->output_overlap_V(ofs, overlap_V);		// Peize Lin add 2020.04.23

        if (MY_RANK==0) ofs.close();
    }
    return;
}

ComplexArray Numerical_Basis::cal_overlap_Q(
    const int &ik,
    const int &np,
    const ComplexMatrix &psi,
	const int derivative_order) const
{
    TITLE("Numerical_Basis","cal_overlap_Q");
    timer::tick("Numerical_Basis","cal_overlap_Q");

	ofs_running << " OUTPUT THE OVERLAP BETWEEN SPHERICAL BESSEL FUNCTIONS AND BLOCH WAVE FUNCTIONS" << endl;
	ofs_running << " Q = < J_mu, q | Psi_n, k > " << endl;

    ComplexArray overlap_Q(NBANDS, NLOCAL, this->bessel_basis.get_ecut_number() );
    overlap_Q.zero_out();

	const double normalization = (4 * PI) / sqrt(ucell.omega);			// Peize Lin add normalization 2015-12-29

    std::vector<Vector3<double>> gk(np);
    for (int ig=0; ig<np; ig++)
        gk[ig] = wf.get_1qvec_cartesian(ik, ig);

	const realArray flq = this->cal_flq(ik, gk);

    const matrix ylm = Numerical_Basis::cal_ylm(gk);

    ofs_running << "\n " << setw(5)
        << "ik" << setw(8) 
        << "Type1" << setw(8) 
        << "Atom1" << setw(8) 
        << "L" << endl;

    for (int T = 0; T < ucell.ntype; T++)
    {
        //OUT("T",T);
        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            //OUT("I",I);
            complex<double> *sk = wf.get_sk(ik, T, I);
            for (int L=0; L< ucell.atoms[T].nwl+1; L++)
            {
                ofs_running << " " << setw(5) << ik+1
                            << setw(8) << ucell.atoms[T].label
                            << setw(8) << I+1 
							<< setw(8) << L
							<< endl;
                //OUT("l",l);
                complex<double> lphase = normalization * pow(IMAG_UNIT, L);			// Peize Lin add normalization 2015-12-29
                for (int ie=0; ie < this->bessel_basis.get_ecut_number(); ie++)
                {
                    const int N = 0;
                    assert( ucell.nmax == 1);
                    for (int m=0; m<2*L+1; m++)
                    {
                        const int lm = L*L+m;
                        for (int ib=0; ib<NBANDS; ib++)
                        {
                            complex<double> overlap_tmp = ZERO;
                            for (int ig=0; ig<np; ig++)
                            {
//                              const complex<double> local_tmp = lphase * sk[ig] * ylm(lm, ig) * flq[ig];
                                const complex<double> local_tmp = lphase * sk[ig] * ylm(lm, ig) * flq(L,ie,ig) * pow(gk[ig].norm2(),derivative_order);		// Peize Lin add for dpsi 2020.04.23
                                overlap_tmp += conj( local_tmp ) * psi(ib, ig); // psi is bloch orbitals
                            }
                            overlap_Q(ib, this->mu_index[T](I, L, N, m), ie) = overlap_tmp;
                        }
                    }
                }//end ie
            }//end l
            delete[] sk;    sk=nullptr;
        }
    }

    timer::tick("Numerical_Basis","cal_overlap_Q");
    return overlap_Q;
}

ComplexArray Numerical_Basis::cal_overlap_Sq(
    const int &ik,
    const int &np,
	const int derivative_order) const
{
    TITLE("Numerical_Basis","cal_overlap_Sq");
    timer::tick("Numerical_Basis","cal_overlap_Sq");

	ofs_running << " OUTPUT THE OVERLAP BETWEEN SPHERICAL BESSEL FUNCTIONS"  << endl;
	ofs_running << " S = < J_mu,q1 | J_nu,q2 >" << endl; 

    const int enumber = this->bessel_basis.get_ecut_number();
    ComplexArray overlap_Sq( NLOCAL, NLOCAL, enumber, enumber );
    overlap_Sq.zero_out();

	const double normalization = (4 * PI) * (4 * PI) / ucell.omega;			// Peize Lin add normalization 2015-12-29
	
    std::vector<Vector3<double>> gk(np);
    for (int ig=0; ig<np; ig++)
        gk[ig] = wf.get_1qvec_cartesian(ik, ig);

	const realArray flq = this->cal_flq(ik, gk);

    const matrix ylm = Numerical_Basis::cal_ylm(gk);

    ofs_running << "\n " << setw(5)
        << "ik" << setw(8) 
        << "Type1" << setw(8) 
        << "Atom1" << setw(8) 
        << "L1" << setw(8) 
        << "Type2" << setw(8) 
        << "Atom2" << setw(8) 
        << "L2" << endl;

    for (int T1 = 0; T1 < ucell.ntype; T1++) // 1.1
    {
        for (int I1 = 0; I1 < ucell.atoms[T1].na; I1++) // 1.2
        {
            complex<double> *sk1 = wf.get_sk(ik, T1, I1);
            for (int T2=0; T2<ucell.ntype; T2++) // 2.1
            {
                for (int I2=0; I2<ucell.atoms[T2].na; I2++) // 2.2
                {
                    complex<double> *sk2 = wf.get_sk(ik, T2, I2);
                    for (int l1 = 0; l1 < ucell.atoms[T1].nwl+1; l1++) // 1.3
                    {
                        const complex<double> lphase1 = normalization * pow(IMAG_UNIT, l1);			// Peize Lin add normalization 2015-12-29
                        for (int l2 = 0; l2 < ucell.atoms[T2].nwl+1; l2++) // 2.3
                        {
                            ofs_running << " " << setw(5)
                                << ik+1 << setw(8)
                                << ucell.atoms[T1].label << setw(8)
                                << I1+1 << setw(8)
                                << l1 << setw(8)
                                << ucell.atoms[T2].label << setw(8)
                                << I2+1 << setw(8)
                                << l2 << setw(8) << endl;

                            const complex<double> lphase2 = pow(IMAG_UNIT, l2);
                            for (int ic1=0; ic1 < ucell.nmax; ic1++) // 1.5
                            {
                                for (int ic2=0; ic2 < ucell.nmax; ic2++) // 2.5
                                {
                                    for (int m1=0; m1<2*l1+1; m1++) // 1.6
                                    {
                                        const int lm1 = l1*l1+m1;
                                        std::vector<std::complex<double>> about_ig(np, std::complex<double>(0.0,0.0));
                                        for (int ig=0; ig<np; ig++)
                                        {
//                                          about_ig[ig] = conj( lphase1 * sk1[ig] * ylm(lm1, ig) );
                                            about_ig[ig] = conj( lphase1 * sk1[ig] * ylm(lm1, ig) ) * pow(gk[ig].norm2(),derivative_order);		// Peize Lin add for dpsi 2020.04.23
                                        }
                                        for (int m2=0; m2<2*l2+1; m2++) // 2.6
                                        {
                                            const int lm2 = l2*l2+m2;
                                            const int iwt1 = this->mu_index[T1](I1,l1,ic1,m1);
                                            const int iwt2 = this->mu_index[T2](I2,l2,ic2,m2);
                                            for (int ig=0; ig<np; ig++)
                                            {
                                                const complex<double> about_ig3= lphase2 * sk2[ig] * ylm(lm2, ig)
                                                                                 * about_ig[ig];

                                                for (int ie1=0; ie1 < enumber; ie1++) // 1.4
                                                {
                                                    for (int ie2=0; ie2 < enumber; ie2++) // 2.4
                                                    {
                                                        overlap_Sq( iwt1, iwt2, ie1, ie2) += 
                                                            about_ig3 * flq(l1,ie1,ig) * flq(l2,ie2,ig);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    delete[] sk2;   sk2=nullptr;
                }
            }
            delete[] sk1;   sk1=nullptr;
        }
    }

    timer::tick("Numerical_Basis","cal_overlap_Sq");
    return overlap_Sq;
}

// Peize Lin add for dpsi 2020.04.23
matrix Numerical_Basis::cal_overlap_V(
	const ComplexMatrix *psi,
	const int derivative_order)
{
	matrix overlap_V(kv.nks, NBANDS);
	for(int ik=0; ik<kv.nks; ++ik)
	{
		for(int ib=0; ib<NBANDS; ++ib)
		{
			for(int ig=0; ig<kv.ngk[ik]; ++ig)
			{
				overlap_V(ik,ib)+= norm(psi[ik](ib,ig)) * pow(wf.get_1qvec_cartesian(ik,ig).norm2(),derivative_order);
			}
		}
	}
	return overlap_V;
}

realArray Numerical_Basis::cal_flq(const int ik, const std::vector<Vector3<double>> &gk) const
{
	const int np = gk.size();
	const int enumber = this->bessel_basis.get_ecut_number();

    // get flq(G) = \int f(r)jl(G*r) from interpolation table.
    realArray flq(ucell.lmax+1, enumber, np);
    for (int il=0; il<ucell.lmax+1; il++)
        for (int ie=0; ie<enumber; ie++)
            for (int ig=0; ig<np; ig++)
                flq(il,ie,ig) = this->bessel_basis.Polynomial_Interpolation2(il, ie, gk[ig].norm() * ucell.tpiba );
	return flq;	
}

matrix Numerical_Basis::cal_ylm(const std::vector<Vector3<double>> &gk)
{
    const int total_lm = ( ucell.lmax + 1) * ( ucell.lmax + 1);
    matrix ylm(total_lm, gk.size());
    YlmReal::Ylm_Real(total_lm, gk.size(), gk.data(), ylm);
    return ylm;
}

std::vector<IntArray> Numerical_Basis::init_mu_index(void)
{
	ofs_running << " Initialize the mu index" << endl;
    std::vector<IntArray> mu_index_(ucell.ntype);

    int mu = 0;
    for (int it=0; it<ucell.ntype; it++)
    {
        mu_index_[it].create(
            ucell.atoms[it].na,
            ucell.atoms[it].nwl+1,
            ucell.nmax,
            2*(ucell.atoms[it].nwl+1)+1); // m ==> 2*l+1

		// mohan added 2021-01-03
		ofs_running << "Type " << it+1 
            << " number_of_atoms " << ucell.atoms[it].na
            << " number_of_L " << ucell.atoms[it].nwl+1
            << " number_of_n " << ucell.nmax
            << " number_of_m " << 2*(ucell.atoms[it].nwl+1)+1 << endl;

        for (int ia=0; ia<ucell.atoms[it].na; ia++)
        {
            for (int l=0; l< ucell.atoms[it].nwl+1; l++)
            {
                for (int n=0; n< ucell.atoms[it].l_nchi[l]; n++)
                {
                    for (int m=0; m<2*l+1; m++)
                    {
                        mu_index_[it](ia,l,n,m) = mu;
                        mu++;
                    }
                }
            }
        }
    }
    return mu_index_;
}

void Numerical_Basis::numerical_atomic_wfc(
    const int &ik,
    const int &np,
    ComplexMatrix &psi)
{
    TITLE("Numerical_Basis", "numerical_atomic_wfc");

    std::vector<Vector3<double>> gk(np);
    for (int ig=0; ig<np; ig++)
        gk[ig] = wf.get_1qvec_cartesian(ik, ig);

    const int total_lm = ( ucell.lmax + 1) * ( ucell.lmax + 1);
    matrix ylm(total_lm, np);
    YlmReal::Ylm_Real(total_lm, np, gk.data(), ylm);

    std::vector<double> flq(np);
    for (int it = 0; it < ucell.ntype; it++)
    {
        //OUT("it",it);
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            //OUT("ia",ia);
            complex<double> *sk = wf.get_sk(ik, it, ia);
            for (int l = 0; l < ucell.atoms[it].nwl+1; l++)
            {
                //OUT("l",l);
                complex<double> lphase = pow(IMAG_UNIT, l);
                for (int ic=0; ic < ucell.atoms[it].l_nchi[l]; ic++)
                {
                    //OUT("ic",ic);
                    for (int ig=0; ig<np; ig++)
                    {
                        flq[ig] = this->bessel_basis.
                                  Polynomial_Interpolation(it, l, ic, gk[ig].norm() * ucell.tpiba );
                    }

                    for (int m=0; m<2*l+1; m++)
                    {
                        //OUT("m",m);
                        const int lm = l*l+m;
                        for (int ig=0; ig<np; ig++)
                        {
                            psi( this->mu_index[it](ia,l,ic,m), ig) =
                                lphase * sk[ig] * ylm(lm, ig) * flq[ig];
                        }
                    }
                }
            }
            delete[] sk;    sk=nullptr;
        }
    }
}

void Numerical_Basis::output_info(
    ofstream &ofs,
    const Bessel_Basis &bessel_basis)
{
    // only print out to the information by the first processor
    if (MY_RANK==0)
    {
        ofs.precision(10);
        ofs << ucell.lat0 << endl;

        ofs << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
        ofs << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
        ofs << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;

        ofs << ucell.ntype << " ntype" << endl;
        for (int it=0; it<ucell.ntype; it++)
        {
            ofs << ucell.atoms[it].label << " label" << endl; // mohan add 2009-07-23
            ofs << ucell.atoms[it].na << " na" << endl;
            for (int ia=0; ia<ucell.atoms[it].na; ia++)
            {
                ofs << ucell.atoms[it].tau[ia].x << " " 
                    << ucell.atoms[it].tau[ia].y << " " 
                    << ucell.atoms[it].tau[ia].z << endl;
            }
        }
        // ecutwfc_jlq determine the jlq corresponding to plane wave calculation.
        ofs << pw.ecutwfc << " ecutwfc" << endl; // mohan add 2009-09-08

        // this parameter determine the total number of jlq.
        ofs << bessel_basis.get_ecut() << " ecutwfc_jlq" << endl;//mohan modify 2009-09-08
        ofs << bessel_basis.get_rcut() << " rcut_Jlq" << endl;

        // mohan add 'smooth' and 'sigma' 2009-08-28
        ofs << bessel_basis.get_smooth() << " smooth" << endl;
        ofs << bessel_basis.get_sigma() << " sigma" << endl;

        ofs << bessel_basis.get_tolerence() << " tolerence" << endl;

        ofs << ucell.lmax << " lmax" << endl;
    }

    ofs << scientific;

    ofs << setprecision(8);
    // NOTICE: ofs_warning << "\n The precison may affect the optimize result.";
    
    if (MY_RANK==0)
    {
        ofs << kv.nkstot << " nks" << endl;
        ofs << NBANDS << " nbands" << endl;
        ofs << NLOCAL << " nwfc" << endl;
        ofs << bessel_basis.get_ecut_number() << " ne " << endl;        
    }
}

void Numerical_Basis::output_k(
    ofstream &ofs)
{
    // (1)
    if (MY_RANK==0)
    {
        ofs << "<WEIGHT_OF_KPOINTS>";
    }

    // (2)
    for (int ik=0; ik<kv.nkstot; ik++)
    {
        double kx, ky, kz, wknow;
#ifdef __MPI
        const int pool = Pkpoints.whichpool[ik];
        const int iknow = ik - Pkpoints.startk_pool[MY_POOL];
        if (RANK_IN_POOL==0)
        {
            if (MY_POOL==0)
            {
                if (pool==0)
                {
                    kx = kv.kvec_c[ik].x;
                    ky = kv.kvec_c[ik].y;
                    kz = kv.kvec_c[ik].z;
                    wknow = kv.wk[ik];
                }
                else
                {
                    MPI_Status ierror;
                    MPI_Recv(&kx, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&ky, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4+1, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&kz, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4+2, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&wknow, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4+3, MPI_COMM_WORLD,&ierror);
                }
            }
            else
            {
                if (MY_POOL == pool)
                {
                    MPI_Send(&kv.kvec_c[iknow].x, 1, MPI_DOUBLE, 0, ik*4, MPI_COMM_WORLD);
                    MPI_Send(&kv.kvec_c[iknow].y, 1, MPI_DOUBLE, 0, ik*4+1, MPI_COMM_WORLD);
                    MPI_Send(&kv.kvec_c[iknow].z, 1, MPI_DOUBLE, 0, ik*4+2, MPI_COMM_WORLD);
                    MPI_Send(&kv.wk[iknow], 1, MPI_DOUBLE, 0, ik*4+3, MPI_COMM_WORLD);
                }
            }
        }
        // this barrier is very important
        MPI_Barrier(MPI_COMM_WORLD);
#else
        if (MY_RANK==0)
        {
            kx = kv.kvec_c[ik].x;
            ky = kv.kvec_c[ik].y;
            kz = kv.kvec_c[ik].z;
            wknow = kv.wk[ik];
        }
#endif

        if (MY_RANK==0)
        {
            ofs << "\n" << kx << " " << ky << " " << kz;
            ofs << " " << wknow * 0.5;
        }
    }
    
    if (MY_RANK==0)
    {
        ofs << "\n</WEIGHT_OF_KPOINTS>" << endl;
    }
}

void Numerical_Basis::output_overlap_Q(
    ofstream &ofs,
    const std::vector<ComplexArray> &overlap_Q)
{
    // (3)
    if (MY_RANK==0)
    {
        ofs << "\n<OVERLAP_Q>";
    }

    // (4)
    /*
    if(MY_RANK==0)
    {
    //    	for( int i=0; i<overlap_Q1.getSize(); i++)
    //    	{
    //    		if( i%2==0 ) ofs << "\n";
    //    		ofs << " " << overlap_Q1.ptr[i] << " " << overlap_Q2.ptr[i];
    //    	}
    }
    */

    // Copy to overlap_Q_k for Pkpoints.pool_collection temporaly.
    // It's better to refactor to Pkpoints.pool_collection(overlap_Q) in the future.
    // Peize Lin comments 2021.07.25
    assert(kv.nks>0);
    ComplexArray overlap_Q_k(kv.nks, overlap_Q[0].getBound1(), overlap_Q[0].getBound2(), overlap_Q[0].getBound3());
    for(int ik=0; ik<kv.nks; ++ik)
    {
        std::memcpy(
            overlap_Q_k.ptr + ik*overlap_Q[ik].getSize(),
            overlap_Q[ik].ptr,
            overlap_Q[ik].getSize() * sizeof(std::complex<double>));
    }

    int count = 0;
    for (int ik=0; ik<kv.nkstot; ik++)
    {
        ComplexArray Qtmp(overlap_Q[ik].getBound1(), overlap_Q[ik].getBound2(), overlap_Q[ik].getBound3());
        Qtmp.zero_out();
        Pkpoints.pool_collection(Qtmp.ptr, overlap_Q_k, ik);
        if (MY_RANK==0)
        {
    //        ofs << "\n ik=" << ik;
            // begin data writing.
            const int dim = Qtmp.getSize();
            for (int i=0; i<dim; i++)
            {
                if ( count%4==0 ) ofs << "\n";
                ofs << " " << Qtmp.ptr[i].real() << " " << Qtmp.ptr[i].imag();
                ++count;
            }
            // end data writing.
        }
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    // (5)
    if (MY_RANK==0)
    {
        ofs << "\n</OVERLAP_Q>" << endl;
    }
}

void Numerical_Basis::output_overlap_Sq(
    const string &name,
    ofstream &ofs,
    const std::vector<ComplexArray> &overlap_Sq)
{
    if (MY_RANK==0)
    {
        ofs << "\n<OVERLAP_Sq>";
        ofs.close();
    }

    int count = 0;
    for (int ik=0; ik< kv.nkstot; ik++)
    {
        if ( MY_POOL == Pkpoints.whichpool[ik] )
        {
            if ( RANK_IN_POOL == 0)
            {
                ofs.open(name.c_str(), ios::app);
                const int ik_now = ik - Pkpoints.startk_pool[MY_POOL];

                const int size = overlap_Sq[ik_now].getSize();
                for (int i=0; i<size; i++)
                {
                    if (count%2==0) ofs << "\n";
                    ofs << " " << overlap_Sq[ik_now].ptr[i].real() << " " << overlap_Sq[ik_now].ptr[i].imag();
                    ++count;
                }

                ofs.close();
            }
#ifdef __MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif
        }
        else
        {
#ifdef __MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif
        }

        /*
        if(MY_RANK==0)
        for(int i=0; i< Sq_real[ik].getSize(); i++)
        {
        	if(i%2==0) ofs << "\n";
        	ofs << " " << Sq_real[ik].ptr[i] << " " << Sq_imag[ik].ptr[i];
        }
        */

    }
    if (MY_RANK==0)
    {
        ofs.open(name.c_str(), ios::app);
        ofs << "\n</OVERLAP_Sq>" << endl;
    }
}

// Peize Lin add 2020.04.23
void Numerical_Basis::output_overlap_V(
    ofstream &ofs,
	const matrix &overlap_V)
{
	if (MY_RANK==0)
    {
        ofs << "\n<OVERLAP_V>" <<endl;;
		ofs << overlap_V;
		ofs << "</OVERLAP_V>" <<endl;
	}
}
