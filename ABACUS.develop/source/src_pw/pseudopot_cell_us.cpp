#include "global.h"
#include "pseudopot_cell_us.h"
#include "tools.h"

pseudopot_cell_us::pseudopot_cell_us(void)
{
}

pseudopot_cell_us::~pseudopot_cell_us(void)
{
}

/***************************************************************
//void start_clock (char *);
//void stop_clock (char *);
//void reduce ( int , ::real *qq );
void divide (int nqxq, int startq, int lastq);	// ?
void setqf(::real , ::real , ::real r, int nqf,int , int ilast);	// ?
// void invmat(int llx, matrix ylm, matrix mly, ::real dum); ?

template <class T>
void copyarr(T *orig, T *dest, int len);

pseudopot_cell_us::pseudopot_cell_us(void)
{

	nlx = (lmaxx+1) * (lmaxx+1);
						// maximum number of combined angular momentum
	mx = 2*lqmax-1;		// maximum magnetic angular momentum of Q
	lpx = matrix(nlx, nlx);	// (nlx,nlx); for each pair of combined momenta lm(1),lm(2):
							// maximum combined angular momentum LM
    lpl = realArray (mx,nlx,nlx);	//(nlx,nlx,mx)  list of combined angular momenta  LM

	ap = realArray(nlx,nlx,lqmax*lqmax);	//(lqmax*lqmax,nlx,nlx)
	// Clebsch-Gordan coefficients for spherical harmonics
//	qq = realArray(ntyp,nhm,nhm);		//(:,:,:), the q functions in the solid
//	qrad = realArray(nbeta*(nbeta+1),nbrx,nbrx);		//(:,:,:,:), radial FT of Q functions

//	int nkbus;        // as above, for US-PP only

//	realArray deeq;		//(:,:,:,:), the integral of V_eff and Q_{nm}
//	realArray becsum;	//(:,:,:), \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
//	realArray beta;		//(:,:,:), beta functions for CP (without struct.factor) //not sure!

	// from spin_orb
	lspinorb = true;
	domag = true;			// if .TRUE. this is a spin-robit calculation
	// from us
	okvan = false;;			// if .TRUE. at least one pseudo is Vanderbilt
	rot_ylm = ComplexMatrix(2*lmaxx+1,2*lmaxx+1);
				//(2*lmaxx+1,2*lmaxx+1), transform real
				// spherical harmonics into complex ones
//	fcoef = ComplexArray(ntyp,nhm,nhm,2,2);		// (:,:,:,:,:), function needed to
							// account for spinors.

	// from uspp
//	dvan_so = ComplexArray(ntyp,nhm,nhm,4);	// (:,:,:,:), D_{nm}
//	qq_so = ComplexArray(ntyp,nhm,nhm,4);	// (:,:,:,:), Q_{nm}
}

pseudopot_cell_us::~pseudopot_cell_us(void)
{
}

//----------------------------------------------------------------------
// from init_us_1.f90
void pseudopot_cell_us::init_us_1(PW_Basis pw)
{
	//------------------------------------------------------------------
	//   This routine performs the following tasks:
	//   a) For each non vanderbilt pseudopotential it computes the D and
	//      the betar in the same form of the Vanderbilt pseudopotential.
	//   b) It computes the indices indv which establish the correspondence
	//      nh <-> beta in the atom
	//   c) It computes the indices nhtol which establish the correspondence
	//      nh <-> angular momentum of the beta function
	//   d) It computes the indices nhtolm which establish the correspondence
	//      nh <-> combined (l,m) index for the beta function.
	//   e) It computes the coefficients c_{LM}^{nm} which relates the
	//      spherical harmonics in the Q expansion
	//   f) It computes the radial fourier transform of the Q function on
	//      all the g vectors
	//   g) It computes the q terms which define the S matrix.
	//   h) It fills the interpolation table for the beta functions

	// USE parameters, ONLY : lmaxx, nbrx, lqmax
	// USE constants,  ONLY : fpi
	// USE atom,       ONLY : r, rab
	// USE ions_base,  ONLY : ntyp => nsp
	// USE cell_base,  ONLY : omega, tpiba
	// USE gvect,      ONLY : g, gg
	// USE pseud,      ONLY : lloc, lmax
	// USE lsda_mod,   ONLY : nspin
	// USE us,         ONLY : okvan, nqxq, dq, nqx, tab, qrad
	// USE uspp,       ONLY : nhtol, nhtoj, nhtolm, dvan, qq, indv, ap, aainit,
	//						  qq_so, dvan_so
	// USE uspp_param, ONLY : lmaxq, dion, betar, qfunc, qfcoef, rinner, nbeta,
	//						  kkbeta, nqf, nqlc, lll, jjj, lmaxkb, nh, tvanp, nhm
	// USE spin_orb,   ONLY : lspinorb, rot_ylm, fcoef

	//     here a few local variables
	int nt, ih, jh, nb, mb, nmb, l, m, ir, iq, is, startq,
		lastq, ilast, ndm;
	// various counters

	::real *aux , *aux1 , *besr ;
	realArray qtot ;		// (:,:,:)
	// various work space

	::real prefr, pref, q, qi;
	// the prefactor of the q functions
	// the prefactor of the beta functions
	// the modulus of g for each shell
	// q-point grid for interpolation

	matrix ylmk0;
	// the spherical harmonics

	::real *vll;		// (0:lmaxx),
	::real vqint, sqrt2, j;
	// the denominator in KB case
	// interpolated value

	int n1, m0, m1, n, li, mi, vi, vj, ijs, is1, is2,
		lk, mk, vk, kh, lh;
	realArray sph_ind;
	complex <::real> coeff, qgm;		//(1);?
	realArray spinor;
	::real ji, jk;

	::real x;
	complex <::real> cx;
	complex <::real> c0,c10;
	c0 = complex <::real>(0,0);
	c10 = complex <::real>(1,0);

	int ntyp,nbeta,kkbeta;
	int nh;
	::real omega;
	int *lll;
	omega = pw.omega;
	ntyp = pw.ntype;
//	start_clock ("init_us_1");
	//    Initialization of the variables

	// ndm = MAXVAL (kkbeta),(1:ntyp));
	ndm = 0;
	for(n=0;n<ntyp;n++){
		kkbeta = pw.atoms[n].kkbeta;
		ndm = (ndm > kkbeta)?ndm:kkbeta;
	}
	// allocate
	aux = new real[ndm];
	aux1 = new real[ndm];
	besr = new real [ndm];
	qtot = realArray(nbrx, nbrx, ndm);
	ylmk0 = matrix( lmaxq , lmaxq);
	ap.zero_out();	// (:,:,:)   = 0.0;
	if (lmaxq > 0)
		qrad.zero_out();	//(:,:,:,:)= 0.0;
	prefr = fpi / omega;
	if (lspinorb){
		//  In the spin-orbit case we need the unitary matrix u which
		// rotates the real spherical harmonics and yields the complex ones.

		sqrt2=1.0/sqrt(2.0);	// dsqrt()
		rot_ylm.zero_out();		//(0.0,0.0)
		l=lmaxx;
		rot_ylm(l+1,1)=c10;		//(1.0,0.0);
		for(n1=2;n1<2*l+1;n1+=2){ // do n1=2,2*l+1,2
			m=n1/2;
			n=l+1-m;
			rot_ylm(n,n1)=std::complex <real>(pow(-1.0, m)*sqrt2,0.0);
			rot_ylm(n,n1+1)=std::complex <real>(0.0,pow(-(-1.0),m)*sqrt2);
			n=l+1+m;
			rot_ylm(n,n1)=std::complex <real>(sqrt2,0.0);
			rot_ylm(n,n1+1)=std::complex <real>(0.0, sqrt2);
		} //  enddo
		fcoef.zero_out();
		dvan_so.zero_out();
		qq_so.zero_out();
	}else{
		qq.zero_out();	// (:,:,:)   = 0.d0
		dvan.zero_out();
	} // endif
	//   For each pseudopotential we initialize the indices nhtol, nhtolm,
	//   nhtoj, indv, and if the pseudopotential is of KB type we initialize
	// the atomic D terms
	for(nt=0;nt<ntyp;nt++){ // do nt = 1, ntyp
		ih = 1;
		nbeta = pw.atoms[nt].nbeta;
		nh = pw.atoms[nt].nh;
		for(nb=0;nb<nbeta;nb++){ // do nb = 1, nbeta (nt)
			l = pw.atoms[nt].lll [nb];
// ?			j = pw.atoms[nt].jjj [nb];
			for(m=0;m<2*l+1;m++){ // do m = 1, 2 * l + 1
				nhtol (ih, nt) = l;
				nhtolm(ih, nt) = l*l+m;
				nhtoj (ih, nt) = j;
				indv  (ih, nt) = nb;
				ih = ih + 1;
			} // enddo
		} // enddo

		//    From now on the only difference between KB and US pseudopotentials
		//    is in the presence of the q and Q functions.
		//    Here we initialize the D of the solid
		if (lspinorb){
			//  first calculate the fcoef coefficients
			for(ih=0;ih<nh;ih++){ // do ih = 1, nh (nt)
				li = nhtol(ih, nt);
				ji = nhtoj(ih, nt);
				mi = nhtolm(ih, nt)-li*li;
				vi = indv (ih, nt);
				for(kh=0;kh<nh;kh++){ // do kh=1,nh(nt)
					lk = nhtol(kh, nt);
					jk = nhtoj(kh, nt);
					mk = nhtolm(kh, nt)-lk*lk;
					vk = indv (kh, nt);
					if (li == lk && fabs(ji-jk) < 1.e-7){
						for(is1=0;is1<2;is1++){ // do is1=1,2
							for(is2=0;is2<2;is2++){ // do is2=1,2
								coeff = c0;		//(0.d0, 0.d0)
								for(m=li-1;m<=li;m++){ // do m=-li-1, li ?
									m0= sph_ind.valueof(li,ji,m,is1) + lmaxx + 1;
									m1= sph_ind.valueof(lk,jk,m,is2) + lmaxx + 1;
									coeff=coeff + rot_ylm(m0,mi)*spinor.valueof(li,ji,m,is1)*
										conjg(rot_ylm(m1,mk))*spinor.valueof(lk,jk,m,is2);
								} //  enddo
								fcoef.assign(nt,ih,kh,is1,is2, coeff);
							} // enddo
						} // enddo
					} // endif
				} // enddo
			} // enddo
			//   and calculate the bare coefficients

			for(ih=0;ih<nh;ih++){ // do ih = 1, nh (nt)
				vi = indv (ih, nt);
				for(jh=0;jh<nh;jh++){ // do jh = 1, nh (nt)
					vj = indv (jh, nt);
					ijs=0;
					for(is1=0;is1<2;is1++){ // do is1=1,2
						for(is2=0;is2<2;is2++){ //do is2=1,2
							ijs=ijs+1;
							dvan_so.assign(nt,ih,jh,ijs,
								fcoef.valueof(nt,ih,jh,is1,is2));
							if (vi != vj)
								fcoef.assign(nt,ih,jh,is1,is2, c0);	//=(0.d0,0.d0)
						} // enddo
					} // enddo
				} // enddo
			} // enddo
		} else{
			for(ih=0;ih<nh;ih++){ // do ih = 1, nh (nt)
				for(jh=0;jh<nh;jh++){ // do jh = 1, nh (nt)
					if (nhtol (ih, nt) == nhtol (jh, nt) &&
						nhtolm(ih, nt) == nhtolm(jh, nt) ){
						ir = indv (nt,ih);
						is = indv (nt,jh);
						dvan.assign (nt,ih, jh,
							pw.atoms[nt].dion (ir, is));
					} // endif
				} // enddo
			} // enddo
		} // endif
	} // enddo

	//  compute Clebsch-Gordan coefficients
	if (okvan){
		aainit (lmaxkb + 1);
	//   here for the US types we compute the Fourier transform of the
	//   Q functions.

	divide (nqxq, startq, lastq);
	for(nt=0;nt<ntyp;nt++){ // do nt = 1, ntyp
		nbeta = pw.atoms[nt].nbeta;
		kkbeta = pw.atoms[nt].kkbeta;
		copyarr(pw.atoms[nt].lll,lll,nbeta);
		if (pw.atoms[nt].tvanp ){
			for(l=0;l<pw.atoms[nt].ppus.nqlc-1;l++){ // do l = 0, nqlc (nt) - 1
				// first we build for each nb,mb,l the total Q(|r|) function
				// note that l is the true angular momentum, and the arrays
				// have dimensions 1..l+1
				for(nb=0;nb<nbeta;nb++){ // do nb = 1, nbeta (nt)
					for(mb=nb;mb<nbeta;mb++){ // do mb = nb, nbeta (nt)
						if ( (l >= fabs (lll [nb] - lll [mb]) ) &&
							(l <= lll [nb] + lll [mb]) &&
							(((int)(l + lll [nb] + lll [mb]) % 2) == 0) ){
							// mod (l + lll [nb] + lll [mb], 2)
							for(ir=0;ir<kkbeta;ir++){ // do ir = 1, kkbeta (nt)
								if (pw.atoms[nt].r [ir] >=
									pw.atoms[nt].ppus.rinner [l + 1]){
									qtot.assign (nb, mb, ir,
										pw.atoms[nt].ppus.qfunc.valueof (ir, nb, mb));
								}else{
									ilast = ir;
								} // endif
							} // enddo
							if (pw.atoms[nt].ppus.rinner [l + 1] > 0.0)
								setqf(pw.atoms[nt].ppus.qfcoef.valueof (1, l+1, nb, mb),
								qtot.valueof(nb,mb,1), pw.atoms[nt].r[1],
								pw.atoms[nt].ppus.nqf,l,ilast);
						} // endif
					} // enddo
				} // enddo
				// here we compute the spherical bessel function for each |g|
				for(iq=startq;iq<lastq;iq++){ // do iq = startq, lastq
					q = (iq - 1) * dq * pw.tpiba;
					sph_bes (kkbeta, pw.atoms[nt].r ,
						q, l, aux);
					//   and then we integrate with all the Q functions
					for(nb=0;nb<nbeta;nb++){ // do nb = 1, nbeta (nt)
						//    the Q are symmetric with respect to indices
						for(mb=nb;mb<nbeta;mb++){ // do mb = nb, nbeta (nt)
							nmb = mb * (mb - 1) / 2 + nb;
							if ( (l >= fabs (lll [nb] - lll [mb]) ) &&
								(l <= lll [nb] + lll [mb]) &&
								(((int)(l + lll[nb] + lll[mb]) % 2) == 0) ){
								// mod (l + lll[nb] + lll[mb], 2)
								for(ir=0;ir<kkbeta;ir++){ // do ir = 1, kkbeta (nt)
									aux1 [ir] = aux [ir] * qtot.valueof (ir, nb, mb);
								} // enddo
								x = qrad.valueof(nt,iq,nmb,l + 1);
								simpson (kkbeta, aux1,
									pw.atoms[nt].rab, x);
								qrad.assign(nt,iq,nmb,l + 1,x);
							} // endif
						} // enddo
					} // enddo
					// igl
				} // enddo
			} // l
		} // enddo
//        qrad.assign (, , , nt, qrad.valueof (, :, :, nt)*prefr);
#ifdef __PARA
//        reduce (nqxq * nbrx * (nbrx + 1) / 2 * lmaxq, qrad.valueof (1, 1, 1, nt) );
#endif
	} // endif
	// ntyp
	} // enddo

	//   and finally we compute the qq coefficients by integrating the Q.
	//   q are the g=0 components of Q.
#ifdef __PARA
//?	if (pw.gg [1] > 1.0e-8)
//?		goto 100;
#endif
	ylmr2 (lmaxq * lmaxq, 1, pw.g, pw.gg, ylmk0);
	for(nt=0;nt<ntyp;nt++){ // do nt = 1, ntyp
		nh = pw.atoms[nt].nh;
		if (pw.atoms[nt].tvanp) {
			if (lspinorb){
				for(ih=0;ih<nh;ih++){ // do ih=1,nh(nt)
					for(jh=0;jh<nh;jh++){ // do jh=1,nh(nt)
						qvan2 (1, ih, jh, nt, pw.gg, qgm, ylmk0);
						for(kh=0;kh<nh;kh++){ // do kh=1,nh(nt)
							for(lh=0;lh<nh;lh++){ // do lh=1,nh(nt)
								ijs=0;
								for(is1 =0;is1<2;is1++){ // do is1=1,2
									for(is2=0;is2<2;is2++){ // do is2=1,2
										ijs=ijs+1;
										for(is=0;is<2;is++){ // do is=1,2
											cx = qq_so.valueof(nt,kh,lh,ijs) +
												omega * qgm.real() *
												fcoef.valueof(nt,kh,ih,is1,is) *
												fcoef.valueof(nt,jh,lh,is,is2);
											qq_so.assign(nt,kh,lh,ijs,cx);
										} // enddo
									} // enddo
								} // enddo
							} // enddo
						} // enddo
					} // enddo
				} // enddo
			}else{
				for(ih=0;ih<nh;ih++){ // do ih = 1, nh (nt)
					for(jh=ih;jh<nh;jh++){ // do jh = ih, nh (nt)
						qvan2 (1, ih, jh, nt, pw.gg, qgm, ylmk0);
						x = omega * qgm.real() ;
						qq.assign (nt,ih, jh, x);
						x = qq.valueof (nt,ih, jh);
						qq.assign (nt,jh, ih, x);
					} // enddo
				} // enddo
			} // endif
		} // endif
	} // enddo
#ifdef __PARA
//	100 continue;
	if (lspinorb){
		reduce ( nhm * nhm * ntyp * 8, qq_so );
	}else{
		reduce ( nhm * nhm * ntyp, qq );
	} // endif
#endif
//	stop_clock ("init_us_1");
**********************************************************/
/**********************************************************
// He Lixin: this block is used for non-local potential
//     fill the interpolation table tab
************************************************************/

/*
	pref = fpi / sqrt (omega);
	divide (nqx, startq, lastq);
	tab.zero_out();
	for(nt=0;nt<ntyp;nt++){ // do nt = 1, ntyp
		nbeta = pw.atoms[nt].nbeta;
		kkbeta = pw.atoms[nt].kkbeta;
		for(nb=0;nb<nbeta;nb++){ // do nb = 1, nbeta (nt)
			l = pw.atoms[nt].lll [nb];
			for(iq=startq;iq<lastq;iq++){ // do iq = startq, lastq
				qi = (iq - 1) * dq;
				sph_bes (kkbeta, pw.atoms[nt].r , qi, l, besr);
				for(ir=0;ir<kkbeta;ir++){ // do ir = 1, kkbeta (nt)
					aux [ir] = pw.atoms[nt].betar (ir, nb) *
						besr [ir] * pw.atoms[nt].r [ir];
				} // enddo
				simpson (kkbeta , aux, pw.atoms[nt].rab ,
					vqint);
				tab.assign (iq, nb, nt, vqint * pref);
			} // enddo
		} // enddo
	} // enddo

#ifdef __PARA
	reduce (nqx * nbrx * ntyp, tab);
#endif
	// deallocate (ylmk0)
	// deallocate (qtot)
	// deallocate (besr)
	// deallocate (aux1)
	// deallocate (aux)
	return;

} // end subroutine init_us_1
*/
/*
// -----------------------------------------------------------------------
// from uspp.f90
void pseudopot_cell_us::aainit(int lli)
{
	//-----------------------------------------------------------------------
    //
    // this routine computes the coefficients of the expansion of the product
    // of two real spherical harmonics into real spherical harmonics.
    //
    //     Y_limi(r) * Y_ljmj(r) = \sum_LM  ap(LM,limi,ljmj)  Y_LM(r)
    //
    // On output:
    // ap     the expansion coefficients
    // lpx    for each input limi,ljmj is the number of LM in the sum
    // lpl    for each input limi,ljmj points to the allowed LM
    //
    // The indices limi,ljmj and LM assume the order for real spherical
    // harmonics given in routine ylmr2
    //
    // input: the maximum li considered
    // integer :: lli

    // local variables
	int llx, l, li, lj;
    Vector3 <::real> *r;	//(:,:),
	::real *rr;	//(:),
	matrix ylm;	//(:,:),
	matrix mly;	//(:,:)
    // an array of random vectors: r(3,llx)
    // the norm of r: rr(llx)
    // the real spherical harmonics for array r: ylm(llx,llx)
    // the inverse of ylm considered as a matrix: mly(llx,llx)
    ::real dum;
	int lli2;
	lli2 = lli * lli;

    if (lli < 0)
		cout << "\n error,aainit,lli not allowed,"
		<< lli;

    if (lli2 > nlx)
		cout << "\n error, aainit, nlx is too small,"
		<< lli2;

    llx = (2*lli-1);
	llx *= llx;
    if (2*lli-1 > lqmax)
		cout << "\n error, aainit, ap leading dimension is too small,"
		<< llx;

    rr = new ::real[ llx ];
    r = new Vector3 <::real>[llx];
    ylm = matrix( llx, llx );
    mly = matrix( llx, llx );

//    r.zero_out();
    ylm.zero_out();
    mly.zero_out();
    ap.zero_out();	//(:,:,:)= 0.d0

    // - generate an array of random vectors (uniform deviate on unitary sphere)
    gen_rndm_r(llx,r,rr);

    // - generate the real spherical harmonics for the array: ylm(ir,lm)
    ylmr2(llx,llx,r,rr,ylm);

    //-  store the inverse of ylm(ir,lm) in mly(lm,ir)
    // invmat(llx, ylm, mly, dum);
	mly = ylm.Inverse();	//? dum

	::real x;
    //-  for each li,lj compute ap(l,li,lj) and the indices, lpx and lpl
	for(li=0;li<lli2;li++){ // do li = 1, lli*lli
		for(lj=0;lj<lli2;lj++){ // do lj = 1, lli*lli
			lpx(li,lj)=0;
			for(l=0;l<llx;l++){ // do l = 1, llx
				x = compute_ap(l,li,lj,llx,ylm,mly);
				ap.assign(li,lj,l,x);
				if (fabs(ap.valueof(li,lj,l)) > 1.e-3){
					lpx(li,lj) = lpx(li,lj) + 1;
					if (lpx(li,lj) > mx)
						cout << "\n error, aainit, mx dimension too small,"
						<< lpx(li,lj);
					lpl.assign(li,lj,lpx(li,lj),l);
				} // end if
			} // end do
		} // end do
	} // end do

    //deallocate(mly)
    //deallocate(ylm)
    //deallocate(rr)
    //deallocate(r)
    return;
} // end subroutine aainit

// -----------------------------------------------------------------------
// from uspp.f90
void gen_rndm_r(int llx,Vector3 <::real> *r,::real *rr)	// Vector3 *r ?
{
    //-----------------------------------------------------------------------
    // - generate an array of random vectors (uniform deviate on unitary sphere)
    //
    // USE constants, ONLY: tpi

    // first the I/O variables
    // int llx;				// input: the dimension of r and rr
    // matrix r(3,llx),		// output: an array of random vectors
    // real *rr(llx)		// output: the norm of r

    // here the local variables
    int ir;
    ::real costheta, sintheta, phi;

    for(ir=0;ir<llx;ir++){ // do ir = 1, llx
		costheta = 2.0 * rand() - 1.0;	// rndm()
		sintheta = sqrt ( 1.0 - costheta*costheta);
		phi = tpi * rand();	//rndm();
		r[ir].x = sintheta * cos(phi);
		r[ir].y = sintheta * sin(phi);
		r[ir].z = costheta;
		rr[ir]   = 1.0;
	} // end do
//	printrm(ofs, "r :",r);
    return;
} // end subroutine gen_rndm_r

// -----------------------------------------------------------------------
// from uspp.f90
::real compute_ap(int l,int li,int lj,int llx,matrix ylm,matrix mly)
{
    //-----------------------------------------------------------------------
    //-  given an l and a li,lj pair compute ap(l,li,lj)
    //
    // first the I/O variables
    //
	// int llx,         // the dimension of ylm and mly
    // int l,li,lj       // the arguments of the array ap

    // real(kind=DP) ::
    // function compute_ap,  // this function
    //   ylm(llx,llx), // the real spherical harmonics for array r
    //   mly(llx,llx)  // the inverse of ylm considered as a matrix

    // here the local variables
    int ir;
	::real apx;

    apx = 0;
	for(ir=0;ir<llx;ir++){ // do ir = 1,llx
		apx = apx + mly(l,ir) * ylm(ir,li) * ylm(ir,lj);
	} // end do

    return apx;
} // end function compute_ap


void printus(ofstream &ofs)
{
	// which data is need printed ?
}

template <class T>
void copyarr(T *orig, T *dest, int len)
{
	for(int i=0;i<len;i++){
		dest[i] = orig[i];
	}
	return;
}
*/
