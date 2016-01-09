#include "../src_pw/global.h"
#include "myfunc.h"
#include "mymath.h"


// from gradcorr.f90
void gradcorr(double &etxc, double &vtxc, matrix &v)
{
	BLOCK_HERE("gradcorr");
    if (xcf.igcx == 0  &&  xcf.igcc == 0)
    {
        return;
    }

    double zeta, rh, grh2;
    int k, ipol, is;

    double *dh;
    double grho2[2], sx, sc, v1x, v2x, v1c, v2c, v1xup, v1xdw,
    v2xup, v2xdw, v1cup, v1cdw , etxcgc, vtxcgc, segno, arho, fac;

    etxcgc = 0.0;
    vtxcgc = 0.0;

    realArray h(NSPIN, pw.nrxx, 3);
    realArray grho(NSPIN, pw.nrxx, 3);

    // calculate the gradient of rho+rho_core in real space
    Vector3<double> *grho_v3 = new Vector3<double>[pw.nrxx];
    Vector3<double> *h_v3 = new Vector3<double>[pw.nrxx];

    fac = 1.0 / NSPIN;

    for (is=0; is<NSPIN; is++)
    {
        for (k = 0; k <pw.nrxx ;k++)
        {
            grho_v3[k].x = grho(is, k, 0);
            grho_v3[k].y = grho(is, k, 1);
            grho_v3[k].z = grho(is, k, 2);
        }
        daxpy(pw.nrxx , fac, chr.rho_core, 1, chr.rho1, 1);
        gradient(chr.rho1, grho_v3);
    }

    for (k=0; k<pw.nrxx; k++)
    {
        for (is = 0;is < NSPIN;is++)
        {
            grho2 [is] = grho(is, k, 0) * grho(is, k, 0) +
                         grho(is, k, 1) * grho(is, k, 1) +
                         grho(is, k, 2) * grho(is, k, 2);
        }

        if (NSPIN == 1)
        {
            // This is the spin-unpolarised case
            arho = abs( chr.rho1[k] );
            segno = arho;		// sign (1.0, rho (k, 1) );

            if (arho > epsr && grho2 [1]  > epsg)
            {
                gcxc(arho, grho2[1], sx, sc, v1x, v2x, v1c, v2c);  //?grho2[0]
                // first term of the gradient correction : D(rho*Exc)/D(rho)
                v(0, k) = v(0, k) + e2 * (v1x + v1c);
                // h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|

                for (ipol = 0;ipol < 3;ipol++)  //  do ipol = 1, 3
                {
                    h(0, k, ipol) =  e2 * (v2x + v2c) * grho(0, k, ipol);
                } // enddo

                vtxcgc = vtxcgc + e2 * (v1x + v1c) * ( chr.rho1[k] - chr.rho_core[k]);

                etxcgc = etxcgc + e2 * (sx + sc) * segno;
            }
            else
            {
                for (ipol = 0;ipol < 3;ipol++)  // do ipol = 1, 3
                {
                    h(0, k, ipol) = 0.0;
                } //  enddo
            } // endif
        }
        else
        {
            // spin-polarised case
            gcx_spin( chr.rho1[k], chr.rho2[k], grho2 [0], grho2 [1],
                     sx, v1xup, v1xdw, v2xup, v2xdw);
            rh = chr.rho1[k] + chr.rho2[k];

            if (rh > epsr)
            {
                zeta = ( chr.rho1[k] - chr.rho2[k] ) / rh;
                Vector3 < double> v3;
                v3.x = grho(0, k, 0) + grho(1, k, 0);
                v3.y = grho(0, k, 1) + grho(1, k, 1);
                v3.z = grho(0, k, 2) + grho(1, k, 2);
                grh2 = v3 * v3;
                //		grh2 = (grho (k, 1, 1) + grho (k, 1, 2) ) **2 +
                //			   (grho (k, 2, 1) + grho (k, 2, 2) ) **2 +
                //			   (grho (k, 3, 1) + grho (k, 3, 2) ) **2;
                gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
            }
            else
            {
                sc = 0.0;
                v1cup = 0.0;
                v1cdw = 0.0;
                v2c = 0.0;
            } //  endif

            // first term of the gradient correction : D(rho*Exc)/D(rho)
            v(0, k) = v(0, k) + e2 * (v1xup + v1cup);

            v(1, k) = v(1, k) + e2 * (v1xdw + v1cdw);

            // h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
            for (ipol = 0;ipol < 3;ipol++)  // do ipol = 1, 3
            {
                h(0, k, ipol) =  e2 * ((v2xup + v2c) * grho(0, k, ipol)
                                       + v2c * grho(1, k, ipol));
                h(1, k, ipol) =  e2 * ((v2xdw + v2c) * grho(1, k, ipol)
                                       + v2c * grho(0, k, ipol));
            } //  enddo

            vtxcgc = vtxcgc + e2 * (v1xup + v1cup) * ( chr.rho1[k] -
                     chr.rho_core[k] * fac);

            vtxcgc = vtxcgc + e2 * (v1xdw + v1cdw) * ( chr.rho2[k] -
                     chr.rho_core[k] * fac);

            etxcgc = etxcgc + e2 * (sx + sc);
        } //  endif
    } //  enddo


	daxpy(pw.nrxx, -fac, chr.rho_core, 1, chr.rho1, 1);
	if(NSPIN==2)
	{
		daxpy(pw.nrxx, -fac, chr.rho_core, 1, chr.rho2, 1);
	}
	
    dh = new double [pw.nrxx];

    // second term of the gradient correction :
    // \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )

    for (is = 0;is < NSPIN;is++)  // do is = 1, n_spin
    {
        for (k = 0; k < pw.nrxx; k++)
        {
            h_v3[k].x = h(is, k, 0);
            h_v3[k].y = h(is, k, 1);
            h_v3[k].z = h(is, k, 2);
        }

        grad_dot(pw.ncx, pw.ncy, pw.ncz, pw.ncxyz, h_v3, pw.ngmc, pw.g, pw.ig2fftc, ucell.lat0, dh);

        for (k = 0;k < pw.nrxx;k++)  // do k = 1, ncxyz
        {
            v(is, k) = v(is, k) - dh [k];
        }

		if(is==0)
		{
        	for (k = 0;k < pw.nrxx;k++)
			{
            	vtxcgc = vtxcgc - dh [k] * chr.rho1[k];
			}
		}
		else if(is==1)
		{
        	for (k = 0;k < pw.nrxx;k++)
			{
            	vtxcgc = vtxcgc - dh [k] * chr.rho2[k];
			}
		}
		
    } //  enddo

    double vtxcaux = vtxc;
    double etxcaux = etxc;

#ifdef __MPI
    MPI_Allreduce(&vtxcaux,&vtxc,1,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
    MPI_Allreduce(&etxcaux,&etxc,1,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
#endif

    vtxc += ucell.omega * vtxcgc / pw.ncxyz;//(nr1 * nr2 * nr3);
    etxc += ucell.omega * etxcgc / pw.ncxyz;//(nr1 * nr2 * nr3);

    delete [] grho_v3;
    delete [] h_v3;
    delete [] dh;

    return;
} // end subroutine gradcorr

// from gradcorr.f90
void gradient( double *a, Vector3<double> *ga)
{
    //-------------------------------------------------------
    // Calculates ga = \grad a in R-space (a is also in R-space)
    // use gvect, only: int *nlm;
    // use wvfct, only: bool gamma_only;
    double tpiba = TWO_PI / lat0;

    int n=0;
    int ipol=0;
    matrix aux(2, pw.nrxx);
    matrix gaux(2, pw.nrxx);

    double *aux_1d= new double[pw.nrxx];
    double *gaux_1d= new double[pw.nrxx];
    double *ga_1d= new double[pw.nrxx];

    matrix gx(ngmc, 3);

    for (n = 0;n < ngmc;n++)
    {
        gx(n, 0) = g[n].x;
        gx(n, 1) = g[n].y;
        gx(n, 2) = g[n].z;
    }

    // copy a(r) to complex array...
    // aux(2,:) = 0.0l; ?
    for (n = 0;n < ncxyz;n++)
    {
        aux(0, n) = 0.01;
        aux(1, n) = 0.01;
    }

    for (n = 0;n < ncxyz;n++)
    {
        aux_1d[n] = aux(1, n);
    }

    dcopy(ncxyz, a, 1, aux_1d, 1);
    // bring a(r) to G-space, a(G) ...

    // cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1);
    //	setupFFT( ncx, ncy,ncz);
#ifdef __MPI
#else
//	fftchg.FFT3D(aux, -1);
#endif
    // multiply by (iG) to get (\grad_ipol a)(G) ...

    // ga(:,:) = 0.0;
    //	ga.zero_out();

    for (ipol = 0;ipol < 3;ipol++)
    {
        // gaux(:,:) = 0.0;
        gaux.zero_out();

        for (n = 0;n < ngmc;n++)  // do n = 1, ngm
        {
            gaux(0, ig2fftc[n]) = - gx(n, 0) * aux(1, ig2fftc[n]);
            gaux(1, ig2fftc[n]) =   gx(n, 1) * aux(0, ig2fftc[n]);
        } //  enddo

        if (wf.gamma_only)
        {
//			for(n=0;n<ngm;n++){ //	do n = 1, ngm
//				gaux (0, nlm[n] ) =   gaux (0, nl[n] );
//				gaux (1, nlm[n] ) = - gaux (1, nl[n] );
//			} // enddo
        } // end if

        // bring back to R-space, (\grad_ipol a)(r) ...
        // cft3 (gaux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1);
//		setupFFT3D( ncx, ncy,ncz);
//
#ifdef __MPI
#else
        //	fftchg.FFT3D(gaux, 1);
#endif

        // ...and add the factor 2\pi/a  missing in the definition of G
        // call DAXPY (nrxx, tpiba, gaux, 2, ga (ipol, 1), 3)
        //if(ipol = 0){
        if (ipol == 0) //mohan modify 2008-01-23
        {
            for (n = 0; n < ncxyz; n++)
            {
                ga_1d[n] = ga[n].x;
            }

            //}else if(ipol = 1){
        }
        else if (ipol == 1) //mohan modify 2008-01-23
        {
            for (n = 0; n < ncxyz; n++)
            {
                ga_1d[n] = ga[n].y;
            }
        }
        else
        {
            for (n = 0; n < ncxyz; n++)
            {
                ga_1d[n] = ga[n].z;
            }
        }

        //call DAXPY (nrxx, tpiba, gaux, 2, ga (ipol, 1), 3)
        daxpy(ncxyz, tpiba, gaux_1d, 1, ga_1d, 1);
    } //  enddo

    delete [] aux_1d;
    aux_1d = 0;
    delete [] gaux_1d;
    gaux_1d = 0;
    delete [] ga_1d;
    ga_1d = 0;

    return;
} // end subroutine gradient

// from gradcorr.f90
//---------------------------------------------------------------
//void grad_dot (int nrx1, int nrx2, int nrx3, int ncx, int ncy, int ncz,
//			   int ncxyz, matrix a, int ngm, Vector3 < double> *g, int *nl, double alat,
//			   double *da)
void grad_dot(int ncx, int ncy, int ncz, int ncxyz, Vector3 < double> *a, int ngmc,
              Vector3 < double> *g, int *ig2fftc, double lat0, double *da)
{
    //-----------------------------------------------------------

    // Calculates da = \sum_i \grad_i a_i in R-space

    // use gvect, only: nlm
    // use wvfct, only: gamma_only

    // integer :: nrx1, nrx2, nrx3, nr1, nr2, nr3, ncxyz, ngm, nl (ngm);
    // real(kind=DP) :: a (3, ncxyz), g (3, ngm), da (ncxyz), alat;
    int n, ipol;
    matrix aux;		// (:,:),
    matrix gaux;	// (:,:);
    aux.create(2, ncxyz);	// allocate (aux( 2,ncxyz));
    gaux.create(2, ncxyz);	// allocate (gaux(2,ncxyz));

    double tpiba= TWO_PI / lat0;
    matrix gx;
    gx.create(ngmc, 3);

    for (n = 0;n < ngmc;n++)
    {
        gx(n, 0) = g[n].x;
        gx(n, 1) = g[n].y;
        gx(n, 2) = g[n].z;
    }

    for (ipol = 0;ipol < 3;ipol++)  // do ipol = 1, 3
    {
        // copy a(ipol,r) to a complex array...

        for (n = 0;n < ncxyz;n++)
        {
            aux(0, n) = 0.0;
        }

        //dcopy (ncxyz, a (ipol, 1), 3, aux, 2);
        //if(ipol = 0){
        if (ipol == 0) //mohan modify 2007-01-23
        {
            for (n = 0; n < ncxyz; n++)
            {
                aux(0, n) = a[n].x;
            }

            //}else if(ipol = 1){
        }
        else if (ipol == 1) //mohan modify 2007-01-23
        {
            for (n = 0; n < ncxyz; n++)
            {
                aux(0, n) = a[n].y;
            }
        }
        else
        {
            for (n = 0; n < ncxyz; n++)
            {
                aux(0, n) = a[n].z;
            }
        }

        // bring a(ipol,r) to G-space, a(G) ...
        // cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1);
//		setupFFT3D( ncx, ncy,ncz);
//
#ifdef __MPI
#else
        //	fftchg.FFT3D(aux, -1);
#endif

        // multiply by (iG) to get (\grad_ipol a)(G) ...
        for (n = 0;n < ngmc;n++)  // do n = 1, ngm
        {
            gaux(0, ig2fftc[n]) -= gx(n, ipol) * aux(1, ig2fftc[n]);
            gaux(1, ig2fftc[n]) += gx(n, ipol) * aux(0, ig2fftc[n]);
        } //  enddo
    } //  enddo

    if (wf.gamma_only)
    {
//		for(n=0;n<ngm;n++){ // do n = 1, ngm
//			gaux ( nlm[n],0 ) =   gaux ( ig2fftc[n],0 );
//			gaux ( nlm[n],1 ) = - gaux ( ig2fftc[n],1 );
//		} // enddo
    } // end if

    //  bring back to R-space, (\grad_ipol a)(r) ...
    // cft3 (gaux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1);
//	setupFFT( ncx, ncy,ncz);
//

#ifdef __MPI
#else
//	fftchg.FFT3D(gaux, 1);
#endif

    // ...add the factor 2\pi/a  missing in the definition of G and sum
//	tpiba = tpi / alat;
    for (n = 0;n < ncxyz;n++)  // do n=1,ncxyz
    {
        da[n] = gaux(n, 0) * tpiba;
    } // end do

    return;
} // end subroutine grad_dot

