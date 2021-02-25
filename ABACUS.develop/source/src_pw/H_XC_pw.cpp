void potential::v_xc
(
    double **rho_in,
    double &etxc,
    double &vtxc,
    matrix &v)
{
    TITLE("potential","v_xc");
    timer::tick("potential","v_xc");
    //Exchange-Correlation potential Vxc(r) from n(r)
    etxc = 0.0;
    vtxc = 0.0;

    // the square of the e charge
    // in Rydeberg unit, so * 2.0.
    double e2 = 2.0;

    double rhox = 0.0;
    double arhox = 0.0;
    double zeta = 0.0;
    double ex = 0.0;
    double ec = 0.0;
    double vx[2];
    double vc[2];

    int ir, is;
    int neg [3];

    double vanishing_charge = 1.0e-10;
    if (NSPIN == 1 || ( NSPIN ==4 && !DOMAG && !DOMAG_Z))
    {
        // spin-unpolarized case
        // for parallel : ncxyz change ==> nrxx
        // 2008-06-01 mohan
        for (int ir = 0;ir < pw.nrxx;ir++)
        {
            rhox = rho_in[0][ir] + CHR.rho_core[ir];
            arhox = abs(rhox);
            if (arhox > vanishing_charge)
            {
                // call
                XC_Functional::xc(arhox, ex, ec, vx[0], vc[0]);
                //if(ir<10)
                //{
                //    cout << "\n ir = " << ir << " ex = " << ex << " ec = " << ec;
                //}
                v(0,ir) = e2 * (vx[0] + vc[0]);
                etxc += e2 * (ex + ec) * rhox;
                vtxc += v(0, ir) * rho_in[0][ir];
            } // endif
        } //enddo
    }
    else if(NSPIN ==2)
    {
        // spin-polarized case
        neg [0] = 0;
        neg [1] = 0;
        neg [2] = 0;

        // 2008-06-01 mohan
        if (test_potential>0) cout<<"\n Begin calculate Exc(r) and Vxc(r)";
        for (ir = 0;ir < pw.nrxx;ir++)
        {
            rhox = rho_in[0][ir] + rho_in[1][ir] + CHR.rho_core[ir]; //HLX(05-29-06): bug fixed
            arhox = abs(rhox);

            if (arhox > vanishing_charge)
            {
                zeta = (rho_in[0][ir] - rho_in[1][ir]) / arhox; //HLX(05-29-06): bug fixed

                if (abs(zeta)  > 1.0)
                {
                    ++neg[2];
                    zeta = (zeta > 0.0) ? 1.0 : (-1.0);
                }

                if (rho_in[0][ir] < 0.0)
                {
                    ++neg[0];
                }

                if (rho_in[1][ir] < 0.0)
                {
                    ++neg[1];
                }

                // call
                XC_Functional::xc_spin(arhox, zeta, ex, ec, vx[0], vx[1], vc[0], vc[1]);

                //if(ir<10)
                //{
                //    cout << "\n ir = " << ir << " ex = " << ex << " ec = " << ec;
                //}

                for (is = 0;is < NSPIN;is++)
                {
                    v(is, ir) = e2 * (vx[is] + vc[is]);
                }

                etxc += e2 * (ex + ec) * rhox;

                vtxc += v(0, ir) * rho_in[0][ir] + v(1, ir) * rho_in[1][ir];
            } // endif
        } //   enddo
        if (test_potential>0) cout<<"\n End calculate Exc(r) and Vxc(r) with SPIN == 2";

    } // nspin 2
    else if(NSPIN == 4)//noncollinear case added by zhengdy
    {
        for( ir = 0;ir<pw.nrxx; ir++)
        {
            double amag = sqrt( pow(rho_in[1][ir],2) + pow(rho_in[2][ir],2) + pow(rho_in[3][ir],2) );

            rhox = rho_in[0][ir] + CHR.rho_core[ir];

            if ( rho_in[0][ir] < 0.0 )  neg[0] -= rho_in[0][ir];

            arhox = abs( rhox );

            if ( arhox > vanishing_charge )
            {
                zeta = amag / arhox;

                if ( abs( zeta ) > 1.0 )
                {
                    neg[1] += 1.0 / ucell.omega;

                    zeta = (zeta > 0.0) ? 1.0 : (-1.0);
                }//end if

                XC_Functional::xc_spin( arhox, zeta, ex, ec, vx[0], vx[1], vc[0], vc[1] );
				
                etxc += e2 * ( ex + ec ) * rhox;

                v(0, ir) = e2*( 0.5 * ( vx[0] + vc[0] + vx[1] + vc[1] ) );
                vtxc += v(0,ir) * rho_in[0][ir];

                double vs = 0.5 * ( vx[0] + vc[0] - vx[1] - vc[1] );
                if ( amag > vanishing_charge )
                {
                    for(int ipol = 1;ipol< 4;ipol++)
                    {
                        v(ipol, ir) = e2 * vs * rho_in[ipol][ir] / amag;
                        vtxc += v(ipol,ir) * rho_in[ipol][ir];
                    }//end do
                }//end if
            }//end if
        }//end do
    }//end if
    // energy terms, local-density contributions

    // add gradient corrections (if any)
    // mohan modify 2009-12-15
    GGA_PW::gradcorr(etxc, vtxc, v);

    // parallel code : collect vtxc,etxc
    // mohan add 2008-06-01
    Parallel_Reduce::reduce_double_pool( etxc );
    Parallel_Reduce::reduce_double_pool( vtxc );

    etxc *= ucell.omega / pw.ncxyz;
    vtxc *= ucell.omega / pw.ncxyz;
    if (test_potential > 1)
    {
        OUT(ofs_running,"etxc",etxc);
        OUT(ofs_running,"vtxc",vtxc);
    }

    if (test_potential >0 ) cout<<"\n End calculate gradient";

    if (test_potential > 2)
    {
        cout<<"\n After gradcorr : vtxc = "<<vtxc;
        cout<<"\n After gradcorr : etxc = "<<etxc;
    }
    timer::tick("potential","v_xc");
    return;
}
