#include "tools.h"
#include "global.h"
#include "potential.h"
#include "xc_functional.h"
#include "xc_gga_pw.h"
#include "efield.h"
#include "math.h"
#include "potential_libxc.h"

potential::potential()
{
    vltot = new double[1];
    vrs1 = new double[1];
    this->out_potential = 0;
}

potential::~potential()
{
    delete[] vltot;
    delete[] vrs1;
}

void potential::init(const int nrxx)
{
    if (test_potential) TITLE("potential","init");
    assert(nrxx>0);

    delete[] this->vltot;
    this->vltot = new double[nrxx];
    Memory::record("potential","vltot",nrxx,"double");

    this->vr.create(NSPIN,nrxx);
    this->vrs.create(NSPIN,nrxx);
    Memory::record("potential","vr",NSPIN*nrxx,"double");
    Memory::record("potential","vrs",NSPIN*nrxx,"double");

    delete[] this->vrs1;
    this->vrs1 = new double[nrxx];//mohan add 2007-11-12
    Memory::record("potential","vrs1",nrxx,"double");

    this->vnew.create(NSPIN,nrxx);
    Memory::record("potential","vnew",NSPIN*nrxx,"double");

    return;
}

//----------------------------------------------------------
//  EXPLAIN :
//  Initializes the self consistent potential in the array vr
//  if 'delta_vh=0' & 'vna=1', then 
//  the final potential: atom Hartree + Local pp = Vna
//  else if 'delta_vh=1' & 'vna=0', then
//  the final potential: the remaining Hartree potential.
//----------------------------------------------------------
void potential::init_pot(const int &istep, const bool delta_vh, const bool vna)
{
    TITLE("potential","init_pot");
    timer::tick("potential","init_pot");

    assert(istep>=0);

    //ofs_running << " istep=" << istep << " delta_vh=" << delta_vh << " vna=" << vna << endl;

    vrs.zero_out();

    // mohan fix bug 2011-07-08
    // the vltot should and must be zero here.
    ZEROS(this->vltot, pw.nrxx);

    // vna use this line to get delta_vh and V_Efield to do grid integration.
    // if detal_vh is set to 1, then local pseudopotential is not added.
    if(delta_vh) 
    {
        if(EFIELD && !DIPOLE)
        {
            Efield EFID;
            // in fact, CHR.rho is not used here.
            // if charge correction due to Efield is considered,
            // the structure here need to be updated.

            static bool first = true;
            if(first)
            {
                cout << " ADD THE EFIELD (V/A) : " << Efield::eamp*51.44 << endl;
                first = false;
            }
            EFID.add_efield(CHR.rho[0], this->vltot);	
        }
    }
    else
    {
        // (1) local part of pseudopotentials.
        // set vltot
        this->set_local(this->vltot);

        // mohan fix bug 2011-07-07
        // set pseudopotentials.
        int nspin0=NSPIN;//zhengdy-soc, pauli matrix, just index 0 has vlocal term.
        if(NSPIN==4) nspin0=1;
        for(int is=0; is<nspin0; ++is)
        {
            for(int ir=0; ir<pw.nrxx; ++ir)
            {
                this->vrs(is,ir) = this->vltot[ir];	
            }
        }
    	// (2) core correction potential.
        CHR.set_rho_core( pw.strucFac );

        //if(vna==1)return; // tmp by mohan
    }

	// before ion relaxation
	// need reconstruction -- mohan add 2021-02-09
    if(istep==0)
    {
        OUT(ofs_running,"start_pot",start_pot);

        cout << " START POTENTIAL      : " << start_pot << endl;
        if (this->start_pot == "atomic")//mohan add 2007-10-17
        {
            start_from_atomic:
            CHR.atomic_rho(NSPIN, CHR.rho);
        }
        else if (this->start_pot == "file")
        {
            ofs_running << " try to start potential from file : ";
            for(int is=0; is<NSPIN; is++)
            {
                stringstream ssc;
                ssc << global_out_dir << "SPIN" << is + 1 << "_CHG";
                ofs_running << ssc.str() << endl;
                // mohan update 2012-02-10
                if(CHR.read_rho( is, ssc.str() )) 
                {
                    ofs_running << " Read in the charge density: " << ssc.str() << endl;
				}
				else if(is>0 && NSPIN==4)
				{
					if(PRENSPIN == 1)
					{//read only up+down , others set to zero.
						ofs_running << " Didn't read in the charge density but autoset it for spin " <<is+1<< endl;
						for(int ir=0;ir<pw.nrxx;ir++)
							CHR.rho[is][ir] = 0.0;
					}
					else if(PRENSPIN == 2)
					{//read up and down , then rearrange them.
						if(is==1) 
						{
							WARNING_QUIT("potential::init_pot","Incomplete charge density file!");
						}
						else if(is==2) 
						{
							ofs_running << " Didn't read in the charge density but would rearrange it later. "<< endl;
						}
						else if(is==3)
						{
							ofs_running << " rearrange charge density " << endl;
							for(int ir=0;ir<pw.nrxx;ir++)
							{
								CHR.rho[3][ir] = CHR.rho[0][ir] - CHR.rho[1][ir];
								CHR.rho[0][ir] = CHR.rho[0][ir] + CHR.rho[1][ir];
								CHR.rho[1][ir] = 0.0;
								CHR.rho[2][ir] = 0.0;
							}
						}
					}
					else
					{
						WARNING_QUIT("potential::init_pot","Incomplete charge density file!");
					}
				}
				else
                {
                    ofs_running << " Start charge density from atomic charge density." << endl;
                    goto start_from_atomic;
                }
            }
        }
        else
        {
            WARNING_QUIT("potential::init_pot","start_pot is wrong!");
        }
		
		// Peize Lin add 2020.04.04
		if(restart.info_load.load_charge && !restart.info_load.load_charge_finish)
		{
			for(int is=0; is<NSPIN; ++is)
			{
				restart.load_disk("charge", is);
			}
			restart.info_load.load_charge_finish = true;
		}
    }
    else
    {
        // the extrapolation part moves to ions.cpp.
    }

    CHR.renormalize_rho();

    //-------------------------------------------------------
    // Here we computer the potential which correspond
    // to the charge density
    // if(vna=1)
    //   get the Hartree potential from atomic charge density.
    // else
    //   delta_vh = 1: get the delta Hartree potential.
    //   delta_vh = 0: get the full Hartree potential.
    //--------------------------------------------------------
    this->v_of_rho( CHR.rho, en.ehart, en.etxc, en.vtxc, vr, delta_vh, vna);

    //----------------------------------------------------------
    // Define the total local potential (external+scf) in DFT
	// Define TDDFT potential, by Fuxiang He
    //----------------------------------------------------------
    if(vext == 0) 
	{
		this->set_vrs();
	}
    else 
	{
		this->set_vrs_tddft(istep);
	}

    //figure::picture(this->vrs1,pw.ncx,pw.ncy,pw.ncz);
    timer::tick("potential","init_pot");
    return;
}


//==========================================================
// This routine computes the local potential in real space
// vltot(ir)
//==========================================================
void potential::set_local(double* vl_pseudo)const
{
    TITLE("potential","set_local");
    timer::tick("potential","set_local");

    complex<double> *vg = new complex<double>[pw.ngmc];

    ZEROS( vg, pw.ngmc );

    for (int it=0; it<ucell.ntype; it++)
    {
        for (int ig=0; ig<pw.ngmc; ig++)
        {
            vg[ig] += ppcell.vloc(it, pw.ig2ngg[ig]) * pw.strucFac(it,ig);
        }
    }

    UFFT.ToRealSpace(vg, vl_pseudo); 

    delete[] vg;

    if(EFIELD && !DIPOLE)
    {
        Efield EFID;
        // in fact, CHR.rho is not used here.
        // if charge correction due to Efield is considered,
        // the structure here need to be updated.

        static bool first = true;
        if(first)
        {
            cout << " ADD THE EFIELD (V/A) : " << Efield::eamp*51.44 << endl;
            first = false;
        }
        EFID.add_efield(CHR.rho[0], vl_pseudo);	
    }

    //ofs_running <<" set local pseudopotential done." << endl;
    timer::tick("potential","set_local");
    return;
}

//==========================================================
// This routine computes the Hartree and Exchange and Correlation
// potential and energies which corresponds to a given charge density
// The XC potential is computed in real space, while the
// Hartree potential is computed in reciprocal space.
//==========================================================
void potential::v_of_rho
(
    double **rho_in,
    double &ehart,
    double &etxc,
    double &vtxc,
    matrix &v_in,
    const bool delta_vh,
    const bool vna
)
{
    TITLE("potential","v_of_rho");
    v_in.zero_out();

    if(vna)
    {
        timer::tick("potential","vna_h");
        double** rho_atom = new double*[NSPIN];
        for(int is=0; is<NSPIN; is++)
        {
            rho_atom[is] = new double[pw.nrxx];
            ZEROS(rho_atom[is], pw.nrxx);
        }

        CHR.atomic_rho(NSPIN,rho_atom);

        this->v_h(NSPIN, ehart, v_in, rho_atom);

        for(int is=0; is<NSPIN; is++)
        {
            delete[] rho_atom[is];
        }
        delete[] rho_atom;

        timer::tick("potential","vna_h");
        return;
    }

    timer::tick("potential","v_of_rho",'E');

//----------------------------------------------------------
//  calculate the exchange-correlation potential
//----------------------------------------------------------
	
	#ifdef TEST_LIBXC
    Potential_Libxc::v_xc(rho_in, etxc, vtxc, v_in);
	#else
    this->v_xc(rho_in, etxc, vtxc, v_in);
	#endif

//----------------------------------------------------------
//  calculate the Hartree potential
//----------------------------------------------------------
    if(delta_vh)
    {
        //--------------------------------------------
        // get the atomic charge on real space grid.	
        //--------------------------------------------
        double** rho_atom = new double*[NSPIN];
        for(int is=0; is<NSPIN; is++)
        {
            rho_atom[is] = new double[pw.nrxx];
            ZEROS(rho_atom[is], pw.nrxx);
        }

        CHR.atomic_rho(NSPIN,rho_atom);

        //--------------------------------------------
        // get the delta atomic charge on real space
        // grid.
        //--------------------------------------------
        for(int is=0; is<NSPIN; is++)
        {
            for(int ir=0; ir<pw.nrxx; ir++)
            {
                rho_atom[is][ir] = rho_in[is][ir] - rho_atom[is][ir];
            }
        }

        //--------------------------------------------
        // get the potential from the delta charge
        //--------------------------------------------
        this->v_h(NSPIN, ehart, v_in, rho_atom);

        for(int is=0; is<NSPIN; is++)
        {
            delete[] rho_atom[is];
        }
        delete[] rho_atom;
    }
    else
    {
        this->v_h(NSPIN, ehart, v_in, rho_in);
    }

    // mohan add 2011-06-20
    if(EFIELD && DIPOLE)
    {
        Efield EFID;
        for (int is = 0;is < NSPIN;is++)
        {
            EFID.add_efield(rho_in[is], &v_in.c[is*pw.nrxx]);
        }
    }
    timer::tick("potential","v_of_rho",'E');
    return;
} //end subroutine v_of_rho

//--------------------------------------------------------------------
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

//--------------------------------------------------------------------
// Transform charge density to hartree potential.
//--------------------------------------------------------------------
void potential::v_h(int NSPIN,double &ehart, matrix &v, double** rho)
{
    TITLE("potential","v_h");
    timer::tick("potential","v_hartree");

    complex<double> *Porter = UFFT.porter;

    //  Hartree potential VH(r) from n(r)
    ZEROS( Porter, pw.nrxx );
    int nspin0 = 1;
    if(NSPIN==2)nspin0 = NSPIN;
    for(int is=0; is<nspin0; is++)
    {
        for (int ir=0; ir<pw.nrxx; ir++) 
        {
            Porter[ir] += complex<double>( rho[is][ir], 0.0 );
        }
    }

    //=============================
    //  bring rho (aux) to G space
    //=============================
    pw.FFT_chg.FFT3D(Porter, -1);

    //double charge;
    //if (pw.gstart == 1)
    //{
    //    charge = ucell.omega * Porter[pw.ig2fftc[0]].real();
    //}
    //OUT(ofs_running, "v_h charge", charge);

    //=======================================================
    // calculate hartree potential in G-space (NB: V(G=0)=0 )
    //=======================================================
    ehart = 0.0;

    complex<double> *vh_g  = new complex<double>[pw.ngmc];
    ZEROS(vh_g, pw.ngmc);

    for (int ig = pw.gstart; ig<pw.ngmc; ig++)
    {
        const int j = pw.ig2fftc[ig];
        if(pw.gg[ig] >= 1.0e-12) //LiuXh 20180410
        {
            const double fac = e2 * FOUR_PI / (ucell.tpiba2 * pw.gg [ig]);

            ehart += ( conj( Porter[j] ) * Porter[j] ).real() * fac;
            vh_g[ig] = fac * Porter[j];
        }
    }

    Parallel_Reduce::reduce_double_pool( ehart );
    ehart *= 0.5 * ucell.omega;

    //cout << " ehart=" << ehart << endl;

    ZEROS(Porter, pw.nrxx);

    for (int ig = 0;ig < pw.ngmc;ig++)
    {
        Porter[pw.ig2fftc[ig]] = vh_g[ig];
    }

    //==========================================
    //transform hartree potential to real space
    //==========================================
    pw.FFT_chg.FFT3D(Porter, 1);
    //==========================================
    //Add hartree potential to the xc potential
    //==========================================

	if(NSPIN==4)
		for (int ir = 0;ir < pw.nrxx;ir++)
		{
			v(0, ir) += Porter[ir].real();
		}
	else
    	for (int is = 0;is < NSPIN;is++)
    	{
        	for (int ir = 0;ir < pw.nrxx;ir++)
        	{
            		v(is, ir) += Porter[ir].real();
        	}
    	}

	//-------------------------------------------
	// output the Hartree potential into a file.
	//-------------------------------------------
	if(out_potential==-2)
	{
		cout << " output VH" << endl;
		int is = 0;
		int iter = 0;
		int precision = 3;
		string fn = "VH.dat";
		stringstream ss;
		ss << global_out_dir << fn;
		matrix v;
		v.create(1,pw.nrxx);
		for(int ir=0; ir<pw.nrxx; ++ir)
		{
			v(0,ir) = Porter[ir].real();
		}
		this->write_potential( is, iter, ss.str(), v, precision, 1 );
	}

    timer::tick("potential","v_hartree");
    delete[] vh_g;
    return;
} // end subroutine v_h


//==========================================================
// MEMBER FUNCTION : set_vrs
// set the total local potential vrs on the smooth mesh to
// be used in h_psi, adding the (spin dependent) scf (H+xc)
// part and the sum of all the local pseudopotential
// contributions.
//==========================================================
// from set_vrs.f90
void potential::set_vrs(void)
{
    TITLE("potential","set_vrs");
    timer::tick("potential","set_vrs");

    for (int is = 0;is < NSPIN;is++)
    {
        //=================================================================
        // define the total local potential (external + scf) for each spin
        //=================================================================
		if(NSPIN==4&&is>0)
			for (int i = 0;i < pw.nrxx;i++)
			{
				this->vrs(is, i) = this->vr(is, i);
			}
		else        
			for (int i = 0;i < pw.nrxx;i++)
	        {
	            this->vrs(is, i) = this->vltot[i] + this->vr(is, i);
	//	    	cout <<"i: "<< i <<"	"<< "vrs: " << vrs(is,i) <<endl;
			}

        //====================================================
        // add external linear potential, fuxiang add in 2017/05
        //====================================================
/*
        if (istep >= 0)
        {
            for (int i = 0;i < pw.nrxx;i++)
            {
                this->vrs(is, i) = this->vltot[i] + this->vr(is, i);
            }
        }
        else
        {
            this->vext = new double[pw.nrxx];
            const int yz = pw.ncy*pw.nczp;
            int index, i, j, k;

            for(int ir=0; ir<pw.nrxx; ++ir)
            {
                index = ir;
                i     = index / yz; // get the z, z is the fastest
                index = index - yz * i;// get (x,y)
                j     = index / pw.nczp;// get y
                k     = index - pw.nczp*j + pw.nczp_start;// get x

                //if (k<pw.ncx*0.1) this->vext[ir] = -0.4*k/pw.ncx+0.05;
                //else if (k>=pw.ncx*0.1 && k<pw.ncx*0.9) this->vext[ir] = 0.1*k/pw.ncx;
                //else if (k>=pw.ncx*0.9) this->vext[ir] = -0.4*(1.0*k/pw.ncx-1)+0.05;

                if (k<pw.ncx*0.05) this->vext[ir] = (0.019447*k/pw.ncx-0.001069585)*ucell.lat0;
                else if (k>=pw.ncx*0.05 && k<pw.ncx*0.95) this->vext[ir] = -0.0019447*k/pw.ncx*ucell.lat0;
                else if (k>=pw.ncx*0.95) this->vext[ir] = (0.019447*(1.0*k/pw.ncx-1)-0.001069585)*ucell.lat0;

                this->vrs(is,ir) = this->vltot[ir] + this->vr(is, ir) + this->vext[ir];

                //cout << "x: " << k <<"	" << "y: " << j <<"	"<< "z: "<< i <<"	"<< "ir: " << ir << endl;
                //cout << "vext: " << this->vext[ir] << endl;
                //cout << "vrs: " << vrs(is,ir) <<endl;
            }
        }
*/

    }


    timer::tick("potential","set_vrs");
    return;
} //end subroutine set_vrs

// ----------------------------------------------------------------------
void potential::newd(void)
{
    if (test_potential) TITLE("potential","newd");

    // distringuish non-local pseudopotential in REAL or RECIPROCAL space.
    // if in real space, call new_r
    // if in reciprocal space, call new_g

    // new g:
    //----------------------------------------------------------------------
    //  This routine computes the integral of the effective potential with
    //  the Q function and adds it to the bare ionic D term which is used
    //  to compute the non-local term in the US scheme.

	// no ultrasoft potentials: use bare coefficients for projectors
	// if( spin_orbital) ....
	// else if(noncolin) ....
	for (int iat=0; iat<ucell.nat; iat++)
	{
		const int it = ucell.iat2it[iat];
		const int nht = ucell.atoms[it].nh;
		// nht: number of beta functions per atom type
		for (int is = 0; is < NSPIN; is++)
		{
			for (int ih=0; ih<nht; ih++)
			{
				for (int jh=ih; jh<nht; jh++)
				{
					if(LSPINORB)
					{
						ppcell.deeq_nc(is , iat , ih , jh)= ppcell.dvan_so(is , it , ih , jh);
						ppcell.deeq_nc(is , iat , jh , ih)= ppcell.dvan_so(is , it , jh , ih);
					}
					else if( NSPIN==4 )
					{
						if(is==0)
						{
							ppcell.deeq_nc(is, iat, ih, jh) = ppcell.dvan(it, ih, jh);
							ppcell.deeq_nc(is, iat, jh, ih) = ppcell.dvan(it, ih, jh);
						}
						else if(is==1)
						{
							ppcell.deeq_nc(is, iat, ih, jh) = complex<double>(0.0 , 0.0);
							ppcell.deeq_nc(is, iat, jh, ih) = complex<double>(0.0 , 0.0);
						}
						else if(is==2)
						{
							ppcell.deeq_nc(is, iat, ih, jh) = complex<double>(0.0 , 0.0);
							ppcell.deeq_nc(is, iat, jh, ih) = complex<double>(0.0 , 0.0);
						}
						else if(is==3)
						{
							ppcell.deeq_nc(is, iat, ih, jh) = ppcell.dvan(it, ih, jh);
							ppcell.deeq_nc(is, iat, jh, ih) = ppcell.dvan(it, ih, jh);
						}
					}
					else{
						ppcell.deeq(is, iat, ih, jh) = ppcell.dvan(it, ih, jh);
						ppcell.deeq(is, iat, jh, ih) = ppcell.dvan(it, ih, jh);
					}
				}
			}
		}
	}
	return;
} // end subroutine newd

void potential::print_pot(ofstream &ofs)const
{
    //  ofs << "\n potinit() : ";
    ofs << "\n vtxc       = " << en.vtxc;
    ofs << "\n etxc       = " << en.etxc;
    ofs << "\n ehart      = " << en.ehart;
    //   ofs << "\n charge     = " << charge;
    out.printr1_d(ofs, " vltot : ", vltot, pw.nrxx);
    out.printrm(ofs, "  vrs : ", vrs);
}


//==========================================================
// this function aims to add external time-dependent potential 
// (eg: linear potential) used in tddft
// fuxiang add in 2017-05
//==========================================================
void potential::set_vrs_tddft(const int istep)
{
    TITLE("potential","set_vrs_tddft");
    timer::tick("potential","set_vrs_tddft");

    for (int is = 0;is < NSPIN;is++)
    {
        //====================================================
        // add external linear potential, fuxiang add in 2017/05
        //====================================================

        const int timescale = 1;  // get the time that vext influences;
        if (istep >= timescale)
        {
            for (int i = 0;i < pw.nrxx;i++)
            {
                this->vrs(is, i) = this->vltot[i] + this->vr(is, i);
            }
            cout << "vext = 0! " << endl;
        }
        else
        {
            this->vextold = new double[pw.nrxx];
            this->vext = new double[pw.nrxx];
            const int yz = pw.ncy*pw.nczp;
            int index, i, j, k;

            for(int ir=0; ir<pw.nrxx; ++ir)
            {
                index = ir;
                i     = index / yz; // get the z, z is the fastest
                index = index - yz * i;// get (x,y)
                j     = index / pw.nczp;// get y
                k     = index - pw.nczp*j + pw.nczp_start;// get x

                if(vext_dire == 1)
                {
                    if (k<pw.ncx*0.05) this->vextold[ir] = (0.019447*k/pw.ncx-0.001069585)*ucell.lat0;
                    else if (k>=pw.ncx*0.05 && k<pw.ncx*0.95) this->vextold[ir] = -0.0019447*k/pw.ncx*ucell.lat0;
                    else if (k>=pw.ncx*0.95) this->vextold[ir] = (0.019447*(1.0*k/pw.ncx-1)-0.001069585)*ucell.lat0;
                }
                else if(vext_dire == 2)
                {
                    if (j<pw.ncx*0.05) this->vextold[ir] = (0.019447*j/pw.ncx-0.001069585)*ucell.lat0;
                    else if (j>=pw.ncx*0.05 && j<pw.ncx*0.95)	this->vextold[ir] = -0.0019447*j/pw.ncx*ucell.lat0;
                    else if (j>=pw.ncx*0.95) this->vextold[ir] = (0.019447*(1.0*j/pw.ncx-1)-0.001069585)*ucell.lat0;
                }
                else if(vext_dire == 3)
                {
                    if (i<pw.ncx*0.05) this->vextold[ir] = (0.019447*i/pw.ncx-0.001069585)*ucell.lat0;
                    else if (i>=pw.ncx*0.05 && i<pw.ncx*0.95) this->vextold[ir] = -0.0019447*i/pw.ncx*ucell.lat0;
                    else if (i>=pw.ncx*0.95) this->vextold[ir] = (0.019447*(1.0*i/pw.ncx-1)-0.001069585)*ucell.lat0;
                }

                // Gauss
/*
                const double w = 22.13;    // eV
                const double sigmasquare = 6836;
                const double timecenter = 700;
                const double timenow = (istep-timecenter)*INPUT.md_dt*41.34;
                this->vext[ir] = this->vextold[ir]*cos(w/27.2116*timenow)*exp(-timenow*timenow*0.5/(sigmasquare))*0.25;  //0.1 is modified in 2018/1/12
*/

                //HHG of H atom
/*
                if(istep < 1875)
                {
                    this->vext[ir] = this->vextold[ir]*2.74*istep/1875*cos(0.0588*istep*INPUT.md_dt*41.34);	// 2.75 is equal to E0;
                }
                else if(istep < 5625)
                {
                    this->vext[ir] = this->vextold[ir]*2.74*cos(0.0588*istep*INPUT.md_dt*41.34);
                }
                else if(istep < 7500)
                {
                    this->vext[ir] = this->vextold[ir]*2.74*(7500-istep)/1875*cos(0.0588*istep*INPUT.md_dt*41.34);
                }
*/

                //HHG of H2

                //const double timenow = (istep)*INPUT.md_dt*41.34;
                //this->vext[ir] = this->vextold[ir]*2.74*cos(0.856*timenow)*sin(0.0214*timenow)*sin(0.0214*timenow);
                //this->vext[ir] = this->vextold[ir]*2.74*cos(0.856*timenow)*sin(0.0214*timenow)*sin(0.0214*timenow)*0.01944;
                //this->vext[ir] = this->vextold[ir]*2.74*cos(0.0428*timenow)*sin(0.00107*timenow)*sin(0.00107*timenow);

                this->vrs(is,ir) = this->vltot[ir] + this->vr(is, ir) + this->vext[ir];

                //cout << "x: " << k <<"	" << "y: " << j <<"	"<< "z: "<< i <<"	"<< "ir: " << ir << endl;
                //cout << "vext: " << this->vext[ir] << endl;
                //cout << "vrs: " << vrs(is,ir) <<endl;
            }
            cout << "vext is existed!" << endl;

            delete[] this->vextold;
            delete[] this->vext;
        }
    }


    timer::tick("potential","set_vrs_tddft");
    return;
} //end subroutine set_vrs_tddft
