#include "tools.h"
#include "global.h"
#include "potential.h"
#include "xc_functional.h"
#include "gga_pw.h"
#include "efield.h"
#include "math.h"

potential::potential()
{
    this->test = 0;
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

// from PW/potinit.f90
//----------------------------------------------------------
//  EXPLAIN :
//  This routine initializes the self consistent potential
//  in the array vr
//  if 'delta_vh=0' & 'vna=1', then 
//  the final potential: atom Hartree + Local pp = Vna
//  else if 'delta_vh=1' & 'vna=0', then
//  the final potential: remain Hartree potential.
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
            // in fact, chr.rho is not used here.
            // if charge correction due to Efield is considered,
            // the structure here need to be updated.

            static bool first = true;
            if(first)
            {
                cout << " ADD THE EFIELD (V/A) : " << Efield::eamp*51.44 << endl;
                first = false;
            }
            EFID.add_efield(chr.rho[0], this->vltot);	
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
        chr.set_rho_core( pw.strucFac );

        //if(vna==1)return; // tmp by mohan
    }

    if(istep==0)//not begin to do ion relaxation.
    {
        OUT(ofs_running,"start_pot",start_pot);

        cout << " START POTENTIAL      : " << start_pot << endl;
        if (this->start_pot == "atomic")//mohan add 2007-10-17
        {
            start_from_atomic:
            chr.atomic_rho(NSPIN, chr.rho);
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
                if(chr.read_rho( is, ssc.str() )) 
                {
                    ofs_running << " Read in the charge density: " << ssc.str() << endl;
                }
		else if(is>0 && NSPIN==4)
		{
			if(PRENSPIN == 1)
			{//read only up+down , others set to zero.
				ofs_running << " Didn't read in the charge density but autoset it for spin " <<is+1<< endl;
				for(int ir=0;ir<pw.nrxx;ir++)
					chr.rho[is][ir] = 0.0;
			}
			else if(PRENSPIN == 2)
			{//read up and down , then rearrange them.
				if(is==1) WARNING_QUIT("potential::init_pot","Incomplete charge density file!");
				if(is==2) 
				{
					ofs_running << " Didn't read in the charge density but would rearrange it later. "<< endl;
				}
				if(is==3)
				{
					ofs_running << " rearrange charge density " << endl;
					for(int ir=0;ir<pw.nrxx;ir++)
					{
						chr.rho[3][ir] = chr.rho[0][ir] - chr.rho[1][ir];
						chr.rho[0][ir] = chr.rho[0][ir] + chr.rho[1][ir];
						chr.rho[1][ir] = 0.0;
						chr.rho[2][ir] = 0.0;
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
				restart.load_disk("charge", is);
			restart.info_load.load_charge_finish = true;
		}
    }
    else
    {
        // the extrapolation part moves to ions.cpp.
    }

    chr.renormalize_rho();

    //-------------------------------------------------------
    // Here we computer the potential which correspond
    // to the charge density
    // if(vna=1)
    //   get the Hartree potential from atomic charge density.
    // else
    //   delta_vh = 1: get the delta Hartree potential.
    //   delta_vh = 0: get the full Hartree potential.
    //--------------------------------------------------------
    this->v_of_rho( chr.rho, en.ehart, en.etxc, en.vtxc, vr, delta_vh, vna);

    //----------------------------------------------------------
    // Define the total local potential (external+scf)
    //----------------------------------------------------------
    if(vext == 0) this->set_vrs(pw.doublegrid);
    else this->set_vrs_tddft(pw.doublegrid, istep);

    //figure::picture(this->vrs1,pw.ncx,pw.ncy,pw.ncz);
    timer::tick("potential","init_pot");
    return;
}

// from setlocal.f90
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
        // in fact, chr.rho is not used here.
        // if charge correction due to Efield is considered,
        // the structure here need to be updated.

        static bool first = true;
        if(first)
        {
            cout << " ADD THE EFIELD (V/A) : " << Efield::eamp*51.44 << endl;
            first = false;
        }
        EFID.add_efield(chr.rho[0], vl_pseudo);	
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

        chr.atomic_rho(NSPIN,rho_atom);

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
    this->v_xc(rho_in, etxc, vtxc, v_in);

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

        chr.atomic_rho(NSPIN,rho_atom);

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
            rhox = rho_in[0][ir] + chr.rho_core[ir];
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
            rhox = rho_in[0][ir] + rho_in[1][ir] + chr.rho_core[ir]; //HLX(05-29-06): bug fixed
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

            rhox = rho_in[0][ir] + chr.rho_core[ir];

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

                double vs = 0.5 * ( vx[0] + vc[0] - vx[1] - vc[1] );
                v(0, ir) = e2*( 0.5 * ( vx[0] + vc[0] + vx[1] + vc[1] ) );

                if ( amag > vanishing_charge )
                {
                    for(int ipol = 1;ipol< 4;ipol++)
                    {
                        v(ipol, ir) = e2 * vs * rho_in[ipol][ir] / amag;
                        vtxc += v(ipol,ir) * rho_in[ipol][ir];
                    }//end do
                }//end if
                etxc += e2 * ( ex + ec ) * rhox;
                vtxc += v(0,ir) * rho_in[0][ir];
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

// translate from write_rho in charge.cpp.
void potential::write_potential(const int &is, const int &iter, const string &fn, 
const matrix &v, const int &precision, const int &hartree)const
{
    TITLE("potential","write_potential");

    if(out_potential==0) 
    {
        return;
    }
    else if(out_potential<0)
    {
        if(hartree==0) return;
    }
    else if(iter % out_potential != 0)
    {
        return;
    }
    timer::tick("potential","write_potential");

    ofstream ofs;

    if(MY_RANK==0)
    {
        ofs.open( fn.c_str() );

        ofs << ucell.latName << endl;//1
        ofs << " " << ucell.lat0 * 0.529177 << endl;
        ofs << " " << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
        ofs << " " << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
        ofs << " " << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;

        for(int it=0; it<ucell.ntype; it++)
        {
            ofs << " " << ucell.atoms[it].label;
        }
        ofs << endl;
        for(int it=0; it<ucell.ntype; it++)
        {
            ofs << " " << ucell.atoms[it].na;
        }
        ofs << endl;
        ofs << "Direct" << endl;

        for(int it=0; it<ucell.ntype; it++)
        {
            for(int ia=0; ia<ucell.atoms[it].na; ia++)
            {
                ofs << " " << ucell.atoms[it].taud[ia].x
                    << " " << ucell.atoms[it].taud[ia].y
                    << " " << ucell.atoms[it].taud[ia].z << endl;
            }
        }
        ofs << pw.ncx << " " << pw.ncy << " " << pw.ncz;
        ofs << setprecision(precision);
        ofs << scientific; 
        if(!ofs)
        {
            WARNING("potential::write_potential","Can't create VHartree File!");
        }
    }	

#ifndef __MPI
    int count=0;
    for(int k=0; k<pw.ncz; k++)
    {
        ofs << "\n" << k << " iz";
        double value = 0.0;
        double ave = 0.0;
        for(int j=0; j<pw.ncy; j++)
        {
            for(int i=0; i<pw.ncx; i++)
            {
                if(count%8==0) ofs << "\n";
                value = v(is, i*pw.ncy*pw.ncz + j*pw.ncz + k);
                ofs << " " << value;
                ave += value;
                ++count;
            }
        }
        ofs << "\n" << ave/pw.ncx/pw.ncy << " average";
    }
#else
    MPI_Barrier(MPI_COMM_WORLD);
    // only do in the first pool.
    if(MY_POOL==0)
    {
        // num_z: how many planes on processor 'ip'
        int *num_z = new int[NPROC_IN_POOL];
        ZEROS(num_z, NPROC_IN_POOL);
        for (int iz=0;iz<pw.ncz;iz++)
        {
            int ip = iz % NPROC_IN_POOL;
            num_z[ip]++;
        }

        // start_z: start position of z in
        // processor ip.
        int *start_z = new int[NPROC_IN_POOL];
        ZEROS(start_z, NPROC_IN_POOL);
        for (int ip=1;ip<NPROC_IN_POOL;ip++)
        {
            start_z[ip] = start_z[ip-1]+num_z[ip-1];
        }

        // which_ip: found iz belongs to which ip.
        int *which_ip = new int[pw.ncz];
        ZEROS(which_ip, pw.ncz);
        for(int iz=0; iz<pw.ncz; iz++)
        {
            for(int ip=0; ip<NPROC_IN_POOL; ip++)
            {
                if(iz>=start_z[NPROC_IN_POOL-1])
                {
                    which_ip[iz] = NPROC_IN_POOL-1;
                    break;
                }
                else if(iz>=start_z[ip] && iz<start_z[ip+1])
                {
                    which_ip[iz] = ip;
                    break;
                }
            }
            //ofs_running << "\n iz=" << iz << " ip=" << which_ip[iz];
        }
        int count=0;
        int nxy = pw.ncx * pw.ncy;
        double* zpiece = new double[nxy];
        // save the rho one z by one z.
        for(int iz=0; iz<pw.ncz; iz++)
        {
            //ofs_running << "\n" << iz << " iz"; //LiuXh modify 20200624
            // tag must be different for different iz.
            ZEROS(zpiece, nxy);
            int tag = iz;
            MPI_Status ierror;

            // case 1: the first part of rho in processor 0.
            if(which_ip[iz] == 0 && RANK_IN_POOL ==0)
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v(is, ir*pw.nczp+iz-start_z[RANK_IN_POOL] );
                    //ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*pw.nczp+iz=" << ir*pw.nczp+iz;
                }
            }
            // case 2: > first part rho: send the rho to
            // processor 0.
            else if(which_ip[iz] == RANK_IN_POOL )
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v(is, ir*pw.nczp+iz-start_z[RANK_IN_POOL]);
                }
                MPI_Send(zpiece, nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
            }

            // case 2: > first part rho: processor 0 receive the rho
            // from other processors
            else if(RANK_IN_POOL==0)
            {
                MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
                //ofs_running << "\n Receieve First number = " << zpiece[0];
            }

            // write data
            if(MY_RANK==0)
            {
                //ofs << "\niz=" << iz;
                double ave = 0.0;
                for(int ir=0; ir<nxy; ir++)
                {
                    if(count%8==0) ofs << "\n";
                    ofs << " " << zpiece[ir];
                    ave += zpiece[ir];
                    ++count;
                }
                ofs << "\n" << ave/nxy << " average"; 
            }
        }
        delete[] zpiece;
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(MY_RANK==0) ofs.close();
    timer::tick("potential","write_potential");
    return;
}

//==========================================================
// MEMBER FUNCTION : set_vrs
// set the total local potential vrs on the smooth mesh to
// be used in h_psi, adding the (spin dependent) scf (H+xc)
// part and the sum of all the local pseudopotential
// contributions.
//==========================================================
// from set_vrs.f90
#ifdef __FP
#include "../src_lcao/bfield.h"
#endif
void potential::set_vrs(const bool doublegrid)
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
            //cout <<"i: "<< i <<"	"<< "vrs: " << vrs(is,i) <<endl;
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

        //====================================================
        // And interpolate it on the smooth mesh if necessary
        //====================================================
        /*
        if (doublegrid)
        {
            double *vrs_1d = new double[pw.nrxx]();
            assert(vrs_1d != 0);
            dcopy(vrs, is, vrs_1d);
            //interpolate (vrs (1, is), vrs (1, is), - 1);
            interpolate(vrs_1d, vrs_1d, -1);
            //(vrs ( is, 0), vrs ( is, 0), - 1);
            delete [] vrs_1d;
        }
        */
    }

#ifdef __FP
    if(BFIELD)
    {
        if(NSPIN==2)
        {
            bfid.add_zeeman();
        }
    }
#endif

    timer::tick("potential","set_vrs");
    return;
} //end subroutine set_vrs

// ----------------------------------------------------------------------
void potential::newd()
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
    if (!ppcell.okvan)
    {
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
    }

    //====================
    // if use PAW method
    //====================
    /*   if (okpaw) then
     ! Add paw contributions to deeq (computed in paw_potential)
     do na=1,nat
        nt = ityp(na)
        IF (.not.upf(nt)%tpawp) cycle
        ijh=0
        do ih=1,nh(nt)
           do jh=ih,nh(nt)
              ijh=ijh+1
              deeq(ih,jh,na,1:nspin) = deeq(ih,jh,na,1:nspin) &
                                     + ddd_paw(ijh,na,1:nspin)
              deeq(jh,ih,na,1:nspin) = deeq(ih,jh,na,1:nspin)
           end do
        end do
     end do
    end IF
    */

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

/*

// from mix_potential.f90
//-----------------------------------------------------------------------
void potential::mix_potential (int ndim, double *vout, double *vin, double alphamix,
                           double &dr2, double tr2, int iter, int n_iter, string filename,
                           bool &conv)
{
    cout << "\n ==== mix_potential() ==== ";
    //-----------------------------------------------------------------------
    //
    // Modified Broyden"s method for potential/charge density mixing
    //             D.D.Johnson, PRB 38, 12807 (1988)
    // On input :
    //    ndim      dimension of arrays vout, vin
    //    vout      output potential/rho at current iteration
    //    vin       potential/rho at previous iteration
    //    alphamix  mixing factor (0 < alphamix <= 1)
    //    tr2       threshold for selfconsistency
    //    iter      current iteration number
    //    n_iter    number of iterations used in the mixing
    //    filename  if present save previous iterations on file "filename"
    //              otherwise keep everything in memory
    // On output:
    //    dr2       [(vout-vin)/ndim]^2
    //    vin       mixed potential
    //    vout      vout-vin
    //    conv      true if dr2 <= tr2

    //   First the dummy variables
    // character (len=42) :: filename
    // integer :: ndim, iter, n_iter
    // real(kind=DP) :: vout (ndim), vin (ndim), alphamix, dr2, tr2
    // logical :: conv

    //   Here the local variables

    // max number of iterations used in mixing: n_iter must be .le. maxter
    // parameter :
    int maxter = 8;

    int n, i, j, *iwork, info, iter_used,
    ipos, inext, ndimtot;
    iwork = new int[maxter];
    assert(iwork != 0);
    // work space containing info from previous iterations:
    // must be kept in memory and saved between calls if filename=' '

    // allocatable, save ::
    matrix  df, dv ;
    double *vinsave, *ework;
    ework = new double[maxter];
    assert(ework != 0);
    matrix beta (maxter, maxter);
    double gamma, norm;

    // adjustable parameters as suggested in the original paper
    double *w, w0;
    w = new double[maxter];
    // data :
    w0 = 0.010 ;
    // w / maxter * 1.0 / ?

    double *work;
    work = new double[maxter];//?
    assert(work != 0);
    start_clock ("mix_pot");

    if (iter < 1)
        cout << "\n mix_potential, iter is wrong";
    if (nmix > maxter)
        cout << "\n mix_potential, nmix too big";
    if (ndim <= 0)
        cout << "\n mix_potential, ndim <= 0";

    for(n=0;n<ndim;n++){ //	do n = 1, ndim
        vout [n] -= vin [n];
    } //  enddo
    // dr2 = dnrm2 (ndim, vout, 1) **2;
    dr2 = dnrm2 (ndim, vout, 1) ;
    dr2 *= dr2;
    ndimtot = ndim;
#ifdef __PARA
    reduce (1, dr2);
    ireduce (1, ndimtot);
#endif
    dr2 = (sqrt (dr2) / ndimtot);	// **2;
    dr2 *= dr2;
    conv_elec = dr2 < tr2;
    if (iter == 1){
        matrix df( nmix, ndim );
        matrix dv( nmix, ndim );
    } //  endif
    if (conv_elec) {
        dv.freemem();
        df.freemem();
        stop_clock ("mix_pot");
        return;
    } // endif
    vinsave = new double [ndim];

    // iter_used = iter-1  if iter <= n_iter
    // iter_used = n_iter  if iter >  n_iter
    iter_used = (iter - 1 < nmix)?(iter-1):nmix;	// min
    // ipos is the position in which results from the present iteraction
    // are stored. ipos=iter-1 until ipos=n_iter, then back to 1,2,...

    ipos = iter - 1 - ( (iter - 2) / nmix) * nmix;
    if (iter > 1){
        for(n=0;n<ndim;n++){ //  do n = 1, ndim
            df (ipos, n) = vout [n] - df (ipos, n);
            dv (ipos, n) = vin  [n] - dv (ipos, n);
        } //  enddo
        // norm = (dnrm2 (ndim, df (1, ipos), 1) ) **2;
        norm = dnrm2 (df,ipos);
        norm *= norm;
#ifdef __PARA
        reduce (1, norm);
#endif
        norm = sqrt (norm);
        dscal (1.0 / norm, df,ipos);	// dscal(ndim,1.0/norm,df (1, ipos), 1);
        dscal (1.0 / norm, dv,ipos);	// dscal(ndim,1.0/norm,dv (1, ipos), 1);
    } //  endif

    dcopy (ndim, vin, 1, vinsave, 1);

    for(i=0;i<iter_used;i++){ // do i = 1, iter_used
        for(j=i+1;j<iter_used;j++){ //  do j = i + 1, iter_used
            beta (i, j) = w [i] * w [j] * ddot (df,j, df,i);
                            //ddot(ndim, df (1, j), 1, df (1, i), 1);
#ifdef __PARA
            reduce (1, beta (i, j) );
#endif
        } // enddo
        beta (i, i) = w0 * w0 + w [i] * w[i];
    } // enddo

    dsytrf ('U', iter_used, beta, maxter, iwork, work, maxter, info);
    cout << "\n broyden, factorization, info: "
        << info;
    dsytri ('U', iter_used, beta, maxter, iwork, work, info);
    cout << "\n broyden, dsytri, info : "
        << info;
    for(i=0;i<iter_used;i++){ // do i = 1, iter_used
        for(j=i+1;j<iter_used;j++){ // do j = i + 1, iter_used
            beta (j, i) = beta (i, j);
        } //enddo
    } //  enddo
    for(i=0;i<iter;i++){ // do i = 1, iter_used;
        work [i] = ddot (df,i,vout);	//ddot(ndim,df (1, i), 1, vout, 1);
    } // enddo
#ifdef __PARA
    reduce (iter_used, work);
#endif
    for(n=0;n<ndim;n++){ // do n = 1, ndim
        vin [n] = vin [n] + mixing_beta * vout [n];
    } // enddo
    for(i=0;i<iter_used;i++){ // do i = 1, iter_used
        gamma = 0.0;
        for(j=0;j<iter_used;j++){ // do j = 1, iter_used
            gamma = gamma + beta (i, j) * w [j] * work [j];
        } // enddo
        for(n=0;n<ndim;n++){ // do n = 1, ndim
            vin [n] = vin [n] - w [i] * gamma * (mixing_beta * df (i, n) +
                dv (i, n) );
        } // enddo
    } // enddo

    inext = iter - ( (iter - 1) / nmix) * nmix;
    //dcopy (ndim, vout,    df,inext);	//df(1, inext), 1);
    dcopy (vout,    df, inext);	//df(1, inext), 1);
    //dcopy (ndim, vinsave, dv,inext);	//dv(1, inext), 1);
    dcopy (vinsave, dv, inext);	//dv(1, inext), 1);

    delete [] vinsave;
    stop_clock ("mix_pot");
    return;
} //end subroutine mix_potential
*/

double potential::vr_ave(const int n,const int size,const double *p)
{
    // return the average of potential here <V^n>
    // [sum_i v(i)^n] / size of the array
    double ave = 0.0;
    if (size < 0)
    {
        cout << "\n The size of the array should > = 1, error";
        exit(0);
    }
    for (int i =0; i< size; i++) ave += pow(p[i], n);
    return ( ave / static_cast<double>(size) ) ;
}


//==========================================================
// this function aims to add external time-dependent potential (eg: linear potential)
// used in tddft
// fuxiang add in 2017-05
//==========================================================

void potential::set_vrs_tddft(const bool doublegrid, const int istep)
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

        //====================================================
        // And interpolate it on the smooth mesh if necessary
        //====================================================
        /*
        if (doublegrid)
        {
            double *vrs_1d = new double[pw.nrxx]();
            assert(vrs_1d != 0);
            dcopy(vrs, is, vrs_1d);
            //interpolate (vrs (1, is), vrs (1, is), - 1);
            interpolate(vrs_1d, vrs_1d, -1);
            //(vrs ( is, 0), vrs ( is, 0), - 1);
            delete [] vrs_1d;
        }
        */
    }

#ifdef __FP
    if(BFIELD)
    {
        if(NSPIN==2)
        {
            bfid.add_zeeman();
        }
    }
#endif

    timer::tick("potential","set_vrs_tddft");
    return;
} //end subroutine set_vrs_tddft

void potential::write_elecstat_pot(const string &fn, const string &fn_ave)
{
    TITLE("potential","write_elecstat_pot");
    timer::tick("potential","write_elecstat_pot");

    double *v_elecstat;
    v_elecstat = new double[pw.nrxx];
    ZEROS(v_elecstat, pw.nrxx);

    complex<double> *Porter = UFFT.porter;
    ZEROS( Porter, pw.nrxx );
    
    int nspin0 = 1;
    if(NSPIN==2) nspin0 = NSPIN;
    for(int is=0; is<nspin0; is++)
    {
        for(int ir=0; ir<pw.nrxx; ir++)
        {
            Porter[ir] += complex<double>( chr.rho[is][ir], 0.0 );
        }
    }

    //=============================
    //  bring rho (aux) to G space
    //=============================
    pw.FFT_chg.FFT3D(Porter, -1);

    //=======================================================
    // calculate hartree potential in G-space (NB: V(G=0)=0 )
    //=======================================================
    complex<double> *vh_g  = new complex<double>[pw.ngmc];
    ZEROS(vh_g, pw.ngmc);

    for(int ig = pw.gstart; ig<pw.ngmc; ig++)
    {
        const int j = pw.ig2fftc[ig];
        if(pw.gg[ig] >= 1.0e-12) //LiuXh 20180410
        {
            const double fac = e2 * FOUR_PI / (ucell.tpiba2 * pw.gg [ig]);
            vh_g[ig] = fac * Porter[j];
        }
    }

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
    //Add hartree potential and local pseudopot
    //==========================================
    for (int ir = 0;ir < pw.nrxx;ir++)
    {
        v_elecstat[ir] = Porter[ir].real() + this->vltot[ir];
    }

    //-------------------------------------------
    // output the electrostatic potential into a file.
    //-------------------------------------------
    ofstream ofs;
    ofstream ofs_ave;

    if(MY_RANK==0)
    {
        ofs.open( fn.c_str() );
        ofs_ave.open( fn_ave.c_str() );

        ofs << ucell.latName << endl;//1
        ofs << " " << ucell.lat0 * 0.529177 << endl;
        ofs << " " << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
        ofs << " " << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
        ofs << " " << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;

        ofs_ave << ucell.latName << endl;//1
        ofs_ave << " " << ucell.lat0 * 0.529177 << endl;
        ofs_ave << " " << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
        ofs_ave << " " << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
        ofs_ave << " " << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;

        for(int it=0; it<ucell.ntype; it++)
        {
            ofs << " " << ucell.atoms[it].label;
            ofs_ave << " " << ucell.atoms[it].label;
        }
        ofs << endl;
        ofs_ave << endl;
        for(int it=0; it<ucell.ntype; it++)
        {
            ofs << " " << ucell.atoms[it].na;
            ofs_ave << " " << ucell.atoms[it].na;
        }
        ofs << endl;
        ofs << "Direct" << endl;

        ofs_ave << endl;
        ofs_ave << "Direct" << endl;

        for(int it=0; it<ucell.ntype; it++)
        {
            for(int ia=0; ia<ucell.atoms[it].na; ia++)
            {
                ofs << " " << ucell.atoms[it].taud[ia].x
                    << " " << ucell.atoms[it].taud[ia].y
                    << " " << ucell.atoms[it].taud[ia].z << endl;

                ofs_ave << " " << ucell.atoms[it].taud[ia].x
                    << " " << ucell.atoms[it].taud[ia].y
                    << " " << ucell.atoms[it].taud[ia].z << endl;
            }
        }

        ofs << pw.ncx << " " << pw.ncy << " " << pw.ncz;
        ofs_ave << pw.ncx << " " << pw.ncy << " " << pw.ncz;

        int precision = 9;
        ofs << setprecision(precision);
        ofs << scientific; 
        ofs_ave << setprecision(precision);
        ofs_ave << scientific; 
        if(!ofs)
        {
            WARNING("potential::write_potential","Can't create VHartree File!");
        }
    }	

#ifndef __MPI
    int count=0;
    for(int k=0; k<pw.ncz; k++)
    {
        ofs << "\n" << k << " iz";
        double value = 0.0;
        double ave = 0.0;
        for(int j=0; j<pw.ncy; j++)
        {
            for(int i=0; i<pw.ncx; i++)
            {
                //if(count%8==0) ofs << "\n";
                if(count%5==0) ofs << "\n";
                value = v_elecstat[i*pw.ncy*pw.ncz + j*pw.ncz + k];
                ofs << " " << value;
                ave += value;
                ++count;
            }
        }
        //ofs << "\n" << ave/pw.ncx/pw.ncy << " average";
        if(k==0) ofs_ave << "iz" << "\taverage";
        ofs_ave << "\n" << k << "\t" << ave/pw.ncx/pw.ncy;
    }
#else
    MPI_Barrier(MPI_COMM_WORLD);
    // only do in the first pool.
    if(MY_POOL==0)
    {
        // num_z: how many planes on processor 'ip'
        int *num_z = new int[NPROC_IN_POOL];
        ZEROS(num_z, NPROC_IN_POOL);
        //for (int iz=0;iz<pw.ncz;iz++)
        //{
        //    int ip = iz % NPROC_IN_POOL;
        //    num_z[ip]++;
        //}
        for (int iz=0;iz<pw.nbz;iz++)
        {
            int ip = iz % NPROC_IN_POOL;
            num_z[ip] += pw.bz;
        }

        // start_z: start position of z in
        // processor ip.
        int *start_z = new int[NPROC_IN_POOL];
        ZEROS(start_z, NPROC_IN_POOL);
        for (int ip=1;ip<NPROC_IN_POOL;ip++)
        {
            start_z[ip] = start_z[ip-1]+num_z[ip-1];
        }

        // which_ip: found iz belongs to which ip.
        int *which_ip = new int[pw.ncz];
        ZEROS(which_ip, pw.ncz);
        for(int iz=0; iz<pw.ncz; iz++)
        {
            for(int ip=0; ip<NPROC_IN_POOL; ip++)
            {
                if(iz>=start_z[NPROC_IN_POOL-1])
                {
                    which_ip[iz] = NPROC_IN_POOL-1;
                    break;
                }
                else if(iz>=start_z[ip] && iz<start_z[ip+1])
                {
                    which_ip[iz] = ip;
                    break;
                }
            }
            //ofs_running << "\n iz=" << iz << " ip=" << which_ip[iz];
        }
        int count=0;
        int nxy = pw.ncx * pw.ncy;
        double* zpiece = new double[nxy];
        // save the rho one z by one z.
        for(int iz=0; iz<pw.ncz; iz++)
        {
            //ofs_running << "\n" << iz << " iz";
            // tag must be different for different iz.
            ZEROS(zpiece, nxy);
            int tag = iz;
            MPI_Status ierror;

            // case 1: the first part of rho in processor 0.
            if(which_ip[iz] == 0 && RANK_IN_POOL ==0)
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v_elecstat[ir*pw.nczp+iz-start_z[RANK_IN_POOL] ];
                    //ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*pw.nczp+iz=" << ir*pw.nczp+iz;
                }
            }
            // case 2: > first part rho: send the rho to
            // processor 0.
            else if(which_ip[iz] == RANK_IN_POOL )
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v_elecstat[ir*pw.nczp+iz-start_z[RANK_IN_POOL]];
                }
                MPI_Send(zpiece, nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
            }

            // case 2: > first part rho: processor 0 receive the rho
            // from other processors
            else if(RANK_IN_POOL==0)
            {
                MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
                //ofs_running << "\n Receieve First number = " << zpiece[0];
            }

            // write data
            if(MY_RANK==0)
            {
                //ofs << "\niz=" << iz;
                double ave = 0.0;
                /*
                for(int ir=0; ir<nxy; ir++)
                {
                    //if(count%8==0) ofs << "\n";
                    if(count%5==0) ofs << "\n";
                    //ofs << " " << zpiece[ir];
                    ofs << setw(17) << zpiece[ir];
                    ave += zpiece[ir];
                    ++count;
                }
                //ofs << "\n" << ave/nxy << " average"; 
                */
                for(int iy=0; iy<pw.ncy; iy++)
                {
                    for(int ix=0; ix<pw.ncx; ix++)
                    {
                        //if(count%8==0) ofs << "\n";
                        if(count%5==0) ofs << "\n";
                        //ofs << " " << zpiece[ir];
                        ofs << setw(17) << zpiece[ix*pw.ncy+iy];
                        ave += zpiece[ix*pw.ncy+iy];
                        ++count;
                    }
                }
                if(iz==0) ofs_ave << "\niz" << "\taverage";
                ofs_ave << "\n" << iz << "\t" << ave/pw.ncx/pw.ncy;
            }
        }
        delete[] num_z;
        delete[] start_z;
        delete[] which_ip;
        delete[] zpiece;
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(MY_RANK==0) ofs.close();

    delete[] v_elecstat;
    delete[] vh_g;

    timer::tick("potential","write_potential");
    return;
}
