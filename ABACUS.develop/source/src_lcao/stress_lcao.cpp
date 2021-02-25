#include "stress_lcao.h"
#include "../src_pw/global.h"
#include "../src_pw/xc_functional.h"
#include "../src_pw/xc_gga_pw.h"
#include "src_pw/myfunc.h"
// new
#include "src_pw/H_Hartree_pw.h"

double Stress_LCAO::stress_invalid_threshold_ev = 0.00;

Stress_LCAO::Stress_LCAO ()
{
    allocate_flag = false;
	output_acc = 1.0e-8;
}

Stress_LCAO::~Stress_LCAO ()
{
	this->destroy();
}

void Stress_LCAO::destroy (void)
{
    if (allocate_flag)
    {

        allocate_flag = false;
    }
}

void Stress_LCAO::allocate(void)
{
    TITLE("Stress_LCAO","init");

    if (allocate_flag)
    {
        this->destroy();
    }

    allocate_flag = true;
    return;
}

#include "../src_pw/efield.h"
#include "../src_pw/stress.h"
// be called in : Local_Orbital_Ions::force_stress
void Stress_LCAO::start_stress(double overlap[][3],double tvnl_dphi[][3],double vnl_dbeta[][3],double vl_dphi[][3], const matrix& stress_vdw)
{
    TITLE("Stress_LCAO","start_stress");
	timer::tick("Stress_LCAO","start_stress",'E');

    for(int i=0;i<3;i++)
	{
       for(int j=0;j<3;j++)
		{
          scs[i][j] = 0.0;
          soverlap[i][j] = overlap[i][j];
          stvnl_dphi[i][j] = tvnl_dphi[i][j];
          svnl_dbeta[i][j] = vnl_dbeta[i][j];
          svl_dphi[i][j] = vl_dphi[i][j];
          sigmacc[i][j] = 0.0;
          sigmadvl[i][j] = 0.0;
          sigmaewa[i][j] = 0.0;
          sigmaxc[i][j] = 0.0;
          sigmahar[i][j] = 0.0;
       }
    }
	//--------------------------------------------------------
	// local pseudopotential stress: 
	// use charge density; plane wave; local pseudopotential;
	//--------------------------------------------------------
    this->cal_stress_loc ();
 
	//--------------------------------------------------------
	//hartree term
	//--------------------------------------------------------
	this->cal_stress_har ();
    
	//--------------------------------------------------------
	// ewald stress: use plane wave only.
	//--------------------------------------------------------
    this->cal_stress_ew (); //remain problem


	//--------------------------------------------------------
	// stress due to core correlation. 
	//--------------------------------------------------------
	this->cal_stress_cc();

	//--------------------------------------------------------
	// stress due to self-consistent charge.
	//--------------------------------------------------------
	for(int i=0;i<3;i++)
	{
		// sigmaxc[i][i] = - (en.etxc-en.vtxc) / ucell.omega;
		sigmaxc[i][i] =  -(en.etxc) / ucell.omega;

		// sigmahar[i][i] = en.ehart /ucell.omega;

	}
	//Exchange-correlation for PBE
	cal_stress_gradcorr();

	//--------------------------------------------------------
	// need to move atom positions here.
	//--------------------------------------------------------
	if(GAMMA_ONLY_LOCAL)
	{
//    	this->stable_gamma();
	}
	else
	{
//		this->stable_k();
	}

/*	if(vdwd2.vdwD2)									//Peize Lin add 2014-04-04, update 2019-04-26
	{
		vdwd2.stress();
	}
*/	
/*	if(vdwd3.vdwD3)									//jiyy add 2019-05-18
	{
		vdwd3.stress();
	}
*/	
/*	matrix sefield;
	if(EFIELD)
	{
		sefield.create(3, 3);
		Efield::compute_stress(sefield);
	}*/

	for(int i=0; i<3; i++)
	{
		for (int j=0;j<3;j++)
		{
			scs[i][j] += soverlap[i][j]
				+ stvnl_dphi[i][j] 
				+ svnl_dbeta[i][j] 
				+ svl_dphi[i][j] 
				+ sigmadvl[i][j] // derivative of local potential stress (pw)
				+ sigmaewa[i][j] // ewald stress (pw)
				+ sigmacc[i][j] //nonlinear core correction stress (pw)
				+ sigmaxc[i][j]//exchange corretion stress 
				+ sigmahar[i][j];// hartree stress 

			if(vdwd2.vdwD2)			// Peize Lin update 2019-04-26
			{
				scs[i][j] += stress_vdw(i , j);
			}
			if(vdwd3.vdwD3)			// jiyy add 2019-05-18
			{
				scs[i][j] += stress_vdw(i , j);
			}   
			/*		if(EFIELD)
					{
					scs(iat, i) = scs(iat, i) + sefield(iat, i);
					}
			 */

		}

		//	if(OUT_LEVEL != "m") ofs_running << " correction stress for each atom along direction " 
		// << i+1 << " is " << sum/ucell.nat << endl;
	}

/*if(SYMMETRY)     
{
    double *pos;
    double d1,d2,d3;
    pos = new double[ucell.nat*3];
    ZEROS(pos, ucell.nat*3);
    int iat = 0;
    for(int it = 0;it < ucell.ntype;it++)
        {
                //Atom* atom = &ucell.atoms[it];
                for(int ia =0;ia< ucell.atoms[it].na;ia++)
                {
                        pos[3*iat  ] = ucell.atoms[it].taud[ia].x ;
                        pos[3*iat+1] = ucell.atoms[it].taud[ia].y ;
                        pos[3*iat+2] = ucell.atoms[it].taud[ia].z;
                        for(int k=0; k<3; ++k)
                        {
                                symm.check_translation( pos[iat*3+k], -floor(pos[iat*3+k]));
                                symm.check_boundary( pos[iat*3+k] );
                        }
                        iat++;

                }
        }

    for(int iat=0; iat<ucell.nat; iat++)
    {
                             Mathzone::Cartesian_to_Direct(fcs(iat,0),fcs(iat,1),fcs(iat,2),
                                        ucell.a1.x, ucell.a1.y, ucell.a1.z,
                                        ucell.a2.x, ucell.a2.y, ucell.a2.z,
                                        ucell.a3.x, ucell.a3.y, ucell.a3.z,
                                        d1,d2,d3);
                             fcs(iat,0) = d1;fcs(iat,1) = d2;fcs(iat,2) = d3;

    }
    symm.stress_symmetry(fcs , pos);
    for(int iat=0; iat<ucell.nat; iat++)
    {
                             Mathzone::Direct_to_Cartesian(fcs(iat,0),fcs(iat,1),fcs(iat,2),
                                        ucell.a1.x, ucell.a1.y, ucell.a1.z,
                                        ucell.a2.x, ucell.a2.y, ucell.a2.z,
                                        ucell.a3.x, ucell.a3.y, ucell.a3.z,
                                        d1,d2,d3);
                             fcs(iat,0) = d1;fcs(iat,1) = d2;fcs(iat,2) = d3;

    }
    //cout << "nrotk =" << symm.nrotk << endl;
   delete [] pos;
}
*/
	// test
    double svlocal[3][3];
    for (int i = 0; i<3; i++)
    {
		for(int j=0; j<3; j++)
		{
			svlocal[i][j] = 0.0;
			svlocal[i][j] = svl_dphi[i][j] + sigmadvl[i][j];
		}
	}

	// test
	double stvnl[3][3];
	for (int i = 0; i < 3; i++)
	{
		for(int j=0; j<3; j++)
		{
			stvnl[i][j] = 0.0;
			stvnl[i][j] = stvnl_dphi[i][j] + svnl_dbeta[i][j];
		}
	}

    if(Symmetry::symm_flag)
    {
		symm.stress_symmetry(scs);
    }//end symmetry

	// print Rydberg stress or not
	bool ry = false;
  
//  int TEST_STRESS = 1;
	if(TEST_STRESS)
	{
		ofs_running << "\n PARTS OF STRESS: " << endl;
		ofs_running << setiosflags(ios::showpos);
		ofs_running << setiosflags(ios::fixed) << setprecision(8) << endl;
		this->print_stress("OVERLAP    STRESS",soverlap,TEST_STRESS,ry);
		//test
		this->print_stress("T      STRESS",stvnl_dphi,TEST_STRESS,ry);
		this->print_stress("VNL      STRESS",svnl_dbeta,TEST_STRESS,ry);	

		this->print_stress("T_VNL      STRESS",stvnl,TEST_STRESS,ry);

		this->print_stress("VL_dPHI    STRESS",svl_dphi,TEST_STRESS,ry);
		this->print_stress("VL_dVL     STRESS",sigmadvl,TEST_STRESS,ry);
		this->print_stress("HAR     STRESS",sigmahar,TEST_STRESS,ry);

		this->print_stress("EWALD      STRESS",sigmaewa,TEST_STRESS,ry);
		this->print_stress("cc      STRESS",sigmacc,TEST_STRESS,ry);
		//		this->print_stress("NLCC       STRESS",sigmacc,TEST_STRESS,ry);
		this->print_stress("XC        STRESS",sigmaxc,TEST_STRESS,ry);
		this->print_stress("TOTAL        STRESS",scs,TEST_STRESS,ry);
	}


/*	if(EFIELD) 
	{
		STRESS::print("EFIELD     STRESS", sefield);
	}
*/
	ofs_running << setiosflags(ios::left);  
	
	this->printstress_total(ry);

	ofs_running << resetiosflags(ios::showpos);
	
	if(TEST_STRESS)
	{
		ofs_running << "\n STRESS INVALID TABLE." << endl;
	//	ofs_running << " " << setw(8) << "atom" << setw(5) << "x" << setw(5) << "y" << setw(5) << "z" << endl;
		for(int i=0;i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				if( abs( scs[i][j]*Ry_to_eV/0.529177 ) < Stress_LCAO::stress_invalid_threshold_ev)
				{
					scs[i][j] = 0.0;
					ofs_running << setw(5) << "1";
				}
				else
				{
					ofs_running << setw(5) << "0";
				}
			}
			ofs_running << endl;
		}
	}
	timer::tick("Stress_LCAO","start_stress",'E');
	return;
}

void Stress_LCAO::print_stress(const string &name, double f[][3], const bool screen, bool ry)const
{
	ofs_running << " --------------------------- " << name << " ----------------------------" << endl;
	
	
	double fac = 1.0;
	
	if(!ry)
	{
	 //	fac = Ry_to_eV / 0.529177;
	}

	cout << setprecision(5);
	cout << setiosflags(ios::showpos);

	if(screen)
	{
		cout << " ------------------- " << name << " --------------------" << endl;
	
	}

	for (int i=0;i<3;i++)
	{
		ofs_running << setw(15)<< " ";
		if( abs(f[i][0]) >output_acc) ofs_running << setw(15) << f[i][0]*fac;
		else ofs_running << setw(15) << "0";
		if( abs(f[i][1]) >output_acc) ofs_running << setw(15) << f[i][1]*fac;
		else ofs_running << setw(15) << "0";
		if( abs(f[i][2]) >output_acc) ofs_running << setw(15) << f[i][2]*fac;
		else ofs_running << setw(15) << "0";
		ofs_running << endl;

		if(screen)
		{
			if( abs(f[i][0]) >output_acc) cout << setw(15) << f[i][0]*fac;
			else cout << setw(15) << "0";
			if( abs(f[i][1]) >output_acc) cout << setw(15) << f[i][1]*fac;
			else cout << setw(15) << "0";
			if( abs(f[i][2]) >output_acc) cout << setw(15) << f[i][2]*fac;
			else cout << setw(15) << "0";
			cout << endl;
		}	
	}


	cout << resetiosflags(ios::showpos);

    return;
}

void Stress_LCAO::printstress_total(bool ry)
{
// zhengdy update 2016-10-08
	double unit_transform = 1;

	if(!ry)
	{
		unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * 1.0e-8;
	}
//	cout.setf(ios::fixed);


	//ofs_running << setiosflags(ios::right);
 	ofs_running << setprecision(6) << setiosflags(ios::showpos) << setiosflags(ios::fixed) << endl;
	NEW_PART("TOTAL-STRESS (KBAR)");//Ryd/(a.u.)^3
    cout << " ><><><><><><><><><><><><><><><><><><><><><><" << endl;
    cout << " TOTAL-STRESS (KBAR):" << endl;
    cout << " ><><><><><><><><><><><><><><><><><><><><><><" << endl;

//        if(INPUT.stress_set == 1)
//        int TEST_STRESS = 1;

 	if(TEST_STRESS) 
	{
		cout << setiosflags(ios::fixed) << setprecision(6);
		cout << setiosflags(ios::showpos);
		cout << " ------------------- TOTAL      STRESS --------------------" << endl;
    	cout << " " << setw(8) << "STRESS" << endl;
    	ofs_running << " " << setw(12) << "STRESS" << endl;
	}

    
    for (int i=0; i<3; i++)
    {

			//if(TEST_STRESS)
            cout << " " << setw(15) << scs[i][0]*unit_transform << setw(15)
                 << scs[i][1]*unit_transform << setw(15) << scs[i][2]*unit_transform << endl;

            ofs_running << " " << setw(15) << scs[i][0]*unit_transform << setw(15)
            << scs[i][1]*unit_transform << setw(15) << scs[i][2]*unit_transform << endl;

     
    }
	ofs_running << setiosflags(ios::left);
	cout << resetiosflags(ios::showpos);

    return;
}

void Stress_LCAO::cal_stress_loc(void)
{
    timer::tick("Stress_LCAO","cal_stress_loc",'F');

    double *dvloc;
    double evloc,fact=1.0;
    int ng,nt,l,m,is;

    dvloc = new double[pw.ngmc];

	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			sigmadvl[l][m]=0;
		}
	}

     complex<double> *Porter = UFFT.porter;

	 ZEROS( Porter, pw.nrxx );
	 for(int is=0; is<NSPIN; is++)
	 {
		 for (int ir=0; ir<pw.nrxx; ir++)
		 {
			 Porter[ir] += complex<double>(CHR.rho[is][ir], 0.0 );
		 }
	 }
	 pw.FFT_chg.FFT3D(Porter, -1);

//    if(INPUT.gamma_only==1) fact=2.0;
//    else fact=1.0;

    evloc=0.0;
    double g[3]={0,0,0};


    complex<double> *vg = new complex<double>[pw.ngmc];
    ZEROS( vg, pw.ngmc );
    for (int it=0; it<ucell.ntype; it++)
    {
        if (pw.gstart==1) evloc += ppcell.vloc(it, pw.ig2ngg[0]) * (pw.strucFac(it,0) * conj(Porter[pw.ig2fftc[0]])).real();
        for (int ig=pw.gstart; ig<pw.ngmc; ig++)
        {
            const int j = pw.ig2fftc[ig];
            evloc += ppcell.vloc(it, pw.ig2ngg[ig]) * (pw.strucFac(it,ig) * conj(Porter[j]) * fact).real();
        }
    }
    for( nt = 0;nt< ucell.ntype; nt++)
	{
        const Atom* atom = &ucell.atoms[nt];
        //mark by zhengdy for check
       // if ( ppcell.vloc == NULL ){
        if(0){
        //
        // special case: pseudopotential is coulomb 1/r potential
        //
        str.dvloc_coul (atom->zv, dvloc);
        //
        }
        else{
        //
        // normal case: dvloc contains dV_loc(G)/dG
        //
        str.dvloc_of_g ( atom->msh, atom->rab, atom->r,
          atom->vloc_at, atom->zv, dvloc);
        //
        }

        for( ng = 0;ng< pw.ngmc;ng++)
		{
            const int j = pw.ig2fftc[ng];
            g[0]=pw.gcar[ng].x;
            g[1]=pw.gcar[ng].y;
            g[2]=pw.gcar[ng].z;
            for (l = 0;l< 3;l++){
               for (m = 0; m<l+1;m++){
                  sigmadvl[l][m] = sigmadvl[l][m] +  ( conj( Porter[j] )
                                   * pw.strucFac (nt, ng) ).real() * 2.0 * dvloc [ pw.ig2ngg[ng] ]
                                   * ucell.tpiba2 * g[l] * g[m] * fact;
               }
            }
         }
    }

    for( l = 0;l< 3;l++)
	{
//     sigmadvl [l][l] = sigmadvl [l][l] +  evloc;//modified 2017/3/21
         for (m = 0; m<l+1; m++)
               sigmadvl [m][l] = sigmadvl [l][m];
    }
    for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigmadvl[l][m] );
		}
	}
    delete[] dvloc;
    delete[] vg;


    timer::tick("Stress_LCAO","cal_stress_loc");
    return;
}


#include "src_pw/H_Ewald_pw.h"
void Stress_LCAO::cal_stress_ew(void)
{
    timer::tick("Stress_lo","cal_stress_ew",'E');

    int i,j,l,m;
    double g[3];

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			sigmaewa[i][j]=0;
		}
	}

    double charge=0;
    for(int it=0; it < ucell.ntype; it++)
	{
		for(int i=0; i<ucell.atoms[it].na; i++)
		{
			charge = charge + ucell.atoms[it].zv;
		}
	}
    //choose alpha in order to have convergence in the sum over G
    //upperbound is a safe upper bound for the error ON THE ENERGY

    double alpha=2.9;
    double upperbound;
    do{
       alpha-=0.1;
       if(alpha==0.0)
          WARNING_QUIT("stres_ew", "optimal alpha not found");
       upperbound =e2 * pow(charge,2) * sqrt( 2 * alpha / (TWO_PI)) * erfc(sqrt(ucell.tpiba2 * pw.ggchg / 4.0 / alpha));
    }
    while(upperbound>1e-7);

    //G-space sum here
    //Determine if this processor contains G=0 and set the constant term 
    double sdewald;
    if(pw.gstart == 1){
       sdewald = (TWO_PI) * e2 / 4.0 / alpha * pow(charge/ucell.omega,2);
    }
    else {
       sdewald = 0.0;
    }

    //sdewald is the diagonal term 

    double fact=1.0;
//    if (INPUT.gamma_only) fact=2.0;
//    else fact=1.0;

    double g2,g2a;
    double arg;
    complex<double> rhostar;
    double sewald;
    for(int ng=pw.gstart;ng<pw.ngmc;ng++){
       g2 = pw.gg[ng]* ucell.tpiba2;
       g2a = g2 /4.0/alpha;
       rhostar=(0.0,0.0);
       for(int it=0; it < ucell.ntype; it++){
                for(int i=0; i<ucell.atoms[it].na; i++){
                       arg = (pw.gcar[ng].x * ucell.atoms[it].tau[i].x + pw.gcar[ng].y * ucell.atoms[it].tau[i].y + pw.gcar[ng].z * ucell.atoms[it].tau[i].z) * (TWO_PI);
                       rhostar = rhostar + complex<double>(ucell.atoms[it].zv * cos(arg),ucell.atoms[it].zv * sin(arg));
                }
       }
       rhostar /= ucell.omega;
       sewald = fact* (TWO_PI) * e2 * exp(-g2a) / g2 * pow(abs(rhostar),2);
       sdewald = sdewald - sewald;
       g[0]=pw.gcar[ng].x;
       g[1]=pw.gcar[ng].y;
       g[2]=pw.gcar[ng].z;
       for(l=0;l<3;l++){
          for(m=0;m<l+1;m++){
             sigmaewa[l][m]+=sewald * ucell.tpiba2 * 2.0 * g[l] * g[m] / g2 * (g2a+1);
          }
       }
    }

    for(l=0;l<3;l++){
       sigmaewa[l][l]+=sdewald;
    }
    //R-space sum here (only for the processor that contains G=0) 
    int mxr = 50;
    int *irr;
    Vector3<double> *r;
    double *r2;
    r  = new Vector3<double>[mxr];
    r2 = new double[mxr];
    irr = new int[mxr];
    double rr;
    Vector3<double> d_tau;
    double r0[3];
    int rmax , nrm=0;
    double fac;
    if(pw.gstart==1){
       rmax = 4.0/sqrt(alpha)/ucell.lat0;
       //with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2*10^-1)
       for(int it=0; it < ucell.ntype; it++){
         for(int i=0; i<ucell.atoms[it].na; i++){
           for(int jt=0; jt < ucell.ntype; jt++){
             for(int j=0; j<ucell.atoms[jt].na; j++){
               //calculate tau[na]-tau[nb]
               d_tau = ucell.atoms[it].tau[i] - ucell.atoms[jt].tau[j];
               //generates nearest-neighbors shells 
               H_Ewald_pw::rgen(d_tau, rmax, irr, ucell.latvec, ucell.G, r, r2, nrm);
               for(int nr=0 ; nr<nrm ; nr++){
                rr=sqrt(r2[nr]) * ucell.lat0;
                fac = -e2/2.0/ucell.omega*pow(ucell.lat0,2)*ucell.atoms[it].zv * ucell.atoms[it].zv / pow(rr,3) * (erfc(sqrt(alpha)*rr)+rr * sqrt(8 * alpha / (TWO_PI)) * exp(-alpha * pow(rr,2)));
                for(l=0; l<3; l++){
                   for(m=0; m<l+1; m++){
                      r0[0] = r[nr].x;
                      r0[1] = r[nr].y;
                      r0[2] = r[nr].z;
                      sigmaewa[l][m] += fac * r0[l] * r0[m];
                   }//end m
                }//end l
             }//end nr
          }//end j
        }//end jt
       }//end i
      }//end it
    }//end if
	for(l=0;l<3;l++){
		for(m=0;m<l+1;m++){
			sigmaewa[m][l]=sigmaewa[l][m];
		}
	}
	for(l=0;l<3;l++){
		for(m=0;m<3;m++){
			sigmaewa[l][m]=-sigmaewa[l][m];
			Parallel_Reduce::reduce_double_pool( sigmaewa[l][m] );
		}
	}

	delete[] r;
	delete[] r2;
	delete[] irr;
	// this->print(ofs_running, "ewald stress", stression);
	timer::tick("Force_lo","cal_stress_ew");

	return;
}

void Stress_LCAO::cal_stress_cc(void)
{
	timer::tick("Stress_LCAO","cal_stress_cc",'E');
        
	int nt,ng,l,m,ir;
	double fact=1.0;
	complex<double> sigmadiag;
	double* rhocg;
	double g[3];


	int judge=0;
	for(nt=0;nt<ucell.ntype;nt++){
		if(ucell.atoms[nt].nlcc) judge++;
	}
	if(judge==0) return;

	//recalculate the exchange-correlation potential
	matrix vxc(NSPIN, pw.nrxx);
	pot.v_xc(CHR.rho, en.etxc, en.vtxc, vxc);

	complex<double> * psic = new complex<double> [pw.nrxx];

	ZEROS(psic, pw.nrxx);

	if(NSPIN==1||NSPIN==4){
		for(ir=0;ir<pw.nrxx;ir++){
			// psic[ir] = vxc(0,ir);
			psic[ir] = complex<double>(vxc(0, ir),  0.0);
		}
	}
	else{
		for(ir=0;ir<pw.nrxx;ir++){
			psic[ir] = 0.5 * (vxc(0, ir) + vxc(1, ir));
		}
	}
	// to G space
	pw.FFT_chg.FFT3D(psic, -1);

	//psic cantains now Vxc(G)
	rhocg= new double [pw.nggm];
	ZEROS(rhocg, pw.nggm);

	sigmadiag=0.0;
	for(nt=0;nt<ucell.ntype;nt++){
		if(ucell.atoms[nt].nlcc){
			//drhoc();
			CHR.non_linear_core_correction(
				ppcell.numeric,
				ucell.atoms[nt].msh,
				ucell.atoms[nt].r,
				ucell.atoms[nt].rab,
				ucell.atoms[nt].rho_atc,
				rhocg);


			//diagonal term 
			if (pw.gstart==1) sigmadiag += conj(psic [pw.ig2fftc[0]] ) * pw.strucFac (nt, 0) * rhocg [pw.ig2ngg[0] ];
			for( ng = pw.gstart;ng< pw.ngmc;ng++){
				sigmadiag +=  conj(psic[pw.ig2fftc[ng]] ) *
					pw.strucFac (nt, ng) * rhocg [pw.ig2ngg[ng] ] * fact;
			}
			str.deriv_drhoc (
				ppcell.numeric,
				ucell.atoms[nt].msh,
				ucell.atoms[nt].r,
				ucell.atoms[nt].rab,
				ucell.atoms[nt].rho_atc,
				rhocg);
			// non diagonal term (g=0 contribution missing)
			for( ng = pw.gstart;ng< pw.ngmc;ng++){
				g[0] = pw.gcar[ng].x;
				g[1] = pw.gcar[ng].y;
				g[2] = pw.gcar[ng].z;
				for( l = 0;l< 3;l++){
					for (m = 0;m< 3;m++){
						const complex<double> t = conj(psic[pw.ig2fftc[ng]] )
							* pw.strucFac (nt, ng) * rhocg [pw.ig2ngg[ng] ] * ucell.tpiba *
							g [l] * g [m] / pw.gcar[ng].norm() * fact;
						sigmacc [l][ m] += t.real();
					}//end m
				}//end l
			}//end ng
		}//end if
	}//end nt
	for( l = 0;l< 3;l++){
		sigmacc [l][ l] += sigmadiag.real();
	}
	for( l = 0;l< 3;l++){
		for (m = 0;m< 3;m++){
			Parallel_Reduce::reduce_double_pool( sigmacc[l][m] );
		}
	}

	delete[] rhocg;
	delete[] psic;

	timer::tick("Stress_LCAO","cal_stress_cc",'E');
	return;
}

void Stress_LCAO::cal_stress_har(void)
{
	TITLE("Stress_LCAO","cal_stress_har");

	double shart,g2;
	const double eps=1e-8;
	int is,ig,l,m,nspin0;

	complex<double> *Porter = UFFT.porter;

	//  Hartree potential VH(r) from n(r)
	ZEROS( Porter, pw.nrxx );
	for(int is=0; is<NSPIN; is++)
	{
		for (int ir=0; ir<pw.nrxx; ir++)
		{
			Porter[ir] += complex<double>( CHR.rho[is][ir], 0.0 );
		}
	}
	//=============================
	//  bring rho (aux) to G space
	//=============================
	pw.FFT_chg.FFT3D(Porter, -1);

	complex<double> *psic = new complex<double> [pw.nrxx];
	double *psic0 = new double[pw.nrxx];
	ZEROS( psic0, pw.nrxx);
	for(int is=0; is<NSPIN; is++)
	{
		daxpy (pw.nrxx, 1.0, CHR.rho[is],1,psic0,2 );
		for (int ir=0; ir<pw.nrxx; ir++)
		{
			psic[ir] = complex<double>(psic0[ir], 0.0);
		}
	}

	pw.FFT_chg.FFT3D(psic, -1) ;

	double charge;
	if (pw.gstart == 1)
	{
		charge = ucell.omega * Porter[pw.ig2fftc[0]].real();
	}

	complex<double> *vh_g  = new complex<double>[pw.ngmc];
	ZEROS(vh_g, pw.ngmc);

	double g[3];

	// test
	// int i=pw.gstart;
	// cout<< "gstart " <<pw.gstart<<endl;

	double ehart=0;
	for (int ig = pw.gstart; ig<pw.ngmc; ig++)
	{
		const int j = pw.ig2fftc[ig];
		const double fac = e2 * FOUR_PI / (ucell.tpiba2 * pw.gg [ig]);

		ehart += ( conj( Porter[j] ) * Porter[j] ).real() * fac;
		//            vh_g[ig] = fac * Porter[j];
		shart= ( conj( Porter[j] ) * Porter[j] ).real()/(ucell.tpiba2 * pw.gg [ig]);
		g[0]=pw.gcar[ig].x;
		g[1]=pw.gcar[ig].y;
		g[2]=pw.gcar[ig].z;
		//test

		//   cout<<"g "<<g[0]<<" "<<g[1]<<" "<<g[2]<<endl;
		//   cout<<"gg "<<pw.gg[i]<<" "<<pw.gcar[i].norm2()<<endl;
		//   cout<<"Porter "<<Porter[j]<<endl;
		//   cout<<"tpiba2 "<<ucell.tpiba2<<endl;

		for(l=0;l<3;l++)
		{
			for(m=0;m<l+1;m++)
			{
				sigmahar[l][m]=sigmahar[l][m]+shart *2*g[l]*g[m]/pw.gg[ig];
			}
		}

	}
	Parallel_Reduce::reduce_double_pool( ehart );
	ehart *= 0.5 * ucell.omega;
	//cout<<"ehart "<<ehart<<" en.ehart "<< en.ehart<<endl;

	for(l=0;l<3;l++)
	{
		for(m=0;m<l+1;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigmahar[l][m] );
		}
	}

	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			sigmahar[l][m] *= 0.5 * e2 * FOUR_PI;
		}
	}

	for(l=0;l<3;l++)
	{
		sigmahar[l][l] += H_Hartree_pw::hartree_energy /ucell.omega;
	}

	for(l=0;l<3;l++)
	{
		for(m=0;m<l;m++)
		{
			sigmahar[m][l]=sigmahar[l][m];
		}
	}

	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			sigmahar[l][m]*=-1;
		}
	}

	delete[] vh_g;
	delete[] psic;
	delete[] psic0;
	return;
}

void Stress_LCAO::cal_stress_gradcorr() 
{
     
     if (xcf.igcx == 0  &&  xcf.igcc == 0)
     {
        return;
     } 
     double sigma_gradcorr[3][3];
     double* p= &sigma_gradcorr[0][0];
     for(int i=0;i<9;i++)
        *p++ = 0;

     bool igcc_is_lyp = false;
     if( xcf.igcc == 3 || xcf.igcc == 7)
     {
                igcc_is_lyp = true;
     }

     assert(NSPIN>0);
     const double fac = 1.0/ NSPIN;

     // doing FFT to get rho in G space: rhog1 
     CHR.set_rhog(CHR.rho[0], CHR.rhog[0]);
     if(NSPIN==2)//mohan fix bug 2012-05-28
     {
          CHR.set_rhog(CHR.rho[1], CHR.rhog[1]);
     }
     CHR.set_rhog(CHR.rho_core, CHR.rhog_core);
    
	double* rhotmp1;
	double* rhotmp2;
	complex<double>* rhogsum1;
	complex<double>* rhogsum2;
	Vector3<double>* gdr1;
	Vector3<double>* gdr2;
 
     rhotmp1 = new double[pw.nrxx];
     rhogsum1 = new complex<double>[pw.ngmc];
     ZEROS(rhotmp1, pw.nrxx);
     ZEROS(rhogsum1, pw.ngmc);
     for(int ir=0; ir<pw.nrxx; ir++) rhotmp1[ir] = CHR.rho[0][ir] + fac * CHR.rho_core[ir];
     for(int ig=0; ig<pw.ngmc; ig++) rhogsum1[ig] = CHR.rhog[0][ig] + fac * CHR.rhog_core[ig];
	gdr1 = new Vector3<double>[pw.nrxx];
	ZEROS(gdr1, pw.nrxx);

        GGA_PW::grad_rho( rhogsum1 , gdr1 );

        if(NSPIN==2)
        {
                rhotmp2 = new double[pw.nrxx];
                rhogsum2 = new complex<double>[pw.ngmc];
                ZEROS(rhotmp2, pw.nrxx);
                ZEROS(rhogsum2, pw.ngmc);
                for(int ir=0; ir<pw.nrxx; ir++) rhotmp2[ir] = CHR.rho[1][ir] + fac * CHR.rho_core[ir];
                for(int ig=0; ig<pw.ngmc; ig++) rhogsum2[ig] = CHR.rhog[1][ig] + fac * CHR.rhog_core[ig];

                        gdr2 = new Vector3<double>[pw.nrxx];
                        ZEROS(gdr2, pw.nrxx);

                GGA_PW::grad_rho( rhogsum2 , gdr2 );
        }
        
        const double epsr = 1.0e-6;
        const double epsg = 1.0e-10;

        double grho2a = 0.0;
        double grho2b = 0.0;
        double sx = 0.0;
        double sc = 0.0;
        double v1x = 0.0;
        double v2x = 0.0;
        double v1c = 0.0;
        double v2c = 0.0;
        double vtxcgc = 0.0;
        double etxcgc = 0.0;

        if(NSPIN==1||NSPIN==4)
        {
                double segno;
                for(int ir=0; ir<pw.nrxx; ir++)
                {
                        const double arho = std::abs( rhotmp1[ir] );
			if(arho > epsr)
			{
				grho2a = gdr1[ir].norm2();
				if( grho2a > epsg )
				{
					if( rhotmp1[ir] >= 0.0 ) segno = 1.0;
					if( rhotmp1[ir] < 0.0 ) segno = -1.0;

					XC_Functional::gcxc( arho, grho2a, sx, sc, v1x, v2x, v1c, v2c);
					double tt[3];
					tt[0] = gdr1[ir].x;
					tt[1] = gdr1[ir].y;
					tt[2] = gdr1[ir].z;
					for(int l = 0;l< 3;l++){
						for(int m = 0;m< l+1;m++){
							sigma_gradcorr[l][m] += tt[l] * tt[m] * e2 * (v2x + v2c);
						}
					}
				}
			}
		} 
	}
	else if(NSPIN==2)
	{
		double v1cup = 0.0;
		double v1cdw = 0.0;
		double v2cup = 0.0;
		double v2cdw = 0.0;
		double v1xup = 0.0;
		double v1xdw = 0.0;
		double v2xup = 0.0;
		double v2xdw = 0.0;
		double v2cud = 0.0;
		double v2c = 0.0;
		for(int ir=0; ir<pw.nrxx; ir++)
		{
			double rh = rhotmp1[ir] + rhotmp2[ir];
			grho2a = gdr1[ir].norm2();;
			grho2b = gdr2[ir].norm2();;
			//XC_Functional::gcx_spin();
			gcx_spin(rhotmp1[ir], rhotmp2[ir], grho2a, grho2b,
				sx, v1xup, v1xdw, v2xup, v2xdw);

			if(rh > epsr)
			{
				if(igcc_is_lyp)
				{
					WARNING_QUIT("stress","igcc_is_lyp is not available now.");
				}
				else
				{
					double zeta = ( rhotmp1[ir] - rhotmp2[ir] ) / rh;
					double grh2 = (gdr1[ir]+gdr2[ir]).norm2();
					//XC_Functional::gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
					gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
					v2cup = v2c;
					v2cdw = v2c;
					v2cud = v2c;
				}
			}
			else
			{
				sc = 0.0;
				v1cup = 0.0;
				v1cdw = 0.0;
				v2c = 0.0;
				v2cup = 0.0;
				v2cdw = 0.0;
				v2cud = 0.0;
			}
			double tt1[3],tt2[3];
			{
				tt1[0] = gdr1[ir].x;
				tt1[1] = gdr1[ir].y;
				tt1[2] = gdr1[ir].z;
				tt2[0] = gdr2[ir].x;
				tt2[1] = gdr2[ir].y;
				tt2[2] = gdr2[ir].z;
			}
			for(int l = 0;l< 3;l++){
			    for(int m = 0;m< l+1;m++){
				//    exchange
				sigma_gradcorr [l][m] += tt1[l] * tt1[m] * e2 * v2xup + 
							tt2[l] * tt2[m] * e2 * v2xdw;
				//    correlation
				sigma_gradcorr [l][m] += ( tt1[l] * tt1[m] * v2cup + 
							tt2[l] * tt2[m] * v2cdw + 
							(tt1[l] * tt2[m] +
							tt2[l] * tt1[m] ) * v2cud ) * e2;
			    }
			}
		}
	}

	for(int l = 0;l< 3;l++){
		for(int m = 0;m< l;m++){
			sigma_gradcorr[m][l] = sigma_gradcorr[l][m];
		}
	}
	for(int l = 0;l<3;l++){
		for(int m = 0;m<3;m++){
			Parallel_Reduce::reduce_double_pool( sigma_gradcorr[l][m] );
		}
	}
	p= &sigma_gradcorr[0][0];
	double* p1 = &sigmaxc[0][0];
	for(int i=0;i<9;i++){
		*p /= pw.ncxyz ;
		*p1++ += *p++;  
	}

	delete[] rhotmp1;
	delete[] rhogsum1;
	delete[] gdr1;
	if(NSPIN==2)
	{
		delete[] rhotmp2;
		delete[] rhogsum2;
		delete[] gdr2;
	}
	return;
}
