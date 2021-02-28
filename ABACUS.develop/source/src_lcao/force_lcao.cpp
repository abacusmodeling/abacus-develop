#include "force_lcao.h"
#include "src_pw/global.h"
#include "src_pw/potential_libxc.h"
#include "dftu.h"  //Quxin add for DFT+U on 20201029
// new
#include "src_pw/H_XC_pw.h"

double Force_LCAO::force_invalid_threshold_ev = 0.00;

Force_LCAO::Force_LCAO ()
{
    allocate_flag = false;
	output_acc = 1.0e-8;
}

Force_LCAO::~Force_LCAO ()
{
	this->destroy();
}

void Force_LCAO::destroy (void)
{
    if (allocate_flag)
    {
        for (int iat = 0; iat < ucell.nat; iat++)
        {
			//-------------------------------
			// use energy density matrix
			//-------------------------------
            delete [] foverlap[iat];
			//-------------------------------
			// use density matrix
			//-------------------------------
            delete [] ftvnl_dphi[iat];
            delete [] fvnl_dbeta[iat];
			//-------------------------------
			// use grid integration
			//-------------------------------
            delete [] fvl_dphi[iat];
			//-------------------------------
			// use plane wave
			//-------------------------------
            delete [] fvl_dvl[iat];
            delete [] fewalds[iat];
			delete [] fcc[iat];
            delete [] fscc[iat];
        }

        delete [] foverlap;
        delete [] ftvnl_dphi;
        delete [] fvnl_dbeta;
		delete [] fvl_dphi;        
		delete [] fvl_dvl;
        delete [] fewalds;
		delete [] fcc;
        delete [] fscc;

        allocate_flag = false;
    }
}

void Force_LCAO::allocate(void)
{
    TITLE("Force_LCAO","init");

    // reduce memory occupy by vlocal
    delete[] ParaO.sender_local_index;
    delete[] ParaO.sender_size_process;
    delete[] ParaO.sender_displacement_process;
    delete[] ParaO.receiver_global_index;
    delete[] ParaO.receiver_size_process;
    delete[] ParaO.receiver_displacement_process;

    ParaO.sender_local_index = new int[1];
    ParaO.sender_size_process = new int[1];
    ParaO.sender_displacement_process = new int[1];
    ParaO.receiver_global_index = new int[1];
    ParaO.receiver_size_process = new int[1];
    ParaO.receiver_displacement_process = new int[1];

    if (allocate_flag)
    {
        this->destroy();
    }

	const int nat = ucell.nat;
	
    this->fcs.create (nat, 3);

	// part of total force
    foverlap = new double*[nat];
    ftvnl_dphi = new double*[nat];
    fvnl_dbeta = new double*[nat];
    fvl_dphi = new double*[nat];
    fvl_dvl = new double*[nat];
    fewalds = new double*[nat];
	fcc = new double*[nat];
    fscc = new double*[nat];

    for (int iat = 0; iat < nat; iat++)
    {
        foverlap[iat] = new double[3];
        ftvnl_dphi[iat] = new double[3];
        fvnl_dbeta[iat] = new double[3];
        fvl_dphi[iat] = new double[3];
		fvl_dvl[iat] = new double[3];
        fewalds[iat] = new double[3];
		fcc[iat] = new double[3];
        fscc[iat] = new double[3];

        ZEROS (foverlap[iat], 3);
        ZEROS (ftvnl_dphi[iat], 3);
        ZEROS (fvnl_dbeta[iat], 3);
        ZEROS (fvl_dphi[iat],3);
		ZEROS (fvl_dvl[iat],3);
        ZEROS (fewalds[iat], 3);
		ZEROS (fcc[iat], 3);
        ZEROS (fscc[iat],3);
    }

    allocate_flag = true;

    Memory::record("Force_LCAO","Paremeters",nat*3*8,"double");
    return;
}

#include "../src_pw/efield.h"
#include "../src_pw/forces.h"
// be called in : Local_Orbital_Ions::force_stress
void Force_LCAO::start_force(void)
{
    TITLE("Force_LCAO","start_force");
	timer::tick("Force_LCAO","start_force",'E');


	//--------------------------------------------------------
	// local pseudopotential force: 
	// use charge density; plane wave; local pseudopotential;
	//--------------------------------------------------------
    this->cal_force_loc ();

	//--------------------------------------------------------
	// ewald force: use plane wave only.
	//--------------------------------------------------------
    this->cal_force_ew (); //remain problem


	//--------------------------------------------------------
	// force due to core correlation. 
	//--------------------------------------------------------
	this->cal_force_cc();

	//--------------------------------------------------------
	// force due to self-consistent charge.
	//--------------------------------------------------------
    this->cal_force_scc();

	//--------------------------------------------------------
	// need to move atom positions here.
	//--------------------------------------------------------
	if(GAMMA_ONLY_LOCAL)
	{
    	this->ftable_gamma();
	}
	else
	{
		this->ftable_k();
	}

	// clear the data.
	this->fcs.zero_out();

	// zhengdy added in 2018-10-29
	stress_vdw.create(3,3);
	// Peize Lin add 2014-04-04, update 2019-04-26
	if(vdwd2.vdwD2)
	{
		vdwd2.force(stress_vdw, STRESS);
	}
	// jiyy add 2019-05-18
	else if(vdwd3.vdwD3)
	{
		vdwd3.force(stress_vdw, STRESS);
	}																			 
		
	matrix fefield;
	if(EFIELD)
	{
		fefield.create(ucell.nat, 3);
		Efield::compute_force(fefield);
	}

	for(int i=0; i<3; i++)
    {
		double sum = 0.0;

    	for (int iat = 0; iat < ucell.nat; iat++)
		{
        	fcs(iat, i) += foverlap[iat][i]
			+ ftvnl_dphi[iat][i] 
			+ fvnl_dbeta[iat][i] 
			+ fvl_dphi[iat][i] 
			+ fvl_dvl[iat][i] // derivative of local potential force (pw)
			+ fewalds[iat][i] // ewald force (pw)
			+ fcc[iat][i] //nonlinear core correction force (pw)
			+ fscc[iat][i];//self consistent corretion force (pw)

			// Force contribution from DFT+U, Quxin add on 20201029
            if(INPUT.dft_plus_u) 
			{
				fcs(iat, i) += dftu.force_dftu.at(iat).at(i);
			}
	
			// Peize Lin add 2014-04-04, update 2019-04-261
			if(vdwd2.vdwD2)
			{
				switch(i)
				{
					case 0:	fcs(iat,i) += vdwd2.force_result[iat].x;	break;
					case 1:	fcs(iat,i) += vdwd2.force_result[iat].y;	break;
					case 2:	fcs(iat,i) += vdwd2.force_result[iat].z;	break;
				}
				
			}
			// jiyy add 2019-05-18
	        if(vdwd3.vdwD3)
			{
				switch(i)
				{
					case 0:	fcs(iat,i) += vdwd3.force_result[iat][0];	break;
					case 1:	fcs(iat,i) += vdwd3.force_result[iat][1];	break;
					case 2:	fcs(iat,i) += vdwd3.force_result[iat][2];	break;
				}
				
			}
			
			if(EFIELD)
			{
				fcs(iat, i) = fcs(iat, i) + fefield(iat, i);
			}

			sum += fcs(iat, i);
		}

		for(int iat=0; iat<ucell.nat; ++iat)
		{
			fcs(iat, i) -= sum/ucell.nat;
		}

 		//xiaohui add "OUT_LEVEL", 2015-09-16
		if(OUT_LEVEL != "m") 
		{
			ofs_running << " correction force for each atom along direction " 
				<< i+1 << " is " << sum/ucell.nat << endl;
		}
    }
	

	// pengfei 2016-12-20
	if(Symmetry::symm_flag)
	{
		double *pos;
		double d1,d2,d3;
		pos = new double[ucell.nat*3];
		ZEROS(pos, ucell.nat*3);
		int iat = 0;
		for(int it = 0;it < ucell.ntype;it++)
		{
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
		symm.force_symmetry(fcs , pos);
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
		delete[] pos;
	}

	// test
	double** fvlocal = new double*[ucell.nat];
    for (int iat = 0; iat < ucell.nat; iat++)
    {
		fvlocal[iat] = new double[3];
		ZEROS(fvlocal[iat], 3);
		for(int i=0; i<3; i++)
		{
        	fvlocal[iat][i] = fvl_dphi[iat][i] + fvl_dvl[iat][i];
		}
    }

	// test
	double** ftvnl = new double*[ucell.nat];
    for (int iat = 0; iat < ucell.nat; iat++)
    {
		ftvnl[iat] = new double[3];
		ZEROS(ftvnl[iat], 3);
		for(int i=0; i<3; i++)
		{
        	ftvnl[iat][i] = ftvnl_dphi[iat][i] + fvnl_dbeta[iat][i];
		}
    }

	// print Rydberg force or not
	bool ry = false;
	if(TEST_FORCE)
	{
		ofs_running << "\n PARTS OF FORCE: " << endl;
		ofs_running << setiosflags(ios::showpos);
		ofs_running << setiosflags(ios::fixed) << setprecision(8) << endl;
		this->print_force("OVERLAP    FORCE",foverlap,TEST_FORCE,ry);
		//  this->print_force("TVNL_DPHI  force",ftvnl_dphi,TEST_FORCE);
		//  this->print_force("VNL_DBETA  force",fvnl_dbeta,TEST_FORCE);
		this->print_force("T_VNL      FORCE",ftvnl,TEST_FORCE,ry);

		this->print_force("VL_dPHI    FORCE",fvl_dphi,TEST_FORCE,ry);
		this->print_force("VL_dVL     FORCE",fvl_dvl,TEST_FORCE,ry);
		// 	this->print_force("VLOCAL     FORCE",fvlocal,TEST_FORCE);

		this->print_force("EWALD      FORCE",fewalds,TEST_FORCE,ry);
		this->print_force("NLCC       FORCE",fcc,TEST_FORCE,ry);
		this->print_force("SCC        FORCE",fscc,TEST_FORCE,ry);
	}

	// mohan fix bug 2012-03-22
	for(int iat=0; iat<ucell.nat; ++iat)
	{
		delete[] ftvnl[iat];
		delete[] fvlocal[iat];
	}
	delete[] ftvnl;
	delete[] fvlocal; 


	if(EFIELD) 
	{
		Forces::print("EFIELD     FORCE", fefield);
	}

	ofs_running << setiosflags(ios::left);  
	
	this->printforce_total(ry);

	ofs_running << resetiosflags(ios::showpos);
	
	if(TEST_FORCE)
	{
		ofs_running << "\n FORCE INVALID TABLE." << endl;
		ofs_running << " " << setw(8) << "atom" << setw(5) << "x" << setw(5) << "y" << setw(5) << "z" << endl;
		for(int iat=0; iat<ucell.nat; iat++)
		{
			ofs_running << " " << setw(8) << iat;
			for(int i=0; i<3; i++)
			{
				if( abs( fcs(iat,i)*Ry_to_eV/0.529177 ) < Force_LCAO::force_invalid_threshold_ev)
				{
					fcs(iat,i) = 0.0;
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
	timer::tick("Force_LCAO","start_force",'E');
	return;
}

void Force_LCAO::print_force(const string &name, double** f, const bool screen, bool ry)const
{
	ofs_running << " --------------------------- " << name << " ----------------------------" << endl;
	ofs_running << " " << setw(8) << "atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
	
	double fac = 1.0;
	
	if(!ry)
	{
	 	fac = Ry_to_eV / 0.529177;
	}

	cout << setprecision(5);
	cout << setiosflags(ios::showpos);

	if(screen)
	{
		cout << " ------------------- " << name << " --------------------" << endl;
		cout << " " << setw(8) << "atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
	}

    int iat = 0;
    for (int it = 0;it < ucell.ntype;it++)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;ia++)
        {
			stringstream ss;
			ss << ucell.atoms[it].label << ia+1;

			ofs_running << " " << setw(8) << ss.str();
			if( abs(f[iat][0]) >output_acc) ofs_running << setw(15) << f[iat][0]*fac;
			else ofs_running << setw(15) << "0";
			if( abs(f[iat][1]) >output_acc) ofs_running << setw(15) << f[iat][1]*fac;
			else ofs_running << setw(15) << "0";
			if( abs(f[iat][2]) >output_acc) ofs_running << setw(15) << f[iat][2]*fac;
			else ofs_running << setw(15) << "0";
			ofs_running << endl;

			if(screen)
			{
				cout << " " << setw(8) << ss.str();
				if( abs(f[iat][0]) >output_acc) cout << setw(15) << f[iat][0]*fac;
				else cout << setw(15) << "0";
				if( abs(f[iat][1]) >output_acc) cout << setw(15) << f[iat][1]*fac;
				else cout << setw(15) << "0";
				if( abs(f[iat][2]) >output_acc) cout << setw(15) << f[iat][2]*fac;
				else cout << setw(15) << "0";
				cout << endl;
			}	
				
            iat++;
        }
    }


	cout << resetiosflags(ios::showpos);

    return;
}

void Force_LCAO::printforce_total (bool ry)
{
// mohan update 2011-03-15
	double unit_transform = 1;

	if(!ry)
	{
		unit_transform = Ry_to_eV / 0.529177;
	}
//	cout.setf(ios::fixed);

    int iat=0;

	//ofs_running << setiosflags(ios::right);
 	ofs_running << setprecision(6) << setiosflags(ios::showpos) << setiosflags(ios::fixed) << endl;
	NEW_PART("TOTAL-FORCE (eV/Angstrom)");

        if(INPUT.force_set == 1)
        {
           ofstream ofs("FORCE.dat");
           if(!ofs)
           {
              cout << "open FORCE.dat error !" <<endl;
           }

           for(int iat=0; iat<ucell.nat; iat++)
           {
               ofs << "   " << fcs(iat,0)*Ry_to_eV / 0.529177 << "   " << fcs(iat,1)*Ry_to_eV / 0.529177 << "   " << fcs(iat,2)*Ry_to_eV / 0.529177 << endl;
           }

           ofs.close();
        }

 	if(TEST_FORCE) 
	{
		cout << setiosflags(ios::fixed) << setprecision(6);
		cout << setiosflags(ios::showpos);
		cout << " ------------------- TOTAL      FORCE --------------------" << endl;
    	cout << " " << setw(8) << "Atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
    	ofs_running << " " << setw(12) << "Atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
	}

    iat=0;
    for (int it=0; it<ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            stringstream ss;
            ss << ucell.atoms[it].label << ia+1;

			if(TEST_FORCE)
            cout << " " << setw(8) << ss.str() << setw(15) << fcs(iat,0)*unit_transform << setw(15)
                 << fcs(iat,1)*unit_transform << setw(15) << fcs(iat,2)*unit_transform << endl;

            ofs_running << " " << setw(12) << ss.str() << setw(15) << fcs(iat,0)*unit_transform << setw(15)
            << fcs(iat,1)*unit_transform << setw(15) << fcs(iat,2)*unit_transform << endl;

            ++iat;
        }
    }
	ofs_running << setiosflags(ios::left);
	cout << resetiosflags(ios::showpos);

    return;
}

void Force_LCAO::cal_force_loc(void)
{
    timer::tick("Force_LCAO","cal_force_loc",'F');

    complex<double> *aux = new complex<double>[pw.nrxx];
    ZEROS(aux, pw.nrxx);

    // now, in all pools , the charge are the same,
    // so, the force calculated by each pool is equal.

    for (int is=0; is<NSPIN; is++)
    {
        for (int ir=0; ir<pw.nrxx; ir++)
        {
            aux[ir] += complex<double>( CHR.rho[is][ir], 0.0 );
        }
    }

    pw.FFT_chg.FFT3D(aux, -1);

    int gstart_here = pw.gstart;
    if (pw.ggs[0] != 0) gstart_here = 0;

    int iat = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
			fvl_dvl[iat][0]=0.0;
			fvl_dvl[iat][1]=0.0;
			fvl_dvl[iat][2]=0.0;

            for (int ig = gstart_here; ig < pw.ngmc; ig++)
            {
                const double phase = TWO_PI * (pw.gcar[ig] * ucell.atoms[it].tau[ia]);
                const double factor = ppcell.vloc(it, pw.ig2ngg[ig]) *
                                      ( cos(phase) * aux[ pw.ig2fftc[ig] ].imag()
                                        + sin(phase) * aux[ pw.ig2fftc[ig] ].real());

                fvl_dvl[iat][0] += pw.gcar[ig].x * factor;
                fvl_dvl[iat][1] += pw.gcar[ig].y * factor;
                fvl_dvl[iat][2] += pw.gcar[ig].z * factor;
            }
            for (int ipol = 0; ipol < 3; ipol++)
            {
                this->fvl_dvl[iat][ipol] *= (ucell.tpiba * ucell.omega);
            }
            ++iat;

        }
    }

    //this->print(ofs_running, "local forces", forcelc);
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        Parallel_Reduce::reduce_double_pool(this->fvl_dvl[iat], 3);
    }

    delete[] aux;

    timer::tick("Force_LCAO","cal_force_loc");
    return;
}


#include "src_pw/H_Ewald_pw.h"
void Force_LCAO::cal_force_ew(void)
{
    timer::tick("Force_lo","cal_force_ew",'E');

    double fact = 2.0;
    complex<double> *aux = new complex<double> [pw.ngmc];
    ZEROS(aux, pw.ngmc);

    int gstart = pw.gstart;

    if (pw.ggs[0] != 0) gstart = 0;

    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ig = gstart; ig < pw.ngmc; ig++)
        {
            aux[ig] += static_cast<double>(ucell.atoms[it].zv) * conj(pw.strucFac(it, ig));
        }
    }

    for (int ig = gstart; ig < pw.ngmc; ig++)
    {
        aux[ig] *= exp(-1.0 * pw.gg[ig] * ucell.tpiba2 / H_Ewald_pw::alpha / 4.0) / (pw.gg[ig] * ucell.tpiba2);
    }

    int iat = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
			ZEROS(fewalds[iat],3);
            for (int ig = gstart; ig < pw.ngmc; ig++)
            {
                const double arg = TWO_PI * (pw.gcar[ig] * ucell.atoms[it].tau[ia]);
                double sumnb =  -cos(arg) * aux[ig].imag() + sin(arg) * aux[ig].real();
                fewalds[iat][0] += pw.gcar[ig].x * sumnb;
                fewalds[iat][1] += pw.gcar[ig].y * sumnb;
                fewalds[iat][2] += pw.gcar[ig].z * sumnb;
            }
            for (int ipol = 0; ipol < 3; ipol++)
            {
                fewalds[iat][ipol] *= ucell.atoms[it].zv * e2 * ucell.tpiba * TWO_PI / ucell.omega * fact;
            }
            iat++;
        }
    }
    delete [] aux;

    for (int iat = 0; iat < ucell.nat; iat++)
    {
        Parallel_Reduce::reduce_double_pool(this->fewalds[iat], 3);
    }

    // mohan fix bug 2010-07-15
    // means that the processor contains G=0 term.
    if (gstart == 1)
    {
        double rmax = 5.0 / (sqrt(H_Ewald_pw::alpha) * ucell.lat0);
        int nrm = 0;

        //output of rgen: the number of vectors in the sphere
        const int mxr = 50;
        // the maximum number of R vectors included in r
        Vector3<double> *r  = new Vector3<double>[mxr];
        double *r2 = new double[mxr];
        int *irr = new int[mxr];
		ZEROS(r2, mxr);
		ZEROS(irr, mxr);
        // the square modulus of R_j-tau_s-tau_s'

        int iat1 = 0;
        for (int it1 = 0; it1 < ucell.ntype; it1++)
        {
            for (int ia1 = 0; ia1 < ucell.atoms[it1].na; ia1++)
            {
        		int iat2 = 0; // mohan fix serious bug 2011-06-07
                for (int it2 = 0; it2 < ucell.ntype; it2++)
                {
                    for (int ia2 = 0; ia2 < ucell.atoms[it2].na; ia2++)
                    {
                        if (iat1 != iat2)
                        {
                            Vector3<double> d_tau = ucell.atoms[it1].tau[ia1] - ucell.atoms[it2].tau[ia2];
                            H_Ewald_pw::rgen(d_tau, rmax, irr, ucell.latvec, ucell.G, r, r2, nrm);
                            // return r2[n], r(n, ipol)
//							cout << "nrm = " << nrm << endl;

                            for (int n = 0; n < nrm; n++)
                            {
                                assert (H_Ewald_pw::alpha >= 0.0);

                                const double rr = sqrt(r2[n]) * ucell.lat0;
                                double factor = ucell.atoms[it1].zv * ucell.atoms[it2].zv * e2 / (rr * rr)
                                                * (erfc(sqrt(H_Ewald_pw::alpha) * rr) / rr
                                                   + sqrt(8.0 * H_Ewald_pw::alpha / TWO_PI) 
								* exp(-1.0 * H_Ewald_pw::alpha * rr * rr)) * ucell.lat0;

                                //remian problem
                                //fix bug iat -> iat1
                                fewalds[iat1][0] -= factor * r[n].x;
                                fewalds[iat1][1] -= factor * r[n].y;
                                fewalds[iat1][2] -= factor * r[n].z;
                            }
                        }
                        iat2++;
                    }//ia2
                }//atom b
                iat1++;
            }//ia1
        }//atom a

		delete[] r;
		delete[] r2;
		delete[] irr;
    }//gstart

#ifdef __MPI
    // (1) only used force to do BFGS in processor 0.
    // or (2) ewalds bcast from processor 0, which is
    // happen to be gstart=1.
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        Parallel_Common::bcast_double(this->fewalds[iat], 3);
    }
#endif

    // this->print(ofs_running, "ewald forces", forceion);
    timer::tick("Force_lo","cal_force_ew");

    return;
}

void Force_LCAO::cal_force_cc(void)
{
	timer::tick("Force_LCAO","cal_force_cc",'E');
	// recalculate the exchange-correlation potential.
    matrix vxc(NSPIN, pw.nrxx);
	#ifdef USE_LIBXC
	Potential_Libxc::v_xc(CHR.rho, en.etxc, en.vtxc, vxc);
	#else
    H_XC_pw::v_xc(pw.nrxx, pw.ncxyz, ucell.omega, CHR.rho, CHR.rho_core, vxc);
	#endif

    complex<double> * psiv = new complex<double> [pw.nrxx];
    ZEROS(psiv, pw.nrxx);
    if (NSPIN == 1 || NSPIN == 4)
    {
        for (int ir = 0;ir < pw.nrxx;ir++)
        {
            psiv[ir] = complex<double>(vxc(0, ir),  0.0);
        }
    }
    else
    {
        for (int ir = 0;ir < pw.nrxx;ir++)
        {
            psiv[ir] = 0.5 * (vxc(0 ,ir) + vxc(1, ir));
        }
    }

	// to G space
    pw.FFT_chg.FFT3D(psiv, -1);

    //psiv contains now Vxc(G)
    double * rhocg = new double [pw.nggm];
    ZEROS(rhocg, pw.nggm);

    for (int T1 = 0;T1 < ucell.ntype;T1++)
    {
        if (ucell.atoms[T1].nlcc)
        {
            //call drhoc
            CHR.non_linear_core_correction(
                ppcell.numeric,
                ucell.atoms[T1].msh,
                ucell.atoms[T1].r,
                ucell.atoms[T1].rab,
                ucell.atoms[T1].rho_atc,
                rhocg);
            int iat = 0;


			complex<double> ipol0, ipol1, ipol2;
            for (int T2 = 0;T2 < ucell.ntype;T2++)
            {
            	for (int I2 = 0;I2 < ucell.atoms[T2].na;I2++)
                {
                    if (T2 == T1)
                    {
                        for (int ig = pw.gstart; ig < pw.ngmc; ig++)
                        {
                            const double arg = TWO_PI * (pw.gcar[ig].x * ucell.atoms[T2].tau[I2].x
                                                      + pw.gcar[ig].y * ucell.atoms[T2].tau[I2].y
                                                      + pw.gcar[ig].z * ucell.atoms[T2].tau[I2].z);

                            ipol0 = ucell.tpiba * ucell.omega * rhocg[pw.ig2ngg[ig]]
                                                    * pw.gcar[ig].x * conj(psiv[pw.ig2fftc[ig]]) 
													* complex<double>(sin(arg), cos(arg))  ;
                            this->fcc[iat][0] +=  ipol0.real();

                            ipol1 = ucell.tpiba * ucell.omega * rhocg[pw.ig2ngg[ig]]
                                                    * pw.gcar[ig].y * conj(psiv[pw.ig2fftc[ig]]) 
													* complex<double>(sin(arg), cos(arg)) ;
                            this->fcc[iat][1] += ipol1.real();

                            ipol2 = ucell.tpiba * ucell.omega * rhocg[pw.ig2ngg[ig]]
                                                    * pw.gcar[ig].z * conj(psiv[pw.ig2fftc[ig]]) 
													* complex<double>(sin(arg), cos(arg)) ;

							this->fcc[iat][2] += ipol2.real();
                        }
                        ++iat;
                    }
                }
            }
        }
    }
    delete [] rhocg;
	delete [] psiv; // mohan fix bug 2012-03-22

	// need to be improved here, otherwise the calling time is
	// to many!
	for(int iat=0; iat<ucell.nat; iat++)
	{
    	Parallel_Reduce::reduce_double_all(this->fcc[iat], 3);
	}
	timer::tick("Force_LCAO","cal_force_cc",'E');
	return;
}


void Force_LCAO::cal_force_scc(void)
{
    TITLE ("Force_LCAO", "cal_force_scc");
    timer::tick ("Force_LCAO", "cal_force_scc",'E');

//	cout << " calculate scc force" << endl;

    complex<double>* psic = new complex<double> [pw.nrxx];
    ZEROS(psic, pw.nrxx);

    if (NSPIN == 1 || NSPIN == 4)
    {
        for (int i = 0; i < pw.nrxx; i++)
        {
            psic[i] = pot.vnew(0,i);
        }
    }
    else
    {
        int isup = 0;
        int isdw = 1;
        for (int i = 0; i < pw.nrxx; i++)
        {
            psic[i] = (pot.vnew(isup, i) + pot.vnew(isdw, i)) * 0.5;
        }
    }

    int ndm = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        if (ndm < ucell.atoms[it].msh)
        {
            ndm = ucell.atoms[it].msh;
        }
    }

    //work space
    double* aux = new double[ndm];
    ZEROS(aux, ndm);

    double* rhocgnt = new double[pw.nggm];
    ZEROS(rhocgnt, pw.nggm);

    pw.FFT_chg.FFT3D(psic, -1);

    double fact = 2.0;
    for (int nt = 0; nt < ucell.ntype; nt++)
    {
//		Here we compute the G.ne.0 term
        const int mesh = ucell.atoms[nt].msh;
        for (int ig = pw.gstart; ig < pw.nggm; ig++)
        {
            const double gx = sqrt(pw.ggs[ig]) * ucell.tpiba;
            for (int ir = 0; ir < mesh; ir++)
            {
                if (ucell.atoms[nt].r[ir] < 1.0e-8)
                {
                    aux[ir] = ucell.atoms[nt].rho_at[ir];
                }
                else
                {
                    const double gxx = gx * ucell.atoms[nt].r[ir];
                    aux[ir] = ucell.atoms[nt].rho_at[ir] * sin(gxx) / gxx;
                }
            }
            Mathzone::Simpson_Integral(mesh , aux, ucell.atoms[nt].rab , rhocgnt [ig]);
        }

        int iat = 0;
        for (int it = 0; it < ucell.ntype; it++)
        {
			Atom* atom = &ucell.atoms[it];
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                if (nt == it)
                {
                    for (int ig = pw.gstart; ig < pw.ngmc; ig++)
                    {
                        const double arg = TWO_PI * (pw.gcar[ig].x * atom->tau[ia].x
                                                     + pw.gcar[ig].y * atom->tau[ia].y
                                                     + pw.gcar[ig].z * atom->tau[ia].z);

                        const complex<double> cpm = complex<double>(sin(arg), cos(arg)) * conj(psic[pw.ig2fftc[ig] ]);

                        fscc[iat][0] += fact * rhocgnt [pw.ig2ngg[ig] ] * ucell.tpiba * pw.gcar[ig].x * cpm.real();
                        fscc[iat][1] += fact * rhocgnt [pw.ig2ngg[ig] ] * ucell.tpiba * pw.gcar[ig].y * cpm.real();
                        fscc[iat][2] += fact * rhocgnt [pw.ig2ngg[ig] ] * ucell.tpiba * pw.gcar[ig].z * cpm.real();
                    }
                }
	//			cout << " fscc=" << fscc[iat][0] << " " << fscc[iat][1] << " " << fscc[iat][2] << endl;
                iat++;
            }
        }
    }

    for (int iat = 0; iat < ucell.nat; iat++)
    {
        Parallel_Reduce::reduce_double_pool(this->fscc[iat], 3);
    }

	delete[] psic; // mohan fix bug 2012-03-22
	delete[] aux; // mohan fix bug 2012-03-22
	delete[] rhocgnt; // mohan fix bug 2012-03-22
    timer::tick ("Force_LCAO", "cal_force_scc",'E');
    return;
}


void Force_LCAO::cal_stress(matrix &stress)
{
	TITLE("Force_LCAO","cal_stress");

	Stress_LCAO SS;
	SS.allocate();
	SS.start_stress(this->soverlap, this->stvnl_dphi, this->svnl_dbeta, this->svl_dphi, this->stress_vdw);

	double unit_transform = 0.0;
	unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * 1.0e-8;
	double external_stress[3] = {PRESS1,PRESS2,PRESS3};

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			stress(i,j) = SS.scs[i][j];

            //quxin added for DFT+U; stress contribution from DFT+U
            if(INPUT.dft_plus_u) if(i!=j) stress(i,j) += dftu.stress_dftu.at(i).at(j);          
		}
		stress(i,i) = SS.scs[i][i] - external_stress[i]/unit_transform;

        if(INPUT.dft_plus_u) stress(i,i) += dftu.stress_dftu.at(i).at(i);
	}
    PRESSURE = (SS.scs[0][0]+SS.scs[1][1]+SS.scs[2][2])/3;

    if(INPUT.dft_plus_u) 
	{
		PRESSURE += (dftu.stress_dftu.at(0).at(0) + dftu.stress_dftu.at(1).at(1) + dftu.stress_dftu.at(2).at(2))/3.0;
	}

    return;
}
