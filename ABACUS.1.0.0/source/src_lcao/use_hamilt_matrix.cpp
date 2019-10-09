#include "../src_pw/global.h"
#include "use_hamilt_matrix.h"
#include "build_st_pw.h"
#include "sltk_atom_arrange.h"
#include "bfield.h"
#include "lcao_vna.h"

Use_Hamilt_Matrix::Use_Hamilt_Matrix()
{ 
    init_s = false;
}

Use_Hamilt_Matrix::~Use_Hamilt_Matrix()
{
	if(test_deconstructor)
	{
		cout << " ~Use_Hamilt_Matrix()" << endl;	
	}
}

void Use_Hamilt_Matrix::set_ion(void)
{
	TITLE("Use_Hamilt_Matrix","set_ion");
	timer::tick("Use_Hamilt_Matrix","set_ion",'D');
	if(GAMMA_ONLY_LOCAL)
	{
		// mohan add 2012-03-29
		if(GRID_SPEED==1)
		{
			// calculate the grid integration of 'Vl' matrix for gamma algorithms.
			this->GG.prepare(ucell.latvec, ucell.lat0);
		}
		else if(GRID_SPEED==2)
		{
			this->GG.prepare(ucell.latvec, ucell.lat0);
			this->GS.prepare(ucell.latvec, ucell.lat0);
		}
		else
		{
			WARNING_QUIT("Use_Hamilt_Matrix::set_ion","GRID_SPEED should be 1 or 2");
		}
	
		// calulate the 'S', 'T' and 'Vnl' matrix for gamma algorithms.
		if(!BFIELD)
		{
	    	this->calculate_STNR_gamma();	
		}
		else
		{
			this->calculate_STNR_gamma_B();
		}

		// vna only need to do once, you can check this in lcao_vna.cpp
		if(VNA>=1)
		{
			LCAO_Vna lv;
			const char matrix_type='T';
			lv.dense_vna(matrix_type); 
		}
	}
	else
	{
		// calculate the 'S', 'T' and 'Vnl' matrix for k-points algorithms.
		this->calculate_STNR_k();

		// calculate the grid integration of 'Vl' matrix for l-points algorithms.
		this->GK.init(pw.nbx, pw.nby, pw.nbzp, pw.nbzp_start, pw.ncxyz);

		// do real space grid integration.
		if(VNA>=1)
		{
			LCAO_Vna lv;
			const char matrix_type='T';
			lv.dense_vna(matrix_type);
		}
		else if(VNA==-1)
		{
			LCAO_Vna lv;
			lv.two_center_vna();
		}
	}

	// initial the overlap matrix is done.	
    this->init_s = true;
	//cout << " init_s=" << init_s << endl; //delete 2015-09-06, xiaohui
//	OUT(ofs_running,"init_s",init_s);

	timer::tick("Use_Hamilt_Matrix","set_ion",'D');
	return;
}

void Use_Hamilt_Matrix::calculate_Hgamma(void)
{
	TITLE("Use_Hamilt_Matrix","calculate_Hgamma");
	timer::tick("Use_Hamilt_Matrix","cal_Hgamma",'F');

	// Set the matrix 'H' to zero.
	if(BFIELD)
	{
		LM.zeros_HSk('H'); // 3 stands for Hloc.
	}
	else
	{
		LM.zeros_HSgamma('H'); // 3 stands for Hloc.
	}

	bool local_pw = false;

	if(local_pw)
	{
		cout << "\n Call build_H in plane wave basis!" << endl;
	 	// Use plane wave basis to calculate 'Vl' matrix. 
		Build_ST_pw bsp;
		// 0 stands for, 0 stands for k point.
		bsp.set_local(0);
	}
	else
	{
		time_t time_vlocal_start = time(NULL);
	
		// calculate the 'Vl' matrix using gamma-algorithms.
		if(VL_IN_H)
		{	
			if(VNA==0)
			{
				if(GRID_SPEED==1)
				{
					this->GG.cal_vlocal(pot.vrs1);
				}
				else if(GRID_SPEED==2)
				{
					this->GS.cal_vlocal(pot.vrs1);
				}
			}
			else if(VNA>1)
			{
				//------------------------------------------------
				// calculate the Vna from dense grid integration.
				//------------------------------------------------
				LCAO_Vna lv;
				lv.smooth_vl2();
				//LCAO_Vna::real_space_vna();
			}
			else if(VNA==1)
			{
				//------------------------------------------------
				// calculat the Vna from normal real space grid.
				//------------------------------------------------
				LCAO_Vna::smooth_vl1();
			}
		}

		time_t time_vlocal_end = time(NULL);
		OUT_TIME("vlocal integration",time_vlocal_start,time_vlocal_end);
	}

	//add T+VNL+Vl matrix.
	if(BFIELD)
	{
		LM.update_Hloc2();
	}
	else
	{
		LM.update_Hloc();
	}

	//@@@@@@@
	//test
	//@@@@@@@
	if(NURSE)
	{
		if(BFIELD)
		{
			LM.print_HSk('S');
			LM.print_HSk('T');
			LM.print_HSk('H','R'); // 
			LM.print_HSk('H','I'); // 
		}
		else
		{
			LM.print_HSgamma('S'); // S
			LM.print_HSgamma('T');
			LM.print_HSgamma('H');
		}
	//	WARNING_QUIT("Use_Hamilt_Matrix::calculate_Hgamma","print the H,S matrix");
//		QUIT();
	}


	timer::tick("Use_Hamilt_Matrix","cal_Hgamma",'F');
	return;
}


void Use_Hamilt_Matrix::calculate_STNR_gamma_B(void)
{
	TITLE("Use_Hamilt_Matrix","calculate_STNR_gamma_B");

	// s matrix
	bfid.cal_A_of_Atom();
	bfid.make_table();
	LM.zeros_HSk('S');

	// calculate overlap 
	this->GG.cal_S_T_AP('S', GridT);
	//this->UOM.calculate_S_no();	
//	LM.print_HSk('S','R');

	LM.zeros_HSk('T');

	time_t time_vnl_start = time(NULL);
	if(VNL_IN_H)
	{
		bfid.cal_A_of_Atom();
		// calculate A 
		this->GG.cal_S_T_AP('A', GridT);
		// calculate nonlocal
		bfid.calculate_NL_B();
	}
	time_t time_vnl_end = time(NULL);

	time_t time_t_start = time(NULL);
	if(T_IN_H)
	{
		// calculate T
		if(T_IN_H==1)
		{
			UOM.calculate_T_no();
		}
		else if(T_IN_H==2)
		{
			this->GG.cal_S_T_AP('T', GridT);
		}
//		LM.print_HSk('T','R');
	}
	time_t time_t_end = time(NULL);
	return;
}


void Use_Hamilt_Matrix::calculate_STNR_gamma(void)
{
	TITLE("Use_Hamilt_Matrix","calculate_fixed");

	OUT(ofs_running,"gamma_only_local",GAMMA_ONLY_LOCAL);


	// must be done after "setup_this_ion_iter"
	// because some basic parameters should be initialized
	// in UHM.GG.init();

	LM.zeros_HSgamma('S');    	
	this->UOM.calculate_S_no();	
	//LM.print_HSgamma('S');


	//-------------------------------------
	// test using plane wave calculations.
	// all the matrixs are stored in LM.
	// LM.allocate_HS_k(ParaO.nloc);
	// Build_ST_pw bsp;
	// bsp.set_ST(0, 'S');
	// LM.print_HSk('S','R',1.0e-5);
	//-------------------------------------

	// set T and Vnl matrix to zero.
	// 2 stands for LM.Hloc_fixed matrix.
	LM.zeros_HSgamma('T'); 

	//add nonlocal pseudopotential matrix element
	time_t time_vnl_start = time(NULL);
	if(VNL_IN_H)
	{
		UOM.calculate_NL_no();
	}
	time_t time_vnl_end = time(NULL);

//	OUT(ofs_running, "Time to calculate <psi|Vnl|psi>", std::difftime(time_vnl_end, time_vnl_start));
	
	//add kinetic energy matrix element
	time_t time_t_start = time(NULL);
	if(T_IN_H)
	{
		UOM.calculate_T_no();
//		LM.print_HSgamma('T');
	}
	time_t time_t_end = time(NULL);

//	ofs_running << " T+Vnl matrix" << endl;
	//LM.print_HSgamma('T');

	OUT_TIME("kinetical matrix",time_t_start, time_t_end);
	OUT_TIME("vnl matrix",time_vnl_start, time_vnl_end);

	return;
}

#include "lcao_nnr.h"
// be called in Local_Orbital_Elec::cal_bands(). 
void Use_Hamilt_Matrix::calculate_Hk(const int &ik)
{
	TITLE("Use_Hamilt_Matrix","calculate_Hk");
	timer::tick("Use_Hamilt_Matrix","calculate_Hk",'F');

	// whether you want to calculate the local potential
	// or not, you need to set this matrix to 0.
	LM.zeros_HSk('H');

	if(VL_IN_H)
	{
		//-------------------------
		// set the local potential
		// in plane wave basis.
		//-------------------------
//		Build_ST_pw bsp;
//		bsp.set_local(ik);	
//		LM.print_HSk('H','C',1.0e-5);

		//--------------------------
		// set the local potential
		// in LCAO basis.
		//--------------------------
		LM.zeros_HSR('H', LNNR.nnr);
		if(!NONCOLIN) this->GK.folding_vl_k(ik);
		else this->GK.folding_vl_k_nc(ik);
	}


	//-----------------------------------------
	// folding matrix here: S(k) (SlocR->Sloc2)
	// folding matrix here: T(k)+Vnl(k)
	// (Hloc_fixed->Hloc_fixed2)
	//-----------------------------------------
	LM.zeros_HSk('S');
	LM.zeros_HSk('T');
//	cout << " after folding Hfixed k." << endl;
	LNNR.folding_fixedH(ik);

	//------------------------------------------
	// Add T(k)+Vnl(k)+Vlocal(k)
	// (Hloc2 += Hloc_fixed2), (complex matrix)
	//------------------------------------------
//	cout << " Folding matrix here." << endl;
	LM.update_Hloc2();

/*
	if(NURSE)
	{
		LM.print_HSk('H','R',1.0e-5);
//		LM.print_HSk('S','R',1.0e-5);
	}
	*/
	
	timer::tick("Use_Hamilt_Matrix","calculate_Hk",'F');
	return;
}

// only need to do the first time.
// available for all k points.
void Use_Hamilt_Matrix::calculate_STNR_k(void)
{
    TITLE("Hamilt_Linear","calculate_STBR_k");

	//--------------------------------------------
	// set S(R) to zero.
	// the total value of S(R) in this processor
	// is LNNR.nnr.
	// and store in LM.SlocR.
	//--------------------------------------------
	LM.zeros_HSR('S', LNNR.nnr);
    this->UOM.calculate_S_no();	

	//------------------------------
	// set T(R) and Vnl(R) to zero.
	// and then calculate it
	// and store in LM.Hloc_fixedR.
	//------------------------------
	LM.zeros_HSR('T', LNNR.nnr);
	


	if(T_IN_H)
	{
		this->UOM.calculate_T_no();	
	}


	if(VNL_IN_H)
	{
		this->UOM.calculate_NL_no();
	}


	return;

	//-----------------------------------
	// this part is used for checking
	// the consistent between LCAO
	// and plane wave basis.
	//-----------------------------------
	
	// check in plane wave basis.	
	Build_ST_pw bsp;
	for(int ik=0; ik<kv.nks; ik++)
	{
		cout << " ik=" << ik << " ------------------------------------------" << endl;
	
		//----------------------------------------------
		// Check the matrix in plane wave basis.
		//----------------------------------------------
		LM.zeros_HSk('S');
		LM.zeros_HSk('T');

		bsp.set_ST(ik, 'S');
		bsp.set_ST(ik, 'T');

		cout << " --> PW S" << endl;
		LM.print_HSk('S','R',1.0e-5);
		cout << " --> PW T" << endl;
		LM.print_HSk('T','R',1.0e-5);

		string fn = "Sloc2pw.dat";
		LM.output_HSk('S', fn);
		
		//------------------------------------------
		// folding the SlocR and Hloc_fixedR matrix
		// into Sloc2 and Hloc_fixed2 matrix.
		//------------------------------------------
		LM.zeros_HSk('S');
		LM.zeros_HSk('T');
		LNNR.folding_fixedH(ik);
		cout << " --> LCAO S" << endl;
		LM.print_HSk('S','R',1.0e-5);	
		cout << " --> LCAO T+Vnl" << endl;
		LM.print_HSk('T','R',1.0e-5);	

		string fn2 = "Sloc2lcao.dat";
		LM.output_HSk('S',fn2);

		//----------------
		// test gamma Vnl	
		//----------------
//		GAMMA_ONLY_LOCAL = true;
//		LM.allocate_HS_gamma(ParaO.nloc);
//		LM.zeros_HSgamma('H');
//		UHM.UOM.calculate_NL_no( nstart );
//		GAMMA_ONLY_LOCAL = false;
//		cout << " Correct LCAO Vnl " << endl;
//		LM.print_HSgamma('H');		
//		UHM.UOM.calculate_NL_no( nstart );
//		GAMMA_ONLY_LOCAL = false;
//		cout << " Correct LCAO Vnl " << endl;
//		LM.print_HSgamma('H');		
		
	}
	
    return;
}


void Use_Hamilt_Matrix::calculate_STN_R(void)
{
    TITLE("Use_Hamilt_Matrix","calculate_STN_R");

    int iat = 0;
    int index = 0;
    Vector3<double> dtau, tau1, tau2;
    Vector3<double> dtau1, dtau2, tau0;

    LM.allocate_Hloc_fixedR_tr();
    LM.allocate_HR_tr();
    LM.allocate_SlocR_tr();

    double R_minX = GridD.getD_minX();
    double R_minY = GridD.getD_minY();
    double R_minZ = GridD.getD_minZ();

    int R_x;
    int R_y;
    int R_z;

    for(int T1 = 0; T1 < ucell.ntype; ++T1)
    {
        Atom* atom1 = &ucell.atoms[T1];
        for(int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            //GridD.Find_atom(tau1);
            GridD.Find_atom(tau1, T1, I1);
            Atom* atom1 = &ucell.atoms[T1];
            const int start = ucell.itiaiw2iwt(T1,I1,0) * NPOL;

            for(int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = GridD.getType(ad);
                const int I2 = GridD.getNatom(ad);
                Atom* atom2 = &ucell.atoms[T2];

                tau2 = GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * ucell.lat0;
                double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

                bool adj = false;

                if(distance < rcut) adj = true;
                else if(distance >= rcut)
                {
                    for(int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
                    {
                        const int T0 = GridD.getType(ad0);
                        const int I0 = GridD.getNatom(ad0);
                        const int iat0 = ucell.itia2iat(T0, I0);
                        const int start0 = ucell.itiaiw2iwt(T0, I0, 0) * NPOL;

                        tau0 = GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * ucell.lat0;
                        double distance2 = dtau2.norm() * ucell.lat0;

                        double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
                        double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

                        if( distance1 < rcut1 && distance2 < rcut2 )
                        {
                            adj = true;
                            break;
                        }
                    }
                }

                if(adj)
                {
                    const int start2 = ucell.itiaiw2iwt(T2,I2,0) * NPOL;

                    Vector3<double> dR(GridD.getBox(ad).x, GridD.getBox(ad).y, GridD.getBox(ad).z);
                    R_x = (int) (dR.x - R_minX);
                    R_y = (int) (dR.y - R_minY);
                    R_z = (int) (dR.z - R_minZ);

                    for(int ii=0; ii<atom1->nw*NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = ParaO.trace_loc_row[iw1_all];

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = ParaO.trace_loc_col[iw2_all];

                            if(nu<0)continue;

                            int iic;
                            if(KS_SOLVER=="genelpa")  // save the matrix as column major format
                            {
                                iic=mu+nu*ParaO.nrow;
                            }
                            else
                            {
                                iic=mu*ParaO.ncol+nu;
                            }

                            if(!NONCOLIN)
                            {
                                LM.SlocR_tr[R_x][R_y][R_z][iic] = LM.SlocR[index];
                                LM.Hloc_fixedR_tr[R_x][R_y][R_z][iic] = LM.Hloc_fixedR[index];
                            }
                            else
                            {
                                LM.SlocR_tr_soc[R_x][R_y][R_z][iic] = LM.SlocR_soc[index];
                                LM.Hloc_fixedR_tr_soc[R_x][R_y][R_z][iic] = LM.Hloc_fixedR_soc[index];
                            }

                            ++index;
                        }
                    }
                }
            }
        }
    }

    return;
}
