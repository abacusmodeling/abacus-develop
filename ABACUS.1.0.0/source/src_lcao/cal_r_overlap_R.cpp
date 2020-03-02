#include "cal_r_overlap_R.h"

cal_r_overlap_R::cal_r_overlap_R(){}

cal_r_overlap_R::~cal_r_overlap_R()
{
	if(this->allocate_psi_r_psi)
	{
		for(int ix = 0; ix < R_x_num; ix++)
		{			
			for(int iy = 0; iy < R_y_num; iy++)
			{				
				for(int iz = 0; iz < R_z_num; iz++)
				{				
					delete[] psi_r_psi[ix][iy][iz];				
				}
				delete[] psi_r_psi[ix][iy];
			}
			delete[] psi_r_psi[ix];
		}
		delete[] psi_r_psi;
	}
}

void cal_r_overlap_R::init()
{
	this->R_x_num = GridD.getCellX();
    this->R_y_num = GridD.getCellY();
    this->R_z_num = GridD.getCellZ();
	this->R_minX = (int)GridD.getD_minX();
	this->R_minY = (int)GridD.getD_minY();
	this->R_minZ = (int)GridD.getD_minZ();
	
	
	psi_r_psi = new Vector3<double> ***[R_x_num];
	for(int ix = 0; ix < R_x_num; ix++)
	{
		psi_r_psi[ix] = new Vector3<double> **[R_y_num];
		for(int iy = 0; iy < R_y_num; iy++)
		{
			psi_r_psi[ix][iy] = new Vector3<double> *[R_z_num];
			for(int iz = 0; iz < R_z_num; iz++)
			{				
				psi_r_psi[ix][iy][iz] = new Vector3<double> [ParaO.nloc];				
			}
		}
	}
	
	this->allocate_psi_r_psi = true;
	
	int Lmax_used, Lmax;
	
	MOT.allocate(
		ORB.get_ntype(),// number of atom types
		ORB.get_lmax(),// max L used to calculate overlap
		ORB.get_kmesh(), // kpoints, for integration in k space
		ORB.get_Rmax(),// max value of radial table
		ORB.get_dR(),// delta R, for making radial table
		ORB.get_dk()); // delta k, for integration in k space
		
	MOT.init_Table_Spherical_Bessel (2, 3, Lmax_used, Lmax);

	Ylm::set_coefficients();

	MGT.init_Gaunt_CH( Lmax );
	MGT.init_Gaunt( Lmax );	

	const int T = 0;  // arbitrary atom type 
	orb_r.set_orbital_info(
	ORB.Phi[T].PhiLN(0,0).getLabel(),  //atom label
	T,    //atom type
	1,    //angular momentum L
	1,    //number of orbitals of this L , just N
	ORB.Phi[T].PhiLN(0,0).getNr(),  //number of radial mesh
	ORB.Phi[T].PhiLN(0,0).getRab(), //the mesh interval in radial mesh 
	ORB.Phi[T].PhiLN(0,0).getRadial(),  // radial mesh value(a.u.)
	Numerical_Orbital_Lm::Psi_Type::Psi,
	ORB.Phi[T].PhiLN(0,0).getRadial(),  // radial wave function
	ORB.Phi[T].PhiLN(0,0).getNk(),
	ORB.Phi[T].PhiLN(0,0).getDk(),
	ORB.Phi[T].PhiLN(0,0).getDruniform(),
	false,
	true);
	
	for(int TA = 0; TA < ucell.ntype; TA++)
	{
		for (int TB = 0;  TB < ucell.ntype; TB++)
		{
			for (int LA=0; LA <= ORB.Phi[TA].getLmax() ; LA++)
			{
				for (int NA = 0; NA < ORB.Phi[TA].getNchi(LA); ++NA)
				{
					for (int LB = 0; LB <= ORB.Phi[TB].getLmax(); ++LB)
					{
						for (int NB = 0; NB < ORB.Phi[TB].getNchi(LB); ++NB)
						{
							center2_orb21_r[TA][TB][LA][NA][LB].insert( 
								make_pair(NB, Center2_Orb::Orb21(
									ORB.Phi[TA].PhiLN(LA,NA),	
									orb_r,									
									ORB.Phi[TB].PhiLN(LB,NB),
									MOT, MGT)));
						}
					}
				}
			}
		}
	}
	
	for( auto &co1 : center2_orb21_r )
		for( auto &co2 : co1.second )
			for( auto &co3 : co2.second )
				for( auto &co4 : co3.second )
					for( auto &co5 : co4.second )
						for( auto &co6 : co5.second )
							co6.second.init_radial_table();
						
	
}


void cal_r_overlap_R::out_r_overlap_R(const int nspin)
{
	Vector3<double> tau1, tau2, dtau;
	Vector3<double> origin_point(0.0,0.0,0.0);

    int R_x;
    int R_y;
    int R_z;
	
	for (int T1 = 0; T1 < ucell.ntype; ++T1)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			GridD.Find_atom(tau1, T1, I1);			
			
			for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
			{
				const int T2 = GridD.getType(ad);
				const int I2 = GridD.getNatom(ad);
				Atom* atom2 = &ucell.atoms[T2];
				const int R_direct_x = GridD.getBox(ad).x;
				const int R_direct_y = GridD.getBox(ad).y;
				const int R_direct_z = GridD.getBox(ad).z;
				
				R_x = R_direct_x - R_minX;
				R_y = R_direct_y - R_minY;
				R_z = R_direct_z - R_minZ;
				
			
				tau2 = GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
				if(distance < rcut - 1.0e-15)
				{
					Vector3<double> R_car = R_direct_x * ucell.a1 + 
											R_direct_y * ucell.a2 +
											R_direct_z * ucell.a3;
					
					for(int iw1 = 0; iw1 < atom1->nw; iw1++)
					{			
						for(int iw2 = 0; iw2 < atom2->nw; iw2++)
						{
							int orb_index_in_NLOCAL_1 = ucell.itiaiw2iwt( T1, I1, iw1 );
							int orb_index_in_NLOCAL_2 = ucell.itiaiw2iwt( T2, I2, iw2 );
							
							int irow = ParaO.trace_loc_row[orb_index_in_NLOCAL_1];
							int icol = ParaO.trace_loc_col[orb_index_in_NLOCAL_2];
							
							if(irow >= 0 && icol >= 0)
							{
								int atomType1 = iw2it(orb_index_in_NLOCAL_1);  int ia1 = iw2ia(orb_index_in_NLOCAL_1);  int N1 = iw2iN(orb_index_in_NLOCAL_1);  int L1 = iw2iL(orb_index_in_NLOCAL_1);  int m1 = iw2im(orb_index_in_NLOCAL_1); 
								int atomType2 = iw2it(orb_index_in_NLOCAL_2);  int ia2 = iw2ia(orb_index_in_NLOCAL_2);  int N2 = iw2iN(orb_index_in_NLOCAL_2);  int L2 = iw2iL(orb_index_in_NLOCAL_2);  int m2 = iw2im(orb_index_in_NLOCAL_2);

								Vector3<double> r_distance = ( ucell.atoms[atomType2].tau[ia2] - ucell.atoms[atomType1].tau[ia1] + R_car ) * ucell.lat0;							
															
								double overlap_x = -1 * sqrt(FOUR_PI/3.0) * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 1, m2 ); // m = 1
								double overlap_y = -1 * sqrt(FOUR_PI/3.0) * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 2, m2 ); // m = -1
								double overlap_z =      sqrt(FOUR_PI/3.0) * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap( origin_point, r_distance, m1, 0, m2 ); // m =0
								
								int icc = irow + icol * ParaO.nrow;
								psi_r_psi[R_x][R_y][R_z][icc] = Vector3<double>( overlap_x,overlap_y,overlap_z );
								
							}
							
						} // end iw2
						
					} // end iw1
					
				}
				
			} //  end ad
			
		} // end I1
		
	} // end T1
	
	
	// out r_overlap_R file
	ofstream out_r;
	stringstream ssh;
	ssh << global_out_dir << "data-rR-tr_SPIN" << nspin;
	if(DRANK == 0)
	{
		out_r.open(ssh.str().c_str());
		out_r << "Matrix Dimension of vector r(R): " << NLOCAL <<endl;
	}
	
	for(int ix = 0; ix < R_x_num; ix++)
    {
        int dRx = ix + R_minX;
        for(int iy = 0; iy < R_y_num; iy++)
        {
            int dRy = iy + R_minY;
            for(int iz = 0; iz < R_z_num; iz++)
            {
                int dRz = iz + R_minZ;				
				
				int ir,ic;
				for(int i = 0; i < NLOCAL; i++)
				{
					double *liner_x, *liner_y, *liner_z;
					liner_x = new double[NLOCAL];
					liner_y = new double[NLOCAL];
					liner_z = new double[NLOCAL];
					ZEROS(liner_x,NLOCAL);
					ZEROS(liner_y,NLOCAL);
					ZEROS(liner_z,NLOCAL);
					
					ir = ParaO.trace_loc_row[i];
					
					if(ir >= 0)
					{
						for(int j = 0; j < NLOCAL; j++)
						{
							ic = ParaO.trace_loc_col[j];
							if(ic >= 0)
							{
								int iic = ir + ic * ParaO.nrow;
								liner_x[j] = psi_r_psi[ix][iy][iz][iic].x;
								liner_y[j] = psi_r_psi[ix][iy][iz][iic].y;
								liner_z[j] = psi_r_psi[ix][iy][iz][iic].z;
							}
						}
					}
					
					Parallel_Reduce::reduce_double_all(liner_x,NLOCAL);
					Parallel_Reduce::reduce_double_all(liner_y,NLOCAL);
					Parallel_Reduce::reduce_double_all(liner_z,NLOCAL);
					
					if(DRANK == 0)
					{
						for(int j = 0; j < NLOCAL; j++)
						{
							if(i==0 && j==0)
							{
								out_r << dRx << " " << dRy << " " << dRz  << "    //R vector(R2 - R1,unit: lattice vector)" <<endl;
							}
							
							out_r << "(" << setw(20) << setprecision(9) << setiosflags(ios::scientific) << liner_x[j]
							      << "," << setw(20) << setprecision(9) << setiosflags(ios::scientific) << liner_y[j]
								  << "," << setw(20) << setprecision(9) << setiosflags(ios::scientific) << liner_z[j]
								  << ") "; 							
						}
						
						out_r << endl;
					}
					
					delete[] liner_x;
					delete[] liner_y;
					delete[] liner_z;
					
				
				}
				
				
			}
		}
	}
	
	if(DRANK == 0) out_r.close();
	
}

int cal_r_overlap_R::iw2it(int iw)
{
    int ic, type;
    ic = 0;
    for(int it = 0; it < ucell.ntype; it++)
	{
        for(int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i=0; i<(2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           type = it;
                        }
                        ic++;
					}
                }
			}
        }
	}
    return type;
}

int cal_r_overlap_R::iw2ia(int iw)
{
    int ic, na;
    ic = 0;
    for(int it = 0; it < ucell.ntype; it++)
	{
        for(int ia = 0; ia<  ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i = 0; i < (2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           na = ia;
                        }
                        ic++;
                    } 
				}
			}				
		}
    }
    return na;
}

int cal_r_overlap_R::iw2iL(int iw)
{
	int ic, iL;
    ic = 0;
    for(int it = 0; it < ucell.ntype; it++)
	{
        for(int ia = 0; ia<  ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i = 0; i < (2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           iL = L;
                        }
                        ic++;
                    } 
				}
			}				
		}
    }
    return iL;
}

int cal_r_overlap_R::iw2iN(int iw)
{
	int ic, iN;
    ic = 0;
    for(int it = 0; it < ucell.ntype; it++)
	{
        for(int ia = 0; ia<  ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i = 0; i < (2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           iN = N;
                        }
                        ic++;
                    } 
				}
			}				
		}
    }
    return iN;
	
}

int cal_r_overlap_R::iw2im(int iw)
{
	int ic, im;
    ic = 0;
    for(int it = 0; it < ucell.ntype; it++)
	{
        for(int ia = 0; ia<  ucell.atoms[it].na; ia++)
        {
            for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
                for(int N = 0; N < ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i = 0; i < (2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           im = i;
                        }
                        ic++;
                    } 
				}
			}				
		}
    }
    return im;
}