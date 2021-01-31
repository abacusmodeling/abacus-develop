#include "grid_integral.h"
#include "../src_pw/global.h"
#include "ylm.h"
#include "lcao_orbitals.h"


// use density kernel(density matrix) to calcualte the density 
// on real space grid.
void Grid_Integral::cal_rho(void)
{
	if(test>0)TITLE("Grid_Integral","cal_rho");
	timer::tick("Grid_Integral","cal_rho");
	
	// there are two job available: cal_charge or cal_local.
	this->job = cal_charge;
	
	//=========================================
	// calculte chrge density from local basis
	//=========================================
	this->deal_region_atoms();
	
	timer::tick("Grid_Integral","cal_rho");
	return;
}

void Grid_Integral::cal_vlocal( 
		const double* vlocal_in, 
		SparseMatrix &SM_in)
{
	if(test>0)TITLE("Grid_Integral","cal_vlocal");
	timer::tick("Grid_Integral","cal_vlocal");
	this->test = 0;
	this->job = cal_local;

	// a pointer point to the sparse matrix.
	this->SM = &SM_in;

	// a pointer point to the real space grid vlocal.
	this->vlocal = vlocal_in;

	//=================================================
	// calculte vlocal grid integral from local basis
	//=================================================
	this->deal_region_atoms();
	timer::tick("Grid_Integral","cal_vlocal");
   	return;
}

// be called by : Grid_Integral::cal_rho & Grid_Integral::cal_vlocal
// this program is in the core of the whole LCAO program.
// it's very important to make this part efficient.
void Grid_Integral::deal_region_atoms(void)
{
//	if(test>0)TITLE("Grid_Integral","cal_region_alpha");
	for(int T1=0; T1<ucell.ntype; T1++)
	{
		for(int I1=0; I1<ucell.atoms[T1].na; I1++)
		{
			//GridD.Find_atom(ucell.atoms[T1].tau[I1]);
			GridD.Find_atom(ucell.atoms[T1].tau[I1], T1, I1);
			if(test==1) cout << " T1=" << T1 << " I1=" << I1 << " adjacent number = " << GridD.getAdjacentNum()+1 << endl;
			
			for (int ad = 0; ad < GridD.getAdjacentNum()+1; ad++)
			{
				const int T2 = GridD.getType(ad);
				const int I2 = GridD.getNatom(ad);

				if( T2 >= T1 )
				{
					Vector3<double> R1 = ucell.atoms[T1].tau[I1];
					Vector3<double> R2 = ucell.atoms[T2].tau[I2];
					const double distance = ( R1 - R2 ).norm();
					if( distance*lat0 > ( Rcut_max[T1] + Rcut_max[T2] ) ) continue;

					// (3) get the constrained grid direct coordinate.
					Vector3<double> max_direct_coordinate;
					Vector3<double> min_direct_coordinate;
					Vector3<double> tau_max, tau_min, tau_max2, tau_min2;
					
					this->get_small_box( R1, T1, tau_max, tau_min );
					this->get_small_box( R2, T2, tau_max2, tau_min2 );

					max_direct_coordinate.x = std::min( tau_max.x, tau_max2.x );
					max_direct_coordinate.y = std::min( tau_max.y, tau_max2.y );
					max_direct_coordinate.z = std::min( tau_max.z, tau_max2.z );
	
					min_direct_coordinate.x = std::max( tau_min.x, tau_min2.x );
					min_direct_coordinate.y = std::max( tau_min.y, tau_min2.y );
					min_direct_coordinate.z = std::max( tau_min.z, tau_min2.z );

					if( min_direct_coordinate.x<0 ||
						min_direct_coordinate.y<0 ||
						min_direct_coordinate.z<0 ||
						max_direct_coordinate.x>1 ||
						max_direct_coordinate.y>1 ||
						max_direct_coordinate.z>1)
					{
						ofs_warning << "\n" << "min_direct_coordinate = " << min_direct_coordinate.x
						<< " " <<min_direct_coordinate.y
						<< " " << min_direct_coordinate.z << endl;

						ofs_warning << "\n" << "max_direct_coordinate = " << max_direct_coordinate.x
						<< " " << max_direct_coordinate.y
						<< " " << max_direct_coordinate.z << endl;
			
						WARNING_QUIT("Grid_Base::constrain_grid","out of range(<0 or >1)!");
					}
	
					if( min_direct_coordinate.z > max_direct_coordinate.z )
					{
						ofs_warning << "\n min_direct_z=" << min_direct_coordinate.z; 
						ofs_warning << "\n max_direct_z=" << max_direct_coordinate.z; 
						WARNING_QUIT("Grid_Base::constrain_grid","min_direct_coordinate.z > max_direct_coordinate.z");
					}

					this->lmax1 = ucell.atoms[T1].nwl + 1;
					this->lmax2 = ucell.atoms[T2].nwl + 1;
					this->n1 = lmax1 * lmax1;
					this->n2 = lmax2 * lmax2;	

					// (4) 8 Edge Points In Which Includes the Two Orbitals in Direct Coordinates
					this->edge_grid_points( R1*this->lat0, R2*this->lat0, max_direct_coordinate, min_direct_coordinate );
					this->iw1_all = ucell.itiaiw2iwt(T1, I1, 0);

					//cout << "\n T1=" << T1 << " I1=" << I1 << " T2=" << T2 << " I2=" << I2; 
            		this->deal_region_orbitals(
               		R1, //atomic position of the first atom
               		T1,
               		I1,
               		R2, //atomic position of the second atom
               		T2,
               		I2);
				}
			}
        }
    }
	return;
}

// be called by cal_region_atoms
// FUNCTION : deal things happen between two atoms.
void Grid_Integral::deal_region_orbitals
(
    const Vector3<double> &R1, //atomic position of the first atom
    const int &T1,
    const int &I1,
    const Vector3<double> &R2, //atomic position of the second atom
    const int &T2,
    const int &I2
)
{
	this->Rcut1 = ORB.Phi[T1].getRcut();
	this->Rcut2 = ORB.Phi[T2].getRcut();
	for (int L1 = 0; L1 < ucell.atoms[T1].nwl + 1; L1++)
	{
		for (int N1 = 0; N1 < ucell.atoms[T1].l_nchi[L1]; N1++)
		{
			for (int m1 = 0; m1 < 2*L1 + 1; m1++)
			{
				// calculate iw2_all and index1
				this->iw2_all = ucell.itiaiw2iwt(T2, I2, 0);
				this->index1 = L1*L1+m1;
				
				for (int L2 = 0; L2 < ucell.atoms[T2].nwl+1; L2++)
				{
					for (int N2 = 0; N2 < ucell.atoms[T2].l_nchi[L2]; N2++)
					{
						for (int m2 = 0; m2 < 2*L2 + 1; m2++)
						{
							
							if( iw1_all > iw2_all) 
							{
								++iw2_all;
								continue;
							}
							
							this->index2 = L2*L2+m2;
							this->pointer1 = &ORB.Phi[T1].PhiLN(L1, N1);
							this->pointer2 = &ORB.Phi[T2].PhiLN(L2, N2);
							
							switch( job )
							{
								case cal_local: 
								this->vlocal_in_small_box();
								break;

								case cal_charge: 
								this->rho_in_small_box(); 
								break;

								default: 
								WARNING_QUIT("Grid_Integral::cal_region","check job!");
							}
							++iw2_all;
						}// m2
					}// N2
				}// L2
				++iw1_all;
			}// m1
		}// N1
	}// L1
	return;
}


void Grid_Integral::rho_in_small_box(void)
{
	int count = 0;

	double factor;
	if(iw1_all == iw2_all) factor=1;
	else factor = 2;

	double* yy1_lm = this->yy1[index1];
	double* yy2_lm = this->yy2[index2];
	double dk = LOC.DM[0][iw1_all][iw2_all];

//	double sum = 0.0;
	for (int i = this->edge_min.x; i < this->edge_max.x; i++)
	{
		for (int j = this->edge_min.y; j < this->edge_max.y; j++)
		{
			for (int k = this->edge_min.z; k < this->edge_max.z; k++)
			{
				if ( this->norm1[count] < this->Rcut1 && this->norm2[count] < this->Rcut2)
				{			
					CHR.rho[0][ ijk_index[count] ] += factor *
						Mathzone::Polynomial_Interpolation( 
						pointer1->getPsiuniform(),				// Peize Lin update 2016-05-14
						pointer1->getNruniform(),
						pointer1->getDruniform(),
						norm1[count])
						*						
						Mathzone::Polynomial_Interpolation( 
						pointer2->getPsiuniform(),				// Peize Lin update 2016-05-14
						pointer2->getNruniform(),
						pointer2->getDruniform(),
						norm2[count])
						*
						yy1_lm[count]
						*
						yy2_lm[count]
						*
						dk;
				
//					sum += psi1 * yy1( this->index1, count)* psi2 * yy2( this->index2, count);
				}		
				count++;
			
			}
		}
	}
//	if( sum > 1.0e-5 )
//	cout << "\n" << " iw1_all = " << iw1_all << " iw2_all = " << iw2_all << " sum = " << sum * pw.omega / pw.ncxyz;
	return;
}

void Grid_Integral::vlocal_in_small_box(void)
{
	double v = 0.0;
	int count=0;

	double* yy1_lm = this->yy1[index1];
	double* yy2_lm = this->yy2[index2];

	for (int i = this->edge_min.x; i < this->edge_max.x; i++)
	{
		for (int j = this->edge_min.y; j < this->edge_max.y; j++)
		{
			for (int k = this->edge_min.z; k < this->edge_max.z; k++)
			{
				//timer::tick("Grid_Integral","interpo");
										
				if ( this->norm1[count] < this->Rcut1 && this->norm2[count] < this->Rcut2)
				{								
					
					v += 	Mathzone::Polynomial_Interpolation( 
							pointer1->getPsiuniform(),				// Peize Lin update 2016-05-14
							pointer1->getNruniform(),
							pointer1->getDruniform(),
							norm1[count])
							*
							Mathzone::Polynomial_Interpolation( 
							pointer2->getPsiuniform(),				// Peize Lin update 2016-05-14
							pointer2->getNruniform(),
							pointer2->getDruniform(),
							norm2[count])
							*
							yy1_lm[count]
							*
							yy2_lm[count]
							*
							this->vlocal[ ijk_index[count] ];
				}
				count++;
			}
		}
	}
	// unit : (a.u.)
	SM->set_add(iw1_all, iw2_all, v * latvec0.Det() / nxyz );
	
	if(iw1_all!=iw2_all) 
	{
		SM->reset( iw2_all, iw1_all, SM[0](iw1_all,iw2_all)  );
	}
				
	return;					
}



