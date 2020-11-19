#include "gint_speed.h"
#include "grid_technique.h"
#include "lcao_orbitals.h"
#include "../src_pw/global.h"
#include "local_orbital_elec.h" //mohan add 2012-03-29

#include <stdexcept>

Gint_Speed::Gint_Speed()
{
	nov = 0;
}
Gint_Speed::~Gint_Speed(){}

void Gint_Speed::cal_vlocal(
    const double* vlocal_in)
{
    TITLE("Gint_Speed","cal_vlocal");
    timer::tick("Gint_Speed","cal_vlocal",'J');

	if(Local_Orbital_Elec::iter==1)
	{
		save_phi();
	}

    this->job = cal_local;
    this->vlocal = vlocal_in;
	this->max_size = GridT.max_atom;

	assert(GridT.ncxyz>0);
	this->vfactor = std::abs(this->latvec0.Det())/GridT.ncxyz;

	// call function
	this->gamma_vlocal();

    timer::tick("Gint_Speed","cal_vlocal",'J');
	
	// Peize Lin add 2016-12-03
	if(5==xcf.iexch_now)
		throw logic_error("Exx unfinished in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	else if(6==xcf.iexch_now)
		throw logic_error("Exx unfinished in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	
    return;
}


// evaluate the <phi | V | phi>
void Gint_Speed::gamma_vlocal(void)
{
    TITLE("Gint_Speed","gamma_vlocal");
	timer::tick("Gint_Speed","gamma_vlocal",'K');

	const int nwmax = ucell.nwmax;
	const int nat = ucell.nat;

	bool perform_gint = true;


	//--------------------------------------------
	// allocate the matrix <phi_i | V | phi_j>
	//--------------------------------------------
	double** GridVlocal;
	const int lgd_now = GridT.lgd;
	if(lgd_now > 0)
	{
		GridVlocal = new double*[lgd_now];
		for (int i=0; i<lgd_now; i++)
		{
			GridVlocal[i] = new double[lgd_now];
			ZEROS(GridVlocal[i], lgd_now);
		}
		Memory::record("Gint_Speed","GridVlocal",lgd_now*lgd_now,"double");
	}
	else if(lgd_now <= 0)
	{
		perform_gint = false;
	}


//	ofs_running << " lgd_now = " << lgd_now << endl;

	if(nov==0)
	{

	}
	else
	{
//		goto jump;

		int incx=1;
		int incy=1;
		double **p1 = new double*[nwmax];
		double **p2 = new double*[nwmax];
		for(int iw=0; iw<nwmax; ++iw) p1[iw] = new double[this->nov];
		for(int iw=0; iw<nwmax; ++iw) p2[iw] = new double[this->nov];

		for(int iat=0; iat<nat; ++iat)
		{
			const int it = ucell.iat2it[ iat ];
			const int ia = ucell.iat2ia[ iat ];
			const int start1 = ucell.itiaiw2iwt(it, ia, 0);
			const int iw_lo = GridT.trace_lo[start1];
			if(iw_lo<0)continue;


			const int nw1 = ucell.atoms[it].nw;

			//cout << " iat=" << iat << endl;
			for(int iat2=iat; iat2<nat; ++iat2)
			{
				const int it2 = ucell.iat2it[ iat2 ];
				const int ia2 = ucell.iat2ia[ iat2 ];
				const int start2 = ucell.itiaiw2iwt(it2, ia2, 0);
				const int iw2_lo = GridT.trace_lo[start2];
				if(iw2_lo<0)continue;

				const int nw2 = ucell.atoms[it2].nw;

				// for example
				// psi1: 1 3 3 5 7
				// psi2: 2 3 3 5 9
				int is=0;
				int is_rec=0;
				int ir2_start=0; //start ir of atom2
				int ee=0;
				bool same=false;
				for(int ir=0; ir<nr[iat]; ++ir)
				{
					if( ir>0 )
					{
						if( grid_label[iat][ir] == grid_label[iat][ir-1] )
						{
							same = true;
							ir2_start = is_rec;	
						}
						else
						{
							same = false;
							ir2_start = is;
						}
					}

					is_rec = ir2_start;

					for(int jr=ir2_start; jr<nr[iat2]; ++jr)
					{
						if( grid_label[iat2][jr] > grid_label[iat][ir] )
						{
							is=jr;
							break;
						}
						else if( grid_label[iat2][jr] == grid_label[iat][ir])
						{
							assert(ee<nov);
							const int dd = this->grid_label[iat2][jr];
							double v1 = this->vlocal[vindex[dd]] * vfactor;
							for(int iw=0; iw<nw1; ++iw)
							{
								p1[iw][ee] = phiylm[iat][iw][ir] * v1;
							}
							for(int iw=0; iw<nw2; ++iw)
							{
								p2[iw][ee] = phiylm[iat2][iw][jr];
							}
							++ee;
						}
					}// jr
				}// ir

				
				if(ee==0) continue;

				for(int iw=0; iw<nw1; ++iw)
				{
					int iw1_all = iw_lo + iw;
					for(int iw2=0; iw2<nw2; ++iw2)
					{
						int iw2_all = iw2_lo + iw2;
						if(iw1_all > iw2_all)
						{
							continue;
						}
						GridVlocal[iw1_all][iw2_all] = ddot_(&ee,p1[iw],&incx,p2[iw2],&incy);
					}
				}

				//---------------------------
				// very slow,but correct!
				//---------------------------
				/*
				   int is=0;
				   for(int i=0; i<nr[iat]; ++i)
				   {
				   for(int j=is; j<nr[iat2]; ++j)
				   {
				   if( grid_label[iat2][j] > grid_label[iat][i] )
				   {
				   is=j;
				   break;	
				   } 
				   else if( grid_label[iat2][j] == grid_label[iat][i])
				   {
				   const int dd = grid_label[iat2][j];
				   is=j+1;
				   double v1 =  this->vlocal[vindex[dd]] * vfactor; 
				   for(int iw=0; iw<nw; ++iw)
				   {
				   int iw1_all = iw_lo+iw;
				   double v2 = phiylm[iat][iw][i] * v1;
				   for(int iw2=0; iw2<nw; ++iw2)
				   {
				   int iw2_all = iw2_lo+iw2;
				   if(iw1_all>iw2_all)
				   {
				   continue;
				   }

				   GridVlocal[iw1_all][iw2_all] += phiylm[iat2][iw2][j] * v2;
				   }//end iw2
				   }//end iw
				   break;
				   }//same grid
				   }//j
				   }//i
				 */

			}//iat2
		}//iat


		for(int iw=0; iw<nwmax; ++iw) delete[] p1[iw];
		for(int iw=0; iw<nwmax; ++iw) delete[] p2[iw];
		delete[] p1;
		delete[] p2;



	}//end NOV

//jump:			// Peize Lin delete 2019-05-01

	//----------------------------
	// PART4
	// distribute the data.
	//----------------------------
#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	timer::tick("Gint_Speed","gamma_vlocal",'K');
	timer::tick("Gint_Speed","distri_vl",'K');
    if (job==cal_local)
    {
        double* tmp;
        for (int i=0; i<NLOCAL; i++)
        {
            tmp = new double[NLOCAL];
            ZEROS(tmp, NLOCAL);
            const int mu = GridT.trace_lo[i];
            if (mu >= 0)
            {
                for (int j=0; j<NLOCAL; j++)
                {
                    const int nu = GridT.trace_lo[j];
                    if (nu >=0)
                    {
                        if (mu <= nu)
                        {
                            tmp[j] = GridVlocal[mu][nu];
                        }
                        else
                        {
							//-------------------------------
							// origin:
                            // tmp[i] = GridVlocal[nu][mu];
                            // mohan fix bug 
							// 2011-01-13
							//-------------------------------
							tmp[j] = GridVlocal[nu][mu];
                        }
                    }
                }
            }
            Parallel_Reduce::reduce_double_grid( tmp, NLOCAL );
            Parallel_Reduce::reduce_double_diag( tmp, NLOCAL );
            for (int j=0; j<NLOCAL; j++)
            {
                if (!ParaO.in_this_processor(i,j))
                {
                    continue;
                }
			
				// mohan update 2011-04-15	
				if(BFIELD)
				{
					LM.set_HSk(i,j,complex<double>(tmp[j],0.0),'L');
				}
				else
				{
                	LM.set_HSgamma(i,j,tmp[j],'L');
				}
            }
            delete[] tmp;
        }
    }


	/*
	cout << "GridVlocal" << endl;
	for(int i=0; i<GridT.lgd; ++i)
	{
		for(int j=0; j<GridT.lgd; ++j)
		{
			cout << setw(15) << GridVlocal[i][j];
		}
		cout << endl;
	}
	*/


	// mohan update 2010-09-07
	if(GridT.lgd>0)
	{
		for(int i=0; i<GridT.lgd; i++)
		{
			delete[] GridVlocal[i];
		}
		delete[] GridVlocal;
	}





	timer::tick("Gint_Speed","distri_vl",'K');
    return;
}

