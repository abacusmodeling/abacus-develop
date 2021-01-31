#include "gint_speed.h"
#include "grid_technique.h"
#include "lcao_orbitals.h"
#include "../src_pw/global.h"
#include "local_orbital_elec.h" //mohan add 2012-03-29

#include "global_fp.h" // mohan add 2021-01-30


double Gint_Speed::cal_rho(void)
{
    TITLE("Gint_Speed","cal_rho");
    timer::tick("Gint_Speed","cal_rho",'J');

    this->job = cal_charge;

	assert(GridT.ncxyz>0);
	this->vfactor = std::abs(this->latvec0.Det())/GridT.ncxyz;

	// call function
	double ne = this->gamma_charge();

    timer::tick("Gint_Speed","cal_rho",'J');
	return ne;

}


double Gint_Speed::gamma_charge(void)
{
    TITLE("Gint_Speed","gamma_charge");
	timer::tick("Gint_Speed","gamma_charge",'K');

//	cout << " gamma_charge" << endl;


	const int nwmax = ucell.nwmax;
	const int nat = ucell.nat;

	// PART3: evaluate the <phi | V | phi>
	assert(allocate_phiylm);
	assert(nov>=0);


	if(nov > 0)
	{
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
				int ir2_start=0; //start position of atom 2
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
							const int dd = grid_label[iat2][jr];
							const int ig = vindex[dd];
							double sum = 0.0;
							for(int iw=0; iw<nw1; ++iw)
							{
								const int iw1_all = iw_lo + iw;
								double pp1 = phiylm[iat][iw][ir] + phiylm[iat][iw][ir];
								for(int iw2=0; iw2<nw2; ++iw2)
								{
									const int iw2_all = iw2_lo + iw2;
									double dm = LOC.DM[CURRENT_SPIN][iw1_all][iw2_all];

									if(iw1_all > iw2_all)
									{
										continue;
									}

									double tmp1 = pp1 * phiylm[iat2][iw2][jr];
									double tmp2 = tmp1 * dm;

									if(iw1_all < iw2_all)
									{
										sum += tmp2; 
									}
									else
									{
										sum += tmp2*0.5;
									}

									/*
									   if( ig > pw.nrxx )
									   {
									   cout << "afsdfasdfasdf" << endl;
									   WARNING_QUIT("test","test");
									   }
									   double sum = phiylm[iat][iw][i] * phiylm[iat2][iw2][jr] * dm;
									   */
								}
							}
							CHR.rho[CURRENT_SPIN][ig] += sum;
						}
					}// jr
				}// i
			}//iat2
		}//iat


		for(int iw=0; iw<nwmax; ++iw) delete[] p1[iw];
		for(int iw=0; iw<nwmax; ++iw) delete[] p2[iw];
		delete[] p1;
		delete[] p2;


	}//only if nov>0


	double sum = 0.0;
	for(int is=0; is<NSPIN; is++)
	{
		for (int ir=0; ir<pw.nrxx; ir++)
		{
			sum += CHR.rho[is][ir];
		}
	}
#ifdef __MPI
	Parallel_Reduce::reduce_double_pool( sum );
#endif
	OUT(ofs_warning,"charge density sumed from grid", sum * ucell.omega/ pw.ncxyz);

	double ne = sum * ucell.omega / pw.ncxyz;
	timer::tick("Gint_Speed","gamma_charge",'K');
	return ne;
}

