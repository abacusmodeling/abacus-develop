#include "gint_gamma.h"
#include "grid_technique.h"
#include "lcao_orbitals.h"
#include "../src_pw/global.h"
#include "lcao_nnr.h"

#include "global_fp.h" // mohan add 2021-01-30

void Gint_k::cal_vna(const double* vrs1, const Grid_Technique &gt)
{

	TITLE("Gint_k","cal_vna");

    if(!pvnapR_alloc_flag)
    {
        WARNING_QUIT("Gint_k::destroy_pvpR","pvnapR has not been allocated yet!");
    }
    else
    {
        // reduce the dimension of array, which
        // is used to save <phi | Vl | phi>
        if(this->reduced)
        {
            ZEROS(this->pvnapR_reduced, LNNR.nnrg);
        }
    }		

    timer::tick("Gint_k","cal_vna");

    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    double delta_r = ORB.dr_uniform;
    // possible max atom number in real space grid.
    const int max_size = gt.max_atom;
    // how many meshcells in bigcell.
    const int bxyz = pw.bxyz;

    double*** dr;// vectors between atom and grid: [bxyz, maxsize, 3]
    double** distance; // distance between atom and grid: [bxyz, maxsize]
    double*** psir_ylm; // psi(r-R) * Ylm(r-R)
    bool** cal_flag;
    double* ylma;

    if(max_size!=0)
    {
        // save the small box information for a big box.
        dr = new double**[bxyz];
        distance = new double*[bxyz];
        psir_ylm = new double**[bxyz];
        cal_flag = new bool*[bxyz];

        // mohan fix bug 2011-05-02
        // possible number of atom configureation (l,m)
        int nn = 0;
        for(int it=0; it<ucell.ntype; it++)
        {
            nn = std::max(nn, (ucell.atoms[it].nwl+1)*(ucell.atoms[it].nwl+1));
        }

        // ylma
        ylma = new double[nn];
        ZEROS(ylma, nn);

        for(int i=0; i<bxyz; i++)
        {
            // possible max atom number in a big box.
            dr[i] = new double*[max_size];
            psir_ylm[i] = new double*[max_size];
            distance[i] = new double[max_size];
            cal_flag[i] = new bool[max_size];

            ZEROS(distance[i], max_size);
            ZEROS(cal_flag[i], max_size);

            for(int j=0; j<max_size; j++)
            {
                dr[i][j] = new double[3];
                psir_ylm[i][j] = new double[ucell.nwmax];
                ZEROS(dr[i][j],3);
                ZEROS(psir_ylm[i][j],ucell.nwmax);
            }
        }
    }

	int ncxyz=0;
	if(VNA>0)
	{
		ncxyz = pw.ncxyz * VNA * VNA * VNA;
	}
	else if(VNA==0)
	{
		ncxyz = pw.ncxyz;
	}
    const double dv = ucell.omega/(double)ncxyz;
    int vl_index=0;

    // array to store local potential for each small box in
    // a big box.
    double* vldr3 = new double[bxyz];
    ZEROS(vldr3, bxyz);

    for(int i=0; i<nbx; i++)
    {
        for(int j=0; j<nby; j++)
        {
            // count the z according to big box.
            for(int k=nbz_start; k<nbz_start+nbz; k++)
            {
                const int grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
                const int size = gt.how_many_atoms[ grid_index ];
                if(size==0) continue;

                //---------------------------------
                // get the wave functions in this
                // grid.
                //---------------------------------
                this->set_ijk_atom(grid_index, size,
                psir_ylm, dr, cal_flag,
                distance, ylma, delta_r);

                int bindex = 0;
                // z is the fastest,
                for(int ii=0; ii<pw.bx; ii++)
                {
                    for(int jj=0; jj<pw.by; jj++)
                    {
                        for(int kk=0; kk<pw.bz; kk++)
                        {
                            const int iii = i*pw.bx + ii;
                            const int jjj = j*pw.by + jj;
                            const int kkk = k*pw.bz + kk;
                            vl_index = (kkk-pw.nczp_start) + jjj*pw.nczp + iii*pw.ncy*pw.nczp;
                            vldr3[bindex] = pot.vrs1[ vl_index ] * dv;
                            ++bindex;
                        }
                    }
                }

                if(this->reduced)
                {
                    this->evaluate_pvpR_reduced(this->pvnapR_reduced,grid_index, size,i,j,k,
                        psir_ylm, cal_flag, vldr3, distance, gt);
                }

            }// int k
        }// int j
    } // int i

	delete[] vldr3;
    if(max_size!=0)
    {
        for(int i=0; i<pw.bxyz; i++)
        {
            for(int j=0; j<max_size; j++)
            {
                delete[] dr[i][j];
                delete[] psir_ylm[i][j];
            }
            delete[] dr[i];
            delete[] distance[i];
            delete[] psir_ylm[i];
            delete[] cal_flag[i];
        }
        delete[] dr;
        delete[] distance;
        delete[] psir_ylm;
        delete[] cal_flag;

        delete[] ylma;
    }
    timer::tick("Gint_k","cal_vna");
	
	return;
}




void Gint_k::evaluate_pvnapR_reduced(const int &grid_index, const int &size, const int &i, const int &j, const int &k,
	double*** psir_ylm, bool** cal_flag, double* vldr3, double** distance)
{
	double *psi1, *psi2;
	double *iw1p, *iw2p;
	double *end1, *end2;
	double *pvp;
	int iw1_lo, iw2_lo;
	int iwi, iww;
	double vpsir1;

	bool *all_out_of_range = new bool[size];
	for(int ia=0; ia<size; ia++)
	{
		all_out_of_range[ia] = true;
		for(int ib=0; ib<pw.bxyz; ib++)
		{
			if(cal_flag[ib][ia])
			{
				all_out_of_range[ia] = false;
			}
		}
	}


	for (int ia1=0; ia1<size; ia1++)
	{
		if(all_out_of_range[ia1]) continue;

        const int mcell_index1 = GridT.bcell_start[grid_index] + ia1;
		const int iat = GridT.which_atom[mcell_index1];
        const int T1 = ucell.iat2it[iat];
        const int I1 = ucell.iat2ia[iat];
        const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
		const int iw1_start = GridT.trace_lo[start1];
        Atom *atom1 = &ucell.atoms[T1];
	
        //~~~~~~~~~~~~~~~~
        // get cell R1.
        //~~~~~~~~~~~~~~~~
        const int id1 = GridT.which_unitcell[mcell_index1];
        const int R1x = GridT.ucell_index2x[id1];
        const int R1y = GridT.ucell_index2y[id1];
        const int R1z = GridT.ucell_index2z[id1];
        const int DM_start = LNNR.nlocstartg[iat];

        // get (j,beta,R2)
        for (int ia2=0; ia2<size; ia2++)
        {
			if(all_out_of_range[ia2]) continue;

			//---------------------------------------------
			// check if we need to calculate the big cell.
			//---------------------------------------------
			bool same_flag = false;
			for(int ib=0; ib<pw.bxyz; ib++)
			{
				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
				{
				//	cout << " same flag is = " << ib << endl;
					same_flag = true;
					break;
				}
			} 

			if(!same_flag) continue;

            const int bcell2 = GridT.bcell_start[grid_index] + ia2;
            const int T2 = ucell.iat2it[ GridT.which_atom[bcell2]];
			const int iat2 = GridT.which_atom[bcell2];

            if (T2 >= T1)
            {
                Atom *atom2 = &ucell.atoms[T2];
                const int I2 = ucell.iat2ia[iat2];
                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
				const int iw2_start = GridT.trace_lo[start2];

	            //~~~~~~~~~~~~~~~~
                // get cell R2.
                //~~~~~~~~~~~~~~~~
                const int id2 = GridT.which_unitcell[bcell2];
                const int R2x = GridT.ucell_index2x[id2];
                const int R2y = GridT.ucell_index2y[id2];
                const int R2z = GridT.ucell_index2z[id2];

				//------------------------------------------------
				// calculate the 'offset': R2 position relative
				// to  R1 atom.
				//------------------------------------------------
                const int dRx = R1x - R2x;
                const int dRy = R1y - R2y;
                const int dRz = R1z - R2z;
	
				const int index = LNNR.cal_RindexAtom(dRx, dRy, dRz, iat2);
                int offset = -1;

				int* find_start = LNNR.find_R2[iat];
				int* findend = LNNR.find_R2[iat] + LNNR.nad[iat];
				
				// the nad should be a large expense of time.
				for(int* find=find_start; find < findend; find++)
				{
					if( find[0] == index )
					{
						offset = find - find_start;
						break;
					}
				}

                assert(offset < LNNR.nad[iat]);

				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				// what I do above is to get 'offset' for atom pair (iat1, iat2)
				// if I want to simplify this searching for offset,
				// I should take advantage of GridT.which_unitcell.
				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				const int iatw = DM_start + LNNR.find_R2st[iat][offset];
				
				for(int ib=0; ib<pw.bxyz; ib++)
				{
					if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
					{
						psi1 = psir_ylm[ib][ia1];
						psi2 = psir_ylm[ib][ia2];
						end1 = psi1 + atom1->nw;
						end2 = psi2 + atom2->nw;
						iw1_lo = iw1_start;	
						//------------------------------------
						// circle for wave functions of atom 1.
						//------------------------------------
						iwi = 0;
						for (iw1p=psi1; iw1p < end1; ++ iw1p)
						{
							vpsir1 = iw1p[0] * vldr3[ib];
							iw2_lo = iw2_start;
							iww = iatw + iwi;// -1 because ++iww from below.
							iwi += atom2->nw;

							//note! here is pvnapR, not pvpR
							//mohan add 2012-06-22
							pvp = &pvnapR_reduced[iww]; 

							//------------------------------------
							// circle for wave functions of atom 2.
							//------------------------------------
							for(iw2p=psi2; iw2p < end2; ++ iw2p)
							{
								if( iw1_lo > iw2_lo)
								{
									++iw2_lo;
									++pvp;
									continue;
								}
								pvp[0] += vpsir1 * iw2p[0];
								++iw2_lo;
								++pvp;
							}
							++iw1_lo;
						}// iw
					}//end flag
				}//end ib
            }// T
        }// ia2
	}

	delete[] all_out_of_range;
	return;
}
