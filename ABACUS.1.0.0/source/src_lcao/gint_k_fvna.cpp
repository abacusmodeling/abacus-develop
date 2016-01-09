#include "gint_k.h"
#include "gint_gamma.h"
#include "grid_technique.h"
#include "lcao_orbitals.h"
#include "../src_pw/global.h"
#include "lcao_nnr.h"

void Gint_k::fvna_k_RealSpace(const Grid_Technique &gt, double** fvna_dphi)
{
	TITLE("Gint_k","fvna_k_RealSpace");
    timer::tick("Gint_k","fvna_k_RealSpace",'i');


    if(!this->reduced)
    {
        WARNING_QUIT("Gint_k::fvna_k_RealSpace","The force with k can only with reduced H.");
    }


    int nnrg = LNNR.nnrg;
    ofs_running << " LNNR.nnrg in fvna_k_RealSpace = " << LNNR.nnrg << endl;
    assert(nnrg>=0);

    // just because to make thea arrys meaningful.
    if(LNNR.nnrg == 0)
    {
        nnrg = 1;
    }


    // to store < phi | vlocal | dphi>
    double* pvdpx = new double[nnrg];
    double* pvdpy = new double[nnrg];
    double* pvdpz = new double[nnrg];
    ZEROS(pvdpx, nnrg);
    ZEROS(pvdpy, nnrg);
    ZEROS(pvdpz, nnrg);



    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    double delta_r = ORB.dr_uniform;
    // possible max atom number in real space grid.
    const int max_size = gt.max_atom;
    // how many meshcells in bigcell.

    const int bxyz = gt.bxyz;
    const int bx = gt.bx;
    const int by = gt.by;
    const int bz = gt.bz;


    double*** dr;// vectors between atom and grid: [bxyz, maxsize, 3]
    double** distance; // distance between atom and grid: [bxyz, maxsize]
    double*** psir_ylm; // psi(r-R) * Ylm(r-R)
    bool** cal_flag;
    double* ylma;

    double*** dphi_x;
    double*** dphi_y;
    double*** dphi_z;



    if(max_size!=0)
    {
        // save the small box information for a big box.
        dr = new double**[bxyz];
        distance = new double*[bxyz];
        psir_ylm = new double**[bxyz];
        cal_flag = new bool*[bxyz];
        dphi_x = new double**[bxyz];
        dphi_y = new double**[bxyz];
        dphi_z = new double**[bxyz];


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
            dphi_x[i] = new double*[max_size];
            dphi_y[i] = new double*[max_size];
            dphi_z[i] = new double*[max_size];


            ZEROS(distance[i], max_size);
            ZEROS(cal_flag[i], max_size);

            for(int j=0; j<max_size; j++)
            {
                dr[i][j] = new double[3];
                psir_ylm[i][j] = new double[ucell.nwmax];
                dphi_x[i][j] = new double[ucell.nwmax];
                dphi_y[i][j] = new double[ucell.nwmax];
                dphi_z[i][j] = new double[ucell.nwmax];

                ZEROS(dr[i][j],3);
                ZEROS(psir_ylm[i][j],ucell.nwmax);
                ZEROS(dphi_x[i][j],ucell.nwmax);
                ZEROS(dphi_y[i][j],ucell.nwmax);
                ZEROS(dphi_z[i][j],ucell.nwmax);

            }
        }
    }

	// mohan fix bug 2012-06-18	
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
    double* vna3d = new double[bxyz];
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
                ZEROS(vna3d, bxyz);

                this->set_ijk_atom_fvna(grid_index, size,
                psir_ylm, dr, cal_flag,
                distance, ylma, delta_r, 
				dphi_x, dphi_y, dphi_z,
				gt, vna3d);


				// for test, use original vl(G)->vl(r),
				// remember to turn on the potential in force_lcao_k.cpp

/*
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
                            //vldr3[bindex] = dv; // for overlap test

                            ++bindex;
                        }
                    }
                }
				*/

				for(int ii=0; ii<bxyz; ++ii)
                {
                    vldr3[ii] = vna3d[ii] * dv;
                }

                if(this->reduced)
                {
 	                this->evaluate_vl_force(grid_index, size,i,j,k,
                    psir_ylm, cal_flag, vldr3, distance,
                    dphi_x, dphi_y, dphi_z,
                    pvdpx, pvdpy, pvdpz, gt);
                }

            }// int k
        }// int j
    } // int i

    delete[] vldr3;
	delete[] vna3d;


    if(max_size!=0)
    {
        for(int i=0; i<bxyz; i++)
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

	//LM.DHloc_fixedR_x
	this->folding_force(fvna_dphi,pvdpx, pvdpy, pvdpz);


	delete[] pvdpx;
	delete[] pvdpy;
	delete[] pvdpz;

    timer::tick("Gint_k","fvna_k_RealSpace",'i');
	return;
}


