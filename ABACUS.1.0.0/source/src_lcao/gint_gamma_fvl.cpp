#include "gint_gamma.h"
#include "grid_technique.h"
#include "lcao_orbitals.h"
#include "../src_pw/global.h"

void Gint_Gamma::cal_force(const double* vlocal_in)
{
	timer::tick("Gint_Gamma","cal_force",'H');
	this->vlocal = vlocal_in;
	this->save_atoms_on_grid(GridT);
	this->gamma_force();
	timer::tick("Gint_Gamma","cal_force",'H');
}

// this subroutine lies in the heart of LCAO algorithms.
// so it should be done very efficiently, very carefully.
// I might repeat again to emphasize this: need to optimize
// this code very efficiently, very carefully.
void Gint_Gamma::gamma_force(void)
{
    TITLE("Grid_Integral","gamma_force");
    timer::tick("Gint_Gamma","gamma_force",'I');
	// GridT.lgd: local grid dimension (sub-FFT-mesh).
    double** DGridV_x = new double*[GridT.lgd];
    double** DGridV_y = new double*[GridT.lgd];
    double** DGridV_z = new double*[GridT.lgd];
    for (int i=0; i<GridT.lgd; i++)
    {
        DGridV_x[i] = new double[GridT.lgd];
        DGridV_y[i] = new double[GridT.lgd];
        DGridV_z[i] = new double[GridT.lgd];
        ZEROS(DGridV_x[i], GridT.lgd);
        ZEROS(DGridV_y[i], GridT.lgd);
        ZEROS(DGridV_z[i], GridT.lgd);
    }
    Memory::record("Gint_Gamma","DGridV",3*GridT.lgd*GridT.lgd,"double");


	//-------------------------------------------------------------
	// first part of Simpson integration.
	//-------------------------------------------------------------
	/*
	cout << " GridT.lgd = " << GridT.lgd << endl;
	cout << " simd " << GridT.lgd*GridT.lgd << endl;
	cout << " pw.ncxyz = " << pw.ncxyz << endl;
	cout << " total dim = " << GridT.lgd*GridT.lgd*pw.ncxyz << endl;
	cout << " Mem = " << Memory::record("Gint_Gamma","gamma_force",GridT.lgd*GridT.lgd*pw.ncxyz,"double") << " MB" << endl;
	int simd = GridT.lgd*GridT.lgd;
	double **simp = new double*[simd];
	for(int i=0; i<simd; i++)
	{
		simp[i] = new double[pw.ncxyz];
		ZEROS(simp[i], pw.ncxyz);
	}
	*/


    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    const double delta_r = ORB.dr_uniform;
	const Numerical_Orbital_Lm* pointer;

    double*** dr; // vectors between atom and grid: [bxyz, maxsize, 3]
    double** distance; // distance between atom and grid: [bxyz, maxsize]
    double*** psir_ylm;
	double*** dphix;
	double*** dphiy;
	double*** dphiz;
    bool** cal_flag;
    if(max_size!=0)
    {
        dr = new double**[pw.bxyz];
        distance = new double*[pw.bxyz];
        psir_ylm = new double**[pw.bxyz];
		dphix = new double**[pw.bxyz];
		dphiy = new double**[pw.bxyz];
		dphiz = new double**[pw.bxyz];
        cal_flag = new bool*[pw.bxyz];
        for(int i=0; i<pw.bxyz; i++)
        {
            dr[i] = new double*[max_size];
            distance[i] = new double[max_size];
            psir_ylm[i] = new double*[max_size];
			dphix[i] = new double*[max_size];
			dphiy[i] = new double*[max_size];
			dphiz[i] = new double*[max_size];
            cal_flag[i] = new bool[max_size];

            ZEROS(distance[i], max_size);
            ZEROS(cal_flag[i], max_size);

            for(int j=0; j<max_size; j++)
            {
                dr[i][j] = new double[3];
                psir_ylm[i][j] = new double[ucell.nwmax];
                dphix[i][j] = new double[ucell.nwmax];
                dphiy[i][j] = new double[ucell.nwmax];
                dphiz[i][j] = new double[ucell.nwmax];
                ZEROS(dr[i][j],3);
                ZEROS(psir_ylm[i][j],ucell.nwmax);
				ZEROS(dphix[i][j],ucell.nwmax);
				ZEROS(dphiy[i][j],ucell.nwmax);
				ZEROS(dphiz[i][j],ucell.nwmax);
            }
        }
    }

    int nnnmax=0;
    for(int T=0; T<ucell.ntype; T++)
    {
        nnnmax = max(nnnmax, nnn[T]);
    }

    //array to store spherical harmonics and its derivatives
	assert(nnnmax<400);
    double rly[400];
    double grly[400][3];

    double mt[3]={0,0,0};
    double *vldr3 = new double[pw.bxyz];
    ZEROS(vldr3, pw.bxyz);
	int *vindex = new int[pw.bxyz];
	ZEROS(vindex, pw.bxyz);

    for (int i=0; i< GridT.nbx; i++)
    {
        for (int j=0; j< GridT.nby; j++)
        {
            for (int k= GridT.nbzp_start; k< GridT.nbzp_start+GridT.nbzp; k++)
            {
                this->grid_index = (k-GridT.nbzp_start) + j * GridT.nbzp + i * GridT.nby * GridT.nbzp;
                const int size = GridT.how_many_atoms[ this->grid_index ];
				if(size==0)continue;
                for (int id=0; id<size; id++)
                {
                    const int mcell_index = GridT.bcell_start[grid_index] + id;
                    const int imcell = GridT.which_bigcell[mcell_index];
                    int iat = GridT.which_atom[mcell_index];
                    const int it = ucell.iat2it[ iat ];
                    const int ia = ucell.iat2ia[ iat ];
                    Atom *atom = &ucell.atoms[it];

					mt[0] = GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0];
					mt[1] = GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1];
					mt[2] = GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2];

					for(int ib=0; ib<pw.bxyz; ib++)
                    {
                        dr[ib][id][0] = GridT.meshcell_pos[ib][0] + mt[0];
                        dr[ib][id][1] = GridT.meshcell_pos[ib][1] + mt[1];
                        dr[ib][id][2] = GridT.meshcell_pos[ib][2] + mt[2];

                        distance[ib][id] = std::sqrt(dr[ib][id][0]*dr[ib][id][0] 
						+ dr[ib][id][1]*dr[ib][id][1] 
						+ dr[ib][id][2]*dr[ib][id][2]);
                        
						if(distance[ib][id] <= ORB.Phi[it].getRcut())
                        {
                            cal_flag[ib][id]=true;
                        }
                        else
                        {
                            cal_flag[ib][id]=false;
                            continue;
                        }


                    	// >>>	the old method
                    	// ylma[id] = new double[nnn[it]]; // liaochen found this bug 2010/03/29
                    	// Ylm::get_ylm_real(ucell.atoms[it].nwl+1, this->dr[id], ylma[id]);
                    	// <<<
                    	// Ylm::rlylm(ucell.atoms[it].nwl+1, dr[id].x, dr[id].y, dr[id].z, rly, grly);
                    	// Ylm::rlylm(ucell.atoms[it].nwl+1, dr[id].x, dr[id].y, dr[id].z, rly, grly);
                    	Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl, dr[ib][id][0], dr[ib][id][1], dr[ib][id][2], rly, grly);

                    	// 1E-7 is necessary in case of R is just on one grid
						// the following code is about interpolation,
						// here we use Polynomia interpolation,
						// iq, x0, x1, x3, x3, x12, x03 are all the variables 
						// needed in this interpolation method.
						// It's the same as the function in mathzone.cpp:
						// Mathzone::Polynomial_Interpolation
						// BUT THE IMPORTANT THING IS,
						// in order to minimiz the total operations,
						// we save some variable in advance here,
						// such as x12 and x03, which is only calculated 
						// once and saved here.
						//
						// 'distane' variable here is the distance between
						// this grid and the atom.
                    	const double position = distance[ib][id] / delta_r;
                    	this->iq[id] = static_cast<int>(position);
                    	this->x0[id] = position - static_cast<double>(iq[id]);
                    	this->x1[id] = 1.0 - x0[id];
                    	this->x2[id] = 2.0 - x0[id];
                    	this->x3[id] = 3.0 - x0[id];
                    	this->x12[id] = x1[id]*x2[id] / 6.0;
						this->x03[id] = x0[id]*x3[id] / 2.0;

						double tmp, dtmp;

						// circle around the localized wave functions
						// for this atom, here localized wave functions
						// has the index (n,l,m), n is the multiplicity,
						// l is the angular momentum and m is the magnetic momentum,
						// for each localized wave function 'iw', we can know
						// the n according to atom->iw2n[iw] function
						// and know the l according to atom->iw2l[iw] function.
						// also the most important thing I always emphasis,
						// is to minimize the operations!!
						// so for a set of differnt m of a specific l,
						// we don't need to do interpolation for each set of (n,l,m),
						// we only need to do interpolation for each (n,l),
						// for differemtn m, only the Spherical Harmonical functions
						// are differnt, what's why we have atom->iw2_new[iw],
						// if the localized wave funcction 'iw' reprensents a new 'l',
						// we do interpolation.
						for (int iw=0; iw< atom->nw; iw++)
						{
							// this is a new 'l', we need 1D orbital wave
							// function from interpolation method.
							if ( atom->iw2_new[iw] )
							{
								pointer = &ORB.Phi[it].PhiLN(
										atom->iw2l[iw],
										atom->iw2n[iw]);

								if ( iq[id] >= pointer->nr_uniform-4)
								{
									tmp = dtmp = 0.0;
								}
								else
								{
									// use Polynomia Interpolation method to get the 
									// wave functions
									tmp = x12[id]*(pointer->psi_uniform[iq[id]]*x3[id]
											+pointer->psi_uniform[iq[id]+3]*x0[id])
										+ x03[id]*(pointer->psi_uniform[iq[id]+1]*x2[id]
												-pointer->psi_uniform[iq[id]+2]*x1[id]);

									dtmp = x12[id]*(pointer->dpsi_uniform[iq[id]]*x3[id]
											+pointer->dpsi_uniform[iq[id]+3]*x0[id])
											+ x03[id]*(pointer->dpsi_uniform[iq[id]+1]*x2[id]
												-pointer->dpsi_uniform[iq[id]+2]*x1[id]);

					//				dtmp = x12[id]*(pointer->psi_uniform[iq[id]]*x3[id]
					//						+pointer->psi_uniform[iq[id]+3]*x0[id])
					//					+ x03[id]*(pointer->psi_uniform[iq[id]+1]*x2[id]
					//							-pointer->psi_uniform[iq[id]+2]*x1[id]);
								}
							}//new l is used.

							// get the 'l' of this localized wave function
							int ll = atom->iw2l[iw];
							int idx_lm = atom->iw2_ylm[iw];

							//special case for distance[id] -> 0
							//Problems Remained
							//You have to add this two lines
							double rr = distance[ib][id];
							if (rr < 1e-9)
							{
								// if l==0 for localized orbital,
								// here is how we choose the 3D wave function psir_ylm,
								// and the derivative at r=0 point.
								if (ll == 0)
								{
									// psir_ylm is the three dimensional localized wave functions
									// (n,l,m), which is the multiply between 1D wave function 'tmp' and
									// spherical harmonic function rly.
									psir_ylm[ib][id][iw] = tmp * rly[idx_lm];
									// the derivative of psir_ylm with respect to atom position R,
									// it's a vector, so it has x,y,z components.
									dphix[ib][id][iw] = dphiy[ib][id][iw] = dphiz[ib][id][iw] = 0.0;
								}
								// if l>0, how do we choose 3D localized wave function
								// and the derivative.
								else
								{
									pointer = &ORB.Phi[it].
										PhiLN(atom->iw2l[iw], atom->iw2n[iw]);

									double Zty = pointer->zty;
									psir_ylm[ib][id][iw] = Zty * rly[idx_lm];
									dphix[ib][id][iw] = Zty * grly[idx_lm][0];
									dphiy[ib][id][iw] = Zty * grly[idx_lm][1];
									dphiz[ib][id][iw] = Zty * grly[idx_lm][2];
								}
							}
							// if r >0, here is how we choose the 3D wave function psir_ylm,
							// and the derivative at r=0 point.
							else
							{
								double rl;
								if(ll==0)
								{
									rl = 1.0;
								}
								else if(ll==1)
								{
									rl = rr;
								}
								else 
								{
									rl = pow(rr, ll);
								}

								// 3D wave functions
								psir_ylm[ib][id][iw] = tmp * rly[idx_lm] / rl;

								// derivative of wave functions with respect to atom positions.
								double tmpdphi_rly = (dtmp  - tmp * ll / rr) / rl * rly[idx_lm] / rr;
								double tmprl = tmp/rl;

								dphix[ib][id][iw] = tmpdphi_rly * dr[ib][id][0]  + tmprl * grly[idx_lm][0];
								dphiy[ib][id][iw] = tmpdphi_rly * dr[ib][id][1]  + tmprl * grly[idx_lm][1];
								dphiz[ib][id][iw] = tmpdphi_rly * dr[ib][id][2]  + tmprl * grly[idx_lm][2];
							}
						}
					}// ib
				}//!id //finish loop of calc pre-info for each adjacent atom

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
							vindex[bindex] = (kkk-pw.nczp_start) + jjj*pw.nczp + iii*pw.ncy*pw.nczp;
							vldr3[bindex] = this->vlocal[ vindex[bindex] ] * this->vfactor;
							//vldr3[bindex] = vfactor;//for overlap test

							assert(vindex[bindex] < pw.nrxx);
							++bindex;
						}
					}
				}

				// make pair among adjacent atoms of grid with number size
				//(ia1, ia2)
				for (int ia1=0; ia1<size; ia1++)
				{
					const int mcell_index1 = GridT.bcell_start[grid_index] + ia1;
					const int iat1 = GridT.which_atom[mcell_index1];
					const int T1 = ucell.iat2it[iat1];
					const int I1 = ucell.iat2ia[iat1];
					const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
					Atom *atom1 = &ucell.atoms[T1];

					for (int ia2=0; ia2<size; ia2++)
					{
						const int mcell_index2 = GridT.bcell_start[grid_index] + ia2;
						const int iat2 = GridT.which_atom[mcell_index2];
						const int T2 = ucell.iat2it[iat2];
						const int I2 = ucell.iat2ia[iat2];
						const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
						const int iw2_lo = GridT.trace_lo[start2];
						Atom *atom2 = &ucell.atoms[T2];

						for(int ib=0; ib<pw.bxyz; ib++)
						{
							if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
							{
								double* psi1 = &psir_ylm[ib][ia1][0];
								double* dphi2x = dphix[ib][ia2];
								double* dphi2y = dphiy[ib][ia2];
								double* dphi2z = dphiz[ib][ia2];
//								double* dphi2x = &psir_ylm[ib][ia2][0]; // for overlap test
//								double* dphi2y = &psir_ylm[ib][ia2][0]; // for overlap test
//								double* dphi2z = &psir_ylm[ib][ia2][0]; // for overlap test

								int iw1_lo = GridT.trace_lo[start1];
								// how many orbitals in this type: SZ or DZP or TZP...
								for (int iw=0; iw< atom1->nw; ++iw, ++psi1)
								{
									const double vv = vldr3[ib] * psi1[0];

									double* dpvpx = &DGridV_x[iw1_lo][iw2_lo];
									double* dpvpy = &DGridV_y[iw1_lo][iw2_lo];
									double* dpvpz = &DGridV_z[iw1_lo][iw2_lo];

									//-------------------------
									// dphi: derivative of phi
									// x,y,z: direction
									// p: pointer
									//-------------------------
									double* dphi2xp = dphi2x;
									double* dphi2yp = dphi2y;
									double* dphi2zp = dphi2z;
									
									double* phi2_end = dphi2xp + atom2->nw;

									for (; dphi2xp<phi2_end; 
									++dpvpx, ++dpvpy, ++dpvpz, 
									++dphi2xp, ++dphi2yp, ++dphi2zp)
									{
										dpvpx[0] -= dphi2xp[0] * vv;
										dpvpy[0] -= dphi2yp[0] * vv;	
										dpvpz[0] -= dphi2zp[0] * vv;

	//									simp[iw1_lo*GridT.lgd+iw2_lo][vindex[ib]] -= vz;
									}
									iw1_lo++;
								}//iw1
							}// cal_flag
						}//ib
					}//!ia
				}// ia1
			}// k
		}// j
	}// i


	//-------------------------------------------------------------
	// first part of Simpson integration.
	//-------------------------------------------------------------
	/*
	int simz, simy, simx;
	if(pw.ncz%2==0) simz = pw.ncz + 1;
	else simz = pw.ncz;
	if(pw.ncy%2==0) simy = pw.ncy + 1;
	else simy = pw.ncy;
	if(pw.ncx%2==0) simx = pw.ncx + 1;
	else simx = pw.ncx;

	double* vvvz = new double[simz];
	double* rabx = new double[simx]; 
	double* raby = new double[simy]; 
	double* rabz = new double[simz]; 

	ZEROS(vvvz, simz);
	ZEROS(rabx, simx);
	ZEROS(raby, simy);
	ZEROS(rabz, simz);
	//for(int i=0; i<simx; i++) rabx[i] = ucell.a1.x / (double)pw.ncx * ucell.lat0;
	//for(int j=0; j<simy; j++) raby[j] = ucell.a2.y / (double)pw.ncy * ucell.lat0;
	//for(int k=0; k<simz; k++) rabz[k] = ucell.a3.z / (double)pw.ncz * ucell.lat0;

	for(int i=0; i<simx; i++) rabx[i] = 1 ;
	for(int j=0; j<simy; j++) raby[j] = 1 ;
	for(int k=0; k<simz; k++) rabz[k] = 1; 
	
	double* result1 = new double[simx * simy];
	double* result2 = new double[simx];
	ZEROS(result1, simx * simy);
	ZEROS(result2, simx);

	cout << setprecision(6) << endl;
	for(int io=0; io<NLOCAL; io++)
	{
		for(int jo=0; jo<NLOCAL; jo++)
		{
			int orbital = io*NLOCAL+jo;

			ZEROS(result2, simx);
			for(int i=0; i<simx; i++)
			{
				const int ii = i%pw.ncx;
				ZEROS(&result1[i*simy], simy);
				for(int j=0; j<simy; j++)
				{	
					const int jj=j%pw.ncy;
					ZEROS(vvvz, simz);

					const int index = i*simy + j;
					for(int k=0; k<simz; k++)
					{
						const int kk=k%pw.ncz;
						vvvz[k] = simp[orbital][kk + jj * pw.ncz + ii * pw.ncy * pw.ncz]; 
			//			cout << " k=" << k << " vvvz=" << simp[orbital][kk + jj * pw.ncz + ii * jj * kk] << endl;
					}
					

//					for(int k=0; k<pw.ncz; k++) result1[index]+=vvvz[k];
					Mathzone::Simpson_Integral(simz,vvvz,rabz,result1[index]);
			//				cout << " result1[" << index << "]=" << result1[index] << endl;
			//				int ok; cin >> ok;
				}

				//for(int j=0; j<pw.ncy; j++) result2[i]+=result1[i*simy+j];
				Mathzone::Simpson_Integral(simy,&result1[i*simy],raby,result2[i]);
			}
			double v = 0.0;
			//for(int i=0; i<pw.ncx; i++) v += result2[i];
			Mathzone::Simpson_Integral(simx,result2,rabx,v);
		
			double t = 0.0;
			for(int m=0; m<pw.ncxyz; m++)
			{
				t += simp[orbital][m];
			}	

			//DGridV_z[io][jo] = v;
			
			//cout << " v=" << v << " gridV=" << DGridV_z[io][jo] << " t=" << t << endl;
//			int ok; cin >> ok;
		}
	}

	delete[] vvvz;
	delete[] rabx;
	delete[] raby;
	delete[] rabz;
	delete[] result1;
	delete[] result2;


	for(int i=0; i<simd; i++)
	{
		delete[] simp[i];
	}
	delete[] simp;
	*/

	delete[] vindex;

	
    if(max_size!=0)
    {
        for(int i=0; i<pw.bxyz; i++)
        {
            for(int j=0; j<max_size; j++)
            {
                delete[] dr[i][j];
                delete[] psir_ylm[i][j];
				delete[] dphix[i][j];
				delete[] dphiy[i][j];
				delete[] dphiz[i][j];
            }
            delete[] dr[i];
            delete[] distance[i];
            delete[] psir_ylm[i];
            delete[] cal_flag[i];
			delete[] dphix[i];
			delete[] dphiy[i];
			delete[] dphiz[i];
        }
        delete[] dr;
        delete[] distance;
        delete[] psir_ylm;
        delete[] cal_flag;
		delete[] dphix;
		delete[] dphiy;
		delete[] dphiz;
    }
	delete[] vldr3;

    timer::tick("Gint_Gamma","gamma_force",'I');
    timer::tick("Gint_Gamma","gamma_force2",'I');


    for (int i=0; i<NLOCAL; i++)
    {
        double* tmpx = new double[NLOCAL];
        double* tmpy = new double[NLOCAL];
        double* tmpz = new double[NLOCAL];
        
		ZEROS(tmpx, NLOCAL);
		ZEROS(tmpy, NLOCAL);
		ZEROS(tmpz, NLOCAL);
        
		const int mu = GridT.trace_lo[i];
		// mohan fix bug 2010-09-05
		// lack mu>=0 and nu>=0 in previous version.
		if(mu >=0)
		{

        	for (int j=0; j<NLOCAL; j++)
       		{
				const int nu = GridT.trace_lo[j];
				if(nu>=0)
				{
            		tmpx[j] = DGridV_x[mu][nu];
            		tmpy[j] = DGridV_y[mu][nu];
            		tmpz[j] = DGridV_z[mu][nu];
				}
        	}
		}

		// There may be overlaps of tmpx,y,z between different
		// processors, however, the true value is the sum of it.
		// so it would be totally correct.
        Parallel_Reduce::reduce_double_pool( tmpx, NLOCAL );
        Parallel_Reduce::reduce_double_pool( tmpy, NLOCAL );
        Parallel_Reduce::reduce_double_pool( tmpz, NLOCAL );

        for (int j=0; j<NLOCAL; j++)
        {
            if (!ParaO.in_this_processor(i,j))
            {
                continue;
            }
            LM.set_force (i,j,tmpx[j], tmpy[j], tmpz[j],'N');
        }

        delete[] tmpx;
		delete[] tmpy;
		delete[] tmpz;
    }

	//delete DGridV_x,y,z
	for (int i = 0; i < GridT.lgd; i++)
	{
		delete [] DGridV_x[i];
		delete [] DGridV_y[i];
		delete [] DGridV_z[i];
	}

	delete [] DGridV_x;
	delete [] DGridV_y;
	delete [] DGridV_z;

    timer::tick("Gint_Gamma","gamma_force2",'I');
    return;
}

