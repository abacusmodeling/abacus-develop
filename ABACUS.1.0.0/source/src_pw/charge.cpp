//==========================================================
// AUTHOR : Lixin He, mohan , fangwei
// DATE : 2008-11-10
//==========================================================
//----------------------------------------------------------
// EXPLAIN : This routine calculates rhoa as the
// superposition of atomic charges.
//
// nspina is the number of spin components to be calculated
//
// nspina = 1 the total atomic charge density is calculated
// nspina = 2 the spin up and spin down atomic charge
// densities are calculated assuming an uniform atomic
// spin-polarization equal to starting_magnetization(nt)
// nspina = 4 noncollinear case. The total density is
// calculated in the first component and the magnetization
// vector in the other three.
//
// NB: nspina may not be equal to nspin because in some cases
// (as in update) the total charge only could be needed,
// even in a LSDA calculation.
//----------------------------------------------------------
#include "tools.h"
#include "global.h"
#include "charge.h"
#include "magnetism.h"
#include "../src_parallel/parallel_grid.h"

Charge::Charge()
{
	allocate_rho = false;
}

Charge::~Charge()
{
	if(allocate_rho)
	{
		for(int i=0; i<NSPIN; i++)
		{
			delete[] rho[i];
			delete[] rhog[i];
			delete[] rho_save[i];
			delete[] rhog_save[i];
		}
		delete[] rho;
		delete[] rhog;
		delete[] rho_save;
		delete[] rhog_save;
    	delete[] rho_core;
		delete[] rhog_core;
	}
}

void Charge::init()
{
    if (test_charge) TITLE("Charge","init");

	assert(allocate_rho == false);

    if (test_charge > 1)
    {
        cout << "\n spin_number = " << NSPIN
             << " real_point_number = " << pw.nrxx;
    }

	// allocate memory
	rho = new double*[NSPIN];
	rhog = new complex<double>*[NSPIN];
	rho_save = new double*[NSPIN];
	rhog_save = new complex<double>*[NSPIN];

	for(int is=0; is<NSPIN; is++)
	{
		rho[is] = new double[pw.nrxx];
		rhog[is] = new complex<double>[pw.ngmc];
		rho_save[is] = new double[pw.nrxx];
		rhog_save[is] = new complex<double>[pw.ngmc];			
		ZEROS(rho[is], pw.nrxx);
		ZEROS(rhog[is], pw.ngmc);
		ZEROS(rho_save[is], pw.nrxx);
		ZEROS(rhog_save[is], pw.ngmc);
	}

    Memory::record("Charge","rho",NSPIN*pw.nrxx,"double");
    Memory::record("Charge","rho_save",NSPIN*pw.nrxx,"double");
    Memory::record("Charge","rhog",NSPIN*pw.ngmc,"double");
    Memory::record("Charge","rhog_save",NSPIN*pw.ngmc,"double");

    this->rho_core = new double[pw.nrxx]; // core charge in real space
    ZEROS( rho_core, pw.nrxx);

	this->rhog_core = new complex<double>[pw.ngmc]; // reciprocal core charge
	ZEROS( rhog_core, pw.ngmc);

    Memory::record("Charge","rho_core",pw.nrxx,"double");
    Memory::record("Charge","rhog_core",pw.ngmc,"double");

	this->allocate_rho = true;
    return;
}

double Charge::sum_rho(void) const
{
    double sum_rho = 0.0;
    	
	for(int is=0; is<NSPIN; is++)
	{
		for(int ir=0; ir<pw.nrxx; ir++) 
		{
			sum_rho += this->rho[is][ir];
		}
	}

	// Here multiply the sum of charge density by a factor
    sum_rho *= ucell.omega / static_cast<double>( pw.ncxyz );
    Parallel_Reduce::reduce_double_pool( sum_rho );

	// mohan fixed bug 2010-01-18, 
	// sum_rho may be smaller than 1, like Na bcc.
    if (sum_rho <= 0.1)
    {
		ofs_warning << " sum_rho=" << sum_rho << endl;
        WARNING_QUIT("Charge::renormalize_rho","Can't find even an electron!");
    }
    return sum_rho;
}


void Charge::renormalize_rho(void)
{
    TITLE("Charge","renormalize_rho");

    const double sr = this->sum_rho();
	ofs_warning << setprecision(15);
	OUT(ofs_warning,"charge before normalized",sr);
    const double normalize_factor = ucell.nelec / sr;

	for(int is=0; is<NSPIN; is++)
	{
    	for(int ir=0; ir<pw.nrxx; ir++) 
		{
			rho[is][ir] *= normalize_factor;
		}
	}

	OUT(ofs_warning,"charge after normalized",this->sum_rho());

	ofs_running << setprecision(6);
    return;
}

//-------------------------------------------------------
// superposition of atomic charges contained in the array
// rho_at (read from pseudopotential files)
// allocate work space (psic must already be allocated)
//-------------------------------------------------------
void Charge::atomic_rho(const int spin_number_need, double** rho_in)const
{
    TITLE("Charge","atomic_rho");
    timer::tick("Charge","atomic_rho");

    double rho_at = 0.0;
    double x = 0.0;

	assert(ucell.meshx>0);
    double *rho1d = new double[ucell.meshx];
    
	// one dimension of charge in G space.
	double *rho_lgl= new double[ pw.nggm ];
    ZEROS(rho1d, ucell.meshx);
    ZEROS(rho_lgl, pw.nggm);

	// use interpolation to get three dimension charge density.
    ComplexMatrix rho_g3d( spin_number_need, pw.ngmc);


	// check the start magnetization
	int startmag_type = 1;
	for(int it=0; it<ucell.ntype; it++)
	{
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			if(ucell.atoms[it].mag[ia]!=0.0)
			{
				startmag_type = 2;
				break;
			}
		}
	}
	OUT(ofs_warning,"startmag_type",startmag_type);


    for (int it = 0;it < ucell.ntype;it++)
    {
		Atom* atom = &ucell.atoms[it];

        if (test_charge>0) cout<<" To get charge in G for atom type : " << it <<endl;
		// mesh point of this element.
        const int mesh = atom->msh;



        //----------------------------------------------------------
        // Here we check the electron number 
        //----------------------------------------------------------
		double* rhoatm = new double[mesh];
		for(int ir=0; ir<mesh; ++ir)
		{
			double r2=atom->r[ir]*atom->r[ir];
			rhoatm[ir]=atom->rho_at[ir]/FOUR_PI/r2;
		}

		double charge = 0.0;
		Mathzone::Simpson_Integral(atom->msh,atom->rho_at,atom->rab,charge);

		OUT(ofs_warning,"charge from rho_at",charge);
		assert(charge!=0.0);
		double scale=1.0;
		if(charge!=atom->zv)
		{
			OUT(ofs_warning,"charge should be",atom->zv);
			scale = atom->zv/charge;
		}

		for(int ir=0; ir<mesh; ++ir)
		{
			rhoatm[ir] *= scale;
			rhoatm[ir] *= (FOUR_PI*atom->r[ir]*atom->r[ir]);
		}

        //----------------------------------------------------------
        // Here we compute the G=0 term
        //----------------------------------------------------------

        if (pw.gstart == 1)
        {
            for (int ir = 0;ir < mesh;ir++)
            {
//                rho1d [ir] = atom->rho_at[ir];
				rho1d[ir] = rhoatm[ir];
            }
            Mathzone::Simpson_Integral(mesh, rho1d, atom->rab , rho_lgl[0]);
        }


        if (test_charge>0) cout<<"\n |G|=0 term done." <<endl;
        //----------------------------------------------------------
        // Here we compute the G<>0 term
        // But if in parallel case
        // G=0 term only belong to 1 cpu.
        // Other processors start from '0'
        //----------------------------------------------------------
        for (int ig = pw.gstart; ig < pw.nggm ;ig++)
        {
            const double gx = sqrt(pw.ggs [ig]) * ucell.tpiba;
            for (int ir = 0; ir < mesh;ir++)
            {
                if ( atom->r[ir] < 1.0e-8 )
                {
                    rho1d[ir] = rhoatm[ir];
                    //rho1d[ir] = atom->rho_at[ir];
                }
                else
                {
                    const double gxx = gx * atom->r[ir];
                    rho1d[ir] = rhoatm[ir] * sin(gxx) / gxx;
                    rho1d[ir] = rhoatm[ir] * sin(gxx) / gxx;
                }
            }
            Mathzone::Simpson_Integral(mesh , rho1d, atom->rab , rho_lgl [ig]);
        }
		
		
		delete[] rhoatm;
        
		
		
		
		
		
		
		
		
		
		if (test_charge>0) cout<<" |G|>0 term done." <<endl;
        //----------------------------------------------------------
        // EXPLAIN : Complete the transfer of rho from real space to
        // reciprocal space
        //----------------------------------------------------------
        for (int ig=0; ig< pw.nggm ; ig++)
        {
            rho_lgl[ig] /= ucell.omega;
        }
        //----------------------------------------------------------
        // EXPLAIN : compute the 3D atomic charge in reciprocal space
        //----------------------------------------------------------
        if(spin_number_need==1)
        {
            for (int ig=0; ig< pw.ngmc ;ig++)
            {
                rho_g3d(0, ig) += pw.strucFac(it, ig) * rho_lgl[ pw.ig2ngg[ig] ];
            }
		}
		// mohan add 2011-06-14, initialize the charge density according 
		// to each atom.
		else if(spin_number_need==2)
		{
			if(startmag_type==1)
			{
				for (int ig = 0; ig < pw.ngmc ; ig++)
				{
					const complex<double> swap = pw.strucFac(it, ig)* rho_lgl[pw.ig2ngg[ig]];
					rho_g3d(0, ig) += swap * mag.nelup_percent(it);
					rho_g3d(1, ig) += swap * mag.neldw_percent(it);
				}
			}
			// mohan add 2011-06-14
			else if(startmag_type==2)
			{
				complex<double> swap = ZERO;
				complex<double> ci_tpi = NEG_IMAG_UNIT * TWO_PI;
				for (int ia = 0; ia < atom->na; ia++)
				{
					const double up = 0.5 * ( 1 + atom->mag[ia] );
					const double dw = 0.5 * ( 1 - atom->mag[ia] );
					//cout << " atom " << ia << " up=" << up << " dw=" << dw << endl;

					for (int ig = 0; ig < pw.ngmc ; ig++)
					{
						const double Gtau = 
							pw.gcar[ig].x * atom->tau[ia].x
							+ pw.gcar[ig].y * atom->tau[ia].y
							+ pw.gcar[ig].z * atom->tau[ia].z; 

						swap = exp(ci_tpi * Gtau) * rho_lgl[pw.ig2ngg[ig]];

						rho_g3d(0, ig) += swap * up;
						rho_g3d(1, ig) += swap * dw;
					}
				}
			}
		}
		else
		{
			WARNING_QUIT("Charge::spin_number_need"," Either 1 or 2, check SPIN number !");
		}
	}

    delete [] rho_lgl;
    delete [] rho1d;;

    if (test_charge>0) cout<<"\n Charge in G space is done."<<endl;

	assert( spin_number_need > 0 );
	double* ne = new double[spin_number_need];
	ZEROS( ne, spin_number_need);
    for (int is = 0; is < spin_number_need;is++)
    {
        UFFT.ToRealSpace( is, rho_g3d, rho_in[is]);

		for(int ir=0; ir<pw.nrxx; ++ir)
		{
			ne[is] += rho_in[is][ir];
		}
		ne[is] *= ucell.omega/(double)pw.ncxyz; 
		Parallel_Reduce::reduce_double_pool( ne[is] );

        // we check that everything is correct
        double neg = 0.0;
        double rea = 0.0;
        double ima = 0.0;
		double sumrea = 0.0;
        for (int ir=0;ir < pw.nrxx; ir++)
        {
            rea = UFFT.porter[ir].real();
			sumrea += rea;
            neg += std::min(0.0, rea);
            ima += abs(UFFT.porter[ir].imag());
        }

		Parallel_Reduce::reduce_double_pool( neg );	
		Parallel_Reduce::reduce_double_pool( ima );	
		Parallel_Reduce::reduce_double_pool( sumrea );	

		// mohan fix bug 2011-04-03
        neg = neg / (double)pw.ncxyz * ucell.omega;
        ima = ima / (double)pw.ncxyz * ucell.omega;
		sumrea = sumrea / (double)pw.ncxyz * ucell.omega;

        if ((neg < -1.0e-4) && (is == 0 || NSPIN == 2) || ima > 1.0e-4)
        {
            ofs_warning << " Warning: negative or imaginary starting charge : " ;
            ofs_warning << " neg = " << neg
                 << " ima = " << ima
                 << " SPIN = " << is << endl;
        }

//		cout << " sum rho for spin " << is << " = " << sumrea << endl;
//		cout << " sum rho for spin " << is << " = " << sumrea << endl;

    }//end is

//	for(int it=0; it<ucell.ntype; it++)
//	{
		//cout << " nelup_percent = " << mag.nelup_percent(it) << endl;
		//cout << " neldw_percent = " << mag.neldw_percent(it) << endl;
//	}


	double ne_tot = 0.0;
	for(int is=0; is<spin_number_need; ++is)
	{
		ofs_warning << "\n SETUP ATOMIC RHO FOR SPIN " << is+1 << endl;
		OUT(ofs_warning,"Electron number from rho",ne[is]);
		ne_tot += ne[is];
	}
	OUT(ofs_warning,"total electron number from rho",ne_tot);
	OUT(ofs_warning,"should be",ucell.nelec);
	for(int is=0; is<spin_number_need; ++is)
	{
		for(int ir=0; ir<pw.nrxx; ++ir)
		{
			rho_in[is][ir] = rho_in[is][ir] / ne_tot * ucell.nelec;
		}
	}


	
	// if TWO_EFEMI, 
	// the total magnetism will affect the calculation of
	// occupations.
	// mag.compute_magnetization();

	//ofs_running << " Superposition of atomic wave function as First-Charge done." << endl;
	//2014-06-22
	delete[] ne;

    timer::tick("Charge","atomic_rho");
    return;
}

//==========================================================
// MEMBER FUNCTION :
// set_rhoc : computes the core charge on the real space
// 3D mesh.
// from set_rhoc.f90
//==========================================================
void Charge::set_rho_core(
    const ComplexMatrix &structure_factor
)
{
    TITLE("Charge","set_rho_core");
    timer::tick("Charge","set_rho_core");
    double eps = 1.e-10;
    en.etxcc = 0.0;
//----------------------------------------------------------
// LOCAL VARIABLES :
// counter on mesh points
// counter on atomic types
// counter on g vectors
//----------------------------------------------------------
    int ir = 0;
    int it = 0;
    int ig = 0;

    bool bl = false;
    for (it = 0; it<ucell.ntype; it++)
    {
//----------------------------------------------------------
// USE GLOBAL CLASS VARIABLE :
// NAME : nlcc( Non linear core corrections )
//----------------------------------------------------------
        if (ucell.atoms[it].nlcc)
        {
            bl = true;
            break;
        }
    }

    if (!bl)
    {
        ZEROS( this->rho_core, pw.nrxx);
    	timer::tick("Charge","set_rho_core");
        return;
    }

    double *rhocg = new double[pw.nggm];
    ZEROS(rhocg, pw.nggm );

	// three dimension.
    complex<double> *vg = new complex<double>[pw.ngmc];	

    for (int it = 0; it < ucell.ntype;it++)
    {
        if (ucell.atoms[it].nlcc)
        {
//----------------------------------------------------------
// EXPLAIN : drhoc compute the radial fourier transform for
// each shell of g vec
//----------------------------------------------------------
            this->non_linear_core_correction(
                ppcell.numeric,
                ucell.atoms[it].msh,
                ucell.atoms[it].r,
                ucell.atoms[it].rab,
                ucell.atoms[it].rho_atc,
                rhocg);
//----------------------------------------------------------
// EXPLAIN : multiply by the structure factor and sum
//----------------------------------------------------------
            for (int ig = 0; ig < pw.ngmc ; ig++)
            {
                vg[ig] += structure_factor(it, ig) * rhocg[pw.ig2ngg[ig]];
            }
        }
    }

	// for tmp use.
	for(int ig=0; ig< pw.ngmc; ig++)
	{
		this->rhog_core[ig] = vg[ig];
	}

    UFFT.ToRealSpace(vg, this->rho_core);

    // test on the charge and computation of the core energy
    double rhoima = 0.0;
    double rhoneg = 0.0;
    for (ir = 0; ir < pw.nrxx; ir++)
    {
        rhoneg += min(0.0, UFFT.porter[ir].real());
        rhoima += abs(UFFT.porter[ir].imag());
        // NOTE: Core charge is computed in reciprocal space and brought to real
        // space by FFT. For non smooth core charges (or insufficient cut-off)
        // this may result in negative values in some grid points.
        // Up to October 1999 the core charge was forced to be positive definite.
        // This induces an error in the force, and probably stress, calculation if
        // the number of grid points where the core charge would be otherwise neg
        // is large. The error disappears for sufficiently high cut-off, but may be
        // rather large and it is better to leave the core charge as it is.
        // If you insist to have it positive definite (with the possible problems
        // mentioned above) uncomment the following lines.  SdG, Oct 15 1999
    }

	// mohan fix bug 2011-04-03
	Parallel_Reduce::reduce_double_pool( rhoneg );
	Parallel_Reduce::reduce_double_pool( rhoima );

	// mohan changed 2010-2-2, make this same as in atomic_rho.
	// still lack something......
    rhoneg /= pw.ncxyz * ucell.omega;
    rhoima /= pw.ncxyz * ucell.omega;

    // calculate core_only exch-corr energy etxcc=E_xc[rho_core] if required
    // The term was present in previous versions of the code but it shouldn't
    delete [] rhocg;
    delete [] vg;
    timer::tick("Charge","set_rho_core");
    return;
} // end subroutine set_rhoc

// from drhoc.f90
void Charge::non_linear_core_correction
(
    const bool &numeric,
    const int mesh,
    const double *r,
    const double *rab,
    const double *rhoc,
    double *rhocg)
{
    if (test_charge) TITLE("charge","drhoc");
    double gx = 0.0;
    double rhocg1 = 0.0;
    double *aux;

    //    here we compute the fourier transform is the charge in numeric form
    if (numeric)
    {
        aux = new double [mesh];
        // G=0 term

        int igl0 = 0;
        if (pw.ggs [0] < 1.0e-8)
        {
            for (int ir = 0;ir < mesh; ir++)
            {
                aux [ir] = r [ir] * r [ir] * rhoc [ir];
            }
            Mathzone::Simpson_Integral(mesh, aux, rab, rhocg1);
            //rhocg [1] = fpi * rhocg1 / omega;
            rhocg [0] = FOUR_PI * rhocg1 / ucell.omega;//mohan modify 2008-01-19
            igl0 = 1;
        }

        // G <> 0 term
        for (int igl = igl0; igl < pw.nggm;igl++) 
        {
            gx = sqrt(pw.ggs [igl] * ucell.tpiba2);
            Mathzone::Spherical_Bessel(mesh, r, gx, 0, aux);
            for (int ir = 0;ir < mesh; ir++) 
            {
                aux [ir] = r[ir] * r[ir] * rhoc [ir] * aux [ir];
            } //  enddo
            Mathzone::Simpson_Integral(mesh, aux, rab, rhocg1);
            rhocg [igl] = FOUR_PI * rhocg1 / ucell.omega;
        } //  enddo
        delete [] aux;
    }
    else
    {
        // here the case where the charge is in analytic form,
        // check old version before 2008-12-9 or see pwscf
    }
    return;
}


//----------------------------------------------------------
// MEMBER FUNCTION :
// NAME : sum_band
//----------------------------------------------------------
void Charge::sum_band(void)
{
    TITLE("Charge","sum_band");
    timer::tick("Charge","sum_band");
//----------------------------------------------------------
// Calculates the symmetrized charge density and sum of
// occupied eigenvalues.
//----------------------------------------------------------

	for(int is=0; is<NSPIN; is++)
	{
		ZEROS(rho[is], pw.nrxx);
	}
	
    sum_band_k();

    // Symmetrization of the charge density (and local magnetization)
    timer::tick("Charge","sum_band");
    return;
}

void Charge::sum_band_k(void)
{
    TITLE("Charge","sum_band_k");
    en.eband = 0.0;

    complex<double>* porter = UFFT.porter;

    for (int ik = 0;ik < kv.nks;ik++)
    {
		//cout << "\n ik=" << ik;
        if (NSPIN==2) CURRENT_SPIN = kv.isk[ik];
		
        //  here we compute the band energy: the sum of the eigenvalues
        for (int ibnd = 0;ibnd < NBANDS;ibnd++)
        {
            en.eband += wf.ekb[ik][ibnd] * wf.wg(ik, ibnd);
			//cout << "\n ekb = " << wf.ekb[ik][ibnd] << " wg = " << wf.wg(ik, ibnd);

            ZEROS( porter, pw.nrxx );
            for (int ig = 0;ig < kv.ngk[ik] ; ig++)
            {
                porter[ pw.ig2fftw[wf.igk(ik, ig)] ] = wf.evc[ik](ibnd, ig);
            }
            pw.FFT_wfc.FFT3D(UFFT.porter, 1);

            const double w1 = wf.wg(ik, ibnd) / ucell.omega;

            if (w1 != 0.0)
            {
                for (int ir=0; ir<pw.nrxx; ir++) 
				{
					rho[CURRENT_SPIN][ir]+=w1* norm( porter[ir] );
				}
            }
        }
    } // END DO k_loop

#ifdef __MPI
    this->rho_mpi();
    //==================================
    // Reduce all the Energy in each cpu
    //==================================
    en.eband /= NPROC_IN_POOL;
    Parallel_Reduce::reduce_double_all( en.eband );
#endif
	
	// check how many electrons on this grid.
	/*
	double sum = 0.0;
	for(int ir=0; ir<pw.nrxx; ir++)
	{
		sum += rho1[ir];
	}
	cout << "\n sum=" << sum * ucell.omega / pw.nrxx << endl;
	*/

    return;
}

#ifdef __MPI
void Charge::rho_mpi(void)
{
	TITLE("Charge","rho_mpi");
    if (NPROC==1) return;
    timer::tick("Charge","rho_mpi");
    int ir;//counters on real space mesh point.
    int iz;//counters on z direction of fft grid.
    int ip;//counters on processors

    //=========================================
    // There are two steps to do before getting
    // the final charge:
    // (1) sum up the plane rhos in each pool.
    // (2) sum up all the rhos from all pools.
    //=========================================

    //=================================================
    // Searching in all planes, and calculated each
    // plane belong to which processor.
    // Count number of planes for each cpu in this pool
	// num_z: how many planes on processor 'ip'
    //=================================================
    int *num_z = new int[NPROC_IN_POOL];
    ZEROS(num_z, NPROC_IN_POOL);
    for (iz=0;iz<pw.nbz;iz++)
    {
        ip = iz % NPROC_IN_POOL;
        num_z[ip]++;
    }

	// mohan update 2011-04-26
	for(int ip=0; ip<NPROC_IN_POOL; ip++)
	{
		num_z[ip]*=pw.bz;
	}

    //=======================================
    // Find current number of planes (nz)
	// start_z: start position of z in 
	// processor ip.
    //=======================================
    int *start_z = new int[NPROC_IN_POOL];
    ZEROS(start_z, NPROC_IN_POOL);
    for (ip=1;ip<NPROC_IN_POOL;ip++)
    {
        start_z[ip] = start_z[ip-1]+num_z[ip-1];
    }

    //====================================================
    // Find "number of data" in each processor in each pool
    //====================================================
    int *rec = new int[NPROC_IN_POOL];
	ZEROS(rec, NPROC_IN_POOL);
    const int ncxy = pw.ncx * pw.ncy;
    for (ip=0;ip<NPROC_IN_POOL;ip++)
    {
        rec[ip] = num_z[ip]*ncxy;
    }

    //======================================================
    // Find current "index of data" in each cpu in this pool
	// also, we mean start position of data.
    //======================================================
    int *dis = new int[NPROC_IN_POOL];
	ZEROS(dis, NPROC_IN_POOL);
    for (ip=1;ip<NPROC_IN_POOL;ip++)
    {
        dis[ip]=dis[ip-1]+rec[ip-1];
    }

    //==========================================
    // Collection of rho in each pool
    // ( according to different k distribution,
    // so the rho in each pool is different
    //==========================================
    double *rho_tmp = new double[pw.nrxx];
    double *rho_tot = new double[pw.ncxyz];
    double *rho_tot_aux = new double[pw.ncxyz];
	ZEROS(rho_tot_aux, pw.ncxyz);

    for (int is=0; is< NSPIN; is++)
    {
        ZEROS(rho_tot, pw.ncxyz);

		for (ir=0;ir<pw.nrxx;ir++)
		{
			rho_tmp[ir] = this->rho[is][ir] / static_cast<double>(NPROC_IN_POOL);
		}

        MPI_Allgatherv(rho_tmp, pw.nrxx, MPI_DOUBLE, rho_tot, rec, dis, MPI_DOUBLE, POOL_WORLD);
        //=================================================================
        // Change the order of rho_tot in each pool , make them consistent
        // this is the most complicated part !!
        //=================================================================
        ZEROS(rho_tot_aux, pw.ncxyz);
        for (ip=0;ip<NPROC_IN_POOL;ip++)
        {
            for (ir=0;ir<ncxy;ir++)
            {
                for (iz=0;iz<num_z[ip];iz++)
                {
					// -------------------------------------------------
					// very carefully with the order of charge density.
					// the data (ir,iz) is now in processor 'ip'.
					// different POOL has different ordering.
					// we want to collect them in each processor
					// in a unit format,
					// and then reduce among all POOLS to yield
					// the correct charge density.
					// we know the division of 'z' is indipendent
					// in each processor, so the 'unit format'
					// must have no relationship with 'z' divide method.
					// -------------------------------------------------
					// rot_tot_aux : suitable among all pools.
					// (1) the data save along z direction.
					// (2) and each element number of group 'z data' 
					// is 'pw.ncz'
					// (3) however, the data rearrange is occured
					// between [ start_z[ip], start_z[ip]+num_z[ip] )
					// (4) start_z[ip] + iz yields correct z coordiante.
					// -------------------------------------------------
					// rot_tot: suitable for local pool.
					// (1) the data save along z direction, only 
					// in a small distance.
					// (2) however, the number of z in each processor 
					// 'ip' is num_z[ip]
					// (3) the index of data increases with the ip,
					// so the data on large 'ip' processor must 
					// have large 'start position', which we label
					// start_z[ip] * ncxy.
					// -------------------------------------------------
                    rho_tot_aux[pw.ncz*ir    + start_z[ip]      + iz]
                      = rho_tot[num_z[ip]*ir + start_z[ip]*ncxy + iz];
                }
            }
        }
        //==================================
        // Reduce all the rho in each cpu
        //==================================
        MPI_Allreduce(rho_tot_aux,rho_tot,pw.ncxyz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        
		//=====================================
        // Change the order of rho in each cpu
        //=====================================
        for (ir =0;ir<ncxy;ir++)
        {
            for (iz=0;iz<num_z[RANK_IN_POOL];iz++)
            {
                this->rho[is][num_z[RANK_IN_POOL]*ir+iz] = rho_tot[pw.ncz*ir + start_z[RANK_IN_POOL] + iz ];
            }
        }
    }
    delete[] rho_tot_aux;
    delete[] rho_tot;
    delete[] rho_tmp;

    delete[] rec;
    delete[] dis;

    delete[] num_z;
    delete[] start_z;
    timer::tick("Charge","rho_mpi");
    return;
}
#endif

void Charge::save_rho_before_sum_band(void)
{
	for(int is=0; is<NSPIN; is++)
	{
    	DCOPY( rho[is], rho_save[is], pw.nrxx);
    }
    return;
}


void Charge::write_rho(const int &is, const int &iter, const string &fn, const int &precision, const bool for_plot)
{
    TITLE("Charge","write_rho");
    if (out_charge==0) 
	{
		return;
	}
	else if(iter % out_charge != 0) 
	{
		return; // mohan add 2010-05-22
	}
	
	time_t start, end;
	ofstream ofs;
	
	if(MY_RANK==0)
	{
		start = time(NULL);
    	
		ofs.open(fn.c_str());
    	if (!ofs)
    	{
        	WARNING("Charge::write_rho","Can't create Charge File!");
    	}	

		//ofs_running << "\n Output charge file." << endl;

		ofs << ucell.latName << endl;//1
		ofs << " " << ucell.lat0 * 0.529177 << endl;
		ofs << " " << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
		ofs << " " << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
		ofs << " " << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;
		for(int it=0; it<ucell.ntype; it++)
		{
			ofs << " " << ucell.atoms[it].label;
		}
		ofs << endl;
		for(int it=0; it<ucell.ntype; it++)
		{
			ofs << " " << ucell.atoms[it].na;
		}
		ofs << endl;
		ofs << "Direct" << endl;

		for(int it=0; it<ucell.ntype; it++)
		{
			for(int ia=0; ia<ucell.atoms[it].na; ia++)
			{
				ofs << " " << ucell.atoms[it].taud[ia].x
					<< " " << ucell.atoms[it].taud[ia].y
					<< " " << ucell.atoms[it].taud[ia].z << endl;
			}
		}


		if(for_plot)
		{

		}
		else
		{
			ofs << "\n  " << NSPIN;
			if(NSPIN==1)
			{
				ofs << "\n " << en.ef << " (fermi energy)";
			}
			else if(NSPIN==2)
			{
				if(is==0)ofs << "\n " << en.ef_up << " (fermi energy for spin=1)"; 
				else if(is==1)ofs << "\n " << en.ef_dw << " (fermi energy for spin=2)";
			}
			else
			{
				WARNING_QUIT("write_rho","check nspin!");
			}
		}
		ofs << "\n  " << pw.ncx << " " << pw.ncy << " " << pw.ncz << endl;

		ofs << setprecision(precision);
		ofs << scientific;

	}

	
#ifndef __MPI
	int count=0;
	for(int k=0; k<pw.ncz; k++)
	{
		for(int j=0; j<pw.ncy; j++)
		{
			for(int i=0; i<pw.ncx; i++)
			{
				if(count%8==0) ofs << "\n";
				ofs << " " << rho_save[is][i*pw.ncy*pw.ncz + j*pw.ncz + k];
				++count;
			}
		}
	}
#else
//	for(int ir=0; ir<pw.nrxx; ir++) chr.rho[0][ir]=1; // for testing
//	ofs_running << "\n RANK_IN_POOL = " << RANK_IN_POOL;
	
	// only do in the first pool.
	if(MY_POOL==0)
	{
		// num_z: how many planes on processor 'ip'
    	int *num_z = new int[NPROC_IN_POOL];
    	ZEROS(num_z, NPROC_IN_POOL);
    	for (int iz=0;iz<pw.nbz;iz++)
    	{
        	int ip = iz % NPROC_IN_POOL;
        	num_z[ip] += pw.bz;
    	}	

		// start_z: start position of z in 
		// processor ip.
    	int *start_z = new int[NPROC_IN_POOL];
    	ZEROS(start_z, NPROC_IN_POOL);
    	for (int ip=1;ip<NPROC_IN_POOL;ip++)
    	{
        	start_z[ip] = start_z[ip-1]+num_z[ip-1];
    	}	

		// which_ip: found iz belongs to which ip.
		int *which_ip = new int[pw.ncz];
		ZEROS(which_ip, pw.ncz);
		for(int iz=0; iz<pw.ncz; iz++)
		{
			for(int ip=0; ip<NPROC_IN_POOL; ip++)
			{
				if(iz>=start_z[NPROC_IN_POOL-1]) 
				{
					which_ip[iz] = NPROC_IN_POOL-1;
					break;
				}
				else if(iz>=start_z[ip] && iz<start_z[ip+1])
				{
					which_ip[iz] = ip;
					break;
				}
			}
//			ofs_running << "\n iz=" << iz << " ip=" << which_ip[iz];
		}

		
		int count=0;
		int nxy = pw.ncx * pw.ncy;
		double* zpiece = new double[nxy];

		// save the rho one z by one z.
		for(int iz=0; iz<pw.ncz; iz++)
		{
			//	cout << "\n iz=" << iz << endl;
			// tag must be different for different iz.
			ZEROS(zpiece, nxy);
			int tag = iz;
			MPI_Status ierror;

			// case 1: the first part of rho in processor 0.
			if(which_ip[iz] == 0 && RANK_IN_POOL ==0)
			{
				for(int ir=0; ir<nxy; ir++)
				{
					// mohan change to rho_save on 2012-02-10
					// because this can make our next restart calculation lead
					// to the same dr2 as the one saved.
					zpiece[ir] = rho_save[is][ir*pw.nczp+iz-start_z[RANK_IN_POOL]];
					//						ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*pw.nczp+iz=" << ir*pw.nczp+iz;
				}
			}
			// case 2: > first part rho: send the rho to 
			// processor 0.
			else if(which_ip[iz] == RANK_IN_POOL )
			{
				for(int ir=0; ir<nxy; ir++)
				{
					//						zpiece[ir] = rho[is][ir*num_z[RANK_IN_POOL]+iz];
					zpiece[ir] = rho_save[is][ir*pw.nczp+iz-start_z[RANK_IN_POOL]];
					//						ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*pw.nczp+iz=" << ir*pw.nczp+iz;
				}
				MPI_Send(zpiece, nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
			}

			// case 2: > first part rho: processor 0 receive the rho
			// from other processors
			else if(RANK_IN_POOL==0)
			{
				MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
				//					ofs_running << "\n Receieve First number = " << zpiece[0];
			}

			// write data	
			if(MY_RANK==0)
			{
				//	ofs << "\niz=" << iz;
				// mohan update 2011-03-30
				for(int iy=0; iy<pw.ncy; iy++)
				{
					for(int ix=0; ix<pw.ncx; ix++)
					{
						if(count%8==0) ofs << "\n";
						ofs << " " << zpiece[ix*pw.ncy+iy];
						++count;
					}
				}
			}
		}// end iz
		delete[] zpiece;
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(MY_RANK==0) 
	{
		end = time(NULL);
		OUT_TIME("write_rho",start,end);
		ofs.close();
	}

    return;
}

bool Charge::read_rho(const int &is, const string &fn) //add by dwan
{
    TITLE("Charge","read_rho");
    ifstream ifs(fn.c_str());
    if (!ifs) 
	{
		ofs_running << " !!! Couldn't find the charge file !!!" << endl;
		return false;
	}
	else
	{
    	ofs_running << " Find the file, try to read charge from file." << endl;
	}

	bool quit=false;

    string name;
	ifs >> name;
    
	// check lattice constant, unit is Angstrom
	CHECK_DOUBLE(ifs,ucell.lat0 * 0.529177,quit);
    CHECK_DOUBLE(ifs,ucell.latvec.e11,quit);
    CHECK_DOUBLE(ifs,ucell.latvec.e12,quit);
    CHECK_DOUBLE(ifs,ucell.latvec.e13,quit);
    CHECK_DOUBLE(ifs,ucell.latvec.e21,quit);
    CHECK_DOUBLE(ifs,ucell.latvec.e22,quit);
    CHECK_DOUBLE(ifs,ucell.latvec.e23,quit);
    CHECK_DOUBLE(ifs,ucell.latvec.e31,quit);
    CHECK_DOUBLE(ifs,ucell.latvec.e32,quit);
    CHECK_DOUBLE(ifs,ucell.latvec.e33,quit);

	for(int it=0; it<ucell.ntype; it++)
	{
		CHECK_STRING(ifs,ucell.atoms[it].label,quit);
	}

	for(int it=0; it<ucell.ntype; it++)
	{
		CHECK_DOUBLE(ifs,ucell.atoms[it].na,quit);
	}

	string coordinate;
	ifs >> coordinate;

	for(int it=0; it<ucell.ntype; it++)
	{
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			CHECK_DOUBLE(ifs,ucell.atoms[it].taud[ia].x,quit);
			CHECK_DOUBLE(ifs,ucell.atoms[it].taud[ia].y,quit);
			CHECK_DOUBLE(ifs,ucell.atoms[it].taud[ia].z,quit);
		}
	}

	CHECK_INT(ifs, NSPIN);
	if(NSPIN == 1)
	{
		READ_VALUE(ifs, en.ef);
		ofs_running << " read in fermi energy = " << en.ef << endl;
	}
	else if(NSPIN == 2)
	{
		if(is==0)READ_VALUE(ifs, en.ef_up);
		else if(is==1)READ_VALUE(ifs, en.ef_dw);
	}
	else 
	{
		WARNING_QUIT("read_rho","check nspin!");
	}
	CHECK_INT(ifs, pw.ncx);	
	CHECK_INT(ifs, pw.ncy);	
	CHECK_INT(ifs, pw.ncz);	

#ifndef __MPI
	ofs_running << " Read SPIN = " << is+1 << " charge now." << endl;
	for(int k=0; k<pw.ncz; k++)
	{
		// consistent with the write_rho, something is
		// wrong.... but it works now.
		for(int j=0; j<pw.ncy; j++)
		{
			for(int i=0; i<pw.ncx; i++)
			{
				ifs >> rho[is][i*pw.ncy*pw.ncz + j*pw.ncz +k];
			}
		}
	}
#else
	
	const int nxy = pw.ncx * pw.ncy;
	double *zpiece = new double[nxy];
	for(int iz=0; iz<pw.ncz; iz++)
	{
		ZEROS(zpiece, nxy);
		if(MY_RANK==0)
		{
			//				ofs_running << " Read charge density iz=" << iz << endl;
			for(int j=0; j<pw.ncy; j++)
			{
				for(int i=0; i<pw.ncx; i++)
				{
					ifs >> zpiece[ i*pw.ncy + j ];
				}
			}
		}
		Pgrid.zpiece_to_all(zpiece, iz, rho[is]);
	}// iz
	delete[] zpiece;
#endif

    if(MY_RANK==0) ifs.close();
    return true;
}

double Charge::check_ne(const double *rho_in) const 
{
	double ne= 0.0;
	for(int ir=0; ir<pw.nrxx; ir++)
	{
		ne += rho_in[ir];
	}
	Parallel_Reduce::reduce_double_pool( ne );
	ne = ne * ucell.omega / (double)pw.ncxyz;
	cout << setprecision(10);
	cout << " check the electrons number from rho, ne =" << ne << endl;
	cout << setprecision(6);
	return ne;
}
