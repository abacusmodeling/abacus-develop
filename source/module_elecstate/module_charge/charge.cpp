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
// std::vector in the other three.
//
// NB: nspina may not be equal to nspin because in some cases
// (as in update) the total charge only could be needed,
// even in a LSDA calculation.
//----------------------------------------------------------
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "charge.h"
#include "module_elecstate/magnetism.h"
#include "module_elecstate/energy.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"
#include <vector>
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_base/libm/libm.h"
#include "module_io/rho_io.h"

Charge::Charge()
{
	allocate_rho = false;
	allocate_rho_final_scf = false; //LiuXh add 20180619
}


Charge::~Charge()
{
    this->destroy();
}

void Charge::destroy()
{
    if(allocate_rho || allocate_rho_final_scf) //LiuXh add 20180619
    {
        for(int i=0; i<GlobalV::NSPIN; i++)
        {
            delete[] rho[i];
            delete[] rhog[i];
            delete[] rho_save[i];
            delete[] rhog_save[i];
            if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                delete[] kin_r[i];
                delete[] kin_r_save[i];
            }
        }
        delete[] rho;
        delete[] rhog;
        delete[] rho_save;
        delete[] rhog_save;
        delete[] rho_core;
        delete[] rhog_core;
        if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            delete[] kin_r;
            delete[] kin_r_save;
        }
    }
}

void Charge::allocate(const int &nspin_in, const int &nrxx_in, const int &ngmc_in)
{
    ModuleBase::TITLE("Charge","allocate");

    if(allocate_rho == true)
    {
        this->destroy();
        allocate_rho = false;
    }

	assert(allocate_rho == false);

	//  mohan add 2021-02-20
	this->nspin = nspin_in;
	this->nrxx = nrxx_in;
	this->ngmc = ngmc_in;

    if (GlobalV::test_charge > 1)
    {
        std::cout << "\n spin_number = " << nspin
             << " real_point_number = " << nrxx;
    }

	// allocate memory
	rho = new double*[nspin];
	rhog = new std::complex<double>*[nspin];
	rho_save = new double*[nspin];
	rhog_save = new std::complex<double>*[nspin];
	if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
	{
		kin_r = new double*[nspin];
		kin_r_save = new double*[nspin];
	}

	for(int is=0; is<nspin; is++)
	{
		rho[is] = new double[nrxx];
		rhog[is] = new std::complex<double>[ngmc];
		rho_save[is] = new double[nrxx];
		rhog_save[is] = new std::complex<double>[ngmc];			
		ModuleBase::GlobalFunc::ZEROS(rho[is], nrxx);
		ModuleBase::GlobalFunc::ZEROS(rhog[is], ngmc);
		ModuleBase::GlobalFunc::ZEROS(rho_save[is], nrxx);
		ModuleBase::GlobalFunc::ZEROS(rhog_save[is], ngmc);
		if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
		{
			kin_r[is] = new double[nrxx];
			ModuleBase::GlobalFunc::ZEROS(kin_r[is], nrxx);
			kin_r_save[is] = new double[nrxx];
			ModuleBase::GlobalFunc::ZEROS(kin_r_save[is], nrxx);
		}
	}

    ModuleBase::Memory::record("Chg::rho", sizeof(double) * nspin*nrxx);
    ModuleBase::Memory::record("Chg::rho_save", sizeof(double) * nspin*nrxx);
    ModuleBase::Memory::record("Chg::rhog", sizeof(double) * nspin*ngmc);
    ModuleBase::Memory::record("Chg::rhog_save", sizeof(double) * nspin*ngmc);
	if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
	{
		ModuleBase::Memory::record("Chg::kin_r", sizeof(double) * nspin*ngmc);
		ModuleBase::Memory::record("Chg::kin_r_save", sizeof(double) * nspin*ngmc);
	}

    this->rho_core = new double[nrxx]; // core charge in real space
    ModuleBase::GlobalFunc::ZEROS( rho_core, nrxx);

	this->rhog_core = new std::complex<double>[ngmc]; // reciprocal core charge
	ModuleBase::GlobalFunc::ZEROS( rhog_core, ngmc);

    ModuleBase::Memory::record("Chg::rho_core", sizeof(double) * nrxx);
    ModuleBase::Memory::record("Chg::rhog_core", sizeof(double) * ngmc);

	this->allocate_rho = true;
    return;
}

void Charge::init_rho()
{
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "init_chg", GlobalV::init_chg);

    std::cout << " START CHARGE      : " << GlobalV::init_chg << std::endl;
    if (GlobalV::init_chg == "atomic") // mohan add 2007-10-17
    {
        this->atomic_rho(GlobalV::NSPIN, GlobalC::ucell.omega, rho, GlobalC::rhopw);
    }
    else if (GlobalV::init_chg == "file")
    {
	GlobalV::ofs_running << " try to read charge from file : ";
	for (int is=0; is<GlobalV::NSPIN; ++is)
	{
		std::stringstream ssc;
		ssc << GlobalV::global_readin_dir << "SPIN" << is + 1 << "_CHG.cube";
		GlobalV::ofs_running << ssc.str() << std::endl;
		double& ef_tmp = GlobalC::en.get_ef(is,GlobalV::TWO_EFERMI);
		if (ModuleIO::read_rho(
#ifdef __MPI
			&(GlobalC::Pgrid),
#endif
			is,
			GlobalV::NSPIN,
			ssc.str(),
			this->rho[is],
			GlobalC::rhopw->nx,
			GlobalC::rhopw->ny,
			GlobalC::rhopw->nz,
			ef_tmp,
			&(GlobalC::ucell),
			this->prenspin))
		{
			GlobalV::ofs_running << " Read in the charge density: " << ssc.str() << std::endl;
		}
		else if(is > 0)
		{
			if (prenspin == 1)
			{
			    GlobalV::ofs_running << " Didn't read in the charge density but autoset it for spin " << is + 1
			                         << std::endl;
			    for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
			    {
			        this->rho[is][ir] = 0.0;
			    }
			}
			//
			else if (prenspin == 2)
			{ // read up and down , then rearrange them.
			    if (is == 1)
			    {
			        ModuleBase::WARNING_QUIT("Charge::init_rho", "Incomplete charge density file!");
			    }
			    else if (is == 2)
			    {
			        GlobalV::ofs_running << " Didn't read in the charge density but would rearrange it later. "
			                             << std::endl;
			    }
			    else if (is == 3)
			    {
			        GlobalV::ofs_running << " rearrange charge density " << std::endl;
			        for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
			        {
			            this->rho[3][ir] = this->rho[0][ir] - this->rho[1][ir];
			            this->rho[0][ir] = this->rho[0][ir] + this->rho[1][ir];
			            this->rho[1][ir] = 0.0;
			            this->rho[2][ir] = 0.0;
			        }
			    }
			}
		}
		else
		{
			ModuleBase::WARNING_QUIT("init_rho","!!! Couldn't find the charge file !!! The default directory \n of SPIN1_CHG.cube is OUT.suffix, or you must set read_file_dir \n to a specific directory. ");
		}
	}


        
	if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
	{
			for (int is = 0; is < GlobalV::NSPIN; is++)
			{
				std::stringstream ssc;
				ssc << GlobalV::global_readin_dir << "SPIN" << is + 1 << "_TAU.cube";
				GlobalV::ofs_running << " try to read kinetic energy density from file : " << ssc.str() << std::endl;
				// mohan update 2012-02-10, sunliang update 2023-03-09
				if (ModuleIO::read_rho(
#ifdef __MPI
							&(GlobalC::Pgrid),
#endif
							is,
							GlobalV::NSPIN,
							ssc.str(),
							this->kin_r[is],
							GlobalC::rhopw->nx,
							GlobalC::rhopw->ny,
							GlobalC::rhopw->nz,
							GlobalC::en.ef,
							&(GlobalC::ucell),
							this->prenspin))
				{
					GlobalV::ofs_running << " Read in the kinetic energy density: " << ssc.str() << std::endl;
				}
			}
		}
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge::init_rho", "init_chg is wrong!");
    }

    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_load.load_charge && !GlobalC::restart.info_load.load_charge_finish)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            GlobalC::restart.load_disk("charge", is, rho);
        }
        GlobalC::restart.info_load.load_charge_finish = true;
    }
}


double Charge::sum_rho(void) const
{
	ModuleBase::TITLE("Charge","sum_rho");

    double sum_rho = 0.0;
    int nspin0 = 1;
    if(nspin==2) 
	{
		nspin0 = 2;
	}
	
	for(int is=0; is<nspin0; is++)
	{
		for(int ir=0; ir<nrxx; ir++) 
		{
			sum_rho += this->rho[is][ir];
		}
	}

	// multiply the sum of charge density by a factor
    sum_rho *= GlobalC::ucell.omega / static_cast<double>( GlobalC::rhopw->nxyz );
    Parallel_Reduce::reduce_double_pool( sum_rho );

	// mohan fixed bug 2010-01-18, 
	// sum_rho may be smaller than 1, like Na bcc.
    if (sum_rho <= 0.1)
    {
		GlobalV::ofs_warning << " sum_rho=" << sum_rho << std::endl;
        ModuleBase::WARNING_QUIT("Charge::renormalize_rho","Can't find even an electron!");
    }

    return sum_rho;
}


void Charge::renormalize_rho(void)
{
    ModuleBase::TITLE("Charge","renormalize_rho");

    const double sr = this->sum_rho();
	GlobalV::ofs_warning << std::setprecision(15);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"charge before normalized",sr);
    const double normalize_factor = GlobalV::nelec / sr;

	for(int is=0; is<nspin; is++)
	{
    	for(int ir=0; ir<nrxx; ir++)
		{
			rho[is][ir] *= normalize_factor;
		}
	}

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"charge after normalized",this->sum_rho());

	GlobalV::ofs_running << std::setprecision(6);
    return;
}

//-------------------------------------------------------
// superposition of atomic charges contained in the array
// rho_at (read from pseudopotential files)
// allocate work space (psic must already be allocated)
//-------------------------------------------------------
void Charge::atomic_rho(const int spin_number_need, const double& omega, double** rho_in, ModulePW::PW_Basis* rho_basis)const		// Peize Lin refactor 2021.04.08
{
    ModuleBase::TITLE("Charge","atomic_rho");
    ModuleBase::timer::tick("Charge","atomic_rho");

	ModuleBase::ComplexMatrix rho_g3d = [&]()->ModuleBase::ComplexMatrix
	{
		// use interpolation to get three dimension charge density.
		ModuleBase::ComplexMatrix rho_g3d( spin_number_need, rho_basis->npw);

		for (int it = 0;it < GlobalC::ucell.ntype;it++)
		{
			// check the start magnetization
			const int startmag_type = [&]()->int
			{
				if( GlobalC::ucell.magnet.start_magnetization[it] != 0.0) return 1;
				return 2;
			}();
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"startmag_type",startmag_type);

			const Atom* const atom = &GlobalC::ucell.atoms[it];

			if(!atom->flag_empty_element)		// Peize Lin add for bsse 2021.04.07
			{		
				const std::vector<double> rho_lgl = [&]()->std::vector<double>
				{
					// one dimension of charge in G space.
					std::vector<double> rho_lgl(rho_basis->ngg,0);

					// mesh point of this element.
					const int mesh = atom->ncpp.msh;

					//----------------------------------------------------------
					// Here we check the electron number 
					//----------------------------------------------------------
					const std::vector<double> rhoatm = [&]()->std::vector<double>
					{
						std::vector<double> rhoatm(mesh);		
						for(int ir=0; ir<mesh; ++ir)
						{
							double r2=atom->ncpp.r[ir]*atom->ncpp.r[ir];
							rhoatm[ir]=atom->ncpp.rho_at[ir]/ModuleBase::FOUR_PI/r2;
						}
						rhoatm[0] = pow( (rhoatm[2]/rhoatm[1]), 1./(atom->ncpp.r[2]-atom->ncpp.r[1]) );//zws add
						rhoatm[0] = pow(rhoatm[0], atom->ncpp.r[1]);
						rhoatm[0] = rhoatm[1] / rhoatm[0];

						double charge = 0.0;
						ModuleBase::Integral::Simpson_Integral(atom->ncpp.msh,atom->ncpp.rho_at,atom->ncpp.rab,charge);
						ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"charge from rho_at",charge);
						assert(charge!=0.0 || charge==atom->ncpp.zv);		// Peize Lin add charge==atom->zv for bsse 2021.04.07

						double scale=1.0;
						if(charge!=atom->ncpp.zv)
						{
							ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"charge should be",atom->ncpp.zv);
							scale = atom->ncpp.zv/charge;
						}

						for(int ir=0; ir<mesh; ++ir)
						{
							rhoatm[ir] *= scale;
							rhoatm[ir] *= (ModuleBase::FOUR_PI*atom->ncpp.r[ir]*atom->ncpp.r[ir]);
						}
						return rhoatm;
					}();

					assert(GlobalC::ucell.meshx>0);
					//----------------------------------------------------------
					// Here we compute the G=0 term
					//----------------------------------------------------------
					int gstart = 0;
					if(rho_basis->gg_uniq[0] < 1e-8)
					{
						std::vector<double> rho1d(GlobalC::ucell.meshx);
						for (int ir = 0;ir < mesh;ir++)
						{
			//              rho1d [ir] = atom->rho_at[ir];
							rho1d[ir] = rhoatm[ir];
						}
						ModuleBase::Integral::Simpson_Integral(mesh, rho1d.data(), atom->ncpp.rab, rho_lgl[0]);
						gstart = 1;
					}
					if (GlobalV::test_charge>0) std::cout<<"\n |G|=0 term done." <<std::endl;
					//----------------------------------------------------------
					// Here we compute the G<>0 term
					// But if in parallel case
					// G=0 term only belong to 1 cpu.
					// Other processors start from '0'
					//----------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel
{
#endif
					std::vector<double> rho1d(GlobalC::ucell.meshx);

#ifdef _OPENMP
#pragma omp for
#endif
					for (int igg = gstart; igg < rho_basis->ngg ;++igg)
					{
						const double gx = sqrt(rho_basis->gg_uniq[igg]) * GlobalC::ucell.tpiba;
						for (int ir = 0; ir < mesh;ir++)
						{
							if ( atom->ncpp.r[ir] < 1.0e-8 )
							{
								rho1d[ir] = rhoatm[ir];
								//rho1d[ir] = atom->rho_at[ir];
							}
							else
							{
								const double gxx = gx * atom->ncpp.r[ir];
								rho1d[ir] = rhoatm[ir] * ModuleBase::libm::sin(gxx) / gxx;
							}
						}
						ModuleBase::Integral::Simpson_Integral(mesh, rho1d.data(), atom->ncpp.rab, rho_lgl[igg]);
					}
#ifdef _OPENMP
#pragma omp single
#endif
					{ if (GlobalV::test_charge>0) std::cout<<" |G|>0 term done." <<std::endl; }
					//----------------------------------------------------------
					// EXPLAIN : Complete the transfer of rho from real space to
					// reciprocal space
					//----------------------------------------------------------
#ifdef _OPENMP
#pragma omp for
#endif
					for (int igg=0; igg< rho_basis->ngg ; igg++)
						rho_lgl[igg] /= omega;
#ifdef _OPENMP
}
#endif
					return rho_lgl;
				}();
				//----------------------------------------------------------
				// EXPLAIN : compute the 3D atomic charge in reciprocal space
				//----------------------------------------------------------
				if(spin_number_need==1)
				{
#ifdef _OPENMP
#pragma omp parallel for
#endif
					for (int ig=0; ig< rho_basis->npw ;ig++)
					{
						rho_g3d(0, ig) += GlobalC::sf.strucFac(it, ig) * rho_lgl[ rho_basis->ig2igg[ig] ];
					}
				}
				// mohan add 2011-06-14, initialize the charge density according to each atom 
				else if(spin_number_need==2)
				{
					if(startmag_type==1)
					{
#ifdef _OPENMP
#pragma omp parallel for
#endif
						for (int ig = 0; ig < rho_basis->npw ; ig++)
						{
							const std::complex<double> swap = GlobalC::sf.strucFac(it, ig)* rho_lgl[rho_basis->ig2igg[ig]];
							const double up = 0.5 * ( 1 + GlobalC::ucell.magnet.start_magnetization[it] / atom->ncpp.zv );
							const double dw = 0.5 * ( 1 - GlobalC::ucell.magnet.start_magnetization[it] / atom->ncpp.zv );
							rho_g3d(0, ig) += swap * up;
							rho_g3d(1, ig) += swap * dw;
						}
					}
					// mohan add 2011-06-14
					else if(startmag_type==2)
					{
						std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;
						for (int ia = 0; ia < atom->na; ia++)
						{
							//const double up = 0.5 * ( 1 + atom->mag[ia] );
							//const double dw = 0.5 * ( 1 - atom->mag[ia] );
							const double up = 0.5 * ( 1 + atom->mag[ia] / atom->ncpp.zv );
							const double dw = 0.5 * ( 1 - atom->mag[ia] / atom->ncpp.zv );
							//std::cout << " atom " << ia << " up=" << up << " dw=" << dw << std::endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
							for (int ig = 0; ig < rho_basis->npw ; ig++)
							{
								const double Gtau =
									rho_basis->gcar[ig][0] * atom->tau[ia].x + 
									rho_basis->gcar[ig][1] * atom->tau[ia].y + 
									rho_basis->gcar[ig][2] * atom->tau[ia].z;

								std::complex<double> swap = ModuleBase::libm::exp(ci_tpi * Gtau) * rho_lgl[rho_basis->ig2igg[ig]];

								rho_g3d(0, ig) += swap * up;
								rho_g3d(1, ig) += swap * dw;
							}
						}
					}
				}
				else if(spin_number_need==4)
				{
					//noncolinear case
					if(startmag_type == 1)
					{
						double sin_a1, sin_a2, cos_a1, cos_a2;
						if(GlobalV::DOMAG)
						{//will not be used now, will be deleted later
							ModuleBase::libm::sincos(atom->angle1[0], &sin_a1, &cos_a1);
							ModuleBase::libm::sincos(atom->angle2[0], &sin_a2, &cos_a2);
						}
#ifdef _OPENMP
#pragma omp parallel for
#endif
						for (int ig = 0; ig < rho_basis->npw ; ig++)
						{
							const std::complex<double> swap = GlobalC::sf.strucFac(it, ig)* rho_lgl[rho_basis->ig2igg[ig]];
							rho_g3d(0, ig) += swap ;
							if(GlobalV::DOMAG)
							{//will not be used now, will be deleted later
								rho_g3d(1, ig) += swap * (GlobalC::ucell.magnet.start_magnetization[it] / atom->ncpp.zv) 
								* sin_a1 * cos_a2;
								rho_g3d(2, ig) += swap * (GlobalC::ucell.magnet.start_magnetization[it] / atom->ncpp.zv) 
								* sin_a1 * sin_a2;
								rho_g3d(3, ig) += swap * (GlobalC::ucell.magnet.start_magnetization[it] / atom->ncpp.zv) 
								* cos_a1;
							}
							else if(GlobalV::DOMAG_Z)
							{
								rho_g3d(1, ig) = 0.0;
								rho_g3d(2, ig) = 0.0;
								rho_g3d(3, ig) += swap * (GlobalC::ucell.magnet.start_magnetization[it] / atom->ncpp.zv);
							}
						}
					}
					else if(startmag_type == 2)
					{//zdy-warning-not-available
						std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;
						for(int ia = 0;ia<atom->na;ia++)
						{
							double sin_a1, sin_a2, cos_a1, cos_a2;
							if(GlobalV::DOMAG || GlobalV::DOMAG_Z)
							{
								ModuleBase::libm::sincos(atom->angle1[ia], &sin_a1, &cos_a1);
							}
							if(GlobalV::DOMAG)
							{
								ModuleBase::libm::sincos(atom->angle2[ia], &sin_a2, &cos_a2);
							}
#ifdef _OPENMP
#pragma omp parallel for
#endif
							for (int ig = 0; ig < rho_basis->npw ; ig++)
							{
								const double Gtau =
									rho_basis->gcar[ig][0] * atom->tau[ia].x + 
									rho_basis->gcar[ig][1] * atom->tau[ia].y + 
									rho_basis->gcar[ig][2] * atom->tau[ia].z;

								std::complex<double> swap = exp(ci_tpi * Gtau) * rho_lgl[rho_basis->ig2igg[ig]];

								//calculate rho_total
								rho_g3d(0, ig) += swap;
								//calculate mag_z
								if(GlobalV::DOMAG || GlobalV::DOMAG_Z)
								{
									rho_g3d(3, ig) += swap * (atom->mag[ia] / atom->ncpp.zv) 
										* cos_a1;
								}
								//calculate mag_x and mag_y
								if(GlobalV::DOMAG)
								{
									rho_g3d(1, ig) += swap * (atom->mag[ia] / atom->ncpp.zv) 
										* sin_a1 * cos_a2;
									rho_g3d(2, ig) += swap * (atom->mag[ia] / atom->ncpp.zv) 
										* sin_a1 * sin_a2;
								}
								else
								{
									rho_g3d(1, ig) = 0.0;
									rho_g3d(2, ig) = 0.0;
								}
							}
						}
					}
				}
				else
				{
					ModuleBase::WARNING_QUIT("Charge::spin_number_need"," Either 1 or 2 or 4, check SPIN number !");
				}
			}
		}
		return rho_g3d;
	}();

	assert( spin_number_need > 0 );
	std::vector<double> ne(spin_number_need);
    for (int is = 0; is < spin_number_need;is++)
    {
        rho_basis->recip2real( &rho_g3d(is,0), rho_in[is]);

		for(int ir=0; ir<rho_basis->nrxx; ++ir)
			ne[is] += rho_in[is][ir];
		ne[is] *= omega/(double)rho_basis->nxyz; 
		Parallel_Reduce::reduce_double_pool( ne[is] );

        // we check that everything is correct
        double neg = 0.0;
        double rea = 0.0;
        double ima = 0.0;
		double sumrea = 0.0;
        for (int ir=0;ir < rho_basis->nrxx; ir++)
        {
            rea = rho_basis->ft.get_auxr_data<double>()[ir].real();
			sumrea += rea;
            neg += std::min(0.0, rea);
            ima += abs(rho_basis->ft.get_auxr_data<double>()[ir].imag());
        }

		Parallel_Reduce::reduce_double_pool( neg );	
		Parallel_Reduce::reduce_double_pool( ima );	
		Parallel_Reduce::reduce_double_pool( sumrea );	

		// mohan fix bug 2011-04-03
        neg = neg / (double)rho_basis->nxyz * omega;
        ima = ima / (double)rho_basis->nxyz * omega;
		sumrea = sumrea / (double)rho_basis->nxyz * omega;

        if( ((neg<-1.0e-4) && (is==0||GlobalV::NSPIN==2)) || ima>1.0e-4)
        {
            GlobalV::ofs_warning << " Warning: negative or imaginary starting charge : " ;
            GlobalV::ofs_warning << " neg = " << neg
                 << " ima = " << ima
                 << " SPIN = " << is << std::endl;
        }

//		std::cout << " sum rho for spin " << is << " = " << sumrea << std::endl;
//		std::cout << " sum rho for spin " << is << " = " << sumrea << std::endl;

    }//end is


	double ne_tot = 0.0;
	int spin0=1;
	if(spin_number_need == 2) spin0 = spin_number_need;
	for(int is=0; is<spin0; ++is)
	{
		GlobalV::ofs_warning << "\n SETUP ATOMIC RHO FOR SPIN " << is+1 << std::endl;
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"Electron number from rho",ne[is]);
		ne_tot += ne[is];
	}
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"total electron number from rho",ne_tot);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"should be",GlobalV::nelec);
	for(int is=0; is<spin_number_need; ++is)
		for(int ir=0; ir<rho_basis->nrxx; ++ir)
			rho_in[is][ir] = rho_in[is][ir] / ne_tot * GlobalV::nelec;

	//wenfei 2021-7-29 : initial tau = 3/5 rho^2/3, Thomas-Fermi
	if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
	{
		const double pi = 3.141592653589790;
		double fact = (3.0/5.0)*pow(3.0*pi*pi,2.0/3.0);
		int nspin = spin_number_need;
		//ofstream test_tau0("tau0");
		for(int is=0; is<spin_number_need; ++is)
			for(int ir=0; ir<rho_basis->nrxx; ++ir)
			{
				kin_r[is][ir] = fact * pow(abs(rho_in[is][ir])*nspin,5.0/3.0)/nspin;
				//test_tau0 << rho_in[is][ir] << " " << kin_r[is][ir] << endl;
			}
	}

    ModuleBase::timer::tick("Charge","atomic_rho");
    return;
}


//==========================================================
// computes the core charge on the real space 3D mesh.
//==========================================================
void Charge::set_rho_core(
    const ModuleBase::ComplexMatrix &structure_factor
)
{
    ModuleBase::TITLE("Charge","set_rho_core");
    ModuleBase::timer::tick("Charge","set_rho_core");

    //double eps = 1.e-10;
    GlobalC::en.etxcc = 0.0;
//----------------------------------------------------------
// LOCAL VARIABLES :
// counter on mesh points
// counter on atomic types
// counter on g vectors
//----------------------------------------------------------
    //int ir = 0;
    //int it = 0;
    //int ig = 0;

    bool bl = false;
    for (int it = 0; it<GlobalC::ucell.ntype; it++)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {
            bl = true;
            break;
        }
    }

    if (!bl)
    {
        ModuleBase::GlobalFunc::ZEROS( this->rho_core, GlobalC::rhopw->nrxx);
    	ModuleBase::timer::tick("Charge","set_rho_core");
        return;
    }

    double *rhocg = new double[GlobalC::rhopw->ngg];
    ModuleBase::GlobalFunc::ZEROS(rhocg, GlobalC::rhopw->ngg );

	// three dimension.
    std::complex<double> *vg = new std::complex<double>[GlobalC::rhopw->npw];	

    for (int it = 0; it < GlobalC::ucell.ntype;it++)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {
//----------------------------------------------------------
// EXPLAIN : drhoc compute the radial fourier transform for
// each shell of g vec
//----------------------------------------------------------
            this->non_linear_core_correction(
                GlobalC::ppcell.numeric,
                GlobalC::ucell.atoms[it].ncpp.msh,
                GlobalC::ucell.atoms[it].ncpp.r,
                GlobalC::ucell.atoms[it].ncpp.rab,
                GlobalC::ucell.atoms[it].ncpp.rho_atc,
                rhocg,
				GlobalC::rhopw);
//----------------------------------------------------------
// EXPLAIN : multiply by the structure factor and sum
//----------------------------------------------------------
            for (int ig = 0; ig < GlobalC::rhopw->npw ; ig++)
            {
                vg[ig] += structure_factor(it, ig) * rhocg[GlobalC::rhopw->ig2igg[ig]];
            }
        }
    }

	// for tmp use.
	for(int ig=0; ig< GlobalC::rhopw->npw; ig++)
	{
		this->rhog_core[ig] = vg[ig];
	}

    GlobalC::rhopw->recip2real(vg, this->rho_core);

    // test on the charge and computation of the core energy
    double rhoima = 0.0;
    double rhoneg = 0.0;
    for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
    {
        rhoneg += min(0.0, GlobalC::rhopw->ft.get_auxr_data<double>()[ir].real());
        rhoima += abs(GlobalC::rhopw->ft.get_auxr_data<double>()[ir].imag());
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
    rhoneg /= GlobalC::rhopw->nxyz * GlobalC::ucell.omega;
    rhoima /= GlobalC::rhopw->nxyz * GlobalC::ucell.omega;

    // calculate core_only exch-corr energy etxcc=E_xc[rho_core] if required
    // The term was present in previous versions of the code but it shouldn't
    delete [] rhocg;
    delete [] vg;
    ModuleBase::timer::tick("Charge","set_rho_core");
    return;
} // end subroutine set_rhoc



void Charge::non_linear_core_correction
(
    const bool &numeric,
    const int mesh,
    const double *r,
    const double *rab,
    const double *rhoc,
    double *rhocg,
	ModulePW::PW_Basis* rho_basis) const
{
    ModuleBase::TITLE("charge","drhoc");

	// use labmda instead of repeating codes
	const auto kernel = [&](int num_threads, int thread_id)
	{

	double gx = 0.0;
    double rhocg1 = 0.0;
    double *aux;

    // here we compute the fourier transform is the charge in numeric form
    if (numeric)
    {
        aux = new double [mesh];
        // G=0 term

        int igl0 = 0;
        if (rho_basis->gg_uniq [0] < 1.0e-8)
        {
			// single thread term
			if (thread_id == 0)
			{
				for (int ir = 0;ir < mesh; ir++)
				{
					aux [ir] = r [ir] * r [ir] * rhoc [ir];
				}
				ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
				//rhocg [1] = fpi * rhocg1 / omega;
				rhocg [0] = ModuleBase::FOUR_PI * rhocg1 / GlobalC::ucell.omega;//mohan modify 2008-01-19
			}
            igl0 = 1;
        }

		int igl_beg, igl_end;
		// exclude igl0
		ModuleBase::TASK_DIST_1D(num_threads, thread_id, rho_basis->ngg - igl0, igl_beg, igl_end);
		igl_beg += igl0;
		igl_end += igl_beg;

        // G <> 0 term
        for (int igl = igl_beg; igl < igl_end;igl++) 
        {
            gx = sqrt(rho_basis->gg_uniq[igl] * GlobalC::ucell.tpiba2);
            ModuleBase::Sphbes::Spherical_Bessel(mesh, r, gx, 0, aux);
            for (int ir = 0;ir < mesh; ir++) 
            {
                aux [ir] = r[ir] * r[ir] * rhoc [ir] * aux [ir];
            } //  enddo
            ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
            rhocg [igl] = ModuleBase::FOUR_PI * rhocg1 / GlobalC::ucell.omega;
        } //  enddo
        delete [] aux;
    }
    else
    {
        // here the case where the charge is in analytic form,
        // check old version before 2008-12-9
    }

	}; // end kernel

	// do not use omp parallel when this function is already in parallel block
	// 
	// it is called in parallel block in Forces::cal_force_cc,
	// but not in other funtcion such as Stress_Func::stress_cc.
	ModuleBase::TRY_OMP_PARALLEL(kernel);

    return;
}


#ifdef __MPI
void Charge::rho_mpi(void)
{
	ModuleBase::TITLE("Charge","rho_mpi");
    if (GlobalV::NPROC==1) return;
	if(GlobalV::ESOLVER_TYPE == "sdft" && GlobalV::NPROC_IN_STOGROUP==1) 
		return;//qinarui add it temporarily.
    ModuleBase::timer::tick("Charge","rho_mpi");
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
    int *num_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
    for (iz=0;iz<GlobalC::bigpw->nbz;iz++)
    {
        ip = iz % GlobalV::NPROC_IN_POOL;
        num_z[ip]++;
    }

	// mohan update 2011-04-26
	for(int ip=0; ip<GlobalV::NPROC_IN_POOL; ip++)
	{
		num_z[ip]*=GlobalC::bigpw->bz;
	}

    //=======================================
    // Find current number of planes (nz)
	// start_z: start position of z in 
	// processor ip.
    //=======================================
    int *start_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(start_z, GlobalV::NPROC_IN_POOL);
    for (ip=1;ip<GlobalV::NPROC_IN_POOL;ip++)
    {
        start_z[ip] = start_z[ip-1]+num_z[ip-1];
    }

    //====================================================
    // Find "number of data" in each processor in each pool
    //====================================================
    int *rec = new int[GlobalV::NPROC_IN_POOL];
	ModuleBase::GlobalFunc::ZEROS(rec, GlobalV::NPROC_IN_POOL);
    const int ncxy = GlobalC::rhopw->nx * GlobalC::rhopw->ny;
    for (ip=0;ip<GlobalV::NPROC_IN_POOL;ip++)
    {
        rec[ip] = num_z[ip]*ncxy;
    }

    //======================================================
    // Find current "index of data" in each cpu in this pool
	// also, we mean start position of data.
    //======================================================
    int *dis = new int[GlobalV::NPROC_IN_POOL];
	ModuleBase::GlobalFunc::ZEROS(dis, GlobalV::NPROC_IN_POOL);
    for (ip=1;ip<GlobalV::NPROC_IN_POOL;ip++)
    {
        dis[ip]=dis[ip-1]+rec[ip-1];
    }

    //==========================================
    // Collection of rho in each pool
    // ( according to different k distribution,
    // so the rho in each pool is different
    //==========================================
    double *rho_tmp = new double[GlobalC::rhopw->nrxx];
    double *rho_tot = new double[GlobalC::rhopw->nxyz];
    double *rho_tot_aux = new double[GlobalC::rhopw->nxyz];
	ModuleBase::GlobalFunc::ZEROS(rho_tot_aux, GlobalC::rhopw->nxyz);

	double *tau_tmp;
	double *tau_tot;
	double *tau_tot_aux;

	if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
	{
    	tau_tmp = new double[GlobalC::rhopw->nrxx];
	    tau_tot = new double[GlobalC::rhopw->nxyz];
    	tau_tot_aux = new double[GlobalC::rhopw->nxyz];
		ModuleBase::GlobalFunc::ZEROS(tau_tot_aux, GlobalC::rhopw->nxyz);
	}

    for (int is=0; is< GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(rho_tot, GlobalC::rhopw->nxyz);
		if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) ModuleBase::GlobalFunc::ZEROS(tau_tot, GlobalC::rhopw->nxyz);

		for (ir=0;ir<GlobalC::rhopw->nrxx;ir++)
		{
			rho_tmp[ir] = this->rho[is][ir] / static_cast<double>(GlobalV::NPROC_IN_POOL);
			if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
			{
				tau_tmp[ir] = this->kin_r[is][ir] / static_cast<double>(GlobalV::NPROC_IN_POOL);
			}
		}

        MPI_Allgatherv(rho_tmp, GlobalC::rhopw->nrxx, MPI_DOUBLE, rho_tot, rec, dis, MPI_DOUBLE, POOL_WORLD);
		if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
		{
        	MPI_Allgatherv(tau_tmp, GlobalC::rhopw->nrxx, MPI_DOUBLE, tau_tot, rec, dis, MPI_DOUBLE, POOL_WORLD);
		}
        //=================================================================
        // Change the order of rho_tot in each pool , make them consistent
        // this is the most complicated part !!
        //=================================================================
        ModuleBase::GlobalFunc::ZEROS(rho_tot_aux, GlobalC::rhopw->nxyz);
		if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
		{
        	ModuleBase::GlobalFunc::ZEROS(tau_tot_aux, GlobalC::rhopw->nxyz);
		}

        for (ip=0;ip<GlobalV::NPROC_IN_POOL;ip++)
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
					// is 'GlobalC::rhopw->nz'
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
                    rho_tot_aux[GlobalC::rhopw->nz*ir    + start_z[ip]      + iz]
                      = rho_tot[num_z[ip]*ir + start_z[ip]*ncxy + iz];
					if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
					{
                    	tau_tot_aux[GlobalC::rhopw->nz*ir    + start_z[ip]      + iz]
                      	  = tau_tot[num_z[ip]*ir + start_z[ip]*ncxy + iz];
					}	
                }
            }
        }
        //==================================
        // Reduce all the rho in each cpu
        //==================================
		if(GlobalV::ESOLVER_TYPE == "sdft") //qinarui add it temporarily.
		{
			MPI_Allreduce(rho_tot_aux,rho_tot,GlobalC::rhopw->nxyz,MPI_DOUBLE,MPI_SUM,STO_WORLD);
			if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
			{
				MPI_Allreduce(tau_tot_aux,tau_tot,GlobalC::rhopw->nxyz,MPI_DOUBLE,MPI_SUM,STO_WORLD);
			}
		}
		else
        MPI_Allreduce(rho_tot_aux,rho_tot,GlobalC::rhopw->nxyz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
		{
   	    	MPI_Allreduce(tau_tot_aux,tau_tot,GlobalC::rhopw->nxyz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		}
        
		//=====================================
        // Change the order of rho in each cpu
        //=====================================
        for (ir =0;ir<ncxy;ir++)
        {
            for (iz=0;iz<num_z[GlobalV::RANK_IN_POOL];iz++)
            {
                this->rho[is][num_z[GlobalV::RANK_IN_POOL]*ir+iz] = rho_tot[GlobalC::rhopw->nz*ir + GlobalC::rhopw->startz_current + iz ];
				if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
				{
                	this->kin_r[is][num_z[GlobalV::RANK_IN_POOL]*ir+iz] = tau_tot[GlobalC::rhopw->nz*ir + GlobalC::rhopw->startz_current + iz ];
				}	
            }
        }
    }
    delete[] rho_tot_aux;
    delete[] rho_tot;
    delete[] rho_tmp;
    
	if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
	{
		delete[] tau_tot_aux;
    	delete[] tau_tot;
	    delete[] tau_tmp;
	}
    delete[] rec;
    delete[] dis;

    delete[] num_z;
    delete[] start_z;
    ModuleBase::timer::tick("Charge","rho_mpi");
    return;
}
#endif


void Charge::save_rho_before_sum_band(void)
{
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
    	ModuleBase::GlobalFunc::DCOPY( rho[is], rho_save[is], GlobalC::rhopw->nrxx);
    	if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) ModuleBase::GlobalFunc::DCOPY( kin_r[is], kin_r_save[is], GlobalC::rhopw->nrxx);
    }
    return;
}


double Charge::check_ne(const double *rho_in) const 
{
	double ne= 0.0;
	for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	{
		ne += rho_in[ir];
	}
	Parallel_Reduce::reduce_double_pool( ne );
	ne = ne * GlobalC::ucell.omega / (double)GlobalC::rhopw->nxyz;
	std::cout << std::setprecision(10);
	std::cout << " check the electrons number from rho, ne =" << ne << std::endl;
	std::cout << std::setprecision(6);
	return ne;
}


//LiuXh add 20180619
void Charge::init_final_scf()
{
    ModuleBase::TITLE("Charge","init_after_scf");

	assert(allocate_rho_final_scf == false);

    if (GlobalV::test_charge > 1)
    {
        std::cout << "\n spin_number = " << GlobalV::NSPIN
             << " real_point_number = " << GlobalC::rhopw->nrxx;
    }

	// allocate memory
	rho = new double*[GlobalV::NSPIN];
	rhog = new std::complex<double>*[GlobalV::NSPIN];
	rho_save = new double*[GlobalV::NSPIN];
	rhog_save = new std::complex<double>*[GlobalV::NSPIN];

	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		rho[is] = new double[GlobalC::rhopw->nrxx];
		rhog[is] = new std::complex<double>[GlobalC::rhopw->npw];
		rho_save[is] = new double[GlobalC::rhopw->nrxx];
		rhog_save[is] = new std::complex<double>[GlobalC::rhopw->npw];			
		ModuleBase::GlobalFunc::ZEROS(rho[is], GlobalC::rhopw->nrxx);
		ModuleBase::GlobalFunc::ZEROS(rhog[is], GlobalC::rhopw->npw);
		ModuleBase::GlobalFunc::ZEROS(rho_save[is], GlobalC::rhopw->nrxx);
		ModuleBase::GlobalFunc::ZEROS(rhog_save[is], GlobalC::rhopw->npw);
	}

    ModuleBase::Memory::record("Chg::rho", sizeof(double) * GlobalV::NSPIN*GlobalC::rhopw->nrxx);
    ModuleBase::Memory::record("Chg::rho_save", sizeof(double) * GlobalV::NSPIN*GlobalC::rhopw->nrxx);
    ModuleBase::Memory::record("Chg::rhog", sizeof(double) * GlobalV::NSPIN*GlobalC::rhopw->npw);
    ModuleBase::Memory::record("Chg::rhog_save", sizeof(double) * GlobalV::NSPIN*GlobalC::rhopw->npw);

    this->rho_core = new double[GlobalC::rhopw->nrxx]; // core charge in real space
    ModuleBase::GlobalFunc::ZEROS( rho_core, GlobalC::rhopw->nrxx);

	this->rhog_core = new std::complex<double>[GlobalC::rhopw->npw]; // reciprocal core charge
	ModuleBase::GlobalFunc::ZEROS( rhog_core, GlobalC::rhopw->npw);

    ModuleBase::Memory::record("Chg::rho_core", sizeof(double) * GlobalC::rhopw->nrxx);
    ModuleBase::Memory::record("Chg::rhog_core", sizeof(double) * GlobalC::rhopw->npw);

	this->allocate_rho_final_scf = true;
    return;
}

//=========================================================
// calculate total number of electrons (GlobalV::nelec) and default
// number of bands (GlobalV::NBANDS).
//=========================================================
#include "module_elecstate/occupy.h"
void Charge::cal_nelec(void)
{
	ModuleBase::TITLE("UnitCell","cal_nelec");
	//=======================================================
	// calculate the total number of electrons in the system
	// if GlobalV::nelec <>0; use input number (setup.f90)
	//=======================================================

	GlobalV::ofs_running << "\n SETUP THE ELECTRONS NUMBER" << std::endl;

	if (GlobalV::nelec == 0)
	{
		for (int it = 0; it < GlobalC::ucell.ntype;it++)
		{
			std::stringstream ss1, ss2;
			ss1 << "electron number of element " << GlobalC::ucell.atoms[it].label;
			const int nelec_it = GlobalC::ucell.atoms[it].ncpp.zv * GlobalC::ucell.atoms[it].na;
			GlobalV::nelec += nelec_it;
			ss2 << "total electron number of element " << GlobalC::ucell.atoms[it].label; 
			
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,ss1.str(),GlobalC::ucell.atoms[it].ncpp.zv);
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,ss2.str(),nelec_it);
		}
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "AUTOSET number of electrons: ", GlobalV::nelec);
	}
	return;
}
