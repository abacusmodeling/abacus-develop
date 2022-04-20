#include "charge_pulay.h"
#include "global.h"
#include "../module_base/inverse_matrix.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/memory.h"

Charge_Pulay::Charge_Pulay()
{
    rstep = 0;
	dstep = rstep - 1;//alway like this.
    initp = false;
}
Charge_Pulay::~Charge_Pulay()
{
    if (initp)
    {
		// delete: Rrho[i] = rho_out[i] - rho_in[i];
        for (int is=0; is<GlobalV::NSPIN; is++)
        {
            for (int i=0; i<rstep; i++)
            {
                delete[] Rrho[is][i];
            }
            delete[] Rrho[is];
        }
        delete[] Rrho;

		// delete: dRrho[i] = Rrho[i+1] - Rrho[i]
       	for (int is=0; is<GlobalV::NSPIN; is++)
		{
			for (int i=0; i<dstep; i++)
            {
                delete[] dRrho[is][i];
				delete[] drho[is][i];
            }
            delete[] dRrho[is];
			delete[] drho[is];
		}
        delete[] dRrho;
		delete[] drho;

		// dimension: dstep
		delete[] dRR;
		delete[] alpha;

		// dimension of rho_save2(GlobalV::NSPIN, GlobalC::pw.nrxx)
		for (int is=0; is<GlobalV::NSPIN; is++)
		{
			delete[] rho_save2[is];
		}
		delete[] rho_save2;	
    }
}


void Charge_Pulay::Pulay_mixing(void)
{
//  ModuleBase::TITLE("Charge_Pulay","Pulay_mixing");
	rstep = this->mixing_ndim;
	dstep = this->mixing_ndim - 1;
	assert(dstep>0);

	const int scheme = 1;
	
	// scheme 2 will only work correctly for the first
	// time to calculate residual std::vector norm,
	// which is one way to check if scheme 1 right.
	// scheme 1 is correct to provide the final charge
	// density.

	// (1) allocate 
    if(GlobalV::FINAL_SCF && totstep==0) initp = false;
	this->allocate_pulay(scheme);

	// irstep: iteration step for rstep (Rrho)
	// idstep: iteration step for dstep (dRrho)
	// totstep only used for the first few iterations.
	// At the beginning of each ion iteration, reset the three variables.
	// mohan add 2010-07-16
	if(this->new_e_iteration)
	{
		irstep = 0;
		idstep = 0;
		totstep = 0;
	}

    if (irstep==rstep) irstep=0;
	if (idstep==dstep) idstep=0;

	if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"irstep",irstep);
	if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"idstep",idstep);
	if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"totstep",totstep);

	//----------------------------------------------
	// calculate "dR^{i} = R^{i+1} - R^{i}"
	// calculate "drho^{i} = rho^{i+1} - rho^{i}"
	//----------------------------------------------
	this->generate_datas(irstep, idstep, totstep);

	// not enough steps, not full matrix.
	if(totstep < dstep) 
	{
		int premix = 1;
		
		// mohan update 2011-06-14
		if(totstep==0) premix = 2;

		if(premix == 1)
		{
			this->generate_Abar(scheme,Abar);
			
			//-----------------------------
			// inverse part of the matrix.
			//-----------------------------
			const int predim = idstep+1;
			ModuleBase::matrix preA(predim, predim);
			for(int i=0; i<predim; i++)
			{
				for(int j=0; j<predim; j++)
				{
					preA(i,j) = Abar(i,j);
				}
			}
			this->inverse_preA(predim,preA);
			
			
			//----------------------------------------------
			// get the information from part of the matrix. 
			//----------------------------------------------
			Abar.zero_out();
			for(int i=0; i<predim; i++)
			{
				for(int j=0; j<predim; j++)
				{
					Abar(i,j) = preA(i,j);
				}			
			}

			this->generate_dRR(irstep);

			this->generate_alpha(scheme);

			for(int is=0; is<GlobalV::NSPIN; is++)
			{
				this->generate_new_rho(is,irstep);
			}
		}

		if(premix == 2)
		{
			// if not enough step, take kerker mixing method.	
			for(int is=0; is<GlobalV::NSPIN; is++)
			{
				//this->Kerker_mixing( this->rho[is], this->rhog[is], this->rho_save[is] );
				this->plain_mixing( this->rho[is], this->rho_save[is]);
			}
		}

		if(totstep>0)
		{
			++idstep;
		}
		++irstep;
		++totstep;
	}
	else
	{
		// generate matrix A = <dR|dR>
		this->generate_Abar(scheme,Abar);
		ModuleBase::matrix A(Abar);

		// inverse A matrix to become Abar.
		this->inverse_real_symmetry_matrix(scheme,Abar);

		if(scheme==1)
		{
			this->generate_dRR(irstep);
		}

		this->generate_alpha(scheme);

		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			this->generate_new_rho(is,irstep);
		}

		++irstep;
		++idstep;
		++totstep;	
	}

    // output file name
    /*
    std::stringstream ofname;
    ofname << GlobalV::global_out_dir << "PULAY_MIXING_" << counter << ".dat";
    this->write_rho(ofname.str());

    for(int i=0; i<counter; i++)
    {
    	std::stringstream ifname;
    	ifname << GlobalV::global_out_dir << "PULAY_MIXING_" << counter << ".dat";
    	this->read_rho(ifname.str());
    }
    */
	#if TEST_EXX_LCAO==1
		std::cout<<"Charge_Pulay::Pulay_mixing\t"<<__FILE__<<"\t"<<__LINE__<<std::endl;
		std::cout<<"irstep:\t"<<irstep<<std::endl;
		std::cout<<"idstep:\t"<<idstep<<std::endl;
		std::cout<<"rstep:\t"<<rstep<<std::endl;
		std::cout<<"dstep:\t"<<dstep<<std::endl;
		std::cout<<"totstep:\t"<<totstep<<std::endl;
	#elif TEST_EXX_LCAO==-1
		#error
	#endif
	return;		
}

void Charge_Pulay::set_new_e_iteration( const bool new_e_iteration_in )		// Peize Lin add 2018-11-01
{
	this->new_e_iteration = new_e_iteration_in;
	
	if(this->new_e_iteration)
	{
		irstep = 0;
		idstep = 0;
		totstep = 0;
	}	
}

void Charge_Pulay::allocate_pulay(const int &scheme)
{
	//std::cout << "\n initp = " << initp << std::endl;
	if(!this->initp)
	{
		ModuleBase::TITLE("Charge_Pulay","allocate_pulay");
		ModuleBase::GlobalFunc::NOTE("rstep is used to record Rrho");
		if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"rstep",rstep);
		ModuleBase::GlobalFunc::NOTE("dstep is used to record dRrho, drho");
		if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"dstep",dstep);
		assert(rstep>1);
		dstep = rstep - 1;
		
		// (1) allocate Rrho[i]: rho_out[i] - rho_in[i]
    	this->Rrho = new double**[GlobalV::NSPIN];
    	for (int is=0; is<GlobalV::NSPIN; is++)
    	{
        	Rrho[is] = new double*[rstep];
        	for (int i=0; i<rstep; i++)
        	{
            	Rrho[is][i] = new double[GlobalC::pw.nrxx];
				ModuleBase::GlobalFunc::ZEROS( Rrho[is][i], GlobalC::pw.nrxx );
        	}
		}
    	ModuleBase::Memory::record("Charge_Pulay","Rrho", GlobalV::NSPIN*rstep*GlobalC::pw.nrxx,"double");

		// (2) allocate "dRrho[i] = Rrho[i+1] - Rrho[i]" of the last few steps.
		// allocate "drho[i] = rho[i+1] - rho[i]" of the last few steps.
		// allocate rho_save2: rho[i] of the last few steps.
		this->dRrho = new double**[GlobalV::NSPIN];
		this->drho = new double**[GlobalV::NSPIN];
		this->rho_save2 = new double*[GlobalV::NSPIN];

		for (int is=0; is<GlobalV::NSPIN; is++)
		{
			dRrho[is] = new double*[dstep];
			drho[is] = new double*[dstep];
			rho_save2[is] = new double[GlobalC::pw.nrxx];
			ModuleBase::GlobalFunc::ZEROS( rho_save2[is], GlobalC::pw.nrxx);

			for (int i=0; i<dstep; i++)
			{
				dRrho[is][i] = new double[GlobalC::pw.nrxx];	
				drho[is][i] = new double[GlobalC::pw.nrxx];
				ModuleBase::GlobalFunc::ZEROS( dRrho[is][i], GlobalC::pw.nrxx );
				ModuleBase::GlobalFunc::ZEROS( drho[is][i], GlobalC::pw.nrxx);
			}
		}
    	ModuleBase::Memory::record("Charge_Pulay","dRrho", GlobalV::NSPIN*dstep*GlobalC::pw.nrxx,"double");
    	ModuleBase::Memory::record("Charge_Pulay","drho", GlobalV::NSPIN*dstep*GlobalC::pw.nrxx,"double");
    	ModuleBase::Memory::record("Charge_Pulay","rho_save2", GlobalV::NSPIN*GlobalC::pw.nrxx,"double");

		ModuleBase::GlobalFunc::NOTE("Allocate Abar = <dRrho_j | dRrho_i >, dimension = dstep.");
		if(scheme==1)
		{
			this->Abar.create(dstep, dstep);
    		ModuleBase::Memory::record("Charge_Pulay","Abar", dstep*dstep,"double");
		}
		else if(scheme==2)
		{
			this->Abar.create(rstep, rstep);
			ModuleBase::Memory::record("Charge_Pulay","Abar", rstep*rstep,"double");
		}

		// (4) allocate dRR = <delta R|R>
		ModuleBase::GlobalFunc::NOTE("Allocate dRR = < dR | R >, dimension = dstep");
		this->dRR = new double[dstep];
		ModuleBase::GlobalFunc::ZEROS(dRR, dstep);

		// (5) allocate alpha
		ModuleBase::GlobalFunc::NOTE("Allocate alpha, dimension = dstep");
		this->alpha = new double[dstep];
		ModuleBase::GlobalFunc::ZEROS(alpha, dstep);

		this->initp = true;
    }

	// mohan add 2010-07-16
	if(this->new_e_iteration)
	{
		for(int i=0; i<GlobalV::NSPIN; i++)
		{
			for(int j=0; j<rstep; j++)
			{
				ModuleBase::GlobalFunc::ZEROS(Rrho[i][j], GlobalC::pw.nrxx);
			}
		}
		
		for(int i=0; i<GlobalV::NSPIN; i++)
		{
			ModuleBase::GlobalFunc::ZEROS(rho_save2[i], GlobalC::pw.nrxx);
			for(int j=0; j<dstep; j++)
			{
				ModuleBase::GlobalFunc::ZEROS( dRrho[i][j], GlobalC::pw.nrxx );
				ModuleBase::GlobalFunc::ZEROS( drho[i][j], GlobalC::pw.nrxx);
			}
		}

		ModuleBase::GlobalFunc::ZEROS(dRR, dstep);
		ModuleBase::GlobalFunc::ZEROS(alpha, dstep);
	}
	return;
}


// calculate < dR | dR >
// if spin is considered, double the size.
// < dR1,dR2 | dR1,dR2 > = < dR1 | dR1 > + < dR2 | dR2 >
void Charge_Pulay::generate_Abar(const int &scheme, ModuleBase::matrix &A)const
{
	int step = 0;

	// default is scheme = 1.
	if(scheme==1) step=dstep;
	if(scheme==2) step=rstep;

	A.zero_out();

	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int i=0; i<step; i++)
		{
			for(int j=0; j<=i; j++)
			{
				if(scheme==1)
				{
					A(i,j) += this->calculate_residual_norm( this->dRrho[is][j], this->dRrho[is][i] );
				}
				else if(scheme==2)
				{
					A(i,j) += this->calculate_residual_norm( this->Rrho[is][j], this->Rrho[is][i] );
				}
				A(j,i) = A(i,j);
			}
		}
	}

	return;
}

#include "../module_base/complexmatrix.h"
void Charge_Pulay::inverse_preA(const int &dim, ModuleBase::matrix &preA)const
{
	ModuleBase::ComplexMatrix B(dim, dim);
	ModuleBase::ComplexMatrix C(dim, dim);
	for(int i=0; i<dim; i++)
	{
		for(int j=0; j<dim; j++)
		{
			B(i,j) = std::complex<double> (preA(i,j), 0.0);
		}
	}
	ModuleBase::Inverse_Matrix_Complex IMC;
	IMC.init(dim);
	IMC.using_zheev(B,C);
//	out.printcm("Inverse B", C);
	for(int i=0; i<dim; i++)
	{
		for(int j=0; j<dim; j++)
		{
			preA(i,j) = C(i,j).real();
		}
	}
	return;
}

void Charge_Pulay::inverse_real_symmetry_matrix(const int &scheme, ModuleBase::matrix &A)const // indicate the spin.
{
//	ModuleBase::TITLE("Charge_Pulay","inverse_Abar");

	int step = 0;

	if(scheme==1) step=dstep;
	if(scheme==2) step=rstep;

	// Notice that it's a symmetry matrix!!!
//	out.printrm("Abar",Abar);	
	ModuleBase::ComplexMatrix B(step,step);
	ModuleBase::ComplexMatrix C(step,step);
	for(int i=0; i<step; i++)
	{
		for(int j=0; j<step; j++)
		{
			B(i,j) = std::complex<double> (A(i,j),0.0);
		}
	}
		
	ModuleBase::Inverse_Matrix_Complex IMC;
	IMC.init(step);
	IMC.using_zheev(B,C);

//	out.printcm("Inverse B", C);

	for(int i=0; i<step; i++)
	{
		for(int j=0; j<step; j++)
		{
			A(i,j) = C(i,j).real();
		}
	}

//	out.printrm("Inverse Abar",A);	
	
	return;
}

void Charge_Pulay::generate_dRR(const int &m)
{
	ModuleBase::GlobalFunc::ZEROS(dRR, dstep);
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for(int i=0; i<dstep; i++)
		{
			this->dRR[i] += this->calculate_residual_norm(this->dRrho[is][i], this->Rrho[is][m]);	
		}
	}

	return;
}

// scheme1 : use dstep to genearte Abar(dstep, dstep)
// scheme2 : use rstep to generate Abar(rstep, rstep)
void Charge_Pulay::generate_alpha(const int &scheme)
{
//	ModuleBase::TITLE("Charge_Pulay","generate_alpha");

	ModuleBase::GlobalFunc::ZEROS(alpha, dstep);
	if(scheme==1)
	{
		for(int i=0; i<dstep; i++)
		{
			for(int j=0; j<dstep; j++)
			{
				// Abar is a symmetry matrix.
				this->alpha[i] -= this->Abar(j,i) * this->dRR[j];
			}	
		}
	}	
	else if(scheme==2)
	{
		double sum = 0.0;
		for(int i=0; i<rstep; i++)
		{
			for(int j=0; j<rstep; j++)
			{
				sum += this->Abar(i,j);	
			} 
		}
		assert(sum!=0.0);
		for(int i=0; i<rstep; i++)
		{
			for(int j=0; j<rstep; j++)
			{
				this->alpha[i] += this->Abar(j,i);
			}
			this->alpha[i] /= sum;
			//std::cout << "\n alpha[" << i << "]=" << alpha[i];
		}

		// test if sum(alpha)=1;
		double suma = 0.0;
		for(int i=0; i<rstep; i++)
		{
			suma += alpha[i];
		}
//		std::cout << "\n suma = " << suma << std::endl;
	}

	return;
}

void Charge_Pulay::generate_new_rho(const int &is, const int &m)
{
//	ModuleBase::TITLE("Charge_Pulay","generate_new_rho");

//	std::cout << " generate_new_rho " << std::endl;
//	this->check_ne(rho[is]);
//	this->check_ne(rho_save[is]);
	
	double mixp = this->mixing_beta;
	
	// rho tmp
	double* rhonew = new double[GlobalC::pw.nrxx];
	ModuleBase::GlobalFunc::ZEROS(rhonew, GlobalC::pw.nrxx);
	
	for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
	{
		rhonew[ir] = this->rho_save[is][ir] + mixp * this->Rrho[is][m][ir];
		for(int i=0; i<dstep; i++)
		{
			rhonew[ir] += this->alpha[i] * ( this->drho[is][i][ir] + mixp * this->dRrho[is][i][ir] );
		}
	}

	ModuleBase::GlobalFunc::DCOPY(rhonew, rho[is], GlobalC::pw.nrxx);
	
	// this is done in save_rho_before_sum_bands
//	ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is], GlobalC::pw.nrxx);


	delete[] rhonew;

	return;
}

void Charge_Pulay::generate_residual_vector(double *residual, const double* rho_out, const double* rho_in)const
{
	for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
	{
		residual[ir]= rho_out[ir] - rho_in[ir];
	}
	//std::cout << "\n Calculate residual norm = " << calculate_residual_norm(residual, residual) << std::endl;
	return;
}

double Charge_Pulay::calculate_residual_norm(double *residual1, double* residual2)const
{
	// calculate the norm of the residual std::vector:
	// (the target to minimize in Pulay's algorithm)
	double rnorm = 0.0;
	for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
	{
		rnorm += residual1[ir]*residual2[ir];
	}
	Parallel_Reduce::reduce_double_pool(rnorm); // mohan fix bug 2010-07-22
	//GlobalV::ofs_running << " rnorm = " << rnorm << std::endl;
	return rnorm;
}


// calculate "dR^{i} = R^{i+1} - R^{i}"
// calculate "drho^{i} = rho^{i+1} - rho^{i}"
void Charge_Pulay::generate_datas(const int &irstep, const int &idstep, const int &totstep)
{
	//===============================================
	// calculate the important "Rrho".
	// init the "Rrho = rho - rho_save"
	// which Rrho to be update now? answer: irstep
	//===============================================

	ModuleBase::GlobalFunc::NOTE("Generate Residual std::vector from rho and rho_save.");
    for (int is=0; is<GlobalV::NSPIN; is++)
    {
//		std::cout << " generate datas , spin=" << is << std::endl;
//		double c1=check_ne(rho[is]);
//		double c2=check_ne(rho_save[is]);
		this->generate_residual_vector( this->Rrho[is][irstep], this->rho[is], this->rho_save[is]);

		if(this->mixing_gg0 > 0.0)
		{
			std::complex<double> *kerpulay = new std::complex<double>[GlobalC::pw.ngmc];
			double* kerpulayR = new double[GlobalC::pw.nrxx];
			
			set_rhog(Rrho[is][irstep], kerpulay);

			const double fac = this->mixing_gg0;
			const double gg0 = std::pow(fac * 0.529177 /GlobalC::ucell.tpiba, 2);
			double* filter_g = new double[GlobalC::pw.ngmc];
			for(int ig=0; ig<GlobalC::pw.ngmc; ig++)
			{
				double gg = GlobalC::pw.get_NormG_cartesian(ig);
				filter_g[ig] = max(gg / (gg + gg0), 0.1);
				kerpulay[ig] = (1 - filter_g[ig]) * kerpulay[ig];
			}

			set_rhor(kerpulay, kerpulayR);
			for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
			{
				Rrho[is][irstep][ir] = Rrho[is][irstep][ir] - kerpulayR[ir];
			}

			delete[] kerpulay;
			delete[] kerpulayR;
			delete[] filter_g;
		}
    }

	if(totstep==0)
	{
		// don't need to calculate 'dRrho' and 'drho'.
	}
	else if(totstep>0)
	{	
		// which dRrho to be update now? answer: dstep=istep-1;
		// irstep range: [0, rstep)
		const int nowR = irstep;
		int lastR = irstep - 1;

		if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "now irstep", nowR);
		if(GlobalV::test_charge)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "last irstep",lastR); 

		if(lastR < 0) lastR += rstep;
		//std::cout << "\n nowR(irstep) = " << nowR;
		//std::cout << "\n lastR(irstep-1) = " << lastR;
		for (int is=0; is<GlobalV::NSPIN; is++)
		{
			for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
			{
				this->dRrho[is][idstep][ir] = this->Rrho[is][nowR][ir] - this->Rrho[is][lastR][ir];
				this->drho[is][idstep][ir] = this->rho_save[is][ir] - this->rho_save2[is][ir];
			}
			//std::cout << "\n Calculate <dR|dR> norm = " 
			//<< calculate_residual_norm(dRrho[is][idstep], dRrho[is][idstep]);
			//std::cout << "\n Calculate <drho|drho> norm = " 
			//<< calculate_residual_norm(drho[is][idstep],drho[is][idstep]) << std::endl;
		}
	}

	// mohan fix the bug 2010/03/26.
	// save 'rho_save2' in order to calculate drho in
	// the next iteration.
	ModuleBase::GlobalFunc::NOTE("Calculate drho = rho_{in}^{i+1} - rho_{in}^{i}");
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		ModuleBase::GlobalFunc::DCOPY(this->rho_save[is], this->rho_save2[is], GlobalC::pw.nrxx);
	}
	return;
}

