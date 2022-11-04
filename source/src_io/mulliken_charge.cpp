/**********************************************************************
  Mulliken_Charge.cpp:

     Mulliken_Charge.cpp is a subrutine to calculate Mulliken charge.
 
  Log of Mulliken_Charge.cpp:

     12/Oct/2018  Released by Feng Qi

***********************************************************************/

#include "mulliken_charge.h"
#include "../src_pw/global.h"
#include "../src_pw/wavefunc.h"
#ifdef __LCAO
#include "../src_lcao/local_orbital_charge.h"
#include "../src_lcao/LCAO_gen_fixedH.h"
#include "../src_lcao/LCAO_matrix.h"
#include "../src_lcao/global_fp.h"
#endif
#include "../module_base/lapack_connector.h"
#include "../module_base/scalapack_connector.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include <vector>
#ifdef __MPI
#include<mpi.h>
#endif
#include "../module_neighbor/sltk_atom_arrange.h"



Mulliken_Charge::Mulliken_Charge(const psi::Psi<double> *psi_gamma_in,
    const psi::Psi<std::complex<double>> *psi_k_in)
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
        this->psi_gamma = psi_gamma_in;
	}
	else 
	{
        this->psi_k = psi_k_in;
	}

	mug = new  std::complex<double>   [GlobalV::NLOCAL];


	DecMulP = new double* [GlobalV::NSPIN];
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		DecMulP[is] = new double[GlobalV::NLOCAL];
		ModuleBase::GlobalFunc::ZEROS(DecMulP[is], GlobalV::NLOCAL);
	}
	MecMulP = new double* [GlobalV::NSPIN];
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		MecMulP[is] = new double[GlobalV::NLOCAL];
		ModuleBase::GlobalFunc::ZEROS(MecMulP[is], GlobalV::NLOCAL);
	}


	ADecMulP = new double**[GlobalV::NSPIN];
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		ADecMulP[is] = new double*[GlobalC::ucell.nat];

		for (int i=0; i<GlobalC::ucell.nat; i++)
		{
			ADecMulP[is][i] = new double[(2*GlobalC::ucell.lmax+1)*(2*GlobalC::ucell.lmax+1)*GlobalC::ucell.nmax];
			ModuleBase::GlobalFunc::ZEROS(ADecMulP[is][i], (2*GlobalC::ucell.lmax+1)*(2*GlobalC::ucell.lmax+1)*GlobalC::ucell.nmax);
		}
	}
}

Mulliken_Charge::~Mulliken_Charge()
{  

	delete[] mug;

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		delete[] DecMulP[is];
	}	
	delete[] DecMulP;
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		delete[] MecMulP[is];
	}	
	delete[] MecMulP;

	for (int is=0; is<GlobalV::NSPIN; is++)
	{
		for (int i=0; i<GlobalC::ucell.nat; i++)
		{
			delete[] ADecMulP[is][i]; 
		}
		delete[] ADecMulP[is];
	}
	delete[] ADecMulP;

}

  
void Mulliken_Charge::cal_mulliken(LCAO_Hamilt &uhm)
{
    ModuleBase::TITLE("Mulliken_Charge", "cal_mulliken");
    const Parallel_Orbitals* pv = uhm.LM->ParaV;

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		if(GlobalV::GAMMA_ONLY_LOCAL)
		{
			std::vector<ModuleBase::matrix>   mud;
			mud.resize(1);
			mud[0].create(pv->ncol,pv->nrow);

			this->psi_gamma->fix_k(is);
			const double* ppsi = this->psi_gamma->get_pointer();
			for (int i=0; i<GlobalV::NBANDS; ++i)		  
			{     
				ModuleBase::GlobalFunc::ZEROS(mug, GlobalV::NLOCAL);
				const int NB= i+1;

				const double one_float=1.0, zero_float=0.0;
				const int one_int=1;

			#ifdef __MPI
				const char T_char='T';		
				pdgemv_(
						&T_char,
						&GlobalV::NLOCAL,&GlobalV::NLOCAL,
						&one_float,
						uhm.LM->Sloc.data(), &one_int, &one_int, pv->desc,
						ppsi, &one_int, &NB, pv->desc, &one_int,
						&zero_float,
						mud[0].c, &one_int, &NB, pv->desc,
						&one_int);
			#endif

				for (int j=0; j<GlobalV::NLOCAL; ++j)
				{

					if ( pv->in_this_processor(j,i) )
					{

						const int ir = pv->trace_loc_row[j];
						const int ic = pv->trace_loc_col[i];

						mug[j] = mud[0](ic,ir)*this->psi_gamma[0](ic,ir);

						const double x = mug[j].real();

						MecMulP[is][j] +=x*GlobalC::wf.wg(0,i);
					}
				} 
			}//ib
		}//if
		else
		{   
			std::vector<ModuleBase::ComplexMatrix> mud;
			mud.resize(1);
			mud[0].create(pv->ncol,pv->nrow);

			GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
				GlobalV::ofs_running,
				GlobalV::OUT_LEVEL,
				GlobalC::ORB.get_rcutmax_Phi(), 
				GlobalC::ucell.infoNL.get_rcutmax_Beta(), 
				GlobalV::GAMMA_ONLY_LOCAL);

			atom_arrange::search(
				GlobalV::SEARCH_PBC,
				GlobalV::ofs_running,
				GlobalC::GridD, 
				GlobalC::ucell, 
				GlobalV::SEARCH_RADIUS,
				GlobalV::test_atom_input);//qifeng-2019-01-21


			uhm.LM->allocate_HS_R(pv->nnr);
			uhm.LM->zeros_HSR('S');
			uhm.genH.calculate_S_no(uhm.LM->SlocR.data());
			uhm.genH.build_ST_new('S', false, GlobalC::ucell, uhm.LM->SlocR.data());

			for(int ik=0;ik<GlobalC::kv.nks;ik++)
			{
				if(is == GlobalC::kv.isk[ik])
				{
					uhm.LM->allocate_HS_k(pv->nloc);
					uhm.LM->zeros_HSk('S');
					uhm.LM->folding_fixedH(ik);

					this->psi_k->fix_k(ik);
					psi::Psi<std::complex<double>> Dwfc(this->psi_k[0], 1);
					std::complex<double>* p_dwfc = Dwfc.get_pointer();
					for(int index = 0; index < Dwfc.size(); ++index)
					{
						p_dwfc[index] = conj(p_dwfc[index]);
					}

					for (int i=0; i<GlobalV::NBANDS; ++i)		  
					{     
						ModuleBase::GlobalFunc::ZEROS(mug, GlobalV::NLOCAL);

						const int NB= i+1;

						const double one_float[2]={1.0, 0.0}, zero_float[2]={0.0, 0.0};
						const int one_int=1;
						//   const int two_int=2;

						const char T_char='T';		// N_char='N',U_char='U'

					#ifdef __MPI
						pzgemv_(
								&T_char,
								&GlobalV::NLOCAL,&GlobalV::NLOCAL,
								&one_float[0],
								uhm.LM->Sloc2.data(), &one_int, &one_int, pv->desc,
								p_dwfc, &one_int, &NB, pv->desc, &one_int,
								&zero_float[0],
								mud[0].c, &one_int, &NB, pv->desc,
								&one_int);
					#endif

						for (int j=0; j<GlobalV::NLOCAL; ++j)
						{

							if ( pv->in_this_processor(j,i) )
							{

								const int ir = pv->trace_loc_row[j];
								const int ic = pv->trace_loc_col[i];

								mug[j] = mud[0](ic,ir)*this->psi_k[0](ic,ir);
								const double x = mug[j].real();
								MecMulP[is][j] +=x*GlobalC::wf.wg(ik,i);
								// std::cout <<   wavog[j] << std::endl; 
							}
						}                             
					}//ib
				}//if                       
			}//ik
#ifdef __MPI
			atom_arrange::delete_vector(
				GlobalV::ofs_running,
				GlobalV::SEARCH_PBC,
				GlobalC::GridD, 
				GlobalC::ucell, 
				GlobalV::SEARCH_RADIUS, 
				GlobalV::test_atom_input);
#endif

		}//else          
	#ifdef __MPI           
		MPI_Reduce(MecMulP[is], DecMulP[is] , GlobalV::NLOCAL , MPI_DOUBLE , MPI_SUM, 0, MPI_COMM_WORLD);
	#endif
		if(GlobalV::MY_RANK == 0)
		{
			for (int i=0; i<GlobalC::ucell.nat; i++)
			{   
				int a = GlobalC::ucell.iat2ia[i];
				int t = GlobalC::ucell.iat2it[i];
				Atom* atom1 = &GlobalC::ucell.atoms[t];
				for(int j=0; j<atom1->nw; ++j)
				{
					int k = GlobalC::ucell.itiaiw2iwt(t,a,j);
					ADecMulP[is][i][j] = DecMulP[is][k];

				}
			}

		}

	}//is

	return;                									
}				   

void Mulliken_Charge::stdout_mulliken(LCAO_Hamilt &uhm)
{                    
	this->cal_mulliken(uhm);
	
	if(GlobalV::MY_RANK == 0)
	{
		ModuleBase::TITLE("Dos","calculate_Mulliken");
		std::ofstream fout;
		std::stringstream ss;
		ss << GlobalV::global_out_dir <<"mulliken.txt";
		fout.open(ss.str().c_str());
		// std::ofstream fout;
		// std::string wordqf="mulliken.txt";
		// wordqf += char(GlobalV::MY_RANK + 48);
		//  wordqf += ".txt";
		//   fout.open(wordqf.c_str(),ios::app);     
		int num,l,m,mul;
		double Tcharge;
		double*   sum_l = new double[2];
		double*   sum_mul = new double[2];

		std::string Name_Angular[5][11];
		/* decomposed Mulliken charge */

		Name_Angular[0][0] = "s          ";
		Name_Angular[1][0] = "px         ";
		Name_Angular[1][1] = "py         ";
		Name_Angular[1][2] = "pz         ";
		Name_Angular[2][0] = "d3z^2-r^2  ";
		Name_Angular[2][1] = "dxy        ";
		Name_Angular[2][2] = "dxz        ";
		Name_Angular[2][3] = "dx^2-y^2   ";
		Name_Angular[2][4] = "dyz        ";
		Name_Angular[3][0] = "f5z^2-3r^2 ";
		Name_Angular[3][1] = "f5xz^2-xr^2";
		Name_Angular[3][2] = "f5yz^2-yr^2";
		Name_Angular[3][3] = "fzx^2-zy^2 ";
		Name_Angular[3][4] = "fxyz       ";
		Name_Angular[3][5] = "fx^3-3*xy^2";
		Name_Angular[3][6] = "f3yx^2-y^3 ";
		Name_Angular[4][0] = "g1         ";
		Name_Angular[4][1] = "g2         ";
		Name_Angular[4][2] = "g3         ";
		Name_Angular[4][3] = "g4         ";
		Name_Angular[4][4] = "g5         ";
		Name_Angular[4][5] = "g6         ";
		Name_Angular[4][6] = "g7         ";
		Name_Angular[4][7] = "g8         ";
		Name_Angular[4][8] = "g9         ";
		fout << "\n CALCULATE THE MULLIkEN ANALYSIS FOR EACH ATOM" << std::endl;

		// calculate the total charge of the system.
		double sch = 0.0;
		fout << std::setprecision(8);
		for(int is=0; is<GlobalV::NSPIN; ++is)
		{
			double sss = 0.0;
			for(int iw=0; iw<GlobalV::NLOCAL; ++iw)
			{
				sch += DecMulP[is][iw];
				sss += DecMulP[is][iw];


			}
			fout << sss << " (Total charge all spin " << is+1 << ")" << std::endl;
		}
		fout << sch << " (Total charge of the system)" << std::endl;
		fout << "  Decomposed Mulliken populations" << std::endl;

		/*for (int i=0; i<GlobalC::ucell.nat; i++)
		  {
		  int   t = GlobalC::ucell.iat2it[i];
		  Atom* atom1 = &GlobalC::ucell.atoms[t];  
		  fout << i << "         " << GlobalC::ucell.atoms[t].label << "			" << "Up spin" << "                    " << "Down spin" << std::endl;
		  fout << "            " << "multiple" << std::endl;
		//num = 0;
		for(int j=0; j<atom1->nw; ++j)
		{
		const int L1 = atom1->iw2l[j];
		const int N1 = atom1->iw2n[j];
		const int m1 = atom1->iw2m[j];
		fout << Name_Angular[L1][m1] << "      " << N1 << "          " << ADecMulP[0][i][j] << "               " << ADecMulP[0][i][j] << std::endl;
		}
		}*/
		for (int i=0; i<GlobalC::ucell.nat; i++)
		{
			Tcharge = 0.0;
			int t = GlobalC::ucell.iat2it[i];
			if (GlobalV::NSPIN==1 || GlobalV::NSPIN==2)
			{
				fout << i << std::setw(25) << GlobalC::ucell.atoms[t].label 
					<< std::setw(30) << "Up spin" 
					<< std::setw(30) << "Down spin" 
					<< std::setw(30) << "Sum"<< std::setw(30) << "Diff"<< std::endl;
			}
			fout << std::setw(29) << "multiple" <<std::endl;

			num = 0;
			int lm = GlobalC::ucell.atoms[t].nwl;

			for (l=0; l<=lm; l++)
			{             
				if (GlobalV::NSPIN==1)
				{
					sum_l[0] = 0.0;
				}
				else if (GlobalV::NSPIN==2)
				{
					sum_l[0] = 0.0;
					sum_l[1] = 0.0;
				}
				/* else if (GlobalV::NSPIN==3){
				   sum_l[0] = 0.0;
				   sum_l[1] = 0.0;
				   }*/
				for (mul=0; mul<GlobalC::ucell.atoms[t].l_nchi[l]; mul++)
				{                 
					if (GlobalV::NSPIN==1){
						sum_mul[0] = 0.0;
					}
					else if (GlobalV::NSPIN==2){
						sum_mul[0] = 0.0;
						sum_mul[1] = 0.0;
					}
					/* else if (SpinP_switch==3){
					   sum_mul[0] = 0.0;
					   sum_mul[1] = 0.0;
					   }*/
					for (m=0; m<(2*l+1); m++)
					{     
						if (GlobalV::NSPIN==1){
							ADecMulP[0][i][num] = 0.5*ADecMulP[0][i][num]; 
							fout << Name_Angular[l][m] << std::setw(14) 
							<< mul << std::setw(32)<< ADecMulP[0][i][num] 
							<< std::setw(30)<< ADecMulP[0][i][num] 
							<< std::setw(30)<< ADecMulP[0][i][num]
							+ ADecMulP[0][i][num]
							<< std::setw(28) << ADecMulP[0][i][num]-ADecMulP[0][i][num]<< std::endl;
							sum_mul[0] += ADecMulP[0][i][num];
						}
						else if (GlobalV::NSPIN==2)
						{
							fout << Name_Angular[l][m] 
								<< std::setw(14) << mul 
								<< std::setw(32) << ADecMulP[0][i][num] 
								<< std::setw(30) << ADecMulP[1][i][num] 
								<< std::setw(30) << ADecMulP[0][i][num]
								+ ADecMulP[1][i][num]<< std::setw(28) 
								<< ADecMulP[0][i][num]-ADecMulP[1][i][num]<< std::endl;
							sum_mul[0] += ADecMulP[0][i][num];
							sum_mul[1] += ADecMulP[1][i][num];
						}
						num++;
					}

					if (GlobalV::NSPIN==1)
					{
						fout <<"  sum over m "<< std::setw(43) 
							<<sum_mul[0]<< std::setw(30) <<sum_mul[0]
							<< std::setw(35) <<sum_mul[0]+sum_mul[0]
							<< std::setw(25) <<sum_mul[0]-sum_mul[0]<<std::endl;
						sum_l[0] += sum_mul[0];
					}
					else if (GlobalV::NSPIN==2)
					{
						fout <<"  sum over m "<< std::setw(43) <<sum_mul[0]
						<< std::setw(30) <<sum_mul[0]
						<< std::setw(35) <<sum_mul[0]+sum_mul[1]<< std::setw(25) <<sum_mul[0]-sum_mul[1]<<std::endl;
						sum_l[0] += sum_mul[0];
						sum_l[1] += sum_mul[1];
					}

				}

				if (GlobalC::ucell.atoms[t].l_nchi[l]!=0)
				{
					if (GlobalV::NSPIN==1)
					{
						fout <<"  sum over m+mul "<< std::setw(36) 
							<<sum_l[0]<< std::setw(30) <<sum_l[0]
							<< std::setw(33) <<sum_l[0]+sum_l[0]
							<< std::setw(29) <<sum_l[0]-sum_l[0]<<std::endl;
						Tcharge =  Tcharge+ sum_l[0]+sum_l[0];
					}
					else if (GlobalV::NSPIN==2)
					{
						fout <<"  sum over m+mul "<< std::setw(36) 
							<<sum_l[0]<< std::setw(30) <<sum_l[1]
							<< std::setw(33) <<sum_l[0]+sum_l[1]
							<< std::setw(29) <<sum_l[0]-sum_l[1]<<std::endl;
						Tcharge =  Tcharge+sum_l[0]+sum_l[1];
					}
				}
			}
			fout <<"Total Charge on atom  "<< GlobalC::ucell.atoms[t].label <<  std::setw(20) << Tcharge <<std::endl<<std::endl<<std::endl;
		}
		fout.close();
		delete[] sum_l;
		delete[] sum_mul;
	}

	return;
}




