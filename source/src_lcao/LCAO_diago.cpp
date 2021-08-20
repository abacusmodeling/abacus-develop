#include "LCAO_diago.h"
#include "../src_pw/global.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_io/write_HS.h"
//xiaohui modified 2013-03-23
//#include "../src_develop/src_siao/selinv.h"

Diago_LCAO_Matrix::Diago_LCAO_Matrix(){}
Diago_LCAO_Matrix::~Diago_LCAO_Matrix(){}

void Diago_LCAO_Matrix::solve_complex_matrix(
	const int &ik, 
	std::complex<double>** wfc, 
	ModuleBase::ComplexMatrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","solve_complex_matrix");
	time_t time_start = time(NULL);

	if(GlobalV::KS_SOLVER=="lapack")
	{
		this->using_LAPACK_complex(ik, wfc);
	}
	else
	{
#ifdef __MPI
		this->using_HPSEPS_complex(ik, wfc, wfc_2d);
#else
		WARNING_QUIT("Diago_LCAO_Matrix::solve_complex_matrix","only lapack is available!");
#endif
	}

	time_t time_end = time(NULL);

	OUT_TIME("diago(std::complex)", time_start, time_end);

	return;
}


void Diago_LCAO_Matrix::solve_double_matrix(
	const int &ik, 
	double** wfc, 
	matrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","solve_double_matrix");
	timer::tick("Diago_LCAO_Matrix","solve_double_matrix");
	time_t time_start = time(NULL);

	
	if(GlobalV::KS_SOLVER=="lapack")
	{
		this->using_LAPACK(ik, wfc);
	}
#ifdef __MPI
	else if(GlobalV::KS_SOLVER=="hpseps" || GlobalV::KS_SOLVER=="genelpa"|| GlobalV::KS_SOLVER=="scalapack_gvx")
	{
		this->using_HPSEPS_double(ik, wfc, wfc_2d);
	}
#endif
	else
	{
		std::cout << " Diago_LCAO_Matrix, diago_type = " << GlobalV::KS_SOLVER << std::endl; 
		WARNING_QUIT("Diago_LCAO_Matrix::init","Check GlobalV::KS_SOLVER.");
	}

	time_t time_end = time(NULL);

	OUT_TIME("diago(double)",time_start, time_end);

	timer::tick("Diago_LCAO_Matrix","solve_double_matrix");
	return;
}

#ifdef __MPI
void Diago_LCAO_Matrix::using_HPSEPS_double(const int &ik, double**wfc, matrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","using_HPSEPS_double");

	// save H and S matrix to disk.
	bool bit = false;
	HS_Matrix::saving_HS(GlobalC::LM.Hloc, GlobalC::LM.Sloc, bit, GlobalC::ParaO.out_hs);
	GlobalV::ofs_running << std::setprecision(6);

	// Distribution of matrix for 
	// prallel eigensolver.
	GlobalC::ParaO.diago_double_begin(ik, wfc, wfc_2d, GlobalC::LM.Hloc, GlobalC::LM.Sloc, GlobalC::wf.ekb[ik]);

	return;
}


void Diago_LCAO_Matrix::using_HPSEPS_complex(const int &ik, std::complex<double>** wfc, ModuleBase::ComplexMatrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","using_HPSEPS_complex");

	//GlobalC::ParaO.out_hs=1;//zhengdy-soc-test
	bool bit = false; //LiuXh, 2017-03-21
	//if set bit = true, there would be error in soc-multi-core calculation, noted by zhengdy-soc
	HS_Matrix::saving_HS_complex(GlobalC::LM.Hloc2, GlobalC::LM.Sloc2, bit, GlobalC::ParaO.out_hs); //LiuXh, 2017-03-21
	GlobalV::ofs_running << std::setprecision(6); //LiuXh, 2017-03-21

	GlobalC::ParaO.diago_complex_begin(ik, wfc, wfc_2d, GlobalC::LM.Hloc2, GlobalC::LM.Sloc2, GlobalC::wf.ekb[ik]);

	//added by zhengdy-soc, rearrange the WFC_K from [up,down,up,down...] to [up,up...down,down...], 
	if(GlobalV::NSPIN==4)
	{
		int row = GlobalC::GridT.lgd;
		std::vector<std::complex<double>> tmp(row);
		for(int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			for(int iw=0; iw<row / GlobalV::NPOL; iw++)
			{
				tmp[iw] = GlobalC::LOWF.WFC_K[ik][ib][iw * GlobalV::NPOL];
				tmp[iw + row / GlobalV::NPOL] = GlobalC::LOWF.WFC_K[ik][ib][iw * GlobalV::NPOL + 1];
			}
			for(int iw=0; iw<row; iw++)
			{
				GlobalC::LOWF.WFC_K[ik][ib][iw] = tmp[iw];
			}
		}
	}

	return;
}
#endif

void Diago_LCAO_Matrix::using_LAPACK_complex(const int &ik, std::complex<double> **wfc)const
{
	TITLE("Diago_LCAO_Matrix","using_LAPACK_complex");

	assert(GlobalV::NPROC = 1);

	ModuleBase::ComplexMatrix Htmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	ModuleBase::ComplexMatrix Stmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			Htmp(i,j) = GlobalC::LM.Hloc2[i*GlobalV::NLOCAL+j];
			Stmp(i,j) = GlobalC::LM.Sloc2[i*GlobalV::NLOCAL+j];
		}
	}

//----------------------------
// keep this for tests
//    out.printcm_norm("Lapack_H", Htmp, 1.0e-5);
//    out.printcm_norm("Lapack_S", Stmp, 1.0e-5);
//----------------------------

	double* en = new double[GlobalV::NLOCAL];
	ZEROS(en, GlobalV::NLOCAL);

	ModuleBase::ComplexMatrix hvec(GlobalV::NLOCAL, GlobalV::NBANDS);
	GlobalC::hm.diagH_LAPACK(GlobalV::NLOCAL, GlobalV::NBANDS, Htmp, Stmp, GlobalV::NLOCAL, en, hvec);

	if(GlobalV::NSPIN!=4)
	{
		for(int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			for(int iw=0; iw<GlobalV::NLOCAL; iw++)
			{
				GlobalC::LOWF.WFC_K[ik][ib][iw] = hvec(iw,ib);
			}
		}
	}
	else
	{
		for(int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			for(int iw=0; iw<GlobalV::NLOCAL / GlobalV::NPOL; iw++)
			{
				GlobalC::LOWF.WFC_K[ik][ib][iw] = hvec(iw * GlobalV::NPOL, ib);
				GlobalC::LOWF.WFC_K[ik][ib][iw + GlobalV::NLOCAL / GlobalV::NPOL] = hvec(iw * GlobalV::NPOL + 1, ib);
			}
		}
	}

	// energy for k-point ik
	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		GlobalC::wf.ekb[ik][ib] = en[ib];
	}

	return;
}


void Diago_LCAO_Matrix::using_LAPACK(const int &ik, double** wfc)const
{
	TITLE("Diago_LCAO_Matrix","using_LAPACK");
	assert(GlobalV::NLOCAL>0);

	// save H and S matrix to disk.
//	bool bit = false;
	bool bit = true;//zhengdy-soc
	HS_Matrix::saving_HS(GlobalC::LM.Hloc, GlobalC::LM.Sloc, bit, GlobalC::ParaO.out_hs);

	matrix Htmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	matrix Stmp(GlobalV::NLOCAL,GlobalV::NLOCAL);
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			Htmp(i,j) = GlobalC::LM.Hloc[i*GlobalV::NLOCAL+j];
			Stmp(i,j) = GlobalC::LM.Sloc[i*GlobalV::NLOCAL+j];
		}
	}

	// @@@@@@@
	// test
	// @@@@@@@
//    out.printrm("Lapack_H", Htmp);
//    out.printrm("Lapack_S", Stmp);
		
	int itype=1;
	int lwork=3*GlobalV::NLOCAL-1;// tmp
	double* w = new double[GlobalV::NLOCAL];
	double* work = new double[lwork];
	ZEROS(w, GlobalV::NLOCAL);
	ZEROS(work, lwork);
	int info;

	clock_t clock_start, clock_end;
	clock_start = std::clock();
	LapackConnector::dsygv(itype,'V','U',GlobalV::NLOCAL,Htmp,GlobalV::NLOCAL,Stmp,GlobalV::NLOCAL,w,work,lwork,&info);
	clock_end = std::clock();
	double duration = (double)(clock_end - clock_start) / CLOCKS_PER_SEC;

	GlobalV::ofs_running << std::setiosflags(ios::fixed) << std::setprecision(20);
//	GlobalV::ofs_running << " clock_start = " << clock_start << std::endl;
//	GlobalV::ofs_running << " clock_end = " << clock_end << std::endl;
	GlobalV::ofs_running << " Time using dsygv in LAPACK (seconds) is " << duration << std::endl;
	GlobalV::ofs_running << std::resetiosflags(ios::fixed) << std::setprecision(10);

	for(int i=0; i<GlobalV::NBANDS; i++)
	{
		// eigenvalues
		GlobalC::wf.ekb[ik][i] = w[i]; 
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			wfc[i][j] = Htmp(j,i);
		}
	}
	

	// @@@@@@@
	// test
	// @@@@@@@
	/*
	std::cout << "\n Lapack, wfc after diago:" << std::endl;
	for(int i=0; i<GlobalV::NBANDS; i++)
	{
		std::cout << " Eigenvalue from LAPACK : " << std::setw(5) << std::setw(12) << GlobalC::wf.ekb[ik][i] << std::endl;
		std::cout << " Eigenfunctions" << std::endl;
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			std::cout << std::setw(12) << wfc[i][j];
		}
		std::cout << std::endl;
	}
	//exit(0);
	*/

	delete[] w;
	delete[] work;
	return;
}
