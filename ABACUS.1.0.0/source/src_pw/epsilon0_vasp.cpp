//===============================================================================
//   AUTHOR : Pengfei Li
//   DATE : 2016-12-30
//===============================================================================
//-------------------------------------------------------------------------------
//this part is for calculating the macroscopic dielectric tensor  using the linear
//response theory .
//-------------------------------------------------------------------------------
#include "global.h"
#include "algorithms.h"
#include "hamilt_pw.h"
#include "../src_pw/wavefunc_in_pw.h"
#include "optical.h"
#include "epsilon0_vasp.h"
#include <iostream>
#include <math.h>
using namespace std;

Epsilon0_vasp epsilon0_vasp;

Epsilon0_vasp::Epsilon0_vasp()
               :init_finish(false)
{
	epsilon = false;
	domega = 0.01;
	nomega = 300;
	eta = 0.01;	
}

Epsilon0_vasp::~Epsilon0_vasp()
{	
}

void Epsilon0_vasp::cal_epsilon0()
{
	cout << "domega = "<<domega<<endl;
	cout << "nomega = "<<nomega<<endl;
	cout << "eta = "<<eta<<endl;
	
	double occupied_bands = static_cast<double>(ucell.nelec/DEGSPIN);
	if( (occupied_bands - std::floor(occupied_bands)) > 0.0 )
	{
		occupied_bands = std::floor(occupied_bands) + 1.0;
	}
	
	oband = int(occupied_bands + 0.5);
	uband = NBANDS - oband;
	
	cout << "oband = "<<oband<<endl;
	cout << "uband = "<<uband<<endl;
	
	if( !init_finish )
	{
		Init();
		init_finish = true;
	}
	
	Cal_epsilon0s();
	
	Cal_T();
	Cal_epsilon0();
	
	ofs_running<<endl; 
	ofs_running<<" The real part of the macroscopic dielectric constant:"<<endl;
	ofs_running<<setw(15)<<"omega"<<setw(15)<<"XX"<<setw(15)<<"XY"<<setw(15)<<"XZ"<<setw(15)<<"YX"<<setw(15)<<"YY"<<setw(15)<<"YZ"<<setw(15)<<"ZX"<<setw(15)<<"ZY"<<setw(15)<<"ZZ"<<endl;  
	for(int i=0; i<nomega; i++)
	{
		ofs_running<<setprecision(2)<<setw(15)<<(i*domega)<<setprecision(2)<<setw(15)<<eps0[0][i].real()+1.0<<setprecision(2)<<setw(15)<<eps0[1][i].real()<<setprecision(2)<<setw(15)<<eps0[2][i].real()<<setprecision(2)<<setw(15)<<eps0[3][i].real()<<setprecision(2)<<setw(15)<<eps0[4][i].real()+1.0<<setprecision(2)<<setw(15)<<eps0[5][i].real()<<setprecision(2)<<setw(15)<<eps0[6][i].real()<<setprecision(2)<<setw(15)<<eps0[7][i].real()<<setprecision(2)<<setw(15)<<eps0[8][i].real()+1.0<<endl;
	}
	
	ofs_running<<endl; 
	ofs_running<<" The imag part of the macroscopic dielectric constant:"<<endl;
	ofs_running<<setw(15)<<"omega"<<setw(15)<<"XX"<<setw(15)<<"XY"<<setw(15)<<"XZ"<<setw(15)<<"YX"<<setw(15)<<"YY"<<setw(15)<<"YZ"<<setw(15)<<"ZX"<<setw(15)<<"ZY"<<setw(15)<<"ZZ"<<endl;  
	for(int i=0; i<nomega; i++)
	{
		ofs_running<<setprecision(2)<<setw(15)<<(i*domega)<<setprecision(2)<<setw(15)<<-eps0[0][i].imag()<<setprecision(2)<<setw(15)<<-eps0[1][i].imag()<<setprecision(2)<<setw(15)<<-eps0[2][i].imag()<<setprecision(2)<<setw(15)<<-eps0[3][i].imag()<<setprecision(2)<<setw(15)<<-eps0[4][i].imag()<<setprecision(2)<<setw(15)<<-eps0[5][i].imag()<<setprecision(2)<<setw(15)<<-eps0[6][i].imag()<<setprecision(2)<<setw(15)<<-eps0[7][i].imag()<<setprecision(2)<<setw(15)<<-eps0[8][i].imag()<<endl;
	}
	
	ofs_running<<" Macroscopic dielectric constant matrix :"<<endl;
	ofs_running<<setprecision(2)<<setw(15)<<eps0[0][0].real()+1.0<<setprecision(2)<<setw(15)<<eps0[1][0].real()<<setprecision(2)<<setw(15)<<eps0[2][0].real()<<endl;
	ofs_running<<setprecision(2)<<setw(15)<<eps0[3][0].real()<<setprecision(2)<<setw(15)<<eps0[4][0].real()+1.0<<setprecision(2)<<setw(15)<<eps0[5][0].real()<<endl;
	ofs_running<<setprecision(2)<<setw(15)<<eps0[6][0].real()<<setprecision(2)<<setw(15)<<eps0[7][0].real()<<setprecision(2)<<setw(15)<<eps0[8][0].real()+1.0<<endl;	
	
	//cout <<"Macroscopic dielectric constant matrix :"<<endl;
	//cout << eps0[0][0].real()+1.0 <<"  "<<eps0[1][0].real() <<"  "<<eps0[2][0].real() <<"  "<<endl;
	//cout << eps0[3][0].real() <<"  "<<eps0[4][0].real()+1.0 <<"  "<<eps0[5][0].real() <<"  "<<endl;
	//cout << eps0[6][0].real() <<"  "<<eps0[7][0].real() <<"  "<<eps0[8][0].real()+1.0 <<"  "<<endl;
	
	Delete();
	
	return;
}

void Epsilon0_vasp:: Init()
{
	b = new complex<double> **[oband];
	for(int ib1=0; ib1<oband; ib1++)
	{
		b[ib1] = new complex<double> *[uband];
		for(int ib2=0; ib2<uband; ib2++)
		{
			b[ib1][ib2] = new complex<double>[3];
		}
	}
	
	eps0s = new complex<double> *[9];
	for(int i=0; i<9; i++)
	{
		eps0s[i] = new complex<double>[nomega];
	}
	
	T = new complex<double> *[nomega];
	for(int i=0; i<nomega; i++)
	{
		T[i] = new complex<double>[nomega];
	}

	eps0 = new complex<double> *[9];
	for(int i=0; i<9; i++)
	{
		eps0[i] = new complex<double>[nomega];
	}

	return;
}

void Epsilon0_vasp:: Delete()
{
	for(int ib1=0; ib1<oband; ib1++)
	{
		for(int ib2=0; ib2<uband; ib2++)
		{
			delete[] b[ib1][ib2];
		}
		delete[] b[ib1];
	}
	delete[] b;
	
	for(int i=0; i<9; i++)
	{
		delete eps0s[i];
	}
	delete[] eps0s;
	
	for(int i=0; i<nomega; i++)
	{
		delete[] T[i];
	}
	delete[] T;
	
	for(int i=0; i<9; i++)
	{
		delete eps0[i];
	}
	delete[] eps0;
	
	return;
}

void Epsilon0_vasp:: Cal_b(int ik)
{
	complex<double> b_core[oband][uband][3];
	
	for(int ib1=0; ib1<oband; ib1++)
		for(int ib2=0; ib2<uband; ib2++)
			for(int i=0; i<3; i++)
			{
				b_core[ib1][ib2][i] = complex<double>(0.0,0.0);
			}
			
	for(int ib1=0; ib1<oband; ib1++)
		for(int ib2=0; ib2<uband; ib2++)
			for(int ig =0; ig< kv.ngk[ik]; ig++)
			{
				/*b_core[ib1][ib2][0] += conj(wf.evc[ik](ib1,ig)) * ((pw.gcar[ig].x)*(TWO_PI/ucell.lat0)) * wf.evc[ik](oband+ib2,ig);
				b_core[ib1][ib2][1] += conj(wf.evc[ik](ib1,ig)) * ((pw.gcar[ig].y)*(TWO_PI/ucell.lat0)) * wf.evc[ik](oband+ib2,ig);
				b_core[ib1][ib2][2] += conj(wf.evc[ik](ib1,ig)) * ((pw.gcar[ig].z)*(TWO_PI/ucell.lat0)) * wf.evc[ik](oband+ib2,ig);*/
				b_core[ib1][ib2][0] += conj(wf.evc[ik](ib1,ig)) * ((kv.kvec_c[ik].x + pw.gcar[ig].x)*(TWO_PI/ucell.lat0)) * wf.evc[ik](oband+ib2,ig);
				b_core[ib1][ib2][1] += conj(wf.evc[ik](ib1,ig)) * ((kv.kvec_c[ik].y + pw.gcar[ig].y)*(TWO_PI/ucell.lat0)) * wf.evc[ik](oband+ib2,ig);
				b_core[ib1][ib2][2] += conj(wf.evc[ik](ib1,ig)) * ((kv.kvec_c[ik].z + pw.gcar[ig].z)*(TWO_PI/ucell.lat0)) * wf.evc[ik](oband+ib2,ig);
			}
			
	double b_core_R[oband][uband][3];
	double b_core_I[oband][uband][3];
	double b_R[oband][uband][3];
	double b_I[oband][uband][3];
	
	for(int ib1=0; ib1<oband; ib1++)
		for(int ib2=0; ib2<uband; ib2++)
			for(int i=0; i<3; i++)
			{
				b_core_R[ib1][ib2][i] = b_core[ib1][ib2][i].real();
				b_core_I[ib1][ib2][i] = b_core[ib1][ib2][i].imag();
			}
			
	MPI_Allreduce(b_core_R,b_R,3*oband*uband,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
	MPI_Allreduce(b_core_I,b_I,3*oband*uband,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
	
	for(int ib1=0; ib1<oband; ib1++)
		for(int ib2=0; ib2<uband; ib2++)
			for(int i=0; i<3; i++)	
			{
				b[ib1][ib2][i] = complex<double>( b_R[ib1][ib2][i], b_I[ib1][ib2][i]);
			}
		
	return;
}

void Epsilon0_vasp:: Cal_epsilon0s()
{
	for(int i=0; i<9; i++)
		for(int j=0; j<nomega; j++)
		{
			eps0s[i][j] = complex<double>(0.0, 0.0);
		}
		
	for(int ik=0; ik<kv.nks; ik++)
	{
		Cal_b(ik);
		/*for(int ib1=0; ib1<oband; ib1++)
			for(int ib2=0; ib2<uband; ib2++)
			{
				cout<<"ik = "<<ik<<" b["<<ib1<<"]["<<ib2<<"][0]"<<b[ib1][ib2][0]<<endl;
			}*/
			
		for(int ib1=0; ib1<oband; ib1++)
			for(int ib2=0; ib2<uband; ib2++)
			{
				double delta_e = wf.ekb[ik][oband+ib2] - wf.ekb[ik][ib1];
				if ((delta_e > 0 || delta_e == 0) && delta_e < ((nomega-1) * domega) )
				{
					int n = int(delta_e/domega);
					double e1 = double(n) * domega;
					double e2 = double(n+1) * domega;
					complex<double> weight1 = complex<double>(32 * (wf.wg(ik,ib1) - wf.wg(ik,oband+ib2)) * (e2 - delta_e)/domega/ucell.omega * PI * PI /delta_e/delta_e/domega, 0.0);
					complex<double> weight2 = complex<double>(32 * (wf.wg(ik,ib1) - wf.wg(ik,oband+ib2)) * (delta_e - e1)/domega/ucell.omega * PI * PI /delta_e/delta_e/domega, 0.0);
					
					eps0s[0][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][0]);
					eps0s[1][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][1]);
					eps0s[2][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][2]);
					eps0s[3][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][0]);
					eps0s[4][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][1]);
					eps0s[5][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][2]);
					eps0s[6][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][0]);
					eps0s[7][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][1]);
					eps0s[8][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][2]);
					
					eps0s[0][n+1] += weight2 * b[ib1][ib2][0] * conj(b[ib1][ib2][0]);
					eps0s[1][n+1] += weight2 * b[ib1][ib2][0] * conj(b[ib1][ib2][1]);
					eps0s[2][n+1] += weight2 * b[ib1][ib2][0] * conj(b[ib1][ib2][2]);
					eps0s[3][n+1] += weight2 * b[ib1][ib2][1] * conj(b[ib1][ib2][0]);
					eps0s[4][n+1] += weight2 * b[ib1][ib2][1] * conj(b[ib1][ib2][1]);
					eps0s[5][n+1] += weight2 * b[ib1][ib2][1] * conj(b[ib1][ib2][2]);
					eps0s[6][n+1] += weight2 * b[ib1][ib2][2] * conj(b[ib1][ib2][0]);
					eps0s[7][n+1] += weight2 * b[ib1][ib2][2] * conj(b[ib1][ib2][1]);
					eps0s[8][n+1] += weight2 * b[ib1][ib2][2] * conj(b[ib1][ib2][2]);					
								
				}
				else if((delta_e > ((nomega-1) * domega) || delta_e == ((nomega-1) * domega)) && delta_e < (nomega * domega))
				{
					int n = int(delta_e/domega);
					double e1 = double(n) * domega;
					double e2 = double(n+1) * domega;					
					complex<double> weight1 = complex<double>(32 * (wf.wg(ik,ib1) - wf.wg(ik,oband+ib2)) * (e2 - delta_e)/domega/ucell.omega * PI * PI /delta_e/delta_e/domega, 0.0);

					eps0s[0][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][0]);
					eps0s[1][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][1]);
					eps0s[2][n] += weight1 * b[ib1][ib2][0] * conj(b[ib1][ib2][2]);
					eps0s[3][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][0]);
					eps0s[4][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][1]);
					eps0s[5][n] += weight1 * b[ib1][ib2][1] * conj(b[ib1][ib2][2]);
					eps0s[6][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][0]);
					eps0s[7][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][1]);
					eps0s[8][n] += weight1 * b[ib1][ib2][2] * conj(b[ib1][ib2][2]);					
					
				}
				
			}
	}
	
	return;
}

void Epsilon0_vasp:: Cal_T()
{
	complex<double> M1, M2;
	for(int n1=0; n1<nomega; n1++)
		for(int n=0; n<nomega; n++)
		{
			double n_e = double(n) * domega;
			double n1_e = double(n1) * domega;
			
			M1 = complex<double>(n1_e - n_e, -eta);
			M2 = complex<double>(-n1_e - n_e, eta);
			
			T[n][n1] = (1.0/M1 + 1.0/M2) * (-0.5) * 2.0/PI * domega;
		}
		
	return;
}

void Epsilon0_vasp:: Cal_epsilon0()
{
	CMatrixMul(9,nomega,nomega,eps0s,T,eps0);
	
	return;
}

void Epsilon0_vasp:: CMatrixMul(int m, int n, int l, complex<double>** A, complex<double>** B, complex<double>** C)
{
	for(int i=0; i<m; i++)
		for(int j=0; j<l; j++)
		{
			C[i][j] = complex<double>(0.0, 0.0);
			for(int k=0; k<n; k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
		
	return;
}
