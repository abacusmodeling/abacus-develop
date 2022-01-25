#include "./stress_func.h"
#include "./xc_functional.h"
#include "./myfunc.h"
#include "./xc_gga_pw.h"
#include "../module_base/timer.h"

//calculate the GGA stress correction in PW and LCAO
void Stress_Func::stress_gga(ModuleBase::matrix& sigma) 
{
	ModuleBase::timer::tick("Stress_Func","stress_gga");
     
	if (GlobalC::xcf.igcx == 0  &&  GlobalC::xcf.igcc == 0)
	{
		ModuleBase::timer::tick("Stress_Func","stress_gga");
		return;
	} 
	double sigma_gradcorr[3][3];
	double* p= &sigma_gradcorr[0][0];
	for(int i=0;i<9;i++)
		*p++ = 0;

	bool igcc_is_lyp = false;
	if( GlobalC::xcf.igcc == 3 || GlobalC::xcf.igcc == 7)
	{
		igcc_is_lyp = true;
	}
	const int nspin_in = GlobalV::NSPIN;
	assert(nspin_in>0);
	const double fac = 1.0/ nspin_in;

	// doing FFT to get rho in G space: rhog1 
	GlobalC::CHR.set_rhog(GlobalC::CHR.rho[0], GlobalC::CHR.rhog[0]);
	if(nspin_in==2)//mohan fix bug 2012-05-28
	{
		GlobalC::CHR.set_rhog(GlobalC::CHR.rho[1], GlobalC::CHR.rhog[1]);
	}
	GlobalC::CHR.set_rhog(GlobalC::CHR.rho_core, GlobalC::CHR.rhog_core);
		
	double* rhotmp1;
	double* rhotmp2;
	std::complex<double>* rhogsum1;
	std::complex<double>* rhogsum2;
	ModuleBase::Vector3<double>* gdr1;
	ModuleBase::Vector3<double>* gdr2;
 
	rhotmp1 = new double[GlobalC::pw.nrxx];
	rhogsum1 = new std::complex<double>[GlobalC::pw.ngmc];
	ModuleBase::GlobalFunc::ZEROS(rhotmp1, GlobalC::pw.nrxx);
	ModuleBase::GlobalFunc::ZEROS(rhogsum1, GlobalC::pw.ngmc);
	for(int ir=0; ir<GlobalC::pw.nrxx; ir++) rhotmp1[ir] = GlobalC::CHR.rho[0][ir] + fac * GlobalC::CHR.rho_core[ir];
	for(int ig=0; ig<GlobalC::pw.ngmc; ig++) rhogsum1[ig] = GlobalC::CHR.rhog[0][ig] + fac * GlobalC::CHR.rhog_core[ig];
	gdr1 = new ModuleBase::Vector3<double>[GlobalC::pw.nrxx];
	ModuleBase::GlobalFunc::ZEROS(gdr1, GlobalC::pw.nrxx);

	GGA_PW::grad_rho( rhogsum1 , gdr1 );

	if(nspin_in==2)
	{
		rhotmp2 = new double[GlobalC::pw.nrxx];
		rhogsum2 = new std::complex<double>[GlobalC::pw.ngmc];
		ModuleBase::GlobalFunc::ZEROS(rhotmp2, GlobalC::pw.nrxx);
		ModuleBase::GlobalFunc::ZEROS(rhogsum2, GlobalC::pw.ngmc);
		for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			rhotmp2[ir] = GlobalC::CHR.rho[1][ir] + fac * GlobalC::CHR.rho_core[ir];
		}
		for(int ig=0; ig<GlobalC::pw.ngmc; ig++)
		{
			rhogsum2[ig] = GlobalC::CHR.rhog[1][ig] + fac * GlobalC::CHR.rhog_core[ig];
		}
		
		gdr2 = new ModuleBase::Vector3<double>[GlobalC::pw.nrxx];
		ModuleBase::GlobalFunc::ZEROS(gdr2, GlobalC::pw.nrxx);

		GGA_PW::grad_rho( rhogsum2 , gdr2 );
	}
        
	const double epsr = 1.0e-6;
	const double epsg = 1.0e-10;

	double grho2a = 0.0;
	double grho2b = 0.0;
	double atau = 0.0;
	double sx = 0.0;
	double sc = 0.0;
	double v1x = 0.0;
	double v2x = 0.0;
	double v3x = 0.0;
	double v1c = 0.0;
	double v2c = 0.0;
	double v3c = 0.0;

	if(nspin_in==1||nspin_in==4)
	{
		//double segno;
		for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			const double arho = std::abs( rhotmp1[ir] );
			if(arho > epsr)
			{
				grho2a = gdr1[ir].norm2();
				if( grho2a > epsg )
				{
					//if( rhotmp1[ir] >= 0.0 ) segno = 1.0;
					//if( rhotmp1[ir] < 0.0 ) segno = -1.0;
					if(GlobalV::DFT_META)
					{
#ifdef USE_LIBXC
						atau = GlobalC::CHR.kin_r[0][ir]/2.0;
						XC_Functional::tau_xc( arho, grho2a, atau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c);
#endif
					}
					else
					{
						XC_Functional::gcxc( arho, grho2a, sx, sc, v1x, v2x, v1c, v2c);
					}
					double tt[3];
					tt[0] = gdr1[ir].x;
					tt[1] = gdr1[ir].y;
					tt[2] = gdr1[ir].z;
					for(int l = 0;l< 3;l++)
					{
						for(int m = 0;m< l+1;m++)
						{
							sigma_gradcorr[l][m] += tt[l] * tt[m] * ModuleBase::e2 * (v2x + v2c);
						}
					}
				}
			}
		} 
	}
	else if(nspin_in==2)
	{
		double v1cup = 0.0;
		double v1cdw = 0.0;
		double v2cup = 0.0;
		double v2cdw = 0.0;
		double v1xup = 0.0;
		double v1xdw = 0.0;
		double v2xup = 0.0;
		double v2xdw = 0.0;
		double v2cud = 0.0;
		double v2c = 0.0;
		for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			if(GlobalV::DFT_META)
			{
				ModuleBase::WARNING_QUIT("stress_gga","stress mGGA not ready for nspin=2");
			}
			double rh = rhotmp1[ir] + rhotmp2[ir];
			grho2a = gdr1[ir].norm2();;
			grho2b = gdr2[ir].norm2();;
			//XC_Functional::gcx_spin();
			gcx_spin(rhotmp1[ir], rhotmp2[ir], grho2a, grho2b,
				sx, v1xup, v1xdw, v2xup, v2xdw);

			if(rh > epsr)
			{
				if(igcc_is_lyp)
				{
					ModuleBase::WARNING_QUIT("stress","igcc_is_lyp is not available now.");
				}
				else
				{
					double zeta = ( rhotmp1[ir] - rhotmp2[ir] ) / rh;
					double grh2 = (gdr1[ir]+gdr2[ir]).norm2();
					//XC_Functional::gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
					gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
					v2cup = v2c;
					v2cdw = v2c;
					v2cud = v2c;
				}
			}
			else
			{
				sc = 0.0;
				v1cup = 0.0;
				v1cdw = 0.0;
				v2c = 0.0;
				v2cup = 0.0;
				v2cdw = 0.0;
				v2cud = 0.0;
			}
			double tt1[3],tt2[3];
			{
				tt1[0] = gdr1[ir].x;
				tt1[1] = gdr1[ir].y;
				tt1[2] = gdr1[ir].z;
				tt2[0] = gdr2[ir].x;
				tt2[1] = gdr2[ir].y;
				tt2[2] = gdr2[ir].z;
			}
			for(int l = 0;l< 3;l++)
			{
			    for(int m = 0;m< l+1;m++)
				{
				//    exchange
				sigma_gradcorr [l][m] += tt1[l] * tt1[m] * ModuleBase::e2 * v2xup + 
							tt2[l] * tt2[m] * ModuleBase::e2 * v2xdw;
				//    correlation
				sigma_gradcorr [l][m] += ( tt1[l] * tt1[m] * v2cup + 
							tt2[l] * tt2[m] * v2cdw + 
							(tt1[l] * tt2[m] +
							tt2[l] * tt1[m] ) * v2cud ) * ModuleBase::e2;
				}
			}
		}
	}

	for(int l = 0;l< 3;l++)
	{
		for(int m = 0;m< l;m++)
		{
			sigma_gradcorr[m][l] = sigma_gradcorr[l][m];
		}
	}
	for(int l = 0;l<3;l++)
	{
		for(int m = 0;m<3;m++)
		{
			Parallel_Reduce::reduce_double_pool( sigma_gradcorr[l][m] );
		}
	}
	
/*	p= &sigma_gradcorr[0][0];
	double* p1 = &sigmaxc[0][0];
	for(int i=0;i<9;i++){
		*p /= GlobalC::pw.ncxyz ;
		*p1++ += *p++;  
	}*/
	
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			sigma(i,j) += sigma_gradcorr[i][j] / GlobalC::pw.ncxyz;
		}
	}

	delete[] rhotmp1;
	delete[] rhogsum1;
	delete[] gdr1;
	if(nspin_in==2)
	{
		delete[] rhotmp2;
		delete[] rhogsum2;
		delete[] gdr2;
	}
	ModuleBase::timer::tick("Stress_Func","stress_gga");
	return;
}
