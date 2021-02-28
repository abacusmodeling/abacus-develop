#include"./stress_lcao.h"
#include "../src_pw/H_XC_pw.h"

void Stress_LCAO::start_stress
(
const double overlap[][3],
const double tvnl_dphi[][3],
const double vnl_dbeta[][3],
const double vl_dphi[][3], 
const matrix& stress_vdw,
matrix& scs
)
{
    TITLE("Stress_LCAO","start_stress");
	timer::tick("Stress_LCAO","start_stress",'E');

	matrix soverlap;
	soverlap.create(3,3);
	matrix stvnl_dphi;
	stvnl_dphi.create(3,3);
	matrix svnl_dbeta;
	svnl_dbeta.create(3,3);
	matrix svl_dphi;
	svl_dphi.create(3,3);
	matrix sigmacc;
	sigmacc.create(3,3);
	matrix sigmadvl;
	sigmadvl.create(3,3);
	matrix sigmaewa;
	sigmaewa.create(3,3);
	matrix sigmaxc;
	sigmaxc.create(3,3);
	matrix sigmahar;
	sigmahar.create(3,3);
	
    for(int i=0;i<3;i++)
	{
       for(int j=0;j<3;j++)
		{
          scs(i,j) = 0.0;
          soverlap(i,j) = overlap[i][j];
          stvnl_dphi(i,j) = tvnl_dphi[i][j];
          svnl_dbeta(i,j) = vnl_dbeta[i][j];
          svl_dphi(i,j) = vl_dphi[i][j];
          sigmacc(i,j) = 0.0;
          sigmadvl(i,j) = 0.0;
          sigmaewa(i,j) = 0.0;
          sigmaxc(i,j) = 0.0;
          sigmahar(i,j) = 0.0;
       }
    }
	//--------------------------------------------------------
	// local pseudopotential stress: 
	// use charge density; plane wave; local pseudopotential;
	//--------------------------------------------------------
    this->stress_loc (sigmadvl, 0);
 
	//--------------------------------------------------------
	//hartree term
	//--------------------------------------------------------
	this->stress_har (sigmahar, 0);
    
	//--------------------------------------------------------
	// ewald stress: use plane wave only.
	//--------------------------------------------------------
    this->stress_ewa (sigmaewa, 0); //remain problem


	//--------------------------------------------------------
	// stress due to core correlation. 
	//--------------------------------------------------------
	this->stress_cc(sigmacc, 0);

	//--------------------------------------------------------
	// stress due to self-consistent charge.
	//--------------------------------------------------------
	for(int i=0;i<3;i++)
	{
		sigmaxc(i,i) =  -(H_XC_pw::etxc) / ucell.omega;
	}
	//Exchange-correlation for PBE
	stress_gga(sigmaxc);

/*	if(vdwd2.vdwD2)									//Peize Lin add 2014-04-04, update 2019-04-26
	{
		vdwd2.stress();
	}
*/	
/*	if(vdwd3.vdwD3)									//jiyy add 2019-05-18
	{
		vdwd3.stress();
	}
*/	
/*	matrix sefield;
	if(EFIELD)
	{
		sefield.create(3, 3);
		Efield::compute_stress(sefield);
	}*/

	for(int i=0; i<3; i++)
	{
		for (int j=0;j<3;j++)
		{
			scs(i,j) += soverlap(i,j)
				+ stvnl_dphi(i,j) 
				+ svnl_dbeta(i,j) 
				+ svl_dphi(i,j) 
				+ sigmadvl(i,j) // derivative of local potential stress (pw)
				+ sigmaewa(i,j) // ewald stress (pw)
				+ sigmacc(i,j) //nonlinear core correction stress (pw)
				+ sigmaxc(i,j)//exchange corretion stress 
				+ sigmahar(i,j);// hartree stress 

			if(vdwd2.vdwD2)			// Peize Lin update 2019-04-26
			{
				scs(i,j) += stress_vdw(i , j);
			}
			if(vdwd3.vdwD3)			// jiyy add 2019-05-18
			{
				scs(i,j) += stress_vdw(i , j);
			}   
			/*		if(EFIELD)
					{
					scs(iat, i) = scs(iat, i) + sefield(iat, i);
					}
			 */

		}

		//	if(OUT_LEVEL != "m") ofs_running << " correction stress for each atom along direction " 
		// << i+1 << " is " << sum/ucell.nat << endl;
	}

	// test
	matrix svlocal;
	svlocal.create(3,3);
	for (int i = 0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			svlocal(i,j) = 0.0;
			svlocal(i,j) = svl_dphi(i,j) + sigmadvl(i,j);
		}
	}

	// test
	matrix stvnl;
	stvnl.create(3,3);
	for (int i = 0; i < 3; i++)
	{
		for(int j=0; j<3; j++)
		{
			stvnl(i,j) = 0.0;
			stvnl(i,j) = stvnl_dphi(i,j) + svnl_dbeta(i,j);
		}
	}

    if(Symmetry::symm_flag)
    {
		symm.stress_symmetry(scs);
    }//end symmetry

	// print Rydberg stress or not
	bool ry = false;
  
//  int TEST_STRESS = 1;
	if(TEST_STRESS)
	{
		ofs_running << "\n PARTS OF STRESS: " << endl;
		ofs_running << setiosflags(ios::showpos);
		ofs_running << setiosflags(ios::fixed) << setprecision(8) << endl;
		this->print_stress("OVERLAP    STRESS",soverlap,TEST_STRESS,ry);
		//test
		this->print_stress("T      STRESS",stvnl_dphi,TEST_STRESS,ry);
		this->print_stress("VNL      STRESS",svnl_dbeta,TEST_STRESS,ry);	

		this->print_stress("T_VNL      STRESS",stvnl,TEST_STRESS,ry);

		this->print_stress("VL_dPHI    STRESS",svl_dphi,TEST_STRESS,ry);
		this->print_stress("VL_dVL     STRESS",sigmadvl,TEST_STRESS,ry);
		this->print_stress("HAR     STRESS",sigmahar,TEST_STRESS,ry);

		this->print_stress("EWALD      STRESS",sigmaewa,TEST_STRESS,ry);
		this->print_stress("cc      STRESS",sigmacc,TEST_STRESS,ry);
		//		this->print_stress("NLCC       STRESS",sigmacc,TEST_STRESS,ry);
		this->print_stress("XC        STRESS",sigmaxc,TEST_STRESS,ry);
		this->print_stress("TOTAL        STRESS",scs,TEST_STRESS,ry);
	}


/*	if(EFIELD) 
	{
		STRESS::print("EFIELD     STRESS", sefield);
	}
*/
	ofs_running << setiosflags(ios::left);  
	
	this->printstress_total(scs, ry);

	
	timer::tick("Stress_LCAO","start_stress",'E');
	return;
}
