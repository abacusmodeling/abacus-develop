#include "lcao_vna.h"
#include "../src_pw/global.h"

LCAO_Vna::LCAO_Vna()
{}

LCAO_Vna::~LCAO_Vna()
{}

void LCAO_Vna::two_center_vna(void)
{
	cout << " two center integration for vna" << endl;

	WARNING_QUIT("LCAO_Vna::two_center_vna","not ready yet");
	return;
}


// be called in Use_Hamilt_Matrix::set_ion
void LCAO_Vna::dense_vna(const char &matrix_type)
{
	TITLE("LCAO_Vna","dense_vna");
	timer::tick("LCAO_Vna","dense_vna",'E');

	ofs_running << "\n SETUP VNA." << endl;

	const int istep = 1;
	const int delta_vh = 0;
	const int cal_vna = 1; // tmp by mohan

	const int dense = VNA;

	//--------------------------------------
	// set up a dense grid for doing dense
	// grid integration
	// note that for Si the length of Vna
	// may be as long as 10.0 Bohr, so if
	// cal_vna is turned on, the data in
	// Grid_Technique will be different.

	// for example, bx=2, by=2, bz=2
	// nbx, nby, nbz stay as unchanged.
	//--------------------------------------
	Grid_Technique gtf;
	gtf.set_pbc_grid(
	pw.ncx*dense,pw.ncy*dense,pw.ncz*dense,
	pw.bx*dense,pw.by*dense,pw.bz*dense,
	pw.nbx,pw.nby,pw.nbz,
	pw.nbxx,pw.nbzp_start,pw.nbzp,
	cal_vna);
	//0);

	//---------------------------------------------
	// delta_vh = 0 means don't need to calculate
	// delta Vhartree.
	// but calculate the neutral potential
	//---------------------------------------------

	//--------------------------------------------------------
	// The neutral potential can be only used when vna = 1,
	// then to be used to check if the real space neutral
	// potential is correct
	//--------------------------------------------------------

	// for test
	bool test_G_vl = false;

	// TEST !
	if(test_G_vl)
	{
		pot.init_pot(istep, delta_vh, cal_vna);
		for(int ir=0; ir<pw.nrxx; ir++)
		{
			pot.vrs1[ir] = pot.vrs(CURRENT_SPIN, ir);
		}
	}
	
	if(GAMMA_ONLY_LOCAL)
	{
		// pot.vrs1 is not used any more, but it can
		// exist for test.
	    UHM.GG.cal_vna_d(gtf,pot.vrs1,matrix_type);
		// for test
		//UHM.GG.cal_vna(pot.vrs1);
	}
	else
	{
		UHM.GK.cal_vna_d(gtf,pot.vrs1,matrix_type);
		// for test
	//	UHM.GK.cal_vna(pot.vrs1, gtf);
	}
	timer::tick("LCAO_Vna","dense_vna",'E');
	return;
}


// if VNA>1, do the local part except the Vna.
void LCAO_Vna::smooth_vl2(void)
{
	TITLE("LCAO_Vna","smooth_vl2");
	timer::tick("LCAO_Vna","smooth_vl2",'I');


	//cout << " calculate the <phi|delta_Vh+Vxc|phi>" << endl;
	const int istep = 1;
	// calculate the (delta_Vh + Vxc + Vefield)
	bool delta_vh = 1;
	bool vna = 0;
	//--------------------------------------------
	// NOTE1 : istep must set to 1, otherwise the 
	// charge density may be changed!
	// NOTE2 : the Vna should include the efield
	// contribution. The unsymmetry potential
	// due to efield can be only evaluated 
	// using grid integration.
	//--------------------------------------------
	pot.init_pot(istep, delta_vh, vna);

	for(int ir=0; ir<pw.nrxx; ir++)
	{
		pot.vrs1[ir] = pot.vrs(CURRENT_SPIN, ir);
	}

//	cout << " delta_Vh+Vxc is not calculated." << endl;
	if(GAMMA_ONLY_LOCAL)
	{
		if(GRID_SPEED==1)
		{
			UHM.GG.cal_vlocal(pot.vrs1);
		}
		else if(GRID_SPEED==2)
		{
			UHM.GS.cal_vlocal(pot.vrs1);
		}
	}
	else
	{
		UHM.GK.cal_vlocal_k(pot.vrs1, GridT);
	}

// for check, in fact, vna only need to do
// once for each ion iteration.
//	char matrix_type = 'L';
//	this->dense_vna(matrix_type);	

	//-------------------------------------
	// very important step!
	// turn back the full local potential.
	//-------------------------------------
	pot.init_pot(istep);


	timer::tick("LCAO_Vna","smooth_vl2",'I');
	return;
}



void LCAO_Vna::smooth_vl1(void)
{
	TITLE("LCAO_Vna","smooth_vl1");

	//cout << " calculate the <phi|delta_Vh+Vxc|phi>" << endl;
	const int istep=1;
	// calculate the (delta_Vh + Vxc)
	bool delta_vh = 1;
	bool vna = 0;
	// istep must set to 1, otherwise the charge
	// density may be changed!
	pot.init_pot(istep, delta_vh, vna);

	for(int ir=0; ir<pw.nrxx; ir++)
	{
		pot.vrs1[ir] = pot.vrs(CURRENT_SPIN, ir);
	}

	//cout << " delta_Vh+Vxc is not calculated." << endl;

	if(GAMMA_ONLY_LOCAL)
	{
		UHM.GG.cal_vlocal(pot.vrs1);
	}
	else
	{
		UHM.GK.cal_vlocal_k(pot.vrs1, GridT);
	}

	// calculate the vna part.
	//cout << " calculate the <phi|Vna|phi>" << endl;
	/*
	delta_vh = 0;
	vna = 1;
	pot.init_pot(istep, delta_vh, vna);

	for(int ir=0; ir<pw.nrxx; ir++)
	{
		pot.vrs1[ir] = pot.vrs(CURRENT_SPIN, ir);
	}

	UHM.GG.cal_vna(pot.vrs1);
	*/

	pot.init_pot(istep);

	return;
}
