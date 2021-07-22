#include "../src_pw/charge.h"
#include "../src_pw/global.h"

bool Charge::read_rho(const int &is, const string &fn, double* rho) //add by dwan
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

	if(NSPIN != 4) CHECK_INT(ifs, NSPIN);
	else
	{
		READ_VALUE(ifs, PRENSPIN);
	}
	if(NSPIN == 1||NSPIN == 4)
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
				ifs >> rho[i*pw.ncy*pw.ncz + j*pw.ncz +k];
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
		Pgrid.zpiece_to_all(zpiece, iz, rho);
	}// iz
	delete[] zpiece;
#endif

    if(MY_RANK==0) ifs.close();
    return true;
}
