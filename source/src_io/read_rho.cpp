#include "../src_pw/charge.h"
#include "../src_pw/global.h"

bool Charge::read_rho(const int &is, const string &fn, double* rho) //add by dwan
{
    TITLE("Charge","read_rho");
    ifstream ifs(fn.c_str());
    if (!ifs) 
	{
		GlobalV::ofs_running << " !!! Couldn't find the charge file !!!" << endl;
		return false;
	}
	else
	{
    	GlobalV::ofs_running << " Find the file, try to read charge from file." << endl;
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

	if(GlobalV::NSPIN != 4) CHECK_INT(ifs, GlobalV::NSPIN);
	else
	{
		READ_VALUE(ifs, GlobalV::PRENSPIN);
	}
	if(GlobalV::NSPIN == 1||GlobalV::NSPIN == 4)
	{
		READ_VALUE(ifs, GlobalC::en.ef);
		GlobalV::ofs_running << " read in fermi energy = " << GlobalC::en.ef << endl;
	}
	else if(GlobalV::NSPIN == 2)
	{
		if(is==0)READ_VALUE(ifs, GlobalC::en.ef_up);
		else if(is==1)READ_VALUE(ifs, GlobalC::en.ef_dw);
	}
	else 
	{
		WARNING_QUIT("read_rho","check nspin!");
	}
	CHECK_INT(ifs, GlobalC::pw.ncx);	
	CHECK_INT(ifs, GlobalC::pw.ncy);	
	CHECK_INT(ifs, GlobalC::pw.ncz);	

#ifndef __MPI
	GlobalV::ofs_running << " Read SPIN = " << is+1 << " charge now." << endl;
	for(int k=0; k<GlobalC::pw.ncz; k++)
	{
		// consistent with the write_rho, something is
		// wrong.... but it works now.
		for(int j=0; j<GlobalC::pw.ncy; j++)
		{
			for(int i=0; i<GlobalC::pw.ncx; i++)
			{
				ifs >> rho[i*GlobalC::pw.ncy*GlobalC::pw.ncz + j*GlobalC::pw.ncz +k];
			}
		}
	}
#else
	
	const int nxy = GlobalC::pw.ncx * GlobalC::pw.ncy;
	double *zpiece = new double[nxy];
	for(int iz=0; iz<GlobalC::pw.ncz; iz++)
	{
		ZEROS(zpiece, nxy);
		if(GlobalV::MY_RANK==0)
		{
			//				GlobalV::ofs_running << " Read charge density iz=" << iz << endl;
			for(int j=0; j<GlobalC::pw.ncy; j++)
			{
				for(int i=0; i<GlobalC::pw.ncx; i++)
				{
					ifs >> zpiece[ i*GlobalC::pw.ncy + j ];
				}
			}
		}
		Pgrid.zpiece_to_all(zpiece, iz, rho);
	}// iz
	delete[] zpiece;
#endif

    if(GlobalV::MY_RANK==0) ifs.close();
    return true;
}
