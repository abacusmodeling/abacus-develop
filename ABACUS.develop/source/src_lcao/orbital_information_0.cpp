//=========================================================
//AUTHOR : mohan
//DATE : 2008-11-29
//=========================================================
#include "orbital_information_0.h"

Orbital_Information_0::Orbital_Information_0()
{
}

Orbital_Information_0::~Orbital_Information_0()
{
	if (init_done)
	{
		delete[] label;
		delete[] enviroment;
		delete[] transparence;
		delete[] lmax;

		for (int i = 0; i < ntype; i++)
		{
			delete[] nchi[i];
		}

		delete[] nchi;
	}
}

void Orbital_Information_0::Read(const string &fn, ofstream &ofs)
{
	if (test_overlap)TITLE("Orbital_Information_0", "Read");

	ifstream ifs(fn.c_str());

	if (!ofs)
	{
		WARNING_QUIT("Orbital_Information_0::Read","Can not find the file.");
	}
	else
	{
		cout << "\n Read in file : " << fn;
	}

	CHECK_NAME(ifs, "INPUT_ORBITAL_INFORMATION_0");

	CHECK_NAME(ifs, "ntype");
	ifs >> ntype;
	OUT(ofs, "ntype", ntype);

	if (ntype <= 0 || ntype > 50)
	{
		WARNING_QUIT("Orbital_Information_0::Read", "ntype<=0 || ntype>50");
	}

	this->init_done = false;

	this->label = new string[ntype];
	this->enviroment  = new string[ntype];
	this->transparence = new double[ntype];
	this->lmax = new int[ntype];
	this->nchi = new int*[ntype];

	for (int i = 0; i < ntype; i++)
	{
		CHECK_NAME(ifs, "label");
		ifs >> label[i];
		OUT(ofs, "label", label[i]);

		CHECK_NAME(ifs, "enviroment");
		ifs >> enviroment[i];
		OUT(ofs, "enviroment", enviroment[i]);

		CHECK_NAME(ifs, "transparence");
		ifs >> transparence[i];
		OUT(ofs, "transparence", transparence[i]);

		CHECK_NAME(ifs, "Lmax");
		ifs >> lmax[i];
		OUT(ofs, "Lmax", lmax[i]);
		nchi[i] = new int[ lmax[i] ];

		for (int j = 0;j <= lmax[i];j++)
		{
			CHECK_NAME(ifs, "L");
			OUT(ofs, "L", j);
			int k;
			ifs >> k;
			assert(k == j);

			CHECK_NAME(ifs, "nchi");
			int m;
			ifs >> nchi[i][j];
			OUT(ofs, "nchi", nchi[i][j]);
		}
	}

	this->init_done = true;

	return;
}
