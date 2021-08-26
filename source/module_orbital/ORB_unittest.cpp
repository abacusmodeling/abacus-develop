#include "ORB_unittest.h"
#include <fstream>
#include <iostream>
#include <string>
#include "../module_base/global_function.h"	

test_orb::test_orb()
{}

test_orb::~test_orb()
{}

void test_orb::set_files()
{
	std::cout << "read names of atomic basis set files" << std::endl;
	std::ifstream ifs("STRU",std::ios::in);

	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs,"NUMERICAL_ORBITAL");
	GlobalC::ORB.read_in_flag = true;

	for(int it=0;it<ntype;it++)
	{
		std::string ofile;
		ifs >> ofile;
		GlobalC::ORB.orbital_file.push_back(ofile);

		std::cout << "Numerical orbital file : " << ofile << std::endl;
	}
	
	return;
}

void test_orb::count_ntype()
{
	std::cout << "count number of atom types" << std::endl;
	std::ifstream ifs("STRU",std::ios::in);

	if (!ifs)
	{
		std::cout << "ERROR : file STRU does not exist" <<std::endl;
		exit(1);
	}

	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs,"ATOMIC_SPECIES");	
	
	ntype = 0;

	std::string x;
	ifs.rdstate();
	while ( ifs.good() )
	{
//read a line
		std::getline(ifs,x);

//trim white space
		const char* typeOfWhitespaces = " \t\n\r\f\v";
		x.erase(x.find_last_not_of(typeOfWhitespaces) + 1);
		x.erase(0,x.find_first_not_of(typeOfWhitespaces));

		if(x=="LATTICE_CONSTANT" || x=="NUMERICAL_ORBITAL" || x=="LATTICE_VECTORS" || x=="ATOMIC_POSITIONS") break;

		std::string tmpid=x.substr(0,1);
		if(!x.empty() && tmpid!="#") ntype++;
	}

	std::cout << "ntype : "<< ntype << std::endl;
	ifs.close();

	return;
}

