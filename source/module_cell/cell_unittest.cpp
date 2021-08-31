#include "cell_unittest.h"

#ifdef __LCAO
test_cell_orb::test_cell_orb()
{}
test_cell_orb::~test_cell_orb()
{}
#else
test_cell::test_cell()
{}
test_cell::~test_cell()
{}
#endif

#ifdef __LCAO
void test_cell_orb::set_parameters()
#else
void test_cell::set_parameters()
#endif
{
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::global_pseudo_type= "auto";
	GlobalV::PSEUDORCUT = 15.0;
	GlobalV::global_out_dir="./";
	GlobalV::ofs_warning.open("warning.log");
	ofs_running.open("log.txt");

	ucell.latName = "test";
	ucell.ntype = ntype;
	return;
}

#ifdef __LCAO
void test_cell_orb::setup_cell()
#else
void test_cell::setup_cell()
#endif
{
	ucell.setup_cell(
#ifdef __LCAO
	ORB,
#endif
	"./",
	out,
	"STRU", 
	ofs_running);

	return;
}

#ifndef __LCAO
void test_cell::count_ntype()
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
#endif
