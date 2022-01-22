#include "../tool_check.h"
#include "gtest/gtest.h"
#include <fstream>
#include <cstdio>

/************************************************
*  unit test of class ToolCheckTest
***********************************************/

/**
 * - Tested Function
 *   - ModuleBase::CHECK_NAME
 *	- check the next input from ifs is std::string
 *
 *   - ModuleBase::CHECK_INT
 *	- check the next input from ifs is int
 *
 *   - ModuleBase::CHECK_DOUBLE
 *	- check the next input from ifs is double
 *
 *   - ModuleBase::CHECK_STRING
 *	- check the next input from ifs is string
 *
 */


class ToolCheckTest : public testing::Test
{
    protected:
        std::ofstream ofs;
	std::ifstream ifs;
	// define std::string, int, double variables
	std::string name="abaqus";
	int ecut = 100;
	double occupation = 0.23;
	std::string caltype="nscf";
	// quit is the swith to control performance of function
	bool quit=false;
};

TEST_F(ToolCheckTest,Name){
	ofs.open("tmp");
	// double input to check continus check function
	ofs << name << std::endl;
	ofs << name << std::endl;
	ofs.close();
	ifs.open("tmp");
	// non-quit check
        ModuleBase::CHECK_NAME(ifs,"abacus",quit);
	// quit check
	EXPECT_DEATH(ModuleBase::CHECK_NAME(ifs,"abacus"),"");
	ifs.close();
}

TEST_F(ToolCheckTest,Int){
	ofs.open("tmp");
	// double input to check continus check function
	ofs << ecut << std::endl;
	ofs << ecut << std::endl;
	ofs.close();
	ifs.open("tmp");
	// non-quit check
        ModuleBase::CHECK_INT(ifs,80,quit);
	// quit check
	EXPECT_DEATH(ModuleBase::CHECK_INT(ifs,80),"");
	ifs.close();
}

TEST_F(ToolCheckTest,Double){
	ofs.open("tmp");
	// double input to check continus check function
	ofs << occupation << std::endl;
	ofs << occupation << std::endl;
	ofs.close();
	ifs.open("tmp");
	// non-quit check
        ModuleBase::CHECK_DOUBLE(ifs,0.23002,quit);
	// quit check
	EXPECT_DEATH(ModuleBase::CHECK_DOUBLE(ifs,0.23002),"");
	ifs.close();
}

TEST_F(ToolCheckTest,String){
	ofs.open("tmp");
	// double input to check continus check function
	ofs << caltype << std::endl;
	ofs << caltype << std::endl;
	ofs.close();
	ifs.open("tmp");
	// non-quit check
        ModuleBase::CHECK_STRING(ifs,"scf",quit);
	// quit check
	EXPECT_DEATH(ModuleBase::CHECK_STRING(ifs,"scf"),"");
	ifs.close();
	// remove file tmp
	remove("tmp");
}

