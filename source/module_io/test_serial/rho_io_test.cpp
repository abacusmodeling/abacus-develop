#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_io/rho_io.h"
#include "module_base/global_variable.h"
#include "prepare_unitcell.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
LCAO_Orbitals::LCAO_Orbitals(){}
LCAO_Orbitals::~LCAO_Orbitals(){}
#endif
Magnetism::Magnetism()
{
	this->tot_magnetization = 0.0;
	this->abs_magnetization = 0.0;
	this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
	delete[] this->start_magnetization;
}

/************************************************
 *  unit test of read_rho and write_rho
 ***********************************************/

/**
 * - Tested Functions:
 *   - read_rho()
 *     - the function to read_rho from file
 *     - the serial version without MPI
 *   - write_rho()
 *     - the function to write_rho to file
 *     - the serial version without MPI
 */

class RhoIOTest : public ::testing::Test
{
protected:
	int nspin = 1;
	int nrxx = 36*36*36;
	int prenspin = 1;
	double **rho;
	UnitCell* ucell;
	void SetUp()
	{
		rho = new double*[nspin];
		ucell = new UnitCell;
		for(int is=0; is<nspin; ++is)
		{
			rho[is] = new double[nrxx];
		}
	}
	void TearDown()
	{
		for(int is=0; is<nspin; ++is)
		{
			delete[] rho[is];
		}
		delete[] rho;
		delete ucell;
	}
};

TEST_F(RhoIOTest,Read)
{
	int is = 0;
	std::string fn = "./support/SPIN1_CHG.cube";
	int nx = 36;
	int ny = 36;
	int nz = 36;
	double ef;
	UcellTestPrepare utp = UcellTestLib["Si"];
	ucell = utp.SetUcellInfo();
	ModuleIO::read_rho(is,nspin,fn,rho[is],nx,ny,nz,ef,ucell,prenspin);
	EXPECT_DOUBLE_EQ(ef,0.461002);
	EXPECT_DOUBLE_EQ(rho[0][0],1.27020863940e-03);
	EXPECT_DOUBLE_EQ(rho[0][46655],1.33581335706e-02);
}

TEST_F(RhoIOTest,Write)
{
	int is = 0;
	std::string fn = "./support/SPIN1_CHG.cube";
	int nx = 36;
	int ny = 36;
	int nz = 36;
	double ef;
	UcellTestPrepare utp = UcellTestLib["Si"];
	ucell = utp.SetUcellInfo();
	//first read
	ModuleIO::read_rho(is,nspin,fn,rho[is],nx,ny,nz,ef,ucell,prenspin);
	EXPECT_DOUBLE_EQ(ef,0.461002);
	EXPECT_DOUBLE_EQ(rho[0][0],1.27020863940e-03);
	EXPECT_DOUBLE_EQ(rho[0][46655],1.33581335706e-02);
	//then write
	std::string ssc = "SPIN1_CHG.cube";
	GlobalV::MY_RANK = 0;
	GlobalV::out_chg = 1;
	ModuleIO::write_rho(rho[is],is,nspin,0,ssc,nx,ny,nz,ef,ucell);
	std::ifstream ifs;
	ifs.open("SPIN1_CHG.cube");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("1 (nspin) 0.461002 (fermi energy, in Ry)"));
	ifs.close();
	remove("SPIN1_CHG.cube");
}
