#include "../xc_functional.h"
#include "gtest/gtest.h"

//a mock function of WARNING_QUIT, to avoid the uncorrected call by matrix.cpp at line 37.
namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {return ;}
}

namespace GlobalV
{
    std::string BASIS_TYPE = "";
    bool STRESS = 0;
    int NSPIN = 1;
}

class XCTest_PBE : public testing::Test
{
    protected:
        void SetUp()
        {
            XC_Functional::set_xc_type("PBE");
            std::cout << XC_Functional::get_func_type() << std::endl;
        }
};

TEST_F(XCTest_PBE, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
}