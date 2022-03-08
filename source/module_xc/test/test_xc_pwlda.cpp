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

class XCTest_PWLDA : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;

        void SetUp()
        {
            XC_Functional::set_xc_type("PWLDA");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};

            for(int i=0;i<5;i++)
            {
                double e,v;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
            }
        }
};

TEST_F(XCTest_PWLDA, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),1);
    std::vector<double> e_lda_ref  = {-0.9570906378,-0.9570906378,-0.9200171474,-0.3808291731,-9.124862983};
    std::vector<double> v_lda_ref  = {-1.259358955 , -1.259358955, -1.210234884,-0.4975768336,-12.12952188};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
    }
}