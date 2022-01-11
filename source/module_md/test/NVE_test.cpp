#include "setcell.cpp"
#include "gtest/gtest.h"
#include "../NVE.h"

UnitCell_pseudo ucell;
class NVE_test : public testing::Test
{
protected:
	Verlet *verlet;

	static void SetUpTestCase()
	{
		setupcell(ucell);
        parameters();
	}

	static void TearDownTestCase(){}
};

TEST_F(NVE_test, setup)
{
    verlet = new NVE(INPUT.mdp, ucell);
	verlet->setup();
    
    EXPECT_DOUBLE_EQ(verlet->temperature_, 299.99999999999665);
    EXPECT_DOUBLE_EQ(verlet->stress(0,0), 6.0100555286436806e-06);
	EXPECT_DOUBLE_EQ(verlet->stress(0,1), -1.4746713013791574e-06);
	EXPECT_DOUBLE_EQ(verlet->stress(0,2), 1.5039983732220751e-06);
	EXPECT_DOUBLE_EQ(verlet->stress(1,0), -1.4746713013791574e-06);
	EXPECT_DOUBLE_EQ(verlet->stress(1,1), 3.4437172989317909e-06);
	EXPECT_DOUBLE_EQ(verlet->stress(1,2), -1.251414906590483e-06);
	EXPECT_DOUBLE_EQ(verlet->stress(2,0), 1.5039983732220751e-06);
	EXPECT_DOUBLE_EQ(verlet->stress(2,1), -1.251414906590483e-06);
	EXPECT_DOUBLE_EQ(verlet->stress(2,2), 1.6060561926126463e-06);
}

TEST_F(NVE_test, first_half)
{
    verlet->first_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.9945454470992772);
	EXPECT_DOUBLE_EQ(verlet->pos[0].y, 0.0029590658162135359);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 9.9994204767196599);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2063192793031234);
	EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1939342598421803);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 0.0039873402383468889);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0944648036273872);
	EXPECT_DOUBLE_EQ(verlet->pos[2].y, 9.9998836061438716);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.997763438399228);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.0046704699702541427);
	EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.303223068197739);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9988287446427613);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00013193932519649473);
	EXPECT_DOUBLE_EQ(verlet->vel[0].y, 7.1576379239356465e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -1.40179977966e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, 0.00015285605661129458);
	EXPECT_DOUBLE_EQ(verlet->vel[1].y, -0.00014672323796402785);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, 9.6449148069800003e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -0.00013388999749840003);
	EXPECT_DOUBLE_EQ(verlet->vel[2].y, -2.8154327428808153e-06);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, -5.4099838013700003e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 0.0001129732660846002);
	EXPECT_DOUBLE_EQ(verlet->vel[3].y, 7.7962291467652202e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -2.83313122596e-05);
}

TEST_F(NVE_test, second_half)
{
    verlet->second_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.9945454470992772);
	EXPECT_DOUBLE_EQ(verlet->pos[0].y, 0.0029590658162135359);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 9.9994204767196599);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2063192793031234);
	EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1939342598421803);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 0.0039873402383468889);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0944648036273872);
	EXPECT_DOUBLE_EQ(verlet->pos[2].y, 9.9998836061438716);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.997763438399228);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.0046704699702541427);
	EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.303223068197739);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9988287446427613);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00013179791402898947);
	EXPECT_DOUBLE_EQ(verlet->vel[0].y, 7.1809815495212933e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -1.40179977966e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, 0.00015267223469058917);
	EXPECT_DOUBLE_EQ(verlet->vel[1].y, -0.00014691320932005571);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, 9.6449148069800003e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -0.00013399051477080008);
	EXPECT_DOUBLE_EQ(verlet->vel[2].y, -2.5857616745616307e-06);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, -5.4099838013700003e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 0.0001131161941102004);
	EXPECT_DOUBLE_EQ(verlet->vel[3].y, 7.7689155499504408e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -2.83313122596e-05);

    delete verlet;
}

int main(int argc, char **argv) 
{
    srand(2000);
#ifdef __MPI
	MPI_Init(&argc,&argv);
#endif

    testing::InitGoogleTest(&argc, argv);
	int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}