#include "setcell.cpp"
#include "gtest/gtest.h"
#include "../Langevin.h"

UnitCell_pseudo ucell;
class Langevin_test : public testing::Test
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

TEST_F(Langevin_test, setup)
{
    verlet = new Langevin(INPUT.mdp, ucell);
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

    EXPECT_DOUBLE_EQ(verlet->force[0].x, -0.023875294316943713);
	EXPECT_DOUBLE_EQ(verlet->force[0].y, -0.48819483738869529);
	EXPECT_DOUBLE_EQ(verlet->force[0].z, 0.29384542562861765);
	EXPECT_DOUBLE_EQ(verlet->force[1].x, -0.5737192187058936);
	EXPECT_DOUBLE_EQ(verlet->force[1].y, 0.29063407548222864);
	EXPECT_DOUBLE_EQ(verlet->force[1].z, -0.54217115712378094);
	EXPECT_DOUBLE_EQ(verlet->force[2].x, 0.24793819074617376);
	EXPECT_DOUBLE_EQ(verlet->force[2].y, 0.026039528498451155);
	EXPECT_DOUBLE_EQ(verlet->force[2].z, 0.54931240592385677);
	EXPECT_DOUBLE_EQ(verlet->force[3].x, -0.38147565167794611);
	EXPECT_DOUBLE_EQ(verlet->force[3].y, -0.23254957368548659);
	EXPECT_DOUBLE_EQ(verlet->force[3].z, -0.16806427435447177);
}

TEST_F(Langevin_test, first_half)
{
    verlet->first_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.994259423387355);
	EXPECT_DOUBLE_EQ(verlet->pos[0].y, 9.9972204283509285);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 0.0028687651122257405);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.1995942596625513);
	EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1973527165828584);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 9.9976249389728018);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0973785242155092);
	EXPECT_DOUBLE_EQ(verlet->pos[2].y, 0.00017968618542555382);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 5.0042096425189211);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.00018792835923789908);
	EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.3005053810142835);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9968565033352039);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00013885790793427856);
	EXPECT_DOUBLE_EQ(verlet->vel[0].y, -6.7234622963379829e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, 6.9392109663872021e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, -9.8143894288626359e-06);
	EXPECT_DOUBLE_EQ(verlet->vel[1].y, -6.40347236553997e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, -5.7449978931817646e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -6.3410466874160769e-05);
	EXPECT_DOUBLE_EQ(verlet->vel[2].y, 4.346400976153433e-06);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, 0.00010182638309204409);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 4.5457696266635515e-06);
	EXPECT_DOUBLE_EQ(verlet->vel[3].y, 1.2224582143593012e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -7.6037548128964855e-05);
}

TEST_F(Langevin_test, second_half)
{
    verlet->second_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.994259423387355);
	EXPECT_DOUBLE_EQ(verlet->pos[0].y, 9.9972204283509285);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 0.0028687651122257405);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.1995942596625513);
	EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1973527165828584);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 9.9976249389728018);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0973785242155092);
	EXPECT_DOUBLE_EQ(verlet->pos[2].y, 0.00017968618542555382);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 5.0042096425189211);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.00018792835923789908);
	EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.3005053810142835);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9968565033352039);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00019020019349660282);
	EXPECT_DOUBLE_EQ(verlet->vel[0].y, -0.00030594937616533209);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, 0.00017443321209239919);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, -9.7187804913322887e-05);
	EXPECT_DOUBLE_EQ(verlet->vel[1].y, 0.0001082901135078902);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, -0.00031318370325007474);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -9.9014035573776046e-05);
	EXPECT_DOUBLE_EQ(verlet->vel[2].y, -7.279491639258441e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, 0.00013041787625213766);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, -0.00021454814645268522);
	EXPECT_DOUBLE_EQ(verlet->vel[3].y, -0.0001458202419857561);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -1.6354503527938183e-05);

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