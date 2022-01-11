#include "gtest/gtest.h"
#include "setcell.h"
#include "module_md/Langevin.h"

class Langevin_test : public testing::Test
{
protected:
	Verlet *verlet;
    UnitCell_pseudo ucell;

	void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters();
        verlet = new Langevin(INPUT.mdp, ucell);
    }

    void TearDown()
    {
        delete verlet;
    }
};

TEST_F(Langevin_test, setup)
{
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
    verlet->setup();
    verlet->first_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.9925629752939198);
	EXPECT_DOUBLE_EQ(verlet->pos[0].y, 9.9959547635280881);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 0.0020388759545583515);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2060734518595);
	EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1993531037385479);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 9.9965961574993294);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0915381206686439);
	EXPECT_DOUBLE_EQ(verlet->pos[2].y, 9.9965416906633902);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.9957223479323636);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 9.9978511429749233);
	EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.2980537676192396);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 5.0022822443252535);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00017989302497367973);
	EXPECT_DOUBLE_EQ(verlet->vel[0].y, -9.7849590989155838e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, 4.9318050901691171e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, 0.00014690977510726564);
	EXPECT_DOUBLE_EQ(verlet->vel[1].y, -1.5647672277038056e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, -8.233501274764825e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -0.00020468307287393493);
	EXPECT_DOUBLE_EQ(verlet->vel[2].y, -8.3652502505295395e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, -0.00010347145540647848);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, -5.1978365778580354e-05);
	EXPECT_DOUBLE_EQ(verlet->vel[3].y, -4.7077110015593379e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, 5.5204850276109106e-05);
}

TEST_F(Langevin_test, second_half)
{
    verlet->setup();
    verlet->first_half();
    verlet->second_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 0.0020913145423326999);
	EXPECT_DOUBLE_EQ(verlet->pos[0].y, 0.0039920695760518807);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 0.00014322232695976257);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2006470094688808);
	EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1992070870097464);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 9.9982468488089733);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0956952366253319);
	EXPECT_DOUBLE_EQ(verlet->pos[2].y, 0.004957501008778049);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.9966542578184852);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 9.9981783389512984);
	EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.2972719686259842);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 5.0033597142197355);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, 0.00014155502230549501);
	EXPECT_DOUBLE_EQ(verlet->vel[0].y, 0.00011214575426567897);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -1.9603979984090928e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, -0.00020698448409339689);
	EXPECT_DOUBLE_EQ(verlet->vel[1].y, 0.0002442859363436786);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, -8.5606654125729555e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -2.180241323111902e-05);
	EXPECT_DOUBLE_EQ(verlet->vel[2].y, 5.5556410086806277e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, 1.3294522338114925e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, -0.00026189121215375625);
	EXPECT_DOUBLE_EQ(verlet->vel[3].y, -9.4115580823576737e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, 0.00014873968993418975);
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