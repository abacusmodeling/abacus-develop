#include "gtest/gtest.h"
#include "setcell.h"
#include "module_md/MSST.h"
#include "module_esolver/esolver_lj.h"

#define doublethreshold 1e-12

class MSST_test : public testing::Test
{
protected:
    MDrun *mdrun;
    UnitCell_pseudo ucell;

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters();

        ModuleESolver::ESolver *p_esolver = new ModuleESolver::ESolver_LJ();
        p_esolver->Init(INPUT, ucell);

        mdrun = new MSST(INPUT.mdp, ucell);
        mdrun->setup(p_esolver);
    }

    void TearDown()
    {
        delete mdrun;
    }
};

TEST_F(MSST_test, setup)
{
    EXPECT_NEAR(mdrun->vel[0].x, -0.0001314186733659715, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.0985331994796372e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.3947731701005279e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015227275651566311, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014579875939315496, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.5965690649087203e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013311885204189453, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -3.0298400368294885e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.3828659173134662e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011226476889319793, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7843267435287586e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8189299775046767e-05, doublethreshold);
    
    EXPECT_NEAR(mdrun->stress(0,0), 5.9579909955800075e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0,1), -1.4582038138067117e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0,2), 1.4889583894898544e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1,0), -1.4582038138067117e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1,1), 3.4199108345556597e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1,2), -1.2389007575245785e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2,0), 1.4889583894898544e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2,1), -1.2389007575245785e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2,2), 1.5964231736437884e-06, doublethreshold);
}

TEST_F(MSST_test, first_half)
{
    mdrun->first_half();

    EXPECT_NEAR(ucell.lat0, 1.0, doublethreshold);
    EXPECT_NEAR(ucell.lat0_angstrom, 0.52917700000000001, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e11, 10.0, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e12, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e13, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e21, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e22, 10.0, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e23, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e31, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e32, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e33, 9.9959581179144905, doublethreshold);
    EXPECT_NEAR(ucell.omega, 999.59581179144902, doublethreshold);
    
    EXPECT_NEAR(mdrun->pos[0].x, 9.9945728176928519, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.0029442816868202821, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 9.9953814995781549, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 5.2062875654254492, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 5.1939646253791674, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0039673531204680191, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 5.0944925283175158, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 9.9998842371692671, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 4.9957537084439876, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.0046470885642230716, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 5.3032068557647492, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 4.9968136746863676, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013127726219846624, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.121876825065284e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.3947730561390963e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.0001520889345949577, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014598873074918282, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.596568280810794e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013321936931429457, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.8001689685103039e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.3828654775006574e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011240769691879813, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7570131467139791e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8189297471809918e-05, doublethreshold);
}

TEST_F(MSST_test, second_half)
{
    mdrun->first_half();
    mdrun->second_half();

    EXPECT_NEAR(ucell.lat0, 1.0, doublethreshold);
    EXPECT_NEAR(ucell.lat0_angstrom, 0.52917700000000001, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e11, 10.0, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e12, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e13, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e21, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e22, 10.0, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e23, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e31, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e32, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e33, 9.9959581179144905, doublethreshold);
    EXPECT_NEAR(ucell.omega, 999.59581179144902, doublethreshold);
    
    EXPECT_NEAR(mdrun->pos[0].x, 9.9945728176928519, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.0029442816868202821, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 9.9953814995781549, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 5.2062875654254492, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 5.1939646253791674, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0039673531204680191, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 5.0944925283175158, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 9.9998842371692671, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 4.9957537084439876, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.0046470885642230716, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 5.3032068557647492, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 4.9968136746863676, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013113585103096098, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1452204506509308e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.3953371489538059e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015190511267425228, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014617870210521068, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.600449453585996e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013331988658669462, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.5704979001911192e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.3850424881082548e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011255062494439833, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7296995498991997e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8200698165338931e-05, doublethreshold);
}