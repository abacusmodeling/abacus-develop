#include "gtest/gtest.h"
#include "setcell.h"
#include "module_md/verlet.h"
#include "module_esolver/esolver_lj.h"

#define doublethreshold 1e-12

class Verlet_test : public testing::Test
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

        mdrun = new Verlet(INPUT.mdp, ucell);
        mdrun->setup(p_esolver);
    }

    void TearDown()
    {
        delete mdrun;
    }
};

TEST_F(Verlet_test, setup)
{   
    EXPECT_NEAR(mdrun->temperature_, 299.99999999999665, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0,0), 6.0100555286436806e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0,1), -1.4746713013791574e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0,2), 1.5039983732220751e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1,0), -1.4746713013791574e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1,1), 3.4437172989317909e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1,2), -1.251414906590483e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2,0), 1.5039983732220751e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2,1), -1.251414906590483e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2,2), 1.6060561926126463e-06, doublethreshold);
}

TEST_F(Verlet_test, first_half)
{
    mdrun->first_half();
    
    EXPECT_NEAR(mdrun->pos[0].x, 9.9945454470992772, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.0029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 9.9994204767196599, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 5.2063192793031234, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 5.1939342598421803, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0039873402383468889, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 5.0944648036273872, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 9.9998836061438716, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 4.997763438399228, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.0046704699702541427, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 5.303223068197739, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 4.9988287446427613, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013193932519649473, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1576379239356465e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.40179977966e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015285605661129458, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014672323796402785, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6449148069800003e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013388999749840003, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.8154327428808153e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4099838013700003e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.0001129732660846002, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7962291467652202e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.83313122596e-05, doublethreshold);
}

TEST_F(Verlet_test, NVE)
{
    mdrun->first_half();
    mdrun->second_half();
    
    EXPECT_NEAR(mdrun->pos[0].x, 9.9945454470992772, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.0029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 9.9994204767196599, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 5.2063192793031234, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 5.1939342598421803, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0039873402383468889, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 5.0944648036273872, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 9.9998836061438716, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 4.997763438399228, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.0046704699702541427, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 5.303223068197739, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 4.9988287446427613, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013179791402898947, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1809815495212933e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.40179977966e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015267223469058917, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014691320932005571, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6449148069800003e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013399051477080008, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.5857616745616307e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4099838013700003e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.0001131161941102004, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7689155499504408e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.83313122596e-05, doublethreshold);
}

TEST_F(Verlet_test, Anderson)
{
    mdrun->first_half();
    mdrun->mdp.md_thermostat = "Anderson";
    mdrun->second_half();
    
    EXPECT_NEAR(mdrun->pos[0].x, 9.9945454470992772, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.0029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 9.9994204767196599, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 5.2063192793031234, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 5.1939342598421803, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0039873402383468889, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 5.0944648036273872, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 9.9998836061438716, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 4.997763438399228, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.0046704699702541427, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 5.303223068197739, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 4.9988287446427613, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013179791402898947, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1809815495212933e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.40179977966e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 6.9452562329904563e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, 7.321611395307015e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, -8.133446733603267e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, 0.00013239881096711222, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, 0.00030862680563211305, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -0.00012925479702246553, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.0001131161941102004, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7689155499504408e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.83313122596e-05, doublethreshold);
}

TEST_F(Verlet_test, Berendsen)
{
    mdrun->first_half();
    mdrun->mdp.md_thermostat = "Berendsen";
    mdrun->second_half();
    
    EXPECT_NEAR(mdrun->pos[0].x, 9.9945454470992772, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.0029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 9.9994204767196599, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 5.2063192793031234, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 5.1939342598421803, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0039873402383468889, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 5.0944648036273872, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 9.9998836061438716, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 4.997763438399228, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.0046704699702541427, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 5.303223068197739, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 4.9988287446427613, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013179175250738632, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1806458403162173e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.4017342458487154e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015266509729938552, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.0001469063411619389, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6444639094723906e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013398425074562592, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.5856407908091386e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4097308858947404e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011311090595462667, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7685523549685863e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8329987777389342e-05, doublethreshold);
}

TEST_F(Verlet_test, Rescaling)
{
    mdrun->first_half();
    mdrun->mdp.md_thermostat = "Rescaling";
    mdrun->second_half();
    
    EXPECT_NEAR(mdrun->pos[0].x, 9.9945454470992772, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.0029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 9.9994204767196599, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 5.2063192793031234, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 5.1939342598421803, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0039873402383468889, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 5.0944648036273872, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 9.9998836061438716, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 4.997763438399228, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.0046704699702541427, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 5.303223068197739, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 4.9988287446427613, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013178559069770653, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1803101154153484e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.4016687089734539e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015265795957447931, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.0001468994726827073, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6440129908834563e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013397798642758268, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.5855199014048311e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4094779585946356e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.0001131056175518098, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7681891430058639e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8328663233253657e-05, doublethreshold);
}

TEST_F(Verlet_test, Rescale_v)
{
    mdrun->first_half();
    mdrun->mdp.md_thermostat = "Rescale_v";
    mdrun->second_half();
    
    EXPECT_NEAR(mdrun->pos[0].x, 9.9945454470992772, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.0029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 9.9994204767196599, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 5.2063192793031234, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 5.1939342598421803, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0039873402383468889, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 5.0944648036273872, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 9.9998836061438716, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 4.997763438399228, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.0046704699702541427, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 5.303223068197739, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 4.9988287446427613, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013178559069770653, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1803101154153484e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.4016687089734539e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015265795957447931, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.0001468994726827073, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6440129908834563e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013397798642758268, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.5855199014048311e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4094779585946356e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.0001131056175518098, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7681891430058639e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8328663233253657e-05, doublethreshold);
}