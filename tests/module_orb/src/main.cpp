#include "gtest/gtest.h"
int main(int argc, char** argv)
{
// test if the result of Center_2_Orb::Orb11:cal_overlap
// is equal to the result of ORB_gen_table::snap_psipsi
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
