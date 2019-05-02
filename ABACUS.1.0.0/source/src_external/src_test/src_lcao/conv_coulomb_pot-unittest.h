#ifndef CONV_COULOMB_POT_UNITTEST_H
#define CONV_COULOMB_POT_UNITTEST_H

#include "src_lcao/conv_coulomb_pot.h"
#include "src_external/src_test/test_function.h"

static void test_ccp()
{
	Conv_Coulomb_Pot ccp(ORB.Phi[0].PhiLN(0,0));
	ccp.cal_conv_coulomb_pot();
	stringstream ss;
	for( int ir=0; ir<=4*ORB.Phi[0].PhiLN(0,0).getNr(); ++ir )
	{
		ss<<ccp.get_conv_coulomb_pot(ir)<<endl;
	}
	MPI_RANK_OFSTREAM( "Conv_Coulomb_Pot", ss );
}

#endif