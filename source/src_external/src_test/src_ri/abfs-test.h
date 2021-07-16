#include "../../../src_ri/abfs.h"
#include <iostream>

static void test_adjs( ostream & os )
{
	for( size_t it=0; it!=ucell.ntype; ++it )
		for( size_t ia=0; ia!=ucell.atoms[it].na; ++ia )
		{
			const int iat = ucell.itia2iat(it,ia);
			const map<size_t,vector<Abfs::Vector3_Order<int>>> adjs = Abfs::get_adjs(iat);
			os<<"@@@\t"<<iat<<endl;
			for( const auto & atom2 : adjs )
			{
				const size_t iat2 = atom2.first;
				for( const Abfs::Vector3_Order<int> &box2 : atom2.second )
					os<<iat2<<"\t"<<box2<<endl;
			}
		}
}
