//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-06-02
//==========================================================

#include "element_basis_index_unittest.h"

Element_Basis_Index::Range construct_range()
{
	Element_Basis_Index::Range range;
	range.resize(2);
	range[0].resize(3);
	for( size_t L=0; L!=range[0].size(); ++L )
		range[0][L].M=2*L+1;
	range[0][0].N = range[0][1].N = 2;
	range[0][2].N = 1;
	return range;
}

int main()
{
	const Element_Basis_Index::Range &&range 
		= construct_range();
	cout<<range<<endl;
	
	const Element_Basis_Index::Index &&index
		= Element_Basis_Index::construct_index( range );
	cout<<index<<endl;

	return 0;
}

/*
0
        0       1       2
        1       3       2
        2       5       1
1

0:      13
        0
                0
                        0       0
                        1       1
        1
                0
                        0       2
                        1       3
                1
                        0       4
                        1       5
                2
                        0       6
                        1       7
        2
                0
                        0       8
                1
                        0       9
                2
                        0       10
                3
                        0       11
                4
                        0       12
1:      0

14
*/