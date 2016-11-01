//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-06-02
//==========================================================

#ifndef ELEMENT_BASIS_INDEX_H
#define ELEMENT_BASIS_INDEX_H

#include <cstddef>
#include <vector>
using std::vector;

class Element_Basis_Index
{
private:
	
	struct MN
	{
		public:
		size_t M;
		size_t N;
	};
	
	class Index_T: public vector<vector<vector<size_t>>>
	{
		public:
		size_t count_size;
	};	
	
public:
	
	typedef vector<vector<MN>> Range; 
	typedef vector<Index_T> Index;
	
	static Index construct_index( const Range &range );
};

#endif