#ifndef EXX_OPT_ORB_H
#define EXX_OPT_ORB_H

#include "../module_base/matrix.h"
#include "../module_base/element_basis_index.h"
#include <vector>
#include <map>
#include <set>
using namespace std;

class Exx_Opt_Orb
{
public:
	void generate_matrix() const;
private:
	std::vector<std::vector<matrix>> cal_I( 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,matrix>>>> &ms, 
		const size_t TA, const size_t IA, const size_t TB, const size_t IB ) const;
	matrix cal_proj( 
		const matrix & m_big, 
		const std::vector<matrix> & m_left, 
		const std::vector<std::vector<matrix>> & m_middle, 
		const std::vector<matrix> & m_right ) const;
	void print_matrix(
		const std::string &file_name,
		const std::vector<matrix> &matrix_Q, 
		const std::vector<std::vector<matrix>> &matrix_S,
		const matrix &matrix_V,
		const size_t TA, const size_t IA, const size_t TB, const size_t IB,
		const Element_Basis_Index::Range &range_jles, 
		const Element_Basis_Index::IndexLNM &index_jles) const;
	std::map<size_t,std::map<size_t,std::set<double>>> get_radial_R() const;
		
	int kmesh_times = 4;
};

#endif