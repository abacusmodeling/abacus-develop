//==========================================================
// Author: ywcui ,mohan
// Last Upstae : 2009-3-8
//==========================================================
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <cassert>
#include <vector>
#include <complex>
#include "matrix.h"

using namespace std;

// define a row in the matrix
struct Row
{
	Row():n_col(0){};
	~Row(){};

	vector< double > member;
	vector<int> j_col;//column number j of a member
	int n_col;//non_zero member number in a row
	// Constructors and destructor
};


class SparseMatrix
{
public:
	// default constructor
	SparseMatrix(const int row_num=0, const int col_num=0);
	~SparseMatrix();

	// matrix size
	int row; //row number
	int col; //column number

	// rows in the matrix
	Row * RowSet;
	
	//allocate new space for a sparse matrix
	void create (const int row_num, const int col_num);

	//add an new element
	void set(
		const int row_wanted,
		const int column_wanted, 
		const double &value);

	void set_add(
		const int row_wanted,
		const int column_wanted,
		const double &value);

	void reset(
		const int row_wanted,
		const int column_wanted,
		const double &value);
	
	// An important function
	void multiply_vector(
		const double* psi, 
		double* h_psi
	)const;
	
	//operators
	double &operator()(const int row,const int column);
	const double &operator()(const int row,const int column)const ;
	SparseMatrix& operator=(const SparseMatrix &m);

	const matrix convert_to_matrix(void)const;

	double rate(void) const;

private:
	void init(const int row_num, const int col_num);
	int size_sparse;
	int size;

	static double zero;

};

#endif
