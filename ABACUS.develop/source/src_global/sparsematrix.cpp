#include "sparsematrix.h"
#include <string>
#include <iostream>
#include <cstdlib>
using namespace std;

double SparseMatrix::zero = 0.0;

SparseMatrix::SparseMatrix(const int row_num, const int col_num)
{
	this->init(row_num, col_num);
}
SparseMatrix::~SparseMatrix()
{
	delete[] this->RowSet;
}

void SparseMatrix::init(const int row_num, const int col_num)
{
	/*Set Size of the Matrix*/
	assert(row_num >= 0);
	assert(col_num >= 0);

	this->row = row_num;
	this->col = col_num;
	this->size = row_num * col_num;
	this->size_sparse = 0;

	/*Allocate Space for the Matrix*/
	this->RowSet = new Row[row];
	for(int i=0;i<row_num; i++)
	{
		RowSet[i].n_col = 0;
	}
	return;
}

void SparseMatrix::create(const int row_num, const int col_num)
{
	assert(row_num >= 0);
	assert(col_num >= 0);
	delete[] this->RowSet;
	this->init(row_num, col_num);
	return;
}

// Set a member
void SparseMatrix::set(const int row_wanted,const int column_wanted,const double &value)
{
	this->RowSet[row_wanted].member.push_back(value);
	this->RowSet[row_wanted].j_col.push_back(column_wanted);
	this->RowSet[row_wanted].n_col++;
	this->size_sparse++;

//	cout<<"\n set()"<<" row="<<row_wanted<<" col="<<column_wanted<<" value="<<value;
//	cout<<" n_col="<<RowSet[row_wanted].n_col;
	return;
}

void SparseMatrix::reset(
		const int row_wanted,
		const int column_wanted,
		const double &value)
{
	for(int ic=0; ic<RowSet[row_wanted].n_col; ic++)
	{
		if(this->RowSet[ row_wanted].j_col[ic] == column_wanted)
		{
			this->RowSet[row_wanted].member[ic] = value;	
			return;
		}
	}
	this->set(row_wanted, column_wanted, value);
	return;
}



void SparseMatrix::set_add(
		const int row_wanted,
		const int column_wanted,
		const double &value)
{
//	static int count = 0;
//	count++;
//	if(count==20) exit(0);
//	cout<<"\n set_add() "<<row_wanted<<" "<<column_wanted<<endl;

	assert(row_wanted < this->row);
	assert(row_wanted >= 0 );
	assert(column_wanted < this->col);
	assert(column_wanted >= 0);

	for(int ic=0; ic<RowSet[row_wanted].n_col; ic++)
	{
		if(this->RowSet[ row_wanted].j_col[ic] == column_wanted)
		{
			this->RowSet[row_wanted].member[ic] += value;	
			return;
		}
	}

	this->set(row_wanted, column_wanted, value);
	return;
}

SparseMatrix& SparseMatrix::operator=(const SparseMatrix &m)
{
	if(m.row!=this->row || m.col!=this->col)
	{
		cout << "\n row/col number can't match in ComplexMatrix '=' operator\n";
		cout << " this row = " << this->row;
		cout << " this col = " << this->col;
		cout << " in row = " << m.row;
		cout << " in col = " << m.col;
		exit(0);
	}
	else {
		this->create(m.row, m.col);
		this->size_sparse=0;
	}

	for(int ir=0; ir<m.row; ir++)
	{
		for(int ic=0; ic< m.RowSet[ir].n_col; ic++)
		{
//			cout << "\n ir=" << ir << " j_col=" << m.RowSet[ir].j_col[ic]
//				<< " memeber=" << m.RowSet[ir].member[ic];
			this->set( ir, m.RowSet[ir].j_col[ic], m.RowSet[ir].member[ic]);
		}
	}


	return *this;
}


// Operators 
double& SparseMatrix::operator()(const int row_wanted,const int column_wanted)
{
	// searching from all column numbers in this row
	for(int ic=0; ic<RowSet[row_wanted].n_col; ic++)
	{
		if(this->RowSet[row_wanted].j_col[ic] == column_wanted)
		{
			return RowSet[row_wanted].member[ic];
		}
	}
	// if we can't find the member, we are sure the value is zero.
	return zero;
}

const double& SparseMatrix::operator()(const int row_wanted,const int column_wanted)const
{
	for(int ic=0; ic<RowSet[row_wanted].n_col; ic++)
	{
		if(this->RowSet[row_wanted].j_col[ic] == column_wanted)
		{
			return RowSet[row_wanted].member[ic];
		}
	}
	return zero;
}

const matrix SparseMatrix::convert_to_matrix(void)const
{
//	cout << "\n ==> SparseMatrix::convert_to_matrix";
	matrix cm(this->row, this->col);
//	cout << "\n row = " << this->row << " col = " << this->col;
	for(int i=0;i<this->row;i++)
	{
		for(int j=0;j<this->RowSet[i].n_col;j++)
		{	
	//		cout<<"\n i="<<i<<" "<<" j="<<this->RowSet[i].j_col[j]<<" "<<this->RowSet[i].member[j];
			cm(i, this->RowSet[i].j_col[j]) = this->RowSet[i].member[j];
		}
	}
	return cm;
}

void SparseMatrix::multiply_vector(
		const double* psi, // 
		double* h_psi)const
{
	for(int i=0;i<this->row;i++)
	{
		h_psi[i] = 0.0;
	}

	for(int i=0;i<this->row;i++)
	{
		for(int j=0;j<this->RowSet[i].n_col;j++)
		{

//			if( abs( this->RowSet[i].member[j].real() ) > 1.0e-5)
//			cout<<"\n i="<<i
//				<<" j="<<this->RowSet[i].j_col[j]
//				<<" "<<this->RowSet[i].member[j].real()
//				<<" "<<psi[ this->RowSet[i].j_col[j] ].real();

			h_psi[i] += this->RowSet[i].member[j] * psi[ this->RowSet[i].j_col[j] ];
		}
//		cout<<"\n result ==> "<< h_psi[i];
	}
	return;
}

double SparseMatrix::rate(void) const 
{  
	assert(size!=0);
//	cout << "\n size_sparse = " << size_sparse;
	const double occupy_rate = static_cast<double>(size_sparse)/
		static_cast<double>(size);
	return occupy_rate;
}	
