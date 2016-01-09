#ifndef INVERSE_MATRIX_H
#define INVERSE_MATRIX_H

#include "../src_pw/tools.h"

class Inverse_Matrix_Complex
{
	public:

	Inverse_Matrix_Complex();
	~Inverse_Matrix_Complex();
	
	ComplexMatrix A;
	void using_zpotrf( const ComplexMatrix &Sin);
	void using_zheev(const ComplexMatrix &in, ComplexMatrix &out);
	void init( const int &dim_in);

	private:
	int dim;
	double *e;
	int lwork;
	complex<double> *work2;
	double* rwork;
	int info;
	bool allocate; //mohan add 2012-04-02

	ComplexMatrix EA;
};

class Inverse_Matrix_Real
{
	public:
	
	Inverse_Matrix_Real(){};
	~Inverse_Matrix_Real(){};

	int using_spotri(matrix &A, const int dim);

};


#endif
