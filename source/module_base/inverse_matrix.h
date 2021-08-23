#ifndef INVERSE_MATRIX_H
#define INVERSE_MATRIX_H

#include "../src_pw/tools.h"
namespace ModuleBase
{

class Inverse_Matrix_Complex
{
	public:

	Inverse_Matrix_Complex();
	~Inverse_Matrix_Complex();
	
	ModuleBase::ComplexMatrix A;
	void using_zpotrf( const ModuleBase::ComplexMatrix &Sin);
	void using_zheev(const ModuleBase::ComplexMatrix &in, ModuleBase::ComplexMatrix &out);
	void init( const int &dim_in);

	private:
	int dim;
	double *e;
	int lwork;
	std::complex<double> *work2;
	double* rwork;
	int info;
	bool allocate; //mohan add 2012-04-02

	ModuleBase::ComplexMatrix EA;
};

//not been used yet!
class Inverse_Matrix_Real
{
	public:
	
	Inverse_Matrix_Real(){};
	~Inverse_Matrix_Real(){};

	int using_spotri(matrix &A, const int dim);

};

}
#endif
