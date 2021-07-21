#include "inverse_matrix.h"

Inverse_Matrix_Complex::Inverse_Matrix_Complex()
{
	allocate=false;	
}
Inverse_Matrix_Complex::~Inverse_Matrix_Complex()
{
	if(allocate)
	{
		delete[] e; //mohan fix bug 2012-04-02
		delete[] work2;
		delete[] rwork;
		allocate=false;
	}	
}

void Inverse_Matrix_Complex::init(const int &dim_in)
{
//	GlobalV::ofs_running << " allocate=" << allocate << endl;
	if(allocate)
	{
		delete[] e; //mohan fix bug 2012-04-02
		delete[] work2;
		delete[] rwork;
		allocate=false;
	}	

	this->dim = dim_in;

	assert(dim>0);
	this->e = new double[dim];
	this->lwork = 2*dim;

	assert(lwork>0);
	this->work2 = new complex<double>[lwork];

	assert(3*dim-2>0);
	this->rwork = new double[3*dim-2];
	this->info = 0;
	this->A.create(dim, dim);
	this->EA.create(dim, dim);

	this->allocate = true;

	return;
}

void Inverse_Matrix_Complex::using_zheev( const ComplexMatrix &Sin, ComplexMatrix &Sout)
{
	timer::tick("Inverse","using_zheev");
	this->A = Sin;

    //LapackConnector::zhegv( 1, 'V', 'U', nwan , A ,  nwan , B , nwan , e, work2 , 80 , rwork , &info );
    LapackConnector::zheev('V', 'U', dim, this->A, dim, e, work2, lwork, rwork, &info);
	
	for(int i=0; i<dim; i++)
	{
		for(int j=0; j<dim; j++)
		{
			EA(i,j)= conj( this->A(j,i) ) / e[i] ;
		}
	}

    Sout = this->A * this->EA;
	timer::tick("Inverse","using_zheev");
    return;
}

void Inverse_Matrix_Complex::using_zpotrf( const ComplexMatrix &Sin)
{
//	timer::tick("Inverse","using_zpotrf");

	for(int i=0; i<dim; i++)
	{
		for(int j=i; j<dim; j++)
		{
			A(i,j) = Sin(i,j);
		}
	}

	LapackConnector::zpotrf('U',dim,A,dim,&info);

	if(info!=0)
	{
		cout << "\n info_zpotrf = " << info;
		QUIT();
	}
	
	LapackConnector::zpotri('U',dim,A,dim,&info);
	
	if(info!=0)
	{
		cout << "\n info_zpotri = " << info;
		QUIT();
	}
//	timer::tick("Inverse","using_zpotrf");
	return;
}

int Inverse_Matrix_Real::using_spotri(matrix &A, const int dim)
{
	int info = 0;
	LapackConnector::spotrf('U',dim,A,dim,&info);
	cout << "\n info_spotrf = " << info;
	LapackConnector::spotri('U',dim,A,dim,&info);
	cout << "\n info_spotri = " << info;
	return info;
}
