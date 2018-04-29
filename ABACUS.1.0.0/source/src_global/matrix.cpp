/* matrix.cpp file */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>

#include <cmath>
#include <cstdlib>
#include<limits>

using namespace std;
#include "matrix.h"
#include "lapack_connector.h"

//*********************************************************
// The init() function is the main initialization routine.
// Sets up sizes and allocates memory for matrix class.
// All constructors call init()
// ********************************************************

int matrix::mCount = 0;

void matrixAlloc()
{
	cout << "\n Allocation error for Matrix " << endl;
	abort();
}

void matrix::init(const int nrows,const int ncols)
{
	nr = nrows;
	nc = ncols;
	c = NULL;
//	hermetian = 0;
	
	if(nr*nc == 0) c = NULL;
	else
	{
		auto handler_old = set_new_handler(matrixAlloc);
		c = new double[nrows*ncols]();
		set_new_handler(handler_old);
		
		this->zero_out();
		assert(c!=0);// terminate if memory not allocated
	}

	mCount++;
}

// Free up memory for matrix
void matrix::freemem(void)
{
	if(c)					// Peize Lin add 2016-08-05
	{
		delete [] c;
		c = NULL;
	}
}

// constructor with sizes
matrix::matrix(const int nrows,const int ncols)
{
	this->init(nrows, ncols);
}

//******************
//
// Copy constructor
//
//******************
matrix::matrix(const matrix &m1)
{
	init(m1.nr, m1.nc);
	//hermetian = m1.hermetian;
	// Copy over m1 contents
	for (int i = 0; i < nr*nc; i++) c[i] = m1.c[i];
}

// Peize Lin add 2016-08-05
matrix::matrix( matrix && m1 )
{
	nr=m1.nr; nc=m1.nc;
	c=m1.c;
	m1.init(1,1);
}

//*************
//
// destructor
//
//*************
matrix::~matrix()
{
	this->freemem();
}

//******************************
// reallocate memory for matrix
//******************************
void matrix::create(const int nrow,const int ncol)
{
	delete [] c;
	nr = nrow;
	nc = ncol;
	c = new double[nr * nc];
	zero_out();				// Peize Lin change 2018-03-12
}

// Peize Lin change 2018-03-12
matrix& matrix::operator=( const matrix & m1 )
{
	if( nr*nc != m1.nr*m1.nc )
		this->create( m1.nr, m1.nc );
	for( int i=0; i<nr*nc; ++i ) 
		this->c[i] = m1.c[i];
	return *this;
}

// Peize Lin add 2016-08-05
matrix& matrix::operator=( matrix && m1 )
{
	nr = m1.nr;		nc = m1.nc;
	delete[] c;		c  = m1.c;
	m1.init(1,1);
	return *this;
}

/* Adding matrices, as a friend */
matrix operator+(const matrix &m1, const matrix &m2)
{
	assert(m1.nr == m2.nr);
	assert(m2.nc == m2.nc);

	matrix tm(m1);
	for (int i = 0; i < m1.nr*m1.nc; i++) tm.c[i] += m2.c[i];
	return tm;
}

/* Subtracting matrices, as a friend */
matrix operator-(const matrix &m1, const matrix &m2)
{
	assert(m1.nr == m2.nr);
	assert(m2.nc == m2.nc);

	matrix tm(m1);
	for(int i = 0; i < m1.nr*m1.nc; i++) tm.c[i] -= m2.c[i];
	return tm;
}

//***************************************
//
// Multiplying matrices, as a friend 
//
// *************************************
matrix operator*(const matrix &m1, const matrix &m2)
{
	// fixed bug 2010-01-26
	assert(m1.nc == m2.nr);
	
    // allocate the result and zero it out
    matrix mprod(m1.nr, m2.nc);

    // do the multiply and return
//    for (int i = 0;i < m1.nr;i++)
//    {
//        for (int j = 0;j < m2.nc;j++)
//        {
//            for (int k = 0;k < m1.nc;k++)
//            {
//                //mprod(i, j) += m2(i, k) * m1(k, j);
//                mprod(i, j) += m1(i, k) * m2(k, j);
//            }
//        }
//    }
	
	// Peize Lin accelerate 2017-10-27
	LapackConnector::gemm('N', 'N', m1.nr, m2.nc, m1.nc,
		1, m1.c, m1.nc, m2.c, m2.nc, 
		0, mprod.c, mprod.nc);

	return mprod;
}

/* Scale a matrix */
matrix operator*(const double &s, const matrix &m)
{
	matrix sm(m);
	for (int i = 0; i < m.nr*m.nc; i++) sm.c[i] *= s;
	return sm;
}

/* matrix * double */
matrix operator*(const matrix &m,const double &s)
{
	matrix sm(m);
	for (int i = 0; i < m.nr*m.nc; i++)sm.c[i] *= s;
	return sm;
}

/* Scale a matrix in place */
void matrix::operator*=(const double &s)
{
	for (int i = 0; i < nc*nr; i++) c[i] *= s;
}

/* Accumulate to a matrix in place */
void matrix::operator+=(const matrix & m)
{
	if (nr != m.nr || nc != m.nc)
	{
		cout << " void matrix::operator+=(const matrix &m) has size mismatch\n";
		abort();
	}
	for (int i = 0; i < nc*nr; i++) c[i] += m.c[i];
}


/* decumulate to a matrix in place */
void matrix::operator-=(const matrix & m)
{
	if (nr != m.nr || nc != m.nc)
	{
		cout << "void matrix::operator-=(const matrix &m) has size mismatch\n";
		abort();
	}
	for (int i = 0; i < nc*nr; i++) c[i] -= m.c[i];
}

/* zero out the matrix */
void matrix::zero_out(void)
{
	for(int i = 0; i < nr*nc; i++)
	{
		c[i] = 0.0;
	}
	return;
}


matrix transpose(matrix m)
{
	matrix tm(m.nc, m.nr);

	for (int i = 0;i < m.nr;i++)
	{
		for (int j = 0;j < m.nc;j++)
		{
			tm(j, i) = m(i, j);
		}
	}
	return tm;
}

double matrix::trace_on(void) const
{
    assert(nr == nc);
    int inch = nc + 1;
    int size = nr * nc;
    double tr = 0.0;
    for (int i = 0; i < size; i += inch)
    {
        tr += c[i];
    }
    return tr;
}

void matrix::get_extreme_eigen_values(double &ev_lower, double &ev_upper)const
{
    double *a = new double[nr];
    double *b = new double[nr];
    for (int i = 0; i < nr; ++i)
    {
        double sum = 0.0;
        for(int j = 0; j < nc; ++j)
        {
            sum += fabs(c[i * nc + j]);
        }
        sum -= fabs(c[i * nc + i]);
        a[i] = c[i * nc + i] - sum;
        b[i] = c[i * nc + i] + sum;
    }

    ev_lower = a[0];
    ev_upper = b[0];

    for (int i = 1; i < nr; ++i)
    {
        if (a[i] < ev_lower) ev_lower = a[i];
        if (b[i] > ev_upper) ev_upper = b[i];
    }
    delete[] a;
    delete[] b;
}

// Peize Lin add 2017-05-27
void matrix::reshape( const double nr_new, const double nc_new )
{
	assert( nr*nc == nr_new*nc_new );
	nr=nr_new;
	nc=nc_new;
}

double trace_on(const matrix &A, const matrix &B)
{
    double tr = 0.0;
    for (int i = 0; i < A.nr; ++i)
    {
        for (int k = 0; k < A.nc; ++k)
        {
            tr += A(i,k) * B(k, i);
        }
    }
    return tr;
}

double mdot(const matrix &A, const matrix &B)
{
    assert (A.nr == B.nr);
    assert (A.nc == B.nc);
    int size = A.nr * A.nc;

    double sum = 0.0;
    for (int i = 0; i < size; ++i)
    {
        sum += A.c[i] * B.c[i];
    }
    return sum;
}

// Peize Lin add 2016-09-08
double max( const matrix & m )
{
	double value = std::numeric_limits<double>::min();
	for( int ir=0; ir!=m.nr; ++ir )
		for( int ic=0; ic!=m.nc; ++ic )
			value = std::max( value, m(ir,ic) );
	return value;
}

// Peize Lin add 2016-09-08
double min( const matrix & m )
{
	double value = std::numeric_limits<double>::max();
	for( int ir=0; ir!=m.nr; ++ir )
		for( int ic=0; ic!=m.nc; ++ic )
			value = std::min( value, m(ir,ic) );
	return value;
}

// Peize Lin add 2016-09-08
std::ostream & operator<<( std::ostream & os, const matrix & m )
{
	for( int ir=0; ir!=m.nr; ++ir )
	{
		for( int ic=0; ic!=m.nc; ++ic )
			os<<m(ir,ic)<<"\t";
		os<<std::endl;
	}	
	return os;
}
