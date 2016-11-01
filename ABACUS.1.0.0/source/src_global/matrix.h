#ifndef MATRIX_H
#define MATRIX_H

#ifdef _MCD_CHECK
#include "../src_parallel/mcd.h"
#endif

class matrix
{
	/* data */
public:
	double *c;    /* Holds the data */
	int nr;
	int nc;   /* Number of rows and columns */
	int hermetian;
	static int mCount;

	/* Constructors and destructor */
	matrix(const int nrows = 1,const int ncols = 1);
	matrix(const matrix &m1); /* copy constructor */
	matrix( matrix && m1 );			// Peize Lin add 2016-08-05
	~matrix();

	void create(const int nrow,const int ncol);
	void operator=(const matrix &m1); /* Nonstandard: returns void */
	matrix& operator=( matrix && m1 );	// Peize Lin add 2016-08-05

	inline double &operator()(const int &ir,const int &ic)
	{ return c[nc*ir + ic]; }

	inline const double &operator()(const int ir,const int ic)const
	{ return c[nc*ir + ic]; }

//	inline double &operator()(int i,int j) const
//	    { return c[nc*i+j]; }

	friend matrix operator+(const matrix &m1, const matrix &m2);
	friend matrix operator-(const matrix &m1, const matrix &m2);
	friend matrix operator*(const matrix &m1, const matrix &m2);
	friend matrix operator*(const double &s, const matrix &m);
	friend matrix operator*(const matrix &m, const double &s);
	void operator*=(const double &s);
	void operator+=(const matrix &m);
	void operator-=(const matrix &m);

	/* member function: */
	matrix Inverse(void);

	double det(void);

	// mohan add 2011-01-13
	double trace_on(void) const;

	/* Does the memory allocations of the constructor */
	void init(int nrows, int ncols);

	/* Free up memory */
	void freemem(void);

	/* zero out all the entries */
	void zero_out(void);

	// mohan add 2011-01-13
	void get_extreme_eigen_values(double &ev_lower, double &ev_upper) const;
};

matrix transpose(matrix m);

// mohan add 2011-01-13
double trace_on(const matrix &A, const matrix &B);

// mohan add 2011-01-13
double mdot(const matrix &A, const matrix &B);

// Peize Lin add 2016-09-08
double max( const matrix & m );
double min( const matrix & m );

#endif // MATRIX_H
