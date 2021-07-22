#include "polint.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;

Polint::Polint()
{}

Polint::~Polint()
{}

/*************************************
 * (N-1)th order 
 * Lagrange Polynomial Interpolation
 * N-1 = 3
 * **********************************/
double Polint::Lagrange3
(
 	const double* xa,
	const double* ya,
	const int& n,
	const double& xp
)
{
	double p0, p1, p2, p3;
	double x0, x1, x2, x3;
	double c0, c1, c2, c3;
	double d01, d02, d03, d12, d13, d23;

	int ns = 0;				// Peize Lin add initialization for compiler warning at 2020.01.31
	for(int i = 0; i < n-1; i++)
		if(xp >= xa[i] && xp <= xa[i+1])
			ns = i;
	if(ns > n-4)
		ns = n-4;
	
	x0 = xa[ns];
	x1 = xa[ns+1];
	x2 = xa[ns+2];
	x3 = xa[ns+3];

	c0 = xp-x0;
	c1 = xp-x1;
	c2 = xp-x2;
	c3 = xp-x3;
	
	d01 = x1-x0;
	d02 = x2-x0;
	d03 = x3-x0;
	d12 = x2-x1;
	d13 = x3-x1;
	d23 = x3-x2;

	if(d01 == 0.0 || d02 == 0.0 || d03 == 0.0
			|| d12 == 0.0 || d13 == 0.0 || d23 == 0.0)
	{
		cout << "In Polint::Lagrange3, Two XA's are EQUAL!" << endl;
		exit(0);
	}
	
	p0 = -c1 * c2 * c3 / d01 / d02 / d03 * ya[ns];
	p1 = c0 * c2 * c3 / d01 / d12 / d13 * ya[ns+1];
	p2 = -c0 * c1 * c3 / d02 / d12 / d23 * ya[ns+2];
	p3 = c0 * c1 * c2 / d03 / d13 / d23 * ya[ns+3];
	
	return p0+p1+p2+p3;
}

/************************************************
 * Interpolation for Numerical Orbitals
 * *********************************************/
double Polint::RadialF
(
	const double* rad,
	const double* rad_f,
	const int& msh,
	const int& l,
	const double& R
)
{
  
  int mp_min, mp_max, m;
  double h1, h2, h3, f1, f2, f3, f4;
  double g1, g2, x1, x2, y1, y2, f;
  double c, result;

  mp_min = 0;
  mp_max = msh - 1;

  //assume psir behaves like r**l
  if (R < rad[0])
  {
    if (l == 0) f = rad_f[0];
//    else f = 0.0;
	else 
	{
		c = rad_f[0] / pow(rad[0], l);
		f = pow(R, l) * c;
	}
  }
  else if (rad[mp_max] < R) f = 0.0;
  else
  {
    do
	{
      m = (mp_min + mp_max)/2;
      if (rad[m] < R) mp_min = m;
      else mp_max = m;
    }
    while((mp_max-mp_min)!=1);
    m = mp_max;

    if (m < 2) m = 2;
    else if (msh <= m) m = msh - 2;

    /****************************************************
                   Spline like interpolation
    ****************************************************/

    if (m == 1)
	{
      h2 = rad[m] - rad[m-1];
      h3 = rad[m+1] - rad[m];

      f2 = rad_f[m-1];
      f3 = rad_f[m];
      f4 = rad_f[m+1];

      h1 = -(h2+h3);
      f1 = f4;
    }
    else if (m == (msh-1))
	{
      h1 = rad[m-1] - rad[m-2];
      h2 = rad[m] - rad[m-1];

      f1 = rad_f[m-2];
      f2 = rad_f[m-1];
      f3 = rad_f[m];

      h3 = -(h1+h2);
      f4 = f1;
    }
    else
	{
      h1 = rad[m-1] - rad[m-2];
      h2 = rad[m]   - rad[m-1];
      h3 = rad[m+1] - rad[m];

      f1 = rad_f[m-2];
      f2 = rad_f[m-1];
      f3 = rad_f[m];
      f4 = rad_f[m+1];
    }

    /****************************************************
                Calculate the value at R
    ****************************************************/

    g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
    g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

    x1 = R - rad[m-1];
    x2 = R - rad[m];
    y1 = x1/h2;
    y2 = x2/h2;

    f =  y2*y2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       + y1*y1*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);
  }
  
  result = f;

  return result;
}

