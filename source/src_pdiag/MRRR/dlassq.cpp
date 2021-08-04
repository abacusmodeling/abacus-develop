#include"mrrr.h"
/* Subroutine */int pdlassq(int *n, double *x, int *incx, double *scale,
		double *sumsq) {
	/* System generated locals */
	int i__1, i__2;
	double d__1;

	/* Local variables */
	int ix;
	double absxi;

	/*  -- LAPACK auxiliary routine (version 3.2) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
	/*     November 2006 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLASSQ  returns the values  scl  and  smsq  such that */

	/*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq, */

	/*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is */
	/*  assumed to be non-negative and  scl  returns the value */

	/*     scl = max( scale, fabs( x( i ) ) ). */

	/*  scale and sumsq must be supplied in SCALE and SUMSQ and */
	/*  scl and smsq are overwritten on SCALE and SUMSQ respectively. */

	/*  The routine makes only one pass through the std::vector x. */

	/*  Arguments */
	/*  ========= */

	/*  N       (input) int */
	/*          The number of elements to be used from the std::vector X. */

	/*  X       (input) DOUBLE PRECISION array, dimension (N) */
	/*          The std::vector for which a scaled sum of squares is computed. */
	/*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n. */

	/*  INCX    (input) int */
	/*          The increment between successive values of the std::vector X. */
	/*          INCX > 0. */

	/*  SCALE   (input/output) DOUBLE PRECISION */
	/*          On entry, the value  scale  in the equation above. */
	/*          On exit, SCALE is overwritten with  scl , the scaling factor */
	/*          for the sum of squares. */

	/*  SUMSQ   (input/output) DOUBLE PRECISION */
	/*          On entry, the value  sumsq  in the equation above. */
	/*          On exit, SUMSQ is overwritten with  smsq , the basic sum of */
	/*          squares from which  scl  has been factored out. */

	/* ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	--x;

	/* Function Body */
	if (*n > 0) {
		i__1 = (*n - 1) * *incx + 1;
		i__2 = *incx;
		for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
			if (x[ix] != 0.) {
				absxi = (d__1 = x[ix], fabs(d__1));
				if (*scale < absxi) {
					/* Computing 2nd power */
					d__1 = *scale / absxi;
					*sumsq = *sumsq * (d__1 * d__1) + 1;
					*scale = absxi;
				} else {
					/* Computing 2nd power */
					d__1 = absxi / *scale;
					*sumsq += d__1 * d__1;
				}
			}
			/* L10: */
		}
	}
	return 0;

	/*     End of DLASSQ */

} /* dlassq_ */
