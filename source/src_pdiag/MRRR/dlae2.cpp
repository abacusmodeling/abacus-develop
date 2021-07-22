#include"mrrr.h"
int pdlae2(double *a, double *b, double *c__, double *rt1, double *rt2) {
	/* System generated locals */
	double d__1;

	/* Local variables */
	double ab, df, tb, sm, rt, adf, acmn, acmx;

	/*  -- LAPACK auxiliary routine (version 3.2) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
	/*     November 2006 */

	/*     .. Scalar Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix */
	/*     [  A   B  ] */
	/*     [  B   C  ]. */
	/*  On return, RT1 is the eigenvalue of larger absolute value, and RT2 */
	/*  is the eigenvalue of smaller absolute value. */

	sm = *a + *c__;
	df = *a - *c__;
	adf = fabs(df);
	tb = *b + *b;
	ab = fabs(tb);
	if (fabs(*a) > fabs(*c__)) {
		acmx = *a;
		acmn = *c__;
	} else {
		acmx = *c__;
		acmn = *a;
	}
	if (adf > ab) {
		/* Computing 2nd power */
		d__1 = ab / adf;
		rt = adf * sqrt(d__1 * d__1 + 1.);
	} else if (adf < ab) {
		/* Computing 2nd power */
		d__1 = adf / ab;
		rt = ab * sqrt(d__1 * d__1 + 1.);
	} else {

		/*        Includes case AB=ADF=0 */

		rt = ab * sqrt(2.);
	}
	if (sm < 0.) {
		*rt1 = (sm - rt) * .5;

		/*        Order of execution important. */
		/*        To get fully accurate smaller eigenvalue, */
		/*        next line needs to be executed in higher precision. */

		*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
	} else if (sm > 0.) {
		*rt1 = (sm + rt) * .5;

		/*        Order of execution important. */
		/*        To get fully accurate smaller eigenvalue, */
		/*        next line needs to be executed in higher precision. */

		*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
	} else {

		/*        Includes case RT1 = RT2 = 0 */

		*rt1 = rt * .5;
		*rt2 = rt * -.5;
	}
	return 0;

	/*     End of DLAE2 */

} /* dlae2_ */
