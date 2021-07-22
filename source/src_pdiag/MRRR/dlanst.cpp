#include"mrrr.h"

static int c__1 = 1;

double pdlanst(char *norm, int *n, double *d__, double *e) {
	/* System generated locals */
	int i__1;
	double ret_val, d__1, d__2, d__3, d__4, d__5;


	/* Local variables */
	int i__;
	double sum, scale;
	double anorm;

	/*  -- LAPACK auxiliary routine (version 3.2) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
	/*     November 2006 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLANST  returns the value of the one norm,  or the Frobenius norm, or */
	/*  the  infinity norm,  or the  element of  largest absolute value  of a */
	/*  real symmetric tridiagonal matrix A. */

	/*  Description */
	/*  =========== */

	/*  DLANST returns the value */

	/*     DLANST = ( max(fabs(A(i,j))), NORM = 'M' or 'm' */
	/*              ( */
	/*              ( norm1(A),         NORM = '1', 'O' or 'o' */
	/*              ( */
	/*              ( normI(A),         NORM = 'I' or 'i' */
	/*              ( */
	/*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */

	/*  where  norm1  denotes the  one norm of a matrix (maximum column sum), */
	/*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
	/*  normF  denotes the  Frobenius norm of a matrix (square root of sum of */
	/*  squares).  Note that  max(fabs(A(i,j)))  is not a consistent matrix norm. */

	/*  Arguments */
	/*  ========= */

	/*  NORM    (input) CHARACTER*1 */
	/*          Specifies the value to be returned in DLANST as described */
	/*          above. */

	/*  N       (input) int */
	/*          The order of the matrix A.  N >= 0.  When N = 0, DLANST is */
	/*          set to zero. */

	/*  D       (input) DOUBLE PRECISION array, dimension (N) */
	/*          The diagonal elements of A. */

	/*  E       (input) DOUBLE PRECISION array, dimension (N-1) */
	/*          The (n-1) sub-diagonal or super-diagonal elements of A. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	--e;
	--d__;

	/* Function Body */
	if (*n <= 0) {
		anorm = 0.;
	} else if (plsame(norm, "M")) {

		/*        Find max(fabs(A(i,j))). */

		anorm = (d__1 = d__[*n], fabs(d__1));
		i__1 = *n - 1;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* Computing MAX */
			d__2 = anorm, d__3 = (d__1 = d__[i__], fabs(d__1));
			anorm = max(d__2, d__3);
			/* Computing MAX */
			d__2 = anorm, d__3 = (d__1 = e[i__], fabs(d__1));
			anorm = max(d__2, d__3);
			/* L10: */
		}
	} else if (plsame(norm, "O") || *(unsigned char *) norm == '1'
			|| plsame(norm, "I")) {

		/*        Find norm1(A). */

		if (*n == 1) {
			anorm = fabs(d__[1]);
		} else {
			/* Computing MAX */
			d__3 = fabs(d__[1]) + fabs(e[1]), d__4 = (d__1 = e[*n - 1], fabs(
					d__1)) + (d__2 = d__[*n], fabs(d__2));
			anorm = max(d__3, d__4);
			i__1 = *n - 1;
			for (i__ = 2; i__ <= i__1; ++i__) {
				/* Computing MAX */
				d__4 = anorm, d__5 = (d__1 = d__[i__], fabs(d__1)) + (d__2 =
						e[i__], fabs(d__2)) + (d__3 = e[i__ - 1], fabs(d__3));
				anorm = max(d__4, d__5);
				/* L20: */
			}
		}
	} else if (plsame(norm, "F") || plsame(norm, "E")) {

		/*        Find normF(A). */

		scale = 0.;
		sum = 1.;
		if (*n > 1) {
			i__1 = *n - 1;
			pdlassq(&i__1, &e[1], &c__1, &scale, &sum);
			sum *= 2;
		}
		pdlassq(n, &d__[1], &c__1, &scale, &sum);
		anorm = scale * sqrt(sum);
	}

	ret_val = anorm;
	return ret_val;

	/*     End of DLANST */

} /* dlanst_ */
