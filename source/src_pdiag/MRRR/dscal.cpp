#include"mrrr.h"
int pdscal(int *n, double *da, double *dx, int *incx) {
	/* System generated locals */
	int i__1, i__2;

	/* Local variables */
	int i__, m, mp1, nincx;

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */
	/* * */
	/*     scales a vector by a constant. */
	/*     uses unrolled loops for increment equal to one. */
	/*     jack dongarra, linpack, 3/11/78. */
	/*     modified 3/93 to return if incx .le. 0. */
	/*     modified 12/3/93, array(1) declarations changed to array(*) */

	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/* Parameter adjustments */
	--dx;

	/* Function Body */
	if (*n <= 0 || *incx <= 0) {
		return 0;
	}
	if (*incx == 1) {
		goto L20;
	}

	/*        code for increment not equal to 1 */

	nincx = *n * *incx;
	i__1 = nincx;
	i__2 = *incx;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
		dx[i__] = *da * dx[i__];
		/* L10: */
	}
	return 0;

	/*        code for increment equal to 1 */

	/*        clean-up loop */

	L20: m = *n % 5;
	if (m == 0) {
		goto L40;
	}
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
		dx[i__] = *da * dx[i__];
		/* L30: */
	}
	if (*n < 5) {
		return 0;
	}
	L40: mp1 = m + 1;
	i__2 = *n;
	for (i__ = mp1; i__ <= i__2; i__ += 5) {
		dx[i__] = *da * dx[i__];
		dx[i__ + 1] = *da * dx[i__ + 1];
		dx[i__ + 2] = *da * dx[i__ + 2];
		dx[i__ + 3] = *da * dx[i__ + 3];
		dx[i__ + 4] = *da * dx[i__ + 4];
		/* L50: */
	}
	return 0;
} /* dscal_ */
