#include"mrrr.h"

int i_nint(double *x) {
	return (int) (*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
}

