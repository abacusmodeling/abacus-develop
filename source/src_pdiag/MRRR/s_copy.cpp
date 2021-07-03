#include"mrrr.h"

/* assign strings:  a = b */

void s_copy(char *a, char *b, int la, int lb)

{
	char *aend, *bend;

	aend = a + la;

	if (la <= lb)
#ifndef NO_OVERWRITE
		if (a <= b || a >= b + la)
#endif
			while (a < aend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else
			for (b += la; a < aend;)
				*--aend = *--b;
#endif

	else {
		bend = b + lb;
#ifndef NO_OVERWRITE
		if (a <= b || a >= bend)
#endif
			while (b < bend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else {
			a += lb;
			while (b < bend)
				*--a = *--bend;
			a += lb;
		}
#endif
		while (a < aend)
			*a++ = ' ';
	}
}
