#include"mrrr.h"
int sort_w(double *w, int n, int il, int iu, int nsplit, int *blk_begins,
		int *blk_sizes, int *a) {
	int i, nb;
	nb = (nsplit % 2 == 0) ? (nsplit + 1) : nsplit;
	for (i = 0; i < n; ++i) {
		a[i] = i;
	}
	int *a_tmp = (int *) malloc(sizeof(int) * n);
	double *w_tmp = (double *) malloc(sizeof(double) * n);

	int size = blk_sizes[0];

	double *src = w, *des = w_tmp, *tmp;
	int *asrc = a, *ades = a_tmp, *tmpa;
	for (i = 1; i < nsplit; i++) {
		merge(src, src + size, &w[blk_begins[i]],
				&w[blk_begins[i] + blk_sizes[i]], des, asrc, &a[blk_begins[i]],
				ades);
		size += blk_sizes[i];
		tmp = src;
		src = des;
		des = tmp;
		tmpa = asrc;
		asrc = ades;
		ades = tmpa;
	}
	if (nb > nsplit) {
		int ONE = 1;
		pdcopy(&n, w_tmp, &ONE, w, &ONE);
		picopy(&n, a_tmp, &ONE, a, &ONE);
	}

	//xiaohui add 2015-04-07
	free(a_tmp);
	free(w_tmp);

	return 0;
}

void merge(double *x, double *xe, double *y, double *ye, double *z, int *a,
		int *b, int *a1) {
	while (x != xe && y != ye) {
		if (*x <= *y) {
			*z = *x;
			*a1 = *a;
			x++;
			a++;
		} else {
			*z = *y;
			*a1 = *b;
			y++;
			b++;
		}
		z++;
		a1++;
	}
	if (x == xe) {
		while (y != ye) {
			*z = *y;
			*a1 = *b;
			++z;
			++y;
			++a1;
			++b;
		}
	} else {
		while (x != xe) {
			*z = *x;
			*a1 = *a;
			++a1;
			++a;
			++z;
			++x;
		}
	}
}

