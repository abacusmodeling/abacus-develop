#include "mr_interface.h"
#include"mrrr.h"

static int c__1 = 1;
static double c_b18 = .001;

int pdstemr_mpi(MPI_Comm comm1D, char *jobz, char *range, int *n, double * d__,
		double *e, double *vl, double *vu, int *il, int *iu, int *m, double *w,
		double *z__, int *ldz, int *nzc, int *isuppz, bool *tryrac,
		double *work, int *lwork, int *iwork, int *liwork, int *info) {
	/* System generated locals */
	int z_dim1, z_offset, i__1, i__2;
	double d__1, d__2;
	int ix, jx;
	/* Local variables */
	int i__, j;
	double r1 = 0.0, r2 = 0.0;
	int jj;
	double cs;
	int in;
	double sn, wl, wu;
	int iil, iiu;
	double eps, tmp;
	int indd, iend, jblk, wend;
	double rmin, rmax;
	int itmp;
	double tnrm;
	int inde2, itmp2;
	double rtol1, rtol2;
	double scale;
	int indgp;
	int iinfo, iindw, ilast;
	int lwmin;
	bool wantz;
	bool alleig;
	int ibegin;
	bool indeig;
	int iindbl;
	bool valeig;
	int wbegin;
	double safmin;
	double bignum;
	int inderr, iindwk, indgrs, offset;

	double thresh;
	int iinspl, ifirst, indwrk, liwmin, nzcmin;
	double pivmin;
	int nsplit;
	double smlnum;
	bool lquery, zquery;

	int tasks, myid, iiil, iiiu, loc_size, size;

	/* Parameter adjustments */
	--d__;
	--e;
	--w;
	z_dim1 = *ldz;
	z_offset = 1 + z_dim1;
	z__ -= z_offset;
	--isuppz;
	--work;
	--iwork;

	/* Function Body */
	wantz = plsame(jobz, "V");
	alleig = plsame(range, "A");
	valeig = plsame(range, "V");
	indeig = plsame(range, "I");

	/*add  begin*********************************/
	if (indeig && *il == 1 && *iu == *m) {
		*range = 'A';
		alleig = true;
		indeig = false;
	}
	/************************************end add*/

	lquery = *lwork == -1 || *liwork == -1;
	zquery = *nzc == -1;
	/*     DSTEMR needs WORK of size 6*N, IWORK of size 3*N. */
	/*     In addition, DLARRE needs WORK of size 6*N, IWORK of size 5*N. */
	/*     Furthermore, DLARRV needs WORK of size 12*N, IWORK of size 7*N. */
	if (wantz) {
		lwmin = *n * 18;
		liwmin = *n * 10;
	} else {
		/*        need less workspace if only the eigenvalues are wanted */
		lwmin = *n * 12;
		liwmin = *n << 3;
	}
	wl = 0.;
	wu = 0.;
	iil = 1;
	iiu = *n;
	if (valeig) {
		/*        We do not reference VL, VU in the cases RANGE = 'I','A' */
		/*        The interval (WL, WU] contains all the wanted eigenvalues. */
		/*        It is either given by the user or computed in DLARRE. */
		wl = *vl;
		wu = *vu;
	} else if (indeig) {
		/*        We do not reference IL, IU in the cases RANGE = 'V','A' */
		iil = *il;
		iiu = *iu;
	}

	*info = 0;
	if (!(wantz || plsame(jobz, "N"))) {
		*info = -1;
	} else if (!(alleig || valeig || indeig)) {
		*info = -2;
	} else if (*n < 0) {
		*info = -3;
	} else if (valeig && *n > 0 && wu <= wl) {
		*info = -7;
	} else if (indeig && (iil < 1 || iil > *n)) {
		*info = -8;
	} else if (indeig && (iiu < iil || iiu > *n)) {
		*info = -9;
	} else if (*ldz < 1 || wantz && *ldz < *n) {
		*info = -13;
	} else if (*lwork < lwmin && !lquery) {
		*info = -17;
	} else if (*liwork < liwmin && !lquery) {
		*info = -19;
	}

	/*     Get machine constants. */

	safmin = dlamch("Safe minimum");
	eps = dlamch("Precision");
	smlnum = safmin / eps;
	bignum = 1. / smlnum;
	rmin = sqrt(smlnum);
	/* Computing MIN */
	d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
	rmax = min(d__1, d__2);

	if (*info == 0) {
		work[1] = (double) lwmin;
		iwork[1] = liwmin;

		if (wantz && alleig) {
			nzcmin = *n;
		} else if (wantz && valeig) {
			pdlarrc("T", n, vl, vu, &d__[1], &e[1], &safmin, &nzcmin, &itmp,
					&itmp2, info);
		} else if (wantz && indeig) {
			nzcmin = iiu - iil + 1;
		} else {
			/*           WANTZ .EQ. FALSE. */
			nzcmin = 0;
		}
		if (zquery && *info == 0) {
			z__[z_dim1 + 1] = (double) nzcmin;
		} else if (*nzc < nzcmin && !zquery) {
			*info = -14;
		}
	}
	if (*info != 0) {

		i__1 = -(*info);
		xerbla("PDSTEMR", &i__1);

		return 0;
	} else if (lquery || zquery) {
		return 0;
	}

	/*     Handle N = 0, 1, and 2 cases immediately */
	/*add  begin*********************************/
	MPI_Comm_size(comm1D, &tasks);
	MPI_Comm_rank(comm1D, &myid);
	size = *iu - *il + 1;
	loc_size = size / tasks;

	if (myid < size % tasks) {
		++loc_size;
		iiil = loc_size * myid + 1;
	} else {
		iiil = loc_size * myid + size % tasks + 1;
	}
	iiiu = iiil + loc_size - 1;
	//printf("%d proc from %d to %d\n", myid, iiil + iil - 1, iiiu + iil - 1);
	/************************************end add*/
	if (loc_size == 0) {
		return 0;
	}

	*m = 0;

	if (*n == 0) {
		return 0;
	}

	if (*n == 1) {

		if (myid == 0) {/*add*/
			if (alleig || indeig) {
				*m = 1;
				w[1] = d__[1];
			} else {
				if (wl < d__[1] && wu >= d__[1]) {
					*m = 1;
					w[1] = d__[1];
				}
			}
			if (wantz && !zquery) {
				z__[z_dim1 + 1] = 1.;
				isuppz[1] = 1;
				isuppz[2] = 1;
			}
		}
		return 0;
	}

	//printf("n = %d\n", *n);
	if (*n == 2) {
		if (myid == 0) {/*add*/
			if (!wantz) {
				pdlae2(&d__[1], &e[1], &d__[2], &r1, &r2);
			} else if (wantz && !zquery) {
				pdlaev2(&d__[1], &e[1], &d__[2], &r1, &r2, &cs, &sn);
			}
			if (alleig || valeig && r2 > wl && r2 <= wu || indeig && iil == 1) {
				++(*m);
				w[*m] = r2;
				if (wantz && !zquery) {
					z__[*m * z_dim1 + 1] = -sn;
					z__[*m * z_dim1 + 2] = cs;

					if (sn != 0.) {
						if (cs != 0.) {
							isuppz[(*m << 1) - 1] = 1;
							isuppz[(*m << 1) - 1] = 2;
						} else {
							isuppz[(*m << 1) - 1] = 1;
							isuppz[(*m << 1) - 1] = 1;
						}
					} else {
						isuppz[(*m << 1) - 1] = 2;
						isuppz[*m * 2] = 2;
					}
				}
			}
			if (alleig || valeig && r1 > wl && r1 <= wu || indeig && iiu == 2) {
				++(*m);
				w[*m] = r1;
				if (wantz && !zquery) {
					z__[*m * z_dim1 + 1] = cs;
					z__[*m * z_dim1 + 2] = sn;
					if (sn != 0.) {
						if (cs != 0.) {
							isuppz[(*m << 1) - 1] = 1;
							isuppz[(*m << 1) - 1] = 2;
						} else {
							isuppz[(*m << 1) - 1] = 1;
							isuppz[(*m << 1) - 1] = 1;
						}
					} else {
						isuppz[(*m << 1) - 1] = 2;
						isuppz[*m * 2] = 2;
					}
				}
			}
		}
		return 0;
	}
	/*     Continue with general N */
	indgrs = 1;
	inderr = (*n << 1) + 1;
	indgp = *n * 3 + 1;
	indd = (*n << 2) + 1;
	inde2 = *n * 5 + 1;
	indwrk = *n * 6 + 1;

	iinspl = 1;
	iindbl = *n + 1;
	iindw = (*n << 1) + 1;
	iindwk = *n * 3 + 1;

	/*     Scale matrix to allowable range, if necessary. */
	/*     The allowable range is related to the PIVMIN parameter; see the */
	/*     comments in DLARRD.  The preference for scaling small values */
	/*     up is heuristic; we expect users' matrices not to be close to the */
	/*     RMAX threshold. */

	scale = 1.;
	tnrm = pdlanst("M", n, &d__[1], &e[1]);
	if (tnrm > 0. && tnrm < rmin) {
		scale = rmin / tnrm;
	} else if (tnrm > rmax) {
		scale = rmax / tnrm;
	}
	if (scale != 1.) {
		pdscal(n, &scale, &d__[1], &c__1);
		i__1 = *n - 1;
		pdscal(&i__1, &scale, &e[1], &c__1);
		tnrm *= scale;
		if (valeig) {
			/*           If eigenvalues in interval have to be found, */
			/*           scale (WL, WU] accordingly */
			wl *= scale;
			wu *= scale;
		}
	}

	/*     Compute the desired eigenvalues of the tridiagonal after splitting */
	/*     into smaller subblocks if the corresponding off-diagonal elements */
	/*     are small */
	/*     THRESH is the splitting parameter for DLARRE */
	/*     A negative THRESH forces the old splitting criterion based on the */
	/*     size of the off-diagonal. A positive THRESH switches to splitting */
	/*     which preserves relative accuracy. */

	if (*tryrac) {
		/*        Test whether the matrix warrants the more expensive relative approach. */
		pdlarrr(n, &d__[1], &e[1], &iinfo);
	} else {
		/*        The user does not care about relative accurately eigenvalues */
		iinfo = -1;
	}
	/*     Set the splitting criterion */
	if (iinfo == 0) {
		thresh = eps;
	} else {
		thresh = -eps;
		/*        relative accuracy is desired but T does not guarantee it */
		*tryrac = false
		;
	}

	if (*tryrac) {
		/*        Copy original diagonal, needed to guarantee relative accuracy */
		pdcopy(n, &d__[1], &c__1, &work[indd], &c__1);
	}
	/*     Store the squares of the offdiagonal values of T */
	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
		/* Computing 2nd power */
		d__1 = e[j];
		work[inde2 + j - 1] = d__1 * d__1;
		/* L5: */
	}
	/*     Set the tolerance parameters for bisection */
	if (!wantz) {
		/*        DLARRE computes the eigenvalues to full precision. */
		rtol1 = eps * 4.;
		rtol2 = eps * 4.;
	} else {
		/*        DLARRE computes the eigenvalues to less than full precision. */
		/*        DLARRV will refine the eigenvalue approximations, and we can */
		/*        need less accurate initial bisection in DLARRE. */
		/*        Note: these settings do only affect the subset case and DLARRE */
		rtol1 = sqrt(eps);
		/* Computing MAX */
		d__1 = sqrt(eps) * .005, d__2 = eps * 4.;
		rtol2 = max(d__1, d__2);
	}
	time_t ebegin = clock();

	pdlarre(range, n, &wl, &wu, &iil, &iiu, &d__[1], &e[1], &work[inde2],
			&rtol1, &rtol2, &thresh, &nsplit, &iwork[iinspl], m, &w[1],
			&work[inderr], &work[indgp], &iwork[iindbl], &iwork[iindw],
			&work[indgrs], &pivmin, &work[indwrk], &iwork[iindwk], &iinfo);
//	for (i__ = 1; i__ <= *n; i__++) {
//		printf("w[%d]=%lf\t", i__, w[i__] + e[*n]);
//	}
	//printf("\n");
	time_t eend = clock();
//	if (myid == 0)
//		printf("dlarre time is %d and m = %d\n", eend - ebegin, *m);
	if (iinfo != 0) {
		*info = abs(iinfo) + 10;
		return 0;
	}
	iiu = iiu - iil + 1;
	iil = 1;

	/*     Note that if RANGE .NE. 'V', DLARRE computes bounds on the desired */
	/*     part of the spectrum. All desired eigenvalues are contained in */
	/*     (WL,WU] */
	if (wantz) {

		/*parallel this segment*/

		/*        Compute the desired eigenvectors corresponding to the computed */
		/*        eigenvalues */
		/*
		 printf("%d %lf %lf %d %lf %lf\n", *n, wl, wu, *m, rtol1, rtol2);

		 for (ix = 1; ix <= *n; ix++)
		 printf("%lf\t", d__[ix]);
		 printf("\n");
		 for (ix = 1; ix <= *n; ix++)
		 printf("%lf\t", e[ix]);
		 printf("\n");
		 */
		pdlarrv_mpi(&iil, &iiu, n, &wl, &wu, &d__[1], &e[1], &pivmin,
				&iwork[iinspl], m, &iiil, &iiiu, &c_b18, &rtol1, &rtol2, &w[1],
				&work[inderr], &work[indgp], &iwork[iindbl], &iwork[iindw],
				&work[indgrs], &z__[z_offset], ldz, &isuppz[1], &work[indwrk],
				&iwork[iindwk], &iinfo);
		//printf("dlarrv over %d\n", myid);
//		pdlarrv(n, &wl, &wu, &d__[1], &e[1], &pivmin, &iwork[iinspl], m, &c__1,
//				m, &c_b18, &rtol1, &rtol2, &w[1], &work[inderr], &work[indgp],
//				&iwork[iindbl], &iwork[iindw], &work[indgrs], &z__[z_offset],
//				ldz, &isuppz[1], &work[indwrk], &iwork[iindwk], &iinfo);
		if (iinfo != 0) {
			printf("%d err: %d\n",myid,iinfo);
			*info = abs(iinfo) + 20;
			return 0;
		}
	} else {
		/*        DLARRE computes eigenvalues of the (shifted) root representation */
		/*        DLARRV returns the eigenvalues of the unshifted matrix. */
		/*        However, if the eigenvectors are not desired by the user, we need */
		/*        to apply the corresponding shifts from DLARRE to obtain the */
		/*        eigenvalues of the original matrix. */
		i__1 = *m;
		for (j = 1; j <= i__1; ++j) {
			itmp = iwork[iindbl + j - 1];
			w[j] += e[iwork[iinspl + itmp - 1]];
			/* L20: */
		}
	}

	if (*tryrac) {
		/*        Refine computed eigenvalues so that they are relatively accurate */
		/*        with respect to the original matrix T. */
		ibegin = 1;
		wbegin = 1;
		i__1 = iwork[iindbl + *m - 1];
		for (jblk = 1; jblk <= i__1; ++jblk) {
			iend = iwork[iinspl + jblk - 1];
			in = iend - ibegin + 1;
			wend = wbegin - 1;
			/*           check if any eigenvalues have to be refined in this block */
			L36: if (wend < *m) {
				if (iwork[iindbl + wend] == jblk) {
					++wend;
					goto L36;
				}
			}
			if (wend < wbegin) {
				ibegin = iend + 1;
				goto L39;
			}
			offset = iwork[iindw + wbegin - 1] - 1;
			ifirst = iwork[iindw + wbegin - 1];
			ilast = iwork[iindw + wend - 1];
			rtol2 = eps * 4.;
			pdlarrj(&in, &work[indd + ibegin - 1], &work[inde2 + ibegin - 1],
					&ifirst, &ilast, &rtol2, &offset, &w[wbegin],
					&work[inderr + wbegin - 1], &work[indwrk], &iwork[iindwk],
					&pivmin, &tnrm, &iinfo);
			ibegin = iend + 1;
			wbegin = wend + 1;
			L39: ;
		}
	}

	/*     If matrix was scaled, then rescale eigenvalues appropriately. */

	if (scale != 1.) {
		d__1 = 1. / scale;
		pdscal(m, &d__1, &w[1], &c__1);
	}

	/*     If eigenvalues are not in increasing order, then sort them, */
	/*     possibly along with eigenvectors. */

	MPI_Barrier(comm1D);
/*
	if (nsplit > 1) {

		if (!wantz) {
			pdlasrt("I", m, &w[1], &iinfo);
			if (iinfo != 0) {
				*info = 3;
				return 0;
			}
		} else {
			double *z_send = (double *) malloc(
					sizeof(double) * loc_size * z_dim1);
			double *z_recv = (double *) malloc(
					sizeof(double) * loc_size * z_dim1);
			int *blk_begins = (int *) malloc(sizeof(int) * nsplit);
			int *blk_sizes = (int *) malloc(sizeof(int) * nsplit);
			int blkcnt = 0;
			blk_begins[blkcnt] = 0;
			int blksize = 1;
			int *ibls = &iwork[iindbl];
			for (i__ = 1; i__ < *m; i__++) {
				if (ibls[i__] != ibls[i__ - 1]) {
					blk_sizes[blkcnt] = blksize;
					blksize = 1;
					blkcnt++;
					blk_begins[blkcnt] = i__;
				} else {
					blksize++;
				}
			}
			blk_sizes[blkcnt] = blksize;
			int *a = (int *) malloc(sizeof(int) * *m);
			int *sendsize = (int *) malloc(sizeof(int) * tasks);
			int *sendpos = (int *) malloc(sizeof(int) * tasks);
			int *sendpostmp = (int *) malloc(sizeof(int) * tasks);
			int *recvsize = (int *) malloc(sizeof(int) * tasks);
			int *recvpos = (int *) malloc(sizeof(int) * tasks);
			int *recvpostmp = (int *) malloc(sizeof(int) * tasks);
			for (i__ = 0; i__ < tasks; ++i__) {
				sendsize[i__] = 0;
				recvsize[i__] = 0;
			}
			sort_w(&w[1], *m, iiil, iiiu, nsplit, blk_begins, blk_sizes, a);
			iiil--;
			iiiu--;
			int pre, post, sep, PID, npre;
			npre = *m % tasks;
			if (npre == 0) {
				pre = *m / tasks;
				sep = *m;
			} else {
				pre = *m / tasks + 1;
				post = pre - 1;
				sep = pre * npre;
			}

			for (i__ = 0; i__ < *m; ++i__) {
				if (a[i__] <= iiiu && a[i__] >= iiil) {
					PID = (i__ < sep) ? (i__ / pre) : ((i__ - npre) / post);
					sendsize[PID]++;
				}
			}

			for (i__ = iiil; i__ <= iiiu; ++i__) {
				j = a[i__];
				PID = (j < sep) ? (j / pre) : ((j - npre) / post);
				recvsize[PID]++;
			}
			sendpos[0] = 0;
			recvpos[0] = 0;
			for (i__ = 1; i__ < tasks; ++i__) {
				sendpos[i__] = sendpos[i__ - 1] + sendsize[i__ - 1];
			}
			picopy(&tasks, sendpos, &c__1, sendpostmp, &c__1);
			for (i__ = 1; i__ < tasks; ++i__) {
				recvpos[i__] = recvpos[i__ - 1] + recvsize[i__ - 1];
			}
			picopy(&tasks, recvpos, &c__1, recvpostmp, &c__1);

			for (i__ = 0; i__ < *m; ++i__) {
				j = a[i__];
				if (j >= iiil && j <= iiiu) {
					PID = (i__ < sep) ? (i__ / pre) : ((i__ - npre) / post);
					if (PID == myid) {
						pdcopy(n, &z__[(j - iiil + 1) * z_dim1 + 1], &c__1,
								&z_recv[recvpostmp[myid] * z_dim1], &c__1);
						recvpostmp[myid]++;
					} else {
						pdcopy(n, &z__[(j - iiil + 1) * z_dim1 + 1], &c__1,
								&z_send[sendpostmp[PID] * z_dim1], &c__1);

						sendpostmp[PID]++;
					}
				}
			}
			MPI_Status ierr;
			for (i__ = 0; i__ < tasks; ++i__) {
				if (myid == i__) {

				} else {
					if (sendsize[i__] != 0) {
						MPI_Send(&z_send[sendpos[i__] * z_dim1],
								sendsize[i__] * z_dim1, MPI_DOUBLE, i__, myid,
								comm1D);
					}
					if (recvsize[i__] != 0) {
						MPI_Recv(&z_recv[recvpos[i__] * z_dim1],
								recvsize[i__] * z_dim1, MPI_DOUBLE, i__, i__,
								comm1D, &ierr);

					}
				}

			}
			picopy(&tasks, recvpos, &c__1, recvpostmp, &c__1);
			for (i__ = iiil; i__ <= iiiu; ++i__) {
				j = a[i__];
				PID = (j < sep) ? (j / pre) : ((j - npre) / post);
				pdcopy(n, &z_recv[recvpostmp[PID] * z_dim1], &c__1,
						&z__[(i__ - iiil + 1) * z_dim1 + 1], &c__1);

				recvpostmp[PID]++;
			}

			free(z_send);
			free(z_recv);
			free(blk_begins);
			free(blk_sizes);
			free(a);
			free(sendsize);
			free(sendpos);
			free(sendpostmp);
			free(recvsize);
			free(recvpos);
			free(recvpostmp);
		}
	}
*/
	if (nsplit > 1) {

		if (!wantz) {
			pdlasrt("I", m, &w[1], &iinfo);
			if (iinfo != 0) {
				*info = 3;
				return 0;
			}
		} else {
			int *blk_begins = (int *) malloc(sizeof(int) * nsplit);
			int *blk_sizes = (int *) malloc(sizeof(int) * nsplit);
			int blkcnt = 0;
			blk_begins[blkcnt] = 0;
			int blksize = 1;
			int *ibls = &iwork[iindbl];
			for (i__ = 1; i__ < *m; i__++) {
				if (ibls[i__] != ibls[i__ - 1]) {
					blk_sizes[blkcnt] = blksize;
					blksize = 1;
					blkcnt++;
					blk_begins[blkcnt] = i__;
				} else {
					blksize++;
				}
			}
			blk_sizes[blkcnt] = blksize;
			/*add*/blkcnt++;

			if (blkcnt > 1) {
				double *z_send = (double *) malloc(
						sizeof(double) * loc_size * z_dim1);
				double *z_recv = (double *) malloc(
						sizeof(double) * loc_size * z_dim1);
				int *a = (int *) malloc(sizeof(int) * *m);
				int *sendsize = (int *) malloc(sizeof(int) * tasks);
				int *sendpos = (int *) malloc(sizeof(int) * tasks);
				int *sendpostmp = (int *) malloc(sizeof(int) * tasks);
				int *recvsize = (int *) malloc(sizeof(int) * tasks);
				int *recvpos = (int *) malloc(sizeof(int) * tasks);
				int *recvpostmp = (int *) malloc(sizeof(int) * tasks);
				for (i__ = 0; i__ < tasks; ++i__) {
					sendsize[i__] = 0;
					recvsize[i__] = 0;
				}
				sort_w(&w[1], *m, iiil, iiiu, blkcnt, blk_begins, blk_sizes, a);
				iiil--;
				iiiu--;
				int pre, post, sep, PID, npre;
				npre = *m % tasks;
				if (npre == 0) {
					pre = *m / tasks;
					sep = *m;
				} else {
					pre = *m / tasks + 1;
					post = pre - 1;
					sep = pre * npre;
				}

				for (i__ = 0; i__ < *m; ++i__) {
					if (a[i__] <= iiiu && a[i__] >= iiil) {
						PID = (i__ < sep) ? (i__ / pre) : ((i__ - npre) / post);
						sendsize[PID]++;
					}
				}

				for (i__ = iiil; i__ <= iiiu; ++i__) {
					j = a[i__];
					PID = (j < sep) ? (j / pre) : ((j - npre) / post);
					recvsize[PID]++;
				}
				sendpos[0] = 0;
				recvpos[0] = 0;
				for (i__ = 1; i__ < tasks; ++i__) {
					sendpos[i__] = sendpos[i__ - 1] + sendsize[i__ - 1];
				}
				picopy(&tasks, sendpos, &c__1, sendpostmp, &c__1);
				for (i__ = 1; i__ < tasks; ++i__) {
					recvpos[i__] = recvpos[i__ - 1] + recvsize[i__ - 1];
				}
				picopy(&tasks, recvpos, &c__1, recvpostmp, &c__1);

				for (i__ = 0; i__ < *m; ++i__) {
					j = a[i__];
					if (j >= iiil && j <= iiiu) {
						PID = (i__ < sep) ? (i__ / pre) : ((i__ - npre) / post);
						if (PID == myid) {
							pdcopy(n, &z__[(j - iiil + 1) * z_dim1 + 1], &c__1,
									&z_recv[recvpostmp[myid] * z_dim1], &c__1);
							recvpostmp[myid]++;
						} else {
							pdcopy(n, &z__[(j - iiil + 1) * z_dim1 + 1], &c__1,
									&z_send[sendpostmp[PID] * z_dim1], &c__1);

							sendpostmp[PID]++;
						}
					}
				}
				MPI_Status ierr;
				for (i__ = 0; i__ < tasks; ++i__) {
					if (myid == i__) {

					} else {
						if (sendsize[i__] != 0) {
							MPI_Send(&z_send[sendpos[i__] * z_dim1],
									sendsize[i__] * z_dim1, MPI_DOUBLE, i__,
									myid, comm1D);
						}
						if (recvsize[i__] != 0) {
							MPI_Recv(&z_recv[recvpos[i__] * z_dim1],
									recvsize[i__] * z_dim1, MPI_DOUBLE, i__,
									i__, comm1D, &ierr);

						}
					}

				}
				picopy(&tasks, recvpos, &c__1, recvpostmp, &c__1);
				for (i__ = iiil; i__ <= iiiu; ++i__) {
					j = a[i__];
					PID = (j < sep) ? (j / pre) : ((j - npre) / post);
					pdcopy(n, &z_recv[recvpostmp[PID] * z_dim1], &c__1,
							&z__[(i__ - iiil + 1) * z_dim1 + 1], &c__1);

					recvpostmp[PID]++;
				}
				free(z_send);
				free(z_recv);
				free(a);
				free(sendsize);
				free(sendpos);
				free(sendpostmp);
				free(recvsize);
				free(recvpos);
				free(recvpostmp);
			}
			free(blk_begins);
			free(blk_sizes);
		}
	}
	return 0;

	/*     End of DSTEMR */

} /* dstemr_ */
