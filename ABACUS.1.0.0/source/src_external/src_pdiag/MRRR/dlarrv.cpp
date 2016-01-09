#include"mrrr.h"
static double c_b5 = 0.;
static int c__1 = 1;
static int c__2 = 2;

int pdlarrv_mpi(int *il, int *iu, int *n, double *vl, double *vu, double *d__,
		double *l, double *pivmin, int *isplit, int *m, int *dol, int *dou,
		double *minrgp, double *rtol1, double *rtol2, double *w, double *werr,
		double *wgap, int *iblock, int *indexw, double *gers, double *z__,
		int *ldz, int *isuppz, double *work, int *iwork, int *info) {
	/* System generated locals */

	//printf("On entry dlarrv_mpi, il=%d, iu=%d, dol=%d, dou=%d\n", *il, *iu,*dol, *dou);
	int z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
	double d__1, d__2;
	bool L__1;

	/* Local variables */
	int minwsize, i__, j, k, p, q, miniwsize, ii;
	double gl;
	int im, in;
	double gu, gap, eps, tau, tol, tmp;
	int zto;
	double ztz;
	int iend, jblk;
	double lgap;
	int done;
	double rgap, left;
	int wend, iter;
	double bstw;
	int itmp1;

	int indld;
	double fudge;
	int idone;
	double sigma;
	int iinfo, iindr;
	double resid;
	bool eskip;
	double right;
	int nclus, zfrom;
	double rqtol;
	int iindc1, iindc2;

	bool stp2ii;
	double lambda;
	int ibegin, indeig;
	bool needbs;
	int indlld;
	double sgndef, mingma;

	int oldien, oldncl, wbegin;
	double spdiam;
	int negcnt;

	int oldcls;
	double savgap;
	int ndepth;
	double ssigma;

	bool usedbs;
	int iindwk, offset;
	double gaptol;
	int newcls, oldfst, indwrk, windex, oldlst;
	bool usedrq;
	int newfst, newftt, parity, windmn, windpl, isupmn, newlst, zusedl;
	double bstres;
	int newsiz, zusedu, zusedw;
	double nrminv, rqcorr;
	bool tryrqc;
	int isupmx;

	int ne;
	double *dfst, *lfst, *dlst, *llst;
	dfst = (double *) malloc(sizeof(double) * *n);
	lfst = (double *) malloc(sizeof(double) * *n);
	dlst = (double *) malloc(sizeof(double) * *n);
	llst = (double *) malloc(sizeof(double) * *n);

	--d__;
	--l;
	--isplit;
	--w;
	--werr;
	--wgap;
	--iblock;
	--indexw;
	--gers;
	z_dim1 = *ldz;
	z_offset = 1 + z_dim1 * *dol;
	z__ -= z_offset;
	--isuppz;
	--work;
	--iwork;

	/* Function Body */
	indld = *n + 1;
	indlld = (*n << 1) + 1;
	indwrk = *n * 3 + 1;
	minwsize = *n * 12;
	i__1 = minwsize;
	for (i__ = 1; i__ <= i__1; ++i__) {
		work[i__] = 0.;
		/* L5: */
	}
	/*     IWORK(IINDR+1:IINDR+N) hold the twist indices R for the */
	/*     factorization used to compute the FP vector */
	iindr = 0;
	/*     IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current */
	/*     layer and the one above. */
	iindc1 = *n;
	iindc2 = *n << 1;
	iindwk = *n * 3 + 1;
	miniwsize = *n * 7;
	i__1 = miniwsize;
	for (i__ = 1; i__ <= i__1; ++i__) {
		iwork[i__] = 0;
		/* L10: */
	}

	zusedl = *dol;
	zusedu = *dou;
	/*     The width of the part of Z that is used */
	zusedw = zusedu - zusedl + 1;
	pdlaset("Full", n, &zusedw, &c_b5, &c_b5, &z__[zusedl * z_dim1 + 1], ldz);
	eps = dlamch("Precision");
	rqtol = eps * 2.;

	/*     Set expert flags for standard code. */
	tryrqc = true;
	if (*il == 1 && *iu == *m) {
	} else {
		/*        Only selected eigenpairs are computed. Since the other evalues */
		/*        are not refined by RQ iteration, bisection has to compute to full */
		/*        accuracy. */
		*rtol1 = eps * 4.;
		*rtol2 = eps * 4.;
	}
	/*     The entries WBEGIN:WEND in W, WERR, WGAP correspond to the */
	/*     desired eigenvalues. The support of the nonzero eigenvector */
	/*     entries is contained in the interval IBEGIN:IEND. */
	/*     Remark that if k eigenpairs are desired, then the eigenvectors */
	/*     are stored in k contiguous columns of Z. */
	/*     DONE is the number of eigenvectors already computed */
	done = 0;
//printf("%d %d %d %d %d %d %d %d\n",iblock[1],iblock[2],iblock[3],iblock[4],iblock[5],iblock[6],iblock[7],iblock[8]);
	int *blk_begin = (int *)malloc(sizeof(int) * *m);
	for (i__ = 0; i__ < *m; ++i__) {
		blk_begin[i__] = 0;
	}
	int *blk_end = (int *)malloc(sizeof(int) * *m);
	int block_num = iblock[1];
	blk_begin--;
	blk_end--;
	blk_begin[1] = 1;
	blk_end[1] = 1;
	for (i__ = 2; i__ <= *m; ++i__) {
		if (iblock[i__] != iblock[i__ - 1]) {
			blk_end[block_num] = i__ - 1;
			block_num = iblock[i__];
			blk_begin[block_num] = i__;
		}
	}
	blk_end[block_num] = *m;

	if (iblock[*dol] == 1) {
		ibegin = 1;
		//wbegin = 1;

	} else {
		ibegin = isplit[iblock[*dol] - 1] + 1;
		//wbegin = isplit[iblock[*dol] - 1] + 1;
	}
	//wbegin = blk_begin[iblock[*dol]];
	//ibegin=1;
	//wbegin = 1;
	//i__1 = iblock[*m];
	i__1 = iblock[*dou];
	for (jblk = iblock[*dol]; jblk <= i__1; ++jblk) {
		if (blk_begin[jblk] == 0) {
			continue;
		}
		//printf("block[%d]\n", jblk);
		iend = isplit[jblk];
		sigma = l[iend];

		if (iblock[*dol] == 1) {
			ibegin = 1;
		} else {
			ibegin = isplit[jblk - 1] + 1;
		}
		wbegin = blk_begin[jblk];
		wend = blk_end[jblk];

		//2014-07-02
		//free(blk_begin);
		//free(blk_end);

		/*        Find local spectral diameter of the block */
		gl = gers[(ibegin << 1) - 1];
		gu = gers[ibegin * 2];
		i__2 = iend;
		for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
			/* Computing MIN */
			d__1 = gers[(i__ << 1) - 1];
			gl = min(d__1, gl);
			/* Computing MAX */
			d__1 = gers[i__ * 2];
			gu = max(d__1, gu);
			/* L20: */
		}
		spdiam = gu - gl;
		/*        OLDIEN is the last index of the previous block */
		oldien = ibegin - 1;
		/*        Calculate the size of the current block */
		in = iend - ibegin + 1;
		/*        The number of eigenvalues in the current block */
		im = wend - wbegin + 1;

		ne = (*dou < wend ? *dou : wend) - (*dol > wbegin ? *dol : wbegin) + 1;

		/*        This is for a 1x1 block */
		if (ibegin == iend) {
			++done;
			z__[ibegin + wbegin * z_dim1] = 1.;
			isuppz[(wbegin << 1) - 1] = ibegin;
			isuppz[wbegin * 2] = ibegin;
			w[wbegin] += sigma;
			work[wbegin] = w[wbegin];
			ibegin = iend + 1;
			++wbegin;
			goto L170;
		}
		/*        The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND) */
		/*        Note that these can be approximations, in this case, the corresp. */
		/*        entries of WERR give the size of the uncertainty interval. */
		/*        The eigenvalue approximations will be refined when necessary as */
		/*        high relative accuracy is required for the computation of the */
		/*        corresponding eigenvectors. */
		pdcopy(&im, &w[wbegin], &c__1, &work[wbegin], &c__1);
		/*        We store in W the eigenvalue approximations w.r.t. the original */
		/*        matrix T. */
		i__2 = im;
		for (i__ = 1; i__ <= i__2; ++i__) {
			w[wbegin + i__ - 1] += sigma;
			/* L30: */
		}
		/*        NDEPTH is the current depth of the representation tree */
		ndepth = 0;
		/*        PARITY is either 1 or 0 */
		parity = 1;
		/*        NCLUS is the number of clusters for the next level of the */
		/*        representation tree, we start with NCLUS = 1 for the root */
		nclus = 1;
		iwork[iindc1 + 1] = wbegin;
		iwork[iindc1 + 2] = wend;
		/*        IDONE is the number of eigenvectors already computed in the current */
		/*        block */
		idone = 0;
		/*        loop while( IDONE.LT.IM ) */
		/*        generate the representation tree for the current block and */
		/*        compute the eigenvectors */
		L40: if (idone < ne) {
			/*           This is a crude protection against infinitely deep trees */
			if (ndepth > *m) {

//	if(*dol==5)	printf("ndepth=%d\n",ndepth);
				*info = -2;
				//xiaohui add 2015-04-07
				free(dfst);
				free(lfst);
				free(dlst);
				free(llst);

				blk_begin++;
				blk_end++;
				free(blk_begin);
				free(blk_end);

				return 0;
			}
			/*           breadth first processing of the current level of the representation */
			/*           tree: OLDNCL = number of clusters on current level */
			oldncl = nclus;
			/*           reset NCLUS to count the number of child clusters */
			nclus = 0;

			parity = 1 - parity;
			if (parity == 0) {
				oldcls = iindc1;
				newcls = iindc2;
			} else {
				oldcls = iindc2;
				newcls = iindc1;
			}
			/*           Process the clusters on the current level */
			i__2 = oldncl;
//if(*dol==5)printf("%d  %d\n",parity,i__2);
			for (i__ = 1; i__ <= i__2; ++i__) {

				j = oldcls + (i__ << 1);
				/*              OLDFST, OLDLST = first, last index of current cluster. */
				/*                               cluster indices start with 1 and are relative */
				/*                               to WBEGIN when accessing W, WGAP, WERR, Z */
				oldfst = iwork[j - 1];
				oldlst = iwork[j];
				if (ndepth > 0) {

					if (wbegin + oldfst - 1 < *dol) {

						pdcopy(&in, dfst, &c__1, &d__[ibegin], &c__1);
						i__3 = in - 1;
						pdcopy(&i__3, lfst, &c__1, &l[ibegin], &c__1);
						sigma = lfst[iend - 1];
						for (j = ibegin - 1; j < iend; j++) {
							dfst[j] = 0.0;
							lfst[j] = 0.0;
						}
					} else if (wbegin + oldlst - 1 > *dou) {
						pdcopy(&in, dlst, &c__1, &d__[ibegin], &c__1);
						i__3 = in - 1;
						pdcopy(&i__3, llst, &c__1, &l[ibegin], &c__1);
						sigma = llst[iend - 1];
						for (j = ibegin - 1; j < iend; j++) {
							dlst[j] = 0.0;
							llst[j] = 0.0;
						}
					} else {
						j = wbegin + oldfst - 1;
						pdcopy(&in, &z__[ibegin + j * z_dim1], &c__1,
								&d__[ibegin], &c__1);
						i__3 = in - 1;
						pdcopy(&i__3, &z__[ibegin + (j + 1) * z_dim1], &c__1,
								&l[ibegin], &c__1);
						sigma = z__[iend + (j + 1) * z_dim1];
						/*                 Set the corresponding entries in Z to zero */
						pdlaset("Full", &in, &c__2, &c_b5, &c_b5,
								&z__[ibegin + j * z_dim1], ldz);
					}

				}
				/*              Compute DL and DLL of current RRR */
				i__3 = iend - 1;
				for (j = ibegin; j <= i__3; ++j) {
					tmp = d__[j] * l[j];
					work[indld - 1 + j] = tmp;
					work[indlld - 1 + j] = tmp * l[j];
					/* L50: */
				}
				if (ndepth > 0) {
					/*                 P and Q are index of the first and last eigenvalue to compute */
					/*                 within the current block */
					p = indexw[wbegin - 1 + oldfst];
					q = indexw[wbegin - 1 + oldlst];
					/*                 Offset for the arrays WORK, WGAP and WERR, i.e., th P-OFFSET */
					/*                 thru' Q-OFFSET elements of these arrays are to be used. */
					/*                  OFFSET = P-OLDFST */
					offset = indexw[wbegin] - 1;
					/*                 perform limited bisection (if necessary) to get approximate */
					/*                 eigenvalues to the precision needed. */
					pdlarrb(&in, &d__[ibegin], &work[indlld + ibegin - 1], &p,
							&q, rtol1, rtol2, &offset, &work[wbegin],
							&wgap[wbegin], &werr[wbegin], &work[indwrk],
							&iwork[iindwk], pivmin, &spdiam, &in, &iinfo);
					if (iinfo != 0) {
						*info = -1;
						//xiaohui add 2015-04-07
						free(dfst);
						free(lfst);
						free(dlst);
						free(llst);

						blk_begin++;
						blk_end++;
						free(blk_begin);
						free(blk_end);

						return 0;
					}
					/*                 We also recompute the extremal gaps. W holds all eigenvalues */
					/*                 of the unshifted matrix and must be used for computation */
					/*                 of WGAP, the entries of WORK might stem from RRRs with */
					/*                 different shifts. The gaps from WBEGIN-1+OLDFST to */
					/*                 WBEGIN-1+OLDLST are correctly computed in DLARRB. */
					/*                 However, we only allow the gaps to become greater since */
					/*                 this is what should happen when we decrease WERR */
					if (oldfst > 1) {
						/* Computing MAX */
						d__1 = wgap[wbegin + oldfst - 2], d__2 = w[wbegin
								+ oldfst - 1] - werr[wbegin + oldfst - 1]
								- w[wbegin + oldfst - 2]
								- werr[wbegin + oldfst - 2];
						wgap[wbegin + oldfst - 2] = max(d__1, d__2);
					}
					if (wbegin + oldlst - 1 < wend) {
						/* Computing MAX */
						d__1 = wgap[wbegin + oldlst - 1], d__2 = w[wbegin
								+ oldlst] - werr[wbegin + oldlst]
								- w[wbegin + oldlst - 1]
								- werr[wbegin + oldlst - 1];
						wgap[wbegin + oldlst - 1] = max(d__1, d__2);
					}
					/*                 Each time the eigenvalues in WORK get refined, we store */
					/*                 the newly found approximation with all shifts applied in W */
					i__3 = oldlst;
					for (j = oldfst; j <= i__3; ++j) {
						w[wbegin + j - 1] = work[wbegin + j - 1] + sigma;
						/* L53: */
					}
				}
				/*              Process the current node. */
				newfst = oldfst;
				i__3 = oldlst;
				for (j = oldfst; j <= i__3; ++j) {
//if(*dol==5)printf("ndepth=%d oldfst=%d oldlst=%d\n",ndepth,oldfst,oldlst);
					if (j == oldlst) {
						newlst = j;
					} else if (wgap[wbegin + j - 1]
							>= *minrgp
									* (d__1 = work[wbegin + j - 1], fabs(d__1))) {
						newlst = j;
					} else {
						continue;
					}
					if (newlst < *dol) {
						newfst = newlst + 1;
						continue;
					} else if (newfst > *dou) {
						break;
					}
//if(*dol==5)printf("newfst=%d newlst=%d\n",newfst,newlst);
					newsiz = newlst - newfst + 1;
					newftt = wbegin + newfst - 1;
					if (newsiz > 1) {
//if(*dol==5)printf("%d\n",newsiz);
						if (newfst == 1) {
							d__1 = 0., d__2 = w[wbegin] - werr[wbegin] - *vl;
							lgap = max(d__1, d__2);
						} else {
							lgap = wgap[wbegin + newfst - 2];
						}
						rgap = wgap[wbegin + newlst - 1];

						/*                    Compute left- and rightmost eigenvalue of child */
						/*                    to high precision in order to shift as close */
						/*                    as possible and obtain as large relative gaps */
						/*                    as possible */

						for (k = 1; k <= 2; ++k) {
							if (k == 1) {
								p = indexw[wbegin - 1 + newfst];
							} else {
								p = indexw[wbegin - 1 + newlst];
							}
							offset = indexw[wbegin] - 1;
							pdlarrb(&in, &d__[ibegin],
									&work[indlld + ibegin - 1], &p, &p, &rqtol,
									&rqtol, &offset, &work[wbegin],
									&wgap[wbegin], &werr[wbegin], &work[indwrk],
									&iwork[iindwk], pivmin, &spdiam, &in,
									&iinfo);
							/* L55: */
						}

//						if(*dol==5)printf("begin=%d\tfrom %d to %d\n",*dol,newfst,newlst);
						if (wbegin + newfst - 1 < *dol) {
							//printf("1\n");
							pdlarrf(&in, &d__[ibegin], &l[ibegin],
									&work[indld + ibegin - 1], &newfst, &newlst,
									&work[wbegin], &wgap[wbegin], &werr[wbegin],
									&spdiam, &lgap, &rgap, pivmin, &tau,
									&dfst[ibegin - 1], &lfst[ibegin - 1],
									&work[indwrk], &iinfo);
							//2014-07-02
							//free(dfst);
						} else if (wbegin + newlst - 1 > *dou) {
							//printf("2\n");
							pdlarrf(&in, &d__[ibegin], &l[ibegin],
									&work[indld + ibegin - 1], &newfst, &newlst,
									&work[wbegin], &wgap[wbegin], &werr[wbegin],
									&spdiam, &lgap, &rgap, pivmin, &tau,
									&dlst[ibegin - 1], &llst[ibegin - 1],
									&work[indwrk], &iinfo);
							//2014-07-02
							//free(dlst);
						} else {
							pdlarrf(&in, &d__[ibegin], &l[ibegin],
									&work[indld + ibegin - 1], &newfst, &newlst,
									&work[wbegin], &wgap[wbegin], &werr[wbegin],
									&spdiam, &lgap, &rgap, pivmin, &tau,
									&z__[ibegin + newftt * z_dim1],
									&z__[ibegin + (newftt + 1) * z_dim1],
									&work[indwrk], &iinfo);
						}
						if (iinfo == 0) {
//if(*dol==5)printf("new clus\n");
							/*                       a new RRR for the cluster was found by DLARRF */
							/*                       update shift and store it */
							ssigma = sigma + tau;
							if (wbegin + newfst - 1 < *dol) {
								lfst[iend - 1] = ssigma;
								//2014-07-02
								//free(lfst);
							} else if (wbegin + newlst - 1 > *dou) {
								llst[iend - 1] = ssigma;
								//2014-07-02
								//free(llst);
							} else {
								z__[iend + (newftt + 1) * z_dim1] = ssigma;
							}
							/*                       WORK() are the midpoints and WERR() the semi-width */
							/*                       Note that the entries in W are unchanged. */
							i__4 = newlst;
							for (k = newfst; k <= i__4; ++k) {
								fudge = eps * 3.
										* (d__1 = work[wbegin + k - 1], fabs(
												d__1));
								work[wbegin + k - 1] -= tau;
								fudge += eps * 4. * (d__1 =
										work[wbegin + k - 1], fabs(d__1));
								/*                          Fudge errors */
								werr[wbegin + k - 1] += fudge;
								/*                          Gaps are not fudged. Provided that WERR is small */
								/*                          when eigenvalues are close, a zero gap indicates */
								/*                          that a new representation is needed for resolving */
								/*                          the cluster. A fudge could lead to a wrong decision */
								/*                          of judging eigenvalues 'separated' which in */
								/*                          reality are not. This could have a negative impact */
								/*                          on the orthogonality of the computed eigenvectors. */
								/* L116: */
							}
							++nclus;
							k = newcls + (nclus << 1);
							iwork[k - 1] = newfst;
							iwork[k] = newlst;
						} else {
							//printf("dlarrf \n");
							*info = -2;
							//xiaohui add 2015-04-07
							free(dfst);
							free(lfst);
							free(dlst);
							free(llst);

							blk_begin++;
							blk_end++;
							free(blk_begin);
							free(blk_end);

							return 0;
						}
					} else {

						/*                    Compute eigenvector of singleton */

						iter = 0;

						tol = log((double) in) * 4. * eps;

						k = newfst;
						windex = wbegin + k - 1;
						/* Computing MAX */
						i__4 = windex - 1;
						windmn = max(i__4, *il);
						//windmn = max(i__4, 1);
						/* Computing MIN */
						i__4 = windex + 1;
						windpl = min(i__4, *iu);
						//windpl = min(i__4, *m);
						lambda = work[windex];
						++done;
						/*                    Check if eigenvector computation is to be skipped */
						if (windex < *dol || windex > *dou) {
							eskip = true
							;
							goto L125;
						} else {
							eskip = false
							;
						}
						left = work[windex] - werr[windex];
						right = work[windex] + werr[windex];
						indeig = indexw[windex];
						/*                    Note that since we compute the eigenpairs for a child, */
						/*                    all eigenvalue approximations are w.r.t the same shift. */
						/*                    In this case, the entries in WORK should be used for */
						/*                    computing the gaps since they exhibit even very small */
						/*                    differences in the eigenvalues, as opposed to the */
						/*                    entries in W which might "look" the same. */
						if (k == 1) {
							/*                       In the case RANGE='I' and with not much initial */
							/*                       accuracy in LAMBDA and VL, the formula */
							/*                       LGAP = MAX( ZERO, (SIGMA - VL) + LAMBDA ) */
							/*                       can lead to an overestimation of the left gap and */
							/*                       thus to inadequately early RQI 'convergence'. */
							/*                       Prevent this by forcing a small left gap. */
							/* Computing MAX */
							d__1 = fabs(left), d__2 = fabs(right);
							lgap = eps * max(d__1, d__2);
						} else {
							lgap = wgap[windmn];
						}
						if (k == im) {
							/*                       In the case RANGE='I' and with not much initial */
							/*                       accuracy in LAMBDA and VU, the formula */
							/*                       can lead to an overestimation of the right gap and */
							/*                       thus to inadequately early RQI 'convergence'. */
							/*                       Prevent this by forcing a small right gap. */
							/* Computing MAX */
							d__1 = fabs(left), d__2 = fabs(right);
							rgap = eps * max(d__1, d__2);
						} else {
							rgap = wgap[windex];
						}
						gap = min(lgap, rgap);
						if (k == 1 || k == im) {
							/*                       The eigenvector support can become wrong */
							/*                       because significant entries could be cut off due to a */
							/*                       large GAPTOL parameter in LAR1V. Prevent this. */
							gaptol = 0.;
						} else {
							gaptol = gap * eps;
						}
						isupmn = in;
						isupmx = 1;
						/*                    Update WGAP so that it holds the minimum gap */
						/*                    to the left or the right. This is crucial in the */
						/*                    case where bisection is used to ensure that the */
						/*                    eigenvalue is refined up to the required precision. */
						/*                    The correct value is restored afterwards. */
						savgap = wgap[windex];
						wgap[windex] = gap;
						/*                    We want to use the Rayleigh Quotient Correction */
						/*                    as often as possible since it converges quadratically */
						/*                    when we are close enough to the desired eigenvalue. */
						/*                    However, the Rayleigh Quotient can have the wrong sign */
						/*                    and lead us away from the desired eigenvalue. In this */
						/*                    case, the best we can do is to use bisection. */
						usedbs = false
						;
						usedrq = false
						;
						/*                    Bisection is initially turned off unless it is forced */
						needbs = !tryrqc;
						L120:
						/*                    Check if bisection should be used to refine eigenvalue */
						if (needbs) {
							/*                       Take the bisection as new iterate */
							usedbs = true
							;
							itmp1 = iwork[iindr + windex];
							offset = indexw[wbegin] - 1;
							d__1 = eps * 2.;
							pdlarrb(&in, &d__[ibegin],
									&work[indlld + ibegin - 1], &indeig,
									&indeig, &c_b5, &d__1, &offset,
									&work[wbegin], &wgap[wbegin], &werr[wbegin],
									&work[indwrk], &iwork[iindwk], pivmin,
									&spdiam, &itmp1, &iinfo);
							if (iinfo != 0) {
								*info = -3;
								//xiaohui add 2015-04-07
								free(dfst);
								free(lfst);
								free(dlst);
								free(llst);

								blk_begin++;
								blk_end++;
								free(blk_begin);
								free(blk_end);

								return 0;
							}
							lambda = work[windex];
							/*                       Reset twist index from inaccurate LAMBDA to */
							/*                       force computation of true MINGMA */
							iwork[iindr + windex] = 0;
						}
						/*                    Given LAMBDA, compute the eigenvector. */
						L__1 = !usedbs;
						pdlar1v(&in, &c__1, &in, &lambda, &d__[ibegin],
								&l[ibegin], &work[indld + ibegin - 1],
								&work[indlld + ibegin - 1], pivmin, &gaptol,
								&z__[ibegin + windex * z_dim1], &L__1, &negcnt,
								&ztz, &mingma, &iwork[iindr + windex],
								&isuppz[(windex << 1) - 1], &nrminv, &resid,
								&rqcorr, &work[indwrk]);
						if (iter == 0) {
							bstres = resid;
							bstw = lambda;
						} else if (resid < bstres) {
							bstres = resid;
							bstw = lambda;
						}
						/* Computing MIN */
						i__4 = isupmn, i__5 = isuppz[(windex << 1) - 1];
						isupmn = min(i__4, i__5);
						/* Computing MAX */
						i__4 = isupmx, i__5 = isuppz[windex * 2];
						isupmx = max(i__4, i__5);
						++iter;
						/*                    sin alpha <= |resid|/gap */
						/*                    Note that both the residual and the gap are */
						/*                    proportional to the matrix, so ||T|| doesn't play */
						/*                    a role in the quotient */

						/*                    Convergence test for Rayleigh-Quotient iteration */
						/*                    (omitted when Bisection has been used) */

						if (resid > tol * gap
								&& fabs(rqcorr) > rqtol * fabs(lambda)
								&& !usedbs) {
							/*                       We need to check that the RQCORR update doesn't */
							/*                       move the eigenvalue away from the desired one and */
							/*                       towards a neighbor. -> protection with bisection */
							if (indeig <= negcnt) {
								/*                          The wanted eigenvalue lies to the left */
								sgndef = -1.;
							} else {
								/*                          The wanted eigenvalue lies to the right */
								sgndef = 1.;
							}
							/*                       We only use the RQCORR if it improves the */
							/*                       the iterate reasonably. */
							if (rqcorr * sgndef >= 0.
									&& lambda + rqcorr <= right
									&& lambda + rqcorr >= left) {
								usedrq = true
								;
								/*                          Store new midpoint of bisection interval in WORK */
								if (sgndef == 1.) {
									/*                             The current LAMBDA is on the left of the true */
									/*                             eigenvalue */
									left = lambda;
									/*                             We prefer to assume that the error estimate */
									/*                             is correct. We could make the interval not */
									/*                             as a bracket but to be modified if the RQCORR */
									/*                             chooses to. In this case, the RIGHT side should */
									/*                             be modified as follows: */
									/*                              RIGHT = MAX(RIGHT, LAMBDA + RQCORR) */
								} else {
									/*                             The current LAMBDA is on the right of the true */
									/*                             eigenvalue */
									right = lambda;
									/*                             See comment about assuming the error estimate is */
									/*                             correct above. */
									/*                              LEFT = MIN(LEFT, LAMBDA + RQCORR) */
								}
								work[windex] = (right + left) * .5;
								/*                          Take RQCORR since it has the correct sign and */
								/*                          improves the iterate reasonably */
								lambda += rqcorr;
								/*                          Update width of error interval */
								werr[windex] = (right - left) * .5;
							} else {
								needbs = true
								;
							}
							if (right - left < rqtol * fabs(lambda)) {
								/*                             The eigenvalue is computed to bisection accuracy */
								/*                             compute eigenvector and stop */
								usedbs = true
								;
								goto L120;
							} else if (iter < 10) {
								goto L120;
							} else if (iter == 10) {
								needbs = true
								;
								goto L120;
							} else {
								*info = 5;
								//xiaohui add 2015-04-07
								free(dfst);
								free(lfst);
								free(dlst);
								free(llst);

								blk_begin++;
								blk_end++;
								free(blk_begin);
								free(blk_end);

								return 0;
							}
						} else {
							stp2ii = false
							;
							if (usedrq && usedbs && bstres <= resid) {
								lambda = bstw;
								stp2ii = true
								;
							}
							if (stp2ii) {
								/*                          improve error angle by second step */
								L__1 = !usedbs;
								pdlar1v(&in, &c__1, &in, &lambda, &d__[ibegin],
										&l[ibegin], &work[indld + ibegin - 1],
										&work[indlld + ibegin - 1], pivmin,
										&gaptol, &z__[ibegin + windex * z_dim1],
										&L__1, &negcnt, &ztz, &mingma,
										&iwork[iindr + windex],
										&isuppz[(windex << 1) - 1], &nrminv,
										&resid, &rqcorr, &work[indwrk]);
							}
							work[windex] = lambda;
						}

						/*                    Compute FP-vector support w.r.t. whole matrix */

						isuppz[(windex << 1) - 1] += oldien;
						isuppz[windex * 2] += oldien;
						zfrom = isuppz[(windex << 1) - 1];
						zto = isuppz[windex * 2];
						isupmn += oldien;
						isupmx += oldien;
						/*                    Ensure vector is ok if support in the RQI has changed */
						if (isupmn < zfrom) {
							i__4 = zfrom - 1;
							for (ii = isupmn; ii <= i__4; ++ii) {
								z__[ii + windex * z_dim1] = 0.;
								/* L122: */
							}
						}
						if (isupmx > zto) {
							i__4 = isupmx;
							for (ii = zto + 1; ii <= i__4; ++ii) {
								z__[ii + windex * z_dim1] = 0.;
								/* L123: */
							}
						}
						i__4 = zto - zfrom + 1;
						pdscal(&i__4, &nrminv, &z__[zfrom + windex * z_dim1],
								&c__1);
						L125:
						/*                    Update W */
						w[windex] = lambda + sigma;
						/*                    Recompute the gaps on the left and right */
						/*                    But only allow them to become larger and not */
						/*                    smaller (which can only happen through "bad" */
						/*                    cancellation and doesn't reflect the theory */
						/*                    where the initial gaps are underestimated due */
						/*                    to WERR being too crude.) */
						if (!eskip) {
							if (k > 1) {
								/* Computing MAX */
								d__1 = wgap[windmn], d__2 = w[windex]
										- werr[windex] - w[windmn]
										- werr[windmn];
								wgap[windmn] = max(d__1, d__2);
							}
							if (windex < wend) {
								/* Computing MAX */
								d__1 = savgap, d__2 = w[windpl] - werr[windpl]
										- w[windex] - werr[windex];
								wgap[windex] = max(d__1, d__2);
							}
						}
						++idone;
					}
					/*                 here ends the code for the current child */

					L139:
					/*                 Proceed to any remaining child nodes */
					newfst = j + 1;
					if (newfst > *dou) {
						break;
					}
					L140: ;
				}
				/* L150: */
			}
//if(*dol==5)printf("dp+\n");
			++ndepth;
			goto L40;
		}
		ibegin = iend + 1;
		wbegin = wend + 1;
		L170: ;
	}
	
	//xiaohui add 2015-04-07
	free(dfst);
	free(lfst);
	free(dlst);
	free(llst);

	blk_begin++;
	blk_end++;
	free(blk_begin);
	free(blk_end);

	return 0;

	/*     End of DLARRV */

} /* dlarrv_ */
