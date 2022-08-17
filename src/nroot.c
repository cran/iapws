/*
 * iapws - IAPWS formulations for the properties of water and steam
 * Copyright (C) 2022 Jonathan Debove
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef NDEBUG
#	include <stdio.h>
#endif

#ifdef _WIN32
#	include <malloc.h>
#	define alloca _alloca
#else
#	include <alloca.h>
#endif

#include <math.h>

#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "nroot.h"

int nroot(root_fun fun, double *x, void *prms, double *errf,
		double tol, int maxiter)
{
	int iter;
	double fx, dfx;
	for (iter = 0; iter < maxiter; ++iter) {
		fun(x, prms, &fx, &dfx);
		*errf = fabs(fx);
#ifndef NDEBUG
		fprintf(stderr, "nroot: i=%d x=%.8f fx=%.8e df/dx=%.8e\n",
				iter, *x, fx, dfx);
#endif
		if (*errf <= tol) {
			break;
		}
		if (dfx == 0) {
			iter = -99;
			break;
		}
		*x -= fx / dfx;
	}
	return iter;
}

int nroot2(root_fun fun, double *x, void *prms,
		double *errf, double tol, int maxiter)
{
	int iter;
	double fx[2], dfx[4];
	double det;

	for (iter = 0; iter < maxiter; ++iter) {
		fun(x, prms, fx, dfx);
		*errf = fabs(fx[0]) + fabs(fx[1]);
#ifndef NDEBUG
		fprintf(stderr, "nroot: iter=%d errf=%.8e\n", iter, *errf);
#endif
		if (*errf <= tol) {
			break;
		}
		det = dfx[0] * dfx[3] - dfx[2] * dfx[1];
		if (det == 0.0) {
			iter = -99;
			break;
		}
		x[0] -= (dfx[3] * fx[0] - dfx[2] * fx[1]) / det;
		x[1] += (dfx[1] * fx[0] - dfx[0] * fx[1]) / det;
	}
	return iter;
}

int nrootn(int n, root_fun fun, double *x, void *prms,
		double *errf, double tol, int maxiter)
{
	int iter = -88;
	double *fx, *dfx;
	/* Blas & Lapack */
	int i1 = 1;
	double d1 = 1.0;
	double dm1 = -1.0;
	int *ipiv;
	int info;

	fx = alloca(sizeof(fx) * n);
	dfx = alloca(sizeof(dfx) * n * n);
	ipiv = alloca(sizeof(ipiv) * n);

	for (iter = 0; iter < maxiter; ++iter) {
		fun(x, prms, fx, dfx);
		*errf = F77_NAME(dasum)(&n, fx, &i1);
#ifndef NDEBUG
		fprintf(stderr, "nroot: iter=%d errf=%.8e\n", iter, *errf);
#endif
		if (*errf <= tol) {
			break;
		}
		F77_NAME(dscal)(&n, &dm1, fx, &i1);
		F77_NAME(dgesv)(&n, &i1, dfx, &n, ipiv, fx, &n, &info);
		if (info != 0) {
			iter = -99;
			break;
		}
		F77_NAME(daxpy)(&n, &d1, fx, &i1, x, &i1);
	}

	return iter;
}
