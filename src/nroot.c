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

int nroot(root_fun fun, double *x, void *prms,
		double *tolf, double *tolx, int *maxiter)
{
	const double epsf = *tolf, epsx = *tolx;
	double fx, dfx, dx;
	while ((*maxiter)--) {
		fun(x, prms, &fx, &dfx);
		*tolf = fabs(fx);
		if (*tolf <= epsf) return 0;

		if (dfx == 0.0) return -2;

		dx = -fx / dfx;
		*tolx = fabs(dx);
		if (*tolx <= epsx) return 0;
#ifndef NDEBUG
		fprintf(stderr, "nroot: iter=%d epsf=%.8e epsx=%.8e\n",
				*maxiter, *tolf, *tolx);
#endif
		*x += dx;
	}
	return -1;
}

int nroot2(root_fun fun, double *x, void *prms,
		double *tolf, double *tolx, int *maxiter)
{
	const double epsf = *tolf, epsx = *tolx;
	double fx[2], dfx[4], dx[2];
	double det;

	while ((*maxiter)--) {
		fun(x, prms, fx, dfx);
		*tolf = fabs(fx[0]) + fabs(fx[1]);
		if (*tolf <= epsf) return 0;

		det = dfx[0] * dfx[3] - dfx[2] * dfx[1];
		if (det == 0.0) return -2;

		dx[0] = -(dfx[3] * fx[0] - dfx[2] * fx[1]) / det;
		dx[1] = +(dfx[1] * fx[0] - dfx[0] * fx[1]) / det;
		*tolx = fabs(dx[0]) + fabs(dx[1]);
		if (*tolx <= epsx) return 0;
#ifndef NDEBUG
		fprintf(stderr, "nroot: iter=%d epsf=%.8e epsx=%.8e\n",
				*maxiter, *tolf, *tolx);
#endif
		x[0] += dx[0];
		x[1] += dx[1];
	}
	return -1;
}

int nrootn(int n, root_fun fun, double *x, void *prms,
		double *tolf, double *tolx, int *maxiter)
{
	const double epsf = *tolf, epsx = *tolx;
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

	while ((*maxiter)--) {
		fun(x, prms, fx, dfx);
		*tolf = F77_NAME(dasum)(&n, fx, &i1);
		if (*tolf <= epsf) return 0;

		F77_NAME(dscal)(&n, &dm1, fx, &i1);
		F77_NAME(dgesv)(&n, &i1, dfx, &n, ipiv, fx, &n, &info);
		if (info != 0) return -2;

		*tolx = F77_NAME(dasum)(&n, fx, &i1);
		if (*tolx <= epsx) return 0;

		F77_NAME(daxpy)(&n, &d1, fx, &i1, x, &i1);
	}

	return -1;
}
