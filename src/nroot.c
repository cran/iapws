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

#include <Rconfig.h> // for HAVE_ALLOCA_H (R >= 3.2.2)
#ifdef __GNUC__
#	undef alloca
#	define alloca(x) __builtin_alloca((x))
#elif defined(HAVE_ALLOCA_H)
#	include <alloca.h>
#endif

#include <assert.h>
#include <math.h>

#include <R_ext/Print.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "nroot.h"

struct nroot_control nroot_default = {
	.trace = 0,
	.maxit = 100,
	.abstol = 1e-9,
	.reltol = 1e-9
};

void nroot_log(const struct nroot_control *ctrl)
{
	if (ctrl->trace > 0)
		REprintf("nroot: iter=%d epsf=%.8e epsx=%.8e\n",
			ctrl->maxit, ctrl->abstol, ctrl->reltol);
}

enum nroot_exit nroot1(nroot_fun fun, double *x, void *prms,
		struct nroot_control *ctrl)
{
	assert(ctrl != &nroot_default);
	const double tolf = ctrl->abstol;
	const double tolx = ctrl->reltol;
	double fx, dfx, dx;
	while (ctrl->maxit--) {
		fun(x, prms, &fx, &dfx);
		ctrl->abstol = fabs(fx);
		if (ctrl->abstol <= tolf) return NROOT_SUCCESS;

		if (dfx == 0.0) return NROOT_ZERODET;

		dx = -fx / dfx;
		ctrl->reltol = fabs(dx) / fabs(*x);
		if (ctrl->reltol <= tolx) return NROOT_SUCCESS;

		nroot_log(ctrl);

		*x += dx;
	}
	return NROOT_MAXITER;
}

enum nroot_exit nroot2(nroot_fun fun, double *x, void *prms,
		struct nroot_control *ctrl)
{
	assert(ctrl != &nroot_default);
	const double tolf = ctrl->abstol;
	const double tolx = ctrl->reltol;
	double fx[2], dfx[4], dx[2];
	double det;

	while (ctrl->maxit--) {
		fun(x, prms, fx, dfx);
		ctrl->abstol = fabs(fx[0]) + fabs(fx[1]);
		if (ctrl->abstol <= tolf) return NROOT_SUCCESS;

		det = dfx[0] * dfx[3] - dfx[2] * dfx[1];
		if (det == 0.0) return NROOT_ZERODET;

		det = 1.0 / det;
		dx[0] = -(dfx[3] * fx[0] - dfx[2] * fx[1]) * det;
		dx[1] = +(dfx[1] * fx[0] - dfx[0] * fx[1]) * det;
		ctrl->reltol = (fabs(dx[0]) + fabs(dx[1])) /
			(fabs(x[0]) + fabs(x[1]));
		if (ctrl->reltol <= tolx) return NROOT_SUCCESS;

		nroot_log(ctrl);

		x[0] += dx[0];
		x[1] += dx[1];
	}
	return NROOT_MAXITER;
}

enum nroot_exit nrootn(int n, nroot_fun fun, double *x, void *prms,
		struct nroot_control *ctrl)
{
	assert(ctrl != &nroot_default);
	const double tolf = ctrl->abstol;
	const double tolx = ctrl->reltol;
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

	while (ctrl->maxit--) {
		fun(x, prms, fx, dfx);
		ctrl->abstol = F77_NAME(dasum)(&n, fx, &i1);
		if (ctrl->abstol <= tolf) return NROOT_SUCCESS;

		F77_NAME(dscal)(&n, &dm1, fx, &i1);
		F77_NAME(dgesv)(&n, &i1, dfx, &n, ipiv, fx, &n, &info);
		if (info != 0) return NROOT_ZERODET;

		ctrl->reltol = F77_NAME(dasum)(&n, fx, &i1) /
			F77_NAME(dasum)(&n, x, &i1);
		if (ctrl->reltol <= tolx) return NROOT_SUCCESS;

		nroot_log(ctrl);

		F77_NAME(daxpy)(&n, &d1, fx, &i1, x, &i1);
	}

	return NROOT_MAXITER;
}

enum nroot_exit sroot(nroot_fun fun, double *x, void *prms,
		struct nroot_control *ctrl)
{
	assert(ctrl != &nroot_default);
	const double tolf = ctrl->abstol;
	const double tolx = ctrl->reltol;
	double fx, dfx, dx, fx0;

	fun(x, prms, &fx0, &dfx);
	ctrl->abstol = fabs(fx0);
	if (ctrl->abstol <= tolf) return NROOT_SUCCESS;
	dx = tolx;
	*x += dx;

	while (ctrl->maxit--) {
		fun(x, prms, &fx, &dfx);
		ctrl->abstol = fabs(fx);
		if (ctrl->abstol <= tolf) return NROOT_SUCCESS;

		dfx = fx - fx0;
		if (dfx == 0.0) return NROOT_ZERODET;

		dx *= -fx / dfx;
		ctrl->reltol = fabs(dx) / fabs(*x);
		if (ctrl->reltol <= tolx) return NROOT_SUCCESS;

		nroot_log(ctrl);

		*x += dx;
		fx0 = fx;
	}
	return NROOT_MAXITER;
}

