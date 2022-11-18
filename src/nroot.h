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

#ifndef IAPWS_NROOT_H
#define IAPWS_NROOT_H

enum nroot_exit {
	NROOT_MAXITER = -3,
	NROOT_ZERODET = -2,
	NROOT_FAILURE = -1,
	NROOT_SUCCESS =  0,
};

struct nroot_control {
	int trace;
	int maxit;
	double abstol;
	double reltol;
};
extern struct nroot_control nroot_default;

typedef void nroot_fun(double *x, void *prms, double *fx, double *dfx);

enum nroot_exit nroot1(nroot_fun fun, double *x, void *prms,
		struct nroot_control *ctrl);
enum nroot_exit nroot2(nroot_fun fun, double *x, void *prms,
		struct nroot_control *ctrl);
enum nroot_exit nrootn(int n, nroot_fun fun, double *x, void *prms,
		struct nroot_control *ctrl);

enum nroot_exit sroot(nroot_fun fun, double *x, void *prms,
		struct nroot_control *ctrl);

#endif
