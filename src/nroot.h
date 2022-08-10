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

#define USE_LAPACK

typedef void root_fun(double *x, void *prms, double *fx, double *dfx);
int nroot(root_fun fun, double *x, void *prms,
		double *errf, double tol, int maxiter);
int nroot2(root_fun fun, double *x, void *prms,
		double *errf, double tol, int maxiter);
#ifdef USE_LAPACK
int nrootn(int n, root_fun fun, double *x, void *prms,
		double *errf, double tol, int maxiter);
#endif

#endif
