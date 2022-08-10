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

#include "R_interface.h"

#include "iapws.h"
#include "iapws95.h"
#include "if97.h"
#include "visc.h"
#include "cond.h"

#if 0
static double (*iapws_a[15])(iapws_phi *phi) = {
	iapws_f, iapws_g,
	iapws_u, iapws_h,
	iapws_s, iapws_t,
	iapws_p, iapws_v,
	iapws_cp, iapws_cv,
	iapws_w, iapws_rho,
	iapws_chit,
	NULL, /* eta */
	NULL, /* lambda */
};
#endif

/* IF97 */

static double (*if97_a[18])(iapws_phi *phi) = {
	iapws_f, iapws_g,
	iapws_u, iapws_h,
	iapws_s, iapws_t,
	iapws_p, iapws_v,
	iapws_cp, iapws_cv,
	iapws_w, iapws_rho,
	iapws_alpha, iapws_beta,
	iapws_chit,
	if97_eta, if97_lambda,
	iapws_sigma,
};

MAKE_R_FUN_2(if97_state, double, double, int)
MAKE_R_FUN_2(if97_region, double, double, int)
MAKE_R_FUN_1(if97_tsat, double, double)
MAKE_R_FUN_1(if97_psat, double, double)

SEXP R_if97(SEXP w, SEXP p, SEXP t, SEXP s)
{
	int iw, nw = length(w);
	int ip, np = length(p);
	int it, nt = length(t);
	int is, ns = length(s);
	int i, n = (np > nt ? np : nt);
	if (n < ns) n = ns;

	SEXP d = PROTECT(allocMatrix(REALSXP, n, nw));
	int *xw = INTEGER(w);
	double *xp = REAL(p);
	double *xt = REAL(t);
	int *xs = INTEGER(s);
	double *xd = REAL(d);

	iapws_phi gamma;
	int err, offset;
	MOD_ITERATE3(n, np, nt, ns, i, ip, it, is,
		if ((err = if97_gamma(xp[ip], xt[it], xs[is], &gamma)) == 0) {
			for (iw = 0, offset = 0; iw < nw; iw++, offset += n) {
				xd[i + offset] = if97_a[xw[iw]](&gamma);
			}
		} else {
			for (iw = 0, offset = 0; iw < nw; iw++, offset += n) {
				xd[i + offset] = NA_REAL;
			}
		}
	);

	UNPROTECT(1);
	return d;
}

/* IAPWS95 */

static double (*iapws95_a[18])(iapws_phi *phi) = {
	iapws_f, iapws_g,
	iapws_u, iapws_h,
	iapws_s, iapws_t,
	iapws_p, iapws_v,
	iapws_cp, iapws_cv,
	iapws_w, iapws_rho,
	iapws_alpha, iapws_beta,
	iapws_chit,
	iapws95_eta, iapws95_lambda,
	iapws_sigma,
};

MAKE_R_FUN_2(iapws95_state, double, double, int)

SEXP R_iapws95(SEXP w, SEXP r, SEXP t)
{
	int iw, nw = length(w);
	int ir, nr = length(r);
	int it, nt = length(t);
	int i, n = (nr > nt ? nr : nt);

	SEXP d = PROTECT(allocMatrix(REALSXP, n, nw));
	int *xw = INTEGER(w);
	double *xr = REAL(r);
	double *xt = REAL(t);
	double *xd = REAL(d);

	iapws_phi phi;
	int err, offset;
	MOD_ITERATE2(n, nr, nt, i, ir, it,
		if ((err = iapws95_phi(xr[ir], xt[it], &phi)) == 0) {
			for (iw = 0, offset = 0; iw < nw; iw++, offset += n) {
				xd[i + offset] = iapws95_a[xw[iw]](&phi);
			}
		} else {
			for (iw = 0, offset = 0; iw < nw; iw++, offset += n) {
				xd[i + offset] = NA_REAL;
			}
		}
	);

	UNPROTECT(1);
	return d;
}

SEXP R_iapws95_pt(SEXP w, SEXP p, SEXP t, SEXP s)
{
	int iw, nw = length(w);
	int ip, np = length(p);
	int it, nt = length(t);
	int is, ns = length(s);
	int i, n = (np > nt ? np : nt);
	if (n < ns) n = ns;

	SEXP d = PROTECT(allocMatrix(REALSXP, n, nw));
	int *xw = INTEGER(w);
	double *xp = REAL(p);
	double *xt = REAL(t);
	int *xs = INTEGER(s);
	double *xd = REAL(d);

	iapws_phi phi;
	int err, offset;
	MOD_ITERATE3(n, np, nt, ns, i, ip, it, is,
		if ((err = iapws95_phi_pt(xp[ip], xt[it], xs[is], &phi)) == 0) {
			for (iw = 0, offset = 0; iw < nw; iw++, offset += n) {
				xd[i + offset] = iapws95_a[xw[iw]](&phi);
			}
		} else {
			for (iw = 0, offset = 0; iw < nw; iw++, offset += n) {
				xd[i + offset] = NA_REAL;
			}
		}
	);

	UNPROTECT(1);
	return d;
}

SEXP R_iapws95_sat(SEXP w, SEXP t)
{
	int iw, nw = length(w);
	int it, nt = length(t);
	int ntw = nt * nw;  /* overflow */

	SEXP d = PROTECT(alloc3DArray(REALSXP, nt, nw, 2));
	int *xw = INTEGER(w);
	double *xt = REAL(t);
	double *xd = REAL(d);

	iapws_phi phil, phig;
	int err, offset;
	R_ITERATE(nt, it,
		if ((err = iapws95_sat(xt[it], &phil, &phig)) == 0) {
			for (iw = 0, offset = 0; iw < nw; iw++, offset += nt) {
				xd[it + offset] = iapws95_a[xw[iw]](&phil);
				xd[it + offset + ntw] = iapws95_a[xw[iw]](&phig);
			}
		} else {
			for (iw = 0, offset = 0; iw < nw; iw++, offset += nt) {
				xd[it + offset] = NA_REAL;
				xd[it + offset + ntw] = NA_REAL;
			}
		}
	);

	UNPROTECT(1);
	return d;
}

SEXP R_iapws95_sat_p(SEXP w, SEXP t)
{
	int iw, nw = length(w);
	int it, nt = length(t);
	int ntw = nt * nw;  /* overflow */

	SEXP d = PROTECT(alloc3DArray(REALSXP, nt, nw, 2));
	int *xw = INTEGER(w);
	double *xt = REAL(t);
	double *xd = REAL(d);

	iapws_phi phil, phig;
	int err, offset;
	R_ITERATE(nt, it,
		if ((err = iapws95_sat_p(xt[it], &phil, &phig)) == 0) {
			for (iw = 0, offset = 0; iw < nw; iw++, offset += nt) {
				xd[it + offset] = iapws95_a[xw[iw]](&phil);
				xd[it + offset + ntw] = iapws95_a[xw[iw]](&phig);
			}
		} else {
			for (iw = 0, offset = 0; iw < nw; iw++, offset += nt) {
				xd[it + offset] = NA_REAL;
				xd[it + offset + ntw] = NA_REAL;
			}
		}
	);

	UNPROTECT(1);
	return d;
}
