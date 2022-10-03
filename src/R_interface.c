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

#include <R_ext/Rdynload.h>

#include "R_interface.h"
#include "iapws.h"
#include "iapws95.h"
#include "if97.h"
#include "ice06.h"
#include "visc.h"
#include "cond.h"
#include "surf.h"
#include "dielec.h"
#include "nroot.h"

/* IF97 */

static double (*if97_a[])(const iapws_phi *phi) = {
	iapws_f, iapws_g,
	iapws_u, iapws_h,
	iapws_s, iapws_t,
	iapws_p, iapws_v,
	iapws_cp, iapws_cv,
	iapws_w, iapws_rho,
	iapws_alpha, iapws_beta,
	iapws_kappat,
	if97_eta, if97_lambda,
	iapws_sigma, iapws_epsilon,
};

#define REGISTER_IF97_REGION 0

#if REGISTER_IF97_REGION
MAKE_R_FUN_2(if97_region_pt, double, double, int)
MAKE_R_FUN_2(if97_region_ph, double, double, int)
#endif
MAKE_R_FUN_2(if97_state_pt, double, double, int)
MAKE_R_FUN_1(if97_tsat, double, double)
MAKE_R_FUN_1(if97_psat, double, double)

SEXP R_if97_pt(SEXP w, SEXP p, SEXP t, SEXP s)
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
	int offset;
	MOD_ITERATE3_CHECK(NINTERRUPT, n, np, nt, ns, i, ip, it, is,
		if (if97_gamma_pt(xp[ip], xt[it], xs[is], &gamma) == 0) {
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

SEXP R_if97_ph(SEXP w, SEXP p, SEXP h)
{
	int iw, nw = length(w);
	int ip, np = length(p);
	int ih, nh = length(h);
	int i, n = (np > nh ? np : nh);

	SEXP d = PROTECT(allocMatrix(REALSXP, n, nw));
	int *xw = INTEGER(w);
	double *xp = REAL(p);
	double *xh = REAL(h);
	double *xd = REAL(d);

	iapws_phi gamma;
	int offset;
	MOD_ITERATE2_CHECK(NINTERRUPT, n, np, nh, i, ip, ih,
		if (if97_gamma_ph(xp[ip], xh[ih], &gamma) == 0) {
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

static double (*iapws95_a[])(const iapws_phi *phi) = {
	iapws_f, iapws_g,
	iapws_u, iapws_h,
	iapws_s, iapws_t,
	iapws_p, iapws_v,
	iapws_cp, iapws_cv,
	iapws_w, iapws_rho,
	iapws_alpha, iapws_beta,
	iapws_kappat,
	iapws95_eta, iapws95_lambda,
	iapws_sigma, iapws_epsilon,
};

MAKE_R_FUN_2(iapws95_state_pt, double, double, int)
MAKE_R_FUN_2(iapws95_state_rhot, double, double, int)

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
	int offset;
	MOD_ITERATE2_CHECK(NINTERRUPT, n, nr, nt, i, ir, it,
		iapws95_phi(xr[ir], xt[it], &phi);
		for (iw = 0, offset = 0; iw < nw; iw++, offset += n) {
			xd[i + offset] = iapws95_a[xw[iw]](&phi);
		}
	);

	UNPROTECT(1);
	return d;
}

SEXP R_iapws95_rhot(SEXP w, SEXP r, SEXP t, SEXP s)
{
	int iw, nw = length(w);
	int ir, nr = length(r);
	int it, nt = length(t);
	int is, ns = length(s);
	int i, n = (nr > nt ? nr : nt);
	if (n < ns) n = ns;

	SEXP d = PROTECT(allocMatrix(REALSXP, n, nw));
	int *xw = INTEGER(w);
	double *xr = REAL(r);
	double *xt = REAL(t);
	int *xs = INTEGER(s);
	double *xd = REAL(d);

	iapws_phi phi;
	int offset;
	MOD_ITERATE3_CHECK(NINTERRUPT, n, nr, nt, ns, i, ir, it, is,
		if (iapws95_phi_rhot(xr[ir], xt[it], xs[is], &phi) == 0) {
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
	int offset;
	MOD_ITERATE3_CHECK(NINTERRUPT, n, np, nt, ns, i, ip, it, is,
		if (iapws95_phi_pt(xp[ip], xt[it], xs[is], &phi) == 0) {
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

SEXP R_iapws95_ph(SEXP w, SEXP p, SEXP h)
{
	int iw, nw = length(w);
	int ip, np = length(p);
	int ih, nh = length(h);
	int i, n = (np > nh ? np : nh);

	SEXP d = PROTECT(allocMatrix(REALSXP, n, nw));
	int *xw = INTEGER(w);
	double *xp = REAL(p);
	double *xh = REAL(h);
	double *xd = REAL(d);

	iapws_phi phi;
	int offset;
	MOD_ITERATE2_CHECK(NINTERRUPT, n, np, nh, i, ip, ih,
		if (iapws95_phi_ph(xp[ip], xh[ih], &phi) == 0) {
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

SEXP R_iapws95_sat_t(SEXP w, SEXP t)
{
	int iw, nw = length(w);
	int it, nt = length(t);
	int ntw = nt * nw;  /* overflow */

	SEXP d = PROTECT(alloc3DArray(REALSXP, nt, nw, 2));
	int *xw = INTEGER(w);
	double *xt = REAL(t);
	double *xd = REAL(d);

	iapws_phi phil, phig;
	int offset;
	R_ITERATE_CHECK(NINTERRUPT, nt, it,
		if (iapws95_sat_t(xt[it], &phil, &phig) == 0) {
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
	int offset;
	R_ITERATE_CHECK(NINTERRUPT, nt, it,
		if (iapws95_sat_p(xt[it], &phil, &phig) == 0) {
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

/* ICE06 */

static double (*ice06_a[])(const iapws_phi *phi) = {
	iapws_f, iapws_g,
	iapws_u, iapws_h,
	iapws_s, iapws_t,
	iapws_p, iapws_v,
	iapws_cp, iapws_cv,
	iapws_w, iapws_rho,
	iapws_alpha, iapws_beta,
	iapws_kappat,
	NULL, NULL, NULL, NULL,
};

SEXP R_ice06_pt(SEXP w, SEXP p, SEXP t)
{
	int iw, nw = length(w);
	int ip, np = length(p);
	int it, nt = length(t);
	int i, n = (np > nt ? np : nt);

	SEXP d = PROTECT(allocMatrix(REALSXP, n, nw));
	int *xw = INTEGER(w);
	double *xp = REAL(p);
	double *xt = REAL(t);
	double *xd = REAL(d);

	iapws_phi gamma;
	int offset;
	MOD_ITERATE2_CHECK(NINTERRUPT, n, np, nt, i, ip, it,
		ice06_gamma(xp[ip], xt[it], &gamma);
		for (iw = 0, offset = 0; iw < nw; iw++, offset += n) {
			if (ice06_a[xw[iw]] != NULL) {
				xd[i + offset] = ice06_a[xw[iw]](&gamma);
			} else {
				xd[i + offset] = NA_REAL;
			}
		}
	);

	UNPROTECT(1);
	return d;
}

SEXP R_nroot_control(SEXP trace, SEXP maxit, SEXP abstol, SEXP reltol)
{
	nroot_default.trace = asInteger(trace);
	nroot_default.maxit = asInteger(maxit);
	nroot_default.abstol = asReal(abstol);
	nroot_default.reltol = asReal(reltol);
	return R_NilValue;
}

#define ADDENTRY(f, n) {"C_" # f, (DL_FUNC) &R_ ## f, n}

static const R_CallMethodDef CallEntries[] = {
	ADDENTRY(if97_pt,	4),
	ADDENTRY(if97_ph,	3),
	ADDENTRY(if97_tsat,	1),
	ADDENTRY(if97_psat,	1),
	ADDENTRY(if97_state_pt,	2),
#if REGISTER_IF97_REGION
	ADDENTRY(if97_region_pt,	2),
	ADDENTRY(if97_region_ph,	2),
#endif
	ADDENTRY(iapws95,	3),
	ADDENTRY(iapws95_rhot,	4),
	ADDENTRY(iapws95_pt,	4),
	ADDENTRY(iapws95_ph,	3),
	ADDENTRY(iapws95_sat_t,	2),
	ADDENTRY(iapws95_sat_p,	2),
	ADDENTRY(iapws95_state_pt,	2),
	ADDENTRY(iapws95_state_rhot,	2),
	ADDENTRY(ice06_pt,	3),
	ADDENTRY(nroot_control,	4),
	{ NULL, NULL, 0 },
};

void R_init_iapws(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}
