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
#include "heavy17.h"
#include "trans.h"
#include "surf.h"
#include "elec.h"
#include "nroot.h"

#define REGISTER_EXTRA 0
#if REGISTER_EXTRA
#	include "sat86.h"
#endif

/* dummy function when property is not available */
static double iapws_na(const struct iapws_phi *phi)
{
	return NA_REAL;
}

/* IF97 */

static double (*if97_a[])(const struct iapws_phi *phi) = {
	iapws_f, iapws_g,
	iapws_u, iapws_h,
	iapws_s, iapws_t,
	iapws_p, iapws_v,
	iapws_cp, iapws_cv,
	iapws_w, iapws_rho,
	iapws_alpha, iapws_beta,
	iapws_kappat,
	if97_eta, if97_lambda,
};

MAKE_R_PHI_CHECK_4(if97_pt, if97_gamma_pt, if97_a)
MAKE_R_PHI_CHECK_3(if97_ph, if97_gamma_ph, if97_a)

MAKE_R_FUN_1(double, if97_tsat, double)
MAKE_R_FUN_1(double, if97_psat, double)

MAKE_R_FUN_2(int, if97_state_pt, double, double)

#if REGISTER_EXTRA
MAKE_R_FUN_2(int, if97_region_pt, double, double)
MAKE_R_FUN_2(int, if97_region_ph, double, double)
#endif

/* IAPWS95 */

static double (*iapws95_a[])(const struct iapws_phi *phi) = {
	iapws_f, iapws_g,
	iapws_u, iapws_h,
	iapws_s, iapws_t,
	iapws_p, iapws_v,
	iapws_cp, iapws_cv,
	iapws_w, iapws_rho,
	iapws_alpha, iapws_beta,
	iapws_kappat,
	iapws95_eta, iapws95_lambda,
};

MAKE_R_PHI_3(iapws95, iapws95_phi, iapws95_a)

MAKE_R_PHI_CHECK_4(iapws95_rhot, iapws95_phi_rhot, iapws95_a)
MAKE_R_PHI_CHECK_4(iapws95_pt, iapws95_phi_pt, iapws95_a)
MAKE_R_PHI_CHECK_3(iapws95_ph, iapws95_phi_ph, iapws95_a)

MAKE_R_PHI_CHECK_2(iapws95_sat_t, iapws95_sat_t, iapws95_a)
MAKE_R_PHI_CHECK_2(iapws95_sat_p, iapws95_sat_p, iapws95_a)

MAKE_R_FUN_2(int, iapws95_state_pt, double, double)
MAKE_R_FUN_2(int, iapws95_state_rhot, double, double)

/* SAT-86 */

#if REGISTER_EXTRA
MAKE_R_FUN_1(double, sat86_p, double)
MAKE_R_FUN_1(double, sat86_t, double)
MAKE_R_FUN_1(double, sat86_rhol, double)
MAKE_R_FUN_1(double, sat86_rhog, double)
#endif

/* ELEC */

MAKE_R_FUN_2(double, iapws_epsilon, double, double)
MAKE_R_FUN_3(double, iapws_n, double, double, double)
MAKE_R_FUN_2(double, iapws_pk, double, double)
MAKE_R_FUN_1(double, iapws_sigma, double)
MAKE_R_FUN_1(double, heavy17_sigma, double)

/* ICE-06 */

static double (*ice06_a[])(const struct iapws_phi *phi) = {
	iapws_f, iapws_g,
	iapws_u, iapws_h,
	iapws_s, iapws_t,
	iapws_p, iapws_v,
	iapws_cp, iapws_cv,
	iapws_w, iapws_rho,
	iapws_alpha, iapws_beta,
	iapws_kappat,
	iapws_na, iapws_na,
};

MAKE_R_PHI_3(ice06_pt, ice06_gamma, ice06_a)

/* HEAVY-17 */

static double (*heavy17_a[])(const struct iapws_phi *phi) = {
	iapws_f, iapws_g,
	iapws_u, iapws_h,
	iapws_s, iapws_t,
	iapws_p, iapws_v,
	iapws_cp, iapws_cv,
	iapws_w, iapws_rho,
	iapws_alpha, iapws_beta,
	iapws_kappat,
	heavy17_eta, heavy17_lambda,
};

MAKE_R_PHI_3(heavy17, heavy17_phi, heavy17_a)

MAKE_R_PHI_CHECK_4(heavy17_rhot, heavy17_phi_rhot, heavy17_a)
MAKE_R_PHI_CHECK_4(heavy17_pt, heavy17_phi_pt, heavy17_a)

MAKE_R_PHI_CHECK_2(heavy17_sat_t, heavy17_sat_t, heavy17_a)
MAKE_R_PHI_CHECK_2(heavy17_sat_p, heavy17_sat_p, heavy17_a)

MAKE_R_FUN_2(int, heavy17_state_pt, double, double)
MAKE_R_FUN_2(int, heavy17_state_rhot, double, double)

#if REGISTER_EXTRA
MAKE_R_FUN_1(double, heavy17_psat, double)
MAKE_R_FUN_1(double, heavy17_tsat, double)
MAKE_R_FUN_1(double, heavy17_rhol, double)
MAKE_R_FUN_1(double, heavy17_rhog, double)
#endif

/* NROOT */

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
	ADDENTRY(if97_pt,		4),
	ADDENTRY(if97_ph,		3),
	ADDENTRY(if97_tsat,		1),
	ADDENTRY(if97_psat,		1),
	ADDENTRY(if97_state_pt,		2),
	ADDENTRY(iapws95,		3),
	ADDENTRY(iapws95_rhot,		4),
	ADDENTRY(iapws95_pt,		4),
	ADDENTRY(iapws95_ph,		3),
	ADDENTRY(iapws95_sat_t,		2),
	ADDENTRY(iapws95_sat_p,		2),
	ADDENTRY(iapws95_state_pt,	2),
	ADDENTRY(iapws95_state_rhot,	2),
	ADDENTRY(ice06_pt,		3),
	ADDENTRY(heavy17,		3),
	ADDENTRY(heavy17_rhot,		4),
	ADDENTRY(heavy17_pt,		4),
	ADDENTRY(heavy17_sat_t,		2),
	ADDENTRY(heavy17_sat_p,		2),
	ADDENTRY(heavy17_state_pt,	2),
	ADDENTRY(heavy17_state_rhot,	2),
	ADDENTRY(iapws_epsilon,		2),
	ADDENTRY(iapws_n,		3),
	ADDENTRY(iapws_pk,		2),
	ADDENTRY(iapws_sigma,		1),
	ADDENTRY(heavy17_sigma,		1),
	ADDENTRY(nroot_control,		4),
#if REGISTER_EXTRA
	ADDENTRY(if97_region_pt,	2),
	ADDENTRY(if97_region_ph,	2),
	ADDENTRY(sat86_p,		1),
	ADDENTRY(sat86_t,		1),
	ADDENTRY(sat86_rhol,		1),
	ADDENTRY(sat86_rhog,		1),
	ADDENTRY(heavy17_psat,		1),
	ADDENTRY(heavy17_tsat,		1),
	ADDENTRY(heavy17_rhol,		1),
	ADDENTRY(heavy17_rhog,		1),
#endif
	{ NULL, NULL, 0 },
};

void R_init_iapws(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}
