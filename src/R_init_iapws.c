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

#include <Rinternals.h>

#define ADDENTRY(f, n) {#f, (DL_FUNC) &f, n}

SEXP R_if97_region(SEXP, SEXP);
SEXP R_if97_state(SEXP, SEXP);
SEXP R_if97_tsat(SEXP);
SEXP R_if97_psat(SEXP);
SEXP R_if97(SEXP, SEXP, SEXP, SEXP);

SEXP R_iapws95_state(SEXP, SEXP);
SEXP R_iapws95_sat(SEXP, SEXP);
SEXP R_iapws95_sat_p(SEXP, SEXP);
SEXP R_iapws95(SEXP, SEXP, SEXP);
SEXP R_iapws95_pt(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
	ADDENTRY(R_if97_region,	2),
	ADDENTRY(R_if97_state,	2),
	ADDENTRY(R_if97_tsat,	1),
	ADDENTRY(R_if97_psat,	1),
	ADDENTRY(R_if97,	4),
	ADDENTRY(R_iapws95_state,	2),
	ADDENTRY(R_iapws95_sat,	2),
	ADDENTRY(R_iapws95_sat_p,	2),
	ADDENTRY(R_iapws95,	3),
	ADDENTRY(R_iapws95_pt,	4),
	{NULL, NULL, 0},
};

void R_init_iapws(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
}
