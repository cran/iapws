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

#ifndef R_INTERFACE_H
#define R_INTERFACE_H

#include <Rinternals.h>
#include <R_ext/Itermacros.h>

#define NINTERRUPT	1000000U

#define RDAT_(x)	RDAT_##x
#define RDAT_int	INTEGER
#define RDAT_double	REAL

#define RSXP_(x)	RSXP_##x
#define RSXP_int	INTSXP
#define RSXP_double	REALSXP

#define MAKE_R_FUN_1(ta, FUN, t1)  \
	SEXP R_ ## FUN (SEXP s1) {  \
		R_xlen_t i, n = xlength(s1);  \
		SEXP ans = PROTECT(allocVector(RSXP_(ta), n));  \
		t1 *d1 = RDAT_(t1)(s1);  \
		ta *da = RDAT_(ta)(ans);  \
		R_ITERATE_CHECK(NINTERRUPT, n, i, da[i] = FUN(d1[i]););  \
		UNPROTECT(1);  \
		return ans;  \
	}

#define MAKE_R_FUN_2(ta, FUN, t1, t2)  \
	SEXP R_ ## FUN (SEXP s1, SEXP s2) {  \
		R_xlen_t i1, n1 = xlength(s1);  \
		R_xlen_t i2, n2 = xlength(s2);  \
		R_xlen_t i, n = (n1 > n2 ? n1 : n2);  \
		SEXP ans = PROTECT(allocVector(RSXP_(ta), n));  \
		t1 *d1 = RDAT_(t1)(s1);  \
		t2 *d2 = RDAT_(t2)(s2);  \
		ta *da = RDAT_(ta)(ans);  \
		MOD_ITERATE2_CHECK(NINTERRUPT, n, n1, n2, i, i1, i2,  \
				da[i] = FUN(d1[i1], d2[i2]););  \
		UNPROTECT(1);  \
		return ans;  \
	}

#define MAKE_R_FUN_3(ta, FUN, t1, t2, t3)  \
	SEXP R_ ## FUN (SEXP s1, SEXP s2, SEXP s3) {  \
		R_xlen_t i1, n1 = xlength(s1);  \
		R_xlen_t i2, n2 = xlength(s2);  \
		R_xlen_t i3, n3 = xlength(s3);  \
		R_xlen_t i, n = (n1 > n2 ? n1 : n2);  \
		if (n < n3) n = n3;  \
		SEXP ans = PROTECT(allocVector(RSXP_(ta), n));  \
		t1 *d1 = RDAT_(t1)(s1);  \
		t2 *d2 = RDAT_(t2)(s2);  \
		t3 *d3 = RDAT_(t3)(s3);  \
		ta *da = RDAT_(ta)(ans);  \
		MOD_ITERATE3_CHECK(NINTERRUPT, n, n1, n2, n3, i, i1, i2, i3, \
				da[i] = FUN(d1[i1], d2[i2], d3[i3]););  \
		UNPROTECT(1);  \
		return ans;  \
	}

#define MAKE_R_PHI_3(FUN, PHI, ARR)  \
	SEXP R_ ## FUN (SEXP s1, SEXP s2, SEXP s3) {  \
		R_xlen_t i1, n1 = xlength(s1);  \
		R_xlen_t i2, n2 = xlength(s2);  \
		R_xlen_t i3, n3 = xlength(s3);  \
		R_xlen_t i, n = (n2 > n3 ? n2 : n3);  \
		SEXP ans = PROTECT(allocMatrix(REALSXP, n, n1));  \
		int *d1 = INTEGER(s1);  \
		double *d2 = REAL(s2);  \
		double *d3 = REAL(s3);  \
		double *da = REAL(ans);  \
		struct iapws_phi phi;  \
		int offset;  \
		MOD_ITERATE2_CHECK(NINTERRUPT, n, n2, n3, i, i2, i3,  \
			PHI(d2[i2], d3[i3], &phi);  \
			for (i1 = 0, offset = 0; i1 < n1; i1++, offset += n) {  \
				da[i + offset] = ARR[d1[i1]](&phi);  \
			}  \
		);  \
		UNPROTECT(1);  \
		return ans;  \
	}

#define MAKE_R_PHI_CHECK_2(FUN, PHI, ARR)  \
	SEXP R_ ## FUN (SEXP s1, SEXP s2) {  \
		int i1, n1 = xlength(s1);  \
		int i2, n2 = xlength(s2);  \
		int n12 = n1 * n2;  \
		SEXP ans = PROTECT(alloc3DArray(REALSXP, n2, n1, 2));  \
		int *d1 = INTEGER(s1);  \
		double *d2 = REAL(s2);  \
		double *da = REAL(ans);  \
		struct iapws_phi phi1, phi2;  \
		int offset;  \
		R_ITERATE_CHECK(NINTERRUPT, n2, i2,  \
			if (PHI(d2[i2], &phi1, &phi2) == 0) {  \
			for (i1 = 0, offset = 0; i1 < n1; i1++, offset += n2) {  \
				da[i2 + offset] = ARR[d1[i1]](&phi1);  \
				da[i2 + offset + n12] = ARR[d1[i1]](&phi2);  \
			}  \
			} else {  \
			for (i1 = 0, offset = 0; i1 < n1; i1++, offset += n2) {  \
				da[i2 + offset] = NA_REAL;  \
				da[i2 + offset + n12] = NA_REAL;  \
			}  \
			}  \
		);  \
		UNPROTECT(1);  \
		return ans;  \
	}

#define MAKE_R_PHI_CHECK_3(FUN, PHI, ARR)  \
	SEXP R_ ## FUN (SEXP s1, SEXP s2, SEXP s3) {  \
		R_xlen_t i1, n1 = xlength(s1);  \
		R_xlen_t i2, n2 = xlength(s2);  \
		R_xlen_t i3, n3 = xlength(s3);  \
		R_xlen_t i, n = (n2 > n3 ? n2 : n3);  \
		SEXP ans = PROTECT(allocMatrix(REALSXP, n, n1));  \
		int *d1 = INTEGER(s1);  \
		double *d2 = REAL(s2);  \
		double *d3 = REAL(s3);  \
		double *da = REAL(ans);  \
		struct iapws_phi phi;  \
		int offset;  \
		MOD_ITERATE2_CHECK(NINTERRUPT, n, n2, n3, i, i2, i3,  \
			if (PHI(d2[i2], d3[i3], &phi) == 0) {  \
			for (i1 = 0, offset = 0; i1 < n1; i1++, offset += n) {  \
				da[i + offset] = ARR[d1[i1]](&phi);  \
			}  \
			} else {  \
			for (i1 = 0, offset = 0; i1 < n1; i1++, offset += n) {  \
				da[i + offset] = NA_REAL;  \
			}  \
			}  \
		);  \
		UNPROTECT(1);  \
		return ans;  \
	}

#define MAKE_R_PHI_CHECK_4(FUN, PHI, ARR)  \
	SEXP R_ ## FUN (SEXP s1, SEXP s2, SEXP s3, SEXP s4) {  \
		R_xlen_t i1, n1 = xlength(s1);  \
		R_xlen_t i2, n2 = xlength(s2);  \
		R_xlen_t i3, n3 = xlength(s3);  \
		R_xlen_t i4, n4 = xlength(s4);  \
		R_xlen_t i, n = (n2 > n3 ? n2 : n3);  \
		if (n < n4) n = n4;  \
		SEXP ans = PROTECT(allocMatrix(REALSXP, n, n1));  \
		int *d1 = INTEGER(s1);  \
		double *d2 = REAL(s2);  \
		double *d3 = REAL(s3);  \
		int *d4 = INTEGER(s4);  \
		double *da = REAL(ans);  \
		struct iapws_phi phi;  \
		int offset;  \
		MOD_ITERATE3_CHECK(NINTERRUPT, n, n2, n3, n4, i, i2, i3, i4,  \
			if (PHI(d2[i2], d3[i3], d4[i4], &phi) == 0) {  \
			for (i1 = 0, offset = 0; i1 < n1; i1++, offset += n) {  \
				da[i + offset] = ARR[d1[i1]](&phi);  \
			}  \
			} else {  \
			for (i1 = 0, offset = 0; i1 < n1; i1++, offset += n) {  \
				da[i + offset] = NA_REAL;  \
			}  \
			}  \
		);  \
		UNPROTECT(1);  \
		return ans;  \
	}

#endif
