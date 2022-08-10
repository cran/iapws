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

#define RDAT_(x)	RDAT_##x
#define RDAT_int	INTEGER
#define RDAT_double	REAL

#define RSXP_(x)	RSXP_##x
#define RSXP_int	INTSXP
#define RSXP_double	REALSXP

#define MAKE_R_FUN_1(FUN, tx1, ty)  \
	SEXP R_ ## FUN (SEXP x1) {  \
		R_xlen_t i, n = xlength(x1);  \
		SEXP y = PROTECT(allocVector(RSXP_(ty), n));  \
		tx1 *xx1 = RDAT_(tx1)(x1);  \
		ty *xy = RDAT_(ty)(y);  \
		R_ITERATE(n, i, xy[i] = FUN(xx1[i]););  \
		UNPROTECT(1);  \
		return y;  \
	}

#define MAKE_R_FUN_2(FUN, tx1, tx2, ty)  \
	SEXP R_ ## FUN (SEXP x1, SEXP x2) {  \
		R_xlen_t i1, n1 = xlength(x1);  \
		R_xlen_t i2, n2 = xlength(x2);  \
		R_xlen_t i, n = (n1 > n2 ? n1 : n2);  \
		SEXP y = PROTECT(allocVector(RSXP_(ty), n));  \
		tx1 *xx1 = RDAT_(tx1)(x1);  \
		tx2 *xx2 = RDAT_(tx2)(x2);  \
		ty *xy = RDAT_(ty)(y);  \
		MOD_ITERATE2(n, n1, n2, i, i1, i2,  \
				xy[i] = FUN(xx1[i1], xx2[i2]););  \
		UNPROTECT(1);  \
		return y;  \
	}

#define MAKE_R_FUN_3(FUN, tx1, tx2, tx3, ty)  \
	SEXP R_ ## FUN (SEXP x1, SEXP x2, SEXP x3) {  \
		R_xlen_t i1, n1 = xlength(x1);  \
		R_xlen_t i2, n2 = xlength(x2);  \
		R_xlen_t i3, n3 = xlength(x3);  \
		R_xlen_t i, n = (n1 > n2 ? n1 : n2);  \
		if (n < n3) n = n3;  \
		SEXP y = PROTECT(allocVector(RSXP_(ty), n));  \
		tx1 *xx1 = RDAT_(tx1)(x1);  \
		tx2 *xx2 = RDAT_(tx2)(x2);  \
		tx3 *xx3 = RDAT_(tx3)(x3);  \
		ty *xy = RDAT_(ty)(y);  \
		MOD_ITERATE3(n, n1, n2, n3, i, i1, i2, i3, \
				xy[i] = FUN(xx1[i1], xx2[i2], xx3[i3]););  \
		UNPROTECT(1);  \
		return y;  \
	}

#endif
