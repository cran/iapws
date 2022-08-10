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

#ifndef IAPWS_H
#define IAPWS_H

#include <Rmath.h>
#define pow_di R_pow_di
#define POW2(x) ((x) * (x))

#define IAPWS_TC	647.096		/* K	   */
#define IAPWS_RHOC	322.0		/* kg/m3   */
#define IAPWS_PC	22.064		/* MPa	   */

extern int iapws_nroot_max_iter;
extern double iapws_nroot_tol;

typedef enum {
	IAPWS_UNDEF = -1,
	IAPWS_SOLID,
	IAPWS_LIQUID,
	IAPWS_GAS,
	IAPWS_CRIT,
} iapws_state_id;

typedef struct {
	enum { IAPWS_PHI, IAPWS_GAMMA } type;
	double d00;
	double d10;
	double d01;
	double d11;
	double d20;
	double d02;
	double p, rho, t;
	double R;
} iapws_phi;

double iapws_rho(iapws_phi *phi);
double iapws_t(iapws_phi *phi);
double iapws_p(iapws_phi *phi);
double iapws_v(iapws_phi *phi);
double iapws_f(iapws_phi *phi);
double iapws_g(iapws_phi *phi);
double iapws_u(iapws_phi *phi);
double iapws_h(iapws_phi *phi);
double iapws_s(iapws_phi *phi);
double iapws_cv(iapws_phi *phi);
double iapws_cp(iapws_phi *phi);
double iapws_w(iapws_phi *phi);
double iapws_alpha(iapws_phi *phi);
double iapws_beta(iapws_phi *phi);
double iapws_chit(iapws_phi *phi);
double iapws_sigma(iapws_phi *phi);

#endif
