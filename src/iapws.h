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

/* Critical point */
#define IAPWS_TC	647.096		/* K     */
#define IAPWS_RHOC	322.0		/* kg/m3 */
#define IAPWS_PC	22.064		/* MPa   */

/* Triple point */
#define IAPWS_TT	273.16		/* K     */
#define IAPWS_PT	611.657e-6	/* MPa   */

#define IAPWS_P0	0.101325	/* MPa   */

typedef enum {
	IAPWS_UNDEF = -1,
	IAPWS_SOLID,
	IAPWS_LIQUID,
	IAPWS_GAS,
	IAPWS_CRIT,
	IAPWS_SAT,
} iapws_state_id;

typedef struct {
	enum { IAPWS_PHI, IAPWS_GAMMA } type;
	double d00;
	double d10;
	double d01;
	double d11;
	double d20;
	double d02;
	double p, rho, t, h;
	double R;
} iapws_phi;

double iapws_rho(const iapws_phi *phi);
double iapws_t(const iapws_phi *phi);
double iapws_p(const iapws_phi *phi);
double iapws_v(const iapws_phi *phi);
double iapws_f(const iapws_phi *phi);
double iapws_g(const iapws_phi *phi);
double iapws_u(const iapws_phi *phi);
double iapws_h(const iapws_phi *phi);
double iapws_s(const iapws_phi *phi);
double iapws_cv(const iapws_phi *phi);
double iapws_cp(const iapws_phi *phi);
double iapws_w(const iapws_phi *phi);
double iapws_alpha(const iapws_phi *phi);
double iapws_beta(const iapws_phi *phi);
double iapws_kappat(const iapws_phi *phi);

#endif
