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

#ifndef IAPWS_IAPWS_H
#define IAPWS_IAPWS_H

/*
 * Constants
 */

/* Molar mass */
#define IAPWS_M		18.015268	/* g/mol */

/* Critical point */
#define IAPWS_TC	647.096		/* K     */
#define IAPWS_PC	22.064		/* MPa   */
#define IAPWS_RHOC	322.0		/* kg/m3 */

/* Triple point */
#define IAPWS_TT	273.16		/* K     */
#define IAPWS_PT	611.657e-6	/* MPa   */

/* Normal condition */
#define IAPWS_PN	0.101325	/* MPa   */

enum iapws_state {
	IAPWS_UNDEF = -1,
	IAPWS_SOLID,
	IAPWS_LIQUID,
	IAPWS_GAS,
	IAPWS_CRIT,
	IAPWS_SAT,
};

enum iapws_ice {
	ICE_IH  = 1,
	ICE_III = 3,
	ICE_V   = 5,
	ICE_VI  = 6,
	ICE_VII = 7,
};

struct iapws_phi {
	enum { IAPWS_PHI, IAPWS_GAMMA } type;
	double d00, d10, d01, d11, d20, d02;
	double p, rho, t, h;
	double R;
};

double iapws_rho(const struct iapws_phi *phi);
double iapws_t(const struct iapws_phi *phi);
double iapws_p(const struct iapws_phi *phi);
double iapws_v(const struct iapws_phi *phi);
double iapws_f(const struct iapws_phi *phi);
double iapws_g(const struct iapws_phi *phi);
double iapws_u(const struct iapws_phi *phi);
double iapws_h(const struct iapws_phi *phi);
double iapws_s(const struct iapws_phi *phi);
double iapws_cv(const struct iapws_phi *phi);
double iapws_cp(const struct iapws_phi *phi);
double iapws_w(const struct iapws_phi *phi);
double iapws_alpha(const struct iapws_phi *phi);
double iapws_beta(const struct iapws_phi *phi);
double iapws_kappat(const struct iapws_phi *phi);

struct iapws_phi_call {
	void (*iapws_phi)(double rho, double t, struct iapws_phi *phi);
	struct iapws_phi *phi;
};

void get_phi_pt(double *x, void *xcall, double *fx, double *dfx);
void get_phi_ph(double *x, void *xcall, double *fx, double *dfx);
void get_sat_t(double *x, void *xcall, double *fx, double *dfx);
void get_sat_p(double *x, void *xcall, double *fx, double *dfx);
void get_gamma_ph(double *x, void *xcall, double *fx, double *dfx);

#endif
