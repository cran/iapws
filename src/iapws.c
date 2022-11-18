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

#include <math.h>

#include "iapws.h"

double iapws_t(const struct iapws_phi *phi)  /* K */
{
	return phi->t;
}

double iapws_rho(const struct iapws_phi *phi)  /* kg/m3 */
{
	switch (phi->type) {
		case IAPWS_PHI:
			return phi->rho;
		case IAPWS_GAMMA:
			return 1.0 / iapws_v(phi);
		default:
			return 0.0;
	}
}

double iapws_p(const struct iapws_phi *phi)  /* MPa */
{
	switch (phi->type) {
		case IAPWS_PHI:
			return phi->d10 * phi->rho * phi->t * phi->R * 1e-3;
		case IAPWS_GAMMA:
			return phi->p;
		default:
			return 0.0;
	}
}

double iapws_v(const struct iapws_phi *phi)  /* m3/kg */
{
	switch (phi->type) {
		case IAPWS_PHI:
			return 1.0 / phi->rho;
		case IAPWS_GAMMA:
			return phi->d10 / phi->p * phi->t * phi->R * 1e-3;
		default:
			return 0.0;
	}
}

double iapws_f(const struct iapws_phi *phi)  /* kJ/kg */
{
	double ans;
	switch (phi->type) {
		case IAPWS_PHI:
			ans = phi->d00;
			break;
		case IAPWS_GAMMA:
			ans = phi->d00 - phi->d10;
			break;
		default:
			return 0.0;
	}
	return ans * phi->t * phi->R;
}

double iapws_g(const struct iapws_phi *phi)  /* kJ/kg */
{
	double ans;
	switch (phi->type) {
		case IAPWS_PHI:
			ans = phi->d00 + phi->d10;
			break;
		case IAPWS_GAMMA:
			ans = phi->d00;
			break;
		default:
			return 0.0;
	}
	return ans * phi->t * phi->R;
}

double iapws_u(const struct iapws_phi *phi)  /* kJ/kg */
{
	double ans;
	switch (phi->type) {
		case IAPWS_PHI:
			ans = phi->d01;
			break;
		case IAPWS_GAMMA:
			ans = phi->d01 - phi->d10;
			break;
		default:
			return 0.0;
	}
	return ans * phi->t * phi->R;
}

double iapws_h(const struct iapws_phi *phi)  /* kJ/kg */
{
	double ans;
	switch (phi->type) {
		case IAPWS_PHI:
			ans = phi->d01 + phi->d10;
			break;
		case IAPWS_GAMMA:
			ans = phi->d01;
			break;
		default:
			return 0.0;
	}
	return ans * phi->t * phi->R;
}

double iapws_s(const struct iapws_phi *phi)  /* kJ/kg/K */
{
	return (phi->d01 - phi->d00) * phi->R;
}

double iapws_cv(const struct iapws_phi *phi)  /* kJ/kg/K */
{
	double ans;
	switch (phi->type) {
		case IAPWS_PHI:
			ans = -phi->d02;
			break;
		case IAPWS_GAMMA:
			ans = -phi->d02 + (phi->d10 - phi->d11) *
				(phi->d10 - phi->d11) / phi->d20;
			break;
		default:
			return 0.0;
	}
	return ans * phi->R;
}

double iapws_cp(const struct iapws_phi *phi)  /* kJ/kg/K */
{
	double ans;
	switch (phi->type) {
		case IAPWS_PHI:
			ans = -phi->d02 + (phi->d10 - phi->d11) *
				(phi->d10 - phi->d11) /
				(phi->d10 * 2.0 + phi->d20);
			break;
		case IAPWS_GAMMA:
			ans = -phi->d02;
			break;
		default:
			return 0.0;
	}
	return ans * phi->R;
}

double iapws_w(const struct iapws_phi *phi)  /* m/s */
{
	double ans;
	switch (phi->type) {
		case IAPWS_PHI:
			ans = phi->d10 * 2.0 + phi->d20 -
				(phi->d10 - phi->d11) *
				(phi->d10 - phi->d11) / phi->d02;
			break;
		case IAPWS_GAMMA:
			ans = phi->d10 * phi->d10 /
				((phi->d10 - phi->d11) *
				 (phi->d10 - phi->d11) /
				 phi->d02 - phi->d20);
			break;
		default:
			return 0.0;
	}
	return sqrt(ans * phi->R * phi->t * 1e3);
}

/* alpha = beta * kappat */
/* -1/v dv/dt */
double iapws_alpha(const struct iapws_phi *phi)  /* 1/K */
{
	double ans = (phi->d10 - phi->d11) / phi->t;
	switch (phi->type) {
		case IAPWS_PHI:
			return ans / (phi->d20 + phi->d10 * 2.0);
		case IAPWS_GAMMA:
			return ans / phi->d10;
		default:
			return 0.0;
	}
}

/* dp/dt */
double iapws_beta(const struct iapws_phi *phi)  /* MPa/K */
{
	double ans = phi->d10 - phi->d11;
	switch (phi->type) {
		case IAPWS_PHI:
			return ans * phi->R * phi->rho * 1e-3;
		case IAPWS_GAMMA:
			return -ans * phi->p / (phi->d20 * phi->t);
		default:
			return 0.0;
	}
}

/* -1/v dv/dp */
double iapws_kappat(const struct iapws_phi *phi)  /* 1/MPa */
{
	switch (phi->type) {
		case IAPWS_PHI:
			return 1e3 / ((phi->d10 * 2.0 + phi->d20) *
					phi->rho * phi->R * phi->t);
		case IAPWS_GAMMA:
			return -phi->d20 / (phi->p * phi->d10);
		default:
			return 0.0;
	}
}

/*
 * wrappers for nroot
 */

void get_phi_pt(double *x, void *xcall, double *fx, double *dfx)
{
	double rho = *x;

	struct iapws_phi_call *call = (struct iapws_phi_call *)(xcall);
	struct iapws_phi *phi = call->phi;

	call->iapws_phi(rho, phi->t, phi);

	double fact = phi->R * 1e-3 * phi->t;
	*fx = phi->d10 * rho * fact - phi->p;
	*dfx = (phi->d10 * 2.0 + phi->d20) * fact;
}

void get_phi_ph(double *x, void *xcall, double *fx, double *dfx)
{
	double rho = x[0];
	double t = x[1];

	struct iapws_phi_call *call = (struct iapws_phi_call *)(xcall);
	struct iapws_phi *phi = call->phi;

	call->iapws_phi(rho, t, phi);

	fx[0] = phi->d10 * rho * phi->R * t * 1e-3 - phi->p;
	fx[1] = (phi->d10 + phi->d01) * phi->R * t - phi->h;

	dfx[0] = (phi->d10 * 2.0 + phi->d20) * phi->R * t * 1e-3;
	dfx[1] = (phi->d10 + phi->d20 + phi->d11) / rho * phi->R * t;
	dfx[2] = (phi->d10 - phi->d11) * rho * phi->R * 1e-3;
	dfx[3] = (phi->d10 - phi->d11 - phi->d02) * phi->R;
}

void get_sat_t(double *x, void *xcall, double *fx, double *dfx)
{
	double rhol = x[0];
	double rhog = x[1];

	struct iapws_phi_call *call = (struct iapws_phi_call *)(xcall);
	struct iapws_phi *phil = call[0].phi;
	struct iapws_phi *phig = call[1].phi;

	call[0].iapws_phi(rhol, phil->t, phil);
	call[1].iapws_phi(rhog, phig->t, phig);

	fx[0] = phil->d10 * rhol - phig->d10 * rhog;
	fx[1] = phil->d00 + phil->d10 - phig->d00 - phig->d10;

	dfx[0] = (phil->d10 * 2.0 + phil->d20);
	dfx[1] = dfx[0] / rhol;
	dfx[2] = -(phig->d10 * 2.0 + phig->d20);
	dfx[3] = dfx[2] / rhog;
}

void get_sat_p(double *x, void *xcall, double *fx, double *dfx)
{
	double rhol = x[0];
	double rhog = x[1];
	double t = x[2];

	struct iapws_phi_call *call = (struct iapws_phi_call *)(xcall);
	struct iapws_phi *phil = call[0].phi;
	struct iapws_phi *phig = call[1].phi;

	call[0].iapws_phi(rhol, t, phil);
	call[1].iapws_phi(rhog, t, phig);

	fx[0] = phil->d10 * rhol - phig->d10 * rhog;
	fx[1] = phil->d00 + phil->d10 - phig->d00 - phig->d10;
	fx[2] = phil->d10 * rhol * phil->R * t * 1e-3 - phil->p;

	dfx[0] = (phil->d10 * 2.0 + phil->d20);
	dfx[1] = dfx[0] / rhol;
	dfx[2] = dfx[0] * phil->R * t * 1e-3;
	dfx[3] = -(phig->d10 * 2.0 + phig->d20);
	dfx[4] = dfx[3] / rhog;
	dfx[5] = 0.0;
	dfx[6] = (-phil->d11 * rhol + phig->d11 * rhog) / t;
	dfx[7] = (-phil->d01 - phil->d11 + phig->d01 + phig->d11) / t;
	dfx[8] = (phil->d10 - phil->d11) * rhol * phil->R * 1e-3;
}

void get_gamma_ph(double *x, void *xcall, double *fx, double *dfx)
{
	double t = *x;

	struct iapws_phi_call *call = (struct iapws_phi_call *)(xcall);
	struct iapws_phi *gamma = call->phi;

	call->iapws_phi(gamma->p, t, gamma);

	*fx = gamma->d01 * gamma->R * gamma->t - gamma->h;
	*dfx = -gamma->d02 * gamma->R;
}

