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

int	iapws_nroot_max_iter = 30;
double	iapws_nroot_tol = 1e-9;

double iapws_t(const iapws_phi *phi)  /* K */
{
	return phi->t;
}

double iapws_rho(const iapws_phi *phi)  /* kg/m3 */
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

double iapws_p(const iapws_phi *phi)  /* MPa */
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

double iapws_v(const iapws_phi *phi)  /* m3/kg */
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

double iapws_f(const iapws_phi *phi)  /* kJ/kg */
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

double iapws_g(const iapws_phi *phi)  /* kJ/kg */
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

double iapws_u(const iapws_phi *phi)  /* kJ/kg */
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

double iapws_h(const iapws_phi *phi)  /* kJ/kg */
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

double iapws_s(const iapws_phi *phi)  /* kJ/kg/K */
{
	return (phi->d01 - phi->d00) * phi->R;
}

double iapws_cv(const iapws_phi *phi)  /* kJ/kg/K */
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

double iapws_cp(const iapws_phi *phi)  /* kJ/kg/K */
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

double iapws_w(const iapws_phi *phi)  /* m/s */
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

/* alpha = beta * chit */
/* -1/v dv/dt */
double iapws_alpha(const iapws_phi *phi)  /* 1/K */
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
double iapws_beta(const iapws_phi *phi)  /* MPa/K */
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
double iapws_chit(const iapws_phi *phi)  /* 1/MPa */
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

double iapws_sigma(const iapws_phi *phi)  /* mN/m */
{
	double t  = 1.0 - phi->t / IAPWS_TC;
	return pow(t, 1.256) * (1.0 - t * 0.625) * 235.8;
}
