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

/*
 * International Association for the Properties of Water and Steam,
 * IAPWS R16-17(2018), Revised Release on the IAPWS Formulation 2017 for the
 * Thermodynamic Properties of Heavy Water (2018)
 */

#include <assert.h>
#include <math.h>

#include "heavy17.h"
#include "nroot.h"
#include "phi.h"
#include "pow.h"

void heavy17_phi(double rho, double t, struct iapws_phi *phi)
{
	const struct iapws_phi_coef0 coef0[7] = {
		{ -8.67099402264600,	0.0 },
		{  6.96033578458778,	0.0 },
		{  3.0,			0.0 },
		{  0.10633e-1,		308.0	/ HEAVY17_TC },
		{  0.99787,		1695.0	/ HEAVY17_TC },
		{  0.21483e1,		3949.0	/ HEAVY17_TC },
		{  0.3549,		10317.0	/ HEAVY17_TC },
	};
	const struct iapws_phi_coef1 coef1[6 + 6] = {
		{ 0,	1,	0.6555,	 0.29695687e1	* M_E },
		{ 0,	1,	0.9369,	-0.37900454e1	* M_E },
		{ 0,	2,	0.5610,	 0.94108960	* M_E },
		{ 0,	2,	0.7017,	-0.92246625	* M_E },
		{ 0,	3,	1.0672,	-0.13960419e-1	* M_E },
		{ 0,	4,	1.0000,	 0.12208206e-1	* M_E },
		{ 1,	1,	3.9515,	-0.12520357	},
		{ 1,	2,	0.2000,	-0.35947024e-1	},
		{ 2,	1,	4.6000,	-0.55539150e1	},
		{ 2,	1,	2.3660,	-0.69183515	},
		{ 2,	2,	5.4644,	-0.93617287e1	},
		{ 2,	3,	5.1590,	-0.49300974e1	},
	};
	const struct iapws_phi_coef2 coef2[12] = {
		{ 1,	1.5745,	 0.86000607e1,	1.5305,	1.2888,	1.7385,	0.5803 },
		{ 1,	1.9000,	 0.42109740,	18.738,	1177.0,	1.0491,	0.9488 },
		{ 1,	2.6200,	 0.51965200e1,	1.2148,	0.9465,	2.5617,	0.9683 },
		{ 1,	3.4553,	-0.45611060e-1,	0.6014,	0.4200,	1.5414,	1.8663 },
		{ 1,	3.8106,	 0.16447690e2,	1.3086,	0.3673,	2.7242,	0.6815 },
		{ 1,	4.3200,	-0.39192110,	18.677,	1167.0,	1.0486,	0.9487 },
		{ 1,	4.8950,	 0.27039336e1,	1.3528,	0.9504,	3.5321,	0.9495 },
		{ 2,	1.4300,	 0.37563747e2,	3.4456,	7.8318,	2.4552,	1.1158 },
		{ 2,	1.5870,	-0.17760776e1,	1.2645,	3.3281,	0.8319,	0.1607 },
		{ 2,	3.7900,	 0.22092464e1,	2.5547,	7.1753,	1.3500,	0.4144 },
		{ 3,	1.4150,	-0.22451330e1,	1.4723,	2.4318,	1.3794,	0.2895 },
		{ 3,	3.4540,	-0.24841042e1,	2.4297,	8.2710,	1.3045,	0.2236 },
	};

	phi->type = IAPWS_PHI;
	phi->rho = rho;
	phi->t = t;
	phi->R = HEAVY17_R;

	iapws_phi(rho / HEAVY17_RHOC, HEAVY17_TC / t,
			coef0, ARRAY_SIZE(coef0),
			coef1, ARRAY_SIZE(coef1),
			coef2, ARRAY_SIZE(coef2), phi);
}

#define PEPS	1.0001
#define REPS	1.01

int heavy17_phi_rhot(double rho, double t, enum iapws_state state,
		struct iapws_phi *phi)
{
	switch (state) {
		case IAPWS_LIQUID:
		case IAPWS_GAS:
		case IAPWS_CRIT:
			heavy17_phi(rho, t, phi);
			return 0;
		default:
			return -1;
	}
	assert(0);
}

int heavy17_phi_pt(double p, double t, enum iapws_state state,
		struct iapws_phi *phi)
{
	struct nroot_control ctrl = nroot_default;

	/* Find suitable start value for rho */
	double rho;
	if (state == IAPWS_LIQUID) {
		/* Rackett equation */
		rho = (HEAVY17_PC * 1.0e3) / (HEAVY17_TC * HEAVY17_R);
		rho *= POW(HEAVY17_RHOC / rho, 1.0 +
				POW(fabs(1.0 - t / HEAVY17_TC), 2.0 / 7.0));
	} else if (state == IAPWS_GAS) {
		/* Ideal gas */
		rho = (p * 1.0e3) / (t * HEAVY17_R);
	} else if (state == IAPWS_CRIT) {
		/* Good enough default value */
		rho = HEAVY17_RHOC * 2.0;
	} else {
		return -1;
	}

	phi->p = p;
	phi->t = t;
	struct iapws_phi_call call = { heavy17_phi, phi };
	return nroot1(get_phi_pt, &rho, &call, &ctrl) == NROOT_SUCCESS ? 0 : -1;
}

static double sumpow6(const double x, const double coef[][2])
{
	return  POW(x, coef[0][1]) * coef[0][0] +
		POW(x, coef[1][1]) * coef[1][0] +
		POW(x, coef[2][1]) * coef[2][0] +
		POW(x, coef[3][1]) * coef[3][0] +
		POW(x, coef[4][1]) * coef[4][0] +
		POW(x, coef[5][1]) * coef[5][0];
}

static double sat_logpi(double theta)
{
	const double coef[6][2] = {
		{ -0.794440e1, 1.00 },
		{  0.194340e1, 1.50 },
		{ -0.243530e1, 2.44 },
		{ -0.342000e1, 5.30 },
		{  0.355000e2, 14.0 },
		{ -0.302000e3, 20.0 },
	};
	return sumpow6(1.0 - theta, coef) / theta;
}

double heavy17_psat(double t)
{
	if (t < HEAVY17_TT || t > HEAVY17_TC) return 0.0;
	return exp(sat_logpi(t / HEAVY17_TC)) * HEAVY17_PC;
}

static void get_thetasat(double *theta, void *logpi, double *fx, double *dfx)
{
	double *lp = (double *)(logpi);
	*fx = sat_logpi(*theta) - (*lp);
}

double heavy17_tsat(double p)
{
	if (p < HEAVY17_PT || p > HEAVY17_PC) return 0.0;

	struct nroot_control ctrl = nroot_default;

	double logpi = log(p / HEAVY17_PC);
	double theta = 1.0 / (1.0 - 1.401228e-01 * logpi -
			1.207096e-03 * POW2(logpi));
	if (sroot(get_thetasat, &theta, &logpi, &ctrl) != NROOT_SUCCESS)
		return 0.0;
	theta *= HEAVY17_TC;
	if (theta > HEAVY17_TC) theta = HEAVY17_TC;
	else if (theta < HEAVY17_TT) theta = HEAVY17_TT;
	return theta;
}

double heavy17_rhol(double t)
{
	const double coef[6][2] = {
		{  0.166200e1, 0.29 },
		{  0.901130e1, 1.00 },
		{ -0.154210e2, 1.30 },
		{  0.115760e2, 1.77 },
		{ -0.516940e1, 2.50 },
		{ -0.236240e3, 16.0 },
	};
	if (t < HEAVY17_TT || t > HEAVY17_TC) return 0.0;
	return (sumpow6(1.0 - t / HEAVY17_TC, coef) + 1.0) * HEAVY17_RHOC;
}

double heavy17_rhog(double t)
{
	const double coef[6][2] = {
		{ -0.247140e1, 0.33 },
		{ -0.266744e2, 1.29 },
		{  0.531080e2, 1.68 },
		{ -0.480150e2, 2.09 },
		{ -0.576230e2, 6.10 },
		{ -0.371720e3, 17.0 },
	};
	if (t < HEAVY17_TT || t > HEAVY17_TC) return 0.0;
	return exp(sumpow6(1.0 - t / HEAVY17_TC, coef)) * HEAVY17_RHOC;
}

int heavy17_sat_t(double t, struct iapws_phi *phil, struct iapws_phi *phig)
{
	struct nroot_control ctrl = nroot_default;

	double rhol = heavy17_rhol(t);
	if (rhol == 0.0) return -1;  /* test if outside of saturation line */
	double x[2] = { rhol * REPS, heavy17_rhog(t) / REPS};

	phil->t = t;
	phig->t = t;
	struct iapws_phi_call call[2] = {
		{ heavy17_phi, phil },
		{ heavy17_phi, phig },
	};

	return nroot2(get_sat_t, x, call, &ctrl) == NROOT_SUCCESS ? 0 : -1;
}

int heavy17_sat_p(double p, struct iapws_phi *phil, struct iapws_phi *phig)
{
	struct nroot_control ctrl = nroot_default;

	double t = heavy17_tsat(p);
	if (t == 0.0) return -1;  /* test if outside of saturation line */
	double x[3] = { heavy17_rhol(t) * REPS, heavy17_rhog(t) / REPS, t };

	phil->p = p;
	phig->p = p;
	struct iapws_phi_call call[2] = {
		{ heavy17_phi, phil },
		{ heavy17_phi, phig },
	};

	return nrootn(3, get_sat_p, x, call, &ctrl) == NROOT_SUCCESS ? 0 : -1;
}

static const double tmelt[5] = {
	HEAVY17_TT,
	254.415,
	258.661,
	275.748,
	315.0,
};
static const double pmelt[4] = {
	HEAVY17_PT,
	222.41,
	352.19,
	634.53,
};
static double melt_p(double t, enum iapws_ice phase)
{
	switch (phase) {
		case ICE_IH:
			if (t < tmelt[1] || t > tmelt[0]) return 0.0;
			t /= tmelt[0];
			return pmelt[0] *
				(1.0 - 0.30153e5 * (1.0 - POW(t, 5.5)) +
				 0.692503e6 * (1.0 - POW(t, 8.2)));
		case ICE_III:
			if (t < tmelt[1] || t > tmelt[2]) return 0.0;
			t /= tmelt[1];
			return pmelt[1] *
				(1.0 - 0.802871 * (1.0 - POW(t, 33)));
		case ICE_V:
			if (t < tmelt[2] || t > tmelt[3]) return 0.0;
			t /= tmelt[2];
			return pmelt[2] *
				(1.0 - 0.1280388e1 * (1.0 - POW(t, 7.6)));
		case ICE_VI:
			if (t < tmelt[3] || t > tmelt[4]) return 0.0;
			t /= tmelt[3];
			return pmelt[3] *
				(1.0 - 0.1276026e1 * (1.0 - POW(t, 4)));
		default:
			return 0.0;
	}
	assert(0);
}

static double sub_p(double t)
{
	if (t < 210.0 || t > HEAVY17_TT) return 0.0;
	t /= HEAVY17_TT;
	return HEAVY17_PT * exp(-0.1314226e2 * (1.0 - POW(t, -1.73)) +
			0.3212969e2 * (1.0 - POW(t, -1.42)));
}

static enum iapws_state melt_sub_state(double p, double t)
{
	//assert((t < tmelt[0] && p < pmelt[1]) || p > pmelt[1]);

	if (p < sub_p(210.0)) {
		if (t >= 210.0) return IAPWS_GAS;
		else return IAPWS_UNDEF;
	} else if (p < pmelt[0]) {
		if (t < 210.0) return IAPWS_SOLID;
		else if (t > tmelt[0]) return IAPWS_GAS;
		else if (p <= sub_p(t)) return IAPWS_GAS;
		else return IAPWS_SOLID;
	} else if (p < pmelt[1]) {
		if (t < tmelt[1]) return IAPWS_SOLID;
		else if (t > tmelt[0]) return IAPWS_LIQUID;
		else if (p >= melt_p(t, ICE_IH)) return IAPWS_LIQUID;
		else return IAPWS_SOLID;
	} else if (p < pmelt[2]) {
		if (t < tmelt[1]) return IAPWS_SOLID;
		else if (t > tmelt[2]) return IAPWS_LIQUID;
		else if (p <= melt_p(t, ICE_III)) return IAPWS_LIQUID;
		else return IAPWS_SOLID;
	} else if (p < pmelt[3]) {
		if (t < tmelt[2]) return IAPWS_SOLID;
		else if (t > tmelt[3]) return IAPWS_LIQUID;
		else if (p <= melt_p(t, ICE_V)) return IAPWS_LIQUID;
		else return IAPWS_SOLID;
	} else if (p < melt_p(tmelt[4], ICE_VI)) {
		if (t < tmelt[3]) return IAPWS_SOLID;
		else if (t > tmelt[4]) return IAPWS_LIQUID;
		else if (p <= melt_p(t, ICE_VI)) return IAPWS_LIQUID;
		else return IAPWS_SOLID;
	} else {  //if (p > pmelt[4])
		if (t < tmelt[4]) return IAPWS_SOLID;
		else return IAPWS_UNDEF;
	}
}

enum iapws_state heavy17_state_pt(double p, double t)
{
	double ps;
	struct iapws_phi phil, phig;

	if (t >= HEAVY17_TT && t < HEAVY17_TC && p < 640.0) {
		/* Try with fast approximation */
		ps = heavy17_psat(t);
		if (p > ps * PEPS) return IAPWS_LIQUID;
		if (p * PEPS < ps) return IAPWS_GAS;

		heavy17_sat_t(t, &phil, &phig);
		ps = iapws_p(&phig);
		if (p > ps) return IAPWS_LIQUID;
		if (p < ps) return IAPWS_GAS;

		return IAPWS_SAT;
	} else if (t >= HEAVY17_TC) {
		if (p < HEAVY17_PC) return IAPWS_GAS;
		else return IAPWS_CRIT;
	}
	return melt_sub_state(p, t);
}

enum iapws_state heavy17_state_rhot(double rho, double t)
{
	double rhol, rhog;
	struct iapws_phi phil, phig;

	if (t >= HEAVY17_TT && t < HEAVY17_TC) {
		/* Try with fast approximation */
		rhol = heavy17_rhol(t);
		rhog = heavy17_rhog(t);
		if (rho > rhol * REPS) return IAPWS_LIQUID;
		else if (rho * REPS < rhog) return IAPWS_GAS;
		else if (rho * REPS < rhol && rho > rhog * REPS)
			return IAPWS_SAT;

		heavy17_sat_t(t, &phil, &phig);
		rhol = iapws_rho(&phil);
		rhog = iapws_rho(&phig);
		if (rho > rhol) return IAPWS_LIQUID;
		else if (rho < rhog) return IAPWS_GAS;
		else return IAPWS_SAT;
	} else if (t >= HEAVY17_TC) {
		//heavy17_phi(rho, t, &phig);
		//if (iapws_p(&phig) < IAPWS_PC) return IAPWS_GAS;
		//else return IAPWS_CRIT;
		if (rho < HEAVY17_RHOC) return IAPWS_GAS;
		else return IAPWS_CRIT;
	}
	return IAPWS_UNDEF;  /* FIXME */
}

