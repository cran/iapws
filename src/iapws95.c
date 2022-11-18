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
 * IAPWS R6-95(2018), Revised Release on the IAPWS Formulation 1995 for the
 * Thermodynamic Properties of Ordinary Water Substance for General and
 * Scientific Use (2018)
 */

#include <assert.h>
#include <math.h>

#include "iapws95.h"
#include "if97.h"
#include "sat86.h"
#include "melt.h"
#include "nroot.h"
#include "phi.h"
#include "pow.h"

static const struct iapws_phi_coef0 coef0[3 + 5] = {
	{ -8.3204464837497,	0.0 },
	{  6.6832105275932,	0.0 },
	{  3.00632,		0.0 },
	{  0.012436,	 1.28728967 },
	{  0.97315,	 3.53734222 },
	{  1.2795,	 7.74073708 },
	{  0.96956,	 9.24437796 },
	{  0.24873,	 27.5075105 },
};
static const struct iapws_phi_coef1 coef1[7 + 44] = {
	{ 0,	1,	-0.5,	 0.12533547935523e-1 * M_E },	
	{ 0,	1,	0.875,	 0.78957634722828e1  * M_E },	
	{ 0,	1,	1,	-0.87803203303561e1  * M_E },	
	{ 0,	2,	0.5,	 0.31802509345418    * M_E },	
	{ 0,	2,	0.75,	-0.26145533859358    * M_E },	
	{ 0,	3,	0.375,	-0.78199751687981e-2 * M_E },	
	{ 0,	4,	1,	 0.88089493102134e-2 * M_E },	
	{ 1,	1,	4,	-0.66856572307965	},	
	{ 1,	1,	6,	 0.20433810950965	},	
	{ 1,	1,	12,	-0.66212605039687e-4	},	
	{ 1,	2,	1,	-0.19232721156002	},	
	{ 1,	2,	5,	-0.25709043003438	},	
	{ 1,	3,	4,	 0.16074868486251	},	
	{ 1,	4,	2,	-0.40092828925807e-1	},	
	{ 1,	4,	13,	 0.39343422603254e-6	},	
	{ 1,	5,	9,	-0.75941377088144e-5	},	
	{ 1,	7,	3,	 0.56250979351888e-3	},	
	{ 1,	9,	4,	-0.15608652257135e-4	},	
	{ 1,	10,	11,	 0.11537996422951e-8	},	
	{ 1,	11,	4,	 0.36582165144204e-6	},	
	{ 1,	13,	13,	-0.13251180074668e-11	},	
	{ 1,	15,	1,	-0.62639586912454e-9	},	
	{ 2,	1,	7,	-0.10793600908932	},	
	{ 2,	2,	1,	 0.17611491008752e-1	},	
	{ 2,	2,	9,	 0.22132295167546	},	
	{ 2,	2,	10,	-0.40247669763528	},	
	{ 2,	3,	10,	 0.58083399985759	},	
	{ 2,	4,	3,	 0.49969146990806e-2	},	
	{ 2,	4,	7,	-0.31358700712549e-1	},	
	{ 2,	4,	10,	-0.74315929710341	},	
	{ 2,	5,	10,	 0.47807329915480	},	
	{ 2,	6,	6,	 0.20527940895948e-1	},	
	{ 2,	6,	10,	-0.13636435110343	},	
	{ 2,	7,	10,	 0.14180634400617e-1	},	
	{ 2,	9,	1,	 0.83326504880713e-2	},	
	{ 2,	9,	2,	-0.29052336009585e-1	},	
	{ 2,	9,	3,	 0.38615085574206e-1	},	
	{ 2,	9,	4,	-0.20393486513704e-1	},	
	{ 2,	9,	8,	-0.16554050063734e-2	},	
	{ 2,	10,	6,	 0.19955571979541e-2	},	
	{ 2,	10,	9,	 0.15870308324157e-3	},	
	{ 2,	12,	8,	-0.16388568342530e-4	},	
	{ 3,	3,	16,	 0.43613615723811e-1	},	
	{ 3,	4,	22,	 0.34994005463765e-1	},	
	{ 3,	4,	23,	-0.76788197844621e-1	},	
	{ 3,	5,	23,	 0.22446277332006e-1	},	
	{ 4,	14,	10,	-0.62689710414685e-4	},	
	{ 6,	3,	50,	-0.55711118565645e-9	},	
	{ 6,	6,	44,	-0.19905718354408	},	
	{ 6,	6,	46,	 0.31777497330738	},	
	{ 6,	6,	50,	-0.11841182425981	},	
};
static const struct iapws_phi_coef2 coef2[3] = {
	{ 3, 0, -0.31306260323435e2, 20.0, 150.0, 1.21, 1.0 },
	{ 3, 1,  0.31546140237781e2, 20.0, 150.0, 1.21, 1.0 },
	{ 3, 4, -0.25213154341695e4, 20.0, 250.0, 1.25, 1.0 },
};
static const struct {
	double a;
	double b;
	double B;
	double n;
	double C;
	double D;
	double A;
	double beta;
} coef3[2] = {
	{ 3.5, 0.85, 0.2, -0.14874640856724, 28.0, 700.0, 0.32, 0.3 },
	{ 3.5, 0.95, 0.2,  0.31806110878444, 32.0, 800.0, 0.32, 0.3 },
};

static void calc_Deltab(int i, double delta, double tau, struct iapws_phi *ans)
{
	assert(delta != 1.0 || tau != 1.0);
	const double dm1 = delta - 1.0;
	const double dm2 = POW2(dm1);
	const double invbeta = 1.0 / coef3[i].beta;
	double Adb, Bda, Db;
	struct iapws_phi theta, Delta;

	/* theta */
	Adb = POW(dm2, invbeta * 0.5 - 1.0) * coef3[i].A;
	theta.d00 = Adb * dm2 + 1.0 - tau;
	theta.d10 = Adb * dm1 * invbeta;
	theta.d20 = Adb * invbeta * (invbeta - 1.0);
	theta.d01 = -1.0;
	theta.d02 = 0.0;
	theta.d11 = 0.0;

	/* Delta */
	Bda = POW(dm2, coef3[i].a - 1.0) * coef3[i].B;
	Delta.d00 = Bda * dm2 + POW2(theta.d00);
	Delta.d10 = (Bda * dm1 * coef3[i].a + theta.d00 * theta.d10) * 2.0;
	Delta.d20 = (Bda * coef3[i].a * (coef3[i].a * 2.0 - 1.0) +
			POW2(theta.d10) + theta.d00 * theta.d20) * 2.0;
	Delta.d01 = theta.d00 * (-2.0);
	Delta.d02 = 2.0;
	Delta.d11 = theta.d10 * (-2.0);

	/* Deltab */
	assert(Delta.d00 > 0);  /* FIXME */
	Db = POW(Delta.d00, coef3[i].b - 2.0);
	ans->d00 = Db * POW2(Delta.d00);
	ans->d10 = Db * Delta.d00 * Delta.d10 * coef3[i].b;
	ans->d20 = Db * (POW2(Delta.d10) * (coef3[i].b - 1.0) +
			Delta.d00 * Delta.d20) * coef3[i].b;
	ans->d01 = Db * Delta.d00 * Delta.d01 * coef3[i].b;
	ans->d02 = Db * (POW2(Delta.d01) * (coef3[i].b - 1.0) +
			Delta.d00 * Delta.d02) * coef3[i].b;
	ans->d11 = Db * (Delta.d00 * Delta.d11 + Delta.d10 * Delta.d01 *
			(coef3[i].b - 1.0)) * coef3[i].b;
}

static void calc_psi(int i, double delta, double tau, struct iapws_phi *psi)
{
	const double dm1 = delta - 1.0;
	const double tm1 = tau - 1.0;
	const double dm2 = POW2(dm1);
	const double tm2 = POW2(tm1);
	const double psi0 = exp(-coef3[i].C * dm2 - coef3[i].D * tm2);

	psi->d00 = psi0;
	psi->d10 = psi0 * dm1 * coef3[i].C * (-2.0);
	psi->d20 = psi0 * (dm2 * coef3[i].C * 4.0 - 2.0) * coef3[i].C;
	psi->d01 = psi0 * tm1 * coef3[i].D * (-2.0);                       
	psi->d02 = psi0 * (tm2 * coef3[i].D * 4.0 - 2.0) * coef3[i].D;
	psi->d11 = psi0 * dm1 * tm1 * coef3[i].C * coef3[i].D * 4.0;
}

void iapws95_phi(double rho, double t, struct iapws_phi *phi)
{
	int i;
	const double delta = rho / IAPWS_RHOC;
	const double tau = IAPWS_TC / t;
	double xn;
	struct iapws_phi Db, psi;

	phi->type = IAPWS_PHI;
	phi->rho = rho;
	phi->t = t;
	phi->R = IAPWS95_R;

	iapws_phi(delta, tau, coef0, ARRAY_SIZE(coef0),
			coef1, ARRAY_SIZE(coef1),
			coef2, ARRAY_SIZE(coef2), phi);

	if (delta != 1.0 || tau != 1.0) {  /* FIXME not totally valid */
		for (i = 0; i < ARRAY_SIZE(coef3); ++i) {
			calc_Deltab(i, delta, tau, &Db);
			calc_psi(i, delta, tau, &psi);
			xn = coef3[i].n * delta;
			phi->d00 += xn * Db.d00 * psi.d00;
			phi->d10 += xn * (Db.d00 * (psi.d00 + psi.d10 * delta) +
					Db.d10 * psi.d00 * delta);
			phi->d01 += xn * (Db.d00 * psi.d01 + Db.d01 * psi.d00) * tau;
			phi->d11 += xn * (Db.d00 * (psi.d01 + psi.d11* delta) +
					Db.d10 * psi.d01 * delta +
					Db.d01 * (psi.d00 + psi.d10 * delta) +
					Db.d11 * psi.d00 * delta) * tau;
			phi->d20 += xn * (Db.d00 * (psi.d10 * 2 + psi.d20 * delta) +
					Db.d10 * (psi.d00 + psi.d10 * delta) * 2 +
					Db.d20 * psi.d00 * delta) * delta;
			phi->d02 += xn * (Db.d00 * psi.d02 +
					Db.d01 * psi.d01 * 2 +
					Db.d02 * psi.d00) * POW2(tau);
		}
	}
}

#define PEPS	1.0001
#define REPS	1.01

int iapws95_phi_rhot(double rho, double t, enum iapws_state state,
		struct iapws_phi *phi)
{
	switch (state) {
		case IAPWS_LIQUID:
		case IAPWS_GAS:
		case IAPWS_CRIT:
			iapws95_phi(rho, t, phi);
			return 0;
		default:
			return -1;
	}
	assert(0);
}

int iapws95_phi_pt(double p, double t, enum iapws_state state,
		struct iapws_phi *phi)
{
	struct nroot_control ctrl = nroot_default;

	/* Find suitable start value for rho */
	double rho;
	if (if97_gamma_pt(p, t, state, phi) == 0) {
		rho = iapws_rho(phi);
		if (state == IAPWS_LIQUID) {
			rho *= REPS;
		} else if (state == IAPWS_GAS) {
			rho /= REPS;
		}
	} else {
		if (state == IAPWS_LIQUID) {
			/* Rackett equation */
			rho = (IAPWS_PC * 1.0e3) / (IAPWS_TC * IAPWS95_R);
			rho *= POW(IAPWS_RHOC / rho, 1.0 +
					POW(fabs(1.0 - t / IAPWS_TC), 2.0 / 7.0));
		} else if (state == IAPWS_GAS) {
			/* Ideal gas */
			rho = (p * 1.0e3) / (t * IAPWS95_R);
		} else if (state == IAPWS_CRIT) {
			/* Good enough default value */
			rho = IAPWS_RHOC * 2.0;
		} else {
			return -10;
		}
	}

	phi->p = p;
	phi->t = t;
	struct iapws_phi_call call = { iapws95_phi, phi };
	return nroot1(get_phi_pt, &rho, &call, &ctrl);
}

int iapws95_phi_ph(double p, double h, struct iapws_phi *phi)
{
	struct nroot_control ctrl = nroot_default;

	/* Find suitable start value for rho, t */
	double rhot[2];
	double p0 = p, h0 = h;
	if (p > 1.0e2) p0 = 1.0e2;
	if (h < 0.0) h0 = 0.0;
	else if (h > 4.0e3) h0 = 4.0e3;
	if (if97_gamma_ph(p0, h0, phi) != 0) return -10;
	rhot[0] = iapws_rho(phi);
	rhot[1] = iapws_t(phi);

	phi->p = p;
	phi->h = h;
	struct iapws_phi_call call = { iapws95_phi, phi };
	return nroot2(get_phi_ph, rhot, &call, &ctrl);
}

int iapws95_sat_t(double t, struct iapws_phi *phil, struct iapws_phi *phig)
{
	struct nroot_control ctrl = nroot_default;

	double p = if97_psat(t);
	if (p == 0.0) return -1;  /* test if outside of saturation line */
	if (if97_gamma_pt(p, t, IAPWS_LIQUID, phil) != 0) return -11;
	if (if97_gamma_pt(p, t, IAPWS_GAS, phig) != 0) return -12;
	double x[2] = { iapws_rho(phil) * REPS, iapws_rho(phig) / REPS };

	//phil->t = t;
	//phig->t = t;
	struct iapws_phi_call call[2] = {
		{ iapws95_phi, phil },
		{ iapws95_phi, phig },
	};
	return nroot2(get_sat_t, x, call, &ctrl);
}

int iapws95_sat_p(double p, struct iapws_phi *phil, struct iapws_phi *phig)
{
	struct nroot_control ctrl = nroot_default;

	double t = if97_tsat(p);
	if (t == 0.0) return -1;  /* test if outside of saturation line */
	if (if97_gamma_pt(p, t, IAPWS_LIQUID, phil) != 0) return -11;
	if (if97_gamma_pt(p, t, IAPWS_GAS, phig) != 0) return -12;
	double x[3] = { iapws_rho(phil) * REPS, iapws_rho(phig) / REPS, t };

	//phil->p = p;
	//phig->p = p;
	struct iapws_phi_call call[2] = {
		{ iapws95_phi, phil },
		{ iapws95_phi, phig },
	};
	return nrootn(3, get_sat_p, x, call, &ctrl);
}

enum iapws_state iapws95_state_pt(double p, double t)
{
	double ps;
	struct iapws_phi phil, phig;

	if (t >= IAPWS_TT && t < IAPWS_TC && p < 620.0) {
		/* Try with fast approximation */
		ps = sat86_p(t);
		if (p > ps * PEPS) return IAPWS_LIQUID;
		if (p * PEPS < ps) return IAPWS_GAS;

		iapws95_sat_t(t, &phil, &phig);
		ps = iapws_p(&phig);
		if (p > ps) return IAPWS_LIQUID;
		if (p < ps) return IAPWS_GAS;

		return IAPWS_SAT;
	} else if (t >= IAPWS_TC) {
		if (p < IAPWS_PC) return IAPWS_GAS;
		else return IAPWS_CRIT;
	}
	return melt_sub_state(p, t);
}

enum iapws_state iapws95_state_rhot(double rho, double t)
{
	double rhol, rhog;
	struct iapws_phi phil, phig;

	if (t >= IAPWS_TT && t < IAPWS_TC) {
		/* Try with fast approximation */
		rhol = sat86_rhol(t);
		rhog = sat86_rhog(t);
		if (rho > rhol * REPS) return IAPWS_LIQUID;
		else if (rho * REPS < rhog) return IAPWS_GAS;
		else if (rho * REPS < rhol && rho > rhog * REPS)
			return IAPWS_SAT;

		iapws95_sat_t(t, &phil, &phig);
		rhol = iapws_rho(&phil);
		rhog = iapws_rho(&phig);
		if (rho > rhol) return IAPWS_LIQUID;
		else if (rho < rhog) return IAPWS_GAS;
		else return IAPWS_SAT;
	} else if (t >= IAPWS_TC) {
		//iapws95_phi(rho, t, &phig);
		//if (iapws_p(&phig) < IAPWS_PC) return IAPWS_GAS;
		//else return IAPWS_CRIT;
		if (rho < IAPWS_RHOC) return IAPWS_GAS;
		else return IAPWS_CRIT;
	}
	return IAPWS_UNDEF;  /* FIXME */
}

