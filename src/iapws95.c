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

/* International Association for the Properties of Water and Steam,
 * IAPWS R6-95(2018), Revised Release on the IAPWS Formulation 1995 for the
 * Thermodynamic Properties of Ordinary Water Substance for General and
 * Scientific Use (2018)
 */

#include <assert.h>
#include <math.h>

#include "iapws95.h"
#include "if97.h"
#include "nroot.h"

enum {
	SIZE0 = 8,
	SIZE1 = 7,
	SIZE2 = 44,
	SIZE3 = 3,
	SIZE4 = 2,
};

static const struct {
	double n;
	double gamma;
} coef0[SIZE0] = {
	{ -8.3204464837497,  0.00000000 },
	{  6.6832105275932,  0.00000000 },
	{  3.0063200000000,  0.00000000 },
	{  0.0124360000000,  1.28728967 }, 
	{  0.9731500000000,  3.53734222 },
	{  1.2795000000000,  7.74073708 },
	{  0.9695600000000,  9.24437796 },
	{  0.2487300000000, 27.50751050 },
};
static const struct {
	int d;
	double t;
	double n;
} coef1[SIZE1] = {
	{ 1, -0.5,	 0.12533547935523e-1	},
	{ 1,  0.875,	 0.78957634722828e1	},
	{ 1,  1.0,	-0.87803203303561e1	},
	{ 2,  0.5,	 0.31802509345418	},
	{ 2,  0.75,	-0.26145533859358	},
	{ 3,  0.375,	-0.78199751687981e-2	},
	{ 4,  1.0,	 0.88089493102134e-2	},
};
static const struct {
	int c;
	int d;
	int t;
	double n;
} coef2[SIZE2] = {
	{ 1,  1,  4,	-0.66856572307965	},
	{ 1,  1,  6,	 0.20433810950965	},
	{ 1,  1, 12,	-0.66212605039687e-4	},
	{ 1,  2,  1,	-0.19232721156002	},
	{ 1,  2,  5,	-0.25709043003438	},
	{ 1,  3,  4,	 0.16074868486251	},
	{ 1,  4,  2,	-0.40092828925807e-1	},
	{ 1,  4, 13,	 0.39343422603254e-6	},
	{ 1,  5,  9,	-0.75941377088144e-5	},
	{ 1,  7,  3,	 0.56250979351888e-3	},
	{ 1,  9,  4,	-0.15608652257135e-4	},
	{ 1, 10, 11,	 0.11537996422951e-8	},
	{ 1, 11,  4,	 0.36582165144204e-6	},
	{ 1, 13, 13,	-0.13251180074668e-11	},
	{ 1, 15,  1,	-0.62639586912454e-9	},
	{ 2,  1,  7,	-0.10793600908932	},
	{ 2,  2,  1,	 0.17611491008752e-1	},
	{ 2,  2,  9,	 0.22132295167546	},
	{ 2,  2, 10,	-0.40247669763528	},
	{ 2,  3, 10,	 0.58083399985759	},
	{ 2,  4,  3,	 0.49969146990806e-2	},
	{ 2,  4,  7,	-0.31358700712549e-1	},
	{ 2,  4, 10,	-0.74315929710341	},
	{ 2,  5, 10,	 0.47807329915480	},
	{ 2,  6,  6,	 0.20527940895948e-1	},
	{ 2,  6, 10,	-0.13636435110343	},
	{ 2,  7, 10,	 0.14180634400617e-1	},
	{ 2,  9,  1,	 0.83326504880713e-2	},
	{ 2,  9,  2,	-0.29052336009585e-1	},
	{ 2,  9,  3,	 0.38615085574206e-1	},
	{ 2,  9,  4,	-0.20393486513704e-1	},
	{ 2,  9,  8,	-0.16554050063734e-2	},
	{ 2, 10,  6,	 0.19955571979541e-2	},
	{ 2, 10,  9,	 0.15870308324157e-3	},
	{ 2, 12,  8,	-0.16388568342530e-4	},
	{ 3,  3, 16,	 0.43613615723811e-1	},
	{ 3,  4, 22,	 0.34994005463765e-1	},
	{ 3,  4, 23,	-0.76788197844621e-1	},
	{ 3,  5, 23,	 0.22446277332006e-1	},
	{ 4, 14, 10,	-0.62689710414685e-4	},
	{ 6,  3, 50,	-0.55711118565645e-9	},
	{ 6,  6, 44,	-0.19905718354408	},
	{ 6,  6, 46,	 0.31777497330738	},
	{ 6,  6, 50,	-0.11841182425981	},
};
static const struct {
	int d;
	int t;
	double n;
	double alpha;
	double beta;
	double gamma;
	double eps;
} coef3[SIZE3] = {
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
} coef4[SIZE4] = {
	{ 3.5, 0.85, 0.2, -0.14874640856724, 28.0, 700.0, 0.32, 0.3 },
	{ 3.5, 0.95, 0.2,  0.31806110878444, 32.0, 800.0, 0.32, 0.3 },
};

static void calc_Deltab(int i, double delta, double tau, iapws_phi *ans)
{
	assert(delta != 1.0 || tau != 1.0);
	double dm1 = delta - 1.0;
	double invbeta = 1.0 / coef4[i].beta;
	double Adb, Bda, Db;
	iapws_phi theta, Delta;

	/* theta */
	Adb = pow(dm1 * dm1, invbeta * 0.5 - 1.0) * coef4[i].A;
	theta.d00 = Adb * POW2(dm1) + 1.0 - tau;
	theta.d10 = Adb * dm1 * invbeta;
	theta.d20 = Adb * invbeta * (invbeta - 1.0);
	theta.d01 = -1.0;
	theta.d02 = 0.0;
	theta.d11 = 0.0;

	/* Delta */
	Bda = pow(dm1 * dm1, coef4[i].a - 1.0) * coef4[i].B;
	Delta.d00 = Bda * POW2(dm1) + POW2(theta.d00);
	Delta.d10 = (Bda * dm1 * coef4[i].a + theta.d00 * theta.d10) * 2.0;
	Delta.d20 = (Bda * coef4[i].a * (coef4[i].a * 2.0 - 1.0) +
			POW2(theta.d10) + theta.d00 * theta.d20) * 2.0;
	Delta.d01 = theta.d00 * (-2.0);
	Delta.d02 = 2.0;
	Delta.d11 = theta.d10 * (-2.0);

	/* Deltab */
	assert(Delta.d00 > 0);  /* FIXME */
	Db = pow(Delta.d00, coef4[i].b - 2.0);
	ans->d00 = Db * POW2(Delta.d00);
	ans->d10 = Db * Delta.d00 * Delta.d10 * coef4[i].b;
	ans->d20 = Db * (POW2(Delta.d10) * (coef4[i].b - 1.0) +
			Delta.d00 * Delta.d20) * coef4[i].b;
	ans->d01 = Db * Delta.d00 * Delta.d01 * coef4[i].b;
	ans->d02 = Db * (POW2(Delta.d01) * (coef4[i].b - 1.0) +
			Delta.d00 * Delta.d02) * coef4[i].b;
	ans->d11 = Db * (Delta.d00 * Delta.d11 + Delta.d10 * Delta.d01 *
			(coef4[i].b - 1.0)) * coef4[i].b;
}

static void calc_psi(int i, double delta, double tau, iapws_phi *psi)
{
	double dm1 = delta - 1.0;
	double tm1 = tau - 1.0;
	double psi0 = exp(-coef4[i].C * POW2(dm1) - coef4[i].D * POW2(tm1));

	psi->d00 = psi0;
	psi->d10 = psi0 * dm1 * coef4[i].C * (-2.0);
	psi->d20 = psi0 * (POW2(dm1) * coef4[i].C * 4.0 - 2.0) * coef4[i].C;
	psi->d01 = psi0 * tm1 * coef4[i].D * (-2.0);                       
	psi->d02 = psi0 * (POW2(tm1) * coef4[i].D * 4.0 - 2.0) * coef4[i].D;
	psi->d11 = psi0 * dm1 * tm1 * coef4[i].C * coef4[i].D * 4.0;
}

int iapws95_phi(double rho, double t, iapws_phi *phi)
{
	int i;
	double delta = rho / IAPWS_RHOC;
	double tau = IAPWS_TC / t;
	double xn, egt, dc, tc;
	iapws_phi Db, psi;

	phi->type = IAPWS_PHI;
	phi->rho = rho;
	phi->t = t;
	phi->R = IAPWS95_R;

	/* phi0 */
	phi->d00 = log(delta) + coef0[0].n + coef0[1].n * tau + coef0[2].n * log(tau);
	phi->d10 = 1.0;
	phi->d01 = coef0[1].n * tau + coef0[2].n;
	phi->d11 = 0.0;
	phi->d20 = -1.0;
	phi->d02 = -coef0[2].n;
	for (i = 3; i < SIZE0; ++i) {
		egt = expm1(-coef0[i].gamma * tau);
		phi->d00 += coef0[i].n * log(-egt);
		phi->d01 -= coef0[i].n * coef0[i].gamma * (tau + tau / egt);
		phi->d02 -= coef0[i].n * POW2(coef0[i].gamma * tau) *
			(egt + 1) / POW2(egt);
	}

	/* phir */
	for (i = 0; i < SIZE1; ++i) {
		xn = coef1[i].n *
			pow_di(delta, coef1[i].d) * pow(tau, coef1[i].t);
		phi->d00 += xn;
		phi->d10 += xn * coef1[i].d; 
		phi->d01 += xn * coef1[i].t;
		phi->d11 += xn * coef1[i].d * coef1[i].t;
		phi->d20 += xn * coef1[i].d * (coef1[i].d - 1);
		phi->d02 += xn * coef1[i].t * (coef1[i].t - 1);
	}
	for (i = 0; i < SIZE2; ++i) {
		dc = pow_di(delta, coef2[i].c);
		xn = coef2[i].n * exp(-dc) *
			pow_di(delta, coef2[i].d) * pow_di(tau, coef2[i].t);
		phi->d00 += xn;
		phi->d10 += xn * (coef2[i].d - coef2[i].c * dc);
		phi->d01 += xn * coef2[i].t;
		phi->d11 += xn * (coef2[i].d - coef2[i].c * dc) * coef2[i].t;
		phi->d20 += xn * ((coef2[i].d - coef2[i].c * dc) *
				(coef2[i].d - coef2[i].c * dc - 1) -
				coef2[i].c * coef2[i].c * dc);
		phi->d02 += xn * coef2[i].t * (coef2[i].t - 1);
	}
	for (i = 0; i < SIZE3; ++i) {
		dc = delta - coef3[i].eps;
		tc = tau - coef3[i].gamma;
		xn = coef3[i].n *
			pow_di(delta, coef3[i].d) * pow_di(tau, coef3[i].t) *
			exp(-coef3[i].alpha * dc * dc -coef3[i].beta * tc * tc);
		phi->d00 += xn;
		phi->d10 += xn * (coef3[i].d - delta * 2 * coef3[i].alpha * dc);
		phi->d01 += xn * (coef3[i].t - tau * 2 * coef3[i].beta * tc);
		phi->d11 += xn * (coef3[i].d - delta * 2 * coef3[i].alpha * dc) *
			(coef3[i].t - tau * 2 * coef3[i].beta * tc);
		phi->d20 += xn * (POW2(coef3[i].d - 2 * coef3[i].alpha * delta * dc) -
				coef3[i].d - 2 * coef3[i].alpha * POW2(delta));
		phi->d02 += xn * (POW2(coef3[i].t - 2 * coef3[i].beta * tau * tc) -
				coef3[i].t - 2 * coef3[i].beta * POW2(tau));
	}
	if (delta != 1.0 || tau != 1.0) {  /* FIXME not totally valid */
		for (i = 0; i < SIZE4; ++i) {
			calc_Deltab(i, delta, tau, &Db);
			calc_psi(i, delta, tau, &psi);
			xn = coef4[i].n * delta;
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
	return 0;
}

static void get_phi_pt(double *rho, void *xphi, double *p, double *dp)
{
	iapws_phi *phi = (iapws_phi *)(xphi);
	iapws95_phi(*rho, phi->t, phi);
	double fact = phi->R * phi->t * 1e-3;
	*p = phi->d10 * (*rho) * fact - phi->p;
	*dp = (phi->d10 * 2.0 + phi->d20) * fact;
	/*
	*p = iapws_p(phi) - phi->p;
	*dp = 1.0 / (iapws_chit(phi) * (*rho));
	*/
}

#define FTOL (1.01)

int iapws95_phi_pt(double p, double t, iapws_state_id state, iapws_phi *phi)
{
	int iter;
	int maxiter = iapws_nroot_max_iter;
	double tol = iapws_nroot_tol;
	double err;

	/* Find suitable start value for rho */
	/* May be more effective if test t0 before p0 */
	double rho;
	double p0 = p <= 100.0 ? p : 100.0;
	double t0 = t >= 273.15 ? t : 273.15;
	if (p0 <= 50.0) {
		if (t0 > 2273.15) t0 = 2273.15;
	} else {
		if (t0 > 1073.15) t0 = 1073.15;
	}
	if (if97_gamma(p0, t0, state, phi) != 0) return -10;
	rho = iapws_v(phi);
	rho = rho != 0.0 ? 1.0 / rho : IAPWS_RHOC;
	if (state == IAPWS_LIQUID) {
		rho *= FTOL;
	} else if (state == IAPWS_GAS) {
		rho /= FTOL;
	}

	phi->p = p;
	phi->t = t;
	iter = nroot(get_phi_pt, &rho, phi, &err, tol, maxiter);
	if (iter >= 0 && iter < maxiter) return 0;
	else return -1;
}

static void get_sat(double *x, void *xphi, double *fx, double *dfx)
{
	double rhol = x[0];
	double rhog = x[1];

	iapws_phi *phil = *((iapws_phi **)(xphi));
	iapws_phi *phig = *((iapws_phi **)(xphi) + 1);
	iapws95_phi(rhol, phil->t, phil);
	iapws95_phi(rhog, phig->t, phig);

	fx[0] = phil->d10 * rhol - phig->d10 * rhog;
	fx[1] = phil->d00 + phil->d10 - phig->d00 - phig->d10;

	dfx[0] = (phil->d10 * 2.0 + phil->d20);
	dfx[1] = dfx[0] / rhol;
	dfx[2] = -(phig->d10 * 2.0 + phig->d20);
	dfx[3] = dfx[2] / rhog;
}

int iapws95_sat(double t, iapws_phi *phil, iapws_phi *phig)
{
	int iter;
	int maxiter = iapws_nroot_max_iter;
	double tol = iapws_nroot_tol;
	double err;
	iapws_phi *phi[2] = { phil, phig };
	double p = if97_psat(t);
	if (p == 0.0) return -1;  /* test if outside of saturation line */
	if (if97_gamma(p, t, IAPWS_LIQUID, phil) != 0) return -11;
	if (if97_gamma(p, t, IAPWS_GAS, phig) != 0) return -12;
	double x[2] = { iapws_rho(phil) * FTOL, iapws_rho(phig) / FTOL };
	iter = nroot2(get_sat, x, phi, &err, tol, maxiter);
	if (iter >= 0 && iter < maxiter) return 0;
	else return -1;
}

static void get_sat_p(double *x, void *xphi, double *fx, double *dfx)
{
	double rhol = x[0];
	double rhog = x[1];
	double t = x[2];

	iapws_phi *phil = *((iapws_phi **)(xphi));
	iapws_phi *phig = *((iapws_phi **)(xphi) + 1);
	iapws95_phi(rhol, t, phil);
	iapws95_phi(rhog, t, phig);

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

int iapws95_sat_p(double p, iapws_phi *phil, iapws_phi *phig)
{
	int iter;
	int maxiter = iapws_nroot_max_iter;
	double tol = iapws_nroot_tol;
	double err;
	iapws_phi *phi[2] = { phil, phig };
	double t = if97_tsat(p);
	if (t == 0.0) return -1;  /* test if outside of saturation line */
	if (if97_gamma(p, t, IAPWS_LIQUID, phil) != 0) return -11;
	if (if97_gamma(p, t, IAPWS_GAS, phig) != 0) return -12;
	double x[3] = { iapws_rho(phil) * FTOL, iapws_rho(phig) / FTOL, t };
	iter = nrootn(3, get_sat_p, x, phi, &err, tol, maxiter);
	if (iter >= 0 && iter < maxiter) return 0;
	else return -1;
}

iapws_state_id iapws95_state(double p, double t)
{
	double ps;
	iapws_phi phil, phig;

	if (t < 273.16) return IAPWS_SOLID;  /* FIXME */
	if (t >= IAPWS_TC && p >= IAPWS_PC) return IAPWS_CRIT;
	if (t >= IAPWS_TC) return IAPWS_GAS;
	iapws95_sat(t, &phil, &phig);
	ps = iapws_p(&phig);
	if (p > ps) return IAPWS_LIQUID;
	if (p < ps) return IAPWS_GAS;
	return IAPWS_UNDEF;
}

