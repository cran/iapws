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
 * IAPWS R12-08(2008), Release on the IAPWS Formulation 2008
 * for the Viscosity of Ordinary Water Substance (2008)
 *
 * International Association for the Properties of Water and Steam,
 * IAPWS R15-11(2011), Release on the IAPWS Formulation 2011 for the Thermal
 * Conductivity of Ordinary Water Substance (2011)
 *
 * International Association for the Properties of Water and Steam,
 * IAPWS R17-20, Release on the IAPWS Formulation 2020
 * for the Viscosity of Heavy Water (2020)
 *
 * International Association for the Properties of Water and Steam,
 * IAPWS R18-21, Release on the IAPWS Formulation 2021 for the Thermal
 * Conductivity of Heavy Water (2021)
 */

#include <assert.h>
#include <math.h>

#include "iapws.h"
#include "iapws95.h"
#include "heavy17.h"
#include "coef.h"
#include "pow.h"

static const double nu = 0.630;
static const double gam = 1.239;
static const double xi0 = 0.13;	/* nm */
static const double gam0 = 0.06;
static const double tr = 1.5;

/* Viscosity */

static double eta0w(double tau, const double *coef, const int n)
{
	int i;
	double ti = 1.0;
	double ans = 0.0;
	for (i = 0; i < n; ++i, ti *= tau) ans += coef[i] * ti;
	return 1.0 / sqrt(tau) / ans;
}

static double eta0h(double tau, const double *coef, const int n)
{
	tau = 1.0 / tau;
	int i;
	double ti = 1.0;
	double ans[2] = { 0.0, 0.0 };
	for (i = 0; i < n; ++i, ti *= tau) {
		ans[0] += coef[i] * ti;
		ans[1] += coef[i + n] * ti;
	}
	return sqrt(tau) * ans[0] / ans[1];
}

static double eta1(double delta, double tau,
		const struct Nij *coef, const int n)
{
	double ans = 0.0;

	tau -= 1.0;
	delta -= 1.0;

	COEF_ITERATE2(coef, n, tau, delta, ans += _ans_;);

	return exp((delta + 1.0) * ans);
}

static double eta2(const double dchi, const double qd, const double xisep)
{
	const double x = 0.068;
	const double qc = 1.0 / 1.9;  /* 1/nm */

	assert(dchi > 0);
	assert(qd > 0);

	const double xi = xi0 * POW(dchi / gam0, nu / gam);
	const double xic = qc * xi;
	const double xid = qd * xi;

	double y;
	if (xi <= xisep) {
		y = 0.2 * xic * POWINT(xid, 5) * (1.0 - xic + POW2(xic) -
				765.0 / 504.0 * POW2(xid));
	} else {
		const double psid = acos(1.0 / sqrt(1.0 + POW2(xid)));
		const double w = sqrt(fabs((xic - 1.0) / (xic + 1.0))) *
			tan(psid * 0.5);
		const double Lw = xic > 1 ?
			log((1.0 + w) / (1.0 - w)) :
			atan(fabs(w)) * 2.0;
		const double invxic = 1.0 / xic;
		y = sin(psid * 3.0) / 12.0 -
			sin(psid * 2.0) * .25 * invxic +
			sin(psid) * (POW2(invxic) - 1.25) -
			invxic * (psid * (POW2(invxic) - 1.5) -
					Lw * POW(fabs(POW2(xic) - 1.0), 1.5) *
					POW2(invxic));
	}

	return exp(x * y);
}

double if97_eta(const struct iapws_phi *phi)  /* µPa.s */
{
	const double delta = iapws_rho(phi) / IAPWS_RHOC;
	const double tau = IAPWS_TC / iapws_t(phi);
	const double coef0[4] = {
		 1.67752e-2,
		 2.20462e-2,
		 0.6366564e-2,
		-0.241605e-2,
	};
	const struct Nij coef1[21] = {
		{ 0,	0,	 5.20094e-1	},
		{ 1,	0,	 8.50895e-2	},
		{ 2,	0,	-1.08374	},
		{ 3,	0,	-2.89555e-1	},
		{ 0,	1,	 2.22531e-1	},
		{ 1,	1,	 9.99115e-1	},
		{ 2,	1,	 1.88797	},
		{ 3,	1,	 1.26613	},
		{ 5,	1,	 1.20573e-1	},
		{ 0,	2,	-2.81378e-1	},
		{ 1,	2,	-9.06851e-1	},
		{ 2,	2,	-7.72479e-1	},
		{ 3,	2,	-4.89837e-1	},
		{ 4,	2,	-2.57040e-1	},
		{ 0,	3,	 1.61913e-1	},
		{ 1,	3,	 2.57399e-1	},
		{ 0,	4,	-3.25372e-2	},
		{ 3,	4,	 6.98452e-2	},
		{ 4,	5,	 8.72102e-3	},
		{ 3,	6,	-4.35673e-3	},
		{ 5,	6,	-5.93264e-4	},
	};
	return eta0w(tau, coef0, ARRAY_SIZE(coef0)) *
		eta1(delta, tau, coef1, ARRAY_SIZE(coef1));
}

double iapws95_eta(const struct iapws_phi *phi)  /* µPa.s */
{
	const double delta = iapws_rho(phi) / IAPWS_RHOC;
	const double tau = IAPWS_TC / iapws_t(phi);
	double eta = if97_eta(phi);

	const double qd = 1.0 / 1.1;  /* 1/nm */
	const double xisep = 0.3817016416;

	struct iapws_phi phir;
	iapws95_phi(delta * IAPWS_RHOC, tr * IAPWS_TC, &phir);
	double dchi = (iapws_kappat(phi) - iapws_kappat(&phir) * tr * tau) *
		IAPWS_PC * POW2(delta);

	if (dchi > 0.0) eta *= eta2(dchi, qd, xisep);

	return eta;
}

double heavy17_eta(const struct iapws_phi *phi)  /* µPa.s */
{
	const double delta = iapws_rho(phi) / HEAVY17_RHOC;
	const double tau = HEAVY17_TC / iapws_t(phi);

	const double coef0[10] = {
		0.889754, 61.22217, -44.8866, 111.5812, 3.547412,
		0.79637, 2.38127, -0.33463, 2.669, 2.11366e-4,
	};
	const struct Nij coef1[25] = {
		{ 0,	0,	 0.510953	},
		{ 2,	0,	-0.558947	},
		{ 3,	0,	-2.718820	},
		{ 4,	0,	 0.480990	},
		{ 5,	0,	 2.404510	},
		{ 6,	0,	-1.824320	},
		{ 0,	1,	 0.275847	},
		{ 1,	1,	 0.762957	},
		{ 3,	1,	 1.760340	},
		{ 4,	1,	 0.0819086	},
		{ 6,	1,	 1.417750	},
		{ 0,	2,	-0.228148	},
		{ 1,	2,	-0.321497	},
		{ 5,	2,	-2.302500	},
		{ 0,	3,	 0.0661035	},
		{ 1,	3,	 0.0449393	},
		{ 2,	3,	 1.466670	},
		{ 5,	3,	 0.938984	},
		{ 6,	3,	-0.108354	},
		{ 0,	4,	-0.00481265	},
		{ 2,	4,	-1.545710	},
		{ 3,	4,	-0.0570938	},
		{ 5,	4,	-0.0753783	},
		{ 2,	5,	 0.553080	},
		{ 2,	6,	-0.0650201	},
	};
	double eta = eta0h(tau, coef0, ARRAY_SIZE(coef0) / 2) *
		eta1(delta, tau, coef1, ARRAY_SIZE(coef1));

	const double qd = 1.0 / 0.4;  /* 1/nm */
	const double xisep = 0.03021806692;

	struct iapws_phi phir;
	heavy17_phi(delta * HEAVY17_RHOC, tr * HEAVY17_TC, &phir);
	double dchi = (iapws_kappat(phi) - iapws_kappat(&phir) * tr * tau) *
		HEAVY17_PC * POW2(delta);

	if (dchi > 0.0) eta *= eta2(dchi, qd, xisep);

	return eta;
}

/*
 * Thermal conductivity
 */

static double lambda01(double delta, double tau)	/* mW/K/m */
{
	const double coef0[5] = {
		2.443221e-3,
		1.323095e-2,
		6.770357e-3,
		-3.454586e-3,
		4.096266e-4,
	};
	const struct Nij coef1[28] = {
		{ 0,	0,	 1.60397357	},
		{ 0,	1,	-0.646013523	},
		{ 0,	2,	 0.111443906	},
		{ 0,	3,	 0.102997357	},
		{ 0,	4,	-0.0504123634	},
		{ 0,	5,	 0.00609859258	},
		{ 1,	0,	 2.33771842	},
		{ 1,	1,	-2.78843778	},
		{ 1,	2,	 1.53616167	},
		{ 1,	3,	-0.463045512	},
		{ 1,	4,	 0.0832827019	},
		{ 1,	5,	-0.00719201245	},
		{ 2,	0,	 2.19650529	},
		{ 2,	1,	-4.54580785	},
		{ 2,	2,	 3.55777244	},
		{ 2,	3,	-1.40944978	},
		{ 2,	4,	 0.275418278	},
		{ 2,	5,	-0.0205938816	},
		{ 3,	0,	-1.21051378	},
		{ 3,	1,	 1.60812989	},
		{ 3,	2,	-0.621178141	},
		{ 3,	3,	 0.0716373224	},
		{ 4,	0,	-2.7203370	},
		{ 4,	1,	 4.57586331	},
		{ 4,	2,	-3.18369245	},
		{ 4,	3,	 1.1168348	},
		{ 4,	4,	-0.19268305	},
		{ 4,	5,	 0.012913842	},
	};
	return eta0w(tau, coef0, ARRAY_SIZE(coef0)) *
		eta1(delta, tau, coef1, ARRAY_SIZE(coef1));
}

static double lambda2(double delta, double tau, double cp, double cv,
		double dchi, double eta, double lam, double qd)	/* mW/K/m */
{
	if (delta == 0.0 || dchi <= 0.0) return 0.0;

	const double y = qd * xi0 * POW(dchi / gam0, nu / gam);
	if (y < 1.2e-7) return 0.0;

	const double invk = cv / cp;
	return lam * delta * cp * 2.0 * M_1_PI / (y * eta * tau) *
		((1.0 - invk) * atan(y) + invk * y +
		 expm1(-y / (1.0 + y*y*y / (POW2(delta) * 3.0))));
}

double if97_lambda(const struct iapws_phi *gamma)	/* mW/K/m */
{
	const double delta = iapws_rho(gamma) / IAPWS_RHOC;
	const double tau = IAPWS_TC / iapws_t(gamma);

	const double A[6][5] = {
		{  6.53786807199516,  6.52717759281799,  5.35500529896124, 
		   1.55225959906681,  1.11999926419994   },
		{ -5.61149954923348, -6.30816983387575, -3.96415689925446, 
		   0.464621290821181,  0.595748562571649 },
		{  3.39624167361325,  8.08379285492595,  8.91990208918795, 
		   8.93237374861479,  9.88952565078920   },
		{ -2.27492629730878, -9.82240510197603, -12.0338729505790, 
		  -11.0321960061126, -10.3255051147040   },
		{  10.2631854662709,  12.1358413791395,  9.19494865194302, 
		   6.16780999933360,  4.66861294457414   },
		{  1.97815050331519, -5.54349664571295, -2.16866274479712, 
		  -0.965458722086812, -0.503243546373828 },
	};
	int i, j;
	double invzr = 0.0;
	double ri;
	if (delta <= 0.310559006) j = 0;
	else if (delta <= 0.776397516) j = 1;
	else if (delta <= 1.242236025) j = 2;
	else if (delta <= 1.863354037) j = 3;
	else j = 4;
	for (i = 0, ri = 1.0; i < 6; ++i, ri *= delta) {
		invzr += A[i][j] * ri;
	}

	double dchi = (iapws_kappat(gamma) * IAPWS_PC * delta -
			tr * tau / invzr) * delta;
	const double lam = 177.8514;
	const double qd = 1.0 / 0.4;
	return lambda01(delta, tau) + lambda2(delta, tau,
			iapws_cp(gamma) / IAPWS95_R,
			iapws_cv(gamma) / IAPWS95_R,
			dchi, if97_eta(gamma), lam, qd);
}

double iapws95_lambda(const struct iapws_phi *phi)	/* mW/K/m */
{
	const double delta = iapws_rho(phi) / IAPWS_RHOC;
	const double tau = IAPWS_TC / iapws_t(phi);
	struct iapws_phi phir;
	iapws95_phi(delta * IAPWS_RHOC, tr * IAPWS_TC, &phir);
	double dchi = (iapws_kappat(phi) - iapws_kappat(&phir) * tr * tau) *
		IAPWS_PC * POW2(delta);

	const double lam = 177.8514;
	const double qd = 1.0 / 0.4;

	return lambda01(delta, tau) + lambda2(delta, tau,
			iapws_cp(phi) / IAPWS95_R, iapws_cv(phi) / IAPWS95_R,
			dchi, iapws95_eta(phi), lam, qd);
}

double heavy17_lambda(const struct iapws_phi *phi)	/* mW/K/m */
{
	const double delta = iapws_rho(phi) / HEAVY17_RHOC;
	const double tau = HEAVY17_TC / iapws_t(phi);

	const double coef0[8] = {
		1.0, 3.3620798, -1.0191198, 2.8518117,
		0.10779213, -0.034637234, 0.036603464, 0.0091018912,
	};
	const struct Nij coef1[30] = {
		{ 0,	0,	 1.50933576	},
		{ 0,	1,	-0.65831078	},
		{ 0,	2,	 0.111174263	},
		{ 0,	3,	 0.140185152	},
		{ 0,	4,	-0.0656227722	},
		{ 0,	5,	 0.00785155213	},
		{ 1,	0,	 2.8414715	},
		{ 1,	1,	-2.9826577	},
		{ 1,	2,	 1.34357932	},
		{ 1,	3,	-0.599233641	},
		{ 1,	4,	 0.28116337	},
		{ 1,	5,	-0.0533292833	},
		{ 2,	0,	 4.86095723	},
		{ 2,	1,	-6.19784468	},
		{ 2,	2,	 2.20941867	},
		{ 2,	3,	 0.224691518	},
		{ 2,	4,	-0.322191265	},
		{ 2,	5,	 0.0596204654	},
		{ 3,	0,	 2.06156007	},
		{ 3,	1,	-3.48612456	},
		{ 3,	2,	 1.47962309	},
		{ 3,	3,	 0.625101458	},
		{ 3,	4,	-0.56123225	},
		{ 3,	5,	 0.0974446139	},
		{ 4,	0,	-2.06105687	},
		{ 4,	1,	 0.416240028	},
		{ 4,	2,	 2.92524513	},
		{ 4,	3,	-2.81703583	},
		{ 4,	4,	 1.00551476	},
		{ 4,	5,	-0.127884416	},
	};
	double lambda01 = eta0h(tau, coef0, ARRAY_SIZE(coef0) / 2) *
		eta1(delta, tau, coef1, ARRAY_SIZE(coef1));

	struct iapws_phi phir;
	heavy17_phi(delta * HEAVY17_RHOC, tr * HEAVY17_TC, &phir);
	double dchi = (iapws_kappat(phi) - iapws_kappat(&phir) * tr * tau) *
		HEAVY17_PC * POW2(delta);

	const double lam = 175.9870;
	const double qd = 1.0 / 0.36;

	return lambda01 + lambda2(delta, tau, iapws_cp(phi) / HEAVY17_R,
			iapws_cv(phi) / HEAVY17_R, dchi, heavy17_eta(phi),
			lam, qd);
}

