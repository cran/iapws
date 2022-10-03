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
 * IAPWS R15-11(2011), Release on the IAPWS Formulation 2011 for the Thermal
 * Conductivity of Ordinary Water Substance (2011)
 */

#include <math.h>

#include "iapws.h"
#include "iapws95.h"
#include "pow.h"
#include "visc.h"

static double lambda01(double rho, double t)	/* mW/K/m */
{
	enum {
		SIZE0 = 5,
		SIZE1 = 5,
		SIZE2 = 6,
	};
	const double coef0[SIZE0] = {
		2.443221e-3,
		1.323095e-2,
		6.770357e-3,
		-3.454586e-3,
		4.096266e-4,
	};
	const double coef1[SIZE1][SIZE2] = {
		{  1.60397357,	-0.646013523,	 0.111443906,
		   0.102997357,	-0.0504123634,	 0.00609859258	},
		{  2.33771842,	-2.78843778,	 1.53616167,
		  -0.463045512,	 0.0832827019,	-0.00719201245	},
		{  2.19650529,	-4.54580785,	 3.55777244,
		  -1.40944978,	 0.275418278,	-0.0205938816	},
		{ -1.21051378,	 1.60812989,	-0.621178141,
		   0.0716373224, 0.0,		 0.0		},
		{ -2.7203370,	 4.57586331,	-3.18369245,
		   1.1168348,	-0.19268305,	 0.012913842	},
	};

	const double tau = IAPWS_TC / t;
	const double delta = rho / IAPWS_RHOC;
	double lambda0 = 0.0;
	double lambda1 = 0.0;
	int i, j;
	double ti, dj;

	for (i = 0, ti = 1.0; i < SIZE0; ++i, ti *= tau) {
		lambda0 += coef0[i] * ti;
	}
	for (i = 0, ti = 1.0; i < SIZE1; ++i, ti *= tau - 1.0) {
		for (j = 0, dj = 1.0; j < SIZE2; ++j, dj *= delta - 1.0) {
			lambda1 += coef1[i][j] * ti * dj;
		}
	}
	return 1.0 / sqrt(tau) / lambda0 * exp(delta * lambda1);
}

static double lambda2(double rho, double t, double cp, double cv,
		double dchi, double eta)	/* mW/K/m */
{
	const double lam = 177.8514;
	const double qd = 1.0 / 0.4;	/* 1/nm */
	const double nu = 0.630;
	const double gam = 1.239;
	const double xi0 = 0.13;	/* nm */
	const double gam0 = 0.06;

	if (rho == 0.0 || dchi <= 0.0) return 0.0;

	const double y = qd * xi0 * POW(dchi / gam0, nu / gam);
	if (y < 1.2e-7) return 0.0;

	const double delta = rho / IAPWS_RHOC;
	const double invk = cv / cp;
	return lam * delta * t * cp * 2.0 * M_1_PI *
		((1.0 - invk) * atan(y) + invk * y +
		 expm1(-y / (1.0 + y*y*y / (POW2(delta) * 3.0)))) /
		(y * eta * IAPWS95_R * IAPWS_TC);
}

double if97_lambda(const iapws_phi *gamma)	/* mW/K/m */
{
	const double rho = iapws_rho(gamma);
	const double t = iapws_t(gamma);

	const double tr = 1.5 * IAPWS_TC;
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
	const double rhob = rho / IAPWS_RHOC;
	double invzr = 0.0;
	double ri;
	if (rhob <= 0.310559006) j = 0;
	else if (rhob <= 0.776397516) j = 1;
	else if (rhob <= 1.242236025) j = 2;
	else if (rhob <= 1.863354037) j = 3;
	else j = 4;
	for (i = 0, ri = 1.0; i < 6; ++i, ri *= rhob) {
		invzr += A[i][j] * ri;
	}

	double dchi = (iapws_kappat(gamma) * IAPWS_PC * rhob -
			tr / (t * invzr)) * rhob;
	return lambda01(rho, t) + lambda2(rho, t, iapws_cp(gamma),
			iapws_cv(gamma), dchi, if97_eta(gamma));
}

double iapws95_lambda(const iapws_phi *phi)	/* mW/K/m */
{
	const double rho = iapws_rho(phi);
	const double t = iapws_t(phi);
	const double tr = 1.5 * IAPWS_TC;
	iapws_phi phir;
	iapws95_phi(rho, tr, &phir);
	double dchi = (iapws_kappat(phi) - iapws_kappat(&phir) * tr / t) *
		IAPWS_PC * POW2(rho / IAPWS_RHOC);
	return lambda01(rho, t) + lambda2(rho, t, iapws_cp(phi),
			iapws_cv(phi), dchi, iapws95_eta(phi));
}

