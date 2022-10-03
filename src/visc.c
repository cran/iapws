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
 * IAPWS R12-08(2008), Release on the IAPWS Formulation 2008
 * for the Viscosity of Ordinary Water Substance (2008)
 */

#include <math.h>

#include "iapws.h"
#include "iapws95.h"
#include "pow.h"

static double eta01(double rho, double t)
{
	enum {
		SIZE0 = 4,
		SIZE1 = 21,
	};
	const double coef0[SIZE0] = {
		1.67752,
		2.20462,
		0.6366564,
		-0.241605,
	};
	const struct {
		int i;
		int j;
		double H;
	} coef1[SIZE1] = {
		{ 0, 0,  5.20094e-1 },
		{ 1, 0,  8.50895e-2 },
		{ 2, 0, -1.08374    },
		{ 3, 0, -2.89555e-1 },
		{ 0, 1,  2.22531e-1 },
		{ 1, 1,  9.99115e-1 },
		{ 2, 1,  1.88797    },
		{ 3, 1,  1.26613    },
		{ 5, 1,  1.20573e-1 },
		{ 0, 2, -2.81378e-1 },
		{ 1, 2, -9.06851e-1 },
		{ 2, 2, -7.72479e-1 },
		{ 3, 2, -4.89837e-1 },
		{ 4, 2, -2.57040e-1 },
		{ 0, 3,  1.61913e-1 },
		{ 1, 3,  2.57399e-1 },
		{ 0, 4, -3.25372e-2 },
		{ 3, 4,  6.98452e-2 },
		{ 4, 5,  8.72102e-3 },
		{ 3, 6, -4.35673e-3 },
		{ 5, 6, -5.93264e-4 },
	};

	const double tau = IAPWS_TC / t;
	const double delta = rho / IAPWS_RHOC;
	double eta0 = 0.0;
	double eta1 = 0.0;
	int i;

	for (i = 0; i < SIZE0; ++i) {
		eta0 += coef0[i] * powint(tau, i);
	}
	for (i = 0; i < SIZE1; ++i) {
		eta1 += coef1[i].H *
			powint(tau - 1.0, coef1[i].i) *
			powint(delta - 1.0, coef1[i].j);
	}
	return 100.0 / sqrt(tau) / eta0 * exp(delta * eta1);
}

double if97_eta(const iapws_phi *phi)  /* µPa.s */
{
	return eta01(iapws_rho(phi), iapws_t(phi));
}

double iapws95_eta(const iapws_phi *phi)  /* µPa.s */
{
	const double x = 0.068;
	const double qc = 1.0 / 1.9;  /* 1/nm */
	const double qd = 1.0 / 1.1;  /* 1/nm */
	const double nu = 0.630;
	const double gam = 1.239;
	const double xi0 = 0.13;  /* nm */
	const double gam0 = 0.06;
	const double tr = 1.5 * IAPWS_TC;

	const double rho = iapws_rho(phi);
	const double t = iapws_t(phi);
	double eta = eta01(rho, t);

	iapws_phi phir;
	iapws95_phi(rho, tr, &phir);
	double dchi = (iapws_kappat(phi) - iapws_kappat(&phir) * tr / t) *
		IAPWS_PC * POW2(rho / IAPWS_RHOC);

	if (dchi > 0.0) {
		double xi = xi0 * POW(dchi / gam0, nu / gam);
		double xic = qc * xi;
		double xid = qd * xi;

		double y;
		if (xi <= 0.3817016416) {
			y = 0.2 * xic * POWINT(xid, 5) * (1.0 - xic + POW2(xic) -
					765.0 / 504.0 * POW2(xid));
		} else {
			double psid = acos(1.0 / sqrt(1.0 + POW2(xid)));
			double w = sqrt(fabs((xic - 1.0) / (xic + 1.0))) *
				tan(psid * 0.5);
			double Lw = xic > 1 ?
				log((1.0 + w) / (1.0 - w)) :
				atan(fabs(w)) * 2.0;
			double invxic = 1.0 / xic;
			y = sin(psid * 3.0) / 12.0 -
				sin(psid * 2.0) * .25 * invxic +
				sin(psid) * (POW2(invxic) - 1.25) -
				invxic * (psid * (POW2(invxic) - 1.5) -
				Lw * POW(fabs(POW2(xic) - 1.0), 1.5) *
				POW2(invxic));
		}
		eta *= exp(x * y);
	}

	return eta;
}

