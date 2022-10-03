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
 * IAPWS R8-97, Release on the Static Dielectric Constant of Ordinary
 * Water Substance for Temperatures from 238 K to 873 K and Pressures
 * up to 1000 MPa (1997)
 */

#include "iapws.h"
#include "pow.h"

static double dielec(double rho, double t)
{
	const double Na = 6.0221367e23 / 0.018015268;
	const double inve0 = 4e-7 * M_PI * POW2(299792458.0);
	const double alpha = 1.636e-40;
	const double mu2 = POW2(6.138e-30);
	const double k = 1.380658e-23;

	enum { SIZE = 12 };
	const struct {
		double N;
		int i;
		double j;
	} coef[SIZE] = {
		{  0.978224486826e-0,	 1,	0.25	},
		{ -0.957771379375e-0,	 1,	1.00	},
		{  0.237511794148e-0,	 1,	2.50	},
		{  0.714692244396e-0,	 2,	1.50	},
		{ -0.298217036956e-0,	 3,	1.50	},
		{ -0.108863472196e-0,	 3,	2.50	},
		{  0.949327488264e-1,	 4,	2.00	},
		{ -0.980469816509e-2,	 5,	2.00	},
		{  0.165167634970e-4,	 6,	5.00	},
		{  0.937359795772e-4,	 7,	0.50	},
		{ -0.123179218720e-9,	10,	10.0	},
		{  0.196096504426e-2,	 1,	-1.2	},
	};

	const double delta = rho / IAPWS_RHOC;
	const double tau = IAPWS_TC / t;

	int i;
	double a, b;
	double g = 1.0 + coef[11].N * delta * POW(t / 228.0 - 1.0, coef[11].j);
	for (i = 0; i < SIZE - 1; ++i) {
		g += coef[i].N *
			POWINT(delta, coef[i].i) *
			POW(tau, coef[i].j);
	}

	a = Na * inve0 * mu2 * rho * g / (k * t);
	b = Na * inve0 * alpha * rho / (3.0);

	return (1.0 + a + b * 5.0 + sqrt(9.0 + a * 2.0 + b * 18.0 +
				a * (a + b * 10.0) + b * b * 9.0)) /
			(4.0 - b * 4.0);
}

double iapws_epsilon(const iapws_phi *phi)  /* adimensional */
{
	return dielec(iapws_rho(phi), iapws_t(phi));
}
