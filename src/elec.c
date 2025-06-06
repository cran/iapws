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
 * IAPWS R8-97, Release on the Static Dielectric Constant of Ordinary
 * Water Substance for Temperatures from 238 K to 873 K and Pressures
 * up to 1000 MPa (1997)
 *
 * International Association for the Properties of Water and Steam,
 * IAPWS R9-97, Release on the Refractive Index of Ordinary Water
 * Substance as Function of Wavelength, Temperature and Pressure (1997)
 *
 * International Association for the Properties of Water and Steam,
 * IAPWS R11-24, Revised Release on the Ionization Constant of H2O
 * (2024)
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

	const struct {
		double N;
		int i;
		double j;
	} coef[12] = {
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
	for (i = 0; i < ARRAY_SIZE(coef) - 1; ++i) {
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

static double rind(double rho, double t, double lambda)
{
	const double a[8] = {
		 0.244257733,
		 9.74634476e-3,
		-3.73234996e-3,
		 2.68678472e-4,
		 1.58920570e-3,
		 2.45934259e-3,
		 0.900704920,
		-1.66626219e-2,
	};
	const double luv2 = POW2(0.2292020);
	const double lir2 = POW2(5.432937);

	/* reduced parameters */
	rho *= 1.0e-3;
	t /= 273.15;
	lambda /= 0.589;

	const double l2 = POW2(lambda);
	const double ar = (a[0] + a[1] * rho + a[2] * t + a[3] * t * l2 +
			a[4] / l2 + a[5] / (l2 - luv2) + a[6] / (l2 - lir2) +
			a[7] * POW2(rho)) * rho;
	return sqrt((ar * 2.0 + 1.0) / (1.0 - ar));
}

static double pk(double rho, double t)
{
	int const n = 6;
	double const gamma[4] = {
		 6.141500e-1,
		 4.825133e4,
		-6.770793e4,
		 1.010210e7,
	};
	double const alpha[3] = {
		-0.702132,
		+8681.05,
		-24145.1,
	};
	double const beta[3] = {
		+0.813876,
		-51.4471,
		-0.469920,
	};
	double const lMG = log10(IAPWS_M) - 3.0;

	rho = rho * 1e-3;  /* g/cm3 */
	t = 1.0 / t;
	double const q = rho * exp(alpha[0] + t * (alpha[1] + t * alpha[2] *
			POW(rho, 2.0 / 3.0)));
	return gamma[0] + t * (gamma[1] + t * (gamma[2] + t * gamma[3])) +
		2.0 * (lMG - n * (log10(1.0 + q) - q / (1.0 + q) * rho *
			(beta[0] + beta[1] * t + beta[2] * rho)));
}

double iapws_epsilon(double rho, double t)
{
	return dielec(rho, t);
}

double iapws_n(double rho, double t, double lambda)
{
	return rind(rho, t, lambda);
}

double iapws_pk(double rho, double t)
{
	return pk(rho, t);
}

