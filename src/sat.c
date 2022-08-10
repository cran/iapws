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
 * IAPWS SR1-86(1992), Revised Supplementary Release on Saturation of Ordinary
 * Water Substance
 */

#include <math.h>

#include "iapws.h"

static inline double sumpow(double x, int n, const double a[], const int I[])
{
	double ans = 0.0;
	int i;
	for (i = 0; i < n; ++i) ans += a[i] * pow_di(x, I[i]);
	return ans;
}

double sat_p(double t)
{
	const int I[6] = { 2, 3, 6, 7, 8, 15 };
	const double a[6] = {
		-7.85951783,  1.84408259, -11.7866497,
		 22.6807411, -15.9618719,  1.80122502,
	};
	double theta = t / IAPWS_TC;
	if (t < 273.15 || t > IAPWS_TC) return 0.0;
	return exp(sumpow(sqrt(1.0 - theta), 6, a, I) / theta) * IAPWS_PC;
}

double sat_rhol(double t)
{
	const int I[6] = { 1, 2, 5, 16, 43, 110 };
	const double b[6] = {
		 1.99274064,  1.09965342, -0.510839303,
		-1.75493479, -45.5170352, -6.74694450e5,
	};
	double theta = t / IAPWS_TC;
	if (t < 273.15 || t > IAPWS_TC) return 0.0;
	return (sumpow(cbrt(1.0 - theta), 6, b, I) + 1.0) * IAPWS_RHOC;
}

double sat_rhog(double t)
{
	const int I[6] = { 2, 4, 8, 18, 37, 71 };
	const double c[6] = {
		-2.03150240, -2.68302940, -5.38626492,
		-17.2991605, -44.7586581, -63.9201063,
	};
	double theta = t / IAPWS_TC;
	if (t < 273.15 || t > IAPWS_TC) return 0.0;
	return exp(sumpow(pow(1.0 - theta, 1.0 / 6.0), 6, c, I)) * IAPWS_RHOC;
}

/*
static const int d[5] = {
	-5.65134998e-8,
	 2690.66631,
	 127.287297,
	-135.003439,
	 0.981825814,
};

static double alpha(double theta)
{
	double xt = sqrt(theta);
	return -1135.905627715 +
		d[0] * pow_di(theta, -19) +
		d[1] * theta +
		d[2] * pow_di(theta, 4) * xt +
		d[3] * pow_di(theta, 5) +
		d[4] * pow_di(theta, 54) * xt;
}

static double phi(double theta)
{
	double xt = sqrt(theta);
	return 2319.5246 +
		d[0] * (19.0/20.0) * pow_di(theta, -20) +
		d[1] * log(theta) +
		d[2] * (9.0/7.0) * pow_di(theta, 3) * xt +
		d[3] * (5.0/4.0) * pow_di(theta, 4) +
		d[4] * (109.0/107.0) * pow_di(theta, 53) * xt;
}
*/

