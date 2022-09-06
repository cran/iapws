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
 * IAPWS SR1-86(1992), Revised Supplementary Release on Saturation of
 * Ordinary Water Substance
 */

#include <math.h>

#include "iapws.h"

static inline double sum6pow(const double x, const double a[], const int I[])
{
	return  a[0] * POWINT(x, I[0]) + a[1] * POWINT(x, I[1]) +
		a[2] * POWINT(x, I[2]) + a[3] * POWINT(x, I[3]) +
		a[4] * POWINT(x, I[4]) + a[5] * POWINT(x, I[5]);
}

double sat_p(double t)
{
	const int I[6] = { 2, 3, 6, 7, 8, 15 };
	const double a[6] = {
		-7.85951783,  1.84408259, -11.7866497,
		 22.6807411, -15.9618719,  1.80122502,
	};
	double theta = t / IAPWS_TC;
	if (t < 273.16 || t > IAPWS_TC) return 0.0;
	return exp(sum6pow(sqrt(1.0 - theta), a, I) / theta) * IAPWS_PC;
}

double sat_rhol(double t)
{
	const int I[6] = { 1, 2, 5, 16, 43, 110 };
	const double b[6] = {
		 1.99274064,  1.09965342, -0.510839303,
		-1.75493479, -45.5170352, -6.74694450e5,
	};
	double theta = t / IAPWS_TC;
	if (t < 273.16 || t > IAPWS_TC) return 0.0;
	return (sum6pow(cbrt(1.0 - theta), b, I) + 1.0) * IAPWS_RHOC;
}

double sat_rhog(double t)
{
	const int I[6] = { 2, 4, 8, 18, 37, 71 };
	const double c[6] = {
		-2.03150240, -2.68302940, -5.38626492,
		-17.2991605, -44.7586581, -63.9201063,
	};
	double theta = t / IAPWS_TC;
	if (t < 273.16 || t > IAPWS_TC) return 0.0;
	return exp(sum6pow(POW(1.0 - theta, 1.0 / 6.0), c, I)) * IAPWS_RHOC;
}

