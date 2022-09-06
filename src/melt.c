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
 * IAPWS R14-08(2011), Revised Release on the Pressure along the Melting
 * and Sublimation Curves of Ordinary Water Substance
 */

#include <math.h>

#include "iapws.h"  /* for POW and POWINT */

static const double tmelt[] = {
	273.16,
	251.165,
	256.164,
	273.31,
	355.0,
	715.0,
};

static const double pmelt[] = {
	611.657e-6,	/* does not match IAPWS95_PT */
	208.566,
	350.1,
	632.4,
	2216.0,
};

double melt_p1h(double t)
{
	const double a[3] = { 0.119539337e7, 0.808183159e5, 0.333826860e4 };
	const double b[3] = { 0.300000e1, 0.257500e2, 0.103750e3 };

	double theta = t / tmelt[0];
	double ans;

	if (t < tmelt[1] || t > tmelt[0]) return 0.0;
	ans = 1.0 + a[0] * (1.0 - POW(theta, b[0])) +
		    a[1] * (1.0 - POW(theta, b[1])) +
		    a[2] * (1.0 - POW(theta, b[2]));
	return ans * pmelt[0];
}

double melt_p3(double t)
{
	double theta = t / tmelt[1];
	if (t < tmelt[1] || t > tmelt[2]) return 0.0;
	return (1.0 - POWINT(theta, 60) * 0.299948) * pmelt[1];
}

double melt_p5(double t)
{
	double theta = t / tmelt[2];
	if (t < tmelt[2] || t > tmelt[3]) return 0.0;
	return (1.0 - POWINT(theta, 8) * 1.18721) * pmelt[2];
}

double melt_p6(double t)
{
	double theta = t / tmelt[3];
	if (t < tmelt[3] || t > tmelt[4]) return 0.0;
	return (1.0 - POW(theta, 4.6) * 1.07476) * pmelt[3];
}

double melt_p7(double t)
{
	const double a[3] = { 0.173683e1, -0.544606e-1, 0.806106e-7 };
	const int b[3] = { -1, 5, 22 };

	double theta = t / tmelt[4];
	double ans;

	if (t < tmelt[4] || t > tmelt[5]) return 0.0;
	ans =   a[0] * (1.0 - POWINT(theta, b[0])) +
		a[1] * (1.0 - POWINT(theta, b[1])) +
		a[2] * (1.0 - POWINT(theta, b[2]));
	return exp(ans) * pmelt[4];
}

double sub_p(double t)
{
	const double a[3] = { -0.212144006e2, 0.273203819e2, -0.610598130e1 };
	const double b[3] = { 0.333333333e-2, 0.120666667e1, 0.170333333e1 };

	double theta = t / tmelt[0];
	double ans;

	if (t < 50.0 || t > tmelt[0]) return 0.0;
	ans =   a[0] * POW(theta, b[0]) +
		a[1] * POW(theta, b[1]) +
		a[2] * POW(theta, b[2]);
	return exp(ans / theta) * pmelt[0];
}

