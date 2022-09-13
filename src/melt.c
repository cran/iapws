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

//#include <assert.h>
#include <math.h>

#include "iapws.h"
#include "melt.h"

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

static double melt_p1h(double t)
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

static double melt_p3(double t)
{
	double theta = t / tmelt[1];
	if (t < tmelt[1] || t > tmelt[2]) return 0.0;
	return (1.0 - (1.0 - POWINT(theta, 60)) * 0.299948) * pmelt[1];
}

static double melt_p5(double t)
{
	double theta = t / tmelt[2];
	if (t < tmelt[2] || t > tmelt[3]) return 0.0;
	return (1.0 - (1.0 - POWINT(theta, 8)) * 1.18721) * pmelt[2];
}

static double melt_p6(double t)
{
	double theta = t / tmelt[3];
	if (t < tmelt[3] || t > tmelt[4]) return 0.0;
	return (1.0 - (1.0 - POW(theta, 4.6)) * 1.07476) * pmelt[3];
}

static double melt_p7(double t)
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

double melt_p(double t, ice_phase_id phase)
{
	switch (phase) {
		case ICE_1H: return melt_p1h(t);
		case ICE_3:  return melt_p3(t);
		case ICE_5:  return melt_p5(t);
		case ICE_6:  return melt_p6(t);
		case ICE_7:  return melt_p7(t);
		default: return 0.0;
	}
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

iapws_state_id melt_sub_state(double p, double t)
{
	//assert((t < tmelt[0] && p < pmelt[1]) || p > pmelt[1]);

	if (p < 2.0e-46) { /* p < sub_p(50.0) */
		if (t >= 50.0) return IAPWS_GAS;
		else return IAPWS_UNDEF;
	} else if (p < pmelt[0]) {
		if (t < 50.0) return IAPWS_SOLID;
		else if (t > tmelt[0]) return IAPWS_GAS;
		else if (p <= sub_p(t)) return IAPWS_GAS;
		else return IAPWS_SOLID;
	} else if (p < pmelt[1]) {
		if (t < tmelt[1]) return IAPWS_SOLID;
		else if (t > tmelt[0]) return IAPWS_LIQUID;
		else if (p >= melt_p1h(t)) return IAPWS_LIQUID;
		else return IAPWS_SOLID;
	} else if (p < pmelt[2]) {
		if (t < tmelt[1]) return IAPWS_SOLID;
		else if (t > tmelt[2]) return IAPWS_LIQUID;
		else if (p <= melt_p3(t)) return IAPWS_LIQUID;
		else return IAPWS_SOLID;
	} else if (p < pmelt[3]) {
		if (t < tmelt[2]) return IAPWS_SOLID;
		else if (t > tmelt[3]) return IAPWS_LIQUID;
		else if (p <= melt_p5(t)) return IAPWS_LIQUID;
		else return IAPWS_SOLID;
	} else if (p < pmelt[4]) {
		if (t < tmelt[3]) return IAPWS_SOLID;
		else if (t > tmelt[4]) return IAPWS_LIQUID;
		else if (p <= melt_p6(t)) return IAPWS_LIQUID;
		else return IAPWS_SOLID;
	} else {  //if (p < pmelt[5]) {
		if (t < tmelt[4]) return IAPWS_SOLID;
		else if (t > tmelt[5]) return IAPWS_LIQUID;
		else if (p <= melt_p7(t)) return IAPWS_LIQUID;
		else return IAPWS_SOLID;
	}
}

