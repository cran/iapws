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
 * IAPWS SR1-86(1992), Revised Supplementary Release on Saturation of
 * Ordinary Water Substance
 */

#include <math.h>

#include "iapws.h"
#include "nroot.h"
#include "coef.h"
#include "pow.h"

static double sumpow6(const double x, const struct Ni coef[6])
{
	double ans = 0.0;
	COEF_ITERATE1(coef, 6, x, ans += _ans_);
	return ans;
}

static double sat_logpi(double theta)
{
	const struct Ni a[6] = {
		{ 2,	-7.85951783 },
		{ 3,	 1.84408259 },
		{ 6,	-11.7866497 },
		{ 7,	 22.6807411 },
		{ 8,	-15.9618719 },
		{ 15,	 1.80122502 },
	};
	return sumpow6(sqrt(1.0 - theta), a) / theta;
}

double sat86_p(double t)
{
	if (t < IAPWS_TT || t > IAPWS_TC) return 0.0;
	return exp(sat_logpi(t / IAPWS_TC)) * IAPWS_PC;
}

static void get_thetasat(double *theta, void *logpi, double *fx, double *dfx)
{
	double *lp = (double *)(logpi);
	*fx = sat_logpi(*theta) - (*lp);
}

double sat86_t(double p)
{
	if (p < IAPWS_PT || p > IAPWS_PC) return 0.0;

	struct nroot_control ctrl = nroot_default;

	double logpi = log(p / IAPWS_PC);
	double theta = 1.0 / (1.0 - 1.416488e-01 * logpi -
			1.047873e-03 * POW2(logpi));
	if (sroot(get_thetasat, &theta, &logpi, &ctrl) != NROOT_SUCCESS)
		return 0.0;
	theta *= IAPWS_TC;
	if (theta > IAPWS_TC) theta = IAPWS_TC;
	else if (theta < IAPWS_TT) theta = IAPWS_TT;
	return theta;
}

double sat86_rhol(double t)
{
	const struct Ni b[6] = {
		{ 1,	 1.99274064	},
		{ 2,	 1.09965342	},
		{ 5,	-0.510839303	},
		{ 16,	-1.75493479	},
		{ 43,	-45.5170352	},
		{ 110,	-6.74694450e5	},
	};
	if (t < IAPWS_TT || t > IAPWS_TC) return 0.0;
	return (sumpow6(cbrt(1.0 - t / IAPWS_TC), b) + 1.0) * IAPWS_RHOC;
}

double sat86_rhog(double t)
{
	const struct Ni c[6] = {
		{ 2,	-2.03150240 },
		{ 4,	-2.68302940 },
		{ 8,	-5.38626492 },
		{ 18,	-17.2991605 },
		{ 37,	-44.7586581 },
		{ 71,	-63.9201063 },
	};
	if (t < IAPWS_TT || t > IAPWS_TC) return 0.0;
	return exp(sumpow6(POW(1.0 - t / IAPWS_TC, 1.0 / 6.0), c)) * IAPWS_RHOC;
}

