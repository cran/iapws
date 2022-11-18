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
 * IAPWS R10-06(2009), Revised Release on the Equation of State 2006
 * for H2O Ice Ih (2009)
 */

#include <complex.h>
#include "iapws.h"
#include "pow.h"  /* for ARRAY_SIZE */

void ice06_gamma(double p, double t, struct iapws_phi *gamma)
{
	const double coef0[5] = {
		-0.632020233335886e6	/ IAPWS_TT,
		 0.655022213658955	/ IAPWS_TT,
		-0.189369929326131e-7	/ IAPWS_TT,
		 0.339746123271053e-14	/ IAPWS_TT,
		-0.556464869058991e-21	/ IAPWS_TT,
	};
	const double s0 = -0.332733756492168e4;	/* IAPWS-95 */
	//const double s0 = 0.18913e3;		/* absolute */
	const double complex coeft[2] = {
		0.368017112855051e-1 + I * 0.510878114959572e-1,
		0.337315741065416    + I * 0.335449415919309,
	};
	const double complex r1 = 0.447050716285388e2 + I * 0.656876847463481e2;
	const double complex coefr[3] = {
		-0.725974574329220e2   - I * 0.781008427112870e2,
		-0.557107698030123e-4  + I * 0.464578634580806e-4,
		 0.234801409215913e-10 - I * 0.285651142904972e-10,
	};

	const double pi = p / IAPWS_PT;
	const double p0 = IAPWS_PN / IAPWS_PT;
	const double tau = IAPWS_TT / t;
	const double theta = 1.0 / tau;

	double pik;
	double complex r2[3] = { 0 };
	double complex lt[2][3] = {
		{ clog(coeft[0]), clog(coeft[0]-theta), clog(coeft[0]+theta) },
		{ clog(coeft[1]), clog(coeft[1]-theta), clog(coeft[1]+theta) },
	};
	double complex xt[2];
	struct iapws_phi gamma_r;
	int i;

	/* gamma = g/RT */
	gamma->type = IAPWS_GAMMA;
	gamma->p = p;
	gamma->t = t;
	gamma->R = 1.0e-3;  /* J -> kJ */

	gamma->d00 = 0.0;
	gamma->d10 = 0.0;
	gamma->d01 = 0.0;
	gamma->d11 = 0.0;
	gamma->d20 = 0.0;
	gamma->d02 = 0.0;

	/* g0/RT */
	for (i = 0, pik = 1.0; i < ARRAY_SIZE(coef0); ++i, pik *= pi - p0) {
		gamma->d00 += coef0[i] * pik;
		if (i < ARRAY_SIZE(coef0) - 1)
			gamma->d10 += coef0[i + 1] * (i + 1) * pik;
		if (i < ARRAY_SIZE(coef0) - 2)
			gamma->d20 += coef0[i + 2] * (i + 2) * (i + 1) * pik;
	}
	gamma->d00 *= tau;
	gamma->d10 *= tau * pi;
	gamma->d20 *= tau * pi * pi;
	gamma->d01 = gamma->d00;
	gamma->d11 = gamma->d10;

	/* s0 */
	gamma->d00 -= s0;

	/* r2 */
	for (i = 0, pik = 1.0; i < ARRAY_SIZE(coefr); ++i, pik *= pi - p0) {
		r2[0] += coefr[i] * pik;
		if (i < ARRAY_SIZE(coefr) - 1)
			r2[1] += coefr[i + 1] * (i + 1) * pik;
		if (i < ARRAY_SIZE(coefr) - 2)
			r2[2] += coefr[i + 2] * (i + 2) * (i + 1) * pik;
	}
	r2[1] *= pi;
	r2[2] *= pi * pi;

	/* gr/RTt */
	xt[0] = -theta*theta / coeft[0] - coeft[0] * lt[0][0] * 2.0 +
		(coeft[0] - theta) * lt[0][1] + (coeft[0] + theta) * lt[0][2];
	xt[1] = -theta*theta / coeft[1] - coeft[1] * lt[1][0] * 2.0 +
		(coeft[1] - theta) * lt[1][1] + (coeft[1] + theta) * lt[1][2];

	gamma_r.d00 = creal(r1 * xt[0] + r2[0] * xt[1]);
	gamma_r.d10 = creal(r2[1] * xt[1]);
	gamma_r.d20 = creal(r2[2] * xt[1]);

	xt[0] = -theta / coeft[0] * 2.0 - lt[0][1] + lt[0][2];
	xt[1] = -theta / coeft[1] * 2.0 - lt[1][1] + lt[1][2];

	gamma_r.d01 = creal(r1 * xt[0] + r2[0] * xt[1]);
	gamma_r.d11 = creal(r2[1] * xt[1]);

	xt[0] = -2.0 / coeft[0] + 1.0 / (coeft[0] - theta)
		+ 1.0 / (coeft[0] + theta);
	xt[1] = -2.0 / coeft[1] + 1.0 / (coeft[1] - theta)
		+ 1.0 / (coeft[1] + theta);

	gamma_r.d02 = creal(r1 * xt[0] + r2[0] * xt[1]);

	/* gamma */
	gamma->d00 += gamma_r.d00 * tau;
	gamma->d10 += gamma_r.d10 * tau;
	gamma->d01 += gamma_r.d00 * tau - gamma_r.d01;
	gamma->d11 += gamma_r.d10 * tau - gamma_r.d11;
	gamma->d20 += gamma_r.d20 * tau;
	gamma->d02 += gamma_r.d02 * theta;
}

