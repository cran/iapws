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

#include <assert.h>
#include <math.h>

#include "iapws.h"
#include "phi.h"
#include "pow.h"

void iapws_phi(const double delta, const double tau,
		const struct iapws_phi_coef0 *coef0, const int n0,
		const struct iapws_phi_coef1 *coef1, const int n1,
		const struct iapws_phi_coef2 *coef2, const int n2,
		struct iapws_phi *phi)
{
	assert(n0 >= 2);

	int i;
	double xd, xt, xn, dc, edc, tc;

	phi->type = IAPWS_PHI;

	/* phi0 */
	phi->d00 = log(delta) + coef0[0].n + coef0[1].n * tau +
		coef0[2].n * log(tau);
	phi->d10 = 1.0;
	phi->d01 = coef0[1].n * tau + coef0[2].n;
	phi->d11 = 0.0;
	phi->d20 = -1.0;
	phi->d02 = -coef0[2].n;
	for (i = 3; i < n0; ++i) {
		const double gt = coef0[i].gamma * tau;
		const double egt = expm1(-gt);
		phi->d00 += coef0[i].n * log(-egt);
		phi->d01 -= coef0[i].n * gt * (1.0 + 1.0 / egt);
		phi->d02 -= coef0[i].n * POW2(gt) * (egt + 1) / POW2(egt);
	}

	/* phir */
	for (i = 0; i < n1; ++i) {
		if (i != 0) {
			dc *= POW(delta, coef1[i].c - coef1[i - 1].c);
			xd *= POW(delta, coef1[i].d - coef1[i - 1].d);
			xt *= POW(tau, coef1[i].t - coef1[i - 1].t);
			if (coef1[i].c != coef1[i - 1].c) edc = exp(-dc);
		} else {
			dc = POW(delta, coef1[i].c);
			xd = POW(delta, coef1[i].d);
			xt = POW(tau, coef1[i].t);
			edc = exp(-dc);
		}
		xn = coef1[i].n * edc * xd * xt;
		const double dcdc = coef1[i].d - coef1[i].c * dc;
		phi->d00 += xn;
		phi->d10 += xn * dcdc;
		phi->d01 += xn * coef1[i].t;
		phi->d11 += xn * dcdc * coef1[i].t;
		phi->d20 += xn * (dcdc * (dcdc - 1.0) - POW2(coef1[i].c) * dc);
		phi->d02 += xn * coef1[i].t * (coef1[i].t - 1.0);
	}
	for (i = 0; i < n2; ++i) {
		if (i != 0) {
			xd *= POW(delta, coef2[i].d - coef2[i - 1].d);
			xt *= POW(tau, coef2[i].t - coef2[i - 1].t);
		} else {
			xd = POW(delta, coef2[i].d);
			xt = POW(tau, coef2[i].t);
		}
		dc = delta - coef2[i].eps;
		tc = tau - coef2[i].gamma;
		xn = coef2[i].n * xd * xt *
			exp(-coef2[i].alpha * dc * dc -coef2[i].beta * tc * tc);
		dc = coef2[i].d - 2 * coef2[i].alpha * delta * dc;
		tc = coef2[i].t - 2 * coef2[i].beta * tau * tc;
		phi->d00 += xn;
		phi->d10 += xn * dc;
		phi->d01 += xn * tc;
		phi->d11 += xn * dc * tc;
		phi->d20 += xn * (POW2(dc) - coef2[i].d -
				2 * coef2[i].alpha * POW2(delta));
		phi->d02 += xn * (POW2(tc) - coef2[i].t -
				2 * coef2[i].beta * POW2(tau));
	}
}

