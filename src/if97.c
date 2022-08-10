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
 * IAPWS R7-97(2012), Revised Release on the IAPWS Industrial Formulation 1997
 * for the Thermodynamic Properties of Water and Steam (2007)
 */

#include <assert.h>
#include <math.h>

#include "if97.h"
#include "sat.h"
#include "nroot.h"

/* Static function forward declarations */
static void gamma_r1(double p, double t, iapws_phi *gamma);
static void gamma_r2(double p, double t, int meta, iapws_phi *gamma);
static void phi_r3(double delta, double tau, iapws_phi *phi);
static void get_delta_r3(double p, double tau, double *delta);
static void gamma_r5(double p, double t, iapws_phi *gamma);
static double pi_r4(double theta);
static double theta_r4(double pi);
static double pi_b23(double theta);
/*
static double theta_b23(double pi);
*/

typedef struct {
	int I;
	int J;
	double n;
} if97_coef;

/* Exported functions */

iapws_state_id if97_state(double p, double t)
{
	double ps;
	if (t < 273.15) return IAPWS_UNDEF;
	if (t >= IAPWS_TC && p >= IAPWS_PC) return IAPWS_CRIT;
	if (t >= IAPWS_TC) return IAPWS_GAS;
	ps = if97_psat(t);
	if (p > ps) return IAPWS_LIQUID;
	if (p < ps) return IAPWS_GAS;
	return IAPWS_UNDEF;
}

if97_region_id if97_region(double p, double t)
{
	double ps;
	if (t >= 273.15 && t <= 623.15) {
		ps = if97_psat(t);
		if (p > 0.0 && p <= ps) {
			return IF97_STEAM;
		} else if (p >= ps && p <= 100.0) {
			return IF97_WATER;
		}
	} else if (t >= 623.15 && t <= 863.15) {
		ps = pi_b23(t);
		if (p > 0.0 && p <= ps) {
			return IF97_STEAM;
		} else if (p >= ps && p <= 100.0) {
			return IF97_SUPER;
		}
	} else if (t >= 863.15 && t <= 1073.15) {
		if (p > 0.0 && p <= 100.0) {
			return IF97_STEAM;
		}
	} else if (t >= 1073.15 && t <= 2273.15) {
		if (p > 0.0 && p <= 50.0) {
			return IF97_GAS;
		}
	}
	return IF97_UNDEF;
}

double if97_tsat(double p)
{
	if (p >= 611.213e-6 && p <= IAPWS_PC) {
		return theta_r4(p);
	}
	return 0.0;
}

double if97_psat(double t)
{
	if (t >= 273.15 && t <= IAPWS_TC) {
		return pi_r4(t);
	}
	return 0.0;
}

int if97_gamma(double p, double t, iapws_state_id state, iapws_phi *gamma)
{
	if97_region_id reg = if97_region(p, t);
	int meta = 0;
	double tau, delta = 2.0;

	if (state == IAPWS_LIQUID) {
		if (reg == IF97_WATER) {
		} else if (reg == IF97_STEAM) {
			reg = IF97_WATER;
		} else if (reg == IF97_SUPER) {
			delta = sat_rhol(t) / IAPWS_RHOC * 1.3;
			if (delta == 0.0) delta = 2.0;
		} else {
			reg = IF97_UNDEF;
		}
	} else if (state == IAPWS_GAS) {
		if (reg == IF97_STEAM) {
		} else if (reg == IF97_WATER) {
			meta = (p < 10.0 ? 1 : 0);
			reg = IF97_STEAM;
		} else if (reg == IF97_SUPER) {
			delta = sat_rhog(t) / IAPWS_RHOC;
			if (delta == 0.0) delta = 0.5;
		} else if (reg == IF97_GAS) {
		} else {
			reg = IF97_UNDEF;
		}
	} else if (state == IAPWS_CRIT) {
		if (reg == IF97_STEAM) {
		} else if (reg == IF97_SUPER) {
		} else if (reg == IF97_GAS) {
		} else {
			reg = IF97_UNDEF;
		}
	} else {
		reg = IF97_UNDEF;
	}

	gamma->p = p;
	gamma->t = t;
	gamma->R = IF97_R;

	switch (reg) {
		case IF97_WATER:
			gamma_r1(p, t, gamma);
			break;
		case IF97_STEAM:
			gamma_r2(p, t, meta, gamma);
			break;
		case IF97_SUPER:
			tau = IAPWS_TC / t;
			get_delta_r3(p, tau, &delta);
			gamma->rho = delta * IAPWS_RHOC;
			phi_r3(delta, tau, gamma);
			break;
		case IF97_GAS:
			gamma_r5(p, t, gamma);
			break;
		default:
			return -1;
	}
	return 0;
}

/* static functions */

#define XI(i)	(XCOEF[i].I)
#define XJ(i)	(XCOEF[i].J)
#define XN(i)	(XCOEF[i].n)

#define XCOEF	(coef_r1)
#define XSIZE	(34)
static const if97_coef XCOEF[XSIZE] = {
	{  0,	 -2,	 0.14632971213167     },
	{  0,	 -1,	-0.84548187169114     },
	{  0,	  0,	-0.37563603672040e1   },
	{  0,	  1,	 0.33855169168385e1   },
	{  0,	  2,	-0.95791963387872     },
	{  0,	  3,	 0.15772038513228     },
	{  0,	  4,	-0.16616417199501e-1  },
	{  0,	  5,	 0.81214629983568e-3  },
	{  1,	 -9,	-0.28319080123804e-3  },
	{  1,	 -7,	 0.60706301565874e-3  },
	{  1,	 -1,	 0.18990068218419e-1  },
	{  1,	  0,	 0.32529748770505e-1  },
	{  1,	  1,	 0.21841717175414e-1  },
	{  1,	  3,	 0.52838357969930e-4  },
	{  2,	 -3,	-0.47184321073267e-3  },
	{  2,	  0,	-0.30001780793026e-3  },
	{  2,	  1,	 0.47661393906987e-4  },
	{  2,	  3,	-0.44141845330846e-5  },
	{  2,	 17,	-0.72694996297594e-15 },
	{  3,	 -4,	 0.31679644845054e-4  },
	{  3,	  0,	 0.28270797985312e-5  },
	{  3,	  6,	 0.85205128120103e-9  },
	{  4,	 -5,	-0.22425281908000e-5  },
	{  4,	 -2,	-0.65171222895601e-6  },
	{  4,	 10,	-0.14341729937924e-12 },
	{  5,	 -8,	 0.40516996860117e-6  },
	{  8,	-11,	-0.12734301741641e-8  },
	{  8,	 -6,	-0.17424871230634e-9  },
	{ 21,	-29,	 0.68762131295531e-18 },
	{ 23,	-31,	-0.14478307828521e-19 },
	{ 29,	-38,	-0.26335781662795e-22 },
	{ 30,	-39,	-0.11947622640071e-22 },
	{ 31,	-40,	-0.18228094581404e-23 },
	{ 32,	-41,	-0.93537087292458e-25 }
};
static void gamma_r1(double p, double t, iapws_phi *gamma)
{
	int i;
	double pi = p / 16.53;
	double tau = 1386.0 / t;
	double xp = pi - 7.1;
	double xt = tau - 1.222;
	double xxp = pi / xp;
	double xxt = tau / xt;
	double xn;
	gamma->type = IAPWS_GAMMA;
	gamma->d00 = 0.0;
	gamma->d10 = 0.0;
	gamma->d01 = 0.0;
	gamma->d11 = 0.0;
	gamma->d20 = 0.0;
	gamma->d02 = 0.0;
	for (i = 0; i < XSIZE; ++i) {
		xn = XN(i) * pow_di(xp, XI(i)) * pow_di(xt, XJ(i));
		gamma->d00 += xn;
		gamma->d10 += xn * XI(i);
		gamma->d01 += xn * XJ(i);
		gamma->d11 += xn * XI(i) * XJ(i);
		gamma->d20 += xn * XI(i) * (XI(i) - 1);
		gamma->d02 += xn * XJ(i) * (XJ(i) - 1);
	}
	gamma->d10 *= xxp;
	gamma->d01 *= xxt;
	gamma->d11 *= xxp * xxt;
	gamma->d20 *= xxp * xxp;
	gamma->d02 *= xxt * xxt;
}
#undef XSIZE
#undef XCOEF

#define XCOEF	(coef_r2)
#define XSIZE1	(9)
#define XSIZE2	(XSIZE1 + 43)
#define XSIZE3	(XSIZE2 + 9)
#define XSIZE4	(XSIZE3 + 13)
static const if97_coef XCOEF[XSIZE4] = {
	/* gamma_0 */
	{ 0,	 0,	-0.96927686500217e1  },
	{ 0,	 1,	 0.10086655968018e2  },
	{ 0,	-5,	-0.56087911283020e-2 },
	{ 0,	-4,	 0.71452738081455e-1 },
	{ 0,	-3,	-0.40710498223928    },
	{ 0,	-2,	 0.14240819171444e1  },
	{ 0,	-1,	-0.43839511319450e1  },
	{ 0,	 2,	-0.28408632460772    },
	{ 0,	 3,	 0.21268463753307e-1 },
	/* gamma_r */
	{  1,	 0,	-0.17731742473213e-2  },
	{  1,	 1,	-0.17834862292358e-1  },
	{  1,	 2,	-0.45996013696365e-1  },
	{  1,	 3,	-0.57581259083432e-1  },
	{  1,	 6,	-0.50325278727930e-1  },
	{  2,	 1,	-0.33032641670203e-4  },
	{  2,	 2,	-0.18948987516315e-3  },
	{  2,	 4,	-0.39392777243355e-2  },
	{  2,	 7,	-0.43797295650573e-1  },
	{  2,	36,	-0.26674547914087e-4  },
	{  3,	 0,	 0.20481737692309e-7  },
	{  3,	 1,	 0.43870667284435e-6  },
	{  3,	 3,	-0.32277677238570e-4  },
	{  3,	 6,	-0.15033924542148e-2  },
	{  3,	35,	-0.40668253562649e-1  },
	{  4,	 1,	-0.78847309559367e-9  },
	{  4,	 2,	 0.12790717852285e-7  },
	{  4,	 3,	 0.48225372718507e-6  },
	{  5,	 7,	 0.22922076337661e-5  },
	{  6,	 3,	-0.16714766451061e-10 },
	{  6,	16,	-0.21171472321355e-2  },
	{  6,	35,	-0.23895741934104e2   },
	{  7,	 0,	-0.59059564324270e-17 },
	{  7,	11,	-0.12621808899101e-5  },
	{  7,	25,	-0.38946842435739e-1  },
	{  8,	8 ,	 0.11256211360459e-10 },
	{  8,	36,	-0.82311340897998e1   },
	{  9,	13,	 0.19809712802088e-7  },
	{ 10,	 4,	 0.10406965210174e-18 },
	{ 10,	10,	-0.10234747095929e-12 },
	{ 10,	14,	-0.10018179379511e-8  },
	{ 16,	29,	-0.80882908646985e-10 },
	{ 16,	50,	 0.10693031879409     },
	{ 18,	57,	-0.33662250574171     },
	{ 20,	20,	 0.89185845355421e-24 },
	{ 20,	35,	 0.30629316876232e-12 },
	{ 20,	48,	-0.42002467698208e-5  },
	{ 21,	21,	-0.59056029685639e-25 },
	{ 22,	53,	 0.37826947613457e-5  },
	{ 23,	39,	-0.12768608934681e-14 },
	{ 24,	26,	 0.73087610595061e-28 },
	{ 24,	40,	 0.55414715350778e-16 },
	{ 24,	58,	-0.94369707241210e-6  },
	/* gamma_0m */
	{ 0,	 0,	-0.96937268393049e1  },
	{ 0,	 1,	 0.10087275970006e2  },
	{ 0,	-5,	-0.56087911283020e-2 },
	{ 0,	-4,	 0.71452738081455e-1 },
	{ 0,	-3,	-0.40710498223928    },
	{ 0,	-2,	 0.14240819171444e1  },
	{ 0,	-1,	-0.43839511319450e1  },
	{ 0,	 2,	-0.28408632460772    },
	{ 0,	 3,	 0.21268463753307e-1 },
	/* gamma_rm */
	{ 1,	 0,	-0.73362260186506e-2 },
	{ 1,	 2,	-0.88223831943146e-1 },
	{ 1,	 5,	-0.72334555213245e-1 },
	{ 1,	11,	-0.40813178534455e-2 },
	{ 2,	 1,	 0.20097803380207e-2 },
	{ 2,	 7,	-0.53045921898642e-1 },
	{ 2,	16,	-0.76190409086970e-2 },
	{ 3,	 4,	-0.63498037657313e-2 },
	{ 3,	16,	-0.86043093028588e-1 },
	{ 4,	 7,	 0.75321581522770e-2 },
	{ 4,	10,	-0.79238375446139e-2 },
	{ 5,	 9,	-0.22888160778447e-3 },
	{ 5,	10,	-0.26456501482810e-2 },
};
static void gamma_r2(double p, double t, int meta, iapws_phi *gamma)
{
	assert(meta == 0 || meta == 1);

	int i, n0, n1, n2;
	double pi = p;
	double tau = 540.0 / t;
	double xt = tau - 0.5;
	double xxt = tau / xt;
	double xn;

	if (meta == 0) {
		n0 = 0;
		n1 = XSIZE1;
		n2 = XSIZE2;
	} else {
		n0 = XSIZE2;
		n1 = XSIZE3;
		n2 = XSIZE4;
	}

	/* gamma_0 */
	gamma->type = IAPWS_GAMMA;
	gamma->d00 = log(pi);
	gamma->d10 = 1.0;
	gamma->d01 = 0.0;
	gamma->d11 = 0.0;
	gamma->d20 = -1.0;
	gamma->d02 = 0.0;
	for (i = n0; i < n1; ++i) {
		xn = XN(i) * pow_di(tau, XJ(i));
		gamma->d00 += xn;
		gamma->d01 += xn * XJ(i);
		gamma->d02 += xn * XJ(i) * (XJ(i) - 1);
	}

	/* gamma_r */
	iapws_phi gamma_r = { IAPWS_GAMMA };
	for (i = n1; i < n2; ++i) {
		xn = XN(i) * pow_di(pi, XI(i)) * pow_di(xt, XJ(i));
		gamma_r.d00 += xn;
		gamma_r.d10 += xn * XI(i);
		gamma_r.d01 += xn * XJ(i);
		gamma_r.d11 += xn * XI(i) * XJ(i);
		gamma_r.d20 += xn * XI(i) * (XI(i) - 1);
		gamma_r.d02 += xn * XJ(i) * (XJ(i) - 1);
	}
	gamma_r.d01 *= xxt;
	gamma_r.d11 *= xxt;
	gamma_r.d02 *= xxt * xxt;

	/* gamma */
	gamma->d00 += gamma_r.d00;
	gamma->d10 += gamma_r.d10;
	gamma->d01 += gamma_r.d01;
	gamma->d11 += gamma_r.d11;
	gamma->d20 += gamma_r.d20;
	gamma->d02 += gamma_r.d02;
}
#undef XSIZE4
#undef XSIZE3
#undef XSIZE2
#undef XSIZE1
#undef XCOEF

#define XCOEF	(coef_r3)
#define XSIZE	(40)
static const if97_coef XCOEF[XSIZE] = {
	{ 0,	0,	 0.10658070028513e1  },
	{ 0,	0,	-0.15732845290239e2  },
	{ 0,	1,	 0.20944396974307e2  },
	{ 0,	2,	-0.76867707878716e1  },
	{ 0,	7,	 0.26185947787954e1  },
	{ 0,	10,	-0.28080781148620e1  },
	{ 0,	12,	 0.12053369696517e1  },
	{ 0,	23,	-0.84566812812502e-2 },
	{ 1,	2,	-0.12654315477714e1  },
	{ 1,	6,	-0.11524407806681e1  },
	{ 1,	15,	 0.88521043984318    },
	{ 1,	17,	-0.64207765181607    },
	{ 2,	0,	 0.38493460186671    },
	{ 2,	2,	-0.85214708824206    },
	{ 2,	6,	 0.48972281541877e1  },
	{ 2,	7,	-0.30502617256965e1  },
	{ 2,	22,	 0.39420536879154e-1 },
	{ 2,	26,	 0.12558408424308    },
	{ 3,	0,	-0.27999329698710    },
	{ 3,	2,	 0.13899799569460e1  },
	{ 3,	4,	-0.20189915023570e1  },
	{ 3,	16,	-0.82147637173963e-2 },
	{ 3,	26,	-0.47596035734923    },
	{ 4,	0,	 0.43984074473500e-1 },
	{ 4,	2,	-0.44476435428739    },
	{ 4,	4,	 0.90572070719733    },
	{ 4,	26,	 0.70522450087967    },
	{ 5,	1,	 0.10770512626332    },
	{ 5,	3,	-0.32913623258954    },
	{ 5,	26,	-0.50871062041158    },
	{ 6,	0,	-0.22175400873096e-1 },
	{ 6,	2,	 0.94260751665092e-1 },
	{ 6,	26,	 0.16436278447961    },
	{ 7,	2,	-0.13503372241348e-1 },
	{ 8,	26,	-0.14834345352472e-1 },
	{ 9,	2,	 0.57922953628084e-3 },
	{ 9,	26,	 0.32308904703711e-2 },
	{ 10,	0,	 0.80964802996215e-4 },
	{ 10,	1,	-0.16557679795037e-3 },
	{ 11,	26,	-0.44923899061815e-4 },
};
static void phi_r3(double delta, double tau, iapws_phi *phi)
{
	int i;
	double xn;

	phi->type = IAPWS_PHI;
	phi->d00 = XN(0) * log(delta);
	phi->d10 = XN(0);
	phi->d01 = 0.0;
	phi->d11 = 0.0;
	phi->d20 = -XN(0);
	phi->d02 = 0.0;

	for (i = 1; i < XSIZE; ++i) {
		xn = XN(i) * pow_di(delta, XI(i)) * pow_di(tau, XJ(i));
		phi->d00 += xn;
		phi->d10 += xn * XI(i);
		phi->d01 += xn * XJ(i);
		phi->d11 += xn * XI(i) * XJ(i);
		phi->d20 += xn * XI(i) * (XI(i) - 1);
		phi->d02 += xn * XJ(i) * (XJ(i) - 1);
	}
}
#undef XCOEF
#undef XSIZE

static void get_pdp_r3(double *delta, void *ptau, double *p, double *dp)
{
	double *dptau = (double *)(ptau);
	double p0 = dptau[0];
	double tau = dptau[1];
	iapws_phi phi;
	phi_r3(*delta, tau, &phi);
	double fact = (IF97_R * IAPWS_TC * IAPWS_RHOC * 1e-3) / tau;
	*p = phi.d10 * (*delta) * fact - p0;
	*dp = (phi.d10 * 2. + phi.d20) * fact;
}

static void get_delta_r3(double p, double tau, double *delta)
{
	int iter;
	int maxiter = iapws_nroot_max_iter;
	double tol = iapws_nroot_tol;
	double ptau[2] = { p, tau };
	double err;
	iter = nroot(get_pdp_r3, delta, ptau, &err, tol, maxiter);
	assert(iter > 0 && iter < maxiter);  /* FIXME */
}

#define XCOEF	(coef_r5)
#define XSIZE1	(6)
#define XSIZE2	(XSIZE1 + 6)
static const if97_coef XCOEF[XSIZE2] = {
	/* gamma_0 */
	{ 0,	 0,	-0.13179983674201e2  },
	{ 0,	 1,	 0.68540841634434e1  },
	{ 0,	-3,	-0.24805148933466e-1 },
	{ 0,	-2,	 0.36901534980333    },
	{ 0,	-1,	-0.31161318213925e1  },
	{ 0,	 2,	-0.32961626538917    },
	/* gamma_r */
	{1,	1,	 0.15736404855259e-2 },
	{1,	2,	 0.90153761673944e-3 },
	{1,	3,	-0.50270077677648e-2 },
	{2,	3,	 0.22440037409485e-5 },
	{2,	9,	-0.41163275453471e-5 },
	{3,	7,	 0.37919454822955e-7 },
};
static void gamma_r5(double p, double t, iapws_phi *gamma)
{
	int i;
	double pi = p;
	double tau = 1000.0 / t;
	double xn;

	/* gamma_0 */
	gamma->type = IAPWS_GAMMA;
	gamma->d00 = log(pi);
	gamma->d10 = 1.0;
	gamma->d01 = 0.0;
	gamma->d11 = 0.0;
	gamma->d20 = -1.0;
	gamma->d02 = 0.0;
	for (i = 0; i < XSIZE1; ++i) {
		xn = XN(i) * pow_di(tau, XJ(i));
		gamma->d00 += xn;
		gamma->d01 += xn * XJ(i);
		gamma->d02 += xn * XJ(i) * (XJ(i) - 1);
	}

	/* gamma_r */
	for (i = XSIZE1; i < XSIZE2; ++i) {
		xn = XN(i) * pow_di(pi, XI(i)) * pow_di(tau, XJ(i));
		gamma->d00 += xn;
		gamma->d10 += xn * XI(i);
		gamma->d01 += xn * XJ(i);
		gamma->d11 += xn * XI(i) * XJ(i);
		gamma->d20 += xn * XI(i) * (XI(i) - 1);
		gamma->d02 += xn * XJ(i) * (XJ(i) - 1);
	}
}
#undef XSIZE2
#undef XSIZE1
#undef XCOEF

/* boundaries */

#define XSIZE	(10)
#define XCOEF	(coef_r4)
static double XCOEF[XSIZE] = {
	0.11670521452767e4, -0.72421316703206e6,
	-0.17073846940092e2, 0.12020824702470e5,
	-0.32325550322333e7, 0.14915108613530e2,
	-0.48232657361591e4, 0.40511340542057e6,
	-0.23855557567849, 0.65017534844798e3,
};

static double pi_r4(double theta)
{
	theta += XCOEF[8] / (theta - XCOEF[9]);
	double theta2 = theta * theta;
	double a = theta2 + XCOEF[0] * theta + XCOEF[1];
	double b = XCOEF[2] * theta2 + XCOEF[3] * theta + XCOEF[4];
	double c = XCOEF[5] * theta2 + XCOEF[6] * theta + XCOEF[7];
	double ans = 2.0 * c / (-b + sqrt(b * b - 4.0 * a * c));
	ans *= ans;
	return ans * ans;
}

static double theta_r4(double pi)
{
	double beta = pow(pi, .25);
	double beta2 = beta * beta;
	double e = beta2 + XCOEF[2] * beta + XCOEF[5];
	double f = XCOEF[0] * beta2 + XCOEF[3] * beta + XCOEF[6];
	double g = XCOEF[1] * beta2 + XCOEF[4] * beta + XCOEF[7];
	double d = -2.0 * g / (f + sqrt(f * f - 4.0 * e * g));
	return (XCOEF[9] + d - sqrt((XCOEF[9] + d) * (XCOEF[9] + d) -
				4.0 * (XCOEF[8] + XCOEF[9] * d))) * .5;
}
#undef XSIZE
#undef XCOEF

#define XSIZE	5
#define XCOEF	(coef_b23)
static double XCOEF[XSIZE] = {
	.34805185628969e+3, -.11671859879975e+1,
	.10192970039326e-2, .57254459862746e+3,
	.13918839778870e+2,
};

static double pi_b23(double theta)
{
	return XCOEF[0] + XCOEF[1] * theta + XCOEF[2] * theta * theta;
}

/*
static double theta_b23(double pi)
{
	return XCOEF[3] + sqrt((pi - XCOEF[4]) / XCOEF[3]);
}
*/
#undef XSIZE
#undef XCOEF

