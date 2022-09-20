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
#include "melt.h"
#include "nroot.h"

/* Static function forward declarations */
static void gamma_r1(double p, double t, iapws_phi *gamma);
static int gamma_r1_ph(double p, double h, iapws_phi *gamma);
static double t_r1_ph(double p, double h);

static void gamma_r2(double p, double t, int meta, iapws_phi *gamma);
static int gamma_r2_ph(double p, double h, iapws_phi *gamma);
static double pi_b2bc(double eta);
static double t_r2_ph(double p, double h);

static void phi_r3(double rho, double t, iapws_phi *phi);
static int phi_r3_pt(double p, double t, iapws_phi *phi);
static int phi_r3_ph(double p, double h, iapws_phi *phi);
static double eta_b3ab(double pi);
static double t_r3_ph(double p, double h);
static double v_r3_ph(double p, double h);

static double pi_r4(double theta);
static double theta_r4(double pi);

static void gamma_r5(double p, double t, iapws_phi *gamma);

static double pi_b23(double theta);
static double theta_b23(double pi);

typedef struct {
	int I;
	int J;
	double n;
} if97_coef;

/* Exported functions */

iapws_state_id if97_state(double p, double t)
{
	double ps;
	if (t >= 273.16 && t < IAPWS_TC && p < 620.0) {
		ps = if97_psat(t);
		if (p > ps) return IAPWS_LIQUID;
		if (p < ps) return IAPWS_GAS;
		return IAPWS_SAT;
	} else if (t >= IAPWS_TC) {
		if (p < IAPWS_PC) return IAPWS_GAS;
		else return IAPWS_CRIT;
	}
	return melt_sub_state(p, t);
}

if97_region_id if97_region(double p, double t)
{
	double ps;
	if (t >= 273.15 && t <= 623.15) {
		ps = if97_psat(t);
		if (p > 0.0 && p <= ps) {
			return IF97_REGION_2;
		} else if (p >= ps && p <= 100.0) {
			return IF97_REGION_1;
		}
	} else if (t >= 623.15 && t <= 863.15) {
		ps = pi_b23(t);
		if (p > 0.0 && p <= ps) {
			return IF97_REGION_2;
		} else if (p >= ps && p <= 100.0) {
			return IF97_REGION_3;
		}
	} else if (t >= 863.15 && t <= 1073.15) {
		if (p > 0.0 && p <= 100.0) {
			return IF97_REGION_2;
		}
	} else if (t >= 1073.15 && t <= 2273.15) {
		if (p > 0.0 && p <= 50.0) {
			return IF97_REGION_5;
		}
	}
	return IF97_REGION_UNDEF;
}

if97_region_id if97_region_ph(double p, double h)
{
	int err;
	double ts;
	iapws_phi gamma;
	gamma.R = IF97_R;

	if (p <= pi_r4(273.15)) {
		gamma_r2(p, 273.15, 0, &gamma);
		if (h >= iapws_h(&gamma)) return IF97_REGION_2;
		else return IF97_REGION_UNDEF;
	} else if (p <= pi_r4(623.15)) {
		gamma_r1(p, 273.15, &gamma);
		if (h < iapws_h(&gamma)) return IF97_REGION_UNDEF;
		ts = theta_r4(p);
		gamma_r1(p, ts, &gamma);
		if (h <= iapws_h(&gamma)) return IF97_REGION_1;
		gamma_r2(p, ts, 0, &gamma);
		if (h >= iapws_h(&gamma)) return IF97_REGION_2;
		return IF97_REGION_4;
	} else if (p <= 100.0) {
		gamma_r1(p, 273.15, &gamma);
		if (h < iapws_h(&gamma)) return IF97_REGION_UNDEF;
		gamma_r1(p, 623.15, &gamma);
		if (h < iapws_h(&gamma)) return IF97_REGION_1;
		gamma_r2(p, theta_b23(p), 0, &gamma);
		if (h > iapws_h(&gamma)) return IF97_REGION_2;
		if (p < IAPWS_PC) {
			ts = theta_r4(p);
			if (h <= 2087.546845117154) { /* HC */
				gamma.rho = 650.0;
				if ((err = phi_r3_pt(p, ts, &gamma)) != 0)
					return IF97_REGION_UNDEF;
				if (h > iapws_h(&gamma)) return IF97_REGION_4;
				else return IF97_REGION_3;
			} else {
				gamma.rho = 150.0;
				if ((err = phi_r3_pt(p, ts, &gamma)) != 0)
					return IF97_REGION_UNDEF;
				if (h < iapws_h(&gamma)) return IF97_REGION_4;
				else return IF97_REGION_3;
			}
		}
		return IF97_REGION_3;
	}
	return IF97_REGION_UNDEF;
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
	int err = 0;
	double rho0 = 650.0;

	if (state == IAPWS_LIQUID) {
		if (reg == IF97_REGION_1) {
		} else if (reg == IF97_REGION_2) {
			reg = IF97_REGION_1;
		} else if (reg == IF97_REGION_3) {
		} else {
			reg = IF97_REGION_UNDEF;
		}
	} else if (state == IAPWS_GAS) {
		if (reg == IF97_REGION_2) {
		} else if (reg == IF97_REGION_1) {
			meta = (p < 10.0 ? 1 : 0);
			reg = IF97_REGION_2;
		} else if (reg == IF97_REGION_3) {
			rho0 = 150.0;
		} else if (reg == IF97_REGION_5) {
		} else {
			reg = IF97_REGION_UNDEF;
		}
	} else if (state == IAPWS_CRIT) {
		if (reg == IF97_REGION_3) {
		} else if (reg == IF97_REGION_2) {
		} else if (reg == IF97_REGION_5) {
		} else {
			reg = IF97_REGION_UNDEF;
		}
	} else {
		reg = IF97_REGION_UNDEF;
	}

	gamma->R = IF97_R;

	switch (reg) {
		case IF97_REGION_1:
			gamma_r1(p, t, gamma);
			break;
		case IF97_REGION_2:
			gamma_r2(p, t, meta, gamma);
			break;
		case IF97_REGION_3:
			gamma->rho = rho0;  /* initial guess */
			err = phi_r3_pt(p, t, gamma);
			break;
		case IF97_REGION_5:
			gamma_r5(p, t, gamma);
			break;
		default:
			err = -1;
	}
	return err;
}

int if97_gamma_ph(double p, double h, iapws_phi *gamma)
{
	int err;
	if97_region_id reg = if97_region_ph(p, h);

	gamma->R = IF97_R;

	switch (reg) {
		case IF97_REGION_1:
			err = gamma_r1_ph(p, h, gamma);
			break;
		case IF97_REGION_2:
			err = gamma_r2_ph(p, h, gamma);
			break;
		case IF97_REGION_3:
			err = phi_r3_ph(p, h, gamma);
		default:
			err = -1;
	}
	return err;
}

/* Region 1 */

static void gamma_r1(double p, double t, iapws_phi *gamma)
{
	enum { SIZE = 34 };
	const if97_coef coef[SIZE] = {
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
	const double pi = p / 16.53;
	const double tau = 1386.0 / t;
	const double xp = pi - 7.1;
	const double xt = tau - 1.222;
	const double xxp = pi / xp;
	const double xxt = tau / xt;

	int i;
	double xn;

	gamma->type = IAPWS_GAMMA;
	gamma->d00 = 0.0;
	gamma->d10 = 0.0;
	gamma->d01 = 0.0;
	gamma->d11 = 0.0;
	gamma->d20 = 0.0;
	gamma->d02 = 0.0;
	gamma->p = p;
	gamma->t = t;
	for (i = 0; i < SIZE; ++i) {
		xn = coef[i].n * POWINT(xp, coef[i].I) * POWINT(xt, coef[i].J);
		gamma->d00 += xn;
		gamma->d10 += xn * coef[i].I;
		gamma->d01 += xn * coef[i].J;
		gamma->d11 += xn * coef[i].I * coef[i].J;
		gamma->d20 += xn * coef[i].I * (coef[i].I - 1);
		gamma->d02 += xn * coef[i].J * (coef[i].J - 1);
	}
	gamma->d10 *= xxp;
	gamma->d01 *= xxt;
	gamma->d11 *= xxp * xxt;
	gamma->d20 *= xxp * xxp;
	gamma->d02 *= xxt * xxt;
}

static double t_r1_ph(double p, double h)
{
	enum { SIZE = 20 };
	const if97_coef coef[SIZE] = {
		{ 0,	 0,	-0.23872489924521e+3  },
		{ 0,	 1,	 0.40421188637945e+3  },
		{ 0,	 2,	 0.11349746881718e+3  },
		{ 0,	 6,	-0.58457616048039e+1  },
		{ 0,	22,	-0.15285482413140e-3  },
		{ 0,	32,	-0.10866707695377e-5  },
		{ 1,	 0,	-0.13391744872602e+2  },
		{ 1,	 1,	 0.43211039183559e+2  },
		{ 1,	 2,	-0.54010067170506e+2  },
		{ 1,	 3,	 0.30535892203916e+2  },
		{ 1,	 4,	-0.65964749423638e+1  },
		{ 1,	10,	 0.93965400878363e-2  },
		{ 1,	32,	 0.11573647505340e-6  },
		{ 2,	10,	-0.25858641282073e-4  },
		{ 2,	32,	-0.40644363084799e-8  },
		{ 3,	10,	 0.66456186191635e-7  },
		{ 3,	32,	 0.80670734103027e-10 },
		{ 4,	32,	-0.93477771213947e-12 },
		{ 5,	32,	 0.58265442020601e-14 },
		{ 6,	32,	-0.15020185953503e-16 },
	};

	int i;
	double ans = 0.0;

	h = h * 4.0e-4 + 1.0;
	for (i = 0; i < SIZE; ++i) {
		ans += coef[i].n * POWINT(p, coef[i].I) * POWINT(h, coef[i].J);
	}
	return ans;
}

static void get_gamma_r1_ph(double *t, void *xphi, double *h, double *dh)
{
	iapws_phi *gamma = (iapws_phi *)(xphi);
	gamma_r1(gamma->p, *t, gamma);
	*h = gamma->d01 * gamma->R * gamma->t - gamma->h;
	*dh = -gamma->d02 * gamma->R;
}

static int gamma_r1_ph(double p, double h, iapws_phi *gamma)
{
	int maxiter = nroot_maxiter;
	double tolf = nroot_tolf;
	double tolx = nroot_tolx;
	gamma->p = p;
	gamma->h = h;
	gamma->t = t_r1_ph(p, h);
	return nroot1(get_gamma_r1_ph, &gamma->t, gamma,
			&tolf, &tolx, &maxiter, nroot_verbose);
}

/* Region 2 */

static void gamma_r2(double p, double t, int meta, iapws_phi *gamma)
{
	assert(meta == 0 || meta == 1);

	enum {
		SIZE1 = 9,
		SIZE2 = SIZE1 + 43,
		SIZE3 = SIZE2 + 9,
		SIZE4 = SIZE3 + 13,
	};
	const if97_coef coef[SIZE4] = {
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

	const double pi = p;
	const double tau = 540.0 / t;
	const double xt = tau - 0.5;
	const double xxt = tau / xt;

	int i, n0, n1, n2;
	double xn;

	if (meta == 0) {
		n0 = 0;
		n1 = SIZE1;
		n2 = SIZE2;
	} else {
		n0 = SIZE2;
		n1 = SIZE3;
		n2 = SIZE4;
	}

	/* gamma_0 */
	gamma->type = IAPWS_GAMMA;
	gamma->d00 = log(pi);
	gamma->d10 = 1.0;
	gamma->d01 = 0.0;
	gamma->d11 = 0.0;
	gamma->d20 = -1.0;
	gamma->d02 = 0.0;
	gamma->p = p;
	gamma->t = t;
	for (i = n0; i < n1; ++i) {
		xn = coef[i].n * POWINT(tau, coef[i].J);
		gamma->d00 += xn;
		gamma->d01 += xn * coef[i].J;
		gamma->d02 += xn * coef[i].J * (coef[i].J - 1);
	}

	/* gamma_r */
	iapws_phi gamma_r = { IAPWS_GAMMA };
	for (i = n1; i < n2; ++i) {
		xn = coef[i].n * POWINT(pi, coef[i].I) * POWINT(xt, coef[i].J);
		gamma_r.d00 += xn;
		gamma_r.d10 += xn * coef[i].I;
		gamma_r.d01 += xn * coef[i].J;
		gamma_r.d11 += xn * coef[i].I * coef[i].J;
		gamma_r.d20 += xn * coef[i].I * (coef[i].I - 1);
		gamma_r.d02 += xn * coef[i].J * (coef[i].J - 1);
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

static double t_r2_ph(double p, double h)
{
	enum {
		SIZE1 = 34,
		SIZE2 = SIZE1 + 38,
		SIZE3 = SIZE2 + 23,
       	};
	const if97_coef coef[SIZE3] = {
		/* region 2a */
		{ 0,	 0,	 0.10898952318288e4	},
		{ 0,	 1,	 0.84951654495535e3	},
		{ 0,	 2,	-0.10781748091826e3	},
		{ 0,	 3,	 0.33153654801263e2	},
		{ 0,	 7,	-0.74232016790248e1	},
		{ 0,	20,	 0.11765048724356e2	},
		{ 1,	 0,	 0.18445749355790e1	},
		{ 1,	 1,	-0.41792700549624e1	},
		{ 1,	 2,	 0.62478196935812e1	},
		{ 1,	 3,	-0.17344563108114e2	},
		{ 1,	 7,	-0.20058176862096e3	},
		{ 1,	 9,	 0.27196065473796e3	},
		{ 1,	11,	-0.45511318285818e3	},
		{ 1,	18,	 0.30919688604755e4	},
		{ 1,	44,	 0.25226640357872e6	},
		{ 2,	 0,	-0.61707422868339e-2	},
		{ 2,	 2,	-0.31078046629583	},
		{ 2,	 7,	 0.11670873077107e2	},
		{ 2,	36,	 0.12812798404046e9	},
		{ 2,	38,	-0.98554909623276e9	},
		{ 2,	40,	 0.28224546973002e10	},
		{ 2,	42,	-0.35948971410703e10	},
		{ 2,	44,	 0.17227349913197e10	},
		{ 3,	24,	-0.13551334240775e5	},
		{ 3,	44,	 0.12848734664650e8	},
		{ 4,	12,	 0.13865724283226e1	},
		{ 4,	32,	 0.23598832556514e6	},
		{ 4,	44,	-0.13105236545054e8	},
		{ 5,	32,	 0.73999835474766e4	},
		{ 5,	36,	-0.55196697030060e6	},
		{ 5,	42,	 0.37154085996233e7	},
		{ 6,	34,	 0.19127729239660e5	},
		{ 6,	44,	-0.41535164835634e6	},
		{ 7,	28,	-0.62459855192507e2	},
		/* region 2b */
		{ 0,	 0,	 0.14895041079516e4	},
		{ 0,	 1,	 0.74307798314034e3	},
		{ 0,	 2,	-0.97708318797837e2	},
		{ 0,	12,	 0.24742464705674e1	},
		{ 0,	18,	-0.63281320016026	},
		{ 0,	24,	 0.11385952129658e1	},
		{ 0,	28,	-0.47811863648625	},
		{ 0,	40,	 0.85208123431544e-2	},
		{ 1,	 0,	 0.93747147377932	},
		{ 1,	 2,	 0.33593118604916e1	},
		{ 1,	 6,	 0.33809355601454e1	},
		{ 1,	12,	 0.16844539671904	},
		{ 1,	18,	 0.73875745236695	},
		{ 1,	24,	-0.47128737436186	},
		{ 1,	28,	 0.15020273139707	},
		{ 1,	40,	-0.21764114219750e-2	},
		{ 2,	 2,	-0.21810755324761e-1	},
		{ 2,	 8,	-0.10829784403677	},
		{ 2,	18,	-0.46333324635812e-1	},
		{ 2,	40,	 0.71280351959551e-4	},
		{ 3,	 1,	 0.11032831789999e-3	},
		{ 3,	 2,	 0.18955248387902e-3	},
		{ 3,	12,	 0.30891541160537e-2	},
		{ 3,	24,	 0.13555504554949e-2	},
		{ 4,	 2,	 0.28640237477456e-6	},
		{ 4,	12,	-0.10779857357512e-4	},
		{ 4,	18,	-0.76462712454814e-4	},
		{ 4,	24,	 0.14052392818316e-4	},
		{ 4,	28,	-0.31083814331434e-4	},
		{ 4,	40,	-0.10302738212103e-5	},
		{ 5,	18,	 0.28217281635040e-6	},
		{ 5,	24,	 0.12704902271945e-5	},
		{ 5,	40,	 0.73803353468292e-7	},
		{ 6,	28,	-0.11030139238909e-7	},
		{ 7,	 2,	-0.81456365207833e-13	},
		{ 7,	28,	-0.25180545682962e-10	},
		{ 9,	 1,	-0.17565233969407e-17	},
		{ 9,	40,	 0.86934156344163e-14	},
		/* region 2c */
		{-7,	 0,	-0.32368398555242e13	},
		{-7,	 4,	 0.73263350902181e13	},
		{-6,	 0,	 0.35825089945447e12	},
		{-6,	 2,	-0.58340131851590e12	},
		{-5,	 0,	-0.10783068217470e11	},
		{-5,	 2,	 0.20825544563171e11	},
		{-2,	 0,	 0.61074783564516e6	},
		{-2,	 1,	 0.85977722535580e6	},
		{-1,	 0,	-0.25745723604170e5	},
		{-1,	 2,	 0.31081088422714e5	},
		{ 0,	 0,	 0.12082315865936e4	},
		{ 0,	 1,	 0.48219755109255e3	},
		{ 1,	 4,	 0.37966001272486e1	},
		{ 1,	 8,	-0.10842984880077e2	},
		{ 2,	 4,	-0.45364172676660e-1	},
		{ 6,	 0,	 0.14559115658698e-12	},
		{ 6,	 1,	 0.11261597407230e-11	},
		{ 6,	 4,	-0.17804982240686e-10	},
		{ 6,	10,	 0.12324579690832e-6	},
		{ 6,	12,	-0.11606921130984e-5	},
		{ 6,	16,	 0.27846367088554e-4	},
		{ 6,	20,	-0.59270038474176e-3	},
		{ 6,	22,	 0.12918582991878e-2	},
	};

	int i, n;
	double ans = 0.0;

	if (p <= 4.0) {
		i = 0;
		n = SIZE1;
		h = h * 5.0e-4 - 2.1;
	} else if (p < pi_b2bc(h)) {
		i = SIZE1;
		n = SIZE2;
		p -= 2.0;
		h = h * 5.0e-4 - 2.6;
	} else {
		i = SIZE2;
		n = SIZE3;
		p += 25.0;
		h = h * 5.0e-4 - 1.8;
	}

	for (; i < n; ++i) {
		ans += coef[i].n * POWINT(p, coef[i].I) * POWINT(h, coef[i].J);
	}
	return ans;
}

static void get_gamma_r2_ph(double *t, void *xphi, double *h, double *dh)
{
	iapws_phi *gamma = (iapws_phi *)(xphi);
	gamma_r2(gamma->p, *t, 0, gamma);
	*h = gamma->d01 * gamma->R * gamma->t - gamma->h;
	*dh = -gamma->d02 * gamma->R;
}

static int gamma_r2_ph(double p, double h, iapws_phi *gamma)
{
	int maxiter = nroot_maxiter;
	double tolf = nroot_tolf;
	double tolx = nroot_tolx;
	gamma->p = p;
	gamma->h = h;
	gamma->t = t_r2_ph(p, h);
	return nroot1(get_gamma_r2_ph, &gamma->t, gamma,
			&tolf, &tolx, &maxiter, nroot_verbose);
}


/* Region 3 */

static void phi_r3(double rho, double t, iapws_phi *phi)
{
	enum { SIZE = 40 };
	const if97_coef coef[SIZE] = {
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

	const double delta = rho / IAPWS_RHOC;
	const double tau = IAPWS_TC / t;

	int i;
	double xn;

	phi->type = IAPWS_PHI;
	phi->d00 = coef[0].n * log(delta);
	phi->d10 = coef[0].n;
	phi->d01 = 0.0;
	phi->d11 = 0.0;
	phi->d20 = -coef[0].n;
	phi->d02 = 0.0;
	phi->rho = rho;
	phi->t = t;

	for (i = 1; i < SIZE; ++i) {
		xn = coef[i].n * POWINT(delta, coef[i].I) * POWINT(tau, coef[i].J);
		phi->d00 += xn;
		phi->d10 += xn * coef[i].I;
		phi->d01 += xn * coef[i].J;
		phi->d11 += xn * coef[i].I * coef[i].J;
		phi->d20 += xn * coef[i].I * (coef[i].I - 1);
		phi->d02 += xn * coef[i].J * (coef[i].J - 1);
	}
}

static void get_phi_r3_pt(double *rho, void *xphi, double *p, double *dp)
{
	iapws_phi *phi = (iapws_phi *)(xphi);
	phi_r3(*rho, phi->t, phi);
	double fact = phi->R * phi->t * 1e-3;
	*p = phi->d10 * (*rho) * fact - phi->p;
	*dp = (phi->d10 * 2.0 + phi->d20) * fact;
}

static int phi_r3_pt(double p, double t, iapws_phi *phi)
{
	int maxiter = nroot_maxiter;
	double tolf = nroot_tolf;
	double tolx = nroot_tolx;
	phi->p = p;
	phi->t = t;
	return nroot1(get_phi_r3_pt, &phi->rho, phi,
			&tolf, &tolx, &maxiter, nroot_verbose);
}

static double t_r3_ph(double p, double h)
{
	enum {
		SIZE1 = 31,
		SIZE2 = SIZE1 + 33,
       	};
	const if97_coef coef[SIZE2] = {
		/* region 3a */
		{ -12,	0,	-0.133645667811215e-6	},
		{ -12,	1,	 0.455912656802978e-5	},
		{ -12,	2,	-0.146294640700979e-4	},
		{ -12,	6,	 0.639341312970080e-2	},
		{ -12,	14,	 0.372783927268847e3	},
		{ -12,	16,	-0.718654377460447e4	},
		{ -12,	20,	 0.573494752103400e6	},
		{ -12,	22,	-0.267569329111439e7	},
		{ -10,	1,	-0.334066283302614e-4	},
		{ -10,	5,	-0.245479214069597e-1	},
		{ -10,	12,	 0.478087847764996e2	},
		{ -8,	0,	 0.764664131818904e-5	},
		{ -8,	2,	 0.128350627676972e-2	},
		{ -8,	4,	 0.171219081377331e-1	},
		{ -8,	10,	-0.851007304583213e1	},
		{ -5,	2,	-0.136513461629781e-1	},
		{ -3,	0,	-0.384460997596657e-5	},
		{ -2,	1,	 0.337423807911655e-2	},
		{ -2,	3,	-0.551624873066791	},
		{ -2,	4,	 0.729202277107470	},
		{ -1,	0,	-0.992522757376041e-2	},
		{ -1,	2,	-0.119308831407288	},
		{ 0,	0,	 0.793929190615421	},
		{ 0,	1,	 0.454270731799386	},
		{ 1,	1,	 0.209998591259910	},
		{ 3,	0,	-0.642109823904738e-2	},
		{ 3,	1,	-0.235155868604540e-1	},
		{ 4,	0,	 0.252233108341612e-2	},
		{ 4,	3,	-0.764885133368119e-2	},
		{ 10,	4,	 0.136176427574291e-1	},
		{ 12,	5,	-0.133027883575669e-1	},
		/* region 3b */
		{ -12,	0,	 0.323254573644920e-4	},
		{ -12,	1,	-0.127575556587181e-3	},
		{ -10,	0,	-0.475851877356068e-3	},
		{ -10,	1,	 0.156183014181602e-2	},
		{ -10,	5,	 0.105724860113781	},
		{ -10,	10,	-0.858514221132534e2	},
		{ -10,	12,	 0.724140095480911e3	},
		{ -8,	0,	 0.296475810273257e-2	},
		{ -8,	1,	-0.592721983365988e-2	},
		{ -8,	2,	-0.126305422818666e-1	},
		{ -8,	4,	-0.115716196364853	},
		{ -8,	10,	 0.849000969739595e2	},
		{ -6,	0,	-0.108602260086615e-1	},
		{ -6,	1,	 0.154304475328851e-1	},
		{ -6,	2,	 0.750455441524466e-1	},
		{ -4,	0,	 0.252520973612982e-1	},
		{ -4,	1,	-0.602507901232996e-1	},
		{ -3,	5,	-0.307622221350501e1	},
		{ -2,	0,	-0.574011959864879e-1	},
		{ -2,	4,	 0.503471360939849e1	},
		{ -1,	2,	-0.925081888584834	},
		{ -1,	4,	 0.391733882917546e1	},
		{ -1,	6,	-0.773146007130190e2	},
		{ -1,	10,	 0.949308762098587e4	},
		{ -1,	14,	-0.141043719679409e7	},
		{ -1,	16,	 0.849166230819026e7	},
		{ 0,	0,	 0.861095729446704	},
		{ 0,	2,	 0.323346442811720	},
		{ 1,	1,	 0.873281936020439	},
		{ 3,	1,	-0.436653048526683	},
		{ 5,	1,	 0.286596714529479	},
		{ 6,	1,	-0.131778331276228	},
		{ 8,	1,	 0.676682064330275e-2	},
	};

	int i, n;
	double ts;
	double ans = 0.0;

	if (h < eta_b3ab(p)) {
		i = 0;
		n = SIZE1;
		ts = 760.0;
		p = p * 0.01 + 0.24;
		h = h / 2300.0 - 0.615;
	} else {
		i = SIZE1;
		n = SIZE2;
		ts = 860.0;
		p = p * 0.01 + 0.298;
		h = h / 2800.0 - 0.72;
	}

	for (; i < n; ++i) {
		ans += coef[i].n * POWINT(p, coef[i].I) * POWINT(h, coef[i].J);
	}
	return ans * ts;
}

static double v_r3_ph(double p, double h)
{
	enum {
		SIZE1 = 32,
		SIZE2 = SIZE1 + 30,
       	};
	const if97_coef coef[SIZE2] = {
		/* region 3a */
		{ -12,	 6,	 0.529944062966028e-2	},
		{ -12,	 8,	-0.170099690234461	},
		{ -12,	12,	 0.111323814312927e2	},
		{ -12,	18,	-0.217898123145125e4	},
		{ -10,	 4,	-0.506061827980875e-3	},
		{ -10,	 7,	 0.556495239685324	},
		{ -10,	10,	-0.943672726094016e1	},
		{ -8,	 5,	-0.297856807561527	},
		{ -8,	12,	 0.939353943717186e2	},
		{ -6,	 3,	 0.192944939465981e-1	},
		{ -6,	 4,	 0.421740664704763	},
		{ -6,	22,	-0.368914126282330e7	},
		{ -4,	 2,	-0.737566847600639e-2	},
		{ -4,	 3,	-0.354753242424366	},
		{ -3,	 7,	-0.199768169338727e1	},
		{ -2,	 3,	 0.115456297059049e1	},
		{ -2,	16,	 0.568366875815960e4	},
		{ -1,	 0,	 0.808169540124668e-2	},
		{ -1,	 1,	 0.172416341519307	},
		{ -1,	 2,	 0.104270175292927e1	},
		{ -1,	 3,	-0.297691372792847	},
		{  0,	 0,	 0.560394465163593	},
		{  0,	 1,	 0.275234661176914	},
		{  1,	 0,	-0.148347894866012	},
		{  1,	 1,	-0.651142513478515e-1	},
		{  1,	 2,	-0.292468715386302e1	},
		{  2,	 0,	 0.664876096952665e-1	},
		{  2,	 2,	 0.352335014263844e1	},
		{  3,	 0,	-0.146340792313332e-1	},
		{  4,	 2,	-0.224503486668184e1	},
		{  5,	 2,	 0.110533464706142e1	},
		{  8,	 2,     -0.408757344495612e-1	},
		/* region 3b */
		{ -12,	 0,	-0.225196934336318e-8	},
		{ -12,	 1,	 0.140674363313486e-7	},
		{ -8,	 0,	 0.233784085280560e-5	},
		{ -8,	 1,	-0.331833715229001e-4	},
		{ -8,	 3,	 0.107956778514318e-2	},
		{ -8,	 6,	-0.271382067378863	},
		{ -8,	 7,	 0.107202262490333e1	},
		{ -8,	 8,	-0.853821329075382	},
		{ -6,	 0,	-0.215214194340526e-4	},
		{ -6,	 1,	 0.769656088222730e-3	},
		{ -6,	 2,	-0.431136580433864e-2	},
		{ -6,	 5,	 0.453342167309331	},
		{ -6,	 6,	-0.507749535873652	},
		{ -6,	10,	-0.100475154528389e3	},
		{ -4,	 3,	-0.219201924648793	},
		{ -4,	 6,	-0.321087965668917e1	},
		{ -4,	10,	 0.607567815637771e3	},
		{ -3,	 0,	 0.557686450685932e-3	},
		{ -3,	 2,	 0.187499040029550	},
		{ -2,	 1,	 0.905368030448107e-2	},
		{ -2,	 2,	 0.285417173048685	},
		{ -1,	 0,	 0.329924030996098e-1	},
		{ -1,	 1,	 0.239897419685483	},
		{ -1,	 4,	 0.482754995951394e1	},
		{ -1,	 5,	-0.118035753702231e2	},
		{  0,	 0,	 0.169490044091791	},
		{  1,	 0,	-0.179967222507787e-1	},
		{  1,	 1,	 0.371810116332674e-1	},
		{  2,	 2,	-0.536288335065096e-1	},
		{  2,	 6,	 0.160697101092520e1	},
	};

	int i, n;
	double vs;
	double ans = 0.0;

	if (h < eta_b3ab(p)) {
		i = 0;
		n = SIZE1;
		vs = 2.8e-3;
		p = p * 0.01 + 0.128;
		h = h / 2100.0 - 0.727;
	} else {
		i = SIZE1;
		n = SIZE2;
		vs = 8.8e-3;
		p = p * 0.01 + 0.0661;
		h = h / 2800.0 - 0.72;
	}

	for (; i < n; ++i) {
		ans += coef[i].n * POWINT(p, coef[i].I) * POWINT(h, coef[i].J);
	}
	return ans * vs;
}

static void get_phi_r3_ph(double *rhot, void *xphi, double *ph, double *dph)
{
	double rho = rhot[0];
	double t = rhot[1];

	iapws_phi *phi = (iapws_phi *)(xphi);
	phi_r3(rho, t, phi);

	ph[0] = phi->d10 * rho * phi->R * t * 1e-3 - phi->p;
	ph[1] = (phi->d10 + phi->d01) * phi->R * t - phi->h;

	dph[0] = (phi->d10 * 2.0 + phi->d20) * phi->R * t * 1e-3;
	dph[1] = (phi->d10 + phi->d20 + phi->d11) / rho * phi->R * t;
	dph[2] = (phi->d10 - phi->d11) * rho * phi->R * 1e-3;
	dph[3] = (phi->d10 - phi->d11 - phi->d02) * phi->R;
}

static int phi_r3_ph(double p, double h, iapws_phi *phi)
{
	int maxiter = nroot_maxiter;
	double tolf = nroot_tolf;
	double tolx = nroot_tolx;
	double rhot[2] = { 1.0 / v_r3_ph(p, h), t_r3_ph(p, h) };
	phi->p = p;
	phi->h = h;
	return nroot2(get_phi_r3_ph, rhot, phi,
			&tolf, &tolx, &maxiter, nroot_verbose);
}


/* Region 4 */

static const double coef_r4[10] = {
	0.11670521452767e4, -0.72421316703206e6,
	-0.17073846940092e2, 0.12020824702470e5,
	-0.32325550322333e7, 0.14915108613530e2,
	-0.48232657361591e4, 0.40511340542057e6,
	-0.23855557567849, 0.65017534844798e3,
};

static double pi_r4(double theta)
{
	theta += coef_r4[8] / (theta - coef_r4[9]);
	double theta2 = theta * theta;
	double a = theta2 + coef_r4[0] * theta + coef_r4[1];
	double b = coef_r4[2] * theta2 + coef_r4[3] * theta + coef_r4[4];
	double c = coef_r4[5] * theta2 + coef_r4[6] * theta + coef_r4[7];
	double ans = 2.0 * c / (-b + sqrt(b * b - 4.0 * a * c));
	ans *= ans;
	return ans * ans;
}

static double theta_r4(double pi)
{
	double beta = POW(pi, .25);
	double beta2 = beta * beta;
	double e = beta2 + coef_r4[2] * beta + coef_r4[5];
	double f = coef_r4[0] * beta2 + coef_r4[3] * beta + coef_r4[6];
	double g = coef_r4[1] * beta2 + coef_r4[4] * beta + coef_r4[7];
	double d = -2.0 * g / (f + sqrt(f * f - 4.0 * e * g));
	return (coef_r4[9] + d - sqrt((coef_r4[9] + d) * (coef_r4[9] + d) -
				4.0 * (coef_r4[8] + coef_r4[9] * d))) * .5;
}

/* Region 5 */

static void gamma_r5(double p, double t, iapws_phi *gamma)
{
	enum {
		SIZE1 = 6,
		SIZE2 = SIZE1 + 6,
	};
	const if97_coef coef[SIZE2] = {
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

	const double pi = p;
	const double tau = 1000.0 / t;

	int i;
	double xn;

	/* gamma_0 */
	gamma->type = IAPWS_GAMMA;
	gamma->d00 = log(pi);
	gamma->d10 = 1.0;
	gamma->d01 = 0.0;
	gamma->d11 = 0.0;
	gamma->d20 = -1.0;
	gamma->d02 = 0.0;
	gamma->p = p;
	gamma->t = t;
	for (i = 0; i < SIZE1; ++i) {
		xn = coef[i].n * POWINT(tau, coef[i].J);
		gamma->d00 += xn;
		gamma->d01 += xn * coef[i].J;
		gamma->d02 += xn * coef[i].J * (coef[i].J - 1);
	}

	/* gamma_r */
	for (i = SIZE1; i < SIZE2; ++i) {
		xn = coef[i].n * POWINT(pi, coef[i].I) * POWINT(tau, coef[i].J);
		gamma->d00 += xn;
		gamma->d10 += xn * coef[i].I;
		gamma->d01 += xn * coef[i].J;
		gamma->d11 += xn * coef[i].I * coef[i].J;
		gamma->d20 += xn * coef[i].I * (coef[i].I - 1);
		gamma->d02 += xn * coef[i].J * (coef[i].J - 1);
	}
}

/* Boundary between regions 2 and 3 */

static const double coef_b23[5] = {
	0.34805185628969e+3, -0.11671859879975e+1, 0.10192970039326e-2,
	0.57254459862746e+3,  0.13918839778870e+2,
};

static double pi_b23(double theta)
{
	return coef_b23[0] + coef_b23[1] * theta + coef_b23[2] * POW2(theta);
}

static double theta_b23(double pi)
{
	return coef_b23[3] + sqrt((pi - coef_b23[4]) / coef_b23[2]);
}

/* Boundary between regions 2b and 2c */

static const double coef_b2bc[5] = {
	0.90584278514723e3, -0.67955786399241e0, 0.12809002730136e-3,
	0.26526571908428e4,  0.45257578905948e1
};

static double pi_b2bc(double eta)
{
	return coef_b2bc[0] + coef_b2bc[1] * eta + coef_b2bc[2] * POW2(eta);
}

//static double eta_b2bc(double pi)
//{
//	return coef_b2bc[3] + sqrt((pi - coef_b2bc[4]) / coef_b2bc[2]);
//}

/* Boundary between regions 3a and 3b */

static double eta_b3ab(double pi)
{
	static double n[4] = {
		 0.201464004206875e+4,
		 0.374696550136983e+1,
		-0.219921901054187e-1,
		 0.875131686009950e-4,
	};
	double pi2 = pi * pi;
	return n[0] + n[1] * pi + n[2] * pi2 + n[3] * pi * pi2;
}

#if 0
#include <stdio.h>
int main()
{
	int err, reg;
	iapws_phi gamma;
	gamma.R = IF97_R;

	phi_r3(IAPWS_RHOC, IAPWS_TC, &gamma);
	printf("%.8e %.8e %.8e %.8e %.8e\n",
			iapws_p(&gamma),
			iapws_t(&gamma),
			iapws_rho(&gamma),
			iapws_s(&gamma),
			iapws_h(&gamma));
}
#endif
