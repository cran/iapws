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
#include "melt08.h"
#include "nroot.h"
#include "pow.h"

/* Static function forward declarations */
static void gamma_r1(double p, double t, iapws_phi *gamma);
static int gamma_r1_ph(double p, double h, iapws_phi *gamma);
static double t_r1_ph(double p, double h);

static void gamma_r2(double p, double t, int meta, iapws_phi *gamma);
static int gamma_r2_ph(double p, double h, iapws_phi *gamma);
static double t_r2_ph(double p, double h);

static void phi_r3(double rho, double t, iapws_phi *phi);
static int phi_r3_pt(double p, double t, iapws_phi *phi);
static int phi_r3_ph(double p, double h, iapws_phi *phi);
static double v_r3_pt(double p, double t);
static double t_r3_ph(double p, double h);
static double v_r3_ph(double p, double h);

static double pi_r4(double theta);
static double theta_r4(double pi);

static void gamma_r5(double p, double t, iapws_phi *gamma);

static double pi_b23(double theta);
static double theta_b23(double pi);
static double p_b34_h(double h);

static double pi_b2bc(double eta);
static double eta_b3ab(double pi);

typedef struct {
	int I;
	int J;
	double n;
} if97_coef;

/* Exported functions */

iapws_state_id if97_state_pt(double p, double t)
{
	double ps;
	if (t >= IAPWS_TT && t < IAPWS_TC && p < 620.0) {
		ps = pi_r4(t);
		if (p > ps) return IAPWS_LIQUID;
		if (p < ps) return IAPWS_GAS;
		return IAPWS_SAT;
	} else if (t >= IAPWS_TC) {
		if (p < IAPWS_PC) return IAPWS_GAS;
		else return IAPWS_CRIT;
	}
	return melt_sub_state(p, t);
}

if97_region_id if97_region_pt(double p, double t)
{
	double ps;
	if (t >= 273.15 && t <= 623.15) {
		ps = pi_r4(t);
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
	double ts;
	iapws_phi gamma;
	gamma.R = IF97_R;

	if (p <= pi_r4(273.15)) {
		gamma_r2(p, 273.15, 0, &gamma);
		if (h >= iapws_h(&gamma)) return IF97_REGION_2;
		return IF97_REGION_UNDEF;
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
		if (h <= iapws_h(&gamma)) return IF97_REGION_1;
		gamma_r2(p, theta_b23(p), 0, &gamma);
		if (h >= iapws_h(&gamma)) return IF97_REGION_2;
		if (p >= IAPWS_PC || p >= p_b34_h(h)) return IF97_REGION_3;
		return IF97_REGION_4;
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

//#define PEPS	1.001
#define PEPS	1.0001

int if97_gamma_pt(double p, double t, iapws_state_id state, iapws_phi *gamma)
{
	if97_region_id reg = if97_region_pt(p, t);
	int meta = 0;
	int err = 0;
	//double rho0 = 650.0;
	double p0 = p;

	if (state == IAPWS_LIQUID) {
		if (reg == IF97_REGION_1) {
		} else if (reg == IF97_REGION_2) {
			reg = IF97_REGION_1;
		} else if (reg == IF97_REGION_3) {
			if (t <= IAPWS_TC && (p0 = pi_r4(t) * PEPS) < p) {
				p0 = p;
			}
		} else {
			reg = IF97_REGION_UNDEF;
		}
	} else if (state == IAPWS_GAS) {
		if (reg == IF97_REGION_2) {
		} else if (reg == IF97_REGION_1) {
			meta = (p < 10.0 ? 1 : 0);
			reg = IF97_REGION_2;
		} else if (reg == IF97_REGION_3) {
			//rho0 = 150.0;
			if (t <= IAPWS_TC && (p0 = pi_r4(t) / PEPS) > p) {
				p0 = p;
			}
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
			//gamma->rho = rho0;
			gamma->rho = 1.0 / v_r3_pt(p0, t);
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
	if97_region_id reg = if97_region_ph(p, h);
	gamma->R = IF97_R;
	switch (reg) {
		case IF97_REGION_1:
			return gamma_r1_ph(p, h, gamma);
		case IF97_REGION_2:
			return gamma_r2_ph(p, h, gamma);
		case IF97_REGION_3:
			return phi_r3_ph(p, h, gamma);
		default:
			return -1;
	}
	assert(0);
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
	double xi, xj, xn;

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
		if (i != 0) {
			xi *= powint(xp, coef[i].I - coef[i - 1].I);
			xj *= powint(xt, coef[i].J - coef[i - 1].J);
		} else {
			xi = powint(xp, coef[i].I);
			xj = powint(xt, coef[i].J);
		}
		xn = coef[i].n * xi * xj;
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
	double xi, xj;
	double ans = 0.0;

	h = h * 4.0e-4 + 1.0;
	for (i = 0; i < SIZE; ++i) {
		if (i != 0) {
			xi *= powint(p, coef[i].I - coef[i - 1].I);
			xj *= powint(h, coef[i].J - coef[i - 1].J);
		} else {
			xi = powint(p, coef[i].I);
			xj = powint(h, coef[i].J);
		}
		ans += coef[i].n * xi * xj;
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
	nroot_control ctrl = nroot_default;
	gamma->p = p;
	gamma->h = h;
	gamma->t = t_r1_ph(p, h);
	return nroot1(get_gamma_r1_ph, &gamma->t, gamma, &ctrl);
}

/* Region 2 */

static void gamma_r2(double p, double t, int meta, iapws_phi *gamma)
{
	assert(meta == 0 || meta == 1);

	enum {
		SIZE0 = 9,
		SIZE1 = SIZE0 + 43,
		SIZE2 = SIZE0 + 13,
	};
	const if97_coef coef1[SIZE1] = {
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
	};
	const if97_coef coef2[SIZE2] = {
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

	int i, n;
	const if97_coef *coef;
	double xi, xj, xn;

	if (meta == 0) {
		n = SIZE1;
		coef = coef1;
	} else {
		n = SIZE2;
		coef = coef2;
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
	for (i = 0; i < SIZE0; ++i) {
		xn = coef[i].n * powint(tau, coef[i].J);
		gamma->d00 += xn;
		gamma->d01 += xn * coef[i].J;
		gamma->d02 += xn * coef[i].J * (coef[i].J - 1);
	}

	/* gamma_r */
	iapws_phi gamma_r = { IAPWS_GAMMA };
	for (i = SIZE0; i < n; ++i) {
		if (i != SIZE0) {
			xi *= powint(pi, coef[i].I - coef[i - 1].I);
			xj *= powint(xt, coef[i].J - coef[i - 1].J);
		} else {
			xi = powint(pi, coef[i].I);
			xj = powint(xt, coef[i].J);
		}
		xn = coef[i].n * xi * xj;
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
		SIZEA = 34,
		SIZEB = 38,
		SIZEC = 23,
	};
	const if97_coef coefa[SIZEA] = {
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
	};
	const if97_coef coefb[SIZEB] = {
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
	};
	const if97_coef coefc[SIZEC] = {
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
	const if97_coef *coef;
	double xi, xj;
	double ans = 0.0;

	if (p <= 4.0) {
		n = SIZEA;
		coef = coefa;
		h = h * 5.0e-4 - 2.1;
	} else if (p < pi_b2bc(h)) {
		n = SIZEB;
		coef = coefb;
		p -= 2.0;
		h = h * 5.0e-4 - 2.6;
	} else {
		n = SIZEC;
		coef = coefc;
		p += 25.0;
		h = h * 5.0e-4 - 1.8;
	}

	for (i = 0; i < n; ++i) {
		if (i != 0) {
			xi *= powint(p, coef[i].I - coef[i - 1].I);
			//xj *= powint(h, coef[i].J - coef[i - 1].J);
			xj = powint(h, coef[i].J);
		} else {
			xi = powint(p, coef[i].I);
			xj = powint(h, coef[i].J);
		}
		ans += coef[i].n * xi * xj;
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
	nroot_control ctrl = nroot_default;
	gamma->p = p;
	gamma->h = h;
	gamma->t = t_r2_ph(p, h);
	return nroot1(get_gamma_r2_ph, &gamma->t, gamma, &ctrl);
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
	double xi, xj, xn;

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
		if (i != 1) {
			xi *= powint(delta, coef[i].I - coef[i - 1].I);
			//xj *= powint(tau, coef[i].J - coef[i - 1].J);
			xj = powint(tau, coef[i].J);
		} else {
			xi = powint(delta, coef[i].I);
			xj = powint(tau, coef[i].J);
		}
		xn = coef[i].n * xi * xj;
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
	nroot_control ctrl = nroot_default;
	phi->p = p;
	phi->t = t;
	return nroot1(get_phi_r3_pt, &phi->rho, phi, &ctrl);
}

static double t_r3_ph(double p, double h)
{
	enum {
		SIZEA = 31,
		SIZEB = 33,
	};
	const if97_coef coefa[SIZEA] = {
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
	};
	const if97_coef coefb[SIZEB] = {
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
	const if97_coef *coef;
	double xi, xj, ts;
	double ans = 0.0;

	if (h < eta_b3ab(p)) {
		n = SIZEA;
		coef = coefa;
		ts = 760.0;
		p = p * 0.01 + 0.24;
		h = h / 2300.0 - 0.615;
	} else {
		n = SIZEB;
		coef = coefb;
		ts = 860.0;
		p = p * 0.01 + 0.298;
		h = h / 2800.0 - 0.72;
	}

	for (i = 0; i < n; ++i) {
		if (i != 0) {
			xi *= powint(p, coef[i].I - coef[i - 1].I);
			//xj *= powint(h, coef[i].J - coef[i - 1].J);
			xj = powint(h, coef[i].J);
		} else {
			xi = powint(p, coef[i].I);
			xj = powint(h, coef[i].J);
		}
		ans += coef[i].n * xi * xj;
	}
	return ans * ts;
}

static double v_r3_ph(double p, double h)
{
	enum {
		SIZEA = 32,
		SIZEB = 30,
	};
	const if97_coef coefa[SIZEA] = {
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
	};
	const if97_coef coefb[SIZEB] = {
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
	const if97_coef *coef;
	double xi, xj, vs;
	double ans = 0.0;

	if (h < eta_b3ab(p)) {
		n = SIZEA;
		coef = coefa;
		vs = 2.8e-3;
		p = p * 0.01 + 0.128;
		h = h / 2100.0 - 0.727;
	} else {
		n = SIZEB;
		coef = coefb;
		vs = 8.8e-3;
		p = p * 0.01 + 0.0661;
		h = h / 2800.0 - 0.72;
	}

	for (i = 0; i < n; ++i) {
		if (i != 0) {
			xi *= powint(p, coef[i].I - coef[i - 1].I);
			//xj *= powint(h, coef[i].J - coef[i - 1].J);
			xj = powint(h, coef[i].J);
		} else {
			xi = powint(p, coef[i].I);
			xj = powint(h, coef[i].J);
		}
		ans += coef[i].n * xi * xj;
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
	nroot_control ctrl = nroot_default;
	double rhot[2] = { 1.0 / v_r3_ph(p, h), t_r3_ph(p, h) };
	phi->p = p;
	phi->h = h;
	return nroot2(get_phi_r3_ph, rhot, phi, &ctrl);
}

/* Region 4 */

static const double coef_r4[10] = {
	 0.11670521452767e4, -0.72421316703206e6,
	-0.17073846940092e2,  0.12020824702470e5,
	-0.32325550322333e7,  0.14915108613530e2,
	-0.48232657361591e4,  0.40511340542057e6,
	-0.23855557567849,    0.65017534844798e3,
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
		SIZE0 = 6,
		SIZER = 6,
	};
	const if97_coef coef0[SIZE0] = {
		{ 0,	 0,	-0.13179983674201e2  },
		{ 0,	 1,	 0.68540841634434e1  },
		{ 0,	-3,	-0.24805148933466e-1 },
		{ 0,	-2,	 0.36901534980333    },
		{ 0,	-1,	-0.31161318213925e1  },
		{ 0,	 2,	-0.32961626538917    },
	};
	const if97_coef coefr[SIZER] = {
		{ 1,	1,	 0.15736404855259e-2 },
		{ 1,	2,	 0.90153761673944e-3 },
		{ 1,	3,	-0.50270077677648e-2 },
		{ 2,	3,	 0.22440037409485e-5 },
		{ 2,	9,	-0.41163275453471e-5 },
		{ 3,	7,	 0.37919454822955e-7 },
	};

	const double pi = p;
	const double tau = 1000.0 / t;

	int i;
	double xi, xj, xn;

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
	for (i = 0; i < SIZE0; ++i) {
		xn = coef0[i].n * powint(tau, coef0[i].J);
		gamma->d00 += xn;
		gamma->d01 += xn * coef0[i].J;
		gamma->d02 += xn * coef0[i].J * (coef0[i].J - 1);
	}

	/* gamma_r */
	for (i = 0; i < SIZER; ++i) {
		if (i != 0) {
			xi *= powint(pi, coefr[i].I - coefr[i - 1].I);
			xj *= powint(tau, coefr[i].J - coefr[i - 1].J);
		} else {
			xi = powint(pi, coefr[i].I);
			xj = powint(tau, coefr[i].J);
		}
		xn = coefr[i].n * xi * xj;
		gamma->d00 += xn;
		gamma->d10 += xn * coefr[i].I;
		gamma->d01 += xn * coefr[i].J;
		gamma->d11 += xn * coefr[i].I * coefr[i].J;
		gamma->d20 += xn * coefr[i].I * (coefr[i].I - 1);
		gamma->d02 += xn * coefr[i].J * (coefr[i].J - 1);
	}
}

/* Boundary between regions 2 and 3 */

static const double coef_b23[5] = {
	 0.34805185628969e+3,
	-0.11671859879975e+1,
	 0.10192970039326e-2,
	 0.57254459862746e+3,
	 0.13918839778870e+2,
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
	 0.90584278514723e3,
	-0.67955786399241e0,
	 0.12809002730136e-3,
	 0.26526571908428e4,
	 0.45257578905948e1,
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

static double p_b34_h(double h)
{
	enum { SIZE = 14 };
	const if97_coef coef[SIZE] = {
		{  0,	 0,	 0.600073641753024    },
		{  1,	 1,	-0.936203654849857e1  },
		{  1,	 3,	 0.246590798594147e2  },
		{  1,	 4,	-0.107014222858224e3  },
		{  1,	36,	-0.915821315805768e14 },
		{  5,	 3,	-0.862332011700662e4  },
		{  7,	 0,	-0.235837344740032e2  },
		{  8,	24,	 0.252304969384128e18 },
		{ 14,	16,	-0.389718771997719e19 },
		{ 20,	16,	-0.333775713645296e23 },
		{ 22,	 3,	 0.356499469636328e11 },
		{ 24,	18,	-0.148547544720641e27 },
		{ 28,	 8,	 0.330611514838798e19 },
		{ 36,	24,	 0.813641294467829e38 },
	};

	int i;
	double xi, xj;
	double ans = 0.0;

	const double etai = h / 2600.0 - 1.02;
	const double etaj = h / 2600.0 - 0.608;

	for (i = 0; i < SIZE; ++i) {
		if (i != 0) {
			xi *= powint(etai, coef[i].I - coef[i - 1].I);
			xj *= powint(etaj, coef[i].J - coef[i - 1].J);
		} else {
			xi = powint(etai, coef[i].I);
			xj = powint(etaj, coef[i].J);
		}
		ans += coef[i].n * xi * xj;
	}
	return ans * 22.0;
}

enum if97_b3 {
	B3AB,
	B3CD,
	B3EF,
	B3GH,
	B3IJ,
	B3JK,
	B3MN,
	B3OP,
	B3QU,
	B3RX,
	B3UV,
	B3WX,
};
static double t_b3(double p, enum if97_b3 b)
{
	enum {
		SIZEAB = 5,
		SIZECD = 4,
		SIZEGH = 5,
		SIZEIJ = 5,
		SIZEJK = 5,
		SIZEMN = 4,
		SIZEOP = 5,
		SIZEQU = 4,
		SIZERX = 4,
		SIZEUV = 4,
		SIZEWX = 5,
	};
	const double coefab[SIZEAB] = {
		 0.918419702359447e3,
		-0.191887498864292e4,
		 0.154793642129415e4,
		-0.187661219490113e3,
		 0.213144632222113e2,
	};
	const double coefcd[SIZECD] = {
		 0.585276966696349e3,
		 0.278233532206915e1,
		-0.127283549295878e-1,
		 0.159090746562729e-3,
	};
	const double coefgh[SIZEGH] = {
		-0.249284240900418e5,
		 0.428143584791546e4,
		-0.269029173140130e3,
		 0.751608051114157e1,
		-0.787105249910383e-1,
	};
	const double coefij[SIZEIJ] = {
		 0.584814781649163e3,
		-0.616179320924617,
		 0.260763050899562,
		-0.587071076864459e-2,
		 0.515308185433082e-4,
	};
	const double coefjk[SIZEJK] = {
		 0.617229772068439e3,
		-0.770600270141675e1,
		 0.697072596851896,
		-0.157391839848015e-1,
		 0.137897492684194e-3,
	};
	const double coefmn[SIZEMN] = {
		 0.535339483742384e3,
		 0.761978122720128e1,
		-0.158365725441648,
		 0.192871054508108e-2,
	};
	const double coefop[SIZEOP] = {
		-0.152313732937084e4,
		 0.773845935768222e3,
		 0.969461372400213e3,
		-0.332500170441278e3,
		 0.642859598466067e2,
	};
	const double coefqu[SIZEQU] = {
		 0.565603648239126e3,
		 0.529062258221222e1,
		-0.102020639611016,
		 0.122240301070145e-2,
	};
	const double coefrx[SIZERX] = {
		 0.584561202520006e3,
		-0.102961025163669e1,
		 0.243293362700452,
		-0.294905044740799e-2,
	};
	const double coefuv[SIZEUV] = {
		 0.528199646263062e3,
		 0.890579602135307e1,
		-0.222814134903755,
		 0.286791682263697e-2,
	};
	const double coefwx[SIZEWX] = {
		0.873371668682417e3,
		0.329196213998375e3,
		0.728052609145380e1,
		0.973505869861952e2,
		0.147370491183191e2,
	};

	int i, n;
	const double *coef;
	double pi = 1.0;
	double ans = 0.0;

	switch (b) {
		case B3AB:
			n = SIZEAB;
			coef = coefab;
			p = log(p);
			pi = 1.0 / (p * p);
			break;
		case B3CD:
			n = SIZECD;
			coef = coefcd;
			break;
		case B3EF:
			return 3.727888004 * (p - IAPWS_PC) + IAPWS_TC;
		case B3GH:
			n = SIZEGH;
			coef = coefgh;
			break;
		case B3IJ:
			n = SIZEIJ;
			coef = coefij;
			break;
		case B3JK:
			n = SIZEJK;
			coef = coefjk;
			break;
		case B3MN:
			n = SIZEMN;
			coef = coefmn;
			break;
		case B3OP:
			n = SIZEOP;
			coef = coefop;
			p = log(p);
			pi = 1.0 / (p * p);
			break;
		case B3QU:
			n = SIZEQU;
			coef = coefqu;
			break;
		case B3RX:
			n = SIZERX;
			coef = coefrx;
			break;
		case B3UV:
			n = SIZEUV;
			coef = coefuv;
			break;
		case B3WX:
			n = SIZEWX;
			coef = coefwx;
			p = log(p);
			pi = 1.0 / (p * p);
			break;
		default:
			assert(0);
	}
	for (i = 0; i < n; ++i, pi *= p) ans += coef[i] * pi;
	return ans;
}

static double v_r3_pt(double p, double t)
{
	enum {
		SIZEA = 30,
		SIZEB = 32,
		SIZEC = 35,
		SIZED = 38,
		SIZEE = 29,
		SIZEF = 42,
		SIZEG = 38,
		SIZEH = 29,
		SIZEI = 42,
		SIZEJ = 29,
		SIZEK = 34,
		SIZEL = 43,
		SIZEM = 40,
		SIZEN = 39,
		SIZEO = 24,
		SIZEP = 27,
		SIZEQ = 24,
		SIZER = 27,
		SIZES = 29,
		SIZET = 33,
		SIZEU = 38,
		SIZEV = 39,
		SIZEW = 35,
		SIZEX = 36,
		SIZEY = 20,
		SIZEZ = 23,
	};
	const if97_coef coefa[SIZEA] = {
		{ -12,	  5,	 0.110879558823853e-02	},
		{ -12,	 10,	 0.572616740810616e+03	},
		{ -12,	 12,	-0.767051948380852e+05	},
		{ -10,	  5,	-0.253321069529674e-01	},
		{ -10,	 10,	 0.628008049345689e+04	},
		{ -10,	 12,	 0.234105654131876e+06	},
		{  -8,	  5,	 0.216867826045856e+00	},
		{  -8,	  8,	-0.156237904341963e+03	},
		{  -8,	 10,	-0.269893956176613e+05	},
		{  -6,	  1,	-0.180407100085505e-03	},
		{  -5,	  1,	 0.116732227668261e-02	},
		{  -5,	  5,	 0.266987040856040e+02	},
		{  -5,	 10,	 0.282776617243286e+05	},
		{  -4,	  8,	-0.242431520029523e+04	},
		{  -3,	  0,	 0.435217323022733e-03	},
		{  -3,	  1,	-0.122494831387441e-01	},
		{  -3,	  3,	 0.179357604019989e+01	},
		{  -3,	  6,	 0.442729521058314e+02	},
		{  -2,	  0,	-0.593223489018342e-02	},
		{  -2,	  2,	 0.453186261685774e+00	},
		{  -2,	  3,	 0.135825703129140e+01	},
		{  -1,	  0,	 0.408748415856745e-01	},
		{  -1,	  1,	 0.474686397863312e+00	},
		{  -1,	  2,	 0.118646814997915e+01	},
		{   0,	  0,	 0.546987265727549e+00	},
		{   0,	  1,	 0.195266770452643e+00	},
		{   1,	  0,	-0.502268790869663e-01	},
		{   1,	  2,	-0.369645308193377e+00	},
		{   2,	  0,	 0.633828037528420e-02	},
		{   2,	  2,	 0.797441793901017e-01	},
	};
	const if97_coef coefb[SIZEB] = {
		{ -12,	 10,	-0.827670470003621e-01	},
		{ -12,	 12,	 0.416887126010565e+02	},
		{ -10,	  8,	 0.483651982197059e-01	},
		{ -10,	 14,	-0.291032084950276e+05	},
		{  -8,	  8,	-0.111422582236948e+03	},
		{  -6,	  5,	-0.202300083904014e-01	},
		{  -6,	  6,	 0.294002509338515e+03	},
		{  -6,	  8,	 0.140244997609658e+03	},
		{  -5,	  5,	-0.344384158811459e+03	},
		{  -5,	  8,	 0.361182452612149e+03	},
		{  -5,	 10,	-0.140699677420738e+04	},
		{  -4,	  2,	-0.202023902676481e-02	},
		{  -4,	  4,	 0.171346792457471e+03	},
		{  -4,	  5,	-0.425597804058632e+01	},
		{  -3,	  0,	 0.691346085000334e-05	},
		{  -3,	  1,	 0.151140509678925e-02	},
		{  -3,	  2,	-0.416375290166236e-01	},
		{  -3,	  3,	-0.413754957011042e+02	},
		{  -3,	  5,	-0.506673295721637e+02	},
		{  -2,	  0,	-0.572212965569023e-03	},
		{  -2,	  2,	 0.608817368401785e+01	},
		{  -2,	  5,	 0.239600660256161e+02	},
		{  -1,	  0,	 0.122261479925384e-01	},
		{  -1,	  2,	 0.216356057692938e+01	},
		{   0,	  0,	 0.398198903368642e+00	},
		{   0,	  1,	-0.116892827834085e+00	},
		{   1,	  0,	-0.102845919373532e+00	},
		{   1,	  2,	-0.492676637589284e+00	},
		{   2,	  0,	 0.655540456406790e-01	},
		{   3,	  2,	-0.240462535078530e+00	},
		{   4,	  0,	-0.269798180310075e-01	},
		{   4,	  1,	 0.128369435967012e+00	},
	};
	const if97_coef coefc[SIZEC] = {
		{ -12,	  6,	 0.311967788763030e+01	},
		{ -12,	  8,	 0.276713458847564e+05	},
		{ -12,	 10,	 0.322583103403269e+08	},
		{ -10,	  6,	-0.342416065095363e+03	},
		{ -10,	  8,	-0.899732529907377e+06	},
		{ -10,	 10,	-0.793892049821251e+08	},
		{  -8,	  5,	 0.953193003217388e+02	},
		{  -8,	  6,	 0.229784742345072e+04	},
		{  -8,	  7,	 0.175336675322499e+06	},
		{  -6,	  8,	 0.791214365222792e+07	},
		{  -5,	  1,	 0.319933345844209e-04	},
		{  -5,	  4,	-0.659508863555767e+02	},
		{  -5,	  7,	-0.833426563212851e+06	},
		{  -4,	  2,	 0.645734680583292e-01	},
		{  -4,	  8,	-0.382031020570813e+07	},
		{  -3,	  0,	 0.406398848470079e-04	},
		{  -3,	  3,	 0.310327498492008e+02	},
		{  -2,	  0,	-0.892996718483724e-03	},
		{  -2,	  4,	 0.234604891591616e+03	},
		{  -2,	  5,	 0.377515668966951e+04	},
		{  -1,	  0,	 0.158646812591361e-01	},
		{  -1,	  1,	 0.707906336241843e+00	},
		{  -1,	  2,	 0.126016225146570e+02	},
		{   0,	  0,	 0.736143655772152e+00	},
		{   0,	  1,	 0.676544268999101e+00	},
		{   0,	  2,	-0.178100588189137e+02	},
		{   1,	  0,	-0.156531975531713e+00	},
		{   1,	  2,	 0.117707430048158e+02	},
		{   2,	  0,	 0.840143653860447e-01	},
		{   2,	  1,	-0.186442467471949e+00	},
		{   2,	  3,	-0.440170203949645e+02	},
		{   2,	  7,	 0.123290423502494e+07	},
		{   3,	  0,	-0.240650039730845e-01	},
		{   3,	  7,	-0.107077716660869e+07	},
		{   8,	  1,	 0.438319858566475e-01	},
	};
	const if97_coef coefd[SIZED] = {
		{ -12,	  4,	-0.452484847171645e-09	},
		{ -12,	  6,	 0.315210389538801e-04	},
		{ -12,	  7,	-0.214991352047545e-02	},
		{ -12,	 10,	 0.508058874808345e+03	},
		{ -12,	 12,	-0.127123036845932e+08	},
		{ -12,	 16,	 0.115371133120497e+13	},
		{ -10,	  0,	-0.197805728776273e-15	},
		{ -10,	  2,	 0.241554806033972e-10	},
		{ -10,	  4,	-0.156481703640525e-05	},
		{ -10,	  6,	 0.277211346836625e-02	},
		{ -10,	  8,	-0.203578994462286e+02	},
		{ -10,	 10,	 0.144369489909053e+07	},
		{ -10,	 14,	-0.411254217946539e+11	},
		{  -8,	  3,	 0.623449786243773e-05	},
		{  -8,	  7,	-0.221774281146038e+02	},
		{  -8,	  8,	-0.689315087933158e+05	},
		{  -8,	 10,	-0.195419525060713e+08	},
		{  -6,	  6,	 0.316373510564015e+04	},
		{  -6,	  8,	 0.224040754426988e+07	},
		{  -5,	  1,	-0.436701347922356e-05	},
		{  -5,	  2,	-0.404213852833996e-03	},
		{  -5,	  5,	-0.348153203414663e+03	},
		{  -5,	  7,	-0.385294213555289e+06	},
		{  -4,	  0,	 0.135203700099403e-06	},
		{  -4,	  1,	 0.134648383271089e-03	},
		{  -4,	  7,	 0.125031835351736e+06	},
		{  -3,	  2,	 0.968123678455841e-01	},
		{  -3,	  4,	 0.225660517512438e+03	},
		{  -2,	  0,	-0.190102435341872e-03	},
		{  -2,	  1,	-0.299628410819229e-01	},
		{  -1,	  0,	 0.500833915372121e-02	},
		{  -1,	  1,	 0.387842482998411e+00	},
		{  -1,	  5,	-0.138535367777182e+04	},
		{   0,	  0,	 0.870745245971773e+00	},
		{   0,	  2,	 0.171946252068742e+01	},
		{   1,	  0,	-0.326650121426383e-01	},
		{   1,	  6,	 0.498044171727877e+04	},
		{   3,	  0,	 0.551478022765087e-02	},
	};
	const if97_coef coefe[SIZEE] = {
		{ -12,	 14,	 0.715815808404721e+09	},
		{ -12,	 16,	-0.114328360753449e+12	},
		{ -10,	  3,	 0.376531002015720e-11	},
		{ -10,	  6,	-0.903983668691157e-04	},
		{ -10,	 10,	 0.665695908836252e+06	},
		{ -10,	 14,	 0.535364174960127e+10	},
		{ -10,	 16,	 0.794977402335603e+11	},
		{  -8,	  7,	 0.922230563421437e+02	},
		{  -8,	  8,	-0.142586073991215e+06	},
		{  -8,	 10,	-0.111796381424162e+07	},
		{  -6,	  6,	 0.896121629640760e+04	},
		{  -5,	  6,	-0.669989239070491e+04	},
		{  -4,	  2,	 0.451242538486834e-02	},
		{  -4,	  4,	-0.339731325977713e+02	},
		{  -3,	  2,	-0.120523111552278e+01	},
		{  -3,	  6,	 0.475992667717124e+05	},
		{  -3,	  7,	-0.266627750390341e+06	},
		{  -2,	  0,	-0.153314954386524e-03	},
		{  -2,	  1,	 0.305638404828265e+00	},
		{  -2,	  3,	 0.123654999499486e+03	},
		{  -2,	  4,	-0.104390794213011e+04	},
		{  -1,	  0,	-0.157496516174308e-01	},
		{   0,	  0,	 0.685331118940253e+00	},
		{   0,	  1,	 0.178373462873903e+01	},
		{   1,	  0,	-0.544674124878910e+00	},
		{   1,	  4,	 0.204529931318843e+04	},
		{   1,	  6,	-0.228342359328752e+05	},
		{   2,	  0,	 0.413197481515899e+00	},
		{   2,	  2,	-0.341931835910405e+02	},
	};
	const if97_coef coeff[SIZEF] = {
		{   0,	 -3,	-0.251756547792325e-07	},
		{   0,	 -2,	 0.601307193668763e-05	},
		{   0,	 -1,	-0.100615977450049e-02	},
		{   0,	  0,	 0.999969140252192e+00	},
		{   0,	  1,	 0.214107759236486e+01	},
		{   0,	  2,	-0.165175571959086e+02	},
		{   1,	 -1,	-0.141987303638727e-02	},
		{   1,	  1,	 0.269251915156554e+01	},
		{   1,	  2,	 0.349741815858722e+02	},
		{   1,	  3,	-0.300208695771783e+02	},
		{   2,	  0,	-0.131546288252539e+01	},
		{   2,	  1,	-0.839091277286169e+01	},
		{   3,	 -5,	 0.181545608337015e-09	},
		{   3,	 -2,	-0.591099206478909e-03	},
		{   3,	  0,	 0.152115067087106e+01	},
		{   4,	 -3,	 0.252956470663225e-04	},
		{   5,	 -8,	 0.100726265203786e-14	},
		{   5,	  1,	-0.149774533860650e+01	},
		{   6,	 -6,	-0.793940970562969e-09	},
		{   7,	 -4,	-0.150290891264717e-03	},
		{   7,	  1,	 0.151205531275133e+01	},
		{  10,	 -6,	 0.470942606221652e-05	},
		{  12,	-10,	 0.195049710391712e-12	},
		{  12,	 -8,	-0.911627886266077e-08	},
		{  12,	 -4,	 0.604374640201265e-03	},
		{  14,	-12,	-0.225132933900136e-15	},
		{  14,	-10,	 0.610916973582981e-11	},
		{  14,	 -8,	-0.303063908043404e-06	},
		{  14,	 -6,	-0.137796070798409e-04	},
		{  14,	 -4,	-0.919296736666106e-03	},
		{  16,	-10,	 0.639288223132545e-09	},
		{  16,	 -8,	 0.753259479898699e-06	},
		{  18,	-12,	-0.400321478682929e-12	},
		{  18,	-10,	 0.756140294351614e-08	},
		{  20,	-12,	-0.912082054034891e-11	},
		{  20,	-10,	-0.237612381140539e-07	},
		{  20,	 -6,	 0.269586010591874e-04	},
		{  22,	-12,	-0.732828135157839e-10	},
		{  24,	-12,	 0.241995578306660e-09	},
		{  24,	 -4,	-0.405735532730322e-03	},
		{  28,	-12,	 0.189424143498011e-09	},
		{  32,	-12,	-0.486632965074563e-09	},
	};
	const if97_coef coefg[SIZEG] = {
		{ -12,	  7,	 0.412209020652996e-04	},
		{ -12,	 12,	-0.114987238280587e+07	},
		{ -12,	 14,	 0.948180885032080e+10	},
		{ -12,	 18,	-0.195788865718971e+18	},
		{ -12,	 22,	 0.496250704871300e+25	},
		{ -12,	 24,	-0.105549884548496e+29	},
		{ -10,	 14,	-0.758642165988278e+12	},
		{ -10,	 20,	-0.922172769596101e+23	},
		{ -10,	 24,	 0.725379072059348e+30	},
		{  -8,	  7,	-0.617718249205859e+02	},
		{  -8,	  8,	 0.107555033344858e+05	},
		{  -8,	 10,	-0.379545802336487e+08	},
		{  -8,	 12,	 0.228646846221831e+12	},
		{  -6,	  8,	-0.499741093010619e+07	},
		{  -6,	 22,	-0.280214310054101e+31	},
		{  -5,	  7,	 0.104915406769586e+07	},
		{  -5,	 20,	 0.613754229168619e+28	},
		{  -4,	 22,	 0.802056715528378e+32	},
		{  -3,	  7,	-0.298617819828065e+08	},
		{  -2,	  3,	-0.910782540134681e+02	},
		{  -2,	  5,	 0.135033227281565e+06	},
		{  -2,	 14,	-0.712949383408211e+19	},
		{  -2,	 24,	-0.104578785289542e+37	},
		{  -1,	  2,	 0.304331584444093e+02	},
		{  -1,	  8,	 0.593250797959445e+10	},
		{  -1,	 18,	-0.364174062110798e+28	},
		{   0,	  0,	 0.921791403532461e+00	},
		{   0,	  1,	-0.337693609657471e+00	},
		{   0,	  2,	-0.724644143758508e+02	},
		{   1,	  0,	-0.110480239272601e+00	},
		{   1,	  1,	 0.536516031875059e+01	},
		{   1,	  3,	-0.291441872156205e+04	},
		{   3,	 24,	 0.616338176535305e+40	},
		{   5,	 22,	-0.120889175861180e+39	},
		{   6,	 12,	 0.818396024524612e+23	},
		{   8,	  3,	 0.940781944835829e+09	},
		{  10,	  0,	-0.367279669545448e+05	},
		{  10,	  6,	-0.837513931798655e+16	},
	};
	const if97_coef coefh[SIZEH] = {
		{ -12,	  8,	 0.561379678887577e-01	},
		{ -12,	 12,	 0.774135421587083e+10	},
		{ -10,	  4,	 0.111482975877938e-08	},
		{ -10,	  6,	-0.143987128208183e-02	},
		{ -10,	  8,	 0.193696558764920e+04	},
		{ -10,	 10,	-0.605971823585005e+09	},
		{ -10,	 14,	 0.171951568124337e+14	},
		{ -10,	 16,	-0.185461154985145e+17	},
		{  -8,	  0,	 0.387851168078010e-16	},
		{  -8,	  1,	-0.395464327846105e-13	},
		{  -8,	  6,	-0.170875935679023e+03	},
		{  -8,	  7,	-0.212010620701220e+04	},
		{  -8,	  8,	 0.177683337348191e+08	},
		{  -6,	  4,	 0.110177443629575e+02	},
		{  -6,	  6,	-0.234396091693313e+06	},
		{  -6,	  8,	-0.656174421999594e+07	},
		{  -5,	  2,	 0.156362212977396e-04	},
		{  -5,	  3,	-0.212946257021400e+01	},
		{  -5,	  4,	 0.135249306374858e+02	},
		{  -4,	  2,	 0.177189164145813e+00	},
		{  -4,	  4,	 0.139499167345464e+04	},
		{  -3,	  1,	-0.703670932036388e-02	},
		{  -3,	  2,	-0.152011044389648e+00	},
		{  -2,	  0,	 0.981916922991113e-04	},
		{  -1,	  0,	 0.147199658618076e-02	},
		{  -1,	  2,	 0.202618487025578e+02	},
		{   0,	  0,	 0.899345518944240e+00	},
		{   1,	  0,	-0.211346402240858e+00	},
		{   1,	  2,	 0.249971752957491e+02	},
	};
	const if97_coef coefi[SIZEI] = {
		{   0,	  0,	 0.106905684359136e+01	},
		{   0,	  1,	-0.148620857922333e+01	},
		{   0,	 10,	 0.259862256980408e+15	},
		{   1,	 -4,	-0.446352055678749e-11	},
		{   1,	 -2,	-0.566620757170032e-06	},
		{   1,	 -1,	-0.235302885736849e-02	},
		{   1,	  0,	-0.269226321968839e+00	},
		{   2,	  0,	 0.922024992944392e+01	},
		{   3,	 -5,	 0.357633505503772e-11	},
		{   3,	  0,	-0.173942565562222e+02	},
		{   4,	 -3,	 0.700681785556229e-05	},
		{   4,	 -2,	-0.267050351075768e-03	},
		{   4,	 -1,	-0.231779669675624e+01	},
		{   5,	 -6,	-0.753533046979752e-12	},
		{   5,	 -1,	 0.481337131452891e+01	},
		{   5,	 12,	-0.223286270422356e+22	},
		{   7,	 -4,	-0.118746004987383e-04	},
		{   7,	 -3,	 0.646412934136496e-02	},
		{   8,	 -6,	-0.410588536330937e-09	},
		{   8,	 10,	 0.422739537057241e+20	},
		{  10,	 -8,	 0.313698180473812e-12	},
		{  12,	-12,	 0.164395334345040e-23	},
		{  12,	 -6,	-0.339823323754373e-05	},
		{  12,	 -4,	-0.135268639905021e-01	},
		{  14,	-10,	-0.723252514211625e-14	},
		{  14,	 -8,	 0.184386437538366e-08	},
		{  14,	 -4,	-0.463959533752385e-01	},
		{  14,	  5,	-0.992263100376750e+14	},
		{  18,	-12,	 0.688169154439335e-16	},
		{  18,	-10,	-0.222620998452197e-10	},
		{  18,	 -8,	-0.540843018624083e-07	},
		{  18,	 -6,	 0.345570606200257e-02	},
		{  18,	  2,	 0.422275800304086e+11	},
		{  20,	-12,	-0.126974478770487e-14	},
		{  20,	-10,	 0.927237985153679e-09	},
		{  22,	-12,	 0.612670812016489e-13	},
		{  24,	-12,	-0.722693924063497e-11	},
		{  24,	 -8,	-0.383669502636822e-03	},
		{  32,	-10,	 0.374684572410204e-03	},
		{  32,	 -5,	-0.931976897511086e+05	},
		{  36,	-10,	-0.247690616026922e-01	},
		{  36,	 -8,	 0.658110546759474e+02	},
	};
	const if97_coef coefj[SIZEJ] = {
		{   0,	 -1,	-0.111371317395540e-03	},
		{   0,	  0,	 0.100342892423685e+01	},
		{   0,	  1,	 0.530615581928979e+01	},
		{   1,	 -2,	 0.179058760078792e-05	},
		{   1,	 -1,	-0.728541958464774e-03	},
		{   1,	  1,	-0.187576133371704e+02	},
		{   2,	 -1,	 0.199060874071849e-02	},
		{   2,	  1,	 0.243574755377290e+02	},
		{   3,	 -2,	-0.177040785499444e-03	},
		{   4,	 -2,	-0.259680385227130e-02	},
		{   4,	  2,	-0.198704578406823e+03	},
		{   5,	 -3,	 0.738627790224287e-04	},
		{   5,	 -2,	-0.236264692844138e-02	},
		{   5,	  0,	-0.161023121314333e+01	},
		{   6,	  3,	 0.622322971786473e+04	},
		{  10,	 -6,	-0.960754116701669e-08	},
		{  12,	 -8,	-0.510572269720488e-10	},
		{  12,	 -3,	 0.767373781404211e-02	},
		{  14,	-10,	 0.663855469485254e-14	},
		{  14,	 -8,	-0.717590735526745e-09	},
		{  14,	 -5,	 0.146564542926508e-04	},
		{  16,	-10,	 0.309029474277013e-11	},
		{  18,	-12,	-0.464216300971708e-15	},
		{  20,	-12,	-0.390499637961161e-13	},
		{  20,	-10,	-0.236716126781431e-09	},
		{  24,	-12,	 0.454652854268717e-11	},
		{  24,	 -6,	-0.422271787482497e-02	},
		{  28,	-12,	 0.283911742354706e-10	},
		{  28,	 -5,	 0.270929002720228e+01	},
	};
	const if97_coef coefk[SIZEK] = {
		{  -2,	 10,	-0.401215699576099e+09	},
		{  -2,	 12,	 0.484501478318406e+11	},
		{  -1,	 -5,	 0.394721471363678e-14	},
		{  -1,	  6,	 0.372629967374147e+05	},
		{   0,	-12,	-0.369794374168666e-29	},
		{   0,	 -6,	-0.380436407012452e-14	},
		{   0,	 -2,	 0.475361629970233e-06	},
		{   0,	 -1,	-0.879148916140706e-03	},
		{   0,	  0,	 0.844317863844331e+00	},
		{   0,	  1,	 0.122433162656600e+02	},
		{   0,	  2,	-0.104529634830279e+03	},
		{   0,	  3,	 0.589702771277429e+03	},
		{   0,	 14,	-0.291026851164444e+14	},
		{   1,	 -3,	 0.170343072841850e-05	},
		{   1,	 -2,	-0.277617606975748e-03	},
		{   1,	  0,	-0.344709605486686e+01	},
		{   1,	  1,	 0.221333862447095e+02	},
		{   1,	  2,	-0.194646110037079e+03	},
		{   2,	 -8,	 0.808354639772825e-15	},
		{   2,	 -6,	-0.180845209145470e-10	},
		{   2,	 -3,	-0.696664158132412e-05	},
		{   2,	 -2,	-0.181057560300994e-02	},
		{   2,	  0,	 0.255830298579027e+01	},
		{   2,	  4,	 0.328913873658481e+04	},
		{   5,	-12,	-0.173270241249904e-18	},
		{   5,	 -6,	-0.661876792558034e-06	},
		{   5,	 -3,	-0.395688923421250e-02	},
		{   6,	-12,	 0.604203299819132e-17	},
		{   6,	-10,	-0.400879935920517e-13	},
		{   6,	 -8,	 0.160751107464958e-08	},
		{   6,	 -5,	 0.383719409025556e-04	},
		{   8,	-12,	-0.649565446702457e-14	},
		{  10,	-12,	-0.149095328506000e-11	},
		{  12,	-10,	 0.541449377329581e-08	},
	};
	const if97_coef coefl[SIZEL] = {
		{ -12,	 14,	 0.260702058647537e+10	},
		{ -12,	 16,	-0.188277213604704e+15	},
		{ -12,	 18,	 0.554923870289667e+19	},
		{ -12,	 20,	-0.758966946387758e+23	},
		{ -12,	 22,	 0.413865186848908e+27	},
		{ -10,	 14,	-0.815038000738060e+12	},
		{ -10,	 24,	-0.381458260489955e+33	},
		{  -8,	  6,	-0.123239564600519e-01	},
		{  -8,	 10,	 0.226095631437174e+08	},
		{  -8,	 12,	-0.495017809506720e+12	},
		{  -8,	 14,	 0.529482996422863e+16	},
		{  -8,	 18,	-0.444359478746295e+23	},
		{  -8,	 24,	 0.521635864527315e+35	},
		{  -8,	 36,	-0.487095672740742e+55	},
		{  -6,	  8,	-0.714430209937547e+06	},
		{  -5,	  4,	 0.127868634615495e+00	},
		{  -5,	  5,	-0.100752127917598e+02	},
		{  -4,	  7,	 0.777451437960990e+07	},
		{  -4,	 16,	-0.108105480796471e+25	},
		{  -3,	  1,	-0.357578581169659e-05	},
		{  -3,	  3,	-0.212857169423484e+01	},
		{  -3,	 18,	 0.270706111085238e+30	},
		{  -3,	 20,	-0.695953622348829e+33	},
		{  -2,	  2,	 0.110609027472280e+00	},
		{  -2,	  3,	 0.721559163361354e+02	},
		{  -2,	 10,	-0.306367307532219e+15	},
		{  -1,	  0,	 0.265839618885530e-04	},
		{  -1,	  1,	 0.253392392889754e-01	},
		{  -1,	  3,	-0.214443041836579e+03	},
		{   0,	  0,	 0.937846601489667e+00	},
		{   0,	  1,	 0.223184043101700e+01	},
		{   0,	  2,	 0.338401222509191e+02	},
		{   0,	 12,	 0.494237237179718e+21	},
		{   1,	  0,	-0.198068404154428e+00	},
		{   1,	 16,	-0.141415349881140e+31	},
		{   2,	  1,	-0.993862421613651e+02	},
		{   4,	  0,	 0.125070534142731e+03	},
		{   5,	  0,	-0.996473529004439e+03	},
		{   5,	  1,	 0.473137909872765e+05	},
		{   6,	 14,	 0.116662121219322e+33	},
		{  10,	  4,	-0.315874976271533e+16	},
		{  10,	 12,	-0.445703369196945e+33	},
		{  14,	 10,	 0.642794932373694e+33	},
	};
	const if97_coef coefm[SIZEM] = {
		{   0,	  0,	 0.811384363481847e+00	},
		{   0,	 14,	 0.185135446828337e+09	},
		{   0,	 24,	-0.821698160721956e+15	},
		{   0,	 28,	 0.233415869478510e+18	},
		{   0,	 36,	-0.329421923951460e+22	},
		{   1,	  5,	-0.814568209346872e+05	},
		{   1,	  6,	 0.458384828593949e+06	},
		{   1,	 14,	-0.170451090076385e+12	},
		{   1,	 18,	 0.157890366037614e+15	},
		{   1,	 20,	-0.202530509748774e+16	},
		{   1,	 32,	 0.188813911076809e+22	},
		{   2,	  7,	 0.453735800004273e+08	},
		{   2,	 10,	-0.547578313899097e+10	},
		{   2,	 22,	 0.170215539458936e+18	},
		{   2,	 36,	-0.137570282536696e+26	},
		{   3,	  0,	-0.568199310990094e+04	},
		{   3,	  5,	-0.659774567602874e+08	},
		{   3,	 12,	 0.185007245563239e+13	},
		{   3,	 28,	-0.600079934586803e+23	},
		{   3,	 36,	 0.181508996303902e+28	},
		{   4,	  5,	-0.152861148659302e+11	},
		{   4,	  8,	 0.939454935735563e+12	},
		{   4,	 28,	 0.594584382273384e+25	},
		{   4,	 36,	-0.346865122768353e+30	},
		{   5,	  5,	-0.560165667510446e+12	},
		{   5,	 10,	 0.200725701112386e+15	},
		{   5,	 24,	-0.795260241872306e+24	},
		{   6,	  6,	-0.385754000383848e+14	},
		{   8,	  0,	-0.178657198172556e+11	},
		{   8,	 32,	 0.111052244098768e+36	},
		{   8,	 36,	-0.211961148774260e+38	},
		{  12,	 28,	 0.189461279349492e+40	},
		{  14,	  8,	 0.266572856432938e+28	},
		{  14,	 32,	 0.291133958602503e+46	},
		{  14,	 36,	-0.128617899887675e+49	},
		{  16,	 22,	 0.639234909918741e+42	},
		{  16,	 28,	-0.810093428842645e+46	},
		{  20,	  2,	 0.795537657613427e+32	},
		{  24,	 36,	 0.479817895699239e+65	},
		{  28,	 20,	 0.368193926183570e+60	},
	};
	const if97_coef coefn[SIZEN] = {
		{   0,	-12,	 0.280967799943151e-38	},
		{   0,	-10,	-0.135031446451331e-31	},
		{   0,	 -8,	 0.177274872361946e-25	},
		{   0,	 -2,	 0.240560808321713e-06	},
		{   0,	 -1,	-0.890763306701305e-03	},
		{   0,	  1,	 0.159158748314599e+04	},
		{   0,	  2,	-0.792681207132600e+06	},
		{   0,	  5,	 0.354542769185671e+12	},
		{   1,	  0,	-0.302807107747776e+03	},
		{   1,	  1,	 0.232534272709876e+06	},
		{   1,	  4,	-0.869871364662769e+11	},
		{   1,	  6,	 0.400849240129329e+15	},
		{   2,	 -6,	 0.541276911564176e-13	},
		{   2,	 -5,	-0.493111362030162e-10	},
		{   3,	-12,	 0.614869006573609e-30	},
		{   3,	-10,	-0.607246643970893e-23	},
		{   3,	 -8,	-0.334952758812999e-18	},
		{   3,	 -6,	 0.705412100773699e-11	},
		{   3,	 -3,	-0.643064132636925e-02	},
		{   3,	 -1,	-0.440209599407714e+04	},
		{   4,	-12,	 0.582238667048942e-27	},
		{   4,	 -6,	 0.258585887897486e-08	},
		{   4,	 -5,	-0.158649699894543e-05	},
		{   4,	 -4,	 0.220019901729615e-02	},
		{   5,	-10,	-0.402352115234494e-18	},
		{   5,	 -3,	 0.629154149015048e+02	},
		{   6,	-12,	 0.390628369238462e-22	},
		{   6,	-10,	-0.744938506925544e-16	},
		{   6,	 -3,	 0.135147318617061e+03	},
		{   7,	-12,	 0.821445758255119e-20	},
		{   7,	 -8,	-0.421537726098389e-08	},
		{   7,	 -5,	-0.525037427886100e+00	},
		{   8,	-10,	 0.189917206526237e-12	},
		{  10,	-12,	 0.402137961842776e-14	},
		{  12,	-12,	 0.651718171878301e-12	},
		{  12,	-10,	 0.364975183508473e-05	},
		{  12,	 -8,	-0.391048167929649e-01	},
		{  14,	-12,	-0.211773355803058e-07	},
		{  18,	-12,	 0.264953354380072e-02	},
	};
	const if97_coef coefo[SIZEO] = {
		{   0,	-12,	 0.128746023979718e-34	},
		{   0,	 -4,	-0.735234770382342e-11	},
		{   0,	 -1,	 0.289078692149150e-02	},
		{   2,	 -1,	 0.244482731907223e+00	},
		{   3,	-10,	 0.141733492030985e-23	},
		{   4,	-12,	-0.354533853059476e-28	},
		{   4,	 -8,	-0.594539202901431e-17	},
		{   4,	 -5,	-0.585188401782779e-08	},
		{   4,	 -4,	 0.201377325411803e-05	},
		{   4,	 -1,	 0.138647388209306e+01	},
		{   5,	 -4,	-0.173959365084772e-04	},
		{   5,	 -3,	 0.137680878349369e-02	},
		{   6,	 -8,	 0.814897605805513e-14	},
		{   7,	-12,	 0.425596631351839e-25	},
		{   8,	-10,	-0.387449113787755e-17	},
		{   8,	 -8,	 0.139814747930240e-12	},
		{   8,	 -4,	-0.171849638951521e-02	},
		{  10,	-12,	 0.641890529513296e-21	},
		{  10,	 -8,	 0.118960578072018e-10	},
		{  14,	-12,	-0.155282762571611e-17	},
		{  14,	 -8,	 0.233907907347507e-07	},
		{  20,	-12,	-0.174093247766213e-12	},
		{  20,	-10,	 0.377682649089149e-08	},
		{  24,	-12,	-0.516720236575302e-10	},
	};
	const if97_coef coefp[SIZEP] = {
		{   0,	 -1,	-0.982825342010366e-04	},
		{   0,	  0,	 0.105145700850612e+01	},
		{   0,	  1,	 0.116033094095084e+03	},
		{   0,	  2,	 0.324664750281543e+04	},
		{   1,	  1,	-0.123592348610137e+04	},
		{   2,	 -1,	-0.561403450013495e-01	},
		{   3,	 -3,	 0.856677401640869e-07	},
		{   3,	  0,	 0.236313425393924e+03	},
		{   4,	 -2,	 0.972503292350109e-02	},
		{   6,	 -2,	-0.103001994531927e+01	},
		{   7,	 -5,	-0.149653706199162e-08	},
		{   7,	 -4,	-0.215743778861592e-04	},
		{   8,	 -2,	-0.834452198291445e+01	},
		{  10,	 -3,	 0.586602660564988e+00	},
		{  12,	-12,	 0.343480022104968e-25	},
		{  12,	 -6,	 0.816256095947021e-05	},
		{  12,	 -5,	 0.294985697916798e-02	},
		{  14,	-10,	 0.711730466276584e-16	},
		{  14,	 -8,	 0.400954763806941e-09	},
		{  14,	 -3,	 0.107766027032853e+02	},
		{  16,	 -8,	-0.409449599138182e-06	},
		{  18,	 -8,	-0.729121307758902e-05	},
		{  20,	-10,	 0.677107970938909e-08	},
		{  22,	-10,	 0.602745973022975e-07	},
		{  24,	-12,	-0.382323011855257e-10	},
		{  24,	 -8,	 0.179946628317437e-02	},
		{  36,	-12,	-0.345042834640005e-03	},
	};
	const if97_coef coefq[SIZEQ] = {
		{ -12,	 10,	-0.820433843259950e+05	},
		{ -12,	 12,	 0.473271518461586e+11	},
		{ -10,	  6,	-0.805950021005413e-01	},
		{ -10,	  7,	 0.328600025435980e+02	},
		{ -10,	  8,	-0.356617029982490e+04	},
		{ -10,	 10,	-0.172985781433335e+10	},
		{  -8,	  8,	 0.351769232729192e+08	},
		{  -6,	  6,	-0.775489259985144e+06	},
		{  -5,	  2,	 0.710346691966018e-04	},
		{  -5,	  5,	 0.993499883820274e+05	},
		{  -4,	  3,	-0.642094171904570e+00	},
		{  -4,	  4,	-0.612842816820083e+04	},
		{  -3,	  3,	 0.232808472983776e+03	},
		{  -2,	  0,	-0.142808220416837e-04	},
		{  -2,	  1,	-0.643596060678456e-02	},
		{  -2,	  2,	-0.428577227475614e+01	},
		{  -2,	  4,	 0.225689939161918e+04	},
		{  -1,	  0,	 0.100355651721510e-02	},
		{  -1,	  1,	 0.333491455143516e+00	},
		{  -1,	  2,	 0.109697576888873e+01	},
		{   0,	  0,	 0.961917379376452e+00	},
		{   1,	  0,	-0.838165632204598e-01	},
		{   1,	  1,	 0.247795908411492e+01	},
		{   1,	  3,	-0.319114969006533e+04	},
	};
	const if97_coef coefr[SIZER] = {
		{  -8,	  6,	 0.144165955660863e-02	},
		{  -8,	 14,	-0.701438599628258e+13	},
		{  -3,	 -3,	-0.830946716459219e-16	},
		{  -3,	  3,	 0.261975135368109e+00	},
		{  -3,	  4,	 0.393097214706245e+03	},
		{  -3,	  5,	-0.104334030654021e+05	},
		{  -3,	  8,	 0.490112654154211e+09	},
		{   0,	 -1,	-0.147104222772069e-03	},
		{   0,	  0,	 0.103602748043408e+01	},
		{   0,	  1,	 0.305308890065089e+01	},
		{   0,	  5,	-0.399745276971264e+07	},
		{   3,	 -6,	 0.569233719593750e-11	},
		{   3,	 -2,	-0.464923504407778e-01	},
		{   8,	-12,	-0.535400396512906e-17	},
		{   8,	-10,	 0.399988795693162e-12	},
		{   8,	 -8,	-0.536479560201811e-06	},
		{   8,	 -5,	 0.159536722411202e-01	},
		{  10,	-12,	 0.270303248860217e-14	},
		{  10,	-10,	 0.244247453858506e-07	},
		{  10,	 -8,	-0.983430636716454e-05	},
		{  10,	 -6,	 0.663513144224454e-01	},
		{  10,	 -5,	-0.993456957845006e+01	},
		{  10,	 -4,	 0.546491323528491e+03	},
		{  10,	 -3,	-0.143365406393758e+05	},
		{  10,	 -2,	 0.150764974125511e+06	},
		{  12,	-12,	-0.337209709340105e-09	},
		{  14,	-12,	 0.377501980025469e-08	},
	};
	const if97_coef coefs[SIZES] = {
		{ -12,	 20,	-0.532466612140254e+23	},
		{ -12,	 24,	 0.100415480000824e+32	},
		{ -10,	 22,	-0.191540001821367e+30	},
		{  -8,	 14,	 0.105618377808847e+17	},
		{  -6,	 36,	 0.202281884477061e+59	},
		{  -5,	  8,	 0.884585472596134e+08	},
		{  -5,	 16,	 0.166540181638363e+23	},
		{  -4,	  6,	-0.313563197669111e+06	},
		{  -4,	 32,	-0.185662327545324e+54	},
		{  -3,	  3,	-0.624942093918942e-01	},
		{  -3,	  8,	-0.504160724132590e+10	},
		{  -2,	  4,	 0.187514491833092e+05	},
		{  -1,	  1,	 0.121399979993217e-02	},
		{  -1,	  2,	 0.188317043049455e+01	},
		{  -1,	  3,	-0.167073503962060e+04	},
		{   0,	  0,	 0.965961650599775e+00	},
		{   0,	  1,	 0.294885696802488e+01	},
		{   0,	  4,	-0.653915627346115e+05	},
		{   0,	 28,	 0.604012200163444e+50	},
		{   1,	  0,	-0.198339358557937e+00	},
		{   1,	 32,	-0.175984090163501e+58	},
		{   3,	  0,	 0.356314881403987e+01	},
		{   3,	  1,	-0.575991255144384e+03	},
		{   3,	  2,	 0.456213415338071e+05	},
		{   4,	  3,	-0.109174044987829e+08	},
		{   4,	 18,	 0.437796099975134e+34	},
		{   4,	 24,	-0.616552611135792e+46	},
		{   5,	  4,	 0.193568768917797e+10	},
		{  14,	 24,	 0.950898170425042e+54	},
	};
	const if97_coef coeft[SIZET] = {
		{   0,	  0,	 0.155287249586268e+01	},
		{   0,	  1,	 0.664235115009031e+01	},
		{   0,	  4,	-0.289366236727210e+04	},
		{   0,	 12,	-0.385923202309848e+13	},
		{   1,	  0,	-0.291002915783761e+01	},
		{   1,	 10,	-0.829088246858083e+12	},
		{   2,	  0,	 0.176814899675218e+01	},
		{   2,	  6,	-0.534686695713469e+09	},
		{   2,	 14,	 0.160464608687834e+18	},
		{   3,	  3,	 0.196435366560186e+06	},
		{   3,	  8,	 0.156637427541729e+13	},
		{   4,	  0,	-0.178154560260006e+01	},
		{   4,	 10,	-0.229746237623692e+16	},
		{   7,	  3,	 0.385659001648006e+08	},
		{   7,	  4,	 0.110554446790543e+10	},
		{   7,	  7,	-0.677073830687349e+14	},
		{   7,	 20,	-0.327910592086523e+31	},
		{   7,	 36,	-0.341552040860644e+51	},
		{  10,	 10,	-0.527251339709047e+21	},
		{  10,	 12,	 0.245375640937055e+24	},
		{  10,	 14,	-0.168776617209269e+27	},
		{  10,	 16,	 0.358958955867578e+29	},
		{  10,	 22,	-0.656475280339411e+36	},
		{  18,	 18,	 0.355286045512301e+39	},
		{  20,	 32,	 0.569021454413270e+58	},
		{  22,	 22,	-0.700584546433113e+48	},
		{  22,	 36,	-0.705772623326374e+65	},
		{  24,	 24,	 0.166861176200148e+53	},
		{  28,	 28,	-0.300475129680486e+61	},
		{  32,	 22,	-0.668481295196808e+51	},
		{  32,	 32,	 0.428432338620678e+69	},
		{  32,	 36,	-0.444227367758304e+72	},
		{  36,	 36,	-0.281396013562745e+77	},
	};
	const if97_coef coefu[SIZEU] = {
		{ -12,	 14,	 0.122088349258355e+18	},
		{ -10,	 10,	 0.104216468608488e+10	},
		{ -10,	 12,	-0.882666931564652e+16	},
		{ -10,	 14,	 0.259929510849499e+20	},
		{  -8,	 10,	 0.222612779142211e+15	},
		{  -8,	 12,	-0.878473585050085e+18	},
		{  -8,	 14,	-0.314432577551552e+22	},
		{  -6,	  8,	-0.216934916996285e+13	},
		{  -6,	 12,	 0.159079648196849e+21	},
		{  -5,	  4,	-0.339567617303423e+03	},
		{  -5,	  8,	 0.884387651337836e+13	},
		{  -5,	 12,	-0.843405926846418e+21	},
		{  -3,	  2,	 0.114178193518022e+02	},
		{  -1,	 -1,	-0.122708229235641e-03	},
		{  -1,	  1,	-0.106201671767107e+03	},
		{  -1,	 12,	 0.903443213959313e+25	},
		{  -1,	 14,	-0.693996270370852e+28	},
		{   0,	 -3,	 0.648916718965575e-08	},
		{   0,	  1,	 0.718957567127851e+04	},
		{   1,	 -2,	 0.105581745346187e-02	},
		{   2,	  5,	-0.651903203602581e+15	},
		{   2,	 10,	-0.160116813274676e+25	},
		{   3,	 -5,	-0.510254294237837e-08	},
		{   5,	 -4,	-0.152355388953402e+00	},
		{   5,	  2,	 0.677143292290144e+12	},
		{   5,	  3,	 0.276378438378930e+15	},
		{   6,	 -5,	 0.116862983141686e-01	},
		{   6,	  2,	-0.301426947980171e+14	},
		{   8,	 -8,	 0.169719813884840e-07	},
		{   8,	  8,	 0.104674840020929e+27	},
		{  10,	 -4,	-0.108016904560140e+05	},
		{  12,	-12,	-0.990623601934295e-12	},
		{  12,	 -4,	 0.536116483602738e+07	},
		{  12,	  4,	 0.226145963747881e+22	},
		{  14,	-12,	-0.488731565776210e-09	},
		{  14,	-10,	 0.151001548880670e-04	},
		{  14,	 -6,	-0.227700464643920e+05	},
		{  14,	  6,	-0.781754507698846e+28	},
	};
	const if97_coef coefv[SIZEV] = {
		{ -10,	 -8,	-0.415652812061591e-54	},
		{  -8,	-12,	 0.177441742924043e-60	},
		{  -6,	-12,	-0.357078668203377e-54	},
		{  -6,	 -3,	 0.359252213604114e-25	},
		{  -6,	  5,	-0.259123736380269e+02	},
		{  -6,	  6,	 0.594619766193460e+05	},
		{  -6,	  8,	-0.624184007103158e+11	},
		{  -6,	 10,	 0.313080299915944e+17	},
		{  -5,	  1,	 0.105006446192036e-08	},
		{  -5,	  2,	-0.192824336984852e-05	},
		{  -5,	  6,	 0.654144373749937e+06	},
		{  -5,	  8,	 0.513117462865044e+13	},
		{  -5,	 10,	-0.697595750347391e+19	},
		{  -5,	 14,	-0.103977184454767e+29	},
		{  -4,	-12,	 0.119563135540666e-47	},
		{  -4,	-10,	-0.436677034051655e-41	},
		{  -4,	 -6,	 0.926990036530639e-29	},
		{  -4,	 10,	 0.587793105620748e+21	},
		{  -3,	 -3,	 0.280375725094731e-17	},
		{  -3,	 10,	-0.192359972440634e+23	},
		{  -3,	 12,	 0.742705723302738e+27	},
		{  -2,	  2,	-0.517429682450605e+02	},
		{  -2,	  4,	 0.820612048645469e+07	},
		{  -1,	 -2,	-0.188214882341448e-08	},
		{  -1,	  0,	 0.184587261114837e-01	},
		{   0,	 -2,	-0.135830407782663e-05	},
		{   0,	  6,	-0.723681885626348e+17	},
		{   0,	 10,	-0.223449194054124e+27	},
		{   1,	-12,	-0.111526741826431e-34	},
		{   1,	-10,	 0.276032601145151e-28	},
		{   3,	  3,	 0.134856491567853e+15	},
		{   4,	 -6,	 0.652440293345860e-09	},
		{   4,	  3,	 0.510655119774360e+17	},
		{   4,	 10,	-0.468138358908732e+32	},
		{   5,	  2,	-0.760667491183279e+16	},
		{   8,	-12,	-0.417247986986821e-18	},
		{  10,	 -2,	 0.312545677756104e+14	},
		{  12,	 -3,	-0.100375333864186e+15	},
		{  14,	  1,	 0.247761392329058e+27	},
	};
	const if97_coef coefw[SIZEW] = {
		{ -12,	  8,	-0.586219133817016e-07	},
		{ -12,	 14,	-0.894460355005526e+11	},
		{ -10,	 -1,	 0.531168037519774e-30	},
		{ -10,	  8,	 0.109892402329239e+00	},
		{  -8,	  6,	-0.575368389425212e-01	},
		{  -8,	  8,	 0.228276853990249e+05	},
		{  -8,	 14,	-0.158548609655002e+19	},
		{  -6,	 -4,	 0.329865748576503e-27	},
		{  -6,	 -3,	-0.634987981190669e-24	},
		{  -6,	  2,	 0.615762068640611e-08	},
		{  -6,	  8,	-0.961109240985747e+08	},
		{  -5,	-10,	-0.406274286652625e-44	},
		{  -4,	 -1,	-0.471103725498077e-12	},
		{  -4,	  3,	 0.725937724828145e+00	},
		{  -3,	-10,	 0.187768525763682e-38	},
		{  -3,	  3,	-0.103308436323771e+04	},
		{  -2,	  1,	-0.662552816342168e-01	},
		{  -2,	  2,	 0.579514041765710e+03	},
		{  -1,	 -8,	 0.237416732616644e-26	},
		{  -1,	 -4,	 0.271700235739893e-14	},
		{  -1,	  1,	-0.907886213483600e+02	},
		{   0,	-12,	-0.171242509570207e-36	},
		{   0,	  1,	 0.156792067854621e+03	},
		{   1,	 -1,	 0.923261357901470e+00	},
		{   2,	 -1,	-0.597865988422577e+01	},
		{   2,	  2,	 0.321988767636389e+07	},
		{   3,	-12,	-0.399441390042203e-29	},
		{   3,	 -5,	 0.493429086046981e-07	},
		{   5,	-10,	 0.812036983370565e-19	},
		{   5,	 -8,	-0.207610284654137e-11	},
		{   5,	 -6,	-0.340821291419719e-06	},
		{   8,	-12,	 0.542000573372233e-17	},
		{   8,	-10,	-0.856711586510214e-12	},
		{  10,	-12,	 0.266170454405981e-13	},
		{  10,	 -8,	 0.858133791857099e-05	},
	};
	const if97_coef coefx[SIZEX] = {
		{  -8,	 14,	 0.377373741298151e+19	},
		{  -6,	 10,	-0.507100883722913e+13	},
		{  -5,	 10,	-0.103363225598860e+16	},
		{  -4,	  1,	 0.184790814320773e-05	},
		{  -4,	  2,	-0.924729378390945e-03	},
		{  -4,	 14,	-0.425999562292738e+24	},
		{  -3,	 -2,	-0.462307771873973e-12	},
		{  -3,	 12,	 0.107319065855767e+22	},
		{  -1,	  5,	 0.648662492280682e+11	},
		{   0,	  0,	 0.244200600688281e+01	},
		{   0,	  4,	-0.851535733484258e+10	},
		{   0,	 10,	 0.169894481433592e+22	},
		{   1,	-10,	 0.215780222509020e-26	},
		{   1,	 -1,	-0.320850551367334e+00	},
		{   2,	  6,	-0.382642448458610e+17	},
		{   3,	-12,	-0.275386077674421e-28	},
		{   3,	  0,	-0.563199253391666e+06	},
		{   3,	  8,	-0.326068646279314e+21	},
		{   4,	  3,	 0.397949001553184e+14	},
		{   5,	 -6,	 0.100824008584757e-06	},
		{   5,	 -2,	 0.162234569738433e+05	},
		{   5,	  1,	-0.432355225319745e+11	},
		{   6,	  1,	-0.592874245598610e+12	},
		{   8,	 -6,	 0.133061647281106e+01	},
		{   8,	 -3,	 0.157338197797544e+07	},
		{   8,	  1,	 0.258189614270853e+14	},
		{   8,	  8,	 0.262413209706358e+25	},
		{  10,	 -8,	-0.920011937431142e-01	},
		{  12,	-10,	 0.220213765905426e-02	},
		{  12,	 -8,	-0.110433759109547e+02	},
		{  12,	 -5,	 0.847004870612087e+07	},
		{  12,	 -4,	-0.592910695762536e+09	},
		{  14,	-12,	-0.183027173269660e-04	},
		{  14,	-10,	 0.181339603516302e+00	},
		{  14,	 -8,	-0.119228759669889e+04	},
		{  14,	 -6,	 0.430867658061468e+07	},
	};
	const if97_coef coefy[SIZEY] = {
		{   0,	 -3,	-0.525597995024633e-09	},
		{   0,	  1,	 0.583441305228407e+04	},
		{   0,	  5,	-0.134778968457925e+17	},
		{   0,	  8,	 0.118973500934212e+26	},
		{   1,	  8,	-0.159096490904708e+27	},
		{   2,	 -4,	-0.315839902302021e-06	},
		{   2,	 -1,	 0.496212197158239e+03	},
		{   2,	  4,	 0.327777227273171e+19	},
		{   2,	  5,	-0.527114657850696e+22	},
		{   3,	 -8,	 0.210017506281863e-16	},
		{   3,	  4,	 0.705106224399834e+21	},
		{   3,	  8,	-0.266713136106469e+31	},
		{   4,	 -6,	-0.145370512554562e-07	},
		{   4,	  6,	 0.149333917053130e+28	},
		{   5,	 -2,	-0.149795620287641e+08	},
		{   5,	  1,	-0.381881906271100e+16	},
		{   8,	 -8,	 0.724660165585797e-04	},
		{   8,	 -2,	-0.937808169550193e+14	},
		{  10,	 -5,	 0.514411468376383e+10	},
		{  12,	 -8,	-0.828198594040141e+05	},
	};
	const if97_coef coefz[SIZEZ] = {
		{  -8,	  3,	 0.244007892290650e-10	},
		{  -6,	  6,	-0.463057430331240e+07	},
		{  -5,	  6,	 0.728803274777710e+10	},
		{  -5,	  8,	 0.327776302858860e+16	},
		{  -4,	  5,	-0.110598170118410e+10	},
		{  -4,	  6,	-0.323899915729960e+13	},
		{  -4,	  8,	 0.923814007023250e+16	},
		{  -3,	 -2,	 0.842250080413710e-12	},
		{  -3,	  5,	 0.663221436245510e+12	},
		{  -3,	  6,	-0.167170186672140e+15	},
		{  -2,	  2,	 0.253749358701390e+04	},
		{  -1,	 -6,	-0.819731559610520e-20	},
		{   0,	  3,	 0.328380587890660e+12	},
		{   1,	  1,	-0.625004791171540e+08	},
		{   2,	  6,	 0.803197957462020e+21	},
		{   3,	 -6,	-0.204397011338350e-10	},
		{   3,	 -2,	-0.378391047055940e+04	},
		{   6,	 -6,	 0.972876545938620e-02	},
		{   6,	 -5,	 0.154355721681460e+02	},
		{   6,	 -4,	-0.373962862928640e+04	},
		{   6,	 -1,	-0.682859011374570e+11	},
		{   8,	 -8,	-0.248488015614540e-03	},
		{   8,	 -4,	 0.394536049497070e+07	},
	};
	enum {
		R3A,
		R3B,
		R3C,
		R3D,
		R3E,
		R3F,
		R3G,
		R3H,
		R3I,
		R3J,
		R3K,
		R3L,
		R3M,
		R3N,
		R3O,
		R3P,
		R3Q,
		R3R,
		R3S,
		R3T,
		R3U,
		R3V,
		R3W,
		R3X,
		R3Y,
		R3Z
	} r3;

	int i, n, e;
	const if97_coef *coef;
	double vs, xi, xj;
	double ans = 0.0;

	assert(p < 100.0 && p > pi_r4(623.15));
	if (p > 40.0) {
		if (t <= t_b3(p, B3AB)) r3 = R3A;
		else r3 = R3B;
	} else if (p > 25.0) {
		if (t <= t_b3(p, B3CD)) r3 = R3C;
		else if (t <= t_b3(p, B3AB)) r3 = R3D;
		else if (t <= t_b3(p, B3EF)) r3 = R3E;
		else r3 = R3F;
	} else if (p > 23.5) {
		if (t <= t_b3(p, B3CD)) r3 = R3C;
		else if (t <= t_b3(p, B3GH)) r3 = R3G;
		else if (t <= t_b3(p, B3EF)) r3 = R3H;
		else if (t <= t_b3(p, B3IJ)) r3 = R3I;
		else if (t <= t_b3(p, B3JK)) r3 = R3J;
		else r3 = R3K;
	} else if (p > 23.0) {
		if (t <= t_b3(p, B3CD)) r3 = R3C;
		else if (t <= t_b3(p, B3GH)) r3 = R3L;
		else if (t <= t_b3(p, B3EF)) r3 = R3H;
		else if (t <= t_b3(p, B3IJ)) r3 = R3I;
		else if (t <= t_b3(p, B3JK)) r3 = R3J;
		else r3 = R3K;
	} else if (p > 22.5) {
		if (t <= t_b3(p, B3CD)) r3 = R3C;
		else if (t <= t_b3(p, B3GH)) r3 = R3L;
		else if (t <= t_b3(p, B3MN)) r3 = R3M;
		else if (t <= t_b3(p, B3EF)) r3 = R3N;
		else if (t <= t_b3(p, B3OP)) r3 = R3O;
		else if (t <= t_b3(p, B3IJ)) r3 = R3P;
		else if (t <= t_b3(p, B3JK)) r3 = R3J;
		else r3 = R3K;
	} else if (p > 22.11) {
		if (t <= t_b3(p, B3CD)) r3 = R3C;
		else if (t <= t_b3(p, B3QU)) r3 = R3Q;
		else if (t <= t_b3(p, B3UV)) r3 = R3U;
		else if (t <= t_b3(p, B3EF)) r3 = R3V;
		else if (t <= t_b3(p, B3WX)) r3 = R3W;
		else if (t <= t_b3(p, B3RX)) r3 = R3X;
		else if (t <= t_b3(p, B3JK)) r3 = R3R;
		else r3 = R3K;
	} else if (p > IAPWS_PC) {
		if (t <= t_b3(p, B3CD)) r3 = R3C;
		else if (t <= t_b3(p, B3QU)) r3 = R3Q;
		else if (t <= t_b3(p, B3UV)) r3 = R3U;
		else if (t <= t_b3(p, B3EF)) r3 = R3Y;
		else if (t <= t_b3(p, B3WX)) r3 = R3Z;
		else if (t <= t_b3(p, B3RX)) r3 = R3X;
		else if (t <= t_b3(p, B3JK)) r3 = R3R;
		else r3 = R3K;
	} else if (p > 21.93161551) {
		if (t <= t_b3(p, B3CD)) r3 = R3C;
		else if (t <= t_b3(p, B3QU)) r3 = R3Q;
		else if (t <= t_b3(p, B3UV)) r3 = R3U;
		else if (t <= theta_r4(p)) r3 = R3Y;
		else if (t <= t_b3(p, B3WX)) r3 = R3Z;
		else if (t <= t_b3(p, B3RX)) r3 = R3X;
		else if (t <= t_b3(p, B3JK)) r3 = R3R;
		else r3 = R3K;
	} else if (p > 21.90096265) {
		if (t <= t_b3(p, B3CD)) r3 = R3C;
		else if (t <= t_b3(p, B3QU)) r3 = R3Q;
		else if (t <= theta_r4(p)) r3 = R3U;
		else if (t <= t_b3(p, B3WX)) r3 = R3Z;
		else if (t <= t_b3(p, B3RX)) r3 = R3X;
		else if (t <= t_b3(p, B3JK)) r3 = R3R;
		else r3 = R3K;
	} else if (p > pi_r4(643.15)) {
		if (t <= t_b3(p, B3CD)) r3 = R3C;
		else if (t <= t_b3(p, B3QU)) r3 = R3Q;
		else if (t <= theta_r4(p)) r3 = R3U;
		else if (t <= t_b3(p, B3RX)) r3 = R3X;
		else if (t <= t_b3(p, B3JK)) r3 = R3R;
		else r3 = R3K;
	} else if (p > 20.5) {
		if (t <= t_b3(p, B3CD)) r3 = R3C;
		else if (t <= theta_r4(p)) r3 = R3S;
		else if (t <= t_b3(p, B3JK)) r3 = R3R;
		else r3 = R3K;
	} else if (p > 19.00881189173929) {
		if (t <= t_b3(p, B3CD)) r3 = R3C;
		else if (t <= theta_r4(p)) r3 = R3S;
		else r3 = R3T;
	} else if (p > pi_r4(623.15)) {
		if (t <= theta_r4(p)) r3 = R3C;
		else r3 = R3T;
	} else {
		assert(0);
	}

	switch (r3) {
		case R3A:
			vs = 2.4e-3;
			p = p / 100.0 - 0.085;
			t = t / 760.0 - 0.817;
			e = 1;
			n = SIZEA;
			coef = coefa;
			break;
		case R3B:
			vs = 4.1e-3;
			p = p / 100.0 - 0.280;
			t = t / 860.0 - 0.779;
			e = 1;
			n = SIZEB;
			coef = coefb;
			break;
		case R3C:
			vs = 2.2e-3;
			p = p / 40.0 - 0.259;
			t = t / 690.0 - 0.903;
			e = 1;
			n = SIZEC;
			coef = coefc;
			break;
		case R3D:
			vs = 2.9e-3;
			p = p / 40.0 - 0.559;
			t = t / 690.0 - 0.939;
			e = 4;
			n = SIZED;
			coef = coefd;
			break;
		case R3E:
			vs = 3.2e-3;
			p = p / 40.0 - 0.587;
			t = t / 710.0 - 0.918;
			e = 1;
			n = SIZEE;
			coef = coefe;
			break;
		case R3F:
			vs = 6.4e-3;
			p = sqrt(p / 40.0 - 0.587);
			t = t / 730.0 - 0.891;
			e = 4;
			n = SIZEF;
			coef = coeff;
			break;
		case R3G:
			vs = 2.7e-3;
			p = p / 25.0 - 0.872;
			t = t / 660.0 - 0.971;
			e = 4;
			n = SIZEG;
			coef = coefg;
			break;
		case R3H:
			vs = 3.2e-3;
			p = p / 25.0 - 0.898;
			t = t / 660.0 - 0.983;
			e = 4;
			n = SIZEH;
			coef = coefh;
			break;
		case R3I:
			vs = 4.1e-3;
			p = sqrt(p / 25.0 - 0.910);
			t = t / 660.0 - 0.984;
			e = 4;
			n = SIZEI;
			coef = coefi;
			break;
		case R3J:
			vs = 5.4e-3;
			p = sqrt(p / 25.0 - 0.875);
			t = t / 670.0 - 0.964;
			e = 4;
			n = SIZEJ;
			coef = coefj;
			break;
		case R3K:
			vs = 7.7e-3;
			p = p / 25.0 - 0.802;
			t = t / 680.0 - 0.935;
			e = 1;
			n = SIZEK;
			coef = coefk;
			break;
		case R3L:
			vs = 2.6e-3;
			p = p / 24.0 - 0.908;
			t = t / 650.0 - 0.989;
			e = 4;
			n = SIZEL;
			coef = coefl;
			break;
		case R3M:
			vs = 2.8e-3;
			p = p / 23.0 - 1.000;
			t = POW(t / 650.0 - 0.997, 0.25);
			e = 1;
			n = SIZEM;
			coef = coefm;
			break;
		case R3N:
			vs = 3.1e-3;
			p = p / 23.0 - 0.976;
			t = t / 650.0 - 0.997;
			e = 1;
			n = SIZEN;
			coef = coefn;
			break;
		case R3O:
			vs = 3.4e-3;
			p = sqrt(p / 23.0 - 0.974);
			t = t / 650.0 - 0.996;
			e = 1;
			n = SIZEO;
			coef = coefo;
			break;
		case R3P:
			vs = 4.1e-3;
			p = sqrt(p / 23.0 - 0.972);
			t = t / 650.0 - 0.997;
			e = 1;
			n = SIZEP;
			coef = coefp;
			break;
		case R3Q:
			vs = 2.2e-3;
			p = p / 23.0 - 0.848;
			t = t / 650.0 - 0.983;
			e = 4;
			n = SIZEQ;
			coef = coefq;
			break;
		case R3R:
			vs = 5.4e-3;
			p = p / 23.0 - 0.874;
			t = t / 650.0 - 0.982;
			e = 1;
			n = SIZER;
			coef = coefr;
			break;
		case R3S:
			vs = 2.2e-3;
			p = p / 21.0 - 0.886;
			t = t / 640.0 - 0.990;
			e = 4;
			n = SIZES;
			coef = coefs;
			break;
		case R3T:
			vs = 8.8e-3;
			p = p / 20.0 - 0.803;
			t = t / 650.0 - 1.020;
			e = 1;
			n = SIZET;
			coef = coeft;
			break;
		case R3U:
			vs = 2.6e-3;
			p = p / 23.0 - 0.902;
			t = t / 650.0 - 0.988;
			e = 1;
			n = SIZEU;
			coef = coefu;
			break;
		case R3V:
			vs = 3.1e-3;
			p = p / 23.0 - 0.960;
			t = t / 650.0 - 0.995;
			e = 1;
			n = SIZEV;
			coef = coefv;
			break;
		case R3W:
			vs = 3.9e-3;
			p = p / 23.0 - 0.959;
			t = t / 650.0 - 0.995;
			e = 4;
			n = SIZEW;
			coef = coefw;
			break;
		case R3X:
			vs = 4.9e-3;
			p = p / 23.0 - 0.910;
			t = t / 650.0 - 0.988;
			e = 1;
			n = SIZEX;
			coef = coefx;
			break;
		case R3Y:
			vs = 3.1e-3;
			p = p / 22.0 - 0.996;
			t = t / 650.0 - 0.994;
			e = 4;
			n = SIZEY;
			coef = coefy;
			break;
		case R3Z:
			vs = 3.8e-3;
			p = p / 22.0 - 0.993;
			t = t / 650.0 - 0.994;
			e = 4;
			n = SIZEZ;
			coef = coefz;
			break;
		default:
			assert(0);
	}

	for (i = 0; i < n; ++i) {
		if (i != 0) {
			xi *= powint(p, coef[i].I - coef[i - 1].I);
			//xj *= powint(t, coef[i].J - coef[i - 1].J);
			xj = powint(t, coef[i].J);
		} else {
			xi = powint(p, coef[i].I);
			xj = powint(t, coef[i].J);
		}
		ans += coef[i].n * xi * xj;
	}

	return (r3 != R3N ? powint(ans, e) : exp(ans)) * vs;
}

#if 0
#include <stdio.h>
int main()
{
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
