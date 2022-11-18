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

#ifndef IAPWS_POW_H
#define IAPWS_POW_H

#define ARRAY_SIZE(x) (sizeof(x) / sizeof(x[0]))

#define MATHLIB_STANDALONE
#include <Rmath.h>
#define POW2(x)		((x) * (x))
#define POWINT(x, y)	powint((x), (y))
#define POW(x, y)	((y) == ((int)(y)) ? powint(x, (int)(y)) : R_pow(x, y))

inline double powint(double x, int i)
{
	double y;
	switch (i) {
		case 0: return 1.0;
		case 1: return x;
		case 2: return x*x;
		case 3: return x*x*x;
		case 4: x*=x; return x*x;
		case 5: y=x*x; return y*y*x;
		case 6: x*=x; return x*x*x;
		case 7: y=x*x; return y*y*y*x;
		case 8: x*=x; x*=x; return x*x;
		case 9: x*=x*x; return x*x*x;
		default: return R_pow_di(x, i);
	}
}

#endif
