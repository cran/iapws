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
 * IAPWS R1-76(2014), Revised Release on Surface Tension of Ordinary Water
 * Substance (2014)
 *
 * International Association for the Properties of Water and Steam,
 * IAPWS R5-85(1994), IAPWS Release on Surface Tension of Heavy Water
 * Substance (1994)
 */

#include "iapws.h"
#include "heavy17.h"
#include "pow.h"

static double surf(double theta, const double B, const double b, const double mu)
{
	theta = 1.0 - theta;
	return theta >= 0.0 ? POW(theta, mu) * (1.0 + theta * b) * B : 0.0;
}

double iapws_sigma(double t)  /* mN/m */
{
	return surf(t / IAPWS_TC, 235.8, -0.625, 1.256);
}

double heavy17_sigma(double t)  /* mN/m */
{
	return surf(t / HEAVY17_TC, 238.0, -0.639, 1.25);
}

