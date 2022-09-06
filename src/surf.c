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
 * IAPWS R1-76(2014), Revised Release on Surface Tension of Ordinary Water
 * Substance (2014)
 */

#include "iapws.h"

static double surf(double t)
{
	t  = 1.0 - t / IAPWS_TC;
	return POW(t, 1.256) * (1.0 - t * 0.625) * 235.8;
}

double iapws_sigma(const iapws_phi *phi)  /* mN/m */
{
	return surf(iapws_t(phi));
}
