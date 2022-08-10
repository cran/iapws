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

#ifndef IAPWS_IF97_H
#define IAPWS_IF97_H

#include "iapws.h"

#define IF97_R		0.461526	/* kJ/kg/K */

typedef enum {
	IF97_UNDEF = 0,
	IF97_WATER = 1,
	IF97_STEAM = 2,
	IF97_SUPER = 3,
	IF97_SAT   = 4,
	IF97_GAS   = 5,
} if97_region_id;

iapws_state_id if97_state(double p, double t);
if97_region_id if97_region(double p, double t);
double if97_tsat(double p);
double if97_psat(double t);
int if97_gamma(double p, double t, iapws_state_id state, iapws_phi *gamma);
//void if97_gamma_rhot(double rho, double t, iapws_phi *gamma);

#endif
