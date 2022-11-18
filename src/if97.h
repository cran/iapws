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
#include "nroot.h"

#define IF97_R		0.461526	/* kJ/kg/K */
#define IF97_PT		611.657e-06	/* MPa */

enum if97_region {
	IF97_REGION_UNDEF = 0,
	IF97_REGION_1 = 1,
	IF97_REGION_2 = 2,
	IF97_REGION_3 = 3,
	IF97_REGION_4 = 4,
	IF97_REGION_5 = 5,
};

enum iapws_state if97_state_pt(double p, double t);
enum if97_region if97_region_pt(double p, double t);
enum if97_region if97_region_ph(double p, double h);
double if97_tsat(double p);
double if97_psat(double t);
int if97_gamma_pt(double p, double t, enum iapws_state state,
		struct iapws_phi *gamma);
int if97_gamma_ph(double p, double h, struct iapws_phi *gamma);

#endif
