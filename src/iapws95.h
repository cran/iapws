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

#ifndef IAPWS_IAPWS95_H
#define IAPWS_IAPWS95_H

#include "iapws.h"

#define IAPWS95_R	0.46151805	/* kJ/kg/K */
#define IAPWS95_PT	611.654771e-06	/* MPa */

void iapws95_phi(double rho, double t, struct iapws_phi *phi);
int iapws95_phi_rhot(double rho, double t, enum iapws_state state, struct iapws_phi *phi);
int iapws95_phi_pt(double p, double t, enum iapws_state state, struct iapws_phi *phi);
int iapws95_phi_ph(double p, double h, struct iapws_phi *phi);
int iapws95_sat_t(double t, struct iapws_phi *phil, struct iapws_phi *phig);
int iapws95_sat_p(double p, struct iapws_phi *phil, struct iapws_phi *phig);
enum iapws_state iapws95_state_pt(double p, double t);
enum iapws_state iapws95_state_rhot(double rho, double t);

#endif
