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

#ifndef IAPWS_HEAVY17_H
#define IAPWS_HEAVY17_H

#include "iapws.h"

#define HEAVY17_M	20.027508		/* g/mol   */
#define HEAVY17_R	(8.3144598 / HEAVY17_M)	/* kJ/kg/K */

#define HEAVY17_PC	21.661831		/* MPa     */
#define HEAVY17_TC	643.847			/* K       */
#define HEAVY17_RHOC	(17.77555 * HEAVY17_M)	/* kg/m3   */

#define HEAVY17_PT	661.59e-6		/* MPa     */
#define HEAVY17_TT	276.969			/* K       */

void heavy17_phi(double rho, double t, struct iapws_phi *phi);
int heavy17_phi_rhot(double rho, double t, enum iapws_state state,
		struct iapws_phi *phi);
int heavy17_phi_pt(double p, double t, enum iapws_state state,
		struct iapws_phi *phi);
int heavy17_sat_t(double t, struct iapws_phi *phil, struct iapws_phi *phig);
int heavy17_sat_p(double p, struct iapws_phi *phil, struct iapws_phi *phig);
enum iapws_state heavy17_state_pt(double p, double t);
enum iapws_state heavy17_state_rhot(double rho, double t);

double heavy17_psat(double t);
double heavy17_tsat(double p);
double heavy17_rhol(double t);
double heavy17_rhog(double t);

#endif
