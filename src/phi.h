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

#ifndef IAPWS_PHI_H
#define IAPWS_PHI_H

struct iapws_phi_coef0 {
	double n;
	double gamma;
};

struct iapws_phi_coef1 {
	int c;
	int d;
	double t;
	double n;
};

struct iapws_phi_coef2 {
	int d;
	double t;
	double n;
	double alpha;
	double beta;
	double gamma;
	double eps;
};

void iapws_phi(const double delta, const double tau,
		const struct iapws_phi_coef0 *coef0, const int n0,
		const struct iapws_phi_coef1 *coef1, const int n1,
		const struct iapws_phi_coef2 *coef2, const int n2,
		struct iapws_phi *phi);

#endif
