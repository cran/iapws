# iapws - IAPWS formulations for the properties of water and steam
# Copyright (C) 2022 Jonathan Debove
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

iapws_epsilon <- function(rho, t)
	callWrapper(C_iapws_epsilon, r = rho, t = t)

iapws_n <- function(rho, t, lambda)
	callWrapper(C_iapws_n, r = rho, t = t, l = lambda)

iapws_pk <- function(rho, t)
	callWrapper(C_iapws_pk, r = rho, t = t)

iapws_sigma <- function(t)
	callWrapper(C_iapws_sigma, t = t)

heavy17_sigma <- function(t)
	callWrapper(C_heavy17_sigma, t = t)

