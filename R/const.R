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

IAPWS <- c(R    = 8.3144598,	# J/K/mol
	   M    = 18.015268,	# g/mol
	   Pc   = 22.064,	# MPa
	   Tc   = 647.096,	# K
	   RHOc = 322.0,	# kg/m3
	   Pt   = 611.657e-6,	# MPa
	   Tt   = 273.16)	# K

IAPWS95 <- IAPWS
IAPWS95["Pt"] <- 611.654771e-06

HEAVY17 <- c(R    = unname(IAPWS["R"]),
	     M    = 20.027508,	# g/mol
	     Pc   = 21.661831,	# MPa
	     Tc   = 643.847,	# K
	     RHOc = 356.0,	# kg/m3
	     Pt   = 661.59e-6,	# MPa
	     Tt   = 276.969)	# K
