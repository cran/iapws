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

.nroot_control <- function(trace = 0L, maxit = 100L,
			   abstol = 1e-9, reltol = 1e-9)
{
	invisible(callWrapper(C_nroot_control, trace, maxit, abstol, reltol,
			      what = c("integer", "integer",
				       "double", "double")))
}

