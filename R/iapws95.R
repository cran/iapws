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

iapws95 <- function(what, rho, t)
{
	w <- .check_what(what)
	x <- callWrapper(R_iapws95, w = w, r = rho, t = t,
			 what = c("integer", "double", "double"))
	colnames(x) <- what
	x
}

iapws95_pt <- function(what, p, t, state = iapws95_state(p, t))
{
	w <- .check_what(what)
	s <- .check_state(state)
	x <- callWrapper(R_iapws95_pt, w = w, p = p, t = t, s = s,
			 what = c("integer", "double", "double", "integer"))
	colnames(x) <- what
	x
}

iapws95_sat <- function(what, t)
{
	w <- .check_what(what)
	x <- callWrapper(R_iapws95_sat, w = w, t = t,
			 what = c("integer", "double"))
	dimnames(x) <- list(NULL, what, c("liquid", "gas"))
	x
}

iapws95_sat_p <- function(what, p)
{
	w <- .check_what(what)
	x <- callWrapper(R_iapws95_sat_p, w = w, p = p,
			 what = c("integer", "double"))
	dimnames(x) <- list(NULL, what, c("liquid", "gas"))
	x
}

iapws95_state <- function(p, t) {
	s <- callWrapper(R_iapws95_state, p = p, t = t)
	names(.IAPWS_STATE)[match(s, .IAPWS_STATE)]
}

iapws95_state_rhot <- function(rho, t) {
	s <- callWrapper(R_iapws95_state_rhot, rho = rho, t = t)
	names(.IAPWS_STATE)[match(s, .IAPWS_STATE)]
}
