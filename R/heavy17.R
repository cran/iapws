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

.heavy17 <- function(what, rho, t)
{
	w <- .check_what(what)
	x <- callWrapper(C_heavy17, w = w, r = rho, t = t,
			 what = c("integer", "double", "double"))
	colnames(x) <- what
	x
}

heavy17 <- function(what, p, t, rho, state = NULL)
{
	w <- .check_what(what)
	x <- if (!missing(p) && !missing(t)) {
		if (is.null(state)) {
			state <- heavy17_state(p = p, t = t)
		}
		s <- .check_state(state)
		callWrapper(C_heavy17_pt, w = w, p = p, t = t, s = s,
			    what = c("integer", "double", "double", "integer"))
	} else if (!missing(rho) && !missing(t)) {
		if (is.null(state)) {
			state <- heavy17_state(rho = rho, t = t)
		}
		s <- .check_state(state)
		callWrapper(C_heavy17_rhot, w = w, r = rho, t = t, s = s,
			    what = c("integer", "double", "double"))
	} else {
		stop("invalid combination of arguments")
	}
	colnames(x) <- what
	x
}

heavy17_sat <- function(what, p, t)
{
	w <- .check_what(what)
	x <- if (!missing(t)) {
		callWrapper(C_heavy17_sat_t, w = w, t = t,
			    what = c("integer", "double"))
	} else if (!missing(p)) {
		callWrapper(C_heavy17_sat_p, w = w, p = p,
			    what = c("integer", "double"))
	} else {
		stop("invalid combination of arguments")
	}
	dimnames(x) <- list(NULL, what, c("liquid", "gas"))
	x
}

heavy17_state <- function(p, t, rho)
{
	s <- if (!missing(p) && !missing(t)) {
		callWrapper(C_heavy17_state_pt, p = p, t = t)
	} else if (!missing(rho) && !missing(t)) {
		callWrapper(C_heavy17_state_rhot, rho = rho, t = t)
	} else {
		stop("invalid combination of arguments")
	}
	names(.IAPWS_STATES)[match(s, .IAPWS_STATES)]
}

