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

if97 <- function(what, p, t, h, state = NULL)
{
	w <- .check_what(what)
	x <- if (!missing(p) && !missing(t)) {
		if (is.null(state)) {
			state <- if97_state(p, t)
		}
		s <- .check_state(state)
		callWrapper(C_if97_pt, w = w, p = p, t = t, s = s,
			    what = c("integer", "double", "double", "integer"))
	} else if (!missing(p) && !missing(h)) {
		callWrapper(C_if97_ph, w = w, p = p, h = h,
			    what = c("integer", "double", "double"))
	} else {
		stop("invalid combination of arguments")
	}
	colnames(x) <- what
	x
}

if97_psat <- function(t)
{
	x <- callWrapper(C_if97_psat, t = t)
	is.na(x) <- x == 0
	x
}

if97_tsat <- function(p)
{
	x <- callWrapper(C_if97_tsat, p = p)
	is.na(x) <- x == 0
	x
}

if97_state <- function(p, t)
{
	s <- callWrapper(C_if97_state_pt, p = p, t = t)
	names(.IAPWS_STATES)[match(s, .IAPWS_STATES)]
}

