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

#IF97_REGION <- c(WATER = 1L, STEAM = 2L, SUPER = 3L, SAT = 4L, GAZ = 5L)

if97_region <- function(p, t)
	callWrapper(R_if97_region, p = p, t = t)

if97_state <- function(p, t) {
	s <- callWrapper(R_if97_state, p = p, t = t)
	names(.IAPWS_STATE)[match(s, .IAPWS_STATE)]
}

if97_psat <- function(t)
	callWrapper(R_if97_psat, t = t)

if97_tsat <- function(p)
	callWrapper(R_if97_tsat, p = p)

if97 <- function(what, p, t, state = if97_state(p, t))
{
	w <- .check_what(what)
	s <- .check_state(state)
	x <- callWrapper(R_if97, w = w, p = p, t = t, s = s,
			 what = c("integer", "double", "double", "integer"))
	colnames(x) <- what
	x
}

