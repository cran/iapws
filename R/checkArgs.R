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

.IAPWS_STATES <- c(undef  = -1L,
		  solid  = 0L,
		  liquid = 1L,
		  gas    = 2L,
		  supercritical = 3L,
		  saturated = 4L)

.check_state <- function(state)
{
	.call <- format(sys.call(sys.parent()))
	if (!is.character(state)) {
		stop(.call, ": 'state' is not a character vector", call. = FALSE)
	} else if (!all(state %in% names(.IAPWS_STATES))) {
		stop(.call, ": invalid 'state'", call. = FALSE)
	}
	.IAPWS_STATES[state]
}

.IAPWS_PROPERTIES <- c(f  =  0L, g  =  1L,
		 u  =  2L, h  =  3L,
		 s  =  4L, t  =  5L,
		 p  =  6L, v  =  7L,
		 cp =  8L, cv =  9L,
		 w  = 10L, rho = 11L,
		 alpha = 12L, beta = 13L,
		 chit = 14L,
		 eta = 15L, lambda = 16L,
		 sigma = 17L, epsilon = 18L)

.check_what <- function(what)
{
	.call <- format(sys.call(sys.parent()))
	if (!is.character(what)) {
		stop(.call, ": 'what' is not a character vector", call. = FALSE)
	} else if (!all(what %in% names(.IAPWS_PROPERTIES))) {
		stop(.call, ": invalid 'what'", call. = FALSE)
	}
	.IAPWS_PROPERTIES[what]
}

