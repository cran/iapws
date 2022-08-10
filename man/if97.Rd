\name{if97}
\alias{if97}
\title{
	IAPWS-IF97 Basic Equations for Regions 1, 2, 3 and 5
}
\description{
	Compute thermodynamic and transport properties of water and steam
	using the IAPWS industrial formulation 1997.
}
\usage{
if97(what, p, t, state = if97_state(p, t))
}
\arguments{
  \item{what}{
	  a character vector listing the output properties. See Details for
	  available properties.
}
  \item{p}{
	  a numeric vector giving the pressure values in MPa.
}
  \item{t}{
	  a numeric vector giving the temperature values in K.
}
  \item{state}{
	  a character vector giving the physical state. One of
	  \code{"liquid"}, \code{"gas"} or \code{"supercritical"}.
}
}
\details{
	The currently available properties are:
	\itemize{
		\item \code{"f"}: the specific free energy in kJ/kg.
		\item \code{"g"}: the specific Gibbs enthalpy in kJ/kg.
		\item \code{"u"}: the specific internal energy in kJ/kg.
		\item \code{"h"}: the specific enthalpy in kJ/kg.
		\item \code{"s"}: the specific entropy in kJ/K/kg.
		\item \code{"t"}: the temperature in K.
		\item \code{"p"}: the pressure in MPa.
		\item \code{"v"}: the specific volume in m\eqn{^3}/kg.
		\item \code{"rho"}: the mass density in kg/m\eqn{^3}.
		\item \code{"cv"}: the specific isochoric heat capacity in kJ/K/kg.
		\item \code{"cp"}: the specific isobaric heat capacity in kJ/K/kg.
		\item \code{"alpha"}: the expansion coefficient in 1/K
		\item \code{"beta"}: the pressure coefficient in MPa/K.
		\item \code{"chit"}: the isothermal compressibility in 1/MPa.
		\item \code{"w"}: the speed of sound in m/s.
		\item \code{"eta"}: the dynamic viscosity in \eqn{10^{-6}}Pa.s.
		\item \code{"lambda"}: the thermal conductivity in mW/K/m.
	}
}
\value{
	A numeric matrix of dimension \code{c(n, length(what))} with \code{n}
	the maximum length of either \code{p}, \code{t} or \code{state}.
}
\references{
	International Association for the Properties of Water and Steam,
	IAPWS R7-97(2012), Revised Release on the IAPWS Industrial Formulation
	1997 for the Thermodynamic Properties of Water and Steam.

	International Association for the Properties of Water and Steam,
	IAPWS R12-08(2008), Release on the IAPWS Formulation 2008 for the
	Viscosity of Ordinary Water Substance.

	International Association for the Properties of Water and Steam,
	IAPWS R15-11(2011), Release on the IAPWS Formulation 2011 for the
	Thermal Conductivity of Ordinary Water Substance.
}
\author{
	Jonathan Debove
}
\note{
	Computing the needed properties simultaneously (in a single call)
	shoud be more efficient than separately because most of the computation
	is shared between the different properties.
}
\examples{
what <- c("v", "h", "u", "s", "cp", "w")

# Region 1
format(if97(what, t = c(300, 300, 500), p = c(3, 80, 3)),
       scientific = TRUE, digits = 9)

# Region 2
format(if97(what, t = c(300, 700, 700), p = c(3.5e-3, 3.5e-3, 30)),
       scientific = TRUE, digits = 9)

# Region 2 metastable
format(if97(what, t = c(450, 440, 450), p = c(1, 1, 1.5), state = "gas"),
       scientific = TRUE, digits = 9)

# Region 3
format(if97(what, t = c(650, 650, 750),
	    p = c(0.255837018e2, 0.222930643e2, 0.783095639e2)),
       scientific = TRUE, digits = 9)

# Region 5
format(if97(what, t = c(1500, 1500, 2000), p = c(0.5, 30, 30)),
       scientific = TRUE, digits = 9)

# Viscosity and thermal conductivity
what <- c("lambda", "eta")
## Region 1
format(if97(what, p = c(20, 50), t = c(620, 620)),
       scientific = TRUE, digits = 9)

## Region 2
format(if97(what, p = c(0.3, 50), t = c(650, 800)),
       scientific = TRUE, digits = 9)

## Region 3
#format(if97_rhot(c("lambda", "rho", "cp", "cv", "eta"),
#	    rho = c(0.3, 50), t = c(222, 322)),
#       scientific = TRUE, digits = 9)
}
