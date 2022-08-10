\name{iapws95}
\alias{iapws95}
\alias{iapws95_pt}
\title{
	IAPWS Formulation 1995 in the Single Phase Regions
}
\description{
	Compute thermodynamic and transport properties of water and steam
	using the IAPWS formulation 1995.
}
\usage{
iapws95(what, rho, t)
iapws95_pt(what, p, t, state = iapws95_state(p, t))
}
\arguments{
  \item{what}{
	  a character vector listing the output properties. See Details for
	  available properties.
}
  \item{rho}{
	  a numeric vector giving the values of density in kg/m^3.
}
  \item{t}{
	  a numeric vector giving the temperature values in K.
}
  \item{p}{
	  a numeric vector giving the pressure values in MPa.
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
	the maximum length of either \code{rho}, \code{p}, \code{t} or
	\code{state}.
}
\references{
	International Association for the Properties of Water and Steam,
	IAPWS R6-95(2018), Revised Release on the IAPWS Formulation 1995 for
	the Thermodynamic Properties of Ordinary Water Substance for General
	and Scientific Use.

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
# IAPWS95
format(iapws95(c("p", "cv", "w", "s"),
	       t = c(300, 300, 300,
		     500, 500, 500, 500,
		     647,
		     900, 900, 900),
	       rho = c(0.9965560e3, 0.1005308e4, 0.1188202e4,
		       0.4350000e0, 0.4532000e1, 0.8380250e3, 0.1084564e4,
		       0.3580000e3,
		       0.2410000e0, 0.5261500e2, 0.8707690e3)),
       scientific = TRUE, digits = 9)

# Viscosity
format(iapws95("eta",
	       t = c(298.15, 298.15, 373.15,
		     433.15, 433.15,
		     873.15, 873.15, 873.15,
		     1173.15, 1173.15, 1173.15),
	       rho = c(998, 1200, 1000,
		       1, 1000,
		       1, 100, 600,
		       1, 100, 400)),
       digits = 8)
format(iapws95("eta", t = 647.35,
	       rho = c(122, 222, 272, 322, 372, 422)),
       digits = 8)

# Thermal conductivity
format(iapws95("lambda",
	       t = c(298.15, 298.15, 298.15, 873.15),
	       rho = c(0, 998, 1200, 0)),
       digits = 9)
format(iapws95("lambda", t = 647.35,
	       rho = c(1, 122, 222, 272, 322, 372, 422, 750)),
       digits = 9)
}
