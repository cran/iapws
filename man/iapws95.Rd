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
		\item \code{"sigma"}: the surface tension in mN/m.
		\item \code{"epsilon"}: the dielectric constant.
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

	International Association for the Properties of Water and Steam,
	IAPWS R1-76(2014), Revised Release on Surface Tension of Ordinary Water
	Substance

	International Association for the Properties of Water and Steam,
	IAPWS R8-97, Release on the Static Dielectric Constant of Ordinary
	Water Substance for Temperatures from 238 K to 873 K and Pressures
	up to 1000 MPa
}
\author{
	Jonathan Debove
}
\note{
	Computing several properties in a single call may be more efficient
	than separately because most of the computation time is shared.
}
\examples{
# IAPWS95 Tab. 7
tab7 <- iapws95(c("p", "cv", "w", "s"),
                t = c(300, 300, 300,
                      500, 500, 500, 500,
                      647,
                      900, 900, 900),
                rho = c(0.9965560e3, 0.1005308e4, 0.1188202e4,
                        0.4350000e0, 0.4532000e1, 0.8380250e3, 0.1084564e4,
                        0.3580000e3,
                        0.2410000e0, 0.5261500e2, 0.8707690e3))
format(tab7, scientific = TRUE, digits = 9)
\dontshow{
stopifnot(all.equal(tab7[, "p", drop = FALSE],
	  matrix(c(0.992418352e-1, 0.200022515e2, 0.700004704e3,
		   0.999679423e-1, 0.999938125e0, 0.100003858e2,
		   0.700000405e3, 0.220384756e2,
		   0.100062559e0, 0.200000690e2, 0.700000006e3),
		 11, 1, dimnames = list(NULL, c("p"))), tolerance = 1e-9))
}

# Viscosity Tab. 4
tab4 <- iapws95("eta",
                t = c(298.15, 298.15, 373.15,
                      433.15, 433.15,
                      873.15, 873.15, 873.15,
                      1173.15, 1173.15, 1173.15),
                rho = c(998, 1200, 1000,
                        1, 1000,
                        1, 100, 600,
                        1, 100, 400))
format(tab4, digits = 8)
\dontshow{
stopifnot(all.equal(tab4[, "eta", drop = FALSE],
	  matrix(c(889.735100, 1437.649467, 307.883622,
		   14.538324, 217.685358,
		   32.619287, 35.802262, 77.430195,
		   44.217245, 47.640433, 64.154608),
		 11, 1, dimnames = list(NULL, c("eta"))), tolerance = 1e-8))
}

# Viscosity Tab. 5
tab5 <- iapws95("eta", t = 647.35,
               rho = c(122, 222, 272, 322, 372, 422))
format(tab5, digits = 8)
\dontshow{
stopifnot(all.equal(tab5[, "eta", drop = FALSE],
	  matrix(c(25.520677, 31.337589, 36.228143,
		   42.961579, 45.688204, 49.436256),
		 6, 1, dimnames = list(NULL, c("eta"))), tolerance = 1e-8))
}

# Thermal conductivity Tab. 4
tab4 <- iapws95("lambda",
                t = c(298.15, 298.15, 298.15, 873.15),
                rho = c(0, 998, 1200, 0))
format(tab4, digits = 9)
\dontshow{
stopifnot(all.equal(tab4[, "lambda", drop = FALSE],
	  matrix(c(18.4341883, 607.712868, 799.038144, 79.1034659),
		 4, 1, dimnames = list(NULL, c("lambda"))), tolerance = 1e-9))
}

# Thermal conductivity Tab. 5
tab5 <- iapws95("lambda", t = 647.35,
                rho = c(1, 122, 222, 272, 322, 372, 422, 750))
format(tab5, digits = 9)
\dontshow{
stopifnot(all.equal(tab5[, "lambda", drop = FALSE],
	  matrix(c(51.9298924, 130.922885, 367.787459, 757.959776,
		   1443.75556, 650.319402, 448.883487, 600.961346),
		 8, 1, dimnames = list(NULL, c("lambda"))), tolerance = 1e-9))
}

# Dielectric constant Tab. 4
tab4 <- iapws95_pt("epsilon",
                   p = c(0.1013125, 0.1013125, 10, 1000,
                         10, 100, 500, 10, 100, 500),
                   t = c(240, 300, 300, 300,
                         650, 650, 650, 870, 870, 870),
		   state = c(rep("liquid", 4),
			     "gas", rep("supercritical", 2),
			     "gas", rep("supercritical", 2)))
format(tab4, digits = 6)
\dontshow{
stopifnot(all.equal(tab4[, "epsilon", drop = FALSE],
	  matrix(c(104.34982, 77.74735, 78.11269, 103.69632, 1.26715,
		   17.71733, 26.62132, 1.12721, 4.98281, 15.09746),
		 10, 1, dimnames = list(NULL, c("epsilon"))), tolerance = 1e-7))
}
}
