\name{heavy17}
\alias{heavy17}
\title{
  IAPWS-17 Formulation in Single Phase Regions
}
\description{
  Compute thermodynamic and transport properties of heavy water and
  steam using the IAPWS formulation 2017.
}
\usage{
heavy17(what, p, t, rho, state = NULL)
}
\arguments{
  \item{what}{
    a character vector listing the output properties.
    See Details for available properties.
  }
  \item{p}{
    a numeric vector giving the pressure values in MPa.
  }
  \item{t}{
    a numeric vector giving the temperature values in K.
  }
  \item{rho}{
    a numeric vector giving the density values in kg/m\eqn{^3}.
  }
  \item{state}{
    \code{NULL} or a character vector giving the physical state.
    One of \code{"liquid"}, \code{"gas"} or \code{"supercritical"}.
    In case of \code{NULL}, the state is determined by the function
    \code{heavy17_state()}.
  }
}
\details{
  The available properties for \code{what} are:
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
    \item \code{"kappat"}: the isothermal compressibility in 1/MPa.
    \item \code{"w"}: the speed of sound in m/s.
    \item \code{"eta"}: the dynamic viscosity in \eqn{\mu}Pa.s.
    \item \code{"lambda"}: the thermal conductivity in mW/K/m.
    %\item \code{"sigma"}: the surface tension in mN/m.
    %\item \code{"epsilon"}: the dielectric constant.
  }

  A valid pair of \code{p}, \code{rho} and \code{t} is needed:
  either \code{(p, t)} or \code{(rho, t)}.
}
\value{
  A numeric matrix of dimension \code{c(n, length(what))} with \code{n}
  the maximum length of either \code{rho}, \code{p}, \code{t}
  or \code{state}.
}
\references{
  International Association for the Properties of Water and Steam,
  IAPWS R16-17(2018), \emph{Revised Release on the IAPWS Formulation 2017
  for the Thermodynamic Properties of Heavy Water.}

  International Association for the Properties of Water and Steam,
  IAPWS R17-20, \emph{Release on the IAPWS Formulation 2020
  for the Viscosity of Heavy Water.}

  International Association for the Properties of Water and Steam,
  IAPWS R18-21, \emph{Release on the IAPWS Formulation 2021
  for the Thermal Conductivity of Heavy Water.}
}
\author{
  Jonathan Debove
}
\note{
  Computing several properties in a single call may be more efficient
  than separately because most of the computation time is shared.
}
\examples{
# IAPWS-17 Tab. 7
M <- 20.027508  # g/mol
tab <- heavy17(c("p", "cv", "w", "s"),
               t = c(300, 300, 300,
                     500, 500, 500, 500, 500,
                     643.8,
                     800, 800),
               rho = c(0.55126e2, 0.6e2, 0.65e2,
                       0.5e-1, 0.5, 0.4626e2, 0.5e2, 0.6e2,
                       0.2e2,
                       0.1e-1, 0.25) * M)
tab[, c("cv", "s")] <- tab[, c("cv", "s")] / M
format(tab, scientific = TRUE, digits = 9)
\dontshow{
stopifnot(all.equal(tab[, "p", drop = FALSE],
          matrix(c(0.529123711e-1, 0.238222326e3, 0.626176781e3, 0.206052588,
		   0.188967446e1, 0.835329492e1, 0.107462884e3, 0.721798322e3,
		   0.216503820e2, 0.664864175e-1, 0.164466177e1),
                 11, 1, dimnames = list(NULL, c("p"))), tolerance = 1e-9))
}

# Viscosity Tab. 3
tab <- heavy17("eta", t = c(298.15, 298.15, 298.15, 373.15,
                            775.0, 775.0, 775.0),
               rho = c(0, 1105, 1130, 1064, 1, 100, 400))
format(tab, digits = 8)

# Viscosity Tab. 4
tab <- heavy17("eta", t = 644.101, rho = c(145, 245, 295, 345, 395, 445))
format(tab, digits = 8)
\dontshow{
stopifnot(all.equal(tab[, "eta", drop = FALSE],
          matrix(c(26.640959, 32.119967, 36.828275,
		   43.225017, 47.193530, 50.241640),
                 6, 1, dimnames = list(NULL, c("eta"))), tolerance = 1e-7))
}

# Thermal conductivity Tab. 3
tab <- heavy17("lambda", t = c(298.15, 298.15, 298.15, 825),
	       rho = c(0, 1104.5, 1200, 0))
format(tab, digits = 6)
\dontshow{
stopifnot(all.equal(tab[, "lambda", drop = FALSE],
          matrix(c(17.7498, 599.557, 690.421, 76.4492),
                 4, 1, dimnames = list(NULL, c("lambda"))), tolerance = 1e-6))
}

# Thermal conductivity Tab. 4
tab <- heavy17("lambda", t = 644.10,
	       rho = c(1, 106, 256, 306, 356, 406, 456, 750))
format(tab, digits = 6)
\dontshow{
stopifnot(all.equal(tab[, "lambda", drop = FALSE],
          matrix(c(52.4527, 103.342, 394.612, 801.382,
		   1278.423, 670.833, 423.603, 454.846),
                 8, 1, dimnames = list(NULL, c("lambda"))), tolerance = 1e-6))
}
}
