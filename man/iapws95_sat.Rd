\name{iapws95_sat}
\alias{iapws95_sat}
\title{
  IAPWS-95 Formulation in the Saturation Region
}
\description{
  Compute thermodynamic and transport properties of water and steam
  along the saturated line according to the IAPWS formulation 1995.
}
\usage{
iapws95_sat(what, p, t)
}
\arguments{
  \item{what}{
    a character vector listing the output properties.
    See Details for available properties.
  }
  \item{t}{
    a numeric vector giving the temperature values in K.
  }
  \item{p}{
    a numeric vector giving the pressure values in MPa.
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

  Only one of \code{t} or \code{p} is needed.
  If both arguments are given, \code{p} is ignored.
}
\value{
  A numeric array of dimension \code{c(n, length(what), 2L)} with
  \code{n} the length of either \code{p} or \code{t}. The last dimension
  indicate the physical state (\code{"liquid"} or \code{"gas"}).
}
\references{
  International Association for the Properties of Water and Steam,
  IAPWS R6-95(2018), \emph{Revised Release on the IAPWS Formulation 1995
  for the Thermodynamic Properties of Ordinary Water Substance
  for General and Scientific Use.}
}
\author{
  Jonathan Debove
}
\note{
  Computing several properties in a single call may be more efficient
  than separately because most of the computation time is shared.
}
\examples{
# IAPWS-95 Tab. 8
tab <- iapws95_sat(c("p", "rho", "h", "s"), t = c(275, 450, 625))
format(tab, scientific = TRUE, digits = 9)
\dontshow{
stopifnot(all.equal(tab[, "p", "gas", drop = FALSE],
          array(c(0.698451167e-3, 0.932203564e0, 0.169082693e2),
                c(3, 1, 1), dimnames = list(NULL, c("p"), c("gas"))),
          tolerance = 2e-9))
}
}
