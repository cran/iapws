\name{ice06}
\alias{ice06}
\title{
  Equation of State 2006 for Ice Ih
}
\description{
  Compute thermodynamic properties of ice Ih using the Equation of State 2006.
}
\usage{
ice06(what, p, t)
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
  }
}
\value{
  A numeric matrix of dimension \code{c(n, length(what))} with \code{n}
  the maximum length of either \code{p} or \code{t}.
}
\references{
  International Association for the Properties of Water and Steam,
  IAPWS R10-06(2009), \emph{Revised Release on the Equation of State 2006
  for H2O Ice Ih.}
}
\author{
  Jonathan Debove
}
\note{
  Computing several properties in a single call may be more efficient
  than separately because most of the computation time is shared.
}

\examples{
tab <- ice06(c("g", "h", "f", "u", "s", "cp", "rho", "alpha", "beta", "kappat"),
              p = c(611.657e-6, .101325, 100), t = c(273.16, 273.152519, 100))
format(tab, scientific = TRUE, digits = 11)
\dontshow{
stopifnot(all.equal(tab[, "f", drop = FALSE],
          matrix(c(-0.55446875e-4, -0.918701567e-2, -0.328489902347e3),
                 3, 1, dimnames = list(NULL, c("f"))), tolerance = 1e-8))
}
}
