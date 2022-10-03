\name{iapws95_state}
\alias{iapws95_state}
\title{
  IAPWS-95 Physical States
}
\description{
  Compute the physical state of water according the IAPWS formulation 1995.
}
\usage{
iapws95_state(p, t, rho)
}
\arguments{
  \item{p}{a numeric vector giving the pressure values in MPa.}
  \item{t}{a numeric vector giving the temperature values in K.}
  \item{rho}{a numeric vector giving the values of density in kg/m^3.}
}
\details{
  A valid pair of \code{p}, \code{rho}, \code{t} is needed:
  either \code{(p, t)} or \code{(rho, t)}.
}
\value{
  A character vector giving the physical state.
  One of \code{"solid"}, \code{"liquid"}, \code{"gas"},
  \code{"supercritical"}, \code{"saturated"} or \code{"undef"}.
}
\author{
  Jonathan Debove
}
\examples{
iapws95_state(p = c(.1, .1, 23), t = c(293.15, 373.15, 650))
iapws95_state(rho = c(1000, .1, 500, 500),
              t = c(293.15, 373.15, 650, 373.15))
}
