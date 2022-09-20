\name{iapws95_state}
\alias{iapws95_state}
\alias{iapws95_state_rhot}
\title{
  IAPWS-95 Physical States
}
\description{
  Compute the physical state of water according the IAPWS formulation 1995.
}
\usage{
iapws95_state(p, t)
iapws95_state_rhot(rho, t)
}
\arguments{
  \item{rho}{a numeric vector giving the values of density in kg/m^3.}
  \item{p}{a numeric vector giving the pressure values in MPa.}
  \item{t}{a numeric vector giving the temperature values in K.}
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
iapws95_state(c(.1, .1, 23), c(293.15, 373.15, 650))
iapws95_state_rhot(c(1000, .1, 500, 500), c(293.15, 373.15, 650, 373.15))
}
