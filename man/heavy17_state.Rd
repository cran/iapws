\name{heavy17_state}
\alias{heavy17_state}
\title{
  IAPWS-17 Physical States
}
\description{
  Identify the physical state of heavy water according to the IAPWS formulation
  2017.
}
\usage{
heavy17_state(p, t, rho)
}
\arguments{
  \item{p}{a numeric vector giving the pressure values in MPa.}
  \item{t}{a numeric vector giving the temperature values in K.}
  \item{rho}{a numeric vector giving the density values in kg/m\eqn{^3}.}
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
heavy17_state(p = c(.1, .1, 23), t = c(293.15, 393.15, 650))
heavy17_state(rho = c(1200, .1, 500, 500),
              t = c(293.15, 373.15, 650, 373.15))
}
