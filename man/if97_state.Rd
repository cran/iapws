\name{if97_state}
\alias{if97_state}
\title{
	IAPWS-IF97 Physical State
}
\description{
	Compute the physical state of water according the IAPWS industrial
	formulation 1997.
}
\usage{
if97_state(p, t)
}
\arguments{
  \item{p}{
    a numeric vector giving the pressure values in MPa.
  }
  \item{t}{
    a numeric vector giving the temperature values in K.
  }
}
\value{
	A character vector giving the physical state. One of \code{"solid"},
	\code{"liquid"}, \code{"gas"}, \code{"supercritical"} or \code{"undef"}.
}
\author{
	Jonathan Debove
}
\examples{
if97_state(.1, 293.15)
if97_state(.1, 373.15)
if97_state(23, 650)
}
