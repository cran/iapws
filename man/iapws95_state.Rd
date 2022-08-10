\name{iapws95_state}
\alias{iapws95_state}
\title{
	IAPWS Formulation 1995 Physical State
}
\description{
	Compute the physical state of water according the IAPWS formulation 1995.
}
\usage{
iapws95_state(p, t)
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
iapws95_state(.1, 293.15)
iapws95_state(.1, 373.15)
iapws95_state(23, 650)
}
