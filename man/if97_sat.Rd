\name{if97_sat}
\alias{if97_psat}
\alias{if97_tsat}
\title{
	IAPWS-IF97 Equations in Region 4
}
\description{
	Compute the pressure and the temperature along the saturated line
	according to the IAPWS industrial formulation 1997.
}
\usage{
if97_psat(t)
if97_tsat(p)
}
\arguments{
	\item{t}{a numeric vector giving the temperature values in K.}
	\item{p}{a numeric vector giving the pressure values in MPa.}
}
\value{
	A numeric vector containing the saturation-pressure (MPa) or the
	saturation-temperature (K). Return \code{NA} for inputs outside of the
	range of validity.
}
\references{
	International Association for the Properties of Water and Steam,
	IAPWS R7-97(2012), Revised Release on the IAPWS Industrial Formulation
	1997 for the Thermodynamic Properties of Water and Steam.
}
\author{
	Jonathan Debove
}
\examples{
# Region 4
format(if97_psat(t = c(300, 500, 600)), scientific = TRUE, digits = 9)
format(if97_tsat(p = c(0.1, 1.0, 10.0)), scientific = TRUE, digits = 9)
}
