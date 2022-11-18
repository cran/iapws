\name{iapws_sigma}
\alias{iapws_sigma}
\alias{heavy17_sigma}
\title{
  IAPWS Formulations for Surface Tension
}
\description{
  Compute the surface tension of ordinary and heavy water.
}
\usage{
iapws_sigma(t)
heavy17_sigma(t)
}
\arguments{
  \item{t}{a numeric vector giving the temperature values in K.}
}
\value{
  A numeric vector giving the values of surface tension in mN/m.
}
\references{
  International Association for the Properties of Water and Steam,
  IAPWS R1-76(2014), \emph{Revised Release on Surface Tension of Ordinary Water
  Substance.}

  International Association for the Properties of Water and Steam,
  IAPWS R5-85(1994), \emph{IAPWS Release on Surface Tension of Heavy Water
  Substance.}
}
\author{
  Jonathan Debove
}
\examples{
print(iapws_sigma(seq(5, 370, by = 5) + 273.15), digits = 2)
print(heavy17_sigma(seq(5, 370, by = 5) + 273.15), digits = 1)
}
