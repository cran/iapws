\name{IAPWS}
\alias{IAPWS}
\alias{IAPWS95}
\alias{HEAVY17}
\docType{data}
\title{Constants built into iapws}
\description{
  \code{IAPWS} and \code{IAPWS95} contain constants relative to
  ordinary water substance.
  \code{HEAVY17} is the same but for heavy water.
}
\usage{
IAPWS
IAPWS95
HEAVY17
}
\format{
  \code{IAPWS}, \code{IAPWS95} and \code{HEAVY17} are numeric verctors
  with with the following named elements:
  \describe{
    \item{\code{"R"}}{molar gas constant in J/K/mol.}
    \item{\code{"M"}}{molar mass in g/mol.}
    \item{\code{"Pc"}}{critical pressure in MPa.}
    \item{\code{"Tc"}}{critical temperature in K.}
    \item{\code{"RHOc"}}{critical density in kg/m\eqn{^3}.}
    \item{\code{"Pt"}}{triple-point pressure in MPa.}
    \item{\code{"Tt"}}{triple-point temperature in K.}
  }
}
\examples{
# Convert critical density of water in mol/L
IAPWS["RHOc"] / IAPWS["M"]
}
