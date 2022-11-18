\name{if97}
\alias{if97}
\title{
  IAPWS-IF97 Equations for Regions 1, 2, 3 and 5
}
\description{
  Compute thermodynamic and transport properties of water and steam
  using the IAPWS industrial formulation 1997.
}
\usage{
if97(what, p, t, h, state = NULL)
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
  \item{h}{
    a numeric vector giving the specific enthalpy values in kJ/kg.
  }
  \item{state}{
    \code{NULL} or a character vector giving the physical state.
    One of \code{"liquid"}, \code{"gas"} or \code{"supercritical"}.
    In case of \code{NULL}, the state is determined by the function
    \code{if97_state()}.
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

  A valid pair of \code{p}, \code{t} and \code{h} is needed:
  either \code{(p, t)} or \code{(p, h)}.
}
\value{
  A numeric matrix of dimension \code{c(n, length(what))} with \code{n}
  the maximum length of either \code{p}, \code{t}, \code{h} or \code{state}.
}
\references{
  International Association for the Properties of Water and Steam,
  IAPWS R7-97(2012), \emph{Revised Release on the IAPWS Industrial
  Formulation 1997 for the Thermodynamic Properties of Water and Steam.}

  International Association for the Properties of Water and Steam,
  IAPWS R12-08(2008), \emph{Release on the IAPWS Formulation 2008
  for the Viscosity of Ordinary Water Substance.}

  International Association for the Properties of Water and Steam,
  IAPWS R15-11(2011), \emph{Release on the IAPWS Formulation 2011
  for the Thermal Conductivity of Ordinary Water Substance.}
}
\author{
  Jonathan Debove
}
\note{
  Computing several properties in a single call may be more efficient
  than separately because most of the computation time is shared.
}
\examples{
what <- c("v", "h", "u", "s", "cp", "w")

# Region 1 Tab. 5
tab5 <- if97(what, t = c(300, 300, 500), p = c(3, 80, 3))
format(tab5, scientific = TRUE, digits = 9)

# Region 2 Tab. 15
tab15 <- if97(what, t = c(300, 700, 700), p = c(3.5e-3, 3.5e-3, 30))
format(tab15, scientific = TRUE, digits = 9)

# Region 2 metastable Tab. 18
tab18 <- if97(what, t = c(450, 440, 450), p = c(1, 1, 1.5), state = "gas")
format(tab18, scientific = TRUE, digits = 9)

# Region 3 Tab. 33
tab33 <- if97(what, t = c(650, 650, 750),
              p = c(0.255837018e2, 0.222930643e2, 0.783095639e2))
format(tab33, scientific = TRUE, digits = 9)

# Region 5 Tab. 42
tab42 <- if97(what, t = c(1500, 1500, 2000), p = c(0.5, 30, 30))
format(tab42, scientific = TRUE, digits = 9)

\dontshow{
stopifnot(all.equal(tab5[, "v", drop = FALSE],
          matrix(c(0.100215168e-2, 0.971180894e-3, 0.120241800e-2),
                 3, 1, dimnames = list(NULL, "v")),
          tolerance = 2e-9),
          all.equal(tab15[, "v", drop = FALSE],
          matrix(c(0.394913866e2, 0.923015898e2, 0.542946619e-2),
                 3, 1, dimnames = list(NULL, "v")),
          tolerance = 2e-9),
          all.equal(tab18[, "v", drop = FALSE],
          matrix(c(0.192516540, 0.186212297, 0.121685206),
                 3, 1, dimnames = list(NULL, "v")),
          tolerance = 2e-9),
          all.equal(tab33[, "v", drop = FALSE],
          matrix(1 / c(500, 200, 500),
                 3, 1, dimnames = list(NULL, "v")),
          tolerance = 2e-8),
          all.equal(tab42[, "v", drop = FALSE],
          matrix(c(0.138455090e1, 0.230761299e-1, 0.311385219e-1),
                 3, 1, dimnames = list(NULL, "v")),
          tolerance = 2e-9))
}

# Viscosity and thermal conductivity
what <- c("lambda", "eta")

# Region 1 Tab. 7
tab7 <- if97(what, p = c(20, 50), t = c(620, 620))
format(tab7, scientific = TRUE, digits = 9)

# Region 2 Tab. 8
tab8 <- if97(what, p = c(0.3, 50), t = c(650, 800))
format(tab8, scientific = TRUE, digits = 9)

## Region 3 Tab. 9
#tab9 <- if97_rhot(c("lambda", "rho", "cp", "cv", "eta"),
#      rho = c(0.3, 50), t = c(222, 322)),
#format(tab9, scientific = TRUE, digits = 9)

\dontshow{
stopifnot(all.equal(tab7[, "lambda", drop = FALSE],
          matrix(c(0.481485195e3, 0.545038940e3),
                 2, 1, dimnames = list(NULL, "lambda")),
          tolerance = 2e-9),
          all.equal(tab8[, "lambda", drop = FALSE],
          matrix(c(0.522311024e2, 0.177709914e3),
                 2, 1, dimnames = list(NULL, "lambda")),
          tolerance = 3e-9))
}
}
