\name{iapws_epsilon}
\alias{iapws_epsilon}
\alias{iapws_n}
\alias{iapws_pk}
\title{
  IAPWS Formulations for Electrostatic Properties of Water
}
\description{
  Compute the dielectric constant, refractive index and ionization
  constant of water.
}
\usage{
iapws_epsilon(rho, t)
iapws_n(rho, t, lambda)
iapws_pk(rho, t)
}
\arguments{
  \item{rho}{a numeric vector giving the density values in kg/m\eqn{^3}.}
  \item{t}{a numeric vector giving the temperature values in K.}
  \item{lambda}{a numeric vector giving the wavelenth values in \eqn{\mu}m.}
}
\details{
  For calculating the properties as function of pressure, the function
  \code{iapws95()} should be employed.
}
\value{
  A numeric vector.
}
\references{
  International Association for the Properties of Water and Steam,
  IAPWS R8-97, \emph{Release on the Static Dielectric Constant
  of Ordinary Water Substance for Temperatures from 238 K to 873 K
  and Pressures up to 1000 MPa.}

  International Association for the Properties of Water and Steam,
  IAPWS R9-97, \emph{Release on the Refractive Index of Ordinary Water
  Substance as Function of Wavelength, Temperature and Pressure.} 

  International Association for the Properties of Water and Steam,
  IAPWS R11-2024, \emph{Revised Release on the Ionization Constant
  of H2O.}
}
\author{
  Jonathan Debove
}
\examples{
# Dielectric constant Tab. 4
t <- c(240, 300, 300, 300, 650, 650, 650, 870, 870, 870)
rho <- iapws95("rho", t = t,
               p = c(0.1013125, 0.1013125, 10, 1000,
                     10, 100, 500, 10, 100, 500),
               state = c(rep("liquid", 4),
                         "gas", rep("supercritical", 2),
                         "gas", rep("supercritical", 2)))
tab <- iapws_epsilon(rho, t)
print(tab, digits = 6)
\dontshow{
stopifnot(all.equal(tab, c(104.34982, 77.74735, 78.11269, 103.69632, 1.26715,
			   17.71733, 26.62132, 1.12721, 4.98281, 15.09746),
		    tolerance = 1e-7))
}

# Refractive index Tab. 3
d <- expand.grid(t = c(0, 100, 200, 500) + 273.15,
		 p = c(0.1, 1, 10, 100),
                 lambda = c(0.2265, 0.589, 1.01398))
d$state <- iapws95_state(p = d$p, t = d$t)
d$state[d$state == "solid"] <- "liquid"
d$rho <- drop(iapws95("rho", t = d$t, p = d$p, state = d$state))
tab <- iapws_n(d$rho, d$t, d$lambda)
dim(tab) <- c(4, 4, 3)
print(tab, digits = 8)

# Ionization constant Tab. 3
tab <- iapws_pk(rho = c(1, 0.07, 0.7, 0.2, 1.2, 0.0) * 1e3,
                t = c(300, 600, 600, 800, 800, 1270.0))
print(tab, digits = 7)
\dontshow{
stopifnot(all.equal(tab, c(13.906672, 20.161651, 11.147093, 14.487671,
			   6.4058649, 35.081557), tolerance = 1e-7))
}
}
