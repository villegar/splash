# R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
#
# solar.R
#
# VERSION: 1.0-r2
# LAST UPDATED: 2016-08-19
#
# ~~~~~~~~
# license:
# ~~~~~~~~
# Copyright (C) 2016 Prentice Lab
#
# This file is part of the SPLASH model.
#
# SPLASH is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 2.1 of the License, or
# (at your option) any later version.
#
# SPLASH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SPLASH.  If not, see <http://www.gnu.org/licenses/>.
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
# algorithms for simulating habitats (SPLASH): Robust indices of radiation
# evapo-transpiration and plant-available moisture, Geoscientific Model
# Development, 2016 (in progress)
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script contains functions to calculate daily radiation, i.e.:
#   berger_tls(double n, double N)
#   density_h2o(double tc, double pa)
#   dcos(double d)
#   dsin(double d)
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# - fixed Cooper's and Spencer's declination angle equations [14.11.25]
# - replaced simplified_kepler with full_kepler method [14.11.25]
# - added berger_tls function [15.01.13]
# - updated evap function (similar to stash.py EVAP class) [15.01.13]
# - updated some documentation [16.05.27]
# - fixed HN- equation (iss#13) [16.08.19]
#
#### IMPORT SOURCES ##########################################################
# source("const.R")
# const()


#### DEFINE FUNCTIONS ########################################################

# ************************************************************************
# Name:     berger_tls
# Inputs:   - double, day of year (n)
#           - double, days in year (N)
# Returns:  numeric list, true anomaly and true longitude
# Features: Returns true anomaly and true longitude for a given day.
# Depends:  - ke ............. eccentricity of earth's orbit, unitless
#           - komega ......... longitude of perihelion
#  Ref:     Berger, A. L. (1978), Long term variations of daily insolation
#             and quaternary climatic changes, J. Atmos. Sci., 35, 2362-2367.
# ************************************************************************
#' Calculate true anomaly and true longitude
#'
#' @param n Numeric, day of year.
#' @param N Numeric, days in a year.
#' @inheritParams calc_daily_evap
#' @param pir \eqn{\pi} (\eqn{~ `r round(pi / 180, 6)`}) in radians.
#'
#' @return True anomaly and true longitude for a given day.
#' @keywords internal
#'
#' @references
#' Berger, A.L., 1978. Long-term variations of daily insolation and Quaternary
#' climatic changes. Journal of Atmospheric Sciences, 35(12), pp.2362-2367.
#' \doi{10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2}
berger_tls <- function(n,
                       N,
                       ke = 0.01670,
                       keps = 23.44,
                       komega = 283,
                       pir = pi / 180) {
  # Variable substitutes:
  xee <- ke ^ 2
  xec <- ke ^ 3
  xse <- sqrt(1 - ke ^ 2)

  # Mean longitude for vernal equinox:
  xlam <- (ke / 2.0 + xec / 8.0) * (1 + xse) * sin(komega * pir) -
    xee / 4.0 * (0.5 + xse) * sin(2.0 * komega * pir) +
    xec / 8.0 * (1.0 / 3.0 + xse) * sin(3.0 * komega * pir)
  xlam <- 2.0 * xlam / pir

  # Mean longitude for day of year:
  dlamm <- xlam + (n - 80.0) * (360.0 / N)

  # Mean anomaly:
  anm <- dlamm - komega
  ranm <- anm * pir

  # True anomaly (uncorrected):
  ranv <- ranm + (2.0 * ke - xec / 4.0) * sin(ranm) +
    5.0 / 4.0 * xee * sin(2.0 * ranm) +
    13.0 / 12.0 * xec * sin(3.0 * ranm)
  anv <- ranv / pir

  # True longitude:
  my_tls <- anv + komega
  if (my_tls < 0) {
    my_tls <- my_tls + 360
  } else if (my_tls > 360) {
    my_tls <- my_tls - 360
  }

  # True anomaly:
  my_nu <- my_tls - komega
  if (my_nu < 0) {
    my_nu <- my_nu + 360
  }
  return (c(my_nu, my_tls))
}


# ************************************************************************
# Name:     dcos
# Inputs:   double (d), angle in degrees
# Returns:  double, cosine of angle
# Features: This function calculates the cosine of an angle (d) given
#           in degrees.
# Depends:  pir
# Ref:      This script is based on the Javascript function written by
#           C Johnson, Theoretical Physicist, Univ of Chicago
#           - 'Equation of Time' URL: http://mb-soft.com/public3/equatime.html
#           - Javascript URL: http://mb-soft.com/believe/txx/astro22.js
# ************************************************************************
#' Calculate cosine of an angle
#'
#' Calculates the cosine of an angle (d) given in degrees.
#'
#' @param d Numeric, angle in degrees.
#'
#' @return Cosine of an angle.
#' @keywords internal
#'
#' @references
#' C. Johnson, Theoretical Physicist, Univ of Chicago
#'
#' - 'Equation of Time' URL: \url{https://mb-soft.com/public3/equatime.html}
#'
#' - Javascript URL: \url{https://mb-soft.com/believe/txx/astro22.js}
dcos <- function(d, pir = pi / 180) {
  cos(d * pir)
}


# ************************************************************************
# Name:     dsin
# Inputs:   double (d), angle in degrees
# Returns:  double, sine of angle
# Features: This function calculates the sine of an angle (d) given
#           in degrees.
# Depends:  pir
# ************************************************************************
#' Calculate sine of an angle
#'
#' Calculates the sine of an angle (d) given in degrees.
#'
#' @param d Numeric, angle in degrees.
#'
#' @return Sine of an angle.
#' @keywords internal
dsin <- function(d, pir = pi / 180) {
  sin(d * pir)
}


# ************************************************************************
# Name:     calc_daily_solar
# Inputs:   - double, latitude, degrees (lat)
#           - double, day of year (n)
#           - double, elevation (elv)  *optional
#           - double, year (y)         *optional
#           - double, fraction of sunshine hours (sf)        *optional
#           - double, mean daily air temperature, deg C (tc) *optional
# Returns:  list object (et.srad)
#             $nu_deg ............ true anomaly, degrees
#             $lambda_deg ........ true longitude, degrees
#             $dr ................ distance factor, unitless
#             $delta_deg ......... declination angle, degrees
#             $hs_deg ............ sunset angle, degrees
#             $ra_j.m2 ........... daily extraterrestrial radiation, J/m^2
#             $tau ............... atmospheric transmittivity, unitless
#             $ppfd_mol.m2 ....... daily photosyn. photon flux density, mol/m^2
#             $hn_deg ............ net radiation hour angle, degrees
#             $rn_j.m2 ........... daily net radiation, J/m^2
#             $rnn_j.m2 .......... daily nighttime net radiation, J/m^2
# Features: This function calculates daily radiation fluxes.
# Depends:  - kalb_sw ........ shortwave albedo
#           - kalb_vis ....... visible light albedo
#           - kb ............. empirical constant for longwave rad
#           - kc ............. empirical constant for shortwave rad
#           - kd ............. empirical constant for shortwave rad
#           - ke ............. eccentricity
#           - keps ........... obliquity
#           - kfFEC .......... from-flux-to-energy conversion, umol/J
#           - kGsc ........... solar constant
#           - berger_tls() ... calc true anomaly and longitude
#           - dcos() ......... cos(x*pi/180), where x is in degrees
#           - dsin() ......... sin(x*pi/180), where x is in degrees
#           - julian_day() ... date to julian day
# ************************************************************************
#' Calculate daily solar radiation fluxes
#'
#' This function calculates daily solar radiation fluxes.
#'
#' @inheritParams calc_daily_evap
#' @param kA double, empirical constant, degrees Celsius.
#'     Default: \eqn{107} (Monteith and Unsworth, 1990).
#' @param kalb_sw double, shortwave albedo.
#'     Default: \eqn{0.17} (Federer, 1968).
#' @param kalb_vis double, visible light albedo.
#'     Default: \eqn{0.03} (Sellers, 1985).
#' @param kb double, empirical constant.
#'     Default: \eqn{0.20} (Linacre, 1968).
#' @param kc double, cloudy transmittivity.
#'     Default: \eqn{0.25} (Linacre, 1968).
#' @param kd double, angular coefficient of transmittivity.
#'     Default: \eqn{0.50} (Linacre, 1968).
#' @param kfFEC double, flux-to-energy conversion, umol/J.
#'     Default: \eqn{2.04} (Meek et al., 1984).
#' @param kGsc double, solar constant, W/m^2.
#'     Default: \eqn{1360.8} (Kopp and Lean, 2011).
#'
#' @return Returns a \code{list} object with the following variables:
#' \itemize{
#'  \item nu_deg ............ true anomaly, degrees
#'  \item lambda_deg ........ true longitude, degrees
#'  \item dr ................ distance factor, unitless
#'  \item delta_deg ......... declination angle, degrees
#'  \item hs_deg ............ sunset angle, degrees
#'  \item ra_j.m2 ........... daily extraterrestrial radiation, J/m^2
#'  \item tau ............... atmospheric transmittivity, unitless
#'  \item ppfd_mol.m2 ....... daily photosyn. photon flux density, mol/m^2
#'  \item hn_deg ............ net radiation hour angle, degrees
#'  \item rn_j.m2 ........... daily net radiation, J/m^2
#'  \item rnn_j.m2 .......... daily nighttime net radiation, J/m^2
#' }
#' @export
#'
#' @references
#' Berger, A.L., 1978. Long-term variations of daily insolation and Quaternary
#' climatic changes. Journal of Atmospheric Sciences, 35(12), pp.2362-2367.
#' \doi{10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2}
#'
#' Federer, C.A., 1968. Spatial variation of net radiation, albedo and surface
#' temperature of forests. Journal of Applied Meteorology and Climatology, 7(5),
#' pp.789-795. \doi{10.1175/1520-0450(1968)007<0789:SVONRA>2.0.CO;2}
#'
#' Kopp, G. and Lean, J.L., 2011. A new, lower value of total solar irradiance:
#' Evidence and climate significance. Geophys. Res. Lett. 38, L01706.
#' \doi{10.1029/2010GL045777}
#'
#' Linacre, E.T., 1968. Estimating the net-radiation flux. Agricultural
#' meteorology, 5(1), pp.49-63. \doi{10.1016/0002-1571(68)90022-8}
#'
#' Meek, D.W., Hatfield, J.L., Howell, T.A., Idso, S.B. and Reginato, R.J.,
#' 1984. A generalized relationship between photosynthetically active radiation
#' and solar radiation 1. Agronomy journal, 76(6), pp.939-945.
#' \doi{10.2134/agronj1984.00021962007600060018x}
#'
#' Monteith, J., and Unsworth, M., 1990. Principles of Environmental Physics,
#' Butterworth-Heinemann, Oxford.
#'
#' Sellers, P.J., 1985. Canopy reflectance, photosynthesis and transpiration,
#' International Journal of Remote Sensing, 6:8, 1335-1372,
#' \doi{10.1080/01431168508948283}
#'
#' @examples
#' solar <- splash::calc_daily_solar(lat = 37.7,
#'                                   n = 172,
#'                                   elv = 142,
#'                                   y = 2000,
#'                                   sf = 1,
#'                                   tc = 23.0)
#' cat(sprintf("Solar values:\n"))
#' cat(sprintf("  kn: %d\n", solar$kN))
#' cat(sprintf("  nu: %0.6f degrees\n", solar$nu_deg))
#' cat(sprintf("  lambda: %0.6f degrees\n", solar$lambda_deg))
#' cat(sprintf("  rho: %0.6f\n", solar$rho))
#' cat(sprintf("  dr: %0.6f\n", solar$dr))
#' cat(sprintf("  delta: %0.6f degrees\n", solar$delta_deg))
#' cat(sprintf("  ru: %0.6f\n", solar$ru))
#' cat(sprintf("  rv: %0.6f\n", solar$rv))
#' cat(sprintf("  rw: %0.6f\n", solar$rw))
#' cat(sprintf("  hs: %0.6f degrees\n", solar$hs_deg))
#' cat(sprintf("  hn: %0.6f degrees\n", solar$hn_deg))
#' cat(sprintf("  tau_o: %0.6f\n", solar$tau_o))
#' cat(sprintf("  tau: %0.6f\n", solar$tau))
#' cat(sprintf("  Qn: %0.6f mol/m^2\n", solar$ppfd_mol.m2))
#' cat(sprintf("  Rnl: %0.6f w/m^2\n", solar$rnl_w.m2))
#' cat(sprintf("  Ho: %0.6f MJ/m^2\n", (1.0e-6) * solar$ra_j.m2))
#' cat(sprintf("  Hn: %0.6f MJ/m^2\n", (1.0e-6) * solar$rn_j.m2))
#' cat(sprintf("  Hnn: %0.6f MJ/m^2\n", (1.0e-6) * solar$rnn_j.m2))
calc_daily_solar <- function(lat,
                             n,
                             elv = 0,
                             y = 0,
                             sf = 1,
                             tc = 23.0,
                             ke = 0.01670,
                             keps = 23.44,
                             komega = 283,
                             kA = 107,
                             kalb_sw = 0.17,
                             kalb_vis = 0.03,
                             kb = 0.20,
                             kc = 0.25,
                             kd = 0.50,
                             kfFEC = 2.04,
                             kGsc = 1360.8) {
  # Local bindings
  pir <- pi / 180

  # ~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (lat > 90 || lat < -90) {
    stop("Warning: Latitude outside range of validity (-90 to 90)!")
  }
  if (n < 1 || n > 366) {
    stop("Warning: Day outside range of validity (1 to 366)!")
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  solar <- list()

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate the number of days in yeark (kN), days
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (y == 0) {
    kN <- 365
  } else {
    kN <- (julian_day(y + 1, 1, 1) - julian_day(y, 1, 1))
  }
  solar$kN <- kN

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Calculate heliocentric longitudes (nu and lambda), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  my_helio <- berger_tls(n = n, N = kN, ke = ke, keps = keps, komega = komega)
  nu <- my_helio[1]
  lam <- my_helio[2]
  solar$nu_deg <- nu
  solar$lambda_deg <- lam

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Calculate distance factor (dr), unitless
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Berger et al. (1993)
  kee <- ke ^ 2
  rho <- (1 - kee) / (1 + ke * dcos(nu))
  dr <- (1 / rho) ^ 2
  solar$rho <- rho
  solar$dr <- dr

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 04. Calculate the declination angle (delta), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Woolf (1968)
  delta <- asin(dsin(lam) * dsin(keps))
  delta <- delta / pir
  solar$delta_deg <- delta

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 05. Calculate variable substitutes (u and v), unitless
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ru <- dsin(delta) * dsin(lat)
  rv <- dcos(delta) * dcos(lat)
  solar$ru <- ru
  solar$rv <- rv

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 06. Calculate the sunset hour angle (hs), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Note: u/v equals tan(delta) * tan(lat)
  if (ru/rv >= 1.0) {
    hs <- 180  # Polar day (no sunset)
  } else if (ru / rv <= -1.0) {
    hs <- 0 # Polar night (no sunrise)
  } else {
    hs <- acos(-1.0 * ru / rv)
    hs <- hs / pir
  }
  solar$hs_deg <- hs

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 07. Calculate daily extraterrestrial radiation (ra_d), J/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq. 1.10.3, Duffy & Beckman (1993)
  ra_d <- (86400 / pi) * kGsc * dr * (ru * pir * hs + rv * dsin(hs))
  solar$ra_j.m2 <- ra_d

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 08. Calculate transmittivity (tau), unitless
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref:  Eq. 11, Linacre (1968); Eq. 2, Allen (1996)
  tau_o <- (kc + kd * sf)
  tau <- tau_o * (1 + (2.67e-5) * elv)

  solar$tau_o <- tau_o
  solar$tau <- tau

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 09. Calculate daily photosynthetic photon flux density (ppfd_d), mol/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ppfd_d <- (1e-6) * kfFEC * (1 - kalb_vis) * tau * ra_d
  solar$ppfd_mol.m2 <- ppfd_d

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 10. Estimate net longwave radiation (rnl), W/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rnl <- (kb + (1 - kb) * sf) * (kA - tc)
  solar$rnl_w.m2 <- rnl

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 11. Calculate variable substitue (rw), W/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rw <- (1 - kalb_sw) * tau * kGsc * dr
  solar$rw <- rw

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 12. Calculate net radiation cross-over angle (hn), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ((rnl - rw * ru)/(rw * rv) >= 1.0) {
    hn <- 0  # Net radiation is negative all day
  } else if ((rnl - rw * ru) / (rw * rv) <= -1.0) {
    hn <- 180 # Net radiation is positive all day
  } else {
    hn <- acos((rnl - rw * ru) / (rw * rv))
    hn <- hn / pir
  }
  solar$hn_deg <- hn

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 13. Calculate daytime net radiation (rn_d), J/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rn_d <- (86400 / pi) * (hn * pir * (rw * ru - rnl) + rw * rv * dsin(hn))
  solar$rn_j.m2 <- rn_d

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 14. Calculate nighttime net radiation (rnn_d), J/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # fixed iss#13
  rnn_d <- (86400 / pi) * (rw * rv * (dsin(hs) - dsin(hn)) +
                             rw * ru * (hs - hn) * pir -
                             rnl * (pi - hn * pir))
  solar$rnn_j.m2 <- rnn_d

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(solar)
}


# ************************************************************************
# Name:     julian_day
# Inputs:   - double, year (y)
#           - double, month (m)
#           - double, day of month (i)
# Returns:  double, Julian day
# Features: This function converts a date in the Gregorian calendar
#           to a Julian day number (i.e., a method of consecutative
#           numbering of days---does not have anything to do with
#           the Julian calendar!)
#           * valid for dates after -4712 January 1 (i.e., jde >= 0)
# Ref:      Eq. 7.1 J. Meeus (1991), Chapter 7 "Julian Day", Astronomical
#             Algorithms
# ************************************************************************
#' Calculate Julian day
#'
#' This function converts a date in the Gregorian calendar
#' to a Julian day number (i.e., a method of consecutive
#' numbering of days---does not have anything to do with
#' the Julian calendar!)
#'
#' * valid for dates after -4712 January 1 (i.e., jde >= 0)
#'
#' @param y double, year.
#' @param m double, month.
#' @param i double, day of month.
#'
#' @return double, Julian day.
#' @export
#'
#' @references
#' Meeus, J. 1991. Chapter 7 "Julian Day". Astronomical Algorithms.
#' Willmann-Bell.
julian_day <- function(y, m, i) {
  if (m <= 2) {
    y <- y - 1
    m <- m + 12
  }
  a <- floor(y/100)
  b <- 2 - a + floor(a/4)

  jde <- floor(365.25*(y + 4716)) + floor(30.6001*(m + 1)) + i + b - 1524.5
  return(jde)
}
