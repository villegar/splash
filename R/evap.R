# R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
#
# evap.R
#
# VERSION: 1.0-r2
# LAST UPDATED: 2016-09-11
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
# This script contains functions to calculate daily radiation, condensation,
# and evapotranspiration, i.e.:
#   berger_tls(double n, double N)
#   density_h2o(double tc, double pa)
#   dcos(double d)
#   dsin(double d)
#   elv2pres(double z)
#   enthalpy_vap(double tc)
#   evap(double lat, double n, double elv=0, double y=0, double sf=1,
#        double tc=23.0, double sw=1.0)
#   julian_day(double y, double m, double i)
#   psychro(double tc, double pa)
#   sat_slope(double tc)
#   specific_heat(double tc)
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# - fixed Cooper's and Spencer's declination angle equations [14.11.25]
# - replaced simplified_kepler with full_kepler method [14.11.25]
# - added berger_tls function [15.01.13]
# - updated evap function (similar to stash.py EVAP class) [15.01.13]
# - added missing variables to evap list [16.02.17]
# - updated documentation [16.05.27]
# - addressed specific heat limitation [16.09.11]
#
#### IMPORT SOURCES ##########################################################
# source("const.R")
# source("solar.R")

#### DEFINE FUNCTIONS ########################################################
# ************************************************************************
# Name:     density_h2o
# Inputs:   - double (tc), air temperature, degrees C
#           - double (pa), atm pressure, Pa
# Returns:  double, kg/m^3
# Features: This function calculates the temperature and pressure
#           dependent density of pure water
# * Ref:    Chen, C.T., R.A. Fine, and F.J. Millero (1977), The equation
#             of state of pure water determined from sound speeds, The
#             Journal of Chemical Physics 66, 2142;
#             \doi{10.1063/1.434179}
# ************************************************************************
#' Calculate density of water at 1 atm, g/cm^3
#'
#' This function calculates the temperature and pressure dependent density of
#' pure water.
#'
#' @param tc double, air temperature, degrees C.
#' @param pa double, atm pressure, Pa.
#'
#' @return double, kg/m^3.
#' @keywords internal
#'
#' @references
#' Chen, C.T., Fine, R.A. and Millero, F.J., 1977. The equation of state of
#' pure water determined from sound speeds. The Journal of Chemical Physics,
#' 66(5), pp.2142-2144. \doi{10.1063/1.434179}
density_h2o <- function(tc, pa) {
  # Calculate density of water at 1 atm, g/cm^3
  po <- 0.99983952 +
    (6.788260e-5)  * tc +
    -(9.08659e-6)  * tc * tc +
    (1.022130e-7)  * tc * tc * tc +
    -(1.35439e-9)  * tc * tc * tc * tc +
    (1.471150e-11) * tc * tc * tc * tc * tc +
    -(1.11663e-13) * tc * tc * tc * tc * tc * tc +
    (5.044070e-16) * tc * tc * tc * tc * tc * tc * tc +
    -(1.00659e-18) * tc * tc * tc * tc * tc * tc * tc * tc

  # Calculate the bulk modulus of water at 1 atm, atm
  ko <- 19652.17  +
    148.1830      * tc +
    -2.29995      * tc * tc +
    0.01281       * tc * tc * tc +
    -(4.91564e-5) * tc * tc * tc * tc +
    (1.035530e-7) * tc * tc * tc * tc * tc

  # Calculate temperature-dependend coefficients
  ca <- 3.26138 +
    (5.223e-4)  * tc +
    (1.324e-4)  * tc * tc +
    -(7.655e-7) * tc * tc * tc +
    (8.584e-10) * tc * tc * tc * tc

  cb <- (7.2061e-5) +
    -(5.8948e-6)    * tc +
    (8.69900e-8)    * tc * tc +
    -(1.0100e-9)    * tc * tc * tc +
    (4.3220e-12)    * tc * tc * tc * tc

  # Convert pressure to bar (1 bar = 100000 Pa)
  pbar <- (1e-5) * pa

  pw <-
    (1e3) * po * (ko + ca * pbar + cb * pbar ^ 2) / (ko + ca * pbar + cb *
                                                       pbar ^ 2 - pbar)
  return(pw)
}


# ************************************************************************
# Name:     elv2pres
# Inputs:   double (z), meters
# Returns:  double, Pa
# Features: Calculates atmospheric pressure for a given elevation
# Depends:  - kPo ............ base pressure, Pa
#           - kTo ............ base temperature, C
#           - kL ............. temperature lapse rate, K/m
#           - kR ............. universal gas constant, J/mol/K
#           - kMa ............ molecular weight of dry air, kg/mol
#           - kG ............. gravity, m/s^2
# Ref:      Allen et al. (1998)
# ************************************************************************
#' Elevation to pressure
#'
#' Calculates atmospheric pressure for a given elevation.
#'
#' @param z Elevation, m.
#' @param kG Gravitational acceleration, m/s^2.
#'     Default: \eqn{9.80665} (Allen, 1973)
#' @param kL Adiabatic lapse rate, K/m.
#'     Default: \eqn{0.0065} (Allen, 1973)
#' @param kMa Molecular weight of dry air, kg/mol.
#'     Default: \eqn{0.028963} (Tsilingiris, 2008)
#' @param kPo Standard atmosphere, Pa.
#'     Default: \eqn{101325} (Allen, 1973)
#' @param kR Universal gas constant, J/mol/K.
#'     Default: \eqn{8.31447} (Moldover et al., 1988)
#' @param kTo Base temperature, K.
#'     Default: \eqn{288.15} (Berberan-Santos et al., 1997)
#'
#' @return Atmospheric pressure for the given elevation, Pa.
#' @keywords internal
#'
#' @references
#' Allen, R.G., Pereira, L.S., Raes, D. and Smith, M., 1998. Crop
#' evapotranspiration-Guidelines for computing crop water requirements-FAO
#' Irrigation and drainage paper 56. Food and Agriculture Organization of the
#' United Nations, Rome, 300(9), p.D05109. Available:
#' http://www.fao.org/docrep/x0490e/x0490e07.htm
#'
#' Berberan-Santos, M.N., Bodunov, E.N. and Pogliani, L., 1997. On the
#' barometric formula. American Journal of Physics, 65(5), pp.404-412.
#' \doi{10.1119/1.18555}
#'
#' Moldover, M.R., Trusler, J.M., Edwards, T.J., Mehl, J.B. and Davis, R.S.,
#' 1988. Measurement of the universal gas constant R using a spherical acoustic
#' resonator. Physical review letters, 60(4), p.249.
#' \doi{10.1103/PhysRevLett.60.249}
#'
#' Tsilingiris, P.T., 2008. Thermophysical and transport properties of humid air
#' at temperature range between 0 and 100 C. Energy Conversion and Management,
#' 49(5), pp.1098-1110. \doi{10.1016/j.enconman.2007.09.015}
elv2pres <- function(z,
                     kG = 9.80665,
                     kL = 0.0065,
                     kMa = 0.028963,
                     kPo = 101325,
                     kR = 8.31447,
                     kTo = 288.15) {
  kPo * (1 - kL * z / kTo) ^ (kG * kMa / (kR * kL))
}


# ************************************************************************
# Name:     enthalpy_vap
# Inputs:   double (tc), air temperature, degrees C
# Returns:  double, J/kg
# Features: This function calculates the temperature-dependent enthalpy
#           of vaporization (latent heat of vaporization)
# Ref:      Eq. 8, Henderson-Sellers (1984), A new formula for latent heat
#             of vaporization of water as a function of temperature, Quarterly
#             Journal of the Royal Meteorological Society, vol. 110, pp. 1186--
#             1190.
# ************************************************************************
#' Calculate enthalpy of vaporization
#'
#' This function calculates the temperature-dependent enthalpy of vaporization
#' (latent heat of vaporization).
#'
#' @inheritParams density_h2o
#'
#' @return double, J/kg.
#' @keywords internal
#'
#' @references
#' Eq. 8, Henderson‐Sellers, B., 1984. A new formula for latent heat of
#' vaporization of water as a function of temperature. Quarterly Journal of the
#' Royal Meteorological Society, 110(466), pp.1186-1190.
#' \doi{10.1002/qj.49711046626}
enthalpy_vap <- function(tc) {
  1.91846e6 * ((tc + 273.15) / (tc + 273.15 - 33.91)) ^ 2
}


# ************************************************************************
# Name:     psychro
# Inputs:   - double (tc), air temperature, degrees C
#           - double (pa), atm pressure, Pa
# Returns:  double, Pa/K
# Features: This function calculates the temperature and pressure
#           dependent psychrometric constant
# Depends:  - enthalpy_vap
#           - specific_heat
# Ref:      Allen, R.G., L.S. Pereira, D. Raes, M. Smith (1998),
#           'Meteorological data,' Crop evapotranspiration - Guidelines
#           for computing crop water requirements - FAO Irrigation and
#           drainage paper 56, Food and Agriculture Organization of the
#           United Nations, Available:
#           http://www.fao.org/docrep/x0490e/x0490e07.htm
# ************************************************************************
#' Calculate psychrometric constant
#'
#' This function calculates the temperature and pressure dependent
#' psychrometric constant.
#'
#' @inheritParams density_h2o
#' @param kMa double, molecular weight of dry air, kg/mol.
#'     Default: \eqn{0.028963} (Tsilingiris, 2008)
#' @param kMv double, molecular weight of water vapor, kg/mol.
#'     Default: \eqn{0.01802} (Tsilingiris, 2008)
#' @return double, Pa/K.
#' @keywords internal
#'
#' @references
#' Allen, R.G., Pereira, L.S., Raes, D. and Smith, M., 1998. Crop
#' evapotranspiration-Guidelines for computing crop water requirements-FAO
#' Irrigation and drainage paper 56. Food and Agriculture Organization of the
#' United Nations, Rome, 300(9), p.D05109. Available:
#' http://www.fao.org/docrep/x0490e/x0490e07.htm
psychro <- function(tc, pa, kMa = 0.028963, kMv = 0.01802) {
  # Calculate the specific heat capacity of water, J/kg/K
  cp <- specific_heat(tc)

  # Calculate latent heat of vaporization, J/kg
  lv <- enthalpy_vap(tc)

  # Calculate psychrometric constant, Pa/K
  return(cp * kMa * pa / (kMv * lv))
}


# ************************************************************************
# Name:     specific_heat
# Inputs:   double (tc), air temperature, degrees C
# Returns:  double, specific heat of moist air, J/kg/K
# Features: This function calculates the spefic heat of moist air
# Ref:      Tsilingris (2008), Thermophysical and transport properties of
#           humid air at temperature range between 0 and 100 °C, Energy
#           Conversion and Management, vol. 49, pp. 1098--1110.
# ************************************************************************
#' Calculate specific heat
#'
#' This function calculates the specific heat of moist air.
#'
#' @inheritParams density_h2o
#'
#' @return double, specific heat of moist air, J/kg/K.
#' @keywords internal
#'
#' @references
#' Tsilingiris, P.T., 2008. Thermophysical and transport properties of humid air
#' at temperature range between 0 and 100 C. Energy Conversion and Management,
#' 49(5), pp.1098-1110. \doi{10.1016/j.enconman.2007.09.015}
specific_heat <- function(tc) {
  if (tc < 0) {
    tc <- 0
  } else if (tc > 100) {
    tc <- 100
  }
  cp <- 1.0045714270  +
    (2.050632750e-3)  * tc -
    (1.631537093e-4)  * tc * tc +
    (6.212300300e-6)  * tc * tc * tc -
    (8.830478888e-8)  * tc * tc * tc * tc +
    (5.071307038e-10) * tc * tc * tc * tc * tc
  cp <- (1e3) * cp
  return(cp)
}


# ************************************************************************
# Name:     sat_slope
# Inputs:   double (tc), degrees C
# Returns:  double, Pa/K
# Features: This function calculates the temperature-dependent slope of
#           the saturation pressure temperature curve using the
#           methodology presented in the eMast energy.cpp script
# Ref:      - Eq. 6, Prentice et al. (1993);
#           - Eq. 13, Allen et al. (1998)
# ************************************************************************
#' Calculate the temperature-dependent slope
#'
#' This function calculates the temperature-dependent slope of the saturation
#' pressure temperature curve using the methodology presented in the eMast
#' energy.cpp script.
#' @inheritParams density_h2o
#'
#' @return double, Pa/K.
#' @keywords internal
#'
#' @references
#' Allen, R.G., Pereira, L.S., Raes, D. and Smith, M., 1998. Crop
#' evapotranspiration-Guidelines for computing crop water requirements-FAO
#' Irrigation and drainage paper 56. Food and Agriculture Organization of the
#' United Nations, Rome, 300(9), p.D05109. Available:
#' http://www.fao.org/docrep/x0490e/x0490e07.htm
#'
#' Prentice, I.C., Sykes, M.T. and Cramer, W., 1993. A simulation model for the
#' transient effects of climate change on forest landscapes. Ecological
#' modelling, 65(1-2), pp.51-70. \doi{10.1016/0304-3800(93)90126-D}
sat_slope <- function(tc) {
  (17.269) * (237.3) * (610.78) *
    exp(17.269 * tc / (237.3 + tc)) / (237.3 + tc) ^ 2
}

# ************************************************************************
# Name:     calc_daily_evap
# Inputs:   - double, latitude, degrees (lat)
#           - double, day of year (n)
#           - double, elevation (elv)  *optional
#           - double, year (y)         *optional
#           - double, fraction of sunshine hours (sf)        *optional
#           - double, mean daily air temperature, deg C (tc) *optional
#           - double, evaporative supply rate, mm/hr (sw)    *optional
# Returns:  list object (evap)
#             $nu_deg ............ true anomaly, degrees
#             $lambda_deg ........ true longitude, degrees
#             $dr ................ distance factor, unitless
#             $delta_deg ......... declination angle, degrees
#             $hs_deg ............ sunset angle, degrees
#             $ra_j.m2 ........... daily extraterrestrial radiation, J/m^2
#             $tau ............... atmospheric transmittivity, unitless
#             $ppfd_mol.m2 ....... daily photosyn photon flux density, mol/m^2
#             $hn_deg ............ net radiation hour angle, degrees
#             $rn_j.m2 ........... daily net radiation, J/m^2
#             $rnn_j.m2 .......... daily nighttime net radiation, J/m^2
#             $econ_m3.j ......... water to energy conversion, m^3/J
#             $cond_mm ........... daily condensation, mm
#             $eet_mm ............ daily equilibrium evapotranspiration, mm
#             $pet_mm ............ daily potential evapotranspiration, mm
#             $hi_deg ............ intersection hour angle, degrees
#             $aet_mm ............ daily actual evapotranspiration, mm
# Features: This function calculates daily radiation, condensation, and
#           evaporation fluxes.
# Depends:  - kw ............. entrainment factor for PET
#           - calc_daily_solar daily radiation fluxes
#           - dcos() ......... cos(x*pi/180), where x is in degrees
#           - dsin() ......... sin(x*pi/180), where x is in degrees
#           - density_h2o() .. density of water
#           - elv2pres() ..... elevation dependent atm. pressure
#           - enthalpy_vap() . latent heat of vaporization
#           - psychro() ...... psychrometric constant
#           - sat_slope() .... slope of sat. pressure temp curve
# ************************************************************************
#' Calculate daily evaporation fluxes
#'
#' This function calculates daily radiation, condensation, and evaporation
#' fluxes.
#'
#' @param lat double, decimal degrees.
#' @param n double, day of year.
#' @param elv double, elevation, m A.S.L.
#'     Default: \eqn{0}.
#' @param y double, year.
#'     Default: \eqn{0}.
#' @param sf double, fraction of sunshine hours.
#'     Default: \eqn{1}.
#' @param tc double, mean daily air temperature, degrees C.
#'     Default: \eqn{23.0}.
#' @param sw double, evaporative supply rate, mm/hr.
#'     Default: \eqn{1.0}.
#' @param ke double, eccentricity of earth's orbit.
#'     Default: \eqn{0.01670}, 2000CE (Berger, 1978).
#' @param keps double, obliquity of earth's elliptic.
#'     Default: \eqn{23.44}, 2000CE (Berger, 1978).
#' @param komega double, lon. of perihelion, degrees
#'     Default: \eqn{283}, 2000CE (Berger, 1978).
#' @param kw double, PET entrainment, \eqn{(1 + kw) * EET}
#'     Default: \eqn{0.26} (Priestley-Taylor, 1972)
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
#'  \item ppfd_mol.m2 ....... daily photosyn photon flux density, mol/m^2
#'  \item hn_deg ............ net radiation hour angle, degrees
#'  \item rn_j.m2 ........... daily net radiation, J/m^2
#'  \item rnn_j.m2 .......... daily nighttime net radiation, J/m^2
#'  \item econ_m3.j ......... water to energy conversion, m^3/J
#'  \item cond_mm ........... daily condensation, mm
#'  \item eet_mm ............ daily equilibrium evapotranspiration, mm
#'  \item pet_mm ............ daily potential evapotranspiration, mm
#'  \item hi_deg ............ intersection hour angle, degrees
#'  \item aet_mm ............ daily actual evapotranspiration, mm
#' }
#' @export
#'
#' @references
#' Berger, A.L., 1978. Long-term variations of daily insolation and Quaternary
#' climatic changes. Journal of Atmospheric Sciences, 35(12), pp.2362-2367.
#' \doi{10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2}
#'
#' Priestley, C.H.B. and Taylor, R.J., 1972. On the assessment of surface heat
#' flux and evaporation using large-scale parameters. Monthly weather review,
#' 100(2), pp.81-92. \doi{10.1175/1520-0493(1972)100<0081:OTAOSH>2.3.CO;2}
#'
#' @examples
#' evap <- splash::calc_daily_evap(lat = 37.7,
#'                                 n = 172,
#'                                 elv = 142,
#'                                 y = 2000,
#'                                 sf = 1,
#'                                 tc = 23.0,
#'                                 sw = 0.9)
#' cat(sprintf("Evaporation values:\n"))
#' cat(sprintf("  s: %0.6f Pa/K\n", evap$s_pa.k))
#' cat(sprintf("  Lv: %0.6f MJ/kg\n", (1e-6) * evap$lv_j.kg))
#' cat(sprintf("  Patm: %0.6f bar\n", (1e-5) * evap$patm_pa))
#' cat(sprintf("  pw: %0.6f kg/m^3\n", evap$pw_kg.m3))
#' cat(sprintf("  gamma: %0.6f Pa/K\n", evap$gam_pa.k))
#' cat(sprintf("  Econ: %0.6f mm^3/J\n", (1e9) * evap$econ_m3.j))
#' cat(sprintf("  Cn: %0.6f mm\n", evap$cond_mm))
#' cat(sprintf("  rx: %0.6f\n", evap$rx))
#' cat(sprintf("  hi: %0.6f degrees\n", evap$hi_deg))
#' cat(sprintf("  EET: %0.6f mm\n", evap$eet_mm))
#' cat(sprintf("  PET: %0.6f mm\n", evap$pet_mm))
#' cat(sprintf("  AET: %0.6f mm\n", evap$aet_mm))
calc_daily_evap <-function(lat,
                           n,
                           elv = 0,
                           y = 0,
                           sf = 1,
                           tc = 23.0,
                           sw = 1.0,
                           ke = 0.01670,
                           keps = 23.44,
                           komega = 283,
                           kw = 0.26) {
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
  evap <- list()

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 1. Calculate radiation fluxes
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  solar <- calc_daily_solar(lat = lat,
                            n = n,
                            elv = elv,
                            y = y,
                            sf = sf,
                            tc = tc,
                            ke = ke,
                            keps = keps,
                            komega = komega)
  ru <- solar$ru;
  rv <- solar$rv;
  rw <- solar$rw;
  rnl <- solar$rnl_w.m2;
  hn <- solar$hn_deg;
  rn_d <- solar$rn_j.m2;
  rnn_d <- solar$rnn_j.m2;
  evap$ra_j.m2  <- solar$ra_j.m2
  evap$rn_j.m2 <- solar$rn_j.m2
  evap$ppfd_mol.m2 <- solar$ppfd_mol.m2

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2. Calculate water-to-energy conversion (econ), m^3/J
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Atmospheric pressure, Pa
  patm <- elv2pres(elv)
  evap$patm_pa <- patm

  # Slope of saturation vap press temp curve, Pa/K
  s <- sat_slope(tc)
  evap$s_pa.k <- s

  # Enthalpy of vaporization, J/kg
  lv <- enthalpy_vap(tc)
  evap$lv_j.kg <- lv

  # Density of water, kg/m^3
  pw <- density_h2o(tc, patm)
  evap$pw_kg.m3 <- pw

  # Psychrometric constant, Pa/K
  gam <- psychro(tc, patm)
  evap$gam_pa.k <- gam

  econ <- s / (lv * pw * (s + gam))
  evap$econ_m3.j <- econ

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 3. Calculate daily condensation (cn), mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cn <- (1e3) * econ * abs(rnn_d)
  evap$cond_mm <- cn

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 4. Estimate daily equilibrium evapotranspiration (eet_d), mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eet_d <- (1e3) * econ * rn_d
  evap$eet_mm <- eet_d

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 5. Estimate daily potential evapotranspiration (pet_d), mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pet_d <- (1 + kw) * eet_d
  evap$pet_mm <- pet_d

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 6. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rx <- (3.6e6) * (1 + kw) * econ
  evap$rx <- rx

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 7. Calculate the intersection hour angle (hi), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cos_hi <- sw / (rw * rv * rx) + rnl / (rw * rv) - ru / rv
  #     print( paste( "in evap: sw =", sw))
  #     print( paste( "in evap: ru =", ru))
  #     print( paste( "in evap: rv =", rv))
  #     print( paste( "in evap: rw =", rw))
  #     print( paste( "in evap: rx =", rx))
  #     print( paste( "in evap: rnl =", rnl))
  #     print( paste( "in evap: cos_hi =", cos_hi))
  if (cos_hi >= 1.0) {
    hi <- 0.0       # supply exceeds demand
  } else if (cos_hi <= -1.0) {
    hi <- 180.0     # supply limits demand everywhere
  } else {
    hi <- acos(cos_hi)
    hi <- hi / pir
  }
  evap$hi_deg <- hi

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 8. Estimate daily actual evapotranspiration (aet_d), mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  aet_d <- (24 / pi) * (
    sw * hi * pir +
      rx * rw * rv * (dsin(hn) - dsin(hi)) +
      (rx * rw * ru - rx * rnl) * (hn - hi) * pir
  )
  evap$aet_mm <- aet_d

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(evap)
}
