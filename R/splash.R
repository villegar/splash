# R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
#
# splash.R
#
# VERSION: 1.0
# LAST UPDATED: 2016-02-19
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
# This script contains functions for running SPLASH for point-based data, i.e.:
#   spin_up(list mdat, list dtot)
#   quick_run(double lat, double elv, double n, double y, double wn, double sf,
#             double tc, double pn)
#   run_one_day(double lat, double elv, double n, double y, double wn,
#               double sf, double tc, double pn)
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# - added fix to daily soil moisture when n and ny are 365 [15.01.27]
#
#### IMPORT SOURCES ##########################################################
# source("const.R")
# source("evap.R")


#### DEFINE FUNCTIONS ########################################################
# ************************************************************************
# Name:     spin_up
# Inputs:   - list, meteorological data (mdat)
#               $num_lines ..... double, length of meteorol. variable lists
#               $lat_deg ....... double latitude (degrees)
#               $elv_m ......... double, elevation (m)
#               $year .......... double, year
#               $sf ............ list, fraction of sunshine hours
#               $tair .......... list, mean daily air temperature (deg. C)
#               $pn ............ list, precipitation (mm/d)
#           - list, daily totals (dtot)
#               $wm ............ list, daily soil moisture (mm)
# Returns:  list, daily totals
# Features: Updates the soil moisture in daily totals until equilibrium
# Depends:  quick_run
# ************************************************************************
#' Calculate daily totals
#'
#' Calculate daily totals updating the soil moisture until equilibrium.
#'
#' @details The list with meteorological data, \code{mdat}, should have the
#' following components:
#' \itemize{
#'  \item num_lines ..... double, length of meteorol. variable lists
#'  \item lat_deg ....... double latitude (degrees)
#'  \item elv_m ......... double, elevation (m)
#'  \item year .......... double, year
#'  \item sf ............ list, fraction of sunshine hours
#'  \item tair .......... list, mean daily air temperature (deg. C)
#'  \item pn ............ list, precipitation (mm/d)
#' }
#'
#' The list with daily totals, \code{dtot}, should have the following component:
#' \itemize{
#'  \item wm ............ list, daily soil moisture (mm)
#' }
#'
#' @param mdat list with meteorological data (see the details section).
#' @param dtot list with daily totals (see the details section).
#'
#' @return list, daily totals
#' @export
#'
#' @examples
#' daily_totals <- matrix(data = rep(0, 366), nrow = 366, ncol = 1)
#' daily_totals <- as.data.frame(daily_totals)
#' names(daily_totals) <- c("wn")
#' my_file <- system.file("extdata/example_data.csv", package = "splash")
#' my_data <- splash::read_csv(my_file, 2000)
#' my_data$lat_deg <- 37.7
#' my_data$elv_m <- 142
#' daily_totals <- splash::spin_up(my_data, daily_totals)
#' cat(sprintf("Spin-Up:\n"))
#' for (i in seq(from = 1, to = my_data$num_lines, by = 1)) {
#'   if (i == 1) cat(sprintf("Day\tWn (mm)\n"))
#'   cat(sprintf("%d\t%0.6f\n", i, daily_totals$wn[i]))
#' }
spin_up <- function(mdat, dtot) {
  # Run one year:
  for (i in seq(from = 1, to = mdat$num_lines, by = 1)) {
    if (i == 1) {
      wn <- dtot$wn[mdat$num_lines]
    } else {
      wn <- dtot$wn[i - 1]
    }
    my_vals <- quick_run(mdat$lat_deg,
                         mdat$elv_m,
                         i,
                         mdat$year,
                         wn,
                         mdat$sf[i],
                         mdat$tair[i],
                         mdat$pn[i])
    dtot$wn[i] <- my_vals$sm
  }

  # Calculate the change:
  start_sm <- dtot$wn[1]
  end_vals <- quick_run(mdat$lat_deg, mdat$elv_m, 1, mdat$year,
                        dtot$wn[mdat$num_lines], mdat$sf[1], mdat$tair[1],
                        mdat$pn[1])
  diff_sm <- abs(end_vals$sm - start_sm)

  # Equilibrate:
  spin_count <- 1
  while (diff_sm > 1.0) {
    for (i in seq(from = 1, to = mdat$num_lines, by = 1)) {
      if (i == 1) {
        wn <- dtot$wn[mdat$num_lines]
      } else {
        wn <- dtot$wn[i - 1]
      }
      my_vals <- quick_run(mdat$lat_deg, mdat$elv_m, i, mdat$year, wn,
                           mdat$sf[i], mdat$tair[i], mdat$pn[i])
      dtot$wn[i] <- my_vals$sm
    }

    # Calculate the change:
    start_sm <- dtot$wn[1]
    end_vals <- quick_run(mdat$lat_deg, mdat$elv_m, 1, mdat$year,
                          dtot$wn[mdat$num_lines], mdat$sf[1],
                          mdat$tair[1], mdat$pn[1])
    diff_sm <- abs(end_vals$sm - start_sm)
    spin_count <- spin_count + 1
  }
  cat(paste("Spun", spin_count, "years\n"))
  return(dtot)
}


# ************************************************************************
# Name:     quick_run
# Inputs:   - double, latitude, deg (lat)
#           - double, elevation, m (elv)
#           - double, day of year (n)
#           - double, year (y)
#           - double, daily soil moisture content, mm (wn)
#           - double, daily fraction of bright sunshine (sf)
#           - double, daily air temperature, deg C (tc)
#           - double, daily precipitation, mm (pn)
# Returns:  list
#             $sm - soil moisture, mm
#             $ro - runoff, mm
# Features: Returns daily soil moisture and runoff
# Depends:  evap
# ************************************************************************
#' Calculate daily soil moisture and runoff
#'
#' @inheritParams calc_daily_evap
#' @param pn double, daily precipitation, mm/day.
#' @param kCw double, supply constant, mm/hr.
#'     Default: \eqn{1.05} (Federer, 1982)
#' @param kWm double, soil moisture capacity, mm.
#'     Default: \eqn{150} (Cramer-Prentice, 1988)
#'
#' @return Returns daily soil moisture and runoff.
#' @keywords internal
#'
#' @references
#' Cramer, W. and Prentice, I.C., 1988. Simulation of regional soil moisture
#' deficits on a European scale. Norsk Geografisk Tidsskrift - Norwegian Journal
#' of Geography, 42(2-3), pp.149–151. doi:10.1080/00291958808552193
#'
#' Federer, C.A., 1982. Transpirational supply and demand: plant, soil, and
#' atmospheric effects evaluated by simulation. Water Resources Research, 18(2),
#' pp.355-362. doi:10.1029/WR018i002p00355
quick_run <- function(lat, elv, n, y, wn, sf, tc, pn, kCw = 1.05, kWm = 150) {
  # Calculate evaporative supply (mm/hr)
  sw <- kCw * wn / kWm

  # Compute daily radiation and evaporations values:
  ET <- calc_daily_evap(lat, n, elv, y, sf, tc, sw)

  # Update daily soil moisture:
  sm <- wn + pn + ET$cond_mm - ET$aet_mm

  if (sm > kWm) {
    # Bucket is full:
    # - set soil moisture to capacity
    # - add remaining water to runoff
    ro <- sm - kWm
    sm <- kWm
  } else if (sm < 0) {
    # Bucket is empty:
    # - set runoff and soil moisture equal to zero
    ro <- 0
    sm <- 0
  } else {
    ro <- 0
  }

  rval <- list()
  rval$sm <- sm
  rval$ro <- ro
  return(rval)
}


# ************************************************************************
# Name:     run_one_day
# Inputs:   - double, latitude, deg (lat)
#           - double, elevation, m (elv)
#           - double, day of year (n)
#           - double, year (y)
#           - double, daily soil moisture content, mm (wn)
#           - double, daily fraction of bright sunshine (sf)
#           - double, daily air temperature, deg C (tc)
#           - double, daily precipitation, mm (pn)
# Returns:  list
#             $ho - daily solar irradiation, J/m2
#             $hn - daily net radiation, J/m2
#             $ppfd - daily PPFD, mol/m2
#             $cond - daily condensation water, mm
#             $eet - daily equilibrium ET, mm
#             $pet - daily potential ET, mm
#             $aet - daily actual ET, mm
#             $wn - daily soil moisture, mm
#             $ro - daily runoff, mm
# Features: Runs SPLASH at a single location for one day.
# Depends:  evap
# ************************************************************************
#' Runs SPLASH at a single location for one day
#'
#' @param wn double, daily soil moisture content, mm (wn).
#' @inheritParams quick_run
#'
#' @return List with the following components:
#' \itemize{
#'  \item ho .......... daily solar irradiation, J/m2
#'  \item hn .......... daily net radiation, J/m2
#'  \item ppfd ........ daily PPFD, mol/m2
#'  \item cond ........ daily condensation water, mm
#'  \item eet ......... daily equilibrium ET, mm
#'  \item pet ......... daily potential ET, mm
#'  \item aet ......... daily actual ET, mm
#'  \item wn .......... daily soil moisture, mm
#'  \item ro .......... daily runoff, mm
#' }
#'
#' @export
#'
#' @references
#' Cramer, W. and Prentice, I.C., 1988. Simulation of regional soil moisture
#' deficits on a European scale. Norsk Geografisk Tidsskrift - Norwegian Journal
#' of Geography, 42(2-3), pp.149–151. doi:10.1080/00291958808552193
#'
#' Federer, C.A., 1982. Transpirational supply and demand: plant, soil, and
#' atmospheric effects evaluated by simulation. Water Resources Research, 18(2),
#' pp.355-362. doi:10.1029/WR018i002p00355
#'
#' @examples
#' soil <- run_one_day(lat = 37.7,
#'                     elv = 142,
#'                     n = 172,
#'                     y = 2000,
#'                     wn = 75,
#'                     sf = 1,
#'                     tc = 23,
#'                     pn = 5)
#' cat(sprintf("Soil moisture (run one day):\n"))
#' cat(sprintf("  Ho: %0.6f J/m2\n", soil$ho))
#' cat(sprintf("  Hn: %0.6f J/m2\n", soil$hn))
#' cat(sprintf("  PPFD: %0.6f mol/m2\n", soil$ppfd))
#' cat(sprintf("  EET: %0.6f mm/d\n", soil$eet))
#' cat(sprintf("  PET: %0.6f mm/d\n", soil$pet))
#' cat(sprintf("  AET: %0.6f mm/d\n", soil$aet))
#' cat(sprintf("  Cn: %0.6f mm/d\n", soil$cond))
#' cat(sprintf("  Wn: %0.6f mm\n", soil$wn))
#' cat(sprintf("  RO: %0.6f mm\n", soil$ro))
run_one_day <- function(lat, elv, n, y, wn, sf, tc, pn, kCw = 1.05, kWm = 150) {
  # Return values
  rvals <- list()

  # Calculate evaporative supply (mm/hr)
  sw <- kCw*wn/kWm

  # Compute daily radiation and evaporations values:
  ET <- calc_daily_evap(lat, n, elv, y, sf, tc, sw)
  rvals$ho <- ET$ra_j.m2
  rvals$hn <- ET$rn_j.m2
  rvals$ppfd <- ET$ppfd_mol.m2
  rvals$cond <- ET$cond_mm
  rvals$eet <- ET$eet_mm
  rvals$pet <- ET$pet_mm
  rvals$aet <- ET$aet_mm

  # Update daily soil moisture:
  sm <- wn + pn + ET$cond_mm - ET$aet_mm

  print( paste( "in run_one_day: pn =", pn))

  if (sm > kWm) {
    # Bucket is full:
    # - set soil moisture to capacity
    # - add remaining water to runoff
    ro <- sm - kWm
    sm <- kWm
  } else if (sm < 0) {
    # Bucket is empty:
    # - reduce actual ET by discrepancy amount
    # - set runoff and soil moisture equal to zero
    rvals$aet <- rvals$aet + sm
    ro <- 0
    sm <- 0
  } else {
    ro <- 0
  }

  rvals$wn <- sm
  rvals$ro <- ro
  return(rvals)
}
