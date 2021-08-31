# R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
#
# data.R
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
# This script contains functions to handle the file IO for reading and writing
# data, i.e.:
#   read_csv(character fname, double y=-1)
#   read_txt(list my_data, character fname, character var, double y=-1)
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# - created read_csv and read_txt functions [15.02.23]
#
#### DEFINE FUNCTIONS ########################################################
# ************************************************************************
# Name:     read_csv
# Inputs:   - character, file name (fname)
#           - double, year (y)
# Returns:  list object (data)
#             $file_name ............ file name
#             $sf ................... sunshine fraction
#             $tair ................. air temperature
#             $pn ................... precipitation
#             $num_lines ............ number of data points
#             $year ................. year of data
# Features: Reads all three daily input variables (sf, tair, and pn)
#           for a single year from a CSV file that includes a header
# ************************************************************************
#' Read CSV file
#'
#' Reads all three daily input variables (sf, tair, and pn)
#' for a single year from a CSV file that includes a header.
#'
#' @importFrom utils read.csv
#'
#' @param fname String, file name.
#' @param y Numeric, year.
#'
#' @return List with the following properties:
#' \describe{
#'   \item{$file_name}{File name.}
#'   \item{$sf}{Sunshine fraction.}
#'   \item{$tair}{Air temperature.}
#'   \item{$pn}{Precipitation.}
#'   \item{$num_lines}{Number of data points.}
#'   \item{$year}{Year of data.}
#' }
#' @export
read_csv <- function(fname, y = -1) {
  my_data <- list()
  my_data$file_name <- fname

  # Read data from CSV file and save to return list:
  DATA <- read.csv(fname)
  my_data$sf <- DATA$sf
  my_data$tair <- DATA$tair
  my_data$pn <- DATA$pn
  my_data$num_lines <- dim(DATA)[1]

  my_year <- y
  if (y == -1) {
    if (dim(DATA)[1] == 366) {
      my_year <- 2000
    } else if (dim(DATA)[1] == 365) {
      my_year <- 2001
    }
  }
  my_data$year <- my_year
  return (my_data)
}


# ************************************************************************
# Name:     read_txt
# Inputs:   - list object (my_data)
#           - character, file name (fname)
#           - character, variable name (var)
#           - double, year (y)
# Returns:  list object (my_data)
# Features: Reads plain text file (no header) of one of the input
#           arrays
# ************************************************************************
#' Read plain text file
#'
#' Reads plain text file (no header) of one of the input arrays.
#'
#' @param my_data List same as the output from \code{\link{read_csv}}.
#' @param fname String, file name.
#' @param var String, variable name.
#' @param y Numeric, year.
#'
#' @return List with the following properties:
#' \describe{
#'   \item{$file_name}{File name.}
#'   \item{$sf}{Sunshine fraction.}
#'   \item{$tair}{Air temperature.}
#'   \item{$pn}{Precipitation.}
#'   \item{$num_lines}{Number of data points.}
#'   \item{$year}{Year of data.}
#' }
#' @export
read_txt <- function(my_data, fname, var, y = -1) {
  my_data$file_name <- c(my_data$file_name, fname)
  DATA <- scan(fname)
  if (var == "sf") {
    my_data$sf <- DATA
  } else if (var == "tair") {
    my_data$tair <- DATA
  } else if (var == "pn") {
    my_data$pn <- DATA
  }
  my_data$num_lines <- length(DATA)

  my_year <- y
  if (y == -1) {
    if (length(DATA) == 366) {
      my_year <- 2000
    } else if (length(DATA) == 365) {
      my_year <- 2001
    }
  }
  my_data$year <- my_year
  return (my_data)
}
