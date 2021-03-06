% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{read_csv}
\alias{read_csv}
\title{Read CSV file}
\usage{
read_csv(fname, y = -1)
}
\arguments{
\item{fname}{String, file name.}

\item{y}{Numeric, year.}
}
\value{
List with the following properties:
\describe{
\item{$file_name}{File name.}
\item{$sf}{Sunshine fraction.}
\item{$tair}{Air temperature.}
\item{$pn}{Precipitation.}
\item{$num_lines}{Number of data points.}
\item{$year}{Year of data.}
}
}
\description{
Reads all three daily input variables (sf, tair, and pn)
for a single year from a CSV file that includes a header.
}
