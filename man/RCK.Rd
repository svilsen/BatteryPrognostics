% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCKLookupTable.R
\name{RCK}
\alias{RCK}
\title{Simulate RC-K model}
\usage{
RCK(current, extractedTables, dt = 1, SOCStart = 0, trace = T,
  traceLimit = 10000)
}
\arguments{
\item{current}{A current profile.}

\item{extractedTables}{A list containing the possible 'SOC' and 'I' values, and the following look-up tables: 'R0', 'Rk', 'Ck', 'Cap', and 'OCV', where 'Rk' and 'Ck' are lists containing a look-up table for each RC-unit.}

\item{dt}{The simulation step-size.}

\item{SOCStart}{The starting SOC value.}

\item{trace}{TRUE/FALSE: Show trace?}

\item{traceLimit}{Only show the trace for the iterations of specified 'traceLimit'.}

\item{voltage}{A current profile reflecting the change in current (not the measured current).}
}
\value{
A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
}
\description{
Simulate battery voltage using look-up tables to approximate the parameters and assuming the Electrical Equivalent Circuit (EEC) model has 'K' RC-units.
}
