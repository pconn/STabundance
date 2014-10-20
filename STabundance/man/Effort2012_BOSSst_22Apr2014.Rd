\docType{data}
\name{Effort2012_BOSSst_22Apr2014.Rda}
\alias{Effort2012_BOSSst_22Apr2014.Rda}
\title{Bering Sea sampling effort as used in spotted seal analyses}
\description{
  Bering Sea sampling effort as used in spotted seal analyses (dataset last compiled 22 Apr    2014).  Surveys were conducted over U.S. portions of the Bering Sea during April and May of 2012, and are summarized here from April 10-May 8 (with the convention that April 10 is `Day 1' and May 8 is `Day 29.'  The R list object `Effort' includes four objects: 
  \itemize{
    \item{\code{Area.hab}}The proportion of each sample unit that is in the ocean (sample units that have some land receive a value less than 1.0)
    \item{\code{Area.trans}}For each combination of sample unit and day surveyed, 
    \code{Area.trans} indicates the proportion of the sample unit that was surveyed
    \item{\code{Count.data}}A data.frame providing the sample unit (column \code{Cell}) surveyed, the day (column \code{Time}) surveyed, the proportion of the grid cell that was surveyed (column \code{AreaSurveyed}), and the spotted seal count (column \code{Count}).
    \item{\code{Mapping}}A matrix indicating which sample unit (first column) was surveyed on which day (second column) and at which hour of day (UTC; third column).
  }
}
\author{
  Paul Conn \email{paul.conn@noaa.gov}
}
\keyword{data}

