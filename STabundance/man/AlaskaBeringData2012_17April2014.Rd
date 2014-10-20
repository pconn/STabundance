\docType{data}
\name{AlaskaBeringData2012_17April2014}
\alias{AlaskaBeringData2012_17April2014}
\title{Alaska covariate data for spotted seal analyses}
\description{
  Bering Sea covariate data from 4 April 2012 through 14 May 2012 for use with spotted seal   analyses.  The list \code{Data} includes four objects:
  \itemize{
    \item{\code{Adj}} An adjacency matrix describing the spatial proximity of sampling units (sample units that are queen's move neighbors receive a 1, those that are not neighbors receive a 0)
    \item{\code{Adj2}} An adjacency matrix for a larger spatial extent (not used)
    \item{\code{Grid}} A list vector, each element of which provides a spatial polygons dataframe encoding the spatial position and covariates associated with sample units on each day.  For example, \code{Data$Grid[[1]]} holds covariate data for day 1 (4 April 2012).  Covariates provided include proportion land cover (\code{land_cover}), standardized distance from sample unit centroid to mainland (\code{dist_mainland}), standardized distance from land (including islands; \code{dist_land}), standardized distance from 1000m shelf break (\code{dist_shelf}), standardized depth (\code{depth}), proportion sea ice coverage (obtained from NSIDC; \code{ice_conc}), standardized distance from 90 percent sea ice contour (\code{dist_contour}), standardized distance from the southern sea ice edge (\code{dist_edge}), and Piatt and Springer ecoregion \code{Ecoregion}.
    \item{\code{Meta}} Provides some basic metadata
  }
}
\author{
  Paul Conn \email{paul.conn@noaa.gov}
}
\keyword{data}

