% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forestSuccession.R
\name{createRaster}
\alias{createRaster}
\title{createRaster:
This function creates a rasterized version of LiDAR data using the raster function
from the "raster" package}
\usage{
createRaster(data, resolution = 1)
}
\arguments{
\item{data}{LiDAR point cloud data}

\item{resolution}{Resolution of the raster}
}
\value{
A RasterLayer object
}
\description{
createRaster:
This function creates a rasterized version of LiDAR data using the raster function
from the "raster" package
}
\examples{
raster_data <- createRaster(myLiDARData, resolution = 1)
}
