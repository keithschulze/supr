% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segmentation.R
\name{remove_small_objects}
\alias{remove_small_objects}
\title{Remove connected regions/objects smaller than a specified size}
\usage{
remove_small_objects(image, size = 50, output_mask = FALSE)
}
\arguments{
\item{image}{binary image (instance of \code{\link{im}} object)}

\item{size}{Optional. Defaults to 50 units.
Exclude any objects smaller than this size. Units default to units of
the input \code{im}.}

\item{output_mask}{Optional. Logical value to indicate whether a mask of objects should
be returned. If \code{FALSE} a labelled image will be returned.}
}
\value{
labelled image (\code{\link{im}} object) or mask containing only
connected components > \code{size}.
}
\description{
Remove connected regions/objects smaller than a specified size
}

