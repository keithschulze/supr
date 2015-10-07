#' Measure properties of different regions in a labelled image.
#'
#' Given a labelled image (e.g., result of \code{\link{connected}}),
#' extract properties of different regions.
#'
#' @seealso \code{\link{connected}}
#'          \code{\link{im}}
#' @param label_im \code{\link{im}} 
#'  Label input image for whcih properties are extracted for each label.
#' @param intensity_im \code{\link{im}} 
#'  Intensity image with same dimensions as label image.
#' @param points \code{\link{ppp}} 
#'  Point pattern with a window matching the dimensions of the input label image.
#' @export
#' @return list
#'  List of region properties. The following properties can be accessed as methods
#'  for each region.
#'  \itemize{
#'    \item \emph{area}: Numeric
#'      Numeric value representing the area of the region.
#'    \item \emph{bbox}: \code{\link{owin}}
#'      Rectangular window representing the smallest rectangle that encloses the
#'      region.
#'    \item \emph{centroid}: list
#'      List containing x and y coordinates of the region.
#'    \item \emph{diameter}: Numeric
#'      Numeric value representing the size of the largest line that can be
#'      drawn across the region.
#'    \item \emph{perimeter}: Numeric
#'      Numeric value representing the perimeter around the region.
#'  }
#'  Note: the values are only computed when the method is called for a specific
#'  region i.e., they are lazy. This means you can filter on a specific property
#'  without needing to compute all properties.
region_props <- function(label_im, intensity_im = NULL, points = NULL) {
    t <- spatstat::tiles(spatstat::tess(image = label_im))

    # Creates a list of functions to each different properties for each
    # tile/region
    lapply(t, function(tile) {
      list(
        area = function() {
          spatstat::area.owin(tile)
        },
        bbox = function() {
          spatstat::boundingbox(tile)
        },
        centroid = function() {
          spatstat::centroid.owin(tile)
        },
        diameter = function() {
          spatstat::diameter.owin(tile)
        },
        perimeter = function() {
          spatstat::perimeter(tile)
        }
      )
    })
}

#' Calculate r value at L(r) - r peaks
#'
#' @seealso \code{\link{Lest}}
#' @param lest fv output from \code{\link{Lest}} function
#' @param variable names of the curve for which to calculate peak r
#' @return num peak r-value
#' @export
peak_r <- function(lest, variable) {
  df <- as.data.frame(lest)
  lminusr <- df[,variable] - df$r
  r_peak <- which.max(lminusr)
  df$r[r_peak]
}
