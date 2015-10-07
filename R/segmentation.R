#' Clear connected regions touching edges
#'
#' Clear connect regions of a binary image that touch either the edges of the
#' image or edges in specified binary mask (or window).
#'
#' @param image binary image (instance of \code{\link{im}} object)
#' @param win Optional.
#'      Window/mask that defines edges (instance of \code{\link{owin}} object)
#' @param r Optional.
#'      Radius of morphological erosion operation to use to define the edge. This
#'      will determine the 'thickness' of the edge. Defaults to 2px.
#' @export
#' @return binary image (\code{\link{im}} object) excluding any connected 
#'      components that touch the edge.
clear_border <- function(image, win = NULL, r = image$xstep*2) {
    label_im <- spatstat::connected(image)

    if (is.null(win))
        win <- owin(image$xrange, image$yrange, unitname=image$units)

    edge <- spatstat::setminus.owin(spatstat::as.mask(win, eps=label_im$xstep),
                                    spatstat::erosion.owin(spatstat::as.mask(win, eps=label_im$xstep), r = r))
    masked <- label_im[edge, drop=FALSE]
    label_im$v[label_im$v %in% levels(factor(masked$v))] <- NA
    spatstat::as.mask(label_im)
}

#' Remove connected regions/objects smaller than a specified size
#'
#' @param image binary image (instance of \code{\link{im}} object)
#' @param size Optional. Defaults to 50 units.
#'      Exclude any objects smaller than this size. Units default to units of
#'      the input \code{image}.
#' @export
#' @return binary image (\code{\link{im}} object) containing only connected 
#'      components > \code{size}.
remove_small_objects <- function(image, size = 50) {
    label_im <- connected(image)
    rp <- region_props(label_im)
    exclude <- names(Filter(function(r) { r$area() <= size }, rp))
    label_im$v[label_im$v %in% exclude] <- NA
    spatstat::as.mask(label_im)
}
