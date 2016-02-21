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
clear_border <- function(image, win = NULL, r = image$xstep*2, output_mask = FALSE) {
    if (spatstat::is.mask(image)) image <- spatstat::connected(image)

    if (is.null(win))
        win <- spatstat::owin(image$xrange, image$yrange, unitname=image$units)

    edge <- spatstat::setminus.owin(spatstat::as.mask(win, eps=image$xstep),
                                    spatstat::erosion.owin(spatstat::as.mask(win, eps=image$xstep), r = r))
    masked <- image[edge, drop=FALSE]
    image$v[image$v %in% levels(factor(masked$v))] <- NA
    if (output_mask) {
      return(spatstat::as.mask(image))
    } else {
      return(image)
    }
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
remove_small_objects <- function(image, size = 50, output_mask = FALSE) {
    if (spatstat::is.mask(image)) image <- spatstat::connected(image)
    rp <- region_props(image)
    exclude <- names(Filter(function(r) { r$area() <= size }, rp))
    image$v[image$v %in% exclude] <- NA
    if (output_mask) {
      return(spatstat::as.mask(image))
    } else {
      return(image)
    }
}

#' Watershed segmentation of spatstat im object. Requires imager package
#' \link{\code{watershed}} function to perform watershed segmentation.
#'
#' @param im image object (grayscale or binary)
#' @param seeds planar point pattern containing seed/marker coordinates
#' @export
#' @return label image.
watershed <- function(mask, seeds = NULL, fill_lines = TRUE) {
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("imager package needed for this function to work. Please install it.",
      call. = FALSE)
  }

  cmask <- imager::im2cimg(spatstat::as.im(mask, na.replace = 0))

  # Create seed image
  seed_img <- spatstat::pixellate(seeds, W = mask, padzero = TRUE)
  seed_img$v <- as.matrix(as.data.frame.matrix(seed_img$v)) #Hack to get rid of table annot
  seed_cimg <- imager::im2cimg(seed_img)

  # Creat distance transform and seed labels
  seed_cimg_dm <- 1 - imager::distance_transform(sign(seed_cimg), 1)
  # seed_cimg_dm <- imager::distance_transform(cmask, value = 0)
  seed_cimg_dm <- imager::mult(list(seed_cimg_dm-min(seed_cimg_dm), cmask))
  seed_cimg_lbl <- imager::label(seed_cimg)

  # Do watershed
  seed_cimg_ws <- imager::watershed(seed_cimg_lbl, seed_cimg_dm, fill_lines = fill_lines)
  seed_cimg_ws <- imager::mult(list(seed_cimg_ws, cmask))

  # Convert back to im
  im_ws <- imager::cimg2im(seed_cimg_ws, W = spatstat::as.rectangle(mask))

  #Convert fill line to background
  im_ws$v[im_ws$v == 0] <- NA

  return(im_ws)
}
