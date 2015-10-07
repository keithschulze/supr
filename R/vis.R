#' Fortify method for \code{\link{ppp}} objects
#'
#' This method extracts coordinate and marks data from
#' \code{\link{ppp}}objects to a data.frame, so that it can be
#' plotted with ggplot.
#'
#' @method fortify ppp
#' @seealso \code{\link{ppp}}
#' @param model ppp object
#' @param data not used by this method
#' @param ... ignored by this method
#' @return data.frame containing at least two columns: \code{x} and \code{y}
#'  which represent the \code{x} and \code{y} coordinates of each point in a
#'  \code{\link{ppp}} object. Returned data.frame can also contain
#'  multiple other columns depending of the type of marks associated with the
#'  \code{\link{ppp}} object (see \code{\link{ppp}}). If mark type is:
#'  \itemize{
#'    \item \emph{vector}: a single column with heading \code{marks} will be added.
#'    \item \emph{data.frame}: the marks data.frame with be column bound to the
#'      output data.frame, where the column headings from the marks data.frame
#'      are preserved.
#'  }
#' @export fortify.ppp
fortify.ppp <- function(model, data, ...) {
  spatstat::as.data.frame.ppp(model)
}

#' Fortify method for \code{\link{psp}} objects
#'
#' This method extracts coordinate and marks data from
#' \code{\link{psp}} objects to a data.frame, so that it can be
#' plotted with ggplot.
#'
#' @method fortify psp
#' @seealso \code{\link{psp}}
#' @param model psp object
#' @param data not used by this method
#' @param ... ignored by this method
#' @return data.frame containing at least four columns: \code{x0}, \code{x1},
#'  \code{y0} and \code{y1} with represent the coordinates of a line segement.
#'  Returned data.frame can also contain multiple other columns depending of
#'  the type of marks associated with the \code{\link{ppp}} object
#'  (see \code{\link{ppp}}). If mark type is:
#'  \itemize{
#'    \item \emph{vector}: a single column with heading \code{marks} will be added.
#'    \item \emph{data.frame}: the marks data.frame with be column bound to the
#'      output data.frame, where the column headings from the marks data.frame
#'      are preserved.
#'  }
#' @export fortify.psp
fortify.psp <- function(model, data, ...) {
  spatstat::as.data.frame.psp(model)
}

#' Fortify method for \code{\link{im}} objects
#'
#' This method extracts coordinate and intensity data from a
#' \code{\link{im}} to a data frame, so that it can be plotted
#' with ggplot.
#'
#' @method fortify im
#' @seealso \code{\link{im}}
#' @param model im object
#' @param data not used by this method
#' @param ... ignored by this method
#' @return data.frame with 3 columns: \code{x}, \code{y} and \code{value}
#'  where \code{x} & \code{y} are the pixel coordinates, \code{value} is the
#'  intensity value at that pixel coordinate.
#' @export fortify.im
fortify.im <- function(model, data, ...) {
  spatstat::as.data.frame.im(model)
}

#' Fortify method for unpacking \code{\link{owin}} objects.
#'
#' This method extracts relevant coordinates for plotting a
#' \link{\code{owin}} object to allow it to be plotted with ggplot.
#'
#' @method fortify owin
#' @seealso \code{\link{owin}}
#' @param model owin object
#' @param data not used by this method
#' @param ... ignored
#' @return data.frame containing windows information. Information in the
#'  data.frame depends on the type of the \code{\link{owin}}:
#'  \itemize{
#'    \item \emph{rectangular windows}: data.frame with a single row and
#'      4 columns that contain xmin, xmax, ymin, ymax extents of the recangle.
#'    \item \emph{polygonal window}: data.frame with three columns: \code{id},
#'      \code{x}, \code{y} where \code{id} represents a unique id for each
#'      polygon when there are multiple, and \code{x} & \code{y} represent the
#'      x and y coordinates of each vertex in the polygon.
#'    \item \emph{mask window}: data.frame with three columns: \code{x},
#'      \code{y} and \code{value} where \code{x} & \code{y} are pixel
#'      coordinates, and \code{value} is a logical value denoting whether the
#'      pixel is part of the window.
#'  }
#' @export fortify.owin
fortify.owin <- function(model, data, ...) {
  if (model$type == "rectangle") {
    xmin <- model$xrange[1]
    xmax <- model$xrange[2]
    ymin <- model$yrange[1]
    ymax <- model$yrange[2]
    data.frame(xmin, xmax, ymin, ymax)
  } else if (model$type == "polygonal") {
    df <- ldply(mapply(function(win, index) {
      data.frame(id=index, x=win$x, y=win$y)
    }, model$bdry, 1:length(model$bdry), SIMPLIFY = FALSE),
    data.frame)
    df$id <- factor(df$id)
    df
  } else if (model$type == "mask") {
    df <- spatstat::as.data.frame.owin(model, drop=FALSE)
    df$inside <- factor(df$inside)
    setNames(df, c("x", "y", "value"))
  } else {
    stop("Unrecognised window type")
  }
}
