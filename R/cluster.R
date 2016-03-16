#' Summary statistics for all types in a Multitype Point Pattern
#'
#' Wrapper for \code{\link{alltypes}} that that allows for a
#' subset of the input \code{\link{ppp}} \code{marks} to be used,
#' as specified by the parameter \code{which.marks}.
#' @seealso \code{\link{alltypes}}
#' @param points multitype \code{\link{ppp}} object
#' @param which.marks string specifying column of the \code{\link{ppp}}
#'              \code{marks} data.frame to use as marks for cross-
#'              type analysis
#' @param ... argument passed to \code{\link{alltypes}}
#' @return Object of class "fv" - see \code{\link{fv.object}}
#' @export
all_types <- function(points, which.marks=NULL, ...) {
  if (!spatstat::is.marked(points))
    stop("Multitype point pattern is required.")

  if (spatstat::markformat(points) == "dataframe") {
    if (is.null(which.marks)) {
      stop("which.marks must be specified for a ppp with multiple marks.")
    } else if (!is.character(which.marks) || length(which.marks) != 1) {
      stop("which.marks should be a string specifying a single column of marks.")
    }
    points <- subset(points, select = which.marks)
  }

  spatstat::alltypes(points, ...)
}

#' Local L-function analysis
#'
#' Wrapper for \code{\link{local_l}} that allows local L-function to
#' be computed on each type of a multitype \code{\link{ppp}} object
#' where type is specidied by the \code{marks} column denoted by the parameter
#' \code{which.marks}.
#'
#' If an r-value is specified, \code{\link{local_l}} values can be appended as
#' marks on the original \code{\link{ppp}} object. Otherwise, a
#' \code{\link{fv}} object is returned.
#'
#' @param points multitype \code{\link{ppp}} object
#' @param which.marks string specifies the column of the \code{\link{ppp}}
#'              \code{marks} data.frame that denotes points types.
#' @param rvalue radius of to determine L-function at.
#' @param mark.original if rvalue is specified, local L-function values can be
#'    appended to the original \code{\link{ppp}} object as a new column
#'    in \code{marks}.
#' @param ... argument passed to \code{\link{local_l}}
#' @export
#' @return If \code{rvalue} and \code{mark.original} are specified,
#'    a vector of L(r) values is appended to the original \code{\link{ppp}}
#'    marks. If mark.original is not specified then the L(r) vector is simply returned.
#'    If \code{rvalue} is not specified, \code{\link{fv}} object is returned.
multi_local_l <- function(points, which.marks=spatstat::marks(points),
                            rvalue=NULL, mark.original=FALSE, ...) {
  if (mark.original && is.null(rvalue))
    stop("An r-value must be specified to mark the original PPP.")
  if (!spatstat::is.marked(points))
    stop("Multitype point pattern is required.")

  spoints <- split(points, f = which.marks)

  ll <- lapply(spoints,
                function(p, rvalue, ...) {
                  local_l(p, rvalue=rvalue, ...)
                }, rvalue, ...)
  ll_labels <- names(ll)

  if (mark.original) {
    out.points <- Reduce(spatstat::superimpose, mapply(function(p, ll) {
      spatstat::marks(p) <- cbind(spatstat::marks(p), ll)
      p
    }, spoints, ll, SIMPLIFY=FALSE))

    # Fix marks data.frame headings and re-factorise split column.
    marx <- spatstat::marks(out.points)
    if (spatstat::markformat(points) == "vector") {
      colnames(marx) <- c(spatstat::short.deparse(substitute(points)),
                          paste("L", rvalue, sep=""))
      marx[,1] <- factor(marx[,1], labels = ll_labels)
    } else {
      colnames(marx) <- c(colnames(marx)[-length(colnames(marx))],
                          paste("L", rvalue, sep=""))
      marx[, colnames(marx) %in% which.marks] <- factor(marx[, colnames(marx) %in% which.marks], labels = ll_labels)
    }
    spatstat::marks(out.points) <- marx

    return(out.points)
  } else {
    return(spatstat::listof(ll))
  }
}

#' Cross-type Local L-function analysis
#'
#' Wrapper for \code{\link{local_l_cross}} that allows cross-type
#' local L-function to be computed for specified multitype
#' \code{\link{ppp}} object using \code{marks} specified by the
#' parameter \code{which.marks}.
#'
#' If an r-value is specified, localL values can be appended as marks on the
#' original \code{\link{ppp}} object.
#'
#' @seealso \code{\link{local_k_cross}}
#'          \code{\link{local_l_cross}}
#' @param points multitype \code{\link{ppp}} object
#' @param which.marks string specifies the column of the \code{\link{ppp}}
#'              \code{marks} data.frame that denotes points types.
#' @param rvalue radius of to determine L-function at.
#' @param mark.original if rvalue is specified, local L-function values can be
#'    appended to the original \code{\link{ppp}} object as a new column
#'    in \code{marks}.
#' @param ... argument passed to \code{\link{local_l}}
#' @export
#' @return If \code{rvalue} and \code{mark.original} are specified,
#'    a vector of L(r)cross values is appended to the original \code{\link{ppp}}
#'    marks. If mark.original is not specified then the L(r)cross vector is simply returned.
#'    If \code{rvalue} is not specified, \code{\link{fv}} object is returned.
multi_local_lcross <- function(points, which.marks=NULL,
                                 rvalue=NULL, mark.original=FALSE, ...) {
  if (mark.original && is.null(rvalue))
    stop("An r-value must be specified to mark the original PPP.")
  if (!spatstat::is.marked(points))
    stop("Multitype point pattern is required.")

  if (spatstat::markformat(points) == "dataframe") {
    if (is.null(which.marks)) {
      stop("which.marks must be specified for a ppp with multiple marks.")
    } else if (!is.character(which.marks) || length(which.marks) != 1) {
      stop("which.marks should be a string specifying a single column of marks.")
    }
    points <- subset(points, select = which.marks)
  }

  if (mark.original && length(levels(spatstat::marks(points))) > 2) {
    stop(paste("Adding localLcross values back to original multitype ppp with",
      "greater than two mark types is not supported."))
  }

  llc <- local_l_cross(points, rvalue=rvalue, ...)
  if (mark.original) {
    # Fix marks data.frame headings and re-factorise split column.
    marx <- cbind(spatstat::marks(points), llc)
    if (spatstat::markformat(points) == "vector") {
      colnames(marx) <- c(spatstat::short.deparse(substitute(points)),
                          paste("L", rvalue, sep=""))
    } else {
      colnames(marx) <- c(colnames(marx)[-length(colnames(marx))],
                          paste("L", rvalue, sep=""))
    }
    spatstat::marks(points) <- marx

    return(points)
  } else {
    return(llc)
  }
}

#' Co-cluster analysis using self- and cross-type local L-function analysis
#'
#' Performs self- and cross-type local L-function analysis at a specified
#' radius on a multi-type PPP object. The column of marks to use for types
#' is specified using the \code{which.marks}. Users have the option to return
#' either a data.frame containing 2 columns with the self- and cross-type
#' local L-function values or the original PPP where the self- and cross-type
#' local L-function value are added to the marks.
#'
#' @param points A point pattern (object of class \code{ppp})
#' @param rvalue A single value of the distance argument r at which
#'  the L-function should be computed.
#' @param which.marks Specifies the column of marks to use as types. Must be
#'  specified is the PPP has multiple mark categories.
#' @param mark.original Logical value specifying whether the original point
#'  pattern should be marked with the output and returned.
#' @param ... Extra arguments are passed to the \code{local_l}
#'  and \code{local_l_cross} functions.
#' @return Self- and cross-type local L-functions values for each point are
#'  returned as either:
#'    \itemize{
#'      \item data.frame (if \code{mark.original} is \code{FALSE})
#'      \item marks in the original points (\code{ppp}) object passed as
#'            as argument (if \code{mark.original} is \code{TRUE})
#'    }
#' @export
co_cluster_l_cross <- function(points, rvalue, which.marks=NULL,
                      mark.original=FALSE, ...) {
  if (missing(rvalue))
    stop("An r-value must be specified.")
  if (!spatstat::is.marked(points))
    stop("Multitype point pattern is required.")

  if (spatstat::markformat(points) == "dataframe") {
    if (is.null(which.marks)) {
      stop("which.marks must be specified for a ppp with multiple marks.")
    } else if (!is.character(which.marks) || length(which.marks) != 1) {
      stop("which.marks should be a string specifying a single column of marks.")
    }
    ppp <- subset(points, select = which.marks)
  } else {
    ppp <- points
  }

  spoints <- split(points, f = spatstat::marks(ppp))
  mark.names <- names(spoints)
  if(mark.original && length(mark.names) > 2)
    stop("Marking an original point pattern with greater than 2 types is not
         supported. Please disable mark.original")

  ll <- mapply(function(p, n, rvalue, ...) {
                  ll <- local_l(p, rvalue=rvalue, ...)
                  llc <- local_l_cross(ppp, i = n,
                                       j = mark.names[mark.names != n],
                                       rvalue = rvalue, ...)
                  df <- cbind(ll, llc)
                  df
                }, spoints, mark.names, rvalue, ..., SIMPLIFY=FALSE)

  if (mark.original) {
    out.points <- Reduce(spatstat::superimpose, mapply(function(p, ll) {
      spatstat::marks(p) <- cbind(spatstat::marks(p), ll)
      p
    }, spoints, ll, SIMPLIFY=FALSE))

    # Fix marks data.frame headings and re-factorise split column.
    marx <- spatstat::marks(out.points)
    if (spatstat::markformat(points) == "vector") {
      colnames(marx) <- c(spatstat::short.deparse(substitute(points)),
                          paste("L", rvalue, sep=""),
                          paste("L", rvalue, "c", sep=""))
      marx[,1] <- factor(marx[,1], labels = mark.names)
    } else {
      colnames(marx) <- c(colnames(spatstat::marks(points)),
                          paste("L", rvalue, sep=""),
                          paste("L", rvalue, "c", sep=""))
      marx[, colnames(marx) %in% which.marks] <- factor(marx[, colnames(marx) %in% which.marks],
                                                        labels = mark.names)
    }
    spatstat::marks(out.points) <- marx

    return(out.points)
  } else {
    out <- Reduce(rbind, mapply(function(p, ll) {
      df <- cbind(spatstat::marks(p), ll)
      colnames(df) <- c("type",
                        paste("L", rvalue, sep=""),
                        paste("L", rvalue, "c", sep=""))
      df
    }, spoints, ll, SIMPLIFY=FALSE))

    out[,1] <- factor(out[,1], labels = mark.names)
    return(out)
  }
}


#' Quadrant analysis of self- vs. cross local L-function values
#'
#' @description Divides a relationship between self- and cross-type local
#' L-function values into quadrants based on a sepcified threshold. Determines
#' the quadrant number in which each point in a multi-type point pattern lies,
#' where quadrants are labelled:
#' \enumerate{
#'    \item lower left
#'    \item lower right
#'    \item upper left
#'    \item upper right
#'  }
#'
#' @param points point pattern (multi-type \code{ppp} object)
#' @param threshold threshold that splits x and y axis of the quadrant plot
#' @param ll Column name of the self-type local L values
#' @param llcross Column names of the cross-type local L values
#' @param mark.original Logical value indicating whether labels for each point should
#'    be added back to the original point pattern as marks or return as a vector.
#' @return vector or point pattern depending on \code{mark.original}
#' @export
co_cluster_quadrant <- function(points, threshold, ll, llcross, mark.original = TRUE) {
  quad <- function(ll, llcross, threshold) {
    if(ll > threshold && llcross > threshold) {
      return("4")
    } else if(ll > threshold && llcross <= threshold) {
      return("2")
    } else if(ll <= threshold && llcross > threshold) {
      return("3")
    } else if(ll <= threshold && llcross <= threshold) {
      return("1")
    }
  }

  spatstat::verifyclass(points, "ppp")
  if(!spatstat::is.marked(points)) stop("Point pattern must be marked.")
  if(spatstat::markformat(points) != "dataframe")
    stop("Point pattern marks must be a data.frame which includes columns for LocalL and LocalLCross values.")

  marx <- spatstat::marks(points)
  if(all(colnames(marx) %in% c(ll, llcross)))
    stop(paste(ll, " and/or ", llcross, "were not found in mark columns."))

  group <- mapply(quad, ll=marx[,ll], llcross=marx[,llcross],
                                  threshold=threshold)
  group <- factor(group)
  if (mark.original) {
    n <- paste("Quad", threshold, sep="")
    marx[,n] <- group
    spatstat::marks(points) <- marx
    return(points)
  } else {
    return(group)
  }
}

#' Interpolate a surface over a spatial point pattern
#'
#' Grid interpolation to create a surface map over a spatial point pattern
#' where the marks of the point pattern are interpolated between the points.
#' Function currently uses matlab \code{griddata} function for interpolation
#' using the \code{R.matlab} package.
#'
#' @param matlab \code{Matlab} instance
#' @param data marked point pattern (\code{\link{ppp}} object) with
#'    with a single vector of marks to be interpolated as a surface
#' @param xsize x pixel size of output image
#' @param ysize y pixel size of output image
#' @param method interpolation method for matlab gridata function, choices are
#'    "nearest", "linear", "natural", "cubic" and "v4" (default).
#' @export
#' @return \code{\link{im}} of interpolated surface with pixel sizes
#'    specified by \code{xsize} and \code{ysize}.
interpolate_surface <- function(matlab, data, xsize = 10, ysize = 10, method = "v4") {
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("R.matlab needed for this function to work. Please install it.",
         call. = FALSE)
  }

  stopifnot(R.matlab::isOpen(matlab))
  stopifnot(spatstat::is.ppp(data))
  if (!spatstat::is.marked(data)) {
    stop("Marked point pattern is required.")
  } else if (spatstat::markformat(data) != "vector") {
    stop("Marks must be a vector")
  }

  data.mat <- as.matrix(as.data.frame(data))
  win.mat <- as.matrix(data.frame(xrange = data$window$xrange,
                                  yrange = data$window$yrange))

  R.matlab::setVariable(matlab, data = data.mat, win = win.mat, xsize = xsize,
                        ysize = ysize, m = method)
  R.matlab::evaluate(matlab, "xlin=linspace(win(1,1), win(2,1), (win(2,1)-win(1,1))/xsize);",
                     "ylin=linspace(win(1,2), win(2,2), (win(2,2)-win(1,2))/ysize);",
                     "[xq, yq] = meshgrid(xlin, ylin);",
                     "vq = griddata(data(:,1), data(:,2), data(:,3), xq, yq, m);")

  out <- R.matlab::getVariable(matlab, "vq")
  spatstat::im(out$vq, xrange = data$window$xrange, yrange = data$window$yrange,
               unitname = spatstat::unitname(data))
}


