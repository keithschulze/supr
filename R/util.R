#' Exchange \code{x} and \code{y} coords for others stored as marks
#'
#' Utility function that allows the \code{x} and/or \code{y} coords
#' of a PPP object to be swapped with another set of coords stored as marks.
#' The previous \code{x} and \code{y} coords are stored in the marks under
#' names \code{orig_x} and \code{orig_y} respectively.
#'
#' @param points ppp object containing marks
#' @param x_name name of the marks column to be used as the new \code{x} coords
#' @param y_name name of the marks column to be used as the new \code{y} coords
#' @export
#' @return \code{\link{ppp}} where x and y coords have been exchanged for those
#'   those specified as marks.
exchange_xy_for_marks <- function(points, x_name = "x", y_name = "y") {
  df <- as.data.frame(points)
  x <- df[, x_name]
  y <- df[, y_name]
  marks <- df[, !colnames(df) %in% c(x_name, y_name)]
  names(marks)[names(marks) == "x"] <- "orig_x"
  names(marks)[names(marks) == "y"] <- "orig_y"
  return(spatstat::ppp(x, y, xrange = c(min(x), max(x)),
    yrange = c(min(y), max(y)), marks = marks,
    units = spatstat::unitname(points)))
}

#' Simple function to establish whether things can be run in parallel using
#' \code{\link{foreach}}
#'
#' Checks whether a parallel backend is setup for the \code{\link{foreach}}
#' package.
#'
#' @param verbose specify whether to show output messages.
#' @return boolean specifying whether a parallel backend has been setup.
check_parallel <- function(verbose = TRUE) {
  if (foreach::getDoParWorkers() == 1) {
    if (verbose)
      message("No parallel backend registered. Execution will be sequential.")
    return(FALSE)
  } else {
    if (verbose)
      message("Progress is disabled when running in parallel.")
    return(TRUE)
  }
}

#' Utility function to start matlab server using R.matlab.
#'
#' Simple wrapper function for starting the matlab server using R.matlab
#' @param matlab_path path to the matlab executable
#' @param port port to start matlab server on
#' @param remote If TRUE, all data to and from the MATLAB server will be 
#'    transferred through the socket connection, otherwise the data will be 
#'    transferred via a temporary file.
#' @param interval interval at which to poll server to check for results.
#' @param maxTries maximum number of times to poll server for results.
#' @export
#' @return \code{\link{Matlab}} instance
start_matlab <- function(matlab_path=NULL, port=9999, remote=FALSE, interval=NULL,
                         maxTries=NULL) {
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("R.matlab needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!is.null(matlab_path))
    options(matlab=matlab_path)

  R.matlab::Matlab$startServer(port = port)
  matlab <- R.matlab::Matlab(remote = remote, port = port)
  if (!is.null(interval))
    R.matlab::setOption(matlab, "readResult/interval", interval)
  if (!is.null(maxTries))
    R.matlab::setOption(matlab, "readResult/maxTries", maxTries)

  return(matlab)
}

#' Run matlab code using R.matlab 
#'
#' Utility that starts a matlab server, initialises an R.matlab 
#' \code{\link{Matlab}} instance connected to the Matlab server and passes it to 
#' a single arity input function. Therefore code inside the input function has
#' access to the \code{\link{Matlab}} instance. Once the function has completed 
#' execution or if a failure occurs, the \code{\link{Matlab}} instance is 
#' disconnected and the Matlab server shutdown. 
#'
#' @param fn Single arity function that accepts an initialised/connected 
#'  \code{\link{Matlab}} object as input.
#' @param ... passed as arguments to \code{\link{start_matlab}}
#' @export
#' @return results of the input \code{fn} function
with_matlab <- function(fn, ...) {
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("R.matlab needed for this function to work. Please install it.",
         call. = FALSE)
  }

  mlab <- start_matlab(...)

  R.matlab::open.Matlab(mlab)
  result <- tryCatch(fn(mlab),
                     finally = { R.matlab::close.Matlab(mlab) })

  return(result)
}
