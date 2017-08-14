#' Thunderstorm molecule list reader
#'
#' Reader function for comma separated molecule list output of the
#' Thunderstorm plugin for ImageJ/Fiji.
#'
#' @param filepath string denoting the path to the file to be read.
#' @return \code{\link{ppp}} object containing the coordinates for
#'          single molecule localisations and other parameters/columns in
#'          original csv file attached as marks.
#' @export
read.moleculelist_thunderstorm <-
    function(filepath) {
  data <- read.csv(filepath, header=TRUE)
  names(data) <- c("frame", "x", "y", "sigma", "intensity", "offset", "bkgstd", "uncertainty")
  fd_ppp <- spatstat::ppp(data$x, data$y,
    xrange=c(min(data$x), max(data$x)),
    yrange=c(min(data$y), max(data$y)),
    marks = data[, !colnames(data) %in% c("x", "y")],
    unitname="nm")
  return(fd_ppp)
}


#' Rapidstorm molecule list reader
#'
#' Reader function for comma separated molecule list output of the
#' RapidStorm.
#'
#' @param filepath string denoting the path to the file to be read.
#' @return \code{\link{ppp}} object containing the coordinates for
#'          single molecule localisations and other parameters/columns in
#'          original csv file attached as marks.
#' @export
read.moleculelist_rapidstorm <-
    function(filepath, hdrs = c("x", "y", "frame", "amplitude", "chisq", "bkgd")) {
  data <- read.table(file = filepath, header = FALSE)
  names(data) <- hdrs
  fd_ppp <- spatstat::ppp(data$x, data$y,
    xrange=c(min(data$x), max(data$x)),
    yrange=c(min(data$y), max(data$y)),
    marks = data[, !colnames(data) %in% c("x", "y")],
    unitname="nm")
  return(fd_ppp)
}

#' Nikon N-STORM molecule list reader
#'
#' Reader function for comma separated molecule list output of the
#' NIS N-STORM module.
#'
#' @param filepath string denoting the path to the file to be read.
#' @param nonspecific include 'Non Specific Activation' data (default = TRUE)
#' @return \code{\link{ppp}} object containing the coordinates for
#'          single molecule localisations and other parameters/columns in
#'          original csv file attached as marks.
#' @export
read.moleculelist_nstorm <- function(filepath, nonspecific = TRUE) {
  data <- read.delim(filepath, sep="\t")
  names(data) <- tolower(names(data))
  if(!nonspecific) {
    data <- data[!(data$channel.name %in% c("Non Specific Activation")), ]
  }
  fd_ppp <- spatstat::ppp(data$x, data$y,
    xrange = c(min(data$x), max(data$x)),
    yrange = c(min(data$y), max(data$y)),
    marks = data[, !colnames(data) %in% c("x", "y")],
    unitname="nm")
  return(fd_ppp)
}

#' Function to read multiple datasets using a specific reader
#'
#' @param file_paths list of file_paths to read
#' @param rdr reader function that takes a file path as input and returns
#'            a ppp object.
#' @param ... extra arguments passed to \code{rdr} function.
#' @return \code{list(\link{ppp})}
#' @export
read.moleculelists <- function(file_paths, rdr = read.moleculelist_rapidstorm, ...) {
  return(lapply(file_paths, rdr, ...))
}

#' Read a polygon owin from a tab-delimited txt file
#'
#' Reads x and y coordinates from a tab delimited files where each XY coord is
#' on a separate row and X and Y coords are in columns 1 and 2, respectively.
#'
#' @param file_path path to txt file
#' @param xsize x pixel size
#' @param ysize y pixel size
#' @return \code{\link{owin}}
#' @export
read.polygon_roi <- function(file_path, xsize = 1, ysize = 1) {
  win <- read.table(file_path, header=FALSE)

  asxy <- function(xy) {
    list(x = xy[, 1]*xsize, y = xy[, 2]*ysize)
  }
  win <- asxy(win)

  w.area <- spatstat.utils::Area.xypolygon(win)
  if (sum(w.area) < 0) win <- spatstat.utils::reverse.xypolygon(win)

  return(spatstat::owin(poly=win, unitname="nm"))
}

#' Function to read multiple windows using a specified reader function
#'
#' @param file_paths list of file_paths to read
#' @param rdr reader function that takes a file path as input and returns
#'            an owin object.
#' @param ... passed to reader function
#' @return \code{list(\link{owin})}
#' @export
read.rois <- function(file_paths, rdr = read.polygon_roi, ...) {
  return(lapply(file_paths, rdr, ...))
}
