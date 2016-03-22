#' Feature detection functions.

#' Peak local max
#'
#' Determine local maxima of marked point pattern. Indentifies
#' local maxima of marks and returns the point pattern with only
#' those marks.
#'
#' @param points Marked planar point pattern (ppp)
#' @param min_distance Min distance around a local max in which other maxima
#'  cannot occur
#' @param threshold Minimal value that maxima can have (default = 50).
#' @param which.marks If ppp is marked with multiple types, this specifies
#'  which marks to calculated peaks from.
#' @export
#' @return ppp containing only the detected maxima.
peak_local_max <- function(points, min_distance, threshold = 50, which.marks = NULL) {

  # Determine local max from a subset of a planar point pattern
  local_max <- function(df, mask) {
    if (!"index" %in% colnames(df)) df$index <- rownames(df)
    if (!"mask" %in% colnames(df)) df$mask <- TRUE

    masked_df <- df[mask, ]
    masked_df_min_idx <- masked_df[-which.max(masked_df$marks),]$index
    df[masked_df_min_idx,'mask'] <- FALSE
    return(df)
  }

  if (!spatstat::is.ppp(points)) stop("points must be a planar point pattern (ppp.object).")

  if (!spatstat::is.marked(points)) stop("Marked point pattern is required.")

  if (spatstat::markformat(points) == "dataframe") {
    if (is.null(which.marks)) {
      stop("which.marks must be specified for a point pattern with multiple marks.")
    } else if (!is.character(which.marks) || length(which.marks) != 1) {
      stop("which.marks should be a string specifying a single column of marks.")
    }
    points <- subset(points, select = which.marks)
  }

  marx <- spatstat::marks(points)
  if (!is.numeric(marx)) stop("point pattern marks must be numeric")

  # exclude any points below threshold
  points <- subset(points, marks >= threshold)
  n <- spatstat::npoints(points)
  if (n == 0) return(ppp(c(), c(), window = spatstat::as.rectangle(points$window)))

  # calculate pairwise distance matrix between remaining marks
  dm <- spatstat::pairdist(points)

  # create a mask where only pairs within the minimum distance are TRUE
  dm_mask <- dm < min_distance

  # convert point pattern to a data.frame
  df <- spatstat::as.data.frame.ppp(points)

  # Create a mask of local maxima only using the local_max function. For each
  # column in distance matrix mask i.e. each point, determine the local maxima
  # of all points within the minimum distance e.g. those that are visible
  # through the mask.
  res <- Reduce(f=local_max, x=split(dm_mask, col(dm_mask)), init=df)
  points <- points[res$mask]

  points
}

#' Calculate r value at which L(r) - r peaks
#'
#' @seealso \code{\link{Lest}}
#' @param lest fv output from \code{\link{Lest}} function
#' @param variable names of the curve for which to calculate peak r
#' @return num peak r-value
#' @export
peak_r <- function(lest, variable, span=25) {
  peaks <- function(series, span) {
    z <- embed(series, span)
    s <- span%/%2
    v<- max.col(z) == 1 + s
    result <- c(rep(FALSE,s),v)
    result <- result[1:(length(result)-s)]
    result
  }

  r <- lest[, "r", drop=TRUE]
  obs <- lest[, variable, drop=TRUE]

  which(peaks(obs-r, span))[1]
  # df <- as.data.frame(lest)
  # lminusr <- df[,variable] - df$r
  # r_peak <- which.max(lminusr)
  # df$r[r_peak]
}
