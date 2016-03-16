#' Neighbourhood density function
#'
#' Computes the neighbourhood density function, a local version of
#' the \eqn{K}-function or \eqn{L}-function, defined by Getis and Franklin (1987).
#' Note: Equivalent to \code{\link{localK}} and \code{\link{localL}}
#' except that they can take advantage of parallel computation \eqn{K}-function
#' or \eqn{L}-function via \code{\link{foreach}} package.
#'
#' The command \code{localL} computes the \emph{neighbourhood density function},
#' a local version of the \eqn{L}-function (Besag's transformation of Ripley's
#' \eqn{K}-function) that was proposed by Getis and Franklin (1987).
#' The command \code{localK} computes the corresponding
#' local analogue of the K-function.
#'
#' Given a spatial point pattern \code{X}, the neighbourhood density function
#' \eqn{L_i(r)}{L[i](r)} associated with the \eqn{i}th point
#' in \code{X} is computed by
#' \deqn{
#'   L_i(r) = \sqrt{\frac a {(n-1) \pi} \sum_j e_{ij}}
#' }{
#'   L[i](r) = sqrt( (a/((n-1)* pi)) * sum[j] e[i,j])
#' }
#' where the sum is over all points \eqn{j \neq i}{j != i} that lie
#' within a distance \eqn{r} of the \eqn{i}th point,
#' \eqn{a} is the area of the observation window, \eqn{n} is the number
#' of points in \code{X}, and \eqn{e_{ij}}{e[i,j]} is an edge correction
#' term (as described in \code{\link{Kest}}).
#' The value of \eqn{L_i(r)}{L[i](r)} can also be interpreted as one
#' of the summands that contributes to the global estimate of the L
#' function.
#'
#' By default, the function \eqn{L_i(r)}{L[i](r)} or
#' \eqn{K_i(r)}{K[i](r)} is computed for a range of \eqn{r} values
#' for each point \eqn{i}. The results are stored as a function value
#' table (object of class \code{"fv"}) with a column of the table
#' containing the function estimates for each point of the pattern
#' \code{X}.
#'
#' Alternatively, if the argument \code{rvalue} is given, and it is a
#' single number, then the function will only be computed for this value
#' of \eqn{r}, and the results will be returned as a numeric vector,
#' with one entry of the vector for each point of the pattern \code{X}.
#'
#' Inhomogeneous counterparts of \code{local_k} and \code{local_l}
#' are computed by \code{local_k_inhom} and \code{local_l_inhom}.
#'
#' Computation can be done in parallel by registering a parallel backend for
#' the \code{\link{foreach}} package.
#'
#' @aliases local_l
#' @seealso \code{\link{Kest}}
#'          \code{\link{Lest}}
#'          \code{\link{localKinhom}}
#'          \code{\link{localLinhom}}
#'          \code{\link{local_k_inhom}}
#'          \code{\link{local_l_inhom}}
#' @param X A point pattern (object of class \code{"ppp"}).
#' @param ... ignored.
#' @param correction tring specifying the edge correction to be applied.
#'  Options are \code{"none"}, \code{"translate"}, \code{"translation"},
#'  \code{"Ripley"},
#'  \code{"isotropic"} or \code{"best"}.
#'  Only one correction may be specified.
#' @param verbose Logical flag indicating whether to print progress
#'  reports during the calculation.
#' @param rvalue Optional. A \emph{single} value of the distance argument
#'  \eqn{r} at which the function L or K should be computed.
#' @export
#' @return If \code{rvalue} is given, the result is a numeric vector
#'  of length equal to the number of points in the point pattern.
#'
#'  If \code{rvalue} is absent, the result is
#'  an object of class \code{"fv"}, see \code{\link{fv.object}},
#'  which can be plotted directly using \code{\link{plot.fv}}.
#'  Essentially a data frame containing columns
#'  \item{r}{the vector of values of the argument \eqn{r}
#'    at which the function \eqn{K} has been  estimated
#'  }
#'  \item{theo}{the theoretical value \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}
#'   or \eqn{L(r)=r} for a stationary Poisson process
#'  }
#'  together with columns containing the values of the
#'  neighbourhood density function for each point in the pattern.
#'  Column \code{i} corresponds to the \code{i}th point.
#'  The last two columns contain the \code{r} and \code{theo} values.
#' @references Getis, A. and Franklin, J. (1987)
#'      Second-order neighbourhood analysis of mapped point patterns.
#'      \emph{Ecology} \bold{68}, 473--477.
#' @examples
#'  X <- spatstat::ponderosa
#'
#'  # compute all the local L functions
#'  L <- local_l(X)
#'
#'  # All local functions can also be executed in parallel. Simply register a
#'  # parallel backend for the foreach package. For example using the
#'  # doParallel (this needs to be installed) backend:
#'
#'  cl <- parallel::makeCluster(2)
#'  doParallel::registerDoParallel(cl)
#'  L <- local_l(X)
#'  foreach::registerDoSEQ()
#'  parallel::stopCluster(cl)
#'
#'  # plot all the local L functions against r
#'  plot(L, main="local L functions for ponderosa", legend=FALSE)
local_k <- function(X, ..., correction="Ripley", verbose=TRUE, rvalue=NULL) {
  spatstat::verifyclass(X, "ppp")
  local_k_engine(X, ..., correction=correction, verbose=verbose, rvalue=rvalue)
}

#' @export
local_l <- function(X, ..., correction="Ripley", verbose=TRUE, rvalue=NULL) {
  local_k(X, wantL=TRUE, correction=correction, verbose=verbose, rvalue=rvalue)
}

#' Inhomogeneous Neighbourhood density function
#'
#' Computes spatially-weighted versions of the local \eqn{K}-function or \eqn{L}-function.
#' Note: Equivalent to \code{\link{localKinhom}} and \code{\link{localLinhom}}
#' except that they can take advantage of parallel computation \eqn{K}-function
#' or \eqn{L}-function via \code{\link{foreach}} package.
#'
#' The functions \code{local_k_inhom} and \code{local_l_inhom}
#' are inhomogeneous or weighted versions of the
#' neighbourhood density function implemented in
#' \code{\link{local_k}} and \code{\link{local_l}}.
#'
#' Given a spatial point pattern \code{X}, the inhomogeneous neighbourhood
#' density function \eqn{L_i(r)}{L[i](r)} associated with the \eqn{i}th point
#' in \code{X} is computed by
#' \deqn{
#'   L_i(r) = \sqrt{\frac 1 \pi \sum_j \frac{e_{ij}}{\lambda_j}}
#' }{
#'   L[i](r) = sqrt( (1/pi) * sum[j] e[i,j]/lambda[j])
#' }
#' where the sum is over all points \eqn{j \neq i}{j != i} that lie
#' within a distance \eqn{r} of the \eqn{i}th point,
#' \eqn{\lambda_j}{\lambda[j]} is the estimated intensity of the
#' point pattern at the point \eqn{j}, and \eqn{e_{ij}}{e[i,j]} is an edge correction
#' term (as described in \code{\link{Kest}}).
#' The value of \eqn{L_i(r)}{L[i](r)} can also be interpreted as one
#' of the summands that contributes to the global estimate of the inhomogeneous L
#' function (see \code{\link{Linhom}}).
#'
#' By default, the function \eqn{L_i(r)}{L[i](r)} or
#' \eqn{K_i(r)}{K[i](r)} is computed for a range of \eqn{r} values
#' for each point \eqn{i}. The results are stored as a function value
#' table (object of class \code{"fv"}) with a column of the table
#' containing the function estimates for each point of the pattern
#' \code{X}.
#'
#' Alternatively, if the argument \code{rvalue} is given, and it is a
#' single number, then the function will only be computed for this value
#' of \eqn{r}, and the results will be returned as a numeric vector,
#' with one entry of the vector for each point of the pattern \code{X}.
#'
#' Computation can be done in parallel by registering a parallel backend for
#' the \code{\link{foreach}} package.
#'
#' @aliases local_l_inhom
#' @seealso \code{\link{Kinhom}}
#'          \code{\link{Linhom}}
#'          \code{\link{localKinhom}}
#'          \code{\link{localLinhom}}
#'          \code{\link{local_k}}
#'          \code{\link{local_l}}
#' @param X A point pattern (object of class \code{"ppp"}).
#' @param lambda Optional.
#'  Values of the estimated intensity function.
#'  Either a vector giving the intensity values
#'  at the points of the pattern \code{X},
#'  a pixel image (object of class \code{"im"}) giving the
#'  intensity values at all locations, a fitted point process model
#'  (object of class \code{"ppm"}) or a \code{function(x,y)} which
#'  can be evaluated to give the intensity value at any location.
#' @param ... ignored.
#' @param correction tring specifying the edge correction to be applied.
#'  Options are \code{"none"}, \code{"translate"}, \code{"translation"},
#'  \code{"Ripley"},
#'  \code{"isotropic"} or \code{"best"}.
#'  Only one correction may be specified.
#' @param verbose Logical flag indicating whether to print progress
#'  reports during the calculation.
#' @param rvalue Optional. A \emph{single} value of the distance argument
#'  \eqn{r} at which the function L or K should be computed.
#' @param sigma
#'  Optional arguments passed to \code{\link{density.ppp}} to control
#'  the kernel smoothing procedure for estimating \code{lambda},
#'  if \code{lambda} is missing.
#'  @param varcov
#'  Optional arguments passed to \code{\link{density.ppp}} to control
#'  the kernel smoothing procedure for estimating \code{lambda},
#'  if \code{lambda} is missing.
#' @export
#' @return If \code{rvalue} is given, the result is a numeric vector
#'  of length equal to the number of points in the point pattern.
#'
#'  If \code{rvalue} is absent, the result is
#'  an object of class \code{"fv"}, see \code{\link{fv.object}},
#'  which can be plotted directly using \code{\link{plot.fv}}.
#'  Essentially a data frame containing columns
#'  \item{r}{the vector of values of the argument \eqn{r}
#'    at which the function \eqn{K} has been  estimated
#'  }
#'  \item{theo}{the theoretical value \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}
#'   or \eqn{L(r)=r} for a stationary Poisson process
#'  }
#'  together with columns containing the values of the
#'  neighbourhood density function for each point in the pattern.
#'  Column \code{i} corresponds to the \code{i}th point.
#'  The last two columns contain the \code{r} and \code{theo} values.
#' @references Getis, A. and Franklin, J. (1987)
#'      Second-order neighbourhood analysis of mapped point patterns.
#'      \emph{Ecology} \bold{68}, 473--477.
local_k_inhom <- function(X, lambda=NULL, ..., correction="Ripley", verbose=TRUE,
                        rvalue=NULL, sigma=NULL, varcov=NULL) {
  spatstat::verifyclass(X, "ppp")

  if(is.null(lambda)) {
    # No intensity data provided
    # Estimate density by leave-one-out kernel smoothing
    lambda <- density(X, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=TRUE)
    lambda <- as.numeric(lambda)
  } else {
    # validate
    if(spatstat::is.im(lambda))
      lambda <- spatstat::safelookup(lambda, X)
    else if(spatstat::is.ppm(lambda))
      lambda <- predict(lambda, locations=X, type="trend")
    else if(is.function(lambda))
      lambda <- lambda(X$x, X$y)
    else if(is.numeric(lambda) && is.vector(as.numeric(lambda)))
      spatstat::check.nvector(lambda, spatstat::npoints(X))
    else stop(paste(sQuote("lambda"),
                    "should be a vector, a pixel image, or a function"))
  }
  local_k_engine(X, lambda=lambda, ...,
               correction=correction, verbose=verbose, rvalue=rvalue)
}

#' @export
local_l_inhom <- function(X, lambda=NULL, ..., correction="Ripley", verbose=TRUE,
                        rvalue=NULL, sigma=NULL, varcov=NULL) {
  local_k_inhom(X, lambda=lambda, wantL=TRUE, ...,
              correction=correction, verbose=verbose, rvalue=rvalue,
              sigma=sigma, varcov=varcov)
}

local_k_engine <- function(X, ..., wantL=FALSE, lambda=NULL,
                         correction="Ripley", verbose=TRUE, rvalue=NULL) {
  npts <- spatstat::npoints(X)
  W <- X$window
  areaW <- spatstat::area(W)
  lambda.ave <- npts/areaW
  lambda1.ave <- (npts - 1)/areaW

  weighted <- !is.null(lambda)

  if(is.null(rvalue))
    rmaxdefault <- spatstat::rmax.rule("K", W, lambda.ave)
  else {
    stopifnot(is.numeric(rvalue))
    stopifnot(length(rvalue) == 1)
    stopifnot(rvalue >= 0)
    rmaxdefault <- rvalue
  }
  breaks <- spatstat::handle.r.b.args(NULL, NULL, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max

  correction.given <- !missing(correction)
  correction <- spatstat::pickoption("correction", correction,
                                     c(none="none",
                                       isotropic="isotropic",
                                       Ripley="isotropic",
                                       trans="translate",
                                       translate="translate",
                                       translation="translate",
                                       best="best"),
                                     multi=FALSE)

  correction <- spatstat::implemented.for.K(correction, W$type, correction.given)

  # recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))

  # identify all close pairs
  rmax <- max(r)
  close <- spatstat::closepairs(X, rmax)
  DIJ <- close$d
  XI <- spatstat::ppp(close$xi, close$yi, window=W, check=FALSE)
  I <- close$i
  if(weighted) {
    J <- close$j
    lambdaJ <- lambda[J]
    weightJ <- 1/lambdaJ
  }

  iseq <- 1:npts
  icode <- sprintf("%02d", iseq)
  bkt <- function(x) { paste("[", x, "]", sep="") }

  parallel <- check_parallel(verbose)
  fe <- foreach::foreach(i = iseq, .export = c("local_k_calc"))

  switch(correction,
         none={
           # uncorrected! For demonstration purposes only!
           df <- foreach::`%dopar%`(fe, local_k_calc(i, closeI=I, distIJ=DIJ, breaks=breaks,
                                                   weight=if (weighted) weightJ else NULL,
                                                   pb = if(verbose && !parallel) function(i) spatstat::progressreport(i, n=npts) else NULL))

           # Hack to quickly convert list to data.frame
           attributes(df) <- list(row.names=c(NA_integer_, length(r)),
                                  class="data.frame",
                                  names=make.names(paste("un", icode, sep=""),
                                                   unique=TRUE))
           labl <- paste("%s", bkt(icode), "(r)", sep="")
           desc <- paste("uncorrected estimate of %s",
                         "for point", icode)

           if (!weighted) df <- df / lambda1.ave
         },
         translate={
           # Translation correction
           XJ <- spatstat::ppp(close$xj, close$yj, window=W, check=FALSE)
           edgewt <- spatstat::edge.Trans(XI, XJ, paired=TRUE)
           if (weighted)
             edgewt <- edgewt * weightJ

           df <- foreach::`%dopar%`(fe, local_k_calc(i, closeI=I, distIJ=DIJ,
                                                   breaks=breaks, weight=edgewt,
                                                   pb = if(verbose && !parallel) function(i) spatstat::progressreport(i, n=npts) else NULL))

           # Hack to quickly convert list to data.frame
           attributes(df) <- list(row.names=c(NA_integer_, length(r)),
                                  class="data.frame",
                                  names=make.names(paste("trans", icode, sep=""),
                                                   unique=TRUE))
           labl <- paste("%s", bkt(icode), "(r)", sep="")
           desc <- paste("translation-corrected estimate of %s",
                         "for point", icode)

           if(!weighted) df <- df/lambda1.ave
           h <- spatstat::diameter(W)/2
           df[r >= h, ] <- NA
         },
         isotropic={
           # Ripley isotropic correction
           edgewt <- spatstat::edge.Ripley(XI, matrix(DIJ, ncol=1))
           if(weighted)
             edgewt <- edgewt * weightJ

           df <- foreach::`%dopar%`(fe, local_k_calc(i, closeI=I, distIJ=DIJ,
                                                   breaks=breaks, weight=edgewt,
                                                   pb = if(verbose && !parallel) function(i) spatstat::progressreport(i, n=npts) else NULL))

           # Hack to quickly convert list to data.frame
           attributes(df) <- list(row.names=c(NA_integer_, length(r)),
                                  class="data.frame",
                                  names=make.names(paste("iso", icode, sep=""),
                                                   unique=TRUE))
           labl <- paste("%s", bkt(icode), "(r)", sep="")
           desc <- paste("Ripley isotropic correction estimate of %s",
                         "for point", icode)

           if(!weighted) df <- df/lambda1.ave
           h <- spatstat::diameter(W)/2
           df[r >= h, ] <- NA
         })

  # transform values if L required
  if (wantL)
    df <- sqrt(df / pi)

  # return vector of values at r=rvalue, if desired
  if(!is.null(rvalue)) {
    nr <- length(r)
    if(r[nr] != rvalue)
      stop("Internal error - rvalue not attained")
    return(as.numeric(df[nr,]))
  }

  # function value table required
  # add r and theo
  if(!wantL) {
    df <- cbind(df, data.frame(r=r, theo=pi * r^2))
    if(!weighted) {
      ylab <- quote(K[loc](r))
      fnam <- "K[loc][',']"
    } else {
      ylab <- quote(Kinhom[loc](r))
      fnam <- "Kinhom[loc][',']"
    }
  } else {
    df <- cbind(df, data.frame(r=r, theo=r))
    if(!weighted) {
      ylab <- quote(L[loc](r))
      fnam <- "L[loc][',']"
    } else {
      ylab <- quote(Linhom[loc](r))
      fnam <- "Linhom[loc][',']"
    }
  }
  desc <- c(desc, c("distance argument r", "theoretical Poisson %s"))
  labl <- c(labl, c("r", "%s[pois](r)"))
  # create fv object
  K <- spatstat::fv(df, "r", ylab, "theo", "", alim, labl, desc, fname=fnam)
  # default is to display them all
  spatstat::formula(K) <- (. ~ r)
  spatstat::unitname(K) <- spatstat::unitname(X)
  attr(K, "correction") <- correction
  return(K)
}

# Determine Ripley's K value for each point
local_k_calc <- function(i, closeI, distIJ, breaks, weight=NULL, pb=NULL) {
  if (!is.null(pb))
    pb(i)
  ii <- (closeI == i)
  wh <- spatstat::whist(distIJ[ii], breaks$val,
                        if(!is.null(weight)) weight[ii] else NULL)
  cs <- cumsum(wh)

  return(cs)
}
