#' Multitype (cross-type) neighbourhood density function
#'
#' Computes the cross-type neighbourhood density function, a local version of
#' the \eqn{K}-function or \eqn{L}-function, defined by Getis and Franklin (1987).
#'
#' The command \code{local_l_cross} computes the \emph{neighbourhood density function},
#' a local version of the \eqn{L}-function (Besag's transformation of Ripley's
#' \eqn{K}-function) proposed by Getis and Franklin (1987), for 2 types in a multitype
#' spatial point pattern. The command \code{local_k_cross} computes the corresponding
#' local analogue of the cross-type K-function.
#'
#' Given a multitype spatial point pattern \code{X} with types \code{i} and \code{j},
#' the neighbourhood density function
#' \eqn{L_{ij}(r)}{L[ij](r)} associated with the \eqn{i}th point
#' in \code{X} is computed by
#' \deqn{
#'   L_{ij}(r) = \sqrt{\frac a {(n_j) \pi} \sum_j e_{ij}}
#' }{
#'   L[ij](r) = sqrt( (a/((n[j])* pi)) * sum[j] e[i,j])
#' }
#' where the sum is over all points \eqn{j}{j} that lie
#' within a distance \eqn{r} of the \eqn{i}th point,
#' \eqn{a} is the area of the observation window, \eqn{n_j}{n[j]} is the number
#' of \code{j} points in \code{X}, and \eqn{e_{ij}}{e[i,j]} is an edge correction
#' term (as described in \code{\link{Kest}}).
#' The value of \eqn{L_{ij}(r)}{L[ij](r)} can also be interpreted as one
#' of the summands that contributes to the global estimate of the cross-type
#' \eqn{L}-function.
#'
#' By default, the function \eqn{L_{ij}(r)}{L[ij](r)} or
#' \eqn{K_{ij}(r)}{K[ij](r)} is computed for a range of \eqn{r} values
#' for each point \eqn{i}. The results are stored as a function value
#' table (object of class \code{"fv"}) with a column of the table
#' containing the function estimates for each point \eqn{i} of the pattern
#' \code{X}.
#'
#' Alternatively, if the argument \code{rvalue} is given, and it is a
#' single number, then the function will only be computed for this value
#' of \eqn{r}, and the results will be returned as a numeric vector,
#' with one entry of the vector for each point \eqn{i} of the pattern \code{X}.
#'
#' Inhomogeneous counterparts of \code{local_k_cross} and \code{local_l_cross}
#' are computed by \code{local_k_cross_inhom} and \code{local_l_cross_inhom}.
#'
#' Computation can be done in parallel by registering a parallel backend for
#' the \code{\link{foreach}} package.
#'
#' @aliases local_l_cross
#' @seealso \code{\link{localL}}
#'          \code{\link{localK}}
#'          \code{\link{Lcross}}
#'          \code{\link{Kcross}}
#' @param X Multitype point pattern (\code{\link{ppp}} object). It must be a multitype
#'      point pattern (or marked point pattern).
#' @param i The type (mark value) of the points in X from which distances are
#'      measured. A character string (or something that will be converted to a
#'      character string). Defaults to the first level of \code{\link{marks}(X)}.
#' @param j The type (marks value) of the points in X to which distances are
#'      measured. A character string (or something that will be converted to a
#'      character string'). Defaults to the second level of \code{\link{marks}(X)}.
#' @param ... Ignored.
#' @param correction String sprecifying the edge correction to be applied.
#'      Option are \code{"none"}, \code{"translate"}, \code{"translation"},
#'      \code{"Ripley"}, \code{"isotropic"} or \code{"best"}. Only one correction
#'      may be given.
#' @param verbose Logical flag indicating whether to print progress reports
#'      during the calculation.
#' @param rvalue Optional. A \emph{single} value of the distance argument \eqn{r}
#'      at which the function L or K should be computed.
#' @return If \code{rvalue} is given, the result is a numeric vector of equal
#'      length to the number of points in X_i.
#'
#'      If the code{rvalue} is absent, the result is an object of class \code{"fv"},
#'      see \code{\link{fv.object}}, which can be plotted directly using
#'      \code{\link{plot.fv}}. Essentially a data frame containing columns:
#'      \item{r}{the vector of values of the argument \eqn{r}
#'      at which the function \eqn{K} has been  estimated}
#'      \item{theo}{the theoretical value \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}
#'      or \eqn{L(r)=r} for a stationary Poisson process
#'      }
#'      together with columns containing the values of the
#'      neighbourhood density function for each point in the pattern.
#'      Column \code{i} corresponds to the \code{i}th point.
#'      The last two columns contain the \code{r} and \code{theo} values.
#' @references Getis, A. and Franklin, J. (1987)
#'      Second-order neighbourhood analysis of mapped point patterns.
#'      \emph{Ecology} \bold{68}, 473--477.
local_k_cross <- function(X, i, j, ..., correction="Ripley", verbose=TRUE, rvalue=NULL) {
  spatstat::verifyclass(X, "ppp")
  if(!spatstat::is.multitype(X, dfok=FALSE))
    stop("Point pattern must be multitype")
  marx <- spatstat::marks(X)
  if(missing(i))
    i <- levels(marx)[1]
  if(missing(j))
    j <- levels(marx)[2]

  I <- (marx == i)
  if(!any(I))
    stop(paste("No points have mark i =", i))

  J <- (marx == j)
  if(!any(J))
    stop(paste("No points have mark j =", j))

  local_k_multi_engine(X=X, I=I, J=J, ..., correction=correction,
                    verbose=verbose, rvalue=rvalue)
}

#' @export
local_l_cross <- function(X, i, j, ..., correction="Ripley", verbose=TRUE, rvalue=NULL) {
  local_k_cross(X=X, i=i, j=j, ..., wantL=TRUE, correction=correction,
              verbose=verbose, rvalue=rvalue)
}

#' Inhomogeneous multitype (cross-type) neighbourhood density function
#'
#' Computes spatially-weighted versions of the local \eqn{K}-function or
#' \eqn{L}-function, defined by Getis and Franklin (1987).
#'
#' The command \code{local_l_cross_inhom} computes the \emph{neighbourhood density function},
#' a local version of the \eqn{L}-function (Besag's transformation of Ripley's
#' \eqn{K}-function) proposed by Getis and Franklin (1987), for 2 types in a multitype
#' spatial point pattern. The command \code{local_k_cross} computes the corresponding
#' local analogue of the cross-type K-function.
#'
#' Given a multitype spatial point pattern \code{X} with types \code{i} and \code{j},
#' the neighbourhood density function
#' \eqn{L_{ij}(r)}{L[ij](r)} associated with the \eqn{i}th point
#' in \code{X} is computed by
#' \deqn{
#'   L_{ij}(r) = \sqrt{\frac a {(n_j) \pi} \sum_j e_{ij}}
#' }{
#'   L[ij](r) = sqrt( (a/((n[j])* pi)) * sum[j] e[i,j])
#' }
#' where the sum is over all points \eqn{j}{j} that lie
#' within a distance \eqn{r} of the \eqn{i}th point,
#' \eqn{\lambda_j}{\lambda[j]} is the estimated intensity of type \eqn{j} of the
#' point pattern at the point \eqn{j}, and \eqn{e_{ij}}{e[i,j]} is an edge correction
#' term (as described in \code{\link{Kest}}).
#' The value of \eqn{L_{ij}(r)}{L[ij](r)} can also be interpreted as one
#' of the summands that contributes to the global estimate of the cross-type \eqn{L}-
#' function.
#'
#' By default, the function \eqn{L_{ij}(r)}{L[ij](r)} or
#' \eqn{K_{ij}(r)}{K[ij](r)} is computed for a range of \eqn{r} values
#' for each point \eqn{i}. The results are stored as a function value
#' table (object of class \code{"fv"}) with a column of the table
#' containing the function estimates for each point \eqn{i} of the pattern
#' \code{X}.
#'
#' Alternatively, if the argument \code{rvalue} is given, and it is a
#' single number, then the function will only be computed for this value
#' of \eqn{r}, and the results will be returned as a numeric vector,
#' with one entry of the vector for each point \eqn{i} of the pattern \code{X}.
#'
#' Computation can be done in parallel by registering a parallel backend for
#' the \code{\link{foreach}} package.
#'
#' @aliases local_l_cross_inhom
#' @seealso \code{\link{localLinhom}}
#'          \code{\link{localKinhom}}
#'          \code{\link{Lcross.inhom}}
#'          \code{\link{Kcross.inhom}}
#' @param X Multitype point pattern (\code{\link{ppp}} object). It must be a multitype
#'      point pattern (or marked point pattern).
#' @param i The type (mark value) of the points in X from which distances are
#'      measured. A character string (or something that will be converted to a
#'      character string). Defaults to the first level of \code{\link{marks}(X)}.
#' @param j The type (marks value) of the points in X to which distances are
#'      measured. A character string (or something that will be converted to a
#'      character string'). Defaults to the second level of \code{\link{marks}(X)}.
#' @param lambdaI Optional. Values of the estimated intensity of the sub-process
#'      of points of type \code{i}. Either a pixel image (object of class \code{"im"}),
#'      a numeric vector containing the intensity values at each of the type
#'      \code{i} points in \code{X}, or a \code{function(x,y)} which can be
#'      evaluated to give the intensity value at any location.
#' @param lambdaJ Optional. Values of the estimated intensity of the sub-process
#'      of points of type \code{j}. Either a pixel image (object of class \code{"im"}),
#'      a numeric vector containing the intensity values at each of the type
#'      \code{j} points in \code{X}, or a \code{function(x,y)} which can be
#'      evaluated to give the intensity value at any location.
#' @param ... Ignored.
#' @param correction String sprecifying the edge correction to be applied.
#'      Option are \code{"none"}, \code{"translate"}, \code{"translation"},
#'      \code{"Ripley"}, \code{"isotropic"} or \code{"best"}. Only one correction
#'      may be given.
#' @param verbose Logical flag indicating whether to print progress reports
#'      during the calculation.
#' @param rvalue Optional. A \emph{single} value of the distance argument \eqn{r}
#'      at which the function L or K should be computed.
#' @return If \code{rvalue} is given, the result is a numeric vector of equal
#'      length to the number of points in X_i.
#'
#'      If the code{rvalue} is absent, the result is an object of class \code{"fv"},
#'      see \code{\link{fv.object}}, which can be plotted directly using
#'      \code{\link{plot.fv}}. Essentially a data frame containing columns:
#'      \item{r}{the vector of values of the argument \eqn{r}
#'      at which the function \eqn{K} has been  estimated }
#'      \item{theo}{the theoretical value \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}
#'      or \eqn{L(r)=r} for a stationary Poisson process
#'      }
#'      together with columns containing the values of the
#'      neighbourhood density function for each point in the pattern.
#'      Column \code{i} corresponds to the \code{i}th point.
#'      The last two columns contain the \code{r} and \code{theo} values.
#' @references Getis, A. and Franklin, J. (1987)
#'      Second-order neighbourhood analysis of mapped point patterns.
#'      \emph{Ecology} \bold{68}, 473--477.
local_k_cross_inhom <-
  function(X, I, J, lambdaI=NULL, lambdaJ=NULL, ..., wantL=FALSE,
           correction="Ripley", verbose=TRUE, rvalue=NULL, lambdaIJ= NULL,
           sigma=NULL, varcov=NULL) {
    spatstat::verifyclass(X, "ppp")

    if(!spatstat::is.multitype(X, dfok=FALSE))
      stop("Point pattern must be multitype")
    marx <- spatstat::marks(X)
    if(missing(i))
      i <- levels(marx)[1]
    if(missing(j))
      j <- levels(marx)[2]

    I <- (marx == i)
    if(!any(I))
      stop(paste("No points have mark i =", i))

    J <- (marx == j)
    if(!any(J))
      stop(paste("No points have mark j =", j))

    dangerous <- c("lambdaI", "lambdaJ", "lambdaIJ")
    dangerI <- dangerJ <- TRUE

    ## intensity data
    if(is.null(lambdaI)) {
      # estimate intensity
      dangerI <- FALSE
      lambdaI <- density(X[I], ..., sigma=sigma, varcov=varcov,
                         at="points", leaveoneout=TRUE)
    } else if(spatstat::is.im(lambdaI)) {
      # look up intensity values
      lambdaI <- spatstat::safelookup(lambdaI, X[I])
    } else if(is.function(lambdaI)) {
      # evaluate function at locations
      XI <- X[I]
      lambdaI <- lambdaI(XI$x, XI$y)
    } else if(is.numeric(lambdaI) && is.vector(as.numeric(lambdaI))) {
      # validate intensity vector
      if(length(lambdaI) != nI)
        stop(paste("The length of", sQuote("lambdaI"),
                   "should equal the number of", Iname))
    } else
      stop(paste(sQuote("lambdaI"), "should be a vector or an image"))

    if(is.null(lambdaJ)) {
      # estimate intensity
      dangerJ <- FALSE
      lambdaJ <- density(X[J], ..., sigma=sigma, varcov=varcov,
                         at="points", leaveoneout=TRUE)
    } else if(spatstat::is.im(lambdaJ)) {
      # look up intensity values
      lambdaJ <- spatstat::safelookup(lambdaJ, X[J])
    } else if(is.function(lambdaJ)) {
      # evaluate function at locations
      XJ <- X[J]
      lambdaJ <- lambdaJ(XJ$x, XJ$y)
    } else if(is.numeric(lambdaJ) && is.vector(as.numeric(lambdaJ))) {
      # validate intensity vector
      if(length(lambdaJ) != nJ)
        stop(paste("The length of", sQuote("lambdaJ"),
                   "should equal the number of", Jname))
    } else
      stop(paste(sQuote("lambdaJ"), "should be a vector or an image"))

    # Weight for each pair
    if(!is.null(lambdaIJ)) {
      dangerIJ <- TRUE
      if(!is.matrix(lambdaIJ))
        stop("lambdaIJ should be a matrix")
      if(nrow(lambdaIJ) != nI)
        stop(paste("nrow(lambdaIJ) should equal the number of", Iname))
      if(ncol(lambdaIJ) != nJ)
        stop(paste("ncol(lambdaIJ) should equal the number of", Jname))
    } else dangerIJ <- FALSE

    danger <- dangerI || dangerJ || dangerIJ
    if(danger) {
      d <- dangerous[c(dangerI, dangerJ, dangerIJ)]
      warning(paste("User specified: ", d))
    }

    local_k_multi_engine(X=X, I=I, J=J, lambdaI=lambdaI, lambdaJ=lambdaJ,
                      ..., wantL=wantL, correction=correction, lambdaIJ=lambdaIJ,
                      verbose=verbose, rvalue=rvalue)
  }

#' @export
local_l_cross_inhom <-
  function(X, i, j, lambdaI=NULL, lambdaJ=NULL, ..., correction="Ripley",
           verbose=TRUE, rvalue=NULL, lambdaIJ=NULL, sigma=NULL, varcov=NULL) {
    local_k_cross_inhom(X=X, i=i, j=j, lambdaI=lambdaI, lambdaJ=lambdaJ, ...,
                     wantL=TRUE, correction=correction, verbose=verbose, rvalue=rvalue,
                     lambdaIJ=lambdaIJ, sigma=sigma, varcov=varcov)
  }

local_k_multi_engine <-
  function(X, I, J, lambdaI=NULL, lambdaJ=NULL, ..., wantL=FALSE,
           correction="Ripley", lambdaIJ=NULL, verbose=TRUE, rvalue=NULL)
  {
    spatstat::verifyclass(X, "ppp")

    npts <- spatstat::npoints(X)
    x <- X$x
    y <- X$y
    W <- X$window
    areaW <- spatstat::area(W)

    I <- spatstat::ppsubset(X, I)
    J <- spatstat::ppsubset(X, J)
    if(is.null(I) || is.null(J))
      stop("I and J must be valid subset indices")

    if(!any(I)) stop("no points belong to subset I")
    if(!any(J)) stop("no points belong to subset J")

    nI <- sum(I)
    nJ <- sum(J)
    lambdaI.ave <- nI/areaW
    lambdaJ.ave <- nJ/areaW

    if(is.null(rvalue))
      rmaxdefault <- spatstat::rmax.rule("K", W, lambdaJ.ave)
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

    #find close pairs
    XI <- X[I]
    XJ <- X[J]
    inpts <- spatstat::npoints(XI)
    close <- spatstat::crosspairs(XI, XJ, rmax)
    # close$i and close$j are serial numbers in XI and XJ respectively;
    # map them to original serial numbers in X
    orig <- seq_len(npts)
    imap <- orig[I]
    jmap <- orig[J]
    iX <- imap[close$i]
    jX <- jmap[close$j]
    # eliminate any identical pairs
    if(any(I & J)) {
      ok <- (iX != jX)
      if(!all(ok)) {
        close$i  <- close$i[ok]
        close$j  <- close$j[ok]
        close$xi <- close$xi[ok]
        close$yi <- close$yi[ok]
        close$xj <- close$xj[ok]
        close$yj <- close$yj[ok]
        close$dx <- close$dx[ok]
        close$dy <- close$dy[ok]
        close$d  <- close$d[ok]
      }
    }
    # extract information for these pairs (relative to orderings of XI, XJ)
    dcloseIJ <- close$d
    icloseI  <- close$i
    jcloseJ  <- close$j

    weighted <- FALSE
    # Form weight for each pair
    if(!is.null(lambdaI) && !is.null(lambdaJ)) {
      weight <- 1/(lambdaI[icloseI] * lambdaJ[jcloseJ])
      weighted <- TRUE
    } else if(!is.null(lambdaIJ)) {
      weight <- 1/lambdaIJ[cbind(icloseI, jcloseJ)]
      weighted <- TRUE
    }

    # initialise
    iseq <- 1:inpts
    icode <- sprintf("%03d", iseq)
    bkt <- function(x) { paste("[", x, "]", sep="") }

    parallel <- check_parallel(verbose)
    fe <- foreach::foreach(i = iseq, .export = c("local_k_calc"))

    switch(correction,
           none={
             # uncorrected! For demonstration purposes only!
             df <- foreach::`%dopar%`(fe, local_k_calc(i, closeI=icloseI, distIJ=dcloseIJ, breaks=breaks,
                                                     weight=if (weighted) weight else NULL,
                                                     pb = if(verbose && !parallel) function(i) spatstat::progressreport(i, n=npts) else NULL))
             # Hack to quickly convert list to data.frame
             attributes(df) <- list(row.names=c(NA_integer_, length(r)),
                                    class="data.frame",
                                    names=make.names(paste("un", icode, sep=""),
                                                     unique=TRUE))
             labl <- paste("%s", bkt(icode), "(r)", sep="")
             desc <- paste("uncorrected estimate of %s",
                           "for point", icode)

             if (!weighted) df <- df / lambdaJ.ave
           },
           translate={
             edgewt <- spatstat::edge.Trans(XI[icloseI], XJ[jcloseJ], paired=TRUE)
             if(weighted) {
               edgewt <- edgewt * weight
             }

             df <- foreach::`%dopar%`(fe, local_k_calc(i, closeI=icloseI, distIJ=dcloseIJ,
                                                     breaks=breaks, weight=edgewt,
                                                     pb = if(verbose && !parallel) function(i) spatstat::progressreport(i, n=npts) else NULL))

             # Hack to quickly convert list to data.frame
             attributes(df) <- list(row.names=c(NA_integer_, length(r)),
                                    class="data.frame",
                                    names=make.names(paste("trans", icode, sep=""),
                                                     unique=TRUE))

             desc <- paste("translation-corrected estimate of %s",
                           "for point", icode)
             labl <- paste("%s", bkt(icode), "(r)", sep="")

             if (!weighted) df <- df / lambdaJ.ave
             h <- spatstat::diameter(W) / 2
             df[r >= h, ] <- NA
           },
           isotropic={
             edgewt <- spatstat::edge.Ripley(XI[icloseI], matrix(dcloseIJ, ncol=1))
             if (weighted) {
               edgewt <- edgewt * weight
             }

             df <- foreach::`%dopar%`(fe, local_k_calc(i, closeI=icloseI, distIJ=dcloseIJ,
                                                     breaks=breaks, weight=edgewt,
                                                     pb = if(verbose && !parallel) function(i) spatstat::progressreport(i, n=npts) else NULL))

             # Hack to quickly convert list to data.frame
             attributes(df) <- list(row.names=c(NA_integer_, length(r)),
                                    class="data.frame",
                                    names=make.names(paste("iso", icode, sep=""),
                                                     unique=TRUE))
             desc <- paste("Ripley isotropic correction estimate of %s",
                           "for point", icode)
             labl <- paste("%s", bkt(icode), "(r)", sep="")
             if (!weighted) df <- df / lambdaJ.ave
             h <- spatstat::diameter(W) / 2
             df[r >= h, ] <- NA
           })

    if (wantL) {
      df <- sqrt(df / pi)
    }

    if(!is.null(rvalue)) {
      nr <- length(r)
      if(r[nr] != rvalue)
        stop("Internal error - rvalue not attained")
      return(as.numeric(df[nr,]))
    }

    # function value table required
    # add r and theo
    if (!wantL) {
      df <- cbind(df, data.frame(r=r, theo=pi * r^2))
      if (!weighted) {
        ylab <- quote(K[ij](r))
        fnam <- "K[ij][',']"
      } else {
        ylab <- quote(Kinhom[ij](r))
        fnam <- "Kinhom[ij][',']"
      }
    } else {
      df <- cbind(df, data.frame(r=r, theo=r))
      if (!weighted) {
        ylab <- quote(L[ij](r))
        fnam <- "L[ij][',']"
      } else {
        ylab <- quote(Linhom[ij](r))
        fnam <- "Linhom[ij][',']"
      }
    }
    desc <- c(desc, c("distance argument r", "theoretical Poisson %s"))
    labl <- c(labl, c("r", "%s[pois](r)"))
    # create fv object
    K <- spatstat::fv(df, "r", ylab, "theo", , alim, labl, desc, fname=fnam)
    # default is to display them all
    spatstat::formula(K) <- . ~ r
    spatstat::unitname(K) <- spatstat::unitname(X)
    attr(K, "correction") <- correction
    return(K)
  }
