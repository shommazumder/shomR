#================================================================
## Conley Standard Errors
## Date: August 21, 2017
## Vectorized Code using Rcpp
#================================================================

pkgs <- c("data.table", "lfe", "geosphere", "Rcpp", "RcppArmadillo")
invisible(sapply(pkgs, require, character.only = TRUE))
sourceCpp("cpp/cpp-functions.cpp")

#' Returns conley SEs corrected for spatial autocorrelation
#'
#' @param reg a regression-like object
#' @param unit character indicating unit
#' @param time character indicating time
#' @param lat character indicating latitude
#' @param lon character indicating longitude
#' @param kernel kernel for weights. defualts to bartlett
#' @param dist_fn distance function. defaults to haversine
#' @param dist_cutoff distance for cutoff. defaults to 500
#' @param lag_cutoff lag of cutoff. defaults to 5
#' @param lat_scale scale. defaults to 111
#' @param verbose verbose output. defaults to FALSE
#' @param cores number of cores to use. defaults to 1
#' @param balanced_pnl logical for whether balanced panel. defaults to FALSE
#'
#'
#' @examples
#' ConleySEs()
#'
#' @export
ConleySEs <- function(reg,
    unit, time, lat, lon,
    kernel = "bartlett", dist_fn = "Haversine",
    dist_cutoff = 500, lag_cutoff = 5,
    lat_scale = 111, verbose = FALSE, cores = 1, balanced_pnl = FALSE) {

    Fac2Num <- function(x) {as.numeric(as.character(x))}
    #source("iterate-obs-function.R", local = TRUE)
    if(cores > 1) {invisible(library(parallel))}

    if(class(reg) == "felm") {
        Xvars <- rownames(reg$coefficients)
        dt = data.table(reg$cY, reg$cX,
            fe1 = Fac2Num(reg$fe[[1]]),
            fe2 = Fac2Num(reg$fe[[2]]),
            coord1 = Fac2Num(reg$clustervar[[1]]),
            coord2 = Fac2Num(reg$clustervar[[2]]))
        setnames(dt,
            c("fe1", "fe2", "coord1", "coord2"),
            c(names(reg$fe), names(reg$clustervar)))
        dt = dt[, e := as.numeric(reg$residuals)]

    } else {
        message("Model class not recognized.")
        break
    }

    n <- nrow(dt)
    k <- length(Xvars)

    # Renaming variables:
    orig_names <- c(unit, time, lat, lon)
    new_names <- c("unit", "time", "lat", "lon")
    setnames(dt, orig_names, new_names)

    # Empty Matrix:
    XeeX <- matrix(nrow = k, ncol = k, 0)

    #================================================================
    # Correct for spatial correlation:
    timeUnique <- unique(dt[, time])
    Ntime <- length(timeUnique)
    setkey(dt, time)

    if(verbose){message("Starting to loop over time periods...")}

    if(balanced_pnl){
        sub_dt <- dt[time == timeUnique[1]]
        lat <- sub_dt[, lat]; lon <- sub_dt[, lon]; rm(sub_dt)

        if(balanced_pnl & verbose){message("Computing Distance Matrix...")}

        d <- DistMat(cbind(lat, lon), cutoff = dist_cutoff, kernel, dist_fn)
        rm(list = c("lat", "lon"))
    }

    if(cores == 1) {
        XeeXhs <- lapply(timeUnique, function(t) iterateObs(sub_index = t,
            type = "spatial", cutoff = dist_cutoff))
    } else {
        XeeXhs <- mclapply(timeUnique, function(t) iterateObs(sub_index = t,
            type = "spatial", cutoff = dist_cutoff), mc.cores = cores)
    }

    if(balanced_pnl){rm(d)}

    # First Reduce:
	XeeX <- Reduce("+",  XeeXhs)

    # Generate VCE for only cross-sectional spatial correlation:
    X <- as.matrix(dt[, eval(Xvars), with = FALSE])
    invXX <- solve(t(X) %*% X) * n

    V_spatial <- invXX %*% (XeeX / n) %*% invXX / n

    V_spatial <- (V_spatial + t(V_spatial)) / 2

    if(verbose) {message("Computed Spatial VCOV.")}

    #================================================================
    # Correct for serial correlation:
    panelUnique <- unique(dt[, unit])
    Npanel <- length(panelUnique)
    setkey(dt, unit)

    if(verbose){message("Starting to loop over units...")}

    if(cores == 1) {
        XeeXhs <- lapply(panelUnique, function(t) iterateObs(sub_index = t,
            type = "serial", cutoff = lag_cutoff))
    } else {
        XeeXhs <- mclapply(panelUnique,function(t) iterateObs(sub_index = t,
            type = "serial", cutoff = lag_cutoff), mc.cores = cores)
    }

	XeeX_serial <- Reduce("+",  XeeXhs)

	XeeX <- XeeX + XeeX_serial

    V_spatial_HAC <- invXX %*% (XeeX / n) %*% invXX / n
    V_spatial_HAC <- (V_spatial_HAC + t(V_spatial_HAC)) / 2

    return_list <- list(
        "OLS" = reg$vcv,
        "Spatial" = V_spatial,
        "Spatial_HAC" = V_spatial_HAC)
    return(return_list)
}




