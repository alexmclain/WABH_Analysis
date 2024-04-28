# Calculation of false-discovery rates


#' Calculate the science-wise FDR (swfdr)
#' 
#' @param pValues Numerical vector of p-values
#' @param truncated Vector of 0s and 1s with indices corresponding to those in pValues; 1 indicates that the p-values is truncated, 0 that it is not truncated
#' @param rounded Vector of 0s and 1s with indices corresponding to those in pValues; 1 indicates that the p-values is rounded, 0 that it is not rounded
#' @param pi0 Initial prior probability that a hypothesis is null (default is 0.5)
#' @param alpha Initial value of parameter alpha from Beta(alpha, beta) true positive distribution (default is 1)
#' @param beta Initial value of parameter beta from Beta(alpha, beta) true positive distribution (default is 50)
#' @param numEmIterations The number of EM iterations (default is 100)
#' 
#' @return pi0 Final value of prior probability - estimated from EM - that a hypothesis is null, i.e. estimated swfdr
#' @return alpha Final value of parameter alpha - estimated from EM - from Beta(alpha, beta) true positive distribution
#' @return beta Final value of parameter beta - estimated from EM - from Beta(alpha, beta) true positive distribution
#' @return z Vector of expected values of the indicator of whether the p-value is null or not - estimated from EM - for the non-rounded p-values (values of NA represent the rounded p-values)
#' @return n0 Expected number of rounded null p-values - estimated from EM - between certain cutpoints (0.005, 0.015, 0.025, 0.035, 0.045, 0.05)
#' @return n Number of rounded p-values between certain cutpoints (0.005, 0.015, 0.025, 0.035, 0.045, 0.05)
#' 
#' @importFrom stats4 mle coef
#' @importFrom stats dbeta lsfit pbeta smooth.spline
#' 
#' @examples
#' pVals <- runif(100)
#' tt <- rr <- rep(0, 100)
#' resSwfdr <- calculateSwfdr(pValues = pVals, truncated = tt, rounded = rr, numEmIterations=100)
#' 
#' @export
calculateSwfdr = function(pValues,truncated,rounded,pi0 = 0.5,alpha=1,beta=50,numEmIterations=100){
  pp = pValues
  tt = truncated
  rr = rounded
  
  
  ll = function(a,b){
    tmp1 = rep(0,length(pp))
    tmp1[tt==0 & rr==0] = dbeta(pp[tt==0 & rr==0],a,b,log=TRUE) - pbeta(0.05,a,b,log.p = TRUE)
    tmp1[tt > 0 & rr==0] = pbeta(pp[tt > 0 & rr==0],a,b,log.p = TRUE) - pbeta(0.05,a,b,log.p = TRUE)
    tmp1 = -sum((1-z)*tmp1, na.rm=TRUE)
    
    probvec = log(pbeta(c(0.05,0.045,0.035,0.025,0.015,0.005),a,b) - pbeta(c(0.045,0.035,0.025,0.015,0.005,0),a,b)) - pbeta(0.05,a,b,log.p = TRUE)
    probvec =  sort(rev(probvec), decreasing = TRUE)
    tmp2 = sum(-n1*c(probvec))
    
    return(tmp1 + tmp2)
  }
  
  ##make vector of rounded p-values
  ppR <- pp[rr > 0]
  ##change values of 0 to 10^-10 and values of 0.05 to 0.05 - 10^-10 
  ##(so we can take the lowest point in "cut" as 0 and highest point as 0.05)
  ppR[ppR == 0] <- 10^-10
  ppR[ppR == 0.05] <- 0.05-10^10-10
  n = table(cut(ppR,c(0,0.005,0.015,0.025,0.035,0.045,0.050)))
  
  
  for(i in 1:numEmIterations){
    
    ## E-step
    
    probvec1 = (pbeta(c(0.05,0.045,0.035,0.025,0.015,0.005),alpha,beta) - pbeta(c(0.045,0.035,0.025,0.015,0.005,0),alpha,beta))/pbeta(0.05,alpha,beta)
    probvec1 = rev(probvec1)
    probvec0 = c(0.005,0.01,0.01,0.01,0.01,0.005)*20
    
    pij0 = pi0*probvec0/(probvec0*pi0 + probvec1*(1-pi0))
    n0 = n*pij0
    n1 = n - n0
    
    z = rep(NA,length(pp))
    z[tt == 0 & rr ==0] <- pi0*20/(pi0*20 + (1-pi0)*dbeta(pp[tt==0 & rr == 0],alpha,beta)/pbeta(0.05,alpha,beta))
    z[tt > 0 & rr ==0] <- pi0*20*pp[tt > 0 & rr ==0]/(pi0*20*pp[tt > 0 & rr==0] + (1-pi0)*pbeta(pp[tt > 0 & rr==0],alpha,beta)/pbeta(0.05,alpha,beta))
    
    ## M-step
    
    pi0 = (sum(n0) + sum(z, na.rm=TRUE))/(sum(n) + sum(rr == 0))
    tmp = stats4::mle(ll,start=list(a=0.05,b=100),lower=c(0.001,1),upper=c(1,500),method="L-BFGS-B")
    alpha = stats4::coef(tmp)[1]
    beta = stats4::coef(tmp)[2]
  }
  return(list(pi0 = pi0, alpha=alpha, beta = beta, z=z,n0=n0,n=n))
}



# Functions for validation of inputs for lm_pi0 and lm_qvalue
#
# All these functions are for internal use for the package only


#' validate pvalues. They must be finite, in range [0,1]
#'
#' @keywords internal
#' @noRd
#' @param p numeric vector of p-values
#'
#' @return numeric vector of length 2 with (min(p), max(p))
check_p <- function(p) {
  if (missing(p)) {
    stop("p is a required argument\n", call. = FALSE)
  }
  if (!is(p, "numeric")) {
    stop("p must be a numeric vector\n", call. = FALSE)
  }
  prange <- range(p, na.rm=TRUE)
  if (prange[1] < 0 | prange[2] > 1) {
    stop("p must be numeric in range [0, 1]\n", call. = FALSE)
  }
  if (any(is.na(p))) {
    stop("p must not have missing values\n", call. = FALSE)
  }
  prange
}


#' validate lambda. They must be numeric, finite, sorted with unique value
#' 
#' @keywords internal
#' @noRd
#' @param x vector of lambda values
#' @param pmax numeric, maximal pvalue
#'
#' @return numeric vector of sorted unique x,
#' or stop if x does not satisfy criteria
check_lambda <- function(x, pmax) {
  if (!is(x, "numeric")) {
    stop("lambda must be a numeric vector \n", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop("lambda must not contain NAs, NULL, or non-finite elements\n",
         call. = FALSE)
  }
  x <- sort(unique(x))
  if (length(x)<4) {
    stop("lambda must be of length >=4\n", call. = FALSE)
  }
  if (min(x)<0 | max(x)>=1) {
    stop("lambda values must all be in range [0, 1)\n", call. = FALSE)
  }
  if (pmax < max(x)) {
    warning("maximal p is smaller than maximal lambda", call. = FALSE)
  }
  x
}


#' validate degrees of freedom. 
#'
#' @keywords internal
#' @noRd
#' @param x expect a single number
#' @param max.value numeric, maximal value allowed for x
#'
#' @return integer derived from x
check_df <- function(x, max.value) {
  if (!is(x, "numeric") & !is(x, "integer")) {
    stop("df must be a number")
  }
  if (length(x) != 1 | any(!is.finite(x))) {
    stop("df must be a single finite number", call. = FALSE)
  }
  x <- round(x)
  if (x <= 1 | x > max.value) {
    stop("df must be in range 1 < df < length(lambda)\n", call. = FALSE)
  }
  x
}


#' validate matrix of covariates.
#'
#' The matrix must be compatible with a vector of pvalues
#'
#' @keywords internal
#' @noRd
#' @param X vector or matrix of covariates
#' @param p vector of p-values
#'
#' @return matrix
check_X <- function(X, p) {
  # allow for null input (no covariates)
  if (missing(X)) {
    X <- NULL
  }
  if (is.null(X)) {
    warning("X is missing or NULL - modeling will have no effect",
            call.=FALSE)
    X <- cbind(rep(1, length(p)))
    rownames(X) <- names(p)
  }
  # allow for a single covariates specified as a vector
  if (is.null(dim(X))) {
    X <- cbind(X=X)
  }
  # ensure that X and pvalues are compatible
  if (!is(X, "matrix")) {
    if (is(X, "data.frame")) {
      X <- as.matrix(X)
    }
  }
  if (length(p) != nrow(X)) {
    stop("incompatible X and p - different lengths\n", call. = FALSE)
  }
  if (!is.null(names(p)) & !identical(rownames(X), names(p))) {
    stop("X and p have different names", call. = FALSE)
  }
  # ensure that all columns in X are numeric
  is.number.class = function(z) {
    is(z, "numeric") | is(z, "integer") | is(z, "factor")
  }
  if (!all(apply(X, 2, is.number.class))) {
    stop("X must be a numeric vector or numeric matrix\n", call. = FALSE)
  }
  # ensure that all values are set
  num.bad <- sum(!is.finite(X))
  if (num.bad>0) {
    stop("X must not contain missing or non-finite values\n", call. = FALSE)
  }
  X
}


#' check if an object is of a certain class
#'
#' @keywords internal
#' @noRd
#' @param x object
#' @param classname character
#' 
#' @return nothing, emit error if check not satisfied
check_class <- function(x, classname) {
  if (!is(x, classname)) {
    stop(paste0("object is not of class '", classname, "'\n"), call. = FALSE)
  }
}


#' Get number of decimals (i.e. return total number of digits after decimal point) for any vector of numbers in [0,1) if number of decimals <= 6
#' 
#' @param x Numerical vector where all elements are in [0,1)
#' 
#' @return Vector giving the number of decimals for each element in x if the number is <= 6; otherwise return 7 with a warning
#'
#' @examples
#' get_number_decimals(c(0.0006, 0.0750, 0.0420, 0.0031, 0.0001, 0.0100))
#' get_number_decimals(c(6*10^-4, 7.5*10^-2, 4.2*10^-2, 3.1*10^-3, 10^-4, 10^-2))
#' get_number_decimals(c(6.5*10^-4, 0.0100)) 
#' get_number_decimals(c(6.5e-4, 0.0100))
#' get_number_decimals(c(0.00065, 0.0100))
#' get_number_decimals(c(10^-7, 10e-7, 10e-3))
#'
#' @export
get_number_decimals <- function(x)
{
  if((any(x<0))|(any(x>=1)))
  {
    stop("All elements of x should be in [0,1)")
  }
  
  ##get the maximum number of digits
  max_digits <- 6
  
  ##get all numbers spaced 10^-k apart from 0 to 1
  list_grid <- lapply(1:max_digits, function(k){(1:10^k)/(10^k)})
  
  ##function for a single number
  n_dec_single <- function(x_single){
    ##round to get rid of funny numerical issues
    x_single <- round(x_single, 12)
    
    if(x_single < 10^-max_digits)
    {
      n_dec <- max_digits + 1
      warning(paste(max_digits + 1, " is a place-holder. Number of decimals seems to be >", 
                    max_digits, ". This case is not implemented. Beware floating point arithmetic!",
                    sep=""))
    } else {
      ##get which vector the query number is in, which corresponds to the "number of digits"
      grid_x_is_in <- sapply(list_grid, function(l,a){min(abs(l-a)) < 10^-12}, x_single)
      
      if(x_single==0)
      {
        n_dec <- 0
      } else {
        if(sum(grid_x_is_in) >= 1)
        {
          n_dec <- min((1:max_digits)[grid_x_is_in])
        } else {
          n_dec <- max_digits + 1
          warning(paste(max_digits + 1, " is a place-holder. Number of decimals seems to be >", 
                        max_digits, ". This case is not implemented. Beware floating point arithmetic!",
                        sep=""))        
        }
      }
    }
    n_dec
  }
  n_d <- sapply(x, n_dec_single) 
  n_d
}

# Computation of pi0 from pvalues and matrix of covariates
#


#' Estimation of pi0, proportion of p-values consistent with a null hypothesis
#'
#' @param p numeric vector, p-values
#' @param lambda numeric vector, thresholds used to bin pvalues,
#' must be in [0,1).
#' @param X numeric matrix, covariates that might be related to p values
#' (one test per row, one variable per column). 
#' @param type character, type of regression used to fit features to pvalues
#' @param smooth.df integer, degrees of freedom when estimating pi0(x) with
#' a smoother.
#' @param threshold logical, if TRUE, all estimates are thresholded into unit
#' interval; if FALSE, all estimates are left as they are computed
#' @param smoothing character, type of smoothing used to fit pi0
#'
#' @return pi0 numerical vector of smoothed estimate of pi0(x).
#' The length is the number of rows in X.
#' @return pi0.lambda numeric matrix of estimated pi0(x) for each value of
#' lambda. The number of columns is the number of tests, the number of rows is
#' the length of lambda.
#' @return lambda numeric vector of thresholds used in calculating pi0.lambda
#' @return pi0.smooth (only output with smoothing="smooth.spline") Matrix of
#' fitted values from the smoother fit to the pi0(x) estimates at each value
#' of lambda (same number of rows and columns as pi0.lambda)
#'
#' @importFrom stats binomial glm
#'
#' @examples
#' # define a covariate
#' X <- seq(-1,2,length=1000)
#' # set probability of being null
#' pi0 <- 1/4*X + 1/2
#' # generate null/alternative p-values
#' nullI <- rbinom(1000,prob=pi0,size=1)> 0
#' # vector of p-values
#' pValues <- rep(NA,1000) 
#' pValues[nullI] <- runif(sum(nullI)) # from U(0,1)
#' pValues[!nullI] <- rbeta(sum(!nullI),1,2) # from Beta
#' pi0x <- lm_pi0(pValues, X=X)
#'
#' @export
lm_pi0 <- function(p, lambda = seq(0.05, 0.95, 0.05), X,
                   type=c("logistic", "linear"),
                   smooth.df=3,
                   threshold=TRUE,
                   smoothing=c("unit.spline", "smooth.spline")) {
  
  # check validity of inputs
  type <- match.arg(type)
  smoothing <- match.arg(smoothing)
  prange <- check_p(p)
  lambda <- check_lambda(lambda, prange[2])
  n.lambda <- length(lambda)
  smooth.df <- check_df(smooth.df, n.lambda)
  X <- check_X(X=X, p=p)
  
  # pick a modeling function, fit odels for each lambda
  available.regressions <- list(logistic=fit_logistic, linear=fit_linear)
  fit.function <- available.regressions[[type]]
  pi0.lambda <- matrix(NA, nrow=nrow(X), ncol=n.lambda)
  for (i in 1:n.lambda) {
    y <- (p >= lambda[i])
    pi0.lambda[, i] <- fit.function(y, X)/(1-lambda[i])
  }
  if (threshold) {
    pi0.lambda <- regularize.interval(pi0.lambda)
  }
  
  # smooth over values of lambda (for each p-value/ row in X)
  # (instead of taking limit lambda->1, use largest available lambda)
  available.smoothings <- list(smooth.spline=smooth.spline.pi0,
                               unit.spline=unit.spline.pi0)
  smoothing.function <- available.smoothings[[smoothing]]
  pi0 <- smoothing.function(lambda, pi0.lambda, smooth.df)
  if (threshold) {
    pi0 <- regularize.interval(pi0)
  }
  
  result <- c(list(call=match.call(), lambda=lambda, X.names = colnames(X),
                   pi0.lambda=pi0.lambda), pi0)
  
  class(result) <- "lm_pi0"
  result
}





# #############################################################################
# modeling of a response vector with covariates


#' Fit response values using a binomial/logit model
#'
#' @keywords internal
#' @noRd
#' @param y numeric vector
#' @param X numeric matrix (covariates)
#'
#' @return numeric vector
#'
#' @importFrom stats binomial glm
fit_logistic <- function(y, X) {
  glm(y ~ X, family=binomial)$fitted.values
}


#' Fit response values using a linear regression
#' 
#' (This could be implemnted via glm(... family=gaussian) but the
#' implementation with lsfit is much faster
#'
#' @keywords internal
#' @noRd
#' @param y numeric vector
#' @param X numeric matrix (covariates)
#'
#' @return numeric vector
#'
#' @importFrom stats lsfit
fit_linear <- function(y, X) {
  regFit <- lsfit(X, y)$coefficients
  regFit[1] + (X %*% matrix(regFit[-1], ncol=1))
}




# #############################################################################
# other helper functions


#' Force a set of values in a regular interval
#'
#' (This is a simpler/faster implementation than with ifelse)
#'
#' @keywords internal
#' @noRd
#' @param x numeric vector or matrix
#' @param interval numeric vector of length 2 with min/max values
#'
#' @return same object like x, with values truncated by [0,1]
regularize.interval <- function(x, interval = c(0, 1)) {
  if (is(x, "list")) {
    return (lapply(x, regularize.interval, interval=interval))
  }
  x[x < interval[1]] <- interval[1]
  x[x > interval[2]] <- interval[2]
  x
}
# Esimation of qvalues from pvalues and matrix of covariates


#' Compute qvalues taking into account a matrix of covariates
#'
#' The recipe for turning pvalues into qvalues is adapted from package
#' 'qvalue' and articles by Storey, Tibshirani, Taylor, Siegmund.
#' 
#' @param p numeric vector of p-values
#' @param X matrix of covariates (can be missing if pi0 is specified instead)
#' @param pfdr logical, making estimates robust for small p-values and a small
#' sample size
#' @param pi0 list with pi0 estimates from lm_pi0
#' @param smoothing character, type of smoothing used to fit pi0. Note the
#' default in this function is different than in lm_pi0.
#' @param ... other parameters (passed on to lm_pi0 if pi0 is not provided)
#' 
#' @return list
#'
#' @examples
#' # define a covariate
#' X <- rep(c(0, 1), each=1000)
#' # generate p-values, randomly for group 0 and with low values for group 1
#' pVal <- c(runif(1000), rbeta(1000, 0.2, 1))
#' # compute qvalues, using the covariate
#' qVal <- lm_qvalue(pVal, X=X)
#' 
#' @export
lm_qvalue <- function(p, X, pfdr=FALSE, pi0=NULL, 
                      smoothing=c("unit.spline", "smooth.spline"), ...) {
  
  # check inputs
  prange <- check_p(p)
  smoothing <- match.arg(smoothing)
  
  # pi0 is required - compute it if not provided
  if (is.null(pi0)) {
    pi0 <- lm_pi0(p, X=X, smoothing=smoothing, ...)
  }
  
  # block to compute qvalues (adapted from package qvalue)
  # This is almost the same as pi0s*p.adjust(p, method="BH")
  # However, this implementation allows setting pfdr=TRUE similarly
  # as in qvalue::qvalue()
  n <- length(p)
  i <- n:1
  o <- order(p, decreasing=TRUE)
  ro = order(o, decreasing=FALSE)
  if (pfdr) {
    q <- pmin(1, cummin( p[o]*n/ (i*(1- (1-p[o])^n))))
  } else {
    q <- pmin(1, cummin( p[o]*n/ i ))
  }
  q <- (q* pi0$pi0[o])[ro]
  
  # create output
  if ("call" %in% names(pi0)) {
    pi0 <- pi0[-which(names(pi0)=="call")]
  }
  result <- c(list(call=match.call()), pi0,
              list(pvalues = p, qvalues = q))
  class(result) <- "lm_qvalue"
  result
}

# Functions for display of package objects via print()


#' Display a summary of object lm_pi0
#'
#' @keywords internal
#' @noRd
#' @param x object of class lm_pi0
#' @param ... other arguments ignored
#'
#' @method print lm_pi0
#' @export
print.lm_pi0 <- function(x, ...) {
  check_class(x, "lm_pi0")
  components <- get.components(c("call", "lambda", "X", "pi0"),
                               list(...)[["components"]])
  compound.message(x, components)
  invisible(x)
}


#' Display a summary of an lm_qvalue object
#'
#' @keywords internal
#' @noRd
#' @param x lm_qvalue object
#' @param ... ignored
#'
#' @method print lm_qvalue
#' @export
print.lm_qvalue <- function(x, ...) {
  check_class(x, "lm_qvalue")
  components <- get.components(c("call", "lambda", "X", "pi0", "hits"),
                               list(...)[["components"]])
  compound.message(x, components)
  invisible(x)
}


#' Display a summary of object lm_pi0
#'
#' (This is the same as print.lm_pi0)
#'
#' @keywords internal
#' @noRd
#' @param object object of class lm_pi0
#' @param ... other arguments ignored
#'
#' @method summary lm_pi0
#' @export
summary.lm_pi0 <- function(object, ...) {
  invisible(print.lm_pi0(object))
}


#' Display a summary of object lm_qvalue
#'
#' (This is the same as print.lm_qvalue)
#'
#' @keywords internal
#' @noRd
#' @param object object of class lm_qvalue
#' @param ... other arguments ignored
#'
#' @method summary lm_qvalue
#' @export
summary.lm_qvalue <- function(object, ...) {
  invisible(print.lm_qvalue(object))
}




# #############################################################################
# Some helper functions to the print() and summary()


#' Helper to create a single string by concatenating items from a vector
#'
#' @keywords internal
#' @noRd
#' @param x vector of things
#' @param width vector of character widths for each item in x
#'
#' @return a single string
v2s <- function(x, width=8) {
  empty <- paste(rep(" ", max(width)), collapse="")
  xlen <- length(x)
  if (length(width)<xlen) {
    width <- rep(width, length=xlen)[1:xlen]
  }
  if (is(x, "numeric")) {
    x <- as.character(round(x, 4))
  }
  result <- as.character(x)
  for (i in seq_along(x)) {
    ichars <- nchar(x[i])
    if (ichars < width[i]) {
      result[i] <- paste0(substr(empty, 1, width[i]-ichars), result[i])
    }
  }
  paste(result, collapse=" ")
}


#' Get a set intersection, but when second set is null default to first set
#'
#' This is meant to identify a subset of supported features that are requested
#'
#' @keywords internal
#' @noRd
#' @param supported vector of supported feature names
#' @param requested vector of requested feature names
#'
#' @return character vector with an intersection, or all supported features
get.components <- function(supported, requested) {
  if (is.null(requested)) {
    return(supported)
  }
  intersect(supported, requested)
}


#' Compose a two line report about a numeric vector
#'
#' @keywords internal
#' @noRd
#' @param v numeric vector
#'
#' @return vector with two strings a header line and a data line
compose.stats <- function(v) {
  header <- c("(Length)", "Min", "Mean", "Median", "Max")
  data <- c(length(v), min(v), mean(v), median(v), max(v))
  c(v2s(header), v2s(data))
}


#' Compose and output a compound message and output
#'
#' @keywords internal
#' @noRd
#' @param x list object of type lm_qvalue or lm_pi0 (not checked)
#' @param components character vector, identifiers suggesting what to include
#' in output
#'
#' @return character vector
compound.message <- function(x, components=c("call", "lambda",
                                             "X", "pi0", "hits")) {
  comps = components
  output = setNames(vector("list", length=length(comps)), comps)
  if ("call" %in% comps) {
    output$call = c("Call:", x$call)
  }
  if ("lambda" %in% comps) {
    output$lambda =   c("lambda:", compose.stats(x$lambda))
  }
  if ("X" %in% comps) {
    output$X = c("covariates:", "(Length)", v2s(length(x$X.names)))
  }
  if ("pi0" %in% comps) {
    output$pi0 = c("pi0:", compose.stats(x$pi0))
  }
  if ("hits" %in% comps) {
    hits.header <- c(" ", "<1e-4", "<1e-3", "<0.01", "<0.05", "<0.1", "<1")
    hits.widths = c(9, rep(7, 6))
    thresholds = c(1e-4, 1e-3, 1e-2, 0.05, 0.1, 1)
    hits.p <- c("p-value",
                vapply(thresholds,
                       function(t) { sum(x$pvalues<t) },
                       integer(1)))
    hits.q <- c("q-value",
                vapply(thresholds,
                       function(t) { sum(x$qvalues<t) },
                       integer(1)
                ))
    output$hits <- c("Cumulative number of significant calls:",
                     v2s(hits.header, hits.widths),
                     v2s(hits.p, hits.widths),
                     v2s(hits.q, hits.widths))
  }
  # add a separator line
  output = lapply(output, function(x) { c(x, "") })
  message(paste(c("", unlist(output)), collapse="\n")) 
}

# Estimation of a pi0 curves using splines 
# This is used in lm_pi0 to estimate pi0 based on estimates at various lambda


#' Fit smoothing splines to several curves 
#'
#' This implementation uses smooth.spline. It gives a slight optimization
#' by not reporting original data (keep.data=FALSE) and by setting numeric
#' tolerance to a fixed number (tol=1e-7). In the latter, the precise number
#' does not matter, but it is important that it is not computed from x
#' at each iteration in the loop.
#'
#' @keywords internal
#' @noRd
#' @param x numeric vector, position of knots
#' @param ymat numeric matrix, ncol(ymat) should match length(x)
#'
#' @return list with components pi0.smooth and pi0
smooth.spline.pi0 <- function(x, ymat, df=3) {
  nx <- length(x)
  pi0.smooth = matrix(NA, nrow=nrow(ymat), ncol=nx)
  pi0 <- rep(NA, nrow(ymat))
  for (i in seq_len(nrow(ymat))) {
    yfit <- smooth.spline(x, ymat[i,], df=df, tol=1e-7, keep.data=FALSE)$y
    pi0.smooth[i,] = yfit
    pi0[i] <- yfit[nx]
  }
  list(pi0.smooth=pi0.smooth, pi0=pi0)
}


#' Fit smoothing unit splines to several curves
#'
#' In contrast to smooth.spline.pi0, this function uses the bs() bases
#' and uses boundary knots at x=(0, 1)
#'
#' @keywords internal
#' @noRd
#' @param x numeric vector, position of knots
#' @param ymat numeric matrix, ncol(ymat) should match length(x)
#'
#' @return list with component pi0 (numeric vector)
unit.spline.pi0 <- function(x, ymat, df=3) {
  transform.matrix <- unit.spline.matrix(x, df=df)
  ## make prediction only from the last row of the transformation matrix
  transform.last <- transform.matrix[nrow(transform.matrix), ]
  pi0 <- apply(ymat, 1, function(z) { sum(transform.last * z) })
  list(pi0=pi0)
}


#' Construct a transformation matrix for a unit spline
#'
#' This function uses library::bs to construct a b-spline basis
#' from knots, using [0,1] as the boundary knots.
#'
#' @keywords internal
#' @noRd
#' @param x numeric vector, position of internal knots within [0,1]
#'
#' @importFrom splines bs
#'
#' @return matrix of dimension c(length(x), length(x))
unit.spline.matrix <- function(x, df=3) {
  # construct a splines::bs basis from x
  x.bs = bs(x, knots=x, Boundary.knots=c(0, 1))
  x.bsTbs = t(x.bs) %*% x.bs
  n.bs = ncol(x.bs)
  n.x = length(x)
  
  # construct a smoothing component
  omega.sqrt = matrix(0, ncol=n.bs, nrow=n.bs)
  for (i in seq(3, n.bs)) {
    omega.sqrt[i, c(i-2, i-1, i)] = c(1,-2,1)
  }
  omega = t(omega.sqrt) %*% omega.sqrt
  
  # optimize smoothing weight to get number of degrees of freedom right
  multiplier = optimize.multiplier(x.bsTbs, omega, target=2+df)
  
  # output a transformation matrix
  x.bs %*% solve(x.bsTbs + multiplier*omega) %*% t(x.bs)
}


#' Fit smoothing unit spline to one curve
#'
#' @keywords internal
#' @noRd
#' @param x numeric vector
#' @param y numeric vector
#'
#' @return vector of same length as y
unit.spline <- function(x, y, df=3) {
  transform.matrix = unit.spline.matrix(x, df=df)
  as.numeric(transform.matrix %*% cbind(y))
}


#' recursive internal function to find an optimal value of lambda/multiplier
#' so that trace(X+lambda*omega) = target
#'
#' This implementation relies on the knowledge that a larger lagrange
#' multiplier leads to lower trace(X+lambda*omega).
#'
#' The implementation starts in an interval (0, 1, 2), which are lower,
#' middle, and upper bounds, and then expands the interval or zooms in to
#' find a reasonable multiplier
#'
#' @keywords internal
#' @noRd
#' @param X matrix
#' @param omega matrix
#' @param target numeric, target for trace(X+lambda*omega)
#' @param interval numeric vector of length 2, current range for lambda
#' (internal use)
#' @param values numeric vector of length 2, current target estimates
#' (internal use)
#' @param tol numeric, numerical tolerance, does not need to be very small
#'
#' @return numeric, lagrange multiplier that brings criterion close to its
#' target
optimize.multiplier <- function(X, omega, target=5,
                                interval=c(1e-6, 1, 2), values=c(NA, NA, NA),
                                tol=1e-4) {
  # exit early if the interval is very narrow
  interval.width <- interval[3]-interval[1]
  if (interval.width/interval[2] < tol) {
    return(interval[2])
  }
  # estimate values for the target based on the test interval range
  for (i in c(1,2,3)) {
    if (is.na(values[i])) {
      values[i] <- sum(diag(solve(X + interval[i]*omega)))
    }
  }
  # exit if satisfactory solution is found
  if (abs(values[1] - values[3]) < tol) {
    return(mean(interval))
  }
  # adjust the interval range and look further
  if (target < values[3]) {
    interval.new <- interval[3] + c(0, 2*interval.width, 4*interval.width)
    values.new <- c(values[2], NA, NA)
  } else {
    if (target < values[2]) {
      interval.new <- c(interval[2], mean(interval[2:3]), interval[3])
      values.new <- c(values[2], NA, values[3])      
    } else {
      interval.new <- c(interval[1], mean(interval[1:2]), interval[2])
      values.new <- c(values[1], NA, values[2])
    }
  }
  optimize.multiplier(X, omega, target, interval.new, values.new)
}

