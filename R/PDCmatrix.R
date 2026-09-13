

##### PDC function based on a multivariate arima model  ######
# By Guillermo Granados
# Department of Mathematics and Statistics Lancaster University

#' Quadratic multiplication of a matrix M with respect to a matrix A:
#'  Conj(M)  A M, where Conj() is the complex conjugate function

#' @param M A Matrix of dimension P by P
#' @param A A squared Matrix of dimension P by P
#' @return The squared root of the absolute values of the matrix result of the 
#' operation Conj(M) A M 
#'  
#' @export
#' @examples
#' M=matrix( rnorm(100), ncol=10  )
#' A=matrix( rnorm(100), ncol=10  )
#' sqnorms(M, A)
sqnorms<-function(M,A){ sqrt(abs( Conj( M) %*% A %*% M ) )  }




#' Fit a vector autoregressive model of order p by ordinary least squares
#'
#' Internal replacement for the AR-only estimation previously delegated to
#' the archived \pkg{marima} package (archived from CRAN 2026-05-18 because
#' its own dependency, 'marima', was itself archived).
#'
#' \strong{Provenance / citation.} This is not a from-scratch guess at what
#' marima computed -- it was derived by reading marima's actual estimation
#' code, obtained from the official read-only CRAN-on-GitHub mirror:
#' Spliid, H. (2017) \emph{marima: Multivariate ARIMA and ARIMA-X Analysis},
#' version 2.2, source: \url{https://github.com/cran/marima}, file
#' \code{R/marima.R} (as retrieved 2026-09; the corresponding CRAN archive
#' copy is \url{https://cran.r-project.org/src/contrib/Archive/marima/marima_2.2.tar.gz}).
#'
#' In \code{marima()} (see \code{R/marima.R}, roughly lines 517-620), when no
#' \code{ma.pattern} is supplied the model reduces to, for each variable
#' \code{i}, an OLS call \code{lm(x_i ~ -1 + <lagged regressors>)}; the
#' fitted coefficients are then stored as
#' \code{ar.estimates[i, ] <- -MOD$coefficients} (line ~613 negates the
#' regression coefficient). This sign flip exists because marima's internal
#' "indicator form" is \code{A_0 x_t + A_1 x_(t-1) + ... = e_t} with
#' \code{A_0 = I}, so a regression of \code{x_t[i]} on the lagged terms
#' naturally yields \code{-A_lag[i,j]}, which marima negates back to
#' \code{A_lag[i,j]}.
#'
#' This package's \code{matrix_PDC()} then applied a second negation,
#' \code{arest <- -short.form(fit$ar.estimates, leading = FALSE)}. Composing
#' both sign flips: \code{arest[i,j,lag] = -A_lag[i,j]}, which is exactly the
#' plain OLS coefficient of \code{x_(t-lag)[j]} on \code{x_t[i]} -- i.e. the
#' standard VAR coefficient \code{Phi_lag[i,j]} in
#' \code{x_t = Phi_1 x_(t-1) + ... + Phi_p x_(t-p) + e_t}. Because the
#' \code{ar.pattern} built by \code{define.model(kvar, ar)} (as used by
#' \code{matrix_PDC}) is a full, unrestricted pattern -- every variable
#' regressed on every lagged variable, same regressors for every equation
#' \code{i} -- marima's per-equation \code{lm()} calls are mathematically
#' identical to one joint multivariate least-squares solve. That is what
#' this function computes directly, with no external dependency, and it
#' returns \code{Phi} in the exact orientation \code{matrix_PDC} previously
#' read off \code{-marima::short.form(fit$ar.estimates, leading = FALSE)}.
#'
#' @param unit A matrix with time in rows and each variable (time series) in
#' a column, matching the convention used throughout this package.
#' @param order Integer AR order p (the model uses contiguous lags 1..p).
#' @return A list with \code{Phi}, an array of dimension
#' \code{c(kvar, kvar, order)} where \code{Phi[,,lag]} is the coefficient
#' matrix for that lag such that
#' \code{x_t = Phi_1 x_(t-1) + ... + Phi_p x_(t-p) + e_t}; and
#' \code{resid}, the regression residuals.
#' @keywords internal
fit_var_ols <- function(unit, order) {
  unit <- as.matrix(unit)
  n <- nrow(unit)
  kvar <- ncol(unit)

  # marima's default (means = 1) mean-adjusts all variables before fitting
  mu <- colMeans(unit)
  X <- scale(unit, center = mu, scale = FALSE)

  neff <- n - order
  # design matrix of lagged regressors: [x_(t-1) | x_(t-2) | ... | x_(t-p)]
  Z <- matrix(0, nrow = neff, ncol = kvar * order)
  for (lag in 1:order) {
    Z[, ((lag - 1) * kvar + 1):(lag * kvar)] <-
      X[(order - lag + 1):(n - lag), , drop = FALSE]
  }
  Y <- X[(order + 1):n, , drop = FALSE]

  # multivariate OLS: B = (Z'Z)^-1 Z'Y, B is (kvar*order) x kvar
  B <- solve(t(Z) %*% Z, t(Z) %*% Y)

  Phi <- array(0, dim = c(kvar, kvar, order))
  for (lag in 1:order) {
    block <- B[((lag - 1) * kvar + 1):(lag * kvar), , drop = FALSE]
    Phi[, , lag] <- t(block)
  }

  list(Phi = Phi, resid = Y - Z %*% B)
}


#' Partial directed coherence matrix    

#' @param unit A Matrix containing the multivariate time series. Each column 
#' represents a univariate time series.
#' @param ar  Integer vector containing all the lags considered for the
#' vector autoregressive model
#' @return An real array of dimensions, ncol(unit), ncol(unit), n, where n is the number of 
#' frequencies at which the PDC is estimated.   
#'  
#' @export
#' @examples
#' X=matrix( rnorm(2000), ncol=10  )
#' ar=c(1, 2)
#' matrix_PDC(X, ar)
matrix_PDC=function(  unit, ar ){
  kvar=ncol(unit)
  order=length(ar)
  # Previously: model = marima::define.model(kvar=kvar, ar=ar)
  #             fit   = marima::marima(unit, ar.pattern = model$ar.pattern, penalty = 2)
  #             arest = -marima::short.form(fit$ar.estimates, leading = FALSE)
  # 'marima' was archived from CRAN 2026-05-18 (its own dependency was archived
  # first) and its author was unreachable, so it is replaced here with an
  # in-package OLS VAR(p) fit. See fit_var_ols() above for the full derivation
  # and citation of marima's original source code showing this is an exact
  # reproduction (not an approximation) of marima's pure-AR estimation path.
  fit <- fit_var_ols(unit, order = order)
  arest <- fit$Phi
  omegas=TSA::periodogram(unit[,1], plot = F)$freq  
  #PDC in complex values with a certain level of precision 
  dimens<-dim(arest)[1]
  AF<-array(1:length(omegas), c(dimens,dimens,length(omegas) ))  
  for(j in 1:length(omegas)){
    AF[,,j]<- diag(rep(1,dimens )) 
    for(i in 1:order){
      mycom<- complex( real = cos(2*pi*omegas[j]), imaginary = sin(2*pi*omegas[j]) )^i
      AF[,,j]=AF[,,j] - arest[,,i]*mycom
    }# end i
  }#end j
  # PDC in real values
  dimens<-dim(AF)[1]
  PDCmat<-array(1:length(omegas), c(dimens,dimens,length(omegas) ))
  for(j in 1:length(omegas)){  
    #solve(fit$resid.cov)   diag( rep(1,dimens) )
    sq<- apply ( AF[,,j], 2, sqnorms, A=diag( rep(1,dimens) ) )   
    PDCmat[,,j]<-  t( t(AF[,,j] )/ sq )
  }
  return(PDCmat)
}
