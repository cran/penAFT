\name{penAFT.plot}
\alias{penAFT.plot}
\title{Plot cross-validation curves}
\description{A function for plotting the cross-validation curves for the regularized semiparametric accelerated failure time model estimator.}
\usage{
penAFT.plot(fit)
}

\arguments{
\item{fit}{A fitted model from \code{penAFT.cv}.}
}

\value{No return value; prints a plot of cross-validation curves as described in Details.}

\details{
	This function returns a plot with the cross-validation curves for the regularized Gehan estimator.  The vertical blue line indicates the tuning parameter which minimized cross-validated linear predictor scores and the vertical black line indicates the tuning parameter which minimized the cross-validated Gehan loss according to the one-standard-error rule. The vertical axis (and blue line) denotes the cross-validated linear predictor scores whereas the right vertical axis (and black line) denotes cross-validated Gehan loss and standard errors). To make matters simple, we do not allow for customization of the plot: please refer to the source code if extensive customization is desired. 
}

\examples{
\donttest{
  # --------------------------------------
  # Generate data  
  # --------------------------------------
  set.seed(1)
  genData <- genSurvData(n = 50, p = 100, s = 10, mag = 2, cens.quant = 0.6)
  X <- genData$X
  logY <- genData$logY
  delta <- genData$status
  
  
  # ------------------------------------------
  # Fit elastic net penalized estimator with CV
  # -------------------------------------------
  fit.cv <- penAFT.cv(X = X, logY = logY, delta = delta,
                    nlambda = 50,
                    penalty = "EN", tol.rel = 1e-5,
                    alpha = 1, nfolds = 10)
  penAFT.plot(fit.cv)
}

}