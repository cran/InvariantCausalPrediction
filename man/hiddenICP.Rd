\name{hiddenICP}
\alias{hiddenICP}
\title{
Invariant Causal Prediction with hidden variables
}
\description{
Confidence intervals for causal effects in a regression setting with possible confounders.
}
\usage{
hiddenICP(X, Y, ExpInd, alpha = 0.1, mode = "asymptotic", intercept=TRUE)
}
\arguments{
  \item{X}{
  A matrix (or data frame) with the predictor variables for all experimental settings
}
  \item{Y}{
 The response or target variable of interest.   Can be numeric for regression or a factor with two levels for binary classification.
}
  \item{ExpInd}{
  Indicator of the experiment or the intervention type an observation belongs to. 
  Can be a numeric vector of the same length as \code{Y} with \code{K}
  unique entries if there are \code{K} different experiments 
  (for example entry 1 for all observational data and entry 2 for
  intervention data). 
  Can also be a list, where each element of the list contains the
  indices of all observations that belong to the corresponding grouping
  in the data (for example two elements: first element is a vector that contains indices
  of observations that are observational data and 
  second element is a vector that contains indices 
  of all observations that are of interventional type).
}
  \item{alpha}{
  The level of the test procedure. Use the default \code{alpha}=0.1 
  to obtain 90\% confidence intervals.
}
  \item{mode}{
  Currently only mode "asymptotic" is implemented; the argument is thus in the current version without effect.
}
  \item{intercept}{
  Boolean variable; if \code{TRUE}, an intercept is added to the design matrix (but coefficients are returned without intercept term).
}
}

\value{
A list with elements
  \item{betahat}{The point estimator for the causal effects}
  \item{maximinCoefficients}{The value in the confidence interval for each variable effects that is 
    closest to 0. Is hence non-zero for variables with significant effects.}
}
\references{
    none yet...
}
\author{
Nicolai Meinshausen <meinshausen@stat.math.ethz.ch>
}

\examples{
 ##########################################
 ####### 1st example:
 ####### Simulate data with interventions
      set.seed(1)
    ## sample size n
      n <- 2000
    ## 3 predictor variables
      p  <- 3
    ## simulate as independent Gaussian variables
      X <- matrix(rnorm(n*p),nrow=n)
    ## divide data into observational (ExpInd=1) and interventional (ExpInd=2)
      ExpInd <- c(rep(1,n/2),rep(2,n/2))
    ## for interventional data (ExpInd==2): change distribution
      X[ExpInd==2,] <- sweep(X[ExpInd==2,],1, 6*rnorm(p) ,FUN="*")
      ## add hidden variables
      W <- rnorm(n) * 10
      X <- X + outer(W, rexp(p))
      
      ## first two variables are the causal predictors of Y
      beta <- c(1,1,0)
    ## response variable Y
      Y <- as.numeric(X\%*\%beta + W + rnorm(n))
      

####### Compute "hidden Invariant Causal Prediction" Confidence Intervals
      icp <- hiddenICP(X,Y,ExpInd,alpha=0.01)

 ###### Print point estimates and points in the confidence interval closest to 0
      print(icp$betahat)
      print(icp$maximinCoefficients)
      cat("true coefficients are:", beta)

 #### compare with coefficients from a linear model
      cat("coefficients from linear model:")
      print(coef( lm(Y ~ X-1)))



      
##########################################
####### 2nd example:
####### Simulate model X -> Y -> Z with hidden variables, trying to
######  estimate causal effects from (X,Z) on Y
      set.seed(21)
    ## sample size n
      n <- 2000
    ## simulate as independent Gaussian variables
      W <- rnorm(n)
      noiseX <- rnorm(n)
      noiseY <- rnorm(n)
      noiseZ <- rnorm(n)
    ## divide data into observational (ExpInd=1) and interventional (ExpInd=2)
      ExpInd <- c(rep(1,n/2),rep(2,n/2))
      noiseX[ which(ExpInd==2)] <- noiseX[ which(ExpInd==2)] * 5
      noiseZ[ which(ExpInd==2)] <- noiseZ[ which(ExpInd==2)] * 3

    ## simulate equilibrium data
      beta <- -0.5
      alpha <- 0.9
      X <- noiseX + 3*W
      Y <- beta* X + noiseY + 3*W
      Z <- alpha*Y + noiseZ
    

      
 ####### Compute "Invariant Causal Prediction" Confidence Intervals
      icp <- hiddenICP(cbind(X,Z),Y,ExpInd,alpha=0.01)

 ###### Print/plot/show summary of output (truth here is (beta,0))
      print(signif(icp$betahat,3))
      print(signif(icp$maximinCoefficients,3))
      cat("true coefficients are:", beta,0)

 #### compare with coefficients from a linear model
      cat("coefficients from linear model:")
      print(signif(coef( lm(Y ~ X + Z -1)),3))

      
  
}

\keyword{Causality }
\keyword{Regression }
\keyword{Hidden Variables }
