\name{ICP}
\alias{ICP}
\title{
Invariant Causal Prediction
}
\description{
Confidence intervals for causal effects in a regression setting.
}
\usage{
ICP(X, Y, ExpInd, alpha=0.01, test = "normal",
              selection = c("lasso", "all", "stability", "boosting")
                                      [if (ncol(X) <=8) 2 else 4],
              maxNoVariables = 8, maxNoVariablesSimult = 8,
              maxNoObs = 200,
              showAcceptedSets = TRUE,
              showCompletion = TRUE,
              stopIfEmpty = FALSE,
              gof = max(0.01,alpha))
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
  in the data (for example two elements: first element is a vector that contains indices of observations that are observational data and 
  second element is a vector that contains indices 
  of all observations that are of interventional type).
}
\item{alpha}{
     Determines coverage of the confidence regions. \code{alpha} is the probability of non-coverage, so the default of 0.01 produces 99\% confidence intervals. If no model has a p-value of at least \code{alpha}, the model assumptions seem violated and \code{modelReject} will return \code{TRUE}.
}
  \item{test}{
  Use "exact" for an exact test in a regression setting, especially if sample size is small. However, this test is computationally demanding if sample size is high. The default "normal" tests for a shift in mean and variance of the residuals between different populations/environments. Using "correlation" tests additionally for vanishing correlation between predictor variables and residuals in each environment. Other options are "ranks" that uses a rank-based alternative and "ks" for a Kolmogorov-Smirnov test to detect differences in the distributions of the residuals between different environments. It is also possible to supply a function of the form \code{function(x,z)}  that takes samples \code{x} and \code{z} of two populations as arguments and whose return value is the p-value for the null hypothesis that the two underlying distribution are identical.
}
  \item{selection}{
  The method for pre-selection of variables to save computational resources. 
  Can use "all" for no pre-selection (which guarantees coverage but might take longer to compute), 
  "boosting" for a boosting-type, "lasso" for Lasso cross-validated or 
  "stability" for a stability-selection-type pre-selection.
   Default is "all" if p does not exceed 10 and "boosting" otherwise.    
}
\item{maxNoVariables}{
  The maximal number of variables to pre-select 
  (choosing smaller values saves computational resources 
  but increases approximation error).
}
  \item{maxNoVariablesSimult}{
  The maximal size of sets of variables considered in the procedure 
  (same comment as for \code{maxNoVariables}).
}
  \item{maxNoObs}{
  The maximal number of observations used for the "exact" test 
  (same comment as for \code{maxNoVariables}).
}
\item{showAcceptedSets}{
  If \code{TRUE}, print out information about accepted sets of variables.
}
  \item{showCompletion}{
  If \code{TRUE}, print out information about progress of computation.
}
 \item{stopIfEmpty}{
  If \code{TRUE}, the procedure will stop computing confidence intervals
  if the empty set has been accepted (and hence no variable can have a
  signicificant causal effect). Setting to \code{TRUE} will save
  computational time in these cases, but means that the confidence intervals lose their
  coverage properties for values different to 0.
}
 \item{gof}{
 If no set of variables (including the empty set) leads to a p-value larger than the goodness-of-fit cutoff \code{gof}, 
the whole model will be rejected. If the model is correct, this will happen with a probability of \code{gof} 
and this option protects again making statements when the model is obviously not suitable for the data.
}
}

\value{
A list with elements
  \item{ConfInt}{The matrix with confidence intervals for the causal
    coefficient of all variables. 
    First row is the upper bound and second row the lower bound.}
  \item{maximinCoefficients}{The value in the confidence interval
    closest to 0. Is hence non-zero for variables with significant effects.}
  \item{colnames}{The names of the variables (replaced with generic "Variable 1" etc. if not available).}
   \item{factor}{Logical indicating whether the response is a factor or not.}
   \item{dimX}{The dimensions of the matrix with predictor variables.}
   \item{Coeff}{A list which contains for all variables the vector with point-estimates among all accepted sets where the variable was part of the set.}
  \item{CoeffVar}{Same as \code{Coeff} but with the standard deviation of the point-estimate.}
  \item{modelReject}{Logical indicating if the whole model was rejected (the p-value of the best fitting model is too low).}
  \item{acceptedSets}{A list with one element per accepted set. For each accepted model the list entry is a vector that contains the indices of the variables in the accepted set.}
  \item{usedvariables}{The pre-selected variables if not all variables are used for the analysis; otherwise all variables.}
  \item{pvalues}{The p-values of all variables.}
  \item{stopIfEmpty}{A boolean value indicating whether computations stop as soon as intersection of accepted sets is empty.}
  \item{noEnv}{The number of distinct environments.}
  \item{gof}{The goodness-of-fit cutoff \code{gof} used.}
  \item{bestModel}{The largest p-value across all tested sets of variables.}
}
\references{
 Jonas Peters, Peter Buhlmann, Nicolai Meinshausen (2015):

 Causal inference using invariant prediction: identification and confidence intervals

 arxiv preprint \url{http://arxiv.org/abs/1501.01332}
}
\author{
Nicolai Meinshausen <meinshausen@stat.math.ethz.ch>
}
\seealso{  \code{\link{hiddenICP}} for reconstructing the parents of a variable in the presence of  hidden variables (but assuming shift/additive interventions), which also allows construction of confidence intervals for the causal coefficients. See package "backShift" for constructing point estimates of causal cyclic models in the presence of hidden variables (again under shift interventions). 
}
\examples{
 ##########################################
 ####### 1st example:
 ####### Simulate data with interventions
set.seed(1)
    ## sample size n
      n <- 4000
    ## 5 predictor variables
      p  <- 5
    ## simulate as independent Gaussian variables
      X <- matrix(rnorm(n*p),nrow=n)
    ## divide data into observational (ExpInd=1) and interventional (ExpInd=2)
      ExpInd <- c(rep(1,n/2),rep(2,n/2))
    ## for interventional data (ExpInd==2): change distribution
      X[ExpInd==2,] <- sweep(X[ExpInd==2,],2, 5*rnorm(p) ,FUN="*")
    ## first two variables are the causal predictors of Y
      beta <- c(1,1,rep(0,p-2))
    ## response variable Y
      Y <- as.numeric(X\%*\%beta + rnorm(n))
   ## optinal: make last variable a child of Y (so last variable is non-causal for Y)
      X[,p] <- 0.3*Y + rnorm(n)

 ####### Compute "Invariant Causal Prediction" Confidence Intervals
      icp <- ICP(X,Y,ExpInd)

 ###### Print/plot/show summary of output
      print(icp)
      plot(icp)

#### compare with linear model 
      cat("\n compare with linear model  \n")            
      print(summary(lm(Y~X)))


##########################################
 ####### 2nd example:
 ####### Simulate a DAG where X1 -> Y, Y -> X2 and Y -> X3
 ####### noise interventions on second half of data on X1
 ####### structure of DAG (at Y -> X2) is changing under interventions
     n1 <- 400
     n2 <- 500
     ExpInd <- c(rep(1,n1), rep(2,n2))
   ## index for observational (ExpInd=1) and intervention data (ExpInd=2)
     X1 <- c(rnorm(n1), 2 * rnorm(n2) + 1)
     Y <- 0.5 * X1 + 0.2 * rnorm(n1 + n2)
     X2 <- c(1.5 * Y[1:n1]      + 0.4 * rnorm(n1),
            -0.3 * Y[(n1+1):n2] + 0.4 * rnorm(n2))
     X3 <-  -0.4 * Y            + 0.2 * rnorm(n1 + n2)
     X <- cbind(X1, X2, X3)

 ### Compute "Invariant Causal Prediction" Confidence Intervals
  ## use a rank-based test to detect shift in distribution of residuals
     icp <- ICP(X, Y, ExpInd,test="ranks")
  ## use a Kolmogorov-Smirnov test to detect shift in distribution of residuals
     icp <- ICP(X, Y, ExpInd,test="ks")
  ## can also supply test as a function
  ## here chosen to be equivalent to option "ks" above
     icp <- ICP(X, Y, ExpInd,test=function(x,z) ks.test(x,z)$p.val)
 ## use a test based on normal approximations
     icp <- ICP(X, Y, ExpInd, test="normal")


 ### Print/plot/show summary of output
     print(icp)
     plot(icp)

#### compare with linear model 
     cat("\n compare with linear model \n")            
     print(summary(lm(Y~X)))




\dontrun{
 ##########################################
 ####### 3rd example:
 ####### College Distance data 
     library(AER)
     data("CollegeDistance")
     CD <- CollegeDistance

  ##  define two experimental settings by
  ##  distance to closest 4-year college
     ExpInd <- list()
     ExpInd[[1]] <- which(CD$distance < quantile(CD$distance,0.5))
     ExpInd[[2]] <- which(CD$distance >= quantile(CD$distance,0.5))

  ## target variable is binary (did education lead at least to BA degree?)
     Y <- as.factor(CD$education>=16)
  ## use these predictors
     X <- CD[,c("gender","ethnicity","score","fcollege","mcollege","home",
        "urban","unemp","wage","tuition","income","region")]

  ## searching all subsets (use selection="lasso" or selection="stability"
  ##     to select a subset of subsets to search)
  ## with selection="all" the function will take several minutes
    icp <- ICP(X,Y,ExpInd,selection="all",alpha=0.1)

  ## Print/plot/show summary of output
     print(icp)
     summary(icp)
     plot(icp)

}


}

\keyword{Causality }
\keyword{Regression }
\keyword{Classification }
