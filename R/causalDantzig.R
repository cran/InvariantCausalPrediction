# Dependencies: parallel, linprog, RColorBrewer

# Test: coverage for unregularized causal Dantzig
# Test: many environments for unregularized causal Dantzig
# Test: paths and chosen lambda for unregularized causal Dantzig
# Test: many environments for regularized causal Dantzig
# Test: p=1 for regularized causal Dantzig
# Test: p=1 for unregularized causal Dantzig

# Last changes:
# how to deal with many environments in the unregularized case

## Main function; serves as a wrapper for regularized and unregularized causal Dantzig
#' Causal Dantzig
#'
#' Causal inference in linear structural equation models with hidden variables under additive interventions. This function implements both regularized and unregularized causal Dantzig. For unregularized causal Dantzig, the function returns point estimates and asymptotic confidence intervals. For regularized causal Dantzig, point estimates and a ranking via stability selection are given. 
#' 
#' @param X A nxp matrix with the predictor variables for all experimental settings, where n is the total number of samples and p ist the number of predictors.
#' @param Y The numeric response or target variable of interest of length n. 
#' @param ExpInd Indicator of the experiment or the intervention type an observation belongs to. A discrete vector of the same length as Y with K unique entries if there are K different experiments (for example entry 1 for all observational data and entry 2 for intervention data).
#' @param regularization Use \code{regularization=FALSE} for unregularized causal Dantzig and \code{regulariation=TRUE} for regularized causal Dantzig.
#' @param alpha The significance level of the procedure. This parameter is ignored for regularized causal Dantzig. It defaults to \code{alpha=0.05} for two-sided 95\% confidence intervals.
#' @param lambda User supplied regularization parameter for regularized causal Dantzig. If \code{lambda=NULL}, then the regularization parameter is chosen by cross-validation. The option \code{lambda=NULL} is experimental.
#' @param scale The default \code{scale=TRUE} scales the variables before running regularized causal Dantzig. This parameter is ignored for unregularized causal Dantzig.
#' @param mc.cores Number of CPU cores used for cross-validation. Defaults to \code{mc.cores=1}. This parameter is ignored for unregularized causal Dantzig.
#' @seealso \code{summary} and \code{plot} methods for "causalDantzig" objects. Furthermore, see \code{ICP} in the package \code{InvariantCausalPrediction} for causal inference in a setting without hidden variables.
#' @references Dominik Rothenhausler, Peter Buhlmann, Nicolai Meinshausen (2017):
#' 
#' Causal Dantzig: fast inference in linear structrual equation models with hidden variables under additive interventions
#' 
#' arxiv preprint \url{https://arxiv.org/abs/1706.06159}
#' @export
# #' @importFrom base abs apply array as.character as.matrix as.numeric as.vector by c cat cbind class colMeans colnames do.call exp floor ifelse is.null lapply length list log match.call max mean min names ncol nrow order paste pmax pmin print rbind Reduce rep replicate return rev round rownames rowSums sample sapply simplify2array solve sqrt stop sum t unique unlist warning which which.min
#' @importFrom graphics abline lines plot
#' @importFrom stats pnorm printCoefmat qnorm var
#' @importFrom RColorBrewer brewer.pal
#' @importFrom parallel mclapply
#' @importFrom linprog solveLP
#' @examples
#' set.seed(1)
#' H <- rnorm(2000)
#' X2 <- H+c(rnorm(1000),rnorm(1000,sd=2))
#' Y <- X2 + H + rnorm(2000)
#' X1 <- Y + X2 + c(rnorm(1000),rnorm(1000,sd=2))
#' X3 <- X1 + H + c(rnorm(1000),rnorm(1000,sd=2))
#' ExpInd <- c(rep(0,1000),rep(1,1000))
#' 
#' X <- cbind(X1,X2,X3)
#' 
#' fit <- causalDantzig(X,Y,ExpInd,regularization=FALSE)
#' summary(fit)
#' 
#' fit <- causalDantzig(X,Y,ExpInd,regularization=TRUE)
#' # for parallel processing:
#' # mc.cores <- 2
#' # fit <- causalDantzig(X,Y,ExpInd,regularization=TRUE,mc.cores=mc.cores)
#' summary(fit)
#' plot(fit)
#' 




causalDantzig <- function(X,Y,ExpInd, regularization=FALSE,alpha=0.05,lambda=NULL,scale=TRUE, mc.cores=1) {
    
    X <- as.matrix(X)
    ExpInd <- as.character(ExpInd)
                                        ## warnings
    if(is.null(lambda) && regularization) message("Choosing lambda by cross-validation can be slow. Specifying lambda alleviates this issue. Consult the documentation for more information.")
    if(mc.cores <= 1 && regularization) message("Choosing lambda by cross-validation can be slow. Parallel processing can alleviate this issue. Consult the documentation for more information.")
    if( !(length(ExpInd) == length(Y))) stop("Y and ExpInd do not have equal length.")
    if( !(length(Y) == nrow(X))) stop("The number of rows of X must equal the length of Y.")
    if (length(unique(ExpInd)) <= 1 ) {  stop("ExpInd must have at least two levels.")}

                                        ## normalize to overall mean zero
    foo <- by(data=cbind(X,Y),INDICES=ExpInd,FUN=colMeans)
    overall_mean <- colMeans(do.call(rbind,foo))
    X <- X-overall_mean[1:ncol(X)]
    Y <- Y-overall_mean[ncol(X)+1]
    X <- as.matrix(X)


                                        ## regularized causal Dantzig
    if (regularization==TRUE) {
    if (ncol(X) <=1) stop("The number of columns of X must exceed 1")        
                                        ## prepare scaling
        if (scale) {
                                        ## split of one against all
                                        ## weights <- lapply(unique(ExpInd),FUN= function(e) sqrt(diag(t(X[ExpInd==e,]) %*% X[ExpInd==e,]*(1/sum(ExpInd==e)-1/length(ExpInd))^2+(1/length(ExpInd))^2*diag(t(X[ExpInd!=e,]) %*% X[ExpInd!=e,]))))

                                        ## split of one against the rest
            foo <- by(X,INDICES=ExpInd,FUN= function(x) colMeans(x^2))
            weights <- lapply(unique(ExpInd), FUN= function(e) sqrt(foo[[e]]+(Reduce("+",foo)-foo[[e]])/(length(foo)-1)^2))
            names(weights) <- unique(ExpInd)
        } else {
            weights <- lapply(rep(1,length(unique(ExpInd))),FUN=rep,times=ncol(X))
            names(weights) <- unique(ExpInd)
        }

                                        ## Prepare lambdasequence
        
        if(is.null(lambda)) {

                                        ## Define a lambdasequence if lambda is not given
            XxY <- X * Y
            foo <- by(XxY, INDICES=ExpInd,FUN=function(x) colMeans(as.matrix(x)))
            cov_XY_diff_abs_weighted <- lapply(unique(ExpInd), FUN= function(e) abs(foo[[e]]-(Reduce("+",foo)-foo[[e]])/(length(foo)-1))/weights[[e]])
            maxi <- do.call(max,cov_XY_diff_abs_weighted)
            mini <- maxi*0.001
            lambdasequence <- mini*exp(0:29*log(maxi/mini)/29)
                                        ##        some choices of lambdasequence can lead to numerical instabilities in linprog
                                        ##        lambdasequence <- 0.002*exp(0:99*log(max/0.002)/99)
            
                                        ## compute loss of cross-validation
            if (mc.cores < 2) {
                crossv <- simplify2array(lapply(lambdasequence,compute_loss_cDantzig,X,Y,ExpInd,weights=weights))
            } else {
                crossv <- simplify2array(parallel::mclapply(lambdasequence,compute_loss_cDantzig,X,Y,ExpInd,weights=weights,mc.cores =mc.cores))
            }
            lambda_selected<- lambdasequence[which.min(crossv)]
            names(crossv) <- lambdasequence
            # some of the crossv might be infinite; indicating that linprog returned numeric(0)

                                        ## compute causaldantzig path on the sequence of lambdas
            if (mc.cores < 2) {
                pathlapply <- lapply(lambdasequence,cDantzig,X,Y,ExpInd,weights=weights)
            } else {
                pathlapply <- parallel::mclapply(lambdasequence,cDantzig,X,Y,ExpInd,weights=weights,mc.cores = mc.cores)
            }
            names(pathlapply) <- lambdasequence
            ## Make sure that if cDantzig returns numeric(0) that the corresponding lambda is dropped
            ## both from lambdasequence and the path
            for (i in names(pathlapply)) if (length(pathlapply[[i]])==0){ pathlapply[[i]] <- NULL}
            lambdasequence <- as.numeric(names(pathlapply))
            
            
            path <- as.matrix(t(simplify2array(pathlapply)))
        } else {
            lambda_selected <- lambda
            path <- NULL
            crossv <- NULL
            lambdasequence <- NULL
            
        }

        coefficients <- cDantzig(X,Y,ExpInd,weights=weights,lambda=lambda_selected)

                                        ## stability selection
        stability_selection <- stability_selection(lambda_selected,X,Y,ExpInd,weights=weights) 
        names(coefficients) <- colnames(X)
        
        fit <- list(coefficients=coefficients, path=path, lambdasequence = lambdasequence, crossv=crossv,stability_selection=stability_selection,call=match.call(),regularization=TRUE,lambda_selected = lambda_selected,weights=weights)
        class(fit) <- "causalDantzig"
        return(fit)
    }
    else {
        if (length(unique(ExpInd))>20) warning("Too many environments may lead to long runtimes and very conservative p-values and confidence intervals. Pooling environments can alleviate this issue.")
        if (scale != TRUE) warning("Option 'scale' ignored for unregularized causalDantzig (causalDantzig(...,regularization=FALSE))")
        if (!is.null(lambda)) warning("Option 'lambda' ignored for unregularized causalDantzig (causalDantzig(...,regularization=FALSE))")

        CI <- t(replicate(ncol(X),c(-Inf,Inf)))
        p_value <- rep(1,ncol(X))
        coefficients <- rep(0,ncol(X))
        coefmax <- rep(0,ncol(X))
        coefmin <- rep(0,ncol(X))
          
        for (l in unique(ExpInd)) { 
          Ep <- ifelse(ExpInd == l,0,1) 
          fit_unregularized <- unreg_cDantzig(X,Y,Ep,mult_adjust = ifelse(length(unique(ExpInd))==2,1,length(unique(ExpInd))),alpha=alpha)
          CI[,1] <- pmax(fit_unregularized$CI[,1],CI[,1])
          CI[,2] <- pmin(fit_unregularized$CI[,2],CI[,2])
          for (j in 1:nrow(CI)){
          if (CI[j,1]==fit_unregularized$CI[j,1]) {coefmin <- fit_unregularized$coefficients}
          if (CI[j,2]==fit_unregularized$CI[j,2]) {coefmax <- fit_unregularized$coefficients}
          }      
          p_value <- pmin(p_value,fit_unregularized$p_value)
        }
        ## compute hat beta
        # this is problematic: sometimes CIs are bonkers but coefficients are fine
        #coefficients <- (CI[,2]+CI[,1])/2
        # this is still problematic... wtf
        coefficients <- (coefmax + coefmin)/2
        names(coefficients) <- colnames(X)
        rownames(CI) <- colnames(X)
        colnames(CI) <- c("lower bound", "upper bound")
        names(p_value) <- colnames(X)
                                        ## prepare return object
        fit <- list(coefficients= coefficients, CI = CI, p_value = p_value,alpha=alpha, call=match.call(), regularization = FALSE)
        class(fit) <- "causalDantzig"
        return(fit)
    }
}


#' @keywords internal
                                        ## unregularized causal Dantzig - called within function causalDantzig 
unreg_cDantzig <- function(X,Y,ExpInd,mult_adjust=1,alpha=0.05) {
    X1 <-  as.matrix(X[ExpInd==1,])
    Y1 <-  as.matrix(Y[ExpInd==1])
    X0 <- as.matrix(X[ExpInd==0,])
    Y0 <- as.matrix(Y[ExpInd==0])
    if (nrow(X1) < ncol(X1) | nrow(X0) < ncol(X0)) stop("System is singular. Unregularized causalDantzig does not have a unique solution. Call causalDantzig(...,regularization=TRUE) for regularized causalDantzig.")
                                        ## Compute hat G^{-1} and hat Z
    Ginv <- solve(t(X1) %*% X1 / nrow(X1) - t(X0) %*% X0 / nrow(X0))
    Z <-  t(X1) %*% Y1 / nrow(X1) - t(X0) %*% Y0 / nrow(X0)
                                        ## compute the estimator
    causalDantzig <-   Ginv %*% Z
    
    
                                        ## Compute samples for approximate variance of betahat
    samples1 <- lapply(1:nrow(X1),function(i)  -Ginv %*% X1[i,] %*% X1[i,] %*% Ginv %*% Z +Ginv  %*% X1[i,] * Y1[i])
    samples0 <- lapply(1:nrow(X0),function(i)  -Ginv %*% X0[i,] %*% X0[i,] %*% Ginv %*% Z +Ginv  %*% X0[i,] * Y0[i])
    list2array <- function(list) array(unlist(list), dim = c(nrow(list[[1]]), ncol(list[[1]]), length(list)))
    var <- rep(0,length(causalDantzig))
                                        ## compute asymptotic variance for each component of betahat
    for (i in 1:length(causalDantzig)) {
        samples1transf <- list2array(samples1)[i,1,]
        samples0transf <- list2array(samples0)[i,1,]
        var[i] <- var(samples1transf)/length(samples1transf)+var(samples0transf)/length(samples0transf)
    }

    CI <- cbind(causalDantzig-qnorm(1-0.5*alpha/mult_adjust)*sqrt(var),causalDantzig+qnorm(1-0.5*alpha/mult_adjust)*sqrt(var))
    colnames(CI) <- c("lower bound","upper bound")
    p_value <- round( 2*(1-pnorm(abs(causalDantzig/(mult_adjust*sqrt(var))))),3)

    return(list(coefficients = causalDantzig, CI=CI,sd=sqrt(var),p_value=p_value))
}

#' @keywords internal
                                        ## regularized causal Dantzig - called within function causalDantzig and compute_loss_cDantzig
cDantzig <- function (lambda,X,Y,ExpInd,weights=NULL) {
    ExpInd <- as.character(ExpInd)

    foo <- by(data=as.matrix(X),INDICES=ExpInd,FUN= function(x) t(as.matrix(x)) %*% as.matrix(x) / nrow(as.matrix(x)))
    cov_XX_diff <- lapply(unique(ExpInd), FUN= function(e) foo[[e]]-(Reduce("+",foo)-foo[[e]])/(length(foo)-1))
                                        ### one environment vs unweighted rest
                                        ## cov_XX_diff<- function(e) {return(t(X[ExpInd==e,]) %*% X[ExpInd==e,]/nrow(X[ExpInd==e,]) -  t(X) %*% X / nrow(X)  )}
    diffX <- do.call(rbind,cov_XX_diff)

                                        ##cov_XY_diff<- function(e) {return(t(X[ExpInd==e,]) %*% Y[ExpInd==e] / length(Y[ExpInd==e]) - t(X) %*% Y / nrow(X))}
    XxY <- X * Y
    foo <- by(XxY, INDICES=ExpInd,FUN=function(x) colMeans(as.matrix(x)))
    cov_XY_diff<- lapply(unique(ExpInd), FUN= function(e) foo[[e]]-(Reduce("+",foo)-foo[[e]])/(length(foo)-1))
    diffY <- do.call("c",cov_XY_diff)

    Atop <- cbind( - diffX , diffX )
    A <- rbind(Atop, -Atop)
    btop <- - diffY
    rownames(A) <- NULL

    if (is.null(weights)) {
      b <- c(btop+ rep(lambda,length(btop)),-btop+rep(lambda,length(btop)))
    } else {
      expanded_weights <- do.call('c',weights)
      weightsclipped <- sapply(expanded_weights,max,0.1)
      b <- c(btop+ lambda*weightsclipped,-btop+lambda*weightsclipped)
    }   

    rownames(b) <- NULL
    c <- rep(1,2*ncol(X))
    
    LPsolution <- linprog::solveLP(c,b,A,lpSolve=TRUE)
    coefficients <- LPsolution$solution[1:ncol(X)] - LPsolution$solution[(ncol(X)+1):(2*ncol(X))]
    
    return(coefficients)
}

#######################################################################
### functions for stability selection and crossvalidation
#######################################################################

#' @keywords internal
                                        ## Performs stability selection for regularized causal Dantzig - called within function causalDantzig
stability_selection <- function(lambda,X,Y,ExpInd,weights=NULL,n=30) {
    sample_and_run <- function() {
        sample_indices <- do.call("c",lapply(unique(ExpInd), function(e) sample(x=which(ExpInd==e),size=sum(ExpInd==e),replace=TRUE) ))
        Xp <- X[sample_indices,]
        Yp <- Y[sample_indices]
        Ep <- ExpInd[sample_indices]
        coefficients <- cDantzig(lambda=lambda,X=Xp,Y=Yp,ExpInd=Ep,weights=weights)
        return(coefficients)
    }
    m <- replicate(n=n,sample_and_run())
    ## get rid of numeric(0) entries
    i <-1
    while (i <= length(m)) { 
      if (length(m[[i]])== 0) { m[[i]] <- NULL }
      else { i <- i+1}
    }
    m <- simplify2array(m)
    
    
    active_set <- abs(m) >= 10e-10
    probability <- apply(active_set, 1,FUN=mean)
    names(probability) <- colnames(X)
    average_strength <- apply(m,1,FUN=mean)
    names(average_strength) <- colnames(X)
    
    return(list(probability=probability,average_strength=average_strength))
}


## Estimates loss of causalDantzig for a given lambda via crossvalidation
## Called within function causalDantzig()

#' @keywords internal
compute_loss_cDantzig <- function(lambda,X,Y,ExpInd,weights=NULL) {
    ExpInd <- as.character(ExpInd)

    folds <- 10
    residuals <- rep(0,folds)
    i <- 0
    for (i in 0:(folds-1)) {
        select <- rep(TRUE,nrow(X))
        for (e in unique(ExpInd)) {
            samples <- sum(ExpInd==e)
            selectionregion <- (floor(i*samples/folds)+1):floor((i+1)*samples/folds)

            select[ExpInd==e][selectionregion] <- FALSE
        }

        gamma <- cDantzig(X[select,],Y[select],ExpInd[select],lambda=lambda,weights=weights)
        ## gamma <- cDantzig(X[select,],Y[select],ExpInd[select],lambda=lambda,weights=NULL)

        Xp <- X[!select,]
        Yp <- Y[!select]
        Ep <- ExpInd[!select]

        if(is.null(weights)){
          weights <- lapply(rep(1,length(unique(ExpInd))),FUN=rep,times=ncol(X))
          names(weights) <- unique(ExpInd)
        }

        if (length(gamma) != 0) {

        XxY <- Xp * as.vector(Yp- Xp %*% gamma) 
        foo <- by(XxY, INDICES=Ep,FUN=function(x) colMeans(as.matrix(x)))

        cov_XY_diff_abs_weighted <- lapply(unique(Ep), FUN= function(e) abs(foo[[e]]-(Reduce("+",foo)-foo[[e]])/(length(foo)-1))/weights[[e]])
        ### without weights
        ## cov_XY_diff_abs_weighted <- lapply(unique(Ep), FUN= function(e) abs(foo[[e]]-(Reduce("+",foo)-foo[[e]])/(length(foo)-1)))
        maxi <- do.call(max,cov_XY_diff_abs_weighted)
                
        residuals[i+1] <- maxi
        
        } else { residuals[i+1] <- Inf}

    }
    return(mean(residuals))
}

#######################################################################
#### Plot, print and summary functions for causalDantzig objects
#######################################################################

#' Plot coefficients paths from a "causalDantzig" object with regularization
#'
#' Gives a plot of the coefficient paths for a fitted "causalDantzig" object with regularization.
#'
#' @param x fitted "causalDantzig" model with regularization
#' @param ... Additional inputs to generic plot function (not used).
#' @seealso \code{causalDantzig} and \code{summary} method.
#' @export



plot.causalDantzig <-  function(x,...) {
    
    if(x$regularization == FALSE | is.null(x$lambdasequence)) stop("Object does not contain the solution path of causalDantzig. Call causalDantzig(...,lambda=NULL, regularization=TRUE) to compute the solution path of causalDantzig with regularization")
    
    cols <- RColorBrewer::brewer.pal(9,"Set1")
    cols <- c(cols, RColorBrewer::brewer.pal(8,"Set2"))
    cols <- c(cols, RColorBrewer::brewer.pal(12,"Set3"))
    
    path <- x$path
    crossv <- x$crossv
    
    # path may have deleted different lambdas as opposed to crossv (i.e. linprog returned numeric(0))
    lambda_index <- which(rownames(path)==names(which.min(crossv)))
    
                                        ## reverse path and normalize
    normalized_path <- rev(rowSums(abs(path))/max(rowSums(abs(path))))
    argmin <- normalized_path[length(normalized_path)+1-lambda_index]
    
                                        ## stabilize coefficients path - sometimes values with higher lambda have a smaller l1 norm    
    partialmax <- function(i) {return(max(normalized_path[1:i]))}
    normalized_path<-sapply(1:length(normalized_path),partialmax)
    
    plot(normalized_path,normalized_path,type="n", xlim=c(min(normalized_path),max(normalized_path)),ylim=c(min(x$path),max(x$path)),xlab="|beta|/|beta_max|",ylab="coefficients")
    for (i in 1:length(x$coefficients)) {
        lines(normalized_path,rev(path[,i]),col=cols[(i %% length(cols))+1],lwd=2.2)
    }
    abline(v=argmin)
}

#' Summarizes "causalDantzig" fits
#' 
#' \code{summary} method for class "causalDantzig".
#' 
#' @param object fitted "causalDantzig" model with or without regularization
#' @param ... Additional inputs to generic summary function (not used).
#' @seealso \code{causalDantzig} and \code{plot} method.
#' @export


summary.causalDantzig <- function(object,...) {
    
    if(object$regularization  == TRUE) {
        cat("Regularized causal Dantzig \n")
        if (!is.null(object$lambdasequence)) {cat("Regularization parameter lambda is chosen by crossvalidation\n \n")}
        else {cat("\n")}
        cat("Call:\n")
        print(object$call)
        cat("\n")
        cat("Stability selection:\n")
        cat("Variables with largest probability of being in the active set under subsampling\n \n")
        if(is.null(names(object$coefficients))) {names(object$coefficients) <- paste("X",1:length(object$coefficients),sep="")}
        TAB <- do.call(cbind, object$stability_selection)
        rownames(TAB) <- names(object$coefficients)
        colnames(TAB) <- c("probability", "average_strength")
        TAB <- TAB[order(-TAB[,1]),][1:min(nrow(TAB),10),]
        TAB <- round(TAB,3)
        printCoefmat(TAB)
        cat("\n")
    } else {
        TAB <- cbind(Estimate = round(object$coefficients,3),
                     "Lower bound" = round(object$CI[,1],3),
                     "Upper bound" = round(object$CI[,2],3),
                     p.value = object$p_value)
        if(is.null(names(object$coefficients))) {rownames(TAB) <- paste("X",1:length(object$coefficients),sep="")}
        colnames(TAB) <- c("Estimate","Lower bound", "Upper bound", "p.value")
        cat("Unregularized causal Dantzig at level",object$alpha, "\n")
        cat("Call:\n")
        print(object$call)
        cat("\n")
        printCoefmat(TAB,has.Pvalue = TRUE)
    }
    
}


#' @export

print.causalDantzig <- function(x,...) summary.causalDantzig(x)
