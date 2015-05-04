hiddenICE <- function(X, ExpInd, covariance=TRUE, alpha=0.1, threshold =0.75, nsim=100, sampleSettings=1/sqrt(2), sampleObservations=1/sqrt(2), nodewise=TRUE ){

    if(!is.matrix(X) & !is.data.frame(X)) stop("'X' must be a matrix of data frame")

    X <- as.matrix(X)
    p <- ncol(X)
    
    if(!is.vector(ExpInd)) stop("ExpInd must be a vector")
    if( length(ExpInd) != nrow(X)) stop(" 'ExpInd' must have as many entries as X has rows")
    if(alpha<0) stop("alpha must be non-negative")
    if( threshold <=0.5 | threshold >1) stop("threshold must be between 0.5 and 1")
    if(nsim < 2) stop("'nsim' must be at least 2 (but usually at least 20-100")
    if(length((settings <- unique(ExpInd)))<2) stop("need at least three different settings")
     if(sampleSettings<=0) stop("sampleSettings needs to be positive")
    if(sampleObservations<=0) stop("sampleObservations needs to be positive")
    if(sampleSettings>1) stop("sampleSettings needs to be at most 1")
    if(sampleObservations>1) stop("sampleObservations needs to be at most 1")
    q <- if(!nodewise) sqrt(alpha*(2*threshold-1)*(p^2-p)) else sqrt(alpha*(2*threshold-1))
    subs <- sampleSettings* length(settings)
    drawE <- function(x){
        z <- floor(x)
        if( runif(1) <  x-z) z <- z+1
        z <- max(3,z)
        return(z)
    }        
    
    ## point estimator
    
    Deltalist <- computeDelta(X,  ExpInd, covariance=covariance)$Delta
    Delta <- array(unlist(Deltalist), dim = c(nrow(Deltalist[[1]]), ncol(Deltalist[[1]]), length(Deltalist)))
    estimatedB <- ffdiag(Delta,eps=10^(-4))$B
    rescaledB <- diag(1/diag(estimatedB))%*% estimatedB
    Ahat <- diag(p)- t(rescaledB)
    
    if(alpha>0){

        AhatList <- vector("list", nsim)
        
        for(i in 1:nsim){
            useSettings <- sample( settings, drawE(subs))
            ind <- which(  ExpInd %in% useSettings)
            useSamples <- sort(sample(ind, round(length(ind)*sampleObservations)))
        
            Xcurrent <- X[useSamples,]
            ExpIndcurrent <- ExpInd[useSamples]
            
            ## difference between the covariance matrices of X under the different interventions (Perturb)
            
            Deltalist <- computeDelta(Xcurrent,  ExpIndcurrent, covariance=covariance)$Delta
            Delta <- array(unlist(Deltalist), dim = c(nrow(Deltalist[[1]]), ncol(Deltalist[[1]]), length(Deltalist)))
            estimatedB <- ffdiag(Delta,eps=10^(-4))$B
            rescaledB <- diag(1/diag(estimatedB))%*% estimatedB
            AhatJD <- diag(p)- t(rescaledB)
            
            diag(AhatJD) <- 0
            AhatList[[i]] <- edgeSelection(AhatJD, q,  nodewise=nodewise)
        }
        
        AhatAdjacency <- edgeRetention(AhatList, threshold, p)
    }else{
        AhatAdjacency <- NULL
    }
    list(Ahat=Ahat, AhatAdjacency = AhatAdjacency)
}
