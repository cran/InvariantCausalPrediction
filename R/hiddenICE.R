hiddenICE <- function(X, ExpInd, covariance=TRUE, alpha=0.1, threshold =0.75, nsim=100, sampleScenarios=0.5, nodewise=TRUE ){


    X <- as.matrix(X)
    p <- ncol(X)
    
    if(!is.vector(ExpInd)) stop("ExpInd must be a vector")
    if( length(ExpInd) != nrow(X)) stop(" 'ExpInd' must have as many entries as X has rows")
    if(alpha<0) stop("alpha must be non-negative")
    if( threshold <=0.5 | threshold >1) stop("threshold must be between 0.5 and 1")
    if(nsim < 2) stop("'nsim' must be at least 2 (but usually at least 20-100")
    if(sampleScenarios<=0 | sampleScenarios >1) stop(" 'sampleScenarios' must be between 0 and 1")
    if(length((settings <- unique(ExpInd)))<2) stop("need at least three different settings")
    
    q <- if(!nodewise) sqrt(alpha*(2*threshold-1)*(p^2-p)) else sqrt(alpha*(2*threshold-1))
    noS <- min(length(unique(ExpInd)),max(3, ceiling( sampleScenarios * length(settings))))
                   
    
    ## point estimator
    
    Deltalist <- computeDelta(X,  ExpInd, covariance=covariance)$Delta
    Delta <- array(unlist(Deltalist), dim = c(nrow(Deltalist[[1]]), ncol(Deltalist[[1]]), length(Deltalist)))
    estimatedB <- ffdiag(Delta,eps=10^(-4))$B
    rescaledB <- diag(1/diag(estimatedB))%*% estimatedB
    Ahat <- diag(p)- t(rescaledB)
    
    if(alpha>0){

        AhatList <- vector("list", nsim)
        
        for(i in 1:nsim){
            selD <- sort(sample( 1: length(settings), noS))
            eligible <- which( ExpInd %in% settings[selD])
            indices <- sample(eligible, size = length(eligible)/2)
            Xcurrent <- X[indices,]
            ExpIndcurrent <- ExpInd[indices]
            
            ## difference between the covariance matrices of X under the different interventions (Perturb)
            
            Deltalist <- computeDelta(Xcurrent,  ExpIndcurrent, covariance=covariance)$Delta
            Delta <- array(unlist(Deltalist), dim = c(nrow(Deltalist[[1]]), ncol(Deltalist[[1]]), length(Deltalist)))
            estimatedB <- ffdiag(Delta,eps=10^(-4))$B
                                        # B should have ones on the diagonal
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
