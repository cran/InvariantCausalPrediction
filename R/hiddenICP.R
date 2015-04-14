
hiddenICP <-
function(X,Y,ExpInd,alpha=0.1,mode="asymptotic",intercept=TRUE){

    if(nrow(X) <= ncol(X)) stop(paste( "hiddenICP not suitable for high-dimensional data (at the moment) \n -- need nrow(X)>ncol(X) but have nrow(X)=",nrow(X)," and ncol(X)=",ncol(X)))

    if(intercept){
        X <- cbind(rep(1,nrow(X)),X)
    }
    
    if(!is.list(ExpInd)){
        if(length(ExpInd)!=length(Y)) stop("if `ExpInd' is a vector, it needs to have the same length as `Y'")
        uni <- unique(ExpInd)
        K <- length(uni)
        ExpIndNEW <- list()
        for (uc in 1:K){
            ExpIndNEW[[uc]] <- which(ExpInd==uni[uc])
            attr(ExpIndNEW[[uc]],"value") <- uni[uc]
        }
        ExpInd <- ExpIndNEW
        rm(ExpIndNEW)
    }else{
        ran <- range(unlist(ExpInd))
        if(ran[1]<1) stop(paste("if `ExpInd' is a list with indicies of observations, \n minimal entry has to be at least 1 but is",ran[1]))
        if(ran[2]>length(Y)) stop(paste("if `ExpInd' is a list with indicies of observations, \n maximal entry has to be at most equal \n to the length",length(Y),"of the observations but is",ran[2]))
    }
    X <- as.matrix(X)
    K <- length(ExpInd)
    p <- ncol(X)
    n <- nrow(X)

    kc <- 1
    KC <- if(K>2) K else 1
    for (kc in 1: KC){
        ins <- ExpInd[[kc]]
        out <- (1:n)[-ins]
        DS <- (t(X[ins,])%*%X[ins,])/length(ins) - (t(X[out,])%*%X[out,])/length(out)
        Drho <- (t(X[ins,])%*%Y[ins])/length(ins) - (t(X[out,])%*%Y[out])/length(out)
        DSI <- solve(DS)
    
        
        betahat <- as.numeric(solve( DS,Drho))
        if(kc==1) betahatall <- betahat else betahatall <- betahatall+betahat
        Zin <- matrix(nrow=length(ins),ncol=p)
        Zout <- matrix(nrow=length(out),ncol=p)
        for (i in 1:length(ins)){
            tmp <- DSI %*% X[ins[i],]
            Zin[i,] <- as.numeric(- tmp * sum(tmp*Drho) + Y[ins[i]] *tmp)
        }
        for (i in 1:length(out)){
            tmp <- DSI %*% X[out[i],]
            Zout[i,] <- as.numeric(- tmp * sum(tmp*Drho) + Y[out[i]] *tmp)
        }
        sigmavec <- sqrt(diag((cov(Zin)/length(ins)+cov(Zout)/length(out))))
        
        maximineffectsN <- sign(betahat) * pmax( 0, abs(betahat) - qnorm(max(0.5,1-alpha/(p*K))) * sigmavec)
        if(kc==1){
            maximineffects <- maximineffectsN
        }else{
            for (varc in 1:p){
                if(abs(maximineffectsN[varc]) > abs(maximineffects[varc])) maximineffects[varc] <- maximineffectsN[varc]
            }
        }
    }
    betahat <- betahatall/KC
    maximinCoefficients <- maximineffects
    if(intercept){
        betahat <- betahat[-1]
        maximinCoefficients <- maximinCoefficients[-1]
    }
    retobj <- list(betahat=betahat,maximinCoefficients=maximinCoefficients)
    class(retobj) <- "hiddenInvariantCausalPrediction"
    return(retobj)
                       
}



