ICP <-
function(X,Y,ExpInd,alpha=0.1, test="approximate", selection = c("lasso","all","stability","boosting")[if(ncol(X)<=10) 2 else 4], maxNoVariables=10, maxNoVariablesSimult=10, maxNoObs=200, showAcceptedSets=TRUE, showCompletion=TRUE, stopIfEmpty = FALSE){
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

    getblanket <- getblanketall
    if(selection=="lasso"){
        getblanket <- getblanketlasso
    }
    if(selection=="stability"){
        getblanket <- getblanketstability
    }
    if(selection=="boosting"){
        getblanket <- getblanketboosting
    }

    if(is.data.frame(X)){
        if(any(sapply(X,class)=="factor")){
            Z <- X
            X <- matrix(nrow=nrow(Z),ncol=sum(sapply(X,function(x) if(is.factor(x)) length(levels(x)) else 1)))
            cc <- 0
            colX <- character(0)
            for (k in 1:ncol(Z)){
                if(is.numeric(Z[,k])){
                    cc <- cc+1
                    X[,cc] <- Z[,k]
                    colX <- c(colX, colnames(Z)[k])
                }else{
                    nf <- length(lev <- levels(Z[,k]))
                    for (nfc in 1:(nf-1)){
                        cc <- cc+1
                        X[,cc] <- as.numeric(Z[,k]==lev[nfc])
                        colX <- c(colX, paste(colnames(Z)[k],"_",lev[nfc],sep=""))
                    }
                }
            }
            X <- X[,1:cc]
            colnames(X) <- colX
            X <- as.matrix(X)
        }
    }

    if( length(unique(Y))==2 & !is.factor(Y)){
        warning("\n Y only has 2 unique values -- using classification")
        Y <- as.factor(Y)
    }
            
    K <- length(ExpInd)
    
    n <- nrow(X)
    p <- ncol(X)
    X <- cbind(rep(1,nrow(X)),X)
    
    Clist <- matrix(NA,nrow=2,ncol=p)
    Coeff <- list()
    CoeffVar <- list()
    for (k in 1:p){
        Coeff[[k]] <- numeric(0)
        CoeffVar[[k]] <- numeric(0)
    }
    pvall <- numeric(0)
    
    cont <- TRUE
    pvalempty <- getpval(Y,X[,1,drop=FALSE],ExpInd,test=test,maxNoObs=maxNoObs)$pval
    
    
    if(pvalempty > alpha){
        
        Clist <- matrix(0,nrow=2,ncol=p)
        pvall <- c(pvall, pvalempty)
        if(showAcceptedSets) cat(paste("\n  accepted empty set"))
        if(stopIfEmpty) cont <- FALSE
    }
    testsets <- if( any(!is.null(c(maxNoVariables,maxNoVariablesSimult))))  getblanket(X[,-1,drop=FALSE],Y,maxNoVariables=maxNoVariables,maxNoVariablesSimult=maxNoVariablesSimult) else getblanket(X[,-1,drop=FALSE],Y)
    len <- sapply(testsets,length)
    lcsingle <- sum(len==1)
    usedvariables <- unique(unlist(testsets))

    
    lc <- 0

    printoutat <- 2^(1:ceiling(log2(length(testsets))))
    while(cont && lc < length(testsets)){
        if(showCompletion){
            if( lc %in% printoutat){
                cat(paste("\n *** ",round(100*lc/length(testsets)),"% complete: tested ", lc," of ", length(testsets), " sets of variables ",sep=""))
            }
        }
        lc <- lc+1
        usevariab <- testsets[[lc]]
        notusevariab <- (1:p)[-testsets[[lc]]]
        
        tmp <- getpval(Y,X[, c(1,1+usevariab),drop=FALSE],ExpInd,test=test,maxNoObs=maxNoObs)
        pval <- tmp$pval
        if(pval > alpha){
            if(showAcceptedSets) cat(paste("\n accepted set of variables ", paste(usevariab,collapse=","),sep=""))
            Clist[1,usevariab] <- pmax(Clist[1,usevariab,drop=FALSE],tmp$coefficients + qnorm(1-alpha/4) * tmp$coefficientsvar,na.rm=TRUE)
            Clist[2,usevariab] <- pmin(Clist[2,usevariab,drop=FALSE],tmp$coefficients - qnorm(1-alpha/4) * tmp$coefficientsvar,na.rm=TRUE)
            for (kc in usevariab){
                Coeff[[kc]] <- c(Coeff[[kc]],tmp$coefficients[ which(usevariab==kc)])
                CoeffVar[[kc]] <- c(CoeffVar[[kc]],tmp$coefficientsvar[ which(usevariab==kc)])
            }
            if(length(notusevariab)>=1){
                Clist[1,notusevariab] <- pmax(Clist[1,notusevariab,drop=FALSE],0,na.rm=TRUE)
                Clist[2,notusevariab] <- pmin(Clist[2,notusevariab,drop=FALSE],0,na.rm=TRUE)
            }
            
            
            pvall <- c(pvall,pval)
        }        
    }
    colnames(Clist) <- colnames(X[,-1,drop=FALSE])
    if(is.null(colnames(Clist))) colnames(Clist) <- paste("Variable", 1:ncol(Clist))

    sig <- apply(sign(Clist[2:1,]),2,function(x) prod(x))
    sigo <- sign(Clist[1,])
    maximin <- sigo*apply(abs(Clist[2:1,]),2,min) * (sig>=0)
    
    retobj <- list(ConfInt=Clist,maximinCoefficients=maximin,alpha=alpha,colnames=colnames(Clist),factor=is.factor(Y),dimX=dim(X),Coeff=Coeff,CoeffVar=CoeffVar,modelReject=(max(sapply(Coeff,length))==0),usedvariables=usedvariables)
    class(retobj) <- "InvariantCausalPrediction"
    cat("\n\n")
    return(retobj)
}
