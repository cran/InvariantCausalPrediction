summary.InvariantCausalPrediction <-
function(object,maxshow=50,...){
    cat(paste("\n Invariant Linear Causal", if(object$factor)" Classification " else " Regression ", "at level ", object$alpha,sep=""))
    usedvariables <- object$usedvariables
    if(length(usedvariables)< length(object$colnames)) cat(paste("\n\n ",length(usedvariables)," have been pre-selected for analysis -- all others \n    have been dropped by the pre-selection; \n use selection='all' in ICP-function to re-run with all variables if desired",sep=""))
    
    if( (p <- length(usedvariables))>maxshow){
        cat(paste("\n Number of variables (",p,") exceeds chosen limit (",maxshow,") \n ... only showing first ",maxshow," variables "  ,sep=""))
    }else{
        maxshow <- p
    }
    usedv <- object$usedvariables[1:maxshow]
    sig <- apply(sign(object$ConfInt[2:1,usedv]),2,function(x) prod(x))
    sigo <- sign(object$ConfInt[1,usedv])
    dis <- rbind(object$ConfInt[2:1,usedv], sigo*apply(abs(object$ConfInt[2:1,usedv]),2,min) * (sig>=0))
    colnames(dis) <- object$colnames[usedv]
    rownames(dis) <- c( "LOWER BOUND", "UPPER BOUND", " MAXIMIN EFFECT")
    cat("\n \n")
    if(object$modelReject){
        cat(paste("Model has been rejected due to presence of non-linearities or hidden variables \n\n",sep=""))
    }else{
        printCoefmat( t(dis), digits=3, signif.stars=FALSE,P.values=FALSE,has.Pvalue=FALSE)
    }
    cat("\n \n")

}
