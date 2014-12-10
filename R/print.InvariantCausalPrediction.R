print.InvariantCausalPrediction <-
function(x,...){
    ns <- sum(x$maximinCoefficients!=0)
    wh <- which( x$maximinCoefficients!=0)
    cat(paste("\n Invariant Linear Causal", if(x$factor)" Classification " else " Regression ", "at level ", x$alpha,sep=""))
    if(x$modelReject){
        cat(paste("\n Model has been rejected, for example due to presence of non-linearities or hidden variables \n\n",sep=""))
    }else{
        if(ns>0){
            cat(paste("\n ",if(ns>1) "Variables: " else "Variable ",paste(x$colnames[wh[1:min(10,length(wh))]],collapse=", ")  ,if(length(wh)>10) paste("... (and ", length(wh)-10, "more) ",sep="") else ""," ", if(ns>1) "show" else "shows"," a significant causal effect",sep=""))
            sig <- apply(sign(x$ConfInt[2:1,]),2,function(x) prod(x))
            sigo <- sign(x$ConfInt[1,])
            dis <- rbind(x$ConfInt[2:1,], sigo*apply(abs(x$ConfInt[2:1,]),2,min) * (sig>=0))
            colnames(dis) <- x$colnames
            rownames(dis) <- c( "LOWER BOUND", "UPPER BOUND", " MAXIMIN EFFECT")
            dis <- dis[ ,wh,drop=FALSE]
            cat("\n \n ")
            printCoefmat( t(dis), digits=5, signif.stars=FALSE,P.values=FALSE,has.Pvalue=FALSE)

        }else{
            cat(paste("\n \n No variables show a significant causal effect"))
        }
    }
    cat("\n\n")
}
