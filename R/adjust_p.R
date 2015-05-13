############################
# compute ajusted p-values #
############################

adjust.p=function(p, pi0.method=1, alpha=0.05, nbins = 20, pz=0.05){
  #library(multtest)
  if ((pi0.method=="bky")==0){
    if (is.numeric(pi0.method)==TRUE){pi0=pi0.method;
                                      if (pi0.method==1){cat("Procedure of Benjamini-Hochberg is used. pi0 is fixed to 1.");}
    }
    if (is.numeric(pi0.method)==FALSE){
      r=estim.pi0(p, pi0.method = pi0.method,  nbins = nbins, pz=pz);
      pi0=r$pi0
      cat(pi0.method,"method is used. pi0 is estimated to",as.numeric(pi0),".");
    }
    qa=mt.rawp2adjp(p, proc = "BH");
    adjp=data.frame(qa$adjp[order(qa$index),1],qa$adjp[order(qa$index),2]*as.numeric(pi0));
    colnames(adjp)=c("rawp","adjusted.p");
    return(list(adjp=adjp,pi0=pi0));
  }
  else{
    qa=mt.rawp2adjp(p, proc = "TSBH", alpha=alpha);
    adjp=data.frame(qa$adjp[order(qa$index),1:2]);
    pi0=as.numeric(qa$h0.TSBH)/length(p);
    cat("Procedure of BKY with a FDR of",alpha,"is used. pi0 is estimated to",pi0,".");
    colnames(adjp)=c("rawp","adjusted.p");
    return(list(adjp=adjp,pi0=pi0));
  }
}
