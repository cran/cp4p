
############################################################
# plot to verify the assumptions of FDR control procedures #
############################################################

calibration.plot=function(p, pi0.method="pounds",  nbins = 20, pz=0.05){
  #
  c=max(p);
  Fr=ecdf(1-p);
  minabs=max(0, min(1-p)-0.05);
  abs=seq(minabs,1,by=0.001);
  fc=Fr(abs);
  AUC=NULL;
  AUC2=NULL;
  plot(abs,fc,ty="l",xlab="1-p.value",ylab="Cumulative Distribution Function of 1-p.value",lwd=2,ylim=c(0,1.05));
  if (is.numeric(pi0.method)==TRUE){
    if (pi0.method<=1 && pi0.method>=0){
      pi0=pi0.method;
      title(main=paste("Calibration Plot - pi0 =",pi0.method))

      #straight line y=pi0*x
      dr=as.numeric(pi0)*(abs-1+c)/c;
      dr[dr<0]=0;
      dif=(abs-1+c)/c-fc;
      dif[dif<0]=0;
      dif[dif==0]=0.0001;
      f=fc[dif>0]-dr[dif>0];
      f[f<0]=0;
      f=f+dr[dif>0];
      #Gray area
      id=1;
      while (f[id]==0){id=id+1;}
      x.bons=c(abs[dif>0],abs[dif>0][id-1]);
      y.bons=c(f,f[id-1]);
      polygon(x = x.bons, y = y.bons, col = "gray93",border="white");
    
      dif=fc-dr;
      lfc=length(fc);
      i=lfc;
      if (sign(dif[lfc])>0){
        while (sign(dif[i-1])==sign(dif[i])){i=i-1;}
      }
      #Green area
      x.bons=c(abs[i:lfc],1,abs[i],abs[i]);
      y.bons=c(fc[i:lfc],pi0,pi0*(abs[i]-1+c)/c,fc[i]);
      polygon(x = x.bons, y = y.bons, col = "darkseagreen1",border="darkseagreen1");  
      #Compute DA protein concentration
      if (i<lfc){auc1=MESS::auc(abs[i:lfc],dif[i:lfc]);}else{auc1=0;}
      if (c*(1-pi0)>2*auc1){
        AUC=(c*(1-pi0)/2-auc1)/(c*(1-pi0)/2);
      }else{AUC=0;}
      text(minabs,0.9,paste("DA protein concentration: ",floor(AUC*1000)/10,"%"),col="green4",pos=4);
      #Compute divergence to uniformity assumption
      dif2=dif;
      dif2[dif2<0]=0;
      fc2=dr+dif2;
      AUC2=MESS::auc(abs[1:(i-1)],dif2[1:(i-1)]);
      text(minabs,0.85,paste("Uniformity underestimation: ",floor(AUC2*100000)/1000),col="red",pos=4);
      #Red area
      x.bons=c(abs[id:(i-1)],abs[id]);
      y.bons=c(fc2[id:(i-1)],fc2[id]);
      polygon(x = x.bons, y = y.bons, col = "red",border="white");
    
      lines(abs,dr,lwd=2,col="blue",lty=1);
      lines(abs,fc,lwd=2);
      text(minabs,0.95,paste("Non-DA protein proportion: ",floor(pi0*1000)/10, "%"),col="blue",pos=4);
    }else{warning("\n Error in input pi0.method:\n Please write a numeric value between 0 and 1 or the name of an estimation method among st.spline, st.boot, jiang, histo, langaas, pounds, abh, slim or ALL.\n");}
  }
  if (is.numeric(pi0.method)==FALSE){
    if (pi0.method=="ALL"||pi0.method=="st.spline"||pi0.method=="st.boot"||pi0.method=="jiang"||pi0.method=="histo"||pi0.method=="langaas"||pi0.method=="pounds"||pi0.method=="abh"||pi0.method=="slim"){
      r=cp4p::estim.pi0(p=p, pi0.method = pi0.method,  nbins = nbins, pz=pz);
      pi0=r$pi0;
      if (pi0.method=="ALL"){
        for (i in 1:length(pi0)){
          if (i==1){dro=as.numeric(pi0[i])*(abs-1+c)/c;dro[dro<0]=0;lines(abs,dro,ty="l",col="chartreuse1", lty=i,lwd=2);name="st.spline";}
          if (i==2){dro=as.numeric(pi0[i])*(abs-1+c)/c;dro[dro<0]=0;lines(abs,dro,ty="l",col="chartreuse3", lty=i,lwd=2);name=c(name,"st.boot");}
          if (i==3){dro=as.numeric(pi0[i])*(abs-1+c)/c;dro[dro<0]=0;lines(abs,dro,ty="l",col="darkgreen", lty=i,lwd=2);name=c(name,"jiang");}
          if (i==4){dro=as.numeric(pi0[i])*(abs-1+c)/c;dro[dro<0]=0;lines(abs,dro,ty="l",col="darkorange", lty=i,lwd=2);name=c(name,"histo");}
          if (i==5){dro=as.numeric(pi0[i])*(abs-1+c)/c;dro[dro<0]=0;lines(abs,dro,ty="l",col="darkorange3", lty=i,lwd=2);name=c(name,"langaas");}
          if (i==6){dro=as.numeric(pi0[i])*(abs-1+c)/c;dro[dro<0]=0;lines(abs,dro,ty="l",col="red", lty=i,lwd=2);name=c(name,"pounds");}
          if (i==7){dro=as.numeric(pi0[i])*(abs-1+c)/c;dro[dro<0]=0;lines(abs,dro,ty="l",col="deepskyblue", lty=i,lwd=2);name=c(name,"abh");}
          if (i==8){dro=as.numeric(pi0[i])*(abs-1+c)/c;dro[dro<0]=0;lines(abs,dro,ty="l",col="blue", lty=i,lwd=2);name=c(name,"slim");}
        }
        title(main="Calibration Plot - All methods");
        legend("topleft",name,col=c("chartreuse1","chartreuse3","darkgreen","darkorange","darkorange3","red","deepskyblue","blue"),lty=1:8,lwd=rep(2,8),title="pi0.method");
      }else{
        title(main=paste("Calibration Plot -",pi0.method,"method"));
        
        #straight line y=pi0*x
        dr=as.numeric(pi0)*(abs-1+c)/c;
        dr[dr<0]=0;
        dif=(abs-1+c)/c-fc;
        dif[dif<0]=0;
        dif[dif==0]=0.0001;
        f=fc[dif>0]-dr[dif>0];
        f[f<0]=0;
        f=f+dr[dif>0];
        #Gray area
        id=1;
        while (f[id]==0){id=id+1;}
        x.bons=c(abs[dif>0],abs[dif>0][id-1]);
        y.bons=c(f,f[id-1]);
        polygon(x = x.bons, y = y.bons, col = "gray93",border="white");
        
        dif=fc-dr;
        lfc=length(fc);
        i=lfc;
        if (sign(dif[lfc])>0){
          while (sign(dif[i-1])==sign(dif[i])){i=i-1;}
        }
        #Green area
        x.bons=c(abs[i:lfc],1,abs[i],abs[i]);
        y.bons=c(fc[i:lfc],pi0,pi0*(abs[i]-1+c)/c,fc[i]);
        polygon(x = x.bons, y = y.bons, col = "darkseagreen1",border="darkseagreen1");  
        #Compute DA protein concentration
        if (i<lfc){auc1=MESS::auc(abs[i:lfc],dif[i:lfc]);}else{auc1=0;}
        if (c*(1-pi0)>2*auc1){
          AUC=(c*(1-pi0)/2-auc1)/(c*(1-pi0)/2);
        }else{AUC=0;}
        text(minabs,0.9,paste("DA protein concentration: ",floor(AUC*1000)/10,"%"),col="green4",pos=4);
        #Compute divergence to uniformity assumption
        dif2=dif;
        dif2[dif2<0]=0;
        fc2=dr+dif2;
        AUC2=MESS::auc(abs[1:(i-1)],dif2[1:(i-1)]);
        text(minabs,0.85,paste("Uniformity underestimation: ",floor(AUC2*100000)/1000),col="red",pos=4);
        #Red area
        x.bons=c(abs[id:(i-1)],abs[id]);
        y.bons=c(fc2[id:(i-1)],fc2[id]);
        polygon(x = x.bons, y = y.bons, col = "red",border="white");
        
        lines(abs,dr,lwd=2,col="blue",lty=1);
        lines(abs,fc,lwd=2);
        text(minabs,0.95,paste("Non-DA protein proportion: ",floor(pi0*1000)/10, "%"),col="blue",pos=4);
      }
    }else{warning("\n Error in input pi0.method:\n Please write a numeric value between 0 and 1 or the name of an estimation method among st.spline, st.boot, jiang, histo, langaas, pounds, abh, slim or ALL.\n");}
  }
  return(list(pi0=pi0,h1.concentration=AUC,unif.under=AUC2*100));
}
