getdata<-function(status,tlabel){
  address = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
  ext = '.csv';
  filename = paste(address,'time_series_covid19_',status,'_global',ext,sep='')
  if (status=='confirmed'){
    results =read.csv(filename);
    
  }else if (status=='deaths'){
    results =read.csv(filename);
    
  }else if (status=='recovered'){
    results =read.csv(filename);
  }
  names(results)<-c('region','country','lat','long',tlabel)
  results$region <- paste('-',results$region)
  results$region[results$region=='- ']=''
  results$region <- paste(results$country,results$region)
  return(results)
}

helper<-function(Y,A,M){
  return(A%*%Y+M)
}

epidemic<-function(para,Npop,E0,I0,Q0,R0,D0,t,lambdaFun,kappaFun,fitflag=0){
  alpha = abs(para[1])
  beta = abs(para[2])
  gamma = abs(para[3])
  delta = abs(para[4])
  lambda0 = abs(para[5:7])
  kappa0 = abs(para[8:10])
  
  N = length(t);
  Y = matrix(0L, nrow = 7, ncol=N)
  Y[1,1] = Npop-Q0-E0-R0-D0-I0;
  Y[2,1] = E0;
  Y[3,1] = I0;
  Y[4,1] = Q0;
  Y[5,1] = R0;
  Y[6,1] = D0;
  if (abs(sum(Y[,1])-Npop)>0.5){
    print('Error: the sum must be equal to the total population')
    return()
  }
  dt = median(diff(t));
  lambda = lambdaFun(lambda0,t);
  kappa = kappaFun(kappa0,t);
  
  for (i in 1:(N-1)){
    
    A=matrix(0L, nrow = 7, ncol=7)
    A[1,1] = -alpha;
    A[2,2] = -gamma;
    A[3,2:3] = c(gamma,-delta);
    A[4,3:4] = c(delta,-kappa[i]-lambda[i]);
    A[5,4] = lambda[i];
    A[6,4] = kappa[i];
    A[7,1] = alpha;
    
    SI = Y[1,i]*Y[3,i];
    M = matrix(0L, nrow = 7, ncol=1);
    M[1:2,1] = c(-beta/Npop,beta/Npop)*SI;
    
    k_1 = helper(Y[,i],A,M);
    k_2 = helper(Y[,i]+0.5*dt*k_1,A,M);
    k_3 = helper(Y[,i]+0.5*dt*k_2,A,M);
    k_4 = helper(Y[,i]+k_3*dt,A,M);
    Y[,i+1] = Y[,i] + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    
  }
  if (fitflag==0){
   res=list(
    S = Y[1,1:N],
    E = Y[2,1:N],
    I = Y[3,1:N],
    Q = Y[4,1:N],
    R = Y[5,1:N],
    D = Y[6,1:N],
    P = Y[7,1:N],
    t=t)
   return(res)
  }else{
    t0=seq(1,length(t),by=10)
   res=Y[4:6,t0]
   res=matrix(res,ncol=1)
   return(res)
  }
}

getKappaFun<-function(tTarget,Q,D,guess){
  if (max(D)<10){
    kappaFun<-function(a,t){
      res=a[1]/(exp(a[2]*(t-a[3])) + exp(-a[2]*(t-a[3])));
      return(res)
    }
  } else{
    result<-tryCatch({
    helper1<-function(a,t){
      res=a[1]/(exp(a[2]*(t-a[3])) + exp(-a[2]*(t-a[3])));
      return(res)
    }
    kappaFun = helper1;
    helper2<-function(a,t){
      res=a[1]*exp(-a[2]*(t-a[3])^2);
      return(res)
    }
    helper3<-function(a,t){
      res=a[1] + exp(-a[2]*(t+a[3]));
      return(res)
    }
    rate = (diff(D)/median(diff(tTarget)))/Q[2:length(Q)];
    rate[abs(rate)>3]=NULL;
    dl=list(y=rate,t=tTarget[2:length(tTarget)])
    fit1<-nls(y~helper1(para,t),data=dl,start=list(para=guess[8:10]),lower=c(0,0,0),upper=c(1,1,100),algorithm = 'port',na.action = na.omit)
    coeff1<-coef(fit1)
    guess[8:10]<-coeff1;
    r1<-sigma(fit1)
    minR=r1;
    fit2<-nls(y~helper2(para,t),data=dl,start=list(para=guess[8:10]),lower=c(0,0,0),upper=c(1,1,100),algorithm = 'port',na.action = na.omit)
    coeff2<-coef(fit2)
    r2<-sigma(fit2)
    if (r2<minR){
      minR=r2;
      guess[8:10]<-coeff2;
      kappaFun = helper2;
      }
    fit3<-nls(y~helper3(para,t),data=dl,start=list(para=guess[8:10]),lower=c(0,0,0),upper=c(1,1,100),algorithm = 'port',na.action = na.omit)
    coeff3<-coef(fit3)
    r3<-sigma(fit3)
    if (r3<minR){
      guess[8:10]<-coeff3;
      kappaFun = helper3;
    }
    return(list(guess=guess,kappaFun=kappaFun))
    },error=function(err){
      kappaFun =  helper1;
      return(list(guess=guess,kappaFun=kappaFun))
    })
  }
  return(list(guess=guess,kappaFun=kappaFun))
}


getLambdaFun<-function(tTarget,Q,R,guess){
  if (max(R)<20){
    lambdaFun<-function(a,t){
      res=a[1]/(1+exp(-a[2]*(t-a[3])));
      return(res)
    }
    result=list(guess=guess,lambdaFun=lambdaFun)
  } else{
    result<-tryCatch({
      helper1<-function(a,t){
        res=a[1]/(1+exp(-a[2]*(t-a[3])));
        return(res)
      }
      lambdaFun = helper1;
      helper2<-function(a,t){
        res=a[1] + exp(-a[2]*(t+a[3]));
        return(res)
      }
      
      rate = (diff(R)/median(diff(tTarget)))/Q[2:length(Q)];
      rate[abs(rate)>1]=NULL;
      t=tTarget[2:length(tTarget)]
      dl=list(y=rate,t=t)
      fit1<-nls(y~helper1(para,t),data=dl,start=list(para=guess[5:7]),lower=c(0,0,0),upper=c(1,1,100),algorithm = 'port',na.action = na.omit)
      coeff1<-coef(fit1)
      guess[5:7] = coeff1;
      
      r1<-sigma(fit1)
      fit2<-nls(y~helper2(para,t),data=dl,start=list(para=guess[5:7]),lower=c(0,0,0),upper=c(1,1,100),algorithm = 'port',na.action = na.omit)
      coeff2<-coef(fit2)
      r2<-sigma(fit2)
      if ((r1<r2)|(coeff2[1]>0.99)|(coeff2[2]>4.9)){
        lambdaGuess = coeff1;
        lambdaFun = helper1;
        
      }else{
        lambdaGuess = coeff2;
        lambdaFun = helper2;
      }
      guess[5:7] = lambdaGuess;
      return(list(guess=guess,lambdaFun=lambdaFun))
      
    },error=function(err){
      #print(err)
      lambdaFun =  helper1;
      return(list(guess=guess,lambdaFun=lambdaFun))
    })
  }
  return(list(guess=guess,lambdaFun=lambdaFun))
}


fit_epidemic<-function(Q,R,D,Npop,E0,I0,tlabel,guess,maxiter){
  Q[Q<0]=0;
  R[R<0]=0;
  D[D<0]=0;
  input = rbind(Q,R,D);
  dt=0.1
  tol=1e-5
  fs = 1./dt;
  #print(tlabel)
  tTarget = 0:(length(tlabel)-1);
  t = seq(0, length(tlabel)-1, by = dt);
  result = getLambdaFun(tTarget,Q,R,guess);
  guess=result$guess;
  lambdaFun=result$lambdaFun; 
  result= getKappaFun(tTarget,Q,D,guess);
  guess=result$guess;
  kappaFun=result$kappaFun; 
  lambdaMax = c(1,5,100); 
  lambdaMin = c(0,0.000001,0);
  kappaMax = c(0.5,1,100);
  kappaMin = guess[8:10]/3;
  ub = c(1, 5, 1, 1, lambdaMax, kappaMax); 
  lb = c(0, 0, 0, 0, lambdaMin, kappaMin);
  print(guess)
  input=matrix(input,ncol=1)
  #names(guess)<-c('a0','a1','a2','a3','a4','a5','a6','a7','a8','a9')
  dl=list(y=input,t0=t,Npop=Npop,E0=E0,I0=I0,Q0=Q[1],R0=R[1],D0=D[1],lambdaFun=lambdaFun,kappaFun=kappaFun,fitflag=1)
  c=nls.control(maxiter = maxiter, tol = 1e-06, minFactor = 1/1024,printEval = FALSE, warnOnly = TRUE)
  #(alpha,beta,gamma,delta,lambda0,kappa0,Npop,E0,I0,Q0,R0,D0,t,lambdaFun,kappaFun,fitflag)
  fit<-nls(y~epidemic(para,Npop,E0,I0,Q0,R0,D0,t0,lambdaFun,kappaFun,fitflag),dl,start=list(para=guess),lower=lb,upper=ub,algorithm='port',trace = T,control=c)
  print(summary(fit))
  Coeff<-coef(fit)
  result=list(
    alpha1 = abs(Coeff[1]),
    beta1 = abs(Coeff[2]),
    gamma1 = abs(Coeff[3]),
    delta1 = abs(Coeff[4]),
    Lambda1 = abs(Coeff[5:7]),
    Kappa1 = abs(Coeff[8:10]),
    lambdaFun=lambdaFun,
    kappaFun=kappaFun
  )
}

checkRates<-function(time,Q,R,D,kappaFun,lambdaFun,Kappa1,Lambda1){
  t=0:(length(time)-1)
  rateD = (diff(D)/diff(t))/Q[2:length(Q)];
  if(isTRUE(max(abs(rateD)>3))){
  rateD[abs(rateD)>3] = NULL;}
  rateR = (diff(R)/diff(t))/Q[2:length(Q)];
  if(isTRUE(max(abs(rateR)>3))){
  rateR[abs(rateR)>3] = NULL;}
  x=1:(length(time)-1)
  xf=seq(1,length(time),0.1)
  yf=kappaFun(Kappa1,xf)
  d <- data.frame(x = x, y = rateD)
  df<-data.frame(xf=xf,yf=yf)
  colors <- c("Measured" = "blue", "Fitted" = "red")
  linetypes<-c('Fitted'=1)
  p1<-ggplot()+geom_point(data = d, mapping = aes(x = x, y = y,color='Measured')) + 
    geom_line(data = df,mapping = aes(x = xf, y = yf,linetype='Fitted'),color='red')+
    labs(x = "Time (days)",
         y = "Recovery rate (day^{-1})",
         color = "Legend",
         linetype = "Legend") +
    scale_color_manual('',values = colors)+
    scale_linetype_manual('',values = linetypes)
  
  yf=lambdaFun(Lambda1,xf)
  d <- data.frame(x = x, y = rateR)
  df<-data.frame(xf=xf,yf=yf)
  p2<-ggplot() + geom_point(data = d, mapping = aes(x = x, y = y,color='Measured'))+
    geom_line(data = df,mapping = aes(x = xf, y = yf,linetype='Fitted'),color='red')+
    labs(x = "Time (days)",
       y = "Recovery rate (day^{-1})",
       color = "Legend",
       linetype = "Legend") +
    scale_color_manual('',values = colors)+
    scale_linetype_manual('',values = linetypes)
  
  return(list(p1=p1,p2=p2))
  
}

fitregion<-function(region,tlabel,confirmed,deaths,recovered,maxiter=50,Npop=6e6){
  maxiter=as.numeric(maxiter)
  #print(region)
  indR=grep(region,recovered$region);
  indC=grep(region,confirmed$region);
  indD=grep(region,deaths$region);
  print(c(indR,indD,indC))
  N=length(recovered);
  R=as.matrix(recovered[indR,5:N]);
  C=as.matrix(confirmed[indC,5:N]);
  D=as.matrix(deaths[indD,5:N]);
  
  minNum= round(0.20*max(C));
  indRemoved = c(which(C<=minNum),which(is.na(C)))
  R<-R[-indRemoved];
  D<-D[-indRemoved];
  tlabel<-tlabel[-indRemoved];
  C<-C[-indRemoved];
  
  Npop= as.numeric(Npop)
  # Definition of the first estimates for the parameters
  alpha_guess = 0.2; # protection rate
  beta_guess = 1.0; # Infection rate
  LT_guess = 5; # latent time in days
  Q_guess = 0.1; # rate at which infectious people enter in quarantine
  lambda_guess = c(0.01,0.001,0.01); # recovery rate
  kappa_guess = c(0.001,0.001,10); # death rate
  
  guess = c(alpha_guess,beta_guess,1/LT_guess,Q_guess,lambda_guess,kappa_guess);
  
  # Initial conditions
  E0 = C[1]; # Initial number of exposed cases. Unknown but unlikely to be zero.
  I0 = C[1]; # Initial number of infectious cases. Unknown but unlikely to be zero.
  Q0 = C[1]-R[1]-D[1];
  R0 = R[1];
  D0 = D[1];
  A = C-R-D;
  A[A<0] = 0;
  #result=fit_epidemic(A,R,D,Npop,E0,I0,tlabel,guess,maxiter);
  result<-tryCatch({
    suppressWarnings(fit_epidemic(A,R,D,Npop,E0,I0,tlabel,guess,maxiter))
    },
    error=function(cond) {
      return(NULL)
    },
    warning=function(cond) {
    },
    finally={
    }
  )
  if(is.null(result)){
    return(NULL)
  }
  
  dt = 1/24;
  N = length(tlabel)+40;
  t = seq(0,N-1,by=dt);
  simulate = epidemic(c(result$alpha1,result$beta1,result$gamma1,result$delta1,result$Lambda1,result$Kappa1),Npop,E0,I0,Q0,R0,D0,t,result$lambdaFun,result$kappaFun);
  simulate=data.frame(simulate)
  realdata=data.frame(A=A,R=R,D=D,t=1:length(tlabel))
  
  d0=tlabel[1]
  
  pics=checkRates(tlabel,A,R,D,result$kappaFun,result$lambdaFun,result$Kappa1,result$Lambda1)
  p1<-pics$p1
  p2<-pics$p2
  #p1
  #p2
  #ggarrange(p1,p2,width=c(1.5,2))
  colors <- c("Active" = "red", "Recovered" = "green","Deceased"="blue")
  p3<-ggplot(data = realdata) +geom_point(aes(x = t, y = A,color='Active'),shape=1,size=2)+
    geom_point(aes(x = t, y = R,color='Recovered'),shape=1,size=2) +
    geom_point(aes(x = t, y = D,color='Deceased'),shape=1,size=2)+    
    geom_line(data = simulate,mapping = aes(x = t, y = Q),color='red')+
    geom_line(data = simulate,mapping = aes(x = t, y = R),color='green')+ 
    geom_line(data = simulate,mapping = aes(x = t, y = D),color='blue')+
    labs(x = paste("Time ( days since",d0,')'),y = "Number of cases",color = "Legend")+
    scale_color_manual(name='',values = colors)
  
  param=c(round(result$alpha1,2),round(result$beta1,2),round(result$gamma1,2),round(result$delta1,2),round(result$Kappa1,2),round(result$Lambda1,2))
  names(param)<-c('1','2','3','4','5','6','3','4','5','6')
  param=table(param)
  allpics=list(
    p1=p1,
    p2=p2,
    p3=p3,
    alpha=as.character(round(result$alpha1,2)),
    beta=as.character(round(result$beta1,2)),
    gammainv=as.character(round(result$gamma1,2)),
    delta=as.character(round(result$delta1,2)),
    kappa=as.character(round(result$Kappa1,2)),
    lambda=as.character(round(result$Lambda1,2)),
    param=param
  )

  return(allpics)
}
