library("quantmod")
library("stats")
library("hypergeo")
library("microbenchmark")


dBM = function(t,n){
  return(matrix(rnorm(t*n),t,n)/sqrt(t))
}

HesVol = function(sigma,eta,dW){
  t = length(dW[,1])
  n = length(dW[1,])
  V = matrix(sigma,t,n)
  for(i1 in 1:t){
    for(i2 in 1:n){
      V[i1,i2] = V[i1,i2] + eta*sqrt(V[i1,i2])*dW[i1,i2]
    }
  }
  return(V)
}

rBerVol = function(sigma,eta,H,t,n){
  dW = dBM(t,n)
  V = matrix(sigma,t,n)
  v = matrix(sigma,t,n)
  for(i1 in 1:t){
    for(i2 in 1:i1){
      gamma = eta*sqrt(2*H)*((i1-i2+1)/t)^(H-1/2)
      for(i3 in 1:n){
        V[i1,i3] = V[i1,i3] + V[i1,i3]*gamma*dW[i2,i3]
      }
    }
  }
  for(i2 in 1:t){
    gamma = eta*sqrt(2*H)*((t-i2+1)/t)^(H-1/2)
    for(i3 in 1:n){
      v[i2,i3] = v[i2-(i2>1),i3] + v[i2-(i2>1),i3]*gamma*dW[i2,i3]
    }
  }
  return(list(V,v,dW))
}

BSV = function(S0,mu,dt,sigma,dZ){
  t = length(dZ[,1])
  n = length(dZ[1,])
  S = matrix(0,t,n)
  for (i1 in 1:t){
    for (i2 in 1:n){
      S[i1,i2] = S[i1-(i1>1),i2] + (mu-sigma[i1,i2]/2)*dt + sqrt(sigma[i1,i2])*dZ[i1,i2]
    }
  }
  S = S0*exp(S)
  return(S)
}

CoVT = function(X,Y){
  cov = rep(0,length(X[,1]))
  for(i1 in 1:length(X[,1])){
    cov[i1] = t(X[i1,])%*%Y[i1,]/length(X[1,])
  }
  return(cov)
}

PriceSurface = function(S,t,K){
  Surface = matrix(0,length(t),length(K))
  for(i1 in 1:length(t)){
    for(i2 in 1:length(K)){
      Surface[i1,i2] = mean((abs(S[t[i1],]-K[i2])+S[t[i1],]-K[i2])/2)
    }
  }
  return(Surface)
}

FitPrice = function(S0,K,mu,P,tau,sigma0){
  if (sigma0 == 0){
    d1 = Inf
    d2 = Inf
  }
  else{
    d1 = (log(S0/K)+tau*(mu+sigma0^2/2))/(sigma0*sqrt(tau))
    d2 = (log(S0/K)+tau*(mu-sigma0^2/2))/(sigma0*sqrt(tau))
  }
  price = pnorm(d1)*S0-pnorm(d2)*K*exp(-mu*tau)
  return((price-P)^2)
}

ImpliedVol = function(S,tau,K,Surface,mu){
  VolSurface = matrix(0,length(tau),length(K))
  for (i1 in 1:length(tau)){
    for (i2 in 1:length(K)){
      VolSurface[i1,i2] = optimize(FitPrice,c(0,10), S0 = S, K = K[i2], P = Surface[i1,i2], mu = mu, tau = tau[i1])$minimum
    }
  }
  return(VolSurface)
}

FitPowerLaw = function(Surface,start,end,tau,titel = ""){
  if(end == 0){ y = Surface }
  else{ y = Surface[,end] - Surface[,start] }
  
  plot(tau,y,xlab = expression(tau),ylab=expression(psi(tau)),main = titel)
  s = sign(mean(exp(-c(1:length(y))/2)*y))
  Y = y*s
  X = cbind(1,log(tau[(Y>0)]))
  lines(tau[(Y>0)],s*exp(X%*%solve(t(X)%*%X)%*%t(X)%*%log(Y[(Y>0)])))
  b = solve(t(X)%*%X)%*%t(X)%*%log(Y[(Y>0)])
  return(c(exp(b[1])*s,b[2]))
}


ProcessOptions = function(D,date){
  n = 0
  for (i in 1:length(D)){
    if (length(D[[i]]$calls[,1])==0){next}
    if (difftime(D[[i]]$calls[1,4],paste(date," 02:00:00"))[[1]]>1000){next}
    if (difftime(D[[i]]$calls[1,4],paste(date," 02:00:00"))[[1]] == 0){next}
    n = n + length(D[[i]]$calls[,1])
  }
  Data = matrix(0,n,5)
  i = 0
  for (i1 in 1:length(D)){
    if (length(D[[i1]]$calls[,1])==0){next}
    if (difftime(D[[i1]]$calls[1,4],paste(date," 02:00:00"))[[1]]>1000){next}
    if (difftime(D[[i1]]$calls[1,4],paste(date," 02:00:00"))[[1]] == 0){next}
    for (i2 in 1:length(D[[i1]]$calls[,1])){
      i = i + 1
      Data[i,] = c(difftime(D[[i1]]$calls[i2,4],paste(date," 02:00:00"))[[1]],D[[i1]]$calls[i2,5],D[[i1]]$calls[i2,9],D[[i1]]$calls[i2,14],D[[i1]]$calls[i2,12])
      Data[is.na(Data)] = 0
    }
  }
  return(Data)
}

AddImpliedVol = function(Data,X,r){
  Data = cbind(Data,0) 
  for (i in 1:length(Data[,1])){
    Data[i,length(Data[1,])] = optimize(FitPrice,c(0,1), S0 = X, K = Data[i,2],mu = r, P = Data[i,3],tau = Data[i,1])$minimum
  }
  return(Data)
}

Uniques = function(X){
  Xs = c()
  for(i in X){
    if (!(i %in% Xs)){ Xs = c(Xs,i)}
  }
  return(sort(Xs))
}

InterpolateSurface = function(X,Y,Z){
  Xs = Uniques(X)
  Ys = Uniques(Y)
  Zs = matrix(0,length(Xs),length(Ys))
  Zw = matrix(0,length(Xs),length(Ys))
  W = matrix(0,2*length(Xs)-1,2*length(Ys)-1)
  for (i1 in 1:(2*length(Xs)-1)){ for (i2 in 1:(2*length(Ys)-1)){ W[i1,i2] = exp(-(i1-length(Xs))^2/2-(i2-length(Ys))^2/2) } }
  W[length(Xs),length(Ys)] = 2
  for(i in 1:length(X)){
    Xind = match(X[i],Xs)
    Yind = match(Y[i],Ys)
    Zw = Zw + W[(length(Xs)-Xind+1):(2*length(Xs)-Xind),(length(Ys)-Yind+1):(2*length(Ys)-Yind)]
    Zs = Zs + Z[i]*W[(length(Xs)-Xind+1):(2*length(Xs)-Xind),(length(Ys)-Yind+1):(2*length(Ys)-Yind)]
  }
  persp(Xs,Ys,Zs/Zw,theta=120,phi=30)
  return(list(Xs,Ys,Zs/Zw))
}

AssetVolatility = function(Asset,Source,X,date){
  D = getOptionChain(Asset,Exp=NULL,src=Source)
  
  Data = ProcessOptions(D,date)
  #Data = AddImpliedVol(Data,X,r/365)
  Surface = InterpolateSurface(Data[,1],Data[,2],Data[,3])
  Surface = InterpolateSurface(Data[,1],Data[,2],Data[,4])
  
  s = 1
  for(i in (1:length(Surface[[2]]))){
    if(Surface[[2]][i] > X){
      s = i
      break
    }
  }
  print(FitPowerLaw(Surface[[3]],s-1+12,s+12,Surface[[1]]))
  return(list(Data,cbind(Surface[[1]],Surface[[3]][,s+12]-Surface[[3]][,s-1+12]),X))
}

rBerVol2 = function(sigma,eta,H,t,n,k){
  Y = matrix(sigma,t,n)
  Sigma = matrix(1/t,k+1,k+1)
  for(i in 2:(k+1)){
    Sigma[i,i] = ((i-1)^(2*H)-(i-2)^(2*H))/((2*H)*t^(2*H))
    Sigma[1,i] = ((i-1)^(H+1/2)-(i-2)^(H+1/2))/((H+1/2)*t^(H+1/2))
    Sigma[i,1] = Sigma[1,i]
    if(i == 2){next}
    for(j in 2:(i-1)){
      Sigma[i,j] = ((j-1)^(H+1/2)*(i-1)^(H-1/2)*Re(hypergeo(1/2-H,1,H+3/2,(j-1)/(i-1))))/((H+1/2)*t^(2*H))
      Sigma[i,j] = Sigma[i,j] - ((j-2)^(H+1/2)*(i-2)^(H-1/2)*Re(hypergeo(1/2-H,1,H+3/2,(j-2)/(i-2))))/((H+1/2)*t^(2*H))
      Sigma[j,i] = Sigma[i,j]
    }
  }
  Sigma = t(chol(Sigma))
  set.seed(2)
  dW = Sigma%*%matrix(rnorm((k+1)*t*n),k+1,t*n)
  g = c(rep(0,k), (((k+1):t)/t)^(H-1/2),rep(0,t))
  G = fft(g)
  dWW = dW[2:(k+1),]
  dW = dW[1,]
  for(i in 2:k){
    dWW[i,] = c(rep(0,i-1),dWW[i,1:(t*n-i+1)])
  }
  for(i1 in 1:(n-1)){
    for(i2 in 2:k){
      for(i3 in 1:(i2-1)){
        dWW[i2,i1*t+i3] = 0
      }
    }
  }
  V = matrix(rep(1,k)%*%dWW,t,n)
  for(i in 1:n){
    DW = fft(c(rep(0,t),dW[t*i-t+(1:t)]))
    V[,i] = V[,i] + round(Re(fft(DW*G,inverse=TRUE)[(t+1):(2*t)])/t/2,9)
  }
  
  
  V = sigma*exp(eta*sqrt(2*H)*V-eta^2/2*(((1:t)/t)^(2*H))%*%t(rep(1,n)))
  return(list(V,matrix(dW,t,n)))
}

ModelDelta = function(S0,K,mS,tau){
  n = length(mS[1,])
  t = length(mS[,1])
  delta = mean(0.0002*mS[t,]/mS[tau,]+abs((S0+0.0001)*mS[t,]/mS[tau,]-K)-abs((S0-0.0001)*mS[t,]/mS[tau,]-K))/0.0004
  return(delta)
}

CallOpt = function(S0,K,mS,tau){
  n = length(mS[1,])
  t = length(mS[,1])
  return(mean(S0*S[t,]/S[tau,]-K+abs(S0*S[t,]/S[tau,]-K))/2)
}

Hedging = function(S,K,mS,tau){
  n = length(mS[1,])
  t = length(mS[,1])
  Positions = rep(0,n)
  Cost = rep(-CallOpt(mean(S[1,]),K,mS,1),n)
  for(i1 in tau){
    M = (0:99)*(max(S[i1,])-min(S[i1,]))/99+min(S[i1,])
    Delta = rep(0,100)
    for(i2 in 1:100){
      Delta[i2] = ModelDelta(M[i2],K,mS,i1)
    }
    for(i2 in 1:n){
      for(i3 in 1:100){ if(M[i3] > S[i1,i2]){ break } }
      delta = Delta[i3-1] + (M[i3] - S[i1,i2])/(M[i3] - M[i3-1]) * (Delta[i3] - Delta[i3-1])
      Cost[i2] = Cost[i2] - S[i1,i2]*(delta + Positions[i2])
      Positions[i2] = -delta
    }
  }
  for(i2 in 1:n){
    Cost[i2] = Cost[i2] - S[t,i2]*Positions[i2] + max(0,S[t,i2]-K)
  }
  return(Cost)
}

rBerSurface = function(H,t,n,Tau,K){
  sigma0 = H[1]
  eta = H[2]
  rho = H[4]
  r = H[5]
  S0 = H[6]
  H = H[3]
  B = rBerVol2(sigma0,eta,H,t,n,3)
  dW = B[[2]]
  B = B[[1]]
  dZ = rho*dW+sqrt(1-rho^2)*dBM(t,n)
  S = BSV(S0,r,1/t,B,dZ)
  PSurface = PriceSurface(S,Tau,K)
  return(PSurface)
}

BSSurface = function(H,t,n,Tau,K){
  sigma = H[1]
  r = H[2]
  S0 = H[3]
  dZ = dBM(t,n)
  S = BSV(S0,r,1/t,matrix(sigma,t,n),dZ)
  PSurface = PriceSurface(S,Tau,K)
  return(PSurface)
}

PriceVal = function(S,t,K){
  Surface = rep(0,length(t))
  for(i in 1:length(t)){
    Surface[i] = mean((abs(S[t[i],]-K[i])+S[t[i],]-K[i])/2)
  }
  return(Surface)
}

CleanData = function(Fit,SData){
  Residuals = Fit-SData[,3]
  lims = boxplot(Residuals)$stats[c(1,5)]
  lims = c(min(lims[1],sort(Residuals)[floor(length(Residuals)*0.025)]),max(lims[2],sort(Residuals)[ceiling(length(Residuals)*0.975)]))
  exclude = c()
  for(i in 1:length(Residuals)){
    if(Residuals[i] < lims[1]){
      exclude = c(exclude,i)
    }
    if(Residuals[i] > lims[2]){
      exclude = c(exclude,i)
    }
  }
  return(SData[-exclude,])
}

OptBS = function(H,S0,t,n,SData,Est=FALSE){
  print(H)
  set.seed(2)
  dZ = dBM(t,n)
  S = BSV(S0,H[2],1/t,matrix(exp(H[1]),t,n),dZ)
  PEst = PriceVal(S,SData[,1],SData[,2])
  if(Est){return(PEst)}
  else{return(sum((PEst-SData[,3])^2))}
}

OptrBer = function(H,S0,t,n,SData,Est=FALSE){
  print(H)
  if(abs(H[3])>6){ return(Inf) }
  B = rBerVol2(exp(H[1]),H[2],0.5/(1+exp(H[3])),t,n,3)
  dW = B[[2]]
  B = B[[1]]
  G = 2*1/(1+exp(H[4]))-1
  dZ = G*dW+sqrt(1-G^2)*dBM(t,n)
  S = BSV(S0,H[5],1/t,B,dZ)
  
  PEst = PriceVal(S,SData[,1],SData[,2])
  if(Est){return(PEst)}
  else{return(sum((PEst-SData[,3])^2))}
}

FitModel = function(H,S0,t,n,SData,model){
  if (model == "BS"){
    Model = optim(H,OptBS,S0=S0,t=t,n=n,SData=SData,control = list(maxit=100))
    Fit = OptBS(Model$par,S0,t,n,SData,TRUE)
    print(1)
    SData.clean = CleanData(Fit,SData)
    Model = optim(Model$par,OptBS,S0=S0,t=t,n=n,SData=SData.clean)
    Fit = OptBS(Model$par,S0,t,n,SData.clean,TRUE)
  }
  else if (model == "rBer"){
    Model = optim(H,OptrBer,S0=S0,t=t,n=n,SData=SData,control = list(maxit=100))
    Fit = OptrBer(Model$par,S0,t,n,SData,TRUE)
    print(1)
    SData.clean = CleanData(Fit,SData)
    Model = optim(Model$par,OptrBer,S0=S0,t=t,n=n,SData=SData.clean)
    Fit = OptBS(Model$par,S0,t,n,SData.clean,TRUE)
  }
  else if (model == "Hes"){
    Model = optim(H,OptHes,S0=S0,t=t,n=n,SData=SData,control = list(maxit=100))
    Fit = OptHes(Model$par,S0,t,n,SData,TRUE)
    print(1)
    SData.clean = CleanData(Fit,SData)
    Model = optim(Model$par,OptHes,S0=S0,t=t,n=n,SData=SData.clean)
    Fit = OptBS(Model$par,S0,t,n,SData.clean,TRUE)
  }
  print(sqrt(Model$value/length(SData[,3]))/S0)
  return(list(Model$par,Fit,SData.clean))
}
