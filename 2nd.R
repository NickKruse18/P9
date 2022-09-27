library("quantmod")

dBM = function(t,n){
  return(matrix(rnorm(t*n),t,n)/sqrt(t))
}


rBerVol = function(sigma,eta,H,dW){
  t = length(dW[,1])
  n = length(dW[1,])
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
  return(list(V,v))
}

BSV = function(S0,mu,dt,sigma,dZ){
  t = length(dZ[,1])
  n = length(dZ[1,])
  S = matrix(S0,t,n)
  for (i1 in 1:t){
    for (i2 in 1:n){
      S[i1,i2] = S[i1-(i1>1),i2] + S[i1-(i1>1),i2]*(mu*dt + sqrt(sigma[i1,i2])*dZ[i1,i2])
    }
  }
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

FitPowerLaw = function(Surface,start,end,tau){
  y = Surface[,end] - Surface[,start]
  plot(tau,y)
  s = sign(mean(exp(-c(1:length(y))/2)*y))
  Y = y*s
  X = cbind(1,log(tau[(Y>0)]))
  lines(tau[(Y>0)],s*exp(X%*%solve(t(X)%*%X)%*%t(X)%*%log(Y[(Y>0)])))
  b = solve(t(X)%*%X)%*%t(X)%*%log(Y[(Y>0)])
  return(c(exp(b[1])*s,b[2]))
}


ProcessOptions = function(D,date){
  n = 0
  for (i1 in 1:length(D)){
    n = n + length(D[[i1]]$calls[,1])
  }
  Data = matrix(0,n,3)
  i = 0
  for (i1 in 1:length(D)){
    for (i2 in 1:length(D[[i1]]$calls[,1])){
      i = i + 1
      Data[i,] = c(difftime(D[[i1]]$calls[i2,4],paste(date," 02:00:00"))[[1]],D[[i1]]$calls[i2,5],(D[[i1]]$calls[i2,9]+D[[i1]]$calls[i2,10])/2)
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
