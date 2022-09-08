
FBM = function(t,H){
  Sigma = matrix(0,nrow=length(t),ncol=length(t))
  a = t^(2*H)
  for (i1 in 1:length(t)){
    Sigma[i1,i1] = a[i1]
    if(i1 == length(t)){
      next
    }
    for (i2 in (i1+1):length(t)){
      Sigma[i1,i2] = 0.5*(a[i1]+a[i2]-a[abs(i1-i2)])
      Sigma[i2,i1] = Sigma[i1,i2]
    }
  }
  C = t(chol(Sigma))
  return(C)
}

FOU = function(t,x0,nu,alpha,m,W){
  X = c(x0,rep(0,length(t)))
  for(i in 1:length(t)){
    X[1+i] = X[i] + nu*(W[1+i]-W[i])-alpha*(X[i]-m)*(t[2]-t[1])
  }
  return(X)
}

RSFV = function(t,x0,mu,V,dZ){
  X = c(x0,rep(0,length(t)))
  for(i in 1:length(t)){
    X[1+i] = X[i] + X[i]*(mu*(t[2]-t[1])+V[i+1]*dZ[i])
  }
  return(X)
}

rVol = function(n,t,x0,y0,mu,nu,alpha,m,H){
  s = Sys.time()
  Sigma = FBM(t,H)
  N = rep(0,n)
  for (i in 1:n){
    W = c(0,Sigma%*%rnorm(length(t)))
    Y = FOU(t,y0,nu,alpha,m,W)
    X = RSFV(t,x0,mu,exp(Y),rnorm(length(t))*sqrt(diff(c(0,t))))
    N[i] = X[length(X)]
  }
  print(Sys.time()-s)
  return(N)
}

set.seed(2)
hist(rVol(1000,1:1000/1000,1,0,0,0.25,0,0,0.1))

Price = function(n,S0,K,sigma0,mu,tau,nu,dt,dZ,dW){
  price = rep(0,n)
  for (i1 in 1:n){
    S = S0
    sigma = sigma0
    for (i2 in (floor(length(dZ)/n)-floor(tau/dt)+1):floor(length(dZ)/n)){
      sigma = sigma + nu*dW[i1,i2]
      S = S + S*(mu*dt + exp(sigma)*dZ[i1,i2])
    }
    price[i1] = max(0,S-K)
  }
  return(price)
}

mean(Price(1000,1,1,0,0,1,0.01,0.01,matrix(rnorm(100000),1000,100)/sqrt(100),t(FBM(1:100/100,0.1)%*%matrix(rnorm(100000),100,1000))))

(pnorm(1/2)-pnorm(-1/2))

#S exp(mu tau) N(d1) - K N(d2)
#d1 = (log(S / K) + (mu + 1/2 sigma^2) tau) / sigma / sqrt(tau)
#d2 = (log(S / K) + (mu - 1/2 sigma^2) tau) / sigma / sqrt(tau)

D = getOptionChain("AAPL",Exp=NULL,src="yahoo")


n = 0
for (i1 in 1:length(D)){
  n = n + length(D[[i1]]$calls[,1])
}

#time K/S price/S
Data = matrix(0,n,3)
i = 0

for (i1 in 1:length(D)){
  S = 155.4
  for (i2 in 1:length(D[[i1]]$calls[,1])){
    i = i + 1
    Data[i,] = c(difftime(D[[i1]]$calls[i2,4],"2022-09-07 02:00:00")[[1]],D[[i1]]$calls[i2,5]/S,(D[[i1]]$calls[i2,9]+D[[i1]]$calls[i2,10])/2/S)
  }
}


Price = function(n,S0,K,P,sigma,mu,tau,dt,dZ){
  price = rep(0,n)
  for (i1 in 1:n){
    S = S0
    for (i2 in (floor(length(dZ)/n)-floor(tau/dt)+1):floor(length(dZ)/n)){
      S = S + S*(mu*dt + sigma*dZ[i1,i2])
    }
    price[i1] = max(0,S-K)
  }
  return((mean(price)-P)^2)
}
deltaZ = matrix(rnorm(1306000),1000,1306)/sqrt(20)
Data = cbind(Data,0)
for (i in 1:length(Data[,1])){
  print(i)
  Data[i,4] = optimize(Price,c(0,1), tol=0.01, n = 1000,S0 = 155.4,K = Data[i,2]*155.4,P = Data[i,3]*155.4,mu = 0,tau = Data[i,1], dt = 0.5, dZ = deltaZ)$minimum
}


Price(n = 1000,S0 = 155.4,K = Data[1,2]*155.4,P = Data[1,3]*155.4,sigma = -2,mu = 0,tau = Data[1,1], dt = 0.1, dZ = deltaZ)


plot(Data[Data[,1]==135,2],Data[Data[,1]==135,4])



