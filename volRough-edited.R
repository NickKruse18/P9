
library("somebm")
library("rumidas")

data=rv5

m.proces=function(q,delta,data){
  m=0
  for (i in 1:(length(data)-delta)) {
    m=m+abs(log(data[[i]])-log(data[[i+delta]]))^q
  }
  return(m/(length(data)-delta))
}

m.delta=function(q,data,delta){
  m=numeric(length(delta))
  for (i in 1:length(delta)) {
    m[i]=m.proces(q,delta[i],data)
  }
  return(m)
}

#Fractional brownian motion
fBm = fbm(hurst=0.1,2000)
fm = m.delta(2,exp(fBm),1:20)

#volatility process
m=m.delta(2,1:1000,1:20)

#Heston model
HestonVol = function(V0,sigma,n){
  V = numeric(n); V[1] = V0
  for(i in 2:n){
    V[i] = V[i-1] + sigma*sqrt(V[i-1])*rnorm(1)/sqrt(n)
  }
  return(V)
}

hv = HestonVol(1,0.5,20000)
hm = m.delta(2,hv,1:20)

###Plots and fitting###
 
plot(log(1:20),log(fm))
plot(log(1:20),log(hm))
plot(log(1:20),log(m))


points(log(1:20),log(fm),col=4)
lines(log(1:20),log(1:20)*0.2513-1.6275)


lm(log(m)~log(1:20))
lm(log(hm)~log(1:20))
lm(log(fm)~log(1:20))

data[[1:365]]
#the scaling of the slope using Heston model is approximately 3*H



