
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

fBm = fbm(0.15,20000)
warnings()
fm = m.delta(3,exp(fBm/2),1:20)
#-0.5761       0.4594 
plot(log(1:20),log(hm))
m=m.delta(3,data,1:20)
points(log(1:20),log(m),col=3)
lines(log(1:20),log(1:20)*0.3761-0.6955)
lm(log(hm)~log(1:20))
plot(log(1:20),log(m),col="blue")
data[[1:365]]
3*H

hv = HestonVol(1,0.5,20000)
hm = m.delta(3,hv,1:20)

HestonVol = function(V0,sigma,n){
  V = numeric(n); V[1] = V0
  for(i in 2:n){
    V[i] = V[i-1] + sigma*sqrt(V[i-1])*rnorm(1)/sqrt(n)
  }
  return(V)
}


