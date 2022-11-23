
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
fBm = fbm(hurst=0.15,20000)*3
fm.1 = m.delta(1,exp(fBm),1:20)
fm.2 = m.delta(1.5,exp(fBm),1:20)
fm.3 = m.delta(2,exp(fBm),1:20)
fm.4 = m.delta(2.5,exp(fBm),1:20)
fm.5 = m.delta(3,exp(fBm),1:20)

#volatility process
mvol.1=m.delta(1,data,1:20)
mvol.2=m.delta(1.5,data,1:20)
mvol.3=m.delta(2,data,1:20)
mvol.4=m.delta(2.5,data,1:20)
mvol.5=m.delta(3,data,1:20)

#Heston model
HestonVol = function(V0,sigma,n){
  V = numeric(n); V[1] = V0
  for(i in 2:n){
    V[i] = V[i-1] + sigma*sqrt(V[i-1])*rnorm(1)/sqrt(n)
  }
  return(V)
}

hv = HestonVol(1,0.5,20000)
hm.1 = m.delta(1,hv,1:20)
hm.2 = m.delta(1.5,hv,1:20)
hm.3 = m.delta(2,hv,1:20)
hm.4 = m.delta(2.5,hv,1:20)
hm.5 = m.delta(3,hv,1:20)


########################################### Plots and fitting ############################################
plot(log(1:20),log(mvol.5))

points(log(1:20),log(mvol.1),col=5)
points(log(1:20),log(mvol.2),col=4)
points(log(1:20),log(mvol.3),col=3)
points(log(1:20),log(mvol.4),col=2)
#points(log(1:20),log(mvol.5),col=6)

lines(log(1:20),log(fm.1),col=2)
lines(log(1:20),log(fm.2),col=2)
lines(log(1:20),log(fm.3),col=2)
lines(log(1:20),log(fm.4),col=2)
lines(log(1:20),log(fm.5),col=2)
#In the plot above it looks like the the best fit of the fbm to the volatility process is for q equal to 2

#best fit of fbm to volatility process
plot(log(1:20),log(mvol.3))
lines(log(1:20),log(fm.3),col=2)
lm(log(mvol.3)~log(1:20)) # -0.7540       0.2961
lm(log(fm.3)~log(1:20)) #-0.7564       0.2914 
 

#fit of Heston model for q=2 
plot(log(1:20),log(hm.3))
plot(log(1:20),log(mvol.3))
lm(log(hm)~log(1:20)) # -11.232        1.022

#Heston model has a much steeper slope than the fbm and the intercept is also much more off than the fbm.
#the scaling of the slope using Heston model is approximately 3.35*H

#plotted together (not the best)
plot(log(1:20),log(mvol.3),ylim=c(-12,2))
lines(log(1:20),log(fm.3),col='red')
lines(log(1:20),log(hm.3),col='blue')

######## slope and concavity effect ################################
q=c(1,1.5,2,2.5,3)
fm.Hq.points = 0.15*q
slopes=c(0.1449,0.2192,0.2961,0.3761,0.4594)

plot(q,fm.Hq.points)
points(q,slopes,col=2)
lines(q,fm.Hq.points,col='blue')
lines(q,slopes,col='red')

hm.Hq.points=0.5*q
plot(q,hm.Hq.points)

#plotted togetger (looks horrible)
plot(q,fm.Hq.points,ylim=c(0.2,1.8))
points(q,slopes,col=2)
lines(q,fm.Hq.points,col='blue')
lines(q,slopes,col='red')
lines(q,hm.Hq.points)


#h øger slope og multiplikatoren øger intercept i fBm

#### volatility skew and volatility smile #####



