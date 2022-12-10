
library("somebm")
library("rumidas")

data=rv5
plot(data,type="l",main="S&P 500 realized variance at 5-minutes",xlab="Days",ylab="Realized variance")
#the high realized variance around 2200-2300 days was in the middle of the financial crisis.
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

plot(fm.3)
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
    V[i] = abs(V[i-1] + sigma*sqrt(V[i-1])*rnorm(1)/sqrt(n))
  }
  return(V)
}

hv = HestonVol(1,3,20000)
hm.1 = m.delta(1,hv,1:20)
hm.2 = m.delta(1.5,hv,1:20)
hm.3 = m.delta(2,hv,1:20)
hm.4 = m.delta(2.5,hv,1:20)
hm.5 = m.delta(3,hv,1:20)

#The Heston model is very non-consistent and do changes a lot when adjusting the paremeters.

########################################### Plots and fitting ############################################
## The fBm model
plot(log(1:20),log(mvol.5),main="The log-volatility processes with fitted fBm",ylim=c(-0.8,0.8),  xlab = expression(log(Delta)), ylab = expression(log(m(q,Delta))))

points(log(1:20),log(mvol.1),col="purple")
points(log(1:20),log(mvol.2),col="blue")
points(log(1:20),log(mvol.3),col="red")
points(log(1:20),log(mvol.4),col="green")
#points(log(1:20),log(mvol.5),col=6)

#lines(log(1:20),log(1:20)*0.219-0.7398)
#lines(log(1:20),log(1:20)*0.2961-0.7540)
#lines(log(1:20),log(1:20)*0.376-0.695)
#lines(log(1:20),log(1:20)*0.459-0.576)

lines(log(1:20),log(fm.1),col="purple")
lines(log(1:20),log(fm.2),col="blue")
lines(log(1:20),log(fm.3),col="red")
lines(log(1:20),log(fm.4),col="green")
lines(log(1:20),log(fm.5),col="black")

legend("topleft", inset=.05, title="q values", c("q=1","q=1.5","q=2", "q=2.5", "q=3"), 
       col=c("purple","blue","red","green","black"),  lty=1, cex=0.8)

## Heston model
plot(log(1:20),log(mvol.5),main="The log-volatility processes and fitted Heston model", ylim=c(-0.8,0.8),  xlab = expression(log(Delta)), ylab = expression(log(m(q,Delta))))

points(log(1:20),log(mvol.1),col="purple")
points(log(1:20),log(mvol.2),col="blue")
points(log(1:20),log(mvol.3),col="red")
points(log(1:20),log(mvol.4),col="green")
#points(log(1:20),log(mvol.5),col=6)

lines(log(1:20),log(hm.1),col="purple")
lines(log(1:20),log(hm.2),col="blue")
lines(log(1:20),log(hm.3),col="red")
lines(log(1:20),log(hm.4),col="green")
lines(log(1:20),log(hm.5),col="black")

legend("topleft", inset=.05, title="q values", c("q=1","q=1.5","q=2", "q=2.5", "q=3"), 
       col=c("purple","blue","red","green","black"),  lty=1, cex=0.8)


#In the plot above it looks like the the best fit of the fbm to the volatility process is for q equal to 2

########## LM test ##################
mlm1=lm(log(mvol.1)~log(1:20))$coefficients # -0.6352747   0.1449026 
mlm2=lm(log(mvol.2)~log(1:20))$coefficients # -0.7398497   0.2192481 
mlm3=lm(log(mvol.3)~log(1:20))$coefficients # -0.7539825   0.2961300 
mlm4=lm(log(mvol.4)~log(1:20))$coefficients # -0.6954881   0.3760548 
mlm5=lm(log(mvol.5)~log(1:20))$coefficients # -0.5760885   0.4594152

lm(log(fm.1)~log(1:20))$coefficients # -0.6046370   0.1455551
lm(log(fm.2)~log(1:20))$coefficients # -0.7182084   0.2183343
lm(log(fm.3)~log(1:20))$coefficients # -0.7563779   0.2913825
lm(log(fm.4)~log(1:20))$coefficients # -0.7361408   0.3646410 
lm(log(fm.5)~log(1:20))$coefficients # -0.6683289   0.4379853

lm(log(hm.1)~log(1:20))$coefficients # -1.7451924   0.3543493
lm(log(hm.2)~log(1:20))$coefficients # -1.7242139   0.4586978
lm(log(hm.3)~log(1:20))$coefficients # -1.4654457   0.5440321 
lm(log(hm.4)~log(1:20))$coefficients # -1.0612765   0.6150094 
lm(log(hm.5)~log(1:20))$coefficients # -0.5487761   0.6734308


################ Plot for best fit which is q=2 #################
#best fit of fbm to volatility process
plot(log(1:20),log(mvol.3))
lines(log(1:20),log(fm.3),col=2)

 

#fit of Heston model for q=2 
plot(log(1:20),log(hm.3))
lines(log(1:20),log(mvol.3))
lm(log(hm)~log(1:20)) # -11.232        1.022

#Heston model has a much steeper slope than the fbm and the intercept is also much more off than the fbm.
#the scaling of the slope using Heston model is approximately 3.35*H

#plotted together (not the best)
plot(log(1:20),log(mvol.3), main="Log-volatility process for q=2",xlab=expression(log(Delta)),ylab=expression(log(m(q,Delta))))
lines(log(1:20),log(fm.3),col='red')
lines(log(1:20),log(hm.3),col='blue')
legend("topleft", inset=.05, title="model", c("fBm","Heston"), 
       col=c("red","blue"),lty=1, cex=0.8)
#We cna only change the Hetson model in the way it cus the vertical axis. We cannot change or manipulate with its slope.

######## slope and concavity effect ################################
q=c(1,1.5,2,2.5,3)
fm.Hq.points =0.15*q
slopes=c(0.1449,0.2192,0.2961,0.3761,0.4594)

lm(slopes~q)

plot(q,slopes,main="Slope of line for log-volatility process",ylab=expression(zeta))
points(q,slopes,col='red')
lines(q,-0.01522+0.15*q, col='blue') #just a try 
lines(q,slopes,col='red')

plot(q,slopes,main="Slope of line for log-volatility and fBm",ylab = "slope")
points(q,fm.Hq.points,col='blue')
lines(q,fm.Hq.points,col='blue')
lines(q,slopes,col='red')

hm.Hq.points=0.5*q
plot(q,hm.Hq.points,main="Slope of line for log-volatility and Heston",ylab="slope",ylim = c(0.15,1.5))
points(q,slopes,col='red')
points(q,hm.Hq.points,col='blue')
lines(q,hm.Hq.points,col='blue')
lines(q,slopes,col='red')



#plotted togetger (looks horrible)
plot(q,fm.Hq.points,ylim=c(0.1,1.5))
points(q,slopes,col=2)
lines(q,fm.Hq.points,col='blue')
lines(q,slopes,col='red')
lines(q,hm.Hq.points)


#Hurst parameter øger slope og multiplikatoren øger intercept i fBm



plot(fbm(0.15,1000))

