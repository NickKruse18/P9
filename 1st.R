


plot((331:720)/60+4,Asset[331:720],type="l",xlab = "hour",ylab="price",main="SPY on the 29. of August 2022")
plot(((331:720)/60+4)[D[1951:2340,1]!=0],D[1951:2340,1][D[1951:2340,1]!=0],ylim=c(2,7),col="green",type="l",xlab = "hour",ylab="price",main="Call Options for SPY with maturity on the 31. of August 2022")
lines(((331:720)/60+4)[D[3901:4290,1]!=0],D[3901:4290,1][D[3901:4290,1]!=0],col="red")
legend(13, 5, legend=c("K = 400","K = 405"), col=3:2,lty=1)

set.seed(2)

t = 365
n = 10000
H = 0.1
rho = 0.9
S0 = mean(Asset[331:720])
S = SampleOU(S0,-2,0,0,1,0.15,-0.9,10000,t,Tau)
S = SampleOU(S0,-2,0,0,1.23,0.25,-0.9,10000,t,Tau)
S = SampleOU(S0,-1,0,0,1.2,0.22,-0.9,10000,t,Tau)
S = SampleHes(1,0.02,4,0.02,0.7,-0.9,n,1000,Tau)
S = SampleOU(S0,-1.689767327,0,0,1.155426937,0.22,-0.9,n,t,Tau)
set.seed(2)
S = SampleOU(S0,model[1],model[2],model[3],model[4],model[5],model[6],n,t,Tau)
hist(S[100,])
Tau = 1:100
K = S0*exp((-2:2)*0.001)
P = matrix(0,length(Tau),length(K))

for(i in 1:length(Tau)){ P[i,] = Price(S[i,],K) }
V = ImpliedVol(S0,Tau,K,P,0,t)

FitPowerLaw(Tau/t,abs((V[,5]-V[,1]))/(log(K[5])-log(K[1])),"Skew for an RFSV Model")[[1]]

FitPowerLaw(Tau/t,abs((V[,5]-V[,1]))/(log(K[5])-log(K[1])),"Skew for Fitted Model for the 28. of January 2022")[[1]]

as.Date("2022-01-28")+Tau[c(3,9,16,20)]

Tau[c(3,9,16,20)]/365

plot(log(K/S0),V.Data[,3],xlab = "log-moneyness",ylab = "implied volatility",main = "Smile for the 5. of February 2022,  tau = 0.020")
lines(log(K/S0),V[3,],col="red")
plot(log(K/S0),V.Data[,9],xlab = "log-moneyness",ylab = "implied volatility",main = "Smile for the 18. of February 2022,  tau = 0.058")
lines(log(K/S0),V[9,],col="red")
plot(log(K[-c(1,2,3,4,7,8,9,10)]/S0),V.Data[-c(1,2,3,4,7,8,9,10),16],xlab = "log-moneyness",ylab = "implied volatility",main = "Smile for the 7. of March 2022,  tau = 0.104")
lines(log(K/S0),V[16,],col="red")
plot(log(K[-c(1,2)]/S0),V.Data[-c(1,2),20],xlab = "log-moneyness",ylab = "implied volatility",main = "Smile for the 4. of April 2022,  tau = 0.208")
lines(log(K/S0),V[20,],col="red")


plot(log(K/S0),V[3,],xlab = "log-moneyness",ylab = "implied volatility",main = "Model Smile for the 6. of September 2022,  tau = 0.022",type = "l",col = "red")
plot(log(K/S0),V[9,],xlab = "log-moneyness",ylab = "implied volatility",main = "Model Smile for the 19. of September 2022,  tau = 0.058",type = "l",col = "red")
plot(log(K/S0),V[16,],xlab = "log-moneyness",ylab = "implied volatility",main = "Model Smile for the 7. of October 2022,  tau = 0.107",type = "l",col = "red")
plot(log(K/S0),V[18,],xlab = "log-moneyness",ylab = "implied volatility",main = "Model Smile for the 18. of November 2022,  tau = 0.222",type = "l",col = "red")
plot(log(K[-1]/S0),V.Data[-1,3],xlab = "log-moneyness",ylab = "implied volatility",main = "Smile for the 6. of September 2022,  tau = 0.022")
lines(log(K/S0),V[3,],col="red")
plot(log(K[-1]/S0),V.Data[-1,9],xlab = "log-moneyness",ylab = "implied volatility",main = "Smile for the 19. of September 2022,  tau = 0.058")
lines(log(K/S0),V[9,],col="red")
plot(log(K[c(-1,-3,-5)]/S0),V.Data[c(-1,-3,-5),16],xlab = "log-moneyness",ylab = "implied volatility",main = "Smile for the 7. of October 2022,  tau = 0.107")
lines(log(K/S0),V[16,],col="red")
plot(log(K[c(-7,-8,-9,-10)]/S0),V.Data[c(-7,-8,-9,-10),18],xlab = "log-moneyness",ylab = "implied volatility",main = "Smile for the 18. of November 2022,  tau = 0.222")
lines(log(K/S0),V[18,],col="red")
model = c(-1.7413394,0.1505670,1.3808247,1.0972277,0.2441113,-0.8499624)
FitPowerLaw(Tau[-c(19:21)]/t,abs((V[-c(19:21),5]-V[-c(19:21),1]))/(log(K[5])-log(K[1])),"Skew for The Fitted Model")[[1]]
X = Skew(D[(5*390+1):(6*390),],D[(10*390+1):(11*390),],Asset[331:720],Tau,t,K.Data[c(6,11)])

FitPowerLaw(Tau[-c(16,17,24)]/365,abs(X[-c(16,17,24)]),"Approximate Volatility Skew for the 28. of January 2022")[[1]]
#FitPowerLaw(Tau[-13]/365,abs(X[-13]),"Approximate Volatility Skew for the 29. of August 2022")[[1]]

model = FitModel(V.Data,model,Tau,K.Data,S0,10000,t,c(0.6033132, -0.2483098))
#model = FitModel(V.Data,model,Tau,K.Data,S0,10000,t,c(0.3719807,-0.2821460))
model

MC.Error = PriceError(model[1],model[2],model[3],model[4],model[5],model[6],t)

MC.Error

plot(log(MC.Error[[1]]),MC.Error[[4]],xlab="log repetitions", ylab = "error",main = "Monte Carlo Error in IV for Various Repetitions")
lines(log(MC.Error[[1]]),0.58*MC.Error[[1]]^-0.56,col="red")


V.Data = EstimateVol(D,Asset[331:720],Tau,K.Data,365)
persp(log(K.Data),Tau,V.Data,phi=30,theta = 120)

K.Data = c(405,410,415,420,425,430,431,432,433,434,435,440,445,450,455,460)

weekdays(as.Date("2022-01-29"))
D = GetData("SPY",as.Date("2022-01-28"),180,0)
Asset = D[[1]]
Tau = D[[2]]
D = D[[3]]

#model jan 28: c(-1.705766673,-0.008114986,1.523440189,1.752889018,0.336198337,-0.926887688)
#K.Data = c(405,410,415,420,425,430,431,432,433,434,435,440,445,450,455,460)
#model aug 29: c(-1.7413394,0.1505670,1.3808247,1.0972277,0.2441113,-0.8499624)
#K.Data = c(375,380,385,390,395,400,401,402,403,404,405,410,415,420,425,430)
Tau = readRDS("SPYAug29")[[2]]
saveRDS(list(Asset,Tau,D),"SPYJan28")
D
#Simulation Errors

#Fractional Brownian motion

fBm.Error = ComparefBm(0.1,1000,c(10,20,50,100,200,500,1000,2000,5000,10000))
fBm.Error[4,] = fBm.Error[4,]/sqrt(c(10,20,50,100,200,500,1000,2000,5000,10000))
plot(log(fBm.Error[1,]),fBm.Error[2,],type = "l",xlab="log steps",ylab="time, seconds",main="Computation time for fBM Simulations")
lines(log(fBm.Error[1,]),fBm.Error[4,],col="red")
plot(log(fBm.Error[1,]),fBm.Error[3,],type = "l",xlab="log steps",ylab="RMSE",main="RMSE for fBM Simulations")
lines(log(fBm.Error[1,]),fBm.Error[5,],col="red")
fBm.Error

#Ornstein-Uhlenbeck

OU.Error = CompareOU(0,1,1,1,0.1,100000,1000,c(10,20,50,100,200,500,1000,2000,5000,10000))
OU.Error
lm(log(OU.Error[2,-1])~log(OU.Error[1,-1]))
lm(log(OU.Error[4,-1])~log(OU.Error[1,-1]))
lm(log(OU.Error[3,-1])~log(OU.Error[1,-1]))
lm(log(OU.Error[5,-1])~log(OU.Error[1,-1]))
OU.Error[2,] = OU.Error[2,]^2
plot(log(OU.Error[1,-1]),OU.Error[2,-1],type="l",xlab="log steps",ylab="time, seconds",main="Computation time for fOU Approximations")
lines(log(OU.Error[1,-1]),OU.Error[4,-1],col="red")
plot(log(OU.Error[1,-1]),OU.Error[3,-1],type="l",xlab="log steps",ylab="RMSE",main="RMSE for fOU Approximations")
lines(log(OU.Error[1,-1]),OU.Error[5,-1],col="red")

#Rough Fractional Volatility

MC.Sim.Error = PriceError(-1,0,0,1,0.1,0,512)
MC.Sim.Error
plot(log(MC.Sim.Error[[1]]),MC.Sim.Error[[4]],xlab="log repetitions", ylab = "error",main = "Monte Carlo Error in Price for Various Repetitions")
lines(log(MC.Sim.Error[[1]]),0.32*MC.Sim.Error[[1]]^-0.44,col="red")


plot(log(MC.Sim.Error[[1]]),MC.Sim.Error[[5]],xlab="log repetitions", ylab = "time",main = "Computation time for Monte Carlo estimate for various Repetitions")
lines(log(MC.Sim.Error[[1]]),0.001157*MC.Sim.Error[[1]],col="red")

CompareRFSV(1,-2,0,0,1,0.1,0,c(256),1000)
CompareRFSV(1,-2,0,0,1,0.1,0.9,c(512,256,128,64),10000)
RFSV.Error = CompareRFSV(1,-1,0,0,1,0.1,0,c(4096,2048,1024,512,256,128,64),10000)
plot(RFSV.Error[2,])
plot(7:1,RFSV.Error[1,],xlab="log steps",ylab="time, seconds",main = "Computation time for RFSV Simulations",type="l")


legend(0.8, 0.2, legend=c("n=4096","n=2048","n=1024","n=512","n=256","n=128","n=64"), col=1:7,lty=1)
N256 = MCDist(-1,0,0,1,0.1,0,256)
N2048 = MCDist(-1,0,0,1,0.1,0,2048)
N4096 = MCDist(-1,0,0,1,0.1,0,4096)
hist(N256,30,xlab="MC price",main="Histogram of 100 Monte Carlo RFSV, n = 256")
hist(N2048,30,xlab="MC price",main="Histogram of 100 Monte Carlo RFSV, n = 2048")
hist(N4096,50,xlab="MC price",main="Histogram of 100 Monte Carlo RFSV, n = 4096")



BSPrice = function(S,K,tau,sigma){
  dp = (log(S/K)+sigma^2*tau/2)/(sigma*sqrt(tau));  dn = dp - sigma*sqrt(tau)
  return(pnorm(dp)*S-pnorm(dn)*K)
}

BSDelta = function(S,K,tau,sigma){
  dp = (log(S/K)+sigma^2*tau/2)/(sigma*sqrt(tau))
  return(pnorm(dp))
}

Hedge = function(S0,model,K,Tau,t,n,Data){
  Delta = matrix(0,length(Tau),length(K))
  set.seed(2)
  S = SampleOU(1,model[1],model[2],model[3],model[4],model[5],model[6],n,t,Tau)
  for(i in 1:length(Tau)){ 
    for(j in 1:length(K)){
      Delta[i,j] = (Price(K[j]*S[i,]*1.0001,S0)-Price(K[j]*S[i,]*0.9999,S0))/(K[j]*0.0002)
    }
  }
  S = S0*S
  Hist = -Price(S[length(Tau),],S0)
  Hedge = numeric(n)
  S = SampleOU(S0,model[1],model[2],model[3],model[4],model[5],model[6],n,t,Tau)
  Hist = Hist + (S[length(Tau),]-S0)*(S[length(Tau),]>S0)
  for(i in 2:length(Tau)){ Hedge = Hedge + (S[i-1,]-S[i,])*Delta[length(Tau)-i+2,round(S[i-1,])] }
  return(list(Delta,Hist+Hedge,Hist))
}
K.Data
Delta = Hedge(50,model,1:400,1:144,365,10000)
Delta[[2]] = 405/50*Delta[[2]]
Delta[[3]] = 405/50*Delta[[3]]

Delta[[2]]
hist(Delta[[2]],xlab="return",main="Histogram for Delta Hedging")
lines(c(mean(Delta[[2]]),mean(Delta[[2]])),c(0,10000),col="red")
lines(c(mean(Delta[[2]])+sqrt(var(Delta[[2]])),mean(Delta[[2]])+sqrt(var(Delta[[2]]))),c(0,10000),col="blue")
lines(c(mean(Delta[[2]])-sqrt(var(Delta[[2]])),mean(Delta[[2]])-sqrt(var(Delta[[2]]))),c(0,10000),col="blue")
lines(c(sort(Delta[[2]])[9750],sort(Delta[[2]])[9750]),c(0,10000))
lines(c(sort(Delta[[2]])[250],sort(Delta[[2]])[250]),c(0,10000))

hist(Delta[[3]],xlab="return",main="Histogram without Delta Hedging")
lines(c(mean(Delta[[3]]),mean(Delta[[3]])),c(0,10000),col="red")
lines(c(mean(Delta[[3]])+sqrt(var(Delta[[3]])),mean(Delta[[3]])+sqrt(var(Delta[[3]]))),c(0,10000),col="blue")
lines(c(mean(Delta[[3]])-sqrt(var(Delta[[3]])),mean(Delta[[3]])-sqrt(var(Delta[[3]]))),c(0,10000),col="blue")
lines(c(sort(Delta[[3]])[9750],sort(Delta[[3]])[9750]),c(0,10000))
lines(c(sort(Delta[[3]])[250],sort(Delta[[3]])[250]),c(0,10000))
sqrt(var(Delta[[2]]))
sqrt(var(Delta[[3]]))
sort(Delta[[2]])[250]
sort(Delta[[2]])[9750]
sort(Delta[[3]])[250]
sort(Delta[[3]])[9750]
mean(Delta[[2]])
mean(Delta[[3]])
persp(1:144,370:430,Hedge(405,model,370:430,1:144,365,10000),theta=30,phi=30,xlab="tau",ylab="S",zlab="Delta",main="Delta for K = 405 on the 29. of August")


