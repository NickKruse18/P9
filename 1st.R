
set.seed(2)

t = 365
n = 10000
H = 0.1
rho = 0.9
S0 = mean(Asset[331:720])
S = SampleOU(S0,-1,0,0,1.2,0.22,-0.9,n,t,Tau)
S = SampleOU(S0,-1.689767327,0,0,1.155426937,0.22,-0.9,n,t,Tau)
set.seed(2)
S = SampleOU(S0,model[1],model[2],model[3],model[4],model[5],model[6],n,t,Tau)

P = matrix(0,length(Tau),5)
K = S0*exp((-2:2)*0.001)

for(i in 1:length(Tau)){
  P[i,] = Price(S[i,],K)
}
V = ImpliedVol(S0,Tau,K,P,0,t)
Tau/365
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
FitPowerLaw(Tau/t,abs((V[,5]-V[,1]))/(log(K[5])-log(K[1])),"Skew for The Fitted Model")[[1]]
X = Skew(D[(5*390+1):(6*390),],D[(10*390+1):(11*390),],Asset[331:720],Tau,t,c(400,405))
FitPowerLaw(Tau[-13]/365,abs(X[-13]),"Approximate Volatility Skew for the 29. of August 2022")[[1]]

model = FitModel(V.Data,model,Tau,K.Data,S0,1000,t,c(0.3719807,-0.2821460))


MC.Error = PriceError(model[1],model[2],model[3],model[4],model[5],model[6],t)
model
MC.Error

plot(log(MC.Error[[1]]),MC.Error[[4]],xlab="log repetitions", ylab = "error",main = "Monte Carlo error for various Repetitions")
lines(log(MC.Error[[1]]),0.58*MC.Error[[1]]^-0.56,col="red")


V.Data = EstimateVol(D,Asset[331:720],Tau,K.Data,365)
persp(log(K.Data),Tau,V.Data,phi=30,theta = 120)

K.Data = c(375,380,385,390,395,400,401,402,403,404,405,410,415,420,425,430)

as.Date("2022-08-29")+Tau
D = GetData("SPY",as.Date("2022-08-29"),180,0)
Asset = D[[1]]
Tau = D[[2]]
D = D[[3]]

saveRDS(list(Asset,Tau,D),"SPYAug29")

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

MC.Sim.Error = PriceError(-1,0,0,1,0.1,0,1024)

RFSV.Error = CompareRFSV(1,-1,0,0,1,0.1,0,c(4096,2048,1024,512,256,128,64),1000)
plot(RFSV.Error[2,])
plot(RFSV.Error[1,])
RFSV.Error















