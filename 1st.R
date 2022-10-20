#Dette er hoved filen hvor alle ting bliver udregnet.
#Denne fil benytter funktionerne skrevet i 2nd.R, så hele den fil skal være loadet ind i memory før denne fil virker.
#Tryk ctrl + A og ctrl + enter, i 2nd.R for at loade det hele.
#Efter 2nd.R er loaded er det ikke nødvendigt at kigge på den mere.


#---------------------------------------------------------------------------------------------------------------

#Simuler rBergomi med Euler Maruyama og Hybrid scheme

t = 400 # mængden af skridt i simuleringen (mere akkurat simulering) 
n = 10000 # mængden af simuleringer (mere akkurat forventet værdi E[(S-K)^+])
eta = 0.4 # rBergomi parameter der øger volatiliteten af volatiliteten
H = 0.1 # rBergomi parameter ændrer powerlawen på volatity skew: tau^(H-1/2)
rho = 0.85 # korrelationen mellem den brownske bevægelse Z der hører til prisen S og den brownske bevægelse der høre til volatiliteten Xi
sigma0 = 0.235^2 #Begyndelses volatiliteten
S0 = 1 #Begyndelses pris


B = rBerVol2(sigma0,eta,H,t,n,3)

dW = B[[2]]
B = B[[1]]
plot(B%*%rep(1,n)/n)
dZ = rho*dW+sqrt(1-rho^2)*dBM(t,n)

S = BSV(S0,-0.327,1/t,B,dZ) # simulerer n paths af prisen S for en rBergomi model. kompexitet af O(n t)
mean(S%*%rep(1,n)/n) # et problem jeg ikke har løst er at S har en smule uønsket bias afhængigt af parameterne. så -0.263 i BSV skal justeres indtil mean er lig 1
plot(S%*%rep(1,n)/n,type = "l") # gennemsnits path for S_tau

Tau = (1:20)*(t/20) # de expiration tider for de simulerede options
K = S0*exp((-50:50)*0.01) # de strike priser for de simulerede options

PSurface = PriceSurface(S,Tau,K) #Pris overfladen for forskellige options med rBergomi som underlying. Tau er expiration tider og K er strike priser
persp(Tau,log(K/S0),PSurface,theta=120,phi=30,xlab="tau",ylab="k",zlab="P",main="Monte Carlo Price Surface") #3D plot af pris overfalden. theta og phi kontrollerer syns vinklen. log er taget til K for at få log moneyness
VolSurface = ImpliedVol(S0,(1:20)/20,K,PSurface,0) #implied volatilitets overflade for option priserne.
persp(Tau,log(K/S0),VolSurface,theta=120,phi=30,xlab="tau",ylab="k",zlab="V",main="Monte Carlo Implied Volatility Surface") #3D plot af volatilitets overfalden. theta og phi kontrollerer syns vinklen. log er taget til K for at få log moneyness
FitPowerLaw(VolSurface/0.01,50,51,Tau) # volatility skewet of log(K/S0) = 0 og det bedste powerlaw fit (a tau^b). (50,51 er valgt til at matche den søjle i VolSurface der svarer til log(K/S) = 0 samt den lige ved siden af.)
lines(Tau,rho*eta/2/(H+1/2)/(H+3/2)*sqrt(2*H)*(Tau/t)^(H-0.5),col="red")
# Den røde linje er teoretisk volatility skew og de sorte prikker er simuleret. De matcher næsten.


#---------------------------------------------------------------------------------------------------------------


# udregner den implied volatility surface og skew for en given aktie.

Asset = "^SPX" # Det underliggende aktiv for optionerne. Apple i dette tilfælde
Source = "yahoo"

Data = AssetVolatility(Asset,Source,3712,as.Date("2022-10-19")) #3712 og as.Date skal være pris og dato for tidspunktet dataen var taget
S0 = Data[[3]]
Skew = Data[[2]]
Data = Data[[1]]
plot(Skew) 

# Volatility Smile
for(i in Skew[,1]){
  plot(Data[(Data[,1] == i),2],Data[(Data[,1] == i),4],main = i)
}

SData = Data[(Data[,5]>=100)&(Data[,1]<=800),] # Fjerner data med mindre end 100 volumen og mere end 800 dage til maturity
Fit = FitModel(c(-1.36,0.7,0,0,0.31),S0,800,1000,SData,"rBer") # Fit med rBergomi model
Fit.BS = FitModel(c(-1.36,0.02),S0,800,1000,SData,"BS") # Fit med Black-Scholes model
#SPX (old) rBer 0.002844473. BS 0.005176576
#AAPL rBer 0.003353342, BS 0.005057349. Resultat for AAPL data Oct 19
#SPX rBer 0.002139801, BS 0.003452251, pris 3712. Resultat for ^SPX data Oct 19

Data = readRDS("SPXData")
Fit = readRDS("SPXFit")
Fit.BS = readRDS("SPXFitBS")

Tau = (1:20)*40 # Maturities hvor prisoverfladen er evalueret
K = S0*exp((-50:50)*0.01) # Strikes hvor prisoverfladen er evalueret
# Prisoverflade for rBergomi modellen
PSurface = rBerSurface(c(exp(Fit[[1]][1]),Fit[[1]][2],0.5/(1+exp(Fit[[1]][3])),2*1/(1+exp(Fit[[1]][4]))-1,Fit[[1]][5],S0),800,1000,Tau,K)
# Prisoverflade for Black-Scholes modellen
PSurface.BS = BSSurface(c(exp(Fit.BS[[1]][1]),Fit.BS[[1]][2],S0),800,1000,Tau,K)
# Prisoverflade i rød og faktiske priser (Cleaned) i sort for rBergomi
points(trans3d(Fit[[3]][,1],log(Fit[[3]][,2]/S0),Fit[[3]][,3],pmat = persp(Tau,log(K/S0),PSurface,theta=120,phi=30,xlab="tau",ylab="k",zlab="P",main="rBergomi Price Surface (Cleaned Data)",border="red")))
# Prisoverflade i rød og faktiske priser (Cleaned) i sort for Black-Scholes
points(trans3d(Fit.BS[[3]][,1],log(Fit.BS[[3]][,2]/S0),Fit.BS[[3]][,3],pmat = persp(Tau,log(K/S0),PSurface.BS,theta=120,phi=30,xlab="tau",ylab="k",zlab="P",main="Black-Scholes Price Surface (Cleaned Data)",border="red")))
# Prisoverflade i rød og alle faktiske priser (SData) i sort for rBergomi
points(trans3d(SData[,1],log(SData[,2]/S0),SData[,3],pmat = persp(Tau,log(K/S0),PSurface,theta=120,phi=30,xlab="tau",ylab="k",zlab="P",main="rBergomi Price Surface (Uncleaned Data)",border="red")))


#---------------------------------------------------------------------------------------------------------------
#Heston Model
t = 400
n = 100000
eta = 0.4
H = 0.1
rho = 0.85
sigma0 = 0.235^2
S0 = 1


dZ = dBM(t,n)
dW = rho*dZ+sqrt(1-rho^2)*dBM(t,n)

H = HesVol(sigma0,eta,dW)

plot(H%*%rep(1,n)/n,type = "l")

S = BSV(S0,-0.168,1/t,H,dZ)
mean(S%*%rep(1,n)/n)
plot(S%*%rep(1,n)/n,type = "l")

Tau = (1:20)*(t/20)
K = S0*exp((-50:50)*0.01)

PSurface = PriceSurface(S,Tau,K)
persp(Tau,log(K/S0),PSurface,theta=120,phi=30,xlab="tau",ylab="k",zlab="P",main="Heston Model Price Surface")
VolSurface = ImpliedVol(S0,(1:20)/20,K,PSurface,0)
persp(Tau,log(K/S0),VolSurface,theta=120,phi=30,xlab="tau",ylab="k",zlab="V",main="Heston Model Implied Volatility Surface")
FitPowerLaw(VolSurface/0.01,50,51,Tau,"Heston Model Volatility Skew")

#---------------------------------------------------------------------------------------------------------------
#Delta Hedging

t = 800 # mængden af skridt i simuleringen (mere akkurat simulering) 
n = 1000 # mængden af simuleringer (mere akkurat forventet værdi E[(S-K)^+])
eta = 0.4 # rBergomi parameter der øger volatiliteten af volatiliteten
H = 0.1 # rBergomi parameter ændrer powerlawen på volatity skew: tau^(H-1/2)
rho = 0.85 # korrelationen mellem den brownske bevægelse Z der hører til prisen S og den brownske bevægelse der høre til volatiliteten Xi
sigma0 = 0.235^2 #Begyndelses volatiliteten
S0 = 1 #Begyndelses pris

dW = dBM(t,n)
dZ = dBM(t,n)
S = BSV(S0,0,1/t,matrix(sigma0,t,n),dZ)

V = Hedging(S,1,S,c())
sort(V)[n*0.95]
sort(V)[n*0.05]
hist(V)
mean(V)
var(V)

#---------------------------------------------------------------------------------------------------------------
#Benchmark

t = 400 # mængden af skridt i simuleringen (mere akkurat simulering) 
n = 1000 # mængden af simuleringer (mere akkurat forventet værdi E[(S-K)^+])
eta = 0.4 # rBergomi parameter der øger volatiliteten af volatiliteten
H = 0.1 # rBergomi parameter ændrer powerlawen på volatity skew: tau^(H-1/2)
rho = 0.85 # korrelationen mellem den brownske bevægelse Z der hører til prisen S og den brownske bevægelse der høre til volatiliteten Xi
sigma0 = 0.235^2 #Begyndelses volatiliteten
S0 = 1 #Begyndelses pris

mark = microbenchmark("HS t = 100" = {rBerVol2(sigma0,eta,H,100,n,3)},"EM t = 100" = {rBerVol(sigma0,eta,H,100,n)},
                      "HS t = 400" = {rBerVol2(sigma0,eta,H,400,n,3)},"EM t = 400" = {rBerVol(sigma0,eta,H,400,n)},
                      "HS t = 1000" = {rBerVol2(sigma0,eta,H,1000,n,3)},"EM t = 1000" = {rBerVol(sigma0,eta,H,1000,n)})
microbenchmark("HS t = 100" = {rBerVol2(sigma0,eta,H,100,n,3)},
               "HS t = 400" = {rBerVol2(sigma0,eta,H,400,n,3)},
               "HS t = 1000" = {rBerVol2(sigma0,eta,H,1000,n,3)})
mark

#---------------------------------------------------------------------------------------------------------------

#Simuler rBergomi med Euler Maruyama

t = 400 # mængden af skridt i simuleringen (mere akkurat simulering) 
n = 100 # mængden af simuleringer (mere akkurat forventet værdi E[(S-K)^+])
eta = 0.4 # rBergomi parameter der øger volatiliteten af volatiliteten
H = 0.1 # rBergomi parameter ændrer powerlawen på volatity skew: tau^(H-1/2)
rho = 0.85 # korrelationen mellem den brownske bevægelse Z der hører til prisen S og den brownske bevægelse der høre til volatiliteten Xi
sigma0 = 0.235^2 #Begyndelses volatiliteten
S0 = 1 #Begyndelses pris




b = rBerVol(sigma0,eta,H,t,n) # simulerer n paths af volatiliten Xi for en rBergomi model. Meste af processer tiden: kompexitet af O(n t^2)
B = b[[1]] # Xi_tau(tau) for tau = 1/t, 2/t, ..., (t-1)/t, 1
dW = b[[3]]
b = b[[2]] # Xi_tau(1)
dZ = rho*dW+sqrt(1-rho^2)*dBM(t,n) # samme som dZ, dog med en korrelation til dZ på rho
plot(B%*%rep(1,n)/n,type = "l") # gennemsnits path for Xi_tau(tau)
plot(b%*%rep(1,n)/n,type = "l") # gennemsnits path for Xi_tau(1)
S = BSV(S0,-0.198,1/t,B,dZ) # simulerer n paths af prisen S for en rBergomi model. kompexitet af O(n t)
mean(S%*%rep(1,n)/n) # et problem jeg ikke har løst er at S har en smule uønsket bias afhængigt af parameterne. så -0.263 i BSV skal justeres indtil mean er lig 1
plot(S%*%rep(1,n)/n,type = "l") # gennemsnits path for S_tau

plot(CoVT(diff(log(S)),diff(b)),xlab = expression(tau),ylab=expression(paste("E[d",x[tau],"d", xi[tau](u),"]/dt")),main=expression(paste("Monte Carlo vs Theoretical E[d",x[tau],"d", xi[tau](u),"]/dt"))) 
lines(rho*eta*sqrt(2*H)*sqrt(B[1:(t-1),]%*%rep(1,n)/n)*b[1:(t-1),]%*%rep(1,n)/n*((t-1:(t-1))/t)^(H-0.5)/t, col = "red")
# Et egenskaberne for rBergomi er at de sorte prikker (simuleret) skal match den røde linje (teoretisk), (ligningen kan ses på s. 898 ligningen over 7.1. Special case...)

Tau = (1:20)*(t/20) # de expiration tider for de simulerede options
K = S0*exp((-50:50)*0.01) # de strike priser for de simulerede options

PSurface = PriceSurface(S,Tau,K) #Pris overfladen for forskellige options med rBergomi som underlying. Tau er expiration tider og K er strike priser
persp(Tau,log(K/S0),PSurface,theta=120,phi=30,xlab="tau",ylab="k",zlab="P",main="Monte Carlo Price Surface") #3D plot af pris overfalden. theta og phi kontrollerer syns vinklen. log er taget til K for at få log moneyness
VolSurface = ImpliedVol(S0,(1:20)/20,K,PSurface,0) #implied volatilitets overflade for option priserne.
persp(Tau,log(K/S0),VolSurface,theta=120,phi=30,xlab="tau",ylab="k",zlab="V",main="Monte Carlo Implied Volatility Surface") #3D plot af volatilitets overfalden. theta og phi kontrollerer syns vinklen. log er taget til K for at få log moneyness
FitPowerLaw(VolSurface/0.01,50,51,Tau) # volatility skewet of log(K/S0) = 0 og det bedste powerlaw fit (a tau^b). (50,51 er valgt til at matche den søjle i VolSurface der svarer til log(K/S) = 0 samt den lige ved siden af.)
lines(Tau,rho*eta/2/(H+1/2)/(H+3/2)*sqrt(2*H)*(Tau/t)^(H-0.5),col="red")
# Den røde linje er teoretisk volatility skew og de sorte prikker er simuleret. De matcher næsten.


