#Dette er hoved filen hvor alle ting bliver udregnet.
#Denne fil benytter funktionerne skrevet i 2nd.R, så hele den fil skal være loadet ind i memory før denne fil virker.
#Tryk ctrl + A og ctrl + enter, i 2nd.R for at loade det hele.
#Efter 2nd.R er loaded er det ikke nødvendigt at kigge på den mere.

#---------------------------------------------------------------------------------------------------------------

#Simuler rBergomi med Euler Maruyama

t = 400 # mængden af skridt i simuleringen (mere akkurat simulering) 
n = 1000 # mængden af simuleringer (mere akkurat forventet værdi E[(S-K)^+])
eta = 0.4 # rBergomi parameter der øger volatiliteten af volatiliteten
H = 0.1 # rBergomi parameter ændrer powerlawen på volatity skew: tau^(H-1/2)
rho = 0.85 # korrelationen mellem den brownske bevægelse Z der hører til prisen S og den brownske bevægelse der høre til volatiliteten Xi
sigma0 = 0.235^2 #Begyndelses volatiliteten
S0 = 1 #Begyndelses pris


dZ = dBM(t,n) #t X n matrice af normalt fordelte værdier til at simulere Z_{tau/t}-Z_{(tau-1)/t}
dW = rho*dZ+sqrt(1-rho^2)*dBM(t,n) # samme som dZ, dog med en korrelation til dZ på rho

b = rBerVol(sigma0,eta,H,dW) # simulerer n paths af volatiliten Xi for en rBergomi model. Meste af processer tiden: kompexitet af O(n t^2)
B = b[[1]] # Xi_tau(tau) for tau = 1/t, 2/t, ..., (t-1)/t, 1
b = b[[2]] # Xi_tau(1)
plot(B%*%rep(1,n)/n,type = "l") # gennemsnits path for Xi_tau(tau)
plot(b%*%rep(1,n)/n,type = "l") # gennemsnits path for Xi_tau(1)
S = BSV(S0,-0.263,1/t,B,dZ) # simulerer n paths af prisen S for en rBergomi model. kompexitet af O(n t)
mean(S%*%rep(1,n)/n) # et problem jeg ikke har løst er at S har en smule uønsket bias afhængigt af parameterne. så -0.263 i BSV skal justeres indtil mean er lig 1
plot(S%*%rep(1,n)/n,type = "l") # gennemsnits path for S_tau

plot(CoVT(diff(log(S)),diff(b))) 
lines(rho*eta*sqrt(2*H)*sqrt(B[1:(t-1),]%*%rep(1,n)/n)*b[1:(t-1),]%*%rep(1,n)/n*((t-1:(t-1))/t)^(H-0.5)/t, col = "red")
# Et egenskaberne for rBergomi er at de sorte prikker (simuleret) skal match den røde linje (teoretisk), (ligningen kan ses på s. 898 ligningen over 7.1. Special case...)

Tau = (1:20)*(t/20) # de expiration tider for de simulerede options
K = S0*exp((-50:50)*0.01) # de strike priser for de simulerede options

PSurface = PriceSurface(S,Tau,K) #Pris overfladen for forskellige options med rBergomi som underlying. Tau er expiration tider og K er strike priser
persp(Tau,log(K/S0),PSurface,theta=120,phi=30) #3D plot af pris overfalden. theta og phi kontrollerer syns vinklen. log er taget til K for at få log moneyness
VolSurface = ImpliedVol(S0,(1:20)/20,K,PSurface,0) #implied volatilitets overflade for option priserne.
persp(Tau,log(K/S0),VolSurface,theta=120,phi=30) #3D plot af volatilitets overfalden. theta og phi kontrollerer syns vinklen. log er taget til K for at få log moneyness
FitPowerLaw(VolSurface/0.01,50,51,Tau) # volatility skewet of log(K/S0) = 0 og det bedste powerlaw fit (a tau^b). (50,51 er valgt til at matche den søjle i VolSurface der svarer til log(K/S) = 0 samt den lige ved siden af.)
lines(Tau,rho*eta/2/(H+1/2)/(H+3/2)*sqrt(2*H)*(Tau/t)^(H-0.5),col="red")
# Den røde linje er teoretisk volatility skew og de sorte prikker er simuleret. De matcher næsten.

#---------------------------------------------------------------------------------------------------------------

# udregner den implied volatility surface og skew for en given aktie.

Asset = "AAPL" # Det underliggende aktiv for optionerne. Apple i dette tilfælde

X = getSymbols(Asset,src="yahoo",auto.assign = F)[length(getSymbols(Asset,src="yahoo",auto.assign = F)[,1]),6] # handelsprisen og dato
date = index(X) # Isolere dato
X = X[[1]] # Isolere pris
D = getOptionChain(Asset,Exp=NULL,src="yahoo") # Hent optionerne for aktivet
r = 0.00000 # risk neutral rate


Data = ProcessOptions(D,date) # Konverterer det til en matrice med de nødvendige data, (dage til expiration, strike pris, og options pris)
Data = AddImpliedVol(Data,X,r) # Tilføjer implied vol. 

Surface = InterpolateSurface(Data[,1],Data[,2],Data[,3]) # Pris surface for optionerne. For at få en ordenlig surface er det nødvendigt at interpolere nogle værdier, så jeg lavede en hurtig funktion til det.
Surface = InterpolateSurface(Data[,1],Data[,2],Data[,4]) # Implied vol surface
FitPowerLaw(Surface[[3]],44,45,Surface[[1]]) # volatility skewet of log(K/S0) = 0 og det bedste powerlaw fit (a tau^b). (44,45 er valgt til at matche den søjle i VolSurface der svarer til log(K/S) = 0 samt den lige ved siden af. Se næste linjes Surface[[2]] for at finde de bedste indeks)
Surface[[2]] #Strike pris K for forskellige søjler i implied Vol surface

