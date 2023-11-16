## Set working directory
setwd("~/Gavin's Dropbox/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")
#setwd("~/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")

## Read in optimization functions
source("profitOptimTest.R")
source("RunDynamicOptimizationEnfTest.R")

## Enforcement Parameters
enfVariable = (23.3+180)*8
enfFixed = 60000

enfVec=seq(0,1,.1)
enfDet=c(0,.04,.25,.38,.47,.54,.59,.64,.68,.72,.75)

fine = 74.07  # $ / IUU trip
fine = fine / 63 # $/kg IUU @ x0

## Tourism Parameters
alpha0 = 9.6448
alpha1 = -.003
alpha2 = .00004

# Biological parameters
r = 0.79787
K = 2.66e5
x0 = K*.5

## Economic Parameters
price = 8.954
cost = 97.79 # $/trip
cost = cost / 63 # $/kg @ x0
discount = 0.05


# Preference parameters
beta = 0 # cost of illegal fishing, $/kg @ x0
T = 50

## Financing parameters
v =rep(price*.1,T) # Tax per unit catch each time step
LF = rep(18.52*34,T) # Licensing fee each time step

## Optimization parameters
bvec = seq(0,2,0.1) # Define vector of B/BMSY for optimal harvest and enforcement functions
tolF = .01 # Tolerance value for fishing effort optimization
tolE = .01 # Tolerance value for enforcement effort optimization

## Functions
costEnforcementFunc = function(enfEffort)
{
  cost = enfVariable * enfEffort
  return(cost)
}

fineProbFunc = function(enfEffort)
{
  prob = approx(enfVec,enfDet,enfEffort)$y
  return(prob)
}

profitPrivateLegalFunc = function (QLegal,x,licFee,tax)
{
  if(x==0) profit = -licFee else profit = price * QLegal - cost * QLegal * (x0 / x) - licFee - tax * QLegal
  return(profit)
}

profitPrivateIUUFunc = function (QIUU,x,enfEffort)
{
  if(x==0) profit = 0 else profit = price * QIUU - (cost + beta) * QIUU * (x0 / x) - fine * fineProbFunc(enfEffort) * QIUU * (x0 / x)
  return(profit)
}

tourismRevenueOptimal = function (xt)
{
  revenue = ((alpha0 + alpha2 * xt) / (-2*alpha1)) * (alpha0 + (alpha0 + alpha2 * xt) / -2 + alpha2 * xt)
  return (revenue)
}

profitSocialFunc = function(QLegal, QIUU, x, enfEffort, licFee, tax)
{
  profit = profitPrivateLegalFunc(QLegal,x,licFee,tax) + 
    tourismRevenueOptimal(x) + 
#    fineProbFunc(enfEffort) * fine * QIUU * (x0 / x)- 
    - costEnforcementFunc(enfEffort) +
    licFee +
    tax * QLegal
  return (profit)
}

stockGrowthFunc <- function (x)
{
  x1 = x + r * x - r * x^2 / K
  if (x<0) x1=0
  return (x1)
}

NPV = function(discount,P)
{
  pv=vector()
  for (i in 1:length(P))
  {
    pv[i] = P[i] / (1+discount)^i
  }
  return(sum(pv))
}

## Optimization to determine vector of optimal fishing effort and optimal enforcement effort based on B/BMSY
## Choose Policy
policy = 3 # 1: FMSY at all time steps; 2: 0 until biomass reaches BMSY, FMSY afterwards; 3: Optimal policy
optimal = RunDynamicOptEnf(K,r,price,cost,discount,bvec,tolF,tolE, LF[1], v[1],pol=policy)

# Run forward projection model
biomass = vector()
QLegal = vector()
QIUU = vector()
profitLegal = vector()
profitIUU = vector()
profitSocial = vector()
enfCost = vector()
tourismRev = vector()
internalFinancing = vector()
enfEffort = vector()
npvCurrent = vector()
pv=vector()

for (i in 1:T)
{
  if (i == 1) biomass[i] = x0 else biomass[i] = ifelse(biomass[i-1] - QLegal[i-1] - QIUU[i-1] >= 0, stockGrowthFunc(biomass[i-1] - QLegal[i-1] - QIUU[i-1]), 0)
  
  enfEffort[i] = approx(bvec,optimal$Enforcement,biomass[i]/(K/2))$y
  QIUU[i] = approx(bvec,optimal$fIUU,(biomass[i])/(K/2))$y * (biomass[i])
  QLegal[i] = approx(bvec,optimal$Policy,biomass[i]/(K/2))$y * (biomass[i] - QIUU[i])
  
    
  tourismRev[i] = tourismRevenueOptimal(biomass[i])
  
  if (i==1) enfCost[i] = enfFixed
  else enfCost[i] = costEnforcementFunc(enfEffort[i])
  
  if (biomass[i]==0) internalFinancing[i] = tourismRev[i] + LF[i] else internalFinancing[i] = tourismRev[i] + LF[i] + v[i] * QLegal [i] + QIUU[i] * fineProbFunc(enfEffort[i]) * fine * (x0 / biomass[i])
  
  npvCurrent[i] = NPV(discount,internalFinancing-enfCost)
  
  pv[i] = profitSocialFunc (QLegal[i], QIUU[i], biomass[i], enfEffort[i], LF[i], v[i])
  
  #  browser()
}
breakEven=min(which(npvCurrent>0))
npv= NPV(discount,pv)

## Plot Results

par(mfrow=c(3,3),oma = c(0, 0, 6, 0))
plot(optimal$Policy~bvec,xlab="B/BMSY",type="l",ylab="Optimal fLegal",col="red")
plot(optimal$Enforcement~bvec,xlab="B/BMSY",type="l",ylab="Optimal enfEffort",col="red")
plot(optimal$fIUU~bvec,xlab="B/BMSY",type="l",ylab="Optimal fIUU",col="red")
curve(fineProbFunc,from=0,to=1,xlab="Enforcement Effort",col="blue")
curve(costEnforcementFunc,from=0,to=1,xlab="Enforcement Effort",col="blue")
plot(enfEffort,xlab="Years",type="l")
plot(biomass/(K/2),xlab="Years",type="l",ylab="B/BMSY")
plot(QLegal/biomass,xlab="Years",type="l",ylab="f_Legal")
plot(QIUU/biomass,xlab="Years",type="l",ylab="f_Illegal")

mtext(paste("Policy",policy,"; NPV =",round(npv,0),"\n B/B0 = ",signif(x0/K,3),"; Beta = ",beta,"\n Price =",signif(price,3),"[$/kg] ; Cost =",signif(cost,3),"[$/kg] [x0/x] ; Fine = ",signif(fine,3),"[$/kg] [x0/x]"), side=3, line=1, outer=TRUE, cex=1, font=2)

par(mfrow=c(1,1),mar = c(5,4,0,4) + 0.1)
plot(enfCost,type="l",xlab="Years",ylab="$",col="red",lwd=2,ylim=range(c(internalFinancing,enfCost)))
par(new=T)
plot(internalFinancing,type="l",col="blue",lwd=2,ylim=range(c(internalFinancing,enfCost)),ylab = "",xlab="")
par(new=T)
plot(tourismRev,type="l",col="green",lwd=2,ylim=range(c(internalFinancing,enfCost)),ylab = "",xlab="")
par(new=T)
plot(v*QLegal,type="l",col="orange",lwd=2,ylim=range(c(internalFinancing,enfCost)),ylab = "",xlab="")

mtext(paste("Policy",policy,"\n Break Even point =",breakEven,"years"), side=3, line=1, cex=1, font=1)
legend("topright",legend=c("Total Internal Financing Revenue","Tourism Revenue","Fishery Tax Revenue","Enforcement Costs"),col=c("blue","green","orange","red"),"l",lwd=2)
abline(v=breakEven)