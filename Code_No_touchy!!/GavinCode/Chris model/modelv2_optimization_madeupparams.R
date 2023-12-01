## Set working directory
setwd("~/Gavin's Dropbox/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")
#setwd("~/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")

## Read in optimization functions
source("profitOptim.R")
source("RunDynamicOptimizationEnf.R")

## Enforcement Parameters
ECFitPar1 = 2.31*500
ECFitPar2 = 0
ECFitPar3 = 0.29

EPFitPar1 = 0.89
EPFitPar2 = 0.19
EPFitPar3 = 0.045

fine = 10

## Economic Parameters
price = 10
cost = 5
discount = 0.05

## Tourism Parameters
alpha0 = 9.6448
alpha1 = -.003
alpha2 = .00004
fixedTourismPrice = 10

# Biological parameters
r = 0.2
K = 1000000
x0 = 150000


# Preference parameters
beta = 1
T = 50

## Financing parameters
v = rep(1,T) # Tax per unit catch each time step
LF = rep(50,T) # Licensing fee each time step

## Optimization parameters
bvec = seq(0,2,0.1) # Define vector of B/BMSY for optimal harvest and enforcement functions
tolF = .01 # Tolerance value for fishing effort optimization
tolE = .01 # Tolerance value for enforcement effort optimization

## Functions
costEnforcementFunc = function(enfEffort)
{
  cost = ECFitPar1 * enfEffort ^ ECFitPar3 + ECFitPar2
  return(cost)
}

fineProbFunc = function(enfEffort)
{
  prob = EPFitPar1 / (1+exp ((EPFitPar2-enfEffort)/EPFitPar3))
  return(prob)
}

profitPrivateLegalFunc = function (QLegal,x,licFee,tax)
{
  if(x==0) profit = -licFee else profit = price * QLegal - cost * QLegal ^2 / x - licFee - tax * QLegal
  return(profit)
}

profitPrivateIUUFunc = function (QIUU,x,enfEffort)
{
  profit = beta * (price * QIUU - cost * QIUU^2 / x) - fine * fineProbFunc(enfEffort) * QIUU
  return(profit)
}

tourismRevenueFixed = function (xt)
{
  revenue = fixedTourismPrice * (alpha0 + alpha2 * xt) / (-2*alpha1)
  return (revenue)
}

tourismRevenueOptimal = function (xt)
{
  revenue = ((alpha0 + alpha2 * xt) / (-2*alpha1)) * (alpha0 + (alpha0 + alpha2 * xt) / -2 + alpha2 * xt)
  return (revenue)
}

profitSocialFunc = function(QLegal, QIUU, x, enfEffort, licFee, tax)
{
  profit = profitPrivateLegalFunc(QLegal,x,licFee,tax) + 
#    tourismRevenueFixed(x) - 
#    fineProbFunc(enfEffort) * fine * QIUU - 
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

# ## Assume static enforcement over time
# enfEffortIntercept = .5
# enfEffortSlope = -.001
# 
# ## Assume static QLegal over time
# QLegalIntercept = 2500
# QLegalSlope = 250
# 
# QLegal=vector()
# enfEffort = vector()
# for (i in 1:T)
# {
#   QLegal[i] = QLegalIntercept + QLegalSlope * (i-1)
#   enfEffort[i] = enfEffortIntercept + enfEffortSlope * (i-1)
# }

## Run optimization to determine vector of optimal fishing effort and optimal enforcement effort based on B/BMSY
## Choose Policy
policy = 3 # 1: FMSY at all time steps; 2: 0 until biomass reaches BMSY, FMSY afterwards; 3: Optimal policy
optimal = RunDynamicOptEnf(K,r,price,cost,discount,bvec,tolF,tolE, LF[1], v[1],pol=policy)

## Plot optimal paths

par(mfrow=c(2,2),oma = c(0, 0, 4, 0))

plot(optimal$Policy~bvec,xlab="B/BMSY",type="l",ylab="Optimal fLegal")
plot(optimal$Enforcement~bvec,xlab="B/BMSY",type="l",ylab="Optimal enfEffort")
plot(optimal$fIUU~bvec,xlab="B/BMSY",type="l",ylab="Optimal fIUU")
mtext(paste("Policy",policy), side=3, line=1, outer=TRUE, cex=2, font=2)

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

for (i in 1:T)
{
  if (i == 1) biomass[i] = x0 else biomass[i] = ifelse(biomass[i-1] - QLegal[i-1] - QIUU[i-1] >= 0, stockGrowthFunc(biomass[i-1] - QLegal[i-1] - QIUU[i-1]), 0)
  
  enfEffort[i] = approx(bvec,optimal$Enforcement,biomass[i]/(K/2))$y
  QLegal[i] = approx(bvec,optimal$Policy,biomass[i]/(K/2))$y * biomass[i]
  QIUU[i] = approx(bvec,optimal$fIUU,biomass[i]/(K/2))$y * (biomass[i] - QLegal[i])
    
  profitLegal[i]=profitPrivateLegalFunc(QLegal[i],biomass[i],LF[i],v[i])
  profitIUU[i]=profitPrivateIUUFunc(QIUU[i],biomass[i],enfEffort[i])
  profitSocial[i]=profitSocialFunc(QLegal[i],QIUU[i],biomass[i],enfEffort[i],LF[i],v[i])
  enfCost[i] = costEnforcementFunc(enfEffort[i])
  tourismRev[i] = tourismRevenueFixed(biomass[i])
  internalFinancing[i] = tourismRev[i] + LF[i] + v[i] * QLegal [i] + QIUU[i] * fineProbFunc(enfEffort[i]) * fine
  
}

## Plot Results

par(mfrow=c(3,4),oma = c(0, 0, 4, 0))

curve(fineProbFunc,from=0,to=1,xlab="Enforcement Effort")
curve(costEnforcementFunc,from=0,to=1,xlab="Enforcement Effort")
plot(enfEffort,xlab="Years",type="l")
plot(QLegal,xlab="Years",type="l")
plot(QIUU,xlab="Years",type="l")
plot(profitLegal,xlab="Years",type="l")
plot(profitIUU,xlab="Years",type="l")
plot(profitSocial,xlab="Years",type="l")
plot(biomass,xlab="Years",type="l")
plot(tourismRev,xlab="Years",type="l")
plot(enfCost,xlab="Years",type="l")
NPV(discount,profitSocial)

par(mfrow=c(1,1),mar = c(5,4,4,4) + 0.1)
plot(enfCost,type="l",xlab="Years",ylab="Enforcement Costs",col="red",lwd=2,ylim=range(c(internalFinancing,enfCost)))
par(new=T)
plot(internalFinancing,type="l",col="blue",lwd=2,ylim=range(c(internalFinancing,enfCost)),ylab = "",xlab="")
axis(side = 4)
mtext("Internal Financing Revenue", side = 4, line = 3, cex = par("cex.lab"))
legend("center",legend=c("Enforcement Costs","Internal Financing Revenue"),col=c("red","blue"),"l",lwd=2)
