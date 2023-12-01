#setwd("~/Gavin's Dropbox/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")
setwd("~/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")
## Enforcement Parameters
ECFitPar1 = 2.31 * 30000
ECFitPar2 = 0
ECFitPar3 = 0.29

EPFitPar1 = 0.89
EPFitPar2 = 0.19
EPFitPar3 = 0.045

fine = 10

## Economic Parameters
price = 10
cost = 5
r = 0.05

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
  profit = price * QLegal - cost * QLegal ^2 / x - licFee - tax * QLegal
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
    tourismRevenueFixed(x) + 
    fineProbFunc(enfEffort) * fine * QIUU - 
    costEnforcementFunc(enfEffort) +
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

NPV = function(r,P)
{
  pv=vector()
  for (i in 1:length(P))
  {
    pv[i] = P[i] / (1+r)^i
  }
  return(sum(pv))
}

## Assume static enforcement over time
enfEffortIntercept = .5
enfEffortSlope = -.001

## Assume static QLegal over time
QLegalIntercept = 2500
QLegalSlope = 250

QLegal=vector()
enfEffort = vector()
for (i in 1:T)
{
  QLegal[i] = QLegalIntercept + QLegalSlope * (i-1)
  enfEffort[i] = enfEffortIntercept + enfEffortSlope * (i-1)
}

# Run model
biomass = vector()
QIUU = vector()
profitLegal = vector()
profitIUU = vector()
profitSocial = vector()
enfCost = vector()
tourismRev = vector()
internalFinancing = vector()

bioEconModel = function(QLeg,Enf)
{
  for (i in 1:T)
  {
    if (i == 1) biomass[i] = x0 else biomass[i]=stockGrowthFunc(biomass[i-1])
    if (biomass[i] == 0) break
    QIUU[i]=optimize(profitPrivateIUUFunc,x=biomass[i],enfEffort=Enf[i],c(0,biomass[i]),maximum=TRUE)$objective
    biomass[i] = ifelse(biomass[i] - QLeg[i] - QIUU[i] >= 0, biomass[i] - QLeg[i] - QIUU[i], 0)
  }
  return(NPV(r,profitSocial))
}

#QOpt=optimx(par=QLegal,fn=bioEconModel,Enf=enfEffort,upper=K,lower=0,method = c("L-BFGS-B"), control = list( maximize = TRUE ))
#QLegal = coef(QOpt)

for (i in 1:T)
{
  if (i == 1) biomass[i] = x0 else biomass[i]=stockGrowthFunc(biomass[i-1])
  if (biomass[i] == 0) break
  QIUU[i]=optimize(profitPrivateIUUFunc,x=biomass[i],enfEffort=enfEffort[i],c(0,biomass[i]),maximum=TRUE)$objective
  profitLegal[i]=profitPrivateLegalFunc(QLegal[i],biomass[i],LF[i],v[i])
  profitIUU[i]=profitPrivateIUUFunc(QIUU[i],biomass[i],enfEffort[i])
  profitSocial[i]=profitSocialFunc(QLegal[i],QIUU[i],biomass[i],enfEffort[i],LF[i],v[i])
  enfCost[i] = costEnforcementFunc(enfEffort[i])
  tourismRev[i] = tourismRevenueFixed(biomass[i])
  internalFinancing[i] = tourismRev[i] + LF[i] + v[i] * QLegal [i] + QIUU[i] * fineProbFunc(enfEffort[i]) * fine
  biomass[i] = ifelse(biomass[i] - QLegal[i] - QIUU[i] >= 0, biomass[i] - QLegal[i] - QIUU[i], 0)
}

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
#plot(enfCost,xlab="Years",type="l")
#NPV(r,profitSocial)

plot(profitSocial,xlab="Years",type="l",col="blue",lwd=2,ylim=range(c(profitSocial,enfCost)))
par(new=T)
plot(enfCost,type="l",col="red",lwd=2,ylim=range(c(profitSocial,enfCost)),axes=FALSE )
par(new=T)
plot(internalFinancing,type="l",col="green",lwd=2,ylim=range(c(profitSocial,enfCost)),axes=FALSE)