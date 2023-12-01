## Model for "Fisheries Financing for Enforcement" paper
## Gavin McDonald, Lennon Thomas, Chris Costello, Katie Peterson

setwd("~/Gavin's Dropbox/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model")

######### Parameterize Model ##########
##                                            
## Enforcement Parameters
eta = 0.5   # Detection Probability
a = 100     # Enforcement Cost Parameter
f = 3.5       # Size of enforcement fine

## Biological Parameters
r = 0.2   # Logistic Growth Parameter r
K = 1000  # Logistic Growth Parameter K
B0 = 1000 # Starting biomass

MSY = K*r/4 # MSY

## Economic Paramters
p = 1       # Ex-vessel Price
c = 0.7     # Variable Fishing Cost
FK = 1    # Fixed Fishing Cost
r = 0.05    # Social Planner Discount Rate

## Fishery Parameters


## Temporal Parameters
yearsCS = 100 # Time horizon of catch share

##
########################################

######### Simulation Functions #########

## stockEvolve
## Evolve stock to find new stock size (x1) based on old stock size (x0) 
## and insider and outsider catch (q_insider and q_outsider). 
## Use alpha and beta growth parameters.
stockEvolve = function (x0,q_insider,q_outsider)
{
  x1 = x0 + r * x0 - r * x0^2 / K - q_insider - q_outsider
  if (x1<0) x1=0
  return (x1)
}

## penalityProbability1
## Find the probability of penalty (detection and prosecution)
## given enforcement effort (e) and detection probability (eta)
penaltyProbability1 = function (e,eta)
{
  #PI = e / (eta+e)
  PI = e
  return(PI)
}

## penaltyFunction1
## Find the penalty that would be paid
## given actual catch (q) and allowed catch (q_star)
penaltyFunction1 = function(f,q,q_star)
{
  if (q>q_star) F = f * (q - q_star) else F=0
  return(F)
}

## enforcementCost1
## Find the cost of enforcement 
## given eneforcement effort (e) and enforcement cost parameter(a)
enforcementCost1 = function (e,a)
{
  C = e* a
  return(C)
}

## profitInsider
## Find insider fisher profit 
## given old stock size (x0), enforcement effort, actual catch (q), and allowed catch (q_star). 
profitInsider = function (q,x0,q_star,e)
{
  B = p * q - c * q^2 / x0 - FK - penaltyProbability1(e,eta) * penaltyFunction1(f,q,q_star)
  return(B)
}

## profitOutsider
## Find outsider fisher profit 
## given old stock size (x0), enforcement effort, actual catch (q), and allowed catch (q_star). 
profitOutsider = function (x0,q,e)
{
  B = p * q - c * q^2 / x0 - FK - penaltyProbability1(e,eta) * penaltyFunction1(f,q,0)
  return(B)
}

## NPV
## Find NPV given discount rate (r) and vector of profits (P) over time period T
NPV = function(r,P)
{
  pv=vector()
  for (i in 1:length(P))
  {
    pv[i] = P[i] / (1+r)^i
  }
  return(sum(pv))
}

OptHarvest = function(h)
{
  P={}
  x0 = B0
  for (i in 1:yearsCS)
  {
    if (i==1){
      P[i] = profitInsider(B0,B0*h,B0*h,0)
      x0[i] = B0
    }
    else {
      P[i] = profitInsider(x0[i-1]*h,x0[i-1],x0[i-1]*h,0)
      x0[i] = stockEvolve(x0[i-1],h*x0[i-1],0)
  }
  }
  profit = NPV(r,P)
  return(profit)
}

findOpt = optimize(OptHarvest,c(0,1),maximum=TRUE)
u_opt=findOpt$maximum

enf=1
#quota = MSY * 0.75
quota=MSY*.25
P=vector()
x0=vector()
q=vector()
for (i in 1:yearsCS)
{
  if (i==1){
    x0[i]=B0
    harvestOpt = optimize(profitInsider,q_star=quota,e=enf,x0=x0[i],c(0,K),maximum=TRUE)
    if (x0[i]<=0) q[i]=0 else q[i]=harvestOpt$maximum
    #if (x0[i]<=0) q[i]=0 else q[i]=uniroot(profitInsider,x0=x0[i],q_star=quota,e=enf,upper=K,lower=0)$root
    if (x0[i]<=0) P[i]=0 else P[i]=profitInsider(q[i],x0[i],quota,enf)}
  else {
    x0[i] = stockEvolve(x0[i-1],q[i-1],0)
    harvestOpt = optimize(profitInsider,q_star=quota,e=enf,x0=x0[i],c(0,K),maximum=TRUE)
    if (x0[i]<=0) q[i]=0 else q[i]=harvestOpt$maximum
    #if (x0[i]<=0) q[i]=0 else q[i]=uniroot(profitInsider,x0=x0[i],q_star=quota,e=enf,upper=K,lower=0)$root
    if (x0[i]<=0) P[i]=0 else P[i]=profitInsider(q[i],x0[i],quota,enf)}
}
plot(P)
plot(x0)
plot(q)
NPV(r,P)