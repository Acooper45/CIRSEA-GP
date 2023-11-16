## Set working directory
#setwd("~/Gavin's Dropbox/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")
#setwd("~/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")

## Global Functions
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

NPV = function(discount,P)
{
  pv=vector()
  for (i in 1:length(P))
  {
    pv[i] = P[i] / (1+discount)^i
  }
  return(sum(pv))
}

breakevenNPV = function(discount,P)
{
  for (i in 1:length(P))
  {
    npvCurrent = NPV(discount,P[1:i])
    if (npvCurrent >= 0) pointer = i-1
    if (npvCurrent >= 0) break
    if (i==length(P) & npvCurrent <= 0) pointer=NA
  }
  return(pointer)
}

numDives = function(xt)
{
  nd = ((alpha0 + alpha2 * xt) / (-2*alpha1))
  return(nd)
}

OP = function(xt)
{
  op = (alpha0 + (alpha0 + alpha2 * xt) / -2 + alpha2 * xt)
  return(op)
}

tourismRevenueOptimal = function (xt)
{
  revenue = ((alpha0 + alpha2 * xt) / (-2*alpha1)) * (alpha0 + (alpha0 + alpha2 * xt) / -2 + alpha2 * xt)
  return (revenue)
}

stockGrowthFunc <- function (x)
{
  x1 = x + r * x - r * x^2 / K
  if (x<0) x1=0
  return (x1)
}

## Read in optimization functions

wrapperFunction = function(pol,arch,B_B0,Beta,Fine,Tax,Fee,wTax){
  
## Define Model Run Parameters
x0 = K*B_B0
v =rep(price*Tax,T) # Tax per unit catch each time step
w = rep(wTax,T) # Tax per unit tourism revenue
LF = rep(Fee*34,T) # Licensing fee each time step
#beta=Beta
tax=Tax
beta=cost*Beta
fine=Fine
archetype=arch
pol=pol

fineExp = function(enfEffort,fine,Q,Qstar)
{
  if (Q>Qstar) {
    expect = fineProbFunc(enfEffort) * (Q-Qstar) * fine
  }
  else{
    expect = 0
  }
  return(expect)
}

profitPrivateFunc = function (Qstar,Q,x,licFee,tax,enfEffort,fine)
{
#  taxRev = min(tax*Qstar,tax*Q)
#  if (Q>Qstar) taxRev = fineProbFunc(enfEffort) * (Q-Qstar)* tax + Qstar*tax else taxRev = tax*price*Q
  taxRev = min(tax*price*Qstar,tax*price*Q)
  profit = price * Q - cost * Q^2/x - licFee - taxRev - fineExp(enfEffort,fine,Q,Qstar)
  if (x==0) profit = - licFee
  return(profit)
}

profitPrivateOpt = function (Q,x,licFee,tax)
{
  #  taxRev = min(tax*Qstar,tax*Q)
  #  if (Q>Qstar) taxRev = fineProbFunc(enfEffort) * (Q-Qstar)* tax + Qstar*tax else taxRev = tax*price*Q
#  browser()
  if (archetype == 1) profit = price * Q - cost * Q^2/x - licFee - tax*Q
  if (archetype == 2) profit = tourismRevenueOptimal(x)
  if (archetype == 3) profit = price * Q - cost * Q^2/x - licFee - tax*Q + tourismRevenueOptimal(x)
  if (x==0) profit = - licFee
  return(profit)
}

QFunc = function(enfE,x,q_star,tax)
{

  #  if (archetype == 2) q = max((price - fineProbFunc(enfE)*fine)*x/(2*(cost+beta)),q_star) else q = max((price - tax - fineProbFunc(enfE)*fine)*x/(2*(cost+beta)),q_star)

  if (archetype == 2) q = max((price - fineProbFunc(enfE)*fine)*x/(2*(cost+beta)),q_star) else q = max((price - tax - fineProbFunc(enfE)*fine)*x/(2*(cost+beta)),q_star)
  return(q)
}

QFuncZero = function(enfE,x,q_star,tax)
{
  
  #  if (archetype == 2) q = max((price - fineProbFunc(enfE)*fine)*x/(2*(cost+beta)),q_star) else q = max((price - tax - fineProbFunc(enfE)*fine)*x/(2*(cost+beta)),q_star)
  
  if (archetype == 2) q = (price - fineProbFunc(enfE)*fine)*x/(2*(cost+beta)) - q_star else q = (price - tax - fineProbFunc(enfE)*fine)*x/(2*(cost+beta)) - q_star
  return(q)
}

# QFuncOpt = function(x)
# {
#   q = (price - tax)*x/(2*cost)
#   return(q)
# }


if (archetype ==1) {
  profitSocialFunc = function(Q, Qstar, x, enfE, licFee, tax,fine)
  {
    taxRev = min(tax*price*Qstar,tax*price*Q)
    profit = profitPrivateFunc(Qstar,Q,x,licFee,tax,enfE,fine) + 
      fineExp(enfE,fine,Q,Qstar) - 
      costEnforcementFunc(enfE) +
      licFee +
      taxRev
    return (profit)
  }
}

if (archetype ==2) {
  profitSocialFunc = function(Q, Qstar, x, enfE, licFee, tax,fine)
  {
    profit = tourismRevenueOptimal(x) + 
      fineExp(enfE,fine,Q,Qstar) -  
      costEnforcementFunc(enfE)
    return (profit)
  }
}

if (archetype ==3) {
  profitSocialFunc = function(Q, Qstar, x, enfE, licFee, tax,fine)
  {
    taxRev = min(tax*price*Qstar,tax*price*Q)
    profit = profitPrivateFunc(Qstar,Q,x,licFee,tax,enfE,fine) + 
      tourismRevenueOptimal(x) + 
      fineExp(enfE,fine,Q,Qstar) -  
      costEnforcementFunc(enfE) +
      licFee +
      taxRev
    return (profit)
  }
}

## Optimization to determine vector of optimal fishing effort and optimal enforcement effort based on B/BMSY

profitOptimEnf = function(enfEffort,b,p,K,c,r,V,bvec,delta, f, licFee, tax)
{  
  
  bcurrent = b*K/2
  
  if (b==0) qstar = 0 else qstar = f*bcurrent
  
  if (b==0) q = 0 else q = QFunc(enfEffort,bcurrent,qstar,tax)
    
  profit = profitSocialFunc (q,qstar, bcurrent, enfEffort,licFee, tax, fine)
  
  bnext=max(min(bvec),stockGrowthFunc(bcurrent - q)/(K/2))
  bnext=min(max(bvec),bnext)
  
  out= approx(bvec,t(V),bnext)
  Vnext= out$y
  
  negout= -(profit + delta*Vnext)
  return(negout)
}

profitOptimQ = function(f,b,p,K,c,r,V,bvec,delta,licFee, tax)
{  
  
  bcurrent = b*K/2
  
  if (b==0) q = 0 else q = f*bcurrent
  
  profit = profitPrivateOpt(q,bcurrent,licFee,tax)
  
  bnext=max(min(bvec),stockGrowthFunc(bcurrent - q)/(K/2))
  bnext=min(max(bvec),bnext)
  
  out= approx(bvec,t(V),bnext)
  Vnext= out$y
  
  negout= -(profit + delta*Vnext)
  return(negout)
}

profitOptimTwo = function(params=c(f,enfEffort),b,p,K,c,r,V,bvec,delta,licFee, tax)
{  
  
  bcurrent = b*K/2
  
  if (b==0) qstar = 0 else qstar = params[1]*bcurrent
  
  if (b==0) q = 0 else q = QFunc(params[2],bcurrent,qstar,tax)
    
  profit = profitSocialFunc (q,qstar, bcurrent, params[2],licFee, tax, fine)
  
  bnext=max(min(bvec),stockGrowthFunc(bcurrent - q)/(K/2))
  bnext=min(max(bvec),bnext)
  
  out= approx(bvec,t(V),bnext)
  Vnext= out$y
  negout= -(profit + delta*Vnext)
  return(negout)
}

RunDynamicOptEnf= function(K,r,p,c,disc,bvec,tolF,tolE,licFee, tax,pol)
{
  
  # K<- K
  # r<- r
  # p<- price
  # c<- cost
  # disc<- discount
  # bvec<- bvec
  # tol<- tol
  #   
  
  delta= 1/(1+disc) #Discount parameter
  t=0
  
  if (pol==0) f1 = rep(0,length(bvec))
  if (pol==1) f1 = rep(r/2,length(bvec))
  if (pol==2) f1 = c(rep(0,floor(length(bvec)/2)),rep(r/2,ceiling(length(bvec)/2)))
  if (pol==3) f1 = matrix(0,length(bvec),1)
  if (pol==4) f1 = matrix(0,length(bvec),1)
  
  Vnew= matrix(0,length(bvec),1)
  e1 = matrix(1,length(bvec),1)
  
  diffF= 10*tolF
  diffE= 10*tolE
  
  while (t<4 | diffF>tolF | diffE>tolE)
  {
    t= t+1
    V= Vnew
    oldf1= f1
    olde1 = e1
    for (i in 1:length(bvec))
    {
      b= bvec[i]
      if(i==1 | i==2)
      {guess= f1[1]
       guessEnf=e1[1]}
      else
      {guess= f1[i-1]
       guessEnf=e1[i-1]}
      
      if (pol==0 | pol==1 | pol==2)
      {
        FishOut= optim(par=guessEnf,fn=profitOptimEnf,lower=0.0001,upper=.9999,b=b,p=p,K=K,c=c,r=r,V=V,bvec=bvec,delta=delta,f=f1[i],licFee=licFee, tax=tax,method="L-BFGS-B")
        e1[i]= FishOut$par
        Vnew[i]= -FishOut$value
      }
      if (pol==3)
      {
        FishOut= optim(par=c(guess,guessEnf),fn=profitOptimTwo,lower=0.0001,upper=.9999,b=b,p=p,K=K,c=c,r=r,V=V,bvec=bvec,delta=delta,licFee=licFee, tax=tax,method="L-BFGS-B")
        f1[i]= FishOut$par[1]
        e1[i]= FishOut$par[2]
        Vnew[i]= -FishOut$value
      }
      if (pol==4)
      {
        FishOut= optim(par=guess,fn=profitOptimQ,lower=0.0001,upper=.9999,b=b,p=p,K=K,c=c,r=r,V=V,bvec=bvec,delta=delta,licFee=licFee, tax=tax,method="L-BFGS-B")
        f1[i]= FishOut$par
        e1[i]= 0
        Vnew[i]= -FishOut$value
      }
    } #Close bvec loop
    
    diffF= sum(abs(f1-oldf1))
    diffE= sum(abs(e1-olde1))
    
  }# Close while loop

if (pol==4) {
  
  delta= 1/(1+disc) #Discount parameter
  t=0
  if (pol==4) f1 = f1
  
  Vnew= matrix(0,length(bvec),1)
  e1 = matrix(1,length(bvec),1)
  
  diffF= 10*tolF
  diffE= 10*tolE
while (t<4 | diffF>tolF | diffE>tolE)
{
  t= t+1
  V= Vnew
  oldf1= f1
  olde1 = e1
  for (i in 1:length(bvec))
  {
    b= bvec[i]
    if(i==1 | i==2)
    {guess= f1[1]
     guessEnf=e1[1]}
    else
    {guess= f1[i-1]
     guessEnf=e1[i-1]}
    
      FishOut= optim(par=guessEnf,fn=profitOptimEnf,lower=0.0001,upper=.9999,b=b,p=p,K=K,c=c,r=r,V=V,bvec=bvec,delta=delta,f=f1[i],licFee=licFee, tax=tax,method="L-BFGS-B")
      e1[i]= FishOut$par
      Vnew[i]= -FishOut$value
  } #Close bvec loop
  
  diffF= sum(abs(f1-oldf1))
  diffE= sum(abs(e1-olde1))
  
}# Close while loop
}  

  eZero = vector()
  for (i in 1:length(bvec))
  {
    eZero[i] = uniroot(function(enfEffort) QFuncZero(enfEffort,x=bvec[i]*K/2,q_star=f1[i]*bvec[i]*K/2,tax=tax),c(0.0001,0.9999))$root
  }

  fOpt = vector()
  for (i in 1:length(bvec))
  {
    if (bvec[i]==0) fOpt[i] = 0 else fOpt[i] = QFunc(e1[i],bvec[i]*K/2,f1[i]*bvec[i]*K/2,tax)/(bvec[i]*K/2)
  }

  fOptZero = vector()
  for (i in 1:length(bvec))
  {
    if (bvec[i]==0) fOptZero[i] = 0 else fOptZero[i] = QFunc(eZero[i],bvec[i]*K/2,f1[i]*bvec[i]*K/2,tax)/(bvec[i]*K/2)
  }

  return(list(Policy=f1,Enforcement=e1,fOpt=fOpt,eZero=eZero,fOptZero=fOptZero))
} #Close function

optimal = RunDynamicOptEnf(K,r,price,cost,discount,bvec,tolF,tolE, LF[1], v[1],pol=pol)

# Run forward projection model
biomass = vector()
Qstar = vector()
Q = vector()
profitFishing = vector()
profitTourism = vector()
profitIndustry = vector()
enfCost = vector()
tourismRev = vector()
enfRev=vector()
fishingRev=vector()
licensingRev=vector()
internalFinancing = vector()
enfEffort = vector()
npvCurrent = vector()
pv=vector()
socialRev=vector()



for (i in 1:T)
{
  if (i == 1) biomass[i] = x0 else biomass[i] = ifelse(biomass[i-1] - Q[i-1] >= 0, stockGrowthFunc(biomass[i-1] -Q[i-1]), 0)
  
  enfEffort[i] = approx(bvec,optimal$Enforcement,biomass[i]/(K/2))$y
  Qstar[i] = approx(bvec,optimal$Policy,biomass[i]/(K/2))$y * biomass[i]
  Q[i] = approx(bvec,optimal$fOpt,biomass[i]/(K/2))$y * biomass[i]

  if (archetype == 2) profitFishing[i] =0 else profitFishing[i]=profitPrivateFunc(Qstar[i],Q[i],biomass[i],LF[i],v[i]/price,enfEffort[i],fine)
  if (archetype == 1) profitTourism[i] =0  else profitTourism[i]=(1-w[i])*tourismRevenueOptimal(biomass[i] - Q[i])
  profitIndustry[i]=profitFishing[i] + profitTourism[i]
    
  if (archetype == 1) tourismRev[i]=0 else tourismRev[i] = w[i]*tourismRevenueOptimal(biomass[i] - Q[i])
  fishingRev[i] = min(Qstar[i]*v[i],Q[i]*v[i])
  
  
  
#  if (archetype == 2) fishingRev[i]=0 else fishingRev[i] = min(v[i] * Qstar[i], v[i] * Q[i])
  if (archetype == 2) licensingRev[i]=0 else licensingRev[i] = LF[i]
  enfRev[i] = fineExp(enfEffort[i],fine,Q[i],Qstar[i])
  
  socialRev[i] = fishingRev[i] + tourismRev[i] + enfRev[i] + licensingRev[i]

  if (i==1) enfCost[i] = enfFixed + costEnforcementFunc(enfEffort[i])
  else enfCost[i] = costEnforcementFunc(enfEffort[i])
  
  if (biomass[i]==0) internalFinancing[i] = tourismRev[i] + LF[i] else internalFinancing[i] = tourismRev[i] + fishingRev[i] + enfRev[i] + licensingRev[i]
  
  npvCurrent[i] = NPV(discount,internalFinancing-enfCost)
  
}
npv= NPV(discount,internalFinancing-enfCost)
return(list(biomass=biomass,
            Qstar=Qstar,
            Q=Q,
            fOpt=optimal$fOpt,
            profitFishing=profitFishing,
            profitTourism=profitTourism,
            profitIndustry=profitIndustry,
            enfEffort=enfEffort,
            enfCost=enfCost,
            tourismRev=tourismRev,
            fishingRev=fishingRev,
            licensingRev=licensingRev,
            enfRev=enfRev,
            internalFinancing=internalFinancing,
            enfEffort=enfEffort,
            npvCurrent=npvCurrent,
            pv=pv,
            breakEven=breakEven,
            npv=npv,
            eZero=optimal$eZero,
            fOptZero=optimal$fOptZero,
            Policy=optimal$Policy,
            Enforcement=optimal$Enforcement,
            cost=cost,
            socialRev=socialRev,
            breakEven=breakEven))
}
