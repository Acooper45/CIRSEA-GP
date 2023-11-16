## Set working directory
setwd("~/Gavin's Dropbox/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")
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
    if (npvCurrent > 0) pointer = i-1
    if (npvCurrent > 0) break
  }
  return(pointer)
}
## Read in optimization functions

wrapperFunction = function(pol,arch,B_B0,Beta,Fine,Tax,Fee,wTax){
  
## Define Model Run Parameters
x0 = K*B_B0
v =rep(price*Tax,T) # Tax per unit catch each time step
w = rep(wTax,T) # Tax per unit tourism revenue
LF = rep(Fee*34,T) # Licensing fee each time step
#beta=Beta
beta=cost*Beta
fine=Fine
archetype=arch
pol=pol

 
profitPrivateLegalFunc = function (QLegal,x,licFee,tax)
{
  if(x==0) profit = -licFee else profit = price * QLegal - cost * (QLegal/x)^2 - licFee - tax * QLegal
  return(profit)
}

#  profitPrivateIUUFunc = function (QIUU,x=387275,enfEffort=1)
#  {
#    if(x==0) profit = 0 else profit = price * QIUU - (cost + beta) * (QIUU)^2 - fineExp(enfEffort,fine,QIUU,x)
#    return(profit)
#  }
#  curve(profitPrivateIUUFunc,from=0,to=K)


# fineExp = function(enfEffort,fine,QIUU,x)
# {
#   if (x==0) expect=0 else expect = fineProbFunc(enfEffort) * fine * QIUU
#   return(expect)
# }

# profitPrivateIUUFunc = function (QIUU,x,enfEffort)
# {
#   if(x==0) profit = 0 else profit = price * QIUU - (cost + beta) * (QIUU/x)^2 - fineExp(enfEffort,fine,QIUU,x)
#   return(profit)
# }
# 
# QIUUFunc = function(enfEffort,x)
# {
#   if (x==0) qIUU = 0 else qIUU = ((price - fineProbFunc(enfEffort)*fine)*x^2)/(2*(cost+beta))
#   if (qIUU<0) qIUU=0
#   if(qIUU>x) qIUU=x
#   return(qIUU)
# }

fineExp = function(enfEffort,fine,QIUU,x)
{
  if (x==0) expect=0 else expect = fineProbFunc(enfEffort) * fine * (QIUU/x)^2
  return(expect)
}

QIUUFunc = function(enfEffort,x)
{
  if (x==0) qIUU = 0 else qIUU = (price*x^2)/(2*(cost+beta+fineProbFunc(enfEffort)*fine))
  if (qIUU<0) qIUU=0
  if(qIUU>x) qIUU=x
  return(qIUU)
}

tourismRevenueOptimal = function (xt)
{
  revenue = ((alpha0 + alpha2 * xt) / (-2*alpha1)) * (alpha0 + (alpha0 + alpha2 * xt) / -2 + alpha2 * xt)
  return (revenue)
}

if (archetype ==1) {
profitSocialFunc = function(QLegal, QIUU, x, enfEffort, licFee, tax)
  {
    profit = profitPrivateLegalFunc(QLegal,x,licFee,tax) - 
#    fineExp(enfEffort,fine,QIUU,x) - 
    costEnforcementFunc(enfEffort) +
    licFee +
    tax * QLegal
    return (profit)
  }
}

if (archetype ==2) {
  profitSocialFunc = function(QLegal, QIUU, x, enfEffort, licFee, tax)
  {
      profit = tourismRevenueOptimal(x) - 
#      fineExp(enfEffort,fine,QIUU,x) -  
      costEnforcementFunc(enfEffort)
      return (profit)
  }
}

if (archetype ==3) {
  profitSocialFunc = function(QLegal, QIUU, x, enfEffort, licFee, tax)
  {
      profit = profitPrivateLegalFunc(QLegal,x,licFee,tax) + 
      tourismRevenueOptimal(x) - 
#      fineExp(enfEffort,fine,QIUU,x) -  
      costEnforcementFunc(enfEffort) +
      licFee +
      tax * QLegal
      return (profit)
  }
}

stockGrowthFunc <- function (x)
{
  x1 = x + r * x - r * x^2 / K
  if (x<0) x1=0
  return (x1)
}

## Optimization to determine vector of optimal fishing effort and optimal enforcement effort based on B/BMSY

profitOptimEnf = function(enfEffort,b,p,K,c,r,V,bvec,delta, f, licFee, tax)
{  
  
  bcurrent = b*K/2
  
  QLegal = f*(bcurrent)
  
  if (b==0) QIUU =0 else QIUU = QIUUFunc(enfEffort,bcurrent-QLegal)
  
  profit = profitSocialFunc (QLegal, QIUU, bcurrent, enfEffort, licFee, tax)
  
  bnext=max(min(bvec),stockGrowthFunc(bcurrent - QIUU - QLegal)/(K/2))
  bnext=min(max(bvec),bnext)
  
  out= approx(bvec,t(V),bnext)
  Vnext= out$y
  
  negout= -(profit + delta*Vnext)
  return(negout)
}

profitOptimTwo = function(params=c(f,enfEffort),b,p,K,c,r,V,bvec,delta,licFee, tax)
{  
  
  bcurrent = b*K/2
  
  QLegal = params[1]*(bcurrent)
  
  if (b==0) QIUU =0 else QIUU = QIUUFunc(params[2],bcurrent-QLegal)
  
  profit = profitSocialFunc (QLegal, QIUU, bcurrent, params[2], licFee, tax)
  
  bnext=max(min(bvec),stockGrowthFunc(bcurrent - QLegal - QIUU)/(K/2))
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
  
  f1IUU= matrix(0,length(bvec),1)
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
        FishOut= optim(par=guessEnf,fn=profitOptimEnf,lower=0.0001,upper=0.9999,b=b,p=p,K=K,c=c,r=r,V=V,bvec=bvec,delta=delta,f=f1[i],licFee=licFee, tax=tax,method="L-BFGS-B")
        e1[i]= FishOut$par
        Vnew[i]= -FishOut$value
        #          if (b==0) f1IUU[i]=0 else f1IUU[i] = max(0,optim(par=0.5,fn=profitPrivateIUUFunc,lower=0,upper=b*K/2-f1[i]*b*K/2,x=b*K/2-f1[i]*b*K/2,enfEffort=e1[i],method="L-BFGS-B",control=list(fnscale=-1))$par/(b*K/2-f1[i]*b*K/2))
        #        if (b==0) f1IUU[i]=0 else f1IUU[i] = min(1,max(0,QIUUFunc(e1[i],b*K/2-f1[i]*b*K/2)))/(b*K/2-f1[i]*b*K/2)
      }
      if (pol==3)
      {
        FishOut= optim(par=c(guess,guessEnf),fn=profitOptimTwo,lower=0.0001,upper=0.9999,b=b,p=p,K=K,c=c,r=r,V=V,bvec=bvec,delta=delta,licFee=licFee, tax=tax,method="L-BFGS-B")
        f1[i]= FishOut$par[1]
        e1[i]= FishOut$par[2]
        Vnew[i]= -FishOut$value
        #          if (b==0) f1IUU[i]=0 else f1IUU[i] =  max(0,optim(par=0.5,fn=profitPrivateIUUFunc,lower=0,upper=b*K/2-f1[i]*b*K/2,x=b*K/2-f1[i]*b*K/2,enfEffort=e1[i],method="L-BFGS-B",control=list(fnscale=-1))$par/(b*K/2-f1[i]*b*K/2))
        #        if (b==0) f1IUU[i]=0 else f1IUU[i] = min(1,max(0,QIUUFunc(e1[i],b*K/2-f1[i]*b*K/2)))/(b*K/2-f1[i]*b*K/2)
      }
    } #Close bvec loop
    
    diffF= sum(abs(f1-oldf1))
    diffE= sum(abs(e1-olde1))
    
  }# Close while loop
  
  
  #  return(list(Policy=f1,Enforcement=e1,Value=V,fIUU=f1IUU))
  return(list(Policy=f1,Enforcement=e1))
} #Close function

optimal = RunDynamicOptEnf(K,r,price,cost,discount,bvec,tolF,tolE, LF[1], v[1],pol=pol)


fIUU = vector()
for (i in 1:length(bvec))
{
  if (bvec[i]==0) fIUU[i] = 0 else fIUU[i] = QIUUFunc(optimal$Enforcement[i],bvec[i]*K/2)/(bvec[i]*K/2)
}


# Run forward projection model
biomass = vector()
QLegal = vector()
QIUU = vector()
profitLegal = vector()
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



for (i in 1:T)
{
  if (i == 1) biomass[i] = x0 else biomass[i] = ifelse(biomass[i-1] - QLegal[i-1] - QIUU[i-1] >= 0, stockGrowthFunc(biomass[i-1] - QLegal[i-1] - QIUU[i-1]), 0)
  
  enfEffort[i] = approx(bvec,optimal$Enforcement,biomass[i]/(K/2))$y
  QLegal[i] = approx(bvec,optimal$Policy,biomass[i]/(K/2))$y * biomass[i]
  QIUU[i] = approx(bvec,fIUU,(biomass[i] - QLegal[i])/(K/2))$y * (biomass[i] - QLegal[i])

  if (archetype == 2) profitLegal[i]=0 else profitLegal[i]=profitPrivateLegalFunc(QLegal[i],biomass[i],LF[i],v[i])
  if (archetype == 1) profitTourism[i] =0  else profitTourism[i]=tourismRevenueOptimal(biomass[i] - QLegal[i] - QIUU[i])
  profitIndustry[i]=profitLegal[i] + profitTourism[i]
    
  if (archetype == 1) tourismRev[i]=0 else tourismRev[i] = w[i]*tourismRevenueOptimal(biomass[i] - QLegal[i] - QIUU[i])
  if (archetype == 2) fishingRev[i]=0 else fishingRev[i] = v[i] * QLegal [i]
  if (archetype == 2) licensingRev[i]=0 else licensingRev[i] = LF[i]
  enfRev[i] = fineExp(enfEffort[i],fine,QIUU[i],(biomass[i]-QLegal[i]))
  
  if (i==1) enfCost[i] = enfFixed + costEnforcementFunc(enfEffort[i])
  else enfCost[i] = costEnforcementFunc(enfEffort[i])
  
  if (biomass[i]==0) internalFinancing[i] = tourismRev[i] + LF[i] else internalFinancing[i] = tourismRev[i] + fishingRev[i] + enfRev[i] + licensingRev[i]
  
  npvCurrent[i] = NPV(discount,internalFinancing-enfCost)
  
}
breakEven=min(which(npvCurrent>0))
npv= NPV(discount,internalFinancing-enfCost)
return(list(biomass=biomass,
            QLegal=QLegal,
            profitLegal=profitLegal,
            profitTourism=profitTourism,
            profitIndustry=profitIndustry,
            enfEffort=enfEffort,
            QIUU=QIUU,
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
            Policy=optimal$Policy,
            Enforcement=optimal$Enforcement,
            fIUU=fIUU,
            cost=costTotal0/Q0,
            breakEven=breakEven))
}
