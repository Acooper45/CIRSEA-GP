path <- "~/Gavin's Dropbox/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/COBECOS/"
#path<-("~/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/COBECOS/")
setwd(path)
source("COBECOScode_v1.r")
library(optimx)

## Initialize new COBECOS object called fish
fish<-new("COBECOS", path = path, fitinfo = TRUE, graphics = FALSE)

## Initial Bioeconomic Parameters
fish@r <- 0.2 #Set value of carrying capacity
fish@K <- 1000 #Set value of intrinsic growth rate
fish@Price <- 1  	#  price of harvest    (p)
fish@FCost <- 0.5	#  fishing cost parameter  (c)
B0 = 150 #  the initial biomass (x)
fish@Biomass <- B0	#  the initial biomass (x)
fish@Fine <- 1		#  the level of fining  (f)
fish@ShadowVB <- 0.2	# The shadow value of biomass
fish@TAC <- 0.25 * fish@K * fish@r / 4 #Set value of TAC at fraction MSY

r = 0.05 # Discount rate
#fish@Effort <- 0.1246 # Set enforcement effort

## NPV
NPV = function(r,P)
{
  pv=vector()
  for (i in 1:length(P))
  {
    pv[i] = P[i] / (1+r)^i
  }
  return(sum(pv))
}


### Dynamics
horizon = 50

# Define stock growth function
fish@stockGrowth <- function (x0,q)
{
  x1 = x0 + fish@r * x0 - fish@r * x0^2 / fish@K - q
  if (x1<0) x1=0
  return (x1)
}

# Define benefit functions

fish@PrivateBFunc <- function(Harvest,Price,FCost,Biomass,ExpFine,UnitTax,TAC)
{
  if (Harvest>TAC)
    {
      Price*Harvest-FCost*((Harvest*Harvest)/Biomass)-(ExpFine*(Harvest-TAC))
    }
  else
  {
    Price*Harvest-FCost*((Harvest*Harvest)/Biomass)
  }
}

fish@SocialBFunc <- function(Price,ShadowVB,Harvest,FCost,Biomass,TotalCost,ExpFine,UnitTax,TAC)
  (Price-ShadowVB)*Harvest-FCost*((Harvest*Harvest)/Biomass)-TotalCost

######################
## Static Enforcement
######################


enfVec = seq(from=0.05,to=0.5,by=0.05)
enfNum=length(enfVec)

biomass = matrix(data=NA,nrow=horizon,ncol=enfNum)
harvest = matrix(data=NA,nrow=horizon,ncol=enfNum)
socialBenefit = matrix(data=NA,nrow=horizon,ncol=enfNum)
privateBenefit = matrix(data=NA,nrow=horizon,ncol=enfNum)
effort = matrix(data=NA,nrow=horizon,ncol=enfNum)
expfine = matrix(data=NA,nrow=horizon,ncol=enfNum)
npvPrivate = vector()
npvSocial = vector()

for (j in 1:length(enfVec))
{
  fish@Biomass = B0
  for (i in 1:horizon)
  {
    fish@Effort = enfVec[j] # Set enforcement effort
    fish <- BCalc(fish,Optimize = FALSE,stoch = 0,graphics = FALSE)
    biomass[i,j] = fish@stockGrowth(fish@Biomass,fish@Harvest)
    fish@Biomass = biomass[i,j]
    harvest[i,j] = fish@Harvest
    socialBenefit[i,j] = fish@SocialB
    privateBenefit[i,j] = fish@PrivateB
    effort[i,j] = fish@Effort
    expfine[i,j] = fish@ExpFine
  }
  npvPrivate[j]=NPV(r,privateBenefit[,j])
  npvSocial[j]=NPV(r,socialBenefit[,j])
}

par(mfrow=c(3,3),oma = c(0, 0, 4, 0))

eff<-seq(fish@Effort_Min,fish@Effort_Max,length.out=100)
plot(eff,fish@EPFitFunc(eff,fish),xlab = "Enforcement effort",ylab="Probability",type = "l", col = "red", pch=19,xlim=c(fish@Effort_Min,fish@Effort_Max),ylim=c(0,1),lwd=2)
try(points(fish@EPData$Effort,fish@EPData$Prob,pch=21),silent=TRUE)
title(paste(fish@Ename,"EP relationship"))

plot(eff,fish@ECFitFunc(eff,fish),xlab = "Enforcement effort",ylab="Cost",type = "l", col = "red", pch=19,xlim=c(fish@Effort_Min,fish@Effort_Max),lwd=2)
try(points(fish@ECData$Effort,fish@ECData$Cost,pch=21),silent=TRUE)
title(paste(fish@Ename,"EC relationship"))      

matplot(effort,type="l",xlab="Time (years)",ylab = "Enforcement Effort")
matplot(biomass,type="l",xlab="Time (years)",ylab = "Biomass")
matplot(socialBenefit,type="l",xlab="Time (years)",ylab = "Social Benefit")
matplot(privateBenefit,type="l",xlab="Time (years)",ylab = "Private Benefit")
matplot(harvest,type="l",xlab="Time (years)",ylab = "Harvest")
plot(npvSocial~enfVec,type="l",xlab="Enforcement Effort",ylab = "Social NPV",col="red",lwd=2.5)
plot(npvPrivate~enfVec,type="l",col="blue",xlab="Enforcement Effort",ylab = "Private NPV",lwd=2.5)
mtext(paste("r = ",fish@r,";  K =",fish@K,";  B0 =",B0, ";  ShadowVB =",fish@ShadowVB,"\n Price =",fish@Price,";  Cost =",fish@FCost,";  Fine =",fish@Fine,";  TAC =",fish@TAC,"; Discount Rate =",r,"\n Static Enforcement"),outer=TRUE)

######################
## Dynamic Enforcement
######################

enfVec = seq(from=0.05,to=0.5,by=0.05)
enfNum=length(enfVec)
enfSlope=-.005

biomass = matrix(data=NA,nrow=horizon,ncol=enfNum)
harvest = matrix(data=NA,nrow=horizon,ncol=enfNum)
socialBenefit = matrix(data=NA,nrow=horizon,ncol=enfNum)
privateBenefit = matrix(data=NA,nrow=horizon,ncol=enfNum)
effort = matrix(data=NA,nrow=horizon,ncol=enfNum)
expfine = matrix(data=NA,nrow=horizon,ncol=enfNum)
npvPrivate = vector()
npvSocial = vector()

for (j in 1:length(enfVec))
{
  fish@Biomass = B0
  for (i in 1:horizon)
  {
    if (enfVec[j] + enfSlope*i < 0) fish@Effort = 0 else fish@Effort = enfVec[j] + enfSlope*i # Set enforcement effort
    fish <- BCalc(fish,Optimize = FALSE,stoch = 0,graphics = FALSE)
    biomass[i,j] = fish@stockGrowth(fish@Biomass,fish@Harvest)
    fish@Biomass = biomass[i,j]
    harvest[i,j] = fish@Harvest
    socialBenefit[i,j] = fish@SocialB
    privateBenefit[i,j] = fish@PrivateB
    effort[i,j] = fish@Effort
    expfine[i,j] = fish@ExpFine
  }
  npvPrivate[j]=NPV(r,privateBenefit[,j])
  npvSocial[j]=NPV(r,socialBenefit[,j])
}

par(mfrow=c(3,3),oma = c(0, 0, 4, 0))

eff<-seq(fish@Effort_Min,fish@Effort_Max,length.out=100)
plot(eff,fish@EPFitFunc(eff,fish),xlab = "Enforcement effort",ylab="Probability",type = "l", col = "red", pch=19,xlim=c(fish@Effort_Min,fish@Effort_Max),ylim=c(0,1),lwd=2)
try(points(fish@EPData$Effort,fish@EPData$Prob,pch=21),silent=TRUE)
title(paste(fish@Ename,"EP relationship"))

plot(eff,fish@ECFitFunc(eff,fish),xlab = "Enforcement effort",ylab="Cost",type = "l", col = "red", pch=19,xlim=c(fish@Effort_Min,fish@Effort_Max),lwd=2)
try(points(fish@ECData$Effort,fish@ECData$Cost,pch=21),silent=TRUE)
title(paste(fish@Ename,"EC relationship"))      

matplot(effort,type="l",xlab="Time (years)",ylab = "Enforcement Effort")
matplot(biomass,type="l",xlab="Time (years)",ylab = "Biomass")
matplot(socialBenefit,type="l",xlab="Time (years)",ylab = "Social Benefit")
matplot(privateBenefit,type="l",xlab="Time (years)",ylab = "Private Benefit")
matplot(harvest,type="l",xlab="Time (years)",ylab = "Harvest")
plot(npvSocial~enfVec,type="l",xlab="Enforcement Effort",ylab = "Social NPV",col="red",lwd=2.5)
plot(npvPrivate~enfVec,type="l",col="blue",xlab="Enforcement Effort",ylab = "Private NPV",lwd=2.5)
mtext(paste("r = ",fish@r,";  K =",fish@K,";  B0 =",B0, ";  ShadowVB =",fish@ShadowVB,"\n Price =",fish@Price,";  Cost =",fish@FCost,";  Fine =",fish@Fine,";  TAC =",fish@TAC,"; Discount Rate =",r,"\n Dynamic Enforcement Negative"),outer=TRUE)

enfVec = seq(from=0.05,to=0.5,by=0.05)
enfNum=length(enfVec)
enfSlope=.005

biomass = matrix(data=NA,nrow=horizon,ncol=enfNum)
harvest = matrix(data=NA,nrow=horizon,ncol=enfNum)
socialBenefit = matrix(data=NA,nrow=horizon,ncol=enfNum)
privateBenefit = matrix(data=NA,nrow=horizon,ncol=enfNum)
effort = matrix(data=NA,nrow=horizon,ncol=enfNum)
expfine = matrix(data=NA,nrow=horizon,ncol=enfNum)
npvPrivate = vector()
npvSocial = vector()

for (j in 1:length(enfVec))
{
  fish@Biomass = B0
  for (i in 1:horizon)
  {
    if (enfVec[j] + enfSlope*i < 0) fish@Effort = 0 else fish@Effort = enfVec[j] + enfSlope*i # Set enforcement effort
    fish <- BCalc(fish,Optimize = FALSE,stoch = 0,graphics = FALSE)
    biomass[i,j] = fish@stockGrowth(fish@Biomass,fish@Harvest)
    fish@Biomass = biomass[i,j]
    harvest[i,j] = fish@Harvest
    socialBenefit[i,j] = fish@SocialB
    privateBenefit[i,j] = fish@PrivateB
    effort[i,j] = fish@Effort
    expfine[i,j] = fish@ExpFine
  }
  npvPrivate[j]=NPV(r,privateBenefit[,j])
  npvSocial[j]=NPV(r,socialBenefit[,j])
}

par(mfrow=c(3,3),oma = c(0, 0, 4, 0))

eff<-seq(fish@Effort_Min,fish@Effort_Max,length.out=100)
plot(eff,fish@EPFitFunc(eff,fish),xlab = "Enforcement effort",ylab="Probability",type = "l", col = "red", pch=19,xlim=c(fish@Effort_Min,fish@Effort_Max),ylim=c(0,1),lwd=2)
try(points(fish@EPData$Effort,fish@EPData$Prob,pch=21),silent=TRUE)
title(paste(fish@Ename,"EP relationship"))

plot(eff,fish@ECFitFunc(eff,fish),xlab = "Enforcement effort",ylab="Cost",type = "l", col = "red", pch=19,xlim=c(fish@Effort_Min,fish@Effort_Max),lwd=2)
try(points(fish@ECData$Effort,fish@ECData$Cost,pch=21),silent=TRUE)
title(paste(fish@Ename,"EC relationship"))      

matplot(effort,type="l",xlab="Time (years)",ylab = "Enforcement Effort")
matplot(biomass,type="l",xlab="Time (years)",ylab = "Biomass")
matplot(socialBenefit,type="l",xlab="Time (years)",ylab = "Social Benefit")
matplot(privateBenefit,type="l",xlab="Time (years)",ylab = "Private Benefit")
matplot(harvest,type="l",xlab="Time (years)",ylab = "Harvest")
plot(npvSocial~enfVec,type="l",xlab="Enforcement Effort",ylab = "Social NPV",col="red",lwd=2.5)
plot(npvPrivate~enfVec,type="l",col="blue",xlab="Enforcement Effort",ylab = "Private NPV",lwd=2.5)
mtext(paste("r = ",fish@r,";  K =",fish@K,";  B0 =",B0, ";  ShadowVB =",fish@ShadowVB,"\n Price =",fish@Price,";  Cost =",fish@FCost,";  Fine =",fish@Fine,";  TAC =",fish@TAC,"; Discount Rate =",r,"\n Dynamic Enforcement Positive"),outer=TRUE)


# enfOpt = function(enf)
# {
#   biomass = vector()
#   harvest = vector()
#   socialBenefit = vector()
#   privateBenefit = vector()
#   effort = vector()
#   expfine = vector()
#   fish@Biomass = B0  
# 
#   for (i in 1:horizon)
#   {
#     fish@Effort = enf[i] # Set enforcement effort
#     fish <- BCalc(fish,Optimize = FALSE,stoch = 0,graphics = FALSE)
#     biomass[i] = fish@stockGrowth(fish@Biomass,fish@Harvest)
#     fish@Biomass = biomass[i]
#     harvest[i] = fish@Harvest
#     socialBenefit[i] = fish@SocialB
#     privateBenefit[i] = fish@PrivateB
#     effort[i] = fish@Effort
#     expfine[i] = fish@ExpFine
#   }
#   npvPrivate=NPV(r,privateBenefit[i])
#   npvSocial=NPV(r,socialBenefit[i])
#   return(npvSocial)
# }
# 
# e=rep(0.1,horizon)
# enfMax = optimx(par=e, fn=enfOpt , method="L-BFGS-B", lower=0, upper=1, control=list(maximize = TRUE ))
# enfVec = coef(enfMax)
