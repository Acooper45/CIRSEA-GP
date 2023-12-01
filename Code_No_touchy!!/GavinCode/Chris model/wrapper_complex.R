## Set working directory
setwd("~/Gavin's Dropbox/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")
#setwd("~/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")

## Read in optimization functions
source("model1217.R")

## Define Constant Parameters
##
## Enforcement Parameters
enfVariable = 551113
enfFixed = 60000
enfVec=seq(0,1,.1)
enfDet=c(0,.04,.25,.38,.47,.54,.59,.64,.68,.72,.75)

## Tourism Parameters
alpha0 = 9.6448
alpha1 = -.003
alpha2 = .00004 * 2

# Biological parameters
r = 0.79787
K = 2.66e5
Q0 = 63*120*34
K = 4 * Q0/r # Assume f0 is MSY

## Economic Parameters
price = 8.954
cost0 = 97.79 # $/trip at time 0
costTotal0 = cost0 * 120 * 34 # $ Total cost at time 0
discount = 0.05
T = 50

## Optimization parameters
bvec = seq(0,2,0.1) # Define vector of B/BMSY for optimal harvest and enforcement functions
tolF = .01 # Tolerance value for fishing effort optimization
tolE = .01 # Tolerance value for enforcement effort optimization


######################
## SINGLE MODEL RUN ##
######################

# ## Model Run Parameters
# policy = 3 # 0: No fishing; 1: FMSY at all time steps; 2: 0 until biomass reaches BMSY, FMSY afterwards; 3: Optimal policy
# archetype = 3 # 1: Only fishing profit considered; 2: Only tourism profit considered; 3: Fishing and Tourism Profit Considered
# if (archetype==2) policy =0; # If only considering tourism, do optimization routine for policy 1
# B_B0=.5 # Starting biomass
# Beta=0 # Social preference to fish illegally parameter
# Fine=74.07/63*10 # $/kg @ x0
# Tax=0.15 # Percentage of price/kg
# Fee=18 # Annual fee per boat
# 
# ## Run Model
# results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee)
# 
# ## Plot Results
# 
# par(mfrow=c(3,3),oma = c(0,0,6,0) + 0.1)
# plot(results$Policy~bvec,xlab="B/BMSY",type="l",ylab="Optimal fLegal",col="red")
# plot(results$Enforcement~bvec,xlab="B/BMSY",type="l",ylab="Optimal enfEffort",col="red")
# plot(results$fIUU~bvec,xlab="B/BMSY",type="l",ylab="fIUU",col="red")
# curve(fineProbFunc,from=0,to=1,xlab="Enforcement Effort",col="blue")
# curve(costEnforcementFunc,from=0,to=1,xlab="Enforcement Effort",col="blue")
# plot(results$enfEffort,xlab="Years",type="l")
# plot(results$biomass/(K/2),xlab="Years",type="l",ylab="B/BMSY")
# plot(results$QLegal/results$biomass,xlab="Years",type="l",ylab="f_Legal")
# plot(results$QIUU/results$biomass,xlab="Years",type="l",ylab="f_Illegal")
# mtext(paste("Policy",policy,"; Archetype",archetype,"; NPV =",round(results$npv,0),"\n B/B0 = ",signif(B_B0,3),"; Beta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1,outer=TRUE)
# 


######################
## LOOPED MODEL RUN ##
######################

## Model Run Parameters

Run = "archetype"

if Run = "archetpye"{
  run=c(1,2,3)
  pol=3
  B_B0=0.5
  Beta=0
  Fine=74.07/63*10
  Tax=0.15
  Fee=18
} 

if Run = "policy" {
  run = c(1,2,3)
  arch = 3
  B_B0=0.5
  Beta=0
  Fine=74.07/63*10
  Tax=0.15
  Fee=18
}

if Run = ""
if Run = "tax" run = seq(0.05,0.25,0.05)

modelRuns=length(run)
bvecMat=matrix(bvec)
tMat=matrix(seq(1,T))
fLegalOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
enfOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
fIUUOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))

fLegal = matrix(0,nrow=modelRuns,ncol=T)
fIUU = matrix(0,nrow=modelRuns,ncol=T)
biomass = matrix(0,nrow=modelRuns,ncol=T)
internalFinancing = matrix(0,nrow=modelRuns,ncol=T)
enfCost = matrix(0,nrow=modelRuns,ncol=T)

for (i in 3:modelRuns){
  print(i)
  pol = 3
  if Run = "archetpye" arch = i else arch = arch
  if (arch==2) policy =0
  B_B0=i/10
  Beta=0
  Fine=74.07/63*10
  Tax=0.15
  Fee=18
  results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee)
  fLegalOpt[i,] = results$Policy
  enfOpt[i,] = results$Enforcement
  fIUUOpt[i,] = results$fIUU
  
  fLegal[i,] = results$QLegal/results$biomass
  fIUU[i,] = results$QIUU/results$biomass
  biomass[i,] = results$biomass/(K/2)
  internalFinancing[i,] = results$internalFinancing
  enfCost[i,] = results$enfCost
}

par(mar=c(4,4,5,7))
matplot(bvec,t(fLegalOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Fishing Effort")
legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)

par(mar=c(4,4,5,7))
matplot(bvec,t(enfOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Optimal Enforcement Effort")
legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)

par(mar=c(4,4,5,7))
matplot(bvec,t(fIUUOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Illegal Fishing Effort")
legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)

par(mar=c(4,4,5,7))
matplot(tMat,t(biomass),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="B/BMSY")
legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)

par(mar=c(4,4,5,7))
matplot(tMat,t(fLegal),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Legal Fishing Effort")
legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)

par(mar=c(4,4,5,7))
matplot(tMat,t(fIUU),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Illegal Fishing Effort")
legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)

par(mar=c(4,4,5,7))
matplot(tMat,t(internalFinancing-enfCost),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Social Planner Profit")
legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)


## Loop over various B/B0
# 
# for (i in 3:modelRuns){
#   print(i)
#   pol = 3
#   arch = 3
#   if (arch==2) policy =0
#   B_B0=i/10
#   Beta=0
#   Fine=74.07/63*10
#   Tax=0.15
#   Fee=18
#   results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee)
#   fLegalOpt[i,] = results$Policy
#   enfOpt[i,] = results$Enforcement
#   fIUUOpt[i,] = results$fIUU
#   
#   fLegal[i,] = results$QLegal/results$biomass
#   fIUU[i,] = results$QIUU/results$biomass
#   biomass[i,] = results$biomass/(K/2)
#   internalFinancing[i,] = results$internalFinancing
#   enfCost[i,] = results$enfCost
# }
# 
# par(mar=c(4,4,5,7))
# matplot(bvec,t(fLegalOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(bvec,t(enfOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Optimal Enforcement Effort")
# legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(bvec,t(fIUUOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Illegal Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(biomass),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="B/BMSY")
# legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(fLegal),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Legal Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(fIUU),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Illegal Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(internalFinancing-enfCost),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Social Planner Profit")
# legend("right",inset=-.2,xpd=NA,title="B/B0",legend=seq(.1,1,.1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,4))
# plot(results$enfCost,type="l",xlab="Years",ylab="$",col="red",lwd=4,ylim=range(c(results$internalFinancing,results$enfCost,results$fishingRev,results$tourismRev)))
# par(new=T)
# plot(results$internalFinancing,type="l",col="blue",lwd=4,ylim=range(c(results$internalFinancing,results$enfCost,results$fishingRev,results$tourismRev)),ylab = "",xlab="")
# par(new=T)
# plot(results$tourismRev,type="l",col="green",lwd=2,ylim=range(c(results$internalFinancing,results$enfCost,results$fishingRev,results$tourismRev)),ylab = "",xlab="")
# par(new=T)
# plot(results$fishingRev,type="l",col="orange",lwd=2,ylim=range(c(results$internalFinancing,results$enfCost,results$fishingRev,results$tourismRev)),ylab = "",xlab="")
# par(new=T)
# plot(results$enfRev,type="l",col="purple",lwd=2,ylim=range(c(results$internalFinancing,results$enfCost,results$fishingRev,results$tourismRev)),ylab = "",xlab="")
# 
# mtext(paste("Policy",pol,"; Archetype",arch,"\n Break Even point =",results$breakEven,"years"), side=3, line=1, cex=1, font=1)
# legend("topright",legend=c("Total Internal Financing Revenue","Tourism Revenue","Fishery Revenue","Enforcement Revenue","Enforcement Costs"),col=c("blue","green","orange","purple","red"),"l",lwd=c(4,2,2,2,4))
# abline(v=results$breakEven)

## Loop over various archetypes
# 
# for (i in 1:modelRuns){
#   print(i)
#   pol = 1
#   arch = i
#   if (arch==2) pol =0
#   B_B0=0.5
#   Beta=0
#   Fine=74.07/63*10
#   Tax=0.15
#   Fee=18
#   results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee)
#   fLegalOpt[i,] = results$Policy
#   enfOpt[i,] = results$Enforcement
#   fIUUOpt[i,] = results$fIUU
#   
#   fLegal[i,] = results$QLegal/results$biomass
#   fIUU[i,] = results$QIUU/results$biomass
#   biomass[i,] = results$biomass/(K/2)
#   internalFinancing[i,] = results$internalFinancing
#   enfCost[i,] = results$enfCost
# }
# 
# par(mar=c(4,4,5,7))
# matplot(bvec,t(fLegalOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="Archetype",legend=seq(1,3,1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(bvec,t(enfOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Optimal Enforcement Effort")
# legend("right",inset=-.2,xpd=NA,title="Archetype",legend=seq(1,3,1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(bvec,t(fIUUOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Illegal Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="Archetype",legend=seq(1,3,1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(biomass),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="B/BMSY")
# legend("right",inset=-.2,xpd=NA,title="Archetype",legend=seq(1,3,1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(fLegal),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Legal Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="Archetype",legend=seq(1,3,1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(fIUU),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Illegal Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="Archetype",legend=seq(1,3,1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(internalFinancing-enfCost),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Social Planner Profit")
# legend("right",inset=-.2,xpd=NA,title="Archetype",legend=seq(1,3,1),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"\nBeta = ",Beta,"; Tax =",signif(Tax*price,3),"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)

# # Loop over various landing taxes
# 
# for (i in 1:modelRuns){
#   print(i)
#   pol = 3
#   arch = 3
#   if (arch==2) pol =0
#   B_B0=0.5
#   Beta=0
#   Fine=74.07/63*10
#   Tax=i*0.05
#   Fee=18
#   results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee)
#   fLegalOpt[i,] = results$Policy
#   enfOpt[i,] = results$Enforcement
#   fIUUOpt[i,] = results$fIUU
#   
#   fLegal[i,] = results$QLegal/results$biomass
#   fIUU[i,] = results$QIUU/results$biomass
#   biomass[i,] = results$biomass/(K/2)
#   internalFinancing[i,] = results$internalFinancing
#   enfCost[i,] = results$enfCost
# }
# 
# par(mar=c(4,4,5,7))
# matplot(bvec,t(fLegalOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="Tax",legend=seq(0.05,0.25,0.05),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(bvec,t(enfOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Optimal Enforcement Effort")
# legend("right",inset=-.2,xpd=NA,title="Tax",legend=seq(0.05,0.25,0.05),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(bvec,t(fIUUOpt),col=rainbow(modelRuns),type="l",lwd=2,xlab="B/BMSY",ylab="Illegal Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="Tax",legend=seq(0.05,0.25,0.05),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(biomass),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="B/BMSY")
# legend("right",inset=-.2,xpd=NA,title="Tax",legend=seq(0.05,0.25,0.05),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(fLegal),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Legal Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="Tax",legend=seq(0.05,0.25,0.05),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(fIUU),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Illegal Fishing Effort")
# legend("right",inset=-.2,xpd=NA,title="Tax",legend=seq(0.05,0.25,0.05),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)
# 
# par(mar=c(4,4,5,7))
# matplot(tMat,t(internalFinancing-enfCost),col=rainbow(modelRuns),type="l",lwd=2,xlab="Time",ylab="Social Planner Profit")
# legend("right",inset=-.2,xpd=NA,title="Tax",legend=seq(0.05,0.25,0.05),col=rainbow(modelRuns),lwd=2)
# mtext(paste("Policy",pol,"; Archetype",arch,"\nBeta = ",Beta,"[$/kg] \n Price =",signif(price,3),"[$/kg] ; Cost =",signif(results$cost,3),"[$/kg @ x0] ; Fine = ",signif(Fine,3),"[$/kg @ x0]"), side=3, line=1, cex=1, font=1)

