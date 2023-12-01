## Set working directory
setwd("~/Gavin's Dropbox/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")
#setwd("~/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")
library(ggplot2)
library(reshape)
## Read in optimization functions
source("model19.R")

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
Q0 = 63*120*34
K = 4 * Q0/r # Assume Q0 is MSY

## Economic Parameters
price = 8.954
cost0 = 97.79 # $/trip at time 0
costTotal0 = cost0 * 120 * 34 # $ Total cost at time 0
cost = costTotal0 / (Q0/(K/2))^2 # $
# fine0 = 74.07 # $/trip at time 0
# fineTotal0 = fine0 * 120 * 34 # $ Total fine at time 0
# fine=fineTotal0 / (Q0/(K/2))^2 # $
discount = 0.05
T = 50

## Paramter Plots

# Cost of Enforcement
par(mfrow=c(1,1),mar=c(5,5,5,5))
curve(costEnforcementFunc,from=0,to=1,xlab="Enforcement Effort",ylab="Variable Enforcement Cost (USD)",lwd=2)
## Tourism Revenue
par(mfrow=c(1,1),mar=c(5,5,5,5))
curve(tourismRevenueOptimal,from=0,to=K,xlab="Biomass (kg)",ylab="Tourism Revenue (USD)",lwd=2)
## Probability of receiving a fine
par(mfrow=c(1,1),mar=c(5,5,5,5))
curve(fineProbFunc,from=0,to=1,xlab="Enforcement Effort",ylab="Probability of Receiving a Fine",lwd=2)

## Optimization parameters
bvec = seq(0,2,0.1) # Define vector of B/BMSY for optimal harvest and enforcement functions
tolF = .01 # Tolerance value for fishing effort optimization
tolE = .01 # Tolerance value for enforcement effort optimization

######################
## LOOPED MODEL RUNS ##
######################

# Loop over various archetypes
modelRuns=3
bvecMat=matrix(bvec)
tMat=matrix(seq(0,T-1))
fLegalOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
enfOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
fIUUOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))

fLegal = matrix(0,nrow=modelRuns,ncol=T)
fIUU = matrix(0,nrow=modelRuns,ncol=T)
biomass = matrix(0,nrow=modelRuns,ncol=T)
internalFinancing = matrix(0,nrow=modelRuns,ncol=T)
enfCost = matrix(0,nrow=modelRuns,ncol=T)
tourismRev= matrix(0,nrow=modelRuns,ncol=T)
fishingRev= matrix(0,nrow=modelRuns,ncol=T)
enfRev= matrix(0,nrow=modelRuns,ncol=T)
licensingRev=matrix(0,nrow=modelRuns,ncol=T)
enfEffort= matrix(0,nrow=modelRuns,ncol=T)
QIUU= matrix(0,nrow=modelRuns,ncol=T)
QLegal= matrix(0,nrow=modelRuns,ncol=T)
breakEven = rep(0,modelRuns)
profitIndustry=matrix(0,nrow=modelRuns,ncol=T)

for (i in 1:modelRuns){
  print(i)
  pol = 3
  arch = i
  if (arch==2) pol =0
  B_B0=0.5
  Beta=0
  Fine=cost*10
  Tax=0.15
  wTax=0.25
  Fee=18
  results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee,wTax=wTax)
  fLegalOpt[i,] = results$Policy
  enfOpt[i,] = results$Enforcement
  enfEffort[i,]=results$enfEffort
  fIUUOpt[i,] = results$fIUU
  
  fLegal[i,] = results$QLegal/results$biomass
  biomass[i,] = results$biomass/(K/2)
  fIUU[i,] = results$QIUU/(results$biomass-results$QLegal)
  internalFinancing[i,] = results$internalFinancing
  enfCost[i,] = results$enfCost
  tourismRev[i,]=results$tourismRev
  fishingRev[i,]=results$fishingRev
  licensingRev[i,]=results$licensingRev
  enfRev[i,]=results$enfRev
  QIUU[i,]=results$QIUU
  QLegal[i,]=results$QLegal
  profitIndustry[i,]=results$profitIndustry
  breakEven[i]=min(which(results$npvCurrent>0))
}

# B/BMSY Space Plots
par(mfrow=c(3,2),mar=c(5,5,4,8))
matplot(bvec,t(fLegalOpt),col="black",type="l",lwd=2,xlab="B/BMSY",ylab="Optimal Fishing Effort")
mtext("a) Optimal Fishing Effort in B/BMSY Space",line=1)
legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))

plot.new()
matplot(bvec,t(enfOpt),col="black",type="l",lwd=2,xlab="B/BMSY",ylab="Optimal Enforcement Effort")
mtext("b) Optimal Enforcement Effort in B/BMSY Space",line=1)
legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))

plot.new()
matplot(bvec,t(fIUUOpt),col="black",type="l",lwd=2,xlab="B/BMSY",ylab="Illegal Fishing Effort")
mtext("c) Illegal Fishing Effort in B/BMSY Space",line=1)
legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))

## Projection Plots
par(mfrow=c(3,2),mar=c(5,5,4,8))
matplot(tMat,t(biomass),col="black",type="l",lwd=2,xlab="Time",ylab="B/BMSY")
mtext("Projected B/BMSY",line=1)
legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))

matplot(tMat,t(enfEffort),col="black",type="l",lwd=2,xlab="Time",ylab="Enforcement Effort")
mtext("Projected Enforcement Effort",line=1)
legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))

matplot(tMat,t(fLegal),col="black",type="l",lwd=2,xlab="Time",ylab="Legal Fishing Effort",ylim=range(fLegal,fIUU))
mtext("Projected Legal Fishing Effort",line=1)
legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))

matplot(tMat,t(fIUU),col="black",type="l",lwd=2,xlab="Time",ylab="Illegal Fishing Effort",ylim=range(fLegal,fIUU))
mtext("Projected Illegal Fishing Effort",line=1)
legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))

matplot(tMat,t(profitIndustry),col="black",type="l",lwd=2,xlab="Time",ylab="Industry Profit (USD)",ylim=range(profitIndustry,internalFinancing-enfCost))
mtext("Projected Industry Profit",line=1)
legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))
abline(xpd=FALSE,h=0)

matplot(tMat,t(internalFinancing-enfCost),col="black",type="l",lwd=2,xlab="Time",ylab="Social Planner Profit (USD)",ylim=range(profitIndustry,internalFinancing-enfCost))
mtext("Projected Social Planner Profit",line=1)
legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))
abline(xpd=FALSE,h=0)

financingMatrix=matrix(0,nrow=modelRuns,ncol=5)
breakEvenNPVSocial=vector()
for (i in 1:modelRuns)
{
  financingMatrix[i,1] = -NPV(discount,enfCost[i,])
  financingMatrix[i,2] = NPV(discount,enfRev[i,])
  financingMatrix[i,3] = NPV(discount,fishingRev[i,])
  financingMatrix[i,4] = NPV(discount,tourismRev[i,])
  financingMatrix[i,5] = NPV(discount,licensingRev[i,])
  breakEvenNPVSocial[i]=breakevenNPV(discount,internalFinancing[i,]-enfCost[i,])
}

par(mfrow=c(1,1),mar=c(4,4,4,14),xpd=FALSE)
matplot(tMat,t(internalFinancing-enfCost),xlim=c(0,20),col="black",type="l",lwd=3,xlab="Time",ylab="Social Planner Profit (USD)")
mtext("Projected Social Planner Profits \n and Break Even Points",line=1)
legend("right",inset=-.5,xpd=NA,title="Archetype",legend=c(paste(1,"; Break Even at",breakEvenNPVSocial[1],"Years"),paste(2,"; Break Even at",breakEvenNPVSocial[2],"Years"),paste(3,"; Break Even at",breakEvenNPVSocial[3],"Years")),col="black",lwd=3,lty=c(1,2,3))
abline(v=breakEvenNPVSocial[1],lty=1)
abline(v=breakEvenNPVSocial[2],lty=2)
abline(v=breakEvenNPVSocial[3],lty=3)
abline(h=0,lty=6)

par(mfrow=c(1,1),mar=c(4,4,4,14))
barplot(t(financingMatrix),xlab="Archetype",beside=TRUE,ylab="NPV (USD)")
mtext("NPV of Revenues and Costs",line=1)
axis(side=1,at=c(3,9,15),labels=c(1,2,3))
legend("right",inset=-.45,xpd=NA,fill=gray.colors(5),legend=c("Enforcement Cost","Enforcement Fine Revenue","Landings Tax Revenue","Tourism Tax Revenue","Licensning Fee Revenue"))

for (i in 1:3)
{
  d=cbind(enfCost[i,],tourismRev[i,],fishingRev[i,],enfRev[i,],licensingRev[i,])
  colnames(d)=c("Enforcement Cost","Tourism Revenue","Landings Tax Revenue","Enforcement Fine Revenue","Licensning Fee Revenue")
  md=melt(d,id=1)
  colnames(md)=c("X1","Financing_Stream","value")
  p0=ggplot(md, aes(x=X1,y=value,group=Financing_Stream,fill=Financing_Stream)) + 
    geom_area(position="fill") +
    scale_fill_grey() +
#    scale_colour_manual(c("Enforcement Cost"=gray.colors(5)[1],"Tourism Revenue"=gray.colors(5)[2],"Landings Tax Revenue"=gray.colors(5)[3],"Enforcement Fine Revenue"=gray.colors(5)[4],"Licensning Fee Revenue"=gray.colors(5)[5]))
    theme_bw() +
    geom_hline(aes(yintercept=0.5)) + 
    labs(x="Year",y="Percentage of Total Cost and Revenue",title=paste("Cost and Revenue Streams for Arcehtype",i)) +
    xlim(0,15)
  plot(p0)
}

# Loop over various landing taxes

modelRuns=10
bvecMat=matrix(bvec)
tMat=matrix(seq(0,T-1))
fLegalOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
enfOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
fIUUOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))

fLegal = matrix(0,nrow=modelRuns,ncol=T)
fIUU = matrix(0,nrow=modelRuns,ncol=T)
biomass = matrix(0,nrow=modelRuns,ncol=T)
internalFinancing = matrix(0,nrow=modelRuns,ncol=T)
enfCost = matrix(0,nrow=modelRuns,ncol=T)
tourismRev= matrix(0,nrow=modelRuns,ncol=T)
fishingRev= matrix(0,nrow=modelRuns,ncol=T)
licensingRev=matrix(0,nrow=modelRuns,ncol=T)
enfRev= matrix(0,nrow=modelRuns,ncol=T)
breakEven = rep(0,modelRuns)
enfEffort= matrix(0,nrow=modelRuns,ncol=T)
npvSocialPlanner=vector()
npvIndustry=vector()
breakEvenSocialPlanner=vector()

interval = 0.025

for (i in 1:modelRuns){
  print(i)
  pol = 3
  arch = 3
  if (arch==2) pol =0
  B_B0=0.5
  Beta=0
  Fine=cost*10
  Tax=interval*i
  wTax=.15
  Fee=18
  results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee,wTax=wTax)
  fLegalOpt[i,] = results$Policy
  enfOpt[i,] = results$Enforcement
  fIUUOpt[i,] = results$fIUU
  
  fLegal[i,] = results$QLegal/results$biomass
  fIUU[i,] = results$QIUU/(results$biomass-results$QLegal)
  biomass[i,] = results$biomass/(K/2)
  internalFinancing[i,] = results$internalFinancing
  enfCost[i,] = results$enfCost
  enfEffort[i,]=results$enfEffort
  tourismRev[i,]=results$tourismRev
  fishingRev[i,]=results$fishingRev
  licensingRev[i,]=results$licensingRev
  enfRev[i,]=results$enfRev
  breakEven[i]=min(which(results$npvCurrent>0))
  npvIndustry[i]=NPV(discount,results$profitIndustry)
  npvSocialPlanner[i]=NPV(discount,results$internalFinancing-results$enfCost)
  breakEvenSocialPlanner[i]=breakevenNPV(discount,results$internalFinancing-results$enfCost)
}

par(mfrow=c(2,1),mar=c(4,4,3,13))
vec=seq(interval,interval*modelRuns,interval)*100
plot(vec,npvIndustry,ylab="NPV (USD)",xlab="Landings Tax (Percentage of Revenue)",type="l",lty=1,lwd=2,ylim=range(npvIndustry,npvSocialPlanner))
par(new=T)
plot(vec,npvSocialPlanner,ylab="",xlab="",type="l",lty=2,lwd=2,ylim=range(npvIndustry,npvSocialPlanner))
legend("right",inset=-.4,xpd=NA,legend=c("Industry NPV","Social Planner NPV"),col="black",lwd=c(2,2),lty=c(1,2))
plot(vec,breakEvenSocialPlanner,ylab="Social Planner Break Even Point (years)",xlab="Landings Tax (Percentage of Revenue)",type="l",lty=1,lwd=2)

# Loop over various enforcement fines

modelRuns=10
bvecMat=matrix(bvec)
tMat=matrix(seq(0,T-1))
fLegalOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
enfOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
fIUUOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))

fLegal = matrix(0,nrow=modelRuns,ncol=T)
fIUU = matrix(0,nrow=modelRuns,ncol=T)
biomass = matrix(0,nrow=modelRuns,ncol=T)
internalFinancing = matrix(0,nrow=modelRuns,ncol=T)
enfCost = matrix(0,nrow=modelRuns,ncol=T)
tourismRev= matrix(0,nrow=modelRuns,ncol=T)
fishingRev= matrix(0,nrow=modelRuns,ncol=T)
licensingRev=matrix(0,nrow=modelRuns,ncol=T)
enfRev= matrix(0,nrow=modelRuns,ncol=T)
breakEven = rep(0,modelRuns)
enfEffort= matrix(0,nrow=modelRuns,ncol=T)
npvSocialPlanner=vector()
npvIndustry=vector()
breakEvenSocialPlanner=vector()

interval = 1.5

for (i in 1:modelRuns){
  print(i)
  pol = 3
  arch = 3
  if (arch==2) pol =0
  B_B0=0.5
  Beta=0
  Fine=i*interval*cost
  Tax=.15
  wTax=.15
  Fee=18
  results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee,wTax=wTax)
  fLegalOpt[i,] = results$Policy
  enfOpt[i,] = results$Enforcement
  fIUUOpt[i,] = results$fIUU
  
  fLegal[i,] = results$QLegal/results$biomass
  fIUU[i,] = results$QIUU/(results$biomass-results$QLegal)
  biomass[i,] = results$biomass/(K/2)
  internalFinancing[i,] = results$internalFinancing
  enfCost[i,] = results$enfCost
  enfEffort[i,]=results$enfEffort
  tourismRev[i,]=results$tourismRev
  fishingRev[i,]=results$fishingRev
  licensingRev[i,]=results$licensingRev
  enfRev[i,]=results$enfRev
  breakEven[i]=min(which(results$npvCurrent>0))
  npvIndustry[i]=NPV(discount,results$profitIndustry)
  npvSocialPlanner[i]=NPV(discount,results$internalFinancing-results$enfCost)
  breakEvenSocialPlanner[i]=breakevenNPV(discount,results$internalFinancing-results$enfCost)
}
par(mfrow=c(2,1),mar=c(4,4,3,13))
vec=seq(interval,modelRuns*interval,interval)
plot(vec,npvIndustry,ylab="NPV (USD)",xlab="Enforcement Fine (Scalar of Fishing Cost)",type="l",lty=1,lwd=2,ylim=range(npvIndustry,npvSocialPlanner))
par(new=T)
plot(vec,npvSocialPlanner,ylab="",xlab="",type="l",lty=2,lwd=2,ylim=range(npvIndustry,npvSocialPlanner))
legend("right",inset=-.4,xpd=NA,legend=c("Industry NPV","Social Planner NPV"),col="black",lwd=c(2,2),lty=c(1,2))
plot(vec,breakEvenSocialPlanner,ylab="Social Planner Break Even Point (years)",xlab="Enforcement Fine (Scalar of Fishing Cost)",type="l",lty=1,lwd=2)

# Loop over various tourism taxes

modelRuns=10
bvecMat=matrix(bvec)
tMat=matrix(seq(0,T-1))
fLegalOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
enfOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
fIUUOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))

fLegal = matrix(0,nrow=modelRuns,ncol=T)
fIUU = matrix(0,nrow=modelRuns,ncol=T)
biomass = matrix(0,nrow=modelRuns,ncol=T)
internalFinancing = matrix(0,nrow=modelRuns,ncol=T)
enfCost = matrix(0,nrow=modelRuns,ncol=T)
tourismRev= matrix(0,nrow=modelRuns,ncol=T)
fishingRev= matrix(0,nrow=modelRuns,ncol=T)
licensingRev=matrix(0,nrow=modelRuns,ncol=T)
enfRev= matrix(0,nrow=modelRuns,ncol=T)
breakEven = rep(0,modelRuns)
enfEffort= matrix(0,nrow=modelRuns,ncol=T)
npvSocialPlanner=vector()
npvIndustry=vector()
breakEvenSocialPlanner=vector()

interval = 0.025

for (i in 1:modelRuns){
  print(i)
  pol = 3
  arch = 3
  if (arch==2) pol =0
  B_B0=0.5
  Beta=0
  Fine=cost*10
  Tax=0.15
  wTax=interval*i
  Fee=18
  results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee,wTax=wTax)
  fLegalOpt[i,] = results$Policy
  enfOpt[i,] = results$Enforcement
  fIUUOpt[i,] = results$fIUU
  
  fLegal[i,] = results$QLegal/results$biomass
  fIUU[i,] = results$QIUU/(results$biomass-results$QLegal)
  biomass[i,] = results$biomass/(K/2)
  internalFinancing[i,] = results$internalFinancing
  enfCost[i,] = results$enfCost
  enfEffort[i,]=results$enfEffort
  tourismRev[i,]=results$tourismRev
  fishingRev[i,]=results$fishingRev
  licensingRev[i,]=results$licensingRev
  enfRev[i,]=results$enfRev
  breakEven[i]=min(which(results$npvCurrent>0))
  npvIndustry[i]=NPV(discount,results$profitIndustry)
  npvSocialPlanner[i]=NPV(discount,results$internalFinancing-results$enfCost)
  breakEvenSocialPlanner[i]=breakevenNPV(discount,results$internalFinancing-results$enfCost)
}
par(mfrow=c(2,1),mar=c(4,4,3,13))
vec=seq(interval,interval*modelRuns,interval)*100
plot(vec,npvIndustry,ylab="NPV (USD)",xlab="Tourism Tax (Percentage of Revenue)",type="l",lty=1,lwd=2,ylim=range(npvIndustry,npvSocialPlanner))
par(new=T)
plot(vec,npvSocialPlanner,ylab="",xlab="",type="l",lty=2,lwd=2,ylim=range(npvIndustry,npvSocialPlanner))
legend("right",inset=-.4,xpd=NA,legend=c("Industry NPV","Social Planner NPV"),col="black",lwd=c(2,2),lty=c(1,2))
plot(vec,breakEvenSocialPlanner,ylab="Social Planner Break Even Point (years)",xlab="Tourism Tax (Percentage of Revenue)",type="l",lty=1,lwd=2)

# Loop over various license fees

modelRuns=10
bvecMat=matrix(bvec)
tMat=matrix(seq(0,T-1))
fLegalOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
enfOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
fIUUOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))

fLegal = matrix(0,nrow=modelRuns,ncol=T)
fIUU = matrix(0,nrow=modelRuns,ncol=T)
biomass = matrix(0,nrow=modelRuns,ncol=T)
internalFinancing = matrix(0,nrow=modelRuns,ncol=T)
enfCost = matrix(0,nrow=modelRuns,ncol=T)
tourismRev= matrix(0,nrow=modelRuns,ncol=T)
fishingRev= matrix(0,nrow=modelRuns,ncol=T)
licensingRev=matrix(0,nrow=modelRuns,ncol=T)
enfRev= matrix(0,nrow=modelRuns,ncol=T)
breakEven = rep(0,modelRuns)
enfEffort= matrix(0,nrow=modelRuns,ncol=T)
npvSocialPlanner=vector()
npvIndustry=vector()
breakEvenSocialPlanner=vector()

interval = 10

for (i in 1:modelRuns){
  print(i)
  pol = 3
  arch = 3
  if (arch==2) pol =0
  B_B0=0.5
  Beta=0
  Fine=cost*10
  Tax=0.15
  wTax=0.15
  Fee=price*i*interval
  results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee,wTax=wTax)
  fLegalOpt[i,] = results$Policy
  enfOpt[i,] = results$Enforcement
  fIUUOpt[i,] = results$fIUU
  
  fLegal[i,] = results$QLegal/results$biomass
  fIUU[i,] = results$QIUU/(results$biomass-results$QLegal)
  biomass[i,] = results$biomass/(K/2)
  internalFinancing[i,] = results$internalFinancing
  enfCost[i,] = results$enfCost
  enfEffort[i,]=results$enfEffort
  tourismRev[i,]=results$tourismRev
  fishingRev[i,]=results$fishingRev
  licensingRev[i,]=results$licensingRev
  enfRev[i,]=results$enfRev
  breakEven[i]=min(which(results$npvCurrent>0))
  npvIndustry[i]=NPV(discount,results$profitIndustry)
  npvSocialPlanner[i]=NPV(discount,results$internalFinancing-results$enfCost)
  breakEvenSocialPlanner[i]=breakevenNPV(discount,results$internalFinancing-results$enfCost)
}

par(mfrow=c(2,1),mar=c(4,4,3,13))
vec=seq(interval,interval*modelRuns,interval)
plot(vec,npvIndustry,ylab="NPV (USD)",xlab="License Fee (Scalar of Landings Price)",type="l",lty=1,lwd=2,ylim=range(npvIndustry,npvSocialPlanner))
par(new=T)
plot(vec,npvSocialPlanner,ylab="",xlab="",type="l",lty=2,lwd=2,ylim=range(npvIndustry,npvSocialPlanner))
legend("right",inset=-.4,xpd=NA,legend=c("Industry NPV","Social Planner NPV"),col="black",lwd=c(2,2),lty=c(1,2))
plot(vec,breakEvenSocialPlanner,ylab="Social Planner Break Even Point (years)",xlab="License Fee (Scalar of Landings Price)",type="l",lty=1,lwd=2)