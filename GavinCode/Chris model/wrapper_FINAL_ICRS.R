## Set working directory
setwd("~/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")
#setwd("~/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/Chris model")
library(ggplot2)
library(reshape)
## Read in optimization functions
source("model_FINAL.R")

## Define Constant Parameters
##
## Enforcement Parameters
enfVariable = 551113
enfFixed = 60000
enfVec=seq(0,1,.1)
enfDet=c(0,.04,.25,.38,.47,.54,.59,.64,.68,.72,.75)/3

## Tourism Parameters
alpha0 = 9.6448
alpha1 = -.003
alpha2 = .00004 * 3.05

# Biological parameters
r = 0.79787
Q0 = 63*120*34
K = 4 * Q0/r # Assume Q0 is MSY

## Economic Parameters
price = 10
cost0 = 97.79 # $/trip at time 0
costTotal0 = cost0 * 120 * 34 # $ Total cost at time 0
#cost = costTotal0 / (Q0^2/(K/2)) # $
cost = 7
# fine0 = 74.07 # $/trip at time 0
# fineTotal0 = fine0 * 120 * 34 # $ Total fine at time 0
# fine=fineTotal0 / (Q0/(K/2))^2 # $
discount = 0.05
T = 20

## Paramter Plots

# Cost of Enforcement
par(mfrow=c(1,1),mar=c(5,5,5,5))
curve(costEnforcementFunc,from=0,to=1,xlab="Enforcement Effort",ylab="Variable Enforcement Cost (USD)",lwd=2)
## Tourism Revenue
par(mfrow=c(1,1),mar=c(5,5,5,5))
curve(tourismRevenueOptimal,from=0,to=K,xlab="Biomass (kg)",ylab="Tourism Revenue (USD)",lwd=2)
## Probability of receiving a fine
jpeg('ICRS_Figures/fineProbability.jpg')
par(mfrow=c(1,1),mar=c(5,5,5,5))
curve(fineProbFunc,from=0,to=1,xlab="Enforcement Effort",ylab="Probability of Receiving a Fine",lwd=4,
      cex.lab=2.25, cex.axis=2.25, cex.main=2.25, cex.sub=2.25,,col=c("#0072B2"))
dev.off()

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
fStarOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
enfOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
fOpt = matrix(0,nrow=modelRuns,ncol=length(bvec))
enfOptZero = matrix(0,nrow=modelRuns,ncol=length(bvec))
fOptZero = matrix(0,nrow=modelRuns,ncol=length(bvec))

fStar = matrix(0,nrow=modelRuns,ncol=T)
fActual = matrix(0,nrow=modelRuns,ncol=T)
biomass = matrix(0,nrow=modelRuns,ncol=T)
internalFinancing = matrix(0,nrow=modelRuns,ncol=T)
enfCost = matrix(0,nrow=modelRuns,ncol=T)
tourismRev= matrix(0,nrow=modelRuns,ncol=T)
fishingRev= matrix(0,nrow=modelRuns,ncol=T)
enfRev= matrix(0,nrow=modelRuns,ncol=T)
licensingRev=matrix(0,nrow=modelRuns,ncol=T)
enfEffort= matrix(0,nrow=modelRuns,ncol=T)
QstarRun= matrix(0,nrow=modelRuns,ncol=T)
QRun= matrix(0,nrow=modelRuns,ncol=T)
breakEven = rep(0,modelRuns)
profitIndustry=matrix(0,nrow=modelRuns,ncol=T)

for (i in 1:modelRuns){
  print(i)
  pol = 4
  arch = i
  if (arch==2) pol =0
  B_B0=0.25
  Beta=0
  Fine=price*10
  Tax=0.05
  wTax=0.05
  Fee=18
  results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee,wTax=wTax)
  fStarOpt[i,] = results$Policy
  enfOpt[i,] = results$Enforcement
  fOpt[i,] = results$fOpt
  enfOptZero[i,] = results$eZero
  fOptZero[i,] = results$fOptZero
  
  enfEffort[i,]=results$enfEffort
  biomass[i,] = results$biomass/(K/2)
  fStar[i,] = results$Qstar/results$biomass
  fStar[i,][is.nan(fStar[i,])] = 0
  fStar[i,][is.infinite(fStar[i,])] = 0
  fActual[i,] = results$Q/results$biomass
  fActual[i,][is.nan(fActual[i,])] = 0
  fActual[i,][is.infinite(fActual[i,])] = 0
  internalFinancing[i,] = results$internalFinancing
  enfCost[i,] = results$enfCost
  tourismRev[i,]=results$tourismRev
  fishingRev[i,]=results$fishingRev
  licensingRev[i,]=results$licensingRev
  enfRev[i,]=results$enfRev
  QRun[i,]=results$Q
  QstarRun[i,]=results$Qstar
  profitIndustry[i,]=results$profitIndustry
  if(sum(results$npvCurrent>=0) == 0) breakEven[i]=NA else breakEven[i] = min(which(results$npvCurrent>=0))
}

# B/BMSY Space Plots
jpeg('ICRS_Figures/BvBMSYPlots.jpg',height=600,width=1000)
par(mfrow=c(2,3),mar=c(5,5,5,5))
matplot(bvec,t(fStarOpt),type="l",lwd=4,xlab="B/BMSY",ylab="Legal Effort Quota",ylim=range(fOpt),lty=1,col=c("#0072B2", "#D55E00", "#CC79A7"),
        cex.lab=2.25, cex.axis=2.25, cex.main=2.25, cex.sub=2.25)
mtext("Optimal Legal\nFishing Effort Quota",line=1,cex=1.5)

matplot(bvec,t(enfOpt),type="l",lwd=4,xlab="B/BMSY",ylab="Enforcement Effort",lty=1,col=c("#0072B2", "#D55E00", "#CC79A7"),
        cex.lab=2.25, cex.axis=2.25, cex.main=2.25, cex.sub=2.25)
mtext("Optimal Enforcement Effort",line=1,cex=1.5)

matplot(bvec,t(fOpt-fStarOpt),type="l",lwd=4,xlab="B/BMSY",ylab="Illegal Fishing Effort",ylim=range(fOpt),lty=1,col=c("#0072B2", "#D55E00", "#CC79A7"),
        cex.lab=2.25, cex.axis=2.25, cex.main=2.25, cex.sub=2.25)
mtext("Illegal Fishing Effort with\nOptimal Enforcement Effort",line=1,cex=1.5)

matplot(bvec,t(enfOptZero),type="l",lwd=4,xlab="B/BMSY",ylab="Enforcement Effort",lty=1,col=c("#0072B2", "#D55E00", "#CC79A7"),
        cex.lab=2.25, cex.axis=2.25, cex.main=2.25, cex.sub=2.25)
mtext("Enforcement Effort to\nEliminate Illegal Fishing",line=1,cex=1.5)

matplot(bvec,t(enfOptZero-enfOpt)*enfVariable,type="l",lwd=4,xlab="B/BMSY",ylab="Cost Differential [USD]",lty=1,col=c("#0072B2", "#D55E00", "#CC79A7"),
        cex.lab=2.25, cex.axis=2.25, cex.main=2.25, cex.sub=2.25)
mtext("Cost Differential to\nEliminate Illegal Fishing",line=1,cex=1.5)

plot.new()
legend("center",cex=2.25,xpd=NA,title="Archetype",legend=c("1 - Fishing Only","2 - Dive Tourism Only","3 - Fishing and Dive Tourism"),col=c("#0072B2", "#D55E00", "#CC79A7"),lwd=4,lty=1)
dev.off()

resultsTable = rbind(fStarOpt[,c(5,9,13,17,21)],
                              enfOpt[,c(5,9,13,17,21)],
                              (fOpt-fStarOpt)[,c(5,9,13,17,21)],
                              enfOptZero[,c(5,9,13,17,21)],
                              (enfOptZero-enfOpt)[,c(5,9,13,17,21)]*enfVariable)
write.csv(resultsTable,file="resultsTable.csv")

# ## Projection Plots
# par(mfrow=c(3,2),mar=c(5,5,4,8))
# matplot(tMat,t(biomass),col="black",type="l",lwd=2,xlab="Time",ylab="B/BMSY")
# mtext("Projected B/BMSY",line=1)
# legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))
# 
# matplot(tMat,t(enfEffort),col="black",type="l",lwd=2,xlab="Time",ylab="Enforcement Effort")
# mtext("Projected Enforcement Effort",line=1)
# legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))
# 
# matplot(tMat,t(fStar),col="black",type="l",lwd=2,xlab="Time",ylab="Fishing Effort Quota",ylim=range(fStar,fActual))
# mtext("Projected Fishing Effort Quota",line=1)
# legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))
# 
# matplot(tMat,t(fActual),col="black",type="l",lwd=2,xlab="Time",ylab="Actual Fishing Effort",ylim=range(fStar,fActual))
# mtext("Projected Actual Fishing Effort",line=1)
# legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))
# 
# matplot(tMat,t(profitIndustry),col="black",type="l",lwd=2,xlab="Time",ylab="Industry Profit (USD)",ylim=range(profitIndustry,internalFinancing-enfCost))
# mtext("Projected Industry Profit",line=1)
# legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))
# abline(xpd=FALSE,h=0)
# 
# matplot(tMat,t(internalFinancing-enfCost),col="black",type="l",lwd=2,xlab="Time",ylab="Social Planner Profit (USD)",ylim=range(profitIndustry,internalFinancing-enfCost))
# mtext("Projected Social Planner Profit",line=1)
# legend("right",inset=-.3,xpd=NA,title="Archetype",legend=seq(1,3,1),col="black",lwd=2,lty=c(1,2,3))
# abline(xpd=FALSE,h=0)

financingMatrix=matrix(0,nrow=modelRuns,ncol=5)
breakEvenNPVSocial=vector()
for (i in 1:modelRuns)
{
  financingMatrix[i,1] = -NPV(discount,enfCost[i,])
  financingMatrix[i,2] = NPV(discount,enfRev[i,])
  financingMatrix[i,3] = NPV(discount,fishingRev[i,])
  financingMatrix[i,4] = NPV(discount,tourismRev[i,])
  financingMatrix[i,5] = NPV(discount,licensingRev[i,])
  breakEvenNPVSocial[i] = breakevenNPV(discount,internalFinancing[i,]-enfCost[i,])
}
# 
# par(mfrow=c(1,1),mar=c(4,4,4,14),xpd=FALSE)
# matplot(tMat,t(internalFinancing-enfCost),xlim=c(0,20),col="black",type="l",lwd=3,xlab="Time",ylab="Social Planner Profit (USD)")
# mtext("Projected Social Planner Profits \n and Break Even Points",line=1)
# legend("right",inset=-.5,xpd=NA,title="Archetype",legend=c(paste(1,"; Break Even at",breakEvenNPVSocial[1],"Years"),paste(2,"; Break Even at",breakEvenNPVSocial[2],"Years"),paste(3,"; Break Even at",breakEvenNPVSocial[3],"Years")),col="black",lwd=3,lty=c(1,2,3))
# abline(v=breakEvenNPVSocial[1],lty=1)
# abline(v=breakEvenNPVSocial[2],lty=2)
# abline(v=breakEvenNPVSocial[3],lty=3)
# abline(h=0,lty=6)
# 
# par(mfrow=c(1,1),mar=c(4,4,4,14))
# barplot(t(financingMatrix),xlab="Archetype",beside=TRUE,ylab="NPV (USD)")
# mtext("NPV of Revenues and Costs",line=1)
# axis(side=1,at=c(3,9,15),labels=c(1,2,3))
# legend("right",inset=-.45,xpd=NA,fill=gray.colors(5),legend=c("Enforcement Cost","Enforcement Fine Revenue","Landings Tax Revenue","Tourism Tax Revenue","Licensning Fee Revenue"))
# 
# for (i in 1:3)
# {
#   d=cbind(enfCost[i,],tourismRev[i,],fishingRev[i,],enfRev[i,],licensingRev[i,])
#   colnames(d)=c("Enforcement Cost","Tourism Revenue","Landings Tax Revenue","Enforcement Fine Revenue","Licensning Fee Revenue")
#   md=melt(d,id=1)
#   colnames(md)=c("X1","Financing_Stream","value")
#   p0=ggplot(md, aes(x=X1,y=value,group=Financing_Stream,fill=Financing_Stream)) + 
#     geom_area(position="fill") +
#     scale_fill_grey() +
# #    scale_colour_manual(c("Enforcement Cost"=gray.colors(5)[1],"Tourism Revenue"=gray.colors(5)[2],"Landings Tax Revenue"=gray.colors(5)[3],"Enforcement Fine Revenue"=gray.colors(5)[4],"Licensning Fee Revenue"=gray.colors(5)[5]))
#     theme_bw() +
#     geom_hline(aes(yintercept=0.5)) + 
#     labs(x="Year",y="Percentage of Total Cost and Revenue",title=paste("Cost and Revenue Streams for Archetype",i)) +
#     xlim(0,15)
#   plot(p0)
# }

## Loop over archetype
bioRuns = 10
leverRuns = 3
NPVsocialRev = matrix(0,nrow=leverRuns,ncol=bioRuns)
NPVSocialMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
NPVenfCostMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
NPVtourismRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
NPVfishingRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
NPVenfRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
NPVlicenseRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
breakEvenMat = matrix(0,nrow=leverRuns,ncol=bioRuns)

for (i in 1:bioRuns){
  for (j in 1:leverRuns){
    print(paste("B0/BMSY run",i,"; Archetype run",j))
    pol = 4
    arch = j
    if (arch==2) pol =0
    B_B0=i/10
    Beta=0
    Fine=price*10
    Tax=0.05
    wTax=0.05
    Fee=18
    results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee,wTax=wTax)
    NPVsocialRev[j,i] = NPV(discount,results$socialRev)
    NPVenfCostMat[j,i] = NPV(discount,results$enfCost)
    NPVtourismRevMat[j,i] = NPV(discount,results$tourismRev)
    NPVfishingRevMat[j,i] = NPV(discount,results$fishingRev)
    NPVenfRevMat[j,i] = NPV(discount,results$enfRev)
    NPVlicenseRevMat[j,i] = NPV(discount,results$licensingRev)
    breakEvenMat[j,i] = breakevenNPV(discount,results$internalFinancing-results$enfCost)
  }
}

jpeg('ICRS_Figures/revenueCost.jpg',height=600,width=1000)
par(mfrow=c(2,2),mar=c(4,8,4,4))
for (i in 1:leverRuns){
  financingCost=matrix(0,nrow=bioRuns,ncol=1)
  financingRev=matrix(0,nrow=bioRuns,ncol=4)
  breakEvenVec=vector()
  for (j in 1:bioRuns){
    financingCost[j,1] = -NPVenfCostMat[i,j]
    financingRev[j,1] = NPVenfRevMat[i,j]
    financingRev[j,2] = NPVfishingRevMat[i,j]
    financingRev[j,3] = NPVtourismRevMat[i,j]
    financingRev[j,4] =NPVlicenseRevMat[i,j]
    breakEvenVec[i] = breakEvenMat[i,j]
  }
  b = barplot((t(financingCost)/-min(financingCost))[c(2,4,6,8,10)],ylim=range(-1,sum(financingRev[10,])/-min(financingCost)),col=c("#D55E00"),xlab="B0/BMSY",ylab="Revenue and Cost NPV\n(Normalized)",
              cex.lab=2.25, cex.axis=2.25, cex.main=2.25, cex.sub=2.25)
  par(new=TRUE)
  c = barplot((t(financingRev)/-min(financingCost))[,c(2,4,6,8,10)],ylim=range(-1,sum(financingRev[10,])/-min(financingCost)),col=c("#E69F00", "#56B4E9", "#009E73","#0072B2"),
              cex.lab=2.25, cex.axis=2.25, cex.main=2.25, cex.sub=2.25)
  abline(h=0)
  if (i ==1) mtext(paste("Archetype 1 - Fishing Only"),line=1,cex=2.25)
  if (i ==2) mtext(paste("Archetype 2 - Dive Tourism Only"),line=1,cex=2.25)
  if (i ==3) mtext(paste("Archetype 3 - Fishing and Dive Tourism"),line=1,cex=2.25)
  axis(side=1,at=b,labels=seq(0.4,2,.4),cex.axis=2.25)
  box()
}


plot.new()
legend("center",cex=2.25,xpd=NA,fill=rev(c("#D55E00","#E69F00", "#56B4E9", "#009E73","#0072B2")),legend=rev(c("Enforcement Cost","Enforcement Fine Revenue","Landings Tax Revenue","Tourism Tax Revenue","Licensning Fee Revenue")))
dev.off()

jpeg('ICRS_Figures/breakEven.jpg',width=800,height=400)
par(mfrow=c(1,1),mar=c(4,6,4,25))
matplot(seq(1,bioRuns)/5,t(breakEvenMat),type="l",col=c("#0072B2", "#D55E00", "#CC79A7"),lwd=4,lty=1,xlab="B0/BMSY",ylab="Break-Even Point [Years]",
        cex.lab=2.25, cex.axis=2.25, cex.main=2.25, cex.sub=2.25)
legend("right",inset=-1,cex=1.5,xpd=NA,title="Archetype",legend=c("1 - Fishing Only","2 - Dive Tourism Only","3 - Fishing and Dive Tourism"),col=c("#0072B2", "#D55E00", "#CC79A7"),lwd=4,lty=1)
dev.off()
## Financing Plots

jpeg('ICRS_Figures/financing.jpg',height=600,width=1000)
for (aRun in 3:3) # Start financing plot loop
{
  
  ## Loop over landings tax
  
  bioRuns = 5
  leverRuns = 10
  interval= 0.025
  NPVIndustryMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVSocialMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVenfCostMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVtourismRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVfishingRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVenfRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVlicenseRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  breakEvenMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  
  for (i in 1:bioRuns){
    for (j in 1:leverRuns){
      print(paste("B0/BMSY run",i,"; Financing Lever run",j))
      pol = 4
      arch = aRun
      if (arch==2) pol =0
      B_B0=i/5
      Beta=0
      Fine=price*10
      Tax=j*interval
      wTax=0.05
      Fee=18
      results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee,wTax=wTax)
      NPVIndustryMat[j,i] = NPV(discount,results$profitIndustry)
      NPVSocialMat[j,i]=NPV(discount,results$internalFinancing-results$enfCost)
      NPVenfCostMat[j,i] = NPV(discount,results$enfCost)
      NPVtourismRevMat[j,i] = NPV(discount,results$tourismRev)
      NPVfishingRevMat[j,i] = NPV(discount,results$fishingRev)
      NPVenfRevMat[j,i] = NPV(discount,results$enfRev)
      NPVlicenseRevMat[j,i] = NPV(discount,results$licensingRev)
      breakEvenMat[j,i] = breakevenNPV(discount,results$internalFinancing-results$enfCost)
    }
  }
  
  par(mfrow=c(3,3),mar=c(5,5,5,5))
  vec=seq(interval,leverRuns*interval,interval)*100
  matplot(vec,NPVIndustryMat,type="l",lty=1,lwd=4,col=c("#D55E00","#E69F00", "#56B4E9", "#009E73","#0072B2"),ylab="Industry NPV (USD)",xlab="Landings Tax (Percentage of Revenue)",
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  mtext("a) Industry NPV versus\nLandings Tax", line=1,cex=1.5)
  
  matplot(vec,breakEvenMat,type="l",lty=1,lwd=4,col=c("#D55E00","#E69F00", "#56B4E9", "#009E73","#0072B2"),ylab="Break Even Point (Years)",xlab="Landings Tax (Percentage of Revenue)",
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  mtext("b) Break-even point versus\nLandings Tax", line=1,cex=1.5)
  
  ## Loop over enforcement fines
  
  bioRuns = 5
  leverRuns = 10
  interval=5
  NPVIndustryMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVSocialMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVenfCostMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVtourismRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVfishingRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVenfRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVlicenseRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  breakEvenMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  
  for (i in 1:bioRuns){
    for (j in 1:leverRuns){
      print(paste("B0/BMSY run",i,"; Financing Lever run",j))
      pol = 4
      arch = aRun
      if (arch==2) pol =0
      B_B0=i/5
      Beta=0
      Fine=price*j*interval
      Tax=0.05
      wTax=0.05
      Fee=18
      results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee,wTax=wTax)
      NPVIndustryMat[j,i] = NPV(discount,results$profitIndustry)
      NPVSocialMat[j,i]=NPV(discount,results$internalFinancing-results$enfCost)
      NPVenfCostMat[j,i] = NPV(discount,results$enfCost)
      NPVtourismRevMat[j,i] = NPV(discount,results$tourismRev)
      NPVfishingRevMat[j,i] = NPV(discount,results$fishingRev)
      NPVenfRevMat[j,i] = NPV(discount,results$enfRev)
      NPVlicenseRevMat[j,i] = NPV(discount,results$licensingRev)
      breakEvenMat[j,i] = breakevenNPV(discount,results$internalFinancing-results$enfCost)
    }
  }
  
  vec=seq(interval,leverRuns*interval,interval)
  matplot(vec,NPVIndustryMat,type="l",lty=1,lwd=4,col=c("#D55E00","#E69F00", "#56B4E9", "#009E73","#0072B2"),ylab="Industry NPV (USD)",xlab="Enforcement Fine (Scalar of Landings Price)",
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  mtext("c) Industry NPV versus\nEnforcement Fine", line=1,cex=1.5)
  
  matplot(vec,breakEvenMat,type="l",lty=1,lwd=4,col=c("#D55E00","#E69F00", "#56B4E9", "#009E73","#0072B2"),ylab="Break Even Point (Years)",xlab="Enforcement Fine (Scalar of Landings Price)",
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  mtext("d) Break-even point versus\nEnforcement Fine", line=1,cex=1.5)
  
  ## Loop over tourism tax
  
  bioRuns = 5
  leverRuns = 10
  interval= 0.025
  NPVIndustryMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVSocialMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVenfCostMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVtourismRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVfishingRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVenfRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVlicenseRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  breakEvenMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  
  for (i in 1:bioRuns){
    for (j in 1:leverRuns){
      print(paste("B0/BMSY run",i,"; Financing Lever run",j))
      pol = 4
      arch = aRun
      if (arch==2) pol =0
      B_B0=i/5
      Beta=0
      Fine=price*10
      Tax=0.05
      wTax=j*interval
      Fee=18
      results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee,wTax=wTax)
      NPVIndustryMat[j,i] = NPV(discount,results$profitIndustry)
      NPVSocialMat[j,i]=NPV(discount,results$internalFinancing-results$enfCost)
      NPVenfCostMat[j,i] = NPV(discount,results$enfCost)
      NPVtourismRevMat[j,i] = NPV(discount,results$tourismRev)
      NPVfishingRevMat[j,i] = NPV(discount,results$fishingRev)
      NPVenfRevMat[j,i] = NPV(discount,results$enfRev)
      NPVlicenseRevMat[j,i] = NPV(discount,results$licensingRev)
      breakEvenMat[j,i] = breakevenNPV(discount,results$internalFinancing-results$enfCost)
    }
  }
  
  vec=seq(interval,leverRuns*interval,interval)*100
  matplot(vec,NPVIndustryMat,type="l",lty=1,lwd=4,col=c("#D55E00","#E69F00", "#56B4E9", "#009E73","#0072B2"),ylab="Industry NPV (USD)",xlab="Tourism Tax (Percentage of Revenue)",
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  mtext("e) Industry NPV versus\nTourism Tax", line=1,cex=1.5)
  
  matplot(vec,breakEvenMat,type="l",lty=1,lwd=4,col=c("#D55E00","#E69F00", "#56B4E9", "#009E73","#0072B2"),ylab="Break Even Point (Years)",xlab="Tourism Tax (Percentage of Revenue)",
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  mtext("f) Break-even point versus\nTourism Tax", line=1,cex=1.5)
  
  ## Loop over licsense fee
  
  bioRuns = 5
  leverRuns = 10
  interval= 10
  NPVIndustryMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVSocialMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVenfCostMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVtourismRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVfishingRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVenfRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  NPVlicenseRevMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  breakEvenMat = matrix(0,nrow=leverRuns,ncol=bioRuns)
  
  for (i in 1:bioRuns){
    for (j in 1:leverRuns){
      print(paste("B0/BMSY run",i,"; Financing Lever run",j))
      pol = 4
      arch = aRun
      if (arch==2) pol =0
      B_B0=i/5
      Beta=0
      Fine=price*4
      Tax=0.05
      wTax=0.05
      Fee=price*j*interval
      results=wrapperFunction(pol=pol,arch=arch,B_B0=B_B0,Beta=0,Fine=Fine,Tax=Tax,Fee=Fee,wTax=wTax)
      NPVIndustryMat[j,i] = NPV(discount,results$profitIndustry)
      NPVSocialMat[j,i]=NPV(discount,results$internalFinancing-results$enfCost)
      NPVenfCostMat[j,i] = NPV(discount,results$enfCost)
      NPVtourismRevMat[j,i] = NPV(discount,results$tourismRev)
      NPVfishingRevMat[j,i] = NPV(discount,results$fishingRev)
      NPVenfRevMat[j,i] = NPV(discount,results$enfRev)
      NPVlicenseRevMat[j,i] = NPV(discount,results$licensingRev)
      breakEvenMat[j,i] = breakevenNPV(discount,results$internalFinancing-results$enfCost)
    }
  }
  
  vec=seq(interval,leverRuns*interval,interval)
  matplot(vec,NPVIndustryMat,type="l",lty=1,lwd=4,col=c("#D55E00","#E69F00", "#56B4E9", "#009E73","#0072B2"),ylab="Industry NPV (USD)",xlab="Fishing Licsense Fee (scalar of landings price)",
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  mtext("g) Industry NPV versus\nLicense Fee", line=1,cex=1.5)
  
  matplot(vec,breakEvenMat,type="l",lty=1,lwd=4,col=c("#D55E00","#E69F00", "#56B4E9", "#009E73","#0072B2"),ylab="Break Even Point (Years)",xlab="Fishing Licsense Fee (scalar of landings price)",
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  mtext("h) Break-even point versus\nLicense Fee", line=1,cex=1.5)
  
  plot.new()
  legend("center",cex=2.25,xpd=NA,title="B0/BMSY",legend=rev(seq(1,bioRuns)/2.5),col=rev(c("#D55E00","#E69F00", "#56B4E9", "#009E73","#0072B2")),lwd=4,lty=1)
  
} # End financing plot loop
dev.off()
