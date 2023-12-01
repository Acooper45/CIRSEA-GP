

cat("","\n")
cat("","\n")
cat("","\n")
cat("----------------------------------------------------------------------------","\n")
cat("COBECOS   COsts and BEnefits of COntrol Strategies", "\n")
cat("Alpha v1.3","\n")
cat("April 2008","\n")
cat("","\n")
cat("T.Carruthers & C.Edwards","\n")
cat("Imperial College","\n")
cat("","\n")
cat("----------------------------------------------------------------------------","\n")


# METHOD 1
# -----------------------------------  INITIALISATION  ---------------------------------------------------------------------------
# A number of methods are defined which allow the user to investigate cost-benefit
# create a class of control object called "COBECOS"

library(graphics)
library(MASS)
setClass("COBECOS",representation(nE = "numeric", Enames = "list", EPdata = "data.frame",
ECdata = "data.frame",EPFitList = "list", ECFitList = "list",EPDataList = "list", 
ECDataList = "list", EPFitFuncs = "list", ECFitFuncs = "list", controller = "data.frame", 
Effort = "numeric", Fine = "numeric", Costs = "numeric", TotalCost = "numeric", Probs = "numeric", 
Biomass = "numeric",FCost = "numeric",Harvest = "numeric",Price = "numeric",ShadowVB= "numeric",
ExpFine = "numeric", SocialBFunc = "character", PrivateBFunc = "character",FishingRFunc = "character",
PrivateB = "numeric", SocialB = "numeric", SocialOBJFunc = "function", PrivateOBJFunc = "function",
Stochastic = "logical", Optimized = "list", StochN = "numeric", ErrCV = "numeric"), prototype(EPdata = NULL, 
ECdata = NULL, nE = 1, Biomass =100, FCost = 0.5, Price = 1, Harvest = 1, ShadowVB = 0.2, 
SocialBFunc ="default", PrivateBFunc = "default", FishingRFunc = "default",StochN = 2000,ErrCV=0.1,
Stochastic=FALSE))

# Create an initialisation method of the COBECOS control object

getobjname<-function(.Object){
    deparse(substitute(.Object))
}

setMethod("initialize", "COBECOS", function(.Object, path, fitinfo, graphics){


cat("Method 1. Initialisation of COBECOS object","\n")
flush.console()

# Read in the observed data for which functions and parameter are to be fit
.Object@controller <- read.csv(file = paste(path, "Control.csv",sep = ""),head=TRUE,sep=",")
.Object@EPdata <-read.csv(file = paste(path, "ProbData.csv",sep = ""),head=TRUE,sep=",")
.Object@EPdata<-.Object@EPdata[order(.Object@EPdata$Effort),]
.Object@ECdata <- read.csv(file = paste(path, "CostData.csv",sep = ""),head=TRUE,sep=",")
.Object@ECdata<-.Object@ECdata[order(.Object@ECdata$Effort),]
.Object@Enames <- as.list(.Object@controller$Type)
.Object@nE <- max(.Object@controller$Type.no)  # find maximum number of enforcement types

if(graphics == TRUE){
  quartz()
  par(mfrow = c(.Object@nE,2) , mai = c(0.55,0.55,0.35,0.05))  # set up the plot array
}

for(i in 1:.Object@nE){

    # Some unrelated stuff
    .Object@Effort[i] <- 0.5                           # set the default effort levels for later
    .Object@Fine[i] <- .Object@Price/.Object@nE        # set default fines for later  (this is just a silly first guess that will provide not ridiculous or error prone problems for first time users that haven't specified fines correctly yet!)
  
    cat("","\n")
    cat("","\n")                                                                     # output to screen
    cat(paste("-------------- Enforcement type: " ,.Object@Enames[[i]], " -----------------------------------"),"\n")    # output to screen
    cat("----- Probability of detection versus enforcement effort -----","\n")           # output to screen
    
    ######################
    # Effort-Prob models #
    ######################
    if(is.na(.Object@controller$EPMod[i])){  # Fit effort - prob models - only do this if the user hasn't specified a model for this effort probability model ie they'v stuck an "NA" in the EPMod column of the control.csv excel file
        cat("Model fitted to data","\n")        # output to screen
        # Fit models models to the EPdata
        # select dist for probability of detection modelling on the basis of AIC
        .Object@EPDataList[[i]]<- subset(.Object@EPdata,Type.no==.Object@controller$Type.no[[i]],select=c(Effort,Prob))

        n<-length(.Object@EPDataList[[i]]$Prob)
        
        logistic <- function(params){
            pred <- params[1]/(1+exp((params[2]-.Object@EPDataList[[i]]$Effort)/params[3]))
            mahalanobis(.Object@EPDataList[[i]]$Prob,pred,diag(length(pred)))
        }
        
        exponenti <- function(params){
            pred <- (1-exp(-params[1]*.Object@EPDataList[[i]]$Effort))*(1-params[2])+params[2]
            mahalanobis(.Object@EPDataList[[i]]$Prob,pred,diag(length(pred)))
        }
        
        exponent <- function(params){
            pred <- 1-exp(-params[1]*.Object@EPDataList[[i]]$Effort)
            mahalanobis(.Object@EPDataList[[i]]$Prob,pred,diag(length(pred)))
        }
        
        params<-c(0.95,0.5,0.05)
        templogistic <- optim(params, logistic ,method="Nelder-Mead", hessian = TRUE)
        if(templogistic$convergence!=0){      # try some alternative methods if it doesn't work
             params<-c(0.95,0.5,0.05)
             templogistic <- optim(params, logistic ,method="BFGS", hessian = TRUE)
             cat("NOTE: The optimiser could not converge satisfactorily using the Nelder-Mead method. It will now try the BFGS method", "\n")   
             if(templogistic$convergence!=0) cat("ERROR: The optimiser could not converge")
        }
        AIClogistic <- n*log(templogistic$value/n)+6
        
        params<-c(mean(-(log(1-.Object@EPDataList[[i]]$Prob)/.Object@EPDataList[[i]]$Effort)),min(.Object@EPDataList[[i]]$Prob))
           # this is the initial guess at the exponential parameter 'a' (not nearly an MLE but not so far off!)
        tempexponenti <- optim(params, exponenti ,method="BFGS", hessian = TRUE)        # use the BGFS optim function because apparently the default 'Nelder Mead' method is unreliable for single parameter models
        if(tempexponenti$convergence!=0){      # try some alternative methods if it doesn't work
              params<-c(mean(-(log(1-.Object@EPDataList[[i]]$Prob)/.Object@EPDataList[[i]]$Effort)),min(.Object@EPDataList[[i]]$Prob))    # this is the initial guess at the exponential parameter a  (not an MLE but not far off!)
              tempexponenti <- optim(params, exponenti ,method="Nelder-Mead", hessian = TRUE)
              cat("NOTE: The optimiser could not converge satisfactorily using the BFGS method. It will now try the Nelder-Mead method", "\n")   
              if(tempexponenti$convergence!=0) cat("ERROR: The optimiser could not converge")
        } 
        AICexponenti <- n*log(tempexponenti$value/n)+4        # calcuate the AIC for this fitted model
        
            # this is the initial guess at the exponential parameter 'a' (not nearly an MLE but not so far off!)
        params<-mean(-(log(1-.Object@EPDataList[[i]]$Prob)/.Object@EPDataList[[i]]$Effort))
        tempexponent <- optim(params, exponent ,method="BFGS", hessian = TRUE)        # use the BGFS optim function because apparently the default 'Nelder Mead' method is unreliable for single parameter models
        if(tempexponent$convergence!=0){      # try some alternative methods if it doesn't work
              params<-mean(-(log(1-.Object@EPDataList[[i]]$Prob)/.Object@EPDataList[[i]]$Effort))    # this is the initial guess at the exponential parameter a  (not an MLE but not far off!)
              tempexponent <- optim(params, exponent ,method="Nelder-Mead", hessian = TRUE)        
              cat("NOTE: The optimiser could not converge satisfactorily using the BFGS method. It will now try the Nelder-Mead method", "\n")   
              if(tempexponent$convergence!=0) cat("ERROR: The optimiser could not converge")
        }
        AICexponent <- n*log(tempexponent$value/n)+2        # calcuate the AIC for this fitted model

        if (AIClogistic < AICexponent & AIClogistic < AICexponenti) {
            .Object@EPFitList[[i]]<-templogistic
            .Object@EPFitFuncs[[i]] <- function(Effort ,Type.no,.Object) (.Object@EPFitList[[Type.no]]$par[1]/(1+exp((.Object@EPFitList[[Type.no]]$par[2]-Effort)/.Object@EPFitList[[Type.no]]$par[3])))
            cat("logistic model", "\n")
        }
        
        else if (AICexponenti<AICexponent){
            .Object@EPFitList[[i]]<-tempexponenti
            .Object@EPFitFuncs[[i]] <- function(Effort ,Type.no,.Object) (1-exp(-.Object@EPFitList[[Type.no]]$par[1]*Effort))*(1-.Object@EPFitList[[Type.no]]$par[2])+.Object@EPFitList[[Type.no]]$par[2]
            cat("exponential model with intercept", "\n")
        }
        
        else{
            .Object@EPFitList[[i]]<-tempexponent
            .Object@EPFitFuncs[[i]] <- function(Effort ,Type.no,.Object) (1-exp(-.Object@EPFitList[[Type.no]]$par[1]*Effort))
            cat("exponential model without intercept", "\n")
        }
        
        if (fitinfo == TRUE){
          print(.Object@EPFitFuncs[[i]])
          print(.Object@EPFitList[[i]])
        }
    }
    else{
        cat("User defined","\n")
        # create a new "list" and populate it with the user defined parameter values - this is so that the functions of EPFitList are identical and can be used by other generic functions (as is the case in BCalc)
        .Object@EPFitList[[i]]<-new("list")
        .Object@EPFitList[[i]]$par[1] <- .Object@controller$sym[i]
        .Object@EPFitList[[i]]$par[2] <- .Object@controller$xmd[i]
        .Object@EPFitList[[i]]$par[3] <- .Object@controller$sc[i]
        if(.Object@controller$EPMod[i] == "logistic"){
            
            .Object@EPFitFuncs[[i]] <- function(Effort,Type.no,.Object) .Object@EPFitList[[Type.no]]$par[1]/(1+exp((.Object@EPFitList[[Type.no]]$par[2]-Effort)/.Object@EPFitList[[Type.no]]$par[3]))
            cat("logistic model","\n")
        }
        else if(.Object@controller$EPMod[i] == "exponentiali"){

            .Object@EPFitFuncs[[i]] <- function(Effort, Type.no,.Object) (1-exp(-.Object@EPFitList[[Type.no]]$par[1]*Effort))*(1-.Object@EPFitList[[Type.no]]$par[2])+.Object@EPFitList[[Type.no]]$par[2]
            cat("exponential model with intercept","\n")
        }
        else if(.Object@controller$EPMod[i] == "exponential"){

            .Object@EPFitFuncs[[i]] <- function(Effort, Type.no,.Object) 1-exp(-.Object@EPFitList[[Type.no]]$par[1]*Effort)
            cat("exponential model without intercept","\n")
        }
        else{

            .Object@EPFitFuncs[[i]] <- function(Effort, Type.no,.Object) {ifelse(Effort>.Object@EPFitList[[Type.no]]$par[2],.Object@EPFitList[[Type.no]]$par[3],.Object@EPFitList[[Type.no]]$par[1]) }
            .Object@Effort[i]<-as.integer(.Object@Effort[i])
            cat("step function: note that the BCalc optimisation of effort does not function correctly with this user defined relationship. ","\n")
            # On trying a very steep logistic model (below) the surface of social benefits becomes distorted because costs and prob increase unevenly with effort over the space of the inflection point
            #{.Object@EPFitList[[Type.no]]$par[1]/(1+exp((.Object@EPFitList[[Type.no]]$par[2]+0.03-Effort)/0.005)) }
        }
        
        if (fitinfo == TRUE){
            print(.Object@EPFitFuncs[[i]])
        }
    }

    ######################
    # Effort-Cost models #
    ######################
    cat("","\n")
    cat("----- Cost of enforcement versus enforcement effort -----","\n")
    if(is.na(.Object@controller$ECMod[[i]])){ # where models are fitted to data
        cat("Model fitted to data", "\n")
       .Object@ECDataList[[i]] <-subset(.Object@ECdata,Type.no==.Object@controller$Type.no[[i]],select=c(Effort,Cost))
       
       n<-length(.Object@ECDataList[[i]]$Cost)

        # The functions for Cost of enforcement versus effort
        nonlineari <- function(params){
          pred <- params[2]+params[1]*.Object@ECDataList[[i]]$Effort^params[3]
          mahalanobis(.Object@ECDataList[[i]]$Cost,pred,diag(length(pred)))
        }
        
        nonlinear <- function(params){
          pred <- params[1]*.Object@ECDataList[[i]]$Effort^params[2] 
          mahalanobis(.Object@ECDataList[[i]]$Cost,pred,diag(length(pred)))
        }
        
        lineari <- function(params){
           pred <- params[2]+params[1]*.Object@ECDataList[[i]]$Effort
           mahalanobis(.Object@ECDataList[[i]]$Cost,pred,diag(length(pred)))
        }
        
        linear <- function(params){
           pred <- params[1]*.Object@ECDataList[[i]]$Effort
           mahalanobis(.Object@ECDataList[[i]]$Cost,pred,diag(length(pred)))
        }
        
        # some initial guesses to get the optimisers going - note that these guesses are not very good on purpose - the optimiser has to search over a reasonable number of iterations to get a good representation of the parameter covariance matrix
        startslope = (max(.Object@ECDataList[[i]]$Cost)-min(.Object@ECDataList[[i]]$Cost))/(max(.Object@ECDataList[[i]]$Effort)-min(.Object@ECDataList[[i]]$Effort))
        startintercept = min(.Object@ECDataList[[i]]$Cost)

        params<-c(startslope,startintercept)
        templineari<- optim(params, lineari ,method="Nelder-Mead",  hessian = TRUE)
        if(templineari$convergence!=0){   # try some alternative methods if it doesn't work
              params<-c(startslope,startintercept)                  
              templineari<- optim(params, lineari ,method="BFGS",  hessian = TRUE)
              cat("NOTE: The optimiser could not converge satisfactorily using the Nelder-Mead method. It will now try the BFGS method", "\n")
              if(templineari$convergence!=0) cat("ERROR: The optimiser could not converge")
        }
        AIClineari <-  n*log(templineari$value/n) + 4

        params<-c(startslope)
        templinear<- optim(params, linear ,method="BFGS",  hessian = TRUE)
        if(templinear$convergence!=0) cat("ERROR: The optimiser could not converge")
        AIClinear <-  n*log(templinear$value/n) + 2


        params<-c(startslope,startintercept,1)
        tempnonlineari<- optim(params, nonlineari ,method="Nelder-Mead",  hessian = TRUE)
        if(tempnonlineari$convergence!=0){      # try some alternative methods if it doesn't work
             params<-c(startslope,startintercept,1)
             tempnonlineari<- optim(params, nonlineari ,method="BFGS",  hessian = TRUE)
             cat("NOTE: The optimiser could not converge satisfactorily using the Nelder-Mead method. It will now try the BFGS method", "\n")
             if(tempnonlineari$convergence!=0) cat("ERROR: The optimiser could not converge")
        }        
        AICnonlineari <- n*log(tempnonlineari$value/n)+6

        params<-c(startslope,1) ###
        tempnonlinear<- optim(params, nonlinear ,method="Nelder-Mead",  hessian = TRUE)
        if(tempnonlinear$convergence!=0){      # try some alternative methods if it doesn't work
             params<-c(startslope,1)    ###
             tempnonlinear<- optim(params, nonlinear ,method="BFGS",  hessian = TRUE)
             cat("NOTE: The optimiser could not converge satisfactorily using the Nelder-Mead method. It will now try the BFGS method", "\n")
             if(tempnonlinear$convergence!=0) cat("ERROR: The optimiser could not converge")
        }
        AICnonlinear <- n*log(tempnonlinear$value/n)+4

        # Pick the 'best' model ----
        AICs<-c(AIClineari,AIClinear,AICnonlineari,AICnonlinear)    # make the AIC's into a list
        selectAIC <- match(min(AICs,na.rm=T),AICs)                  # find the smallest that returns a value

        
        #if(inherits(try(tempnonlineari),"try-error")){     # if a non-linear model won't work then use a linear one

        if(selectAIC == 1){
            .Object@ECFitList[[i]]<-templineari
            .Object@ECFitFuncs[[i]]<-function(Effort,Type.no,.Object)ifelse(Effort<0.001,0,.Object@ECFitList[[Type.no]]$par[1]*Effort+.Object@ECFitList[[Type.no]]$par[2])
            cat("linear model with intercept", "\n")

        }
        else if(selectAIC==2) {
            .Object@ECFitList[[i]]<-templinear
            .Object@ECFitFuncs[[i]]<-function(Effort,Type.no,.Object).Object@ECFitList[[Type.no]]$par[1]*Effort
            cat("linear model without intercept", "\n")
        }
        else if(selectAIC==3) {
            .Object@ECFitList[[i]]<-tempnonlineari
            .Object@ECFitFuncs[[i]]<-function(Effort,Type.no,.Object)ifelse(Effort<0.001,0,.Object@ECFitList[[Type.no]]$par[1]*Effort^.Object@ECFitList[[Type.no]]$par[3]+.Object@ECFitList[[Type.no]]$par[2])
            cat("nonlinear model with intercept", "\n")
        }
        else{
            .Object@ECFitList[[i]]<-tempnonlinear
            .Object@ECFitFuncs[[i]]<-function(Effort,Type.no,.Object).Object@ECFitList[[Type.no]]$par[1]*Effort^.Object@ECFitList[[Type.no]]$par[2]   ###
            cat("nonlinear model without intercept", "\n")
        }

        # now that the model is selected, some warning messages are in order just in case the optimiser failed to run properly for the nonlinear models:
        if(is.na(AICnonlineari))cat("NOTE: The optimiser could not converge satisfactorily on the parameters of the non-linear model with intercept.", "\n")
        if(is.na(AICnonlinear))cat("NOTE: The optimiser could not converge satisfactorily on the parameters of the non-linear model without intercept.", "\n")
        if(is.na(AICnonlinear) & is.na(AICnonlinear))cat("NOTE: The failure of the optimiser to converge satisfactorily on the parameters of either non-linear models has lead it to accept a linear model.", "\n")

        # this is one of the arguments of the new COBECOS object function: do you want to see the detailed information about fitting?
        if(fitinfo==TRUE){
            print(.Object@ECFitFuncs[[i]])
            print(.Object@ECFitList[[i]])    # print summary of fitting
        }
    }

    else{         # where the user defines the relationship
        cat("User defined", "\n")
        # create a new "list" and populate it with the user defined parameter values - this is so that the functions of EPFitList are identical and can be used by other generic functions (as is the case in BCalc)
        .Object@ECFitList[[i]]<-new("list")
        .Object@ECFitList[[i]]$par[1]<-.Object@controller$slope[i]
        .Object@ECFitList[[i]]$par[2]<-.Object@controller$intercept[i]
        .Object@ECFitList[[i]]$par[3]<-.Object@controller$pow[i]


        if(.Object@controller$ECMod[i] == "nonlinear"){
            .Object@ECFitFuncs[[i]]<-function(Effort,Type.no,.Object) .Object@ECFitList[[Type.no]]$par[1]*Effort^.Object@ECFitList[[Type.no]]$par[3]+.Object@ECFitList[[Type.no]]$par[2]
            cat("nonlinear model","\n")
        }

        else if(.Object@controller$ECMod[i] == "linear"){
        
            .Object@ECFitFuncs[[i]] <- function(Effort,Type.no,.Object) .Object@ECFitList[[Type.no]]$par[1]*Effort+.Object@ECFitList[[Type.no]]$par[2]
            cat("linear model","\n")
        }
        else{
        
             .Object@ECFitFuncs[[i]] <- function(Effort, Type.no,.Object) ifelse(Effort>.Object@ECFitList[[Type.no]]$par[2],.Object@ECFitList[[Type.no]]$par[3],.Object@ECFitList[[Type.no]]$par[1])
             cat("step function. Note that the BCalc optimisation of effort may not work correctly with this user defined relationship","\n")
        }
        if(fitinfo==TRUE){
            print(.Object@ECFitFuncs[[i]])
        }


    }
    
    # Now to plot the fitted / user defined relationships and the data.
    if(graphics == TRUE){
    
        eff<-c(1:100)/100
        
        if(is.na(.Object@controller$EPMod[[i]])){        # plots the fitted prob-effort relationship and observed data
            plot(eff,.Object@EPFitFuncs[[i]](eff,i,.Object),xlab = "Enforcement effort",type = "l", col = "red", pch=19,ylab="Prob",xlim=c(0,1))
            points(.Object@EPDataList[[i]],pch=21,xlab = "Enforcement effort")
            title(paste(.Object@Enames[[i]]," (fitted)"))
        }
        else{     # plots the user defined prob-effort relationship
            plot(eff,.Object@EPFitFuncs[[i]](eff,i,.Object),xlab = "Enforcement effort",type = "l", col = "red", pch=19,ylab="Prob",xlim=c(0.01,1))
            title(paste(.Object@Enames[[i]]," (user defined)"))
        }
        
                
        if(is.na(.Object@controller$ECMod[[i]])){    # plots the fitted cost-effort relationship and observed data
            plot(eff,.Object@ECFitFuncs[[i]](eff,i,.Object),xlab = "Enforcement effort",type = "l", col = "red", pch=19,ylab="Cost",xlim=c(0.01,1))
            points(.Object@ECDataList[[i]],xlab = "Enforcement effort",pch=21)
            title(paste(.Object@Enames[[i]]," (fitted)"))
        }
        else{                  #plots the user defined cost-effort relationship
            plot(eff,.Object@ECFitFuncs[[i]](eff,i,.Object),xlab = "Enforcement effort",type = "l", col = "red", pch=19,ylab="Cost",xlim=c(0.01,1))
            title(paste(.Object@Enames[[i]]," (user defined)"))
        }
        
    }

}


cat("","\n")
cat("","\n")

# requirement of S4 class initialisation - ultimately the function returns the new Object
.Object

})

# METHOD 2
#------------------------------------ BENEFIT CALCULATOR / OPTIMISER ----------------------------------------------

BCalc <- function(.Object,Optimize,stoch,graphics){
     
    cat("","\n")
    cat("","\n")
    cat("","\n")
    cat("----------------------------------------------------------------------------","\n")
    cat("Method 2. Calculate benefits","\n")
    cat("","\n")
    flush.console()
    
    tempEffort<-new("numeric")
    tempLower<-new("numeric")
    tempUpper<-new("numeric")
    tempLower[1:.Object@nE]<-0.001
    tempUpper[1:.Object@nE]<-0.999
    
    .Object@Stochastic <- FALSE # set default - resets a stochastic object that is run deterministically
    
    # reset to a non-vector:
    .Object@Harvest<-1
    .Object@SocialB<-1
    .Object@PrivateB<-1
     
    # the private benefit function
    .Object@PrivateOBJFunc <- function(Eff,rp=1){
        
        for(i in 1:.Object@nE){
          .Object@Costs[i] <- .Object@ECFitFuncs[[i]](Eff[i],i,.Object)
          .Object@Probs[i] <- .Object@EPFitFuncs[[i]](Eff[i],i,.Object)
        }

        .Object@ExpFine <-  sum(.Object@Probs*.Object@Fine)
        .Object@TotalCost <-sum(.Object@Costs)
    
        # Fishing response function
        
        if(.Object@FishingRFunc == "default"){
        
            .Object@Harvest[rp] <- ((.Object@Price-.Object@ExpFine)*.Object@Biomass)/(2*.Object@FCost) 
        }
        
        if(.Object@PrivateBFunc == "default"){
            
            # Private benefit function           
            .Object@Price*.Object@Harvest[rp]-.Object@FCost*((.Object@Harvest[rp]*.Object@Harvest[rp])/.Object@Biomass)-(.Object@ExpFine*.Object@Harvest[rp])
                            
        }
    }    
    
    #the social benefit function    
    .Object@SocialOBJFunc <- function(Eff,rp=1){
        
        for(i in 1:.Object@nE){
          .Object@Costs[i] <- .Object@ECFitFuncs[[i]](Eff[i],i,.Object)
          .Object@Probs[i] <- .Object@EPFitFuncs[[i]](Eff[i],i,.Object)
        }
        
        #slot(.Object,"ExpFine") <- sum(.Object@Probs*.Object@Fine)
        .Object@ExpFine <-  sum(.Object@Probs*.Object@Fine)
        .Object@TotalCost <- sum(.Object@Costs)
        
        if(.Object@FishingRFunc == "default"){
            # Fishing response function
            .Object@Harvest[rp] <- ((.Object@Price-.Object@ExpFine)*.Object@Biomass)/(2*.Object@FCost)
        } 
         
        if(.Object@SocialBFunc == "default"){
             
            # Social benefit function           
            (.Object@Price-.Object@ShadowVB)*.Object@Harvest[rp]-.Object@FCost*((.Object@Harvest[rp]*.Object@Harvest[rp])/.Object@Biomass)-.Object@TotalCost

        }
        
    }
                  
    if(Optimize == TRUE){
        
        #tempEffort[1:.Object@nE] <- 0.005   # initial values for optimization
        tempEffort[1:.Object@nE] <- 0.5
        
        .Object@Optimized<-optim(tempEffort, .Object@SocialOBJFunc, rp=1 ,method="L-BFGS-B", lower = tempLower, upper =tempUpper ,control = c(fnscale = -100))    #control = c(fnscale = -1)
        .Object@Effort <- .Object@Optimized$par  # set new effort to optimised values
        tempOptEffort<-.Object@Optimized$par # record the effort levels at the mean parameter values for the cost-effort, prob-effort relationships (they are used below to replace effort that is optimised on the last sample of parameters)
    }
        
    for(i in 1:.Object@nE){
        .Object@Costs[i] <- .Object@ECFitFuncs[[i]](.Object@Effort[i],i,.Object)
        .Object@Probs[i] <- .Object@EPFitFuncs[[i]](.Object@Effort[i],i,.Object)
    }
        
    .Object@ExpFine <-  sum(.Object@Probs*.Object@Fine)
    .Object@TotalCost <-sum(.Object@Costs)
    
    # Fishing response function
    if(.Object@FishingRFunc == "default"){
        
        .Object@Harvest[1] <- ((.Object@Price-.Object@ExpFine)*.Object@Biomass)/(2*.Object@FCost)
    } 
            
    if(graphics == TRUE){     # generic graphics (output for both stochastic and deterministic runs)
        
        quartz()    # new quartz graphics device (plot window)
        par(mfrow = c(.Object@nE,3),mai=c(0.55,0.55,0.35,0.05))  # set up the plot array
        
        
        for(i in 1:.Object@nE){     # loop through enforcement types
            
            Effort<-c(1:100)/100
            if(is.na(.Object@controller$EPMod[[i]])){   # plots the fitted prob-effort relationship and observed data
                plot(Effort,.Object@EPFitFuncs[[i]](Effort,i,.Object),xlab = "Enforcement effort",type = "l", col = "red", pch=19,ylab="Prob",xlim=c(0,1))
                abline(v=.Object@Effort[i], col = "blue")
                points(.Object@EPDataList[[i]],pch=21,xlab = "Enforcement effort")
                title(paste(.Object@Enames[[i]],"  (fitted)"))
            }
            else{                             # plots the user defined prob-effort relationship
                plot(Effort,.Object@EPFitFuncs[[i]](Effort,i,.Object),xlab = "Enforcement effort",type = "l", col = "red", pch=19,ylab="Prob",xlim=c(0,1))
                abline(v=.Object@Effort[i], col = "blue")
                title(paste(.Object@Enames[[i]],"  (user defined)"))
            }
                                          
            if(is.na(.Object@controller$ECMod[[i]])){    # plots the fitted cost-effort relationship and observed data
                plot(Effort,.Object@ECFitFuncs[[i]](Effort,i,.Object),xlab = "Enforcement effort",type = "l", col = "red", pch=19,ylab="Cost",xlim=c(0,1))
                abline(v=.Object@Effort[i], col = "blue")
                points(.Object@ECDataList[[i]],xlab = "Enforcement effort",pch=21,xlim=c(0,100))
                title(paste(.Object@Enames[[i]],"  (fitted)"))
            }
            else{                   # plots the user defined cost-effort relationship
                Effort<-1:100/100
                plot(Effort,.Object@ECFitFuncs[[i]](Effort,i,.Object),xlab = "Enforcement effort",type = "l", col = "red", pch=19,ylab="Cost",xlim=c(0,1))
                abline(v=.Object@Effort[i], col = "blue")
                title(paste(.Object@Enames[[i]],"  (user defined)"))
            }
            
            # set up temporary effort array for social benefit profiles
            npoints<-100    # resolution of the social benefit profile
            Ef<-matrix(.Object@Effort[1:length(.Object@Effort)],npoints,length(.Object@Effort),byrow = TRUE)   # set up the temporary effort array (changes for each enforcement type)
            SocLHP<-1:npoints    # set up the initial vector of social benefits 
            Ef[1:npoints,i]<-1:npoints/npoints   # alter the ith (enforcement type) array of effort
            for(w in 1:npoints){
               SocLHP[w]<-.Object@SocialOBJFunc(Eff=Ef[w,1:.Object@nE],rp=1) # calculate the array of social benefits for each instance of efforts (note only the ith enforcement type has effort that changes)
            }
            plot(Ef[,i],SocLHP,type="l",xlab = "Enforcement effort",ylab="Social benefits") # plot the marginal social benefits
            abline(v=.Object@Effort[i],col="blue")    # plot the current level of effort (user defined or optimized)
            if(Optimize == FALSE) abline(v=Ef[match(max(SocLHP),SocLHP),i],col="red",lty="dotted") # plot the marginal effort where social benefits are maximised if optimization has not been undertaken
            title(paste(.Object@Enames[[i]],": benefit profile"))
            
          
          
        } # End of loop of enforcement types
    }   # End of if graphics = true loop  
     
    
    
    if(stoch != 0){      # Extra stuff where stochasticity is employed 
      
      .Object@Stochastic <- TRUE   
      
      if(stoch ==1 & .Object@ErrCV==0){.Object@Stochastic<-FALSE}    # where CV is zero of Implemenation error and this is the only form of stochasticity, there is effectively no stochasticity
      
      if(stoch !=1){     # just stoch 2 & 3
        cat("    % complete","\n")
        cat("   ")
        PerCentComplete<-1:10*(.Object@StochN/10)
        counter = 1
           
        
        oldEPparams <- new("list")
        oldECparams <- new("list")
        EPtemp<-new("list")
        ECtemp<-new("list")
                   
        for(i in 1:.Object@nE){
           oldEPparams[[i]] <- .Object@EPFitList[[i]]$par
           oldECparams[[i]] <- .Object@ECFitList[[i]]$par
           if(is.na(.Object@controller$EPMod[[i]])){  #multivariate samples are taken at once to improve running speed
               EPtemp[[i]]<-mvrnorm(n=.Object@StochN,mu=oldEPparams[[i]],Sigma=solve(.Object@EPFitList[[i]]$hessian),tol=1e-4) 
           }
           if(is.na(.Object@controller$ECMod[[i]])){
              ECtemp[[i]]<- mvrnorm(n=.Object@StochN,mu=oldECparams[[i]],Sigma=solve(.Object@ECFitList[[i]]$hessian),tol=1e-4)  
           }     
        } 
                       
        for(k in 1:.Object@StochN){
                             
            for(i in 1:.Object@nE){
                if(is.na(.Object@controller$EPMod[[i]])){                                   
                    .Object@EPFitList[[i]]$par<- EPtemp[[i]][k,]
                }
                if(is.na(.Object@controller$ECMod[[i]])){
                    .Object@ECFitList[[i]]$par<-ECtemp[[i]][k,]
                }
            }
            if(Optimize == TRUE){
      
                tempEffort[1:.Object@nE] <- 0.5   # initial values for optimization
                .Object@Optimized<-optim(tempEffort, .Object@SocialOBJFunc,rp=k ,method="L-BFGS-B", lower = tempLower, upper = tempUpper,control = c(fnscale = -100))    #control = c(fnscale = -1)
                .Object@Effort <- .Object@Optimized$par  # set new effort to optimised values
        
            }
                             
            for(i in 1:.Object@nE){
                .Object@Costs[i] <- .Object@ECFitFuncs[[i]](.Object@Effort[i],i,.Object)
                .Object@Probs[i] <- .Object@EPFitFuncs[[i]](.Object@Effort[i],i,.Object)                 
            }
                
            .Object@ExpFine <-  sum(.Object@Probs*.Object@Fine)
            .Object@TotalCost <-sum(.Object@Costs)  
                                  
            if(.Object@FishingRFunc == "default"){
                # Fishing response function
                .Object@Harvest[k] <- ((.Object@Price-.Object@ExpFine)*.Object@Biomass)/(2*.Object@FCost)
            } 
      
            .Object@PrivateB[k] <- .Object@PrivateOBJFunc(.Object@Effort,k)
            .Object@SocialB[k] <- .Object@SocialOBJFunc(.Object@Effort,k)
                
                
            if(k==PerCentComplete[counter]){  # percent complete meter.. updates every 10 per cent
                cat(paste(" ",as.character(counter*10),"...",sep=""))
                counter=counter+1
                flush.console() 
            }
                
        }  # end of loop of stochastic samples 
        cat("\n")   
        for(i in 1:.Object@nE){ # replace the old cost-effort prob-effort relationship parameters (otherwise they would be set to the last sample)
            .Object@EPFitList[[i]]$par<-oldEPparams[[i]]
            .Object@ECFitList[[i]]$par<-oldECparams[[i]]
        }  
        if(Optimize ==TRUE){   # replace the old optimised effort levels (that are a product of a sample of cost-effort, prob-effort parameters) with the optimised effort from the mean parameter values
            .Object@Effort <- tempOptEffort #reset to old optimised effort levels
        }
      
      } 
     
        
      # add implementation error      (stoch = 1 or 3)
      if(.Object@ErrCV>0){     # only apply method if error is set above zero
          if(stoch!=2){  
          
              if(stoch==3){
              
                  for(rp in 1 : .Object@StochN){    # add a rlnorm error to the harvests already found
                         .Object@Harvest[rp] <- .Object@Harvest[rp]* rlnorm(1,meanlog=0,sdlog=.Object@ErrCV)
                         
                  }
                     
              }
              else{    # stoch = 1
                 tempharv<-.Object@Harvest       # temporarily store harvest
                 # creat a rlnorm error around the single estimated value of harvest
                 .Object@Harvest<- rlnorm(.Object@StochN,log(tempharv),sdlog=.Object@ErrCV)
              }
              
              for(rp in 1 : .Object@StochN){    # add a rlnorm error to the harvests already found
                         
                  if(.Object@SocialBFunc == "default"){
                      .Object@SocialB[rp]<- (.Object@Price-.Object@ShadowVB)*.Object@Harvest[rp]-.Object@FCost*((.Object@Harvest[rp]*.Object@Harvest[rp])/.Object@Biomass)-.Object@TotalCost
                  }
                  if(.Object@PrivateBFunc == "default"){
                      # Private benefit function           
                      .Object@PrivateB[rp]<- .Object@Price*.Object@Harvest[rp]-.Object@FCost*((.Object@Harvest[rp]*.Object@Harvest[rp])/.Object@Biomass)-(.Object@ExpFine*.Object@Harvest[rp])
                  }
              }
          
          
          }
       }
    }
    if(.Object@Stochastic==TRUE){    
      # stochastic console output         
      cat("","\n")
      cat("Private profit = ","\n")
      cat(paste("mean = ", format(mean(.Object@PrivateB),digits=4), "   S.D. = ", format(sd(.Object@PrivateB),digits=4)),"\n")
      cat("","\n")
      cat("Social benefit = ", "\n")
      cat(paste("mean = ", format(mean(.Object@SocialB),digits=4), "   S.D. = ", format(sd(.Object@SocialB),digits=4)),"\n")
      cat("","\n")
      cat("Harvest = ", "\n")
      cat(paste("mean = ", format(mean(.Object@Harvest),digits=4), "   S.D. = ", format(sd(.Object@Harvest),digits=4)),"\n")
      cat("","\n")    
      if (Optimize == TRUE){
        cat("Optimized effort levels:", "\n")
      }
      else{
        cat("User defined effort levels:", "\n")
      }
          
      cat(.Object@Effort,"\n")
          
      if (Optimize == TRUE){
         cat("Note that 0.001 to 0.999 is the range of enforcement effort that can be evaluated by the optimizer","\n")
      }
          
      if(graphics == TRUE){ 
            
          quartz()
          par(mfrow = c(3,1),mai=c(0.55,0.55,0.35,0.05))  # set up the plot array
                        
          r <- density(.Object@PrivateB,adjust=1.1,n=800)
          xl <- c(quantile(.Object@PrivateB,0.001),quantile(.Object@PrivateB,0.999))
          plot(r,xlab = "Private profits",main="",xlim=xl,col="red",ylab="Relative Frequency")
         
          r <- density(.Object@SocialB,adjust=1.1,n=800)
          xl <- c(quantile(.Object@SocialB,0.001),quantile(.Object@SocialB,0.999))
          plot(r,xlab = "Social benefits",main="",xlim=xl,col="red",ylab="Relative Frequency")
              
          
          r <- density(.Object@Harvest,adjust=1.1,n=800)
          xl <- c(quantile(.Object@Harvest,0.001),quantile(.Object@Harvest,0.999))
          plot(r,xlab = "Harvest",main="",xlim=xl,col="red",ylab="Relative Frequency")
    
     }
    }
    else{  # deterministic console output
                
        cat("Private profit = ","\n")
        .Object@PrivateB <- .Object@PrivateOBJFunc(.Object@Effort,1)
        cat(format(.Object@PrivateB,digits=4),"\n")
        cat("","\n")
        cat("Social benefit = ", "\n")
        .Object@SocialB<- .Object@SocialOBJFunc(.Object@Effort,1)
        cat(format(.Object@SocialB,digits=4),"\n")
        cat("","\n")
        cat("Harvest level = ", "\n")
        cat(format(.Object@Harvest,digits=4),"\n")
        cat("","\n")    
        if (Optimize == TRUE){
            cat("Optimized effort levels:", "\n")
        }
        else{
            cat("User defined effort levels:", "\n")
        }
            
        cat(.Object@Effort,"\n")
        
        if (Optimize == TRUE){
            cat("\n")
            cat("(Note that 0.001 to 0.999 is the range of enforcement effort that can be evaluated by the optimizer)","\n")
        }
    }
    
    
    .Object   #The function returns the new object. 
}

# METHOD 3
#------------------------------------ COMPARE BENEFITS ----------------------------------------------


BComp<-function(a){
    
    cat("","\n")
    cat("","\n")
    cat("","\n")
    cat("----------------------------------------------------------------------------","\n")
    cat("Method 3. Compare COBECOS objects","\n")
    cat("","\n")    
    flush.console()
        
    #string<-new("character")
    #temp<-new("character")   
    Bcol<-colours()[c(28,552,50,24,148,367,52,464,652,636,401,139,573,609,30,90,204,130,137,153,552,456,139,451,525,137,628,588,460,573,631)]   # define the colours for the lines of graphs
    
    string <- deparse(substitute(a[[]]),width.cutoff = 500,.deparseOpts(control=c("warnIncomplete")))   # take the name of the list a
    n <- nchar(string)                       # take its length
    string<-substr(string,3,n-6)            # remove the "c(" from the start and the "([[i]]" from the end
    temp <-strsplit(string, split=", ")      # turn the string into a list of items
    nams <- as.character(temp[[1]])          # turn the list into a character vector of names
    nEnf<-new("numeric")
    maxe<-0  # upper bounds for the y-axis of the effort bar plots   
    
    # the header for the console output
    dns<-new("list")
    dns[[1]]<-c(" ",nams)
    print(dns)
    dns[[2]]<-c("Stochastic","Harvest","","Social ","benefits","Private ","profits","Effort","","","","","","","","","","","","","","","","","")
    output<-matrix("",nrow = length(a)+1, ncol = 25, dimnames=dns)

    output[1,2]<-"mean"
    output[1,3]<-"sd"
    output[1,4]<-"mean"
    output[1,5]<-"sd"
    output[1,6]<-"mean"
    output[1,7]<-"sd"
    for(i in 1:a[[1]]@nE){
      output[1,7+i]<-as.character(a[[1]]@Enames[[i]])
    }
      
    # initial bounds for the x-axis of density plots
    harvxlim <-c(10000000000,-10000000000)   
    socxlim <-c(10000000000,-10000000000)
    prixlim <-c(10000000000,-10000000000)

    # initial upper bound for the y-axis of density plots
    harvylim<--1000000
    socylim<--1000000
    priylim<--1000000
    temp <- 0
    
    for(i in 1:length(a)){
               
        nEnf[i]<-a[[i]]@nE
        output[i+1,1]<-a[[i]]@Stochastic  
        
        # limits for the effort plots    
        if(max(a[[i]]@Effort) > maxe)  {maxe<- max(a[[i]]@Effort)}   
        
        if(a[[i]]@Stochastic ==TRUE){
            
            # populate the output table
            
            output[i+1,2]<-format(mean(a[[i]]@Harvest),digits=4)
            output[i+1,3]<-format(sd(a[[i]]@Harvest),digits=4)
            output[i+1,4]<-format(mean(a[[i]]@SocialB),digits=4)
            output[i+1,5]<-format(sd(a[[i]]@SocialB),digits=4)
            output[i+1,6]<-format(mean(a[[i]]@PrivateB),digits=4)
            output[i+1,7]<-format(sd(a[[i]]@PrivateB),digits=4)
            for(k in 1:a[[1]]@nE){
                output[i+1,7+k]<-format(a[[i]]@Effort[k],digits=4)
            }
            
            # limits for the effort plots    
            #if(max(a[[i]]@Effort) > maxe)  {maxe<- max(a[[i]]@Effort)}
                    
            # limits for the density plots (if at least one object is stochastic
            if(max(density(a[[i]]@Harvest)$y)>harvylim)   {harvylim<-max(density(a[[i]]@Harvest)$y)}
            if(max(density(a[[i]]@SocialB)$y)>socylim)    {socylim<-max(density(a[[i]]@SocialB)$y)}
            if(max(density(a[[i]]@PrivateB)$y)>priylim)   {priylim<-max(density(a[[i]]@PrivateB)$y)}
            
            
            # set upper and lower bounds for the harvest density plot
            if(quantile(a[[i]]@Harvest,0.001)< harvxlim[1]){
                harvxlim[1]<-quantile(a[[i]]@Harvest,0.001)
            }
            if(quantile(a[[i]]@Harvest,0.999)> harvxlim[2]){
                harvxlim[2]<-quantile(a[[i]]@Harvest,0.999)
            }
            
            # set upper and lower bounds for the social benefits density plot
            if(quantile(a[[i]]@SocialB,0.001)< socxlim[1]){
                socxlim[1]<-quantile(a[[i]]@SocialB,0.001)
            }
            if(quantile(a[[i]]@SocialB,0.999)> socxlim[2]){
                socxlim[2]<-quantile(a[[i]]@SocialB,0.999)
            }
            
            # set upper and lower bounds for the private profits density plot
            if(quantile(a[[i]]@PrivateB,0.001)< prixlim[1]){
                prixlim[1]<-quantile(a[[i]]@PrivateB,0.001)
            }
            if(quantile(a[[i]]@PrivateB,0.999)> prixlim[2]){
                prixlim[2]<-quantile(a[[i]]@PrivateB,0.999)
            }
            
            if (a[[i]]@Stochastic == TRUE){
                if(temp == 0){
                    temp<-i      # the first object in the list to be stochastic
                }
            }
        }
        else{
        
        
            # populate table
            output[i+1,2]<-format(a[[i]]@Harvest,digits=4)
            output[i+1,4]<-format(a[[i]]@SocialB,digits=4)
            output[i+1,6]<-format(a[[i]]@PrivateB,digits=4)
            for(k in 1:a[[1]]@nE){
                output[i+1,7+k]<-format(a[[i]]@Effort[k],digits=4)
            }
        
            # set upper and lower bounds for the harvest density plot
            if(a[[i]]@Harvest< harvxlim[1]){
                harvxlim[1]<-a[[i]]@Harvest
            }
            if(a[[i]]@Harvest> harvxlim[2]){
                harvxlim[2]<-a[[i]]@Harvest
            }
            
            # set upper and lower bounds for the social benefits density plot
            if(a[[i]]@SocialB< socxlim[1]){
                socxlim[1]<-a[[i]]@SocialB
            }
            if(a[[i]]@SocialB> socxlim[2]){
                socxlim[2]<-a[[i]]@SocialB
            }
            
            # set upper and lower bounds for the private profits density plot
            if(a[[i]]@PrivateB< prixlim[1]){
                prixlim[1]<-a[[i]]@PrivateB
            }
            if(a[[i]]@PrivateB> prixlim[2]){
                prixlim[2]<-a[[i]]@PrivateB
            }
        
        
        }
        
    }
    print(output, quote=FALSE, zero.print = "") #print output table
    
    eff <-array(1,c(nEnf[1],length(a)))
    
    if(all(nEnf == nEnf[1])){    # only does it if all elements have the same number of enforcement types - a condition of comparability
        
        effort<-new("list")
        
        quartz()    # new quartz graphics device (plot window)
        par(mfrow = c(nEnf[1],1),mai=c(0.55,0.55,0.35,0.05))  # set up the plot array  
        
        for(k in 1:nEnf[1]){
            for(i in 1:length(a)){
                eff[k,i]<-a[[i]]@Effort[k]   # takes the effort levels for each ith object and kth enforcement type
            }    
                  
            if(k==nEnf[1]){
                par(mar=c(3,4,2,0.1))
                barplot(eff[k,],main = as.character(a[[1]]@Enames[[k]][1]),ylim=c(0,maxe),names.arg=nams,ylab="Effort",col=Bcol)   
                         
            }
            else{
                par(mar=c(0.7,4,2,0.1))
                barplot(eff[k,],main = as.character(a[[1]]@Enames[[k]][1]),ylim=c(0,maxe),xlab = "",ylab="Effort",col=Bcol)   
                axis(1, labels = F)
            }
          
        }
                
        quartz()    # new quartz graphics device (plot window)
        par(mfrow = c(3,1),mai=c(0.55,0.55,0.35,0.05))  # set up the plot array 
             
        if(all(temp == 0)){    # only does loop (a series of bar charts) for a full array of deterministic COBECOS objects

            harv <-new("numeric")
            SB <-new("numeric")
            PP <-new("numeric")
            
            for(i in 1:length(a)){
                harv[i]<-a[[i]]@Harvest   # takes the effort levels for each ith object and kth enforcement type
                SB[i] <-a[[i]]@SocialB
                PP[i] <-a[[i]]@PrivateB
            }
            par(mar=c(3,4,2,0.2))
            barplot(harv,main = "Harvest",xlab = "",ylab="",col=Bcol)   
            axis(1, labels = F)
            
            par(mar=c(3,4,2,0.2))
            barplot(SB,main = "Social benefits",xlab = "",ylab="",col=Bcol)   
            axis(1, labels = F)
            
            par(mar=c(3,4,2,0.2))
            barplot(PP,main = "Private profits",names.arg=nams,ylab="",col=Bcol)   
                        
        }
        else{               # if there is at least one stochastic object (a density plot with vertical lines for deterministic objects)
                   
            # harvest graph
            par(mar=c(4,4,1,0.2))
            r <- density(a[[temp]]@Harvest,adjust=1.5)
            plot(r,xlab = "Harvest",main="",xlim=harvxlim,ylim=c(0,harvylim),col=Bcol[temp],ylab="Relative Frequency")  # plot the first object to be stochastic
        
            for(i in 1:length(a)){  # loop through objects
            
                if(i!=temp){       # if the object is not the first to be stochastic which has already been plotted (and thus requires a density plot)
                    if(a[[i]]@Stochastic == TRUE){
                        r <- density(a[[i]]@Harvest,adjust=1.3,n=800)
                        lines(r,col=Bcol[i])
                    }
                    else{
                        abline(v=a[[i]]@Harvest,  col = Bcol[i])
                    }
                }
            }
            legend("topright", legend = nams,text.col=Bcol)

            # social benefit graph
            par(mar=c(4,4,1,0.2))
            r <- density(a[[temp]]@SocialB,adjust=1.5)
            plot(r,xlab = "Social benefits",xlim = socxlim,ylim=c(0,socylim),main="",col=Bcol[temp],ylab="Relative Frequency")  # plot the first object to be stochastic
        
            for(i in 1:length(a)){  # loop through objects
                if(i!=temp){       # if the object is not the first to be stochastic which has already been plotted (and thus require a density plot)
                    if(a[[i]]@Stochastic == TRUE){
                        r <- density(a[[i]]@SocialB,adjust=1.3,n=800)
                        lines(r,col=Bcol[i])
                    }
                    else{
                        abline(v=a[[i]]@SocialB,  b = 10000000, col = Bcol[i])
                    }
                }
            }            
            
            # Private profits
            par(mar=c(4,4,1,0.2))
            r <- density(a[[temp]]@PrivateB,adjust=1.5)
            plot(r,xlab = "Private profits",xlim = prixlim,ylim=c(0,priylim), main="",col=Bcol[temp],ylab="Relative Frequency")  # plot the first object to be stochastic
        
            for(i in 1:length(a)){  # loop through objects
                if(i!=temp){       # if the object is not the first to be stochastic which has already been plotted (and thus require a density plot)
                    if(a[[i]]@Stochastic == TRUE){
                        r <- density(a[[i]]@PrivateB,adjust=1.3,n=800)
                        lines(r,col=Bcol[i])
                    }
                    else{
                        abline(v=a[[i]]@PrivateB,  b = 10000000, col = Bcol[i])
                    }
                }
            }
        
        }

     } # end of check for equal enforcment numbers


     
}

cat("","\n")
cat("Source code loaded successfully", "\n")
cat("","\n")
cat("","\n")