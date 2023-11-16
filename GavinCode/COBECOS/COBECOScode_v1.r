

cat("\n ---------------------------------------------------------------------------- \n")
cat(" COBECOS (COsts and BEnefits of COntrol Strategies) \n")
cat(" Alpha v.1 \n")
cat(" October 2008 \n")
cat("\n")
cat(" T.Carruthers & C.Edwards \n")
cat(" Imperial College \n")
cat("\n ---------------------------------------------------------------------------- \n")


# METHOD 1
# -----------------------------------  INITIALISATION  ---------------------------------------------------------------------------
# A number of methods are defined which allow the user to investigate cost-benefit
# create a class of control object called "COBECOS"

setClass("COBECOS",
representation( 
r = "numeric", #GM
K = "numeric", #GM
stockGrowth = "function", #GM
Ename = "character", 
EPData = "data.frame",
ECData = "data.frame",
EPFitPar = "numeric", 
ECFitPar = "numeric", 
EPFitFunc = "function", 
ECFitFunc = "function", 
control = "data.frame",
Effort="numeric", 
StochEffort = "numeric",
StochProb = "numeric",
StochCost = "numeric", 
Fine = "numeric", 
Cost = "numeric",  
Prob = "numeric",
UnitTax = "numeric", 
TAC = "numeric", 
Biomass = "numeric",
FCost = "numeric",
Harvest = "numeric",
StochExpFine = "numeric",
StochHarvest = "numeric",
Price = "numeric",
ShadowVB= "numeric",
ExpFine = "numeric", 
SocialBFunc = "function", 
PrivateBFunc = "function", 
FishingRFunc = "function",
HarvestErr="numeric",
Effort_Min = "numeric", 
Effort_Max = "numeric", 
PrivateB = "numeric",
StochPrivateB = "numeric", 
SocialB = "numeric",
StochSocialB = "numeric", 
SocialOBJFunc = "function",  
Stochastic = "logical", 
ImpError = "logical", 
EstError = "logical", 
Optimized = "list", 
StochN = "numeric", 
ImpErrorSE = "numeric", 
EPEstErrorSE = "numeric",
ECEstErrorSE = "numeric",
EPErr = "numeric",
ECErr = "numeric"), 
prototype(
StochN = 1000,
ImpErrorSE = 0.01, 
EPEstErrorSE = 0.01, 
ECEstErrorSE = 0.01, 
Stochastic = FALSE, 
ImpError = FALSE, 
EstError = FALSE, 
PrivateBFunc = NULL, 
SocialBFunc = NULL))

# Create an initialisation method of the COBECOS control object

setMethod("initialize", "COBECOS", function(.Object, path, fitinfo, graphics) {

    # ------------------------------------------------------------------------ #
    # ------- PRELIMINARIES -------------------------------------------------- #
    # ------------------------------------------------------------------------ #
    # Read in the observed data for which functions and parameter are to be fit
    .Object@control <- read.csv(file = paste(path, "Control.csv",sep = ""),head=TRUE,sep=",")
    
    readEPdata <- FALSE
    if(is.na(.Object@control$EPMod)) {
      readEPdata <- TRUE
    } else {
      if(.Object@control$EPMod=="logistic") {
        if(any(is.na(c(.Object@control$EPpar1,.Object@control$EPpar2,.Object@control$EPpar3))))
          readEPdata <- TRUE
      }
      if(.Object@control$EPMod=="exponential") {
        if(any(is.na(c(.Object@control$EPpar1,.Object@control$EPpar2))))
          readEPdata <- TRUE
      }
    }
    if(readEPdata) {
      .Object@EPData <-read.csv(file = paste(path, "ProbData.csv",sep = ""),head=TRUE,sep=",")
      .Object@EPData<-.Object@EPData[order(.Object@EPData$Effort),]
    } else {
      .Object@EPData <- data.frame(NA)
    }
    
    readECdata <- FALSE
    if(is.na(.Object@control$ECMod)) {
      readECdata <- TRUE
    } else {
      if(.Object@control$ECMod=="nonlinear") {
        if(any(is.na(c(.Object@control$ECpar1,.Object@control$ECpar2,.Object@control$ECpar3))))
          readECdata <- TRUE
      }
      if(.Object@control$ECMod=="linear") {
        if(any(is.na(c(.Object@control$ECpar1,.Object@control$ECpar2))))
          readECdata <- TRUE
      }
    }
    if(readECdata) {
      .Object@ECData <- read.csv(file = paste(path, "CostData.csv",sep = ""),head=TRUE,sep=",")
      .Object@ECData<-.Object@ECData[order(.Object@ECData$Effort),]
    } else {
      .Object@ECData <- data.frame(NA)
    }
    
    .Object@Ename <- as.character(.Object@control$Type)
        
    # load function parameters
    .Object@EPFitPar<-as.numeric(rep(NA,3))  
    if(!is.na(.Object@control$EPpar1)) {.Object@EPFitPar[1] <- .Object@control$EPpar1 }
    if(!is.na(.Object@control$EPpar2)) {.Object@EPFitPar[2] <- .Object@control$EPpar2 }
    if(!is.na(.Object@control$EPpar3)) {.Object@EPFitPar[3] <- .Object@control$EPpar3 }
    
    .Object@ECFitPar<-as.numeric(rep(NA,3))    
    if(!is.na(.Object@control$ECpar1)) {.Object@ECFitPar[1] <- .Object@control$ECpar1 }
    if(!is.na(.Object@control$ECpar2)) {.Object@ECFitPar[2] <- .Object@control$ECpar2 }
    if(!is.na(.Object@control$ECpar3)) {.Object@ECFitPar[3] <- .Object@control$ECpar3 }
    
    # load Effort Min and Effort Max
    .Object@Effort <- as.numeric(NA)
    if(is.na(.Object@control$Effort_Min)) stop("Assign Effort_Min value in control file")
    .Object@Effort_Min <- .Object@control$Effort_Min
    if(is.na(.Object@control$Effort_Max)) stop("Assign Effort_Max value in control file")
    .Object@Effort_Max <- .Object@control$Effort_Max
    
    cat("\n Initialisation of COBECOS object for enforcement type:",.Object@Ename,"\n \n")
    flush.console()

    # ------------------------------------------------------------------------ #
    # ------- DEFINE EFFORT-PROBABILITY MODELS ------------------------------- #
    # ------------------------------------------------------------------------ # 
    logistic <- function(params){
        if(!is.na(.Object@EPFitPar[1])) { params[1]<-.Object@EPFitPar[1] }
        if(!is.na(.Object@EPFitPar[2])) { params[2]<-.Object@EPFitPar[2] }
        if(!is.na(.Object@EPFitPar[3])) { params[3]<-.Object@EPFitPar[3] }
        pred <- params[1]/(1+exp((params[2]-.Object@EPData$Effort)/params[3]))
        mahalanobis(.Object@EPData$Prob,pred,diag(length(pred)))
    }
        
    exponenti <- function(params){
        if(!is.na(.Object@EPFitPar[1])) { params[1]<-.Object@EPFitPar[1] }
        if(!is.na(.Object@EPFitPar[2])) { params[2]<-.Object@EPFitPar[2] }
        pred <- (1-exp(-params[1]*.Object@EPData$Effort))*(1-params[2])+params[2]
        mahalanobis(.Object@EPData$Prob,pred,diag(length(pred)))
        }
        
    exponent <- function(params){
        if(!is.na(.Object@EPFitPar[1])) { params[1]<-.Object@EPFitPar[1] }
        pred <- 1-exp(-params[1]*.Object@EPData$Effort)
        mahalanobis(.Object@EPData$Prob,pred,diag(length(pred)))
        }
        
    # ------------------------------------------------------------------------ #
    # ------- DEFINE EFFORT-COST MODELS -------------------------------------- #
    # ------------------------------------------------------------------------ # 
    nonlineari <- function(params){
        if(!is.na(.Object@ECFitPar[1])) { params[1]<-.Object@ECFitPar[1] }
        if(!is.na(.Object@ECFitPar[2])) { params[2]<-.Object@ECFitPar[2] }
        if(!is.na(.Object@ECFitPar[3])) { params[3]<-.Object@ECFitPar[3] }
        pred <- params[2]+params[1]*.Object@ECData$Effort^params[3]
        mahalanobis(.Object@ECData$Cost,pred,diag(length(pred)))
    }
        
    nonlinear <- function(params){
        if(!is.na(.Object@ECFitPar[1])) { params[1]<-.Object@ECFitPar[1] }
        if(!is.na(.Object@ECFitPar[3])) { params[2]<-.Object@ECFitPar[3] }
        pred <- params[1]*.Object@ECData$Effort^params[2] 
        mahalanobis(.Object@ECData$Cost,pred,diag(length(pred)))
    }
        
    lineari <- function(params){
        if(!is.na(.Object@ECFitPar[1])) { params[1]<-.Object@ECFitPar[1] }
        if(!is.na(.Object@ECFitPar[2])) { params[2]<-.Object@ECFitPar[2] }
        pred <- params[2]+params[1]*.Object@ECData$Effort
        mahalanobis(.Object@ECData$Cost,pred,diag(length(pred)))
    }
        
    linear <- function(params){
        if(!is.na(.Object@ECFitPar[1])) { params[1]<-.Object@ECFitPar[1] }
        pred <- params[1]*.Object@ECData$Effort
        mahalanobis(.Object@ECData$Cost,pred,diag(length(pred)))
    }
    
    cat(" ------------------ \n")    
    # ------------------------------------------------------------------------ #
    # ------- SELECT EP MODEL ACCORDING TO DATA ------------------------------ #
    # ------------------------------------------------------------------------ # 
    if(is.na(.Object@control$EPMod)){  # Fit effort - prob models - only do this if the user hasn't specified a model for this effort probability model ie they'v stuck an "NA" in the EPMod column of the control.csv excel file
        cat("\n AIC selected EP model: ")      
        
        # stop if there is no data
        if(all(is.na(.Object@EPData))) stop(paste("\n No prob. data for enforcement type: ",.Object@Ename,"\n",sep=""))
        
        # re-set all parameters
        .Object@EPFitPar<-as.numeric(rep(NA,3))  
        
        n<-length(.Object@EPData$Prob)
        
        #######################
        # FIT COMPLETE MODELS #
        #######################
        params<-c(0.95,0.5,0.05)
        templogistic <- optim(params, logistic ,method="Nelder-Mead")
        if(templogistic$convergence!=0){      # try some alternative methods if it doesn't work
             params<-c(0.95,0.5,0.05)
             templogistic <- optim(params, logistic ,method="BFGS")
             cat("NOTE: The optimiser could not converge satisfactorily using the Nelder-Mead method. It will now try the BFGS method", "\n")   
             if(templogistic$convergence!=0) stop("The optimiser could not converge")
        }
        AIClogistic <- n*log(templogistic$value/n)+6
        
        params<-c(mean(-(log(1-.Object@EPData$Prob)/.Object@EPData$Effort)),min(.Object@EPData$Prob)) # this is the initial guess at the exponential parameter 'a' (not nearly an MLE but not so far off!)
        tempexponenti <- optim(params, exponenti ,method="BFGS")        # use the BGFS optim function because apparently the default 'Nelder Mead' method is unreliable for single parameter models
        if(tempexponenti$convergence!=0){      # try some alternative methods if it doesn't work
              params<-c(mean(-(log(1-.Object@EPData$Prob)/.Object@EPData$Effort)),min(.Object@EPData$Prob))    # this is the initial guess at the exponential parameter a  (not an MLE but not far off!)
              tempexponenti <- optim(params, exponenti ,method="Nelder-Mead")
              cat("NOTE: The optimiser could not converge satisfactorily using the BFGS method. It will now try the Nelder-Mead method", "\n")   
              if(tempexponenti$convergence!=0) stop("The optimiser could not converge")
        } 
        AICexponenti <- n*log(tempexponenti$value/n)+4        # calcuate the AIC for this fitted model
              
        params<-mean(-(log(1-.Object@EPData$Prob)/.Object@EPData$Effort))      # this is the initial guess at the exponential parameter 'a' (not nearly an MLE but not so far off!)
        tempexponent <- optim(params, exponent ,method="BFGS")        # use the BGFS optim function because apparently the default 'Nelder Mead' method is unreliable for single parameter models
        if(tempexponent$convergence!=0){      # try some alternative methods if it doesn't work
              params<-mean(-(log(1-.Object@EPData$Prob)/.Object@EPData$Effort))    # this is the initial guess at the exponential parameter a  (not an MLE but not far off!)
              tempexponent <- optim(params, exponent ,method="Nelder-Mead")        
              cat("NOTE: The optimiser could not converge satisfactorily using the BFGS method. It will now try the Nelder-Mead method", "\n")   
              if(tempexponent$convergence!=0) stop("The optimiser could not converge")
        }
        AICexponent <- n*log(tempexponent$value/n)+2        # calcuate the AIC for this fitted model
        
        ##############################
        # SELECT MODEL USING THE AIC #
        ##############################
        AICs<-c(AIClogistic,AICexponenti,AICexponent)    # make the AIC's into a list
        selectAIC <- match(min(AICs,na.rm=T),AICs) 
        
        if (selectAIC == 1) {
            .Object@EPFitPar<-templogistic$par
            names(.Object@EPFitPar)<-c("asymptote","location","slope")
            .Object@EPFitFunc <- function(Effort,.Object) (.Object@EPFitPar[1]/(1+exp((.Object@EPFitPar[2]-Effort)/.Object@EPFitPar[3])))
            #names(.Object@EPFitFunc)<-"logistic"
            cat("logistic", "\n")
        }
        
        else if (selectAIC == 2){
            .Object@EPFitPar<-tempexponenti$par
            names(.Object@EPFitPar)<-c("par","intercept")
            .Object@EPFitFunc <- function(Effort,.Object) (1-exp(-.Object@EPFitPar[1]*Effort))*(1-.Object@EPFitPar[2])+.Object@EPFitPar[2]
            #names(.Object@EPFitFunc)<-"exponentiali"
            cat("exponential with intercept", "\n")
        }
        
        else{
            .Object@EPFitPar<-tempexponent$par
            names(.Object@EPFitPar)<-c("par")
            .Object@EPFitFunc <- function(Effort,.Object) (1-exp(-.Object@EPFitPar[1]*Effort))
            #names(.Object@EPFitFunc)<-"exponential"
            cat("exponential without intercept", "\n")
        }
        
        if (fitinfo == TRUE){
          cat("\n Parameter values: \n")
          print(.Object@EPFitPar)
        }
        cat(" ------------------ \n") 
    }
    
    # ------------------------------------------------------------------------ #
    # ------- USER DEFINED EP MODEL WITH SOME OR ----------------------------- #
    # ------- ALL PARAMETERS ESTIMATED FROM DATA ----------------------------- #
    # ------------------------------------------------------------------------ #
    if(!is.na(.Object@control$EPMod)){
        cat("\n User defined EP model: ")
        
        if(.Object@control$EPMod == "logistic"){
        
          cat("logistic \n")
          
          if (any(is.na(.Object@EPFitPar[1:3]))) {   #if any parameters need estimating
            if(!all(is.na(.Object@EPData))) {     #if there is some data
              params<-c(0.95,0.5,0.05)
              templogistic <- optim(params, logistic ,method="BFGS")
              if(templogistic$convergence==0) {
                if(is.na(.Object@EPFitPar[1])) {.Object@EPFitPar[1] <- templogistic$par[1] }
                if(is.na(.Object@EPFitPar[2])) {.Object@EPFitPar[2] <- templogistic$par[2] }
                if(is.na(.Object@EPFitPar[3])) {.Object@EPFitPar[3] <- templogistic$par[3] } 
              } else stop("The optimiser did not converge when fitting logistic function to prob. data \n")
            } else stop("No prob. data to fit logistic function \n")
          }
          names(.Object@EPFitPar)<-c("asymptote","location","slope")  
          .Object@EPFitFunc <- function(Effort,.Object) .Object@EPFitPar[1]/(1+exp((.Object@EPFitPar[2]-Effort)/.Object@EPFitPar[3]))
          #names(.Object@EPFitFunc)<-"logistic"
        }
        else if(.Object@control$EPMod == "exponential"){
           
           cat("exponential \n")
          
          .Object@EPFitPar[3]<-NA
          if (any(is.na(.Object@EPFitPar[1:2]))) {   #if any parameters need estimating
            if(!all(is.na(.Object@EPData))) {     #if there is some data
              params<-c(mean(-(log(1-.Object@EPData$Prob)/.Object@EPData$Effort)),min(.Object@EPData$Prob))
              tempexponent <- optim(params, exponenti ,method="BFGS")
              if(tempexponent$convergence==0) {            #if optimisation was successful
                if(is.na(.Object@EPFitPar[1])) {.Object@EPFitPar[1] <- tempexponent$par[1] }
                if(is.na(.Object@EPFitPar[2])) {.Object@EPFitPar[2] <- tempexponent$par[2] } 
              } else stop("The optimiser did not converge when fitting exponential (with intercept) function to prob. data \n")
            } else stop("No prob. data to fit exponential (with intercept) function \n")
          }
          names(.Object@EPFitPar)<-c("par","intercept")
          .Object@EPFitFunc <- function(Effort,.Object) (1-exp(-.Object@EPFitPar[1]*Effort))*(1-.Object@EPFitPar[2])+.Object@EPFitPar[2]          
          #names(.Object@EPFitFunc)<-"exponential"
        }
        else if(.Object@control$EPMod == "step"){

            cat("step function \n Note that the BCalc optimisation of effort does not function correctly with this relationship. \n")
            
            if (any(is.na(.Object@EPFitPar[1:3]))) stop("Step function parameters not defined")
            .Object@EPFitFunc <- function(Effort,.Object) {ifelse(Effort>.Object@EPFitPar[2],.Object@EPFitPar[3],.Object@EPFitPar[1]) }            
        }
        else stop("\n User defined function not recognised")
        
        if (fitinfo == TRUE){
          cat("\n Parameter values: \n ")
          print(.Object@EPFitPar)
        }
        cat(" ------------------ \n")
    }
    
    # ------------------------------------------------------------------------ #
    # ------- SELECT EC MODEL ACCORDING TO DATA ------------------------------ #
    # ------------------------------------------------------------------------ #
    if(is.na(.Object@control$ECMod)){
      cat("\n AIC selected EC model: ")
        
      # stop if there is no data
      if(all(is.na(.Object@ECData))) stop("No cost data \n")
        
      # re-set all parameters
      .Object@ECFitPar<-as.numeric(rep(NA,3))  
        
      n<-length(.Object@ECData$Cost)
 
      #######################
      # FIT COMPLETE MODELS #
      ####################### 
      # some initial guesses to get the optimisers going - note that these guesses are not very good on purpose - the optimiser has to search over a reasonable number of iterations to get a good representation of the parameter covariance matrix
      startslope = (max(.Object@ECData$Cost)-min(.Object@ECData$Cost))/(max(.Object@ECData$Effort)-min(.Object@ECData$Effort))
      startintercept = min(.Object@ECData$Cost)

      params<-c(startslope,startintercept)
      templineari<- optim(params, lineari ,method="Nelder-Mead",  hessian = TRUE)
      if(templineari$convergence!=0){   # try some alternative methods if it doesn't work
          params<-c(startslope,startintercept)                  
          templineari<- optim(params, lineari ,method="BFGS",  hessian = TRUE)
          cat("NOTE: The optimiser could not converge satisfactorily using the Nelder-Mead method. It will now try the BFGS method", "\n")
          if(templineari$convergence!=0) stop("The optimiser could not converge")
      }
      AIClineari <-  n*log(templineari$value/n) + 4

      params<-c(startslope)
      templinear<- optim(params, linear ,method="BFGS",  hessian = TRUE)
      if(templinear$convergence!=0) stop("The optimiser could not converge")
      AIClinear <-  n*log(templinear$value/n) + 2


      params<-c(startslope,startintercept,1)
      tempnonlineari<- optim(params, nonlineari ,method="Nelder-Mead",  hessian = TRUE)
      if(tempnonlineari$convergence!=0){      # try some alternative methods if it doesn't work
          params<-c(startslope,startintercept,1)
          tempnonlineari<- optim(params, nonlineari ,method="BFGS",  hessian = TRUE)
          cat("NOTE: The optimiser could not converge satisfactorily using the Nelder-Mead method. It will now try the BFGS method", "\n")
          if(tempnonlineari$convergence!=0) stop("The optimiser could not converge")
      }        
      AICnonlineari <- n*log(tempnonlineari$value/n)+6

      params<-c(startslope,1)
      tempnonlinear<- optim(params, nonlinear ,method="Nelder-Mead",  hessian = TRUE)
      if(tempnonlinear$convergence!=0){      # try some alternative methods if it doesn't work
          params<-c(startslope,1)   
          tempnonlinear<- optim(params, nonlinear ,method="BFGS",  hessian = TRUE)
          cat("NOTE: The optimiser could not converge satisfactorily using the Nelder-Mead method. It will now try the BFGS method", "\n")
          if(tempnonlinear$convergence!=0) stop("The optimiser could not converge")
      }
      AICnonlinear <- n*log(tempnonlinear$value/n)+4

      ##############################
      # SELECT MODEL USING THE AIC #
      ##############################
      AICs<-c(AIClineari,AIClinear,AICnonlineari,AICnonlinear)    # make the AIC's into a list
      selectAIC <- match(min(AICs,na.rm=T),AICs)                  # find the smallest that returns a value

      if(selectAIC == 1){
          .Object@ECFitPar<-templineari$par
          names(.Object@ECFitPar)<-c("slope","intercept") 
          .Object@ECFitFunc<-function(Effort,.Object) .Object@ECFitPar[1]*Effort+.Object@ECFitPar[2]
          #names(.Object@ECFitFunc)<-"lineari"
          cat("linear with intercept", "\n")
      }
      else if(selectAIC==2) {
          .Object@ECFitPar<-templinear$par
          names(.Object@ECFitPar)<-c("slope") 
          .Object@ECFitFunc<-function(Effort,.Object) .Object@ECFitPar[1]*Effort
          #names(.Object@ECFitFunc)<-"linear"
          cat("linear without intercept", "\n")
      }
      else if(selectAIC==3) {
          .Object@ECFitPar<-tempnonlineari$par
          names(.Object@ECFitPar)<-c("slope","intercept","power") 
          .Object@ECFitFunc<-function(Effort,.Object) .Object@ECFitPar[1]*Effort^.Object@ECFitPar[3]+.Object@ECFitPar[2]
          #names(.Object@ECFitFunc)<-"nonlineari"
          cat("nonlinear with intercept", "\n")
      }
      else{
          .Object@ECFitPar[1]<-tempnonlinear$par[1]
          .Object@ECFitPar[3]<-tempnonlinear$par[2]
          names(.Object@ECFitPar)<-c("slope","","power") 
          .Object@ECFitFunc<-function(Effort,.Object) .Object@ECFitPar[1]*Effort^.Object@ECFitPar[3]  
          #names(.Object@ECFitFunc)<-"nonlinear"
          cat("nonlinear without intercept", "\n")
      }

      # now that the model is selected, some warning messages are in order just in case the optimiser failed to run properly for the nonlinear models:
      if(is.na(AICnonlineari))cat("NOTE: The optimiser could not converge satisfactorily on the parameters of the non-linear model with intercept.", "\n")
      if(is.na(AICnonlinear))cat("NOTE: The optimiser could not converge satisfactorily on the parameters of the non-linear model without intercept.", "\n")
      if(is.na(AICnonlinear) & is.na(AICnonlinear))cat("NOTE: The failure of the optimiser to converge satisfactorily on the parameters of either non-linear models has lead it to accept a linear model.", "\n")

      # this is one of the arguments of the new COBECOS object function: do you want to see the parameter values?
      if(fitinfo==TRUE){
          cat("\n Parameter values: \n")
          print(.Object@ECFitPar) 
      }
      cat(" ------------------ \n")
    }

    # ------------------------------------------------------------------------ #
    # ------- USER DEFINED EC MODEL WITH SOME OR ----------------------------- #
    # ------- ALL PARAMETERS ESTIMATED FROM DATA ----------------------------- #
    # ------------------------------------------------------------------------ #
    if(!is.na(.Object@control$ECMod)){      
        cat("\n User defined EC model: ")
        
        if(!all(is.na(.Object@ECData))) {
          startslope = (max(.Object@ECData$Cost)-min(.Object@ECData$Cost))/(max(.Object@ECData$Effort)-min(.Object@ECData$Effort))
          startintercept = min(.Object@ECData$Cost)
        }

        if(.Object@control$ECMod == "nonlinear"){
        
          cat("nonlinear \n")
          
          if (any(is.na(.Object@ECFitPar[1:3]))) {   #if any parameters need estimating
            if(!all(is.na(.Object@ECData))) {     #if there is some data
              params<-c(startslope,startintercept,1)
              tempnonlineari<- optim(params, nonlineari ,method="BFGS",  hessian = TRUE)
              if(tempnonlineari$convergence==0) {
                if(is.na(.Object@ECFitPar[1])) {.Object@ECFitPar[1] <- tempnonlineari$par[1] }
                if(is.na(.Object@ECFitPar[2])) {.Object@ECFitPar[2] <- tempnonlineari$par[2] }
                if(is.na(.Object@ECFitPar[3])) {.Object@ECFitPar[3] <- tempnonlineari$par[3] } 
              } else stop("The optimiser did not converge when fitting nonlinear function to cost data")
            } else stop("No cost data to fit nonlinear function \n")
          }
          names(.Object@ECFitPar)<-c("slope","intercept","power") 
          .Object@ECFitFunc<-function(Effort,.Object) .Object@ECFitPar[1]*Effort^.Object@ECFitPar[3]+.Object@ECFitPar[2]          
          #names(.Object@ECFitFunc)<-"nonlinear"
        }

        else if(.Object@control$ECMod == "linear"){
        
          cat("linear \n")
          
          .Object@ECFitPar[3]<-NA
          if (any(is.na(.Object@ECFitPar[1:2]))) {   #if any parameters need estimating
            if(!all(is.na(.Object@ECData))) {     #if there is some data
              params<-c(startslope,startintercept)
              templineari<- optim(params, lineari ,method="BFGS",  hessian = TRUE)
              if(templineari$convergence==0) {
                if(is.na(.Object@ECFitPar[1])) {.Object@ECFitPar[1] <- templineari$par[1] }
                if(is.na(.Object@ECFitPar[2])) {.Object@ECFitPar[2] <- templineari$par[2] } 
              } else stop("The optimiser did not converge when fitting linear function to cost data")
            } else stop("No cost data to fit linear function \n")
          }
          names(.Object@ECFitPar)<-c("slope","intercept") 
          .Object@ECFitFunc <- function(Effort,.Object) .Object@ECFitPar[1]*Effort+.Object@ECFitPar[2]          
          #names(.Object@ECFitFunc)<-"linear"
        }
        else if(.Object@control$ECMod == "step"){
        
             cat("step function \n Note that the BCalc optimisation of effort may not work correctly with this user defined relationship","\n")
             
             .Object@ECFitFunc <- function(Effort,.Object) ifelse(Effort>.Object@ECFitPar[2],.Object@ECFitPar[3],.Object@ECFitPar[1])             
        } 
        else stop("User defined function not recognised")
        
        if(fitinfo==TRUE){
            cat("\n Parameter values: \n ")
          print(.Object@ECFitPar)
        }
        cat(" ------------------ \n")
    }
    
    ########
    # PLOT #
    ########
    if(graphics == TRUE){
    
      windows(width=12)
      par(mfrow = c(1,2)) 
      
      eff<-seq(.Object@Effort_Min,.Object@Effort_Max,length.out=100)
        
      plot(eff,.Object@EPFitFunc(eff,.Object),xlab = "Enforcement effort",ylab="Probability",type = "l", col = "red", pch=19,xlim=c(.Object@Effort_Min,.Object@Effort_Max),ylim=c(0,1),lwd=2)
      try(points(.Object@EPData$Effort,.Object@EPData$Prob,pch=21),silent=TRUE)
      title(paste(.Object@Ename,"EP relationship"))
                
      plot(eff,.Object@ECFitFunc(eff,.Object),xlab = "Enforcement effort",ylab="Cost",type = "l", col = "red", pch=19,xlim=c(.Object@Effort_Min,.Object@Effort_Max),lwd=2)
      try(points(.Object@ECData$Effort,.Object@ECData$Cost,pch=21),silent=TRUE)
      title(paste(.Object@Ename,"EC relationship"))        
    }

    cat("\n \n")
                                                  
  # requirement of S4 class initialisation - ultimately the function returns the new Object
  return(.Object)

})

#------------------------------------------------------------------------------------------------------------------
# METHOD 2
#------------------------------------ BENEFIT CALCULATOR / OPTIMISER ----------------------------------------------

BCalc <- function(.Object,Optimize,stoch,graphics){
     
    cat("\n---------------------------------------------------------------------------- \n")
    cat(" Calculate benefits \n \n")
    flush.console()
    
    # ------------------------------------------------------------------------ #
    # ------- INITIALISE ----------------------------------------------------- #
    # ------------------------------------------------------------------------ #
    
    .Object@StochEffort <- as.numeric(rep(NA,.Object@StochN))
    .Object@StochCost <- as.numeric(rep(NA,.Object@StochN))
    .Object@StochProb <- as.numeric(rep(NA,.Object@StochN))
    .Object@StochExpFine <- as.numeric(rep(NA,.Object@StochN))
    .Object@StochHarvest <- as.numeric(rep(NA,.Object@StochN))
    .Object@StochSocialB <- as.numeric(rep(NA,.Object@StochN))
    .Object@StochPrivateB <- as.numeric(rep(NA,.Object@StochN))
    .Object@HarvestErr <- rep(1,.Object@StochN)
    .Object@EPErr <- rep(0,.Object@StochN)
    .Object@ECErr <- rep(0,.Object@StochN)
     
    # ------------------------------------------------------------------------ #
    # ------- DEFINE FUNCTIONS ----------------------------------------------- #
    # ------------------------------------------------------------------------ # 
    if (is.null(.Object@PrivateBFunc))
    {
      cat("Defined default Private Benefit Function \n")
      ####################################
      # DEFAULT PRIVATE BENEFIT FUNCTION #
      ####################################  
      .Object@PrivateBFunc <- function(Harvest,Price,FCost,Biomass,ExpFine,UnitTax,TAC)
        Price*Harvest-FCost*((Harvest*Harvest)/Biomass)-(ExpFine*Harvest)
    }
    
    #############################
    # FISHING RESPONSE FUNCTION #
    #############################
    .Object@FishingRFunc <- function(Price,FCost,Biomass,ExpFine,UnitTax,TAC)
        optimise(.Object@PrivateBFunc,interval=c(0,Biomass),Price,FCost,Biomass,ExpFine,UnitTax,TAC,maximum=TRUE)$maximum
    
    if (is.null(.Object@SocialBFunc))
    {
      cat("Defined default Social Benefit Function \n \n")
      ###################################
      # DEFAULT SOCIAL BENEFIT FUNCTION #
      ###################################   
      .Object@SocialBFunc <-function(Price,ShadowVB,Harvest,FCost,Biomass,TotalCost,ExpFine,UnitTax,TAC)
        (Price-ShadowVB)*Harvest-FCost*((Harvest*Harvest)/Biomass)-TotalCost
    }
        
    #############################
    # SOCIAL OBJECTIVE FUNCTION #
    #############################   
    .Object@SocialOBJFunc <- function(Eff,rp=1){
        
        
          Cost <- .Object@ECFitFunc(Eff,.Object)+.Object@ECErr[rp]
          Prob <- .Object@EPFitFunc(Eff,.Object)+.Object@EPErr[rp]
          if(Cost<0) Cost<-0
          if(Prob<0) Prob<-0 
         
        ExpFine <-  Prob*.Object@Fine
         
        Harvest <- .Object@FishingRFunc(.Object@Price,.Object@FCost,.Object@Biomass,ExpFine,.Object@UnitTax,.Object@TAC)*.Object@HarvestErr[rp]
        
        .Object@SocialBFunc(.Object@Price,.Object@ShadowVB,Harvest,.Object@FCost,.Object@Biomass,Cost,ExpFine,.Object@UnitTax,.Object@TAC)
    }
    # ------------------------------------------------------------------------ #
    # ------- DEFINE STOCHASTICITY ------------------------------------------- #
    # ------------------------------------------------------------------------ #
    if(!any(stoch==c(0:3))) stop("stoch!=c(0:3)")
    
    .Object@Stochastic <- FALSE
    .Object@ImpError <- FALSE
    .Object@EstError <- FALSE
    
    if (stoch!=0) {
      cat("Stochastic analysis with \n")
      .Object@Stochastic <- TRUE
      if (stoch == 1 || stoch == 3) { .Object@ImpError <- TRUE
      cat("- Implementation error \n") }
      if (stoch == 2 || stoch == 3) { .Object@EstError <- TRUE
      cat("- Estimation error \n") }
      cat("\n")
    }
    # ------------------------------------------------------------------------ #
    # ------- ESTIMATE STD. ERR. AROUND FITTED RELATIONSHIPS ----------------- #
    # ------------------------------------------------------------------------ #
    if(.Object@EstError==TRUE) {
        if(!all(is.na(.Object@ECData))) {     #if there is some data (otherwise it defaults to user inputted value)
          .Object@ECEstErrorSE <- sqrt(mean((.Object@ECFitFunc(.Object@ECData$Effort,.Object)-.Object@ECData$Cost)^2)) #sqrt of error mean squared (biased estimate of SE assuming normality)
        } #else .Object@ECEstErrorSE[i] <- 0
        if(!all(is.na(.Object@EPData))) {    
          .Object@EPEstErrorSE <- sqrt(mean((.Object@EPFitFunc(.Object@EPData$Effort,.Object)-.Object@EPData$Prob)^2))
        } #else .Object@EPEstErrorSE[i] <- 0
    }
    # ------------------------------------------------------------------------ #
    # ------- DETERMINISTIC ANALYSIS ----------------------------------------- #
    # ------------------------------------------------------------------------ #
    if (!.Object@Stochastic) {
                  
      if(Optimize == TRUE){
        cat("Optimisation of Enforcement effort\n")
        .Object@Optimized<-optimise(.Object@SocialOBJFunc,interval=c(.Object@Effort_Min,.Object@Effort_Max),rp=1,maximum=TRUE)
        .Object@Effort <- .Object@Optimized$maximum
      }
     
        .Object@Cost <- .Object@ECFitFunc(.Object@Effort,.Object)
        .Object@Prob <- .Object@EPFitFunc(.Object@Effort,.Object)
        
      .Object@ExpFine <-  .Object@Prob*.Object@Fine
            
      .Object@Harvest <- .Object@FishingRFunc(.Object@Price,.Object@FCost,.Object@Biomass,.Object@ExpFine,.Object@UnitTax,.Object@TAC)
       
      .Object@PrivateB <- .Object@PrivateBFunc(.Object@Harvest,.Object@Price,.Object@FCost,.Object@Biomass,.Object@ExpFine,.Object@UnitTax,.Object@TAC)
      .Object@SocialB <- .Object@SocialBFunc(.Object@Price,.Object@ShadowVB,.Object@Harvest,.Object@FCost,.Object@Biomass,.Object@Cost,.Object@ExpFine,.Object@UnitTax,.Object@TAC)
            
    }
    
    # ------------------------------------------------------------------------ #  
    # ------- STOCHASTIC ANALYSIS -------------------------------------------- #
    # ------------------------------------------------------------------------ #
    if(.Object@Stochastic){  
      
      if(.Object@ImpError & .Object@ImpErrorSE<=0) stop("Stochastic analysis with ImpErrorSE<=0")
      if(.Object@EstError & .Object@EPEstErrorSE<=0) stop("Stochastic analysis with all EPEstErrorSE<=0")
      if(.Object@EstError & .Object@ECEstErrorSE<=0) stop("Stochastic analysis with all ECEstErrorSE<=0")
      
      if (Optimize == TRUE) {
        cat("Optimisation of Enforcement effort \n    % complete \n   ")
        PerCentComplete<-1:10*(.Object@StochN/10)
        counter = 1
      }
      
      if (.Object@ImpError)
      {
        ###################################
        # LOG-NORMAL ERROR AROUND HARVEST #
        ################################### 
        .Object@HarvestErr <- rlnorm(n=.Object@StochN,sdlog=.Object@ImpErrorSE)
      }
      
      if (.Object@EstError)
      {
        #######################################################
        # NORMAL ERROR AROUND PREDICTED PROBILITIES AND COSTS #
        #######################################################                   
           .Object@EPErr <- rnorm(n=.Object@StochN,sd=.Object@EPEstErrorSE)
           .Object@ECErr <- rnorm(n=.Object@StochN,sd=.Object@ECEstErrorSE) 
             
      }
                       
      for(k in 1:.Object@StochN){  # loop through stochastic samples
                                     
        if(Optimize == TRUE){
                   
          .Object@Optimized<-optimise(.Object@SocialOBJFunc,interval=c(.Object@Effort_Min,.Object@Effort_Max),rp=k,maximum=TRUE)
          .Object@StochEffort[k] <- .Object@Optimized$maximum 
              
            .Object@StochCost[k] <- .Object@ECFitFunc(.Object@StochEffort[k],.Object)+.Object@ECErr[k]
            .Object@StochProb[k] <- .Object@EPFitFunc(.Object@StochEffort[k],.Object)+.Object@EPErr[k]
            if(.Object@StochCost[k]<0) .Object@StochCost[k]<-0
            if(.Object@StochProb[k]<0) .Object@StochProb[k]<-0                 
                
          .Object@StochExpFine[k] <-  .Object@StochProb[k]*.Object@Fine 
                        
          .Object@StochHarvest[k] <- .Object@FishingRFunc(.Object@Price,.Object@FCost,.Object@Biomass,.Object@StochExpFine[k],.Object@UnitTax,.Object@TAC)*.Object@HarvestErr[k]

          .Object@StochSocialB[k] <- .Object@Optimized$objective                 
          .Object@StochPrivateB[k] <- .Object@PrivateBFunc(.Object@StochHarvest[k],.Object@Price,.Object@FCost,.Object@Biomass,.Object@StochExpFine[k],.Object@UnitTax,.Object@TAC)
          
          if(k==PerCentComplete[counter]){  # percent complete meter.. updates every 10 per cent
                cat(paste(" ...",as.character(counter*10)," ",sep=""))
                counter=counter+1
                flush.console()
          } 
        }                              
        else {
            
            .Object@StochEffort[k] <- .Object@Effort
            
                .Object@StochCost[k] <- .Object@ECFitFunc(.Object@StochEffort[k],.Object)+.Object@ECErr[k]
                .Object@StochProb[k] <- .Object@EPFitFunc(.Object@StochEffort[k],.Object)+.Object@EPErr[k]
                if(.Object@StochCost[k]<0) .Object@StochCost[k]<-0
                if(.Object@StochProb[k]<0) .Object@StochProb[k]<-0                  
                 
            .Object@StochExpFine[k] <- .Object@StochProb[k]*.Object@Fine 
                                  
            .Object@StochHarvest[k] <- .Object@FishingRFunc(.Object@Price,.Object@FCost,.Object@Biomass,.Object@StochExpFine[k],.Object@UnitTax,.Object@TAC)*.Object@HarvestErr[k]

            .Object@StochSocialB[k] <- .Object@SocialBFunc(.Object@Price,.Object@ShadowVB,.Object@StochHarvest[k],.Object@FCost,.Object@Biomass,.Object@StochCost[k],.Object@StochExpFine[k],.Object@UnitTax,.Object@TAC)
            .Object@StochPrivateB[k] <- .Object@PrivateBFunc(.Object@StochHarvest[k],.Object@Price,.Object@FCost,.Object@Biomass,.Object@StochExpFine[k],.Object@UnitTax,.Object@TAC)            
          }     
        }  # end of loop of stochastic samples
        
      .Object@Effort <- median(.Object@StochEffort,na.rm=TRUE)  #set to median of stochastic samples  (NA's indicate that optimisation failed)
      .Object@Cost <- median(.Object@StochCost,na.rm=TRUE)
      .Object@Prob <- median(.Object@StochProb,na.rm=TRUE)
      .Object@ExpFine <- median(.Object@StochExpFine,na.rm=TRUE)
      .Object@Harvest <- median(.Object@StochHarvest,na.rm=TRUE)
      .Object@PrivateB <- median(.Object@StochPrivateB,na.rm=TRUE)
      .Object@SocialB <- median(.Object@StochSocialB,na.rm=TRUE)
       
      cat("\n")
    }
    # ------------------------------------------------------------------------ #    
    # ------- CONSOLE OUTPUT ------------------------------------------------- #
    # ------------------------------------------------------------------------ #
    if(.Object@Stochastic){             
      cat("\n")
      cat("Private benefit: \n")
      cat(paste("Median = ", round(.Object@PrivateB,digits=1), "   95% quantiles = ",round(quantile(na.omit(.Object@StochPrivateB),0.025),digits=1),",",round(quantile(na.omit(.Object@StochPrivateB),0.975),digits=1),"\n \n",sep=""))
      
      cat("Social benefit: \n")
      cat(paste("Median = ", round(.Object@SocialB,digits=1), "   95% quantiles = ",round(quantile(na.omit(.Object@StochSocialB),0.025),digits=1),",",round(quantile(na.omit(.Object@StochSocialB),0.975),digits=1),"\n \n",sep=""))
      
      cat(paste("Harvest: \n",sep=""))
      cat(paste("Median = ", round(.Object@Harvest,digits=1),"   95% quantiles = ",round(quantile(na.omit(.Object@StochHarvest),0.025),digits=1),",",round(quantile(na.omit(.Object@StochHarvest),0.975),digits=1),"\n \n",sep=""))
          
      if (Optimize == TRUE){
          cat("Optimized effort levels [",.Object@Effort_Min,",",.Object@Effort_Max,"]: \n",sep="")      
          cat("Median = ",signif(.Object@Effort,digits=3),"   95% quantiles = ", round(quantile(na.omit(.Object@StochEffort),0.025),digits=3),",",round(quantile(na.omit(.Object@StochEffort),0.975),digits=3),"\n",sep="")
      }
      else {
          cat("User defined effort levels: \n")
          cat(round(.Object@Effort,digits=3),"\n")
      }
    } 
    
    if(!.Object@Stochastic){
        cat("\n")            
        cat("Private benefit: \n")
        cat(round(.Object@PrivateB,digits=1),"\n \n")
        
        cat("Social benefit: \n")
        cat(round(.Object@SocialB,digits=1),"\n \n")
        
        cat("Harvest: \n",sep="")
        cat(round(.Object@Harvest,digits=1),"\n \n")
            
        if (Optimize == TRUE){
            cat("Optimized effort levels [",.Object@Effort_Min,",",.Object@Effort_Max,"]: \n",sep="")
        }
        else{
            cat("User defined effort levels:\n")
        }
        cat(round(.Object@Effort,digits=3),"\n")
    }
    cat("\n")
    
    # ------------------------------------------------------------------------ #
    # ------- GRAPHICS OUTPUT ------------------------------------------------ #
    # ------------------------------------------------------------------------ #
    if (graphics==TRUE) {
      if (.Object@Stochastic) {
          
          windows()
          par(mfrow = c(3,1),mai=c(0.55,0.55,0.35,0.05))  # set up the plot array
                        
          r <- density(na.omit(.Object@StochPrivateB),adjust=1.1,n=800)
          xl <- c(quantile(na.omit(.Object@StochPrivateB),0.001),quantile(na.omit(.Object@StochPrivateB),0.999))
          x2 <- c(quantile(na.omit(.Object@StochPrivateB),0.025),quantile(na.omit(.Object@StochPrivateB),0.975))
          x3 <- quantile(na.omit(.Object@StochPrivateB),0.5)
          plot(r,main = "Private profits",xlab="",xlim=xl,col="red",ylab="Relative Frequency",lwd=2)
          abline(v=x2,col="blue")
          abline(v=x3,col="blue",lwd=2)
          
          r <- density(na.omit(.Object@StochSocialB),adjust=1.1,n=800)
          xl <- c(quantile(na.omit(.Object@StochSocialB),0.001),quantile(na.omit(.Object@StochSocialB),0.999))
          x2 <- c(quantile(na.omit(.Object@StochSocialB),0.025),quantile(na.omit(.Object@StochSocialB),0.975))
          x3 <- quantile(na.omit(.Object@StochSocialB),0.5)
          plot(r,main = "Social benefits",xlab="",xlim=xl,col="red",ylab="Relative Frequency",lwd=2)
          abline(v=x2,col="blue")
          abline(v=x3,col="blue",lwd=2)    
          
          r <- density(na.omit(.Object@StochHarvest),adjust=1.1,n=800)
          xl <- c(quantile(na.omit(.Object@StochHarvest),0.001),quantile(na.omit(.Object@StochHarvest),0.999))
          x2 <- c(quantile(na.omit(.Object@StochHarvest),0.025),quantile(na.omit(.Object@StochHarvest),0.975))
          x3 <- quantile(na.omit(.Object@StochHarvest),0.5)
          plot(r,main = "Harvest",xlab="",xlim=xl,col="red",ylab="Relative Frequency",lwd=2)
          abline(v=x2,col="blue")
          abline(v=x3,col="blue",lwd=2)
      }   
        
      #######################################
      # PLOT OF SOCIAL BENEFIT OPTIMISATION #
      #######################################
      if (.Object@Stochastic) {
          windows(width=6,height=18)
          par(mfcol=c(3,1))
      }
      if (!.Object@Stochastic) {
          windows(width=6,height=12)
          par(mfcol=c(2,1))
      }
                      
            # set up temporary effort array for social benefit profiles
            npoints<-100    # resolution of the social benefit profile
            plot_sb<-1:npoints    # set up the initial vectors
            plot_hvst<-1:npoints 
            Ef<- seq(.Object@Effort_Min,.Object@Effort_Max,length.out=npoints)   # alter the ith (enforcement type) array of effort
            for(w in 1:npoints){         # loop over effort
              if(!.Object@Stochastic) {            
        
                  Cost <- .Object@ECFitFunc(Ef[w],.Object)
                  Prob <- .Object@EPFitFunc(Ef[w],.Object)
                
         
                ExpFine <-  Prob*.Object@Fine
         
                plot_hvst[w] <- .Object@FishingRFunc(.Object@Price,.Object@FCost,.Object@Biomass,ExpFine,.Object@UnitTax,.Object@TAC)
                plot_sb[w] <- .Object@SocialBFunc(.Object@Price,.Object@ShadowVB,plot_hvst[w],.Object@FCost,.Object@Biomass,Cost,ExpFine,.Object@UnitTax,.Object@TAC)              
              }
              if(.Object@Stochastic) {       
               
                Harvest <- as.numeric(rep(NA,.Object@StochN))
                SocialB <- as.numeric(rep(NA,.Object@StochN))
                for(k in 1:.Object@StochN){  # loop through stochastic samples
                                            
        
                    Cost <- .Object@ECFitFunc(Ef[w],.Object)+.Object@ECErr[k]
                    Prob <- .Object@EPFitFunc(Ef[w],.Object)+.Object@EPErr[k]
                    if(Cost<0) Cost<-0
                    if(Prob<0) Prob<-0 
         
                  ExpFine <-  Prob*.Object@Fine
         
                  Harvest[k] <- .Object@FishingRFunc(.Object@Price,.Object@FCost,.Object@Biomass,ExpFine,.Object@UnitTax,.Object@TAC)*.Object@HarvestErr[k]
                  SocialB[k] <- .Object@SocialBFunc(.Object@Price,.Object@ShadowVB,Harvest[k],.Object@FCost,.Object@Biomass,Cost,ExpFine,.Object@UnitTax,.Object@TAC)
                }  # end of loop over stochastic samples
              plot_hvst[w]<-median(Harvest)  # integrate over stochastic uncertainty
              plot_sb[w]<-median(SocialB)
              }
            }
             
            plot(Ef,plot_sb,type="l",xlab = "Enforcement effort",ylab="Social benefits",lwd=2) # plot the marginal social benefits
            if(Optimize==TRUE) {
              txt.x<-.Object@Effort
              txt.y<-min(plot_sb,na.rm=TRUE)+(max(plot_sb,na.rm=TRUE)-min(plot_sb,na.rm=TRUE))/2
              text(txt.x,txt.y,labels="Optimised Effort level",srt=90,col=4,adj=c(0.5,-0.5))
            }
            abline(v=Ef[match(max(plot_sb),plot_sb)],col="red",lty="dotted") # plot the marginal effort where social benefits are maximised if optimization has not been undertaken
            abline(v=.Object@Effort,col="blue",lwd=2)
            if(.Object@Stochastic) {
              abline(v=c(quantile(na.omit(.Object@StochEffort),0.025),quantile(na.omit(.Object@StochEffort),0.975)),col="blue")
            }
            title("Social benefit profile")
            
            plot(Ef,plot_hvst,type="l",xlab = "Enforcement effort",ylab="Harvest",lwd=2)
            abline(v=Ef[match(max(plot_sb),plot_sb)],col="red",lty="dotted")
            abline(v=.Object@Effort,col="blue",lwd=2)
            if(.Object@Stochastic) {
              abline(v=c(quantile(na.omit(.Object@StochEffort),0.025),quantile(na.omit(.Object@StochEffort),0.975)),col="blue")
            }            
            title("Harvest response")
          
          if (.Object@Stochastic) {   
            x1 <- c(quantile(na.omit(.Object@StochEffort),0.025),quantile(na.omit(.Object@StochEffort),0.975))
            x2 <- quantile(na.omit(.Object@StochEffort),0.5)
            hist(na.omit(.Object@StochEffort),breaks=seq(0,1,by=0.01),xlim=c(.Object@Effort_Min,.Object@Effort_Max),main = "Distribution of Effort",xlab="Enforcement Effort",ylab="Frequency",yaxt="n")
            abline(v=x1,col="blue")
            abline(v=x2,col="blue",lwd=2)
           } 
            
          # End of loop of enforcement type
                     
          ########################################
          # PLOT OF PRIVATE BENEFIT OPTIMISATION #
          # (FISHING RESPONSE FUNCTION)          #
          ########################################
          hvst<-seq(0,.Object@Biomass,length.out=npoints)
          windows(height=12)
          par(mfrow=c(2,1))
          plot_sb<-1:npoints    # set up the initial vectors
          plot_pb<-1:npoints
          for(w in 1:npoints){
            if (!.Object@Stochastic) {
              plot_pb[w]<-.Object@PrivateBFunc(hvst[w],.Object@Price,.Object@FCost,.Object@Biomass,.Object@ExpFine,.Object@UnitTax,.Object@TAC)
              plot_sb[w] <- .Object@SocialBFunc(.Object@Price,.Object@ShadowVB,hvst[w],.Object@FCost,.Object@Biomass,.Object@Cost,.Object@ExpFine,.Object@UnitTax,.Object@TAC)
            }
            if (.Object@Stochastic) {
            
              PrivateB <- as.numeric(rep(NA,.Object@StochN))
              SocialB <- as.numeric(rep(NA,.Object@StochN))
              for(k in 1:.Object@StochN){  # loop through stochastic samples
                                         
                    Cost <- .Object@ECFitFunc(.Object@Effort,.Object)+.Object@ECErr[k]
                    Prob <- .Object@EPFitFunc(.Object@Effort,.Object)+.Object@EPErr[k]
                    if(Cost<0) Cost<-0
                    if(Prob<0) Prob<-0
         
                  ExpFine <-  Prob*.Object@Fine
                  
                  PrivateB[k] <- .Object@PrivateBFunc(hvst[w],.Object@Price,.Object@FCost,.Object@Biomass,ExpFine,.Object@UnitTax,.Object@TAC)
                  SocialB[k] <- .Object@SocialBFunc(.Object@Price,.Object@ShadowVB,hvst[w],.Object@FCost,.Object@Biomass,Cost,ExpFine,.Object@UnitTax,.Object@TAC)
              }
              plot_pb[w]<-median(PrivateB)   # integrate over stochastic uncertainty     
              plot_sb[w]<-median(SocialB)
            }
          }
                      
          plot(hvst,plot_pb,type="l",xlab="Harvest",ylab="Private Benefit",main="Private benefit profile",col=2,lwd=2)
          txt.x<-.Object@Harvest
          txt.y<-min(plot_pb,na.rm=TRUE)+(max(plot_pb,na.rm=TRUE)-min(plot_pb,na.rm=TRUE))/2
          text(txt.x,txt.y,labels="Optimised Harvest level",srt=90,col=4,adj=c(0.5,-0.5))
          abline(v=hvst[match(max(plot_pb),plot_pb)],col="red",lty="dotted")
          abline(v=.Object@Harvest,col="blue",lwd=2)
          if(.Object@Stochastic) {
              abline(v=c(quantile(na.omit(.Object@StochHarvest),0.025),quantile(na.omit(.Object@StochHarvest),0.975)),col="blue")
          }

          plot(hvst,plot_sb,type="l",xlab = "Harvest",ylab="Social benefits",main="Social benefit profile",lwd=2,col=3) # plot the marginal social benefits
          abline(v=hvst[match(max(plot_sb),plot_sb)],col="red",lty="dotted")
          abline(v=.Object@Harvest,col="blue",lwd=2)
          if(.Object@Stochastic) {
              abline(v=c(quantile(na.omit(.Object@StochHarvest),0.025),quantile(na.omit(.Object@StochHarvest),0.975)),col="blue")
          }
    }
    # ------------------------------------------------------------------------ #
    
    return(.Object)   #The function returns the new object. 
}

cat("\n Source code loaded successfully \n \n")