setwd("~/Gavin's Dropbox/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/COBECOS")

# -------------------------------------------------------------------------------------------- #
# Tutorial: Hake
# -------------------------------------------------------------------------------------------- #
# Notes: PLEASE ENSURE YOU HAVE INSTALLED THE LATEST VERSION OF R BEFORE RUNNING 
# THIS TUTORIAL; tutorial will work with data provided in previous releases.
# -------------------------------------------------------------------------------------------- #

#  first load the source code (adjust file path as necessary)
source("COBECOScode_v1_3.r")

#  ------------------------- Data -------------------------------------------------------- 
#
#  Before continuing, it is a good idea to define a variable which represents the current 
#  working directory (the location on disk where files and objects will be loaded to and from):

HakePath <- "~/Gavin's Dropbox/Dropbox/SFG/Lennon/Fisheries Funding Paper/Model/COBECOS/"

#  Examine the control file ('Control.csv') and datafiles 
#  ('Cost.csv' and 'Prob.csv'). 
#
#  The optimisers employed in this package are relatively robust. However, if the 
#  optimiser fails to fit a non-linear function it will default to a linear relationship 
#  between enforcement effort and enforcement cost that may not be a satisfactory description
#  of the observed relationship. 

#############################################################
# Tutorial Part I:                                          #
# ----------------                                          #
# INTRODUCTION TO METHODS 1 AND 2                           #
#############################################################

#  ----------------------- Method 1 ------------------------------------------------------ 
#
#  Create an instance of a COBECOS object (which calls the initialisation method)
#
#  There are three arguments: 
#  (i)   path      (the location on your computer of the control file and datafiles 
#                  (currently reprsented by 'HakePath')) 
#  (ii)  fitinfo   (whether you would like to see detailed information about the selected models
#                  and their fit to the observed data - assuming you do not specify all of them)
#  (iii) graphics  (whether you would like to see a graphical output of the fit of the models
#                  to the data (where fitted) or the specified relationship (where user  
#                  defined))
#
#  As input the initialisation method will read the control file ('/COBECOS_prog/Tutorial/Hake/Control.csv') 
#  and datafiles ('/COBECOS_prog/Tutorial/Hake/Cost.csv' and '/COBECOS_prog/Tutorial/Hake/Prob.csv'). Using
#  the control file, functional forms and parameter values for the effort-cost (EC) and effort-probability (EP)
#  relationships can be specified by the user. 
#  EP relationships can be one of:
#   - logistic
#   - exponential
#  EC relationships can be one of:
#   - linear  
#   - nonlinear
#  If the functional form is specified (e.g. by writing 'logisitic' or 'exponential' 
#  under EPMod/ECMod in the control file) then any parameters marked as NA will be
#  estimated from the data, whilst other parameters will remain fixed. If the model 
#  type is marked as NA then user-specified parameters will be ignored and the model 
#  selected (by the AIC) using the data with all parameters estimated. Although it is possible to input
#  a step function this is not advised as it may lead to problems with the optimisation.

#  Initialise the COBECOS object:

Hake<-new("COBECOS", path = HakePath, fitinfo = TRUE, graphics = TRUE)

#  The COBECOS program has fitted EP and EC relationships dependent on the specifications
#  of the control file. To familiarise yourself with its workings, it is advisable to change
#  some of the control file settings and examine the initialisation results (with fitinfo = TRUE and 
#  graphics = TRUE). 
#
#  Take a look at your new COBECOS object:

Hake

#  A number of names (these represent 'slots' in the COBECOS object Hake) are listed along with 
#  their corresponding value and class(including numeric values, a vectors of numeric values, 
# 'data.frames', lists and character strings). 
#
#  For example, the Hake@control slot contains all of the information from the Control.csv 
#  control file:

Hake@control

#  Many of the slots are parameters of benefit functions or input data which are initialised 
#  to a default setting.
#
#  To obtain a more concise list of the slots of your COBECOS object you can use the function:

slotNames(Hake)

#  Note: clearly the functional form of effort-cost and effort-prob relationships will strongly 
#  affect the subsequent analyses. Take some time to check these. 
#
#  Where applicable, examine the fit of the graphs to the observed data. Subjectively, does the 
#  fitted line follow the approximate pattern of the data?  If not, you can always fit your own 
#  models to data and specify them in the control file (/Control.csv).
#
#  Where applicable, study your user defined effort-cost and effort-prob relationships. Are they
#  consistent with your desired specification? 
#
#  If you are satisfied with the fit of the models to the data (and/or that the relationships 
#  reflect those specified in Control.csv), before going further, it is worth taking a look at 
#  the effort-cost, effort-prob graphs for each enforcement type to get an idea of their 
#  relative cost efficacy. 
# 
#
#
#
#  ---------------------- Method 2 ------------------------------------------------------- 
#
#  Calculate the social benefits, private profits and compliance (harvesting) given the fitted 
#  or specified relationships between enforcement effort and cost of enforcement and probabality
#  of detecting infringements.
#
#  There are four arguments:
#  (i)   the name of an initialised COBECOS object
#  (ii)  Optimize (whether you would like to allow the program to find the enforcement efforts
#                 which maximise social benefits. Where this is FALSE, the program will use 
#                 those specified (in this case Hake@Effort)) (note Optimize is case sensitive)
#  (iii) stoch    0:  no stochasticity 
#                 1:  lognormal error in the fishers behavioural function 
#                 2:  normal error around effort-cost and effort-prob relationships
#                 3:  1+2 combined 
#  (iv)  graphics (whether you would like to see graphical outputs, such as the user-defined or
#                 fitted effort (where Optimize = TRUE) and the distribution of Social benefits
#                 private profits and harvesting (where stoch = TRUE)
#
#

# User values must be specified for the following bioeconomic parameters for the case study
Hake@Price <- 1		#  price of harvest    (p)
Hake@FCost <- 0.5	#  fishing cost parameter  (c)
Hake@Biomass <- 1000	#  the biomass available to fishing   (x)
Hake@Fine <- 1		#  the level of fining  (f)
Hake@ShadowVB <- 0.2	# The shadow value of biomass

#  In the first instance we will simply use the fitted effort-cost and effort-prob 
#  relationships of our Hake COBECOS object to make a deterministic calculation of harvest, 
#  social benefits and private profits at the default levels of effort, shadow value of biomass,
#  fining, harvest price and cost of fishing. 

# By default the effort value is half-way between Effort_Min and Effort_Max
Hake@Effort_Min
Hake@Effort_Max
Hake@Effort

# Now call BCalc (without optimisation)
Hake<-BCalc(Hake,Optimize = FALSE,stoch = 0,graphics = TRUE)

#  Two graphical output windows are produced. One shows the relationship between the harvest and the 
#  private and social benefits. The relationship between harvest and private benefit is important since
#  it illustrates the first stage optimisation that calculates the fishers response to their economic
#  surroundings. Rather than being defined analytically the fishing response function is
#  a simple optimisation (over harvest) of the private benefit function (see manual equation 9). 

Hake@PrivateBFunc
Hake@FishingRFunc 

#  The economic model assumes that fishers are rationale and omnipotent and will
#  behave in a way that maximises their benefit. Harvest should therefore be equal to the level that 
#  gives the maximum private benefit. This level is shown by the blue vertical line of the graphical output. The user should check 
#  that it does indeed equate to the maximum (shown by the red dotted line). Next to the private benefit function is 
#  plotted the social benefit profile. This is simply for comparison purposes and may prove insightful during 
#  optimisation of the enforcement level.
#  The second graphical output primarily illustrates the relationship between enforcement effort and social benefit
#  for each enforcement type. On the top 
#  row is shown the social benefit profile against effort. The utility of this plot will beocome apparent 
#  when enforcement effort is optimised. The harvest response function (i.e. how harvest levels
#  will change at different enforcement levels) is also shown for illustrative purposes. The blue vertical 
#  line represents the user defined (or optimised) levels of enforcement effort.  The red dotted 
#  vertical line represents the marginal level of enforcement effort that corresponds with the 
#  maximum social benefit. 
#
#  The level of harvesting, social benefits and private profits are stored in the slots:

Hake@Harvest
Hake@SocialB
Hake@PrivateB


#  You can change the levels of these enforcement effort manually and re-run the same BCalc 
#  function:

Hake@Effort<- 0.2
Hake<-BCalc(Hake,Optimize = FALSE,stoch = 0,graphics = TRUE)

#  The BCalc function includes an optimisation routine to maximize social benefits by changing 
#  the levels of enforcement effort.

Hake<-BCalc(Hake,Optimize = TRUE,stoch = 0,graphics = TRUE)

#  Notice that the graphs have updated to reflect the socially optimal effort levels.  If the 
#  optimiser has worked correctly, these optimised levels of effort should correspond with the 
#  maxima of each marginal social benefit profile.  
#
#  Changing the value of slots which represent parameters of the enforcement model has strong 
#  implications for optimal effort in this example. E.g. shadow value of biomass (lambda):
#
#  This will return the shadow value of biomass user set earlier (line 151):

Hake@ShadowVB

#  You can easily examine the consequences for socially optimal effort of changes in this 
#  parameter: 

Hake@ShadowVB <- 0.8
Hake<-BCalc(Hake,Optimize = TRUE,stoch = 0,graphics = TRUE)

#  The graphical output may cause concern as it may be clear that the optimisation was
#  not successful. Optimisation can be dependent on the value of the starting parameters.
#  A solution therefore is to re-set the effort levels and re-run the optimisation.

#  Note that every time that BCalc is called with the argument graphics = TRUE, a new graphics 
#  window is created allowing fitted efforts to be compared among sequential runs (unfortunately 
#  there is no way to return the name of an object as it is being created so you will have to 
#  keep track of these.  You could paste the images into a word document for example).
#
#  Clearly, there are other key parameters that are worth investigating:

Hake@Price
Hake@Fine 
Hake@FCost
Hake@Biomass

#############################################################
# Tutorial Part II:                                         #
# -----------------                                         #
# SPECIFICATION OF THE SOCIAL AND PRIVATE BENEFIT FUNCTIONS #
#############################################################

#  Private benefit and social benefit functions are
#  defined explicitly in the COBECOS object and initialised when BCalc is called
#  provided a functional form has not already been defined.

Hake@PrivateBFunc
Hake@SocialBFunc

#  Rember that the fishing response function is
#  a simple optimisation (over harvest) of the private benefit function.

Hake@FishingRFunc

#  Any Private benefit function can be defined by the user and FishingRFunc
#  will search over harvest levels to maximise the private benefit. It is
#  important that the private benefit function is smooth and that the optimisation
#  does not get trapped on a local maxima. The plots make it easy to assess
#  whether these condistions have been me. Two slots delimit the harvest levels
#  over which the optimiser searches.

#  The private benefit function is defined in terms of the following input values:
#  Harvest; Price; FCost; Biomass; ExpFine; UnitTax; TAC
#  It is important that any user defined function includes all these variables as
#  inputs (in the correct order) although not all have to be used.
#  In its current form for example, the private benefit function does not make use
#  of UnitTax or TAC.

#  To illustrate the potential versatility of the program we begin by 
#  introducing a simple tax on catches

Hake@PrivateBFunc <- function(Harvest,Price,FCost,Biomass,ExpFine,UnitTax,TAC)    
    Price*Harvest-FCost*((Harvest*Harvest)/Biomass)-(ExpFine*Harvest)-(UnitTax*Harvest)
        
#  Remember to set the UnitTax to a level >0

Hake@UnitTax
Hake@UnitTax<-0.3

#  Now re-optimise the social benefit creating a new COBECOS object in the process

Hake_TAX<-BCalc(Hake,Optimize = TRUE,stoch = 0,graphics = TRUE)
   
#  As might have been expected, introducing a tax of this form (effectively reducing the value
#  of the resource) results in a smaller private profit and social benefit, with less effort
#  expended on enforcement. But surely taxation should be of overall social benefit? Let's
#  change the social benefit function to reflect this. Note that all function input variables are retained
#  in the correct order in the new specification.

Hake_TAX@SocialBFunc <- function(Price,ShadowVB,Harvest,FCost,Biomass,TotalCost,ExpFine,UnitTax,TAC)
        (Price-ShadowVB)*Harvest-FCost*((Harvest*Harvest)/Biomass)-TotalCost+(Harvest*UnitTax)
        
Hake_TAX<-BCalc(Hake_TAX,Optimize = TRUE,stoch = 0,graphics = TRUE)

#  Illegal harvest now has an additional social benefit (through taxation) and the optimal
#  enforcment is at the minimum to reflect this. The model is clearly flawed, since
#  illegal catch is unlikely to be taxed and so no social benefit would be accrued. Lets
#  update the model again so that taxation is only on the legal catch, whilst illegal
#  catches are of no wider social benefit.
#  User must now also specify the TAC which is now used in private benefit function:
Hake_TAX@TAC <- 500

Hake_TAX@PrivateBFunc <- function(Harvest,Price,FCost,Biomass,ExpFine,UnitTax,TAC)
        if (Harvest <= TAC) { Price*Harvest-FCost*((Harvest*Harvest)/Biomass)-(UnitTax*Harvest)
        } else  Price*Harvest-FCost*((Harvest*Harvest)/Biomass)-(ExpFine*(Harvest-TAC))
        
Hake_TAX@SocialBFunc <- function(Price,ShadowVB,Harvest,FCost,Biomass,TotalCost,ExpFine,UnitTax,TAC)
        if (Harvest <= TAC) { (Price-ShadowVB)*Harvest-FCost*((Harvest*Harvest)/Biomass)-TotalCost+(Harvest*UnitTax)
        } else  (Price-ShadowVB)*Harvest-FCost*((Harvest*Harvest)/Biomass)-TotalCost+(TAC*UnitTax)

#  Adding such a disconintuity to the Social or Private benefit functions should only
#  be done with great care and be accompanied by close inspection of the graphical
#  outputs. This is necessary to ensure that there is no unusual behaviour and that
#  optimisation was successful.

Hake_TAX<-BCalc(Hake_TAX,Optimize = TRUE,stoch = 0,graphics = TRUE) 

#  The Harvest now reflects total legal and illegal catch.

Hake_TAX@Biomass
Hake_TAX@TAC

#  The optimal effort levels have been returned 
#  showing that the illegal catch is now less (Harvest minus the TAC) and the 
#  social benefit has increased. 

#############################################################
# Tutorial Part III:                                        #
# ------------------                                        #
# STOCHASTIC ANALYSIS                                       #
#############################################################

#  It is advisable to read the COBECOS manual (Section 7) before proceeding with
#  with this section of the tutorial.
# First we initialise a new object
Hake<-new("COBECOS", path = HakePath, fitinfo = TRUE, graphics = FALSE)

#  There are a number of slots in the COBECOS object for the stochastic analysis.

Hake@Stochastic
Hake@ImpError
Hake@EstError

#  These are logical values for the stochastic analysis specifying whether or not there has been any
#  implementation error (a log-normal distribution around harvest) or estimation
#  error (a normal distribution around the predicted probability and cost values).

Hake@ImpErrorSE
Hake@EPEstErrorSE
Hake@ECEstErrorSE

#  These standard errors define error distributions. The implementation error SE is 
#  defined by the user. For the estimation error SE, if there is data then it will
#  be estimated, otherwise the user defined values are used. 

Hake@HarvestErr
Hake@EPErr
Hake@ECErr

#  Contain error values for each sample drawn during a stochastic analysis. Harvest
#  error is multiplicative and sampled from a log-normal distribution. Estimation
#  error is additive and sampled from a normal distribution ensuring non-negative
#  values.

Hake@StochEffort
Hake@StochHarvest
Hake@StochPrivateB
Hake@StochSocialB

#  First run a stochastic analysis with implementation error. The default SE for the 
#  error distribution is in

Hake@ImpErrorSE

#  And the number of stochastic samples in

Hake@StochN
   
#  During the run, the HarvestErr, StochHarvest and StochEffort slots will 
#  be populated. Harvest and Effort are calculated from the median of these 
#  sample distributions.

#  Again we must specify the following bioeconomic parameters for this new object:
Hake@Price <- 1
Hake@FCost <- 0.5
Hake@Biomass <- 1000
Hake@Fine <- 1
Hake@ShadowVB <- 0.2
Hake@Effort <- 0.5

# And call BCalc, this time with the stochasticity argument set to 1
Hake_STOCH1<-BCalc(Hake,Optimize = FALSE,stoch = 1,graphics = FALSE)

Hake_STOCH1@StochHarvest
Hake_STOCH1@HarvestErr
Hake_STOCH1@StochEffort
Hake_STOCH1@Harvest
Hake_STOCH1@Effort

#  Note that StochEffort is constant because optimisation is switched off.
#  The slots used to define the stochasticity have also been updated

Hake_STOCH1@Stochastic
Hake_STOCH1@ImpError
Hake_STOCH1@EstError

#  Running BCalc with Estimation error first estimates the standared error
#  around any fitted EC or EP reltionship using the data. These standard errors
#  are stored in EPEstErrSE and ECEstErrSE. If there is no data 
#  then the standard errors are set equal to user defined values (default is zero). 
#  During the run, the EPErr, 
#  ECErr and StochEffort slots will be populated. StochHarvest is 
#  also recorded. Again, Harvest and Effort are calculated 
#  from the median of these sample distributions.

Hake_STOCH2<-BCalc(Hake,Optimize = FALSE,stoch = 2,graphics = FALSE)

Hake_STOCH2@StochHarvest
Hake_STOCH2@EPErr
Hake_STOCH2@ECErr
Hake_STOCH2@StochEffort
Hake_STOCH2@Harvest
Hake_STOCH2@Effort

Hake_STOCH2@Stochastic
Hake_STOCH2@ImpError
Hake_STOCH2@EstError

#  Now let's do a complete run including both types of stochasticity and with 
#  optimisation switched on

Hake_STOCH3<-BCalc(Hake,Optimize = TRUE,stoch = 3,graphics = TRUE)

#  A new row has been added to the bottom of the graphics window showing the relationship
#  between social benefit/harvest and enforcement effort. This shows the distribution
#  of effort values around the median. 95% quantiles of this distribution are shown as
#  thin blue lines. A third graphical window is also now produced, showing the distribution
#  of the harvest, social and private benefits across stochastic samples.

#  To ensure successful optimisation it is important to examine the graphical outputs.
#  Particularly the social benefit curves. For each stochastic sample the effort is
#  optimised, giving a distribution of optimised effort values from which a median and
#  quantile values can be obtained. 95% quantiles are shown in the graphical output. If 
#  optimisation was unsuccessful it may be necessary to run the program a second time. A
#  second optimisation will automatically be started from the current Effort levels. If the 
#  error distribution SE's are high it may be necessary to increase the stochastic sample
#  size
