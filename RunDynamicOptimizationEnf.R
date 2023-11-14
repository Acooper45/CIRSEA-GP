######################################
# Run Dynamic Optimization --------------------------------------------------
# Solves for optimal policy function f (as function of bvec) given model parameters.
# last input argument "tol" is convergence tolerance (use tol=.01)
######################################

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
      t= t+1 # t is the equivalent of i in a for loop; we need to define the state changes each iteration for basically every variable
      V= Vnew
      oldf1= f1
      olde1 = e1
      for (i in 1:length(bvec)) #for-loop in the while-loop which is pretty fucked
      {
        b= bvec[i] #bvec is a vector of biomass that we're going to iterate, from 0.1 to K
        if(i==1 | i==2) # if i equals one, it says do this particular guess (avoid zero values). otherwise run the function below
          {guess= f1[1]
          guessEnf=e1[1]}
        else
        {guess= f1[i-1] # a good starting guess for this period is the value from last period; iterate!
        guessEnf=e1[i-1]}

          f1 =  rep(r/2,length(bvec))
          FishOut= optim(par=guessEnf, # guess vec (xo = ...)
                         fn=profitOptimEnf, # fn is f.eval (your function)
                         lower=0.0001, # lb
                         upper=0.9999, # ub
                         b=b, # parameters that we pass into our profit optimization function
                         p=p,
                         K=K,
                         c=c,
                         r=r,
                         V=V,
                         bvec=bvec,
                         delta=delta,
                         f=f1[i],
                         licFee=licFee,
                         tax=tax,
                         method="L-BFGS-B") #WE NEED TO USE THIS METHOD. Effectively, "opts = "
          e1[i]= FishOut$par # have to get our parameter values ($solution)
          Vnew[i]= -FishOut$value # have to get our objective values $objective
          # don't set b == 0 because then shit goes bad; bvec[1] = 10% biomass
          
         # NGGrimes recommends using optim instead of nloptr; a simple maximization problem where we don't need to worry about constraints

#      browser()
          
      } #Close bvec loop
      
      # Spline should be somewhere around here(?)
      
      diffF= sum(abs(f1-oldf1))
      diffE= sum(abs(e1-olde1))
      
    }# Close while loop
    

  return(list(Policy=f1,Enforcement=e1,Value=V,fIUU=f1IUU))
  
} #Close function