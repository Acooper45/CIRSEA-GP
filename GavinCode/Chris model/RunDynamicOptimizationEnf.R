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
                  
#        FishOutF= optim(par=guess,fn=profitOptimFish,lower=0.0001,upper=0.9999,b=b,p=p,K=K,c=c,r=r,V=V,bvec=bvec,delta=delta,enfEffort=guessEnf, fIUU=f1IUU[i],licFee=licFee, tax=tax,method="L-BFGS-B")
        
#        f1[i]= FishOutF$par
                  
#        FishOutE= optim(par=guessEnf,fn=profitOptimEnf,lower=0.0001,upper=0.9999,b=b,p=p,K=K,c=c,r=r,V=V,bvec=bvec,delta=delta,f=f1[i], fIUU=f1IUU[i],licFee=licFee, tax=tax,method="L-BFGS-B")
          
#        e1[i]= FishOutE$par
  
#        e1[i]=0.5
  
#browser()
        if (pol==1)
        {
          f1 =  rep(r/2,length(bvec))
          FishOut= optim(par=guessEnf,fn=profitOptimEnf,lower=0.0001,upper=0.9999,b=b,p=p,K=K,c=c,r=r,V=V,bvec=bvec,delta=delta,f=f1[i],licFee=licFee, tax=tax,method="L-BFGS-B")
          e1[i]= FishOut$par
          Vnew[i]= -FishOut$value
          if (b==0) f1IUU[i]=0 else f1IUU[i] = max(0,optim(par=0.5,fn=profitPrivateIUUFunc,lower=0,upper=b*K/2-f1[i]*b*K/2,x=b*K/2-f1[i]*b*K/2,enfEffort=e1[i],method="L-BFGS-B",control=list(fnscale=-1))$par/(b*K/2-f1[i]*b*K/2))
        }

        if (pol==2)
        {
          f1 =  c(rep(0,floor(length(bvec)/2)),rep(r/2,ceiling(length(bvec)/2)))
          FishOut= optim(par=guessEnf,fn=profitOptimEnf,lower=0.0001,upper=0.9999,b=b,p=p,K=K,c=c,r=r,V=V,bvec=bvec,delta=delta,f=f1[i],licFee=licFee, tax=tax,method="L-BFGS-B")
          e1[i]= FishOut$par
          Vnew[i]= -FishOut$value
          if (b==0) f1IUU[i]=0 else f1IUU[i] =  max(0,optim(par=0.5,fn=profitPrivateIUUFunc,lower=0,upper=b*K/2-f1[i]*b*K/2,x=b*K/2-f1[i]*b*K/2,enfEffort=e1[i],method="L-BFGS-B",control=list(fnscale=-1))$par/(b*K/2-f1[i]*b*K/2))
#        browser()
        }

        if (pol==3)
        {
          FishOut= optim(par=c(guess,guessEnf),fn=profitOptimTwo,lower=0.0001,upper=0.9999,b=b,p=p,K=K,c=c,r=r,V=V,bvec=bvec,delta=delta,licFee=licFee, tax=tax,method="L-BFGS-B")
          f1[i]= FishOut$par[1]
          e1[i]= FishOut$par[2]
          Vnew[i]= -FishOut$value
          if (b==0) f1IUU[i]=0 else f1IUU[i] =  max(0,optim(par=0.5,fn=profitPrivateIUUFunc,lower=0,upper=b*K/2-f1[i]*b*K/2,x=b*K/2-f1[i]*b*K/2,enfEffort=e1[i],method="L-BFGS-B",control=list(fnscale=-1))$par/(b*K/2-f1[i]*b*K/2))
#          browser()
        }

#      browser()

      } #Close bvec loop
    
      diffF= sum(abs(f1-oldf1))
      diffE= sum(abs(e1-olde1))
      
    }# Close while loop
    

  return(list(Policy=f1,Enforcement=e1,Value=V,fIUU=f1IUU))
  
} #Close function