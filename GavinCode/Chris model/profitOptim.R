######################################
# Internal Optimization Functions
# To be used with RunDynamicOptimizationEnf.R
# Gives (negative) value function value for value function iteration code--------------------------------------------------
######################################

profitOptimEnf = function(enfEffort,b,p,K,c,r,V,bvec,delta, f, licFee, tax)
{  
  
  bcurrent = b*K/2
  
  QLegal = f*(bcurrent)
  
  #  if (bcurrent<=0) QIUU=0 else QIUU = min(bcurrent, optimize(profitPrivateIUUFunc,x=bcurrent,enfEffort=params[2],c(0,bcurrent),maximum=TRUE)$objective)
  
  if (bcurrent - QLegal<=0) QIUU=0 else QIUU = optim(par=0.5,fn=profitPrivateIUUFunc,x=bcurrent - QLegal,enfEffort=enfEffort,lower=0,upper=bcurrent - QLegal,method="L-BFGS-B",control=list(fnscale=-1))$par
  
  profit = profitSocialFunc (QLegal, QIUU, bcurrent, enfEffort, licFee, tax)
  
  bnext=max(min(bvec),stockGrowthFunc(bcurrent - QIUU - QLegal)/(K/2))
  bnext=min(max(bvec),bnext)
  
  out= approx(bvec,t(V),bnext)
  Vnext= out$y
  
  negout= -(profit + delta*Vnext)
#    browser()
  return(negout)
}

profitOptimTwo = function(params=c(f,enfEffort),b,p,K,c,r,V,bvec,delta,licFee, tax)
{  
    
  bcurrent = b*K/2
  
  QLegal = params[1]*(bcurrent)
  
#  if (bcurrent<=0) QIUU=0 else QIUU = min(bcurrent, optimize(profitPrivateIUUFunc,x=bcurrent,enfEffort=params[2],c(0,bcurrent),maximum=TRUE)$objective)
 
  if (bcurrent - QLegal<=0) QIUU=0 else QIUU = optim(par=0.5,fn=profitPrivateIUUFunc,x=bcurrent - QLegal,enfEffort=params[2],lower=0,upper=bcurrent - QLegal,method="L-BFGS-B",control=list(fnscale=-1))$par
  
  profit = profitSocialFunc (QLegal, QIUU, bcurrent, params[2], licFee, tax)
  
  bnext=max(min(bvec),stockGrowthFunc(bcurrent - QLegal - QIUU)/(K/2))
  bnext=min(max(bvec),bnext)
  
  out= approx(bvec,t(V),bnext)
  Vnext= out$y
  
  negout= -(profit + delta*Vnext)
#  browser()
  return(negout)
}