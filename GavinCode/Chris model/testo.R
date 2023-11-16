fTest=vector()
for (j in 1:length(bvec))
{
  b=bvec[j]
#  browser()
  if (b==0) fTest[j]=0 else fTest[j] = optim(par=0.5,fn=profitPrivateIUUFunc,lower=0,upper=b*K/2-optimal$Policy[j]*b*K/2,x=b*K/2-optimal$Policy[j]*b*K/2,enfEffort=optimal$Enforcement[j],method="L-BFGS-B",control=list(fnscale=-1))$par/(b*K/2-optimal$Policy[j]*b*K/2)
}
fTest