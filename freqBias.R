#Utility Function for counting cases
instances <- function(x,cases)
{
  x <- c(x,cases)
  return(table(x) - 1)
}

freqBias=function(N,N0=N[1],mu0=mu[1],mu,b0=b[1],b=0,warm=1000,ts=50)
{
  variants = 1:N0
  tcounter = N0
  if (length(mu)==1){mu=rep(mu,ts)}
  if (length(N)==1){N=rep(N,ts)}
  if (length(b)==1){b=rep(b,ts)}
  
  for (i in 1:warm)
  {
    variants.cnt=table(variants)
    trans.prb= ((variants.cnt/N0)^(1-b0)) / sum((variants.cnt/N0)^(1-b0))
    sample.pool = as.numeric(names(variants.cnt))
    variants = sample(sample.pool,size=N0,prob=trans.prb,replace=TRUE)
    
    index = which(runif(N0)<=mu)
    if (length(index)>0)
    {
      variants[index]=(tcounter+1):c(tcounter+length(index))
      tcounter = tcounter + length(index) 
    }
  }
  
  res = vector('list',length=ts)
  res[[1]]=variants
  for (i in 2:ts)
  {
    variants.cnt=table(res[[i-1]])
    trans.prb= ((variants.cnt/N[i-1])^(1-b[i])) / sum((variants.cnt/N[i-1])^(1-b[i]))
    sample.pool = as.numeric(names(variants.cnt))
    variants = sample(sample.pool,size=N0,prob=trans.prb,replace=TRUE)
    
    res[[i]] = sample(sample.pool,size=N[i],prob=trans.prb,replace=TRUE)
    
    index = which(runif(N[i])<=mu[i])
    if (length(index)>0)
    {
      res[[i]][index]=(tcounter+1):c(tcounter+length(index))
      tcounter = tcounter + length(index) 
    }
  }
  
  res=sapply(res,instances,cases=unique((unlist(res))))
  return(res)
}