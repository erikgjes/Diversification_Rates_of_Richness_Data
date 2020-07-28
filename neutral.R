#Utility Function for counting cases
instances <- function(x,cases)
{
  x <- c(x,cases)
  return(table(x) - 1)
}

neutral=function(N,N0=N[1],mu0=mu[1],mu,warm=1000,ts=50)
{
  variants = 1:N0
  tcounter = N0
  if (length(mu)==1){mu=rep(mu,ts)}
  if (length(N)==1){N=rep(N,ts)}
  
  for (i in 1:warm)
  {
    variants = sample(variants,size=N0,replace=TRUE)
    index = which(runif(N0)<=mu0)
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
    res[[i]] = sample(res[[i-1]],size=N[i],replace=TRUE)
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
