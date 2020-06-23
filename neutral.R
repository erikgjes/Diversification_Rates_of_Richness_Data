#Utility Function for counting cases
instances <- function(x,cases)
{
  x <- c(x,cases)
  return(table(x) - 1)
}

neutral=function(N,mu,warm=1000,ts=50)
{
  variants = 1:N
  tcounter = N
  for (i in 1:warm)
  {
    variants = sample(variants,size=N,replace=TRUE)
    index = which(runif(N)<=mu)
    if (length(index)>0)
    {
      variants[index]=(tcounter+1):c(tcounter+length(index))
      tcounter = tcounter + length(index) 
    }
  }
  
  res = matrix(NA,nrow=ts,ncol=N)
  res[1,]=variants
  for (i in 2:ts)
  {
    res[i,] = sample(res[i-1,],size=N,replace=TRUE)
    index = which(runif(N)<=mu)
    if (length(index)>0)
    {
      res[i,index]=(tcounter+1):c(tcounter+length(index))
      tcounter = tcounter + length(index) 
    }
  }
  
  res=apply(res,1,instances,cases=unique(as.numeric(res)))
  return(res)
}