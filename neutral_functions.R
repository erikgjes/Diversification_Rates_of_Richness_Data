

instances <- function(x,cases)
{
  x <- c(x,cases)
  return(table(x) - 1)
}


neutral=function(N,N0=N[1],mu0=mu[1],mu,warm,ts)
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

neutral_replicates<-function(N,mu,warm,ts,reps){
  res=list()
  for(i in 1:reps){
    temp_res<-neutral(N = N,mu = mu,warm = warm, ts=ts)
    res[[i]]<-temp_res
  }
  return(res)
}

summary_reps<-function(x){
  res=matrix(nrow=ncol(x),ncol=6)
  for(j in 1:ncol(x)){
    stan_div <- colSums(x != 0)[j]
    new_vars<- table(x[,j-1]==0 & x[,j]>=1)[2]
    loss_vars<- table(x[,j-1]>=1 & x[,j]==0)[2]
    Innov_rate <- table(x[,j-1]==0 & x[,j]>=1)[2] / colSums(x)[j]
    orig_rate <- table(x[,j-1]==0 & x[,j]>=1)[2] / colSums(x != 0)[j]
    ex_rate <- table(x[,j-1]>=1 & x[,j]==0)[2] / colSums(x != 0)[j]
    res[j,1:6]<-c(stan_div,new_vars,loss_vars,Innov_rate,orig_rate,ex_rate)
  }
  colnames(res)<-c("Richness","New_Vars","Loss_Vars","Innov_Rate","Orig_Rate","Ex_Rate")
  return(res)
}

#Function for calculating the expexted richness based on model parameters
exp.std.div = function(N,mu,n=N)
{
  theta=2*N*mu
  EK = 0
  for (i in 0:(n-1))
  {
    EK = EK + theta/(theta+i)
  }
  return(EK)
}


plot_res<-function(x,color){
  plot(x=c(1:ncol(x)),y=x[ts/2,],type="n",xlim=c(ts,0),
       xlab="Time",ylab=paste0(deparse(substitute(x))),
       ylim=c(min(x,na.rm=TRUE)*0.8,max(x,na.rm=TRUE)*1.1))
  apply(x,MAR=1,lines,col="grey80",lwd=0.5,x=seq(from=ts,to=1,by=-1))
  lines(x=seq(from=ts,to=1,by=-1),y=colMeans(x),col=color,lwd=3)
  
}
