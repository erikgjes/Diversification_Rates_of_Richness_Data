#' @title Taphonomy Model
#' @params x Variant frequency matrix generated using the model_cultChange() function
#' @params r Sampling fraction to be retained from objects surviving taphonomic loss. Default is 1.
#' @params lambda Magnitude of taphonomic loss. Default is 0.1
#' @details The model removes a propotion 1-exp(-lambda * t), where t is the number of generations till present (i.e. the last generation in the simulation is assumed to be the present, and hence not affected by taphonomic loss). 
#' @examples 
#' step_nb=50
#' r=0.3
#' N1=1000
#' pop_size = rep(N1,step_nb+1)
#' mutation_rate = rep(0.01,,step_nb+1)
#' b = rep(0,step_nb+1)
#' fullInfo=1
#' cb = 0.3
#' set.seed(123)
#' x = model_cultChange(b,r, pop=initial_pop, step_nb, mutation_rate, pop_size, fullInfo, N1, cb)
#' par(mfrow=c(1,2))
#' plot(0,0,xlim=c(0,ncol(x)),ylim=c(1,nrow(x)),xlab='N.Generations',ylab='Variants',axes=FALSE,type='n',main='Full Sample')
#' for(i in 1:nrow(x)){lines(which(x[i,]>0),rep(i,length(which(x[i,]>0))))}
#' axis(1)
#' x2 = taphonomic_loss(x,lambda=0.1)
#' plot(0,0,xlim=c(0,ncol(x2)),ylim=c(1,nrow(x2)),xlab='N.Generations',ylab='Variants',axes=FALSE,type='n', main='Taphonomic Loss with Lambda=0.1')
#' for(i in 1:nrow(x2)){lines(which(x2[i,]>0),rep(i,length(which(x2[i,]>0))))}
#' axis(1)

taphonomic_loss=function(x,r=1,lambda=0.1)
{
  ngen = ncol(x)
  k = nrow(x)
  N = apply(x,2,sum)[1]
  
  # Internal utility function
  instances<-function(x,variants)    
  {
    x=c(x,variants)
    res<-table(x)-1
    return(res)    
  }
  
  
  
  varmat = matrix(NA,nrow=ngen,ncol=N)
  for (i in 1:nrow(varmat)){varmat[i,] = rep(1:nrow(x),x[,i])}

  recovered_sample_size = round(N*exp(-lambda*(ngen-1):0)*r)
  
  for (i in 1:ngen)
  {
    index = sample(1:N,size=N-recovered_sample_size[i])
    varmat[i,index]=NA
  }
  
  output = matrix(NA,nrow=k,ncol=ngen)
  
  for (i in 1:ngen)
  {
    output[,i] = as.numeric(instances(varmat[i,],1:k))
  }
  
  
  return(output)
}










