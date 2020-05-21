#' @title Rare Variant Model
#' @params x Variant frequency matrix generated using the model_cultChange() function
#' @params cutoff Cutoff value below which the variant is eliminated
#' @params perGeneration Logical value. When set to TRUE observations are eliminated per generation (i.e. cells with counts below cutoff within each generation are converted to 0), when set to FALSE entire variants with a total frequency below cutoff are eliminated. Default is FALSE.  


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
#' par(mfrow=c(1,3))
#' plot(0,0,xlim=c(0,ncol(x)),ylim=c(1,nrow(x)),xlab='N.Generations',ylab='Variants',axes=FALSE,type='n',main='Full Sample')
#' for(i in 1:nrow(x)){lines(which(x[i,]>0),rep(i,length(which(x[i,]>0))))}
#' axis(1)
#' x2 = rare_variant(x,cutoff=5,perGeneration=FALSE)
#' plot(0,0,xlim=c(0,ncol(x2)),ylim=c(1,nrow(x2)),xlab='N.Generations',ylab='Variants',axes=FALSE,type='n', main='b')
#' for(i in 1:nrow(x2)){lines(which(x2[i,]>0),rep(i,length(which(x2[i,]>0))))}
#' axis(1)
#' x3 = rare_variant(x,cutoff=5,perGeneration=TRUE)
#' plot(0,0,xlim=c(0,ncol(x3)),ylim=c(1,nrow(x3)),xlab='N.Generations',ylab='Variants',axes=FALSE,type='n', main='b')
#' for(i in 1:nrow(x3)){lines(which(x3[i,]>0),rep(i,length(which(x3[i,]>0))))}
#' axis(1)

rare_variant=function(x,cutoff,perGeneration=FALSE)
{
  ngen = ncol(x)
  k = nrow(x)
  N = apply(x,2,sum)[1]

  if (perGeneration){x[which(x<cutoff,arr.ind = TRUE)]=0}
  if (!perGeneration)
  {
    i = which(apply(x,1,sum)<cutoff)
    x = x[-i,]
  }
  return(x)
}



