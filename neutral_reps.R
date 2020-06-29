
#1. Create 1000 replicates from neutral model based on parameter settings
#2. Calculate "macro" measures for each time step (each list will consists of 50 rows, 8 columns)
#3. Plot "macro" measures and the expected value for each of these measures
#4. Calculate the mean for each measure in each time step across all the replicates (across the list), these statistics are really only meaningful when evaluating settings with constant innovation (mu) and constant population sizes  

library(tidyverse)
library(magrittr)
library(rlist)
source("neutral_functions.R")

#1
ts=100
N=seq(from=1000,to=1000,length.out = ts)
mu=seq(from=0.01,to=0.05,length.out = ts)
warm=500
rep=1000

ls<-neutral_replicates(N=N,mu=mu,warm=warm,ts=ts,rep=rep)

#2
ls_sum<-lapply(ls,summary_reps)

#Richness
Richness<-lapply(ls_sum,function(x)(cbind(x[,1]))) %>% list.cbind() %>% 
  t() %>% set_colnames(paste0("t_",seq(ts,1,by=-1))) %>% as_tibble()

New_Vars<-lapply(ls_sum,function(x)(cbind(x[,2]))) %>% list.cbind() %>% 
  t() %>% set_colnames(paste0("t_",seq(ts,1,by=-1))) %>% as_tibble()

Lost_Vars<-lapply(ls_sum,function(x)(cbind(x[,3]))) %>% list.cbind() %>% 
  t() %>% set_colnames(paste0("t_",seq(ts,1,by=-1))) %>% as_tibble()

Innov_Rate<-lapply(ls_sum,function(x)(cbind(x[,4]))) %>% list.cbind() %>% 
  t() %>% set_colnames(paste0("t_",seq(ts,1,by=-1))) %>% as_tibble()

Orig_Rate<-lapply(ls_sum,function(x)(cbind(x[,5]))) %>% list.cbind() %>% 
  t() %>% set_colnames(paste0("t_",seq(ts,1,by=-1))) %>% as_tibble()

Ex_Rate<-lapply(ls_sum,function(x)(cbind(x[,6]))) %>% list.cbind() %>% 
  t() %>% set_colnames(paste0("t_",seq(ts,1,by=-1))) %>% as_tibble()


#3
png("plot_mu_change_0.1_0.5.png",units="in",height=10,width=8,res=300)

par(mfrow=c(3,2))
plot_res(Richness,"goldenrod4")
lines(x=seq(from=ts,to=1,by=-1),y=exp.std.div(N,mu),lwd=3,col="brown",lty=2)

plot_res(Innov_Rate,"seagreen2")
lines(x=seq(from=ts,to=1,by=-1),y=mu,lwd=3,col="springgreen4",lty=2)

plot_res(New_Vars,"dodgerblue")
lines(x=seq(from=ts,to=1,by=-1),y=N*mu,lwd=3,col="navyblue",lty=2)

plot_res(Orig_Rate,"dodgerblue")
lines(x=seq(from=ts,to=1,by=-1),y=N*mu / exp.std.div(N,mu),lwd=3,col="navyblue",lty=2)

plot_res(Lost_Vars,"pink")
lines(x=seq(from=ts,to=1,by=-1),y=N*mu,lwd=3,col="firebrick",lty=2)

plot_res(Ex_Rate,"pink")
dev.off()

#4
ls_means<-lapply(ls_sum,colMeans,na.rm=TRUE)
means_df<-data.frame(matrix(unlist(ls_means), nrow=length(ls_means), byrow=T))
colnames(means_df)<-c("Richness","Mean_New_Vars","Mean_Loss_Var","Mean_Mu","Mean_Orig_Rate","Mean_Ex_Rate")

