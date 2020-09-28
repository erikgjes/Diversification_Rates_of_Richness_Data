
######### Assigning random letters sequences to varaints #############
random_names <- function(n = 5000) {
  a <- do.call(paste0, replicate(7, sample(LETTERS, n, TRUE), FALSE))
  return(a)
}

######### Population Growth Scenarios ###############
sig_growth<- function(timesteps, assympote, lower_bound, midpoint, scale) {
  
  midpoint = midpoint
  scale = scale
  assympote = assympote
  x_mid = lower_bound - assympote
  
  est = assympote + ( x_mid / ( 1 + exp( scale * (timesteps - midpoint) ) ) )
  
  return(rev(est))
}

expo_growth <- function(timesteps, upper_bound, lower_bound, scale) {
  scale = scale
  x_mid = lower_bound - upper_bound
  est = upper_bound + (x_mid * (1 - exp(-scale*timesteps)))
  return(est)
}


###### Innovation Rate Scenarios #############

innovation_rate<-function(mu_shifts,mu_values,ts)
{
  shifts_nb=length(mu_shifts)
  innovation_rate=NULL
  for(i in 1:shifts_nb){
    innovation_rate[mu_shifts[i]:(ts+1)]=mu_values[i]*rep(1,ts+1-mu_shifts[i]+1)
  }
  innovation_rate <- innovation_rate[-ts]
  return(innovation_rate)  
}

###### Transmission Rate Scenarios #############

trans_rate<-function(cult_shifts,b_values,ts)
{
  shifts_nb=length(cult_shifts)
  b=NULL
  for(i in 1:shifts_nb){
    b[cult_shifts[i]:(ts+1)]=b_values[i]*rep(1,ts+1-cult_shifts[i]+1)
  }
  b <- b[-ts]
  return(b)  
}


############ Function for creating data based on neutral model parameters ############


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
    variants.cnt=table(res[[i-1]])
    trans.prb= ((variants.cnt/N[i-1])^(1-b[i])) / sum((variants.cnt/N[i-1])^(1-b[i]))
    sample.pool = as.numeric(names(variants.cnt))
    res[[i]] = sample(sample.pool,size=N[i],prob=trans.prb,replace=TRUE)
    
    index = which(runif(N[i])<=mu[i])
    if (length(index)>0)
    {
      res[[i]][index]=(tcounter+1):c(tcounter+length(index))
      tcounter = tcounter + length(index) 
    }
  }
  
  res=sapply(res,instances,cases=unique((unlist(res))))
  rownames(res)<-random_names(nrow(res))
  colnames(res)<-paste0("t_",seq(from=ts,to=1,by=-1))
  return(res)
}


freqBias_replicates<-function(N,mu,b,warm,ts,reps){
  res=list()
  for(i in 1:reps){
    temp_res<-freqBias(N = N,mu = mu,warm = warm, ts=ts, b = b)
    res[[i]]<-temp_res
    rownames(res[[i]])<-random_names(nrow(res[[i]]))
    colnames(res[[i]])<-paste0("t_",seq(from=ts,to=1,by=-1))
  }
  
  return(res)
}


########### Function for Formatting Data Matrix into LiteRate############


literate_rep_function<-function(x,replicates,new_directory_name){
  dir.create(new_directory_name)
  wd <- paste0(getwd(),"/",new_directory_name)
  index <- sample(1:length(x),replicates)
  list <- x[c(index)]
  names(list) <- paste0("r",index)
  x_name <- deparse(substitute(x))
  for (i in 1:length(list)) {
    r1 <- matrix_to_literate(list[[i]])
    rep_nb <- names(list[i])
    write.table(r1,paste0(wd,"/",x_name,"_",rep_nb,".txt"),quote=FALSE,sep="\t",row.names=FALSE)
    cat(paste0("python3 LiteRateForward.py -d ",paste0(wd, "/",x_name,"_",rep_nb,".txt")," -TBP"),
        append=TRUE,file=paste0(wd,"/",x_name,"_commands.txt"),sep="\n")
  }
}


matrix_to_literate<-function(x){
  x_t<-t(x)
  time_step <- seq(from=nrow(x_t),to=1,by=-1)
  x_t<- as.data.frame(cbind(x_t,time_step))
  x_long <- x_t %>% pivot_longer(-time_step,names_to="variants",values_to="count")
  x_long <- uncount(x_long,count)
  x_ts<- x_long %>% group_by(variants) %>% summarise(ts=max(time_step), .groups = 'drop')
  x_te<- x_long %>% group_by(variants) %>% summarise(te=min(time_step), .groups = 'drop')
  x_literate<-data.frame(Species=x_ts$variants,
                         max_age=x_ts$ts,
                         min_age=x_te$te)
  colnames(x_literate)<-c("Species","max_age","min_age")
  return(x_literate)
}

### Time-average and sampled version of functions above

matrix_to_literate_ta_samp<-function(x,base,s_frac){
  x_t<-t(x)
  time_step <- seq(from=nrow(x_t),to=1,by=-1)
  x_t<- as.data.frame(cbind(x_t,time_step))
  x_long <- x_t %>% pivot_longer(-time_step,names_to="variants",values_to="count")
  x_long <- uncount(x_long,count)
  x_ts<- x_long %>% group_by(variants) %>% summarise(ts=max(time_step), .groups = 'drop')
  x_te<- x_long %>% group_by(variants) %>% summarise(te=min(time_step), .groups = 'drop')
  x_literate<-data.frame(Species=x_ts$variants,
                         max_age=x_ts$ts,
                         min_age=x_te$te)
  colnames(x_literate)<-c("Species","max_age","min_age")
  x_literate$max_age <- base*round(x_literate$max_age/base) #rounding to nearest multiple of 5
  x_literate$max_age <- x_literate$max_age / base
  x_literate$min_age <- base*round(x_literate$min_age/base) #rounding to nearest multiple of 5
  x_literate$min_age <- x_literate$min_age / base
  x_literate<-sample_frac(x_literate,s_frac) #Randomly selecting 10% of the cultural variants
  return(x_literate)
}

literate_rep_function_ta_samp<-function(x,replicates,new_directory_name,base,s_frac){
  dir.create(new_directory_name)
  wd <- paste0(getwd(),"/",new_directory_name)
  index <- sample(1:length(x),replicates)
  list <- x[c(index)]
  names(list) <- paste0("r",index)
  x_name <- deparse(substitute(x))
  for (i in 1:length(list)) {
    r1 <- matrix_to_literate_ta_samp(list[[i]],base,s_frac)
    rep_nb <- names(list[i])
    write.table(r1,paste0(wd,"/",x_name,"_",rep_nb,".txt"),quote=FALSE,sep="\t",row.names=FALSE)
    cat(paste0("python3 LiteRateForward.py -d ",paste0(wd, "/",x_name,"_",rep_nb,".txt")," -TBP"),
        append=TRUE,file=paste0(wd,"/",x_name,"_commands.txt"),sep="\n")
  }
}


######## Summarizing the Replicated Empirical Values ##############

empirical_values<-function(x){
  res=matrix(nrow=ncol(x),ncol=9)
  for(j in 1:ncol(x)){
    stan_div <- colSums(x != 0)[j]
    new_vars<- table(x[,j-1]==0 & x[,j]>=1)[2]
    loss_vars<- table(x[,j-1]>=1 & x[,j]==0)[2]
    innov_rate <- table(x[,j-1]==0 & x[,j]>=1)[2] / colSums(x)[j]
    orig_rate <- table(x[,j-1]==0 & x[,j]>=1)[2] / colSums(x != 0)[j]
    ex_rate <- table(x[,j-1]>=1 & x[,j]==0)[2] / colSums(x != 0)[j]
    net_rate <- orig_rate - ex_rate
    longevity <- 1 / ex_rate
    turnover_rate <- orig_rate / ex_rate
    res[j,1:9]<-c(stan_div,new_vars,loss_vars,innov_rate,orig_rate,ex_rate,net_rate,longevity,turnover_rate)
  }
  colnames(res)<-c("Richness","New_Vars","Loss_Vars","Innov_Rate","Orig_Rate","Ex_Rate","Net_Rate","Longevity","Turnover_Rate")
  res<-as.data.frame(res)
  return(res)
}

empirical_summary<-function(x){
list.do<-function(.data,fun,...){do.call(what = fun, args = as.list(.data), ...)}
list.cbind<-function (x) {list.do(x, "cbind")}
temp_list<-list(9)
for (i in 1:9){
  temp_list[[i]]<-lapply(x,function(x)(cbind(x[,i]))) %>% list.cbind() %>% as_tibble()
}
  names(temp_list)<-c("Richness","New_Vars","Lost_Vars","Innov_Rate","Orig_Rate","Ex_Rate","Net_Rate","Longevity","Turnover_Rate")
  return(temp_list)
}


