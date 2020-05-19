calculate_b<-function(switch_cult,b_value,step_nb)
{
  switch_nb=length(switch_cult)
  b=NULL
  for(i in 1:switch_nb){
    b[switch_cult[i]:(step_nb+1)]=b_value[i]*rep(1,step_nb+1-switch_cult[i]+1)
  }
  return(b)  
}

calculate_pop_size<-function(N1,N2,step_nb){
  pop_size=NULL
{
  if(N1==N2){
    pop_size=N1*rep(1,step_nb+1)
  }
  else{
    #Linear interpolation to find the size at each time step
    pop_size[1]=N1
    pop_size[step_nb+1]=N2
    t=c(1,step_nb+1)
    lm=lm(c(pop_size[1],pop_size[step_nb+1])~t)
    pop_size=round(predict(lm, data.frame(t=seq(from = 1, to = step_nb+1, by = 1))))
  }
  return(pop_size)
}
}

calculate_mutation_rate<-function(switch_mut,mut_value,step_nb)
{
  switch_nb=length(switch_mut)
  mutation_rate=NULL
  for(i in 1:switch_nb){
    mutation_rate[switch_mut[i]:(step_nb+1)]=mut_value[i]*rep(1,step_nb+1-switch_mut[i]+1)
  }
  return(mutation_rate)  
}

calculate_initial_pop<-function(N1)
  {
    a1=round(N1/4)
    a2=round(N1/4)
    a3=round(N1/4)
    a4=N1-(a1+a2+a3)
    initial_pop=c(a1,a2,a3,a4)
    return(initial_pop)
}

r_to_pyrate<-function(x,y) {
  comb=NULL
  for (i in 1:ncol(x)) {
    cv.1<-data.frame(a=x[,i],b=seq(from=1,to=y,by=1),c=x[,i]/sum(x[,i]))
    cv.2<-filter(cv.1,a!=0)
    cv.3<-arrange(cv.2,b)
    cv.i.pyrate<-data.frame(clade=1, species=i, ts=y-first(cv.3$b), te=y-last(cv.3$b),trait=mean(cv.3$c),stringsAsFactors = FALSE)
    comb<-rbind(comb,cv.i.pyrate)
  }
  comb_final<<-comb
}