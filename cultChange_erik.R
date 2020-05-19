library("gtools")
library("Matrix")
library("parallel")

model_cultChange<-function(b,r, pop=initial_pop, step_nb, mutation_rate, pop_size, fullInfo,N1,cb)
  
{ 
  #------------------------------------------------------------------------------------------------------------
  # Burn-in period to achieve stationary instantaneous distribution  
  #------------------------------------------------------------------------------------------------------------
  for(t in 1:500){

    relFreq=pop/sum(pop)
    type_nb=length(pop)
    
    #Total number of variants removed (u) and added (v) in each time step (kept constant as we want to generate a stationary state)
    u=ceiling(r*N1);
    v=u;
    
    # Remove u variants randomly, probability that a variant of type i is removed is given by its relative frequency
    l=sample(1:type_nb,u,replace=TRUE,prob=relFreq) 
    index=1
    while(index<(u+1)){
      if(pop[l[index]]>=1){
        pop[l[index]]=pop[l[index]]-1
        index=index+1
      }
      else{
        bool=TRUE
        while(bool){
          a=sample(1:type_nb,1,prob=relFreq)
          if(pop[a]>0){
            pop[a]=pop[a]-1
            bool=FALSE
            index=index+1
          }
        }
      }
    }
    x=pop/sum(pop)
    
    #Determination of the probabilities that a variant of type j is added
    copyProb_h=rep(0,type_nb)
    copyProb=rep(0,type_nb)
    for(j in 1:type_nb){
      
      copyProb_h[j]=x[j]+b[1]*x[j]-b[1]*cb
      if (copyProb_h[j]<0 || x[j]==0 ){
        copyProb_h[j]=0
      }
    }
    #normalisation
    copyProb=copyProb_h/sum(copyProb_h)
    
    # Adding u variants according to the considered process of cultural transmission
    l=sample(1:type_nb,v,replace = T,prob=copyProb)
    
    for(i in 1:v){
      #Mutation
      if(runif(1)<mutation_rate[1]){
        pop=c(pop,1)
      }
      else{
        pop[l[i]]=pop[l[i]]+1
      }
    }
    pop=pop[pop != 0]
  } #End of burn in time
  #------------------------------------------------------------------------------------------------------------


  #------------------------------------------------------------------------------------------------------------
  #Analysis period
  #------------------------------------------------------------------------------------------------------------

  #Determination of the total number of variants removed (u) and added (v) in each time step 
  #This allows for temporal changes in population size and therefore in u and v
  u=rep(0,step_nb)
  v=rep(0,step_nb) 
  for(i in 1:step_nb){
    if(pop_size[i]<=pop_size[i+1]){
      u[i]=round(r*pop_size[i])
      v[i]=pop_size[i+1]-(pop_size[i]-u[i])
    } else{
      u[i]=round(r*pop_size[i+1]+pop_size[i]-pop_size[i+1])
      v[i]=pop_size[i+1]-(pop_size[i]-u[i])
    }
  }
  
  
  if(fullInfo == 2){
    pop=c(pop,0)
    type_nb=length(pop)
  }
  
  for(t in 1:step_nb){
    
    relFreq=pop/sum(pop)
    if(fullInfo == 1){type_nb=length(pop)
    }
    
    # Remove u variants randomly, probability that a variant of type i is removed is given by its relative frequency
    l=sample(1:type_nb,u[t],replace=TRUE,prob=relFreq) 
    index=1
    while(index<(u[t]+1)){
      if(pop[l[index]]>=1){
        pop[l[index]]=pop[l[index]]-1
        index=index+1
      }
      else{
        bool=TRUE
        while(bool){
          a=sample(1:type_nb,1,prob=relFreq)
          if(pop[a]>0){
            pop[a]=pop[a]-1
            bool=FALSE
            index=index+1
          }
        }
      }
    }
    x=pop/sum(pop)
    
    #Determination of the probabilities that a variant of type j is added
    copyProb_h=rep(0,type_nb)
    copyProb=rep(0,type_nb)
    for(j in 1:type_nb){
      
      copyProb_h[j]=x[j]+b[t]*x[j]-b[t]*cb
      
      if (copyProb_h[j]<0 || x[j]==0 ){
        copyProb_h[j]=0
      }
    }
    #normalisation
    copyProb=copyProb_h/sum(copyProb_h)
    
    # Adding u variants according to the considered process of cultural transmission
    l=sample(1:type_nb,v[t],replace = T,prob=copyProb)
    
    for(i in 1:v[t]){
      #Mutation
      if(runif(1)<mutation_rate[t]){
        if(fullInfo == 1){pop=c(pop,1)
        }
        else{pop[type_nb]=pop[type_nb]+1
        }
      }
      else{
        pop[l[i]]=pop[l[i]]+1
      }
    }

    #saving each time step in the matrix data
    data_h=matrix(0,length(pop),t)
    if(t>1){
      data_h[1:data_n,1:t-1]=data
      data_n=length(pop)
      data_h[1:data_n,t]=t(pop)
    }
    else{
      data_n=length(pop)
      data_h[1:data_n,t]=t(pop)
    }
    data=data_h
    
  } #End of time loop
  
  return(data)

}

