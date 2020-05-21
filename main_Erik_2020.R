source("cultChange_erik.R")
source("inputData_erik.R")
source("pyrate_utilities.r")
library(tidyverse)
#number of analysis time steps
step_nb=50
#fraction of the population that is replaced in each year
r=0.3

#---------------------------------------------------------------------------------------------------------------
#Population size
#---------------------------------------------------------------------------------------------------------------
#population size at the beginning (N1) and end (N2) of consideration
N1=1000
N2=1000
b_value=c(0) #(-0.03,-0.01,0,0.01,0.03)
mut_value=0.01 #or 0.1

#population size for each time step (linear interpolation between N1 and N2)
pop_size=calculate_pop_size(N1,N2,step_nb)

#Initial condition (meaningless for analysis as we allow for burn-in time)
a1=round(N1/4)
a2=round(N1/4)
a3=round(N1/4)
a4=N1-(a1+a2+a3)
initial_pop=c(a1,a2,a3,a4)
#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#Parameter for cultural transmission process
#---------------------------------------------------------------------------------------------------------------
#Strength of frequency-dependent selection in each time step
#switch_cult contains the time steps when a switch occurs and b_value the corresponding selection strength
#IMPORTANT: switch_cult has to contain 1 as this defines the initial b
#IMPORTANT: for constant selection strength choose switch_cult=c(1) and define the value in b_value
switch_cult=c(1)
#switch_cult=c(1) #constant
#b_value=c(-0.03)
#b_value=0.03 #adjust accordingly (-0.03,-0.01,0.01,0.03)
b=calculate_b(switch_cult,b_value,step_nb)
#commoness threshold
cb=0.3
#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#Mutation
#---------------------------------------------------------------------------------------------------------------
#Mutation rate in each time step
#switch_mut contains the time steps when a switch occurs and mut_value the corresponding mutation rate
#IMPORTANT: switch_mut has to contain 1 as this defines the initial mutation_rate
#IMPORTANT: for a constant mutation process choose switch_mut=c(1) and define the value in mut_value
#switch_mut=c(1,20)
switch_mut=c(1)
#mut_value=c(0.01,0.1)
#mut_value=0.01 #or 0.01
mutation_rate=calculate_mutation_rate(switch_mut,mut_value,step_nb)
#---------------------------------------------------------------------------------------------------------------

#fullInfo=1: every mutation is recorded separately
#fullInfo=2: all mutation are summarised in type k+1
fullInfo=1

#Output: a matrix where the columns correspond to the time steps and the rows row to the different variant types; each entry gives the frequency 
#of a specific variant at time t
tempData=model_cultChange(b,r, pop=initial_pop, step_nb, mutation_rate, pop_size, 
                          fullInfo, N1, cb)


################
#Taking the output of the matrix above and wrangling it into pyrate format
#Clade is a generic variable, if additional information is known about each variant this can be adjusted to reflect this information
#species is the cultural variant number
#ts is the time of variant first apperance (50 is oldest as there are 50 time steps)
#te is the time of variant last apperance (0 is most recent)
#trait is the average frequency compared to all other variants over the life of the variant, this can really be an discrete value that is meaningful

#Currently data manipulation only takes into account the first and last apperance. It does not currently take into account if a variant 
#emerges goes extinct and then re-emerges. Currently, that scenario would be marked by the very first and very last apperance and not by 
#apperance in between these (although I will fix this in the near future)
#library(dplyr)


####Sampling the data for PyRate
#tempData <- as.data.frame(tempData)
#tempData <- tempData %>% mutate(Status=ifelse(.[[,step_nb]]!=0,"Extinct","Extant"))
random_names <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  return(a)
}
sample_names<-random_names(nrow(tempData))
rownames(tempData)<-sample_names

status <- data.frame(variants=rownames(tempData),
                     Status=ifelse(tempData[,step_nb]==0,"Extinct","Extant"))
#tempData <- cbind(tempData,status)

tempData_t<-t(tempData) #transposing the data matrix

#tempData_t<-as.data.frame(tempData_t)

#sample_names<-random_names(ncol(tempData_t))
#colnames(tempData_t)<-sample_names
time_step <- seq(from=nrow(tempData_t),to=1,by=-1)
tempData_t<- cbind(tempData_t,time_step)
#rownames(tempData_t)<-tempData_t[,ncol(tempData_t)]

tempData_t <- as_tibble(tempData_t)

#tempData_t %>% add_row(ifelse(tempData_t$time_step==1,"Extant","Extinct"))

t2 <- tempData_t %>% pivot_longer(-time_step,names_to="variants",values_to="count")
t3 <- uncount(t2,count)
t4 <- left_join(t3,status,by="variants")

sample_data_pyrate <- data.frame(Species=t4$variants,
                          Status=t4$Status,
                          min_age=t4$time_step-0.5,
                          max_age=t4$time_step+0.5)

sample_data_literate <- sample_data_pyrate %>% group_by(Species) %>% 
  summarise(min_age=min(min_age),max_age=max(max_age))

write.table(sample_data_pyrate,"A3_5.txt",
            quote=FALSE,row.names = FALSE,sep="\t")

write.table(sample_data_literate,"A3_5_Lite.txt",
            quote=FALSE,row.names = FALSE,sep="\t")

##Only Necessary For PyRate
extract.ages("A3_5.txt",replicates=10)


