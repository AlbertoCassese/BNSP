US_simulation2 <- function(NEIGHBOR,sd,M){
  states = state.abb[c(-2,-11)]
  region1 = c("WA","OR","ID","MT","WY")
  region2 = c("CA","NV","UT","AZ","CO","NM")
  region3 = c("ND","SD","NE","KS","MN","IA","MO","WI","IL","MI","IN","OH")
  region4 = c("OK","TX","AR","LA","KY","TN","MS","AL")
  region5 = c("WV","VA","DE","MD","NC","SC","GA","FL")
  region6 = c("PA","NY","NJ","CT","RI","MA","VT","NH","ME")

  #region = state.division[-c(2,11)] #excluding Alaska and Hawaii
  n = length(states)
  regions = rep(0,n)
  regions[match(region1,states)] = 1
  regions[match(region2,states)] = 2
  regions[match(region3,states)] = 3
  regions[match(region4,states)] = 4
  regions[match(region5,states)] = 5
  regions[match(region6,states)] = 6
  
  L=4
  a11=0.8
  A=matrix(c(0.9,1-a11,0.1,a11),nr=2)
  pi=rep(1,L)/L ; #probability allocating among L components
  mu0=0
  mu=rep(NA,L)
  mu[1]=-1
  mu[2]=-0.5
  mu[3]=0.5
  mu[4]=1
  SD=sd
  
  #replications=500
  state1=rep(NA,M)
  state2=rep(NA,M)
  state3=rep(NA,M)
  state4=rep(NA,M)
  state5=rep(NA,M)
  state6=rep(NA,M)
 
  mu1 = mu2 = mu3 = mu4 = mu5 = mu6 = rep(0,M)
  STATES=seq(0,L,1)
  state1[1]=sample(STATES,1,prob=c(1,rep(0,L)))
  state2[1]=sample(STATES,1,prob=c(1,rep(0,L)))
  state3[1]=sample(STATES,1,prob=c(1,rep(0,L)))
  state4[1]=sample(STATES,1,prob=c(1,rep(0,L)))
  state5[1]=sample(STATES,1,prob=c(1,rep(0,L)))
  state6[1]=sample(STATES,1,prob=c(1,rep(0,L)))

  
  mu1[1]=0; mu2[1]=0; mu3[1]=0; mu4[1]=0
  mu5[1]=0; mu6[1]=0; 
  
  for(m in 2:M){
    if(state1[m-1]==0){temp=A[1,]}else{temp=A[2,]}
    state1[m]=sample(c(0,1),1,prob=temp)
    #if(state[m]==0){y[m]=rnorm(1,0,1)}else{y[m]=rnorm(1,mu1[state[m]],1)}
    if(state2[m-1]==0){temp=A[1,]}else{temp=A[2,]}
    state2[m]=sample(c(0,1),1,prob=temp)
    if(state3[m-1]==0){temp=A[1,]}else{temp=A[2,]}
    state3[m]=sample(c(0,1),1,prob=temp)
    if(state4[m-1]==0){temp=A[1,]}else{temp=A[2,]}
    state4[m]=sample(c(0,1),1,prob=temp)
    if(state5[m-1]==0){temp=A[1,]}else{temp=A[2,]}
    state5[m]=sample(c(0,1),1,prob=temp)
    #if(state[m]==0){y[m]=rnorm(1,0,1)}else{y[m]=rnorm(1,mu1[state[m]],1)}
    if(state6[m-1]==0){temp=A[1,]}else{temp=A[2,]}
    state6[m]=sample(c(0,1),1,prob=temp)
 
  }
  
  ones1 = sum(state1!=0)
  ones2 = sum(state2!=0)
  ones3 = sum(state3!=0)
  ones4 = sum(state4!=0)
  ones5 = sum(state5!=0)
  ones6 = sum(state6!=0)

  
  
  state1[state1!=0] = sample.int(L,size=ones1,replace=TRUE,prob=pi) # the four states
  state2[state2!=0] = sample.int(L,size=ones2,replace=TRUE,prob=pi)
  state3[state3!=0] = sample.int(L,size=ones3,replace=TRUE,prob=pi)
  state4[state4!=0] = sample.int(L,size=ones4,replace=TRUE,prob=pi)
  state5[state5!=0] = sample.int(L,size=ones5,replace=TRUE,prob=pi) # the four states
  state6[state6!=0] = sample.int(L,size=ones6,replace=TRUE,prob=pi)
 
  mu1[state1!=0] = mu[state1] # give y the corresponding mean value
  mu2[state2!=0] = mu[state2]
  mu3[state3!=0] = mu[state3]
  mu4[state4!=0] = mu[state4]
  mu5[state5!=0] = mu[state5] # give y the corresponding mean value
  mu6[state6!=0] = mu[state6]
 
  
  US.cluster = matrix(NA,n,M)
  US.series = matrix(0,n,M)
  
  US.cluster[regions==1,] = rep(state1,each=length(which(regions==1)))
  US.cluster[regions==2,] = rep(state2,each=length(which(regions==2)))
  US.cluster[regions==3,] = rep(state3,each=length(which(regions==3)))
  US.cluster[regions==4,] = rep(state4,each=length(which(regions==4)))
  US.cluster[regions==5,] = rep(state5,each=length(which(regions==5)))
  US.cluster[regions==6,] = rep(state6,each=length(which(regions==6)))
  
  US.series[regions==1,] = rep(mu1,each=length(which(regions==1)))
  US.series[regions==2,] = rep(mu2,each=length(which(regions==2)))
  US.series[regions==3,] = rep(mu3,each=length(which(regions==3)))
  US.series[regions==4,] = rep(mu4,each=length(which(regions==4)))
  US.series[regions==5,] = rep(mu5,each=length(which(regions==5)))
  US.series[regions==6,] = rep(mu6,each=length(which(regions==6)))
  
  
  US.series = US.series + rnorm(n*M,mean=0,sd=SD)
  
  US.cluster2 = US.series2 = array(0,dim=c(1,n,M))
  US.cluster2[1,,]=US.cluster; US.series2[1,,]=US.series
  
  return(list(US.series2,US.cluster2))
}
