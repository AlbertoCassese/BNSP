US_simulation1 <- function(NEIGHBOR,sd,M){
  L=4
  a11=0.8
  A=matrix(c(0.9,1-a11,0.1,a11),nr=2)
  pi=rep(1,L)/L ; #probability allocating among L components
  mu0=0
  mu1=rep(NA,L)
  mu1[1]=-1
  mu1[2]=-0.5
  mu1[3]=0.5
  mu1[4]=1
  SD=sd
  
  state=rep(NA,M)
  mu = rep(0,M)
  STATES=seq(0,L,1)
  state[1]=sample(STATES,1,prob=c(1,rep(0,L)))
  mu[1]=0
  
  for(m in 2:M){
    if(state[m-1]==0){temp=A[1,]}else{temp=A[2,]}
    state[m]=sample(c(0,1),1,prob=temp)
  }
  
  ones = sum(state)
  state[state!=0] = sample.int(L,size=ones,replace=TRUE,prob=pi) # the four states
  mu[state!=0] = mu1[state] # give y the corresponding mean value
  
  n = dim(NEIGHBOR)[1] # No. of states
  US.series = US.01 = US.cluster = array(0,dim=c(1,n,M))
  dimnames(US.series) <- list(c(""),dimnames(NEIGHBOR)[[1]])
  US.series[1,"Washington",] = mu
  state.count = 1
  lag = rep(NA,n)
  current.state = c(45) # 45 index of Washington
  lag[45] = 0
  itration = 0
  while(state.count < n){
    neighbor.state = NULL
    temp = length(current.state)
    for(i in 1:temp){
      neighbor.state = c(neighbor.state,which(NEIGHBOR[current.state[i],]!=0))
      # all the neighbors of current.state
    }
    # remove the repated neighbors
    repeated = which(!is.na(lag[neighbor.state]))
    if(any(repeated)){neighbor.state = neighbor.state[-repeated]}
    neighbor.state = unique(neighbor.state)
    
    state.count = state.count+length(neighbor.state)
    itration = itration+1
    lag[neighbor.state] = itration
    current.state = neighbor.state
  }
  
  for(j in 1:n){
    US.cluster[1,j,(lag[j]+1):M] = state[1:(M-lag[j])]
    US.series[1,j,(lag[j]+1):M] = mu[1:(M-lag[j])] 
  }
  
  US.series = US.series + rnorm(prod(dim(US.series)),0,SD)
  return(list(US.series,US.cluster))
}
