rm(list = ls())

######################################################################################################################
################### Constructing neighbor matrix of the United States, excluding Alaska and Hawaii ###################
######################################################################################################################
path=getwd()
source(paste(path,'/neighbor.R',sep=""))
NEIGHBOR=neighbor()

###############################################################
################### First simulation method ###################
###############################################################
source(paste(path,'/sim1.R',sep=""))

################################################################
################### Second simulation method ###################
################################################################
source(paste(path,'/sim2.R',sep=""))

###############################################################################
################### Simulates data under the first scenario ################### 
###############################################################################
data.simu = US_simulation1(NEIGHBOR,sd=0.2,M=200)

####################################################################
################### Employ our modeling approach ###################
####################################################################
library(RcppArmadillo)
library(Rcpp)
sourceCpp(paste(path,"/Ccode4_tau.cpp",sep=""))

################### Initialize parameters and hyperparameters settings###################
data = data.simu[[1]]
nrow=dim(data)[1]; ncol=dim(data)[2]; n=dim(data)[3];
C=array(NA,dim=c(nrow,ncol,n))#Pairing labels
for(i in 1:nrow){
  for(j in 1:ncol){
    C[i,j,]=1:n
  }
}
IS0 = array(1,dim=c(nrow,ncol,n))#Binary variable selection
mu0=0; sigma0=5;#Base measure
an =1:n; bn=rep(1,n)#Beta-GOS hyperparam
d=-1; e=0 #MRF hyperparam
ALPHA=2.001; BETA=0.04;#Invers-Gamma hyperparam
dimension = dim(C)

B=10000#Number of MCMC iterations
bi=1000#Number of iterations discarded as burn-in

t1 = Sys.time()
simulation2=Gibbs2D(data,C,IS0,B,an,bn,mu0=mu0,sigma0=sigma0,alpha0=ALPHA,beta0=BETA,d,e,dimension,NEIGHBOR)
runtime = Sys.time()-t1
print(c(runtime,"minutes"))

#######################################################
################### Post processing ###################
#######################################################

posmean2 = rowMeans(simulation2[[3]][,bi:B])
Pos.mean.2 = array(posmean2,dim=dimension); 

#....... PPI .......#
PPI_vec = rowMeans(simulation2[[2]][,bi:B])
PPI_mat = array(PPI_vec,dim=dimension);

##Save file
#filename = paste(path,"/","simulation1d",d,"e",e,"alp",ALPHA,"bet",BETA,".RData",sep="")
#save.image(filename)


###########################################################################################
################### Bayesian FDR, sensitivity, specificity and accuracy ###################
###########################################################################################
grid01.2 = data.simu[[2]]; grid01.2[grid01.2!=0] = 1 #True values of variable selection
alpha = 0.05

source(paste(path,'/BAyesianFDR.R',sep=""))

threshold = apply(PPI_mat[1,,], 1, BayFDR, alpha) # vector of threshold for each state
M.2 = array(0, dim = dimension)

for (i in 1:48){
  M.2[1,i,][PPI_mat[1,i,]>threshold[i]] = 1
  M.2[1,i,][PPI_mat[1,i,]<threshold[i]] = 0
}

sen3 = sum(M.2[grid01.2!=0])/sum(grid01.2)
spe3 = sum(M.2[grid01.2==0]==0)/sum(grid01.2==0)
acc3 = (sum(M.2[grid01.2!=0])+sum(M.2[grid01.2==0]==0))/prod(dim(grid01.2))
out=c(sen3,spe3,acc3)
names(out)<-c("sen","spe","acc")
