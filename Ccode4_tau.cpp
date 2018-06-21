//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

void GibbsW(int* mC,double* man,double* mbn,vec *Wout,int N){//Gibbs sampling for the SS weights
  ivec C(mC,N,false);
  vec an(man,N,false);
  vec bn(mbn,N,false);
  vec W(N);W.zeros();
  double An,Bn;
  for(int n=0;n<(N-1);n++){//Rcpp::Rcout << "\n n: " << n << std::endl;
    An=an(n);
    Bn=bn(n);
    for(int i=(n+1);i<N;i++){
      if(C(i)<(n+1)||C(i)==(i+1)){An+=1;}//Rcpp::Rcout << 1;}else{Rcpp::Rcout << 0;}//Here (n+1) and (i+1)
      if(C(i)==(n+1)){Bn+=1;}//Rcpp::Rcout << 1;}else{Rcpp::Rcout << 0;}//Here (n+1)
    }//Rcpp::Rcout << "An: " << An << std::endl;Rcpp::Rcout << "Bn: " << Bn << std::endl;
    W(n)=Rf_rbeta(An,Bn);
  }
  W(N-1)=Rf_rbeta(an(N-1),bn(N-1));
  *Wout=W;
}

ivec Seq2Cluster(int* mC,int N){//From SS sequence to cluster labels
  ivec C(mC,N,false);
  ivec lb=C;
  int max=0;
  for(int i=0;i<N;i++){
    if(C(i)<(i+1))
      lb(i)=lb(C(i)-1);
    else if(C(i)==(i+1)){
      max+=1;lb(i)=max;
    }
  }
  return lb;
}


vec GibbsX(int* mC,double* mY,int* mIS0,double *tauOut,double tau,double mu0,double sigma0,double alpha0,double beta0,int N){//Gibbs sampling for theta and tau_s
  ivec C(mC,N,false);
  vec Y(mY,N,false);
  ivec IS0(mIS0,N,false);
  ivec lb=Seq2Cluster(mC,N);//Rcpp::Rcout << "lb: " << lb << std::endl;
  
  int L=lb.max();
  vec Xstar(L); Xstar.zeros(); 
  vec Xvec(N); Xvec.zeros();
  vec tauJ(L); tauJ.zeros();
  ivec nJ(L); nJ.zeros();
  double mui,sigmai,tau2,sigma02,alpha,beta,temp;
  tau2=pow(tau,2);
  sigma02=pow(sigma0,2);
  for(int i=0;i<L;i++){
    vec YJ=Y(find(lb==(i+1)));//Rcpp::Rcout << "YJ: " << YJ << std::endl;
    temp=YJ.n_elem;nJ(i)=temp;
    if(sum(IS0(find(lb==(i+1))))!=0){
      mui=(temp*mean(YJ)/tau2+mu0/sigma02)/(1/sigma02+temp/tau2);//Rcpp::Rcout << "mui: " << mui << std::endl;
      sigmai=1/(1/sigma02+temp/tau2);//Rcpp::Rcout << "sigmai: " << sigmai << std::endl;
      Xstar(i)=Rf_rnorm(mui,sigmai);
    }
    alpha=alpha0+temp/2;//Rcpp::Rcout << "alpha: " << alpha << std::endl;
    vec vectemp=YJ-Xstar(i);//Rcpp::Rcout << "vectemp: " << vectemp << std::endl;
    beta=beta0+0.5*dot(vectemp,vectemp);//Rcpp::Rcout << "beta: " << beta << std::endl;
    tauJ(i)=sqrt(1/(Rf_rgamma(alpha,1/beta)));
  }
  if((sum(nJ-1))!=0){tau=sqrt(sum(tauJ%tauJ%(nJ-1))/sum(nJ-1));}
  for(int j=0;j<N;j++){
    Xvec(j) = Xstar(lb(j)-1);
  }
  *tauOut=tau;
  return Xvec;
  
}

ivec updateState(int *mC,int *mIS0,int* mclSizeAcc,ivec *clSizeAccout,double* mYacc,vec *Yaccout,ivec clStarters,double piupdate,double *piupdateout,int n, int i,double mu0,double sigma0,double tau,int N){//Function employed by GibbsC to update the states and other parameters
  ivec C(mC,N,false);
  ivec IS0(mIS0,N,false);
  ivec clSizeAcc(mclSizeAcc,N,false);
  vec Yacc(mYacc,N,false);
  
  double mu02,sigma02,tau2;
  mu02=mu0*mu0;
  tau2=tau*tau;
  sigma02=sigma0*sigma0;
  if(C(n-1)==n){
    clStarters.shed_row(as_scalar(find(clStarters==n,1)));//This works ONLY because one and only one element of clStarters==n
    if(IS0(n-1)!=0){
      piupdate+=mu02/(2*sigma02)+0.5*log(1+clSizeAcc(n-1)*sigma02/tau2)-0.5*pow((mu0/sigma02+Yacc(n-1)/tau2),2)/(1/sigma02+clSizeAcc(n-1)/tau2);
    }
  }
  if(i==n){
    int temp=clStarters.n_elem;
    clStarters.resize(temp+1);
    clStarters(temp)=n;
    piupdate=piupdate-0.5*mu02/sigma02-0.5*log(1+clSizeAcc(n-1)*sigma02/tau2)+0.5*pow((mu0/sigma02+Yacc(n-1)/tau2),2)/(1/sigma02+clSizeAcc(n-1)/tau2);
  }
  int clSizeN=clSizeAcc(n-1);//Rcpp::Rcout << "clSizeN: " << clSizeN << std::endl;
  double YaccN=Yacc(n-1);//Rcpp::Rcout << "YaccN: " << YaccN << std::endl;
  int pivot=i;
  if(pivot!=n && pivot!=(n+1)){
    while(C(pivot-1)<pivot){
      clSizeAcc(pivot-1)+=clSizeN;
      Yacc(pivot-1)+=YaccN;
      pivot=C(pivot-1);
    }
    if(IS0(pivot-1)!=0){
      piupdate+=0.5*log(1+clSizeAcc(pivot-1)*sigma02/tau2)-0.5*pow((mu0/sigma02+Yacc(pivot-1)/tau2),2)/(1/sigma02+clSizeAcc(pivot-1)/tau2);
      clSizeAcc(pivot-1)+=clSizeN;
      Yacc(pivot-1)+=YaccN;
      piupdate+=0.5*pow((mu0/sigma02+Yacc(pivot-1)/tau2),2)/(1/sigma02+clSizeAcc(pivot-1)/tau2)-0.5*log(1+clSizeAcc(pivot-1)*sigma02/tau2);
    }
  }
  pivot=C(n-1);
  if(pivot!=n){
    while(C(pivot-1)<pivot){
      clSizeAcc(pivot-1)-=clSizeN;
      Yacc(pivot-1)-=YaccN;
      pivot=C(pivot-1);
    }
    if(IS0(pivot-1)!=0){
      piupdate+=0.5*log(1+clSizeAcc(pivot-1)*sigma02/tau2)-0.5*pow((mu0/sigma02+Yacc(pivot-1)/tau2),2)/(1/sigma02+clSizeAcc(pivot-1)/tau2);
      clSizeAcc(pivot-1)-=clSizeN;
      Yacc(pivot-1)-=YaccN;
      piupdate+=0.5*pow((mu0/sigma02+Yacc(pivot-1)/tau2),2)/(1/sigma02+clSizeAcc(pivot-1)/tau2)-0.5*log(1+clSizeAcc(pivot-1)*sigma02/tau2);
    }
  }
  *clSizeAccout=clSizeAcc;
  *Yaccout=Yacc;
  *piupdateout=piupdate;
  return clStarters;
}

int samplep(double* mp,int P){//Function employed by GibbsC to identify the position to be updated
  vec p(mp,P,false);
  vec pcumsum=cumsum(p);//Rcpp::Rcout << "pcumsum: " << pcumsum.t();
  double rand=::Rf_runif(0,1);//Rcpp::Rcout << "rand: " << rand << std::endl;
  uvec temp=find(pcumsum<rand);
  int out=temp.n_elem+1;
  return out;
}

void GibbsC(vec W,int* mC,int* mIS0,double* mY,ivec *Cout,ivec *IS0out,double tau,double mu0,double sigma0,int* mneighbor,double d,double e,int N){//Gibbs sampling for the values of C and Gamma
  ivec C(mC,N,false);
  vec Y(mY,N,false);
  ivec IS0(mIS0,N,false);
  ivec neighbor(mneighbor,N,false);
  
  W = log(W);
  double const0 = -mu0*mu0/(2*sigma0*sigma0);
  double const1 = sigma0*sigma0/(tau*tau);
  double const2 = mu0/(sigma0*sigma0);
  double const3 = 1/(sigma0*sigma0);
  double tau2=tau*tau;
  double mu02=mu0*mu0;
  double sigma02=sigma0*sigma0;
  int nextCn;
  
  ivec lb=Seq2Cluster(mC,N);//Rcpp::Rcout << "lb: " << lb << std::endl;
  int clusters=lb.max();//Rcpp::Rcout << "clusters: " << clusters << std::endl;
  ivec clStarters0(clusters);
  for(int cl=0;cl<clusters;cl++){
    clStarters0(cl) = as_scalar(find(lb==(cl+1),1))+1; 
  }//Rcpp::Rcout << "clStarters0: " << clStarters0 << std::endl;
  
  ivec clSizeAcc0 = zeros<ivec>(N);int* mclSizeAcc0;mclSizeAcc0=clSizeAcc0.memptr();
  vec Yacc0=zeros<vec>(N);double* mYacc0;mYacc0=Yacc0.memptr();
  int temp;
  for(int i=(N-1);i>(-1);i--){
    temp=C(i)-1;
    if(temp!=i){
      clSizeAcc0(temp)+=clSizeAcc0(i)+1;
      Yacc0(temp)+=Yacc0(i)+Y(i);
    }else{
      clSizeAcc0(temp)=clSizeAcc0(i);
      Yacc0(temp)=Yacc0(i);
    }
  }
  clSizeAcc0+=1;
  Yacc0+=Y;
  //Rcpp::Rcout << "clSizeAcc0: " << clSizeAcc0 << std::endl;
  //Rcpp::Rcout << "Yacc0: " << Yacc0 << std::endl;
  
  double pi_update0=0;
  for(int cl=0;cl<clusters;cl++){
    temp = clStarters0(cl)-1;
    if(IS0(temp)!=0){
      pi_update0 += const0-0.5*log(1+clSizeAcc0(temp)*const1)+0.5*pow((const2+Yacc0(temp)/tau2),2)/(const3+clSizeAcc0(temp)/tau2);
    }
  }//Rcpp::Rcout << "pi_update0: " << pi_update0 << std::endl;
  
  C(1)=1;
  double  sumW1n=W(0);
  vec omega(N);
  double temp2;
  for(int i=0;i<N;i++){
    temp2=d+e*neighbor(i);
    omega(i)=temp2-log(1+exp(temp2));
  }//Rcpp::Rcout << "omega: " << omega << std::endl;
  double prob1,prob2;
  prob1 = (1-exp(omega(0)))*::Rf_dnorm4(Y(0),0,tau,0);//Rcpp::Rcout << "prob1: " << prob1 << std::endl;
  prob2 = exp(omega(0))*::Rf_dnorm4(Y(0),mu0,sqrt(tau2+sigma02),0);//Rcpp::Rcout << "prob2: " << prob2 << std::endl;
  //Random
  if(::Rf_runif(0,prob1+prob2)<prob1){
    IS0(0)=0;
  }else{
    IS0(0)=1;
  }//Rcpp::Rcout << "IS0: " << IS0 << std::endl;
  //End Random
  vec p(N+1); double* mp;mp=p.memptr();
  ivec clSizeAcc(N);
  vec Yacc(N);
  double piupdate,YaccN;piupdate=0;
  int clSizeN,pivot;
  
  for(int n=1;n<N;n++){
    p.zeros();
    if(n>1){
      p(0)=sumW1n-W(0)+log(1-exp(W(0)));
      if(n>2){
        for(int i=1;i<(n-1);i++){
          p(i)=p(i-1)-W(i)+log(1-exp(W(i)))-log(1-exp(W(i-1)));
        }
      }
    }//}Rcpp::Rcout << "p: \n" << p << std::endl;
    
    p(n-1) = log(1-exp(W(n-1)));
    p(n) = sumW1n;
    vec omega_temp(p.n_elem); omega_temp.fill(omega(n));
    p+=omega_temp;
    p(n+1)=sumW1n+log(1-exp(omega(n)));
    sumW1n += W(n);//}Rcpp::Rcout << "p: \n" << p << std::endl;Rcpp::Rcout << "sumW1n: " << sumW1n << std::endl;
    for(int i=(n+1);i>(-1);i--){
      ivec clSizeAcc(clSizeAcc0);//Rcpp::Rcout << "clSizeAcc: " << clSizeAcc.t() << std::endl;
      vec Yacc(Yacc0);//Rcpp::Rcout << "Yacc: " << Yacc.t() << std::endl;
      ivec clStarters(clStarters0);//Rcpp::Rcout << "clStarters: " << clStarters.t() << std::endl;
      piupdate = pi_update0;//Rcpp::Rcout << "piupdate: " << piupdate << std::endl;
      if(C(n)==(n+1)){
        clStarters.shed_row(as_scalar(find(clStarters==(n+1))));
        if(IS0(n)!=0){
          piupdate+=mu02/(2*sigma02)+0.5*log(1+clSizeAcc(n)*sigma02/tau2)-0.5*pow((mu0/sigma02+Yacc(n)/tau2),2)/(1/sigma02+clSizeAcc(n)/tau2);
        }
      }//Rcpp::Rcout << "clStarters: " << clStarters.t() << std::endl;Rcpp::Rcout << "piupdate: " << piupdate << std::endl;
      if(i==n){
        temp=clStarters.n_elem;
        clStarters.resize(temp+1);
        clStarters(temp)=(n+1);//Here n+1
        piupdate=piupdate-0.5*mu02/sigma02-0.5*log(1+clSizeAcc(n)*sigma02/tau2)+0.5*pow((mu0/sigma02+Yacc(n)/tau2),2)/(1/sigma02+clSizeAcc(n)/tau2);
      }//Rcpp::Rcout << "clStarters: " << clStarters.t() << std::endl;Rcpp::Rcout << "piupdate: " << piupdate << std::endl;
      clSizeN=clSizeAcc(n);//Rcpp::Rcout << "clSizeN: " << clSizeN << std::endl;
      YaccN=Yacc(n);//Rcpp::Rcout << "YaccN: " << YaccN << std::endl;
      
      pivot=i;//Rcpp::Rcout << "i: " << i << std::endl;
      if(pivot!=n && pivot!=(n+1)){//Rcpp::Rcout << "pivot: " << C(pivot) << std::endl;
        while(C(pivot)<(pivot+1)){
          clSizeAcc(pivot)+=clSizeN;
          Yacc(pivot)+=YaccN;
          pivot=C(pivot)-1;
        }//Rcpp::Rcout << "clSizeAcc: " << clSizeAcc.t();
        //Rcpp::Rcout << "Yacc: " << Yacc.t();Rcpp::Rcout << "pivot: " << pivot << std::endl;
        
        if(IS0(pivot)!=0){
          piupdate+=0.5*log(1+clSizeAcc(pivot)*sigma02/tau2)-0.5*pow((mu0/sigma02+Yacc(pivot)/tau2),2)/(1/sigma02+clSizeAcc(pivot)/tau2);
          clSizeAcc(pivot)+=clSizeN;
          Yacc(pivot)+=YaccN;
          piupdate+=0.5*pow((mu0/sigma02+Yacc(pivot)/tau2),2)/(1/sigma02+clSizeAcc(pivot)/tau2)-0.5*log(1+clSizeAcc(pivot)*sigma02/tau2);
        }
      }//Rcpp::Rcout << "piupdate: " << piupdate << std::endl;
      //Rcpp::Rcout << "clSizeAcc: " << clSizeAcc.t();Rcpp::Rcout << "Yacc: " << Yacc.t();
      pivot=C(n)-1;//C(n)-1
      if(pivot!=n){
        while(C(pivot)<(pivot+1)){
          clSizeAcc(pivot)-=clSizeN;
          Yacc(pivot)-=YaccN;
          pivot=C(pivot)-1;
        }//Rcpp::Rcout << "clSizeAcc: " << clSizeAcc.t();
        //Rcpp::Rcout << "Yacc: " << Yacc.t();Rcpp::Rcout << "pivot: " << pivot << std::endl;
        if(IS0(pivot)!=0){
          piupdate+=0.5*log(1+clSizeAcc(pivot)*sigma02/tau2)-0.5*pow((mu0/sigma02+Yacc(pivot)/tau2),2)/(1/sigma02+clSizeAcc(pivot)/tau2);
          clSizeAcc(pivot)-=clSizeN;
          Yacc(pivot)-=YaccN;
          piupdate+=0.5*pow((mu0/sigma02+Yacc(pivot)/tau2),2)/(1/sigma02+clSizeAcc(pivot)/tau2)-0.5*log(1+clSizeAcc(pivot)*sigma02/tau2);
        }//Rcpp::Rcout << "piupdate: " << piupdate << std::endl;
        //Rcpp::Rcout << "clSizeAcc: " << clSizeAcc.t();Rcpp::Rcout << "Yacc: " << Yacc.t();
      }
      p(i) += piupdate;//Rcpp::Rcout << "p: " << p.t();
    }
    for(unsigned int hhh=n+2;hhh<p.n_elem;hhh++){p(hhh)=0;}
    //Rcpp::Rcout << "p: " << p.t();
    //double maxp=max(p.rows(0,n+1));Rcpp::Rcout << "maxp: " << maxp << std::endl;
    //p.rows(0,n+1)=p.rows(0,n+1)-max(p.rows(0,n+1));Rcpp::Rcout << "p: " << p.t();
    p.rows(0,n+1) = exp(p.rows(0,n+1)-max(p.rows(0,n+1)));//Rcpp::Rcout << "p: " << p.t();
    p = p/sum(p);//Rcpp::Rcout << "p: " << p.t();
    //From now on based on a random result
    nextCn = samplep(mp,n+1);//Rcpp::Rcout << "nextCn: " << nextCn << std::endl;
    if(nextCn==(n+2)){
      IS0(n)=0;
    }else if(nextCn==(n+1)){
      IS0(n)=1;
    }else{
      IS0(n)=IS0(nextCn-1);
    }
    
    if(nextCn==(n+2)){
      nextCn = n+1;
    }//Rcpp::Rcout << "n: " << n << std::endl;Rcpp::Rcout << "nextCn: " << nextCn << std::endl;
    //Rcpp::Rcout << "IS0: " << IS0.t();
    clStarters0=updateState(mC,mIS0,mclSizeAcc0,&clSizeAcc0,mYacc0,&Yacc0,clStarters0,piupdate,&pi_update0,n+1,nextCn,mu0,sigma0,tau,N);
    
    C[n] = nextCn;//Rcpp::Rcout << "C: " << C.t();
  }
  
  *Cout=C;
  *IS0out=IS0;
  
}

//[[Rcpp::export]]
Rcpp::List Gibbs2D(vec dataarray,ivec Carray,ivec IS0array,int B, vec an,vec bn,double mu0, double sigma0,double alpha0,double beta0,double d, double e,ivec dim,imat NEIGHBOR){//This is the main
  Rcpp::List out(4);//Rcpp::Rcout << "IS0array: " << IS0array << std::endl;
  int x = dim[0];//Rcpp::Rcout << "x: " << x << std::endl;
  int y = dim[1];//Rcpp::Rcout << "y: " << y << std::endl;
  int t = dim[2];//Rcpp::Rcout << "t: " << t << std::endl;
  int xxxx,yyyy;
  xxxx=(x==0)?1:x;//Rcpp::Rcout << "xxxx: " << xxxx << std::endl;
  yyyy=(y==0)?1:y;//Rcpp::Rcout << "yyyy: " << yyyy << std::endl;
  ivec C(t);ivec IS0(t);vec Y(t);ivec neighbor(t);
  int *mC,*mIS0,*mneighbor;
  mC=C.memptr();mneighbor=neighbor.memptr();
  double *mY;
  mY=Y.memptr();
  mIS0=IS0.memptr();
  double* man;man=an.memptr();
  double* mbn;mbn=bn.memptr();
  int xyt=dataarray.n_elem;//Rcpp::Rcout << "xyt: " << xyt << std::endl;
  unsigned int index;
  imat matC(xyt,B);
  imat matIS0(xyt,B);
  mat matW(xyt,B);
  mat matX(xyt,B);
  mat matTau(x*y,B);
  vec Warray(xyt);Warray.zeros();vec W(t);
  vec Xarray(xyt);Xarray.zeros();
  uvec timeStart(t);uvec indexVec(t);
  
  for(int k=0;k<t;k++){
    timeStart[k]=k*(x*y);
    for(int j=0;j<y;j++){
      for(int i=0;i<x;i++){
        index=k*(x*y)+j*(x)+i;Warray[index]=Rf_rbeta(an[index],bn[index]);//Rcpp::Rcout << "index: " << index << std::endl;
      }
    }
  }//Rcpp::Rcout << "Warray: " << Warray.t();
  //Rcpp::Rcout << "timeStart: " << timeStart.t();
  
  vec tau(x*y);
  for(int j=0;j<y;j++){
    for(int i=0;i<x;i++){
      index=j*x+i;tau[index]=sqrt(1/(Rf_rgamma(alpha0,1/beta0)));//Rcpp::Rcout << "index: " << index << std::endl;
    }
  }
  for(int b=0;b<B;b++){if (b%1000==0 && b!=0){Rcpp::Rcout << "Iteration " << b << " done."<< std::endl;}
    for(int i=0;i<x;i++){
      for(int j=0;j<y;j++){
        index=j*x+i;indexVec=index+timeStart;//Rcpp::Rcout << "indexVec: " << indexVec.t();
        uvec neighborlocation = find(NEIGHBOR.row(j));//Rcpp::Rcout << "neighborlocation: " << neighborlocation << std::endl;
        for(int kkkk=0;kkkk<t;kkkk++){//Rcpp::Rcout << "check: " << IS0array(kkkk*xxxx*yyyy+neighborlocation) << std::endl;
          neighbor[kkkk]=sum(IS0array(kkkk*xxxx*yyyy+neighborlocation));//Rcpp::Rcout << "neighbor: " << neighbor << std::endl;
        }//Rcpp::Rcout << "neighbor: " << neighbor << std::endl;	
        C=Carray(indexVec);//Rcpp::Rcout << "C: " << C.t();
        IS0=IS0array(indexVec);//Rcpp::Rcout << "IS0: " << IS0.t();
        Y=dataarray(indexVec);//Rcpp::Rcout << "Y: " << Y.t();
        W=Warray(indexVec);//Rcpp::Rcout << "W: " << W.t();
        GibbsC(W,mC,mIS0,mY,&C,&IS0,tau[index],mu0,sigma0,mneighbor,d,e,t);//Update C and Gamma
        Carray(indexVec)=C;
        IS0array(indexVec)=IS0;
        vec Xnew;
        Xnew = GibbsX(mC,mY,mIS0,&tau[index],tau[index],mu0,sigma0,alpha0,beta0,t);//Update theta and tau_s
        Xarray(indexVec)=Xnew;
        GibbsW(mC,man,mbn,&W,t);//Update W
        Warray(indexVec)=W;
      }
    }//Rcpp::Rcout << "Carray: " << Carray.t();
    matTau.col(b) = tau;
    matC.col(b)=Carray;
    matIS0.col(b)=IS0array;
    matW.col(b)=Warray;
    matX.col(b)=Xarray;
  }
  
  out[0]=matC;
  out[1]=matIS0;
  out[2]=matX;
  out[3]=matTau;
  return out;
}

