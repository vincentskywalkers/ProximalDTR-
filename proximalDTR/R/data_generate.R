data_generate <- function(M,N,ini_M){
  
  X_1 <- function(d,s_1,s_2){
    val = 0.15*100 +  0.9*s_1 - 0.1*s_2 +
      0*I(d==0) + 
      -4*I(d==1)+
      -6*I(d==2)+
      6*I(d==3)+
      4*I(d==4)+ 
      10*I(d==5)-
      8*I(d==6)-
      12*I(d==7)-
      14*I(d==8)-
      2*I(d==9)-
      4*I(d==10)+
      2*I(d==11)+
      rnorm(1,0,2)
    
    return(val)
  }
  
  
  X_2 <- function(d,s_1,s_2){
    
    val = 0.15*12 + 0.9*s_2 + 0.01*s_1 +
      0*I(d==0) + 
      2*I(d==1)+
      -2*I(d==2)+
      3*I(d==3)-
      1*I(d==4)+ 
      1*I(d==5)-
      0*I(d==6) + 
      2*I(d==7)+
      -2*I(d==8)+
      3*I(d==9)-
      1*I(d==10)+ 
      1*I(d==11)-
      rnorm(1,0,1)
    
    return(val)
  }
  
  X_3 <- function(d,s_3){
    
    val =  0.1*70 + 0.9*s_3 +
      0*I(d==0)- 
      4*I(d==1)+
      8*I(d==2)+
      0*I(d==3)+
      12*I(d==4)+ 
      4*I(d==5)-
      0*I(d==6)- 
      4*I(d==7)+
      8*I(d==8)+
      0*I(d==9)+
      12*I(d==10)+ 
      4*I(d==11)-
      rnorm(1,0,1)
    
    return(val)
  }
  
  
  #Utility function
  Reward_func <- function(s_1,s_2,s_3){
    
    val = 3*I(70 <= s_1 & s_1 <= 120) - 3*I(70 > s_1 | s_1 > 150) - 
      I(120 < s_1 & s_1 <= 150) + 
      2*I(5 <= s_2 & s_2 <= 23) - I(s_2 <5 | s_2 > 23) + 
      2*I(66 <= s_3 & s_3 <= 80 ) - I(s_3 < 66 | s_3 > 80)
    
    return(val)
  }
  
  
  X_01 <- matrix(rnorm(N,100,5),N,1)
  X_02 <- matrix(rnorm(N,12,2),N,1)
  X_03 <- matrix(rnorm(N,72,2),N,1)
  X_0 = cbind(X_01,X_02,X_03)
  
  A=matrix(sample(0:11, size = N*M,replace = T,prob = rep(1/12,12)),N,M)
  d=matrix(0,N,M);
  X=array(NA,dim=c(N,M+1,3))
  
  
  for (j in 1:N){
    X[j,1,1] = X_1(A[j,1],X_01[j],X_02[j])
    X[j,1,2] = X_2(A[j,1],X_01[j],X_02[j])
    X[j,1,3] = X_3(A[j,1],X_03[j])
  }
  
  
  for (i in 2:(M+1)){
    for (j in 1:N){
      X[j,i,1] = X_1(A[j,i-1],X[j,i-1,1],X[j,i-1,2])
      X[j,i,2] = X_2(A[j,i-1],X[j,i-1,1],X[j,i-1,2])
      X[j,i,3] = X_3(A[j,i-1],X[j,i-1,3])
    }
  }
  
  R=matrix(0,N,(M));
  for (i in 1:(M)){
    for (j in 1:N){
      R[j,i] = Reward_func(X[j,i,1],X[j,i,2],X[j,i,3])
    }
  }
  
  
  
  X_vec = matrix(NA,(M+1)*N,3)
  for (i in 1:(M+1)){
    for (j in 1:3){
      X_vec[c(1:N)+(i-1)*N,j] = X[,i,j]
    }
  }
  
  A_vec = matrix(NA,(M)*N,1)
  for (i in 1:(M)){
    A_vec[c(1:N)+(i-1)*N,1] = A[,i]
  }
  
  R_vec = matrix(NA,(M)*N,1)
  for (i in 1:(M)){
    R_vec[c(1:N)+(i-1)*N,1] = R[,i]
  }
  
  
  Xmat = matrix(NA,N*M,3)
  Xmat = X_vec[1:(N*(M-1)),]
  Xmat = rbind(X_0,Xmat)
  
  Xmat_prime = matrix(NA,N*M,3)
  Xmat_prime = X_vec[1:(N*(M)),]
  
  train_data = data.frame(Xmat,A_vec,R_vec,Xmat_prime)
  train_data = train_data[((ini_M)*100+1):dim(train_data)[1],]
  names(train_data) = c("X1","X2","X3","A","R","Xprime_1","Xprime_2","Xprime_3")
  
  return(train_data)
  
}


################################################################

data_generate_nopair <- function(stage,N,ini_M){
  
  M = stage-1
  
  X_1 <- function(d,s_1,s_2){
    val = 0.15*100 +  0.9*s_1 - 0.1*s_2 +
      0*I(d==0) + 
      -4*I(d==1)+
      -6*I(d==2)+
      6*I(d==3)+
      4*I(d==4)+ 
      10*I(d==5)-
      8*I(d==6)-
      12*I(d==7)-
      14*I(d==8)-
      2*I(d==9)-
      4*I(d==10)+
      2*I(d==11)+
      rnorm(1,0,2)
    
    return(val)
  }
  
  
  X_2 <- function(d,s_1,s_2){
    
    val = 0.15*12 + 0.9*s_2 + 0.01*s_1 +
      0*I(d==0) + 
      2*I(d==1)+
      -2*I(d==2)+
      3*I(d==3)-
      1*I(d==4)+ 
      1*I(d==5)-
      0*I(d==6) + 
      2*I(d==7)+
      -2*I(d==8)+
      3*I(d==9)-
      1*I(d==10)+ 
      1*I(d==11)-
      rnorm(1,0,1)
    
    return(val)
  }
  
  X_3 <- function(d,s_3){
    
    val =  0.1*70 + 0.9*s_3 +
      0*I(d==0)- 
      4*I(d==1)+
      8*I(d==2)+
      0*I(d==3)+
      12*I(d==4)+ 
      4*I(d==5)-
      0*I(d==6)- 
      4*I(d==7)+
      8*I(d==8)+
      0*I(d==9)+
      12*I(d==10)+ 
      4*I(d==11)-
      rnorm(1,0,1)
    
    return(val)
  }
  
  
  #Utility function
  Reward_func <- function(s_1,s_2,s_3){
    
    val = 3*I(70 <= s_1 & s_1 <= 120) - 3*I(70 > s_1 | s_1 > 150) - 
      I(120 < s_1 & s_1 <= 150) + 
      2*I(5 <= s_2 & s_2 <= 23) - I(s_2 <5 | s_2 > 23) + 
      2*I(66 <= s_3 & s_3 <= 80 ) - I(s_3 < 66 | s_3 > 80)
    
    return(val)
  }
  
  
  X_01 <- matrix(rnorm(N,100,5),N,1)
  X_02 <- matrix(rnorm(N,12,2),N,1)
  X_03 <- matrix(rnorm(N,72,2),N,1)
  X_0 = cbind(X_01,X_02,X_03)
  
  A=matrix(sample(0:11, size = N*M,replace = T,prob = rep(1/12,12)),N,M)
  d=matrix(0,N,M);
  X=array(NA,dim=c(N,M+1,3))
  
  
  for (j in 1:N){
    X[j,1,1] = X_1(A[j,1],X_01[j],X_02[j])
    X[j,1,2] = X_2(A[j,1],X_01[j],X_02[j])
    X[j,1,3] = X_3(A[j,1],X_03[j])
  }
  
  
  for (i in 2:(M+1)){
    for (j in 1:N){
      X[j,i,1] = X_1(A[j,i-1],X[j,i-1,1],X[j,i-1,2])
      X[j,i,2] = X_2(A[j,i-1],X[j,i-1,1],X[j,i-1,2])
      X[j,i,3] = X_3(A[j,i-1],X[j,i-1,3])
    }
  }
  
  R=matrix(0,N,(M));
  for (i in 1:(M)){
    for (j in 1:N){
      R[j,i] = Reward_func(X[j,i,1],X[j,i,2],X[j,i,3])
    }
  }
  
  
  
  X_vec = matrix(NA,(M+1)*N,3)
  for (i in 1:(M+1)){
    for (j in 1:3){
      X_vec[c(1:N)+(i-1)*N,j] = X[,i,j]
    }
  }
  
  A_vec = matrix(NA,(M)*N,1)
  for (i in 1:(M)){
    A_vec[c(1:N)+(i-1)*N,1] = A[,i]
  }
  
  R_vec = matrix(NA,(M)*N,1)
  for (i in 1:(M)){
    R_vec[c(1:N)+(i-1)*N,1] = R[,i]
  }
  
  
  Xmat = matrix(NA,N*M,3)
  Xmat = X_vec[1:(N*(M-1)),]
  Xmat = rbind(X_0,Xmat)
  
  Xmat_prime = matrix(NA,N*M,3)
  Xmat_prime = X_vec[1:(N*(M)),]
  
  
  Xmat = rbind(Xmat,Xmat_prime[(N*(M-1)+1):(N*M),])
  
  return(list(X =  Xmat, A= A_vec, R = R_vec ))
  
}


###############################################################
test.reward <- function(train,beta,seed.num,test.M,lambda,x1_max,x1_min,x2_max,x2_min,x3_max,x3_min,node1,node2,node3,node4,node5,node6,A_type){
  
  set.seed(seed.num) 
  
  X_01 <- matrix(rnorm(N,100,5),N,1)
  X_02 <- matrix(rnorm(N,12,2),N,1)
  X_03 <- matrix(rnorm(N,72,2),N,1)
  X_0 = cbind(X_01,X_02,X_03)
  
  M = test.M
  
  est.A = matrix(NA,N,M)
  
  betahat.tmp = matrix(beta,1)
  
  prob_list_ini = matrix(NA,N,12)
  
  for (j in 1:N){
    prob_a1 = c()
    for (a in A_type){
      prob_a1 =  c(prob_a1,min(1,Pi_policy(a,self_norm(X_0[j,],x1_max,x1_min,x2_max,x2_min,x3_max,x3_min),
                                           A_type,betahat.tmp,lambda,node1,node2,node3,node4,node5,node6)))
    }
    prob_list_ini[j,] = prob_a1 
    if (sum(prob_a1) == 0){
      est.A[j,1] = sample(0:11, size = 1,prob = rep(1/12,12))
    }else{
      est.A[j,1] = sample(0:11, size = 1,prob = prob_a1)
    }
  }
  
  est.X=array(NA,dim=c(N,M,3))
  
  for (j in 1:N){
    est.X[j,1,1] = X_1(est.A[j,1],X_01[j],X_02[j])
    est.X[j,1,2] = X_2(est.A[j,1],X_01[j],X_02[j])
    est.X[j,1,3] = X_3(est.A[j,1],X_03[j])
  }
  
  prob_list_all = matrix(NA,N*(M-1),12)
  jwj = 0
  for (i in 2:(M)){
    for (k in 1:N){
      jwj = jwj + 1
      prob_a1 = c()
      for (a in A_type){
        prob_a1 =  c(prob_a1,min(1,Pi_policy(a,self_norm(est.X[k,i-1,],x1_max,x1_min,x2_max,x2_min,x3_max,x3_min),A_type,betahat.tmp,lambda,node1,node2,node3,node4,node5,node6)))
      }
      prob_list_all[jwj,] = prob_a1
      if (sum(prob_a1) == 0){
        est.A[k,i] = sample(0:11, size = 1,prob = rep(1/12,12))
      }else{
        est.A[k,i] = sample(0:11, size = 1,prob = prob_a1)
      }
      est.X[k,i,1] = X_1(est.A[k,i],est.X[k,i-1,1],est.X[k,i-1,2])
      est.X[k,i,2] = X_2(est.A[k,i],est.X[k,i-1,1],est.X[k,i-1,2])
      est.X[k,i,3] = X_3(est.A[k,i],est.X[k,i-1,3])
      
    }
  }
  
  est.R=matrix(0,N,(M));
  for (i in 1:(M)){
    for (j in 1:N){
      est.R[j,i] = Reward_func(est.X[j,i,1],est.X[j,i,2],est.X[j,i,3])
    }
  }
  
  
  est.R_vec = matrix(NA,(M)*N,1)
  for (i in 1:(M)){
    est.R_vec[c(1:N)+(i-1)*N,1] = est.R[,i]
  }
  
  
  datainfo = list(outcome.tmp = mean(est.R_vec), outcome.sd =sd(est.R_vec),
                  A.output = est.A, prob_list_ini = prob_list_ini, prob_list_all = prob_list_all)
  
  return(datainfo)
  
}


normalize <- function(x) { 
  x <- as.matrix(x)
  minAttr=apply(x, 2, min)
  maxAttr=apply(x, 2, max)
  x <- sweep(x, 2, minAttr, FUN="-") 
  x=sweep(x, 2,  maxAttr-minAttr, "/") 
  attr(x, 'normalized:min') = minAttr
  attr(x, 'normalized:max') = maxAttr
  return (x)
} 

self_norm <- function(x,x1_max,x1_min,x2_max,x2_min,x3_max,x3_min){
  
  self_norm.value = (x - c(x1_min,x2_min,x3_min))/(c(x1_max,x2_max,x3_max) - c(x1_min,x2_min,x3_min))
  return(self_norm.value)
  
}

silverman <- function(d, n)
{
  return((4/(d+2))^(1/(d+4))*n^(-1/(d+4)))
}


# 
# X_1 <- function(d,s_1,s_2){
#   val = 0.15*100 +  0.9*s_1 - 0.1*s_2 +
#     0*I(d==0) + 
#     -4*I(d==1)+
#     -6*I(d==2)+
#     6*I(d==3)+
#     4*I(d==4)+ 
#     10*I(d==5)-
#     8*I(d==6)-
#     12*I(d==7)-
#     14*I(d==8)-
#     2*I(d==9)-
#     4*I(d==10)+
#     2*I(d==11)+
#     rnorm(1,0,2)
#   
#   return(val)
# }
# 
# X_2 <- function(d,s_1,s_2){
#   
#   val = 0.15*12 + 0.9*s_2 + 0.01*s_1 +
#     0*I(d==0) + 
#     2*I(d==1)+
#     -2*I(d==2)+
#     3*I(d==3)-
#     1*I(d==4)+ 
#     1*I(d==5)-
#     0*I(d==6) + 
#     2*I(d==7)+
#     -2*I(d==8)+
#     3*I(d==9)-
#     1*I(d==10)+ 
#     1*I(d==11)-
#     rnorm(1,0,1)
#   
#   return(val)
# }
# 
# 
# X_3 <- function(d,s_3){
#   
#   val =  0.1*70 + 0.9*s_3 +
#     0*I(d==0)- 
#     4*I(d==1)+
#     8*I(d==2)+
#     0*I(d==3)+
#     12*I(d==4)+ 
#     4*I(d==5)-
#     0*I(d==6)- 
#     4*I(d==7)+
#     8*I(d==8)+
#     0*I(d==9)+
#     12*I(d==10)+ 
#     4*I(d==11)-
#     rnorm(1,0,1)
#   
#   return(val)
# }
# 
# 
# Reward_func <- function(s_1,s_2,s_3){
# 
#   val = 3*I(70 <= s_1 & s_1 <= 120) - 3*I(70 > s_1 | s_1 > 150) - 
#     I(120 < s_1 & s_1 <= 150) + 
#     2*I(5 <= s_2 & s_2 <= 23) - I(s_2 <5 | s_2 > 23) + 
#     2*I(66 <= s_3 & s_3 <= 80 ) - I(s_3 < 66 | s_3 > 80)
#   
#   return(val)
# }
# 
