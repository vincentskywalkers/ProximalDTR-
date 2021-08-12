#' @title proximalDTR 
#' @name proximalDTR
#' @description The proximal Bellman embedding model for estiamting dynamic treatment regime
#' @param X The longitudinal covariates matrix
#' @param A A is the longitudinal assigned treatment  
#' @param R R is the longitudinal assigned treatment 
#' @param stage The number of stages
#' @param sgd.fold.num The batch size folder for stochastic gradient descent
#' @param lambda.set The tunning candidate set for sparisty and smoothness parameter
#' @param step.beta The learning rate for main parameter 
#' @param step.omega The learning rate for Lagrange parameter 
#' @param desc.rate The decreasing rate for the learning rate
#' @param max.iter The maximum iterations number for main algorithm
#' @param max.iter.cv The maximum iterations number for corss validation algorithm of lambda
#' @param bw The kernel bandwidth for the U-statistic type estimator
#' @param maxitr Maximum number of iterations
#' @param stop.threshold The stoppting criterion for algorithm, which measure the L^2 distance of iterative beta
#' @param cv Should the cross validation be implemented
#' @param norm_scale Should the normalization of state variables be implemented
#' @param node The nodes of B-spline
#' @param trace Should the information be displayed
#' @return A object; a list consisting of
#' \item{betahat}{The optimal \code{beta} value}
#' \item{obj.value}{The final objective value}
#' \item{beta.diff}{The L2 distances of iterative beta}
#' \item{para.value}{The final lower bound of parametric value}
#' \item{norm_scale}{Should the normalization of state variables be implemented}
#' \item{lambda}{The final cross validation selected lambda}
#' \item{X_dim}{The dimension of state space}

#' @references Zhou, W., Zhu, R., Qu, A., "A Proximal Approach for Estimating Dynamic treatment Regime in Infinite Horizon".
#' 


proximalDTR <- function(X, A, R, n_ID, stage, sgd.fold.num = 4, lambda.set = c(0.1,0.5,1),
                        step.beta = 0.01, step.omega = 0.01, desc.rate = 0.005,
                        max.iter = 3000, max.iter.cv = 1000, bw = 0.5, gamma = 0.9,
                        stop.threshold = 0.01, cv=FALSE, trace =T,
                        norm_scale = T, node = c(0.15,0.3,0.4,0.5,0.8,0.9)){
  
  ################################################################
  ####################---Initialization---########################
  ################################################################
  
  X_dim = ncol(X)
  
  train = data.frame(X[1:((stage-1)*n_ID),],A,R,X[(n_ID+1):((stage)*n_ID),])
  colnames(train) = c(paste0("X",1:X_dim),"A","R",  paste0("Xprime_",1:X_dim))
  
  train.org = train 
  
  n_train = dim(train)[1]
  p_train = dim(train)[2]
  M = stage - 1
  A_type = sort(unique(train$A))
  
  node1 = node[1]
  node2 = node[2]
  node3 = node[3]
  node4 = node[4]
  node5 = node[5]
  node6 = node[6]
  
  if (norm_scale == T){
    
    for (i in 1:X_dim){
      
      train[,i] = as.numeric(normalize(train[,i]))
      train[,p_train-i+1] = as.numeric(normalize(train[,p_train-i+1]))
      
    }
    
  }
  
  ################################################################
  ################---Kernel Matrix---#############################
  ################################################################
  
  train_mat = as.matrix(train)
  
  train_list = list()
  kernel_list = list()
  
  bw = bw
  for (kk in 1:N){
    train_list[[kk]] = as.matrix(train_mat[seq(kk,N*M,N),])
    kernel_list[[kk]]  = KernelDist_single(as.matrix(train_list[[kk]][,1:X_dim]), train_list[[kk]][,(X_dim+1)],0,bw,bw)
  }
  
  omega = matrix(rep(0,6),1)
  
  
  ####################################################################
  #########################  Tuning lambda  ##########################
  #################################################################### 
  
  value.lambda <- c()
  
  if (cv != FALSE){
    m = 0
    for (lambda in lambda.set){
      m = m + 1
      obj.tmp = c()
      beta_tmp.list = matrix(NA,200,(X_dim*9+1)*length(A_type))
      for (k in 1:200){
        
        if(trace == T){
          print(paste("set.order: ", m, "search: ",k))
        }
        set.seed(250*k) 
        beta = runif((X_dim*9+1)*length(A_type),-10,10)
        beta_tmp.list[k,] = beta
        
        total.obj.tmp = c()
        for (kj in 1:N){
          total.obj.tmp = c(total.obj.tmp,kernel_obj(train_list[[kj]],matrix(beta,1),omega,A_type,gamma,lambda,kernel_list[[kj]],node1,node2,node3,node4,node5,node6,X_dim))
        }
        obj.tmp[k] = mean(total.obj.tmp)
      }
      
      index1 = order(obj.tmp)[1]
      beta = beta_tmp.list[index1,]
      
      ####################################################################
      #########################  Cross.Validate SGD ######################
      ####################################################################
      
      iters = max.iter.cv
      
      value.func.value <- c()
      obj.value <- c()
      grad.value <- c()
      reward.tmp <- c()
      reward.sd <- c()
      beta_diff <- c()
      
      eta_b0 <- step.beta
      eta_omg <- step.omega
      rate_des <- desc.rate
      
      foldnum = sgd.fold.num
      indexfold = createFolds(c(1:N),k=foldnum)
      
      betalist <- matrix(NA,iters,length(beta))
      
      jj=0
      for(j in 1:iters){
        
        eta_b <- eta_b0/(1+rate_des*sqrt(j))
        
        betalist[j,] <- beta
        subset <- unlist(indexfold[sample(1:foldnum,1)])
        
        beta = matrix(beta,1)
        
        grad.loss = loss_grad(train_list[[subset[1]]],beta,omega,A_type,eta_b,eta_omg,gamma,lambda,kernel_list[[subset[1]]],node1,node2,node3,node4,node5,node6,X_dim)
        grad.beta = grad.loss$G.beta/length(subset)     
        grad.omega = grad.loss$G.omega/length(subset) 
        
        for (kjj in 2:length(subset)){
          
          grad.tmp = loss_grad(train_list[[subset[kjj]]],beta,omega,A_type,eta_b,eta_omg,gamma,lambda,kernel_list[[subset[kjj]]],node1,node2,node3,node4,node5,node6,X_dim)
          grad.beta = grad.beta + grad.tmp$G.beta/length(subset)      
          grad.omega = grad.omega + grad.tmp$G.omega/length(subset)   
          
        }
        
        beta = beta - eta_b*grad.beta
        omega = omega - eta_omg*grad.omega
        
        
        if (j %% 100 == 0 | j==1){
          jj = jj + 1
          beta_diff[jj] = sqrt(sum((betalist[j,] - betalist[j-1,])^2))
          
          final.obj.tmp = c()
          value.total = c()
          
          for (kj in 1:N){
            value.total = c(value.total, V_func(train_mat[kj,1:X_dim],A_type,matrix(beta,1),lambda,node1,node2,node3,node4,node5,node6))
            
            final.obj.tmp = c(final.obj.tmp,kernel_obj(train_list[[kj]],matrix(beta,1),omega,A_type,gamma,lambda,kernel_list[[kj]],node1,node2,node3,node4,node5,node6,X_dim))
          }
          obj.value[jj] = mean(final.obj.tmp)
          grad.value[jj] = sum(grad.beta^2)
          value.func.value[jj] = mean(  value.total)
          
          if(length(beta_diff) > 2 & beta_diff[jj] <= 0.01){
            break
          }
        }
        
        if(trace == T){
          print(paste("set.order: ",m,"iter: ",j))
        }
        
      }
      
      
      lower_bd_value <- c()
      lower_bd_value = value.func.value[length(value.func.value)] - (1/(1-gamma))*0.5*lambda
      
      value.lambda = c(value.lambda, lower_bd_value )
      
    }
    
    lambda = lambda.set[which.max(value.lambda)]
  }else{
    lambda = lambda.set[1]
  }
  
  ####################################################################
  ##########################  Main SGD ###############################
  ####################################################################
  
  obj.tmp = c()
  beta_tmp.list = matrix(NA,200,(X_dim*9+1)*length(A_type))
  for (k in 1:200){
    if(trace == T){
      print(paste("main search: ",k))
    }
    set.seed(250*k) 
    beta = runif((X_dim*9+1)*length(A_type),-10,10)
    beta_tmp.list[k,] = beta
    
    total.obj.tmp = c()
    for (kj in 1:N){
      total.obj.tmp = c(total.obj.tmp,kernel_obj(train_list[[kj]],matrix(beta,1),omega,A_type,gamma,lambda,kernel_list[[kj]],node1,node2,node3,node4,node5,node6,X_dim))
    }
    obj.tmp[k] = mean(total.obj.tmp)
  }
  
  index1 = order(obj.tmp)[1]
  beta = beta_tmp.list[index1,]
  
  ################################################################
  ####################---SGD Iteration---#########################
  ################################################################
  
  iters = max.iter
  
  obj.value <- c()
  grad.value <- c()
  reward.tmp <- c()
  reward.sd <- c()
  beta_diff <- c()
  value.func.value <- c()
  
  eta_b0 <- step.beta
  eta_omg <- step.omega
  rate_des <- desc.rate
  
  foldnum = sgd.fold.num
  indexfold = createFolds(c(1:N),k=foldnum)
  
  betalist <- matrix(NA,iters,length(beta))
  
  jj=0
  for(j in 1:iters){
    
    eta_b <- eta_b0/(1+rate_des*sqrt(j))
    
    betalist[j,] <- beta
    subset <- unlist(indexfold[sample(1:foldnum,1)])
    
    beta = matrix(beta,1)
    
    
    grad.loss = loss_grad(train_list[[subset[1]]],beta,omega,A_type,eta_b,eta_omg,gamma,lambda,kernel_list[[subset[1]]],node1,node2,node3,node4,node5,node6,X_dim)
    grad.beta = grad.loss$G.beta/length(subset)     
    grad.omega = grad.loss$G.omega/length(subset) 
    
    for (kjj in 2:length(subset)){
      
      grad.tmp = loss_grad(train_list[[subset[kjj]]],beta,omega,A_type,eta_b,eta_omg,gamma,lambda,kernel_list[[subset[kjj]]],node1,node2,node3,node4,node5,node6,X_dim)
      grad.beta = grad.beta + grad.tmp$G.beta/length(subset)      
      grad.omega = grad.omega + grad.tmp$G.omega/length(subset)   
      
    }
    
    beta = beta - eta_b*grad.beta
    omega = omega - eta_omg*grad.omega
    
  
    if (j %% 50 == 0 | j==1){
      jj = jj + 1
      beta_diff[jj] = sqrt(sum((betalist[j,] - betalist[j-1,])^2))
      #  reward.tmp[jj] = test.reward(train, beta,rr+5-100,36,lambda,X1_max,X1_min,X2_max,X2_min,X3_max,X3_min,node1,node2,node3,node4,node5,node6,A_type)$outcome.tmp
      
      final.obj.tmp = c()
      value.total = c()
      
      for (kj in 1:N){
        value.total = c(value.total, V_func(train_mat[kj,1:X_dim],A_type,matrix(beta,1),lambda,node1,node2,node3,node4,node5,node6))
        final.obj.tmp = c(final.obj.tmp,kernel_obj(train_list[[kj]],matrix(beta,1),omega,A_type,gamma,lambda,kernel_list[[kj]],node1,node2,node3,node4,node5,node6,X_dim))
      }
      obj.value[jj] = mean(final.obj.tmp)
      grad.value[jj] = sum(grad.beta^2)
      value.func.value[jj] = mean(value.total)
      
      if(length(beta_diff) > 2 & beta_diff[jj] <= 0.01){
        break
      }
      
    }
    
    if(trace == T){
      print(paste("iter: ",j))
    }
    
  }
  
  lower_bd_value <- c()
  lower_bd_value = value.func.value - (1/(1-gamma))*0.5*lambda
  
  
  return(list(betahat = beta, obj.value = obj.value,
              beta.diff = beta_diff, para.value = lower_bd_value, data=train.org,
              norm_scale = norm_scale, node = node, 
              lambda = lambda, X_dim = X_dim))
  
}

