#' @title predict_proximalDTR
#' @description Internal prediction function for proximal DTR estimation
#' @keywords internal
#' @param object fitted object
#' @param testx Testing data
#' @param ... ...
#' @return The predicted object
#'

predict_rl <- function(object, testx, ...){
  
  
  X_dim = object$X_dim
  betahat = object$betahat
  train = object$data
  node = object$node
  norm_scale = object$norm_scale
  lambda = object$lambda
  
  
  node1 = node[1]
  node2 = node[2]
  node3 = node[3]
  node4 = node[4]
  node5 = node[5]
  node6 = node[6]
  
  if (norm_scale == T){
    
  self_norm <- function(x,train){
    
  
    X_max = c()
    X_min = c()
    for (i in 1:X_dim){
      
      X_max[i] = max(train[,i])
      X_min[i] = min(train[,i])
      
    }
    
    self_norm.value = (x - X_min)/(X_max -  X_min)
    
    return(self_norm.value)
    
  }
  
  }
  
  
  A_type = sort(unique(train$A))

  prob_a = c()
  for (a in A_type){
    prob_a =  c(prob_a, Pi_policy(a,self_norm(testx,train),A_type,betahat,lambda,
                                  node1,node2,node3,node4,node5,node6))
  }
  
  est.A = sample(0:(length(A_type)-1), size = 1, prob = prob_a)
  
  return(list(recommend.trt = est.A, prob = prob_a))
}
  
  





