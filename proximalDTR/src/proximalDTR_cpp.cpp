#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//    ----------------------------------------------------------------
//    Kernel Matrix Computation 
//    ----------------------------------------------------------------

// [[Rcpp::export]]


arma::mat KernelDist_single(const arma::mat& X, 
                            const arma::vec& A,
                            double diag, 
                            double bw,
                            double bw_A)
{
  int N = X.n_rows;
  arma::mat kernel_matrix(N, N);
  
  for (int i = 0; i < N; i++)
  {
    kernel_matrix(i, i) = diag;
    for (int j = 0; j < i; j ++)
    {
      kernel_matrix(j,i) = exp(-sum(pow(X.row(i)-X.row(j),2))/bw/bw-sum(pow(A(i)-A(j),2))/bw_A/bw_A);
      kernel_matrix(i,j) = kernel_matrix(j,i);
    }
  }
  
  return(kernel_matrix);
}

// [[Rcpp::export]]


arma::mat EpanKernelDist_single(const arma::mat& X, double diag)
{
  int N = X.n_rows;
  arma::mat kernel_matrix(N, N, arma::fill::zeros);
  double u;
  
  for (int i = 0; i < N; i++)
  {
    kernel_matrix(i, i) = diag;
    for (int j = 0; j < i; j ++)
    {
      u = sum(pow(X.row(i)-X.row(j),2));
      
      if (u > -1 && u < 1)
        kernel_matrix(j,i) = pow((1-u*u), 3);
      
      kernel_matrix(i,j) = kernel_matrix(j,i);
    }
  }
  return(kernel_matrix);
}


//    ----------------------------------------------------------------
//    Q Function for (s,a)
//    ----------------------------------------------------------------

double Qfunc(const arma::rowvec& X_tmp,
             const int A_tmp,
             const arma::mat& beta,
             const arma::colvec& A_type,
             double node1,
             double node2,double node3,double node4,double node5,double node6)
{
  
  int X_size = X_tmp.n_elem;
  //int Num_center = center.n_elem;
  int len_basis = X_size*9+1;
  int A_size = A_type.n_elem;
  
  //Rcout << "basis "   << std::endl;
  
  arma::mat Phi(1,len_basis*A_size);
  Phi.fill(0);
  
  for (int i = 0; i < A_size; i++) {
    
    if (A_tmp == A_type(i)){
      
      Phi(0,i*(len_basis)) = 1;
      
      for (int j = 0; j < X_size; j++) {
        
        Phi(0,i*(len_basis)+j*9+1) = X_tmp(j);
        Phi(0,i*(len_basis)+j*9+2) = pow(X_tmp(j),2);
        Phi(0,i*(len_basis)+j*9+3) = pow(X_tmp(j),3);
        Phi(0,i*(len_basis)+j*9+4) = pow(std::max(0.0,X_tmp(j)-node1),  3);
        Phi(0,i*(len_basis)+j*9+5) = pow(std::max(0.0,X_tmp(j)-node2),  3);
        Phi(0,i*(len_basis)+j*9+6) = pow(std::max(0.0,X_tmp(j)-node3),  3);
        Phi(0,i*(len_basis)+j*9+7) = pow(std::max(0.0,X_tmp(j)-node4),  3);
        Phi(0,i*(len_basis)+j*9+8) = pow(std::max(0.0,X_tmp(j)-node5),  3);
        Phi(0,i*(len_basis)+j*9+9) = pow(std::max(0.0,X_tmp(j)-node6),  3);
        
        
      }
      
    }
  }
  
  
  double Q_value = arma::as_scalar(Phi * beta.t());
  
  return(Q_value);
  
}


//    ----------------------------------------------------------------
//    Basis Function for (s,a)
//    ----------------------------------------------------------------



// [[Rcpp::export]]

arma::rowvec Basis_func(const arma::rowvec& X_tmp,
                        const int A_tmp,
                        const arma::colvec& A_type,
                        double node1,
                        double node2,double node3,double node4,double node5,double node6)
{
  
  int X_size = X_tmp.n_elem;
  // int Num_center = center.n_elem;
  int len_basis = X_size*9+1;
  int A_size = A_type.n_elem;
  
  arma::mat Phi(1,len_basis*A_size);
  Phi.fill(0);
  
  for (int i = 0; i < A_size; i++) {
    
    if (A_tmp == A_type(i)){
      
      Phi(0,i*(len_basis)) = 1;
      
      for (int j = 0; j < X_size; j++) {
        
        Phi(0,i*(len_basis)+j*9+1) = X_tmp(j);
        Phi(0,i*(len_basis)+j*9+2) = pow(X_tmp(j),2);
        Phi(0,i*(len_basis)+j*9+3) = pow(X_tmp(j),3);
        Phi(0,i*(len_basis)+j*9+4) = pow(std::max(0.0,X_tmp(j)-node1),  3);
        Phi(0,i*(len_basis)+j*9+5) = pow(std::max(0.0,X_tmp(j)-node2),  3);
        Phi(0,i*(len_basis)+j*9+6) = pow(std::max(0.0,X_tmp(j)-node3),  3);
        Phi(0,i*(len_basis)+j*9+7) = pow(std::max(0.0,X_tmp(j)-node4),  3);
        Phi(0,i*(len_basis)+j*9+8) = pow(std::max(0.0,X_tmp(j)-node5),  3);
        Phi(0,i*(len_basis)+j*9+9) = pow(std::max(0.0,X_tmp(j)-node6),  3);
        
      }
      
    }
  }
  
  arma::rowvec basis = Phi;
  
  //Rcout << "basis " << basis  << std::endl;
  
  return(basis);
  
}


//    ----------------------------------------------------------------
//    Value Function
//    ----------------------------------------------------------------


// [[Rcpp::export]]


double V_func(const arma::rowvec& X_tmp,
              const arma::colvec& A_type,
              const arma::mat& beta,
              double lambda,
              double node1,
              double node2,double node3,double node4,double node5,double node6)
{
  
  int A_size = A_type.n_elem;
  
  arma::vec Qtmp(A_size); 
  
  // Rcout << "train,   F = " <<  1 << std::endl;
  
  for (int i = 0; i < A_size; i++) {
    
    Qtmp(i) = Qfunc(X_tmp,A_type(i),beta,A_type,node1,node2,node3,node4,node5,node6);
    
  }
  
  
  //Rcout << "train,   F = " <<  Qtmp << std::endl;
  //Rcout << "beta " <<  beta << std::endl;
  
  // Rcout << "train,   F = " <<  Qtmp << std::endl;
  
  arma::ucolvec order_index = sort_index(Qtmp,"descend");
  arma::vec Q_sort = sort(Qtmp,"descend");
  
  //Rcout << "index  F = " <<  order_index << std::endl;
  
  // Rcout << "q_sort,   F = " << Q_sort << std::endl;
  
  //Rcout << "train,   F = " <<  2 << std::endl;
  
  double ii = 0;
  while (lambda + (ii+1)*Q_sort(ii) >= sum(Q_sort.head(ii+1))){
    ii = ii + 1;
    if (ii ==  A_size){
      break;
    }
  }
  ii = ii - 1;
  
  // Rcout << "ii   F = " << ii << std::endl;
  
  arma::vec target_A_set = A_type(order_index.head(ii+1));
  
  //  Rcout << "target = " << target_A_set << std::endl;
  
  double G_value = (sum(Qtmp(order_index.head(ii+1)))/lambda-1)/(ii+1);  
  //Rcout << "G vale   F = " << G_value << std::endl;
  
  double V_value = 0.5*lambda*(1+ sum(pow(Qtmp(order_index.head(ii+1))/lambda,2) 
                                        - pow(G_value,2)));
  // Rcout << "v val   F = " << V_value  << std::endl;
  
  return(V_value);
}  

//    ----------------------------------------------------------------
//    Gradient Value Function
//    ----------------------------------------------------------------


// [[Rcpp::export]]

arma::rowvec V_grad(const arma::rowvec& X_tmp,
                    const arma::colvec& A_type,
                    const arma::mat& beta,
                    double lambda,
                    double node1,
                    double node2,double node3,double node4,double node5,double node6,
                    int X_size)
{
  
  
  // int X_size = X_tmp.n_elem;
  int len_basis = X_size*9+1;
  int A_size = A_type.n_elem;
  
  arma::vec Qtmp(A_size); 
  
  
  for (int i = 0; i < A_size; i++) {
    
    Qtmp(i) = Qfunc(X_tmp,A_type(i),beta,A_type,node1,node2,node3,node4,node5,node6);
    
  }
  
  //Rcout << "what??? 1 " << 321 << std::endl;
  
  arma::ucolvec order_index = sort_index(Qtmp,"descend");
  arma::vec Q_sort = sort(Qtmp,"descend");
  
  //Rcout << "what??? 2 " << Qtmp << std::endl;
  
  double ii = 0;
  while (lambda + (ii+1)*Q_sort(ii) >= sum(Q_sort.head(ii+1))){
    ii = ii + 1;
    if (ii ==  A_size){
      break;
    }
  }
  ii = ii - 1;
  //Rcout << "what??? 2 " << ii << std::endl;
  
  //Rcout << "what " << 123 << std::endl;
  
  arma::vec target_A_set = A_type(order_index.head(ii+1));
  
  // Rcout << "target " << target_A_set << std::endl;
  
  arma::rowvec basis_sum(len_basis*A_size);
  basis_sum.fill(0);
  
  ///Rcout << "basis_sum " << basis_sum << std::endl;
  
  for (int j = 0; j < ii+1; j++) {
    basis_sum += Basis_func(X_tmp,target_A_set(j),A_type,node1,node2,node3,node4,node5,node6);
  }
  
  // Rcout << "Initial value,   F = " << basis_sum << std::endl;
  
  arma::rowvec V_grad_value(len_basis*A_size);
  
  V_grad_value.fill(0);
  
  for (int j = 0; j < ii+1; j++) {
    
    V_grad_value += Qfunc(X_tmp,target_A_set(j),beta,A_type,node1,node2,node3,node4,node5,node6)*Basis_func(X_tmp,target_A_set(j),A_type,node1,node2,node3,node4,node5,node6)/lambda - 
      basis_sum*(sum(Qtmp(order_index.head(j+1)))/lambda-1)/pow(ii+1,2);
  }
  
  return(V_grad_value);
} 


//    ----------------------------------------------------------------
//    Sparse Policy Function
//    ----------------------------------------------------------------
// [[Rcpp::export]]



double Pi_policy(const int A_tmp,
                 const arma::rowvec& X_tmp,
                 const arma::colvec& A_type,
                 const arma::mat& beta,
                 double lambda,
                 double node1,
                 double node2,double node3,double node4,double node5,double node6)
{
  
  //int X_size = X_tmp.n_elem;
  //int len_basis = X_size*3+1;
  int A_size = A_type.n_elem;
  
  arma::vec Qtmp(A_size); 
  
  for (int i = 0; i < A_size; i++) {
    
    Qtmp(i) = Qfunc(X_tmp,A_type(i),beta,A_type,node1,node2,node3,node4,node5,node6);
    
  }
  
  //Rcout << "Qtmp   F = " << Qtmp<< std::endl;
  
  
  arma::ucolvec order_index = sort_index(Qtmp,"descend");
  arma::vec Q_sort = sort(Qtmp,"descend");
  
  //Rcout << "33= " << 1<< std::endl;
  
  double ii = 0;
  while (lambda + (ii+1)*Q_sort(ii) >= sum(Q_sort.head(ii+1))){
    ii = ii + 1;
    if (ii ==  A_size){
      break;
    }
  }
  ii = ii - 1;
  
  arma::vec target_A_set = A_type(order_index.head(ii+1));
  
  double G_value = (sum(Qtmp(order_index.head(ii+1)))/lambda-1)/(ii+1); 
  double policy_value = std::max(0.0,Qfunc(X_tmp,A_tmp,beta,A_type,node1,node2,node3,node4,node5,node6)/lambda - G_value);
  
  return(policy_value);
  
}  


//    ----------------------------------------------------------------
//    Gradient Policy Function
//    ----------------------------------------------------------------

// [[Rcpp::export]]


arma::rowvec Pi_grad(const int A_tmp,
                     const arma::rowvec& X_tmp,
                     const arma::colvec& A_type,
                     const arma::mat& beta,
                     double lambda,
                     double node1,
                     double node2,double node3,double node4,double node5,double node6,
                     int X_size)
{
  
  //int X_size = X_tmp.n_elem;
  int len_basis = X_size*9+1;
  int A_size = A_type.n_elem;
  
  
  
  arma::vec Qtmp(A_size); 
  
  for (int i = 0; i < A_size; i++) {
    
    Qtmp(i) = Qfunc(X_tmp,A_type(i),beta,A_type,node1,node2,node3,node4,node5,node6);
    
  }
  
  
  arma::ucolvec order_index = sort_index(Qtmp,"descend");
  
  
  arma::vec Q_sort = sort(Qtmp,"descend");
  
  
  int ii = 0;
  while (lambda + (ii+1)*Q_sort(ii) >= sum(Q_sort.head(ii+1))){
    ii = ii + 1;
    if (ii == A_size){
      break;
    }
  }
  ii = ii - 1;
  
  
  arma::vec target_A_set = A_type(order_index.head(ii+1));
  
  //Rcout << "Initial value,   F = " <<target_A_set << std::endl;
  
  arma::rowvec Pi_gradient_value(len_basis*A_size);
  
  //  Rcout << "Initial value,   F = " << (std::find(target_A_set.begin(), target_A_set.end(), A_tmp) != target_A_set.end()) << std::endl;
  // Rcout << "Initial value,   F = " <<  target_A_set.end() << std::endl;
  
  if (std::find(target_A_set.begin(), target_A_set.end(), A_tmp) != target_A_set.end()){
    
    arma::rowvec basis_sum(len_basis*A_size);
    basis_sum.fill(0);
    
    for (int j = 0; j < ii+1; j++) {
      basis_sum += Basis_func(X_tmp,target_A_set(j),A_type,node1,node2,node3,node4,node5,node6)/lambda;
    }
    
    Pi_gradient_value = Basis_func(X_tmp,A_tmp,A_type,node1,node2,node3,node4,node5,node6)/lambda- basis_sum/(ii+1); 
    
  }else{
    
    Pi_gradient_value.fill(0);
    
  }
  
  
  return(Pi_gradient_value);
}



//    ----------------------------------------------------------------
//    Lagrange Conjugate Policy Function
//    ----------------------------------------------------------------

// [[Rcpp::export]]


double mu(const int A_tmp,
          const arma::rowvec& X_tmp,
          const arma::colvec& A_type,
          const arma::mat& beta,
          double lambda,
          double node1,
          double node2,double node3,double node4,double node5,double node6)
{
  
  //int X_size = X_tmp.n_elem;
  //int len_basis = X_size*3+1;
  int A_size = A_type.n_elem;
  
  arma::vec Qtmp(A_size); 
  
  for (int i = 0; i < A_size; i++) {
    
    Qtmp(i) = Qfunc(X_tmp,A_type(i),beta,A_type,node1,node2,node3,node4,node5,node6);
    
  }
  
  arma::ucolvec order_index = sort_index(Qtmp,"descend");
  arma::vec Q_sort = sort(Qtmp,"descend");
  
  
  double ii = 0;
  while (lambda + (ii+1)*Q_sort(ii) >= sum(Q_sort.head(ii+1))){
    ii = ii + 1;
    if (ii ==  A_size){
      break;
    }
  }
  ii = ii - 1;
  
  arma::vec target_A_set = A_type(order_index.head(ii+1));
  
  double G_value = (sum(Qtmp(order_index.head(ii+1)))/lambda-1)/(ii+1); 
  double mu_value = std::max(0.0,G_value-Qfunc(X_tmp,A_tmp,beta,A_type,node1,node2,node3,node4,node5,node6)/lambda);
  
  return(mu_value);
  
}  


//    ----------------------------------------------------------------
//    Gradient Lagrange Conjugate Policy Function
//    ----------------------------------------------------------------

// [[Rcpp::export]]


arma::rowvec mu_grad_beta(const int A_tmp,
                          const arma::rowvec& X_tmp,
                          const arma::colvec& A_type,
                          const arma::mat& beta,
                          double lambda,
                          double node1,
                          double node2,double node3,double node4,double node5,double node6,
                          int X_size)
{
  
  //int X_size = X_tmp.n_elem;
  int len_basis = X_size*9+1;
  int A_size = A_type.n_elem;
  
  
  
  arma::vec Qtmp(A_size); 
  
  for (int i = 0; i < A_size; i++) {
    
    Qtmp(i) = Qfunc(X_tmp,A_type(i),beta,A_type,node1,node2,node3,node4,node5,node6);
    
  }
  
  
  arma::ucolvec order_index = sort_index(Qtmp,"descend");
  
  
  arma::vec Q_sort = sort(Qtmp,"descend");
  
  
  int ii = 0;
  while (lambda + (ii+1)*Q_sort(ii) >= sum(Q_sort.head(ii+1))){
    ii = ii + 1;
    if (ii == A_size){
      break;
    }
  }
  ii = ii - 1;
  
  
  arma::vec target_A_set = A_type(order_index.head(ii+1));
  
  //Rcout << "Initial value,   F = " <<  target_A_set << std::endl;
  
  //Rcout << "Initial value,   F = " <<target_A_set << std::endl;
  
  arma::rowvec mu_gradient_value(len_basis*A_size);
  
  //  Rcout << "Initial value,   F = " << (std::find(target_A_set.begin(), target_A_set.end(), A_tmp) != target_A_set.end()) << std::endl;
  // Rcout << "Initial value,   F = " <<  target_A_set.end() << std::endl;
  
  if (std::find(target_A_set.begin(), target_A_set.end(), A_tmp) == target_A_set.end()){
    
    arma::rowvec basis_sum(len_basis*A_size);
    basis_sum.fill(0);
    
    for (int j = 0; j < ii+1; j++) {
      basis_sum += Basis_func(X_tmp,target_A_set(j),A_type,node1,node2,node3,node4,node5,node6)/lambda;
    }
    
    //Rcout << "Initial value,   F = " <<  basis_sum/(ii+1) << std::endl;
    
    
    mu_gradient_value =  basis_sum/(ii+1)-Basis_func(X_tmp,A_tmp,A_type,node1,node2,node3,node4,node5,node6)/lambda; 
    
  }else{
    
    mu_gradient_value.fill(0);
    
  }
  
  return(mu_gradient_value);
}



//    ----------------------------------------------------------------
//    Lagrange Lambda Function (-1/2<Lambda<0)
//    ----------------------------------------------------------------

// [[Rcpp::export]]


double lagr_lambda(const arma::rowvec& X_tmp,
                   double lambda,
                   const arma::mat& omega)
{
  
  int X_size = X_tmp.n_elem;
  //int len_basis = 5+1;
  
  arma::rowvec lar_Phi(6);
  lar_Phi.fill(1);
  
  arma::rowvec lar_center(X_size);
  
  lar_Phi(0,1) = exp(-8.0*sum(pow(X_tmp - lar_center.fill(0),2)));
  lar_Phi(0,2) = exp(-8.0*sum(pow(X_tmp - lar_center.fill(0.25),2)));
  lar_Phi(0,3) = exp(-8.0*sum(pow(X_tmp - lar_center.fill(0.5),2)));
  lar_Phi(0,4) = exp(-8.0*sum(pow(X_tmp - lar_center.fill(0.75),2)));
  lar_Phi(0,5) = exp(-8.0*sum(pow(X_tmp - lar_center.fill(1),2)));
  
  
  // for (int j = 0; j < X_size; j++) {
  //   
  //   lar_Phi(j*5+1) = exp(-8.0*pow(X_tmp(j)-0,2));
  //   lar_Phi(j*5+2) = exp(-8.0*pow(X_tmp(j)-0.25,2));
  //   lar_Phi(j*5+3) = exp(-8.0*pow(X_tmp(j)-0.5,2));
  //   lar_Phi(j*5+4) = exp(-8.0*pow(X_tmp(j)-0.75,2));
  //   lar_Phi(j*5+5) = exp(-8.0*pow(X_tmp(j)-1.0,2));
  //   
  // }
  
  //Rcout << "Initial value,   F = " <<  (2*(1+exp(-( arma::as_scalar(lar_Phi * omega.t()))  )))<< std::endl;
  
  double  lagr_value = -lambda/(2.0*(1+exp(-( arma::as_scalar(lar_Phi * omega.t()))  )));
  
  return(lagr_value);
  
}

//    ----------------------------------------------------------------
//    Gradient Lagrange Lambda Function (-1/2<Lambda<0)
//    ----------------------------------------------------------------

// [[Rcpp::export]]


arma::rowvec lagr_lambda_grad(const arma::rowvec& X_tmp,
                              double lambda,
                              const arma::mat& omega)
{
  
  int X_size = X_tmp.n_elem;
  //int len_basis = 5+1;
  
  arma::rowvec lar_Phi(6);
  lar_Phi.fill(1);
  
  arma::rowvec lar_center(X_size);
  
  lar_Phi(0,1) = exp(-8.0*sum(pow(X_tmp - lar_center.fill(0),2)));
  lar_Phi(0,2) = exp(-8.0*sum(pow(X_tmp - lar_center.fill(0.25),2)));
  lar_Phi(0,3) = exp(-8.0*sum(pow(X_tmp - lar_center.fill(0.5),2)));
  lar_Phi(0,4) = exp(-8.0*sum(pow(X_tmp - lar_center.fill(0.75),2)));
  lar_Phi(0,5) = exp(-8.0*sum(pow(X_tmp - lar_center.fill(1),2)));
  
  // 
  // for (int j = 0; j < X_size; j++) {
  //   
  //   lar_Phi(j*5+1) = exp(-8.0*pow(X_tmp(j)-0,2));
  //   lar_Phi(j*5+2) = exp(-8.0*pow(X_tmp(j)-0.25,2));
  //   lar_Phi(j*5+3) = exp(-8.0*pow(X_tmp(j)-0.5,2));
  //   lar_Phi(j*5+4) = exp(-8.0*pow(X_tmp(j)-0.75,2));
  //   lar_Phi(j*5+5) = exp(-8.0*pow(X_tmp(j)-1,2));
  //   
  // }
  // 
  arma::rowvec lagr_grad_value = (-lambda*lar_Phi/2.0)*(1.0-1.0/(1+exp(-arma::as_scalar(lar_Phi*omega.t()))))*(1.0/(1.0+exp(-arma::as_scalar(lar_Phi*omega.t()))));
  
  return(lagr_grad_value);
  
}


//    ----------------------------------------------------------------
//    Gradient Kernel Loss 
//    ----------------------------------------------------------------

// [[Rcpp::export]]


List loss_grad(const arma::mat& train,
               const arma::mat& beta,
               const arma::mat& omega,
               const arma::colvec& A_type,
               double eta_b,
               double eta_omg,
               double gamma,
               double lambda,
               const arma::mat& kernel_mat,
               double node1,
               double node2,double node3,double node4,double node5,double node6,
               int X_size)
{
  
  int N = train.n_rows;
  int beta_size = beta.n_elem;
  int omega_size = omega.n_elem;
  //int subset_size = subset.n_elem;
  // int X_size = (omega_size - 1)/5;
  
  arma::rowvec error_vec(N); 
  arma::mat  grad_beta_vec(N,beta_size);
  arma::mat  grad_omega_vec(N,omega_size); 
  
  for (int i = 0; i < N; i++) {
    
    
    //Rcout << "count " <<  k << std::endl;
    
    //Rcout << "V = " <<  train(i,X_size+1) << std::endl;
    //Rcout << "V = " <<  V_func(train.row(i).tail(X_size),A_type,beta,lambda,center,bandw,Num_center) << std::endl;
    //Rcout << "pi = " <<  Pi_policy(train(i,X_size),train.row(i).head(X_size),A_type,beta,lambda,center,bandw,Num_center) << std::endl;
    //Rcout << "mu = " <<  mu(train(i,X_size),train.row(i).head(X_size),A_type,beta,lambda,center,bandw,Num_center) << std::endl;
    //Rcout << "larg = " <<  lagr_lambda(train.row(i).head(X_size),lambda,omega) << std::endl;
    //Rcout << "V2 = " <<      V_func(train.row(i).head(X_size),A_type,beta,lambda,center,bandw,Num_center) << std::endl;
    
    error_vec(i) = train(i,X_size+1) + 
      gamma*V_func(train.row(i).tail(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6) -
      lambda*Pi_policy(train(i,X_size),train.row(i).head(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6) + 
      lambda/2.0 + 
      mu(train(i,X_size),train.row(i).head(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6) - 
      lagr_lambda(train.row(i).head(X_size),lambda,omega) - 
      V_func(train.row(i).head(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6);
    
    //Rcout << "err = " <<  error_vec(k) << std::endl;
    
    grad_beta_vec.row(i) = gamma*V_grad(train.row(i).tail(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6,X_size)-
      V_grad(train.row(i).head(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6,X_size) -
      (lambda*Pi_grad(train(i,X_size),train.row(i).head(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6,X_size) -
      mu_grad_beta(train(i,X_size),train.row(i).head(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6,X_size));
    
    
    grad_omega_vec.row(i) =  -lagr_lambda_grad(train.row(i).head(X_size),lambda,omega);
    
    
    //Rcout << "omega= " <<  grad_omega_vec.row(k) << std::endl;
    
  }
  
  // Rcout << "V = " <<  gamma*V_func(train.row(subset(0)).tail(X_size),A_type,beta,lambda) << std::endl;
  
  arma::rowvec total_grad_beta = (1.0/((N-1)*N))*((kernel_mat * error_vec.t()).t() * grad_beta_vec);
  arma::rowvec total_grad_omega = (1.0/((N-1)*N))*((kernel_mat * error_vec.t()).t() * grad_omega_vec);
  
  
  List ret;
  ret["G.beta"] = total_grad_beta ;
  ret["G.omega"] = total_grad_omega;
  //ret["G.beta"] =  grad_beta_vec;
  return(ret);
  
}

//    ----------------------------------------------------------------
//    Objective Function
//    ----------------------------------------------------------------

// [[Rcpp::export]]


double kernel_obj(const arma::mat& train,
                  const arma::mat& beta,
                  const arma::mat& omega,
                  const arma::colvec& A_type,
                  double gamma,
                  double lambda,
                  const arma::mat& kernel_mat,
                  double node1,
                  double node2,double node3,double node4,double node5,double node6,
                  int X_size)
{
  
  int N = train.n_rows;
  arma::rowvec error_vec(N); 
  
  for (int k = 0; k < N; k++) {
    
 
    error_vec(k) = train(k,X_size+1) + 
      gamma*V_func(train.row(k).tail(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6) -
      lambda*Pi_policy(train(k,X_size),train.row(k).head(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6) + 
      lambda/2.0 + 
      mu(train(k,X_size),train.row(k).head(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6) - 
      lagr_lambda(train.row(k).head(X_size),lambda,omega) - 
      V_func(train.row(k).head(X_size),A_type,beta,lambda,node1,node2,node3,node4,node5,node6);
    
  }
  
  //Rcout << "V = " <<  1 << std::endl;
  
  //
  double obj_value =arma::as_scalar((1.0/((N-1)*N))*((kernel_mat * error_vec.t()).t() * error_vec.t()));
  
  return(obj_value);
  
}
