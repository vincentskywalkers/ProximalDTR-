###########################################################
#            A Example for proximal DTR                   #
#Thi example contains total 12 different treatment choices#
###########################################################

#Stage Num
stage = 25
#Subject Num 
N = 75
#Burn in stages Num (set zeor not brun in stage)
ini.M = 0

#Genearate Dataset
set.seed(123) 
train = data_generate_nopair(stage,N,ini.M)

#X is the longitudinal covariates matrix:
#eg. row1 : rowN is subjects covariates matrix at stage 1;
#eg. row(N+1) : row(2*N) is subjects covariates matrix at stage 2
#etc...
X = train$X

#A is the longitudinal assigned treatment 
A = train$A
#R is the longitudinal reward 
R = train$R

#Note that, X has 1 more stage than A and R, since at the end of stage, we only observe X but not A or R. 

#Fit.model
#n_ID is the number of subjects attending experiments 

fit = proximalDTR(X = X, A=A, R=R, n_ID = N, stage=stage, gamma=0.9,
                  lambda.set = c(0.1,0.5,1,2.5,5),step.beta = 0.01, step.omega = 0.01,
                  desc.rate = 0.005, max.iter = 3000, max.iter.cv = 1000, bw = 0.5,
                  cv=F, trace =TRUE)


#predict prob and recommend trt

#generate one sample point 
X_01 <- matrix(rnorm(1,100,5),1,1)
X_02 <- matrix(rnorm(1,12,2),1,1)
X_03 <- matrix(rnorm(1,72,2),1,1)
X_0 = c(X_01,X_02,X_03)

#start prediction
pred.fit <- predict_rl(fit, X_0)

#obatin prediction and recommendation
recommend.trt = pred.fit$recommend.trt
recommend.prob = pred.fit$prob





