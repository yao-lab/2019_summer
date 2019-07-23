load('dgData_lm.RData')
#### Read the workspace, which includes 
##  X.train 2,400,000 * 339,730
##  Y.train 2,400,000 * 1
##  X.test  6,000,000 * 339,730
##  The feature contains 339,730 = 265(baseline) + 1250*265(gene mutation) + 31*265(cancer type)
library(Matrix)
library(foreach)
library(glmnet)

#########################LinearModel########################
#### penalty factor for elastic-net penalty
p.fac = rep(1,(265+1250*265+31*265))
p.fac[1:265] = 0
alpha = c(0.5,0.75,1)  ## alpha is weight parameter between L_1 penalty and L_2 penalty

##### 2-dim cross-validation 
fit.lm.cv = list()
mse = numeric()
foldid = sample(1:10,size=length(Y.train),replace=T) ## fix the fold before cross validation

for(ii in 1:length(alpha)){
## for fix alpha , cross validation on lambda
fit.lm.cv[[ii]] = cv.glmnet(X.train,Y.train,penalty.factor=p.fac,intercept=T,nfold=10,foldid=foldid,alpha=alpha[ii],lambda.min.ratio = 1e-6,nlambda=100)
# this takes a while...

## record the minimum mean square error 
mse[ii] = min(fit.lm.cv[[ii]]$cvm) 
}

## choose optimal alpha
opt.alp.index = which(mse == min(mse))  

#### prediction
fit = fit.lm.cv[[opt.alp.index]] 
opt.lam.index = which(fit$cvm == min(fit$cvm))
Y.test = predict(fit, newx=X.test,s=fit$lambda[opt.lam.index] )
