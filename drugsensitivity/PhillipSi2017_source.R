# By Phillip Si, Summer 2017, HKUST
# 

# Download the dataset to a current directory
# set working directory to the current directory
# For example:
# > setwd("~/Dropbox/2019_summer/drugsensitivity-3/")
# > dir()
# [1] "cell_line_features.xlsx"     "submission.6M_drug_pair.csv" "train3M_drug_pair.csv"      

library(Matrix)
library(foreach)
library(parallel)
library(doParallel)
#library(Rfast)
library(glmnet)

##### Construct the dataset which includes 
##  X.train 2,400,000 * 339,730
##  Y.train 2,400,000 * 1
##  X.test  6,000,000 * 339,730
##  The feature contains 339,730 = 265(baseline) + 1250*265(gene mutation) + 31*265(cancer type)

## Open the file "cell_line_features.xlsx" and save it to "cell_line_features.csv"
#
cellfea.DF = read.csv('cell_line_features.csv',header=T,as.is=T) #Data frame and as.is keeps the characters  instead of converting into factor
gene.names = cellfea.DF[-1,1]
cell.gene = t(apply(as.matrix(cellfea.DF[-1,-1]),2,as.numeric))
colnames(cell.gene) = gene.names
cell.cancer = as.character(cellfea.DF[1,-1])

#dgComp = read.csv('train3M_drug_pairMAIN.csv',header=T)
dgComp = read.csv('train3M_drug_pair.csv',header=T)

#relabel names of drug by 1-265 
drug.names = unique(c(as.vector(dgComp$Drug1),as.vector(dgComp$Drug2)))
ll = factor(c(as.vector(dgComp$Drug1),as.vector(dgComp$Drug2)),levels=drug.names,labels=1:265)
dgComp$Drug1 = ll[1:3000000]
dgComp$Drug2 = ll[3000001:6000000]

#relabel names of cell lines by 1-990 corresponding to order in origin data
cell.names = rownames(cell.gene)
sample.names = factor(as.character(dgComp[,1]),levels=cell.names,labels=1:990)
dgComp[,1] = as.numeric(sample.names)

#relabel names of gene by 1-1250 corresponding to order in origin data
colnames(cell.gene) = 1:1250

#relabel names of cancer silghtly adjusted from origin data
cancer.names = unique(cell.cancer)
cancer.names = c(cancer.names[-3],'NOT CLASSIFIED') 
cell.cancer.l = factor(cell.cancer,levels=cancer.names,labels=1:31)

## Sparse design matrix
# 1st (drug baseline)
ii = as.numeric(dgComp[,2])
jj = as.numeric(dgComp[,3])
X1 = sparseMatrix(rep(1:3000000,2),c(ii,jj),x=rep(c(1,-1),each=3000000),dims=c(3000000,265))
# 2nd (gene mutation 3,000,000 * (1250 * 265))
cell.gene.s = Matrix(cell.gene,sparse=T)

index.cell.dg1 = (as.numeric(dgComp[,2]) - 1) * 990 + dgComp[,1]
index.cell.dg2 = (as.numeric(dgComp[,3]) - 1) * 990 + dgComp[,1]
index.cell.1 = sparseMatrix(i=1:3000000,j=index.cell.dg1, x=rep(1,3000000),dims=c(3000000,990*265))
index.cell.2 = sparseMatrix(i=1:3000000,j=index.cell.dg2, x=rep(-1,3000000),dims=c(3000000,990*265))
index.cell = index.cell.1 + index.cell.2 #3,000,000 * (990*265)

cell.gene.block = bdiag(rep(list(A=cell.gene.s),265)) #(990*265) * (1250*265)
X2 = index.cell %*% cell.gene.block


# 3rd (cancer type 3,000,000 * (265*31))
index.cancer = rep(1,3000000)
for(ii in 1:3000000){
  index.cancer[ii] = cell.cancer.l[dgComp[ii,1]]
}
index.cancer.dg1 = (as.numeric(dgComp[,2]) - 1)*31 + index.cancer
index.cancer.dg2 = (as.numeric(dgComp[,3]) - 1)*31 + index.cancer
X3 = sparseMatrix(1:3000000,index.cancer.dg1,x=rep(1,3000000),dims=c(3000000,31*265)) - sparseMatrix(1:3000000, index.cancer.dg2,x=rep(1,3000000),dims=c(3000000,31*265))
# Train sample index and Test sample index
Y = dgComp[,4]
index.train = rep(append(rep(T,1),rep(F,99)),30000)#rep(T,3000000), change to this if running full data
#index.train2 = rep(append(rep(F,99),rep(T,1)),30000) for testing data
index.0 = which(dgComp[,4] == 0)
index.train[index.0] = F
index.notna = (is.na(dgComp[,4])==F) #not NA index
index.train = index.train & index.notna
index.test= (is.na(dgComp[,4]) ==T)

X.train = cbind(X1[index.train,],X2[index.train,],X3[index.train,])
X.test = cbind(X1[index.test,],X2[index.test,],X3[index.test,])
Y.train = Y[index.train]
#Y.test = Y[index.test] for testing on partitioned

######  Model fitting
submission = read.csv('submission.6M_drug_pair.csv',head=T)

## Coefficient estimation
p.fac = rep(0.6,(265+1250*265+31*265))
p.fac[1:265] = 0

# Be careful! The following run may take a long time if you have bad parameters...
ptm <- proc.time()
fit = glmnet(X.train,Y.train,family='binomial',penalty.factor=p.fac,intercept=T)
proc.time()-ptm
#> proc.time()-ptm
#   user  system elapsed 
#  4.028   0.369   4.421 

Y.test.logistic = predict(fit,newx=X.test,type = "response", s=c(fit$lambda[90]))
Y.test.class = matrix(NA,length(Y.test.logistic[,1]),1)
Y.test.class[Y.test.logistic[,1]<0.5,1] = -1
Y.test.class[Y.test.logistic[,1]>0.5,1] = 1
submission[,2] = Y.test.class[,1]
write.csv(submission,file='submission_phillipsi_ElasticLog.csv',row.names=F)

## for feature rendering in Excel file
coefficientMatrix = coef(fit)
featureMatrix = matrix(NA,265,1282)
foreach(i=1:1282)%do%
{
  featureMatrix[,i]=coefficientMatrix[(265*i-264):(265*i),90]
  print(i)
}
rownames(featureMatrix)=drug.names
cancer.names = unique(cell.cancer)
colnames(featureMatrix)=c("Common Ranking", gene.names, cancer.names)
featureVectors = data.frame(featureMatrix)
write.csv(featureVectors, file='featureCoefficients.csv', row.names=T)
