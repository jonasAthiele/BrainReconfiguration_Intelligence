# Generate a g-factor from the behavioral data
# Code adopted from Dubois et al. 2018: http://dx.doi.org/10.1098/rstb.2017.0284
# Original code: https://github.com/adolphslab/HCP_MRI-behavior

library(ggplot2)
library(psych)
library(lavaan)
library(Hmisc)
library(corrplot)
library(semPlot)
library(colorRamps)
# Helpers functions
# compute Comparative Fit Index for a factor analysis 
CFI <-function(x){
  return((1-((x$STATISTIC-x$dof))/(x$null.chisq-x$null.dof)))
}
# compute Comparative Fit Index for a bifactor analysis 
CFI_biv <-function(x){
  return((1-((x$stats$STATISTIC-x$stats$dof))/(x$stats$null.chisq-x$stats$null.dof)))
}
# compute implied matrix for a factor analysis
impliedMatrix<-function(x){
  if (dim(x$loadings)[2]==1) {
    imp <- x$loadings %*% t(x$loadings)
    } 
  else {imp <- x$loadings %*% x$Phi %*% t(x$loadings)}
  diag(imp)<- diag(imp) + x$uniquenesses
  return(imp)
  }
# compute implied matrix for a bifactor analysis
impliedMatrix_biv<-function(x){
  Gloadings     <- x$schmid$sl[,1]
  Floadings     <- x$schmid$sl[,2:(ncol(x$schmid$sl)-3)]
  uniquenesses  <- x$schmid$sl[,ncol(x$schmid$sl)-1]
  imp           <- Gloadings %*% t(Gloadings) + Floadings %*% t(Floadings)
  diag(imp)     <- diag(imp) + uniquenesses
  return(imp)
}

cogdf <- read.csv(file = 'cogScores_1186Subjetcs.csv') # read cognitive scores

out = fa.parallel(cogdf,plot=F)#error.bars=T,se.bars=F,
faValues = out$fa.values
faSim    = out$fa.sim
faSimR   = out$fa.simr

fm     <- "minres"       # use maximum likelihood estimator
rotate <- "oblimin"   # use oblimin factor rotation

fitInds <- matrix(, nrow = 2, ncol = 9)
rownames(fitInds) <- c('s1','b4')
colnames(fitInds) <- c('CFI','RMSEA','SRMR','BIC','om_h','om_s1','om_s2','om_s3','om_s4')
# observed covariance matrices
obs       <-  cov(cogdf)
lobs      <-  obs[!lower.tri(obs)]

#SINGLE FACTOR
model = 1
f1     <- fa(cogdf,nfactors=1)
imp    <-  impliedMatrix(f1)
limp   <-  imp[!lower.tri(imp)]
fitInds[model,1] <-  CFI(f1)
fitInds[model,2] <-  f1$RMSEA[1]
fitInds[model,3] <-  sqrt(mean((limp - lobs)^2))
fitInds[model,4] <-  f1$BIC

# Parallel analysis to compute numbers of factors for bi-factor model
fa.parallel(cogdf,n.obs=NULL,fm="minres",fa="both",nfactors=1, 
            main="Parallel Analysis Scree Plots",
            n.iter=20,error.bars=FALSE,se.bars=FALSE,SMC=FALSE,ylabel=NULL,show.legend=TRUE,
            sim=TRUE,quant=.95,cor="cor",use="pairwise",plot=TRUE,correct=.5)




# BI-FACTOR MODEL
model = 2
b4      <- omega(cogdf,nfactors=4,fm=fm,key=NULL,flip=FALSE,
                 digits=3,title="Omega",sl=TRUE,labels=NULL, plot=FALSE,
                 n.obs=NA,rotate=rotate,Phi = NULL,option="equal",covar=FALSE)
imp     <-  impliedMatrix_biv(b4)
limp    <-  imp[!lower.tri(imp)]
fitInds[model,1] <-  CFI_biv(b4)
fitInds[model,2] <-  b4$schmid$RMSEA[1]
fitInds[model,3] <-  sqrt(mean((limp - lobs)^2))
fitInds[model,4] <-  b4$stats$BIC
fitInds[model,5] <-  b4$omega_h
fitInds[model,6:9] <-  b4$omega.group[-1,3]

print(fitInds,digits=3)

print(b4)

# %%R -w 500 -h 500 -o b4Scores

diagram(b4,digits=3,cut=0.2)

# export scores
b4Scores    <- factor.scores(cogdf,b4$schmid$sl[,1:5])$scores
write.csv(b4Scores[,1], 'gfac_bi_explore_1186Subjects.csv')





