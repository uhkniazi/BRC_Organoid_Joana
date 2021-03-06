# File: 11_de_analysis_MCMC.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 27/01/2017
# Desc: DE analysis for the count data with Stan

## set variables and source libraries
source('header.R')
## libraries to load
library('RMySQL')

### space for internal functions
## gamma shape function
## this function is from the book: Bayesian data analysis
gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}
### end of functions


##### connect to mysql database to get count matrix
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 11) AND (MetaFile.comment like "%Expression set object%")')
dfCountMatrix = dbGetQuery(db, q)
dfCountMatrix
# load the expression set object
dbDisconnect(db)

## load the count matrix
n = paste0(dfCountMatrix$location, dfCountMatrix$name)
load(n)

dim(oExp.norm)
temp = exprs(oExp.norm)
# perform independent filtering 
iCutoff = 4
fKeep = rowMeans(temp) > iCutoff
table(fKeep)
oExp.norm = oExp.norm[fKeep,]
dim(oExp.norm)

## prepare data to fit model
# fCondition = as.character(oExp.norm$fCondition)
# table(fCondition)
# # drop the si
# f = which(fCondition == 'SI')
# f
# oExp.norm = oExp.norm[,-f]
# dim(oExp.norm)
# # keep the others
# fCondition = as.character(oExp.norm$fCondition)
# table(fCondition)
# fCondition[fCondition != 'CLP_SI'] = 'Media'
# fCondition = factor(fCondition, levels = c('Media', 'CLP_SI'))
mDat = exprs(oExp.norm)
fCondition = oExp.norm$fCondition
levels(fCondition)

#### define function to optimize
library(LearnBayes)
library(MASS)
library(numDeriv)
## define log posterior
mylogpost = function(theta, data){
  size = exp(theta['size'])
  beta0 = theta['beta0']
  beta1 = theta['beta1']
  beta2 = theta['beta2']
  # hyperparameters for the hierarchichal standard deviation parameter
  betaSigma = exp(theta['betaSigma'])# + 0.1 # maybe add a one to prevent zero variance
  ##if (betaSigma <= 0) return(-Inf)
  # hyper-hyperparameters for the hierarchichal standard deviation parameter
  ivBetaParam = data$betaParam
  #ivSizeParam = data$sizeParam
  x = data$pred
  y = data$resp
  ## create model matrix to get fitted value
  mModMatrix = model.matrix(y ~ x)
  ## likelihood function
  lf = function(dat, pred){
    return(log(dnbinom(dat, size = size, mu = pred)))
  }
  mCoef = matrix(c(beta0, beta1, beta2), nrow = 3)
  mu = exp(mModMatrix %*% mCoef)
  val = sum(lf(y, mu))
  ## notice that we share betaSigma hyperparameter with all 3 coefficients, if we desire even more shrinkage
  ## then set beta0 the baseline to a constant e.g. 10^2, this will free up this parameter to shrink even further
  ## and hence shrink the beta1 and beta2. Furthermore, the size parameter is not restricted and has a very large range
  ## this allows the model to assign extra unexplained variance to size while shrinking the coefficients. If size is fixed or too 
  ## restrictive, then you may get the correct betas but their standard deviations will be high i.e. the betaSigma won't shrink enough.
  ret = val + dnorm(beta0, 0, betaSigma, log=T) + dnorm(beta1, 0, betaSigma, log=T) + dnorm(beta2, 0, betaSigma, log=T) + 
    dunif(size, 1, 1e10, log=T) +
    dgamma(betaSigma, shape = ivBetaParam['shape'], rate = ivBetaParam['rate'], log=T)
  return(ret)
}

getOptimizedPValue = function(obj){
  se = sqrt(abs(diag(obj$var)))[c('beta0', 'beta1', 'beta2')]
  m = obj$mode[c('beta0', 'beta1', 'beta2')]
  pnorm(-abs(m/se))*2
}

getOptimizedSummary = function(obj){
  se = sqrt(abs(diag(obj$var)))['beta2']
  m = obj$mode['beta2']
  p = pnorm(-abs(m/se))*2
  ret  = c('Coef' = round(m, 3), 'SE'=round(se, 3), 'P-Value'=signif(p, 3))
  names(ret) = c('Coef', 'SE', 'P-Value')
  return(ret)
}

############################### fit the model on multiple genes
## calculate the distribution for sizes
# ivSizes = apply(mDat, 1, function(x){
#   getsize = function(X) fitdistr(as.integer(X), 'negative binomial')$estimate['size'];
#   tryCatch(expr = getsize(x), error=function(e) NULL)
# })
# ivSizes = unlist(ivSizes)
# 
# iSizeParam = unlist(gammaShRaFromModeSD(sd(ivSizes)/2, 2*sd(ivSizes)))
mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1), data=data)
  # calculate hessian
  fit$hessian = (hessian(logpost, fit$par, data=data))
  colnames(fit$hessian) = names(mode)
  rownames(fit$hessian) = names(mode)
  options(warn = 0)
  mode = fit$par
  h = -solve(fit$hessian)
  stuff = list(mode = mode, var = h, converge = fit$convergence == 
                 0)
  return(stuff)
}


modelFunction = function(dat){
  ## prepare data before fitting
  dfData = data.frame(resp=as.integer(mDat[dat,]), pred=fCondition)
  temp = fitdistr(dfData$resp, 'negative binomial')$estimate['size']
  names(temp) = NULL
  # set starting values for optimizer and 
  # set parameters for optimizer
  lData = list(resp=dfData$resp, pred=dfData$pred)
  # hyper-hyperparameters for deflection variance, calculated using the data
  #lData$betaParam = unlist(gammaShRaFromModeSD(mode = log(sd(dfData$resp+0.5)/2), sd = log(2*sd(dfData$resp+0.5))))
  # use a jeffery's prior 
  lData$betaParam = c('shape'=0.5, 'rate'=0.0001)
  # # gamma hyperparameter for the size 
  # lData$sizeParam = iSizeParam
  start = c('size'=log(max(temp, 1)), 'beta0'=log(mean(dfData$resp)), 'beta1'=0, 'beta2' = 0, 
            'betaSigma'=log(log(sd(dfData$resp))))
  #op = optim(start, mylogpost, control = list(fnscale = -1, maxit=1000), data=lData)
  # see results of optimiser
  #if (op$convergence) cat('Optimizer convergence failed')
  ## you can see the starting values
  #start2 = op$par
  #names(start2) = names(start)
  mylaplace(mylogpost, start, lData)
}

mDat = exprs(oExp.norm)
iIndex = 1:nrow(mDat)

lGlm = lapply(iIndex, function(iIndexSub) {
  tryCatch(modelFunction(iIndexSub), error=function(e) NULL)
})

names(lGlm) = rownames(mDat)

table(sapply(lGlm, is.null))
# remove the elements of list that are empty
lGlm.sub = lGlm[!sapply(lGlm, is.null)]
length(lGlm.sub)

table(sapply(lGlm.sub, function(x) x$converge))

mSummary = t(sapply(lGlm.sub, getOptimizedSummary))
table(fCondition)

dfResults = data.frame('ENTREZID' = rownames(mSummary), logFC=mSummary[,'Coef'], SE=mSummary[,'SE'], P.Value=mSummary[,'P-Value'])
rownames(dfResults) = dfResults$ENTREZID

# ## the genes on which the model did not fit, use MCMC to fit the model
# ## prepare stan dso object
# library(rstan)
# stanDso = rstan::stan_model(file='nb_glm.stan')
# 
# 
# ## functions to extract information out of stan object
# getStanSummary = function(obj){
#   return(apply(extract(obj)$beta, 2, sd))
# }
# getStanSD = function(obj){
#   return(apply(extract(obj)$beta, 2, sd))
# }
# getStanMean = function(obj){
#   return(apply(extract(obj)$beta, 2, mean))
# }
# getStanPValue = function(obj){
#   pnorm(-abs(getStanMean(obj)/getStanSD(obj)))*2
# }
# 
# ## get the names of genes on which optimiser did not converge
# table(sapply(lGlm, is.null))
# f = sapply(lGlm, is.null)
# cvFailed = names(lGlm[f])
# 
# ## write a model function
# modelFunction2 = function(dat){
#   ## prepare data before fitting
#   dfData = data.frame(resp=as.integer(mDat[dat,]+1), pred=fCondition)
#   lData = list(resp=dfData$resp, pred=dfData$pred)
#   lData$betaParam = unlist(gammaShRaFromModeSD(mode = log(sd(dfData$resp)/2), sd = log(2*sd(dfData$resp))))
#   lStanData = list(Ntotal=nrow(dfData), Ncols=2, modmatrix= model.matrix(resp ~ pred, data=dfData), y=dfData$resp,
#                    betaShape=lData$betaParam['shape'], betaRate=lData$betaParam['rate'])
#   fit.stan = sampling(stanDso, data = lStanData, pars = c('beta'),
#                       iter = 5000, chains = 4)
#   return(c('Coef'=getStanMean(fit.stan)[2], 'SE'=getStanSD(fit.stan)[2], 'P-Value'=getStanPValue(fit.stan)[2]))
# }
# 
# lGlm.stan = lapply(cvFailed, function(cvIndexSub) {
#   tryCatch(modelFunction2(cvIndexSub), error=function(e) NULL)
# })

dfResults$adj.P.Val = p.adjust(dfResults$P.Value, 'BH')

### create some volcano plots with annotations
library(org.Mm.eg.db)
# add annotation to the data set after selecting comparison
rn = rownames(dfResults)
df = select(org.Mm.eg.db, as.character(rn), c('SYMBOL'), 'ENTREZID')
head(df); head(dfResults)
df = df[!duplicated(df$ENTREZID),]
dim(df); dim(dfResults)
rownames(df) = df$ENTREZID
rn = rownames(dfResults)
dfResults = cbind(dfResults[rn,], df[rn,'SYMBOL'])
colnames(dfResults)[6] = 'SYMBOL'
dfPlot = na.omit(dfResults)
dim(dfPlot); dim(dfResults)

## write csv file

write.csv(dfPlot, file='Results/DEAnalysis_CLP_SI.vs.Media_sizePrior_3coef_mar01.xls')

dfGenes = data.frame(P.Value=dfPlot$P.Value, logFC=dfPlot$logFC, adj.P.Val = dfPlot$adj.P.Val, SYMBOL=dfPlot$SYMBOL)

f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim)
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values, ignore this set probs=0
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}

f_plotVolcano(dfGenes, 'CLP_SI vs Media')

dfResults = na.omit(dfResults)

par(mfrow=c(2,2))
hist(dfResults$logFC)
hist(dfResults$P.Value)
hist(dfResults$adj.P.Val)

table(dfResults$adj.P.Val < 0.01)

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$p.adj < p.cut] = 'red'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$logfc, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC', ylim=c(-2, 2), cex=0.6)
}

dfGenes = data.frame(p.adj=dfResults$adj.P.Val, logfc=dfResults$logFC)
iMean = rowMeans(mDat[as.character(dfResults$ENTREZID),])

plotMeanFC(log(iMean), dfGenes, 0.05, 'MA Plot - Bayes')


### make some heatmaps
fSamples = fCondition
rn = rownames(dfResults[dfResults$adj.P.Val < 0.01 , ])
dim(mDat)
m1 = log(mDat[rn,]+1)
dim(m1)

fGroups = fSamples
colnames(m1) = oExp.norm$phenotype
m1 = m1[,order(fGroups)]
fGroups = fGroups[order(fGroups)]

## download the cgraph library
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CGraphClust/master/CGraphClust.R'
download(url, 'CGraphClust.R')
m1 = log(m1+1)
# load the required packages
source('CGraphClust.R')
# delete the file after source
unlink('CGraphClust.R')
library('NMF')
# ignore step if stabalization not required
m1 = t(apply(m1, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(m1) = fGroups

# scale across rows
dfGenes = select(org.Mm.eg.db, rn, c('SYMBOL', 'GENENAME'), 'ENTREZID')
rownames(dfGenes) = dfGenes$ENTREZID
rownames(m1) = dfGenes[rownames(m1), 'SYMBOL']
mCounts = t(m1)
mCounts = scale(mCounts)
mCounts = t(mCounts)
# threshhold the values
mCounts[mCounts < -3] = -3
mCounts[mCounts > 3] = 3

# draw the heatmap  color='-RdBu:50'
aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=TRUE)

## save this object in the database for future use
## NOTE: don't run this segment of code again as object is already saved
## commenting for safety
# n = make.names(paste('list of negative binomial glm objects for organoids project with 3 coef and some shrinkage rds'))
# n2 = paste0('~/Data/MetaData/', n)
# save(lGlm, file=n2)
# 
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='list of negative binomial glm objects for organoids project fit using custom laplace function with 3 coefficients and some shrinkage')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

head(dfResults)
dfResults = dfResults[order(dfResults$adj.P.Val, decreasing = F),]

dfResults.plot = dfResults[dfResults$logFC > 0,]
dfResults.plot = dfResults.plot[1:5,]

fSamples = fCondition
rn = rownames(dfResults.plot)
dim(mDat)
m1 = log(mDat[rn,]+3)
fGroups = fSamples
colnames(m1) = oExp.norm$phenotype
m1 = m1[,order(fGroups)]
fGroups = fGroups[order(fGroups)]
head(m1)
dim(m1)

matplot(scale(t(m1)), type='l', lty=1, lwd=2, col=1:5, xlab='Samples', ylab='Z Scale Expression', main='Positive Fold Change - Bayes', xaxt='n',
        ylim=c(-2, 2))
axis(side = 1, labels = fGroups, at = 1:9, cex.axis=0.8, las=2)
abline(h = 0, lty=2, lwd=0.8)
gn = select(org.Mm.eg.db, keys = rownames(m1), keytype = 'ENTREZID', columns = 'SYMBOL')
legend('topleft', legend = gn$SYMBOL, fill=1:5)

### negative fold change
dfResults.plot = dfResults[dfResults$logFC < 0,]
dfResults.plot = dfResults.plot[1:5,]

fSamples = fCondition
rn = rownames(dfResults.plot)
dim(mDat)
m1 = log(mDat[rn,]+3)
fGroups = fSamples
colnames(m1) = oExp.norm$phenotype
m1 = m1[,order(fGroups)]
fGroups = fGroups[order(fGroups)]
head(m1)
dim(m1)

matplot(scale(t(m1)), type='l', lty=1, lwd=2, col=1:5, xlab='Samples', ylab='Z Scale Expression', main='Negative Fold Change - Bayes', xaxt='n',
        ylim=c(-2, 2))
axis(side = 1, labels = fGroups, at = 1:9, cex.axis=0.8, las=2)
abline(h = 0, lty=2, lwd=0.8)
gn = select(org.Mm.eg.db, keys = rownames(m1), keytype = 'ENTREZID', columns = 'SYMBOL')
legend('topright', legend = gn$SYMBOL, fill=1:5)

