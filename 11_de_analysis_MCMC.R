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
  if ( mode <=0 ) stop("mode must be > 0")
  if ( sd <=0 ) stop("sd must be > 0")
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
iCutoff = 2
fKeep = rowMeans(temp) > iCutoff
table(fKeep)
oExp.norm = oExp.norm[fKeep,]
dim(oExp.norm)

## prepare data to fit model
fCondition = as.character(oExp.norm$fCondition)
fCondition[fCondition != 'CLP_SI'] = 'Media'
fCondition = factor(fCondition, levels = c('Media', 'CLP_SI'))
mDat = exprs(oExp.norm)


#### define function to optimize
library(LearnBayes)
library(MASS)
## define log posterior
mylogpost = function(theta, data){
  size = theta['size']
  beta0 = theta['beta0']
  beta1 = theta['beta1']
  # hyperparameters for the hierarchichal standard deviation parameter
  betaSigma = exp(theta['betaSigma']) # maybe add a one to prevent zero variance
  ##if (betaSigma <= 0) return(-Inf)
  # hyper-hyperparameters for the hierarchichal standard deviation parameter
  ivBetaParam = data$betaParam
  x = data$pred
  y = data$resp
  ## create model matrix to get fitted value
  mModMatrix = model.matrix(y ~ x)
  ## likelihood function
  lf = function(dat, pred){
    return(log(dnbinom(dat, size = size, mu = pred)))
  }
  mCoef = matrix(c(beta0, beta1), nrow = 2)
  mu = exp(mModMatrix %*% mCoef)
  val = sum(lf(y, mu))
  ret = val + dnorm(beta0, 0, betaSigma, log=T) + dnorm(beta1, 0, betaSigma, log=T) + dunif(size, 1, 1e+3, log=T) +
    dgamma(betaSigma, shape = ivBetaParam['shape'], rate = ivBetaParam['rate'], log=T)
  return(ret)
}

getOptimizedPValue = function(obj){
  se = sqrt(diag(obj$var))[c('beta0', 'beta1')]
  m = obj$mode[c('beta0', 'beta1')]
  pnorm(-abs(m/se))*2
}

getOptimizedSummary = function(obj){
  se = sqrt(diag(obj$var))['beta1']
  m = obj$mode['beta1']
  p = pnorm(-abs(m/se))*2
  ret  = c('Coef' = round(m, 3), 'SE'=round(se, 3), 'P-Value'=signif(p, 3))
  names(ret) = c('Coef', 'SE', 'P-Value')
  return(ret)
}

############################### fit the model on multiple genes
modelFunction = function(dat){
  ## prepare data before fitting
  dfData = data.frame(resp=as.integer(mDat[dat,]), pred=fCondition)
  temp = fitdistr(dfData$resp, 'negative binomial')$estimate['size']
  names(temp) = NULL
  # set starting values for optimizer and 
  # set parameters for optimizer
  lData = list(resp=dfData$resp, pred=dfData$pred)
  lData$betaParam = unlist(gammaShRaFromModeSD(mode = log(sd(dfData$resp+0.5)/2), sd = log(2*sd(dfData$resp+0.5))))
  start = c('size'=temp, 'beta0'=log(mean(dfData$resp)), 'beta1'=0, 
            'betaSigma'=log(mean(rgamma(1000, shape = lData$betaParam['shape'], rate = lData$betaParam['rate']))))
  
  laplace(mylogpost, start, lData)
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
write.csv(dfPlot, file='Results/DEAnalysis_CLP_SI.vs.Media_2.xls')

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

hist(dfResults$logFC)
hist(dfResults$P.Value)
hist(dfResults$adj.P.Val)

table(dfResults$adj.P.Val < 0.1)

### make some heatmaps
fSamples = fCondition
rn = rownames(dfResults[dfResults$adj.P.Val < 0.1 , ])
dim(mDat)
m1 = mDat[rn,]
dim(m1)

fGroups = fSamples
colnames(m1) = fGroups
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
         annColors=NA, Colv=NA)


