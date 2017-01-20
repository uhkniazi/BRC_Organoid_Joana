# File: 10_de_analysis.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 20/01/2017
# Desc: DE analysis for the count data

## set variables and source libraries
source('header.R')
## libraries to load
library('RMySQL')

### space for internal functions

### end of functions


##### connect to mysql database to get count matrix
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 11) AND (MetaFile.comment like "%count%")')
dfCountMatrix = dbGetQuery(db, q)
dfCountMatrix
# load the sample meta data
q = paste0('select Sample.id as sid, Sample.group1 as description, Sample.group2 as phenotype, Sample.title, File.* from Sample, File
           where (Sample.idData = 11) AND (File.idSample = Sample.id AND File.type like "%duplicates removed%")')
dfSample.names = dbGetQuery(db, q)
dim(dfSample.names)
dfSample.names
# close connection after getting data
dbDisconnect(db)

## load the count matrix
n = paste0(dfCountMatrix$location, dfCountMatrix$name)
load(n)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)

# reorder the count matrix columns according to the order in samples table
i = match(dfSample.names$name, colnames(mCounts))
mCounts = mCounts[,i]
# sanity check
identical(dfSample.names$name, colnames(mCounts))

### factors for comparisons
# condition represents the biological condition
f2 = dfSample.names$phenotype
table(f2)
fCondition = factor(f2, levels=c('SI_compM', 'SI', 'CLP_SI'))

dfSample.names$fCondition = fCondition
str(dfSample.names)
# put all the data - matrix and annotation in expression set object
rownames(dfSample.names) = dfSample.names$title
colnames(mCounts) = dfSample.names$title
identical(rownames(dfSample.names), colnames(mCounts))

library(Biobase)
oExp = ExpressionSet(mCounts)
pData(oExp) = dfSample.names
dim(oExp)

library(DESeq2)
# get the size factors for the librar
sf = estimateSizeFactorsForMatrix(exprs(oExp))
oExp$LibrarySizeFactor = sf
# drop any unused levels
pData(oExp) = droplevels.data.frame(pData(oExp))

## create the DESeq object
dfDesign = data.frame(condition=oExp$fCondition, row.names = colnames(oExp))
dfDesign
## DE analysis
# call deseq2 constructor
oDseq = DESeqDataSetFromMatrix(exprs(oExp), dfDesign, design = ~ condition)
oDseq = DESeq(oDseq)
plotDispEsts(oDseq)

## get the normalized data
oExp.norm = oExp
exprs(oExp.norm) = counts(oDseq, normalized=T)

## save this object in the database for future use
## NOTE: don't run this segment of code again as object is already saved
# n = make.names(paste('Expression set object for normalised joana organoids mouse rds'))
# n2 = paste0('~/Data/MetaData/', n)
# save(oExp.norm, file=n2)
# 
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='Normalised Expression set object with comparison factor list from joana organoids mouse sequencing run with quality 10 duplicates removed and strict settings on trimmomatic')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

# get the results/contrasts for each comparison
levels(oExp.norm$fCondition)
oRes.SI.vs.SI_compM = results(oDseq, contrast = c('condition', 'SI', 'SI_compM'))
oRes.CLP_SI.vs.SI_compM = results(oDseq, contrast = c('condition', 'CLP_SI', 'SI_compM'))
oRes.CLP_SI.vs.SI = results(oDseq, contrast = c('condition', 'CLP_SI', 'SI'))

plotMA(oRes.SI.vs.SI_compM, main='SI vs SI_compM')
plotMA(oRes.CLP_SI.vs.SI_compM, main='CLP_SI vs SI_compM')
plotMA(oRes.CLP_SI.vs.SI, main='CLP_SI vs SI')

# plot histograms of adjusted p values
hist(oRes.CLP_SI.vs.SI$log2FoldChange)
hist(oRes.CLP_SI.vs.SI$padj)

hist(oRes.CLP_SI.vs.SI_compM$log2FoldChange)
hist(oRes.CLP_SI.vs.SI_compM$padj)

hist(oRes.SI.vs.SI_compM$log2FoldChange)
hist(oRes.SI.vs.SI_compM$padj)

table(oRes.CLP_SI.vs.SI$padj < 0.1)
table(oRes.CLP_SI.vs.SI_compM$padj < 0.1)
table(oRes.SI.vs.SI_compM$padj < 0.1)


# get results with significant p-values
dfD.vs.C = as.data.frame(oRes.D.vs.C[which(oRes.D.vs.C$padj < 0.1),])
dfH.vs.C = as.data.frame(oRes.H.vs.C[which(oRes.H.vs.C$padj < 0.1),])
dfD.vs.H = as.data.frame(oRes.D.vs.H[which(oRes.D.vs.H$padj < 0.1),])

nrow(dfD.vs.C)
nrow(dfD.vs.H)
nrow(dfH.vs.C)









# remove count matrix objects
rm(list=ls()[grepl('count', ls(), ignore.case = T)])
print(ls()[grepl('^f', ls(), ignore.case = T)])
rm(list=ls()[grepl('^f', ls(), ignore.case = T)])
print(ls()[grepl('^df', ls(), ignore.case = T)])
rm(list=ls()[grepl('^df', ls(), ignore.case = T)])


########## DE Analysis using glmer
library(multcomp)
#library(lme4)
library(glmmADMB)

## use the data matrix  a matrix of data
mDat = round(exprs(oExp),0)
str(mDat)

# remove low expression features
i = rowMeans(mDat)
i = which(i > 3)
mDat = mDat[i,]
str(mDat)

# fit glm to each feature
index = 1:nrow(mDat)

# check clock time
ptm = proc.time()

modelFunction = function(dat){
  df = data.frame(resp=mDat[dat,], cond.time=oExp$fCondition.t, patient=oExp$fTitle)
  return(tryCatch(glmmADMB::glmmadmb(resp ~ 0 + cond.time + (1 | patient), data=df, family = 'nbinom', link='log'), error=function(e) NULL))
}

# modelFunction = function(dat){
#   df = data.frame(resp=mDat[dat,], cond.time=oExp$fCondition.t, patient=oExp$fTitle)
#   #return(tryCatch(glmer.nb(resp ~ 0 + cond.time + (1 | patient), data=df), warning=function(w) NULL, error=function(e) NULL))
#   return(tryCatch(glmer.nb(resp ~ 0 + cond.time + (1 | patient), data=df), error=function(e) NULL))
# }

# lGlm = lapply(index, function(dat){
#   return(glmer.nb(mDat[dat,] ~ 0 + cond.time + (1 | patient)))
# })

lGlm = lapply(index, modelFunction)

names(lGlm) = rownames(mDat)

ptm.end = proc.time()

# # total time took
# > ptm.end - ptm
# user   system  elapsed 
# 6426.412  440.156 6758.919 ~ 2 hours

## save this object in the database for future use
## NOTE: don't run this segment of code again as object is already saved
## commenting for safety
# n = make.names(paste('list of glmer nb objects for keloids q10 rdup using glmmADMB rds'))
# n2 = paste0('~/Data/MetaData/', n)
# save(lGlm, file=n2)
# 
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=2, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='list of glmer negative binomial objects for keloids q10 rdup using glmmADMB - resp ~ 0 + cond.time + (1 | patient) for keloids S014 S021 and S032 sequencing runs with quality 10 duplicates removed')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

# remove the elements of list that are empty
lGlm.sub = lGlm[!sapply(lGlm, is.null)]

# extract a contrast at a time
mContrasts = rbind('Control:2 vs Control:1' = c(-1, 1, 0, 0),
                   'Keloid:1 vs Control:1' = c(-1, 0, 1, 0),
                   #'Keloid:2 vs Control:1' = c(-1, 0, 0, 1),
                   'Keloid:2 vs Keloid:1' = c(0, 0, -1, 1),
                   'Keloid:2 vs Control:2' = c(0, -1, 0, 1))


## perform contrasts tests
index = 1:length(lGlm.sub)

lContrast1 = lapply(index, function(dat){
  s = summary(glht(lGlm.sub[[dat]], t(mContrasts[1,])))
  ret = c(s$test$coefficients[1], s$test$pvalues[1])
  names(ret) = c('logfc', 'p.value')
  return(ret)
})

dfContrast1 = data.frame(do.call(rbind, lContrast1))
dfContrast1$p.adj = p.adjust(dfContrast1$p.value, method = 'BH')
rownames(dfContrast1) = names(lGlm.sub)

# second contrast
lContrast2 = mclapply(index, function(dat){
  s = summary(glht(lGlm.sub[[dat]], t(mContrasts[2,])))
  ret = c(s$test$coefficients[1], s$test$pvalues[1])
  names(ret) = c('logfc', 'p.value')
  return(ret)
})

dfContrast2 = data.frame(do.call(rbind, lContrast2))
dfContrast2$p.adj = p.adjust(dfContrast2$p.value, method = 'BH')
rownames(dfContrast2) = names(lGlm.sub)

# third contrast
lContrast3 = mclapply(index, function(dat){
  s = summary(glht(lGlm.sub[[dat]], t(mContrasts[3,])))
  ret = c(s$test$coefficients[1], s$test$pvalues[1])
  names(ret) = c('logfc', 'p.value')
  return(ret)
})

dfContrast3 = data.frame(do.call(rbind, lContrast3))
dfContrast3$p.adj = p.adjust(dfContrast3$p.value, method = 'BH')
rownames(dfContrast3) = names(lGlm.sub)

# fourth contrast
lContrast4 = mclapply(index, function(dat){
  s = summary(glht(lGlm.sub[[dat]], t(mContrasts[4,])))
  ret = c(s$test$coefficients[1], s$test$pvalues[1])
  names(ret) = c('logfc', 'p.value')
  return(ret)
})

dfContrast4 = data.frame(do.call(rbind, lContrast4))
dfContrast4$p.adj = p.adjust(dfContrast4$p.value, method = 'BH')
rownames(dfContrast4) = names(lGlm.sub)


## assign annotation to genes
library(org.Hs.eg.db)

df = select(org.Hs.eg.db, as.character(rownames(dfContrast1)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast1$SYMBOL = df$SYMBOL

df = select(org.Hs.eg.db, as.character(rownames(dfContrast2)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast2$SYMBOL = df$SYMBOL

df = select(org.Hs.eg.db, as.character(rownames(dfContrast3)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast3$SYMBOL = df$SYMBOL

df = select(org.Hs.eg.db, as.character(rownames(dfContrast4)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast4$SYMBOL = df$SYMBOL

## estimate dispersions and some quality plots
par(mfrow=c(2,2))
calculateDispersion = function(fm){
  n = length(resid(fm))
  sqrt(sum(c(as.numeric(resid(fm)), as.numeric(fm$U[[1]]))^2)/n)
}

iDispersion = sapply(lGlm.sub, calculateDispersion)

plotDispersion = function(dis, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$p.adj < p.cut] = 'red'
  plot(dis, dat$logfc, col=col, pch=20, main=title, xlab='Dispersion', ylab='logFC')
}

plotDispersion(iDispersion, dfContrast1, 0.01, 'Control:2 vs Control:1')
plotDispersion(iDispersion, dfContrast2, 0.1, 'Keloid:1 vs Control:1')
plotDispersion(iDispersion, dfContrast3, 0.01, 'Keloid:2 vs Keloid:1')
plotDispersion(iDispersion, dfContrast4, 0.1, 'Keloid:2 vs Control:2')

iMean = rowMeans(mDat[names(lGlm.sub),])

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$p.adj < p.cut] = 'red'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$logfc, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC')
}

par(mfrow=c(2,2))
plotMeanFC(log(iMean), dfContrast1, 0.01, 'Control:2 vs Control:1')
plotMeanFC(log(iMean), dfContrast2, 0.1, 'Keloid:1 vs Control:1')
plotMeanFC(log(iMean), dfContrast3, 0.01, 'Keloid:2 vs Keloid:1')
plotMeanFC(log(iMean), dfContrast4, 0.1, 'Keloid:2 vs Control:2')


### grouping of genes
dfContrast1.sub = na.omit(dfContrast1[dfContrast1$p.adj < 0.01,])
dfContrast2.sub = na.omit(dfContrast2[dfContrast2$p.adj < 0.1,])
dfContrast3.sub = na.omit(dfContrast3[dfContrast3$p.adj < 0.01,])
dfContrast4.sub = na.omit(dfContrast4[dfContrast4$p.adj < 0.1,])

cvCommonGenes = unique(c(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub), rownames(dfContrast4.sub)))
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=4)
mCommonGenes[,1] = cvCommonGenes %in% rownames(dfContrast1.sub)
mCommonGenes[,2] = cvCommonGenes %in% rownames(dfContrast2.sub)
mCommonGenes[,3] = cvCommonGenes %in% rownames(dfContrast3.sub)
mCommonGenes[,4] = cvCommonGenes %in% rownames(dfContrast4.sub)
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = gsub(' ', '', rownames(mContrasts))

## there appears to be a lot of genes coming up as significant
## perform a sensitivity analysis to see if certain samples are causing this problem
mCounts = mDat[cvCommonGenes,]
dim(mCounts)
str(mCounts)

## functions to perform sensitivity analysis
## posterior check
nbPosterior = function(x, prior=c(1/2, 1/2)){
  # calculate r i.e. alpha or size and p
  est = c('size'= abs(mean(x)^2/(var(x)-mean(x))), 'mu' = mean(x))
  est = c(est, est['size']/(est['size']+est['mu']))
  names(est)[3] = 'prob'
  # If the likelihood function for an observation x is negative binomial(r, p) and
  # p is distributed a priori as Beta(a, b) then the posterior distribution for p is
  # Beta(a + r, b + x). Note that this is the same as having observed r successes
  # and x failures with a binomial(r + x, p) likelihood. All that matters from a
  # Bayesian perspective is that r successes were observed and x failures.
  post = rbeta(1000, est['size']+prior[1], est['mu']+prior[2])
}

mSensitivityCheck = function(x){
  # create matrix to hold data
  mRet = matrix(NA, nrow=1000, ncol=length(x)+1)
  ## get posterior for full data
  mRet[,1] = nbPosterior(x)
  ## repeat with drop one observation
  for (i in 1:length(x)){
    mRet[,i+1] = nbPosterior(x[-i])
  }
  return(mRet)
}

mSensitivityCheckPvalues = function(x){
  m = mSensitivityCheck(x)
  #p.adjust(apply(m[,-1], 2, function(x) ks.test(m[,1], x)$p.value),method = 'bonf')
  apply(m[,-1], 2, function(x) ks.test(m[,1], x)$p.value)
}

lSenPvalues = lapply(seq_along(1:nrow(mCounts)), function(x){
  return(tryCatch(mSensitivityCheckPvalues(mCounts[x,]), error=function(e) NULL))
})

table(sapply(lSenPvalues, is.null))

mSenPvalues = do.call(rbind, lSenPvalues)
dim(mSenPvalues)
colnames(mSenPvalues) = colnames(mCounts)
rownames(mSenPvalues) = rownames(mCounts)
str(mSenPvalues)

## check the average to see if there is a pattern
cm = colMeans(mSenPvalues)
barplot(cm, las=2)
boxplot(mSenPvalues, las=2, ylab='p values distribution', main='bayesian sensitivity analysis for samples')
plot(cm, type='l', las=2, xlab='Samples', ylab='P values')

## is there a correlation between size factor and this
logit = function(p) log(p/(1-p))
sf = oExp$LibrarySizeFactor
plot(sf, logit(cm), main='Scatter plot for Sensitivity Analysis Average vs Library Size', xlab='Size Factor', ylab='Average P-Value for Sample')
fm = lm(logit(cm) ~ sf)
summary(fm)
iDrop = which.max(hatvalues(fm))
summary(update(fm, logit(cm[-iDrop]) ~ (sf[-iDrop])))
plot((sf[-iDrop]), logit(cm[-iDrop]), main='Scatter plot for Sensitivity Analysis Average vs Library Size', xlab='Size Factor', ylab='Average P-Value for Sample')

cor.test(sf[-iDrop], cm[-iDrop])

##### repeat the analysis for the specific genes after sensitivity analysis
mSenPvalues = t(apply(mSenPvalues, 1, p.adjust, 'bonf'))
fRepeat = apply(mSenPvalues, 1, function(x) any(x < 0.05))
table(fRepeat)
cvRepeat = rownames(mCounts[fRepeat,])
table(names(lGlm.sub) %in% cvRepeat)

# how many samples have an outlier value per gene
mRepeat = t(apply(mSenPvalues, 1, function(x) (x < 0.05)))
i = rowSums(mRepeat)
hist(i)
table(i > 5)

# fit glm to these genes
ptm = proc.time()

modelFunction2 = function(dat){
  df = data.frame(resp=mDat[dat,], cond.time=oExp$fCondition.t, patient=oExp$fTitle, pvalues=mSenPvalues[dat,])
  # drop the observations with low p-values
  i = which(df$pvalues < 0.05)
  df = df[-i,]
  df = droplevels.data.frame(df)
  return(tryCatch(glmmADMB::glmmadmb(resp ~ 0 + cond.time + (1 | patient), data=df, family = 'nbinom', link='log'), error=function(e) NULL))
}

lGlm.rep = lapply(cvRepeat, modelFunction2)

names(lGlm.rep) = cvRepeat

ptm.end = proc.time()

## save this object in the database for future use
## NOTE: don't run this segment of code again as object is already saved
## commenting for safety
# n = make.names(paste('list of glmer nb objects for keloids q10 rdup using glmmADMB after doing sensitivity analysis rds'))
# n2 = paste0('~/Data/MetaData/', n)
# save(lGlm.rep, file=n2)
# 
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=2, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='list of glmer negative binomial objects for keloids q10 rdup using glmmADMB - resp ~ 0 + cond.time + (1 | patient) for keloids S014 S021 and S032 sequencing runs with quality 10 duplicates removed, repeated after doing sensitivity analysis')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)


## remove any null elements
table(sapply(lGlm.rep, is.null))
f = sapply(lGlm.rep, is.null)
lGlm.rep[f] = NULL

i = match(names(lGlm.rep), names(lGlm.sub))
head(names(lGlm.sub[i]),3)
head(names(lGlm.rep), 3)
lGlm.sub[i] = lGlm.rep

##################################################################
### repeat the contrasts 
# extract a contrast at a time
mContrasts = rbind('Control:2 vs Control:1' = c(-1, 1, 0, 0),
                   'Keloid:1 vs Control:1' = c(-1, 0, 1, 0),
                   #'Keloid:2 vs Control:1' = c(-1, 0, 0, 1),
                   'Keloid:2 vs Keloid:1' = c(0, 0, -1, 1),
                   'Keloid:2 vs Control:2' = c(0, -1, 0, 1))


## perform contrasts tests
index = 1:length(lGlm.sub)

lContrast1 = lapply(index, function(dat){
  tryCatch({
    s = summary(glht(lGlm.sub[[dat]], t(mContrasts[1,])))
    ret = c(s$test$coefficients[1], s$test$pvalues[1])
    names(ret) = c('logfc', 'p.value')
    return(ret)
  }, error=function(e) {
    ret = c('logfc'=NA, 'p.value'=NA)
    return(ret)
  })
})

dfContrast1 = data.frame(do.call(rbind, lContrast1))
dfContrast1$p.adj = p.adjust(dfContrast1$p.value, method = 'BH')
rownames(dfContrast1) = names(lGlm.sub)

# second contrast
lContrast2 = mclapply(index, function(dat){
  tryCatch({
    s = summary(glht(lGlm.sub[[dat]], t(mContrasts[2,])))
    ret = c(s$test$coefficients[1], s$test$pvalues[1])
    names(ret) = c('logfc', 'p.value')
    return(ret)
  }, error=function(e) {
    ret = c('logfc'=NA, 'p.value'=NA)
    return(ret)
  })
})

dfContrast2 = data.frame(do.call(rbind, lContrast2))
dfContrast2$p.adj = p.adjust(dfContrast2$p.value, method = 'BH')
rownames(dfContrast2) = names(lGlm.sub)

# third contrast
lContrast3 = mclapply(index, function(dat){
  tryCatch({
    s = summary(glht(lGlm.sub[[dat]], t(mContrasts[3,])))
    ret = c(s$test$coefficients[1], s$test$pvalues[1])
    names(ret) = c('logfc', 'p.value')
    return(ret)
  }, error=function(e) {
    ret = c('logfc'=NA, 'p.value'=NA)
    return(ret)
  })
})

dfContrast3 = data.frame(do.call(rbind, lContrast3))
dfContrast3$p.adj = p.adjust(dfContrast3$p.value, method = 'BH')
rownames(dfContrast3) = names(lGlm.sub)

# fourth contrast
lContrast4 = mclapply(index, function(dat){
  tryCatch({
    s = summary(glht(lGlm.sub[[dat]], t(mContrasts[4,])))
    ret = c(s$test$coefficients[1], s$test$pvalues[1])
    names(ret) = c('logfc', 'p.value')
    return(ret)
  }, error=function(e) {
    ret = c('logfc'=NA, 'p.value'=NA)
    return(ret)
  })
})

dfContrast4 = data.frame(do.call(rbind, lContrast4))
dfContrast4$p.adj = p.adjust(dfContrast4$p.value, method = 'BH')
rownames(dfContrast4) = names(lGlm.sub)

## assign annotation to genes
library(org.Hs.eg.db)

df = select(org.Hs.eg.db, as.character(rownames(dfContrast1)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast1$SYMBOL = df$SYMBOL
dfContrast1$GENENAME = df$GENENAME

df = select(org.Hs.eg.db, as.character(rownames(dfContrast2)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast2$SYMBOL = df$SYMBOL
dfContrast2$GENENAME = df$GENENAME

df = select(org.Hs.eg.db, as.character(rownames(dfContrast3)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast3$SYMBOL = df$SYMBOL
dfContrast3$GENENAME = df$GENENAME

df = select(org.Hs.eg.db, as.character(rownames(dfContrast4)), c('SYMBOL', 'GENENAME'), 'ENTREZID')
table(duplicated(df$ENTREZID))
dfContrast4$SYMBOL = df$SYMBOL
dfContrast4$GENENAME = df$GENENAME

## estimate dispersions and some quality plots
par(mfrow=c(2,2))
calculateDispersion = function(fm){
  n = length(resid(fm))
  sqrt(sum(c(as.numeric(resid(fm)), as.numeric(fm$U[[1]]))^2)/n)
}

iDispersion = sapply(lGlm.sub, calculateDispersion)

plotDispersion = function(dis, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$p.adj < p.cut] = 'red'
  plot(dis, dat$logfc, col=col, pch=20, main=title, xlab='Dispersion', ylab='logFC')
}

plotDispersion(iDispersion, dfContrast1, 0.01, 'Control:2 vs Control:1')
plotDispersion(iDispersion, dfContrast2, 0.1, 'Keloid:1 vs Control:1')
plotDispersion(iDispersion, dfContrast3, 0.01, 'Keloid:2 vs Keloid:1')
plotDispersion(iDispersion, dfContrast4, 0.1, 'Keloid:2 vs Control:2')

iMean = rowMeans(mDat[names(lGlm.sub),])

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$p.adj < p.cut] = 'red'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$logfc, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC')
}

par(mfrow=c(2,2))
plotMeanFC(log(iMean), dfContrast1, 0.01, 'Control:2 vs Control:1')
plotMeanFC(log(iMean), dfContrast2, 0.1, 'Keloid:1 vs Control:1')
plotMeanFC(log(iMean), dfContrast3, 0.01, 'Keloid:2 vs Keloid:1')
plotMeanFC(log(iMean), dfContrast4, 0.1, 'Keloid:2 vs Control:2')

# add dispersion parameter to ecah gene
dfContrast1$Dispersion = iDispersion[rownames(dfContrast1)]
dfContrast2$Dispersion = iDispersion[rownames(dfContrast2)]
dfContrast3$Dispersion = iDispersion[rownames(dfContrast3)]
dfContrast4$Dispersion = iDispersion[rownames(dfContrast4)]

# save the gene lists
write.csv(na.omit(dfContrast1), file='Results/Control:2vsControl:1.xls')
write.csv(na.omit(dfContrast2), file='Results/Keloid:1vsControl:1.xls')
write.csv(na.omit(dfContrast3), file='Results/Keloid:2vsKeloid:1.xls')
write.csv(na.omit(dfContrast4), file='Results/Keloid:2vsControl:2.xls')

### grouping of genes
dfContrast1.sub = na.omit(dfContrast1[dfContrast1$p.adj < 0.01 & dfContrast1$Dispersion > 0.4,])
dfContrast2.sub = na.omit(dfContrast2[dfContrast2$p.adj < 0.1 & dfContrast2$Dispersion > 0.4,])
dfContrast3.sub = na.omit(dfContrast3[dfContrast3$p.adj < 0.01 & dfContrast3$Dispersion > 0.4,])
dfContrast4.sub = na.omit(dfContrast4[dfContrast4$p.adj < 0.1 & dfContrast4$Dispersion > 0.4,])

cvCommonGenes = unique(c(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub), rownames(dfContrast4.sub)))
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=4)
mCommonGenes[,1] = cvCommonGenes %in% rownames(dfContrast1.sub)
mCommonGenes[,2] = cvCommonGenes %in% rownames(dfContrast2.sub)
mCommonGenes[,3] = cvCommonGenes %in% rownames(dfContrast3.sub)
mCommonGenes[,4] = cvCommonGenes %in% rownames(dfContrast4.sub)
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = gsub(' ', '', rownames(mContrasts))


############

#### analysis by grouping genes
# create groups in the data based on 4^2-1 combinations
mCommonGenes.grp = mCommonGenes
set.seed(123)
dm = dist(mCommonGenes.grp, method='binary')
hc = hclust(dm)
plot(hc)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mCommonGenes.grp = cbind(mCommonGenes.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mCommonGenes.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## make venn diagrams of comparisons
library(VennDiagram)
par(p.old)
# create a list for overlaps
lVenn = list(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub), rownames(dfContrast4.sub))
names(lVenn) = c('C2vsC1', 'K1vsC1', 'K2vsK1', 'K2vsC2')
# calculate overlaps
#lVenn.overlap = calculate.overlap(lVenn)
venn.diagram(lVenn, filename = 'Results/venn_all_contrasts.tif')
venn.diagram(lVenn[c(1,3)], filename = 'Results/venn_time_contrasts.tif')

## save the genes in the overlaps of interest
# genes common between C2vsC1 and K2vsK1
rn = apply(mCommonGenes, 1, function(x) all(x[c(1, 3)] == c(T, T)))
rn = names(rn[rn])
length(rn)

df.rn = select(org.Hs.eg.db, keys = rn, columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
dim(df.rn)
str(df.rn)
# write csv to look at gene list
write.csv(df.rn, file=paste('Temp/', 'common.between.c2vsc1.and.k2vsk1', '.xls', sep=''))

# genes unique to C2vsC1 VS K2vsK1
rn = apply(mCommonGenes, 1, function(x) all(x[c(1, 3)] == c(T, F)))
rn = names(rn[rn])
length(rn)

df.rn = select(org.Hs.eg.db, keys = rn, columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
dim(df.rn)
str(df.rn)
# write csv to look at gene list
write.csv(df.rn, file=paste('Temp/', 'unique.to.c2vsc1.VS.k2vsk1', '.xls', sep=''))

# genes unique to K2vsK1 VS C2vsC1
rn = apply(mCommonGenes, 1, function(x) all(x[c(1, 3)] == c(F, T)))
rn = names(rn[rn])
length(rn)

df.rn = select(org.Hs.eg.db, keys = rn, columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
dim(df.rn)
str(df.rn)
# write csv to look at gene list
write.csv(df.rn, file=paste('Temp/', 'unique.to.k2vsk1.VS.c2vsc1', '.xls', sep=''))

# genes present in other 2 remaining contrasts
df.rn = select(org.Hs.eg.db, keys = rownames(dfContrast2.sub), columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
dim(df.rn)
str(df.rn)
write.csv(df.rn, file='Temp/K1vsC1.xls')

df.rn = select(org.Hs.eg.db, keys = rownames(dfContrast4.sub), columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
dim(df.rn)
str(df.rn)
write.csv(df.rn, file='Temp/K2vsC2.xls')

## insert this gene list at innate db to see pathways overrepresented
## import innate db result
dfCommonAcrossTime = read.csv(file.choose(), header = T, sep='\t', stringsAsFactors = F)
dfCommonAcrossTime = dfCommonAcrossTime[,-10]
# sort the pathway table on adjusted p-values
dfCommonAcrossTime = dfCommonAcrossTime[order(dfCommonAcrossTime$Pathway.p.value),]

## make a barplot
plot.bar = function(ivBar, title='', ...){
  p.old = par(mar=c(8,4,2,2)+0.1)
  l = barplot(ivBar, beside=T, xaxt='n', main=title, ...)
  axis(side = 1, l[,1], labels=F)
  text(l[,1], y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=names(ivBar), srt=45, adj=1, xpd=TRUE, cex=0.6)
  par(p.old)
}

i = -1*log10(dfCommonAcrossTime$Pathway.p.value[1:7])
names(i) = dfCommonAcrossTime$Pathway.Name[1:7]
# make on of the names shorter
names(i)[2] = 'Extracellular matrix degradation'
pdf('Temp/injury_response.pdf')
plot.bar(i, title='Common Response to Injury', ylab='-log10 PValue')
dev.off(dev.cur())

## recycling the code above 
## import innate db result for unique to keloids across time
dfCommonAcrossTime = read.csv(file.choose(), header = T, sep='\t', stringsAsFactors = F)
dfCommonAcrossTime = dfCommonAcrossTime[,-10]
# sort the pathway table on adjusted p-values
dfCommonAcrossTime = dfCommonAcrossTime[order(dfCommonAcrossTime$Pathway.p.value),]

## make a barplot
plot.bar = function(ivBar, title='', ...){
  p.old = par(mar=c(8,4,2,2)+0.1)
  l = barplot(ivBar, beside=T, xaxt='n', main=title, ...)
  axis(side = 1, l[,1], labels=F)
  text(l[,1], y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=names(ivBar), srt=45, adj=1, xpd=TRUE, cex=0.6)
  par(p.old)
}

i = -1*log10(dfCommonAcrossTime$Pathway.p.value[1:7])
names(i) = dfCommonAcrossTime$Pathway.Name[1:7]
# make on of the names shorter
names(i)[7] = 'Molecular Transport'
pdf('Temp/unique_keloids_injury_response.pdf')
plot.bar(i, title='Unique in Keloids Response to Injury', ylab='-log10 PValue')
dev.off(dev.cur())


### make some heatmaps
## all overexpressed genes if interested in
fSamples = oExp$fCondition.t
i = 1:nrow(mCommonGenes)

m1 = as.matrix(mCommonGenes[i,])
m1 = mDat[rownames(m1),]

fGroups = fSamples
colnames(m1) = fGroups
m1 = m1[,order(fGroups)]
fGroups = fGroups[order(fGroups)]

## download the cgraph library
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CGraphClust/master/CGraphClust.R'
download(url, 'CGraphClust.R')

# load the required packages
source('CGraphClust.R')
# delete the file after source
unlink('CGraphClust.R')
library('NMF')
m1 = log(m1+0.5)
# ignore step if stabalization not required
m1 = t(apply(m1, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(m1) = fGroups

# scale across rows
dfGenes = select(org.Hs.eg.db, cvCommonGenes, c('SYMBOL', 'GENENAME'), 'ENTREZID')
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


####################################################################################################
########### pathway analysis using CGraph library
# uniprot annotation for data
cvGenes = rownames(mCommonGenes)
dfGenes = select(org.Hs.eg.db, keys = cvGenes, columns = c('ENTREZID', 'SYMBOL', 'GENENAME', 'UNIPROT'), keytype = 'ENTREZID')
dfGenes = na.omit(dfGenes)

# get reactome data
url = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
dir.create('Data_external', showWarnings = F)
csReactomeFile = 'Data_external/UniProt2Reactome_All_Levels.txt'
# download the reactome file if it doesnt exist
if (!file.exists(csReactomeFile)) download(url, csReactomeFile)
# map reactome pathways
dfReactome = read.csv(csReactomeFile, header = F, stringsAsFactors=F, sep='\t')
x = gsub('.+\\w+-\\w+-(\\d+)', replacement = '\\1', x = dfReactome$V2, perl = T)
dfReactome$V2 = x
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfGenes$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfGenes$UNIPROT)
dfReactome.sub$ENTREZID = dfGenes$ENTREZID[i]
dfGraph = dfReactome.sub[,c('ENTREZID', 'V2')]
dfGraph = na.omit(dfGraph)

# get expression data
mCounts = mDat[unique(dfGenes$ENTREZID),]
# ## optional adjusting the data for repeated measurements
# temp = sapply(rownames(mCounts), function(x){
#   return(fitted(lGlm.sub[[x]]))
# })
# mCounts = t(mCounts)

fGroups = fSamples
names(fGroups) = oExp$fTitle
colnames(mCounts) = fGroups
# reorder on grouping factor
mCounts = mCounts[,order(fGroups)]
fGroups = fGroups[order(fGroups)]
mCounts = t(mCounts)

# select genes that have a reactome term attached
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[,n]
print(paste('Total number of genes with Reactome terms', length(n)))
levels(fGroups)

# [1] "Total number of genes with Reactome terms 910"
# > levels(fGroups)
# [1] "Control:1" "Control:2" "Keloid:1"  "Keloid:2" 
# > 

# create a correlation matrix to decide cor cutoff
mCor = cor(log(mCounts+0.5))

# check distribution 
hist(mCor, prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# create the graph cluster object
# using absolute correlation vs actual values lead to different clusters
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.6, bSuppressPlots = F)

## general graph structure
set.seed(1)
plot.final.graph(oGr)
ecount(getFinalGraph(oGr))
vcount(getFinalGraph(oGr))
## community structure
# plot the main communities and centrality graphs
ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 2)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, mark.groups=NULL, edge.color='lightgrey')
set.seed(1)
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 2)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, edge.color='darkgrey')

# get community sizes
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
iSizes = sort(table(dfCluster$cluster))
# remove communities smaller than 5 members or choose a size of your liking
i = which(iSizes <= 5)
if (length(i) > 0) {
  cVertRem = as.character(dfCluster[dfCluster$cluster %in% names(i),'gene'])
  iVertKeep = which(!(V(getFinalGraph(oGr))$name %in% cVertRem))
  oGr = CGraphClust.recalibrate(oGr, iVertKeep)
}

## centrality diagnostics
## centrality parameters should not be correlated significantly and the location of the central
## genes can be visualized
# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
# plot the genes or vertex sizes by fold change
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 2)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
par(p.old)

## the diagnostic plots show the distribution of the centrality parameters
# these diagnostics plots should be looked at in combination with the centrality graphs
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)

## top vertices based on centrality scores
## get a table of top vertices 
dfTopGenes.cent = dfGetTopVertices(oGr, iQuantile = 0.90)
rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
# assign metadata annotation to these genes and clusters
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
dfTopGenes.cent = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)

####### NOTE: This section of code is very slow, use only if you need data from genbank
# loop and get the data from genbank
n = rep(NA, length=nrow(dfTopGenes.cent))
names(n) = as.character(dfTopGenes.cent$VertexID)
for (i in seq_along(n)){
  n[i] = f_csGetGeneSummaryFromGenbank(names(n)[i])
  # this wait time is required to stop sending queries to ncbi server very quickly
  Sys.sleep(time = 3)
}
cvSum.2 = as.character(dfTopGenes.cent$VertexID)
dfTopGenes.cent$Summary = n[cvSum.2]
####### Section ends

write.csv(dfTopGenes.cent, file='Temp/Top_Centrality_Genes.xls')

## if we want to look at the expression profiles of the top genes
# plot a heatmap of these top genes
library(NMF)
m1 = log(mCounts[,as.character(dfTopGenes.cent$VertexID)]+0.5)
m1 = scale(m1)
m1 = t(m1)
# threshhold the values
m1[m1 < -3] = -3
m1[m1 > 3] = 3
rownames(m1) = as.character(dfTopGenes.cent$SYMBOL)
# draw the heatmap  color='-RdBu:50'
aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)

m1 = log(mCounts[,as.character(dfTopGenes.cent$VertexID)]+0.5)
m1 = apply(m1, 2, f_ivStabilizeData, fGroups)
rownames(m1) = fGroups
m1 = scale(m1)
m1 = t(m1)
# threshhold the values
m1[m1 < -3] = -3
m1[m1 > 3] = 3
rownames(m1) = as.character(dfTopGenes.cent$SYMBOL)
# draw the heatmap  color='-RdBu:50'
aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)

## in addition to heatmaps the graphs can be plotted
# plot a graph of these top genes
# plot for each contrast i.e. base line vs other level
lev = levels(fGroups)[-1]
m = mCounts
m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = induced_subgraph(getFinalGraph(oGr), vids = as.character(dfTopGenes.cent$VertexID))
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=5)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.2, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
       main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}


### Looking at the largest clique can be informative in the graph
# plot the graph with location of the clique highlighted
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 2)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

# plot the largest clique at each grouping contrast
lev = levels(fGroups)[-1]
m = mCounts
#m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
#rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = induced_subgraph(getFinalGraph(oGr), vids = unlist(getLargestCliques(oGr)))
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=2)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, layout=layout_with_fr, main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}

## instead of looking at individual genes we can look at clusters
## we can look at the problem from the other direction and look at clusters instead of genes
# some sample plots
# mean expression of groups in every cluster
par(p.old)
plot.mean.expressions(oGr, log(t(mCounts)+0.5), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster', cex.axis=0.7)
# only significant clusters
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, log(t(mCounts)+0.5), fGroups, main='Significant Clusters', lwd=1, bStabalize = T, cex.axis=0.7)#, p.cut=0.05 )
# principal component plots
pr.out = plot.components(oGr, log(t(mCounts)+0.5), fGroups, bStabalize = T)#, p.cut=0.05)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
# plot summary heatmaps
# marginal expression level in each cluster
plot.heatmap.significant.clusters(oGr, log(t(mCounts)+0.5), fGroups, bStabalize = F)#, p.cut=0.05)
#plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T, p.cut=0.05)
# plot variance of cluster
m = getSignificantClusters(oGr, log(t(mCounts)+0.5), fGroups)$clusters
m = getClusterMarginal(oGr, log(t(mCounts)+0.5))
csClust = rownames(m)
length(csClust)
pdf('Temp/cluster_variance.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust), las=2)
for (i in seq_along(csClust)){
  temp = t(as.matrix(m[csClust[i],]))
  rownames(temp) = csClust[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
}
dev.off(dev.cur())


# print the labels of the clusters from reactome table
i = which(dfReactome.sub$V2 %in% csClust)
dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
rownames(dfCluster.name) = NULL
dfCluster.name

#### plot a graph of clusters 
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
write.csv(dfCluster, file='Temp/dfCluster_members.xls')
# how many genes in each cluster
data.frame(sort(table(dfCluster$cluster)))
#csClust = rownames(m$clusters)
csClust = as.character(unique(dfCluster$cluster))

mMarginal = getClusterMarginal(oGr, log(t(mCounts)+0.5))

### create lattice plots for each cluster and group i.e. keloids and normals/controls
library(lattice)
dfData = data.frame(scale(t(mMarginal)))
#colnames(dfData) = gsub('X', '', colnames(dfData))
dfData$groups = fGroups
dfData$patient = factor(names(fGroups))
# extract patient conditions from the fGroups
fTime = factor(gsub('^Control:|Keloid:', '', as.character(fGroups)))
fCondition = factor(gsub(':1|:2', '', as.character(fGroups)))
dfData$time = fTime
dfData$condition = fCondition

str(dfData)
## stack the dataframe expression values
dfStack = stack(dfData)
str(dfStack)
dfStack$time = dfData$time
dfStack$condition = dfData$condition
dfStack$patient = dfData$patient

xyplot(values ~ time | ind, data=dfStack, type=c('n','smooth'), groups=condition, auto.key=T)

## xyplots with population average
mMarginal = getClusterMarginal(oGr, log(t(mCounts)+0.5))
mMarginal = apply(mMarginal, 1, function(x) tapply((x), fGroups, mean))
mMarginal = scale(mMarginal)
dfData = data.frame(mMarginal)
cn = colnames(mMarginal)
str(dfData)
rownames(dfCluster.name) = dfCluster.name$V2
dfCluster.name[cn,]
cn = c('Transmembrane transport', 'Amino A Met', 'Organelle biogenesis', 'Neutrophil degran/GTPase',
       'Epigenetics', 'ECM Degradation', 'Semaphorin interactions', 'Membrane Transport', 'ECM Proteoglycans',
       'Class I MHC', 'GPCR binding', 'Muscle contraction', 'Sphingolipid metabolism', 'Keratinization',
       'Wnt Signaling', 'Deubiquitination', 'Neuronal System', 'Cell communication', 'Stress response', 'Cytokine Signaling',
       'FA metabolism', 'rRNA processing', 'Lipid metabolism', 'p53 regulation', 'Olf Rec Path')
colnames(dfData) = cn
temp = dfCluster.name[colnames(mMarginal),]
temp$cn = cn
temp
#colnames(dfData) = gsub('X', '', colnames(dfData))
dfData$groups = factor(rownames(dfData))
dfData$time = factor(gsub('^Control:|Keloid:', '', as.character(dfData$groups)))
dfData$condition = factor(gsub(':1|:2', '', as.character(dfData$groups)))

dfStack = stack(dfData)
dfStack$time = dfData$time
dfStack$condition = dfData$condition
pdf('Temp/module_summary.pdf')
xyplot(values ~ time | ind, data=dfStack, type=c('b'), groups=condition, auto.key=T, par.strip.text=list(cex=0.6),
       scales=list(cex=0.6), xlab='Time', ylab='Scaled Module Average')#, layout=c(4,7))
dev.off(dev.cur())
xyplot(values ~ groups | ind, data=dfStack, type='o')