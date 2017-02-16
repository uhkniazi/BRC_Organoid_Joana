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
#oRes.SI.vs.SI_compM = results(oDseq, contrast = c('condition', 'SI', 'SI_compM'))
oRes.CLP_SI.vs.SI_compM = results(oDseq, contrast = c('condition', 'CLP_SI', 'SI_compM'))
#oRes.CLP_SI.vs.SI = results(oDseq, contrast = c('condition', 'CLP_SI', 'SI'))

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$p.adj < p.cut] = 'red'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$logfc, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC', ylim=c(-2, 2), cex=0.6)
}

# perform independent filtering 
iCutoff = 2
fKeep = rowMeans(counts(oDseq, normalized=T)) > iCutoff
table(fKeep)

dim(oRes.CLP_SI.vs.SI_compM)
oRes.CLP_SI.vs.SI_compM = oRes.CLP_SI.vs.SI_compM[fKeep,]
dim(oRes.CLP_SI.vs.SI_compM)
oRes.CLP_SI.vs.SI_compM$padj = p.adjust(oRes.CLP_SI.vs.SI_compM$pvalue, method = 'BH')

dfGenes = data.frame(p.adj=oRes.CLP_SI.vs.SI_compM$padj, logfc=log(2^oRes.CLP_SI.vs.SI_compM$log2FoldChange))
iMean = oRes.CLP_SI.vs.SI_compM$baseMean

plotMeanFC(log(iMean), dfGenes, 0.1, 'MA Plot')

plotMA(oRes.CLP_SI.vs.SI_compM, main='CLP_SI vs SI_compM')


# plot histograms of adjusted p values

hist(oRes.CLP_SI.vs.SI_compM$log2FoldChange)
hist(oRes.CLP_SI.vs.SI_compM$padj)

table(oRes.CLP_SI.vs.SI_compM$padj < 0.1)

# get results with significant p-values
dfCLP_SI.vs.SI_compM = as.data.frame(oRes.CLP_SI.vs.SI_compM[which(oRes.CLP_SI.vs.SI_compM$padj < 0.1),])
dim(dfCLP_SI.vs.SI_compM)

### create some volcano plots with annotations
## choose the comparison for plotting
library(org.Mm.eg.db)
# add annotation to the data set after selecting comparison
res = as.data.frame(oRes.CLP_SI.vs.SI_compM)

rn = rownames(res)
df = select(org.Mm.eg.db, as.character(rn), c('SYMBOL'), 'ENTREZID')
head(df); head(res)
df = df[!duplicated(df$ENTREZID),]
rownames(df) = df$ENTREZID
dfPlot = res
dfPlot = cbind(dfPlot[rn,], df[rn,])
dfPlot = na.omit(dfPlot)

## write csv file
write.csv(dfPlot, file='Results/DEAnalysis_CLP_SI.vs.Media_with_DESeq.xls')

dfGenes = data.frame(P.Value=dfPlot$pvalue, logFC=dfPlot$log2FoldChange, adj.P.Val = dfPlot$padj, SYMBOL=dfPlot$SYMBOL)

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
  y.cut = quantile(p.val[c], probs=0)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}

f_plotVolcano(dfGenes, 'CLP_SI vs SI')

