# File: 11.2_de_analysis_edgeR.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 20/02/2017
# Desc: DE analysis for the count data using EdgeR

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

library(edgeR)
# drop any unused levels
pData(oExp) = droplevels.data.frame(pData(oExp))

fCondition = oExp$fCondition
table(fCondition)
levels(fCondition)
# # drop the si
# f = which(fCondition == 'SI')
# f
# oExp = oExp[,-f]
dim(oExp)
# # keep the others
# fCondition = as.character(oExp$fCondition)
# table(fCondition)
# fCondition[fCondition != 'CLP_SI'] = 'Media'
# fCondition = factor(fCondition, levels = c('Media', 'CLP_SI'))
# oExp$fCondition = fCondition

#fCondition=oExp$fCondition
oEdge = DGEList(exprs(oExp), group=fCondition)
# mCpm = cpm(oEdge)
# fKeep = rowMeans(mCpm) > 1
# table(fKeep)
# oEdge = DGEList(exprs(oExp)[fKeep,], group=fCondition)
## remove low counts
oEdge = oEdge[rowSums(1e+06 * oEdge$counts/expandAsMatrix(oEdge$samples$lib.size, dim(oEdge)) > 1) >= 3, ]
dim(oEdge)
oEdge = calcNormFactors(oEdge)
levels(fCondition)
mModMatrix = model.matrix(~ fCondition)
colnames(mModMatrix) = levels(fCondition)
#oEdge = estimateDisp(oEdge, mModMatrix, )
oEdge = estimateCommonDisp(oEdge)
oEdge = estimateTagwiseDisp(oEdge, prior.n = 10)
oEdge.glm = glmFit(oEdge, mModMatrix)
oLrt = glmLRT(oEdge.glm, coef=3)

dfResults = oLrt$table
# fKeep = dfResults$logCPM > 1
# table(fKeep)
# dfResults = dfResults[fKeep,]
table(dfResults$PValue < 0.01)
dfResults$padj = p.adjust(dfResults$PValue, method='BH')
table(dfResults$padj < 0.1)
# cCont = sprintf("%s-%s", "CLP_SI", "SI_compM") # Samples to be compared
# oContrast = makeContrasts(contrasts=cCont, levels=mModMatrix) # Make contrast
# oRes = glmLRT(oEdge, contrast = oContrast)
# oRes$table$QValue = p.adjust(oRes$table$PValue, method = "BH")
# 
# dfRes = oRes$table

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$p.adj < p.cut] = 'red'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$logfc, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC', ylim=c(-2, 2), cex=0.6)
}

dfGenes = data.frame(p.adj=dfResults$padj, logfc=dfResults$logFC)
iMean = dfResults$logCPM

plotMeanFC(iMean, dfGenes, 0.1, 'MA Plot')

oExact = exactTest(oEdge, pair=c("CLP_SI", "Media"))


### create some volcano plots with annotations
## choose the comparison for plotting
library(org.Mm.eg.db)
# add annotation to the data set after selecting comparison
res = dfResults

rn = rownames(res)
df = select(org.Mm.eg.db, as.character(rn), c('SYMBOL'), 'ENTREZID')
head(df); head(res)
df = df[!duplicated(df$ENTREZID),]
rownames(df) = df$ENTREZID
dfPlot = res
dfPlot = cbind(dfPlot[rn,], df[rn,])
dfPlot = na.omit(dfPlot)

## write csv file
write.csv(dfPlot, file='Results/DEAnalysis_CLP_SI.vs.SI_compM_edgeR.xls')

## what is the concordance between the 3 results 
dfEdgeR = dfPlot
dfDeseq = read.csv(file.choose(), header=T, row.names=1)
dfBayes = read.csv(file.choose(), header=T, row.names=1)

## select common genes 
cvGenes = unique(c(rownames(dfEdgeR), rownames(dfDeseq), rownames(dfBayes)))
dfGenes = data.frame(genes=cvGenes)
rownames(dfGenes) = cvGenes

## add the data e.g. p-values and fold changes
dfPv = data.frame(gene=dfGenes[cvGenes,], deseq=dfDeseq[cvGenes,'pvalue'], edger=dfEdgeR[cvGenes, 'PValue'], bayes=dfBayes[cvGenes, 'P.Value'])
dfPv = na.omit(dfPv)
dim(dfPv)
plot(dfPv$deseq, dfPv$edger, pch=20, xlab='deseq', ylab='edger', main='pvalue')
plot(dfPv$deseq, dfPv$bayes, pch=20, xlab='deseq', ylab='bayes', main='pvalue')
plot(dfPv$edger, dfPv$bayes, pch=20, xlab='edger', ylab='bayes', main='pvalue')

dfPv = data.frame(gene=dfGenes[cvGenes,], deseq=dfDeseq[cvGenes,'log2FoldChange'], edger=dfEdgeR[cvGenes, 'logFC'], bayes=dfBayes[cvGenes, 'logFC'])
dfPv = na.omit(dfPv)
dim(dfPv)
plot(dfPv$deseq, dfPv$edger, pch=20, xlab='deseq', ylab='edger', main='fc')
plot(dfPv$deseq, dfPv$bayes, pch=20, xlab='deseq', ylab='bayes', main='fc')
plot(dfPv$edger, dfPv$bayes, pch=20, xlab='edger', ylab='bayes', main='fc')
