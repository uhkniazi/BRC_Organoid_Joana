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
  betaSigma = exp(theta['betaSigma'])+1
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
  lData$betaParam = unlist(gammaShRaFromModeSD(mode = sd(dfData$resp)/2, sd = 2*sd(dfData$resp)))
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
write.csv(dfPlot, file='Results/DEAnalysis_CLP_SI.vs.Media.xls')

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
## all overexpressed genes if interested in
fSamples = oExp.norm$fCondition
i = 1:nrow(mCommonGenes)

m1 = as.matrix(mCommonGenes[i,])
mDat = counts(oDseq, normalized=T)
dim(mDat)
m1 = mDat[rownames(m1),]
dim(m1)

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
dfGenes = select(org.Mm.eg.db, cvCommonGenes, c('SYMBOL', 'GENENAME'), 'ENTREZID')
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
dfGenes = select(org.Mm.eg.db , keys = cvGenes, columns = c('ENTREZID', 'SYMBOL', 'GENENAME', 'UNIPROT'), keytype = 'ENTREZID')
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
dim(dfGraph)
head(dfGraph)
# get expression data
mCounts = mDat[unique(dfGenes$ENTREZID),]
# ## optional adjusting the data for repeated measurements
# temp = sapply(rownames(mCounts), function(x){
#   return(fitted(lGlm.sub[[x]]))
# })
# mCounts = t(mCounts)

fGroups = fSamples
names(fGroups) = oExp$title
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