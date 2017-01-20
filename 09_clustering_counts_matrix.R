# File: 09_clustering_counts_matrix.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: cluster the samples based on count table
# Date: 20/01/2017


## set variables and source libraries
source('header.R')

######################### functions
# Name: f_Plot3DPCA
# Args: mComp = n X 3 matrix with first 3 components as the 3 column vectors
#       color = colours for the points
#       ... additional arguments to the plot function
# Rets: none
# Desc: takes the first 3 components and plots the data in a 3d plot
f_Plot3DPCA = function(mComp, color, ...) {
  x = mComp[,1]
  y = mComp[,2]
  z = mComp[,3]
  if (!require(scatterplot3d)) stop('scatterplot3d library required')
  scatterplot3d(x, y, z, color, ...)
}

## connect to mysql database to get find path to appropriate file
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 11) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
dfSample
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2 as phenotype, Sample.title, File.* from Sample, File
           where (Sample.idData = 11) AND (File.idSample = Sample.id AND File.type like "%duplicates removed%")')
dfSample.names = dbGetQuery(db, q)
dim(dfSample.names)
dfSample.names
# close connection after getting data
dbDisconnect(db)

n = paste0(dfSample$location, dfSample$name)

load(n)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)

# reorder the count matrix columns according to the order in samples table
i = match(dfSample.names$name, colnames(mCounts))
mCounts = mCounts[,i]
# sanity check
identical(dfSample.names$name, colnames(mCounts))

#### scale each variable vector i.e. gene
## add a normal jitter to each cell to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
mCounts = mCounts + rnorm(n)

#### standardize samples first
s = apply(mCounts, 2, sd)
mCounts.s = sweep(mCounts, 2, s, '/')

## PCA with strandardizing samples
mCounts.s = t(mCounts.s)
# set scaling to TRUE to scale variables i.e. genes in columns
pr.out = prcomp(mCounts.s, scale = T)

# set the factor for colours
fSamples = factor(dfSample.names$phenotype)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('topleft', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

## normalize these datasets with DESeq2 
library(DESeq2)

n = paste0(dfSample$location, dfSample$name)

load(n)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)

# reorder the count matrix columns according to the order in samples table
i = match(dfSample.names$name, colnames(mCounts))
mCounts = mCounts[,i]
# sanity check
identical(dfSample.names$name, colnames(mCounts))

dfDesign = data.frame(condition=factor(dfSample.names$phenotype), row.names = colnames(mCounts))

oDseq = DESeqDataSetFromMatrix(mCounts, dfDesign, design = ~ condition)
mCounts.rlog = assays(rlog(oDseq))
oDseq = DESeq(oDseq)
mCounts.norm = counts(oDseq, normalized=T)
mCounts.rlog = mCounts.rlog@listData[[1]]

# adding some noise to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
mCounts.s = t(mCounts.rlog + rnorm(n))

# set scaling to TRUE to scale columns 
pr.out = prcomp(mCounts.s, scale = T)

# set the factor for colours
fSamples = factor(dfSample.names$phenotype)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, log transformed')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('bottomleft', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

## for normalized data
# adding some noise to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
# adding noise to avoid negative numbers
mCounts.s = t(mCounts.norm + rnorm(n))

# set scaling to TRUE to scale columns 
pr.out = prcomp(mCounts.s, scale = T)

# set the factor for colours
fSamples = factor(dfSample.names$phenotype)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]

plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalised')
text(pr.out$x[,1:2], labels = dfSample.names$title, pos = 1, cex=0.6)
legend('bottomleft', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], cex=0.8)

############## some additional QC checks with distribution of gene expressions
## and sample types
n = paste0(dfSample$location, dfSample$name)
load(n)
## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)

# reorder the count matrix columns according to the order in samples table
i = match(dfSample.names$name, colnames(mCounts))
mCounts = mCounts[,i]
# sanity check
identical(dfSample.names$name, colnames(mCounts))

# adding some noise to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
mCounts.s = mCounts + rnorm(n)

#### standardize the genes first 
s = apply(mCounts.s, 1, sd)
mCounts.s = sweep(mCounts.s, 1, s, '/')

# get the mean vector and total vector for each sample
ivMean = colMeans(mCounts.s)
ivTotal = colSums(mCounts.s)

dfData = data.frame(ivMean, ivTotal, condition = dfSample.names$phenotype)
library(lattice)

densityplot(~ ivMean, data=dfData, groups=condition, auto.key=TRUE, main='Average Gene Expression Density in Each Phenotype',
            xlab='Mean Gene Expression')

## for normalized data
# adding some noise to remove zeros
n = dim(mCounts)[1] * dim(mCounts)[2]
# adding noise to avoid negative numbers
mCounts.s = mCounts.norm + rnorm(n)

# get the mean vector and total vector for each sample
ivMean = colMeans(mCounts.s)
ivTotal = colSums(mCounts.s)

dfData = data.frame(ivMean, ivTotal, condition = dfSample.names$phenotype)

densityplot(~ ivMean , data=dfData, groups=condition, auto.key=TRUE, main='Average Gene Expression Density in Each Phenotype, Normalised',
            xlab='Mean Gene Expression')

dotplot(ivMean ~ condition, data=dfData, auto.key=TRUE, main='Average Gene Expression in Each Sample, Normalised',
        xlab='Mean Gene Expression', pch=20, cex.axis=0.7)





