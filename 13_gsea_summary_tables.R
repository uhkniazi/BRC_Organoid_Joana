# File: 13_gsea_summary_tables.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 02/02/2017
# Desc: summary tables after merging gsea data

## set variables and source libraries
source('header.R')

lFiles = list.files('Results/', pattern='*mSigDb_all.xls', full.names = T, ignore.case = T)

# load the files
ldfData = lapply(lFiles, function(x) as.data.frame(read.csv(x, row.names=1)))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = c('CLP_SI_vs_Media') 
names(ldfData.up) = sn

sn = c('CLP_SI_vs_Media') 
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
cn = c('CLP_SI_vs_Media') 
colnames(mMerged.up) = paste(cn, 'up', sep='-')
colnames(mMerged.down) = paste(cn, 'down', sep='-')

mMerged.c2 = cbind(mMerged.up, mMerged.down)

# remove na sections
dim(mMerged.c2)
mMerged.c2 = na.omit(mMerged.c2)
dim(mMerged.c2)
head(mMerged.c2)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c2.bin = getBinaryMatrix(mMerged.c2)

## group this matrix into combinations
mMerged.c2.bin.grp = mMerged.c2.bin
set.seed(123)
dm = dist(mMerged.c2.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c2.bin.grp = cbind(mMerged.c2.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c2.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c2.bin)
dfMerged.c2 = data.frame(round(mMerged.c2, 3), sig.pvals, groups, DB='mSigDB-Mouse')
str(dfMerged.c2)
head(dfMerged.c2)
tail(dfMerged.c2)

#write.csv(dfMerged.c2, file='Results/gsea_msigdb_c2_merged.xls')
#write.csv(dfMerged.c2, file='Results/gsea_msigdb_all_merged.xls')
## keep only the ones with significant p-values
table(dfMerged.c2$groups)
dfMerged.c2.sub = dfMerged.c2[dfMerged.c2$groups != 2,]
# 
# table(dfMerged.c5$groups)
# dfMerged.c5.sub = dfMerged.c5[dfMerged.c5$groups != 2,]
# 
# dfMerged = rbind(dfMerged.c2.sub, dfMerged.c5.sub)
# dfMerged = droplevels.data.frame(dfMerged)
# 
write.csv(dfMerged.c2.sub, file='Results/gsea_msigdb_significant_all_merged.csv')

### heatmaps
dfMerged = dfMerged.c2
df = dfMerged
head(df)
mMat = as.matrix(df[,c(1:2)])
head(mMat)
mMat = -10*log10(mMat+1e-16)
g1 = df[,'groups']
g1 = factor(as.character(g1))
levels(g1)
g2 = df[,'DB']
g2 = factor(as.character(g2))
levels(g2)

ann = data.frame(DB=g2, Group=g1 )
range(mMat)
mMat[mMat < 15] = 0 
mMat[mMat > 100] = 100

library(NMF)
library(RColorBrewer)
aheatmap(mMat, annRow = ann, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

pdf('Results/gsea_msigdb_significant_merged.pdf')
aheatmap(mMat, annRow = ann, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))
dev.off(dev.cur())