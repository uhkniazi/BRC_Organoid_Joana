# File: 04_fastQC_trimmed.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the fastq files after trimming
# Date: 17/01/2017


## set variables and source libraries
source('header.R')
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CFastqQuality/master/CFastqQuality.R'
download(url, 'CFastqQuality.R')

# load the required packages
source('CFastqQuality.R')
# delete the file after source
unlink('CFastqQuality.R')

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'Sample')
# another way to get the query, preferred
g_did
dfSample = dbGetQuery(db, "select title, group1, group2 from Sample where idData=11;")
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
head(dfSample)
#### get the names of the fastq files for first sequencing run
setwd('Data_external/fastqTrimmed/')
csFiles = list.files('.', pattern = '*.gz')

# match the files to samples by title
i = sapply(dfSample$title, function(x) grep(x, csFiles))

# index of list with no matching files
i2 = sapply(i, length)
# skip this step if i2 has no 0 values
# dfSample = dfSample[which(i2 != 0),]

# match the files to their respective index in the table
lFilesIndex = lapply(dfSample$title, function(x) grep(x, csFiles))
names(lFilesIndex) = dfSample$title


## perform the analysis one sample at a time
## function to write the qa files
write.qa = function(fls, indir, title){
  wd = getwd()
  setwd(indir)
  ob = CFastqQuality(fls, title)
  setwd(wd)
  cat(paste('done', title, '\n'))
  return(ob)
}

ivFilesIndex = unlist(lFilesIndex)

lOb = lapply(ivFilesIndex, function(x){
  write.qa(csFiles[x], getwd(), csFiles[x])
})

setwd(gcswd)
n = make.names(paste('CFastqQuality organoids rna_seq trimmed with less restrictions rds'))
lOb$meta = dfSample
lOb$desc = paste('CFastqQuality object from organoids rna_seq trimmed run with less restrictions', date())
n2 = paste0('~/Data/MetaData/', n)
save(lOb, file=n2)

## note: comment out as this entry has been made in db
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/', comment='organoids project FASTQ file quality data after trimming with less restrictions on trimmomatic')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

### create the plots of interest
getwd()
lOb$desc = NULL
lOb$meta = NULL
pdf(file='Results/qa.fastq_trimmed_less_restrictions.pdf')

iReadCount = sapply(lOb, CFastqQuality.getReadCount)
iReadCount = iReadCount/1e+6

barplot(iReadCount, las=2, main='Trimmed Read Count', ylab = 'No. of Reads in Millions', cex.names =0.6, col=grey.colors(2))

mQuality = sapply(lOb, function(x){
  m = mGetReadQualityByCycle(x)
  m = colMeans(m, na.rm = T)
  return(m)
})

matplot(mQuality, type='l', main='Trimmed base quality', ylab = 'Mean Score', xlab='Position in Read')

lReadWidth = lapply(lOb, iGetReadWidth)
boxplot(lReadWidth, las=2, main='Trimmed Read Width', ylab = 'Read Width', col=grey.colors(2), outline=F, xaxt='n')
axis(1, at=1:length(lReadWidth), labels = names(lReadWidth), cex.axis=0.7, las=2)
dev.off(dev.cur())

## some samples may have quality drops which can be seen by clustering
hc = hclust(dist(t(mQuality)))
plot(hc, main='Clustering of Trimmed data based on Per base quality', cex=0.6, xlab='', sub='')
abline(h = 3, col=2)

c = cutree(hc, h=3)
table(c)
# draw for each of the clusters
sapply(unique(c), function(x){
m = mQuality[,c == x]
matplot(m, type='l', main='Trimmed Per base quality Clusters', ylab = 'Mean Score', xlab='Position in Read', col=1:ncol(m))
legend('bottomleft', legend = c(colnames(m)), fill = 1:ncol(m))
})
#n = colnames(m)



