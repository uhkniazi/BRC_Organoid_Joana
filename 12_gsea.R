# File: 12_gsea.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 02/02/2017
# Desc: gene set enrichment analysis for the datasets

## set variables and source libraries
source('header.R')
## libraries to load
library(gage)

## load gene expression data
dfContrast = read.csv(file='Results/DEAnalysis_CLP_SI.vs.Media_3.xls', stringsAsFactors = F)
dfContrast = dfContrast[,-1]

## load msig db data
oMsigGS = readList(file.choose())
## OR try other datasets
# library(gageData)
# data("kegg.sets.hs")
# oMsigGS = kegg.sets.hs
## GO
# data(go.sets.hs)
# data(go.subs.hs)
# names(go.subs.hs)
# oMsigGS=go.sets.hs[go.subs.hs$MF]
# oMsigGS=go.sets.hs[go.subs.hs$BP]
# oMsigGS=go.sets.hs[go.subs.hs$CC]

#dfContrast = dfContrast4
# for a contrats of choice create the list
iContFc = dfContrast$logFC
names(iContFc) = as.character(dfContrast[,1])
head(iContFc)
head(dfContrast)

oGage = gage(iContFc, oMsigGS)

dfGreater = data.frame(oGage$greater)
str(dfGreater)
i = which(dfGreater$p.val < 0.01)
rownames(dfGreater[i,])

dfLess = data.frame(oGage$less)
str(dfLess)
i = which(dfLess$p.val < 0.01)
rownames(dfLess[i,])

write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file='Results/upregulated_pathways_msigdb_all.xls')
write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file='Results/downregulated_pathways_msigdb_all.xls')



## make a barplot
plot.bar = function(ivBar, title='', ...){
  p.old = par(mar=c(8,4,2,2)+0.1)
  l = barplot(ivBar, beside=T, xaxt='n', main=title, ...)
  axis(side = 1, l[,1], labels=F)
  text(l[,1], y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=names(ivBar), srt=45, adj=1, xpd=TRUE, cex=0.6)
  par(p.old)
}

dfBar = dfGreater[order(dfGreater$p.val), ]
dfBar = dfLess[order(dfLess$p.val), ]

i = -1*log10(dfBar$p.val[1:7])
names(i) = rownames(dfBar)[1:7]

plot.bar(i, title='Contrast 4 - down regulated GO BP', ylab='-log10 PValue')

# #### test with kegg pathways data
# library(pathview)
# dir.create('Results/Kegg_figures')
# pv.out.list <- sapply(substr(rownames(dfBar)[1], 1, 8), function(pid) pathview(
#   gene.data = iContFc, pathway.id = pid,
#   species = "hsa", out.suffix='Contrast 1 Down', kegg.dir = 'Results/Kegg_figures/'))
# 
# 
# ### protein enrichment analysis using interaction database 
# # ## some enrichment using string db package
# library(STRINGdb)
# # find database id for human
# temp = get_STRING_species(version = '10')
# temp[grep('homo', temp$official_name, ignore.case = T, perl=TRUE),]
# 
# dbString = STRINGdb$new(version='10', species=9606, input_directory='')
# 
# df.rn = dfContrast4
# str(df.rn)
# df.rn$ENTREZ.ID = as.character(df.rn$ENTREZ.ID)
# rownames(df.rn) = df.rn$ENTREZ.ID
# head(df.rn)
# 
# string.mapped = dbString$map(df.rn, 'ENTREZ.ID', removeUnmappedRows = T)
# 
# string.mapped = string.mapped[order(string.mapped$p.value),]
# head(string.mapped)
# 
# # plot enrichment
# dbString$plot_ppi_enrichment(string.mapped$STRING_id, quiet = T, title = 'Contrast 4 All')
# dbString$plot_ppi_enrichment(string.mapped$STRING_id[1:700], quiet = T, title = 'Contrast 4')
# #ig = dbString$get_graph()
# #enrichment = ppi_enrichment_full(string.mapped$STRING_id, ig)
# #enrichment = ppi_enrichment(string.mapped$STRING_id, ig)
# #enrichment = dbString$get_ppi_enrichment(string.mapped$STRING_id[1:1000])
# par(p.old)
# hits = string.mapped$STRING_id[250:430]
# #dbString$plot_network(hits)
# pdf('Temp/string_contrast4_250-430.pdf')
# dbString$plot_network(hits)
# dev.off(dev.cur())
# 
# 
# # category for which to compute the enrichment (i.e. "Process", "Component",
# # "Function", "KEGG", "Pfam", "InterPro"). The default category is "Process".
# enrichmentProcess = dbString$get_enrichment(string.mapped$STRING_id[1:500], category = 'Process', iea=F)
# enrichmentKEGG = dbString$get_enrichment(string.mapped$STRING_id[1:500], category = 'KEGG', iea=F)
# 
# ## make bar plots and save data
# write.csv(enrichmentProcess, file='Results/Keloid:2vsControl:2_protein_pathways_GO_BP.xls')
# write.csv(enrichmentKEGG, file='Results/Keloid:2vsControl:2_protein_pathways_KEGG.xls')
# 
# dfBar = enrichmentProcess[order(enrichmentProcess$pvalue), ]
# dfBar = enrichmentKEGG[order(enrichmentKEGG$pvalue), ]
# 
# i = -1*log10(dfBar$pvalue[1:7])
# names(i) = dfBar$term_description[1:7]
# 
# plot.bar(i, title='Contrast 4 - Proteins KEGG', ylab='-log10 PValue')
# 
# 
# # #################### compare 2 contrasts
# # df.rn = dfContrast4
# # str(df.rn)
# # df.rn$ENTREZ.ID = as.character(df.rn$ENTREZ.ID)
# # rownames(df.rn) = df.rn$ENTREZ.ID
# # head(df.rn)
# # 
# # string.mapped = dbString$map(df.rn, 'ENTREZ.ID', removeUnmappedRows = T)
# # 
# # string.mapped = string.mapped[order(string.mapped$p.value),]
# # 
# # string.mapped.1 = string.mapped
# # string.mapped.2 = string.mapped
# # 
# # hits.contrast1 = string.mapped.1$STRING_id[1:600]
# # hits.contrast3 = string.mapped.2$STRING_id[1:600]
# # 
# # eh = dbString$enrichment_heatmap(list(hits.contrast1, hits.contrast3), vectorNames = list('Contrast1', 'Contrast3'), output_file = 'Temp/hm.pdf',
# #                                  enrichmentType = 'Process',
# #                                 title = 'Contrats 1 and Contrast 2 GO BP')
# 
