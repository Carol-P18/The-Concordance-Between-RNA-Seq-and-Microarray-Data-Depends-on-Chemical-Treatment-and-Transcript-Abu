
#reads counts csv
ds1<- read.csv("countsfinal.csv",sep=',')
#ds1<- read.csv("sample1_n.csv",sep=',') #use for real data(3 samples)
#countsfinal csv generated from part 4
#

head (ds1)

#contol counts
ds2<- read.csv("control_counts.csv",sep=',')  # remains the same
ds1<- ds1[!colnames(ds1)%in%ds2_col]
#ds2<- read.csv("control_counts.csv",sep=',')
ds2_col <- colnames(ds2)
ds2_col <- ds2_col[ds2_col!='Geneid']
total<-merge(ds1,ds2,by="Geneid")
head(total)

library(DESeq2)
# load counts

cnts <- total
# filter out rows that have any zeros for funzies
cnts <- subset(cnts,rowSums(cnts==0)==0)

# READS TOX_GROUP FILES
#info <- read.csv("toxgroup_1_rna_info.csv")
info <- read.csv("toxgroup_EX_rna_info.csv")

###MEtoxgroup_EX_rna_info.csv CHANGE TO TOX GROUP 3 WHEN YOU GET CURATORS RESULSTS
#toxgroup_EX_rna_info.csv is precomputed data

info_col <- as.vector(info['Run'])
#Filter counts by sample names from meta data file
row.names(cnts) <- cnts[,1]
cnts <- cnts[,-1]
cnts_filter <- cnts
cnts_filter <- apply(cnts_filter,1:2,round)
rownames(info) <- info$Run
# Check for Gene IDS
all(row.names(info) %in% colnames(cnts_filter))
all(row.names(info) == colnames(cnts_filter))
info_row <- intersect(rownames(info),colnames(cnts_filter))
cnts_filter <- subset.matrix(cnts_filter,select = info_row)
# create the DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = cnts_filter,
  colData = info,
  design= ~ mode_of_action
)
# relevel mode_of_action as factor
dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')
# run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c('mode_of_action','HMGCOA','Control'))
res <- lfcShrink(dds, coef=2)
# write out DE results
write.csv(res,'deseq_results.csv')
# write out matrix of normalized counts
write.csv(counts(dds,normalized=TRUE),'deseq_norm_counts.csv')
#ds3<-read.csv("deseq_results.csv",sep=',')

library(dplyr)
res<-arrange(res,res$pvalue)
head(res,c=10)
sigDE <- subset(res,padj <0.05)
write.csv(sigDE,'sigDE.csv')
hist(sigDE$log2FoldChange,breaks=100)
plot(-log(sigDE$pvalue)~sigDE$log2FoldChange)



###Volcan PLOT original
library(EnhancedVolcano)
#SDE SIG GENES #RES=ALL GENES
EnhancedVolcano(res, 
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-5, 8), pCutoff =0.05,FCcutoff = 1.5, width = 7 * aspect_ratio)
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install("EnhancedVolcano")

#,pCutoff =0.05, FCcutoff = 1.5, pCutoff =0.05,FCcutoff = 1.5, width = 7 * aspect_ratio, width = 7 * aspect_ratio,
##res is for nonsignificant genes

####ADJUSTED
#EnhancedVolcano(sigDE,
               # lab = rownames(sigDE),
                #x = 'log2FoldChange',
                #y = 'pvalue',
                #xlim = c(-8, 8),
               # title = 'N061011 versus N61311',
               # pCutoff = 10e-16,
               # FCcutoff = 1.5,
               # pointSize = 4.0,
               # labSize = 3.0,
               # shape = 8,
                #colAlpha = 1)

