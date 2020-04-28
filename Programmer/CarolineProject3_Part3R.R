library("edgeR")
library("magrittr")
library("stringr")

#args <- commandArgs(trailingOnly = TRUE)
#stopifnot(length(args) > 0, file.exists(args))
#f_counts <- args
# For testing:

f_counts <- Sys.glob("SRR*txt") # FROM FEATURE COUNTS

raw <- readDGE(f_counts, columns = c(1,7), comment.char = "#",header=TRUE)
counts <- as.matrix(raw)

# Change column names from filenames to sample names ------------------------

#COMBINES THE 9 COUNT FILES
#names(dimnames(counts)) <- c("Gene_id", "")
counts=data.frame(counts)

##combines the count files and R SAMPLE COUNT FILES
counts=cbind(rownames(counts),counts)
colnames(counts)[1]="Geneid"
rownames(counts)=NULL

counts <- counts[order(rownames(counts)), order(as.numeric(colnames(counts)))]


# Write to standard out --------------------------------------------------------
# head(counts)
write.csv(counts, file = "countsfinal.csv", quote = TRUE,row.names=FALSE)
par(mar=c(8,2,1,1))
counts.nozero=counts[which (rowSums(counts[,-1])>0),]
boxplot(log2(counts.nozero[,-1]),las=2,col="green", frame.plot = FALSE)
#dev.copy(png,'box plots results for each sample.png')
#dev.off()
