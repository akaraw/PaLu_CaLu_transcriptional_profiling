####################################### ~ Loading libraries ~ #####
library("DESeq2")
library(dplyr)
library("RColorBrewer")
library("gplots")
library( "genefilter" )
library("pheatmap")
library(tools)
library("ggplot2")
library(magrittr)
library("biomaRt")
library(apeglm)
library("genefilter")
library(tximport)
library(pcaExplorer)
library(GenomicFeatures)
library(goseq)
library(GO.db)
library(fgsea)
library('EnhancedVolcano')
library(PCAtools)
library(dplyr)

####################################### ~ DECLARE YOUR VARIABLES HERE ~ #####
workdir <- "C:/Users/kar131/OneDrive - CSIRO/"
projdir <- "M:/work/BatProject_MB/DEseq2/"
htseqdir <- paste0(projdir, "htseq/")
filemeta <- "M:/work/BatProject_MB/DEseq2/ID_and_labels.txt"
resultsdir <- "Human_bat_results"
(resultsdir <- paste0(workdir, resultsdir))
dir.create(resultsdir)

setwd(resultsdir)
sample_table <- read.table(paste0(filemeta), sep = "\t", header=F)
length(sample_table)
head(sample_table)
(sample_table)

rownames(sample_table) <- paste0(sample_table$V1, "ortho.tsv")
rownames(sample_table)

sampleFiles <- paste0(htseqdir, rownames(sample_table))

input_path <- htseqdir
file_list <- list.files(input_path, pattern = "*ortho.tsv", full.names = T)
file_list

for (file in file_list) {
  df <- read.csv(file, header = F, stringsAsFactors = FALSE, sep = "\t")
  t <- tail(df, 4)
  print(t)
}

all(file.exists(sampleFiles))
sampleFiles
sampleTable <- data.frame(sampleName = sample_table$V1, 
                          fileName = rownames(sample_table),
                          sp = as.factor(sample_table$V3),
                          time = as.factor(sample_table$V4),
                          infection = as.factor(sample_table$V5))

filedir <- htseqdir
setwd(htseqdir)
sampleTable

# create a data frame with the sample names and file paths
samples <- data.frame(
  sp = sampleTable$sp,
  time = sampleTable$time,
  samplename = sampleTable$fileName,
  infection = sampleTable$infection)

countList <- lapply(samples$samplename, read.table, header=F, row.names=1)
countList

countMatrix <- do.call(cbind, countList)
countMatrix
class(countMatrix)
tail(countMatrix)
nrow(countMatrix)

#remove last fine lines
countMatrix <- countMatrix[1:(nrow(countMatrix)-5),]
tail(countMatrix)

dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = samples,
                              design = ~ sp + time + infection)
nrow(counts(dds))
rownames(rowData(dds))
colnames(colData(dds))
dds$sp
dds$time
dds$infection

baselevel <- "mock"
dds$infection<- relevel(dds$infection, ref = baselevel)

baselevel <- "Pteropus alecto"
dds$sp <- relevel(dds$sp, ref = baselevel)

keep <- rowSums(counts(dds)) > 10
ddsKEEP <- dds[keep,]

ddsKEEP <- DESeq(ddsKEEP) #DESeq analysis
colData(ddsKEEP)
resultsNames(ddsKEEP)

nonnormalized_counts <- counts(dds)
resultsdir
write.table(nonnormalized_counts, file=paste(resultsdir, "/sarscov2_HandB_non_normalized_counts.txt", sep = ""), sep="\t", quote=F, col.names=NA)

normalized_counts <- counts(ddsKEEP, normalized=TRUE)
head(normalized_counts)
t <- as.data.frame(normalized_counts)
head(t)
colnames(t)
resultsdir
class(normalized_counts)
write.table(t, file=paste(resultsdir, "/handb_normalized_counts.txt", sep = ""), sep="\t", quote=F, col.names=NA)

vsdata <- vst(ddsKEEP, blind=FALSE)
plotPCA(vsdata, intgroup=c("sp", "infection", "time"))
resname <- 'sp_Homo.sapiens_vs_Pteropus.alecto'
res <- results(ddsKEEP, name = resname)
res
colData(ddsKEEP)
table(is.na(res$padj))
sum(res$padj<0.01, na.rm = T)
res[which(res$padj<0.01),]
sum(res$log2FoldChange>2, na.rm = T)
res[which(res$log2FoldChange>2),]
resSig <- res[which(res$padj < 0.01),]
resSig

write.csv( as.data.frame(resSig), file=paste(resultsdir, "/", resname, "_DEGS.csv", sep = "") )

assayed.genes <- rownames(ddsKEEP)
topGene <- rownames(res)[which.min(res$padj)]
topGene
resname

plotCounts(ddsKEEP, gene=topGene, intgroup=c("sp"))
mygene <- topGene
d <- plotCounts(ddsKEEP,gene = mygene, intgroup = "sp", main = paste( "gene -", mygene), returnData=TRUE)
p <- ggplot(d, aes(x=sp, y=count)) + 
  #geom_violin()  +
  geom_boxplot(colour = "red", fill = "orange", alpha = 0.2, 
               outlier.colour="black", outlier.shape=8, outlier.size=2, notch=F, notchwidth = 0.5) + 
  #geom_dotplot(binwidth = 50, binaxis='y', stackdir='center', dotsize=1)
  geom_point(position=position_jitter(w=0.1,h=0), colour = 'purple', size = 1) + 
  scale_y_log10(breaks=c(25,100,400)) + 
  theme(
    #panel background elements
    panel.background = element_rect(
      fill = "grey90",
      colour = "black",
      size = 1,
    ),
    legend.position= "bottom",
    plot.title = element_text(color="black", size=12, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=11, face="bold"),
    axis.title.y = element_text(color="#993333", size=11, face="bold")
  ) + 
  ggtitle(paste0("Condition: ", " gene - ", mygene)) + xlab("Species") + 
  ylab("Noramlized gene count") +
  labs(fill = d$sp) +
  stat_summary(fun=mean, geom="point", shape=23, size=4) + scale_color_grey() +
  scale_fill_manual(values=c("#999999", "#E69F00"))
print(p)

(resSig_up= resSig[which(resSig$log2FoldChange >= 2),])
resSig_up=resSig_up[order(resSig_up$log2FoldChange),]
head(resSig_up)
write.csv( as.data.frame(resSig_up), file=paste(resultsdir, "/", resname, "_UP_DEGS.csv", sep = "" ))


dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = samples,
                              design = ~ sp + time + infection + sp:infection)

baselevel <- "mock"
dds$infection

dds$infection<- relevel(dds$infection, ref = baselevel)

dds$sp
baselevel <- "Pteropus alecto"
dds$sp <- relevel(dds$sp, ref = baselevel)

keep <- rowSums(counts(dds)) > 10
ddsKEEP <- dds[keep,]

ddsKEEP <- DESeq(ddsKEEP) 
coef(ddsKEEP)[resultsNames(ddsKEEP) == "infection_Infe_vs_mock"]
#DESeq analysis
#colData(ddsKEEP)
resultsNames(ddsKEEP)

#resname <- 'spPteropus.alecto.infectionmock'
resname <- 'sp_Homo.sapiens_vs_Pteropus.alecto'
res <- results(ddsKEEP, name = resname)

colData(ddsKEEP)
table(is.na(res$padj))
sum(res$padj<0.01, na.rm = T)
res[which(res$padj<0.01),]
sum(res$log2FoldChange>2, na.rm = T)
res[which(res$log2FoldChange>2),]
resSig <- res[which(res$padj < 0.01),]
resSig

write.csv( as.data.frame(resSig), file=paste(resultsdir, "/", resname, "_DEGS.csv", sep = "") )

resultsNames(ddsKEEP)
resname <- 'sp_Homo.sapiens_vs_Pteropus.alecto'
res <- results(ddsKEEP, name = resname)

colData(ddsKEEP)
table(is.na(res$padj))
sum(res$padj<0.01, na.rm = T)
res[which(res$padj<0.01),]
sum(res$log2FoldChange>2, na.rm = T)
res[which(res$log2FoldChange>2),]
resSig <- res[which(res$padj < 0.01),]
resSig

write.csv( as.data.frame(resSig), file=paste(resultsdir, "/", resname, "_DEGS.csv", sep = "") )
