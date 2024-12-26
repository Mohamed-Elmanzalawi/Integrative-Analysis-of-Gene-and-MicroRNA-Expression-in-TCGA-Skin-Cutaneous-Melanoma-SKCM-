#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#install_github("vqv/ggbiplot")
library(DESeq2)
library(ggplot2)
library(ggbiplot)
library(ggsci)
library(RColorBrewer)
library(NMF)
library(pheatmap)
library(devtools)
library(ggfortify)
library(gtools)
#load the mirna data as matrix to use it in deseq2
data = as.matrix(read.table("F:/Integerative project/mirna", quote="\"", comment.char="" , row.names=1))
#load the melanoma data (metadata)
melanoma <- read.delim("F:/Integerative project/melanoma", header=T ,row.names = 1 )
#as.matrix(meta)
table(melanoma$vital_status)
#show the dimension of the data as how many columns that present samples and how many rows that present genes in this data
dim(data)
dim(melanoma)
#check if is there any missing value of expression 
sum(is.na(data))
#sum(is.na(melanoma))
#remove NA from condition column (vital_status)
meta <- melanoma[!is.na(melanoma$vital_status), ] 

#patients (samples) name are dotted in data and dashed in meta data so make the same dashed in both
#data1 = gsub(".", "-", colnames(data))
#colnames(data) 
colnames(data) <- gsub (".", "-", colnames(data), fixed = TRUE) 
colnames(data)
#intersect data and meta data to make cols of data equal raws of meta data to apply deseq2
data3=data[,intersect(colnames(data),rownames(meta))]
dim(data3)
dim(meta)
#plot histogram to visualize the distribution of the data 
hist(log2(data3+1), col = "orange", main="Histogram") 
#histogram result shows that the data is right skewed 
#scale the data by log2 transformation for better visualization, the +1 at the end of command is to avoid the infinity at log the values equal to zero 
boxplot(log(data3[1:5,]+1))
#draw QQ plot for checking the normality of data
qqnorm(data3[1,])
qqline(data3[1,])
# the qq plot results show that the data is not normally distributed as not all the points are on the line of qq norm 
#in deseq2, it is critical that the columns of the data and the rows of the metadata (melanoma) that gives information about samples, to be in the same order as DESeq2 will not know which column of the count matrix belongs to which row of the column of the data, these must be provided to deseq2 in a consistent order so the following code will do that
#meta=meta[sample(1:nrow(data3)),]
#data3 <- data3[,rownames(meta)]
meta=meta[colnames(data3),]
dim(data3)
dim(meta)
#check the order to know if they are in a consistent order or not
all(colnames(data3) == rownames(meta))
all(colnames(data3) %in% rownames(meta))
#mat <- data3[,-1]
#rownames(mat) <- data3[,1]
identical(rownames(meta), colnames(data3))
######################################
#The deseq2 package require the count data values to be integers 
#save the gene names in a variable
genes=row.names(data3)
#convert the data values to integers
data3=apply(round(data3),2,as.integer)
#rename the rows of the data
row.names(data3)=genes
###### DO the differential EXP analysis using DeSeq2 
#creat a deseq dataset object
dds= DESeqDataSetFromMatrix( countData = data3 , colData = meta, design = ~ vital_status )
# Run pipeline for differential expression steps
dds.run = DESeq(dds)
#specify how many conditions do you want to compare according to the phenotypic table
#meta$X_OS_IND <- gsub ("0", "LIVING", meta$X_OS_IND) 
#meta$X_OS_IND <- gsub ("1", "DEAD", meta$X_OS_IND) 
cond1="DECEASED"
cond2="LIVING"
#as.factor(meta$vital_status)
#meta$vital_status = as.factor(meta$vital_status)
str(meta$vital_status)
#specifying the contrast (to make a res object based on two specific conditions)
res=results(dds.run, contrast = c("vital_status",cond1,cond2))

# remove nulls
res=as.data.frame(res[complete.cases(res), ])

#choose the statstical significant differentaily expressed transcripts of mirnas (DETs) based on the p adjusted value less than 0.05 and biological significance  based on the fold change more than 2
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)> 1.2,]
#export the Dets into your current folder for further analysthis
write.csv(as.matrix(deseq.deg),file="deseq.det.csv", quote=F,row.names=T)
#convert values from log-ratios to fold changes
#FC = logratio2foldchange(res$log2FoldChange)
#results = cbind(FC , res)
#DEGS = results[results$padj<0.05 & abs (results$FC) > 0.5,]
#export the Dets into your current folder for further analysthis
#write.csv(as.matrix(DEGS),file="dets.csv", quote=F,row.names=T)

#drow DETs volcano plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="DECEASED vs LIVING"))
with(subset(res, padj<.05 & (log2FoldChange)> 1.2 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & (log2FoldChange)< -1.2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=-1.5,y=4,c("upregulated","downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=19)

####drow heatmap####
#normalize the data
dds2 <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds2, normalized=TRUE))

#extract counts values of DETs(diffrentially expressed transcripts(micrornas)) only for each stage
exp.degs=as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg), ])
aheatmap(log2(exp.degs+1), annCol =meta$vital_status, col = rev(brewer.pal(9,"RdBu")))

##############################
#PCA 
#transpose rows with columns 
data4 <- t(data3)
dim(data4)
dim(data3)
PCA <- prcomp(data4)
plot(PCA$x[,1], PCA$x[,2])
plot(PCA)
plot(PCA ,type="l")
biplot(PCA , scale = 0)
# Extract PC scores
pcaData <- as.data.frame(PCA$x[, 1:2])
pcaData <- cbind(pcaData, meta$vital_status) 
colnames(pcaData) <- c("PC1", "PC2", "Species")
# plot with ggplot 
library(ggplot2)
ggplot(pcaData) +
  aes(PC1, PC2, color = Species, shape = Species) + # define plot area
  geom_point(size = 2) 
#adding data points
percentVar <- round(100 * summary(PCA)$importance[2, 1:2], 0)
ggplot(pcaData, aes(PC1, PC2, color = Species, shape = Species)) + # starting ggplot2
  geom_point(size = 2) + # add data points
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + # x label
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + # y label
  ggtitle("Principal component analysis (PCA)") + # title
  theme(aspect.ratio = 1) # width and height ratio

ggplot(pcaData, aes(PC1, PC2, col = Species, fill = Species)) +
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")
#####################################################################
#PCA (gender)
#transpose rows with columns 
data4 <- t(data3)
dim(data4)
dim(data3)
PCA <- prcomp(data4)
plot(PCA$x[,1], PCA$x[,2])
plot(PCA)
plot(PCA ,type="l")
biplot(PCA , scale = 0)
# Extract PC scores
pcaData <- as.data.frame(PCA$x[, 1:2])
pcaData <- cbind(pcaData, meta$gender) 
colnames(pcaData) <- c("PC1", "PC2", "Gender")
# plot with ggplot 
library(ggplot2)
ggplot(pcaData) +
  aes(PC1, PC2, color = Gender, shape = Gender) + # define plot area
  geom_point(size = 2) 
#adding data points
percentVar <- round(100 * summary(PCA)$importance[2, 1:2], 0)
ggplot(pcaData, aes(PC1, PC2, color = Gender, shape = Gender)) + # starting ggplot2
  geom_point(size = 2) + # add data points
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + # x label
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + # y label
  ggtitle("Principal component analysis (PCA)") + # title
  theme(aspect.ratio = 1) # width and height ratio

ggplot(pcaData, aes(PC1, PC2, col = Gender, fill = Gender)) +
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")
