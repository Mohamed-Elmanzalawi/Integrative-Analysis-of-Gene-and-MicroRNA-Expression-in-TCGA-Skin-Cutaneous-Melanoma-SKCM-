#BiocManager::install("DESeq2", force = TRUE)#one time
library(DESeq2)#each time you open your R
library(RColorBrewer)
library(NMF)
library(ggplot2)


########################### Setting Up our data ##################################
#loading the gene expression data and as matrix to allow mathematical calculations
data <- as.matrix(read.csv("F:/Files/Diploma/Integrative/Project/Data/exp",row.names=1, header = T, sep=""))
#loading the meta data 

melanoma <- read.delim("F:/Files/Diploma/Integrative/Project/Data/melanoma")
#Knowing how many deceased and living.
table(melanoma$vital_status)
#explore the data, dim function give you the dimension of your data; how 
#many columns(samples) do you have and how many row(genes) do you have
dim(data)
dim(melanoma)
#explore if is there any missing expression value (empty cell)
sum(is.na(data))
sum(is.na(melanoma$vital_status))
#Changing names of patients in melanoma data to meta data
#The difference was that there was a dot in patient names in melanoma  data
#while meta data had - .
for (i in 1:nrow(melanoma)){
  melanoma[i,1]=gsub("-",".",melanoma[i,1])
}
#setting first column as rownames
rownames(melanoma)=melanoma[,1]
#removing NA values
melanoma <- melanoma[complete.cases(melanoma[ , c("vital_status")]), ]
#Getting the intersection between both melanoma and expression data.
melanoma3 = melanoma[intersect(colnames(data),row.names(melanoma)),]
data2=data[,intersect(colnames(data),row.names(melanoma))]

#####explore the data distribution##############################################
#Using histogram
hist(data, col = "orange", main="Histogram", breaks = 100)
#scaling the data using log2 transformation to better visulization
# we use (+1) to avoid the infinity character when we log zero valus 
hist(log2(data+1), col = "orange", main="Histogram")
# Drawing a box plot
boxplot(log2(data[1:5,]+1))
# QQ plot for the normality
qqnorm(data[30,])
qqline(data[30,])

################################################################################
#It is absolutely critical that the columns of the "data" and the rows of 
#the "melanoma" (information about samples) are in the same order.DESeq2 will
#not make guesses as to which column of the count matrix belongsto which row
#of the column data, these must be provided to DESeq2 already in consistent order
melanoma3=melanoma3[colnames(data2),]
dim(melanoma3)
dim(data2)

#The deseq2 package require the count data values to be integers 
#save the gene names in a variable
genes=row.names(data2)
#convert the data values to integers
data2=apply(round(data2),2,as.integer)
#view the data
head(data2)
#rename the rows of the data
row.names(data2)=genes

###### DO the differential EXP analysis using DeSeq2#############################

#specify how many conditions do you want to compare according to 
cond1="DECEASED" 
cond2="LIVING"

#creat a deseq dataset object
dds= DESeqDataSetFromMatrix( countData = data2 , colData = melanoma3, design = ~vital_status)
#run the deseq2 worflow
dds.run = DESeq(dds)
#specifying the contrast (to make a res object based on two specific conditions)
res=results(dds.run, contrast = c("vital_status",cond1 ,cond2))
# remove nulls
res=as.data.frame(res[complete.cases(res), ])

#chose the statstical significant differentaily expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 1.2
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>1.5,]
#export the Degs into your current folder for further analysthis
write.csv(as.matrix(deseq.deg),file="F:/Files/Diploma/Integrative/Project/DEGs.csv", quote=F,row.names=T)

###### Exporting Degs in the form required for GSEA ############################
Deceased=c()
Living=c()
for (i in 1:nrow(melanoma3)) {
  if ( melanoma3[i,c("vital_status")] == "DECEASED" ){
    Deceased=append(Deceased,melanoma3[i,1])
  }else{
    Living=append(Living,melanoma3[i,1])
  }
}
Deceased_df=data.frame(Deceased)
Living_df=data.frame(Living)

Deceased_df=data[,intersect(Deceased_df$Deceased,colnames(data))]
Living_df=data[,intersect(Living_df$Living,colnames(data))]

row.names(Deceased_df)=gsub(".*\\.", "", row.names(Deceased_df))
row.names(Living_df)=gsub(".*\\.", "", row.names(Living_df))

write.csv(as.matrix(Deceased_df),file="F:/Files/Diploma/Integrative/Project/DEGs_GSEA_D.csv", quote=F,row.names=T)
write.csv(as.matrix(Living_df),file="F:/Files/Diploma/Integrative/Project/DEGs_GSEA_L.csv", quote=F,row.names=T)

###### Exporting Degs in the form required for Gprofiler #######################
row.names(deseq.deg)=gsub("\\..*", "", row.names(deseq.deg))
write.csv(as.matrix(deseq.deg),file="F:/Files/Diploma/Integrative/Project/DEGs_Gprofiler.csv", quote=F,row.names=T)

###### Drawing the volcano plot ################################################

par(mfrow=c(1,1))

with(res, plot(log2FoldChange, -log10(padj), pch=20, 
                   main = "DECEASED vs LIVING DEGs", 
                   col="gray", xlim = c(-3,3), ylim = c(0,15)))

with(subset(res, padj < 0.05 & log2FoldChange > 1.5), 
     points(log2FoldChange, -log10(padj), pch=20, 
            col="red", xlim = c(-3,3), ylim = c(0,15)))

with(subset(res, padj < 0.05 & log2FoldChange < 1.5), 
     points(log2FoldChange, -log10(padj), pch=20, 
            col="blue", xlim = c(-3,3), ylim = c(0,15)))
legend("topright", legend = c('Upregulated', 'Downegulated'), 
       col = c('red', 'blue'), pch = c(20, 20))
lines(c(0.6,0.6), c(0,15), lwd = 2, lty = 2)
lines(c(-0.6,-0.6), c(0,15), lwd = 2, lty = 2)
lines(c(-3.5,3.5), c(1.3,1.3), lwd = 2, lty = 2)


###### Drawing the Heat map ###################################################
#normalize the data
dds2 <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds2, normalized=TRUE))
#extract counts values of DEGs only for each stage
exp.degs=as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg), ])
aheatmap(log2(exp.degs+1), annCol =melanoma3$vital_status, col = rev(brewer.pal(9,"RdBu")), main="mRNA DECEASED vs LIVING")
heatmap(log2(exp.degs+1))


######### PCA Based on vital status #############################################

#transpose rows with columns 
data4 <- t(data2)
dim(data2)
dim(data4)
PCA <- prcomp(data4)
plot(PCA$x[,1], PCA$x[,2])
plot(PCA)
plot(PCA ,type="l")
biplot(PCA , scale = 0)
# Extract PC scores
pcaData <- as.data.frame(PCA$x[, 1:2])
pcaData <- cbind(pcaData, melanoma3$vital_status) 
colnames(pcaData) <- c("PC1", "PC2", "Status")
# plot with ggplot 
ggplot(pcaData) +
  aes(PC1, PC2, color = Status, shape = Status) + # define plot area
  geom_point(size = 2) 
#adding data points
percentVar <- round(100 * summary(PCA)$importance[2, 1:2], 0)
ggplot(pcaData, aes(PC1, PC2, color = Status, shape = Status)) + # starting ggplot2
  geom_point(size = 2) + # add data points
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + # x label
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + # y label
  ggtitle("Principal component analysis (PCA)") + # title
  theme(aspect.ratio = 1) # width and height ratio

ggplot(pcaData, aes(PC1, PC2, col = Status, fill = Status)) +
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")

######### PCA Based on Gender ##################################################
# Extract PC scores
pcaData2 <- as.data.frame(PCA$x[, 1:2])
pcaData2 <- cbind(pcaData2, melanoma3$gender)
colnames(pcaData2) <- c("PC1", "PC2", "Gender")
# plot with ggplot 
ggplot(pcaData2) +
  aes(PC1, PC2, color = Gender, shape = Gender) + # define plot area
  geom_point(size = 2) 
#adding data points
percentVar <- round(100 * summary(PCA)$importance[2, 1:2], 0)
ggplot(pcaData2, aes(PC1, PC2, color = Gender, shape = Gender)) + # starting ggplot2
  geom_point(size = 2) + # add data points
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + # x label
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + # y label
  ggtitle("Principal component analysis (PCA)") + # title
  theme(aspect.ratio = 1) # width and height ratio

ggplot(pcaData2, aes(PC1, PC2, col = Gender, fill = Gender)) +
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")


