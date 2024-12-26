
library(MOFA2)
library(MOFAdata)

#load data
#methyl
methyl = read.csv("F:/Integerative project/methy", row.names=1, sep="")
#convert B values to M values
methyl[] <- apply(methyl, 2, function(x) log2(x / (1 - x)))
hist(as.matrix((methyl)))


#exp
exp <- read.csv("F:/Integerative project/exp", row.names = 1, sep="")
#exploration of RNA Expression distribution
hist(as.matrix(log2(exp+1)))
dim(exp)
#mirna
mirna <- read.csv("F:/Integerative project/mirna", row.names=1, sep="")
#exploration of miRNA Expression distribution
hist(as.matrix(log2(mirna+1)))
dim(mirna)
#intersect
new_methyl=methyl[,intersect(intersect(colnames(methyl),colnames(exp)),colnames(mirna))]
dim(new_methyl)
#creating a new matrix contains the intersecting samples between exp matrix and mirna matrix
new_exp=exp[,intersect(colnames(new_methyl),colnames(exp))]
dim(new_exp)
#creating a new matrix contains the intersecting samples between mirna matrix and exp matrix
new_mirna=mirna[,intersect(colnames(mirna),colnames(new_methyl))]
dim(new_mirna)

###################


#loading meta_data
melanoma <- read.delim("F:/Integerative project/melanoma")
row.names(melanoma)=melanoma$sampleID
row.names(melanoma) <- gsub ("-", ".", row.names(melanoma), fixed = TRUE)
melanoma$sampleID <- gsub ("-", ".", melanoma$sampleID, fixed = TRUE)
#creating new meta-date matrix that contain the intersecting samples with new_exp and new_mirna
new_melanoma=melanoma[intersect (row.names(melanoma),colnames(new_mirna)), ]
colnames(new_melanoma)[1]="sample"
#creating the meta-date of MOFAobject
Project_data <- make_example_data(
  n_views = 3, 
  n_samples = 200, 
  n_features = 1000, 
  n_factors = 10
)[[1]]
lapply(Project_data,dim)
names(Project_data)= c("Transcriptomic", "miRNA", "methyl")
Project_data [["Transcriptomic"]]<- as.matrix(new_exp, rownames.force = NA)
Project_data [["miRNA"]]<- as.matrix(new_mirna, rownames.force = NA)
Project_data [["methyl"]]<- as.matrix(new_methyl, rownames.force = NA)

#creat MOFA opject
MOFAobject <- create_mofa(Project_data, extract_metadata = F)
MOFAobject
samples_metadata(MOFAobject)=new_melanoma

#overview of training data ## The rows are the different views and columns are samples. Missing values are indicated by a grey bar.
plot_data_overview(MOFAobject)

#Fit the MOFA model
#Define options
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views = TRUE


#Define model options
model_opts <- get_default_model_options(MOFAobject)
model_opts

#Define training options
train_opts <- get_default_training_options(MOFAobject)
train_opts
#TrainOptions$tolerance=0.01

#prepareMOFA internally performs a set of sanity checks and fills the DataOptions, TrainOptions and ModelOptions slots of the MOFAobject
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#Run MOFA This step can take some time (around 15 min with default parameters
MOFAobject <- run_mofa(MOFAobject,use_basilisk = T)

#Analyse a trained MOFA model 
#Disentangling the heterogeneity: calculation of variance explained by each factor in each view
# Calculate the variance explained (R2) per factor in each view 
r2 <- get_variance_explained(MOFAobject)
r2$r2_total

# Variance explained by each factor in each  
head(r2$r2_per_factor)

# Plot it
plot_variance_explained(MOFAobject)

##############################################

#Characterisation of individual factors
#To get an overview of the weights across all factors in a given view
plot_weights_heatmap(
  MOFAobject, 
  view = "Transcriptomic", 
  factors = 1:5,
  show_colnames = FALSE
)
plot_weights_heatmap(
  MOFAobject, 
  view = "methyl", 
  factors = 1:5,
  show_colnames = FALSE
)

plot_weights_heatmap(
  MOFAobject, 
  view = "miRNA", 
  factors = 1:5,
  show_colnames = FALSE
)
#To explore a given factor in more detail we can plot all weights for a single factor using the plotWeights function.
plot_weights(
  MOFAobject, 
  view = "miRNA", 
  factor = 1, 
  nfeatures = 5
)
#To explore a given factor in more detail we can plot all weights for a single factor using the plotWeights function.
plot_weights(
  MOFAobject, 
  view = "Transcriptomic", 
  factor = 2, 
  nfeatures = 5
)
plot_weights(
  MOFAobject, 
  view = "Transcriptomic", 
  factor = 3, 
  nfeatures = 5
)
plot_weights(
  MOFAobject, 
  view = "methyl", 
  factor = 1, 
  nfeatures = 5
)

#If you are only interested in looking at only the top features you can use the plotTopWeights function.
plot_top_weights(
  MOFAobject, 
  view="Transcriptomic", 
  factor=2
)
plot_top_weights(
  MOFAobject, 
  view="methyl", 
  factor=1
)
plot_top_weights(
  MOFAobject, 
  view="Transcriptomic", 
  factor=3
)
plot_top_weights(
  MOFAobject, 
  view="miRNA", 
  factor=1
)



#using of meta-data to annotate the MOFA factors
plot_data_heatmap(
  MOFAobject, 
  view = "miRNA", 
  factor = 1, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "gender" , show_colnames = FALSE
  
)
#using of meta-data to annotate the MOFA factors
plot_data_heatmap(
  MOFAobject, 
  view = "methyl", 
  factor = 1, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "gender" , show_colnames = FALSE
  
)

plot_data_heatmap(
  MOFAobject, 
  view = "miRNA", 
  factor = 1, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "sample_type",  show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "methyl", 
  factor = 1, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "sample_type",  show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "miRNA", 
  factor = 1, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "tumor_tissue_site", show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "methyl", 
  factor = 1, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "tumor_tissue_site", show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "miRNA", 
  factor = 1, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "vital_status", show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "Transcriptomic", 
  factor = 2, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "vital_status", show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "methyl", 
  factor = 1, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "vital_status", show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "Transcriptomic", 
  factor = 3, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "vital_status", show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "Transcriptomic", 
  factor = 2, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "sample_type", show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "Transcriptomic", 
  factor = 3, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "sample_type", show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "Transcriptomic", 
  factor = 2, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "gender", show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "Transcriptomic", 
  factor = 3, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "gender", show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "Transcriptomic", 
  factor = 2, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "tumor_tissue_site", show_colnames = FALSE
)
plot_data_heatmap(
  MOFAobject, 
  view = "Transcriptomic", 
  factor = 3, 
  features = 1:5, 
  show_rownames = TRUE, denoise = T, annotation_samples = "tumor_tissue_site", show_colnames = FALSE
)


#######################################################################################################################################
