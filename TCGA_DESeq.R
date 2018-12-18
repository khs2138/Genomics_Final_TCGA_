################################################################################
# Kayla Schiffer
# Intro to Genonimc Information Science and Technology
# Final Project
# December 17, 2018
################################################################################

################################################################################
# TCGA Biolinks
# Download raw count data 
################################################################################

# Install/update TCGA biolinks
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

# Load required libraries
library(TCGAbiolinks)
library(dplyr)
library(DT)

# Get gene expression data aligned against hg19.
query <- GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts",
  legacy = FALSE)

# Download the data
GDCdownload(query, method = "api")
data <- GDCprepare(query)

# Gather expression and column data
ge <- assay(data)
col_data <- colData(data)

################################################################################

# Normalize gene expression data

################################################################################

# Load required lirarbies
library(DESeq2)

# Add a column smoking, which indicates T/F, extracted from where  packs per day column > 0
col_data$smoking <- !(is.na(col_data$cigarettes_per_day))

dds <- DESeqDataSetFromMatrix(
  countData = ge, 
  colData = col_data,
  design = ~ smoking)

# Filter out low expressed genes, defined by where the total count of that gene is < 500
keep <- rowSums(counts(dds)) >= 500
dds <- dds[keep,]

# Normalize by library size 
dds <- DESeq(dds)

################################################################################

# Identify cluseters in gene expression data

################################################################################

# Update/install if needed
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")

# Load required library
library(ConsensusClusterPlus)

# Find the top 1,000 most variable genes
ge <- counts(dds, normalized = T )
gene_vars <- apply(ge, 1, var)
top_var_genes_idx <- order(gene_vars, decreasing = TRUE)[1:1000]

# Get the sub matrix of the expression matrix by the top most variable genes 
x <- ge[top_var_genes_idx, ]

# Transpose matrix for row-by-row PCA on top most variable genes
pca_res <- PCA(t(x), graph=FALSE)
plot.PCA(title = "PCA without K Clusters", pca_res, label="none")

# Run Consensus Clustering to determine the best choice of k for cluster analysis
title="consensusClustering"
results = ConsensusClusterPlus(
  d = x, maxK=15, reps=50, pItem=0.8, pFeature=1,
  title=title, clusterAlg="hc", distance="pearson",
  seed=1262118388.71279, plot="png")
# See report for figures from consensus clustering
# Best result: K = 10

# Run K-means clustering
km <- kmeans(t(x), centers=10)
names(km)
km

# Redo PCA with K clusters
plot.PCA(title = "PCA with K Clusters", pca_res, label="none", col.ind=km$cluster)

###################################################################################################

# Heatmap

###################################################################################################

# Load required libraries
library(gplots)
library(RColorBrewer)

# Sort km$cluster vector by clusters 
km_order <- order(km$cluster)

# Reorder x columns by km_sort, i.e. by cluster
clust <- x[,km_order]

# Set coloring using RColorBrewer
my_palette <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(n = 100)
cluster_cols <- rainbow(G)[group]

# Generate a heatmap using all 1,000 genes
heatmap.2(as.matrix(clust), scale = "row",
          main = "Clustered Heatmap",  # heat map title
          density.info = "none",       # turns off density plot inside color legend
          trace = "none",              # turns off trace lines inside the heat map
          margins = c(12,9),           # widens margins around plot
          dendrogram = "none",
          Colv = F,       
          Rowv = T, 
          col = my_palette, 
          ColSideColors = cluster_cols)

# Identify genes for clusters
G <- 10
N <- ncol(clust)
M <- 1000
group <- factor((km$cluster[order(km$cluster)]))

# Perform a T test to find genes significant for each cluster vs. all others
pv.sig <- lapply(levels(group), function(g){
  ingroup <- names(group)[group %in% g]
  outgroup <- names(group)[!(group %in% g)]
  pv <- sapply(1:M, function(i) {
    #t.test(clust[i,ingroup], clust[i,outgroup], alternative='greater')$p.value
    t.test(clust[i,ingroup], clust[i,outgroup])$p.value
  })
  names(pv) <- rownames(clust)
  pv.sig <- names(pv)[pv < 0.05/M/G] 
})

# Regenerate heatmap 
heatmap.2(as.matrix(log10(clust[unique(unlist(pv.sig)),] + 1)), 
          scale = "row",
          main = "Clustered Heatmap Globally Distinct Genes",
          density.info = "none",
          trace = "none", 
          #margins = c(12,9),
          dendrogram = "none",
          Colv = F,       
          Rowv = F, 
          labCol=FALSE, labRow=FALSE,
          col = my_palette, 
          ColSideColors = cluster_cols)

################################################################################

# Survival Analysis

################################################################################

# Load required libraries
library(survival)
library(survminer)
library(dplyr)

# Generate survival data matrix
surv_data <- data.frame(
  id = col_data$sample,
  cluster = km$cluster,
  time = rep(0, length(col_data$sample)),
  event = rep(0, length(col_data$sample))
)
for(i in 1:nrow(col_data)){ # time to event
  if(is.na(col_data$days_to_death[i])){
    surv_data$time[i] <- col_data$days_to_last_follow_up[i]
  } else {
    surv_data$time[i] <- col_data$days_to_death[i]
  }
}
surv_data$event <- 0 + !(is.na(col_data$days_to_death))   # event 1=death


# Get overall survival for all groups
surv_object <- Surv(time = surv_data$time, event = surv_data$event)
fit <- survfit(surv_object ~ cluster, data = surv_data)

# Plot survival
ggsurvplot(fit, data = surv_data, pval = TRUE)

# Plot only top and bottom cluster
top_surv_data <- surv_data[surv_data$cluster %in% c(10,3),]
surv_object_2 <- Surv(time = top_surv_data$time, event = top_surv_data$event)
fit2 <- survfit(surv_object_2 ~ cluster, data = top_surv_data)
ggsurvplot(fit2, data = top_surv_data, pval = TRUE)
