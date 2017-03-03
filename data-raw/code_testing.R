
library(cytofWorkflow)

data_dir <- system.file("extdata", package = "cytofWorkflow")

## Load metadata
md <- read.table(file.path(data_dir, "cytofWrokflow_metadata.txt"), header = TRUE, as.is = TRUE)

## Add colors 
md$colors <- interaction(md$condition2, md$condition1, sep = "_", lex.order = TRUE)
md$colors <- factor(md$colors, labels = c("#B17BA6", "#882E72", "#7BAFDE", "#1965B0"))

color_samples <- as.character(md$color)
names(color_samples) <- md$short_name

## Load panel

panel <- read.table(file.path(data_dir, "cytofWrokflow_panel.txt"), header = TRUE, as.is = TRUE)

## Read the raw FCS files in
library(flowCore)

file_names <- file.path(data_dir, md$file_name)

fcs_raw <- lapply(file_names, read.FCS)

## Find which markers to use
fcs_colnames <- colnames(fcs_raw[[1]])
fcs_isotope <- as.numeric(gsub("[[:alpha:]]", "", fcs_colnames))

mm <- match(panel$Isotope, fcs_isotope)

panel$fcs_colnames <- fcs_colnames[mm]

keep <- panel$fcs_colnames[panel$Use == 1]

## Extract expression + arcsineh normalization

expr <- lapply(fcs_raw, function(x){
  e <- exprs(x)[, keep]
  asinh(e / 5)
})

expr <- do.call("rbind", expr)

colnames(expr) <- panel$Antigen[panel$Use == 1]

## For plotting the heatmaps normalize the expression data to 0-1

expr01 <- expr

rng <- apply(expr01, 2, quantile, p = c(0.01, 0.99))

for(i in 1:ncol(expr01)){
  expr01[,i] <- (expr01[, i] - rng[1, i]) / (rng[2, i] - rng[1, i])
}

expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1


## Generate the sample information

sample_ids <- rep(md$short_name, sapply(fcs_raw, nrow))

## Plot distributions of marker expression 
library(ggplot2)
library(reshape2)

ggdf <- data.frame(sample_id = sample_ids, expr)
ggdf <- melt(ggdf, id.var = "sample_id", value.name = "expression", variable.name = "antigen")

ggp <- ggplot(ggdf, aes(x = expression, color = sample_id)) +
  geom_density() +
  facet_wrap(~ antigen, nrow = 4, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank(), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 4)) +
  scale_color_manual(values = color_samples)

ggp

## Plot number of cells per sample

cell_table <- table(sample_ids)

ggdf <- data.frame(sample_id = names(cell_table), cell_counts = as.numeric(cell_table))

ggp <- ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = sample_id)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = cell_counts), hjust=0.5, vjust=-0.5, size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = color_samples, drop = FALSE) +
  scale_x_discrete(drop = FALSE)

ggp


## The non-redundancy score (NRS).


NRS <- function(x, ncomp = 3){

  pr <- prcomp(x, center = TRUE, scale. = FALSE) 
  
  score <- rowSums(outer(rep(1, ncol(x)), pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
  
  return(score)
  
}


## Split expression per sample
expr_split <- split(data.frame(expr), sample_ids)

nrs_sample <- sapply(expr_split, NRS)

nrs <- rowMeans(nrs_sample, na.rm = TRUE)


## Plot NRS scores for ordered markers
antigens_ordered <- names(sort(nrs, decreasing = TRUE))

ggdf <- data.frame(antigen = factor(names(nrs), levels = antigens_ordered), nrs_score = nrs)

ggplot(ggdf, aes(x = antigen, y = nrs_score)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



## Cell population identification with FlowSOM

library(FlowSOM)
library(ConsensusClusterPlus)

## Run SOM with 100 clusters
set.seed(1234)
fsom <- SOM(expr)

## Meta clustering into 20 clusters with ConsensusClusterPlus
codes <- fsom$codes
nmc <- 20
plot_outdir <- "consensus_plots"

mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100, pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png", clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", distance = "euclidean", seed = 1234)

## Get cluster ids for each cell
fsom_mc <- mc[[nmc]]$consensusClass
clustering <- fsom_mc[fsom$mapping[,1]]

## Plot the SOM codes

## Get code sizes; sometimes not all the codes have mapped cells so they will have size 0
code_sizes <- rep(0, nrow(codes))
names(code_sizes) <- 1:nrow(codes)
code_table <- table(fsom$mapping[, 1]) 
code_sizes[names(code_table)] <- as.numeric(code_table)

## Dimension reduction 
library(cytofkit)

dr <- cytof_dimReduction(codes, method = "tsne", tsneSeed = 1234)
colnames(dr) <- c("dim1", "dim2")

ggdf <- data.frame(dr)
ggdf$cluster <- factor(fsom_mc)
ggdf$size <- code_sizes

ggplot(ggdf,  aes(x = dim1, y = dim2, color = cluster, size = size)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 5)))


## Plotting heatmaps

# ggplot color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

# Get the median expression
expr_median <- aggregate(expr, by = list(clustering), FUN = median)

# Get cluster frequencies
clustering_table <- table(clustering)

# This clustering is based on the markers that were used for the main clustering, and it is used in all the heatmaps
expr_heat <- as.matrix(expr_median[, antigens_ordered])
rownames(expr_heat) <- expr_median[, 1]

cluster_rows <- hclust(dist(expr_heat), method = "average")

labels_row <- paste0(rownames(expr_heat), " (", round(as.numeric(clustering_table)/sum(clustering_table)*100, 2), "%)")
labels_col <- colnames(expr_heat)

# Row annotation for the heatmap
annotation_row <- data.frame(cluster = rownames(expr_heat))
rownames(annotation_row) <- rownames(expr_heat)

color_clusters <- gg_color_hue(nrow(expr_heat))
names(color_clusters) <- rownames(expr_heat)
annotation_colors <- list(cluster = color_clusters)


library(RColorBrewer)
color <- colorRampPalette(brewer.pal(n = 8, name = "YlGnBu"))(100)

library(pheatmap)

pheatmap(expr_heat, color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 8, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors)



































