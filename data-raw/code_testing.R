
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













