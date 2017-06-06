# CyTOF workflow: differential discovery in high-throughput high-dimensional cytometry datasets

## Installation 

First, you need to install `devtools` if you haven't already:

```
install.packages("devtools")
```

Bioconductor dependencies need to be installed manually (since `install_github` can only install CRAN dependencies automatically):

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocStyle", "flowCore", "FlowSOM", "ConsensusClusterPlus", "limma"))
```


After the Bioconductor dependencies are installed, the package and its CRAN dependencies can be installed from GitHub:

```
devtools::install_github("gosianow/cytofWorkflow")
```

The workflow can be executed by following the instructions in the cytofWorkflow.Rmd file available in the  [vignettes directory](https://github.com/gosianow/cytofWorkflow/blob/master/vignettes/cytofWorkflow.Rmd).







