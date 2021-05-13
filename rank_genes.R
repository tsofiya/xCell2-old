library(singscore)
library(GSEABase)

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("singscore")

# You would need to install 'devtools' package first.
install.packages("devtools")

# And install the 'singscore' package from the GitHub repository
# 'singscore' requires these packages to be installed: methods, stats, graphics, ggplot2, ggsci, grDevices,
#  ggrepel, plotly, tidyr, plyr, magrittr, reshape, edgeR, RColorBrewer, Biobase, GSEABase, BiocParallel
devtools::install_github('DavisLaboratory/singscore')
# Set build_vignette = TRUE if would like to browseVignette()


rankData <- rankGenes(tgfb_expr_10_se)
