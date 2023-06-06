# Louise CHORFI, 28/05/23
# MDS and PCA for phenotyping data 
setwd("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - rpz graph/")
rm(list = ls(all = TRUE))


#colors palette for  MDS function
color_palette <- c("#f54242", "#b83030", "#912626",
                            "#f5ad42", "#ab792e","#b37e32",
                            "#7bf542", "#57b02e", "#33661b",
                            "#42d7f5", "#3abcd6","#2e9ab0", "#1a5561",
                            "#7842f5", "#542dad", "#341b6e",
                            "#f542c2", "#d439a8", "#b0308c")
                            
#### loading datas and functions ####
data <- read.delim("./02_data/matrice pheno.txt", sep = "\t", check.names = FALSE)
samples <- paste0(data$treatment, "_", data$Visit, "_", data$Sample_name)
data$treatment <- NULL
data$Visit <- NULL
data$Sample_name <- NULL
data <- t(data)
data <- data.frame(data, check.names = FALSE)
colnames(data) <- samples

source("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell stim/Covivac stim/R/MDS_PCA_fun.R")


# Plots and saves
plot_MDS <- MDS(data)
plot_MDS

pdf(paste0("./03_figures/phenotyping_MDS.pdf"), height = 11.69, width = 8.27)
plot_MDS <- MDS(data) 
plot_MDS <- plot_MDS+facet_wrap(~treatment)
plot(plot_MDS)
plot_MDS <- MDS(data, vec = c(1, 2))
plot(plot_MDS)
dev.off()

pdf(paste0("./03_figures/phenotyping_PCA.pdf"), height = 11.69, width = 8.27)
pca <- PCA(data)
plot(pca$ind)
plot(pca$var)
dev.off()

