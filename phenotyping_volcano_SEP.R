rm(list = ls(all = TRUE))
graphics.off()
setwd("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - rpz graph/")


library(tidyverse)
library("colorspace")
library("ggplot2")
library("ggrepel")
library("ggh4x")

#load volcano plot function : 
source("//contrebasse.transimmunom.org/Users/lch/Documents/Covivac stim/R/HM_VP_fun.R") 

# load data and plot VP

# ------ All treatments merged together (fig. not shown) ------------
# txt file of pvalues
data <- read.delim("./02_data/pvalues_all.txt", sep = "\t")
#keep comparison from baseline
data <- data[data$visit_1 == "J0", ]

#set column for data to plot
data$Time_point = data$visit_2
data$log2fc <- log2(data$foldchange)
data$log10pb <- -log10(data$pvalues)
data$param = gsub("P", "%", data$param)

#plot
plot_HM <- HM(data)
plot_HM
plot_VP <- VP(data)
  #look when changes are occuring with facet
plot_VP <- plot_VP + facet_wrap(~Time_point, ncol=1)
plot_VP

#save
pdf(paste0("./03_figures/pvalues_all_HM.pdf"), height = 11.69, width = 8.27)
plot(plot_HM)
dev.off()
pdf(paste0("./03_figures/pvalues_all_VP.pdf"), height = 11.69, width = 8.27)
plot(plot_VP)
dev.off()


# ------ Patients grouped by treatments -------
data <- read.delim("./02_data/pvalues_traitements_all.txt", sep = "\t")
data <- data[data$visit_1 == "J0", ]
data$Time_point = data$visit_2
data$log2fc <- log2(data$foldchange)
data$log10pb <- -log10(data$pvalues)
data$param = gsub("P", "%", data$param)
plot_HM <- HM(data)
plot_HM <- plot_HM + facet_wrap(~treatment, scale = "free")
plot_VP <- VP(data)
plot_VP <- plot_VP + facet_grid(treatment~Time_point)
pdf(paste0("./03_figures/pvalues_traitements_all_HM.pdf"), height = 11.69, width = 8.27)
plot(plot_HM)
dev.off()
pdf(paste0("./03_figures/pvalues_traitements_all_VP.pdf"), height = 11.69, width = 8.27)
plot(plot_VP)
dev.off()

# ----------- Patients under anti-CD20 --------------
data = read.delim("./02_data/SEP - volcano_Wilcox Adjpvalues_Bdep.txt",sep="\t")
data <- data[data$visit_1 == "J0", ]
data$Time_point = data$visit_2
data$log2fc <- log2(data$foldchange)
data$log10pb <- -log10(data$pvalues)
data$param = gsub("P", "%", data$param)
plot_HM <- HM(data)
plot_HM
plot_VP <- VP(data)
plot_VP <- plot_VP + facet_wrap(~Time_point, ncol=1)
plot_VP
pdf(paste0("./03_figures/pheno_Bdep_HM.pdf"), height = 11.69, width = 8.27)
plot(plot_HM)
dev.off()
pdf(paste0("./03_figures/pheno_Bdep_VP.pdf"), height = 11.69, width = 8.27)
plot(plot_VP)
dev.off()

#------------ Paitent with other treatments ----------
data = read.delim("./02_data/SEP - volcano_Wilcox Adjpvalues_Other DMT.txt",sep="\t")
data <- data[data$visit_1 == "J0", ]
data$Time_point = data$visit_2
data$log2fc <- log2(data$foldchange)
data$log10pb <- -log10(data$pvalues)
data$param = gsub("P", "%", data$param)
plot_HM <- HM(data)
plot_HM
plot_VP <- VP(data)
plot_VP <- plot_VP + coord_cartesian(xlim = c(-2,6), ylim = c(0,5))+
  facet_wrap(~Time_point, ncol=1)
plot_VP
pdf(paste0("./03_figures/pheno_Other_DMT_HM.pdf"), height = 11.69, width = 8.27)
plot(plot_HM)
dev.off()
pdf(paste0("./03_figures/pheno_Other_DMT_VP.pdf"), height = 11.69, width = 8.27)
plot(plot_VP)
dev.off()