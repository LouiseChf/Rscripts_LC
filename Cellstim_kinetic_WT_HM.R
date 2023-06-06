#Louise CHORFI, 12/05/23
# Stats and Heat maps of WT peptide T cell stimulation 
rm(list = ls(all = TRUE))
setwd("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell stim")
library(tidyverse)
library(colorspace)
library(ggpubr)


#------Importer les données de excel : ------
source("./Covivac stim/R/load_cell_stim_datas.R")


#----- select data of interest ------
data <- data.frame()
data <- subset(x = data_P,
              subset = (condition == "wt"))

data <- data[c("SAMPLE_ID","time_point","patient_ID","condition","treatment",
                           "P_IFNg.1.TNFa.0_CD4","P_IFNg.1.TNFa.1_CD4",
               "P_IFNg.1.TNFa.0_CD8","P_IFNg.1.TNFa.1_CD8")]

data$P_IFNg.1_CD4 <- data$P_IFNg.1.TNFa.0_CD4 + data$P_IFNg.1.TNFa.1_CD4
data$P_IFNg.1_CD8 <- data$P_IFNg.1.TNFa.0_CD8 + data$P_IFNg.1.TNFa.1_CD8

View(data)

rm(TEMP, TEMP_N, TEMP_P, data_P, data_N, TEMP_char,j, exploitable)
data <- data[!grepl("newag", data$SAMPLE_ID), ]


#------ normalizer à J0 ------
numeric_cols <- sapply(data, is.numeric)
data[, numeric_cols] <- data[, numeric_cols] + 1
rm(numeric_cols)

norm_data = data_frame()
data <- data %>% group_by(patient_ID)
for(n in unique(data$patient_ID)){
  rowJ0_exists <- any(data$time_point == "J0" & data$patient_ID == n)
  if ( rowJ0_exists == FALSE ){
    next
  }
  J0_CD4 <-data$P_IFNg.1_CD4 [data$patient_ID == n & data$time_point == "J0"]
  J0_CD8 <-data$P_IFNg.1_CD8 [data$patient_ID == n & data$time_point == "J0"]
  temp <- data_frame(patient_ID = n, treatment = data$treatment[data$patient_ID == n],
                     time_point = data$time_point[data$patient_ID == n],
                     P_IFNg.1_CD4  = data$P_IFNg.1_CD4 [data$patient_ID == n]/J0_CD4,
                     P_IFNg.1_CD8  = data$P_IFNg.1_CD8 [data$patient_ID == n]/J0_CD8)
  norm_data <- rbind(norm_data, temp)
}

rm(J0_CD4, J0_CD8, rowJ0_exists)

# ---------------- stats ----------------------------
my_comparisons <- list(c("J0", "J56"), c("J0", "J28"), c("J28","J56"), c("J0","J84"), c("J28","J84"), c("J56","J84"))
treatments <- c("healthy","Ocrelizumab","Teriflunomide","Fingolimod","Rituximab","Natalizumab","Dimethylfumarate")
param <- c("P_IFNg.1_CD4", "P_IFNg.1_CD8" )

for (t in treatments){
  assign(x = t, value = (norm_data %>% subset(subset = norm_data$treatment == t )))
}

t_pvalues = data.frame()

#### Treatment grouped####
for (t in treatments){
  for (p in param){
    for ( v in my_comparisons){
      
      temp_v1 <- subset(get(t), subset = time_point == v[1], select = c("patient_ID", p))
      names(temp_v1) <- c("Sample_name", "Values")
      
      temp_v2 <- subset(get(t), time_point== v[2], select = c("patient_ID", p))
      names(temp_v2) <- c("Sample_name", "Values")
      
      temp <- merge(temp_v1, temp_v2, by = "Sample_name", incomparables = NA)
      temp <- na.omit(temp)
      
      names(temp) <- c("Sample_name", v)
      #print(head(temp))
      # stat
      if (nrow(temp) == 0){
        wilcoxon <- data.frame("p.value" = NA)
      } else {
        wilcoxon <- wilcox.test(x = as.numeric(temp[[2]]), y = as.numeric(temp[[3]]),paired = T, correct = T)
      }
      temp_pvalues <- data.frame("treatment" = t,
                                 "param" = p, "time_point_1" = v[[1]] , 
                                 "time_point_2" = v[[2]], "pvalues"= wilcoxon$p.value,
                                 "foldchange"= mean(temp_v2$Values, na.rm = T)/mean(temp_v1$Values, na.rm =T)) 
      t_pvalues <- rbind(t_pvalues, temp_pvalues)
      
    }
  }
}
t_pvalues$adj_pBH <- p.adjust(t_pvalues$pvalues, method = "BH" )
t_pvalues <- na.omit(t_pvalues)
write.table( t_pvalues, file = "./pvalues_cinetique_traitements.txt", sep = "\t", row.names = FALSE)

#### All merged ####
for (p in param){
  for ( v in my_comparisons){

    temp_v1 <- subset(norm_data, subset = time_point == v[1], select = c("patient_ID", p))
    names(temp_v1) <- c("Sample_name", "Values")

    temp_v2 <- subset(norm_data, time_point== v[2], select = c("patient_ID", p))
    names(temp_v2) <- c("Sample_name", "Values")

    temp <- merge(temp_v1, temp_v2, by = "Sample_name", incomparables = NA)
    temp <- na.omit(temp)

    names(temp) <- c("Sample_name", v)
    #print(head(temp))
    # stat
    if (nrow(temp) == 0){
      wilcoxon <- data.frame("p.value" = NA)
    } else {
      wilcoxon <- wilcox.test(x = as.numeric(temp[[2]]), y = as.numeric(temp[[3]]),paired = T, correct = T)
    }
    temp_pvalues <- data.frame("treatment" = "all",
                               "param" = p, "time_point_1" = v[[1]] ,
                               "time_point_2" = v[[2]], "pvalues"= wilcoxon$p.value,
                               "foldchange"= mean(temp_v2$Values, na.rm = T)/mean(temp_v1$Values, na.rm =T))
    t_pvalues <- rbind(t_pvalues, temp_pvalues)

  }
}

t_pvalues$adj_pBH <- p.adjust(t_pvalues$pvalues, method = "BH" )
t_pvalues <- na.omit(t_pvalues)
write.table( t_pvalues, file = "./pvalues_cinetique_merged.txt", sep = "\t", row.names = FALSE)

rm(n, p, t, v,temp_v1, temp_v2,temp_pvalues, temp, wilcoxon)

#### Bdep / Other ####

treatments = c("aCD20","Other")
Other <- norm_data %>% subset(norm_data$treatment != "Ocrelizumab" & norm_data$treatment != "Rituximab" & norm_data$treatment != "healthy")
aCD20 <- norm_data %>% subset(norm_data$treatment == "Ocrelizumab" | norm_data$treatment == "Rituximab")

for (t in treatments){
  for (p in param){
    for ( v in my_comparisons){
      
      temp_v1 <- subset(get(t), subset = time_point == v[1], select = c("patient_ID", p))
      names(temp_v1) <- c("Sample_name", "Values")
      
      temp_v2 <- subset(get(t), time_point== v[2], select = c("patient_ID", p))
      names(temp_v2) <- c("Sample_name", "Values")
      
      temp <- merge(temp_v1, temp_v2, by = "Sample_name", incomparables = NA)
      temp <- na.omit(temp)
      
      names(temp) <- c("Sample_name", v)
      #print(head(temp))
      # stat
      if (nrow(temp) == 0){
        wilcoxon <- data.frame("p.value" = NA)
      } else {
        wilcoxon <- wilcox.test(x = as.numeric(temp[[2]]), y = as.numeric(temp[[3]]),paired = T, correct = T)
      }
      temp_pvalues <- data.frame("treatment" = t,
                                 "param" = p, "time_point_1" = v[[1]] , 
                                 "time_point_2" = v[[2]], "pvalues"= wilcoxon$p.value,
                                 "foldchange"= mean(temp_v2$Values, na.rm = T)/mean(temp_v1$Values, na.rm =T)) 
      t_pvalues <- rbind(t_pvalues, temp_pvalues)
      
    }
  }
}
t_pvalues$adj_pBH <- p.adjust(t_pvalues$pvalues, method = "BH" )
t_pvalues <- na.omit(t_pvalues)
write.table( t_pvalues, file = "./pvalues_cinetique_Bdep_Other.txt", sep = "\t", row.names = FALSE)

# ------------- Heatmap ---------------------
source("./Covivac stim/R/HM_VP_fun.R")

####HM kinetics per treatments####
dataHM <- read.delim("./pvalues_cinetique_traitements.txt",sep="\t")
dataHM <- dataHM[dataHM$time_point_1 == "J0", ]
dataHM$Time_point = dataHM$time_point_2
dataHM$log2fc <- log2(dataHM$foldchange)
dataHM$log10pb <- -log10(dataHM$pvalues)

plotHM_t <- HM(dataHM)
plotHM_t
plotHM_t <- plotHM_t + facet_wrap(~dataHM$treatment, scale = "free")
plotHM_t
pdf(paste0("./pvalues_cinetique_traitements_HM.pdf"), height = 11.69, width = 8.27)
plot(plot_HM)
dev.off()


#### HM kinetics Bdep vs Other ####
dataHM <- read.delim("./pvalues_cinetique_Bdep_Other.txt",sep="\t")
dataHM <- dataHM[dataHM$time_point_1 == "J0", ]
dataHM$Time_point = dataHM$time_point_2
dataHM$log2fc <- log2(dataHM$foldchange)
dataHM$log10pb <- -log10(dataHM$pvalues)

plotHM_t <- HM(dataHM)
plotHM_t
plotHM_t <- plotHM_t + facet_wrap(~treatment, scale = "free")
plotHM_t
pdf(paste0("./pvalues_cinetique_Bdep_Other_HM.pdf"), height = 11.69, width = 8.27)
plot(plot_HM)
dev.off()

####HM kinetics all merged####
dataHM <- read.delim("./pvalues_cinetique_merged.txt",sep="\t")
dataHM <- dataHM[dataHM$time_point_1 == "J0", ]
dataHM$Time_point = dataHM$time_point_2
dataHM$log2fc <- log2(dataHM$foldchange)
dataHM$log10pb <- -log10(dataHM$pvalues)

plotHM_t <- HM(dataHM)
plotHM_t
plotHM_t <- plotHM_t + facet_wrap(~dataHM$treatment, scale = "free")
plotHM_t
pdf(paste0("./pvalues_cinetique_merged_HM.pdf"), height = 11.69, width = 8.27)
plot(plot_HM)
dev.off()

