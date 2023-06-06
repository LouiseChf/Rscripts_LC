# Louise CHORFI, 23/05/23
#  Generate statistics tables

setwd("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell numeration panels")
rm(list = ls(all = TRUE))
library("plyr")       
library("dplyr")        
library("readr")
library("stringr")
library(Hmisc)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(rstatix)
library(gridExtra)

input_dir="//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell numeration panels/SEP - tableau par tube"
output_dir="//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell numeration panels/"

# --------------  Load datas ---------------- 
source("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell numeration panels/load_pheno_data.R")
rm(Other, Bdep)

# --------------- Normalizing to J0 ---------------
numeric_cols <- sapply(total, is.numeric)
total[, numeric_cols] <- total[, numeric_cols] + 1

norm_total = data_frame()
total <- total %>% group_by(Sample_name)
for(n in unique(total$Sample_name)){
  rowJ0_exists <- any(total$time_point == "J0" & total$Sample_name == n)
  if ( rowJ0_exists == FALSE ){
    next
  }
  J0 <-total[total$Sample_name == n & total$time_point == "J0",4:54]
  temp <- data_frame(Sample_name = n, 
                     treatment = total$treatment[total$Sample_name == n],
                     time_point = total$time_point[total$Sample_name == n])
  nbrow <- nrow(total[total$Sample_name == n, 4:54])
  J0_repeated <- do.call(rbind, replicate(nbrow, J0, simplify = FALSE))
  
  temp <- cbind(temp, total[total$Sample_name == n, 4:54]/J0_repeated)
  
  norm_total <- rbind(norm_total, temp)
}



write.table( total, file = paste0(output_dir,"matrice pheno.txt"), sep = "\t", row.names = FALSE)
write.table( norm_total, file = paste0(output_dir,"matrice pheno_norm.txt"), sep = "\t", row.names = FALSE)

# ---------------- stats ----------------------------
#------------ pvalues diff between each time point, all treatments merged -------------
all_pvalues <- data.frame()
my_comparisons <- list(c("J0", "J56"), c("J0", "J28"), c("J28","J56"), c("J0","J84"), c("J28","J84"), c("J56","J84"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

for ( p in param){
  
  for ( v in my_comparisons){
    
    temp_v1 <- subset(total, time_point == v[1], select = c("Sample_name", p))
    names(temp_v1) <- c("Sample_name", "Values")
    
    temp_v2 <- subset(total, time_point == v[2], select = c("Sample_name", p))
    names(temp_v2) <- c("Sample_name", "Values")
    
    temp <- merge(temp_v1, temp_v2, by = "Sample_name", incomparables = NA)
    temp <- na.omit(temp)
    
    names(temp) <- c("Sample_name", v)
    print(head(temp))
    # stat
    wilcoxon <- wilcox.test(x = as.numeric(temp[[2]]), y = as.numeric(temp[[3]]),paired = T, correct = T)
    temp_pvalues <-data.frame("param" = p, "time_point_1" = v[[1]] , 
                              "time_point_2" = v[[2]], "pvalues"= wilcoxon$p.value,
                              "foldchange"= mean(temp_v2$Values, na.rm = T)/mean(temp_v1$Values, na.rm =T)) 
    all_pvalues <- rbind(all_pvalues, temp_pvalues)
    
  }
}

all_pvalues$adj_pBH <- p.adjust(all_pvalues$pvalues, method = "BH" )
write.table(all_pvalues, file = paste0(output_dir,"pvalues_all.txt"), sep = "\t", row.names = FALSE)

#------------- pvalues  diff between each time point,per treatment group ----------
#reminder object "treatment" = list os treatments
my_comparisons <- list(c("J0", "J56"), c("J0", "J28"), c("J28","J56"), c("J0","J84"), c("J28","J84"), c("J56","J84"))

for (t in treatments){
  assign(x = t, value = (norm_total %>% subset(subset = norm_total$treatment == t )))
}
treatments <- c("Ocrelizumab","Teriflunomide","Fingolimod","Rituximab","Natalizumab","Dimethylfumarate")
t_pvalues = data.frame()

for (t in treatments){
  
  for ( p in param){
    
    for ( v in my_comparisons){
      
      temp_v1 <- subset(get(t), subset = time_point == v[1], select = c("Sample_name", p))
      names(temp_v1) <- c("Sample_name", "Values")
      
      temp_v2 <- subset(get(t), time_point == v[2], select = c("Sample_name", p))
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
write.table( t_pvalues, file = paste0(output_dir,"pvalues_traitements_all.txt"), sep = "\t", row.names = FALSE)


# ------------ pvalues diff between each traetment, per time point -------------

my_comparisons <- combn(treatments, 2)
my_comparisons <- lapply(seq(ncol(my_comparisons)), function(i) {
  my_comparisons[, i]
})
for (t in time_pointes){
  assign(x = t, value = (norm_total %>% subset(subset = norm_total$time_point == t )))
}
trt_pvalues = data.frame()

for (t in time_pointes){
  
  for ( p in param){
    
    for ( v in my_comparisons){
      
      temp_v1 <- subset(get(t), subset = treatment == v[1], select = c("Sample_name", p))
      names(temp_v1) <- c("Sample_name", "Values")
      temp_v1 <- na.omit(temp_v1)
      
      temp_v2 <- subset(get(t), treatment == v[2], select = c("Sample_name", p))
      names(temp_v2) <- c("Sample_name", "Values")
      temp_v2<- na.omit(temp_v2)
      
      
      if (nrow(temp_v1) == 0 | nrow(temp_v2) == 0 ){
        wilcoxon <- data.frame("p.value" = NA)
      } else {wilcoxon <- wilcox.test(x = as.numeric(temp_v1[[2]]), y = as.numeric(temp_v2[[2]]),paired = F, correct = T)
      }
      
      temp_pvalues <- data.frame("Time_point" = t,
                                 "param" = p, "traitement_1" = v[[1]] , 
                                 "traitement_2" = v[[2]], "pvalues"= wilcoxon$p.value,
                                 "foldchange"= mean(temp_v2$Values, na.rm = T)/mean(temp_v1$Values, na.rm =T)) 
      trt_pvalues <- rbind(trt_pvalues, temp_pvalues)
      
    }
  }
}
trt_pvalues$adj_pBH <- p.adjust(trt_pvalues$pvalues, method = "BH" )
trt_pvalues <- na.omit(trt_pvalues)
write.table( trt_pvalues, file = paste0(output_dir,"pvalues_timepoints_all.txt"), sep = "\t", row.names = FALSE)

