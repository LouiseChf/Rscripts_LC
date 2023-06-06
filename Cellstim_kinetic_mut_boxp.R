# Louise CHORFI, 05.05.23
# Boxplot representation, kinetics of T cell stimulation with Mut peptides

setwd("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell stim")
rm(list = ls(all = TRUE))
library(tidyverse)
library(colorspace)
library(ggpubr)


#Import data
source("./Covivac stim/R/load_cell_stim_datas.R")

#-------------- Subset data of interest --------------
graph_data <- data.frame()
graph_data <- subset(x = T1_P,
                     subset = condition == "mut")
names(graph_data)

graph_data <- graph_data[c("SAMPLE_ID","time_point","patient_ID","condition","treatment",
                           "P_IFNg.1.TNFa.0_CD8","P_IFNg.1.TNFa.1_CD8")]

graph_data$P_IFNg.1.TNFa.1_CD8 <- graph_data$P_IFNg.1.TNFa.0_CD8 + graph_data$P_IFNg.1.TNFa.1_CD8
View(graph_data)

rm(TEMP, TEMP_N, TEMP_P, T1_P, T1_N)
graph_data <- graph_data[!grepl("newag", graph_data$SAMPLE_ID), ]


#----------------- Normalize to J0 ------------

numeric_cols <- sapply(graph_data, is.numeric)
graph_data[, numeric_cols] <- graph_data[, numeric_cols] + 1

norm_graph_data = data_frame()
graph_data <- graph_data %>% group_by(patient_ID)
for(n in unique(graph_data$patient_ID)){
  rowJ0_exists <- any(graph_data$time_point == "J0" & graph_data$patient_ID == n)
  if ( rowJ0_exists == FALSE ){
    next
  }
  J0 <-graph_data$P_IFNg.1.TNFa.1_CD8[graph_data$patient_ID == n & graph_data$time_point == "J0"]
  temp <- data_frame(patient_ID = n, treatment = graph_data$treatment[graph_data$patient_ID == n],
                     condition = graph_data$condition[graph_data$patient_ID == n],
                     time_point = graph_data$time_point[graph_data$patient_ID == n],
                     increase = graph_data$increase[graph_data$patient_ID == n],
                     P_IFNg.1.TNFa.1_CD8 = graph_data$P_IFNg.1.TNFa.1_CD8[graph_data$patient_ID == n]/J0)
  norm_graph_data <- rbind(norm_graph_data, temp)
}


#-------------- Split anti-CD20 p. vs Other p. ----------------
norm_Bdep <- norm_graph_data %>% 
  subset(norm_graph_data$treatment == "Ocrelizumab" | norm_graph_data$treatment == "Rituximab")


norm_Other <- norm_graph_data %>%
  subset(norm_graph_data$treatment != "Ocrelizumab" & norm_graph_data$treatment != "Rituximab" ) 

# ------------- generic stat & useful vectors ---------------
summary(graph_data)
# Variances homogeneity : var homogeneous if p > 0.5
bartlett.test(graph_data$P_IFNg.1.TNFa.1_CD8 ~ graph_data$time_point)
#normal distribution : normal if p > 0.5
shapiro.test(graph_data$P_IFNg.1.TNFa.1_CD8)
##  No homogeneity, not normal

my_comparison <- list(c("J0", "J56"), c("J0", "J28"), c("J28","J56"), c("J0","J84"), c("J28","J84"), c("J56","J84"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

#---------------All patients merged --------------- 
# ----- Wilcox -----

p_all = data_frame()
for (comp in my_comparison){
  print(comp)
  statdf <- norm_graph_data %>%
    subset (subset = norm_graph_data$time_point == comp[1] | norm_graph_data$time_point == comp[2]) 
  print((statdf))
  
  for (n in unique(statdf$patient_ID)){
    row1_exists <- sum(statdf$time_point == comp[1] & statdf$patient_ID == n)
    row2_exists <- sum(statdf$time_point == comp[2]& statdf$patient_ID == n)
    if ( row1_exists != row2_exists ){
      statdf <- subset(x=statdf, subset = patient_ID != n)
    }
  }
  row1_exists <- any(statdf$time_point == comp[1])
  row2_exists <- any(statdf$time_point == comp[2])
  if (row1_exists == FALSE | row2_exists == FALSE){
    next
  }
  
  
  statres <- compare_means(P_IFNg.1.TNFa.1_CD8 ~ time_point,  
                           method = "wilcox.test", paired = T,
                           data=statdf,
                           rm.na = T,
                           symnum.args = symnum.args)
  
  temp <- data_frame(t1 = comp[1], t2 = comp[2],
                     pvalue = statres$p, signif = statres$p.signif)
  print(temp)
  p_all = rbind(p_all, temp)
}
# ----- Graph -----
merged_graph <- norm_graph_data %>%
  group_by(patient_ID)%>%
  ggplot(aes(x = time_point, y = P_IFNg.1.TNFa.1_CD8))+ 
  #geom_boxplot(alpha = 0.5, outlier.alpha = 1)+
  geom_line(aes(group = patient_ID, color = increase), linewidth = 1, alpha = 0.3)+
  ggtitle(" Response to mut peptide stimulation at each time point \n all Patients")+
  xlab("Time point")+
  ylab("% of IFN+/ TNFa+ CD8 Tcells")+
  geom_dotplot(binaxis='y', stackdir = "center", stackratio = 0.3, alpha = 0.4, dotsize = 0.8)+
  stat_summary(fun=mean, geom="line", lwd = 1, linetype ='longdash',color = "blue", aes(group=1))+
  stat_summary(fun=mean, geom="point", size= 3, color = "blue", aes(group=1))+
  geom_bracket(xmin = 1, xmax = 2, y.position = 1.2,
               label = p_all$signif[2])+
  geom_bracket(xmin = 1, xmax = 3, y.position = 1.30,
               label = p_all$signif[1])+
  geom_bracket(xmin = 1, xmax = 4, y.position = 1.40,
               label = p_all$signif[4])+
  geom_bracket(xmin = 2, xmax = 3, y.position = 1.55,
               label = p_all$signif[3])+
  geom_bracket(xmin = 2, xmax = 4, y.position = 1.65,
               label = p_all$signif[5])+
  geom_bracket(xmin = 3, xmax = 4, y.position = 1.80,
               label = p_all$signif[6])  

merged_graph

ggsave(filename = "mut, doublepos, CD8, all.png",
       path = "./plots/",
       plot = merged_graph)



# -------------- anti CD20 treated patients , Bdep (fig not shown) ---------
# ----- Wilcox -----

p_Bdep = data_frame()
for (comp in my_comparison){
  print(comp)
  statdf <- norm_Bdep %>%
    subset (subset = norm_Bdep$time_point == comp[1] | norm_Bdep$time_point == comp[2]) 
  print((statdf))
  
  for (n in unique(statdf$patient_ID)){
    row1_exists <- sum(statdf$time_point == comp[1] & statdf$patient_ID == n)
    row2_exists <- sum(statdf$time_point == comp[2]& statdf$patient_ID == n)
    if ( row1_exists != row2_exists ){
      statdf <- subset(x=statdf, subset = patient_ID != n)
    }
  }
  row1_exists <- any(statdf$time_point == comp[1])
  row2_exists <- any(statdf$time_point == comp[2])
  if (row1_exists == FALSE | row2_exists == FALSE){
    next
  }
  
  
  statres <- compare_means(P_IFNg.1.TNFa.1_CD8 ~ time_point,  
                           method = "wilcox.test", paired = T,
                           data=statdf,
                           rm.na = T,
                           symnum.args = symnum.args)
  
  temp <- data_frame(t1 = comp[1], t2 = comp[2],
                     pvalue = statres$p, signif = statres$p.signif)
  print(temp)
  p_Bdep = rbind(p_Bdep, temp)
}
# ----- Graph -----
B_graph <- norm_Bdep %>%
  group_by(patient_ID)%>%
  ggplot(aes(x = time_point, y = P_IFNg.1.TNFa.1_CD8))+ 
  #geom_boxplot(alpha = 0.5, outlier.alpha = 1)+
  geom_line(aes(group = patient_ID, color = increase), linewidth = 1, alpha = 0.3)+
  ggtitle(" Response to mut peptide stimulation at each time point \n Patients treated with aCD20 treatment")+
  xlab("Time point")+
  ylab("% of IFN+/TNFa+ CD8 Tcells")+
  geom_dotplot(binaxis='y', stackdir = "center", stackratio = 0.3, alpha = 0.4, dotsize = 0.8)+
  stat_summary(fun=mean, geom="line", lwd = 1, linetype ='longdash',color = "blue", aes(group=1))+
  stat_summary(fun=mean, geom="point", size= 3, color = "blue", aes(group=1))+
  geom_bracket(xmin = 1, xmax = 2, y.position = 1.2,
               label = p_Bdep$signif[1])+
  geom_bracket(xmin = 1, xmax = 3, y.position = 1.30,
               label = p_Bdep$signif[2])+
  geom_bracket(xmin = 1, xmax = 4, y.position = 1.40,
               label = p_Bdep$signif[3])+
  geom_bracket(xmin = 2, xmax = 3, y.position = 1.55,
               label = "NA")+
  geom_bracket(xmin = 2, xmax = 4, y.position = 1.65,
               label = p_Bdep$signif[4])+
  geom_bracket(xmin = 3, xmax = 4, y.position = 1.80,
               label = p_Bdep$signif[5])  

B_graph

ggsave(filename = "mut, IFNg, CD8, aCD20.png",
       path = "./plots/",
       plot = B_graph)



# -------------- other treatment patients, Other (fig not shown) ---------
# ----- Wilcox -----
p_Other = data_frame()
for (comp in my_comparison){
  print(comp)
  statdf <- norm_Other %>%
    subset (subset = norm_Other$time_point == comp[1] | norm_Other$time_point == comp[2]) 
  print((statdf))
  
  for (n in unique(statdf$patient_ID)){
    row1_exists <- sum(statdf$time_point == comp[1] & statdf$patient_ID == n)
    row2_exists <- sum(statdf$time_point == comp[2]& statdf$patient_ID == n)
    if ( row1_exists != row2_exists ){
      statdf <- subset(x=statdf, subset = patient_ID != n)
    }
  }
  row1_exists <- any(statdf$time_point == comp[1])
  row2_exists <- any(statdf$time_point == comp[2])
  if (row1_exists == FALSE | row2_exists == FALSE){
    next
  }
  
  
  statres <- compare_means(P_IFNg.1.TNFa.1_CD8 ~ time_point,  
                           method = "wilcox.test", paired = T,
                           data=statdf,
                           rm.na = T,
                           symnum.args = symnum.args)
  
  temp <- data_frame(t1 = comp[1], t2 = comp[2],
                     pvalue = statres$p, signif = statres$p.signif)
  print(temp)
  p_Other = rbind(p_Other, temp)
}


# ----- graph  ---------
O_graph <- norm_Other %>%
  group_by(patient_ID)%>%
  ggplot(aes(x = time_point, y = P_IFNg.1.TNFa.1_CD8))+ 
  #geom_boxplot(alpha = 0.5, outlier.alpha = 1)+
  geom_line(aes(group = patient_ID, color = increase), linewidth = 1, alpha = 0.3)+
  ggtitle(" Response to mut peptide stimulation at each time point \n Patients treated with OtherDMT")+
  xlab("Time point")+
  ylab("% of IFN+/TNFa+ CD8 Tcells")+
  geom_dotplot(binaxis='y', stackdir = "center", stackratio = 0.3, alpha = 0.4, dotsize = 0.8)+
  stat_summary(fun=mean, geom="line", lwd = 1, linetype ='longdash',color = "blue", aes(group=1))+
  stat_summary(fun=mean, geom="point", size= 3, color = "blue", aes(group=1))+
  geom_bracket(xmin = 1, xmax = 2, y.position = 1.2, bracket.shorten = 0.1,
               label = p_Other$signif[1])+
  geom_bracket(xmin = 2, xmax = 3, y.position = 1.30,bracket.shorten = 0.1,
               label = p_Other$signif[3])+
  geom_bracket(xmin = 1, xmax = 3, y.position = 1.40,bracket.shorten = 0.1,
               label = p_Other$signif[2])

O_graph

ggsave(filename = "mut,IFNg, CD8, OtherDMT.png",
       path = "./plots/",
       plot = O_graph)
