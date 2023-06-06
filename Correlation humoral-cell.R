# Louise CHORFI, 17.05.23
# correlation humoral - cellular responses

setwd("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell stim")
rm(list = ls(all = TRUE))
library(tidyverse)
library(colorspace)
library(ggpubr)


#Import data
source("./Covivac stim/R/load_cell_stim_datas.R")
source("./Covivac stim/R/load_humoral_datas.R")

data_P_wt <- subset(data_P, data_P$condition == "wt")
rm(TEMP,TEMP_N, TEMP_P,TEMP_char,treatments, exploitable,data_P, data_N )

#----------- select data of interest ----------
graph_data_wt =  merge(humoral, data_P_wt)
graph_data_wt$P_IFNg.1_CD8 <- graph_data_wt$P_IFNg.1.TNFa.0_CD8 + graph_data_wt$P_IFNg.1.TNFa.1_CD8
graph_data_wt$P_IFNg.1_CD4 <- graph_data_wt$P_IFNg.1.TNFa.0_CD4 + graph_data_wt$P_IFNg.1.TNFa.1_CD4

#J0 not taken in consideration as they should not have humoral nor cellular response
graph_data_wt <- graph_data_wt %>% subset(graph_data_wt$time_point != "J0")

# ----- correlation S UA------
correlation = cor(x = graph_data_wt$Titre_anti_S_UA, y = graph_data_wt$P_IFNg.1_CD8, method = "spearman")
print(correlation)

graph_corr <- graph_data_wt %>%
  ggplot(aes(x = Titre_anti_S_UA, y = P_IFNg.1_CD8))+ 
  geom_point(pch = 21, alpha = 0.9, aes(fill = time_point), size = 2)+
  geom_smooth()+
  geom_text(x =1000, y = 1.95, hjust = 0, nudge_x = 0.05,
            label = paste0("Spearman correlation coef: r=", round(correlation, 2))) +
  ggtitle("Correlation between CD8+ activated cells after vaccination and humoral response")+
  xlab("Anti Spike Ab AU titer")+
  ylab("% of IFNg+ CD8 Tcells")+
  scale_fill_brewer("Time point", palette = "Blues")

graph_corr = graph_corr + coord_cartesian(ylim = c(-1,2))  

graph_corr

ggsave(filename = "Corr spearman, IFNg - spike titer UA, CD8, WT.png",
       path = "./plots/",
       plot = graph_corr)


# ------- differential analysis T cell resp between seropos and seroneg -------
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
pj = position_jitter(width = 0.2, height=0, seed = 9)
graph_data_wt = graph_data_wt [-31,]
time_point = c("J0", "J28","J56", "J84")

p_graph = data_frame()


for (t in time_point) {
  statres <-  compare_means( P_IFNg.1_CD4 ~ COV19_IgG_anti_S,  
                             method = "wilcox.test", paired = F,
                             data=graph_data_wt[graph_data_wt$time_point == t,],  
                             rm.na = T,
                             symnum.args = symnum.args)
  statres  
  
  boxp <- graph_data_wt [graph_data_wt$time_point == t,] %>%
    ggplot(aes(x = COV19_IgG_anti_S, y = P_IFNg.1_CD4))+ 
    geom_point(alpha = 0.4, position = pj)+
    geom_boxplot(alpha = 0.5, outlier.alpha = 1)+
    stat_summary(fun=mean, geom="point", size= 3, 
                 aes(group= COV19_IgG_anti_S, pch = "mean"))+
    ggtitle("% of IFNg+ CD4+ Tcell in patients with positive or negative Antibody response")+
    xlab("Détection d'IGg Anti Spike")+
    ylab("% of IFNg+ CD4 Tcells") +
    coord_cartesian(ylim = c(0,3)) +
    geom_bracket(xmin = 1, xmax = 2, y = 2.5,label = statres$p.signif)
    
    boxp
  # ggsave(filename = paste0(t," - Boxplot, CD4 stim IFNg- IgG resp.png"),
  #        path = "./plots/",
  #        plot = boxp)
  
  pdf(file = paste0(t,"CD4, plot.pdf"), width = 3, height = 6)
  plot(boxp)
  dev.off()
}
