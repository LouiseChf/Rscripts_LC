rm(list = ls(all = TRUE))
library("plyr")       
library("dplyr")  
library(ggplot2)
library(ggpubr)
library(tidyverse)
library("ggrepel")

setwd("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell numeration panels")
input_dir="//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell numeration panels/VIH - tableau"
output_dir="//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell numeration panels/VIH - analysis"
options(stringsAsFactors = FALSE)
options(error = traceback)
options(show.error.locations = TRUE)

# ------------- Load data ---------------
####numerations ####
num=data.frame()
num <- read.table("./VIH - tableau/extract_2021-07-15_142429.csv", sep=";",dec=",",stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                  check.names = FALSE, fileEncoding = "UTF-8")
names(num)[names(num)=="ID"]="Sample_name"

num <- num[,c(1,2,3,5,7,9,11)]
num$Visit=gsub("M1", "J56", num$Visit)
num$Sample_name <- substr(num$Sample_name, 1, 8)

#### phenotyping data ####
pheno=data.frame()
pheno <- read.table("./VIH - tableau/tableau_007.csv", sep=";",dec=",",stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                 check.names = FALSE, fileEncoding = "UTF-8")
names(pheno)[names(pheno)=="Sample name"]="Sample_name"
pheno$Visit=gsub("M1", "J56", pheno$Visit)

#### total ####
total=merge(num, pheno, by= c("Sample_name", "Visit"), all=TRUE)

total=total[-grep(c("F"),names(total))]
total=total[-grep(c("N"),names(total))]
total=total[-grep(c(".x"), names(total))]
total=total[-grep(c(".y"), names(total))]

#### organizing data####
#Vectors listing stuff
param=names(total)[3:length(names(total))]
param=param[param !="Cal.Factor"]
patients=as.character(unique(total$Sample_name))
visites=as.character(unique(total$Visit))

nomscols=names(total)[-c(1:2)]
nomscols=gsub("P ","% ",nomscols)
nomscols=gsub("P_","% ",nomscols)
nomscols=gsub("N ","Abs. Number of ",nomscols)

# ---------------- normalizing to J0 --------------------
numeric_cols <- sapply(total, is.numeric)
total[, numeric_cols] <- total[, numeric_cols] + 1

norm_total = data_frame()
total <- total %>% group_by(Sample_name)
for(n in unique(total$Sample_name)){
  rowJ0_exists <- any(total$Visit == "J0" & total$Sample_name == n)
  if ( rowJ0_exists == FALSE ){
    next
  }
  J0 <-total[total$Sample_name == n & total$Visit == "J0",3:53]
  temp <- data_frame(Sample_name = n, 
                     Visit = total$Visit[total$Sample_name == n],
                     )
  nbrow <- nrow(total[total$Sample_name == n, 3:53])
  J0_repeated <- do.call(rbind, replicate(nbrow, J0, simplify = FALSE))
  
  temp <- cbind(temp, total[total$Sample_name == n, 3:53]/J0_repeated)

  norm_total <- rbind(norm_total, temp)
}


write.table( total, file = paste0(output_dir,"VIH matrice pheno.txt"), sep = "\t", row.names = FALSE)
write.table( norm_total, file = paste0(output_dir,"VIH matrice pheno_norm.txt"), sep = "\t", row.names = FALSE)


# ---------------- stats ----------------------------
#### pvalues diff between each time point, all treatments merged####
all_pvalues <- data.frame()
my_comparisons <- list(c("J0", "J56"), c("J0", "J28"), c("J28","J56"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

for ( p in param){
  
  for ( v in my_comparisons){
    
    temp_v1 <- subset(total, Visit == v[1], select = c("Sample_name", p))
    names(temp_v1) <- c("Sample_name", "Values")
    
    temp_v2 <- subset(total, Visit == v[2], select = c("Sample_name", p))
    names(temp_v2) <- c("Sample_name", "Values")
    
    temp <- merge(temp_v1, temp_v2, by = "Sample_name", incomparables = NA)
    temp <- na.omit(temp)
    
    names(temp) <- c("Sample_name", v)
    print(head(temp))
    # stat
    wilcoxon <- wilcox.test(x = as.numeric(temp[[2]]), y = as.numeric(temp[[3]]),paired = T, correct = T)
    temp_pvalues <-data.frame("param" = p, "visit_1" = v[[1]] , 
                              "visit_2" = v[[2]], "pvalues"= wilcoxon$p.value,
                              "foldchange"= mean(temp_v2$Values, na.rm = T)/mean(temp_v1$Values, na.rm =T)) 
    all_pvalues <- rbind(all_pvalues, temp_pvalues)
    
  }
}

all_pvalues$adj_pBH <- p.adjust(all_pvalues$pvalues, method = "BH" )
write.table(all_pvalues, file = paste0(output_dir,"VIH pvalues_all.txt"), sep = "\t", row.names = FALSE)

# -------------------- VP --------------------
VP <- function(data, pv.th = 0.05) {
  
  data$dir <- "NS"
  data$dir[data$log2fc > 0.14 & data$pvalues < pv.th] <- "up"
  data$dir[data$log2fc < -0.14 & data$pvalues < pv.th] <- "down"
  
  plot <- ggplot() +
    geom_point(data = data, aes(x = log2fc, y = log10pb, fill = dir), shape = 21, size = 3) +
    geom_text_repel(data = data[data$dir != "NS", ], aes(x = log2fc, y = log10pb, label = param), 
                    size = 3) +
    geom_hline(yintercept = -log10(pv.th), linetype = "dashed", size = 0.2) +
    geom_vline(xintercept = c(-0.14,0.14), linetype = "dashed", size = 0.2) +
    scale_fill_manual(values = c("down" = "green", "NS" = "gray", "up" = "red")) +
    xlab("log2(fold-change)") +
    ylab("-log10(p-value)") +
    scale_x_continuous(breaks=-100:100) +
    theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.background = element_rect(fill = "white"),
      plot.title = element_text(size = 7, hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      axis.text = element_text(size = 6),
      panel.border = element_rect(size = 0.1),
      axis.ticks = element_line(size = 0.05),
      strip.background = element_rect(size = 0.1),
      strip.text = element_text(hjust = 0.5, size = 5, vjust = 1),
      legend.text = element_text(size = 4),
      legend.title = element_blank()
    )
  
  return(plot)
}

data = read.delim("./VIH - analysis/VIH pvalues_all.txt",sep="\t")
data <- data[data$visit_1 == "J0", ]
data$Time_point = data$visit_2
data$log2fc <- log2(data$foldchange)
data$log10pb <- -log10(data$pvalues)
data$param = gsub("P", "%", data$param)
# plot_HM <- HM(data)
# plot_HM
plot_VP <- VP(data)
plot_VP <- plot_VP + facet_wrap(~Time_point, ncol=1)
plot_VP
# pdf(paste0(".pheno_VIH_HM.pdf"), height = 11.69, width = 8.27)
# plot(plot_HM)
# dev.off()
pdf(paste0("./pheno_VIH_VP.pdf"), height = 11.69, width = 8.27)
plot(plot_VP)
dev.off()

