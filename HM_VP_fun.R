#Heatmap and Volcano plot functions
# To load functions use :
#source("//contrebasse.transimmunom.org/Users/lch/Documents/Covivac stim/R/HM_VP_fun.R") 

library(tidyverse)
library("colorspace")
library("ggplot2")
library("ggrepel")
library("ggh4x")

HM <- function(data, pv.th = 0.05, gradient.th = 1) {

#select only parameters in which at least one comparison  is significant  
  selected <- unique(data[data$pvalues < pv.th, ]$param)
  data <- data[data$param %in% selected, ]

  plot <- ggplot() +
    geom_tile(data = data, aes(x = Time_point, y = param, fill = log2fc), color = "black") +
    geom_point(data = data[data$pvalues < 0.05, ], aes(x = Time_point, y = param), shape = 21, fill = "red") +
    scale_fill_gradient2(low = "orange", mid = "black", high = "blue", breaks = seq(-gradient.th, gradient.th, 1), limits = c(-gradient.th - 0.001, gradient.th + 0.001), guide = guide_colourbar(title.theme = element_text(size = 5, margin = margin(0.15, 0, 0, 0, "cm")), label.theme = element_text(size = 4), title = "log2(fold-change)", order = 1, title.position = "top", title.hjust = 0.5, title.vjust = 0.5, direction = "horizontal"), oob = scales::squish) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    xlab("") +
    ylab("") +
    theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.background = element_rect(fill = "black"),
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

VP <- function(data, pv.th = 0.05) {

#evolution of the parameters. At least 10% of change considered significant : log2(1.1)=0.14  
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