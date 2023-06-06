graphics.off()
library("tidyverse")
library("ggplot2")
library("ggrepel")
library("ggh4x")
library("MASS")
library("ggarchery")
#To load functions use :
# source("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell stim/Covivac stim/R/MDS_PCA_fun.R")

#------------ MDS function ---------------
MDS <- function(data, vec = 1) {

  #load data and replace NA and 0 to 0.0001
       d <- dist(t(data))
       d[is.na(d)] <- 0.0001
       d[d == 0] <- 0.0001
  #create distance matrix
       fit <- isoMDS(d, k = 2)

       x <- fit$points[, 1]
       y <- fit$points[, 2]

       min_lim <- min(min(x), min(y)) * 1.1
       max_lim <- max(max(x), max(y)) * 1.1

       datai <- data.frame(
              x = x,
              y = y
       )

       datai$id <- colnames(data)
       datai$cond <- as.factor(as.vector(sapply(datai$id, function(x) paste(strsplit(x, "_")[[1]][vec], collapse = "_"))))
       datai$treatment <- as.factor(as.vector(sapply(datai$id, function(x) paste(strsplit(x, "_")[[1]][1], collapse = "_"))))
       datai$timepoint <- as.factor(as.vector(sapply(datai$id, function(x) paste(strsplit(x, "_")[[1]][2], collapse = "_"))))


       png <- ggplot() +
              geom_hline(yintercept = (min_lim + max_lim) / 2, linetype = "dashed", size = 0.2) +
              geom_vline(xintercept = (min_lim + max_lim) / 2, linetype = "dashed", size = 0.2) +
              xlim(min_lim, max_lim) +
              ylim(min_lim, max_lim) +
              geom_point(data = datai, aes(x = x, y = y, fill = cond), size = 2, shape = 21, stroke = 0.2, color = "black") +
              scale_fill_manual(values = color_palette,
                             guide = guide_legend(title.theme = element_blank())) +
              #geom_text_repel(data = datai, aes(x = x, y = y, label = id), size = 2, min.segment.length = 0, seed = 42, segment.color = "black", segment.size = 0.1) +
              coord_fixed() +
              coord_cartesian(xlim = c(min_lim, max_lim), ylim = c(min_lim, max_lim)) +
              xlab(paste0("MDS1")) +
              ylab(paste0("MDS2")) +
              theme_bw() +
              theme(
                     panel.background = element_blank(),
                     panel.border = element_rect(size = 0.1),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     aspect.ratio = 1,
                     axis.title.x = element_text(size = 7),
                     axis.title.y = element_text(size = 7),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.ticks = element_blank(),
                     legend.position = "bottom",
                     plot.title = element_text(size = 7, hjust = 0.5),
                     plot.subtitle = element_text(size = 7, hjust = 0.5, face = "italic")
              )

       png

       return(png)
}

#-------------------- PCA function -------------------
PCA <- function(data, axe1 = 1, axe2 = 2) {
       circleFun <- function(center = c(0, 0), diameter = 1, npoints = 100) {
              r <- diameter / 2
              t <- seq(0, 2 * pi, length.out = npoints)
              x <- center[1] + r * cos(t)
              y <- center[2] + r * sin(t)
              return(data.frame(x = x, y = y))
       }

       res.PCA <- FactoMineR::PCA((t(data)), scale.unit = TRUE, graph = FALSE)

       var_explained <- res.PCA$eig[, 2]

       datai1 <- data.frame((res.PCA$ind$coord[, c(axe1, axe2)]))
       datai1$id <- colnames(data)
       datai1$x <- datai1[, 1]
       datai1$y <- datai1[, 2]
       datai1$cond <- as.factor(as.vector(sapply(datai1$id, function(x) paste(strsplit(x, "_")[[1]][vec], collapse = "_"))))

       data2i <- data.frame((res.PCA$var$coord[, c(axe1, axe2)]))
       data2i$id <- rownames(data)
       data2i$x <- data2i[, 1]
       data2i$y <- data2i[, 2]
       data2i$cond <- as.factor(as.vector(sapply(data2i$id, function(x) paste(strsplit(x, "_")[[1]][vec], collapse = "_"))))

       max <- max(max(abs(datai1$x)), max(abs(datai1$y)))

       pca1 <- ggplot() +
              labs(title = paste0("PCA at the disease level")) +
              geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
              geom_vline(xintercept = 0, linetype = "dashed", size = 0.2) +
              geom_point(data = datai1, aes_string(x = "x", y = "y", fill = datai1$cond), shape = 21, size = 2, stroke = 0.1, color = "black") +
              scale_fill_discrete(#values = color_palette,
                                guide = guide_legend(title.theme = element_blank())) +
              #geom_text_repel(data = datai1, aes(x = x, y = y, label = id), size = 1.5, color = "black", min.segment.length = 0, seed = 42, segment.color = "black", segment.size = 0.1) +
              scale_x_continuous(limits = c(-max, max)) +
              scale_y_continuous(limits = c(-max, max)) +
              xlab(paste0("PC1 (", round(var_explained[1], 2), "%)")) +
              ylab(paste0("PC2 (", round(var_explained[2], 2), "%)")) +
              theme_bw() +
              theme(
                     panel.background = element_blank(),
                     panel.border = element_rect(size = 0.1),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     aspect.ratio = 1,
                     axis.title.x = element_text(size = 7),
                     axis.title.y = element_text(size = 7),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = c(0.995, 0.999),
                     legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                     legend.justification = c(1, 1),
                     legend.text = element_text(size = 5),
                     legend.box = "horizontal",
                     legend.key.height = unit(0.5, "cm"),
                     legend.key.size = unit(0, "cm"),
                     legend.background = element_rect(fill = NA),
                     legend.spacing.x = unit(0.1, "cm"),
                     legend.spacing.y = unit(0, "cm"),
                     legend.box.spacing = unit(0, "cm"),
                     legend.margin = margin(t = 0, r = 0, b = 0, l = 0.15, unit = "cm"),
                     axis.ticks = element_blank(),
                     plot.title = element_text(size = 7, hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5, face = "italic")
              )

       rad_th <- 0
       circle1 <- circleFun(c(0, 0), 2, npoints = 100)
       circle2 <- circleFun(c(0, 0), rad_th + rad_th, npoints = 100)
       data2ic <- data2i[sqrt(data2i$x^2 + data2i$y^2) > rad_th, ]
       data2ic$id <- rownames(data2ic)
       pca2 <- ggplot() +
              labs(title = paste0("PCA correlation circle of cell parameters")) +
              geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
              geom_vline(xintercept = 0, linetype = "dashed", size = 0.2) +
              geom_path(data = circle1, aes(x = -x, y = y), size = 0.2, color = "gray") +
              geom_path(data = circle2, aes(x = -x, y = y), size = 0.2, linetype = "dashed", color = "gray") +
              geom_segment(data = data2ic, aes(x = 0, y = 0, xend = x, yend = y), size = 0.01, linetype = "solid", color = "gray", arrow = grid::arrow(length = unit(0.05, "inches"), type = "closed"), position = ggarchery::position_attractsegment(
                     start_shave = 0.00,
                     end_shave = 0.02,
                     type_shave = "distance"
              )) +
              geom_text_repel(data = data2ic, aes(x = x, y = y, label = id), size = 1.5, color = "black", min.segment.length = 0, seed = 42, segment.color = "black", segment.size = 0.1) +
              geom_point(data = data2ic, aes(x = x, y = y), fill = "grey", shape = 21, size = 1.75, stroke = 0.1, color = "black") +
              scale_fill_discrete(#values = color_palette,
                           guide = guide_legend(title.theme = element_blank())) +
              coord_fixed() +
              scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +
              scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
              xlab(paste0("PC1 (", round(var_explained[1], 2), "%)")) +
              ylab(paste0("PC2 (", round(var_explained[2], 2), "%)")) +
              theme_bw() +
              theme(
                     panel.border = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     aspect.ratio = 1,
                     axis.title.x = element_text(size = 7),
                     axis.title.y = element_text(size = 7),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.background = element_rect(color = NA),
                     axis.ticks = element_blank(),
                     plot.title = element_text(size = 7, hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5, face = "italic"),
                     legend.position = "bottom",
                     legend.box.spacing = unit(0, "cm"),
                     legend.margin = margin(t = 0, r = 0, b = 0, l = 0.15, unit = "cm"),
                     legend.text = element_text(size = 4),
                     legend.title = element_blank(),
                     legend.key.height = unit(0.25, "cm"),
                     legend.key.width = unit(0.25, "cm")
              )


       return(list("ind" = pca1, "var" = pca2, "res" = res.PCA))
}
