##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: GxE-Compare stability in 
##                 Midwest vs. non-Midwest
## Date: 2018-09-20
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggpubr)

#### Read in, format, and organize data ####
midwest_slope <- read.csv("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/stability_in_midwest/selectedTraitsSlopes_traitsStandardized.csv", header = T)
midwest_MSE <- read.csv("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/stability_in_midwest/selectedTraitsMSE_traitsStandardized.csv", header = T)
midwest_GGE <- read.csv("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/stability_in_midwest/GGE_stability.csv", header = T)

non_midwest_slope <- read.csv("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/stability_in_non_midwest/selectedTraitsSlopes_traitsStandardized.csv", header = T)
non_midwest_MSE <- read.csv("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/stability_in_non_midwest/selectedTraitsMSE_traitsStandardized.csv", header = T)
non_midwest_GGE <- read.csv("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/stability_in_non_midwest/GGE_stability.csv", header = T)

colnames(midwest_slope) <- c("Genotype", "Anthesis", "Silking",
                             "Plant height", "Ear height", "Plot grain weight",
                             "Ear length", "Ear width", "Kernels per row", "Kernel row number",
                             "Kernel weight", "Kernel area", "Kernel length", "Kernel width", "Kernel thickness")
colnames(midwest_MSE) <- c("Genotype", "Anthesis", "Silking",
                           "Plant height", "Ear height", "Plot grain weight",
                           "Ear length", "Ear width", "Kernels per row", "Kernel row number",
                           "Kernel weight", "Kernel area", "Kernel length", "Kernel width", "Kernel thickness")
colnames(midwest_GGE) <- c("Genotype", "Anthesis", "Silking",
                           "Plant height", "Ear height", "Plot grain weight",
                           "Ear length", "Ear width", "Kernels per row", "Kernel row number",
                           "Kernel weight", "Kernel area", "Kernel length", "Kernel width", "Kernel thickness")

colnames(non_midwest_slope) <- c("Genotype", "Anthesis", "Silking",
                                 "Plant height", "Ear height", "Plot grain weight",
                                 "Ear length", "Ear width", "Kernels per row", "Kernel row number",
                                 "Kernel weight", "Kernel area", "Kernel length", "Kernel width", "Kernel thickness")
colnames(non_midwest_MSE) <- c("Genotype", "Anthesis", "Silking",
                               "Plant height", "Ear height", "Plot grain weight",
                               "Ear length", "Ear width", "Kernels per row", "Kernel row number",
                               "Kernel weight", "Kernel area", "Kernel length", "Kernel width", "Kernel thickness")
colnames(non_midwest_GGE) <- c("Genotype", "Anthesis", "Silking",
                               "Plant height", "Ear height", "Plot grain weight",
                               "Ear length", "Ear width", "Kernels per row", "Kernel row number",
                               "Kernel weight", "Kernel area", "Kernel length", "Kernel width", "Kernel thickness")
#### Compare stability in Midwest vs. non-Midwest (West/East) ####

### Slope
midwest_slope_melt <- melt(midwest_slope)
colnames(midwest_slope_melt) <- c("genotype", "trait", "slope")
midwest_slope_melt$region <- "Midwest"

non_midwest_slope_melt <- melt(non_midwest_slope)
colnames(non_midwest_slope_melt) <- c("genotype", "trait", "slope")
non_midwest_slope_melt$region <- "West/East"

slopes <- bind_rows(midwest_slope_melt, non_midwest_slope_melt)

midwest_vs_non_slope <- ggplot(slopes, aes(x = trait, y = slope, fill = region)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test",
                     label.y = -0.5,
                     aes(label = ..p.signif..)
  ) +
  labs(x = "Trait", y = "Slope", fill = "Region") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black")
        # , legend.position = c(0.15,0.9)
  )
ggsave(plot = midwest_vs_non_slope, filename = "midwest_vs_non_slope.tif", device = "tiff",
       width = 7, height = 5, units = "in")


### Type II
typeII <- slopes

ggplot(typeII, aes(x = trait, y = slope, fill = region)) +
  geom_boxplot() +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
# plant height
typeII_plant_height_midwest <- subset(typeII, typeII$trait == "Plant height" & typeII$region == "Midwest")
median(typeII_plant_height_midwest$slope, na.rm = T)

typeII_plant_height_other <- subset(typeII, typeII$trait == "Plant height" & typeII$region == "West/East")
median(typeII_plant_height_other$slope, na.rm = T)

# ear height
typeII_ear_height_midwest <- subset(typeII, typeII$trait == "Ear height" & typeII$region == "Midwest")
median(typeII_ear_height_midwest$slope, na.rm = T)

typeII_ear_height_other <- subset(typeII, typeII$trait == "Ear height" & typeII$region == "West/East")
median(typeII_ear_height_other$slope, na.rm = T)

# plot weight
typeII_plot_weight_midwest <- subset(typeII, typeII$trait == "Plot grain weight" & typeII$region == "Midwest")
median(typeII_plot_weight_midwest$slope, na.rm = T)

typeII_plot_weight_other <- subset(typeII, typeII$trait == "Plot grain weight" & typeII$region == "West/East")
median(typeII_plot_weight_other$slope, na.rm = T)

# ear length
typeII_ear_length_midwest <- subset(typeII, typeII$trait == "Ear length" & typeII$region == "Midwest")
median(typeII_ear_length_midwest$slope, na.rm = T)

typeII_ear_length_other <- subset(typeII, typeII$trait == "Ear length" & typeII$region == "West/East")
median(typeII_ear_length_other$slope, na.rm = T)


typeII$slope <- abs(typeII$slope - 1)
colnames(typeII)[3] <- "value"
typeII$stability_measure <- "Type II"

midwest_vs_non_typeII <- ggplot(typeII, aes(x = trait, y = value, fill = region)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test",
                     label.y = -0.5,
                     aes(label = ..p.signif..)
  ) +
  labs(x = "Trait", y = "Abs(Slope-1)", fill = "Region") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black")
        # , legend.position = c(0.15,0.9),
  )
ggsave(plot = midwest_vs_non_typeII, filename = "midwest_vs_non_typeII.tif", device = "tiff",
       width = 7, height = 5, units = "in")

### MSE
midwest_MSE_melt <- melt(midwest_MSE)
colnames(midwest_MSE_melt) <- c("genotype", "trait", "MSE")
midwest_MSE_melt$region <- "Midwest"

non_midwest_MSE_melt <- melt(non_midwest_MSE)
colnames(non_midwest_MSE_melt) <- c("genotype", "trait", "MSE")
non_midwest_MSE_melt$region <- "West/East"

MSEs <- bind_rows(midwest_MSE_melt, non_midwest_MSE_melt)

midwest_vs_non_MSE <- ggplot(MSEs, aes(x = trait, y = MSE, fill = region)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test",
                     label.y = 0,
                     aes(label = ..p.signif..)
  ) +
  labs(x = "Trait", y = "MSE", fill = "Region") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black")
        #, legend.position = c(0.15,0.9)
  )
ggsave(plot = midwest_vs_non_MSE, filename = "midwest_vs_non_MSE.tif", device = "tiff",
       width = 7, height = 5, units = "in")


### GGE
midwest_GGE_melt <- melt(midwest_GGE)
colnames(midwest_GGE_melt) <- c("genotype", "trait", "GGE")
midwest_GGE_melt$region <- "Midwest"

non_midwest_GGE_melt <- melt(non_midwest_GGE, id.vars = "Genotype")
colnames(non_midwest_GGE_melt) <- c("genotype", "trait", "GGE")
non_midwest_GGE_melt$region <- "West/East"

GGEs <- bind_rows(midwest_GGE_melt, non_midwest_GGE_melt)

midwest_vs_non_GGE <- ggplot(GGEs, aes(x = trait, y = GGE, fill = region)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test",
                     label.y = -0.1,
                     aes(label = paste0(..p.signif..))
  ) +
  labs(x = "Trait", y = "GGE", fill = "Region") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black")
        #, legend.position = c(0.15,0.9)
  )
ggsave(plot = midwest_vs_non_GGE, filename = "midwest_vs_non_GGE.tif", device = "tiff",
       width = 7, height = 5, units = "in")

colnames(slopes)[3] <- "value"
colnames(MSEs)[3] <- "value"
colnames(GGEs)[3] <- "value"

slopes$stability_measure <- "Slope"
MSEs$stability_measure <- "MSE"
GGEs$stability_measure <- "GGE"
typeII$stability_measure <- "Type II"

stability <- bind_rows(slopes, typeII, MSEs, GGEs)
stability$stability_measure <- factor(stability$stability_measure,
                                      levels = c("Slope", "Type II", "MSE", "GGE"))

midwest_vs_non_all <-
  ggplot(stability, aes(x = trait, y = value, fill = region)) +
  geom_boxplot() +
  scale_fill_manual(values = c("steelblue", "darkred")) +
  stat_compare_means(method = "t.test",
                     label.y = -0.5,
                     aes(label = paste0(..p.signif..), family = "Courier")) +
  labs(x = "Trait", y = "", fill = "Region") +
  facet_grid(stability_measure~., scales = "free_y") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9),
        strip.text = element_text(size = 10), strip.background = element_rect(color = "white", fill = "white")
  )

ggsave(plot = midwest_vs_non_all, filename = "midwest_vs_non_all.tif", device = "tiff",
       width = 9, height = 10, units = "in")


stability2 <- bind_rows(typeII, MSEs, GGEs)
stability2$stability_measure <- factor(stability2$stability_measure,
                                       levels = c("Type II", "MSE", "GGE"))

midwest_vs_non_all2 <-
  ggplot(stability2, aes(x = trait, y = value, fill = region)) +
  geom_boxplot() +
  scale_fill_manual(values = c("steelblue", "darkred")) +
  stat_compare_means(method = "t.test",
                     label.y = -0.5,
                     aes(label = paste0(..p.signif..), family = "Courier")) +
  labs(x = "Trait", y = "Stability value", fill = "Region") +
  facet_grid(stability_measure~., scales = "free_y") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9),
        strip.text = element_text(size = 10), strip.background = element_rect(color = "white", fill = "white")
  )

ggsave(plot = midwest_vs_non_all2, filename = "midwest_vs_non_all2.tif", device = "tiff",
       width = 9, height = 10, units = "in")


stability3 <- bind_rows(typeII, MSEs)
stability3$stability_measure <- factor(stability3$stability_measure,
                                       levels = c("Type II", "MSE"))

midwest_vs_non_all3 <-
  ggplot(stability3, aes(x = trait, y = value, fill = region)) +
  geom_boxplot() +
  scale_fill_manual(values = c("steelblue", "darkred")) +
  stat_compare_means(method = "t.test",
                     label.y = -0.5,
                     aes(label = paste0(..p.signif..), family = "Courier")) +
  labs(x = "Trait", y = "", fill = "Region") +
  facet_grid(stability_measure~., scales = "free_y") +
  theme_bw() +
  theme(
    # text = element_text(family = "Courier"),
    axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
    axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 10),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
    panel.border = element_blank(), axis.line = element_line(color = "black"),
    # legend.position = c(0.15,0.9),
    strip.text = element_text(size = 10), strip.background = element_rect(color = "white", fill = "white")
  )

ggsave(plot = midwest_vs_non_all3, filename = "midwest_vs_non_all3.tif", device = "tiff",
       width = 9, height = 6.5, units = "in")


