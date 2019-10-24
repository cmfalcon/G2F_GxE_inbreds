##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: GxE-GGE-Stability analysis
## Date: 2018-09-20
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)
library(reshape2)
library(ggplot2)
library(forcats)
library(ggpubr)

rowMedian <- function(x, na.rm = FALSE) apply(x, 1, median, na.rm = na.rm)
median_ <- function(...) median(..., na.rm = T)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

#### Read in, format, and organize data ####
stabilityRank_allTraits <- read.csv("stability_rank_all_traits.csv", header = T)

selectedPhenoData <- read.csv("./selectedPhenoData.csv", header = T)
colnames(selectedPhenoData) <- c("Plot", "Enviro", "Rep", "Genotype", "Anthesis (GDU)", "Silking (GDU)", 
                                 "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                 "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                 "Kernel length", "Kernel width", "Kernel thickness")


#### Stability ####
stabilityRank_allTraits2 <- stabilityRank_allTraits[,seq(1, ncol(stabilityRank_allTraits), 2)]
stabilityRank_allTraits2$medianRank <- as.numeric(rowMedian(stabilityRank_allTraits2[,2:ncol(stabilityRank_allTraits2)], na.rm = T))
stabilityRank_allTraits2[stabilityRank_allTraits2 == "NaN"] <- NA_character_
write.csv(stabilityRank_allTraits2, "stabilityRank_allTraits2.csv", row.names = F)

stabilityRank_allTraits2_melt <- melt(stabilityRank_allTraits2)

stability <- ggplot(stabilityRank_allTraits2_melt, mapping = aes(x = reorder(genos, -value, median, na.rm = T), y = value)) +
  geom_boxplot() +
  # geom_jitter() +
  labs(x = "Genotype", y = "Stability rank") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = c(0.15,0.9), legend.title = element_text(" "))
ggsave(plot = stability, filename = "medianStability.tif", device = "tiff",
       width = 5, height = 7.75, units = "in")


stabilityValue_allTraits <- stabilityRank_allTraits[,c(1, seq(2, ncol(stabilityRank_allTraits), 2))]
# stabilityValue_allTraits$median <- as.numeric(rowMedian(stabilityValue_allTraits[,2:ncol(stabilityValue_allTraits)], na.rm = T))
stabilityValue_allTraits[stabilityValue_allTraits == "NaN"] <- NA_character_
traits <- colnames(selectedPhenoData[5:ncol(selectedPhenoData)])
colnames(stabilityValue_allTraits) <- c("Genotype", traits)
stabilityValue_allTraits <- subset(stabilityValue_allTraits, Genotype != "PHG83")

write.csv(stabilityValue_allTraits, "stability_value_all_traits.csv", row.names = F)

stabilityValue_allTraits_withMed <- stabilityValue_allTraits
stabilityValue_allTraits_withMed$median <- as.numeric(rowMedian(stabilityValue_allTraits_withMed[,2:ncol(stabilityValue_allTraits_withMed)], na.rm = T))

write.csv(stabilityValue_allTraits_withMed, "stability_value_all_traits_with_median.csv", row.names = F)

stabilityValue_allTraits_melt <- melt(stabilityValue_allTraits)

stabilityValue <- ggplot(stabilityValue_allTraits_melt, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot(outlier.color = "white") +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  coord_flip() +
  labs(x = "Genotypes", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue, filename = "medianStabilityValue.tif", device = "tiff",
       width = 5, height = 7.75, units = "in")

stabilityValue2 <- ggplot(stabilityValue_allTraits_melt, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = variable, shape = variable)) +
  # scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  # scale_shape_manual(values = c(rep(0:6, 3))) +
  coord_flip() +
  labs(x = "Genotype", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue2, filename = "medianStabilityValue_noJitter.tif", device = "tiff",
       width = 5, height = 7.75, units = "in")


stabilityValue_allTraits_melt_traitType <- stabilityValue_allTraits_melt
stabilityValue_allTraits_melt_traitType$traitType <- NA
for (i in 1:nrow(stabilityValue_allTraits_melt_traitType)) {
  if (stabilityValue_allTraits_melt_traitType[i, "variable"] == "Anthesis (GDU)" |
      stabilityValue_allTraits_melt_traitType[i, "variable"] == "Silking (GDU)" |
      stabilityValue_allTraits_melt_traitType[i, "variable"] == "Plant height" |
      stabilityValue_allTraits_melt_traitType[i, "variable"] == "Ear height") {
    stabilityValue_allTraits_melt_traitType[i, "traitType"] <- "Non"
  } else {
    stabilityValue_allTraits_melt_traitType[i, "traitType"] <- "Yield-component"
  }
}

stabilityValue_allTraits_melt_traitType$traitType2 <- NA
for (i in 1:nrow(stabilityValue_allTraits_melt_traitType)) {
  if (stabilityValue_allTraits_melt_traitType[i, "variable"] == "Anthesis (GDU)" |
      stabilityValue_allTraits_melt_traitType[i, "variable"] == "Silking (GDU)") {
    stabilityValue_allTraits_melt_traitType[i, "traitType2"] <- "Flowering date"
  } else if (stabilityValue_allTraits_melt_traitType[i, "variable"] == "Plant height" |
             stabilityValue_allTraits_melt_traitType[i, "variable"] == "Ear height") {
    stabilityValue_allTraits_melt_traitType[i, "traitType2"] <- "Height"
  } else {
    stabilityValue_allTraits_melt_traitType[i, "traitType2"] <- "Yield-component"
  }
}

stabilityValue_traitType <- ggplot(stabilityValue_allTraits_melt_traitType, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = traitType, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:14, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_traitType, filename = "medianStabilityValue_traitType.tif", device = "tiff")

stabilityValue_traitType2 <- ggplot(stabilityValue_allTraits_melt_traitType, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = traitType2, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:14, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_traitType2, filename = "medianStabilityValue_traitType2.tif", device = "tiff")

stabilityValue_melt_yieldComponents <- subset(stabilityValue_allTraits_melt_traitType, stabilityValue_allTraits_melt_traitType$traitType == "Yield-component")
stabilityValue_melt_nonYieldComponents <- subset(stabilityValue_allTraits_melt_traitType, stabilityValue_allTraits_melt_traitType$traitType == "Non")

stabilityValue_yieldComponents <- ggplot(stabilityValue_melt_yieldComponents, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_yieldComponents, filename = "medianStabilityValue_yieldComponents.tif", device = "tiff")

stabilityValue_nonYieldComponents <- ggplot(stabilityValue_melt_nonYieldComponents, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "darkgreen", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_nonYieldComponents, filename = "medianstabilityValue_nonYieldComponents.tif", device = "tiff")

stabilityValue_melt_flowering <- subset(stabilityValue_allTraits_melt_traitType, stabilityValue_allTraits_melt_traitType$traitType2 == "Flowering date")
stabilityValue_melt_height <- subset(stabilityValue_allTraits_melt_traitType, stabilityValue_allTraits_melt_traitType$traitType2 == "Height")

stabilityValue_flowering <- ggplot(stabilityValue_melt_flowering, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_flowering, filename = "medianStabilityValue_flowering.tif", device = "tiff")

stabilityValue_height <- ggplot(stabilityValue_melt_height, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("darkgreen", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_height, filename = "medianStabilityValue_height.tif", device = "tiff")

genos_ordered <- reorder(x = stabilityValue_melt_yieldComponents$Genotype, X = -stabilityValue_melt_yieldComponents$value, FUN = median_)
genos_ordered <- levels(genos_ordered)
genos_ordered <- genos_ordered[1:31]
genos_ordered2 <- c(genos_ordered, genos_ordered)

stabilityValue_flowering <- ggplot(stabilityValue_melt_flowering, mapping = aes(x = fct_inorder(genos_ordered2), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  scale_x_discrete(limits = genos_ordered) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "stabilityinabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_flowering, filename = "medianStabilityValue_flowering_orderedYC.tif", device = "tiff")

stabilityValue_height <- ggplot(stabilityValue_melt_height, mapping = aes(x = fct_inorder(genos_ordered2), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("darkgreen", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "stabilityinabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_height, filename = "medianStabilityValue_height_orderedYC.tif", device = "tiff")

genos_ordered3 <- c(genos_ordered, genos_ordered, genos_ordered, genos_ordered, 
                    genos_ordered, genos_ordered, genos_ordered, genos_ordered,
                    genos_ordered, genos_ordered, genos_ordered, genos_ordered,
                    genos_ordered, genos_ordered, genos_ordered)

stabilityValue_allTraits_melt_traitType$value2 <- NA
for (i in 1:nrow(stabilityValue_allTraits_melt_traitType)) {
  if (stabilityValue_allTraits_melt_traitType[i,"traitType2"] == "Yield-component") {
    stabilityValue_allTraits_melt_traitType[i,"value2"] <- stabilityValue_allTraits_melt_traitType[i,"value"]
  }
}

stabilityValue_panels <- 
  ggplot(stabilityValue_allTraits_melt_traitType, mapping = aes(x = reorder(Genotype, value2, median_), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = variable, shape = variable)) +
  # scale_color_manual(values = c("red", "goldenrod", rep(c("darkgreen", "blue"), 20))) +
  # scale_shape_manual(values = c(rep(0:6, 20))) +
  coord_flip() +
  facet_grid(~ traitType2) +
  labs(x = "Genotype", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16), strip.background = element_rect(color = "white", fill = "white"),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(plot = stabilityValue_panels, filename = "medianStabilityValue_orderedYC_panels.tif", device = "tiff",
       width = 7.5, height = 6.5, units = "in")

stabilityValue_panels_eachTrait <- 
  ggplot(stabilityValue_allTraits_melt_traitType, mapping = aes(x = reorder(Genotype, value2, median_), y = value)) +
  geom_col() +
  coord_flip() +
  facet_grid(~ variable) +
  labs(x = "Genotype", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12, angle = 90), strip.background = element_rect(color = "white", fill = "white"),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(plot = stabilityValue_panels_eachTrait, filename = "medianStabilityValue_panelForEachTrait.tif", device = "tiff",
       width = 9, height = 6.5, units = "in")

stabilityValue3 <- ggplot(stabilityValue_allTraits_melt, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.title = element_text(" "))
ggsave(plot = stabilityValue3, filename = "medianStabilityValue3.tif", device = "tiff",
       width = 5, height = 7.75, units = "in")


combinedStabilityBoxplots <- ggarrange(stabilityValue3, stabilityValue_panels,
                                       widths = c(1,1.95),
                                       labels = c("A", "B"),
                                       ncol = 2, nrow = 1)

ggsave(combinedStabilityBoxplots, filename = "combined_stability_boxplots.tif", device = "tiff",
       height = 8.5, width = 11, units = "in")

## stability by release year
lineMetadata <- read.csv("gxeInbredLines-Metadata.csv", header = T)
lineMetadata$Year <- as.factor(lineMetadata$Year)
colnames(stabilityRank_allTraits2)[1] <- "Genotype"
stabilityRank_andMetadata <- merge(stabilityRank_allTraits2, lineMetadata, by = "Genotype")

stabilityRank_andMetadata_melt <- melt(stabilityRank_andMetadata)
stabilityRank_andMetadata_melt$Year <- as.numeric.factor(stabilityRank_andMetadata_melt$Year)

stability2 <- ggplot(stabilityRank_andMetadata_melt, mapping = aes(x = reorder(stabilityRank_andMetadata_melt$Genotype, stabilityRank_andMetadata_melt$Year), y = value)) +
  geom_boxplot() +
  # geom_jitter() +
  labs(x = "Genotype", y = "Stability rank") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = c(0.15,0.9), legend.title = element_text(" "))
ggsave(plot = stability2, filename = "medianStability2.tif", device = "tiff")


