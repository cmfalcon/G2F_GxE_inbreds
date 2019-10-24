##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: GxE-Compare types of stability
## Date: 2018-09-20
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)
library(ggplot2)
library(ggcorrplot)
library(corrplot)
library(agricolae)
library(data.table)
library(ggpubr)
library(dplyr)
library(broman)

median_ <- function(...) median(..., na.rm = T)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

#### Read in, format, and organize data ####
allSlopesWithMedian <- read.csv("all_slopes_with_median.csv", header = T)
colnames(allSlopesWithMedian) <- c("Genotype", "Anthesis", "Silking", 
                                        "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                        "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                        "Kernel length", "Kernel width", "Kernel thickness", "Median")

allMSEWithMedian <- read.csv("all_MSE_with_median.csv", header = T)
colnames(allMSEWithMedian) <- c("Genotype", "Anthesis", "Silking", 
                                        "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                        "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                        "Kernel length", "Kernel width", "Kernel thickness", "Median")

stabilityValue_allTraits <- read.csv("stability_value_all_traits.csv", header = T)
colnames(stabilityValue_allTraits) <- c("Genotype", "Anthesis", "Silking", 
                        "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                        "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                        "Kernel length", "Kernel width", "Kernel thickness")
stabilityValue_allTraits_melt <- melt(stabilityValue_allTraits)

stabilityValue_allTraits_withMed <- read.csv("stability_value_all_traits_with_median.csv", header = T)
colnames(stabilityValue_allTraits_withMed) <- c("Genotype", "Anthesis", "Silking", 
                                        "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                        "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                        "Kernel length", "Kernel width", "Kernel thickness", "Median")

allBlups <- as.data.frame(read.csv("gxe20142015_allTraits_BLUPs.csv", header = T))
colnames(allBlups) <- c("Genotype", "Anthesis", "Silking", 
                                 "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                 "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                 "Kernel length", "Kernel width", "Kernel thickness")



#### Rank correlation of stability of different types of traits/different traits ####
data_stability <- stabilityValue_allTraits[,-1]
cormat <- cor(data_stability, use = "p", method = "pearson")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(data_stability)
head(p.mat[, 1:5])

tiff("stability-correlationPlot.tif")
corrplot(cormat, method = "color", 
         # addCoef.col = "gray40",
         tl.col = "black",  
         type = "upper", diag = F, p.mat = p.mat, 
         mar = c(0,3,0,0), insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white")
dev.off()

ggcorrplot_spearman <- ggcorrplot(cor(data_stability, use = "complete.obs", method = "spearman"), p.mat = p.mat, type = "upper",
                                  colors = c("darkred", "white", "steelblue"),
                                  ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_spearman, filename = "stability-ggCorrPlot-spearman.tif", device = "tiff", width = 6, height = 4.5, units = "in")

ggcorrplot_pearson <- ggcorrplot(cor(data_stability, use = "complete.obs", method = "pearson"), p.mat = p.mat, type = "upper",
                                 colors = c("darkred", "white", "steelblue"),
                                 ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_pearson, filename = "stability-ggCorrPlot-pearson.tif", device = "tiff", width = 6, height = 4.5, units = "in")


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

combinedStabilityBoxplotAndCorrplot <- ggarrange(stabilityValue3, ggcorrplot_spearman,
                                                 widths = c(1,1.8),
                                                 labels = c("A", "B"),
                                                 ncol = 2, nrow = 1)

ggsave(combinedStabilityBoxplotAndCorrplot, filename = "combined_stability_boxplots_and_correlation.tif", device = "tiff",
       height = 8, width = 11, units = "in")

#### Combine BLUPs with line metadata and do boxplots by category ####
lineMetadata <- read.csv("gxeInbredLines-Metadata.csv", header = T)
allMetadataAndBlups <- merge(lineMetadata, allBlups, by = "Genotype")
allMetadataAndBlups$X <- NULL; allMetadataAndBlups$X.1 <- NULL; allMetadataAndBlups$X.2 <- NULL; allMetadataAndBlups$X.3 <- NULL; allMetadataAndBlups$X.4 <- NULL; allMetadataAndBlups$X.5 <- NULL
rm(lineMetadata)

boxplot_theme <- theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())

## Boxplots Silking
boxplot_silking_allInbreds <- ggplot(allMetadataAndBlups, mapping = aes(x = "", y = Silking)) +
  geom_boxplot() +
  labs(x = "All inbreds") +
  coord_cartesian(ylim = c(1275, 1550)) +
  boxplot_theme

boxplot_silking_group <- ggplot(allMetadataAndBlups, mapping = aes(x = Group, y = Silking)) +
  geom_boxplot() +
  labs(y = "") +
  boxplot_theme +
  coord_cartesian(ylim = c(1275, 1550)) +
  theme(axis.text.y = element_blank())

boxplot_silking_exPVP <- ggplot(allMetadataAndBlups, mapping = aes(x = factor(Ex.PVP,levels = c("Non", "Ex-PVP")), y = Silking)) +
  geom_boxplot() +
  labs(x = "Ex-PVP", y = "") +
  boxplot_theme +
  coord_cartesian(ylim = c(1275, 1550)) +
  theme(axis.text.y = element_blank())

boxplot_silking_combined <- ggarrange(boxplot_silking_allInbreds, boxplot_silking_group, boxplot_silking_exPVP, 
                                      nrow = 1, ncol = 3)

ggsave(boxplot_silking_combined, filename = "boxplot-silking-by_categories.tif", device = "tiff",
       height = 3.5, width = 4.5, units = "in")

## Boxplots Plant height
boxplot_plant_height_allInbreds <- ggplot(allMetadataAndBlups, mapping = aes(x = "", y = `Plant height`)) +
  geom_boxplot() +
  labs(x = "All inbreds") +
  coord_cartesian(ylim = c(120, 191)) +
  boxplot_theme

boxplot_plant_height_group <- ggplot(allMetadataAndBlups, mapping = aes(x = Group, y = `Plant height`)) +
  geom_boxplot() +
  labs(y = "") +
  coord_cartesian(ylim = c(120, 191)) +
  boxplot_theme +
  theme(axis.text.y = element_blank())

boxplot_plant_height_exPVP <- ggplot(allMetadataAndBlups, mapping = aes(x = factor(Ex.PVP,levels = c("Non", "Ex-PVP")), y = `Plant height`)) +
  geom_boxplot() +
  labs(x = "Ex-PVP", y = "") +
  coord_cartesian(ylim = c(120, 191)) +
  boxplot_theme +
  theme(axis.text.y = element_blank())


boxplot_plant_height_combined <- ggarrange(boxplot_plant_height_allInbreds, boxplot_plant_height_group, boxplot_plant_height_exPVP, 
                                           nrow = 1, ncol = 3)

ggsave(boxplot_plant_height_combined, filename = "boxplot-plant_height-by_categories.tif", device = "tiff",
       height = 3.5, width = 4.5, units = "in")

## Boxplots Kernel row number
boxplot_kernel_row_number_allInbreds <- ggplot(allMetadataAndBlups, mapping = aes(x = "", y = `Kernel row number`)) +
  geom_boxplot() +
  labs(x = "All inbreds") +
  coord_cartesian(ylim = c(10.5, 17.4)) +
  boxplot_theme

boxplot_kernel_row_number_group <- ggplot(allMetadataAndBlups, mapping = aes(x = Group, y = `Kernel row number`)) +
  geom_boxplot() +
  labs(y = "") +
  coord_cartesian(ylim = c(10.5, 17.4)) +
  boxplot_theme +
  theme(axis.text.y = element_blank())

boxplot_kernel_row_number_exPVP <- ggplot(allMetadataAndBlups, mapping = aes(x = factor(Ex.PVP,levels = c("Non", "Ex-PVP")), y = `Kernel row number`)) +
  geom_boxplot() +
  labs(x = "Ex-PVP", y = "") +
  coord_cartesian(ylim = c(10.5, 17.4)) +
  boxplot_theme +
  theme(axis.text.y = element_blank()) 


boxplot_kernel_row_number_combined <- ggarrange(boxplot_kernel_row_number_allInbreds, boxplot_kernel_row_number_group, boxplot_kernel_row_number_exPVP, 
                                                nrow = 1, ncol = 3)

ggsave(boxplot_kernel_row_number_combined, filename = "boxplot-kernel_row_number-by_categories.tif", device = "tiff",
       height = 3.5, width = 4.5, units = "in")

## Boxplots Kernel weight
boxplot_kernel_weight_allInbreds <- ggplot(allMetadataAndBlups, mapping = aes(x = "", y = `Kernel weight`)) +
  geom_boxplot() +
  labs(x = "All inbreds") +
  coord_cartesian(ylim = c(0.185, 0.31)) +
  boxplot_theme

boxplot_kernel_weight_group <- ggplot(allMetadataAndBlups, mapping = aes(x = Group, y = `Kernel weight`)) +
  geom_boxplot() +
  labs(y = "") +
  boxplot_theme +
  coord_cartesian(ylim = c(0.185, 0.31)) +
  theme(axis.text.y = element_blank())

boxplot_kernel_weight_exPVP <- ggplot(allMetadataAndBlups, mapping = aes(x = factor(Ex.PVP,levels = c("Non", "Ex-PVP")), y = `Kernel weight`)) +
  geom_boxplot() +
  labs(x = "Ex-PVP", y = "") +
  boxplot_theme +
  coord_cartesian(ylim = c(0.185, 0.31)) +
  theme(axis.text.y = element_blank())

boxplot_kernel_weight_combined <- ggarrange(boxplot_kernel_weight_allInbreds, boxplot_kernel_weight_group, boxplot_kernel_weight_exPVP, 
                                            nrow = 1, ncol = 3)

ggsave(boxplot_kernel_weight_combined, filename = "boxplot-kernel_weight-by_categories.tif", device = "tiff",
       height = 3.5, width = 4.5, units = "in")

boxplot_all_traits_combined <- ggarrange(boxplot_silking_combined, boxplot_plant_height_combined, boxplot_kernel_row_number_combined, boxplot_kernel_weight_combined,
                                         nrow = 4, ncol = 1,
                                         labels = c("A", "B", "C", "D"))
ggsave(boxplot_all_traits_combined, filename = "boxplot-all_traits-combined.tif", device = "tiff",
       height = 12, width = 5, units = "in")


boxplot_all_traits_combined_no_categories <- ggarrange(boxplot_silking_allInbreds, boxplot_plant_height_allInbreds, boxplot_kernel_row_number_allInbreds, boxplot_kernel_weight_allInbreds,
                                                       nrow = 2, ncol = 2,
                                                       labels = c("A", "B", "C", "D"))
ggsave(boxplot_all_traits_combined_no_categories, filename = "boxplot-all_traits-combined_no_categories.tif", device = "tiff",
       height = 5, width = 4, units = "in")


#### Combine slope, MSE, and GGE stability with line metadata ####
lineMetadata <- read.csv("gxeInbredLines-Metadata.csv", header = T)

slopeAndMetadata <- merge(lineMetadata, allSlopesWithMedian, by = "Genotype")

MSEandMetadata <- merge(lineMetadata, allMSEWithMedian, by = "Genotype")

GGEstabilityAndMetadata <- merge(lineMetadata, stabilityValue_allTraits, by = "Genotype")

#### Calculate Type II stability ####
typeIIandMetadata <- slopeAndMetadata %>%
  mutate(Anthesis_type2 = abs(Anthesis - 1)) %>%
  mutate(Silking_type2 = abs(Silking - 1)) %>%
  mutate(PlantHeight_type2 = abs(`Plant height` - 1)) %>%
  mutate(EarHeight_type2 = abs(`Ear height` - 1)) %>%
  mutate(PlotWeight_type2 = abs(`Plot grain weight` - 1)) %>%
  mutate(EarLength_type2 = abs(`Ear length` - 1)) %>%
  mutate(EarWidth_type2 = abs(`Ear width` - 1)) %>%
  mutate(KernelsPerRow_type2 = abs(`Kernels per row` - 1)) %>%
  mutate(KRN_type2 = abs(`Kernel row number` - 1)) %>%
  mutate(KernelWeight_type2 = abs(`Kernel weight` - 1)) %>%
  mutate(KernelArea_type2 = abs(`Kernel area` - 1)) %>%
  mutate(KernelLength_type2 = abs(`Kernel length` - 1)) %>%
  mutate(KernelWidth_type2 = abs(`Kernel width` - 1)) %>%
  mutate(KernelThickness_type2 = abs(`Kernel thickness` - 1)) %>%
  mutate(median_type2 = abs(Median - 1)) %>%
  select(Genotype, Pedigree, Derivation, Ex.PVP, Origin, Group, Company, Year, Release, Accession_Number.Origin,
         Anthesis_type2, Silking_type2, PlantHeight_type2, EarHeight_type2, PlotWeight_type2,
         EarLength_type2, EarWidth_type2, KernelsPerRow_type2, KRN_type2, KernelWeight_type2,
         KernelArea_type2, KernelLength_type2, KernelWidth_type2, KernelThickness_type2, median_type2)

#### Compare NSS vs SSS ####
# dir.create("./NSSvsSSS")
# dir.create("./NSSvsSSS/boxplotSlope")
# dir.create("./NSSvsSSS/boxplotTypeII")
# dir.create("./NSSvsSSS/boxplotMSE")
# dir.create("./NSSvsSSS/boxplotGGEstability")

#### Slope
sigDiffSlope_NSSvsSSS <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffSlope_NSSvsSSS) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "Slope", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./NSSvsSSS/boxplotSlope/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(slopeAndMetadata) - 1)) {
  for (j in 6) {
    # i = 20
    # j = 6
    
    subpopMeans <- aggregate(slopeAndMetadata[,i],by = list(slopeAndMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(slopeAndMetadata$Genotype,by = list(slopeAndMetadata[,j]),FUN = length)
    
    NSS <- slopeAndMetadata[which(slopeAndMetadata$Group == "NSS"),]
    SSS <- slopeAndMetadata[which(slopeAndMetadata$Group == "SSS"),]
    tTest <- t.test(NSS[,i], SSS[,i])
    pval <- paste0("p = ", myround(tTest$p.value,2))
    yLabPos <- max(slopeAndMetadata[,i]) + (0.00000001 * mean(slopeAndMetadata[,i]))
    
    myBoxplot(slopeAndMetadata, slopeAndMetadata[,j], slopeAndMetadata[,i])
    
    sigDiffSlope_NSSvsSSS <- rbind(sigDiffSlope_NSSvsSSS, myround(tTest$p.value,2))
  }
}

rownames(sigDiffSlope_NSSvsSSS) <- colnames(slopeAndMetadata[,11:(ncol(slopeAndMetadata) - 1)])


#### Type II Stability (abs(slope-1))
sigDiffTypeII_NSSvsSSS <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffTypeII_NSSvsSSS) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "|slope-1|", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./NSSvsSSS/boxplotTypeII/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(typeIIandMetadata) - 1)) {
  for (j in 6) {
    # i = 20
    # j = 6
    
    subpopMeans <- aggregate(typeIIandMetadata[,i],by = list(typeIIandMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(typeIIandMetadata$Genotype,by = list(typeIIandMetadata[,j]),FUN = length)
    
    NSS <- typeIIandMetadata[which(typeIIandMetadata$Group == "NSS"),]
    SSS <- typeIIandMetadata[which(typeIIandMetadata$Group == "SSS"),]
    tTest <- t.test(NSS[,i], SSS[,i])
    pval <- paste0("p = ", myround(tTest$p.value,2))
    yLabPos <- max(typeIIandMetadata[,i]) + (0.00000001 * mean(typeIIandMetadata[,i]))
    
    myBoxplot(typeIIandMetadata, typeIIandMetadata[,j], typeIIandMetadata[,i])
    
    sigDiffTypeII_NSSvsSSS <- rbind(sigDiffTypeII_NSSvsSSS, myround(tTest$p.value,2))
  }
}

rownames(sigDiffTypeII_NSSvsSSS) <- colnames(typeIIandMetadata[,11:(ncol(typeIIandMetadata) - 1)])


#### MSE
sigDiffMSE_NSSvsSSS <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffMSE_NSSvsSSS) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "MSE", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./NSSvsSSS/boxplotMSE/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(MSEandMetadata) - 1)) {
  for (j in 6) {
    # i = 12
    # j = 6
    subpopMeans <- aggregate(MSEandMetadata[,i],by = list(MSEandMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(MSEandMetadata$Genotype,by = list(MSEandMetadata[,j]),FUN = length)
    
    NSS <- MSEandMetadata[which(MSEandMetadata$Group == "NSS"),]
    SSS <- MSEandMetadata[which(MSEandMetadata$Group == "SSS"),]
    tTest <- t.test(NSS[,i], SSS[,i])
    pval <- paste0("p = ", myround(tTest$p.value, 2))
    yLabPos <- max(MSEandMetadata[,i]) + (0.00000001 * mean(MSEandMetadata[,i]))
    
    myBoxplot(MSEandMetadata, MSEandMetadata[,j], MSEandMetadata[,i])
    
    sigDiffMSE_NSSvsSSS <- rbind(sigDiffMSE_NSSvsSSS, myround(tTest$p.value,2))
  }
}

rownames(sigDiffMSE_NSSvsSSS) <- colnames(MSEandMetadata[,11:(ncol(MSEandMetadata) - 1)])

#### GGE Stability
sigDiffGGEstability_NSSvsSSS <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffGGEstability_NSSvsSSS) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "GGE Stability Value", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./NSSvsSSS/boxplotGGEstability/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:ncol(GGEstabilityAndMetadata)) {
  for (j in 6) {
    # i = 12
    # j = 6
    subpopMeans <- aggregate(GGEstabilityAndMetadata[,i],by = list(GGEstabilityAndMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(GGEstabilityAndMetadata$Genotype,by = list(GGEstabilityAndMetadata[,j]),FUN = length)
    
    NSS <- GGEstabilityAndMetadata[which(GGEstabilityAndMetadata$Group == "NSS"),]
    SSS <- GGEstabilityAndMetadata[which(GGEstabilityAndMetadata$Group == "SSS"),]
    tTest <- t.test(NSS[,i], SSS[,i])
    pval <- paste0("p = ", myround(tTest$p.value, 2))
    yLabPos <- max(GGEstabilityAndMetadata[,i]) + (0.00000001 * mean(GGEstabilityAndMetadata[,i]))
    
    myBoxplot(GGEstabilityAndMetadata, GGEstabilityAndMetadata[,j], GGEstabilityAndMetadata[,i])
    
    sigDiffGGEstability_NSSvsSSS <- rbind(sigDiffGGEstability_NSSvsSSS, myround(tTest$p.value,2))
  }
}

rownames(sigDiffGGEstability_NSSvsSSS) <- colnames(GGEstabilityAndMetadata[,11:ncol(GGEstabilityAndMetadata)])

#### Compare exPVP vs Non ####
# dir.create("./exPVPvsNon")
# dir.create("./exPVPvsNon/boxplotSlope")
# dir.create("./exPVPvsNon/boxplotTypeII")
# dir.create("./exPVPvsNon/boxplotMSE")
# dir.create("./exPVPvsNon/boxplotGGEstability")

#### Slope
sigDiffSlope_exPVPvsNon <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffSlope_exPVPvsNon) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "Slope", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./exPVPvsNon/boxplotSlope/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(slopeAndMetadata) - 1)) {
  for (j in 4) {
    # i = 11
    # j = 4
    
    subpopMeans <- aggregate(slopeAndMetadata[,i],by = list(slopeAndMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(slopeAndMetadata$Genotype,by = list(slopeAndMetadata[,j]),FUN = length)
    
    exPVP <- slopeAndMetadata[which(slopeAndMetadata$Ex.PVP == "Ex-PVP"),]
    Non <- slopeAndMetadata[which(slopeAndMetadata$Ex.PVP == "Non"),]
    tTest <- t.test(exPVP[,i], Non[,i])
    pval <- paste0("p = ", myround(tTest$p.value,2))
    yLabPos <- max(slopeAndMetadata[,i]) + (0.00000001 * mean(slopeAndMetadata[,i]))
    
    myBoxplot(slopeAndMetadata, slopeAndMetadata[,j], slopeAndMetadata[,i])
    
    sigDiffSlope_exPVPvsNon <- rbind(sigDiffSlope_exPVPvsNon, myround(tTest$p.value,2))
  }
}

rownames(sigDiffSlope_exPVPvsNon) <- colnames(slopeAndMetadata[,11:(ncol(slopeAndMetadata) - 1)])

#### Type II Stability (abs(slope-1))
sigDiffTypeII_exPVPvsNon <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffTypeII_exPVPvsNon) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "|slope-1|", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./exPVPvsNon/boxplotTypeII/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(typeIIandMetadata) - 1)) {
  for (j in 4) {
    # i = 11
    # j = 4
    
    subpopMeans <- aggregate(typeIIandMetadata[,i],by = list(typeIIandMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(typeIIandMetadata$Genotype,by = list(typeIIandMetadata[,j]),FUN = length)
    
    exPVP <- typeIIandMetadata[which(typeIIandMetadata$Ex.PVP == "Ex-PVP"),]
    Non <- typeIIandMetadata[which(typeIIandMetadata$Ex.PVP == "Non"),]
    tTest <- t.test(exPVP[,i], Non[,i])
    pval <- paste0("p = ", myround(tTest$p.value,2))
    yLabPos <- max(typeIIandMetadata[,i]) + (0.00000001 * mean(typeIIandMetadata[,i]))
    
    myBoxplot(typeIIandMetadata, typeIIandMetadata[,j], typeIIandMetadata[,i])
    
    sigDiffTypeII_exPVPvsNon <- rbind(sigDiffTypeII_exPVPvsNon, myround(tTest$p.value,2))
  }
}
rownames(sigDiffTypeII_exPVPvsNon) <- colnames(typeIIandMetadata[,11:(ncol(typeIIandMetadata) - 1)])


#### MSE
sigDiffMSE_exPVPvsNon <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffMSE_exPVPvsNon) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "MSE", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./exPVPvsNon/boxplotMSE/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(MSEandMetadata) - 1)) {
  for (j in 4) {
    # i = 12
    # j = 4
    subpopMeans <- aggregate(MSEandMetadata[,i],by = list(MSEandMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(MSEandMetadata$Genotype,by = list(MSEandMetadata[,j]),FUN = length)
    
    exPVP <- MSEandMetadata[which(MSEandMetadata$Ex.PVP == "Ex-PVP"),]
    Non <- MSEandMetadata[which(MSEandMetadata$Ex.PVP == "Non"),]
    tTest <- t.test(exPVP[,i], Non[,i])
    pval <- paste0("p = ", myround(tTest$p.value, 2))
    yLabPos <- max(MSEandMetadata[,i]) + (0.00000001 * mean(MSEandMetadata[,i]))
    
    myBoxplot(MSEandMetadata, MSEandMetadata[,j], MSEandMetadata[,i])
    
    sigDiffMSE_exPVPvsNon <- rbind(sigDiffMSE_exPVPvsNon, myround(tTest$p.value,2))
  }
}

rownames(sigDiffMSE_exPVPvsNon) <- colnames(MSEandMetadata[,11:(ncol(MSEandMetadata) - 1)])

#### GGE Stability
sigDiffGGEstability_exPVPvsNon <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffGGEstability_exPVPvsNon) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "GGE Stability Value", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./exPVPvsNon/boxplotGGEstability/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:ncol(GGEstabilityAndMetadata)) {
  for (j in 4) {
    # i = 12
    # j = 4
    subpopMeans <- aggregate(GGEstabilityAndMetadata[,i],by = list(GGEstabilityAndMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(GGEstabilityAndMetadata$Genotype,by = list(GGEstabilityAndMetadata[,j]),FUN = length)
    
    exPVP <- GGEstabilityAndMetadata[which(GGEstabilityAndMetadata$Ex.PVP == "Ex-PVP"),]
    Non <- GGEstabilityAndMetadata[which(GGEstabilityAndMetadata$Ex.PVP == "Non"),]
    tTest <- t.test(exPVP[,i], Non[,i])
    pval <- paste0("p = ", myround(tTest$p.value, 2))
    yLabPos <- max(GGEstabilityAndMetadata[,i]) + (0.00000001 * mean(GGEstabilityAndMetadata[,i]))
    
    myBoxplot(GGEstabilityAndMetadata, GGEstabilityAndMetadata[,j], GGEstabilityAndMetadata[,i])
    
    sigDiffGGEstability_exPVPvsNon <- rbind(sigDiffGGEstability_exPVPvsNon, myround(tTest$p.value,2))
  }
}

rownames(sigDiffGGEstability_exPVPvsNon) <- colnames(GGEstabilityAndMetadata[,11:ncol(GGEstabilityAndMetadata)])


sigDiffSlope_NSSvsSSS <- as.data.frame(sigDiffSlope_NSSvsSSS)
sigDiffTypeII_NSSvsSSS <- as.data.frame(sigDiffTypeII_NSSvsSSS)
sigDiffMSE_NSSvsSSS <- as.data.frame(sigDiffMSE_NSSvsSSS)
sigDiffGGEstability_NSSvsSSS <- as.data.frame(sigDiffGGEstability_NSSvsSSS)

sigDiffSlope_exPVPvsNon <- as.data.frame(sigDiffSlope_exPVPvsNon)
sigDiffTypeII_exPVPvsNon <- as.data.frame(sigDiffTypeII_exPVPvsNon)
sigDiffMSE_exPVPvsNon <- as.data.frame(sigDiffMSE_exPVPvsNon)
sigDiffGGEstability_exPVPvsNon <- as.data.frame(sigDiffGGEstability_exPVPvsNon)

sigDiff <- cbind(sigDiffSlope_NSSvsSSS, sigDiffTypeII_NSSvsSSS, sigDiffMSE_NSSvsSSS, sigDiffGGEstability_NSSvsSSS,
                 sigDiffSlope_exPVPvsNon, sigDiffTypeII_exPVPvsNon, sigDiffMSE_exPVPvsNon, sigDiffGGEstability_exPVPvsNon)

colnames(sigDiff) <- c("slope_NSSvsSSS", "typeII_NSSvsSSS", "MSE_NSSvsSSS", "GGE_NSSvsSSS",
                       "slope_exPVPvsNon", "typeII_exPVPvsNon", "MSE_exPVPvsNon", "GGE_exPVPvsNon")

sigDiff[,1] <- as.numeric.factor(sigDiff[,1])
sigDiff[,2] <- as.numeric.factor(sigDiff[,2])
sigDiff[,3] <- as.numeric.factor(sigDiff[,3])
sigDiff[,4] <- as.numeric.factor(sigDiff[,4])
sigDiff[,5] <- as.numeric.factor(sigDiff[,5])
sigDiff[,6] <- as.numeric.factor(sigDiff[,6])
sigDiff[,7] <- as.numeric.factor(sigDiff[,7])
sigDiff[,8] <- as.numeric.factor(sigDiff[,8])

write.csv(sigDiff, "sigDiffs_stability_genotypeCategories.csv")

#### Rank correlation of measures of stability ####
colnames(allSlopesWithMedian)[2:ncol(allSlopesWithMedian)] <- paste0("slope_", colnames(allSlopesWithMedian)[2:ncol(allSlopesWithMedian)])
colnames(allMSEWithMedian)[2:ncol(allMSEWithMedian)] <- paste0("MSE_", colnames(allMSEWithMedian)[2:ncol(allMSEWithMedian)])
colnames(stabilityValue_allTraits_withMed)[2:ncol(stabilityValue_allTraits_withMed)] <- paste0("GGE_", colnames(stabilityValue_allTraits_withMed)[2:ncol(stabilityValue_allTraits_withMed)])

typeIIStability <- typeIIandMetadata[,c(1,11:25)]

stabilityMeasures <- merge(allSlopesWithMedian, allMSEWithMedian, by = "Genotype")
stabilityMeasures <- merge(stabilityMeasures, stabilityValue_allTraits_withMed, by = "Genotype")
stabilityMeasures <- merge(stabilityMeasures, typeIIStability, by = "Genotype")

corr_slope_MSE <- cor(stabilityMeasures$slope_Median, stabilityMeasures$MSE_Median, method = "spearman")
corr_slope_GGE <- cor(stabilityMeasures$slope_Median, stabilityMeasures$GGE_Median, method = "spearman")

corr_typeII_MSE <- cor(stabilityMeasures$median_type2, stabilityMeasures$MSE_Median, method = "spearman")
corr_typeII_GGE <- cor(stabilityMeasures$median_type2, stabilityMeasures$GGE_Median, method = "spearman")

corr_MSE_GGE <- cor(stabilityMeasures$MSE_Median, stabilityMeasures$GGE_Median, method = "spearman")


cormat <- cor(stabilityMeasures[,-1], use = "p", method = "spearman")
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(stabilityMeasures[,-1])

p.mat2 <- p.mat %>%
  as.data.frame() %>%
  dplyr::select(MSE_Median,GGE_Median,median_type2)
  
ggcorrplot <- 
  ggcorrplot(cor(stabilityMeasures[,-1]), p.mat = cor.mtest(stabilityMeasures[,-1]), type = "upper",
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot, filename = "ggCorrPlot-stabilityMeasures.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")

allSlopes <- allSlopesWithMedian[,-16]
colnames(allSlopes) <- gsub(pattern = "slope_", replacement = "", x = colnames(allSlopes))
colnames(allSlopes) <- gsub(pattern = "_st", replacement = "", x = colnames(allSlopes))

ggcorrplot_slopeAmongTraits <- 
  ggcorrplot(cor(allSlopes[,-1]), p.mat = cor.mtest(allSlopes[,-1]), type = "upper",
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_slopeAmongTraits, filename = "ggCorrPlot-slopeAmongTraits.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")

typeIIStability <- typeIIStability[,-ncol(typeIIStability)]
colnames(typeIIStability) <- c("Genotype", "Anthesis", "Silking", "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                               "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                               "Kernel length", "Kernel width", "Kernel thickness")

ggcorrplot_typeIIAmongTraits <- 
  ggcorrplot(cor(typeIIStability[,-1]), p.mat = cor.mtest(typeIIStability[,-1]), type = "upper",
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_typeIIAmongTraits, filename = "ggCorrPlot-typeIIAmongTraits.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")

MSEstability <- allMSEWithMedian[,-16]
colnames(MSEstability) <- gsub(pattern = "MSE_", replacement = "", x = colnames(MSEstability))
colnames(MSEstability) <- gsub(pattern = "_st", replacement = "", x = colnames(MSEstability))

ggcorrplot_MSEAmongTraits <- 
  ggcorrplot(cor(MSEstability[,-1]), p.mat = cor.mtest(MSEstability[,-1]), type = "upper",
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_MSEAmongTraits, filename = "ggCorrPlot-MSEAmongTraits.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")

GGEstability <- stabilityValue_allTraits_withMed[-16]
colnames(GGEstability) <- gsub(pattern = "GGE_", replacement = "", x = colnames(GGEstability))
colnames(GGEstability) <- gsub(pattern = "\\(GDU)", replacement = "", x = colnames(GGEstability))

ggcorrplot_GGEAmongTraits <- 
  ggcorrplot(cor(GGEstability[,-1]), p.mat = cor.mtest(GGEstability[,-1]), type = "upper",
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_GGEAmongTraits, filename = "ggCorrPlot-GGEAmongTraits.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")


combinedStabilityCorrAmongTraits <- ggarrange(ggcorrplot_typeIIAmongTraits, ggcorrplot_MSEAmongTraits,
                                              nrow = 1, ncol = 2,
                                              labels = c("A", "B"))
ggsave(combinedStabilityCorrAmongTraits, filename = "combinedStabilityCorrAmongTraits.tif", device = "tiff",
       width = 9, height = 3.5, units = "in")


#### Plot of genotypes that are in 5 most or 5 least stable for slope, MSE, GGE ####
traits <- c("Anthesis", "Silking", 
            "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
            "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
            "Kernel length", "Kernel width", "Kernel thickness")

colnames(stabilityMeasures) <- c("Genotype", "slope_Anthesis", "slope_Silking", "slope_Plant height", "slope_Ear height", 
                                 "slope_Plot grain weight", "slope_Ear length", "slope_Ear width", "slope_Kernels per row", "slope_Kernel row number", 
                                 "slope_Kernel weight", "slope_Kernel area", "slope_Kernel length", "slope_Kernel width", "slope_Kernel thickness",
                                 "slope_medianSlope", "MSE_Anthesis", "MSE_Silking", "MSE_Plant height", "MSE_Ear height",     
                                 "MSE_Plot grain weight", "MSE_Ear length", "MSE_Ear width", "MSE_Kernels per row", "MSE_Kernel row number",           
                                 "MSE_Kernel weight", "MSE_Kernel area", "MSE_Kernel length", "MSE_Kernel width", "MSE_Kernel thickness",  
                                 "MSE_medianMSE", "GGE_Anthesis", "GGE_Silking", "GGE_Plant height", "GGE_Ear height",          
                                 "GGE_Plot grain weight", "GGE_Ear length", "GGE_Ear width", "GGE_Kernels per row", "GGE_Kernel row number",   
                                 "GGE_Kernel weight", "GGE_Kernel area", "GGE_Kernel length", "GGE_Kernel width", "GGE_Kernel thickness",    
                                 "GGE_median", "type2_Anthesis", "type2_Silking", "type2_Plant height", "type2_Ear height", 
                                 "type2_Plot grain weight", "type2_Ear length", "type2_Ear width", "type2_Kernels per row", "type2_Kernel row number", 
                                 "type2_Kernel weight", "type2_Kernel area", "type2_Kernel length", "type2_Kernel width", "type2_Kernel thickness",
                                 "type2_median")
stabilityMeasures$MSE_medianMSE <- NULL
stabilityMeasures$GGE_median <- NULL
stabilityMeasures$slope_medianSlope <- NULL
stabilityMeasures$type2_median <- NULL

stabilityMeasures_melt <- melt(stabilityMeasures, id.vars = "Genotype")
stabilityMeasures_melt$variable1 <- gsub(pattern = "_.*", replacement = "", stabilityMeasures_melt$variable)
stabilityMeasures_melt$variable2 <- gsub(pattern = "slope_", replacement = "", stabilityMeasures_melt$variable)
stabilityMeasures_melt$variable2 <- gsub(pattern = "MSE_", replacement = "", stabilityMeasures_melt$variable2)
stabilityMeasures_melt$variable2 <- gsub(pattern = "GGE_", replacement = "", stabilityMeasures_melt$variable2)
stabilityMeasures_melt$variable <- NULL
stabilityMeasures_melt <- stabilityMeasures_melt[,c(1,3,4,2)]
colnames(stabilityMeasures_melt) <- c("Genotype", "stability measure", "trait", "value")

stabilityMeasures_melt$`stability measure` <- as.factor(stabilityMeasures_melt$`stability measure`)
stabilityMeasures_melt$trait <- as.factor(stabilityMeasures_melt$trait)

####

least_and_most_stable_genos <- matrix(NA, nrow = 0, ncol = 4)
colnames(least_and_most_stable_genos) <- c("measure", "trait", "least", "most")

m = "GGE"
measureStability <- stabilityMeasures_melt %>%
  filter(`stability measure` == m)

for (t in 1:length(traits)) {
  # t = 1
  
  traitStability <- measureStability %>%
    filter(measureStability$trait == traits[t])
  
  # most stable = smallest distance from AEC
  mostStable <- traitStability %>%
    top_n(-5) %>%
    .$Genotype
  
  # least stable = greatest distance from AEC
  leastStable <- traitStability %>%
    top_n(5) %>%
    .$Genotype
  
  genos_rows <- c(rep(m, 5), rep(traits[t], 5), as.character(leastStable), as.character(mostStable))
  genos_matrix <- matrix(genos_rows, nrow = 5, ncol = 4)
  least_and_most_stable_genos <- rbind(least_and_most_stable_genos, genos_matrix)
}

m = "MSE"
measureStability <- stabilityMeasures_melt %>%
  filter(`stability measure` == m)

for (t in 1:length(traits)) {
  # t = 1
  
  traitStability <- measureStability %>%
    filter(measureStability$trait == traits[t])
  
  # most stable = smallest distance from AEC
  mostStable <- traitStability %>%
    top_n(-5) %>%
    .$Genotype
  
  # least stable = greatest distance from AEC
  leastStable <- traitStability %>%
    top_n(5) %>%
    .$Genotype
  
  genos_rows <- c(rep(m, 5), rep(traits[t], 5), as.character(leastStable), as.character(mostStable))
  genos_matrix <- matrix(genos_rows, nrow = 5, ncol = 4)
  least_and_most_stable_genos <- rbind(least_and_most_stable_genos, genos_matrix)
}

m = "slope"
measureStability <- stabilityMeasures_melt %>%
  filter(`stability measure` == m)
measureStability$value2 <- abs(measureStability$value - 1)

for (t in 1:length(traits)) {
  # t = 1
  
  traitStability <- measureStability %>%
    filter(measureStability$trait == traits[t])
  
  # most stable = smallest distance from AEC
  mostStable <- traitStability %>%
    top_n(-5) %>%
    .$Genotype
  
  # least stable = greatest distance from AEC
  leastStable <- traitStability %>%
    top_n(5) %>%
    .$Genotype
  
  genos_rows <- c(rep(m, 5), rep(traits[t], 5), as.character(leastStable), as.character(mostStable))
  genos_matrix <- matrix(genos_rows, nrow = 5, ncol = 4)
  least_and_most_stable_genos <- rbind(least_and_most_stable_genos, genos_matrix)
}

least_and_most_stable_genos <- as.data.frame(least_and_most_stable_genos)
least_and_most_stable_genos_melt <- melt(least_and_most_stable_genos, id.var = c("measure", "trait"))

write.csv(least_and_most_stable_genos, "least_and_most_stable_genos.csv", row.names = F)

myTable <- table(least_and_most_stable_genos_melt$measure, least_and_most_stable_genos_melt$variable, least_and_most_stable_genos_melt$value)
myTable_melt <- melt(myTable)

for (i in 1:nrow(myTable_melt)) {
  if (myTable_melt[i,"Var2"] == "least") {
    myTable_melt[i,"value"] <- -myTable_melt[i,"value"]
  }
}

colnames(myTable_melt) <- c("measure", "stability", "Genotype", "value")

stable_genos <-
  ggplot(myTable_melt, aes(x = Genotype, y = value)) +
  geom_bar(stat = "identity", fill = "gray50") +
  labs(x = "Genotype", y = "") +
  scale_y_continuous(breaks = c(-8, -4, 0, 4),
                     labels = c("8", "4", "0", "4")) +
  geom_abline(slope = 0, intercept = 0, size = 1.2) +
  annotate("text", label = "Most stable", x = 16, y = 8, size = 4, family = "Courier") +
  annotate("text", label = "Least stable", x = 16, y = -9, size = 4, family = "Courier") +
  facet_grid(measure~.) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = c(0.225,0.9), legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))
ggsave(plot = stable_genos, filename = paste("./most_least_stable.tif",sep = ""),
       device = "tiff", width = 5, height = 8, units = "in")


myTable_melt_with_metadata <- merge(lineMetadata, myTable_melt)

stable_genos_byYear <-
  ggplot(myTable_melt_with_metadata, aes(x = reorder(Genotype, Year), y = value)) +
  geom_bar(stat = "identity", fill = "gray50") +
  labs(x = "Genotype", y = "") +
  scale_y_continuous(breaks = c(-8, -4, 0, 4),
                     labels = c("8", "4", "0", "4")) +
  geom_abline(slope = 0, intercept = 0, size = 1.2) +
  annotate("text", label = "Most stable", x = 16, y = 8, size = 4, family = "Courier") +
  annotate("text", label = "Least stable", x = 16, y = -9, size = 4, family = "Courier") +
  facet_grid(measure~.) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = c(0.225,0.9), legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))
ggsave(plot = stable_genos_byYear, filename = paste("./most_least_stable_byYear.tif",sep = ""),
       device = "tiff", width = 5, height = 8, units = "in")

stable_genos_byExPVP <-
  ggplot(myTable_melt_with_metadata, aes(x = reorder(Genotype, Year), y = value, fill = Ex.PVP)) +
  geom_bar(stat = "identity") +
  labs(x = "Genotype", y = "") +
  scale_y_continuous(breaks = c(-8, -4, 0, 4),
                     labels = c("8", "4", "0", "4")) +
  geom_abline(slope = 0, intercept = 0, size = 1.2) +
  annotate("text", label = "Most stable", x = 16, y = 8, size = 4, family = "Courier") +
  annotate("text", label = "Least stable", x = 16, y = -9, size = 4, family = "Courier") +
  facet_grid(measure~.) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))
ggsave(plot = stable_genos_byExPVP, filename = paste("./most_least_stable_byExPVP.tif",sep = ""),
       device = "tiff", width = 6, height = 8, units = "in")


myTable_melt_with_metadata$Genotype <- factor(myTable_melt_with_metadata$Genotype,
                                              levels = unique(myTable_melt_with_metadata[order(myTable_melt_with_metadata$Group),"Genotype"]))


stable_genos_byGroup <-
  ggplot(myTable_melt_with_metadata, aes(x = reorder(Genotype, Group), y = value, fill = Group)) +
  geom_bar(stat = "identity") +
  labs(x = "Genotype", y = "") +
  scale_y_continuous(breaks = c(-8, -4, 0, 4),
                     labels = c("8", "4", "0", "4")) +
  geom_abline(slope = 0, intercept = 0, size = 1.2) +
  annotate("text", label = "Most stable", x = 16, y = 8, size = 4, family = "Courier") +
  annotate("text", label = "Least stable", x = 16, y = -9, size = 4, family = "Courier") +
  facet_grid(measure~.) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))
ggsave(plot = stable_genos_byGroup, filename = paste("./most_least_stable_byGroup.tif",sep = ""),
       device = "tiff", width = 6, height = 8, units = "in")


#### Correlation of performance (BLUPs) and stability ####
BLUPs_stabilityMeasures <- merge(allBlups, stabilityMeasures, by = "Genotype")

i = 2
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_anthesis <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_anthesis$typeII <- abs(BLUPs_stabilityMeasures_anthesis$slope_Anthesis - 1)

i = 3
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_silking <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_silking$typeII <- abs(BLUPs_stabilityMeasures_silking$slope_Silking - 1)

i = 4
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_plantHeight <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_plantHeight$typeII <- abs(BLUPs_stabilityMeasures_plantHeight$`slope_Plant height` - 1)

i = 5
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_earHeight <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_earHeight$typeII <- abs(BLUPs_stabilityMeasures_earHeight$`slope_Ear height` - 1)

i = 6
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_plotWeight <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_plotWeight$typeII <- abs(BLUPs_stabilityMeasures_plotWeight$`slope_Plot grain weight` - 1)

i = 7
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_earLength <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_earLength$typeII <- abs(BLUPs_stabilityMeasures_earLength$`slope_Ear length` - 1)

i = 8
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_earWidth <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_earWidth$typeII <- abs(BLUPs_stabilityMeasures_earWidth$`slope_Ear width` - 1)

i = 9
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelsPerRow <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelsPerRow$typeII <- abs(BLUPs_stabilityMeasures_kernelsPerRow$`slope_Kernels per row` - 1)

i = 10
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelRowNumber <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelRowNumber$typeII <- abs(BLUPs_stabilityMeasures_kernelRowNumber$`slope_Kernel row number` - 1)

i = 11
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelWeight <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelWeight$typeII <- abs(BLUPs_stabilityMeasures_kernelWeight$`slope_Kernel weight` - 1)

i = 12
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelArea <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelArea$typeII <- abs(BLUPs_stabilityMeasures_kernelArea$`slope_Kernel area` - 1)

i = 13
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelLength <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelLength$typeII <- abs(BLUPs_stabilityMeasures_kernelLength$`slope_Kernel length` - 1)

i = 14
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelWidth <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelWidth$typeII <- abs(BLUPs_stabilityMeasures_kernelWidth$`slope_Kernel width` - 1)

i = 15
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelThickness <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelThickness$typeII <- abs(BLUPs_stabilityMeasures_kernelThickness$`slope_Kernel thickness` - 1)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

ggcorrplot_anthesis_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_anthesis[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_anthesis[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_anthesis_BLUPs_stability, filename = "ggcorrplot_anthesis_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_silking_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_silking[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_silking[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_silking_BLUPs_stability, filename = "ggcorrplot_silking_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_plantHeight_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_plantHeight[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_plantHeight[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_plantHeight_BLUPs_stability, filename = "ggcorrplot_plantHeight_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_earHeight_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_earHeight[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_earHeight[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_earHeight_BLUPs_stability, filename = "ggcorrplot_earHeight_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_plotWeight_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_plotWeight[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_plotWeight[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_plotWeight_BLUPs_stability, filename = "ggcorrplot_plotWeight_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_earLength_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_earLength[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_earLength[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_earLength_BLUPs_stability, filename = "ggcorrplot_earLength_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_earWidth_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_earWidth[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_earWidth[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_earWidth_BLUPs_stability, filename = "ggcorrplot_earWidth_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_kernelsPerRow_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelsPerRow[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelsPerRow[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelsPerRow_BLUPs_stability, filename = "ggcorrplot_kernelsPerRow_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_kernelRowNumber_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelRowNumber[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelRowNumber[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelRowNumber_BLUPs_stability, filename = "ggcorrplot_kernelRowNumber_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_kernelWeight_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelWeight[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelWeight[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelWeight_BLUPs_stability, filename = "ggcorrplot_kernelWeight_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_kernelArea_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelArea[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelArea[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelArea_BLUPs_stability, filename = "ggcorrplot_kernelArea_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_kernelLength_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelLength[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelLength[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelLength_BLUPs_stability, filename = "ggcorrplot_kernelLength_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")


ggcorrplot_kernelWidth_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelWidth[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelWidth[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelWidth_BLUPs_stability, filename = "ggcorrplot_kernelWidth_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")


ggcorrplot_kernelThickness_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelThickness[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelThickness[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelThickness_BLUPs_stability, filename = "ggcorrplot_kernelThickness_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")




