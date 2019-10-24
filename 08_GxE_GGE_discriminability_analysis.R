##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: GxE-GGE-Discriminability analysis
## Date: 2018-09-20
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)
library(forcats)
library(ggplot2)
library(ggpubr)
library(broman)
library(reshape2)
library(ggcorrplot)
library(dplyr)
library(broom)

source("../ggThemes.R")

rowMedian <- function(x, na.rm = FALSE) apply(x, 1, median, na.rm = na.rm)
median_ <- function(...) median(..., na.rm = T)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

#### Read in, format, and organize data ####
discrimRank_allTraits <- read.csv("discrim_rank_all_traits.csv", header = T)

#### Discriminability ####
discrimRank_allTraits2 <- discrimRank_allTraits[,seq(1, ncol(discrimRank_allTraits), 2)]
discrimRank_allTraits2$medianRank <- as.numeric(rowMedian(discrimRank_allTraits2[,2:ncol(discrimRank_allTraits2)], na.rm = T))
discrimRank_allTraits2[discrimRank_allTraits2 == "NaN"] <- NA_character_
write.csv(discrimRank_allTraits2, "discrimRank_allTraits2_st.csv", row.names = F)

discrimRank_allTraits3 <- discrimRank_allTraits2[,1:(ncol(discrimRank_allTraits2) - 1)]
discrimRank_allTraits3_melt <- melt(discrimRank_allTraits3)
discrimRank_allTraits3_melt2 <- discrimRank_allTraits3_melt
discrimRank_allTraits3_melt2$site <- discrimRank_allTraits3_melt2$enviros
discrimRank_allTraits3_melt2$site <- gsub("_14", "", discrimRank_allTraits3_melt2$site)
discrimRank_allTraits3_melt2$site <- gsub("_15", "", discrimRank_allTraits3_melt2$site)

discrimValue_allTraits <- discrimRank_allTraits[,c(1, seq(2, ncol(discrimRank_allTraits), 2))]
discrimValue_allTraits[discrimValue_allTraits == "NaN"] <- NA_character_
traits <- c("Anthesis", "Silking", "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
            "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
            "Kernel length", "Kernel width", "Kernel thickness")
colnames(discrimValue_allTraits) <- c("enviros", traits)

discrimValue_allTraits_withMed <- discrimValue_allTraits
discrimValue_allTraits_withMed$median <- as.numeric(rowMedian(discrimValue_allTraits_withMed[,2:ncol(discrimValue_allTraits_withMed)], na.rm = T))

write.csv(discrimValue_allTraits_withMed, "discriminability_values.csv", row.names = F)

discrimValue_allTraits_melt <- melt(discrimValue_allTraits)

discrimValue <- ggplot(discrimValue_allTraits_melt, mapping = aes(x = reorder(x = enviros, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = discrimValue, filename = "medianDiscrimValue.tif", device = "tiff",
       width = 5, height = 7.75, units = "in")

discrimValue_melt_traitType <- discrimValue_allTraits_melt
discrimValue_melt_traitType$traitType <- NA
for (i in 1:nrow(discrimValue_melt_traitType)) {
  if (discrimValue_melt_traitType[i, "variable"] == "Anthesis" |
      discrimValue_melt_traitType[i, "variable"] == "Silking" |
      discrimValue_melt_traitType[i, "variable"] == "Plant height" |
      discrimValue_melt_traitType[i, "variable"] == "Ear height") {
    discrimValue_melt_traitType[i, "traitType"] <- "Non"
  } else {
    discrimValue_melt_traitType[i, "traitType"] <- "Yield-component"
  }
}

discrimValue_melt_traitType$traitType2 <- NA
for (i in 1:nrow(discrimValue_melt_traitType)) {
  if (discrimValue_melt_traitType[i, "variable"] == "Anthesis" |
      discrimValue_melt_traitType[i, "variable"] == "Silking") {
    discrimValue_melt_traitType[i, "traitType2"] <- "Flowering"
  } else if (discrimValue_melt_traitType[i, "variable"] == "Plant height" |
             discrimValue_melt_traitType[i, "variable"] == "Ear height") {
    discrimValue_melt_traitType[i, "traitType2"] <- "Height"
  } else {
    discrimValue_melt_traitType[i, "traitType2"] <- "Yield-component"
  }
}

discrimValue_traitType <- ggplot(discrimValue_melt_traitType, mapping = aes(x = reorder(x = enviros, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = traitType, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:14, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = discrimValue_traitType, filename = "medianDiscrimValue_traitType.tif", device = "tiff")

discrimValue_traitType2 <- ggplot(discrimValue_melt_traitType, mapping = aes(x = reorder(x = enviros, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = traitType2, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:14, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = discrimValue_traitType2, filename = "medianDiscrimValue_traitType2.tif", device = "tiff")

discrimValue_melt_yieldComponents <- subset(discrimValue_melt_traitType, discrimValue_melt_traitType$traitType == "Yield-component")
discrimValue_melt_nonYieldComponents <- subset(discrimValue_melt_traitType, discrimValue_melt_traitType$traitType == "Non")

discrimValue_yieldComponents <- ggplot(discrimValue_melt_yieldComponents, mapping = aes(x = reorder(x = enviros, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = discrimValue_yieldComponents, filename = "medianDiscrimValue_yieldComponents.tif", device = "tiff")

discrimValue_nonYieldComponents <- ggplot(discrimValue_melt_nonYieldComponents, mapping = aes(x = reorder(x = enviros, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "darkgreen", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = discrimValue_nonYieldComponents, filename = "medianDiscrimValue_nonYieldComponents.tif", device = "tiff")

discrimValue_melt_flowering <- subset(discrimValue_melt_traitType, discrimValue_melt_traitType$traitType2 == "Flowering")
discrimValue_melt_height <- subset(discrimValue_melt_traitType, discrimValue_melt_traitType$traitType2 == "Height")

discrimValue_flowering <- ggplot(discrimValue_melt_flowering, mapping = aes(x = reorder(x = enviros, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = discrimValue_flowering, filename = "medianDiscrimValue_flowering.tif", device = "tiff")

discrimValue_height <- ggplot(discrimValue_melt_height, mapping = aes(x = reorder(x = enviros, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("darkgreen", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = discrimValue_height, filename = "medianDiscrimValue_height.tif", device = "tiff")

enviros_ordered <- reorder(x = discrimValue_melt_yieldComponents$enviros, X = -discrimValue_melt_yieldComponents$value, FUN = median_)
enviros_ordered <- levels(enviros_ordered)
enviros_ordered2 <- c(enviros_ordered, enviros_ordered)

discrimValue_flowering <- ggplot(discrimValue_melt_flowering, mapping = aes(x = fct_inorder(enviros_ordered2), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  scale_x_discrete(limits = enviros_ordered) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = discrimValue_flowering, filename = "medianDiscrimValue_flowering_orderedYC.tif", device = "tiff")

discrimValue_height <- ggplot(discrimValue_melt_height, mapping = aes(x = fct_inorder(enviros_ordered2), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("darkgreen", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = discrimValue_height, filename = "medianDiscrimValue_height_orderedYC.tif", device = "tiff")

enviros_ordered3 <- c(enviros_ordered, enviros_ordered, enviros_ordered, enviros_ordered, 
                      enviros_ordered, enviros_ordered, enviros_ordered, enviros_ordered,
                      enviros_ordered, enviros_ordered, enviros_ordered, enviros_ordered,
                      enviros_ordered, enviros_ordered, enviros_ordered)

discrimValue_melt_traitType$value2 <- NA
for (i in 1:nrow(discrimValue_melt_traitType)) {
  if (discrimValue_melt_traitType[i,"traitType2"] == "Yield-component") {
    discrimValue_melt_traitType[i,"value2"] <- discrimValue_melt_traitType[i,"value"]
  }
}

discrimValue_panels <- 
  ggplot(discrimValue_melt_traitType, mapping = aes(x = reorder(enviros, value2, median_), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = variable, shape = variable)) +
  # scale_color_manual(values = c("red", "goldenrod", rep(c("darkgreen", "blue"), 20))) +
  # scale_shape_manual(values = c(rep(0:6, 20))) +
  coord_flip() +
  facet_grid(~ traitType2) +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16), strip.background = element_rect(color = "white", fill = "white"),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(plot = discrimValue_panels, filename = "medianDiscrimValue_orderedYC_panels.tif", device = "tiff",
       width = 7.5, height = 6.5, units = "in")

discrimValue_panels_eachTrait <- 
  ggplot(discrimValue_melt_traitType, mapping = aes(x = reorder(enviros, value2, median_), y = value)) +
  geom_col() +
  coord_flip() +
  facet_grid(~ variable) +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12, angle = 90), strip.background = element_rect(color = "white", fill = "white"),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(plot = discrimValue_panels_eachTrait, filename = "medianDiscrimValue_panelForEachTrait.tif", device = "tiff",
       width = 9, height = 6.5, units = "in")

discrimValue2 <- ggplot(discrimValue_allTraits_melt, mapping = aes(x = reorder(x = enviros, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  # coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Discriminabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = discrimValue2, filename = "medianDiscrimValue2.tif", device = "tiff",
       width = 5, height = 7.75, units = "in")


combinedDiscrimBoxplots <- ggarrange(discrimValue2, discrimValue_panels,
                                     widths = c(1,1.95),
                                     labels = c("A", "B"),
                                     ncol = 2, nrow = 1)

ggsave(combinedDiscrimBoxplots, filename = "combined_discriminability_boxplots.tif", device = "tiff",
       height = 8.5, width = 11, units = "in")

#### Rank correlation of discriminabilty of different types of traits/different traits ####
data_discrim <- discrimValue_allTraits[,-1]

cormat <- cor(data_discrim, use = "complete.obs", method = "spearman")

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
p.mat <- cor.mtest(data_discrim)
head(p.mat[, 1:5])


ggcorrplot_spearman <- ggcorrplot(cor(data_discrim, use = "complete.obs", method = "spearman"), p.mat = p.mat, type = "upper",
                         colors = c("darkred", "white", "steelblue"),
                         ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot, filename = "discriminability-ggCorrPlot-spearman.tif", device = "tiff", width = 6, height = 4.5, units = "in")

combinedDiscrimBoxplotAndCorrplot <- ggarrange(discrimValue2, ggcorrplot_spearman,
                                               widths = c(1,1.8),
                                               labels = c("A", "B"),
                                               ncol = 2, nrow = 1)

ggsave(combinedDiscrimBoxplotAndCorrplot, filename = "combined_discriminability_boxplots_and_correlation.tif", device = "tiff",
       height = 8, width = 11, units = "in")

#### Count instances of enviros showing in 5 least and 5 most discriminable for a trait ####
least_and_most_discrim_enviros <- matrix(NA, nrow = 0, ncol = 3)
colnames(least_and_most_discrim_enviros) <- c("trait", "least", "most")

for (trait in traits) {
  # trait = "Kernel width"
  
  traitDiscrim <- discrimValue_allTraits[,c("enviros",trait)]
  
  mostDiscrim <- traitDiscrim %>%
    top_n(5) %>%
    .$enviros
  
  leastDiscrim <- traitDiscrim %>%
    top_n(-5) %>%
    .$enviros
  
  enviros_rows <- c(rep(trait, 5), as.character(leastDiscrim), as.character(mostDiscrim))
  enviros_matrix <- matrix(enviros_rows, nrow = 5, ncol = 3)
  least_and_most_discrim_enviros <- rbind(least_and_most_discrim_enviros, enviros_matrix)
}

least_and_most_discrim_enviros <- as.data.frame(least_and_most_discrim_enviros)
write.csv(least_and_most_discrim_enviros, "least_and_most_discrim_enviros.csv", row.names = F)

least_and_most_discrim_enviros <- read.csv("least_and_most_discrim_enviros.csv")

least_and_most_discrim_enviros_melt <- melt(least_and_most_discrim_enviros, id.var = "trait")

myTable <- table(least_and_most_discrim_enviros_melt$variable, least_and_most_discrim_enviros_melt$value)
myTable_melt <- melt(myTable)

for (i in 1:nrow(myTable_melt)) {
  if (myTable_melt[i,"Var1"] == "least") {
    myTable_melt[i,"value"] <- -myTable_melt[i,"value"]
  }
}

discrim_enviros <- ggplot(myTable_melt, aes(x = reorder(Var2, value), y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "Environment", y = "") +
  scale_fill_manual(name = "Discriminability",
                    values = c("darkred", "steelblue")) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = c(0.225,0.9), legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  coord_flip(ylim = c(-12, 12))
ggsave(plot = discrim_enviros, filename = paste("./most_vs_least_discrim/discrim_enviros.tif",sep = ""),
       device = "tiff", width = 4, height = 4, units = "in")

discrim_enviros_byLocation <- ggplot(myTable_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "Environment", y = "") +
  scale_fill_manual(name = "Discriminability",
                    values = c("darkred", "steelblue")) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = c(0.225,0.9), legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  coord_flip(ylim = c(-12, 12))
ggsave(plot = discrim_enviros_byLocation, filename = paste("./most_vs_least_discrim/discrim_enviros_byLocation.tif",sep = ""),
       device = "tiff", width = 4, height = 4, units = "in")

myTable_melt$year <- myTable_melt$Var2
myTable_melt$year <- as.numeric(gsub(".*_", "", x = myTable_melt$year))

discrim_enviros_byYear <- ggplot(myTable_melt, aes(x = reorder(Var2, myTable_melt$year), y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "Environment", y = "") +
  scale_fill_manual(name = "Discriminability",
                    values = c("darkred", "steelblue")) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = c(0.225,0.9), legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  coord_flip(ylim = c(-12, 12))
ggsave(plot = discrim_enviros_byYear, filename = paste("./most_vs_least_discrim/discrim_enviros_byYear.tif",sep = ""),
       device = "tiff", width = 4, height = 4, units = "in")

lat_long_plDens <- read.csv("latitude_longitude_plantingDensity.csv", header = T)
lat_long_plDens <- lat_long_plDens[,-2]

colnames(lat_long_plDens) <- c("Var2", "Latitude", "Longitude", "Planting density")

myTable_melt <- merge(myTable_melt, lat_long_plDens)

discrim_enviros_byLatitude <- ggplot(myTable_melt, aes(x = reorder(Var2, myTable_melt$Latitude), y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "Environment", y = "") +
  scale_fill_manual(name = "Discriminability",
                    values = c("darkred", "steelblue")) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = c(0.225,0.9), legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  coord_flip(ylim = c(-12, 12))
ggsave(plot = discrim_enviros_byLatitude, filename = paste("./most_vs_least_discrim/discrim_enviros_byLatitude.tif",sep = ""),
       device = "tiff", width = 4, height = 4, units = "in")

discrim_enviros_byLongitude <- ggplot(myTable_melt, aes(x = reorder(Var2, myTable_melt$Longitude), y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "Environment", y = "") +
  scale_fill_manual(name = "Discriminability",
                    values = c("darkred", "steelblue")) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = c(0.225,0.9), legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  coord_flip(ylim = c(-12, 12))
ggsave(plot = discrim_enviros_byLongitude, filename = paste("./most_vs_least_discrim/discrim_enviros_byLongitude.tif",sep = ""),
       device = "tiff", width = 4, height = 4, units = "in")

discrim_enviros_byPlantDens <- ggplot(myTable_melt, aes(x = reorder(Var2, myTable_melt$`Planting density`), y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "Environment", y = "") +
  scale_fill_manual(name = "Discriminability",
                    values = c("darkred", "steelblue")) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = c(0.225,0.9), legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  coord_flip(ylim = c(-12, 12))
ggsave(plot = discrim_enviros_byPlantDens, filename = paste("./most_vs_least_discrim/discrim_enviros_byPlantDens.tif",sep = ""),
       device = "tiff", width = 4, height = 4, units = "in")

discrim_enviros_combined <- ggarrange(discrim_enviros_byLocation, discrim_enviros_byYear, discrim_enviros_byLongitude, discrim_enviros_byLatitude,
                                      ncol = 2, nrow = 2,
                                      labels = c("A", "B", "C", "D"))
ggsave(plot = discrim_enviros_combined, filename = paste("./most_vs_least_discrim/discrim_enviros_combined.tif",sep = ""),
       device = "tiff", width = 16, height = 16, units = "in")

colnames(least_and_most_discrim_enviros_melt) <- c("Trait", "Discrim", "Experiment")
lat_long_plDens <- read.csv("latitude_longitude_plantingDensity.csv", header = T)

least_most_discrim_enviros_melt_plus_lat_long_plDens <- merge(least_and_most_discrim_enviros_melt, lat_long_plDens)
least_most_discrim_enviros_melt_plus_lat_long_plDens$Discrim_value <- least_most_discrim_enviros_melt_plus_lat_long_plDens$Discrim
least_most_discrim_enviros_melt_plus_lat_long_plDens$Discrim_value <- gsub("most", "1", least_most_discrim_enviros_melt_plus_lat_long_plDens$Discrim_value)
least_most_discrim_enviros_melt_plus_lat_long_plDens$Discrim_value <- gsub("least", "-1", least_most_discrim_enviros_melt_plus_lat_long_plDens$Discrim_value)
least_most_discrim_enviros_melt_plus_lat_long_plDens$Discrim_value <- as.numeric(least_most_discrim_enviros_melt_plus_lat_long_plDens$Discrim_value)

least_most_discrim_enviros_melt_plus_lat_long_plDens$Trait <- as.character(least_most_discrim_enviros_melt_plus_lat_long_plDens$Trait)
least_most_discrim_enviros_melt_plus_lat_long_plDens$Trait <- factor(least_most_discrim_enviros_melt_plus_lat_long_plDens$Trait,
                                                                     levels = c("Anthesis", "Silking", 
                                                                                "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                                                                "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                                                                "Kernel length", "Kernel width", "Kernel thickness"))

#### Fig. 5
discrim_enviros_byLongitudeValue <-
  ggplot(least_most_discrim_enviros_melt_plus_lat_long_plDens, aes(x = reorder(Experiment, Longitude), y = Discrim_value)) + 
  geom_bar(stat = "identity", fill = "gray50") +
  scale_y_continuous(breaks = c(-10, -5, 0, 5, 10),
                     labels = c("10", "5", "0", "5", "10")) +
  labs(x = "Environment", y = "") +
  geom_abline(slope = 0, intercept = 0, size = 1.2) +
  annotate("text", label = "Most discriminating", x = 17, y = 10.5, size = 4, family = "Courier") +
  annotate("text", label = "Least discriminating", x = 17, y = -10.5, size = 4, family = "Courier") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 0), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))
ggsave(plot = discrim_enviros_byLongitudeValue, filename = "./most_vs_least_discrim/barplot_discrim_enviros_byLongitude.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")

#### Fig. 5 - bars colored by Midwest/non
least_most_discrim_enviros_melt_plus_lat_long_plDens$Region <- NA

for (i in 1:nrow(least_most_discrim_enviros_melt_plus_lat_long_plDens)) {
  if (least_most_discrim_enviros_melt_plus_lat_long_plDens[i, "Longitude"] > -94 &&
      least_most_discrim_enviros_melt_plus_lat_long_plDens[i, "Longitude"] < -84) {
    least_most_discrim_enviros_melt_plus_lat_long_plDens[i, "Region"] <- "Midwest"  
  }  else {
    least_most_discrim_enviros_melt_plus_lat_long_plDens[i, "Region"] <- "West/East"
  }
}

discrim_enviros_byLongitudeValue_colorMidwestVsOther <-
  ggplot(least_most_discrim_enviros_melt_plus_lat_long_plDens, aes(x = reorder(Experiment, Longitude), y = Discrim_value, fill = Region)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("steelblue", "darkred")) +
  scale_y_continuous(breaks = c(-10, -5, 0, 5, 10),
                     labels = c("10", "5", "0", "5", "10")) +
  labs(x = "Environment", y = "Number of traits") +
  geom_abline(slope = 0, intercept = 0, size = 1.2) +
  annotate("text", label = "Most discriminating", x = 17, y = 10.5, size = 4, family = "Courier") +
  annotate("text", label = "Least discriminating", x = 17, y = -10.5, size = 4, family = "Courier") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        legend.position = c(0.15, 0.85))
ggsave(plot = discrim_enviros_byLongitudeValue_colorMidwestVsOther, filename = "./most_vs_least_discrim/barplot_discrim_enviros_byLongitude_colorByMidwestVsOther.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")

#### Look at genotype performance in least/most discriminating enviros ####
selectedPhenoData <- read.csv("./selectedPhenoData.csv", header = T)
colnames(selectedPhenoData) <- c("Plot", "Enviro", "Rep", "Genotype", "Anthesis (GDU)", "Silking (GDU)", 
                                 "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                 "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                 "Kernel length", "Kernel width", "Kernel thickness")
selectedPhenoData_standardized <- selectedPhenoData
selectedPhenoData_standardized <- selectedPhenoData_standardized %>% transmute(
  Anthesis_st = (`Anthesis (GDU)` - mean(`Anthesis (GDU)`, na.rm = T)) / sd(`Anthesis (GDU)`, na.rm = T),
  Silking_st = (`Silking (GDU)` - mean(`Silking (GDU)`, na.rm = T)) / sd(`Silking (GDU)`, na.rm = T),
  PlantHeight_st = (`Plant height` - mean(`Plant height`, na.rm = T)) / sd(`Plant height`, na.rm = T),
  EarHeight_st = (`Ear height` - mean(`Ear height`, na.rm = T)) / sd(`Ear height`, na.rm = T),
  PlotWeight_st = (`Plot grain weight` - mean(`Plot grain weight`, na.rm = T)) / sd(`Plot grain weight`, na.rm = T),
  EarLength_st = (`Ear length` - mean(`Ear length`, na.rm = T)) / sd(`Ear length`, na.rm = T),
  EarWidth_st = (`Ear width` - mean(`Ear width`, na.rm = T)) / sd(`Ear width`, na.rm = T),
  KernelsPerRow_st = (`Kernels per row` - mean(`Kernels per row`, na.rm = T)) / sd(`Kernels per row`, na.rm = T),
  KRN_st = (`Kernel row number` - mean(`Kernel row number`, na.rm = T)) / sd(`Kernel row number`, na.rm = T),
  KernelWeight_st = (`Kernel weight` - mean(`Kernel weight`, na.rm = T)) / sd(`Kernel weight`, na.rm = T),
  KernelArea_st = (`Kernel area` - mean(`Kernel area`, na.rm = T)) / sd(`Kernel area`, na.rm = T),
  KernelLength_st = (`Kernel length` - mean(`Kernel length`, na.rm = T)) / sd(`Kernel length`, na.rm = T),
  KernelWidth_st = (`Kernel width` - mean(`Kernel width`, na.rm = T)) / sd(`Kernel width`, na.rm = T),
  KernelThickness_st = (`Kernel thickness` - mean(`Kernel thickness`, na.rm = T)) / sd(`Kernel thickness`, na.rm = T)
)  

selectedPhenoData_standardized <- cbind(selectedPhenoData[,1:4], selectedPhenoData_standardized)

lineMetadata <- read.csv("gxeInbredLines-Metadata.csv", header = T)
colnames(selectedPhenoData_standardized)[5:ncol(selectedPhenoData_standardized)] <- c("Anthesis", "Silking", 
                                                                                      "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                                                                      "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                                                                      "Kernel length", "Kernel width", "Kernel thickness")
# dir.create("./stability_in_most_least_discrim")

inbreds <- unique(selectedPhenoData_standardized$Genotype)
traits <- colnames(selectedPhenoData_standardized[ ,5:ncol(selectedPhenoData_standardized)])

slopes_intercepts_grouped_all_traits <- tibble()

#### Environmental index (means per enviro)
means_by_exp <- aggregate(selectedPhenoData_standardized,
                          by = list(selectedPhenoData_standardized$Enviro), FUN = mean, na.rm = T)
means_by_exp$Enviro <- means_by_exp$Group.1
means_by_exp <- means_by_exp[,-c(1, 2, 4, 5)]


for (trait in traits) {
  # trait = "Kernel row number"
  
  #### Most and least discrim enviros
  traitDiscrim <- discrimValue_allTraits[,c("enviros",trait)]
  colnames(traitDiscrim)[1] <- "Enviro"
  
  mostDiscrim <- traitDiscrim %>%
    top_n(5) %>%
    mutate(Category = "most")
  
  leastDiscrim <- traitDiscrim %>%
    top_n(-5) %>%
    mutate(Category = "least")
  
  most_least_discrim <- bind_rows(mostDiscrim, leastDiscrim)
  
  means_by_exp_trait <- data.frame(Enviro = means_by_exp[,"Enviro"], trait = means_by_exp[,trait])
  colnames(means_by_exp_trait) <- c("Enviro", paste0("EI_", trait))
  
  most_least_discrim <- left_join(most_least_discrim, means_by_exp_trait)
  
  data_most_least_discrim <- selectedPhenoData_standardized[which(selectedPhenoData_standardized$Enviro %in% most_least_discrim$Enviro),]
  data_most_least_discrim_trait <- cbind(data_most_least_discrim[,c(1:4)], data_most_least_discrim[,trait])
  colnames(data_most_least_discrim_trait)[5] <- trait
  
  dir.create(paste0("./stability_in_most_least_discrim/", trait))
  
  inb_loc_means_all <- matrix(NA, nrow = 0, ncol = 4)
  
  for (inbred in inbreds) {
    # inbred = "740"
    inbredDat <- data_most_least_discrim_trait[data_most_least_discrim_trait$Genotype == inbred, c("Enviro", trait)]
    
    if (all(is.na(inbredDat[,trait]))) {
      reg[inbred, ] <- NA
      next
    }
    inbredDat <- inbredDat[!is.na(inbredDat[,trait]), ]
    
    inb_loc_means <- aggregate(inbredDat, by = list(inbredDat$Enviro), FUN = mean, na.rm = T)
    inb_loc_means <- inb_loc_means[,-2]
    colnames(inb_loc_means)[1]  <- "Enviro"
    
    inb_loc_means$Inbred <- inbred
    
    inb_loc_means <- left_join(inb_loc_means, most_least_discrim, by = "Enviro")
    colnames(inb_loc_means) <- c("Enviro", trait, "Genotype", "Discriminability value", "Category", "Enviro_Index")
    inb_loc_means <- inb_loc_means[,c(1,5,3,4,6,2)]
    
    inb_loc_means_all <- rbind(inb_loc_means_all, inb_loc_means)
  }
  
  inb_loc_means_all_metadata <- left_join(inb_loc_means_all,lineMetadata, by = "Genotype")
  
  stability_regression <-
    ggplot(inb_loc_means_all_metadata, aes(x = `Enviro_Index`, y = get(trait))) +
    geom_point(aes(color = Ex.PVP), shape = 1) +
    # geom_smooth(method = lm, se = F, fullrange = F, aes(color = Genotype)) +
    geom_smooth(method = lm, se = F, fullrange = T, aes(color = Ex.PVP)) +
    # geom_smooth(method = lm, se = F, fullrange = T, color = "black") +
    coord_cartesian(ylim = c(-3, 4)) +
    facet_grid(.~inb_loc_means_all_metadata$Category
               # , scales = "free_x"
    ) +
    labs(x = "Environmental Index", y = trait) +
    theme_bw() +
    theme(text = element_text(family = "Courier"),
          axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10), axis.ticks.y = element_blank(),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_rect(color = "white", fill = "white"))
  ggsave(plot = stability_regression, filename = paste0("./stability_in_most_least_discrim/", trait, "/stability_regression_plot.tif"),
         device = "tiff", width = 6.5, height = 5, units = "in")
  
  fitted_models_grouped <- inb_loc_means_all_metadata %>%
    group_by(Category, Ex.PVP) %>%
    do(model = lm(get(trait) ~ Enviro_Index, data = .))
  
  slopes_intercepts_grouped <- tidy(x = fitted_models_grouped, object = model)
  
  # fitted_models_ungrouped <- inb_loc_means_all_metadata %>%
  #   group_by(Category, Genotype) %>%
  #   do(model = lm(get(trait) ~ Enviro_Index, data = .))
  # slopes_intercepts_ungrouped <- tidy(x = fitted_models_ungrouped, object = model)
  
  slopes_intercepts_grouped$Trait <- trait 
  slopes_intercepts_grouped_all_traits <- bind_rows(slopes_intercepts_grouped_all_traits, slopes_intercepts_grouped)
}
write.csv(slopes_intercepts_grouped_all_traits, "slopes_intercepts_grouped_all_traits.csv", row.names = F)

#### Look at genotype performance in lowest/highest EI enviros ####

lineMetadata <- read.csv("gxeInbredLines-Metadata.csv", header = T)
colnames(selectedPhenoData_standardized)[5:ncol(selectedPhenoData_standardized)] <- c("Anthesis", "Silking", 
                                                                                      "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                                                                      "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                                                                      "Kernel length", "Kernel width", "Kernel thickness")
# dir.create("./stability_in_highest_lowest_EI")

inbreds <- unique(selectedPhenoData_standardized$Genotype)
traits <- colnames(selectedPhenoData_standardized[ ,5:ncol(selectedPhenoData_standardized)])

slopes_intercepts_grouped_all_traits <- tibble()

#### Environmental index (means per enviro)
means_by_exp <- aggregate(selectedPhenoData_standardized,
                          by = list(selectedPhenoData_standardized$Enviro), FUN = mean, na.rm = T)
means_by_exp$Enviro <- means_by_exp$Group.1
means_by_exp <- means_by_exp[,-c(1, 2, 4, 5)]


for (trait in traits) {
  # trait = "Anthesis (GDU)"
  
  means_by_exp_trait <- data.frame(Enviro = means_by_exp[,"Enviro"], trait = means_by_exp[,trait])
  colnames(means_by_exp_trait) <- c("Enviro", paste0("EI_", trait))
  
  #### highest lowest EI enviros
  
  highestEI <- means_by_exp_trait %>%
    top_n(5) %>%
    mutate(Category = "highest")
  
  lowestEI <- means_by_exp_trait %>%
    top_n(-5) %>%
    mutate(Category = "lowest")
  
  highest_lowest_EI <- bind_rows(highestEI, lowestEI)
  
  data_highest_lowest_EI <- selectedPhenoData_standardized[which(selectedPhenoData_standardized$Enviro %in% highest_lowest_EI$Enviro),]
  data_highest_lowest_EI_trait <- cbind(data_highest_lowest_EI[,c(1:4)], data_highest_lowest_EI[,trait])
  colnames(data_highest_lowest_EI_trait)[5] <- trait
  
  dir.create(paste0("./stability_in_highest_lowest_EI/", trait))
  
  inb_loc_means_all <- matrix(NA, nrow = 0, ncol = 4)
  
  for (inbred in inbreds) {
    # inbred = "740"
    inbredDat <- data_highest_lowest_EI_trait[data_highest_lowest_EI_trait$Genotype == inbred, c("Enviro", trait)]
    
    if (all(is.na(inbredDat[,trait]))) {
      reg[inbred, ] <- NA
      next
    }
    inbredDat <- inbredDat[!is.na(inbredDat[,trait]), ]
    
    inb_loc_means <- aggregate(inbredDat, by = list(inbredDat$Enviro), FUN = mean, na.rm = T)
    inb_loc_means <- inb_loc_means[,-2]
    colnames(inb_loc_means)[1]  <- "Enviro"
    
    inb_loc_means$Inbred <- inbred
    
    inb_loc_means <- left_join(inb_loc_means, highest_lowest_EI, by = "Enviro")
    colnames(inb_loc_means) <- c("Enviro", trait, "Genotype", "Enviro_Index", "Category")
    inb_loc_means <- inb_loc_means[,c(1,5,3,4,2)]
    
    inb_loc_means_all <- rbind(inb_loc_means_all, inb_loc_means)
  }
  
  inb_loc_means_all_metadata <- left_join(inb_loc_means_all,lineMetadata, by = "Genotype")
  
  inb_loc_means_all_metadata$Category <- factor(inb_loc_means_all_metadata$Category, levels = c("lowest", "highest"))
  
  stability_regression <-
    ggplot(inb_loc_means_all_metadata, aes(x = Enviro_Index, y = get(trait))) +
    geom_point(aes(color = Ex.PVP), shape = 1) +
    geom_smooth(method = lm, se = F, fullrange = T, aes(color = Ex.PVP)) +
    # geom_smooth(method = lm, se = F, fullrange = F, aes(color = Genotype)) +
    # geom_smooth(method = lm, se = F, fullrange = T, color = "black") +
    coord_cartesian(ylim = c(-3,4)) +
    facet_grid(.~Category, scales = "free_x") +
    labs(x = "Environmental Index", y = trait) +
    theme_bw() +
    theme(text = element_text(family = "Courier"),
          axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10), axis.ticks.y = element_blank(),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_rect(color = "white", fill = "white"))
  
  ggsave(plot = stability_regression, filename = paste0("./stability_in_highest_lowest_EI/", trait, "/stability_regression_plot.tif"),
         device = "tiff", width = 6.5, height = 5, units = "in")
  
  fitted_models_grouped <- inb_loc_means_all_metadata %>%
    group_by(Category, Ex.PVP) %>%
    do(model = lm(get(trait) ~ Enviro_Index, data = .))
  
  slopes_intercepts_grouped <- tidy(x = fitted_models_grouped, object = model)
  
  # fitted_models_ungrouped <- inb_loc_means_all_metadata %>%
  #   group_by(Category, Genotype) %>%
  #   do(model = lm(get(trait) ~ Enviro_Index, data = .))
  # slopes_intercepts_ungrouped <- tidy(x = fitted_models_ungrouped, object = model)
  
  slopes_intercepts_grouped$Trait <- trait 
  slopes_intercepts_grouped_all_traits <- bind_rows(slopes_intercepts_grouped_all_traits, slopes_intercepts_grouped)
}
write.csv(slopes_intercepts_grouped_all_traits, "slopes_intercepts_grouped_all_traits_EI.csv", row.names = F)

#### Compare season accumulated rain in most and least discriminable enviros ####
rainSummary <- read.csv("./weatherData/rainSeasonSummary.csv", header = T)

sigDiff_seasonalRain <- matrix(NA, nrow = 0, ncol = 4)
colnames(sigDiff_seasonalRain) <- c("trait", "mean_mostDiscrim", "mean_leastDiscrim", "pvalue")

seasonalRainData_compare_allTraits <- matrix(NA, nrow = 0, ncol = 4)
colnames(seasonalRainData_compare_allTraits) <- c("experiment", "seasonal_accumulated_rainfall",
                                                  "category", "trait")

myBoxplot <- function(data, category, weather, trait, experiment) {
  # plot <-
  ggplot(data = data, aes(factor(trait), weather)) + 
    geom_boxplot(aes(fill = category)) +
    # geom_jitter(aes(color = experiment, shape = experiment)) +
    # scale_color_manual(values = rep(c("red", "goldenrod", "darkgreen", "blue"), 10)) +
    # scale_shape_manual(values = c(rep(0:6, 20))) +
    # geom_text(data = data, aes(category, weather, label = experiment), position = position_jitter(w=1, h=.2)) +
    # geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Category", y = colnames(data)[i], title = pval) +
    theme_bw() +
    theme(text = element_text(family = "Courier"),
          axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12), 
          axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  # ggsave(plot = plot, filename = paste("./most_vs_least_discrim/", trait, "_", colnames(data)[i], ".tif",sep = ""), 
  # device = "tiff", width = 7, height = 4, units = "in")
}

for (trait in traits) {
  # trait = "Anthesis (GDU)"
  
  traitDiscrim <- discrimValue_allTraits[,c("enviros",trait)]
  
  mostDiscrim <- traitDiscrim %>%
    top_n(5) %>%
    .$enviros
  
  leastDiscrim <- traitDiscrim %>%
    top_n(-5) %>%
    .$enviros
  
  rainData_mostDiscrim <- rainSummary %>%
    filter(experiment %in% mostDiscrim)
  
  rainData_leastDiscrim <- rainSummary %>%
    filter(experiment %in% leastDiscrim)
  
  rainData_mostDiscrim$category <- "Most discriminating"
  rainData_leastDiscrim$category <- "Least discriminating"
  
  rainData_compare <- rbind(rainData_mostDiscrim, rainData_leastDiscrim)
  rainData_compare$trait <- trait
  
  tTest <- t.test(rainData_mostDiscrim$total_season_rain, rainData_leastDiscrim$total_season_rain)
  pval <- paste0("p = ", myround(tTest$p.value,3))
  
  row <- c(trait, mean(rainData_mostDiscrim$total_season_rain, na.rm = T), mean(rainData_leastDiscrim$total_season_rain, na.rm = T) , tTest$p.value)  
  sigDiff_seasonalRain <- rbind(sigDiff_seasonalRain, row)
  
  seasonalRainData_compare_allTraits <- rbind(seasonalRainData_compare_allTraits, rainData_compare)
} 

seasonalRainData_compare_allTraits2 <- unique(seasonalRainData_compare_allTraits) # remove identical rows

seasonalRainData_compare_allTraits_melt <- melt(seasonalRainData_compare_allTraits2,
                                                id.vars = c("trait", "category", "experiment"))

seasonalRainData_compare_allTraits_melt$trait <- factor(seasonalRainData_compare_allTraits_melt$trait,
                                                        levels = unique(seasonalRainData_compare_allTraits_melt$trait))

seasonal_rain <- ggplot(data = seasonalRainData_compare_allTraits_melt, aes(x = factor(category), y = value, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 150,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ variable, scales = "free_x", switch = "y")
ggsave(plot = seasonal_rain, filename = paste("./most_vs_least_discrim/seasonal_rain.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")

#### Critical period rain #### 
silking_rain_summary <- read.csv("weatherData/silking_rain_summary.csv", header = T)

planting_to_silking_rain_summary <- read.csv("weatherData/planting_to_silking_rain_summary.csv", header = T)

sigDiff_silking_rain <- matrix(NA, nrow = 0, ncol = 5)
colnames(sigDiff_silking_rain) <- c("days_before_after_flowering", "trait", "mean_mostDiscrim", "mean_leastDiscrim", "pvalue")

silking_rain_compare_allTraits <- matrix(NA, nrow = 0, ncol = 4)
colnames(silking_rain_compare_allTraits) <- c("experiment", "silking_accumulated_rainfall",
                                              "category", "trait")

sigDiff_planting_to_silking_rain <- matrix(NA, nrow = 0, ncol = 5)
colnames(sigDiff_planting_to_silking_rain) <- c("days_before_after_flowering", "trait", "mean_mostDiscrim", "mean_leastDiscrim", "pvalue")

planting_to_silking_rain_compare_allTraits <- matrix(NA, nrow = 0, ncol = 4)
colnames(planting_to_silking_rain_compare_allTraits) <- c("experiment", "silking_accumulated_rainfall",
                                                          "category", "trait")

days_before_after_flowering <- c(7, 10, 14)

for (trait in traits) {
  # trait = "Anthesis (GDU)"
  
  traitDiscrim <- discrimValue_allTraits[,c("enviros",trait)]
  
  mostDiscrim <- traitDiscrim %>%
    top_n(5) %>%
    .$enviros
  
  leastDiscrim <- traitDiscrim %>%
    top_n(-5) %>%
    .$enviros
  
  for (d in days_before_after_flowering) {
    # d = 14
    ## Silking period
    silking_rain_summary_period <- subset(silking_rain_summary, silking_rain_summary$days_before_after_silking == d)
    
    rainDataSilking_mostDiscrim <- silking_rain_summary_period %>%
      filter(experiment %in% mostDiscrim)
    
    rainDataSilking_leastDiscrim <- silking_rain_summary_period %>%
      filter(experiment %in% leastDiscrim)
    
    rainDataSilking_mostDiscrim$category <- "Most discriminating"
    rainDataSilking_leastDiscrim$category <- "Least discriminating"
    
    rainDataSilking_compare <- rbind(rainDataSilking_mostDiscrim, rainDataSilking_leastDiscrim)
    rainDataSilking_compare$trait <- trait
    
    silking_rain_compare_allTraits <- rbind(silking_rain_compare_allTraits, rainDataSilking_compare)
    
    tTest <- t.test(rainDataSilking_mostDiscrim[,3], rainDataSilking_leastDiscrim[,3])
    
    row <- c(d, trait, 
             mean(rainDataSilking_mostDiscrim[,3], na.rm = T), 
             mean(rainDataSilking_leastDiscrim[,3], na.rm = T) , 
             tTest$p.value)  
    sigDiff_silking_rain <- rbind(sigDiff_silking_rain, row)
    
    ## Planting to silking period
    planting_to_silking_rain_summary_period <- subset(planting_to_silking_rain_summary, planting_to_silking_rain_summary$days_before_after_silking == d)
    
    rainDataPlantingToSilking_mostDiscrim <- planting_to_silking_rain_summary_period %>%
      filter(experiment %in% mostDiscrim)
    
    rainDataPlantingToSilking_leastDiscrim <- planting_to_silking_rain_summary_period %>%
      filter(experiment %in% leastDiscrim)
    
    rainDataPlantingToSilking_mostDiscrim$category <- "Most discriminating"
    rainDataPlantingToSilking_leastDiscrim$category <- "Least discriminating"
    
    rainDataPlantingToSilking_compare <- rbind(rainDataPlantingToSilking_mostDiscrim, rainDataPlantingToSilking_leastDiscrim)
    rainDataPlantingToSilking_compare$trait <- trait
    
    planting_to_silking_rain_compare_allTraits <- rbind(planting_to_silking_rain_compare_allTraits, rainDataPlantingToSilking_compare)
    
    tTest <- t.test(rainDataPlantingToSilking_mostDiscrim[,3], rainDataPlantingToSilking_leastDiscrim[,3])
    
    row <- c(d, trait, 
             mean(rainDataPlantingToSilking_mostDiscrim[,3], na.rm = T), 
             mean(rainDataPlantingToSilking_leastDiscrim[,3], na.rm = T) , 
             tTest$p.value)  
    sigDiff_planting_to_silking_rain <- rbind(sigDiff_planting_to_silking_rain, row)
  } 
} 

silking_rain <- ggplot(data = silking_rain_compare_allTraits, aes(x = factor(category), y = total_rain, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 100,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = silking_rain, filename = paste("./most_vs_least_discrim/silking_rain.tif",sep = ""),
       device = "tiff", width = 9, height = 6.5, units = "in")

planting_to_silking_rain <- ggplot(data = planting_to_silking_rain_compare_allTraits, aes(x = factor(category), y = total_rain, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 500,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = planting_to_silking_rain, filename = paste("./most_vs_least_discrim/planting_to_silking_rain.tif",sep = ""),
       device = "tiff", width = 9, height = 6.5, units = "in")




#### Compare min/median/max temperature in most and least discriminable enviros ####
tempSummary <- read.csv("./weatherData/tempSeasonSummary.csv", header = T)

sigDiff_seasonaltemp <- matrix(NA, nrow = 0, ncol = 5)
colnames(sigDiff_seasonaltemp) <- c("trait", "weather", "mean_mostDiscrim", "mean_leastDiscrim", "pvalue")

seasonaltempData_compare_allTraits <- matrix(NA, nrow = 0, ncol = 6)
colnames(seasonaltempData_compare_allTraits) <- c("experiment", "min_temp", "median_temp", "max_temp",
                                                  "category", "trait")

for (trait in traits) {
  # trait = "Anthesis (GDU)"
  
  traitDiscrim <- discrimValue_allTraits[,c("enviros",trait)]
  
  mostDiscrim <- traitDiscrim %>%
    top_n(5) %>%
    .$enviros
  
  leastDiscrim <- traitDiscrim %>%
    top_n(-5) %>%
    .$enviros
  
  tempData_mostDiscrim <- tempSummary %>%
    filter(experiment %in% mostDiscrim)
  
  tempData_leastDiscrim <- tempSummary %>%
    filter(experiment %in% leastDiscrim)
  
  tempData_mostDiscrim$category <- "Most discriminating"
  tempData_leastDiscrim$category <- "Least discriminating"
  
  tempData_compare <- rbind(tempData_mostDiscrim, tempData_leastDiscrim)
  tempData_compare$trait <- trait
  
  for (i in 2:4) {
    # i = 4
    
    tTest <- t.test(tempData_mostDiscrim[,i], tempData_leastDiscrim[,i])
    
    row <- c(trait, colnames(tempData_compare)[i], mean(tempData_mostDiscrim[,i], na.rm = T), mean(tempData_leastDiscrim[,i], na.rm = T), tTest$p.value)  
    sigDiff_seasonaltemp <- rbind(sigDiff_seasonaltemp, row)
  } 
  seasonaltempData_compare_allTraits <- rbind(seasonaltempData_compare_allTraits, tempData_compare)
} 

seasonaltempData_compare_allTraits_melt <- melt(seasonaltempData_compare_allTraits,
                                                id.vars = c("trait", "category", "experiment"))

seasonaltempData_compare_allTraits_melt$trait <- factor(seasonaltempData_compare_allTraits_melt$trait,
                                                        levels = unique(seasonaltempData_compare_allTraits_melt$trait))

seasonal_temp <- ggplot(data = seasonaltempData_compare_allTraits_melt, aes(x = factor(category), y = value, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 25,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ variable, switch = "y")
ggsave(plot = seasonal_temp, filename = paste("./most_vs_least_discrim/seasonal_temp.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")



#### Critical period temp #### 
silking_temp_summary <- read.csv("weatherData/silking_temp_summary.csv", header = T, na.strings = c("NA", "Inf", "-Inf"))

sigDiff_silking_temp <- matrix(NA, nrow = 0, ncol = 6)
colnames(sigDiff_silking_temp) <- c("days_before_after_flowering", "weather", "trait", "mean_mostDiscrim", "mean_leastDiscrim", "pvalue")

silking_temp_compare_allTraits <- matrix(NA, nrow = 0, ncol = 7)
colnames(silking_temp_compare_allTraits) <- c("experiment", "days_before_after_flowering", 
                                              "min_temp", "median_temp", "max_temp",
                                              "category", "trait")

planting_to_silking_temp_summary <- read.csv("weatherData/planting_to_silking_temp_summary.csv", header = T, na.strings = c("NA", "Inf", "-Inf"))

sigDiff_planting_to_silking_temp <- matrix(NA, nrow = 0, ncol = 6)
colnames(sigDiff_planting_to_silking_temp) <- c("days_before_after_flowering", "weather", "trait", "mean_mostDiscrim", "mean_leastDiscrim", "pvalue")

planting_to_silking_temp_compare_allTraits <- matrix(NA, nrow = 0, ncol = 7)
colnames(planting_to_silking_temp_compare_allTraits) <- c("experiment", "days_before_after_flowering", 
                                                          "min_temp", "median_temp", "max_temp",
                                                          "category", "trait")
days_before_after_flowering <- c(7, 10, 14)

for (trait in traits) {
  # trait = "Anthesis (GDU)"
  
  traitDiscrim <- discrimValue_allTraits[,c("enviros",trait)]
  
  mostDiscrim <- traitDiscrim %>%
    top_n(5) %>%
    .$enviros
  
  leastDiscrim <- traitDiscrim %>%
    top_n(-5) %>%
    .$enviros
  
  for (d in days_before_after_flowering) {
    # d = 14
    
    ## Silking
    silking_temp_summary_period <- subset(silking_temp_summary, silking_temp_summary$days_before_after_silking == d)
    
    tempDataSilking_mostDiscrim <- silking_temp_summary_period %>%
      filter(experiment %in% mostDiscrim)
    
    tempDataSilking_leastDiscrim <- silking_temp_summary_period %>%
      filter(experiment %in% leastDiscrim)
    
    tempDataSilking_mostDiscrim$category <- "Most discriminating"
    tempDataSilking_leastDiscrim$category <- "Least discriminating"
    
    tempDataSilking_compare <- rbind(tempDataSilking_mostDiscrim, tempDataSilking_leastDiscrim)
    tempDataSilking_compare$trait <- trait
    
    for (i in 3:5) {
      # i = 3
      tTest <- t.test(tempDataSilking_mostDiscrim[,i], tempDataSilking_leastDiscrim[,i])
      
      row <- c(d, trait, colnames(tempDataSilking_compare)[i],
               mean(tempDataSilking_mostDiscrim[,i], na.rm = T), 
               mean(tempDataSilking_leastDiscrim[,i], na.rm = T) , 
               tTest$p.value)
      sigDiff_silking_temp <- rbind(sigDiff_silking_temp, row)
    } 
    silking_temp_compare_allTraits <- rbind(silking_temp_compare_allTraits, tempDataSilking_compare)
    
    ## Planting to silking
    planting_to_silking_temp_summary_period <- subset(planting_to_silking_temp_summary, planting_to_silking_temp_summary$days_before_after_silking == d)
    
    tempDataPlantingToSilking_mostDiscrim <- planting_to_silking_temp_summary_period %>%
      filter(experiment %in% mostDiscrim)
    
    tempDataPlantingToSilking_leastDiscrim <- planting_to_silking_temp_summary_period %>%
      filter(experiment %in% leastDiscrim)
    
    tempDataPlantingToSilking_mostDiscrim$category <- "Most discriminating"
    tempDataPlantingToSilking_leastDiscrim$category <- "Least discriminating"
    
    tempDataPlantingToSilking_compare <- rbind(tempDataPlantingToSilking_mostDiscrim, tempDataPlantingToSilking_leastDiscrim)
    tempDataPlantingToSilking_compare$trait <- trait
    
    for (i in 3:5) {
      # i = 3
      tTest <- t.test(tempDataPlantingToSilking_mostDiscrim[,i], tempDataPlantingToSilking_leastDiscrim[,i])
      
      row <- c(d, trait, colnames(tempDataPlantingToSilking_compare)[i],
               mean(tempDataPlantingToSilking_mostDiscrim[,i], na.rm = T), 
               mean(tempDataPlantingToSilking_leastDiscrim[,i], na.rm = T) , 
               tTest$p.value)
      sigDiff_planting_to_silking_temp <- rbind(sigDiff_planting_to_silking_temp, row)
    } 
    planting_to_silking_temp_compare_allTraits <- rbind(planting_to_silking_temp_compare_allTraits, tempDataPlantingToSilking_compare)
  }
} 

silking_temp_min <- ggplot(data = silking_temp_compare_allTraits, aes(x = factor(category), y = min_temp, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 5,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = silking_temp_min, filename = paste("./most_vs_least_discrim/silking_temp_min.tif",sep = ""),
       device = "tiff", width = 9, height = 6.5, units = "in")

silking_temp_median <- ggplot(data = silking_temp_compare_allTraits, aes(x = factor(category), y = median_temp, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 15,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = silking_temp_median, filename = paste("./most_vs_least_discrim/silking_temp_median.tif",sep = ""),
       device = "tiff", width = 9, height = 6.5, units = "in")

silking_temp_max <- ggplot(data = silking_temp_compare_allTraits, aes(x = factor(category), y = max_temp, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 28,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = silking_temp_max, filename = paste("./most_vs_least_discrim/silking_temp_max.tif",sep = ""),
       device = "tiff", width = 9, height = 6.5, units = "in")

## Planting to silking
planting_to_silking_temp_min <- ggplot(data = planting_to_silking_temp_compare_allTraits, aes(x = factor(category), y = min_temp, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 10,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = planting_to_silking_temp_min, filename = paste("./most_vs_least_discrim/planting_to_silking_temp_min.tif",sep = ""),
       device = "tiff", width = 9, height = 6.5, units = "in")

planting_to_silking_temp_median <- ggplot(data = planting_to_silking_temp_compare_allTraits, aes(x = factor(category), y = median_temp, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 15,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = planting_to_silking_temp_median, filename = paste("./most_vs_least_discrim/planting_to_silking_temp_median.tif",sep = ""),
       device = "tiff", width = 9, height = 6.5, units = "in")

planting_to_silking_temp_max <- ggplot(data = planting_to_silking_temp_compare_allTraits, aes(x = factor(category), y = max_temp, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 28,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = planting_to_silking_temp_max, filename = paste("./most_vs_least_discrim/planting_to_silking_temp_max.tif",sep = ""),
       device = "tiff", width = 9, height = 6.5, units = "in")


#### # of hrs: Temperature > 30, 25, 22, < 22 ####

hours_temp_hum_summary <- read.csv("./weatherData/hours_temp_hum_summary.csv", header = T)

sigDiff_hours_temp_hum <- matrix(NA, nrow = 0, ncol = 5)
colnames(sigDiff_hours_temp_hum) <- c("trait", "weather", "mean_mostDiscrim", "mean_leastDiscrim", "pvalue")

hoursTempHum_compare_allTraits <- matrix(NA, nrow = 0, ncol = 17)
colnames(hoursTempHum_compare_allTraits) <- c("experiment", "Tg30", "Tg25", "Tg22", "Tl22", "Hg90", "Hg75", "Hg50", "Hl50",
                                              "Tg25Hl90", "Tg22Hl90", "Tg25Hl75", "Tg22Hl75", "Tg25Hl50", "Tg22Hl50",  
                                              "category", "trait")

for (trait in traits) {
  # trait = "Anthesis (GDU)"
  
  traitDiscrim <- discrimValue_allTraits[,c("enviros",trait)]
  
  mostDiscrim <- traitDiscrim %>%
    top_n(5) %>%
    .$enviros
  
  leastDiscrim <- traitDiscrim %>%
    top_n(-5) %>%
    .$enviros
  
  hours_temp_hum_mostDiscrim <- hours_temp_hum_summary %>%
    filter(experiment %in% mostDiscrim)
  
  hours_temp_hum_leastDiscrim <- hours_temp_hum_summary %>%
    filter(experiment %in% leastDiscrim)
  
  hours_temp_hum_mostDiscrim$category <- "Most discriminating"
  hours_temp_hum_leastDiscrim$category <- "Least discriminating"
  
  hours_temp_hum_compare <- rbind(hours_temp_hum_mostDiscrim, hours_temp_hum_leastDiscrim)
  hours_temp_hum_compare$trait <- trait
  
  for (i in 2:15) {
    # i = 2
    
    tTest <- t.test(hours_temp_hum_mostDiscrim[,i], hours_temp_hum_leastDiscrim[,i])
    
    row <- c(trait, colnames(hours_temp_hum_compare)[i], mean(hours_temp_hum_mostDiscrim[,i], na.rm = T), mean(hours_temp_hum_leastDiscrim[,i], na.rm = T), tTest$p.value)  
    sigDiff_hours_temp_hum <- rbind(sigDiff_hours_temp_hum, row)
  } 
  hoursTempHum_compare_allTraits <- rbind(hoursTempHum_compare_allTraits, hours_temp_hum_compare)
} 

if (all(is.na(hoursTempHum_compare_allTraits))) {
  
}
hoursTempHum_compare_allTraits <- hoursTempHum_compare_allTraits[,colSums(hoursTempHum_compare_allTraits != 0) > 0]

hoursTempHum_compare_allTraits_melt <- melt(hoursTempHum_compare_allTraits,
                                            id.vars = c("trait", "category", "experiment"))

hoursTempHum_compare_allTraits_melt$trait <- factor(hoursTempHum_compare_allTraits_melt$trait,
                                                    levels = unique(hoursTempHum_compare_allTraits_melt$trait))

hoursTempHum <- ggplot(data = hoursTempHum_compare_allTraits_melt, aes(x = factor(category), y = value, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 25,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ variable, scales = "free_x", switch = "y")
ggsave(plot = hoursTempHum, filename = paste("./most_vs_least_discrim/hoursTempHum.tif",sep = ""),
       device = "tiff", width = 12, height = 6.5, units = "in")



#### Critical periods
silking_hours_temp_hum_summary <- read.csv("./weatherData/silking_hours_temp_hum_summary.csv", header = T)
planting_to_silking_hours_temp_hum_summary <- read.csv("./weatherData/planting_to_silking_hours_temp_hum_summary.csv", header = T)

sigDiff_silking_hours_temp_hum <- matrix(NA, nrow = 0, ncol = 6)
colnames(sigDiff_silking_hours_temp_hum) <- c("days_before_after_silking", "trait", "weather", "mean_mostDiscrim", "mean_leastDiscrim", "pvalue")

silking_hoursTempHum_compare_allTraits <- matrix(NA, nrow = 0, ncol = 18)
colnames(silking_hoursTempHum_compare_allTraits) <- c("experiment", "days_before_after_silking", "Tg30", "Tg25", "Tg22", "Tl22", "Hg90", "Hg75", "Hg50", "Hl50",
                                                      "Tg25Hl90", "Tg22Hl90", "Tg25Hl75", "Tg22Hl75", "Tg25Hl50", "Tg22Hl50",  
                                                      "category", "trait")

sigDiff_planting_to_silking_hours_temp_hum <- matrix(NA, nrow = 0, ncol = 6)
colnames(sigDiff_planting_to_silking_hours_temp_hum) <- c("days_before_after_silking", "trait", "weather", "mean_mostDiscrim", "mean_leastDiscrim", "pvalue")

planting_to_silking_hoursTempHum_compare_allTraits <- matrix(NA, nrow = 0, ncol = 18)
colnames(planting_to_silking_hoursTempHum_compare_allTraits) <- c("experiment", "days_before_after_silking", "Tg30", "Tg25", "Tg22", "Tl22", "Hg90", "Hg75", "Hg50", "Hl50",
                                                                  "Tg25Hl90", "Tg22Hl90", "Tg25Hl75", "Tg22Hl75", "Tg25Hl50", "Tg22Hl50",  
                                                                  "category", "trait")

for (trait in traits) {
  # trait = "Anthesis (GDU)"
  
  traitDiscrim <- discrimValue_allTraits[,c("enviros",trait)]
  
  mostDiscrim <- traitDiscrim %>%
    top_n(5) %>%
    .$enviros
  
  leastDiscrim <- traitDiscrim %>%
    top_n(-5) %>%
    .$enviros
  
  for (d in days_before_after_flowering) {
    # d = 14
    
    ## Silking
    silking_hours_temp_hum <- subset(silking_hours_temp_hum_summary, silking_hours_temp_hum_summary$days_before_after_silking == d)
    
    silking_hours_temp_hum_mostDiscrim <- silking_hours_temp_hum %>%
      filter(experiment %in% mostDiscrim)
    
    silking_hours_temp_hum_leastDiscrim <- silking_hours_temp_hum %>%
      filter(experiment %in% leastDiscrim)
    
    silking_hours_temp_hum_mostDiscrim$category <- "Most discriminating"
    silking_hours_temp_hum_leastDiscrim$category <- "Least discriminating"
    
    silking_hours_temp_hum_compare <- rbind(silking_hours_temp_hum_mostDiscrim, silking_hours_temp_hum_leastDiscrim)
    silking_hours_temp_hum_compare$trait <- trait
    
    for (i in 3:16) {
      # i = 3
      
      tTest <- t.test(silking_hours_temp_hum_mostDiscrim[,i], silking_hours_temp_hum_leastDiscrim[,i])
      
      row <- c(d, trait, colnames(silking_hours_temp_hum_compare)[i], mean(silking_hours_temp_hum_mostDiscrim[,i], na.rm = T), mean(silking_hours_temp_hum_leastDiscrim[,i], na.rm = T), tTest$p.value)  
      sigDiff_silking_hours_temp_hum <- rbind(sigDiff_silking_hours_temp_hum, row)
    } 
    silking_hoursTempHum_compare_allTraits <- rbind(silking_hoursTempHum_compare_allTraits, silking_hours_temp_hum_compare)
    
    ## Planting to silking
    planting_to_silking_hours_temp_hum <- subset(planting_to_silking_hours_temp_hum_summary, planting_to_silking_hours_temp_hum_summary$days_before_after_silking == d)
    
    planting_to_silking_hours_temp_hum_mostDiscrim <- planting_to_silking_hours_temp_hum %>%
      filter(experiment %in% mostDiscrim)
    
    planting_to_silking_hours_temp_hum_leastDiscrim <- planting_to_silking_hours_temp_hum %>%
      filter(experiment %in% leastDiscrim)
    
    planting_to_silking_hours_temp_hum_mostDiscrim$category <- "Most discriminating"
    planting_to_silking_hours_temp_hum_leastDiscrim$category <- "Least discriminating"
    
    planting_to_silking_hours_temp_hum_compare <- rbind(planting_to_silking_hours_temp_hum_mostDiscrim, planting_to_silking_hours_temp_hum_leastDiscrim)
    planting_to_silking_hours_temp_hum_compare$trait <- trait
    
    for (i in 3:16) {
      # i = 3
      
      tTest <- t.test(planting_to_silking_hours_temp_hum_mostDiscrim[,i], planting_to_silking_hours_temp_hum_leastDiscrim[,i])
      
      row <- c(d, trait, colnames(planting_to_silking_hours_temp_hum_compare)[i],
               mean(planting_to_silking_hours_temp_hum_mostDiscrim[,i], na.rm = T), 
               mean(planting_to_silking_hours_temp_hum_leastDiscrim[,i], na.rm = T) , 
               tTest$p.value)
      sigDiff_planting_to_silking_hours_temp_hum <- rbind(sigDiff_planting_to_silking_hours_temp_hum, row)
    } 
    planting_to_silking_hoursTempHum_compare_allTraits <- rbind(planting_to_silking_hoursTempHum_compare_allTraits, planting_to_silking_hours_temp_hum_compare)
  }
}

silking_hoursTempHum_compare_allTraits <- silking_hoursTempHum_compare_allTraits[,colSums(silking_hoursTempHum_compare_allTraits != 0) > 0]
planting_to_silking_hoursTempHum_compare_allTraits <- planting_to_silking_hoursTempHum_compare_allTraits[,colSums(planting_to_silking_hoursTempHum_compare_allTraits != 0) > 0]

sigDiff_silking_hours_temp_hum <- as.data.frame(sigDiff_silking_hours_temp_hum)
sigDiff_silking_hours_temp_hum$pvalue <- as.numeric.factor(sigDiff_silking_hours_temp_hum$pvalue)
sig_silking <- sigDiff_silking_hours_temp_hum %>%
  filter(pvalue < 0.05)

sigDiff_planting_to_silking_hours_temp_hum <- as.data.frame(sigDiff_planting_to_silking_hours_temp_hum)
sigDiff_planting_to_silking_hours_temp_hum$pvalue <- as.numeric.factor(sigDiff_planting_to_silking_hours_temp_hum$pvalue)
sig_planting_to_silking <- sigDiff_planting_to_silking_hours_temp_hum %>%
  filter(pvalue < 0.05)

#### Silking
silking_hoursTempHum_Tg30 <- ggplot(data = silking_hoursTempHum_compare_allTraits, aes(x = factor(category), y = Tg30, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 60,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = silking_hoursTempHum_Tg30, filename = paste("./most_vs_least_discrim/silking_hoursTempHum_Tg30.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")


silking_hoursTempHum_Tg25 <- ggplot(data = silking_hoursTempHum_compare_allTraits, aes(x = factor(category), y = Tg25, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 60,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = silking_hoursTempHum_Tg25, filename = paste("./most_vs_least_discrim/silking_hoursTempHum_Tg25.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")

silking_hoursTempHum_Tg22 <- ggplot(data = silking_hoursTempHum_compare_allTraits, aes(x = factor(category), y = Tg22, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 60,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = silking_hoursTempHum_Tg22, filename = paste("./most_vs_least_discrim/silking_hoursTempHum_Tg22.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")

silking_hoursTempHum_Tl22 <- ggplot(data = silking_hoursTempHum_compare_allTraits, aes(x = factor(category), y = Tl22, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 60,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = silking_hoursTempHum_Tl22, filename = paste("./most_vs_least_discrim/silking_hoursTempHum_Tl22.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")


#### Planting to silking
planting_to_silking_hoursTempHum_Tg30 <- ggplot(data = planting_to_silking_hoursTempHum_compare_allTraits, aes(x = factor(category), y = Tg30, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 60,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = planting_to_silking_hoursTempHum_Tg30, filename = paste("./most_vs_least_discrim/planting_to_silking_hoursTempHum_Tg30.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")


planting_to_silking_hoursTempHum_Tg25 <- ggplot(data = planting_to_silking_hoursTempHum_compare_allTraits, aes(x = factor(category), y = Tg25, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 60,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = planting_to_silking_hoursTempHum_Tg25, filename = paste("./most_vs_least_discrim/planting_to_silking_hoursTempHum_Tg25.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")

planting_to_silking_hoursTempHum_Tg22 <- ggplot(data = planting_to_silking_hoursTempHum_compare_allTraits, aes(x = factor(category), y = Tg22, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 60,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = planting_to_silking_hoursTempHum_Tg22, filename = paste("./most_vs_least_discrim/planting_to_silking_hoursTempHum_Tg22.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")

planting_to_silking_hoursTempHum_Tl22 <- ggplot(data = planting_to_silking_hoursTempHum_compare_allTraits, aes(x = factor(category), y = Tl22, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 60,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ days_before_after_silking, scales = "free_x", switch = "y")
ggsave(plot = planting_to_silking_hoursTempHum_Tl22, filename = paste("./most_vs_least_discrim/planting_to_silking_hoursTempHum_Tl22.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")

#### Plots of latitude, longitude, and planting density ####
lat_long_plDens <- read.csv("latitude_longitude_plantingDensity.csv", header = T)

lat_long_plDens$Longitude <- abs(lat_long_plDens$Longitude)
colnames(lat_long_plDens) <- c("Experiment", "Location", "Latitude", "-Longitude", "Planting density")

for (i in 3:4) {
  # i = 5
  plot <- 
    ggplot(lat_long_plDens, mapping = aes(x = reorder(lat_long_plDens[,"Experiment"], lat_long_plDens[,i], FUN = median_), y = lat_long_plDens[,i])) +
    geom_point() +
    labs(x = "Environment", y = colnames(lat_long_plDens)[i]) +
    theme_bw() +
    theme(text = element_text(family = "Courier"),
          axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 70, hjust = 1), 
          axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.line = element_line(color = "black"))
  ggsave(plot = plot, filename = paste0("scatterplot_", colnames(lat_long_plDens)[i], ".tif"), device = "tiff",
         height = 4, width = 8, units = "in")
}

plot_plantingDensity <-
  ggplot(lat_long_plDens, mapping = aes(x = reorder(lat_long_plDens[, "Experiment"], lat_long_plDens[,5], FUN = median_), y = lat_long_plDens[,5])) +
  geom_point() +
  labs(x = "Environment", y = "Planting density\n(plants per ha)") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 70, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(plot = plot_plantingDensity, filename = paste0("scatterplot_", colnames(lat_long_plDens)[5], ".tif"), device = "tiff",
       height = 4, width = 8, units = "in")

#### Compare latitude, longitude, and planting density in most/least discriminating enviros ####
lat_long_plDens <- read.csv("latitude_longitude_plantingDensity.csv", header = T)
lat_long_plDens <- lat_long_plDens[,-2]

colnames(lat_long_plDens) <- c("experiment", "Latitude", "Longitude", "Planting density")

sigDiff_latLongPlDens <- matrix(NA, nrow = 0, ncol = 5)
colnames(sigDiff_latLongPlDens) <- c("trait", "weather", "mean_mostDiscrim", "mean_leastDiscrim", "pvalue")

latLongPlDens_compare_allTraits <- matrix(NA, nrow = 0, ncol = 6)
colnames(latLongPlDens_compare_allTraits) <- c("experiment", "Latitude", "-Longitude", "Planting density",
                                               "category", "trait")

for (trait in traits) {
  # trait = "Anthesis (GDU)"
  
  traitDiscrim <- discrimValue_allTraits[,c("enviros",trait)]
  
  mostDiscrim <- traitDiscrim %>%
    top_n(5) %>%
    .$enviros
  
  leastDiscrim <- traitDiscrim %>%
    top_n(-5) %>%
    .$enviros
  
  latLongPlDens_mostDiscrim <- lat_long_plDens %>%
    filter(experiment %in% mostDiscrim)
  
  latLongPlDens_leastDiscrim <- lat_long_plDens %>%
    filter(experiment %in% leastDiscrim)
  
  latLongPlDens_mostDiscrim$category <- "Most discriminating"
  latLongPlDens_leastDiscrim$category <- "Least discriminating"
  
  latLongPlDens_compare <- rbind(latLongPlDens_mostDiscrim, latLongPlDens_leastDiscrim)
  latLongPlDens_compare$trait <- trait
  
  for (i in 2:4) {
    i = 2
    
    # pvalue <- my.t.test.p.value(silking_humidity_mostDiscrim[,i], silking_humidity_leastDiscrim[,i])
    tTest <- t.test(latLongPlDens_mostDiscrim[,i], latLongPlDens_leastDiscrim[,i])
    
    row <- c(trait, colnames(latLongPlDens_compare)[i], mean(latLongPlDens_mostDiscrim[,i], na.rm = T), mean(latLongPlDens_leastDiscrim[,i], na.rm = T), tTest$p.value)  
    sigDiff_latLongPlDens <- rbind(sigDiff_latLongPlDens, row)
  } 
  latLongPlDens_compare_allTraits <- rbind(latLongPlDens_compare_allTraits, latLongPlDens_compare)
} 

latLongPlDens_compare_allTraits_melt <- melt(latLongPlDens_compare_allTraits,
                                             id.vars = c("trait", "category", "experiment"))

latitude <- latLongPlDens_compare_allTraits_melt %>%
  filter(variable == "Latitude")

latitude_plot <- ggplot(data = latitude, aes(x = factor(category), y = value, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.y = 25,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ variable, switch = "y")
ggsave(plot = latitude_plot, filename = paste("./most_vs_least_discrim/latitude.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")

longitude <- latLongPlDens_compare_allTraits_melt %>%
  filter(variable == "Longitude")

longitude_plot <- ggplot(data = longitude, aes(x = factor(category), y = value, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     label.x = 2,
                     label.y = -100,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ variable, switch = "y")
ggsave(plot = longitude_plot, filename = paste("./most_vs_least_discrim/longitude.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")

plantingDensity <- latLongPlDens_compare_allTraits_melt %>%
  filter(variable == "Planting density")

plantingDensity_plot <- ggplot(data = plantingDensity, aes(x = factor(category), y = value, group = category)) +
  geom_boxplot(aes(fill = category)) +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  stat_compare_means(method = "t.test",
                     
                     label.y = 25,
                     aes(label = paste0("p = ", ..p.format..))) +
  ggTheme_panels +
  coord_flip() +
  facet_grid(trait ~ variable, switch = "y")
ggsave(plot = plantingDensity_plot, filename = paste("./most_vs_least_discrim/plantingDensity.tif",sep = ""),
       device = "tiff", width = 6.5, height = 6.5, units = "in")

#### For each trait, rank correlation of environmental index, discriminability, ####
#### latitude, longitude, planting density, precipitation, and temperature factors ####
# dir.create("./Correlation_EI_discrim_weather_location")

for (i in 5:ncol(selectedPhenoData_standardized)) {
  # i = 5
  means_by_exp <- aggregate(selectedPhenoData_standardized, 
                            by = list(selectedPhenoData_standardized$Enviro), FUN = mean, na.rm = T)
  means_by_exp$Enviro <- means_by_exp$Group.1
  means_by_exp <- means_by_exp[,-1]
}

means_by_exp <- means_by_exp[,-c(1,3,4)]
hours_temp_summary <- hours_temp_hum_summary[,1:5]

colnames(lat_long_plDens)[1] <- "Enviro"
colnames(rainSummary)[1] <- "Enviro"
colnames(tempSummary)[1] <- "Enviro"
colnames(hours_temp_summary)[1] <- "Enviro"

for (i in 2:ncol(means_by_exp)) {
  # i = 2
  
  enviroIndex <- means_by_exp[,c(1,i)]
  colnames(enviroIndex) <- c("Enviro", paste0("EI_", colnames(means_by_exp)[i]))
  
  enviroDiscrim <- discrimValue_allTraits[,c(1,i)]
  colnames(enviroDiscrim) <- c("Enviro", paste0("discrim_", colnames(discrimValue_allTraits)[i]))
  
  enviroCharacteristics <- left_join(enviroIndex, enviroDiscrim)
  enviroCharacteristics <- left_join(enviroCharacteristics, lat_long_plDens)
  enviroCharacteristics <- left_join(enviroCharacteristics, rainSummary)
  enviroCharacteristics <- left_join(enviroCharacteristics, tempSummary)
  enviroCharacteristics <- left_join(enviroCharacteristics, hours_temp_summary)
  
  
  cormat <- cor(enviroCharacteristics[,-1], use = "p", method = "spearman")
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
  p.mat <- cor.mtest(enviroCharacteristics[,-1])
  
  ggcorrplot <- 
    ggcorrplot(cormat, p.mat = p.mat, type = "upper", lab = T,
               colors = c("darkred", "white", "steelblue"),
               ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
    labs(x = NULL, y = NULL) +
    theme(text = element_text(family = "Courier"),
          axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
          axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 8),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.line = element_line(color = "black"))
  ggsave(ggcorrplot, filename = paste0("./Correlation_EI_discrim_weather_location/", colnames(means_by_exp)[i],"_ggCorrPlot-EI,discrim,weather,location.tif"), device = "tiff")
}

#### Plot number of traits with sig diff weather patterns for most least discrim enviros ####

# !!!!! BE SURE TO UPDATE THE CSV FILE BELOW INSIDE EXCEL ####
silking_weather_least_most_discrim_enviros <- read.csv("most_vs_least_discrim/silking_weather_least_most_discrim_enviros.csv")
silking_weather_least_most_discrim_enviros$Trait <- as.character(silking_weather_least_most_discrim_enviros$Trait)
silking_weather_least_most_discrim_enviros$Trait <- gsub(pattern = "Plot weight", "Plot grain weight", silking_weather_least_most_discrim_enviros$Trait)
silking_weather_least_most_discrim_enviros$Trait <- factor(silking_weather_least_most_discrim_enviros$Trait,
                                                           levels = c("Anthesis", "Silking",
                                                                      "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                                                      "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                                                      "Kernel length", "Kernel width", "Kernel thickness"))

silking_weather_least_most_discrim_enviros$Period <- as.character(silking_weather_least_most_discrim_enviros$Period)
silking_weather_least_most_discrim_enviros$Period <- gsub("Planting_to_Silking", "Planting to silking", silking_weather_least_most_discrim_enviros$Period)
silking_weather_least_most_discrim_enviros$Period <- factor(silking_weather_least_most_discrim_enviros$Period,
                                                            levels = c("Season", "Planting to silking", "Silking"))

silking_weather_least_most_discrim_enviros$Weather <- gsub("Tl22", "hrs temp < 22 °C", silking_weather_least_most_discrim_enviros$Weather)
silking_weather_least_most_discrim_enviros$Weather <- gsub("Tg22", "hrs temp > 22 °C", silking_weather_least_most_discrim_enviros$Weather)
silking_weather_least_most_discrim_enviros$Weather <- gsub("Tg25", "hrs temp > 25 °C", silking_weather_least_most_discrim_enviros$Weather)
silking_weather_least_most_discrim_enviros$Weather <- gsub("Tg30", "hrs temp > 30 °C", silking_weather_least_most_discrim_enviros$Weather)
silking_weather_least_most_discrim_enviros$Weather <- gsub("min temp", "minimum temp", silking_weather_least_most_discrim_enviros$Weather)
silking_weather_least_most_discrim_enviros$Weather <- gsub("max temp", "maximum temp", silking_weather_least_most_discrim_enviros$Weather)

silking_weather_least_most_discrim_enviros$Weather <- as.character(silking_weather_least_most_discrim_enviros$Weather)
silking_weather_least_most_discrim_enviros$Weather <- factor(silking_weather_least_most_discrim_enviros$Weather,
                                                             levels = unique(silking_weather_least_most_discrim_enviros$Weather))

count_traits_with_diff_discrim_weather <-
  ggplot(silking_weather_least_most_discrim_enviros, aes(x = Weather, y = Value, fill = Trait)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Trait",
                    values = c("#54278f", "#756bb1", 
                               "#006d2c", "#31a354", 
                               "#fecc5c",
                               "#a50f15", "#de2d26", "#fb6a4a", "#fcae91",
                               "#08306b", "#08519c", "#3182bd", "#6baed6", "#c6dbef"
                    )) +
  labs(y = "Number of traits") +
  geom_abline(slope = 0, intercept = 0) +
  annotate("text", label = "Most discriminating", x = 4, y = 7, size = 3, family = "Courier") +
  annotate("text", label = "Least discriminating", x = 4, y = -3, size = 3, family = "Courier") +
  facet_grid(Period~., switch = "y") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 10, angle = 53, hjust = 1), 
        axis.title.y = element_text(size = 0), axis.text.y = element_text(size = 0),
        axis.ticks.y = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(color = "white", fill = "white"),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))
ggsave(plot = count_traits_with_diff_discrim_weather, filename = "stacked_barplot.tif", device = "tiff",
       width = 4.5, height = 8, units = "in")

silking_weather_least_most_discrim_enviros_no_seasonal <- silking_weather_least_most_discrim_enviros %>%
  filter(., Period != "Season")

count_traits_with_diff_discrim_weather_no_seasonal <-
  ggplot(silking_weather_least_most_discrim_enviros_no_seasonal, aes(x = Weather, y = Value, fill = Trait)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Trait",
                    values = c("#54278f", "#756bb1", 
                               "#006d2c", "#31a354", 
                               "#fecc5c",
                               "#a50f15", "#de2d26", "#fb6a4a", "#fcae91",
                               "#08306b", "#08519c", "#3182bd", "#6baed6", "#c6dbef"
                    )) +
  labs(y = "Number of traits") +
  geom_abline(slope = 0, intercept = 0) +
  annotate("text", label = "Most discriminating", x = 4, y = 7, size = 3, family = "Courier") +
  annotate("text", label = "Least discriminating", x = 4, y = -3, size = 3, family = "Courier") +
  facet_grid(Period~., switch = "y") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 10, angle = 53, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 0),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        strip.background = element_rect(color = "white", fill = "white"),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
ggsave(plot = count_traits_with_diff_discrim_weather_no_seasonal, filename = "stacked_barplot_no_seasonal.tif", device = "tiff",
       width = 4.5, height = 6, units = "in")

