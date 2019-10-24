##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: Load, combine, and clean data
## Date: 2018-09-18
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)

#### Read in, format, and organize EAR IMAGING data for 2014 and 2015 ####

#### 2014
# earData_2014 <- read.csv("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/Ear Photometry Phenotypic Information/gxe2014_PhenoAvgs.csv", header = T)
earData_2014 <- read.csv("./gxe2014_PhenoAvgs.csv", header = T)
earData_2014$Loc = as.factor(earData_2014$Loc)
earData_2014$Genotype = as.factor(earData_2014$Genotype)

earData_2014 <- as.data.frame(earData_2014)
earData_2014$Year <- "14"
earData_2014$Enviro <- paste(earData_2014$Loc, earData_2014$Year, sep = "_")

earData_2014 <- earData_2014[,c(1:4, 19:20, 5:6, 9, 7:8, 10:13, 15, 14)]
colnames(earData_2014) <- c("Plot", "Loc", "Genotype", "Rep", "Year", "Enviro",
                            "Plot Grain Weight", "Cup Weight", "Kernel Width", "Kernel Length",
                            "Kernel Area", "Kernel Weight", "KRN", "Kernel Thickness",
                            "Kernels Per Row", "Ear Length", "Ear Width")

earData_2014 <- earData_2014[-which(earData_2014$Enviro == "MNI2_14" & earData_2014$Rep == "2"),] # MN2_14, keep only reps 1 and 3 (correct planting date)
earData_2014 <- earData_2014[-which(earData_2014$Enviro == "MNI2_14" & earData_2014$Rep > 3),] # MN2_14, keep only reps 1 and 3 (correct planting date)

earData_2014 <- earData_2014[-which(earData_2014$Enviro == "PAI1_14" & earData_2014$Rep == "1"),] # PA1_14, keep only reps 2 and 4 (normal N)
earData_2014 <- earData_2014[-which(earData_2014$Enviro == "PAI1_14" & earData_2014$Rep == "3"),] # PA1_14, keep only reps 2 and 4 (normal N)

#### 2015
# earData_2015 <- read.csv("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/Ear Photometry Phenotypic Information/gxe2015_PhenoAvgs.csv")
earData_2015 <- read.csv("./gxe2015_PhenoAvgs.csv", header = T)
earData_2015$Plot = as.factor(earData_2015$Plot)
earData_2015$Loc = as.factor(earData_2015$Loc)
earData_2015$Genotype = as.factor(earData_2015$Genotype)
earData_2015$Rep = as.factor(earData_2015$Rep)

earData_2015 <- as.data.frame(earData_2015)
earData_2015$Year <- "15"
earData_2015$Enviro <- paste(earData_2015$Loc, earData_2015$Year, sep = "_")

earData_2015 <- earData_2015[,c(1:4, 19:20, 5:6, 9, 7:8, 10:13, 15, 14)]
colnames(earData_2015) <- c("Plot", "Loc", "Genotype", "Rep", "Year", "Enviro",
                            "Plot Grain Weight", "Cup Weight", "Kernel Width", "Kernel Length",
                            "Kernel Area", "Kernel Weight", "KRN", "Kernel Thickness",
                            "Kernels Per Row", "Ear Length", "Ear Width")

earData_2015 <- earData_2015[-which(earData_2015$Loc == "PAI2"),] # get rid of PAI2 (whole exp is low N)

#### Combine 2014 & 2015
earData_20142015 <- rbind(earData_2014,earData_2015)

earData_20142015[which(earData_20142015$Genotype == "W117HT"), "Genotype"] <- "W117"
earData_20142015[which(earData_20142015$Genotype == "W117Ht"), "Genotype"] <- "W117"

earData_20142015$Enviro <- gsub("WI1", "WII1", earData_20142015$Enviro)
earData_20142015$Enviro <- gsub("WI2", "WII2", earData_20142015$Enviro)

earData_20142015$Enviro <- gsub(pattern = "I1_", replacement = "1_", x = earData_20142015$Enviro)
earData_20142015$Enviro <- gsub(pattern = "I2_", replacement = "2_", x = earData_20142015$Enviro)
earData_20142015$Enviro <- gsub(pattern = "I3_", replacement = "3_", x = earData_20142015$Enviro)
earData_20142015$Enviro <- gsub(pattern = "I4_", replacement = "4_", x = earData_20142015$Enviro)

earData_20142015$Plot = as.factor(earData_20142015$Plot)
earData_20142015$Enviro = as.factor(earData_20142015$Enviro)
earData_20142015$Genotype = as.factor(earData_20142015$Genotype)
earData_20142015$Rep = as.factor(earData_20142015$Rep)

write.csv(earData_20142015, "gxe20142015earImagingData.csv", row.names = F)
rm(earData_2014, earData_2015)

#### Read in, format, and organize AGRONOMIC data for 2014 and 2015 ####

#### 2014
agronomicData_2014 <- read.csv("g2f_2014_inbred_data_20161206.csv", header = T)

#### Remove extra MN and PA reps (and other location extra reps)
agronomicData_2014 <- agronomicData_2014[-which(agronomicData_2014$Location == "GXE_inb_MN2" & agronomicData_2014$Rep > 3),]
agronomicData_2014 <- agronomicData_2014[-which(agronomicData_2014$Location == "GXE_inb_MN2" & agronomicData_2014$Rep == 2),]
agronomicData_2014 <- agronomicData_2014[-which(agronomicData_2014$Location == "GXE_inb_PA1" & agronomicData_2014$Rep == 1),]
agronomicData_2014 <- agronomicData_2014[-which(agronomicData_2014$Location == "GXE_inb_PA1" & agronomicData_2014$Rep == 3),]

agronomicData_2014 <- agronomicData_2014[,c(14,2,9,12,7, 17:19,47:52)]
colnames(agronomicData_2014) <- c("Plot", "Location", "Rep", "Genotype", "Planting date", 
                                  "Stand count", "Stalk lodging", "Root lodging",
                                  "Anthesis date (DAP)", "Silking date (DAP)", "Plant height (cm)",
                                  "Ear height (cm)", "Root lodging (%)", "Stalk lodging (%)")
agronomicData_2014$`Root lodging (%)` <- (agronomicData_2014$`Root lodging`/agronomicData_2014$`Stand count`)*100

agronomicData_2014$`Stalk lodging (%)` <- (agronomicData_2014$`Stalk lodging`/agronomicData_2014$`Stand count`)*100
agronomicData_2014 <- agronomicData_2014[,c(1:5, 9:ncol(agronomicData_2014))]

agronomicData_2014$Location <- gsub("GXE_inb_", "", agronomicData_2014$Location)
agronomicData_2014$Location <- gsub("_MAD", "", agronomicData_2014$Location)
agronomicData_2014$Location <- gsub("_ARL", "", agronomicData_2014$Location)
agronomicData_2014$Year <- "14"
agronomicData_2014$Enviro <- paste(agronomicData_2014$Location, agronomicData_2014$Year, sep = "_")
agronomicData_2014$`Planting date YYYY-MM-DD` <- as.Date(agronomicData_2014$`Planting date`, format = "%m/%d/%Y")
agronomicData_2014$`Planting date YYYY-MM-DD` <- gsub(pattern = "00", replacement = "20", x = agronomicData_2014$`Planting date YYYY-MM-DD`)
agronomicData_2014$`Planting date (DoY)` <- as.numeric(strftime(agronomicData_2014$`Planting date YYYY-MM-DD`, format = "%j"))
agronomicData_2014$`Anthesis date (DoY)` <- agronomicData_2014$`Planting date (DoY)` + agronomicData_2014$`Anthesis date (DAP)`
agronomicData_2014$`Silking date (DoY)` <- agronomicData_2014$`Planting date (DoY)` + agronomicData_2014$`Silking date (DAP)`

#### 2015
agronomicData_2015 <- read.csv("g2f_2015_inbred_Celeste.csv", header = T)

agronomicData_2015 <- agronomicData_2015[,c(1, 3, 9:10, 14, 18:24)]
colnames(agronomicData_2015) <- c("Plot", "Location", "Genotype", "Rep", "Planting date",
                                  "Stand count","Pollen (DAP)", "Silk (DAP)", "Plant height (cm)", "Ear height (cm)",
                                  "Root lodging", "Stalk lodging")
agronomicData_2015$rootLodgingPercent <- (agronomicData_2015$`Root lodging`/agronomicData_2015$`Stand count`)*100
agronomicData_2015$stalkLodgingPercent <- (agronomicData_2015$`Stalk lodging`/agronomicData_2015$`Stand count`)*100
agronomicData_2015 <- agronomicData_2015[,c(1:2,4,3,5,7:10,13:14)]
colnames(agronomicData_2015) <- c("Plot", "Location", "Rep", "Genotype","Planting date",
                                  "Anthesis date (DAP)", "Silking date (DAP)", "Plant height (cm)", "Ear height (cm)",
                                  "Root lodging (%)", "Stalk lodging (%)")
agronomicData_2015$Year <- "15"
agronomicData_2015$Enviro <- paste(agronomicData_2015$Location, agronomicData_2015$Year, sep = "_")

agronomicData_2015 <- agronomicData_2015[-which(agronomicData_2015$Location == "PAI2"),] # get rid of PAI2 (whole exp is low N)

agronomicData_2015$`Planting date YYYY-MM-DD` <- as.Date(agronomicData_2015$`Planting date`, format = "%m/%d/%Y")
agronomicData_2015$`Planting date YYYY-MM-DD` <- gsub(pattern = "00", replacement = "20", x = agronomicData_2015$`Planting date YYYY-MM-DD`)
agronomicData_2015$`Planting date (DoY)` <- as.numeric(strftime(agronomicData_2015$`Planting date YYYY-MM-DD`, format = "%j"))
agronomicData_2015$`Anthesis date (DoY)` <- agronomicData_2015$`Planting date (DoY)` + agronomicData_2015$`Anthesis date (DAP)`
agronomicData_2015$`Silking date (DoY)` <- agronomicData_2015$`Planting date (DoY)` + agronomicData_2015$`Silking date (DAP)`

#### Convert flowering data to GDUs ####
GDU_2014 <- read.csv("./GDU_2014.csv")
GDU_2015 <- read.csv("./GDU_2015.csv")

#### 2014
GDU_2014 <- subset(GDU_2014, location == "TXH1  TXI1  TXI2_2014" | location == "TXH2  TXI3_2014" |
                     location == "GAI2_2014" | location == "PAI1_2014" |
                     location == "MOH1 MOI1_2014" | location == "MOH2 MOI2 MOI3_2014" |
                     location == "WII2_2014" | location == "WIH1 WII1_2014" |
                     location == "MNI2_2014" | location == "ILH1 ILI1_2014" |
                     location == "NEH1 NEI1_2014" | location == "INH1 INI1_2014" |
                     location == "NCI1_2014" | location == "DEI1_2014" |
                     location == "IAH1 IAI1_2014" | location == "IAI3_2014" |
                     location == "IAI2_2014" | location == "NYH1 NYI1_2014")
GDU_2014$location <- gsub("WIH1 WII1_2014", "WII1_2014", x = GDU_2014$location)
GDU_2014$location <- gsub("TXH2  TXI3_2014", "TXI3_2014", x = GDU_2014$location)
GDU_2014$location <- gsub("MOH1 MOI1_2014", "MOI1_2014", x = GDU_2014$location)
GDU_2014$location <- gsub("ILH1 ILI1_2014", "ILI1_2014", x = GDU_2014$location)
GDU_2014$location <- gsub("NEH1 NEI1_2014", "NEI1_2014", x = GDU_2014$location)
GDU_2014$location <- gsub("INH1 INI1_2014", "INI1_2014", x = GDU_2014$location)
GDU_2014$location <- gsub("NYH1 NYI1_2014", "NYI1_2014", x = GDU_2014$location)
GDU_2014$location <- gsub("IAH1 IAI1_2014", "IAI1_2014", x = GDU_2014$location)
GDU_2014$location <- gsub("TXH1  TXI1  TXI2_2014", "TXI1_2014  TXI2_2014", x = GDU_2014$location)
GDU_2014$location <- gsub("MOH2 MOI2 MOI3_2014", "MOI2_2014 MOI3_2014", x = GDU_2014$location)

### Fix TXI1 and TXI2
GDU_2014_TXI1_2014 <- subset(GDU_2014, location == "TXI1_2014  TXI2_2014")
GDU_2014_TXI1_2014$location <- gsub("TXI1_2014  TXI2_2014", "TXI1_2014", x = GDU_2014_TXI1_2014$location)
GDU_2014$location <- gsub("TXI1_2014  TXI2_2014", "TXI2_2014", x = GDU_2014$location)
GDU_2014 <- rbind(GDU_2014, GDU_2014_TXI1_2014)

### Fix MOI2 and MOI3
GDU_2014_MOI2_2014 <- subset(GDU_2014, location == "MOI2_2014 MOI3_2014")
GDU_2014_MOI2_2014$location <- gsub("MOI2_2014 MOI3_2014", "MOI2_2014", x = GDU_2014_MOI2_2014$location)
GDU_2014$location <- gsub("MOI2_2014 MOI3_2014", "MOI3_2014", x = GDU_2014$location)
GDU_2014 <- rbind(GDU_2014, GDU_2014_MOI2_2014)

rm(GDU_2014_MOI2_2014, GDU_2014_TXI1_2014)

GDU_2014$location <- gsub("I1_20", "1_", GDU_2014$location)
GDU_2014$location <- gsub("I2_20", "2_", GDU_2014$location)
GDU_2014$location <- gsub("I3_20", "3_", GDU_2014$location)
GDU_2014$location <- gsub("I4_20", "4_", GDU_2014$location)

GDU_2014 <- GDU_2014[with(GDU_2014, order(location, dayOfYear)),]

rownames(GDU_2014) <- 1:nrow(GDU_2014)

agronomicData_2014$`Anthesis (GDU)` <- NA
agronomicData_2014$`Silking (GDU)` <- NA

locations <- unique(agronomicData_2014$Enviro)
agronomicData_2014_withGDU <- data.frame()

for (location in locations) {
  # location = "TX1_14"
  GDU_location <- GDU_2014[which(GDU_2014$location == location),]
  data_location <- agronomicData_2014[which(agronomicData_2014$Enviro == location),]
  
  rownames(GDU_location) <- 1:nrow(GDU_location)
  
  for (i in 1:nrow(data_location)) {
    # i = 1
    plantingDoY <- data_location[i,15]
    anthesisDoY <- data_location[i,16]
    silkingDoY <- data_location[i,17]
    
    plantingDoY_rowNum <- rownames(GDU_location[which(GDU_location$dayOfYear == plantingDoY),])
    anthesisDoY_rowNum <- rownames(GDU_location[which(GDU_location$dayOfYear == anthesisDoY),])
    silkingDoY_rowNum <- rownames(GDU_location[which(GDU_location$dayOfYear == silkingDoY),])
    
    if (is.na(anthesisDoY)) {
      data_location[i,18] <- NA
    } else if (length(plantingDoY_rowNum) == 0) {
      data_location[i,18] <- NA
    } else {
      data_location[i,18] <- sum(GDU_location[plantingDoY_rowNum:anthesisDoY_rowNum, 5], na.rm = T)
    }
    
    
    if (is.na(silkingDoY)) {
      data_location[i,19] <- NA
    } else if (length(plantingDoY_rowNum) == 0) {
      data_location[i,19] <- NA
    } else {
      data_location[i,19] <- sum(GDU_location[plantingDoY_rowNum:silkingDoY_rowNum, 5], na.rm = T)
    }
  }
  agronomicData_2014_withGDU <- rbind(agronomicData_2014_withGDU, data_location)
}


#### 2015
GDU_2015 <- subset(GDU_2015, location == "KSH1  KSI1_2015" | location == "IAI2_2015" |
                     location == "TXH1  TXI1  TXI2_2015" | location == "GAI1_2015" |
                     location == "AZI1  AZI2_2015" | location == "PAI1  PAI2_2015" |
                     location == "MOH1  MOI1  MOH2  MOI2_2015" | location == "WIH2  WII2_2015" |
                     location == "WIH1  WII1_2015" | location == "MNI1_2015" |
                     location == "ILH1  ILI1  ILH2_2015" | location == "INH1  INI1_2015" |
                     location == "NCI1_2015" | location == "DEI1_2015" |
                     location == "IAI1_2015" | location == "IAI3_2015" |
                     location == "NYH3  NYI2_2015" | location == "NYH1  NYI1_2015" |
                     location == "IAI4_2015")


GDU_2015$location <- gsub("KSH1  KSI1_2015", "KSI1_2015", x = GDU_2015$location)
GDU_2015$location <- gsub("TXH1  TXI1  TXI2_2015", "TXI1  TXI2_2015", x = GDU_2015$location)
GDU_2015$location <- gsub("MOH1  MOI1  MOH2  MOI2_2015", "MOI1  MOI2_2015", x = GDU_2015$location)
GDU_2015$location <- gsub("WIH2  WII2_2015", "WII2_2015", x = GDU_2015$location)
GDU_2015$location <- gsub("WIH1  WII1_2015", "WII1_2015", x = GDU_2015$location)
GDU_2015$location <- gsub("KSH1  KSI1_2015", "KSI1_2015", x = GDU_2015$location)
GDU_2015$location <- gsub("ILH1  ILI1  ILH2_2015", "ILI1_2015", x = GDU_2015$location)
GDU_2015$location <- gsub("INH1  INI1_2015", "INI1_2015", x = GDU_2015$location)
GDU_2015$location <- gsub("NYH3  NYI2_2015", "NYI2_2015", x = GDU_2015$location)
GDU_2015$location <- gsub("NYH1  NYI1_2015", "NYI1_2015", x = GDU_2015$location)


### Fix AZI1 and AZI2
GDU_2015_AZI1_2015 <- subset(GDU_2015, location == "AZI1  AZI2_2015")
GDU_2015_AZI1_2015$location <- gsub("AZI1  AZI2_2015", "AZI1_2015", x = GDU_2015_AZI1_2015$location)
GDU_2015$location <- gsub("AZI1  AZI2_2015", "AZI2_2015", x = GDU_2015$location)
GDU_2015 <- rbind(GDU_2015, GDU_2015_AZI1_2015)

### Fix TXI1 and TXI2
GDU_2015_TXI1_2015 <- subset(GDU_2015, location == "TXI1  TXI2_2015")
GDU_2015_TXI1_2015$location <- gsub("TXI1  TXI2_2015", "TXI1_2015", x = GDU_2015_TXI1_2015$location)
GDU_2015$location <- gsub("TXI1  TXI2_2015", "TXI2_2015", x = GDU_2015$location)
GDU_2015 <- rbind(GDU_2015, GDU_2015_TXI1_2015)

### Fix PAI1 and PAI2
GDU_2015_PAI1_2015 <- subset(GDU_2015, location == "PAI1  PAI2_2015")
GDU_2015_PAI1_2015$location <- gsub("PAI1  PAI2_2015", "PAI1_2015", x = GDU_2015_PAI1_2015$location)
GDU_2015$location <- gsub("PAI1  PAI2_2015", "PAI2_2015", x = GDU_2015$location)
GDU_2015 <- rbind(GDU_2015, GDU_2015_PAI1_2015)

### Fix MOI1 and MOI2
GDU_2015_MOI1_2015 <- subset(GDU_2015, location == "MOI1  MOI2_2015")
GDU_2015_MOI1_2015$location <- gsub("MOI1  MOI2_2015", "MOI1_2015", x = GDU_2015_MOI1_2015$location)
GDU_2015$location <- gsub("MOI1  MOI2_2015", "MOI2_2015", x = GDU_2015$location)
GDU_2015 <- rbind(GDU_2015, GDU_2015_MOI1_2015)

rm(GDU_2015_AZI1_2015, GDU_2015_MOI1_2015, GDU_2015_PAI1_2015, GDU_2015_TXI1_2015)

GDU_2015$location <- gsub("I1_20", "I1_", GDU_2015$location)
GDU_2015$location <- gsub("I2_20", "I2_", GDU_2015$location)
GDU_2015$location <- gsub("I3_20", "I3_", GDU_2015$location)
GDU_2015$location <- gsub("I4_20", "I4_", GDU_2015$location)

GDU_2015 <- GDU_2015[with(GDU_2015, order(location, dayOfYear)),]

rownames(GDU_2015) <- 1:nrow(GDU_2015)

agronomicData_2015$`Anthesis (GDU)` <- NA
agronomicData_2015$`Silking (GDU)` <- NA


locations <- unique(agronomicData_2015$Enviro)
agronomicData_2015_withGDU <- data.frame()

#### remove SDI1_15 because missing weather data (therefore no GDUs)
locations <- locations[-which(locations == "SDI1_15")]

for (location in locations) {
  # location = "TXI1_15"
  GDU_location <- GDU_2015[which(GDU_2015$location == location),]
  data_location <- agronomicData_2015[which(agronomicData_2015$Enviro == location),]
  
  rownames(GDU_location) <- 1:nrow(GDU_location)
  
  for (i in 1:nrow(data_location)) {
    # i = 1
    plantingDoY <- data_location[i,15]
    anthesisDoY <- data_location[i,16]
    silkingDoY <- data_location[i,17]
    
    plantingDoY_rowNum <- rownames(GDU_location[which(GDU_location$dayOfYear == plantingDoY),])
    anthesisDoY_rowNum <- rownames(GDU_location[which(GDU_location$dayOfYear == anthesisDoY),])
    silkingDoY_rowNum <- rownames(GDU_location[which(GDU_location$dayOfYear == silkingDoY),])
    
    if (is.na(anthesisDoY)) {
      data_location[i,18] <- NA
    } else if (length(plantingDoY_rowNum) == 0) {
      data_location[i,18] <- NA
    } else {
      data_location[i,18] <- sum(GDU_location[plantingDoY_rowNum:anthesisDoY_rowNum, 5], na.rm = T)
    }
    
    if (is.na(silkingDoY)) {
      data_location[i,19] <- NA
    } else if (length(plantingDoY_rowNum) == 0) {
      data_location[i,19] <- NA
    } else {
      data_location[i,19] <- sum(GDU_location[plantingDoY_rowNum:silkingDoY_rowNum, 5], na.rm = T)
    }
  }
  agronomicData_2015_withGDU <- rbind(agronomicData_2015_withGDU, data_location)
}

#### Combine 2014 & 2015
agronomicData_20142015 <- rbind(agronomicData_2014_withGDU,agronomicData_2015_withGDU)
agronomicData_20142015$uniquePlot <- seq(from = 1, to = nrow(agronomicData_20142015))
agronomicData_20142015 <- agronomicData_20142015[,c(20, 1, 13, 3, 4, 14, 6:9, 15:19)]
colnames(agronomicData_20142015)[6] <- "Planting date"

agronomicData_20142015[which(agronomicData_20142015$Genotype == "W117HT"), "Genotype"] <- "W117"
agronomicData_20142015[which(agronomicData_20142015$Genotype == "W117Ht"), "Genotype"] <- "W117"

agronomicData_20142015$Enviro <- gsub("WI1", "WII1", agronomicData_20142015$Enviro)
agronomicData_20142015$Enviro <- gsub("WI2", "WII2", agronomicData_20142015$Enviro)

agronomicData_20142015$Enviro <- gsub(pattern = "I1_", replacement = "1_", x = agronomicData_20142015$Enviro)
agronomicData_20142015$Enviro <- gsub(pattern = "I2_", replacement = "2_", x = agronomicData_20142015$Enviro)
agronomicData_20142015$Enviro <- gsub(pattern = "I3_", replacement = "3_", x = agronomicData_20142015$Enviro)
agronomicData_20142015$Enviro <- gsub(pattern = "I4_", replacement = "4_", x = agronomicData_20142015$Enviro)

write.csv(agronomicData_20142015, "gxe20142015-agronomicData.csv", row.names = F)
rm(agronomicData_2014, agronomicData_2015, agronomicData_2014_withGDU, agronomicData_2015_withGDU)

#### Combine EAR IMAGING and AGRONOMIC data ####
allPhenoData <- merge(agronomicData_20142015, earData_20142015, by = "Plot", all = T)
selectedPhenoData <- allPhenoData[,c(1,3,20,4,18,5,17,6:8,14:15,9:10,21:31)]

for (i in 1:nrow(selectedPhenoData)) {
  # i=507
  if (is.na(selectedPhenoData[i,"Rep.x"])) {
    selectedPhenoData[i,"Rep.x"] = as.character(selectedPhenoData[i,"Rep.y"])
  }
}

for (i in 1:nrow(selectedPhenoData)) {
  # i=507
  if (is.na(selectedPhenoData[i,"Genotype.x"])) {
    selectedPhenoData[i,"Genotype.x"] = as.character(selectedPhenoData[i,"Genotype.y"])
  }
}

for (i in 1:nrow(selectedPhenoData)) {
  # i=507
  if (is.na(selectedPhenoData[i,"Enviro.x"])) {
    selectedPhenoData[i,"Enviro.x"] = as.character(selectedPhenoData[i,"Enviro.y"])
  }
}

selectedPhenoData <- selectedPhenoData[,-c(3,5,7)]

colnames(selectedPhenoData) <- c("Plot", "Enviro", "Rep", "Genotype", "Planting date", "Anthesis date (DAP)",
                                 "Silking date (DAP)", "Anthesis (GDU)", "Silking (GDU)", 
                                 "Plant height", "Ear height", "Plot grain weight", "Cup weight", "Kernel width", 
                                 "Kernel length", "Kernel area", "Kernel weight", "Kernel row number", 
                                 "Kernel thickness", "Kernels per row", "Ear length", "Ear width")
selectedPhenoData$Plot = as.factor(selectedPhenoData$Plot)
selectedPhenoData$Enviro = as.factor(selectedPhenoData$Enviro)
selectedPhenoData$Genotype = as.factor(selectedPhenoData$Genotype)
selectedPhenoData$Rep = as.factor(selectedPhenoData$Rep)

traits <- colnames(selectedPhenoData[,6:ncol(selectedPhenoData)])

rm(allPhenoData, earData_20142015, agronomicData_20142015)

#### Remove unwanted locations ####
### MO1_2014, MO2_2014, MO3_2014, TX3_2014, WI2_2014, PA1_2015
### TX1_2015 <-- get rid of plant height
### TX1_14, TX2_14, TX1_15, NY1_14, PA1_14, WI1_15 <-- get rid of anthesis and silking (GDU)

selectedPhenoData <- selectedPhenoData[-which(selectedPhenoData$Enviro == "MO1_14"),]
selectedPhenoData <- selectedPhenoData[-which(selectedPhenoData$Enviro == "MO2_14"),]
selectedPhenoData <- selectedPhenoData[-which(selectedPhenoData$Enviro == "MO3_14"),]
selectedPhenoData <- selectedPhenoData[-which(selectedPhenoData$Enviro == "TX3_14"),]
selectedPhenoData <- selectedPhenoData[-which(selectedPhenoData$Enviro == "WI2_14"),]
selectedPhenoData <- selectedPhenoData[-which(selectedPhenoData$Enviro == "AZ1_15"),]
selectedPhenoData <- selectedPhenoData[-which(selectedPhenoData$Enviro == "AZ2_15"),]
selectedPhenoData <- selectedPhenoData[-which(selectedPhenoData$Enviro == "PA1_15"),]


for (r in 1:nrow(selectedPhenoData)) {
  if (selectedPhenoData[r, 2] == "TX1_15") {
    selectedPhenoData[r,"Plant height"] <- NA
  }
}

for (r in 1:nrow(selectedPhenoData)) {
  if (selectedPhenoData[r,2] == "TX1_14" | 
      selectedPhenoData[r,2] == "TX2_14" |
      selectedPhenoData[r,2] == "TX1_15" |
      selectedPhenoData[r,2] == "NY1_14" |
      selectedPhenoData[r,2] == "PA1_14" |
      selectedPhenoData[r,2] == "WI1_15") {
    selectedPhenoData[r,"Anthesis (GDU)"] <- NA
    selectedPhenoData[r,"Silking (GDU)"] <- NA
  }
}

#### Remove PHG83 ####
selectedPhenoData <- selectedPhenoData[-which(selectedPhenoData$Genotype == "PHG83"),]

#### Remove outliers from AGRONOMIC data ####
NA_count_before <- data.frame(trait = colnames(selectedPhenoData)[8:11],
                              NA_count_before = c(sum(is.na(selectedPhenoData$`Anthesis (GDU)`)),
                                                  sum(is.na(selectedPhenoData$`Silking (GDU)`)),
                                                  sum(is.na(selectedPhenoData$`Plant height`)),
                                                  sum(is.na(selectedPhenoData$`Ear height`))))

for (i in 6:11) {
  # i = 8
  allMeas = c(selectedPhenoData[,i])
  allMeasPlots = c(selectedPhenoData[,1])
  
  hist(allMeas,breaks = 100, main = colnames(selectedPhenoData[i]))
  abline(v = ((mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T))),col = "red")
  abline(v = ((mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T))),col = "red")
  
  which(allMeas > (mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T)) | allMeas < (mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T)))
  plots = allMeasPlots[which(allMeas > (mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T)) | 
                               allMeas < (mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T)))]
  plots = unique(plots)
  
  probs = selectedPhenoData[which(rownames(selectedPhenoData) == plots[1]),c("Plot",colnames(selectedPhenoData[i]))]
  
  for (r in 2:length(plots)) {
    probs = rbind(probs,selectedPhenoData[which(selectedPhenoData[,1] == plots[r]),c("Plot",colnames(selectedPhenoData[i]))])
  }
  ((mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T)))
  ((mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T)))
  probs
  
  for (p in 1:nrow(probs)) {
    selectedPhenoData[which(selectedPhenoData[,1] == probs[p,1]), colnames(selectedPhenoData[i])] <- NA
  }
  
  ## check it again
  allMeas = c(selectedPhenoData[,i])
  allMeasPlots = c(selectedPhenoData[,1])
  
  hist(allMeas,breaks = 100, main = colnames(selectedPhenoData[i]))
  abline(v = ((mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T))),col = "red")
  abline(v = ((mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T))),col = "red")
  
  which(allMeas > (mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T)) | allMeas < (mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T)))
  plots = allMeasPlots[which(allMeas > (mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T)) | 
                               allMeas < (mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T)))]
  plots = unique(plots)
  
  probs = selectedPhenoData[which(rownames(selectedPhenoData) == plots[1]),c("Plot",colnames(selectedPhenoData[i]))]
  
  for (r in 2:length(plots)) {
    probs = rbind(probs,selectedPhenoData[which(selectedPhenoData[,1] == plots[r]),c("Plot",colnames(selectedPhenoData[i]))])
  }
  ((mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T)))
  ((mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T)))
  probs
  
  for (p in 1:nrow(probs)) {
    selectedPhenoData[which(selectedPhenoData[,1] == probs[p,1]), colnames(selectedPhenoData[i])] <- NA
    
    ## check it again
    allMeas = c(selectedPhenoData[,i])
    allMeasPlots = c(selectedPhenoData[,1])
    
    hist(allMeas,breaks = 100, main = colnames(selectedPhenoData[i]))
    abline(v = ((mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T))),col = "red")
    abline(v = ((mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T))),col = "red")
    
    which(allMeas > (mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T)) | allMeas < (mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T)))
    plots = allMeasPlots[which(allMeas > (mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T)) | 
                                 allMeas < (mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T)))]
    plots = unique(plots)
    
    probs = selectedPhenoData[which(rownames(selectedPhenoData) == plots[1]),c("Plot",colnames(selectedPhenoData[i]))]
    
    for (r in 2:length(plots)) {
      probs = rbind(probs,selectedPhenoData[which(selectedPhenoData[,1] == plots[r]),c("Plot",colnames(selectedPhenoData[i]))])
    }
    ((mean(allMeas,na.rm = T)) + (5*sd(allMeas,na.rm = T)))
    ((mean(allMeas,na.rm = T)) - (5*sd(allMeas,na.rm = T)))
    probs
    
    for (p in 1:nrow(probs)) {
      selectedPhenoData[which(selectedPhenoData[,1] == probs[p,1]), colnames(selectedPhenoData[i])] <- NA
    }
  }
}

rm(probs, allMeas, allMeasPlots, i, p, plots, r)

NA_count_after <- data.frame(trait = colnames(selectedPhenoData)[8:11],
                             NA_count_after = c(sum(is.na(selectedPhenoData$`Anthesis (GDU)`)),
                                                sum(is.na(selectedPhenoData$`Silking (GDU)`)),
                                                sum(is.na(selectedPhenoData$`Plant height`)),
                                                sum(is.na(selectedPhenoData$`Ear height`))))
outlier_count <- merge(NA_count_before, NA_count_after)
outlier_count$outliers <- outlier_count$NA_count_after - outlier_count$NA_count_before

#### Remove flowering dates as DAP and cup weight ####
selectedPhenoData <- selectedPhenoData[,-c(5:7, 13)]

selectedPhenoData <- selectedPhenoData[,c(1:9, 17:18, 16,14, 13,12,11,10,15)]

write.csv(selectedPhenoData, "selectedPhenoData.csv", row.names = F,quote = F)

