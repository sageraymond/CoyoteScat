#load libraries
library(lme4) ##
library(survival)
library(Hmisc)
library(ggplot2)
library(AICcmodavg)
library(rafalib)
library(MuMIn)
library(MASS)
library(pROC)
library(car)
library(caret)
library(psych)
library(dplyr)
library(tidyr)
library(broom)
library(tidyverse)
library(magrittr)
library(irr)
library(splitstackshape)
library(jtools)
library(ggstance)
library(ggh4x)
library(grid)
library(gridExtra)
library(jtools)
library(interactions)
library(BAMMtools)
library(DescTools)
library(AMR)
library(patchwork)
library(reshape2)
library(AER)
library(glmmTMB)
library(performance)
library(rsq)
library(glmmTMB)
library(ggpubr)
library(corrplot)
library(magick)
library(sjPlot)
#########################################################################

#Part 1 = Infection status among different habitat types

#########################################################################
###Load and summarise data
Deanna <- read.csv("Deanna_Scat2023.csv")
sum(Deanna$season == "winter") #130
sum(Deanna$season == "winter") / 269 #0.4832714

sum(Deanna$season == "spring") #105
sum(Deanna$season == "spring") / 269 #0.3903346

sum(Deanna$season == "summer") #16
sum(Deanna$season == "summer") / 269 #0.05947955

sum(Deanna$season == "autumn") #18
sum(Deanna$season == "autumn") / 269 #0.0669145

#Create new column identifying infection status
Deanna$Infected <- 0
Deanna$Infected <- ifelse(Deanna$em.qpcr == "NEG", 0, 1)

#Determine overall infection rate
sum(Deanna$Infected == "1") #70
sum(Deanna$Infected == "1") / 269 #0.260223

#Does infection status vary among seasons?
tblSeason <- table(Deanna$season, Deanna$Infected) 
SeasonChi <- chisq.test(tblSeason) #X-squared = 5.1519, df = 3, p-value = 0.161
SeasonChi$expected #25% of values < 5, so fisher's test is more appropriate
FTSeason <- fisher.test(tblSeason) 
FTSeason # p-value = 0.1298

#Does infection status vary among ecological types?
#Create subset of scats where ecological type was recorded
SiteTypeSubset <- Deanna %>%
  filter(ecological.type.primary != "")
SiteTypeSubset$InfectionStatus <- ifelse(SiteTypeSubset$Infected == "1",
                                         "Infected", "Not Infected")
SiteTypeSubset$ecological.type.primary <- ifelse(SiteTypeSubset$ecological.type.primary == "Forest",
                                                 "Forested", SiteTypeSubset$ecological.type.primary)

tblType <- table(SiteTypeSubset$ecological.type.primary, SiteTypeSubset$InfectionStatus)
EcoChi <- chisq.test(tblType) 
EcoChi$expected #Chi square is appropriate, not fisher's
EcoChi #X-squared = 14.282, df = 4, p-value = 0.006446

#get standardized residuals
EcoChi$stdres

#Determine how many times more than expected
19/9.710843 #1.956576

#Create a chart showing observed and expected values
#include label to indicate magnitude of standardized residual
Ecoobs <- as.data.frame(EcoChi$observed)
Ecoobs <- Ecoobs %>%
  dplyr::filter(Var2 == "Infected") %>%
  dplyr::select(-(Var2)) %>%
  dplyr::rename("Habitat" = Var1,
                "Count" = Freq)
Ecoobs$OorE <- "Observed"

Ecoexp <- as.data.frame(EcoChi$expected)
Ecoexp$OorE <- "Expected"
Ecoexp$Habitat <- c("Compost", "Developed", "Forested", "Maintained",
                        "Naturalized")
Ecoexp <- Ecoexp %>%
  dplyr::select(Infected, OorE, Habitat) %>%
  dplyr::rename("Count" = Infected)
rownames(Ecoexp) <- NULL

EcoStR <- as.data.frame(EcoChi$stdres)
EcoStR <- EcoStR %>%
  dplyr::filter(Var2 == "Infected") %>%
  dplyr::rename("StRes" = Freq,
                 "Habitat" = Var1) %>%
  dplyr::select(c(Habitat, StRes))

EcoChartData <- rbind(Ecoobs, Ecoexp)

EcoChartData <- left_join(EcoChartData, EcoStR, by = "Habitat")
EcoChartData %>% mutate_if(is.numeric, ~round(., 1))
EcoChartData$stRes2 <- c("3.7", "-0.6", "-1.3", "-1.0", "-0.5", NA, NA, NA, NA, NA)

Ecochart1 <- ggplot(data = EcoChartData, mapping = aes(x = Habitat, y = Count, 
                                                       fill = OorE)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge()) +
  geom_text(aes(label=stRes2), vjust=c(-0.2, -2.6, -5.8, -3.8, -2, 0, 0, 0, 0, 0),
            hjust = 0.55, color="black", size=4.5) +
  theme_bw() +
  scale_fill_manual(values = c("#004D40", "#4DB6AC")) +
  scale_x_discrete(limits = c("Compost", "Forested", "Developed", "Naturalized",
                              "Maintained")) +
  theme(axis.text.x = element_text(colour = "black", face = "plain"),
        axis.text.y = element_text(colour = "black", face = "plain"),
        legend.title = element_blank(),
        legend.position = c(0.75,0.90)) +
  labs(x="Habitat", y="Number of infected scats")

Ecochart2 <- Ecochart1 + theme_classic()

EcoChart3 <- Ecochart2 + theme(axis.text.x = element_text(colour = "black", 
                                                          face = "plain", size = 12),
                               axis.text.y = element_text(colour = "black", 
                                                          face = "plain", size = 12),
                               axis.title.y = element_text(colour = "black", 
                                                           face = "plain", size = 12),
                               axis.title.x = element_text(colour = "black", 
                                                           face = "plain", size = 12),
                               legend.title = element_blank(),
                               legend.position = c(0.9, 0.95),
                               legend.text = element_text(colour = "black", 
                                                          face = "plain", size = 12))


#Does infection status vary among scats that do and don't contain anthro food?
DeannaAnthro <- Deanna %>%
  dplyr::filter(SageAnthro != "Unk") #264 scats had this info


sum(DeannaAnthro$SageAnthro == "Yes") #142
sum(DeannaAnthro$SageAnthro == "Yes") / 264 #0.538

DeannaAnthro$Infected <- ifelse(DeannaAnthro$Infected == 1, "Infected", "Not Infected")
DeannaAnthro$SageAnthro <- ifelse(DeannaAnthro$SageAnthro == "Yes", "Anthro", "No Anthro")

tblAnthro <- table(DeannaAnthro$SageAnthro, DeannaAnthro$Infected) 
AnthroChi <- chisq.test(tblAnthro)
AnthroChi#X-squared (with Yates) = X-squared = X-squared = 8.5154, df = 1, p-value = 0.003521
AnthroChi$expected
AnthroChi$observed
(48 - 37.11364) / 37.11364 #0.293325

#How many times more common is infection when scat has anthro food?
48/ 37.1 #1.293801

#Create a chart showing observed and expected values
#include label to indicate magnitude of standardized residual
Anthroobs <- as.data.frame(AnthroChi$observed)
Anthroobs <- Anthroobs %>%
  dplyr::filter(Var2 == "Infected") %>%
  dplyr::select(-(Var2)) %>%
  dplyr::rename("Anthro" = Var1,
                "Count" = Freq)
Anthroobs$OorE <- "Observed"

Anthroexp <- as.data.frame(AnthroChi$expected)
Anthroexp$OorE <- "Expected"
Anthroexp$Anthro <- c("Anthro", "No Anthro")

Anthroexp <- Anthroexp %>%
  dplyr::select(Infected, OorE, Anthro) %>%
  dplyr::rename("Count" = Infected)
rownames(Anthroexp) <- NULL

AnthroStR <- as.data.frame(AnthroChi$stdres)
AnthroStR <- AnthroStR %>%
  dplyr::filter(Var2 == "Infected") %>%
  dplyr::rename("StRes" = Freq,
                "Anthro" = Var1) %>%
  dplyr::select(c(Anthro, StRes))

AnthroChartData <- rbind(Anthroobs, Anthroexp)

AnthroChartData <- left_join(AnthroChartData, AnthroStR, by = "Anthro")
AnthroChartData$stRes2 <- c("3.1", "-3.1", NA, NA)
AnthroChartData$Anthro <- ifelse(AnthroChartData$Anthro == "Anthro", 
                                 "Anthropogenic food", "No anthropogenic food")

AnthroChartData$Anthro <- ifelse(AnthroChartData$Anthro == "Anthropogenic food", 
                                 "Present", "Not detected")

Anthrochart1 <- ggplot(data = AnthroChartData, mapping = aes(x = Anthro, y = Count,
                                                             fill = OorE)) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge()) +
  geom_text(aes(label=stRes2), vjust = c(-0.5,-6.6,0,0), hjust = 0.55, color="black",
            size=4.5) +
  theme_classic() +
  scale_fill_manual(values = c("#004D40", "#4DB6AC")) +
  scale_x_discrete(limits = c("Present", "Not detected")) +
  theme(axis.text.x = element_text(colour = "black", size = 12, face = "plain"),
        axis.text.y = element_text(colour = "black", size = 12, face = "plain"),
        axis.title.x = element_text(colour = "black", size = 12, face = "plain"),
        axis.title.y = element_text(colour = "black", size = 12, face = "plain"),
        legend.title = element_blank(),
        legend.position = "none") +
  ylim(0,52) +
  labs(x="Anthropogenic content", y="Number of infected scats")

#determine infection intensity among positive scats in different habitats and contents
Pos <- Deanna %>%
  dplyr::select(ecological.type.primary, SageAnthro, Infected, em.cp) %>%
  na.omit(Pos) %>%
  dplyr::rename("Habitat" = ecological.type.primary)

Pos$InfInt <- 1 / Pos$em.cp *100

PosHab <- Pos %>%
  dplyr::filter(Habitat != "")

aggregate(PosHab$em.cp, list(PosHab$Habitat),FUN=mean)
aggregate(PosHab$em.cp, list(PosHab$Habitat),FUN=SD)

PosHab$Habitat <- ifelse(PosHab$Habitat == "Forest", "Forested", PosHab$Habitat)

#Use a GLM approach to asses influence of site type
PosHab1 <- PosHab
PosHab1$Habitat <- factor(PosHab1$Habitat, levels = c("Forest", "Compost", 
                                                      "Developed", "Maintained", 
                                                      "Naturalized"))

HabGLM <- glm(em.cp~Habitat, data = PosHab1)
summary(HabGLM, scale = TRUE, center = TRUE)
summ(HabGLM, confint = TRUE, digits = 4)
summ(HabGLM, exp = TRUE, digits = 4)
summ(HabGLM, exp = TRUE, confint = TRUE, digits = 4)

#Repeat for anthropogenic content
PosHab2 <- PosHab
PosHab2 <- PosHab2 %>%
  dplyr::filter(SageAnthro != "Unk")

AnthGLM <- glm(em.cp~SageAnthro, data = PosHab2)
summary(AnthGLM)
summ(AnthGLM, confint = TRUE, digits = 4)
summ(AnthGLM, exp = TRUE, digits = 4)
summ(AnthGLM, exp = TRUE, confint = TRUE, digits = 4)

#Make charts displaying infection intensity (among infected scats)
Habchart1 <- ggplot(data = PosHab, mapping = aes(x = Habitat, y = InfInt, fill = 
                                                   "black")) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = "#4DB6AC") +
  scale_x_discrete(limits = c("Maintained", "Compost", "Forested", "Developed", 
                              "Naturalized")) +
  theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(x="Habitat", y="Infection intensity (1 / number of PCR cycles\nneeded to identify infection) * 100")


PosAnth <- Pos %>%
  dplyr::filter(SageAnthro != "Unk")
t.test(PosAnth$em.cp ~ PosAnth$SageAnthro) #t = 1.478, df = 30.216, p-value = 0.1498
aggregate(PosAnth$em.cp, list(PosAnth$SageAnthro),FUN=mean)

PosAnth$SageAnthro <- ifelse(PosAnth$SageAnthro == "No", "Not detected", "Present")

AnthIntchart1 <- ggplot(data = PosAnth, mapping = aes(x = SageAnthro, y = InfInt, 
                                                      fill = "black")) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = "#4DB6AC") +
  scale_x_discrete(limits = c("Present", "Not detected")) +
    theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(x="Anthropogenic content", y="Infection intensity (1 / number of PCR cycles\nneeded to identify infection) * 100")

#Is there a difference in content among habitats?
DeannaAnthro2 <- DeannaAnthro %>%
  dplyr::filter(ecological.type.primary != "")
AnthHab <- table(DeannaAnthro2$SageAnthro, DeannaAnthro2$ecological.type.primary)
AnthHabCHi <- chisq.test(AnthHab)
AnthHabCHi #X-squared = 9.7602, df = 4, p-value = 0.04467
AnthHabCHi$stdres #anthro less comon in naturalized and dorested, and more common in compot, developed, and maintained

#Create composite chart with four panels re. infection prevalence and intensity
EcoChart4 <- EcoChart3 + ggtitle("(A) Prevalence by habitat")
AnthChart2 <- Anthrochart1 + ggtitle("(B) Prevalence by content")  
HabChart2 <- Habchart1  + ggtitle("(C) Intensity by habitat")
AnthIntchart2 <- AnthIntchart1 + ggtitle("(D) Intensity by content")

DeannaComp <- ggarrange(EcoChart4, AnthChart2, HabChart2, AnthIntchart2, nrow = 2, ncol = 2, widths = c(9,5))
ggsave("DeannaComp.png", DeannaComp, width = 9, height = 8, dpi = 700)
#########################################################################

#Part 2 = scat distribution on the landscape (unit of rep = scat)

#########################################################################
#Load data
Scat1 <- read.csv("AllScats_April52022.csv")
#remove urine
Scat1 <- Scat1 %>%
  dplyr::filter(Scat1$Pee != "Y")

#Summarise seasonal data
#Winter
sum(Scat1$Month == "December") + sum(Scat1$Month == "January") + sum(Scat1$Month == "February") #748
(sum(Scat1$Month == "December") + sum(Scat1$Month == "January") + sum(Scat1$Month == "February")) / 1263 #0.5922407

#Spring
sum(Scat1$Month == "March") + sum(Scat1$Month == "April") + sum(Scat1$Month == "May") #444
(sum(Scat1$Month == "March") + sum(Scat1$Month == "April") + sum(Scat1$Month == "May")) / 1263 #0.3515439

#Autumn
sum(Scat1$Month == "September") + sum(Scat1$Month == "October") + sum(Scat1$Month == "November") #71
(sum(Scat1$Month == "September") + sum(Scat1$Month == "October") + sum(Scat1$Month == "November")) / 1263 #0.05621536

#Summarise Size
sum(Scat1$Size == "DNR") + sum(Scat1$Size == "Unk") #366
1263-366 #we know size for 897 scats

sum(Scat1$Size == "Med") #510
sum(Scat1$Size == "Med") / 897 #56.9%

sum(Scat1$Size == "Large") #236
sum(Scat1$Size == "Large") / 897 #26.3

sum(Scat1$Size == "Small") #151
sum(Scat1$Size == "Small") / 897 #16.8

#Summarise freshness
sum(Scat1$Freshness == "DNR") + sum(Scat1$Freshness == "Unk") #102
1263-102 #we know size for 1161 scats

sum(Scat1$Freshness == "Vfresh") #105
sum(Scat1$Freshness == "Vfresh") / 1161 #9.0%

sum(Scat1$Freshness == "Fresh") #329
sum(Scat1$Freshness == "Fresh") / 1161 #28.3%

sum(Scat1$Freshness == "Old") #439
sum(Scat1$Freshness == "Old") / 1161 #37.8%

sum(Scat1$Freshness == "Vold") #288
sum(Scat1$Freshness == "Vold") / 1161 #24.8


###########################################################################
#Stitch together the many files created in ArcGIS
stitchU1 <- read.csv('UseScatJun14.csv')
stitchU1 <- stitchU1 %>%
  dplyr::filter(stitchU1$Pee != "Y") %>%
  dplyr::select(FID, DistWater, DistBldg, DistRoad,
                DistEdge, DistNatEdg, DistScat, nateddens2,
                densrd25, densed25, densbldg25, sitetypera,
                curvature1, curvatur_1, curvatur_2, curvatur_3,
                curvatur_4, curvatur_5, curvatur_6, slope_ju_1,
                aspectju_1)
stitchU1 <- stitchU1 %>%
  rename("DistNatEdge" = DistNatEdg) %>%
  rename("NatEdDens" = nateddens2) %>%
  rename("EdDens" = densed25) %>%
  rename("RoadDens" = densrd25) %>%
  rename("BldgDens" = densbldg25) %>%
  rename("LC" = sitetypera) %>%
  rename("Curve1" = curvature1) %>%
  rename("Curve2" = curvatur_1) %>%
  rename("Curve3" = curvatur_2) %>%
  rename("Curve4" = curvatur_3) %>%
  rename("Curve5" = curvatur_4) %>%
  rename("Curve6" = curvatur_5) %>%
  rename("Curve7" = curvatur_6) %>%
  rename("Slope" = slope_ju_1) %>%
  rename("Aspect" = aspectju_1)

stitchA1 <- read.csv('AvaScatJun14.csv')
stitchA1 <- stitchA1 %>%
  dplyr::select(FID, DistWater, DistBldg, DistRoad,
                DistEdge, DistNEdge, DistScat, natedden_1,
                densrd25_1, densed25_1, densbldg_1, sitetype_1,
                curvature_, curvature2, curvature3, curvature4,
                curvature5, curvature6, curvature7, slope_jun1,
                aspectjun1)
stitchA1 <- stitchA1 %>%
  rename("DistNatEdge" = DistNEdge) %>%
  rename("NatEdDens" = natedden_1) %>%
  rename("EdDens" = densed25_1) %>%
  rename("RoadDens" = densrd25_1) %>%
  rename("BldgDens" = densbldg_1) %>%
  rename("LC" = sitetype_1) %>%
  rename("Curve1" = curvature_) %>%
  rename("Curve2" = curvature2) %>%
  rename("Curve3" = curvature3) %>%
  rename("Curve4" = curvature4) %>%
  rename("Curve5" = curvature5) %>%
  rename("Curve6" = curvature6) %>%
  rename("Curve7" = curvature7) %>%
  rename("Slope" = slope_jun1) %>%
  rename("Aspect" = aspectjun1)

#add LAND COVER
LC_U <- read.csv("LC_Use.csv")
LC_U[is.na(LC_U)] <- 0

LC_A <- read.csv("LC_Ava.csv")
LC_A[is.na(LC_A)] <- 0

stitchU2 <- stitchU1 %>%
  left_join(LC_U, by = "FID")

stitchA2 <- stitchA1 %>%
  full_join(LC_A, by = "FID")

SC_U <- read.csv('ScatCount_U.csv')
SC_U <- SC_U %>%
  dplyr::filter(SC_U$Pee != "Y") %>%
  dplyr::select(Join_Count, FID) %>%
  rename("ScatCount" = Join_Count)

SC_A <- read.csv('ScatCount_A.csv')
SC_A <- SC_A %>%
  dplyr::select(Join_Count, FID) %>%
  rename("ScatCount" = Join_Count)

stitchU3 <- stitchU2 %>%
  full_join(SC_U, by = "FID")

stitchA3 <- stitchA2 %>%
  full_join(SC_A, by = "FID")

EdgeU <- read.csv("Edge_U.csv")

stitchU4 <- stitchU3 %>%
  left_join(EdgeU, by = "FID")

EdgeA <- read.csv("Edge_A.csv")
stitchA4 <- stitchA3 %>%
  full_join(EdgeA, by = "FID")

NatEdgeU <- read.csv("NatEdge_U.csv")
stitchU5 <- stitchU4 %>%
  left_join(NatEdgeU, by = "FID")

NatEdgeA <- read.csv("NatEdge_A.csv")
stitchA5 <- stitchA4 %>%
  full_join(NatEdgeA, by = "FID")

TrailJuncU <- read.csv("Trails_U.csv")
stitchU6 <- stitchU5 %>%
  left_join(TrailJuncU, by = "FID")

TrailJuncA <- read.csv("Trails_A.csv")
stitchA6 <- stitchA5 %>%
  left_join(TrailJuncA, by = "FID")

WaterU <- read.csv("Water_U.csv")
stitchU7 <- stitchU6 %>%
  left_join(WaterU, by = "FID")

WaterA <- read.csv("Water_A.csv")
stitchA7 <- stitchA6 %>%
  left_join(WaterA, by = "FID")

CampGardPGU <- read.csv("CampGardPG_U.csv")
stitchU8 <- stitchU7 %>%
  left_join(CampGardPGU, by = "FID")

CampGardPGA <- read.csv("CampGardPG_A.csv")
stitchA8 <- stitchA7 %>%
  left_join(CampGardPGA, by = "FID")

#replace NA with 0
stitchU8[is.na(stitchU8)] <- 0
stitchA8[is.na(stitchA8)] <- 0

#Add a 'use' column
stitchU8$Use <- "1"
stitchA8$Use <- "0"

#subtract one from scat count for use sites
stitchU8$ScatCount <- stitchU8$ScatCount - 1

#Now join to make master dataset
Data <- rbind(stitchA8, stitchU8)

#And add eastness and northness
Data <- Data %>%
  mutate(AspRad = Aspect*0.01745329251) %>%
  mutate(East = sin(AspRad)) %>%
  mutate(North = cos(AspRad))
#write.csv(Data, file = "MasterScatJul26.csv")

####Now read in your data file and start working with it
MScat <- read.csv('MasterScatJul26.csv')
MScat <- MScat %>%
  dplyr::select(-c(X, FID))
MScat$LC <- as.character(MScat$LC)
MScat$LC[is.na(MScat$LC)] <- "0"
MScat$Use <- factor(MScat$Use, levels = c(0,1))
str(MScat)

MScat <- MScat %>%
  dplyr::mutate(Curve1a = abs(Curve1))

#Multiply your densities by 1000
MScat$NatEdDens <- MScat$NatEdDens*1000
MScat$RoadDens <- MScat$RoadDens*1000
MScat$EdDens <- MScat$EdDens*1000
MScat$BldgDens <- MScat$BldgDens*1000

#retain an original copy of each variable to facilitat ecaluclating effects size later
MScat$DistWaterOG <- MScat$DistWater
MScat$DistBldgOG <- MScat$DistBldg
MScat$DistRoadOG <- MScat$DistRoad
MScat$DistEdgeOG <- MScat$DistEdge
MScat$DistNatEdgeOG <- MScat$DistNatEdge
MScat$DistScatOG <- MScat$DistScat

MScat$DistJuncOG <- MScat$DistJunc
MScat$DistTrailOG <- MScat$DistTrail
MScat$DistCampOG <- MScat$DistCamp
MScat$DistPGOG <- MScat$DistPG

MScat$Curve1OG <- MScat$Curve1
MScat$NATOG <- MScat$NAT

#Make all the variables scaled and centered
MScat$DistWater <- scale(MScat$DistWater, scale = TRUE, center = TRUE)
MScat$DistBldg <- scale(MScat$DistBldg, scale = TRUE, center = TRUE)
MScat$DistRoad <- scale(MScat$DistRoad, scale = TRUE, center = TRUE)
MScat$DistEdge <- scale(MScat$DistEdge, scale = TRUE, center = TRUE)
MScat$DistNatEdge <- scale(MScat$DistNatEdge, scale = TRUE, center = TRUE)
MScat$DistScat <- scale(MScat$DistScat, scale = TRUE, center = TRUE)

MScat$NatEdDens <- scale(MScat$NatEdDens, scale = TRUE, center = TRUE)
MScat$RoadDens <- scale(MScat$RoadDens, scale = TRUE, center = TRUE)
MScat$EdDens <- scale(MScat$EdDens, scale = TRUE, center = TRUE)
MScat$BldgDens <- scale(MScat$BldgDens, scale = TRUE, center = TRUE)

MScat$Curve1 <- scale(MScat$Curve1, scale = TRUE, center = TRUE)

MScat$Slope <- scale(MScat$Slope, scale = TRUE, center = TRUE)
MScat$East <- scale(MScat$East, scale = TRUE, center = TRUE)
MScat$Slope <- scale(MScat$Slope, scale = TRUE, center = TRUE)

MScat$ANTH <- scale(MScat$ANTH, scale = TRUE, center = TRUE)
MScat$GRASS <- scale(MScat$GRASS, scale = TRUE, center = TRUE)
MScat$NAT <- scale(MScat$NAT, scale = TRUE, center = TRUE)
MScat$WATER <- scale(MScat$WATER, scale = TRUE, center = TRUE)

MScat$DistJunc <- scale(MScat$DistJunc, scale = TRUE, center = TRUE)
MScat$DistTrail <- scale(MScat$DistTrail, scale = TRUE, center = TRUE)
MScat$DistCamp <- scale(MScat$DistCamp, scale = TRUE, center = TRUE)
MScat$DistPG <- scale(MScat$DistPG, scale = TRUE, center = TRUE)

#change to numeric
MScat$DistWater <- as.numeric(MScat$DistWater)
MScat$DistBldg <- as.numeric(MScat$DistBldg)
MScat$DistRoad <- as.numeric(MScat$DistRoad)
MScat$DistEdge <- as.numeric(MScat$DistEdge)
MScat$DistNatEdge <- as.numeric(MScat$DistNatEdge)
MScat$DistScat <- as.numeric(MScat$DistScat)

MScat$NatEdDens <- as.numeric(MScat$NatEdDens)
MScat$RoadDens <- as.numeric(MScat$RoadDens)
MScat$EdDens <- as.numeric(MScat$EdDens)
MScat$BldgDens <- as.numeric(MScat$BldgDens)

MScat$Curve1 <- as.numeric(MScat$Curve1)

MScat$Slope <- as.numeric(MScat$Slope)
MScat$East <- as.numeric(MScat$East)
MScat$Slope <- as.numeric(MScat$Slope)

MScat$ANTH <- as.numeric(MScat$ANTH)
MScat$GRASS <- as.numeric(MScat$GRASS)
MScat$NAT <- as.numeric(MScat$NAT)
MScat$WATER <- as.numeric(MScat$WATER)

MScat$DistJunc <- as.numeric(MScat$DistJunc)
MScat$DistTrail <- as.numeric(MScat$DistTrail)
MScat$DistCamp <- as.numeric(MScat$DistCamp)
MScat$DistPG <- as.numeric(MScat$DistPG)


#Create univariate GLMs
mod1 <- glm(Use ~ DistWater, family = binomial, data = MScat)
mod2 <- glm(Use ~ DistBldg, family = binomial, data = MScat)
mod3 <- glm(Use ~ DistRoad, family = binomial, data = MScat) #NOT SIG
mod4 <- glm(Use ~ DistEdge, family = binomial, data = MScat)
mod5 <- glm(Use ~ DistNatEdge, family = binomial, data = MScat)
mod6 <- glm(Use ~ DistScat, family = binomial, data = MScat)

mod7 <- glm(Use ~ NatEdDens, family = binomial, data = MScat)
mod8 <- glm(Use ~ RoadDens, family = binomial, data = MScat)
mod9 <- glm(Use ~ EdDens, family = binomial, data = MScat)
mod10 <- glm(Use ~ BldgDens, family = binomial, data = MScat)

mod11 <- glm(Use ~ Curve1, family = binomial, data = MScat)

mod12 <- glm(Use ~ Slope, family = binomial, data = MScat)

mod13 <- glm(Use ~ ANTH, family = binomial, data = MScat)
mod14 <- glm(Use ~ GRASS, family = binomial, data = MScat)
mod15 <- glm(Use ~ NAT, family = binomial, data = MScat)

mod16 <- glm(Use ~ East, family = binomial, data = MScat)

mod17 <- glm(Use ~ DistJunc, family = binomial, data = MScat)
mod18 <- glm(Use ~ DistTrail, family = binomial, data = MScat)

mod19 <- glm(Use ~ Curve1a, family = binomial, data = MScat)

mod20 <- glm(Use ~ DistCamp, family = binomial, data = MScat)
mod21 <- glm(Use ~ DistPG, family = binomial, data = MScat)

summary(mod1) #MORE SCATS closer to water
summary(mod2) #more scats closer to buildings
summary(mod3) #NOT SIG (but more scats closer to roads)
summary(mod4) #More scats closer to edges
summary(mod5) #more scats closer to natural edges
summary(mod6) #more scats closer to other scats
summary(mod7) #NOT SIG of natural edge density
summary(mod8) #more scats at lower road densities
summary(mod9) #more scats at lower edge densities
summary(mod10) #NOT SIG (but more scats at lower byuilding density)
summary(mod11) #curv1 = sig  ####GO WITH CURVE 1 for now
summary(mod12) #more scats on steeper slopes
summary(mod13) #less anthro = more scat
summary(mod14) #less grass = more scat
summary(mod15) #more nat = more scat
summary(mod16) #more east = more scat
summary(mod17) #Dist Junction is not signficant
summary(mod18) #Dist Trail is not significant either
summary(mod19) #negative values are happier, which suggests flatter spots
summary(mod20) #dist to camp is not sig (p = 0.18)
summary(mod21) # dist to pg is not sig (p = 0.09)

#Decay terms
#Our distance terms are: water, bldg, road, natedge, scatdistance, Junc, Trail
mod1a <- glm(Use ~ exp(-0.002*DistWaterOG), family = binomial, data = MScat)
mod2a <- glm(Use ~ exp(-0.002*DistBldgOG), family = binomial, data = MScat)
mod3a <- glm(Use ~ exp(-0.002*DistRoadOG), family = binomial, data = MScat)
mod5a <- glm(Use ~ exp(-0.002*DistNatEdgeOG), family = binomial, data = MScat)
mod6a <- glm(Use ~ exp(-0.002*DistScatOG), family = binomial, data = MScat)
mod30a <- glm(Use ~ exp(-0.002*DistJuncOG), family = binomial, data = MScat)
mod31a <- glm(Use ~ exp(-0.002*DistTrailOG), family = binomial, data = MScat)
mod35a <- glm(Use ~ exp(-0.002*DistCampOG), family = binomial, data = MScat)
mod36a <- glm(Use ~ exp(-0.002*DistPGOG), family = binomial, data = MScat)

mod1b <- glm(Use ~ exp(-0.004*DistWaterOG), family = binomial, data = MScat)
mod2b <- glm(Use ~ exp(-0.004*DistBldgOG), family = binomial, data = MScat)
mod3b <- glm(Use ~ exp(-0.004*DistRoadOG), family = binomial, data = MScat)
mod5b <- glm(Use ~ exp(-0.004*DistNatEdgeOG), family = binomial, data = MScat)
mod6b <- glm(Use ~ exp(-0.004*DistScatOG), family = binomial, data = MScat)
mod30b <- glm(Use ~ exp(-0.004*DistJuncOG), family = binomial, data = MScat)
mod31b <- glm(Use ~ exp(-0.004*DistTrailOG), family = binomial, data = MScat)
mod35b <- glm(Use ~ exp(-0.004*DistCampOG), family = binomial, data = MScat)
mod36b <- glm(Use ~ exp(-0.004*DistPGOG), family = binomial, data = MScat)

mod1c <- glm(Use ~ exp(-0.006*DistWaterOG), family = binomial, data = MScat)
mod2c <- glm(Use ~ exp(-0.006*DistBldgOG), family = binomial, data = MScat)
mod3c <- glm(Use ~ exp(-0.006*DistRoadOG), family = binomial, data = MScat)
mod5c <- glm(Use ~ exp(-0.006*DistNatEdgeOG), family = binomial, data = MScat)
mod6c <- glm(Use ~ exp(-0.006*DistScatOG), family = binomial, data = MScat)
mod30c <- glm(Use ~ exp(-0.006*DistJuncOG), family = binomial, data = MScat)
mod31c <- glm(Use ~ exp(-0.006*DistTrailOG), family = binomial, data = MScat)
mod35c <- glm(Use ~ exp(-0.006*DistCampOG), family = binomial, data = MScat)
mod36c <- glm(Use ~ exp(-0.006*DistPGOG), family = binomial, data = MScat)

mod1d <- glm(Use ~ exp(-0.012*DistWaterOG), family = binomial, data = MScat)
mod2d <- glm(Use ~ exp(-0.012*DistBldgOG), family = binomial, data = MScat)
mod3d <- glm(Use ~ exp(-0.012*DistRoadOG), family = binomial, data = MScat)
mod5d <- glm(Use ~ exp(-0.012*DistNatEdgeOG), family = binomial, data = MScat)
mod6d <- glm(Use ~ exp(-0.012*DistScatOG), family = binomial, data = MScat)
mod30d <- glm(Use ~ exp(-0.012*DistJuncOG), family = binomial, data = MScat)
mod31d <- glm(Use ~ exp(-0.012*DistTrailOG), family = binomial, data = MScat)
mod35d <- glm(Use ~ exp(-0.012*DistCampOG), family = binomial, data = MScat)
mod36d <- glm(Use ~ exp(-0.012*DistPGOG), family = binomial, data = MScat)

mod1e <- glm(Use ~ exp(-0.03*DistWaterOG), family = binomial, data = MScat)
mod2e <- glm(Use ~ exp(-0.03*DistBldgOG), family = binomial, data = MScat)
mod3e <- glm(Use ~ exp(-0.03*DistRoadOG), family = binomial, data = MScat)
mod5e <- glm(Use ~ exp(-0.03*DistNatEdgeOG), family = binomial, data = MScat)
mod6e <- glm(Use ~ exp(-0.03*DistScatOG), family = binomial, data = MScat)
mod30e <- glm(Use ~ exp(-0.03*DistJuncOG), family = binomial, data = MScat)
mod31e <- glm(Use ~ exp(-0.03*DistTrailOG), family = binomial, data = MScat)
mod35e <- glm(Use ~ exp(-0.03*DistCampOG), family = binomial, data = MScat)
mod36e <- glm(Use ~ exp(-0.03*DistPGOG), family = binomial, data = MScat)

mod1f <- glm(Use ~ exp(-0.06*DistWaterOG), family = binomial, data = MScat)
mod2f <- glm(Use ~ exp(-0.06*DistBldgOG), family = binomial, data = MScat)
mod3f <- glm(Use ~ exp(-0.06*DistRoadOG), family = binomial, data = MScat)
mod5f <- glm(Use ~ exp(-0.06*DistNatEdgeOG), family = binomial, data = MScat)
mod6f <- glm(Use ~ exp(-0.06*DistScatOG), family = binomial, data = MScat)
mod30f <- glm(Use ~ exp(-0.06*DistJuncOG), family = binomial, data = MScat)
mod31f <- glm(Use ~ exp(-0.06*DistTrailOG), family = binomial, data = MScat)
mod35f <- glm(Use ~ exp(-0.06*DistCampOG), family = binomial, data = MScat)
mod36f <- glm(Use ~ exp(-0.06*DistPGOG), family = binomial, data = MScat)

mod1g <- glm(Use ~ exp(-0.2*DistWaterOG), family = binomial, data = MScat)
mod2g <- glm(Use ~ exp(-0.2*DistBldgOG), family = binomial, data = MScat)
mod3g <- glm(Use ~ exp(-0.2*DistRoadOG), family = binomial, data = MScat)
mod5g <- glm(Use ~ exp(-0.2*DistNatEdgeOG), family = binomial, data = MScat)
mod6g <- glm(Use ~ exp(-0.2*DistScatOG), family = binomial, data = MScat)
mod30g <- glm(Use ~ exp(-0.2*DistJuncOG), family = binomial, data = MScat)
mod31g <- glm(Use ~ exp(-0.2*DistTrailOG), family = binomial, data = MScat)
mod35g <- glm(Use ~ exp(-0.2*DistCampOG), family = binomial, data = MScat)
mod36g <- glm(Use ~ exp(-0.2*DistPGOG), family = binomial, data = MScat)

AIC(mod1); AIC(mod1a); AIC(mod1b); AIC(mod1c); AIC(mod1d); AIC(mod1e); AIC(mod1f); AIC(mod1g)
#water is best with no decay

AIC(mod2); AIC(mod2a); AIC(mod2b); AIC(mod2c); AIC(mod2d); AIC(mod2e); AIC(mod2f); AIC(mod2g)
#building is best without a decay

AIC(mod3); AIC(mod3a); AIC(mod3b); AIC(mod3c); AIC(mod3d); AIC(mod3e); AIC(mod3f); AIC(mod3g)
#dist road best with alpha = 0.06

AIC(mod5); AIC(mod5a); AIC(mod5b); AIC(mod5c); AIC(mod5d); AIC(mod5e); AIC(mod5f); AIC(mod5g)
#dist nat edge best with alpha = 0.002

AIC(mod6); AIC(mod6a); AIC(mod6b); AIC(mod6c); AIC(mod6d); AIC(mod6e); AIC(mod6f); AIC(mod6g)
#dist scat best with alpha = 0.06

AIC(mod30); AIC(mod30a); AIC(mod30b); AIC(mod30c); AIC(mod30d); AIC(mod30e); AIC(mod30f); AIC(mod30g)
#dist junction best with alpha = 0.002

AIC(mod31); AIC(mod31a); AIC(mod31b); AIC(mod31c); AIC(mod31d); AIC(mod31e); AIC(mod31f); AIC(mod31g)
#dist trail best with alpha = 0.2

AIC(mod35); AIC(mod35a); AIC(mod35b); AIC(mod35c); AIC(mod35d); AIC(mod35e); AIC(mod35f); AIC(mod35g)
#dist camp best at alpha =  0.012

AIC(mod36); AIC(mod36a); AIC(mod36b); AIC(mod36c); AIC(mod36d); AIC(mod36e); AIC(mod36f); AIC(mod36g)
#dist to pg best without decay term

#Add decay terms to your dataframe
MScat <- MScat %>%
  mutate(DecRoad = 1-(exp(-0.06*DistRoadOG))) %>%
  mutate(DecNatEdge = 1-(exp(-0.002*DistNatEdgeOG))) %>%
  mutate(DecScat = 1-(exp(-0.06*DistScatOG))) %>%
  mutate(DecJunc = 1-(exp(-0.002*DistJuncOG))) %>%
  mutate(DecTrail = 1-(exp(-0.2*DistTrailOG)))%>%
  mutate(DecCamp = 1-(exp(-0.012*DistCampOG)))

#Save originals
MScat$DecRoadOG <- MScat$DecRoad
MScat$DecNatEdgeOG <- MScat$DecNatEdge
MScat$DecScatOG <- MScat$DecScat
MScat$DecJuncOG <- MScat$DecJunc
MScat$DecTrailOG <- MScat$DecTrail
MScat$DecCampOG <- MScat$DecCamp

#Scale and center 
MScat$DecRoad <- scale(MScat$DecRoad, scale = TRUE, center = TRUE)
MScat$DecNatEdge <- scale(MScat$DecNatEdge, scale = TRUE, center = TRUE)
MScat$DecScat <- scale(MScat$DecScat, scale = TRUE, center = TRUE)
MScat$DecJunc <- scale(MScat$DecJunc, scale = TRUE, center = TRUE)
MScat$DecTrail <- scale(MScat$DecTrail, scale = TRUE, center = TRUE)
MScat$DecCamp <- scale(MScat$DecCamp, scale = TRUE, center = TRUE)

MScat$DecRoad <- as.numeric(MScat$DecRoad)
MScat$DecNatEdge <- as.numeric(MScat$DecNatEdge)
MScat$DecScat <- as.numeric(MScat$DecScat)
MScat$DecJunc <- as.numeric(MScat$DecJunc)
MScat$DecTrail <- as.numeric(MScat$DecTrail)
MScat$DecCamp <- as.numeric(MScat$DecCamp)


mod38 <- glm(Use ~ DecRoad, family = binomial, data = MScat)
mod39 <- glm(Use ~ DecNatEdge, family = binomial, data = MScat)
mod40 <- glm(Use ~ DecScat, family = binomial, data = MScat)
mod41 <- glm(Use ~ DecJunc, family = binomial, data = MScat)
mod42 <- glm(Use ~ DecTrail, family = binomial, data = MScat)
mod44 <- glm(Use ~ DecCamp, family = binomial, data = MScat)

summary(mod38) #more scats further from roads
summary(mod39) #more scats closert to nat edges
summary(mod40) #more scats closer to ecxising scats
summary(mod41) #more scats closert to junctions
summary(mod42) #more scats further from trails
summary(mod44) #more scats closer to camps

#Now develop univariate quadratics and compare to linear
qod1 <- glm(Use ~ DistWater + I(DistWater^2), family = binomial, data = MScat)
qod2 <- glm(Use ~ DistBldg + I(DistBldg^2), family = binomial, data = MScat)
qod3 <- glm(Use ~ DistRoad + (I(DistRoad^2)^2), family = binomial, data = MScat) #NOT SIG
qod4 <- glm(Use ~ DistEdge + I(DistEdge^2), family = binomial, data = MScat)
qod5 <- glm(Use ~ DistNatEdge + I(DistNatEdge^2), family = binomial, data = MScat)
qod6 <- glm(Use ~ DistScat + I(DistScat^2), family = binomial, data = MScat)

qod7 <- glm(Use ~ NatEdDens + I(NatEdDens^2), family = binomial, data = MScat)
qod8 <- glm(Use ~ RoadDens + I(RoadDens^2), family = binomial, data = MScat)
qod9 <- glm(Use ~ EdDens + I(EdDens^2), family = binomial, data = MScat)
qod10 <- glm(Use ~ BldgDens + I(BldgDens^2), family = binomial, data = MScat)

qod12 <- glm(Use ~ Curve1 + I(Curve1^2), family = binomial, data = MScat)

qod13 <- glm(Use ~ Slope + I(Slope^2), family = binomial, data = MScat)

qod14 <- glm(Use ~ ANTH + I(ANTH^2), family = binomial, data = MScat)
qod15 <- glm(Use ~ GRASS + I(GRASS^2), family = binomial, data = MScat)
qod16 <- glm(Use ~ NAT + I(NAT^2), family = binomial, data = MScat)

qod17 <- glm(Use ~ East + I(East^2), family = binomial, data = MScat)

qod18 <- glm(Use ~ DistJunc + I(DistJunc^2), family = binomial, data = MScat)
qod19 <- glm(Use ~ DistTrail + I(DistTrail^2), family = binomial, data = MScat)

qod20 <- glm(Use ~ DistCamp + I(DistCamp^2), family = binomial, data = MScat)
qod21 <- glm(Use ~ DistPG + I(DistPG^2), family = binomial, data = MScat)

qod38 <- glm(Use ~ DecRoad + I(DecRoad^2), family = binomial, data = MScat)
qod39 <- glm(Use ~ DecNatEdge + I(DecNatEdge^2), family = binomial, data = MScat)
qod40 <- glm(Use ~ DecScat + I(DecScat^2), family = binomial, data = MScat)
qod41 <- glm(Use ~ DecJunc + I(DecJunc^2), family = binomial, data = MScat)
qod42 <- glm(Use ~ DecTrail + I(DecTrail^2), family = binomial, data = MScat)
qod44 <- glm(Use ~ DecCamp + I(DecCamp^2), family = binomial, data = MScat)

AIC(mod1); AIC(qod1)
AIC(mod2); AIC(qod2)
AIC(mod3); AIC(qod3)
AIC(mod4); AIC(qod4)
AIC(mod5); AIC(qod5)
AIC(mod6); AIC(qod6)
AIC(mod7); AIC(qod7)
AIC(mod8); AIC(qod8)
AIC(mod9); AIC(qod9)
AIC(mod10); AIC(qod10)
AIC(mod11); AIC(qod11)
AIC(mod12); AIC(qod12)
AIC(mod13); AIC(qod13)
AIC(mod14); AIC(qod14)
AIC(mod15); AIC(qod15)
AIC(mod16); AIC(qod16)
AIC(mod17); AIC(qod17)
AIC(mod18); AIC(qod18)
AIC(mod19); AIC(qod19)
AIC(mod20); AIC(qod20)
AIC(mod21); AIC(qod21)
AIC(mod38); AIC(qod38)
AIC(mod39); AIC(qod39)
AIC(mod40); AIC(qod40)
AIC(mod41); AIC(qod41)
AIC(mod42); AIC(qod42)
AIC(mod43); AIC(qod43)


#final variable list is as follows:
#Dist water (quad) = qod1
#Dist building (lin) = mod2
#Dist road (quad) = qod3
#Dist Nat Edge (quad) = qod5
#Nat edge dens (quad) = qod7
#Road dens (lin) = mod8
#Curvature (lin) = mod11
#Slope (quad) = qod13
#ANthro (quad) = qod14
#Grass (line) = mod14
#Nat (quad) = qod16
#East (lin) = mod16
#Dist junction (quad) = qod18
#decay scat (quad) = qod40
#decay trail (quad) = qod42
#decay to camp (linear) = mod44

#Make a data frame with only retained variables plus 0/1!
ModDataSet <- MScat %>%
  dplyr::select(Use, DistWater, DistBldg, DistRoad, DistNatEdge,
                NatEdDens, RoadDens, Curve1, Slope, ANTH,
                GRASS, NAT, ScatCount, East, DistJunc, DecScat,
                DecTrail, DecCamp) %>% #and then make quad variables for quad terms
  mutate(DistRoadS = I(DistRoad^2)) %>%
  mutate(DistNatEdgeS = I(DistNatEdge^2)) %>%
  mutate(NatEdDensS = I(NatEdDens^2)) %>%
  mutate(SlopeS = I(Slope^2)) %>%
  mutate(ANTHS = I(ANTH^2)) %>%
  mutate(NATS = I(NAT^2)) %>%
  mutate(ScatCountS = I(ScatCount^2)) %>%
  mutate(DistJuncS = I(DistJunc^2)) %>%
  mutate(DecScatS = I(DecScat^2)) %>%
  mutate(DecTrailS = I(DecTrail^2))

#Summarise this suff for a table with univeriate models
export_summs(qod1, mod2, qod3, qod5, qod7, mod8, mod11, qod13, qod14, mod14, 
             qod16, mod16, qod18, qod40, qod42, mod44,
             ci_level = 0.95,
             stars = NULL,
             error_format = '(p = {p.value}), CI [{conf.low}, {conf.high}]',
             model.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"),
             to.file = "xlsx",
             digits = 4,
            file.name = "UnivMod ResultsAug24.xlsx")
            
#get pseudo r2 for all these models
PseudoR2(mod12, which = "Nagelkerke") #0.01418896 
PseudoR2(mod2, which = "Nagelkerke") #0.01077304 
PseudoR2(mod29, which = "Nagelkerke") #0.009098187 
PseudoR2(mod8, which = "Nagelkerke") #0.002667243 
PseudoR2(qod1, which = "Nagelkerke") #0.007273486 
PseudoR2(qod19, which = "Nagelkerke") #0.006016053 
PseudoR2(qod23, which = "Nagelkerke") #0.02161793 
PseudoR2(qod30, which = "Nagelkerke") #0.009183186 
PseudoR2(qod5, which = "Nagelkerke") #0.0371857 
PseudoR2(qod7, which = "Nagelkerke") #0.006270318 
PseudoR2(mod44, which = "Nagelkerke") #0.0187167 
PseudoR2(qod42, which = "Nagelkerke") #0.006186454 
PseudoR2(qod40, which = "Nagelkerke") #0.282947  


PseudoR2(mod12, which = "McFadden") #0.008559976 
PseudoR2(mod2, which = "McFadden") #0.006491662 
PseudoR2(mod29, which = "McFadden") #0.005479306 
PseudoR2(mod8, which = "McFadden") #0.001602829 
PseudoR2(qod1, which = "McFadden") #0.004377685 
PseudoR2(qod19, which = "McFadden") #0.003619333 
PseudoR2(qod23, which = "McFadden") #0.01307485 
PseudoR2(qod30, which = "McFadden") #0.005530656 
PseudoR2(qod5, which = "McFadden") #0.02261099 
PseudoR2(qod7, which = "McFadden") #0.003772627 
PseudoR2(mod44, which = "McFadden") #0.01130893
PseudoR2(qod42, which = "McFadden") #0.003722063 
PseudoR2(qod40, which = "McFadden") #0.1885238  


#Now, check for correlation.
#We will use r > 0.5 as a cutoff
cor(ModDataSet[sapply(ModDataSet, is.numeric)], method = c("pearson"))

#Correlation between:
#DistBldg and DistRd --> DIST BULDING retained
#Dist water and Dist small water--> WATER retained
#Slope and DIst Riv --> Dist River retained
#ANTH and NAT --> NAT retained
#NAT and GRASS --> NAT retained
#Scat Count and ScatDecay --> SCat decay retained
#DistJunc and DIstRiv --> Dist Junction retained


#Final variable list is now as follows:
#Dist water (quad)
#Dist building (lin)
#Dist Nat Edge (quad)
#Nat edge dens (quad)
#Road dens (lin)
#Curvature (lin)
#Nat (quad)
#East (lin)
#Dist junction (quad)
#decay scat (quad)
#decay trail (quad)
#decay to nearest camp (lin)
#dust to garden (quad)

#Make another (refined) data frame
ModDataSet1 <- MScat %>%
  dplyr::select(Use, DistWater, DistBldg, DistNatEdge,
                NatEdDens, RoadDens, Curve1,
                NAT, East, DistJunc, DecScat,
                DecTrail, Slope, DecCamp)

#Execute dredge
options(na.action = "na.fail")
Dredge3<-glm(Use~ DistWater + DistBldg + DistNatEdge + NatEdDens + RoadDens +
               Curve1 + NAT + East + DistJunc + DecScat + DecTrail + I(NAT^2) + 
               I(DecScat^2) + I(DecTrail^2) + I(DistNatEdge^2)+ I(NatEdDens^2) + 
               I(DistJunc^2)+ I(DistWater^2) + Slope + I(Slope^2) + DecCamp, 
             data = ModDataSet1, family = binomial)

#Create subset specifying that quadratic items must also show up as linear terms
subset1 <- expression(dc(`NAT`, `I(NAT^2)`) &
                        dc(`DecScat`, `I(DecScat^2)`) &
                        dc(`DecTrail`, `I(DecTrail^2)`) &
                        dc(`DistNatEdge`, `I(DistNatEdge^2)`) &
                        dc(`NatEdDens`, `I(NatEdDens^2)`) &
                        dc(`DistJunc`, `I(DistJunc^2)`)&
                        dc(`Slope`, `I(Slope^2)`)&
                        dc(`DistWater`, `I(DistWater^2)`))

Dredge5<-dredge(Dredge3, subset = subset1, trace = 2, rank = "BIC")
subset(Dredge5, delta < 2)

#Get models within 2 BIC
BICmod1a <- get.models(Dredge5, 1)[[1]]
BICmod2a <- get.models(Dredge5, 2)[[1]]
BICmod3a <- get.models(Dredge5, 3)[[1]]
BICmod4a <- get.models(Dredge5, 4)[[1]]
BICmod5a <- get.models(Dredge5, 5)[[1]]
BICmod6a <- get.models(Dredge5, 6)[[1]]
BICmod7a <- get.models(Dredge5, 7)[[1]]
BICmod8a <- get.models(Dredge5, 8)[[1]]

#Create models from call functions so you don't have to run dredge again
BICmod1a$call
BICmod1 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + 
                 DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)

BICmod2a$call
BICmod2 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + 
                 +                  DistBldg + DistNatEdge + NAT + 1, family = binomial, data = ModDataSet1)

BICmod3a$call
BICmod3 <- glm(formula = Use ~ Curve1 + DecScat + I(DecScat^2) + DistBldg + 
                 DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)

BICmod4a$call
BICmod4 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) +
                 DecTrail + DistBldg + DistNatEdge + NAT + 1, family = binomial, 
               data = ModDataSet1)

BICmod5a$call
BICmod5 <- glm(formula = Use ~ Curve1 + DecScat + I(DecScat^2) + DistBldg + 
                 DistNatEdge + NAT + 1, family = binomial, data = ModDataSet1)

BICmod6a$call
BICmod6 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + 
                 DecTrail + DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 
                 1, family = binomial, data = ModDataSet1)

BICmod7a$call
BICmod7 <- glm(formula = Use ~ Curve1 + DecScat + I(DecScat^2) + DecTrail + 
                 DistBldg + DistNatEdge + NAT + 1, family = binomial, data = ModDataSet1)

BICmod8a$call
BICmod8 <- glm(formula = Use ~ Curve1 + DecScat + I(DecScat^2) + DecTrail + 
                 DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)

#Create plot showing beta coefficients of all 8 main effects models
BICPlot1 <- plot_summs(BICmod1, BICmod2, BICmod3, BICmod4,
                       BICmod5, BICmod6, BICmod7, BICmod8,
                       colors = c("#00796B", "#00796B", "#00796B", "#00796B",
                                  "#00796B", "#00796B", "#00796B", "#00796B"),
                       point.size = 7,
                       point.shape = FALSE,
                       coefs = c("Distance to natural patch edge" = "DistNatEdge",
                                 "Decay distance to scat" = "DecScat",
                                 "Natural land cover (10-m)" = "NAT",
                                 "Decay distance to scat (Squared)" = "I(DecScat^2)",
                                 "Distance to building" = "DistBldg",
                                 "Topographic concavity" = "Curve1",
                                 "Distance to natural patch edge (Squared)" = "I(DistNatEdge^2)",
                                 "Decay distance to maintained trail" = "DecTrail",
                                 "Decay distance to camp" = "DecCamp"))
BICPlot1

BICPlot2 <- BICPlot1 + theme_classic()

#Improve Plot
BICPlot <- BICPlot2 + xlim(-1.5, 0.5) +
  theme(legend.position = "none",
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Coefficient estimate") +
  scale_shape_manual(values = c(24,24,24,24,24,24,24,24))

BICPlot

#Save file
ggsave("BICPlot1AVG.png", BICPlot, width = 7, height = 4, dpi = 700)


#Create a table containing this information for the SI
export_summs(BICmod1, BICmod2, BICmod3, BICmod4,
             BICmod5, BICmod6, BICmod7, BICmod8,
             ci_level = 0.95,
             digits = 3,
             stars = NULL,
             error_format = '(p = {p.value}), CI [{conf.low}, {conf.high}]',
             model.names = c("Model 1", "Model 2", "Model 3", "Model 4",
                             "Model 5", "Model 6", "Model 7", "Model 8"),
             to.file = "xlsx",
            # file.name = "BICModelResultsSummary_Aug26.xlsx")


#Use model averaging to make a final model

ScatAvg <- model.avg(BICmod1, BICmod2, BICmod3, BICmod4, BICmod5, BICmod6, BICmod7, BICmod8)

#Get model weights
summary(ScatAvg)[1] #Top model is weighted at ~ 0.9, so use it as final model
#instead of reporting averaged results

#Determine whether biologically plausible interaction terms improve performance
ScatAvgMod <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecTrail +
                    DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
                  data = ModDataSet1)
BIC(ScatAvgMod)

IntMod1 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + DecScat:DistWater + DistWater + I(DecScat^2) + 
                 DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)
BIC(IntMod1) #4580.668
BIC(IntMod1) - BIC(BICmod1) #12.79534
summ(IntMod1, confint = TRUE, scale = TRUE)

IntMod2 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecScat:DistNatEdge +
                 DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)
BIC(IntMod2) #4575.451
BIC(IntMod2) - BIC(BICmod1) #7.57827
summ(IntMod2, confint = TRUE, scale = TRUE)

IntMod3 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + DecScat:NatEdDens + NatEdDens + I(DecScat^2) + 
                 DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)
BIC(IntMod3) #4579.455
BIC(IntMod3) - BIC(BICmod1) #11.58261
summ(IntMod3, confint = TRUE, scale = TRUE)

IntMod4 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecScat:East + East +
                 DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)
BIC(IntMod4) #4580.55
BIC(IntMod4) - BIC(BICmod1) #12.6772
summ(IntMod4, confint = TRUE, scale = TRUE)

IntMod5 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + DecScat:Slope + Slope + I(DecScat^2) + 
                 DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)
BIC(IntMod5) #4582.162
BIC(IntMod5) - BIC(BICmod1) #14.28888
summ(IntMod5, confint = TRUE, scale = TRUE)

IntMod6 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecScat:Curve1 +
                 DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)
BIC(IntMod6) #4571.969
BIC(IntMod6) - BIC(BICmod1) #4.095868
summ(IntMod6, confint = TRUE, scale = TRUE)

IntMod7 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecScat:DistNatEdge +
                 DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)
BIC(IntMod7) #4575.451
BIC(IntMod7) - BIC(BICmod1) #7.57827
summ(IntMod7, confint = TRUE, scale = TRUE)

IntMod8 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecScat:RoadDens + RoadDens +
                 DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)
BIC(IntMod8) #4572.64
BIC(IntMod8) - BIC(BICmod1) #4.767384
summ(IntMod8, confint = TRUE, scale = TRUE)

IntMod9 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecScat:DistBldg +
                 DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
               data = ModDataSet1)
BIC(IntMod9) #4576.386
BIC(IntMod9) - BIC(BICmod1) #8.513476
summ(IntMod9, confint = TRUE, scale = TRUE)

IntMod10 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecScat:DecTrail + DecTrail +
                  DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
                data = ModDataSet1)
BIC(IntMod10) #4572.865
BIC(IntMod10) - BIC(BICmod1) #4.99192
summ(IntMod10, confint = TRUE, scale = TRUE)

IntMod11 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecScat:DistJunc + DistJunc +
                  DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
                data = ModDataSet1)
BIC(IntMod11) #4580.564
BIC(IntMod11) - BIC(BICmod1) #12.69078
summ(IntMod11, confint = TRUE, scale = TRUE)

IntMod12 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecScat:DecCamp + DecCamp +
                  DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
                data = ModDataSet1)
BIC(IntMod12) #4576.198
BIC(IntMod12) - BIC(BICmod1) #8.324761
summ(IntMod12, confint = TRUE, scale = TRUE)

IntMod13 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecScat:NAT +
                  DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
                data = ModDataSet1)
BIC(IntMod13) #4551.863
BIC(IntMod13) - BIC(BICmod1) #-16.00936
summ(IntMod13, confint = TRUE, scale = TRUE)

IntMod14 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DistNatEdge:NatEdDens + NatEdDens +
                  DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
                data = ModDataSet1)
BIC(IntMod14) #4581.599
BIC(IntMod14) - BIC(BICmod1) #13.72578
summ(IntMod14, confint = TRUE, scale = TRUE)

IntMod15 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + East:Curve1 + East +
                  DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
                data = ModDataSet1)
BIC(IntMod15) #4579.771
BIC(IntMod15) - BIC(BICmod1) #11.89791
summ(IntMod15, confint = TRUE, scale = TRUE)

IntMod16 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + RoadDens:DistBldg + RoadDens +
                  DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
                data = ModDataSet1)
BIC(IntMod16) #4582.583
BIC(IntMod16) - BIC(BICmod1) #14.70984
summ(IntMod16, confint = TRUE, scale = TRUE)

IntMod17 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecTrail:DistJunc + DecTrail + DistJunc +
                  DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
                data = ModDataSet1)
BIC(IntMod17) #4584.959
BIC(IntMod17) - BIC(BICmod1) #17.08599
summ(IntMod17, confint = TRUE, scale = TRUE)

IntMod18 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecTrail:DecCamp + DecTrail +
                  DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
                data = ModDataSet1)
BIC(IntMod18) #4577.818
BIC(IntMod18) - BIC(BICmod1) #9.945462
summ(IntMod18, confint = TRUE, scale = TRUE)

IntMod19 <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + DecCamp:DistJunc + DistJunc +
                  DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + 1, family = binomial, 
                data = ModDataSet1)
BIC(IntMod19) #4582.611
BIC(IntMod19) - BIC(BICmod1) #14.73774
summ(IntMod19, confint = TRUE, scale = TRUE)

#One interaction improves performance and CI don't overlap zero
#Add that effect to the final model
FinalScatMod <- glm(formula = Use ~ Curve1 + DecCamp + DecScat + I(DecScat^2) + 
                      DistBldg + DistNatEdge + I(DistNatEdge^2) + NAT + DecScat*NAT,
                    family = binomial, 
                    data = ModDataSet1)

#Get beta coefficients, CI, p values
summ(FinalScatMod, digits = 3, confint = TRUE)

#Check collinearity
check_collinearity(FinalScatMod)

#Get pseudo R2
PseudoR2(FinalScatMod, which = "Nagelkerke") #  0.316632 
PseudoR2(FinalScatMod, which = "McFadden") #0.2138827 

#Get ROC
roc.1 <- roc(ModDataSet1$Use~fitted(FinalScatMod))
auc(roc.1) #0.7842

#Conduct cross validation
#create index matrix
index1 <- createDataPartition(ModDataSet1$Use, p = 0.8, list = FALSE, times = 1)

#create train and test
train1 <- ModDataSet1[index1,]
test1 <- ModDataSet1[-index1,]

#convert outcome variable to type factor
train1$Use <- as.factor(train1$Use)
train1$Use2 <- ifelse(train1$Use == "0", "N", "Y")
test1$Use <- as.factor(test1$Use)
test1$Use2 <- ifelse(test1$Use == "0", "N", "Y")

#specify number of folds and train control
ctrlspecs <- trainControl(method = "cv", number = 5,
                          savePredictions = "all",
                          classProbs = TRUE)

#specify logistic regression model
modelo1 <- train(Use2~Curve1 + DecScat + I(DecScat^2) + DistBldg + 
                   DistNatEdge + I(DistNatEdge^2) + NAT + DecCamp + DecScat:NAT,
                 data = train1,
                 method = "glm", family = binomial,
                 trControl = ctrlspecs)
print(modelo1) #  0.8105435  0.3971579
summary(modelo1)
varImp(modelo1)

#apply model to test_df
predictions1 <- predict(modelo1, newdata = test1)
predictions1

#create confusion matrix
test1$Use2 <- factor(test1$Use2, levels = c("N", "Y"))
confusionMatrix(data = predictions1, test1$Use2) #Accuracy : 0.8105435, Kappa : 0.397 

#Create plot showing beta coefficients for final model (INCLUDE P VALUES)
BICPlot1 <- plot_summs(FinalScatMod,
                       colors = c("#00796B"),
                       point.size = 7,
                       point.shape = FALSE,
                       coefs = c("Distance to natural patch edge\nP < 0.001" = "DistNatEdge",
                                 "Decay distance to scat\nP < 0.001" = "DecScat",
                                 "Natural land cover (10-m)\nP < 0.001" = "NAT",
                                 "Decay distance to scat : Natural land cover\nP < 0.001" = "DecScat:NAT",
                                 "Decay distance to scat (Squared)\nP < 0.001" = "I(DecScat^2)",
                                 "Distance to building\nP < 0.001" = "DistBldg",
                                 "Topographic concavity\nP = 0.001" = "Curve1",
                                 "Distance to natural patch edge (Squared)\nP = 0.004" = "I(DistNatEdge^2)",
                                 "Decay distance to camp\nP = 0.003" = "DecCamp"))
BICPlot1

BICPlot2 <- BICPlot1 + theme_classic()

#Improve Plot
BICPlot <- BICPlot2 + xlim(-1.5, 0.5) +
  theme(legend.position = "none",
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Coefficient estimate") +
  scale_shape_manual(values = c(24,24,24,24,24,24,24,24))

BICPlot


#Create plot showing beta coefficients for final model (EXCLUDE P VALUES)
BICPlot1_NP <- plot_summs(FinalScatMod,
                       colors = c("#00796B"),
                       point.size = 7,
                       point.shape = FALSE,
                       coefs = c("Distance to natural\npatch edge" = "DistNatEdge",
                                 "Decay distance\nto scat" = "DecScat",
                                 "Natural land cover\n(10-m)" = "NAT",
                                 "Decay distance to scat :\nNatural land cover" = "DecScat:NAT",
                                 "Decay distance to scat\n(Squared)" = "I(DecScat^2)",
                                 "Distance to\nbuilding" = "DistBldg",
                                 "Topographic\nconcavity" = "Curve1",
                                 "Distance to natural patch\nedge" = "I(DistNatEdge^2)",
                                 "Decay distance\nto camp" = "DecCamp"))
BICPlot1_NP

BICPlot2_NP <- BICPlot1_NP + theme_classic()

#Improve Plot
BICPlot_NP <- BICPlot2_NP + xlim(-1.5, 0.5) +
  theme(legend.position = "none",
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Coefficient estimate") +
  scale_shape_manual(values = c(24,24,24,24,24,24,24,24))

BICPlot_NP

#Create interaction plot to better understand effect
IntPlot1 <- interact_plot(FinalScatMod, pred = NAT, modx = DecScat, interval = TRUE,
                          int.type = "confidence", int.width = .95, colors = c("#004D40", "#4DB6AD"),
                          legend.main = "Decay distance\nto nearest scat")

IntPlot <- IntPlot1 + 
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y=element_text(colour = "black", face = "plain", size = 12),
        axis.title.x=element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x = "Natural land cover within 10 m", y = "Use")

#Save plot
ggsave("IntPlot.png", IntPlot, dpi = 700, width = 6, height = 4)


#Create marginal plots for important effects

#First create model using predictors that aren't scaled/ centered
FinalScatMod$call

UNSCALEDFinalScatMod <- glm(formula = Use ~ Curve1OG + DecCampOG + DecScatOG +
                              I(DecScatOG^2) + DistBldgOG + DistNatEdgeOG + I(DistNatEdgeOG^2) +
                              NATOG + DecScatOG*NATOG, family = binomial, data = MScat)
plot_model(UNSCALEDFinalScatMod, type = "pred", terms = "DistNatEdgeOG")

#DecCamp
DecCampplot1 <- ggpredict(UNSCALEDFinalScatMod, terms = "DecCampOG")
DecCampplot2 <- ggplot(DecCampplot1, aes(x, predicted)) + 
  geom_line(aes(y = predicted), color="#004D40") +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), 
              fill = "#26A69A", alpha = 0.25)

DecCampplot3 <- DecCampplot2 + theme_classic() + ggtitle("(B) Distance to camp")
DecCampplot4 <- DecCampplot3 +
  theme(axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x = "Distance to camp (m)", y = "Likelihood of scat deposition")
DecCampplot4 

#DecScat
Scatplot1 <- ggpredict(UNSCALEDFinalScatMod, terms = "DecScatOG")
Scatplot2 <- ggplot(Scatplot1, aes(x, predicted)) + 
  geom_line(aes(y = predicted), color="#004D40") +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), 
              fill = "#26A69A", alpha = 0.25)

Scatplot3 <- Scatplot2 + theme_classic() + ggtitle("(A) Distance to scat")
Scatplot4 <- Scatplot3 +
  theme(axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x = "Distance to scat (m)", y = "Likelihood of scat deposition")
Scatplot4 


#DistBuilding
DistBldgplot1 <- ggpredict(UNSCALEDFinalScatMod, terms = "DistBldgOG")
DistBldgplot2 <- ggplot(DistBldgplot1, aes(x, predicted)) + 
  geom_line(aes(y = predicted), color="#004D40") +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), 
              fill = "#26A69A", alpha = 0.25)

DistBldgplot3 <- DistBldgplot2 + theme_classic() + ggtitle("(C) Distance to building")
DistBldgplot4 <- DistBldgplot3 +
  theme(axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x = "Distance to building (m)", y = "Likelihood of scat deposition")
DistBldgplot4 

EffectPlots <- ggarrange(Scatplot4, DecCampplot4, DistBldgplot4, nrow = 3, ncol = 1)

#Save plots
ggsave("EffectPlots.png", EffectPlots, width = 4, height = 10, dpi = 700)


# Check fitted values for the model
DD <- fitted(FinalScatMod)
hist(DD) 

################COMPARE TO URINE############################################
#Build a data frame of scats and urinations
stitchP1 <- read.csv('UseScatJun14.csv')
stitchP1$Type <- ifelse(stitchP1$Pee == "Y", 1, 0) #Code pee as 1 and scat as 0
stitchP1 <- stitchP1 %>%
  dplyr::select(FID, DistWater, DistBldg, DistRoad,
                DistEdge, DistNatEdg, DistScat, nateddens2,
                densrd25, densed25, densbldg25, sitetypera,
                curvature1, slope_ju_1, aspectju_1, Type)
stitchP1 <- stitchP1 %>%
  rename("DistNatEdge" = DistNatEdg) %>%
  rename("NatEdDens" = nateddens2) %>%
  rename("EdDens" = densed25) %>%
  rename("RoadDens" = densrd25) %>%
  rename("BldgDens" = densbldg25) %>%
  rename("LC" = sitetypera) %>%
  rename("Curve1" = curvature1) %>%
  rename("Slope" = slope_ju_1) %>%
  rename("Aspect" = aspectju_1)

#LAND COVER
LC_P <- read.csv("LC_Use.csv")
LC_P[is.na(LC_P)] <- 0

stitchP2 <- stitchP1 %>%
  left_join(LC_P, by = "FID")
# OK. we are good on landcover. Let's add scat count

SC_P <- read.csv('ScatCount_U.csv')
SC_P <- SC_P %>%
  dplyr::select(Join_Count, FID) %>%
  rename("ScatCount" = Join_Count)

#Append to those other files
stitchP3 <- stitchP2 %>%
  full_join(SC_P, by = "FID")

##Edges
EdgeP <- read.csv("Edge_U.csv")
stitchP4 <- stitchP3 %>%
  left_join(EdgeP, by = "FID")

NatEdgeP <- read.csv("NatEdge_U.csv")
stitchP5 <- stitchP4 %>%
  left_join(NatEdgeP, by = "FID")

TrailJuncP <- read.csv("Trails_U.csv")
stitchP6 <- stitchP5 %>%
  left_join(TrailJuncP, by = "FID")

WaterP <- read.csv("Water_U.csv")
stitchP7 <- stitchP6 %>%
  left_join(WaterP, by = "FID")

CampGardPGP <- read.csv("CampGardPG_U.csv")
stitchP8 <- stitchP7 %>%
  left_join(CampGardPGP, by = "FID")

#replace NA with 0
stitchP8[is.na(stitchP8)] <- 0

#subtract one from scat count for SCAT sites, but not for urine sites
stitchUrine <- subset(stitchP8, stitchP8$Type == "1")
stitchScat <- subset(stitchP8, stitchP8$Type == "0")

stitchScat$ScatCount <- stitchScat$ScatCount - 1

#Now join to make combined dataset
UrineData <- rbind(stitchScat, stitchUrine)

#Add east index
UrineData <- UrineData %>%
  mutate(AspRad = Aspect*0.01745329251) %>%
  mutate(East = sin(AspRad)) %>%
  mutate(North = cos(AspRad))

#Change scale for road and building density (multiply by 1000)
UrineData <- UrineData %>%
  mutate(RoadDens2 = RoadDens*1000) %>%
  mutate(BldgDens2 = BldgDens*1000)

#Save file
write.csv(UrineData, file = "UrineDataJu30.csv")

###Read in dataset
Urine <- read.csv("UrineDataJu30.csv")
#Compare average values for continuous predictors at scat vs. urine sites
Urine2 <- Urine %>%
  dplyr::select(-c(X, FID, LC, Aspect, AspRad))
Urine2$Type <- factor(Urine2$Type, levels = c("0", "1"))

#Summarise statistics for all these variables
Urine3 <- Urine2 %>% 
  dplyr::group_by(Type) %>% 
  dplyr::summarise_if(is.numeric, .funs=list(mean, sd, min, max)) %>%
  tidyr::pivot_longer(cols = -Type, 
                      names_to = c('.value', 'variable'), 
                      names_sep = '_fn')
Urine3$variable[Urine3$variable==1] <- "mean"
Urine3$variable[Urine3$variable==2] <- "sd"
Urine3$variable[Urine3$variable==3] <- "min"
Urine3$variable[Urine3$variable==4] <- "max"
Urine3 #this is your summary table

#And summarise t test results for quantititaive things
Urine2$Type2 <- ifelse(Urine2$Type == 1, "Urine", "Scat")
A <- 
  Urine2 %>% #specify dataframe
  dplyr::select(-c(Type)) %>% #select quanttaive variables and grouping variable (Use)
  gather(key = variable, value = value, -Type2) %>% #creates one list
  group_by(Type2, variable) %>% 
  summarise(value = list(value)) %>%
  spread(Type2, value) %>% 
  group_by(variable) %>%
  mutate(p_value = t.test(unlist(Scat), unlist(Urine))$p.value,
         t_value = t.test(unlist(Scat), unlist(Urine))$statistic)
#now, write this to a table
UrineTresults <- A %>%
  dplyr::select(variable, p_value, t_value)
UrineTresults

##OK Now combine T results with previous table
#Start by making a SCAT table
SCAT <- data.frame(t(dplyr::filter(Urine3, Type == "0")))
colnames(SCAT) <- as.character(SCAT[2,])
SCAT <- SCAT[-c(1:2),]
SCAT <- tibble::rownames_to_column(SCAT, var = "variable")
SCAT

#And a urine table
URINE <- data.frame(t(dplyr::filter(Urine3, Type == "1")))
colnames(URINE) <- as.character(URINE[2,])
URINE <- URINE[-c(1:2),]
URINE <- tibble::rownames_to_column(URINE, var = "variable")
URINE

##now join all three tables
Urine4 <- UrineTresults %>%
  left_join(SCAT, by = "variable") %>%
  left_join(URINE, by = "variable") #in table, x represents SCAT, and y represents URINE
Urine4

##reorder columns
Urine4 <- Urine4 %>%
  select(mean.x, sd.x, min.x, max.x, mean.y, sd.y, min.y, max.y, t_value, p_value)
Urine4

##rename columns
Urine4 <- Urine4 %>%
  rename("Variable" = variable, "MeanU" = mean.y, "SDU" = sd.y, "MinU" = min.y, "MaxU" = max.y, "MeanS" = mean.x, "SDS" = sd.x, "MinS" = min.x, "MaxS" = max.x, "t-stat" = t_value, "p-value" = p_value)
Urine4

#reorder from smallest to largest p values
Urine4 <- 
  Urine4 %>%
  arrange(`p-value`)

#replace p values < 0.001 with "<0.001"
Urine4$`p-value`[(Urine4$`p-value` <= 0.001)] <- "<0.001"

#Save file
write.csv(Urine4, file = "Urine_vs_Scat_TtestResults_Jul4.csv")

#use final scat model to predict urine, and see how it does
testUrine <- Urine2 %>%
  dplyr::select(DistNatEdge, DistBldg, DistScat, Curve1, Type2, NAT, Type2, DistCamp) %>%
  mutate(DecScat = 1-(exp(-0.06*DistScat))) %>%
  mutate(DecCamp = 1-(exp(-0.012*DistCamp))) %>%
  dplyr::filter(Type2 == "Urine")
testUrine$Use <- "Y"
testUrine$Use <- factor(testUrine$Use, levels = c("Y", "N"))

ModDataSet1$Use2 <- ifelse(ModDataSet1$Use == "0", "N", "Y")
ModDataSet1$Use2 <- factor(ModDataSet1$Use2, levels = c("N", "Y"))

modelo2 <- train(Use2~Curve1 + DecScat + I(DecScat^2) + DistBldg + 
                   DistNatEdge + I(DistNatEdge^2) + NAT + DecCamp + DecScat:NAT,
                 data = ModDataSet1,
                 method = "glm", family = binomial,
                 trControl = ctrlspecs)
print(modelo2)
summary(modelo2)
varImp(modelo2)

#apply model to urine
predictionsUrine <- predict(modelo2, newdata = testUrine)
predictionsUrine

#create confusion matrix
testUrine$Use <- factor(testUrine$Use, levels = c("Y", "N"))
confusionMatrix(data = predictionsUrine, as.factor(testUrine$Use)) #Accuracy : 37.7 %

###########################################################################

##########PART 3 = Green spaces as unit of replication########################

#########################################################################
#Create dataframe to work with
Parks8 <- read.csv("LC_parksBuff50.csv")
Parks8 <- Parks8 %>%
  dplyr::rename("AG50" = AG) %>%
  dplyr::rename("ANTH50" = ANTH) %>%
  dplyr::rename("WATER50" = WATER) %>%
  dplyr::rename("GRASS50" = GRASS) %>%
  dplyr::rename("NAT50" = NAT)
Parks8$FID <- as.character(Parks8$FID)

Parks1 <- read.csv("LC_parksBuff150AB.csv")
Parks1 <- Parks1 %>%
  dplyr::rename("AG150" = AG) %>%
  dplyr::rename("ANTH150" = ANTH) %>%
  dplyr::rename("WATER150" = WATER) %>%
  dplyr::rename("GRASS150" = GRASS) %>%
  dplyr::rename("NAT150" = NAT)
Parks1$FID <- as.character(Parks1$FID)

Parks5 <- read.csv("LC_Parksbuff200.csv")
Parks5 <- Parks5 %>%
  dplyr::rename("AG200" = AG) %>%
  dplyr::rename("ANTH200" = ANTH) %>%
  dplyr::rename("WATER200" = WATER) %>%
  dplyr::rename("GRASS200" = GRASS) %>%
  dplyr::rename("NAT200" = NAT)

Parks2 <- read.csv("ParksBuff50Roaddens.csv")
Parks2$FID <- as.character(Parks2$FID)

Parks3 <- read.csv("ParksBuff150_ABRoaddens.csv")
Parks3$FID <- as.character(Parks3$FID)

Parks7 <- read.csv("ParksBuff200RoadDens.csv")
Parks7$FID <- as.character(Parks7$FID)

Parks6 <- read.csv("ParksBuff50BldgDens.csv")
Parks6$FID <- as.character(Parks6$FID)

Parks9 <- read.csv("ParksBuff150BldgDens.csv")
Parks9$FID <- as.character(Parks9$FID)

Parks4 <- read.csv('ParksBuff200BldgDens.csv')
Parks4$FID <- as.character(Parks4$FID)

Parks11 <- read.csv("ParkPGCount.csv")
Parks11 <- Parks11 %>%
  dplyr::rename('PGCount' = Count_)
Parks11$FID <- as.character(Parks11$FID)

Parks12 <- read.csv("ParkCampCount200.csv")
Parks12 <- Parks12 %>%
  dplyr::rename("CampCount" = Count_)
Parks12$FID <- as.character(Parks12$FID)

Parks13 <- read.csv('ParksCC.csv')
Parks13$FID <- as.character(Parks13$FID)

Parks10 <- read.csv("ParksInfoJul31.csv")
Parks10 <- Parks10 %>%
  dplyr::select(FID, AREA, ParksSca_8, PathLength, Dist_RVR, Dist_Garde,
                Dist_PG, Dist_Camp, Dist_Water, ) %>%
  dplyr::rename("ScatCount" = ParksSca_8) %>%
  dplyr::rename("Dist_Garden" = Dist_Garde)
Parks10$FID <- as.character(Parks10$FID)

#Join them all together
Parks <- list(Parks10, Parks1, Parks2, Parks3, Parks4, Parks5, Parks6,
              Parks7, Parks8, Parks9, Parks11, Parks12, Parks13) %>% 
  reduce(left_join, by = "FID")

#make  new variables
Parks <- Parks %>%
  dplyr::select(-(PathLength.x)) %>%
  dplyr::rename(PathLength = PathLength.y) %>%
  dplyr::mutate(PGperArea = PGCount / Area50) %>%
  dplyr::mutate(CampperLength = (CampCount / PathLength)*1000) %>%
  dplyr::select(-(c(PGCount, CampCount)))

Parks$RVR_YN = ifelse(Parks$Dist_RVR < 100, "Y", "N")

Parks[is.na(Parks)] <- 0

#Save file
write.csv(Parks, file = "ParksAug15.csv")

#Read in file
Parks <- read.csv("ParksAug15.csv")

#Summarise
str(Parks) #126 green spaces

#Remove green spaces with < 25 m walked
Parks <- Parks %>%
  dplyr::filter(PathLength > 25) #123 green spaces with > 25 m

Parks <- Parks %>%
  dplyr::mutate(ScatPerm = ScatCount / PathLength)

sum(Parks$ScatCount) #1009 scats occurred in green spaces
min(Parks$ScatCount) #zero scats in some parks
max(Parks$ScatCount) #80 was max in a green spaces
mean(Parks$ScatCount) #8.2 was the average, but doesn't account for effort

#Summarise scats/ km
min(Parks$ScatPerm)
max(Parks$ScatPerm) #0.09813609
mean(Parks$ScatPerm) #0.00449815
sd(Parks$ScatPerm) #0.01139749

#summarise distance walked
sum(Parks$PathLength) #scats detected over 301595 m = 301.6 km
min(Parks$PathLength) #min distance = 37.15634
max(Parks$PathLength) #as much as 27013.6m 
mean(Parks$PathLength) #mean 2451.992
sd(Parks$PathLength) #4218.899

#Summarise park area
min(Parks$AREA)#2399 m2
max(Parks$AREA)#3935333 m2
mean(Parks$AREA)#184804.8 m2
sd(Parks$AREA)#412262 m2

#Save original copy of all layers, then scale and center
Parks$Dist_RVROG <- Parks$Dist_RVR
Parks$Dist_PGOG <- Parks$Dist_PG
Parks$Dist_CampOG <- Parks$Dist_Camp
Parks$Dist_WaterOG <- Parks$Dist_Water
Parks$Dist_CCOG <- Parks$Dist_CC
Parks$BldgDens50OG <- Parks$BldgDens50
Parks$CampperLengthOG <- Parks$CampperLength
Parks$NAT50OG <- Parks$NAT50

Parks$Dist_RVR <- scale(Parks$Dist_RVR, scale = TRUE, center = TRUE)
Parks$Dist_PG <- scale(Parks$Dist_PG, scale = TRUE, center = TRUE)
Parks$Dist_Camp <- scale(Parks$Dist_Camp, scale = TRUE, center = TRUE)
Parks$Dist_Water <- scale(Parks$Dist_Water, scale = TRUE, center = TRUE)
Parks$RoadDens50 <- scale(Parks$RoadDens50, scale = TRUE, center = TRUE)
Parks$BldgDens50 <- scale(Parks$BldgDens50, scale = TRUE, center = TRUE)
Parks$ANTH50 <- scale(Parks$ANTH50, scale = TRUE, center = TRUE)
Parks$GRASS50 <- scale(Parks$GRASS50, scale = TRUE, center = TRUE)
Parks$NAT50 <- scale(Parks$NAT50, scale = TRUE, center = TRUE)
Parks$Dist_CC <- scale(Parks$Dist_CC, scale = TRUE, center = TRUE)
Parks$AREA <- scale(Parks$AREA, scale = TRUE, center = TRUE)
Parks$PGperArea <- scale(Parks$PGperArea, scale = TRUE, center = TRUE)
Parks$CampperLength <- scale(Parks$CampperLength, scale = TRUE, center = TRUE)


Parks$Dist_RVR <- as.numeric(Parks$Dist_RVR)
Parks$Dist_PG <- as.numeric(Parks$Dist_PG)
Parks$Dist_Camp <- as.numeric(Parks$Dist_Camp)
Parks$Dist_Water <- as.numeric(Parks$Dist_Water)
Parks$RoadDens50 <- as.numeric(Parks$RoadDens50)
Parks$BldgDens50 <- as.numeric(Parks$BldgDens50)
Parks$ANTH50 <- as.numeric(Parks$ANTH50)
Parks$GRASS50 <- as.numeric(Parks$GRASS50)
Parks$NAT50 <- as.numeric(Parks$NAT50)
Parks$Dist_CC <- as.numeric(Parks$Dist_CC)
Parks$AREA <- as.numeric(Parks$AREA)
Parks$PGperArea <- as.numeric(Parks$PGperArea)
Parks$CampperLength <- as.numeric(Parks$CampperLength)


#Likely overdispersed, so use negative binomial distribution
#Create univariate GLMs
Parkmod1 <- glm.nb(ScatCount ~ Dist_RVR + offset(log(PathLength)), data = Parks)
Parkmod3 <- glm.nb(ScatCount ~ Dist_PG + offset(log(PathLength)), data = Parks)
Parkmod4 <- glm.nb(ScatCount ~ Dist_Camp + offset(log(PathLength)), data = Parks)
Parkmod5 <- glm.nb(ScatCount ~ Dist_Water + offset(log(PathLength)), data = Parks)
Parkmod6 <- glm.nb(ScatCount ~ RVR_YN + offset(log(PathLength)), data = Parks)

Parkmod7 <- glm.nb(ScatCount ~ AG50 + offset(log(PathLength)), data = Parks)
Parkmod8 <- glm.nb(ScatCount ~ ANTH50 + offset(log(PathLength)), data = Parks)
Parkmod9 <- glm.nb(ScatCount ~ WATER50 + offset(log(PathLength)), data = Parks)
Parkmod10 <- glm.nb(ScatCount ~ GRASS50 + offset(log(PathLength)), data = Parks)
Parkmod11 <- glm.nb(ScatCount ~ NAT50 + offset(log(PathLength)), data = Parks)

Parkmod12 <- glm.nb(ScatCount ~ RoadDens50 + offset(log(PathLength)), data = Parks)
Parkmod13 <- glm.nb(ScatCount ~ BldgDens50 + offset(log(PathLength)), data = Parks)

Parkmod14 <- glm.nb(ScatCount ~ PGperArea + offset(log(PathLength)), data = Parks)
Parkmod15 <- glm.nb(ScatCount ~ CampperLength + offset(log(PathLength)), control = glm.control(maxit = 50), data = Parks)

Parkmod16 <- glm.nb(ScatCount ~ Dist_CC + offset(log(PathLength)), data = Parks)

Parkmod16a <- glm.nb(ScatCount ~ AREA + offset(log(PathLength)), data = Parks)


summary(Parkmod1)
summary(Parkmod3)
summary(Parkmod4)
summary(Parkmod5)
summary(Parkmod6)

summary(Parkmod7)
summary(Parkmod8)
summary(Parkmod9)
summary(Parkmod10)
summary(Parkmod11)

summary(Parkmod12)
summary(Parkmod13)
summary(Parkmod14)
summary(Parkmod15)
summary(Parkmod16)
summary(Parkmod16a)

#Decay stuff
ParkDecmod1a <- glm.nb(ScatCount ~ exp(-0.002*Dist_RVROG) + offset(log(PathLength)), data = Parks)
ParkDecmod3a <- glm.nb(ScatCount ~ exp(-0.002*Dist_PGOG) + offset(log(PathLength)), data = Parks)
ParkDecmod4a <- glm.nb(ScatCount ~ exp(-0.002*Dist_CampOG) + offset(log(PathLength)), data = Parks)
ParkDecmod5a <- glm.nb(ScatCount ~ exp(-0.002*Dist_WaterOG) + offset(log(PathLength)), data = Parks)
ParkDecmod16a <- glm.nb(ScatCount ~ exp(-0.002*Dist_CCOG) + offset(log(PathLength)), data = Parks)

ParkDecmod1b <- glm.nb(ScatCount ~ exp(-0.004*Dist_RVROG) + offset(log(PathLength)), data = Parks)
ParkDecmod3b <- glm.nb(ScatCount ~ exp(-0.004*Dist_PGOG) + offset(log(PathLength)), data = Parks)
ParkDecmod4b <- glm.nb(ScatCount ~ exp(-0.004*Dist_CampOG) + offset(log(PathLength)), data = Parks)
ParkDecmod5b <- glm.nb(ScatCount ~ exp(-0.004*Dist_WaterOG) + offset(log(PathLength)), data = Parks)
ParkDecmod16b <- glm.nb(ScatCount ~ exp(-0.004*Dist_CCOG) + offset(log(PathLength)), data = Parks)

ParkDecmod1c <- glm.nb(ScatCount ~ exp(-0.006*Dist_RVROG) + offset(log(PathLength)), data = Parks)
ParkDecmod3c <- glm.nb(ScatCount ~ exp(-0.006*Dist_PGOG) + offset(log(PathLength)), data = Parks)
ParkDecmod4c <- glm.nb(ScatCount ~ exp(-0.006*Dist_CampOG) + offset(log(PathLength)), data = Parks)
ParkDecmod5c <- glm.nb(ScatCount ~ exp(-0.006*Dist_WaterOG) + offset(log(PathLength)), data = Parks)
ParkDecmod16c <- glm.nb(ScatCount ~ exp(-0.006*Dist_CCOG) + offset(log(PathLength)), data = Parks)

ParkDecmod1d <- glm.nb(ScatCount ~ exp(-0.012*Dist_RVROG) + offset(log(PathLength)), data = Parks)
ParkDecmod3d <- glm.nb(ScatCount ~ exp(-0.012*Dist_PGOG) + offset(log(PathLength)), data = Parks)
ParkDecmod4d <- glm.nb(ScatCount ~ exp(-0.012*Dist_CampOG) + offset(log(PathLength)), data = Parks)
ParkDecmod5d <- glm.nb(ScatCount ~ exp(-0.012*Dist_WaterOG) + offset(log(PathLength)), data = Parks)
ParkDecmod16d <- glm.nb(ScatCount ~ exp(-0.012*Dist_CCOG) + offset(log(PathLength)), data = Parks)

ParkDecmod1e <- glm.nb(ScatCount ~ exp(-0.03*Dist_RVROG) + offset(log(PathLength)), data = Parks)
ParkDecmod3e <- glm.nb(ScatCount ~ exp(-0.03*Dist_PGOG) + offset(log(PathLength)), data = Parks)
ParkDecmod4e <- glm.nb(ScatCount ~ exp(-0.03*Dist_CampOG) + offset(log(PathLength)), data = Parks)
ParkDecmod5e <- glm.nb(ScatCount ~ exp(-0.03*Dist_WaterOG) + offset(log(PathLength)), data = Parks)
ParkDecmod16e <- glm.nb(ScatCount ~ exp(-0.03*Dist_CCOG) + offset(log(PathLength)), data = Parks)

ParkDecmod1f <- glm.nb(ScatCount ~ exp(-0.06*Dist_RVROG) + offset(log(PathLength)), data = Parks)
ParkDecmod3f <- glm.nb(ScatCount ~ exp(-0.06*Dist_PGOG) + offset(log(PathLength)), data = Parks)
ParkDecmod4f <- glm.nb(ScatCount ~ exp(-0.06*Dist_CampOG) + offset(log(PathLength)), data = Parks)
ParkDecmod5f <- glm.nb(ScatCount ~ exp(-0.06*Dist_WaterOG) + offset(log(PathLength)), data = Parks)
ParkDecmod16f <- glm.nb(ScatCount ~ exp(-0.06*Dist_CCOG) + offset(log(PathLength)), data = Parks)

ParkDecmod1g <- glm.nb(ScatCount ~ exp(-0.2*Dist_RVROG) + offset(log(PathLength)), data = Parks)
ParkDecmod3g <- glm.nb(ScatCount ~ exp(-0.2*Dist_PGOG) + offset(log(PathLength)), data = Parks)
ParkDecmod4g <- glm.nb(ScatCount ~ exp(-0.2*Dist_CampOG) + offset(log(PathLength)), data = Parks)
ParkDecmod5g <- glm.nb(ScatCount ~ exp(-0.2*Dist_WaterOG) + offset(log(PathLength)), data = Parks)
ParkDecmod16g <- glm.nb(ScatCount ~ exp(-0.2*Dist_CCOG) + offset(log(PathLength)), data = Parks)

AIC(Parkmod1); AIC(ParkDecmod1a); AIC(ParkDecmod1b); AIC(ParkDecmod1c); AIC(ParkDecmod1d); AIC(ParkDecmod1e); AIC(ParkDecmod1f); AIC(ParkDecmod1g)
#dist RVR best without decay

AIC(Parkmod3); AIC(ParkDecmod3a); AIC(ParkDecmod3b); AIC(ParkDecmod3c); AIC(ParkDecmod3d); AIC(ParkDecmod3e); AIC(ParkDecmod3f); AIC(ParkDecmod3g)
#0.2 best for Dist_PG

AIC(Parkmod4); AIC(ParkDecmod4a); AIC(ParkDecmod4b); AIC(ParkDecmod4c); AIC(ParkDecmod4d); AIC(ParkDecmod4e); AIC(ParkDecmod4f); AIC(ParkDecmod4g)
#0.012 best for dist_Camp

AIC(Parkmod5); AIC(ParkDecmod5a); AIC(ParkDecmod5b); AIC(ParkDecmod5c); AIC(ParkDecmod5d); AIC(ParkDecmod5e); AIC(ParkDecmod5f); AIC(ParkDecmod5g)
#decay 0.006 best for dist water

AIC(Parkmod16); AIC(ParkDecmod16a); AIC(ParkDecmod16b); AIC(ParkDecmod16c); AIC(ParkDecmod16d); AIC(ParkDecmod16e); AIC(ParkDecmod16f); AIC(ParkDecmod16g)
#decay 0.002 best for dist_cc

#Make decay terms
Parks <- Parks %>%
  mutate(Dec_PG = 1-(exp(-0.2*Dist_PGOG))) %>%
  mutate(Dec_Camp = 1-(exp(-0.012*Dist_CampOG))) %>%
  mutate(Dec_Water = 1-(exp(-0.006*Dist_WaterOG))) %>%
  mutate(Dec_CC = 1-(exp(-0.002*Dist_CCOG)))

Parks$Dec_PGOG <- Parks$Dec_PG
Parks$Dec_CampOG <- Parks$Dec_Camp
Parks$Dec_WaterOG <- Parks$Dec_Water
Parks$Dec_CCOG <- Parks$Dec_CC

Parks$Dec_PG <- scale(Parks$Dec_PG, scale = TRUE, center = TRUE)
Parks$Dec_Camp <- scale(Parks$Dec_Camp, scale = TRUE, center = TRUE)
Parks$Dec_Water <- scale(Parks$Dec_Water, scale = TRUE, center = TRUE)
Parks$Dec_CC <- scale(Parks$Dec_CC, scale = TRUE, center = TRUE)

Parks$Dec_Camp <- as.numeric(Parks$Dec_Camp)
Parks$Dec_PG <- as.numeric(Parks$Dec_PG)
Parks$Dec_Water <- as.numeric(Parks$Dec_Water)
Parks$Dec_CC <- as.numeric(Parks$Dec_CC)

Parkmod17 <- glm.nb(ScatCount ~ Dec_Camp + offset(log(PathLength)), data = Parks)
Parkmod18 <- glm.nb(ScatCount ~ Dec_PG + offset(log(PathLength)), data = Parks)
Parkmod19 <- glm.nb(ScatCount ~ Dec_Water + offset(log(PathLength)), data = Parks)
Parkmod20 <- glm.nb(ScatCount ~ Dec_CC + offset(log(PathLength)), data = Parks)

summary(Parkmod17)
summary(Parkmod18)
summary(Parkmod19)
summary(Parkmod20)

ParkQmod1 <- glm.nb(ScatCount ~ Dist_RVR + I(Dist_RVR^2) + offset(log(PathLength)), data = Parks)
ParkQmod3 <- glm.nb(ScatCount ~ Dist_PG + I(Dist_PG^2) + offset(log(PathLength)), data = Parks)
ParkQmod4 <- glm.nb(ScatCount ~ Dist_Camp + I(Dist_Camp^2) + offset(log(PathLength)), data = Parks)
ParkQmod5 <- glm.nb(ScatCount ~ Dist_Water + I(Dist_Water^2) + offset(log(PathLength)), data = Parks)

ParkQmod7 <- glm.nb(ScatCount ~ AG50 + I(AG50^2) + offset(log(PathLength)), data = Parks)
ParkQmod8 <- glm.nb(ScatCount ~ ANTH50 + I(ANTH50^2) + offset(log(PathLength)), data = Parks)
ParkQmod9 <- glm.nb(ScatCount ~ WATER50 + I(WATER50^2) + offset(log(PathLength)), data = Parks)
ParkQmod10 <- glm.nb(ScatCount ~ GRASS50 + I(GRASS50^2) + offset(log(PathLength)), data = Parks)
ParkQmod11 <- glm.nb(ScatCount ~ NAT50 + I(NAT50^2) + offset(log(PathLength)), data = Parks)

ParkQmod12 <- glm.nb(ScatCount ~ RoadDens50 + I(RoadDens50^2) + offset(log(PathLength)), data = Parks)
ParkQmod13 <- glm.nb(ScatCount ~ BldgDens50 + I(BldgDens50^2) + offset(log(PathLength)), data = Parks, maxit = 100)

ParkQmod14 <- glm.nb(ScatCount ~ PGperArea + I(PGperArea^2) + offset(log(PathLength)), data = Parks)
ParkQmod15 <- glm.nb(ScatCount ~ CampperLength + I(CampperLength^2) + offset(log(PathLength)), data = Parks)

ParkQmod16 <- glm.nb(ScatCount ~ Dist_CC + I(Dist_CC^2) + offset(log(PathLength)), data = Parks)
ParkQmod16a <- glm.nb(ScatCount ~ AREA + I(AREA^2) + offset(log(PathLength)), data = Parks)

ParkQmod17 <- glm.nb(ScatCount ~ Dec_Camp + I(Dec_Camp^2) + offset(log(PathLength)), data = Parks)
ParkQmod18 <- glm.nb(ScatCount ~ Dec_PG + I(Dec_PG^2) + offset(log(PathLength)), data = Parks)
ParkQmod19 <- glm.nb(ScatCount ~ Dec_Water + I(Dec_Water^2) + offset(log(PathLength)), data = Parks)
ParkQmod20 <- glm.nb(ScatCount ~ Dec_CC + I(Dec_CC^2) + offset(log(PathLength)), data = Parks)


summary(ParkQmod1)
summary(ParkQmod3)
summary(ParkQmod4)
summary(ParkQmod5)

summary(ParkQmod7)
summary(ParkQmod8)
summary(ParkQmod9)
summary(ParkQmod10)
summary(ParkQmod11)

summary(ParkQmod12)
summary(ParkQmod13)

summary(ParkQmod14)
summary(ParkQmod15)

summary(ParkQmod16)
summary(ParkQmod16a)

summary(ParkQmod17)
summary(ParkQmod18)
summary(ParkQmod19)
summary(ParkQmod20)

AIC(Parkmod1); AIC(ParkQmod1) #Distane to RVR NOT RETAINED
AIC(Parkmod3); AIC(ParkQmod3) #don't retain dist to pg
AIC(Parkmod4); AIC(ParkQmod4) #don't retain dist to camp
AIC(Parkmod5); AIC(ParkQmod5) #retain dist water as quad

summary(Parkmod6) #Don't retain RV YN

AIC(Parkmod7); AIC(ParkQmod7) # ag 50 not sig
AIC(Parkmod8); AIC(ParkQmod8) # Anth 50 is sig lin
AIC(Parkmod9); AIC(ParkQmod9) # water 50 is sig lin
AIC(Parkmod10); AIC(ParkQmod10) # grass 50 is sig lin
AIC(Parkmod11); AIC(ParkQmod11) #Retain nat 50 as quad

AIC(Parkmod12); AIC(ParkQmod12) #quad road 50 sig 
AIC(Parkmod13); AIC(ParkQmod13) #lin bldg 50; sig

AIC(Parkmod14); AIC(ParkQmod14) #retain PG per area
AIC(Parkmod15); AIC(ParkQmod15) #quad but a sus model
AIC(Parkmod16); AIC(ParkQmod16) #quad; retain dist cc as quad
AIC(Parkmod16a); AIC(ParkQmod16a) #retain area 
AIC(Parkmod17); AIC(ParkQmod17) #Decay PG retain lin
AIC(Parkmod18); AIC(ParkQmod18) #decay camp as lin
AIC(Parkmod19); AIC(ParkQmod19) #decay water as lin
AIC(Parkmod20); AIC(ParkQmod20) #dec cc as lin


#Check for correlation (use r > 0.6)
PData <- Parks %>%
  dplyr::select(Dist_Garden, ANTH50, WATER50, GRASS50, NAT50, RoadDens50, BldgDens50,
                CampperLength, Dec_PG, Dec_Camp, Dist_Water, Dist_CC, PGperArea, PathLength,
                ScatCount) %>%
  dplyr::filter(PathLength > 25) 

cor(PData[sapply(PData, is.numeric)], method = c("pearson"))
#with 0.6, Grass and Nat can't be in the same model

#Make global model
GlobParkModP <- glm(ScatCount ~ ANTH50 + NAT50 + 
                      RoadDens50 + BldgDens50 + CampperLength + Dec_Camp + 
                      Dist_Water + Dist_CC + I(NAT50^2)+ I(CampperLength^2) +
                      I(RoadDens50^2) + I(Dist_CC^2) + I(Dist_Water^2) +
                      offset(log(PathLength)), 
                    data = PData, family = poisson)
check_overdispersion(GlobParkModP) #disp ratio = 7.230, P < .001, OD detected

GlobParkMod <- glm.nb(ScatCount ~ ANTH50 + NAT50 + 
                        RoadDens50 + BldgDens50 + CampperLength + Dec_Camp + 
                        Dist_Water + Dist_CC + I(NAT50^2)+ I(CampperLength^2) +
                        I(RoadDens50^2) + I(Dist_CC^2) + I(Dist_Water^2) +
                        offset(log(PathLength)), 
                      data = PData)
check_overdispersion(GlobParkMod) #disp ratio = 1.230, P = 0.052 (no od detected)

subsetPark <- expression(dc(`NAT50`, `I(NAT50^2)`) &
                           dc(`RoadDens50`, `I(RoadDens50^2)`) &
                           dc(`CampperLength`, `I(CampperLength^2)`) &
                           dc(`Dist_Water`, `I(Dist_Water^2)`) &
                           dc(`Dist_CC`, `I(Dist_CC^2)`))

#dredge model
options(na.action = "na.fail")
ParksResults2 <- dredge(GlobParkMod, subset = subsetPark, trace = 2, rank = BIC)

#Get BIC results (two models within 2 BIC)
subset(ParksResults2, delta < 2)

#The first 2 models are within 2 BIC
PBICmod1a <- get.models(ParksResults2, 1)[[1]]
PBICmod2a <- get.models(ParksResults2, 2)[[1]]

PBICmod1a$call
PBICmod1 <- glm.nb(formula = ScatCount ~ BldgDens50 + CampperLength + I(CampperLength^2) + 
                     Dist_CC + NAT50 + I(NAT50^2) + offset(log(PathLength)), data = PData)
PBICmod2a$call
PBICmod2 <- glm.nb(formula = ScatCount ~ BldgDens50 + CampperLength + I(CampperLength^2) + 
                     Dist_CC + NAT50 + offset(log(PathLength)), data = PData)

summ(PBICmod1, digits = 4, confint = TRUE)
PseudoR2(PBICmod1, which = "Nagelkerke") #0.3955961 
PseudoR2(PBICmod1, which = "McFadden") #0.09543276 

summ(PBICmod2, digits = 4, confint = TRUE)
PseudoR2(PBICmod2, which = "Nagelkerke") #0.3657552 
PseudoR2(PBICmod2, which = "McFadden") #0.08631392 

#Create Plot (WITHOUT P VALUES)
BICPlotGS_NP <- plot_summs(PBICmod1, PBICmod2,
                        colors = c("#00796B", "#00796B"),
                        point.size = 7,
                        point.shape = FALSE,
                        coefs = c("Camps / km" = "CampperLength",
                                  "Natural land\ncover" = "NAT50",
                                  "Distance to\ncity center" = "Dist_CC",
                                  "Building\ndensity" = "BldgDens50",
                                  "Natural land cover\n(Squared)" = "I(NAT50^2)",
                                  "Camps / km\n(Squared)" = "I(CampperLength^2)"))
BICPlotGS_NP
BICPlotGS1_NP <- BICPlotGS_NP + theme_classic()
BICPlotGS2_NP <- BICPlotGS1_NP + 
  theme(legend.position = "none",
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Coefficient estimate") +
  scale_shape_manual(values = c(24,24))

BICPlotGS2_NP

#Save plot
ggsave("BICPlotGS_NP.png", BICPlotGS2, width = 5, height = 4, dpi = 700)

#Create Plot (WITH P VALUES)
BICPlotGS <- plot_summs(PBICmod1, PBICmod2,
                        colors = c("#00796B", "#00796B"),
                        point.size = 7,
                        point.shape = FALSE,
                        coefs = c("Camps / km\nP = 0.018" = "CampperLength",
                                  "Natural land cover\nP < 0.001" = "NAT50",
                                  "Distance to city center\nP = 0.001" = "Dist_CC",
                                  "Building density\nP = 0.001" = "BldgDens50",
                                  "Natural land cover (Squared)\nP = 0.013" = "I(NAT50^2)",
                                  "Camps / km (Squared)\nP = 0.418" = "I(CampperLength^2)"))
BICPlotGS
BICPlotGS1 <- BICPlotGS + theme_classic()
BICPlotGS2 <- BICPlotGS1 + 
  theme(legend.position = "none",
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Coefficient estimate") +
  scale_shape_manual(values = c(24,24))

BICPlotGS2

#Save plot
ggsave("BICPlotGS.png", BICPlotGS2, width = 5, height = 4, dpi = 700)


#Average models
GreenSpaceAcg <- model.avg(PBICmod1, PBICmod2, rank = BIC)
summary(GreenSpaceAcg)
confint(GreenSpaceAcg, subset = TRUE)

#check collinearity
check_collinearity(PBICmod1)
check_collinearity(PBICmod2)

#Get pseudo r2
PseudoR2(PBICmod1, which = "Nagelkerke") #0.3955961 
PseudoR2(PBICmod2, which = "Nagelkerke") #0.3657552 

PseudoR2(PBICmod1, which = "McFadden") #0.09543276 
PseudoR2(PBICmod2, which = "McFadden") #0.08631392 

#Conduct multiclss ROC AUC
roc.2 <- multiclass.roc(Parks$ScatCount~fitted(PBICmod1))
auc(roc.2) #0.8748

roc.3 <- multiclass.roc(Parks$ScatCount~fitted(PBICmod2))
auc(roc.3) #0.8699


#Cross validation:
data_ctrl <- trainControl(method = "cv", number = 5)
model_caret <- train(ScatCount ~ BldgDens50 + CampperLength + I(CampperLength^2) + 
                       Dist_CC + NAT50 + I(NAT50^2) + offset(log(PathLength)), 
                     data = PData, trControl = data_ctrl,  method = "glm.nb", na.action = na.pass)
print(model_caret)
model_caret
model_caret$finalModel
model_caret$resample
sd(model_caret$resample$Rsquared)


#Do interaction terms improve model performance?
INTmod1 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                    Dist_CCS + NAT50S + I(NAT50S^2) + RVR_YN + RVR_YN*Dist_CCS + 
                    offset(log(PathLength)), data = Parks)
BIC(INTmod1) - BIC(PBICmod1)
confint(INTmod1, digits = 2)


INTmod2 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                    Dist_CCS + NAT50S + I(NAT50S^2) + RVR_YN + RVR_YN*NAT50S + 
                    offset(log(PathLength)), data = Parks)
BIC(INTmod2) - BIC(PBICmod1)
confint(INTmod2, digits = 2)

INTmod3 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                    Dist_CCS + NAT50S + I(NAT50S^2) + Dist_RVRS + Dist_WaterS + 
                    Dist_RVRS:Dist_WaterS + offset(log(PathLength)), data = Parks)
BIC(INTmod3) - BIC(PBICmod1)
confint(INTmod3, digits = 2)


INTmod4 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                    Dist_CCS + NAT50S + I(NAT50S^2) + Dist_RVRS  + 
                    Dist_RVRS:Dist_CCS + offset(log(PathLength)), data = Parks)
BIC(INTmod4) - BIC(PBICmod1)
confint(INTmod4, digits = 2)


INTmod5 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                    Dist_CCS + NAT50S + I(NAT50S^2) + Dist_RVRS  + 
                    Dist_RVRS:NAT50S + offset(log(PathLength)), data = Parks)
BIC(INTmod5) - BIC(PBICmod1)
confint(INTmod5, digits = 2)


INTmod6 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                    Dist_CCS + NAT50S + I(NAT50S^2) + Dist_WaterS  + 
                    Dist_WaterS:NAT50S + offset(log(PathLength)), data = Parks)
BIC(INTmod6) - BIC(PBICmod1)
confint(INTmod6, digits = 2)


INTmod7 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                    Dist_CCS + NAT50S + I(NAT50S^2) + Dist_WaterS  + AREAS + 
                    Dist_WaterS:AREAS + offset(log(PathLength)), data = Parks)
BIC(INTmod7) - BIC(PBICmod1)
confint(INTmod7, digits = 2)


INTmod8 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                    Dist_CCS + NAT50S + I(NAT50S^2) + RoadDens50S + 
                    Dist_CCS:RoadDens50S + offset(log(PathLength)), data = Parks)
BIC(INTmod8) - BIC(PBICmod1)
confint(INTmod8, digits = 2)


INTmod9 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                    Dist_CCS + NAT50S + I(NAT50S^2) + 
                    Dist_CCS:BldgDens50S + offset(log(PathLength)), data = Parks)
BIC(INTmod9) - BIC(PBICmod1)
confint(INTmod9, digits = 2)


INTmod10 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) + 
                     Dist_CCS:CampperLengthS + offset(log(PathLength)), data = Parks)
BIC(INTmod10) - BIC(PBICmod1)
confint(INTmod10, digits = 2)


INTmod11 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) + Dec_CampS + RoadDens50S +
                     Dec_CampS:RoadDens50S + offset(log(PathLength)), data = Parks)
BIC(INTmod11) - BIC(PBICmod1)
confint(INTmod11, digits = 2)


INTmod12 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) + ANTH50S +
                     NAT50S:ANTH50S + offset(log(PathLength)), data = Parks)
BIC(INTmod12) - BIC(PBICmod1)
confint(INTmod12, digits = 2)


INTmod13 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) + Dec_CampS +
                     Dec_CampS:BldgDens50S + offset(log(PathLength)), data = Parks)
BIC(INTmod13) - BIC(PBICmod1)
confint(INTmod13, digits = 2)


INTmod14 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) + Dec_CampS +
                     Dec_CampS:CampperLengthS + offset(log(PathLength)), data = Parks)
BIC(INTmod14) - BIC(PBICmod1)
confint(INTmod14, digits = 2)


INTmod15 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) + AREAS +
                     NAT50S:AREAS + offset(log(PathLength)), data = Parks)
BIC(INTmod15) - BIC(PBICmod1)
confint(INTmod15, digits = 2)


INTmod16 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) + AREAS + RoadDens50S +
                     RoadDens50S:AREAS + offset(log(PathLength)), data = Parks)
BIC(INTmod16) - BIC(PBICmod1)
confint(INTmod16, digits = 2)


INTmod17 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) + AREAS +
                     BldgDens50S:AREAS + offset(log(PathLength)), data = Parks)
BIC(INTmod17) - BIC(PBICmod1)
confint(INTmod17, digits = 2)


INTmod18 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) + AREAS +
                     Dist_CCS:AREAS + offset(log(PathLength)), data = Parks)
BIC(INTmod18) - BIC(PBICmod1)
confint(INTmod18, digits = 2)


INTmod19 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) +
                     Dist_CCS:CampperLengthS + offset(log(PathLength)), data = Parks)
BIC(INTmod19) - BIC(PBICmod1)
confint(INTmod19, digits = 2)


INTmod20 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) + RoadDens50S +
                     CampperLengthS:RoadDens50S + offset(log(PathLength)), data = Parks)
BIC(INTmod20) - BIC(PBICmod1)
confint(INTmod20, digits = 2)


INTmod21 <- glm.nb(ScatCount ~ BldgDens50S + CampperLengthS + I(CampperLengthS^2) + 
                     Dist_CCS + NAT50S + I(NAT50S^2) +
                     CampperLengthS:BldgDens50S + offset(log(PathLength)), data = Parks)
BIC(INTmod21) - BIC(PBICmod1)
confint(INTmod21, digits = 2)

PBICmod1_UNSCALED <- glm.nb(formula = ScatCount ~ BldgDens50OG + CampperLengthOG + I(CampperLengthOG^2) + 
                     Dist_CCOG + NAT50OG + I(NAT50OG^2) + offset(log(PathLength)), data = Parks)
summary(PBICmod1_UNSCALED)

#Combine beta coefficient plot for your two different units of replication (scat, green space)
BICPlot
BICPlotGS2

BICPlot_annot <- BICPlot + ggtitle("(A) Scats") +
  theme(plot.title = element_text(color = "black",
                                  size = 12,
                                  face = "plain",
                                  hjust = -1.35))

BICPlotGS_annot <- BICPlotGS2 + ggtitle("(B) Green spaces") +
  theme(plot.title = element_text(color = "black",
                                  size = 12,
                                  face = "plain",
                                  hjust = -1.45))

BICmods <- ggarrange(BICPlot_annot, BICPlotGS_annot, ncol = 2, widths = c(5,4))

#Save file
ggsave(BICmods, dpi = 600, file = "ScatPredBICModsFeb19.png", width = 11.7, height = 6)


#Repeat for no p value plots
BICPlot_NP
BICPlotGS2_NP

BICPlot_annot_NP <- BICPlot_NP + ggtitle("(A) Scats") +
  theme(plot.title = element_text(color = "black",
                                  size = 12,
                                  face = "plain",
                                  hjust = -0))

BICPlotGS_annot_NP <- BICPlotGS2_NP + ggtitle("(B) Green spaces") +
  theme(plot.title = element_text(color = "black",
                                  size = 12,
                                  face = "plain",
                                  hjust = 0))

BICmods_NP <- ggarrange(BICPlot_annot_NP, BICPlotGS_annot_NP, ncol = 2, widths = c(3.5,3))

#Save file
ggsave(BICmods_NP, dpi = 600, file = "ScatPredBICModsMay21.png", width = 9, height = 6)

###########################################################################
#PART 4 = CONTENT
###########################################################################
#Read in data
Content1 <- read.csv('UseScatJun14.csv')
Content1 <- Content1 %>%
  dplyr::filter(Content1$Pee != "Y") %>%
  dplyr::select(-c(ID, Timestamp, DateTime, Year, Month, Day, LocationDe, Lat, Long,
                   TOSS, Comments_, NEAR_FID, NEAR_DIST, NEAR_X, NEAR_Y, curvatur_1,
                   curvatur_2, curvatur_3, curvatur_4, curvatur_5, curvatur_6, 
                   BUFF_DIST, ORIG_FID))
Content1 <- Content1 %>%
  rename("DistNatEdge" = DistNatEdg) %>%
  rename("NatEdDens" = nateddens2) %>%
  rename("EdDens" = densed25) %>%
  rename("RoadDens" = densrd25) %>%
  rename("BldgDens" = densbldg25) %>%
  rename("LC" = sitetypera) %>%
  rename("Curve1" = curvature1) %>%
  rename("Slope" = slope_ju_1) %>%
  rename("Aspect" = aspectju_1)

#LAND COVER
LC_Content <- read.csv("LC_Use.csv")
LC_Content[is.na(LC_Content)] <- 0

Content2 <- Content1 %>%
  left_join(LC_Content, by = "FID")

SC_Content <- read.csv('ScatCount_U.csv')
SC_Content <- SC_Content %>%
  dplyr::filter(SC_Content$Pee != "Y") %>%
  dplyr::select(Join_Count, FID) %>%
  rename("ScatCount" = Join_Count)

Content3 <- Content2 %>%
  full_join(SC_Content, by = "FID")

##Edges
EdgeContent <- read.csv("Edge_U.csv")

Content4 <- Content3 %>%
  left_join(EdgeContent, by = "FID")

NatEdgeContent <- read.csv("NatEdge_U.csv")
Content5 <- Content4 %>%
  left_join(NatEdgeContent, by = "FID")

TrailJuncContent <- read.csv("Trails_U.csv")
Content6 <- Content5 %>%
  left_join(TrailJuncContent, by = "FID")

WaterContent <- read.csv("Water_U.csv")
Content7 <- Content6 %>%
  left_join(WaterContent, by = "FID")

PGCGCampContent <- read.csv("CampGardPG_U.csv")
Content8 <- Content7 %>%
  left_join(PGCGCampContent, by = "FID")

#replace NA with 0
Content8[is.na(Content8)] <- 0

#subtract one from scat count for use sites
Content8$ScatCount <- Content8$ScatCount - 1

#And remove entries with no content info
Content8 <- Content8 %>%
  dplyr::filter(Content != "DNR")

#And add eastness and northness
Content <- Content8 %>%
  mutate(AspRad = Aspect*0.01745329251) %>%
  mutate(East = sin(AspRad)) %>%
  mutate(North = cos(AspRad))

#Save file
write.csv(Content, file = "ScatContentJu30.csv")

#Read in content data frame
Cont <- read.csv("ScatContentJu30.csv")
str(Cont)
Cont <- Cont %>%
  dplyr::select(-c(X, FID, Pee, Content))


#Summary: what percent of scats contain each type of content?
Fruit <- (sum(Cont$Fruit == "Y") / 1173) * 100 #620 scats
Fruit #52.9%

Anthropogenic <- (sum(Cont$Anthro == "Y") / 1173) * 100 # 431 scats
Anthropogenic #36/7%

Garbage <- (sum(Cont$Garbage == "Y") / 1173) * 100  #37 scats
Garbage #3.2%

NaturalPrey <- (sum(Cont$NatPrey == "Y") / 1173) * 100 #184
NaturalPrey #15.7%

Birdseed <- (sum(Cont$Seed == "Y") / 1173) * 100 #188
Birdseed #16.0%

Vegetation <- (sum(Cont$Veg == "Y") / 1173) * 100 #181
Vegetation #15.4%

Percent <- c(Fruit, Anthropogenic, Garbage, NaturalPrey, Birdseed, Vegetation)
Item <- c("Fruit", "Anthropogenic", "Garbage", "Natural Prey", "Birdseed", "Vegetation")
ContChart1 <- data.frame(Percent, Item)

ContChart1$Item <- ifelse(ContChart1$Item == "Natural Prey", "Natural prey", ContChart1$Item)

#Make a Bar chart containing percent of scats containing each type of content
ContentChart1 <- ggplot(data = ContChart1, mapping = aes(x = Item, y = Percent, fill = Item)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = c("Fruit", "Anthropogenic", "Birdseed", "Natural prey", "Vegetation", "Garbage")) +
  scale_fill_manual(values=c("#00695C", #Anthro
                             "#00796B", #Birdseed
                             "#004D40", #Fruitc("#004D40", "#00695C","#00796B","#00897B","#009688", "#26A69A")
                             "#26A69A", #Garbage
                             "#00897B", #Nat Prey
                             "#009688")) + #vege
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(x="Content item", y="Percent of scats containing item")

ContentChart1

#Save chart
ggsave('ContentChart.png', ContentChart1, dpi = 700, width = 5, height = 4)

###Use logistic regression to compare scats that contain a target thing (1) and those that don't (0)
ContF <- Cont %>%
  dplyr::select(-c(Season, Location, Size, Freshness, Magpie,
                   Unk, EdDens, DistEdge, LC, Aspect, AG, WATER,
                   ScatCount, Ed_Len, NatEdLen, DistSWater, DistRiv,
                   DistGarden, North, AspRad))

ContF <- ContF %>%
  mutate(Fruit = ifelse(Fruit == "N", 0, 1)) %>%
  mutate(Anthro = ifelse(Anthro == "N", 0, 1)) %>%
  mutate(Garbage = ifelse(Garbage == "N", 0, 1)) %>%
  mutate(NatPrey = ifelse(NatPrey == "N", 0, 1)) %>%
  mutate(Seed = ifelse(Seed == "N", 0, 1)) %>%
  mutate(Veg = ifelse(Veg == "N", 0, 1))

ContF$Fruit = factor(ContF$Fruit, levels = c(0,1))
ContF$Anthro = factor(ContF$Anthro, levels = c(0,1))
ContF$Garbage = factor(ContF$Garbage, levels = c(0,1))
ContF$NatPrey = factor(ContF$NatPrey, levels = c(0,1))
ContF$Seed = factor(ContF$Seed, levels = c(0,1))
ContF$Veg = factor(ContF$Veg, levels = c(0,1))

#Make decay terms, MANY, for each distance metric. Then choose the best one
ContF <- ContF %>%
  mutate(Dec_PG0.002 = 1-(exp(-0.002*DistPG))) %>%
  mutate(Dec_PG0.004 = 1-(exp(-0.004*DistPG))) %>%
  mutate(Dec_PG0.006 = 1-(exp(-0.006*DistPG))) %>%
  mutate(Dec_PG0.012 = 1-(exp(-0.012*DistPG))) %>%
  mutate(Dec_PG0.03 = 1-(exp(-0.03*DistPG))) %>%
  mutate(Dec_PG0.06 = 1-(exp(-0.06*DistPG))) %>%
  mutate(Dec_PG0.2 = 1-(exp(-0.2*DistPG))) %>%
  mutate(Dec_Water0.002 = 1-(exp(-0.002*DistWater))) %>%
  mutate(Dec_Water0.004 = 1-(exp(-0.004*DistWater))) %>%
  mutate(Dec_Water0.006 = 1-(exp(-0.006*DistWater))) %>%
  mutate(Dec_Water0.012 = 1-(exp(-0.012*DistWater))) %>%
  mutate(Dec_Water0.03 = 1-(exp(-0.03*DistWater))) %>%
  mutate(Dec_Water0.06 = 1-(exp(-0.06*DistWater))) %>%
  mutate(Dec_Water0.2 = 1-(exp(-0.2*DistWater))) %>%
  mutate(Dec_Bldg0.002 = 1-(exp(-0.002*DistBldg))) %>%
  mutate(Dec_Bldg0.004 = 1-(exp(-0.004*DistBldg))) %>%
  mutate(Dec_Bldg0.006 = 1-(exp(-0.006*DistBldg))) %>%
  mutate(Dec_Bldg0.012 = 1-(exp(-0.012*DistBldg))) %>%
  mutate(Dec_Bldg0.03 = 1-(exp(-0.03*DistBldg))) %>%
  mutate(Dec_Bldg0.06 = 1-(exp(-0.06*DistBldg))) %>%
  mutate(Dec_Bldg0.2 = 1-(exp(-0.2*DistBldg))) %>%
  mutate(Dec_Road0.002 = 1-(exp(-0.002*DistRoad))) %>%
  mutate(Dec_Road0.004 = 1-(exp(-0.004*DistRoad))) %>%
  mutate(Dec_Road0.006 = 1-(exp(-0.006*DistRoad))) %>%
  mutate(Dec_Road0.012 = 1-(exp(-0.012*DistRoad))) %>%
  mutate(Dec_Road0.03 = 1-(exp(-0.03*DistRoad))) %>%
  mutate(Dec_Road0.06 = 1-(exp(-0.06*DistRoad))) %>%
  mutate(Dec_Road0.2 = 1-(exp(-0.2*DistRoad))) %>%
  mutate(Dec_NatEdge0.002 = 1-(exp(-0.002*DistNatEdge))) %>%
  mutate(Dec_NatEdge0.004 = 1-(exp(-0.004*DistNatEdge))) %>%
  mutate(Dec_NatEdge0.006 = 1-(exp(-0.006*DistNatEdge))) %>%
  mutate(Dec_NatEdge0.012 = 1-(exp(-0.012*DistNatEdge))) %>%
  mutate(Dec_NatEdge0.03 = 1-(exp(-0.03*DistNatEdge))) %>%
  mutate(Dec_NatEdge0.06 = 1-(exp(-0.06*DistNatEdge))) %>%
  mutate(Dec_NatEdge0.2 = 1-(exp(-0.2*DistNatEdge))) %>%
  mutate(Dec_Scat0.002 = 1-(exp(-0.002*DistScat))) %>%
  mutate(Dec_Scat0.004 = 1-(exp(-0.004*DistScat))) %>%
  mutate(Dec_Scat0.006 = 1-(exp(-0.006*DistScat))) %>%
  mutate(Dec_Scat0.012 = 1-(exp(-0.012*DistScat))) %>%
  mutate(Dec_Scat0.03 = 1-(exp(-0.03*DistScat))) %>%
  mutate(Dec_Scat0.06 = 1-(exp(-0.06*DistScat))) %>%
  mutate(Dec_Scat0.2 = 1-(exp(-0.2*DistScat))) %>%
  mutate(Dec_Junc0.002 = 1-(exp(-0.002*DistJunc))) %>%
  mutate(Dec_Junc0.004 = 1-(exp(-0.004*DistJunc))) %>%
  mutate(Dec_Junc0.006 = 1-(exp(-0.006*DistJunc))) %>%
  mutate(Dec_Junc0.012 = 1-(exp(-0.012*DistJunc))) %>%
  mutate(Dec_Junc0.03 = 1-(exp(-0.03*DistJunc))) %>%
  mutate(Dec_Junc0.06 = 1-(exp(-0.06*DistJunc))) %>%
  mutate(Dec_Junc0.2 = 1-(exp(-0.2*DistJunc))) %>%
  mutate(Dec_Trail0.002 = 1-(exp(-0.002*DistTrail))) %>%
  mutate(Dec_Trail0.004 = 1-(exp(-0.004*DistTrail))) %>%
  mutate(Dec_Trail0.006 = 1-(exp(-0.006*DistTrail))) %>%
  mutate(Dec_Trail0.012 = 1-(exp(-0.012*DistTrail))) %>%
  mutate(Dec_Trail0.03 = 1-(exp(-0.03*DistTrail))) %>%
  mutate(Dec_Trail0.06 = 1-(exp(-0.06*DistTrail))) %>%
  mutate(Dec_Trail0.2 = 1-(exp(-0.2*DistTrail))) %>%
  mutate(Dec_Camp0.002 = 1-(exp(-0.002*DistCamp))) %>%
  mutate(Dec_Camp0.004 = 1-(exp(-0.004*DistCamp))) %>%
  mutate(Dec_Camp0.006 = 1-(exp(-0.006*DistCamp))) %>%
  mutate(Dec_Camp0.012 = 1-(exp(-0.012*DistCamp))) %>%
  mutate(Dec_Camp0.03 = 1-(exp(-0.03*DistCamp))) %>%
  mutate(Dec_Camp0.06 = 1-(exp(-0.06*DistCamp))) %>%
  mutate(Dec_Camp0.2 = 1-(exp(-0.2*DistCamp)))

Decmods <- list()
for(i in c(1:9)) {
  if(i==1){decs <- c("DistWater", "Dec_Water0.002", "Dec_Water0.004", "Dec_Water0.006", "Dec_Water0.012", "Dec_Water0.03", "Dec_Water0.06", "Dec_Water0.2")}
  if(i==2){decs <- c("DistBldg", "Dec_Bldg0.002", "Dec_Bldg0.004", "Dec_Bldg0.006", "Dec_Bldg0.012", "Dec_Bldg0.03", "Dec_Bldg0.06", "Dec_Bldg0.2")}
  if(i==3){decs <- c("DistRoad", "Dec_Road0.002", "Dec_Road0.004", "Dec_Road0.006", "Dec_Road0.012", "Dec_Road0.03", "Dec_Road0.06", "Dec_Road0.2")}
  if(i==4){decs <- c("DistNatEdge", "Dec_NatEdge0.002", "Dec_NatEdge0.004", "Dec_NatEdge0.006", "Dec_NatEdge0.012", "Dec_NatEdge0.03", "Dec_NatEdge0.06", "Dec_NatEdge0.2")}
  if(i==5){decs <- c("DistScat", "Dec_Scat0.002", "Dec_Scat0.004", "Dec_Scat0.006", "Dec_Scat0.012", "Dec_Scat0.03", "Dec_Scat0.06", "Dec_Scat0.2")}
  if(i==6){decs <- c("DistJunc", "Dec_Junc0.002", "Dec_Junc0.004", "Dec_Junc0.006", "Dec_Junc0.012", "Dec_Junc0.03", "Dec_Junc0.06", "Dec_Junc0.2")}
  if(i==7){decs <- c("DistTrail", "Dec_Trail0.002", "Dec_Trail0.004", "Dec_Trail0.006", "Dec_Trail0.012", "Dec_Trail0.03", "Dec_Trail0.06", "Dec_Trail0.2")}
  if(i==8){decs <- c("DistCamp", "Dec_Camp0.002", "Dec_Camp0.004", "Dec_Camp0.006", "Dec_Camp0.012", "Dec_Camp0.03", "Dec_Camp0.06", "Dec_Camp0.2")}
  if(i==9){decs <- c("DistPG", "Dec_PG0.002", "Dec_PG0.004", "Dec_PG0.006", "Dec_PG0.012", "Dec_PG0.03", "Dec_PG0.06", "Dec_PG0.2")}
  
  
  # Create a data frame to store the results from all the models.
  Decmods[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(decs)){
    
    model <- glm(Fruit ~ ContF[,j], # Create the predictive model
                 data=ContF, family = binomial)
    
    null <- glm(Fruit ~1, # Create the null model
                data=ContF, family = binomial)
    
    
    Decmods[[i]][j,1] <- j # First column is the predictor
    Decmods[[i]][j,2] <- AIC(model) # Second column is AIC
  }
  
  colnames(Decmods[[i]]) <- c("predictor", "AIC")
}
DecFruit <- dplyr::bind_rows(Decmods)
DecFruit$Cont <- "Fruit"



Decmods <- list()
for(i in c(1:9)) {
  if(i==1){decs <- c("DistWater", "Dec_Water0.002", "Dec_Water0.004", "Dec_Water0.006", "Dec_Water0.012", "Dec_Water0.03", "Dec_Water0.06", "Dec_Water0.2")}
  if(i==2){decs <- c("DistBldg", "Dec_Bldg0.002", "Dec_Bldg0.004", "Dec_Bldg0.006", "Dec_Bldg0.012", "Dec_Bldg0.03", "Dec_Bldg0.06", "Dec_Bldg0.2")}
  if(i==3){decs <- c("DistRoad", "Dec_Road0.002", "Dec_Road0.004", "Dec_Road0.006", "Dec_Road0.012", "Dec_Road0.03", "Dec_Road0.06", "Dec_Road0.2")}
  if(i==4){decs <- c("DistNatEdge", "Dec_NatEdge0.002", "Dec_NatEdge0.004", "Dec_NatEdge0.006", "Dec_NatEdge0.012", "Dec_NatEdge0.03", "Dec_NatEdge0.06", "Dec_NatEdge0.2")}
  if(i==5){decs <- c("DistScat", "Dec_Scat0.002", "Dec_Scat0.004", "Dec_Scat0.006", "Dec_Scat0.012", "Dec_Scat0.03", "Dec_Scat0.06", "Dec_Scat0.2")}
  if(i==6){decs <- c("DistJunc", "Dec_Junc0.002", "Dec_Junc0.004", "Dec_Junc0.006", "Dec_Junc0.012", "Dec_Junc0.03", "Dec_Junc0.06", "Dec_Junc0.2")}
  if(i==7){decs <- c("DistTrail", "Dec_Trail0.002", "Dec_Trail0.004", "Dec_Trail0.006", "Dec_Trail0.012", "Dec_Trail0.03", "Dec_Trail0.06", "Dec_Trail0.2")}
  if(i==8){decs <- c("DistCamp", "Dec_Camp0.002", "Dec_Camp0.004", "Dec_Camp0.006", "Dec_Camp0.012", "Dec_Camp0.03", "Dec_Camp0.06", "Dec_Camp0.2")}
  if(i==9){decs <- c("DistPG", "Dec_PG0.002", "Dec_PG0.004", "Dec_PG0.006", "Dec_PG0.012", "Dec_PG0.03", "Dec_PG0.06", "Dec_PG0.2")}
  
  
  # Create a data frame to store the results from all the models.
  Decmods[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(decs)){
    
    model <- glm(Anthro ~ ContF[,j], # Create the predictive model
                 data=ContF, family = binomial)
    
    null <- glm(Anthro ~1, # Create the null model
                data=ContF, family = binomial)
    
    
    Decmods[[i]][j,1] <- j # First column is the predictor
    Decmods[[i]][j,2] <- AIC(model) # Second column is AIC
  }
  
  colnames(Decmods[[i]]) <- c("predictor", "AIC")
}
DecAnthro <- dplyr::bind_rows(Decmods)
DecAnthro$Cont <- "Anthro"


Decmods <- list()
for(i in c(1:9)) {
  if(i==1){decs <- c("DistWater", "Dec_Water0.002", "Dec_Water0.004", "Dec_Water0.006", "Dec_Water0.012", "Dec_Water0.03", "Dec_Water0.06", "Dec_Water0.2")}
  if(i==2){decs <- c("DistBldg", "Dec_Bldg0.002", "Dec_Bldg0.004", "Dec_Bldg0.006", "Dec_Bldg0.012", "Dec_Bldg0.03", "Dec_Bldg0.06", "Dec_Bldg0.2")}
  if(i==3){decs <- c("DistRoad", "Dec_Road0.002", "Dec_Road0.004", "Dec_Road0.006", "Dec_Road0.012", "Dec_Road0.03", "Dec_Road0.06", "Dec_Road0.2")}
  if(i==4){decs <- c("DistNatEdge", "Dec_NatEdge0.002", "Dec_NatEdge0.004", "Dec_NatEdge0.006", "Dec_NatEdge0.012", "Dec_NatEdge0.03", "Dec_NatEdge0.06", "Dec_NatEdge0.2")}
  if(i==5){decs <- c("DistScat", "Dec_Scat0.002", "Dec_Scat0.004", "Dec_Scat0.006", "Dec_Scat0.012", "Dec_Scat0.03", "Dec_Scat0.06", "Dec_Scat0.2")}
  if(i==6){decs <- c("DistJunc", "Dec_Junc0.002", "Dec_Junc0.004", "Dec_Junc0.006", "Dec_Junc0.012", "Dec_Junc0.03", "Dec_Junc0.06", "Dec_Junc0.2")}
  if(i==7){decs <- c("DistTrail", "Dec_Trail0.002", "Dec_Trail0.004", "Dec_Trail0.006", "Dec_Trail0.012", "Dec_Trail0.03", "Dec_Trail0.06", "Dec_Trail0.2")}
  if(i==8){decs <- c("DistCamp", "Dec_Camp0.002", "Dec_Camp0.004", "Dec_Camp0.006", "Dec_Camp0.012", "Dec_Camp0.03", "Dec_Camp0.06", "Dec_Camp0.2")}
  if(i==9){decs <- c("DistPG", "Dec_PG0.002", "Dec_PG0.004", "Dec_PG0.006", "Dec_PG0.012", "Dec_PG0.03", "Dec_PG0.06", "Dec_PG0.2")}
  
  
  # Create a data frame to store the results from all the models.
  Decmods[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(decs)){
    
    model <- glm(NatPrey ~ ContF[,j], # Create the predictive model
                 data=ContF, family = binomial)
    
    null <- glm(NatPrey ~1, # Create the null model
                data=ContF, family = binomial)
    
    
    Decmods[[i]][j,1] <- j # First column is the predictor
    Decmods[[i]][j,2] <- AIC(model) # Second column is AIC
  }
  
  colnames(Decmods[[i]]) <- c("predictor", "AIC")
}
DecNatPrey <- dplyr::bind_rows(Decmods)
DecNatPrey$Cont <- "NatPrey"


Decmods <- list()
for(i in c(1:9)) {
  if(i==1){decs <- c("DistWater", "Dec_Water0.002", "Dec_Water0.004", "Dec_Water0.006", "Dec_Water0.012", "Dec_Water0.03", "Dec_Water0.06", "Dec_Water0.2")}
  if(i==2){decs <- c("DistBldg", "Dec_Bldg0.002", "Dec_Bldg0.004", "Dec_Bldg0.006", "Dec_Bldg0.012", "Dec_Bldg0.03", "Dec_Bldg0.06", "Dec_Bldg0.2")}
  if(i==3){decs <- c("DistRoad", "Dec_Road0.002", "Dec_Road0.004", "Dec_Road0.006", "Dec_Road0.012", "Dec_Road0.03", "Dec_Road0.06", "Dec_Road0.2")}
  if(i==4){decs <- c("DistNatEdge", "Dec_NatEdge0.002", "Dec_NatEdge0.004", "Dec_NatEdge0.006", "Dec_NatEdge0.012", "Dec_NatEdge0.03", "Dec_NatEdge0.06", "Dec_NatEdge0.2")}
  if(i==5){decs <- c("DistScat", "Dec_Scat0.002", "Dec_Scat0.004", "Dec_Scat0.006", "Dec_Scat0.012", "Dec_Scat0.03", "Dec_Scat0.06", "Dec_Scat0.2")}
  if(i==6){decs <- c("DistJunc", "Dec_Junc0.002", "Dec_Junc0.004", "Dec_Junc0.006", "Dec_Junc0.012", "Dec_Junc0.03", "Dec_Junc0.06", "Dec_Junc0.2")}
  if(i==7){decs <- c("DistTrail", "Dec_Trail0.002", "Dec_Trail0.004", "Dec_Trail0.006", "Dec_Trail0.012", "Dec_Trail0.03", "Dec_Trail0.06", "Dec_Trail0.2")}
  if(i==8){decs <- c("DistCamp", "Dec_Camp0.002", "Dec_Camp0.004", "Dec_Camp0.006", "Dec_Camp0.012", "Dec_Camp0.03", "Dec_Camp0.06", "Dec_Camp0.2")}
  if(i==9){decs <- c("DistPG", "Dec_PG0.002", "Dec_PG0.004", "Dec_PG0.006", "Dec_PG0.012", "Dec_PG0.03", "Dec_PG0.06", "Dec_PG0.2")}
  
  
  # Create a data frame to store the results from all the models.
  Decmods[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(decs)){
    
    model <- glm(Seed ~ ContF[,j], # Create the predictive model
                 data=ContF, family = binomial)
    
    null <- glm(Seed ~1, # Create the null model
                data=ContF, family = binomial)
    
    
    Decmods[[i]][j,1] <- j # First column is the predictor
    Decmods[[i]][j,2] <- AIC(model) # Second column is AIC
  }
  
  colnames(Decmods[[i]]) <- c("predictor", "AIC")
}
DecSeed <- dplyr::bind_rows(Decmods)
DecSeed$Cont <- "Seed"


Decmods <- list()
for(i in c(1:9)) {
  if(i==1){decs <- c("DistWater", "Dec_Water0.002", "Dec_Water0.004", "Dec_Water0.006", "Dec_Water0.012", "Dec_Water0.03", "Dec_Water0.06", "Dec_Water0.2")}
  if(i==2){decs <- c("DistBldg", "Dec_Bldg0.002", "Dec_Bldg0.004", "Dec_Bldg0.006", "Dec_Bldg0.012", "Dec_Bldg0.03", "Dec_Bldg0.06", "Dec_Bldg0.2")}
  if(i==3){decs <- c("DistRoad", "Dec_Road0.002", "Dec_Road0.004", "Dec_Road0.006", "Dec_Road0.012", "Dec_Road0.03", "Dec_Road0.06", "Dec_Road0.2")}
  if(i==4){decs <- c("DistNatEdge", "Dec_NatEdge0.002", "Dec_NatEdge0.004", "Dec_NatEdge0.006", "Dec_NatEdge0.012", "Dec_NatEdge0.03", "Dec_NatEdge0.06", "Dec_NatEdge0.2")}
  if(i==5){decs <- c("DistScat", "Dec_Scat0.002", "Dec_Scat0.004", "Dec_Scat0.006", "Dec_Scat0.012", "Dec_Scat0.03", "Dec_Scat0.06", "Dec_Scat0.2")}
  if(i==6){decs <- c("DistJunc", "Dec_Junc0.002", "Dec_Junc0.004", "Dec_Junc0.006", "Dec_Junc0.012", "Dec_Junc0.03", "Dec_Junc0.06", "Dec_Junc0.2")}
  if(i==7){decs <- c("DistTrail", "Dec_Trail0.002", "Dec_Trail0.004", "Dec_Trail0.006", "Dec_Trail0.012", "Dec_Trail0.03", "Dec_Trail0.06", "Dec_Trail0.2")}
  if(i==8){decs <- c("DistCamp", "Dec_Camp0.002", "Dec_Camp0.004", "Dec_Camp0.006", "Dec_Camp0.012", "Dec_Camp0.03", "Dec_Camp0.06", "Dec_Camp0.2")}
  if(i==9){decs <- c("DistPG", "Dec_PG0.002", "Dec_PG0.004", "Dec_PG0.006", "Dec_PG0.012", "Dec_PG0.03", "Dec_PG0.06", "Dec_PG0.2")}
  
  
  # Create a data frame to store the results from all the models.
  Decmods[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(decs)){
    
    model <- glm(Veg ~ ContF[,j], # Create the predictive model
                 data=ContF, family = binomial)
    
    null <- glm(Veg ~1, # Create the null model
                data=ContF, family = binomial)
    
    
    Decmods[[i]][j,1] <- j # First column is the predictor
    Decmods[[i]][j,2] <- AIC(model) # Second column is AIC
  }
  
  colnames(Decmods[[i]]) <- c("predictor", "AIC")
}
DecVeg <- dplyr::bind_rows(Decmods)
DecVeg$Cont <- "Veg"

Decmods <- list()
for(i in c(1:9)) {
  if(i==1){decs <- c("DistWater", "Dec_Water0.002", "Dec_Water0.004", "Dec_Water0.006", "Dec_Water0.012", "Dec_Water0.03", "Dec_Water0.06", "Dec_Water0.2")}
  if(i==2){decs <- c("DistBldg", "Dec_Bldg0.002", "Dec_Bldg0.004", "Dec_Bldg0.006", "Dec_Bldg0.012", "Dec_Bldg0.03", "Dec_Bldg0.06", "Dec_Bldg0.2")}
  if(i==3){decs <- c("DistRoad", "Dec_Road0.002", "Dec_Road0.004", "Dec_Road0.006", "Dec_Road0.012", "Dec_Road0.03", "Dec_Road0.06", "Dec_Road0.2")}
  if(i==4){decs <- c("DistNatEdge", "Dec_NatEdge0.002", "Dec_NatEdge0.004", "Dec_NatEdge0.006", "Dec_NatEdge0.012", "Dec_NatEdge0.03", "Dec_NatEdge0.06", "Dec_NatEdge0.2")}
  if(i==5){decs <- c("DistScat", "Dec_Scat0.002", "Dec_Scat0.004", "Dec_Scat0.006", "Dec_Scat0.012", "Dec_Scat0.03", "Dec_Scat0.06", "Dec_Scat0.2")}
  if(i==6){decs <- c("DistJunc", "Dec_Junc0.002", "Dec_Junc0.004", "Dec_Junc0.006", "Dec_Junc0.012", "Dec_Junc0.03", "Dec_Junc0.06", "Dec_Junc0.2")}
  if(i==7){decs <- c("DistTrail", "Dec_Trail0.002", "Dec_Trail0.004", "Dec_Trail0.006", "Dec_Trail0.012", "Dec_Trail0.03", "Dec_Trail0.06", "Dec_Trail0.2")}
  if(i==8){decs <- c("DistCamp", "Dec_Camp0.002", "Dec_Camp0.004", "Dec_Camp0.006", "Dec_Camp0.012", "Dec_Camp0.03", "Dec_Camp0.06", "Dec_Camp0.2")}
  if(i==9){decs <- c("DistPG", "Dec_PG0.002", "Dec_PG0.004", "Dec_PG0.006", "Dec_PG0.012", "Dec_PG0.03", "Dec_PG0.06", "Dec_PG0.2")}
  
  
  # Create a data frame to store the results from all the models.
  Decmods[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(decs)){
    
    model <- glm(Garbage ~ ContF[,j], # Create the predictive model
                 data=ContF, family = binomial)
    
    null <- glm(Garbage ~1, # Create the null model
                data=ContF, family = binomial)
    
    
    Decmods[[i]][j,1] <- j # First column is the predictor
    Decmods[[i]][j,2] <- AIC(model) # Second column is AIC
  }
  
  colnames(Decmods[[i]]) <- c("predictor", "AIC")
}
DecGarbage <- dplyr::bind_rows(Decmods)
DecGarbage$Cont <- "Garbage"

DecayAICresults <- rbind(DecFruit, DecAnthro, DecSeed, DecNatPrey, DecVeg, DecGarbage)
DecayAICresults2 <- t(DecayAICresults)

#Create a table that shows the AIC for all models with all decay terms to facilitate selecting the right ones
write.csv(DecayAICresults2, file = "DecayAICresults.csv")

#Scale and center everything in the dataframe
ContF[7:87] <- scale(ContF[7:87], scale = TRUE, center = TRUE)


###FRUIT
ContF_Fruit <- ContF %>%
  dplyr::select(Fruit, DistWater, DistBldg, DistRoad, DistNatEdge, DistScat, 
                NatEdDens, RoadDens, BldgDens, Curve1, Slope, ANTH, GRASS, NAT,
                DistJunc, DistTrail, DistCamp, DistPG, East, Dec_Water0.2,
                Dec_Bldg0.2, Dec_Road0.2, Dec_NatEdge0.03, Dec_Scat0.03, Dec_Camp0.004,
                Dec_PG0.06) %>%
  dplyr::rename("Dec_Water" = Dec_Water0.2,
                "Dec_Bldg" = Dec_Bldg0.2,
                "Dec_Road" = Dec_Road0.2,
                "Dec_NatEdge" = Dec_NatEdge0.03,
                "Dec_Scat" = Dec_Scat0.03,
                "Dec_Camp" = Dec_Camp0.004, 
                "Dec_PG" = Dec_PG0.06) %>%
  dplyr::select(Fruit, )

#Get univariate results for FRUIT
#start by getting linear model info in a table
Fruitmods_lin <- list()
for(i in c(1:14)) {
  if(i==1){fruits <- ("DistBldg")}
  if(i==2){fruits <- ("DistRoad")}
  if(i==3){fruits <- ("DistNatEdge")}
  if(i==4){fruits <- ("DistScat")}
  if(i==5){fruits <- ("DistTrail")}
  if(i==6){fruits <- ("DistCamp")}
  if(i==7){fruits <- ("DistPG")}
  if(i==8){fruits <- ("NatEdDens")}
  if(i==9){fruits <- ("RoadDens")}
  if(i==10){fruits <- ("BldgDens")}
  if(i==11){fruits <- ("DistNatEdge")}
  if(i==12){fruits <- ("ANTH")}
  if(i==13){fruits <- ("GRASS")}
  if(i==14){fruits <- ("NAT")}

  # Create a data frame to store the results from all the models.
  Fruitmods_lin[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(fruits)){
    
    model <- glm(Fruit ~ ContF_Fruit[,j], # Create the predictive model
                 data=ContF_Fruit, family = binomial)
    
    Fruitmods_lin[[i]][j,1] <- j # First column is the predictor
    Fruitmods_lin[[i]][j,2] <- AIC(model) # Second column is AIC
    Fruitmods_lin[[i]][j,3:4] <- model$coefficients # thirs and fourth column are beta for intercept and predictor
    Fruitmods_lin[[i]][j,5:6] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Fruitmods_lin[[i]][j,7:10] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Fruitmods_lin[[i]]) <- c("predictor", "AIC", "Intercept_Beta", "Predictor_beta", "Intercept_pval", "Predictor_pval", "Intercept_LowerCI", "Predictor_LowerCI", "Intercept_UpperCI", "Predictor_UpperCI")
}
FruitUnivariates_lin <- dplyr::bind_rows(Fruitmods_lin)

FruitUnivariates_lin <- FruitUnivariates_lin %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictor_pval, Predictor_LowerCI, Predictor_UpperCI) %>%
  dplyr::rename("Beta" = Predictor_beta,
                "pvalue" = Predictor_pval,
                "LowCI" = Predictor_LowerCI,
                "UppCI" = Predictor_UpperCI)

#Now repeat the process for variables modelled as quadratics
Fruitmods_quad <- list()
for(i in c(1:14)) {
  if(i==1){fruits <- ("DistBldg")}
  if(i==2){fruits <- ("DistRoad")}
  if(i==3){fruits <- ("DistNatEdge")}
  if(i==4){fruits <- ("DistScat")}
  if(i==5){fruits <- ("DistTrail")}
  if(i==6){fruits <- ("DistCamp")}
  if(i==7){fruits <- ("DistPG")}
  if(i==8){fruits <- ("NatEdDens")}
  if(i==9){fruits <- ("RoadDens")}
  if(i==10){fruits <- ("BldgDens")}
  if(i==11){fruits <- ("DistNatEdge")}
  if(i==12){fruits <- ("ANTH")}
  if(i==13){fruits <- ("GRASS")}
  if(i==14){fruits <- ("NAT")}
  
  # Create a data frame to store the results from all the models.
  Fruitmods_quad[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(fruits)){
    
    model <- glm(Fruit ~ ContF_Fruit[,j] + I(ContF_Fruit[,j]^2), # Create the predictive model
                 data=ContF_Fruit, family = binomial)
    
    Fruitmods_quad[[i]][j,1] <- j # First column is the predictor
    Fruitmods_quad[[i]][j,2] <- AIC(model) # Second column is AIC
    Fruitmods_quad[[i]][j,3:5] <- model$coefficients # thirs, fourth,fifth column are beta for intercept and predictor, predictor^2
    Fruitmods_quad[[i]][j,6:8] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Fruitmods_quad[[i]][j,9:14] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Fruitmods_quad[[i]]) <- c("predictor", "AIC", "Intercept_Beta", 
                                     "Predictor_beta", "Predictore_beta2", 
                                     "Intercept_pval", "Predictor_pval", 
                                     "Predictor_pval2", "Intercept_LowerCI", 
                                     "Predictor_LowerCI", "Predictor_LowerCI2",
                                     "Intercept_UpperCI", "Predictor_UpperCI",
                                     "Predictor_UpperCI2")
}
FruitUnivariates_quad <- dplyr::bind_rows(Fruitmods_quad)

FruitUnivariates_quad <- FruitUnivariates_quad %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictore_beta2,
                Predictor_pval, Predictor_pval2, Predictor_LowerCI, Predictor_UpperCI,
                Predictor_LowerCI2, Predictor_UpperCI2) %>%
  dplyr::rename("Q_Beta(lin)" = Predictor_beta,
                "Q_Beta(quad)" = Predictore_beta2,
                "Q_pvalue(lin)" = Predictor_pval,
                "Q_pvalue(quad)" = Predictor_pval2,
                "Q_LowCI(lin)" = Predictor_LowerCI,
                "Q_UppCI(lin)" = Predictor_UpperCI,
                "Q_LowCI(quad)" = Predictor_LowerCI2,
                "Q_UppCI(quad)" = Predictor_UpperCI2,
                "AIC_quad" = AIC)

FruitUnivariates <- cbind(FruitUnivariates_lin, FruitUnivariates_quad)
write.csv(FruitUnivariates, file = "FruitUnivariatesMarch2023.csv")

#ANTHROPOGENIC 
ContF_Anthro <- ContF %>%
  dplyr::select(Anthro, DistBldg, DistRoad, DistNatEdge, DistScat, 
                NatEdDens, RoadDens, BldgDens, NAT, ANTH, GRASS,
                DistTrail, DistCamp, DistPG)

#Get univariates results for ANTHRO
#start by getting linear model info in a table
Anthromods_lin <- list()
for(i in c(1:14)) {
  if(i==1){anthros <- ("DistBldg")}
  if(i==2){anthros <- ("DistRoad")}
  if(i==3){anthros <- ("DistNatEdge")}
  if(i==4){anthros <- ("DistScat")}
  if(i==5){anthros <- ("DistTrail")}
  if(i==6){anthros <- ("DistCamp")}
  if(i==7){anthros <- ("DistPG")}
  if(i==8){anthros <- ("NatEdDens")}
  if(i==9){anthros <- ("RoadDens")}
  if(i==10){anthros <- ("BldgDens")}
  if(i==11){anthros <- ("DistNatEdge")}
  if(i==12){anthros <- ("ANTH")}
  if(i==13){anthros <- ("GRASS")}
  if(i==14){anthros <- ("NAT")}
  
  
  # Create a data frame to store the results from all the models.
  Anthromods_lin[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(anthros)){
    
    model <- glm(Anthro ~ ContF_Anthro[,j], # Create the predictive model
                 data=ContF_Anthro, family = binomial)
    
    Anthromods_lin[[i]][j,1] <- j # First column is the predictor
    Anthromods_lin[[i]][j,2] <- AIC(model) # Second column is AIC
    Anthromods_lin[[i]][j,3:4] <- model$coefficients # thirs and fourth column are beta for intercept and predictor
    Anthromods_lin[[i]][j,5:6] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Anthromods_lin[[i]][j,7:10] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Anthromods_lin[[i]]) <- c("predictor", "AIC", "Intercept_Beta", "Predictor_beta", "Intercept_pval", "Predictor_pval", "Intercept_LowerCI", "Predictor_LowerCI", "Intercept_UpperCI", "Predictor_UpperCI")
}
AnthroUnivariates_lin <- dplyr::bind_rows(Anthromods_lin)

AnthroUnivariates_lin <- AnthroUnivariates_lin %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictor_pval, Predictor_LowerCI, Predictor_UpperCI) %>%
  dplyr::rename("Beta" = Predictor_beta,
                "pvalue" = Predictor_pval,
                "LowCI" = Predictor_LowerCI,
                "UppCI" = Predictor_UpperCI)

#Now repeat the process for variables modelled as quadratics
Anthromods_quad <- list()
for(i in c(1:14)) {
  if(i==1){anthros <- ("DistBldg")}
  if(i==2){anthros <- ("DistRoad")}
  if(i==3){anthros <- ("DistNatEdge")}
  if(i==4){anthros <- ("DistScat")}
  if(i==5){anthros <- ("DistTrail")}
  if(i==6){anthros <- ("DistCamp")}
  if(i==7){anthros <- ("DistPG")}
  if(i==8){anthros <- ("NatEdDens")}
  if(i==9){anthros <- ("RoadDens")}
  if(i==10){anthros <- ("BldgDens")}
  if(i==11){anthros <- ("DistNatEdge")}
  if(i==12){anthros <- ("ANTH")}
  if(i==13){anthros <- ("GRASS")}
  if(i==14){anthros <- ("NAT")}
  
  # Create a data frame to store the results from all the models.
  Anthromods_quad[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(anthros)){
    
    model <- glm(Anthro ~ ContF_Anthro[,j] + I(ContF_Anthro[,j]^2), # Create the predictive model
                 data=ContF_Anthro, family = binomial)
    
    Anthromods_quad[[i]][j,1] <- j # First column is the predictor
    Anthromods_quad[[i]][j,2] <- AIC(model) # Second column is AIC
    Anthromods_quad[[i]][j,3:5] <- model$coefficients # thirs, fourth,fifth column are beta for intercept and predictor, predictor^2
    Anthromods_quad[[i]][j,6:8] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Anthromods_quad[[i]][j,9:14] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Anthromods_quad[[i]]) <- c("predictor", "AIC", "Intercept_Beta", 
                                      "Predictor_beta", "Predictore_beta2", 
                                      "Intercept_pval", "Predictor_pval", 
                                      "Predictor_pval2", "Intercept_LowerCI", 
                                      "Predictor_LowerCI", "Predictor_LowerCI2",
                                      "Intercept_UpperCI", "Predictor_UpperCI",
                                      "Predictor_UpperCI2")
}
AnthroUnivariates_quad <- dplyr::bind_rows(Anthromods_quad)

AnthroUnivariates_quad <- AnthroUnivariates_quad %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictore_beta2,
                Predictor_pval, Predictor_pval2, Predictor_LowerCI, Predictor_UpperCI,
                Predictor_LowerCI2, Predictor_UpperCI2) %>%
  dplyr::rename("Q_Beta(lin)" = Predictor_beta,
                "Q_Beta(quad)" = Predictore_beta2,
                "Q_pvalue(lin)" = Predictor_pval,
                "Q_pvalue(quad)" = Predictor_pval2,
                "Q_LowCI(lin)" = Predictor_LowerCI,
                "Q_UppCI(lin)" = Predictor_UpperCI,
                "Q_LowCI(quad)" = Predictor_LowerCI2,
                "Q_UppCI(quad)" = Predictor_UpperCI2,
                "AIC_quad" = AIC)

AnthroUnivariates <- cbind(AnthroUnivariates_lin, AnthroUnivariates_quad)
write.csv(AnthroUnivariates, file = "AnthroUnivariatesMarch2023.csv")

###BIRDSEED###
ContF_Seed <- ContF %>%
  dplyr::select(Seed, DistBldg, DistRoad, DistNatEdge, DistScat, 
                NatEdDens, RoadDens, BldgDens, NAT, ANTH, GRASS,
                DistTrail, DistCamp, DistPG)

#start by getting linear model info in a table
Seedmods_lin <- list()
for(i in c(1:14)) {
  if(i==1){seeds <- ("DistBldg")}
  if(i==2){seeds <- ("DistRoad")}
  if(i==3){seeds <- ("DistNatEdge")}
  if(i==4){seeds <- ("DistScat")}
  if(i==5){seeds <- ("DistTrail")}
  if(i==6){seeds <- ("DistCamp")}
  if(i==7){seeds <- ("DistPG")}
  if(i==8){seeds <- ("NatEdDens")}
  if(i==9){seeds <- ("RoadDens")}
  if(i==10){seeds <- ("BldgDens")}
  if(i==11){seeds <- ("DistNatEdge")}
  if(i==12){seeds <- ("ANTH")}
  if(i==13){seeds <- ("GRASS")}
  if(i==14){seeds <- ("NAT")}
  
  # Create a data frame to store the results from all the models.
  Seedmods_lin[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(seeds)){
    
    model <- glm(Seed ~ ContF_Seed[,j], # Create the predictive model
                 data=ContF_Seed, family = binomial)
    
    Seedmods_lin[[i]][j,1] <- j # First column is the predictor
    Seedmods_lin[[i]][j,2] <- AIC(model) # Second column is AIC
    Seedmods_lin[[i]][j,3:4] <- model$coefficients # thirs and fourth column are beta for intercept and predictor
    Seedmods_lin[[i]][j,5:6] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Seedmods_lin[[i]][j,7:10] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Seedmods_lin[[i]]) <- c("predictor", "AIC", "Intercept_Beta", "Predictor_beta", "Intercept_pval", "Predictor_pval", "Intercept_LowerCI", "Predictor_LowerCI", "Intercept_UpperCI", "Predictor_UpperCI")
}
SeedUnivariates_lin <- dplyr::bind_rows(Seedmods_lin)

SeedUnivariates_lin <- SeedUnivariates_lin %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictor_pval, Predictor_LowerCI, Predictor_UpperCI) %>%
  dplyr::rename("Beta" = Predictor_beta,
                "pvalue" = Predictor_pval,
                "LowCI" = Predictor_LowerCI,
                "UppCI" = Predictor_UpperCI)

#Now repeat the process for variables modelled as quadratics
Seedmods_quad <- list()
for(i in c(1:14)) {
  if(i==1){seeds <- ("DistBldg")}
  if(i==2){seeds <- ("DistRoad")}
  if(i==3){seeds <- ("DistNatEdge")}
  if(i==4){seeds <- ("DistScat")}
  if(i==5){seeds <- ("DistTrail")}
  if(i==6){seeds <- ("DistCamp")}
  if(i==7){seeds <- ("DistPG")}
  if(i==8){seeds <- ("NatEdDens")}
  if(i==9){seeds <- ("RoadDens")}
  if(i==10){seeds <- ("BldgDens")}
  if(i==11){seeds <- ("DistNatEdge")}
  if(i==12){seeds <- ("ANTH")}
  if(i==13){seeds <- ("GRASS")}
  if(i==14){seeds <- ("NAT")}
  
  # Create a data frame to store the results from all the models.
  Seedmods_quad[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(seeds)){
    
    model <- glm(Seed ~ ContF_Seed[,j] + I(ContF_Seed[,j]^2), # Create the predictive model
                 data=ContF_Seed, family = binomial)
    
    Seedmods_quad[[i]][j,1] <- j # First column is the predictor
    Seedmods_quad[[i]][j,2] <- AIC(model) # Second column is AIC
    Seedmods_quad[[i]][j,3:5] <- model$coefficients # thirs, fourth,fifth column are beta for intercept and predictor, predictor^2
    Seedmods_quad[[i]][j,6:8] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Seedmods_quad[[i]][j,9:14] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Seedmods_quad[[i]]) <- c("predictor", "AIC", "Intercept_Beta", 
                                    "Predictor_beta", "Predictore_beta2", 
                                    "Intercept_pval", "Predictor_pval", 
                                    "Predictor_pval2", "Intercept_LowerCI", 
                                    "Predictor_LowerCI", "Predictor_LowerCI2",
                                    "Intercept_UpperCI", "Predictor_UpperCI",
                                    "Predictor_UpperCI2")
}
SeedUnivariates_quad <- dplyr::bind_rows(Seedmods_quad)

SeedUnivariates_quad <- SeedUnivariates_quad %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictore_beta2,
                Predictor_pval, Predictor_pval2, Predictor_LowerCI, Predictor_UpperCI,
                Predictor_LowerCI2, Predictor_UpperCI2) %>%
  dplyr::rename("Q_Beta(lin)" = Predictor_beta,
                "Q_Beta(quad)" = Predictore_beta2,
                "Q_pvalue(lin)" = Predictor_pval,
                "Q_pvalue(quad)" = Predictor_pval2,
                "Q_LowCI(lin)" = Predictor_LowerCI,
                "Q_UppCI(lin)" = Predictor_UpperCI,
                "Q_LowCI(quad)" = Predictor_LowerCI2,
                "Q_UppCI(quad)" = Predictor_UpperCI2,
                "AIC_quad" = AIC)

SeedUnivariates <- cbind(SeedUnivariates_lin, SeedUnivariates_quad)
write.csv(SeedUnivariates, file = "SeedUnivariatesMarch2023.csv")



#NATURAL PREY
ContF_NatPrey <- ContF %>%
  dplyr::select(NatPrey, DistBldg, DistRoad, DistNatEdge, DistScat, 
                NatEdDens, RoadDens, BldgDens, NAT, ANTH, GRASS,
                DistTrail, DistCamp, DistPG)
  
#Get univariates results for Nat prey
#start by getting linear model info in a table
NatPreymods_lin <- list()
  for(i in c(1:14)) {
    if(i==1){natpreys <- ("DistBldg")}
    if(i==2){natpreys <- ("DistRoad")}
    if(i==3){natpreys <- ("DistNatEdge")}
    if(i==4){natpreys <- ("DistScat")}
    if(i==5){natpreys <- ("DistTrail")}
    if(i==6){natpreys <- ("DistCamp")}
    if(i==7){natpreys <- ("DistPG")}
    if(i==8){natpreys <- ("NatEdDens")}
    if(i==9){natpreys <- ("RoadDens")}
    if(i==10){natpreys <- ("BldgDens")}
    if(i==11){natpreys <- ("DistNatEdge")}
    if(i==12){natpreys <- ("ANTH")}
    if(i==13){natpreys <- ("GRASS")}
    if(i==14){natpreys <- ("NAT")}
    
  # Create a data frame to store the results from all the models.
  NatPreymods_lin[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(natpreys)){
    
    model <- glm(NatPrey ~ ContF_NatPrey[,j], # Create the predictive model
                 data=ContF_NatPrey, family = binomial)
    
    NatPreymods_lin[[i]][j,1] <- j # First column is the predictor
    NatPreymods_lin[[i]][j,2] <- AIC(model) # Second column is AIC
    NatPreymods_lin[[i]][j,3:4] <- model$coefficients # thirs and fourth column are beta for intercept and predictor
    NatPreymods_lin[[i]][j,5:6] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    NatPreymods_lin[[i]][j,7:10] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(NatPreymods_lin[[i]]) <- c("predictor", "AIC", "Intercept_Beta", "Predictor_beta", "Intercept_pval", "Predictor_pval", "Intercept_LowerCI", "Predictor_LowerCI", "Intercept_UpperCI", "Predictor_UpperCI")
}
NatPreyUnivariates_lin <- dplyr::bind_rows(NatPreymods_lin)

NatPreyUnivariates_lin <- NatPreyUnivariates_lin %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictor_pval, Predictor_LowerCI, Predictor_UpperCI) %>%
  dplyr::rename("Beta" = Predictor_beta,
                "pvalue" = Predictor_pval,
                "LowCI" = Predictor_LowerCI,
                "UppCI" = Predictor_UpperCI)

#Now repeat the process for variables modelled as quadratics
NatPreymods_quad <- list()
for(i in c(1:14)) {
  if(i==1){natpreys <- ("DistBldg")}
  if(i==2){natpreys <- ("DistRoad")}
  if(i==3){natpreys <- ("DistNatEdge")}
  if(i==4){natpreys <- ("DistScat")}
  if(i==5){natpreys <- ("DistTrail")}
  if(i==6){natpreys <- ("DistCamp")}
  if(i==7){natpreys <- ("DistPG")}
  if(i==8){natpreys <- ("NatEdDens")}
  if(i==9){natpreys <- ("RoadDens")}
  if(i==10){natpreys <- ("BldgDens")}
  if(i==11){natpreys <- ("DistNatEdge")}
  if(i==12){natpreys <- ("ANTH")}
  if(i==13){natpreys <- ("GRASS")}
  if(i==14){natpreys <- ("NAT")}
  
  
  # Create a data frame to store the results from all the models.
  NatPreymods_quad[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(natpreys)){
    
    model <- glm(NatPrey ~ ContF_NatPrey[,j] + I(ContF_NatPrey[,j]^2), # Create the predictive model
                 data=ContF_NatPrey, family = binomial)
    
    NatPreymods_quad[[i]][j,1] <- j # First column is the predictor
    NatPreymods_quad[[i]][j,2] <- AIC(model) # Second column is AIC
    NatPreymods_quad[[i]][j,3:5] <- model$coefficients # thirs, fourth,fifth column are beta for intercept and predictor, predictor^2
    NatPreymods_quad[[i]][j,6:8] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    NatPreymods_quad[[i]][j,9:14] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(NatPreymods_quad[[i]]) <- c("predictor", "AIC", "Intercept_Beta", 
                                       "Predictor_beta", "Predictore_beta2", 
                                       "Intercept_pval", "Predictor_pval", 
                                       "Predictor_pval2", "Intercept_LowerCI", 
                                       "Predictor_LowerCI", "Predictor_LowerCI2",
                                       "Intercept_UpperCI", "Predictor_UpperCI",
                                       "Predictor_UpperCI2")
}
NatPreyUnivariates_quad <- dplyr::bind_rows(NatPreymods_quad)

NatPreyUnivariates_quad <- NatPreyUnivariates_quad %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictore_beta2,
                Predictor_pval, Predictor_pval2, Predictor_LowerCI, Predictor_UpperCI,
                Predictor_LowerCI2, Predictor_UpperCI2) %>%
  dplyr::rename("Q_Beta(lin)" = Predictor_beta,
                "Q_Beta(quad)" = Predictore_beta2,
                "Q_pvalue(lin)" = Predictor_pval,
                "Q_pvalue(quad)" = Predictor_pval2,
                "Q_LowCI(lin)" = Predictor_LowerCI,
                "Q_UppCI(lin)" = Predictor_UpperCI,
                "Q_LowCI(quad)" = Predictor_LowerCI2,
                "Q_UppCI(quad)" = Predictor_UpperCI2,
                "AIC_quad" = AIC)

NatPreyUnivariates <- cbind(NatPreyUnivariates_lin, NatPreyUnivariates_quad)
write.csv(NatPreyUnivariates, file = "NatPreyUnivariatesMarch23.csv")


#VEGETATION
ContF_Veg <- ContF %>%
  dplyr::select(Veg, DistBldg, DistRoad, DistNatEdge, DistScat, 
                NatEdDens, RoadDens, BldgDens, NAT, ANTH, GRASS,
                DistTrail, DistCamp, DistPG)
  
#Get univariates results for VEG
#start by getting linear model info in a table
Vegmods_lin <- list()
for(i in c(1:14)) {
  if(i==1){vegs <- ("DistBldg")}
  if(i==2){vegs <- ("DistRoad")}
  if(i==3){vegs <- ("DistNatEdge")}
  if(i==4){vegs <- ("DistScat")}
  if(i==5){vegs <- ("DistTrail")}
  if(i==6){vegs <- ("DistCamp")}
  if(i==7){vegs <- ("DistPG")}
  if(i==8){vegs <- ("NatEdDens")}
  if(i==9){vegs <- ("RoadDens")}
  if(i==10){vegs <- ("BldgDens")}
  if(i==11){vegs <- ("DistNatEdge")}
  if(i==12){vegs <- ("ANTH")}
  if(i==13){vegs <- ("GRASS")}
  if(i==14){vegs <- ("NAT")}
  
  
  # Create a data frame to store the results from all the models.
  Vegmods_lin[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(vegs)){
    
    model <- glm(Veg ~ ContF_Veg[,j], # Create the predictive model
                 data=ContF_Veg, family = binomial)
    
    Vegmods_lin[[i]][j,1] <- j # First column is the predictor
    Vegmods_lin[[i]][j,2] <- AIC(model) # Second column is AIC
    Vegmods_lin[[i]][j,3:4] <- model$coefficients # thirs and fourth column are beta for intercept and predictor
    Vegmods_lin[[i]][j,5:6] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Vegmods_lin[[i]][j,7:10] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Vegmods_lin[[i]]) <- c("predictor", "AIC", "Intercept_Beta", "Predictor_beta", "Intercept_pval", "Predictor_pval", "Intercept_LowerCI", "Predictor_LowerCI", "Intercept_UpperCI", "Predictor_UpperCI")
}
VegUnivariates_lin <- dplyr::bind_rows(Vegmods_lin)

VegUnivariates_lin <- VegUnivariates_lin %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictor_pval, Predictor_LowerCI, Predictor_UpperCI) %>%
  dplyr::rename("Beta" = Predictor_beta,
                "pvalue" = Predictor_pval,
                "LowCI" = Predictor_LowerCI,
                "UppCI" = Predictor_UpperCI)

#Now repeat the process for variables modelled as quadratics
Vegmods_quad <- list()
for(i in c(1:14)) {
  if(i==1){vegs <- ("DistBldg")}
  if(i==2){vegs <- ("DistRoad")}
  if(i==3){vegs <- ("DistNatEdge")}
  if(i==4){vegs <- ("DistScat")}
  if(i==5){vegs <- ("DistTrail")}
  if(i==6){vegs <- ("DistCamp")}
  if(i==7){vegs <- ("DistPG")}
  if(i==8){vegs <- ("NatEdDens")}
  if(i==9){vegs <- ("RoadDens")}
  if(i==10){vegs <- ("BldgDens")}
  if(i==11){vegs <- ("DistNatEdge")}
  if(i==12){vegs <- ("ANTH")}
  if(i==13){vegs <- ("GRASS")}
  if(i==14){vegs <- ("NAT")}
  
  # Create a data frame to store the results from all the models.
  Vegmods_quad[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(vegs)){
    
    model <- glm(Veg ~ ContF_Veg[,j] + I(ContF_Veg[,j]^2), # Create the predictive model
                 data=ContF_Veg, family = binomial)
    
    Vegmods_quad[[i]][j,1] <- j # First column is the predictor
    Vegmods_quad[[i]][j,2] <- AIC(model) # Second column is AIC
    Vegmods_quad[[i]][j,3:5] <- model$coefficients # thirs, fourth,fifth column are beta for intercept and predictor, predictor^2
    Vegmods_quad[[i]][j,6:8] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Vegmods_quad[[i]][j,9:14] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Vegmods_quad[[i]]) <- c("predictor", "AIC", "Intercept_Beta", 
                                   "Predictor_beta", "Predictore_beta2", 
                                   "Intercept_pval", "Predictor_pval", 
                                   "Predictor_pval2", "Intercept_LowerCI", 
                                   "Predictor_LowerCI", "Predictor_LowerCI2",
                                   "Intercept_UpperCI", "Predictor_UpperCI",
                                   "Predictor_UpperCI2")
}
VegUnivariates_quad <- dplyr::bind_rows(Vegmods_quad)

VegUnivariates_quad <- VegUnivariates_quad %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictore_beta2,
                Predictor_pval, Predictor_pval2, Predictor_LowerCI, Predictor_UpperCI,
                Predictor_LowerCI2, Predictor_UpperCI2) %>%
  dplyr::rename("Q_Beta(lin)" = Predictor_beta,
                "Q_Beta(quad)" = Predictore_beta2,
                "Q_pvalue(lin)" = Predictor_pval,
                "Q_pvalue(quad)" = Predictor_pval2,
                "Q_LowCI(lin)" = Predictor_LowerCI,
                "Q_UppCI(lin)" = Predictor_UpperCI,
                "Q_LowCI(quad)" = Predictor_LowerCI2,
                "Q_UppCI(quad)" = Predictor_UpperCI2,
                "AIC_quad" = AIC)

VegUnivariates <- cbind(VegUnivariates_lin, VegUnivariates_quad)
write.csv(VegUnivariates, file = "VegUnivariatesMarch2023.csv")

#GARBAGE
ContF_Garbage <- ContF %>%
  dplyr::select(Garbage, DistBldg, DistRoad, DistNatEdge, DistScat, 
                NatEdDens, RoadDens, BldgDens, NAT, ANTH, GRASS,
                DistTrail, DistCamp, DistPG)

#Get univariates results for GARBAGE
#start by getting linear model info in a table
Garbagemods_lin <- list()
for(i in c(1:14)) {
  if(i==1){garb <- ("DistBldg")}
  if(i==2){garb <- ("DistRoad")}
  if(i==3){garb <- ("DistNatEdge")}
  if(i==4){garb <- ("DistScat")}
  if(i==5){garb <- ("DistTrail")}
  if(i==6){garb <- ("DistCamp")}
  if(i==7){garb <- ("DistPG")}
  if(i==8){garb <- ("NatEdDens")}
  if(i==9){garb <- ("RoadDens")}
  if(i==10){garb <- ("BldgDens")}
  if(i==11){garb <- ("DistNatEdge")}
  if(i==12){garb <- ("ANTH")}
  if(i==13){garb <- ("GRASS")}
  if(i==14){garb <- ("NAT")}
  
  # Create a data frame to store the results from all the models.
  Garbagemods_lin[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(garb)){
    
    model <- glm(Garbage ~ ContF_Garbage[,j], # Create the predictive model
                 data=ContF_Garbage, family = binomial)
    
    Garbagemods_lin[[i]][j,1] <- j # First column is the predictor
    Garbagemods_lin[[i]][j,2] <- AIC(model) # Second column is AIC
    Garbagemods_lin[[i]][j,3:4] <- model$coefficients # thirs and fourth column are beta for intercept and predictor
    Garbagemods_lin[[i]][j,5:6] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Garbagemods_lin[[i]][j,7:10] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Garbagemods_lin[[i]]) <- c("predictor", "AIC", "Intercept_Beta", "Predictor_beta", "Intercept_pval", "Predictor_pval", "Intercept_LowerCI", "Predictor_LowerCI", "Intercept_UpperCI", "Predictor_UpperCI")
}
GarbageUnivariates_lin <- dplyr::bind_rows(Garbagemods_lin)

GarbageUnivariates_lin <- GarbageUnivariates_lin %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictor_pval, Predictor_LowerCI, Predictor_UpperCI) %>%
  dplyr::rename("Beta" = Predictor_beta,
                "pvalue" = Predictor_pval,
                "LowCI" = Predictor_LowerCI,
                "UppCI" = Predictor_UpperCI)

#Now repeat the process for variables modelled as quadratics
Garbagemods_quad <- list()
for(i in c(1:14)) {
  if(i==1){garb <- ("DistBldg")}
  if(i==2){garb <- ("DistRoad")}
  if(i==3){garb <- ("DistNatEdge")}
  if(i==4){garb <- ("DistScat")}
  if(i==5){garb <- ("DistTrail")}
  if(i==6){garb <- ("DistCamp")}
  if(i==7){garb <- ("DistPG")}
  if(i==8){garb <- ("NatEdDens")}
  if(i==9){garb <- ("RoadDens")}
  if(i==10){garb <- ("BldgDens")}
  if(i==11){garb <- ("DistNatEdge")}
  if(i==12){garb <- ("ANTH")}
  if(i==13){garb <- ("GRASS")}
  if(i==14){garb <- ("NAT")}
  
  # Create a data frame to store the results from all the models.
  Garbagemods_quad[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(garb)){
    
    model <- glm(Garbage ~ ContF_Garbage[,j] + I(ContF_Garbage[,j]^2), # Create the predictive model
                 data=ContF_Garbage, family = binomial)
    
    Garbagemods_quad[[i]][j,1] <- j # First column is the predictor
    Garbagemods_quad[[i]][j,2] <- AIC(model) # Second column is AIC
    Garbagemods_quad[[i]][j,3:5] <- model$coefficients # thirs, fourth,fifth column are beta for intercept and predictor, predictor^2
    Garbagemods_quad[[i]][j,6:8] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Garbagemods_quad[[i]][j,9:14] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Garbagemods_quad[[i]]) <- c("predictor", "AIC", "Intercept_Beta", 
                                       "Predictor_beta", "Predictore_beta2", 
                                       "Intercept_pval", "Predictor_pval", 
                                       "Predictor_pval2", "Intercept_LowerCI", 
                                       "Predictor_LowerCI", "Predictor_LowerCI2",
                                       "Intercept_UpperCI", "Predictor_UpperCI",
                                       "Predictor_UpperCI2")
}
GarbageUnivariates_quad <- dplyr::bind_rows(Garbagemods_quad)

GarbageUnivariates_quad <- GarbageUnivariates_quad %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictore_beta2,
                Predictor_pval, Predictor_pval2, Predictor_LowerCI, Predictor_UpperCI,
                Predictor_LowerCI2, Predictor_UpperCI2) %>%
  dplyr::rename("Q_Beta(lin)" = Predictor_beta,
                "Q_Beta(quad)" = Predictore_beta2,
                "Q_pvalue(lin)" = Predictor_pval,
                "Q_pvalue(quad)" = Predictor_pval2,
                "Q_LowCI(lin)" = Predictor_LowerCI,
                "Q_UppCI(lin)" = Predictor_UpperCI,
                "Q_LowCI(quad)" = Predictor_LowerCI2,
                "Q_UppCI(quad)" = Predictor_UpperCI2,
                "AIC_quad" = AIC)

GarbageUnivariates <- cbind(GarbageUnivariates_lin, GarbageUnivariates_quad)
write.csv(GarbageUnivariates, file = "GarbageUnivariatesMarch2023.csv")



#Proceed with all subsets approach for FRUIT ; I used P = 0.25 as a cutoff for inclusion
#First, create a dataframe with only the variables that will be going in to the model
ContF_Fruit2 <- ContF_Fruit %>%
  dplyr::select(Fruit, DistTrail, NAT, ANTH, DistScat, BldgDens, GRASS, RoadDens)

cor(ContF_Fruit2[sapply(ContF_Fruit2, is.numeric)], method = c("pearson")) #I will use 0.6 as cutoff 
#ANTH and NAT are correlated --> RETAIN NAT
FruitDredge <- glm(Fruit ~ DistTrail + NAT + DistScat + BldgDens + GRASS + RoadDens + 
                     I(NAT^2) + I(DistScat^2) + I(RoadDens^2),
                   family = binomial, data = ContF_Fruit2)

subsetFruit <- expression(dc(`RoadDens`, `I(RoadDens^2)`) &
                            dc(`NAT`, `I(NAT^2)`) &
                            dc(`DistScat`, `I(DistScat^2)`))

options(na.action = "na.fail")
FruitRes <- dredge(FruitDredge, subset = subsetFruit, trace = 2, rank = BIC)

Fruitmod1a <- get.models(FruitRes, 1)[[1]]
Fruitmod2a <- get.models(FruitRes, 2)[[1]]

Fruitmod1 <- glm(formula = Fruit ~ DistTrail + 1, family = binomial, data = ContF_Fruit2)

Fruitmod2 <- glm(formula = Fruit ~ DistScat + DistTrail + 1, family = binomial, 
                 data = ContF_Fruit2)

plot_summs(Fruitmod1, Fruitmod2)

#Export all the model info
export_summs(Fruitmod1, Fruitmod2,
             ci_level = 0.95,
             stars = NULL,
             error_format = '(p = {p.value}), CI [{conf.low}, {conf.high}]',
             model.names = c("1", "2"),
             digits = 3,
             to.file = "xlsx",
             file.name = "March2023FruitBICModResults.xlsx")


#Proceed with all subsets approach for ANTHRO ; I used P = 0.25 as a cutoff for inclusion
#First, create a dataframe with only the variables that will be going in to the model
ContF_Anthro2 <- ContF_Anthro %>%
  dplyr::select(Anthro, DistCamp, ANTH, NAT, RoadDens, DistScat, DistTrail, DistPG, BldgDens, DistNatEdge)

cor(ContF_Anthro2[sapply(ContF_Anthro2, is.numeric)], method = c("pearson")) #I will use 0.6 as cutoff 
#ANTH and NAT are correlated --> RETAIN ANTH
AnthroDredge <- glm(Anthro ~ DistCamp + ANTH + RoadDens + DistScat + DistTrail + DistPG + BldgDens + I(DistCamp^2) + I(DistPG^2) + I(DistCamp^2), 
                    family = binomial, data = ContF_Anthro2)

subsetAnthro <- expression(dc(`DistCamp`, `I(DistCamp^2)`) &
                             dc(`DistNatEdge`, `I(DistNatEdge^2)`) &
                             dc(`DistPG`, `I(DistPG^2)`) &
                             dc(`DistCamp`, `I(DistCamp^2)`))

options(na.action = "na.fail")
AnthroRes <- dredge(AnthroDredge, subset = subsetAnthro, trace = 2, rank = BIC)

Anthromod1a <- get.models(AnthroRes, 1)[[1]]
Anthromod2a <- get.models(AnthroRes, 2)[[1]]

Anthromod1 <- glm(formula = Anthro ~ DistCamp + I(DistCamp^2) + RoadDens + 
                    1, family = binomial, data = ContF_Anthro2)
Anthromod2 <- glm(formula = Anthro ~ DistCamp + I(DistCamp^2) + DistScat + 
                    RoadDens + 1, family = binomial, data = ContF_Anthro2)

plot_summs(Anthromod1, Anthromod2)

AnthroBIC1 <- plot_summs(Anthromod1,
                         colors = "#00695C",
                         point.size = 30,
                         point.shape = FALSE,
                         coefs = c("Decay distance to camp" = "Dec_Camp",
                                   "Decay distance to camp (Squared)" = "I(Dec_Camp^2)",
                                   "Decay distance to scat" = "Dec_Scat",
                                   "Road density" = "RoadDens"))
AnthroBIC <- AnthroBIC1 +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Coefficient estimate") +
  scale_shape_manual(values = c(19)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1))



#Export all the model info
export_summs(Anthromod1, Anthromod2,
             ci_level = 0.95,
             stars = NULL,
             error_format = '(p = {p.value}), CI [{conf.low}, {conf.high}]',
             model.names = c("1", "2"),
             digits = 3,
             to.file = "xlsx",
             file.name = "March2023AnthroBICModResults.xlsx")


#Seed
#First, create a dataframe with only the variables that will be going in to the model
ContF_Seed2 <- ContF_Seed %>%
  dplyr::select(Seed, DistBldg, DistRoad, DistPG, NatEdDens, ANTH, RoadDens, DistCamp)

cor(ContF_Seed2[sapply(ContF_Seed2, is.numeric)], method = c("pearson")) #I will use 0.6 as cutoff 
#Dist building anfdist road are correlated --> retain DistBLDG

SeedDredge <- glm(Seed ~ DistBldg + DistPG + NatEdDens + ANTH + RoadDens + DistCamp + I(DistBldg^2) + I(NatEdDens^2) + 
                    I(RoadDens^2) + I(ANTH^2) + I(DistCamp^2),
                  family = binomial, data = ContF_Seed2)

subsetSeed <- expression(dc(`DistBldg`, `I(DistBldg^2)`) &
                           dc(`NatEdDens`, `I(NatEdDens^2)`) &
                           dc(`RoadDens`, `I(RoadDens^2)`) &
                           dc(`ANTH`, `I(ANTH^2)`) &
                           dc(`DistCamp`, `I(DistCamp^2)`))

options(na.action = "na.fail")
SeedRes <- dredge(SeedDredge, subset = subsetSeed, trace = 2, rank = BIC)

Seedmod1a <- get.models(SeedRes, 1)[[1]]
Seedmod2a <- get.models(SeedRes, 2)[[1]]
Seedmod3a <- get.models(SeedRes, 3)[[1]]

Seedmod1 <- glm(formula = Seed ~ ANTH + DistBldg + I(DistBldg^2) + DistPG + 
                  1, family = binomial, data = ContF_Seed2)
Seedmod2 <- glm(formula = Seed ~ ANTH + DistBldg + I(DistBldg^2) + DistCamp + 
                  I(DistCamp^2) + DistPG + 1, family = binomial, data = ContF_Seed2)
Seedmod3 <- glm(formula = Seed ~ DistBldg + I(DistBldg^2) + DistPG + 1, family = binomial, data = ContF_Seed2)
  
plot_summs(Seedmod1, Seedmod2, Seedmod3)

SeedBIC1 <- plot_summs(Seedmod1, Seedmod2, Seedmod3, Seedmod4,
                       colors = c("#00796B", "#00796B", "#00796B", "#00796B"),
                       point.size = 30,
                       point.shape = FALSE, 
                       coefs = c("Decay distance to water" = "Dec_Water",
                                 "Decay distance to camp" = "Dec_Camp",
                                 "Distance to playground" = "DistPG",
                                 "Decay distance to junction" = "Dec_Junc",
                                 "Distance to maintained trail" = "DistTrail",
                                 "East index" = "East"))
SeedBIC <- SeedBIC1 +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Coefficient estimate") +
  scale_shape_manual(values = c(17,17,17,17)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1))



#Export all the model info
export_summs(Seedmod1, Seedmod2, Seedmod3,
             ci_level = 0.95,
             stars = NULL,
             error_format = '(p = {p.value}), CI [{conf.low}, {conf.high}]',
             model.names = c("1", "2", "3"),
             digits = 3,
             to.file = "xlsx",
             file.name = "March2023SeedBICModResults.xlsx")

#NAT PREY
#First, create a dataframe with only the variables that will be going in to the model
ContF_NatPrey2 <- ContF_NatPrey %>%
  dplyr::select(NatPrey, DistBldg, ANTH, NAT, DistScat, DistPG, RoadDens)

cor(ContF_NatPrey2[sapply(ContF_NatPrey2, is.numeric)], method = c("pearson")) #I will use 0.6 as cutoff 
#ANTH AND NAT correlated --> keep ANTH

NatPreyDredge <- glm(NatPrey ~ DistBldg + ANTH + DistScat + DistPG + RoadDens + I(DistScat^2),
                     family = binomial, data = ContF_NatPrey2)

subsetNatPrey <- expression(dc(`DistScat`, `I(DistScat^2)`))
                           
options(na.action = "na.fail")
NatPreyRes <- dredge(NatPreyDredge, trace = 2, rank = BIC, subset = subsetNatPrey)

NatPreymod1a <- get.models(NatPreyRes, 1)[[1]]

NatPreymod1 <- glm(formula = NatPrey ~ DistBldg + 1, family = binomial, data = ContF_NatPrey2)


plot_summs(NatPreymod1)

NatPreyBIC1 <- plot_summs(NatPreymod1, NatPreymod2,
                          colors = c("#00897B", "#00897B"),
                          point.size = 30,
                          point.shape = FALSE, 
                          coefs = c("Decay distance to camp" = "Dec_Camp",
                                    "Decay distance to scat" = "Dec_Scat",
                                    "East index" = "East",
                                    "Decay distance to building" = "Dec_Bldg"))
NatPreyBIC <- NatPreyBIC1 +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Coefficient estimate") +
  scale_shape_manual(values = c(18,18)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1))



#Export all the model info
export_summs(NatPreymod1,
             ci_level = 0.95,
             stars = NULL,
             error_format = '(p = {p.value}), CI [{conf.low}, {conf.high}]',
             model.names = c("1"),
             digits = 3,
             to.file = "xlsx",
             file.name = "March2023NatPreyBICModResults.xlsx")


#Proceed with all subsets approach for VEGETATION ; I used P = 0.25 as a cutoff for inclusion
#First, create a dataframe with only the variables that will be going in to the model
ContF_Veg2 <- ContF_Veg %>%
  dplyr::select(Veg, DistPG, DistRoad, DistBldg, DistScat, ANTH, NatEdDens, DistCamp, DistTrail)

cor(ContF_Veg2[sapply(ContF_Veg2, is.numeric)], method = c("pearson")) #I will use 0.6 as cutoff 
#Dist bldg and dtsr road are correlated --> retain dist bldg

VegDredge <- glm(Veg ~ DistPG + DistBldg + DistScat + ANTH + NatEdDens + DistCamp+ DistTrail + I(DistBldg^2) + I(DistScat^2) + I(NatEdDens^2) +
                   + I(DistTrail^2) + I(DistCamp^2), family = binomial, data = ContF_Veg2)

subsetVeg <- expression(dc(`NatEdDens`, `I(NatEdDens^2)`) &
                          dc(`DistBldg`, `I(DistBldg^2)`) &
                          dc(`DistScat`, `I(DistScat^2)`) &
                          dc(`DistTrail`, `I(DistTrail^2)`) &
                          dc(`DistCamp`, `I(DistCamp^2)`))

options(na.action = "na.fail")
VegRes <- dredge(VegDredge, subset = subsetVeg, trace = 2, rank = BIC)

Vegmod1a <- get.models(VegRes, 1)[[1]]
Vegmod2a <- get.models(VegRes, 2)[[1]]
Vegmod3a <- get.models(VegRes, 3)[[1]]

Vegmod1 <- glm(formula = Veg ~ DistPG + 1, family = binomial, data = ContF_Veg2)
Vegmod2 <- glm(formula = Veg ~ DistBldg + I(DistBldg^2) + DistPG + 1, family = binomial, 
               data = ContF_Veg2)
Vegmod3 <- glm(formula = Veg ~ DistBldg + DistPG + 1, family = binomial, 
               data = ContF_Veg2)

plot_summs(Vegmod1, Vegmod2, Vegmod3)

VegBIC1 <- plot_summs(Vegmod1, Vegmod2,
                      colors = c("#009688", "#009688"),
                      point.size = 30,
                      point.shape = FALSE,
                      coefs = c("Decay distance to scat" = "Dec_Scat",
                                "Decay distance to water" = "Dec_Water",
                                "Distance to playground" = "DistPG"))
VegBIC <- VegBIC1 +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Coefficient estimate") +
  scale_shape_manual(values = c(25,25)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1))



#Export all the model info
export_summs(Vegmod1, Vegmod2, Vegmod3,
             ci_level = 0.95,
             stars = NULL,
             error_format = '(p = {p.value}), CI [{conf.low}, {conf.high}]',
             model.names = c("1", "2", "3"),
             digits = 3,
             to.file = "xlsx",
             file.name = "March2023VegBICModResults.xlsx")


#Proceed with all subsets approach for GARBAGE ; I used P = 0.25 as a cutoff for inclusion
#First, create a dataframe with only the variables that will be going in to the model
ContF_Garbage2 <- ContF_Garbage %>%
  dplyr::select(Garbage, RoadDens, DistScat, GRASS, NatEdDens, NAT)

cor(ContF_Garbage2[sapply(ContF_Garbage2, is.numeric)], method = c("pearson")) #I will use 0.6 as cutoff 

GarbageDredge <- glm(Garbage ~ RoadDens + DistScat + GRASS + NAT + NatEdDens,
                     family = binomial, data = ContF_Garbage2)

options(na.action = "na.fail")
GarbageRes <- dredge(GarbageDredge, trace = 2, rank = BIC)

Garbagemod1a <- get.models(GarbageRes, 2)[[1]]

Garbagemod1 <- glm(formula = Garbage ~ RoadDens + 1, family = binomial, data = ContF_Garbage2)

plot_summs(Garbagemod1)

GarbageBIC1 <- plot_summs(Garbagemod1, Garbagemod2, Garbagemod3,
                          colors = c("#26A69A", "#26A69A", "#26A69A"),
                          point.size = 30,
                          point.shape = FALSE,
                          coefs = c("Road density" = "RoadDens",
                                    "Decay distance to scat" = "Dec_Scat",
                                    "Decay distance to human trail" = "Dec_Trail"))

GarbageBIC <- GarbageBIC1 +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Coefficient estimate") +
  scale_shape_manual(values = c(21,21,21)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1))



#Export all the model info
export_summs(Garbagemod1,
             ci_level = 0.95,
             stars = NULL,
             error_format = '(p = {p.value}), CI [{conf.low}, {conf.high}]',
             model.names = c("1"),
             digits = 3,
             to.file = "xlsx",
             file.name = "GarbageBICModResults.xlsx")


#Get mean pseudo R2 for all models
A <- PseudoR2(Fruitmod1, which = "McFadden")
B <- PseudoR2(Fruitmod2, which = "McFadden")
mean(A,B)

C <-PseudoR2(Fruitmod1, which = "Nagelkerke")
D <- PseudoR2(Fruitmod2, which = "Nagelkerke")
mean(C,D)

PseudoR2(Anthromod1, which = "McFadden")
PseudoR2(Anthromod1, which = "Nagelkerke")

E <- PseudoR2(Seedmod1, which = "McFadden")
F <- PseudoR2(Seedmod2, which = "McFadden")
G <- PseudoR2(Seedmod3, which = "McFadden")
H <- PseudoR2(Seedmod4, which = "McFadden")
mean(E, F, G, H)

I <- PseudoR2(Seedmod1, which = "Nagelkerke")
J <- PseudoR2(Seedmod2, which = "Nagelkerke")
K <- PseudoR2(Seedmod3, which = "Nagelkerke")
L <- PseudoR2(Seedmod4, which = "Nagelkerke")
mean(I, J, K, L)

M <- PseudoR2(NatPreymod1, which = "McFadden")
N <- PseudoR2(NatPreymod2, which = "McFadden")
mean(M, N)

O <- PseudoR2(NatPreymod1, which = "Nagelkerke")
P <- PseudoR2(NatPreymod2, which = "Nagelkerke")
mean(O, P)

Q <- PseudoR2(Vegmod1, which = "McFadden")
R <- PseudoR2(Vegmod2, which = "McFadden")
mean(Q,R)

S <- PseudoR2(Vegmod1, which = "Nagelkerke")
T <- PseudoR2(Vegmod2, which = "Nagelkerke")
mean(S,T)

U <- PseudoR2(Garbagemod1, which = "McFadden")
V <- PseudoR2(Garbagemod2, which = "McFadden")
W <- PseudoR2(Garbagemod3, which = "McFadden")
mean(U, V, W)

X <- PseudoR2(Garbagemod1, which = "Nagelkerke")
Y <- PseudoR2(Garbagemod2, which = "Nagelkerke")
Z <- PseudoR2(Garbagemod3, which = "Nagelkerke")
mean(X, Y, Z)

#Get VIF for all models
vif(Fruitmod1)
vif(Fruitmod2)

vif(Anthromod1)

vif(Seedmod1)
vif(Seedmod2)
vif(Seedmod3)
vif(Seedmod4)

vif(NatPreymod1)
vif(NatPreymod2)

vif(Vegmod1)
vif(Vegmod2)

vif(Garbagemod1)
vif(Garbagemod2)
vif(Garbagemod3)

#Get ROC AUC for all models
A1 <- auc(roc(ContF_Fruit2$Fruit~fitted(Fruitmod1)))
A2 <- auc(roc(ContF_Fruit2$Fruit~fitted(Fruitmod2)))
mean(A1, A2)

auc(roc(ContF_Anthro2$Anthro~fitted(Anthromod1)))

B1 <- auc(roc(ContF_Seed2$Seed~fitted(Seedmod1)))
B2 <- auc(roc(ContF_Seed2$Seed~fitted(Seedmod2)))
B3 <- auc(roc(ContF_Seed2$Seed~fitted(Seedmod3)))
B4 <- auc(roc(ContF_Seed2$Seed~fitted(Seedmod4)))
mean(B1, B2, B3, B4)

C1 <- auc(roc(ContF_NatPrey2$NatPrey~fitted(NatPreymod1)))
C2 <- auc(roc(ContF_NatPrey2$NatPrey~fitted(NatPreymod2)))
mean(C1, C2)

D1 <- auc(roc(ContF_Veg2$Veg~fitted(Vegmod1)))
D2 <- auc(roc(ContF_Veg2$Veg~fitted(Vegmod2)))
mean(D1, D2)

E1 <- auc(roc(ContF_Garbage2$Garbage~fitted(Garbagemod1)))
E2 <- auc(roc(ContF_Garbage2$Garbage~fitted(Garbagemod2)))
E3 <- auc(roc(ContF_Garbage2$Garbage~fitted(Garbagemod3)))
mean(E1, E2, E3)

#Do cross validation for all models
#Fruit
data_ctrl <- trainControl(method = "cv", number = 5)
Fruit_caret1 <- train(Fruit ~ Dec_Camp + Dec_Scat + DistTrail, data = ContF_Fruit2,
                      trControl = data_ctrl,
                      method = "glm", na.action = na.pass)
Fruit_caret2 <- train(Fruit ~ Dec_Camp + Dec_PG + Dec_Scat + DistTrail, data = ContF_Fruit2,
                      trControl = data_ctrl,
                      method = "glm", na.action = na.pass)

Fruit_caret1
Fruit_caret2

mean(0.564419, 0.5609711)
mean(0.1172385, 0.1118293)

#Anthro
data_ctrlA <- trainControl(method = "cv", number = 5)
Anthro_caret1 <- train(Anthro ~ Dec_Camp + I(Dec_Camp^2) + Dec_Scat + 
                         RoadDens, data = ContF_Anthro2,
                       trControl = data_ctrlA,
                       method = "glm", na.action = na.pass)
Anthro_caret1


#Seed
data_ctrlS <- trainControl(method = "cv", number = 5)
Seed_caret1 <- train(Seed ~ Dec_Camp + Dec_Water, data = ContF_Seed2,
                     trControl = data_ctrlS,
                     method = "glm", na.action = na.pass)
Seed_caret2 <- train(Seed ~ Dec_Camp + Dec_Water + DistPG, data = ContF_Seed2,
                     trControl = data_ctrlS,
                     method = "glm", na.action = na.pass)
Seed_caret3 <- train(Seed ~ Dec_Camp + Dec_Water + DistPG + East, data = ContF_Seed2,
                     trControl = data_ctrlS,
                     method = "glm", na.action = na.pass)
Seed_caret4 <- train(Seed ~ Dec_Junc + Dec_Water + DistTrail, data = ContF_Seed2,
                     trControl = data_ctrlS,
                     method = "glm", na.action = na.pass)

Seed_caret1
Seed_caret2
Seed_caret3
Seed_caret4

#Nat Prey
data_ctrlN <- trainControl(method = "cv", number = 5)
NatPrey_caret1 <- train(NatPrey ~ Dec_Camp + Dec_Scat + East, data = ContF_NatPrey2,
                        trControl = data_ctrlN,
                        method = "glm", na.action = na.pass)
NatPrey_caret2 <- train(NatPrey ~ Dec_Bldg + Dec_Camp + Dec_Scat + East, data = ContF_NatPrey2,
                        trControl = data_ctrlN,
                        method = "glm", na.action = na.pass)

NatPrey_caret1
NatPrey_caret2

#Veg
data_ctrlV <- trainControl(method = "cv", number = 5)
Veg_caret1 <- train(Veg ~ Dec_Scat + Dec_Water + DistPG, data = ContF_Veg2,
                    trControl = data_ctrlV,
                    method = "glm", na.action = na.pass)
Veg_caret2 <- train(Veg ~ Dec_Scat + DistPG, data = ContF_Veg2,
                    trControl = data_ctrlV,
                    method = "glm", na.action = na.pass)

Veg_caret1
Veg_caret2

mean(0.8456956, 0.8439898)
mean(0, -0.003347101)

#Garbage
data_ctrlG <- trainControl(method = "cv", number = 5)
Garbage_caret1 <- train(Garbage ~ RoadDens, data = ContF_Garbage2,
                        trControl = data_ctrlG,
                        method = "glm", na.action = na.pass)
Garbage_caret2 <- train(Garbage ~ Dec_Scat, data = ContF_Garbage2,
                        trControl = data_ctrlG,
                        method = "glm", na.action = na.pass)
Garbage_caret3 <- train(Garbage ~ Dec_Trail, data = ContF_Garbage2,
                        trControl = data_ctrlG,
                        method = "glm", na.action = na.pass)

Garbage_caret1
Garbage_caret2
Garbage_caret3

##Create a large, multi-panel figure
#Start by making a legend
AA <- plot_summs(Garbagemod1, Garbagemod1, Garbagemod1, Garbagemod1, Garbagemod1, Garbagemod1,
                 colors = c("#004D40", "#00695C","#00796B","#00897B","#009688", "#26A69A"),
                 point.shape = FALSE,
                 model.names = c("Fruit", "Anthropogenic", "Birdseed", "Natural prey",
                                 "Vegetation", "Garbage"))

BB <- AA + 
  theme(axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(colour = "black", size = 12, face = "plain"),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
        legend.title = element_text(colour = "black", size = 12, face = "plain")) +
  labs(x = "Coefficient estimate") +
  scale_shape_manual(name = "Scat content",
                     values = c(15, 19, 17, 18, 25, 21),
                     breaks = c("Fruit", "Anthropogenic", "Birdseed",
                                "Natural prey", "Vegetation", "Garbage")) +
  scale_color_manual(name = "Scat content",
                     values = c("#004D40", "#00695C","#00796B","#00897B","#009688", "#26A69A"),
                     breaks = c("Fruit", "Anthropogenic", "Birdseed",
                                "Natural prey", "Vegetation", "Garbage"))

#Get legend
Contlegend <- cowplot::get_legend(BB)
grid.newpage()
grid.draw(Contlegend)

ContCharts <- ggarrange(FruitBIC, AnthroBIC, SeedBIC, NatPreyBIC, VegBIC, GarbageBIC,
                        nrow = 3, ncol = 2,
                        align="v")

#Save file
#ggsave(ContCharts, height = 8, width = 10, dpi = 600, file = "YYY.png", bg = "white")


ContCharts2 <- ggarrange(ContentChart1, Contlegend, nrow = 2, ncol = 1,
                         heights = c(2,1))

#Save file
ggsave(ContCharts2, height = 8, width = 3, dpi = 600, file = "XXX.png", bg = "white")


#Load as photos
ContPart1 <- image_read("XXX.png")
ContPart2 <- image_read("YYY.png")

#Annotate FirstPart
ContPart1_a <- image_annotate(ContPart1, text = "(A)", size = 125, color = "black",
                              gravity = "northwest")
image_browse(ContPart1_a)

#Annotate SecondPart
ContPart2_a <- image_annotate(ContPart2, text = "(B) Fruit", size = 125, color = "black", 
                              gravity = "northwest")
image_browse(ContPart2_a)

ContPart2_b <- image_annotate(ContPart2_a, text = "(C) Anthropogenic", size = 125, color = "black", 
                              gravity = "northeast", location = '+2000+0')
image_browse(ContPart2_b)

ContPart2_c <- image_annotate(ContPart2_b, text = "(D) Birdseed", size = 125, color = "black", 
                              gravity = "northwest", location = '+0+1550')
image_browse(ContPart2_c)

ContPart2_d <- image_annotate(ContPart2_c, text = "(E) Natural prey", size = 125, color = "black", 
                              gravity = "northeast", location = '+2100+1550')
image_browse(ContPart2_d)

ContPart2_e <- image_annotate(ContPart2_d, text = "(F) Vegetation", size = 125, color = "black", 
                              gravity = "northwest", location = '+0+3200')
image_browse(ContPart2_e)

ContPart2_f <- image_annotate(ContPart2_e, text = "(G) Garbage", size = 125, color = "black", 
                              gravity = "northwest", location = '+3050+3200')
image_browse(ContPart2_f)

#create one figure with three panels of photos
imgX <- c(ContPart1_a, ContPart2_f)

XXFig1 <- image_append(imgX)

image_browse(XXFig1)

image_write(XXFig1, path = "LargeContentFigureMov21.png", format = "png")


#####PHOTO FIGURE################
#Make a composite figure of 8 pictures
#Load all photos
PicGarb <- image_read("Garbage.jpg")
PicAnthro <- image_read("PresAnthro.jpg")
PicNatPrey <- image_read("natprey.jpg")
PicEstrus <- image_read("Estrus.jpg")
Picchipbag <- image_read("chipbag.jpg")
PicVeg <- image_read("Veg.jpg")
PicFruit <- image_read("Fruit.jpg")
PicSeed <- image_read("Seed.jpg")

#get info for images
image_info(PicGarb) #H
PicGarb <- image_scale(PicGarb, "2000")

image_info(PicAnthro) #E
PicAnthro <- image_scale(PicAnthro, "2000")

image_info(PicNatPrey) #I
PicNatPrey <- image_scale(PicNatPrey, "2000")

image_info(PicEstrus) #A
PicEstrus <- image_scale(PicEstrus, "2000")

image_info(Picchipbag) #B
Picchipbag <- image_scale(Picchipbag, "2000")

image_info(PicVeg) #G
PicVeg <- image_scale(PicVeg, "2000")

image_info(PicFruit) #D
PicFruit <- image_scale(PicFruit, "2000")

image_info(PicSeed) #F
PicSeed <- image_scale(PicSeed, "2000")

image_info(PicScatMP) #C
PicScatMP <- image_scale(PicScatMP, "2000")

#Annotate Pic1
PicGarba <- image_annotate(PicGarb, text = "G. Garbage", size = 150, color = "black", boxcolor = "white",
                           gravity = "northwest", location = '+150+250')
image_browse(PicGarba)

#Annotate Pic2
PicAnthroa <- image_annotate(PicAnthro, text = "D. Anthropogenic", size = 150, color = "black", boxcolor = "white",
                             gravity = "northwest", location = '+150+250')
image_browse(PicAnthroa)

#Annotate Pic3
PicNatPreya <- image_annotate(PicNatPrey, text = "H. Natural prey", size = 150, color = "black", boxcolor = "white",
                              gravity = "northwest", location = '+150+250')
image_browse(PicNatPreya)


#Annotate Pic3
PicEstrusa <- image_annotate(PicEstrus, text = "A. Double urination", size = 150, color = "black", boxcolor = "white",
                             gravity = "northwest", location = '+150+250')
image_browse(PicEstrusa)

#Annotate Pic3
PicChipbagA <- image_annotate(Picchipbag, text = "B. Conspicuous scat", size = 150, color = "black", boxcolor = "white",
                              gravity = "northwest", location = '+150+250')
image_browse(PicChipbagA)

#Annotate Pic3
PicVegA <- image_annotate(PicVeg, text = "F. Vegetation", size = 150, color = "black", boxcolor = "white",
                          gravity = "northwest", location = '+150+250')
image_browse(PicVegA)


#Annotate Pic3
PicFruitA <- image_annotate(PicFruit, text = "C. Fruit", size = 150, color = "black", boxcolor = "white",
                            gravity = "northwest", location = '+150+250')
image_browse(PicFruitA)


#Annotate Pic3
PicSeeda <- image_annotate(PicSeed, text = "E. Birdseed", size = 150, color = "black", boxcolor = "white",
                           gravity = "northwest", location = '+150+250')
image_browse(PicSeeda)

PicEstrusa <- image_border(PicEstrusa, color = "white", "10 x 10")
PicChipbagA <- image_border(PicChipbagA, color = "white", "10x10")
PicFruitA <- image_border(PicFruitA, color = "white", "10x10")
PicAnthroa <- image_border(PicAnthroa, color = "white", "10x10")
PicSeeda <- image_border(PicSeeda, color = "white", "10x10") 
PicVegA <- image_border(PicVegA, color = "white", "10x10")
PicGarba <- image_border(PicGarba, color = "white", "10x10")
PicNatPreya <- image_border(PicNatPreya, color = "white", "10x10")

#create one figure with four rows of two photos
img1 <- c(PicEstrusa, PicChipbagA)
img2 <- c(PicFruitA, PicAnthroa)
img3 <- c(PicSeeda, PicVegA)
img4 <- c(PicGarba, PicNatPreya)
PhotoFig1 <- image_append(img1)
PhotoFig2 <- image_append(img2)
PhotoFig3 <- image_append(img3)
PhotoFig4 <- image_append(img4)

PhotoFigFin <- c(PhotoFig1, PhotoFig2, PhotoFig3, PhotoFig4)
PhotoFigFin2 <- image_append(PhotoFigFin, stack = TRUE)

image_browse(PhotoFigFin2)

#Save three-paneled photo figure
#image_write(PhotoFigFin2, path = "PhotoCollageNov21.png", format = "png")

#And make a slightly different one for ppt presentations
#create one figure with four rows of two photos
img2_a <- c(PicFruitA, PicAnthroa, PicSeeda)
img3_a <- c(PicVegA, PicGarba, PicNatPreya)
PhotoFig1_a <- image_append(img2_a)
PhotoFig2_a <- image_append(img3_a)

PhotoFigFin_ppt <- c(PhotoFig1_a, PhotoFig2_a)
PhotoFigFin2_ppt <- image_append(PhotoFigFin_ppt, stack = TRUE)

image_browse(PhotoFigFin2_ppt)
image_write(PhotoFigFin2_ppt, path = "PhotoCollageNov21_forppt.png", format = "png")
