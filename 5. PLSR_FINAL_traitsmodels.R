# # # # # # # # #  PLSR - Analysing dunes trait data set # # # # # # # # # # # # # # # 

# # # 1st PhD chapter - Rodríguez-Arias,L.

# 12/2023

# PLSR ~ Reduction dimension strategy 
# An alternative to PCR is the Partial Least Squares (PLS) regression, which identifies
# new principal components that not only summarizes the original predictors (chosen by cross-validation), but also that are related to the outcome. 
# Summarize the data into few non-redundant components --> regression model.

# http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/152-principal-component-and-partial-least-squares-regression-essentials/#partial-least-squares-regression

rm(list=ls())

# Library ----
library(tidyverse)
library(tidylog)
library(pls)
library(caret)
library(matrixStats)
library(ggsci)
library(readr)
library(dplyr)

# Directory----

work.dir=("C:/Users/lucia/Desktop/PhD_Blanes/3_R/Paper1/Paper1_mechanism_stabilization")

data.dir=paste(work.dir,"Data",sep="/")
plots.dir=paste(work.dir,"Plots",sep="/")
model.outs=paste(work.dir,"ModelOuts",sep="/")

# Data.frame: patches & vegetation traits ----
setwd(data.dir)
dir()
M <- read_csv2("manchasR.csv")
Traits <- read_csv2("Traits.csv")

# # Regaridng with the species traits----
Traits$spp #creating a new column
Traits <- Traits%>%
  mutate(specie = case_when(
    spp=="Ammophila"~ "Ammophila arenaria",
    spp == "Cakile" ~ "Cakile maritima",
    spp=="Echinophora" ~"Echinophora spinosa",
    spp =="Euphorbia" ~"Euphorbia paralias",
    spp =="Salicornia" ~ "Salicornia europaea",
    spp =="Sarcocornia" ~ "Sarcocornia fruticosum",
    spp == "Salsola" ~ "Salsola kali",
    spp == "Limonium" ~"Limonium vigoi",
    spp == "Sueda"~"Suaeda maritima",
    spp == "Sporobolus" ~"Sporobolus pungens",
    spp == "Salsola"~"Salsola kali",
    spp == "Phargmites"~"Phargmites australis",
    spp == "Sedum" ~ "Sedum spp.",
    spp == "Elytrigia"~"Elymus farctus",
    spp == "Artemisia" ~"Arteminia maritima",
    spp == "Atriplex" ~ "Atriplex spp.",
    TRUE ~ spp
  ))

# drop the spp & varibles that we won't need
Traits <- Traits %>%
  filter(!acro %in% c("ATRIP", "ELYT","PHAR","ART","SED"))

# Variables for the PCA (I keep the averages, and eliminate the variables I don't want: categorical or meaningless)
Traits<- Traits[,-c(1:8,20,21,23,24,26,27,35)]
#filter row: mean values + columns:variable sof interest

# We'll doing the PLSR with all species but for some of them, we couldn't get the traits
Traits_full <- Traits %>% 
  select_if(~!any(is.na(.))) # This line drops any column containing a missing value
# select(where(~!any(is.na(.)))) 

# To be used as filtering vector in the next section 
SpeciesList <- Traits_full$specie

# # Regarding with the patches ----

M<- M %>%
  filter(!locality %in% c("Ref_Lances", "Guadalete") & !dominantSpecies %in% c("Caulerpa", "Cymodocea", "Cymodocea-PROT", "Rugulopteryx") & habitat != "SEA") 

Manchas <- M%>%
  mutate(habitat = case_when(
    habitat == "COASTLINE" ~ "Strandline",
    habitat == "COAST" ~ "Strandline",
    TRUE ~ habitat
  ))

Manchas$dominantSpecies
Manchas_specie <- Manchas%>%
  mutate(specie = case_when(
    dominantSpecies=="AMM"~ "Ammophila arenaria",
    dominantSpecies=="AMM-PROT"~ "Ammophila arenaria",
    dominantSpecies=="ATRIP_POSTR"~ "Atriplex halimus",
    dominantSpecies=="ATRIP_PORTU"~ "Atriplex halimus",
    dominantSpecies =="ART" ~"Arthemisia maritima",
    dominantSpecies =="ARTH" ~"Arthrocnemum macrostachyum",
    dominantSpecies =="ILEX" ~"Ilex aquifoliu",
    dominantSpecies =="BOLB" ~"Bolboschoenus maritimus",
    dominantSpecies == "CAK" ~ "Cakile maritima",
    dominantSpecies == "CAK-PROT" ~ "Cakile maritima",
    dominantSpecies =="CALY" ~"Calystegia soldanella",
    dominantSpecies =="CARPO" ~"Carpobrotus acinaciformis",
    dominantSpecies =="CENT" ~"Centaurea seridis",
    dominantSpecies =="CRIT" ~"Crithmum maritimum",
    dominantSpecies =="CRUC" ~"Crucianella maritima",
    dominantSpecies =="CIS" ~"Cistaceae",
    dominantSpecies =="ECHIN" ~"Echinophora spinosa",
    dominantSpecies =="ECHINO" ~"Echinophora spinosa",
    dominantSpecies =="ELYM" ~"Elymus farctus",
    dominantSpecies =="EQUI" ~"Equisetum arvense",
    dominantSpecies =="EPH" ~"Ephedra spp.",
    dominantSpecies =="EUPH" ~"Euphorbia paralias",
    dominantSpecies =="EUPH-PROT" ~"Euphorbia parabis",
    dominantSpecies =="ERYNG" ~"Eryngium maritimum",
    dominantSpecies =="GLAU" ~"Glaucium flavum",
    dominantSpecies =="HELY" ~"Helichrysum stoechas",
    dominantSpecies =="JUNCA" ~"Juncus acutus",
    dominantSpecies =="JUNM" ~"Juncus maritimus",
    dominantSpecies =="LOT" ~"Lotus creticus",
    dominantSpecies =="MED" ~"Medicago maritima",
    dominantSpecies =="PANC" ~"Pancratium maritimum",
    dominantSpecies =="PHRAG" ~"Phragmites australis",
    dominantSpecies =="PLANT" ~"Plantago crassifolia",
    dominantSpecies =="POLY" ~"Polygonum maritimum",
    dominantSpecies =="TAM" ~"Tamarix spp.",
    dominantSpecies =="THYM" ~"Thymelaea hirsuta",
    dominantSpecies =="SPO" ~ "Sporobolus pungens",
    dominantSpecies =="SALI" ~ "Salicornia europaea",
    dominantSpecies =="SARCO" ~ "Sarcocornia fruticosum",
    dominantSpecies == "SALS" ~ "Salsola kali",
    dominantSpecies =="SUA" ~"Suaeda maritima",
    dominantSpecies =="XAN" ~ "Xanthium strumarium", #sinonimo de X.echinatum
    dominantSpecies =="GESP" ~ "Xanthium strumarium",
    dominantSpecies =="LAG" ~ "Lagurus ovatus",
    dominantSpecies =="LIM" ~ "Limoniastrum monopetalum",
    dominantSpecies == "LIMO" ~ "Limonium vigoi",
    dominantSpecies == "LOB" ~ "Lobularia app.",
    dominantSpecies =="ROM" ~"Romulea spp.",
    dominantSpecies =="SCIR" ~ "Scirpus maritimus",
    dominantSpecies =="SCOL" ~ "Scolymus hispanicus",
    dominantSpecies =="SED" ~ "Sedum",
    dominantSpecies == "SPERG" ~ "Spergularia marina",
    dominantSpecies == "SPO" ~ "Sporobolus pungens",
    dominantSpecies == "PINUS" ~ "Pinus pinea",
    dominantSpecies =="Caulerpa" ~ "Caulerpa sp.",
    dominantSpecies =="Cymodocea" ~ "Cymodocea nodosa",
    dominantSpecies =="Cymodocea-PROT" ~ "Cymodocea nodosa",
    dominantSpecies == "Rugulopteryx" ~ "Rugulopteryx okamurae",
    TRUE ~ dominantSpecies
  ))

#standarization & outliers! 
Manchas_specie$area <- as.numeric(Manchas_specie$area) #square meters
Manchas_specie$elevation <- as.numeric(Manchas_specie$elevation)
Manchas_specie$elevation <- Manchas_specie$elevation * 100 #elevation data from m to cm

Manchas_specie  <- Manchas_specie  %>%
  filter(elevation<100) %>% #elevation < 100 centimeters
  filter(area <10) #area < 10 meters

# The full data set with traits + patch data----
manchas_traits_all <- Manchas_specie %>% 
  filter(!is.na(elevation)) %>% 
  filter(specie %in% SpeciesList) %>%
  left_join(Traits_full, by = "specie") 

manchas_traits_all1<- manchas_traits_all[,-c(1:6,8,9,12,14,18)]  #drop categorical variables (less spp-habitat)

# Looped version  ----
predictor.var.out <- NULL
outcome.var.out <- NULL
bestTune <- NULL
bestTune.out <- NULL
r2 <- NULL
r2.out <- NULL
varImp <- NULL
imp.out <- NULL
coef1.out <- NULL
coef2.out <- NULL
manchas_traits <- as.data.frame(manchas_traits_all1)
scoresPlsr1.out <- NULL
for(i in 1:100){
  # We’ll randomly split the data into training set (80% for building a predictive model) 
  # and test set (20% for evaluating the model). Make sure to set seed for reproducibility.
  # Split the data into training and test set
  # set.seed(123)
  training.samples <- createDataPartition(manchas_traits$elevation, p = 0.8, list = FALSE)
  train.data  <- manchas_traits[training.samples, ]
  test.data <- manchas_traits[-training.samples, ]
  
  # The R function train() [caret package] provides an easy workflow to compute PCR and PLS by invoking the pls package.
  # It has an option named method, which can take the value pcr or pls.
  # An additional argument is scale = TRUE for standardizing the variables to make them comparable.
  # caret uses cross-validation to automatically identify the optimal number of principal components (ncomp) to be incorporated in the model.
  # Here, we’ll test 10 different values of the tuning parameter ncomp. 
  # This is specified using the option tuneLength. 
  # The optimal number of principal components is selected so that the cross-validation error (RMSE) is minimized.
  
  # Build the model on training set
  # set.seed(123)
  fitControl <- trainControl(method = "cv",  number = 10)
  model.initial <- train(elevation ~ area+canopy+TotalLength+ CAN_H + LeavesOpening +
                           RootOpening + ShootDensity_m2+RootDensity_m2+MaximumRootLength+TotalRootLength_cm_m2+MeanRootLength+
                           AbovegroundBiomass_g_m2+SRL+RWR+RS+RLR+RTD+DLR+
                           AboveBelowBiomass + RootBiomass_g_m2+LeavesBiomass_g_m2+StemBiomass_g_m2,
                         data = train.data, 
                         method = "simpls",
                         scale = TRUE,
                         trControl = fitControl,
                         tuneLength = 27)
  
  # Print the best tuning parameter ncomp that minimizes the cross-validation error, RMSE
  bestTune <- model.initial$bestTune
  
  # Summarize the final model
  summary(model.initial$finalModel)
  predictor.var <- str_extract(capture.output(summary(model.initial))[7], "[0-9]+")
  outcome.var <- str_extract(capture.output(summary(model.initial))[8], "[0-9]+")
  
  model.coefficients <- as_tibble(model.initial$finalModel$coefficients)
  model.coefficients.tibble <- tibble(predictors = rownames(model.initial$finalModel$coefficients), model.coefficients)
  coef1 <- model.coefficients.tibble$`.outcome.1 comps`
  
  # Model loadings
  model.loadings <- model.initial$finalModel$loadings
  model.loadings1.tibble <- tibble(predictors = rownames(model.initial$finalModel$loadings), model.initial$finalModel$loadings[,1])
  loadings1 <- as.numeric(model.loadings1.tibble$`model.initial$finalModel$loadings[, 1]`)
  # if(bestTune > 1){
  #   coef2 <- model.coefficients.tibble$`.outcome.2 comps`
  # }
  
  # Variable Importance on the Projection
  varIMP <- varImp(model.initial, scale = FALSE)
  varIMP
  # Model performance metrics
  predictions <- model.initial %>% predict(test.data)
  r2 <- caret::R2(predictions, test.data$elevation)
  #  Rsquare
  #  0.4289
  
  # # Make predictions an plot them in 1:1 plots to see residual variation
  test.data$predictions <- model.initial %>% 
    predict(test.data)
  p <- ggplot(data = test.data, mapping = aes(x = predictions, y = elevation)) +
    geom_point(aes(colour = specie)) +
    scale_colour_d3("category20", 
                    alpha = 1,
                    guide = guide_legend(override.aes = list(size = 3,
                                                             alpha = 1))) +
    geom_abline() +
    ylab("Observed sediment accumulation") +
    xlab("Predicted sediment accumulation (PLSR)") + 
    theme_bw() +
    theme(
      text = element_text(size = 18),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16),
      legend.box.spacing = unit(0.2, "cm"),  #espaciado leyenda
      #   legend.key.size = unit(1, "lines"),   # tamaño simbolos en la leyenda
      legend.title = element_blank(),
      legend.text = element_text(size = 16, face = "italic")   
    ) 
  
  p + annotate(geom = "text", x = 10, y = 45, label = bquote(R^2 == .(round(r2, 2))))
  
  setwd(plots.dir) 
  dir()
  name<-'7.22PLRS_species-Standarized1' 
  ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
  ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)
  
  
  # Predictions are not super good, but overall the species effect is very obvious. With some intraspecific residual variance
  # that it's escaping our model (we're not managing to predict this intraspecific variance).
  
  if(i == 1){
    imp.out <- varIMP$importance
    loadings1.out <- loadings1
    # if(bestTune > 1){
    #   coef2.out <- coef2
    # }
  }
  else{
    imp.out <- cbind(imp.out, varIMP$importance)
    loadings1.out <- cbind(loadings1.out, loadings1)
    # if(bestTune > 1){
    #   coef2.out <- cbind(coef2.out, coef2)
    # }
  }
  r2.out <- c(r2.out, r2)
  bestTune.out <- c(bestTune.out, bestTune[[1]])
  predictor.var.out <- c(predictor.var.out, predictor.var)
  outcome.var.out <- c(outcome.var.out, outcome.var)
  if(i<20){
    train.data$scoresPlsr1 <- model.initial$finalModel$scores[1:353,1]
    ggplot(data = train.data, mapping = aes(x = scoresPlsr1, y = elevation)) +
      geom_point(aes(colour = specie)) +
      geom_smooth(method = "lm", se=T, colour = "black") + 
      theme_bw() +
      theme(text = element_text(size=14),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    
  }
}

setwd(plots.dir) 
dir()
name<-'7.2 PLSR2' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)

# Variable importance
names(imp.out) <- 1:100
# imp.out <- as.data.frame(t(imp.out))
# imp.means <- sort(colMeans(imp.out), decreasing = T)
# imp.sds <- colSds(as.matrix(imp.out))
imp.means <- rowMeans(imp.out)
imp.sds <- rowSds(as.matrix(imp.out))
imp.summarised <- as.data.frame(cbind(imp.means, imp.sds))
imp.summarised$x <- rownames(imp.summarised)


# Loadings 1st component
loadings1.out <- as.data.frame(loadings1.out)
rownames(loadings1.out) <- rownames(imp.summarised)
names(loadings1.out) <- 1:100
loadings1.means <- rowMeans(loadings1.out)
loadings1.sds <- rowSds(as.matrix(loadings1.out))
loadings1.summarised <- as.data.frame(cbind(loadings1.means, loadings1.sds))
loadings1.summarised$x <- rownames(imp.summarised)


# Best tune
mean(bestTune.out) # 3.8 components
sd(bestTune.out) # 1.263313 components
min(bestTune.out) # 3 components
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(bestTune.out) # 3 components is the most frequent one

bestTune.out <- as.data.frame(bestTune.out)

max(bestTune.out) #9 max

bestTune.out %>% 
  mutate(name = as.factor(bestTune.out)) %>% 
  group_by(name) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = name, y = n)) +
  geom_bar(stat = 'identity')

bestTune.out %>% 
  mutate(name = as.factor(bestTune.out)) %>% 
  group_by(name) %>% 
  summarise(n = n()) %>%
  mutate(percent = (n/sum(n))*100) %>% 
  ggplot(aes(x = name, y = percent)) +
  geom_bar(stat = 'identity') +
  ylim(c(0,100)) + 
  ylab("%") +
  xlab('Number of components') + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# R2
mean(r2.out) # R2 = 0.2329885
sd(r2.out) # 0.08911194
min(r2.out) # 0.08373869
max(r2.out) # 0.4513046

# # Scores PLSR1
names(scoresPlsr1.out) <- 1:100
scoresPlsr1.out <- as.data.frame(scoresPlsr1.out)
scoresPlsr1.means <- rowMeans(scoresPlsr1.out)
scoresPlsr1.sds <- rowSds(as.matrix(scoresPlsr1.out))
scoresPlsr1.summarised <- as.data.frame(cbind(scoresPlsr1.means, scoresPlsr1.sds))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# GRAPHS FOR THE PAPER HIGHLIGHTING MODEL RESULTS FOR SOIL D13C ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# # Variable importance on the projection
imp.summarised <- imp.summarised %>% 
  mutate(fill = as.numeric(imp.means < 1))

ggplot(imp.summarised, aes(x=reorder(x, imp.means, FUN = max), y = imp.means, fill = fill)) + 
  geom_bar(stat="identity", colour = "black", show.legend = F) +
  geom_errorbar(aes(ymin=imp.means, ymax=imp.means+imp.sds), width=.2, position=position_dodge(.9)) + 
  # scale_y_continuous(name = "Variable importance", limits = c(0,25,50,75,100)) +
  xlab(label = "") + 
  ylab(label ="Variable Importance in Projection (%)") +
  coord_flip() + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

setwd(plots.dir) 
dir()
name<-'7.2 PLSR3-Standarized1' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)

# Another plot for Variable importance
imp_plot <- ggplot(data = pivot_longer(as_tibble(t(imp.out)), cols = 1:21), aes(y = value, x = reorder(name, value, FUN = median))) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 1, lty = 2) +
  xlab(label = "") + 
  ylab(label ="Variable Importance in Projection") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16),
        legend.box.spacing = unit(0.2, "cm"),  #espaciado leyenda
        #   legend.key.size = unit(1, "lines"),   # tamaño simbolos en la leyenda
        legend.title = element_blank(),
        legend.text = element_text(size = 16, face = "italic")   
  ) 
imp_plot

setwd(plots.dir) 
dir()
name<-'7.2VarIMpoPLSR4-Standarized1' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)


imp.summarised <- pivot_longer(as_tibble(t(imp.out)), cols = 1:22) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(fill = as.numeric(mean < 1))
imp.summarised

# PLSR Loadings 1st component
load_plot <- loadings1.summarised %>% 
  left_join(imp.summarised, by = c("x" = "name")) %>% 
  # filter(fill == 0) %>% 
  ggplot(aes(x = reorder(x, loadings1.means, FUN = max), y = loadings1.means, fill = fill)) + 
  geom_bar(stat="identity", colour = "black", show.legend = F) +
  geom_errorbar(aes(ymin=loadings1.means - loadings1.sds, ymax = loadings1.means + loadings1.sds), width=.2, position=position_dodge(0.9)) + 
  ylab("PLSR loadings 1st component") +
  xlab("") + 
  # scale_y_continuous(limits = c(-2, 2)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16),
        legend.box.spacing = unit(0.2, "cm"),  #espaciado leyenda
        #   legend.key.size = unit(1, "lines"),   # tamaño simbolos en la leyenda
        legend.title = element_blank(),
        legend.text = element_text(size = 16)   
  ) 
load_plot

setwd(plots.dir) 
dir()
name<-'2load_plot-Standarized1' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)


# Another plot for Loadings importance
ggplot(data = pivot_longer(as_tibble(t(loadings1.out)), cols = 1:22), aes(y = value, x = reorder(name, value, FUN = median))) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab(label = "") + 
  ylab(label ="Variable Importance in Projection") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16),
        legend.box.spacing = unit(0.2, "cm"),  #espaciado leyenda
        #   legend.key.size = unit(1, "lines"),   # tamaño simbolos en la leyenda
        legend.title = element_blank(),
        legend.text = element_text(size = 16)   
  ) 

setwd(plots.dir) 
dir()
name<-'2varIMpo--Standarized1' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)


# Figures PLSR scores vs. RESPONSE VARIABLE (similar to Carrascal et al 2009)
# We use just the last result from the for loop as an example
train.data$scoresPlsr1 <- model.initial$finalModel$scores[1:353,1]
plsr_line_plot <- ggplot(data = train.data, mapping = aes(x = scoresPlsr1, y = elevation)) +
  geom_point(aes(colour = specie)) +
  # scale_colour_simpsons(alpha = 1) +
  scale_colour_d3("category10", 
                  alpha = 0.5,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  geom_smooth(method = "lm", se=T, colour = "black") + 
  xlab("PLSR scores 1st component") + 
  ylab(expression(Delta*'Elevation (cm)')) +
  theme_bw()+
  theme(
    text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.box.spacing = unit(0.2, "cm"),  #espaciado leyenda
    #   legend.key.size = unit(1, "lines"),   # tamaño simbolos en la leyenda
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "italic")   
  ) 


plsr_line_plot

setwd(plots.dir) 
dir()
name<-'7.22speicesPLSRpp--Standarized1' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)


library(cowplot)
panel1 <- plot_grid(imp_plot, 
                    load_plot + geom_hline(yintercept = 0),
                    nrow = 1, labels = "AUTO", rel_widths = c(1,1), align = "h")

panel2 <- plot_grid(NULL,
                    plsr_line_plot,
                    NULL,
                    ncol = 3, rel_widths = c(1,3,1), labels = c("", "C", ""))
plot_grid(panel1,
          panel2,
          nrow = 2, rel_heights = c(3,2), align = "hv")

setwd(plots.dir) 
dir()
name<-'Final2-Standarized1' 
ggsave(paste(name,".png",sep="."),width = 31, height = 30,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 31, height = 30,units = "cm",dpi=1600,useDingbats=FALSE)

# Final model summary results
summary(model.initial)

#Data: 	X dimension: 353 23 
#Y dimension: 353 1
#Fit method: simpls
#Number of components considered: 3
#TRAINING: % variance explained
#  1 comps  2 comps  3 comps
#X          35.34    46.73    55.61
#.outcome    11.39    19.77    21.33