# # # # # # Mixed models analysis - Model 1 # # # # # #

# # # 1st PhD chapter - Rodríguez-Arias,L.

# 12/2023

rm (list=ls()) 

options(contrasts=c(factor="contr.sum", ordered="contr.poly"))

library(easypackages)
libraries("lme4","nlme","lmerTest","ggplot2","sjPlot" ,"MASS","jtools","ggthemes", "tidyverse",
          "performance", "multcomp","car", "dplyr", "tidyr", "readr", "ggstatsplot",
          "statsExpressions", "scales", "MuMIn", "car", "DHARMa", "sandwich", 
          "visreg","fitdistrplus","lmtest", "ggpubr", "corrplot","reghelper", "effects")

Theme2<- theme(
  axis.title = element_text(size = 18, face = "bold"),
  axis.text = element_text(size = 18, face = "plain"), 
  axis.title.x = element_text(size = 18),  
  axis.title.y = element_text(size = 18), 
  legend.text = element_text(size = 18),  
  legend.title = element_text(size = 15),
  strip.text = element_text(size = 15, face = "italic"))  

# Directory ----

work.dir=("C:/Users/lucia/Desktop/PhD_Blanes/3_R/Paper1/Paper1_mechanism_stabilization")

data.dir=paste(work.dir,"Data",sep="/")
plots.dir=paste(work.dir,"Plots",sep="/")
model.outs=paste(work.dir,"ModelOuts",sep="/")

# # Raw data object ----
setwd(data.dir)
dir()
Manchas <- read.csv2("manchasR.csv")
summary(Manchas) 
str(Manchas)
glimpse(Manchas)

#Functions
n_fun <- function(x) {
  return(data.frame(y = max(x), label = paste0("n = ", length(x))))
}

source("Functions.R")

# variable requalification ----

#Both the AREA and the ELEVATION are in square meters and meters respectively in the database
Manchas$area <- as.numeric(Manchas$area) #square meters
Manchas$elevation <- as.numeric(Manchas$elevation)
Manchas$elevation <- Manchas$elevation * 100 #elevation data from m to cm

system=as.factor(Manchas$system)
locality=as.factor(Manchas$locality)
HIC=as.factor(Manchas$HIC)
habitat=as.factor(Manchas$habitat)

ManchasDelta <- Manchas %>%
  filter(!locality %in% c("Ref_Lances", "Guadalete") & !dominantSpecies %in% c("Caulerpa", "Cymodocea", "Cymodocea-PROT", "Rugulopteryx") & habitat != "SEA")

# # Scientific name for all the spp ---- 

ManchasDelta$dominantSpecies  
Manchas_habitat <- ManchasDelta%>%
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
    dominantSpecies =="EUPH-PROT" ~"Euphorbia paralias",
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

Manchas_habitat <- Manchas_habitat%>%
  mutate(habitat = case_when(
    habitat == "COASTLINE" ~ "Strandline",
    habitat == "COAST" ~ "Strandline",
    TRUE ~ habitat
  ))

table(Manchas_habitat$habitat)


# Selection variables of interest for models ----
variablesMorfometria <- dplyr::select(Manchas_habitat, -long, -lat, -alt, -dominantSpecies)
Morfometria <- na.omit(variablesMorfometria) 

# Outliers Consideration
# We filter elevations of more than 1 meter and spots of area more than 10 meters
Morfometria <- Morfometria  %>%
  filter(elevation<100) %>% #elevation < 100 centimetros
  filter(area <10) #area < 10 metros

# We filter at least 8 observations in the record for each plant species to be taken into account in the model
numero <- table(Morfometria$specie)
especies_filtradas <- names(numero[numero>=8])
Morfometria <- subset(Morfometria, specie %in% especies_filtradas)  

Morfometria <- Morfometria %>% 
  mutate(locality2 = gsub("[[:digit:]]", "", locality))


# Colineality. We select the uncorrelated predictor variables (R2>0.7)
#numerical variables
correlacion=cor(Morfometria[, c("area","elevation","distance")]) 
print(correlacion)
corrplot(correlacion, method="circle") 

#categorical variables
predictors <- Morfometria[, c("area", "canopy", "multiMono", "specie", "speciesRichness", "HIC", "habitat", "system", "distance")]
predictors <- model.matrix(~., data = predictors)[, -1] # categorical variables to dummy variables
cor_matrix <- cor(predictors)
print(cor_matrix)
corrplot(cor_matrix, method = "circle")
#multiMono-speciesRichness high collinearity

###########################################################################################
plotdist(Morfometria$elevation,histo = TRUE, demp = TRUE) 
descdist(Morfometria$elevation,discrete = FALSE, boot = 1000)

dists <- c("gamma", "lnorm", "weibull") 
fit <- list()
for (i in 1:length(dists)){
  fit[[i]] <- fitdist(Morfometria$elevation, dists[i])
}

for (i in 1:length(dists)){
  print(summary(fit[[i]]))
}

plot.legend <- dists
denscomp(fit, legendtext = plot.legend, fitlwd = 3) 
cdfcomp(fit, legendtext = plot.legend, datapch = 21, fitlwd = 3)
qqcomp(fit, legendtext = plot.legend) 
ppcomp(fit, legendtext = plot.legend)

gofstat(fit, fitnames = dists)
#We determine that the lnorm is the distribution that best fits the data; in addition to having a lower AIC value

##################################################################################################
#do we really need the random part?
Mod_lineal <- nlme::gls(log(elevation) ~ log(area) + log(canopy) + multiMono + specie + habitat+ system+distance, 
                       data = Morfometria)
Mod_mixto <- nlme::lme(log(elevation) ~ log(area) + log(canopy) + multiMono + specie + habitat+ system+distance,
                        random = ~1 | locality2, data = Morfometria)

#m1 <- lme(log(elevation) ~ log(area) + log(canopy) + multiMono + specie + habitat,
#          random = ~1|locality2, 
#          weights = varIdent(form = ~1|specie),
#          method = "REML",
#          control = lmeControl(maxIter = 1000, msMaxIter = 1000),
#          data=Morfometria)
#summary(m1)
#anova(m1) #just to check 

anova(Mod_lineal, Mod_mixto)
# Preliminary test of mixed model versus linear model: the random effect (lower AIC & p-value <.0001)

##################################################################################################

#checking the fixed-random part of the models, and removing interactions 

Mod_mixto_all <- lmer(log(elevation) ~ log(area) + log(canopy) + multiMono + specie+ system+distance+
                        (1 | locality2), data = Morfometria,  na.action = "na.fail", REML = FALSE)

Mod_mixto_nested <- lmer(log(elevation) ~ log(area) + log(canopy) + multiMono + specie+ system+distance+
                           (1 | locality2/habitat), data = Morfometria,  na.action = "na.fail", REML = FALSE)

Mod_mixto_random <- lmer(log(elevation) ~ log(area) + log(canopy) + multiMono + specie+ system+distance+
                           (1 | locality2) + (1|habitat), data = Morfometria,  na.action = "na.fail", REML = FALSE)

#Checking random model part:
AIC(Mod_mixto_nested, Mod_mixto_random)
AICc(Mod_mixto_nested, Mod_mixto_random)
lrtest(Mod_mixto_nested, Mod_mixto_random)
anova(Mod_mixto_nested, Mod_mixto_random) #"Mod_mixto_nested" has no sense,'cause it was a non-nested experimental design

AIC(Mod_mixto_random, Mod_mixto_all)
AICc(Mod_mixto_random, Mod_mixto_all)
lrtest(Mod_mixto_random, Mod_mixto_all)
anova(Mod_mixto_random, Mod_mixto_all) #they're similar, but
#we chose the mixed model, since there is no improvement in the model with the other, nor is there any difference between them.
#It has more parsimony as it is a simpler model, without losing explanatory power.

#Finally, we choose "Mod_mixto_all" (lower AIC & BIC)

# Checking the fixed part:
Modelo_mixto1 <- lmerTest::lmer(log(elevation) ~ log(area) + log(canopy) + multiMono + specie+ habitat+system+distance+
                        (1 | locality2), data = Morfometria,  na.action = "na.fail", REML = FALSE)
Mod_mixto2 <- lmerTest::lmer(log(elevation) ~ log(area) + multiMono + specie*habitat +system+distance+
                              (1 | locality2) , data = Morfometria,  na.action = "na.fail", REML = FALSE)

AIC(Modelo_mixto1,Mod_mixto2)
AICc(Modelo_mixto1,Mod_mixto2)
lrtest(Modelo_mixto1,Mod_mixto2)
anova(Modelo_mixto1,Mod_mixto2) 

summary(Mod_mixto2)
car::Anova(Mod_mixto2, type = 3)
car::vif(Mod_mixto2)
#Potential issues with multicollinearity: specie, habitat, and the interaction term specie:habitat exhibit high VIF values
check_model(Mod_mixto2) #alternative ggplot
check_normality(Mod_mixto2)
check_heteroscedasticity(Mod_mixto2)
check_autocorrelation(Mod_mixto2)
hist(resid(Mod_mixto2), breaks = 20)

simulateResiduals(fittedModel = Mod_mixto2, plot=T)   

## Mod_mixto2 is good, but it exist multicollinearity problems.
#We selected the model 'Mod_mixto_all' as a result of AIC, AICc and Likelihood Ratio Tests
# + non multicollinearity problems
#Mod_mixto_all: log(elevation) ~ log(area) + multiMono + specie + system + distance + (1 | locality2)

##################################################################################################

Modelo_initial <- lmerTest::lmer(log(elevation) ~ log(area) + log(canopy) + multiMono + specie+ system+distance+
                               (1 | locality2), data = Morfometria,  na.action = "na.fail", REML = FALSE)

summary(Modelo_initial)
car::Anova(Modelo_initial, type = 3)
car::vif(Modelo_initial)
#The correlation matrix between variables is shown as there are no multicollinearity problems.

check_model(Modelo_initial)
step(Modelo_initial) #drop variables!

Morfometria <- tidyr::drop_na(Morfometria)  


Morfometria$specie <- factor(Morfometria$specie, 
                             levels = c("Limoniastrum monopetalum","Ammophila arenaria", "Limonium vigoi",  
                                         "Salsola kali", "Euphorbia paralias","Sporobolus pungens",
                                        "Sarcocornia fruticosum","Echinophora spinosa",
                                         "Cakile maritima","Salicornia europaea","Suaeda maritima","Xanthium strumarium"),
                             labels = c("Limoniastrum\nmonopetalum","Ammophila\narenaria", "Limonium\nvigoi",
                                         "Salsola\nkali","Euphorbia\nparalias","Sporobolus\npungens",
                                         "Sarcocornia\nfruticosum","Echinophora\nspinosa", 
                                        "Cakile\nmaritima","Salicornia\neuropaea","Suaeda\nmaritima","Xanthium\nstrumarium"))


Modelo_final <- lmerTest::lmer(log(elevation) ~ log(area) + multiMono + specie + (1 | locality2),
                                   data = Morfometria,  na.action = "na.fail", REML = TRUE)

#summary(Modelo_final)
print(summary(Modelo_final), correlation = TRUE)

car::Anova(Modelo_final, type = 3)
#Response: log(elevation)
#              Chisq Df Pr(>Chisq)    
#(Intercept) 715.693  1  < 2.2e-16 ***
#  log(area)   169.137  1  < 2.2e-16 ***
#  multiMono   23.210  1  1.452e-06 ***
#  specie      31.327 11  0.0009773 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

car::vif(Modelo_final) #ok
performance::check_model(Modelo_final) 
#it's seems to have some normality problems, but it is not that bad or problematic for the consistency of our results
mcheck(Modelo_final)
check_normality(Modelo_final)
check_heteroscedasticity(Modelo_final) #ok
#Regression models where the variance of the residuals (or errors) is constant across all levels of the independent variables.
#high p-value (usually above the 0.05 threshold) suggests that there is no evidence to reject the null hypothesis of homoscedasticity
check_autocorrelation(Modelo_final)
hist(resid(Modelo_final), breaks = 20)

plot(Modelo_final)
simulateResiduals(fittedModel = Modelo_final, plot=T) 
#No significant problems detected (Q-Q plot residuals & residual vs predicted values)


#R-squared values (Nakagawa and Schielzeth, 2013) for mixed models
#Conditional R-squared (proportion of variance explained by both fixed and random effects) 
#Marginal R-squared (proportion of variance explained by fixed effects only)
performance::r2_nakagawa(Modelo_final)
# R2 for Mixed Models

#Conditional R2: 0.462
#Marginal R2: 0.389

## Visualizing Model Results
interactions <- allEffects(Modelo_final) 
plot(interactions)
axis(side = 1, las = 1) 

######################Graphics - model output ----
visreg(Modelo_final,xvar = "specie", trans=exp,ylab=expression(Delta*'Elevation (cm)'), #~' (log scale)'
       gg = T, partial = T, line=list(col="black", size =.51)) +
  xlab(bquote("")) +
  theme_classic()+ 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 18),
        axis.text.x = element_text(face = "italic", 
                                   size = 18,
                                   angle = 55, 
                                   vjust = 1, 
                                   hjust = 1))+
  Theme2


setwd(plots.dir) 
dir()
name<-'1.PAPER_SuplMat_OutcomeGLMM-spp' 
ggsave(paste(name,".png",sep="."),width = 31, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 31, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)


# # Area 

visreg(Modelo_final,xvar ="area", trans=exp, ylab=expression(Delta*'Elevation (cm)'), gg = T, partial = T, line=list(col="black", size = 0.5)) +
  xlab(bquote("Patch area ("*m^2*")")) +
  theme_classic()+ 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 18),
        axis.text.x = element_text(size = 18))+
  Theme2

setwd(plots.dir) 
dir()
name<-'1.PAPER_SuplMat_Outcome-area' 
ggsave(paste(name,".png",sep="."),width = 31, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 31, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)

# # Multimono ---- 
Morfometria$multiMono <- factor(Morfometria$multiMono, levels = c("Mono", "Multi"),
                                labels = c("Monospecific", "Multispecific"))

visreg(Modelo_final,xvar ="multiMono",trans=exp, ylab=expression(Delta*'Elevation (cm)'), gg = T, partial = T, line=list(col="black", size = 0.5)) +
  xlab(bquote("")) +
  theme_classic()+ 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 18),
        axis.text.x = element_text(size = 18))+
  Theme2

setwd(plots.dir) 
dir()
name<-'1.PAPER_SuplMat_OutcomemultiMono' 
ggsave(paste(name,".png",sep="."),width = 30, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 30, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)

##################################################################################################
# Multiple comparisons with library(multcompar) FOR SPECIES & SPECIES COMPOSTION

#Tuket test for specie (multicomp)
SpecieSignificance <- multcomp::glht(Modelo_final, linfct=multcomp::mcp(specie="Tukey"))
summary_spp <- broom::tidy(summary(SpecieSignificance))
letters_spp <- broom::tidy(multcomp::cld(SpecieSignificance))
# Significant multiple comparisons
Tukey_hab <- summary_spp %>%
  filter(p.value<0.05) %>%
  print(n = Inf)


# A tibble: 2 × 6
#lhs                                       rhs estimate std.error statistic p.value
#<chr>                                   <dbl>    <dbl>     <dbl>     <dbl>   <dbl>
#  1 Cakile maritima - Ammophila arenaria        0   -0.409     0.113     -3.63 1.25e-2
#  2 Salicornia europaea - Ammophila arenar…     0   -0.496     0.116     -4.29 9.90e-4

#Tuket test for mono or mixed patches (multicomp)
MultiSignificance <- multcomp::glht(Modelo_final, linfct=multcomp::mcp(multiMono="Tukey"))
summary_multi <- broom::tidy(summary(MultiSignificance))
letters_multi <- broom::tidy(multcomp::cld(MultiSignificance))
# Significant multiple comparisons
Tukey_multi <- summary_multi %>%
  filter(p.value<0.05) %>%
  print(n = Inf)

# A tibble: 1 × 6
#lhs            rhs estimate std.error statistic  p.value
#<chr>        <dbl>    <dbl>     <dbl>     <dbl>    <dbl>
#  1 Multi - Mono     0    0.498     0.103      4.82 0.00000145

###########################################################################################
######################Graphics - data en bruto ----
library(ggsci)
library(ggplot2)
library(gridExtra)
compare_means(elevation ~area, data = Morfometria, ref.group = ".all.")


compare_means(elevation ~ specie,  data = Morfometria, ref.group = ".all.",
              method = "t.test")
# A tibble: 12 × 8
#.y.       group1 group2                     p  p.adj p.format p.signif method
#<chr>     <chr>  <chr>                  <dbl>  <dbl> <chr>    <chr>    <chr> 
#  1 elevation .all.  "Limoniastrum\nmonopeta… 1.03e- 1 5.3 e- 1 0.1032   ns       T-test
#2 elevation .all.  "Ammophila\narenaria"    1.38e- 3 1.4 e- 2 0.0014   **       T-test
#3 elevation .all.  "Limonium\nvigoi"        3.95e-24 4.70e-23 < 2e-16  ****     T-test
#4 elevation .all.  "Salsola\nkali"          2.96e- 1 1   e+ 0 0.2962   ns       T-test
#5 elevation .all.  "Euphorbia\nparalias"    5.77e- 8 6.30e- 7 5.8e-08  ****     T-test
#6 elevation .all.  "Sporobolus\npungens"    7.65e- 1 1   e+ 0 0.7654   ns       T-test
#7 elevation .all.  "Sarcocornia\nfruticosu… 3.09e- 1 1   e+ 0 0.3085   ns       T-test
# 8 elevation .all.  "Echinophora\nspinosa"   1.95e- 3 1.8 e- 2 0.0019   **       T-test
# 9 elevation .all.  "Cakile\nmaritima"       2.46e- 2 2   e- 1 0.0246   *        T-test
#10 elevation .all.  "Salicornia\neuropaea"   5.10e- 1 1   e+ 0 0.5103   ns       T-test
#11 elevation .all.  "Suaeda\nmaritima"       5.82e- 2 4.1 e- 1 0.0582   ns       T-test
#12 elevation .all.  "Xanthium\nstrumarium"   8.76e- 2 5.3 e- 1 0.0876   ns       T-test       T-test

compare_means(elevation ~ multiMono,  data = Morfometria, ref.group = ".all.",
              method = "t.test")

# A tibble: 2 × 8
#.y.       group1 group2       p  p.adj p.format p.signif method
#<chr>     <chr>  <chr>    <dbl>  <dbl> <chr>    <chr>    <chr> 
#1 elevation .all.  Mono   0.127    0.13    0.12702  ns       T-test
#2 elevation .all.  Multi  0.000110 0.00022 0.00011  ***      T-test

ggboxplot(Morfometria,x="specie", y="elevation", color="specie", add = "jitter",
          facet.by = "habitat", short.panel.labs = TRUE)+
  rotate_x_text(angle = 45)+
  stat_compare_means(method = "anova", label.y = 40)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")              # Pairwise comparison against all

ggboxplot(Morfometria,x="specie", y="elevation", color="specie", add = "jitter",
          facet.by = "multiMono",
          short.panel.labs = TRUE)+
  rotate_x_text(angle = 45)+
  stat_compare_means(method = "anova", label.y = 40)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")              # Pairwise comparison against all

Specie_letters<- Morfometria%>%
  left_join(letters_spp, by = "specie")
Specie_letters<- Specie_letters[,c(4,15,17)]

ggplot(Specie_letters, aes(x = specie, y = elevation)) +
  geom_boxplot(aes(fill = specie, colour = specie)) +
  geom_text(data = Specie_letters, aes(x = specie, y = elevation, label = letters), size = 4) 

elevation_box <- Morfometria %>% 
  group_by(specie) %>% 
  summarise(n = n(),
            mean_stock = mean(elevation, na.rm = TRUE),
            median_stock = median(elevation, na.rm = TRUE),
            std.error = sd(elevation, na.rm = TRUE) / sqrt(n),  # Cálculo del error estándar de la media
            max_stock = max(elevation, na.rm = TRUE)) %>% 
  left_join(letters_spp, by = "specie")
 
#possible final graph
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
               "#1a55FF", "#dbb40c")


plot_species <- Morfometria%>%
  left_join(letters_spp, by = "specie")%>%
  ggplot(aes(x = reorder(specie, elevation, FUN = median), y = elevation)) +
  #geom_boxplot(aes(fill = specie, colour = specie),alpha = 0.6)+
  #geom_jitter(aes(fill = specie, colour = specie),width = 0.2, size = 1.5, alpha = 0.6)+
  #scale_colour_d3("category20") +
  #scale_fill_npg()+
  #scale_fill_d3("category10") +
  geom_boxplot(fill = "#828282", colour = "#828282", alpha = 0.6) +
  stat_summary(fun = median, geom = "crossbar", width = 0.8, lwd = 0.2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 3), size = 5) +
  geom_text(data = elevation_box, aes(x = specie, y = max_stock, label = letters), size = 7.5,position = position_nudge(y = 7.5)) +
  xlab("") +
  ylab(expression(Delta*'Elevation (cm)')) +
  # coord_flip() +
  # facet_wrap(~habitat) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 18),
        axis.text.x = element_text(face = "italic",
                                   color="black",
                                   size=18,
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))
plot_species

setwd(plots.dir) 
dir()
name<-'1.PAPER_BLACK-MulticompSpeciesSignificance' 
ggsave(paste(name,".png",sep="."),width = 25, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 25, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)

library(forcats)


#or this one
my_grey <- c("grey40","grey40","grey40","grey40","grey40","grey40","grey40","grey40","grey40",
             "grey40","grey40","grey40")
my_grey2 <- c("grey70","grey70","grey70","grey70","grey70","grey70","grey70","grey70","grey70",
             "grey70","grey70","grey70")
Morfometria%>%
  left_join(letters_spp, by = "specie")%>%
  ggplot(aes(x = reorder(specie, elevation, FUN = median), y = elevation)) +
  geom_boxplot(aes(fill = specie, colour = specie),alpha = 0.6,outlier.shape = NA, col = NA)+ # fill = NA Sin box
  geom_point(aes(fill = specie, colour = specie),size = 3, shape = 21, alpha = 0.5, 
             position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0.5))+
  theme_few() + 
  scale_colour_manual(values = my_grey2) +
  scale_fill_manual(values= my_grey) +
  stat_summary(fun = median, geom = "crossbar", width = 0.8, lwd = 0.2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 3), size = 5) +
  geom_text(data = elevation_box, aes(x = specie, y = max_stock, label = letters), size = 7.5,position = position_nudge(y = 7.5)) +
  xlab("") +
  ylab(expression(Delta*'Elevation (cm)')) +
  # coord_flip() +
  # facet_wrap(~habitat) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 18),
        axis.text.x = element_text(face = "italic",
                                   color="black",
                                   size=18,
                                   angle = 55, 
                                   vjust = 1, 
                                   hjust = 1))+
  Theme2

setwd(plots.dir) 
dir()
name<-'GREY.2 MulticompSpeciesSignificanceJitt' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)

# Patches composition: mono - mixed
multi_box <- Morfometria %>% 
  group_by(multiMono) %>% 
  summarise(n=n(),
            mean_stock = mean(elevation, na.rm = TRUE),
            median_stock = median(elevation, na.rm = TRUE),
            std.error = sd(elevation, na.rm = TRUE) / sqrt(n),  # Cálculo del error estándar de la media
            max_stock = max(elevation, na.rm = TRUE)) %>% 
  left_join(letters_multi, by = "multiMono")

compare_means(elevation ~ multiMono,  data = Morfometria, ref.group = ".all.",
              method = "t.test")
#A tibble: 2 × 8
#.y.       group1 group2           p  p.adj p.format p.signif method
#<chr>     <chr>  <chr>        <dbl>  <dbl> <chr>    <chr>    <chr> 
# 1 elevation .all.  Mono   0.127    0.13    0.12702  ns       T-test
# 2 elevation .all.  Multi  0.000110 0.00022 0.00011  ***      T-test

Multi_letters<- Morfometria%>%
  left_join(letters_multi, by = "multiMono")

ggplot(Multi_letters, aes(x = multiMono, y = elevation)) +
  geom_boxplot(aes(fill = multiMono, colour = multiMono)) +
  geom_text(data = Multi_letters, aes(x = multiMono, y = elevation, label = letters), size = 4) 

#possible final graph
plot_multi <- Morfometria %>%
  left_join(letters_multi, by = "multiMono") %>%
  ggplot(aes(x = reorder(multiMono, elevation, FUN = median), y = elevation)) +
  geom_boxplot(aes(fill = multiMono, colour = multiMono), alpha = 0.6) +
  scale_fill_manual(values = c("grey70", "grey40")) +
  scale_color_manual(values = c("grey70", "grey40")) +
  stat_summary(fun = median, geom = "crossbar", width = 0.8, lwd = 0.2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 2), size = 5) +
  geom_text(data = multi_box, aes(x = multiMono, y = max_stock, label = letters), 
            size = 7.5, position = position_nudge(y = 7.5)) +
  scale_x_discrete(labels = c("Mono" = "Monospecific", "Multi" = "Multispecific")) +
  scale_y_continuous(breaks = c(0, 30, 60, 90), limits = c(0, 90)) +  
  xlab("") +
  ylab(expression(Delta*'Elevation (cm)')) +
  theme_classic() + 
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 18),
    axis.text.x = element_text(size = 18, color = "black")
  ) +
  Theme2  

plot_multi

setwd(plots.dir) 
dir()
name<-'1.PAPER_GREYMulticompMixedCompositionSignificanceBox' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)


#or that one 

Morfometria%>%
  left_join(letters_multi, by = "multiMono")%>%
  ggplot(aes(x = reorder(multiMono, elevation, FUN = median), y = elevation)) +
  geom_boxplot(aes(fill = multiMono, colour = multiMono),alpha = 0.6,outlier.shape = NA, col = NA)+
  geom_point(aes(fill = multiMono, colour = multiMono),size = 5, shape = 21, alpha = 0.5, 
             position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0.5))+
  theme_few() + 
  scale_color_manual(values = c("grey70", "grey40")) +  #"#9467bd", "#dbb40c"
  scale_fill_manual(values = c("grey70", "grey40")) + 
  stat_summary(fun = median, geom = "crossbar", width = 0.8, lwd = 0.2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 3), size = 5) +
  geom_text(data = multi_box, aes(x = multiMono, y = max_stock, label = letters), size = 7.5,position = position_nudge(y = 7.5)) +
  scale_x_discrete(labels=c("Mono" = "Monospecific", "Multi" = "Multispecific")) +
  xlab("") +
  ylab(expression(Delta*'Elevation (cm)')) +
  # coord_flip() +
  # facet_wrap(~habitat) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 18),
        axis.text.x = element_text(size=18))+
  Theme2

setwd(plots.dir) 
dir()
name<-'GREY-MulticompMixedCompositionSignificanceBox' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)


# Area 
Morfometria %>% 
  ggplot(aes(x = area, y = elevation)) +
  #scale_y_continuous(trans='exp', labels = scales::comma) + 
  #scale_x_continuous(trans='exp', labels = scales::comma) +
  geom_point()+
  geom_smooth()+
  theme_classic()

Morfometria %>%
  ggplot(aes(x = area, y = elevation)) +
  geom_point(color="grey40",alpha = 0.7) +
  geom_smooth(method="lm",color = "grey40", fill = "grey70", se = TRUE) +
  xlab(bquote("Vegetation patch area ("*m^2*")")) +
  ylab(expression(Delta*'Elevation (cm)')) +
  #scale_y_continuous(trans='exp', labels = scales::comma) + # Adding exponential transformation for y-axis
  #scale_x_continuous(trans='exp', labels = scales::comma) + # Adding exponential transformation for x-axis
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 18),
        axis.text.x = element_text(face = "italic", 
                                   size=18, 
                                   vjust = 1, 
                                   hjust = 1))+
  #scale_y_continuous(breaks = c(0, 30, 60, 90)) +
  theme_classic()+
  Theme2

setwd(plots.dir) 
dir()
name<-'1.PAPER_GREYAreaMixedSignificanceEXP' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)


Morfometria %>%
  ggplot(aes(x = area, y = elevation)) + 
  geom_point(aes(colour = multiMono)) +
  geom_smooth(method = "lm", aes(colour = multiMono, fill = multiMono), se = TRUE) +
 # scale_y_continuous(trans='exp', labels = scales::comma) + 
#  scale_x_continuous(trans='exp', labels = scales::comma) +
  scale_colour_manual(values=c("#3182bd", "orange")) +
  scale_fill_manual(values=c("#3182bd", "orange")) + 
  xlab("Vegetation patch area") +
  ylab("Sediment accumulation capacity") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(face = "italic", 
                                   vjust = 1, 
                                   hjust = 1))


setwd(plots.dir) 
dir()
name<-'0. AreaMixedSignificanceGroup1' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
#ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)


###########################################################################################
