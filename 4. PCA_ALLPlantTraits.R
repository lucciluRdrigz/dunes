# # # # # # PCA -  multivariate analysis: functional traits approach  # # # # # #

# # # 1st PhD chapter - Rodríguez-Arias,L.

# 12/2023

#https://rstudio-pubs-static.s3.amazonaws.com/585948_abd70a6fc3e24d4fad8944197bc5dd25.html

rm(list=ls())

# library & theme----
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(ggsci)
library(formattable)
library(vegan)
library(ggbiplot)
library(caret)
library(corrplot)

#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}

# Install ggbiplot from GitHub
#devtools::install_github("vqv/ggbiplot")


# Work directory----

work.dir=("C:/Users/lucia/Desktop/PhD_Blanes/3_R/Paper1/Paper1_mechanism_stabilization")

data.dir=paste(work.dir,"Data",sep="/")
plots.dir=paste(work.dir,"Plots",sep="/")


# Data.frame----
setwd(data.dir)
dir()
Traits <- read_csv2("Traits.csv")


Traits$spp #new column
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

Traits <- Traits %>%
  filter(!acro %in% c("ATRIP", "ELYT","PHAR","ART","SED"))


# Variables for the PCA (I keep the averages, and eliminate the variables I don't want: categorical or meaningless)
Manchas_repl1<- Traits[,-c(1:8,20,21,23,24,26,27)]


row.names(Manchas_repl1)<-Manchas_repl1$specie #names in the lateral of the table (but it is not a colunm)

## We take another column of retention data from another table and add it
#setwd(data.dir)
#dir()
#M <- read_csv2("manchasR.csv")

#Table1<- M[,c(7,14,15)] #I choose only with the three columns of interest (elevation-dominantspp-habitat)
#Table1$dominantSpecies <- sub("-PROT", "", Table1$dominantSpecies) #I remove those that are considered "protected habitats"
#Table2<- Table1[,-3]
#Table3 <- aggregate(elevation ~ dominantSpecies, data = Table2, FUN = mean) # mean of all elevations per specie
#Table4<- Table3[c(1,2,5,6,13:15),] #take the relevant species
#Now order the rows so that they are in the same order as Manchas_repl since I am going to extract the column of values and insert it into that table
#Table5<-Table4[c(2,6,3,4,1,5,7),]
#Table_final<- cbind(Manchas_repl, Table5$elevation) 

#I think I won't consider "elevation" in the analysis,b but...


# Should I scale numerical variables? 
#need to scale; drop categorical variable(spp)
#drop LA, SLA, LDW, 'cause for Salicornia & Sarcocornia - no leaf data
boxplot(scale(Manchas_repl1[,-31])) #Scale variables

# 1. PCA all traits ----
# Tabla_final: drop columns 36(spp) & la 37 (elevacion)
# Manchas_repl: drop the column 36 (spp)

# # 1. collineality: it is not a problem, in fact is good for the Principal Component Analysis
# # 2. Variance visualization 

ManchasPCA<- prcomp((Manchas_repl1[,-31]),scale. = TRUE, center = TRUE)
ManchasPCA$rotation #look at PCs
summary(ManchasPCA) #coefficient of importance
#Importance of components:
#  PC1    PC2    PC3    PC4    PC5     PC6       PC7
#Standard deviation     .9169 2.8166 2.449 1.8569 1.48414 0.94186 0.77828 0.64409 6.581e-16
#Proportion of Variance 0.2836 0.2644 0.200 0.1149 0.07342 0.02957 0.02019 0.01383 0.000e+00
#Cumulative Proportion  0.2836 0.5481 0.748 0.8630 0.93641 0.96598 0.98617 1.00000 1.000e+00


screeplot(ManchasPCA) # PC selection

rownames(ManchasPCA$x)<-Manchas_repl1$specie #para que aparezcan los nombres de las especies en el plot

# Plot final 
biplot1 <- autoplot(ManchasPCA,
         data = na.omit(Manchas_repl1),
         colour = "specie",
         x = 1, y = 2,
         label = F, frame = TRUE,
         label.vjust = -0.5,
         loadings = T,
         loadings.label.size = 2,
         loadings.label = T,
         loadings.label.hjust = 1.1,
         loadings.colour = 'black',
         loadings.label.colour = 'black',
         scale = 0,
         size = 3
) + theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    legend.text = element_text(face = "italic")
  )

biplot1
#PCs explain 54,8% of the variability


setwd(plots.dir) 
dir()
name<-'PCAfinal_allTraits' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)


# 2. PCA traits - assumptions ----
# # 1. collinearity
correl=cor(Manchas_repl1[,-31]) #analizamos las correlaciones entre predictores
print(correl)
corrplot(correl, method = "circle")

highly_correlated_vars <- findCorrelation(correl, cutoff = 0.8)
Manchas_Traits <- Manchas_repl1[, -highly_correlated_vars]

# PCA 
PCA_Traits<- prcomp(Manchas_Traits[, -c(17)],scale. = TRUE, center = TRUE)
PCA_Traits$rotation #vemos todos los Componentes Principales (PC) de cada uno
summary(PCA_Traits) #impoetancia de los componentes
screeplot(PCA_Traits) # elección de ambos componentes principales

rownames(PCA_Traits$x)<-Manchas_Traits$specie #para que aparezcan los nombres de las especies en el plot


# 2. Plot PCA traits
biplot2 <- autoplot(PCA_Traits,
                    data = na.omit(Manchas_Traits),
                    colour = "specie",
                    x = 1, y = 2,
                    label = F, frame = TRUE,
                    label.vjust = -0.5,
                    loadings = T,
                    loadings.label.size = 2,
                    loadings.label = T,
                    loadings.label.hjust = 1.1,
                    loadings.colour = 'black',
                    loadings.label.colour = 'black',
                    scale = 0,
                    size = 3
) + theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    legend.text = element_text(face = "italic")
  )

biplot2
# Around ~50.8%

setwd(plots.dir) 
dir()
name<-'PCAfinal_colineallity_out' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)

# 3. PCA traits selection ----
Manchas<- Manchas_repl1[,-c(3,7,21,24:30)] #length aerial and underground not consider
boxplot(scale(Manchas[,-21])) #ok, scale 

PCA3<- prcomp((Manchas[,-21]),scale. = TRUE, center = TRUE)
PCA3$rotation 
summary(PCA3) # more or less the same than PCA1; 5 PCs
#Importance of components:
#  PC1    PC2    PC3    PC4     PC5     PC6       PC7
#Standard deviation    2.604 2.3289 2.0472 1.37899 0.90532 0.74492 0.50627 0.26880 4.039e-16
#Proportion of Variance 0.339 0.2712 0.2095 0.09508 0.04098 0.02775 0.01282 0.00361 0.000e+00
#Cumulative Proportion  0.339 0.6102 0.8198 0.91485 0.95583 0.98357 0.99639 1.00000 1.000e+00

screeplot(PCA3) 

rownames(PCA3$x)<-Manchas$specie #para que aparezcan los nombres de las especies en el plot

# Plot final 
biplot3 <- autoplot(PCA3,
                    data = na.omit(Manchas),
                    colour = "specie",
                    x = 1, y = 2,
                    label = F, frame = TRUE,
                    label.vjust = -0.5,
                    loadings = T,
                    loadings.label.size = 2,
                    loadings.label = T,
                    loadings.label.hjust = 1.1,
                    loadings.colour = 'black',
                    loadings.label.colour = 'black',
                    scale = 0,
                    size = 3
) + theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    legend.text = element_text(face = "italic")
  )

biplot3
#PCs explain 61,02% of the variability

setwd(plots.dir) 
dir()
name<-'3. PCAfinal_selectTraits1' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)


## Plot semifinal ----
library(ggfortify)
library(ggplot2)

autoplot(
  PCA3,
  data = na.omit(Manchas),
  colour = "specie",
  x = 1, y = 2,
  label = FALSE, #Etiquetas de las especies
  frame = TRUE,
  label.vjust = -2,
  label.colour = "black",
  label.size = 4,
  loadings = TRUE,
  loadings.label.size = 4, # Tamaño del texto de las etiquetas de carga
  loadings.label = TRUE,
  loadings.label.hjust = 1.1,
  loadings.colour = 'grey40', # Color de los vectores de carga
  loadings.label.colour = 'black', # Color del texto de las etiquetas de carga
  scale = 0,
  size = 3) +
  theme_bw() +
  scale_colour_d3("category10", guide = guide_legend(override.aes = list(size = 3))) +
  coord_cartesian(xlim = c(-5,5))+
  theme(
    text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "top",
    legend.box.spacing = unit(0.2, "cm"),  # Espaciado entre elementos de la leyenda
    legend.key.size = unit(1, "lines"),   # Tamaño de los símbolos en la leyenda
    legend.text = element_text(size = 10, face = "italic")   # Tamaño del texto en la leyenda
  )  


setwd(plots.dir) 
dir()
name<-'5.4 PCA_SuuuperSelection' 
ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)

## Plot final ----

autoplot(
  PCA3,
  data = na.omit(Manchas),
  colour = "specie",
  x = 1, y = 2,
  label.vjust = -2,
  label.colour = "black",
  label.size = 7,
  loadings = TRUE, 
  loadings.label = TRUE, 
  loadings.label.hjust = 1.1,
  loadings.colour = 'grey50',
  loadings.label.colour = 'black',
  loadings.label.repel=T,
  geom_text_repel=2,
  size = 7,
  obs.scale=-5,
  var.scale=1)+
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

setwd(plots.dir) 
dir()
name<-'PCA_FinalSelection_Scaled' 
ggsave(paste(name,".png",sep="."),width = 25, height = 20,units = "cm",dpi=1600)
ggsave(paste(name,".pdf",sep="."),width = 25, height = 20,units = "cm",dpi=1600,useDingbats=FALSE)

# Tabla results PCA ----
library(flextable)

PCA3_flextable <- flextable::flextable(Manchas)
PCA3_flextable

# Scores PCA3 final to add in the GLMM ----
scoresPCA3 <- PCA3$x

score_PC1 <- PCA3$x[, 1]
score_PC1
score_PC2 <- PCA3$x[,2]
a<-cbind(Traits[,-c(1:45)],score_PC1)
b<-cbind(a,score_PC2)
rownames(score_PC1) <- a$specie
rownames(score_PC2) <- a$specie

PCAscores <- b
write.csv2(PCAscores, file = "C:/Users/lucia/Desktop/PhD_Blanes/3_R/Paper1/Paper1_mechanism_stabilization/Data/PCAscores.csv", row.names = FALSE)


# # Extrayendo data  ----
#https://rstudio-pubs-static.s3.amazonaws.com/585948_abd70a6fc3e24d4fad8944197bc5dd25.html
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
library(ade4) 
pcaADE <- dudi.pca(Manchas[,-21], scannf = F, nf = 8)
#The $co is the coordinates of variables in PCA space. Equivalent to loadings*sdev as calculated in theory section above for prcomp.
#The $li correspond to the individual, or row, cooridinates 
scatter(pcaADE)
library(grid) # has the viewport function, needed for insetting the scree plot
library(ggpubr) #for making publication ready plots
library(ggforce) #for some nice polygons geoms
library(ggalt) #contains some extra geoms
library(hrbrthemes) #some nice themes for ggplot
library(factoextra)

fviz_eig(PCA3)
# Eigenvalues
eig.val <- get_eigenvalue(PCA3)
eig.val
#        eigenvalue variance.percent cumulative.variance.percent
#Dim.1 6.780314e+00     3.390157e+01                    33.90157
#Dim.2 5.423934e+00     2.711967e+01                    61.02124
#Dim.3 4.191068e+00     2.095534e+01                    81.97658
#Dim.4 1.901616e+00     9.508081e+00                    91.48466
#Dim.5 8.196023e-01     4.098011e+00                    95.58267
#Dim.6 5.549028e-01     2.774514e+00                    98.35719
#Dim.7 2.563111e-01     1.281555e+00                    99.63874


res.var <- get_pca_var(PCA3)
res.var
#res.var$coord          # Coordinates
#res.var$contrib        # Contributions to the PCs
#res.var$cos2           # Quality of representation 
# Results for individuals (i.e rows)
res.ind <- get_pca_ind(PCA3)
res.ind

# Con cuantos PC me quedo?? ----
#visualmente
t.pca <- prcomp(Manchas[,-21], center = TRUE, scale. = TRUE)
screeplot(t.pca, type="lines") #retener 4 PC segun caída/codo de la pendiente
#criterio Kaiser
(t.pca$sdev)^2
#Criterio de Kaiser: sólo debemos retener los componentes principales cuya varianza sea superior a 1 (cuando se aplicó el análisis de componentes principales a datos estandarizados)
# mantengo 4 componentes ; mantener el número de componentes necesarios para explicar al menos una cantidad mínima de la varianza total.
summary(t.pca) #con 4 componentes, 91.48 varianza explicada
#Importance of components:
#PC1    PC2    PC3     PC4     PC5     PC6     PC7     PC8       PC9
#Standard deviation     2.604 2.3289 2.0472 1.37899 0.90532 0.74492 0.50627 0.26880 4.039e-16
#Proportion of Variance 0.339 0.2712 0.2095 0.09508 0.04098 0.02775 0.01282 0.00361 0.000e+00
#Cumulative Proportion  0.339 0.6102 0.8198 0.91485 0.95583 0.98357 0.99639 1.00000 1.000e+00

t.pca <- prcomp(Manchas[,-21], center = TRUE, scale. = TRUE)
site.groups <- c(rep("a", 3), rep("b", 3), rep("c", 3))
ggbiplot(t.pca, labels=rownames(Manchas), groups = site.groups, ellipse = TRUE)


fviz_pca_var(PCA3,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+theme_classic()

fviz_pca_biplot(PCA3, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)+theme_classic()

fviz_pca_biplot(PCA3, 
                col.ind = Manchas$specie, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species") 


fviz_pca_biplot(PCA3, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = Manchas$specie,
                col.ind = "black",
                # Color variable by groups
              #  col.var = factor(c("sepal", "sepal", "petal", "petal")),
                
                legend.title = list(fill = "specie", color = "Clusters"),
                repel = TRUE        # Avoid label overplotting
)+
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("npg")  +
  theme_classic()# Variable colors

fviz_pca_biplot(PCA3, 
                # Individuals
                geom.ind = "point",
                fill.ind = Manchas$specie, col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdYlBu",
                
                legend.title = list(fill = "specie", color = "Contrib",
                                    alpha = "Contrib")
)+
  theme_classic()
