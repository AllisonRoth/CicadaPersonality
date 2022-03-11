#Update R packages
update.packages(ask = FALSE, checkBuilt = TRUE)

# the pacakge needed
library(here) # making the path works for all of people running it from R
library(tidyverse)
library(broom)
library(broom.mixed)
library(lme4)
library(ggplot2) 
library(lmeresampler)
library(dplyr)
library(MuMIn)
library(lmerTest)
library(glmmTMB)
library(ggpubr)
library(emmeans)


#######Novel Arena#######
#Only look at MALES for NA
m<-read.csv(here("data/CicadaRawNAforPCA.csv"),colClasses="character")
head(m)
#change number columns of csv load to numeric
m$Leave<-as.numeric(m$Leave)
m$Quadrats <-as.numeric(m$Quadrats )
m$Trans<-as.numeric(m$Trans)
m$Sides<-as.numeric(m$Sides)
m$ID<-as.factor(m$ID)

#Find out how many unique individuals were tested
length(unique(m[["ID"]]))

#Standardize (Z-transform) variables
m$LeaveNew <- (m$Leave - mean(m$Leave)) / sd(m$Leave)
m$Quadrats <- (m$Quadrats - mean(m$Quadrats)) / sd(m$Quadrats)
m$Trans <- (m$Trans - mean(m$Trans)) / sd(m$Trans)
m$Sides <- (m$Sides - mean(m$Sides)) / sd(m$Sides)

##Run PCA
#Subset dataframe to only keep columns we want
SubAll<-data.frame(m$Leave, m$Quadrats, m$Trans, m$Sides)
head(SubAll)
#PCA
pcAll<-princomp(SubAll,cor=TRUE,scores=TRUE)
summary(pcAll)

biplot(pcAll)
pcAll$loadings
AllScores<-pcAll$scores
AllScores

#See loadings
PCApr <- prcomp(SubAll,center = TRUE,scale. = TRUE) 
print(PCApr)
plot(PCApr, type = "l")
summary(PCApr)
biplot(PCApr)
#write.csv(AllScores, file = "CicadaNAPC1.csv")


##HISTOGRAM and qqplots for PC1 WITHOUT TRANSFORMATION
#Load PC1 data
n<-read.csv(here("data/CicadaNAPC1.csv"),colClasses="character")
head(n)
n$PC1<-as.numeric(n$PC1)
#Switch sign in PC1 so slow explorers have a smaller PC1 and fast explorers have a larger PC
n$PC1<-n$PC1*-1
hist(n$PC1,xlab="PC1",breaks=20,main="")
qqnorm(n$PC1,main="")
qqline(n$PC1)
#Data is NOT normally distributed!!!!So use MODIFIED VERSION OF Nakagawa's code that I used in Roth et al. 2021 JAE 

#Residual Plot for PC1 WITHOUT TRANSFORMATION
#library(lme4)
Naka<-lmer(PC1~(1|Assay)+(1|DaysSinceMay15)+(1|Time)+(1|Temp)+(1|Humidity)+(1|DaysSinceEmerge)+(1|ID),data=n)
summary(Naka)
#library(ggplot2) 
head(fortify(Naka))
residPlot<-ggplot(aes(x=.fitted,y=.resid),data=Naka)+geom_point()+geom_hline(yintercept=0)+labs(x="Fitted Values",y="Residuals")
residPlot

###########Test for REPEATABILITY for NA Assay
###############REPEATABILITIES WITH NAKAGAWA'S *MODIEFIED* CODE###############
#This is the final and finished method!!!!

# remove.packages("lmeresampler")
# install.packages("lmeresampler")
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("aloy/lmeresampler")
#library(lmeresampler)

#---
#title: "Figuring out how to get CI with non-parametric boostraps with linear mixed models"
#author: "Shinichi Nakagawa"
#date: "09/08/2018"
#output: html_document
#---

## Background

#We will need to get confidence intervals (CIs) for repeatability for Tweedie distributed data usimng the `cplm` package. We will use chicken behaviour data wiht adjusted PC scores. Actually, it is easy to get "credible" intervals using the `bcplm` function so we will just do that. However, it turns out moding with `cplm` does not really work (see below). So we will use non-parametric bootstrapping and normal linear models (`lme4::lmer`). In this way, we can violate linear model assumptions but we can still get meaningful CIs from `lmer`. 

## Setting up

#```{r setup, echo=FALSE}
# knitr::opts_chunk$set(
#   message = FALSE,
#   warning = FALSE,
#   tidy = TRUE,
#   echo = TRUE,
#   fig.width = 8
# )
# # clearning up
# rm(list=ls())

# loading packages
#Install pacman
pacman::p_load(tidyverse, purrr, cplm, lme4, lmeresampler, boot)

# getting functions
source(here("R/function2.R"), chdir = TRUE)

#```


## Data preparation

#```{r prep}

#The dataset for the cicadas (males only is "CicadaNAPC1" with the file loaded into R as:

n<-read.csv(here("data/CicadaNAPC1.csv"),colClasses="character")

#Change PC1, Assay, Time, Temp, Humidity, and DaysSinceEmerge to numeric
n$PC1<-as.numeric(n$PC1)
n$Assay<-as.numeric(n$Assay)
n$DaysSinceMay15<-as.numeric(n$DaysSinceMay15)
n$Time<-as.numeric(n$Time)
n$Temp<-as.numeric(n$Temp)
n$Humidity<-as.numeric(n$Humidity)
n$DaysSinceEmerge<-as.numeric(n$DaysSinceEmerge)
head(n)

#Switch sign in PC1 so slow explorers have a smaller PC1 and fast explorers have a larger PC
n$PC1<-n$PC1*-1

### Modeling with lmer

#```{r lmer}
lmer1 <- lmer(PC1 ~ Assay + DaysSinceMay15 + Time + Temp + Humidity + DaysSinceEmerge + (1 | ID), data = n)
lmer100 <- lmer(PC1*100 ~ Assay + DaysSinceMay15 + Time + Temp + Humidity + DaysSinceEmerge + (1 | ID), data = n)

# unlike cpglmm models, t values are invariant
# this is what we want
summary(lmer1)
tidy(lmer1)
summary(lmer100)
tidy(lmer100)

#```


### Getting correct CI by boostrapping (including repeatability)
#Reinstall dplyr (Shinichi's code wasn't working because I was using an old version of dplyr) 
# install.packages("dplyr")
# library(dplyr)

#We will use fully non-boostrap models so that we can violate almost all assumptions of linear models including the normality of residuals

#```{r boot}
# the function to get non-paramertic boostrap confidence intervals for all parameters including repeatablity (rep)
# the number of boostrap should probably be 1000 - 10000 for your publications
# Note we use the lmeresampler package and the boot package
# Warning - this takes a while!
new_res <- nonpara_boot2(model = lmer1, bootN = 1000)

# sigma2 = within-individual variance (residual); tau2 = between-indvidual variance (Individual)
# Note that percentile CI would look very asymmetric for variance components and related values like repetablity
new_res 

#```

## Information for this R session
#```{r session}
#sessionInfo()
#```  

#```{r session}
#sessionInfo()
#``` 



###############Get a single EB Score for each individual###############
n<-read.csv(here("data/CicadaNAPC1.csv"),colClasses="character")
n$PC1<-as.numeric(n$PC1)
n$Assay<-as.numeric(n$Assay)
n$DaysSinceMay15<-as.numeric(n$DaysSinceMay15)
n$Time<-as.numeric(n$Time)
n$Temp<-as.numeric(n$Temp)
n$Humidity<-as.numeric(n$Humidity)
n$DaysSinceEmerge<-as.numeric(n$DaysSinceEmerge)
#Multiply PC1 by -1 so that faster explorers have higher scores.
n$PC1n <- n$PC1*-1

#In order to control for fixed effects we used a LM of PC1 to obtain a single EB score for each individual. 
#Use: Assay, DaysSinceMay15, Time, Temp, Humidity, DaysSinceEmerge, and Individual identity as fixed effects
#USE BELOW METHOD FOR MODEL SELECTION 
install.packages("MuMIn")
#library(MuMIn)
require(MuMIn)
globalmodel <- lm(PC1n~ID+Assay+DaysSinceMay15+Time+Temp+Humidity+DaysSinceEmerge,data=n,na.action = "na.fail")
combinations <- dredge(globalmodel)
print(combinations)
#Now from this printed list find the model with the lowest AICc that also INCLUDES Individual ID (since ID always needs to be kept in the model)
NA1<-lm(PC1n~ID - 1 +Time+Temp,data=n)
summary(NA1)


#Extract the coefficients data frame so we can write the output into a CSV file
results_dfm <-summary.lm(NA1)$coefficients
results_dfm
#write.csv(results_dfm, file = "LmEBCicadas.csv")
#Output shows faster explores as having higher values 














#######Tonic Immobility: MALES#######
tm<-read.csv(here("data/CicadaRawTIforPCAm.csv"),colClasses="character")
head(tm)
#change number columns of csv load to numeric
tm$LatencyToFreeze<-as.numeric(tm$LatencyToFreeze)
tm$TimeFrozen <-as.numeric(tm$TimeFrozen)
tm$TimesFlipped<-as.numeric(tm$TimesFlipped)
tm$ID<-as.factor(tm$ID)

#Find out how many unique individuals were tested
length(unique(tm[["ID"]]))

#Standardize (Z-transform) variables
tm$LatencyToFreeze <- (tm$LatencyToFreeze - mean(tm$LatencyToFreeze)) / sd(tm$LatencyToFreeze)
tm$TimeFrozen <- (tm$TimeFrozen - mean(tm$TimeFrozen)) / sd(tm$TimeFrozen)
tm$TimesFlipped <- (tm$TimesFlipped - mean(tm$TimesFlipped)) / sd(tm$TimesFlipped)

##Run PCA
#Subset dataframe to only keep columns we want
SubAlltm<-data.frame(tm$LatencyToFreeze, tm$TimeFrozen, tm$TimesFlipped)
head(SubAlltm)
#PCA
pcAlltm<-princomp(SubAlltm,cor=TRUE,scores=TRUE)
summary(pcAlltm)

biplot(pcAlltm)
pcAlltm$loadings
AllScorestm<-pcAlltm$scores
AllScorestm

#See loadings
PCAprtm <- prcomp(SubAlltm,center = TRUE,scale. = TRUE) 
print(PCAprtm)
plot(PCAprtm, type = "l")
summary(PCAprtm)
biplot(PCAprtm)

#write.csv(AllScorestm, file = "CicadaTIPC1m.csv")


##HISTOGRAM and qqplots for PC1 WITHOUT TRANSFORMATION
#Load PC1 data
ntm<-read.csv("data/CicadaTIPC1m.csv",colClasses="character")
head(ntm)
ntm$PC1<-as.numeric(ntm$PC1)
#Bolder Individuals Already Have Higher PC1
hist(ntm$PC1,xlab="PC1",breaks=20,main="")
qqnorm(ntm$PC1,main="")
qqline(ntm$PC1)
#Data is NOT normally distributed!!!!So use A MODIFIED VERSION OF Nakagawa's code that I used in Roth et al. 2021 JAE 




###########Test for REPEATABILITY for TI ASSAY MALES
###############REPEATABILITIES WITH NAKAGAWA'S *MODIFIED* CODE###############
#This is the final and finished method!!!!

# remove.packages("lmeresampler")
# install.packages("lmeresampler")
# library(lmeresampler)

#---
#title: "Figuring out how to get CI with non-parametric boostraps with linear mixed models"
#author: "Shinichi Nakagawa"
#date: "09/08/2018"
#output: html_document
#---

## Background

#We will need to get confidence intervals (CIs) for repeatability for Tweedie distributed data usimng the `cplm` package. We will use chicken behaviour data wiht adjusted PC scores. Actually, it is easy to get "credible" intervals using the `bcplm` function so we will just do that. However, it turns out moding with `cplm` does not really work (see below). So we will use non-parametric bootstrapping and normal linear models (`lme4::lmer`). In this way, we can violate linear model assumptions but we can still get meaningful CIs from `lmer`. 


ntm<-read.csv("data/CicadaTIPC1m.csv",colClasses="character")

#Change PC1, Assay, Time, Temp, Humidity, and DaysSinceEmerge to numeric
ntm$PC1<-as.numeric(ntm$PC1)
ntm$PC2<-as.numeric(ntm$PC2)
ntm$Assay<-as.numeric(ntm$Assay)
ntm$DaysSinceMay15<-as.numeric(ntm$DaysSinceMay15)
ntm$Time<-as.numeric(ntm$Time)
ntm$Temp<-as.numeric(ntm$Temp)
ntm$Humidity<-as.numeric(ntm$Humidity)
ntm$DaysSinceEmerge<-as.numeric(ntm$DaysSinceEmerge)
head(ntm)

#Bolder Individuals Already Have Higher PC1

### Modeling with lmer

#```{r lmer}
lmer1tm <- lmer(PC1 ~ Assay + DaysSinceMay15 + Time + Temp + Humidity + Observer + DaysSinceEmerge + (1 | ID), data = ntm)
lmer100tm <- lmer(PC1*100 ~ Assay + DaysSinceMay15 + Time + Temp + Humidity + Observer + DaysSinceEmerge + (1 | ID), data = ntm)

# unlike cpglmm models, t values are invariant
# this is what we want
summary(lmer1tm)
summary(lmer100tm)


#We will use fully non-boostrap models so that we can violate almost all assumptions of linear models including the normality of residuals

#```{r boot}
# the function to get non-paramertic boostrap confidence intervals for all parameters including repeatablity (rep)
# the number of boostrap should probably be 1000 - 10000 for your publications
# Note we use the lmeresampler package and the boot package
# Warning - this takes a while!
new_restm <- nonpara_boot2(model = lmer1tm, bootN = 1000)

# sigma2 = within-individual variance (residual); tau2 = between-indvidual variance (Individual)
# Note that percentile CI would look very asymmetric for variance components and related values like repetablity
new_restm 

#```

## Information for this R session
#```{r session}
#sessionInfo()
#```  

#```{r session}
#sessionInfo()
#``` 



###############Get a single TI Score for each male###############
ntm<-read.csv("data/CicadaTIPC1m.csv",colClasses="character")
ntm$PC1<-as.numeric(ntm$PC1)
ntm$Assay<-as.numeric(ntm$Assay)
ntm$DaysSinceMay15<-as.numeric(ntm$DaysSinceMay15)
ntm$Time<-as.numeric(ntm$Time)
ntm$Temp<-as.numeric(ntm$Temp)
ntm$Humidity<-as.numeric(ntm$Humidity)
ntm$DaysSinceEmerge<-as.numeric(ntm$DaysSinceEmerge)

#Bolder individuals (i.e. those that spent less time in TI, took longer to enter TI, or had to be righted more often) alread have higher scores.

#In order to control for fixed effects we used a LM of PC1 to obtain a single EB score for each individual. 

#Use: Assay, DaysSinceMay15, Time, Temp, Humidity, DaysSinceEmerge, Observer, and Individual identity as fixed effects

#USE BELOW METHOD FOR MODEL SELECTION 
#library(MuMIn)
require(MuMIn)
globalmodel <- lm(PC1~ID+Assay+DaysSinceMay15+Time+Temp+Humidity+DaysSinceEmerge+Observer,data=ntm,na.action = "na.fail")
combinations <- dredge(globalmodel)
print(combinations)

#Now from this printed list find the model with the lowest AICc that also INCLUDES Individual ID (since ID always needs to be kept in the model)
TI1m<-lm(PC1~ID - 1 +Temp,data=ntm)
summary(TI1m)


#Extract the coefficients data frame so we can write the output into a CSV file
results_dfm <-summary.lm(TI1m)$coefficients
results_dfm
#write.csv(results_dfm, file = "LmTICicadasM.csv")
#Output shows Bolder Individuals as having higher values 














#######Tonic Immobility: FEMALES#######
tf<-read.csv("data/CicadaRawTIforPCAf.csv",colClasses="character")
head(tf)
#change number columns of csv load to numeric
tf$LatencyToFreeze<-as.numeric(tf$LatencyToFreeze)
tf$TimeFrozen <-as.numeric(tf$TimeFrozen)
tf$TimesFlipped<-as.numeric(tf$TimesFlipped)
tf$ID<-as.factor(tf$ID)

#Find out how many unique individuals were tested
length(unique(tf[["ID"]]))

#Standardize (Z-transform) variables
tf$LatencyToFreeze <- (tf$LatencyToFreeze - mean(tf$LatencyToFreeze)) / sd(tf$LatencyToFreeze)
tf$TimeFrozen <- (tf$TimeFrozen - mean(tf$TimeFrozen)) / sd(tf$TimeFrozen)
tf$TimesFlipped <- (tf$TimesFlipped - mean(tf$TimesFlipped)) / sd(tf$TimesFlipped)

##Run PCA
#Subset dataframe to only keep columns we want
SubAlltf<-data.frame(tf$LatencyToFreeze, tf$TimeFrozen, tf$TimesFlipped)
head(SubAlltf)
#PCA
pcAlltf<-princomp(SubAlltf,cor=TRUE,scores=TRUE)
summary(pcAlltf)

biplot(pcAlltf)
pcAlltf$loadings
AllScorestf<-pcAlltf$scores
AllScorestf

#See loadings
PCAprtf <- prcomp(SubAlltf,center = TRUE,scale. = TRUE) 
print(PCAprtf)
plot(PCAprtf, type = "l")
summary(PCAprtf)
biplot(PCAprtf)

#write.csv(AllScorestf, file = "CicadaTIPC1f.csv")


##HISTOGRAM and qqplots for PC1 WITHOUT TRANSFORMATION
#Load PC1 data
ntf<-read.csv("C:/Users/Owner/Desktop/CicadaTIPC1f.csv",colClasses="character")
head(ntf)
ntf$PC1<-as.numeric(ntf$PC1)
#Bolder Individuals Already Have Higher PC1
hist(ntf$PC1,xlab="PC1",breaks=20,main="")
qqnorm(ntf$PC1,main="")
qqline(ntf$PC1)
#Data is NOT normally distributed!!!!So use A MODIFIED VERSION OF Nakagawa's code that I used in Roth et al. 2021 JAE 




###########Test for REPEATABILITY for TI ASSAY FEMALES
###############REPEATABILITIES WITH NAKAGAWA'S *MODIEFIED* CODE###############


#The dataset for the cicadas (FEMALEs only is "CicadaTIPC1f" with the file loaded into R as:

ntf<-read.csv("data/CicadaTIPC1f.csv",colClasses="character")

#Change PC1, Assay, Time, Temp, Humidity, and DaysSinceEmerge to numeric
ntf$PC1<-as.numeric(ntf$PC1)
ntf$PC2<-as.numeric(ntf$PC2)
ntf$Assay<-as.numeric(ntf$Assay)
ntf$DaysSinceMay15<-as.numeric(ntf$DaysSinceMay15)
ntf$Time<-as.numeric(ntf$Time)
ntf$Temp<-as.numeric(ntf$Temp)
ntf$Humidity<-as.numeric(ntf$Humidity)
ntf$DaysSinceEmerge<-as.numeric(ntf$DaysSinceEmerge)
head(ntf)

#Bolder Individuals Already Have Higher PC1

### Modeling with lmer

#```{r lmer}
lmer1tf <- lmer(PC1 ~ Assay + DaysSinceMay15 + Time + Temp + Humidity + Observer + DaysSinceEmerge + (1 | ID), data = ntf)
lmer100tf <- lmer(PC1*100 ~ Assay + DaysSinceMay15 + Time + Temp + Humidity + Observer + DaysSinceEmerge + (1 | ID), data = ntf)

# unlike cpglmm models, t values are invariant
# this is what we want
summary(lmer1tf)
summary(lmer100tf)

#```


### Getting correct CI by boostrapping (including repeatability)

#We will use fully non-boostrap models so that we can violate almost all assumptions of linear models including the normality of residuals

#```{r boot}
# the function to get non-paramertic boostrap confidence intervals for all parameters including repeatablity (rep)
# the number of boostrap should probably be 1000 - 10000 for your publications
# Note we use the lmeresampler package and the boot package
# Warning - this takes a while!
new_restf <- nonpara_boot2(model = lmer1tf, bootN = 1000)

# sigma2 = within-individual variance (residual); tau2 = between-indvidual variance (Individual)
# Note that percentile CI would look very asymmetric for variance components and related values like repetablity
new_restf 

#```

## Information for this R session
#```{r session}
#sessionInfo()
#```  

#```{r session}
#sessionInfo()
#```

###############Get a single TI Score for each female###############
ntf<-read.csv("data/CicadaTIPC1f.csv",colClasses="character")
ntf$PC1<-as.numeric(ntf$PC1)
ntf$Assay<-as.numeric(ntf$Assay)
ntf$DaysSinceMay15<-as.numeric(ntf$DaysSinceMay15)
ntf$Time<-as.numeric(ntf$Time)
ntf$Temp<-as.numeric(ntf$Temp)
ntf$Humidity<-as.numeric(ntf$Humidity)
ntf$DaysSinceEmerge<-as.numeric(ntf$DaysSinceEmerge)

#Bolder individuals (i.e. those that spent less time in TI, took longer to enter TI, or had to be righted more often) alread have higher scores.

#In order to control for fixed effects we used a LM of PC1 to obtain a single EB score for each individual. 

#Use: Assay, DaysSinceMay15, Time, Temp, Humidity, DaysSinceEmerge, Observer, and Individual identity as fixed effects

#USE BELOW METHOD FOR MODEL SELECTION 
# library(MuMIn)
# require(MuMIn)
globalmodel <- lm(PC1~ID+Assay+DaysSinceMay15+Time+Temp+Humidity+DaysSinceEmerge+Observer,data=ntf,na.action = "na.fail")
combinations <- dredge(globalmodel)
print(combinations)

#Now from this printed list find the model with the lowest AICc that also INCLUDES Individual ID (since ID always needs to be kept in the model)
TI1f<-lm(PC1~ID - 1 +Humidity,data=ntf)
summary(TI1f)


#Extract the coefficients data frame so we can write the output into a CSV file
results_dfm <-summary.lm(TI1f)$coefficients
results_dfm
#write.csv(results_dfm, file = "LmTICicadasF.csv")
#Output shows Bolder Individuals as having higher values



########Are EB and TI correlated?
BS<-read.csv("data/CicadaBehavSynd.csv",colClasses="character")
BS$TI<-as.numeric(BS$TI)
BS$EB<-as.numeric(BS$EB)
cor.test(BS$EB,BS$TI)
#No significant correlation between EB and TI


#########MATING ANALYSES#########
#library(glmmTMB)
I<-read.csv("data/CicadaCollatedMatingData.csv",colClasses="character")
head(I)
I$Pen[I$Pen == "1"] <- "A"
I$Pen[I$Pen == "2"] <- "B"
I$Pen[I$Pen == "3"] <- "C"
I$Pen[I$Pen == "4"] <- "D"

I$ID<-as.factor(I$ID)
I$Group<-as.factor(I$Group)
I$Died<-as.numeric(I$Died)
I$Mated<-as.numeric(I$Mated)
I$DurationCop<-as.numeric(I$DurationCop)
I$MM<-as.numeric(I$MM)
I$MU<-as.numeric(I$MU)
I$Weight<-as.numeric(I$Weight)
I$EB<-as.numeric(I$EB)
I$TI<-as.numeric(I$TI)


####Check to see EB distributions of males that died during the mating trials compared to those that did NOT die 
alive<-I[I$DiedDuring=="N",]
dead<-I[I$DiedDuring=="Y",]

range(alive$EB, na.rm = TRUE)
range(dead$EB, na.rm = TRUE)

hist(alive$EB, breaks=5,xlab="Exploration Behavior")
hist(dead$EB, breaks=5,xlab="Exploration Behavior")

shapiro.test(alive$EB)
shapiro.test(dead$EB)
#library(ggpubr)
t.test(alive$EB, dead$EB, alternative = "two.sided", var.equal = FALSE)
#NO SIGNIFICANT DIFFERENCE IN EB BETWEEN DEAD AND ALIVE MALES - so it's ok to leave them in the analyses



####Check to see TI distributions of males that died during the mating trials compared to those that did NOT die 
range(alive$TI, na.rm = TRUE)
range(dead$TI, na.rm = TRUE)

hist(alive$TI, breaks=5,xlab="Tonic Immobility")
hist(dead$TI, breaks=5,xlab="Tonic Immobility")

shapiro.test(alive$TI)
shapiro.test(dead$TI)
if(!require(psych)){install.packages("psych")}
if(!require(FSA)){install.packages("FSA")}
if(!require(lattice)){install.packages("lattice")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(coin)){install.packages("coin")}
if(!require(DescTools)){install.packages("DescTools")}
if(!require(effsize)){install.packages("effsize")}
wilcox.test(TI ~ DiedDuring, data=I)
#NO SIGNIFICANT DIFFERENCE IN TI BETWEEN DEAD AND ALIVE MALES - so it's ok to leave them in the analyses




####Check to see Weight distributions of males that died during the mating trials compared to those that did NOT die 
range(alive$Weight, na.rm = TRUE)
range(dead$Weight, na.rm = TRUE)

hist(alive$Weight, breaks=5,xlab="Tonic Immobility")
hist(dead$Weight, breaks=5,xlab="Tonic Immobility")

shapiro.test(alive$Weight)
shapiro.test(dead$Weight)
if(!require(psych)){install.packages("psych")}
if(!require(FSA)){install.packages("FSA")}
if(!require(lattice)){install.packages("lattice")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(coin)){install.packages("coin")}
if(!require(DescTools)){install.packages("DescTools")}
if(!require(effsize)){install.packages("effsize")}
wilcox.test(Weight ~ DiedDuring, data=I)
#NO SIGNIFICANT DIFFERENCE IN WEIGHT BETWEEN DEAD AND ALIVE MALES - so it's ok to leave them in the analyses




#SCALE VARIABLES 
#https://stackoverflow.com/questions/41766181/correct-way-to-scale-for-multilevel-regression-using-lmer-r
I <- transform(I,Weight.sc=scale(Weight),EB.sc=scale(EB),TI.sc=scale(TI))







###########SUCCESSFUL COPULATIONS###########
Test1<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Location + Date +  (1|Pen) + (1|ID) + (1|Group), data=I, family = binomial)
summary(Test1) 
tidy(Test1)

# just assuming this is the best mode
#library(emmeans)

pred_results <- Test1 %>% 
  emmeans(~ EB.sc,
          #by = I$EB.sc,
          at = list(EB.sc = seq(min(I$EB.sc),max(I$EB.sc),length.out = 100)),
          type="response",
          re_formula = NA) %>% as.data.frame()

ggplot() +
  # putting bubbles
  geom_point(data = I, aes(x = EB.sc, y = Mated, color = as.factor(Mated)), size = 4, alpha = 0.8, fill = "grey90" ) +
  # confidence interval
  geom_smooth(data = pred_results, aes(x = EB.sc, y = lower.CL), method =  "loess", formula = y~x, se = FALSE,lty = "dotted", col = "black") +
  geom_smooth(data = pred_results, aes(x = EB.sc, y = upper.CL), method =  "loess", formula = y~x, se = FALSE, lty ="dotted", col = "black") +
  # main line
  geom_smooth(data = pred_results, aes(x = EB.sc, y = prob), method =  "loess", formula = y~x, se = FALSE, col = "black") +
  theme_bw() +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.position = "none")  +
  labs(y = "Successful Copulation (Yes = 1 & No = 0)", x = "Exploration Score (z-transformed)")







########### M-M ATTEMPTED COPULATIONS ###########
hist(I$MM, breaks=20,xlab="M-M Attempted Copulations")
Test2<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Location + Date + (1|Pen) + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
summary(Test2) 
tidy(Test2)


pred_results <- Test2 %>% 
  emmeans(~ EB.sc,
          #by = I$EB.sc,
          at = list(EB.sc = seq(min(I$EB.sc),max(I$EB.sc),length.out = 100)),
          type="response",
          re_formula = NA) %>% as.data.frame()

ggplot() +
  # putting bubbles
  geom_point(data = I, aes(x = EB.sc, y = MM, color = MM), size = 4, alpha = 0.6, fill = "grey90" ) +
  # confidence interval
  geom_smooth(data = pred_results, aes(x = EB.sc, y = lower.CL), method =  "loess", formula = y~x, se = FALSE,lty = "dotted", col = "black") +
  geom_smooth(data = pred_results, aes(x = EB.sc, y = upper.CL), method =  "loess", formula = y~x, se = FALSE, lty ="dotted", col = "black") +
  # main line
  geom_smooth(data = pred_results, aes(x = EB.sc, y = rate), method =  "loess", formula = y~x, se = FALSE, col = "black") +
  theme_bw() +
  scale_color_gradient() +
  theme(legend.position = "none")  +
  labs(y = "Number of Male-Male Copulation Attempts", x = "Exploration Score (z-transformed)")







###########ALL (BOTH SUCCESSFUL AND UNSUCCESSFUL) ATTEMPTED COPULATIONS WITH EITHER SEX###########
I$All<-I$Mated+I$MU+I$MM
Test0<-glmmTMB(All ~ EB.sc + TI.sc + Weight.sc + Location + Date + (1|Pen) + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
summary(Test0)
tidy(Test0)



pred_results <- Test0 %>% 
  emmeans(~ EB.sc,
          #by = I$EB.sc,
          at = list(EB.sc = seq(min(I$EB.sc),max(I$EB.sc),length.out = 100)),
          type="response",
          re_formula = NA) %>% as.data.frame()

ggplot() +
  # putting bubbles
  geom_point(data = I, aes(x = EB.sc, y = All, color = All), size = 4, alpha = 0.6, fill = "grey90" ) +
  # confidence interval
  geom_smooth(data = pred_results, aes(x = EB.sc, y = lower.CL), method =  "loess", formula = y~x, se = FALSE,lty = "dotted", col = "black") +
  geom_smooth(data = pred_results, aes(x = EB.sc, y = upper.CL), method =  "loess", formula = y~x, se = FALSE, lty ="dotted", col = "black") +
  # main line
  geom_smooth(data = pred_results, aes(x = EB.sc, y = rate), method =  "loess", formula = y~x, se = FALSE, col = "black") +
  ylim(c(0,6))+
  theme_bw() +
  scale_color_gradient() +
  theme(legend.position = "none")  +
  labs(y = "Number of All Copulation Attempts", x = "Exploration Score (z-transformed)")







###########M-M ADVERTISEMENT RATE###########
A<-read.csv("data/CicadaCollatedAdvertData.csv",colClasses="character")
head(A)
A$ID<-as.factor(A$ID)
A$DurationCI<-as.numeric(A$DurationCI)
A$NumberAdvert<-as.numeric(A$NumberAdvert)
A$Weight<-as.numeric(A$Weight)
A$EB<-as.numeric(A$EB)
A$TI<-as.numeric(A$TI)
#SCALE VARIABLES 
A <- transform(A,Weight.sc=scale(Weight),EB.sc=scale(EB),TI.sc=scale(TI))

hist(A$DurationCI)
shapiro.test(A$DurationCI)
hist(A$NumberAdvert)

Test3 <- lm(DurationCI ~ EB.sc, data = A)
summary(Test3)










####Summary stats for MS
#Median length of successful male-female copulation attempts
mated<-I[I$Mated==1,]
median(mated$DurationCop)

#Range and mean of unsuccessful male-female copulations per male
range(I$MU)
mean(I$MU)

#Range and mean of male-male copulation atempts per male
range(I$MM)
mean(I$MM)

#Range and mean duration of singing
range(A$DurationCI)
mean(A$DurationCI)

#Range and mean number of calling song bouts
range(A$NumberAdvert)
mean(A$NumberAdvert)



