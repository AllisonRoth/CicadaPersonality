#Update R packages
update.packages(ask = FALSE, checkBuilt = TRUE)
#######Novel Arena#######
#Only look at MALES for NA!
m<-read.csv("C:/Users/Owner/Desktop/CicadaRawNAforPCA.csv",colClasses="character")
head(m)
#change number columns of csv load to numeric
m$Leave<-as.numeric(m$Leave)
m$Quadrats <-as.numeric(m$Quadrats )
m$Trans<-as.numeric(m$Trans)
m$Sides<-as.numeric(m$Sides)
m$ID<-as.factor(m$ID)

#Find out how many unique individuals were tested
length(unique(m[["ID"]]))

#Change Leave so it actually represents 5 minutes (300s) minus leave time so that it lines up in directionality with the other variables. e.g., A larger variable for this new variable which I'll call "LeaveNew" (i.e., time spent moving) will represent a faster exploring individual. Faster exploring individuals will also have higher values for Quadrats, Trans, and Sides.
#m$LeaveNew<-300-(m$Leave)
#head(m)
#NOTE: Looks like this doesn't actually matter for PCA (so didn't actually need to do this!)

#Standardize (Z-transform) variables
#https://stats.stackexchange.com/questions/69157/why-do-we-need-to-normalize-data-before-principal-component-analysis-pca  
#https://stackoverflow.com/questions/15215457/standardize-data-columns-in-r
#The most common normalization is the z-transformation, where you subtract the mean and divide by the standard deviation of your variable. The result will have mean=0 and sd=1.
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
#0.836176 of variance is explained bt PC1 

biplot(pcAll)
pcAll$loadings
AllScores<-pcAll$scores
AllScores

#See loadings
PCApr <- prcomp(SubAll,center = TRUE,scale. = TRUE) 
####IMPORTANT Different variables are loaded strongly and consistently
print(PCApr)
plot(PCApr, type = "l")
summary(PCApr)
biplot(PCApr)

#write.csv(AllScores, file = "CicadaNAPC1.csv")
#Have to add Individual IDs in to the CSV file (Copy and paste "ID" column from "CicadaRawNAforPCA.csv" into "CicadaNAPC1.csv"). Remove everything except PC1. Save file with individuals as CicadaNAPC1 for simplicity.


##HISTOGRAM and qqplots for PC1 WITHOUT TRANSFORMATION
#Load PC1 data
n<-read.csv("C:/Users/Owner/Desktop/CicadaNAPC1.csv",colClasses="character")
head(n)
n$PC1<-as.numeric(n$PC1)
#Switch sign in PC1 so slow explorers have a smaller PC1 and fast explorers have a larger PC
n$PC1<-n$PC1*-1
hist(n$PC1,xlab="PC1",breaks=20,main="")
qqnorm(n$PC1,main="")
qqline(n$PC1)
#Data is NOT normally distributed!!!!So use MODIFIED VERSION OF Nakagawa's code that I used in Roth et al. 2021 JAE 

#Residual Plot for PC1 WITHOUT TRANSFORMATION
library(lme4)
Naka<-lmer(PC1~(1|Assay)+(1|DaysSinceMay15)+(1|Time)+(1|Temp)+(1|Humidity)+(1|DaysSinceEmerge)+(1|ID),data=n)
summary(Naka)
library(ggplot2) 
head(fortify(Naka))
residPlot<-ggplot(aes(x=.fitted,y=.resid),data=Naka)+geom_point()+geom_hline(yintercept=0)+labs(x="Fitted Values",y="Residuals")
residPlot

###########Test for REPEATABILITY for NA Assay
###############REPEATABILITIES WITH NAKAGAWA'S *MODIEFIED* CODE###############
#This is the final and finished method!!!!

remove.packages("lmeresampler")
install.packages("lmeresampler")
library(lmeresampler)

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
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  tidy = TRUE,
  echo = TRUE,
  fig.width = 8
)
# clearning up
rm(list=ls())

# loading packages
#Install pacman
pacman::p_load(tidyverse, purrr, cplm, lme4, lmeresampler, boot)

# getting functions
source("/Users/Owner/Desktop/function2.R", chdir = TRUE)

#```


## Data preparation

#```{r prep}

#The dataset for the cicadas (males only is "CicadaNAPC1" with the file loaded into R as:

n<-read.csv("C:/Users/Owner/Desktop/CicadaNAPC1.csv",colClasses="character")

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
summary(lmer100)

#```


### Getting correct CI by boostrapping (including repeatability)
#Reinstall dplyr (Shinichi's code wasn't working because I was using an old version of dplyr) 
install.packages("dplyr")
library(dplyr)

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
n<-read.csv("C:/Users/Owner/Desktop/CicadaNAPC1.csv",colClasses="character")
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

#######IMPORTANT: See below all of this AICc stuff for an easier way to do model selection than manually figuring out every possible model combination 
#######IMPORTANT: See below all of this AICc stuff for an easier way to do model selection than manually figuring out every possible model combination 

#Install and load AICcmodavg
library(AICcmodavg)

#Use: Assay, DaysSinceMay15, Time, Temp, Humidity, DaysSinceEmerge, and Individual identity as fixed effects

#Full model
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Time+Temp+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

#5 variables (in addition to ID)
#Should be 6 combinations because (6!)/(5!*((6-5)!))= 6
#https://socratic.org/questions/how-many-different-three-member-teams-can-be-formed-from-six-students
AICc((glm(PC1n~ID+DaysSinceMay15+Time+Temp+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+Time+Temp+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Temp+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Time+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Time+Temp+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Time+Temp+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

#4 variables (in addition to ID)
#Should be 15 combinations because (6!)/(4!*((6-4)!))=15
#https://socratic.org/questions/how-many-different-three-member-teams-can-be-formed-from-six-students
AICc((glm(PC1n~ID+Time+Temp+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceMay15+Temp+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceMay15+Time+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceMay15+Time+Temp+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceMay15+Time+Temp+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Assay+Temp+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+Time+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+Time+Temp+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+Time+Temp+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Temp+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Temp+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Time+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Time+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Time+Temp,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)


#3 variables (in addition to ID)
#Should be 20 combinations because (6!)/(3!*((6-3)!))=20
#https://socratic.org/questions/how-many-different-three-member-teams-can-be-formed-from-six-students
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Time,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Temp,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+DaysSinceMay15+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Assay+Time+Temp,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+Time+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+Time+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Assay+Temp+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+Temp+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Assay+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+DaysSinceMay15+Time+Temp,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceMay15+Time+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceMay15+Time+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+DaysSinceMay15+Temp+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceMay15+Temp+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+DaysSinceMay15+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Time+Temp+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Time+Temp+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Time+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Temp+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)


#2 Variables (in addition to ID)
#Should be 15 combinations because (6!)/(2!*((6-2)!))=15
#https://socratic.org/questions/how-many-different-three-member-teams-can-be-formed-from-six-students
AICc((glm(PC1n~ID+Assay+DaysSinceMay15,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+Time,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+Temp,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Assay+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+DaysSinceMay15+Time,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceMay15+Temp,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceMay15+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceMay15+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Time+Temp,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Time+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Time+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Temp+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Temp+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)

AICc((glm(PC1n~ID+Humidity+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)


#1 Variable (in addition to ID)
#Should be 5 combinations because (6!)/(1!*((6-1)!))=5
#https://socratic.org/questions/how-many-different-three-member-teams-can-be-formed-from-six-students
AICc((glm(PC1n~ID+Assay,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceMay15,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Time,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Temp,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+Humidity,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)
AICc((glm(PC1n~ID+DaysSinceEmerge,data=n,family = gaussian)), return.K = FALSE, second.ord = TRUE, nobs = NULL, c.hat = 1)


#######IMPORTANT: There was an easier way to do model selection than manually figuring out every possible model combination 
#######IMPORTANT: There was an easier way to do model selection than manually figuring out every possible model combination 
#USE BELOW METHOD FOR MODEL SELECTION FOR TONIC IMMOBILITY
#Method use the MuMIN package: https://stackoverflow.com/questions/28606549/how-to-run-lm-models-using-all-possible-combinations-of-several-variables-and-a/52300594
#Also note this link to set the na.action argument and make sure there are no errors when using the dredge function: https://stackoverflow.com/questions/25281739/dredge-function-error-r-package-mumln  
install.packages("MuMIn")
library(MuMIn)
require(MuMIn)
globalmodel <- lm(PC1n~ID+Assay+DaysSinceMay15+Time+Temp+Humidity+DaysSinceEmerge,data=n,na.action = "na.fail")
combinations <- dredge(globalmodel)
print(combinations)
#Now from this printed list find the model with the lowest AICc that also INCLUDES Individual ID (since ID always needs to be kept in the model)
#Note this result matches up with what we found foing things the hard way (i.e., manually)
#The model with ID, Temperature, and Time is the best fit model 

####Model with lowest AICc is:
#NA1<-glm(PC1n~ID+Time+Temp,data=n,family = gaussian)
#summary(NA1)
#The number that is listed under each individual should be added to the value listed under "Intercept" to obtain the resulting EB score. 
#For example for Individual 102 the EB score would be 4.640e-02 + 3.051e+00 = 3.0974. 
#Furthermore, we  tested 99 individuals, but the model output only lists 98 or those 99 individuals. 
#Given this, the value found under intercept is the EB score for the missing individual (i.e. the individual who is not listed here).


#SHINICHI says "you probably want to run your model with = ~ ID - 1 - then you can get SE for each individual. You should now see, ID101, ID102, ID103, ID104... (no Intercept)"
NA1<-lm(PC1n~ID - 1 +Time+Temp,data=n)
summary(NA1)


#Extract the coefficients data frame so we can write the output into a CSV file
results_dfm <-summary.lm(NA1)$coefficients
#Note: Command would be "summary.glm" if left in glm with guassian distribution format 
results_dfm
#write.csv(results_dfm, file = "LmEBCicadas.csv")
#DON'T NEED TO DO THIS ANYMORE THANKS TO SHINICHI'S SUGGESTIONS: Add intercept to all individual's estimates to get EB Score for each individual
#DON'T NEED TO DO THIS ANYMORE THANKS TO SHINICHI'S SUGGESTIONS:Find missing individual (it is male 101 here)
#Output shows faster explores as having higher values 














#######Tonic Immobility: MALES#######
tm<-read.csv("C:/Users/Owner/Desktop/CicadaRawTIforPCAm.csv",colClasses="character")
head(tm)
#change number columns of csv load to numeric
tm$LatencyToFreeze<-as.numeric(tm$LatencyToFreeze)
tm$TimeFrozen <-as.numeric(tm$TimeFrozen)
tm$TimesFlipped<-as.numeric(tm$TimesFlipped)
tm$ID<-as.factor(tm$ID)

#Find out how many unique individuals were tested
length(unique(tm[["ID"]]))

#Standardize (Z-transform) variables
#https://stats.stackexchange.com/questions/69157/why-do-we-need-to-normalize-data-before-principal-component-analysis-pca  
#https://stackoverflow.com/questions/15215457/standardize-data-columns-in-r
#The most common normalization is the z-transformation, where you subtract the mean and divide by the standard deviation of your variable. The result will have mean=0 and sd=1.
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
#0.5395744 of variance is explained bt PC1 - originally not sure if this was good enoug - YES IT IS GOOD ENOUGH as it matches up with the amount of variance explained by PC1 for Red Junglefowl Novel Object Assay!!!!

biplot(pcAlltm)
pcAlltm$loadings
AllScorestm<-pcAlltm$scores
AllScorestm

#See loadings
PCAprtm <- prcomp(SubAlltm,center = TRUE,scale. = TRUE) 
####IMPORTANT Different variables are loaded strongly and consistently
print(PCAprtm)
plot(PCAprtm, type = "l")
summary(PCAprtm)
biplot(PCAprtm)

#write.csv(AllScorestm, file = "CicadaTIPC1m.csv")
#Have to add Individual IDs in to the CSV file (Copy and paste "ID" column from "CicadaRawTIforPCAm.csv" into "CicadaTIPC1m.csv"). Remove everything except PC1 (also kept in PC2 here - decided PC2 was not necessary (see above)). Save file with individuals as CicadaTIPC1m for simplicity.


##HISTOGRAM and qqplots for PC1 WITHOUT TRANSFORMATION
#Load PC1 data
ntm<-read.csv("C:/Users/Owner/Desktop/CicadaTIPC1m.csv",colClasses="character")
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

remove.packages("lmeresampler")
install.packages("lmeresampler")
library(lmeresampler)

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
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  tidy = TRUE,
  echo = TRUE,
  fig.width = 8
)
# clearning up
rm(list=ls())

# loading packages
#Install pacman
pacman::p_load(tidyverse, purrr, cplm, lme4, lmeresampler, boot)

# getting functions
source("/Users/Owner/Desktop/function2.R", chdir = TRUE)

#```


## Data preparation

#```{r prep}

#The dataset for the cicadas (males only is "CicadaTIPC1m" with the file loaded into R as:

ntm<-read.csv("C:/Users/Owner/Desktop/CicadaTIPC1m.csv",colClasses="character")

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

#```


### Getting correct CI by boostrapping (including repeatability)
#Reinstall dplyr (Shinichi's code wasn't working because I was using an old version of dplyr) 
install.packages("dplyr")
library(dplyr)

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
ntm<-read.csv("C:/Users/Owner/Desktop/CicadaTIPC1m.csv",colClasses="character")
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

#USE BELOW METHOD FOR MODEL SELECTION FOR TONIC IMMOBILITY
#Method use the MuMIN package: https://stackoverflow.com/questions/28606549/how-to-run-lm-models-using-all-possible-combinations-of-several-variables-and-a/52300594
#Also note this link to set the na.action argument and make sure there are no errors when using the dredge function: https://stackoverflow.com/questions/25281739/dredge-function-error-r-package-mumln  
install.packages("MuMIn")
library(MuMIn)
require(MuMIn)
globalmodel <- lm(PC1~ID+Assay+DaysSinceMay15+Time+Temp+Humidity+DaysSinceEmerge+Observer,data=ntm,na.action = "na.fail")
combinations <- dredge(globalmodel)
print(combinations)
#Now from this printed list find the model with the lowest AICc that also INCLUDES Individual ID (since ID always needs to be kept in the model)
#The model with ID and Temperature is the best fit model

####Model with lowest AICc is:
#TI1m<-lm(PC1~ID+Temp,data=ntm)
#summary(TI1m)
#SHINICHI says "you probably want to run your model with = ~ ID - 1 - then you can get SE for each individual. You should now see, ID101, ID102, ID103, ID104... (no Intercept)"
TI1m<-lm(PC1~ID - 1 +Temp,data=ntm)
summary(TI1m)


#Extract the coefficients data frame so we can write the output into a CSV file
results_dfm <-summary.lm(TI1m)$coefficients
results_dfm
#write.csv(results_dfm, file = "LmTICicadasM.csv")
#Output shows Bolder Individuals as having higher values 














#######Tonic Immobility: FEMALES#######
tf<-read.csv("C:/Users/Owner/Desktop/CicadaRawTIforPCAf.csv",colClasses="character")
head(tf)
#change number columns of csv load to numeric
tf$LatencyToFreeze<-as.numeric(tf$LatencyToFreeze)
tf$TimeFrozen <-as.numeric(tf$TimeFrozen)
tf$TimesFlipped<-as.numeric(tf$TimesFlipped)
tf$ID<-as.factor(tf$ID)

#Find out how many unique individuals were tested
length(unique(tf[["ID"]]))

#Standardize (Z-transform) variables
#https://stats.stackexchange.com/questions/69157/why-do-we-need-to-normalize-data-before-principal-component-analysis-pca  
#https://stackoverflow.com/questions/15215457/standardize-data-columns-in-r
#The most common normalization is the z-transformation, where you subtract the mean and divide by the standard deviation of your variable. The result will have mean=0 and sd=1.
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
#0.5535 of variance is explained by PC1 - originally not sure if this was good enoug - YES IT IS GOOD ENOUGH as it matches up with the amount of variance explained by PC1 for Red Junglefowl Novel Object Assay!!!!

biplot(pcAlltf)
pcAlltf$loadings
AllScorestf<-pcAlltf$scores
AllScorestf

#See loadings
PCAprtf <- prcomp(SubAlltf,center = TRUE,scale. = TRUE) 
####IMPORTANT Different variables are loaded strongly and consistently
print(PCAprtf)
plot(PCAprtf, type = "l")
summary(PCAprtf)
biplot(PCAprtf)

#write.csv(AllScorestf, file = "CicadaTIPC1f.csv")
#Have to add Individual IDs in to the CSV file (Copy and paste "ID" column from "CicadaRawTIforPCAf.csv" into "CicadaTIPC1f.csv"). Remove everything except PC1 (also kept in PC2 here - decided PC2 was not necessary (see above)). Save file with individuals as CicadaTIPC1f for simplicity.


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
#This is the final and finished method!!!!

remove.packages("lmeresampler")
install.packages("lmeresampler")
library(lmeresampler)

#---
#title: "Figuring out how to get CI with non-parametric boostraps with linear mixed models"
#author: "Shinichi Nakagawa"
#date: "09/08/2018"
#output: htfl_document
#---

## Background

#We will need to get confidence intervals (CIs) for repeatability for Tweedie distributed data usimng the `cplm` package. We will use chicken behaviour data wiht adjusted PC scores. Actually, it is easy to get "credible" intervals using the `bcplm` function so we will just do that. However, it turns out moding with `cplm` does not really work (see below). So we will use non-parametric bootstrapping and normal linear models (`lme4::lmer`). In this way, we can violate linear model assumptions but we can still get meaningful CIs from `lmer`. 

## Setting up

#```{r setup, echo=FALSE}
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  tidy = TRUE,
  echo = TRUE,
  fig.width = 8
)
# clearning up
rm(list=ls())

# loading packages
#Install pacman
pacman::p_load(tidyverse, purrr, cplm, lme4, lmeresampler, boot)

# getting functions
source("/Users/Owner/Desktop/function2.R", chdir = TRUE)

#```


## Data preparation

#```{r prep}

#The dataset for the cicadas (FEMALEs only is "CicadaTIPC1f" with the file loaded into R as:

ntf<-read.csv("C:/Users/Owner/Desktop/CicadaTIPC1f.csv",colClasses="character")

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
#Reinstall dplyr (Shinichi's code wasn't working because I was using an old version of dplyr) 
install.packages("dplyr")
library(dplyr)

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
ntf<-read.csv("C:/Users/Owner/Desktop/CicadaTIPC1f.csv",colClasses="character")
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

#USE BELOW METHOD FOR MODEL SELECTION FOR TONIC IMMOBILITY
#Method use the MuMIN package: https://stackoverflow.com/questions/28606549/how-to-run-lm-models-using-all-possible-combinations-of-several-variables-and-a/52300594
#Also note this link to set the na.action argument and make sure there are no errors when using the dredge function: https://stackoverflow.com/questions/25281739/dredge-function-error-r-package-mumln  
install.packages("MuMIn")
library(MuMIn)
require(MuMIn)
globalmodel <- lm(PC1~ID+Assay+DaysSinceMay15+Time+Temp+Humidity+DaysSinceEmerge+Observer,data=ntf,na.action = "na.fail")
combinations <- dredge(globalmodel)
print(combinations)
#Now from this printed list find the model with the lowest AICc that also INCLUDES Individual ID (since ID always needs to be kept in the model)
#The model with ID and Humidity is the best fit model 

####Model with lowest AICc is:
#TI1f<-lm(PC1~ID+Humidity,data=ntf)
#summary(TI1f)
#SHINICHI says "you probably want to run your model with = ~ ID - 1 - then you can get SE for each individual. You should now see, ID101, ID102, ID103, ID104... (no Intercept)"
TI1f<-lm(PC1~ID - 1 +Humidity,data=ntf)
summary(TI1f)


#Extract the coefficients data frame so we can write the output into a CSV file
results_dfm <-summary.lm(TI1f)$coefficients
results_dfm
#write.csv(results_dfm, file = "LmTICicadasF.csv")
#Output shows Bolder Individuals as having higher values







########Are EB and TI correlated?
BS<-read.csv("/Users/Owner/Desktop/CicadaBehavSynd.csv",colClasses="character")
BS$TI<-as.numeric(BS$TI)
BS$EB<-as.numeric(BS$EB)
cor.test(BS$EB,BS$TI)
#No significant correlation between EB and TI







#########MATING ANALYSES#########
#########MATING ANALYSES#########
#########MATING ANALYSES#########
#########MATING ANALYSES#########
library(lmerTest)
install.packages("glmmTMB")
library(glmmTMB)
I<-read.csv("/Users/Owner/Desktop/CicadaCollatedMatingData.csv",colClasses="character")
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
#-8.230178353 to -2.12536741
range(dead$EB, na.rm = TRUE)
#-7.240319661 to -4.04970104

hist(alive$EB, breaks=5,xlab="Exploration Behavior")
hist(dead$EB, breaks=5,xlab="Exploration Behavior")

shapiro.test(alive$EB)
#Normal enough (p > 0.05)
shapiro.test(dead$EB)
#Normal enough (p > 0.05)
##Because EB measures are "normal enough" for both alive and dead males, we can use a upaired 2 tailed t-test (which is a parametric test)
#See Table 3: https://www.biochemia-medica.com/en/journal/20/1/10.11613/BM.2010.004/fullArticle
install.packages("ggpubr")
library(ggpubr)
t.test(alive$EB, dead$EB, alternative = "two.sided", var.equal = FALSE)
#NO SIGNIFICANT DIFFERENCE IN EB BETWEEN DEAD AND ALIVE MALES - so it's ok to leave them in the analyses



####Check to see TI distributions of males that died during the mating trials compared to those that did NOT die 
range(alive$TI, na.rm = TRUE)
#-3.6777105 to 0.2359917
range(dead$TI, na.rm = TRUE)
#-3.8510263 to -0.2984309

hist(alive$TI, breaks=5,xlab="Tonic Immobility")
hist(dead$TI, breaks=5,xlab="Tonic Immobility")

shapiro.test(alive$TI)
#NOT normal enough (p < 0.05)
shapiro.test(dead$TI)
#Normal enough (p > 0.05)
##Because TI measures are NOT "normal enough" for alive males, we will use a Mann-Whitney U Test (which is a NONparametric test)
#See Table 3: https://www.biochemia-medica.com/en/journal/20/1/10.11613/BM.2010.004/fullArticle
#https://rcompanion.org/handbook/F_04.html
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
#524 to 924
range(dead$Weight, na.rm = TRUE)
#552 to 799

hist(alive$Weight, breaks=5,xlab="Tonic Immobility")
hist(dead$Weight, breaks=5,xlab="Tonic Immobility")

shapiro.test(alive$Weight)
#NOT normal enough (p < 0.05)
shapiro.test(dead$Weight)
#Normal enough (p > 0.05)
##Because Weight measures are NOT "normal enough" for alive males, we will use a Mann-Whitney U Test (which is a NONparametric test)
#See Table 3: https://www.biochemia-medica.com/en/journal/20/1/10.11613/BM.2010.004/fullArticle
#https://rcompanion.org/handbook/F_04.html
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
#Weight by whether they died midway through the day (i.e., before the evening) 
#Not sure this is necessary since no significnt difference between alive and dead males in EB or TI (also 3/11 males that died during the mating trial were the ACTORS of some interaction (compared to 78/108 males that did not die))
#Either way though, do NOT resolve by completely removing dead males since 3/11 males that died during the mating trial were the ACTORS of some interaction
#Test1<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, weight = Died, family = binomial)
#summary(Test1) 
#SHINICHI DOESN'T THINK WEIGHTING OUR MODELS BY DIED VS NOT DIED IS NECESSARY 
#Updated unweighted models
Test1<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, family = binomial)
summary(Test1) 

#Currently Pen A is the reference class
#Change Pen B to reference class
I$Pen[I$Pen == "B"] <- "1B"
#Rerun analyses
Test1<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, family = binomial)
summary(Test1) 
#Switch Pen 1B back to Pen B
I$Pen[I$Pen == "1B"] <- "B"
#Change Pen C to reference class
I$Pen[I$Pen == "C"] <- "1C"
#Rerun analyses
Test1<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, family = binomial)
summary(Test1) 
#Switch Pen 1C back to Pen C
I$Pen[I$Pen == "1C"] <- "C"
#Change Pen D to reference class
I$Pen[I$Pen == "D"] <- "1D"
#Rerun analyses
Test1<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, family = binomial)
summary(Test1) 
#Switch Pen 1D back to Pen D
I$Pen[I$Pen == "1D"] <- "D"

#LRT Location
Test1<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, family = binomial)
Test1n<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Date + Pen + (1|ID) + (1|Group), data=I, family = binomial)
anova(Test1, Test1n, test="Chisq") 

#LRT Date
Test1<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, family = binomial)
Test1n<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Location + Pen + (1|ID) + (1|Group), data=I, family = binomial)
anova(Test1, Test1n, test="Chisq") 

#LRT Pen
Test1<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, family = binomial)
Test1n<-glmmTMB(Mated ~  EB.sc + TI.sc + Weight.sc + Location + Date + (1|ID) + (1|Group), data=I, family = binomial)
anova(Test1, Test1n, test="Chisq") 







###########M-M ATTEMPTED COPULATIONS###########
hist(I$MM, breaks=20,xlab="M-M Attempted Copulations")
#Will need to use a zero inflated Poisson
#Weight by whether they died midway through the day
#Test2<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, weight = Died, ziformula=~1, family = poisson)
#summary(Test2) 
#SHINICHI DOESN'T THINK WEIGHTING OUR MODELS BY DIED VS NOT DIED IS NECESSARY 
#Updated unweighted models
Test2<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
summary(Test2) 

#Currently Pen A is the reference class
#Change Pen B to reference class
I$Pen[I$Pen == "B"] <- "1B"
#Rerun analyses
Test2<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
summary(Test2)
#Switch Pen 1B back to Pen B
I$Pen[I$Pen == "1B"] <- "B"
#Change Pen C to reference class
I$Pen[I$Pen == "C"] <- "1C"
#Rerun analyses
Test2<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
summary(Test2)
#Switch Pen 1C back to Pen C
I$Pen[I$Pen == "1C"] <- "C"
#Change Pen D to reference class
I$Pen[I$Pen == "D"] <- "1D"
#Rerun analyses
Test2<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
summary(Test2)
#Switch Pen 1D back to Pen D
I$Pen[I$Pen == "1D"] <- "D"

#LRT Location
Test2<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
Test2n<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
anova(Test2, Test2n, test="Chisq") 

#LRT Date
Test2<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
Test2n<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Location + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
anova(Test2, Test2n, test="Chisq") 

#LRT Pen
Test2<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
Test2n<-glmmTMB(MM ~ EB.sc + TI.sc + Weight.sc + Location + Date + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
anova(Test2, Test2n, test="Chisq") 







###########ALL (BOTH SUCCESSFUL AND UNSUCCESSFUL) ATTEMPTED COPULATIONS WITH EITHER SEX###########
I$All<-I$Mated+I$MU+I$MM
Test0<-glmmTMB(All ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
summary(Test0)

#Currently Pen A is the reference class
#Change Pen B to reference class
I$Pen[I$Pen == "B"] <- "1B"
#Rerun analyses
Test0<-glmmTMB(All ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
summary(Test0) 
#Switch Pen 1B back to Pen B
I$Pen[I$Pen == "1B"] <- "B"
#Change Pen C to reference class
I$Pen[I$Pen == "C"] <- "1C"
#Rerun analyses
#Test0<-glmmTMB(All ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
summary(Test0) 
#Switch Pen 1C back to Pen C
I$Pen[I$Pen == "1C"] <- "C"
#Change Pen D to reference class
I$Pen[I$Pen == "D"] <- "1D"
#Rerun analyses
Test0<-glmmTMB(All ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
summary(Test0) 
#Switch Pen 1D back to Pen D
I$Pen[I$Pen == "1D"] <- "D"

#LRT Location
Test0<-glmmTMB(All ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
Test0n<-glmmTMB(All ~ EB.sc + TI.sc + Weight.sc + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
anova(Test0, Test0n, test="Chisq") 

#LRT Date
Test0<-glmmTMB(All ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
Test0n<-glmmTMB(All ~ EB.sc + TI.sc + Weight.sc + Location + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
anova(Test0, Test0n, test="Chisq") 

#LRT Pen
Test0<-glmmTMB(All ~ EB.sc + TI.sc + Weight.sc + Location + Date + Pen + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
Test0n<-glmmTMB(All ~ EB.sc + TI.sc + Weight.sc + Location + Date + (1|ID) + (1|Group), data=I, ziformula=~1, family = poisson)
anova(Test0, Test0n, test="Chisq") 







###########M-M ADVERTISEMENT RATE###########
A<-read.csv("/Users/Owner/Desktop/CicadaCollatedAdvertData.csv",colClasses="character")
head(A)
A$ID<-as.factor(A$ID)
A$DurationCI<-as.numeric(A$DurationCI)
A$NumberAdvert<-as.numeric(A$NumberAdvert)
A$Weight<-as.numeric(A$Weight)
A$EB<-as.numeric(A$EB)
A$TI<-as.numeric(A$TI)
#SCALE VARIABLES TO TRY TO STOP CONVERGENCE WARNINGS
#https://stackoverflow.com/questions/41766181/correct-way-to-scale-for-multilevel-regression-using-lmer-r
A <- transform(A,Weight.sc=scale(Weight),EB.sc=scale(EB),TI.sc=scale(TI))

hist(A$DurationCI)
shapiro.test(A$DurationCI)
#Normal enough (p > 0.05)
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



