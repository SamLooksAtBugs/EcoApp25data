#packages ####
library(tidyverse)
library(readxl)
library(vegan)
library(lme4)
library(ggstatsplot)
library(MASS)
library(car)
library(gridExtra)
library(ordinal)
library(lmerTest)
library(lmeresampler)

setwd("C:/Users/Sam/OneDrive/Documents")
#Read in data, substitute pathways if needed
genera <- read_excel("ecoBI.xlsx", sheet="Genus2")
family <- read_excel("ecoBI.xlsx", sheet="Family2")
allmech <- read_excel("ArkansasValleyOP.xlsx", sheet = "All")
testdata <- read_excel("indextestdata.xlsx", sheet="Sheet1")

#testdata <- testdata[,-c(22:25)]

#Remove Chironomidae cause they just cause issues
testdata <- testdata %>% filter(Family!="Chironomidae")
hist(testdata$TN)
#Assign 
generaAll <- genera %>% filter(Ecoregion == "All")
generaAll$scoreTN<- as.integer(ntile(generaAll$TKN,10))
generaAll$scoreTP<- as.integer(ntile(generaAll$TP,10))

familyAll <- family %>% filter(Ecoregion == "All")
familyAll$scoreTN<- as.integer(ntile(familyAll$TKN,10))
familyAll$scoreTP<- as.integer(ntile(familyAll$TP,10))

#Combine data
allindexscores <- rbind(familyAll,generaAll)
allindexscores <- allindexscores[,c(1,2,17:ncol(allindexscores))]

#transform environmental variables
testdata$TP <- log10(testdata$TP)
testdata$TN <- log10(testdata$TN)

#Remove outliers
testdata$TN[testdata$TN<(quantile(testdata$TN,0.25,na.rm=T)-(1.5*(quantile(testdata$TN,0.75,na.rm=T)-quantile(testdata$TN,0.25,na.rm=T))))|testdata$TN>(quantile(testdata$TN,0.75,na.rm=T)+(1.5*(quantile(testdata$TN,0.75,na.rm=T)-quantile(testdata$TN,0.25,na.rm=T))))]<-NA
testdata$TP[testdata$TP<(quantile(testdata$TP,0.25,na.rm=T)-(1.5*(quantile(testdata$TP,0.75,na.rm=T)-quantile(testdata$TP,0.25,na.rm=T))))|testdata$TP>(quantile(testdata$TP,0.75,na.rm=T)+(1.5*(quantile(testdata$TP,0.75,na.rm=T)-quantile(testdata$TP,0.25,na.rm=T))))]<-NA

#Combine the test data and the scores
allindexscores <- left_join(testdata,allindexscores, by="Species")

#also good to add the mechanisms to the dataset, but first we create a new
#column, dispersal, which factors in whether they disperse aquatically
#or terrestrially
allmech$Dispersal <- allmech$SFP*ifelse(allmech$DispersalType=="Aquatic",0.5,1)
allindexscores <- left_join(allindexscores, allmech, by="Species")

#Begin index calculation: Score x log-transformed Count
allindexscores <- allindexscores %>% mutate(TNscore = scoreTN*log(1+Count))
allindexscores <- allindexscores %>% mutate(TPscore =scoreTP*log(1+Count))
allindexscores <- allindexscores %>% mutate(logcount = log(1+Count))

#Correct counts based on assembly mechanisms
{allindexscores$logcount.D <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)
  allindexscores$logcount.I <- allindexscores$logcount*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)
  allindexscores$logcount.N<- allindexscores$logcount*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DI <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)
  allindexscores$logcount.DN <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.IN <- allindexscores$logcount*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DNI <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.D.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)
  allindexscores$logcount.I.5 <- allindexscores$logcount*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)
  allindexscores$logcount.N.5<- allindexscores$logcount*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5I <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)
  allindexscores$logcount.DI.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)
  allindexscores$logcount.D.5N <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DN.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.I.5N <- allindexscores$logcount*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.IN.5 <- allindexscores$logcount*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5NI <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DN.5I <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DNI.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5N.5I <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DN.5I.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5NI.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5N.5I.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5I.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)
  allindexscores$logcount.D.5N.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.I.5N.5 <- allindexscores$logcount*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)}

#TP site scores
{allindexscores$TPscore.D <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)
  allindexscores$TPscore.I <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)
  allindexscores$TPscore.N<- allindexscores$TPscore*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DI <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)
  allindexscores$TPscore.DN <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.IN <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DNI <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.D.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)
  allindexscores$TPscore.I.5 <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)
  allindexscores$TPscore.N.5<- allindexscores$TPscore*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5I <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)
  allindexscores$TPscore.DI.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)
  allindexscores$TPscore.D.5N <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DN.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.I.5N <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.IN.5 <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5NI <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DN.5I <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DNI.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5N.5I <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DN.5I.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5NI.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5N.5I.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5I.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)
  allindexscores$TPscore.D.5N.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.I.5N.5 <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
}

#TN site scores
{allindexscores$TNscore.D <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)
  allindexscores$TNscore.I <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)
  allindexscores$TNscore.N<- allindexscores$TNscore*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DI <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)
  allindexscores$TNscore.DN <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.IN <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DNI <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.D.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)
  allindexscores$TNscore.I.5 <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)
  allindexscores$TNscore.N.5<- allindexscores$TNscore*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5I <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)
  allindexscores$TNscore.DI.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)
  allindexscores$TNscore.D.5N <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DN.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.I.5N <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.IN.5 <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5NI <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DN.5I <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DNI.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5N.5I <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DN.5I.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5NI.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5N.5I.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5I.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)
  allindexscores$TNscore.D.5N.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=allindexscores$DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.I.5N.5 <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=allindexscores$CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
}

#Summarize measured values, dates, etc. 

testdata1 <- allindexscores[,c(1,6,7,12:21)] %>% group_by(UID) %>% summarise(across(.fns=mean, na.rm=T))
#unique_data <- distinct(allindexscores[,1:2])
#testdata1 <- left_join(testdata1,unique_data, by="UID")
testdata2 <- allindexscores[,c(1,31:ncol(allindexscores))] %>% group_by(UID) %>% summarise(across(.fns=sum, na.rm=T))

#View(testfinal)
testfinal <- left_join(testdata1,testdata2,by="UID")

#TP final score
{testfinal$TPscore <- testfinal$TPscore/testfinal$logcount
testfinal$TPscore.D <- testfinal$TPscore.D/testfinal$logcount.D
testfinal$TPscore.I <- testfinal$TPscore.I/testfinal$logcount.I
testfinal$TPscore.N<- testfinal$TPscore.N/testfinal$logcount.N
testfinal$TPscore.DI <- testfinal$TPscore.DI/testfinal$logcount.DI
testfinal$TPscore.DN <- testfinal$TPscore.DN/testfinal$logcount.DN
testfinal$TPscore.IN <- testfinal$TPscore.IN/testfinal$logcount.IN
testfinal$TPscore.DNI <- testfinal$TPscore.DNI/testfinal$logcount.DNI
testfinal$TPscore.D.5 <- testfinal$TPscore.D.5/testfinal$logcount.D.5
testfinal$TPscore.I.5 <- testfinal$TPscore.I.5/testfinal$logcount.I.5
testfinal$TPscore.N.5<- testfinal$TPscore.N.5/testfinal$logcount.N.5
testfinal$TPscore.D.5I <- testfinal$TPscore.D.5I/testfinal$logcount.D.5I
testfinal$TPscore.DI.5 <- testfinal$TPscore.DI.5/testfinal$logcount.DI.5
testfinal$TPscore.D.5N <- testfinal$TPscore.D.5N/testfinal$logcount.D.5N
testfinal$TPscore.DN.5 <- testfinal$TPscore.DN.5/testfinal$logcount.DN.5
testfinal$TPscore.I.5N <- testfinal$TPscore.I.5N/testfinal$logcount.I.5N
testfinal$TPscore.IN.5 <- testfinal$TPscore.IN.5/testfinal$logcount.IN.5
testfinal$TPscore.D.5NI <- testfinal$TPscore.D.5NI/testfinal$logcount.D.5NI
testfinal$TPscore.DN.5I <- testfinal$TPscore.DN.5I/testfinal$logcount.DN.5I
testfinal$TPscore.DNI.5 <- testfinal$TPscore.DNI.5/testfinal$logcount.DNI.5
testfinal$TPscore.D.5N.5I <- testfinal$TPscore.D.5N.5I/testfinal$logcount.D.5N.5I
testfinal$TPscore.DN.5I.5 <- testfinal$TPscore.DN.5I.5/testfinal$logcount.DN.5I.5
testfinal$TPscore.D.5NI.5 <- testfinal$TPscore.D.5NI.5/testfinal$logcount.D.5NI.5
testfinal$TPscore.D.5N.5I.5 <- testfinal$TPscore.D.5N.5I.5/testfinal$logcount.D.5N.5I.5
testfinal$TPscore.D.5I.5 <- testfinal$TPscore.D.5I.5/testfinal$logcount.D.5I.5
testfinal$TPscore.D.5N.5 <- testfinal$TPscore.D.5N.5/testfinal$logcount.D.5N.5
testfinal$TPscore.I.5N.5 <- testfinal$TPscore.I.5N.5/testfinal$logcount.I.5N.5
}

#TN final score
{testfinal$TNscore <- testfinal$TNscore/testfinal$logcount
testfinal$TNscore.D <- testfinal$TNscore.D/testfinal$logcount.D
testfinal$TNscore.I <- testfinal$TNscore.I/testfinal$logcount.I
testfinal$TNscore.N<- testfinal$TNscore.N/testfinal$logcount.N
testfinal$TNscore.DI <- testfinal$TNscore.DI/testfinal$logcount.DI
testfinal$TNscore.DN <- testfinal$TNscore.DN/testfinal$logcount.DN
testfinal$TNscore.IN <- testfinal$TNscore.IN/testfinal$logcount.IN
testfinal$TNscore.DNI <- testfinal$TNscore.DNI/testfinal$logcount.DNI
testfinal$TNscore.D.5 <- testfinal$TNscore.D.5/testfinal$logcount.D.5
testfinal$TNscore.I.5 <- testfinal$TNscore.I.5/testfinal$logcount.I.5
testfinal$TNscore.N.5<- testfinal$TNscore.N.5/testfinal$logcount.N.5
testfinal$TNscore.D.5I <- testfinal$TNscore.D.5I/testfinal$logcount.D.5I
testfinal$TNscore.DI.5 <- testfinal$TNscore.DI.5/testfinal$logcount.DI.5
testfinal$TNscore.D.5N <- testfinal$TNscore.D.5N/testfinal$logcount.D.5N
testfinal$TNscore.DN.5 <- testfinal$TNscore.DN.5/testfinal$logcount.DN.5
testfinal$TNscore.I.5N <- testfinal$TNscore.I.5N/testfinal$logcount.I.5N
testfinal$TNscore.IN.5 <- testfinal$TNscore.IN.5/testfinal$logcount.IN.5
testfinal$TNscore.D.5NI <- testfinal$TNscore.D.5NI/testfinal$logcount.D.5NI
testfinal$TNscore.DN.5I <- testfinal$TNscore.DN.5I/testfinal$logcount.DN.5I
testfinal$TNscore.DNI.5 <- testfinal$TNscore.DNI.5/testfinal$logcount.DNI.5
testfinal$TNscore.D.5N.5I <- testfinal$TNscore.D.5N.5I/testfinal$logcount.D.5N.5I
testfinal$TNscore.DN.5I.5 <- testfinal$TNscore.DN.5I.5/testfinal$logcount.DN.5I.5
testfinal$TNscore.D.5NI.5 <- testfinal$TNscore.D.5NI.5/testfinal$logcount.D.5NI.5
testfinal$TNscore.D.5N.5I.5 <- testfinal$TNscore.D.5N.5I.5/testfinal$logcount.D.5N.5I.5
testfinal$TNscore.D.5I.5 <- testfinal$TNscore.D.5I.5/testfinal$logcount.D.5I.5
testfinal$TNscore.D.5N.5 <- testfinal$TNscore.D.5N.5/testfinal$logcount.D.5N.5
testfinal$TNscore.I.5N.5 <- testfinal$TNscore.I.5N.5/testfinal$logcount.I.5N.5
}

#last step: remove sites without a score
test <- testfinal %>% filter(TNscore != "NaN")

#there is an issue
#test$Ecoregion.x <- stringr::str_to_title(test$Ecoregion.x)

# test <- test %>% mutate(region = case_when(
#   Ecoregion.x %in% c("Arkansas Valley", "Boston Mountains", "Ouachita Mountains","Ozark Highlands") ~ "Temperate Forest",
#   Ecoregion.x  %in% c("Central Great Plains", "Central Irregular Plains", "Cross Timbers") ~ "Central Plains",
#   Ecoregion.x  %in% c("South Central Plains", "Southwestern Tablelands") ~ "Southern Plains"))
# test <- test %>% mutate(season =case_when(Month %in% c(4,5) ~ "Spring",Month %in% c(6,7,8) ~"Summer",Month %in% c(9,10,11)~ "Fall"))
# 
# 
# spring <- test %>% filter(season=="Spring")
# summer <- test %>% filter(season=="Summer")
# fall <- test %>% filter(season=="Fall")
# =
# south <- test %>% filter(region=="Southern Plains")
# central <- test %>% filter(region=="Central Plains")
# forest <- test %>% filter(region=="Temperate Forest")
# 
# cor(spring$TN,spring$TNscore, use = "na.or.complete")

#TP models
{
  TPModel.null <- lm(TP~1,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore','TPscore.I')]))
  TPModel <-lm(TP~1+TPscore ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore','TPscore.I')]))
  TPModel.D <-lm(TP~1+TPscore.D ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D','TPscore.I')]))
  TPModel.I <-lm(TP~1+TPscore.I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I')]))
  TPModel.N<-lm(TP~1+TPscore.N,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.N','TPscore.I')]))
  TPModel.DI <-lm(TP~1+TPscore.DI ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DI')]))
  TPModel.DN <-lm(TP~1+TPscore.DN ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN','TPscore.I')]))
  TPModel.IN <-lm(TP~1+TPscore.IN ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.IN')]))
  TPModel.DNI <-lm(TP~1+TPscore.DNI ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DNI')]))
  TPModel.D.5 <-lm(TP~1+TPscore.D.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5','TPscore.I')]))
  TPModel.I.5 <-lm(TP~1+TPscore.I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I.5','TPscore.I')]))
  TPModel.N.5<-lm(TP~1+TPscore.N.5,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.N.5','TPscore.I')]))
  TPModel.D.5I <-lm(TP~1+TPscore.D.5I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5I')]))
  TPModel.DI.5 <-lm(TP~1+TPscore.DI.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DI.5')]))
  TPModel.D.5N <-lm(TP~1+TPscore.D.5N ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N','TPscore.I')]))
  TPModel.DN.5 <-lm(TP~1+TPscore.DN.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN.5','TPscore.I')]))
  TPModel.I.5N <-lm(TP~1+TPscore.I.5N ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I.5N')]))
  TPModel.IN.5 <-lm(TP~1+TPscore.IN.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.IN.5')]))
  TPModel.D.5NI <-lm(TP~1+TPscore.D.5NI ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5NI')]))
  TPModel.DN.5I <-lm(TP~1+TPscore.DN.5I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN.5I')]))
  TPModel.DNI.5 <-lm(TP~1+TPscore.DNI.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DNI.5')]))
  TPModel.D.5N.5I <-lm(TP~1+TPscore.D.5N.5I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N.5I')]))
  TPModel.DN.5I.5 <-lm(TP~1+TPscore.DN.5I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN.5I.5')]))
  TPModel.D.5NI.5 <-lm(TP~1+TPscore.D.5NI.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5NI.5')]))
  TPModel.D.5N.5I.5 <-lm(TP~1+TPscore.D.5N.5I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N.5I.5')]))
  TPModel.D.5I.5 <-lm(TP~1+TPscore.D.5I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5I.5')]))
  TPModel.D.5N.5 <-lm(TP~1+TPscore.D.5N.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N.5','TPscore.I')]))
  TPModel.I.5N.5 <-lm(TP~1+TPscore.I.5N.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I.5N.5')]))
  TPAICBIC <- cbind(AIC(TPModel,TPModel.D ,TPModel.I ,TPModel.N,TPModel.DI ,TPModel.DN ,TPModel.IN ,
                        TPModel.DNI ,TPModel.D.5 ,TPModel.I.5 ,TPModel.N.5,TPModel.D.5I ,TPModel.DI.5 ,
                        TPModel.D.5N ,TPModel.DN.5 ,TPModel.I.5N ,TPModel.IN.5 ,TPModel.D.5NI ,
                        TPModel.DN.5I ,TPModel.DNI.5 ,TPModel.D.5N.5I ,TPModel.DN.5I.5 ,TPModel.D.5NI.5 ,
                        TPModel.D.5N.5I.5 ,TPModel.D.5I.5 ,TPModel.D.5N.5 ,TPModel.I.5N.5 ,TPModel.null),
                    BIC(TPModel,TPModel.D ,TPModel.I ,TPModel.N,TPModel.DI ,TPModel.DN ,TPModel.IN ,
                        TPModel.DNI ,TPModel.D.5 ,TPModel.I.5 ,TPModel.N.5,TPModel.D.5I ,TPModel.DI.5 ,
                        TPModel.D.5N ,TPModel.DN.5 ,TPModel.I.5N ,TPModel.IN.5 ,TPModel.D.5NI ,
                        TPModel.DN.5I ,TPModel.DNI.5 ,TPModel.D.5N.5I ,TPModel.DN.5I.5 ,TPModel.D.5NI.5 ,
                        TPModel.D.5N.5I.5 ,TPModel.D.5I.5 ,TPModel.D.5N.5 ,TPModel.I.5N.5 ,TPModel.null))}

#TN models
{TNModel.null <- lm(TN~1,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore','TNscore.I')]))
  TNModel <-lm(TN~1+TNscore ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore','TNscore.I')]))
  TNModel.D <-lm(TN~1+TNscore.D ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D','TNscore.I')]))
  TNModel.I <-lm(TN~1+TNscore.I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I')]))
  TNModel.N<-lm(TN~1+TNscore.N,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.N','TNscore.I')]))
  TNModel.DI <-lm(TN~1+TNscore.DI ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DI')]))
  TNModel.DN <-lm(TN~1+TNscore.DN ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN','TNscore.I')]))
  TNModel.IN <-lm(TN~1+TNscore.IN ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.IN')]))
  TNModel.DNI <-lm(TN~1+TNscore.DNI ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DNI')]))
  TNModel.D.5 <-lm(TN~1+TNscore.D.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5','TNscore.I')]))
  TNModel.I.5 <-lm(TN~1+TNscore.I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I.5','TNscore.I')]))
  TNModel.N.5<-lm(TN~1+TNscore.N.5,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.N.5','TNscore.I')]))
  TNModel.D.5I <-lm(TN~1+TNscore.D.5I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5I')]))
  TNModel.DI.5 <-lm(TN~1+TNscore.DI.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DI.5')]))
  TNModel.D.5N <-lm(TN~1+TNscore.D.5N ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N','TNscore.I')]))
  TNModel.DN.5 <-lm(TN~1+TNscore.DN.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN.5','TNscore.I')]))
  TNModel.I.5N <-lm(TN~1+TNscore.I.5N ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I.5N')]))
  TNModel.IN.5 <-lm(TN~1+TNscore.IN.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.IN.5')]))
  TNModel.D.5NI <-lm(TN~1+TNscore.D.5NI ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5NI')]))
  TNModel.DN.5I <-lm(TN~1+TNscore.DN.5I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN.5I')]))
  TNModel.DNI.5 <-lm(TN~1+TNscore.DNI.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DNI.5')]))
  TNModel.D.5N.5I <-lm(TN~1+TNscore.D.5N.5I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N.5I')]))
  TNModel.DN.5I.5 <-lm(TN~1+TNscore.DN.5I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN.5I.5')]))
  TNModel.D.5NI.5 <-lm(TN~1+TNscore.D.5NI.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5NI.5')]))
  TNModel.D.5N.5I.5 <-lm(TN~1+TNscore.D.5N.5I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N.5I.5')]))
  TNModel.D.5I.5 <-lm(TN~1+TNscore.D.5I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5I.5')]))
  TNModel.D.5N.5 <-lm(TN~1+TNscore.D.5N.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N.5','TNscore.I')]))
  TNModel.I.5N.5 <-lm(TN~1+TNscore.I.5N.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I.5N.5')]))
  TNAICBIC <- cbind(AIC(TNModel,TNModel.D ,TNModel.I ,TNModel.N,TNModel.DI ,TNModel.DN ,TNModel.IN ,
                        TNModel.DNI ,TNModel.D.5 ,TNModel.I.5 ,TNModel.N.5,TNModel.D.5I ,TNModel.DI.5 ,
                        TNModel.D.5N ,TNModel.DN.5 ,TNModel.I.5N ,TNModel.IN.5 ,TNModel.D.5NI ,
                        TNModel.DN.5I ,TNModel.DNI.5 ,TNModel.D.5N.5I ,TNModel.DN.5I.5 ,TNModel.D.5NI.5 ,
                        TNModel.D.5N.5I.5 ,TNModel.D.5I.5 ,TNModel.D.5N.5 ,TNModel.I.5N.5 ,TNModel.null),
                    BIC(TNModel,TNModel.D ,TNModel.I ,TNModel.N,TNModel.DI ,TNModel.DN ,TNModel.IN ,
                        TNModel.DNI ,TNModel.D.5 ,TNModel.I.5 ,TNModel.N.5,TNModel.D.5I ,TNModel.DI.5 ,
                        TNModel.D.5N ,TNModel.DN.5 ,TNModel.I.5N ,TNModel.IN.5 ,TNModel.D.5NI ,
                        TNModel.DN.5I ,TNModel.DNI.5 ,TNModel.D.5N.5I ,TNModel.DN.5I.5 ,TNModel.D.5NI.5 ,
                        TNModel.D.5N.5I.5 ,TNModel.D.5I.5 ,TNModel.D.5N.5 ,TNModel.I.5N.5 ,TNModel.null))}


View(TNAICBIC)


#plot colors

region_colors <- c("Central Plains" = "orange", "Temperate Forest" = "forestgreen", "Southern Plains" = "red4")
season_colors <- c("Fall" = "brown4", "Spring" = "forestgreen", "Summer" = "yellow3")

??ggeffects::plot
#Nitrogen Model Plots
?geom_smooth


TPmodel.initial <- ggpredict(TPModel,terms = "TPscore [0:10]")
TPmodel.corrected <-ggpredict(TPModel.D.5,terms = "TPscore.D.5 [0:10]")
initial <- ggplot()+  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),
                                 aes(x=TPscore,y=TP, color="red3"),alpha = 0.8)+
  geom_ribbon(data=TPmodel.initial, aes(x=x,ymin = conf.low,ymax=conf.high),fill = "red3",alpha=0.3)+
  stat_smooth(data=TPmodel.initial, aes(x=x,y=predicted),
                                se=F,method = "glm",color = "red3")+ 
  xlab("Initial index score")+ylab("Log10-transformed total phosphorus (ug/L)")+

  # geom_smooth(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),method="glm",
  #             se=F,aes(x=TPscore,y=TP),size=0.1,linetype="dashed")+
  # geom_smooth(data=TPmodel.initial, se=F,aes(x=x,y=predicted),color="red3")+
  labs(title="")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
corrected<- ggplot()+  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),
                                  aes(x=TPscore.D.5,y=TP), color="green3",alpha = 0.8)+
  geom_ribbon(data=TPmodel.corrected ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "green3",alpha=0.3)+
  stat_smooth(data=TPmodel.corrected, aes(x=x,y=predicted),
              se=F,method = "glm",color = "green3")+ 
  xlab("Corrected index score")+ylab('')+

  xlim(0,10)+labs(title="")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()

gridExtra::grid.arrange(initial,corrected,ncol=2)


TNmodel.initial <- ggpredict(TNModel,terms = "TNscore [0:10]")
TNmodel.corrected <-ggpredict(TNModel.D,terms = "TNscore.D [0:10]")
initial <- ggplot()+  geom_point(data=na.omit(test[,c("TNscore","TNscore.D.5","TN")]),
                                 aes(x=TNscore,y=TN, color="red3"),alpha = 0.8)+
  geom_ribbon(data=TNmodel.initial, aes(x=x,ymin = conf.low,ymax=conf.high),fill = "red3",alpha=0.3)+
  stat_smooth(data=TNmodel.initial, aes(x=x,y=predicted),
              se=F,method = "glm",color = "red3")+ 
  xlab("Initial index score")+ylab("Log10-transformed total nitrogen (ug/L)")+

  # geom_smooth(data=na.omit(test[,c("TNscore","TNscore.D.5","TN")]),method="glm",
  #             se=F,aes(x=TNscore,y=TN),size=0.1,linetype="dashed")+
  # geom_smooth(data=TNmodel.initial, se=F,aes(x=x,y=predicted),color="red3")+
  labs(title="")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
corrected<- ggplot()+  geom_point(data=na.omit(test[,c("TNscore","TNscore.D.5","TN")]),
                                  aes(x=TNscore.D.5,y=TN), color="green3",alpha = 0.8)+
  geom_ribbon(data=TNmodel.corrected ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "green3",alpha=0.3)+
  stat_smooth(data=TNmodel.corrected, aes(x=x,y=predicted),
              se=F,method = "glm",color = "green3")+ 
  xlab("Corrected index score")+ylab('')+

  xlim(0,10)+labs(title="")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()

gridExtra::grid.arrange(initial,corrected,ncol=2)




initial2 <- plot(tnmodel.initial,line_size = 1,limit_range = F)+ 
  xlab("Initial index score")+ylab("Log10-transformed total nitrogen (mg/L)")+
  geom_point(data=na.omit(test[,c("TNscore","TNscore.D","TN")]),
             aes(x=TNscore,y=TN,group=region,color=season),alpha = 0.7)+
  scale_color_manual(values = season_colors)+xlim(0,10)+labs(title="",color="Season")+scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
corrected2<- plot(tnmodel.corrected,line_size = 1,limit_range = F)+ 
  xlab("Corrected index score")+ylab("")+
  geom_point(data=na.omit(test[,c("TNscore","TNscore.D","TN","region","season")]),
             aes(x=TNscore.D,y=TN,group=region,color=season),alpha = 0.7)+
  scale_color_manual(values = season_colors)+xlim(0,10)+labs(title="",color="Season")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
gridExtra::grid.arrange(initial,corrected,ncol=2)
gridExtra::grid.arrange(initial2,corrected2,ncol=2)


TPmodel.initial <- ggpredict(TPModel,terms = "TPscore")
TPmodel.corrected <-ggpredict(TPModel.D.5,terms = "TPscore.D.5")
initial <- plot(TPmodel.initial,line_size = 1,limit_range = T,colors = "red3")+ 
  xlab("Initial index score")+ylab("Log10-transformed total phosphorus (ug/L)")+
  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),
             aes(x=TPscore,y=TP, color="red3"),alpha = 0.7)+
  # geom_smooth(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),method="glm",
  #             se=F,aes(x=TPscore,y=TP),size=0.1,linetype="dashed")+
  # geom_smooth(data=TPmodel.initial, se=F,aes(x=x,y=predicted),color="red3")+
labs(title="")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
corrected<- plot(TPmodel.corrected,line_size = 1,limit_range = F, color="green3")+ 
  xlab("Corrected index score")+ylab('')+
  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),
             aes(x=TPscore.D.5,y=TP, color="green3"),alpha = 0.7)+
  xlim(0,10)+labs(title="")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
initial2 <- plot(TPmodel.initial,line_size = 1,limit_range = F)+ 
  xlab("Initial index score")+ylab("Log10-transformed total phosphorus (ug/L)")+
  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TP","region","season")]),
             aes(x=TPscore,y=TP,group=region,color=season),alpha = 0.7)+
  scale_color_manual(values = season_colors)+xlim(0,10)+labs(title="",color="Season")+scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
corrected2<- plot(TPmodel.corrected,line_size = 1,limit_range = F)+ 
  xlab("Corrected index score")+ylab("")+
  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TP","region","season")]),
             aes(x=TPscore.D.5,y=TP,group=region,color=season),alpha = 0.7)+
  scale_color_manual(values = season_colors)+xlim(0,10)+labs(title="",color="Season")+scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')
grid.arrange(initial,corrected,ncol=2)
grid.arrange(initial2,corrected2,ncol=2)

dataTN <- na.omit(test[,c('TNscore.DNI','Year','TN','TNscore','TNscore.D','TNscore.I',"TNscore.DN.5")])
dataTP <- na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5','TPscore.I',"TPscore")])

initialse <- NULL
initial2se <- NULL
correctse <- NULL
initialr <- NULL
initial2r <- NULL
correctr <- NULL
for (i in 1:1000) {
  #Creating a resampled dataset from the sample data
  sample_d = dataTN[sample(1:nrow(dataTN), nrow(dataTN), replace = TRUE), ]
  
  #Running the regression on these data
  initialTN <- summary(lm(TN~1+TNscore ,data=sample_d))
  initial2TN <- summary(lm(TN~1+TNscore.DN.5 ,data=sample_d))
  correctTN <-summary(lm(TN~1+TNscore.D ,data=sample_d))
  
  
  initial2se <- c(initial2se, initial2TN$coefficients[2,2])
  initialse <- c(initialse, initialTN$coefficients[2,2])
  correctse <- c(correctse,correctTN$coefficients[2,2])
  
  initialr <- c(initialr,initialTN$adj.r.squared)
  initial2r <- c(initial2r,initial2TN$adj.r.squared)
  correctr <- c(correctr, correctTN$adj.r.squared)

}

tncoefse <- cbind(initialse,correctse)
tncoefse2 <- pivot_longer(data.frame(tncoefse),cols=1:2, names_to = "Index",values_to = "Coefficent SE")

tnr <- cbind(initialr, correctr)
tnr2 <- pivot_longer(data.frame(tnr),cols=1:2, names_to = "Index",
                     values_to = "Rsq")

ggstatsplot::ggbetweenstats(data=tnr2 ,x=Index,y=Rsq,plot.type = "box", ylab = "Adjusted r-squared", 
                            results.subtitle = F, bf.message = F, 
                            centrality.plotting = F,xlab=c(" "),ggsignif.args = list(textsize =3, tip_length=0.001),
                            point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), 
                                              alpha = 0.4, size = 2, stroke = 0),
                            violin.args = list(fill = NA, color = "transparent"), boxplot.args= list(alpha=0),
)
t.test(Rsq~Index,data=tnr2)
ggstatsplot::ggbetweenstats(data=tncoefse2 ,x=Index,y=`Coefficent SE`,plot.type = "box", ylab = "SE of coefficient", 
                            results.subtitle = F, bf.message = F, 
                            centrality.plotting = F,xlab=c(" "),ggsignif.args = list(textsize =3, tip_length=0.001),
                            point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), 
                                              alpha = 0.4, size = 2, stroke = 0),
                            violin.args = list(fill = NA, color = "transparent"), boxplot.args= list(alpha=0),
)
t.test(`Coefficent SE`~Index,data=tncoefse2)
