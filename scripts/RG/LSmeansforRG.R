#LSmeans for RG
#Nicole E Soltis
#05/02/16
#-----------------------------------
rm(list=ls())
setwd("~/mydir")
#read in file 
ModDat <- read.csv("SlBcDATAFRAME.csv")

#----------------------------------------------------------------------
library(lme4); library(car); library(lmerTest)
#lsmeans calculations
#run model per isolate WITHIN each plant genotype
#so include no species terms or plant genotype terms
attach(ModDat)
out <- split( ModDat , f = ModDat$Pgeno)
head(out[[1]]) 

#Using a for loop, iterate over the list of data frames in out[[]]
sink(file="LSMeans032116.txt")
for (i in c(1:12)) {
  print(unique(out[[i]]$Pgeno))
  Lesion.lm <- lmer(Scale.LS ~ Igeno + Plant/Leaf/AorB + (1|Exp), data=out[[i]])
  Lesion.lsm <- lsmeans(Lesion.lm, "Igeno")
  print(Lesion.lsm)
}
sink()
