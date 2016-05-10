#Preliminary data analysis Eudicot spp. x Botrytis cinerea lesions
#011916
#-------------------------------------------------------------
#clear memory
rm(list=ls())
#set working directory
setwd("S:/FacultyData/KLIEBENSTEIN/KLIEBENSTEINShared/Gongjun Shi/Inoculation/CiAnalysis")
#read-in your data files
ModDat <- read.csv("CiModelData.csv")
attach(ModDat)
#--------------------------------------------------------
#STATISTICS!
#check data structure
xtabs(~ Rep + ImageName, ModDat)
xtabs(~ PlantID + PnumID, ModDat)
#ImageName is nested within Rep
#and both are random
#PnumID is nested within PlantID
#and both are random
#And we can consider PlantID to be nested within Domestication
#IsolateID, Domestication, PlantID, AorB are fixed(What's A and B)

#Variance output of summary(Mod) gives you SS for the random factors
#rand(Mod) gives Chi-sq  and P values for random factors in packages lmerTest
install.packages("MASS")
library(MASS)
library(lmerTest)

library(lme4)


#first check normality of Scale.LS
install.packages("car")
library(car)

require(car)


#Data transformation
ModDat$Scale.LS.t <- ModDat$Scale.LS + 1
#fairly normal
qqp(ModDat$Scale.LS.t, "norm")

#definitely not log-normal
qqp(ModDat$Scale.LS.t, "lnorm")
#data must be integers for the following...
ModDat$Scale.LS.i <- ModDat$Scale.LS*100 + 100
ModDat$Scale.LS.i <- round(ModDat$Scale.LS.i)
#negative binomial
nbinom <- fitdistr(ModDat$Scale.LS.i, "Negative Binomial")
qqp(ModDat$Scale.LS.i, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
#poisson
poisson <- fitdistr(ModDat$Scale.LS.i, "Poisson")
qqp(ModDat$Scale.LS.i, "pois", poisson$estimate)

#------------------------------------------------------------------------------------------------
#full model
#nesting terms are already included- don't need to add Species as a separate term BUT for random effects (ExpBlock alone) do need a separate term

#only 2 ExpBlocks so does it make more sense to include them as fixed effects?
#reasonable to consider AgFlat as a random effect -- have 16 x 3 per exp
#but maybe include a term for "bench"?? random or fixed?
#ExpBlock/Bench/AgFlat

#PlGenoNm is a term nested within Species (but not CODED as if nested within Species = it's not an implicitly nested factor)
#AgFlat IS implicitly nested within ExpBlock -- let's fix this
#non-numeric factors to use: Igeno, PlGenoNm, Species, ExpBlock, AgFlat

#optional to fix: coding of AgFlat so that it is not implicitly nested

#the lme4 mixed model:
#WARNING: this can take ~hours to run (Nicole)
fullmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)

#Raoni's model, model B with interactions
fullmodA <- lmer(Scale.LS ~ IsolateID + Domestication + Domestication/PlantID + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/ImageName) + (1|Exp/Rep/PnumID) , data = ModDat)

fullmodB <- lmer(Scale.LS ~ IsolateID + Domestication + Domestication/PlantID + IsolateID:Domestication + IsolateID:Domestication/PlantID + (1|Exp) + (1|Exp/Rep) + (1|Exp/Rep/ImageName) + (1|Exp/Rep/PnumID) , data = ModDat)

anova(fullmodA)
rand(fullmodA)

anova(fullmodB)
rand(fullmodB)

anova(fullmodA, fullmodB)
#end of Raoni's script

#alternate: writing the whole thing out to better understand nesting:
#Does it make sense to include Igeno:PlGenoNm separately? I think not.
Sys.time()
fullmodALT <- lmer(Scale.LS ~ Igeno + Species + Species:PlGenoNm + Igeno:Species + Igeno:Species:PlGenoNm + (1|ExpBlock) + (1|ExpBlock:AgFlat) + (1|IndPlant) + AorB , data = ModDat)
save.image(file = "model.RData") #model is saved as .RData
#use ?load to reload it

#print system time to keep track fo how long processes take
Sys.time()

#copy output to a text file
#other model options- glmm.object? glmm package? nlme?
sink(file='output.txt')
Sys.time()
summary(fullmodALT) # the code generating output
Sys.time()
rand(fullmodALT)
Anova(fullmodALT, type=2)
anova(fullmodALT)
Sys.time()
sink()

library(sjPlot)

#can also minimize model to avoid overfitting

sjp.glmer(fullmod, type = "fe.cor")
sjp.glmer(fullmod, type = "re.qq")

#also look at whether experiment ALTERS the effect of isolate or plant genotype:
#fullmod including Isolate*Exp and Plant*Exp
#add a factor by isolate * experiment
names(ModDat)
ModDat$IbyX <- paste(ModDat$Igeno, ModDat$ExpBlock, sep='') 
sort(unique(ModDat$IbyX))
#remove isolates with no match across experiments:
#94.1, Gallo3
ModDat2 <- ModDat[ModDat$Igeno!="94.1",]
ModDat2 <- ModDat2[ModDat2$Igeno!="Gallo3",]

library(lme4); library(lmerTest)
library(car)
fullmod5 <- lmer(Scale.LS ~ Igeno + Species/Pgeno + Igeno:Species/Pgeno + Species:Igeno + (1|ExpBlock/AgFlat) + (1|Plant) + AorB + (1|ExpBlock) + Igeno:ExpBlock + Species:ExpBlock + ExpBlock:Species/Pgeno + ExpBlock:Igeno:Species/Pgeno, data = ModDat2)
#copy output to a text file
sink(file='output.txt')
summary(fullmod5) # the code generating output
anova(fullmod5) # the code generating output
Anova(fullmod5, type=2)
rand(fullmod5)
sink()

library(car)
qqnorm(residuals(fullmod))
qqline(residuals(fullmod))
plot(fullmod)
print(fullmod)


aovmod1 <- aov(lm(Lesion.Size ~ Pexp + Pexp/PImage + PPlant*Isolate, data=SrtDat3))
summary(aovmod1)

#an alternate model
mod1 <- lm(Lesion.Size ~ PPlant*Isolate + PPlant/PInPlant/PInLeaf + PInLflt + Pexp , data=SrtDat3)
#lme for lsmeans
library(lme4)
mod1 <- lmer(Lesion.Size ~ PPlant*Isolate + (1|PImage) + (1|Pexp), data=SrtDat3)
summary(mod1) 
#lsmeans for GWAS
#doesn't work yet
library(lsmeans)
#myLSmeans <- lsmeans(mod1, Lesion.Size ~ PPlant|Isolate, adjust="none")