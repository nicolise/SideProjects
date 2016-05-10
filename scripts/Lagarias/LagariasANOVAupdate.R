#Nicole E Soltis
#03/17/16
#ANOVA for Lagarias RNAseq

#-------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/SideProjects/data/Lagarias")
myData <- read.csv("for_pathwayANOVAed.csv")
names(myData)

#is data normal?
attach(myData)
#graphically...
hist(RPKM)
#more graphs
require(car); require(MASS)
myData$RPKM.t <- myData$RPKM + 1
#is it more normal or log-normal?
#very long right tail, log-normal is also a bad estimate
qqp(myData$RPKM.t, "norm")
qqp(myData$RPKM.t, "lnorm")
#statistically...
#dataset is too large (37980 obs. instead of 5000 max)
#shapiro.test(myData$RPKM.t) 
#try transformations
transf <- log10(RPKM)
hist(transf)
#randomly select a subset of values from transf to test
set.seed(100)
sample.shapiro <- sample( 1:nrow(myData) , size=1e3 , replace=TRUE )
sample.RPKM <- myData$RPKM[ sample.shapiro ]
df.shapiro <- data.frame(matrix(unlist(sample.RPKM), nrow=1e3, byrow=T))
df.shapiro$RPKM <- df.shapiro$matrix.unlist.sample.RPKM...nrow...1000..byrow...T.
df.shapiro$transf <- (log(df.shapiro$RPKM+1))
hist(df.shapiro$transf)
shapiro.test(df.shapiro$transf) #still significantly non-normal
qqp(df.shapiro$transf) #but it looks pretty good

#at any rate I'll log-transform
myData$RPKM.t <- (log(myData$RPKM+1))

#next check assumption of homoscedasticity
#graphically...
attach(myData)
boxplot(RPKM.t~GenotypeID*Time,
        ylab="YTITLE", main="PLOTTITLE", las=3) #looks good yay!
#statistically...
bartlett.test(RPKM.t ~ GenotypeID, data=myData) #meh not homoscedastic
bartlett.test(RPKM.t ~ Time, data=myData) #also iffy
leveneTest(RPKM.t~GenotypeID) #also not homoscedastic

#go ahead with ANOVA anyway
head(myData)
#30 pathways, 979 genes
as.data.frame(table(unique(myData[])$GeneID)) 

#small model
Model.lm <- lm(RPKM.t ~ PathwayID + GenotypeID + Time)
#big model with RPKM.t
#gives error with / or : or %in%
#gives error with aov
Model.lm <- aov(RPKM.t ~ PathwayID/GeneIDV55)

#maybe I can subset to find out what component in the big model is making it so large
#currently 3.8 x 1e4 observations
#randomly select a subset of values to test
set.seed(100)
sample.mini <- sample( 1:nrow(myData) , size=1e3 , replace=TRUE )
df.mini <- myData[sample.mini,]
minimod.lm <- lm(RPKM.t ~ PathwayID/GeneIDV55, data=df.mini)

#try running within each pathway
out <- split( myData , f = myData$PathwayID )
head(out[[1]]) #30 elements

#Using a for loop, iterate over the list of data frames in out[[]]
sink(file='ModelsBYpathway_032316.txt')
#skip 13: blank, 59: 94.1, 68: blank, 77: Gallo3, 99: blank
for (i in c(1:30)) {
  print(unique(out[[i]]$PathwayID))
  Mod <- lm(RPKM.t ~ GeneIDV55 + GenotypeID + Time + GeneIDV55:GenotypeID + GeneIDV55:Time + GenotypeID:Time, data=out[[i]])
  result <- anova(Mod)
  print(result)
}
sink()


# more analysis
MY.ANOVA <- anova(Model.lm)
summary(MY.ANOVA)
MY.ANOVA #sig fx of genotype and time
interaction.plot(GenotypeID,Time,RPKM.t)
TukeyHSD(MY.ANOVA)
MY.aov <- aov(RPKM.t~GenotypeID*Time)
summary(MY.aov)

#adjusted p-values give 6 sig. pairs and 1 marg. sig.
TukeyHSD(MY.aov)