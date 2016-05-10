rm(list=ls())
setwd("~/Projects/SideProjects/data/labInterns")
myDataX <- read.csv("BR_selectMedia_0407.csv")
myData <- read.csv("BR_Hyphae_imagedata.csv")

names(myData)

#add a new variable for our digital estimate of waviness
myData$DigWave <- myData$TimesCrossing / myData$Length
#standardize it
myData$DigWave.s <- (myData$DigWave - mean(myData$DigWave))/sd(myData$DigWave)
myData$DigWave.sp <- myData$DigWave.s+2

library(plyr)
attach(myData)
myFigDat <- myData[complete.cases(myData),]
myFigDat <- ddply(myFigDat, "Isolate", summarise, #na.rm=T
                  N    = sum(!is.na(DigWave.sp)) ,
                  mean = mean(DigWave.sp, na.rm=T),
                  sd   = sd(DigWave.sp, na.rm=T),
                  se   = sd / sqrt(N))

library(ggplot2)
limits <- aes(ymax = mean + se, ymin=mean - se)
ggplot(myFigDat, aes(x = factor(Isolate), y = mean))+
  geom_bar(stat="identity", fill="dodgerblue3")+
  theme_bw()+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))+
labs(y=expression(Mean ~ Wavy ~ Index), x=element_blank())+
  geom_errorbar(limits, width=0.25)
#from here is for plots with media concentration only
+
  facet_grid(.~Date, scales="free")+ 
  scale_y_continuous(limits = c(0,4.2)) 


#compare to manual scale
IsolateNm <- c("DavisNavel", "Gallo1", "Grape", "Pepper")
DigMean <- c(2.15359, 1.786423, 2.0566, 1.972669)
ManMean <- c(7.45, 6.33333, 7.72, 6.55)
dnew <- data.frame( Isolate=IsolateNm , ManuMean=ManMean , DigiMean=DigMean )
plot(dnew$ManuMean, dnew$DigiMean)
m0 <- lm( DigiMean ~ ManuMean, dnew )
abline( m0  )

1 - var(resid(m0))/var(dnew$DigiMean)
