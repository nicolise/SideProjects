rm(list=ls())
setwd("~/Projects/SideProjects/data/labInterns")
setwd("~/Documents/GitRepos/SideProjects/data/LabIntern")
myDatam <- read.csv("BR_selectMedia_0509.csv")

names(myDatam)
newData <- myDatam[(myDatam$PDAConc)=="0.75",]
#on March 07 we started the 0:10 numbering (previously ~0:5)
newData <- newData[newData$Date %in% c("4-Apr","7-Mar","28-Mar","21-Mar", "28-Apr"), ]

#data for PDA conc comparison
mediaData <- myDatam[myDatam$Date %in% c("1-Feb","22-Feb","7-Mar"), ]

names(mediaData)
mediaData$IsoBYpda <- paste(mediaData$Isolate, mediaData$PDAConc, sep='.') 
attach(mediaData)
Media.lm <- lm(Count ~ Isolate * PDAConc, data=mediaData)
Media.lm <- lm(Count ~ Isolate * PDAConc, data=newData)
Media.lm2 <- lm(Count ~ IsoBYpda + Date, data=myDatam)
anova(Media.lm)
anova(Media.lm2)
pairwise.t.test(Count, PDAConc, p.adj="none")

#remove rows with NA
myFigDatm <- newData[complete.cases(newData),]

mediaFig <- ddply(mediaData, c("Isolate", "PDAConc"), summarise,
                 N    = sum(!is.na(Count)),
                 mean = mean(Count, na.rm=T),
                 sd   = sd(Count, na.rm=T),
                 se   = sd / sqrt(N))

myFigDat2 <- ddply(myFigDatm, c("Isolate"), summarise,
                  N    = sum(!is.na(Count)) ,
                  mean = mean(Count, na.rm=T),
                  sd   = sd(Count, na.rm=T),
                  se   = sd / sqrt(N))

#reorder bars by mean waviness
myFigDat2 <- transform(myFigDat2, 
                          Isolate = reorder(Isolate, mean))

limits <- aes(ymax = mean + se, ymin=mean - se)
ggplot(myFigDat2, aes(x = factor(Isolate), y = mean))+
  geom_bar(stat="identity", fill="dodgerblue3")+
  theme_bw()+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))+
labs(y=expression(Mean ~ Wavy ~ Index), x=element_blank())+
  geom_errorbar(limits, width=0.25)


#from here is for plots with media concentration only
mediaFig <- transform(mediaFig, 
                       Isolate = reorder(Isolate, mean))

medialimits <- aes(ymax = mean + se, ymin=mean - se)
ggplot(mediaFig, aes(x = factor(Isolate), y = mean))+
  geom_bar(stat="identity", fill="dodgerblue3")+
  theme_bw()+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y=expression(Mean ~ Wavy ~ Index), x=element_blank())+
  geom_errorbar(medialimits, width=0.25)+
  facet_grid(.~PDAConc, scales="free")+ 
  scale_y_continuous(limits = c(0,9)) 