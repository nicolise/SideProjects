rm(list=ls())
setwd("~/Projects/SideProjects/data/labInterns")
myDatam <- read.csv("BR_selectMedia_0509.csv")

names(myDatam)
newData <- myDatam[(myDatam$PDAConc)=="0.75",]
#on March 07 we started the 0:10 numbering (previously ~0:5)
newData <- newData[newData$Date %in% c("4-Apr","7-Mar","28-Mar","21-Mar", "28-Apr"), ]

names(myData)
myData$IsoBYpda <- paste(myData$Isolate, myData$PDAConc, sep='.') 
attach(myData)
Media.lm <- lm(Count ~ Isolate * PDAConc, data=myData)
Media.lm <- lm(Count ~ Isolate * PDAConc, data=newData)
Media.lm2 <- lm(Count ~ IsoBYpda + Date, data=myData)
anova(Media.lm)
anova(Media.lm2)
pairwise.t.test(Count, PDAConc, p.adj="none")

#remove rows with NA
myFigDatm <- newData[complete.cases(newData),]

myFigDat <- ddply(myFigDat, c("Isolate", "Date"), na.rm=T, summarise,
                 N    = length(Count),
                 mean = mean(Count),
                 sd   = sd(Count),
                 se   = sd / sqrt(N))

myFigDat2 <- ddply(myFigDatm, c("Isolate"), summarise,
                  N    = sum(!is.na(Count)) ,
                  mean = mean(Count, na.rm=T),
                  sd   = sd(Count, na.rm=T),
                  se   = sd / sqrt(N))

limits <- aes(ymax = mean + se, ymin=mean - se)
ggplot(myFigDat2, aes(x = factor(Isolate), y = mean))+
  geom_bar(stat="identity", fill="dodgerblue3")+
  theme_bw()+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))+
labs(y=expression(Mean ~ Wavy ~ Index), x=element_blank())+
  geom_errorbar(limits, width=0.25)
#from here is for plots with media concentration only
+
  facet_grid(.~Date, scales="free")+ 
  scale_y_continuous(limits = c(0,4.2)) 