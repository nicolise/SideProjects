 #Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
#clear memory
rm(list=ls())
#set working directory
setwd("C:/Users/Raoni/Desktop/Example")
getwd()
#read-in your data files
ModDat <- read.csv("ModelDataDT2.csv")
attach(ModDat)
ModDat$SpLabs <- factor(ModDat$origin, labels = c("Domesticated", "Plant Introduction"))
#-------------------------------------------------
#plot my data: scatterplot!
names(ModDat)
attach(ModDat)
head(ModDat)
names(ModDat)
#if you do not have these packages, first install them

library(beanplot); library(ggplot2); library(RColorBrewer); library(plyr)
library(dplyr)
#add PlantNum as an integer sorted by mean lesion size
###NES added mean(Scale.ls, na.rm=T)
FigDat3 <- ddply(ModDat, c("genotypename", "isolatename", "SpLabs"), summarise,
                 mLS   = mean(Scale.LS, na.rm=T))
head(FigDat3)
attach(FigDat3)
names(FigDat3)

#add a column of mmLS (mean of mean lesion size) per isolate
#sort dataframe by mmLS 
#then color by the new factor mmLS
###NES added ave(mLS, isolate, na.rm=T)
FigDat3$mmLS <- ave(FigDat3$mLS, FigDat3$isolatename, na.rm=T)

attach(FigDat3)
FigDat3 <- FigDat3[order(FigDat3$mmLS),]
ggplot(FigDat3, aes(x = genotypename, y = mLS))+
  geom_point()+
  theme_bw()+
 # geom_line(size=1, aes(color=isolatename, group=factor(isolatename)), show.legend=F)+
 geom_line(size=1, aes(color=factor(mmLS), group=factor(isolatename)), show_guide=F)+
   theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())


#-------------------------------------------------------
#barplot with average lesion size by genotype
#and SE bars
FigDat <- ModDat
names(FigDat)
###NES added na.rm=T
FigDat2 <- ddply(FigDat, "genotypename", summarise, mLS = mean(Scale.LS, na.rm=T))

library(sciplot)
FigDat3 <- ddply(FigDat, "genotypename", summarise, seLS = se(Scale.LS, na.rm=T))
head(FigDat2)

###NES modified to deal with NA values
FigDat2 <- ddply(FigDat, c("genotypename", "SpLabs"), summarise,
               N    = sum(!is.na(Scale.LS)),
               mean = mean(Scale.LS, na.rm=T),
               sd   = sd(Scale.LS, na.rm=T),
               se   = sd / sqrt(N))
limits <- aes(ymax = mean + se, ymin=mean - se)
ggplot(FigDat2, aes(x = factor(genotypename), y = mean))+
  geom_bar(stat="identity", fill="dodgerblue3")+
  theme_bw()+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y=expression(Mean ~ Lesion ~ Area ~ (cm^{2})), x=element_blank())+
  geom_errorbar(limits, width=0.25)+
  facet_grid(.~SpLabs, scales="free")
names(ModDat)
#---------------------------------------------------------
#WE STOPPED HERE - NICOLE
#scatterplot: mean of each isolate on each host // domestication level
ModDat$SbyI <- paste(ModDat$genotypename, ModDat$isolatename, sep='') 
FigDat4 <- ddply(ModDat, c("genotypename", "isolatename","origin", "SbyI","PbyI"), summarise,PlMean   = mean(Scale.LS))
FDmeans <- ddply(ModDat, c("SbyI"), summarise, SpMean=mean(Scale.LS))
FigDat4 <- merge(FigDat4, FDmeans, by="SbyI")
head(FigDat4)
attach(FigDat4)
ggplot(FigDat4, aes(x = PlMean, y = SpMean, color=origin, group=genotypename))+
  geom_point()+
  geom_line()+
  theme_bw()+
  geom_line(size=1, aes(color=factor(isolatename)), show_guide=F)+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x="Plant Genotype")+
  theme(axis.text.x = element_blank())

#or more simply...
# mean lesion size per domestication status
names(FigDat4)
FigDat5 <- FigDat4[!duplicated(FigDat4$PbyI), ]
#color lines by positive or negative
unique(unlist(FigDat5$origin))
FigDat5w <- filter(FigDat5, origin == "Wl") 
FigDat5d <- filter(FigDat5, origin == "Dm") 
FigDat5d <-FigDat5d[,c(1,3,4,7)]
FigDat5w <-FigDat5w[,c(1,3,4,7)]
unique(FigDat5w$isolatename)
FigDat5d <- FigDat5d[!duplicated(FigDat5d$isolatename), ]
FigDat5w <- FigDat5w[!duplicated(FigDat5w$isolatename), ]
FigDat6 <- merge(FigDat5w, FigDat5d, by="isolatename")
#get mean difference (Dm - Wl)
head(FigDat6)
FigDat6 <- transform(FigDat6, MeanDiff=(SpMean.y - SpMean.x))
FigDat6$Direction <- ifelse(FigDat6$MeanDiff > 0,1,0)
FigDat7 <- FigDat6[,c(1,9)]
FigDat5 <- merge(FigDat5, FigDat7, by="isolatename")

names(FigDat5)
ggplot(FigDat5, aes(x = origin, y = SpMean, group=isolatename))+
  geom_point()+
  geom_line(size=1, aes(color=factor(Direction)), show.legend = F)+
  theme_bw()+
  labs(y=expression(Mean ~ Lesion ~ Area ~ per ~ Isolate ~ (cm^{2})), x=element_blank())+
  scale_x_discrete(labels=c("Domesticated","Wild"))

#---------------------------------------------------------------------------------------------

library(ggplot2)
ggplot (data = ModDat, 
  aes(x=GenFx, y=bigH))+
  geom_violin(adjust = 0.5, scale = "width", fill="#E6F598")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  stat_summary(fun.y="median", geom="point")+
  labs(y=expression(Heritability~(H^2)), x=NULL)+
  ylim(c(0,NA))+ #because anything less than 0 is uninteresting
  #and just means that variance for some genos > total variance
  scale_fill_manual(
    values = c("#E6F598","#E6F598", "#E6F598"))+
#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)
geom_jitter(shape=16, position=position_jitter(0.2))


library(vioplot)
#use complete cases only
MyPlotcc <- MyPlot[complete.cases(MyPlot),]
#set H<0 to 0
MyPlotcc$bigH[MyPlotcc$bigH<0] <- 0
x1 <- MyPlotcc$bigH[MyPlotcc$GenFx=="Isolate"]
x2 <- MyPlotcc$bigH[MyPlotcc$GenFx=="Plant"]
x3 <- MyPlotcc$bigH[MyPlotcc$GenFx=="PbyI"]
vioplot(x1, x2, x3, names=c("Isolate", "Plant", "PbyI"), 
        col="green")