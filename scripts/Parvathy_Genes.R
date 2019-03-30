
#Nicole E Soltis
#03/29/19


#--------------------------------------------------------------------------------------
#get SNPs +- 2 kb of these genes
#first, read in gene list
#then, from gtf get gene start and end positions
#then, from XXXX get SNPs for gene +- 2kb

rm(list=ls())
setwd("~/Projects/SideProjects/data")
mydata <- read.csv("Parvathy_GeneList.csv")
names(mydata)[1] <- "GeneID"

#from gtf, get gene start and end
#gtf location: 
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
my.gtf <- my.gtf[,c(1:8,10,12,14)]

library(plyr)
summ.gtf <- ddply(my.gtf, c("V10"), summarise,
                                   geneMin = min(V4),
                                   geneMax = max(V5))
mydata$V10 <- paste("transcript:",mydata$GeneID,sep="")
mydata2 <- merge(mydata, summ.gtf, by="V10")
setwd("~/Projects/SideProjects/data")
#write.csv(mydata2, "Parvathy_GeneRegions.csv")
mydat3 <- read.csv("Parvathy_GeneRegions_ed.csv")

#SNP data (binary, not bases) is at C:\Users\nesol\Documents\Projects\BcAt_RNAGWAS\data\GEMMA_eachAt_Bc\01_PLINK\OriginalSNPdata.csv

#try vcf
setwd("~/Projects/BcSolGWAS/data/BcGenome/WGS/big_set_v97iso_SNPs_filtered_qual30_dp6_maf20.recode.vcf")
library(vcfR)
vcf <- read.vcfR("big_set_v97iso_SNPs_filtered_qual30_dp6_maf20.vcf.recode.vcf", verbose = FALSE)
my.vcf <- vcf
my.vcf.df <- cbind(as.data.frame(getFIX(my.vcf)), INFO2df(my.vcf))
my.vcf.df$POS <- as.numeric(paste(my.vcf.df$POS))

setwd("~/Projects/SideProjects/data")
#1_g01260
my.vcf_1_g01260 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome1",]
my.vcf_1_g01260 <- my.vcf_1_g01260[my.vcf_1_g01260$POS > 508198,]
my.vcf_1_g01260 <- my.vcf_1_g01260[my.vcf_1_g01260$POS < 509296,]

#3_g04480
#NONE-- the max for Chr3 is 173095
my.vcf_3_g04480 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome3",]
max(my.vcf_3_g04480$POS)
my.vcf_3_g04480 <- my.vcf_3_g04480[my.vcf_3_g04480$POS > 1493741,]
my.vcf_3_g04480 <- my.vcf_3_g04480[my.vcf_3_g04480$POS < 1495275,]

#5g04960
my.vcf_5_g04960 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome5",]
max(my.vcf_5_g04960$POS)
my.vcf_5_g04960 <- my.vcf_5_g04960[my.vcf_5_g04960$POS > 1768689,]
my.vcf_5_g04960 <- my.vcf_5_g04960[my.vcf_5_g04960$POS < 1770143,]

#6g00026
my.vcf_6_g00026 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome6",]
max(my.vcf_6_g00026$POS)
my.vcf_6_g00026 <- my.vcf_6_g00026[my.vcf_6_g00026$POS > 22572,]
my.vcf_6_g00026 <- my.vcf_6_g00026[my.vcf_6_g00026$POS < 28113,]

#6g05230
my.vcf_6_g05230 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome6",]
max(my.vcf_6_g05230$POS)
my.vcf_6_g05230 <- my.vcf_6_g05230[my.vcf_6_g05230$POS > 1772070,]
my.vcf_6_g05230 <- my.vcf_6_g05230[my.vcf_6_g05230$POS < 1779326,]

#8g03460
#none, Chr8 ends at 66714
my.vcf_8_g03460 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome8",]
max(my.vcf_8_g03460$POS)
my.vcf_8_g03460 <- my.vcf_8_g03460[my.vcf_8_g03460$POS > 1327135,]
my.vcf_8_g03460 <- my.vcf_8_g03460[my.vcf_8_g03460$POS < 1332202,]

#9g01110
#none, Chr9 ends at 151640
my.vcf_9_g01110 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome9",]
max(my.vcf_9_g01110$POS)
my.vcf_9_g01110 <- my.vcf_9_g01110[my.vcf_9_g01110$POS > 407026,]
my.vcf_9_g01110 <- my.vcf_9_g01110[my.vcf_9_g01110$POS < 414107,]

#11g01310
my.vcf_11_g01310 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome11",]
max(my.vcf_11_g01310$POS)
my.vcf_11_g01310 <- my.vcf_11_g01310[my.vcf_11_g01310$POS > 447537,]
my.vcf_11_g01310 <- my.vcf_11_g01310[my.vcf_11_g01310$POS < 453290,]

#12g05690
#none, Chr12 ends at 937014
my.vcf_12_g05690 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome12",]
max(my.vcf_12_g05690$POS)
my.vcf_12_g05690 <- my.vcf_12_g05690[my.vcf_12_g05690$POS > 1962879,]
my.vcf_12_g05690 <- my.vcf_12_g05690[my.vcf_12_g05690$POS < 1968150,]

#12g06180
#none, Chr12 ends at 937014
my.vcf_12_g06180 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome12",]
max(my.vcf_12_g06180$POS)
my.vcf_12_g06180 <- my.vcf_12_g06180[my.vcf_12_g06180$POS > 2136114,]
my.vcf_12_g06180 <- my.vcf_12_g06180[my.vcf_12_g06180$POS < 2141607,]

#13g00710
my.vcf_13_g00710 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome13",]
max(my.vcf_13_g00710$POS)
my.vcf_13_g00710 <- my.vcf_13_g00710[my.vcf_13_g00710$POS > 274010,]
my.vcf_13_g00710 <- my.vcf_13_g00710[my.vcf_13_g00710$POS < 282827,]

#15g05080
#none, Chr15 ends at 114409
my.vcf_15_g05080 <- my.vcf.df[my.vcf.df$CHROM=="Chromosome15",]
max(my.vcf_15_g05080$POS)
my.vcf_15_g05080 <- my.vcf_15_g05080[my.vcf_15_g05080$POS > 1711546,]
my.vcf_15_g05080 <- my.vcf_15_g05080[my.vcf_15_g05080$POS < 1718383,]


setwd("~/Projects/SideProjects/data")
write.csv(my.vcf_1_g01260, "P_c1g01260_vcf.csv")
write.csv(my.vcf_11_g01310, "P_c11g01310_vcf.csv")
write.csv(my.vcf_13_g00710, "P_c13g00710_vcf.csv")
write.csv(my.vcf_3_g04480, "P_c3g04480_vcf.csv")
write.csv(my.vcf_5_g04960, "P_c5g04960_vcf.csv")
write.csv(my.vcf_6_g00026, "P_c6g00026_vcf.csv")
write.csv(my.vcf_6_g05230, "P_c6g05230_vcf.csv")