#
#  Copyright (C) 2016 Kliebenstein Lab(http://www.plantsciences.ucdavis.edu/kliebenstein/)
#  


wei.BotrytisANOVA <- function(){
  source("R/util.R")
  source("R/GlmNbAll.R")
 
  x <- read.csv("results/reads/botrytis.reads/selected.gene.counts.csv",row.names = 1)
  x <- data.frame(t(x))

  factors <- read.csv("keys/Select.SampleKey.csv",row.names = 1)
  factors <- factors[,c("Experiment","GrowingFlat","AgarFlat","Isolate","HostGenotype")]
  factors$Isolate <- as.factor(factors$Isolate)
  
  formula4glm <- "Experiment + Experiment/GrowingFlat + Experiment/GrowingFlat/AgarFlat + Isolate + HostGenotype + HostGenotype*Isolate"
  formula4lsm <- "~ Isolate | HostGenotype"
  
  dir.create(file.path("results/anova/botrytis.anova/"), showWarnings = FALSE, recursive = TRUE )
  
  octopus.glm.nb.all(x,factors,formula4glm,formula4lsm,dist = "results/anova/botrytis.anova/")
  
  octopus.info("adjust result.lsm remove duplicated column.")
  result.lsm <- read.csv("results/anova/botrytis.anova/result.lsm.csv",row.names = 1)
  if(length(grep("Isolate.",colnames(result.lsm))) >0 ) {
    result.lsm <- result.lsm[,-c(grep("Isolate.",colnames(result.lsm)),grep("HostGenotype.",colnames(result.lsm)))]
    write.csv(result.lsm,file= "results/anova/botrytis.anova/result.lsm.csv" )
  }
  octopus.info("result.lsm size [",dim(result.lsm)[1],",",dim(result.lsm)[2],"]")
  
  octopus.info("adjust result.se remove duplicated column.")
  result.se <- read.csv("results/anova/botrytis.anova/result.se.csv",row.names = 1)
  if(length(grep("Isolate.",colnames(result.se))) >0 ) {
    result.se <- result.se[,-c(grep("Isolate.",colnames(result.se)),grep("HostGenotype.",colnames(result.se)))]
    write.csv(result.se,file= "results/anova/botrytis.anova/result.se.csv" )
  }
  octopus.info("result.se size [",dim(result.se)[1],",",dim(result.se)[2],"]")
  
  #results <- read.csv("results/anova/botrytis.anova/results.csv",row.names = 1)
  
  octopus.info("BotrytisANOVA done!!!")
}

