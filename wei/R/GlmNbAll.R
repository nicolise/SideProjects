#
#  Copyright (C) 2016 Kliebenstein Lab(http://www.plantsciences.ucdavis.edu/kliebenstein/) and Guocai Chen
#
#


##
#
# x       : counts matrix, expecting sample names to be rownames, gene names to be colnames
# faxtors : factor matrix, expecting sample names to be rownames, factor names to be colnames
# formula4glm : formula for glm.nb
# formula4lsm : formula for lsmeans
#
# pack glm.nb,anova,lsmeans in one founction
#
#
#
##
octopus.glm.nb.all <- function(x,factors,formula4glm,formula4lsm,dist="tmp/") {
  source("R/util.R")

  octopus.info("start run model with formula4glm -> [",formula4glm,"]")
  octopus.info("start run model with formula4lsm -> [",formula4lsm,"]")
  formula_glm <- formula( paste("x[,i] ~ ",formula4glm))
  formula_lsm <- formula( formula4lsm)

  octopus.library("MASS")
  octopus.library("car")
  octopus.library("lsmeans")

  result <- NULL
  result.lsm <- NULL
  result.se <- NULL
  skipped.genes <- NULL

  ## functions
  run_slice <- function(i){
    Gene.name <- colnames(x)[i]
    is_skip <- FALSE

    Gene.glm <- suppressWarnings(tryCatch({  glm.nb(formula_glm , data = x) }, error=function(e) {
      octopus.info(e)
      is_skip <<- TRUE
      "skip"
    }))

    if(is_skip) {  #Skip and record if a bad gene.
      octopus.info("SKIP! [",Gene.name,"] if you are getting a lots of this, better check you SampleKey file, probably caused by mismatching sample names, otherwise you are just fine :)")
      skipped.genes <<- c(skipped.genes, Gene.name)
      return(FALSE)
    }

    ##Pull p.values/variances/DFs from the model. Can't pull variances/deviances with Type II SS, so stuck with Type I
    Gene.aov <- suppressWarnings(tryCatch({ anova(Gene.glm) }, error=function(e) {
      octopus.info(e)
      is_skip <<- TRUE
      "skip"
    }))

    if(is_skip) {  #Skip and record if a bad gene.
      octopus.info("SKIP!  [",Gene.name,"]")
      skipped.genes <<- c(skipped.genes, Gene.name)
      return(FALSE)
    }

    Gene.res <- c(Gene.aov[,1], Gene.aov[,2], Gene.aov[,3], Gene.aov[,4], Gene.aov[,5])
    Gene.lsm <- summary(lsmeans(Gene.glm,  formula_lsm)) ##Calc LSMeans

    if(is.null(result)) {
      result <<- apply(expand.grid(rownames(Gene.aov), colnames(Gene.aov)),1,paste,collapse="_") #set up name column for p.values and variances to result
      result.lsm <<- Gene.lsm[1:(length(names(Gene.lsm))-5)] #set up factor column for lsmeans result
      result.se <<- Gene.lsm[1:(length(names(Gene.lsm))-5)]
    }

    #Get p.values and variances to result
    result <<- cbind(result,Gene.res)
    colnames(result)[length(colnames(result))] <<- Gene.name
    #Add lsmeans and se to result
    result.lsm <<- cbind(result.lsm,Gene.lsm$lsmean)
    colnames(result.lsm)[length(colnames(result.lsm))] <<- Gene.name
    result.se <<- cbind(result.se,Gene.lsm$SE)
    colnames(result.se)[length(colnames(result.se))] <<- Gene.name
    TRUE
  }

  run.normalize <- function(){
    #Normalize using edgeR
    octopus.library("edgeR")
    x <<- data.frame(t(x))
    x <<- DGEList(x) ###Transform to DGEobject
    octopus.info("TMM Normalization ...")
    x <<- calcNormFactors(x, method ="TMM")
    x <<- estimateCommonDisp(x, verbose=TRUE)

    octopus.info("write out TMM Normalization info to norm.counts.sample.info.csv")
    write.csv(x$samples,file = paste0(dist,"/norm.counts.sample.info.csv"))
    x <<- x$pseudo.counts
    x <<- data.frame(t(x))
    x <<- merge(factors, x, by = "row.names")
    colnames(x)[1] <<- "SampleName"
  }

  tmp_dist <- "tmp/glm.nb.all/"
  dir.create(file.path(tmp_dist), showWarnings = FALSE, recursive = TRUE )
  get_slice_file <- function(i){
    c(paste0(tmp_dist,i,"_result.csv"),
      paste0(tmp_dist,i,"_result.lsm.csv"),
      paste0(tmp_dist,i,"_result.se.csv"),
      paste0(tmp_dist,i,"_result.skipped.genes.csv")
    )
  }

  save_slice <- function(i){
    octopus.info("write out temporary results NO.",i)
    fs <- get_slice_file(i)
    write.csv(result, file = fs[1])
    write.csv(result.lsm,file=fs[2])
    write.csv(result.se,file=fs[3])
    write.csv(skipped.genes,file = fs[4])
    result <<- NULL
    result.lsm <<- NULL
    result.se <<- NULL
    skipped.genes <<- NULL
  }

  # make slices.
  slice_size <- 100
  slices <- NULL
  iStart = (dim(factors)[2] + 2)
  iEnd = dim(x)[2] + dim(factors)[2]
  for(i in iStart:iEnd) {  if(i%%slice_size == 0){ slices <- c(slices,i)} }
  if(iEnd%%slice_size != 0){ slices <- c(slices,iEnd)}

  merge_slices <- function(){
    octopus.info("merge results...")
    result <<- NULL
    result.lsm <<- NULL
    result.se <<- NULL
    skipped.genes <<- NULL

    for(i in slices) {
      sf <- get_slice_file(i)
      if(sum(file.exists(sf)) !=4 ){  octopus.info("missing file for slice NO.", i) ; return(0) }

      partial.result <- read.csv(sf[1] ,row.names = 1 )
      if(is.null(result)) {  result <<- partial.result
      }else{  result <<- cbind(result,partial.result[,2:dim(partial.result)[2]]) }

      partial.result.lsm <- read.csv(sf[2], row.names = 1)
      if(is.null(result.lsm)) { result.lsm <<- partial.result.lsm
      }else{ result.lsm <<- cbind(result.lsm,partial.result.lsm) }

      partial.result.se <- read.csv(sf[3] , row.names = 1)
      if(is.null(result.se)) { result.se <<- partial.result.se
      }else{ result.se <<- cbind(result.se,partial.result.se) }

      partial.skipped.genes <- tryCatch({
        read.csv(sf[4],row.names = 1  ) },
        error=function(e) { octopus.info("skip ",paste0(tmp_dist,i,"_result.skipped.genes.csv")) ; NULL
        })

      skipped.genes <<- rbind(skipped.genes,partial.skipped.genes)
    }
    octopus.info("result size [",dim(result)[1],",",dim(result)[2],"]")
    write.csv(result, file = paste0(dist,"results.csv"))
    write.csv(result.lsm,file= paste0(dist,"result.lsm.csv"))
    write.csv(result.se,file= paste0(dist,"result.se.csv"))
    octopus.info("skipped.genes size [",dim(skipped.genes)[1],"]")
    write.csv(skipped.genes,file = paste0(dist,"skipped.genes.csv"))
    octopus.info("merge results done!!!")
    return(0)
  }

  # try to reuse slices.
  pos <- iStart-1
  for(i in slices) {
    j <- sum(file.exists(get_slice_file(i)))
    if(j == 0) { next }
    if(j==4 && i-slice_size <= pos){  octopus.info("reusing file for slice NO.", i) ; pos <- i ; next }
    octopus.info("found tainted file for slice NO.", i )
  }
  if(pos == iEnd) { merge_slices() ; return(0) } # have got all slices

  # continue at pos+1
  run.normalize()
  sp <- octopus.speed.init(iEnd-pos+1)
  for(i in (pos+1):iEnd) {
    sp <- octopus.speed.tick(sp)
    run_slice(i)
    if(i%%slice_size == 0){ save_slice(i) }
  }
  if(iEnd%%slice_size != 0){ save_slice(iEnd) }
  merge_slices()
  return(0)
}




