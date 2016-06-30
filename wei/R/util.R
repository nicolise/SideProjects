#
#  Copyright (C) 2015 Kliebenstein Lab(http://www.plantsciences.ucdavis.edu/kliebenstein/) and Guocai Chen
#  
#


octopus.DEBUGGING <- TRUE

octopus.debug <- function(...) {
  if(!octopus.DEBUGGING) { return }
  print(paste("DEBUG" , Sys.time() , ": " ,  ... ) )
}

octopus.info <- function(...) {
  print(paste("INFO ",Sys.time() , ": " , ... ) )
}


octopus.library <- function(name){
  if(name %in% rownames(installed.packages()) == FALSE) {
    octopus.info("missing library:-->",name)
    octopus.info("Trying to install it from http://cran.rstudio.com/")
    install.packages(name, repos = "http://cran.rstudio.com/")
    # source("https://www.bioconductor.org/biocLite.R") 
    # biocLite("edgeR")
  }
  library(name,character.only = TRUE)
}

octopus.library("tools")

octopus.speed.init <- function( max_count = 0 ) {
  c(c=0,m=max_count,lb=Sys.time()-100,b=Sys.time())
}

octopus.speed.tick <- function( speed ,... ) {
  speed["c"] <- speed["c"] + 1
  if(Sys.time() -  speed["lb"] < 5){ return(speed)	}
  speed["lb"] <- Sys.time()
  calendar_time <- as.integer(Sys.time() - speed["b"])
  octopus.info(...)
  cps <- as.integer(speed["c"])/calendar_time
  cpn <- paste0(round(cps),"/s ")  
  cpm <- (60 * as.integer(speed["c"]))/calendar_time
  if(cps <1) { cpn <- paste0(cpn , round(cpm) ,"/m ") }
  cph <- (3600 * as.integer(speed["c"]))/calendar_time
  if(cps < 1) { cpn <- paste0(cpn , round(cph) ,"/h ") } 
  
  octopus.info("ETA:",speed["c"],"/",speed["m"],", speed:", cpn , " , calendar time:",calendar_time," seconds ...")
  speed
}

octopus.2Number <- function(x,names) {
  xx <- x[,names]
  xx <- apply(xx, c(1,2),function(i){
    if(is.na(i)){ 0 }else{ as.numeric(i) }
  })
  x[,names] <- xx
  x
}

##
# ram disk path.
# setup this to use a ramdisk as temporary storage.
# this could avoid lots of disk IO, speed things up and protect you SSD.
# Example for linux : 
#	  create ram disk in console : mount -t tmpfs -o size=2g tmpfs /mnt/ramdisk
#   setup octopus.ramdisk in R script : octopus.ramdisk.set("/mnt/ramdisk")
##
octopus.ramdisk <- NULL

##
# try use ramdisk
##
octopus.ramdisk.use <- function(f) {
  if(!is.null(octopus.ramdisk)){
    paste0(octopus.ramdisk,"/",f)
  }else{
    f
  }
}

###
# set octopus.ramdisk  path
###
octopus.ramdisk.set <- function(f) {
  octopus.ramdisk <<- f
  dir.create(file.path(f), showWarnings = FALSE, recursive = TRUE )
}
octopus.ramdisk.set("tmp")

###
# split file path 
###
octopus.dir_reverse <- function(x) {
  if (dirname(x)==x) {
    x
  }else{
    c(basename(x),octopus.dir_reverse(dirname(x)))
  } 
}

##
# return data frame
# group data by x[,i]  , apply data group on func
# data group could be list(1 row) or frame(multiple rows)
#
#
# sample: 
# func <- function(x,gvalue) {
#   if(is.null(dim(x))) return(x) # one row
#   ret <- apply(x,c(1,2),as.numeric) 
#   ret <- apply(ret, 2, sum, na.rm = FALSE) # do sum on column
#   c(gvalue,ret[-1]) # re-attach name 
# }
#
# sort x by x[,i] before call this function
# x <- x[with(x, order(colname)), ] # sort 
# a <- gapply(x ,1,func)
##

octopus.gapply <- function(x,i,func){
  
  rows <- NULL
  gvalue <- ""
  sp <- octopus.speed.init(dim(x)[1])
  
  octopus.info("octopus.gapply start!!!")
  group <- function(xx){
    gvalue_new <- xx[i]
    sp <<- octopus.speed.tick(sp,gvalue_new)
    ret <- NULL
    if(gvalue_new == gvalue ){# row with same gvalue
      rows <<- rbind(rows,xx)
    }else{#row with new gvalue
      if(!is.null(rows)) ret <- func(rows,gvalue)
      
      #if(!is.null(dim(rows))) { # one row
      #  ret <- apply(rows,c(1,2),as.numeric) 
      #  ret <- apply(ret, 2, sum, na.rm = FALSE)
      # ret <- c(gvalue,ret[-1]) # re-attach name 
      #}
      gvalue <<- gvalue_new
      rows <<- xx
    }
    return(ret)
  }
  
  a <- apply(x,1,group)
  a <- a[!sapply(a,is.null)] # remove NULL
  a <- t(data.frame(a))
  if(!is.null(rows)) { #don't forget the last line
    a <- rbind(a,func(rows,gvalue)) # join the last line
  }
  data.frame(a)  # data.frame(a,row.names = 1)
}

##
# bowtie build
#
##
octopus.bowtie_build <- function(src_dir="GenomeSeq/cds/",dist_dir="GenomeSeq/index/") {
  octopus.library("tools")
	files <- list.files(src_dir, pattern="", recursive = TRUE , all.files = FALSE)
	for(f in files){
		inf = paste0(src_dir , "/" ,f)
		outf = paste0(dist_dir , "/" , file_path_sans_ext(f))
		dir.create(outf, showWarnings = FALSE, recursive = TRUE )
		outf = paste0(outf , "/" , file_path_sans_ext(f))
		cmd_line <- paste0("bowtie-build " , inf , " " , outf )
		octopus.debug(cmd_line)
		system(cmd_line)
	}
}


