source("../Mergeomics.R")

setwd("~/Desktop/Mergeomics/")

job.kda <- list()
job.kda$label<-"wKDA"
job.kda$folder<-"Results" ## parent folder for results
## Input a network
## columns: TAIL HEAD WEIGHT
job.kda$netfile<-"bayesian_network.txt"

## Gene sets derived from ModuleMerge, containing two columns, MODULE, 
## NODE, delimited by tab 
job.kda$modfile<-"merged_modules.txt"

## "0" means we do not consider edge weights while 1 is opposite.
job.kda$edgefactor<-1
## The searching depth for the KDA
job.kda$depth<-1
## 0 means we do not consider the directions of the regulatory interactions
## while 1 is opposite.
job.kda$direction <- 1
job.kda$nperm <- 2000

moddata <- tool.read(job.kda$modfile)
mod.names <- unique(moddata$MODULE)
moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))),]
## save this to a temporary file and set its path as new job.kda$modfile:
tool.save(moddata, "subsetof.supersets.txt")
job.kda$modfile <- "subsetof.supersets.txt"

## Let's run KDA!
job.kda <- kda.configure(job.kda)
job.kda <- kda.start(job.kda)
job.kda <- kda.prepare(job.kda)
job.kda <- kda.analyze(job.kda)
job.kda <- kda.finish(job.kda)

job.kda <- kda2cytoscape(job.kda)
