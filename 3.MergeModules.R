source("../Mergeomics.R")

## below, i assume that you have your 1st msea results in the object "plan":
#merging modules
plan = c()
plan$folder = ""
plan$label = "" 
plan$modfile=""
plan$inffile=""
#=====================================================
## first pick the significantly enriched modules for 1st ssea:
FDR=0.05
pool=c()
file.name=paste0(plan$folder, "/",plan$label, ".txt")
if(file.exists(file.name)){
  aa=(read.table(file.name, header=T, sep='\t', check.names=F, quote=NULL)) 
  if(length(which(aa[,"MODULE"] == "_ctrlA")) >0) aa=aa[-which(aa[,"MODULE"] == "_ctrlA"),]
  if(length(which(aa[,"MODULE"] == "_ctrlB")) >0) aa=aa[-which(aa[,"MODULE"] == "_ctrlB"),]
  aa=aa[which( (as.numeric(aa[,"FDR"])< FDR)  ), ] 
  if (nrow(aa) > 0) pool=unique(c(pool, as.character(aa[,"MODULE"])))
}

#=====================================================
#=== Merge the modules before 2nd SSEA
if (length(pool)>0){
  meg.mods<- tool.read(plan$modfile)
  merged.modules <- pool
  moddata <- meg.mods[which(!is.na(match(meg.mods[,1], merged.modules))),]
  # Merge and trim overlapping modules.
  rmax <- 0.33
  moddata$OVERLAP <- moddata$MODULE
  moddata <- tool.coalesce(items=moddata$GENE, groups=moddata$MODULE, rcutoff=rmax) #, ncore=500)
  moddata$MODULE <- moddata$CLUSTER
  moddata$GENE <- moddata$ITEM
  moddata$OVERLAP <- moddata$GROUPS
  moddata <- moddata[,c("MODULE", "GENE", "OVERLAP")]
  moddata <- unique(moddata)
  
  moddatainfo <- tool.read(plan$inffile)
  moddatainfo <- moddatainfo[which(!is.na(match(moddatainfo[,1], moddata[,1]))), ]
  # Mark modules with overlaps.
  for(j in which(moddata$MODULE != moddata$OVERLAP)){
    moddatainfo[which(moddatainfo[,"MODULE"] == moddata[j,"MODULE"]), "MODULE"] <- paste(moddata[j,"MODULE"], "..", sep=",")
    moddata[j,"MODULE"] <- paste(moddata[j,"MODULE"], "..", sep=",")
  }
  # Save module info for 2nd SSEA and KDA.
  moddata <- unique(moddata)
  moddata[, 4] <- moddata[, 2];  names(moddata)[4] <- c("NODE")
  mdfile="modules.txt"; mifile="modules.info.txt" #change output names
  
  write.table(moddata, paste0(plan$folder,"/merged_", mdfile),  
              sep='\t', col.names=T, row.names=F, quote=F) 
  
  write.table(moddatainfo, paste0(plan$folder,"/merged_", mifile),  
              sep='\t', col.names=T, row.names=F, quote=F) 
  
  #==================================================================# 
  #    Apply 2nd SSEA with the merged files				             #
  #==================================================================#
  
}
