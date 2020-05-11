############################################################
# Example Meta-MSEA of two GWAS studies                    #
############################################################

# set working directory to where you have your results files
source("Mergeomics.R")

# The below set up is assuming you have your results files set up as:
# trait.mapping_method.results.txt
# Example:
# DIAGRAMstage1_T2D.Artery_Aorta_sQTL.results.txt

GWAS_studies = c("_gwas1197","_gwas1931")
# eqTL methods
# tissues = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta","Artery_Tibial",
#             "Brain_Cerebellar_Hemisphere","Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", 
#             "Brain_Hypothalamus", "Colon_Sigmoid", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Left_Ventricle", 
#             "Liver", "Muscle_Skeletal","Nerve_Tibial", "Pancreas", "Pituitary", "Spleen", "Stomach", "Thyroid", "Whole_Blood")

#Step1. Compile result files into one

# outside loop is the item that is not being meta-analyzed
# for(g in GWAS_studies){ # tissues is eQTL methods
  res=data.frame()
  for(study in GWAS_studies){
  name=paste("Leprosy_MSEA_pathways",study,".results.txt",sep="") 
    if (file.exists(paste0("./Results/msea/",name))){ 
      print("yes");
      temp=read.delim(name,header=T,as.is=T)[,c(1,2,7)]  #columns 1,2,7 are MODULE,p-value,FDR
      print(temp);
      if (dim(temp)[1]>0){
         temp$MAP=study;temp$Platform="Leprosy";res=rbind(res,temp)
      }
    }
  }

  res=res[!grepl("_ctrlA",res$MODULE)&!grepl("_ctrlB",res$MODULE),] #Remove internal controls
  write.table(res,paste("./combined/combined.", t, ".txt", sep = ""),row.names=F,quote=F,sep="\t") 
  
  #Step2a. Meta-MSEA
  options(stringsAsFactors = FALSE)
  res=read.delim(paste("./combined/combined.", t, ".txt", sep = ""))
  fullRes = res
  fullData = data.frame()
  count = 0
  addFlag = TRUE
  for(mod in unique(res$MODULE))
  {
    count = count +1
    
    res=fullRes[fullRes$MODULE==mod,c(1,2,3,4)]

    out.meta=data.frame(MODULE=unique(res$MODULE))
    index=1
    
    #To create a wide format MSEA result file
    for(study in GWAS_studies){ #this mapping is the step 1. 
      temp=res[res$MAP==study,]
      index=index+1
      out.meta=merge(out.meta,temp[,c(1,2,3)],by="MODULE",all.x=T)
      names(out.meta)[2*index-2]=paste("P", study, sep=".")
      names(out.meta)[2*index-1]=paste("FDR", study, sep=".")
    }
    
    #To calculate Meta-MSEA p value and FDR
    out.meta$METAP=as.numeric(unlist(apply(out.meta[2*1:(index-1)],1,function(x){
      if (sum(is.na(x))==0){ #Important, only modules with p-values in every mapping will be considered in MSEA
        x=x[!is.na(x)] #Remove empty p-values
        z=qnorm(x)
        z=sum(z,na.rm=T)/sqrt(length(z))
        return(pnorm(z))
      } else {
        print("addFlag false")
        return(NA)}
    })))
  
    if(addFlag){
      print(paste("Adding ", mod))
      fullData <- rbind(fullData, out.meta)
    }
    
    addFlag = TRUE
  }
  
  fullData$FDR=p.adjust(fullData$METAP,method="fdr")
  fullDataSorted = fullData[order(as.numeric(fullData$FDR), decreasing = FALSE),]
  write.table(fullDataSorted,paste("./meta_results/meta.Leprosy_", t, ".txt", sep = ""),row.names=F,quote=F,sep="\t") 
  #This file can assigned to the plan$label item as "meta.mega___" in the "merge" script
#}





