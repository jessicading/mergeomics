#source("/u/project/xyang123/xyang123-NOBACKUP/jading/Mergeomics.R")
#setwd("/u/flashscratch/j/jading/10X_Single_Cell/AD_HYP/")
#library(ggplot2) # used for GWAS_MSEA_Screen

# Marker Dependency Filtering
# 
# Filters markers based on dependency (such as linkage disequilibrium) file given. Prepares
# files for MSE (mapping and loci files)
#
# @param MARFILE: association file with 'LOCUS' and 'VALUE' headers. Value reflects strength of association.
#                 (for example, can use -log10 pvalues)
# @param GENFILE: file maping loci to genes - 'GENE' 'LOCUS'
# @param LNKFILE: marker dependency file 'LOCUSa' 'LOCUSb' 'WEIGHT'
# @param trait_name: optional, becomes part of the output directory name with mapping_name
# @param mapping_name: optional, becomes part of the output directory name with trait_name
# @param nmax: percent of associations to consider
# @param ld_threshold: corresponds to LD threshold used - must be changed if not 50%, 
#                     used to name output files
# @param ldprune: path to ldprune program
#
# @return creates directory with name that combines trait and mapping information and that 
#         contains the mapping (-.g.txt) and loci (-.l.txt) file with percent associations
#         and ld threshold information appended (ex. 50.50.g.txt)
#
# @examples
# runMDF(MARFILE = "./GWAS/Kunkle_AD.txt",
#        GENFILE = "./mapping/Brain_Hippocampus.eQTL.txt", 
#        LNKFILE = "./linkage/LD50.1000G.CEU.txt", 
#        output_dir = "./MSEA/Data/",
#        ldprune = "./MDPRUNE/ldprune")
#
runMDF <-function(MARFILE,
                  GENFILE,
                  LNKFILE,
                  trait_name=NULL, 
                  mapping_name=NULL, 
                  output_dir="./MSEA/Data/",
                  nmax=0.5,
                  ld_threshold=50,
                  editFiles=FALSE,
                  mdprune="mdprune"){
  
  # change to appropriate headers
  if(editFiles){
    system(paste0("sed -i \"1s/.*/GENE\\tMARKER/\" ",GENFILE))
    system(paste0("sed -i \"1s/.*/MARKERa\\tMARKERb\\tWEIGHT/\" ",LNKFILE))
  }
  
  if(is.null(trait_name) & is.null(mapping_name)){
    trait_name = unlist(strsplit(MARFILE,"/"))[length(unlist(strsplit(MARFILE,"/")))]
    trait_name = gsub(".txt","",trait_name)
    
    mapping_name = unlist(strsplit(GENFILE,"/"))[length(unlist(strsplit(GENFILE,"/")))]
    mapping_name = gsub(".txt","",mapping_name)
  }
  
  label=paste(output_dir,trait_name,'.',mapping_name,sep="")
  ifelse(!dir.exists(label), dir.create(label, recursive = TRUE),FALSE)
  
  bash_file <- file(paste0(label,".bash"))
  writeLines(c(paste('MARFILE="',MARFILE,'"',sep=''),
               paste('GENFILE="',GENFILE,'"',sep=''),
               paste('LNKFILE="', LNKFILE,'"',sep=""),
               paste('OUTPATH="',output_dir,trait_name,'.',mapping_name,'/"',sep=""),
               paste('NTOP=',nmax,sep=""),
               paste("echo -e \"MARKER\\tVALUE\" > /tmp/header.txt"),
               paste("nice sort -r -g -k 2 $MARFILE > /tmp/sorted.txt"),
               paste("NMARKER=$(wc -l < /tmp/sorted.txt)"),
               paste("NMAX=$(echo \"($NTOP*$NMARKER)/1\" | bc)"),
               paste("nice head -n $NMAX /tmp/sorted.txt > /tmp/top.txt"),
               paste("cat /tmp/header.txt /tmp/top.txt > /tmp/subset.txt"),
               paste('nice ',mdprune," /tmp/subset.txt $GENFILE $LNKFILE $OUTPATH",sep="")),
             bash_file)
  close(bash_file)
  
  ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive = TRUE),FALSE)
  
  name=paste(as.character(nmax*100), ld_threshold, sep=".")
  system(paste0("bash ", paste0(label,".bash")))
  system(paste("mv ",label,"/genes.txt ",label,"/",name,".g.txt",sep=""))
  system(paste("mv ",label,"/marker.txt ",label,"/",name,".m.txt",sep=""))
}

# Marker Set Enrichment Analysis
# 
# Finds significantly enriched marker sets in marker associations
#
# @param MDF_output_dir: output directory created from MDF that contains the corrected
#                        mapping and association files. association_file and mapping_file
#                        do not need to be specified if this is set.
# @param association_file: marker association file (gwas, transcriptome, proteome, 
#                          epigenome,metabolome) 'LOCUS' 'VALUE'
# @param mapping_file: needed if GWAS file used in association_file. not needed otherwise
# @param marker_set: marker sets to test for enrichment 'MODULE' 'GENE'
# @param info: information for marker_set file 'MODULE' 'SOURCE' 'DESCR'
# @param permtype: permutation type. automatically set as "locus" when doing gene level 
#                  enrichment analysis (not from GWAS)
# @param maxoverlap: corresponds to LD threshold used - must be changed if not 50%, 
#                     used to name output files
# @param ldprune: path to ldprune program
#
# @return creates directory with name that combines trait and mapping information and that 
#         contains the mapping (-.g.txt) and loci (-.l.txt) file with percent associations
#         and ld threshold information appended (ex. 50.50.g.txt)
#
# @examples
# runMDF(MARFILE = "./GWAS/Kunkle_AD.txt",
#        GENFILE = "./mapping/Brain_Hippocampus.eQTL.txt", 
#        LNKFILE = "./linkage/LD50.1000G.CEU.txt", 
#        output_dir = "./MSEA/Data/",
#        ldprune = "./MDPRUNE/ldprune")
#
runMSEA <- function(MDF_output_dir = NULL,
                    association_file,
                    mapping_file=NULL,
                    marker_set, 
                    output_dir="./MSEA/Results", 
                    label=NULL, 
                    info=NULL,
                    permtype = "gene",
                    maxoverlap = .33,
                    nperm=10000,
                    max_module_genes=5000,
                    trim = 0.002,
                    return_job = FALSE,
                    perc_top_associations=NULL){
  
  # either gene level enrichment or skipped MDF
  
  
  if(is.null(MDF_output_dir)){
    association_orig = association_file
    association <- read.delim(association_file, stringsAsFactors = FALSE)
    association = association[order(association$VALUE, decreasing = TRUE),]
    write.table(association,file = association_file, 
                row.names = FALSE, quote = FALSE, sep = "\t")
    if(!is.null(perc_top_associations)){ # if skipped MDF and doing gwas enrichment
      num = round(nrow(association)*perc_top_associations)
      association = association[1:num,]
      mapping_orig = mapping_file
      mapping <- read.delim(mapping_file, stringsAsFactors = FALSE)
      association = association[association$LOCUS %in% mapping$LOCUS,]
      
      association_file = gsub(pattern="//",replacement = "/", association_file)
      association_lab <- unlist(strsplit(association_file,"/"))[length(unlist(strsplit(association_file,"/")))]
      association_lab <- gsub(".txt", "", association_lab)
      write.table(association,
                  file = paste0(association_lab,"_for_MSEA_top",perc_top_associations*100,".txt"), 
                  row.names = FALSE, quote = FALSE, sep = "\t")
      association_file <- paste0(association_lab,"_for_MSEA_top",perc_top_associations*100,".txt")
      
      # make new mapping file
      mapping = mapping[mapping$LOCUS %in% association$LOCUS,]
      mapping_file = gsub(pattern="//",replacement = "/", mapping_file)
      mapping_lab <- unlist(strsplit(mapping_file,"/"))[length(unlist(strsplit(mapping_file,"/")))]
      mapping_lab <- gsub(".txt", "", mapping_lab)
      write.table(mapping,
                  file = paste0(mapping_lab,"_for_MSEA.txt"), 
                  row.names = FALSE, quote = FALSE, sep = "\t")
      mapping_file <- paste0(mapping_lab,"_for_MSEA.txt")
    }
    if(is.null(mapping_file)){ # if doing gene level enrichment analysis
      marker_associations <- read.delim(association_file, stringsAsFactors = FALSE)
      mapping_file = data.frame("GENE"=marker_associations$LOCUS, 
                                "LOCUS" = marker_associations$LOCUS, 
                                stringsAsFactors = FALSE)
      write.table(mapping_file, "./genfile_for_geneEnrichment.txt", 
                  row.names = FALSE, quote = FALSE, sep = "\t")
      mapping_file <- "./genfile_for_geneEnrichment.txt"
      permtype = "locus"
      maxoverlap = 1
      association_file = gsub(pattern="//",replacement = "/", association_file)
      association_lab <- unlist(strsplit(association_file,"/"))[length(unlist(strsplit(association_file,"/")))]
      association_lab <- gsub(".txt", "", association_lab)
      label = association_lab
    }else{
      if(is.null(label)){
        association_orig = gsub(pattern="//",replacement = "/", association_orig)
        association_lab <- unlist(strsplit(association_orig,"/"))[length(unlist(strsplit(association_orig,"/")))]
        association_lab <- gsub(".txt", "", association_lab)
        mapping_orig = gsub(pattern="//",replacement = "/", mapping_orig)
        mapping_lab <- unlist(strsplit(mapping_orig,"/"))[length(unlist(strsplit(mapping_orig,"/")))]
        mapping_lab <- gsub(".txt", "", mapping_lab)
        label = paste(association_lab, mapping_lab, sep = ".")
      }
    }
  }else{
    association_file = list.files(MDF_output_dir, 
                                  full.names = TRUE)[grep(".m.txt",
                                                          list.files(MDF_output_dir, 
                                                                     full.names = TRUE))]
    mapping_file = list.files(MDF_output_dir,
                              full.names = TRUE)[grep(".g.txt",
                                                      list.files(MDF_output_dir,
                                                                 full.names = TRUE))]
    if(is.null(label)){
      association_file = gsub(pattern="//",replacement = "/", association_file)
      association_lab <- unlist(strsplit(association_file,"/"))[length(unlist(strsplit(association_file,"/")))-1]
      association_lab <- gsub(".txt", "", association_lab)
      label = association_lab
    }
  }
  
  ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive = TRUE),FALSE)
  
  ass <- read.delim(association_file, stringsAsFactors = FALSE)
  colnames(ass) <- c("MARKER","VALUE")
  write.table(ass, file = association_file, row.names = FALSE, quote = FALSE, sep = "\t")
  
  map <- read.delim(mapping_file, stringsAsFactors = FALSE)
  colnames(map) <- c("GENE","MARKER")
  write.table(map, file = mapping_file, row.names = FALSE, quote = FALSE, sep = "\t")
  
  job.ssea <- list()
  job.ssea$label <- label
  job.ssea$folder <- output_dir
  job.ssea$genfile <- mapping_file
  job.ssea$marfile <- association_file #marfile in package, locfile in Mergeomics.R functions
  job.ssea$modfile <- marker_set
  if(!is.null(info)){
    job.ssea$inffile <- info
  }
  job.ssea$permtype <- permtype	#optional
  job.ssea$maxgenes <- max_module_genes	#optional
  job.ssea$nperm <- nperm	#optional
  job.ssea$maxoverlap <- maxoverlap
  job.ssea <- ssea.start(job.ssea)
  job.ssea <- ssea.prepare(job.ssea)
  job.ssea <- ssea.control(job.ssea)
  job.ssea <- ssea.analyze(job.ssea,trim_start=trim,trim_end=(1-trim))
  job.ssea <- ssea.finish(job.ssea)
  
  if(return_job) return(job.ssea)
}

# nodes_file can either be .results.txt file from MSEA or a 'MODULE' 'NODE' file
runKDA <- function(MSEA_results_dir = NULL,
                   nodes_file, # can either be .results.txt file from MSEA or a 'MODULE' 'NODE' file
                   marker_set=NULL,
                   network, 
                   name = NULL,
                   trim = c(".results.txt",".50.20",".50.50",".txt",".results"), 
                   merge = FALSE, 
                   rcutoff=.33, 
                   info=NULL,
                   edgefactor = 0,
                   depth = 1,
                   direction = 0,
                   nperm = 10000,
                   nKDs_show = NULL,
                   return_job = FALSE
){
  
  job.kda <- list()
  job.kda$label<-"wKDA"
  
  if(!is.null(MSEA_results_dir)){
    results = list.files(MSEA_results_dir, 
                                  full.names = TRUE)[grep(".results.txt",
                                                          list.files(MSEA_results_dir, 
                                                                     full.names = TRUE))]
    if(length(results)>1) cat("Error: More than one MSEA results files detected.\n")
  } else{
    results = nodes_file
  }

  
  result <- read.delim(results, stringsAsFactors = FALSE)
  if(ncol(result)>4){ # result file from MSEA
    sig_mods = result$MODULE[result$FDR<0.05]
    sig_modfile = data.frame(stringsAsFactors = FALSE)
    
    modfile <- read.delim(marker_set, stringsAsFactors = FALSE)
    
    for(mod in sig_mods){
      sig_modfile = rbind(sig_modfile, modfile[modfile$MODULE==mod,])
    }
    colnames(sig_modfile) <- c("MODULE","NODE")
    write.table(sig_modfile, "Significant_module_nodes.txt", row.names = FALSE, quote = FALSE, sep = "\t")
    job.kda$modfile <- "./Significant_module_nodes.txt"
    forMerging = data.frame("MODULE"=unique(sig_modfile$MODULE), stringsAsFactors = FALSE)
  }
  else if(ncol(result)==4){ # output from merge modules
    result <- result[,c("MODULE","NODE")]
    write.table(result, "nodes_file_forKDA.txt", row.names = FALSE, quote = FALSE, sep = "\t")
    job.kda$modfile <- "nodes_file_forKDA.txt"
  }
  else{
    if(length(intersect(colnames(result),c("MODULE","NODE")))!=2){
      colnames(result) <- c("MODULE","NODE")
      write.table(result, "nodes_file_forKDA.txt", row.names = FALSE, quote = FALSE, sep = "\t")
      job.kda$modfile <- "nodes_file_forKDA.txt"
    }else{
      job.kda$modfile <- nodes_file
    }
    
    forMerging <- read.delim(nodes_file, stringsAsFactors = FALSE)
    forMerging = data.frame("MODULE"=unique(forMerging$MODULE), stringsAsFactors = FALSE)
  }
  
  if(merge){
    if(ncol(result)==2){
      modfile <- result
      colnames(modfile) <- c("MODULE","GENE")
      write.table(modfile,"modfile_forMerging.txt",
                  row.names = FALSE, quote = FALSE, sep = "\t")
      marker_set = "modfile_forMerging.txt"
    }else{
      modfile <- read.delim(marker_set, stringsAsFactors = FALSE)
    }
    
    if(is.null(info)){
      infofile = data.frame("MODULE"=unique(modfile$MODULE), "SOURCE"=unique(modfile$MODULE), "DESCR"=unique(modfile$MODULE))
      write.table(infofile,"infofile_forMerging.txt", row.names = FALSE, quote = FALSE, sep = "\t")
      merge_modules(name = "nodes", modules_df = forMerging, rcutoff = rcutoff, output_Dir = "./", 
                    modfile_path = marker_set, infofile_path = "./infofile_forMerging.txt")
    }
    else{
      merge_modules(name = "nodes", modules_df = forMerging, rcutoff = rcutoff, 
                    output_Dir = "./", modfile_path = marker_set, infofile_path = info)
    }
    job.kda$modfile <- "merged_nodes.mod.txt"
  }
  
  nodes_name = nodes_file
  nodes_name = unlist(strsplit(nodes_name, split = "/"))[length(unlist(strsplit(nodes_name, split = "/")))]
  for(i in trim){
    nodes_name = gsub(i,"",nodes_name)
  }
  
  network_name = gsub(".txt","",unlist(strsplit(network, split = "/"))[length(unlist(strsplit(network, split = "/")))])
  
  cat("\nNow analyzing:",nodes_name, "with", network_name, "\n")
  
  if(is.null(name)){
    name = paste(nodes_name, network_name,sep = "_")
  }
  
  job.kda$folder<- name ## parent folder for results
  
  ## Input a network
  ## columns: TAIL HEAD WEIGHT
  job.kda$netfile <- network
  
  
  ## "0" means we do not consider edge weights while 1 is opposite.
  job.kda$edgefactor <- edgefactor
  ## The searching depth for the KDA
  job.kda$depth <- depth
  ## 0 means we do not consider the directions of the regulatory interactions
  ## while 1 is opposite.
  job.kda$direction <- direction
  job.kda$nperm <- nperm
  
  moddata <- tool.read(job.kda$modfile)
  mod.names <- unique(moddata$MODULE)
  moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))),]
  ## save this to a temporary file and set its path as new job.kda$modfile:
  tool.save(moddata, "subsetof.supersets.txt")
  job.kda$modfile <- "subsetof.supersets.txt"
  
  ## Run KDA
  job.kda <- kda.configure(job.kda)
  job.kda <- kda.start(job.kda)
  job.kda <- kda.prepare(job.kda)
  job.kda <- kda.analyze(job.kda)
  job.kda <- kda.finish(job.kda)
  
  if(!is.null(nKDs_show)){
    job.kda <- kda2cytoscape(job.kda, ndrivers = nKDs_show)
  }
  else{
    job.kda <- kda2cytoscape(job.kda)
  }
  
  if(return_job){
    saveRDS(job.kda, file = paste0(name, "_kda_job.rds"))
  }
}

GWAS_MSEA_screen <- function(GWAS_dir="/u/project/xyang123/xyang123-NOBACKUP/icestrik/public_gwas/MSEA_Gtex_v7/",
                             SNP_map_method = "Adipose_Subcutaneous",
                             geneset="", 
                             info = NULL,
                             ntop = 50,
                             ld_threshold = 50,
                             output_dir = "./GWAS_Screen_Results/",
                             output_name = "Results",
                             max_module_genes = 5000,
                             trim = .005,
                             nperm = 10000,
                             permtype = "gene",
                             nGWASshow = 20){
  traits = list.files(GWAS_dir)[grep(SNP_map_method, list.files(GWAS_dir))]
  trim_traits = c()
  for(trait in traits){
    trim_traits = append(trim_traits, unlist(strsplit(trait, split = ".", fixed = TRUE))[1])
  }
  traits = trim_traits
  rm(trim_traits)
  para = paste0(ld_threshold,".",ntop)
  mapping = SNP_map_method
  
  ifelse(!dir.exists(output_dir), dir.create(output_dir), FALSE)
  setwd(output_dir)
  for(trait.item in traits){
    for (item.1 in mapping){
      for (item.2 in para){
        job.ssea <- list()
        job.ssea$label <- paste(trait.item,item.1,item.2,sep=".")
        job.ssea$folder <- output_dir
        job.ssea$genfile <- paste(GWAS_dir,trait.item,".",item.1,"/",item.2,".g.txt",sep="")
        job.ssea$locfile <- paste(GWAS_dir,trait.item,".",item.1,"/",item.2,".l.txt",sep="")		
        job.ssea$modfile <- geneset
        if(!is.null(info)){
          job.ssea$inffile <- info
        }
        job.ssea$permtype <- permtype	#optional
        job.ssea$maxgenes <- max_module_genes	#optional
        job.ssea$nperm <- nperm	#optional
        job.ssea <- ssea.start(job.ssea)
        job.ssea <- ssea.prepare(job.ssea)
        job.ssea <- ssea.control(job.ssea)
        job.ssea <- ssea.analyze(job.ssea,trim_start=trim,trim_end=(1-trim))
        job.ssea <- ssea.finish(job.ssea)
      }
    }
  }
  
  # make plot of results
  result_files <- list.files(paste0(output_dir, "msea/"))[grep("results", list.files(paste0(output_dir, "msea/")))]
  results_df = data.frame(stringsAsFactors = FALSE)
  for(result in result_files){
    msea <- read.delim(paste0(output_dir,"msea/", result), stringsAsFactors = FALSE)
    study = unlist(strsplit(result, split = ".", fixed = TRUE))[1]
    temp = data.frame("Study" = study, 
                      "Module" = msea$MODULE,
                      "FDR" = msea$FDR)
    results_df = rbind(results_df, temp)
  }
  
  results_df = results_df[results_df$Module!="_ctrlA",]
  results_df = results_df[results_df$Module!="_ctrlB",]
  
  results_df$log10FDR <- -log10(results_df$FDR)
  results_df <- results_df[order(results_df$log10FDR, decreasing = TRUE),]
  # get top 10 GWAS studies 
  show_GWAS = unique(results_df$Study)[1:nGWASshow]
  
  toPlot = data.frame(stringsAsFactors = FALSE)
  for(gwas in show_GWAS){
    toPlot = rbind(toPlot, results_df[results_df$Study==gwas,])
  }
  
  toPlot$FDR = toPlot$FDR/nGWASshow
  toPlot$log10FDR <- -log10(toPlot$FDR)
  
  heat <- ggplot(toPlot,aes(x=Module,y = Study, fill=log10FDR)) + geom_tile(colour="gray") + 
    scale_fill_gradientn(colors = c("white", "red")) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 15, vjust = 0, hjust = 0, 
                                     margin = margin(t=0,b=0,r=0,l=0)), 
          axis.text.y =  element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.border = element_blank(),
          axis.title = element_text(size=15, face = "bold")) + 
    labs(fill=expression(-log[10]~FDR)) +
    scale_x_discrete(position = "top") 
  
  pdf(paste0(output_name,"_heatmap.pdf"), width = 10, height = 7)
  print(heat)
  dev.off()
  
  write.table(results_df,paste0(output_name,"_data.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
}

merge_modules <- function(msea_res, 
                          rcutoff=0.33, 
                          fdr_cutoff=NULL,
                          output_Dir="Merged_modules/", 
                          modfile_path, 
                          infofile_path=NULL,
                          label=""){
  
  if(length(msea_res)==1){ # msea result file
    res <- read.delim(msea_res)
    if(!is.null(fdr_cutoff)){
      res <- res[res$FDR<fdr_cutoff,]
    }
    pool <- as.character(res[,"MODULE"])
  } else{ # vector of modules
    pool=msea_res
  }
  pool <- pool[!(pool=="_ctrlA" | pool=="_ctrlB")]
  
  #=====================================================
  #=== Merge the modules before 2nd SSEA
  if (length(pool)>0){
    meg.mods<- tool.read(modfile_path)
    merged.modules <- pool
    moddata <- meg.mods[which(!is.na(match(meg.mods[,1], merged.modules))),]
    all_mod <- moddata
    
    # Merge and trim overlapping modules.
    moddata$OVERLAP <- moddata$MODULE
    moddata <- tool.coalesce(items=moddata$GENE, groups=moddata$MODULE, rcutoff=rcutoff)
    moddata$MODULE <- moddata$CLUSTER
    moddata$GENE <- moddata$ITEM
    moddata$OVERLAP <- moddata$GROUPS
    moddata <- moddata[,c("MODULE", "GENE", "OVERLAP")]
    moddata <- unique(moddata)
    
    if(!is.null(infofile_path)){
      moddatainfo <- tool.read(infofile_path)
      moddatainfo <- moddatainfo[which(!is.na(match(moddatainfo[,1], moddata[,1]))), ]
    }
    # Mark modules with overlaps.
    for(j in which(moddata$MODULE != moddata$OVERLAP)){
      if(!is.null(infofile_path)){
        moddatainfo[which(moddatainfo[,"MODULE"] == moddata[j,"MODULE"]), "MODULE"] <- paste(moddata[j,"MODULE"], "..", sep=",")
      }
      moddata[j,"MODULE"] <- paste(moddata[j,"MODULE"], "..", sep=",")
    }
    # Save module info for 2nd SSEA and KDA.
    moddata <- unique(moddata)
    moddata[, 4] <- moddata[, 2];  names(moddata)[4] <- c("NODE")
    if(label!="") label <- paste0(label,".")
    mdfile=paste0(label, "mod.txt"); mifile=paste0(label, "info.txt")
    
    if(length(setdiff(all_mod$GENE, moddata$GENE))>0){
      addBack <- data.frame("GENE"=unique(setdiff(all_mod$GENE, moddata$GENE)))
      cat("There were ", nrow(addBack), " genes lost from merging. Adding back...\n")
      modules <- c()
      for(gene in addBack$GENE){
        modules <- c(modules, do.call("paste",c(all_mod$MODULE[all_mod$GENE==gene],
                                                list("sep"=", "))))
      }
      addBack$Module <- modules
      merged_modules <- unique(moddata$OVERLAP)
      names(merged_modules) <- unique(moddata$MODULE)
      for(iter in 1:nrow(addBack)){
        # get module name from merged
        mod <- names(merged_modules)[grepl(unlist(strsplit(addBack$Module[iter], split = ", "))[1], merged_modules)]
        temp <- data.frame("MODULE"=mod, 
                           "GENE"=addBack$GENE[iter], 
                           "OVERLAP"=merged_modules[names(merged_modules)==mod],
                           "NODE"=addBack$GENE[iter], stringsAsFactors = FALSE)
        moddata <- rbind(moddata, temp)
      }
    }
    
    if(!dir.exists(output_Dir)) dir.create(output_Dir)
    
    write.table(moddata, paste0(output_Dir,"/merged_", mdfile),  
                sep='\t', col.names=T, row.names=F, quote=F) 
    if(!is.null(infofile_path)){
      write.table(moddatainfo, paste0(output_Dir,"/merged_", mifile),  
                  sep='\t', col.names=T, row.names=F, quote=F) 
    }
  }
}

