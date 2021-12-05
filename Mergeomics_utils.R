# Wrapper functions for Mergeomics pipeline

# Marker Dependency Filtering
# 
# Filters markers based on dependency (such as linkage disequilibrium). 
# Prepares files for MSEA (marker association and mapping files)
#
# @param MARFILE: association file with 'MARKER' and 'VALUE' headers. Value 
#                 reflects strength of association. (e.g., -log10 pvalues)
# @param GENFILE: file maping loci to genes - 'GENE' 'MARKER'
# @param LNKFILE: marker dependency file 'MARKERa' 'MARKERb' 'WEIGHT'
# @param label: directory name containing result files
# @param NTOP: percent of associations to consider
# @param md_threshold: corresponds to dependency (correlation) threshold, used
#                      to name output files (optional)
# @param mdprune: path to mdprune program
#
# @return void, creates directory containing the marker (_.m.txt) and mapping 
#         (_.g.txt) file with percent associations and dependency threshold (if 
#         provided) appended (ex. top50.md50.g.txt)
#
# @examples
# runMDF(MARFILE = "./GWAS/Kunkle_AD.txt",
#        GENFILE = "./mapping/Brain_Hippocampus.eQTL.txt", 
#        LNKFILE = "./linkage/LD50.1000G.CEU.txt", 
#        output_dir = "./MSEA/Data/",
#        mdprune = "./MDPRUNE/ldprune")
#
runMDF <-function(marker_associations,
                  marker_mapping,
                  marker_dependency,
                  label=NULL, 
                  output_dir="Data/",
                  n_top=0.5,
                  md_threshold=NULL, # correlation cutoff for output file name
                  edit_files=FALSE,
                  mdprune="mdprune"){
  
  # change to appropriate headers
  if(edit_files){
    system(paste0("sed -i \"1s/.*/MARKER\\VALUE/\" ",marker_associations))
    system(paste0("sed -i \"1s/.*/GENE\\tMARKER/\" ",marker_mapping))
    system(paste0("sed -i \"1s/.*/MARKERa\\tMARKERb\\tWEIGHT/\" ",marker_dependency))
  }
  
  if(is.null(label)){
    trait_name = unlist(strsplit(marker_associations,"/"))[length(unlist(strsplit(marker_associations,"/")))]
    trait_name = gsub(".txt","",trait_name)
    
    mapping_name = unlist(strsplit(marker_mapping,"/"))[length(unlist(strsplit(marker_mapping,"/")))]
    mapping_name = gsub(".txt","",mapping_name)
    
    label=paste(output_dir, "/",trait_name, '.', mapping_name, sep="")
  } else {
    label=paste(output_dir, "/",label)
  }
  
  ifelse(!dir.exists(label), dir.create(label, recursive = TRUE),FALSE)
  
  bash_file <- file(paste0(label,".bash"))
  writeLines(c(paste('MARFILE="',marker_associations,'"',sep=''),
               paste('GENFILE="',marker_mapping,'"',sep=''),
               paste('LNKFILE="', marker_dependency,'"',sep=""),
               paste('OUTPATH="',output_dir,trait_name,'.',mapping_name,'/"',sep=""),
               paste('NTOP=',n_top,sep=""),
               paste("echo -e \"MARKER\\tVALUE\" > /tmp/header.txt"),
               paste("nice sort -r -g -k 2 $MARFILE > /tmp/sorted.txt"),
               paste("NMARKER=$(wc -l < /tmp/sorted.txt)"),
               paste("NMAX=$(echo \"($NTOP*$NMARKER)/1\" | bc)"),
               paste("nice head -n $NMAX /tmp/sorted.txt > /tmp/top.txt"),
               paste("cat /tmp/header.txt /tmp/top.txt > /tmp/subset.txt"),
               paste('nice ',mdprune," /tmp/subset.txt $GENFILE $LNKFILE $OUTPATH",sep="")),
             bash_file)
  close(bash_file)
  
  if(!is.null(md_threshold)){
    name=paste0("top", as.character(n_top*100), ".md",md_threshold)
  } else {
    name=paste0("top", as.character(n_top*100))
  }
  
  marker_associations <- paste0(label,"/",name,".m.txt")
  marker_mapping <- paste0(label,"/",name,".g.txt")
  
  system(paste0("bash ", paste0(label,".bash")))
  system(paste0("mv ",label,"/genes.txt ",marker_mapping))
  system(paste0("mv ",label,"/marker.txt ",marker_associations))
  
  job <- list()
  job$marker_associations <- marker_associations
  job$marker_mapping <- marker_mapping
  return(job)
}

# Marker Set Enrichment Analysis
# 
# Finds significantly enriched marker sets from marker associations 
# (summary statistics of genome-, epigenome-, transcriptome-, proteome-, 
#  metabolome-wide)
#
# @param marker_associations: marker associations file path 
#                             'MARKER' 'VALUE' headed tab delimited .txt file
# @param marker_mapping: marker to gene mapping file path. Required for GWAS
#                        and EWAS. Not needed for TWAS/PWAS
#                        'GENE' 'MARKER' headed tab delimited .txt file
# @param marker_set: marker set file path 
#                    'MODULE' 'GENE' headed tab delimited .txt file
# @param marker_set_info: marker set info file path (optional)
#                         'MODULE' 'DESCR' headed tab delimited .txt file
# @param output_dir: output directory name 
# @param label: prefix appended to result file names
# @param permtype: permutation type. automatically set as "marker" if mapping
#                  file not provide (gene level enrichment analysis, i.e. not GWAS)
# @param nperm:
# @param maxoverlap: overlap ratio threshold for merging genes with shared markers. 
# @param max_module_genes: maximum number of genes for module to be included
# @param min_module_genes: minimum number of genes for module to be included
# @param trim: percentile taken from the beginning and end for trimming away a 
#              defined proportion of genes with significant trait association to 
#              avoid signal inflation of null background in gene permutation
# @param seed: seed for random number generator
# @param return_job: whether to return job list at the end of analysis
#
# @return job list with inputs and outputs, if return_job = TRUE. produces result
#         files, most notably _.results.txt which contains the full results including
#         top genes, corresponding markers and their association values. full list
#         of gene contributing to the marker set enrichment is in _.details.txt.
#
# @examples
# runMSEA(marker_associations = "./GWAS/Kunkle_AD.txt",
#         marker_mapping = "./mapping/Brain_Hippocampus.eQTL.txt", 
#         marker_set = "./linkage/LD50.1000G.CEU.txt")
#
runMSEA <- function(job=NULL,
                    marker_associations,
                    marker_mapping=NULL,
                    marker_set, 
                    marker_set_info=NULL,
                    output_dir="Results", 
                    label=NULL, 
                    permtype="gene",
                    nperm=10000,
                    maxoverlap=.33,
                    max_module_genes=500,
                    min_module_genes=10,
                    trim=0.002,
                    seed=1,
                    return_job=TRUE){
  
  # either gene level enrichment or skipped MDF
  if(!is.null(job)){
    marker_associations <- job$marker_associations
    marker_mapping <- job$marker_mapping
  }
  
  if(is.null(label)){
    label = "msea"
  }
  
  ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive = TRUE),FALSE)
  
  job.ssea <- list()
  job.ssea$marfile <- marker_associations
  if(!is.null(marker_mapping)){
    job.ssea$genfile <- marker_mapping
  } else {
    permtype <- "marker"
    maxoverlap <- 1
  }
  job.ssea$label <- label
  job.ssea$folder <- output_dir
  job.ssea$modfile <- marker_set
  if(!is.null(marker_set_info)){
    job.ssea$inffile <- marker_set_info
  }
  job.ssea$permtype <- permtype
  job.ssea$nperm <- nperm
  job.ssea$maxoverlap <- maxoverlap 
  job.ssea$trim <- trim
  job.ssea$seed <- seed
  job.ssea <- ssea.start(job.ssea)
  job.ssea <- ssea.prepare(job.ssea)
  job.ssea <- ssea.control(job.ssea)
  job.ssea <- ssea.analyze(job.ssea)
  job.ssea <- ssea.finish(job.ssea)
  
  if(return_job) return(job.ssea)
}

# Meta-MSEA
# 
# Finds significantly enriched marker sets from marker associations 
# (summary statistics of genome-, epigenome-, transcriptome-, proteome-, 
#  metabolome-wide)
#
# @param marker_associations: marker associations file path 
#                             'MARKER' 'VALUE' headed tab delimited .txt file
# @param marker_mapping: marker to gene mapping file path. Required for GWAS
#                        and EWAS. Not needed for TWAS/PWAS
#                        'GENE' 'MARKER' headed tab delimited .txt file
# @param marker_set: marker set file path 
#                    'MODULE' 'GENE' headed tab delimited .txt file
# @param marker_set_info: marker set info file path (optional)
#                         'MODULE' 'DESCR' headed tab delimited .txt file
# @param output_dir: output directory name 
# @param label: prefix appended to result file names
# @param permtype: permutation type. automatically set as "marker" if mapping
#                  file not provide (gene level enrichment analysis, i.e. not GWAS)
# @param nperm:
# @param maxoverlap: overlap ratio threshold for merging genes with shared markers. 
# @param max_module_genes: maximum number of genes for module to be included
# @param min_module_genes: minimum number of genes for module to be included
# @param trim: percentile taken from the beginning and end for trimming away a 
#              defined proportion of genes with significant trait association to 
#              avoid signal inflation of null background in gene permutation
# @param seed: seed for random number generator
# @param return_job: whether to return job list at the end of analysis
#
# @return job list with inputs and outputs, if return_job = TRUE. produces result
#         files, most notably _.results.txt which contains the full results including
#         top genes, corresponding markers and their association values. full list
#         of gene contributing to the marker set enrichment is in _.details.txt.
#
# @examples
# runMSEA(marker_associations = "./GWAS/Kunkle_AD.txt",
#         marker_mapping = "./mapping/Brain_Hippocampus.eQTL.txt", 
#         marker_set = "./linkage/LD50.1000G.CEU.txt")
#
runMetaMSEA <- function(msea_input_list, 
                        marker_set, 
                        marker_set_info=NULL, 
                        output_dir="Results",
                        label="meta",
                        return_job=TRUE){
  
  default_params = list(permtype="gene", nperm=10000, maxoverlap=0.33,
                        max_module_genes=500, min_module_genes=10,
                        trim=0.002, seed=1)
  
  meta_job_list <- list()
  
  for(msea in names(msea_input_list)){
    set_params <- names(msea_input_list[[msea]])
    set_params <- setdiff(set_params, c("marker_associations","marker_mapping"))
    params <- default_params
    for(p in set_params){ # overwrite default params
      params[[p]] <- msea_input_list[[msea]][[p]]
    }
    
    job.ssea <- list()
    job.ssea$label <- msea
    job.ssea$folder <- "individual_msea_results_for_meta"
    job.ssea$marfile <- msea_input_list[[msea]][["marker_associations"]]
    
    for(p in names(params)){
      job.ssea[[p]] <- params[[p]]
    }

    if(!is.null(msea_input_list[[msea]][["marker_mapping"]])){
      job.ssea$genfile <- msea_input_list[[msea]][["marker_mapping"]]
    } else {
      job.ssea[["permtype"]] <- "marker"
      job.ssea[["maxoverlap"]] <- 1
    }
    
    job.ssea$modfile <- marker_set
    if(!is.null(msea_input_list[[msea]][["marker_set_info"]])){
      job.ssea$inffile <- marker_set_info
    }
    
    cat("\nRunning MSEA for", msea, "\n")
    job.ssea <- ssea.start(job.ssea)
    job.ssea <- ssea.prepare(job.ssea)
    job.ssea <- ssea.control(job.ssea)
    job.ssea <- ssea.analyze(job.ssea)
    job.ssea <- ssea.finish(job.ssea)
    
    meta_job_list[[msea]] <- job.ssea
  }
  
  meta_job_list <- ssea.meta(jobs = meta_job_list, 
                             label = label, 
                             folder = output_dir)

  if(return_job) return(meta_job_list)
}

# Key Driver Analysis
# 
# Identifies key drivers with network neighbors significantly enriched for
# input module genes
#
# @param job: completed msea job to retrieve significant modules for KDA
# @param MSEA_results: path to MSEA result .txt to extract significant modules for KDA
#                      ('MODULE' and 'FDR' columns required - _.results.txt file)
# @param MSEA_fdr_cutoff: FDR cutoff to include modules in the MSEA
# @param marker_sets: marker sets file path if setting MSEA_results parameter
# @param merge_modules: whether to merge redundant modules - TRUE or FALSE
# @param merge_rcutoff: minimum ratio of overlap to merge modules
# @param nodes: either a 'MODULE' 'NODE' file or vector of genes 
# @param marker_set_file: file to retrieve module genes if MSEA_results is set 
#                         'MODULE' 'GENE' headed tab delimited .txt file
# @param marker_set_info_file: marker sets info file (optional)
#                              'MODULE' 'DESCR' headed tab delimited .txt file
# @param network: path to network file
#                 'HEAD' 'TAIL' 'WEIGHT' headed tab delimited .txt file
# @param label: prefix appended to result file names
# @param output_dir: directory storing result files
# @param edgefactor: influence of node strengths, 0.0 no influence, 1.0 full influence
# @param maxoverlap: maximum allowed overlap between two key driver neighborhoods
# @param minsize: minimum module size
# @param mindegree: minimum node degree to qualify as a hub
# @param maxdegree: maximum node degree to include
# @param depth: search depth for subgraph search
# @param direction: use 0 for undirected, -1 for downstream and 1 for upstream
# @param nperm: number of permutations
# @param seed: seed for random number generator
# @param nKDs_subnetwork: maximum number of drivers per module to generate subnetwork
# @param return_job: whether to return job list at the end of analysis
# @param save_job: whether to save job list at the end of analysis
#
# @return job list with inputs and outputs, if return_job = TRUE. produces result
#         files, most notably _.results.txt which contains the full results of
#         they drivers. produces cytoscape-ready network files of top key drivers
#
# @examples
# kda_job <- runKDA(job=msea_job, 
#                   network="network/bayesian_network_brain.txt")  
#
# nodes_file can either be .results.txt file from MSEA or a 'MODULE' 'NODE' file
runKDA <- function(job=NULL,
                   MSEA_results=NULL,
                   MSEA_fdr_cutoff=0.05,
                   marker_sets=NULL,
                   merge_modules=FALSE, 
                   merge_rcutoff=.33,
                   nodes=NULL, # 'MODULE' 'NODE' file
                   marker_set_file=NULL,
                   marker_set_info_file=NULL,
                   network, 
                   label=NULL,
                   output_dir ="Results",
                   maxoverlap=0.33,
                   minsize=20,
                   mindegree="automatic",
                   maxdegree="automatic",
                   edgefactor=0,
                   depth=1,
                   direction=0,
                   nperm=10000,
                   seed=1,
                   nKDs_subnetwork=5,
                   return_job=TRUE,
                   save_job = TRUE
){
  
  job.kda <- list()
  job.kda$folder <- output_dir
  if(is.null(label)){
    if(!is.null(job)){
      if(label=="msea"){
        label <- "wKDA"
      } else {
        label <- job$label
      }
    } else {
      label <- "wKDA"
    }
  }
  
  if(!is.null(nodes)){
    if(length(nodes)>1){ # if vector of nodes/genes
      modfile <- data.frame("MODULE"="Input nodes",
                            "NODE"=nodes)
      write.table(modfile, "temp/nodes_file_forKDA.txt", row.names = FALSE, quote = FALSE, sep = "\t")
      job.kda$modfile <- "temp/nodes_file_forKDA.txt"
    } else { # if MODULE NODE file
      job.kda$modfile <- nodes
    }
  } else {
    if(!is.null(job)){
      marker_set_file <- job$modfile
      MSEA_results <- job$resultfile
      marker_set_info_file <- job$inffile
      job.kda$msea_results <- job$msea_results
    }
    
    if(is.null(marker_set_file)) stop("Marker set needed!")
    
    modfile <- read.delim(marker_set_file, stringsAsFactors = FALSE)
    
    if(!is.null(MSEA_results)){
      result <- read.delim(MSEA_results, stringsAsFactors = FALSE)
      modules = result$MODULE[result$FDR<MSEA_fdr_cutoff]
      modules = modules[modules!="_ctrlA"]
      modules = modules[modules!="_ctrlB"]
      if(length(modules)==0){
        if(MSEA_fdr_cutoff>=0.25){
          stop("No modules left after applying MSEA FDR cutoff.")
        } else {
          cat("\nNo modules left after applying MSEA FDR cutoff.\n")
          cat("Changing FDR cutoff to 0.25, please interpret results accordingly.\n")
          MSEA_fdr_cutoff = 0.25
          result <- read.delim(MSEA_results, stringsAsFactors = FALSE)
          modules = result$MODULE[result$FDR<MSEA_fdr_cutoff]
          modules = modules[modules!="_ctrlA"]
          modules = modules[modules!="_ctrlB"]
          if(length(modules)==0){
            stop("\nNo modules left after applying a 0.25 MSEA FDR cutoff.\n")
          }
        }
      }
      modfile <- modfile[modfile$MODULE %in% modules,]
      
    } else if(!is.null(marker_sets)) {
      modfile <- modfile[modfile$MODULE %in% marker_sets,]
    } else{
      stop("No input nodes provided!")
    }
    
    if(merge_modules){
      dir.create("merged_modules")
      merge_modules(msea_res = modfile$MODULE, 
                    rcutoff = merge_rcutoff, 
                    output_Dir = "merged_modules/", 
                    label = label,
                    modfile_path = marker_set_file, 
                    infofile_path = marker_set_info_file)
      job.kda$modfile <- paste0("merged_modules/",label,".merged_mod.txt")
      
    } else {
      dir.create("temp")
      colnames(modfile) <- c("MODULE","NODE")
      write.table(modfile, paste0("temp/",label,"_nodes_file_forKDA.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
      job.kda$modfile <- paste0("temp/",label,"_nodes_file_forKDA.txt")
    }
  }
  
  if(!is.null(marker_set_info_file)){
    job.kda$inffile <- marker_set_info_file
  }
  
  job.kda$label <- label
  
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
  
  # moddata <- tool.read(job.kda$modfile)
  # mod.names <- unique(moddata$MODULE)
  # moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))),]
  # ## save this to a temporary file and set its path as new job.kda$modfile:
  # tool.save(moddata, "subsetof.supersets.txt")
  # job.kda$modfile <- "subsetof.supersets.txt"
  
  ## Run KDA
  job.kda <- kda.configure(job.kda)
  job.kda <- kda.start(job.kda)
  job.kda <- kda.prepare(job.kda)
  job.kda <- kda.analyze(job.kda)
  job.kda <- kda.finish(job.kda)
  
  if(save_job){
    saveRDS(job.kda, file = paste0(output_dir,"/",label, ".kda.job.rds"))
  }
  
  job.kda <- kda2cytoscape(job.kda, ndrivers = nKDs_subnetwork)
  job.kda <- kda2cytoscape(job.kda)
  
  if(return_job) return(job.kda)
  
}

merge_modules <- function(msea_res, 
                          rcutoff=0.33, 
                          fdr_cutoff=NULL,
                          output_dir="Merged_modules/", 
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
    
    if(!dir.exists(output_dir)) dir.create(output_dir)
    
    write.table(moddata, paste0(output_dir,"/merged_", mdfile),  
                sep='\t', col.names=T, row.names=F, quote=F) 
    if(!is.null(infofile_path)){
      write.table(moddatainfo, paste0(output_dir,"/merged_", mifile),  
                  sep='\t', col.names=T, row.names=F, quote=F) 
    }
  }
}

