## ----------------- DOWNSTREAM ANALYSIS OF MSEA RESULTS FROM MERGEOMICS ----------------- ##
#                                                                                           #
# Collection of functions dealing with downstream analysis of MSEA results                  #
#                                                                                           #
# Functions:                                                                                #
#                                                                                           #
#     1) directOverlap - outputs exact module (gene set) matches between studies            #
#                        (i.e. "replication")                                               #
#     2) geneOverlap   - does "looser" overlap analysis by % gene overlap between           #
#                        modules                                                            #
#     3) SGMrepTable   - makes table that records what SGM methods are being represented    #
#                        for each module                                                    #
#     4) summaryTable  - makes table for each SGM method detailing supersets and            #
#                        comprising modules                                                 #
#                                                                                           #
# Written by Jessica Ding, 2019                                                             #
#                                                                                           #
# ------------------------------------------------------------------------------------------#

source("~/Desktop/Yang_Lab/T2D_AD/Data/T2D/meta_canonical/Mergeomics.R")

#  Notes:
#  - These sets of functions are highly dependent on file names (case sensitive) 
#  - If functions are done in the logical order that they were created, then this should not be 
#    a problem. Problems may arise if you have file names in the folder that have string matches
#    to the string I am using to query the files to be analyzed. For example, a query for ".mod"
#    are for files like AD_Adipose_Subcutaneous.mod.txt. If you have a file in the folder that 
#    is called "file.modified.txt" - it will capture this file and result in an error.


#-----------------------------------------------------------------------------------------------------#
# Abbreviation: SGM - SNP to gene mapping
# This function does a % gene overlap for two studies for each 
# SGM method (as opposed to "direct overlap" - exact name 
# matching of gene sets/modules)
#
# Inputs: 1) results_Dir         folder with MSEA results for each study with the SGM 
#                                as the common factor(s) between studies. file name format 
#                                should be "Study, SGM, .results.txt" (.result.txt is the 
#                                output extention of MSEA). 
#                                ex. "UKB_T2D.Adipose_Subcutaneous.results.txt"
#         2) output_Dir          folder to store results
#         3) perc_gene_overlap   minimum gene overlap to merge modules
#         4) FDRcutoff           usually 0.05
#         5) studies             vector with study names that are reflected in the msea results
#                                ex. c("DIAGRAMstage1_T2D", "UKB_T2D")
#         6) modfile_path        path to mod file, columns: 'MODULE' 'GENE' 
#         7) infofile_path       path to info file, columns: 'MODULE' 'GENE' 'DESCR'
#
# Outputs: 1) "combined_..."               combined results file for each study 
#          2) "merged_..."                 mod and info file from module merging 
#          3) "analyzed_..."               shows all supersets (does not include unmerged 
#                                          modules) and study of origin. "OVERLAP" indicates 
#                                          module came from all studies
#          4) "analyzed_overlap_summary"   a summary file of sets of modules that were merged with 
#                                          different study representations (indicating new overlaps 
#                                          based on % gene overlap)
#          5) "shared_..."                 'MODULE' 'METHOD' 'DESCR' file that can be used to merge
#                                          modules a second time given new overlaps
#
# Notes: 1) % gene overlap is set with "rmax" - ex. .3 is 30% gene overlap
#        2) This function is written to consider multiple studies but all msea
#           results files must be in the same folder
#        3) Directory path string parameters must end in "/"
#        4) If DESCR unknown, put MODULE name
#
# Version April 2019
#---------------------------------------------------------------------------------------------------#

geneOverlap <- function(results_Dir, output_Dir, perc_gene_overlap, fdr_cutoff, cohorts, modfile_path, infofile_path, trim, study){
  
  # MAKE COMBINED FILES FOR ALL STUDIES 
  files = list.files(results_Dir)
  
  ##get SGM methods. This assumes your results file is formatted - "study.SGM.results.txt" 
  study_subset = files[grep(cohorts[1], files)]
  SGMs = study_subset
  for(a in trim){
    SGMs = gsub(a, "", SGMs)
  }
  cat("The SNP to gene mappings used is/are:\n")
  print(SGMs) #this is to check that the SGMs are being extracted in the desired way
  
  ##combine studies for each SGM
  cat("Now combining studies together in one file for each SGM.\n")
  for(i in SGMs){
    results = files[grep(i,files)]
    all = data.frame()
    for(j in results){
      temp = read.delim(paste(results_Dir, j ,sep = ""), header = TRUE, stringsAsFactors = FALSE)
      temp = temp[temp$FDR<fdr_cutoff,]
      temp=temp[!grepl("_ctrlA",temp$MODULE)&!grepl("_ctrlB",temp$MODULE),]
      if(nrow(temp)>0){
        for(r in cohorts){
          if(grepl(r, j)){
            temp$STUDY = r
          }
        }
      }
      all = rbind(all, temp)
    }
    write.table(all, paste0(output_Dir,"combined_", i, ".txt"), row.names = FALSE, quote = FALSE, sep="\t")
  }
  
  combined_files = list.files(output_Dir)
  combined_files = combined_files[grep("combined_", combined_files)]
  
  cat("Now merging modules for % gene overlap...\n")
  # MERGE MODULES - this is modified from a module merge script 
  for(l in combined_files){
    cat("Now merging:\n")
    cat(l,"\n")
    plan = c()
    plan$folder = output_Dir #this includes ending "/"
    plan$label = l
    plan$modfile = modfile_path
    plan$inffile = infofile_path
    #=====================================================
    ## Get modules
    pool=c()
    file.name=paste0(plan$folder, plan$label)
    if(file.exists(file.name)){
      aa=(read.table(file.name, header=T, sep='\t', check.names=F, quote=NULL)) 
      if(length(aa$MODULE)==0 | length(aa$MODULE)==1) next
      if (nrow(aa) > 0) pool=unique(c(pool, as.character(aa[,"MODULE"])))
    }
    
    #=====================================================
    #=== Merge the modules
    if (length(pool)>0){
      meg.mods<- tool.read(plan$modfile)
      merged.modules <- pool
      moddata <- meg.mods[which(!is.na(match(meg.mods[,1], merged.modules))),]
      # Merge and trim overlapping modules.
      rmax <- perc_gene_overlap #set desired gene overlap
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
      for(q in which(moddata$MODULE != moddata$OVERLAP)){
        moddatainfo[which(moddatainfo[,"MODULE"] == moddata[q,"MODULE"]), "MODULE"] <- paste(moddata[q,"MODULE"], "..", sep=",")
        moddata[q,"MODULE"] <- paste(moddata[q,"MODULE"], "..", sep=",")
      }
    #save results
    moddata <- unique(moddata)
    moddata[, 4] <- moddata[, 2];  names(moddata)[4] <- c("NODE")
    mdfile=paste(gsub(".txt","", gsub("combined_", "",l)), ".mod.txt", sep = "")
    mifile=paste(gsub(".txt","",gsub("combined_", "",l)), ".info.txt", sep = "")
    
    write.table(moddata, paste0(plan$folder,"merged_", mdfile),  
                sep='\t', col.names=T, row.names=F, quote=F) 
    
    write.table(moddatainfo, paste0(plan$folder,"merged_", mifile),  
                sep='\t', col.names=T, row.names=F, quote=F) 
    }
  }
  
  cat("Now analyzing for newly overlapped modules...\n")
  # FIND IF THERE WERE MERGED MODULES WITH DIFFERENT STUDY REPRESENTATION
  files1 = list.files(output_Dir)
  merged_mod_files = files1[grep("merged_",files1)] #works if no SGMs have "mod" in name
  merged_mod_files = merged_mod_files[grep(".mod", merged_mod_files)]

  uncommon_df = data.frame(stringsAsFactors = FALSE)
  for(m in merged_mod_files){
    
    file = read.delim(paste0(output_Dir,m), header = TRUE, stringsAsFactors = FALSE)
    temp = unique(file$OVERLAP[grep(",", file$OVERLAP)])

    SGM = gsub("merged_","",gsub(".mod.txt","", m))
    combined_file_name = combined_files[grep(SGM, combined_files)]
    combined_df = read.delim(paste0(output_Dir,combined_file_name), header = TRUE, stringsAsFactors = FALSE)
    study_reps = c()
    for(n in temp){
      modules = unlist(strsplit(n,","))
      study_rep = ""
      for(o in modules){
        if(length(combined_df$STUDY[combined_df$MODULE==o])>1){
          study_rep = paste(study_rep, "OVERLAP",sep = ", ")
        }
        else{
          study_rep = paste(study_rep, combined_df$STUDY[combined_df$MODULE==o], sep = ", ")
        }
      }
      study_rep = substring(study_rep, 3)
      study_reps = append(study_reps, study_rep)
    }
    analyzed = data.frame("Merged Modules"= temp,
                          "Study Representation"=study_reps)
    print(SGM)
    write.table(analyzed, paste0(output_Dir,"analyzed_", SGM, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
    bool = c()
    for(p in study_reps){
      p = gsub(" ", "", p)
      splitted = unlist(strsplit(p, ","))
      if(sum(!splitted[1]==splitted)>0){
        bool = append(bool, TRUE)
      }
      else{
        bool = append(bool, FALSE)
      }
    }
    if(sum(bool)>0){
      uncommon = analyzed[bool,]
      uncommon$SNPtoGeneMapping = SGM
      uncommon_df = rbind(uncommon_df, uncommon)
    }
  }
  
  #output summary of all different study representation "overlaps"
  write.table(uncommon_df, paste0(output_Dir,"analyzed_overlap_summary.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  
  #make overlap module files including all direct and quasi overlap modules
  for(s in SGMs){
    cat("Now analyzing\n")
    cat(s, "\n")
    #get direct overlaps
    comb_file = read.delim(paste0(output_Dir,combined_files[grep(s, combined_files)]), header = TRUE, stringsAsFactors = FALSE)
    studies1 = unique(comb_file$STUDY)
    num = 1
    shared_mods = c()
    while(num<length(studies1)){
      if(num==1){
          modules = comb_file$MODULE[comb_file$STUDY==studies1[num]]
      }
      else {
        modules = shared_mods
      }
      temp_modules = comb_file$MODULE[comb_file$STUDY==studies1[num+1]]
      shared_mods = intersect(modules,temp_modules)
      num = num + 1 
    }
    if(length(shared_mods>0)){ #if there are any direct overlaps
    shared = data.frame("MODULE" = shared_mods,
                        "METHOD" = "direct")
    }
    
    uncommon_df = read.delim(paste0(output_Dir,"analyzed_overlap_summary.txt"), header = TRUE, stringsAsFactors = FALSE) #to convert back into strings
    if(sum(grepl(s, uncommon_df$SNPtoGeneMapping))>=1){ #if there are any quasi overlaps
      shared_quasi_mods = c()
      SGM_mod_rep = uncommon_df$Merged.Modules[grep(s, uncommon_df$SNPtoGeneMapping)]
      print(SGM_mod_rep)
      SGM_study_rep = uncommon_df$Study.Representation[grep(s, uncommon_df$SNPtoGeneMapping)]
      for(u in 1:length(SGM_mod_rep)){
        SGM_merged_mod_rep = unlist(strsplit(SGM_mod_rep[u], split = ","))
        SGM_merged_study_rep = unlist(strsplit(SGM_study_rep[u], split = ", "))
        not_overlap = !(SGM_merged_study_rep=="OVERLAP")
        shared_quasi_mods = append(shared_quasi_mods, SGM_merged_mod_rep[not_overlap])
      }
      
      shared_quasi = data.frame("MODULE"= shared_quasi_mods,
                                "METHOD"= "quasi")
      if(length(shared_mods)>0){ #if both quasi and direct
        shared = rbind(shared, shared_quasi)
      }
      else{ # if only quasi
        shared = shared_quasi
      }
    }

    if(length(shared$MODULE)>0){
    #annotate
      info = read.delim(infofile_path, header = TRUE, stringsAsFactors = FALSE)
      annot = c()
      for(v in shared$MODULE){
        annot = append(annot, info$DESCR[info$MODULE==v])
      }
      shared$DESCR = annot
      write.table(shared, paste0(output_Dir, "shared_", study, ".", s, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
    }
    else{
      next
    }
    shared = data.frame()
  }
  cat("Finished. shared_ files made.")
}

source("~/Desktop/Yang_Lab/Mergeomics/geneOntology.R")
source("~/Desktop/Yang_Lab/resources/genesets/R-tomfunctions.R")

#results_Dir with your merged mod files
#ontologyfname = "~/Desktop/Yang_Lab/resources/genesets/MSigDB_Canonical_Pathways.txt"
#--------------------------------------------------------------------------#
#
# 
# Notes:If you get the error "Error in cbind(rep(modulesize, dim(fm)[1]), fm) : object 'fm' not found", the 
# offending module is printed before this error and this module must be taken out from the .mod.txt file 
#
annotate_supersets <- function(results_Dir, ontologyfname, study, trim){ #make trim argument that has a vector of things to trim to get the "SGM" and also so that you don't wind up with .mod.txt in subseqeunt result file names
  files = list.files(results_Dir)
  mod_files = files[grep(".mod.txt", files)]

  for(i in mod_files){
    mod_file = read.delim(paste0(results_Dir,i), header = TRUE, stringsAsFactors = FALSE)
    if(sum(grepl(",", mod_file$MODULE))>0){
      mod_file = mod_file[grep(",", mod_file$OVERLAP),]
      if(length(mod_files)>1){
        SGM = i
        for(k in trim){
          SGM = gsub(k, "", SGM)
        }
        SGM = gsub(study, "", SGM)
        name = paste0("supersets.", study, ".", SGM)
      }
    else{
      mod_file = mod_file[grep(",", mod_file$OVERLAP),]
      name = paste0("supersets.", study)
    }
    file_name = paste0(results_Dir,name, ".txt")
    mod_file = mod_file[,c(1,2)]
    write.table(mod_file, file_name, row.names = FALSE, quote = FALSE, sep = "\t")
    geneOntology(work_Dir = results_Dir, output_Dir = results_Dir, ontologyfname = ontologyfname, file = file_name, name = name)
    }
    
    else{
      next
    }
  }
  
  #input for geneOntology is a .txt file

}


summaryTable <- function(results_Dir, output_Dir, study, infofile_path){
  #get SGMs
  files = list.files(results_Dir)
  mod_files = files[grep("mod", files)]
  SGMs = gsub("merged_", "", gsub(".mod.txt","",mod_files)) 
  all_table = data.frame()
  for(i in SGMs){
    table = data.frame(stringsAsFactors = FALSE)
    individual_table = data.frame()
    mod_file = read.delim(paste0(results_Dir,mod_files[grep(i, mod_files)]), header = TRUE, stringsAsFactors = FALSE)
    mods = unique(mod_file$OVERLAP)
    supersets = mods[grep(",", mods)] 
    mod_names = unique(mod_file$MODULE)
    superset_names = mod_names[grep(",", mod_names)]
    superset_annotation = c()
    if(length(supersets)>0){
      SGM_superset_table = data.frame()
      for(j in supersets){
        superset_table = data.frame()
        sub_modules_DESCR = c()
        name = unique(mod_file$MODULE[mod_file$OVERLAP==j])
        sub_modules = unlist(strsplit(j, split = ","))
        for(l in sub_modules){
          info = read.delim(infofile_path, header = TRUE, stringsAsFactors = FALSE)
          sub_modules_DESCR = append(sub_modules_DESCR, info$DESCR[info$MODULE==l])
        }

        annot = read.delim(paste0(results_Dir,"supersets.",study,".", i,"/","a.supersets.", study, ".",i, 
                                       "_Ontology_", name,".xls"), header = TRUE, stringsAsFactors = FALSE) #probably will cause problems
        size = c()
        annotation = c()
        if(is.null(annot$ModuleSize[1])){
          annotation = append(annotation,"Unknown")
          size = append(size, "Error in annotation gave no module size")
        }
        else if(is.na(annot$ModuleSize[1])){
            annotation = append(annotation,"Unknown")
            size = append(size, annot$ModuleSize[1])
         }
        else if(sum(annot$pvalue_corrected<0.05)==0){
            annotation = append(annotation,"Unknown")
            size = append(size, annot$ModuleSize[1])
        }
        else{
            sig_mods = annot[annot$pvalue_corrected<0.05,]
            max_mod_overlap_annot = annot$Gene.Category[annot$ModuleOverlap==max(annot$ModuleOverlap)]
            if(length(max_mod_overlap_annot)>1){
              names = c()
              for(k in max_mod_overlap_annot){
                names = paste(names, k, sep = ", ") #?
              }
              names = substring(names, 3)
              annotation = append(annotation, names)
            }
            else{
              annotation = append(annotation,max_mod_overlap_annot)
            }
            size = append(size, annot$ModuleSize[1])
          
        
        superset_table = data.frame("SGM" = i,
                                    "MODULE"= sub_modules,
                                    "DESCR" = sub_modules_DESCR,
                                    "SUPERSET"= annotation,
                                    "MODULE_SIZE" = size, stringsAsFactors = FALSE)
      }
      SGM_superset_table = rbind(SGM_superset_table, superset_table)
      }
      table = SGM_superset_table
    }
    
    individual_modules = mods[!grepl(",", mods)]
    if(length(individual_modules)>0){
      size2 = c()
      individ_modules_DESCR = c()
      for(m in individual_modules){
        individ_modules_DESCR = append(individ_modules_DESCR, info$DESCR[info$MODULE==m])
        size2 = append(size2, sum(mod_file$MODULE==m))
      }
      
      individual_table = data.frame("SGM" = i,
                                    "MODULE"= individual_modules,
                                    "DESCR" = individ_modules_DESCR,
                                    "SUPERSET" = "none",
                                    "MODULE_SIZE"=size2, stringsAsFactors = FALSE)
      
      table = rbind(table, individual_table)
    }
   if(grepl("quasi", study) & !(grepl("AD_T2D", study))){ #not applicable 
     shared_files = files[grep("shared_",files)]
     shared_file = shared_files[grep(i, shared_files)]
     shared = read.delim(paste0(results_Dir,shared_file), header=TRUE, stringsAsFactors = FALSE)
     methods = c()
     for(n in table$MODULE){
       methods = append(methods, shared$METHOD[shared$MODULE==n])
     }

     table$OVERLAP_METHOD = methods
   }
   write.table(table,paste0(output_Dir, "table_", study, ".", i, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
   all_table = rbind(all_table, table)
   all_table= rbind(all_table,"")
  }
  
  #include SGMs in which there was only one module overlap (was included above^)
  shared_files = files[grep("shared_", files)] 
  #summary = read.delim(paste0(results_Dir,"/results/table_summary_",study,".txt"), stringsAsFactors = FALSE, header = TRUE)
  for(i in shared_files){
    df = data.frame(stringsAsFactors = FALSE)
    file = read.delim(paste0(results_Dir, i), header = TRUE, stringsAsFactors = FALSE)
    if(length(file$MODULE)==1){
      SGM = gsub(paste0("shared_", study, "."),"",i)
      SGM = gsub(".txt","", SGM)
      if(grepl("quasi",study) & !grepl("AD_T2D",study)){
        df = data.frame("SGM"=SGM,
                        "MODULE"=file$MODULE,
                        "DESCR"=file$DESCR,
                        "SUPERSET" = "none",
                        "MODULE_SIZE"="value",
                        "OVERLAP_METHOD"=file$METHOD, stringsAsFactors = FALSE)
      }
      else{
        info =read.delim("/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.info.txt", header = TRUE, stringsAsFactors = FALSE)
        df = data.frame("SGM"=SGM,
                        "MODULE"=file$MODULE,
                        "DESCR"=info$DESCR[info$MODULE==file$MODULE],
                        "SUPERSET" = "none",
                        "MODULE_SIZE"="value", stringsAsFactors = FALSE
        )
      }
      all_table = rbind(all_table, df)
      all_table = rbind(all_table,"")
    }
    else{
      next
    }
  }
  
  #annotate coexpression modules
  modules = all_table$MODULE
  for(j in modules){
    if(grepl("GTEXv7",j)){
      descrip = c()
      if(j==""){
        descrip = append(descrip, "")
      }
      else{
        annot = read.delim(paste("/Users/jessicading/Desktop/Yang_Lab/resources/coexpression_modules/annotations/a.coexpr.mod_Ontology_", j, ".xls", sep = ""), header = TRUE, stringsAsFactors = FALSE)
        if(is.na(annot$ModuleSize[1])){
          descrip = append(descrip, "Unknown")
        }
        else if(sum(annot$pvalue_corrected<0.05)==0){
          descrip = append(descrip, "Unknown")
        }
        else{
          sig_mods = annot$Gene.Category[annot$pvalue_corrected<0.05]
          if(length(sig_mods)>1){
            names = c()
            for(k in sig_mods){
              names = paste(names, k, sep = ", ")
            }
            names = substring(names, 3)
          }
          else{
            names = sig_mods
          }
          descrip = append(descrip, names)
        }
      }
      all_table$DESCR[all_table$MODULE==j] = descrip
      
    }
  }
  all_table$MODULE_SIZE[is.na(all_table$MODULE_SIZE)]<-""
  
  
  #-------------------------------------------------------#
  write.table(all_table, paste0(output_Dir, "table_summary_",study, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
}



#!!!need to include: 
#For T2D quasi: adipose visceral omentum, adrenal gland, brain frontal cortex, brain hippocampus, esophagus muscularis, spleen
#For AD quasi: pancreas, 
#for T2D direct: adipose visceral omentum, adrenal gland, brain cerebellar hemisphere, brain frontal cortex, brain hippocampus, esophagus muscularis, spleen
#For AD direct: hippocampus, pancreas, 
#For AD_T2D quasi: all_brain_13, hippocampus, !! add nerve tibial at the end
#For AD_T2D direct: all_brain_13, cerebellum, HP, 


#set output_Dir to results folder 
mergedSummaryTable <- function(summary_table, study, output_Dir, modfile_path, infofile_path, ontologyfname){
  summary_tbl = read.delim(paste0(output_Dir,summary_table), header=TRUE, stringsAsFactors = FALSE)
  modules = summary_tbl$MODULE[!(summary_tbl$MODULE=="")]
  modules = unique(modules)
  modules_df = data.frame("MODULE"=modules)
  cat("The study being analyzed is ", study, "\n")
  merge_modules(name = study, modules_df = modules_df, output_Dir = output_Dir, modfile_path = modfile_path, infofile_path = infofile_path)
  #annotate supersets
  annotate_supersets(results_Dir = output_Dir, ontologyfname = ontologyfname, study = study)
  merged = read.delim(paste0(output_Dir, "merged_",study, ".mod.txt"), header = TRUE, stringsAsFactors = FALSE)
  supersets = unique(merged$OVERLAP[grep(",", merged$OVERLAP)])
  all_melted_df = data.frame(stringsAsFactors = FALSE)
  number = 1
  for(i in supersets){ #maybe make function for this? break down supersets into individual modules
    name = unique(merged$MODULE[merged$OVERLAP==i])
    sub_modules = unlist(strsplit(i, split = ","))
    
    melted_df = data.frame("ID" = paste0("S",number),
                           "MODULE"=sub_modules,
                           "SUPERSET"=i,
                           "DESCR"=extract_annotation(study = study, name = name, results_Dir = output_Dir)[,1], stringsAsFactors = FALSE)
    all_melted_df = rbind(all_melted_df, melted_df)
    number = number + 1
  }
  standalone_modules = unique(merged$MODULE[!grepl(",", merged$MODULE)])
  descriptions = c()
  for(j in standalone_modules){
    descriptions = append(descriptions,unique(summary_tbl$DESCR[summary_tbl$MODULE==j]))
  }
  s_melted_df = data.frame("ID" = standalone_modules,
                           "MODULE"=standalone_modules,
                           "SUPERSET"="none",
                           "DESCR"= descriptions)
  all_melted_df = rbind(all_melted_df, s_melted_df)
  write.table(all_melted_df, paste0(output_Dir, "merged_summary_table.",study,".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
}

#input table from output of summaryTable.
#summary_table - path to table 
#1. merge all modules that were significant for all SGMs
#2. annotate modules - pick one annotation in the end
#3. extract comprising modules for each superset 
#4. go through each module for each SGM in summary_table and reference it to its superset (if any) - add new column that corresponds to its superset
#5. do SGMrepTable function as I've written previously. 
#6. only did for quasi
#set to results folder

keys = c("gene2loci.050kb"=1,"all_48esnps"=2,"all_periph_35esnps"=4, "all_brain_13esnps"=3,"select_24esnps"=5, "periph_18esnps"=6,
         "brain_6esnps"=7,"Adipose_Subcutaneous"=8, "Adipose_Visceral_Omentum"=9, "Adrenal_Gland"=10, "Artery_Aorta"=11,
         "Artery_Tibial"=12,"Brain_Cerebellar_Hemisphere"=13,"Brain_Cerebellum"=14, "Brain_Cortex"=15, "Brain_Frontal_Cortex_BA9"=16, 
         "Brain_Hippocampus"=17, "Brain_Hypothalamus"=18, "Colon_Sigmoid"=19, "Esophagus_Mucosa"=20, "Esophagus_Muscularis"=21, 
         "Heart_Left_Ventricle"=22,"Liver"=23, "Muscle_Skeletal"=24,"Nerve_Tibial"=25, "Pancreas"=26, "Pituitary"=27, "Spleen"=28, 
         "Stomach"=29, "Thyroid"=30, "Whole_Blood"=31)

#module_categories is your merged table
#Results directory must have table_summary_`study`_.txt - case sensitive
SGMrepTable <- function(keys, results_Dir){
 
  #do for both by superset and by module 
  files = list.files(results_Dir)
  table_file = files[grep("table_summary",files)]
  merged_tbl_file = files[grep("merged_summary_table", files)]
  
  #read in summary table 
  summary_tbl = read.delim(paste0(results_Dir,table_file), header=TRUE, stringsAsFactors = FALSE)
  summary_tbl = summary_tbl[!(summary_tbl$SGM==""),]
  summary_tbl = summary_tbl[!(summary_tbl$SGM=="gene2loci.010kb"),]
  summary_tbl = summary_tbl[!(summary_tbl$SGM=="gene2loci.020kb"),]
  
  #read in merging of all modules (module_categories)
  mod_categories = read.delim(paste0(results_Dir,merged_tbl_file), header = TRUE, stringsAsFactors = FALSE)
  #for multiple superset names listed, convert to just showing the first one.
  modify = mod_categories$DESCR[mod_categories$SUPERSET!="none"]
  for(i in 1:length(modify)){
    splitted = unlist(strsplit(modify[i], split = ", "))
    modify[i] = splitted[1]
  }
  mod_categories$DESCR[mod_categories$SUPERSET!="none"] = modify
  
  
  mod_descr = unique(mod_categories$DESCR[mod_categories$SUPERSET!="none"]) #SUPERSETS
  
  combined_supersets = c()
  num_mod = c()
  for(j in mod_descr){
    if(grepl(",", mod_categories$SUPERSET[mod_categories$DESCR==j][1])){
      supersets = which(mod_categories$DESCR==j)
      position = 1
      num_mod = append(num_mod, length(supersets))
      for(k in supersets){
        
          if(mod_categories$SUPERSET[k]=="none"){
            temp[position] = mod_categories$MODULE[k]
            supersets[position] = mod_categories$MODULE[k]
          }
          else{
            supersets[position] = mod_categories$SUPERSET[k]
          }
          position = position + 1
        
      }
      supersets = unique(supersets)
      supersets = paste0(supersets, collapse = ",")
    }
    combined_supersets = append(combined_supersets, supersets)
  }
  
  combined_supersets = append(combined_supersets, mod_categories$MODULE[mod_categories$SUPERSET=="none"]) #INDIVIDUAL
  
  num_mod = append(num_mod, rep(1, length(mod_categories$MODULE[mod_categories$SUPERSET=="none"])))
  
  number = 1
  ids = c()
  for(b in combined_supersets){
    name = paste0("S", number)
    if(grepl(",", b)){
      ids = append(ids, name)
      number = number + 1
    }
    else{
      ids = append(ids, b)
    }
  }
  
  description = append(mod_descr, mod_categories$DESCR[mod_categories$SUPERSET=="none"])
  
  new_df = data.frame("ID" = ids,
                      "MODULE"=combined_supersets,
                      "DESCR"=description) #check if was fine
  
  sign_map = c()
  sign_numMap = c()
  num_map = c()
  
  for(l in description){
    #num_sign = 0
    temp2=c()
    for(m in 1:nrow(summary_tbl)){
      if(l==mod_categories$DESCR[mod_categories$MODULE==summary_tbl$MODULE[m]]){
        if(is.na(keys[summary_tbl$SGM[m]])){
          mapping = summary_tbl$SGM[m]
          key = summary_tbl$SGM[m]
        }
        else{
          mapping = summary_tbl$SGM[m]
          key = keys[summary_tbl$SGM[m]]
        }

        temp2 = paste(temp2, key, sep = ", ") 
        #num_sign = num_sign + 1
      }
    }

    temp2 = substring(temp2,3)
    foo2 = unlist(strsplit(temp2,", "))
    foo2 = unique(foo2)
    foo2 = as.integer(foo2)
    foo2 = sort(foo2)
    num_map = append(num_map, length(foo2))
    
    temp = names(keys[foo2])
    
    temp = paste0(temp, collapse = ", ")
    
    temp2 = paste0(foo2, collapse = ", ")
    
    
    sign_map = append(sign_map, temp)
    sign_numMap = append(sign_numMap, temp2)
    
    #num_mod = append(num_mod, num_sign)
    
  }
  
  new_df$"Significant Mapping Description" = sign_map
  new_df$"Significant Mapping"= sign_numMap
  new_df$"Number of Significant Mappings" = num_map
  new_df$"Number of Significant Modules"= num_mod
  
  new_df_sorted = new_df[order(-new_df$"Number of Significant Mappings"),]
  output_file_name = paste0(results_Dir, "SGMrep_table_", gsub("table_summary_","",table_file))
  write.table(new_df_sorted, output_file_name, row.names = FALSE, quote = FALSE, sep="\t")
  
}

#modified merge function that just takes in a list of modules - already filtered for FDR, no ctrl sets, etc. 
#include last "/" in output_Dir
merge_modules <- function(name,modules_df, output_Dir, modfile_path, infofile_path){
  plan = c()
  plan$folder = output_Dir
  plan$label =  name
  plan$modfile= modfile_path
  plan$inffile= infofile_path
  #=====================================================
  pool=c()
    aa = modules_df
    if (nrow(aa) > 0) pool=unique(c(pool, as.character(aa[,"MODULE"])))
  
  #=====================================================
  #=== Merge the modules before 2nd SSEA
  if (length(pool)>0){
    meg.mods<- tool.read(plan$modfile)
    merged.modules <- pool
    moddata <- meg.mods[which(!is.na(match(meg.mods[,1], merged.modules))),]
    # Merge and trim overlapping modules.
    rmax <- 0.20
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
    mdfile=paste0(name, ".mod.txt"); mifile=paste0(name, ".info.txt")
    
    write.table(moddata, paste0(plan$folder,"merged_", mdfile),  #deleted "/" before "merged"
                sep='\t', col.names=T, row.names=F, quote=F) 
    
    write.table(moddatainfo, paste0(plan$folder,"merged_", mifile),  
                sep='\t', col.names=T, row.names=F, quote=F) 
  }
}

#name includes ("supersets.",study,"." and maybe some other thing)
extract_annotation <- function(study, name, results_Dir){
  file_name = paste0(results_Dir,"supersets.",study,"/","a.supersets.", study, "_Ontology_", name,".xls")
  annot = read.delim(file_name, header = TRUE, stringsAsFactors = FALSE)
  size = c()
  annotation = c()
  if(is.null(annot$ModuleSize[1])){
    annotation = append(annotation,"Unknown")
    size = append(size, "Error in annotation gave no module size")
  }
  else if(is.na(annot$ModuleSize[1])){
    annotation = append(annotation,"Unknown")
    size = append(size, annot$ModuleSize[1])
  }
  else if(sum(annot$pvalue_corrected<0.05)==0){
    annotation = append(annotation,"Unknown")
    size = append(size, annot$ModuleSize[1])
  }
  else{
    sig_mods = annot[annot$pvalue_corrected<0.05,]
    max_mod_overlap_annot = annot$Gene.Category[annot$ModuleOverlap==max(annot$ModuleOverlap)]
    if(length(max_mod_overlap_annot)>1){
      names = c()
      for(k in max_mod_overlap_annot){
        names = paste(names, k, sep = ", ") #?
      }
      names = substring(names, 3)
      annotation = append(annotation, names)
    }
    else{
      annotation = append(annotation,max_mod_overlap_annot)
    }
    size = append(size, annot$ModuleSize[1])
  }
    
    superset_table = data.frame("SUPERSET"= annotation,
                                "MODULE_SIZE" = size, stringsAsFactors = FALSE)
    return(superset_table)
}

#cohort is so that you only get one set of SGMs
#trim is a vector of strings to trim off files to get the "SGM" 
#for example, things to trim include the study or any other common information ex. DIAGRAMstage1_T2D.Adipose_Subcutaneous.50.20.results.txt
#the goal of the trim parameter is to reduce the example to "Adipose_Subcutaneous"
#so your trim vector would be c("DIAGRAMstage1_T2D", ".50.20.results.txt")
#because you just need to get SGMs and can get SGM files using "grep" and this is just an "intersect" or "replication" function, you don't need to put all studies
replication <- function(results_Dir, cohort, study, trim=NULL){
  files = list.files(results_Dir)
  #get SGMs
  subset = files[grep(cohort, files)]
  SGMs = subset
  for(i in trim){
    SGMs = gsub(i, "", SGMs)
  }
  shared = c()
  for(j in SGMs){
    results = files[grep(j, files)]
    iter = 1
    for(k in results){
      if(iter==1){
        df = read.delim(paste0(results_Dir, k), header = TRUE, stringsAsFactors = FALSE)
        shared = df$MODULE
      }
      if(iter>1){
        df = read.delim(paste0(results_Dir, k), header = TRUE, stringsAsFactors = FALSE)
        shared_temp = df$MODULE
        shared = intersect(shared,shared_temp)
      }
      iter = iter + 1 
    }
    shared_df = data.frame("MODULE"=shared)
    write.table(shared_df, paste0(results_Dir, "shared_", study, j, ".txt"), row.names = FALSE, quote = FALSE, sep="\t")
  }
}


move_results <- function(study, results_Dir, coexpr_Dir, output_Dir){
  results = list.files(results_Dir)
  files = results[grep(study, results)]
  files = files[grep("results", files)]
  coexpr_files = list.files(coexpr_Dir)
  for(i in files){
    file = read.delim(paste0(results_Dir, i), header = TRUE, stringsAsFactors = FALSE)
    file = file[file$FDR<=0.05,]
    if(sum(grepl(i,coexpr_files))>0){
      coexpr = coexpr_files[grep(i, coexpr_files)]
      coexpr_df = read.delim(paste0(coexpr_Dir,coexpr), header = TRUE, stringsAsFactors = FALSE)
      coexpr_df = coexpr_df[coexpr_df$FDR<=0.05,]
      combined = rbind(file, coexpr_df)
      combined = combined[!(combined$MODULE=="_ctrlA" | combined$MODULE=="_ctrlB"),]
      write.table(combined, paste0(output_Dir, i), row.names = FALSE, quote = FALSE, sep="\t")
      next
    }
    file = file[!(file$MODULE=="_ctrlA" | file$MODULE=="_ctrlB"),]
    write.table(file, paste0(output_Dir, i), row.names = FALSE, quote = FALSE, sep="\t")
  }
}

# nerve = read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/Data/AD/quasi_overlap/supersets.AD_quasi.Nerve_Tibial.txt", header = TRUE, stringsAsFactors = FALSE)
# nerve = nerve[!(nerve$MODULE=="Nerve_Tibial.GTEXv7.MEGENA_126,.."),]  #this module was not able to be annotated!!!
# write.table(nerve, "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/Data/AD/quasi_overlap/supersets.AD_quasi.Nerve_Tibial1.txt", row.names = FALSE, quote = FALSE, sep = "\t")
#"Nerve_Tibial.GTEXv7.MEGENA_126,.." THIS MODULE WAS NOT ABLE TO BE ANNOTATED IN DIRECT OVERLAP EITHER 
#if you get the error: "Error in cbind(rep(modulesize, dim(fm)[1]), fm) : object 'fm' not found" - FIND OFFENDING MODULE AND TAKE IT OUT

#WORKFLOW FOR GENE OVERLAP (QUASI)
#1. have files for both studies (ex. DIAGRAM and UKB) - combine canonical and coexpression
#             for AD, I put them in "canonical_coexpr". for T2D, I put them in "combined" (FOLDERS)
#2. use geneOverlap function for to produce "shared_" files. all other files are used in this function to produce the "shared_" files in the end and have useful information
#3. transfer "shared_" files to new folder. I called this "quasi_overlap"
#4. merge modules in "shared_" files
#5. annotate supersets resulting from step 4 using annotate_supersets function
#6. do summaryTable function to get summary table! 

# setwd("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/Data/AD_T2D/shared_quasi/")
# for(i in files){
#   file = read.delim(i, header = TRUE, stringsAsFactors = FALSE)
#   write.table(file, paste0("AD.",i), row.names = FALSE, quote = FALSE, sep = "\t")
# }
# 
# files = list.files("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/Data/AD_T2D/shared_quasi/")
# subset = files[grep("AD",files)]
# SGMs = gsub("AD.shared_", "", gsub(".50.20.txt", "", subset))
# 
# #T2D does not have colon sigmoid or pituitary
# #AD does not have brain_6esnps, heart left ventricle, or liver
# 
# 
# setwd("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/Data/AD_T2D/shared_quasi/")
# for(i in SGMs){
#   shared_files = files[grep(i, files)]
#   if(length(shared_files)==1){
#     next
#   }
#   first = read.delim(shared_files[1], header = TRUE, stringsAsFactors = FALSE)
#   second = read.delim(shared_files[2], header = TRUE, stringsAsFactors = FALSE)
#   shared_modules = intersect(first$MODULE, second$MODULE)
#   shared_modules_df = data.frame("MODULE"= shared_modules, stringsAsFactors = FALSE)
#   write.table(shared_modules_df, paste0("intersect_AD_T2D.",i,".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
# }

#took out Esophagus_Muscularis.GTEXv7.MEGENA_168,..