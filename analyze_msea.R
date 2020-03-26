# Please excuse the roughness of these functions! I made these a long time ago while I was learning R and I will be making
# efforts to improve the code hopefully soon.

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
# Written by Jessica Ding                                                           #
#                                                                                           #
# ------------------------------------------------------------------------------------------#


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
#         5) cohorts             vector with study names that are reflected in the msea results
#                                ex. c("DIAGRAMstage1_T2D", "UKB_T2D")
#         6) modfile_path        path to mod file, columns: 'MODULE' 'GENE' 
#         7) infofile_path       path to info file, columns: 'MODULE' 'GENE' 'DESCR'
#         8) study               
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

geneOverlap <- function(results_Dir, output_Dir, perc_gene_overlap, fdr_cutoff=0.05, cohorts, modfile_path, infofile_path, study, type=NULL){
  
  # MAKE COMBINED FILES FOR ALL STUDIES 
  files = c()
  for(c in cohorts){
    files = append(files, list.files(results_Dir)[grep(c, list.files(results_Dir))])
  }
  
  ##get SGM methods. This assumes your results file is formatted - "cohort.SGM.results.txt" 
  study_subset = files[grep(cohorts[1], files)]
  SGMs = c()
  for(f in study_subset){
    SGMs = append(SGMs, unlist(strsplit(f, split = ".",fixed = TRUE))[2]) 
  }
  SGMs = unique(SGMs)
  cat("The SNP to gene mappings used is/are:\n")
  print(SGMs) #this is to check that the SGMs are being extracted in the desired way
  
  ##combine studies for each SGM
  cat("Now combining studies together in one file for each SGM.\n")
  for(i in SGMs){
    results = files[grep(i,files)]
    if(length(results)<2){
      cat("The method ", i, " does not have a file for one of the studies.\n")
    }
    all = data.frame()
    for(j in results){
      temp = read.delim(paste(results_Dir, j ,sep = ""), header = TRUE, stringsAsFactors = FALSE)
      if(is.null(type)){
        temp = temp[temp$FDR<fdr_cutoff,]
      }
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
  # MERGE MODULES - this was modified 5.19.2019 and not currently tested 
  for(l in combined_files){
    cat("Now merging:\n")
    cat(l,"\n")
    # unfortunately hard coding this in...
    if(study=="UKB_Wojcik" & grepl("Pituitary_sQTL",l)) next 
    df = read.delim(paste0(output_Dir,l), header = TRUE, stringsAsFactors = FALSE)
    if(sum(is.na(df$MODULE))>0){
      cat("Skipping: ", l, " because there is no significance.\n")
      next
    }
    if(length(df$MODULE)<2){
      cat("Skipping: ", l, " because there is only one or 0 significant modules meaning no overlap at all.\n")
      next
    }
    merge_modules(name = gsub("combined_","",gsub(".txt","",l)), 
                  modules_df = df, 
                  rcutoff = 0.30, 
                  output_Dir = output_Dir, 
                  modfile_path = modfile_path, 
                  infofile_path = infofile_path)
    }
  
  cat("Now analyzing for newly overlapped modules...\n")
  # FIND IF THERE WERE MERGED MODULES WITH DIFFERENT STUDY REPRESENTATION
  files1 = list.files(output_Dir)
  merged_mod_files = files1[grep("merged_",files1)] #works if no SGMs have "mod" in name
  merged_mod_files = merged_mod_files[grep(".mod", merged_mod_files)]

  uncommon_df = data.frame(stringsAsFactors = FALSE)
  for(m in merged_mod_files){
    cat("Now inspecting ", m , "\n")
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
    if(sum(is.na(comb_file$MODULE))>0){
      cat("Skipping: ", l, " because there is no significance.\n")
      next
    }
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
    shared = data.frame(stringsAsFactors = FALSE)
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

#-----------------------------------------------------------------------------------------------------#
# This function makes annotations for supersets that were obtained from module merge. 
#
# Inputs: 1) results_Dir         folder with ".mod.txt" files (output from module merge)
#         2) ontologyfname       path to ontology file
#         3) study               name to append to files 
#         4) trim                vector of strings to "trim" off the mod.txt file name to get desired unique name
#                                ex. UKB_AD.Adipose_Subcutaneous.mod.txt --> trim = c("UKB_AD.",".mod.txt")
#
# Outputs: 1) "supersets.study" folder     contains superset annotations in .xls format
#
# Notes: If you get the error "Error in cbind(rep(modulesize, dim(fm)[1]), fm) : object 'fm' not found", the 
#        offending module is printed before this error and this module must be taken out from the .mod.txt file 

# Version April 2019
#---------------------------------------------------------------------------------------------------#
# ontologyfname = "~/Desktop/Yang_Lab/resources/genesets/MSigDB_Canonical_Pathways.txt"
#--------------------------------------------------------------------------#
#
annotate_supersets <- function(results_Dir, ontologyfname, study="", trim){ #make trim argument that has a vector of things to trim to get the "SGM" and also so that you don't wind up with .mod.txt in subseqeunt result file names
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
   if(grepl("quasi", study)){ 
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
  for(i in shared_files){
    df = data.frame(stringsAsFactors = FALSE)
    file = read.delim(paste0(results_Dir, i), header = TRUE, stringsAsFactors = FALSE)
    if(length(file$MODULE)==1){
      SGM = gsub(paste0("shared_", study, "."),"",i)
      SGM = gsub(".txt","", SGM)
      if(grepl("quasi",study)){
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
                        "MODULE_SIZE"="value", stringsAsFactors = FALSE)
        
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
        annot = read.delim(paste("~/Desktop/Yang_Lab/resources/genesets/reannotate/trim_coexpr/a.trim_coexpr_Ontology_", j, ".xls", sep = ""), header = TRUE, stringsAsFactors = FALSE)
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


#set output_Dir to results folder 
mergedSummaryTable <- function(summary_table, study, output_Dir, modfile_path, infofile_path, ontologyfname, rcutoff){
  mod_file = read.delim(modfile_path, header = TRUE, stringsAsFactors = FALSE)
  summary_tbl = read.delim(paste0(output_Dir,summary_table), header=TRUE, stringsAsFactors = FALSE)
  modules = summary_tbl$MODULE[!(summary_tbl$MODULE=="")]
  modules = unique(modules)
  modules_df = data.frame("MODULE"=modules)
  cat("The study being analyzed is", study, "\n")
  merge_modules(name = study, modules_df = modules_df, output_Dir = output_Dir, modfile_path = modfile_path, 
                infofile_path = infofile_path, rcutoff = rcutoff)
  #annotate supersets
  annotate_supersets(results_Dir = output_Dir, ontologyfname = ontologyfname, study = study, trim = c("merged_", ".mod.txt"))
  merged = read.delim(paste0(output_Dir, "merged_",study, ".mod.txt"), header = TRUE, stringsAsFactors = FALSE)
  supersets = unique(merged$OVERLAP[grep(",", merged$OVERLAP)])
  all_melted_df = data.frame(stringsAsFactors = FALSE)
  number = 1
  index = c()
  aggregated_descr = c()
  for(i in supersets){ #maybe make function for this? break down supersets into individual modules
    name = unique(merged$MODULE[merged$OVERLAP==i])
    sub_modules = unlist(strsplit(i, split = ","))
    descriptions = c()
    sizes = c()
    
    for(k in sub_modules){
      descriptions = append(descriptions,unique(summary_tbl$DESCR[summary_tbl$MODULE==k]))
      sizes = append(sizes, sum(mod_file$MODULE==k))
    }
    
    temp = c()
    for(c in descriptions){
      temp = paste(temp, c, sep = ", ")
    }
    temp = substring(temp, 3)
    aggregated_descr = append(aggregated_descr, temp) #for renaming
    melted_df = data.frame("ID" = paste0("S",number),
                           "MODULE"=sub_modules,
                           "DESCR" = descriptions,
                           "INDIVID_SIZE" = sizes,
                           "SUPERSET_SIZE" = extract_annotation(study = study, name = name, results_Dir = output_Dir)[,2],
                           "SUPERSET"=i,
                           "SUPERSET_DESCR"=extract_annotation(study = study, name = name, results_Dir = output_Dir)[,1], stringsAsFactors = FALSE)
    
    all_melted_df = rbind(all_melted_df, melted_df)
    index = append(index, nrow(all_melted_df))
    number = number + 1
  }
  superset_length = nrow(all_melted_df)
  standalone_modules = unique(merged$MODULE[!grepl(",", merged$MODULE)])
  descriptions = c()
  sizes = c()
  for(j in standalone_modules){
    descriptions = append(descriptions,unique(summary_tbl$DESCR[summary_tbl$MODULE==j]))
    sizes = append(sizes, sum(mod_file$MODULE==j))
  }
  s_melted_df = data.frame("ID" = standalone_modules,
                           "MODULE"=standalone_modules,
                           "DESCR" = descriptions,
                           "INDIVID_SIZE" = sizes,
                           "SUPERSET_SIZE"= "none",
                           "SUPERSET"="none",
                           "SUPERSET_DESCR"= "none")
  
  all_melted_df = rbind(all_melted_df, s_melted_df)
  index = append(index, seq(from = (superset_length+1), to = nrow(all_melted_df), by = 1))
  
  
  # ------- For further refined downstream analysis ------- #
  aggregated_descr = append(aggregated_descr, descriptions)
  renamed = data.frame("ID"= all_melted_df$ID[index],
                       "MODULE" = all_melted_df$SUPERSET[index],
                       "DESCR" = aggregated_descr,
                       "RENAME" = "",
                       "ANNOTATION" = all_melted_df$SUPERSET_DESCR[index])
  write.table(all_melted_df, paste0(output_Dir, "30merged_summary_table.",study,".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(renamed, paste0(output_Dir, "30merged_summary_table_rename.",study,".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
}

summaryTablev2 <- function(merged_table, merged_table_rename, summary_table, output_Dir){
  merged_tbl = read.delim(merged_table, header = TRUE, stringsAsFactors = FALSE)
  merged_tbl_rename = read.delim(merged_table_rename, header = TRUE, stringsAsFactors = FALSE)
  summary_tbl = read.delim(summary_table, header = TRUE, stringsAsFactors = FALSE)
  
  modified_summary_tbl = summary_tbl
  
  # replace SUPERSET names from summary table with merged_tbl SUPERSET
  modules = summary_tbl$MODULE
  
  renames = c()
  ids = c()
  
  for(m in modules){
    if(m==""){
      renames = append(renames, "")
      ids = append(ids, "")
      next
    }
    if(merged_tbl$SUPERSET[merged_tbl$MODULE==m]=="none"){
      if(sum(grepl(m, merged_tbl_rename$ID))==0){
        renames = append(renames, merged_tbl_rename$RENAME[grep(m, merged_tbl_rename$MODULE)])
        ids = append(ids, merged_tbl_rename$ID[grep(m, merged_tbl_rename$MODULE)])
        next
      }
      if(merged_tbl_rename$RENAME[merged_tbl_rename$ID==m]!=""){
        renames = append(renames, merged_tbl_rename$RENAME[merged_tbl_rename$ID==m])
        ids = append(ids, "none")
        next
      }
      renames = append(renames, "none")
      ids = append(ids, "none")
      next
    }
    # get superset ID
    ID = merged_tbl$ID[merged_tbl$MODULE==m]
    ids = append(ids, ID)
    renames = append(renames, merged_tbl_rename$RENAME[merged_tbl_rename$ID==ID])
  }
  
  modified_summary_tbl$SUPERSET_ID = ids
  modified_summary_tbl$RENAME = renames
  
  modified_summary_tbl = modified_summary_tbl[,c("SGM","MODULE","DESCR","SUPERSET_ID","RENAME","MODULE_SIZE","OVERLAP_METHOD")]
  write.table(modified_summary_tbl, paste0(output_Dir,"combined.renamed_summary_table.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
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
SGMrepTable <- function(keys, results_Dir, infofile_path){
 
  #do for both by superset and by module 
  files = list.files(results_Dir)
  table_file = files[grep("renamed_table_summary",files)]
  merged_tbl_file = files[grep("merged_summary_table.AD_quasi", files)]
  
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
            supersets[position] = mod_categories$MODULE[k] # necessarily????
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
                      "DESCR"=description, stringsAsFactors = FALSE) #check if was fine
  
  sign_map = c()
  sign_numMap = c()
  num_map = c()
  
  for(l in description){
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
    
  }
  
  new_df$"Significant Mapping Description" = sign_map
  new_df$"Significant Mapping"= sign_numMap
  new_df$"Number of Significant Mappings" = num_map
  new_df$"Number of Significant Modules"= num_mod
  
  # fix problem of mapping to unknown modules for every module that had no annotation. might run into problems that
  # merged modules with no annotation will have "1" Significant modules when in fact it was a superset with multiple modules
  unknown = new_df[new_df$DESCR=="Unknown",]
  if(length(unknown$MODULE)>0){
    for(d in 1:nrow(unknown)){
      for(e in names(keys)){
        if(grepl(e, unknown$MODULE[d])){
          unknown$"Significant Mapping Description"[d]=e
          unknown$"Significant Mapping"[d] = keys[e]
          unknown$"Number of Significant Mappings"[d] = 1
          unknown$"Number of Significant Modules"[d]= 1
        }
      }
    }
    new_df = new_df[!(new_df$DESCR=="Unknown"),]
    new_df = rbind(new_df, unknown)
  }
  
  
  new_df_sorted = new_df[order(-new_df$"Number of Significant Mappings"),]
  study = gsub("table_summary_","",table_file)
  output_file_name = paste0(results_Dir, "SGMrep_table_", study)
  write.table(new_df_sorted, output_file_name, row.names = FALSE, quote = FALSE, sep="\t")
  
  info = read.delim(infofile_path, header = TRUE, stringsAsFactors = FALSE)
  #Create superset information table
  superset_ids = ids[1:length(mod_descr)]
  superset_info = data.frame(stringsAsFactors = FALSE)
  for(a in 1:length(mod_descr)){
    sub_modules = unlist(strsplit(new_df$MODULE[new_df$DESCR==mod_descr[a]], split = ","))
    temp_df = data.frame("SUPERSET_ID"=superset_ids[a],
                         "MODULE"=sub_modules, stringsAsFactors = FALSE,
                         "SUPERSET_DESCR"= mod_descr[a])
    superset_info = rbind(superset_info, temp_df)
  }
  individual_df = data.frame("SUPERSET_ID"="none",
                             "MODULE"=mod_categories$MODULE[mod_categories$SUPERSET=="none"],
                             "SUPERSET_DESCR"="none")
  all_mod_info = rbind(superset_info, individual_df)
  
  sources = c()
  for(b in all_mod_info$MODULE){
    sources = append(sources, info$SOURCE[info$MODULE==b])
  }
  all_mod_info$SOURCE = sources
  
  descriptions = c()
  for(c in all_mod_info$MODULE){
    descriptions = append(descriptions, info$DESCR[info$MODULE==c])
  }
  
  all_mod_info$DESCR = descriptions
  
  all_mod_info = all_mod_info[, c(1,2,4,5,3)]
  
  write.table(all_mod_info, paste0(results_Dir, "mod.info_",study), row.names = FALSE, quote = FALSE, sep="\t")
}

SGMrepTablev2 <- function(keys, table_summary, results_Dir, study){
  
  tbl_summary = read.delim(paste0(results_Dir, table_summary), header = TRUE, stringsAsFactors = FALSE)
  supersets = unique(tbl_summary$SUPERSET_ID[grep("S",tbl_summary$SUPERSET_ID)])
  all_modules = append(supersets, unique(tbl_summary$MODULE[tbl_summary$SUPERSET=="none"]))
  
  mapping_names = c()
  key_map = c()
  num_map = c()
  descriptions = c()
  superset_size = c()
  
  for(l in all_modules){
    temp = c()
    temp2 = c()
    for(m in 1:nrow(tbl_summary)){
      if(sum(supersets==l)>0){
        if(tbl_summary$SUPERSET_ID[m]==l){
          if(is.na(keys[tbl_summary$SGM[m]])){
            key = tbl_summary$SGM[m]
          }
          else{
            key = keys[tbl_summary$SGM[m]]
          }
          
          temp = paste(temp, key, sep = ", ")
      }
      }
      else {
        if(tbl_summary$MODULE[m]==l){
          if(is.na(keys[tbl_summary$SGM[m]])){
            key = tbl_summary$SGM[m]
          }
          else{
            key = keys[tbl_summary$SGM[m]]
          }
          
          temp = paste(temp, key, sep = ", ")

        }
      }
    }
    if(sum(supersets==l)>0){
      descriptions = append(descriptions, tbl_summary$RENAME[tbl_summary$SUPERSET_ID==l][1])
    }
    else if(unique(tbl_summary$RENAME[tbl_summary$MODULE==l])!="none"){
      descriptions = append(descriptions, tbl_summary$RENAME[tbl_summary$MODULE==l][1])
    }
    else{
      descriptions = append(descriptions, tbl_summary$DESCR[tbl_summary$MODULE==l][1])
    }
    
    if(sum(supersets==l)>0){
      superset_size = append(superset_size, length(unique(tbl_summary$MODULE[tbl_summary$SUPERSET_ID==l])))
    }
    else{
      superset_size = append(superset_size, 1)
    }
    
    temp = substring(temp,3)
    foo = unlist(strsplit(temp,", "))
    foo = unique(foo)
    foo = as.integer(foo)
    foo = sort(foo)
    num_map = append(num_map, length(foo))
    
    temp2 = names(keys[foo]) # keep keys as vector
    
    temp2 = paste0(temp2, collapse = ", ")
    
    temp = paste0(foo, collapse = ", ") # now conjoin keys
    
    
    key_map = append(key_map, temp) 
    mapping_names = append(mapping_names, temp2)
    
  }
  new_df = data.frame("MODULE"= all_modules,
                      "DESCR" = descriptions,
                      "Significant Mapping Description" = mapping_names,
                      "Significant Mapping" = key_map, 
                      "Number_of_Significant_Mappings" = num_map,
                      "Number of Modules" = superset_size,
                      stringsAsFactors = FALSE)
  
  
  new_df_sorted = new_df[order(-new_df$Number_of_Significant_Mappings),]
  unknowns = new_df_sorted[new_df_sorted$DESCR=="Unknown",]
  new_df_sorted = new_df_sorted[!(new_df_sorted$DESCR=="Unknown"),]
  new_df_sorted = rbind(new_df_sorted, unknowns)
  
  # study = gsub("30renamed_table_summary_","",table_summary)
  output_file_name = paste0(results_Dir, "combined.SGMrep_table_", study,".txt")
  write.table(new_df_sorted, output_file_name, row.names = FALSE, quote = FALSE, sep="\t")
  
}

combineSGMrepTable <- function(results1, results2){
  first = read.delim(results1, header = TRUE, stringsAsFactors = FALSE)
  second = read.delim(results2, header = TRUE, stringsAsFactors = FALSE)
   
}

#modified merge function that just takes in a list of modules - already filtered for FDR, no ctrl sets, etc. 
#include last "/" in output_Dir
merge_modules <- function(name, modules_df, rcutoff, output_Dir, modfile_path, infofile_path){
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
    rmax <- rcutoff
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
    #max_mod_overlap_annot = annot$Gene.Category[annot$ModuleOverlap==max(annot$ModuleOverlap)] #before wanted to get highest overlap
    max_mod_overlap_annot = sig_mods$Gene.Category
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
    #return(annotation)
}

# cohort is so that you only get one set of SGMs
# trim is a vector of strings to trim off files to get the "SGM" 
# for example, things to trim include the study or any other common information ex. DIAGRAMstage1_T2D.Adipose_Subcutaneous.50.20.results.txt
# the goal of the trim parameter is to reduce the example to "Adipose_Subcutaneous"
# so your trim vector would be c("DIAGRAMstage1_T2D", ".50.20.results.txt")
# because you just need to get SGMs and can get SGM files using "grep" and this is just an "intersect" or "replication" function, you don't need to put all studies
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
    write.table(shared_df, paste0(results_Dir, "shared_", study,".", j, ".txt"), row.names = FALSE, quote = FALSE, sep="\t")
  }
}


combine_results <- function(folders, output_Dir){
  # just to get traits and mapping
  files = list.files(folders[1]) 
  mapping = c()
  traits = c()
  for(f in files){
    mapping = append(mapping, unlist(strsplit(f, split = ".",fixed = TRUE))[2]) 
    traits = append(traits, unlist(strsplit(f, split = ".",fixed = TRUE))[1]) 
  }
  mapping = unique(mapping)
  traits = unique(traits)
  for(t in traits){
    for(m in mapping){
      combined = data.frame(stringsAsFactors = FALSE)
      files = c()
      for(f in folders){
        file = list.files(f, full.names = TRUE)[grep(t, list.files(f))] 
        file = file[grep(m, file)]
        # should just get one file from each folder 
        if(length(file)>1){
          cat("More than one similar file in ", f, "\n")
        }
        if(length(file)==0){
          next
        }
        files = append(files, file)
      }
      if(is.null(files)) next
      for(f in files){
        df <- read.delim(f, stringsAsFactors = FALSE)
        combined = rbind(combined, df)
      }
      combined = combined[!(combined$MODULE=="_ctrlA" | combined$MODULE=="_ctrlB"),]
      # will output the same name
      write.table(combined, paste0(output_Dir, tail(unlist(strsplit(files[1], split = "/")),n=1)), 
                  row.names = FALSE, quote = FALSE, sep="\t")
    }
  }
  
  # results = list.files(results_Dir)
  # files = results[grep(study, results)]
  # files = files[grep("results", files)]
  # coexpr_files = list.files(coexpr_Dir)
  # for(i in files){
  #   file = read.delim(paste0(results_Dir, i), header = TRUE, stringsAsFactors = FALSE)
  #   file = file[file$FDR<=0.05,]
  #   if(sum(grepl(i,coexpr_files))>0){
  #     coexpr = coexpr_files[grep(i, coexpr_files)]
  #     coexpr_df = read.delim(paste0(coexpr_Dir,coexpr), header = TRUE, stringsAsFactors = FALSE)
  #     coexpr_df = coexpr_df[coexpr_df$FDR<=0.05,]
  #     combined = rbind(file, coexpr_df)
  #     combined = combined[!(combined$MODULE=="_ctrlA" | combined$MODULE=="_ctrlB"),]
  #     write.table(combined, paste0(output_Dir, i), row.names = FALSE, quote = FALSE, sep="\t")
  #     next
  #   }
  #   file = file[!(file$MODULE=="_ctrlA" | file$MODULE=="_ctrlB"),]
  #   write.table(file, paste0(output_Dir, i), row.names = FALSE, quote = FALSE, sep="\t")
  # }
}

combine_details <- function(canon_Dir, coexp_Dir, output_Dir){
  files = list.files(canon_Dir)
  for(file in files){
    if(sum(list.files(coexp_Dir)==file)==0){
      mapping = file.copy(from = paste0(canon_Dir, file), to = paste0(output_Dir,file))
    }
    else{
      canon = read.delim(paste0(canon_Dir, file), stringsAsFactors = FALSE)
      coexp = read.delim(paste0(coexp_Dir, file), stringsAsFactors = FALSE)
      both = rbind(canon, coexp)
      write.table(both, paste0(output_Dir, file), row.names = FALSE, quote = FALSE, sep="\t")
    }
  }
}

runKDA <- function(nodes, network, trim = NULL){
  cat("\nNow analyzing:", gsub(".txt","",nodes), "with", gsub(".txt","",network), "\n")
  job.kda <- list()
  job.kda$label<-"wKDA"
  network_name = gsub(".txt","",unlist(strsplit(network, split = "/"))[length(unlist(strsplit(network, split = "/")))])
  nodes_name = nodes
  for(i in trim){
    nodes_name = gsub(i,"",nodes_name)
  }
  name = paste(nodes_name, network_name,sep = "_")
  
  job.kda$folder<- name ## parent folder for results
  
  ## Input a network
  ## columns: TAIL HEAD WEIGHT
  job.kda$netfile <- network
  ## Input gene list
  ## columns: MODULE NODE
  job.kda$modfile <- nodes
  
  ## "0" means we do not consider edge weights while 1 is opposite.
  job.kda$edgefactor <- 0
  ## The searching depth for the KDA
  job.kda$depth <- 1
  ## 0 means we do not consider the directions of the regulatory interactions
  ## while 1 is opposite.
  job.kda$direction <- 1
  job.kda$nperm <- 10000
  
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
  
  job.kda <- kda2cytoscape(job.kda)
}

summarizeKDA <- function(KDA_folder, protein_descrip, name){
  results = data.frame(stringsAsFactors = FALSE)
  KDs = data.frame(stringsAsFactors = FALSE)
  protein_names = read.delim(protein_descrip, 
                             header = TRUE, stringsAsFactors = FALSE, quote = "")
  for(l in list.files(KDA_folder)){ 
    if(length(list.files(paste0("./",l, "/")))==1){
      cat("No significant key drivers found for", l, "\n")
      next
    }
    else{
      key_drivers = read.delim(paste0(KDA_folder,l,"/kda/wKDA.results.txt"), stringsAsFactors = FALSE)
      valued_indices = which(!is.na(key_drivers$NODE[1:5]))
      edges = read.delim(paste0(KDA_folder,l,"/cytoscape/kda2cytoscape.edges.txt"), stringsAsFactors = FALSE)
      edges = edges[!duplicated(edges[,c(1,2)]),] # makes an edge for each module but don't want this..
      nodes = read.delim(paste0(KDA_folder,l,"/cytoscape/kda2cytoscape.nodes.txt"), stringsAsFactors = FALSE)
      temp = data.frame(stringsAsFactors = FALSE)
      
      # make network summary table
      temp = data.frame("NETWORK"=l,
                        "KDs"=length(key_drivers$NODE),
                        "KDs_p<0.05"=sum(key_drivers$P<0.05),
                        "KDs_fdr<0.05"=sum(key_drivers$FDR<0.05),
                        "topKDs"=concatenate(key_drivers$NODE[valued_indices], mysep = ", "),
                        "topKDs_fdr"=concatenate(key_drivers$FDR[valued_indices], mysep = ", "),
                        "n_nodes" = length(nodes$NODE),
                        "n_edges"=length(edges$TAIL),
                        "avg_degree"= (length(edges$TAIL)/length(unique(append(edges$TAIL, edges$HEAD)))),
                        "perc_member"=(sum(grepl("chart", nodes$URL)))/length(nodes$NODE), stringsAsFactors = FALSE)
      results = rbind(results, temp)
      
      info = c()
      
      indices = which(!is.na(key_drivers$NODE[key_drivers$FDR<0.05]))
      for(m in key_drivers$NODE[indices]){
        if(length(protein_names$annotation[protein_names$preferred_name==m])==0){
          info = append(info, "Not annotated")
        }
        else if(length(protein_names$annotation[protein_names$preferred_name==m])>1){
          info = append(info, 
                        concatenate(protein_names$annotation[protein_names$preferred_name==m],
                                    mysep = ","))
          cat("This protein had more than one annotation:", m, "\n")
        }
        else{
          info = append(info, protein_names$annotation[protein_names$preferred_name==m])
        }
      }
      
      # make KD summary table
      kd_temp = data.frame("NETWORK"=l,
                           "KDs"= key_drivers$NODE[indices],
                           "Degrees" = key_drivers$N.neigh[indices],
                           "P-value" = key_drivers$P[indices],
                           "FDR" = key_drivers$FDR[indices],
                           "INFO" = info,
                           "MEMBER" = key_drivers$MEMBER[indices], stringsAsFactors = FALSE)
      KDs = rbind(KDs, kd_temp)
    }
  }
  
  write.table(results, paste0(name, "_results_summary.txt"), 
              row.names = FALSE, quote = FALSE, sep="\t")
  
  write.table(KDs, paste0(name, "_KDs_summary.txt"), 
              row.names = FALSE, quote = FALSE, sep="\t")
}

# this function outputs a new network in the folder with beginning "mod_"
# this function has repetitive code.. will make better later
# when modify nodes file to add GWAS information, have to delete, overwrite, or put in another folder the original nodes file
trim_network <- function(folder,keep_only_orig_input = FALSE, keep_only_GWAS_hits = FALSE){
  edges_file = list.files(folder)[grep("edges",list.files(folder))]
  edges = read.delim(paste0(folder,edges_file), stringsAsFactors = FALSE)
  nodes = read.delim(paste0(folder, list.files(folder)[grep("nodes",list.files(folder))]), stringsAsFactors = FALSE)
  prefixes = c("^MMT", "^NM_","^XM_","ENSMUST","^ri", "_at")
  # AK genes are long noncoding RNAs
  for(pre in prefixes){ 
    edges = edges[!grepl(pre,edges$HEAD),]
  }
  if(keep_only_orig_input & !keep_only_GWAS_hits){
    # keep only KDs and nodes with URLs
    # make list of KDs and URL nodes 
    # keep edges that connect KD to member node
    keep_nodes = c()
    for(i in nodes$NODE){
      if(nodes$SHAPE[nodes$NODE==i]!="Diamond" & nodes$URL[nodes$NODE==i]!=""){
        keep_nodes = append(keep_nodes, i)
      }
    }
    total = data.frame(stringsAsFactors = FALSE)
    for(r in keep_nodes){
      temp = edges[grep(paste0("^",r,"$"), edges$HEAD ),]# KDs in the TAIL column - keep KDs
      total = rbind(total, temp)
    }
    total = total[!duplicated(total),]
    total = total[!duplicated(total[,c(1,2)]),]
    
    write.table(total, paste0(folder, "mod_",edges_file), row.names = FALSE, quote = FALSE, sep="\t")
    return(total)
  }
  # should do a version that filters ONLY gwas hits so that genes that are not members can be included
  # though this happens rarely
  else if(keep_only_orig_input & keep_only_GWAS_hits){ # must have GWAS_hit_meta column, further filtering after filtering by module
    # keep only KDs and nodes with URLs
    # make list of KDs and URL nodes 
    # keep edges that connect KD to member node
    keep_nodes = c()
    for(i in nodes$NODE){
      if(nodes$SHAPE[nodes$NODE==i]!="Diamond" & nodes$URL[nodes$NODE==i]!="" & nodes$GWAS_hit_meta[nodes$NODE==i]!="no"){
        keep_nodes = append(keep_nodes, i)
      }
    }
    total = data.frame(stringsAsFactors = FALSE)
    for(r in keep_nodes){
      temp = edges[grep(paste0("^",r,"$"), edges$HEAD ),]# KDs in the TAIL column - keep KDs
      total = rbind(total, temp)
    }
    total = total[!duplicated(total),]
    total = total[!duplicated(total[,c(1,2)]),]
    
    write.table(total, paste0(folder, "mod_",edges_file), row.names = FALSE, quote = FALSE, sep="\t")
    return(total)
  }
  else if(!keep_only_orig_input & keep_only_GWAS_hits){
    keep_nodes = c()
    for(i in nodes$NODE){
      if(nodes$SHAPE[nodes$NODE==i]!="Diamond" & nodes$GWAS_hit_meta[nodes$NODE==i]!="no"){
        keep_nodes = append(keep_nodes, i)
      }
    }
    total = data.frame(stringsAsFactors = FALSE)
    for(r in keep_nodes){
      temp = edges[grep(paste0("^",r,"$"), edges$HEAD ),]# KDs in the TAIL column - keep KDs
      total = rbind(total, temp)
    }
    total = total[!duplicated(total),]
    total = total[!duplicated(total[,c(1,2)]),]
    
    write.table(total, paste0(folder, "mod_",edges_file), row.names = FALSE, quote = FALSE, sep="\t")
    return(total)
  }
  else{
    write.table(edges, paste0(folder, "mod_",edges_file), row.names = FALSE, quote = FALSE, sep="\t")
    return(edges)
  }
}

# would be better to make reference that says all the genes that could be mapped from the significant values, not just for
# individual mapping methods
# also want to add the info in it was in both AD and T2D
# ^ can do after this
# do same nodes file on AD and T2D, then combine the T2D column result with AD, do if statements to make new column stating 
# whether there was any "overlap"
add_GWAS_hit_info_network <- function(nodes, SGM, genes_folder, study){
  nodes = read.delim(nodes, stringsAsFactors = FALSE)
  files = list.files(genes_folder)
  file = files[intersect(grep("noDescrip",files), grep(SGM,files))]
  df = read.delim(paste0(genes_folder,file), stringsAsFactors = FALSE)
  GWAS_hit = c()
  for(i in nodes$NODE){
    if(is.character(df$Shared.genes)){
      if(sum(unlist(strsplit(df$Shared.genes, split = ", "))==i)>0){
        GWAS_hit = append(GWAS_hit, paste0("Replicated_GWAS_",study))
        next
      }
      else if(sum(unlist(strsplit(df[,5], split = ", "))==i)>0){ 
        GWAS_hit = append(GWAS_hit, paste0("GWAS_hit_", gsub("_genes","",colnames(df)[5])))
      }
      else if(sum(unlist(strsplit(df[,6], split = ", "))==i)>0){
        GWAS_hit = append(GWAS_hit, paste0("GWAS_hit_", gsub("_genes","",colnames(df)[6])))
      }
      else{
        GWAS_hit = append(GWAS_hit, "no")
      }
    }
    else{
      GWAS_hit = append(GWAS_hit, "no")
    }
  }
  nodes$GWAS_hit = GWAS_hit
  return(nodes)
}


extract_annotationv2 <- function(module, set, results_Dir){
  file_name = paste0(results_Dir, set, module,".xls")
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
    #max_mod_overlap_annot = annot$Gene.Category[annot$ModuleOverlap==max(annot$ModuleOverlap)] #before wanted to get highest overlap
    max_mod_overlap_annot = sig_mods$Gene.Category
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
  #return(annotation)
}

trimSummaryTable <- function(summary_table, study, output_Dir){
  new_summary_tbl = data.frame(stringsAsFactors = FALSE)
  summary_tbl = read.delim(summary_table, header = TRUE, stringsAsFactors = FALSE)
  SGMs = unique(summary_tbl$SGM)
  SGMs = SGMs[SGMs!=""]
  for(i in SGMs){
    all_n_subset = data.frame(stringsAsFactors = FALSE)
    SGM_subset = summary_tbl[summary_tbl$SGM==i,]
    supers = SGM_subset$SUPERSET_ID[SGM_subset$SUPERSET_ID!="none"]
    supersets = unique(supers)
      for(j in supersets){
        subset = SGM_subset[SGM_subset$SUPERSET_ID==j,]
        n_subset = data.frame("SGM"=i, 
                              "MODULE"=concatenate(subset$MODULE, mysep = ", "),
                              "DESCR"=concatenate(subset$DESCR, mysep = ", "),
                              "SUPERSET_ID"=subset$SUPERSET_ID[1],
                              "RENAME"=subset$RENAME[1],
                              "MODULE_SIZE"=subset$MODULE_SIZE[1],
                              "OVERLAP_METHOD"=concatenate(subset$OVERLAP_METHOD, mysep = ", "), stringsAsFactors = FALSE)
        all_n_subset = rbind(all_n_subset, n_subset)  
      }
    all_n_subset = rbind(all_n_subset, SGM_subset[SGM_subset$SUPERSET_ID=="none",])
    new_summary_tbl = rbind(new_summary_tbl, all_n_subset)
    new_summary_tbl = rbind(new_summary_tbl,"")
  }
  write.table(new_summary_tbl, paste0(output_Dir, "trimmed_summary_table_",study,".txt"), row.names = FALSE, quote = FALSE, sep="\t")
}

enrichSupersetTable <- function(name1, sgm_rep_table1, name2, sgm_rep_table2, superset_info, output_Dir){
  
  stbl_1 = read.delim(sgm_rep_table1, header = TRUE, stringsAsFactors = FALSE)
  stbl_2 = read.delim(sgm_rep_table2, header = TRUE, stringsAsFactors = FALSE)
  supersets = read.delim(superset_info, header = TRUE, stringsAsFactors = FALSE)
  
  modules = supersets$MODULE
  
  represent = c()
  # mapping1 = c()
  # mapping2 = c()
  
  for(i in modules){
    if(sum(stbl_1$MODULE==i)>0 & sum(stbl_2$MODULE==i)>0){
      represent = append(represent, paste(name1, name2, sep = ", "))
      # mapping1 = append(mapping1, stbl_1$Significant.Mapping[stbl_1$MODULE==i])
      # mapping2 = append(mapping2, stbl_2$Significant.Mapping[stbl_2$MODULE==i])
    }
    else if(sum(stbl_1$MODULE==i)>0){
      represent = append(represent, name1)
      # mapping1 = append(mapping1, stbl_1$Significant.Mapping[stbl_1$MODULE==i])
      # mapping2 = append(mapping2, "")
    }
    else if(sum(stbl_2$MODULE==i)>0){
      represent = append(represent, name2)
      # mapping1 = append(mapping1, "")
      # mapping2 = append(mapping2, stbl_2$Significant.Mapping[stbl_2$MODULE==i])
    }
    else{
      represent = append(represent, "")
      # mapping1 = append(mapping1, "")
      # mapping2 = append(mapping2, "")
    }
  }
  
  supersets$Representation = represent
  column1 = paste0(name1,"_Mapping")
  column2 = paste0(name2,"_Mapping")
  # supersets[,column1]= mapping1
  # supersets[,column2] = mapping2
  
  # study1 = supersets[supersets$Representation==name1,]
  # study2 = supersets[supersets$Representation==name2,]
  # both_studies = supersets[supersets$Representation==paste(name1,name2,sep = ", "),]

  # supersets_reorder = rbind(both_studies, study1)
  # supersets_reorder = rbind(supersets_reorder, study2)
  
  # supersets_reorder = supersets_reorder[,c(1,6,2,3,4,5,7,8)]
  supersets = supersets[,c(1,2,9,3,4,5,6,7,8)]
  
  write.table(supersets, paste0(output_Dir, "detailed_superset_info_rep.txt"),
              row.names = FALSE, quote = FALSE, sep="\t")
  # write.table(supersets_reorder, paste0(output_Dir, "superset_info_rep_reordered.txt"),
  #            row.names = FALSE, quote = FALSE, sep="\t")
  
}


# get numbers 
# 2 studies
# results_Dir is preprocess_quasi_overlap
quantifyOverlap <- function(results_Dir, studies, quasi_Dir){
  combined_files = list.files(results_Dir)[grep("combined_", list.files(results_Dir))]
  study1_num = c()
  study2_num = c()
  
  shared = c() # number
  shared_modules = c() # actual modules - at the end, unique to get number of unique shared modules
  
  all_study1_modules = c() # unique at end - overall across all mappings
  all_study2_modules = c()
  # at end after "uniquing", get num of shared modules, number of duplicated
  
  for(file in combined_files){
    mapping = read.delim(paste0(results_Dir,file), stringsAsFactors = FALSE)
    study1_num = append(study1_num, length(mapping$STUDY[mapping$STUDY==studies[1]]))
    study2_num = append(study2_num, length(mapping$STUDY[mapping$STUDY==studies[2]]))
    shared = append(shared, sum(duplicated(mapping$MODULE)))
    
    shared_modules = append(shared_modules,mapping$MODULE[duplicated(mapping$MODULE)])
    all_study1_modules = append(all_study1_modules, mapping$MODULE[mapping$STUDY==studies[1]])
    all_study2_modules = append(all_study2_modules, mapping$MODULE[mapping$STUDY==studies[2]])
  }
  
  quasi_overlapped = c()
  mappings = gsub("combined_","",gsub(".txt","",combined_files))
  shared_files = list.files(quasi_Dir)[grep("shared_", list.files(quasi_Dir))]
  for(map in mappings){
    if(sum(grepl(map,shared_files))==0)
      quasi_overlapped = append(quasi_overlapped, 0)
    else{
      file = shared_files[grep(map, shared_files)]
      shared_file = read.delim(paste0(quasi_Dir, file), stringsAsFactors = FALSE)
      quasi_overlapped = append(quasi_overlapped, length(shared_file$MODULE))
    }
  }
  
  superset_overlap = c()
  merged_files = list.files(quasi_Dir)[grep("merged_", list.files(quasi_Dir))]
  merged_info = merged_files[grep("info.txt", merged_files)]
  for(map in mappings){
    if(sum(grepl(map, merged_info))==0)
      superset_overlap = c(superset_overlap, 0)
    else{
      file = merged_info[grep(map, merged_info)]
      merged_file = read.delim(paste0(quasi_Dir,file), stringsAsFactors = FALSE)
      superset_overlap = append(superset_overlap, length(merged_file$MODULE))
    }
  }
  
  quantified = data.frame("Mapping"= mappings,stringsAsFactors = FALSE)
  
  quantified[,studies[1]] = study1_num
  quantified[,studies[2]] = study2_num
  quantified$`Number of Shared Modules` = shared
  quantified$`Number of Shared Modules with Quasi Overlap` = quasi_overlapped
  quantified$`Number of Shared Supersets` = superset_overlap
  
  # address no merged modules because there was only 1 module
  quantified$`Number of Shared Supersets`[quantified$`Number of Shared Modules with Quasi Overlap`==1] <- 1 
  
  return(quantified)
}

# study1num = length(unique(all_study1_modules))
# study2num = length(unique(all_study2_modules))
# all = append(unique(all_study1_modules), unique(all_study2_modules))
# shared_overall = sum(duplicated(all))
# shared_byMapping = length(unique(shared_modules))

# take shared modules and get genes from each study. only direct. find consistent modules and unique ones for each different mappings. 
# method to get SGM is case specific. second to last thing separated by "." in the file
# value_cutoff is association value - only show genes from SNPs that have values greater than the cutoff (lower p-value)
# modify function such that it looks in the kbr.mod.txt file and checks if the shared genes are part of the pathways because in the details file,
# it will output all the genes that are mapped to by that SNP but are not necessarily part of the pathway (but at least one will be)
# ^^ changed this so it doesnt do that in this output! 
# ^^ make sure genes are actually in module - this is what parameter modfile is for
# there is probably a better, modular way I could do this that can do more than 2 studies.
createGenesFiles <- function(results_Dir, genes_Dir, output_Dir, studies, protein_descrip, value_cutoff, modfile){
  modfile <- read.delim(modfile, stringsAsFactors = FALSE)
  shared_files = list.files(results_Dir)[grep("shared_", list.files(results_Dir))]
  all = data.frame(stringsAsFactors = FALSE)
  gene_files = c()
  for(s in studies){
    gene_files = append(gene_files, list.files(genes_Dir)[grep(s, list.files(genes_Dir))])
  }
  for(file in shared_files){
    # direct overlap
    mapping = read.delim(paste0(results_Dir,file), stringsAsFactors = FALSE)
    SGM = unlist(strsplit(file, split = ".", fixed = TRUE))[2]
    cat(SGM,"\n")
    sig_modules = mapping$MODULE
    
    genes_s1 = c()
    genes_s2 = c()
    descrip_s1 = c()
    descrip_s2 = c()
    
    shared_genes = c()
    shared_descrip = c()
    
    all_long = data.frame(stringsAsFactors = FALSE)
    for(mod in sig_modules){
      cat(mod,"\n")
      
      #gene_list1 = c() # to get shared genes
      #gene_list2 =c() # to get shared genes
      
      details = gene_files[grep(SGM, gene_files)]
      
      details_s1 = read.delim(paste0(genes_Dir,details[grep(studies[1], details)]), stringsAsFactors = FALSE)
      # this function will modify genes_s1 and descrip_s1 and return the genes
      gene_set1 <- retrieveGenes(genes_df = details_s1, genes = genes_s1, descrip = descrip_s1, mod = mod, modfile = modfile, value_cutoff = value_cutoff)
      genes_s1 = append(genes_s1, gene_set1[["Concatenated genes"]])
      descrip_s1 = append(descrip_s1, gene_set1[["Concatenated descrip"]])
      # get genes in list 
      # if(sum(details_s1$MODULE==mod)==0){
      #   genes_s1 = append(genes_s1, "none")
      #   descrip_s1 = append(descrip_s1, "none")
      # }
      # else{
      #   if(length(details_s1$GENE[details_s1$MODULE==mod & details_s1$VALUE>value_cutoff])==0){
      #     genes_s1 = append(genes_s1, "Genes not meeting association value threshold: ")
      #     descrip_s1 = append(descrip_s1, "Genes not meeting association value threshold: ")
      #   }
      #   for(cell in details_s1$GENE[details_s1$MODULE==mod & details_s1$VALUE>value_cutoff]){
      #     if(grepl(",",cell)){
      #       gene_list1 = append(gene_list1, unlist(strsplit(cell,split = ",")))
      #     }
      #     else{
      #       gene_list1 = append(gene_list1, cell)
      #     }
      #   }
      #   genes_s1 = append(genes_s1, concatenate(gene_list1, mysep = ", "))
      #   descrip_s1 = append(descrip_s1, concatenate(annotate_gene(gene_list1), mysep = ", "))
      # }  
      details_s2 = read.delim(paste0(genes_Dir,details[grep(studies[2], details)]), stringsAsFactors = FALSE)
      gene_set2 <- retrieveGenes(genes_df = details_s2, genes = genes_s2, descrip = descrip_s2, mod = mod, modfile = modfile, value_cutoff = value_cutoff)
      genes_s2 = append(genes_s2, gene_set2[["Concatenated genes"]])
      descrip_s2 = append(descrip_s2, gene_set2[["Concatenated descrip"]])
      # if(sum(details_s2$MODULE==mod)==0){
      #   genes_s2 = append(genes_s2, "none")
      #   descrip_s2 = append(descrip_s2, "none")
      # }
      # else{
      #   for(cell in details_s2$GENE[details_s2$MODULE==mod & details_s1$VALUE>value_cutoff]){
      #     if(grepl(",",cell)){
      #       gene_list2 = append(gene_list2, unlist(strsplit(cell,split = ",")))
      #     }
      #     else{
      #       gene_list2 = append(gene_list2, cell)
      #     }
      #   }
      #   genes_s2 = append(genes_s2, concatenate(gene_list2, mysep = ", "))
      #   descrip_s2 = append(descrip_s2, concatenate(annotate_gene(gene_list2), mysep = ", "))
      # }
      
      shared_genes = append(shared_genes, concatenate(intersect(gene_set1[["Gene list"]], gene_set2[["Gene list"]]), mysep = ", "))
      shared_descrip = append(shared_descrip, concatenate(annotate_gene(intersect(gene_set1[["Gene list"]], gene_set2[["Gene list"]]))))
      shared = intersect(gene_set1[["Gene list"]], gene_set2[["Gene list"]])
      # for specifically gene2loci T2D, it is duplicating specificaly LDLC_Control 4 times??? (5 sets of the same shared_genes) 
      # tried to debug, could not figure out why
      if(length(intersect(gene_set1[["Gene list"]], gene_set2[["Gene list"]]))>0){
        genes_df_long = data.frame("Mapping"=SGM,
                                   "MODULE" = mod,
                                   "Shared genes" = shared,
                                   "Gene Descriptions" = annotate_gene(shared), stringsAsFactors = FALSE)
        cat(" ", length(shared), "\n")
      }
      all_long = rbind(all_long, genes_df_long[!duplicated(genes_df_long),]) #???? how are things getting duplicated?
      all_long = all_long[!duplicated(all_long),]
    }
    #ugh = all_long[duplicated(all_long),]
    write.table(all_long, paste0(output_Dir,"long_",SGM,"_overlap_genes.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
    genes_df = data.frame("Mapping" = SGM, 
                          "MODULE" = mapping$MODULE,
                          "DESCR"=mapping$DESCR,
                          stringsAsFactors = FALSE)
    
    genes_df$`Overlap Method` = mapping$METHOD
    genes_df[,paste0(studies[1],"_genes")] = genes_s1
    genes_df[,paste0(studies[2],"_genes")] = genes_s2
    genes_df[,paste0(studies[1],"_genes_descrip")] = descrip_s1
    genes_df[,paste0(studies[2],"_genes_descrip")] = descrip_s2
    genes_df$`Shared genes` = shared_genes
    genes_df$`Shared genes description` = shared_descrip
    #**problem: overflow of cell so text file can't handle it, results in losing some columns (?)
    # do a version with no description
    
    write.table(genes_df[,c("Mapping","MODULE","DESCR","Overlap Method",paste0(studies[1],"_genes"), paste0(studies[2],"_genes"), "Shared genes")], 
                paste0(output_Dir,"noDescrip_",SGM,"_overlap_genes.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
    #write.table(genes_df, paste0(output_Dir,SGM,"_overlap_genes.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
    #write.csv(genes_df, paste0(output_Dir,SGM,"_overlap_genes.txt"), row.names = FALSE, quote = FALSE, )
    
    all = rbind(all, genes_df)
    
  }
  write.table(all, paste0(output_Dir,"All_mapping_Genes.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
}

# this function modifies genes_s1!
retrieveGenes <- function(genes_df, genes, descrip, mod, modfile, value_cutoff){
  gene_list1 = c() # will modularize this later....
  set = list()
  if(sum(genes_df$MODULE==mod)==0){
    genes = "none"
    descrip = "none"
    set[["Gene list"]] = ""
  }
  else{
    if(length(genes_df$GENE[genes_df$MODULE==mod & genes_df$VALUE>value_cutoff])==0){
      genes = "Genes not meeting association value threshold: "
      descrip = "Genes not meeting association value threshold: "
      set[["Gene list"]] = ""
    }
    else{
      for(cell in genes_df$GENE[genes_df$MODULE==mod & genes_df$VALUE>value_cutoff]){
        if(grepl(",",cell)){
          gene_list1 = append(gene_list1, unlist(strsplit(cell,split = ",")))
        }
        else{
          gene_list1 = append(gene_list1, cell)
        }
        keep_indices = which(!is.na(match(gene_list1, modfile$GENE[modfile$MODULE==mod])))
        gene_list1 = gene_list1[keep_indices]
      }
      genes = concatenate(gene_list1, mysep = ", ")
      descrip = concatenate(annotate_gene(gene_list1), mysep = ", ")
      set[["Gene list"]] = gene_list1
    }
  }
  set[["Concatenated genes"]] = genes
  set[["Concatenated descrip"]] = descrip
  return(set)
}

# must have protein_names data frame 
# changed vector$GENE to vector
annotate_gene <- function(vector){
  info = c()
  for(m in vector){
    if(length(protein_names$annotation[protein_names$preferred_name==m])==0){
      if(length(protein_names$annotation[protein_names$preferred_name==gsub("[[:digit:]]","",m)])==0){
        info = append(info, "Not annotated")
      }
      else{
        info = append(info, paste0("Possible annotation: ",
                                   concatenate(protein_names$annotation[protein_names$preferred_name==gsub("[[:digit:]]","",m)],
                                               mysep = ",")))
      }
    }
    else if(length(protein_names$annotation[protein_names$preferred_name==m])>1){
      info = append(info, 
                    concatenate(protein_names$annotation[protein_names$preferred_name==m],
                                mysep = ","))
      cat("This protein had more than one annotation:", m, "\n")
    }
    else{
      info = append(info, protein_names$annotation[protein_names$preferred_name==m])
    }
  }
  return(info)
}


replicate <- function(results_Dir, FDR_trim, FDR_cutoff, info_file){
  files = list.files(results_Dir)
  mapping = c()
  traits = c()
  for(f in files){
    mapping = append(mapping, unlist(strsplit(f, split = ".",fixed = TRUE))[2]) 
    traits = append(traits, unlist(strsplit(f, split = ".",fixed = TRUE))[1]) 
  }
  mapping = unique(mapping)
  traits = unique(traits)
  total = data.frame(stringsAsFactors = FALSE)
  for(map in mapping){
    mapping_files = files[grep(map,files)]
    mapping_total = data.frame()
    for(file in mapping_files){
      df = read.delim(paste0(results_Dir,file), stringsAsFactors = FALSE)
      df = df[df$FDR<FDR_trim,]
      df = df[!grepl("_ctrlA",df$MODULE),]
      df = df[!grepl("_ctrlB",df$MODULE),]
      if(length(df$MODULE)==0 | is.na(df$FDR[1])){
        next
      }
      df$Trait = unlist(strsplit(file, split = ".", fixed = TRUE))[1]
      df$Mapping = unlist(strsplit(file, split = ".", fixed = TRUE))[2]
      mapping_total = rbind(mapping_total, df)
    }
    if(dim(mapping_total)[1]==0){
      next
    }
    mapping_ls = list()
    for(trait in 1:length(traits)){
      mapping_ls[[traits[trait]]] = c()
      mapping_ls[[paste0(traits[trait],"_FDR")]] = c()
      for(module in unique(mapping_total$MODULE)){
        if(sum(mapping_total$Trait[mapping_total$MODULE==module]==traits[trait])==0){
          mapping_ls[[traits[trait]]] = append(mapping_ls[[traits[trait]]], "NO")
        }
        else if(mapping_total$FDR[mapping_total$MODULE==module & mapping_total$Trait==traits[trait]]>FDR_cutoff){
          mapping_ls[[traits[trait]]] = append(mapping_ls[[traits[trait]]], "NO")
        }
        else{
          mapping_ls[[traits[trait]]] = append(mapping_ls[[traits[trait]]], "YES")
        }
        if(sum(mapping_total$Trait[mapping_total$MODULE==module]==traits[trait])>0){
          mapping_ls[[paste0(traits[trait],"_FDR")]] = append(mapping_ls[[paste0(traits[trait],"_FDR")]], 
                                                           mapping_total$FDR[mapping_total$MODULE==module & mapping_total$Trait==traits[trait]])
        }
        else{
          mapping_ls[[paste0(traits[trait],"_FDR")]] = append(mapping_ls[[paste0(traits[trait],"_FDR")]], "NS")
        }
      }
    }
    mapping_df = data.frame("Mapping" = map,
                            "MODULE"=unique(mapping_total$MODULE))
    for(item in 1:length(mapping_ls)){
      #mapping_df = cbind(mapping_df, mapping_ls[[item]])
      mapping_df[,names(mapping_ls)[item]] = mapping_ls[[item]]
    }
    total = rbind(total, mapping_df)
  }
  # sum "YES"s in a row
  rownames(total) <- seq(length=nrow(total))
  YESs = c()
  for(i in 1:nrow(total)){
    YESs[i] = sum(total[i,]=="YES")
  }
  total["n_Overlap"] = YESs
  total = total[!(total$n_Overlap<2),] # get rid of no overlaps
  total = total[order(total$n_Overlap, decreasing = TRUE),]
  total = addDESCR(df = total, position_to_add = 3, info_file)
  return(total)
}

# add DESCR column to df at a position position_to_add
# assumes df has MODULE column
addDESCR <- function(df, position_to_add, info_file){
  info <- read.delim(info_file, stringsAsFactors = FALSE)
  descr = c()
  for(i in df$MODULE){
    descr = append(descr, info$DESCR[info$MODULE==i])
  }
  df$DESCR = descr
  col_names = colnames(df)
  post_total = length(colnames(df))+1
  before = 1:(position_to_add-1)
  after = (position_to_add):(post_total-1)
  # indices = c(before, post_total, after)
  df = df[,c(col_names[before], "DESCR",col_names[after])]
  return(df)
}

subset_replicate <- function(df, studies){
  result = data.frame(stringsAsFactors = FALSE)
  for(i in studies){
    temp = df[df[,i]=="YES",]
    result = rbind(result, temp)
  }
}

makeDf <- function(Dir){
  files = list.files(Dir)[grep("shared_", list.files(Dir))]
  all = data.frame(stringsAsFactors = FALSE)
  for(foo in files){
    SGM = unlist(strsplit(foo, split = ".", fixed = TRUE))[2]
    mods = read.delim(paste0(Dir,foo), stringsAsFactors = FALSE)
    mods$Mapping = SGM
    mods = mods[,c(4,1,3,2)]
    all = rbind(all, mods)
  }
  return(all)
}

describePathways <- function(modules_list, df, DESCR){
  mappings = c()
  num = c()
  for(mod in modules_list){
    num = append(num, length(df$Mapping[df$MODULE==mod]))
    mappings = append(mappings, concatenate(df$Mapping[df$MODULE==mod], mysep = ", "))
  }
  describe = c()
  for(mod in modules_list){
    describe = append(describe, all_DESCR$DESCR[all_DESCR$MODULE==mod])
    if(length(all_DESCR$DESCR[all_DESCR$MODULE==mod])==0)
      cat(mod)
  }
  new_df = data.frame("MODULE"=modules_list, 
                      "DESCR"=describe,
                      "Mapping" = mappings, 
                      "nMapping" = num, stringsAsFactors = FALSE)
  new_df = new_df[order(new_df$nMapping, decreasing = TRUE),]
  return(new_df)
}


# find consistent modules and unique ones for each different mappings. Hm, I already did this.

# Purpose: 
# 1) consolidate results from all mappings 
# 2) take genes result file and append to each module in a long format. (each gene separated by comma)
# 3) add gene information 



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
#5. annotatkee supersets resulting from step 4 using annotate_supersets function
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


#took out Esophagus_Muscularis.GTEXv7.MEGENA_168,..
