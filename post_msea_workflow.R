# Workflow from getting MSEA results (multiple mapping methods) --------------
# file format: trait.mapping_method.LDfiltering.topPercentAssociations.results.txt !!! Important for all functions to work!!!!!
#             example: DIAGRAMstage1_T2D.Artery_Aorta_sQTL.50.20.results.txt
#             explanation: DIAGRAMstage1_T2D is the trait (GWAS), Artery_Aorta_sQTL is a tissue-specific splicing QTL 
#                          50% LD filtering, top 20% associations (SNPs)  
# Notes for me; my study:
#       5 studies: DIAGRAMstage1_T2D, UKB_T2D, Wojcik_T2D, Transethnic_AD, UKB_AD 
#       For each study, ran eQTLs and sQTLs for 49 tissues - GTEx v8
#       Ran canonical (Kegg, Biocarta, Reactome) and coexpression (WGCNA, MEGENA)
#       Coexpression - 120 eQTL, 120 sQTL (480 files)
#       Canonical - 

source("/Users/jessicading/Desktop/Yang_Lab/source/analyze_msea.R")
source("/Users/jessicading/Desktop/Yang_Lab/Mergeomics.R") # has tool.read function

# This study combines any separate msea run (for the same mapping method)
# For my study, I want to combine coexpression and canonical together for each mapping method (ran separately for both)
# (for my study, 24 of the tissues have coexpression results)
# Folders must have identical file names (ex. doesn't distinguish b/w canonical and coexpression)
# Can have more than two folders
combine_results(folders = c("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/canonical_results/",
                            "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/coexpression_results/"), 
                output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/all/") # combined file will have the same name
combine_details(canon_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/canonical_genes/", 
                coexp_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/coexpression_genes/", 
                output_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/all_genes/")

# Examine study consistency -------
# (outputs table that shows for each (at least 2) overlapped pathway, which studis (traits) did the pathway show up)
# Required file name: items separated by "." - ex. DIAGRAMstage1_T2D.Artery_Aorta_sQTL.50.20.results.txt
# df <- replicate(results_Dir = "~/Downloads/temp4/", FDR_trim = .25, FDR_cutoff = 0.05, 
#                 info_file = "~/Desktop/Yang_Lab/Resources/genesets/all.info.txt")

All_study_rep <- replicate(results_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/all/", FDR_trim = .25, FDR_cutoff = 0.05,
                info_file = "~/Desktop/Yang_Lab/Resources/genesets/all.info.txt")
setwd("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/")
write.table(All_study_rep, "./Supplement/All_study_representation.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Subset table based on which traits had the pathway enriched ("YES")
# For example, want to see overlapped pathways in AD studies
AD_replication = All_study_rep[All_study_rep$Transethnic_AD=="YES" & All_study_rep$UKB_AD=="YES",]
# I have 3 T2D studies and I want show pathways that are replicated in at least two of the 3 GWAS
T2D_replication = All_study_rep[(All_study_rep$DIAGRAMstage1_T2D=="YES" & All_study_rep$UKB_T2D=="YES") |
                                  (All_study_rep$DIAGRAMstage1_T2D=="YES" & All_study_rep$Wojcik_T2D=="YES") |
                                  (All_study_rep$UKB_T2D=="YES" & All_study_rep$Wojcik_T2D=="YES"),]
#
AD_T2D_replication = All_study_rep[((All_study_rep$DIAGRAMstage1_T2D=="YES" & All_study_rep$UKB_T2D=="YES") |
                                     (All_study_rep$DIAGRAMstage1_T2D=="YES" & All_study_rep$Wojcik_T2D=="YES") |
                                     (All_study_rep$UKB_T2D=="YES" & All_study_rep$Wojcik_T2D=="YES")) & 
                                     (All_study_rep$Transethnic_AD=="YES" & All_study_rep$UKB_AD=="YES"),]

# Quasi Overlap -------
# next do quasi overlap as the results from this function for other functions
# haven't written a version of geneOverlap that finds the quasi overlap between 3 studies so doing it separately
results_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/all/"
output_Dirs = c("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/", "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/", 
                "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_Wojcik/", "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/UKB_Wojcik/")
rcutoff = 0.30
FDR = 0.05
cohorts = list("AD" = c("Transethnic_AD","UKB_AD"), "T2D_1" = c("DIAGRAMstage1_T2D", "UKB_T2D"), 
               "T2D_2" = c("DIAGRAMstage1_T2D", "Wojcik_T2D"), "T2D_3" = c("UKB_T2D","Wojcik_T2D"))
modfile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.mod.txt"
infofile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.info.txt"
studies = c("AD","DIAGRAM_UKB", "DIAGRAM_Wojcik", "UKB_Wojcik")
processed_Dir = "overlap/" # to append to output_Dirs ex. ~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/ (got to make that directory)
for(iter in 3){
  cat("Now interrogating ", studies[iter], "\n")
  geneOverlap(results_Dir = results_Dir,
              #output_Dir = output_Dirs[iter], 
              output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/UKB_Wojcik_new/",
              perc_gene_overlap = 0.30, 
              fdr_cutoff = 0.05, 
              cohorts = cohorts[[iter]], 
              modfile_path = modfile_path, 
              infofile_path = infofile_path, 
              study = studies[iter]) # just to name the files
  files = list.files(output_Dirs[iter]) # choose
  files = files[grep("shared_", files)]
  for(i in files){
    if(iter==3 & grepl("Artery_Tibial_eQTL", i)) next # did not have time to figure out why this wasn't working. error using tool.coalesce 
                                                      # Error:  Error in `$<-.data.frame`(`*tmp*`, "COUNT", value = 0) : replacement has 1 row, data has 0 
                                                      # will have to know later that Artery_Tibial had 2 modules overlapped (for quantifyOverlap)
    name = unlist(strsplit(i, split = ".",fixed = TRUE))[2]
    cat("Merging ", name, "\n")
    modules_df = read.delim(paste0(output_Dirs[iter], i), header = TRUE, stringsAsFactors = FALSE)
    if(length(modules_df$MODULE)==0) next
    if(length(modules_df$MODULE)==1){
      write.table(modules_df, paste0(paste0(output_Dirs[iter], processed_Dir), "merged_",name,".info.txt"), row.names = FALSE, quote = FALSE, sep="\t")
      cat(name, " had only 1 significant module. Making fake merged file.\n") # 10/12/19 reminder - will have to redo this with this part
      next
    }
    merge_modules(name = name, 
                  rcutoff = 0.30,
                  modules_df = modules_df, 
                  output_Dir = paste0(output_Dirs[iter], processed_Dir), 
                  modfile_path = modfile_path, 
                  infofile_path = infofile_path)
  }
}


# Next, do quasi overlap on quasi results from AD and T2D (meta-quasi? :P)
# First, have to combine T2D quasi results as there was three studies (at least studies - considered overlap)
# each folder has for every mapping shared_quasi results
folders <- c("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/UKB_Wojcik/overlap/",
             "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_Wojcik/overlap/",
             "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/overlap/")
# get mapping methods
mapping = c()
for(f in folders){
  files = list.files(f)
  files = files[grep("shared_", files)]
  for(i in files){
    mapping = append(mapping, unlist(strsplit(i, split = ".",fixed = TRUE))[2]) 
  }
}
mapping = unique(mapping)
output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D.2/"
for(m in mapping){
  combined_mapping = data.frame(stringsAsFactors = FALSE)
  for(f in folders){
    files = list.files(f)
    files = files[grep("shared_", files)]
    files = files[grep(m, files)]
    if(length(files)==0){
      next
    }
    if(length(files)>1){
      cat("Error: multiple mappings!\n")
      next
    }
    dtfr = read.delim(paste0(f,files), stringsAsFactors = FALSE)
    combined_mapping = rbind(combined_mapping, dtfr)
  }
  combined_mapping = combined_mapping[!duplicated(combined_mapping$MODULE),]
  write.table(combined_mapping,paste0(output_Dir,"shared_T2D_quasi.",m,".txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")
}
# copy shared_AD_quasi and shared_T2D_quasi results to results_Dir below
shared_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/overlap/"
shared_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D_new/"
geneOverlap(results_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/",
            output_Dir = shared_Dir, 
            perc_gene_overlap = 0.30, 
            fdr_cutoff = 0.05, 
            cohorts = c("shared_AD", "shared_T2D_quasi"), 
            modfile_path = modfile_path, 
            infofile_path = infofile_path,
            type = "meta",
            study = "AD_T2D_quasi")
files = list.files("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/") 
files = files[grep("shared_", files)]
for(i in files){
  name = unlist(strsplit(i, split = ".",fixed = TRUE))[2]
  cat("Merging ", name, "\n")
  modules_df = read.delim(paste0("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/", i), header = TRUE, stringsAsFactors = FALSE)
  if(length(modules_df$MODULE)==0) next
  if(length(modules_df$MODULE)==1){
    write.table(modules_df, paste0("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/", "merged_",name,".info.txt"), row.names = FALSE, quote = FALSE, sep="\t")
    cat(name, " had only 1 significant module. Making fake merged file.\n")
    next
  }
  if(grepl("Muscle_Skeletal_eQTL",i)){
    write.table(modules_df, paste0("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/", "merged_",name,".info.txt"), row.names = FALSE, quote = FALSE, sep="\t")
    cat(name, " had only 1 significant module. Making fake merged file.\n")
    next
  } 
  merge_modules(name = name, 
                rcutoff = 0.30,
                modules_df = modules_df, 
                output_Dir = paste0("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/"), 
                modfile_path = modfile_path, 
                infofile_path = infofile_path)
}


# For T2D, make merged modules for combined shared
shared_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D/"
files = list.files(shared_Dir) 
files = files[grep("shared_", files)]
for(i in files){
  name = unlist(strsplit(i, split = ".",fixed = TRUE))[2]
  cat("Merging ", name, "\n")
  modules_df = read.delim(paste0(shared_Dir, i), header = TRUE, stringsAsFactors = FALSE)
  if(length(modules_df$MODULE)==0 | length(modules_df$MODULE)==1) next
  if(grepl("Muscle_Skeletal_eQTL",i)) next
  merge_modules(name = name, 
                rcutoff = 0.30,
                modules_df = modules_df, 
                output_Dir = shared_Dir, 
                modfile_path = modfile_path, 
                infofile_path = infofile_path)
}


# Select some to run KDA 
Nerve_Tibial_AD = read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/merged_Nerve_Tibial_eQTL.mod.txt", stringsAsFactors = FALSE)
Nerve_Tibial_T2D = read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D/merged_Nerve_Tibial_eQTL.mod.txt", stringsAsFactors = FALSE)
Nerve_tib_total = rbind(Nerve_Tibial_AD, Nerve_Tibial_T2D)
Nerve_tib_total = Nerve_tib_total[!duplicated(Nerve_tib_total),]
write.table(Nerve_tib_total,"./networks/combined_AD_T2D_nerve_tib.txt", row.names = FALSE, quote = FALSE, sep = "\t")
runKDA(nodes = "./networks/combined_AD_T2D_nerve_tib.txt",
       network = "/Users/jessicading/Downloads/nerve_tibial_BN_GTEx_V7_keep_DAG.txt",
       trim = c("merged_",".mod.txt"))

# AD specific pathways, T2D specific pathways, shared pathways -------
# specific overall
# specific by eQTL

# extract GWAS genes


# get AD_T2D shared modules
sharedDir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/"
shared_AD_T2D_files = list.files(sharedDir, full.names = "TRUE")[grep("shared_", list.files(sharedDir))]
AD_T2D = c()
for(foo in shared_AD_T2D_files){
  mods = read.delim(foo, stringsAsFactors = FALSE)
  AD_T2D = append(AD_T2D, mods$MODULE)
}
AD_T2D = unique(AD_T2D)

AD_files = list.files("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/", full.names = "TRUE")[grep("shared_", list.files("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/"))]
AD = c()
for(foo in AD_files){
  mods = read.delim(foo, stringsAsFactors = FALSE)
  AD = append(AD, mods$MODULE)
}
AD = unique(AD)

T2D_files = list.files("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D/", full.names = "TRUE")[grep("shared_", list.files("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D/"))]
T2D = c()
for(foo in T2D_files){
  mods = read.delim(foo, stringsAsFactors = FALSE)
  T2D = append(T2D, mods$MODULE)
}
T2D = unique(T2D)

AD_df = makeDf(Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/") # takes shared files from the quasi overlap
T2D_df = makeDf(Dir="~/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D/")

AD_T2D_df = makeDf(Dir="~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/")
AD_T2D_mod_df = describePathways(modules_list = AD_T2D, df = AD_T2D_df, DESCR = all_DESCR)
write.table(AD_T2D_mod_df,"./Supplement/Top_AD_T2D_pathways.txt",row.names = FALSE, quote = FALSE, sep="\t")

all_DESCR <- read.delim(infofile_path, stringsAsFactors = FALSE)
# make mod file that is inclusive of both disease-specific and shared b/w AD and T2D
AD_all_df = describePathways(modules_list = AD, df=AD_df, DESCR = all_DESCR)
T2D_all_df = describePathways(modules_list = T2D, df=T2D_df, DESCR = all_DESCR)

# from AD and T2D take out pathways that are shared to get disease-specific pathways
AD = setdiff(AD,AD_T2D) # changing AD and T2D for input to describePathways to get AD specific and T2D specific
T2D = setdiff(T2D, AD_T2D)
# numbers don't check out? probably because of quasi overlap... 

# make list of AD-specific pathways and mappings, same for T2D - eQTL specific!
AD_specific_df = describePathways(modules_list = AD, df = AD_df, DESCR = all_DESCR)
T2D_specific_df = describePathways(modules_list = T2D,df = T2D_df, DESCR = all_DESCR)

# take shared pathways but in different mapping methods
AD_T2D_shared_pathways_dif_SGMs = intersect(AD_specific_df$MODULE, T2D_specific_df$MODULE)
AD_mapping = c()
T2D_mapping = c()
for(m in AD_T2D_shared_pathways_dif_SGMs){
  AD_mapping = append(AD_mapping, AD_specific_df$Mapping[AD_specific_df$MODULE==m])
  T2D_mapping = append(T2D_mapping, T2D_specific_df$Mapping[T2D_specific_df$MODULE==m])
}
AD_T2D_shared_pathways_dif_SGMs_df = data.frame("MODULE"=AD_T2D_shared_pathways_dif_SGMs,
                                                "AD_mapping" = AD_mapping, 
                                                "T2D_mapping"=T2D_mapping, stringsAsFactors = FALSE)

AD_very_specific = AD_specific_df[setdiff(1:length(AD_specific_df$MODULE), 
                                          match(AD_T2D_shared_pathways_dif_SGMs, AD_specific_df$MODULE)),]
T2D_very_specific = T2D_specific_df[setdiff(1:length(T2D_specific_df$MODULE), 
                                           match(AD_T2D_shared_pathways_dif_SGMs, T2D_specific_df$MODULE)),]


setwd("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/")
write.table(AD_specific_df, "./Supplement/AD_specific_pathways.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(T2D_specific_df, "./Supplement/T2D_specific_pathways.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(AD_all_df, "./Supplement/AD_all_pathways.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(T2D_all_df, "./Supplement/T2D_all_pathways.txt", row.names = FALSE, quote = FALSE, sep="\t")

# will have to process these modules to create supersets (redundancy of pathways)

# OR write to excel file
library(WriteXLS)
MSEA_results_replication_results = list()
MSEA_results_replication_results[["AD_specific"]] = AD_specific_df
MSEA_results_replication_results[["T2D_specific"]] = T2D_specific_df
MSEA_results_replication_results[["AD_all"]] = AD_all_df
MSEA_results_replication_results[["T2D_all"]] = T2D_all_df
WriteXLS(MSEA_results_replication_results, "top5_all_direction_pathways.xls",names(MSEA_results_replication_results))


# KDA ----
# first prepare merged files (this is to get rid of redundancy)
# merge shared files 
# move shared_ files generated from the geneOverlap function into a new directory (should change script so the user doesn't have to do all these things..)
# first do for AD_T2D
files = list.files("~/Desktop/Yang_Lab/T2D_AD/Data/AD/05192019/quasi_overlap/") # choose
files = files[grep("shared_", files)]
for(i in files){
  name = gsub("shared_AD_quasi.", "", i)
  name = gsub(".txt", "", name)
  modules_df = read.delim(paste0("~/Desktop/Yang_Lab/T2D_AD/Data/AD/05192019/quasi_overlap/", i), header = TRUE, stringsAsFactors = FALSE)
  if(length(modules_df$MODULE)==0 | length(modules_df$MODULE)==1) next
  merge_modules(name = name, 
                rcutoff = 0.30,
                modules_df = modules_df, 
                output_Dir = "~/Desktop/Yang_Lab/T2D_AD/Data/AD/05192019/quasi_overlap/", 
                modfile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.mod.txt", 
                infofile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.info.txt")
}

# want to try method to combine brain_esnps between AD and T2D (not a shared pathway but could be connected in a network?)
# shared b/w AD T2D, then AD-specific, then T2D-specific
brain_AD = read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/merged_brain_esnps_eQTL.mod.txt", stringsAsFactors = FALSE)
brain_T2D = read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D/merged_brain_esnps_eQTL.mod.txt", stringsAsFactors = FALSE)
# shared modules - HDLC_Control and CHB_positive
# take out modules from AD and T2D (just for sake of not duplicate edges)
shared_mods = brain_AD[brain_AD$MODULE=="HDLC_Control" | brain_AD$MODULE=="CHD_positive",]
shared_mods$MODULE <- paste(shared_mods$MODULE, "shared", sep = "_")
brain_AD = brain_AD[brain_AD$MODULE!="HDLC_Control",]
brain_AD = brain_AD[brain_AD$MODULE!="CHD_positive",]
brain_T2D = brain_T2D[brain_T2D$MODULE!="HDLC_Control",]
brain_T2D = brain_T2D[brain_T2D$MODULE!="CHD_positive",]
brain_AD$MODULE <- paste(brain_AD$MODULE, "AD", sep = "_")
brain_T2D$MODULE <- paste(brain_T2D$MODULE,"T2D",sep = "_")
total_brain = rbind(brain_AD, brain_T2D)
total_brain = rbind(total_brain, shared_mods)
# get other brain modules 
AD_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/"
T2D_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D/"
brain_cer_ctx_AD = list.files(AD_folder)[intersect(grep("merged_Brain", list.files(AD_folder)), 
                                                   grep(".mod.txt",list.files(AD_folder)))]
AD_mods = data.frame()
for(file in brain_cer_ctx_AD){
  temp = read.delim(paste0(AD_folder, file), stringsAsFactors = FALSE)
  AD_mods = rbind(AD_mods, temp)
}
AD_mods = AD_mods[!duplicated(AD_mods),]


brain_cer_ctx_T2D = list.files(T2D_folder)[intersect(grep("merged_Brain", list.files(T2D_folder)), 
                                                     grep(".mod.txt",list.files(T2D_folder)))]
T2D_mods = data.frame()
for(file in brain_cer_ctx_T2D){
  temp = read.delim(paste0(T2D_folder, file), stringsAsFactors = FALSE)
  T2D_mods = rbind(T2D_mods, temp)
}
T2D_mods = T2D_mods[!duplicated(T2D_mods),]

shared_T2D_AD = intersect(unique(AD_mods$MODULE), unique(T2D_mods$MODULE))
total = rbind(AD_mods, T2D_mods)
total = total[!duplicated(total),]
total_brain2 = rbind(total_brain, total)
  
write.table(total_brain, "./networks/Brain_AD_T2D_nodes3.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(total_brain2,"./networks/Brain_AD_T2D_nodes2.txt", row.names = FALSE, quote = FALSE, sep="\t")

runKDA(nodes = "./networks/Brain_AD_T2D_nodes.txt", 
       network = "~/Desktop/Yang_Lab/Resources/network/AD_T2D/pan_esnp/networks.hs.all.txt", 
       trim = c("merged_",".mod.txt"))

# add that the module was from AD or T2D

# Add GWAS hit information
total_brain3 <- trim_network(folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/Brain_AD_T2D_nodes.txt_networks.hs.all/cytoscape/", 
                             keep_only_orig_input = TRUE)
nodes <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/Brain_AD_T2D_nodes.txt_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt", 
                                   SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/sig_genes/", study = "AD")
nodes_T2D <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/Brain_AD_T2D_nodes.txt_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt", 
                                       SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/sig_genes/", study = "T2D")
nodes_T2D2 <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/Brain_AD_T2D_nodes.txt_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt", 
                                       SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_Wojcik/sig_genes/", study = "T2D")
# make my own brain network based on KDs I want to include?


# Refine Peripheral sQTL network --------
perip_edges <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/peripheral_ssnps_sQTL_networks.hs.all/cytoscape/kda2cytoscape.edges.txt", 
                          stringsAsFactors = FALSE)
perip_edges <- trim_network(folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/peripheral_ssnps_sQTL_networks.hs.all/cytoscape/", 
                            keep_only_orig_input = TRUE)
# get rid of double edges
perip_edges = perip_edges[!duplicated(cbind(pmin(perip_edges$HEAD, perip_edges$TAIL), pmax(perip_edges$HEAD, perip_edges$TAIL))),]

# pastel colors
# 9DBAD5  darker blue
# FAF3DD beige yellow
# 8FC1A9 darker green
# E0BBE4 purple
# FFDFD3 orange
# FEC8D8 pink/red
# DFBF9F brown 
# CCCCCC gray

nodes <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/peripheral_ssnps_sQTL_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt", 
                    stringsAsFactors = FALSE)

nodes$URL <- gsub("406399","9DBAD5", nodes$URL) # rctm0218,..
nodes$URL <- gsub("409900","FAF3DD", nodes$URL) # rctm0089,..
nodes$URL <- gsub("409999","8FC1A9", nodes$URL) # rctm0415,..
nodes$URL <- gsub("514099","E0BBE4", nodes$URL) # M5669
nodes$URL <- gsub("874099","FFDFD3", nodes$URL) # rctm1335
nodes$URL <- gsub("992500","FEC8D8", nodes$URL) # rctm0863
nodes$URL <- gsub("992559","DFBF9F", nodes$URL) # M15902
nodes$URL <- gsub("996900","CCCCCC", nodes$URL) # rctm0598,..

nodes <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/peripheral_ssnps_sQTL_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt", 
                                   SGM = "peripheral_ssnps_sQTL", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/sig_genes/", study = "AD")
nodes2 <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/peripheral_ssnps_sQTL_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt", 
                                    SGM = "peripheral_ssnps_sQTL", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/sig_genes/", study = "T2D")
nodes$GWAS_T2D = nodes2$GWAS_hit


# cbind last columns
# not adding those with Wojcik because there was no genes meeting association value threshold
GWAS_hit_meta = c()
for(i in 1:length(nodes$NODE)){
  if(nodes$GWAS_hit[i]!="no" & nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "AD_T2D"
  }
  else if(nodes$GWAS_hit[i]!="no"){
    GWAS_hit_meta[i] = "AD"
  }
  else if(nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "T2D"
  }
  else{
    GWAS_hit_meta[i] = "no"
  }
}
nodes$GWAS_hit_meta = GWAS_hit_meta

write.table(nodes, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/peripheral_ssnps_sQTL_networks.hs.all/cytoscape/mod_nodes.txt", row.names = FALSE, quote = FALSE, sep="\t")

# change colors of nodes
nodes$URL <- gsub("514099","99CCFF", nodes$URL) 
nodes$URL <- gsub("992559","FFD7D7", nodes$URL) 
nodes$URL <- gsub("992500","FFCC99", nodes$URL) 
nodes$URL <- gsub("874099","FFFF99", nodes$URL) 
nodes$URL <- gsub("406399","FF9999", nodes$URL) 
# get module info
module_color_mapping <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/Brain_AD_T2D_nodes.txt_networks.hs.all/cytoscape/")

# read in nodes and add column saying disease of origin. edges will retain module information 
# add other brain modules from AD and T2D, all of the coexpression ones as well....

           # SECOND METHOD to see shared AD/T2D casual interactions
# combine merged modules from AD and T2D from the same tissue (just rbind)
# interrogate the network - see if in the network structure, there is some interaction
# to module name, append 

# Get GWAS genes --------
# AD
protein_names = read.delim("~/Desktop/10X_Analysis/9606.protein.info.v11.0.txt", 
                           header = TRUE, stringsAsFactors = FALSE, quote = "")
for(iter in 3){
createGenesFiles(results_Dir = paste0(output_Dirs[iter],"overlap/"), 
                 genes_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/all_genes/", 
                 #output_Dir = paste0(output_Dirs[iter], "sig_genes/"), # make the genes folder
                 output_Dir = paste0(output_Dirs[iter], "genes.p.01/"), # make the genes folder May 2020
                 studies = cohorts[[iter]],
                 value_cutoff = 1.3, 
                 modfile = "~/Desktop/Yang_Lab/Resources/genesets/all.mod.txt") # from above geneOverlap
} # got error for iter 3 but appears to be done though?


# once get genes for AD and T2D, for a single tissue (for now), 
# can do a version of this that only subsets the shared for which the METHOD was "direct" because probably will run into issues with quasi overlap
#   (one study might not have genes for that module at all)

# output 

# gene directories
AD_Dir = paste0(output_Dirs[1], "sig_genes/")
AD_Dir = paste0(output_Dirs[1], "very_sig_genes/")
T2D_Dirs = c(paste0(output_Dirs[2:4], "sig_genes/"))
T2D_Dirs = c(paste0(output_Dirs[2], "very_sig_genes/"))
AD_T2D_shared_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/overlap/"
AD_T2D_shared_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D_new/"

# get tissues that had shared modules from AD_T2D (output from geneOverlap)

files = list.files(AD_T2D_shared_Dir)[grep("shared_", list.files(AD_T2D_shared_Dir))]


total = data.frame(stringsAsFactors = FALSE)
for(f in files){
  mapping = gsub("shared_AD_T2D_quasi.","",gsub(".txt","", f))
  cat(mapping,"\n")
  shared_df = read.delim(paste0(AD_T2D_shared_Dir, f), stringsAsFactors = FALSE)
  AD_files = list.files(AD_Dir)[grep("noDescrip_", list.files(AD_Dir))]
  AD_file = AD_files[grep(mapping, AD_files)]
  AD_genes_df = read.delim(paste0(AD_Dir,AD_file), stringsAsFactors = FALSE)
  
  for(mod in shared_df$MODULE){
    Shared = "possible"
    AD_genes = c()
    T2D_genes = c()
    if(sum(AD_genes_df$MODULE==mod)==0){
      AD_specific = ""
      Shared = ""
    }
    else{
      if(AD_genes_df$Shared[AD_genes_df$MODULE==mod]==""){
        for(i in 5:6){
          if(AD_genes_df[,i][AD_genes_df$MODULE==mod]!="none"){
            AD_genes = AD_genes_df[,i][AD_genes_df$MODULE==mod]
            break
          }
        }
      }
      else{
        AD_genes = AD_genes_df$Shared.genes[AD_genes_df$MODULE==mod]
      }
      AD_genes = unlist(strsplit(AD_genes, split = ", "))
      T2D_genes_df = data.frame(stringsAsFactors = FALSE)
      for(dir in T2D_Dirs){
        T2D_files = list.files(dir)[grep("noDescrip_", list.files(dir))]
        if(sum(grepl(mapping, T2D_files))==0) next
        T2D_file = T2D_files[grep(mapping, T2D_files)]
        T2D_genes_df = read.delim(paste0(dir,T2D_file), stringsAsFactors = FALSE)
        if(sum(T2D_genes_df$MODULE==mod)>0){
          cat("Breaking out of loop\n")
          break
        }
      }
      if(sum(T2D_genes_df$MODULE==mod)==0){
        T2D_specific = ""
        Shared = ""
      }
      else{
        if(T2D_genes_df$Shared.genes[T2D_genes_df$MODULE==mod]=="" | is.na(T2D_genes_df$Shared.genes[T2D_genes_df$MODULE==mod])){
          for(i in 5:6){
            if(T2D_genes_df[,i][T2D_genes_df$MODULE==mod]!="none"){
              T2D_genes = T2D_genes_df[,i][T2D_genes_df$MODULE==mod]
              break
            }
          }
        }
        else{
            T2D_genes = T2D_genes_df$Shared.genes[T2D_genes_df$MODULE==mod]
        }
        T2D_genes = unlist(strsplit(T2D_genes, split = ", "))
      }
      if(Shared!=""){
        Shared = intersect(AD_genes, T2D_genes)
        AD_specific = setdiff(AD_genes, Shared)
        AD_specific = concatenate(AD_specific, mysep = ", ")
        T2D_specific = setdiff(T2D_genes, Shared)
        T2D_specific = concatenate(T2D_specific, mysep = ", ")
        Shared = concatenate(Shared, mysep = ", ")
      }
      if(Shared=="" & length(AD_genes)!=0){
        AD_specific = concatenate(AD_genes, mysep = ", ")
      }
      if(Shared=="" & length(T2D_genes)!=0){
        T2D_specific = concatenate(T2D_genes, mysep = ", ")
      }
    }
    
    temp = data.frame("Mapping"=mapping,
                      "MODULE"=mod,
                      "AD Specific Genes" = AD_specific,
                      "T2D Specific Genes"= T2D_specific,
                      "Shared Genes" = Shared)
    total = rbind(total, temp)
  }
}

total1 <- addDESCR(df = total, position_to_add = 3, info_file = infofile_path)
write.table(total1, "./Supplement/Shared_AD_T2D_GWAS_Genes2.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(total1, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Supplement/Shared_AD_T2D_GWAS_Genes3.txt", row.names = FALSE, quote = FALSE, sep="\t")


# Quantify overlaps -----
AD_overlap <- quantifyOverlap(results_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/", 
                studies = c("Transethnic_AD","UKB_AD"), 
                quasi_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/")
T2D_overlap <- quantifyOverlap(results_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/", 
                               studies = c("DIAGRAMstage1_T2D","UKB_T2D"), 
                               quasi_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/overlap/")
T2D_2_overlap <- quantifyOverlap(results_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_Wojcik/", 
                               studies = c("DIAGRAMstage1_T2D","Wojcik_T2D"), 
                               quasi_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_Wojcik/overlap/")
T2D_3_overlap <- quantifyOverlap(results_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/UKB_Wojcik/", 
                                 studies = c("UKB","Wojcik_T2D"), 
                                 quasi_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/UKB_Wojcik/overlap/")

# just want to quantify shared_T2D_combined 
T2D_4_overlap <- quantifyOverlap(results_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_Wojcik/", 
                                 studies = c("DIAGRAMstage1_T2D","Wojcik_T2D"), 
                                 quasi_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D.2/")

setwd("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/Supplement/")
write.table(AD_overlap, "AD_numbers.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(T2D_overlap, "DIAGRAM_UKB_numbers.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(T2D_4_overlap, "Wojcik_actual_shared_numbers.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(AD_overlap[grep("_sQTL", AD_overlap$Mapping),], "102119_AD_numbers_sQTL.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(T2D_overlap[grep("_sQTL", T2D_overlap$Mapping),], "DIAGRAM_UKB_numbers_sQTL.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(T2D_4_overlap[grep("_sQTL", T2D_4_overlap$Mapping),], "102119_Wojcik_actual_shared_numbers_sQTL.txt", row.names = FALSE, quote = FALSE, sep="\t")

pseudo <- quantifyOverlap(results_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/", 
                              studies = c("Transethnic_AD","UKB_AD"), 
                              quasi_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/")
write.table(pseudo[grep("_sQTL", pseudo$Mapping),], "101219_AD_T2D_numbers_sQTL.txt", row.names = FALSE, quote = FALSE, sep="\t")

files <- c("/Users/jessicading/Downloads/temp6/shared_DIAGRAM_UKB.Adipose_Subcutaneous_eQTL.txt",
           "/Users/jessicading/Downloads/temp6/shared_DIAGRAM_Wojcik.Adipose_Subcutaneous_eQTL.txt",
           "/Users/jessicading/Downloads/temp6/shared_T2D.Adipose_Subcutaneous_eQTL.txt",
           "/Users/jessicading/Downloads/temp6/shared_UKB_Wojcik.Adipose_Subcutaneous_eQTL.txt")
total = data.frame(stringsAsFactors = FALSE)
for(file in files){
  df <- read.delim(file, stringsAsFactors = FALSE)
  df$Study = "T2D"
  total = rbind(total, df)
}
total = total[,c(1,3,4)]
total = total[!duplicated(total),]


# Make stacked bar plots --------
# version with supersets
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
eQTL_prop = read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/Supplement/proportion_figure_sQTL.txt", stringsAsFactors = FALSE)
eQTL_prop_melt2 = melt(eQTL_prop, id.vars = c("Mapping"))

eQTL_prop = read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/Supplement/proportion_figure_sQTL.txt", stringsAsFactors = FALSE)

eQTL_prop$AD_percent = (eQTL_prop$AD/sum(eQTL_prop$AD))*100
eQTL_prop$T2D_percent = (eQTL_prop$T2D/sum(eQTL_prop$T2D))*100
eQTL_prop$Shared_percent = (eQTL_prop$Shared/sum(eQTL_prop$Shared))*100
eQTL_prop_melt1 = melt(eQTL_prop, id.vars = c("Mapping","AD","T2D", "Shared")) # get percent column

eQTL_prop_final = cbind(eQTL_prop_melt2, eQTL_prop_melt1$value)
colnames(eQTL_prop_final) = c("SNP Mapping", "Disease","Numbers","Percentage")


eQTL_prop_final$`SNP Mapping` <- factor(eQTL_prop_final$`SNP Mapping`, levels = rev(eQTL_prop$Mapping))
eQTL_prop_final$Disease <- factor(eQTL_prop_final$Disease, levels = rev(c("AD","T2D","Shared")))

df2 <- eQTL_prop_final %>%
  group_by(Disease) %>%
  arrange(Disease, desc(`SNP Mapping`)) %>%
  mutate(lab_ypos = cumsum(Percentage) - 0.5 * Percentage)
# dont know why the second day I ran this it was shifting T2D numbers to the right by 100
#df2$lab_ypos[df2$Disease=="T2D"] = (df2$lab_ypos[df2$Disease=="T2D"]-100)

df2$Numbers[df2$Numbers==0] <- ""

#df2$`SNP Mapping` = gsub("eQTL","",df2$`SNP Mapping`)

pal <- c(
  "Combined" = "lightblue3",
  "Adipose" = "palevioletred1",
  "Adrenal Gland" = "royalblue4",
  "Artery" = "red1",
  "Brain" = "sienna1",
  "Fibroblasts" = "pink",
  "Colon" = "brown1",
  "Esophagus"="goldenrod1",
  "Heart" = "mediumseagreen",
  "Liver" = "lightskyblue",
  "Lung" = "lightsteelblue",
  "Muscle Skeletal" = "forestgreen", 
  "Nerve Tibial" = "lightseagreen",
  "Pancreas" = "royalblue1",
  "Pituitary" = "slateblue",
  "Skin" = "steelblue3",
  "Spleen" = "mediumpurple1",
  "Thyroid" = "orange4",
  "Whole Blood" = "gray61",
  "Distance" = "black",
  "ENCODE" = "salmon"
)

g <- ggplot(df2, aes(x=Disease, y=Percentage)) + geom_col(aes(fill=`SNP Mapping`)) +
  geom_text(aes(y=lab_ypos, label=Numbers, group = `SNP Mapping`), color = "white") + coord_flip() + 
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "top",
        axis.text.y = element_text(margin = margin(r=5), size=14, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #legend.title = element_text(face = "bold", size=13)) +
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_fill_manual(
    values = pal,
    limits = names(pal))

g


# top pathways for each mapping


pdf("./stacked_bar_superset_sQTL.pdf", height=2.5, width=10)
print(g)
dev.off()

# refine network
network <- read.delim("~/Downloads/kda2cytoscape.edges.txt", stringsAsFactors = FALSE)
nodes <- read.delim("~/Downloads/kda2cytoscape.nodes.txt", stringsAsFactors = FALSE)
prefixes = c("^MMT", "^NM_","^XM_","ENSMUST","^ri", "_at")
for(pre in prefixes){
  network = network[!grepl(pre,network$HEAD),]
}
write.table(network, "./networks/peripheral_ssnps.txt", row.names = FALSE, quote = FALSE, sep = "\t")
# keep only KDs and nodes with URLs
# make list of KDs and URL nodes
keep_nodes = c()
for(i in nodes$NODE){
  if(nodes$SHAPE[nodes$NODE==i]=="Diamond" | nodes$URL[nodes$NODE==i]!=""){
    keep_nodes = append(keep_nodes, i)
  }
}
new_net = network

remove_nodes = setdiff(nodes$NODE, keep_nodes)
for(r in remove_nodes){
  new_net = new_net[!grepl(r, new_net$HEAD),]
}
total = data.frame(stringsAsFactors = FALSE)
for(r in keep_nodes){
  temp = new_net[grep(r, new_net$HEAD),]
  total = rbind(total, temp)
}
total = total[!duplicated(total),]


# change dark colors
nodes$URL <- gsub("514099","99CCFF", nodes$URL) 
nodes$URL <- gsub("992559","FFD7D7", nodes$URL) 
nodes$URL <- gsub("992500","FFCC99", nodes$URL) 
nodes$URL <- gsub("874099","FFFF99", nodes$URL) 
nodes$URL <- gsub("406399","FF9999", nodes$URL) 
total = total[,c(1,2,3)]
total = total[!duplicated(total),]
write.table(total, "./networks/peripheral_ssnps_mod2.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(nodes, "./networks/mod_nodes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

summarizeKDA(KDA_folder = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/", 
             protein_descrip = "/Users/jessicading/Desktop/Yang_Lab/Epicardial/String/9606.protein.info.v11.0.txt", 
             name = "Peripheral_ssnps_network")


# Detailed table for AD_T2D --------
summaryTable(results_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/",
             output_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/", 
             study = "AD_T2D_quasi",
             infofile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.info.txt")








files <- c("/Users/jessicading/Downloads/temp6/shared_DIAGRAM_UKB.Adipose_Subcutaneous_eQTL.txt",
           "/Users/jessicading/Downloads/temp6/shared_DIAGRAM_Wojcik.Adipose_Subcutaneous_eQTL.txt",
           "/Users/jessicading/Downloads/temp6/shared_T2D.Adipose_Subcutaneous_eQTL.txt",
           "/Users/jessicading/Downloads/temp6/shared_UKB_Wojcik.Adipose_Subcutaneous_eQTL.txt")
total = data.frame(stringsAsFactors = FALSE)
for(file in files){
  df <- read.delim(file, stringsAsFactors = FALSE)
  df$Study = "T2D"
  total = rbind(total, df)
}
total = total[,c(1,3,4)]
total = total[!duplicated(total),]

geneOverlap(results_Dir = "~/Downloads/temp5/",
            output_Dir = "~/Downloads/temp9/", 
            perc_gene_overlap = 0.30, 
            fdr_cutoff = 0.05, 
            cohorts = c("Transethnic_AD","UKB_AD"), 
            modfile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.mod.txt", 
            infofile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.info.txt", 
            study = "AD")
df <- read.delim("/Users/jessicading/Downloads/temp6/shared_AD.Adipose_Subcutaneous_eQTL.txt", stringsAsFactors = FALSE)
df$Study = "AD"
df = df[,c(1,3,4)]
total = rbind(total,df)
total2 = total[duplicated(total$MODULE),]
# This function combines coexpression and canonical results for one study and puts them in a new location

# merge modules



runKDA(nodes = "/Users/jessicading/Downloads/temp6/merged_modified2.AD_T2D_Adipose.mod.txt", 
       network = "/Users/jessicading/Desktop/Yang_Lab/Resources/network/AD_T2D/Adipose/networks.hs.adipose.txt", trim = c(".mod.txt", "merged_"))
nodes <- read.delim("/Users/jessicading/Downloads/temp6/merged_AD_T2D_Adipose.mod.txt", stringsAsFactors = FALSE)
nodes = nodes[nodes$MODULE!="rctm0089",]
write.table(nodes, "/Users/jessicading/Downloads/temp6/merged_modified2.AD_T2D_Adipose.mod.txt", row.names = FALSE, quote = FALSE, sep = "\t")

summarizeKDA(KDA_folder = "~/Downloads/temp6/KDA/", 
             protein_descrip = "/Users/jessicading/Desktop/Yang_Lab/Epicardial/String/9606.protein.info.v11.0.txt", 
             name = "Adipose_AD_T2D_KDA")

# get genes for those modules...
protein_names = read.delim("~/Desktop/Yang_Lab/Epicardial/String/9606.protein.info.v11.0.txt", 
                           header = TRUE, stringsAsFactors = FALSE, quote = "")
combineCanonCoexp(canon_Dir = "/Users/jessicading/Downloads/genes_Wojcik/" , 
                  coexp_Dir = "/Users/jessicading/Downloads/genes_Wojcik_coexp/", 
                  combined_Dir = "/Users/jessicading/Downloads/temp9/")
createGenesFiles(results_Dir = "~/Downloads/temp7/", 
                 genes_Dir = "/Users/jessicading/Downloads/temp9/", 
                 output_Dir = "/Users/jessicading/Downloads/temp8/", 
                 studies = c("DIAGRAMstage1_T2D", "Wojcik_T2D"))


# do for sQTLs 
files <- c("/Users/jessicading/Downloads/temp6/shared_DIAGRAM_Wojcik.Adipose_Subcutaneous_sQTL.txt", 
           "/Users/jessicading/Downloads/temp6/shared_T2D.Adipose_Subcutaneous_sQTL.txt")
total = data.frame(stringsAsFactors = FALSE)
for(file in files){
  df <- read.delim(file, stringsAsFactors = FALSE)
  df$Study = "T2D"
  total = rbind(total, df)
}
total = total[,c(1,3,4)]
total = total[!duplicated(total),]

df <- read.delim("/Users/jessicading/Downloads/temp6/shared_AD.Adipose_Subcutaneous_sQTL.txt", stringsAsFactors = FALSE)
df$Study="AD"
df = df[,c(1,3,4)]
total = rbind(total,df)
total3 = total[duplicated(total$MODULE),]



# combine peripheral and brain networks together - shared modules ------
# what is contributing from brain, what from periphery
# from brain, shared with periphery
#   Positive control gene set for HDLC, Positive control gene set for coronary heart disease

# unique to AD and T2D but might be shared

# shared AD_T2D periphery 
AD_T2D_periph <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/overlap/shared_AD_T2D_quasi.peripheral_ssnps_sQTL.txt", 
                            stringsAsFactors = FALSE)
AD_T2D_periph = AD_T2D_periph[AD_T2D_periph$MODULE!="quasi",]

# shared AD_T2D brain
AD_T2D_brain <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/overlap/shared_AD_T2D_quasi.brain_esnps_eQTL.txt",
                           stringsAsFactors = FALSE)

# AD peripheral
AD_periph <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/shared_AD.peripheral_ssnps_sQTL.txt",
                        stringsAsFactors = FALSE)

# T2D peripheral
T2D_periph <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D/shared_T2D_quasi.peripheral_ssnps_sQTL.txt",
                         stringsAsFactors = FALSE)


# AD brain
AD_brain <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/shared_AD.brain_esnps_eQTL.txt", 
                       stringsAsFactors = FALSE)

# T2D brain
T2D_brain <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D/shared_T2D_quasi.brain_esnps_eQTL.txt", 
                       stringsAsFactors = FALSE)

AD_specific_periph = setdiff(AD_periph$MODULE, AD_T2D_periph$MODULE)
T2D_specific_periph = setdiff(T2D_periph$MODULE, AD_T2D_periph$MODULE)

AD_specific_periph_df = AD_periph[match(AD_specific_periph, AD_periph$MODULE),]
T2D_specific_periph_df = T2D_periph[match(T2D_specific_periph, T2D_periph$MODULE),]

AD_specific_brain = setdiff(AD_brain$MODULE, AD_T2D_brain$MODULE)
T2D_specific_brain = setdiff(T2D_brain$MODULE, AD_T2D_brain$MODULE)

AD_specfic_brain_df = AD_brain[match(AD_specific_brain, AD_brain$MODULE),]
T2D_specific_brain_df = T2D_brain[match(T2D_specific_brain, T2D_brain$MODULE),]

# merge AD_specific and T2D_specific modules
merge_modules(name = "AD_specific_periphery", 
              rcutoff = 0.30,
              modules_df = AD_specific_periph_df, 
              output_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/peripheral_not_shared/", 
              modfile_path = modfile_path, 
              infofile_path = infofile_path)
merge_modules(name = "T2D_specific_periphery", 
              rcutoff = 0.30,
              modules_df = T2D_specific_periph_df, 
              output_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/peripheral_not_shared/", 
              modfile_path = modfile_path, 
              infofile_path = infofile_path)
merged_AD_periphery <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/peripheral_not_shared/merged_AD_specific_periphery.mod.txt",
                                  stringsAsFactors = FALSE)
merged_AD_periphery$MODULE <- paste0(merged_AD_periphery$MODULE, "_AD")
merged_T2D_periphery <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/peripheral_not_shared/merged_T2D_specific_periphery.mod.txt",
                                   stringsAsFactors = FALSE)
merged_T2D_periphery$MODULE <- paste0(merged_T2D_periphery$MODULE, "_T2D")
total_periphery = rbind(merged_AD_periphery, merged_T2D_periphery)
write.table(total_periphery,"./networks/peripheral_not_shared/total_periph_mod.txt", row.names = FALSE, quote = FALSE, sep = "\t")
setwd("./networks/peripheral_not_shared/")
runKDA(nodes = "total_periph_mod.txt",
       network = "/Users/jessicading/Desktop/Yang_Lab/Resources/network/AD_T2D/pan_esnp/networks.hs.all.txt",
       trim = c("merged_","mod.txt"))
periph_edges <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/peripheral_not_shared/total_periph__networks.hs.all/cytoscape/kda2cytoscape.edges.txt",
                           stringsAsFactors = FALSE)
perip_edges <- trim_network(folder = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/peripheral_not_shared/total_periph__networks.hs.all/cytoscape/", 
                            keep_only_orig_input = TRUE)
# get rid of double edges
perip_edges = perip_edges[!duplicated(cbind(pmin(perip_edges$HEAD, perip_edges$TAIL), pmax(perip_edges$HEAD, perip_edges$TAIL))),]




merge_modules(name = "AD_specific_periphery", 
              rcutoff = 0.30,
              modules_df = AD_specfic_brain_df, 
              output_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/peripheral_not_shared/", 
              modfile_path = modfile_path, 
              infofile_path = infofile_path)
merge_modules(name = "AD_specific_periphery", 
              rcutoff = 0.30,
              modules_df = T2D_specfic_brain_df, 
              output_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/peripheral_not_shared", 
              modfile_path = modfile_path, 
              infofile_path = infofile_path)


# look at string PPI
# changed to very sig genes - 285 nodes

string_brain <- trim_network(folder = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network/cytoscape/", 
                             keep_only_orig_input = TRUE)
nodes <- add_GWAS_hit_info_network(nodes = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network/cytoscape/kda2cytoscape.nodes.txt", 
                                   SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/very_sig_genes/", study = "AD")
nodes2 <- add_GWAS_hit_info_network(nodes = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network/cytoscape/kda2cytoscape.nodes.txt", 
                                    SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/very_sig_genes/", study = "T2D")
nodes$GWAS_T2D = nodes2$GWAS_hit

# cbind last columns
# not adding those with Wojcik because there was no genes meeting association value threshold
GWAS_hit_meta = c()
for(i in 1:length(nodes$NODE)){
  if(nodes$GWAS_hit[i]!="no" & nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "AD_T2D"
  }
  else if(nodes$GWAS_hit[i]!="no"){
    GWAS_hit_meta[i] = "AD"
  }
  else if(nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "T2D"
  }
  else{
    GWAS_hit_meta[i] = "no"
  }
}
nodes$GWAS_hit_meta = GWAS_hit_meta
write.table(nodes, "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network/cytoscape/mod_nodes_3.txt", row.names = FALSE, quote = FALSE, sep="\t")

string_brain <- trim_network(folder = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network/cytoscape/", 
                             keep_only_orig_input = TRUE, keep_only_GWAS_hits = TRUE)
# I want to trim the network even further there are too many nodes... (presentation purposes)
for(iter in 1:2){
  createGenesFiles(results_Dir = paste0(output_Dirs[iter],"brain_genes/"), 
                   genes_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/all_genes/", 
                   output_Dir = paste0(output_Dirs[iter], "very_sig_genes/"), # make the genes folder
                   studies = cohorts[[iter]],
                   value_cutoff = 2, 
                   modfile = "~/Desktop/Yang_Lab/Resources/genesets/all.mod.txt") # from above geneOverlap
}


string_brain_new <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network/cytoscape/mod_kda2cytoscape.edges.txt",
                                    stringsAsFactors = FALSE)
Module_Trait = c()
for(i in 1:length(string_brain_new$TAIL)){
  if(grepl("AD", string_brain_new$MODULE[i])) Module_Trait[i] = "AD"
  else if(grepl("T2D", string_brain_new$MODULE[i])) Module_Trait[i] = "T2D"
  else if(grepl("shared", string_brain_new$MODULE[i])) Module_Trait[i] = "Shared"
  else Module_Trait[i] = "Unknown"
}
string_brain_new$Module_Trait = Module_Trait
write.table(string_brain_new, "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network/cytoscape/mod_kda2cytoscape.edges.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")



# Trim peripheral network based on Xia's suggestions
peripheral_take_out <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_kda2cytoscape.edges.txt", 
                                  stringsAsFactors = FALSE)
take_out = c("PTPRC", "ITGAL","BIRC5","CDCA8","SPAG5","CCNA2","MKI67")

for(t in take_out){
  peripheral_take_out = peripheral_take_out[!grepl(t, peripheral_take_out$TAIL),]
}
write.table(peripheral_take_out, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_kda2cytoscape.edges_2.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# add in brain 

brain_periph <- trim_network(folder = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/", 
                             keep_only_orig_input = TRUE, keep_only_GWAS_hits = TRUE)
# nodes <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt", 
#                                    SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/sig_genes/", study = "AD")
nodes <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt", 
                                   SGM = "peripheral_ssnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/sig_genes/", study = "AD")
# nodes3 <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt", 
#                                     SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/sig_genes/", study = "T2D")
nodes2 <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt", 
                                    SGM = "peripheral_ssnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/sig_genes/", study = "T2D")
nodes$GWAS_T2D = nodes2$GWAS_hit


including_Brain = nodes

# might just manually color nodes for the brain GWAS hits

# cbind last columns
# not adding those with Wojcik because there was no genes meeting association value threshold
GWAS_hit_meta = c()
for(i in 1:length(nodes$NODE)){
  if(nodes$GWAS_hit[i]!="no" & nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "AD_T2D"
  }
  else if(nodes$GWAS_hit[i]!="no"){
    GWAS_hit_meta[i] = "AD"
  }
  else if(nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "T2D"
  }
  else{
    GWAS_hit_meta[i] = "no"
  }
}
nodes$GWAS_hit_meta = GWAS_hit_meta

# do again including brain mapping
GWAS_hit_meta = c()
for(i in 1:length(nodes$NODE)){
  if(nodes$GWAS_hit[i]!="no" & nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "AD_T2D"
  }
  else if(nodes$GWAS_hit[i]!="no"){
    GWAS_hit_meta[i] = "AD"
  }
  else if(nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "T2D"
  }
  else{
    GWAS_hit_meta[i] = "no"
  }
}
nodes$GWAS_hit_meta = GWAS_hit_meta


# make new column that is for sizes - key drivers are biggest, then GWAS genes, then non-GWAS genes

new_size = c()
for(i in 1:length(nodes$NODE)){
  if(nodes$SHAPE[i]=="Diamond"){
    new_size[i] = "Big"
  }
  else if(nodes$GWAS_hit_meta[i]!="no"){
    new_size[i] = "Medium"
  }
  else{
    new_size[i] = "Small"
  }
}
nodes$New_Size = new_size
write.table(nodes, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_3.txt", 
            row.names = FALSE, quote = FALSE, sep="\t")



GWAS_hit_meta = c()
for(i in 1:length(including_Brain$NODE)){
  if(sum(including_Brain[i,c(9,10,11,12)]!="no")>0){
    if(sum(including_Brain[i,c(9,10,11,12)]!="no")>3){
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(including_Brain$GWAS_hit[i]!="no" & including_Brain$GWAS_brain_T2D[i]!="no"){ # matching brain
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(including_Brain$GWAS_perip_AD[i]!="no" & including_Brain$GWAS_periph_T2D[i]!="no"){ # matching perip
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(including_Brain$GWAS_perip_AD[i]!="no" & including_Brain$GWAS_brain_T2D[i]!="no"){ # not match
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(including_Brain$GWAS_hit[i]!="no" & including_Brain$GWAS_periph_T2D[i]!="no"){ # not match
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(including_Brain$GWAS_hit[i]!="no" & including_Brain$GWAS_periph_T2D[i]!="no"){ # not match
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(sum(including_Brain[i,c(9,10,11,12)]!="no")==2){
      if(including_Brain$GWAS_perip_AD[i]!="no" & including_Brain$GWAS_hit[i]!="no"){ # not match
        GWAS_hit_meta[i] = "AD"
      }
      else if(including_Brain$GWAS_brain_T2D[i]!="no" & including_Brain$GWAS_periph_T2D[i]!="no"){ 
        GWAS_hit_meta[i] = "T2D"
      }
      else{
        GWAS_hit_meta[i] = "warning"
        cat(i, " Warning\n")
      }
    }
    else if(sum(including_Brain[i,c(9,10,11,12)]!="no")==1){
      if(including_Brain$GWAS_hit[i]!="no"){
        GWAS_hit_meta[i] = "AD"
      }
      else if(including_Brain$GWAS_perip_AD[i]!="no"){
        GWAS_hit_meta[i] = "AD"
      }
      else if(including_Brain$GWAS_brain_T2D[i]!="no"){
        GWAS_hit_meta[i] = "T2D"
      }
      else if(including_Brain$GWAS_periph_T2D[i]!="no"){
        GWAS_hit_meta[i] = "T2D"
      }
      else{
        GWAS_hit_meta[i] = "warning"
        cat(i, " Warning\n")
      }
    }
    else{
      GWAS_hit_meta[i] = "warning"
      cat(i, " Warning\n")
    }
  }
  else{
    GWAS_hit_meta[i] = "no"
  }
  
}
including_Brain$GWAS_hit_meta = GWAS_hit_meta
new_size = c()
for(i in 1:length(including_Brain$NODE)){
  if(including_Brain$SHAPE[i]=="Diamond"){
    new_size[i] = "Big"
  }
  else if(including_Brain$GWAS_hit_meta[i]!="no"){
    new_size[i] = "Medium"
  }
  else{
    new_size[i] = "Small"
  }
}
including_Brain$New_Size = new_size
write.table(including_Brain, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_4.txt", 
            row.names = FALSE, quote = FALSE, sep="\t")



# refine brain version 2... ------
nodes <- add_GWAS_hit_info_network(nodes = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network2/cytoscape/kda2cytoscape.nodes.txt",
                                   SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/sig_genes/", study = "AD")
nodes3 <- add_GWAS_hit_info_network(nodes = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network2/cytoscape/kda2cytoscape.nodes.txt",
                                    SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/sig_genes/", study = "T2D")
nodes$GWAS_T2D = nodes3$GWAS_hit

GWAS_hit_meta = c()
for(i in 1:length(nodes$NODE)){
  if(nodes$GWAS_hit[i]!="no" & nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "AD_T2D"
  }
  else if(nodes$GWAS_hit[i]!="no"){
    GWAS_hit_meta[i] = "AD"
  }
  else if(nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "T2D"
  }
  else{
    GWAS_hit_meta[i] = "no"
  }
}
nodes$GWAS_hit_meta = GWAS_hit_meta
write.table(nodes, "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network2/cytoscape/mod_nodes.txt", 
            row.names = FALSE, quote = FALSE, sep="\t")
brain <- trim_network(folder = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network2/cytoscape/", 
                      keep_only_orig_input = TRUE,keep_only_GWAS_hits = TRUE) # need to move other nodes file to another folder
# OR for this module, only keep GWAS hits
remove_modules <- c("rctm0522,.._AD","rctm0917,.._AD")
# HIV infection, Protein folding
# take out M963 entirely? RNA degradation
# take out rctm1360 entirely - OXPHOS (known)



keep <- c("GGH", "GSTM1","GSTP1", "GGH","TUBB4B","CCT2","CCT8")


# account for unique modules for Adipose, Brain, Artery, Colon, and Esophagus (because just added modules up) ----
tissues <- c("Adipose","Brain","Artery","Colon","Esophagus", "Skin", "Heart")
traits <- c("AD", "T2D", "Shared")
folders = c("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/","/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D.2/",
            "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/")

result_list = c()
for(i in 1:length(traits)){
  vector = c()
  folder = folders[i]
  for(t in tissues){
    files = list.files(folder)[intersect(grep("info",list.files(folder)), grep("merged",list.files(folder)))]
    files = files[intersect(grep("eQTL",files),grep(t,files))]
    if(length(files)==0) next
    modules = c()
    for(f in files){
      df = read.delim(paste0(folder, f), stringsAsFactors = FALSE)
      modules = append(modules, df$MODULE)
    }
    vector[t] = length(unique(modules))
  }
  result_list[[traits[i]]] = vector
}

folder = folders[3]
files = list.files(folder)[intersect(grep("info",list.files(folder)), grep("merged",list.files(folder)))]
files = files[grep("sQTL",files)]
for(f in files){
  df = read.delim(paste0(folder, f), stringsAsFactors = FALSE)
  modules = df$MODULE
  cat(gsub("merged_","",gsub(".info.txt","",f)), ": ", length(unique(modules)),"\n")
}

files = list.files(folder)[intersect(grep("mod",list.files(folder)), grep("merged",list.files(folder)))]
files = files[intersect(grep("Brain",files),grep("sQTL",files))]
gene_list = c()
for(f in files){
  df = read.delim(paste0(folder, f), stringsAsFactors = FALSE)
  gene_list = append(gene_list, df$GENE)
}
gene_list = unique(gene_list)


nodes = read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_4.txt",
                   stringsAsFactors = FALSE)
brain_member = c()
for(n in nodes$NODE){
  if(sum(gene_list==n)>0){
    brain_member = append(brain_member,"YES")
  }
  else{
    brain_member = append(brain_member,"NO")
  }
}
nodes$brain_member = brain_member
write.table(nodes, "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_5.txt",
            row.names = FALSE, quote = FALSE, sep="\t")
edges <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_kda2cytoscape.edges_2.txt",
                    stringsAsFactors = FALSE)

brain_member_edge = c()
for(i in 1:length(edges$TAIL)){
  if(sum(gene_list==edges$TAIL[i])>0){
    brain_member_edge = append(brain_member_edge, "YES")
  }
  else if(sum(gene_list==edges$HEADL[i])>0){
    brain_member_edge = append(brain_member_edge, "YES")
  }
  else{
    brain_member_edge = append(brain_member_edge, "NO")
  }
}
edges$Brain_member = brain_member_edge
write.table(edges, "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/edges_5.txt",
            row.names = FALSE, quote = FALSE, sep="\t")
write.table(edges[,6], "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/edges_6.txt",
            row.names = FALSE, quote = FALSE, sep="\t")

# Annotate significant coexpression modules -------
# get significant modules - this does not include the quasi overlapped modules.
# does not annotate Zhang et al. brain modules
source("/Users/jessicading/Desktop/Yang_Lab/pathwayEnrichment.R")
library(WriteXLS)
setwd("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/Supplement/annotate_pathway/") # so that all the files go here
sig_coexpr_modules = append(AD_replication$MODULE[grep("GTEXv7",AD_replication$MODULE)],
                            T2D_replication$MODULE[grep("GTEXv7",T2D_replication$MODULE)])
sig_coexpr_modules = append(sig_coexpr_modules, AD_T2D_replication$MODULE[grep("GTEXv7",AD_T2D_replication$MODULE)]) # redundant but just make sure to get all of them
sig_coexpr_modules = unique(sig_coexpr_modules)
module_df = read.delim("/Users/jessicading/Desktop/Yang_Lab/Resources/genesets/all.mod.txt", stringsAsFactors = FALSE)
sig_coexpr_modules_df = data.frame(stringsAsFactors = FALSE)
for(mod in sig_coexpr_modules){
  temp = module_df[module_df$MODULE==mod,]
  sig_coexpr_modules_df = rbind(sig_coexpr_modules_df, temp) # why is this turning out to be more?
}
# change module names for readability and also because sheet names cannot be more than 31 characters
sig_coexpr_modules_df$MODULE <- gsub(".GTEXv7","",sig_coexpr_modules_df$MODULE)
sig_coexpr_modules_df$MODULE <- gsub("WGCNA","W", sig_coexpr_modules_df$MODULE)
sig_coexpr_modules_df$MODULE <- gsub("MEGENA","M", sig_coexpr_modules_df$MODULE)
sig_coexpr_modules_df$MODULE <- gsub("Brain_","", sig_coexpr_modules_df$MODULE)

pathwayEnrichment(modules_list = sig_coexpr_modules_df, resources_path = "/Users/jessicading/Desktop/Yang_Lab/Resources")



# redo T2D.......
folder = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/"
files = list.files(folder)[grep("shared",list.files(folder))]
for(i in files){
  name = gsub("shared_AD.", "", i)
  name = gsub(".txt", "", name)
  modules_df = read.delim(paste0(folder, i), header = TRUE, stringsAsFactors = FALSE)
  if(length(modules_df$MODULE)==0) next
  if(length(modules_df$MODULE)==1){
    write.table(modules_df, paste0(folder, "merged_",name,".info.txt"), row.names = FALSE, quote = FALSE, sep="\t")
    cat(name, " had only 1 significant module. Making fake merged file.\n")
    next
  }
  merge_modules(name = name, 
                rcutoff = 0.30,
                modules_df = modules_df, 
                output_Dir = folder, 
                modfile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.mod.txt", 
                infofile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.info.txt")
}




# change node file for combined peripheral ssnp and eqtl to make GWAS hit less than 0.001... see how that looks...
for(iter in 1:2){
  createGenesFiles(results_Dir = paste0(output_Dirs[iter],"periph_brain_genes/"), 
                   genes_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/all_genes/", 
                   output_Dir = paste0(output_Dirs[iter], "sig_genes/"), # make the genes folder
                   studies = cohorts[[iter]],
                   value_cutoff = 1.3, 
                   modfile = "~/Desktop/Yang_Lab/Resources/genesets/all.mod.txt") # from above geneOverlap
}

nodes <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_6.txt",
                                   SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/p0.001_sig_genes/", study = "AD")
nodes2 <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_6.txt", 
                                   SGM = "peripheral_ssnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/p0.001_sig_genes/", study = "AD")
nodes3 <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_6.txt",
                                    SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/p0.001_sig_genes/", study = "T2D")
nodes4 <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_6.txt", 
                                    SGM = "peripheral_ssnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/p0.001_sig_genes/", study = "T2D")
nodes$GWAS_perip_AD = nodes2$GWAS_hit
nodes$GWAS_brain_T2D = nodes3$GWAS_hit
nodes$GWAS_periph_T2D = nodes4$GWAS_hit
including_Brain = nodes

GWAS_hit_meta = c()
for(i in 1:length(including_Brain$NODE)){
  if(sum(including_Brain[i,c(9,10,11,12)]!="no")>0){
    if(sum(including_Brain[i,c(9,10,11,12)]!="no")>3){
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(including_Brain$GWAS_hit[i]!="no" & including_Brain$GWAS_brain_T2D[i]!="no"){ # matching brain
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(including_Brain$GWAS_perip_AD[i]!="no" & including_Brain$GWAS_periph_T2D[i]!="no"){ # matching perip
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(including_Brain$GWAS_perip_AD[i]!="no" & including_Brain$GWAS_brain_T2D[i]!="no"){ # not match
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(including_Brain$GWAS_hit[i]!="no" & including_Brain$GWAS_periph_T2D[i]!="no"){ # not match
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(including_Brain$GWAS_hit[i]!="no" & including_Brain$GWAS_periph_T2D[i]!="no"){ # not match
      GWAS_hit_meta[i] = "AD_T2D"
    }
    else if(sum(including_Brain[i,c(9,10,11,12)]!="no")==2){
      if(including_Brain$GWAS_perip_AD[i]!="no" & including_Brain$GWAS_hit[i]!="no"){ # not match
        GWAS_hit_meta[i] = "AD"
      }
      else if(including_Brain$GWAS_brain_T2D[i]!="no" & including_Brain$GWAS_periph_T2D[i]!="no"){ 
        GWAS_hit_meta[i] = "T2D"
      }
      else{
        GWAS_hit_meta[i] = "warning"
        cat(i, " Warning\n")
      }
    }
    else if(sum(including_Brain[i,c(9,10,11,12)]!="no")==1){
      if(including_Brain$GWAS_hit[i]!="no"){
        GWAS_hit_meta[i] = "AD"
      }
      else if(including_Brain$GWAS_perip_AD[i]!="no"){
        GWAS_hit_meta[i] = "AD"
      }
      else if(including_Brain$GWAS_brain_T2D[i]!="no"){
        GWAS_hit_meta[i] = "T2D"
      }
      else if(including_Brain$GWAS_periph_T2D[i]!="no"){
        GWAS_hit_meta[i] = "T2D"
      }
      else{
        GWAS_hit_meta[i] = "warning"
        cat(i, " Warning\n")
      }
    }
    else{
      GWAS_hit_meta[i] = "warning"
      cat(i, " Warning\n")
    }
  }
  else{
    GWAS_hit_meta[i] = "no"
  }
  
}
including_Brain$GWAS_hit_meta = GWAS_hit_meta 
new_size = c()
for(i in 1:length(including_Brain$NODE)){
  if(including_Brain$SHAPE[i]=="Diamond"){
    new_size[i] = "Big"
  }
  else if(including_Brain$GWAS_hit_meta[i]!="no"){
    new_size[i] = "Medium"
  }
  else{
    new_size[i] = "Small"
  }
}
including_Brain$New_Size = new_size
write.table(including_Brain, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_7.txt", 
            row.names = FALSE, quote = FALSE, sep="\t")

# include nodes file that includes both GWAS hit info from 0.05 and 0.001
nodes_p0.05 <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_5.txt",
                          stringsAsFactors = FALSE)

including_Brain$GWAS_hit_meta_0.001 = including_Brain$GWAS_hit_meta
including_Brain$GWAS_hit_meta_0.05 = nodes_p0.05$GWAS_hit_meta
# including_Brain$GWAS_hit_meta will just be for the color

GWAS_significance = c()
for(i in 1:length(including_Brain$NODE)){
  if(including_Brain$GWAS_hit_meta_0.001[i]!="no"){
    GWAS_significance[i] = "0.001"
  }
  else if(including_Brain$GWAS_hit_meta_0.001[i]=="no" & including_Brain$GWAS_hit_meta_0.05[i]!="no"){
    GWAS_significance[i] = "0.05"
  }
  else{
    GWAS_significance[i] = "no"
  }
}
including_Brain$GWAS_significance = GWAS_significance


new_size = c()
for(i in 1:length(including_Brain$NODE)){
  if(including_Brain$SHAPE[i]=="Diamond"){
    new_size[i] = "Big"
  }
  else if(including_Brain$GWAS_hit_meta_0.001[i]!="no"){
    new_size[i] = "Medium Big"
  }
  else if(including_Brain$GWAS_hit_meta_0.05[i]!="no"){
    new_size[i] = "Medium"
  }
  else{
    new_size[i] = "Small"
  }
}
including_Brain$Newest_Size = new_size

including_Brain$GWAS_hit_meta = nodes_p0.05$GWAS_hit_meta # just for the color, based on p<0.05

write.table(including_Brain, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_8.txt", 
            row.names = FALSE, quote = FALSE, sep="\t")


# now change the brain network one....

brain_nodes <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network2/cytoscape/mod_nodes_2.txt",
                                   SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/p0.001_sig_genes/", study = "AD")
brain_nodes3 <- add_GWAS_hit_info_network(nodes = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network2/cytoscape/mod_nodes_2.txt",
                                    SGM = "brain_esnps", genes_folder = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/p0.001_sig_genes/", study = "T2D")
brain_nodes$GWAS_T2D = brain_nodes3$GWAS_hit

GWAS_hit_meta = c()
for(i in 1:length(brain_nodes$NODE)){
  if(brain_nodes$GWAS_hit[i]!="no" & brain_nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "AD_T2D"
  }
  else if(brain_nodes$GWAS_hit[i]!="no"){
    GWAS_hit_meta[i] = "AD"
  }
  else if(brain_nodes$GWAS_T2D[i]!="no"){
    GWAS_hit_meta[i] = "T2D"
  }
  else{
    GWAS_hit_meta[i] = "no"
  }
}
brain_nodes$GWAS_hit_meta_p0.001 = GWAS_hit_meta
orig_brain_nodes <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network2/cytoscape/mod_nodes.txt",
                               stringsAsFactors = FALSE)
brain_nodes$GWAS_hit_meta_p0.05 = orig_brain_nodes$GWAS_hit_meta
write.table(brain_nodes,"~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network2/cytoscape/mod_nodes_2.txt",
            row.names = FALSE, quote = FALSE, sep="\t")
new_size = c()
for(i in 1:length(brain_nodes$NODE)){
  if(brain_nodes$SHAPE[i]=="Diamond"){
    new_size[i] = "Big"
  }
  else if(brain_nodes$GWAS_hit_meta_p0.001[i]!="no"){
    new_size[i] = "Medium Big"
  }
  else if(brain_nodes$GWAS_hit_meta_p0.05[i]!="no"){
    new_size[i] = "Medium"
  }
  else{
    new_size[i] = "Small"
  }
}
brain_nodes$Newest_Size = new_size

brain_nodes$GWAS_hit_meta = orig_brain_nodes$GWAS_hit_meta # just for the color, based on p<0.05
write.table(brain_nodes,"~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network2/cytoscape/mod_nodes_2.txt",
            row.names = FALSE, quote = FALSE, sep="\t")

# count number of genes where the node has edges from both AD and T2D
edges <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/brain_AD_T2D_string_network2/cytoscape/kda2cytoscape.edges.txt", 
                    stringsAsFactors = FALSE)
# get all nodes
nodes <- append(edges$TAIL, edges$HEAD)
nodes = unique(nodes)
connector_nodes = c()
for(n in nodes){
  # get edges - both in head and tail
  subnet = rbind(edges[edges$TAIL==n,], edges[edges$HEAD==n,])
  if(sum(grepl("AD",subnet$MODULE))>0 & sum(grepl("T2D",subnet$MODULE))>0){
    connector_nodes = append(connector_nodes, n)
  }
  else{
    next
  }
}

# these connector nodes probably have membership between AD and T2D modules
# for each of these nodes, get membership

# make version of DESCR.txt (coexpression annotation) that only has the modules that were shared between T2D and AD



# look at GWAS genes from GWAS Catalog from AD and T2D
AD_GWAS <- read.delim("~/Downloads/gwas-association-downloaded_2019-11-02-EFO_0000249-withChildTraits.tsv",
                      stringsAsFactors = FALSE)
T2D_GWAS <- read.delim("~/Downloads/gwas-association-downloaded_2019-10-28-EFO_0001360-withChildTraits.tsv",
                       stringsAsFactors = FALSE)
AD_GWAS = unique(AD_GWAS$REPORTED.GENE.S.)
AD_GWAS_new = c()
for(item in 1:length(AD_GWAS)){
  if(grepl(",",AD_GWAS[item])){
    splitted = unlist(strsplit(AD_GWAS[item], split = ", "))
    AD_GWAS_new = append(AD_GWAS_new, splitted)
  }
  else{
    AD_GWAS_new = append(AD_GWAS_new, AD_GWAS[item])
  }
}
AD_dup <- names(table(AD_GWAS_new))[table(AD_GWAS_new)>1]
AD_GWAS_new = unique(AD_GWAS_new)

T2D_GWAS = unique(T2D_GWAS$REPORTED.GENE.S.)
T2D_GWAS_new = c()
for(item in 1:length(T2D_GWAS)){
  if(grepl(",",T2D_GWAS[item])){
    splitted = unlist(strsplit(T2D_GWAS[item], split = ", "))
    T2D_GWAS_new = append(T2D_GWAS_new, splitted)
  }
  else{
    T2D_GWAS_new = append(T2D_GWAS_new, T2D_GWAS[item])
  }
}
T2D_dup <- names(table(T2D_GWAS_new))[table(T2D_GWAS_new)>1]
T2D_GWAS_new = unique(T2D_GWAS_new)

shared_AD_T2D_GWAS <- intersect(AD_GWAS_new, T2D_GWAS_new)

write.table(shared_AD_T2D_GWAS,"~/Downloads/shared_gwas_ad_t2d.txt", quote = FALSE,row.names = FALSE, sep = "\t")


addColumn <- function(column, df){
  vect = c()
  for(id in df$projid){
    vect = append(vect, clinical[,column][clinical$projid==id])
  }
  return(vect)
}




# 12/2/2019
# Problem: Annotation of coexpression modules but some annotated pathways might NOT have "overlaps" with the GWAS suggestive genes
# Solution: look through pathway enrichment files - only keep the annotations that have overlaps with GWAS suggestive genes
# need: 1) genes folder for AD/T2D 2) pathway enrichment files 
# make data frame: for every, coexpression module, if it was significantly enriched, what genes in the overlap column
# were actually gwas suggestive genes - do for every study
# annotation has to have at least one gene from one study
# COEXPRESSION_MODULE       ANNOTATION       UKB_Overlap        DIAGRAM_Overlap    etc...
genes_folder = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/all_genes"
pathways_folder = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/pathways/"
pathway_files = list.files(pathways_folder, pattern = ".txt")
pathway_files = pathway_files[!grepl("DESCR", pathway_files)]
pathway_files = pathway_files[!grepl("tmp", pathway_files)]

total = data.frame(stringsAsFactors = FALSE)
for(p in pathway_files){
  path_df = read.delim(paste0(pathways_folder, p), stringsAsFactors = FALSE)
  path_df = path_df[path_df$nOverlap>1,]
  path_name = gsub(".txt","",p)
  # look through details files if it was a significant pathway FDR<0.05
  traits = c()
  for(f in list.files(genes_folder)){
    traits = append(traits, unlist(strsplit(f, split = ".",fixed = TRUE))[1]) 
  }
  traits = unique(traits)
  for(annotation in path_df$Pathway){
    temp = data.frame("MODULE"=p, 
                      "Annotation"=annotation,stringsAsFactors = FALSE)
    for(t in traits){
      t_files = list.files(genes_folder)[intersect(grep(t, list.files(genes_folder)), grep("eQTL",list.files(genes_folder)))]
      t_files = t_files[grep(unlist(strsplit(path_name, split = "."))[1], t_files)]
      for(g in t_files){
        df <- read.delim(paste0(genes_folder, "/",g), stringsAsFactors = FALSE)
        df$FDR <- as.numeric(gsub("%","",df$FDR))
        if(is.na(df$FDR[1])){
          temp[,t] = "None"
          next
        }
        df = df[df$FDR<=5,]
        if(nrow(df)==0){
          temp[,t] = "None"
          next
        }
        if(sum(df$MODULE==path_name)>0){
          genes = df$GENE[df$MODULE==path_name]
          single_genes = c()
          for(gene in genes){
            splitted = unlist(strsplit(gene, split = ","))
            single_genes = append(single_genes, splitted)
          }
            path_genes = unlist(strsplit(path_df$Overlap[path_df$Pathway==annotation], split = ","))
            if(length(intersect(path_genes, single_genes))>0){
              temp[,t] = concatenate(intersect(path_genes, single_genes), mysep = ", ")
            }
            else{
              temp[,t] = "None"
            }
          }
        }
    }
    total = rbind(total, temp)
  }
 
}



# DRUG REPOSITIONING -------
# from the peripheral network, get all genes
nodes <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/peripheral_ssnps_sQTL_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt", 
                    stringsAsFactors = FALSE)
nodes = nodes[nodes$SHAPE=="Diamond",]


# run KDA on unique AD modules
AD_periph <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/shared_AD.peripheral_esnps_eQTL.txt",
                        stringsAsFactors = FALSE)

T2D_periph <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D/shared_T2D_quasi.peripheral_esnps_eQTL.txt",
                         stringsAsFactors = FALSE)

shared_AD_T2D_with_quasi <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/shared_AD_T2D_quasi.peripheral_esnps_eQTL.txt",
                                       stringsAsFactors = FALSE)
shared_AD_T2D_periph <- intersect(AD_periph$DESCR, T2D_periph$DESCR)

AD_periph_specific = setdiff(AD_periph$MODULE, T2D_periph$MODULE)
AD_periph_specific = setdiff(AD_periph_specific, shared_AD_T2D_with_quasi$MODULE)
T2D_periph_specific <-setdiff(T2D_periph$MODULE, AD_periph$MODULE)
T2D_periph_specific = setdiff(T2D_periph_specific, shared_AD_T2D_with_quasi$MODULE)

AD_periph_specific_df = data.frame("MODULE" = AD_periph_specific)
T2D_periph_specific_df = data.frame("MODULE" = T2D_periph_specific)

modfile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.mod.txt"
infofile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.info.txt"

merge_modules(name = "AD_periph_specific", modules_df = AD_periph_specific_df, rcutoff = .30,
              output_Dir = "./", modfile_path = modfile_path, infofile_path = infofile_path)
merge_modules(name = "T2D_periph_specific", modules_df = T2D_periph_specific_df, rcutoff = .30, 
              output_Dir = "./",  modfile_path = modfile_path, infofile_path = infofile_path)
merge_modules(name = "AD_T2D_shared_eQTL", modules_df = data.frame("MODULE"=shared_AD_T2D_with_quasi$MODULE), rcutoff = .30, 
              output_Dir = "./",  modfile_path = modfile_path, infofile_path = infofile_path)
merge_modules(name = "AD_selected", modules_df = data.frame("MODULE"=c("rctm0917","M19540","rctm1050","rctm0250","M14933","rctm0647")), rcutoff = .30, 
              output_Dir = "./",  modfile_path = modfile_path, infofile_path = infofile_path)
runKDA(nodes = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/committee_meeting/merged_AD_periph_specific.mod.txt", 
       network = "/Users/jessicading/Desktop/Yang_Lab/Resources/network/AD_T2D/pan_esnp/networks.hs.all.txt", trim = c(".mod.txt",".txt"))
runKDA(nodes = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/committee_meeting/merged_T2D_periph_specific.mod.txt", 
       network = "/Users/jessicading/Desktop/Yang_Lab/Resources/network/AD_T2D/pan_esnp/networks.hs.all.txt", trim = c(".mod.txt",".txt"))
runKDA(nodes = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/committee_meeting/merged_AD_T2D_shared_eQTL.mod.txt", 
       network = "/Users/jessicading/Desktop/Yang_Lab/Resources/network/AD_T2D/pan_esnp/networks.hs.all.txt", trim = c(".mod.txt",".txt"))
runKDA(nodes = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/committee_meeting/merged_AD_selected.mod.txt", 
       network = "/Users/jessicading/Desktop/Yang_Lab/Resources/network/AD_T2D/pan_esnp/networks.hs.all.txt", trim = c(".mod.txt",".txt"))

summarizeKDA(KDA_folder = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/committee_meeting/network_results/", 
             protein_descrip = "/Users/jessicading/Desktop/Yang_Lab/Epicardial/String/9606.protein.info.v11.0.txt", 
             name = "Peripheral_eQTL_AD_T2D")

KDA_AD_specific = read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/committee_meeting/network_results/merged_AD_periph_specific_networks.hs.all/kda/wKDA.results.txt",
                             stringsAsFactors = FALSE)
KDA_AD_specific = KDA_AD_specific[KDA_AD_specific$FDR<0.05,]
KDA_T2D_specific = read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/committee_meeting/network_results/merged_T2D_periph_specific_networks.hs.all/kda/wKDA.results.txt",
                             stringsAsFactors = FALSE)
KDA_T2D_specific = KDA_T2D_specific[KDA_T2D_specific$FDR<0.05,]
KDA_AD_T2D = read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/committee_meeting/network_results/merged_AD_T2D_shared_eQTL_networks.hs.all/kda/wKDA.results.txt",
                             stringsAsFactors = FALSE)
KDA_AD_T2D = KDA_AD_T2D[KDA_AD_T2D$FDR<0.05,]

AD_KDs = setdiff(KDA_AD_specific$NODE, KDA_AD_T2D$NODE)
AD_KDs = setdiff(AD_KDs, KDA_T2D_specific$NODE)

T2D_KDs = setdiff(KDA_T2D_specific$NODE, KDA_AD_T2D$NODE)
T2D_KDs = setdiff(T2D_KDs, KDA_AD_specific$NODE)

AD_T2D_KDs = setdiff(KDA_AD_T2D$NODE, KDA_T2D_specific$NODE)
AD_T2D_KDs = setdiff(AD_T2D_KDs, KDA_AD_specific$NODE)

network_summary = read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/committee_meeting/Peripheral_eQTL_AD_T2D_KDs_summary.txt",
                             stringsAsFactors = FALSE)

AD_KDA = data.frame(stringsAsFactors = FALSE)
for(kd in AD_KDs){
  AD_KDA = rbind(AD_KDA, network_summary[network_summary$KDs==kd & network_summary$NETWORK=="merged_AD_periph_specific_networks.hs.all",])
}

AD_KDA = AD_KDA[!duplicated(AD_KDA$KDs),]
write.table(AD_KDA, "AD_KDA.txt", row.names = FALSE, quote = FALSE, sep = "\t")

T2D_KDA = data.frame(stringsAsFactors = FALSE)
for(kd in T2D_KDs){
  T2D_KDA = rbind(T2D_KDA, network_summary[network_summary$KDs==kd & network_summary$NETWORK=="merged_T2D_periph_specific_networks.hs.all",])
}

AD_T2D_KDA = data.frame(stringsAsFactors = FALSE)
for(kd in AD_T2D_KDs){
  AD_T2D_KDA = rbind(AD_T2D_KDA, network_summary[network_summary$KDs==kd & network_summary$NETWORK=="merged_AD_T2D_shared_eQTL_networks.hs.all",])
}


# change network 
nodes_descrip <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_8.txt",
                            stringsAsFactors = FALSE)
url = "http://chart.apis.google.com/chart?cht=p&chs=200x200&chf=bg,s,00000000&chd=t:"
AD_T2D = "1,1&chco=A4CB97|A9CDE5"
AD = "1&chco=A4CB97"
T2D = "1&chco=A9CDE5"

URL2 = c()
GWAS_hit = c()
for(node in 1:nrow(nodes_descrip)){
  if(nodes_descrip$GWAS_hit_meta[node]=="T2D"){
    URL2[node] = paste0(url, T2D)
  }
  else if(nodes_descrip$GWAS_hit_meta[node]=="AD"){
    URL2[node] = paste0(url, AD)
  }
  else if(nodes_descrip$GWAS_hit_meta[node]=="AD_T2D"){
    URL2[node] = paste0(url, AD_T2D)
  }
  else{
    URL2[node] = ""
  }
  if(sum(nodes_descrip$NODE[node]==AD_GWAS_new)>0 & sum(nodes_descrip$NODE[node]==T2D_GWAS_new)>0){
    GWAS_hit[node] = "AD_T2D"
  }
  else if(sum(nodes_descrip$NODE[node]==AD_GWAS_new)>0){
    GWAS_hit[node] = "AD"
  }
  else if(sum(nodes_descrip$NODE[node]==T2D_GWAS_new)>0){
    GWAS_hit[node] = "T2D"
  }
  else{
    GWAS_hit[node] = "no"
  }
}

nodes_descrip$URL2 = URL2
nodes_descrip$GWAS_actual = GWAS_hit

write.table(nodes_descrip, "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/results/combined_periph_ssnp_brain_esnp_AD_T2D_networks.hs.all/cytoscape/mod_nodes_8.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# 12/10/19 For paper ----
# merge all significant modules to create independent supersets
all_sig_mods = unique(c(AD_T2D, AD, T2D))
merge_modules(name = "All_mods", modules_df = data.frame("MODULE"=all_sig_mods), rcutoff = .30, 
              output_Dir = "/Users/jessicading/Desktop/Yang_Lab/T2D_AD/Paper/", 
              modfile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.mod.txt", 
              infofile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.info.txt")
All_mods_mod <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/All_mods.mod.txt", stringsAsFactors = FALSE)
# 628 sig mods merged into 256 modules
# create superset info

# first want to see modules specific to combined, then combined+peripheral, then just peripheral. 
# (to see if I want to redo peripheral to not include nerve tibial) (decide if i should even consider peripheral)
# just peripheral does have interesting modules like regulation of insulin secretion, etc.


# 1/14/20 
# annotation of significant modules is here:
# also want to annotate supersets 

# split supersets by coexpression vs canonical 
supersets <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/All_mods.mod.txt", stringsAsFactors = FALSE)
coexpr_supersets <- supersets[grep("GTEXv7", supersets$MODULE),]
canon_supersets <- supersets[!grepl("GTEXv7", supersets$MODULE),]

# get all colors (Zhang modules)
zhang_modules <- supersets$MODULE[grepl(" ", supersets$MODULE)]
# only hot pink. Gold, light cyan, etc. were merged in a very large superset

supersets_trim <- supersets[unique(c(grep("GTEX",supersets$MODULE), grep(",..", supersets$MODULE), grep("pink", supersets$MODULE))),]
supersets_trim$MODULE <- gsub(" ","_",supersets_trim$MODULE)
supersets_trim$MODULE <- gsub(".GTEXv7","",supersets_trim$MODULE)
supersets_trim$MODULE <- gsub("WGCNA","W", supersets_trim$MODULE)
supersets_trim$MODULE <- gsub("MEGENA","M", supersets_trim$MODULE)
supersets_trim$MODULE <- gsub("Brain_","", supersets_trim$MODULE)
supersets_trim$MODULE <- gsub("Visceral_Omentum","Visc_Oment", supersets_trim$MODULE)


source("/Users/jessicading/Desktop/Yang_Lab/Resources/genesets/R-tomfunctions.R")
setwd("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/annotate_supersets")
pathwayEnrichment(modules_list = supersets_trim[,c("MODULE","GENE")], resources_path = "/Users/jessicading/Desktop/Yang_Lab/Resources", name = "Annotate")

extract_annotationv2()


makeNewGeneSetFileOnlySNPMappedGenes <- function(geneset_file, details_file, return_only_mod_genes=TRUE){
  details <- read.delim(details_file, header=TRUE, stringsAsFactors = FALSE)
  geneset <- read.delim(geneset_file, header=TRUE, stringsAsFactors = FALSE)
  details = details[!(grepl("_ctrl", details$MODULE)),]
  result_all = data.frame(stringsAsFactors = FALSE)
  result = data.frame(stringsAsFactors = FALSE)
  for(mod in unique(details$MODULE)){
    details %>% 
      filter(MODULE==mod) %>% 
      pull(GENE) -> genes
    multiple_genes <- genes[grep(",", genes)]
    split_genes <- unlist(sapply(multiple_genes, 
                                 FUN = function(x){return(unlist(strsplit(x,
                                                                          split = ",")))}))
    trimmed_genes <- genes[!grepl(",", genes)]
    genes <- c(trimmed_genes, split_genes)
    if(!return_only_mod_genes){
      result_all = rbind(result_all, 
                         data.frame("MODULE"=mod, 
                                    "GENE"=genes, stringsAsFactors = FALSE))
    }
    toTrim = c()
    for(iter in 1:length(genes)){
      if(!(genes[iter] %in% geneset$GENE[geneset$MODULE==mod])) toTrim = c(toTrim,iter)
    }
    if(length(toTrim)>0) genes = genes[-toTrim]
    if(return_only_mod_genes){
      result = rbind(result, 
                    data.frame("MODULE"=mod, 
                               "GENE"=genes, stringsAsFactors = FALSE))
    }
  }
  if(!return_only_mod_genes){
    write.table(result_all,"New_GeneSets_SNP_mapped_all.txt", sep="\t", quote = FALSE, row.names = FALSE)
  }
  else{
    write.table(result,"New_GeneSets_SNP_mapped_t2.txt", sep="\t", quote = FALSE, row.names = FALSE)
  }
}

makeNewGeneSetFileOnlySNPMappedGenes(geneset_file = "~/Desktop/Yang_Lab/Resources/genesets/all.mod.txt",
                                     details_file = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/all_genes/DIAGRAMstage1_T2D.Adipose_Subcutaneous_eQTL.50.20.details.txt")
look <- read.delim("/Users/jessicading/Desktop/10X_Analysis/AD_Metabolic/New_GeneSets_SNP_mapped.txt", stringsAsFactors = FALSE)
orig_file <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/all_genes/DIAGRAMstage1_T2D.Adipose_Subcutaneous_eQTL.50.20.details.txt")

# May 28, 2020

# specific AD pathways ---------
AD_pathways <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Supplement/AD_all_pathways.txt")
T2D_pathways <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Supplement/T2D_all_pathways.txt")
AD_specific = setdiff(unique(AD_pathways$MODULE), unique(T2D_pathways$MODULE))
AD_specific_d = setdiff(unique(AD_pathways$DESCR), unique(T2D_pathways$DESCR))

# merge AD specific
AD_df = makeDf(Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/") # takes shared files from the quasi overlap
infofile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.info.txt"
all_DESCR <- read.delim(infofile_path, stringsAsFactors = FALSE)

# very sig is a value cut off of 0.01. not very high
# make list of AD-specific pathways and mappings, same for T2D - eQTL specific!
AD_specific_df = describePathways(modules_list = AD_specific, df = AD_df, DESCR = all_DESCR)
AD_specific_df <- addGenes(df = AD_specific_df, genes_files = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/very_sig_genes/")

# how many s

# look at all diagram studies
results_files <- list.files("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/all")
results <- results_files[grep("DIAGRAM", results_files)]
all_DIAGRAM = data.frame(stringsAsFactors = FALSE)
for(file in results){
  temp = read.delim(paste0("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/all/",file), stringsAsFactors = FALSE)
  all_DIAGRAM = rbind(all_DIAGRAM, temp[temp$FDR<0.05,])
}
intersect(unique(all_DIAGRAM$MODULE), unique(AD_specific_df$MODULE))

truly_unique_AD = setdiff(unique(AD_specific_df$MODULE), unique(all_DIAGRAM$MODULE))
truly_unique_AD_df = AD_specific_df[AD_specific_df$MODULE %in% truly_unique_AD,]
write.table(AD_specific_df, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_Specific_pathways.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(truly_unique_AD_df, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_Specific_pathways_wo_DIAGRAM.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# 106/249 modules that are "AD" specific are also replicated in T2D


# get all AD genes and all T2D genes

results_files <- list.files("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/all_genes/")
results <- results_files[!grepl("Wojcik", results_files)]
results <- results[!grepl("DIAGRAMstage1_T2D.Kidney_Cortex_eQTL.50.20.details.txt|DIAGRAMstage1_T2D.Kidney_Cortex_sQTL.50.20.details.txt|Transethnic_AD.Kidney_Cortex_eQTL.50.20.details.txt|Transethnic_AD.Kidney_Cortex_sQTL.50.20.details.txt|UKB_AD.Brain_Substantia_nigra_sQTL.50.20.details.txt|UKB_AD.Kidney_Cortex_eQTL.50.20.details.txt|UKB_T2D.Brain_Amygdala_sQTL.50.20.details.txt|UKB_T2D.Kidney_Cortex_eQTL.50.20.details.txt|UKB_T2D.Kidney_Cortex_sQTL.50.20.details.txt", results)]
setwd("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/all_genes/")
truly_all_genes = data.frame(stringsAsFactors = FALSE)
for(file in results){
  temp <- makeNewGeneSetFileOnlySNPMappedGenesMod(geneset_file = "~/Desktop/Yang_Lab/Resources/genesets/all.mod.txt",
                                                  details_file = file)
  if(is.null(temp)) next
  truly_all_genes = rbind(truly_all_genes, temp)
}
write.table(all_genes, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/All_Genes.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(all_genes, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/All_GenesVALUE.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(all_genes, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/All_GenesVALUE.MARKER.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(truly_all_genes, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/Truly_All_Genes.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
all_genes_novalue <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/All_Genes.txt")
all_genes <-read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/All_GenesVALUE.txt", stringsAsFactors = FALSE)
all_genes <-read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/All_GenesVALUE.MARKER.txt", stringsAsFactors = FALSE)

makeNewGeneSetFileOnlySNPMappedGenesMod <- function(geneset_file, details_file, return_only_mod_genes=TRUE){
  details <- read.delim(details_file, header=TRUE, stringsAsFactors = FALSE)
  geneset <- read.delim(geneset_file, header=TRUE, stringsAsFactors = FALSE)
  details$FDR = gsub("%","", details$FDR)
  details$FDR <- as.numeric(details$FDR)
  details = details[details$FDR<10,] # changed just to get all the genes
  details = details[!(grepl("_ctrl", details$MODULE)),]
  if(nrow(details)==0) return(NULL)
  result_all = data.frame(stringsAsFactors = FALSE)
  result = data.frame(stringsAsFactors = FALSE)
  for(mod in unique(details$MODULE)){
    details %>% 
      filter(MODULE==mod) %>% 
      pull(GENE) -> genes
    details %>% 
      filter(MODULE==mod) %>% 
      pull(VALUE) -> values
    details %>% 
      filter(MODULE==mod) %>% 
      pull(MARKER) -> markers

    multiple_genes <- genes[grep(",", genes)]
    multiple_values = values[grep(",", genes)]
    multiple_markers = markers[grep(",", genes)]
    if(length(multiple_genes)!=0){
      split_genes <- lapply(multiple_genes, 
                            FUN = function(x){return(unlist(strsplit(x,split = ",")))})
      for(iter in 1:length(split_genes)){
        take_out <- which(!(split_genes[[iter]] %in% geneset$GENE[geneset$MODULE==mod]))
        if(length(take_out)==0) next
        split_genes[[iter]] = split_genes[[iter]][-take_out]
      }
      split_values = c()
      split_markers = c()
      for(iter in 1:length(multiple_values)){
        split_values = c(split_values, rep(multiple_values[iter],times=length(split_genes[[iter]])))
        split_markers = c(split_markers, rep(multiple_markers[iter],times=length(split_genes[[iter]])))
      }
      split_genes <- unlist(split_genes)
      trimmed_genes <- genes[!grepl(",", genes)]
      trimmed_values <- values[!grepl(",", genes)]
      trimmed_markers <- markers[!grepl(",", genes)]
      genes <- c(trimmed_genes, split_genes)
      values <- c(trimmed_values, split_values)
      markers <- c(trimmed_markers, split_markers)
    }
    if(!return_only_mod_genes){
      result_all = rbind(result_all, 
                         data.frame("MODULE"=mod, 
                                    "GENE"=genes, stringsAsFactors = FALSE))
    }
    if(return_only_mod_genes){
      result = rbind(result, 
                     data.frame("MODULE"=mod, 
                                "GENE"=genes, 
                                "VALUE"=values, 
                                "MARKER"=markers,
                                stringsAsFactors = FALSE))
    }
  }
  if(!return_only_mod_genes){
    write.table(result_all,"New_GeneSets_SNP_mapped_all.txt", sep="\t", quote = FALSE, row.names = FALSE)
  }
  else{
    result$Study = unlist(strsplit(details_file, split = ".", fixed = TRUE))[1]
    result$SGM = unlist(strsplit(details_file, split = ".", fixed = TRUE))[2]
    return(result)
    #write.table(result,"New_GeneSets_SNP_mapped_t2.txt", sep="\t", quote = FALSE, row.names = FALSE)
  }
}

makeNewGeneSetFileOnlySNPMappedGenes(geneset_file = "~/Desktop/Yang_Lab/Resources/genesets/all.mod.txt",
                                     details_file = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/all_genes/DIAGRAMstage1_T2D.Adipose_Subcutaneous_eQTL.50.20.details.txt")



# nrows all = 373916
# nrwos AD
# nrows T2D

# replicated T2D genes
study_genes = c()
for(study in unique(all_genes$Study)){
  study_genes <- c(study_genes, unique(all_genes$GENE[all_genes$Study==study & all_genes$VALUE>4]))
}
names(table(study_genes))[table(study_genes)==4]


tissues <- list.files("~/Desktop/GTEx/")
tissues = gsub(".txt","", tissues)
tissues = gsub("_"," ", tissues)
tissues = paste0("GTEx ", tissues)
tissues = paste0(tissues, " sQTL")
write.table(data.frame("Name"=c(tissues, "GTEx 49 Combined sQTLs")),"~/Downloads/GTEx_mapping_names_sQTL.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


files <- list.files("./")
for(file in files){
  df <- read.delim(file, stringsAsFactors = FALSE)
  colnames(df) <- c("GENE","MARKER")
  write.table(df, file=file, row.names = FALSE, quote = FALSE, sep = "\t")
}

setwd("/home/www/abhatta3-webserver/Data/Pipeline/Resources/GTEx_v8_eQTL")

setwd("/home/www/abhatta3-webserver/Data/Pipeline/Resources/GTEx_v8_sQTL")

# 108 mapping methods
SGMs <- list.files("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/all_genes/")
# SGMs <- SGMs[grep("Transethnic", SGMs)]
# SGMs <- gsub("Transethnic_AD.","",SGMs)

SGMs <- SGMs[grep("DIAGRAM", SGMs)]
SGMs <- gsub("DIAGRAMstage1_T2D.","",SGMs)

SGMs <- gsub(".50.20.details.txt","",SGMs)
SGMs <- gsub(".details.txt","",SGMs)

eQTLs <- SGMs[grep("eQTL", SGMs)]
eQTLs = eQTLs[-c(26,13,40,43)]
sQTLs <- SGMs[grep("sQTL", SGMs)]
sQTLs = sQTLs[-c(26,19,40,43)]

beginning <- SGMs[c(53, 54, 51, 52, 81, 82, 25, 38)]
SGMs <- c(beginning, eQTLs, sQTLs, SGMs[87], SGMs[88])
names(SGMs) <- 1:length(SGMs)

gsub("_"," ",concatenate(paste0(names(SGMs),": ", SGMs), mysep = ", "))

AD_specific_df$Mapping_full <- AD_specific_df$Mapping

for(sgm in names(SGMs)){
  AD_specific_df$Mapping <- gsub(SGMs[[sgm]],sgm,AD_specific_df$Mapping)
}
AD_specific_df$Mapping <- sapply(AD_specific_df$Mapping, 
                                 function(x){return(concatenate(unlist(strsplit(x, ", "))[order(as.numeric(unlist(strsplit(x, ", "))))],
                                                                mysep = ", "))})
AD_specific_df$Mapping <- gsub(", 107", "", AD_specific_df$Mapping)
AD_specific_df$Mapping <- gsub(", 108", "", AD_specific_df$Mapping)

write.table(AD_specific_df, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_Specific_pathways_thesis.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")
# merge modules
merge_modules(name = "AD_specific", 
              output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/merged/", 
              rcutoff = 0.3,
              modules_df = data.frame("MODULE"=AD_specific_df[,"MODULE"], stringsAsFactors = FALSE), 
              modfile_path = modfile_path, infofile_path = infofile_path)
# get merged modules
merged <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/merged/merged_AD_specific.mod.txt")
merged_mods <- unique(merged[grep(",", merged$OVERLAP),"OVERLAP"])
merged_mods[grep("rctm0917", merged_mods)]

merged_coexp <- unique(merged$OVERLAP[grep("GTEXv7", merged$OVERLAP)])
pathwayEnrichment(modules_list = data.frame("MODULE"=merged$MODULE[grep("GTEXv7",merged$OVERLAP)], "GENE"=merged$GENE[grep("GTEXv7",merged$OVERLAP)], stringsAsFactors = FALSE), 
                  resources_path = "~/Desktop/10X_Analysis/Resources/", 
                  output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_Specific_annotate/")
  
# Writing thesis -------
AD_specific_df <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_Specific_pathways.txt", 
                             stringsAsFactors = FALSE)
# how many coexpression and canonical pathways
AD_pathways <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Supplement/AD_all_pathways.txt")

# out of...
all_DESCR <- read.delim(infofile_path, stringsAsFactors = FALSE)
# 9315 coexpression modules, 1891 canonical pathways
length(grep("GTEXv7", all_DESCR$MODULE))  
11206-9315
length(grep("GTEXv7", AD_pathways$MODULE))  


length(grep("GTEXv7", AD_specific_df$MODULE)) 
nrow(AD_specific_df) - length(grep("GTEXv7", AD_specific_df$MODULE)) 

# making table 
all_genes %>% 
  filter(Study=="Transethnic_AD"|Study=="UKB_AD", MODULE=="Nerve_Tibial.GTEXv7.MEGENA_3", VALUE>3) %>% 
  pull(GENE) -> ex_genes
table(ex_genes)
merged_mods[grep("rctm0857", merged_mods)]

all_genes %>% 
  filter(Study=="DIAGRAMstage1_T2D", MODULE=="rctm0857", VALUE>3) %>% 
  pull(GENE) -> ex_genes
table(ex_genes)


# for table describing coexpression modules
all_mods <- read.delim("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/All_mods.mod.txt", stringsAsFactors = FALSE)
all_mods$Tissue <- sapply(all_mods$OVERLAP, function(x){return(unlist(strsplit(x, split = ".", fixed = TRUE))[1])})
all_mods$Merged <- sapply(all_mods$OVERLAP, function(x){return(ifelse(grep(",", x), "merged",""))})
unique(all_mods$OVERLAP)
new_df <- data.frame("MODULE"= unique(all_mods$OVERLAP), stringsAsFactors = FALSE)
new_df$Tissue <- sapply(new_df$MODULE, function(x){return(unlist(strsplit(x, split = ".", fixed = TRUE))[1])})
new_df$Merged <- lapply(new_df$MODULE, function(x){return(ifelse(grep(",", x), "merged","no"))})
new_df$tissue_Merged <- paste(new_df$Tissue, new_df$Merged, sep = "_")
new_df = new_df[grep("GTEXv7", new_df$MODULE),]
table(new_df$tissue_Merged)

all_coexp <- c()
for(row in 1:nrow(new_df)){
  all_coexp <- c(all_coexp, unlist(strsplit(new_df$MODULE[row], split = ",")))
}

AD_specific_df <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_Specific_pathways.txt", 
                             stringsAsFactors = FALSE)
AD_specific_coexp <- AD_specific_df[grep("GTEXv7", AD_specific_df$MODULE),]
AD_specific_coexp$MODULE <- gsub(".MEGENA","M", AD_specific_coexp$MODULE)
AD_specific_coexp$MODULE <- gsub(".WGCNA","W", AD_specific_coexp$MODULE)
AD_specific_coexp$MODULE <- gsub("GTEXv7","", AD_specific_coexp$MODULE)
write.table(AD_specific_coexp, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_Specific_coexp_pathways.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# specific T2D pathways ---------
AD_pathways <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Supplement/AD_all_pathways.txt")
T2D_pathways <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Supplement/T2D_all_pathways.txt")

length(grep("GTEXv7",T2D_pathways$MODULE)) # and 3 zhang, makes 116 significant coexpression modules

T2D_specific = setdiff(unique(T2D_pathways$MODULE), unique(AD_pathways$MODULE))
T2D_specific_d = setdiff(unique(T2D_pathways$DESCR), unique(AD_pathways$DESCR))

# merge AD specific
T2D_df = makeDf(Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/T2D.2/") # takes shared files from the quasi overlap

all_DESCR <- read.delim(infofile_path, stringsAsFactors = FALSE)

# very sig is a value cut off of 0.01. not very high
# make list of AD-specific pathways and mappings, same for T2D - eQTL specific!
T2D_specific_df = describePathways(modules_list = T2D_specific, df = T2D_df, DESCR = all_DESCR)
T2D_specific_df <- addGenes(df = T2D_specific_df, genes_files = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/DIAGRAM_UKB/very_sig_genes/")
genes_df_T2D <- genes_df[genes_df$Study=="DIAGRAMstage1_T2D"| 
                         genes_df$Study=="UKB_T2D",]
T2D_specific_df_genes2 <- addGenesv2(df = T2D_specific_df, 
                                    genes_df = genes_df_T2D, 
                                    value_threshold = 2, 
                                    replicated_studies = c("DIAGRAMstage1_T2D",
                                                           "UKB_T2D"))
T2D_specific_df_genes3 <- addGenesv2(df = T2D_specific_df, 
                                     genes_df = genes_df_T2D, 
                                     value_threshold = 3, 
                                     replicated_studies = c("DIAGRAMstage1_T2D",
                                                            "UKB_T2D"))
T2D_specific_df_genes3 <- T2D_specific_df_genes3[T2D_specific_df_genes3$Mapping!="select_24ssnps_sQTL",]
T2D_specific_df_genes3 <- T2D_specific_df_genes3[T2D_specific_df_genes3$Mapping!="select_24esnps_eQTL",]
# count canonical and coexpression (not merged) after taking out the select_esnps
length(T2D_specific_df_genes3$MODULE[grep("GTEXv7",T2D_specific_df_genes3$MODULE)])
# 43 + 1 coexpression 
nrow(T2D_specific_df_genes3) - length(T2D_specific_df_genes3$MODULE[grep("GTEXv7",T2D_specific_df_genes3$MODULE)]) # 201-2 is 199

merge_modules(name = "T2D_specific", 
              output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/merged/", 
              rcutoff = 0.3,
              modules_df = data.frame("MODULE"=T2D_specific_df_genes3[,"MODULE"], stringsAsFactors = FALSE), 
              modfile_path = modfile_path, infofile_path = infofile_path)

# get merged modules 
merged <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/merged/merged_T2D_specific.mod.txt", stringsAsFactors = FALSE)
T2D_specific_df_genes3 <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/T2D_Specific_pathways_3VALUE2Replicated.txt", stringsAsFactors = FALSE)
T2D_specific_df_genes3 <- addSuperSet(df = T2D_specific_df_genes3, merged_df = merged)

write.table(T2D_specific_df_genes3, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/T2D_Specific_pathways_3VALUE2Replicated_noSelect.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")


# shared AD T2D ------
AD_T2D_df = makeDf(Dir="~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D_new/")
AD_T2D_mod_df = describePathways(modules_list = unique(AD_T2D_df$MODULE), 
                                 df = AD_T2D_df, DESCR = all_DESCR)
genes_df <- all_genes
AD_T2D_mod_df <- addGenesv2(df=AD_T2D_mod_df, 
                            genes_df = genes_df, 
                            value_threshold = 3, 
                            replicated_studies = c("DIAGRAMstage1_T2D",
                                                   "Transethnic_AD",
                                                   "UKB_AD","UKB_T2D"))
write.table(AD_T2D_mod_df,"~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_T2D_pathways_4replicated3value_genes.txt",row.names = FALSE, quote = FALSE, sep="\t")
write.table(AD_T2D_mod_df,"~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_T2D_pathways_3replicated3value_genes.txt",row.names = FALSE, quote = FALSE, sep="\t")
write.table(AD_T2D_mod_df,"~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_T2D_pathways_4replicated2value_genes.txt",row.names = FALSE, quote = FALSE, sep="\t")
write.table(AD_T2D_mod_df,"~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_T2D_pathways_4replicated0value_genes.txt",row.names = FALSE, quote = FALSE, sep="\t")

# merge modules
# take out just select esnps
AD_T2D_mod_df <- AD_T2D_mod_df[AD_T2D_mod_df$Mapping!="select_24ssnps_sQTL",]
AD_T2D_mod_df <- AD_T2D_mod_df[AD_T2D_mod_df$Mapping!="select_24ssnps_eQTL",]
# count canonical and coexpression (not merged) after taking out the select_esnps
length(AD_T2D_mod_df$MODULE[grep("GTEXv7",AD_T2D_mod_df$MODULE)])
# 75 + 2 (light cyan and gold) is 77
nrow(AD_T2D_mod_df) - length(AD_T2D_mod_df$MODULE[grep("GTEXv7",AD_T2D_mod_df$MODULE)]) # 201-2 is 199

merge_modules(name = "AD_T2D", 
              output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/merged/", 
              rcutoff = 0.3,
              modules_df = data.frame("MODULE"=AD_T2D_mod_df[,"MODULE"], stringsAsFactors = FALSE), 
              modfile_path = modfile_path, infofile_path = infofile_path)
# get merged modules 
merged <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/merged/merged_AD_T2D.mod.txt", stringsAsFactors = FALSE)
merged_mods <- unique(merged[grep(",", merged$OVERLAP),"OVERLAP"])
merged_mods[grep("rctm0917", merged_mods)]

# add superset information 
library(stringr)
AD_T2D_mod_df <- addSuperSet(df = AD_T2D_mod_df, merged_df = merged)

# convert mapping methods to numbers
AD_T2D_mod_df$Mapping_full <- AD_T2D_mod_df$Mapping

for(sgm in names(SGMs)){
  AD_T2D_mod_df$Mapping <- gsub(SGMs[[sgm]],sgm,AD_T2D_mod_df$Mapping)
}
AD_T2D_mod_df$Mapping <- sapply(AD_T2D_mod_df$Mapping, 
                                function(x){return(concatenate(unlist(strsplit(x, ", "))[order(as.numeric(unlist(strsplit(x, ", "))))],
                                                               mysep = ", "))})
AD_T2D_mod_df$Mapping <- gsub(", 107", "", AD_T2D_mod_df$Mapping)
AD_T2D_mod_df$Mapping <- gsub(", 108", "", AD_T2D_mod_df$Mapping)


# order coexpression last - 276 total
coexp = AD_T2D_mod_df[grep("GTEXv7|Gold|Light cyan",AD_T2D_mod_df$MODULE),]
AD_T2D_mod_df = AD_T2D_mod_df[!grepl("GTEXv7|Gold|Light cyan",AD_T2D_mod_df$MODULE),]
AD_T2D_mod_df = rbind(AD_T2D_mod_df, coexp)

write.table(AD_T2D_mod_df, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_T2D_pathways_genes4replicated3value.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

merged_coexp <- unique(merged$OVERLAP[grep("GTEXv7", merged$OVERLAP)])
pathwayEnrichment(modules_list = data.frame("MODULE"=merged$MODULE[grep("GTEXv7",merged$OVERLAP)], "GENE"=merged$GENE[grep("GTEXv7",merged$OVERLAP)], stringsAsFactors = FALSE), 
                  resources_path = "~/Desktop/10X_Analysis/Resources/", 
                  output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/AD_Specific_annotate/")
length(AD_T2D_mod_df$MODULE[grep("GTEXv7", AD_T2D_mod_df$MODULE)])

# i didn't run zhang modules with combined brain eQTLs??
# also concerned that Yellow module got filtered out because it is over 1000 genes




# redo numbering 

# new names for SGMs
# 3 is combined
# 4 is peripheral
# 5 is brain
new_SGM_names <- c(1:2, "3e","3s","4e","4s","5e","5s",paste0(6:54,"e"), paste0(6:54,"s"))
old_names_new_names = new_SGM_names
names(old_names_new_names) = 1:(length(SGMs)-2)
SGMs = SGMs[-c(107,108)]
names(SGMs) = new_SGM_names
info = read.delim(infofile_path, stringsAsFactors = FALSE)

changeSNPcodeAddSource <- function(df, info){
  df$Mapping <- sapply(df$Mapping, 
                       function(x){return(concatenate(sapply(unlist(strsplit(x, split = ", ")),
                                                      function(x){return(old_names_new_names[x])}), 
                                                      mysep = ", "))})
  df$Pathway <- sapply(df$Pathway, function(x){return(info$SOURCE[info$MODULE==x])})
  df$Pathway <- gsub("reactome","Reactome", df$Pathway)
  df$Pathway <- gsub("kegg","KEGG", df$Pathway)
  df$Pathway <- gsub("biocarta","BioCarta", df$Pathway)
  #df$Ex..Suggestive.genes <- gsub("_","",df$Ex..Suggestive.genes)
  return(df)
}
trimmed = read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/Trimmed_AD_Specific_Table.txt", 
                     stringsAsFactors = FALSE)
write.table(changeSNPcodeAddSource(df = trimmed, info = info),
            "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/Trimmed_AD_Specific_Table.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)
trimmed = read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/Trimmed_T2D_Specific_Table.txt", 
                     stringsAsFactors = FALSE)
write.table(changeSNPcodeAddSource(df = trimmed, info = info),
            "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/Trimmed_T2D_Specific_Table.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)
trimmed = read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/Trimmed_AD_T2D_Table.txt", 
                     stringsAsFactors = FALSE)
write.table(changeSNPcodeAddSource(df = trimmed, info = info),
            "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/Trimmed_AD_T2D_Table.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)     
            
concatenate(gsub("_"," ",paste0(3:54, ": ",gsub("_eQTL","", SGMs[c(3,5,7,9:57)]))),mysep = ", ")
            
            
            
# make a node file for combined eQTLs
shared_AD_T2D <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/shared_AD_T2D_quasi.combined_49esnps_eQTL.txt", 
                            stringsAsFactors = FALSE)
all <- read.delim(modfile_path, stringsAsFactors = FALSE)
nodes <- data.frame(stringsAsFactors = FALSE)
for(mod in shared_AD_T2D$MODULE){
  nodes <- rbind(nodes, 
                 data.frame("MODULE"=mod,
                            "NODE"=all$GENE[all$MODULE==mod], 
                            stringsAsFactors = FALSE))
}
write.table(nodes, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Combined_eQTLs.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")          
merge_modules(name = "AD_T2D_Combined_eQTLs", 
              output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/", 
              rcutoff = 0.3,
              modules_df = shared_AD_T2D, 
              modfile_path = modfile_path, infofile_path = infofile_path)

mod<-read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/merged_AD_T2D_Combined_eQTLs.mod.txt")
write.table(mod[,c("MODULE","NODE")], "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Combined_eQTLs_merged.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")  


shared_AD_T2D <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/shared_AD_T2D_quasi.peripheral_esnps_eQTL.txt", 
                            stringsAsFactors = FALSE)
all <- read.delim(modfile_path, stringsAsFactors = FALSE)
nodes <- data.frame(stringsAsFactors = FALSE)
for(mod in shared_AD_T2D$MODULE){
  nodes <- rbind(nodes, 
                 data.frame("MODULE"=mod,
                            "NODE"=all$GENE[all$MODULE==mod], 
                            stringsAsFactors = FALSE))
}
write.table(nodes, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Combined_eQTLs.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")          
merge_modules(name = "AD_T2D_Peripheral_eQTLs", 
              output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/", 
              rcutoff = 0.3,
              modules_df = shared_AD_T2D, 
              modfile_path = modfile_path, infofile_path = infofile_path)

mod<-read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/merged_AD_T2D_Peripheral_eQTLs.mod.txt")
write.table(mod[,c("MODULE","NODE")], "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Peripheral_eQTLs_merged.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")  

files <- list.files("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/")
files <- files[intersect(grep("rain_",files), grep("_eQTL",files))]
files <- files[grep("shared_",files)]
shared_AD_T2D <- data.frame(stringsAsFactors = FALSE)
for(file in files){
shared_AD_T2D <- rbind(shared_AD_T2D, read.delim(paste0("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD_T2D/merged/",file)))
}
all <- read.delim(modfile_path, stringsAsFactors = FALSE)
nodes <- data.frame(stringsAsFactors = FALSE)
for(mod in unique(shared_AD_T2D$MODULE)){
  nodes <- rbind(nodes, 
                 data.frame("MODULE"=mod,
                            "NODE"=all$GENE[all$MODULE==mod], 
                            stringsAsFactors = FALSE))
}
write.table(nodes, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Combined_eQTLs.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")          
merge_modules(name = "AD_T2D_Brain_eQTLs", 
              output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/", 
              rcutoff = 0.3,
              modules_df = shared_AD_T2D, 
              modfile_path = modfile_path, infofile_path = infofile_path)

mod<-read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/merged_AD_T2D_Brain_eQTLs.mod.txt")
write.table(mod[,c("MODULE","NODE")], "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Brain_eQTLs_merged.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")  


# quantify input nodes
mod<-read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/merged_AD_T2D_Combined_eQTLs.mod.txt")
length(unique(mod$MODULE))
length(unique(mod$GENE))
net <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Combined_eQTLs_merged_networks.hs.all/cytoscape/kda2cytoscape.edges.txt")
nodes <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Combined_eQTLs_merged_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt")
concatenate(unique(net$TAIL), mysep = ", ")
length(unique(net$TAIL))
length(unique(c(net$TAIL, net$HEAD)))

mod<-read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/merged_AD_T2D_Peripheral_eQTLs.mod.txt") 
mod<-read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/merged_AD_T2D_Brain_eQTLs.mod.txt")

net <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Brain_eQTLs_merged_networks.hs.brain/cytoscape/kda2cytoscape.edges.txt")
nodes <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Brain_eQTLs_merged_networks.hs.brain/cytoscape/kda2cytoscape.nodes.txt")
concatenate(unique(net$TAIL), mysep = ", ")

sum(nodes$URL!="")

combined_KDs <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Combined_eQTLs_merged_networks.hs.all/kda/wKDA.results.txt")
brain_KDs <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Brain_eQTLs_merged_networks.hs.brain/kda/wKDA.results.txt")
concatenate(intersect(combined_KDs$NODE[combined_KDs$FDR<0.05], brain_KDs$NODE[brain_KDs$FDR<0.05]), mysep = ", ")

periph_KDs <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Peripheral_eQTLs_merged_networks.hs.all/kda/wKDA.results.txt")

concatenate(intersect(periph_KDs$NODE[periph_KDs$FDR<0.05], brain_KDs$NODE[brain_KDs$FDR<0.05]), mysep = ", ")

mod <- data.frame("MODULE"="mod","GENE"= intersect(periph_KDs$NODE[periph_KDs$FDR<0.05], brain_KDs$NODE[brain_KDs$FDR<0.05]))
mod <- data.frame("MODULE"="mod","GENE"= setdiff(periph_KDs$NODE[periph_KDs$FDR<0.05], brain_KDs$NODE[brain_KDs$FDR<0.05]))
mod <- data.frame("MODULE"="mod","GENE"= setdiff(brain_KDs$NODE[brain_KDs$FDR<0.05], periph_KDs$NODE[periph_KDs$FDR<0.05]))
periph_KDs_mod <- data.frame("MODULE"="Kds","GENE"=periph_KDs$NODE[periph_KDs$FDR<0.001], stringsAsFactors = FALSE)
pathway_dir = "./forThesisBrainOnly"
pathway_resource = "~/Desktop/10X_Analysis/Resources/"
pathwayEnrichment(modules_list = mod, 
                  resources_path = pathway_resource, output_Dir = pathway_dir)
pathway_files = list.files(pathway_dir)[grep(".txt", list.files(pathway_dir))]
pathway_df = data.frame(stringsAsFactors = FALSE)
for(file in pathway_files){
  temp = read.delim(paste0(pathway_dir,"/",file), stringsAsFactors = FALSE, header = TRUE)
  pathway_df = rbind(pathway_df, temp)
}
pathway_df = pathway_df[order(pathway_df$nOverlap, decreasing = TRUE),]

HYP <- readRDS("~/Downloads/HYP_Neurons_WT_5XFAD_Russ.rds")


net <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Peripheral_eQTLs_merged_networks.hs.all/cytoscape/kda2cytoscape.edges.txt")
nodes <- read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Peripheral_eQTLs_merged_networks.hs.all/cytoscape/kda2cytoscape.nodes.txt")
concatenate(unique(net$TAIL), mysep = ", ")

white_adipo <- all_genes[all_genes$MODULE=="rctm0858",]
white_adipo <- white_adipo[white_adipo$SGM=="Colon_Sigmoid_eQTL",]

genes_study <- all_genes[,c("GENE","Study","VALUE")]
genes_study <- genes_study[genes_study$VALUE>3,]
genes_study <- genes_study[!duplicated(genes_study),]
setdiff(unique(c(genes_study$GENE[genes_study$Study=="DIAGRAMstage1_T2D"], genes_study$GENE[genes_study$Study=="UKB_T2D"])),
        unique(c(genes_study$GENE[genes_study$Study=="Transethnic_AD"], genes_study$GENE[genes_study$Study=="UKB_AD"])))

setdiff(unique(c(genes_study$GENE[genes_study$Study=="Transethnic_AD"], genes_study$GENE[genes_study$Study=="UKB_AD"])),
        unique(c(genes_study$GENE[genes_study$Study=="DIAGRAMstage1_T2D"], genes_study$GENE[genes_study$Study=="UKB_T2D"])))

# genes that are both in DIAGRAMstage1 and Transethnic
genes_study <- genes_study[,c("GENE","Study")]
genes_study <- genes_study[!duplicated(genes_study),]
T2D_replicated = genes_study[genes_study$Study=="DIAGRAMstage1_T2D" | genes_study$Study=="UKB_T2D",]
T2D_replicated = T2D_replicated$GENE[duplicated(T2D_replicated$GENE)]

mod <- T2D_replicated

AD_replicated = genes_study[genes_study$Study=="Transethnic_AD" | genes_study$Study=="UKB_AD",]
AD_replicated = AD_replicated$GENE[duplicated(AD_replicated$GENE)]
length(intersect(T2D_replicated,AD_replicated))
length(setdiff(T2D_replicated, AD_replicated))
length(setdiff(AD_replicated,T2D_replicated))

marker_values <- all_genes[,c("MARKER","VALUE","Study")]
marker_values <- marker_values[!duplicated(marker_values),]
quantile(marker_values$VALUE, probs = c(0.05, 0.95))
sum(marker_values$VALUE<=3)/nrow(marker_values)
hist(marker_values$VALUE[marker_values$VALUE<3.5], xlab = "-log10(p-value)", 
     main = "Histogram of 96% of Marker Association Values")

all_genes$MODULE_Study_SGM <- paste(all_genes$MODULE, all_genes$Study, all_genes$SGM, sep = "_")

all_genes %>% 
  group_by(MODULE_Study_SGM) %>% 
  mutate(Avg_value=mean(VALUE)) -> all_genes_new
hist(all_genes_new$Avg_value)
hist(all_genes_new$Avg_value[all_genes_new$Avg_value<4], xlab = "Average -log10(p-value)",
     main = "Histogram of Gene Module Marker Association Value Averages")
quantile(all_genes_new$Avg_value, probs = c(0.05, 0.95))
sum(all_genes_new$Avg_value<=4)/nrow(all_genes_new)

all_genes_new[all_genes_new$Avg_value>50,]


# get all AD modules
modfile_path = "/Users/jessicading/Desktop/Yang_Lab/resources/genesets/all.mod.txt"
files <- list.files("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/")
files <- files[intersect(grep("rain_",files), grep("_eQTL",files))]
files <- files[grep("shared_AD",files)]
shared_AD_T2D <- data.frame(stringsAsFactors = FALSE)
for(file in files){
  shared_AD_T2D <- rbind(shared_AD_T2D, read.delim(paste0("/Users/jessicading/Desktop/Yang_Lab/T2D_AD/all_results_v8/AD/overlap/",file)))
}
all <- read.delim(modfile_path, stringsAsFactors = FALSE)
nodes <- data.frame(stringsAsFactors = FALSE)
for(mod in unique(shared_AD_T2D$MODULE)){
  nodes <- rbind(nodes, 
                 data.frame("MODULE"=mod,
                            "NODE"=all$GENE[all$MODULE==mod], 
                            stringsAsFactors = FALSE))
}


write.table(nodes, "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/networks/AD_Brain_Nodes.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")          
merge_modules(name = "AD_T2D_Brain_eQTLs", 
              output_Dir = "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/", 
              rcutoff = 0.3,
              modules_df = shared_AD_T2D, 
              modfile_path = modfile_path, infofile_path = infofile_path)

mod<-read.delim("~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/merged_AD_T2D_Brain_eQTLs.mod.txt")
write.table(mod[,c("MODULE","NODE")], "~/Desktop/Yang_Lab/T2D_AD/all_results_v8/Thesis/KDA/AD_T2D_Brain_eQTLs_merged.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")  


# testing Hengjian's mod file (brain, peripheral, shared network genes)
library(stringr)
df <- read.delim("~/Downloads/PathwayEnrichment.txt", stringsAsFactors = FALSE)
df <- df[duplicated(df$GENE),]
df$MODULE <- gsub("perioheral specific","peripheral specific", df$MODULE)
microarray_genes <- df$GENE[str_count(df$GENE, "[0-9]")==6 & str_count(df$GENE,"[A-Z]")==2]
microarray_gene_alpha = unique(gsub("[[:digit:]]","",microarray_genes))
df <- df[!(df$GENE %in% microarray_genes),]
df <- df[!grepl("^MMT|XM_|NM_|ENSMUS|_at|^MMS|^DKF|^mCT|_", df$GENE),]
#df <- df[!grepl("^MMT|XM_|NM_|ENSMUS|_at|^MMS|^AK0|^BC0|^mCT|^BB4|^AW|^DKF|^AY|^BE[[:digit:]]|^AI[[:digit:]]|_", df$GENE),]
table(df$MODULE) # peripheral still has over 20000
# change total genes in makePathwayEnrichmentDf() 


df <- df[1:20000,] # works
df <- df[1:20001,] # doesn't work
df <- df[20000:25000,] # works
# issue of how many genes are being enriched
pathway_df <- makePathwayEnrichmentDf(DEG_df = df, output_Dir = "test_Jane",resources_path = "../Resources/", convertToHuman = FALSE)

pathway_df <- makePathwayEnrichmentDf(DEG_df = PE, output_Dir = "test_Jane",resources_path = "../Resources/", convertToHuman = FALSE)

