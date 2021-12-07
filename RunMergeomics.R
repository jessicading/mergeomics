source("Mergeomics_Version_1.99.0.R")
source("Mergeomics_utils.R")

# See Mergeomics_utils.R to see all possible parameter settings

# GWAS Enrichment
# runMDF must be run on Linux distribution
# alternatively, use the webserver to run MDF (Run Mergeomics > Individual GWAS Enrichment)
msea_job <- runMDF(marker_associations= "GWAS/GWAS_MSEA/Kunkle_AD.txt", # marker association file
                   marker_mapping = "mapping/MONOCYTES_EQTL.txt", # marker to gene mapping file
                   marker_dependency = "linkage/LD50.1000G.CEU.txt", # marker dependency file
                   n_top = 0.2, # top proportion of marker association file to consider
                   md_threshold = 50, # dependency (correlation/R2) cutoff, used to label output
                   mdprune = "mdprune" # path to mdprune program
)
msea_job <- runMSEA(job=msea_job, 
                    marker_set="genesets/KEGG_Reactome_BioCarta.txt")
kda_job <- runKDA(job=msea_job, 
                  merge_modules=TRUE,# whether to merge redundant modules, TRUE recommended
                  network="network/bayesian_network_brain.txt")

# TWAS/PWAS Enrichment
msea_job <- runMSEA(marker_associations="Data/Monocyte_DE_Genes.txt", 
                    marker_set="genesets/KEGG_Reactome_BioCarta.txt")
# KDA can be run as depicted above

# For EWAS enrichment, maxoverlap=1 and permtype="marker" are required.
msea_job <- runMSEA(marker_associations="Data/EWAS_GSE31835.txt", 
                    marker_mapping="Data/EMAP31835.txt",
                    marker_set="genesets/KEGG_Reactome_BioCarta.txt",
                    permtype="marker",
                    maxoverlap=1)
# KDA can be run as depicted above