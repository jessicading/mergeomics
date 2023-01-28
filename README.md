*** Please use Mergeomics_Version_1.99.0.R and report any issues

# Mergeomics
Integrative Network Analysis of Omics Data<br/>
Developed by: Ville-Petteri Makinen, Le Shu, Yuqi Zhao, Zeyneb Kurt, Jessica Ding, Bin Zhang, Xia Yang

#### Table of Contents
- [About](#about)<br/>
- [Tutorial](#tutorial)<br/>
    - [Example Workflow](#example-workflow)
    - [Marker Dependency Filtering (MDF)](#marker-dependency-filtering)
    - [Marker Set Enrichment Analysis (MSEA)](#marker-set-enrichment-analysis)
    - [Weighted Key Driver Analysis (KDA)](#weighted-key-driver-analysis)
    - [Module Merging](#module-merging)
- [Citations](#citations)

## About
Mergeomics is an open source R-based pipeline for multi-dimensional integration of omics data to identify disease-associated pathways and networks. Genes whose network neighborhoods are over-represented with disease associated genes are deemed key drivers and can be targeted in further mechanistic studies.

![abstract figure](https://github.com/jessicading/mergeomics/blob/master/Abstract.jpg?raw=true)

More information can be found in the [methods paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3198-9), [documentation](http://bioconductor.org/packages/release/bioc/html/Mergeomics.html), [web server](http://mergeomics.research.idre.ucla.edu), and [web server paper](https://academic.oup.com/nar/article/49/W1/W375/6287846).

## Tutorial
This tutorial illustrates the workflow for running the Mergeomics pipeline which includes Marker Dependency Filtering (MDF), Marker Set Enrichment Analysis (MSEA), Meta-MSEA, and Key Driver Analysis (KDA) using the R package. For those who prefer a web interface, please use our [web server](http://mergeomics.research.idre.ucla.edu).

A streamlined workflow using wrapper functions in Mergeomics_utils.R is outlined first, followed by descriptions of the parameters, inputs, and outputs for each step. 


#### Notes
- A new version of the Mergeomics R package will be released which fixes minor bugs and improves results output. We suggest using `source(Mergeomics_Version_1.99.0.R)` for the updated functions.
- All input files should be tab delimited .txt files.
- The term 'marker' is used to refer to any type of biological marker including loci, methylation/epigenetic sites, transcripts/genes, proteins, metabolites, etc.


### Example workflow 
Setting of all parameters is displayed in toggle below in "See setting of all possible parameters" section.
```R
source("Mergeomics_Version_1.99.0.R")
source("Mergeomics_utils.R")

# GWAS Enrichment
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

# For EWAS enrichment, maxoverlap=1 and permtype="marker" are required
msea_job <- runMSEA(marker_associations="Data/EWAS_GSE31835.txt", 
                    marker_mapping="Data/EMAP31835.txt",
                    marker_set="genesets/KEGG_Reactome_BioCarta.txt",
                    permtype="marker",
                    maxoverlap=1)
# KDA can be run as depicted above
```
#### Example Meta-MSEA workflow
```R
source("Mergeomics_Version_1.99.0.R")
source("Mergeomics_utils.R")

# Option 1: 
# Running MSEA on multiple of the same type or different types of omics datasets
# followed by Meta-MSEA using runMetaMSEA
job.meta <- runMetaMSEA(msea_input_list=list("gwas"=list(marker_associations="Data/Kunkle_AD.MONOCYTES_EQTL/top50.md50.m.txt",
                                                         marker_mapping="Data/Kunkle_AD.MONOCYTES_EQTL/top50.md50.g.txt"),
                                             "degs"=list(marker_associations="Data/Monocyte_DE_Genes.txt"),
                                             "deps"=list(marker_associations="Data/Monocyte_DE_Proteins.txt")),
                        marker_set = "genesets/KEGG_Reactome_BioCarta.txt")
# then run KDA
kda_job <- runKDA(job=job.meta, 
                  MSEA_fdr_cutoff = c(0.5,0.5,0.5),
                  merge_modules=TRUE, 
                  network="networks/bayesian.hs.blood.txt") 

# Option 2:
# Running Meta-MSEA on finished MSEA runs returned by runMSEA()
msea_job1 <- runMSEA(job=msea_job, 
                     marker_set="genesets/KEGG_Reactome_BioCarta.txt")
msea_job2 <- runMSEA(marker_associations="Data/Monocyte_DE_Genes.txt", 
                     marker_set="genesets/KEGG_Reactome_BioCarta.txt")
msea_job3 <- runMSEA(marker_associations="Data/Data/Monocyte_DE_Proteins.txt",
                     marker_set="genesets/KEGG_Reactome_BioCarta.txt")     
jobs <- list(msea_job1, msea_job2, msea_job3)
job.meta <- ssea.meta(jobs)

job.meta.kda <- list()
job.meta.kda$metamsea <- job.meta
job.meta.kda$modfile <- "genesets/KEGG_Reactome_BioCarta.txt"

# Next, run KDA with completed Meta-MSEA job
kda_job <- runKDA(job=job.meta.kda, 
                  MSEA_fdr_cutoff = c(0.5,0.5,0.5),
                  merge_modules=TRUE, 
                  network="networks/bayesian.hs.blood.txt") 

```
### See setting of all possible parameters
<details>
<summary>Click to see setting of all possible parameters from the functions used above</summary>
	
```R

# MSEA of GWAS followed by KDA

msea_job <- runMDF(marker_associations= "GWAS/GWAS_MSEA/Kunkle_AD.txt", # marker association file
                   marker_mapping = "mapping/MONOCYTES_EQTL.txt", # marker to gene mapping file
                   marker_dependency = "linkage/LD50.1000G.CEU.txt", # marker dependency file
                   NTOP = 0.2, # top proportion of marker association file to consider
                   ld_threshold = 50, # dependency (correlation/R2) cutoff, used to label output
                   mdprune = "mdprune" # path to mdprune program,
                   output_dir="Data/",# where to output results
                   edit_files=FALSE # edit headers of input files to appropriate headere
                   )
                   
msea_job <- runMSEA(job=msea_job, 
                    marker_set="genesets/KEGG_Reactome_BioCarta.txt",
                    marker_set_info="genesets/KEGG_Reactome_BioCarta_info.txt",
                    output_dir="Results", 
                    label="Kunkle_AD.MONOCYTES_EQTL", 
                    permtype="gene",
                    nperm=10000,
                    maxoverlap=.33,
                    max_module_genes=500,
                    min_module_genes=10,
                    trim=0.002,
                    seed=1,
                    return_job=TRUE)
# OR set marker_associations and marker_mapping (no marker_mapping needed for TWAS/PWAS)
msea_job <- runMSEA(marker_associations="Data/Kunkle_AD.MONOCYTES_EQTL/top50.md50.m.txt",
                    marker_mapping="Data/Kunkle_AD.MONOCYTES_EQTL/top50.md50.g.txt",
                    marker_set="genesets/KEGG_Reactome_BioCarta.txt")

# provide complete msea job to retrieve significantly enriched marker sets
kda_job <- runKDA(job=msea_job, 
                  MSEA_fdr_cutoff=0.05,
                  merge_modules=TRUE,
                  merge_rcutoff=.33,
                  network="network/bayesian_network_blood.txt",
                  label="Kunkle_AD.MONOCYTES_EQTL",
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
                  save_job = TRUE)  
                  
# OR provide the MSEA result file and the marker/gene sets file
kda_job <- runKDA(MSEA_results="Results/msea/Kunkle_AD.MONOCYTES_EQTL.results.txt", 
                  marker_set_file="genesets/KEGG_Reactome_BioCarta.txt",
                  network="network/bayesian_network_blood.txt")
                  
# OR provide a vector of marker/gene sets and the marker/gene sets file
kda_job <- runKDA(marker_sets=c("REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"), 
                  marker_set_file="genesets/KEGG_Reactome_BioCarta.txt",
                  network="network/bayesian_network_blood.txt")
                  
# OR provide a vector of genes (no grouping)
kbr <- read.delim("genesets/KEGG_Reactome_BioCarta.txt")
query_genes = kbr$GENE[kbr$MODULE=="REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"]
kda_job <- runKDA(nodes=query_genes, 
                  network="network/bayesian_network_blood.txt")
                  
# OR provide a 'MODULE' 'NODE' header .txt file
kda_job <- runKDA(nodes="Data/KDA_query_nodes.txt", 
                  network="network/bayesian_network_blood.txt")
                  
# by default, subnetworks of 1 depth for the top 5 key drivers of each module 
# are output for network visualization. kda2cytoscape can be rerun on kda_job with
# different depth, number of top key drivers (10, shown beelow), etc.
kda_job <- kda2cytoscape(kda_job, ndrivers = 10)

# Meta-MSEA Example 1
# one MSEA with all possible parameters set
# TWAS/PWAS enrichment does not need permtype and maxoverlap parameters
job.meta <- runMetaMSEA(msea_input_list=list("gwas"=list(marker_associations="Data/Kunkle_AD.MONOCYTES_EQTL/top50.md50.m.txt",
                                                         marker_mapping="Data/Kunkle_AD.MONOCYTES_EQTL/top50.md50.g.txt",
                                                         permtype="gene",
                                                         nperm=10000,
                                                         maxoverlap=0.33,
                                                         max_module_genes=500,
                                                         min_module_genes=10,
                                                         trim=0.002,
                                                         seed=1),
                                             "degs"=list(marker_associations="Data/Monocyte_DE_Genes.txt"),
                                             "deps"=list(marker_associations="Data/Monocyte_DE_Proteins.txt")),
                        marker_set = "genesets/KEGG_Reactome_BioCarta.txt")

# Meta-MSEA Example 2
# EWAS enrichment needs permtype="marker" and maxoverlap=1 parameters
job.meta <- runMetaMSEA(msea_input_list=list("gwas"=list(marker_associations="Data/GWAS_Psoriasis.txt",
                                                         marker_mapping="Data/eQTL_skinblood.txt",
                                                         nperm=50),
                                             "ewas1"=list(marker_associations="Data/EWAS_GSE31835.txt",
                                                          marker_mapping="Data/EMAP31835.txt",
                                                          permtype="marker",
                                                          maxoverlap=1),
                                             "ewas2"=list(marker_associations="Data/EWAS_GSE63315.txt",
                                                          marker_mapping="Data/EMAP63315.txt",
                                                          permtype="marker",
                                                          maxoverlap=1)),
                        marker_set = "genesets/KEGG_Reactome_BioCarta.txt")
```
	
</details>


MSEA and KDA analysis parameter descriptions are detailed below.


### Marker Dependency Filtering

Marker dependency filtering (MDF) removes dependent markers and prepares an optimized marker and gene file for marker set enrichment analysis (MSEA). Dependency between markers can be a result of linkage disequilibrium, for example. Download the mdprune executable in this repository and provide the path to the program in the script below. This is a C++ program and requires a Linux distribution. If you cannot obtain the environment, MDF can be run using our [web server](http://mergeomics.research.idre.ucla.edu) (Run Mergeomics > Individual GWAS Enrichment).

For epigenomics, transcriptomics, proteomics, or metabolomics data, MSEA can be run as a first step (it is recommended to use the `runMSEA` function to simplify the analysis). MDF can be skipped for GWAS as well but it is not recommended.

#### Inputs
1. ```MARFILE```: Disease/Phenotype marker association data (summary statistics) (```marker_associations```)<br/> 
This file must have two columns named 'MARKER' and 'VALUE' where value denotes the strength of the association to the trait. This is usually -log10 p-values but can be other values such as effect size or fold change. Higher values must indicate higher association. The entire marker association file, including nonsignificant associations, should be used.
```
MARKER             VALUE
rs4747841         0.1452
rs4749917         0.1108
rs737656          1.3979
```
2. ```GENFILE```: Mapping file to connect markers from association data to genes (```marker_mapping```)<br/>
(eQTLs, distance based, etc.)
```
GENE              MARKER
TSPAN6            rs1204389
TNMD              rs113630520
SCYL3             rs72691775
```
3. ```LNKFILE```: Marker Dependency File (ex. Linkage Disequilibirum) (```marker_dependency```)<br/>
LD files can be obtained from the [HapMap consortium](https://www.sanger.ac.uk/resources/downloads/human/hapmap3.html). 
```
MARKERa            MARKERb              WEIGHT            
rs143225517       rs3094315           0.962305
rs4475691         rs950122            1
rs4475691         rs3905286           0.921467
```
4. ```OUTPATH```: Output Directory (```output_dir```)
5. ```NTOP```: Top proportion of associations to consider (```n_top```)<br/>
To increase result robustness and conserve memory and time, it is sometimes useful to limit the number of markers. Use 0.5 as default; Try 0.2 for GWAS with high SNP numbers; Try 1.0 for GWAS with low SNP numbers.

#### Run MDF

Download mdprune by cloning this repository using `git clone` or directly from this page.

Give execution permissions to the program.
```bash
chmod +x mdprune
```
Run MDF using function from Mergeomics_utils.R which produces and runs a bash script for MDF. Contents of the bash script that is run is shown below.

```R
source("Mergeomics_utils.R")
msea_job <- runMDF(marker_associations= "GWAS/GWAS_MSEA/Kunkle_AD.txt", # marker association file
                   marker_mapping = "mapping/MONOCYTES_EQTL.txt", # marker to gene mapping file
                   marker_dependency = "linkage/LD50.1000G.CEU.txt", # marker dependency file
                   n_top = 0.2, # top proportion of marker association file to consider
                   md_threshold = 50, # dependency (correlation/R2) cutoff, used to label output
                   mdprune = "mdprune" # path to mdprune program
                   )
```

The default output directory is Data/.

#### MDF Script (bash)
`runMDF` from above code produces and runs the below script.
```bash
#!/bin/bash
#
# Remove dependent markers and prepare an optimized marker and gene file
# for marker set enrichment analysis (MSEA). You must have the MDPrune software
# installed to use this script.
#
# Written by Ville-Petteri Makinen
#
#
MARFILE="GWAS/GWAS_MSEA/Kunkle_AD.txt"
GENFILE="mapping/MONOCYTES_EQTL.txt"
LNKFILE="linkage/LD50.1000G.CEU.txt"
OUTPATH="MSEA/Data/Kunkle_AD.MONOCYTES_EQTL/"
NTOP=0.2
echo -e "MARKER\tVALUE" > /tmp/header.txt
nice sort -r -g -k 2 $MARFILE > /tmp/sorted.txt
NMARKER=$(wc -l < /tmp/sorted.txt)
NMAX=$(echo "($NTOP*$NMARKER)/1" | bc)
nice head -n $NMAX /tmp/sorted.txt > /tmp/top.txt
cat /tmp/header.txt /tmp/top.txt > /tmp/subset.txt

# Remove SNPs in LD and create input files for SSEA.
nice mdprune /tmp/subset.txt $GENFILE $LNKFILE $OUTPATH

```

#### Outputs
These files serve as marker dependency corrected inputs for MSEA.
1. Marker association file (will be named for ex. top20.ld50.m.txt)
```
MARKER             VALUE
rs10000012	  1.9776e+00
rs1000274	  9.4846e-01
rs10003931        1.3696e+00
```

2. Marker to gene mapping file (will be named for ex. top20.ld50.g.txt)
```
GENE              MARKER
RESP18	          rs7600417
RESP18	          rs35083292
ADPRH	          rs60643107
```

The output directory is named by the trait (marker association/GWAS file name) and mapping type. This can be changed with the `label` parameter in `runMDF`. The names of the output files have the top percentage of markers (top20 for ex.) and marker dependency (linkage disequilibrium) (md50 for ex.) appended to the beginning. These are the parameters that can affect the performance of MSEA. Users may want to try different top percentages (100, 50, or 20) and marker dependency thresholds (0.5 or 0.7).

### Marker Set Enrichment Analysis
Marker set enrichment analysis (MSEA) detects pathways and networks affected by multidimensional molecular markers (e.g., SNPs, differential methylation sites) associated with a pathological condition. The pipeline can be concluded after MSEA is run, or the results can be used directly in wKDA. 

MSEA can also be used for gene level enrichment analysis (functional annotation of DEGs, transcription factor target enrichment analysis). This is outlined below. 

#### Inputs
1.```label```: output file name<br/>
2.```folder```: output folder name<br/>
3.```marfile``` and ```genfile```: marker association and marker to gene mapping files (respectively) from MDF (see outputs from MDF section). A `genfile` is not needed for transcriptome-wide and proteome-wide studies as markers are already genes. <br/>
4. ```modfile```: module/pathway file with headers 'MODULE' and 'GENE'
```
MODULE             GENE
rctm0573           APP
M5940              AMPH
Obesity_positive   CD53 
```
5. ```inffile```: description file for modules with headers 'MODULE', 'SOURCE', and 'DESCR' (optional)
```
MODULE             SOURCE           DESCR
rctm0573           reactome         Inflammasomes
M5940              biocarta         Endocytotic role of NDK, Phosphins and Dynamin
Obesity_positive   GWAS Catalog     Positive control gene set for Obesity 
```
6. ```permtype```: permutation type. This is critical. For GWAS enrichment analysis, use “gene”. For gene level enrichment analysis (epigenome-wide, transcriptome-wide, proteome-wide), use "marker". "marker" can be used for GWAS enrichment but may be bias towards genes mapped from many markers.<br/>
7. ```nperm```: number of permutations. Set to 2000 for exploratory analysis, set to 10000-20000 for formal analysis<br/>
8. ```maxoverlap```: overlap ratio threshold for merging genes with shared markers. Default is 0.33. Set to 1 for gene level enrichment analysis (epigenome-wide, transcriptome-wide, proteome-wide).
9. ```trim```: percentile taken from the beginning and end for trimming away a defined proportion of genes with significant trait association to avoid signal inflation of null background in gene permutation. Default is 0.002. Set to 0 to consider all signals in null background.

If no ```genfile``` is provided, transcriptome-wide/proteome-wide enrichment is assumed, and permtype and maxoverlap is automatically set to "marker" and 1, respectively.

#### MSEA Script for GWAS/EWAS Enrichment
<em>See `runMSEA` in Mergeomics_utils.R for a wrapper function of this.</em>
```R
# source functions or load library
source("Mergeomics_Version_1.99.0.R")
job.ssea <- list()
job.ssea$label <- "Kunkle_AD.MONOCYTES_EQTL"
job.ssea$folder <- "../results/"
job.ssea$genfile <- "MSEA/Data/Kunkle_AD.MONOCYTES_EQTL/top20.ld50.g.txt" # marker to gene file
job.ssea$marfile <- "MSEA/Data/Kunkle_AD.MONOCYTES_EQTL/top20.ld50.m.txt" # marker association file		
job.ssea$modfile <- "../resources/genesets/kbr.mod.txt"
job.ssea$inffile <- "../resources/genesets/kbr.info.txt"
job.ssea$permtype <- "gene" # for EWAS, set this to "marker"
job.ssea$nperm <- 10000
job.ssea$maxoverlap <- 0.33 # for EWAS, set this to 1
job.ssea$trim <- 0.002 # default is 0.002, set to 0 for no trimming, users can try increasing to 0.005 to dampen inflation further
job.ssea <- ssea.start(job.ssea)
job.ssea <- ssea.prepare(job.ssea)
job.ssea <- ssea.control(job.ssea)
job.ssea <- ssea.analyze(job.ssea)
job.ssea <- ssea.finish(job.ssea)
```

#### MSEA Script for TWAS, PWAS, MWAS Enrichment
If markers are genes, then a marker to gene mapping file is not needed. If no marker to gene mapping file is provided, TWAS/PWAS/MWAS enrichment will be assumed and the appropriate parameters will be set.

```R
# source functions or load library
source("Mergeomics_Release1.19.0_beta.R")
job.ssea <- list()
job.ssea$label <- "DEGs_enrichment"
job.ssea$folder <- "../results/"
job.ssea$marfile <- "DEGs.txt" # marker association file
job.ssea$modfile <- "../resources/genesets/kbr.mod.txt"
job.ssea$inffile <- "../resources/genesets/kbr.info.txt"
job.ssea$nperm <- 10000
job.ssea <- ssea.start(job.ssea)
job.ssea <- ssea.prepare(job.ssea)
job.ssea <- ssea.control(job.ssea)
job.ssea <- ssea.analyze(job.ssea)
job.ssea <- ssea.finish(job.ssea)
```


#### Outputs

1. Results file<br/>

MODULE | P | FREQ | ZSCORE	| CHI	| CHI_SD | CHI_SE	| NGENES	| NMARKER	| DENSITY	| FDR	| TOPGENES | TOPMARKERS | TOPVALUES
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- 
P1 | 2.44E-18 | 0 | 8.65 | 5.74 |	2.23 | 0.70 | 82 | 1497 | 1 | 1.54E-16 | ITGA6;THBS2;TNR;CD44;ITGB4 | cg24530074;cg04476508;cg26359730;cg15695194;cg20942310 | 10.00;8.89;8.87;8.78;7.89
P2 | 3.36E-17 |	0 |	8.35 | 5.13 |	1.88 | 0.59 | 196 | 2875 | 1 | 1.59E-15 | PIK3CD;BIRC3;ITGA6;THBS2;TNR | cg09706586;cg14296561;cg24530074;cg04476508;cg26359730 | 10.44;10.11;10.00;8.89;8.87
P3 | 7.51E-11 |	0 |	6.40 | 3.71 | 1.43 | 0.45 | 82 | 924 | 1 | 2.84E-09 | LMNA;ITGA6;PRKAG2;SGCA;CACNA1C | cg09820673;cg24530074;cg22615330;cg04582295;cg03047070 | 10.90;10.00;8.79;8.46;8.16

2. Details file<br/>

MODULE | FDR | NMARKER | MARKER | VALUE | DESCR
--- | --- | --- | --- | --- | ---
CHD_positive | 0.00% | ABO | 13 | rs635634 | 6.44 | Positive control gene set for coronary heart disease
rctm0171 | 0.05% | GSTM3 | 7 | rs112479902 | 2.77 | Biological oxidations
M7098 | 0.79% | ITGB6 | 16 | rs1563575 | 6.17 | ECM-receptor interaction

3. Genes file<br/>

GENE | SCORE | NMARKER | MARKER | VALUE 
--- | --- | --- | --- | --- 
CPNE1 | 1.21789059181483 | 16 | rs12480111 | 3.7959 
COX7A1 | -0.307584552536074 | 5 | rs12611259 | 1.9208 
COG4 | -0.358654267820555 | 8 | 16 | rs75076686 | 1.8239

4. Pvalues file<br/>

MODULE | P.Study | FDR.Study | DESCR 
--- | --- | --- | --- 
M14933 | 2.11e-07 | 0.0000 | Steroid hormone biosynthesis
rctm0261 | 1.30e-06 | 0.0001 | Complement cascade
rctm1301 | 1.24e-04 | 0.0043 | Transferrin endocytosis and recycling

See files in example_output of this github.


### Module Merging
This optional step merges redundant pathways (pathways with significant sharing of member genes) into supersets (e.g., "Extracellular matrix" and "Cell adhesion" may be merged). The degree of merging can be modified based on the `rmax` (max gene sharing ratio that modules will remain independent). A recommended `rmax` is 0.33, which means that modules sharing overlap above this ratio will be merged. For more module merging, set the `rmax` lower; for less module merging, set the `rmax` higher (stricter ratio cutoff).

We recommend this step to reduce redundancy in significant modules and to reduce redundancy in input to KDA. When running KDA, the analysis is run on each module in parallel. 

#### Method 1: Use `merge_modules` function using result file from MSEA or vector of modules
This function is in Mergeomics_utils.R.
```R
source("Mergeomics.R")
source("Mergeomics_utils.R")
options(stringsAsFactors = FALSE)

# minimum inputs. msea_res can be a vector of module names or the ".results.txt" result file from MSEA
merge_modules(msea_res="sample_inputs/Sample_Meta.MSEA_modules_full_result.txt", 
              modfile_path = "sample_inputs/Kegg.txt, # MODULE GENE file
	      fdr_cutoff = 0.01) # not required input but highly recommended for downstream analysis with KDA
# OR using a vector of modules
res <- read.delim("sample_inputs/Sample_Meta.MSEA_modules_full_result.txt)
sig_mods <- res$MODULE[res$FDR<0.01]
merge_modules(msea_res=sig_mods, 
              modfile_path = "sample_inputs/Kegg.txt")
	      
# All parameters set
merge_modules(msea_res="~/Downloads/Meta_MSEA_TTK1A14P1X/Sample_Meta.MSEA_modules_full_result.txt", 
              rcutoff = 0.33, # this is the default value
              label = "Psoriasis", 
	      fdr_cutoff = 0.01, 
	      output_dir="Merged_modules/", 
	      modfile_path = "sample_inputs/Kegg.txt",
              infofile_path = "sample_inputs/KEGG_info.txt")

```

#### Method 2: Use `ssea2kda` function using MSEA job with results
The job returned from `ssea.finish` is used in the `ssea2kda` function (see MSEA script above).
```R
# last step of MSEA
job.ssea <- ssea.finish(job.ssea)
# rmax can be adjusted as well
job.kda <- ssea2kda(job.msea, filter=0.05) # filter is the FDR cutoff
```

#### Outputs
Merged modules with ",.." appended to the module name indicates a "superset" and the member modules are detailed in the "OVERLAP" column. This file can be directly used in KDA (KDA nodes input requires "MODULE" and "NODE" columns).

1. Module file ("merged_...mod.txt" from `merge_modules` or "msea2kda.modules.txt" from `ssea2kda`)
```
MODULE		GENE	OVERLAP			NODE	
rctm0089,..	CYBA	rctm0089,rctm0246	CYBA
rctm0089,..	ANAPC1	rctm0089,rctm0246	ANAPC1
rctm0693	ACHE	rctm0693		ACHE

```

2. Info file (if set infofile_path parameter in `merge_modules` function)

```
MODULE		SOURCE		DESCR
rctm0089,..	reactome	Adaptive Immune System
rctm0693	reactome	Metabolism of proteins
```


### Weighted Key Driver Analysis
Weighted key driver analysis (wKDA) selects key regulator genes of the disease related gene sets using gene network topology and edge weight information. wKDA first screens the network for candidate hub genes and then the disease gene-sets are overlaid onto the subnetworks of the candidate hubs to identify key drivers whose neighbors are enriched with disease genes. 

#### Inputs
1. Module file from merge modules or a file containing 'MODULE' and 'NODE' columns. See output #1 from Module Merging. 

2. Network file 
```
TAIL		HEAD        WEIGHT
A1BG		SNHG6       1
A1BG		UNC84A      1
A1CF		KIAA1958    1
```


#### wKDA Script (R)
<em>See `runKDA` in Mergeomics_utils.R for a wrapper function of this.</em>
```R
job.kda <- list()
job.kda$label <- "wKDA"
job.kda$folder <- "Results" # parent folder for results
# Input a network
# columns: TAIL HEAD WEIGHT
job.kda$netfile<-"bayesian_network.txt"
# Tab delimited text file of gene set containing two columns: MODULE, NODE
# Outputs from Module Merge script can be directly used
job.kda$modfile<-"merged_modules.txt"

# 0.0-1.0 - 0.0 means not factoring in edge weights, 0.5 means partial influence,
# and 1.0 means full influence
job.kda$edgefactor<-1
# The searching depth for the KDA
job.kda$depth<-1
# 0 means we do not consider the directions of the regulatory interactions
# while 1 is opposite.
job.kda$direction <- 1
job.kda$nperm <- 2000

## Let's run KDA!
job.kda <- kda.configure(job.kda)
job.kda <- kda.start(job.kda)
job.kda <- kda.prepare(job.kda)
job.kda <- kda.analyze(job.kda)
job.kda <- kda.finish(job.kda)

job.kda <- kda2cytoscape(job.kda)
```

### Outputs
Following completion of the analysis, key driver results are deposited in the "kda" subdirectory and if significant (FDR<0.05) key drivers are found, Cytoscape ready files are deposited in the "cytoscape" subdirectory. See files in example_output of this github.


## Citations
Shu, L., et al., Mergeomics: multidimensional data integration to identify pathogenic perturbations to biological systems. BMC Genomics, 2016. 17(1): p. 874.

Ding, J., et al., Mergeomics 2.0: a web server for multi-omics data integration to elucidate disease networks and predict therapeutics. Nucleic Acids Research, 2021.


