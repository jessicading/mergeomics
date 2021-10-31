*** Please use Mergeomics_Release1.19.0_beta.R and report any issues

# Mergeomics
Integrative Network Analysis of Omics Data<br/>
Developed by: Ville-Petteri Makinen, Le Shu, Yuqi Zhao, Zeyneb Kurt, Bin Zhang, Xia Yang

#### Table of Contents
1. [About](#about)<br/>
2. [Tutorial](#tutorial)<br/>
    1. [Marker Dependency Filtering](#marker-dependency-filtering)
    2. [Marker Set Enrichment Analysis](#marker-set-enrichment-analysis)
    3. [Module Merging](#module-merging)
    4. [Weighted Key Driver Analysis](#weighted-key-driver-analysis)
3. [Workflow using wrapper functions](#workflow-using-wrapper-functions)

## About
Mergeomics is an open source R-based pipeline for multi-dimensional integration of omics data to identify disease-associated pathways and networks. Genes whose network neighborhoods are over-represented with disease associated genes are deemed key drivers and can be targeted in further mechanistic studies.

![abstract figure](https://github.com/jessicading/mergeomics/blob/master/Abstract.jpg?raw=true)

More information can be found in the [methods paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3198-9), [documentation](http://bioconductor.org/packages/release/bioc/html/Mergeomics.html), [web server](http://mergeomics.research.idre.ucla.edu), and [web server paper](https://academic.oup.com/nar/article/49/W1/W375/6287846).

## Tutorial
This tutorial illustrates the workflow for running the Mergeomics pipeline which includes Marker Dependency Filtering (MDF), Marker Set Enrichment Analysis (MSEA), and Key Driver Analysis (KDA). For those who prefer a web interface, please use our [web server](http://mergeomics.research.idre.ucla.edu).

Input and output file examples and individual scripts for each step are given. A more streamlined workflow is displayed at the end which uses wrapper functions in the Mergeomics_utils.R source function.


#### Notes
- A new version of the Mergeomics R package will be released which fixes minor bugs and improves results output. We suggest to using `source(Mergeomics_Release1.19.0_beta.R)` to use the updated functions.
- All input files should be tab delimited .txt files.
- The term 'marker' is used to refer to any type of biological marker including loci, methylation/epigenetic sites, transcripts/genes, proteins, metabolites, etc.

### Marker Dependency Filtering

Marker dependency filtering (MDF) removes dependent markers and prepares an optimized marker and gene file for marker set enrichment analysis (MSEA). Dependency between markers can be a result of linkage disequilibrium, for example. Download the mdprune executable in this repository and provide the path to the program in the script below. This is a C++ program and requires a Linux distribution. If you cannot obtain the environment, MDF can be run using our [web server](http://mergeomics.research.idre.ucla.edu) (Run Mergeomics > Individual GWAS Enrichment).

For epigenomics, transcriptomics, proteomics, or metabolomics data, MSEA can be run as a first step (it is recommended to use the `runMSEA` function to simplify the analysis). MDF can be skipped for GWAS as well but it is not recommended.

#### Inputs
1. ```MARFILE```: Disease/Phenotype association data <br/> 
This file must have two columns named 'MARKER' and 'VALUE' where value denotes the strength of the association to the trait. This is usually -log10 p-values but can be other values such as effect size or fold change. Higher values must indicate higher association.
```
MARKER             VALUE
rs4747841         0.1452
rs4749917         0.1108
rs737656          1.3979
```
2. ```GENFILE```: Mapping file to connect markers from association data to genes<br/>
(eQTLs, distance based, etc.)
```
GENE              MARKER
TSPAN6            rs1204389
TNMD              rs113630520
SCYL3             rs72691775
```
3. ```LNKFILE```: Marker Dependency File (ex. Linkage Disequilibirum)<br/>
LD files can be obtained from the [HapMap consortium](https://www.sanger.ac.uk/resources/downloads/human/hapmap3.html). 
```
MARKERa            MARKERb              WEIGHT            
rs143225517       rs3094315           0.962305
rs4475691         rs950122            1
rs4475691         rs3905286           0.921467
```
4. ```OUTPATH```: Output Directory
5. ```NTOP```: Top proportion of associations to consider <br/>
To increase result robustness and conserve memory and time, it is sometimes useful to limit the number of markers. Use 0.5 as default; Try 0.2 for GWAS with high SNP numbers; Try 1.0 for GWAS with low SNP numbers

#### MDF Script (bash)
<em>See `runMDF` in Mergeomics_utils.R for a wrapper function of this.</em>
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
MARFILE="../GWAS/DIAGRAMstage1_T2D.txt"
GENFILE="../resources/mapping/esnps/Adipose_Subcutaneous.txt"
LNKFILE="../resources/linkage/LD50.1000G.CEU.txt"
OUTPATH="../MSEA/Data/DIAGRAMstage1_T2D.Adipose_Subcutaneous/"
NTOP=0.2
echo -e "MARKER\tVALUE" > /tmp/header.txt
nice sort -r -g -k 2 $MARFILE > /tmp/sorted.txt
NMARKER=$(wc -l < /tmp/sorted.txt)
NMAX=$(echo "($NTOP*$NMARKER)/1" | bc)
nice head -n $NMAX /tmp/sorted.txt > /tmp/top.txt
cat /tmp/header.txt /tmp/top.txt > /tmp/subset.txt

# Remove SNPs in LD and create input files for SSEA.
nice /u/home/j/jading/project-xyang123/GWAS/MDPRUNE/mdprune /tmp/subset.txt $GENFILE $LNKFILE $OUTPATH

```

#### Outputs
These files serve as inputs for MSEA.
1. Gene file
```
GENE              MARKER
RESP18	          rs7600417
RESP18	          rs35083292
ADPRH	          rs60643107
```
2. Loci file
```
MARKER             VALUE
rs10000012	  1.9776e+00
rs1000274	  9.4846e-01
rs10003931        1.3696e+00
```

### Marker Set Enrichment Analysis
Marker set enrichment analysis (MSEA) detects pathways and networks affected by multidimensional molecular markers (e.g., SNPs, differential methylation sites) associated with a pathological condition. The pipeline can be concluded after MSEA is run, or the results can be used directly in wKDA. 

MSEA can also be used for gene level enrichment analysis (functional annotation of DEGs, transcription factor target enrichment analysis). This is outlined below. See the `runMSEA` function from Mergeomics_utils.R for a wrapper function of this script.

#### Inputs
1.```label```: output file name<br/>
2.```folder```: output folder <br/>
3.```genfile``` and ```marfile```: Gene and marker files (respectively) from MDF (see outputs #1 and #2 from MDF section). For gene level enrichment analysis, a "fake" gene (mapping) file can be made. <br/>
4. ```modfile```: module/pathway file with headers 'MODULE' and 'GENE'
```
MODULE             GENE
rctm0573           APP
M5940              AMPH
Obesity_positive   CD53 
```
5. ```inffile```: description file for modules with headers 'MODULE', 'SOURCE', and 'DESCR'
```
MODULE             SOURCE           DESCR
rctm0573           reactome         Inflammasomes
M5940              biocarta         Endocytotic role of NDK, Phosphins and Dynamin
Obesity_positive.  GWAS Catalog     Positive control gene set for Obesity 
```
6. ```permtype```: This is critical. For GWAS enrichment analysis, use “gene”. For gene level enrichment analysis, use "locus"<br/>
7. ```nperm```: Set to 2000 for exploratory analysis, set to 10000 for formal analysis<br/>
8. ```maxoverlap```: Default is 0.33. Set to 1 for gene level enrichment analysis.

#### MSEA Script for GWAS/EWAS Enrichment
<em>See `runMSEA` in Mergeomics_utils.R for a wrapper function of this.</em>
```R
# source functions or load library
source("Mergeomics_Release1.19.0_beta.R")
job.ssea <- list()
job.ssea$label <- "gwas_enrichment"
job.ssea$folder <- "../results/"
job.ssea$genfile <- "eqtls.txt" # marker to gene file
job.ssea$marfile <- "gwas_loci.txt" # marker association file		
job.ssea$modfile <- "../resources/genesets/kbr.mod.txt"
job.ssea$inffile <- "../resources/genesets/kbr.info.txt"
job.ssea$permtype <- "gene"
job.ssea$nperm <- 10000
job.ssea$maxoverlap <- 0.33
job.ssea <- ssea.start(job.ssea)
job.ssea <- ssea.prepare(job.ssea)
job.ssea <- ssea.control(job.ssea)
job.ssea <- ssea.analyze(job.ssea)
job.ssea <- ssea.finish(job.ssea)
```

#### MSEA Script for TWAS, PWAS, MWAS Enrichment
If markers are already functional genes, then a marker to gene mapping file is not needed. If no marker to gene mapping file is provided, TWAS/PWAS/MWAS enrichment will be assumed and the appropriate parameters will be set.

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
1. Details file<br/>

MODULE | FDR | NMARKER | MARKER | VALUE | DESCR
--- | --- | --- | --- | --- | ---
CHD_positive | 0.00% | ABO | 13 | rs635634 | 6.44 | Positive control gene set for coronary heart disease
rctm0171 | 0.05% | GSTM3 | 7 | rs112479902 | 2.77 | Biological oxidations
M7098 | 0.79% | ITGB6 | 16 | rs1563575 | 6.17 | ECM-receptor interaction

2. Genes file<br/>

GENE | SCORE | NMARKER | MARKER | VALUE 
--- | --- | --- | --- | --- 
CPNE1 | 1.21789059181483 | 16 | rs12480111 | 3.7959 
COX7A1 | -0.307584552536074 | 5 | rs12611259 | 1.9208 
COG4 | -0.358654267820555 | 8 | 16 | rs75076686 | 1.8239

3. Pvalues file<br/>

MODULE | P.Study | FDR.Study | DESCR 
--- | --- | --- | --- 
M14933 | 2.11e-07 | 0.0000 | Steroid hormone biosynthesis
rctm0261 | 1.30e-06 | 0.0001 | Complement cascade
rctm1301 | 1.24e-04 | 0.0043 | Transferrin endocytosis and recycling

4. Results file<br/>

MODULE | P | FREQ | NGENES	| NMARKER	| DENSITY	| FDR	| DESCR
--- | --- | --- | --- | --- | --- | --- | ---
rctm1372 | 5.71-23| 0 | 36 | 1337 |	1	| 3.084e-21 |	WNT ligand biogenesis and trafficking
rctm0476 | 8.22e-25 |	0 |	79 | 1486 |	0.99 | 8.87e-23 | GPCR ligand binding
rctm0917 | 0.000175 |	0.00031 |	12 | 435 | 1 | 0.0026 | Protein folding

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
              modfile_path = "sample_inputs/Kegg1.txt,
	      fdr_cutoff = 0.01) # not required input but highly recommended for downstream analysis with KDA
# OR using a vector of modules
res <- read.delim("sample_inputs/Sample_Meta.MSEA_modules_full_result.txt)
sig_mods <- res$MODULE[res$FDR<0.01]
merge_modules(msea_res=sig_mods, 
              modfile_path = "sample_inputs/Kegg1.txt")
	      
# All parameters set
merge_modules(msea_res="~/Downloads/Meta_MSEA_TTK1A14P1X/Sample_Meta.MSEA_modules_full_result.txt", 
              rmax = 0.33, # this is the default value
              modfile_path = "sample_inputs/Kegg1.txt",
              label = "Psoriasis", 
	      fdr_cutoff = 0.01, 
	      output_Dir = "Merged/", 
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
job.kda$label<-"wKDA"
job.kda$folder<-"Results" ## parent folder for results
## Input a network
## columns: TAIL HEAD WEIGHT
job.kda$netfile<-"bayesian_network.txt"
## Tab delimited text file of gene set containing two columns: MODULE, NODE
## Outputs from Module Merge script can be directly used
job.kda$modfile<-"merged_modules.txt"

## "0" means we do not consider edge weights while 1 is opposite.
job.kda$edgefactor<-1
## The searching depth for the KDA
job.kda$depth<-1
## 0 means we do not consider the directions of the regulatory interactions
## while 1 is opposite.
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
If significant key drivers were found, "cytoscape" and "kda" directories are made containing input files for cytoscape and kda results. 


## Workflow using wrapper functions

This functions can be found in Mergeomics_utils.R. The most common values for certain settings are set as default parameters in these wrapper functions. Please review above to see descriptions on these different parameters and they may be tweaked in the following functions. Generally the minimum input is depicted below.

```R
source("Mergeomics.R")
source("Mergeomics_utils.R")
```

### GWAS Enrichment

```R
runMDF(MARFILE = "./GWAS/Kunkle_AD.txt",
       GENFILE = "./mapping/Brain_Hippocampus.eQTL.txt", 
       LNKFILE = "./linkage/LD50.1000G.CEU.txt", 
       output_dir = "./MSEA/Data/", # if more than one directory, must create beforehand such as in this case
       mdprune = "./GWAS/mdprune")
runMSEA(MDF_output_dir = "./MSEA/Data/Kunkle_AD.Brain_Hippocampus.eQTL/",
        marker_set="./HP_MSEA_DEGs.txt")
runKDA(MSEA_results_dir = "./Results/msea/", #contains only one set of results with one "-.results.txt" file
       marker_set="./HP_MSEA_DEGs.txt",
       network="./network/brain_network.txt")
# OR
runKDA(nodes_file = "./Results/msea/Kunkle_AD.Brain_Hippocampus.eQTL.results.txt",
       marker_set="./HP_MSEA_DEGs.txt",
       network="./network/brain_network.txt")
```

### Transcriptome/Proteome/Epigenome/Metabolome Enrichment
```R
runMSEA(association = "./Microglia_DEGs.txt", # has DEGs in LOCUS column and -log10 transformed p values in VALUE column
        marker_set="./genesets/kbr.mod.txt",
	info="./genesets/kbr.info.txt") # info not required
runKDA(MSEA_results_dir = "./Results/msea/Microglia_DEGs.results.txt",
       marker_set="./genesets/kbr.mod.txt",
       network="./network/brain_network.txt")
```

### Run KDA on genes
```R
runKDA(nodes_file = "./HP_MSEA_DEGs.txt",
       network="./network/brain_network.txt")
```

### Skip MDF for GWAS and run MSEA (not recommended for formal analysis) 
```R
runMSEA(association_file = "./GWAS/Kunkle_AD.txt", 
        perc_top_associations = 0.5,
        mapping_file = "./mapping/Brain_Hippocampus.eQTL.txt",
        marker_set="./HP_MSEA_DEGs.txt")
```

