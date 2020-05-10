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

More information can be found in the [paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3198-9), [documentation](http://bioconductor.org/packages/release/bioc/html/Mergeomics.html), and [web server](http://mergeomics.research.idre.ucla.edu).

## Tutorial
This tutorial describes the workflow for running the Mergeomics pipeline with individual scripts (using Mergeomics functions) and requires a minimum level of coding experience. This method allows for more customization and finetuning of the parameters. This is a response to the R package currently being under construction. 

For each step, example input and output files are given, the script, and discussion of the parameters.

There are example files for each input and output. 

#### Notes
All files used are tab delimited text files.

In previous documentations and tutorials, "MARKER" was used as a column name in input files, and this tutorial will still describe associations as "markers" but in input and output files, "LOCUS" will take the place of "MARKER". This is because the particular scripts being used as examples here use "LOCUS" and not "MARKER" column names. This can be modified by the user in their own analysis. 

All scripts except MDF need to source "Mergeomics.R".

### Marker Dependency Filtering (MDF)

MDF removes dependent markers and prepares an optimized marker and gene file for marker set enrichment analysis (MSEA). You must have the ldprune software installed for this script. The inputs are enumerated in the scripts as enumerated below. 

As of writing, MDF is done mostly for GWAS data (to correct for linkage disequilibrium). If starting from transcriptomic, proteomic, epigenomic, or metabolomic data, start from MSEA (it is recommended to use the `runMSEA` function to simplify the analysis). 

#### Inputs
1. ```LOCFILE```: Disease/Phenotype Associated Data <br/> 
This file must have two columns named 'LOCUS' and 'VALUE' where value denotes the association to the trait of interest (p-value). The p-values must be negative log (base 10) transformed (-log P). The higher the value, the stronger the association.
```
LOCUS             VALUE
rs4747841         0.1452
rs4749917         0.1108
rs737656          1.3979
```
2. ```GENFILE```: Mapping File to Connect Markers from Association Data to Genes<br/>
(ex. eQTLs and/or Encode Information)
```
GENE              LOCUS
TSPAN6            rs1204389
TNMD              rs113630520
SCYL3             rs72691775
```
3. ```LNKFILE```: Marker Dependency File (ex. Linkage Disequilibirum)<br/>
LD files can be obtained from the [HapMap consortium](https://www.sanger.ac.uk/resources/downloads/human/hapmap3.html). 
```
LOCUSa            LOCUSb              WEIGHT            
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
LOCFILE="../GWAS/DIAGRAMstage1_T2D.txt"
GENFILE="../resources/mapping/esnps/Adipose_Subcutaneous.txt"
LNKFILE="../resources/linkage/LD50.1000G.CEU.txt"
OUTPATH="../MSEA/Data/DIAGRAMstage1_T2D.Adipose_Subcutaneous/"
NTOP=0.2
echo -e "LOCUS\tVALUE" > /tmp/header.txt
nice sort -r -g -k 2 $LOCFILE > /tmp/sorted.txt
NMARKER=$(wc -l < /tmp/sorted.txt)
NMAX=$(echo "($NTOP*$NMARKER)/1" | bc)
nice head -n $NMAX /tmp/sorted.txt > /tmp/top.txt
cat /tmp/header.txt /tmp/top.txt > /tmp/subset.txt

# Remove SNPs in LD and create input files for SSEA.
nice /u/home/j/jading/project-xyang123/GWAS/MDPRUNE/ldprune /tmp/subset.txt $GENFILE $LNKFILE $OUTPATH

```

#### Outputs
These files serve as inputs for MSEA.
1. Gene file
```
GENE              LOCUS
RESP18	          rs7600417
RESP18	          rs35083292
ADPRH	          rs60643107
```
2. Loci file
```
LOCUS             VALUE
rs10000012	  1.9776e+00
rs1000274	  9.4846e-01
rs10003931        1.3696e+00
```

### Marker Set Enrichment Analysis (MSEA)
MSEA detects pathways and networks affected by multidimensional molecular markers (e.g., SNPs, differential methylation sites) associated with a pathological condition. The pipeline can be concluded after MSEA is run, or the results can be used directly in wKDA. 

MSEA can also be used for gene level enrichment analysis (functional annotation of DEGs, transcription factor target enrichment analysis) with different parameter settings. This is outlined below. Alternatively, you can use the `runMSEA` wrapper function without specifying a mapping_file.

#### Inputs
1.```label```: output file name<br/>
2.```folder```: output folder <br/>
3.```genfile``` and ```locfile```: Gene and loci files (respectively) from MDF (see outputs #1 and #2 from MDF section). For gene level enrichment analysis, a "fake" gene (mapping) file can be made. <br/>
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

#### MSEA Script (R)
<em>See `runMSEA` in Mergeomics_utils.R for a wrapper function of this.</em>
```R
job.ssea <- list()
job.ssea$label <- "DIAGRAMstage2_T2D.Adipose_Subcutaneous"
job.ssea$folder <- "../results/"
job.ssea$genfile <- "./Data/Adipose_Subcutaneous/genes.txt"
job.ssea$locfile <- "./Data/Adipose_Subcutaneous/loci.txt"		
job.ssea$modfile <- "../resources/genesets/kbr.mod.txt"
job.ssea$inffile <- "../resources/genesets/kbr.info.txt"
job.ssea$permtype <- "gene"
job.ssea$nperm <- 10000
job.ssea$maxoverlap <- 0.33
job.ssea <- ssea.start(job.ssea)
job.ssea <- ssea.prepare(job.ssea)
job.ssea <- ssea.control(job.ssea)
job.ssea <- ssea.analyze(job.ssea,trim_start=0.005,trim_end=0.995)
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
This step merges redundant pathways (pathways with significant sharing of member genes) into relatively independent gene sets. The input for this script as exactly written below is the results file from MSEA (see output #4 from Marker Set Enrichment Analysis).

The minimum input for this script is a list of modules. If the user does post-analysis of the results file where the significant modules are already extracted, a text file with only one column 'MODULE' is sufficient for this script and the script can easily be modified to allow this. In the script below, lines with "##" are necessary if inputting a list of already filtered modules. Everything after the section "Merge modules before 2nd SSEA" requires only the list of modules. 

The output module file can be used as input for wKDA. 

This step is optional. Sigificant modules found from MSEA can be used for wKDA (with MODULE and NODE columns where NODE contains the genes).

#### Module Merge Script (R)
<em>See `merge_modules` in Mergeomics_utils.R for a wrapper function of this.</em>
```R
plan = c()
plan$folder = "../msea_results"
plan$label = "DIAGRAMstage1_T2D.Adipose_Subcutaneous" 
plan$modfile="../resources/genesets/kbr.mod.txt"
plan$inffile="../resources/genesets/kbr.info.txt"
#=====================================================
# first pick the significantly enriched modules for 1st ssea:
FDR=0.05
pool=c() ##
file.name=paste0(plan$folder, "/",plan$label, ".txt") ##
if(file.exists(file.name)){
  aa=(read.table(file.name, header=T, sep='\t', check.names=F, quote=NULL)) ##
  if(length(which(aa[,"MODULE"] == "_ctrlA")) >0) aa=aa[-which(aa[,"MODULE"] == "_ctrlA"),] 
  if(length(which(aa[,"MODULE"] == "_ctrlB")) >0) aa=aa[-which(aa[,"MODULE"] == "_ctrlB"),] 
  aa=aa[which( (as.numeric(aa[,"FDR"])< FDR)  ), ] 
  if (nrow(aa) > 0) pool=unique(c(pool, as.character(aa[,"MODULE"]))) ##
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
  mdfile="modules_canonical_plus_coexm_mods.kbr.txt"; mifile="kbr.modules.info.txt" 
  
  write.table(moddata, paste0(plan$folder,"/merged_", mdfile),  
              sep='\t', col.names=T, row.names=F, quote=F) 
  
  write.table(moddatainfo, paste0(plan$folder,"/merged_", mifile),  
              sep='\t', col.names=T, row.names=F, quote=F) 
  
  #=====================================================================================# 
  #    Apply 2nd SSEA with the merged file to check that pathways are still significant #
  #=====================================================================================#
  
}
```

#### Outputs
Merged modules with ",.." appended to the module name indicates a "superset" and the member modules are detailed in the "OVERLAP" column.

1. Module file
```
MODULE		GENE	OVERLAP			NODE	
rctm0089,..	CYBA	rctm0089,rctm0246	CYBA
rctm0089,..	ANAPC1	rctm0089,rctm0246	ANAPC1
rctm0693	ACHE	rctm0693		ACHE

```

2. Info file

```
MODULE		SOURCE		DESCR
rctm0089,..	reactome	Adaptive Immune System
rctm0693	reactome	Metabolism of proteins
```


### Weighted Key Driver Analysis (wKDA)
wKDA selects key regulator genes of the disease related gene sets using gene network topology and edge weight information. wKDA first screens the network for candidate hub genes and then the disease gene-sets are overlaid onto the subnetworks of the candidate hubs to identify key drivers whose neighbors are enriched with disease genes. 

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

moddata <- tool.read(job.kda$modfile)
mod.names <- unique(moddata$MODULE)
moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))),]
## save this to a temporary file and set its path as new job.kda$modfile:
tool.save(moddata, "subsetof.supersets.txt")
job.kda$modfile <- "subsetof.supersets.txt"

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
runMDF(LOCFILE = "./GWAS/Kunkle_AD.txt",
       GENFILE = "./mapping/Brain_Hippocampus.eQTL.txt", 
       LNKFILE = "./linkage/LD50.1000G.CEU.txt", 
       output_dir = "./MSEA/Data/", # if more than one directory, must create beforehand such as in this case
       ldprune = "./GWAS/ldprune")
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

#### Extra notes
We understand the confusion with 'LOCUS' as the column name as opposed to 'MARKER'. In the scripts and ldprune program provided in this github page, 'LOCUS' needs to be the column name as opposed to 'MARKER'. We hope to change to 'MARKER' soon to reflect that SNP, gene, protein, and metabolite information can be used. To clarify, 'LOCUS' can contain <b>any</b> type of information, though LOCUS implies location in genome (SNP). 
