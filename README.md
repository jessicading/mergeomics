# Mergeomics
Integrative Network Analysis of Omics Data

[Custom foo description](#header)

```bash
LOCFILE="/u/home/j/jading/project-xyang123/GWAS/GWAS_MSEA/DIAGRAMstage1_T2D.txt"
GENFILE="/u/home/j/jading/project-xyang123/resources/mapping/esnps/periph_18esnps.txt"
LNKFILE="/u/home/j/jading/project-xyang123/resources/linkage/LD50.1000G.CEU.txt"
OUTPATH="/u/home/j/jading/project-xyang123/MSEA/Data/DIAGRAMstage1_T2D.periph_18esnps/"
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

In line ```R string1 = "string" ```

# Header

