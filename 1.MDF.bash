LOCFILE=""
GENFILE=""
LNKFILE=""
OUTPATH=""
NTOP=0.2
echo -e "LOCUS\tVALUE" > /tmp/header.txt
nice sort -r -g -k 2 $LOCFILE > /tmp/sorted.txt
NMARKER=$(wc -l < /tmp/sorted.txt)
NMAX=$(echo "($NTOP*$NMARKER)/1" | bc)
nice head -n $NMAX /tmp/sorted.txt > /tmp/top.txt
cat /tmp/header.txt /tmp/top.txt > /tmp/subset.txt


# Remove SNPs in LD and create input files for SSEA.
nice /u/home/j/jading/project-xyang123/GWAS/MDPRUNE/ldprune /tmp/subset.txt $GENFILE $LNKFILE $OUTPATH
