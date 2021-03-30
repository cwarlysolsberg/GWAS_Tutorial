FILESET=HapMap_3_r3_1
GENO=0.02 
MIND=0.02 
MAF=0.05

cd HOME/{user}/{path/folder containing your files}
ls

## Visualizing and correcting for missingness
plink -bfile $FILESET --missing

Rscript Histimiss.R

## Filter on missingness 
plink --bfile $FILESET --geno 0.2 --make-bed --out $FILESET.geno.temp
plink --bfile $FILESET.geno.temp --mind 0.2 --make-bed --out $FILESET.mind.temp
plink --bfile $FILESET.mind.temp --geno $GENO --make-bed --out $FILESET.geno
plink --bfile $FILESET.geno --mind $MIND --make-bed --out $FILESET.mind

##Check for sex-discrepancies
plink --bfile $FILESET.mind --check-sex 

Rscript XHomozygosity.R

grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
plink --bfile $FILESET.mind --remove sex_discrepancy.txt --make-bed --out $FILESET.rem

##Limit data to autosomal SNPs only
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' $FILESET.rem.bim > snp_1_22.txt
plink --bfile $FILESET.rem --extract snp_1_22.txt --make-bed --out $FILESET.autosome

## Perform MAF check and filter by MAF Threshold
plink --bfile $FILESET.autosome --freq --out MAF_check

Rscript MafDist.R


plink --bfile $FILESET.autosome --maf $MAF --make-bed --out $FILESET.maf

##Generate a pruned subset of SNPs that are in approximate linkage disequilibrium
plink --bfile $FILESET.maf --indep-pairwise 50 5 0.5


##Compute method-of-moments F coefficient estimates**
plink --bfile $FILESET.maf --extract plink.prune.in --het --out R_check

Rscript Heterozygosity.R

sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}' > het_fail_ind.txt
plink --bfile $FILESET.maf --remove het_fail_ind.txt --make-bed --out $FILESET.het_fix

#mkdir ../2_Population_stratification

cp $FILESET.het_fix.bed ../2_Population_stratification/$FILESET.qcout.bed
cp $FILESET.het_fix.bim ../2_Population_stratification/$FILESET.qcout.bim
cp $FILESET.het_fix.fam ../2_Population_stratification/$FILESET.qcout.fam
cp plink.prune.in ../2_Population_stratification