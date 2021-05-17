FILESET=1000genomes_nomissing.genotypes
GENO=0.02
INDV=0.02
MAF=0.05
FILE_QC=HapMap_3_r3_1.qcout


wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz
md5sum ALL.2of4intersection.20100804.genotypes.vcf.gz
plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --make-bed --out 1000genomes.genotypes

plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --make-bed --out 1000genomes.genotypes
plink --bfile 1000genomes.genotypes --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out 1000genomes_nomissing.genotypes


plink --bfile $FILE_1K --geno 0.2 --make-bed --out $FILE_1K.geno.temp
plink --bfile $FILE_1K.geno.temp --mind 0.2 --allow-no-sex --make-bed --out $FILE_1K.geno.mind.temp
plink --bfile $FILE_1K.geno.mind.temp --geno $GENO --make-bed --out $FILE_1K.geno
plink --bfile $FILE_1K.geno --mind $INDV --allow-no-sex --make-bed --out $FILE_1K.geno.mind

plink --bfile $FILE_1K.geno.mind --maf $MAF --make-bed --out $FILE_1K.geno.mind.maf


awk '{print$2}' "$FILE_QC.bim"> QCFILE_SNPs.txt
plink --bfile $FILE_1K.geno.mind.maf --extract QCFILE_SNPs.txt --make-bed --out $FILE_1K.geno.mind.maf.extract
awk '{print$2}' $FILE_1K.geno.mind.maf.extract.bim > 1kG_SNPs.txt
plink --bfile $FILE_QC --extract 1kG_SNPs.txt --make-bed --out $FILE_QC.extract

## Change build on 1000 Genomes data build to match build of FileSet data
awk '{print$2,$4}' $FILE_QC.extract.bim > buildmap.txt
plink --bfile $FILE_1K.geno.mind.maf.extract --update-map buildmap.txt --make-bed --out $FILE_1K.geno.mind.maf.extract.build

##Merge datasets

#1. set reference allele
awk '{print$2,$5}' $FILE_1K.geno.mind.maf.extract.build.bim > 1kg_ref-list.txt
plink --bfile $FILE_QC.extract --reference-allele 1kg_ref-list.txt --make-bed --out Map-adj

#2. resolve strand issues 
awk '{print$2,$5,$6}' $FILE_1K.geno.mind.maf.extract.build.bim > 1kGMDS_strand_tmp
awk '{print$2,$5,$6}' Map-adj.bim > Map-adj_tmp
sort 1kGMDS_strand_tmp Map-adj_tmp |uniq -u > all_differences.txt
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
plink --bfile Map-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_map

#3. Check for remaining uncorresponding SNPs
awk '{print$2,$5,$6}' corrected_map.bim > corrected_map_tmp
sort 1kGMDS_strand_tmp corrected_map_tmp |uniq -u  > uncorresponding_SNPs.txt
wc -l uncorresponding_SNPs.txt

awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt
plink --bfile corrected_map --exclude SNPs_for_exclusion.txt --make-bed --out $FILE_QC.extract.rem
plink --bfile $FILE_1K.geno.mind.maf.extract.build --exclude SNPs_for_exclusion.txt --make-bed --out $FILE_1K.geno.mind.maf.extract.build.rem

##Perform merge
plink --bfile $FILE_QC.extract.rem --bmerge $FILE_1K.geno.mind.maf.extract.build.rem --allow-no-sex --make-bed --out merged

##throw out ambiguous SNPs
wget https://github.com/eatkinson/Post-QC/raw/master/find_cg_at_snps.py
sed -i -e 's/file(bimfile)/open(bimfile)/g' find_cg_at_snps.py ##for compatibility in python3
sed -i -e 's/ line\[1\]/(line\[1\])/g' find_cg_at_snps.py ##for compatibility in python3
python find_cg_at_snps.py merged.bim > ATCGsites
plink --bfile merged --exclude ATCGsites --make-bed --out MDS_merge2

##Perform PCA 
##**Download the file with population information of the 1000 genomes dataset.**

wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/20100804.ALL.panel

##Convert population codes into superpopulation codes
awk '{print$1,$1,$2}' 20100804.ALL.panel > race_1kG.txt
awk '{print$1,$2,"OWN"}' $FILE_QC.extract.rem.fam > racefile_own.txt
sed 's/JPT/ASN/g' race_1kG.txt>race_1kG2.txt
sed 's/ASW/AFR/g' race_1kG2.txt>race_1kG3.txt
sed 's/CEU/EUR/g' race_1kG3.txt>race_1kG4.txt
sed 's/CHB/ASN/g' race_1kG4.txt>race_1kG5.txt
sed 's/CHD/ASN/g' race_1kG5.txt>race_1kG6.txt
sed 's/YRI/AFR/g' race_1kG6.txt>race_1kG7.txt
sed 's/LWK/AFR/g' race_1kG7.txt>race_1kG8.txt
sed 's/TSI/EUR/g' race_1kG8.txt>race_1kG9.txt
sed 's/MXL/AMR/g' race_1kG9.txt>race_1kG10.txt
sed 's/GBR/EUR/g' race_1kG10.txt>race_1kG11.txt
sed 's/FIN/EUR/g' race_1kG11.txt>race_1kG12.txt
sed 's/CHS/ASN/g' race_1kG12.txt>race_1kG13.txt
sed 's/PUR/AMR/g' race_1kG13.txt>race_1kG14.txt
cat race_1kG14.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt

#Perform PCA
plink --bfile MDS_merge2 --indep-pairwise 50 5 0.5 
plink --bfile MDS_merge2 --extract plink.prune.in  --make-bed --pca 10 'header' --out PCA_MERGE

Rscript PCAPlot.R
Rscript ScreePlot.R
