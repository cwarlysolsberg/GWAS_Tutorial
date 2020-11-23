# Population Stratifacation
## Obtaining the 1000 genome reference data set

We will be using data from the 1000 Genomes Project for the population stratification step.

You can download using the following command 

```bash
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz
```
The data will need to be converted into plink format with unique indentifiers created for each the SNPs with a missing rs-identifier:

!!! note
    You will need `plink` in this section, which can be download from [here](https://www.cog-genomics.org/plink/1.9/).

    Install the program `plink` and include its location in your PATH directory, which allows us to use `plink` instead of `./plink` in the commands below. If PLINK is not in your PATH directory and is instead in your working directory, replace all instances of `plink` in the tutorial with `./plink`.

```bash 
plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --make-bed --out 1000genomes.genotypes
plink --bfile 1000genomes.genotypes --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out 1000genomes_nomissing.genotypes
```

## Define Variables to be used in this section

**input files and variables**
```bash 
FILE_1K=1000genomes_nomissing.genotypes #1000 genomes file
FILE_QC=qcout #file from qc tab
FILE_PRUNEIN=plink.prune.in #snps in approximate linkage equilibrium 
GENO=0.02 #snp missingness filter
INDV=0.02 #individual missingness filter
MAF=0.05 #minor allele frequency filter
HWE_CONTROL=1e-6 #hardy weinburg equilibrium filter
```

**define general file tags**

```bash
sep="_"
tagbed=".bed"
tagbim=".bim"
tagfam=".fam"
tagmap=".map"
```

**tags for filtering 1000genome data**
```bash
TAG_1kG="1KG"
TAG_GENO="geno"
OUT1="$TAG_1kG$sep$TAG_GENO$sep$GENO"
TAG_MIND="mind"
OUT2="$OUT1$sep$TAG_MIND$sep$INDV"
TAG_MAF="maf"
OUT3=$OUT2$sep$TAG_MAF$sep$MAF
TAG_extract="extract"
OUT4=$OUT3$sep$TAG_extract
TAG_BUILD="samebuild"
OUT5=$OUT4$sep$TAG_BUILD
TAG_REM="removeproblem"
OUT6=$OUT5$sep$TAG_REM
```

**tags for filtering QCed data**
```bash
TAG_QC="QCin"
OUT1_QC=$TAG_QC$sep$TAG_extract
OUT2_QC=$OUT1_QC$sep$TAG_REM
tag_euro="euro"
FILEQCEURO=$FILE_QC$sep$tag_euro
```
##QC on 1000 Genomes data.

**Remove variants based on missing genotype data.**
```bash
plink --bfile $FILE_1K --geno $GENO --make-bed --out $OUT1
```
**Remove individuals based on missing genotype data**
```bash
plink --bfile $OUT1 --mind $INDV --allow-no-sex --make-bed --out $OUT2
```
**Remove variants based on MAF**
```bash
plink --bfile $OUT2 --maf $MAF --make-bed --out $OUT3
```
**Extract the variants present in our dataset from the 1000 genomes dataset**
```bash
awk '{print$2}' "$FILE_QC$tagbim"> QCFILE_SNPs.txt
awk '{print$2}' "$OUT3$tagbim"> 1kG_temp.bim
plink --bfile $OUT3 --extract QCFILE_SNPs.txt --make-bed --recode --out $OUT4
```

## Extract the variants present in 1000 Genomes dataset from your dataset.
```bash
awk '{print$2}' $OUT4$tagbim > 1kG_SNPs.txt
plink --bfile $FILE_QC --extract 1kG_SNPs.txt --recode --make-bed --out $OUT1_QC
```
*The datasets now contain the exact same variants.*

## Change build on 1000 Genomes data build to match build of HapMap data

!!! note
    Look at [Liftover tutorial](liftover.md) to see how to move data set to another build. 

```bash 
awk '{print$2,$4}' $OUT1_QC$tagmap > buildmap.txt
# buildmap.txt contains one SNP-id and physical position per line.
plink --bfile $OUT4 --update-map buildmap.txt --make-bed --out $OUT5
```

## Merge the Map and 1000 Genomes data sets

??? note "Prior to merging 1000 Genomes data with the data we want to make sure that the files are mergeable, for this we conduct 3 steps:"

    1) Make sure the reference genome is similar in your data and the 1000 Genomes Project datasets 
    
    2) Resolve strand issues. 
    
    3) Remove the SNPs which after the previous two steps still differ between datasets


**1) set reference genome**
```bash
awk '{print$2,$5}' $OUT5$tagbim > 1kg_ref-list.txt
plink --bfile $OUT1_QC --reference-allele 1kg_ref-list.txt --make-bed --out Map-adj
# The 1kG_MDS6 and the HapMap-adj have the same reference genome for all SNPs.
```

**2) Resolve strand issues**
```bash
awk '{print$2,$5,$6}' $OUT5$tagbim > 1kGMDS_strand_tmp
awk '{print$2,$5,$6}' Map-adj.bim > Map-adj_tmp
sort 1kGMDS_strand_tmp Map-adj_tmp |uniq -u > all_differences.txt
```

*Flip SNPs for resolving strand issues*
```bash
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
plink --bfile Map-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_map
```

*Check for SNPs which are still problematic after they have been flipped.*
```bash
awk '{print$2,$5,$6}' corrected_map.bim > corrected_map_tmp
sort 1kGMDS_strand_tmp corrected_map_tmp |uniq -u  > uncorresponding_SNPs.txt
```

**3) Remove problematic SNPs from your data and from the 1000 Genomes.**
```bash
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt
plink --bfile corrected_map --exclude SNPs_for_exclusion.txt --make-bed --out $OUT2_QC
plink --bfile $OUT5 --exclude SNPs_for_exclusion.txt --make-bed --out $OUT6
```

## Merge outdata with 1000 Genomes Data

```bash
plink --bfile $OUT2_QC --bmerge $OUT6$tagbed $OUT6$tagbim $OUT6$tagfam --allow-no-sex --make-bed --out MDS_merge2
```

## Perform MDS on Map-CEU data anchored by 1000 Genomes data.
# Using a set of pruned SNPs

plink --bfile MDS_merge2 --extract $FILE_PRUNEIN --genome --out MDS_merge2

plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2


### MDS-plot

# Download the file with population information of the 1000 genomes dataset.
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/20100804.ALL.panel

# The file 20130502.ALL.panel contains population codes of the individuals of 1000 genomes.

# Convert population codes into superpopulation codes (i.e., AFR,AMR,ASN, and EUR).
awk '{print$1,$1,$2}' 20100804.ALL.panel > race_1kG.txt
# Create a racefile of your own data.
awk '{print$1,$2,"-"}' $OUT2_QC$tagfam > racefile_own.txt


# Concatenate racefiles.
cat racefile_own.txt race_1kG.txt| sed -e '1i\FID IID race' > MDS_merge2.pop
sed -i -e "1d" MDS_merge2.pop
cut -d " " -f 3- MDS_merge2.pop >temp.txt
#make popfile for admixture script
mv temp.txt MDS_merge2.pop
# conda install -c bioconda admixture
###########run admixture script 
qsub -cwd -pe smp 8 -l mem_free=32G -l scratch=100G -l h_rt=40:20:00 ad.sh
admixture --supervised ./MDS_merge2.bed 12 > log_merge_admixture.out
Rscript --no-save ../R/admixtureplot.R
##european.txt file comes from admixture script
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

# Create a racefile of your own data.

awk '{print$1,$2,"OWN"}' $OUT2_QC$tagfam > racefile_own.txt

# Concatenate racefiles.
cat race_1kG14.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt

# Generate population stratification plot.
Rscript ../R/MDS_merged.R

## Exclude ethnic outliers.
# Select individuals in your own data below cut-off thresholds. The cut-off levels are not fixed ithresholds but have to be determined based on the visualization of the first two dimensions. To exclude ethnic outliers, the thresholds need to be set around the cluster of population of interest.

#awk '{ if ($4 <-0.04 && $5 >0.03) print $1,$2 }' MDS_merge2.mds > EUR_MDS_merge2
##below are the filters used for my own HA data
#awk '{ if ($4 < 0.1 && $4 > -0.05 && $5 > -0.1 && $5 < 0.03) print $1,$2 }' MDS_merge2.mds > EUR_MDS_merge2
# tag_euro="euro"
# FILEQCEURO=$FILE_QC$sep$tag_euro
# plink --bfile $FILE_QC --keep EUR_MDS_merge2 --make-bed --out $FILEQCEURO

## Exclude ethnic outliers.
plink --bfile $FILE_QC --keep europeans.txt --make-bed --out $FILEQCEURO
#plink --bfile $FILE_QC --keep EUR_MDS_merge2 --make-bed --out $FILEQCEURO
#HWE Limitations
TAG_HWE_CONTROL="hwe_control"
OUTQCEURO=$FILEQCEURO$sep$TAG_HWE_CONTROL$sep$HWE_CONTROL
echo "plink --bfile $FILEQCEURO --hwe $HWE_CONTROL --make-bed --out $OUTQCEURO"
plink --bfile $FILEQCEURO --hwe $HWE_CONTROL --make-bed --out $OUTQCEURO
# In this tutorial, we aim to remove all 'relatedness' from our dataset.
# To demonstrate that the majority of the relatedness was due to parent-offspring we only include founders (individuals without parents in the dataset).
found="founder"
plink --bfile $OUTQCEURO --filter-founders --make-bed --out $OUTQCEURO$sep$found

# #Check for cryptic relatedness
##install king and run
plink2 --bfile $OUTQCEURO$sep$found --make-king-table --king-table-filter 0.0884
sed 's/^#//' plink2.kin0 > kin.txt

plink --bfile $OUTQCEURO$sep$found --missing
Rscript ../R/Relatedness_cpw.R
lowcover="unrelated"
OUTCOVER=$lowcover
plink --bfile $OUTQCEURO$sep$found --remove 0.2_low_call_rate.txt --make-bed --out $OUTCOVER

## Create covariates based on MDS.
# Perform an MDS ONLY on qccase data without ethnic outliers. The values of the 10 MDS dimensions are subsequently used as covariates in the association analysis in the third tutorial.
	
plink --bfile $OUTCOVER --extract plink.prune.in --genome --out $OUTCOVER

tag_mds="MDS"
POPSTRATOUT=$OUTCOVER$sep$tag_mds
taggenome=".genome"
echo "plink --bfile $OUTCOVER --read-genome $OUTCOVER$taggenome --cluster --mds-plot 10 --out $POPSTRATOUT"
plink --bfile $OUTCOVER --read-genome $OUTCOVER$taggenome --cluster --mds-plot 10 --out $POPSTRATOUT

tag_dot_mds=".mds"

# Change the format of the .mds file into a plink covariate file.
awk '{print$1, $2, $4, $5, $6, $7,$8,$9,$10,$11,$12,$13}' $POPSTRATOUT$tag_dot_mds > covar_mds.txt

# The values in covar_mds.txt will be used as covariates, to adjust for remaining population stratification, in the third tutorial where we will perform a genome-wide association analysis.

mv $OUTCOVER$tagbed ./popstratout.bed
mv $OUTCOVER$tagbim ./popstratout.bim
mv $OUTCOVER$tagfam ./popstratout.fam
cp popstratout* ../3_Association_GWAS
cp covar_mds.txt ../3_Association_GWAS
cd ../3_Association_GWAS


## \# Sample size
We recommend that users only perform PRS analyses on target data of at least 100 individuals. The sample size of our target data here is 503 individuals. 

## \# File transfer
Usually we do not need to download and transfer the target data file because it is typically generated locally. However, the file should contain an md5sum code in case we send the data file to collaborators who may want to confirm that the file has not changed during the transfer.

??? note "What is the md5sum code for each of the target files?"

    |File|md5sum|
    |:-:|:-:|
    |**EUR.bed**           |98bcef133f683b1272d3ea5f97742e0e|
    |**EUR.bim**           |6b286002904880055a9c94e01522f059|
    |**EUR.cov**           |85ed18288c708e095385418517e9c3bd|
    |**EUR.fam**           |e7b856f0c7bcaffc8405926d08386e97|
    |**EUR.height**        |dd445ce969a81cded20da5c88b82d4df|

## \# Genome build
As stated in the base data section, the genome build for our base and target data is the same, as it should be.

## \# Standard GWAS QC
The target data must be quality controlled to at least the standards 
implemented in GWAS studies, e.g. removing SNPs with low genotyping rate, 
low minor allele frequency, out of Hardy-Weinberg Equilibrium, removing
individuals with low genotyping rate 
(see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).

The following `plink` command applies some of these QC metrics to the target data:


```bash
plink \
    --bfile EUR \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out EUR.QC
```

Each of the parameters corresponds to the following

| Paramter | Value | Description|
|:-:|:-:|:-|
| bfile | EUR | Informs `plink` that the input genotype files should have a prefix of `EUR` |
| maf | 0.01 | Removes all SNPs with minor allele frequency less than 0.05. Genotyping errors typically have a larger influence on SNPs with low MAF. Studies with large sample sizes could apply a lower MAF threshold|
| hwe | 1e-6 | Removes SNPs with low P-value from the Hardy-Weinberg Equilibrium Fisher's exact or chi-squared test. SNPs with significant P-values from the HWE test are more likely affected by genotyping error or the effects of natural selection. Filtering should be performed on the control samples to avoid filtering SNPs that are causal (under selection in cases). When phenotype information is included, plink will automatically perform the filtering in the controls. |
| geno | 0.01 | Excludes SNPs that are missing in a high fraction of subjects. A two-stage filtering process is usually performed (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).|
| mind | 0.01 | Excludes individuals who have a high rate of genotype missingness, since this may indicate problems in the DNA sample or processing. (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/) for more details).|
| make-just-fam | - | Informs `plink` to only generate the QC'ed sample name to avoid generating the .bed file.  |
| write-snplist | - | Informs `plink` to only generate the QC'ed SNP list to avoid generating the .bed file. |
| out | EUR.QC | Informs `plink` that all output should have a prefix of `EUR.QC` |

??? note "How many SNPs and samples were filtered?"
    - `14` samples were removed due to a high rate of genotype missingness
    - `5,353` SNP were removed due to missing genotype data
    -  `944` SNPs were removed due to being out of Hardy-Weinberg Equilibrium
    - `5,061` SNPs were removed due to low minor allele frequency

!!! note
    Normally, we can generate a new genotype file using the new sample list.
    However,  this will use up a lot of storage space. Using `plink`'s
    `--extract`, `--exclude`, `--keep`, `--remove`, `--make-just-fam` and `--write-snplist` functions, we can work 
    solely on the list of samples and SNPs without duplicating the 
    genotype file, reducing the storage space usage.  
    
Very high or low heterozygosity rates in individuals could be due to DNA contamination or to high levels of inbreeding. Therefore, samples with extreme heterozygosity are typically removed prior to downstream analyses. 

First, we perform prunning to remove highly correlated SNPs:

```bash
plink \
    --bfile EUR \
    --keep EUR.QC.fam \
    --extract EUR.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out EUR.QC
```

Each of the parameters corresponds to the following

| Paramter | Value | Description|
|:-:|:-:|:-|
| bfile | EUR | Informs `plink` that the input genotype files should have a prefix of `EUR` |
| keep | EUR.QC.fam | Informs `plink` that we only want to use samples in `EUR.QC.fam` in the analysis |
| extract | EUR.QC.snplist | Informs `plink` that we only want to use SNPs in `EUR.QC.snplist` in the analysis |
|indep-pairwise| 200 50 0.25 | Informs `plink` that we wish to perform prunning with a window size of 200 variants, sliding across the genome with step size of 50 variants at a time, and filter out any SNPs with LD $r^2$ higher than 0.25|
| out | EUR.QC | Informs `plink` that all output should have a prefix of `EUR.QC` |


This will generate two files 1) **EUR.QC.prune.in** and 2) **EUR.QC.prune.out**. All SNPs within **EUR.QC.prune.in** have a pairwise $r^2 < 0.25$. 


Heterozygosity rates can then be computed using `plink`:

```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.fam \
    --het \
    --out EUR.QC
```

This will generate the **EUR.QC.het** file, which contains F coefficient estimates for assessing heterozygosity.
We will remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean, which can be performed using the following `R` command (assuming that you have R downloaded, then you can open an `R` session by typing `R` in your terminal):

=== "Without library"
    ```R
    dat <- read.table("EUR.QC.het", header=T) # Read in the EUR.het file, specify it has header
    m <- mean(dat$F) # Calculate the mean  
    s <- sd(dat$F) # Calculate the SD
    valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
    write.table(valid[,c(1,2)], "EUR.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
    q() # exit R
    ```

=== "With data.table"
    ```R
    library(data.table)
    # Read in file
    dat <- fread("EUR.QC.het")
    # Get samples with F coefficient within 3 SD of the population mean
    valid <- dat[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)] 
    # print FID and IID for valid samples
    fwrite(valid[,c("FID","IID")], "EUR.valid.sample", sep="\t") 
    q() # exit R
    ```

??? note "How many samples were excluded due to high heterozygosity rate?"
    - `2` samples were excluded

## \# Ambiguous SNPs
These were removed during the base data QC.

## \# Mismatching SNPs
SNPs that have mismatching alleles reported in the base and target data may be resolvable by strand-flipping the alleles to their complementary alleles in e.g. the target data, such as for a SNP with A/C in the base data and G/T in the target. This can be achieved with the following steps:

1\. Load the bim file, the summary statistic and the QC SNP list into R

=== "Without data.table"
    ```R
    # Read in bim file
    bim <- read.table("EUR.bim")
    colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
    # Read in QCed SNPs
    qc <- read.table("EUR.QC.snplist", header = F, stringsAsFactors = F)
    # Read in the GWAS data
    height <-
        read.table(gzfile("Height.QC.gz"),
                header = T,
                stringsAsFactors = F, 
                sep="\t")
    # Change all alleles to upper case for easy comparison
    height$A1 <- toupper(height$A1)
    height$A2 <- toupper(height$A2)
    bim$B.A1 <- toupper(bim$B.A1)
    bim$B.A2 <- toupper(bim$B.A2)
    ```

=== "With data.table and magrittr"

    ```R
    # magrittr allow us to do piping, which help to reduce the 
    # amount of intermediate data types
    library(data.table)
    library(magrittr)
    # Read in bim file 
    bim <- fread("EUR.bim") %>%
        # Note: . represents the output from previous step
        # The syntax here means, setnames of the data read from
        # the bim file, and replace the original column names by 
        # the new names
        setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
        # And immediately change the alleles to upper cases
        .[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]
    # Read in summary statistic data (require data.table v1.12.0+)
    height <- fread("Height.QC.gz") %>%
        # And immediately change the alleles to upper cases
        .[,c("A1","A2"):=list(toupper(A1), toupper(A2))]
    # Read in QCed SNPs
    qc <- fread("EUR.QC.snplist", header=F)
    ```


2\. Identify SNPs that require strand flipping 

=== "Without data.table"

    ```R
    # Merge summary statistic with target
    info <- merge(bim, height, by = c("SNP", "CHR", "BP"))
    # Filter QCed SNPs
    info <- info[info$SNP %in% qc$V1,]
    # Function for finding the complementary allele
    complement <- function(x) {
        switch (
            x,
            "A" = "T",
            "C" = "G",
            "T" = "A",
            "G" = "C",
            return(NA)
        )
    }
    # Get SNPs that have the same alleles across base and target
    info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
    # Identify SNPs that are complementary between base and target
    info$C.A1 <- sapply(info$B.A1, complement)
    info$C.A2 <- sapply(info$B.A2, complement)
    info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
    # Update the complementary alleles in the bim file
    # This allow us to match the allele in subsequent analysis
    complement.snps <- bim$SNP %in% info.complement$SNP
    bim[complement.snps,]$B.A1 <-
        sapply(bim[complement.snps,]$B.A1, complement)
    bim[complement.snps,]$B.A2 <-
        sapply(bim[complement.snps,]$B.A2, complement)
    ```

=== "With data.table and magrittr"
    ```R
    # Merge summary statistic with target
    info <- merge(bim, height, by=c("SNP", "CHR", "BP")) %>%
        # And filter out QCed SNPs
        .[SNP %in% qc[,V1]]

    # Function for calculating the complementary allele
    complement <- function(x){
        switch (x,
            "A" = "T",
            "C" = "G",
            "T" = "A",
            "G" = "C",
            return(NA)
        )
    } 
    # Get SNPs that have the same alleles across base and target
    info.match <- info[A1 == B.A1 & A2 == B.A2, SNP]
    # Identify SNPs that are complementary between base and target
    com.snps <- info[sapply(B.A1, complement) == A1 &
                        sapply(B.A2, complement) == A2, SNP]
    # Now update the bim file
    bim[SNP %in% com.snps, c("B.A1", "B.A2") :=
            list(sapply(B.A1, complement),
                sapply(B.A2, complement))]
    ```


3\. Identify SNPs that require recoding in the target (to ensure the coding allele in the target data is the effective allele in the base summary statistic)

=== "Without data.table"

    ```R
    # identify SNPs that need recoding
    info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
    # Update the recode SNPs
    recode.snps <- bim$SNP %in% info.recode$SNP
    tmp <- bim[recode.snps,]$B.A1
    bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
    bim[recode.snps,]$B.A2 <- tmp

    # identify SNPs that need recoding & complement
    info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
    # Update the recode + strand flip SNPs
    com.snps <- bim$SNP %in% info.crecode$SNP
    tmp <- bim[com.snps,]$B.A1
    bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
    bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))

    # Output updated bim file
    write.table(
        bim,
        "EUR.QC.adj.bim",
        quote = F,
        row.names = F,
        col.names = F,
        sep="\t"
    )
    ```

=== "With data.table  and magrittr"

    ```R
    # identify SNPs that need recoding
    recode.snps <- info[B.A1==A2 & B.A2==A1, SNP]
    # Update the bim file
    bim[SNP %in% recode.snps, c("B.A1", "B.A2") :=
            list(B.A2, B.A1)]

    # identify SNPs that need recoding & complement
    com.recode <- info[sapply(B.A1, complement) == A2 &
                        sapply(B.A2, complement) == A1, SNP]
    # Now update the bim file
    bim[SNP %in% com.recode, c("B.A1", "B.A2") :=
            list(sapply(B.A2, complement),
                sapply(B.A1, complement))]
    # Write the updated bim file
    fwrite(bim, "EUR.QC.adj.bim", col.names=F, sep="\t")
    ```


4\. Identify SNPs that have different allele in base and target (usually due to difference in genome build or Indel)

=== "Without data.table"

    ```R
    mismatch <-
        bim$SNP[!(bim$SNP %in% info.match$SNP |
                    bim$SNP %in% info.complement$SNP | 
                    bim$SNP %in% info.recode$SNP |
                    bim$SNP %in% info.crecode$SNP)]
    write.table(
        mismatch,
        "EUR.mismatch",
        quote = F,
        row.names = F,
        col.names = F
    )
    q() # exit R
    ```

=== "With data.table"

    ``` R
    mismatch <- bim[!(SNP %in% info.match |
                        SNP %in% com.snps |
                        SNP %in% recode.snps |
                        SNP %in% com.recode), SNP]
    write.table(mismatch, "EUR.mismatch", quote=F, row.names=F, col.names=F)
    q() # exit R
    ```

5\. Replace **EUR.bim** with **EUR.QC.adj.bim**:

```bash
# Make a back up
mv EUR.bim EUR.bim.bk
ln -s EUR.QC.adj.bim EUR.bim
```
The above commands do the following:

1. Rename **EUR.bim** to **EUR.bim.bk**
2. Soft linking (`ln -s`) **EUR.QC.adj.bim** as **EUR.bim**

!!! note
    Most PRS software will perform strand-flipping automatically, thus this step is usually not required.

## \# Duplicate SNPs
Make sure to remove any duplicate SNPs in your target data (these target data were simulated and so include no duplicated SNPs).

## \# Sex chromosomes 
Sometimes sample mislabelling can occur, which may lead to invalid results. One indication of a mislabelled sample is a difference between reported sex and that indicated by the sex chromosomes. While this may be due to a difference in sex and gender identity, it could also reflect mislabeling of samples or misreporting and, thus, individuals in which there is a mismatch between biological and reported sex are typically removed. A sex check can be performed in PLINK, in which individuals are called as females if their X chromosome homozygosity estimate (F statistic) is < 0.2 and as males if the estimate is > 0.8.

Before performing a sex check, pruning should be performed (see [here](target.md#35-standard-gwas-qc)).
A sex check can then easily be conducted using `plink`
```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.valid.sample \
    --check-sex \
    --out EUR.QC
```

This will generate a file called **EUR.QC.sexcheck** containing the F-statistics for each individual. Individuals are typically called as being biologically male if the F-statistic is > 0.8 and biologically female if F < 0.2.

=== "Without library"

    ```R
    # Read in file
    valid <- read.table("EUR.valid.sample", header=T)
    dat <- read.table("EUR.QC.sexcheck", header=T)
    valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
    write.table(valid[,c("FID", "IID")], "EUR.QC.valid", row.names=F, col.names=F, sep="\t", quote=F) 
    q() # exit R
    ```

=== "With data.table"

    ```R
    library(data.table)
    # Read in file
    valid <- fread("EUR.valid.sample")
    dat <- fread("EUR.QC.sexcheck")[FID%in%valid$FID]
    fwrite(dat[STATUS=="OK",c("FID","IID")], "EUR.QC.valid", sep="\t") 
    q() # exit R
    ```

??? note "How many samples were excluded due mismatched Sex information?"
    - `4` samples were excluded

## \# Sample overlap
Since the target data were simulated there are no overlapping samples between the base and target data here (see the relevant section of [the paper](https://www.nature.com/articles/s41596-020-0353-1) for discussion of the importance of avoiding sample overlap). 

## \# Relatedness
Closely related individuals in the target data may lead to overfitted results, limiting the generalisability of the results. 

Before calculating the relatedness, pruning should be performed (see [here](target.md#35-standard-gwas-qc)).
Individuals that have a first or second degree relative in the sample ($\hat{\pi} > 0.125$) can be removed with the following command:

```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.valid \
    --rel-cutoff 0.125 \
    --out EUR.QC
```

??? note "How many related samples were excluded?"
    - `0` samples were excluded

!!! note
    A greedy algorithm is used to remove closely related individuals in a way that optimizes the size of the sample retained.                However, the algorithm is dependent on the random seed used, which can generate different results. Therefore, to reproduce
    the same result, you will need to specify the same random seed. 
    
    PLINK's algorithm for removing related individuals does not account for the phenotype under study. 
    To minimize the removal of cases of a disease, the following algorithm can be used instead: 
    [GreedyRelated](https://gitlab.com/choishingwan/GreedyRelated).

## Generate final QC'ed target data file
After performing the full analysis, you can generate a QC'ed data set with the following command:
```bash
plink \
    --bfile EUR \
    --make-bed \
    --keep EUR.QC.rel.id \
    --out EUR.QC \
    --extract EUR.QC.snplist \
    --exclude EUR.mismatch
```
