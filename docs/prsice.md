# Background

PRSice-2 is one of the dedicated PRS programs which automates many of the steps from the previous page that used a sequence of PLINK functions (plus some QC steps). 
On this page you will run a PRS analysis using PRSice-2, which implements the standard C+T method.

## Obtaining PRSice-2
`PRSice-2` can be downloaded from:

| Operating System | Link |
| -----------------|:----:|
| Linux 64-bit | [v2.3.3](https://github.com/choishingwan/PRSice/releases/download/2.3.3/PRSice_linux.zip) |
| OS X 64-bit | [v2.3.3](https://github.com/choishingwan/PRSice/releases/download/2.3.3/PRSice_mac.zip) |

and can be directly used after extracting the file. 

In this tutorial, you will only need `PRSice.R` and `PRSice_XXX` where XXX is the operation system

## Required Data

This analysis assumes that you have the following files (or you can download it from [here](https://drive.google.com/file/d/1x_G0Gxk9jFMY-PMqwtg6-vdEyUPp5p5u/view?usp=sharing)): 

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post QC base data file. While PRSice-2 can automatically apply most filtering on the base file, it cannot remove duplicated SNPs|
|**EUR.QC.bed**| This file contains the genotype data that passed the QC steps |
|**EUR.QC.bim**| This file contains the list of SNPs that passed the QC steps |
|**EUR.QC.fam**| This file contains the samples that passed the QC steps |
|**EUR.height**| This file contains the phenotype data of the samples |
|**EUR.cov**| This file contains the covariates of the samples |
|**EUR.eigenvec**| This file contains the principal components (PCs) of the samples |

## Running PRS analysis
To run PRSice-2 we need a single covariate file, and therefore our covariate file and PCs file should be combined. This can be done with `R` as follows:

=== "without data.table"    

    ```R
    covariate <- read.table("EUR.cov", header=T)
    pcs <- read.table("EUR.eigenvec", header=F)
    colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
    cov <- merge(covariate, pcs, by=c("FID", "IID"))
    write.table(cov,"EUR.covariate", quote=F, row.names=F)
    q()
    ```

=== "with data.table"

    ```R
    library(data.table)
    covariate <- fread("EUR.cov")
    pcs <- fread("EUR.eigenvec", header=F)
    colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
    cov <- merge(covariate, pcs)
    fwrite(cov,"EUR.covariate", sep="\t")
    q()
    ```

which generates **EUR.cov**.

PRSice-2 can then be run to obtain the PRS results as follows:

=== "Linux"
    ```bash
    Rscript PRSice.R \
        --prsice PRSice_linux \
        --base Height.QC.gz \
        --target EUR.QC \
        --binary-target F \
        --pheno EUR.height \
        --cov EUR.covariate \
        --base-maf MAF:0.01 \
        --base-info INFO:0.8 \
        --stat OR \
        --or \
        --out EUR
    ```

=== "OS X"
    ```bash
    Rscript PRSice.R \
        --prsice PRSice_mac \
        --base Height.QC.gz \
        --target EUR.QC \
        --binary-target F \
        --pheno EUR.height \
        --cov EUR.covariate \
        --base-maf MAF:0.01 \
        --base-info INFO:0.8 \
        --stat OR \
        --or \
        --out EUR
    ```

=== "Windows"
    ```bash
    Rscript PRSice.R ^
        --prsice PRSice_win64.exe ^
        --base Height.QC.gz ^
        --target EUR.QC ^
        --binary-target F ^
        --pheno EUR.height ^
        --cov EUR.covariate ^
        --base-maf MAF:0.05 ^
        --base-info INFO:0.8 ^
        --stat OR ^
        --or ^
        --out EUR
    ```

The meaning of the parameters are as follow:

| Paramter | Value | Description|
|:-:|:-:|:-|
|prsice|PRSice_xxx| Informs `PRSice.R` that the location of the PRSice binary |
|base| Height.QC.gz| Informs `PRSice` that the name of the GWAS summary statistic |
|target| EUR.QC| Informs `PRSice` that the input genotype files should have a prefix of `EUR.QC` |
|binary-target| F| Indicate if the phenotype of interest is a binary trait. F for no |
|pheno| EUR.height| Provide `PRSice` with the phenotype file |
|cov| EUR.covariate| Provide `PRSice` with the covariate file |
|base-maf| MAF:0.05| Filter out SNPs with MAF < 0.05 in the GWAS summary statistics, using information in the `MAF` column|
|base-info| INFO:0.8| Filter out SNPs with INFO < 0.8 in the GWAS summary statistics, using information in the `INFO` column|
|stat| OR| Column name of the column containing the effect size|
|or|-| Inform `PRSice` that the effect size is an Odd Ratio|
|out | EUR | Informs `PRSice` that all output should have a prefix of `EUR`|

This will automatically perform "high-resolution scoring" and generate the "best-fit" PRS (in **EUR.best**), with associated plots of the results. 
Users should read Section 4.6 of our paper to learn more about issues relating to overfitting in PRS analyses.  

??? note "Which P-value threshold generates the "best-fit" PRS?"
    0.13995

??? note "How much phenotypic variation does the "best-fit" PRS explain?"
    0.166117
    
`lassosum` is one of the dedicated PRS programs which is an `R` package that uses penalised regression (LASSO) in its approach to PRS calculation.

## Installing lassosum

!!! note
    The script used here is based on lassosum version 0.4.4

!!! note
    For more details, please refer to [lassosum's homepage](https://github.com/tshmak/lassosum/)

You can install `lassosum` and its dependencies in `R` with the following command:

```R
install.packages(c("devtools","RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
library(devtools)
install_github("tshmak/lassosum")
```

## Required Data

Again, we assume that we have the following files (or you can download it from [here](https://drive.google.com/file/d/1x_G0Gxk9jFMY-PMqwtg6-vdEyUPp5p5u/view?usp=sharing)): 

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post-QCed summary statistic |
|**EUR.QC.bed**| The genotype file after performing some basic filtering |
|**EUR.QC.bim**| This file contains the SNPs that passed the basic filtering |
|**EUR.QC.fam**| This file contains the samples that passed the basic filtering |
|**EUR.height**| This file contains the phenotype of the samples |
|**EUR.cov**| This file contains the covariates of the samples |
|**EUR.eigenvec**| This file contains the PCs of the samples |

## Running PRS analysis

We can run lassosum as follows: 

``` R
library(lassosum)
# Prefer to work with data.table as it speeds up file reading
library(data.table)
library(methods)
library(magrittr)
# For multi-threading, you can use the parallel package and 
# invoke cl which is then passed to lassosum.pipeline
library(parallel)
# This will invoke 2 threads. 
cl <- makeCluster(2)

sum.stat <- "Height.QC.gz"
bfile <- "EUR.QC"
# Read in and process the covariates
covariate <- fread("EUR.cov")
pcs <- fread("EUR.eigenvec") %>%
    setnames(., colnames(.), c("FID","IID", paste0("PC",1:6)))
# Need as.data.frame here as lassosum doesn't handle data.table 
# covariates very well
cov <- merge(covariate, pcs)

# We will need the EUR.hg19 file provided by lassosum 
# which are LD regions defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
ld.file <- "EUR.hg19"
# output prefix
prefix <- "EUR"
# Read in the target phenotype file
target.pheno <- fread("EUR.height")[,c("FID", "IID", "Height")]
# Read in the summary statistics
ss <- fread(sum.stat)
# Remove P-value = 0, which causes problem in the transformation
ss <- ss[!P == 0]
# Transform the P-values into correlation
cor <- p2cor(p = ss$P,
        n = ss$N,
        sign = log(ss$OR)
        )
fam <- fread(paste0(bfile, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]


# Run the lassosum pipeline
# The cluster parameter is used for multi-threading
# You can ignore that if you do not wish to perform multi-threaded processing
out <- lassosum.pipeline(
    cor = cor,
    chr = ss$CHR,
    pos = ss$BP,
    A1 = ss$A1,
    A2 = ss$A2,
    ref.bfile = bfile,
    test.bfile = bfile,
    LDblocks = ld.file, 
    cluster=cl
)
# Store the R2 results
target.res <- validate(out, pheno = target.pheno, covar=cov)
# Get the maximum R2
r2 <- max(target.res$validation.table$value)^2
```


??? note "How much phenotypic variation does the "best-fit" PRS explain?"
    0.2395471
    
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
