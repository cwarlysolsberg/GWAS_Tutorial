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


We assume that you have the following files (or you can download it from [here](https://drive.google.com/file/d/1x_G0Gxk9jFMY-PMqwtg6-vdEyUPp5p5u/view?usp=sharing)):

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
!!! warning
    While we do provide a rough guide on how to perform LDpred on bed files separated into individual chromosomes, this script is untested and extra caution is required

### 0. Prepare workspace
On some server, you might need to first use the following code in order to run LDpred with multi-thread

=== "prepare workspace and load bigsnpr"

    ```R
    library(bigsnpr)
    options(bigstatsr.check.parallel.blas = FALSE)
    options(default.nproc.blas = NULL)
    ```

### 1. Read in the phenotype and covariate files

=== "read in phenotype and covariates"

    ```R 
    library(data.table)
    library(magrittr)
    phenotype <- fread("EUR.height")
    covariate <- fread("EUR.cov")
    pcs <- fread("EUR.eigenvec")
    # rename columns
    colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
    # generate required table
    pheno <- merge(phenotype, covariate) %>%
        merge(., pcs)
    ``` 
### 2. Obtain HapMap3 SNPs
LDpred2 authors recommend restricting the analysis to only the HapMap3 SNPs

=== "load HapMap3 SNPs"

    ```R
    info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))
    ```
### 3. Load and transform the summary statistic file 

=== "Load summary statistic file"
    
    ``` R
    # Read in the summary statistic file
    sumstats <- bigreadr::fread2("Height.QC.gz") 
    # LDpred 2 require the header to follow the exact naming
    names(sumstats) <-
        c("chr",
        "pos",
        "rsid",
        "a1",
        "a0",
        "n_eff",
        "beta_se",
        "p",
        "OR",
        "INFO",
        "MAF")
    # Transform the OR into log(OR)
    sumstats$beta <- log(sumstats$OR)
    # Filter out hapmap SNPs
    sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
    ```

    !!! Warning
        
        Here, we know the exact ordering of the summary statistics file. 
        However, in many cases, the ordering of the summary statistics differ, 
        thus one must rename the columns according to their actual ordering

### 3. Calculate the LD matrix

=== "Genome Wide bed file"

    ```R
    # Get maximum amount of cores
    NCORES <- nb_cores()
    # Open a temporary file
    tmp <- tempfile(tmpdir = "tmp-data")
    on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
    # Initialize variables for storing the LD score and LD matrix
    corr <- NULL
    ld <- NULL
    # We want to know the ordering of samples in the bed file 
    fam.order <- NULL
    # preprocess the bed file (only need to do once for each data set)
    snp_readBed("EUR.QC.bed")
    # now attach the genotype object
    obj.bigSNP <- snp_attach("EUR.QC.rds")
    # extract the SNP information from the genotype
    map <- obj.bigSNP$map[-3]
    names(map) <- c("chr", "rsid", "pos", "a1", "a0")
    # perform SNP matching
    info_snp <- snp_match(sumstats, map)
    # Assign the genotype to a variable for easier downstream analysis
    genotype <- obj.bigSNP$genotypes
    # Rename the data structures
    CHR <- map$chr
    POS <- map$pos
    # get the CM information from 1000 Genome
    # will download the 1000G file to the current directory (".")
    POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
    # calculate LD
    for (chr in 1:22) {
        # Extract SNPs that are included in the chromosome
        ind.chr <- which(info_snp$chr == chr)
        ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
        # Calculate the LD
        corr0 <- snp_cor(
                genotype,
                ind.col = ind.chr2,
                ncores = NCORES,
                infos.pos = POS2[ind.chr2],
                size = 3 / 1000
            )
        if (chr == 1) {
            ld <- Matrix::colSums(corr0^2)
            corr <- as_SFBM(corr0, tmp)
        } else {
            ld <- c(ld, Matrix::colSums(corr0^2))
            corr$add_columns(corr0, nrow(corr))
        }
    }
    # We assume the fam order is the same across different chromosomes
    fam.order <- as.data.table(obj.bigSNP$fam)
    # Rename fam order
    setnames(fam.order,
            c("family.ID", "sample.ID"),
            c("FID", "IID"))
    ```

=== "Chromosome separated bed files"

    ```R
    # Get maximum amount of cores
    NCORES <- nb_cores()
    # Open a temporary file
    tmp <- tempfile(tmpdir = "tmp-data")
    on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
    # Initialize variables for storing the LD score and LD matrix
    corr <- NULL
    ld <- NULL
    # We want to know the ordering of samples in the bed file 
    info_snp <- NULL
    fam.order <- NULL
    for (chr in 1:22) {
        # preprocess the bed file (only need to do once for each data set)
        # Assuming the file naming is EUR_chr#.bed
        snp_readBed(paste0("EUR_chr",chr,".bed"))
        # now attach the genotype object
        obj.bigSNP <- snp_attach(paste0("EUR_chr",chr,".rds"))
        # extract the SNP information from the genotype
        map <- obj.bigSNP$map[-3]
        names(map) <- c("chr", "rsid", "pos", "a1", "a0")
        # perform SNP matching
        tmp_snp <- snp_match(sumstats[sumstats$chr==chr,], map)
        info_snp <- rbind(info_snp, tmp_snp)
        # Assign the genotype to a variable for easier downstream analysis
        genotype <- obj.bigSNP$genotypes
        # Rename the data structures
        CHR <- map$chr
        POS <- map$pos
        # get the CM information from 1000 Genome
        # will download the 1000G file to the current directory (".")
        POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
        # calculate LD
        # Extract SNPs that are included in the chromosome
        ind.chr <- which(tmp_snp$chr == chr)
        ind.chr2 <- tmp_snp$`_NUM_ID_`[ind.chr]
        # Calculate the LD
        corr0 <- snp_cor(
                genotype,
                ind.col = ind.chr2,
                ncores = NCORES,
                infos.pos = POS2[ind.chr2],
                size = 3 / 1000
            )
        if (chr == 1) {
            ld <- Matrix::colSums(corr0^2)
            corr <- as_SFBM(corr0, tmp)
        } else {
            ld <- c(ld, Matrix::colSums(corr0^2))
            corr$add_columns(corr0, nrow(corr))
        }
        # We assume the fam order is the same across different chromosomes
        if(is.null(fam.order)){
            fam.order <- as.data.table(obj.bigSNP$fam)
        }
    }
    # Rename fam order
    setnames(fam.order,
            c("family.ID", "sample.ID"),
            c("FID", "IID"))
    ```

### 4. Perform LD score regression

===  "Perform LD score regression"

    ```R
    df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
    ldsc <- snp_ldsc(   ld, 
                        length(ld), 
                        chi2 = (df_beta$beta / df_beta$beta_se)^2,
                        sample_size = df_beta$n_eff, 
                        blocks = NULL)
    h2_est <- ldsc[["h2"]]
    ```

### 5. Calculate the null R2

=== "Calculate the null R2 (quantitative trait)"
    
    ```R
    # Reformat the phenotype file such that y is of the same order as the 
    # sample ordering in the genotype file
    y <- pheno[fam.order, on = c("FID", "IID")]
    # Calculate the null R2
    # use glm for binary trait 
    # (will also need the fmsb package to calculate the pseudo R2)
    null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
        paste0("Height~Sex+", .) %>%
        as.formula %>%
        lm(., data = y) %>%
        summary
    null.r2 <- null.model$r.squared
    ```

=== "Calculate the null R2 (binary trait)"
    
    ```R
    library(fmsb)
    # Reformat the phenotype file such that y is of the same order as the 
    # sample ordering in the genotype file
    y <- pheno[fam.order, on = c("FID", "IID")]
    # Calculate the null R2
    # use glm for binary trait 
    # (will also need the fmsb package to calculate the pseudo R2)
    null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
        paste0("Height~Sex+", .) %>%
        as.formula %>%
        glm(., data = y, family=binomial) %>%
        summary
    null.r2 <- fmsb::NagelkerkeR2(null.model)
    ```

!!! important
    Scripts for binary trait analysis only serve as a reference as we have not simulate any binary traits. 
    In addition, Nagelkerke $R^2$ is biased when there are ascertainment of samples. For more information, please refer to [this paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.21614)

### 6. Obtain LDpred adjusted beta
=== "infinitesimal model"
    
    ```R
    beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
    ```

=== "grid model"

    ```R
    # Prepare data for grid model
    p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
    h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
    grid.param <-
        expand.grid(p = p_seq,
                h2 = h2_seq,
                sparse = c(FALSE, TRUE))
    # Get adjusted beta from grid model
    beta_grid <-
        snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)
    ```

=== "auto model"

    ```R
    # Get adjusted beta from the auto model
    multi_auto <- snp_ldpred2_auto(
        corr,
        df_beta,
        h2_init = h2_est,
        vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
        ncores = NCORES
    )
    beta_auto <- sapply(multi_auto, function(auto)
        auto$beta_est)
    ```

### 7. Obtain model PRS
#### Using Genome wide bed file

=== "infinitesimal model"

    ```R
    if(is.null(obj.bigSNP)){
        obj.bigSNP <- snp_attach("EUR.QC.rds")
    }
    genotype <- obj.bigSNP$genotypes
    # calculate PRS for all samples
    ind.test <- 1:nrow(genotype)
    pred_inf <- big_prodVec(    genotype,
                                beta_inf,
                                ind.row = ind.test,
                                ind.col = info_snp$`_NUM_ID_`)
    ```

=== "grid model"

    ```R
    if(is.null(obj.bigSNP)){
        obj.bigSNP <- snp_attach("EUR.QC.rds")
    }
    genotype <- obj.bigSNP$genotypes
    # calculate PRS for all samples
    ind.test <- 1:nrow(genotype)
    pred_grid <- big_prodMat(   genotype, 
                                beta_grid, 
                                ind.col = info_snp$`_NUM_ID_`)
    ```

=== "auto model"

    ```R
    if(is.null(obj.bigSNP)){
        obj.bigSNP <- snp_attach("EUR.QC.rds")
    }
    genotype <- obj.bigSNP$genotypes
    # calculate PRS for all samples
    ind.test <- 1:nrow(genotype)
    pred_auto <-
        big_prodMat(genotype,
                    beta_auto,
                    ind.row = ind.test,
                    ind.col = info_snp$`_NUM_ID_`)
    # scale the PRS generated from AUTO
    pred_scaled <- apply(pred_auto, 2, sd)
    final_beta_auto <-
        rowMeans(beta_auto[,
                    abs(pred_scaled -
                        median(pred_scaled)) <
                        3 * mad(pred_scaled)])
    pred_auto <-
        big_prodVec(genotype,
            final_beta_auto,
            ind.row = ind.test,
            ind.col = info_snp$`_NUM_ID_`)
    ```

#### Using chromosome separated bed files

=== "infinitesimal model"

    ```R
    pred_inf <- NULL
    for(chr in 1:22){
        obj.bigSNP <- snp_attach(paste0("EUR_chr",chr,".rds"))
        genotype <- obj.bigSNP$genotypes
        # calculate PRS for all samples
        ind.test <- 1:nrow(genotype)
        # Extract SNPs in this chromosome
        chr.idx <- which(info_snp$chr == chr)
        ind.chr <- info_snp$`_NUM_ID_`[chr.idx]
        tmp <- big_prodVec(genotype,
                                beta_inf[chr.idx],
                                ind.row = ind.test,
                                ind.col = ind.chr)
        if(is.null(pred_inf)){
            pred_inf <- tmp
        }else{
            pred_inf <- pred_inf + tmp
        }
    }
    ```

=== "grid model"

    ```R
    pred_grid <- NULL
    for(chr in 1:22){
        obj.bigSNP <- snp_attach(paste0("EUR_chr",chr,"_.rds"))
        genotype <- obj.bigSNP$genotypes
        # calculate PRS for all samples
        ind.test <- 1:nrow(genotype)
        # Extract SNPs in this chromosome
        chr.idx <- which(info_snp$chr == chr)
        ind.chr <- info_snp$`_NUM_ID_`[chr.idx]

        tmp <- big_prodMat( genotype, 
                            beta_grid[chr.idx], 
                            ind.col = ind.chr)

        if(is.null(pred_grid)){
            pred_grid <- tmp
        }else{
            pred_grid <- pred_grid + tmp
        }
    }
    ```

=== "auto model"

    ```R
    pred_auto <- NULL
    for(chr in 1:22){
        obj.bigSNP <- snp_attach(paste0("EUR_chr",chr,"_.rds"))
        genotype <- obj.bigSNP$genotypes
        # calculate PRS for all samples
        ind.test <- 1:nrow(genotype)
        # Extract SNPs in this chromosome
        chr.idx <- which(info_snp$chr == chr)
        ind.chr <- info_snp$`_NUM_ID_`[chr.idx]

        tmp <-
            big_prodMat(genotype,
                        beta_auto[chr.idx],
                        ind.row = ind.test,
                        ind.col = ind.chr)
        # scale the PRS generated from AUTO
        pred_scaled <- apply(tmp, 2, sd)
        final_beta_auto <-
            rowMeans(beta_auto[chr.idx,
                        abs(pred_scaled -
                            median(pred_scaled)) <
                            3 * mad(pred_scaled)])
        tmp <-
            big_prodVec(genotype,
                final_beta_auto,
                ind.row = ind.test,
                ind.col = ind.chr)
        if(is.null(pred_auto)){
            pred_auto <- tmp
        }else{
            pred_auto <- pred_auto + tmp
        }
    }
    ```

### 8. Get the final performance of the LDpred models

=== "infinitesimal model"

    ```R
    reg.formula <- paste("PC", 1:6, sep = "", collapse = "+") %>%
        paste0("Height~PRS+Sex+", .) %>%
        as.formula
    reg.dat <- y
    reg.dat$PRS <- pred_inf
    inf.model <- lm(reg.formula, dat=reg.dat) %>%
        summary
    (result <- data.table(
        infinitesimal = inf.model$r.squared - null.r2,
        null = null.r2
    ))
    ```

=== "grid model"

    ```R
    reg.formula <- paste("PC", 1:6, sep = "", collapse = "+") %>%
        paste0("Height~PRS+Sex+", .) %>%
        as.formula
    reg.dat <- y
    max.r2 <- 0
    for(i in 1:ncol(pred_grid)){
        reg.dat$PRS <- pred_grid[,i]
        grid.model <- lm(reg.formula, dat=reg.dat) %>%
            summary  
        if(max.r2 < grid.model$r.squared){
            max.r2 <- grid.model$r.squared
        }
    }
    (result <- data.table(
        grid = max.r2 - null.r2,
        null = null.r2
    ))
    ```

=== "auto model"

    ```R
    reg.formula <- paste("PC", 1:6, sep = "", collapse = "+") %>%
        paste0("Height~PRS+Sex+", .) %>%
        as.formula
    reg.dat <- y
    reg.dat$PRS <- pred_auto
    auto.model <- lm(reg.formula, dat=reg.dat) %>%
        summary
    (result <- data.table(
        auto = auto.model$r.squared - null.r2,
        null = null.r2
    ))
    ```

??? note "How much phenotypic variation does the PRS from each model explain?"
    Infinitesimal = 0.0100
    
    Grid Model = 0.00180

    Auto Model = 0.171
