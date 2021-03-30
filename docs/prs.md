# Calculating and Analysing PRS


## Background
In this section of the tutorial you will learn how to construct a Polygenic Risk Score (PRS) for Alzheimer's disease, using PLINK and ldpred, and use this PRS to predict the likelihood of late-onset alzheimer's disease. As part of this analysis, you will estimate the heritability of alzheimer's disease using ld score regression on the GWAS summary statistics, and using GCTA on the target dataset. 

Whenever evaluating the predictiveness of a PRS, it is of vital importance that the target data set was not included in the original GWAS analysis. When using GWAS results from previously published work, it is important to check the accompanying article for the data sources that the authors used and ensure that your target data was not included in the study. When this is the case, it is advised to search for GWAS results elsewhere, or to contact the original GWAS authors to see whether it is possible to acquire meta-analysed summary statistics that exclude the target population.

The target data that we use in this tutorial is a dataset of 176 cases and 188 controls for late-onset alzheimer's disease ``https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2667989/pdf/main.pdf``. We have cleaned this target data according to our quality control protocol as outlined in tutorials [1](QC.md) and [2](popstrat.md), after which 341 individuals remain (167 cases and 174 controls). `provide download links to cleaned target files` 

The GWAS summary statistics that we will use are from the largest and most recently published GWAS for alzheimer's disease (Jansen et al., 2019). Click [here](https://ctg.cncr.nl/software/summary_statistics) to download the GWAS summary statistics estimated (AD_sumstats_Jansenetal_2019sept.txt.gz), which you will use as your weights for the polygenic score. Note that this is a GWAS, not on late-onset alzheimer's disease, but on a combination of alzheimer's disease and alzheimer's disease by proxy (i.e. through parental diagnoses). The construction of a PRS on a different, but related phenotype may reduce its overall predictiveness, but is still a useful endaveour as long as the phenotypes in the GWAS and target dataset are genetically related. Phenotypes are often coded differently in GWAS to maximize sample size and harmonize different datasets. One of the uses of a PRS is that they can be applied to study the relationships between genotype and phenotype in a smaller dataset in which phenotypes may be recorded in more detail.

`Links to install software here` (GCTA,PLINK,LDpred)

Before you start, make sure the clean target data file in plink binary format (.bed, .bim and .fam) and the GWAS summary statistics (.txt) are in your working directory. We define a TARGETSET and GWASSET variable at the start, to easily refer to our target and GWAS summary data, respectively.

```bash
TARGETSET=adgwas.qcout.clean
GWASSET=AD_sumstats_Jansenetal_2019sept.txt
```
###QC of GWAS summary data
First, let's look at the first 5 lines of the GWAS summary statistics

```bash
head -5 $GWASSET
uniqID.a1a2     CHR     BP      A1      A2      SNP     Z       P       Nsum    Neff    dir     EAF     BETA    SE
1:715265_T_C    1       715265  T       C       rs12184267      2.12197306477   0.03384 359856  359856  ??+?    0.0408069       0.0126426452822487      0.0059579668998386
1:715367_G_A    1       715367  G       A       rs12184277      1.95791489337   0.05024 360230  360230  ??+?    0.0410687       0.0116235111669331      0.00593667845639937
1:717485_A_C    1       717485  A       C       rs12184279      1.91243829548   0.05582 360257  360257  ??+?    0.0405759       0.0114189092819899      0.0059708641627697
1:720381_T_G    1       720381  T       G       rs116801199     2.29540369149   0.02171 360980  360980  ??+?    0.042162        0.0134428918549354      0.00585643906767846

```

The ReadMe of the GWAS summary data clarifies each of these columns. In this tutorial, we will only use the following columns: 
  - **SNP**: the rsID for each tested SNP
  - **BP**: the SNP's base position
  - **A1**: the reference allele
  - **A2**: the alternative allele
  - **Z**: The Z-score
  - **BETA**: the estimated coefficient of the SNP on the diagnosis for alzheimer's disease
  - **P**: the P-value associatied with BETA
  - **Nsum** The sample size for the particular SNP. 

We first define some variables such that we can easily to refer to the columns we need in our code:

```bash 
BPCOLNO=3
REFACOLNO=4
ALTACOLNO=5
SNPCOLNO=6
PVALCOLNO=8
BETACOLNO=13
ZCOLNO=7
NCOLNO=9
```

###Remove duplicate SNPs from GWAS summary statistics

Most GWAS software packages, such as PLINK, do not allow for duplicate SNPs in the GWAS summary data. We will obtain a list of rsIDs using awk, sort these SNPs, and extract the repeated lines using ``uniq -d'', saving these into a new file called duplicate.snp. We next use grep to filter out these duplicate SNPs from our GWAS summary statistics. 

```bash
awk -v c1=$SNPCOLNO '{print $c1}' $GWASSET |\
sort |\
uniq -d > duplicate.snp 
grep -vf duplicate.snp $GWASSET > $GWASSET.nodup
```


```bash
wc -l $GWASSET
13367300 AD_sumstats_Jansenetal_2019sept.txt

wc -l $GWASSET.nodup
13336963 AD_sumstats_Jansenetal_2019sept.txt.nodup
```

This filtering procedure leaves us with over 13 million SNPs that are non-duplicates. Duplicate SNPs occur, for example, due to coding mistakes, or when SNPs are multiallelic. 


###Further cleaning of GWAS summary statistics
When constructing a PRS, it is of importance to ensure that the SNP weights reported in the GWAS summary statistics are correctly matched with the target data. Any SNPs that are incorrectly matched will be assigned a wrong weight, and thus introduce noise into the PRS. Here, we remove strand ambiguous alleles, as well as indels (here marked by "I'/D"). Furthermore, we only keep columns that have rsIDs. 

```bash
awk '!( ($4=="A" && $5=="T") || ($4=="D" && $5=="I")  || \
        ($4=="T" && $5=="A") || ($4=="I" && $5=="D") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' $GWASSET.nodup > $GWASSET.nodup.noambu
        
awk 'NR==1' $GWASSET.nodup.noambu | awk -v c1=$SNPCOLNO '$c1 ~ /^rs/' $GWASSET.nodup.noambu > $GWASSET.nodup.noambu.rsid 
```

```bash
wc -l $GWASSET.nodup.noambu
10331009 AD_sumstats_Jansenetal_2019sept.txt.noambu

wc -l $GWASSET.nodup.noambu.rsid
10227986  AD_sumstats_Jansenetal_2019sept.txt.nodup.noambu.rsid
```

Notes on GWAS cleaning: irrelevant here: but don't include X, Y, or MT chromosome. If available, filter on info. 


###Matching cleaned GWAS summary statistics with target data. 
Our next goal is to take our target data, and only keep the SNPs that are present in the GWAS summary statistics. We first count how many SNPs there are in the GWAS summary and target data, respectively. 

```bash 
###
wc -l $GWASSET.nodup.noambu.rsid
10227986 AD_sumstats_Jansenetal_2019sept.txt.nodup.noambu.rsid

wc -l $TARGETSET.bim 
191069 adgwas.qcout.clean.bim

```

``Resolve strand issues, flip alleles``

Note that only 191069 SNPs are present in our GWAS summary data (``wc -l $TARGETSET.bim``), such that we cannot use any of the 10036917 remaining SNPs. We do not need to observe all SNPs from our GWAS summary statistics in our target data to construct a useful PRS. Even the missingness of causal SNPs for alzheimer's disease in the target data is not a big problem as long as we have still genotyped nearby SNPs that are in linkage disequilibrium with such causal SPNs. Nonetheless, the low SNP density in our target dataset will negatively effect the predictive power of our PRS. The GWAS summary results do have a wide SNP coverage, which is of more crucial importance, such that we can attach PRS weights to each SNP in our limited target data set.

Next, we use awk to get a list of SNP rsIDs and base pair positions from our GWAS summary statistics. We next take our target data, only keep the SNPs that are present in the GWAS summary statistics, and update their base pair positions using the update-map flag in plink. This updates the base pair positions from the build of the target data set (hg36), to the build of the GWAS summary data (hg37) ``check``

```bash
awk -v c1=$SNPCOLNO -v c2=$BPCOLNO '{print $c1,$c2}' $GWASSET.nodup.noambu.rsid > gwassnps.txt
awk -v c1=$SNPCOLNO 'FNR>1{print $c1}' $GWASSET.nodup.noambu.rsid > snplist.txt
plink --bfile $TARGETSET --extract snplist.txt --update-map gwassnps.txt --make-bed --out $TARGETSET.out
```

```bash
wc -l $TARGETSET.out.bim
155949 adgwas.qcout.clean.out.bim
```

155,949 SNPs present in our target data set are also present in our GWAS summary data.

Finally, we use awk  and grep to restrict our gwas summary statistics to the SNPs that are also present in our target data.
```bash
#Use only variables in GWAS summary stats that are also present in genome data. (but keep header of the original summary stats file)
awk '{print $2}' $TARGETSET.out.bim > $TARGETSET.SNPs.txt
head -1 $GWASSET.nodup > GWASanalysis.txt && grep -wFf $TARGETSET.SNPs.txt $GWASSET.nodup >> GWASanalysis.txt
wc -l GWASanalysis.txt
```


Finally, we flip reference alleles and resolve any strand issues.
```bash
awk -v c1=$SNPCOLNO -v c2=$REFACOLNO 'FNR>1{print $c1,$c2}' GWASanalysis.txt > gwasreflist.txt
plink --bfile $TARGETSET.out --reference-allele gwasreflist.txt --make-bed --out $TARGETSET.out.ref

awk -v c1=$SNPCOLNO -v c2=$REFACOLNO -v c3=$ALTACOLNO '{print$c1,$c2,$c3}' GWASanalysis.txt > gwasstrandlist.txt
awk '{print$2,$5,$6}' $TARGETSET.out.ref.bim > $TARGETSET.strandlist
sort gwasstrandlist.txt $TARGETSET.strandlist |uniq -u > all_differences.txt
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
plink --bfile $TARGETSET.out.ref --flip flip_list.txt --reference-allele gwasreflist.txt --make-bed --out $TARGETSET.out.ref.strand
```

Investigate problematic SNPs and throw them out:
```bash
awk '{print$2,$5,$6}' $TARGETSET.out.ref.strand.bim > corrected_map_tmp
sort gwasstrandlist.txt corrected_map_tmp |uniq -u  > uncorresponding_SNPs.txt
wc -l uncorresponding_SNPs.txt

awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt
plink --bfile $TARGETSET.out.ref.strand --exclude SNPs_for_exclusion.txt --make-bed --out $TARGETSET.def
```
####GCTA to estimate SNP-based heritability
Before constructing a PGS on the target data, it is useful to estimate the SNP-based heritability of the phenotype of interest using gcta. Due to measurement error in the PRS weights, the share of explained variance that the PRS can explain in a linear regression is lower than the overall SNP-based heritability. Therefore, estimation of SNP-based heritability gives us a reasonable expectation of the upper bound that we can expect from the predictiveness of the PRSs that we will construct. Genome-based restricted maximum likelihood (GREML) estiamtes the degree of variance explained in the phenotype that is explained by all SNPs in the target data. This is often referred to as SNP-based heritability. GCTA is a software package used to conduct such GREML analyses. It is similar to PLINK in the sense that it is operated through bash commands, using flags to guide the analysis of interest.  

Before using GCTA to estimate SNP-based heritability, it is advisable to conduct a power calculation using the associated power calculator: 
https://shiny.cnsgenomics.com/gctaPower/. To arrive at reasonable power (>80%), a sample size of 2000 individuals is roughly the minimum that is needed for most traits. The target data set here is notably smaller. With our 174 cases and 167 controls, and assuming a disease risk in the population of 0.1, a trait heritability of 0.5, alpha of 0.05, and the default variance of SNP-derived genetic relationships of 0.00002 gives us a dramatically low power of 0.0806. The standard error that we could expect for our analysis is 0.9760, which is incredibly large, given that a heritability estimate is bounded between 0 and 1 by definition. Hence, in our current data set, the estimation of SNP-based heritability using GREML is a useless endevaour. We nonetheless perform the GREML analysis anyway to illustrate how it is done in a dataset that is sufficiently large.

In the following lines of code, we first collect the family and individual IDs, and the phenotype data (case or control) into a new file (phenotype.phen). We next invoke gcta to estimate genome-wide relatedness matrix (using the --make-grm flag), which serves as an input to the heritability estimation. The estimation of the GRM is very sensitive to the inclusion of cryptic related individuals. We use the --grm-cutoff flag to throuw out individuals with relatedness value larger than 0.025. This is a more conservative parameter setting than the one included in our QC pipeline.

We next invoke the --reml flag in gcta64 to estimate the the SNP-based heritability. The necessary inputs are the GRM estimated in the line of code before, using the --grm flag, and the phenotype file using the --pheno flag.

```bash 
awk '{print $1,$2,$6}' $TARGETSET.def.fam > phenotype.phen

##To DO: visualize the phenotype here (write R-script)

##Estimate heritability using GCTA (GREML) (GCTA-GREML power calculator?) --> TO DO: Double-check we exclude close relatives. Check Yang et al. (2017). 
GRMCUT=0.025

gcta64 --bfile $TARGETSET.def --autosome --grm-cutoff $GRMCUT --make-grm --out $TARGETSET.grm
gcta64 --reml --grm $TARGETSET.grm --pheno phenotype.phen --out greml
```

The resulting output:

```
*******************************************************************
* Genome-wide Complex Trait Analysis (GCTA)
* version 1.93.2 beta Linux
* (C) 2010-present, Jian Yang, The University of Queensland
* Please report bugs to Jian Yang <jian.yang.qt@gmail.com>
*******************************************************************
Analysis started at 15:54:59 CEST on Mon Mar 29 2021.
Hostname: int1.bullx

Accepted options:
--reml
--grm adgwas.qcout.clean.grm
--pheno phenotype.phen
--out greml

Note: This is a multi-thread program. You could specify the number of threads by the --thread-num option to speed up the computation if there are multiple processors in your machine.

Reading IDs of the GRM from [adgwas.qcout.clean.grm.grm.id].
341 IDs read from [adgwas.qcout.clean.grm.grm.id].
Reading the GRM from [adgwas.qcout.clean.grm.grm.bin].
GRM for 341 individuals are included from [adgwas.qcout.clean.grm.grm.bin].
Reading phenotypes from [phenotype.phen].
Non-missing phenotypes of 341 individuals are included from [phenotype.phen].
Assuming a disease phenotype for a case-control study: 167 cases and 174 controls
Note: you can specify the disease prevalence by the option --prevalence so that GCTA can transform the variance explained to the underlying liability scale.
341 individuals are in common in these files.

Performing  REML analysis ... (Note: may take hours depending on sample size).
341 observations, 1 fixed effect(s), and 2 variance component(s)(including residual variance).
Calculating prior values of variance components by EM-REML ...
Updated prior values: 0.125147 0.124882
logL: 62.9277
Running AI-REML algorithm ...
Iter.   logL    V(G)    V(e)
1       62.93   0.15331 0.09642
2       63.00   0.21404 0.03527
3       63.06   0.21182 0.03765
4       63.06   0.21194 0.03754
5       63.06   0.21193 0.03754
Log-likelihood ratio converged.

Calculating the logLikelihood for the reduced model ...
(variance component 1 is dropped from the model)
Calculating prior values of variance components by EM-REML ...
Updated prior values: 0.25063
logL: 62.32649
Running AI-REML algorithm ...
Iter.   logL    V(e)
1       62.33   0.25063
Log-likelihood ratio converged.

Summary result of REML analysis:
Source  Variance        SE
V(G)    0.211931        0.176050
V(e)    0.037542        0.174494
Vp      0.249473        0.019152
V(G)/Vp 0.849515        0.699855

Sampling variance/covariance of the estimates of variance components:
3.099372e-02    -3.053760e-02
-3.053760e-02   3.044826e-02

Summary result of REML analysis has been saved in the file [greml.hsq].

Analysis finished at 15:54:59 CEST on Mon Mar 29 2021
Overall computational time: 0.09 sec.
```

The output shows four estimates of interest: V(G), the amount of variance in the phenotype that can be attributed to variance in the SNPs, V(e), the remaining variance that can be attributed to environmental factors, their sum Vp, and estimated SNP-based heritability: V(G) over Vp. The estimated SNP-based heritability is very high: 0.85. However, as expected, the standard error around this estimate (~0.7) is so large that even very low heritabilities can not be ruled out. In sum, the sample size is too small to derive any conclusions about SNP based heritability.

####Create a PRS using PLINK (clumpig and thresholding): 
`Use PRSIce or comment on different thresholds`
We are now ready to estimate our PRS. For each individual, we multiply their reference allele count at each SNP with a SNP weight estimated from the GWAS summary statistics. However, the GWAS coefficients as estimated in GWAS summary data are not corrected for linkage disequilibrium. Constructing a PRS without any correction for LD essentially leads to an overweighting of SNPs that are in dense LD-regions compared to SNPs that are not, resulting in lower predictability of the PRS. One method to deal with this is clumping. Clumping is a form of informed pruning: The R-squared between SNPs that reside within a given kb-window is computed, and one of the SNPs is thrown out of the R-squared is higher than a given threshold. The algorithm differs from pruning because it sorts all SNPs within a window increasingly by p-value, to ensure that SNPs with the lowest p-value are kept. 

The clumping algorithm gives a researcher substantial degrees of freedom when constructing a polygenic score. How to set the optimal parameters? Especially the choice of the p-value threshold has been shown to impact the predictiveness of PRSs, with stark differences in optimal threshold between different traits (Ware et al., 2017). Many researchers optimize over a grid of all possible parameter combinations. Automated tools, such as PRSICE2, are available for this. However, fitting many PRSs and choosing the best in terms of their evaluated predictive value in the target data set comes at the risk of overfitting the data. The smaller the dataset, the larger the risk of overfitting. We strongly suggest researchers to scan the literature first to see whether some consensus on optimal parameter values is available. For alzheimer's disease, many studies suggest that highly restrictive p-value thresholds (such that only the most significant SNPs for alzheimer's disease are included) result in the most predictive polygenic scores. Here, we follow Chaudhury et al. (2019), and set a p-value threshold of 0.000107, a window size of 250 kb, and an r-squared threshold of 0.1. 

```bash
plink --bfile $TARGETSET.def  --clump-p1 0.000107 --clump-p2 0.000107 --clump-r2 0.1 --clump-kb 250 --clump GWASanalysis.txt --clump-snp-field SNP \ 
--clump-field P --out $TARGETSET.clump
#Extract list of clumped SNPs (no header):
awk 'NR!=1{print $3}' $TARGETSET.clump.clumped >  $TARGETSET.clump.snp
```

```bash
wc -l $TARGETSET.clump.snp
76 adgwas.qcout.clean.clump.snp
```

Clumping our target data set using the p-values from gwas summary statistics leaves 76 SNPs in our dataset. We will extract these SNPs from our dataset, and construct a PRS using plink. The --score flag in PLINK next takes the reference allele count of each individual in our target data set, multiplies this value by the coefficient estimated in the GWAS, and sums this result over all included 76 SNPs.  

```bash
awk -v c1=$SNPCOLNO -v c2=$REFACOLNO -v c3=$BETACOLNO 'FNR>1{print $c1,$c2,$c3}' GWASanalysis.txt > score.txt
plink --bfile $TARGETSET.def --extract $TARGETSET.clump.snp --pheno phenotype.phen --score score.txt --out plink_score
```

The output is a .score file, which summarizes the Polygenic score for each individual in the variable ``SCORE'':

Before evaluating our PRSs by checking their predictability against the phenotype of interest, we estimate the first 10 principal components for each individual in our target data. These will be used as control variables. We also extract sex information from the .fam file as additional covariates in our regressions.

```bash
plink --bfile $TARGETSET.def --pca 10 header --out clusters
awk '{print $1,$2,$5}' $TARGETSET.def.fam > temp.txt
awk 'BEGIN{print "FID IID SEX"}1' temp.txt > sex.txt
```



We are now ready to load the score file in R, merge with the principal components, and evaluating the predictiveness of the PRSs. We standardize our PRSs such that they have mean 0 and standard deviation 1, and change the coding of our phenotype such that cases are equal to 1 and controls are equal to zero. We first estimate a null model that regresses the phenotype on the first ten PCs and sex, using both a linear probability model and a logit specification. We next estimate four models: a linear probability regressing the phenotype on the PRS only, one with all controls, and two similar specifications using a logit model. 

=== "Performed in R"

```{r}
library(stargazer)

data <- read.table(file="plink_score.profile",header=TRUE)

clusterdata <- read.table(file="clusters.eigenvec",header=TRUE)
sexdata <- read.table(file="sex.txt",header=TRUE)
mergedata <- merge(data,clusterdata)
mergedata <- merge(mergedata,sexdata)

#standardize data
mergedata$SCORE <- (mergedata$SCORE-mean(mergedata$SCORE))/(sd(mergedata$SCORE))
mergedata$PHENO <- mergedata$PHENO - 1

nullmodel<-lm(PHENO~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+SEX,data=mergedata)
nullr2<-summary(nullmodel)$r.squared

model<-lm(PHENO~SCORE,data=mergedata)
modelcontrols<-lm(PHENO~SCORE+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+SEX,data=mergedata)

model$IncreaseR2<-round(summary(model)$r.squared,digits=3)
modelcontrols$IncreaseR2<-round(summary(modelcontrols)$r.squared - nullr2,digits=3)

logitmodel<-glm(PHENO~SCORE,data=mergedata,family="binomial")
logitmodelcontrols<-glm(PHENO~SCORE+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+SEX,data=mergedata,family="binomial")

nullmodellogit<-glm(PHENO~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+SEX,data=mergedata,family="binomial")
nullpseudor2<-round(1-(nullmodellogit$deviance/nullmodellogit$null.deviance),digits=3)

logitmodel$PseudoR2<-round(1-(logitmodel$deviance/logitmodel$null.deviance),digits=3)
logitmodelcontrols$PseudoR2<-round(1-(logitmodelcontrols$deviance/logitmodelcontrols$null.deviance),digits=3)

logitmodel$IncreaseR2<-logitmodel$PseudoR2
logitmodelcontrols$IncreaseR2<-logitmodelcontrols$PseudoR2 - nullpseudor2

stargazer(model,modelcontrols,logitmodel,logitmodelcontrols,type="text",
          out="PLINKScore.txt",keep = c("SCORE"),
          add.lines = list(c("10 PCs", "No", "Yes", "No", "Yes"),c("Pseudo-R2","","",logitmodel$PseudoR2,logitmodelcontrols$PseudoR2),c("Increase-R2",model$IncreaseR2,modelcontrols$IncreaseR2,logitmodel$IncreaseR2,logitmodelcontrols$IncreaseR2)),
          star.cutoffs = c(0.05, 0.01, 0.001), float=FALSE)
```

##Create a PRS using LDPred 

```bash
python $PREDPATH coord --rs SNP --A1 A1 --A2 A2 --pos BP --eff_type LINREG --chr CHR --pval P --eff BETA --N $NGWAS --ssf GWASanalysis.txt --gf $TARGETSET.out --out pred.coord

## LDpred recommend radius to be Total number of SNPs in target / 3000 (CHECK!)
GSize="$(wc -l <"$TARGETSET.out.bim")"
ldrnum=$(( GSize / 3000 ))
python $PREDPATH gibbs --cf pred.coord --ldr $ldrnum --ldf pred.ld --out pred.weight --N $NGWAS

echo "python3 $PREDPATH score --gf $FILESET$out --rf test.weight --out test.score --pf phenotype.phen --pf-format LSTANDARD"
python $PREDPATH score --gf $TARGETSET.out --rf pred.weight --out test.score --pf phenotype.phen --pf-format LSTANDARD --pcs-file clusters.eigenvec
	
##P+T for reference:
python $PREDPATH p+t --cf pred.coord --ldr $ldrnum --out PTpred
	
python $PREDPATH score --gf $TARGETSET.out --rf PTpred --out PTscore.score --pf phenotype.phen --pf-format \
LSTANDARD --pcs-file clusters.eigenvec
```

=== "Performed in R"

```{r}
library(stargazer)

data <- read.table(file="test.score_LDpred-inf.txt",header=TRUE, sep = ",")

data$PRS <- (data$PRS-mean(data$PRS))/(sd(data$PRS))
data$true_phens <- data$true_phens - 1
#data$true_phens <- (data$true_phens-mean(data$true_phens))/(sd(data$true_phens))

model<-lm(true_phens~PRS,data=data)
logitmodel<-glm(true_phens~PRS,data=data,family="binomial")
modelcontrols<-lm(true_phens~PRS+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=data)
logitmodelcontrols<-glm(true_phens~PRS+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=data,family="binomial")

logitmodel$PseudoR2<-round(1-(logitmodel$deviance/logitmodel$null.deviance),digits=3)
logitmodelcontrols$PseudoR2<-round(1-(logitmodelcontrols$deviance/logitmodelcontrols$null.deviance),digits=3)
stargazer(model,modelcontrols,logitmodel,logitmodelcontrols,
          out="LDPREDstat.tex",keep = c("PRS"), 
          add.lines = list(c("10 PCs", "No", "Yes", "No", "Yes"),c("Pseudo-R2","","",logitmodel$PseudoR2,logitmodelcontrols$PseudoR2)), 
          star.cutoffs = c(0.05, 0.01, 0.001), float=FALSE)
```
##Notes on power: 
``Contrast with some results from a larger dataset (e.g. UKB?)``

##Notes on confouding, within-family analysis, etc.

##Notes on downward bias due to measurement error: Genetic IV?

##Further topics? Genetic correlations, Genomic SEM

##References
Chaudhury, S., Brookes, K. J., Patel, T., Fallows, A., Guetta-Baranes, T., Turton, J. C., ... & Thomas, A. J. (2019). Alzheimer's disease polygenic risk score as a predictor of conversion from mild-cognitive impairment. Translational psychiatry, 9(1), 1-7.
Ware, E. B., Schmitz, L. L., Faul, J., Gard, A., Mitchell, C., Smith, J. A., ... & Kardia, S. L. (2017). Heterogeneity in polygenic scores for common human traits. BioRxiv, 106062.
