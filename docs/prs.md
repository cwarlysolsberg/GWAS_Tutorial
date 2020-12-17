# Calculating and Analysing PRS


## Background
In this section of the tutorial you will learn how to construct a Polygenic Risk Score (PRS) for Alzheimer's disease, using PLINK and ldpred, and use this polygenic score to predict alzheimer's disease on a target data set that was **not included in the original GWAS** `Double-check this!`

Click here to download the [Target data](https://drive.google.com/drive/folders/1ud5F9WN9Xx3oXIkb5xIg1b_zz1nzp3IR)(adgwas.map.zip, adgwas.ped-1.zip, and samples.covar.zip). The target data contains allele frequencies in 364 individuals. We have run this target data through our quality control process as outlined in tutorials [1](QC.md) and [2](popstrat.md), after which 342 individuals remain (166 cases and 176 controls). `Describe cleaning process, provide links to cleaned target files` 

Click [here](https://ctg.cncr.nl/software/summary_statistics) to download GWAS summary statistics for alzheimer's disease as estimated by Iris Jansen et al., 2019 (AD_sumstats_Jansenetal_2019sept.txt.gz), which you will use as your weights for the polygenic score. 

`Links to install software here` (GCTA,PLINK,LDpred)

Before you start, make sure the clean target data file (.bed, .bim and .fam) and the GWAS summary statistics (.txt) are in your working directory. We define a TARGETSET and GWASSET variable at the start, to easily refer to our target and GWAS summary data, respectively.

```bash
TARGETSET=adgwas_geno_0.02_mind_0.05rem_autosomal_only_maf0.05_euro_hwe_control_1e-6_hwe_case_1e-10filterfounderlowcover
GWASSET=AD_sumstats_Jansenetal_2019sept.txt
```
###QC of GWAS summary data
First, let's look at the first 5 lines of the GWAS summary statistics

```bash
head -5 $GWASSET
```
                         
|uniqID.a1a2|CHR|BP|A1|A2|SNP|Z|P|Nsum|Neff|dir|EAF|BETA|SE|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|1:715265_T_C|1|715265|T|C|rs12184267|2.12197306477|0.03384|359856|359856|??+?|0.0408069|0.0126426452822487|0.0059579668998386|
|1:715367_G_A|1|715367|G|A|rs12184277|1.95791489337|0.05024|360230|360230|??+?|0.0410687|0.0116235111669331|0.00593667845639937|
|1:717485_A_C|1|717485|A|C|rs12184279|1.91243829548|0.05582|360257|360257|??+?|0.0405759|0.0114189092819899|0.0059708641627697|
|1:720381_T_G|1|720381|T|G|rs116801199|2.29540369149|0.02171|360980|360980|??+?|0.042162|0.0134428918549354|0.00585643906767846|

The ReadMe of the GWAS summary data clarifies each of these columns. In this tutorial, we will only use the following columns: 
  - **SNP**: the rsID for each tested SNP
  - **BP**: the SNP's base position
  - **A1**: the reference allele
  - **BETA**: the estimated coefficient of the SNP on the diagnosis for alzheimer's disease
  - **P**: the P-value associatied with BETA

We first define some variables such that we can easily to refer to the columns we need in our code:
```bash 
BPCOLNO=3
REFACOLNO=4
SNPCOLNO=6
PvalColNo=8
BETACOLNO=13
```

###Remove duplicate SNPs

Most GWAS software packages, such as PLINK, do not allow for duplicate SNPs in the GWAS summary data. We will obtain a list of the duplicate rsIDs using awk, and next use grep to filter out these duplicate SNPs from our GWAS summary statistics: 

```bash
awk -v c1=$SNPCOLNO '{print $c1}' $GWASSET |\
sort |\
uniq -d > duplicate.snp 
grep -vf duplicate.snp $GWASSET > $GWASSET.nodup
```


```bash
head duplicate.snp
```

10:1000836
10:1002602
10:1003598
10:1003623
10:1003624
10:1007067
10:1008676
10:1011884
10:1015405
10:1016749

`Explain why these are likely to be doubled`

###More elaborate QC of GWAS summary statistics file here
`chip-heritability of GWAS using LD-score regression` 


Next, we take our target data, and only keep the SNPs that are present in the GWAS summary statistics, using the update-map flag in plink. Note that only 190023 SNPs are present in our GWAS summary data, such that we cannot use any of the 13146940 remaining SNPs. We do not need to observe all SNPs from our GWAS summary statistics in our target data to construct a useful PRS, although missing SNPs will negatively effect the predictive power of our PRS. It is more important that our GWAS has wide SNP coverage `reference needed`. 190023 of the 195436 SNPs in our target data are present in the GWAS summary statistics, which is encouraging. 

```bash
###First count how many SNPs there are in the GWAS summary and target data, respectively
wc -l $GWASSET.nodup
wc -l $TARGETSET.bim

awk -v c1=$SNPCOLNO -v c2=$BPCOLNO '{print $c1,$c2}' $GWASSET.nodup > gwassnps.txt
out="out"
plink --bfile $TARGETSET --update-map gwassnps.txt --make-bed --out $TARGETSET.out
wc -l $TARGETSET.out.bim

#Use only variables in GWAS summary stats that are also present in genome data. (but keep header of the original summary stats file)
awk '{print $2}' $TARGETSET.out.bim > $TARGETSET.SNPs.txt
head -1 $GWASSET.nodup > GWASanalysis.txt && grep -wFf $TARGETSET.SNPs.txt $GWASSET.nodup >> GWASanalysis.txt
wc -l GWASanalysis.txt
```

`phenotype is already present in .fam file?`

####GCTA to estimate SNP-based heritability

```bash 
awk '{print $1,$2,$6}' $TARGETSET.out.fam > phenotype.phen

##To DO: visualize the phenotype here (write R-script)

##Estimate heritability using GCTA (GREML) (GCTA-GREML power calculator?) --> TO DO: Double-check we exclude close relatives. Check Yang et al. (2017). 
MAF=0.05
GRMCUT=0.025

gcta64 --bfile $TARGETSET.out --autosome --grm-cutoff $GRMCUT --maf $MAF --make-grm --out $TARGETSET.grm
gcta64 --reml --grm $TARGETSET.grm --pheno phenotype.phen --out greml
```

The resulting output:

```
Accepted options:
--reml
--grm adgwas_geno_0.02_mind_0.05rem_autosomal_only_maf0.05_euro_hwe_control_1e-6_hwe_case_1e-10filterfounderlow
cover.grm
--pheno phenotype.phen
--out greml

Note: This is a multi-thread program. You could specify the number of threads by the --thread-num option to spe
ed up the computation if there are multiple processors in your machine.

Reading IDs of the GRM from [adgwas_geno_0.02_mind_0.05rem_autosomal_only_maf0.05_euro_hwe_control_1e-6_hwe_cas
e_1e-10filterfounderlowcover.grm.grm.id].
342 IDs read from [adgwas_geno_0.02_mind_0.05rem_autosomal_only_maf0.05_euro_hwe_control_1e-6_hwe_case_1e-10fil
terfounderlowcover.grm.grm.id].
Reading the GRM from [adgwas_geno_0.02_mind_0.05rem_autosomal_only_maf0.05_euro_hwe_control_1e-6_hwe_case_1e-10
filterfounderlowcover.grm.grm.bin].
GRM for 342 individuals are included from [adgwas_geno_0.02_mind_0.05rem_autosomal_only_maf0.05_euro_hwe_contro
l_1e-6_hwe_case_1e-10filterfounderlowcover.grm.grm.bin].
Reading phenotypes from [phenotype.phen].
Non-missing phenotypes of 342 individuals are included from [phenotype.phen].
Assuming a disease phenotype for a case-control study: 166 cases and 176 controls
Note: you can specify the disease prevalence by the option --prevalence so that GCTA can transform the variance
 explained to the underlying liability scale.
342 individuals are in common in these files.

Performing  REML analysis ... (Note: may take hours depending on sample size).
342 observations, 1 fixed effect(s), and 2 variance component(s)(including residual variance).
Calculating prior values of variance components by EM-REML ...
Updated prior values: 0.125074 0.125129
logL: 62.7479
Running AI-REML algorithm ...
Iter.   logL    V(G)    V(e)
1       62.75   0.11893 0.13117
2       62.75   0.10559 0.14429
3       62.75   0.10531 0.14457
4       62.75   0.10530 0.14458
Log-likelihood ratio converged.

Calculating the logLikelihood for the reduced model ...
(variance component 1 is dropped from the model)
Calculating prior values of variance components by EM-REML ...
Updated prior values: 0.25052
logL: 62.59235
Running AI-REML algorithm ...
Iter.   logL    V(e)
1       62.59   0.25052
Log-likelihood ratio converged.

Summary result of REML analysis:
Source  Variance        SE
V(G)    0.105304        0.177626
V(e)    0.144584        0.177973
Vp      0.249888        0.019138
V(G)/Vp 0.421406        0.710403

Sampling variance/covariance of the estimates of variance components:
3.155100e-02    -3.142958e-02
-3.142958e-02   3.167444e-02

Summary result of REML analysis has been saved in the file [greml.hsq].

Analysis finished at 19:52:29 CET on Tue Dec 15 2020
Overall computational time: 0.11 sec.
```

####Create a PRS using PLINK: 
`Use PRSIce or comment on different thresholds`

```bash
awk -v c1=$SNPCOLNO -v c2=$REFACOLNO -v c3=$BETACOLNO 'FNR>1{print $c1,$c2,$c3}' GWASanalysis.txt > score.txt
plink --bfile $TARGETSET.out  --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 250 --clump GWASanalysis.txt --clump-snp-field SNP \
--clump-field P --out $TARGETSET.clump
#Extract list of clumped SNPs (no header):
awk 'NR!=1{print $3}' $TARGETSET.clump.clumped >  $TARGETSET.clump.snp

plink --bfile $TARGETSET.out --extract $TARGETSET.clump.snp --pheno phenotype.phen --allow-no-sex --score score.txt --out plink_score
#Calculate PCs as control (should be pruned first?, only use the 1990023 variants?)	
plink --bfile $TARGETSET.out --pca 10 header --out clusters
```

=== "Performed in R"

```{r}
library(stargazer)

data <- read.table(file="plink_score.profile",header=TRUE)
clusterdata <- read.table(file="clusters.eigenvec",header=TRUE)
mergedata <- merge(data,clusterdata)

#standardize data
mergedata$SCORE <- (mergedata$SCORE-mean(mergedata$SCORE))/(sd(mergedata$SCORE))
mergedata$PHENO <- mergedata$PHENO - 1
model<-lm(PHENO~SCORE,data=mergedata)
modelcontrols<-lm(PHENO~SCORE+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=mergedata)
logitmodel<-glm(PHENO~SCORE,data=mergedata,family="binomial")
logitmodelcontrols<-glm(PHENO~SCORE+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=mergedata,family="binomial")

logitmodel$PseudoR2<-round(1-(logitmodel$deviance/logitmodel$null.deviance),digits=3)
logitmodelcontrols$PseudoR2<-round(1-(logitmodelcontrols$deviance/logitmodelcontrols$null.deviance),digits=3)


stargazer(model,modelcontrols,logitmodel,logitmodelcontrols,
          out="PLINKScore.tex",keep = c("SCORE"),
          add.lines = list(c("10 PCs", "No", "Yes", "No", "Yes"),c("Pseudo-R2","","",logitmodel$PseudoR2,logitmodelcontrols$PseudoR2)),
          star.cutoffs = c(0.05, 0.01, 0.001), float=FALSE)
```
