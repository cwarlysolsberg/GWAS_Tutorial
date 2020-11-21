# Home
## Overview 
This tutorial provides a step-by-step guide to performing basic Genome Wide Association Study (GWAS) analyses and accompanies our [GWAS Guide paper](...MAYBE A LINK AT SOME POINT...). The aim of this tutorial is to provide a simple introduction of GWAS analyses while equipping existing users with a better understanding of the processes and implementation of mainstay GWAS tools. 

The tutorial is separated into five main section. 

1. [Quality Control of GWAS Data](QC.md)
2. [Population Stratification](popstrat.md)
3. [Association Analyses](assoc.md)
4. [Visualizing of Association results](visual.md)
5. [Deep dive into underpinnings of associations](assocdive.md)

!!! warning

    Data used in this tutorial are simulated and intended for demonstration purposes only. The results from this tutorial will not reflect the true performance of different software. 

!!! notess

    We assume you have basic knownledges on how to use the terminal, `plink` and `R`. 
    If you are unfamiliar with any of those, you can refer to the following online resources:

    | Software | link |
    |:-:|:-:|
    | terminal (OS X / Linux) |  [1](https://www.digitalocean.com/community/tutorials/basic-linux-navigation-and-file-management), [2](https://linuxconfig.org/bash-scripting-tutorial-for-beginners)
    | terminal (Windows)| [1](https://www.cs.princeton.edu/courses/archive/spr05/cos126/cmd-prompt.html), [2](https://www.tutorialspoint.com/batch_script/index.htm) |
    | plink | [v1.90](https://www.cog-genomics.org/plink/1.9/), [v1.75](http://zzz.bwh.harvard.edu/plink/) |
    | R | [1](https://www.tutorialspoint.com/r/index.htm)


!!! note

    This tutorial is written for Linux and OS X operating systems. 
    Windows users will need to change some commands accordingly.

!!! note
    Throughout the tutorial you will see tabs above some of the code:

    === "A"

        ```bash
        echo "Tab A"
        ```

    === "B"

        ```bash
        echo "Tab B"
        ```

    You can click on the tab to change to an alternative code (eg. to a different operation system)
## Datasets
1. [Base data](base.md)
2. [Target data](https://drive.google.com/file/d/1uhJR_3sn7RA8U5iYQbcmTp6vFdQiF4F2/view?usp=sharing): Simulated data based on the 1000 Genomes Project European samples

## Requirements
To follow the tutorial, you will need the following programs installed:

1. [R](https://www.r-project.org/) (**version 3.2.3+**)
2. [PLINK 1.9](https://www.cog-genomics.org/plink2)

## Citation
If you find this tutorial helpful for a publication, then please consider citing:

!!! important "Citation"

    Choi, S.W., Mak, T.S. & O’Reilly, P.F. Tutorial: a guide to performing polygenic risk score analyses. Nat Protoc (2020). [https://doi.org/10.1038/s41596-020-0353-1](https://doi.org/10.1038/s41596-020-0353-1)
