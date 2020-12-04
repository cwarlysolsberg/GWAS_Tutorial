# Home
## Overview 
This tutorial provides a step-by-step guide to performing basic Genome Wide Association Study (GWAS) analyses and accompanies our [GWAS Guide paper](...MAYBE A LINK AT SOME POINT...). The aim of this tutorial is to provide a simple introduction of GWAS analyses while equipping existing users with a better understanding of the processes and implementation of mainstay GWAS tools. 

The tutorial is separated into five main section. 

1. [Quality Control of GWAS Data](QC.md)
2. [Population Stratification](popstrat.md)
3. [Association Analyses](assoc.md)
4. [Visualizing of Association results](visual.md)
5. [Deep dive into underpinnings of associations](assocdive.md)
6. [Polygenic Risk Score](prs.md)

!!! warning

    Data used in this tutorial are simulated and intended for demonstration purposes only. The results from this tutorial will not reflect the true performance of different software. 

!!! notes

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


## Requirements
To follow the tutorial, you will need the following programs installed:

1. [R](https://www.r-project.org/) (**version 3.2.3+**)
2. [PLINK 1.9](https://www.cog-genomics.org/plink2)

## Citation
If you find this tutorial helpful for a publication, then please consider citing:

!!! important "Citation"



##https://medium.com/tenxor/how-to-generate-a-sequence-diagram-within-markdown-using-js-sequence-diagram-and-mkdocs-91dd4fe0b8fb when we want to make flowchart
