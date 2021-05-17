#PCAPlot.R

options(bitmapType="cairo")

library(dplyr)
data<- read.table(file="PCA_MERGE.eigenvec",header=TRUE)
eigenval<-read.table(file="PCA_MERGE.eigenval",header=F)
race<- read.table(file="racefile.txt",header=TRUE)
datafile<- merge(data,race,by=c("IID","FID"))
datafile=datafile %>% filter(datafile$race %in% c("EUR","ASN","AMR","AFR","OWN"))
head(datafile)
pve = eigenval/sum(eigenval)*100
library(ggplot2)
library(gridExtra)
png("pcagg.png",units="in",width=13,height=7,res=300)
	b <- ggplot(datafile, aes(PC1, PC2, col = race,group=race)) 
	b <- b+ geom_point(size = 2)
	b <-b + xlab(paste0("PC1 (", signif(pve[1,1],3), "%)")) + ylab(paste0("PC2 (", signif(pve[2,1], 3), "%)"))+
    theme(text = element_text(size=13),legend.text=element_text(size=13),legend.title = element_blank(),legend.position='top')+guides(color = guide_legend(override.aes = list(size=8)))
    c <-ggplot(datafile, aes(PC1, PC2, col = race,group=race)) 
    c <-c+ggplot2::stat_ellipse(
      geom = "path",
      position = "identity",
      show.legend = NA,
      size=2,
      inherit.aes = TRUE,
      type = "t",
      level = 0.95,
      segments = 51
    )
    c <- c + xlab(paste0("PC1 (", signif(pve[1,1],3), "%)")) + ylab(paste0("PC2 (", signif(pve[2,1], 3), "%)"))+
    theme(text = element_text(size=13),legend.title = element_blank(),legend.position = "none")
    grid.arrange(b, c, ncol=2)
dev.off()