options(bitmapType="cairo")

library(ggplot2)
png('scree.png',units="in",width=13,height=7,res=300)
val=read.table('PCA_MERGE.eigenval')
val$index=as.numeric(row.names(val)) 
ggplot(val, aes(x=index, y=V1)) +
   geom_point() +
  geom_line()+
  ggtitle("Screeplot of the first 10 PCs")+
  ylab('eigenvalue')+
  xlab('PCs')+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))+
  theme(legend.position = 'none')
dev.off()