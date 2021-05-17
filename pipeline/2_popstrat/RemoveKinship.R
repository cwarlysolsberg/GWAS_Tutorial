options(bitmapType="cairo")

data<- read.table(file="kin.txt",header=TRUE)
data_related=data[,c(1,2,3,4)]
missing <- read.table("plink.imiss", header =TRUE, as.is=T)
FID1=data_related[,c(1,2)]
FID2=data_related[,c(3,4)]
FID1$index=row.names(FID1)
FID2$index=row.names(FID2)
FID1$IID=FID1$IID1
FID2$IID=FID2$IID2
FID1[,1:2]=c()
FID2[,1:2]=c()
FID1=merge(FID1,missing,by="IID")
FID2=merge(FID2,missing,by="IID")
FID1[,4:6]=c()
FID2[,4:6]=c()
FID1$index=as.numeric(FID1$index)
FID2$index=as.numeric(FID2$index)
q=c(setdiff(FID2$index, FID1$index), setdiff(FID1$index, FID2$index))
if (length(q) != 0) {
FID1=FID1[!(FID1$index==q),]
FID2=FID2[!(FID2$index==q),]
}

FID1=FID1[order(FID1[,2]),]
FID2=FID2[order(FID2[,2]),]

##this binds the FIDs,IIDs,missingness for all pairs
bind=cbind(FID1,FID2)
bind$index=c()
bind$index=c()
##this displays just the values for missingness for each pair
bindval=cbind(FID1$F_MISS,FID2$F_MISS)
colnames(bindval)=c(1,2)
max=as.numeric(colnames(bindval)[apply(bindval,1,which.max)])

#finds the corresponding IID and FID (bind1 and bind2) for the individual with the higher missingness so it can be removed
bind1=bind[cbind(seq_along(max*3-1), max*3-1)]
bind2=bind[cbind(seq_along(max*3-2), max*3-2)]
final=cbind(bind1,bind2)
final=unique(final)
write.table(final, 'low_call_rate.txt', append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE,quote=FALSE)