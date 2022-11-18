dat <- read.table("pca_merged.eigenvec",h=F)
colnames(dat)[1] <- "FID" #change column header of first column to FID 
colnames(dat)[2] <- "IID" #change column header of second column to IID 
# add a column in the dat file indicating who is CEU, CHB, JPT, YRI; and who is from your own gwa dataset. This way we can colour the samples in the plot according to their etnicity/dataset. This information is present in the fam file of the merged dataset, in which the last column indicates whether the individual is a case (=='2'), control (=='1'), or which HapMap population (=='3','4','5',or '6'). 
pop <- read.table(file="racefile.txt",header=TRUE) 
colnames(pop)[3] <- "group"
# merge dat and pops
data <- merge(dat, pop, by=c("IID","FID"), all=FALSE)
names(data)
head(data)
dim(data)
str(data)
data$group <- as.factor(data$group)

OWN <- which(data$group=="OWN")
EUR <- which(data$group=="EUR") #note that this indicates that those individuals with '3' in phenotype column are the CEU (european population) individuals
ASN <- which(data$group=="ASN") 
AMR <- which(data$group=="AMR")
AFR <- which(data$group=="AFR")

pdf("pca-ancestry-plot2_corr.pdf")
plot(0,0,pch="",xlim=c(-0.1,0.05),ylim=c(-0.1,0.1),xlab="principal component 1", ylab="principal component 2")
points(data$V3[ASN],data$V4[ASN],pch=20,col="PURPLE")
points(data$V3[AMR],data$V4[AMR],pch=20,col="BLUE")
points(data$V3[AFR],data$V4[AFR],pch=20,col="GREEN")
points(data$V3[EUR],data$V4[EUR],pch=20,col="RED")
par(cex=0.5)
points(data$V3[OWN],data$V4[OWN],pch="+",col="BLACK")
abline(v=0.0075,col="gray32",lty=2)
abline(h=-0.005,col="gray32",lty=2)
legend("topright", pch=c(1,1,1,1,3),c("EUR","ASN","AMR","AFR","OWN"),col=c("RED","PURPLE","BLUE","GREEN","BLACK"),bty="o",cex=1)
dev.off()