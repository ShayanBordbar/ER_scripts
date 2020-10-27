#Normalization and Batch effect removal
Exp.kd.Raw.Dataset.list.Common
#z score normalizition of each gene, separately in each dataset
Exp.kd.geneZscore.Dataset.list.Common <- list()
for(i in 1:length(Exp.kd.Raw.Dataset.list.Common)){
  Exp.kd.geneZscore.Dataset.list.Common[[i]] <- t(scale(t(Exp.kd.Raw.Dataset.list.Common[[i]])))
}
#putting all experiments in one matrix
Exp.kd.geneZscore.Dataset.matrix.Common <- do.call(cbind, Exp.kd.geneZscore.Dataset.list.Common)
Exp.kd.geneZscore.Dataset.matrix.Common <- as.matrix(Exp.kd.geneZscore.Dataset.matrix.Common)

#look at the histogram of variances of genes across columns
a <- apply(Exp.kd.Raw.Dataset.matrix.Common, 1, var)

#removing rows with variance less than 1
Exp.kd.Raw.Dataset.matrix.Common.vgt1 <- Exp.kd.Raw.Dataset.matrix.Common[-c(which(a < 1)),]

#preparing for batch removal via Combat
aa <- unlist(lapply(Exp.kd.Raw.Dataset.list.Common, ncol))
a <- numeric(0)
for(i in 1:length(aa)){
  a <- c(a, rep(i,aa[i]))
}
aa <- data.frame(sample=c(1:ncol(Exp.kd.Raw.Dataset.matrix.Common)), batch=a)
a1 <- model.matrix(~1, data=aa)
#Batch removal
Exp.kd.Batch.Dataset.matrix.Common.vgt1 = ComBat(dat=Exp.kd.Raw.Dataset.matrix.Common.vgt1, batch=aa$batch, mod=a1)
#Checking batch removed
aaa <-  gPCA.batchdetect(t(Exp.kd.Batch.Dataset.matrix.Common.vgt1), batch = aa$batch, filt = NULL, nperm = 1000, center = FALSE, scaleY=FALSE, seed = NULL)
aaaa <- gPCA.batchdetect(t(Exp.kd.Raw.Dataset.matrix.Common.vgt1)  , batch = aa$batch, filt = NULL, nperm = 1000, center = FALSE, scaleY=FALSE, seed = NULL)

aaaaa <- unlist(lapply(Exp.kd.Raw.Dataset.list.Common, ncol))
aaaa <- numeric(0)
for(i in 1:length(aaaaa)){
  aaaa <- c(aaaa, rep(i,aaaaa[i]))
}
#plot pca before batch removal
aaa <- prcomp(t(Exp.kd.Raw.Dataset.matrix.Common.vgt1), scale = F)
aa <- aaa$x
par(mfrow = c(1, 2))
plot(aa[,1], aa[,2], pch = 19,  
     xlab="PC1",ylab="PC2",ylim = range(aaa$x[,2]), xlim = range(aaa$x[,1]), main= "before batch removal")
abline(h=0, v=0, lty = 2)
text(aa[,1], aa[,2], labels=aaaa,
     cex=0.7, pos = 3)
#plot pca after batch removal
a <- prcomp(t(Exp.kd.Batch.Dataset.matrix.Common.vgt1), scale = F)
aa <- a$x
plot(aa[,1], aa[,2], pch = 19,  
     xlab="PC1",ylab="PC2",ylim = range(aaa$x[,2]), xlim = range(aaa$x[,1]), main= "after batch removal")
abline(h=0, v=0, lty = 2)
text(aa[,1], aa[,2], labels=aaaa,
     cex=0.7, pos = 3)
####3D plots
plot3d(aa[,1], aa[,2], aa[,3])
scatter3D(aa[,1], aa[,2], aa[,3], colvar = aa[,3], col = NULL, add = FALSE)
text3D(aa[,1], aa[,2], aa[,3], labels = aaaa, colvar = NULL, add = FALSE)





