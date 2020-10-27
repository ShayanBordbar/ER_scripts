#Estrogen receptor expression data
setwd("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor")
library(stringi)
library(stringr)
ERexp32hrs <- read.delim("Expression data/E-TABM-742-processed-data-1778086214.txt",header = T)
#colnames(ERexp32hrs)
ERexp32hrsNote <- read.delim(file = "Expression data/E-TABM-742.txt",header = T)
# which(as.character( ERexp32hrsNote$Hybridization.Name) =="1842761041_C")
# library(stringr)
# v <- apply(X = colnames(ERexp32hrs),FUN = )
# substr(x = colnames(ERexp32hrs)[2],start = 2, stop = nchar(colnames(ERexp32hrs)[2])-1) 
# nchar(colnames(ERexp32hrs)[2])
unu <- list()
for(i in 1:(nrow(ERexp32hrsNote)-1)){
 unu[[i]]<-which(!is.na(str_match(string = colnames(ERexp32hrs), pattern = as.character(ERexp32hrsNote$Hybridization.Name[i]) ))) 
  }
    
#grep(as.character(ERexp32hrsNote$Hybridization.Name[i]), colnames(ERexp32hrs))
#getting the name of each sample
sampNames <- character()
for(i in 1:53){
  sampNames[unu[[i]]] <- as.character(ERexp32hrsNote$Hybridization.Name[i])
}
#getting the cell type of each sample
celltype <- character()
for (i in 1:53){
  celltype[unu[[i]]] <- as.character(ERexp32hrsNote$Characteristics..CellLine.[i])
}
#getting the condition of each sample
condName <- character()
for (i in 1:53){
  condName[unu[[i]]] <- as.character(ERexp32hrsNote$X.html..head...head..body..pre.style.word.wrap..break.word..white.space..pre.wrap..Source.Name[i])
}
#removing the conditions of ZR-75 cell line
ExpMCF7 <- ERexp32hrs[,-c(which(celltype == "ZR-75.1"))]
#ordering the MCF7 conditions based on time

##getting the time of each condition
condTime <- character()
for (i in 2:length(condName)){
  condTime[i] <- substr(condName[i], nchar(condName[i])-2, nchar(condName[i])-1)
}
condTimeNum <- as.numeric(condTime)
for (i in 1:length(condTime)){
  if (is.na(condTimeNum[i])){
    condTimeNum[i] <- as.numeric(substr(condName[i], nchar(condName[i])-1, nchar(condName[i])-1))
  }
}
#getting the time of just MCF7 cell conditions
condTimeMCF7 <- condTimeNum[-c(which(celltype == "ZR-75.1"))]
condTimeMCF7 <- condTimeMCF7[2:length(condTimeMCF7)]
#sorting the MCF7 expression,names,conds  based on time
rownames(ExpMCF7) <- as.character(ExpMCF7[,1])
ExpMCF7 = ExpMCF7[,2:ncol(ExpMCF7)]
ExpMCF7sorted <- ExpMCF7[,sort(condTimeMCF7,decreasing = F,index.return = T)$ix]

sampNamesMCF7 <- sampNames[-c(which(celltype == "ZR-75.1"))]
sampNamesMCF7 <- sampNamesMCF7[2:length(sampNamesMCF7)]
sampNamesMCF7sorted <- sampNamesMCF7[sort(condTimeMCF7,decreasing = F,index.return = T)$ix]

condNameMCF7 <- condName[-c(which(celltype == "ZR-75.1"))]
condNameMCF7 <- condNameMCF7[2:length(condNameMCF7)]
condNameMCF7sorted <- condNameMCF7[sort(condTimeMCF7,decreasing = F,index.return = T)$ix]
########
ExpMCF7sortedNum <- as.matrix(ExpMCF7sorted)
for (i in 1:ncol(ExpMCF7sortedNum)){
  ExpMCF7sortedNum[,i] <- as.numeric(ExpMCF7sortedNum[,i])
}
#ExpMCF7sortedNum <- as.numeric(ExpMCF7sortedNum)
#range(ExpMCF7sorted[1,])

#ExpMCF7sortedNum[10,20]
#ExpMCF7sorted[10,20]
############
#using the bioconductor package lumi
library(lumi)
x.lumi <- lumiR.batch("Expression data/E-TABM-742-processed-data-1778086214.txt")
summary( x.lumi,'QC')
x.lumi.log2.qnorm <- lumiExpresso(x.lumi,varianceStabilize.param=list(method='log2'))
dataMatrix <- exprs(x.lumi.log2.qnorm )
#To speed up the processing and reduce false positives, remove the unex- pressed genes
presentCount <- detectionCall(x.lumi,Th = 0.05)
selDataMatrix <- dataMatrix[presentCount > 0,]
probeList <- rownames(selDataMatrix)
#getting the name of the condition of each BeadStudio
beadStcorName <- character()
BSN <- condName[2:213]
for(i in 1:53){
  beadStcorName[i] <- BSN[(i-1)*4 +1]
}
#str_match(string =beadStcorName, pattern = "MCF-7" )
#selecting conditions for MCF7
selDataMatrixMCF7 <- selDataMatrix[,!is.na(str_match(string =beadStcorName, pattern = "MCF-7" ))]
#time of each MCF7 exp
condTimeMCF7new <- numeric()
for (i in 1:26){
  condTimeMCF7new[i] <- condTimeMCF7[(i-1)*4 +1]
}
selDataMatrixMCF7sorted <- selDataMatrixMCF7[,sort(condTimeMCF7new,decreasing = F,index.return = T)$ix]
#get the range of expression for each gene
expressionrange <- apply(selDataMatrixMCF7sorted, 1, range)
hist(diff(expressionrange)) 
sum(diff(expressionrange) > 1) #1638 genes

library(ArrayExpress)
#trying to get the differentially expressed genes:
sortedconds <- sort(condTimeMCF7new,decreasing = F)
#try a zscore normalized version as well
#selDataMatrixMCF7sortedNormalized <- scale(selDataMatrixMCF7sorted)
design <- model.matrix(~factor(sortedconds))
fit <- lmFit(selDataMatrixMCF7sorted, design)
fitnormal <- lmFit(scale(selDataMatrixMCF7sorted), design)
ebayes <- eBayes(fit)
ebayesNormal <- eBayes(fitnormal)
#get the genes differentially expressed with adjusted p-values below 0.01. differential expression
#is caculated for each time point as a different factor
tab <- topTable(ebayes, adjust="bonferroni", p.value = 0.01, n = 20000)
tabnormal <- topTable(ebayesNormal, adjust="bonferroni", p.value = 0.000000001, n = 20000)
selDataMatrixMCF7sortedCopy <- selDataMatrixMCF7sorted
colnames(selDataMatrixMCF7sortedCopy) <- sortedconds
heatmap(selDataMatrixMCF7sortedCopy[rownames(tab),])
heatmap(selDataMatrixMCF7sortedCopyNorm[rownames(tabnormal),])

#get the genes differentially expressed with adjusted p-values below 0.01. differential expression
#is caculated for treated vs untreated samples
#design2 <- model.matrix(~factor( c(rep(0,4),rep(1,22))))
#fit2 <- lmFit(selDataMatrixMCF7sorted, design2)
#ebayes2 <- eBayes(fit2)
#tab2 <- topTable(ebayes2, adjust="bonferroni", p.value = 0.01, n = 20000)
#heatmap(selDataMatrixMCF7sortedCopy[rownames(tab2),])
##this led to 69 genes as ooposed to the previous method which led to 4357 genes

probeset <- as.character(rownames(tab)[1:3])

#annotating the differentially expressed genes
# source("http://www.bioconductor.org/getBioC.R")
# getBioC("limma")
# getBioC("affy")
# getBioC("hgu95av2")
# getBioC("estrogen") 
# getBioC("hgu95av2cdf")
# getBioC("simpleaffy")
# getBioC("annotate")
# getBioC("XML")
# library(affy)
# library(limma)
# library(simpleaffy)

source("https://bioconductor.org/biocLite.R")
biocLite("lumiHumanAll.db")
library("lumiHumanAll.db")
biocLite("illuminaHumanv4.db")
library(illuminaHumanv4.db)
library(annotate)
probeset <- as.character(rownames(tab))
probesetNormal <- as.character(rownames(tabnormal))
#getting gene symbols
DiffGeneSymb <- lookUp(probeset, "illuminaHumanv4.db", "SYMBOL")
#
DiffGeneEntrez <- lookUp(probeset, "illuminaHumanv4.db", "ENTREZID")
DiffGeneEntrezNormal <- lookUp(probesetNormal, "illuminaHumanv4.db", "ENTREZID")




#
##
####
########        Reading and Analyzing Chip-seq data
################
################################
################################################################
################################################################################################
#first dataset to be analyzed:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35109
#Hierarchical modularity in ERα transcriptional network is associated with distinct functions and implicates clinical outcomes. Sci Rep 2012
source("http://www.bioconductor.org/biocLite.R")
biocLite("BayesPeak")
library(BayesPeak)
biocLite("ChIPpeakAnno")
library(ChIPpeakAnno)
biocLite("ChIPseeker")
library(ChIPseeker)
biocLite("rGADEM")
library(rGADEM)
#reading the 0 e2 treatment chip data
ERchip0 <- read.bed(filename = "Chip-Data/GSM862178_t0_export.txt_unique.EE_W100_P0.957.bed")
ERchip1 <- read.bed(filename = "Chip-Data/GSM862179_t1_export.txt_unique.EE_W100_P0.951.bed")
ERchip4 <- read.bed(filename = "Chip-Data/GSM862180_t4_export.txt_unique.EE_W100_P0.96.bed")
ERchip24 <- read.bed(filename = "Chip-Data/GSM862181_t24_export.txt_unique.EE_W100_P0.951.bed")

#take a look at width of the peaks
width<-unlist(end(ranges(ERchip0))-start(ranges(ERchip0)))
hist(width, col=3)
#NCBI annotation
data(TSS.human.NCBI36)
TSS.human.NCBI36

#annotatedPeak1 = annotatePeakInBatch(ERchip1, AnnotationData=TSS.human.NCBI36)
#temp1 = as.data.frame(annotatedPeak1)
values(annotatedPeak)
pie(table(temp1$insideFeature))
table(temp1$insideFeature)

x <- org.Hs.egENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
mapped_genes <- mappedkeys(x[as.character(DiffGeneEntrezmod)])
# Convert to a list
xx <- as.list(x[mapped_genes])
hist(temp$distancetoFeature,xlim = c(-20000,20000),breaks = 2000)
sum(temp1$distancetoFeature < 20000)
#map these ensemble ids for the annotated peaks and compare them with the ones obtained by the chipseeker package.




#findOverlappingPeaks
ov1 <- findOverlappingPeaks(toGRanges(ERchip0), toGRanges(ERchip1), NameOfPeaks1="ER0hr", NameOfPeaks2="ER1hr")
summary(ov1)                            
#findOverlapsOfPeaks
summary(ov1$OverlappingPeaks)
table((ov1$OverlappingPeaks)$overlapFeature)
pie(table((ov1$OverlappingPeaks)$overlapFeature))
##############
#annotation using chipseek
#peak4hr <- readPeakFile("Chip-Data/GSM862180_t4_export.txt_unique.EE_W100_P0.96.bed")
peak4hr
covplot(peak4hr, weightCol="V5")
#associate peaks with genes
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb = txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak4hr, weightCol=NULL, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequen")

####
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
orgHsegdb <- org.Hs.eg.db
peakAnno <- annotatePeak(peak4hr, tssRegion=c(-3000, 3000),level = "transcript", TxDb = txdb) 
peakannotation4hr <- as.data.frame(peakAnno)
plotAnnoBar(peakAnno)                         
plotDistToTSS(peakAnno)                         
plotAnnoPie(peakAnno)  

peakAnnoList <- list()
peak0hr <- readPeakFile("Chip-Data/GSM862178_t0_export.txt_unique.EE_W100_P0.957.bed")
peak1hr <- readPeakFile("Chip-Data/GSM862179_t1_export.txt_unique.EE_W100_P0.951.bed")
peak4hr <- readPeakFile("Chip-Data/GSM862180_t4_export.txt_unique.EE_W100_P0.96.bed")
peak24hr <- readPeakFile("Chip-Data/GSM862181_t24_export.txt_unique.EE_W100_P0.951.bed")




peakAnnoList[[1]] <- annotatePeak(peak0hr, tssRegion=c(-3000, 3000),level = "transcript", TxDb = txdb) 
peakAnnoList[[2]] <- annotatePeak(peak1hr, tssRegion=c(-3000, 3000),level = "transcript", TxDb = txdb) 
peakAnnoList[[3]] <- annotatePeak(peak4hr, tssRegion=c(-3000, 3000),level = "transcript", TxDb = txdb) 
peakAnnoList[[4]] <- annotatePeak(peak24hr, tssRegion=c(-3000, 3000),level = "transcript", TxDb = txdb) 


#peakfiles <- c("Chip-Data/GSM862178_t0_export.txt_unique.EE_W100_P0.957.bed", "Chip-Data/GSM862179_t1_export.txt_unique.EE_W100_P0.951.bed","Chip-Data/GSM862180_t4_export.txt_unique.EE_W100_P0.96.bed","Chip-Data/GSM862181_t24_export.txt_unique.EE_W100_P0.951.bed")
#peakAnnoList <- lapply(peakfiles, annotatePeak)


genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
AnnotationChipseekerAllData <- lapply(peakAnnoList, function(i) as.data.frame(i))
par(mfrow = c(2,2), mar = c(3,2,2,2))
for (i in 1:4){
  hist(AnnotationChipseekerAllData[[i]]$distanceToTSS,main = paste(i,"hour"),breaks = 200,xlim= c(-100000,100000))
}
#get the list of the genes that have a peak in their 20kb distance (upstream or downstream)
genes20kb <- list()
for (i in 1:4){
  genes20kb[[i]] <- genes[[i]][which(abs(AnnotationChipseekerAllData[[i]]$distanceToTSS) <20000)]
}
Uniquegenes20kb= lapply(genes20kb, function(i) unique(i))
unionChipgenes20k <- character()
for (i in 1:4){
  unionChipgenes20k <- union(unionChipgenes20k, Uniquegenes20kb[[i]])
}
##########get the union of all genes that are assigned to a chip peak (no distance threshold)
Uniquegenes= lapply(genes, function(i) unique(i))

####################################################################################
####################################################################################
####################################################################################

#get the ensemble id of these genes and compare with the other package results
#should do this sometime, or decide to use one of the packages
#annotatedPeak = annotatePeakInBatch(ERchip1, AnnotationData=TSS.human.NCBI36)
ChIPpeakAnnoList <- list()
ChIPpeakAnnoList[[1]] <- annotatePeakInBatch(ERchip0, AnnotationData=TSS.human.NCBI36)
ChIPpeakAnnoList[[2]] <- annotatePeakInBatch(ERchip1, AnnotationData=TSS.human.NCBI36)
ChIPpeakAnnoList[[3]] <- annotatePeakInBatch(ERchip4, AnnotationData=TSS.human.NCBI36) 
ChIPpeakAnnoList[[4]] <- annotatePeakInBatch(ERchip24, AnnotationData=TSS.human.NCBI36)

genesChIPpeakAnno= lapply(ChIPpeakAnnoList, function(i) as.data.frame(i)$feature)
AllChIPpeakAnno <- lapply(ChIPpeakAnnoList, function(i) as.data.frame(i))

#temp1 = as.data.frame(annotatedPeak1)
#values(annotatedPeak)
#pie(table(temp1$insideFeature))
#table(temp1$insideFeature)

x <- org.Hs.egENSEMBL
xx <- as.list(org.Hs.egENSEMBL2EG)
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
mapped_genes <- mappedkeys(x[as.character(DiffGeneEntrezmod)])
# Convert to a list
xx <- as.list(x[mapped_genes])
entotherapp <- xx[genesChIPpeakAnno[[1]]]
####################################################################################
####################################################################################
####################################################################################
#motif analysis
source("https://bioconductor.org/biocLite.R")  
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
#read the GADEM documentation
#motifs_er=GADEM(RangedData(IRanges(, peak4hr$summit+25), space=space(ets1_pea

#############
#### Get the overlap of genes with ER chip peak nearby them and differentially expressed genes
####################################################
##############################################################################
########################################################################################################
length(DiffGeneEntrez)
sum(is.na(DiffGeneEntrez))
rownames(DiffGeneEntrez)
typeof(DiffGeneEntrez)
#hold the associated illumina id for each gene so the expression can be easily retreived
# the DiffGeneEntrez list has those illumina ids as the names of the list
DiffGeneEntrezmod <- list
#DiffGeneEntrezmodNormal <- list()
#removing the un annotated genes
DiffGeneEntrezmod <- DiffGeneEntrez[-c(which(is.na(DiffGeneEntrez)))]
DiffGeneEntrezmodNormal <- DiffGeneEntrezNormal[-c(which(is.na(DiffGeneEntrezNormal)))]
length(DiffGeneEntrezmod)
#removing the duplicated maps (I should consider looking at the expression profile of these duplicates to if they agree)
which(duplicated.default(DiffGeneEntrezmod))

#names(DiffGeneEntrezmod)
DiffGeneEntrezmod <- DiffGeneEntrezmod[-c(which(duplicated.default(DiffGeneEntrezmod)))]
DiffGeneEntrezmodNormal <- DiffGeneEntrezmodNormal[-c(which(duplicated.default(DiffGeneEntrezmodNormal)))]

length(DiffGeneEntrezmod)
#3276
DiffGeneEntrezmodChar <- as.character(DiffGeneEntrezmod) # <- these are the differentially expressed genes
DiffGeneEntrezmodCharNormal <- as.character(DiffGeneEntrezmodNormal) # <- these are the differentially expressed genes

#getting the union of the genes that were bound by ER at some region

#peakannotation4hr <- as.data.frame(peakAnno)
#peakannotation4hr$geneId
#um(duplicated.default(peakannotation4hr$geneId))
#uniqueBoundgenes <- unique(peakannotation4hr$geneId)

#find the overlap of these ids (uniqueBoundgenes) with differentialy expressed ones


typeof(DiffGeneEntrezmod)
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
Uniquegenes= lapply(genes, function(i) unique(i))
# get the union of bound genes in all chip peaks
unionChipgenes <- character()
for (i in 1:4){
 unionChipgenes <- union(unionChipgenes, Uniquegenes[[i]])
}
#get the intersection of chip associated genes and diff expressed genes
ChipDiffExpressionIntersect <- intersect(unionChipgenes, DiffGeneEntrezmodChar)
ChipDiffExpressionNormIntersect <- intersect(unionChipgenes, DiffGeneEntrezmodCharNormal)

#get the intersection with the genes that have a chip peak in their 20kb
ChipDiffExpressionIntersect20k<- intersect(unionChipgenes20k, DiffGeneEntrezmodChar)
ChipDiffExpressionNormIntersect20k<- intersect(unionChipgenes20k, DiffGeneEntrezmodCharNormal)

length(ChipDiffExpressionIntersect20k) #252
#length(intersect(unionChipgenes, unionChipgenes20k))

length(intersect(unionChipgenes, DiffGeneEntrezmod)) #393
################################################################################################
##################################################
################ Take a look at the expression of the intersection sets
##################################################
################################################################################################

#These are the microarry names of intersection of chip and diff expressed
selDataMatrixMCF7sortedCopyNorm <- scale(selDataMatrixMCF7sortedCopy)
intersectionMicNames <- names(DiffGeneEntrezmod[match(ChipDiffExpressionIntersect,DiffGeneEntrezmod)])
intersectionMicNamesNorm <- names(DiffGeneEntrezmodNormal[match(ChipDiffExpressionNormIntersect,DiffGeneEntrezmodNormal)])
intersectionMicNames20k <- names(DiffGeneEntrezmod[match(ChipDiffExpressionIntersect20k,DiffGeneEntrezmod)])
intersectionMicNames20kNorm <- names(DiffGeneEntrezmodNormal[match(ChipDiffExpressionNormIntersect20k,DiffGeneEntrezmodNormal)])

intersectionExp <- selDataMatrixMCF7sortedCopyNorm[intersectionMicNames,]
intersectionExpNorm <- selDataMatrixMCF7sortedCopyNorm[intersectionMicNamesNorm,]
intersectionExp20k <- selDataMatrixMCF7sortedCopyNorm[intersectionMicNames20k,]
intersectionExpNorm20k <- selDataMatrixMCF7sortedCopyNorm[intersectionMicNames20kNorm,]
par(mfrow = (c(10,10)), mar = c(0.1,0.1,0.1,0.1))
for (i in 1:length(intersectionMicNames)){
  plot(selDataMatrixMCF7sortedCopyNorm[intersectionMicNames[i],], ylim = c(-4,4), type = "l",xlab = "",ylab = "",xaxt = "n",yaxt = "n")
  abline (h = 0, col = 2)
  }
for (i in 1:length(intersectionMicNamesNorm)){
  plot(selDataMatrixMCF7sortedCopyNorm[intersectionMicNamesNorm[i],], ylim = c(-4,4), type = "l",xlab = "",ylab = "",xaxt = "n",yaxt = "n")
  abline (h = 0, col = 2)
}  
for (i in 1:length(intersectionMicNames20kNorm)){
  plot(selDataMatrixMCF7sortedCopyNorm[intersectionMicNames20kNorm[i],], ylim = c(-4,4), type = "l",xlab = "",ylab = "",xaxt = "n",yaxt = "n")
  abline (h = 0, col = 2)
}


ExpCluster$size
nfolds <-15
ExpClusterNorm<- kmeans(intersectionExpNorm ,nfolds , iter.max = 400)
#plot the genes in each group separately
expclusListNorm <- list()

for (eaf in 1:nfolds){
  expclusListNorm[[eaf]]<- which(ExpClusterNorm$cluster %in% eaf)
}
par (mfrow = c(4,4))
for(eaf in 1:nfolds){
  #range(intersectionExp[expclusList[[eaf]],])
  plot(intersectionExpNorm[expclusListNorm[[eaf]][1],], ylim =c(-1,5) , type = "l",xlab = "",ylab = "",xaxt = "n",yaxt = "n")
  abline(h=0, lwd = 2,col=2, lty =2)
  for (j in 2:length(expclusListNorm[[eaf]]))
    
    lines(intersectionExpNorm[expclusListNorm[[eaf]][j],],col = col_vector[j])
}


##################GRO-seq data
####################################################################
######################################################################################################
biocLite("edgeR")
library(edgeR)
biocLite("preprocessCore")
library(preprocessCore)
GROreg <- read.table("Expression data/GRO-seq/GSE27463_RefSeq.reg.tsv",header = T)
par(mfrow =c(1,1), mar = c(3,3,3,3))
pie(table(GROreg$Cluster))
nrow(GROreg) #4691
boxplot.matrix(as.matrix(GROreg[,7:10]), ylim = c(0,400))
GROexp <- as.matrix(GROreg[,7:10])
colnames(GROexp) <- c("0min","10min","40min","160min")
rownames(GROexp) <- GROreg$RefSeqID
heatmap(GROexp,Colv = F)

GROexpQuaNor <- normalize.quantiles(GROexp)
GROexpQuaNorScale <- t(scale(t(GROexpQuaNor)))
GROexplibsizeNor$samples
colMeans(GROexp)
apply(GROexp, 2 ,sd )
GROexpScale <- scale(GROexp)


GROexpRel <- matrix(nrow = nrow(GROexp), ncol = ncol(GROexp))
for(i in 1:ncol(GROexp)){
  GROexpRel[,i] <- GROexp[,i]/GROexp[,1]
}
heatmap(GROexpRel)
#plot expressions in each group

boxplot.matrix(GROexpQuaNorScale[which(GROreg$Cluster == 1),], ylim = c(-2,4))
boxplot.matrix(GROexpQuaNorScale[which(GROreg$Cluster == 2),], ylim = c(-2,4))
boxplot.matrix(GROexpQuaNorScale[which(GROreg$Cluster == 3),], ylim = c(-2,4))
boxplot.matrix(GROexpQuaNorScale[which(GROreg$Cluster == 4),], ylim = c(-2,4))
#plot the quantile normalized data
par(mfrow = c(2,2))
for (i in 1:4){
  plot(GROexpQuaNorScale[which(GROreg$Cluster == i)[1],],ylim = c(-1,2),type = "l")
  mtext(text = paste(length(which(GROreg$Cluster == i)), "transcripts"),side = 3)
  for(j in 2:length(which(GROreg$Cluster == i))){
    lines(GROexpQuaNorScale[which(GROreg$Cluster == i)[j],],col= col_vector[j%%length(col_vector)])
  }
  lines(colMeans(GROexpQuaNorScale[which(GROreg$Cluster == i),],na.rm = T),lwd = 4, col = 1)
}
####look at duplicates:
#make a list of transcripts that map to each gene.
nonDupList <- list()
names(nonDupList)<- unique(as.character(GROreg$GeneSymbol))
# for ( i in 1:length(unique(as.character(GROreg$GeneSymbol)))){
#   nonDupList[[i]] <- GROreg$RefSeqID[as.character(GROreg$GeneSymbol) == unique(as.character(GROreg$GeneSymbol))[i]]
# }

#get the number of duplicates for each gene
sumDup <-numeric()
#nonDupList <- list()
for(i in 1:length(unique(as.character(GROreg$GeneSymbol)))){
  sumDup[i] <- sum(!is.na(match(as.character(GROreg$GeneSymbol),unique(as.character(GROreg$GeneSymbol))[i])))
  nonDupList[[i]] <- which(!is.na(match(as.character(GROreg$GeneSymbol),unique(as.character(GROreg$GeneSymbol))[i])))
}

par(mfrow = c(10,10), mar = c(1,1,1,1))
for(i in 1:length(nonDupList)){
  if (sumDup[i] != 1){
    plot(GROexpQuaNorScale[nonDupList[[i]][1],],type = ("l"),ylim = c(-1,2),main = paste(i),xaxt = "n",yaxt = "n")
    for(j in 2:length(nonDupList[[i]])){
      lines(GROexpQuaNorScale[nonDupList[[i]][j],],col = col_vector[j])
    }
  }
}

#index of the genes with multiple transcripts where their transcripts don't show the same transcriptional profile
#index in the nonDupList
DupNotSim <- c(35,759,983,1333,1475,1505,1589,1742,1751,1820,1853,1973,2002,2119,2363,2413,2648,2693,2697,2798,2857,3014 ,3055)

#get one transcription profile per gene.
GROexpQuaNorScalePerGene <- matrix(nrow = length(nonDupList), ncol = 4)
rownames(GROexpQuaNorScalePerGene) <- names(nonDupList)
colnames(GROexpQuaNorScalePerGene) <- c("0min","10min","40min","160min")
for(i in 1:length(nonDupList)){
  if (sumDup[i] != 1){
    GROexpQuaNorScalePerGene[i,] <- colMeans(GROexpQuaNorScale[nonDupList[[i]],],na.rm = T)
  }else{
    GROexpQuaNorScalePerGene[i,] <- GROexpQuaNorScale[nonDupList[[i]],]
  }
}
sum(is.na(GROexpQuaNorScalePerGene)) #8 

EntrezForDifGRO <- select(x = org.Hs.eg.db, keys = rownames(GROexpQuaNorScalePerGene), columns = c("ENTREZID","GENENAME"), keytype = "SYMBOL")
#be carefull with the genes listed above as "DupNotSim" because those don't have same profiles for 
#their different transcripts.
#intersection of ChipData and GROseq data:
# biocLite("org.Hs.egALIAS2EG")
# library("org.Hs.egALIAS2EG")
# x <- org.Hs.egALIAS2EG
# xvx <- as.list(org.Hs.egALIAS2EG)
# xvx <- xvx[!is.na(xvx)]
# as.list(x[rownames(GROexpQuaNorScalePerGene)])
# typeof(mappedsym)
# mappedsym <- mappedkeys(org.Hs.egALIAS2EG)
# as.list(org.Hs.egALIAS2EG[rownames(GROexpQuaNorScalePerGene)])
# biocLite("org.Hs.egREFSEQ")
# refseqEntrez <- as.list(org.Hs.egREFSEQ)
# mapped_seqs <- mappedkeys(org.Hs.egREFSEQ)
# mapped_seqs[which(mapped_seqs == unionChipgenes20k)]
# unionChipgenes20k
# refseqEntrez[unionChipgenes20k]

unionChipgenes20kSymbol <- getSYMBOL(unionChipgenes20k, data='org.Hs.eg')
length(unique(unionChipgenes20kSymbol))

length(intersect(unique(unionChipgenes20kSymbol),unique(as.character(GROreg$GeneSymbol[which(GROreg$Cluster == 2)])) ))
unique(as.character(GROreg$GeneSymbol[which(GROreg$Cluster == 2)]))

#read a new chip seq dataset
ChipWelboren <- read.table("Chip-Data/GSM365926.peaks.txt",header = T)
rownames(ChipWelboren) <- ChipWelboren$id
ChipWelboren <- ChipWelboren[,2:6]
Chromos <- character()
ChromosIn <- as.character(ChipWelboren$chrom)
for(i in 1:nrow(ChipWelboren)){
  Chromos[i] <- paste("chr",ChromosIn[i],sep = "" )
}
ChipWelboren2 <- ChipWelboren[,2:5]
ChipWelboren2 <- cbind(Chromos, ChipWelboren2)
ChipWelboren2 <- cbind(ChipWelboren2, rep(".",nrow(ChipWelboren2)))
ChipWelboren2$`rep(".", nrow(ChipWelboren2))`[1:10]
ChipWelboren2$Chromos[1:10]
typeof(ChipWelboren2$Chromos)
#after changing the coordinates to hg19
ChipWelboren2hg19gr <- readPeakFile("Chip-Data/ChipWelboren2hg19.bed",as = "GRanges")

write.table(ChipWelboren2,quote = F,file = "ChipWelboren2.bed",row.names = F,col.names = F)


#apply(ChipWelboren$chrom,FUN = paste("Chr",x,sep = ""), ChipWelboren$chrom[x])


ERchipNewSet <- read.bed("Chip-Data/ChipWelboren2.bed")

ChIPpeakAnnoNewDataset<- annotatePeakInBatch(ERchipNewSet, AnnotationData=TSS.human.NCBI36)
ChIPpeakAnnoNewDatasetdf <- as.data.frame(ChIPpeakAnnoNewDataset)
ChIPpeakAnnoNewDatasetdf$feature[1:10]
peakWelboren <- readPeakFile("Chip-Data/ChipWelboren2.bed")




peakAnnoList[[1]] <- annotatePeak(peak0hr, tssRegion=c(-3000, 3000),level = "transcript", TxDb = txdb) 
peakAnnoList[[2]] <- annotatePeak(peak1hr, tssRegion=c(-3000, 3000),level = "transcript", TxDb = txdb) 
peakAnnoList[[3]] <- annotatePeak(peak4hr, tssRegion=c(-3000, 3000),level = "transcript", TxDb = txdb) 
peakAnnoList[[4]] <- annotatePeak(peak24hr, tssRegion=c(-3000, 3000),level = "transcript", TxDb = txdb) 


#peakfiles <- c("Chip-Data/GSM862178_t0_export.txt_unique.EE_W100_P0.957.bed", "Chip-Data/GSM862179_t1_export.txt_unique.EE_W100_P0.951.bed","Chip-Data/GSM862180_t4_export.txt_unique.EE_W100_P0.96.bed","Chip-Data/GSM862181_t24_export.txt_unique.EE_W100_P0.951.bed")
#peakAnnoList <- lapply(peakfiles, annotatePeak)


genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
AnnotationChipseekerAllData <- lapply(peakAnnoList, function(i) as.data.frame(i))
par(mfrow = c(2,2), mar = c(3,2,2,2))
for (i in 1:4){
  hist(AnnotationChipseekerAllData[[i]]$distanceToTSS,main = paste(i,"hour"),breaks = 200,xlim= c(-100000,100000))
}
#get the list of the genes that have a peak in their 20kb distance (upstream or downstream)
genes20kb <- list()
for (i in 1:4){
  genes20kb[[i]] <- genes[[i]][which(abs(AnnotationChipseekerAllData[[i]]$distanceToTSS) <20000)]
}
Uniquegenes20kb= lapply(genes20kb, function(i) unique(i))
unionChipgenes20k <- character()
for (i in 1:4){
  unionChipgenes20k <- union(unionChipgenes20k, Uniquegenes20kb[[i]])
}
##########get the union of all genes that are assigned to a chip peak (no distance threshold)
Uniquegenes= lapply(genes, function(i) unique(i))


#now I want to also add the genes that have an ER chip peak interacting with them at any time point
##################
####################################################################
######################################################################################################
#Hi-C data, Chiapet data

