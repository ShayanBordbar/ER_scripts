#Associating genes with ER chip peak.
##Genes with ER chip peaks in their +/- 10kb or +/- 20kb
##Genes interacting with an ER chip peak in their +/- 10kb or +/- 20kb
#Read the bed files for each ER chip peak from cistrome and the time point
#for all genes check if its ER+ or not
#genes that are ER positive at least in two sets will be considred for further studies
ERpickCistrome <- list()
ERCistromefiles <- list.files(path="Chip-Data/Cistrome/ESR1", pattern="*.txt", full.names=T, recursive=FALSE)
for (i in 1:length(ERCistromefiles)){
  ERpickCistrome[[i]] <- readPeakFile(peakfile =ERCistromefiles[i],as = "GRanges" )
}
names(ERpickCistrome) <- ERCistromefiles
#also add the time point chip data
peak0hr <- readPeakFile("Chip-Data/GSM862178_t0_export.txt_unique.EE_W100_P0.957.bed",as = "GRanges")
peak1hr <- readPeakFile("Chip-Data/GSM862179_t1_export.txt_unique.EE_W100_P0.951.bed",as = "GRanges" )
peak4hr <- readPeakFile("Chip-Data/GSM862180_t4_export.txt_unique.EE_W100_P0.96.bed",as = "GRanges" )
peak24hr <- readPeakFile("Chip-Data/GSM862181_t24_export.txt_unique.EE_W100_P0.951.bed",as = "GRanges" )

ERpickCistromeTimePoint <- c(ERpickCistrome,peak0hr,peak1hr,peak4hr,peak24hr)
names(ERpickCistromeTimePoint)[15:18] <- c("ER0h","ER1h","ER4h","ER24h")
#ERpickCistromeTimePoint is a list containing 18 chip seq data for ER
ERpickCistromeTimePointdf <- sapply(ERpickCistromeTimePoint, as.data.frame)
#Get the overlap of all genes with each of these chip sets
#make the GRange object for tss +/- 10kb
TSSALL20kb <- tssdf
TSSALL20kb$start <- TSSALL20kb$start - 10000
TSSALL20kb$start[which(TSSALL20kb$start < 0 )] <- 1
TSSALL20kb$end <- TSSALL20kb$end + 10000
TSSALL20kb$width <- TSSALL20kb$end - TSSALL20kb$start
TSSALL20kbGR <- makeGRangesFromDataFrame(TSSALL20kb)
#make the GRange object for tss +/- 20kb
TSSALL40kb <- tssdf
TSSALL40kb$start <- TSSALL40kb$start - 20000
TSSALL20kb$start[which(TSSALL20kb$start < 0 )] <- 1
TSSALL40kb$end <- TSSALL40kb$end + 20000
TSSALL40kb$width <- TSSALL40kb$end - TSSALL40kb$start
TSSALL40kbGR <- makeGRangesFromDataFrame(TSSALL40kb)
#getting the overlap of the promoters and each chip peak
Gene20kERpickCistromeTimePointOverlapList <- list()
Gene40kERpickCistromeTimePointOverlapList <- list()

for(i in 1:length(ERpickCistromeTimePoint)){
  Gene20kERpickCistromeTimePointOverlapList[[i]] <- findOverlaps(query = TSSALL20kbGR,subject =ERpickCistromeTimePoint[[i]])
  Gene40kERpickCistromeTimePointOverlapList[[i]] <- findOverlaps(query = TSSALL40kbGR,subject =ERpickCistromeTimePoint[[i]])
}
#making a matrix binary if the chip set has a peak in the region defined for this transcript
TranscriptErChipMat20Kb <- matrix(0L,nrow = nrow(TSSALL20kb),ncol = 18)
TranscriptErChipMat40Kb <- matrix(0L,nrow = nrow(TSSALL40kb),ncol = 18)

for(i in 1:length(ERpickCistromeTimePoint)){
  TranscriptErChipMat20Kb[Gene20kERpickCistromeTimePointOverlapList[[i]]@queryHits,i] <- 1
  TranscriptErChipMat40Kb[Gene40kERpickCistromeTimePointOverlapList[[i]]@queryHits,i] <- 1
  
}
sum(rowSums(TranscriptErChipMat20Kb) > 10 )
colnames(TranscriptErChipMat20Kb)
colnames(TranscriptErChipMat20Kb) <- names(ERpickCistromeTimePoint)
colnames(TranscriptErChipMat40Kb) <- names(ERpickCistromeTimePoint)
###########
#Comparing the cistrome file with the original Welboren file
#make the GRanges object for the welboren file

#welborenGR <- makeGRangesFromDataFrame(ChipWelboren2hg19)
OriginalCistromeOverlap <- findOverlaps(query = ChipWelboren2hg19gr,subject =ERpickCistromeTimePoint[[1]])
length(OriginalCistromeOverlap@queryHits)
#write the chip files to be able to  onvert them to hg19
write.table(x = ChipWelboren2,file ="ChipWelboren2",quote = F,sep = "\t",row.names = F,col.names = F )
for(i in 15:18){
  write.table(x = ERpickCistromeTimePointdf[[i]][,c(1,2,3,6)],file =paste("TimePoint",as.character(i),".bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F )
}
#read the hg19 of the timepoint data
peak0hr <- readPeakFile("Chip-Data/Time0hg19.bed",as = "GRanges")
peak1hr <- readPeakFile("Chip-Data/Time1hg19.bed",as = "GRanges" )
peak4hr <- readPeakFile("Chip-Data/Time4hg19.bed",as = "GRanges" )
peak24hr <- readPeakFile("Chip-Data/Time24hg19.bed",as = "GRanges" )

ERpickCistromeTimePoint <- c(ERpickCistrome,peak0hr,peak1hr,peak4hr,peak24hr)
names(ERpickCistromeTimePoint)[15:18] <- c("ER0h","ER1h","ER4h","ER24h")
ERpickCistromeTimePointdf <- sapply(ERpickCistromeTimePoint, as.data.frame)
#################
#how to choose transcripts?
#1 is 1hr E2 treatment, 2:9 are 45mins E2 treatment, 10,11 are 6 hr E2, 12,13,14,15 control, 16,17,18 are timepoint

sum(rowSums(TranscriptErChipMat20Kb[,1:9]) > 4 )
sum(rowSums(TranscriptErChipMat20Kb[,c(10,11)]) > 1 )
sum(rowSums(TranscriptErChipMat20Kb[,12:15]) > 2 )
sum(rowSums(TranscriptErChipMat20Kb[,16:18]) > 1 )


Un1 <- union(which(rowSums(TranscriptErChipMat20Kb[,1:9]) > 4 ),which(rowSums(TranscriptErChipMat20Kb[,c(10,11)]) > 1 ))
Un2 <- union(which(rowSums(TranscriptErChipMat20Kb[,12:15]) > 2 ),which(rowSums(TranscriptErChipMat20Kb[,16:18]) > 1))
Un3 <-union(Un1,Un2)
length(Un2)
#########
#READ THE XLS FILES OF THE CISTROME DATABASE
Cistrome1 <- read.table("Chip-Data/Human_MCF-7_ESR1_17b-E2-1hr_Stunnenberg.xls",header = T,stringsAsFactors = F)
summary(Cistrome1$X.10.log10.pvalue)
sum(Cistrome1$FDR... < 5 )
################
#Find the regions interacting with each 20kb region around TSS and see if they contain ER chip peaks
#add 2 to the chip matrix of a region for that chip profile if it has an interacting peak


#Get the hits for each 20kb TSS region in the first side of the interactions
OverlapTSS20kbInt1 <- findOverlaps(query = TSSALL20kbGR ,subject =InteractingOneGR)
#Get the hits for each 20kb TSS region in the second side of the interactions
OverlapTSS20kbInt2 <- findOverlaps(query = TSSALL20kbGR ,subject =InteractingTwoGR)
#Get the hits for each chip peak in the first side of the interactions
OverlapChipInt1 <- list()
for (i in 1:length(ERpickCistromeTimePoint)){
  OverlapChipInt1[[i]] <- findOverlaps(query = ERpickCistromeTimePoint[[i]] ,subject =InteractingOneGR)
}
#Get the hits for each chip peak in the second side of the interactions
OverlapChipInt2 <- list()
for (i in 1:length(ERpickCistromeTimePoint)){
  OverlapChipInt2[[i]] <- findOverlaps(query = ERpickCistromeTimePoint[[i]] ,subject =InteractingTwoGR)
}
#Create the matrix each row on 20kb TSS region, each column a Chip condition, 
#Each entry shows the number of interacting chip peaks with the TSS region in this chip condition
OverLap20kbChip <- matrix(0L,nrow = nrow(TSSALL20kb),ncol = length(ERpickCistromeTimePoint))
#Each entry shows the number of interacting chip peaks with the TSS region in this chip condition from intA to intB
OverLap20kbChip1 <- matrix(0L,nrow = nrow(TSSALL20kb),ncol = length(ERpickCistromeTimePoint))
#Each entry shows the number of interacting chip peaks with the TSS region in this chip condition from intB to intA
OverLap20kbChip2 <- matrix(0L,nrow = nrow(TSSALL20kb),ncol = length(ERpickCistromeTimePoint))


for (i in 1:nrow(TSSALL20kb)){
  currInt1 <- OverlapTSS20kbInt1@subjectHits[OverlapTSS20kbInt1@queryHits == i]
  currInt2 <- OverlapTSS20kbInt2@subjectHits[OverlapTSS20kbInt2@queryHits == i]
  #if there is any hit for this 20kb tss region in the first side of the interactions
  #find the chip peaks that are involved in those interactions from the other side
  print(paste("TSS no#",as.character(i),"number of hits in intA:", as.character(length(currInt1))))
  if(length(currInt1) > 0){
    currNumIntChp1 <- matrix(0L,nrow = length(currInt1),ncol = length(OverlapChipInt1))
    for (k in 1:length(currInt1)){
      for(j in 1: length(OverlapChipInt1)){
        currNumIntChp1[k,j] <- length(unique(OverlapChipInt1[[j]]@queryHits[OverlapChipInt1[[j]]@subjectHits == currInt1[k]]))
      }
    }
    OverLap20kbChip1[i,] <- colSums(currNumIntChp1)
  }
  #if there is any hit for this 20kb tss region in the second side of the interactions
  #find the number of chip peaks that are involved in those interactions from the other side
  print(paste("TSS no#",as.character(i),"number of hits in intB:", as.character(length(currInt2))))
  
  if(length(currInt2) > 0){
    currNumIntChp2 <- matrix(0L,nrow = length(currInt2),ncol = length(OverlapChipInt2))
    
    for (k in 1:length(currInt2)){
      for(j in 1: length(OverlapChipInt2)){
        currNumIntChp2[k,j] <- length(unique(OverlapChipInt2[[j]]@queryHits[OverlapChipInt2[[j]]@subjectHits == currInt2[k]]))
      }
    }
    OverLap20kbChip2[i,] <- colSums(currNumIntChp2)
  }

  }

OverLap20kbChip <- OverLap20kbChip1 + OverLap20kbChip2
colnames(OverLap20kbChip) <- names(ERpickCistromeTimePointdf)
#The binary version
OverLap20kbChipBinary <- matrix(0L,nrow = nrow(TSSALL20kb),ncol = length(ERpickCistromeTimePoint))
OverLap20kbChipBinary[OverLap20kbChip > 0] = 1
summary(rowSums(OverLap20kbChipBinary))

which(rowSums(OverLap20kbChipBinary) == 18)
sum(rowSums(OverLap20kbChipBinary[,1:9]) > 4 )
sum(rowSums(OverLap20kbChipBinary[,c(10,11)]) > 1 )
sum(rowSums(OverLap20kbChipBinary[,12:15]) > 2 )
sum(rowSums(OverLap20kbChipBinary[,16:18]) > 1 )


UnInt1 <- union(which(rowSums(OverLap20kbChipBinary[,1:9]) > 4 ),which(rowSums(OverLap20kbChipBinary[,c(10,11)]) > 1 ))
UnInt2 <- union(which(rowSums(OverLap20kbChipBinary[,12:15]) > 2 ),which(rowSums(OverLap20kbChipBinary[,16:18]) > 1))
UnInt3 <-union(UnInt1,UnInt2)
length(UnInt3)

length(union(UnInt3, Un3))
ERchipAssociatedTranscripts <- sort(union(UnInt3, Un3),decreasing = F )#This is the index of all ER associated transcripts

#convert transcript names to gene symbols
mapTranscriptsToGenes <- select(txdb,keys = as.character(TSSALL20kb$tx_id) ,columns = c("TXNAME","GENEID"),keytype = "TXID")
ERchipAssociatedGENEIDs <- mapTranscriptsToGenes$GENEID[ERchipAssociatedTranscripts]
names(ERchipAssociatedGENEIDs) <- as.character(ERchipAssociatedTranscripts)
#for each gene get all of the information about its transcripts and sum them for the gene
UniqueGeneEntrezID <- unique(mapTranscriptsToGenes$GENEID)
UniqueGeneEntrezID <- UniqueGeneEntrezID[-c(3)]
#matrix for gene chip interaction
OverLap20kbChipBinaryGene <- matrix(0L,nrow = length(UniqueGeneEntrezID),ncol = ncol(OverLap20kbChipBinary))
#matrix for gene chip vicinity
TranscriptErChipMat20KbGene <- matrix(0L,nrow = length(UniqueGeneEntrezID),ncol = ncol(OverLap20kbChipBinary))
for(i in 1:length(UniqueGeneEntrezID)){
  CurrTrans <- mapTranscriptsToGenes$TXID[which(mapTranscriptsToGenes$GENEID == UniqueGeneEntrezID[i])]
  if(length(CurrTrans) > 1){
    OverLap20kbChipBinaryGene[i,] <- colSums(OverLap20kbChipBinary[CurrTrans,])
    TranscriptErChipMat20KbGene[i,] <- colSums(TranscriptErChipMat20Kb[CurrTrans,])
  }
  else if(length(CurrTrans) == 1){
    OverLap20kbChipBinaryGene[i,] <- OverLap20kbChipBinary[CurrTrans,]
    TranscriptErChipMat20KbGene[i,] <- TranscriptErChipMat20Kb[CurrTrans,]
  }
 
}
#the Binary versions
OverLap20kbChipBinaryGeneBinary <-  matrix(0L,nrow = length(UniqueGeneEntrezID),ncol = ncol(OverLap20kbChipBinary))
TranscriptErChipMat20KbGeneBinary <- matrix(0L,nrow = length(UniqueGeneEntrezID),ncol = ncol(OverLap20kbChipBinary))
OverLap20kbChipBinaryGeneBinary[OverLap20kbChipBinaryGene > 0] = 1
TranscriptErChipMat20KbGeneBinary[TranscriptErChipMat20KbGene > 0 ] = 1

sum(rowSums(OverLap20kbChipBinaryGeneBinary) > 0)
summary(rowSums(TranscriptErChipMat20KbGeneBinary))
#make the matrix for gene ERchip vicinity and interaction
GeneERchipVicinityPlusInteraction <- OverLap20kbChipBinaryGeneBinary + TranscriptErChipMat20KbGeneBinary
summary(rowSums(GeneERchipVicinityPlusInteraction))
GeneERchipVicinityPlusInteractionBinary <- matrix(0L,nrow = length(UniqueGeneEntrezID),ncol = ncol(OverLap20kbChipBinary))
GeneERchipVicinityPlusInteractionBinary[GeneERchipVicinityPlusInteraction > 0] = 1
#Get the genes with peaks in multple chip datasets
sum(rowSums(GeneERchipVicinityPlusInteractionBinary[,1:9]) > 4 )
sum(rowSums(GeneERchipVicinityPlusInteractionBinary[,c(10,11)]) > 1 )
sum(rowSums(GeneERchipVicinityPlusInteractionBinary[,12:15]) > 2 )
sum(rowSums(GeneERchipVicinityPlusInteractionBinary[,16:18]) > 1 )


GenChipInt1 <- union(which(rowSums(GeneERchipVicinityPlusInteractionBinary[,1:9]) > 4 ),which(rowSums(GeneERchipVicinityPlusInteractionBinary[,c(10,11)]) > 1 ))
GenChipInt2 <- union(which(rowSums(GeneERchipVicinityPlusInteractionBinary[,12:15]) > 2 ),which(rowSums(GeneERchipVicinityPlusInteractionBinary[,16:18]) > 1))
GenChipInt3 <-union(GenChipInt1,GenChipInt2)
length(GenChipInt3)#4039
ERchipAssociatedGENEIDsFinal <- UniqueGeneEntrezID[sort(GenChipInt3,decreasing = F)]


# length(unique(ERchipAssociatedGENEIDs))
# sum(is.na(ERchipAssociatedGENEIDs))
# sum(duplicated(ERchipAssociatedGENEIDs))
# mapTranscriptsToGenesENT <- select(txdb,keys = as.character(TSSALL20kb$tx_id) ,columns = c("TXNAME","GENEID","SYMBOL","ENTREZID"),keytype = "TXID")
# 
# mapTranscriptsToGenes2 <- select(x = org.Hs.eg.db, keys = as.character(TSSALL20kb$tx_id), columns = c("TXNAME","GENEID","SYMBOL","ENTREZID"), keytype = "TXID")



########RNA-seq knock down Datasets
########################################################
################################################################################################
########################################################################################
#THis is microarray not RNA seq : PBXkd
PBXkdRNAseq <- read.delim("Expression data/RNAseq/GSE39417_non-normalized.txt",header = T,sep = "\t")
nrow(PBXkdRNAseq)
PBXkdRNAseqMat <- as.matrix(PBXkdRNAseq[,4:19])
#remove the p-value columns
PBXkdRNAseqMat <- PBXkdRNAseqMat[,seq(1,15,by=2)]
#z score normalize each column 
PBXkdRNAseqMatNorm <- scale(PBXkdRNAseqMat) 
sum(PBXkdRNAseqMatNorm < -1)
boxplot.matrix(PBXkdRNAseqMatNorm,las = 2,ylim = c(-1,2))
PBXkdRNAseq$PROBE_ID[1:10]

########CTCF knock down GSM2257523
GSE53532RNAseq00 <- read.delim(file = "Expression data/RNAseq/GSM2257523_Genes_fpkm_siNT_E0h_cufflinks.txt",header = F,sep = "\t",stringsAsFactors = F)
GSE53532RNAseq01 <- read.delim(file = "Expression data/RNAseq/GSM2257525_Genes_fpkm_siNT_E3h_cufflinks.txt",header = F,sep = "\t",stringsAsFactors = F)
GSE53532RNAseq10 <- read.delim(file = "Expression data/RNAseq/GSM2257524_Genes_fpkm_siCTCF_E0h_cufflinks.txt",header = F,sep = "\t",stringsAsFactors = F)
GSE53532RNAseq11 <- read.delim(file = "Expression data/RNAseq/GSM2257526_Genes_fpkm_siCTCF_E3h_cufflinks.txt",header = F,sep = "\t",stringsAsFactors = F)

sum(GSE53532RNAseq00$V1 == GSE53532RNAseq01$V1)
length(unique(GSE53532RNAseq00$V1))
#make the expression matrix
GSE53532RNAseq <- as.matrix(cbind(as.numeric(GSE53532RNAseq00$V2),as.numeric(GSE53532RNAseq01$V2),as.numeric(GSE53532RNAseq10$V2),as.numeric(GSE53532RNAseq11$V2)))
#remove the line which is not read
GSE53532RNAseq <- GSE53532RNAseq[-c(24667),]
GSE53532RNAseq00 <- GSE53532RNAseq00[-c(24667),]
GSE53532RNAseq01 <- GSE53532RNAseq01[-c(24667),]
GSE53532RNAseq10 <- GSE53532RNAseq10[-c(24667),]
GSE53532RNAseq11 <- GSE53532RNAseq11[-c(24667),]
#Zscore normalize the matrix
GSE53532RNAseqNorm <- scale(GSE53532RNAseq)
rownames(GSE53532RNAseqNorm) <- GSE53532RNAseq01$V1
rownames(GSE53532RNAseq) <- GSE53532RNAseq01$V1

GSE53532RNAseqZNormQuaNorm <- normalize.quantiles(GSE53532RNAseqNorm)

boxplot.matrix(GSE53532RNAseqZNormQuaNorm)
#sum(GSE53532RNAseqNorm[which(GSE53532RNAseqNorm[,1] > 0),2] > 0 )
################
#ERalpha kd :GSE53532
GSE53532_DE_gene <- read.delim("Expression data/RNAseq/GSE53532_DE_gene.tsv",header = T,sep = "\t",stringsAsFactors = F)
GSE53532_DE_gene$EnsemblID[1:10]

################BRD4 kd
GSE55922RNAseq <- read.csv("Expression data/RNAseq/GSE55922_normalized_counts.csv",header = T,stringsAsFactors = F)
GSE55922RNAseqMat <- as.matrix(GSE55922RNAseq[,2:13])
ncol(GSE55922RNAseqMat)
sd(GSE55922RNAseq$cont.e2b)
mean(GSE55922RNAseqMat[1,],na.rm = T)
boxplot.matrix(GSE55922RNAseqMatNorm)
GSE55922RNAseqMatNorm <-scale(GSE55922RNAseqMat)
rownames(GSE55922RNAseqMatNorm) <- GSE55922RNAseq$X
rownames(GSE55922RNAseqMat) <- GSE55922RNAseq$X
GSE55922RNAseq[1,]
#################
#ZNF217 kd : no E2 treatment
GSE58326RNAseq <- read.delim("Expression data/RNAseq/GSE58326_gene_exp.txt",header = T,stringsAsFactors = F)
GSE58326RNAseq$FPKM_siZNF[1:10]

#################
#estradiol- and pro-inflammatory cytokine-treated MCF7 cells.# there is also chip seq and GRO-seq available for this set
GSE67295RNAseq <- read.delim("Expression data/RNAseq/GSE67295GeneExpressionPolyA.txt",header = T,stringsAsFactors = F)
GSE67295RNAseq$Transcript.RepeatID..cmd.analyzeRepeats.pl.rna.hg19..condenseGenes..d..data.hg19.MCF7.RNA.MCF7_RNA_Veh_Rep1_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_Veh_Rep2_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_E2_Rep1_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_E2_Rep2_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_E2_ICI_Rep1_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_E2_ICI_Rep2_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_Il1b_Rep1_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_Il1b_Rep2_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_Il1b_ICI_Rep1_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_Il1b_ICI_Rep2_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_Veh_3h_siCtl_rep1_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_Veh_3h_siCtl_rep2_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_TNFa_3h_siCtl_rep1_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_TNFa_3h_siCtl_rep2_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_TNF_ICI_Rep1_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_TNF_ICI_Rep2_Josh_14_11_03...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_Veh_3h_rep1_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_Veh_3h_rep2_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_E2_3h_rep1_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_E2_3h_rep2_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_E2_TOT_3h_rep1_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_E2_TOT_3h_rep2_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_E2_TOT_Il1b_3h_rep1_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_E2_TOT_Il1b_3h_rep2_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_E2_TOT_TNF_3h_rep1_Josh_14_10_08...data.hg19.MCF7.RNA.MCF7_RNA_RNA_Seq_E2_TOT_TNF_3h_rep2_Josh_14_10_08..[1:10]
colnames(GSE67295RNAseq)[1] <- "TranscriptID"
GSE67295RNAseq$TranscriptID[1:100]
nrow(GSE67295RNAseq)
GSE67295RNAseqMat <- as.matrix(GSE67295RNAseq[,9:34])
rownames(GSE67295RNAseqMat) <- GSE67295RNAseq$TranscriptID
GSE67295RNAseqMatNorm <- scale(GSE67295RNAseqMat)
####################
#read the GSE68358 data #all 24 experiments are E2 or E2 plus something else
GSE68358RNAseq <- list()
GSE68358Files <- list.files(path="Expression data/RNAseq/GSE68358_RAW", pattern=".txt$", full.names=T, recursive=FALSE)

for( i in 1:length(GSE68358Files)){
  GSE68358RNAseq[[i]] <- read.delim(GSE68358Files[[i]],header = F,stringsAsFactors = F)
}
names(GSE68358RNAseq) <-  GSE68358Files

GSE68358RNAseqMat <- matrix(0L,nrow = nrow(GSE68358RNAseq[[1]]), ncol = length(GSE68358Files))
for(i in 1:length(GSE68358RNAseq)){
  GSE68358RNAseqMat[,i] <- GSE68358RNAseq[[i]]$V2
}
colnames(GSE68358RNAseqMat)<- GSE68358Files
rownames(GSE68358RNAseqMat)<- GSE68358RNAseq[[1]]$V1

boxplot.matrix(GSE68358RNAseqMatNorm)
GSE68358RNAseqMatNorm <- scale(GSE68358RNAseqMat)
rownames(GSE68358RNAseqMatNorm)
###############G9 and PHF kd
GSE76507RNAseq <- read.delim("Expression data/RNAseq/GSE76507_shG9a_shPHF20.pool.htseq.txt",header = T,stringsAsFactors = F)
colnames(GSE76507RNAseq) <- c("gene","ensembl","shG9a_NOE2","shG9a_E2","shNT_NOE2","shNT_E2","shNT_NoE2.rep1","shNT_NOE2.rep2","shNT_E2.rep1","shNT_E2.rep2","shPHF20_NOE2.rep1","shPHF20_NOE2.rep2","shPHF20_E2.rep1","shPHF20_E2.rep2")		
GSE76507RNAseqMat <- as.matrix(GSE76507RNAseq[,3:14])
GSE76507RNAseqMatNorm <- scale(GSE76507RNAseqMat)
rownames(GSE76507RNAseqMatNorm) <- GSE76507RNAseq$gene
rownames(GSE76507RNAseqMat) <- GSE76507RNAseq$gene
#################24 hr E2 treatment and before treatment:GSE51403_count.matrix.byGene.round1
GSE51403RNAseq <- read.delim("Expression data/RNAseq/GSE51403_count.matrix.byGene.round1.txt",header = T,stringsAsFactors = F)
#leave this one out for now
#################GSE76453 ZNF143 KO
GSE76453VEH <- read.delim("Expression data/RNAseq/GSE76453_cuffdiff_veh.diff",header = T,stringsAsFactors = F) 
GSE76453siCtl <- read.delim("Expression data/RNAseq/GSE76453_cuffdiff_siCtl.diff",header = T,stringsAsFactors = F) 
GSE76453E2 <- read.delim("Expression data/RNAseq/GSE76453_cuffdiff_E2.diff",header = T,stringsAsFactors = F) 
GSE76453siZNF143 <- read.delim("Expression data/RNAseq/GSE76453_cuffdiff_siZNF143.diff",header = T,stringsAsFactors = F) 
#leave this one out for now
##############GSE78167_RAW.tar timepoint 0, 1, 2, 3, 4, 5, 6, 8, 12 , or 24
GSE78167RNAseq <- list()
GSE78167Files <- list.files(path="Expression data/RNAseq/GSE78169_RAW", pattern=".txt$", full.names=T, recursive=FALSE)

for( i in 1:length(GSE78167Files)){
  GSE78167RNAseq[[i]] <- read.delim(GSE78167Files[[i]],header = T,stringsAsFactors = F)
}
names(GSE78167RNAseq) <-  GSE78167Files

GSE78167RNAseqMat <- matrix(0L,nrow = nrow(GSE78167RNAseq[[1]]), ncol = length(GSE78167Files))
for(i in 1:length(GSE78167RNAseq)){
  GSE78167RNAseqMat[,i] <- GSE78167RNAseq[[i]]$FPKM
}
GSE78167ExpName <- character()
for (i in 1:length(GSE78167Files)){
  GSE78167ExpName[i] <- stri_split_coll(str = colnames(GSE78167RNAseqMat),pattern = "_")[[i]][3]
  
}


colnames(GSE78167RNAseqMat)<- GSE78167ExpName
GSE78167GeneEntrez <- character()
for (i in 1:nrow(GSE78167RNAseq[[1]])){
  GSE78167GeneEntrez[i] <- stri_split_coll(str = GSE78167RNAseq[[1]]$gene_id[i],pattern = "|")[[1]][2]
}
GSE78167GeneSymbol <- character()
for (i in 1:nrow(GSE78167RNAseq[[1]])){
  GSE78167GeneSymbol[i] <- stri_split_coll(str = GSE78167RNAseq[[1]]$gene_id[i],pattern = "|")[[1]][1]
}
  
rownames(GSE78167RNAseqMat)<- GSE78167GeneEntrez
#GSE86316 MEN1 KO
GSE86316RNAseq <- read.delim("Expression data/RNAseq/GSE86316_Expression_transcripts_ShMEN1.txt",header = T,stringsAsFactors = F)
GSE86316RNAseqMat <- as.matrix(GSE86316RNAseq[,3:10])
rownames(GSE86316RNAseqMat) <- GSE86316RNAseq$RefseqID
#GSE80098_MCF7_RNAseq E2, R5020 or both for 12 hours
GSE80098RNAseq <- read.csv("Expression data/RNAseq/GSE80098_MCF7_RNAseq.csv",header = T,stringsAsFactors = F)
GSE80098RNAseqmat <- as.matrix(GSE80098RNAseq[,c(2,3,6,9)])
rownames(GSE80098RNAseqmat) <- GSE80098RNAseq$RefSeq
#GSE73663_inhibitor , GSE73663_timecourse
GSE73663Inh <- read.delim("Expression data/RNAseq/GSE73663_inhibitor.txt",header = T,stringsAsFactors = F)
GSE73663TC <- read.delim("Expression data/RNAseq/GSE73663_timecourse.txt",header = T,stringsAsFactors = F)
GSE73663RNAseqMatInh <- as.matrix(GSE73663Inh[,2:19])
rownames(GSE73663RNAseqMatInh) <- GSE73663Inh$Gene.Symbol
GSE73663RNAseqMatTC <- as.matrix(GSE73663TC[,2:13])
rownames(GSE73663RNAseqMatTC) <- GSE73663TC$Gene.Symbol

#Now all the datasets are imported dataset need to be merged for all genes
#list all the datasets
#GSE53532RNAseqNorm CTCF kd GSE53532RNAseq
#GSE55922RNAseqMatNorm BRD4 kd  GSE55922RNAseqMat
#GSE67295RNAseqMatNorm E2 and other treatments  GSE67295RNAseqMat
#GSE68358RNAseqMatNorm  E2 or E2 plus something else  GSE68358RNAseqMat
#GSE76507RNAseqMatNorm  G9 and PHF kd   GSE76507RNAseqMat
#GSE78167RNAseqMat time points: 0,1,2,3,4,5,6,8,12,24
#GSE86316RNAseqMat MEN1 Kd
#GSE80098RNAseqmat E2, R2050 or both
#GSE73663RNAseqMatInh inhibitors Merdak
#GSE73663RNAseqMatTC time course Merdak

keytypes(org.Hs.eg.db)
ERchipAssociatedGENEIDsFinal #final list of the ER associated ENtrez IDs
mapTranscriptsToGenes <- select(txdb,keys = as.character(TSSALL20kb$tx_id) ,columns = c("TXNAME","GENEID"),keytype = "TXID")
mapTranscriptsToGenes2 <- select(x = org.Hs.eg.db, keys = as.character(TSSALL20kb$tx_id), columns = c("TXNAME","GENEID","SYMBOL","ENTREZID"), keytype = "TXID")

IDrefSample <- select(org.Hs.eg.db,keys = as.character(c(1:1000)) ,columns = c("GENENAME"),keytype = "REFSEQ")

#choose entrezID as the common form and convert everything to that, merge the expressions for the genes associated with ER
GSE53532RNAseqNorm#CTCF kd
mapSymbolToEntrez <- select(x = org.Hs.eg.db, keys = rownames(GSE53532RNAseqNorm), columns = c("ENTREZID"), keytype = "SYMBOL")
GSE53532RNAseqNormCommon<- intersect(mapSymbolToEntrez$ENTREZID,ERchipAssociatedGENEIDsFinal)
GSE53532RNAseqNormMatch <- mapSymbolToEntrez$SYMBOL[match(ERchipAssociatedGENEIDsFinal,mapSymbolToEntrez$ENTREZID )]
names(GSE53532RNAseqNormMatch) <- ERchipAssociatedGENEIDsFinal
GSE53532RNAseqNormMatch <- GSE53532RNAseqNormMatch[-which(is.na(GSE53532RNAseqNormMatch))]
#get the matrix for only the genes common between the chip associated list and dataset
GSE53532commInd <- match(GSE53532RNAseqNormMatch,rownames(GSE53532RNAseqNorm))
GSE53532RNAseqNormComMat <- GSE53532RNAseqNorm[GSE53532commInd,]
GSE53532RNAseqComMat <- GSE53532RNAseq[GSE53532commInd,]

rownames(GSE53532RNAseqNormComMat) <- names(GSE53532RNAseqNormMatch)
rownames(GSE53532RNAseqComMat) <- names(GSE53532RNAseqNormMatch)

####################
GSE55922RNAseqMatNorm #brd4 KD

mapSymbolToEntrezbrd4 <- select(x = org.Hs.eg.db, keys = rownames(GSE55922RNAseqMatNorm), columns = c("ENTREZID"), keytype = "SYMBOL")
GSE55922RNAseqMatNormMatch <- mapSymbolToEntrezbrd4$SYMBOL[match(ERchipAssociatedGENEIDsFinal,mapSymbolToEntrezbrd4$ENTREZID )]
names(GSE55922RNAseqMatNormMatch) <- ERchipAssociatedGENEIDsFinal
GSE55922RNAseqMatNormMatch <- GSE55922RNAseqMatNormMatch[-which(is.na(GSE55922RNAseqMatNormMatch))]
#get the matrix for only the genes common between the chip associated list and dataset
GSE55922commInd <- match(GSE55922RNAseqMatNormMatch,rownames(GSE55922RNAseqMatNorm))
GSE55922RNAseqNormComMat <- GSE55922RNAseqMatNorm[GSE55922commInd,]
GSE55922RNAseqComMat <- GSE55922RNAseqMat[GSE55922commInd,]

rownames(GSE55922RNAseqNormComMat) <- names(GSE55922RNAseqMatNormMatch)
rownames(GSE55922RNAseqComMat) <- names(GSE55922RNAseqMatNormMatch)

##################
GSE67295RNAseqMatNorm #E2 and other treatments

mapNMToEntrez <- select(x = org.Hs.eg.db, keys = rownames(GSE67295RNAseqMatNorm), columns = c("ENTREZID"), keytype = "REFSEQ")
# nrow(GSE67295RNAseqMatNorm)
# sum(is.na(mapNMToEntrez))
GSE67295RNAseqMatNormMatch <- mapNMToEntrez$REFSEQ[match(ERchipAssociatedGENEIDsFinal,mapNMToEntrez$ENTREZID )]
names(GSE67295RNAseqMatNormMatch) <- ERchipAssociatedGENEIDsFinal
GSE67295RNAseqMatNormMatch <- GSE67295RNAseqMatNormMatch[-which(is.na(GSE67295RNAseqMatNormMatch))]
#get the matrix for only the genes common between the chip associated list and dataset
GSE67295commInd <- match(GSE67295RNAseqMatNormMatch,rownames(GSE67295RNAseqMatNorm))
GSE67295RNAseqMatNormComMat <- GSE67295RNAseqMatNorm[GSE67295commInd,]
GSE67295RNAseqMatComMat <- GSE67295RNAseqMat[GSE67295commInd,]

rownames(GSE67295RNAseqMatNormComMat) <- names(GSE67295RNAseqMatNormMatch)
rownames(GSE67295RNAseqMatComMat) <- names(GSE67295RNAseqMatNormMatch)

#rownames(GSE67295RNAseqMatNormComMat)

###################
GSE68358RNAseqMatNorm#  E2 or E2 plus something else
mapEnsembleToEntrez <- select(x = org.Hs.eg.db, keys = rownames(GSE68358RNAseqMatNorm), columns = c("ENTREZID"), keytype = "ENSEMBL")
# nrow(mapEnsembleToEntrez)
# nrow(GSE68358RNAseqMatNorm)
# sum(is.na(mapEnsembleToEntrez))
GSE68358RNAseqMatNormMatch <- mapEnsembleToEntrez$ENSEMBL[match(ERchipAssociatedGENEIDsFinal,mapEnsembleToEntrez$ENTREZID )]
names(GSE68358RNAseqMatNormMatch) <- ERchipAssociatedGENEIDsFinal
GSE68358RNAseqMatNormMatch <- GSE68358RNAseqMatNormMatch[-which(is.na(GSE68358RNAseqMatNormMatch))]
#get the matrix for only the genes common between the chip associated list and dataset
GSE68358commInd <- match(GSE68358RNAseqMatNormMatch,rownames(GSE68358RNAseqMatNorm))
GSE68358RNAseqMatNormComMat <- GSE68358RNAseqMatNorm[GSE68358commInd,]
GSE68358RNAseqMatComMat <- GSE68358RNAseqMat[GSE68358commInd,]

rownames(GSE68358RNAseqMatNormComMat) <- names(GSE68358RNAseqMatNormMatch)
rownames(GSE68358RNAseqMatComMat) <- names(GSE68358RNAseqMatNormMatch)
#nrow(GSE67295RNAseqMatNormComMat)
###################
GSE76507RNAseqMatNorm # G9 and PHF kd

mapG9SymbolToEntrez <- select(x = org.Hs.eg.db, keys = rownames(GSE76507RNAseqMatNorm), columns = c("ENTREZID"), keytype = "SYMBOL")
GSE76507RNAseqMatNormMatch <- mapG9SymbolToEntrez$SYMBOL[match(ERchipAssociatedGENEIDsFinal,mapG9SymbolToEntrez$ENTREZID )]
names(GSE76507RNAseqMatNormMatch) <- ERchipAssociatedGENEIDsFinal
GSE76507RNAseqMatNormMatch <- GSE76507RNAseqMatNormMatch[-which(is.na(GSE76507RNAseqMatNormMatch))]
#get the matrix for only the genes common between the chip associated list and dataset
GSE76507commInd <- match(GSE76507RNAseqMatNormMatch,rownames(GSE76507RNAseqMatNorm))
GSE76507RNAseqMatNormComMat <- GSE76507RNAseqMatNorm[GSE76507commInd,]
GSE76507RNAseqMatComMat <- GSE76507RNAseqMat[GSE76507commInd,]

rownames(GSE76507RNAseqMatNormComMat) <- names(GSE76507RNAseqMatNormMatch)
rownames(GSE76507RNAseqMatComMat) <- names(GSE76507RNAseqMatNormMatch)

##########
#GSE78167RNAseqMat time points: 0,1,2,3,4,5,6,8,12,24

GSE78167CommonEntrez <-  match(ERchipAssociatedGENEIDsFinal,rownames(GSE78167RNAseqMat))
names(GSE78167CommonEntrez) <- ERchipAssociatedGENEIDsFinal
GSE78167CommonEntrez <- GSE78167CommonEntrez[-which(is.na(GSE78167CommonEntrez))]

GSE78167RNAseqMatCommon <- GSE78167RNAseqMat[GSE78167CommonEntrez,]
rownames(GSE78167RNAseqMatCommon) <- names(GSE78167CommonEntrez)


#GSE86316RNAseqMat MEN1 Kd
mapNMToEntrezMEN1Kd <- select(x = org.Hs.eg.db, keys = rownames(GSE86316RNAseqMat), columns = c("ENTREZID"), keytype = "REFSEQ")
# nrow(GSE67295RNAseqMatNorm)
# sum(is.na(mapNMToEntrezMEN1Kd))
GSE86316RNAseqMatMatch <- mapNMToEntrezMEN1Kd$REFSEQ[match(ERchipAssociatedGENEIDsFinal,mapNMToEntrezMEN1Kd$ENTREZID )]
names(GSE86316RNAseqMatMatch) <- ERchipAssociatedGENEIDsFinal
GSE86316RNAseqMatMatch <- GSE86316RNAseqMatMatch[-which(is.na(GSE86316RNAseqMatMatch))]
#get the matrix for only the genes common between the chip associated list and dataset
GSE86316commInd <- match(GSE86316RNAseqMatMatch,rownames(GSE86316RNAseqMat))
GSE86316RNAseqMatComMat <- GSE86316RNAseqMat[GSE86316commInd,]
rownames(GSE86316RNAseqMatComMat) <- names(GSE86316RNAseqMatMatch)
#rownames(GSE67295RNAseqMatNormComMat)

#GSE80098RNAseqmat E2, R2050 or both
mapNMToEntrezE2R2050 <- select(x = org.Hs.eg.db, keys = rownames(GSE80098RNAseqmat), columns = c("ENTREZID"), keytype = "REFSEQ")
GSE80098RNAseqmatMatch <- mapNMToEntrezE2R2050$REFSEQ[match(ERchipAssociatedGENEIDsFinal,mapNMToEntrezE2R2050$ENTREZID )]
names(GSE80098RNAseqmatMatch) <- ERchipAssociatedGENEIDsFinal
GSE80098RNAseqmatMatch <- GSE80098RNAseqmatMatch[-which(is.na(GSE80098RNAseqmatMatch))]
#each rowname has multiple refseq ids : fix this later


#GSE73663RNAseqMatInh inhibitors Merdak
mapSymbolToEntrezMerdakInh <- select(x = org.Hs.eg.db, keys = rownames(GSE73663RNAseqMatInh), columns = c("ENTREZID"), keytype = "SYMBOL")
GSE73663RNAseqMatInhMatch <- mapSymbolToEntrezMerdakInh$SYMBOL[match(ERchipAssociatedGENEIDsFinal,mapSymbolToEntrezMerdakInh$ENTREZID )]
names(GSE73663RNAseqMatInhMatch) <- ERchipAssociatedGENEIDsFinal
GSE73663RNAseqMatInhMatch <- GSE73663RNAseqMatInhMatch[-which(is.na(GSE73663RNAseqMatInhMatch))]
#get the matrix for only the genes common between the chip associated list and dataset
GSE73663InhcommInd <- match(GSE73663RNAseqMatInhMatch,rownames(GSE73663RNAseqMatInh))
GSE73663RNAseqMatInhComMat <- GSE73663RNAseqMatInh[GSE73663InhcommInd,]
rownames(GSE73663RNAseqMatInhComMat) <- names(GSE73663RNAseqMatInhMatch)
a <-  select(x = org.Hs.eg.db, keys = rownames(GSE73663RNAseqMatInh), columns = c("ENTREZID"), keytype = "ALIAS")
aa <- a$ALIAS[which(! is.na(a$ENTREZID))]
a <-  select(x = org.Hs.eg.db, keys = unique(aa), columns = c("ENTREZID"), keytype = "ALIAS")
aaa <- a$ALIAS[!duplicated(a$ALIAS)]
aaaa <- a$ENTREZID[!duplicated(a$ALIAS)]
aa <-  select(x = org.Hs.eg.db, keys = a$ALIAS[which(is.na(a$ENTREZID))] , columns = c("ENTREZID"), keytype = "ENSEMBL")
for (i in keytypes(org.Hs.eg.db)) {
  tryCatch({
    print(i)
    aa <-  select(x = org.Hs.eg.db, keys = c("AC063965.1") , columns = c("ENTREZID"), keytype = i)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for(i in keytypes(org.Hs.eg.db)){
  aa <-  select(x = org.Hs.eg.db, keys = c("AC063965") , columns = c("ENTREZID"), keytype = i)
}


#GSE73663RNAseqMatTC time course Merdak
mapSymbolToEntrezMerdakTC <- select(x = org.Hs.eg.db, keys = rownames(GSE73663RNAseqMatTC), columns = c("ENTREZID"), keytype = "SYMBOL")
GSE73663RNAseqMatTCMatch <- mapSymbolToEntrezMerdakTC$SYMBOL[match(ERchipAssociatedGENEIDsFinal,mapSymbolToEntrezMerdakTC$ENTREZID )]
names(GSE73663RNAseqMatTCMatch) <- ERchipAssociatedGENEIDsFinal
#there are 3255 NAs here, fix this later
GSE73663RNAseqMatTCMatch <- GSE73663RNAseqMatTCMatch[-which(is.na(GSE73663RNAseqMatTCMatch))]
#get the matrix for only the genes common between the chip associated list and dataset
GSE73663TCcommInd <- match(GSE73663RNAseqMatTCMatch,rownames(GSE73663RNAseqMatTC))
GSE73663RNAseqMatTCComMat <- GSE73663RNAseqMatTC[GSE73663TCcommInd,]
rownames(GSE73663RNAseqMatTCComMat) <- names(GSE73663RNAseqMatTCMatch)


length(Reduce(intersect, list(
                        # rownames(GSE73663RNAseqMatTCComMat)
                       rownames(GSE73663RNAseqMatInhComMat)
                       ,rownames(GSE86316RNAseqMatComMat)
                       ,rownames(GSE78167RNAseqMatCommon) 
                       ,rownames(GSE76507RNAseqMatNormComMat) 
                       ,rownames(GSE68358RNAseqMatNormComMat)
                       ,rownames(GSE67295RNAseqMatNormComMat)
                       ,rownames(GSE55922RNAseqNormComMat)
                       ,rownames(GSE53532RNAseqNormComMat) 
                       )))

DatasetMutualGenes <- Reduce(intersect, list(
  # rownames(GSE73663RNAseqMatTCComMat)
  rownames(GSE73663RNAseqMatInhComMat)
  ,rownames(GSE86316RNAseqMatComMat)
  ,rownames(GSE78167RNAseqMatCommon) 
  ,rownames(GSE76507RNAseqMatNormComMat) 
  ,rownames(GSE68358RNAseqMatNormComMat)
  ,rownames(GSE67295RNAseqMatNormComMat)
  ,rownames(GSE55922RNAseqNormComMat)
  ,rownames(GSE53532RNAseqNormComMat) 
))

DatasetMutualGenesKD <- Reduce(intersect, list(
  # rownames(GSE73663RNAseqMatTCComMat)
  rownames(GSE73663RNAseqMatInhComMat)
  ,rownames(GSE86316RNAseqMatComMat)
  #,rownames(GSE78167RNAseqMatCommon) 
  ,rownames(GSE76507RNAseqMatNormComMat) 
  ,rownames(GSE68358RNAseqMatNormComMat)
  ,rownames(GSE67295RNAseqMatNormComMat)
  ,rownames(GSE55922RNAseqNormComMat)
  ,rownames(GSE53532RNAseqNormComMat) 
))
#length(unique(rownames(GSE53532RNAseqNormComMat)))
#removing the batch effects
# plot(prcomp(t(cbind(mutualGSE73663RNA,mutualGSE86316RNA,mutualGSE78167RNA))))
# biplot(prcomp(t(cbind(mutualGSE73663RNA,mutualGSE86316RNA,mutualGSE78167RNA))))
resize.win()
res.pca <- prcomp(t(combat_edata), scale = F)
ind.coord <- res.pca$x
plot(ind.coord[,1], ind.coord[,2], pch = 19,  
     xlab="PC1",ylab="PC2",ylim = range(res.pca2$x[,2]), xlim = range(res.pca2$x[,1]))
abline(h=0, v=0, lty = 2)
text(ind.coord[,1], ind.coord[,2], labels=rownames(ind.coord),
     cex=0.7, pos = 3)
#
resize.win()
res.pca2 <- prcomp(t(mutualAll7KDdatasets), scale = F)
ind.coord <- res.pca2$x
plot(ind.coord[,1], ind.coord[,2], pch = 19,  
     xlab="PC1",ylab="PC2")
abline(h=0, v=0, lty = 2)
text(ind.coord[,1], ind.coord[,2], labels=rownames(ind.coord),
     cex=0.7, pos = 3)
#warnings()
mutualGSE73663RNA <- GSE73663RNAseqMatInhComMat[match(DatasetMutualGenesKD,rownames(GSE73663RNAseqMatInhComMat)),]
mutualGSE86316RNA <- GSE86316RNAseqMatComMat[match(DatasetMutualGenesKD,rownames(GSE86316RNAseqMatComMat)),]
mutualGSE78167RNA <- GSE78167RNAseqMatCommon[match(DatasetMutualGenesKD,rownames(GSE78167RNAseqMatCommon)),]
mutualGSE53532RNA <- GSE53532RNAseqComMat[match(DatasetMutualGenesKD,rownames(GSE53532RNAseqComMat)),] #CTCF kd
mutualGSE55922RNA <- GSE55922RNAseqComMat[match(DatasetMutualGenesKD,rownames(GSE55922RNAseqComMat)),] 
mutualGSE67295RNA <- GSE67295RNAseqMatComMat[match(DatasetMutualGenesKD,rownames(GSE67295RNAseqMatComMat)),] 
mutualGSE68358RNA <- GSE68358RNAseqMatComMat[match(DatasetMutualGenesKD,rownames(GSE68358RNAseqMatComMat)),] 
mutualGSE76507RNA <- GSE76507RNAseqMatComMat[match(DatasetMutualGenesKD,rownames(GSE76507RNAseqMatComMat)),] 

mutualAll7KDdatasets <- cbind(mutualGSE73663RNA #1 Zeynep     ctrl: 11,12 
                            ,mutualGSE86316RNA #2 Men1 kd   ctrl: 1,2
                            #,mutualGSE78167RNA #3 time points 0,1,2,3,...     ctrl: 1, 11, 21
                            ,mutualGSE53532RNA #3 CTCF kd     ctrl: 1
                            ,mutualGSE55922RNA #4 BRD4 kd     cntrl: 5,6
                            ,mutualGSE67295RNA #5 #E2 and other treatments ctrl : 11,12
                            ,mutualGSE68358RNA#6 #E2 or E2 plus something ctrl : 1,5,9,10,11,14,18,24
                            ,mutualGSE76507RNA) #7 G9 and PHF20 kd   ctrl: 3,5,6
#ncol(mutualGSE76507RNA)
#construct the mutual control matrix
mutualControlMat <- cbind(
  mutualGSE73663RNA[,c(11,12)]
  ,mutualGSE86316RNA[,c(1,2)]
  ,mutualGSE78167RNA[,c(1,11,21)]
  ,mutualGSE53532RNA[,1]
  ,mutualGSE55922RNA[,c(5,6)]
  ,mutualGSE67295RNA[,c(11,12)]
  ,mutualGSE68358RNA[,c(1,5,9,10,11,14,18,24)]
  ,mutualGSE76507RNA[,c(3,5,6)]
  )
colnames(mutualAll8datasets) <- c(rep(1,18),rep(2,8),rep(3,30),rep(4,4),rep(5,12),rep(6,26),rep(7,24),rep(8,12))
colnames(mutualAll7KDdatasets) <- c(rep(1,18),rep(2,8),rep(3,4),rep(4,12),rep(5,26),rep(6,24),rep(7,12))

#GSE76507RNAseqMatComMat[1,1]
disEuc <- dist(t(combat_edata))
fitEuc <- cmdscale(disEuc,eig=TRUE, k=2)
fitEuc
# plot solution 
x <- fitEuc$points[,1]
y  <- fitEuc$points[,2]
z <- fitEuc$points[,3]
colCanc <- as.numeric(colnames(mutualAll8datasets))
pchcol<- as.numeric(colnames(mutualAll8datasets))
plot(x,y,col =col_vector[colCanc],pch = pchcol ,main = "Metric MDS",xaxt= "n",yaxt = "n",ylab = "Coordinate 2",xlab = "Coordinate 1",cex = 1.5,cex.lab = 1.4,cex.main = 1.4)

source("https://bioconductor.org/biocLite.R")
biocLite("sva")
library(sva)
# biocLite("bladderbatch")
# library(bladderbatch)
# data(bladderdata)
# pheno = pData(bladderEset)
phenoER <- as.data.frame(cbind(c(1:104),as.numeric(colnames(mutualAll7KDdatasets))))
#rownames(phenoER) <- colnames(mutualAll8datasets)
colnames(phenoER) <- c("sample","batch")
batch = phenoER$batch
edata = mutualAll7KDdatasets
modcombat = model.matrix(~1, data=phenoER)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat,prior.plots=F)


install.packages("gPCA")
library(gPCA)
batchDetect2 <- gPCA.batchdetect(t(combat_edata), batch = batch, filt = NULL, nperm = 1000, center = FALSE, scaleY=FALSE, seed = NULL)
batchDetect2$p.val
gDist(batchDetect)


resize.win()
PCplot(batchDetect2, ug = "unguided", type = "comp", npcs = 2)

#Draw heatmaps for each RNA-seq experiment
par(mfrow = c(1,1),mar = c(8,8,8,8),pin = c(5,5))
heatmap(GSE76507RNAseqMatNormComMat,key = T)
dev.off()
resize.win <- function(Width=6, Height=6)
{
  # works for windows
  #dev.off(); # dev.new(width=6, height=6)
  quartz( width=Width, height=Height)
}
resize.win(8,8)
#dev.new <- function(width = 7, height = 7) 
heatmap.2(GSE76507RNAseqMatNormComMat,key = T)
par(mar = c(1,1,1,1))
# library(lattice)
# levelplot(GSE53532RNAseqNormComMat)
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

heatmap.2(GSE67295RNAseqMatNormComMat,
          #cellnote = GSE76507RNAseqMatNormComMat,  # same data set for cell labels
          #main = " G9 and PHF kd", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(6,6),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
         # breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",
         scale = "row")            # turn off column clustering
resize.win(8,8)
combat_edataQuant <- normalize.quantiles(combat_edata)
combat_edataColScaleQuantRowScale <- t(scale(t(normalize.quantiles(scale(combat_edata)))))
combat_edataColScaleQuant <- normalize.quantiles(scale(combat_edata))
combat_edataQuantScale <- scale(combat_edataQuant)
combat_edataQuantRowscale <- t(scale(t(combat_edataQuant)))
heatmap.2(combat_edataColScaleQuant,col=my_palette,density.info="none",trace="none",Colv = F,dendrogram = "row")  

heatmap.2(GSE55922RNAseqNormComMat,col=my_palette,  scale = "row",density.info="none",trace="none")    
heatmap.2(GSE67295RNAseqMatNormComMat,col=my_palette,  scale = "row",density.info="none",trace="none")    
heatmap.2(GSE68358RNAseqMatNormComMat,col=my_palette,  scale = "row",density.info="none",trace="none")    
heatmap.2(GSE76507RNAseqMatNormComMat,col=my_palette,  scale = "row",density.info="none",trace="none")    


######
read.table("Expression data/DiffExpNames.docx")

####another approach to normalization: look at the correlation between all control (WT with no stimulation)
#which of the 134 conditions are control experiments?


