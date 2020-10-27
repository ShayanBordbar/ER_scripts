# GRO-seq enhancer transcription
# find out which enhancers are transcribed after estrogen treatment

#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         Libraries         ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
library(rtracklayer)

#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         Functions         ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
# list of all functions: 1. RegulatoryCoordinateExtractor
#########################################################################################################
#########################################################################################################
RegulatoryCoordinateExtractor <- function(enhancer.coord=MCFEnhancersGR,
                                          promoter.coor=Promoter1kbGR.gene,
                                          enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                          promoter.promoter.int=Promoter.gene.PromoterIntList,
                                          gene.Names){
  # This function gets as input:
  # 1. enhancer.coord : the genomic coordinates of all enhancers 
  # 2. promoter.coor : the genomic coordinates of all promoters
  # 3. enhancer.per.gene : a list containing list of enhancers for each gene eg: Vicinity100kb.Enhancer.plusPromoter.By.gene
  # 4. promoter.promoter.int : a list containing list of promoters interacting with each gene eg: Promoter.gene.PromoterIntList
  # 5. gene.Names : names of the genes of interest (should be the same as the names used in the privious four inputs)
  # Outputs: a dataframe containing the coordinates of regulatory elements of all genes, rownames: name of the gene + index of the regulatory element
  
  reg.coor.list <- list()
  for(gene in 1:length(gene.Names)){
    gene.Orig.Ind <- match(gene.Names[gene], names(promoter.coor))
    cur.prom <- as.data.frame(promoter.coor[gene.Names[gene]]) #get the coordinates of the promoter
    #if there are any interacting promoters get their coordinates
    
    if (length(promoter.promoter.int[[gene.Orig.Ind]]) > 0){
      cur.prom.int = as.data.frame(promoter.coor[promoter.promoter.int[[gene.Orig.Ind]]])
      cur.prom <- rbind(cur.prom, cur.prom.int)
    }
    # add enhancers if there are any
    if (length(enhancer.per.gene[[gene.Orig.Ind]]) > 0){
      cur.enh <- as.data.frame(enhancer.coord[enhancer.per.gene[[gene.Orig.Ind]]])
      cur.all <- rbind(cur.prom, cur.enh)
    }else{
      cur.all <- cur.prom
    }
    rownames(cur.all) <- paste(gene.Names[gene], c(1:nrow(cur.all)), sep = "_")
    reg.coor.list[[gene]] <- cur.all
  }
  reg.coor.df <- do.call(rbind, reg.coor.list)
  return(reg.coor.df)
}

#########################################################################################################
#########################################################################################################
# example
aa <- RegulatoryCoordinateExtractor(gene.Names = rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2))
#########################################################################################################
#########################################################################################################
# next function starts here


#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         GSE67295          ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################

GSE67295.GROseq <- read.delim("Expression data/GRO-seq/GSE67295GeneExpressionGroSeq.txt",header = T)
colnames(GSE67295.GROseq)[1] <- "Transcript ID"
colnames(GSE67295.GROseq)[2] <- "Chr"
colnames(GSE67295.GROseq)[3] <- "Start"
colnames(GSE67295.GROseq)[4] <- "End"
GSE67295.GROseq.GRanges <- makeGRangesFromDataFrame(GSE67295.GROseq[,2:4])
Exp17.RegulatoryElements.Coordiantes <- makeGRangesFromDataFrame(RegulatoryCoordinateExtractor(gene.Names = rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2)))
aa <- overlapExtractor(Exp17.RegulatoryElements.Coordiantes, GSE67295.GROseq.GRanges)

aaa <- GSE67295.GROseq[unique(unlist(aa)), ]
  
#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         GSE27463          ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################

GSE27463.GROseq <- read.table("Expression data/GRO-seq/GSE27463_RefSeq.all.tsv",header = T)

# find the intersection of enhancers used in ensemble experiment 17 and the ones 

GSE27463.GROseq.GRanges <- makeGRangesFromDataFrame(GSE27463.GROseq[, 1:3])
Exp17.RegulatoryElements.Coordiantes <- makeGRangesFromDataFrame(RegulatoryCoordinateExtractor(gene.Names = rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2)))
aa <- overlapExtractor(Exp17.RegulatoryElements.Coordiantes, GSE27463.GROseq.GRanges)


sum(unlist(lapply(aa, length)) > 0)
unique(unlist(aa))
a <- RegulatoryCoordinateExtractor(gene.Names = rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed.Filtered2))
aaa <- GSE27463.GROseq[unique(unlist(aa)), ]
aaaa <- aaa[aaa$Cluster != "NR",]


aa[[1]]
