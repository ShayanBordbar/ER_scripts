#ER Project:
#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         Directory         ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
setwd("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor")

#Studying transcriptional activity of Estrogen Receptor in MCF7 cell line
#Wrapping up what I have so far:
#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         Libraries         ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
library(rgl)
library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(IRanges)
library(ChIPseeker)
library(stringi)
library(stringr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(sva)
library(gPCA)
library(biomaRt)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19")

library(MotifDb)
library (seqLogo)
library(PWMEnrich)
library(preprocessCore)
library(Biobase)
library(GEOquery)
library(beadarray)
library(illuminaHumanv4.db)
library(hexbin)
library(limma)
library(oligo)
library(affycoretools)
library(clariomdhumantranscriptcluster.db)
library("GEOmetadb")
library("hgu133plus2.db")
library(hgu133a.db)
library("hgu133a2.db")
library(psych)


#install.packages("~/Downloads/creDB")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         Functions         ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
# list of functions in this script:
# 1) StartEndCoordinateCheck , 2) InteractionExtractor , 3) overlapExtractor , 4) InteractionExtender
# 5) ChipPeakInvestigator , 6) OverlapInteractionBatchExtract , 7) geneNameToEntrez , 8)getCommonGenes
# 9) Clustrator 10) clustratorWrapper 11) WriteFastaOfBag 12) AddpsudoRow , 13) RemoveDegenRow 14) PWMtoCount
# 15) MotifWriter 16) ExpressionWriter 17) LogFoldTransform 18) Hal_job_writer 19) bash_directory 20) GEMSTAT_input_constructor
# 21) createDescreteMat 22) CheckAndReplace 23) DiffExpToTopGeneMat 24) CommonDifExpMat 25) QuantileCompare
# 26) TFexpressionDescritizer 27) plotExpression 28) EnhancerChopper_rangeMat 
# 29) EnhancerChopper 30)EnhancerChopper_wrapper 31)my_dist_fuc 32)mergeConnectedRanges

#########################################################################################################
#########################################################################################################
StartEndCoordinateCheck <- function(inputRangedf, ChrLengths=ChromosomeLengthHG38){
  # write a function to check the start and end of a range dataframe, so they are not negative , or larger than chromosome length
  # ChromosomeLengthHG19 <- read.delim(file="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/ChromosomeLengthHG19.txt",
  #                                    header=F,
  #                                    stringsAsFactors = F)
  # ChromosomeLengthHG38 <- read.delim(file="/Users/Shayan/Documents/Shayan/BioInf/hg38_chromosome_sizes.csv",
  #                                    header=F,
  #                                    stringsAsFactors = F, sep = ",")
  # inputRangedf is a dataframe where first column is name of
  # the chromosome (eg chr1), second column is start coordinate, third column is end coordinate
  Newdf <- inputRangedf
  for(i in 1: nrow(ChrLengths)){
    cat(paste0("chr", as.character(i), "\n"))
    currInd <- which(inputRangedf[ , 1] == ChrLengths[i, 1] )
    # Finding the index of the ones larger than chromosome length
    maxOverInd <- currInd[which(inputRangedf[currInd, 3] > ChrLengths[i, 2])]
    cat(paste0("maxOver ", as.character(length(maxOverInd)), "\n"))
    if(length(maxOverInd) > 0){
      # setting the overMax coordinate to max
      Newdf[maxOverInd, 3] <- ChrLengths[i, 2]
    }
    # Finding the index of the ones smaller than zero
    minNegInd  <- currInd[which(inputRangedf[currInd, 2] < 0)]
    if(length(minNegInd) > 0){
      # setting the negative coordinate to one
      Newdf[minNegInd, 2] <- 1
    }
    cat(paste0("minNeg ", as.character(length(minNegInd)), "\n"))
  }
  return(Newdf)
}
#########################################################################################################
#########################################################################################################
#example
aaa <- StartEndCoordinateCheck(inputRangedf = tssdf)
#########################################################################################################
#########################################################################################################
perentOverlapFilter <- function(queryGR, subjectGR, hits, percent){
  # queryGR : GRanges object containing the coordinates of the query sequences
  # subjectGR : GRanges object containing the coordinates of the subject sequences
  # hits : output of findOverlaps() function.
  # percent : threshold on what percent the query sequence to overlap with subject sequence
  # This function modifies the input hits object (output of findOverlap function), and removes all hits with less percent overlap than the specified threshold
  overlaps <- pintersect(queryGR[queryHits(hits)],
                         subjectGR[subjectHits(hits)])
  # percent overlap with the query
  percentOverlap <- (width(overlaps) / width(queryGR[queryHits(hits)])) * 100
  # percent overlap with the subject
  percentOverlap2 <- (width(overlaps) / width(subjectGR[subjectHits(hits)])) * 100
  my_cond <- (percentOverlap > percent) | (percentOverlap2 > percent)
  my_hits <- hits[my_cond]
  return(my_hits)
}
#########################################################################################################
#########################################################################################################
aaa <- findOverlaps(makeGRangesFromDataFrame(Promoter1kbdf.gene), MCFEnhancersGR)
aaaa <- perentOverlapFilter(queryGR = makeGRangesFromDataFrame(Promoter1kbdf.gene),
                            subjectGR = MCFEnhancersGR,
                            hits = aaa,
                            percent = 50)
#########################################################################################################
#########################################################################################################
## Extract clusters from Hits object.
extractClustersFromSelfHits <- function(hits)
{
  stopifnot(is(hits, "Hits"))
  stopifnot(queryLength(hits) == subjectLength(hits))
  hits <- union(hits, t(hits))
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  cid <- seq_len(queryLength(hits))  # cluster ids
  while (TRUE) {
    h <- Hits(qh, cid[sh],
              queryLength(hits), subjectLength(hits))
    cid2 <- pmin(cid, selectHits(h, "first"))
    if (identical(cid2, cid))
      break
    cid <- cid2
  }
  unname(splitAsList(seq_len(queryLength(hits)), cid))
}
## Merge ranges that are "connected" (directly or indirectly)
## via a hit (or several hits) in 'hits'.
mergeConnectedRanges <- function(x, hits)
{
  # x is a GRanges object, 
  #hits is a Hits object where query and subject are both x[it can be a subset of the whole thing]
  stopifnot(is(x, "GenomicRanges"))
  stopifnot(is(hits, "Hits"))
  stopifnot(queryLength(hits) == subjectLength(hits))
  stopifnot(queryLength(hits) == length(x))
  clusters <- extractClustersFromSelfHits(hits)
  ans <- range(extractList(x, clusters))
  if (any(sapply(ans, length) != 1L))
    stop(wmsg("some connected ranges are not on the same ",
              "chromosome and strand, and thus cannot be ",
              "merged"))
  ans <- unlist(ans)
  mcols(ans)$revmap <- clusters
  ans
}
#aad <- mergeConnectedRanges(, aa)

#########################################################################################################
#########################################################################################################
InteractionExtractor <- function(interactingA=Genome4DMCF7[,1:3],
                                 interactingB=Genome4DMCF7[,4:6],
                                 subjectCoord1,
                                 subjectCoord2,
                                 ..minOverlap = 1,
                                 .minOverLapPercent=0,
                                 int_piece_max_length = integer(0)){
  
  # Write a function that gets the coordinates of two sides of
  # interaction + another set of coordinates (e.g. enhancers coordinates) and
  # creates a list of interacting coordinates for each coordinate based on the interaction data 
  
  # interactingA is the coordinates of the first side of interaction matrix as
  # a dataframe with 3 columns. example:Genome4DMCF7[,1:3]
  # interactingB is the coordinates of the second side of interaction matrix as
  # a dataframe with 3 columns. example:Genome4DMCF7[,4:6]
  # subjectCoord1 is the coordinates of the first subject examined 
  # subjectCoord2 is the coordinates of the second subject examined
  # Output is a list of length = nrow(subjectCoord1), each entry is a list containing
  # the index of interacting elements of subjectCoord2 with this element
  # ..minOverlap is the minimum length of overlap between subjectCoord(1 or 2) and interacting(A or B) sequences
  # .minOverLapPercent : is the percent of the length of the query sequence that the subject sequence has to be in overlap with, to be considred a hit
  # int_piece_max_length : maximum length of pieces in the interaction matrix provided
  
  print("extract interactions")
  # create a GRange object for each dataframe
  if (typeof(interactingA) != "S4"){
    InteractingAGR <- makeGRangesFromDataFrame(interactingA)
  }else{
    InteractingAGR <- interactingA
  }
  if (typeof(interactingB) != "S4"){
    InteractingBGR <- makeGRangesFromDataFrame(interactingB)
  }else{
    InteractingBGR <- interactingB
  }
  if (typeof(subjectCoord1) != "S4"){
    subjectCoord1GR <- makeGRangesFromDataFrame(subjectCoord1)
  }else{
    subjectCoord1GR <- subjectCoord1
  }
  if (typeof(subjectCoord2) != "S4"){
    subjectCoord2GR <- makeGRangesFromDataFrame(subjectCoord2)
  }else{
    subjectCoord2GR <- subjectCoord2
  }
  print("past format")
  # filtering the length of the interacting matrix pieces
  if(length(int_piece_max_length) == 0){
    int_piece_max_length <- max(max(unlist(width(InteractingAGR))),
                                max(unlist(width(InteractingBGR))))
  }
  int_a_len <- unlist(width(InteractingAGR))
  int_b_len <- unlist(width(InteractingBGR))
  int_a_b_log <- (int_a_len <= int_piece_max_length) | (int_b_len <= int_piece_max_length)
  InteractingAGR <- InteractingAGR[int_a_b_log]
  InteractingBGR <- InteractingBGR[int_a_b_log]
  
  # finding overlap hits of each subject in each interaction side
  print("findOverlaps_S1_A")
  OverlapSub1IntA <- findOverlaps(query=subjectCoord1GR, 
                                  subject=InteractingAGR,
                                  minoverlap = ..minOverlap)
  print("filterOverlaps_S1_A")
  OverlapSub1IntA <- perentOverlapFilter(queryGR=subjectCoord1GR ,
                                         subjectGR=InteractingAGR ,
                                         hits=OverlapSub1IntA ,
                                         percent=.minOverLapPercent)
  print("findOverlaps_S1_B")
  OverlapSub1IntB <- findOverlaps(query=subjectCoord1GR,
                                  subject=InteractingBGR,
                                  minoverlap = ..minOverlap)
  print("filterOverlaps_S1_B")
  OverlapSub1IntB <- perentOverlapFilter(queryGR=subjectCoord1GR,
                                         subjectGR=InteractingBGR,
                                         hits=OverlapSub1IntB,
                                         percent=.minOverLapPercent )
  print("findOverlaps_S2_A")
  OverlapSub2IntA <- findOverlaps(query=subjectCoord2GR,
                                  subject=InteractingAGR, 
                                  minoverlap = ..minOverlap)
  print("filterOverlaps_S2_A")
  OverlapSub2IntA <- perentOverlapFilter(queryGR=subjectCoord2GR,
                                         subjectGR=InteractingAGR,
                                         hits=OverlapSub2IntA ,
                                         percent=.minOverLapPercent )
  print("findOverlaps_S2_B")
  OverlapSub2IntB <- findOverlaps(query=subjectCoord2GR, 
                                  subject=InteractingBGR,
                                  minoverlap = ..minOverlap)
  print("filterOverlaps_S2_B")
  OverlapSub2IntB <- perentOverlapFilter(queryGR=subjectCoord2GR,
                                         subjectGR=InteractingBGR,
                                         hits=OverlapSub2IntB ,
                                         percent=.minOverLapPercent )
  
  InteractionList <- list()
  for (i in 1:length(subjectCoord1GR )){
    currSubIntA <- OverlapSub1IntA@to[OverlapSub1IntA@from == i]
    # create a numeric vector containing the index of all elements of
    # subject 2 which overlap with the other side of the interaction that
    # overlaps with i th element of subject1
    currSub1Sub2A <- numeric(0)
    if (length(currSubIntA) > 0){
      for(j in 1:length(currSubIntA)){
        currSub1Sub2A <- c(currSub1Sub2A, OverlapSub2IntB@from[OverlapSub2IntB@to == currSubIntA[j]])
      }
    }
    currSubIntB <- OverlapSub1IntB@to[OverlapSub1IntB@from == i]
    # create a numeric vector containing the index of all elements of subject 2 which
    # overlap with the other side of the interaction that overlaps with i th element of subject1
    currSub1Sub2B <- numeric(0)
    if (length(currSubIntB) > 0){
      for(j in 1:length(currSubIntB)){
        currSub1Sub2B <- c(currSub1Sub2B,OverlapSub2IntA@from[OverlapSub2IntA@to == currSubIntB[j]])
      }
    }
    #print("aabeforeuniq")
    InteractionList[[i]] <- union(currSub1Sub2A, currSub1Sub2B)
    #print("aaafteruniq")
  }
  return(InteractionList)
}
#########################################################################################################
#########################################################################################################
#example
aaa <- InteractionExtractor(subjectCoord1=MCFEnhancersdf, subjectCoord2=MCFEnhancersdf)
hist(unlist(lapply(aaa,length)))
hist(sumColumns)
sum(sumColumns == unlist(lapply(aaa,length)))


#########################################################################################################
#########################################################################################################
overlapExtractor <- function(subjectCoord1, subjectCoord2, .minOverlap = 1, minOverlapPercent=0){
  #write a function to get two dataframe or GRanges and output a list of
  #length = nrow(first entry). each entry is the index of entries in
  #the second input which overlap with that element.
  # .minOverlap is the minimum number of positions to overlap
  # minOverlapPercent : is the percent of the length of the query sequence that the subject sequence has to be in overlap with, to be considred a hit
  
  #subjectCoord1 is the coordinates of the first element as a dataframe or a GRange object
  if (typeof(subjectCoord1) != "S4"){
    subjectCoord1 <- makeGRangesFromDataFrame(subjectCoord1)
  }
  if (typeof(subjectCoord2) != "S4"){
    subjectCoord2 <- makeGRangesFromDataFrame(subjectCoord2)
  }
  sub1sub2Overlap <- findOverlaps(query=subjectCoord1, subject=subjectCoord2, minoverlap = .minOverlap)
  sub1sub2Overlap <- perentOverlapFilter(queryGR=subjectCoord1,
                                         subjectGR=subjectCoord2,
                                         hits=sub1sub2Overlap ,
                                         percent=minOverlapPercent)
  sub1sub2OverlapList <- list()
  for(i in 1:length(subjectCoord1)){
    sub1sub2OverlapList[[i]] <- sub1sub2Overlap@to[sub1sub2Overlap@from == i]
  }
  return(sub1sub2OverlapList)
}
#########################################################################################################
#########################################################################################################
#example
aaa <- overlapExtractor(subjectCoord1=Promoter1kbdf.gene, subjectCoord2=MCFEnhancersGR)
#########################################################################################################
#########################################################################################################

InteractionExtender <- function(InteractionList,
                                Subject,
                                interactingA=Genome4DMCF7[,1:3],
                                interactingB=Genome4DMCF7[,4:6],
                                ...minOverlap = 1,
                                ..minOverlapPercent=0,
                                .int_piece_max_length = integer(0)){
  #write a function that gets 
  # InteractionList : a list like the output of overlapExtractor function 
  # Subject : coordinates of the subject that the index in the list is refering to 
  # interactingA = Genome4DMCF7[,1:3],
  # interactingB = Genome4DMCF7[,4:6]
  # ...minOverlap minimum number of overlapping bases
  # ..minOverlapPercent : is the percent of the length of the query sequence that the subject sequence has to be in overlap with, to be considred a hit
  # .int_piece_max_length : maximum length of pieces in the interaction matrix provided
  #outputs a list that for each entry of the input list, adds the indicis of the elements in the subject list that interact with it
  
  #subject is a GRange or dataframe object containing the coordinates of the subject that we want to look at its interactions
  if (typeof(Subject) != "S4"){
    Subject <- makeGRangesFromDataFrame(Subject)
  }
  ExtendedInteractionList <- list()
  Extentions <- list()
  SubjectSelfInteractions <- InteractionExtractor(interactingA=interactingA,
                                                  interactingB=interactingB,
                                                  subjectCoord1=Subject,
                                                  subjectCoord2= Subject,
                                                  ..minOverlap= ...minOverlap,
                                                  .minOverLapPercent=..minOverlapPercent,
                                                  int_piece_max_length = .int_piece_max_length)
  
  for(i in 1:length(InteractionList)){
    if(length(InteractionList[[i]]) > 0){
      NewInt <- SubjectSelfInteractions[InteractionList[[i]]]
      ExtendedInteractionList[[i]] <- unique(c(InteractionList[[i]],unlist(NewInt)))
      Extentions[[i]] <- unique(unlist(NewInt))
    }else{
      ExtendedInteractionList[[i]] <- numeric(0)
      Extentions[[i]] <- numeric(0)
    }
  }
  return(list(Merged =ExtendedInteractionList,Extentions = Extentions))
}
#########################################################################################################
#########################################################################################################
#example
aaa <- InteractionExtender(InteractionList = Promoter.gene.EnhancerIntList,Subject = MCFEnhancersGR)
#########################################################################################################
#########################################################################################################
ChipPeakInvestigator <- function(promoter.promoter = Promoter.gene.PromoterIntList,
                                 gene.enhancer = Vicinity20kb.Enhancer.plusPromoter.By.gene.Extended,
                                 enhancerChipBinary = Enhancer.ERchip.interORover.byEnhancer.Binary,
                                 promoterChipBinary = Promoter.gene.ERchip.interORover.byPromoter.Binary
                                 #,regElemBasedList = F
                                 ){
  #write a function that gets an overlap or interaction list for  a gene-enhancer or a transcript-enhancer (e.g. output of overlapExtractor),
  # an overlap or interaction list for  a gene promoter-gene promoter-gene or a promoter-transcript promoter-transcript (e.g. output of overlapExtractor),
  ##the first two inputs are lists of the same length
  #the binary matrix of enhancer, chip interaction-overlap
  #the binary matrix of promoter, chip interaction-overlap
  #outputs a binary matrix with the nrow = number of genes or transcripts , ncol = number of chip datasets which is 1 if any of the enhancers or promoters associated with that gene or transcript overlap with at least one chip peak in that data set
  #outputs a non-binary version of the above matrix
  # SKIPED FOR NOWif regElemBasedList is True outputs a list with two entries: one entry for promoters, one entry for enhancers, each entry contains list of genes associated with each reg element which has a one in the binary matrix
  AssociatedPeaksBinary <- matrix(0L, nrow=length(promoter.promoter), ncol=ncol(enhancerChipBinary))
  AssociatedPeaksTotal  <- matrix(0L, nrow=length(promoter.promoter), ncol=ncol(enhancerChipBinary))
  for(i in 1:length(promoter.promoter)){
    curPromChip <- promoterChipBinary[unique(c(i,promoter.promoter[[i]])),]
    curEnhChip  <- enhancerChipBinary[unique(gene.enhancer[[i]]),]
    if(ncol(enhancerChipBinary) == 1){
      TotalChip = sum(curPromChip) + sum(curEnhChip)
    }else{
      TotalChip <- colSums(rbind(curPromChip,curEnhChip))
    }
    AssociatedPeaksTotal[i,] <- TotalChip
    AssociatedPeaksBinary[i,(TotalChip > 0)] <- 1
  }
  if(!is.null(names(gene.enhancer))){
    rownames(AssociatedPeaksTotal) <- names(gene.enhancer)
    rownames(AssociatedPeaksBinary) <- names(gene.enhancer)
  }
  # if(regElemBasedList){
  #   regElemBList <- list()
  #   regElemBList[[1]] <- list()
  #   regElemBList[[2]] <- list()
  #   names(regElemBList) <- c("promoters","enhancers")
  #   actproall <- which(rowSums(promoterChipBinary) > 0)
  #   if(length(actproall) > 0){
  #     for(actpro in 1:length(actproall)){
  #       rownames(actproall[actpro])
  #     }
  #   }
  # 
  #   for(actenh in 1:nrow(enhancerChipBinary)){
  #     
  #   }
  # }
  return(list(AssociatedPeaks = AssociatedPeaksTotal,AssociatedPeaksBinary = AssociatedPeaksBinary))
}
#########################################################################################################
#########################################################################################################
#example
aaa <- ChipPeakInvestigator(promoter.promoter=Promoter.gene.PromoterIntList,
                            gene.enhancer=Vicinity20kb.Enhancer.plusPromoter.By.gene.Extended,
                            enhancerChipBinary = Enhancer.ERchip.interORover.byEnhancer.Binary,
                            promoterChipBinary = Promoter.gene.ERchip.interORover.byPromoter.Binary)
#########################################################################################################
#########################################################################################################
OverlapInteractionBatchExtract <- function(interactingA=Genome4DMCF7[,1:3],
                                           interactingB=Genome4DMCF7[,4:6],
                                           subjectCoord1,
                                           subjectCoord2List,
                                           mode = "both",
                                           ...minOverlap = 1,
                                           ...minOverlapPercent=0,
                                           .int_piece_max_length = integer(0)){
  #write a function that gets as input:
  # 1. Enhancer or promoter coordinates (Granges or dataframe)
  # 2. A list of ChIP coordinates
  # 3. interactingA = Genome4DMCF7[,1:3],
  # 4. interactingB = Genome4DMCF7[,4:6]
  # 5. input mode: can be "both","overlap" or "interaction". if mode is overlap: only computes overlap returns the first two outputs, same for interaction
  # 6. ...minOverlap : is the minimum number of overlaping bases
  # 7. ...minOverlapPercent : is the percent of the length of the query sequence that the subject sequence has to be in overlap with, to be considred a hit
  # 8. .int_piece_max_length : maximum length of pieces in the interaction matrix provided
  #outputs a list of Four entires:
  #First one is the overlap list of Input 1 vs input 2
  #Second one is the overlap matrix of nrow(Input1) vs ncol(length of the list input 2) which is the number of overlaps of each enhancer/promoter with the set of ChIP peaks
  #Third is the same as First but for interactions
  #Forth is the same as Second but for interactions
  
  if (typeof(interactingA) != "S4"){
    interactingA <- makeGRangesFromDataFrame(interactingA)
  }
  if (typeof(interactingB) != "S4"){
    interactingB <- makeGRangesFromDataFrame(interactingB)
  }
  if (typeof(subjectCoord1) != "S4"){
    subjectCoord1 <- makeGRangesFromDataFrame(subjectCoord1)
  }
  if (typeof(subjectCoord2List[[1]]) != "S4"){
    for(i in 1:length(subjectCoord2List)){
      subjectCoord2List[[i]] <-  makeGRangesFromDataFrame(subjectCoord2List[[i]])
    }
  }
  if (mode != "interaction"){
    print("calculating overlaps")
    #Get the "overlap" of all "enhancers" with REMAP chip peaks
    Subject1.chip.Overlap.bySubject1 <- list()#a list where each entry is the overlap list of the enhancers with one chip dataset
    for(i in 1:length(subjectCoord2List)){
      Subject1.chip.Overlap.bySubject1[[i]] <- overlapExtractor(subjectCoord1=subjectCoord1,
                                                                subjectCoord2=subjectCoord2List[[i]],
                                                                .minOverlap = ...minOverlap,
                                                                minOverlapPercent = ...minOverlapPercent)
      print(i)
    }
    if (length(names(subjectCoord2List)) > 0){
      names(Subject1.chip.Overlap.bySubject1) <- names(subjectCoord2List)
    }
    #Create a  matrix where each row represents an enhancer, each column represents  a ChIPed factor: entry[i,j] is the number of Chip-peaks of factor j that overlap with the enhancer i
    Subject1.chip.Overlap.bySubject1.Mat <- matrix(0L,
                                                   nrow = length(subjectCoord1),
                                                   ncol = length(Subject1.chip.Overlap.bySubject1))
    for(i in 1:ncol(Subject1.chip.Overlap.bySubject1.Mat)){
      Subject1.chip.Overlap.bySubject1.Mat[ ,i] <- unlist(lapply(Subject1.chip.Overlap.bySubject1[[i]],
                                                                 length))
    }
    if (length(names(subjectCoord2List)) > 0){
      colnames(Subject1.chip.Overlap.bySubject1.Mat) <- names(subjectCoord2List)
    }
  }
  #########################################################################################################
  if(mode != "overlap"){
    print("calculating interactions")
    #Get the "interaction" of  "subject1" with each element of subject2list
    Subject1.chip.Interaction.bySubject1 <- list()#a list where each entry is the interaction list of the enhancers with one chip dataset
    for(i in 1:length(subjectCoord2List)){
      Subject1.chip.Interaction.bySubject1[[i]] <- InteractionExtractor(subjectCoord1=subjectCoord1,
                                                                        subjectCoord2=subjectCoord2List[[i]],
                                                                        interactingA=interactingA,
                                                                        interactingB=interactingB,
                                                                        ..minOverlap = ...minOverlap,
                                                                        .minOverLapPercent = ...minOverlapPercent,
                                                                        int_piece_max_length = .int_piece_max_length)
      print(i)
    }
    if (length(names(subjectCoord2List)) > 0){
      names(Subject1.chip.Interaction.bySubject1) <- names(subjectCoord2List)
    }
    #Create a  matrix where each row represents one entry of subjectcoord1, each column represents  an entry of the subject2list: entry[i,j] is the number of Chip-peaks of factor j that overlap with the enhancer i
    Subject1.chip.Interaction.bySubject1.Mat <- matrix(0L,
                                                       nrow = length(subjectCoord1),
                                                       ncol = length(Subject1.chip.Interaction.bySubject1))
    for(i in 1:ncol(Subject1.chip.Interaction.bySubject1.Mat)){
      Subject1.chip.Interaction.bySubject1.Mat[ ,i] <- unlist(lapply(Subject1.chip.Interaction.bySubject1[[i]],length))
    }
    if (length(names(subjectCoord2List)) > 0){
      colnames(Subject1.chip.Interaction.bySubject1.Mat) <- names(subjectCoord2List)
    }
  }
  if(mode == "both"){
    return(list(OverlapList=Subject1.chip.Overlap.bySubject1,
                OverlapMat=Subject1.chip.Overlap.bySubject1.Mat,
                IntList=Subject1.chip.Interaction.bySubject1,
                IntMat= Subject1.chip.Interaction.bySubject1.Mat))
  }else if(mode == "overlap"){
    return(list(OverlapList=Subject1.chip.Overlap.bySubject1,
                OverlapMat= Subject1.chip.Overlap.bySubject1.Mat))
  }else{
    return(list(IntList=Subject1.chip.Interaction.bySubject1,
                IntMat= Subject1.chip.Interaction.bySubject1.Mat))
  }
}
#########################################################################################################
#########################################################################################################
#example
aa = OverlapInteractionBatchExtract(subjectCoord1=Promoter1kbGR.gene,
                                    subjectCoord2List=ReMapChIP.GRanges.list,
                                    mode = "interaction" )
#########################################################################################################
#########################################################################################################

geneNameToEntrez <- function(inputNames){
  #Write a function that takes in:
  #a charachter array of gene names in these formats (UNIGENE, SYMBOL, ALIAS, ACCNUM, REFSEQ, ENSEMBL) 
  #and returns one-to-one, or "many to one" mapping of those names to ENTREZID:
  #in the format of a dataframe first column are the mappable gene names that we provided, second column is their corresponding ENTREZ ID
  
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  found <- F
  inputNames <- unique(inputNames)
  for (i in c("ALIAS","SYMBOL","ENSEMBL","REFSEQ","UNIGENE","ACCNUM")) {
    sampleCase <- sample(x = c(1:length(inputNames)), size = 10, replace = T) #choose some samples to identify the input type
    tryCatch({
      print(i)
      aa <-  select(x = org.Hs.eg.db,
                    keys = inputNames[sampleCase],
                    columns = c("ENTREZID"),
                    keytype = i)
      print(paste("CHOSEN KEY TYPE: ", i))
      ChosenKey <- i
      found <- T
      break
    }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
  }
  if(!found){
    print("Key was not any of ALIAS,SYMBOL,ENSEMBL,REFSEQ,UNIGENE,ACCNUM types. Can't convert")
    return(c("",""))
  }
  #enter First round
  FirstTry <- select(x = org.Hs.eg.db,
                     keys = inputNames,
                     columns = c("ENTREZID"),
                     keytype = ChosenKey)
  NotConverted <- FirstTry[which(is.na(FirstTry[,2])),1]
  print(paste(length(NotConverted), "names were not converted in this try"))
  FirstTry <- FirstTry[which(!is.na(FirstTry[ ,2])),]
  FirstTry <- FirstTry[(! duplicated(FirstTry[ ,1])),]
  SecondTry = c("","")
  if(length(NotConverted) > 0){
    # print("before merge")
    SecondTry <- rbind(SecondTry, geneNameToEntrez(NotConverted))
    # print("after merge")
    SecondTry <- SecondTry[-c(which(SecondTry[ ,1] == "")), ]
  }
  if (!is.null(dim(SecondTry))){
    print(paste(nrow(SecondTry), "names were converted in the next tries"))
    colnames(SecondTry) <- colnames(FirstTry)
    newRes <- rbind(FirstTry, SecondTry)
    return(newRes)
  }else{
    return(FirstTry)
  }
}
#########################################################################################################
#########################################################################################################
#example
aa = geneNameToEntrez(rownames(GSE55922RNAseqMat))

# 
# library(org.Mm.eg.db)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# aa_to_Entrez <-  select(x = org.Mm.eg.db,
#                         keys = c("NM_001145925", "NM_001145924", "NM_028137", "NM_008126", "NM_001160012"),
#                         columns = c("ENTREZID"),
#                         keytype = "REFSEQ")
# 
# aa_mouse <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
# aa_tssdf <- as.data.frame(resize(aa_mouse, width=1, fix='start'))
# aa_tssdf[aa_tssdf$gene_id %in% aa$ENTREZID,]
#########################################################################################################
#########################################################################################################
geneNameToEntrez_mouse <- function(inputNames){
  #Write a function that takes in:
  #a charachter array of gene names in these formats (UNIGENE, SYMBOL, ALIAS, ACCNUM, REFSEQ, ENSEMBL) 
  #and returns one-to-one, or "many to one" mapping of those names to ENTREZID:
  #in the format of a dataframe first column are the mappable gene names that we provided, second column is their corresponding ENTREZ ID
  
  library(org.Mm.eg.db)
  found <- F
  inputNames <- unique(inputNames)
  for (i in c("ALIAS","SYMBOL","ENSEMBL","REFSEQ","UNIGENE","ACCNUM")) {
    sampleCase <- sample(x = c(1:length(inputNames)), size = 10, replace = T) #choose some samples to identify the input type
    tryCatch({
      print(i)
      aa <-  select(x = org.Mm.eg.db,
                    keys = inputNames[sampleCase],
                    columns = c("ENTREZID"),
                    keytype = i)
      print(paste("CHOSEN KEY TYPE: ", i))
      ChosenKey <- i
      found <- T
      break
    }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
  }
  if(!found){
    print("Key was not any of ALIAS,SYMBOL,ENSEMBL,REFSEQ,UNIGENE,ACCNUM types. Can't convert")
    return(c("",""))
  }
  #enter First round
  FirstTry <- select(x = org.Mm.eg.db,
                     keys = inputNames,
                     columns = c("ENTREZID"),
                     keytype = ChosenKey)
  NotConverted <- FirstTry[which(is.na(FirstTry[,2])),1]
  print(paste(length(NotConverted), "names were not converted in this try"))
  FirstTry <- FirstTry[which(!is.na(FirstTry[ ,2])),]
  FirstTry <- FirstTry[(! duplicated(FirstTry[ ,1])),]
  SecondTry = c("","")
  if(length(NotConverted) > 0){
    # print("before merge")
    SecondTry <- rbind(SecondTry, geneNameToEntrez(NotConverted))
    # print("after merge")
    SecondTry <- SecondTry[-c(which(SecondTry[ ,1] == "")), ]
  }
  if (!is.null(dim(SecondTry))){
    print(paste(nrow(SecondTry), "names were converted in the next tries"))
    colnames(SecondTry) <- colnames(FirstTry)
    newRes <- rbind(FirstTry, SecondTry)
    return(newRes)
  }else{
    return(FirstTry)
  }
}

#########################################################################################################
#########################################################################################################
getCommonGenes <- function(expmatList){
  # Gets a list containing two or more gene expression matrices where
  # each row is one gene and each column is a condition. rownames are the names of genes in some arbitrary format.
  # Converts the names of the genes to ENTREZID and filters each expression matrix to have only the common genes
  # If two names map to the same ENTREZ ID, chooses the one with higher sum of expression values
  
  # remove row duplicates form  the matrices
  for(i in 1:length(expmatList)){
    print(paste("number of duplicate names dataset", i,
                "before processing is",
                sum(duplicated(rownames(expmatList[[i]])))))
    rowdup <- (!duplicated(rownames(expmatList[[i]])))
    currDupAll <- which(duplicated(rownames(expmatList[[i]])))
    if (length(currDupAll) > 0){
      for(j in 1:length(currDupAll)){# go through the duplicates and for each set of duplicates only keep the one with highest sum of expression across all conditions
        currDupSetInd <- which(rownames(expmatList[[i]]) == rownames(expmatList[[i]])[currDupAll[j]])
        currDupSet <- rownames(expmatList[[i]])[currDupSetInd]
        curRowSumsMaxInd <- which.max(rowSums(expmatList[[i]][match(currDupSet, rownames(expmatList[[i]])), ])) 
        rowdup[currDupSetInd[1]] <- 0
        rowdup[currDupSetInd[curRowSumsMaxInd]] <- 1
      }
    }
    expmatList[[i]] <- expmatList[[i]][which(rowdup == 1), ]
    print(paste("number of duplicate names dataset",i,
                "after processing is",
                sum(duplicated(rownames(expmatList[[i]])))))
  }
  ######
  NameList <- list()
  for(i in 1:length(expmatList)){
    NameList[[i]] <- rownames(expmatList[[i]])
  }
  # Convert each entry to entrez id
  ConvertedNameList <- lapply(NameList, geneNameToEntrez)
  if(length(ConvertedNameList) != length(NameList) ){
    print("problem in conversion to ENTREZ (start)")
  }
  for(i in 1:length(ConvertedNameList)){
    ConvertedNameList[[i]] <- cbind(ConvertedNameList[[i]],
                                    matrix(1L, nrow = nrow(ConvertedNameList[[i]]), ncol = 1))
    currDupAll <- which(duplicated(ConvertedNameList[[i]][,2]))
    print(paste("number of duplicated elements of dataset", i, "is", length(currDupAll)))
    if (length(currDupAll) > 0){
      for(j in 1:length(currDupAll)){# go through the duplicates and for each set of duplicates only keep the one with highest sum of expression across all conditions
        currDupSetInd <- which(ConvertedNameList[[i]][ , 2] == ConvertedNameList[[i]][currDupAll[j], 2])
        ConvertedNameList[[i]][currDupSetInd, 3] <- 0
        currDupSet <- ConvertedNameList[[i]][currDupSetInd, 1]
        curRowSumsMaxInd <- which.max(rowSums(expmatList[[i]][match(currDupSet,
                                                                    rownames(expmatList[[i]])),])) 
        ConvertedNameList[[i]][currDupSetInd[curRowSumsMaxInd], 3] <- 1
        }
    }
    print(paste("number of ommited entries of dataset number", i,
                "is", sum(ConvertedNameList[[i]][,3] == 0)))
    ConvertedNameList[[i]] <- ConvertedNameList[[i]][which(ConvertedNameList[[i]][ , 3] == 1), ]
    print(paste("nrow of new ConvertedNameList is",
                nrow(ConvertedNameList[[i]]),
                "ncol of new ConvertedNameList is",
                ncol(ConvertedNameList[[i]] )))
    if(sum(duplicated(ConvertedNameList[[i]][ , 2])) > 0){
      print("Something is wrong, we still have duplicates")
    }
    # remove the rows of expression dataset which are not converted to ENTREZ
    expmatList[[i]] <- expmatList[[i]][which(rownames(expmatList[[i]]) %in% ConvertedNameList[[i]][ , 1]), ]
    # rename the rows to entrez ID
    rownames(expmatList[[i]]) <- ConvertedNameList[[i]][match(rownames(expmatList[[i]]),
                                                              ConvertedNameList[[i]][ , 1]), 2]
  }# end of conversion of each set name to ENTREZ and removing duplicates
  # get the intersect of all converted entrez of all datasets
  ConvertedNameEntrezList <- list()
  for(i in 1:length(ConvertedNameList)){
    ConvertedNameEntrezList[[i]] <- ConvertedNameList[[i]][ , 2]
  }
  print(paste("number of unique converted to Entrez in each dataset is",
              unlist(lapply(ConvertedNameEntrezList,length))))
  MutualEntrez <- Reduce(intersect, ConvertedNameEntrezList)
  print(paste("length of the mutual set is", length(MutualEntrez)))
  MutualDataSetList <- expmatList
  for(i in 1:length(MutualDataSetList)){
    MutualDataSetList[[i]] <- MutualDataSetList[[i]][which(rownames(MutualDataSetList[[i]]) %in% MutualEntrez), ]
    MutualDataSetList[[i]] <- MutualDataSetList[[i]][match(rownames(MutualDataSetList[[i]]), MutualEntrez ), ]
    sortedInd <- sort(rownames(MutualDataSetList[[i]]), index.return = T)$ix
    MutualDataSetList[[i]] <- MutualDataSetList[[i]][sortedInd, ]
  }
  return(MutualDataSetList)
}
#########################################################################################################
#########################################################################################################
#example
aa <- getCommonGenes(list(GSE108308.RNAseq,
                          GSE55922.RNAseq,
                          GSE67295.RNAseq,
                          GSE68358.RNAseq,
                          GSE76507.RNAseq,
                          GSE86316.RNAseq,
                          GSE73663.Inh.RNAseq
                          ))
#########################################################################################################
#########################################################################################################
clustrator <- function(ExpressionMat,
                       nfolds=6,
                       .iter.max=400,
                       plot_mfrow=1,
                       ylim_custum=1,
                       .col_vector=col_vector,
                       Exp.metadata=0,
                       exportPlot=F,
                       feedCluster=numeric(0),
                       filename="GeneCluster.pdf",
                       no.control = F){
  # function to cluster a expression matrix, using kmeans clustering, and plot each clustered group separately
  # inputs: 1.ExpressionMat: expression matrix where each row is a gene, each column is a condition
  # 2. nfolds: number of clusters
  # 3. iter.max: max iteration number of kmeans algorithm
  # 4. plot_mfrow: number of rows and columns of the plot eg: c(4,4)
  # 5. ylim_custum: is the y limits of the plot eg: c(-10,10)
  # 6. col_vector: vector of colors eg: c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F") )
  # 7. Expmetadata: is a dataframe with two columns first: DatasetNo , second: expNo , which store the number of dataset and number of eperiment for each column of gene expression data.eg DataSet.Experiment.Metadata
  # 8. Export : if T plots in a file
  # 9. feedCluster : defualt is numeric(0), causing the function to calculate clusters itself, but can get output of a kmeans clustering as input
  # 10. filename: name of the file for the plots
  # 11. no.control : if the expression data does not include controls set this to True
  # returns a list of two components: first the index of genes in each cluster, second the kmeans cluster object
  if (length(feedCluster) == 0){
    Clustered.DataSet <- kmeans(ExpressionMat, nfolds, iter.max=.iter.max)
    # if the default kmeans fails for some reason: use "MacQueen" algorithm instead
    if (Clustered.DataSet$ifault == 4){
      Clustered.DataSet <- kmeans(ExpressionMat, Clustered.DataSet$centers, algorithm="MacQueen")
    }
  }else{
    Clustered.DataSet <- feedCluster
  }
  
  prev.Par <- par()
  if (length(plot_mfrow) == 1){# if plot_mfrow is not set, set it this way:
    if(nfolds < 7){
      plot_mfrow <- c (nfolds, 1)
    }else if(nfolds < 13){
      plot_mfrow <- c(ceiling(nfolds / 2), 2)
    }else if(nfolds < 19){
      plot_mfrow <- c(ceiling(nfolds / 3), 3)
    }else{
      plot_mfrow <- c(ceiling(sqrt(nfolds)), ceiling(sqrt(nfolds)))
    }
  }
  # plot the genes in each cluster separately
  Clustered.DataSet.index.perCluster <- list()
  for (eaf in 1:nfolds){ # separate the index of genes in each cluster
    Clustered.DataSet.index.perCluster[[eaf]] <- which(Clustered.DataSet$cluster %in% eaf)
  }
  if (exportPlot){ # write to file if exportPlot is true
    pdf(file=filename, width=11, height=8.5)
  }
  par(mfrow=plot_mfrow, mar=c(1,4,1,1))
  if (length(ylim_custum) == 1){ # if ylim_custum is not set, set it to range of data
    ylim_custum <- range(ExpressionMat)
  }
  for (eaf in 1:nfolds){ # iterate through clusters
    # plot the first gene in the cluster
    plot(ExpressionMat[Clustered.DataSet.index.perCluster[[eaf]][1],],
         ylim=ylim_custum, type="l", xlab="", ylab="", xaxt="n")
    # drawing vertical and horizontal lines for grid and metadata
    abline(h=0, lwd=2, col=2, lty=2)
    abline(h=c(-10, -5, 5, 10), lwd=0.2 , col=1, lty=3)
    if (length(Exp.metadata) > 1){
      if(no.control){#if the expression and metadata file do not include controls
        # if we have metadata information create vertical lines separating experiments and datasets
        abline(v=which(Exp.metadata$expNo == 2) - 0.2, col=4, lty=2, lwd=0.5)
      }else{
        # if we have metadata information create vertical lines separating experiments and datasets
        abline(v=which(Exp.metadata$expNo == 1) - 0.2, col=4, lty=2, lwd=0.5)
      }

      a <- numeric(0)
      for (i in 1:length(unique(Exp.metadata$DatasetNo))){
        a <- c(a, which(Exp.metadata$DatasetNo == i)[1])
      }
      abline(v=a - 0.2, col=3, lty=5, lwd=2)
    }
    text((ncol(ExpressionMat) - 5), ylim_custum[1] + (ylim_custum[2] - ylim_custum[1])/20, Clustered.DataSet$size[eaf]) # write the size of the cluster
    for (j in 2:length(Clustered.DataSet.index.perCluster[[eaf]])){ # draw the rest of the genes
      lines(ExpressionMat[Clustered.DataSet.index.perCluster[[eaf]][j],],
            col=.col_vector[(j%%(length(.col_vector)))])
    }
    if(length(Clustered.DataSet.index.perCluster[[eaf]]) > 1){
      #draw the mean of the cluster
      lines(colMeans(ExpressionMat[Clustered.DataSet.index.perCluster[[eaf]],]),
            col=.col_vector[18], lwd=4)
    }
  }
  if(exportPlot){
    dev.off()
  }
  par(mfrow = prev.Par$mfrow, mar = prev.Par$mar)
  return(list(ClusterIndex = Clustered.DataSet.index.perCluster,
              ClusterObject = Clustered.DataSet))
}

#########################################################################################################
#########################################################################################################
#example
aa <- clustrator(DataSet.Experiment.matrix.Batch.vgt1.LFC.ER,
                 nfolds=3,
                 plot_mfrow=c(3, 1),
                 Exp.metadata=DataSet.Experiment.Metadata,
                 exportPlot=F)
#########################################################################################################
#########################################################################################################
#wrapper for clustrator function to feed unique clusters to it
clustratorWrapper <- function(.ExpressionMat,
                              .nfolds=6,
                              ..iter.max=400,
                              .plot_mfrow=1,
                              .ylim_custum=1,
                              ..col_vector=col_vector,
                              .Exp.metadata=DataSet.Experiment.Metadata,
                              .exportPlot=F,
                              #feedCluster=numeric(0),
                              .filename="GeneCluster",
                              RepeatMax=100,
                              harsh=T,
                              harshThreshTimes = 2,
                              dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Expression",
                              .no.control = F
                              ){
  # function to cluster a expression matrix, using kmeans clustering, several times, each time a unique cluster, and plot each clustered group separately
  # inputs: 1.ExpressionMat: expression matrix where each row is a gene, each column is a condition
  # 2. nfolds: number of clusters
  # 3. iter.max: max iteration number of kmeans algorithm
  # 4. plot_mfrow: number of rows and columns of the plot eg: c(4,4)
  # 5. ylim_custum: is the y limits of the plot eg: c(-10,10)
  # 6. col_vector: vector of colors eg: c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F") )
  # 7. Expmetadata: is a dataframe with two columns first: DatasetNo , second: expNo , which store the number of dataset and number of eperiment for each column of gene expression data.eg DataSet.Experiment.Metadata
  # 8. Export : if T plots in a file
  # 9. feedCluster : defualt is numeric(0), causing the function to calculate clusters itself, but can get output of a kmeans clustering as input
  # 10. filename: name of the file for the plots
  # 11. RepeatMax is the number of attempts to create a new unique cluster
  # 12. harsh : if True the new clusters won't be accepted if their size is similar to the previous ones (sum of absolute differences has to be greater than (number of folds * harshThreshTimes)),
  # 12. harsh : if False the new clusters won't be accepted if their size is exactly the same as one of the previous ones
  # 13. harshThreshTimes (sum of absolute differences has to be greater than (number of folds * harshThreshTimes)) if harsh is T
  # 14. dir : directory to store the results if Export is True
  # 15. .no.control : if the expression data does not include controls set this to True

  # returns a list of clustrator outputs which are: lists of two components: first the index of genes in each cluster, second the kmeans cluster object

  # matrix to hold the sorted size of clusters in order to compare with the new ones later
  prev_dir <- getwd()
  if(!dir.exists(paths = dir)){#create the directory if it doesn't exist
    dir.create(path=dir, showWarnings = T,recursive = T)
  }
  setwd(dir)
  ClusterSizes <- matrix(0L,nrow=RepeatMax, ncol=.nfolds)
  Clustered.DataSet.set <- list()
  for(attemp in 1:RepeatMax){
    Clustered.DataSet.set[[attemp]] <- kmeans(.ExpressionMat, .nfolds, iter.max=..iter.max)
    # if the default kmeans fails for some reason: use "MacQueen" algorithm instead
    if (Clustered.DataSet.set[[attemp]]$ifault == 4){
      Clustered.DataSet.set[[attemp]] <- kmeans(.ExpressionMat, Clustered.DataSet.set[[attemp]]$centers,
                                                algorithm="MacQueen")
    }
    ClusterSizes[attemp, ] <- sort(Clustered.DataSet.set[[attemp]]$size)
  }
  if (harsh){# keeping only distinct clusters from the Clustered.DataSet.set
    acceptedIndex <- matrix(1L, nrow=length(Clustered.DataSet.set), ncol=1 ) # stores 1 for unique and 0 for repeated clusters
    for(att in 2:length(Clustered.DataSet.set)){
      for(prev in 1:(att - 1)){# check if the sum of absolute differences between the size of this cluster and any previous cluster is greater than harshThreshTimes*(nfolds)
        sizeDif <- sum(abs(ClusterSizes[prev,] - ClusterSizes[att,]))
        if(sizeDif < (.nfolds * harshThreshTimes)){
          acceptedIndex[att,1] <- 0
          break
        }
      }
    }
    #
    ClusterSizes.Uniq <- ClusterSizes[which(acceptedIndex == 1), ]
    # filer for just the accepted ones
    Clustered.DataSet.set.unique <- Clustered.DataSet.set[which(acceptedIndex == 1)]
    print(sum(acceptedIndex == 1))
    onlyOne = F
    if(sum(acceptedIndex == 1) == 1){
      onlyOne = T
    }
  }else{# if harsh is false just check if the sorted size vector is duplicated or not
    Clustered.DataSet.set.unique <- Clustered.DataSet.set[!duplicated(ClusterSizes)]
    ClusterSizes.Uniq <- ClusterSizes[which(!duplicated(ClusterSizes)), ]
    onlyOne = F
    if(sum(!duplicated(ClusterSizes)) == 1){
      onlyOne = T
    }
  }
  
  par(mforw = c(1, 1), mar = c(4, 4, 4, 4))
  if (!onlyOne){
    boxplot.matrix(ClusterSizes.Uniq, use.cols = T, main = "ClusterSizes")
  }else{
    plot(ClusterSizes.Uniq, type = "l")
  }
  
  # now feed the unique clusters to clustrator
  resultsList <- list()
  for(ClusNo in 1:length(Clustered.DataSet.set.unique)){
    resultsList[[ClusNo]] <- clustrator(ExpressionMat = .ExpressionMat,
                                        nfolds=.nfolds,
                                        .iter.max=..iter.max,
                                        plot_mfrow=.plot_mfrow,
                                        ylim_custum=.ylim_custum,
                                        .col_vector=..col_vector,
                                        Exp.metadata=.Exp.metadata,
                                        exportPlot=.exportPlot,
                                        feedCluster=Clustered.DataSet.set.unique[[ClusNo]],
                                        filename=paste(.filename,"_",as.character(.nfolds),
                                                       "_",as.character(ClusNo),".pdf",sep = ""),
                                        no.control = .no.control)
  }
  setwd(prev_dir)
  return(resultsList)
}

#########################################################################################################
#########################################################################################################
#example
aa <- clustratorWrapper(.ExpressionMat=DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.gt1.ER,
                        .nfolds=4,
                        .exportPlot=T,
                        dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Expression/Clustering_quantile_2334",
                        RepeatMax=50)

#########################################################################################################
#########################################################################################################
WriteFastaOfBag <- function(enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                            Enhancer.Sequence=MCFEnhancers[,"sequ"],
                            EnhancerGR=MCFEnhancersGR,
                            promoter.promoter.int=Promoter.gene.PromoterIntList,
                            promoter.sequence=Promoter.Sequence.Char,
                            PromoterGR=Promoter1kbGR.gene,
                            output.File.Name = "RegulatoryElements",
                            max.enh.length = numeric(0),
                            return.Length = F,
                            returnListNotwrite=F,
                            write_from_list=list()){
  # Function to write the bag of regulatory elements of each gene
  # input:
  # 1) enhancer.per.gene : a list containing list of enhancers for each gene eg: Vicinity100kb.Enhancer.plusPromoter.By.gene: subsetted for the genes of intrest
  # 2) Enhancer.Sequence : charachter vector of all enhancers: MCFEnhancers[,"sequ"]
  # 3) promoter.promoter.int : a list containing list of promoters interacting with each gene eg: Promoter.gene.PromoterIntList: subsetted for the genes of intrest
  # 4) promoter.sequence : sequence of the promoters of genes in the same order as input1 and input 3--->this should be named by genes
  # 5) output.File.Name : name of the output file
  # 6) max.enh.length : max.enh.length
  # 7) return.Length : returns the length of enhancers written
  # 8) returnListNotwrite : returns a list of sequences: one entry per gene: contains charachter vector of reg elements of that 
  #    gene. if True it also returns a list containing a GRAnges object per gene which contains the coordinates of the associated promoters or enhancers of each gene
  # 9) write_from_list : is a named list. instead of looking up the overlapping and interacting enhancers and promoters for each gene, just write the ones specified in the list
  #    list has one entry per gene. the entry contains a charachter vector, with one entry per regulatory element associated with the gene, containing the sequence of that regulatory element.
  # EnhancerGR : GRange object of all enhancers, order is the same as Enhancer.Sequence
  # PromoterGR : GRange object of all promoters, order is the same as promoter.sequence
  # return.GRanges : 
  # IT IS IMPORTANT THAT THE FIRST AND THIRD AND Forth ENTRIES ARE NAMED, (with corresponding genes names)
  # IT IS IMPORTANT THAT THE FIRST AND Third  ENTRIES SHOULD HAVE THE SAME LENGTH, WHICH IS THE NUMBER OF GENES OF INTEREST
  # return.Length : if T it returns a list with length equal to number of genes. each entry is the length of regulatory elments associated with that gene.
  if(length(write_from_list) > 0){
    for(l_cur_gene in 1:length(write_from_list)){
      for(reg.el in 1:length(write_from_list[[l_cur_gene]])){
        if(l_cur_gene == 1 & reg.el == 1){
          write.table(c(paste(">", paste(names(write_from_list)[l_cur_gene],
                                         as.character(reg.el), sep="_"), sep=""),
                        write_from_list[[l_cur_gene]][reg.el]), file=paste0(output.File.Name, ".fa"),
                      sep="\n", row.names=F, col.names=F, quote=F, append=F)
        }else{
          write.table(c(paste(">", paste(names(write_from_list)[l_cur_gene],
                                         as.character(reg.el), sep="_"), sep=""),
                        write_from_list[[l_cur_gene]][reg.el]), file=paste0(output.File.Name, ".fa"),
                      sep="\n", row.names=F, col.names=F, quote=F, append=T)
        }
      }
    }
  }else{
    if (length(max.enh.length) == 0){
      max.enh.length <- max(unlist(lapply(Enhancer.Sequence, nchar))) + 1
    }
    
    Gene.Enhancer.Length <- list()
    SequenceList <- list()
    GRangeList <- list()
    for(gene in 1:length(enhancer.per.gene)){
      # set the size of the gene's regulatory element bag as 1+interacting promoters+enhancers 
      gene.name <- names(enhancer.per.gene)[gene]
      gene.index.orig <- match(gene.name, names(promoter.sequence))
      cur.gene.bag <- character(length = (1+length(promoter.promoter.int[[gene]])
                                          +length(enhancer.per.gene[[gene]])))
      # first element of each bag is its promoter
      cur.gene.bag[1] <- toupper(promoter.sequence[gene.index.orig])
      cur.gene.GRange <- PromoterGR[gene.index.orig, ]
      cnt <- 2
      # add interacting promotes if there are any
      if (length(promoter.promoter.int[[gene]]) > 0){
        for (int.prom in promoter.promoter.int[[gene]]){
          cur.gene.bag[cnt] <- toupper(promoter.sequence[int.prom])
          cur.gene.GRange <- c(cur.gene.GRange, PromoterGR[int.prom, ])
          cnt <- cnt + 1
        }
      }
      # add enhancers if there are any
      if (length(enhancer.per.gene[[gene]]) > 0){
        for (enh in enhancer.per.gene[[gene]]){
          if (nchar(Enhancer.Sequence[enh]) < max.enh.length){
            cur.gene.bag[cnt] <- toupper(Enhancer.Sequence[enh])
            cur.gene.GRange <- c(cur.gene.GRange, EnhancerGR[enh, ])
            cnt <- cnt + 1
          }
        }
      }
      cur.gene.bag <- cur.gene.bag[cur.gene.bag != ""]
      # write the current set to a file, appending to the previous one
      if (names(promoter.sequence)[gene.index.orig] != gene.name){
        print("something is wrong")
      }
      Gene.Enhancer.Length[[gene]] <- nchar(cur.gene.bag)
      SequenceList[[gene]] <- cur.gene.bag
      GRangeList[[gene]] <- cur.gene.GRange
      if(!returnListNotwrite){
        for(reg.el in 1:length(cur.gene.bag)){
          write.table(c(paste(">", paste(names(promoter.sequence)[gene.index.orig],
                                         as.character(reg.el), sep="_"), sep=""),
                        cur.gene.bag[reg.el]), file=paste0(output.File.Name, ".fa"),
                      sep="\n", row.names=F, col.names=F, quote=F, append=T)
        }
      }
    }#loop over genes
    names(Gene.Enhancer.Length) <- names(enhancer.per.gene)
    names(SequenceList) <- names(enhancer.per.gene)
    names(GRangeList) <- names(enhancer.per.gene)
    
    if(returnListNotwrite & return.Length){
      result <- list(sequence = SequenceList, length = Gene.Enhancer.Length, GRange = GRangeList)
      return(result)
    }else if(returnListNotwrite){
      return(SequenceList)
    }else{
      return(Gene.Enhancer.Length)
    }
    
  }
}

#########################################################################################################
#########################################################################################################
#example
aa  <- which(rownames(DataSet.Experiment.matrix.Batch.vgt1) %in% Genes.Associated.REMAP.ER.Entrez)
aaa <- match(rownames(DataSet.Experiment.matrix.Batch.vgt1)[aa], names(Vicinity100kb.Enhancer.plusPromoter.By.gene))
WriteFastaOfBag(enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene[aaa],
                            Enhancer.Sequence=MCFEnhancers[,"sequ"],
                            promoter.promoter.int=Promoter.gene.PromoterIntList[aaa],
                            promoter.sequence=Promoter.Sequence.Char)

WriteFastaOfBag(write_from_list = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_SimAnnFiltered$Chopped_Seq, output.File.Name = "testListSeq")
WriteFastaOfBag(write_from_list = aa_enhancer_list, output.File.Name = "testListSeq")

#########################################################################################################
#########################################################################################################
AddpsudoRow <- function(pwm, pseudo=0.001){
  # write a function to get a motif and output the motif with a
  # pseudocount added so it doesn't have zero entries
  
  # adds a pseudo count to the motif
  # pwm (columns A,C,G,T, rows bases)
  for (i in 1:nrow(pwm)){
    pwm[i, ] <- (pwm[i, ] + pseudo) / (1 + 4 * pseudo)
  }
  return(pwm)
}
#########################################################################################################
#########################################################################################################
RemoveDegenRow <- function(motif,
                           background=c(0.252, 0.249, 0.248, 0.251),
                           threshold=0.3,
                           .pseudo=0.001){
  # write a function to get a motif and output the motif without degenerate rows at beginning and end
  # motif : is the pwm (columns A,C,G,T, rows bases)
  # threshold : is the minimum acceptable information content of a column
  IC <- numeric(length = nrow(motif))
  psudoMotif <- AddpsudoRow(pwm=motif, pseudo=.pseudo)
  for (i in 1:nrow(psudoMotif)){
    for (j in 1:4){
      IC[i] <- IC[i] + psudoMotif[i,j] * log2((psudoMotif[i,j]) / (background[j]))
    }
  }# information content of each column is calculated
  Infmt <- which(IC >= threshold)
  NewMotif <- psudoMotif[Infmt[1]:Infmt[length(Infmt)], ]
  rownames(NewMotif) <- c(1:nrow(NewMotif))
  return(NewMotif)
}
#########################################################################################################
#########################################################################################################
#example
aa <- TF.motifs[[3]]
seqLogo(t(aa))
aaa <- RemoveDegenRow(aa)
seqLogo(t(aaa))
#########################################################################################################
#########################################################################################################
PWMtoCount <- function(PWM,
                       remove.degen=T,
                       times=10000,
                       pseudo=0.0001,
                       .threshold=0.3){
  # input: PWM is matrix where columns are A,C,G,T and rows are bases
  # times is the number to multiply each entry to
  # pseudo is the count being added to PWM before converting to countmatrix
  # .threshold is the threshold on information content of RemoveDegenRow
  # remove.degen if TRUE runs RemoveDegenRow function ro remove the degenarate bases, it also adds a pseudi count
  # outputs: integer count for each entry in the same shape as PWM
  if (remove.degen){
    count.mat <- round(times * RemoveDegenRow(PWM, .pseudo=pseudo, threshold=.threshold))
  }else{
    count.mat <- round(times * PWM)
  }
    return(count.mat)
}
#########################################################################################################
#########################################################################################################
#example
aa <- TF.motifs[[2]]
seqLogo(t(aa))
aaa <- PWMtoCount(aa)
#########################################################################################################
#########################################################################################################
MotifWriter <- function(motif.List, pseudo = 0.01, output.File.Name){
  # inputs a list of motifs and outputs a .wtmx file containg all motifs in the format readable by GEMSTAT
  # the motiflist has to be named with the name of motifs
  # someNumber is the number in wtmx file: it is the psudo count being used in Gemstat
  for (motif.nu in 1:length(motif.List)){
    write.table(paste(">", paste(names(motif.List)[motif.nu],
                                 as.character(nrow(motif.List[[motif.nu]])),
                                 as.character(pseudo) ,sep="\t"), sep = ""),
                file=paste0(output.File.Name, ".wtmx"),
                sep="\n", row.names=F, col.names=F, quote=F, append=!(motif.nu == 1))
    write.table(motif.List[[motif.nu]], file=paste0(output.File.Name, ".wtmx"),
                row.names=FALSE, col.names=FALSE, quote=F, append=T, sep = "\t")
    cat(c("<", "\n"), file = paste0(output.File.Name, ".wtmx"), append=T, sep = "")
    # write.table(paste("<","",sep = "\n"), file=paste0(output.File.Name, ".wtmx"),
    #             row.names=FALSE, col.names=FALSE, quote=F, append=T, sep = "" )
  }
}
#########################################################################################################
#########################################################################################################
#example
MotifWriter(motif.List = TF.motifs.count, output.File.Name = "motifs")
#########################################################################################################
#########################################################################################################
ExpressionWriter <- function(expressionMat,
                             output.File.Name = "ExpressionMat",
                             multi_enhancer_map=integer(0),
                             from_log_fold_change=F){
  # This function gets expression matrix: rows:genes , columns:conditions
  # The rows should be named for genes
  # outputs .tab file containing expressionmat in GEMSTAT readable format
  # multi_enhancer_map is an integer vector containing the mapping of enhancers to genes. It has 
  #  one entry per enhancer which corresponds to its associated gene
  # from_log_fold_change : is a logical. if True the function writes an expression matrix from a
  #  log fold change matrix. it will create two conditions for each fold change condition, where
  #  the first one is the control and the second one is the treatment.
  # normalizes each gene: divides by the max value to have it between 0 and 1
  if(from_log_fold_change){
    # expressionMat[is.na(expressionMat)] <- log2(1234)
    raw_exp_mat <- matrix(nrow = nrow(expressionMat), ncol = 2*ncol(expressionMat))
    rownames(raw_exp_mat) <- rownames(expressionMat)
    for(i in 1:ncol(expressionMat)){
      raw_exp_mat[, (2*i - 1)] <- rep(1, nrow(expressionMat))
      raw_exp_mat[, (2*i)]     <- mapply(expressionMat[, i], FUN = function(x){2^x})
    }
    raw_exp_mat[is.na(raw_exp_mat)] <- 1234
    aa_mx_1234 <- function(x){
      return(max(x[x != 1234]))
    }
    max_row <- apply(raw_exp_mat, 1, aa_mx_1234)
    expressionMat_norm <- raw_exp_mat / max_row
    expressionMat <- raw_exp_mat
  }else{
    max_row <- apply(expressionMat, 1, max, na.rm = T)
    max_row[max_row == 0] <- 1e-6
    expressionMat_norm <- expressionMat / max_row
    expressionMat_norm[is.na(expressionMat_norm)] <- 1234
  }

  if(length(multi_enhancer_map) > 0){
    print("be carefull about indexing, the multi enhancer part uses 0 based index for genes")
    multi_expression_mat <- matrix(nrow = length(multi_enhancer_map), ncol = ncol(expressionMat))
    rownames(multi_expression_mat) <- character(length(multi_enhancer_map))
    for(i in 0:(nrow(expressionMat)-1)){
      cur_ind <- which(multi_enhancer_map %in% i)
      multi_expression_mat[cur_ind, ] <- matrix(rep(expressionMat_norm[i+1,],
                                                 length(cur_ind)),
                                             byrow = T,
                                             nrow = length(cur_ind))
      rownames(multi_expression_mat)[cur_ind] <- paste(rownames(expressionMat)[i+1], c(1:length(cur_ind)), sep = "_")
    }
    expressionMat_norm <- multi_expression_mat
  }
  cat(c("Rows",
        rep(as.character(c(1:(ncol(expressionMat) - 1))), as.integer(ncol(expressionMat) > 1)),
        paste0(as.character(ncol(expressionMat)), "\n")),
      file=paste0(output.File.Name, ".tab"), sep = "\t")
  write.table(expressionMat_norm, file=paste0(output.File.Name, ".tab"),
              row.names=T, col.names=FALSE, quote=F, append=T, sep = "\t")
}
#########################################################################################################
#########################################################################################################
#example
ExpressionWriter(expressionMat = TF.Expression)
ExpressionWriter(expressionMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, output.File.Name = "nolog")
ExpressionWriter(expressionMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, from_log_fold_change = T,output.File.Name = "log")
aa_map <- integer(0)
for(i in 1:nrow(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)){
  aa_map <- c(aa_map, rep(i, sample(c(1:5), 1)))
}
ExpressionWriter(expressionMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, from_log_fold_change = T,output.File.Name = "log_enh", multi_enhancer_map = aa_map)
#########################################################################################################
#########################################################################################################
LogFoldTransform <- function(expressionMatrix,
                             metadataM=DataSet.Experiment.Metadata){
  # write a function to compute log fold change of each condition given the expression mat and metadata
  
  # inputs: 
  # 1.ExpressionMat: expression matrix where each row is a gene, each column is a condition
  # 2. Expmetadata: is a dataframe with two columns first: DatasetNo , second: expNo , which
  # store the number of dataset and number of eperiment for
  # each column of gene expression data.eg DataSet.Experiment.Metadata
  
  myExpMat <- expressionMatrix
  control.index <- which(metadataM$expNo == 1)
  control.index <- c(control.index, (ncol(myExpMat) + 1))
  # set negative values to zero
  myExpMat[myExpMat < 0] <- 0
  # add 1 to all
  myExpMat <- myExpMat + 1
  # divide by control
  for (cnt in 1:(length(control.index) - 1)){
    myExpMat[ ,control.index[cnt]:(control.index[cnt + 1] - 1)] <- myExpMat[ ,control.index[cnt]:(control.index[cnt + 1] - 1)] / myExpMat[ ,control.index[cnt]]
  }
  # compute log2
  myExpMat <- log2(myExpMat)
  return(myExpMat)
}
#########################################################################################################
#########################################################################################################
#example
aaa <- LogFoldTransform(DataSet.Experiment.matrix.Batch.vgt1)
#########################################################################################################
#########################################################################################################
# Function to write hal Jobs for GEMSTAT runs
Hal_job_writer <- function(exp.nu, cluster.nu,
                           seqName="RegulatoryElementsSequence",
                           expressionName="ExpressionMat",
                           Shared_dir="/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/",
                           Ensemble = F,
                           home_dir = character(0),
                           GEMSTAT_dir = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_MultiEnh5/seq2expr",
                           par.exp.nu = exp.nu,
                           .a=character(0),
                           .o="DIRECT",
                           .c=character(0),
                           .i=character(0),
                           .r=character(0),
                           .oo=character(0),
                           .mc=character(0),
                           .p=character(0),
                           .rt=character(0),
                           .na=character(0),
                           .ct=character(0),
                           .oq=character(0),
                           .sigma=character(0),
                           .softmin_groups=character(0),
                           .control_treat_map=character(0)
                           ){
  # Function to write hal Jobs for GEMSTAT runs
  # exp.nu is the number if the experiment
  # cluster.nu is the number of clusters in this experiment
  # seqName is the prefix for sequence file names
  # expressionName is the prefix for expression file names
  # This function assumes a certain hirarchy of directories for Inputs:
  # This is shared for all: Shared_dir/
  # sequence: Shared_dir/Experiment_exp.nu/Cluster_j/Input/Sequence/seqName_i.fa
  # gene expression: Shared_dir/Experiment_exp.nu/Cluster_j/Input/Gene_expression/expressionName_i.tab
  # motifs: Experiment_exp.nu/Cluster_j/Input/Motifs/motifs.wtmx
  # TF expression: Experiment_exp.nu/Cluster_j/Input/TF_Expression/TF_Expression_quantile.tab
  # Shared_dir is the root directory where the data and GEMSTAT program are stored
  # GEMSTAT program to call: Shared_dir/src_MultiEnh5/seq2expr
  # Ensemble : if you want to write jobs for an ensemble of parameters T otherwise false
  # GEMSTAT_dir : address of the gemstat version to be called
  # par.exp.nu : the experiment number that you want to get the parameter files from: by default set to current experiment
  # -a annFile
  # -o modelOption
  # -c coopFile 
  # -i factorInfoFile
  # -r repressionFile 
  # -oo objOption 
  # -mc maxContact 
  # -p parFile 
  # -rt repressionDistThr 
  # -na nAlternations 
  # -ct coopDistThr 
  # -oq if any argument provided for this, there will be a separate qBTM assigned to each CRM otherwise they will share a qBTM
  # -sigma factorIntSigma
  # -softmin_groups address of enhancer to gene mapping file
  # -control_treat_map address of treatment/control meta file
  
  GEMSTAT_call <- GEMSTAT_dir
  for (cl in 1:cluster.nu){
    sequence_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                          "/Cluster_", cl, "/Inputs/Sequence/",
                          seqName, "_", cl, ".fa", sep="")
    gene_expression_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                                 "/Cluster_", cl, "/Inputs/Gene_expression/",
                                 expressionName, "_", cl, ".tab", sep="")
    motif_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                       "/Cluster_", cl,
                       "/Inputs/Motifs/motifs.wtmx", sep="")
    TF_expression_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                               "/Cluster_", cl, 
                               "/Inputs/TF_Expression/TF_Expression_quantile.tab", sep="")
    output_file_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                             "/Cluster_", cl, 
                             "/Outputs/", "Experiment_", exp.nu, "_Cluster_", cl, ".out", sep="")
    log_file_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                          "/Cluster_", cl, 
                          "/Outputs/", "Experiment_", exp.nu, "_Cluster_", cl, ".log", sep="")
    if(Ensemble){
      par_folder_dir <- paste(Shared_dir, "Experiment_", par.exp.nu,
                              "/Cluster_", cl, 
                              "/Inputs/Parameters/", sep="")
      par_folder_dir_home <- paste(home_dir, "Experiment_", par.exp.nu,
                                   "/Cluster_", cl, 
                                   "/Inputs/Parameters/", sep="")
      par_files <- list.files(par_folder_dir_home)
      for (parf in 1:length(par_files)){
        .p <- paste0(par_folder_dir, par_files[parf])
        output_file_dir_par <- paste0(output_file_dir, par_files[parf], ".outEns" )
        log_file_dir_par <- paste0(log_file_dir, par_files[parf], ".logEns" )
        
        cat(c(GEMSTAT_call, "-s", sequence_dir,
              "-e", gene_expression_dir,
              "-m", motif_dir,
              "-f", TF_expression_dir,
              "-fo", output_file_dir_par,
              rep("-a", length(.a)), .a,
              rep("-o", length(.o)), .o,
              rep("-c", (length(.c) > 0)), .c,
              rep("-i", length(.i)), .i,
              rep("-r", length(.r)), .r,
              rep("-oo", length(.oo)), .oo,
              rep("-mc", length(.mc)), .mc,
              rep("-p", length(.p)), .p,
              rep("-rt", length(.rt)), .rt,
              rep("-na", length(.na)), .na,
              rep("-ct", length(.ct)), .ct,
              rep("-oq", length(.oq)), .oq,
              rep("-sigma", length(.sigma)), .sigma,
              rep("-softmin_groups", length(.softmin_groups)), .softmin_groups,
              rep("-control_treat_map", length(.control_treat_map)), .control_treat_map,
              " >> ",log_file_dir_par,  "\n"),
            file=paste0("Experiment_",exp.nu, "job"), sep=" ", append=T)
      }
    }else{
      cat(c(GEMSTAT_call, "-s", sequence_dir,
            "-e", gene_expression_dir,
            "-m", motif_dir,
            "-f", TF_expression_dir,
            "-fo", output_file_dir,
            rep("-a", length(.a)), .a,
            rep("-o", length(.o)), .o,
            rep("-c", (length(.c) > 0)), .c,
            rep("-i", length(.i)), .i,
            rep("-r", length(.r)), .r,
            rep("-oo", length(.oo)), .oo,
            rep("-mc", length(.mc)), .mc,
            rep("-p", length(.p)), .p,
            rep("-rt", length(.rt)), .rt,
            rep("-na", length(.na)), .na,
            rep("-ct", length(.ct)), .ct,
            rep("-oq", length(.oq)), .oq,
            rep("-sigma", length(.sigma)), .sigma,
            rep("-softmin_groups", length(.softmin_groups)), .softmin_groups,
            rep("-control_treat_map", length(.control_treat_map)), .control_treat_map,
            " >> ", log_file_dir, "\n"),
            file=paste0("Experiment_",exp.nu, "job"), sep=" ", append=T)
    }
  }
}
#########################################################################################################
#########################################################################################################
#example
Hal_job_writer(exp.nu = 1, cluster.nu = 5, .oo = "GOLABI")
#########################################################################################################
#########################################################################################################
bash_directory <- function(exp.nu, 
                           root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT",
                           source_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                           Read.output=F, number.Of.Clusters=numeric(0), file_name = "bash_directory_"){
  # write a bash script to create the necessary directories and copy the created files to those directories
  
  # exp.nu is the number of the experiment
  # cl.nu is the number of clusters in this experiment
  # root_dir is the direcroty to mount
  # Read.output : if True, writes bash to copy output from veda to computer instead
  # number.Of.Clusters : only need to provide if reading output
  
  cat(c("#!/bin/bash", "\n"),file=paste0(file_name, "exp_", exp.nu), sep="")
  # mount the directory you want to remote:
  cat(c(paste("sshfs", root_dir,
            "~/remote -oauto_cache,reconnect,defer_permissions,noappledouble,negative_vncache,volname=MySSHFSMount"), "\n"),
      file=paste0(file_name, "exp_", exp.nu), sep="", append=T)
  # cd to ~/remote directory
  cat(c("cd ~/remote", "\n"),
      file=paste0(file_name, "exp_", exp.nu), sep="", append=T)
  # copy the directories:
  if (Read.output){
    for(clus in 1:number.Of.Clusters){
      cat(c(paste("cp ","Experiment_", exp.nu,"/Cluster_",clus,"/Outputs/","* ", source_dir, "/Experiment_", exp.nu,"/Cluster_",clus,"/Outputs/", sep=""), "\n"),
          file=paste0(file_name, "exp_", exp.nu),
          sep="", append=T)
    }
  }else{
    cat(c(paste("cp -r ", source_dir, "/Experiment_", exp.nu," .", sep=""), "\n"),
        file=paste0(file_name, "exp_", exp.nu),
        sep="", append=T)
  }

  cat(c("cd ~", "\n"),
      file=paste0(file_name, "exp_", exp.nu), sep="", append=T)
  cat(c("umount", "remote/", "\n"),
      file=paste0(file_name, "exp_", exp.nu), sep=" ", append=T)
}
#########################################################################################################
#########################################################################################################
#example
bash_directory(exp.nu = 2 )
#example output
bash_directory(exp.nu = 7, Read.output = T, number.Of.Clusters = 1)
#########################################################################################################
#########################################################################################################
GEMSTAT_input_constructor <- function(clusterIndex,
                                      gene_expression_Mat, 
                                      .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                      .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                      .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                      .promoter.sequence=Promoter.Sequence.Char,
                                      .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                                      experiment.nu,
                                      .motifs=TF.motifs.count,
                                      .TF_expression=TF.Expression.Quantile.For.GEMSTAT,
                                      .max.enh.length = numeric(0),
                                      ..a=character(0), ..o="DIRECT", coopertingTFIndex=numeric(0), ..i=character(0), ..r=character(0), ..oo=character(0), ..mc=character(0), ..p=character(0), ..rt=character(0), ..na=character(0), ..ct=character(0), ..sigma=character(0)){
  # get the cluster and create all gemstat inputs and copies them to veda, also creates the job file
  
  #This function gets as inout:
  # 1) clusterIndex <- a list where each entry is index of genes in one cluster
  # 2) gene_expression_Mat <- matrix of gene expression which the indexes in clusterIndex refer to. It has to be named
  #  entries for : Function to write the bag of regulatory elements of each gene
  # 3).enhancer.per.gene : a list containing list of enhancers for each gene eg: Vicinity100kb.Enhancer.plusPromoter.By.gene , for all genes: It will be subsetted using gene names
  # 4) .Enhancer.Sequence : matirx of all enhancers: MCFEnhancers
  # 5).promoter.promoter.int : a list containing list of promoters interacting with each gene eg: Promoter.gene.PromoterIntList: for all genes: It will be subsetted using gene names
  # 6) .promoter.sequence : sequence of the promoters of genes in the same order as input1 and input 3--->this should be named by genes, 
  # IT IS IMPORTANT THAT THE Third AND Fifth AND Forth ENTRIES ARE NAMED, (with corresponding genes names)
  # IT IS IMPORTANT THAT THE Third AND Fifth  ENTRIES SHOULD HAVE THE SAME LENGTH, WHICH IS THE NUMBER OF GENES OF INTEREST
  # 7) dir: the directory in which the files will be stored
  # 8) experiment.nu is the number of experiment for which this data is being created
  # 9) .motifs : motifs in the format to be written 
  # 10).TF_expression : TF experession matrix 
  # 11) .max.enh.length : maximum length of considered enhancers
  # 12) last row are gemstat job parameters used in Hal_job_writer:
  # ..a annFile
  # ..o modelOption
  # coopertingTFIndex is a matrix where nrow = number of interacting pairs, ncol= 2 . Each row has the index of interacting TFs in one interaction (index in TF expression matrix)
  # ..i factorInfoFile
  # ..r repressionFile 
  # ..oo objOption 
  # ..mc maxContact 
  # ..p parFile # if you want to use parameter file it has to be created and its address be added here
  # ..rt repressionDistThr 
  # ..na nAlternations 
  # ..ct coopDistThr 
  # ..sigma factorIntSigma
  
  # This function creates: for each cluster: Input and Output folder, populates the Input folder with sequence, gene and TF expression and motif data.
  # This function creates a bash file for copying the created directory to veda
  # This function creates a job file to be run in hal
  # This function uses five other functions: Hal_job_writer, WriteFastaOfBag, ExpressionWriter, MotifWriter, bash_directory
  prev_dir <- getwd()
  setwd(dir = .dir)
  dir.create(paste0(.dir,"/Experiment_", experiment.nu))
  for(clus in 1:length(clusterIndex)){
    print(paste0("Creating cluster ", clus, " input files ..."))
    # creating directory for this cluster
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus))
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Inputs"))
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Outputs"))
    # writing the expression matrix of this cluster
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Inputs", "/Gene_expression"))
    setwd(paste0(.dir,"/Experiment_", experiment.nu,
                 "/Cluster_", clus, "/Inputs", "/Gene_expression"))
    print("Gene Expression")
    current_exp_mat <- gene_expression_Mat[clusterIndex[[clus]],]
    ExpressionWriter(current_exp_mat,
                     output.File.Name = paste("ExpressionMat_",
                                              as.character(clus), sep=""))
    # write the bag of regulatory elements
    print("Sequence")
    setwd(.dir)
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Inputs", "/Sequence"))
    setwd(paste0(.dir,"/Experiment_", experiment.nu,
                 "/Cluster_", clus, "/Inputs", "/Sequence"))
    Index_in_orig <- match(rownames(current_exp_mat), names(.enhancer.per.gene))
    WriteFastaOfBag(enhancer.per.gene=.enhancer.per.gene[Index_in_orig],
                    Enhancer.Sequence=.Enhancer.Sequence,
                    promoter.promoter.int=.promoter.promoter.int[Index_in_orig],
                    promoter.sequence=.promoter.sequence,
                    output.File.Name=paste("RegulatoryElementsSequence", as.character(clus), sep="_"),
                    max.enh.length=.max.enh.length)
    setwd(.dir)
    print("Motifs")
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Inputs", "/Motifs"))
    setwd(paste0(.dir,"/Experiment_", experiment.nu,
                 "/Cluster_", clus, "/Inputs", "/Motifs"))
    MotifWriter(motif.List=.motifs, output.File.Name = "motifs")
    print("TF Expression")
    dir.create(paste0(.dir,"/Experiment_", experiment.nu,
                      "/Cluster_", clus, "/Inputs", "/TF_Expression"))
    setwd(paste0(.dir,"/Experiment_", experiment.nu,
                 "/Cluster_", clus, "/Inputs", "/TF_Expression"))
    ExpressionWriter(.TF_expression ,output.File.Name = "TF_Expression_quantile")
  }
  # write the script for hal job
  setwd(paste0(.dir,"/Experiment_", experiment.nu))

  if (length(coopertingTFIndex) > 0){
    print("Writing the coop file")
    #create the cooperation file and feed to Hal_job_writer 
    TF.RN <- rownames(.TF_expression)
    for(cptf in 1:nrow(coopertingTFIndex)){
      cat(paste(TF.RN[coopertingTFIndex[cptf,1]], TF.RN[coopertingTFIndex[cptf,2]], sep = "\t"), file = "coop.txt", sep = "\n", append = T)
    }
    ..c <- paste0("/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/","Experiment_", experiment.nu,"/coop.txt")
  }else{
    ..c <- character(0)
  }
    
  print("creating the hal job file ...")
  Hal_job_writer(exp.nu=experiment.nu, cluster.nu=length(clusterIndex),
                 seqName="RegulatoryElementsSequence",
                 expressionName="ExpressionMat",
                 Shared_dir="/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/",
                 .a=..a , .o=..o , .c=..c, .i=..i, .r=..r, .oo=..oo, .mc=..mc, .p=..p, .rt=..rt, .na=..na, .ct=..ct, .sigma=..sigma)
  # create a .submit from the written hal job
  sys_com <- paste0(paste0(paste0(paste0("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl"," ","Experiment_",
                                         experiment.nu),"job "),
                           length(clusterIndex)),
                    paste0(paste0(paste0(paste0(" tmpjob_", experiment.nu),
                                         " >job_exp_"), experiment.nu), ".submit"))
  cat(c("#!/bin/bash", "\n"),file=paste0("hal_sub_creator_", "exp_", experiment.nu), sep="")
  cat(sys_com, file=paste0("hal_sub_creator_", "exp_", experiment.nu), sep="", append = T)
  
  sys_com <- paste0("chmod +x ", paste0(" hal_sub_creator_", "exp_", experiment.nu))
  system(sys_com)
  #write a bash file to copy the created files to shared-mounts
  bash_directory(exp.nu=experiment.nu,
                 root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT",
                 source_dir=.dir)
  #make the created bash file executable
  sys_com <- paste0(paste0("chmod +x ", paste0(paste0(.dir,"/Experiment_", experiment.nu),
                                               "/bash_directory_exp_" )), experiment.nu)
  system(sys_com)
  print("Copying to veda ...")
  sys_com <- paste0(paste0(paste0(.dir,"/Experiment_", experiment.nu),
                           "/bash_directory_exp_" ), experiment.nu)
  system(sys_com)
  # reset the working directory to where it was
  setwd(prev_dir)
}

#########################################################################################################
#########################################################################################################
# example
GEMSTAT_input_constructor(clusterIndex=Experiment1.Chosen.Cluster$ClusterIndex,
                          gene_expression_Mat=DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.gt1.ER,
                          .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                          experiment.nu=1)
#or if there is just one cluster
GEMSTAT_input_constructor(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm))),
                          gene_expression_Mat=GSE78167.RNAseq.Avg.ChosenGenes.Norm, 
                          .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                          .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                          .promoter.promoter.int=Promoter.gene.PromoterIntList,
                          .promoter.sequence=Promoter.Sequence.Char,
                          .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                          experiment.nu = 6,
                          .motifs=TF.motifs.TimeSeries.count,
                         .TF_expression=TF.Expression.TimeSeries,
                         .max.enh.length = numeric(0),
                          coopertingTFIndex = rbind(c(1,1),c(1,2)))

#########################################################################################################
#########################################################################################################
Motif.Enrich.Illustrator <- function(Gene.motif.matrix.list=motifHits.ByGene.RegElement,
                                     LogFoldChangeExpression=DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.ER, 
                                     maxLFCgreaterthan=2,
                                     maxLFClessthan=1,
                                     motifName = motif.map$motifName,
                                     outlier = T){
  # Function to compute enrichment of motifs in subsets of dataset, gets as input:
  # 1) Gene.motif.matrix.list :  a list: one entry per gene, the entry is a matrix: each row represents a
  # motif and each column represents a regulatory element: the entries are the number of
  # motif hits in that reg element 
  # 2) LogFoldChangeExpression : log fold change of expression matrix for genes specified in input 1, nrow should be equal to length of input 1
  # 3) maxLFCgreaterthan : First Group of genes have their max absolute log fold change greater than this number
  # 4) maxLFClessthan    : Second group of genes have ther max absolute log fold change less than this number
  # 5) motifName : charachter vector containing the names of the TFs in the same order that appear in other inputs
  # 6) outlier : if True shows outliers in plot
  # output the count of motif hits for each group of genes as a list with two entries, one for each group of genes: a matrix rows: genes, columns: motifs
  # plots as side by side boxplots
  # takes the maximum number of hits over the regulatory elements of a gene as number of hits for that gene

  if(sum(names(Gene.motif.matrix.list) == rownames(LogFoldChangeExpression)) != nrow(LogFoldChangeExpression)){
    stop("names of The Gene.motif.matrix.list and LogFoldChangeExpression are different, they should match")
  }
  Row.max <- function(x){
    return(apply(x, 1, max))
  }
  LogFoldChangeExpression.abs <- abs(LogFoldChangeExpression)
  max.LFC <- Row.max(LogFoldChangeExpression.abs)
  First.group  <- Gene.motif.matrix.list[max.LFC > maxLFCgreaterthan]
  Second.group <- Gene.motif.matrix.list[max.LFC < maxLFClessthan]
  #gather maximums of each motif count per gene in two matrices, one for each group
  aa <- lapply(X = First.group, FUN = Row.max)
  FirstGroup.matrix <- do.call(rbind, aa)
  aa <- lapply(X = Second.group, FUN = Row.max)
  SecondGroup.matrix <- do.call(rbind, aa)
  #plot side by side boxplots
  total.matrix <- matrix(nrow= max(nrow(FirstGroup.matrix), nrow(SecondGroup.matrix)), ncol = 2 * ncol(FirstGroup.matrix))
  total.matrix[1:nrow(FirstGroup.matrix), seq(1,(ncol(total.matrix)-1), 2)] <- FirstGroup.matrix
  total.matrix[1:nrow(SecondGroup.matrix), seq(2,ncol(total.matrix), 2)] <- SecondGroup.matrix
  par(mfrow = c(1,1), mar = c(5,3,1,1))
  plt.col <- numeric(ncol(total.matrix))
  plt.col[seq(1,(ncol(total.matrix)-1), 2)] <- rep(2, ncol(FirstGroup.matrix))
  plt.col[seq(2, ncol(total.matrix), 2)] <- rep(3, ncol(FirstGroup.matrix))
  boxplot.matrix(total.matrix, main = "motif count per group",
                 col = plt.col, ylab="count", xlab = "TFs", xaxt = "n", outline=outlier)
  axis(side = 1, at = seq(1.5, ncol(total.matrix), 2), labels = motifName , las=2)
  abline(v = seq(0.5, ncol(total.matrix), 2), col = 4, lty = 3, lwd = 0.5)
  result = list(First =FirstGroup.matrix, Second =SecondGroup.matrix)
  return(result)
}

#########################################################################################################
#########################################################################################################
# example
aa <- Motif.Enrich.Illustrator(Gene.motif.matrix.list=motifHits.ByGene.RegElement,
                               LogFoldChangeExpression=DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.ER, 
                               maxLFCgreaterthan=2,
                               maxLFClessthan=1)
#########################################################################################################
#########################################################################################################
AnnotToMotCount <- function(annot.df=Gene.ER.Enhancer.annotation.33TFs,
                            ){
  # count the number of motif hits per gene given a threshold and annotationdf
  # a function take as input the annotation dataframe in this format: each row is
  # a motif hit, 6 columns: first column: "gene" Gene names (factor)),
  # second Column "position" position of the hit(integer),
  # third column "motif" motif name (factor)
  # forth column is "LLRscore" : LLRscore of the motif (float)
  # fifth column is "maxLLR" : max LLRscore of the motif (float)
  # sixth column is "reg.el.number" is the number of regulatory element of the gene in which the hit was detected
  # threshold as fraction of maxLLR for that motif: for example: 0.5 means that LLR of the  hits should be at least 0.5 * maxLLR(motif)
}

#########################################################################################################
#########################################################################################################
# example

#########################################################################################################
#########################################################################################################
GEMSTAT_output_Reader <- function(.exp.nu, .number.Of.Clusters,
                                  .root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT",
                                  .dest_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                                  plot.Res = T){
  # this function mounts to veda (root_dir), copies the output file of the given experiments for all of 
  # its clusters to the directory of that experiment on my laptop (dest_dir)
  # It then reads the output expression file, and outputs the expressio, SSE and plots for each gene if needed
  # Inputs:
  # 1) .exp.nu : number of the experiment
  # 2) .number.Of.Clusters = number of clusters in that experiment
  # 3) .root_dir : the root directory in veda where Experiment folders are located
  # 4) .dest_dir : the destination directory on my laptop where experiment folders are located.
  # 5) if you want the output plotted
  #write a bash file to copy the created files to shared-mounts
  prev_dir <- getwd()
  setwd(paste0(.dest_dir, "/Experiment_", .exp.nu))
  bash_directory(exp.nu=.exp.nu,
                 root_dir=.root_dir,
                 source_dir=.dest_dir,
                 Read.output = T,
                 number.Of.Clusters=.number.Of.Clusters,
                 file_name = "bash_directory_output_")
  #make the created bash file executable
  sys_com <- paste0("chmod +x ","bash_directory_output_exp_", .exp.nu)
  system(sys_com)
  print("Copying from veda ...")
  sys_com <- paste0("./bash_directory_output_exp_", .exp.nu)
  system(sys_com)
  AllResults <- list()
  for(clus in 1:.number.Of.Clusters){
    print(paste0("reading output for cluster_",clus,"..."))
    OutPut_clus <- read.table(paste0(.dest_dir, "/Experiment_",
                                    .exp.nu,"/Cluster_",clus,"/Outputs/Experiment_",
                                    .exp.nu,"_Cluster_",clus,".out" ),
                             header = T, sep = "\t")
    cur_Prediction <- OutPut_clus[seq(2, nrow(OutPut_clus), 2), 2:ncol(OutPut_clus)]
    cur_real <- OutPut_clus[seq(1, nrow(OutPut_clus), 2), 2:ncol(OutPut_clus)]
    rownames(cur_Prediction) <- unique(OutPut_clus[, 1])
    rownames(cur_real) <- unique(OutPut_clus[, 1])
    cur_Prediction <- do.call(cbind, cur_Prediction)
    cur_real <- do.call(cbind, cur_real)
    cur_SSE <- rowSums((cur_real - cur_Prediction)^2)/ncol(cur_real)
    AllResults[[clus]] <- list(prediction=cur_Prediction, TrueValue=cur_real, SSE = cur_SSE)
    prev.Par <- par()
    if(plot.Res){
      print("plotting ...")
      par(mfrow = c(10, 10), mar = c(0.1, 0.1, 0.1, 0.1))
      for(gene in 1:nrow(cur_Prediction)){
        plot(cur_Prediction[gene,],
             ylim = c(0,1), col = 2,
             xaxt = "n", yaxt = "n", main ="", ylab = "", xlab = "",
             type = "l", lwd = 2)
        
        lines(cur_real[gene,], ylim = c(0,1), col = 3)
        aa <- cur_SSE[gene]
        text((ncol(cur_Prediction) - 2), 0.9 , format(round(aa, 2), nsmall = 2), cex = 0.5) # write the size of the SSE
      }#end of loop for plotting genes
    }
  }#end of loop over clusters
  par(mfrow = prev.Par$mfrow, mar = prev.Par$mar)
  setwd(prev_dir)
}
#########################################################################################################
#########################################################################################################
# example
Exp8.GEMSTAT.output <- GEMSTAT_output_Reader(.exp.nu = 8, .number.Of.Clusters = 1)
#########################################################################################################
#########################################################################################################
createDescreteMat <- function(columnNames,
                              conditions,
                              replicates = integer(0),
                              contrasts,
                              expMat,
                              adjpvalThresh = 1e-6,
                              lfcThresh = 0,
                              expfData=numeric(0),
                              exppData=numeric(0),
                              use_arrayWeights = T,
                              weightMethod = "genebygene"){
  # creates the descretized expression matrix. 1 for upregulated, -1 for downregulated and 0 for others
  # columnNames is a character vector containing the name of columns in the output matrix
  # conditions is a chatacter vector containing the types of all treatments in order to be used to create the design matrix : IT SHOULD BE ORDERED CORRESPONDING TO THE expMat
  # replicates is an integer vector indicating which condition is from which replicate experiment (to currect for batch effect).  IT SHOULD BE ORDERED CORRESPONDING TO THE expMat. THIS IS OPTIONAL
  # contrasts is a two column dataframe containing the treatments to be contrasted for diff exp calcs. each row is a contrast : the first column is the control condtion , second column is the treated condition.
  # expMat is the expression matrix which will be the input to diff exp method. rows and columns should be named.
  # adjpvalThresh is the threshold on BH adjusted pvalues.
  # lfcThresh is absolute log fold change threshold
  # expfData is the fData of the expressionset object
  # exppData is the pData of the expressionset object
  # use_arrayWeights : if T it uses arrayWeights function.
  # weightMethod : method used for arrayWeights function. defaults to "genebygene". the other method can be "reml"
  
  #construct the matrix to hold the descrete exp mat
  DescreteMat <- matrix(0L, nrow = nrow(expMat), ncol = length(columnNames))
  colnames(DescreteMat) <- columnNames
  rownames(DescreteMat) <- rownames(expMat)
  
  #sort expMat and conditions based on names
  expMat.ordered <- expMat[, sort(conditions, index.return = T)$ix]
  if (length(replicates) > 0){
    replicates.ordered <- replicates[sort(conditions, index.return = T)$ix]
  }
  conditions.Ordered <- sort(conditions)
  
  #create the design matrix
  conditions.Ordered.fac <- factor(conditions.Ordered)
  if (length(replicates) > 0){
    batch <- factor(replicates.ordered)
    DesignMat <- model.matrix(~0+conditions.Ordered.fac+batch)
    colnames(DesignMat) <- c(levels(conditions.Ordered.fac), paste0("batch", c(2:(length(levels(batch))))))
  }else{
    DesignMat <- model.matrix(~0+conditions.Ordered.fac)
    colnames(DesignMat) <- levels(conditions.Ordered.fac)
  }
  if(use_arrayWeights){
    CondWeight <- arrayWeights(expMat.ordered, DesignMat, method = weightMethod)
    #create the fit
    Myfit <- lmFit(expMat.ordered, DesignMat, weights=CondWeight)
  }else{
    CondWeight = numeric(0)
    Myfit <- lmFit(expMat.ordered, DesignMat)
  }

  #create the contrasts
  ContrastPre <- character(0)
  for(i in 1:nrow(contrasts)){
    ContrastPre <- c(ContrastPre, paste(contrasts[i, 2], "-", contrasts[i, 1]))
  }
  myArgs <- c(as.list(ContrastPre))
  myArgs[[length(myArgs) + 1]] <- DesignMat
  names(myArgs)[length(myArgs)] <- "levels"
  #print(names(myArgs))
  Mycontrasts <- do.call(makeContrasts, myArgs)
  # Mycontrasts <- makeContrasts(unlist(ContrastPre),
  #                              levels=DesignMat)
  # final diff exp step
  Mycontr.fit <- eBayes(contrasts.fit(Myfit, Mycontrasts))
  DiffexpRes <- list()
  for(i in 1:nrow(contrasts)){
    DiffexpRes[[i]] <- topTable(Mycontr.fit,
                                coef=i,
                                adjust.method = "BH",
                                p.value = adjpvalThresh,
                                lfc = lfcThresh,
                                number = Inf)
  }
  # modify the descretized matrix
  for(i in 1:length(DiffexpRes)){
    print(paste("number of diff exp probes in contrast ", ContrastPre[i], " : ", nrow(DiffexpRes[[i]])))
    ResMatch <- match(rownames(DiffexpRes[[i]]), rownames(DescreteMat))
    if(nrow(DiffexpRes[[i]]) > 0){
      for(ii in 1:length(ResMatch)){
        DescreteMat[ResMatch[ii], i] <- ifelse(test = (DiffexpRes[[i]]$logFC[ii] > 0), yes = 1, no = -1)
      }
    }
  }
  log_fold_raw <- matrix(nrow = nrow(DescreteMat),
                         ncol = ncol(DescreteMat))
  rownames(log_fold_raw) <- rownames(DescreteMat)
  colnames(log_fold_raw) <- colnames(DescreteMat)
  for(i in 1:length(DiffexpRes)){
    #print(paste("number of diff exp probes in contrast ", ContrastPre[i], " : ", nrow(DiffexpRes[[i]])))
    ResMatch <- match(rownames(DiffexpRes[[i]]),
                      rownames(DescreteMat))
    if(nrow(DiffexpRes[[i]]) > 0){
      for(ii in 1:length(ResMatch)){
        log_fold_raw[ResMatch[ii], i] <- DiffexpRes[[i]]$logFC[ii]
      }
    }
  }
  
  if(length(exppData) > 0){
    exppData <- exppData[sort(conditions, index.return = T)$ix, ]
  }
  return(list(descrete_expMat = DescreteMat,
              designMatrix = DesignMat,
              FullexpMat_ordered = expMat.ordered,
              fDATA = expfData,
              pDATA_ordered = exppData,
              parameters = c(Adj_pval_thresh = adjpvalThresh, lfc_thresh = lfcThresh),
              FitContrasts = Mycontr.fit,
              ARRAYWeights = CondWeight,
              my_conds_ordered = conditions.Ordered,
              my_contrasts = contrasts,
              raw_lfc = log_fold_raw))
}


#########################################################################################################
#########################################################################################################
#example
GSE37386.exp.soft.filt.Descreteaa <- createDescreteMat(columnNames = c("WT_6h", "WT_24h", "GREB1KD_6h", "GREB1KD_24h"),
                                                     conditions = paste(pData(GSE37386.exp.soft.filt)[,42],pData(GSE37386.exp.soft.filt)[,43], pData(GSE37386.exp.soft.filt)[,44] ,sep = "_"),
                                                     contrasts = cbind(c("siNT_control_vehicle", "siNT_control_vehicle", "siGreb_control_vehicle", "siGreb_control_vehicle"),
                                                                       c("siNT_6hrs_E2", "siNT_24hrs_E2", "siGreb_6hrs_E2", "siGreb_24hrs_E2")),
                                                     expMat=exprs(GSE37386.exp.soft.filt),
                                                     adjpvalThresh = 1e-6,
                                                     replicates = c(rep(1,12), rep(2,11)))
#########################################################################################################
#########################################################################################################
CheckAndReplace <- function(UPSorted, DOWNSorted, nu_gn){
  # this function gets the sorted lists and outputs genes not common between the two lists
  # helper function for DiffExpToTopGeneMat function
  aa <- intersect(UPSorted[1:nu_gn], DOWNSorted[1:nu_gn])
  while(length(aa) > 0){
    UPSorted <- UPSorted[!(UPSorted %in% aa)]
    DOWNSorted <- DOWNSorted[!(DOWNSorted %in% aa)]
    aa <- intersect(UPSorted[1:nu_gn], DOWNSorted[1:nu_gn])
  }
  return(list(up = UPSorted[1:nu_gn], down = DOWNSorted[1:nu_gn]))
}
#########################################################################################################
#########################################################################################################
DiffExpToTopGeneMat <- function(topTable_input,
                                nu_contrasts,
                                nu_top_genes,
                                expMatRowNames,
                                mappedID,
                                is_entrez=T,
                                output_colnames
                                #,sort_based_on = "logFC" # for now only "logFC" is implemented
                                )
                                {
  # This function gets the output of a differential expression analysis (input of topTable) and
  # outputs a matrix with each row a gene, each column a constrast, top 'nu_top_genes' upregulated genes (based 
  #  on 'sort_based_on' (# for now only "logFC" is implemented)) are 1,   
  # top 'nu_top_genes' downregulated genes are -1, the genes with lfc between -0.1 and 0.1 are 0. rest are NA
  
  # topTable_input : is the input for  limma's topTable function (output of eBayes)
  # nu_contrasts : is the number of contrast that were used in the differntial expression
  # nu_top_genes : number of genes to take from top up/down regulated genes.
  # expMatRowNames : is the rownames of the expression matrix, in order to be able to map to ENTREZID
  # mappedID : is a vector where each entry is an ID for each row of topTable: the ID can be geneID, symbol, name, ENTREZID, .. ORDER SHOULD CORRESPOND TO expMatRowNames
  # is_entrez : is a boolean, if True it means that the mappedID is ENTREZID, if they are not ENTREZID it runs geneNameToEntrez() function to convert them to ENTREZID
  # output_colnames : is a character vector containing the column names for the output matrix
  # sort_based_on : sorts the differential expression based on this
  # criteria. (# for now only "logFC" is implemented) it can take one of these values: "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B". these values correspond to colnames of the output of topTable function
  
  #convert to ENTREZ ID if not already in that format
  if(!is_entrez){
    map_MappedID <- geneNameToEntrez(mappedID)
    aa <- match(mappedID, map_MappedID[, 1])
    Entrez_Map_ID <- map_MappedID[aa, 2]
    if(sum(is.na(Entrez_Map_ID)) > 0){
      print("warning: following are the gene names that can't  be mapped to ENTREZID, They will be removed from the diff exp results.")
      print(mappedID[is.na(Entrez_Map_ID)])
      unIdentifieble_rownames <- expMatRowNames[is.na(Entrez_Map_ID)]
      
    }
  }else{
    Entrez_Map_ID <- mappedID
    unIdentifieble_rownames <- expMatRowNames[is.na(Entrez_Map_ID)]
  }
  # construct the output matrix
  OutputMat <- matrix(nrow = length(unique(Entrez_Map_ID)), ncol = nu_contrasts)
  rownames(OutputMat) <- unique(Entrez_Map_ID)
  colnames(OutputMat) <- output_colnames
  # create a list which contains top Upregulated genes for each contrast
  Up_list <- list()
  # create a list which contains top Down-regulated genes for each contrast
  Down_list <- list()
  # create a list which contains non-changing genes for each contrast
  No_Change_List <- list()
  
  #extract up and down genes per contrast
  for(cur_contrast in 1:nu_contrasts){
    print(paste("current contrast being examined:", cur_contrast))
    print(output_colnames[cur_contrast])
    cur_top_Table <- topTable(topTable_input,
                            coef = cur_contrast,
                            number = nrow(topTable_input) )
    row_names <- Entrez_Map_ID[match(rownames(cur_top_Table), expMatRowNames)]
    #up
    cur_uplist <-               row_names[sort(cur_top_Table$logFC,
                                               decreasing = T,
                                               index.return = T)$ix]
    cur_uplist <- cur_uplist[!duplicated(cur_uplist)]
    #down
    cur_downlist <-               row_names[sort(cur_top_Table$logFC,
                                                 decreasing = F,
                                                 index.return = T)$ix]
    cur_downlist <- cur_downlist[!duplicated(cur_downlist)]
    
    ConsistantList <- CheckAndReplace(UPSorted = cur_uplist, 
                                      DOWNSorted = cur_downlist,
                                      nu_gn = nu_top_genes)
    
    #not Changing
    cur_nochangelist <- row_names[abs(cur_top_Table$logFC) < 0.1]
    No_Change_List[[cur_contrast]] <- cur_nochangelist
    
    Up_list[[cur_contrast]] <- ConsistantList$up
    OutputMat[match(Up_list[[cur_contrast]], rownames(OutputMat)), cur_contrast] <- 1

    Down_list[[cur_contrast]] <- ConsistantList$down
    OutputMat[match(Down_list[[cur_contrast]], rownames(OutputMat)), cur_contrast] <- -1
    
    OutputMat[match(No_Change_List[[cur_contrast]], rownames(OutputMat)), cur_contrast] <- 0
  
    
  }
  names(Up_list) <-   output_colnames
  names(Down_list) <- output_colnames
  names(No_Change_List) <- output_colnames
  return(list(output_matrix = OutputMat,
              UP_DOWN_NC_Lists = list(UPgenes = Up_list,
                                   DOWNgenes = Down_list,
                                   NoChangeGenes = No_Change_List)))
}


#########################################################################################################
#########################################################################################################
#example
aa <-       DiffExpToTopGeneMat(topTable_input = GSE37386.exp.soft.filt.Descrete$FitContrasts,
                                nu_contrasts = ncol(GSE37386.exp.soft.filt.Descrete$descrete_expMat),
                                nu_top_genes = 100,
                                expMatRowNames = rownames(GSE37386.exp.soft.filt.Descrete$FullexpMat),
                                mappedID = GSE37386.exp.soft.filt.Descrete$fDATA$Entrez_Gene_ID,
                                is_entrez=T,
                                output_colnames = colnames(GSE37386.exp.soft.filt.Descrete$descrete_expMat)
                                #,sort_based_on = "logFC"
                                )

#########################################################################################################
#########################################################################################################
CommonDifExpMat <- function(DifExpMatList){
  # Gets a list where each entry is an output of DiffExpToTopGeneMat, and creates a matrix where the rows are union of all entries' rows, and columns are concatanated
  # DifExpMatList : a list where each entry is a matrix output of DiffExpToTopGeneMat (first entry in output of DiffExpToTopGeneMat)
  # DifExpMatList should be a named list
  
  for(i in 1:length(DifExpMatList)){
    colnames(DifExpMatList[[i]]) <- paste(names(DifExpMatList)[i], colnames(DifExpMatList[[i]]), sep = "_")
  }
  
  #get the union of the genes and sort them
  union_names <- Reduce(union, lapply(DifExpMatList, rownames))
  union_names <- sort(union_names)
  #get the column names
  col_names <- unlist(lapply(DifExpMatList, colnames))
  
  #construct output matrix
  output_mat <- matrix(nrow = length(union_names), ncol = length(col_names))
  rownames(output_mat) <- union_names
  colnames(output_mat) <- col_names
  
  #this is a counter for the column to start writing to, in the following loop
  start_col <- 1
  for(cur_de in 1:length(DifExpMatList)){
    cur_ind <- match(rownames(DifExpMatList[[cur_de]]), union_names)
    if(!all(colnames(DifExpMatList[[cur_de]]) ==
            col_names[start_col:(start_col + ncol(DifExpMatList[[cur_de]]) - 1)])){
      print("WARNING: column names don't match, something is wrong")
      return(-1)
    }
    output_mat[cur_ind, (start_col:(start_col + ncol(DifExpMatList[[cur_de]]) - 1))] <- DifExpMatList[[cur_de]]
    start_col <- start_col + ncol(DifExpMatList[[cur_de]])
  }
  return(output_mat)
}
#########################################################################################################
#########################################################################################################
#example
my_CommonDifExpMat_16 <- CommonDifExpMat(MyDifExpMatList)
#########################################################################################################
#########################################################################################################
QuantileCompare <- function(exp_vec,
                            quantile_vec,
                            assignedVal = c( 0 , 0.33, 0.66, 1)){
  # This function takes expression vector and a quantile vector and outputs a descretized expression vector based on the quantiles. 0 , 0.33, 0.66, 1
  exp_vec_des <- exp_vec
  stopifnot(length(quantile_vec) == (1 + length(assignedVal)))
  for(i in 2:length(quantile_vec)){
    exp_vec_des[(exp_vec <= quantile_vec[i]) & ((exp_vec >= quantile_vec[i-1])) ] <- assignedVal[i-1]
  }
  exp_vec_des[exp_vec >= quantile_vec[length(quantile_vec)]] <- assignedVal[length(quantile_vec)-1]
  return(exp_vec_des)
}
#########################################################################################################
#########################################################################################################
#example
aa <- QuantileCompare(exp_vec = my_Descretized_List$GSE37386$FullexpMat_ordered[,1],
                      quantile_vec = quantile(my_Descretized_List$GSE37386$FullexpMat_ordered[,1]))
#########################################################################################################
#########################################################################################################

TFexpressionDescritizer <- function(createDescreteMat_output, 
                                    TF_index, 
                                    nu_quant = 4,
                                    return_raw_lfc = F, 
                                    use_raw_lfc=T){
  # This Function takes output of createDescreteMat function, and the index of probes corresponding to a TF or multiple TFs, TF_index SHOULD BE NAMED. 
  # outputs a descretized expression vector for the TF. the number of entries of the vectors are two times
  # the number of contrasts (one entry for control, and one entry for treatment of each contrast)
  # In order to do so it averages the expression of each probes across all replicates of a condition,
  # computes the quantiles of the control experiments and assigns 1 / 2 / 3 / 4 to the expression value of
  # the TF considering the quantiles of its control condition.
  # If there are multiple probes associated with the TF, it outputs one vector corresponding to each probe
  
  # createDescreteMat_output : output of createDescreteMat function
  # TF_index : index of probes associated with the TF -- has to be named
  # nu_quant : the number of quantiles to be created
  # return_raw_lfc : logical, if True return raw log fold change
  # use_raw_lfc : if True it uses quantiles the control condition and uses 
  #  the raw lfc values to get the treatment condition expression. HAS TO BE SET TO TRUE, KEPT AS ARGUMENT FOR HISTORICAL REASONS
  
  stopifnot(length(names(TF_index)) > 0,
            all(colnames(createDescreteMat_output$raw_lfc) == colnames(createDescreteMat_output$descrete_expMat)), 
            use_raw_lfc)
  
  my_conds <- createDescreteMat_output$my_conds_ordered
  my_expMat <- createDescreteMat_output$FullexpMat_ordered
  my_contrast <- createDescreteMat_output$my_contrasts
  my_FC <- 2 ^ createDescreteMat_output$raw_lfc[TF_index, ]
  my_discrete <- matrix(nrow = length(TF_index),
                        ncol = 2 * nrow(my_contrast))
  rownames(my_discrete) <- names(TF_index)
  colnames(my_discrete) <- as.character(c(1:ncol(my_discrete)))
  
  for(cur_cont in 1:nrow(my_contrast)){
    #control
    cur_contMat <- my_expMat[, my_conds %in% my_contrast[cur_cont ,1]]
    if(sum(my_conds %in% my_contrast[cur_cont ,1]) == 1){
      cur_contMat_mean <- cur_contMat
      cur_contMat_mean_quantile <- quantile(cur_contMat_mean,
                                            probs = seq(0, 1, length.out = (nu_quant + 1))) 
#      create a discretized version of expression matrix using the quantiles of control experiment
      cur_contMat_mean_desc <- QuantileCompare(exp_vec = cur_contMat_mean,
                                               quantile_vec = cur_contMat_mean_quantile,
                                               assignedVal = seq((1/nu_quant), 1, length.out = nu_quant))
      cur_contMat_mean_desc_tf <- cur_contMat_mean_desc[TF_index]
      
    }else{ # if there are more than one replicates, assign discrete values to each rep and take the median
      cur_qnt_holder <- matrix(nrow = nrow(cur_contMat), 
                               ncol = ncol(cur_contMat))
      for(cur_rep in 1:ncol(cur_contMat)){
        tmp_qnt <- quantile(cur_contMat[, cur_rep],
                            probs = seq(0, 1, length.out = (nu_quant + 1)))
        cur_qnt_holder[, cur_rep] <- QuantileCompare(exp_vec = cur_contMat[, cur_rep],
                                                     quantile_vec = tmp_qnt, 
                                                     assignedVal = seq((1/nu_quant), 1, length.out = nu_quant))
      }
      cur_qnt_holder_median <- apply(cur_qnt_holder, MARGIN = 1, FUN = median)
      cur_contMat_mean_desc_tf <- cur_qnt_holder_median[TF_index]
#      cur_contMat_mean <- rowMeans(cur_contMat) # compute mean expression of each prob in the control condition
    }
    #treatment
    if(use_raw_lfc){
      cur_treatMat_mean_desc_tf <- cur_contMat_mean_desc_tf * my_FC[, cur_cont]
    }else{
      print("This function doesn't support creating TF expression matrix without using raw fold change values")
      return(0)
      # cur_treatMat <- my_expMat[, my_conds %in% my_contrast[cur_cont ,2]]
      # if(sum(my_conds %in% my_contrast[cur_cont ,2]) == 1){
      #   cur_treatMat_mean <- cur_treatMat
      # }else{
      #   cur_treatMat_mean <- rowMeans(cur_treatMat)
      # }
      # # create a discretized version of expression matrix using the quantiles of control experiment
      # cur_treatMat_mean_desc <- QuantileCompare(exp_vec = cur_treatMat_mean,
      #                                           quantile_vec = cur_contMat_mean_quantile,
      #                                           assignedVal = seq((1/nu_quant), 1, length.out = nu_quant))
      # cur_treatMat_mean_desc_tf <- cur_treatMat_mean_desc[TF_index]
    }
    colnames(my_discrete)[(2 * cur_cont - 1):(2 * cur_cont)] <- my_contrast[cur_cont, ]
    my_discrete[, 2 * cur_cont - 1] <- cur_contMat_mean_desc_tf
    my_discrete[, 2 * cur_cont ] <- cur_treatMat_mean_desc_tf
  }
  if(return_raw_lfc){
    return(list(discrete_TF_mat=my_discrete,
                Raw_LFC = createDescreteMat_output$raw_lfc[TF_index, ]))
  }
  return(my_discrete)
}


#########################################################################################################
#########################################################################################################
#example
aaa <- TF.Index.Shrinked.microarray_ENTREZ$AR$GSE37386
names(aaa) <- c("AR_1", "AR_2")
aa <- TFexpressionDescritizer(createDescreteMat_output = my_Descretized_List$GSE37386, TF_index = aaa)
#########################################################################################################
#########################################################################################################
plotExpression <- function(expMat,
                           exportplot=T,
                           filename = "Expression_heatmap.png",
                           .dendrogram = "row",
                           .Rowv = T,
                           .Colv = F,
                           NAto = -3,
                           colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),
                           .distfun = dist,
                           .ColSideColors = character(0),
                           .RowSideColors = character(0)){
  #colorPoints :  points to vary the colors around, up to five
  #it changes all NAs to "NAto" and then does the heatmap
  if(exportplot){
    png(filename,    # create PNG for the heat map        
        width = 8*300,        # 5 x 300 pixels
        height = 8*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  if(sum(is.na(expMat) > 0 )){
    expMat[is.na(expMat)] <- -3
  }
  #plot the heatmap using gplot heatmap.2
  all_my_cols <- c("white","blue", "green", "yellow", "red","black")
  my_cur_colors <- all_my_cols[1:length(colorPoints)]
  my_palette <- colorRampPalette(my_cur_colors)(n = (100*(length(colorPoints) - 1)- 1))
  col_breaks <- numeric(0)
  for(i in 2:length(colorPoints)){
    col_breaks <- c(col_breaks, seq(colorPoints[i-1], colorPoints[i] - 0.01, length=100))
  }
  # col_breaks = c(seq(colorPoints[1]-0.1, colorPoints[1]+0.1, length=100),     # for green
  #                seq(colorPoints[2]-0.1, colorPoints[2]+0.1, length=100),   # for yellow
  #                seq(colorPoints[3]-0.1, colorPoints[3]+0.1, length=100))         # for red
  if(length(.ColSideColors) == 0){
    .ColSideColors <- rep("white", ncol(expMat))
  }
  if(length(.RowSideColors) == 0){
    .RowSideColors <- rep("white", nrow(expMat))
  }
  heatmap.2(expMat,
            #cellnote = format(round(SSEmat, 2), nsmall = 2),  # same data set for cell labels
            main = "Expression matrix", # heat map title
            notecol="black",      # change font color of cell labels to black
            # density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(9,7),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            breaks=col_breaks,    # enable color transition at specified limits
            dendrogram=.dendrogram,     # only draw a .dendrogram dendrogram
            RowSideColors = .RowSideColors,# grouping row-variables(genes) into different categories
            ColSideColors = .ColSideColors,# grouping column-variables(models) into different categories
            Rowv = .Rowv,
            Colv = .Colv,
            scale = "none",
            na.rm=TRUE,
            distfun = .distfun,
            bg = "grey"
  ) 
  if(exportplot){
    dev.off() 
  }
}


#########################################################################################################
#########################################################################################################
#example
plotExpression(my_CommonDifExpMat_16_ERassoc)
#########################################################################################################
#########################################################################################################
SeqTotalLLR <- function(.motifScore_output,
                        seq_range, 
                        .LLR_Thr = 0,
                        .no_thresh = T,
                        .motif_length, 
                        diff_max_LLR=F,
                        pairwise_mat = F, 
                        pairwise_mat_dist=20,
                        LLR2PVal_list =numeric(0),
                        pval_thresh= numeric(0),
                        length_normalize = T){
  # .motifScore_output is the output of motifScore function given a sequence and list of motifs
  # seq_range is an integer vector of length 2 which indicates the start and end positions of our sub_sequence of interest
  # .LLR_Thr is the LLR threshold
  # .no_thresh : if True it doesn't use any thresholding on LLR, instead uses exp(-delta LLR) as the score in each nucleotide position
  # .motif_length : an integer vector containing the length of all motifs in the same order that they appear in .motifScore_output
  # diff_max_LLR : logical. if Ture it uses exp(-(LLR(max) - LLR(cur_seq))) as the score. if False it uses exp(LLR(cur_seq)) as the score
  # LLR2PVal_list : is a list to convert LLR values to p-values. it has an entry per motif present
  # pairwise_mat: a logical, if True the function also produces a matrix for the given range: nrow = number of
  #  TFs, ncol = number of TFs, it stores the number of pairs of TF in a distance threshold. it uses
  #  the LLR2PVal_list to convert LLRs to pvalues and call sites based on the pval_thresh
  # pairwise_mat_dist : an integer specifying the distance between the sites to be called adjacent and counted in pairwise_mat
  # pval_thresh is a numeric either of length 1 or of length equal to the number of TFs. provides a pvalue cutoff for calling motif hits.
  # This function produces a named numeric vector, with one score per TF for the sub_sequence indicated by seq_range
  
  stopifnot(length(seq_range) == 2, 
            all(seq_range[2] - seq_range[1] + 1 >= .motif_length),
            all(names(.motifScore_output$Forward) == names(.motif_length)),
            is.logical(diff_max_LLR),
            is.logical(.no_thresh),
            ((!pairwise_mat) | (length(LLR2PVal_list) > 0)),
            ((!pairwise_mat) | (length(pval_thresh) > 0))
            )
  if(pairwise_mat){
    stopifnot(length(LLR2PVal_list) == length(.motifScore_output$Forward),
              all(names(LLR2PVal_list) == names(.motif_length)))
  }
  nu_TFs <- length(.motifScore_output$Forward)
  TF_names <- names(.motifScore_output$Forward)
  score_output <- numeric(nu_TFs)
  names(score_output) <- TF_names
  full_tot_score <- list()
  
  for(cur_TF in 1:length(score_output)){
    #change the upper bound if it exceeds the current motifs score matrix (since score matrix doesn't have entries for last part of the sequence which is shorter than motif length)
    upper_bnd <- min(seq_range[2], nrow(.motifScore_output$Forward[[cur_TF]]))
    if(.no_thresh){
      #getting the exponential of delta LLRs
      if(diff_max_LLR){
        forw_score <- .motifScore_output$Forward[[cur_TF]][(seq_range[1]:upper_bnd), 4]
        rev_score <-  .motifScore_output$Reverse[[cur_TF]][seq_range[1]:upper_bnd, 4]
      }else{
        forw_score <- exp(.motifScore_output$Forward[[cur_TF]][(seq_range[1]:upper_bnd), 1])
        rev_score <-  exp(.motifScore_output$Reverse[[cur_TF]][seq_range[1]:upper_bnd, 1])
      }

    }else{
      #getting the LLRs greater than threshold
      forw_score <- .motifScore_output$Forward[[cur_TF]][seq_range[1]:upper_bnd, 1]
      forw_score[forw_score < .LLR_Thr] <- 0
      
      rev_score <-  .motifScore_output$Reverse[[cur_TF]][seq_range[1]:upper_bnd, 1]
      rev_score[rev_score < .LLR_Thr] <- 0
    }
    #getting the maximum of rev and forw
    tot_score <- mapply(FUN = max, forw_score, rev_score)
    
    #setting the parts shorter than motif length to zero
    if(seq_range[2] > upper_bnd){
      tot_score <- c(tot_score, integer(seq_range[2] - upper_bnd))
    }
    tot_score[(length(tot_score) - .motif_length[cur_TF] + 2): length(tot_score)] <- 0
    full_tot_score[[cur_TF]] <- tot_score
    if(length_normalize){
      cur_TF_score <- (sum(tot_score)/length(tot_score)) * 1000
    }else{
      cur_TF_score <- sum(tot_score)
    }
    score_output[cur_TF] <- cur_TF_score
  }
  results <- list(TF_Score = score_output)
  
  if(pairwise_mat){
    pairwise_matrix <- matrix(0L,nrow= nu_TFs, 
                              ncol= nu_TFs)
    overlap_mat <- matrix(0L,nrow= nu_TFs, 
                          ncol= nu_TFs)
    nonoverlap_mat <- matrix(0L,nrow= nu_TFs, 
                             ncol= nu_TFs)
    colnames(pairwise_matrix) <- TF_names
    rownames(pairwise_matrix) <- TF_names
    colnames(overlap_mat) <- TF_names
    rownames(overlap_mat) <- TF_names
    colnames(nonoverlap_mat) <- TF_names
    rownames(nonoverlap_mat) <- TF_names
    min_LLR <- numeric(nu_TFs)
    names(min_LLR) <- TF_names
    if(length(pval_thresh) == 1){
      pval_thresh <- rep(pval_thresh, nu_TFs)
      names(pval_thresh) <- TF_names
    }else if(length(pval_thresh) == nu_TFs){
      names(pval_thresh) <- TF_names
    }else{
      print("the length of pval_thresh provided is not equal to '1' or the number of TFs")
      stop()
    }
    stopifnot(all(names(LLR2PVal_list) == TF_names)
              #,all(.motif_length < (pairwise_mat_dist - 1))
              )
    for(i in 1:nu_TFs){
      aa <- which(LLR2PVal_list[[i]][,2] > pval_thresh[i])[1]
      if(aa == 1){
        min_LLR[i] <- LLR2PVal_list[[i]][aa, 1] / 1000
      }else{
        min_LLR[i] <- LLR2PVal_list[[i]][(aa-1), 1] / 1000
      }
      
    }
    #print("LLR threshold")
    #print(min_LLR)
    
    full_tot_mat <- do.call(cbind, full_tot_score)
    colnames(full_tot_mat) <- TF_names

    for(i in 1:nu_TFs){
      full_tot_mat[(full_tot_mat[, i] < min_LLR[i]), i] <- 0
      full_tot_mat[(full_tot_mat[, i] > 0), i] <- 1
    }
    
    # create a matrix listing the sites that appear, each row corresponds to a site, first column is the
    #  TF name, second column is the starting pos, third column is the ending pos
    site_list_mat <- matrix(nrow =0 , ncol = 3)
    colnames(site_list_mat) <- c("TF", "Start", "End")
    aa_nzrows <- which(! rowSums(full_tot_mat) %in% 0)
    for(crow in aa_nzrows){
      aa_tf <- which(full_tot_mat[crow, ] > 0)
      aa_start <- rep(crow, length(aa_tf))
      aa_end <- aa_start + .motif_length[aa_tf] - 1
      if(length(aa_tf) == 1){
        site_list_mat <- rbind(site_list_mat, 
                               c(aa_tf, aa_start, aa_end))
      }else{
        
        site_list_mat <- rbind(site_list_mat,
                               do.call(cbind, list(aa_tf, aa_start, aa_end)))
      }
    }
    # print(site_list_mat)
    # print(nrow(site_list_mat))
    my_dynamic_mat <- matrix(nrow = 0, ncol = 3)
    colnames(my_dynamic_mat) <- c("TF", "Start", "End")
    for(crow in 1:nrow(site_list_mat)){
      # print("ndim(my_dynamic_mat)_after before_drop")
      # print(dim(my_dynamic_mat))
      if(nrow(site_list_mat) == 0){
        print("There are no sites present for this sequence")
        break()
      }
      cur_drop <- which(my_dynamic_mat[, 3] < (site_list_mat[crow, 2] - pairwise_mat_dist))
      # print("cur_drop")
      # print(cur_drop)
      if(length(cur_drop) > 0 ){
        my_dynamic_mat <- my_dynamic_mat[-cur_drop,]
        if(length(my_dynamic_mat) == 3){
          my_dynamic_mat <- matrix(my_dynamic_mat, nrow = 1)
        }
      }
      
      # print("ndim(my_dynamic_mat)_after cur_drop")
      # print(dim(my_dynamic_mat))
      cur_ol <- my_dynamic_mat[(my_dynamic_mat[, 3] >=  site_list_mat[crow, 2]), 1]
      
      cur_all <- my_dynamic_mat[, 1]
      
      my_dynamic_mat <- rbind(my_dynamic_mat, site_list_mat[crow,])
      # print("my_dynamic_mat")
      # print(my_dynamic_mat)
      # print("cur_ol")
      # print(table(cur_ol))
      # print("cur_all")
      # print(table(cur_all))
      if(length(cur_ol) > 0){
        cur_ol_table <- table(cur_ol)
        overlap_mat[site_list_mat[crow, 1], as.integer(names(cur_ol_table))] <- overlap_mat[site_list_mat[crow, 1], as.integer(names(cur_ol_table))] + cur_ol_table
        overlap_mat[setdiff(as.integer(names(cur_ol_table)), site_list_mat[crow, 1]), site_list_mat[crow, 1]] <- overlap_mat[setdiff(as.integer(names(cur_ol_table)), site_list_mat[crow, 1]), site_list_mat[crow, 1]] + cur_ol_table[! names(cur_ol_table) %in% as.character(site_list_mat[crow, 1])]
      }
      if(length(cur_all) > 0){
        cur_all_table <- table(cur_all)
        pairwise_matrix[site_list_mat[crow, 1], as.integer(names(cur_all_table))] <- pairwise_matrix[site_list_mat[crow, 1], as.integer(names(cur_all_table))] + cur_all_table
        pairwise_matrix[setdiff(as.integer(names(cur_all_table)), site_list_mat[crow, 1]), site_list_mat[crow, 1]] <- pairwise_matrix[setdiff(as.integer(names(cur_all_table)), site_list_mat[crow, 1]), site_list_mat[crow, 1]] + cur_all_table[! names(cur_all_table) %in% as.character(site_list_mat[crow, 1])]
      }
      
    }
    nonoverlap_mat <- pairwise_matrix - overlap_mat
    results[[2]] <- overlap_mat
    results[[3]] <- nonoverlap_mat
    results[[4]] <- pairwise_matrix
    names(results)[2:4] <- c("overlap_mat", "nonoverlap_mat", "pairwise_matrix")
  }
  return(results)
}
#########################################################################################################
#########################################################################################################
#example
aa <- SeqTotalLLR(.motifScore_output = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore[[1]],
                  seq_range = c(400, 1000),
                  .LLR_Thr = 0,
                  .no_thresh = T,
                  .motif_length = unlist(lapply(X = TF.motifs.Shrinked.t, ncol)))


aa <- SeqTotalLLR(.motifScore_output = aa_Positive_set_motifscore[[1]][[1]],
                  seq_range = c(1, nchar(as.character(Positive_set_seq_list[[1]][1]))),
                  .LLR_Thr = 0,
                  .no_thresh = T,
                  .motif_length = unlist(lapply(X = TF.motifs.Expanded_t, ncol)),
                  diff_max_LLR = F, 
                  length_normalize = F,
                  pairwise_mat = T,
                  pairwise_mat_dist = 25,
                  LLR2PVal_list = TF.motifs.Expanded_LLR_to_pval_pseudo,
                  pval_thresh = 0.0001)
# View(aa$overlap_mat)
# aaaaa <- aa$overlap_mat[upper.tri(aa$overlap_mat, diag = T)]
# names(aaaaa) <- aa_names2
# str(aa$overlap_mat[upper.tri(aa$overlap_mat, diag = T)])
# aa_names2 <- character(0)
# for(aai in 1:length(TF.motifs.Expanded_t)){
#   aa_names <- paste(names(TF.motifs.Expanded_t)[aai], names(TF.motifs.Expanded_t)[1:aai],sep = "_vs_")
#   aa_names2 <- c(aa_names2, aa_names)
# }
#########################################################################################################
#########################################################################################################
EnhancerChopper_rangeMat <- function(enh_len, piece_len, .step_size, .must_contain=numeric(0)){
  # creates a range matrix for the chopped sequences
  # enh_len : length of the enhancer
  # piece_len : length of each piece
  # .step_size : size of steps
  # [optional] .must_contain is a matrix that each row corresponds to a site (1st column start pos, 2nd column end pos).
  # we want at least one of the rows to be present in all of the chopped ranges.
  
  stopifnot(enh_len >= piece_len)

  extraNuc <- enh_len - piece_len
  rangeRow <- ceiling(extraNuc / .step_size)
  rangeMat <- matrix(nrow = (rangeRow + 1), ncol = 2)
  rangeMat[1, ] <- c(1, piece_len)
  for(i in 1:(rangeRow)){
    max_len <- min((rangeMat[i, 2] + .step_size), enh_len)
    rangeMat[i+1, ] <- c((rangeMat[i, 1] + .step_size), max_len)
  }
  if(length(.must_contain) > 0){ # if there is a filtering criteria, apply that
    filt <- logical(nrow(rangeMat))
    for(i in 1:nrow(.must_contain)){
      col_1 <- (rangeMat[, 1] <= .must_contain[i, 1])
      col_2 <- (rangeMat[, 2] >= .must_contain[i, 2])
      filt <- filt | (col_1 & col_2)
    }
    print(paste(sum(filt), "out of", length(filt), "ranges remained after .must_contain filtering"))
    rangeMat <- rangeMat[filt, ]
    if(sum(filt) == 1){
      rangeMat <- t(as.matrix(rangeMat))
    }
  }
  return(rangeMat)
}
#########################################################################################################
#########################################################################################################
#example
aa <- EnhancerChopper_rangeMat(1000, 500, 200)
#########################################################################################################
#########################################################################################################
EnhancerChopper <- function(Enhancer_seq,
                            EnhancerGR,
                            motifScore_output,
                            motifList,
                            no_thresh = T,
                            LLR_thresh = 0,
                            piece_length = 1000,
                            step_size = 250,
                            must_contain = numeric(0),
                            .diff_max_LLR=T, 
                            pairwise_closeby_site_freq = F,
                            close_thr=20, 
                            LLR2PVal_lists =numeric(0),
                            pval_threshold= numeric(0),
                            normalize_by_length = T ){
  # Enhancer_seq :  a named character vector containing the sequence of enhancers
  # EnhancerGR : GRanges object of enhancers
  # motifScore_output : named list of outputs of motifScore function given the 'Enhancer_seq' sequences and some motifs.
  # motifList : the named list of motifs used to generate motifScore_output
  # no_thresh : if True doesn't use the LLR threshold and instead adds up exp((LLRmax - LLR))
  # LLR_thresh : is the threshold on LLRs that get summed up for each piece of sequence
  # piece_length : is the length of the chopped pieces 
  # step_size : is the step size of the sliding window creating the choped enhancers.
  # must_contain is a list where each entry corresponds to an enhancer. it cantains a matrix where each row 
  #  corresponds to a site. first column is the start position, last column is the end position. at least one of the rows 
  #  for each enhancer should be present in all of the output chopped enhancers.
  # .diff_max_LLR : logical. if Ture it uses exp(-(LLR(max) - LLR(cur_seq))) as the score. if False it uses exp(LLR(cur_seq)) as the score
  # pairwise_closeby_site_freq: logical. if True it also outputs a nu_motif * nu_motif matrix containing the number of closeby pairs (closer than close_thr) for each sequence
  # close_thr: is an integer specifying the distance used in pairwise_closeby_site_freq.
  # normalize_by_length: if True it normalizes the scores by the length of the sequence
  # pval_threshold, LLR2PVal_lists : read description in SeqTotalLLR function
  
  # it outputs a list with three / five(if pairwise_closeby_site_freq) entries: first one contains the chopped sequences for each input seq and the second one
  # is a list of length equal to length of 'Enhancer_seq'. each entry is
  # a matrix  where each row corresponds to a chopped piece of the parent sequence,  each
  # column corresponds to a TF. entries of the matrix correspond to either "sum of LLRs above LLR_thresh of
  # nucleutides in that piece (forward or reverse)" or "sum of exp((LLRmax - LLR)) for nucleotides". and third entry is
  # the GRanges object associated with the chopped enhancers. 
  # if(pairwise_closeby_site_freq) There will be a  forth entry which is like the second entrybut columns are different. it is a
  #  list of length equal to length of 'Enhancer_seq'. each entry is a matrix  where each row corresponds to a chopped piece of
  #  the parent sequence,  each column corresponds to a pair of TFs, the entry shows the number of sites of the two TFs that are closer than close_thr
  # if(pairwise_closeby_site_freq) There will be a  fifth entry which is the same as 4th entry but only shows overlapping sites.
  
  
  
  stopifnot(is.character(Enhancer_seq),
            length(names(Enhancer_seq)) == length(Enhancer_seq), 
            length(motifScore_output) == length(Enhancer_seq),
            
            all(names(Enhancer_seq) == names(motifScore_output)),
            length(names(motifList)) == length(motifList) ,
            all(unlist(lapply(motifList, nrow)) %in% 4), # insuring that all pwms have 4 rows and number of bases columns
            length(unique(unlist(lapply(lapply(motifScore_output, "[[", 1), length)))) == 1,
            all(names(motifScore_output[[1]][[1]]) == names(motifList))
            )
  # construct the output list
  output_list <- list()
  chopped_Seq_all <- list() #this list holds the sequence of chopped enhancers, each entry corresponds to one entry of Enhancer_seq, and contains a named character vector of all chopped seqs
  chopped_GRanges <- list() # this list has one entry for each input enhancer, entry contains the GRanges object correponding to the chopped sub enhancers of that enhancer
  names(EnhancerGR) <- names(Enhancer_seq)
  output_list_adjacency <- list()
  output_list_overlap <- list()
  if(pairwise_closeby_site_freq){
    comb_names <- character(0)
    for(ctf in 1:length(motifList)){
      aa_names <- paste(names(motifList)[ctf], names(motifList)[1:ctf],sep = "_vs_")
      comb_names <- c(comb_names, aa_names)
    }
  }
  
  for(cur_seq in 1:length(Enhancer_seq)){
    print(paste("seq number", cur_seq, "out of", length(Enhancer_seq), ":"))
    if(nchar(Enhancer_seq[cur_seq]) <= piece_length){
      # if there is no need for chopping since the sequence is shorter than specified piece length
      print("not chopped.")
      # construct the output matrix for this sequence
      output_list[[cur_seq]] <- matrix(nrow = 1, ncol = length(motifList))
      colnames(output_list[[cur_seq]]) <- names(motifList)
      rownames(output_list[[cur_seq]]) <- paste(names(Enhancer_seq)[cur_seq], "1", sep = "_")
      aa_cur_res <- SeqTotalLLR(.motifScore_output = motifScore_output[[cur_seq]],
                                seq_range = c(1, nchar(Enhancer_seq[cur_seq])),
                                .LLR_Thr = LLR_thresh,
                                .no_thresh = no_thresh,
                                .motif_length = unlist(lapply(X = motifList, ncol)),
                                diff_max_LLR = .diff_max_LLR, 
                                length_normalize = normalize_by_length,
                                pairwise_mat = pairwise_closeby_site_freq,
                                pairwise_mat_dist = close_thr,
                                LLR2PVal_list = LLR2PVal_lists,
                                pval_thresh = pval_threshold)
      # "overlap_mat", "nonoverlap_mat", "pairwise_matrix"
      if(pairwise_closeby_site_freq){
        output_list_adjacency[[cur_seq]] <- matrix(nrow = 1, ncol = (length(motifList) * (length(motifList) + 1) / 2))
        output_list_overlap[[cur_seq]] <- matrix(nrow = 1, ncol = (length(motifList) * (length(motifList) + 1) / 2))
        colnames(output_list_adjacency[[cur_seq]] ) <- comb_names
        colnames(output_list_overlap[[cur_seq]] ) <- comb_names
        output_list_adjacency[[cur_seq]][1, ] <-  aa_cur_res$pairwise_matrix[upper.tri(aa_cur_res$pairwise_matrix, diag = T)]
        output_list_overlap[[cur_seq]][1, ] <- aa_cur_res$overlap_mat[upper.tri(aa_cur_res$overlap_mat, diag = T)]
      }
      output_list[[cur_seq]][1, ] <- aa_cur_res[[1]]
      chopped_Seq_all[[cur_seq]] <- Enhancer_seq[cur_seq]
      names(chopped_Seq_all[[cur_seq]]) <- rownames(output_list[[cur_seq]])
      chopped_GRanges[[cur_seq]] <- EnhancerGR[cur_seq]
      names(chopped_GRanges[[cur_seq]]) <- rownames(output_list[[cur_seq]])
    }else{#if the sequence has to be chopped
      # create a range matrix for pieces
      if(length(must_contain) > 0){
        print("using 'must_contain'")
        cur_rangeMat <- EnhancerChopper_rangeMat(enh_len = nchar(Enhancer_seq[cur_seq]),
                                                 piece_len = piece_length,
                                                 .step_size = step_size,
                                                 .must_contain = must_contain[[cur_seq]])
      }else{
        cur_rangeMat <- EnhancerChopper_rangeMat(enh_len = nchar(Enhancer_seq[cur_seq]),
                                                 piece_len = piece_length,
                                                 .step_size = step_size)
      }
      print(paste("chopped into", nrow(cur_rangeMat), "pieces."))
      # construct the output matrix for this sequence
      output_list[[cur_seq]] <- matrix(nrow = nrow(cur_rangeMat), ncol = length(motifList))
      colnames(output_list[[cur_seq]]) <- names(motifList)
      rownames(output_list[[cur_seq]]) <- paste(names(Enhancer_seq)[cur_seq], c(1:nrow(cur_rangeMat)), sep = "_")
      # initialize the chopped sequence list
      chopped_Seq_all[[cur_seq]] <- character(nrow(cur_rangeMat))
      names(chopped_Seq_all[[cur_seq]]) <- rownames(output_list[[cur_seq]])
      if(pairwise_closeby_site_freq){
        output_list_adjacency[[cur_seq]] <- matrix(nrow = nrow(cur_rangeMat), 
                                                   ncol = (length(motifList) * (length(motifList) + 1) / 2))
        output_list_overlap[[cur_seq]] <- matrix(nrow = nrow(cur_rangeMat), 
                                                   ncol = (length(motifList) * (length(motifList) + 1) / 2))
        colnames(output_list_adjacency[[cur_seq]] ) <- comb_names
        colnames(output_list_overlap[[cur_seq]] ) <- comb_names
      }
      
      
      for(cur_piece in 1:nrow(cur_rangeMat)){
        print(paste("piece no", cur_piece))
        # compute the score of the peice for each motif
        aa_cur_Res <-  SeqTotalLLR(.motifScore_output = motifScore_output[[cur_seq]],
                                   seq_range = c(cur_rangeMat[cur_piece, 1], cur_rangeMat[cur_piece, 2]),
                                   .LLR_Thr = LLR_thresh,
                                   .no_thresh = no_thresh,
                                   .motif_length = unlist(lapply(X = motifList, ncol)),
                                   diff_max_LLR = .diff_max_LLR,
                                   length_normalize = normalize_by_length,
                                   pairwise_mat = pairwise_closeby_site_freq,
                                   pairwise_mat_dist = close_thr,
                                   LLR2PVal_list = LLR2PVal_lists,
                                   pval_thresh = pval_threshold)
        if(pairwise_closeby_site_freq){
          output_list_adjacency[[cur_seq]][cur_piece, ] <-  aa_cur_res$pairwise_matrix[upper.tri(aa_cur_res$pairwise_matrix, diag = T)]
          output_list_overlap[[cur_seq]][cur_piece, ] <- aa_cur_res$overlap_mat[upper.tri(aa_cur_res$overlap_mat, diag = T)]
          
          #output_list_adjacency[[cur_seq]][[cur_piece]] <-  aa_cur_Res[2:4]
        }
        
        output_list[[cur_seq]][cur_piece, ] <- aa_cur_Res[[1]]
        if(sum(is.na(output_list[[cur_seq]][cur_piece, ])) > 0){
          print("WARNINING, There are NA scores, find out why")
        }
        # add the sequence of the piece to the list of chopped sequences
        chopped_Seq_all[[cur_seq]][cur_piece] <- substr(Enhancer_seq[cur_seq], cur_rangeMat[cur_piece, 1], cur_rangeMat[cur_piece, 2])
        
        if(cur_piece == 1){      # initialize the chopped GRanges
          cur_Grng <- as.data.frame(EnhancerGR[cur_seq])
          cur_Grng$end[cur_piece] <-  cur_Grng$start[cur_piece] + cur_rangeMat[cur_piece, 2] - 1
        }else{
          cur_Grng <- rbind(cur_Grng, as.data.frame(EnhancerGR[cur_seq]))
          cur_Grng$end[cur_piece] <- cur_Grng$start[cur_piece] + cur_rangeMat[cur_piece, 2] - 1
          cur_Grng$start[cur_piece] <- cur_Grng$start[cur_piece] + cur_rangeMat[cur_piece, 1] - 1
        }
      }
      chopped_GRanges[[cur_seq]] <-  makeGRangesFromDataFrame(cur_Grng)
      names(chopped_GRanges[[cur_seq]]) <- rownames(output_list[[cur_seq]])
    }
  } # end of loop over sequences
  names(chopped_Seq_all) <- names(Enhancer_seq)
  names(output_list) <- names(Enhancer_seq)
  names(chopped_GRanges) <- names(Enhancer_seq)
  if(pairwise_closeby_site_freq){
    names(output_list_adjacency) <- names(Enhancer_seq)
    names(output_list_overlap) <- names(Enhancer_seq)
  }

  return(list(Chopped_Seq = chopped_Seq_all,
              Chopped_Score = output_list,
              Chopped_GRanges=chopped_GRanges,
              Adjan_list = output_list_adjacency,
              overlap_list = output_list_overlap))
}

#########################################################################################################
#########################################################################################################
#example
aa <- EnhancerChopper(Enhancer_seq = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence[1:2],
                      motifScore_output = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore,
                      EnhancerGR = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$GRange[1:2,],
                      motifList = TF.motifs.Shrinked.t,
                      no_thresh = T,
                      LLR_thresh = 0,
                      piece_length = 100,
                      step_size = 50 )
#########################################################################################################
#########################################################################################################
# Needs to be updated.
# has not been modified after adding the adjancent site capability to EnhancerChopper
EnhancerChopper_wrapper <- function(Entrez_IDs,
                                    .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                    .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                    .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                    .promoter.sequence=Promoter.Sequence.Char,
                                    #.Enhancer_seq=character(0),
                                    #.EnhancerGR=character(0),
                                    unique_motifScore,
                                    .motifList,
                                    .no_thresh = T,
                                    .LLR_thresh = 0,
                                    .piece_length = 1000,
                                    .step_size = 250,
                                    ..diff_max_LLR,
                                    dimer_filter_ER=F,
                                    dimer_list=ER.associated.reg.elements_uniq_pair_Concat_sort_un,
                                    dimer_filter_ER_LLR_th=2, 
                                    dimer_filter_corresponding_seq_nonProcessed=ER.associated.reg.elements$sequence,
                                    GRange_overlap=F,
                                    GRange_overlap_object=ReMapChIP.GRanges.list[52],
                                    GRange_overlap_percent=50){

  # this function gets a character vector of entrez ID of the genes, creates the set of enhancers associated with each gene.
  # chops the obtained enhancers into overlapping pieces of specified length and calculates the affinity of each chopped piece for each TF
  # if specified: filters enhancers for overlap and or interaction with specified GRanges, filters chopped pieces for specific positions of ER-ER dimers
  
  # Entrez_IDs : is a character vector of Entrez ID of the genes of interest
  # .enhancer.per.gene= list of enhancers assigned to each gene. list contains one entry per gene.
  # .Enhancer.Sequence : sequence og enhancers, the order is the same as indexed .enhancer.per.gene
  # .promoter.promoter.int a list containing the index of promoters interacting with each promoter.
  # .promoter.sequence : character vector of the sequence of all promoters. ordered the same as .promoter.promoter.int

  # .Enhancer_seq: a named character vector containing the sequence of enhancers
  # .EnhancerGR : GRanges object of enhancers
  # unique_motifScore : named list of outputs of motifScore function given the 'Enhancer_seq' sequences and some motifs.
  # .motifList : the named list of motifs used to generate motifScore_output
  # .no_thresh : if True doesn't use the LLR threshold and instead adds up either exp((LLRmax - LLR)) or exp( LLR)
  # .LLR_thresh : is the threshold on LLRs that get summed up for each piece of sequence
  # .piece_length : is the length of the chopped pieces 
  # .step_size : is the step size of the sliding window creating the choped enhancers.
  # .must_contain is a list where each entry corresponds to an enhancer. it cantains a matrix where each row 
  #   corresponds to a site. first column is the start position, last column is the end position. at least one of the rows 
  #   for each enhancer should be present in all of the output chopped enhancers.
  # ..diff_max_LLR : logical. if Ture it uses exp(-(LLR(max) - LLR(cur_seq))) as the score. if False it uses exp(LLR(cur_seq)) as the score
  # dimer_filter_ER : Filter chopped enhancers based on the position of precalculated dimer ER sites
  # dimer_list : the list containing the ER dimer sites for all enhancers for multiple thresholds e.g. : ER.associated.reg.elements_uniq_pair_Concat_sort_un
  # dimer_filter_ER_LLR_th : the LLR threshold above which the sites are forced to be present in the chopped pieces.
  # GRange_overlap : logical. if True filters the results for only the regulatory elements that overlap with the provided "GRange_overlap_object" at least by "GRange_overlap_Percent"
  # GRange_overlap_object : is a list containing one entry. that entry is a GRanges object which we want to filter the results based in overlap with it.
  
  # Gather the regulatory elements for each gene.
  my_gene_index <- match(Entrez_IDs, names(.enhancer.per.gene))
  
  all_reg_elements <- WriteFastaOfBag(enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene[my_gene_index],
                                      Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                      promoter.promoter.int=Promoter.gene.PromoterIntList[my_gene_index],
                                      promoter.sequence=Promoter.Sequence.Char,
                                      return.Length = T,
                                      returnListNotwrite = T)
  # name each element 
  for(cur_gene in 1:length(all_reg_elements$sequence)){
    names(all_reg_elements$sequence[[cur_gene]]) <- character(length(all_reg_elements$sequence[[cur_gene]]))
    names(all_reg_elements$length[[cur_gene]]) <- character(length(all_reg_elements$sequence[[cur_gene]]))
    names(all_reg_elements$GRange[[cur_gene]]) <- character(length(all_reg_elements$sequence[[cur_gene]]))
    for(cur_el in 1:length(all_reg_elements$sequence[[cur_gene]])){
      names(all_reg_elements$sequence[[cur_gene]])[cur_el] <- as.character(cur_el)
      names(all_reg_elements$length[[cur_gene]])[cur_el] <- as.character(cur_el)
      names(all_reg_elements$GRange[[cur_gene]])[cur_el] <- as.character(cur_el)
    }
  }
  
  # unlisting all gene's entries into one list
  all_reg_elements_unlisted <- list()
  all_reg_elements_unlisted$sequence <- unlist(all_reg_elements$sequence)
  all_reg_elements_unlisted$sequence <- lapply(all_reg_elements_unlisted$sequence, toupper)
  
  ### debug
  if(length(names(all_reg_elements_unlisted$sequence))==0){
    print("names of the sequences disappeared after toupper function, fix it.")
  }
  ### debug
  
  all_reg_elements_unlisted$length <- unlist(all_reg_elements$length)
  aa <- all_reg_elements$GRange[[1]]
  for(i in 2:length(all_reg_elements$GRange)){
    aa <- c(aa, all_reg_elements$GRange[[i]])
  }
  all_reg_elements_unlisted$GRange <- aa
  remove(aa)
  
  # removing duplicate entries before computing the affinity for TFs
  all_reg_elements_unlisted_unique <- list()
  all_reg_elements_unlisted_unique$length <- all_reg_elements_unlisted$length[! duplicated(all_reg_elements_unlisted$sequence)]
  all_reg_elements_unlisted_unique$sequence <- all_reg_elements_unlisted$sequence[! duplicated(all_reg_elements_unlisted$sequence)]
  all_reg_elements_unlisted_unique$GRange <- all_reg_elements_unlisted$GRange[! duplicated(all_reg_elements_unlisted$sequence)]
  
  # scanning for motifs
  all_reg_elements_unlisted_unique$MotifScore <- list()
  if(length(unique_motifScore) > 0){
    stopifnot(all(names(unique_motifScore) == names(all_reg_elements_unlisted_unique$sequence)))
    all_reg_elements_unlisted_unique$MotifScore <- unique_motifScore
  }else{
    for(cur_seq in 1:length(all_reg_elements_unlisted_unique$sequence)){
      print(paste("scanning sequence number", cur_seq, "out of", length(all_reg_elements_unlisted_unique$sequence), "for motifs"))
      all_reg_elements_unlisted_unique$MotifScore[[cur_seq]] <- MotifScore(seq = all_reg_elements_unlisted_unique$sequence[cur_seq],
                                                                           bg = c(0.25, 0.25, 0.25, 0.25),
                                                                           motifList = .motifList)
    }
    names(all_reg_elements_unlisted_unique$MotifScore) <- names(all_reg_elements_unlisted_unique$sequence)
  }

  
  # Filtering for dimer sites
  if(dimer_filter_ER){
    print("Filtering for dimer sites ...")
    aa <- do.call(c,
                  dimer_filter_corresponding_seq_nonProcessed)
    aaa <- aa[!duplicated(aa)]
    aaa2 <- match(all_reg_elements_unlisted_unique$sequence,
                  aaa)
    #get the list of motifhits of ER monomers at most 3 bp apart, for each sequence of interest
    all_reg_elements_unlisted_unique$dimerPos <- dimer_list[[dimer_filter_ER_LLR_th]]$UnsortedOutput$ESR1_ESR1[aaa2]
    names(all_reg_elements_unlisted_unique$dimerPos) <- names(all_reg_elements_unlisted_unique$sequence)
    
    ####   ####   #### small helper function

    ####   ####   #### 
    print(length(all_reg_elements_unlisted_unique$dimerPos))
    all_reg_elements_unlisted_unique$dimerPos_range <- lapply(all_reg_elements_unlisted_unique$dimerPos, aaToRange)
    remove(aa)
    remove(aaa)
    remove(aaa2)
  }else{
    all_reg_elements_unlisted_unique$dimerPos_range <- numeric(0)
  }

  # chop and score unique regulatory elements
  all_reg_elements_unlisted_unique$Chopped <- EnhancerChopper(Enhancer_seq = all_reg_elements_unlisted_unique$sequence,
                                                              EnhancerGR = all_reg_elements_unlisted_unique$GRange,
                                                              motifScore_output = all_reg_elements_unlisted_unique$MotifScore,
                                                              motifList = .motifList,
                                                              no_thresh = .no_thresh,
                                                              LLR_thresh = .LLR_thresh,
                                                              piece_length = .piece_length,
                                                              step_size = .step_size,
                                                              must_contain= all_reg_elements_unlisted_unique$dimerPos_range,
                                                              .diff_max_LLR = ..diff_max_LLR)
  # map unique scores and chopped enhancers to the original non-unique set
  aaa <- match(all_reg_elements_unlisted$sequence,
               all_reg_elements_unlisted_unique$sequence)
  all_reg_elements_unlisted$Chopped <- list()
  all_reg_elements_unlisted$Chopped$Seq <- all_reg_elements_unlisted_unique$Chopped$Chopped_Seq[aaa]
  all_reg_elements_unlisted$Chopped$Score <- all_reg_elements_unlisted_unique$Chopped$Chopped_Score[aaa]
  all_reg_elements_unlisted$Chopped$GRanges <- all_reg_elements_unlisted_unique$Chopped$Chopped_GRanges[aaa]
  names(all_reg_elements_unlisted$Chopped$Seq) <- names(all_reg_elements_unlisted$sequence)
  names(all_reg_elements_unlisted$Chopped$Score) <- names(all_reg_elements_unlisted$sequence)
  names(all_reg_elements_unlisted$Chopped$GRanges) <- names(all_reg_elements_unlisted$sequence)
  for(i in 1:length(all_reg_elements_unlisted$Chopped$Seq)){
    names(all_reg_elements_unlisted$Chopped$Seq[[i]]) <- paste(names(all_reg_elements_unlisted$Chopped$Seq)[i], c(1:length(all_reg_elements_unlisted$Chopped$Seq[[i]])), sep = "_")
    rownames(all_reg_elements_unlisted$Chopped$Score[[i]]) <- paste(names(all_reg_elements_unlisted$Chopped$Seq)[i], c(1:length(all_reg_elements_unlisted$Chopped$Seq[[i]])), sep = "_")
  }
  # map the non-unique set to the gene based structure
  aa <- names(all_reg_elements_unlisted$Chopped$Score)
  aaa <- unlist(lapply(strsplit(aa, split = "\\."), "[[", 1))
  aaaa <- unique(aaa)
  all_reg_elements_perGene <- list()
  all_reg_elements_perGene$Chopped_Seq <- list()
  all_reg_elements_perGene$Chopped_Score <- list()
  all_reg_elements_perGene$Chopped_GRanges <- list()
  for(i in 1:length(aaaa)){
    aaw <- which(aaa %in% aaaa[i])
    aa1 <- do.call(rbind, all_reg_elements_unlisted$Chopped$Score[aaw])
    all_reg_elements_perGene$Chopped_Score[[i]] <- aa1
    aa2 <- do.call(c, all_reg_elements_unlisted$Chopped$Seq[aaw])
    all_reg_elements_perGene$Chopped_Seq[[i]] <- aa2
    aa3 <- lapply(all_reg_elements_unlisted$Chopped$GRanges[aaw], as.data.frame)
    aa3 <- do.call(rbind, aa3)
    aa3 <- makeGRangesFromDataFrame(aa3)
    all_reg_elements_perGene$Chopped_GRanges[[i]] <- aa3
  }
  
  names(all_reg_elements_perGene$Chopped_Seq) <- aaaa
  names(all_reg_elements_perGene$Chopped_Score) <- aaaa
  names(all_reg_elements_perGene$Chopped_GRanges) <- aaaa
  remove(aa)
  remove(aaa)
  remove(aaaa)
  remove(aa1)
  remove(aa2)
  remove(aa3)
  remove(aaw)
  
  # Filter for ChIP peak overlap,or overlap with any other GRanges object
  if(GRange_overlap){
    print("filtering for overlap with the specified GRange")
    all_reg_elements_perGene$interaction <- list()
    for(i in 1:length(all_reg_elements_perGene$Chopped_GRanges)){
      all_reg_elements_perGene$interaction[[i]] <-  OverlapInteractionBatchExtract(subjectCoord1=all_reg_elements_perGene$Chopped_GRanges[[i]],
                                                                                   subjectCoord2List=GRange_overlap_object,
                                                                                   mode = "both"
                                                                                   ,...minOverlapPercent = GRange_overlap_percent,
                                                                                   .int_piece_max_length = 10000)
    }

    all_reg_elements_perGene$Filtered <- list()
    all_reg_elements_perGene$Filtered$Chopped_Seq <- list()
    all_reg_elements_perGene$Filtered$Chopped_Score <- list()
    all_reg_elements_perGene$Filtered$Chopped_GRanges <- list()

    for(i in 1:length(all_reg_elements_perGene$Chopped_GRanges)){
      aaind <- which((all_reg_elements_perGene$interaction[[i]]$OverlapMat) > 0)
      all_reg_elements_perGene$Filtered$Chopped_GRanges[[i]] <- all_reg_elements_perGene$Chopped_GRanges[[i]][aaind]
      all_reg_elements_perGene$Filtered$Chopped_Seq[[i]] <- all_reg_elements_perGene$Chopped_Seq[[i]][aaind]
      if(length(aaind) == 1){
        print(paste("1 out of ", nrow(all_reg_elements_perGene$Chopped_Score[[i]]), "remained"))
        all_reg_elements_perGene$Filtered$Chopped_Score[[i]] <- t(as.matrix(all_reg_elements_perGene$Chopped_Score[[i]][aaind, ]))
      }else if(length(aaind) == 0){
        print(paste("0 left out of ", nrow(all_reg_elements_perGene$Chopped_Score[[i]])))
        all_reg_elements_perGene$Filtered$Chopped_Score[[i]] <- t(as.matrix(all_reg_elements_perGene$Chopped_Score[[i]][1, ]))
      }else{
        print(paste(length(aaind), " out of ", nrow(all_reg_elements_perGene$Chopped_Score[[i]]), "remained"))
        all_reg_elements_perGene$Filtered$Chopped_Score[[i]] <- all_reg_elements_perGene$Chopped_Score[[i]][aaind, ]
      }
    }
    names(all_reg_elements_perGene$Filtered$Chopped_Score) <- names(all_reg_elements_perGene$Chopped_Score)
    names(all_reg_elements_perGene$Filtered$Chopped_Seq) <- names(all_reg_elements_perGene$Chopped_Seq)
    names(all_reg_elements_perGene$Filtered$Chopped_GRanges) <-names(all_reg_elements_perGene$Chopped_GRanges)
  }
  return(all_reg_elements_perGene)
}

aaToRange <- function(x){
  if(nrow(x) > 0){
    aa <- matrix(nrow = nrow(x), ncol = 2)
    aa[, 1] <- x[, 1]
    aa[, 2] <- x[, 2] + 5
    return(aa)
  }else{
    return(numeric(0))
  }
}
#########################################################################################################
#########################################################################################################
#example
aa_example <- EnhancerChopper_wrapper(Entrez_IDs = rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)[3:5],
                                      .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                      .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                      .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                      .promoter.sequence=Promoter.Sequence.Char,
                                      #.Enhancer_seq=character(0),
                                      #.EnhancerGR=character(0),
                                      #.motifScore_output,
                                      .motifList=TF.motifs.Shrinked.t[3],
                                      .no_thresh = T,
                                      .LLR_thresh = 0,
                                      .piece_length = 1000,
                                      .step_size = 250,
                                      ..diff_max_LLR = F,
                                      dimer_filter_ER=T,
                                      dimer_list=ER.associated.reg.elements_uniq_pair_Concat_sort_un,
                                      dimer_filter_ER_LLR_th=2, 
                                      dimer_filter_corresponding_seq_nonProcessed=ER.associated.reg.elements$sequence,
                                      GRange_overlap=F,
                                      GRange_overlap_object=ReMapChIP.GRanges.list[52],
                                      GRange_overlap_percent=50)
#
aa_example_nod <- EnhancerChopper_wrapper(Entrez_IDs = rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)[3:5],
                                      .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                      .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                      .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                      .promoter.sequence=Promoter.Sequence.Char,
                                      #.Enhancer_seq=character(0),
                                      #.EnhancerGR=character(0),
                                      #.motifScore_output,
                                      .motifList=TF.motifs.Shrinked.t[3],
                                      .no_thresh = T,
                                      .LLR_thresh = 0,
                                      .piece_length = 1000,
                                      .step_size = 250,
                                      ..diff_max_LLR = F,
                                      dimer_filter_ER=F,
                                      dimer_list=ER.associated.reg.elements_uniq_pair_Concat_sort_un,
                                      dimer_filter_ER_LLR_th=2, 
                                      dimer_filter_corresponding_seq_nonProcessed=ER.associated.reg.elements$sequence,
                                      GRange_overlap=T,
                                      GRange_overlap_object=ReMapChIP.GRanges.list[52],
                                      GRange_overlap_percent=50)
#
aa_example_nod <- EnhancerChopper_wrapper(Entrez_IDs = rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)[1:2],
                                          .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                          .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                          .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                          .promoter.sequence=Promoter.Sequence.Char,
                                          #.Enhancer_seq=character(0),
                                          #.EnhancerGR=character(0),
                                          #.motifScore_output,
                                          .motifList=TF.motifs.Shrinked.t[3],
                                          .no_thresh = T,
                                          .LLR_thresh = 0,
                                          .piece_length = 1000,
                                          .step_size = 250,
                                          ..diff_max_LLR = F,
                                          dimer_filter_ER=F,
                                          dimer_list=ER.associated.reg.elements_uniq_pair_Concat_sort_un,
                                          dimer_filter_ER_LLR_th=2, 
                                          dimer_filter_corresponding_seq_nonProcessed=ER.associated.reg.elements$sequence,
                                          GRange_overlap=F,
                                          GRange_overlap_object=ReMapChIP.GRanges.list[52],
                                          GRange_overlap_percent=50)

#########################################################################################################
#########################################################################################################

my_dist_fuc <- function(mat){
  # This is a distance function to be used in heatmap hirarchical clustering, it assigns
  #  0 distance to entries that are equal and 1 to the entries that are not equal.
  aa_row_dis <- matrix(nrow = nrow(mat), ncol = nrow(mat))
  for(i in 1:nrow(mat)){
    for(j in i:nrow(mat)){
      aa_row_dis[i,j] <- sum(mat[i, ] != mat[j, ])
      aa_row_dis[j,i] <- aa_row_dis[i,j]
    }
  }
  aa_row_dis <- as.dist(aa_row_dis)
  return(aa_row_dis)
}
#########################################################################################################
#########################################################################################################
# example
#########################################################################################################
#########################################################################################################
GEMSTAT_weight_writer <- function(expmat,
                                  file_name,
                                  multi_enhancer_map=integer(0),
                                  fold_change_to_real=T, 
                                  to_kill=integer(0)){
  # expmat is the fold change expression matrix which contains values -1, 0, 1, or NA
  # This function creates a weight file to input GEMSTAT, it has the same format as the expression matrix
  # The weights are assigned such that all values of -1 , 0, and 1 have the same number*weight
  # fold_change_to_real : if the fold change gene expression matrix is given and the full version (control-treatment) is desired
  # to_kill : is an integer vector containing the labels that I want to assign zero weight to them.
  loss_weight <- matrix(nrow = nrow(expmat), ncol = ncol(expmat))
  rownames(loss_weight) <- rownames(expmat)
  mod_expmat <- expmat
  if(length(to_kill) > 0){
    mod_expmat[expmat  %in% to_kill] <- NA
  }
  
  balance_table <- table(mod_expmat)
  for(i in 1:length(balance_table)){
    loss_weight[mod_expmat %in% as.numeric(names(balance_table)[i])] <- max(balance_table)/balance_table[i]
  }
  loss_weight[is.na(loss_weight)] <- 0
  
  if(length(multi_enhancer_map) > 0){
    print("be carefull about indexing, the multi enhancer part uses 0 based index for genes")
    multi_expression_mat <- matrix(nrow = length(multi_enhancer_map), ncol = ncol(expmat))
    rownames(multi_expression_mat) <- character(length(multi_enhancer_map))
    for(i in 0:(nrow(expmat)-1)){
      cur_ind <- which(multi_enhancer_map %in% i)
      multi_expression_mat[cur_ind, ] <- matrix(rep(loss_weight[i+1,],
                                                    length(cur_ind)),
                                                byrow = T,
                                                nrow = length(cur_ind))
      rownames(multi_expression_mat)[cur_ind] <- paste(rownames(expmat)[i+1], c(1:length(cur_ind)), sep = "_")
    }
    if(fold_change_to_real){
      loss_weight <- matrix(nrow=nrow(multi_expression_mat), ncol = 2*ncol(multi_expression_mat))
      for(exp_col in 1:ncol(multi_expression_mat)){
        loss_weight[, (2*exp_col -1)] <- multi_expression_mat[, exp_col]
        loss_weight[, (2*exp_col)] <- multi_expression_mat[, exp_col]
      }
      rownames(loss_weight) <- rownames(multi_expression_mat)
    }else{
      loss_weight <- multi_expression_mat
    }

  }
  cat(c("Rows", as.character(c(1:ncol(loss_weight))), "\n"),
      file=paste0(file_name, ".weights"), sep = "\t")
  write.table(loss_weight, file=paste0(file_name, ".weights"),
              row.names=T, col.names=FALSE, quote=F, append=T, sep = "\t")
}


#########################################################################################################
#########################################################################################################
# example
aa_map_c <- unlist(lapply(aa_enhancer_list, length))
aa_map = integer(0)
for(i in 0:(length(aa_map_c)-1)){
  aa_map <- c(aa_map, rep(i, aa_map_c[i+1]))
}
aa_map <- c(0:38)
my_exp_mat <- my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3
GEMSTAT_weight_writer(expmat =my_exp_mat , 
                      "first_weights",
                      multi_enhancer_map = aa_map, to_kill = c(-1, 0))
#########################################################################################################
#########################################################################################################
# next function starts here

#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         Enhancers         ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
#EnhancerAtlas Enhancers
MCFEnhancers # is the charachter variable of all enhancers in MCF7, This is obtained from Enhacer Atlas database
MCFEnhancersdf <- as.data.frame(MCFEnhancers,stringsAsFactors = F)
MCFEnhancersGR <- makeGRangesFromDataFrame(MCFEnhancersdf)

#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################           Genes           ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
genesAll <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
tssdf.genes <- as.data.frame(resize(genesAll, width=1, fix='start'))
#make Grange and DataFrame for 
#1: promoter (1 kb upstream of tss)
Promoter1kbdf.gene <- tssdf.genes
Promoter1kbdf.gene$start <- Promoter1kbdf.gene$start - 1000
Promoter1kbdf.gene <- StartEndCoordinateCheck(inputRangedf = Promoter1kbdf.gene)
Promoter1kbGR.gene <- makeGRangesFromDataFrame(Promoter1kbdf.gene)

#2: Vicinity20kb (10 kb upstream and 10 kb downstream of tss)
Vicinity20kbdf.gene <- tssdf.genes
Vicinity20kbdf.gene$start <- Vicinity20kbdf.gene$start - 10000
Vicinity20kbdf.gene$end <- Vicinity20kbdf.gene$start + 10000
Vicinity20kbdf.gene <- StartEndCoordinateCheck(inputRangedf = Vicinity20kbdf.gene)
Vicinity20kbGR.gene <- makeGRangesFromDataFrame(Vicinity20kbdf.gene)

#3: Vicinity100kb (50 kb upstream and 50 kb downstream of tss)
Vicinity100kbdf.gene <- tssdf.genes
Vicinity100kbdf.gene$start <- Vicinity100kbdf.gene$start - 50000
Vicinity100kbdf.gene$end <- Vicinity100kbdf.gene$start + 50000
Vicinity100kbdf.gene <- StartEndCoordinateCheck(inputRangedf = Vicinity100kbdf.gene)
Vicinity100kbGR.gene <- makeGRangesFromDataFrame(Vicinity100kbdf.gene)

#4: Vicinity200kb (100 kb upstream and 100 kb downstream of tss)
Vicinity200kbdf.gene <- tssdf.genes
Vicinity200kbdf.gene$start <- Vicinity200kbdf.gene$start - 100000
Vicinity200kbdf.gene$end <- Vicinity200kbdf.gene$start + 100000
Vicinity200kbdf.gene <- StartEndCoordinateCheck(inputRangedf = Vicinity200kbdf.gene)
Vicinity200kbGR.gene <- makeGRangesFromDataFrame(Vicinity200kbdf.gene)
#########################################################################################################
#Get the enhancers for each gene:
#1.promoter (1kb upstream of tss)
#aaa <- select(txdb,keys = Promoter1kbdf.transcript$tx_name[Overlap.Promoter1kb.Transcript.Enhancer@from],columns = c("GENEID"),keytype = "TXNAME")$GENEID
#make a list: each entry corresponds enhancers associated to 1kb upstream of tss of one gene
Promoter1kb.Enhancer.By.gene <- overlapExtractor(subjectCoord1=Promoter1kbGR.gene, subjectCoord2=MCFEnhancersGR )
names(Promoter1kb.Enhancer.By.gene) <- Promoter1kbdf.gene$gene_id

#2 Vicinity20kb (10 kb upstream and 10 kb downstream of tss)
#make a list: each entry corresponds enhancers associated to 20 kb vicinity of tss of one gene
Vicinity20kb.Enhancer.By.gene <- overlapExtractor(subjectCoord1 = Vicinity20kbGR.gene,subjectCoord2 =MCFEnhancersGR )
names(Vicinity20kb.Enhancer.By.gene) <- Vicinity20kbdf.gene$gene_id

#2 Vicinity100kb (50 kb upstream and 50 kb downstream of tss)
#make a list: each entry corresponds enhancers associated to 20 kb vicinity of tss of one gene
Vicinity100kb.Enhancer.By.gene <- overlapExtractor(subjectCoord1 = Vicinity100kbGR.gene,subjectCoord2 =MCFEnhancersGR )
names(Vicinity100kb.Enhancer.By.gene) <- Vicinity100kbdf.gene$gene_id

#4: Vicinity200kb (100 kb upstream and 100 kb downstream of tss)
#make a list: each entry corresponds enhancers associated to 200 kb vicinity of tss of one gene
Vicinity200kb.Enhancer.By.gene <- overlapExtractor(subjectCoord1 = Vicinity200kbGR.gene,subjectCoord2 =MCFEnhancersGR )
names(Vicinity200kb.Enhancer.By.gene) <- Vicinity200kbdf.gene$gene_id

hist(unlist(lapply(Vicinity100kb.Enhancer.By.gene, length)), breaks = 50)

#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################        Transcripts        ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################


#########################################################################################################
#All genes , getting promoter regions, vicinity regions, map transcript name to gene name and vice versa
#Getting all Transcripts, All genes of hg19
TranscriptsAll <- transcripts(txdb,columns=c("tx_id", "tx_name"))
#make a dataframe of transcriptionstartsites
tssdf <- as.data.frame(resize(TranscriptsAll, width=1, fix='start'))
#map genes to transcripts
mapGenesToTranscripts <- select(txdb,keys = names(genesAll),columns = c("TXID","TXNAME"),keytype = "GENEID")
#map transcripts to genes
mapTranscriptsToGenes <- select(txdb,keys = tssdf$tx_name,columns = c("TXID","GENEID"),keytype = "TXNAME")
#make Grange and DataFrame for 
#1: promoter (1 kb upstream of tss)
Promoter1kbdf.transcript <- tssdf
Promoter1kbdf.transcript$start <- Promoter1kbdf.transcript$start - 1000
Promoter1kbdf.transcript <- StartEndCoordinateCheck(inputRangedf = Promoter1kbdf.transcript)
Promoter1kbGR.transcript <- makeGRangesFromDataFrame(Promoter1kbdf.transcript)

#2: Vicinity20kb (10 kb upstream and 10 kb downstream of tss)
Vicinity20kbdf.transcript <- tssdf
Vicinity20kbdf.transcript$start <- Vicinity20kbdf.transcript$start - 10000
Vicinity20kbdf.transcript$end <- Vicinity20kbdf.transcript$start + 10000
Vicinity20kbdf.transcript <- StartEndCoordinateCheck(inputRangedf = Vicinity20kbdf.transcript)
Vicinity20kbGR.transcript <- makeGRangesFromDataFrame(Vicinity20kbdf.transcript)

#3: Vicinity200kb (100 kb upstream and 100 kb downstream of tss)
Vicinity200kbdf.transcript <- tssdf
Vicinity200kbdf.transcript$start <- Vicinity200kbdf.transcript$start - 100000
Vicinity200kbdf.transcript$end <- Vicinity200kbdf.transcript$start + 100000
Vicinity200kbdf.transcript <- StartEndCoordinateCheck(inputRangedf = Vicinity200kbdf.transcript)
Vicinity200kbGR.transcript <- makeGRangesFromDataFrame(Vicinity200kbdf.transcript)
#########################################################################################################
#Get the enhancers for each transcript:
#1.promoter (1kb upstream of tss)
#make a list: each entry corresponds enhancers associated to 1kb upstream of tss of one transcript
Promoter1kb.Enhancer.By.Transcript <- overlapExtractor(subjectCoord1 = Promoter1kbGR.transcript,subjectCoord2 =MCFEnhancersGR )
names(Promoter1kb.Enhancer.By.Transcript) <- Promoter1kbdf.transcript$tx_name

#2 Vicinity20kb (10 kb upstream and 10 kb downstream of tss)
#make a list: each entry corresponds enhancers associated to 20 kb vicinity of tss of one transcript
Vicinity20kb.Enhancer.By.Transcript <- overlapExtractor(subjectCoord1 = Vicinity20kbGR.transcript,subjectCoord2 =MCFEnhancersGR )
names(Vicinity20kb.Enhancer.By.Transcript) <- Vicinity20kbdf.transcript$tx_name

#3: Vicinity200kb (100 kb upstream and 100 kb downstream of tss)
#make a list: each entry corresponds enhancers associated to 200 kb vicinity of tss of one transcript
Vicinity200kb.Enhancer.By.Transcript <- overlapExtractor(subjectCoord1 = Vicinity200kbGR.transcript,subjectCoord2 =MCFEnhancersGR )
names(Vicinity200kb.Enhancer.By.Transcript) <- Vicinity200kbdf.transcript$tx_name

#sum(unlist(lapply(Vicinity20kb.Enhancer.By.Transcript, length)) == 3 )
#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         Promoters         ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
Promoter1kbdf.transcript#dataframe of promoters: 1kb upstream of TSS of each transcript
Promoter1kbGR.transcript#Granges object of the above dataframe
Promoter1kbdf.gene#dataframe of promoters: 1kb upstream of TSS of each gene
Promoter1kbGR.gene#Granges object of the above dataframe
tss.gene.GR <- makeGRangesFromDataFrame(tssdf.genes)
Promoter.Sequence.DNASTRING <- getPromoterSeq(tss.gene.GR, BSgenome.Hsapiens.UCSC.hg19, upstream=1000, downstream = 0)
Promoter.Sequence.Char <- as.character(Promoter.Sequence.DNASTRING)
names(Promoter.Sequence.Char) <- names(tss.gene.GR)
#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################        interaction        ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#Genome4D
Genome4DMCF7 # is the data frame with all interactions in MCF7 cells
Genome4DMCF7Exp <- Genome4DMCF7[Genome4DMCF7$Detection_Method != "IM-PET",]#JUST experimental interactions:excluding IM-PET predicted interactions
###Creating the Enhancer-Enhancer interaction List
EnhancerEnhancerIntList <- InteractionExtractor(subjectCoord1 =MCFEnhancersdf ,subjectCoord2 =MCFEnhancersdf )
###Creating promoter-Enhancer interaction List --> promoter is per transcript (1 kb upstream of TSS)
Promoter.transcrpt.EnhancerIntList <- InteractionExtractor(subjectCoord1 =Promoter1kbdf.transcript ,subjectCoord2 =MCFEnhancersdf )
hist(unlist(lapply(Promoter.transcrpt.EnhancerIntList,length)),breaks = 200)
###Creating promoter-Enhancer interaction List --> promoter is per gene (1 kb upstream of TSS)
Promoter.gene.EnhancerIntList <- InteractionExtractor(subjectCoord1 =Promoter1kbdf.gene ,subjectCoord2 =MCFEnhancersdf )

###Creating promoter-promoter interaction List --> promoter is per transcript (1 kb upstream of TSS)
Promoter.transcrpt.PromoterIntList <- InteractionExtractor(subjectCoord1 =Promoter1kbdf.transcript ,subjectCoord2 =Promoter1kbdf.transcript )

###Creating promoter-promoter interaction List --> promoter is per gene (1 kb upstream of TSS)
Promoter.gene.PromoterIntList <- InteractionExtractor(subjectCoord1 =Promoter1kbdf.gene ,subjectCoord2 =Promoter1kbdf.gene )

###In order to construct bag of enhancers for each gene:
#1.Get all enhancers in certain genomic range of gene's TSS (20kb or 200kb) and the enhancers which interact with them
#20kb
Vicinity20kb.Enhancer.By.gene
Vicinity20kb.Enhancer.By.gene.Extended <-  InteractionExtender(InteractionList = Vicinity20kb.Enhancer.By.gene,Subject = MCFEnhancersGR)

#100kb
Vicinity100kb.Enhancer.By.gene
Vicinity100kb.Enhancer.By.gene.Extended <-  InteractionExtender(InteractionList = Vicinity100kb.Enhancer.By.gene,Subject = MCFEnhancersGR)
#plots for 100 kb enhancers per gene till this step:
a = unlist(lapply(Vicinity100kb.Enhancer.By.gene.Extended[[1]], length))
aa = unlist(lapply(Vicinity100kb.Enhancer.By.gene, length))
aaa = a - aa
hist(aaa, breaks = 50, xlab = "No. of enhs added.Per.Gene based on enh-enh interaction", main = "")
hist(aaa[aaa > 0], breaks = 50, xlab = "No.of enhs added.Per.Gene based on enh-enh interaction(>0)", main = "")
hist(unlist(lapply(Vicinity100kb.Enhancer.By.gene, length)), breaks=50, xlab="Number of enhancers in 100kb of gene (No interaction)", main="")
hist(unlist(lapply(Vicinity100kb.Enhancer.By.gene.Extended[[1]], length)), breaks=50, xlab="Number of enhancers in 100kb of gene (Plus interaction)", main="")

#200kb
Vicinity200kb.Enhancer.By.gene
Vicinity200kb.Enhancer.By.gene.Extended <-  InteractionExtender(InteractionList = Vicinity200kb.Enhancer.By.gene,Subject = MCFEnhancersGR)

#2.Get all enhancers that interact with the promoter of the gene
Promoter.gene.EnhancerIntList
#3.merge the last two class
#20kb
Vicinity20kb.Enhancer.plusPromoter.By.gene.Extended <-   mapply(c, Promoter.gene.EnhancerIntList, Vicinity20kb.Enhancer.By.gene.Extended[[1]], SIMPLIFY=FALSE)
Vicinity20kb.Enhancer.plusPromoter.By.gene.Extended = mapply(unique,Vicinity20kb.Enhancer.plusPromoter.By.gene.Extended )

#100kb
Vicinity100kb.Enhancer.plusPromoter.By.gene.Extended <-   mapply(c, Promoter.gene.EnhancerIntList, Vicinity100kb.Enhancer.By.gene.Extended[[1]], SIMPLIFY=FALSE)
Vicinity100kb.Enhancer.plusPromoter.By.gene.Extended = mapply(unique,Vicinity100kb.Enhancer.plusPromoter.By.gene.Extended )
Vicinity100kb.Enhancer.plusPromoter.By.gene = mapply(c, Promoter.gene.EnhancerIntList, Vicinity100kb.Enhancer.By.gene, SIMPLIFY=FALSE)
Vicinity100kb.Enhancer.plusPromoter.By.gene = mapply(unique,Vicinity100kb.Enhancer.plusPromoter.By.gene )
names(Vicinity100kb.Enhancer.plusPromoter.By.gene) <- names(genesAll)
hist(unlist(lapply(Vicinity100kb.Enhancer.plusPromoter.By.gene, length)), breaks=50, xlab="Number of enhancers in 100kb of gene (plus prom-enh Int)", main="")
a = unlist(lapply(Vicinity100kb.Enhancer.By.gene.Extended[[1]], length))
aa = unlist(lapply(Vicinity100kb.Enhancer.plusPromoter.By.gene.Extended, length))
aaa = aa - a
hist(aaa, breaks = 50, xlab = "No. of enhs added.Per.Gene based on prom-enh interaction", main = "")
hist(aaa[aaa > 0], breaks = 50, xlab = "No.of enhs added.Per.Gene based on prom-enh interaction(>0)", main = "")

#200kb
Vicinity200kb.Enhancer.plusPromoter.By.gene.Extended <-   mapply(c, Promoter.gene.EnhancerIntList, Vicinity200kb.Enhancer.By.gene.Extended[[1]], SIMPLIFY=FALSE)
Vicinity200kb.Enhancer.plusPromoter.By.gene.Extended = mapply(unique,Vicinity200kb.Enhancer.plusPromoter.By.gene.Extended )

#4.add the promoter of the gene to get the final bag, add the promoters interacting with the genes promoter to the bag
Promoter1kbGR.gene
Promoter.gene.PromoterIntList
hist(unlist(lapply(Promoter.gene.PromoterIntList, length)), breaks=50, xlab="Number of promoters interacting per promoter", main="")
a = unlist(lapply(Promoter.gene.PromoterIntList, length))
hist(a[a > 0], breaks=50, xlab="Number of promoters interacting per promoter(>0)", main="")
#5. Get the size of the bag of enhancers_Plus_promoters for each gene #UPDATE: NOT USING EXTENDED VERSION OF ENHANCERS IN THE BAG:
#USE JUST THE ENHANCERS WITHIN +/-50 KB AND ENHANCERS WHICH INTERACT WITH THE PROMOTER OF THE GENE
aa = unlist(lapply(Vicinity100kb.Enhancer.plusPromoter.By.gene, length))
a = unlist(lapply(Promoter.gene.PromoterIntList, length))
TotalBagSize.PerGene = a + aa + 1
hist(TotalBagSize.PerGene, xlab = "size of the bag per gene", breaks = 52, main ="")
hist(TotalBagSize.PerGene[TotalBagSize.PerGene > 1], xlab = "size of the bag per gene (>1)", breaks = 50, main ="")

###In order to construct bag of enhancers for each transcript:
#1.Get all enhancers in certain genomic range of transcript's TSS (20kb or 200kb) 
#20kb
Vicinity20kb.Enhancer.By.Transcript
Vicinity20kb.Enhancer.By.Transcript.Extended <-  InteractionExtender(InteractionList = Vicinity20kb.Enhancer.By.Transcript,Subject = MCFEnhancersGR)
#200kb
Vicinity200kb.Enhancer.By.Transcript
Vicinity200kb.Enhancer.By.Transcript.Extended <-  InteractionExtender(InteractionList = Vicinity200kb.Enhancer.By.Transcript,Subject = MCFEnhancersGR)

#2.Get all enhancers that interact with the promoter of the transcript
Promoter.transcrpt.EnhancerIntList

#3.merge the last two class
#20kb
Vicinity20kb.Enhancer.plusPromoter.By.Transcript.Extended <-   mapply(c, Promoter.transcrpt.EnhancerIntList, Vicinity20kb.Enhancer.By.Transcript.Extended[[1]], SIMPLIFY=FALSE)
#200kb
Vicinity200kb.Enhancer.plusPromoter.By.Transcript.Extended <-   mapply(c, Promoter.transcrpt.EnhancerIntList, Vicinity200kb.Enhancer.By.Transcript.Extended[[1]], SIMPLIFY=FALSE)

#4.add the promoter of the transcript to get the final bag, add the promoters interacting with the transcript's promoter to the bag
Promoter1kbGR.transcript
Promoter.transcrpt.PromoterIntList
#########################################################################################################

#########################################################################################################
#3Cdb
ThreeCdb <- read.delim(file = "ConformationData/3CDB_MCF7.txt",header = T,stringsAsFactors = F)
#########################################################################################################


#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################          Chip             ########################################
######################################          part1            ########################################
#########################################################################################################
#########################################################################################################
#looking at overlap AND interaction of chip peaks of all datasets with enhncers and genes
ERpickCistromeTimePoint#is the list contining all the chip datasets from different time points after E2 treatments

#Get the overlap of all enhancers with ER chip peaks
Enhancer.ERchip.Overlap.byEnhancer <- list()#a list where each entry is the overlap list of the enhancers with one chip dataset
for(i in 1:length(ERpickCistromeTimePoint)){
  Enhancer.ERchip.Overlap.byEnhancer[[i]] <- overlapExtractor(subjectCoord1 = MCFEnhancersGR,subjectCoord2 =ERpickCistromeTimePoint[[i]] )
}
names(Enhancer.ERchip.Overlap.byEnhancer) <- names(ERpickCistromeTimePoint)
#Create a binary matrix where each row represents an enhancer, each column represents  a chip data set: 1 means that enhancer overlaps with at least one chip peak in that dataset
Enhancer.ERchip.Overlap.byEnhancer.Binary <- matrix(0L,nrow = nrow(MCFEnhancersdf),ncol = length(Enhancer.ERchip.Overlap.byEnhancer))
for(i in 1:ncol(Enhancer.ERchip.Overlap.byEnhancer.Binary)){
  Enhancer.ERchip.Overlap.byEnhancer.Binary[(unlist(lapply(Enhancer.ERchip.Overlap.byEnhancer[[i]],length)) > 0),i] <- 1
}
#plot the histogram of overlaps between Enhancers and chip peaks
hist(rowSums(Enhancer.ERchip.Overlap.byEnhancer.Binary[(rowSums(Enhancer.ERchip.Overlap.byEnhancer.Binary) != 0),]),breaks = 200,main = "No. of overlapping Chip peaks for enhancers")
#Get the overlap of all gene promoters with ER chip peaks
Promoter.gene.ERchip.Overlap.byPromoter <- list()#a list where each entry is the overlap list of the promoters with one chip dataset
for(i in 1:length(ERpickCistromeTimePoint)){
  Promoter.gene.ERchip.Overlap.byPromoter[[i]] <- overlapExtractor(subjectCoord1 = Promoter1kbGR.gene,subjectCoord2 =ERpickCistromeTimePoint[[i]] )
}
names(Promoter.gene.ERchip.Overlap.byPromoter) <- names(ERpickCistromeTimePoint)
#Create a binary matrix where each row represents a gene promoter, each column represents  a chip data set: 1 means that enhancer overlaps with at least one chip peak in that dataset
Promoter.gene.ERchip.Overlap.byPromoter.Binary <- matrix(0L,nrow = nrow(Promoter1kbdf.gene),ncol = length(Promoter.gene.ERchip.Overlap.byPromoter))
for(i in 1:ncol(Promoter.gene.ERchip.Overlap.byPromoter.Binary)){
  Promoter.gene.ERchip.Overlap.byPromoter.Binary[(unlist(lapply(Promoter.gene.ERchip.Overlap.byPromoter[[i]],length)) > 0),i] <- 1
}
#plot the histogram of overlaps between gene promoters and chip peaks
hist(rowSums(Promoter.gene.ERchip.Overlap.byPromoter.Binary[(rowSums(Promoter.gene.ERchip.Overlap.byPromoter.Binary) != 0),]),breaks = 200,main = "No. of overlapping Chip peaks for gene promoters")

#Get the overlap of all transcript promoters with ER chip peaks
Promoter.transcript.ERchip.Overlap.byPromoter <- list()#a list where each entry is the overlap list of the enhancers with one chip dataset
for(i in 1:length(ERpickCistromeTimePoint)){
  Promoter.transcript.ERchip.Overlap.byPromoter[[i]] <- overlapExtractor(subjectCoord1 = Promoter1kbGR.transcript,subjectCoord2 =ERpickCistromeTimePoint[[i]] )
}
names(Promoter.transcript.ERchip.Overlap.byPromoter) <- names(ERpickCistromeTimePoint)
#Create a binary matrix where each row represents a transcript promoter, each column represents  a chip data set: 1 means that enhancer overlaps with at least one chip peak in that dataset
Promoter.transcript.ERchip.Overlap.byPromoter.Binary <- matrix(0L,nrow = nrow(Promoter1kbdf.transcript),ncol = length(Promoter.transcript.ERchip.Overlap.byPromoter))
for(i in 1:ncol(Promoter.transcript.ERchip.Overlap.byPromoter.Binary)){
  Promoter.transcript.ERchip.Overlap.byPromoter.Binary[(unlist(lapply(Promoter.transcript.ERchip.Overlap.byPromoter[[i]],length)) > 0),i] <- 1
}

#Get the interacting enhancers with a chip peak
Enhancer.ERchip.Interaction.byEnhancer <- list()
for(i in 1:length(ERpickCistromeTimePoint)){
  Enhancer.ERchip.Interaction.byEnhancer[[i]] <-   InteractionExtractor(subjectCoord1 = MCFEnhancersGR,subjectCoord2 =ERpickCistromeTimePoint[[i]] )
}
##make the binary matrix just like the overlap binary matrix but this time a 1 indicates an interaction between enhancer and at least one chip peak in that dataset
Enhancer.ERchip.Interaction.byEnhancer.Binary <- matrix(0L,nrow = nrow(MCFEnhancersdf),ncol = length(Enhancer.ERchip.Interaction.byEnhancer))
for(i in 1:ncol(Enhancer.ERchip.Interaction.byEnhancer.Binary)){
  Enhancer.ERchip.Interaction.byEnhancer.Binary[(unlist(lapply(Enhancer.ERchip.Interaction.byEnhancer[[i]],length)) > 0),i] <- 1
}
##Create a binary matrix where each row represents an enhancer, each column represents  a chip data set: 1 means that enhancer either overlaps or interacts with at least one chip peak in that dataset
Enhancer.ERchip.interORover.byEnhancer.Binary  <- Enhancer.ERchip.Interaction.byEnhancer.Binary + Enhancer.ERchip.Overlap.byEnhancer.Binary
Enhancer.ERchip.interORover.byEnhancer.Binary[Enhancer.ERchip.interORover.byEnhancer.Binary > 0] <- 1

#Get the interacting promoters with a chip peak
##gene
Promoter.gene.ERchip.interaction.byPromoter <- list()#a list where each entry is the interaction list of the promoters with one chip dataset
for(i in 1:length(ERpickCistromeTimePoint)){
  Promoter.gene.ERchip.interaction.byPromoter[[i]] <- InteractionExtractor(subjectCoord1 = Promoter1kbGR.gene,subjectCoord2 =ERpickCistromeTimePoint[[i]] )
}
names(Promoter.gene.ERchip.interaction.byPromoter) <- names(ERpickCistromeTimePoint)
###Create a binary matrix where each row represents a gene promoter, each column represents  a chip data set: 1 means that enhancer interacts with at least one chip peak in that dataset
Promoter.gene.ERchip.interaction.byPromoter.Binary <- matrix(0L,nrow = nrow(Promoter1kbdf.gene),ncol = length(Promoter.gene.ERchip.interaction.byPromoter))
for(i in 1:ncol(Promoter.gene.ERchip.interaction.byPromoter.Binary)){
  Promoter.gene.ERchip.interaction.byPromoter.Binary[(unlist(lapply(Promoter.gene.ERchip.interaction.byPromoter[[i]],length)) > 0),i] <- 1
}
###Create a binary matrix where each row represents a gene promoter, each column represents  a chip data set: 1 means that gene promoter either overlaps or interacts with at least one chip peak in that dataset
Promoter.gene.ERchip.interORover.byPromoter.Binary  <- Promoter.gene.ERchip.interaction.byPromoter.Binary + Promoter.gene.ERchip.Overlap.byPromoter.Binary
Promoter.gene.ERchip.interORover.byPromoter.Binary[Promoter.gene.ERchip.interORover.byPromoter.Binary > 0] <- 1

##transcript
Promoter.transcript.ERchip.interaction.byPromoter <- list()#a list where each entry is the interaction list of the enhancers with one chip dataset
for(i in 1:length(ERpickCistromeTimePoint)){
  Promoter.transcript.ERchip.interaction.byPromoter[[i]] <- InteractionExtractor(subjectCoord1 = Promoter1kbGR.transcript,subjectCoord2 =ERpickCistromeTimePoint[[i]] )
}
names(Promoter.transcript.ERchip.interaction.byPromoter) <- names(ERpickCistromeTimePoint)
#Create a binary matrix where each row represents a transcript promoter, each column represents  a chip data set: 1 means that enhancer interacts with at least one chip peak in that dataset
Promoter.transcript.ERchip.interaction.byPromoter.Binary <- matrix(0L,nrow = nrow(Promoter1kbdf.transcript),ncol = length(Promoter.transcript.ERchip.interaction.byPromoter))
for(i in 1:ncol(Promoter.transcript.ERchip.interaction.byPromoter.Binary)){
  Promoter.transcript.ERchip.interaction.byPromoter.Binary[(unlist(lapply(Promoter.transcript.ERchip.interaction.byPromoter[[i]],length)) > 0),i] <- 1
}
#Create a binary matrix where each row represents a transcript promoter, each column represents  a chip data set: 1 means that transcript promoter either overlaps or interacts with at least one chip peak in that dataset
Promoter.transcript.ERchip.interORover.byPromoter.Binary  <- Promoter.transcript.ERchip.interaction.byPromoter.Binary + Promoter.transcript.ERchip.Overlap.byPromoter.Binary
Promoter.transcript.ERchip.interORover.byPromoter.Binary[Promoter.transcript.ERchip.interORover.byPromoter.Binary > 0] <- 1

#Get the overlap of 20kb vicinity of each gene with chip peaks
Vicinity20kb.gene.ERchip.Overlap.byGene <- list()#a list where each entry is the overlap list of the vicinity 20kb with one chip dataset
for(i in 1:length(ERpickCistromeTimePoint)){
  Vicinity20kb.gene.ERchip.Overlap.byGene[[i]] <- overlapExtractor(subjectCoord1 = Vicinity20kbGR.gene,subjectCoord2 =ERpickCistromeTimePoint[[i]] )
}
names(Vicinity20kb.gene.ERchip.Overlap.byGene) <- names(ERpickCistromeTimePoint)
##Create a binary matrix where each row represents a gene 20 kb vicinity, each column represents  a chip data set: 1 means that gene 20 kb vicinity overlaps with at least one chip peak in that dataset
Vicinity20kb.gene.ERchip.Overlap.byGene.Binary <- matrix(0L,nrow = nrow(Vicinity20kbdf.gene),ncol = length(Vicinity20kb.gene.ERchip.Overlap.byGene))
for(i in 1:ncol(Vicinity20kb.gene.ERchip.Overlap.byGene.Binary)){
  Vicinity20kb.gene.ERchip.Overlap.byGene.Binary[(unlist(lapply(Vicinity20kb.gene.ERchip.Overlap.byGene[[i]],length)) > 0),i] <- 1
}
#Get the overlap of 200kb vicinity of each gene with chip peaks
Vicinity200kb.gene.ERchip.Overlap.byGene <- list()#a list where each entry is the overlap list of the vicinity 200kb with one chip dataset
for(i in 1:length(ERpickCistromeTimePoint)){
  Vicinity200kb.gene.ERchip.Overlap.byGene[[i]] <- overlapExtractor(subjectCoord1 = Vicinity200kbGR.gene,subjectCoord2 =ERpickCistromeTimePoint[[i]] )
}
names(Vicinity200kb.gene.ERchip.Overlap.byGene) <- names(ERpickCistromeTimePoint)
##Create a binary matrix where each row represents a gene 200 kb vicinity, each column represents  a chip data set: 1 means that  gene 200 kb vicinity overlaps with at least one chip peak in that dataset
Vicinity200kb.gene.ERchip.Overlap.byGene.Binary <- matrix(0L,nrow = nrow(Vicinity200kbdf.gene),ncol = length(Vicinity200kb.gene.ERchip.Overlap.byGene))
for(i in 1:ncol(Vicinity200kb.gene.ERchip.Overlap.byGene.Binary)){
  Vicinity200kb.gene.ERchip.Overlap.byGene.Binary[(unlist(lapply(Vicinity200kb.gene.ERchip.Overlap.byGene[[i]],length)) > 0),i] <- 1
}
#Get the overlap of 20kb vicinity of each transcript with chip peaks
Vicinity20kb.transcript.ERchip.Overlap.byTranscript <- list()#a list where each entry is the overlap list of the vicinity 20kb with one chip dataset
for(i in 1:length(ERpickCistromeTimePoint)){
  Vicinity20kb.transcript.ERchip.Overlap.byTranscript[[i]] <- overlapExtractor(subjectCoord1 = Vicinity20kbGR.transcript,subjectCoord2 =ERpickCistromeTimePoint[[i]] )
}
names(Vicinity20kb.transcript.ERchip.Overlap.byTranscript) <- names(ERpickCistromeTimePoint)
#Create a binary matrix where each row represents a transcript vicinity 20kb, each column represents  a chip data set: 1 means that vicinity 20kb overlaps with at least one chip peak in that dataset
Vicinity20kb.transcript.ERchip.Overlap.byTranscript.Binary <- matrix(0L,nrow = length(Vicinity20kbGR.transcript),ncol = length(Vicinity20kb.transcript.ERchip.Overlap.byTranscript))
for(i in 1:ncol(Vicinity20kb.transcript.ERchip.Overlap.byTranscript.Binary)){
  Vicinity20kb.transcript.ERchip.Overlap.byTranscript.Binary[(unlist(lapply(Vicinity20kb.transcript.ERchip.Overlap.byTranscript[[i]],length)) > 0),i] <- 1
}
#Get the overlap of 200kb vicinity of each transcript with chip peaks
Vicinity200kb.transcript.ERchip.Overlap.byTranscript <- list()#a list where each entry is the overlap list of the vicinity 200kb with one chip dataset
for(i in 1:length(ERpickCistromeTimePoint)){
  Vicinity200kb.transcript.ERchip.Overlap.byTranscript[[i]] <- overlapExtractor(subjectCoord1 = Vicinity200kbGR.transcript,subjectCoord2 =ERpickCistromeTimePoint[[i]] )
}
names(Vicinity200kb.transcript.ERchip.Overlap.byTranscript) <- names(ERpickCistromeTimePoint)
#Create a binary matrix where each row represents a transcript vicinity 200kb, each column represents  a chip data set: 1 means that vicinity 200kb overlaps with at least one chip peak in that dataset
Vicinity200kb.transcript.ERchip.Overlap.byTranscript.Binary <- matrix(0L,nrow = length(Vicinity200kbGR.transcript),ncol = length(Vicinity200kb.transcript.ERchip.Overlap.byTranscript))
for(i in 1:ncol(Vicinity200kb.transcript.ERchip.Overlap.byTranscript.Binary)){
  Vicinity200kb.transcript.ERchip.Overlap.byTranscript.Binary[(unlist(lapply(Vicinity200kb.transcript.ERchip.Overlap.byTranscript[[i]],length)) > 0),i] <- 1
}
#########################################################################################################
#########################################################################################################
#What genes to consider for modeling? THe ones associated with an ER chip peak in some replicates of each condition
#what are the associations we care about:
##overlap or interact with an enhancer associated with the gene
##overlap or interact with a promoter associated with the gene (gene's promoter and its interacting promoters)
###gene
AssociatedChip.gene.20kb <-  ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList , gene.enhancer = Vicinity20kb.Enhancer.plusPromoter.By.gene.Extended ,enhancerChipBinary = Enhancer.ERchip.interORover.byEnhancer.Binary,promoterChipBinary = Promoter.gene.ERchip.interORover.byPromoter.Binary)
AssociatedChip.gene.200kb <- ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList , gene.enhancer = Vicinity200kb.Enhancer.plusPromoter.By.gene.Extended ,enhancerChipBinary = Enhancer.ERchip.interORover.byEnhancer.Binary,promoterChipBinary = Promoter.gene.ERchip.interORover.byPromoter.Binary)
###transcript
AssociatedChip.transcript.20kb  <-  ChipPeakInvestigator(promoter.promoter = Promoter.transcrpt.PromoterIntList , gene.enhancer = Vicinity20kb.Enhancer.plusPromoter.By.Transcript.Extended ,enhancerChipBinary = Enhancer.ERchip.interORover.byEnhancer.Binary,promoterChipBinary = Promoter.transcript.ERchip.interORover.byPromoter.Binary)
AssociatedChip.transcript.200kb <-  ChipPeakInvestigator(promoter.promoter = Promoter.transcrpt.PromoterIntList , gene.enhancer = Vicinity200kb.Enhancer.plusPromoter.By.Transcript.Extended ,enhancerChipBinary = Enhancer.ERchip.interORover.byEnhancer.Binary,promoterChipBinary = Promoter.transcript.ERchip.interORover.byPromoter.Binary)

#Filter for genes or transcripts that their association with a chip peak has been repeated in experiments of the same conditions
names(ERpickCistromeTimePoint)
 sum(rowSums(AssociatedChip.gene.20kb[[2]][,c(c(1:9),16)]) >7 )
 sum(rowSums(AssociatedChip.gene.20kb[[2]][,c(10,11,17,18)]) > 3 )
 sum(rowSums(AssociatedChip.gene.20kb[[2]][,12:15]) > 3 )
 length(unique(which(rowSums(AssociatedChip.gene.20kb[[2]][,c(c(2:6),c(8,9),16)]) >5),which(rowSums(AssociatedChip.gene.20kb[[2]][,c(11,17,18)]) > 1 ),which(rowSums(AssociatedChip.gene.20kb[[2]][,c(12,15)]) > 0 )))
# sum(rowSums(GeneERchipVicinityPlusInteractionBinary[,c(10,11,17,18)]) > 1 )
# sum(rowSums(GeneERchipVicinityPlusInteractionBinary[,12:15]) > 2 )
-c(1,7,13,14)
 #########################################################################################################
 #########################################################################################################
 ######################################                           ########################################
 ######################################          Chip             ########################################
 ######################################         part II           ########################################
 ######################################          ReMap            ########################################
 ######################################                           ########################################
 #########################################################################################################
 #########################################################################################################
 #Chip data from http://tagc.univ-mrs.fr/remap/celltype.php?CT=MCF7&page=ct
ReMapChIP <- readPeakFile("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Chip-Data/remap2018_MCF7_all_macs2_hg38_v1_2.bed",as = "data.frame")
#Build a dataframe for the name of Chip Experiments containing 1:study number, 2: ChIPed Factor 3: treatment
hist(ReMapChIP$V3 - ReMapChIP$V2, xlim = c(0,2000), breaks = 300)
summary(ReMapChIP$V3 - ReMapChIP$V2)
ReMapChIP_expName = data.frame(studyNumber = character(length = length(levels(x = ReMapChIP$V4))),
                               FactorName=character(length = length(levels(x = ReMapChIP$V4))),
                               Treatment=character(length = length(levels(x = ReMapChIP$V4))),
                               stringsAsFactors = F)

for (i in 1:length(levels(x = ReMapChIP$V4))){
  ReMapChIP_expName[i,] <- strsplit(as.character(levels(x = ReMapChIP$V4)[i]),'\\.')[[1]]
}

#Name of the factors that we have chip for:
ReMapChIP_FactorName_unique <- unique(ReMapChIP_expName$FactorName)

#Group the ReMapChIP into 85 different lists. each list contains the Chip Peaks for one factor
ReMapChIP.df.list <- list()
for(i in 1:length(ReMapChIP_FactorName_unique)){
  aa <- which(as.numeric(ReMapChIP$V4) %in%  which(ReMapChIP_expName$FactorName == ReMapChIP_FactorName_unique[i]))
  ReMapChIP.df.list[[i]] <- ReMapChIP[aa,]
}
for(i in 1:length(ReMapChIP.df.list)){
  names(ReMapChIP.df.list[[i]]) <- c("Chr", "start", "end", "experiment","score","strand", "summit1","summit2","studyNumber" )
}
names(ReMapChIP.df.list) <- ReMapChIP_FactorName_unique
#Create GRanges objects for each of the 85 entires of the list and store in another list:
ReMapChIP.GRanges.list <- mapply(FUN = makeGRangesFromDataFrame, keep.extra.columns = T,ReMapChIP.df.list)
names(ReMapChIP.GRanges.list) <- ReMapChIP_FactorName_unique
ReMapChIP.GRanges.list_merged <- list()
for(i in 1:length(ReMapChIP.GRanges.list)){ # continue from 12
  print(i)
  print(ReMapChIP_FactorName_unique[i])
  aa <- findOverlaps(ReMapChIP.GRanges.list[[i]])
  ReMapChIP.GRanges.list_merged[[i]] <- mergeConnectedRanges(ReMapChIP.GRanges.list[[i]], aa)
}
names(ReMapChIP.GRanges.list_merged) <- ReMapChIP_FactorName_unique
save(list = c("mergeConnectedRanges", "ReMapChIP.GRanges.list", "extractClustersFromSelfHits"),
     file = "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/REMAP_chip.RData")
#########################################################################################################
#########################################################################################################
#########################################################################################################
#find [overlap and interaction] of [Enhancers and promoters] with the grouped ChIP peaks

#########################################################################################################
#############################################Enhancer REMAP ChIP Overlap################################

#Get the "overlap" of all "enhancers" with REMAP chip peaks
Enhancer.ReMapchip.Overlap.byEnhancer <- list()#a list where each entry is the overlap list of the enhancers with one chip dataset
for(i in 1:length(ReMapChIP.GRanges.list)){
  Enhancer.ReMapchip.Overlap.byEnhancer[[i]] <- overlapExtractor(subjectCoord1 = MCFEnhancersGR,subjectCoord2 =ReMapChIP.GRanges.list[[i]] )
}
names(Enhancer.ReMapchip.Overlap.byEnhancer) <- names(ReMapChIP.GRanges.list)
#Create a  matrix where each row represents an enhancer, each column represents  a ChIPed factor: entry[i,j] is the number of Chip-peaks of factor j that overlap with the enhancer i
Enhancer.ReMapchip.Overlap.byEnhancer.Mat <- matrix(0L,nrow = nrow(MCFEnhancersdf),ncol = length(Enhancer.ReMapchip.Overlap.byEnhancer))
for(i in 1:ncol(Enhancer.ReMapchip.Overlap.byEnhancer.Mat)){
  Enhancer.ReMapchip.Overlap.byEnhancer.Mat[,i] <- unlist(lapply(Enhancer.ReMapchip.Overlap.byEnhancer[[i]],length))
}
colnames(Enhancer.ReMapchip.Overlap.byEnhancer.Mat) <- names(ReMapChIP.GRanges.list)
boxplot.matrix(Enhancer.ReMapchip.Overlap.byEnhancer.Mat,las = 2)
#Create a binary matrix where each row represents an enhancer, each column represents  a ChIPed factor: 1 means that enhancer overlaps with at least one chip peak of that factor
Enhancer.ReMapchip.Overlap.byEnhancer.Mat.Binary <- matrix(0L,nrow = nrow(MCFEnhancersdf),ncol = length(Enhancer.ReMapchip.Overlap.byEnhancer))
Enhancer.ReMapchip.Overlap.byEnhancer.Mat.Binary[Enhancer.ReMapchip.Overlap.byEnhancer.Mat > 0] <- 1
colnames(Enhancer.ReMapchip.Overlap.byEnhancer.Mat.Binary) <- names(ReMapChIP.GRanges.list)

#plot the histogram of overlaps between Enhancers and chip peaks
hist(rowSums(Enhancer.ReMapchip.Overlap.byEnhancer.Mat.Binary[(rowSums(Enhancer.ReMapchip.Overlap.byEnhancer.Mat.Binary) != 0),]),breaks = 200,main = "No. of overlapping Chip peaks for enhancers(>0)",xlab = "")
hist(rowSums(Enhancer.ReMapchip.Overlap.byEnhancer.Mat.Binary),breaks = 200,xlab = "No. of ChIPed factors overlapping with each enhancers",main = "")

#########################################################################################################
#############################################Promoter REMAP ChIP Overlap################################
#Get the "overlap" of all gene "promoters" with REMAP chip peaks

Promoter.gene.ReMapchip.Overlap.byPromoter <- list()#a list where each entry is the overlap list of the All promoters with All peaks of one tf
for(i in 1:length(ReMapChIP.GRanges.list)){
  Promoter.gene.ReMapchip.Overlap.byPromoter[[i]] <- overlapExtractor(subjectCoord1 = Promoter1kbGR.gene,subjectCoord2 =ReMapChIP.GRanges.list[[i]] )
  print(i)
}
names(Promoter.gene.ReMapchip.Overlap.byPromoter) <- names(ReMapChIP.GRanges.list)
#Create a  matrix where each row represents a promoter, each column represents  a ChIPed factor: entry[i,j] is the number of Chip-peaks of factor j that overlap with the promoter i
Promoter.gene.ReMapchip.Overlap.byPromoter.Mat <- matrix(0L,nrow = nrow(Promoter1kbdf.gene),ncol = length(Promoter.gene.ReMapchip.Overlap.byPromoter))
for(i in 1:ncol(Promoter.gene.ReMapchip.Overlap.byPromoter.Mat)){
  Promoter.gene.ReMapchip.Overlap.byPromoter.Mat[,i] <- unlist(lapply(Promoter.gene.ReMapchip.Overlap.byPromoter[[i]],length))
}
colnames(Promoter.gene.ReMapchip.Overlap.byPromoter.Mat) <- names(ReMapChIP.GRanges.list)
boxplot.matrix(Promoter.gene.ReMapchip.Overlap.byPromoter.Mat,las = 2)
#Create a binary matrix where each row represents a promoter, each column represents  a ChIPed factor: 1 means that promoter overlaps with at least one chip peak of that factor
Promoter.gene.ReMapchip.Overlap.byPromoter.Mat.Binary <- matrix(0L,nrow = nrow(Promoter1kbdf.gene),ncol = length(Promoter.gene.ReMapchip.Overlap.byPromoter))
Promoter.gene.ReMapchip.Overlap.byPromoter.Mat.Binary[Promoter.gene.ReMapchip.Overlap.byPromoter.Mat > 0] <- 1
colnames(Promoter.gene.ReMapchip.Overlap.byPromoter.Mat.Binary) <- names(ReMapChIP.GRanges.list)

#plot the histogram of overlaps between promoters and chip peaks
hist(rowSums(Promoter.gene.ReMapchip.Overlap.byPromoter.Mat.Binary[(rowSums(Promoter.gene.ReMapchip.Overlap.byPromoter.Mat.Binary) != 0),]),breaks = 200,main = "No. of overlapping Chip peaks for promoters(>0)",xlab = "")
hist(rowSums(Promoter.gene.ReMapchip.Overlap.byPromoter.Mat.Binary),breaks = 200,xlab = "No. of ChIPed factors overlapping with each promoter",main = "")
#########################################################################################################
#############################################Enhancer REMAP ChIP Interaction################################
#Get the interacting enhancers with a chip peak
Enhancer.ReMapchip.Interaction.byEnhancer <- list()#a list where each entry is the interaction list of the All enhancers with All peaks of one factor
for(i in 1:length(ReMapChIP.GRanges.list)){
  Enhancer.ReMapchip.Interaction.byEnhancer[[i]] <- InteractionExtractor(subjectCoord1=MCFEnhancersGR, subjectCoord2=ReMapChIP.GRanges.list[[i]] )
  print(i)
}
names(Enhancer.ReMapchip.Interaction.byEnhancer) <- names(ReMapChIP.GRanges.list)
#Create a  matrix where each row represents an enhancer, each column represents  a ChIPed factor: entry[i,j] is the number of Chip-peaks of factor j that interact with the enhancer i
Enhancer.ReMapchip.Interaction.byEnhancer.Mat <- matrix(0L, nrow=nrow(MCFEnhancersdf), ncol=length(Enhancer.ReMapchip.Interaction.byEnhancer))
for(i in 1:ncol(Enhancer.ReMapchip.Interaction.byEnhancer.Mat)){
  Enhancer.ReMapchip.Interaction.byEnhancer.Mat[,i] <- unlist(lapply(Enhancer.ReMapchip.Interaction.byEnhancer[[i]],length))
}
colnames(Enhancer.ReMapchip.Interaction.byEnhancer.Mat) <- names(ReMapChIP.GRanges.list)
boxplot.matrix(Enhancer.ReMapchip.Interaction.byEnhancer.Mat,las = 2)
#Create a binary matrix where each row represents an enhancer, each column represents  a ChIPed factor: 1 means that enhancer interacts with at least one chip peak of that factor
Enhancer.ReMapchip.Interaction.byEnhancer.Mat.Binary <- matrix(0L,nrow = nrow(MCFEnhancersdf),ncol = length(Enhancer.ReMapchip.Interaction.byEnhancer))
Enhancer.ReMapchip.Interaction.byEnhancer.Mat.Binary[Enhancer.ReMapchip.Interaction.byEnhancer.Mat > 0] <- 1
colnames(Enhancer.ReMapchip.Interaction.byEnhancer.Mat.Binary) <- names(ReMapChIP.GRanges.list)
#plot the histogram of interactions between Enhancers and chip peaks
hist(rowSums(Enhancer.ReMapchip.Interaction.byEnhancer.Mat.Binary[(rowSums(Enhancer.ReMapchip.Interaction.byEnhancer.Mat.Binary) != 0),]),breaks = 200,main = "No. of interacting Chip peaks for enhancers(>0)",xlab = "")
hist(rowSums(Enhancer.ReMapchip.Interaction.byEnhancer.Mat.Binary),breaks = 200,xlab = "No. of ChIPed factors interacting with each enhancers",main = "")

#########################################################################################################
#############################################Promoter REMAP ChIP Interaction################################
aa = OverlapInteractionBatchExtract(subjectCoord1=Promoter1kbGR.gene ,
                                    subjectCoord2List=ReMapChIP.GRanges.list,
                                    mode = "interaction" )
#a list where each entry is the interaction list of the All promoters with All peaks of one ChIPed factor
Promoter.gene.ReMapchip.Interaction.byPromoter <- aa[[1]]
#Create a  matrix where each row represents a promoter, each column represents  a ChIPed factor: entry[i,j] is the number of Chip-peaks of factor j that interact with the promoter i
Promoter.gene.ReMapchip.Interaction.byPromoter.Mat <- aa[[2]]
#Create a binary matrix where each row represents a promoter, each column represents  a ChIPed factor: 1 means that promoter interacts with at least one chip peak of that factor
Promoter.gene.ReMapchip.Interaction.byPromoter.Mat.Binary <- matrix(0L,nrow = nrow(Promoter1kbdf.gene),ncol = length(Promoter.gene.ReMapchip.Interaction.byPromoter))
Promoter.gene.ReMapchip.Interaction.byPromoter.Mat.Binary[Promoter.gene.ReMapchip.Interaction.byPromoter.Mat > 0] <- 1
colnames(Promoter.gene.ReMapchip.Interaction.byPromoter.Mat.Binary) <- names(ReMapChIP.GRanges.list)
#########################################################################################################
#########################################################################################################
#plot the number of enhancers/promoters overlapping/interacting  with at least i ER ChiP peaks
a = matrix(nrow = 50, ncol = 4)
i = 10
for(i in 1:50){
  aa = matrix(0L,nrow = nrow(MCFEnhancersdf),ncol = length(Enhancer.ReMapchip.Overlap.byEnhancer))
  aa[Enhancer.ReMapchip.Overlap.byEnhancer.Mat > (i-1)] <- 1
  a[i,1] <- sum(aa[,52])
  aaa = matrix(0L,nrow = nrow(Promoter1kbdf.gene),ncol = length(Promoter.gene.ReMapchip.Interaction.byPromoter))
  aaa[Promoter.gene.ReMapchip.Overlap.byPromoter.Mat > (i-1)] <- 1
  a[i,2] <- sum(aaa[,52])
  aa = matrix(0L,nrow = nrow(MCFEnhancersdf),ncol = length(Enhancer.ReMapchip.Interaction.byEnhancer))
  aa[Enhancer.ReMapchip.Interaction.byEnhancer.Mat > (i-1)] <- 1
  a[i,3] <- sum(aa[,52])
  aaa = matrix(0L,nrow = nrow(Promoter1kbdf.gene),ncol = length(Promoter.gene.ReMapchip.Interaction.byPromoter))
  aaa[Promoter.gene.ReMapchip.Interaction.byPromoter.Mat > (i-1)] <- 1
  a[i,4] <- sum(aaa[,52])
}
par(mfrow = c(2,2))
plot(c(1:50),a[,1], main = ""
     #, ylim = c(0,21000)
     , xlab = "min no of ER peaks overlapping", ylab="number of Enhancers")
plot(c(1:50),a[,2], main = ""
     #, ylim = c(0,21000)
     , xlab = "min no of ER peaks overlapping", ylab="number of promoters")
plot(c(1:50),a[,3], main = ""
     #, ylim = c(0,21000)
     , xlab = "min no of ER peaks interacting", ylab="number of Enhancers")
plot(c(1:50),a[,4], main = ""
     #, ylim = c(0,21000)
     , xlab = "min no of ER peaks interacting", ylab="number of promoters")

#########################################################################################################
#########################################################################################################
#Create two binary matrices 
#FIRST  for overlap     of bag of regulatory elements of genes and Chiped factors (numberofGenes*numberofChipedfactors)
AssociatedChip.gene.100kb.Overlap <-  ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList,
                                                           gene.enhancer=Vicinity100kb.Enhancer.plusPromoter.By.gene.Extended,
                                                           enhancerChipBinary=Enhancer.ReMapchip.Overlap.byEnhancer.Mat.Binary,
                                                           promoterChipBinary=Promoter.gene.ReMapchip.Overlap.byPromoter.Mat.Binary)
#sum of Binaries
colnames(AssociatedChip.gene.100kb.Overlap[[1]]) <- names(ReMapChIP.GRanges.list)
#Binary
colnames(AssociatedChip.gene.100kb.Overlap[[2]]) <- names(ReMapChIP.GRanges.list)
sum(AssociatedChip.gene.100kb.Overlap[[2]][,52])
#SECOND for interaction of bag of regulatory elements of genes and Chiped factors (numberofGenes*numberofChipedfactors)
#UPDATE: DECIDED NOT TO USE THE INTERACTION OF CHIP AND ENHANCERS OR PROMOTERS SINCE THEY WILL BE CONSIDRED IN THE BAG ANYWAYS
AssociatedChip.gene.100kb.Interaction <-  ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList,
                                                               gene.enhancer=Vicinity100kb.Enhancer.plusPromoter.By.gene.Extended,
                                                               enhancerChipBinary=Enhancer.ReMapchip.Interaction.byEnhancer.Mat.Binary,
                                                               promoterChipBinary=Promoter.gene.ReMapchip.Interaction.byPromoter.Mat.Binary)
#sum of binaries
colnames(AssociatedChip.gene.100kb.Interaction[[1]]) <- names(ReMapChIP.GRanges.list)
#binary
colnames(AssociatedChip.gene.100kb.Interaction[[2]]) <- names(ReMapChIP.GRanges.list)

length(unique(c(which(AssociatedChip.gene.100kb.Overlap[[2]][,52] == 1), which(AssociatedChip.gene.100kb.Interaction[[2]][,52] == 1))))
#############Computing same matrices for at least 10 overlap or interactions
aa = matrix(0L,nrow = nrow(MCFEnhancersdf),ncol = length(Enhancer.ReMapchip.Overlap.byEnhancer))
aa[Enhancer.ReMapchip.Overlap.byEnhancer.Mat > 9] <- 1
aaa = matrix(0L,nrow = nrow(Promoter1kbdf.gene),ncol = length(Promoter.gene.ReMapchip.Interaction.byPromoter))
aaa[Promoter.gene.ReMapchip.Overlap.byPromoter.Mat > 9] <- 1
AssociatedChip.gene.100kb.Overlap_Least10 <-  ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList,
                                                                   gene.enhancer=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                                                   enhancerChipBinary=aa,
                                                                   promoterChipBinary=aaa)
#sum of Binaries
colnames(AssociatedChip.gene.100kb.Overlap_Least10[[1]]) <- names(ReMapChIP.GRanges.list)
#Binary
colnames(AssociatedChip.gene.100kb.Overlap_Least10[[2]]) <- names(ReMapChIP.GRanges.list)

#########################################################################################################
#########################################################################################################
######################################                                 ##################################
######################################       Index of the genes        ##################################
######################################     associated with ER peak     ##################################
######################################                                 ##################################
#########################################################################################################
#########################################################################################################

Genes.Associated.REMAP.ER.index <- which(AssociatedChip.gene.100kb.Overlap_Least10[[2]][,52] > 0)
Genes.Associated.REMAP.ER.Entrez <- tssdf.genes$gene_id[Genes.Associated.REMAP.ER.index ]
#########################################################################################################
#########################################################################################################
#########################################################################################################

aa = matrix(0L,nrow = nrow(MCFEnhancersdf),ncol = length(Enhancer.ReMapchip.Interaction.byEnhancer))
aa[Enhancer.ReMapchip.Interaction.byEnhancer.Mat > 9] <- 1
aaa = matrix(0L,nrow = nrow(Promoter1kbdf.gene),ncol = length(Promoter.gene.ReMapchip.Interaction.byPromoter))
aaa[Promoter.gene.ReMapchip.Interaction.byPromoter.Mat > 9] <- 1
AssociatedChip.gene.100kb.Interaction_Least10 <-  ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList,
                                                                       gene.enhancer=Vicinity100kb.Enhancer.plusPromoter.By.gene.Extended,
                                                                       enhancerChipBinary=aa,
                                                                       promoterChipBinary=aaa)
#sum of binaries
colnames(AssociatedChip.gene.100kb.Interaction_Least10[[1]]) <- names(ReMapChIP.GRanges.list)
#binary
colnames(AssociatedChip.gene.100kb.Interaction_Least10[[2]]) <- names(ReMapChIP.GRanges.list)

length(unique(c(which(AssociatedChip.gene.100kb.Overlap_Least10[[2]][,52] == 1), which(AssociatedChip.gene.100kb.Interaction_Least10[[2]][,52] == 1))))
#########################################################################################################
#########################################################################################################
#######################################    ERa_RARa_Cobinding###############################################
#########################################################################################################
# Cooperative interaction between retinoic acid receptor-a and estrogen receptor in breast cancer: Ross-Innes
#read csv file and write bed file so it can be fed to lift and converted to hg19
RARa.ERa.Binding.hg18 <- read.delim("Chip-Data/ERa_RARa_coBinding.csv",header = T, sep = ",")
write.table(x = RARa.ERa.Binding.hg18, 
            file ="RARa.ERa.Binding.hg18.bed",
            quote = F, sep = "\t",
            row.names = F, col.names = F)
#read the hg19 bed file
RARa.ERa.Binding.hg19 <- read.delim("RARa.ERa.Binding.hg19.bed", header = F, sep = '\t', stringsAsFactors = F)
colnames(RARa.ERa.Binding.hg19) <- c("Chr", "start", "end")
RARa.ERa.Binding.hg19.Granges <- makeGRangesFromDataFrame(RARa.ERa.Binding.hg19)
#Find the enhancers, promoters and subsequently genes associated with these peaks:
aa <- OverlapInteractionBatchExtract(interactingA=Genome4DMCF7[,1:3],
                               interactingB=Genome4DMCF7[,4:6],
                               subjectCoord1 = Promoter1kbGR.gene,
                               subjectCoord2List = list(RARa.ERa.Binding.hg19.Granges),
                               mode = "both")
Promoter.gene.RAR.ER.Interaction.Overlap.byPromoter <- aa
Promoter.gene.RAR.ER.Interaction.Overlap.byPromoter.Mat <- Promoter.gene.RAR.ER.Interaction.Overlap.byPromoter$OverlapMat
aa <- OverlapInteractionBatchExtract(interactingA=Genome4DMCF7[,1:3],
                                     interactingB=Genome4DMCF7[,4:6],
                                     subjectCoord1 = MCFEnhancersGR,
                                     subjectCoord2List = list(RARa.ERa.Binding.hg19.Granges),
                                     mode = "both")
Enhancer.gene.RAR.ER.Interaction.Overlap.byEnhancer <- aa
Enhancer.gene.RAR.ER.Interaction.Overlap.byEnhancer.Mat <- Enhancer.gene.RAR.ER.Interaction.Overlap.byEnhancer$OverlapMat
Enhancer.gene.RAR.ER.Interaction.Overlap.byEnhancer.Mat[Enhancer.gene.RAR.ER.Interaction.Overlap.byEnhancer.Mat > 0] <- 1
AssociatedChip.gene.100kb.RAR.ER<-  ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList,
                                                         gene.enhancer=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                                         enhancerChipBinary=Enhancer.gene.RAR.ER.Interaction.Overlap.byEnhancer.Mat,
                                                         promoterChipBinary=Promoter.gene.RAR.ER.Interaction.Overlap.byPromoter.Mat)
Genes.Associated.RAR.ER.Entrez <- rownames(AssociatedChip.gene.100kb.RAR.ER$AssociatedPeaks)[AssociatedChip.gene.100kb.RAR.ER$AssociatedPeaks > 0]


#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################        expression         ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
#Divided into knock down dataset and time series dataset
#########################################################################################################
######################################        KD dataset         ########################################
#########################################################################################################
#knock down dataset
#1 Zeynep PaPE
GSE73663Inh <- read.delim("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression data/RNAseq/GSE73663_inhibitor.txt",header = T,stringsAsFactors = F)
GSE73663RNAseqMatInh <- as.matrix(GSE73663Inh[,2:19])
rownames(GSE73663RNAseqMatInh) <- GSE73663Inh$Gene.Symbol
GSE73663.Inh.RNAseq <- GSE73663RNAseqMatInh
aa = geneNameToEntrez(rownames(GSE73663.Inh.RNAseq))

#GSE86316 MEN1 KO
GSE86316RNAseq <- read.delim("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression data/RNAseq/GSE86316_Expression_transcripts_ShMEN1.txt",header = T,stringsAsFactors = F)
GSE86316RNAseqMat <- as.matrix(GSE86316RNAseq[,3:10])
rownames(GSE86316RNAseqMat) <- GSE86316RNAseq$RefseqID
GSE86316.RNAseq <- GSE86316RNAseqMat
aa = geneNameToEntrez(rownames(GSE86316.RNAseq))

#GSE80098_MCF7_RNAseq E2, R5020 or both for 12 hours: THis one only has differentially expressed genes, no raw data available
GSE80098RNAseq <- read.csv("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/RNAseq/GSE80098_MCF7_RNAseq.csv",header = T,stringsAsFactors = F)
GSE80098RNAseqmat <- as.matrix(GSE80098RNAseq[,c(2,3,6,9)])
rownames(GSE80098RNAseqmat) <- GSE80098RNAseq$RefSeq
GSE80098.RNAseq <- GSE80098RNAseqmat
aa = geneNameToEntrez(rownames(GSE80098.RNAseq))
match(TF.ENTREZ.Shrinked$ENTREZID, aa$ENTREZID)
###############G9 and PHF kd
GSE76507RNAseq <- read.delim("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression data/RNAseq/GSE76507_shG9a_shPHF20.pool.htseq.txt",header = T,stringsAsFactors = F)
colnames(GSE76507RNAseq) <- c("gene","ensembl","shG9a_NOE2","shG9a_E2","shNT_NOE2","shNT_E2","shNT_NoE2.rep1","shNT_NOE2.rep2","shNT_E2.rep1","shNT_E2.rep2","shPHF20_NOE2.rep1","shPHF20_NOE2.rep2","shPHF20_E2.rep1","shPHF20_E2.rep2")		
GSE76507RNAseqMat <- as.matrix(GSE76507RNAseq[,3:14])
rownames(GSE76507RNAseqMat) <- GSE76507RNAseq$gene
GSE76507.RNAseq<- GSE76507RNAseqMat
aa = geneNameToEntrez(GSE76507RNAseq$ensembl)

###############E2 or E2 plus something else
GSE68358RNAseq <- list()
GSE68358Files <- list.files(path="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression data/RNAseq/GSE68358_RAW", pattern=".txt$", full.names=T, recursive=FALSE)

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
GSE68358.RNAseq<- GSE68358RNAseqMat
aa = geneNameToEntrez(rownames(GSE68358.RNAseq))

a##################estradiol- and pro-inflammatory cytokine-treated MCF7 cells.# there is also chip seq and GRO-seq available for this set

GSE67295RNAseq <- read.delim("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression data/RNAseq/GSE67295GeneExpressionPolyA.txt",header = T,stringsAsFactors = F)
colnames(GSE67295RNAseq)[1] <- "TranscriptID"
GSE67295RNAseq$TranscriptID[1:100]
nrow(GSE67295RNAseq)
GSE67295RNAseqMat <- as.matrix(GSE67295RNAseq[,9:34])
rownames(GSE67295RNAseqMat) <- GSE67295RNAseq$TranscriptID
GSE67295.RNAseq <- GSE67295RNAseqMat
aa = geneNameToEntrez(rownames(GSE67295RNAseqMat))

################BRD4 kd
GSE55922RNAseq <- read.csv("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/RNAseq/GSE55922_normalized_counts.csv",header = T,stringsAsFactors = F)
GSE55922RNAseqMat <- as.matrix(GSE55922RNAseq[,2:13])
rownames(GSE55922RNAseqMat) <- GSE55922RNAseq$X
GSE55922.RNAseq <- GSE55922RNAseqMat
aa = geneNameToEntrez(rownames(GSE55922.RNAseq))

########CTCF knock down GSM2257523
GSE53532RNAseq00 <- read.delim(file = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression data/RNAseq/GSM2257523_Genes_fpkm_siNT_E0h_cufflinks.txt",header = F,sep = "\t",stringsAsFactors = F)
GSE53532RNAseq01 <- read.delim(file = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression data/RNAseq/GSM2257525_Genes_fpkm_siNT_E3h_cufflinks.txt",header = F,sep = "\t",stringsAsFactors = F)
GSE53532RNAseq10 <- read.delim(file = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression data/RNAseq/GSM2257524_Genes_fpkm_siCTCF_E0h_cufflinks.txt",header = F,sep = "\t",stringsAsFactors = F)
GSE53532RNAseq11 <- read.delim(file = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression data/RNAseq/GSM2257526_Genes_fpkm_siCTCF_E3h_cufflinks.txt",header = F,sep = "\t",stringsAsFactors = F)
GSE53532RNAseq <- as.matrix(cbind(as.numeric(GSE53532RNAseq00$V2),as.numeric(GSE53532RNAseq01$V2),as.numeric(GSE53532RNAseq10$V2),as.numeric(GSE53532RNAseq11$V2)))
colnames(GSE53532RNAseq) <- c("siNT_E0h", "siNT_E3h", "siCTCF_E0h", "siCTCF_E3h")
GSM2257523.RNAseq <- GSE53532RNAseq
aa = geneNameToEntrez(rownames(GSM2257523.RNAseq))
########R4K1 stapled peptide
GSE108308.RNAseq <- read.delim(file = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression data/GSE108308_Process.RNASeq.DE.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(GSE108308.RNAseq)[1] <- "TranscriptName"
GSE108308.RNAseq2 <- GSE108308.RNAseq[,9:20]
rownames(GSE108308.RNAseq2) <- GSE108308.RNAseq[,1]
GSE108308.RNAseq <- GSE108308.RNAseq2
GSE108308.RNAseq <- as.matrix(GSE108308.RNAseq)
remove(GSE108308.RNAseq2)
aa = geneNameToEntrez(rownames(GSE108308.RNAseq))
###########Sumo # this is tpm not rpkm #couldn't convert the names to ENTREZ. Drop for now.
GSM2915403RNAseq <- list()
GSM2915403Files <- list.files(path="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression data/Sumo", pattern=".tsv$", full.names=T, recursive=FALSE)

for( i in 1:length(GSM2915403Files)){
  GSM2915403RNAseq[[i]] <- read.delim(GSM2915403Files[i],header = T,stringsAsFactors = F)
}
names(GSM2915403RNAseq) <-  GSM2915403Files
GSM2915403RNAseqMat <- matrix(0L,nrow = nrow(GSM2915403RNAseq[[1]]), ncol = length(GSM2915403Files))
for(i in 1:length(GSM2915403RNAseq)){
  GSM2915403RNAseqMat[,i] <- GSM2915403RNAseq[[i]]$tpm
}
colnames(GSM2915403RNAseqMat)<- GSM2915403Files
rownames(GSM2915403RNAseqMat)<- GSM2915403RNAseq[[1]]$target_id
GSM2915403.RNAseq<- GSM2915403RNAseqMat
colnames(GSM2915403.RNAseq) <- c("E-1","E-2","E-3", "ICI-1","ICI-2","ICI-3", "OHT-1","OHT-2", "OHT-3", "V-1","V-2","V-3")
aa = geneNameToEntrez(rownames(GSM2915403.RNAseq))
#########################################################################################################
#########################################################################################################
#####################################      conversion to entrez      ####################################
#####################################   Extraction of mutual genes   ####################################
#########################################################################################################
#########################################################################################################
#list of used datasets:
Exp.kd.Raw.Dataset.list <- list(GSE73663.Inh.RNAseq  #1 Zeynep PaPe
                               #,GSE86316.RNAseq      #2 MEN1 KO
                               #,GSE55922.RNAseq      #3 BRD4 kd
                               ,GSE67295.RNAseq      #4 estradiol- and pro-inflammatory cytokine-treated
                               ,GSE68358.RNAseq      #5 E2 or E2 plus something else
                               ,GSE76507.RNAseq      #6 G9 and PHF kd
                               ,GSE108308.RNAseq     #7 R4K1 stapled peptide
                                )
names(Exp.kd.Raw.Dataset.list) <- c("GSE73663"
                                    #,"GSE86316"
                                    #,"GSE55922"
                                    ,"GSE67295"
                                    ,"GSE68358"
                                    ,"GSE76507"
                                    ,"GSE108308")

#Omit the genes that are not common between all datasets
Exp.kd.Raw.Dataset.list.Common <- getCommonGenes(Exp.kd.Raw.Dataset.list)
names(Exp.kd.Raw.Dataset.list.Common) <- c("GSE73663"
                                           #,"GSE86316"
                                           #,"GSE55922"
                                           ,"GSE67295"
                                           ,"GSE68358"
                                           ,"GSE76507"
                                           ,"GSE108308")
#put all datasets into one large matrix
Exp.kd.Raw.Dataset.matrix.Common <- do.call(cbind, Exp.kd.Raw.Dataset.list.Common)
Exp.kd.Raw.Dataset.matrix.Common <- as.matrix(Exp.kd.Raw.Dataset.matrix.Common)

#Exp.kd.Raw.Dataset.matrix.Common <- Exp.kd.Raw.Dataset.matrix.Common[-c(which(rowSums(Exp.kd.Raw.Dataset.matrix.Common) == 0)),]

aa <- list(GSE73663.Inh.RNAseq  #1 Zeynep PaPe
           ,GSE86316.RNAseq      #2 MEN1 KO
           ,GSE55922.RNAseq      #3 BRD4 kd
           #,GSE67295.RNAseq      #4 estradiol- and pro-inflammatory cytokine-treated
           #,GSE68358.RNAseq      #5 E2 or E2 plus something else
           ,GSE76507.RNAseq      #6 G9 and PHF kd
           ,GSE108308.RNAseq     #7 R4K1 stapled peptide
)
aaa <- getCommonGenes(aa)

match(rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), rownames(aaa[[1]]))
#########################################################################################################
#########################################################################################################
#####################################             Remove           ######################################
#####################################          Batch effect        ######################################
#########################################################################################################
#########################################################################################################
#done in batch.R
#removed genes which had variance less than 1 across al experiments, performed batch removal using Combat
Exp.kd.Batch.Dataset.matrix.Common.vgt1

#Construct list of raw data after removing all genes with 
#variance less than one across all experiments
#containing the same genes as Exp.kd.Batch.Dataset.matrix.Common.vgt1
a <- apply(Exp.kd.Raw.Dataset.matrix.Common, 1, var)
Exp.kd.Raw.Dataset.list.Common.vgt1 <- list()
for(i in 1:length(Exp.kd.Raw.Dataset.list.Common)){
  Exp.kd.Raw.Dataset.list.Common.vgt1[[i]] <- Exp.kd.Raw.Dataset.list.Common[[i]][-c(which(a < 1)),]
}

#Construct list of batch data after removing all genes with 
#variance less than one across all experiments
#containing the same genes as Exp.kd.Batch.Dataset.matrix.Common.vgt1
Exp.kd.Batch.Dataset.list.Common.vgt1 <- list()
aa <- unlist(lapply(Exp.kd.Raw.Dataset.list.Common, ncol))
a <- numeric(0)
for(i in 1:length(aa)){
  a <- c(a, rep(i,aa[i]))
}
for (i in 1:length(Exp.kd.Raw.Dataset.list.Common.vgt1)){
  Exp.kd.Batch.Dataset.list.Common.vgt1[[i]] <- Exp.kd.Batch.Dataset.matrix.Common.vgt1[,which(a == i)]
}

#########################################################################################################
#########################################################################################################
#####################################             Organize          #####################################
#####################################           Experiments         #####################################
#########################################################################################################
#########################################################################################################
#Done in ExperimentSetup.R
DataSet.Experiment.matrix.Raw.vgt1.ER # raw expression values ER associated genes
DataSet.Experiment.matrix.Raw.vgt1.LFC.ER # raw log fold change ER associated genes
DataSet.Experiment.matrix.Batch.vgt1.ER # batch removed expression values ER associated genes
DataSet.Experiment.matrix.Batch.vgt1.LFC.ER # batch removed log fold change ER associated genes
#########################################################################################################
#########################################################################################################
#write to a file to input GEMSTAT
ExpressionWriter(DataSet.Experiment.matrix.Batch.vgt1.LFC.ER, output.File.Name="Gene_Expression")
#########################################################################################################
#########################################################################################################

#########################################################################################################
#########################################################################################################
#####################################             Clustering              ###############################
#####################################         Expression profiles         ###############################
#########################################################################################################
#########################################################################################################
#cluster the expression profile (log fold change) of ER associated genes
aa = clustrator(DataSet.Experiment.matrix.Batch.vgt1.LFC.ER,nfolds = 3,plot_mfrow = c(3,1),Exp.metadata = DataSet.Experiment.Metadata)
aa <- list()
for(i in 3:15){
  aa[[i]] <- list()
  for(j in 1:10){
    aa[[i]][[j]] = clustrator(DataSet.Experiment.matrix.Batch.vgt1.LFC.ER,nfolds = i,
                              Exp.metadata = DataSet.Experiment.Metadata,exportPlot = T,
                              filename = paste("GeneClustering",i,"_",j,".pdf",sep = ""))
  }
}

#########################################################################################################
#########################################################################################################
#########################################################################################################
######################################        TS dataset         ########################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
# This dataset contains Time series data after estrogen treatment:
# Three independent biological replicates (30 samples: 3 replicates x 10 time points) of MCF7 cells were exposed to
# 10nM Estradiol for 0, 1, 2, 3, 4, 5, 6, 8, 12 , or 24 hours,
# what genes to use? 
# two sets of genes:
# 1. ER-RAR associated and differentially expressed genes
# 2. ER associated, not ER-RAR associated, and differentially expressed

# Genes.Associated.REMAP.ER.Entrez : name of ER associated genes (REMAP)

# Genes.Associated.RAR.ER.Entrez: RAR-ER associated genes names (ROSS-INES)

#read dataset
GSE78167RNAseq <- list()
aa <- list.files(path="Expression data/RNAseq/GSE78169_RAW", pattern=".txt$", full.names=T, recursive=FALSE)

for( i in 1:length(aa)){
  GSE78167RNAseq[[i]] <- read.delim(aa[i], header = T, stringsAsFactors = F)
}
names(GSE78167RNAseq) <-  aa
GSE78167.RNAseq <- matrix(0L,nrow = nrow(GSE78167RNAseq[[1]]), ncol = length(aa))
for(i in 1:length(GSE78167RNAseq)){
  GSE78167.RNAseq[,i] <- GSE78167RNAseq[[i]]$TPM
}

aa <- character()
for (i in 1:nrow(GSE78167.RNAseq)){
  aa[i] <- stri_split_coll(str = GSE78167RNAseq[[1]]$gene_id[i],pattern = "|")[[1]][2]
}
rownames(GSE78167.RNAseq) <- aa
#separate replicates
GSE78167.RNAseq.Rep1 <- GSE78167.RNAseq[, 1:10]
GSE78167.RNAseq.Rep2 <- GSE78167.RNAseq[, 11:20]
GSE78167.RNAseq.Rep3 <- GSE78167.RNAseq[, 21:30]
colnames(GSE78167.RNAseq.Rep1) <- c("0", "1", "2", "3", "4", "5", "6", "8", "12", "24")
colnames(GSE78167.RNAseq.Rep2) <- c("0", "1", "2", "3", "4", "5", "6", "8", "12", "24")
colnames(GSE78167.RNAseq.Rep3) <- c("0", "1", "2", "3", "4", "5", "6", "8", "12", "24")

GSE78167.RNAseq.Avg <- (GSE78167.RNAseq.Rep1 + GSE78167.RNAseq.Rep2 + GSE78167.RNAseq.Rep3)/3
remove(GSE78167RNAseqMat, GSE78167RNAseq, GSE78167ExpName, GSE78167GeneSymbol, GSE78167GeneEntrez)
aa <- apply(GSE78167.RNAseq.Avg, 1, max)
GSE78167.RNAseq.Avg.Norm <- GSE78167.RNAseq.Avg/aa
#########read differentially expressed gene matrix
#########
a <- read.csv("Expression data/RNAseq/TimeSeriesDiffExpNormalized.csv",header = T,stringsAsFactors = F)
aa <- geneNameToEntrez(a$Gene.Name)
aa$ENTREZID[866] <- "57491"
aa$ENTREZID[1300] <- "64699"
aaa <- match(a$Gene.Name, aa$ALIAS)
a <- a[!is.na(aaa),]
aaa <- match(a$Gene.Name, aa$ALIAS)
rownames(a) <- aa$ENTREZID[aaa]
aa <- lapply(X = a, FUN = gsub, pattern =",", replacement = "")
aa <- do.call(cbind, aa)
aaa <- aa[,3:ncol(aa)]
aaa <- apply(X = aaa, MARGIN = 2, FUN = as.numeric)
aa <- cbind(aa[,1:2], aaa)
rownames(aa) <- rownames(a)
GSE78167.RNAseq.Diff <- aa
colnames(GSE78167.RNAseq.Diff) <- c("gene_Alias","cluster","0","1","3","4","5","6","8","12","24")
########################################################################################################################
########################### Find ER associated genes that are differentially expressed in TS dataset ####################
aa <- match(Genes.Associated.REMAP.ER.Entrez, rownames(GSE78167.RNAseq.Diff))
aa <- aa[!is.na(aa)]
Genes.Associated.REMAP.ER.Entrez.TS.DIff <- rownames(GSE78167.RNAseq.Diff)[aa]
################# Find ER-RAR associated genes that are differentially expressed in TS dataset ##########################
aa <- match(Genes.Associated.RAR.ER.Entrez, rownames(GSE78167.RNAseq.Diff))
aa <- aa[!is.na(aa)]
Genes.Associated.RAR.ER.Entrez.TS.DIff <- rownames(GSE78167.RNAseq.Diff)[aa]
########################################################################################################################
#Get expression matrix for total of 210 genes: use Genes.Associated.RAR.ER.Entrez.TS.DIff + 105 of setdiff(Genes.Associated.REMAP.ER.Entrez.TS.DIff, Genes.Associated.RAR.ER.Entrez.TS.DIff)
Genes.Associated.REMAP.ER.Entrez.TS.DIff.NoRAR.Sampled <- sample(x = setdiff(Genes.Associated.REMAP.ER.Entrez.TS.DIff,
                                                                             Genes.Associated.RAR.ER.Entrez.TS.DIff), size = 105, replace = F)
Genes.Associated.REMAP.RAR.ER.Entrez.TS.Diff.Union <- c(Genes.Associated.REMAP.ER.Entrez.TS.DIff.NoRAR.Sampled,
                                                        Genes.Associated.RAR.ER.Entrez.TS.DIff)
########################################################################################################################
#Averaged Expression Matrix of Chosen genes:
aa <- match(Genes.Associated.REMAP.RAR.ER.Entrez.TS.Diff.Union, rownames(GSE78167.RNAseq.Avg))
GSE78167.RNAseq.Avg.ChosenGenes <- GSE78167.RNAseq.Avg[aa, ]
#max normalize the expression matrix:
aa <- apply(GSE78167.RNAseq.Avg.ChosenGenes, 1, max)
GSE78167.RNAseq.Avg.ChosenGenes.Norm <- GSE78167.RNAseq.Avg.ChosenGenes/aa
####################### USE clustrator on this dataset
aa <- clustrator(ExpressionMat =GSE78167.RNAseq.Avg.ChosenGenes.Norm, nfolds = 8 )
aa1 <- aa
####################### Find the genes whose expression drops at the second or third time point significantly and remove them for the new dataset
aa <- GSE78167.RNAseq.Avg.ChosenGenes.Norm[,2] - GSE78167.RNAseq.Avg.ChosenGenes.Norm[,1]
aaa <- GSE78167.RNAseq.Avg.ChosenGenes.Norm[,3] - GSE78167.RNAseq.Avg.ChosenGenes.Norm[,1]
aaaa <- GSE78167.RNAseq.Avg.ChosenGenes.Norm[,3] - GSE78167.RNAseq.Avg.ChosenGenes.Norm[,2]

aaaaa <- (aa < -0.05) + (aaa < -0.05) + (aaaa < -0.05)
GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed <- GSE78167.RNAseq.Avg.ChosenGenes.Norm[aaaaa == 0, ]
aa <- clustrator(ExpressionMat =GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed, nfolds = 5 )


######################
# find other TFs that are associated with these genes
aa = matrix(0L,nrow = nrow(MCFEnhancersdf),ncol = length(Enhancer.ReMapchip.Overlap.byEnhancer))
aa[Enhancer.ReMapchip.Overlap.byEnhancer.Mat > 0] <- 1
aaa = matrix(0L,nrow = nrow(Promoter1kbdf.gene),ncol = length(Promoter.gene.ReMapchip.Interaction.byPromoter))
aaa[Promoter.gene.ReMapchip.Overlap.byPromoter.Mat > 0] <- 1
AssociatedChip.gene.100kb.Overlap_Least1 <-  ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList,
                                                                   gene.enhancer=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                                                   enhancerChipBinary=aa,
                                                                   promoterChipBinary=aaa)
colnames(AssociatedChip.gene.100kb.Overlap_Least1[[1]]) <- colnames(Enhancer.ReMapchip.Overlap.byEnhancer.Mat)
colnames(AssociatedChip.gene.100kb.Overlap_Least1[[2]]) <- colnames(Enhancer.ReMapchip.Overlap.byEnhancer.Mat)
aa <- match(Genes.Associated.REMAP.RAR.ER.Entrez.TS.Diff.Union,
            rownames(AssociatedChip.gene.100kb.Overlap_Least1$AssociatedPeaks))
aaa <- AssociatedChip.gene.100kb.Overlap_Least1$AssociatedPeaksBinary[aa,]
colnames(aaa) <- colnames(Enhancer.ReMapchip.Overlap.byEnhancer.Mat)
sort(colSums(aaa)) #Didn't find anything significantly different from average
##################################################################################
######################################### Write expression of TFs
#only using RAR and ER as TFs
#setting concentration of ER to 1 at every time point except the zero time point
a <- geneNameToEntrez(c("RARA"))
aa <- match(a$ENTREZID, rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm))
TF.Expression.TimeSeries <- rbind(c(0,rep(1,9)), GSE78167.RNAseq.Avg.ChosenGenes.Norm[aa,])
rownames(TF.Expression.TimeSeries) <- c("ESR1","RARA")

##################################################################################
#create motif count matrix for ESR1 and RARA : use half sites for both of them
TF.motifs.TimeSeries <- list()
aa <- list.files(path="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Motifs/Chosen_timeSeries/", full.names=T, recursive=FALSE)
for( i in 1:length(aa)){
  TF.motifs.TimeSeries[[i]] <- read.table(aa[i], header=T, sep ="\t")
  TF.motifs.TimeSeries[[i]] =  TF.motifs.TimeSeries[[i]][,2:5]
}

names(TF.motifs.TimeSeries) <- c("ESR1","RARA")
TF.motifs.TimeSeries.count <- lapply(TF.motifs.TimeSeries, PWMtoCount)
#########
#write the motifs to a file:
MotifWriter(motif.List = TF.motifs.TimeSeries.count, output.File.Name = "motifs_timeseries")




############  EXP6 ####### #with cooperativity: EXP6
GEMSTAT_input_constructor(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm))),
                          gene_expression_Mat=GSE78167.RNAseq.Avg.ChosenGenes.Norm, 
                          .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                          .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                          .promoter.promoter.int=Promoter.gene.PromoterIntList,
                          .promoter.sequence=Promoter.Sequence.Char,
                          .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                          experiment.nu = 6,
                          .motifs=TF.motifs.TimeSeries.count,
                          .TF_expression=TF.Expression.TimeSeries,
                          .max.enh.length = numeric(0),
                          coopertingTFIndex = rbind(c(1,1),c(1,2)))
########Reading the output of GEMSTAT with Coop Sum of enhancer outputs
GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop <- read.table("~/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Experiment_6/Cluster_1/Outputs/Experiment_6_Cluster_1.out", header = T, sep = "\t")
GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop <- GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop[seq(2,nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop), 2),]
rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop) <- GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop[,1]
GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop <- GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop[,2:ncol(GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop)]
GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop <- do.call(cbind, GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop)
######Compare prediction and real values
aa <- match(rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop),
            rownames(GSE78167.RNAseq.Avg.ChosenGenes.Norm))
# GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Sum.Coop.PerGene.SSE <- numeric(length = nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm))
GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Sum.Coop.SSE <- (rowSums((GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop -
                                                                               GSE78167.RNAseq.Avg.ChosenGenes.Norm)^2))/ncol(GSE78167.RNAseq.Avg.ChosenGenes.Norm)
#plot the prediction VS real value for each gene
par(mfrow = c(10,10), mar = c(0.1, 0.1 , 0.1 , 0.1))
for(i in 1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop)){
  print(i)
  plot(GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Output.Sum.Coop[i,],
       ylim = c(0,1), col = 2,
       xaxt = "n", yaxt = "n", main ="", ylab = "", xlab = "",
       type = "l", lwd = 2)

  lines(GSE78167.RNAseq.Avg.ChosenGenes.Norm[i,], ylim = c(0,1), col = 3)
  aa <- GSE78167.RNAseq.Avg.ChosenGenes.Norm.GEMSTAT.Sum.Coop.SSE[i]
  text((ncol(GSE78167.RNAseq.Avg.ChosenGenes.Norm) - 2), 0.9 , format(round(aa, 2), nsmall = 2), cex = 0.5) # write the size of the SSE
}
############  EXP7 ####### same settings as Exp6 but using correlation as obj function
GEMSTAT_input_constructor(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm))),
                          gene_expression_Mat=GSE78167.RNAseq.Avg.ChosenGenes.Norm, 
                          .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                          .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                          .promoter.promoter.int=Promoter.gene.PromoterIntList,
                          .promoter.sequence=Promoter.Sequence.Char,
                          .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                          experiment.nu = 7,
                          .motifs=TF.motifs.TimeSeries.count,
                          .TF_expression=TF.Expression.TimeSeries,
                          .max.enh.length = numeric(0),
                          coopertingTFIndex = rbind(c(1,1),c(1,2)),
                          ..oo = "CORR")
############  reading the output
Exp7.GEMSTAT.output <- GEMSTAT_output_Reader(.exp.nu = 7, .number.Of.Clusters = 1,
                                             .root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT",
                                             .dest_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                                             plot.Res = T)

############  EXP8 ####### same settings as Exp6 but using threshold distance on cooperativity 
GEMSTAT_input_constructor(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm))),
                          gene_expression_Mat=GSE78167.RNAseq.Avg.ChosenGenes.Norm, 
                          .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                          .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                          .promoter.promoter.int=Promoter.gene.PromoterIntList,
                          .promoter.sequence=Promoter.Sequence.Char,
                          .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                          experiment.nu = 8,
                          .motifs=TF.motifs.TimeSeries.count,
                          .TF_expression=TF.Expression.TimeSeries,
                          .max.enh.length = numeric(0),
                          coopertingTFIndex = rbind(c(1,1),c(1,2)),
                          ..oo = "SSE", ..ct = 13)
#read output
Exp8.GEMSTAT.output <- GEMSTAT_output_Reader(.exp.nu = 8, .number.Of.Clusters = 1,
                                             .root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT",
                                             .dest_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                                             plot.Res = T)
############  EXP9 ####### same settings as Exp8 but add more TFs: PGR, AR, LEF1, AP1, NKX3-1

######################################### Write expression of TFs
#only using RAR and ER as TFs
#setting concentration of ER to 1 at every time point except the zero time point
a <- geneNameToEntrez(c("PGR", "AR", "LEF1", "JUN", "NKX3-1"))
aa <- match(a$ENTREZID, rownames(GSE78167.RNAseq.Avg.Norm))
TF.Expression.TimeSeries.Exp9 <- rbind(TF.Expression.TimeSeries, GSE78167.RNAseq.Avg.Norm[aa,])
rownames(TF.Expression.TimeSeries.Exp9) <- c("ESR1","RARA", "PGR", "AR", "LEF1", "JUN_1", "NKX3-1" )

##################################################################################
#create motif count matrix for ESR1 and RARA : use half sites for both of them +  PGR, AR, LEF1, AP1, NKX3-1
aa <- match( c("PGR", "AR", "LEF1", "JUN_1", "NKX3-1" ), names(TF.motifs))
TF.motifs.TimeSeries.Exp9 <- c(TF.motifs.TimeSeries, TF.motifs[aa])
names(TF.motifs.TimeSeries.Exp9) <- c("ESR1","RARA", "PGR", "AR", "LEF1", "JUN_1", "NKX3-1")
TF.motifs.TimeSeries.Exp9.count <- lapply(TF.motifs.TimeSeries.Exp9, PWMtoCount)

GEMSTAT_input_constructor(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm))),
                          gene_expression_Mat=GSE78167.RNAseq.Avg.ChosenGenes.Norm, 
                          .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                          .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                          .promoter.promoter.int=Promoter.gene.PromoterIntList,
                          .promoter.sequence=Promoter.Sequence.Char,
                          .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                          experiment.nu = 9,
                          .motifs=TF.motifs.TimeSeries.Exp9.count,
                          .TF_expression=TF.Expression.TimeSeries.Exp9,
                          .max.enh.length = numeric(0),
                          coopertingTFIndex = rbind(c(1,1),c(1,2)),
                          ..oo = "SSE", ..ct = 13)
#read output
Exp9.GEMSTAT.output <- GEMSTAT_output_Reader(.exp.nu = 9, .number.Of.Clusters = 1,
                                             .root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT",
                                             .dest_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                                             plot.Res = T)
############  EXP10 ####### same settings as Exp8 but removed genes that were reppressed in the early time points
GEMSTAT_input_constructor(clusterIndex = list(c(1:nrow(GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed))),
                          gene_expression_Mat=GSE78167.RNAseq.Avg.ChosenGenes.Norm.NotRepressed, 
                          .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                          .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                          .promoter.promoter.int=Promoter.gene.PromoterIntList,
                          .promoter.sequence=Promoter.Sequence.Char,
                          .dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                          experiment.nu = 10,
                          .motifs=TF.motifs.TimeSeries.count,
                          .TF_expression=TF.Expression.TimeSeries,
                          .max.enh.length = numeric(0),
                          coopertingTFIndex = rbind(c(1,1),c(1,2)),
                          ..oo = "SSE", ..ct = 13)
#read output
Exp10.GEMSTAT.output <- GEMSTAT_output_Reader(.exp.nu = 10, .number.Of.Clusters = 1,
                                             .root_dir="tabebor2@veda.cs.illinois.edu:/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT",
                                             .dest_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT",
                                             plot.Res = T)
#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################            TF             ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
#the list of the TFs, their PWMs, their expression profiles
#List of TFs:
TF.Names <- c(
"FOXA1",
"GATA3" 
,"NFIB"
,"YBX1"
,"AR"
,"NR3C1"  #GR 
,"PGR"
,"PBX1"
#, "CTCF"
,"JUN"     #"AP1"
,"TFAP2C"  #"AP-2g" 
,"SP1"
,"RELA"    #NFkb
,"CEBPB"
,"RARG"
,"RXRA"
#,"PAX6"
,"PAX2"
#,"PITX2"
,"NKX3-1"
,"LEF-1" 
,"NR5A2" #"LRH-1"
,"RUNX1" 
#,"ZNF143"
,"POU5F1" #OCT-4
#,"NR2F1" 
#,"NR4A2" 
,"MYC" 
,"MAX" 
#,"XBP1" 
,"ESR1"
,"PPARD"
,"RARA"
,"RXRB"
,"NR2F2"
,"NR2C1"
,"PPARG")

TF.ENTREZ <- geneNameToEntrez(sort(TF.Names))
#TF.ENTREZ <- TF.ENTREZ$ENTREZID
#TF.ENTREZ <- c(TF.ENTREZ, "5468")
aa <- intersect(TF.ENTREZ$ENTREZID, rownames(DataSet.Experiment.Matrix.Batch.vgt1.logfold) )
aa <- intersect(TF.ENTREZ$ENTREZID, rownames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile))

setdiff(TF.ENTREZ$ENTREZID, aa) #0
#########################################################################################################
#########################################################################################################
#######################################       SHRINKED TFs       ########################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

TF.Names.Shrinked <- c(
  "FOXA1",
  "GATA3" 
#  ,"NFIB"
#  ,"YBX1"
  ,"AR"
  ,"NR3C1"  #GR 
  ,"PGR"
  ,"PBX1"
  #, "CTCF"
  ,"JUN"     #"AP1"
  ,"TFAP2C"  #"AP-2g" 
  ,"SP1"
#  ,"RELA"    #NFkb
  ,"CEBPB"
  ,"RARG"
  ,"RXRA"
  #,"PAX6"
#  ,"PAX2"
  #,"PITX2"
  ,"NKX3-1"
#  ,"LEF-1" 
  ,"NR5A2" #"LRH-1"
  ,"RUNX1" 
  #,"ZNF143"
#  ,"POU5F1" #OCT-4
  #,"NR2F1" 
  #,"NR4A2" 
#  ,"MYC" 
#  ,"MAX" 
  #,"XBP1" 
  ,"ESR1"
#  ,"PPARD"
  ,"RARA"
#  ,"RXRB"
#  ,"NR2F2"
#  ,"NR2C1"
#  ,"PPARG"
)

TF.ENTREZ.Shrinked <- geneNameToEntrez(sort(TF.Names.Shrinked))
TF.ENTREZ.Shrinked$ENTREZID[1] <- "367"
TF.ENTREZ.Shrinked$ENTREZID[17] <- "6667"
#TF.ENTREZ <- TF.ENTREZ$ENTREZID
#TF.ENTREZ <- c(TF.ENTREZ, "5468")
aa <- intersect(TF.ENTREZ.Shrinked$ENTREZID, rownames(DataSet.Experiment.Matrix.Batch.vgt1.logfold) )
aa <- intersect(TF.ENTREZ.Shrinked$ENTREZID, rownames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile))

setdiff(TF.ENTREZ.Shrinked$ENTREZID, aa) #0
#########################################################################################################
#########################################################################################################
#######################################       SHRINKED TFs       ########################################
#######################################        Enrichment        ########################################
#########################################################################################################
#########################################################################################################
#write thr bag of regulatory elements for all genes associated with ER, in order to feed to GEMSTAT and get the annotation.
aa = which(rownames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20) %in% Genes.Associated.REMAP.ER.Entrez)
aaa <- match(rownames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20)[aa], names(Vicinity100kb.Enhancer.plusPromoter.By.gene))
Gene.ER.Enhancer.Length <- WriteFastaOfBag(enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene[aaa],
                                           Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                           promoter.promoter.int=Promoter.gene.PromoterIntList[aaa],
                                           promoter.sequence=Promoter.Sequence.Char,
                                           return.Length = T)
ExpressionWriter(DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.ER,
                 output.File.Name = paste("ExpressionMat_",
                                          0, sep=""))
#read the annotation output
Gene.ER.Enhancer.annotation.19TFs <- readLines("annotation.txt")
aaa <- which(Gene.ER.Enhancer.annotation.19TFs == "19 19 92 92")
Gene.ER.Enhancer.annotation.19TFs <- Gene.ER.Enhancer.annotation.19TFs[-c(1:aaa)]
aaa <- which( Gene.ER.Enhancer.annotation.19TFs== "Parameters for running the program: ")
Gene.ER.Enhancer.annotation.19TFs <- Gene.ER.Enhancer.annotation.19TFs[1:(aaa-1)]
aa <- lapply(X = Gene.ER.Enhancer.annotation.19TFs, FUN = strsplit, split = "\t")
aaa <- lapply(X = aa, FUN = unlist)
Gene.ER.Enhancer.annotation.19TFs <- do.call(rbind, aaa)
Gene.ER.Enhancer.annotation.19TFs[,2] <- as.integer(Gene.ER.Enhancer.annotation.19TFs[,2])
Gene.ER.Enhancer.annotation.19TFs <- data.frame(gene = Gene.ER.Enhancer.annotation.19TFs[,1],
                                                position = as.integer(Gene.ER.Enhancer.annotation.19TFs[,2]),
                                                motif = Gene.ER.Enhancer.annotation.19TFs[,3],
                                                LLRscore = as.numeric(Gene.ER.Enhancer.annotation.19TFs[,4]),
                                                maxLLR = as.numeric(Gene.ER.Enhancer.annotation.19TFs[,5]))
#Get enhancer number for each motif hit and add it to the annotation dataframe
aa = c(1, Gene.ER.Enhancer.annotation.19TFs$position)
aaa = c(Gene.ER.Enhancer.annotation.19TFs$position, 1)
a = aaa - aa
a = a[1:(length(a)-1)]

aa1 <- c(1, as.numeric(Gene.ER.Enhancer.annotation.19TFs$gene))
aaa1 <- c(as.numeric(Gene.ER.Enhancer.annotation.19TFs$gene), 1)
a1 <- (aa1 == aaa1)
a1[1] <- T
a1 <- a1[1:(length(a1) -1)]

a2 <- which(a < 0)
a2 <- c(1, a2)
a2 = c(a2, (length(a) + 1))
aa2 <- which(a1 == 0)
aa2 <- c(1, aa2)
aa3 <- numeric(length(a))
for(i in 1:(length(a2) - 1)){
  if (a2[i] %in% aa2){
    cnt <- 1
  }
  aa3[a2[i]:(a2[i+1] -1)]<- rep(cnt, length(aa3[a2[i]:(a2[i+1] -1)]))
  cnt <- cnt + 1
}

#################### DATA FRAME OF ANNOTATIONS############
Gene.ER.Enhancer.annotation.19TFs <- cbind(Gene.ER.Enhancer.annotation.19TFs, aa3)
colnames(Gene.ER.Enhancer.annotation.19TFs)[6] <- "reg.el.number"
###################
# Now get how many hits of each motif for each enhancer of each gene: 
# create a list: one entry per gene, the entry is a matrix: each row represents a
# motif and each column represents a regulatory element: the entries are the number of
# motif hits in that reg element
motifHits.ByGene.RegElement <- list()
number.of.motifs <- 19
#
for (i in 1:length(levels(Gene.ER.Enhancer.annotation.19TFs$gene))){
  # get the subset of dataset for this gene
  cur.annotation.frame <- subset(Gene.ER.Enhancer.annotation.19TFs, gene == levels(Gene.ER.Enhancer.annotation.19TFs$gene)[i], select = c( "motif", "reg.el.number"))
  # number of regulatory elements
  cur.reg.nu <- length(unique(cur.annotation.frame$reg.el.number))
  motifHits.ByGene.RegElement[[i]] <- matrix(0L, nrow=number.of.motifs, ncol=cur.reg.nu)
  for(reg.el in 1:cur.reg.nu){
    cur.table <- table(cur.annotation.frame$motif[cur.annotation.frame$reg.el.number == reg.el])
    motifHits.ByGene.RegElement[[i]][, reg.el] <- as.numeric(cur.table)
    rownames(motifHits.ByGene.RegElement[[i]]) <- names(cur.table)
  }
}
rm(i, cur.annotation.frame, cur.reg.nu, number.of.motifs, reg.el, cur.table)

names(motifHits.ByGene.RegElement) <- levels(Gene.ER.Enhancer.annotation.19TFs$gene)

motif.map = data.frame(GEMSTAT = sort(paste0("motif",c(0:18))), motifName = names(TF.motifs.Shrinked)[sort(paste0("motif",c(0:18)), index.return = T)$ix])
#################### list of matrces OF ANNOTATIONS############
motifHits.ByGene.RegElement
###############################################################
# now for each motif count its number of occurances in defined sets of genes: 


#########################################################################################################
#########################################################################################################
#######################################       Expanded TFs       ########################################
#######################################        Enrichment        ########################################
#########################################################################################################
#########################################################################################################
#using the set with 30 TFs and adding two other motifs for ESR1: hocomoco and half site 

#read the annotation output
Gene.ER.Enhancer.annotation.33TFs <- readLines("annotation_expanded.txt")
aaa <- which(Gene.ER.Enhancer.annotation.33TFs == "33 33 92 92")
Gene.ER.Enhancer.annotation.33TFs <- Gene.ER.Enhancer.annotation.33TFs[-c(1:aaa)]
Gene.ER.Enhancer.annotation.33TFs <- Gene.ER.Enhancer.annotation.33TFs[-c(1)]
aa <- grepl(pattern = "Annotated sites for", x = Gene.ER.Enhancer.annotation.33TFs)
Gene.ER.Enhancer.annotation.33TFs <- Gene.ER.Enhancer.annotation.33TFs[!aa]
aaa <- which(Gene.ER.Enhancer.annotation.33TFs == "Parameters for running the program: ")
Gene.ER.Enhancer.annotation.33TFs <- Gene.ER.Enhancer.annotation.33TFs[1:(aaa-1)]
aa <- lapply(X = Gene.ER.Enhancer.annotation.33TFs, FUN = strsplit, split = "\t")
aaa <- lapply(X = aa, FUN = unlist)
Gene.ER.Enhancer.annotation.33TFs <- do.call(rbind, aaa)

#aa <- as.integer(Gene.ER.Enhancer.annotation.33TFs[,2])
Gene.ER.Enhancer.annotation.33TFs <- data.frame(gene = Gene.ER.Enhancer.annotation.33TFs[,1],
                                                position = as.integer(Gene.ER.Enhancer.annotation.33TFs[,2]),
                                                motif = Gene.ER.Enhancer.annotation.33TFs[,3],
                                                LLRscore = as.numeric(Gene.ER.Enhancer.annotation.33TFs[,4]),
                                                maxLLR = as.numeric(Gene.ER.Enhancer.annotation.33TFs[,5]))
#Get enhancer number for each motif hit and add it to the annotation dataframe
aa = c(1, Gene.ER.Enhancer.annotation.33TFs$position)
aaa = c(Gene.ER.Enhancer.annotation.33TFs$position, 1)
a = aaa - aa
a = a[1:(length(a)-1)]

aa1 <- c(1, as.numeric(Gene.ER.Enhancer.annotation.33TFs$gene))
aaa1 <- c(as.numeric(Gene.ER.Enhancer.annotation.33TFs$gene), 1)
a1 <- (aa1 == aaa1)
a1[1] <- T
a1 <- a1[1:(length(a1) -1)]

a2 <- which(a < 0)
a2 <- c(1, a2)
a2 = c(a2, (length(a) + 1))
aa2 <- which(a1 == 0)
aa2 <- c(1, aa2)
aa3 <- numeric(length(a))
for(i in 1:(length(a2) - 1)){
  if (a2[i] %in% aa2){
    cnt <- 1
  }
  aa3[a2[i]:(a2[i+1] -1)]<- rep(cnt, length(aa3[a2[i]:(a2[i+1] -1)]))
  cnt <- cnt + 1
}

#################### DATA FRAME OF ANNOTATIONS############
Gene.ER.Enhancer.annotation.33TFs <- cbind(Gene.ER.Enhancer.annotation.33TFs, aa3)
colnames(Gene.ER.Enhancer.annotation.33TFs)[6] <- "reg.el.number"
###################
# Now get how many hits of each motif for each enhancer of each gene: 
# create a list: one entry per gene, the entry is a matrix: each row represents a
# motif and each column represents a regulatory element: the entries are the number of
# motif hits in that reg element
motifHits.ByGene.RegElement.expanded <- list()
number.of.motifs <- 33
#
for (i in 1:length(levels(Gene.ER.Enhancer.annotation.33TFs$gene))){
  # get the subset of dataset for this gene
  cur.annotation.frame <- subset(Gene.ER.Enhancer.annotation.33TFs, gene == levels(Gene.ER.Enhancer.annotation.33TFs$gene)[i], select = c( "motif", "reg.el.number"))
  # number of regulatory elements
  cur.reg.nu <- length(unique(cur.annotation.frame$reg.el.number))
  motifHits.ByGene.RegElement.expanded[[i]] <- matrix(0L, nrow=number.of.motifs, ncol=cur.reg.nu)
  for(reg.el in 1:cur.reg.nu){
    cur.table <- table(cur.annotation.frame$motif[cur.annotation.frame$reg.el.number == reg.el])
    motifHits.ByGene.RegElement.expanded[[i]][, reg.el] <- as.numeric(cur.table)
    rownames(motifHits.ByGene.RegElement.expanded[[i]]) <- names(cur.table)
  }
}
rm(i, cur.annotation.frame, cur.reg.nu, number.of.motifs, reg.el, cur.table)

names(motifHits.ByGene.RegElement.expanded) <- levels(Gene.ER.Enhancer.annotation.33TFs$gene)

motif.map.expanded = data.frame(GEMSTAT = sort(paste0("motif",c(0:32))), motifName = names(TF.motifs.Expanded)[sort(paste0("motif",c(0:32)), index.return = T)$ix])
#################### list of matrces OF ANNOTATIONS############
motifHits.ByGene.RegElement.expanded
###############################################################
# now for each motif count its number of occurances in defined sets of genes: 
aa <- Motif.Enrich.Illustrator(Gene.motif.matrix.list=motifHits.ByGene.RegElement.expanded,
                               LogFoldChangeExpression=DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.ER, 
                               maxLFCgreaterthan=2,
                               maxLFClessthan=1,
                               motifName = motif.map.expanded$motifName)


#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################      TF Expression        ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
TF.Index <- match(TF.ENTREZ$ENTREZID, rownames(DataSet.Experiment.Matrix.Batch.vgt1.logfold))
TF.Expression.LOGFOLD <- DataSet.Experiment.Matrix.Batch.vgt1.logfold[TF.Index,]
TF.Expression <- DataSet.Experiment.matrix.Batch.vgt1[TF.Index,]

TF.Index.Shrinked <- match(TF.ENTREZ.Shrinked$ENTREZID, rownames(DataSet.Experiment.Matrix.Batch.vgt1.logfold))
TF.Expression.LOGFOLD.Shrinked <- DataSet.Experiment.Matrix.Batch.vgt1.logfold[TF.Index.Shrinked, ]
TF.Expression.Shrinked <- DataSet.Experiment.matrix.Batch.vgt1[TF.Index.Shrinked, ]

aaa <- match(TF.ENTREZ$ENTREZID, rownames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile))
TF.Expression.Quantile <- DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile[aaa, ]

aaa <- match(TF.ENTREZ.Shrinked$ENTREZID, rownames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile))
TF.Expression.Quantile.Shrinked <- DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile[aaa, ]

par(mfrow = c(1,1), mar= rep(6,4))
rownames(TF.Expression) <- TF.ENTREZ$ALIAS
rownames(TF.Expression.Quantile) <- TF.ENTREZ$ALIAS

rownames(TF.Expression.Shrinked) <- TF.ENTREZ.Shrinked$ALIAS
rownames(TF.Expression.Quantile.Shrinked) <- TF.ENTREZ.Shrinked$ALIAS

boxplot.matrix(TF.Expression, use.cols = F, las =2 )
boxplot.matrix(TF.Expression.Quantile, use.cols = F, las =2 )

boxplot.matrix(TF.Expression, use.cols = T, las =2 )
#rownames(TF.Expression) <- TF.ENTREZ$ENTREZID
#####REMEMBER TO CHANGE ER EXPRESSION TO 0 AND 1 BASED ON THE CONDITION (E2 OR NOT)

#######Write to file for input to GEMSTAT
#The names and order have to be the same as the names for PWMs, 
#adding one more row for JUN because it has two motifs in motif file
TF.Expression.Quantile.For.GEMSTAT <- rbind(TF.Expression.Quantile[1:6, ], TF.Expression.Quantile[6, ], TF.Expression.Quantile[7:nrow(TF.Expression.Quantile), ])
rownames(TF.Expression.Quantile.For.GEMSTAT) <- c("AR", "CEBPB", "ESR1", "FOXA1", "GATA3", "JUN_1",
                                                  "JUN_2", "LEF1", "MAX", "MYC", "NFIB", "NKX3-1", "NR2C1", "NR2F2",
                                                  "NR3C1", "NR5A2", "PAX2", "PBX1", "PGR", "POU5F1", "PPARD", "PPARG", "RARA",
                                                  "RARG", "RELA", "RUNX1", "RXRA", "RXRB","SP1", "TFAP2C", "YBX1")
ExpressionWriter(TF.Expression.Quantile.For.GEMSTAT ,output.File.Name = "TF_Expression_quantile")
#################shrinked
TF.Expression.Quantile.Shrinked.For.GEMSTAT <- rbind(TF.Expression.Quantile.Shrinked[1:6, ], TF.Expression.Quantile.Shrinked[6, ], TF.Expression.Quantile.Shrinked[7:nrow(TF.Expression.Quantile.Shrinked), ])
rownames(TF.Expression.Quantile.Shrinked.For.GEMSTAT) <- c("AR", "CEBPB", "ESR1", "FOXA1", "GATA3", "JUN_1",
                                                  "JUN_2", "LEF1", "NKX3-1", 
                                                  "NR3C1", "NR5A2", "PBX1", "PGR", "RARA",
                                                  "RARG", "RUNX1", "RXRA","SP1", "TFAP2C")
ExpressionWriter(TF.Expression.Quantile.Shrinked.For.GEMSTAT ,output.File.Name = "TF_Expression_Shrinked_quantile")
#################Expanded
TF.Expression.Quantile.Expanded.For.GEMSTAT <- rbind(TF.Expression.Quantile[1:3, ],
                                                     TF.Expression.Quantile[3, ],
                                                     TF.Expression.Quantile[3:6, ],
                                                     TF.Expression.Quantile[6, ],
                                                     TF.Expression.Quantile[7:nrow(TF.Expression.Quantile), ])
rownames(TF.Expression.Quantile.Expanded.For.GEMSTAT) <- c("AR", "CEBPB", "ESR1_1", "ESR1_2", "ESR1_3","FOXA1", "GATA3", "JUN_1",
                                                           "JUN_2", "LEF1", "MAX", "MYC", "NFIB", "NKX3-1", "NR2C1", "NR2F2",
                                                           "NR3C1", "NR5A2", "PAX2", "PBX1", "PGR", "POU5F1", "PPARD", "PPARG", "RARA",
                                                           "RARG", "RELA", "RUNX1", "RXRA", "RXRB","SP1", "TFAP2C", "YBX1")
ExpressionWriter(TF.Expression.Quantile.Expanded.For.GEMSTAT ,output.File.Name = "TF_Expression_Expanded_quantile")


#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################         TF PWMs           ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################
#read in the PWM of the chosen TFs and write them all to a file for GEMSTAT input.
TF.motifs <- list()
aa <- list.files(path="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Motifs/Chosen", full.names=T, recursive=FALSE)
for( i in 1:length(aa)){
  TF.motifs[[i]] <- read.table(aa[i], header=T, sep ="\t")
  TF.motifs[[i]] =  TF.motifs[[i]][,2:5]
}

names(TF.motifs) <- c("AR", "CEBPB", "ESR1", "FOXA1", "GATA3", "JUN_1",
                      "JUN_2", "LEF1", "MAX", "MYC", "NFIB", "NKX3-1", "NR2C1", "NR2F2",
                      "NR3C1", "NR5A2", "PAX2", "PBX1", "PGR", "POU5F1", "PPARD", "PPARG", "RARA",
                      "RARG", "RELA", "RUNX1", "RXRA", "RXRB","SP1", "TFAP2C", "YBX1")
TF.motifs.count <- lapply(TF.motifs, PWMtoCount)
#########
#write the motifs to a file:
MotifWriter(motif.List = TF.motifs.count, output.File.Name = "motifs")
##########SHRINKED:
TF.motifs.Shrinked <- list()
aa <- list.files(path="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Motifs/Chosen_Shrinked", full.names=T, recursive=FALSE)
for( i in 1:length(aa)){
  TF.motifs.Shrinked[[i]] <- read.table(aa[i], header=T, sep ="\t")
  TF.motifs.Shrinked[[i]] =  TF.motifs.Shrinked[[i]][,2:5]
}

names(TF.motifs.Shrinked) <- c("AR", "CEBPB", "ESR1", "FOXA1", "GATA3", "JUN_1",
                               "JUN_2", "LEF1", "NKX3-1", 
                               "NR3C1", "NR5A2", "PBX1", "PGR", "RARA",
                               "RARG", "RUNX1", "RXRA","SP1", "TFAP2C")
TF.motifs.Shrinked.count <- lapply(TF.motifs.Shrinked, PWMtoCount)
#########
#write the motifs to a file:
MotifWriter(motif.List = TF.motifs.Shrinked.count, output.File.Name = "motifs_shrinked")
##########Expanded:
TF.motifs.Expanded <- list()
aa <- list.files(path="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Motifs/Chosen_expanded", full.names=T, recursive=FALSE)
for( i in 1:length(aa)){
  TF.motifs.Expanded[[i]] <- read.table(aa[i], header=T, sep ="\t")
  TF.motifs.Expanded[[i]] =  TF.motifs.Expanded[[i]][,2:5]
}

names(TF.motifs.Expanded) <- c("AR", "CEBPB", "ESR1_1", "ESR1_2", "ESR1_3","FOXA1", "GATA3", "JUN_1",
                      "JUN_2", "LEF1", "MAX", "MYC", "NFIB", "NKX3-1", "NR2C1", "NR2F2",
                      "NR3C1", "NR5A2", "PAX2", "PBX1", "PGR", "POU5F1", "PPARD", "PPARG", "RARA",
                      "RARG", "RELA", "RUNX1", "RXRA", "RXRB","SP1", "TFAP2C", "YBX1")
TF.motifs.Expanded.count <- lapply(TF.motifs.Expanded, PWMtoCount)
MotifWriter(motif.List = TF.motifs.Expanded.count, output.File.Name = "motifs_expanded")
##########Hocomoco:
TF.motifs.Shrinked.hocomoco <- list()
aa <- list.files(path="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Motifs/HocoMoco_Chosen/", full.names=T, recursive=FALSE)
for( i in 1:length(aa)){
  TF.motifs.Shrinked.hocomoco[[i]] <- read.table(aa[i], header=T, sep ="\t")
  TF.motifs.Shrinked.hocomoco[[i]] <-  TF.motifs.Shrinked.hocomoco[[i]][,2:5]
}

names(TF.motifs.Shrinked.hocomoco) <- c("AR", "CEBPB", "ESR1", "FOXA1", "GATA3", "JUN_1",
                                "LEF1", "NKX3-1", 
                               "NR3C1", "NR5A2", "PBX1", "PGR", "RARA",
                               "RARG", "RUNX1", "RXRA","SP1", "TFAP2C")
# fixing the monomer sites
# AR
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco[[1]]))
aa <- reverseComplement(as.matrix(TF.motifs.Shrinked.hocomoco[[1]] ))
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco[[1]] <- aa[1:6, ]
# ER
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco[[3]]))
aa <- reverseComplement(as.matrix(TF.motifs.Shrinked.hocomoco[[3]]))
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco[[3]] <- aa[4:9, ]
seqLogo::seqLogo(t(aa[4:9,]))
# GR
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco[[9]]))
aa <- reverseComplement(as.matrix(TF.motifs.Shrinked.hocomoco[[9]]))
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco[[9]] <- aa[1:6, ]
seqLogo::seqLogo(t(aa[1:6,]))
# RARA 13
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco[[13]]))
aa <- TF.motifs.Shrinked.hocomoco[[13]]
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco[[13]] <- aa[12:17, ]
seqLogo::seqLogo(t(aa[12:17,]))
# RARG 14
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco[[14]]))
aa <- TF.motifs.Shrinked.hocomoco[[14]]
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco[[14]] <- aa[2:7, ]
seqLogo::seqLogo(t(aa[2:7,]))
# RXRA 16
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco[[16]]))
aa <- TF.motifs.Shrinked.hocomoco[[16]]
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco[[16]] <- aa[8:13, ]
seqLogo::seqLogo(t(aa[8:13,]))

TF.motifs.Shrinked.hocomoco.count <- lapply(TF.motifs.Shrinked.hocomoco, PWMtoCount, remove.degen=F)
names(TF.motifs.Shrinked.hocomoco.count)
#########
#write the motifs to a file:
MotifWriter(motif.List = TF.motifs.Shrinked.hocomoco.count, output.File.Name = "motifs_shrinked_hocomoco", pseudo = 1)
# concatenate the dimer motifs to get the LLR that you want to use to filter the sum of LLR of monomers
TF.motifs.Shrinked.hocomoco.count_dimer_concat <- list()
for(i in c(1, 3, 9, 13, 14, 16)){
  TF.motifs.Shrinked.hocomoco.count_dimer_concat[[i]] <- rbind(TF.motifs.Shrinked.hocomoco.count[[i]], TF.motifs.Shrinked.hocomoco.count[[i]])
}
for(i in setdiff(c(1:18), c(1, 3, 9, 13, 14, 16))){
  TF.motifs.Shrinked.hocomoco.count_dimer_concat[[i]] <- TF.motifs.Shrinked.hocomoco.count[[i]]
}
names(TF.motifs.Shrinked.hocomoco.count_dimer_concat) <- names(TF.motifs.Shrinked.hocomoco.count)
# write motifs one by one for p-value analysis
for(i in 1:length(TF.motifs.Shrinked.hocomoco.count_dimer_concat)){
  MotifWriter(motif.List = TF.motifs.Shrinked.hocomoco.count_dimer_concat[i], 
              output.File.Name = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/MotifLLR2Pvalue/HocoMoco_Shrinked/",
                                        names(TF.motifs.Shrinked.hocomoco.count_dimer_concat)[i]),
              pseudo = 1)
}


# read the pvalue analysis results
aa <- list.files(path = "MotifLLR2Pvalue/HocoMoco_Shrinked/hocomoco_LLR_Pval_res/", pattern = "*.Filter", full.names = T)
TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval <- list()
for(i in 1:length(aa)){
  TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval[[i]] <- read.table(aa[i], header = F)
}
names(TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval) <- names(TF.motifs.Shrinked.hocomoco.count_dimer_concat)
TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_minLLR <- numeric(length(TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval))
for(i in 1:length(TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_minLLR)){
  aa <- which(TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval[[i]][,2] > 0.0004)[1]
  TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_minLLR[i] <- TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval[[i]][(aa-1), 1] / 1000
}
names(TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_minLLR) <- names(TF.motifs.Shrinked.hocomoco.count_dimer_concat)
TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_minLLR
#########################################################################################################
# adding E2F1 to the motifs

TF.motifs.Shrinked.hocomoco_E2F1_added <- list()
aa <- list.files(path="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Motifs/HocoMoco_Chosen_plus_E2F1/", full.names=T, recursive=FALSE)
for( i in 1:length(aa)){
  TF.motifs.Shrinked.hocomoco_E2F1_added[[i]] <- read.table(aa[i], header=T, sep ="\t")
  TF.motifs.Shrinked.hocomoco_E2F1_added[[i]] <-  TF.motifs.Shrinked.hocomoco_E2F1_added[[i]][,2:5]
}

names(TF.motifs.Shrinked.hocomoco_E2F1_added) <- c("AR", "CEBPB", "E2F1",  "ESR1", "FOXA1", "GATA3", "JUN_1",
                                        "LEF1", "NKX3-1", 
                                        "NR3C1", "NR5A2", "PBX1", "PGR", "RARA",
                                        "RARG", "RUNX1", "RXRA","SP1", "TFAP2C")

TF.motifs.Shrinked.hocomoco_E2F1_added_linear <- list()
aa <- list.files(path="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Motifs/HocoMoco_Chosen_plus_E2F1/", full.names=T, recursive=FALSE)
for( i in 1:length(aa)){
  TF.motifs.Shrinked.hocomoco_E2F1_added_linear[[i]] <- read.table(aa[i], header=T, sep ="\t")
  TF.motifs.Shrinked.hocomoco_E2F1_added_linear[[i]] <-  TF.motifs.Shrinked.hocomoco_E2F1_added_linear[[i]][,2:5]
}

names(TF.motifs.Shrinked.hocomoco_E2F1_added_linear) <- c("AR", "CEBPB", "E2F1",  "ESR1", "FOXA1", "GATA3", "JUN_1",
                                                   "LEF1", "NKX3-1", 
                                                   "NR3C1", "NR5A2", "PBX1", "PGR", "RARA",
                                                   "RARG", "RUNX1", "RXRA","SP1", "TFAP2C")
# fixing the monomer sites
# AR
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco_E2F1_added[[1]]))
aa <- reverseComplement(as.matrix(TF.motifs.Shrinked.hocomoco_E2F1_added[[1]] ))
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco_E2F1_added[[1]] <- aa[1:6, ]
# ER
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco_E2F1_added[[4]]))
aa <- reverseComplement(as.matrix(TF.motifs.Shrinked.hocomoco_E2F1_added[[4]]))
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco_E2F1_added[[4]] <- aa[4:9, ]
seqLogo::seqLogo(t(aa[4:9,]))
# GR
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco_E2F1_added[[10]]))
aa <- reverseComplement(as.matrix(TF.motifs.Shrinked.hocomoco_E2F1_added[[10]]))
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco_E2F1_added[[10]] <- aa[1:6, ]
seqLogo::seqLogo(t(aa[1:6,]))
# RARA 14
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco_E2F1_added[[14]]))
aa <- TF.motifs.Shrinked.hocomoco_E2F1_added[[14]]
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco_E2F1_added[[14]] <- aa[12:17, ]
seqLogo::seqLogo(t(aa[12:17,]))
# RARG 15
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco_E2F1_added[[15]]))
aa <- TF.motifs.Shrinked.hocomoco_E2F1_added[[15]]
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco_E2F1_added[[15]] <- aa[2:7, ]
seqLogo::seqLogo(t(aa[2:7,]))
# RXRA 17
seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco_E2F1_added[[17]]))
aa <- TF.motifs.Shrinked.hocomoco_E2F1_added[[17]]
seqLogo::seqLogo(t(aa))
TF.motifs.Shrinked.hocomoco_E2F1_added[[17]] <- aa[8:13, ]
seqLogo::seqLogo(t(aa[8:13,]))

TF.motifs.Shrinked.hocomoco_E2F1_added.count <- lapply(TF.motifs.Shrinked.hocomoco_E2F1_added, PWMtoCount, remove.degen=F)
names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
#########
#write the motifs to a file:
MotifWriter(motif.List = TF.motifs.Shrinked.hocomoco_E2F1_added.count, output.File.Name = "motifs_shrinked_hocomoco_E2F1_added", pseudo = 1)

for(i in 3:3){
  MotifWriter(motif.List = TF.motifs.Shrinked.hocomoco_E2F1_added.count[i], 
              output.File.Name = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/MotifLLR2Pvalue/HocoMoco_Shrinked/",
                                        names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[i]),
              pseudo = 1)
}

# read the pvalue analysis results
aa <- list.files(path = "MotifLLR2Pvalue/HocoMoco_Shrinked_E2F1_added/hocomoco_LLR_Pval_res_2", pattern = "*.Filter", full.names = T)
TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added <- list()
for(i in 1:length(aa)){
  TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added[[i]] <- read.table(aa[i], header = F)
}
names(TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added_minLLR <- numeric(length(TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added))
for(i in 1:length(TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added_minLLR)){
  aa <- which(TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added[[i]][,2] > 0.001)[1]
  TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added_minLLR[i] <- TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added[[i]][(aa-1), 1] / 1000
}
names(TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added_minLLR) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added_minLLR

aa <- cbind(c(0.001, 0.01, 0.01, 0.001, 0.01, 0.01, 0.01,  0.01,  0.01,  0.001, 0.01,  0.01,  0.01, 0.001, 0.001, 0.01, 0.001,  0.01, 0.01),
            c(0.01,  100,  100 , 0.01 , 100,  100,  100,   100,   100,   0.01,  100,   100,   100,  0.01 , 0.01,  100,  0.01,   100,  100))
rownames(aa) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
aa <- c(0.005, 1,   1,   0.005,  1,     1,     1,     1,     1,     0.005, 1,     1,     1,      0.005,  0.005, 1     , 0.005,  1,     1)
names(aa) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)

#########################################################################################################
#########################################################################################################
######################################                           ########################################
######################################          Sequence         ########################################
######################################                           ########################################
#########################################################################################################
#########################################################################################################

#write the sequence of the regulatory elements for all genes of interst to a file
aa  <- which(rownames(DataSet.Experiment.matrix.Batch.vgt1) %in% Genes.Associated.REMAP.ER.Entrez)
aaa <- match(rownames(DataSet.Experiment.matrix.Batch.vgt1)[aa], names(Vicinity100kb.Enhancer.plusPromoter.By.gene))
WriteFastaOfBag(enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene[aaa],
                Enhancer.Sequence=MCFEnhancers[,"sequ"],
                promoter.promoter.int=Promoter.gene.PromoterIntList[aaa],
                promoter.sequence=Promoter.Sequence.Char,
                output.File.Name = "RegulatoryElements" )
#describing length of enhancers
aa = nchar(MCFEnhancers[,"sequ"])
boxplot(aa, outline=F, ylab = "Length of enhancer (bp)", main = "MCF7 enhancers")
hist(aa, breaks = 100, xlab = "length of enhancers", main = "MCF7 enhancers")
a  <- which(rownames(DataSet.Experiment.matrix.Batch.vgt1) %in% Genes.Associated.REMAP.ER.Entrez)
aaa <- match(rownames(DataSet.Experiment.matrix.Batch.vgt1)[a], names(Vicinity100kb.Enhancer.plusPromoter.By.gene))
aaa <- unique(unlist(Vicinity100kb.Enhancer.plusPromoter.By.gene[aaa]))
hist(aa[aaa], breaks = 100, xlab = "length of enhancers", main = "Enhancers of interest")
boxplot(aa[aaa], outline=F, ylab = "Length of enhancer (bp)", main = "Enhancers of interest")

####################################################################################################################################
####################################################################################################################################
############################################       Design        ###################################################################
############################################       GEMSTAT       ###################################################################
############################################     Experiments     ###################################################################
####################################################################################################################################
####################################################################################################################################
# In this part I will generate input data for GEMSTAT to test a hypothesis
#for each experiment you need to provide at least four inputs: Sequence, Motif, TF_expression, Gene expression
################################################ start Experiment #1 ############################################################
########Choose Genes and clusters########

aaa <- apply(X = DataSet.Experiment.matrix.Batch.vgt1,1, mean)
#remove genes with average exp less than one
DataSet.Experiment.matrix.Batch.vgt1.mgt20 <- DataSet.Experiment.matrix.Batch.vgt1[aaa > 20,]
aaaa = which(rownames(DataSet.Experiment.matrix.Batch.vgt1.mgt20) %in% Genes.Associated.REMAP.ER.Entrez)
DataSet.Experiment.matrix.Batch.vgt1.mgt20.ER <- DataSet.Experiment.matrix.Batch.vgt1.mgt20[aaaa, ]
DataSet.Experiment.matrix.Batch.vgt1.mgt20.LFC.ER <- LogFoldTransform(DataSet.Experiment.matrix.Batch.vgt1.mgt20.ER, metadataM = DataSet.Experiment.Metadata)
#make a version which genes with max expression less than 20 are excluded
aaa <- apply(X = DataSet.Experiment.matrix.Batch.vgt1,1, max)
DataSet.Experiment.matrix.Batch.vgt1.maxgt20        <- DataSet.Experiment.matrix.Batch.vgt1[aaa > 20,]
#add the row containing the "5460" TF
DataSet.Experiment.matrix.Batch.vgt1.maxgt20 <- rbind(DataSet.Experiment.matrix.Batch.vgt1.maxgt20, DataSet.Experiment.matrix.Batch.vgt1[which(rownames(DataSet.Experiment.matrix.Batch.vgt1) == "5460"), ])
rownames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20)[nrow(DataSet.Experiment.matrix.Batch.vgt1.maxgt20)] <- "5460"
aaaa <- which(rownames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20) %in% Genes.Associated.REMAP.ER.Entrez)
DataSet.Experiment.matrix.Batch.vgt1.maxgt20.ER     <- DataSet.Experiment.matrix.Batch.vgt1.maxgt20[aaaa, ]
DataSet.Experiment.matrix.Batch.vgt1.maxgt20.LFC.ER <- LogFoldTransform(DataSet.Experiment.matrix.Batch.vgt1.maxgt20.ER,
                                                                        metadataM = DataSet.Experiment.Metadata)


#quantile normalize
DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile   <- normalize.quantiles(DataSet.Experiment.matrix.Batch.vgt1.mgt20)
DataSet.Experiment.matrix.Batch.vgt1.quantile         <- normalize.quantiles(DataSet.Experiment.matrix.Batch.vgt1)
DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile <- normalize.quantiles(DataSet.Experiment.matrix.Batch.vgt1.maxgt20)

rownames(DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile)   <- rownames(DataSet.Experiment.matrix.Batch.vgt1.mgt20)
rownames(DataSet.Experiment.matrix.Batch.vgt1.quantile)         <- rownames(DataSet.Experiment.matrix.Batch.vgt1)
rownames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile) <- rownames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20)
colnames(DataSet.Experiment.matrix.Batch.vgt1.quantile) <- colnames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20)
colnames(DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile) <- colnames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20)
colnames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile) <- colnames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20)

aaaa = which(rownames(DataSet.Experiment.matrix.Batch.vgt1.mgt20) %in% Genes.Associated.REMAP.ER.Entrez)
DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile.ER     <- DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile[aaaa, ]
DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile.LFC.ER <- LogFoldTransform(DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile.ER,
                                                                               metadataM = DataSet.Experiment.Metadata)
aaaa = which(rownames(DataSet.Experiment.matrix.Batch.vgt1) %in% Genes.Associated.REMAP.ER.Entrez)
DataSet.Experiment.matrix.Batch.vgt1.quantile.ER     <- DataSet.Experiment.matrix.Batch.vgt1.quantile[aaaa, ]
DataSet.Experiment.matrix.Batch.vgt1.quantile.LFC.ER <- LogFoldTransform(DataSet.Experiment.matrix.Batch.vgt1.quantile.ER,
                                                                         metadataM = DataSet.Experiment.Metadata)
aaaa = which(rownames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20) %in% Genes.Associated.REMAP.ER.Entrez)
DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.ER     <- DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile[aaaa, ]
DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.ER <- LogFoldTransform(DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.ER,
                                                                                 metadataM = DataSet.Experiment.Metadata)

#remove genes with max absolute fold change less than one
aa <- apply(abs(DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile.LFC.ER), 1, max)
DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile.LFC.gt1.ER <- DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile.LFC.ER[aa > 1, ]
DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile.gt1.ER     <- DataSet.Experiment.matrix.Batch.vgt1.mgt20.quantile.ER[aa > 1, ]

aa <- apply(abs(DataSet.Experiment.matrix.Batch.vgt1.quantile.LFC.ER), 1, max)
DataSet.Experiment.matrix.Batch.vgt1.quantile.LFC.gt1.ER <- DataSet.Experiment.matrix.Batch.vgt1.quantile.LFC.ER[aa > 1, ]
DataSet.Experiment.matrix.Batch.vgt1.quantile.gt1.ER     <- DataSet.Experiment.matrix.Batch.vgt1.quantile.ER[aa > 1, ]

aa <- apply(abs(DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.ER), 1, max)
DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.gt1.ER <- DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.ER[aa > 1, ]
DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.gt1.ER     <- DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.ER[aa > 1, ]

###plot multiple clusterings and next choose the clustering for experiment 1
Experiment1.Cluster.List <- list()
for(i in 4:15){
  Experiment1.Cluster.List[[i]] <- clustratorWrapper(.ExpressionMat=DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.gt1.ER,
                          .nfolds=i,
                          .exportPlot = T,
                          harshThreshTimes = 10,
                          RepeatMax = 25,
                          dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Expression/Experiment_1_Clustering_quantile_2334/")
}

###
#Chose clustering number 
Experiment1.Chosen.Cluster <- Experiment1.Cluster.List[[5]][[3]]
#create all input data for this experiment
#Expression Matrix: DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.gt1.ER
GEMSTAT_input_constructor(clusterIndex = Experiment1.Chosen.Cluster$ClusterIndex,
                          gene_expression_Mat = DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.gt1.ER,
                          dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/GEMSTAT/Experiment_1/Inputs",
                          experiment.nu=1)


################################################  end Experiment #1  ############################################################
################################################ start Experiment #2 ############################################################
#same clustering as experiment 1 but motifs are shrinked to 19 instead of 30
GEMSTAT_input_constructor(clusterIndex = Experiment1.Chosen.Cluster$ClusterIndex,
                          gene_expression_Mat = DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.gt1.ER,
                          experiment.nu=2, .motifs = TF.motifs.Shrinked.count,.TF_expression =TF.Expression.Quantile.Shrinked.For.GEMSTAT)
################################################  end Experiment #2  ############################################################
################################################ start Experiment #3 ############################################################
#same clustering as experiment 1 but motifs are shrinked to 19 instead of 30, and enhancers have to be shorter than 10kb
GEMSTAT_input_constructor(clusterIndex = Experiment1.Chosen.Cluster$ClusterIndex,
                          gene_expression_Mat = DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.gt1.ER,
                          experiment.nu=3,
                          .motifs = TF.motifs.Shrinked.count,
                          .TF_expression =TF.Expression.Quantile.Shrinked.For.GEMSTAT,
                          .max.enh.length = 10000)
################################################  end Experiment #3  ############################################################
################################################ start Experiment #4 ############################################################
Experiment4.Cluster.List <- list()
for(i in 41:60){
  Experiment4.Cluster.List[[i]] <- clustratorWrapper(.ExpressionMat=DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.gt1.ER,
                                                     .nfolds=i,
                                                     .exportPlot = T,
                                                     harshThreshTimes = 4,
                                                     RepeatMax = 20,
                                                     dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Expression/Experiment_4_Clustering_quantile_2334/")
}
################################################  end Experiment #4  ############################################################
################################################ start Experiment #5 ############################################################
#preprocess log fold change data before clustering: 
#set max fold change of each dataset to 1
aa <- DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.ER
a = numeric(length(unique(DataSet.Experiment.Metadata$DatasetNo)))
for(i in 1:length(a)){
  a[i] <- max(abs(aa[,(DataSet.Experiment.Metadata$DatasetNo == i)]))
  aa[,(DataSet.Experiment.Metadata$DatasetNo == i)] <- aa[,(DataSet.Experiment.Metadata$DatasetNo == i)] / a[i]
}
#remove controls, 
aa <- aa[,DataSet.Experiment.Metadata$expNo != 1]
aa.meta <- DataSet.Experiment.Metadata[DataSet.Experiment.Metadata$expNo != 1,]
DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.ER.noControl.DSnorm <- aa
DataSet.Experiment.Metadata.nocontrol <- aa.meta
####perform clustering:
Experiment5.Cluster.List <- list()
for(i in 4:40){
  Experiment5.Cluster.List[[i]] <- clustratorWrapper(.ExpressionMat=DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.ER.noControl.DSnorm,
                                                     .Exp.metadata = DataSet.Experiment.Metadata.nocontrol,
                                                     .no.control = T,
                                                     .nfolds=i,
                                                     .exportPlot = T,
                                                     harshThreshTimes = 10,
                                                     RepeatMax = 20,
                                                     dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/plots/Expression/Experiment_5_Clustering_quantile_4529/")
}
colnames(DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.ER) <- DataSet.Experiment.matrix.Names
#write csv of exp mat for zenvisage

aa <- DataSet.Experiment.matrix.Batch.vgt1.maxgt20.quantile.LFC.ER.noControl.DSnorm
#aaa<- data.frame(Gene = rep(rownames(aa), ncol(aa)), condition, expression)
aaa <- cbind(rownames(aa), aa[,1], rep(colnames(aa)[1], nrow(aa))) 
for(i in 2:ncol(aa)){
  aaa = rbind(aaa, cbind(rownames(aa), aa[,i], rep(colnames(aa)[i], nrow(aa))) )
}
colnames(aaa) <- c("gene", "expression", "condition")
write.table(x = aaa,
            file = "LFC_DSNORM_noControl_4529.csv",
            quote = F,
            sep = ",",
            row.names = F,
            col.names = T)
####################################################################################################################################
####################################################################################################################################
############################################       Analyze        ##################################################################
############################################       GEMSTAT        ##################################################################
############################################       Results        ##################################################################
####################################################################################################################################
####################################################################################################################################
par(mfrow = c(1, 2), mar = c(4,4,4,4))
plot(GSE78167.RNAseq.Avg.Norm[rownames(GSE78167.RNAseq.Avg.Norm) %in% "2099", ], main = "ESR1 expression", xlab = "time points", ylab = "", type = "l")
plot(c(0, rep(1, 9)), main = "ESR1 expression", xlab = "time points", ylab = "", type = "l")

plot(GSE78167.RNAseq.Avg.Norm[rownames(GSE78167.RNAseq.Avg.Norm) %in% "5914", ], main = "RARA expression", xlab = "time points", ylab = "", type = "l")


#########################################################################################################
#########################################################################################################
#########################################################################################################
######################################        Create the +1/0/-1 data set     ###########################
######################################          of purterbation studies       ###########################
#########################################################################################################
#########################################################################################################
#########################################################################################################

# In this part I want to create a new dataset. I will use purterbation studies which include E2 treatment.
# I will construct a matrix where rows are genes and columns are conditions.
# conditions represent knock downs or different treatments, after which E2 treatment is applied 
# Values indicate the comparison of a gene's expression in that condition before and after E2 treatment.
# What are the studies that I will consider?
# for each study 
# 1.download the expression data, specify the address
# 2.specify the conditions  I am interested in
# 3.process everything to assign -1, 0, or 1 to each gene
#########################################################################################################
#########################################################################################################
#########################################################################################################

# 0 #1. GSE60271 : dropped for now: hard to process :GRO-Seq : WT+E2/WT, (siGATA3 + E2)/(siGATA3 + Veh), (shAP2g + E2)/(shAP2g + Veh), (shRARs + E2)/(shRARs + Veh) : hg18

#########################################################################################################
#########################################################################################################
#########################################################################################################

# 0 #2. GSE68358 : dropped because is doesn't have vehicle treatment for progestin conditions. previously ddescribed : progestin

#########################################################################################################
#########################################################################################################
#########################################################################################################

# 1 #3. GSE55923 : previously ddescribed : BRD4/JQ1
GSE55922.RNAseq
colnames(GSE55922.RNAseq)

#########################################################################################################
#########################################################################################################
#########################################################################################################

# 0 #4. GSE67295 : dropped for now: doesn't have the desired conditions : previously described ICI, IL1b, ...

#########################################################################################################
#########################################################################################################
#########################################################################################################

# 0 #5. GSE80098 : dropped for now: doesn't have the desired conditions : previously described R5020

#########################################################################################################
#########################################################################################################
#########################################################################################################

# 2 #6. GSE37386: GREB1 KD : downloaded : 
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE37386_non-normalized_data.txt

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE37386.exp.soft <- getGEO("GSE37386", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE37386.exp.soft <- GSE37386.exp.soft$GSE37386_series_matrix.txt.gz
# convert annotation
GSE37386.exp.soft <- as(GSE37386.exp.soft, "ExpressionSetIllumina")
# mark regular vs negatie control probes
fData(GSE37386.exp.soft)$Status <- ifelse(fData(GSE37386.exp.soft)$PROBEQUALITY=="No match","negative","regular" )
## This is the normalization function which can be done using different methods, my data is already log2 transformed and quantile normalized
# summaryData.norm <- normaliseIllumina(summaryData,method="quantile",
#                                       status=fData(summaryData)$Status, transform = "log2")

# FIlTERATION of bad data
# We recommend removing probes assigned a ‘Bad’ or ‘No match’ quality score after normalization

#getting the quality of probes from the corresponding annotation database. in this case illuminaHumanv4
aa <- unlist(mget(as.character(featureNames(GSE37386.exp.soft)), illuminaHumanv4PROBEQUALITY, ifnotfound=NA))
aa[3574] <- "Good" #changing to keep NKX3-1 in the dataset
table(aa)

#remove probes with bad , nomatch or NA quality
aaa <- aa == "No match" | aa == "Bad" | is.na(aa)
GSE37386.exp.soft.filt <- GSE37386.exp.soft[!aaa,]
dim(GSE37386.exp.soft.filt)
GSE37386.exp.soft.filt <- GSE37386.exp.soft.filt[!is.na(fData(GSE37386.exp.soft.filt)$Entrez_Gene_ID), ]

aa <- apply(X = exprs(GSE37386.exp.soft.filt), MARGIN = 1, FUN = var)
aa[28115] <- 0.1 #in order to keep GATA3 in the dataset
aa[2588]
GSE37386.exp.soft.filt <- GSE37386.exp.soft.filt[aa > quantile(aa)[2], ]
# Differential expression calculations
GSE37386.exp.soft.filt.Descrete <- createDescreteMat(columnNames = c("WT_6h", "WT_24h", "GREB1KD_6h", "GREB1KD_24h"),
                                                     conditions = paste(pData(GSE37386.exp.soft.filt)[,42],pData(GSE37386.exp.soft.filt)[,43], pData(GSE37386.exp.soft.filt)[,44] ,sep = "_"),
                                                     contrasts = cbind(c("siNT_control_vehicle", "siNT_control_vehicle", "siGreb_control_vehicle", "siGreb_control_vehicle"),
                                                                       c("siNT_6hrs_E2", "siNT_24hrs_E2", "siGreb_6hrs_E2", "siGreb_24hrs_E2")),
                                                     expMat=exprs(GSE37386.exp.soft.filt),
                                                     adjpvalThresh = 1,
                                                     lfcThresh = 0,
                                                     expfData = fData(GSE37386.exp.soft.filt),
                                                     exppData = pData(GSE37386.exp.soft.filt)
                                                     )
#adjpval 1e-6:
# [1] "number of diff exp probes in contrast  siNT_6hrs_E2 - siNT_control_vehicle  :  105"
# [1] "number of diff exp probes in contrast  siNT_24hrs_E2 - siNT_control_vehicle  :  537"
# [1] "number of diff exp probes in contrast  siGreb_6hrs_E2 - siGreb_control_vehicle  :  37"
# [1] "number of diff exp probes in contrast  siGreb_24hrs_E2 - siGreb_control_vehicle  :  291"

# lfc > 2:
# [1] "number of diff exp probes in contrast  siNT_6hrs_E2 - siNT_control_vehicle  :  3"
# [1] "number of diff exp probes in contrast  siNT_24hrs_E2 - siNT_control_vehicle  :  15"
# [1] "number of diff exp probes in contrast  siGreb_6hrs_E2 - siGreb_control_vehicle  :  1"
# [1] "number of diff exp probes in contrast  siGreb_24hrs_E2 - siGreb_control_vehicle  :  6"

# lfc > 1:
# [1] "number of diff exp probes in contrast  siNT_6hrs_E2 - siNT_control_vehicle  :  57"
# [1] "number of diff exp probes in contrast  siNT_24hrs_E2 - siNT_control_vehicle  :  217"
# [1] "number of diff exp probes in contrast  siGreb_6hrs_E2 - siGreb_control_vehicle  :  24"
# [1] "number of diff exp probes in contrast  siGreb_24hrs_E2 - siGreb_control_vehicle  :  91"

#TF expression
match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE37386.exp.soft.filt)$Entrez_Gene_ID)
nrow(GSE37386.exp.soft.filt.Descrete$raw_lfc)
ncol(GSE37386.exp.soft.filt.Descrete$raw_lfc)
hist(GSE37386.exp.soft.filt.Descrete$raw_lfc)
par(mfrow = c(1, 1), mar = c(4,4,4,4))
#########################################################################################################
#########################################################################################################
#########################################################################################################

# S #7. GSE25315 : tamoxifen : FOXA1 KD : 2 parts : downloaded
# 3 # part 1 : tamoxifen : GSE25314 : /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE25314_normalized.txt

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE25314.exp.soft <- getGEO("GSE25314", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE25314.exp.soft <- GSE25314.exp.soft$GSE25314_series_matrix.txt.gz
# convert annotation
GSE25314.exp.soft <- as(GSE25314.exp.soft, "ExpressionSetIllumina")
# mark regular vs negatie control probes
fData(GSE25314.exp.soft)$Status <- ifelse(fData(GSE25314.exp.soft)$PROBEQUALITY=="No match","negative","regular" )

match("1869", fData(GSE25314.exp.soft)$Entrez_Gene_ID)
#getting the quality of probes from the corresponding annotation database. in this case illuminaHumanv4
aa <- unlist(mget(as.character(featureNames(GSE25314.exp.soft)), illuminaHumanv3PROBEQUALITY, ifnotfound=NA))
table(aa)
aa[3574] <- "Good" # changing in order to keep NKX3-1
aa[41772]
#remove probes with bad , nomatch or NA quality
aaa <- aa == "No match" | aa == "Bad" | is.na(aa)
GSE25314.exp.soft.filt <- GSE25314.exp.soft[!aaa,]
dim(GSE25314.exp.soft.filt)

match("1869", fData(GSE25314.exp.soft.filt)$Entrez_Gene_ID)
#filter out probes that do not have an ENTREZ ID associated with them and the probes with very low (??) variance, and 
aa <- apply(X = exprs(GSE25314.exp.soft.filt), MARGIN = 1, FUN = var)
aa[6411] <- 0.1 #changing in order to keep CEBPB in the dataset
aa[28285] <- 0.1 #changing in order to keep E2F1 in the dataset
GSE25314.exp.soft.filt <- GSE25314.exp.soft.filt[aa > quantile(aa)[2], ]
GSE25314.exp.soft.filt <- GSE25314.exp.soft.filt[!is.na(fData(GSE25314.exp.soft.filt)$Entrez_Gene_ID), ]
dim(GSE25314.exp.soft.filt)

sum(is.na(fData(GSE25314.exp.soft.filt)$Entrez_Gene_ID))

aaRepicates <- unlist(lapply(strsplit(as.character(pData(GSE25314.exp.soft)$description), split = " "), "[[", 2))
aaCondition <- gsub(" ", "_", as.character(pData(GSE25314.exp.soft)$'treatment:ch1'))
aaCondition <- gsub("\\(", "", aaCondition)
aaCondition <- gsub("\\)", "", aaCondition)
aaColumn <- c("WT_6h", "Tamoxifen_6h")
aaContrasts <- cbind(c("vehicle_control", "tamoxifen"),
                     c("estrogen", "estrogen_and_tamoxifen"))

# Differential expression calculations
GSE25314.exp.soft.filt.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                     conditions = aaCondition,
                                                     replicates = aaRepicates,
                                                     contrasts = aaContrasts,
                                                     expMat=exprs(GSE25314.exp.soft.filt),
                                                     adjpvalThresh = 1,
                                                     lfcThresh = 0,
                                                   expfData = fData(GSE25314.exp.soft.filt),
                                                     exppData = pData(GSE25314.exp.soft.filt) )

# padj 1e-6:
# [1] "number of diff exp probes in contrast  estrogen - vehicle_control  :  156"
# [1] "number of diff exp probes in contrast  estrogen_and_tamoxifen - tamoxifen  :  0"

# lfc > 2:
# [1] "number of diff exp probes in contrast  estrogen - vehicle_control  :  3"
# [1] "number of diff exp probes in contrast  estrogen_and_tamoxifen - tamoxifen  :  0"

# lfc > 1:
# [1] "number of diff exp probes in contrast  estrogen - vehicle_control  :  93"
# [1] "number of diff exp probes in contrast  estrogen_and_tamoxifen - tamoxifen  :  0"
#TF expression
match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE25314.exp.soft.filt)$Entrez_Gene_ID)
match(TF.ENTREZ.Shrinked_2$ENTREZID, fData(GSE25314.exp.soft.filt)$Entrez_Gene_ID)

#########################################################################################################

# 4 # part 2 : FOXA1 KD  : GSE25315 : /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE25315_normalized.txt

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE25315.exp.soft <- getGEO("GSE25315", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE25315.exp.soft <- GSE25315.exp.soft$GSE25315_series_matrix.txt.gz
# convert annotation
GSE25315.exp.soft <- as(GSE25315.exp.soft, "ExpressionSetIllumina")
# mark regular vs negatie control probes
fData(GSE25315.exp.soft)$Status <- ifelse(fData(GSE25315.exp.soft)$PROBEQUALITY=="No match","negative","regular" )

#getting the quality of probes from the corresponding annotation database. in this case illuminaHumanv4
aa <- unlist(mget(as.character(featureNames(GSE25315.exp.soft)), illuminaHumanv4PROBEQUALITY, ifnotfound=NA))
table(aa)
aa[3574] <- "Good"  # changing in order to keep NKX3-1 in the dataset
#remove probes with bad , nomatch or NA quality
aaa <- aa == "No match" | aa == "Bad" | is.na(aa)
GSE25315.exp.soft.filt <- GSE25315.exp.soft[!aaa,]
dim(GSE25315.exp.soft.filt)

aa <- apply(X = exprs(GSE25315.exp.soft.filt), MARGIN = 1, FUN = var)
aa[12289] <- 0.1 #changing in order to keep RARG in dataset
GSE25315.exp.soft.filt <- GSE25315.exp.soft.filt[aa > quantile(aa)[2], ]
GSE25315.exp.soft.filt <- GSE25315.exp.soft.filt[!is.na(fData(GSE25315.exp.soft.filt)$Entrez_Gene_ID), ]
dim(GSE25315.exp.soft.filt)

aaRepicates <- unlist(lapply(strsplit(as.character(pData(GSE25315.exp.soft)$description), split = " "), "[[", 2))
aaCondition <- gsub(" ", "_", as.character(pData(GSE25315.exp.soft)$'treatment:ch1'))
aaCondition <- gsub("\\(", "", aaCondition)
aaCondition <- gsub("\\)", "", aaCondition)
aaColumn <- c("WT_6h", "FoxA1_siRNA_6h")
aaContrasts <- cbind(c("control_siRNA_and_vehicle_control", "FoxA1_siRNA_and_vehicle_control"),c("control_siRNA_and_estrogen", "FoxA1_siRNA_and_estrogen"))

# Differential expression calculations
GSE25315.exp.soft.filt.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                     conditions = aaCondition,
                                                     replicates = aaRepicates,
                                                     contrasts = aaContrasts,
                                                     expMat=exprs(GSE25315.exp.soft.filt),
                                                     adjpvalThresh = 1,
                                                     lfcThresh = 0,
                                                     expfData = fData(GSE25315.exp.soft.filt),
                                                     exppData = pData(GSE25315.exp.soft.filt))
# adjpval 1e-6:
# [1] "number of diff exp probes in contrast  control_siRNA_and_estrogen - control_siRNA_and_vehicle_control  :  0"
# [1] "number of diff exp probes in contrast  FoxA1_siRNA_and_estrogen - FoxA1_siRNA_and_vehicle_control  :  0"

# lfc > 2:
# [1] "number of diff exp probes in contrast  control_siRNA_and_estrogen - control_siRNA_and_vehicle_control  :  0"
# [1] "number of diff exp probes in contrast  FoxA1_siRNA_and_estrogen - FoxA1_siRNA_and_vehicle_control  :  0"

# lfc > 1:
# [1] "number of diff exp probes in contrast  control_siRNA_and_estrogen - control_siRNA_and_vehicle_control  :  1"
# [1] "number of diff exp probes in contrast  FoxA1_siRNA_and_estrogen - FoxA1_siRNA_and_vehicle_control  :  5"
match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE25315.exp.soft.filt)$Entrez_Gene_ID)

colnames(GSE25315.exp.soft.filt.Descrete$raw_lfc)
#########################################################################################################
#########################################################################################################
#########################################################################################################

# 5 #8. GSE102367 : AP1/c-Jun overexpression : downloaded : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE102367_family.soft
# the zipped file is also there, delete the soft file after reading
biocLite("hta20transcriptcluster.db")
library(pd.hta.2.0)
library(hta20transcriptcluster.db)
library(pd.clariom.d.human)


setwd("Expression_data/Descretized_Expression_dataset/GSE102367_RAW")
library(oligo)
library(affycoretools)
library(clariomdhumantranscriptcluster.db)
aadat <- read.celfiles(list.celfiles())
aaeset <- rma(aadat)
aaeset3 <- annotateEset(aaeset, clariomdhumantranscriptcluster.db)
aaeset3 <- aaeset3[!is.na(fData(aaeset3)$ENTREZID), ]
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
pData(aaeset3)
GSE102367.exp <- aaeset3

aa <- apply(GSE102367.exp, MARGIN = 1, var)
aa[match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE102367.exp)$ENTREZID)] <- 0.1
quantile(aa)
GSE102367.exp <- GSE102367.exp[aa > quantile(aa)[2], ]

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE102367.exp.soft <- getGEO("GSE102367", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE102367.exp.soft <- GSE102367.exp.soft$GSE102367_series_matrix.txt.gz

#transfering pDATA to the previous matrix
pData(GSE102367.exp) <- pData(GSE102367.exp.soft)
remove(GSE102367.exp.soft)
aaOVer <- rep(c("Normal", "Normal", "Normal", "AP1_OverExp", "AP1_OverExp", "AP1_OverExp"), 3)
aatreat <- c(rep("vehicle", 6), rep("E2", 6), rep("Tamoxifen", 6))

aaRepicates <- rep(c(1,2,3), 6)
aaCondition <- paste(aaOVer, aatreat, sep = "_")

aaColumn <- c("WT_24h", "AP1_OverExp_24h")
aaContrasts <- cbind(c("Normal_vehicle", "AP1_OverExp_vehicle"),c("Normal_E2", "AP1_OverExp_E2"))

# Differential expression calculations
GSE102367.exp.Descrete <-          createDescreteMat(columnNames = aaColumn,
                                                     conditions = aaCondition,
                                                     replicates = aaRepicates,
                                                     contrasts = aaContrasts,
                                                     expMat=exprs(GSE102367.exp),
                                                     adjpvalThresh = 1,
                                                     lfcThresh = 0,
                                                     expfData = fData(GSE102367.exp),
                                                     exppData = pData(GSE102367.exp))
# padj 1e-6:
# [1] "number of diff exp probes in contrast  Normal_E2 - Normal_vehicle  :  158"
# [1] "number of diff exp probes in contrast  AP1_OverExp_E2 - AP1_OverExp_vehicle  :  179"

# lfc > 2:
# [1] "number of diff exp probes in contrast  Normal_E2 - Normal_vehicle  :  43"
# [1] "number of diff exp probes in contrast  AP1_OverExp_E2 - AP1_OverExp_vehicle  :  47"

# lfc > 1:
# [1] "number of diff exp probes in contrast  Normal_E2 - Normal_vehicle  :  547"
# [1] "number of diff exp probes in contrast  AP1_OverExp_E2 - AP1_OverExp_vehicle  :  499"

match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE102367.exp)$ENTREZID)



#########################################################################################################
#########################################################################################################
#########################################################################################################

# 6 #9. GSE79761 : DEX+E2 or DEX : downloaded : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE79761_family.soft
# the zipped file is also there, delete the soft file after reading

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE79761.exp.soft <- getGEO("GSE79761", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE79761.exp.soft <- GSE79761.exp.soft$GSE79761_series_matrix.txt.gz
# # convert annotation
# GSE79761.exp.soft <- as(GSE79761.exp.soft, "ExpressionSetIllumina")
# # mark regular vs negatie control probes
# fData(GSE79761.exp.soft)$Status <- ifelse(fData(GSE79761.exp.soft)$PROBEQUALITY=="No match","negative","regular" )

# #getting the quality of probes from the corresponding annotation database. in this case illuminaHumanv4
# aa <- unlist(mget(as.character(featureNames(GSE79761.exp.soft)), GPL570PROBEQUALITY, ifnotfound=NA))
# table(aa)

# #remove probes with bad , nomatch or NA quality
# aaa <- aa == "No match" | aa == "Bad" | is.na(aa)
# GSE25315.exp.soft.filt <- GSE25315.exp.soft[!aaa,]
dim(GSE79761.exp.soft)
annotation(GSE79761.exp.soft)
aa <- annotateEset(GSE79761.exp.soft, hgu133plus2.db)

#GSE79761.exp.soft <- GSE79761.exp.soft[!(fData(GSE79761.exp.soft)$ENTREZ_GENE_ID == ""), ]
GSE79761.exp.soft <- annotateEset(GSE79761.exp.soft, hgu133plus2.db)
GSE79761.exp.soft <- GSE79761.exp.soft[!is.na(fData(GSE79761.exp.soft)$ENTREZID), ]

aa <- apply(X = exprs(GSE79761.exp.soft), MARGIN = 1, FUN = var)
quantile(aa)
GSE79761.exp.soft <- GSE79761.exp.soft[aa > quantile(aa)[2], ]

aaRepicates <- rep(c(1,2,3), 4)
aaCondition <- c(rep("Vehicle", 3), rep("E2_4h", 3), rep("Dex_4h", 3), rep("E2_Dex_4h", 3))
aaColumn <- c("WT_4h", "Dex_4h")
aaContrasts <- cbind(c("Vehicle", "Dex_4h"),
                     c("E2_4h", "E2_Dex_4h"))

# Differential expression calculations
GSE79761.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                conditions = aaCondition,
                                                replicates = aaRepicates,
                                                contrasts = aaContrasts,
                                                expMat=exprs(GSE79761.exp.soft),
                                                adjpvalThresh = 1,
                                                lfcThresh = 0,
                                                expfData = fData(GSE79761.exp.soft),
                                                exppData = pData(GSE79761.exp.soft))
# adjpval1e-6:
# [1] "number of diff exp probes in contrast  E2_4h - Vehicle  :  0"
# [1] "number of diff exp probes in contrast  E2_Dex_4h - Dex_4h  :  0"

# lfc > 2:
# [1] "number of diff exp probes in contrast  E2_4h - Vehicle  :  35"
# [1] "number of diff exp probes in contrast  E2_Dex_4h - Dex_4h  :  15"

# lfc > 1:
# [1] "number of diff exp probes in contrast  E2_4h - Vehicle  :  420"
# [1] "number of diff exp probes in contrast  E2_Dex_4h - Dex_4h  :  310"
match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE79761.exp.soft)$ENTREZID)

#########################################################################################################
#########################################################################################################
#########################################################################################################

# 0 # 7 #10. : dropped for now because changes need to be more dramatic: GSE77244 : 0.01, 0.1 or 1 nm E2 : downloaded : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE77244_family.soft
# the zipped file is also there, delete the soft file after reading


#########################################################################################################
#########################################################################################################
#########################################################################################################

# 7 # 8 #11. GSE57935 : SMRT siRNA : downloaded : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE57935_family.soft
# the zipped file is also there, delete the soft file after reading

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE57935.exp.soft <- getGEO("GSE57935", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE57935.exp.soft <- GSE57935.exp.soft$GSE57935_series_matrix.txt.gz
annotation(GSE57935.exp.soft)
library(hugene10sttranscriptcluster.db)
# convert annotation
GSE57935.exp.soft <- annotateEset(GSE57935.exp.soft, hugene10sttranscriptcluster.db)
GSE57935.exp.soft <- normaliseIllumina(GSE57935.exp.soft, method = "quantile", transform = "log2")
GSE57935.exp.soft <- GSE57935.exp.soft[!is.na(fData(GSE57935.exp.soft)$ENTREZID), ]
head(fData(GSE57935.exp.soft))
boxplot(GSE57935.exp.soft)
dim(GSE57935.exp.soft)

aa <- apply(X = exprs(GSE57935.exp.soft), MARGIN = 1, FUN = var)
aa[ match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE57935.exp.soft)$ENTREZID)] <- 0.1
quantile(aa)
GSE57935.exp.soft <- GSE57935.exp.soft[aa > quantile(aa)[2], ]


#differential expression
aa <- strsplit(as.character(pData(GSE57935.exp.soft)$title), split="_" )
aaRepicates <- rep(c(1,2,3), 8)
aaCondition <- paste(unlist(lapply(aa, "[[", 2)),
                     unlist(lapply(aa, "[[", 3)),
                     unlist(lapply(aa, "[[", 4)),
                     sep = "_")

aaColumn <- c("WT_4h", "WT_24h", "siSMRT_4h", "siSMRT_24h")
aaContrasts <- cbind(c("siCtrl_4h_Veh", "siCtrl_24h_Veh", "siSMRT_4h_Veh", "siSMRT_24h_Veh"),
                     c("siCtrl_4h_E2", "siCtrl_24h_E2","siSMRT_4h_E2", "siSMRT_24h_E2"))

# Differential expression calculations
GSE57935.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                conditions = aaCondition,
                                                replicates = aaRepicates,
                                                contrasts = aaContrasts,
                                                expMat=exprs(GSE57935.exp.soft),
                                                adjpvalThresh = 1,
                                                lfcThresh = 0,
                                                expfData = fData(GSE57935.exp.soft),
                                                exppData = pData(GSE57935.exp.soft))

# adjpval 1e-6:
# [1] "number of diff exp probes in contrast  siCtrl_4h_E2 - siCtrl_4h_Veh  :  95"
# [1] "number of diff exp probes in contrast  siCtrl_24h_E2 - siCtrl_24h_Veh  :  261"
# [1] "number of diff exp probes in contrast  siSMRT_4h_E2 - siSMRT_4h_Veh  :  73"
# [1] "number of diff exp probes in contrast  siSMRT_24h_E2 - siSMRT_24h_Veh  :  156"

# lfc > 2:
# [1] "number of diff exp probes in contrast  siCtrl_4h_E2 - siCtrl_4h_Veh  :  2"
# [1] "number of diff exp probes in contrast  siCtrl_24h_E2 - siCtrl_24h_Veh  :  19"
# [1] "number of diff exp probes in contrast  siSMRT_4h_E2 - siSMRT_4h_Veh  :  3"
# [1] "number of diff exp probes in contrast  siSMRT_24h_E2 - siSMRT_24h_Veh  :  5"

# lfc > 1:
# [1] "number of diff exp probes in contrast  siCtrl_4h_E2 - siCtrl_4h_Veh  :  54"
# [1] "number of diff exp probes in contrast  siCtrl_24h_E2 - siCtrl_24h_Veh  :  132"
# [1] "number of diff exp probes in contrast  siSMRT_4h_E2 - siSMRT_4h_Veh  :  47"
# [1] "number of diff exp probes in contrast  siSMRT_24h_E2 - siSMRT_24h_Veh  :  159"

match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE57935.exp.soft)$ENTREZID)

#########################################################################################################
#########################################################################################################
#########################################################################################################

#8 # 9 #12. GSE26740 : si_AP-2g : downloaded : non_normalized
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE26740_non-normalized.txt

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE26740.exp.soft <- getGEO("GSE26740", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE26740.exp.soft <- GSE26740.exp.soft$GSE26740_series_matrix.txt.gz
annotation(GSE26740.exp.soft)
head(fData(GSE26740.exp.soft))
pData(GSE26740.exp.soft)
boxplot(GSE26740.exp.soft)

# Normalization of this expression dataset:
# GSM658200 Raw data was analyzed using GeneSpring GX 11.0 software. Dataset was normalized using percentile shift normalization and baseline transformed to median of control samples.
# Hence I will use a hiuristic approach to make the descrete expression matrix
# I use the difference in expression values, and threshold it at 0.6
GSE26740.exp.soft.descrete <- matrix(0L, nrow = nrow(exprs(GSE26740.exp.soft)), ncol = 2)
colnames(GSE26740.exp.soft.descrete) <- c("WT_12h", "Si_AP2g_12h")
rownames(GSE26740.exp.soft.descrete) <- rownames(exprs(GSE26740.exp.soft))

boxplot(exprs(GSE26740.exp.soft[,2]) - exprs(GSE26740.exp.soft[,1]), outline = F) 
boxplot(exprs(GSE26740.exp.soft[,4]) - exprs(GSE26740.exp.soft[,3]), outline = F) 
boxplot(exprs(GSE26740.exp.soft[,6]) - exprs(GSE26740.exp.soft[,5]), outline = F) 
###
# For comupting fold changes
aa_fc_add <- matrix(nrow = nrow(exprs(GSE26740.exp.soft)), ncol = 2)
aa_fc_add[, 1]<- exprs(GSE26740.exp.soft[,2]) - exprs(GSE26740.exp.soft[,1])
aa_fc_add[, 2] <- (exprs(GSE26740.exp.soft[,4]) - exprs(GSE26740.exp.soft[,3]) +  exprs(GSE26740.exp.soft[,6]) - exprs(GSE26740.exp.soft[,5]))/2
colnames(aa_fc_add) <- c("WT_12h", "Si_AP2g_12h")
GSE26740.exp.soft_fc <- aa_fc_add
###

GSE26740.exp.soft.descrete[(exprs(GSE26740.exp.soft[,2]) - exprs(GSE26740.exp.soft[,1])) >= 0.6, 1] <- 1
GSE26740.exp.soft.descrete[(exprs(GSE26740.exp.soft[,2]) - exprs(GSE26740.exp.soft[,1])) <= -0.6, 1] <- -1
aaa1 <- intersect(which((exprs(GSE26740.exp.soft[,4]) - exprs(GSE26740.exp.soft[,3])) >=  0.595),
                which((exprs(GSE26740.exp.soft[,6]) - exprs(GSE26740.exp.soft[,5])) >=  0.595))
length(unique(fData(GSE26740.exp.soft)$Entrez_Gene_ID[aaa1]))
aaa1 <- unique(fData(GSE26740.exp.soft)$Entrez_Gene_ID[aaa1])

aaa2 <- intersect(which((exprs(GSE26740.exp.soft[,4]) - exprs(GSE26740.exp.soft[,3])) <= -0.495),
                which((exprs(GSE26740.exp.soft[,6]) - exprs(GSE26740.exp.soft[,5])) <= -0.495))
length(unique(fData(GSE26740.exp.soft)$Entrez_Gene_ID[aaa2]))
aaa2 <- unique(fData(GSE26740.exp.soft)$Entrez_Gene_ID[aaa2])

aaa3 <- intersect(which(abs(exprs(GSE26740.exp.soft[,4]) - exprs(GSE26740.exp.soft[,3])) <=  0.1),
                  which(abs(exprs(GSE26740.exp.soft[,6]) - exprs(GSE26740.exp.soft[,5])) <=  0.1))
length(unique(fData(GSE26740.exp.soft)$Entrez_Gene_ID[aaa3]))
aaa3 <- unique(fData(GSE26740.exp.soft)$Entrez_Gene_ID[aaa3])

GSE26740.exp.soft.descrete[aa, 2] <- 1
GSE26740.exp.soft.descrete[ab, 2] <- -1

aa <- list()
aa$output_matrix <- matrix(nrow = length(unique(fData(GSE26740.exp.soft)$Entrez_Gene_ID)),
                           ncol = 2)
colnames(aa$output_matrix ) <- c("WT_12h", "Si_AP2g_12h")
rownames(aa$output_matrix ) <- as.character(unique(fData(GSE26740.exp.soft)$Entrez_Gene_ID))
aajj <- (exprs(GSE26740.exp.soft[,2]) - exprs(GSE26740.exp.soft[,1]))
aaups1 <- as.character(fData(GSE26740.exp.soft)$Entrez_Gene_ID[sort(aajj, 
                                                                    decreasing = T,
                                                                    index.return = T)$ix])
aadowns1 <- as.character(fData(GSE26740.exp.soft)$Entrez_Gene_ID[sort(aajj, 
                                                                      decreasing = F,
                                                                      index.return = T)$ix])
aaups1 <- aaups1[!duplicated(aaups1)]
aadowns1 <- aadowns1[!duplicated(aadowns1)]
aanochange1 <- as.character(fData(GSE26740.exp.soft)$Entrez_Gene_ID[abs(aajj) < 0.1])
aanochange1 <- aanochange1[!duplicated(aanochange1)]

aaChecked1 <- CheckAndReplace(UPSorted = aaups1, DOWNSorted = aadowns1, nu_gn = 100)

aa$output_matrix[match(aaChecked1$up, rownames(aa$output_matrix)), 1] <- 1
aa$output_matrix[match(aaChecked1$down, rownames(aa$output_matrix)), 1] <- -1
aa$output_matrix[match(aanochange1, rownames(aa$output_matrix)), 1] <- 0
# 
aa$output_matrix[match(as.character(aaa1), rownames(aa$output_matrix)), 2] <- 1
aa$output_matrix[match(as.character(aaa2), rownames(aa$output_matrix)), 2] <- -1
aa$output_matrix[match(as.character(aaa3), rownames(aa$output_matrix)), 2] <- 0

aa$UP_DOWN_Lists <- list()
aa$UP_DOWN_Lists$UPgenes <- list()
aa$UP_DOWN_Lists$UPgenes[[1]] <- aaChecked1$up
aa$UP_DOWN_Lists$UPgenes[[2]] <- as.character(aaa1)
names(aa$UP_DOWN_Lists$UPgenes) <- c("WT_12h", "Si_AP2g_12h")
aa$UP_DOWN_Lists$DOWNgenes <- list()
aa$UP_DOWN_Lists$DOWNgenes[[1]] <- aaChecked1$down
aa$UP_DOWN_Lists$DOWNgenes[[2]] <- as.character(aaa2)
names(aa$UP_DOWN_Lists$DOWNgenes) <- c("WT_12h", "Si_AP2g_12h")
aa$UP_DOWN_Lists$NoChangeGenes[[1]] <- aanochange1
aa$UP_DOWN_Lists$NoChangeGenes[[2]] <- as.character(aaa3)
names(aa$UP_DOWN_Lists$NoChangeGenes) <- c("WT_12h", "Si_AP2g_12h")
TopDownMatrixList_addendum <- aa

match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE26740.exp.soft)$Entrez_Gene_ID)

#########################################################################################################
#########################################################################################################
#########################################################################################################

#9 # 10 #13. GSE68918 : si_KDM3A : downloaded : non_normalized
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE68918_non_normalized.txt

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE68918.exp.soft <- getGEO("GSE68918", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE68918.exp.soft <- GSE68918.exp.soft$GSE68918_series_matrix.txt.gz
annotation(GSE68918.exp.soft)
head(fData(GSE68918.exp.soft))
pData(GSE68918.exp.soft)
boxplot(GSE68918.exp.soft)
GSE68918.exp.soft <- GSE68918.exp.soft[!is.na(fData(GSE68918.exp.soft)$Entrez_Gene_ID), ]

aa <- apply(X = exprs(GSE68918.exp.soft), MARGIN = 1, FUN = var)
aa[match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE68918.exp.soft)$Entrez_Gene_ID)[14]] <- 0.1 # to keep PBX1 in dataset

quantile(aa)
GSE68918.exp.soft <- GSE68918.exp.soft[aa > quantile(aa)[2], ]

#differential expression
aa <- gsub(x = as.character(pData(GSE68918.exp.soft)$title), pattern = " - E2 ", replacement = "_vehicle_")
aa <- gsub(x = aa, pattern = " \\+ E2 ", replacement = "_E2_")
aa <- gsub(x = aa, pattern = "_Replicate [[:digit:]]", replacement = "")
aa <- gsub(x = aa, pattern = "-", replacement = "_")
aaRepicates <- c(rep(1,4), rep(2, 4), rep(3, 3))
aaCondition <- aa
aaColumn <- c("WT_4h", "siKDM3A_B_4h")
aaContrasts <- cbind(c("siSCR_vehicle", "siKDM3A_B_vehicle"),
                     c("siSCR_E2", "siKDM3A_B_E2"))

# Differential expression calculations
GSE68918.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                conditions = aaCondition,
                                                replicates = aaRepicates,
                                                contrasts = aaContrasts,
                                                expMat=exprs(GSE68918.exp.soft),
                                                adjpvalThresh = 1,
                                                lfcThresh = 0,
                                                expfData = fData(GSE68918.exp.soft),
                                                exppData = pData(GSE68918.exp.soft))
#padj 1e-6:
# [1] "number of diff exp probes in contrast  siSCR_E2 - siSCR_vehicle  :  19"
# [1] "number of diff exp probes in contrast  siKDM3A_B_E2 - siKDM3A_B_vehicle  :  7"

# lfc > 2:
# [1] "number of diff exp probes in contrast  siSCR_E2 - siSCR_vehicle  :  3"
# [1] "number of diff exp probes in contrast  siKDM3A_B_E2 - siKDM3A_B_vehicle  :  3"

# lfc > 1:
# [1] "number of diff exp probes in contrast  siSCR_E2 - siSCR_vehicle  :  78"
# [1] "number of diff exp probes in contrast  siKDM3A_B_E2 - siKDM3A_B_vehicle  :  85"

match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE68918.exp.soft)$Entrez_Gene_ID)

#########################################################################################################
#########################################################################################################
#########################################################################################################

#10 # 11 #14. GSE24592 : si_ERK1, si_ERK2 : Zeynep : downloaded : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE24592_family.soft
# the zipped file is also there, delete the soft file after reading

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE24592.exp.soft <- getGEO("GSE24592",
                            destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE24592.exp.soft <- GSE24592.exp.soft$GSE24592_series_matrix.txt.gz
annotation(GSE24592.exp.soft)
head(fData(GSE24592.exp.soft))
aa <- annotateEset(GSE24592.exp.soft, hgu133a2.db)
sum(is.na(fData(aa)$ENTREZID))#1199
sum((fData(GSE24592.exp.soft)$ENTREZ_GENE_ID == ""))#1304
GSE24592.exp.soft <- annotateEset(GSE24592.exp.soft, hgu133a2.db)
GSE24592.exp.soft <- GSE24592.exp.soft[!is.na(fData(GSE24592.exp.soft)$ENTREZID), ]
pData(GSE24592.exp.soft)[, c(38:43)]
dim(GSE24592.exp.soft)


aa <- apply(X = exprs(GSE24592.exp.soft), MARGIN = 1, FUN = var)
quantile(aa)
GSE24592.exp.soft <- GSE24592.exp.soft[aa > quantile(aa)[2], ]


#differential expression
aaCondition <- paste(pData(GSE24592.exp.soft)[, 39],
                     pData(GSE24592.exp.soft)[, 42],
                     pData(GSE24592.exp.soft)[, 43], sep = "_")

aaColumn <- c("WT_4h", "WT_24h", "siERK1_4h", "siERK1_24h", "siERK2_4h", "siERK2_24h")
aaContrasts <- cbind(c("EtOH_siCtrl_4h", "EtOH_siCtrl_24h", "EtOH_siERK1_4h", "EtOH_siERK1_24h", "EtOH_siERK2_4h", "EtOH_siERK2_24h"),
                     c("E2_siCtrl_4h",    "E2_siCtrl_24h",  "E2_siERK1_4h",   "E2_siERK1_24h",   "E2_siERK2_4h",   "E2_siERK2_24h"))

# Differential expression calculations
GSE24592.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                conditions = aaCondition,
                                                #replicates = aaRepicates,
                                                contrasts = aaContrasts,
                                                expMat=exprs(GSE24592.exp.soft),
                                                adjpvalThresh = 1,
                                                lfcThresh = 0,
                                                expfData = fData(GSE24592.exp.soft),
                                                exppData = pData(GSE24592.exp.soft),
                                                use_arrayWeights = T)
# adjpval:
# [1] "number of diff exp probes in contrast  E2_siCtrl_4h - EtOH_siCtrl_4h  :  331"
# [1] "number of diff exp probes in contrast  E2_siCtrl_24h - EtOH_siCtrl_24h  :  4"
# [1] "number of diff exp probes in contrast  E2_siERK1_4h - EtOH_siERK1_4h  :  305"
# [1] "number of diff exp probes in contrast  E2_siERK1_24h - EtOH_siERK1_24h  :  105"
# [1] "number of diff exp probes in contrast  E2_siERK2_4h - EtOH_siERK2_4h  :  433"
# [1] "number of diff exp probes in contrast  E2_siERK2_24h - EtOH_siERK2_24h  :  8"

#  lfc > 2
# [1] "number of diff exp probes in contrast  E2_siCtrl_4h - EtOH_siCtrl_4h  :  139"
# [1] "number of diff exp probes in contrast  E2_siCtrl_24h - EtOH_siCtrl_24h  :  678"
# [1] "number of diff exp probes in contrast  E2_siERK1_4h - EtOH_siERK1_4h  :  136"
# [1] "number of diff exp probes in contrast  E2_siERK1_24h - EtOH_siERK1_24h  :  413"
# [1] "number of diff exp probes in contrast  E2_siERK2_4h - EtOH_siERK2_4h  :  204"
# [1] "number of diff exp probes in contrast  E2_siERK2_24h - EtOH_siERK2_24h  :  687"

# lfc > 1
# [1] "number of diff exp probes in contrast  E2_siCtrl_4h - EtOH_siCtrl_4h  :  567"
# [1] "number of diff exp probes in contrast  E2_siCtrl_24h - EtOH_siCtrl_24h  :  3356"
# [1] "number of diff exp probes in contrast  E2_siERK1_4h - EtOH_siERK1_4h  :  657"
# [1] "number of diff exp probes in contrast  E2_siERK1_24h - EtOH_siERK1_24h  :  2075"
# [1] "number of diff exp probes in contrast  E2_siERK2_4h - EtOH_siERK2_4h  :  850"
# [1] "number of diff exp probes in contrast  E2_siERK2_24h - EtOH_siERK2_24h  :  3142"

match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE24592.exp.soft)$ENTREZID)

#########################################################################################################
#########################################################################################################
#########################################################################################################

#11 # 12 #15. GSE53394 : si_ERK5 : downloaded : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE53394_family.soft
# the zipped file is also there, delete the soft file after reading

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE53394.exp.soft <- getGEO("GSE53394", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE53394.exp.soft <- GSE53394.exp.soft$GSE53394_series_matrix.txt.gz
annotation(GSE53394.exp.soft)
head(fData(GSE53394.exp.soft))
aa <- annotateEset(GSE53394.exp.soft, hgu133a.db)
sum(is.na(fData(aa)$ENTREZID))#1205
sum(fData(GSE53394.exp.soft)$ENTREZ_GENE_ID == "")#1310
sum(is.na(fData(GSE53394.exp.soft)$ENTREZ_GENE_ID))
GSE53394.exp.soft <- annotateEset(GSE53394.exp.soft, hgu133a.db)

#GSE53394.exp.soft <- GSE53394.exp.soft[!(fData(GSE53394.exp.soft)$ENTREZ_GENE_ID == ""), ]
GSE53394.exp.soft <- normaliseIllumina(GSE53394.exp.soft, method = "quantile", transform = "log2")
GSE53394.exp.soft <- GSE53394.exp.soft[!is.na(fData(GSE53394.exp.soft)$ENTREZID), ]
pData(GSE53394.exp.soft)
boxplot(GSE53394.exp.soft)
dim(GSE53394.exp.soft)

aa <- apply(X = exprs(GSE53394.exp.soft), MARGIN = 1, FUN = var)
aa[match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE53394.exp.soft)$ENTREZID)[17]] <- 0.1 #to keep SP1 in dataset
quantile(aa)
GSE53394.exp.soft <- GSE53394.exp.soft[aa > quantile(aa)[2], ]

#differential expression
aaCondition <- unlist(lapply(strsplit(x = as.character(pData(GSE53394.exp.soft)$title), split = " "), "[[", 1))
aaRepicates <- rep(c(1,2), 4)
aaColumn <- c("WT_24h", "siERK5_24h")
aaContrasts <- cbind(c("siGL3Veh", "siERK5Veh"),
                     c("siGL3E2",  "siERK5E2"))

# Differential expression calculations
GSE53394.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                conditions = aaCondition,
                                                replicates = aaRepicates,
                                                contrasts = aaContrasts,
                                                expMat=exprs(GSE53394.exp.soft),
                                                adjpvalThresh = 1,
                                                lfcThresh = 0,
                                                expfData = fData(GSE53394.exp.soft),
                                                exppData = pData(GSE53394.exp.soft),
                                                use_arrayWeights = T)

# adjpval 1e-6:
# [1] "number of diff exp probes in contrast  siGL3E2 - siGL3Veh  :  0"
# [1] "number of diff exp probes in contrast  siERK5E2 - siERK5Veh  :  0"

# lfc > 2:
# [1] "number of diff exp probes in contrast  siGL3E2 - siGL3Veh  :  0"
# [1] "number of diff exp probes in contrast  siERK5E2 - siERK5Veh  :  0"

# lfc > 1:
# [1] "number of diff exp probes in contrast  siGL3E2 - siGL3Veh  :  11"
# [1] "number of diff exp probes in contrast  siERK5E2 - siERK5Veh  :  22"

match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE53394.exp.soft)$ENTREZID)


#########################################################################################################
#########################################################################################################
#########################################################################################################

#12 # 13 #16. GSE8597 : protein synthesis inhibition : think if this should be included : downloaded : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE8597_family.soft
# the zipped file is also there, delete the soft file after reading

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE8597.exp.soft <- getGEO("GSE8597", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE8597.exp.soft <- GSE8597.exp.soft$GSE8597_series_matrix.txt.gz
annotation(GSE8597.exp.soft)
head(fData(GSE8597.exp.soft))
GSE8597.exp.soft <- annotateEset(GSE8597.exp.soft, hgu133plus2.db  )
sum(is.na(fData(GSE8597.exp.soft)$ENTREZID))
GSE8597.exp.soft <- GSE8597.exp.soft[!is.na(fData(GSE8597.exp.soft)$ENTREZID), ]
pData(GSE8597.exp.soft)
boxplot(GSE8597.exp.soft)
dim(GSE8597.exp.soft)

aa <- apply(X = exprs(GSE8597.exp.soft), MARGIN = 1, FUN = var)
quantile(aa)
GSE8597.exp.soft <- GSE8597.exp.soft[aa > quantile(aa)[2], ]


#differential expression
aa <- gsub("MCF7_", "", as.character(pData(GSE8597.exp.soft)$title))
aa <- gsub(pattern = "_rep[[:digit:]]", replacement = "", x = aa)
aaCondition <- aa
aaRepicates <- rep(c(1,2,3,4), 4)
aaColumn <- c("WT_24h", "CHX_24h")
aaContrasts <- cbind(c("EtOH_24h", "CHX_EtOH_24h"),
                     c("E2_24h",  "CHX_E2_24h"))

# Differential expression calculations
GSE8597.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                conditions = aaCondition,
                                                replicates = aaRepicates,
                                                contrasts = aaContrasts,
                                                expMat=exprs(GSE8597.exp.soft),
                                                adjpvalThresh = 1,
                                                lfcThresh = 0,
                                                expfData = fData(GSE8597.exp.soft),
                                                exppData = pData(GSE8597.exp.soft),
                                                use_arrayWeights = T)

# adjpval 1e-6:
# [1] "number of diff exp probes in contrast  E2_24h - EtOH_24h  :  911"
# [1] "number of diff exp probes in contrast  CHX_E2_24h - CHX_EtOH_24h  :  337"

# lfc > 2:
# [1] "number of diff exp probes in contrast  E2_24h - EtOH_24h  :  37"
# [1] "number of diff exp probes in contrast  CHX_E2_24h - CHX_EtOH_24h  :  30"

# lfc > 1
# [1] "number of diff exp probes in contrast  E2_24h - EtOH_24h  :  599"
# [1] "number of diff exp probes in contrast  CHX_E2_24h - CHX_EtOH_24h  :  237"
match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE8597.exp.soft)$ENTREZID)


#########################################################################################################
#########################################################################################################
#########################################################################################################

#13 # 14 #17. GSE4006 : +ERbeta : downloaded : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE4006_family.soft
# the zipped file is also there, delete the soft file after reading

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE4006.exp.soft <- getGEO("GSE4006", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE4006.exp.soft <- GSE4006.exp.soft$GSE4006_series_matrix.txt.gz
GSE4006.exp.soft <- normaliseIllumina(GSE4006.exp.soft, method = "quantile", transform = "log2")

annotation(GSE4006.exp.soft)
head(fData(GSE4006.exp.soft))
#look at bioconductor annotation
aa <- annotateEset(GSE4006.exp.soft, hgu133a.db  )
#check which annotation is better
sum(fData(GSE4006.exp.soft)$ENTREZ_GENE_ID == "")
sum(is.na(fData(aa)$ENTREZID))
#set the new annotation
GSE4006.exp.soft <- annotateEset(GSE4006.exp.soft, hgu133a.db  )
sum(is.na(fData(GSE4006.exp.soft)$ENTREZID))

GSE4006.exp.soft <- GSE4006.exp.soft[!is.na(fData(GSE4006.exp.soft)$ENTREZID), ]
pData(GSE4006.exp.soft)
boxplot(GSE4006.exp.soft)
dim(GSE4006.exp.soft)

aa <- apply(X = exprs(GSE4006.exp.soft), MARGIN = 1, FUN = var)
quantile(aa)
GSE4006.exp.soft <- GSE4006.exp.soft[aa > quantile(aa)[2], ]

#differential expression
aa <- as.character(pData(GSE4006.exp.soft)$title)
aa[11] <- "AdERb(moi50)+E2 replicate1"
aa[12] <- "AdERb(moi50)+E2 replicate2"
aa <- unlist(lapply(X = strsplit(x = aa, split = " "), "[[", 1))
aa <-   gsub("\\+", "_", aa)
aa <-   gsub("\\(", "_", aa)
aa <-   gsub("\\)", "", aa)
aaCondition <- aa
aaRepicates <- rep(c(1,2), 6)
aaColumn <- c("WT_24h", "ERb_moi5_24h", "ERb_moi50_24h")
aaContrasts <- cbind(c("Ad_veh", "AdERb_moi5_veh", "AdERb_moi50_veh"),
                     c("Ad_E2",  "AdERb_moi5_E2",  "AdERb_moi50_E2"))

# Differential expression calculations
GSE4006.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                               conditions = aaCondition,
                                               replicates = aaRepicates,
                                               contrasts = aaContrasts,
                                               expMat=exprs(GSE4006.exp.soft),
                                               adjpvalThresh = 1,
                                               lfcThresh = 0,
                                               expfData = fData(GSE4006.exp.soft),
                                               exppData = pData(GSE4006.exp.soft),
                                               use_arrayWeights = T)
#adj pval: 1e-6:
# [1] "number of diff exp probes in contrast  Ad_E2 - Ad_veh  :  8"
# [1] "number of diff exp probes in contrast  AdERb_moi5_E2 - AdERb_moi5_veh  :  0"
# [1] "number of diff exp probes in contrast  AdERb_moi50_E2 - AdERb_moi50_veh  :  0"

# lfc > 2
# [1] "number of diff exp probes in contrast  Ad_E2 - Ad_veh  :  162"
# [1] "number of diff exp probes in contrast  AdERb_moi5_E2 - AdERb_moi5_veh  :  136"
# [1] "number of diff exp probes in contrast  AdERb_moi50_E2 - AdERb_moi50_veh  :  172"

# lfc > 1
# [1] "number of diff exp probes in contrast  Ad_E2 - Ad_veh  :  783"
# [1] "number of diff exp probes in contrast  AdERb_moi5_E2 - AdERb_moi5_veh  :  745"
# [1] "number of diff exp probes in contrast  AdERb_moi50_E2 - AdERb_moi50_veh  :  936"
match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE4006.exp.soft)$ENTREZID)


#########################################################################################################
#########################################################################################################
#########################################################################################################

#14 # 15 #18. GSE74146 : siFGFR2 : downloaded : non_normalized
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE74146_non-normalized.txt

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE74146.exp.soft <- getGEO("GSE74146", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE74146.exp.soft <- GSE74146.exp.soft$GSE74146_series_matrix.txt.gz
boxplot(GSE74146.exp.soft)
annotation(GSE74146.exp.soft)
aa <- annotateEset(GSE74146.exp.soft, illuminaHumanv4.db)

head(fData(GSE74146.exp.soft))
head(fData(aa))

sum(is.na(fData(aa)$ENTREZID)) #11597
sum(is.na(fData(GSE74146.exp.soft)$Entrez_Gene_ID)) #3363
#apparently the original annotation is better, so I'll keep it.
GSE74146.exp.soft <- GSE74146.exp.soft[!is.na(fData(GSE74146.exp.soft)$Entrez_Gene_ID),]

aafunc <- function(x){
  return(any(is.na(x)))
}
aa <- apply(exprs(GSE74146.exp.soft),1,aafunc)
GSE74146.exp.soft <- GSE74146.exp.soft[!aa, ]

aa <- apply(X = exprs(GSE74146.exp.soft), MARGIN = 1, FUN = var)
aa[match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE74146.exp.soft)$Entrez_Gene_ID)] <- 0.1 #to keep TFs in the dataset
quantile(aa)
GSE74146.exp.soft <- GSE74146.exp.soft[aa > quantile(aa)[2], ]
dim(GSE74146.exp.soft)

#differential expression
pData(GSE74146.exp.soft)
aa <- as.character(pData(GSE74146.exp.soft)$title)
aa2 <- paste(unlist(lapply(X = strsplit(x = aa, split = " "), "[[", 2)),
             unlist(lapply(X = strsplit(x = aa, split = " "), "[[", 3)),
             sep = "_")
aaCondition <- aa2
aaRepicates <- as.integer(gsub(pattern = "Rep",
                               replacement = "",
                               x = unlist(lapply(X = strsplit(x = aa, split = " "), "[[", 4))))
aaColumn <- c("WT_6h", "siFGFR2_6h")
aaContrasts <- cbind(c("CTL_Starved", "siFGFR2_Starved"),
                     c("CTL_E2",      "siFGFR2_E2"))

# Differential expression calculations
GSE74146.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                               conditions = aaCondition,
                                               replicates = aaRepicates,
                                               contrasts = aaContrasts,
                                               expMat=exprs(GSE74146.exp.soft),
                                               adjpvalThresh = 1,
                                               lfcThresh = 0,
                                               expfData = fData(GSE74146.exp.soft),
                                               exppData = pData(GSE74146.exp.soft),
                                               use_arrayWeights = T)
# adj p-val 1e-6:
# [1] "number of diff exp probes in contrast  CTL_E2 - CTL_Starved  :  477"
# [1] "number of diff exp probes in contrast  siFGFR2_E2 - siFGFR2_Starved  :  316"

# lfc > 1:
# [1] "number of diff exp probes in contrast  CTL_E2 - CTL_Starved  :  67"
# [1] "number of diff exp probes in contrast  siFGFR2_E2 - siFGFR2_Starved  :  40"

# lfc > 2:
# [1] "number of diff exp probes in contrast  CTL_E2 - CTL_Starved  :  3"
# [1] "number of diff exp probes in contrast  siFGFR2_E2 - siFGFR2_Starved  :  3"
match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE74146.exp.soft)$Entrez_Gene_ID)

#########################################################################################################
#########################################################################################################
#########################################################################################################

#15 # 16 #19. GSE56245 : EGCG (green tea) : downloaded : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE56245_family.soft
# the zipped file is also there, delete the soft file after reading

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE56245.exp.soft <- getGEO("GSE56245", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE56245.exp.soft <- GSE56245.exp.soft$GSE56245_series_matrix.txt.gz
boxplot(GSE56245.exp.soft)
annotation(GSE56245.exp.soft)

sum(is.na(fData(GSE56245.exp.soft)$GeneName))#0
sum((fData(GSE56245.exp.soft)$GeneName == ""))#0
#apparently the original annotation is good. I have to convert to ENTREZ ID though
aa <- geneNameToEntrez(inputNames = fData(GSE56245.exp.soft)$GeneName) # there are 31033 unique gene names, of which 20402 can be converted to ENTREZ ID
aa$ENTREZID[920] <- "367" #correct AR ENTREZID
aa$ENTREZID[9619] <- "6667" #correct SP1 ENTREZID
aaa <- match(fData(GSE56245.exp.soft)$GeneName, aa$ALIAS)
fData(GSE56245.exp.soft)$ENTREZID <- aa$ENTREZID[aaa]
GSE56245.exp.soft <- GSE56245.exp.soft[!is.na(fData(GSE56245.exp.soft)$ENTREZID), ]

#variance check:
aa <- apply(X = exprs(GSE56245.exp.soft), MARGIN = 1, FUN = var)
aa[match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE56245.exp.soft)$ENTREZID)] <- 0.1
quantile(aa)
GSE56245.exp.soft <- GSE56245.exp.soft[aa > quantile(aa)[2], ]
dim(GSE56245.exp.soft)


#differential expression
pData(GSE56245.exp.soft)
aa <- as.character(pData(GSE56245.exp.soft)$title)
aa <-unlist(lapply(X = strsplit(x = aa, split = "-"), "[[", 2))
aa <- gsub(pattern = "\\+",replacement = "_",x = aa)
aaCondition <- aa
aaRepicates <- rep(c(1,2), 4)
aaColumn <- c("WT_24h", "EGCG_24h")
aaContrasts <- cbind(c("EtOH", "EGCG"),
                     c("E2",   "E2_EGCG"))

# Differential expression calculations
GSE56245.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                conditions = aaCondition,
                                                replicates = aaRepicates,
                                                contrasts = aaContrasts,
                                                expMat=exprs(GSE56245.exp.soft),
                                                adjpvalThresh = 1,
                                                lfcThresh = 0,
                                                expfData = fData(GSE56245.exp.soft),
                                                exppData = pData(GSE56245.exp.soft),
                                                use_arrayWeights = T)
#adjpval 1e-6:
# [1] "number of diff exp probes in contrast  E2 - EtOH  :  0"
# [1] "number of diff exp probes in contrast  E2_EGCG - EGCG  :  0"

# lfc > 2
# [1] "number of diff exp probes in contrast  E2 - EtOH  :  151"
# [1] "number of diff exp probes in contrast  E2_EGCG - EGCG  :  238"


# lfc > 1:
# [1] "number of diff exp probes in contrast  E2 - EtOH  :  1231"
# [1] "number of diff exp probes in contrast  E2_EGCG - EGCG  :  1239"

match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE56245.exp.soft)$ENTREZID)

#########################################################################################################
#########################################################################################################
#########################################################################################################

#0 :dropping because diff exp can't be done, also no replicate, shady experiment: #16 # 17 #20. GSE39786 : H2ac knockdown : downloaded : This is raw data soft file
# : maybe not suitable: read the paper: An H2A histone isotype regulates estrogen receptor target genes by mediating enhancer-promoter-3'-UTR interactions in breast cancer cells. Nucleic Acids Res 2014
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE39786_family.soft
# the zipped file is also there, delete the soft file after reading

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE39786.exp.soft <- getGEO("GSE39786", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE39786.exp.soft <- GSE39786.exp.soft$GSE39786_series_matrix.txt.gz
boxplot(GSE39786.exp.soft)
GSE39786.exp.soft <- normaliseIllumina(GSE39786.exp.soft, method = "quantile", transform = "log2")
boxplot(GSE39786.exp.soft)

annotation(GSE39786.exp.soft)
aa <- annotateEset(GSE39786.exp.soft, hgu133plus2.db)

head(fData(GSE39786.exp.soft))
head(fData(aa))

sum(is.na(fData(aa)$ENTREZID)) #10324
sum(fData(GSE39786.exp.soft)$ENTREZ_GENE_ID == "") #10541
#apparently the new annotation is better, so I'll change it.
GSE39786.exp.soft <- annotateEset(GSE39786.exp.soft, hgu133plus2.db)
GSE39786.exp.soft <- GSE39786.exp.soft[!is.na(fData(GSE39786.exp.soft)$ENTREZID),]

#differential expression
pData(GSE39786.exp.soft)
aaCondition <- c("control_veh", "control_E2", "H2acKD_veh", "H2acKD_E2")
aaColumn <- c("WT_4h", "H2acKD_4h")
aaContrasts <- cbind(c("control_veh", "H2acKD_veh"),
                     c("control_E2",   "H2acKD_E2"))

# Differential expression calculations
GSE39786.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                conditions = aaCondition,
                                                #replicates = aaRepicates,
                                                contrasts = aaContrasts,
                                                expMat=exprs(GSE39786.exp.soft),
                                                adjpvalThresh = 1,
                                                lfcThresh = 0,
                                                expfData = fData(GSE39786.exp.soft),
                                                exppData = pData(GSE39786.exp.soft),
                                                use_arrayWeights = T)
#DROPPED

#########################################################################################################
#########################################################################################################
#########################################################################################################

#16 # 18 #21. GSE36586 : siFOS : downloaded : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE36586_family.soft
# the zipped file is also there, delete the soft file after reading

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE36586.exp.soft <- getGEO("GSE36586", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE36586.exp.soft <- GSE36586.exp.soft$GSE36586_series_matrix.txt.gz
boxplot(GSE36586.exp.soft)
GSE36586.exp.soft <- normaliseIllumina(GSE36586.exp.soft, method = "quantile", transform = "log2")
boxplot(GSE36586.exp.soft)

annotation(GSE36586.exp.soft)
#aa <- annotateEset(GSE36586.exp.soft, hgu133plus2.db)

head(fData(GSE36586.exp.soft))
#head(fData(aa))

#sum(is.na(fData(aa)$ENTREZID)) #10324
sum(fData(GSE36586.exp.soft)$GENE_SYMBOL == "") #11969
sum(fData(GSE36586.exp.soft)$GENE_NAME == "") #12572
sum(is.na(fData(GSE36586.exp.soft)$GENE)) #12572
GSE36586.exp.soft <- GSE36586.exp.soft[!is.na(fData(GSE36586.exp.soft)$GENE), ]
colnames(fData(GSE36586.exp.soft))

#apparently the original annotation is good. I have to convert to ENTREZ ID though
aa <- geneNameToEntrez(inputNames = fData(GSE36586.exp.soft)$GENE_SYMBOL)
aa$ENTREZID[665] <- "367" #correcting AR entrez id
aa$ENTREZID[18185] <- "6667" #correcting SP1 entrez id
aaa <- match(fData(GSE36586.exp.soft)$GENE_SYMBOL, aa$ALIAS)
fData(GSE36586.exp.soft)$ENTREZID <- aa$ENTREZID[aaa]
GSE36586.exp.soft <- GSE36586.exp.soft[!is.na(fData(GSE36586.exp.soft)$ENTREZID), ]
dim(GSE36586.exp.soft)

#variance check
aa <- apply(X = exprs(GSE36586.exp.soft), MARGIN = 1, FUN = var)
aa[match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE36586.exp.soft)$ENTREZID)] <- 0.1
quantile(aa)
GSE36586.exp.soft <- GSE36586.exp.soft[aa > quantile(aa)[2], ]
dim(GSE36586.exp.soft)


#differential expression
pData(GSE36586.exp.soft)$title
aa <- as.character(pData(GSE36586.exp.soft)$title)
aa2 <- paste(unlist(lapply(X = strsplit(x = aa, split = "_"), "[[", 2)),
             unlist(lapply(X = strsplit(x = aa, split = "_"), "[[", 3)),
             unlist(lapply(X = strsplit(x = aa, split = "_"), "[[", 4)),
             sep = "_")
aaCondition <- aa2
aaRepicates <- as.integer(gsub(pattern = "rep",
                               replacement = "",
                               x = unlist(lapply(X = strsplit(x = aa, split = "_"), "[[", 5))))
aaColumn <- c("WT_24h", "siFos_24h")
aaContrasts <- cbind(c("siControl_vehicle_24h", "siFos_vehicle_24h"),
                     c("siControl_E2_24h",      "siFos_E2_24h"))

# Differential expression calculations
GSE36586.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                conditions = aaCondition,
                                                replicates = aaRepicates,
                                                contrasts = aaContrasts,
                                                expMat=exprs(GSE36586.exp.soft),
                                                adjpvalThresh = 1,
                                                lfcThresh = 0,
                                                expfData = fData(GSE36586.exp.soft),
                                                exppData = pData(GSE36586.exp.soft),
                                                use_arrayWeights = T)

#adjpval: 1e-6
# [1] "number of diff exp probes in contrast  siControl_E2_24h - siControl_vehicle_24h  :  591"
# [1] "number of diff exp probes in contrast  siFos_E2_24h - siFos_vehicle_24h  :  444"

# lfc > 2
# [1] "number of diff exp probes in contrast  siControl_E2_24h - siControl_vehicle_24h  :  244"
# [1] "number of diff exp probes in contrast  siFos_E2_24h - siFos_vehicle_24h  :  96"

# lfc > 1
# [1] "number of diff exp probes in contrast  siControl_E2_24h - siControl_vehicle_24h  :  1598"
# [1] "number of diff exp probes in contrast  siFos_E2_24h - siFos_vehicle_24h  :  886"
match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE36586.exp.soft)$ENTREZID)
#########################################################################################################
#########################################################################################################
#########################################################################################################

#17 # 19 #22. GSE35428 : 4-hydroxytamoxifen, ICI-182,780, Raloxifene, Bazedoxifene and Lasofoxifene with or without E2 : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE35428_family.soft
# the zipped file is also there, delete the soft file after reading

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE35428.exp.soft <- getGEO("GSE35428", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE35428.exp.soft <- GSE35428.exp.soft$GSE35428_series_matrix.txt.gz
boxplot(GSE35428.exp.soft)
#GSE35428.exp.soft <- normaliseIllumina(GSE35428.exp.soft, method = "quantile", transform = "log2")
#boxplot(GSE35428.exp.soft)

annotation(GSE35428.exp.soft)
aa <- annotateEset(GSE35428.exp.soft, hgu133plus2.db)

head(fData(GSE35428.exp.soft))
head(fData(aa))

sum(is.na(fData(aa)$ENTREZID)) #10324
sum(fData(GSE35428.exp.soft)$ENTREZ_GENE_ID == "") #10541
#apparently the new annotation is better, so I'll change it.
GSE35428.exp.soft <- annotateEset(GSE35428.exp.soft, hgu133plus2.db)
GSE35428.exp.soft <- GSE35428.exp.soft[!is.na(fData(GSE35428.exp.soft)$ENTREZID),]

#variance check
aa <- apply(X = exprs(GSE35428.exp.soft), MARGIN = 1, FUN = var)
aa[match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE35428.exp.soft)$ENTREZID)] <- 0.1
quantile(aa)
GSE35428.exp.soft <- GSE35428.exp.soft[aa > quantile(aa)[2], ]
dim(GSE35428.exp.soft)

#differential expression
pData(GSE35428.exp.soft)$'co-treatment:ch1'
aa <- as.character(pData(GSE35428.exp.soft)$'co-treatment:ch1')
aa <- gsub(pattern = " ",replacement = "_", x = aa)
aa <- gsub(pattern = ",",replacement = "_", x = aa)

aaCondition <- aa

aaColumn <- c("ethanol_24h", "Raloxifene_24h", "4OHT_24h", "Bazedoxifene_24h", "ICI_182_780_24h", "Lasofoxifene_24h")
aaContrasts <- cbind(c("ethanol_ethanol", "ethanol_Raloxifene", "ethanol_4OHT", "ethanol_Bazedoxifene","ethanol_ICI_182_780", "ethanol_Lasofoxifene"),
                     c("E2_ethanol",      "E2_Raloxifene"     , "E2_4OHT"     , "E2_Bazedoxifene",     "E2_ICI_182_780"     , "E2_Lasofoxifene"))

# Differential expression calculations
GSE35428.exp.soft.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                conditions = aaCondition,
                                                #replicates = aaRepicates,
                                                contrasts = aaContrasts,
                                                expMat=exprs(GSE35428.exp.soft),
                                                adjpvalThresh = 1,
                                                lfcThresh = 0,
                                                expfData = fData(GSE35428.exp.soft),
                                                exppData = pData(GSE35428.exp.soft),
                                                use_arrayWeights = T)

# adjpval: 1e-6: 
# [1] "number of diff exp probes in contrast  E2_ethanol - ethanol_ethanol  :  7753"
# [1] "number of diff exp probes in contrast  E2_Raloxifene - ethanol_Raloxifene  :  2"
# [1] "number of diff exp probes in contrast  E2_4OHT - ethanol_4OHT  :  292"
# [1] "number of diff exp probes in contrast  E2_Bazedoxifene - ethanol_Bazedoxifene  :  0"
# [1] "number of diff exp probes in contrast  E2_ICI_182_780 - ethanol_ICI_182_780  :  19"
# [1] "number of diff exp probes in contrast  E2_Lasofoxifene - ethanol_Lasofoxifene  :  7"

# lfc>2:
# [1] "number of diff exp probes in contrast  E2_ethanol - ethanol_ethanol  :  356"
# [1] "number of diff exp probes in contrast  E2_Raloxifene - ethanol_Raloxifene  :  0"
# [1] "number of diff exp probes in contrast  E2_4OHT - ethanol_4OHT  :  9"
# [1] "number of diff exp probes in contrast  E2_Bazedoxifene - ethanol_Bazedoxifene  :  0"
# [1] "number of diff exp probes in contrast  E2_ICI_182_780 - ethanol_ICI_182_780  :  0"
# [1] "number of diff exp probes in contrast  E2_Lasofoxifene - ethanol_Lasofoxifene  :  1"

# lfc>1:
# [1] "number of diff exp probes in contrast  E2_ethanol - ethanol_ethanol  :  1756"
# [1] "number of diff exp probes in contrast  E2_Raloxifene - ethanol_Raloxifene  :  0"
# [1] "number of diff exp probes in contrast  E2_4OHT - ethanol_4OHT  :  74"
# [1] "number of diff exp probes in contrast  E2_Bazedoxifene - ethanol_Bazedoxifene  :  0"
# [1] "number of diff exp probes in contrast  E2_ICI_182_780 - ethanol_ICI_182_780  :  2"
# [1] "number of diff exp probes in contrast  E2_Lasofoxifene - ethanol_Lasofoxifene  :  9"
match(TF.ENTREZ.Shrinked$ENTREZID, fData(GSE35428.exp.soft)$ENTREZID)

#########################################################################################################
#########################################################################################################
#########################################################################################################

#0 : DROP THIS DATASET because its done on mmultiple platforms and diff exp results in 0 diff exp genes possibly because of median centered normalization#18 # 20 #23. GSE13458 : NMNAT1 knockdown : downloaded : This is raw data soft file
# /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Expression_data/Descretized_Expression_dataset/GSE13458_family.soft
# the zipped file is also there, delete the soft file after reading

# import the dataset and convert the annotation to illumianexpressionset so it can be analyzed with bead array
GSE13458.exp.soft <- getGEO("GSE13458", destdir = "Expression_data/Descretized_Expression_dataset/Automatic_download/")
GSE13458.exp.soft_1 <- GSE13458.exp.soft$`GSE13458-GPL570_series_matrix.txt.gz`
GSE13458.exp.soft_2 <- GSE13458.exp.soft$`GSE13458-GPL571_series_matrix.txt.gz`

boxplot(GSE13458.exp.soft_1)
boxplot(GSE13458.exp.soft_2)


annotation(GSE13458.exp.soft_1)
annotation(GSE13458.exp.soft_2)

aa_1 <- annotateEset(GSE13458.exp.soft_1, hgu133plus2.db)
aa_2 <- annotateEset(GSE13458.exp.soft_2, hgu133a2.db)

head(fData(GSE13458.exp.soft_1))
head(fData(aa_1))

head(fData(GSE13458.exp.soft_2))
head(fData(aa_2))

sum(is.na(fData(aa_2)$ENTREZID)) #1199
sum(fData(GSE13458.exp.soft_2)$ENTREZ_GENE_ID == "")#1304


sum(fData(GSE13458.exp.soft_1)$ENTREZ_GENE_ID == "") #10541
sum(is.na(fData(aa_1)$ENTREZID)) #10324

#apparently the new annotation is better, so I'll change it.
GSE13458.exp.soft_1 <- annotateEset(GSE13458.exp.soft_1, hgu133plus2.db)
GSE13458.exp.soft_2 <- annotateEset(GSE13458.exp.soft_2, hgu133a2.db)

GSE13458.exp.soft_1 <- GSE13458.exp.soft_1[!is.na(fData(GSE13458.exp.soft_1)$ENTREZID),]
GSE13458.exp.soft_2 <- GSE13458.exp.soft_2[!is.na(fData(GSE13458.exp.soft_2)$ENTREZID),]


#differential expression: I will just work with GSE13458.exp.soft_2, since the other one has only 4 conditions, which are replicates of these

aaCondition <- c("lucifer_KD_control", "lucifer_KD_control", "lucifer_KD_E2", "lucifer_KD_E2", "NMNAT1_KD_control", "NMNAT1_KD_control", "NMNAT1_KD_E2", "NMNAT1_KD_E2")
aaRepicates <- rep(c(1,3), 4)
aaColumn <- c("WT_3h", "NMNAT1_KD_3h")
aaContrasts <- cbind(c("lucifer_KD_control", "NMNAT1_KD_control"),
                     c("lucifer_KD_E2",      "NMNAT1_KD_E2"))

# Differential expression calculations
GSE13458.exp.soft_2.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                conditions = aaCondition,
                                                replicates = aaRepicates,
                                                contrasts = aaContrasts,
                                                expMat=exprs(GSE13458.exp.soft_2),
                                                adjpvalThresh = 1,
                                                lfcThresh = 2,
                                                expfData = fData(GSE13458.exp.soft_2),
                                                exppData = pData(GSE13458.exp.soft_2),
                                                use_arrayWeights = T)
#########################################################################################################
#########################################################################################################
############################                 Revisit RNA-seq Data                ########################
############################         Can they be used in the new framework?      ########################
#########################################################################################################
#########################################################################################################

######################################################
###################################################### GSE73663 Zeynep
GSE73663.Inh.RNAseq


#differential expression: 
#exp1 : control: 5"AZDVeh4h1" other : 1 "AZD334h1", 3 AZDE24h1
#exp2 : control: 6"AZDVeh4h2" other : 2"AZD334h2",4 AZDE24h2
#exp3 : control: 17 PP242Veh4h1 other: 13 PP242334h1 , 15 PP242E24h1
#exp4 : control: 18 PP242Veh4h2 other: 14 PP242334h2 ,16 PP242E24h2
#exp5 : control: 11 CtrlVeh4h1 other:  7 Ctrl334h1 , 9 CtrlE24h1
#exp6 : control: 12 CtrlVeh4h2 other:  8 Ctrl334h2 , 10 CtrlE24h2
GSE73663_RNA_Seq_ordered <- GSE73663.Inh.RNAseq[, c(5, 3,
                                                    6, 4,
                                                    17, 15,
                                                    18, 16,
                                                    11, 9,
                                                    12, 10)]
GSE73663_RNA_Seq_ordered <- getCommonGenes(list(GSE73663_RNA_Seq_ordered))[[1]]

#adding in the TFs
match("AR", rownames(GSE73663.Inh.RNAseq))
match("SP1",rownames(GSE73663.Inh.RNAseq))
GSE73663_RNA_Seq_ordered <- rbind(GSE73663_RNA_Seq_ordered, GSE73663.Inh.RNAseq[match("AR", rownames(GSE73663.Inh.RNAseq)),  c(5, 3, 6, 4, 17, 15, 18, 16, 11, 9, 12, 10)])
GSE73663_RNA_Seq_ordered <- rbind(GSE73663_RNA_Seq_ordered, GSE73663.Inh.RNAseq[match("SP1", rownames(GSE73663.Inh.RNAseq)), c(5, 3, 6, 4, 17, 15, 18, 16, 11, 9, 12, 10)])

rownames(GSE73663_RNA_Seq_ordered)[nrow(GSE73663_RNA_Seq_ordered)-1] <- TF.ENTREZ.Shrinked$ENTREZID[1]
rownames(GSE73663_RNA_Seq_ordered)[nrow(GSE73663_RNA_Seq_ordered)] <- TF.ENTREZ.Shrinked$ENTREZID[17]



GSE73663_RNA_Seq_ordered <- log2(GSE73663_RNA_Seq_ordered + 1)

# #variance check
# aa <- apply(X = GSE73663_RNA_Seq_ordered, MARGIN = 1, FUN = var)
# aa[match(TF.ENTREZ.Shrinked$ENTREZID, rownames(GSE73663_RNA_Seq_ordered))] <- 1
# quantile(aa)
# GSE73663_RNA_Seq_ordered <- GSE73663_RNA_Seq_ordered[aa > quantile(aa)[2], ]
# dim(GSE73663_RNA_Seq_ordered)



aa <- c("MAPK_Veh_1", "MAPK_E2_1",
        "MAPK_Veh_2", "MAPK_E2_2",
        "mTOR_Veh_1", "mTOR_E2_1",
        "mTOR_Veh_2", "mTOR_E2_2",
        "Veh_Veh_1", "Veh_E2_1",
        "Veh_Veh_2", "Veh_E2_2")
colnames(GSE73663_RNA_Seq_ordered) <- aa
aaCondition <- c("MAPK_control", "MAPK_E2",
                 "MAPK_control", "MAPK_E2",
                 "mTOR_control", "mTOR_E2",
                 "mTOR_control", "mTOR_E2",
                 "WT_control", "WT_E2",
                 "WT_control", "WT_E2")
aaRepicates <- rep(c(1, 1, 2, 2), 3)
aaColumn <- c("WT_4h", "MAPK_4h", "mTOR_4h")
aaContrasts <- cbind(c("WT_control", "MAPK_control", "mTOR_control"),
                     c("WT_E2",      "MAPK_E2",      "mTOR_E2"))

# Differential expression calculations

GSE73663_RNA_Seq_ordered.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                  conditions = aaCondition,
                                                  replicates = aaRepicates,
                                                  contrasts = aaContrasts,
                                                  expMat=GSE73663_RNA_Seq_ordered,
                                                  adjpvalThresh = 1,
                                                  lfcThresh = 0,
                                                  #expfData = fData(GSE13458.exp.soft_2),
                                                  #exppData = pData(GSE13458.exp.soft_2),
                                                  use_arrayWeights = F)
GSE73663_RNA_Seq_ordered.Descrete$designMatrix
GSE73663_RNA_Seq_ordered.Descrete_topdown <- DiffExpToTopGeneMat(topTable_input = GSE73663_RNA_Seq_ordered.Descrete$FitContrasts,
                                              nu_contrasts = ncol( GSE73663_RNA_Seq_ordered.Descrete$descrete_expMat),
                                              nu_top_genes = 100,
                                              expMatRowNames = rownames( GSE73663_RNA_Seq_ordered.Descrete$FullexpMat),
                                              mappedID =  rownames( GSE73663_RNA_Seq_ordered.Descrete$FullexpMat),
                                              is_entrez=T,
                                              output_colnames = colnames( GSE73663_RNA_Seq_ordered.Descrete$descrete_expMat))

aa <- match(rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), rownames(GSE73663_RNA_Seq_ordered.Descrete_topdown$output_matrix))
aa <- GSE73663_RNA_Seq_ordered.Descrete_topdown$output_matrix[aa, ]

######################################################
###################################################### GSE86316 MEN1 KD
rownames(GSE86316.RNAseq)
colnames(GSE86316.RNAseq) #2 Men1 kd   ctrl: 1,2
"kd_smpl_1_s1_treat"               #Veh_ShCtrl_1
"kd_smpl_2_s1_treat"               #Veh_ShCtrl_2
"MCF7_E2_rep1_s1_treat"            #E2_ShCtrl_1 10nM estradiol 45 mins
"MCF7_E2_rep2_s1_treat"            #E2_ShCtrl_2 10nM estradiol 45 mins
"MCF7_Veh_KD_MEN1_rep1_s1_treat"   #Veh_ShMEN1_1
"MCF7_Veh_KD_MEN1_rep2_s1_treat"   #Veh_ShMEN1_2
"MCF7_E2_KD_MEN1_rep1_s1_treat"    #E2_ShMEN1_1 10nM estradiol 45 mins
"kd_smpl_8_s1_treat"               #E2_ShMEN1_2 10nM estradiol 45 mins
####################
#exp1 : control: Veh_ShCtrl_1        other: E2_ShCtrl_1 
#exp2 : control: Veh_ShCtrl_2        other: E2_ShCtrl_2
#exp3 : control: Veh_ShMEN1_1        other: E2_ShMEN1_1
#exp4 : control: Veh_ShMEN1_2        other: E2_ShMEN1_2

GSE86316_RNAseq_ordered <- GSE86316.RNAseq[, c(1, 3,
                                                2, 4,
                                                5, 7,
                                                6, 8)]
aa <- c("Veh_ShCtrl_1", "E2_ShCtrl_1",
        "Veh_ShCtrl_2", "E2_ShCtrl_2",
        "Veh_ShMEN1_1", "E2_ShMEN1_1",
        "Veh_ShMEN1_2", "E2_ShMEN1_2")
colnames(GSE86316_RNAseq_ordered) <- aa
GSE86316_RNAseq_ordered <- getCommonGenes(list(GSE86316_RNAseq_ordered))[[1]]
match(TF.ENTREZ.Shrinked$ENTREZID, rownames(GSE86316_RNAseq_ordered))

GSE86316_RNAseq_ordered <- log2(GSE86316_RNAseq_ordered + 1)


aaCondition <- c("WT_control", "WT_E2",
                 "WT_control", "WT_E2",
                 "MEN1_KD_control", "MEN1_KD_E2",
                 "MEN1_KD_control", "MEN1_KD_E2")
aaRepicates <- rep(c(1, 1, 2, 2), 2)
aaColumn <- c("WT_45min", "MEN1KD_45min")
aaContrasts <- cbind(c("WT_control", "MEN1_KD_control"),
                     c("WT_E2",      "MEN1_KD_E2"))

# Differential expression calculations

GSE86316_RNAseq_ordered.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                       conditions = aaCondition,
                                                       replicates = aaRepicates,
                                                       contrasts = aaContrasts,
                                                       expMat=GSE86316_RNAseq_ordered,
                                                       adjpvalThresh = 1,
                                                       lfcThresh = 0,
                                                       #expfData = fData(GSE13458.exp.soft_2),
                                                       #exppData = pData(GSE13458.exp.soft_2),
                                                       use_arrayWeights = F)
GSE86316_RNAseq_ordered.Descrete$designMatrix
GSE86316_RNAseq_ordered.Descrete_topdown <- DiffExpToTopGeneMat(topTable_input = GSE86316_RNAseq_ordered.Descrete$FitContrasts,
                                                                 nu_contrasts = ncol( GSE86316_RNAseq_ordered.Descrete$descrete_expMat),
                                                                 nu_top_genes = 100,
                                                                 expMatRowNames = rownames( GSE86316_RNAseq_ordered.Descrete$FullexpMat),
                                                                 mappedID =  rownames( GSE86316_RNAseq_ordered.Descrete$FullexpMat),
                                                                 is_entrez=T,
                                                                 output_colnames = colnames( GSE86316_RNAseq_ordered.Descrete$descrete_expMat))

aa <- match(rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), rownames(GSE86316_RNAseq_ordered.Descrete_topdown$output_matrix))
aa <- GSE86316_RNAseq_ordered.Descrete_topdown$output_matrix[aa, ]
length(aa)
sum(is.na(aa))
######################################################
###################################################### GSE55922 BRD4 KD # it was previously dropped because it doesn't have FOXA1
ncol(GSE55922.RNAseq)
colnames(GSE55922.RNAseq) #3 BRD4 kd     cntrl: 5,6
#"brd4.e2a"     #siBRD4 -- Vehicle
#"brd4.e2b"     #siBRD4 -- Vehicle
#"brd4pluse2a"  #siBRD4 -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
#"brd4pluse2b"  #siBRD4 -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
#"cont.e2a"     #sictrl -- Vehicle
#"cont.e2b"     #sictrl -- Vehicle
#"contpluse2a"  #sictrl -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
#"contpluse2b"  #sictrl -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
#"jq1.e2a"      #treated with JQ1 which is a BRD4 inhibitor -- Vehicle
#"jq1.e2b"      #treated with JQ1 which is a BRD4 inhibitor -- Vehicle
#"jq1pluse2a"   #treated with JQ1 which is a BRD4 inhibitor -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
#"jq1pluse2b"  #treated with JQ1 which is a BRD4 inhibitor -- 10 nmol/L 17-β-estradiol (Sigma-Aldrich) for 2 hours.
##############
#exp1: control: brd4.e2a other: brd4pluse2a
#exp2: control: brd4.e2b other: brd4pluse2b
#exp3: control: cont.e2a other: contpluse2a
#exp4: control: cont.e2b other: contpluse2b
#exp5: control: jq1.e2a other: jq1pluse2a
#exp6: control: jq1.e2b other: jq1pluse2b

GSE55922_RNAseq_ordered <- GSE55922.RNAseq[, c(1, 3,
                                               2, 4,
                                               5, 7,
                                               6, 8,
                                               9, 11,
                                               10,12)]
aa <- c("siBRD4_Veh_1", "siBRD4_E2_1",
        "siBRD4_Veh_2", "siBRD4_E2_2",
        "sictrl_Veh_1", "sictrl_E2_1",
        "sictrl_Veh_2", "sictrl_E2_2",
        "JQ1_Veh_1", "JQ1_E2_1",
        "JQ1_Veh_2", "JQ1_E2_2")
colnames(GSE55922_RNAseq_ordered) <- aa
GSE55922_RNAseq_ordered <- getCommonGenes(list(GSE55922_RNAseq_ordered))[[1]]
aa <- geneNameToEntrez(rownames(GSE55922.RNAseq))

# "FOX" %in% toupper(rownames(GSE55922.RNAseq))
# aa <- GSE55922.RNAseq[grep(pattern = "FOX", x = toupper(rownames(GSE55922.RNAseq))), ]
# boxplot.matrix(t(aa), las = 2)
GSE55922_RNAseq_ordered <-rbind(GSE55922_RNAseq_ordered, GSE55922.RNAseq[match("AR" , rownames(GSE55922.RNAseq)), ])
GSE55922_RNAseq_ordered <-rbind(GSE55922_RNAseq_ordered, GSE55922.RNAseq[match("SP1" , rownames(GSE55922.RNAseq)), ])
rownames(GSE55922_RNAseq_ordered)[nrow(GSE55922_RNAseq_ordered)-1] <- TF.ENTREZ.Shrinked$ENTREZID[1]
rownames(GSE55922_RNAseq_ordered)[nrow(GSE55922_RNAseq_ordered)] <- TF.ENTREZ.Shrinked$ENTREZID[17]



match(TF.ENTREZ.Shrinked$ENTREZID, rownames(GSE55922_RNAseq_ordered))

GSE55922_RNAseq_ordered <- log2(GSE55922_RNAseq_ordered + 1)


aaCondition <- c("BRD4_KD_control", "BRD4_KD_E2",
                 "BRD4_KD_control", "BRD4_KD_E2",
                 "WT_control", "WT_E2",
                 "WT_control", "WT_E2",
                 "JQ1_control", "JQ1_E2",
                 "JQ1_control", "JQ1_E2")
aaRepicates <- rep(c(1, 1, 2, 2), 3)
aaColumn <- c("WT_2h", "BRD4KD_2h", "JQ1_2h")
aaContrasts <- cbind(c("WT_control", "BRD4_KD_control", "JQ1_control"),
                     c("WT_E2",      "BRD4_KD_E2", "JQ1_E2"))

# Differential expression calculations

GSE55922_RNAseq_ordered.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                      conditions = aaCondition,
                                                      replicates = aaRepicates,
                                                      contrasts = aaContrasts,
                                                      expMat=GSE55922_RNAseq_ordered,
                                                      adjpvalThresh = 1,
                                                      lfcThresh = 0,
                                                      #expfData = fData(GSE13458.exp.soft_2),
                                                      #exppData = pData(GSE13458.exp.soft_2),
                                                      use_arrayWeights = F)
GSE55922_RNAseq_ordered.Descrete_topdown <- DiffExpToTopGeneMat(topTable_input = GSE55922_RNAseq_ordered.Descrete$FitContrasts,
                                                                nu_contrasts = ncol( GSE55922_RNAseq_ordered.Descrete$descrete_expMat),
                                                                nu_top_genes = 100,
                                                                expMatRowNames = rownames( GSE55922_RNAseq_ordered.Descrete$FullexpMat),
                                                                mappedID =  rownames( GSE55922_RNAseq_ordered.Descrete$FullexpMat),
                                                                is_entrez=T,
                                                                output_colnames = colnames( GSE55922_RNAseq_ordered.Descrete$descrete_expMat))

aa <- match(rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), rownames(GSE55922_RNAseq_ordered.Descrete_topdown$output_matrix))
aa <- GSE55922_RNAseq_ordered.Descrete_topdown$output_matrix[aa, ]
length(aa)
sum(is.na(aa))
######################################################
###################################################### GSE76507 shG9a, shPHF20
colnames(GSE76507.RNAseq)
#1"shG9a_NOE2"       
#2"shG9a_E2"          
#3"shNT_NOE2"         
#4"shNT_E2"          
#5"shNT_NoE2.rep1"  
#6"shNT_NOE2.rep2"    
#7"shNT_E2.rep1"     
#8"shNT_E2.rep2"     
#9"shPHF20_NOE2.rep1" 
#10"shPHF20_NOE2.rep2"
#11"shPHF20_E2.rep1"   
#12"shPHF20_E2.rep2"  
##########
#exp1: control: [1]shG9a_NOE2 other: [2]shG9a_E2
#exp2: control: [3]shNT_NOE2, other: [4]shNT_E2
#exp3: control: [5]shNT_NoE2.rep1, other: [7]shNT_E2.rep1
#exp4: control: [6]shNT_NoE2.rep2, other: [8]shNT_E2.rep2
#exp5: control: [9]shPHF20_NOE2.rep1, other: [11]shPHF20_E2.rep1
#exp6: control: [10]shPHF20_NOE2.rep2, other: [12]shPHF20_E2.rep2
GSE76507_RNAseq_ordered <- GSE76507.RNAseq[, c(1, 2, 3, 4, 5, 7, 6, 8, 9, 11, 10, 12)]
aa <- c("shG9a_Veh", "shG9a_E2",
        "shNT_Veh_1", "shNT_E2_1",
        "shNT_Veh_2", "shNT_E2_2",
        "shNT_Veh_3", "shNT_E2_3",
        "shPHF20_Veh_1", "shPHF20_E2_1",
        "shPHF20_Veh_2", "shPHF20_E2_2")
colnames(GSE76507_RNAseq_ordered) <- aa
GSE76507_RNAseq_ordered <- getCommonGenes(list(GSE76507_RNAseq_ordered))[[1]]
aa <- geneNameToEntrez(rownames(GSE55922.RNAseq))
# checking if all TFs are present
match(TF.ENTREZ.Shrinked$ENTREZID, rownames(GSE76507_RNAseq_ordered))

GSE76507_RNAseq_ordered <-rbind(GSE76507_RNAseq_ordered, GSE76507.RNAseq[match("AR" , rownames(GSE76507.RNAseq)), ])
GSE76507_RNAseq_ordered <-rbind(GSE76507_RNAseq_ordered, GSE76507.RNAseq[match("SP1" , rownames(GSE76507.RNAseq)), ])
rownames(GSE76507_RNAseq_ordered)[nrow(GSE76507_RNAseq_ordered)-1] <- TF.ENTREZ.Shrinked$ENTREZID[1]
rownames(GSE76507_RNAseq_ordered)[nrow(GSE76507_RNAseq_ordered)] <- TF.ENTREZ.Shrinked$ENTREZID[17]

GSE76507_RNAseq_ordered <- log2(GSE76507_RNAseq_ordered + 1)

aaCondition <- c("shG9a_Veh", "shG9a_E2",
                 "shNT_Veh", "shNT_E2",
                 "shNT_Veh", "shNT_E2",
                 "shNT_Veh", "shNT_E2",
                 "shPHF20_Veh", "shPHF20_E2",
                 "shPHF20_Veh", "shPHF20_E2")
aaRepicates <- c(1, 1, 2, 2, rep(c(3,3,4,4), 2))
aaColumn <- c("WT_3h", "G9aKD_3h", "PHF20KD_3h")
aaContrasts <- cbind(c("shNT_Veh", "shG9a_Veh", "shPHF20_Veh"),
                     c("shNT_E2",  "shG9a_E2",  "shPHF20_E2"))

# Differential expression calculations

GSE76507_RNAseq_ordered.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                      conditions = aaCondition,
                                                      replicates = aaRepicates,
                                                      contrasts = aaContrasts,
                                                      expMat=GSE76507_RNAseq_ordered,
                                                      adjpvalThresh = 1,
                                                      lfcThresh = 0,
                                                      #expfData = fData(GSE13458.exp.soft_2),
                                                      #exppData = pData(GSE13458.exp.soft_2),
                                                      use_arrayWeights = F)
GSE76507_RNAseq_ordered.Descrete_topdown <- DiffExpToTopGeneMat(topTable_input = GSE76507_RNAseq_ordered.Descrete$FitContrasts,
                                                                nu_contrasts = ncol( GSE76507_RNAseq_ordered.Descrete$descrete_expMat),
                                                                nu_top_genes = 100,
                                                                expMatRowNames = rownames( GSE76507_RNAseq_ordered.Descrete$FullexpMat),
                                                                mappedID =  rownames( GSE76507_RNAseq_ordered.Descrete$FullexpMat),
                                                                is_entrez=T,
                                                                output_colnames = colnames( GSE76507_RNAseq_ordered.Descrete$descrete_expMat))

aa <- match(rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), rownames(GSE76507_RNAseq_ordered.Descrete_topdown$output_matrix))
aa <- GSE76507_RNAseq_ordered.Descrete_topdown$output_matrix[aa, ]
length(aa)
sum(is.na(aa))
######################################################
###################################################### GSE108308 #R4K1 stapled peptide
colnames(GSE108308.RNAseq)      
#1 Veh_Rep1
#2 Veh_Rep2
#3 E2_Rep1
#4 E2_Rep2
#5 TOT_Rep1
#6 TOT_Rep2
#7 ET_Rep1
#8 ET_Rep2
#9 Peptide_Rep1
#10 Peptide_Rep2
#11 E_Peptide_Rep1
#12 E_Peptide_Rep2
#################
#exp1: control: [1]Veh_Rep1, other: [3]E2_Rep1 ,[5]TOT_Rep1,[7]ET_Rep1, [9]Peptide_Rep1, [11]E_Peptide_Rep1
#exp2: control: [2]Veh_Rep2, other: [4]E2_Rep2 ,[6]TOT_Rep2,[8]ET_Rep2,[10]Peptide_Rep2, [12]E_Peptide_Rep2

GSE108308_RNAseq_ordered <- GSE108308.RNAseq[, c(1, 3,
                                                 2, 4,
                                                 5, 7,
                                                 6, 8,
                                                 9, 11,
                                                 10, 12)]
aa <- c("Veh_1", "E2_1",
        "Veh_2", "E2_2",
        "Tamxfn_1", "E2_Tamxfn_1",
        "Tamxfn_2", "E2_Tamxfn_2",
        "R4K1_1", "E2_R4K1_1",
        "R4K1_2", "E2_R4K1_2")

colnames(GSE108308_RNAseq_ordered) <- aa
GSE108308_RNAseq_ordered <- getCommonGenes(list(GSE108308_RNAseq_ordered))[[1]]
aa <- geneNameToEntrez(rownames(GSE108308.RNAseq))
# checking if all TFs are present
match(TF.ENTREZ.Shrinked$ENTREZID, rownames(GSE108308_RNAseq_ordered))

GSE108308_RNAseq_ordered <- log2(GSE108308_RNAseq_ordered + 1)
aaCondition <- c("Veh", "E2",
                 "Veh", "E2",
                 "Tamxfn", "E2_Tamxfn",
                 "Tamxfn", "E2_Tamxfn",
                 "R4K1", "E2_R4K1",
                 "R4K1", "E2_R4K1")

aaRepicates <-  rep(c(1,1,2,2), 3)
aaColumn <- c("WT_2h", "Tamxfn_2h", "R4K1_2h")
aaContrasts <- cbind(c("Veh", "Tamxfn", "R4K1"),
                     c("E2",  "E2_Tamxfn",  "E2_R4K1"))

# Differential expression calculations

GSE108308_RNAseq_ordered.Descrete <- createDescreteMat(columnNames = aaColumn,
                                                      conditions = aaCondition,
                                                      replicates = aaRepicates,
                                                      contrasts = aaContrasts,
                                                      expMat=GSE108308_RNAseq_ordered,
                                                      adjpvalThresh = 1,
                                                      lfcThresh = 0,
                                                      #expfData = fData(GSE13458.exp.soft_2),
                                                      #exppData = pData(GSE13458.exp.soft_2),
                                                      use_arrayWeights = F)
GSE108308_RNAseq_ordered.Descrete_topdown <- DiffExpToTopGeneMat(topTable_input = GSE108308_RNAseq_ordered.Descrete$FitContrasts,
                                                                nu_contrasts = ncol( GSE108308_RNAseq_ordered.Descrete$descrete_expMat),
                                                                nu_top_genes = 100,
                                                                expMatRowNames = rownames( GSE108308_RNAseq_ordered.Descrete$FullexpMat),
                                                                mappedID =  rownames( GSE108308_RNAseq_ordered.Descrete$FullexpMat),
                                                                is_entrez=T,
                                                                output_colnames = colnames( GSE108308_RNAseq_ordered.Descrete$descrete_expMat))

aa <- match(rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), rownames(GSE108308_RNAseq_ordered.Descrete_topdown$output_matrix))
aa <- GSE108308_RNAseq_ordered.Descrete_topdown$output_matrix[aa, ]
length(aa)
sum(is.na(aa))


#########################################################################################################
#########################################################################################################
#########################################################################################################
######################################### PULL DESCRETIZED MATRICES TOGETHER ############################

####################### NEW APPROACH: USE TOP 100 UP and DOWN regulated genes

my_Descretized_List <- list(GSE37386.exp.soft.filt.Descrete,
                          GSE25314.exp.soft.filt.Descrete,
                          GSE25315.exp.soft.filt.Descrete,
                          GSE102367.exp.Descrete,
                          GSE79761.exp.soft.Descrete,
                          GSE57935.exp.soft.Descrete,
                          GSE68918.exp.soft.Descrete,
                          GSE24592.exp.soft.Descrete,
                          GSE53394.exp.soft.Descrete,
                          GSE8597.exp.soft.Descrete,
                          GSE4006.exp.soft.Descrete,
                          GSE74146.exp.soft.Descrete,
                          GSE56245.exp.soft.Descrete,
                          GSE36586.exp.soft.Descrete,
                          GSE35428.exp.soft.Descrete)
names(my_Descretized_List) <- c("GSE37386",
                                "GSE25314",
                                "GSE25315",
                                "GSE102367",
                                "GSE79761",
                                "GSE57935",
                                "GSE68918",
                                "GSE24592",
                                "GSE53394",
                                "GSE8597",
                                "GSE4006",
                                "GSE74146",
                                "GSE56245",
                                "GSE36586",
                                "GSE35428")
TopDownMatrixList <- list()
for(i in 1:length(my_Descretized_List)){
  print(i)
  if(is.null(my_Descretized_List[[i]]$fDATA$ENTREZID)){
    my_Descretized_List[[i]]$fDATA$ENTREZID <- my_Descretized_List[[i]]$fDATA$Entrez_Gene_ID
  }
  my_Descretized_List[[i]]$fDATA$ENTREZID <- as.character(my_Descretized_List[[i]]$fDATA$ENTREZID)
  
  TopDownMatrixList[[i]] <- DiffExpToTopGeneMat(topTable_input = my_Descretized_List[[i]]$FitContrasts,
                                                nu_contrasts = ncol( my_Descretized_List[[i]]$descrete_expMat),
                                                nu_top_genes = 100,
                                                expMatRowNames = rownames(my_Descretized_List[[i]]$FullexpMat),
                                                mappedID =  my_Descretized_List[[i]]$fDATA$ENTREZID,
                                                is_entrez=T,
                                                output_colnames = colnames( my_Descretized_List[[i]]$descrete_expMat))
}




names(TopDownMatrixList) <- names(my_Descretized_List)
#check if there are any mutual genes in top and down lists
for(i in 1:length(TopDownMatrixList)){
  print(i)
  print(mapply(FUN = intersect,  TopDownMatrixList[[i]]$UP_DOWN_Lists$UPgenes, TopDownMatrixList[[i]]$UP_DOWN_Lists$DOWNgenes ))
}

TopDownMatrixList_addendum # this is for the AP2g KD experiment (GSE26740) which was not included in the above list
TopDownMatrixList[[16]] <- TopDownMatrixList_addendum
names(TopDownMatrixList)[16] <- "GSE26740"

my_Descretized_List[[16]] <- list()
my_Descretized_List[[16]]$FullexpMat <- exprs(GSE26740.exp.soft)
my_Descretized_List[[16]]$fDATA <- fData(GSE26740.exp.soft)
my_Descretized_List[[16]]$pDATA <- pData(GSE26740.exp.soft)
my_Descretized_List[[16]]$raw_lfc <- GSE26740.exp.soft_fc
names(my_Descretized_List)[16] <- "GSE26740"
my_Descretized_List[[16]]$fDATA$ENTREZID <- my_Descretized_List[[16]]$fDATA$Entrez_Gene_ID
#########################################################################################################
#########################################################################################################
#####################################  TopDownMatrixList is the list containing the 100UPDOWN matrices of 16 experiments
#########################################################################################################
#########################################################################################################

MyDifExpMatList <- lapply(TopDownMatrixList, "[[", 1)
names(MyDifExpMatList) <- names(TopDownMatrixList)
my_CommonDifExpMat_16 <- CommonDifExpMat(MyDifExpMatList)
################# creating the same matrix for RNA-seq data:
my_Descretized_List_RNAseq <- list(GSE73663_RNA_Seq_ordered.Descrete,
                                   GSE86316_RNAseq_ordered.Descrete,
                                   GSE55922_RNAseq_ordered.Descrete,
                                   GSE76507_RNAseq_ordered.Descrete,
                                   GSE108308_RNAseq_ordered.Descrete)
names(my_Descretized_List_RNAseq) <- c("GSE73663","GSE86316","GSE55922","GSE76507","GSE108308")

TopDownMatrixList_RNAseq <- list(GSE73663_RNA_Seq_ordered.Descrete_topdown,
                                 GSE86316_RNAseq_ordered.Descrete_topdown,
                                 GSE55922_RNAseq_ordered.Descrete_topdown,
                                 GSE76507_RNAseq_ordered.Descrete_topdown,
                                 GSE108308_RNAseq_ordered.Descrete_topdown)
names(TopDownMatrixList_RNAseq) <- c("GSE73663","GSE86316","GSE55922","GSE76507","GSE108308")
MyDifExpMatList_RNAseq <- lapply(TopDownMatrixList_RNAseq, "[[", 1)
names(MyDifExpMatList_RNAseq) <- names(TopDownMatrixList_RNAseq)
my_CommonDifExpMat_RNAseq <- CommonDifExpMat(MyDifExpMatList_RNAseq)
#########################################################################################################
#########################################################################################################
#####################################  my_CommonDifExpMat_16 is the matrix containing the modified expression profile for union of all genes -- > for microarray experiments
#####################################  my_CommonDifExpMat_RNAseq is the same as above for RNAseq data
#########################################################################################################
#########################################################################################################

#########################################################################################################
#########################################################################################################
             # ALL MICROARRAY DATA  IS SAVED IN: "~/Documents/Shayan/BioInf/EstrogenReceptor/ERPROJECT_microarray_data.RData"
#########################################################################################################
#########################################################################################################

colnames(my_CommonDifExpMat_16)[38] <- "GSE35428_WT_24h"
# Filter for the ER associated genes
length(Genes.Associated.REMAP.ER.Entrez)
length(intersect(Genes.Associated.REMAP.ER.Entrez, rownames(my_CommonDifExpMat_16)))

my_CommonDifExpMat_16_ERassoc <- my_CommonDifExpMat_16[ rownames(my_CommonDifExpMat_16) %in% Genes.Associated.REMAP.ER.Entrez, ]
hist(rowSums(my_CommonDifExpMat_16_ERassoc, na.rm = T), breaks = 100)

aalak <- my_CommonDifExpMat_16_ERassoc
aalak[is.na(aalak)] <- 20

aa <- apply(X = aalak, MARGIN = 1, FUN = (function(x) all(x == 0 | x == 20)))


#aa <- apply(X = my_CommonDifExpMat_16_ERassoc, MARGIN = 1, FUN = (function(x) all(is.na(x))))


my_CommonDifExpMat_16_ERassoc <- my_CommonDifExpMat_16_ERassoc[!aa, ]

aa1 <- apply(X = my_CommonDifExpMat_16_ERassoc, MARGIN = 1, FUN = (function(x) sum(!is.na(x))))
hist(aa1, breaks = 100)
table(aa1)
pie(table(aa1), main = "Genes / no. of non-na entries")
#########################################################################################################
# looking at the correlation of fold change for TFs against the 172 genes chosen


#########################################################################################################
# plot number of genes as a function of nonzero entries
aag <- integer(44)
for(i in 1:44){
  aag[i] <- sum(aa1 >= i)
}
plot(c(1:44), aag, ylab = "number of genes", xlab = "number of non-na entries")
abline(h = seq(0,1200, 50), col = 3, lty = 3)
#########################################################################################################
# excluding the wildtype conditions
my_CommonDifExpMat_16_ERassoc_noWT <- my_CommonDifExpMat_16_ERassoc[, -(which(unlist(lapply(strsplit(colnames(my_CommonDifExpMat_16_ERassoc), split="_"), "[[", 2)) == "WT"))]

aa2 <- apply(X = my_CommonDifExpMat_16_ERassoc_noWT, MARGIN = 1, FUN = (function(x) sum(!is.na(x))))
hist(aa2, breaks = 100)
table(aa2)
pie(table(aa2), main = "Genes/no. of nonWT non_na entries")

# plot number of genes as a function of non_na entries
aag2 <- integer(26)
for(i in 1:26){
  aag2[i] <- sum(aa2 >= i)
}
plot(c(1:26), aag2, ylab = "number of genes", xlab = "number of non-na entries in KD conds")
abline(h = seq(0,1200, 50), col = 3, lty = 3)

sum(aa2 > 3)
table(aa1[which(aa2 > 3)])


aa3 = aa1 - aa2
table(aa3)
pie(table(aa3))

aafakeExp <- my_CommonDifExpMat_16_ERassoc
aafakeExp[is.na(aafakeExp)] <- -100
aa4p <-  apply(X = aafakeExp, MARGIN = 1, FUN = (function(x) sum(x == 1)))
aa4n <-  apply(X = aafakeExp, MARGIN = 1, FUN = (function(x) sum(x == -1)))
aa4z <-  apply(X = aafakeExp, MARGIN = 1, FUN = (function(x) sum(x == 0)))

sum(aa4p > 5 | aa4n > 5 | aa4z > 5)
sum(aa4n > 0)
sum(aa4p > 0 & aa4n > 0) #126
sum(aa4p > 1 & aa4n > 1) #39

table(aa1[which(aa4p > 1 & aa4n > 1)])

plotExpression(my_CommonDifExpMat_16_ERassoc, filename = "All_genes_heatmap.png")

plotExpression(my_CommonDifExpMat_16_ERassoc_noWT, filename = "All_genes_heatmap_noWT.png")


aa1 <- apply(X = aafakeExp, MARGIN = 1, FUN = (function(x) sum(x != 0 & x != -100)))
plotExpression(my_CommonDifExpMat_16_ERassoc[aa1 >=5, ], filename = "genes_gte5nonzero_172_heatmap.png")
sum(aa1 >=5)#172
my_CommonDifExpMat_16_ERassoc_gte5nzAll_172 <- my_CommonDifExpMat_16_ERassoc[aa1 >=5, ]

## Creating the RNA-seq dataset associated with the 172 genes
sum(rownames(my_CommonDifExpMat_RNAseq)  %in% rownames(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172))
my_CommonDifExpMat_RNAseq_5_ERassoc_172 <- my_CommonDifExpMat_RNAseq[match(rownames(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172), 
                                                                           rownames(my_CommonDifExpMat_RNAseq)) ,]
sum(!is.na(my_CommonDifExpMat_RNAseq_5_ERassoc_172))
hist(rowSums(!is.na(my_CommonDifExpMat_RNAseq_5_ERassoc_172)), breaks = 20)
my_CommonDifExpMat_RNAseq_5_ERassoc_172_150 <- my_CommonDifExpMat_RNAseq_5_ERassoc_172[rowSums(!is.na(my_CommonDifExpMat_RNAseq_5_ERassoc_172)) > 0, ] 
plotExpression(my_CommonDifExpMat_RNAseq_5_ERassoc_172_150,
               filename = "RNA_Seq_172_150_notClustered.png", .Rowv = F, .Colv = F, .dendrogram = "none")


##################
aafakeExpKD <- aafakeExp[, -(which(unlist(lapply(strsplit(colnames(my_CommonDifExpMat_16_ERassoc), split="_"), "[[", 2)) == "WT"))]
aa2 <- apply(X = aafakeExpKD, MARGIN = 1, FUN = (function(x) sum(x != 0 & x != -100)))
sum(aa2 >=4)
plotExpression(my_CommonDifExpMat_16_ERassoc[aa2 >=5, ], filename = "genes_gte5nonzeroKD_heatmap.png")
plotExpression(my_CommonDifExpMat_16_ERassoc_noWT[aa2 >=5, ], filename = "genes_gte5nonzeroKD_noWT_heatmap.png")

aa3 <- intersect(which(aa4p > 0 & aa4n > 0), which(aa2 >=4))
plotExpression(my_CommonDifExpMat_16_ERassoc[aa3, ], filename = "genes_gte4nonzeroKD_atleast1p1n_heatmap.png")
#The matrix of gene expression with only genes that have at least 4 nonzero entries in KD conditions, AND have at least one +1 and one -1 
my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52 <- my_CommonDifExpMat_16_ERassoc[aa3, ]
plotExpression(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, filename = "genes_gte4nonzeroKD_atleast1p1n_heatmap.png",.dendrogram = "both", .Colv = T)


plotExpression(my_CommonDifExpMat_16_ERassoc[aa4p >= 1 & aa4n >= 1, ], filename = "genes_atleast1p1n_heatmap.png")
plotExpression(my_CommonDifExpMat_16_ERassoc[aa4p > 1 & aa4n > 1, ], filename = "genes_atleast2p2n_heatmap.png")

############################################################################################################
############################################################################################################
########################################                               #####################################
########################################   Descritized TF expression   #####################################
########################################                               #####################################
############################################################################################################
############################################################################################################
# TF.motifs.Shrinked
# TF.ENTREZ.Shrinked$ENTREZID is the ENTREZ ID for the TFs we are interested in
# should find the TF's index in each experiment's epressionset. 
# assign 0 / 0.33 / 0.66 / 1 to the TF in each condition based on the quantiles of the control experiment
# For each condition in "my_CommonDifExpMat_16" I will have two values for the TF expression corresponding to the contrasting conditions.

length(my_Descretized_List) # GSE26740 is the last entry in this list, but its composition is different from other entries
#find the TFs in each expressionset

################ find entrez ID matches in microarray data:
TF.Index.Shrinked.microarray_ENTREZ <- list()
for(i in 1:length(TF.ENTREZ.Shrinked$ENTREZID)){
  TF.Index.Shrinked.microarray_ENTREZ[[i]] <- list()
  for(j in 1:length(my_Descretized_List)){
    TF.Index.Shrinked.microarray_ENTREZ[[i]][[j]] <- which(my_Descretized_List[[j]]$fDATA$ENTREZID %in%  TF.ENTREZ.Shrinked$ENTREZID[i])
  }
  names(TF.Index.Shrinked.microarray_ENTREZ[[i]]) <- names(my_Descretized_List)
}
names(TF.Index.Shrinked.microarray_ENTREZ) <- TF.ENTREZ.Shrinked$ALIAS
################ find entrez ID matches in microarray data -- added E2F1:
TF.ENTREZ.Shrinked_2 <- TF.ENTREZ.Shrinked
TF.ENTREZ.Shrinked_2 <- rbind(TF.ENTREZ.Shrinked_2[1:2, ], 
                              c("E2F1", "1869"),
                              TF.ENTREZ.Shrinked_2[3:nrow(TF.ENTREZ.Shrinked_2), ])
rownames(TF.ENTREZ.Shrinked_2) <- c(1:nrow(TF.ENTREZ.Shrinked_2))

TF.Index.Shrinked.microarray_ENTREZ_2 <- list()
for(i in 1:length(TF.ENTREZ.Shrinked_2$ENTREZID)){
  TF.Index.Shrinked.microarray_ENTREZ_2[[i]] <- list()
  for(j in 1:length(my_Descretized_List)){
    TF.Index.Shrinked.microarray_ENTREZ_2[[i]][[j]] <- which(my_Descretized_List[[j]]$fDATA$ENTREZID %in%  TF.ENTREZ.Shrinked_2$ENTREZID[i])
  }
  names(TF.Index.Shrinked.microarray_ENTREZ_2[[i]]) <- names(my_Descretized_List)
}
names(TF.Index.Shrinked.microarray_ENTREZ_2) <- TF.ENTREZ.Shrinked_2$ALIAS

################ find entrez ID matches for RNA-seq data:
TF.Index.Shrinked.RNAseq_ENTREZ <- list()
for(i in 1:length(TF.ENTREZ.Shrinked$ENTREZID)){
  TF.Index.Shrinked.RNAseq_ENTREZ[[i]] <- list()
  for(j in 1:length(my_Descretized_List_RNAseq)){
    TF.Index.Shrinked.RNAseq_ENTREZ[[i]][[j]] <- which(rownames(my_Descretized_List_RNAseq[[j]]$descrete_expMat) %in%  TF.ENTREZ.Shrinked$ENTREZID[i])
  }
  names(TF.Index.Shrinked.RNAseq_ENTREZ[[i]]) <- names(my_Descretized_List_RNAseq)
}
names(TF.Index.Shrinked.RNAseq_ENTREZ) <- TF.ENTREZ.Shrinked$ALIAS
### FOXA1 is not present in GSE55922, I'll set the index to 1. and then manually modify it.
TF.Index.Shrinked.RNAseq_ENTREZ$FOXA1$GSE55922 <- 1 

################ find entrez ID matches for RNA-seq data -- added E2F1
TF.Index.Shrinked.RNAseq_ENTREZ_2 <- list()
for(i in 1:length(TF.ENTREZ.Shrinked_2$ENTREZID)){
  TF.Index.Shrinked.RNAseq_ENTREZ_2[[i]] <- list()
  for(j in 1:length(my_Descretized_List_RNAseq)){
    TF.Index.Shrinked.RNAseq_ENTREZ_2[[i]][[j]] <- which(rownames(my_Descretized_List_RNAseq[[j]]$descrete_expMat) %in%  TF.ENTREZ.Shrinked_2$ENTREZID[i])
  }
  names(TF.Index.Shrinked.RNAseq_ENTREZ_2[[i]]) <- names(my_Descretized_List_RNAseq)
}
names(TF.Index.Shrinked.RNAseq_ENTREZ_2) <- TF.ENTREZ.Shrinked_2$ALIAS
### FOXA1 is not present in GSE55922, I'll set the index to 1. and then manually modify it.
TF.Index.Shrinked.RNAseq_ENTREZ_2$FOXA1$GSE55922 <- 1 

##################################################################################################################
TF.Index.Shrinked.microarray_ENTREZ #contains the indices for each TF in each microarray dataset
TF.Index.Shrinked.microarray_ENTREZ_2 # contains the indices for each TF in each microarray dataset after adding E2F1
TF.Index.Shrinked.RNAseq_ENTREZ #contains the indices for each TF in each RNAseq dataset
TF.Index.Shrinked.RNAseq_ENTREZ_2 #contains the indices for each TF in each RNAseq dataset after adding E2F1
##################################################################################################################
# choose one or average for the TFs with more than one index ???
# create the TF expression matrix. in order to do so first:
# create the quantile expression vector for each probe corresponding to a TF:
# Given all expression values in a condition, which quantile of the control experiment does the probe's expression belong to

# create an index list per dataset for microarray data
TF.Index.Shrinked.microarray_ENTREZ_PerDataset <- list()
for(i in 1:length(my_Descretized_List)){
  TF.Index.Shrinked.microarray_ENTREZ_PerDataset[[i]] <- lapply(TF.Index.Shrinked.microarray_ENTREZ, "[[", i)
  names(TF.Index.Shrinked.microarray_ENTREZ_PerDataset[[i]]) <- names(TF.Index.Shrinked.microarray_ENTREZ)
  for(j in 1:length(TF.Index.Shrinked.microarray_ENTREZ_PerDataset[[i]])){
    names(TF.Index.Shrinked.microarray_ENTREZ_PerDataset[[i]][[j]]) <- c(1:length(TF.Index.Shrinked.microarray_ENTREZ_PerDataset[[i]][[j]]))
  }
}
names(TF.Index.Shrinked.microarray_ENTREZ_PerDataset) <- names(my_Descretized_List)
# create an index list per dataset for microarray data after adding E2F1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2 <- list()
for(i in 1:length(my_Descretized_List)){
  TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[i]] <- lapply(TF.Index.Shrinked.microarray_ENTREZ_2, "[[", i)
  names(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[i]]) <- names(TF.Index.Shrinked.microarray_ENTREZ_2)
  for(j in 1:length(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[i]])){
    #if(length(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[i]][[j]]) > 0){
      names(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[i]][[j]]) <- c(1:length(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[i]][[j]]))
    #}
  }
}
names(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2) <- names(my_Descretized_List)
for(i in 1:length(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2)){
  print(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[i]][[3]])
}

# create an index list per dataset for RNAseq data
TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset <- list()
for(i in 1:length(my_Descretized_List_RNAseq)){
  TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset[[i]] <- lapply(TF.Index.Shrinked.RNAseq_ENTREZ, "[[", i)
  names(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset[[i]]) <- names(TF.Index.Shrinked.RNAseq_ENTREZ)
  for(j in 1:length(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset[[i]])){
    names(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset[[i]][[j]]) <- c(1:length(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset[[i]][[j]]))
  }
}
names(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset) <- names(my_Descretized_List_RNAseq)

# create an index list per dataset for RNAseq data -- after adding E2F1
TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset_2 <- list()
for(i in 1:length(my_Descretized_List_RNAseq)){
  TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset_2[[i]] <- lapply(TF.Index.Shrinked.RNAseq_ENTREZ_2, "[[", i)
  names(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset_2[[i]]) <- names(TF.Index.Shrinked.RNAseq_ENTREZ_2)
  for(j in 1:length(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset_2[[i]])){
    names(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset_2[[i]][[j]]) <- c(1:length(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset_2[[i]][[j]]))
  }
}
names(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset_2) <- names(my_Descretized_List_RNAseq)

# create TF expression matrix per microarray dataset --> this is except the last dataset, I have to add that manually
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset <- list()
for(i in 1:(length(my_Descretized_List)-1)){
  aaa <- unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset[[i]])
  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset[[i]] <- TFexpressionDescritizer(createDescreteMat_output = my_Descretized_List[[i]], TF_index = aaa)
}
names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset) <- names(my_Descretized_List)[1:(length(my_Descretized_List)-1)]

# create TF expression matrix per microarray dataset -- after adding E2F1 --> this is except the last dataset, I have to add that manually
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2 <- list()
for(i in 1:(length(my_Descretized_List)-1)){
  aaa <- unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[i]])
  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2[[i]] <- TFexpressionDescritizer(createDescreteMat_output = my_Descretized_List[[i]], TF_index = aaa)
}
names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2) <- names(my_Descretized_List)[1:(length(my_Descretized_List)-1)]


# create TF expression matrix per RNAseq dataset --> remember to fix the values of FOXA1 in GSE55922 dataset
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset <- list()
for(i in 1:(length(my_Descretized_List_RNAseq))){
  aaa <- unlist(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset[[i]])
  TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset[[i]] <- TFexpressionDescritizer(createDescreteMat_output = my_Descretized_List_RNAseq[[i]], TF_index = aaa)
}
names(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset) <- names(my_Descretized_List_RNAseq)
# fix the values of FOXA1 in GSE55922 dataset to 1 since its one in all other exps
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset$GSE55922[4, ] <- rep(1, ncol(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset$GSE55922))


# create TF expression matrix per RNAseq dataset -- after adding E2F1 --> remember to fix the values of FOXA1 in GSE55922 dataset
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_2 <- list()
for(i in 1:(length(my_Descretized_List_RNAseq))){
  aaa <- unlist(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset_2[[i]])
  TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_2[[i]] <- TFexpressionDescritizer(createDescreteMat_output = my_Descretized_List_RNAseq[[i]], TF_index = aaa)
}
names(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_2) <- names(my_Descretized_List_RNAseq)
# fix the values of FOXA1 in GSE55922 dataset to 1 since its one in all other exps
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_2$GSE55922[4, ] <- rep(1, ncol(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_2$GSE55922))

###############################
# Create the discretized TF expression profile for the last data set : GSE26740 and add it to the current set
aaContWT <- my_Descretized_List[[16]]$FullexpMat[, 1]
aaContWT_qu <- quantile(aaContWT)
aacontKD <- rowMeans(my_Descretized_List[[16]]$FullexpMat[, c(3, 5)])
aacontKD_qu <- quantile(aacontKD)
aatreatWT <-  my_Descretized_List[[16]]$FullexpMat[, 2]
aatreatKD <- rowMeans(my_Descretized_List[[16]]$FullexpMat[, c(4, 6)])
aaContWT_dis <- QuantileCompare(exp_vec = aaContWT, quantile_vec = aaContWT_qu, assignedVal = c(0, 0.33, 0.66, 1))
aaContKD_dis <- QuantileCompare(exp_vec = aacontKD, quantile_vec = aacontKD_qu, assignedVal = c(0, 0.33, 0.66, 1))
aatreatWT_dis <- QuantileCompare(exp_vec = aatreatWT, quantile_vec = aaContWT_qu, assignedVal = c(0, 0.33, 0.66, 1))
aatreatKD_dis <- QuantileCompare(exp_vec = aatreatKD, quantile_vec = aacontKD_qu, assignedVal = c(0, 0.33, 0.66, 1))
aaAll_dis <- cbind(aaContWT_dis, aatreatWT_dis, aaContKD_dis, aatreatKD_dis)
colnames(aaAll_dis) <- c("WT_control", "WT_E2_12h", "siAP2g_control", "siAP2g_E2_12h")
aaTFdis <- aaAll_dis[unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset[[16]]), ]
rownames(aaTFdis) <- names(unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset[[16]]))
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset[[16]] <- aaTFdis
names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset)[16] <- names(my_Descretized_List)[16]

# Create the discretized TF expression profile for the last data set -- after adding E2F1 : GSE26740 and add it to the current set
aaContWT <- my_Descretized_List[[16]]$FullexpMat[, 1]
aaContWT_qu <- quantile(aaContWT)
aacontKD <- rowMeans(my_Descretized_List[[16]]$FullexpMat[, c(3, 5)])
aacontKD_qu <- quantile(aacontKD)
aatreatWT <-  my_Descretized_List[[16]]$FullexpMat[, 2]
aatreatKD <- rowMeans(my_Descretized_List[[16]]$FullexpMat[, c(4, 6)])
aaContWT_dis <- QuantileCompare(exp_vec = aaContWT, quantile_vec = aaContWT_qu, assignedVal = c(0, 0.33, 0.66, 1))
aaContKD_dis <- QuantileCompare(exp_vec = aacontKD, quantile_vec = aacontKD_qu, assignedVal = c(0, 0.33, 0.66, 1))
aatreatWT_dis <- QuantileCompare(exp_vec = aatreatWT, quantile_vec = aaContWT_qu, assignedVal = c(0, 0.33, 0.66, 1))
aatreatKD_dis <- QuantileCompare(exp_vec = aatreatKD, quantile_vec = aacontKD_qu, assignedVal = c(0, 0.33, 0.66, 1))
aaAll_dis <- cbind(aaContWT_dis, aatreatWT_dis, aaContKD_dis, aatreatKD_dis)
colnames(aaAll_dis) <- c("WT_control", "WT_E2_12h", "siAP2g_control", "siAP2g_E2_12h")
aaTFdis <- aaAll_dis[unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[16]]), ]
rownames(aaTFdis) <- names(unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[16]]))
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2[[16]] <- aaTFdis
names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2)[16] <- names(my_Descretized_List)[16]

##################################################################################################################
# plot TF expression mats
for(i in 1:length(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset)){
  plotExpression(expMat =  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset[[i]], colorPoints = c(0 , 0.33, 0.66, 1), filename = paste("TF_exp_", names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset)[i], ".png",sep = ""))
}

##################################################################################################################
# GO one by one and choose one probe, or median for each TF in each dataset
# plot without ordering 
for(i in 1:length(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset)){
  plotExpression(expMat =  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset[[i]],
                 colorPoints = c(0 , 0.33, 1),
                 filename = paste("TF_exp_", names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset)[i], ".png",sep = ""), .dendrogram = "none", .Rowv = F)
  
}
# This is the list containing one index per TF in each dataset
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset
# This is the list containing one index per TF in each dataset -- after adding E2F1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2

##################################################################################################################
# 1 : GSE37386
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE37386
#AR.1
#NR3C1.2
#NR5A2.1
#RARA.2
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE37386$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE37386$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE37386$NR3C1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE37386$NR3C1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE37386$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE37386$NR5A2[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE37386$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE37386$RARA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE37386$TFAP2C <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE37386$TFAP2C[1]

# after adding E2F1
# 1 : GSE37386
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE37386
#AR.1
#NR3C1.2
#NR5A2.1
#RARA.2
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE37386$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE37386$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE37386$NR3C1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE37386$NR3C1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE37386$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE37386$NR5A2[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE37386$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE37386$RARA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE37386$TFAP2C <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE37386$TFAP2C[1]

heatmap.2(aa_cor_lfc$GSE37386, Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc$GSE37386), na.rm = T)
aa <- which(rownames(aa_cor_lfc$GSE37386) %in% c("AR.2", "NR3C1.1", "NR5A2.2", "NR5A2.3", "RARA.1", "RARA.3", "TFAP2C.2" )) 
aa_cor_lfc$GSE37386 <- aa_cor_lfc$GSE37386[-aa, ]

aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[1]]) %in% c("AR.2", "NR3C1.1", "NR5A2.2", "NR5A2.3", "RARA.1", "RARA.3", "TFAP2C.2" )) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[1]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[1]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[1]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[1]]) %in% c("AR.2", "NR3C1.1", "NR5A2.2", "NR5A2.3", "RARA.1", "RARA.3", "TFAP2C.2" )) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[1]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[1]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[1]])
##################################################################################################################
# 2 : GSE25314
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE25314
#AR.2
#LEF-1.1
#NR3C1.1
#NR5A2.3
#RARA.1
#RUNX1.2
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$AR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$`LEF-1` <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$`LEF-1`[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$NR3C1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$RUNX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$RUNX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$TFAP2C <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25314$TFAP2C[1]
# after adding E2F1:
# 2 : GSE25314
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE25314
#AR.2
#LEF-1.1
#NR3C1.1
#NR5A2.3
#RARA.1
#RUNX1.2
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$AR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$`LEF-1` <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$`LEF-1`[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$NR3C1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$RUNX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$RUNX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$TFAP2C <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25314$TFAP2C[1]

heatmap.2(aa_cor_lfc$GSE25314, Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc$GSE25314), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[2]]) %in% c("AR.1", "AR.3", "LEF-1.2", "NR3C1.2", "NR5A2.1", "NR5A2.2", "RARA.2", "RUNX1.1", "TFAP2C.2" )) 
aa_cor_lfc[[2]] <- aa_cor_lfc[[2]][-aa, ]

i <- 2
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% c("AR.1", "AR.3", "LEF-1.2", "NR3C1.2", "NR5A2.1", "NR5A2.2", "RARA.2", "RUNX1.1", "TFAP2C.2")) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% c("AR.1", "AR.3", "LEF-1.2", "NR3C1.2", "NR5A2.1", "NR5A2.2", "RARA.2", "RUNX1.1", "TFAP2C.2")) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])

##################################################################################################################
# 3 : GSE25315
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE25315
#AR.2
#NR3C1.2
#NR5A2.3
#RARA.1
#RUNX1.3
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25315$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25315$AR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25315$NR3C1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25315$NR3C1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25315$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25315$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25315$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25315$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25315$RUNX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE25315$RUNX1[3]

# after adding E2F1:
# 3 : GSE25315
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE25315
#AR.2
#NR3C1.2
#NR5A2.3
#RARA.1
#RUNX1.3
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25315$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25315$AR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25315$NR3C1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25315$NR3C1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25315$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25315$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25315$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25315$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25315$RUNX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE25315$RUNX1[3]

i <- 3
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.1", "NR3C1.1", "NR5A2.1", "NR5A2.2", "RARA.2", "RARA.3", "RUNX1.1", "RUNX1.2")) 
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][-aa, ]

i <- 3
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% c("AR.1", "NR3C1.1", "NR5A2.1", "NR5A2.2", "RARA.2", "RARA.3", "RUNX1.1", "RUNX1.2")) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% c("AR.1", "NR3C1.1", "NR5A2.1", "NR5A2.2", "RARA.2", "RARA.3", "RUNX1.1", "RUNX1.2")) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])

##################################################################################################################
# 4 : GSE102367
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE102367
# as is
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE102367

i <- 4
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
#aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.1", "NR3C1.1", "NR5A2.1", "NR5A2.2", "RARA.2", "RARA.3", "RUNX1.1", "RUNX1.2")) 
#aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][-aa, ]

##################################################################################################################
# 5 : GSE79761
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE79761
#AR.3
#ESR1.1
#FOXA1.1
#GATA3.1
#JUN.1
#LEF-1.3
#NKX3-1.1
#NR3C1.1
#NR5A2.4
#PBX1.2
#PGR.2
#RARA.1
#RUNX1.3
#RXRA.2
#SP1.3
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$AR[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$ESR1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$FOXA1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$FOXA1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$GATA3 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$JUN <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$JUN[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$`LEF-1` <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$`LEF-1`[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$`NKX3-1` <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$`NKX3-1`[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$NR3C1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$NR5A2[4]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$PBX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$PGR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$PGR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$RUNX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$RUNX1[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$RXRA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$SP1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$SP1[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$TFAP2C <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE79761$TFAP2C[1]

# after adding E2F1:
# 5 : GSE79761
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE79761
#AR.3
#ESR1.1
#FOXA1.1
#GATA3.1
#JUN.1
#LEF-1.3
#NKX3-1.1
#NR3C1.1
#NR5A2.4
#PBX1.2
#PGR.2
#RARA.1
#RUNX1.3
#RXRA.2
#SP1.3
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$AR       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$AR[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$E2F1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$E2F1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$ESR1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$FOXA1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$FOXA1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$GATA3    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$JUN      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$JUN[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$`LEF-1`  <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$`LEF-1`[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$`NKX3-1` <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$`NKX3-1`[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$NR3C1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$NR5A2    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$NR5A2[4]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$PBX1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$PGR      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$PGR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$RARA     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$RUNX1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$RUNX1[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$RXRA     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$SP1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$SP1[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$TFAP2C   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE79761$TFAP2C[1]

i <- 5
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.1", "AR.2", "AR.4", "E2F1.2","ESR1.3","ESR1.4","ESR1.5","ESR1.6","ESR1.7","ESR1.8","ESR1.2","FOXA1.2", 
                                             "GATA3.2", "GATA3.3", "JUN.2", "JUN.3", "JUN.4","LEF-1.1", "LEF-1.2", "NKX3-1.2", "NKX3-1.3",
                                             "NR3C1.2", "NR3C1.3", "NR3C1.4", "NR3C1.5","NR5A2.1", "NR5A2.2", "NR5A2.3","PBX1.1", "PBX1.3", "PGR.1",
                                             "RARA.2", "RARA.3", "RARA.4","RUNX1.1", "RUNX1.2", "RUNX1.4","RUNX1.5","RUNX1.6","RUNX1.7","RUNX1.8","RUNX1.9",
                                             "RXRA.1", "SP1.1", "SP1.2", 'SP1.4', "TFAP2C.2")) 
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][-aa, ]
aaqq <- c("AR.1", "AR.2", "AR.4", "E2F1.2","ESR1.3","ESR1.4","ESR1.5","ESR1.6","ESR1.7","ESR1.8","ESR1.2","FOXA1.2", 
          "GATA3.2", "GATA3.3", "JUN.2", "JUN.3", "JUN.4","LEF-1.1", "LEF-1.2", "NKX3-1.2", "NKX3-1.3",
          "NR3C1.2", "NR3C1.3", "NR3C1.4", "NR3C1.5","NR5A2.1", "NR5A2.2", "NR5A2.3","PBX1.1", "PBX1.3", "PGR.1",
          "RARA.2", "RARA.3", "RARA.4","RUNX1.1", "RUNX1.2", "RUNX1.4","RUNX1.5","RUNX1.6","RUNX1.7","RUNX1.8","RUNX1.9",
          "RXRA.1", "SP1.1", "SP1.2", 'SP1.4', "TFAP2C.2")
i <- 5
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])


##################################################################################################################
# 6 : GSE57935
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE57935

#ESR1.2
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE57935$ESR1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE57935$ESR1[2]

# after adding E2F1:
# 6 : GSE57935
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE57935

#ESR1.2
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE57935$ESR1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE57935$ESR1[2]

i <- 6
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("ESR1.1")) 
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][-aa, ]
nrow(aa_cor_lfc[[i]])

i <- 6
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% c("ESR1.1")) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% c("ESR1.1")) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])

##################################################################################################################
# 7 : GSE68918
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE68918
#AR.2
#LEF-1.1
#NR3C1.1
#NR5A2.2
#RARA.1
#RUNX1.2
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$AR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$`LEF-1` <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$`LEF-1`[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$NR3C1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$NR5A2[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$RUNX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$RUNX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$TFAP2C <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE68918$TFAP2C[1]

# after adding E2F1:
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE68918
#AR.2
#LEF-1.1
#NR3C1.1
#NR5A2.2
#RARA.1
#RUNX1.2
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$AR         <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$AR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$`LEF-1`    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$`LEF-1`[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$NR3C1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$NR5A2      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$NR5A2[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$RARA       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$RUNX1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$RUNX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$TFAP2C     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE68918$TFAP2C[1]

i <- 7
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.1", "AR.3","LEF-1.2", "NR3C1.2", "NR3C1.3",
                                             "NR5A2.1", "NR5A2.3", "RARA.2",
                                             "RARA.3", "RUNX1.1", "TFAP2C.2")) 
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][-aa, ]
nrow(aa_cor_lfc[[i]])

aaqq <- c("AR.1", "AR.3","LEF-1.2", "NR3C1.2", "NR3C1.3",
          "NR5A2.1", "NR5A2.3", "RARA.2",
          "RARA.3", "RUNX1.1", "TFAP2C.2")
i <- 7
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][-aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])

##################################################################################################################
# 8 : GSE24592
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE24592
#AR.1
#ESR1.1
#GATA3.1
#JUN.1
#LEF-1.3
#NR3C1.3
#NR5A2.3
#PBX1.2
#RARA.1
#RARG.2
#RUNX1.3
#RXRA.2
#TFAP2C.1

TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$AR       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$ESR1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$GATA3    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$JUN      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$JUN[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$`LEF-1`  <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$`LEF-1`[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$NR3C1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$NR3C1[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$NR5A2    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$PBX1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$RARA     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$RARG     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$RARG[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$RUNX1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$RUNX1[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$RXRA     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$TFAP2C   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE24592$TFAP2C[1]
# after adding E2F1:
# 8 : GSE24592
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE24592
#AR.1
#ESR1.1
#GATA3.1
#JUN.1
#LEF-1.3
#NR3C1.3
#NR5A2.3
#PBX1.2
#RARA.1
#RARG.2
#RUNX1.3
#RXRA.2
#TFAP2C.1

TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$AR       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$E2F1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$E2F1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$ESR1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$GATA3    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$JUN      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2$GSE24592$JUN[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$`LEF-1`  <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$`LEF-1`[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$NR3C1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$NR3C1[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$NR5A2    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$PBX1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$RARA     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$RARG     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$RARG[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$RUNX1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2$GSE24592$RUNX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$RXRA     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE24592$TFAP2C   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2$GSE24592$TFAP2C[2]

i <- 8
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.1","CEBPB.1","E2F1.2","ESR1.1","FOXA1.1","GATA3.1","JUN.3","LEF-1.3", "NKX3-1.1",
                                             "NR3C1.3", "NR5A2.3", "PBX1.2", "PGR.1","RARA.1", "RARG.2", "RUNX1.2", "RXRA.2", "SP1.1","TFAP2C.2")) 
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][aa, ]
nrow(aa_cor_lfc[[i]])
aaqq <- c("AR.1","CEBPB.1","E2F1.2","ESR1.1","FOXA1.1","GATA3.1","JUN.3","LEF-1.3", "NKX3-1.1",
          "NR3C1.3", "NR5A2.3", "PBX1.2", "PGR.1","RARA.1", "RARG.2", "RUNX1.2", "RXRA.2", "SP1.1","TFAP2C.2")
i <- 8
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])

##################################################################################################################
# 9 : GSE53394
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE53394
#AR.1
#ESR1.1
#GATA3.1
#JUN.1
#LEF-1.2
#NR3C1.1
#NR5A2.3
#PBX1.2
#RARA.1
#RUNX1.3
#RXRA.2
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$ESR1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$GATA3 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$JUN <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$JUN[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$`LEF-1` <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$`LEF-1`[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$NR3C1 <-  TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$PBX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$RUNX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$RUNX1[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$RXRA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$TFAP2C <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE53394$TFAP2C[1]
# after adding E2F1:
# 9 : GSE53394
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE53394
#AR.1
#ESR1.1
#GATA3.1
#JUN.1
#LEF-1.2
#NR3C1.1
#NR5A2.3
#PBX1.2
#RARA.1
#RUNX1.3
#RXRA.2
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$AR         <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$E2F1       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$E2F1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$ESR1       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$GATA3      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$JUN        <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$JUN[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$`LEF-1`    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$`LEF-1`[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$NR3C1      <-  TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$NR5A2      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$PBX1       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$RARA       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$RUNX1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$RUNX1[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$RXRA       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$TFAP2C     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE53394$TFAP2C[1]

i <- 9
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.1","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.2", "NKX3-1.1",
                                             "NR3C1.1", "NR5A2.3", "PBX1.2", "PGR.1","RARA.1", "RARG.1", "RUNX1.3", "RXRA.2", "SP1.1","TFAP2C.1")) 
length(aa)
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][aa, ]
nrow(aa_cor_lfc[[i]])

aaqq <- c("AR.1","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.2", "NKX3-1.1",
           "NR3C1.1", "NR5A2.3", "PBX1.2", "PGR.1","RARA.1", "RARG.1", "RUNX1.3", "RXRA.2", "SP1.1","TFAP2C.1")
i <- 9
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])


##################################################################################################################
# 10 : GSE8597
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE8597
#AR.1
#ESR1.1
#FOXA1.1
#GATA3.1
#JUN.1
#NR3C1.1
#NR5A2.1
#PBX1.2
#PGR.2
#RARA.3
#RARG.1
#RUNX1.2
#RXRA.2
#SP1.1
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$AR      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$ESR1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$FOXA1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$FOXA1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$GATA3   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$JUN     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$JUN[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$NR3C1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$NR5A2   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$NR5A2[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$PBX1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$PGR     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$PGR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$RARA    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$RARA[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$RARG    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$RARG[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$RUNX1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$RUNX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$RXRA    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$SP1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$SP1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$TFAP2C  <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE8597$TFAP2C[1]

# after adding E2F1:
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE8597
#AR.1
#ESR1.1
#FOXA1.1
#GATA3.1
#JUN.1
#NR3C1.1
#NR5A2.1
#PBX1.2
#PGR.2
#RARA.3
#RARG.1
#RUNX1.2
#RXRA.2
#SP1.1
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$AR      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$E2F1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$E2F1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$ESR1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$FOXA1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$FOXA1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$GATA3   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$JUN     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$JUN[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$NR3C1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$NR5A2   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$NR5A2[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$PBX1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$PGR     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$PGR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$RARA    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$RARA[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$RARG    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$RARG[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$RUNX1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$RUNX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$RXRA    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$SP1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$SP1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$TFAP2C  <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE8597$TFAP2C[1]

i <- 10
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.1","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
                                             "NR3C1.1", "NR5A2.1", "PBX1.2", "PGR.2","RARA.3", "RARG.1", "RUNX1.2", "RXRA.2", "SP1.1","TFAP2C.1")) 
length(aa)
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][aa, ]
nrow(aa_cor_lfc[[i]])

aaqq <- c("AR.1","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
          "NR3C1.1", "NR5A2.1", "PBX1.2", "PGR.2","RARA.3", "RARG.1", "RUNX1.2", "RXRA.2", "SP1.1","TFAP2C.1")
i <- 10
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])


##################################################################################################################
# 11 : GSE4006
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE4006
#AR.1
#ESR1.1
#GATA3.1
#JUN.1
#NR3C1.1
#NR5A2.1
#PBX1.2
#RARA.1
#RUNX1.1
#RXRA.2
#SP1.1
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$AR         <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$ESR1       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$GATA3      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$JUN        <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$JUN[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$NR3C1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$NR5A2      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$NR5A2[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$PBX1       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$RARA       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$RUNX1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$RUNX1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$RXRA       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$SP1        <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$SP1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$TFAP2C     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE4006$TFAP2C[1]

# after adding E2F1:
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE4006
#AR.1
#ESR1.1
#GATA3.1
#JUN.1
#NR3C1.1
#NR5A2.1
#PBX1.2
#RARA.1
#RUNX1.1
#RXRA.2
#SP1.1
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$AR         <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$ESR1       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$GATA3      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$JUN        <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$JUN[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$NR3C1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$NR5A2      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$NR5A2[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$PBX1       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$RARA       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$RUNX1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$RUNX1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$RXRA       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$SP1        <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$SP1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$TFAP2C     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE4006$TFAP2C[1]

i <- 11
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.1","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
                                             "NR3C1.1", "NR5A2.1", "PBX1.2", "PGR.1","RARA.1", "RARG.1", "RUNX1.1", "RXRA.2", "SP1.1","TFAP2C.1")) 
length(aa)
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][aa, ]
nrow(aa_cor_lfc[[i]])

aaqq <- c("AR.1","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
          "NR3C1.1", "NR5A2.1", "PBX1.2", "PGR.1","RARA.1", "RARG.1", "RUNX1.1", "RXRA.2", "SP1.1","TFAP2C.1")
i <- 11
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])


##################################################################################################################
# 12 : GSE74146
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE74146
#AR.1
#LEF-1.1
#NR3C1.1
#NR5A2.3
#RARA.1
#RUNX1.2
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$`LEF-1` <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$`LEF-1`[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$NR3C1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$RUNX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE74146$RUNX1[2]

# after adding E2F1:
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE74146
#AR.1
#LEF-1.1
#NR3C1.1
#NR5A2.3
#RARA.1
#RUNX1.2
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$AR      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$AR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$`LEF-1` <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$`LEF-1`[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$NR3C1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$NR5A2   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$RARA    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$RARA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$RUNX1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE74146$RUNX1[2]

i <- 12
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.1","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
                                             "NR3C1.1", "NR5A2.3", "PBX1.1", "PGR.1","RARA.1", "RARG.1", "RUNX1.2", "RXRA.1", "SP1.1","TFAP2C.1")) 
length(aa)
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][aa, ]
nrow(aa_cor_lfc[[i]])

aaqq <- c("AR.1","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
          "NR3C1.1", "NR5A2.3", "PBX1.1", "PGR.1","RARA.1", "RARG.1", "RUNX1.2", "RXRA.1", "SP1.1","TFAP2C.1")
i <- 12
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])

##################################################################################################################
# 13 : GSE56245
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE56245
#CEBPB.2
#ESR1.1
#FOXA1.1
#LEF-1.2
#NR3C1.1
#NR5A2.3
#PBX1.2
#PGR.2
#RARA.2
#RARG.2
#RUNX1.2
#RXRA.1
#SP1.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$CEBPB     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$CEBPB[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$ESR1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$FOXA1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$FOXA1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$`LEF-1`   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$`LEF-1`[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$NR3C1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$NR5A2     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$PBX1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$PGR       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$PGR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$RARA      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$RARA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$RARG      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$RARG[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$RUNX1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$RUNX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$RXRA      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$RXRA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$SP1       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE56245$SP1[1]

# after adding E2F1:
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE56245
#CEBPB.2
#ESR1.1
#FOXA1.1
#LEF-1.2
#NR3C1.1
#NR5A2.3
#PBX1.2
#PGR.2
#RARA.2
#RARG.2
#RUNX1.2
#RXRA.1
#SP1.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$CEBPB     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$CEBPB[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$E2F1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$E2F1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$ESR1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$FOXA1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$FOXA1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$`LEF-1`   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$`LEF-1`[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$NR3C1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$NR5A2     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$NR5A2[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$PBX1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$PGR       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$PGR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$RARA      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$RARA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$RARG      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$RARG[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$RUNX1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$RUNX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$RXRA      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$RXRA[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$SP1       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE56245$SP1[1]

i <- 13
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.1","CEBPB.2","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.2", "NKX3-1.1",
                                             "NR3C1.1", "NR5A2.3", "PBX1.2", "PGR.2","RARA.2", "RARG.2", "RUNX1.2", "RXRA.1", "SP1.1","TFAP2C.1")) 
length(aa)
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][aa, ]
nrow(aa_cor_lfc[[i]])

aaqq <- c("AR.1","CEBPB.2","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.2", "NKX3-1.1",
         "NR3C1.1", "NR5A2.3", "PBX1.2", "PGR.2","RARA.2", "RARG.2", "RUNX1.2", "RXRA.1", "SP1.1","TFAP2C.1")
i <- 13
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])

##################################################################################################################
# 14 : GSE36586
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE36586

#ESR1.1
#FOXA1.1
#NKX3-1.1
#NR3C1.2
#NR5A2.2
#PBX1.2
#PGR.1
#RARA.2
#RUNX1.3

TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$ESR1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$FOXA1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$FOXA1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$`NKX3-1`  <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$`NKX3-1`[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$NR3C1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$NR3C1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$NR5A2     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$NR5A2[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$PBX1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$PGR       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$PGR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$RARA      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$RARA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$RUNX1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE36586$RUNX1[3]

# after adding E2F1:
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE36586

TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$E2F1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$E2F1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$ESR1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$FOXA1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$FOXA1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$`NKX3-1`  <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$`NKX3-1`[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$NR3C1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$NR3C1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$NR5A2     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$NR5A2[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$PBX1      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$PGR       <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$PGR[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$RARA      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$RARA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$RUNX1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE36586$RUNX1[3]

i <- 14
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.1","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
                                             "NR3C1.2", "NR5A2.2", "PBX1.2", "PGR.1","RARA.2", "RARG.1", "RUNX1.3", "RXRA.1", "SP1.1","TFAP2C.1")) 
length(aa)
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][aa, ]
nrow(aa_cor_lfc[[i]])

aaqq <- c("AR.1","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
          "NR3C1.2", "NR5A2.2", "PBX1.2", "PGR.1","RARA.2", "RARG.1", "RUNX1.3", "RXRA.1", "SP1.1","TFAP2C.1")
i <- 14
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])

##################################################################################################################
# 15 : GSE35428
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE35428
#AR.2
#ESR1.1
#FOXA1.1
#GATA3.1
#JUN.1
#NR3C1.1
#PBX1.2
#PGR.2
#RARA.2
#RARG.1
#RUNX1.2
#RXRA.2
#SP1.3
#TFAP2C.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$AR      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$AR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$ESR1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$FOXA1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$FOXA1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$GATA3   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$JUN     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$JUN[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$NR3C1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$PBX1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$PGR     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$PGR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$RARA    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$RARA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$RARG    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$RARG[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$RUNX1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$RUNX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$RXRA    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$SP1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$SP1[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$TFAP2C  <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE35428$TFAP2C[1]

# after adding E2F1:
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE35428

TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$AR      <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$AR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$E2F1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$E2F1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$ESR1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$ESR1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$FOXA1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$FOXA1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$GATA3   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$GATA3[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$JUN     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$JUN[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$NR3C1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$NR3C1[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$PBX1    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$PBX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$PGR     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$PGR[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$RARA    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$RARA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$RARG    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$RARG[1]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$RUNX1   <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$RUNX1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$RXRA    <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$RXRA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$SP1     <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$SP1[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$TFAP2C  <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE35428$TFAP2C[1]

i <- 15
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.2","CEBPB.1","E2F1.2","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
                                             "NR3C1.1", "NR5A2.1", "PBX1.2", "PGR.2","RARA.2", "RARG.1", "RUNX1.2", "RXRA.2", "SP1.3","TFAP2C.1")) 
length(aa)
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][aa, ]
nrow(aa_cor_lfc[[i]])

aaqq <- c("AR.2","CEBPB.1","E2F1.2","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
          "NR3C1.1", "NR5A2.1", "PBX1.2", "PGR.2","RARA.2", "RARG.1", "RUNX1.2", "RXRA.2", "SP1.3","TFAP2C.1")
i <- 15
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])


##################################################################################################################
# 16 : GSE26740
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset$GSE26740
#AR.3
#NR3C1.2
#NR5A2.2
#RARA.2
#RUNX1.1
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE26740$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE26740$AR[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE26740$NR3C1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE26740$NR3C1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE26740$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE26740$NR5A2[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE26740$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE26740$RARA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE26740$RUNX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique$GSE26740$RUNX1[1]

# after adding E2F1:
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_2$GSE26740
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE26740$AR <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE26740$AR[3]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE26740$NR3C1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE26740$NR3C1[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE26740$NR5A2 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE26740$NR5A2[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE26740$RARA <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE26740$RARA[2]
TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE26740$RUNX1 <- TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2$GSE26740$RUNX1[1]

i <- 16
heatmap.2(aa_cor_lfc[[i]], Rowv = F, Colv = F,
          dendrogram = "none", na.rm = T, trace = "none",
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(n=20))
rowMeans(abs(aa_cor_lfc[[i]]), na.rm = T)
aa <- which(rownames(aa_cor_lfc[[i]]) %in% c("AR.3","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
                                             "NR3C1.2", "NR5A2.2", "PBX1.1", "PGR.1","RARA.2", "RARG.1", "RUNX1.1", "RXRA.1", "SP1.1","TFAP2C.1")) 
length(aa)
aa_cor_lfc[[i]] <- aa_cor_lfc[[i]][aa, ]
nrow(aa_cor_lfc[[i]])

aaqq <- c("AR.3","CEBPB.1","E2F1.1","ESR1.1","FOXA1.1","GATA3.1","JUN.1","LEF-1.1", "NKX3-1.1",
          "NR3C1.2", "NR5A2.2", "PBX1.1", "PGR.1","RARA.2", "RARG.1", "RUNX1.1", "RXRA.1", "SP1.1","TFAP2C.1")
i <- 16
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]])
aa <- which(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) %in% aaqq) 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]][aa, ]
nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]])


##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
#Checking if everything is right with the TF profiles
aa <- 0
for(i in 1:length(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique)){
  print(paste("Dataset number", i))
  print("length of TF index:")
  print(unlist(lapply(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique[[i]], length)))
  print("Sum is.na is:")
  print(sum(unlist(lapply(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique[[i]], is.na))))
  aa <- aa + sum(unlist(lapply(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique[[i]], is.na)))
  print(aa)
}

aa <- 0
for(i in 1:length(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2)){
  print(paste("Dataset number", i))
  print("length of TF index:")
  print(unlist(lapply(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2[[i]], length)))
  print("Sum is.na is:")
  print(sum(unlist(lapply(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2[[i]], is.na))))
  aa <- aa + sum(unlist(lapply(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2[[i]], is.na)))
  print(aa)
}
# Create the unique TF expression profiles
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique <- list()
for(i in 1:(length(my_Descretized_List)-1)){
  aaa <- unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique[[i]])
  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique[[i]] <- TFexpressionDescritizer(createDescreteMat_output = my_Descretized_List[[i]], 
                                                                                      TF_index = aaa)
}
names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique) <- names(my_Descretized_List)[1:(length(my_Descretized_List)-1)]

# Create the unique TF expression profiles -- after adding E2F1
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2 <- list()
for(i in 1:(length(my_Descretized_List)-1)){
  aaa <- unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2[[i]])
  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2[[i]] <- TFexpressionDescritizer(createDescreteMat_output = my_Descretized_List[[i]],
                                                                                        TF_index = aaa)
}
names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2) <- names(my_Descretized_List)[1:(length(my_Descretized_List)-1)]
###############################
# Create the discretized unique TF expression profile for the last data set : GSE26740 and add it to the current set
aaContWT <- my_Descretized_List[[16]]$FullexpMat[, 1]
aaContWT_qu <- quantile(aaContWT)
aacontKD <- rowMeans(my_Descretized_List[[16]]$FullexpMat[, c(3, 5)])
aacontKD_qu <- quantile(aacontKD)
aatreatWT <-  my_Descretized_List[[16]]$FullexpMat[, 2]
aatreatKD <- rowMeans(my_Descretized_List[[16]]$FullexpMat[, c(4, 6)])
aaContWT_dis <- QuantileCompare(exp_vec = aaContWT, quantile_vec = aaContWT_qu, assignedVal = c(0, 0.33, 0.66, 1))
aaContKD_dis <- QuantileCompare(exp_vec = aacontKD, quantile_vec = aacontKD_qu, assignedVal = c(0, 0.33, 0.66, 1))
aatreatWT_dis <- QuantileCompare(exp_vec = aatreatWT, quantile_vec = aaContWT_qu, assignedVal = c(0, 0.33, 0.66, 1))
aatreatKD_dis <- QuantileCompare(exp_vec = aatreatKD, quantile_vec = aacontKD_qu, assignedVal = c(0, 0.33, 0.66, 1))
aaAll_dis <- cbind(aaContWT_dis, aatreatWT_dis, aaContKD_dis, aatreatKD_dis)
colnames(aaAll_dis) <- c("WT_control", "WT_E2_12h", "siAP2g_control", "siAP2g_E2_12h")
aaTFdis <- aaAll_dis[unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique[[16]]), ]
rownames(aaTFdis) <- names(unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique[[16]]))
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique[[16]] <- aaTFdis
names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique)[16] <- names(my_Descretized_List)[16]

# Create the discretized unique TF expression profile for the last data set -- after adding E2F1 : GSE26740 and add it to the current set
aaContWT <- my_Descretized_List[[16]]$FullexpMat[, 1]
aaContWT_qu <- quantile(aaContWT)
aacontKD <- rowMeans(my_Descretized_List[[16]]$FullexpMat[, c(3, 5)])
aacontKD_qu <- quantile(aacontKD)
aatreatWT <-  my_Descretized_List[[16]]$FullexpMat[, 2]
aatreatKD <- rowMeans(my_Descretized_List[[16]]$FullexpMat[, c(4, 6)])
aaContWT_dis <- QuantileCompare(exp_vec = aaContWT, quantile_vec = aaContWT_qu, assignedVal = c(0, 0.33, 0.66, 1))
aaContKD_dis <- QuantileCompare(exp_vec = aacontKD, quantile_vec = aacontKD_qu, assignedVal = c(0, 0.33, 0.66, 1))
aatreatWT_dis <- QuantileCompare(exp_vec = aatreatWT, quantile_vec = aaContWT_qu, assignedVal = c(0, 0.33, 0.66, 1))
aatreatKD_dis <- QuantileCompare(exp_vec = aatreatKD, quantile_vec = aacontKD_qu, assignedVal = c(0, 0.33, 0.66, 1))
aaAll_dis <- cbind(aaContWT_dis, aatreatWT_dis, aaContKD_dis, aatreatKD_dis)
colnames(aaAll_dis) <- c("WT_control", "WT_E2_12h", "siAP2g_control", "siAP2g_E2_12h")
aaTFdis <- aaAll_dis[unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2[[16]]), ]
rownames(aaTFdis) <- names(unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2[[16]]))
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2[[16]] <- aaTFdis
names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2)[16] <- names(my_Descretized_List)[16]

# Concatanate the TF expression profile of all datasets
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat <- do.call(cbind, TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique)
rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat) <- names(TF.Index.Shrinked.microarray_ENTREZ)
plotExpression(expMat =  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat, colorPoints = c(0 , 0.33, 0.66, 1), filename = "TF_Exp_concatenated.png")

# Concatanate the TF expression profile of all datasets -- after adding E2F1
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2 <- do.call(cbind, TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2)
rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2) <- names(TF.Index.Shrinked.microarray_ENTREZ_2)

# TF expression profile for RNAseq data
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat <- do.call(cbind, TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset)
rownames(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat) <- names(TF.Index.Shrinked.microarray_ENTREZ)
plotExpression(expMat =  TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat, filename = "TF_Exp_RNAseq_concatenated.png")

# TF expression profile for RNAseq data --after adding E2F1
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2 <- do.call(cbind, TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_2)
rownames(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2) <- names(TF.Index.Shrinked.microarray_ENTREZ_2)
plotExpression(expMat =  TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat, filename = "TF_Exp_RNAseq_concatenated.png")

##################################################################################################################
# create a TF expression matrix that contains two rows correesponding to JUN because there are tow motifs corresponding to it.
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun <- rbind(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat[1:6, ],
                                                                      TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat[6, ],
                                                                      TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat[7:18, ])
rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun) <- names(TF.motifs.Shrinked)
# same for RNAseq
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun <- rbind(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat[1:6, ],
                                                                  TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat[6, ],
                                                                  TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat[7:18, ])

rownames(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun) <- names(TF.motifs.Shrinked)


# changing the ER expression to 0 for controls and 1 for treatments
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01 <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01[3, seq(1, 90, 2)] <- 0
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01[3, seq(2, 90, 2)] <- 1

TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01 <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01[4, seq(1, 90, 2)] <- 0
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01[4, seq(2, 90, 2)] <- 1

# same for RNAseq
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun_ER01 <- TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun_ER01[3, seq(1, 28, 2)] <- 0
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun_ER01[3, seq(2, 28, 2)] <- 1

# same for RNAseq
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_ER01 <- TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_ER01[4, seq(1, 28, 2)] <- 0
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_ER01[4, seq(2, 28, 2)] <- 1


plotExpression(expMat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01,
               exportplot = T,
               filename = "TF_Exp.png",
               .dendrogram = "none",
               .Rowv = F, .Colv = F, colorPoints = c(-1, 0, 0.33, 0.66, 1))


############################################################################################################################
############################################  Extracting raw fold change values for TFs and Genes  #########################
############################################################################################################################
# get the TF raw fold change values in each dataset
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC <- list()
for(i in 1:(length(my_Descretized_List)-1)){
  aaa <- unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[i]])
  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC[[i]] <- TFexpressionDescritizer(createDescreteMat_output = my_Descretized_List[[i]],
                                                                                            TF_index = aaa, return_raw_lfc = T)
  rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC[[i]]$Raw_LFC) <- names(aaa)
}
names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC) <- names(my_Descretized_List)[1:(length(my_Descretized_List)-1)]
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC <- lapply(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC, "[[", 2)
aaa <- unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_2[[16]])
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC[[16]] <- my_Descretized_List[[16]]$raw_lfc[aaa, ]
rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC[[16]]) <- names(aaa)
names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC)[16] <- names(my_Descretized_List)[16]
# TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC contains raw fold changes for TFs per dataset
# get the 172 gene fold change raw values
aa_gname <- rownames(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172)
aa_raw_gfc <- list()
aa_pgeneEntrez_fc <- list()
for(i in 1:length(my_Descretized_List)){
  if(length(as.character(my_Descretized_List[[i]]$fDATA$Entrez_Gene_ID)) > 0){

    aaa <-  which(as.character(my_Descretized_List[[i]]$fDATA$Entrez_Gene_ID) %in% aa_gname)
    names(aaa) <- as.character(my_Descretized_List[[i]]$fDATA$Entrez_Gene_ID)[aaa]
    
  }else{
    aaa <-  which(as.character(my_Descretized_List[[i]]$fDATA$ENTREZID) %in% aa_gname)
    names(aaa) <- as.character(my_Descretized_List[[i]]$fDATA$ENTREZID)[aaa]
  }
aa_raw_gfc[[i]] <- my_Descretized_List[[i]]$raw_lfc[aaa, ]
  rownames(aa_raw_gfc[[i]]) <- names(aaa)
  aa_pgeneEntrez_fc[[i]] <- matrix(nrow = length(unique(names(aaa))),
                                   ncol = ncol(aa_raw_gfc[[i]]))
  colnames(aa_pgeneEntrez_fc[[i]]) <- colnames((aa_raw_gfc[[i]]))
  for(j in 1:length(unique(names(aaa)))){
    aacur_mat <- aa_raw_gfc[[i]][which(rownames(aa_raw_gfc[[i]]) %in% unique(names(aaa))[j]), ]
    
    if(!is.null(nrow(aacur_mat))){
      aacur_mat[, is.na(MyDifExpMatList[[i]][unique(names(aaa))[j], ])] <- NA
      if(sum(!is.na(aacur_mat)) > 0){
        aamx <- which.max(rowSums(abs(aacur_mat), na.rm = T))
        aacur_mat <- aacur_mat[aamx, ]
      }else{
        aacur_mat <- aacur_mat[1, ]
      }
    }else{
      aacur_mat[is.na(MyDifExpMatList[[i]][unique(names(aaa))[j], ])] <- NA
    }
    aa_pgeneEntrez_fc[[i]][j, ] <- aacur_mat
  }
  rownames(aa_pgeneEntrez_fc[[i]]) <- unique(names(aaa))
}
names(aa_pgeneEntrez_fc) <- names(my_Descretized_List)
#### aa_pgeneEntrez_fc contains the raw fold change 
aa_pgeneEntrez_fc_common <- CommonDifExpMat(aa_pgeneEntrez_fc)
my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_rawLFC <- aa_pgeneEntrez_fc_common
my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_rawFC <- 2 ^ aa_pgeneEntrez_fc_common
nrow(aa_pgeneEntrez_fc_common)
ncol(aa_pgeneEntrez_fc_common)
length(unique(my_Descretized_List[[1]]$fDATA$Entrez_Gene_ID))
# get the 172 gene fold changes descretized
my_CommonDifExpMat_16_ERassoc_gte5nzAll_172

plotExpression(expMat = aa_pgeneEntrez_fc_common, 
               .dendrogram = "none", 
               .Rowv = F,
               filename = "rawexp.png")
plotExpression(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172, 
               .dendrogram = "none", .Rowv = F)
par(mfrow = c(1, 1), mar = c(8,4,4,4))
boxplot.matrix(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172, las= 2)
boxplot.matrix(aa_pgeneEntrez_fc_common, las= 2)
#### calculate the correlation between TF fold change and gene fold change in each of the datasets

TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC <- list()
for(i in 1:length(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC)){
  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC[[i]] <- 2 ^ TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC[[i]]
}

aa_cor_lfc <- list()

for(i in 1:length(aa_pgeneEntrez_fc)){
  aa_cor_lfc[[i]] <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC[[i]]),
                            ncol = nrow(aa_pgeneEntrez_fc[[i]]))
  rownames(aa_cor_lfc[[i]]) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC[[i]])
  colnames(aa_cor_lfc[[i]]) <- rownames(aa_pgeneEntrez_fc[[i]])
  for(j in 1:nrow(aa_cor_lfc[[i]])){
    for(k in 1:ncol(aa_cor_lfc[[i]])){
      aa_cor_lfc[[i]][j, k] <- cor(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC[[i]][j, ], 
                                   aa_pgeneEntrez_fc[[i]][k, ], 
                                   use = "pairwise.complete.obs")
    }
  }
}
names(aa_cor_lfc) <- names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC)
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC

############ IMPORTANT: BEFORE THIS I went through each study and eliminated the extra entries for each TF
for(i in 1:length(aa_cor_lfc)){
  rownames(aa_cor_lfc[[i]]) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
  rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed[[i]]) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
  rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed[[i]]) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
}
############ ############ ############ ############ 

TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common <- CommonDifExpMat(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed)
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed_common <- CommonDifExpMat(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed)

plotExpression(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common, .dendrogram = "none", .Rowv = F)
sum(is.na(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common))

aa_rng <- range(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)
aa_rng <- c(-min(abs(aa_rng)), min(abs(aa_rng)))

heatmap.2(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common,
          Rowv = F, 
          Colv = F, na.rm = T, dendrogram = "none", trace = "none", 
          #breaks = seq(aa_rng[1], aa_rng[2], length.out = 21),
          breaks = seq(-2, 2, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "Raw TF LFC")

aa_rng <- range(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed_common)
#aa_rng <- c(-min(abs(aa_rng)), min(abs(aa_rng)))
heatmap.2(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed_common,
          Rowv = F, 
          Colv = F, na.rm = T, dendrogram = "none", trace = "none", 
          breaks = seq(0, 4, length.out = 21),
          col = colorspace::terrain_hcl(20), main = "Raw TF FC")

aa_rng <- range(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed_common)
#aa_rng <- c(-min(abs(aa_rng)), min(abs(aa_rng)))
heatmap.2(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed_common,
          Rowv = F, 
          Colv = F, na.rm = T, dendrogram = "none", trace = "none", 
          breaks = seq(aa_rng[1], 5, length.out = 101),
          col = colorspace::sequential_hcl(100), main = "TF FC")
############################################################################################################################
#####################################  Raw fold change correlation analysis  ####################################################
############################################################################################################################
#compute the correlation between gene expression (log fold change)or(fold change) values and TFs expression (log fold change)or(fold change) values

aa_tf_gene_lfc_cor <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                             ncol = nrow(aa_pgeneEntrez_fc_common))
rownames(aa_tf_gene_lfc_cor) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)
colnames(aa_tf_gene_lfc_cor) <- rownames(aa_pgeneEntrez_fc_common)

aa_tf_gene_fc_cor <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                             ncol = nrow(aa_pgeneEntrez_fc_common))
rownames(aa_tf_gene_fc_cor) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)
colnames(aa_tf_gene_fc_cor) <- rownames(aa_pgeneEntrez_fc_common)

aafc <- 2 ^ aa_pgeneEntrez_fc_common
for(i in 1:nrow(aa_tf_gene_lfc_cor)){
  for(j in 1:ncol(aa_tf_gene_lfc_cor)){
    aa_tf_gene_lfc_cor[i, j] <- cor(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common[i, ], 
                                    aa_pgeneEntrez_fc_common[j, ], 
                                    use = "pairwise.complete.obs")
    aa_tf_gene_fc_cor[i, j] <- cor(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed_common[i, ], 
                                   aafc[j, ], 
                                    use = "pairwise.complete.obs")
  }
}

heatmap.2(aa_tf_gene_lfc_cor, Rowv = T, 
          Colv = T, na.rm = T, dendrogram = "both", trace = "none", 
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "TF_LFC vs Gene_LFC")

heatmap.2(aa_tf_gene_fc_cor, Rowv = T, 
          Colv = T, na.rm = T, dendrogram = "both", trace = "none", 
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "TF_FC vs Gene_FC")
# compare these to the raw TF expression matrix against discretized exp mat

aa_tf_gene_lfc_cor <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                             ncol = nrow(aa_pgeneEntrez_fc_common))
rownames(aa_tf_gene_lfc_cor) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)
colnames(aa_tf_gene_lfc_cor) <- rownames(aa_pgeneEntrez_fc_common)

aa_tf_gene_fc_cor <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                            ncol = nrow(aa_pgeneEntrez_fc_common))
rownames(aa_tf_gene_fc_cor) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)
colnames(aa_tf_gene_fc_cor) <- rownames(aa_pgeneEntrez_fc_common)

for(i in 1:nrow(aa_tf_gene_lfc_cor)){
  for(j in 1:ncol(aa_tf_gene_lfc_cor)){
    aa_tf_gene_lfc_cor[i, j] <- cor(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common[i, ], 
                                    my_CommonDifExpMat_16_ERassoc_gte5nzAll_172[j, ], 
                                    use = "pairwise.complete.obs")
    aa_tf_gene_fc_cor[i, j] <- cor(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed_common[i, ], 
                                   my_CommonDifExpMat_16_ERassoc_gte5nzAll_172[j, ], 
                                   use = "pairwise.complete.obs")
  }
}

heatmap.2(aa_tf_gene_lfc_cor, Rowv = T, 
          Colv = T, na.rm = T, dendrogram = "both", trace = "none", 
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "TF_LFC vs Gene_discrete")

heatmap.2(aa_tf_gene_fc_cor, Rowv = T, 
          Colv = T, na.rm = T, dendrogram = "both", trace = "none", 
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "TF_FC vs Gene_discrete")

# compare these to the currenlty used TF expression matrix (quantiled) against raw fold changes
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2
aa_tffc <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2),
                  ncol = ncol(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2)/2)
for(i in 1:ncol(aa_tffc)){
  aa_tffc[, i] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2[,2 * i]/TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2[,2 * i - 1]
}
aa_tffc[is.na(aa_tffc)] <- 1
aa_tffc[aa_tffc == Inf] <- 4
rownames(aa_tffc) <- rownames(aa_tf_gene_fc_cor)
colnames(aa_tffc) <- colnames(aa_pgeneEntrez_fc_common)
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_discrete_FC <- aa_tffc
heatmap.2(aa_tffc, Rowv = F, 
          Colv = F, na.rm = T, dendrogram = "none", trace = "none", 
          breaks = seq(0, 4, length.out = 21), margins = c(10,10),
          col = colorspace::terrain_hcl(20), main = "Discrete_TF_FC")

heatmap.2(log2(aa_tffc), Rowv = F, 
          Colv = F, na.rm = T, dendrogram = "none", trace = "none", 
          breaks = seq(-2, 2, length.out = 21), margins = c(10,10),
          col = colorspace::diverge_hcl(20), main = "Discrete_TF_LFC")

aa_tf_gene_lfc_cor <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                             ncol = nrow(aa_pgeneEntrez_fc_common))
rownames(aa_tf_gene_lfc_cor) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)
colnames(aa_tf_gene_lfc_cor) <- rownames(aa_pgeneEntrez_fc_common)

aa_tf_gene_fc_cor <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                            ncol = nrow(aa_pgeneEntrez_fc_common))
rownames(aa_tf_gene_fc_cor) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)
colnames(aa_tf_gene_fc_cor) <- rownames(aa_pgeneEntrez_fc_common)

for(i in 1:nrow(aa_tf_gene_lfc_cor)){
  for(j in 1:ncol(aa_tf_gene_lfc_cor)){
    aa_tf_gene_lfc_cor[i, j] <- cor(aa_tffc[i, ], 
                                    aa_pgeneEntrez_fc_common[j, ], 
                                    use = "pairwise.complete.obs")
    aa_tf_gene_fc_cor[i, j] <- cor(aa_tffc[i, ], 
                                   aafc[j, ], 
                                   use = "pairwise.complete.obs")
  }
}

heatmap.2(aa_tf_gene_lfc_cor[-c(2, 4),], Rowv = T, 
          Colv = T, na.rm = T, dendrogram = "both", trace = "none", 
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "TF_discrete_FC vs Gene_LFC")

heatmap.2(aa_tf_gene_fc_cor[-c(2, 4),], Rowv = T, 
          Colv = T, na.rm = T, dendrogram = "both", trace = "none", 
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "TF_discrete_FC vs Gene_FC")


# compare these to the currenlty used TF expression matrix (quantiled) against discretized exp mat
aa_tf_gene_lfc_cor <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                             ncol = nrow(aa_pgeneEntrez_fc_common))
rownames(aa_tf_gene_lfc_cor) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)
colnames(aa_tf_gene_lfc_cor) <- rownames(aa_pgeneEntrez_fc_common)

aa_tf_gene_fc_cor <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                            ncol = nrow(aa_pgeneEntrez_fc_common))
rownames(aa_tf_gene_fc_cor) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)
colnames(aa_tf_gene_fc_cor) <- rownames(aa_pgeneEntrez_fc_common)

for(i in 1:nrow(aa_tf_gene_lfc_cor)){
  for(j in 1:ncol(aa_tf_gene_lfc_cor)){
    aa_tf_gene_fc_cor[i, j] <- cor(aa_tffc[i, ], 
                                    my_CommonDifExpMat_16_ERassoc_gte5nzAll_172[j, ], 
                                    use = "pairwise.complete.obs")
  }
}

heatmap.2(aa_tf_gene_fc_cor[-c(2, 4),], Rowv = T, 
          Colv = T, na.rm = T, dendrogram = "both", trace = "none", 
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "TF_discrete_FC vs Gene_Discrete")
###############################
# TF discrete vs 2 ^ gene_discrete
aa_tf_gene_fc_cor <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                            ncol = nrow(aa_pgeneEntrez_fc_common))
rownames(aa_tf_gene_fc_cor) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)
colnames(aa_tf_gene_fc_cor) <- rownames(aa_pgeneEntrez_fc_common)

for(i in 1:nrow(aa_tf_gene_lfc_cor)){
  for(j in 1:ncol(aa_tf_gene_lfc_cor)){
    aa_tf_gene_fc_cor[i, j] <- cor(aa_tffc[i, ], 
                                   (2 ^ my_CommonDifExpMat_16_ERassoc_gte5nzAll_172)[j, ], 
                                   use = "pairwise.complete.obs")
  }
}

heatmap.2(aa_tf_gene_fc_cor[-c(2, 4),], Rowv = T, 
          Colv = T, na.rm = T, dendrogram = "both", trace = "none", 
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "TF_discrete_FC vs Gene_FC_Discrete")



###############################
# Look at the correlation between raw and discrete TF fold change
aa_dis_raw_cor <- matrix(nrow = nrow(aa_tffc),
                         ncol = nrow(aa_tffc))
rownames(aa_dis_raw_cor) <- rownames(aa_tffc)
colnames(aa_dis_raw_cor) <- rownames(aa_tffc)
for(i in 1:nrow(aa_dis_raw_cor)){
  for(j in 1:ncol(aa_dis_raw_cor)){
    aa_dis_raw_cor[i, j] <- cor(aa_tffc[i, ], 
                                TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_FC_trimmed_common[j, ], 
                                use = "pairwise.complete.obs")
  }
}
heatmap.2(aa_dis_raw_cor, Rowv = F, 
          Colv = F, na.rm = T, dendrogram = "none", trace = "none", 
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "TF_discrete_FC vs Raw_TF_FC")

aa_dis_raw_cor <- matrix(nrow = nrow(aa_tffc),
                         ncol = nrow(aa_tffc))
rownames(aa_dis_raw_cor) <- rownames(aa_tffc)
colnames(aa_dis_raw_cor) <- rownames(aa_tffc)
for(i in 1:nrow(aa_dis_raw_cor)){
  for(j in 1:ncol(aa_dis_raw_cor)){
    aa_dis_raw_cor[i, j] <- cor(log2(aa_tffc[i, ]), 
                                TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common[j, ], 
                                use = "pairwise.complete.obs")
  }
}
heatmap.2(aa_dis_raw_cor, Rowv = F, 
          Colv = F, na.rm = T, dendrogram = "none", trace = "none", 
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "TF_discrete_LFC vs Raw_TF_LFC")

###############################
# Look at partial correlation of each TF LFC with gene LFC, controlling for 
library(psych)
aa_tf_gene_lfc_pcor <- matrix(nrow = nrow(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                             ncol = nrow(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_rawLFC))
rownames(aa_tf_gene_lfc_pcor) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)
colnames(aa_tf_gene_lfc_pcor) <- rownames(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_rawLFC)
aa_mypartial <- purrr::possibly(.f = partial.r, otherwise = NA)
for(i in 1:nrow(aa_tf_gene_lfc_pcor)){
  for(j in 1:ncol(aa_tf_gene_lfc_pcor)){
    aa_tmp <- cbind(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_rawLFC[j, ], 
                    t(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common))
    aa_tmp <- aa_tmp[!is.na(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_rawLFC[j, ]), ]
    colnames(aa_tmp) <- c(rownames(aa_pgeneEntrez_fc_common)[j], 
                          rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common))
    aatmpres <- aa_mypartial(data = aa_tmp,
                             x = c(rownames(aa_pgeneEntrez_fc_common)[j], 
                                   rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)[i]),
                             y = setdiff(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                                         rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)[i]))
    if(! all(is.na(aatmpres))){
      aa_tf_gene_lfc_pcor[i, j] <- aatmpres[1, 2]
    }
    
  }
}
cor.plot(r = aa_tf_gene_lfc_pcor)
aa_tf_gene_lfc_pcor <- aa_tf_gene_lfc_pcor[, colSums(is.na(aa_tf_gene_lfc_pcor)) == 0]
heatmap.2(aa_tf_gene_lfc_pcor, Rowv = T, 
          Colv = T, na.rm = T, dendrogram = "both", trace = "none", 
          breaks = seq(-1, 1, length.out = 21),
          col = colorspace::diverge_hcl(20), main = "TF_raw_LFC vs Gene_raw_LFC partial")

my_CommonDifExpMat_16_ERassoc_gte5nzAll_172
aa <- partial.r(data = aa_tmp,
                x = c(rownames(aa_pgeneEntrez_fc_common)[j], 
                      rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)[i]),
                y = setdiff(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common),
                            rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_LFC_trimmed_common)[i])
                ,use = "pairwise.complete.obs")
############################################################################################################################
######################################   TF expression matrix based on raw fold change   #######################################################
############################################################################################################################
# In order to create the new expression matrix, I use the descretized version of expression (using quantiles) for
#  control experiment and that times the raw fold change to get the expression in treatment conditions.
# Modify the TFexpressionDescritizer function to get raw fold changes as input and do what is mentioned above
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_BORAW <- list()

for(i in 1:(length(my_Descretized_List)-1)){
  aaa <- unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2[[i]])
  stopifnot(length(aaa) == length(TF.motifs.Shrinked.hocomoco_E2F1_added))
  names(aaa) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added)
  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_BORAW[[i]] <- TFexpressionDescritizer(createDescreteMat_output = my_Descretized_List[[i]],
                                                                                              nu_quant = 5,
                                                                                              use_raw_lfc = T,
                                                                                              TF_index = aaa, 
                                                                                              return_raw_lfc = T)
}
aaContWT <- my_Descretized_List[[16]]$FullexpMat[, 1]
aaContWT_qu <- quantile(aaContWT, probs = seq(0, 1, length.out = 6))
aacontKD1 <- my_Descretized_List[[16]]$FullexpMat[, 3]
aacontKD2 <- my_Descretized_List[[16]]$FullexpMat[, 5]
aacontKD1_qu <- quantile(aacontKD1, probs = seq(0, 1, length.out = 6))
aacontKD2_qu <- quantile(aacontKD2, probs = seq(0, 1, length.out = 6))


aaContWT_dis <- QuantileCompare(exp_vec = aaContWT, quantile_vec = aaContWT_qu, 
                                assignedVal = seq(0.2, 1, length.out = 5))
aatreatWT_dis <- aaContWT_dis * (2 ^ my_Descretized_List[[16]]$raw_lfc[, 1])
aaContKD1_dis <- QuantileCompare(exp_vec = aacontKD1, quantile_vec = aacontKD1_qu, 
                                 assignedVal = seq(0.2, 1, length.out = 5))
aaContKD2_dis <- QuantileCompare(exp_vec = aacontKD2, quantile_vec = aacontKD2_qu,
                                 assignedVal = seq(0.2, 1, length.out = 5))
aaContKD_dis <- rowMeans(cbind(aaContKD1_dis, aaContKD2_dis))
aatreatKD_dis <- aaContKD_dis * (2 ^ my_Descretized_List[[16]]$raw_lfc[, 2])

aaa <- unlist(TF.Index.Shrinked.microarray_ENTREZ_PerDataset_unique_2[[16]])
aaAll_dis <- cbind(aaContWT_dis, aatreatWT_dis, aaContKD_dis, aatreatKD_dis)
colnames(aaAll_dis) <- c("WT_control", "WT_E2_12h", "siAP2g_control", "siAP2g_E2_12h")
aaTFdis <- aaAll_dis[aaa, ]
rownames(aaTFdis) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added)
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_BORAW[[16]] <- list(discrete_TF_mat = aaTFdis,
                                                                          Raw_LFC=my_Descretized_List[[16]]$raw_lfc[aaa,])


names(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_BORAW) <- names(my_Descretized_List)
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW <- do.call(cbind, lapply(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_2_BORAW, "[[", 1))
# Set expression of ER in treatment condition equal to its expression in control condition * 10 
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW
for(i in 1:45){
  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times[4, 2 * i] <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW[4, ((2 * i) - 1)] * 10
}

heatmap.2(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times, 
          Rowv = T, 
          Colv = F, na.rm = T, dendrogram = "row", trace = "none", 
          breaks = seq(0, 3, length.out = 21), margins = c(10,10),
          col = colorspace::terrain_hcl(20), main = "Discrete_TF_expression_ER10x")
heatmap.2(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW, 
          Rowv = T, 
          Colv = F, na.rm = T, dendrogram = "row", trace = "none", 
          breaks = seq(0, 3, length.out = 21), margins = c(10,10),
          col = colorspace::terrain_hcl(20), main = "Discrete_TF_expression")
############################################################################################################################
# RNA-seq TF expression after considering Raw fold changes
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_2_BORAW <- list()
for(i in 1:(length(my_Descretized_List_RNAseq))){
  aaa <- unlist(TF.Index.Shrinked.RNAseq_ENTREZ_PerDataset_2[[i]])
  stopifnot(length(aaa) == length(TF.motifs.Shrinked.hocomoco_E2F1_added))
  names(aaa) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added)
  TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_2_BORAW[[i]] <- TFexpressionDescritizer(createDescreteMat_output = my_Descretized_List_RNAseq[[i]],
                                                                                   TF_index = aaa,
                                                                                   nu_quant = 5,
                                                                                   return_raw_lfc = T, 
                                                                                   use_raw_lfc = T)
  }
names(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_2_BORAW) <- names(my_Descretized_List_RNAseq)

TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_BORAW <- do.call(cbind, lapply(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_2_BORAW, "[[", 1))
# Set expression of ER in treatment condition equal to its expression in control condition * 10 
TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times <- TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_BORAW
for(i in 1:14){
  TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times[4, 2 * i] <- TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times[4, ((2 * i) - 1)] * 10
}

heatmap.2(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times, 
          Rowv = T, 
          Colv = F, na.rm = T, dendrogram = "row", trace = "none", 
          breaks = seq(0, 3, length.out = 21), margins = c(10,10),
          col = colorspace::terrain_hcl(20), main = "Discrete_TF_expression_RNASeq_ER10x")
heatmap.2(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_BORAW, 
          Rowv = T, 
          Colv = F, na.rm = T, dendrogram = "row", trace = "none", 
          breaks = seq(0, 3, length.out = 21), margins = c(10,10),
          col = colorspace::terrain_hcl(20), main = "Discrete_TF_expression_RNASEQ")

##################################################################################################################
######################################   DONE WITH TF EXPRESSION   ##############################################
##################################################################################################################
# TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2 : is TF expression profile after adding E2F1 : no further modification is done on this one
# TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01 : is the above matrix after chaniging ER expression to 0 in control and 1 in treatment

# TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW : is TF expression matrix using the quantiles (0.2, 0.4, 0.6, 0.8 , 1) for control conds and the raw fold changes for treatment conds
# TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times : is the above matrix after chaniging ER expression to quantiles of control and times 10 for treatment.

# TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_BORAW : is TF expression matrix using the quantiles (0.2, 0.4, 0.6, 0.8 , 1) for control conds and the raw fold changes for treatment conds for RNA seq data
# TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times : is the above matrix after chaniging ER expression to quantiles of control and times 10 for treatment.
##################################################################################################################
##################################################################################################################
#####################################                                     ########################################
#####################################   Assign zero to unchanging genes   ########################################
#####################################                                     ########################################
##################################################################################################################
##################################################################################################################
# I assign zero to genes with log fold change ~ 0 and assign NA to the ones that are not important (not in top 100, and not 0 lfc)
# Done in the functions
##################################################################################################################
##################################################################################################################
##########################################                           #############################################
##########################################     First exp mat         #############################################
##########################################                           #############################################
##################################################################################################################
##################################################################################################################
#The matrix of gene expression with only genes that have at least 4 nonzero entries in KD conditions, AND have at least one +1 and one -1 
my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52

#The matrix of gene expression with only genes that have at least 5 nonzero and non-nA entries in all conditions
my_CommonDifExpMat_16_ERassoc_gte5nzAll_172

##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################

##################################################################################################################
##################################################################################################################
##################################                                   #############################################
##################################        Enhancer processing        #############################################
##################################                                   #############################################
##################################################################################################################
##################################################################################################################

# gather all regulatory elements of the genes of interest
aaa <- match(rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), names(Vicinity100kb.Enhancer.plusPromoter.By.gene))
ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52 <- WriteFastaOfBag(enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene[aaa],
                                                                             Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                                                             promoter.promoter.int=Promoter.gene.PromoterIntList[aaa],
                                                                             promoter.sequence=Promoter.Sequence.Char,
                                                                             return.Length = T,
                                                                             returnListNotwrite = T)
aa <- ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$sequence
for(i in 1:length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$sequence)){
  names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$sequence[[i]]) <- character(length = length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$sequence[[i]]))
  for (j in 1:length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$sequence[[i]])){
    names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$sequence[[i]])[j] <- as.character(j)
  }
}
for(i in 1:length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$length)){
  names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$length[[i]]) <- character(length = length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$length[[i]]))
  for (j in 1:length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$length[[i]])){
    names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$length[[i]])[j] <- as.character(j)
  }
}
ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted <- list()
ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence <- unlist(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$sequence)  
ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$length <- unlist(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$length) 

aa <- ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$GRange[[1]]
for(i in 2:length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$GRange)){
  aa <- c(aa, ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52$GRange[[i]])
}

ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$GRange <- aa
aaa <- ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence
aaaGR <- ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$GRange
aaaGR <- aaaGR[!duplicated(aaa), ]
aaa <- aaa[!duplicated(aaa)]
aaa <- lapply(aaa, toupper)

######### ######### ######### #########  Scan for motifs  ######### ######### ######### ######### ######### ######### 
# remove degenerate positions at each end, add pseudo count 0.001
TF.motifs.Shrinked.remdegen <-  lapply(X = TF.motifs.Shrinked, FUN = RemoveDegenRow)
TF.motifs.Shrinked.t <- lapply(TF.motifs.Shrinked.remdegen, t)



aa <- getwd()
setwd("plots/Motif_logo/TF_motifs_shrinked_19")
for(i in 1:length(TF.motifs.Shrinked.t)){
  setEPS()
  postscript(paste0(names(TF.motifs.Shrinked.t)[i], ".eps"), width = 9, height = 6)
  seqLogo::seqLogo(TF.motifs.Shrinked.t[[i]])
  dev.off()
}
setwd(aa)


ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore <- list()
# Following 4 lines are commented out to prevent accidental run. they have been used.
# for(i in 1:length(aaa)){
#   print(i)
#   ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore[[i]] <- MotifScore(seq = aaa[i], bg = c(0.25, 0.25, 0.25, 0.25), motifList = TF.motifs.Shrinked.t)
# }
names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore) <- names(aaa)
ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_uniq_seq <- aaa
ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_uniq_GR <- aaaGR

######### ######### ######### #########  For the 172 genes
# gather all regulatory elements of the genes of interest
aaa <- match(rownames(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172), names(Vicinity100kb.Enhancer.plusPromoter.By.gene))
ER.associated.reg.elements_microarray_gte5nzAll_172 <- WriteFastaOfBag(enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene[aaa],
                                                                       Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                                                       promoter.promoter.int=Promoter.gene.PromoterIntList[aaa],
                                                                       promoter.sequence=Promoter.Sequence.Char,
                                                                       return.Length = T,
                                                                       returnListNotwrite = T)
aa <- ER.associated.reg.elements_microarray_gte5nzAll_172$sequence
for(i in 1:length(ER.associated.reg.elements_microarray_gte5nzAll_172$sequence)){
  names(ER.associated.reg.elements_microarray_gte5nzAll_172$sequence[[i]]) <- character(length = length(ER.associated.reg.elements_microarray_gte5nzAll_172$sequence[[i]]))
  for (j in 1:length(ER.associated.reg.elements_microarray_gte5nzAll_172$sequence[[i]])){
    names(ER.associated.reg.elements_microarray_gte5nzAll_172$sequence[[i]])[j] <- as.character(j)
  }
}
for(i in 1:length(ER.associated.reg.elements_microarray_gte5nzAll_172$length)){
  names(ER.associated.reg.elements_microarray_gte5nzAll_172$length[[i]]) <- character(length = length(ER.associated.reg.elements_microarray_gte5nzAll_172$length[[i]]))
  for (j in 1:length(ER.associated.reg.elements_microarray_gte5nzAll_172$length[[i]])){
    names(ER.associated.reg.elements_microarray_gte5nzAll_172$length[[i]])[j] <- as.character(j)
  }
}
ER.associated.reg.elements_microarray_gte5nzAll_172_unlisted <- list()
ER.associated.reg.elements_microarray_gte5nzAll_172_unlisted$sequence <- unlist(ER.associated.reg.elements_microarray_gte5nzAll_172$sequence)  
ER.associated.reg.elements_microarray_gte5nzAll_172_unlisted$length <- unlist(ER.associated.reg.elements_microarray_gte5nzAll_172$length) 

aa <- ER.associated.reg.elements_microarray_gte5nzAll_172$GRange[[1]]
for(i in 2:length(ER.associated.reg.elements_microarray_gte5nzAll_172$GRange)){
  aa <- c(aa, ER.associated.reg.elements_microarray_gte5nzAll_172$GRange[[i]])
}

ER.associated.reg.elements_microarray_gte5nzAll_172_unlisted$GRange <- aa
aaa <- ER.associated.reg.elements_microarray_gte5nzAll_172_unlisted$sequence
aaaGR <- ER.associated.reg.elements_microarray_gte5nzAll_172_unlisted$GRange
aaaGR <- aaaGR[!duplicated(aaa), ]
aaa <- aaa[!duplicated(aaa)]
aaa <- lapply(aaa, toupper)
# after adding E2F1

TF.motifs.Shrinked.hocomoco_E2F1_added.linear.remdegen <-  lapply(X = TF.motifs.Shrinked.hocomoco_E2F1_added_linear, FUN = RemoveDegenRow)
TF.motifs.Shrinked.hocomoco_E2F1_added.linear.t <- lapply(TF.motifs.Shrinked.hocomoco_E2F1_added.linear.remdegen, t)


ER.associated.reg.elements_microarray_gte5nzAll_172_uniq_MotifScore <- list()


#Following 4 lines are commented out to prevent accidental run. they have been used.
for(i in 1:length(aaa)){
  print(i)
  ER.associated.reg.elements_microarray_gte5nzAll_172_uniq_MotifScore[[i]] <- MotifScore(seq = aaa[i], 
                                                                                         bg = c(0.25, 0.25, 0.25, 0.25), 
                                                                                         motifList = TF.motifs.Shrinked.hocomoco_E2F1_added.linear.t)
}
names(ER.associated.reg.elements_microarray_gte5nzAll_172_uniq_MotifScore) <- names(aaa)
ER.associated.reg.elements_microarray_gte5nzAll_172_unlisted_uniq_seq <- aaa
ER.associated.reg.elements_microarray_gte5nzAll_172_unlisted_uniq_GR <- aaaGR


######### ######### ######### #########  Chop and score motifs per peice  ######### ######### ######### ######### ######### 
# Chop enhancers to 1kb 250 step size and score them for each TF
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece <- EnhancerChopper(Enhancer_seq = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_uniq_seq,
                                                                                                  motifScore_output = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore,
                                                                                                  motifList = TF.motifs.Shrinked.t,
                                                                                                  no_thresh = T,
                                                                                                  LLR_thresh = 0,
                                                                                                  piece_length = 1000,
                                                                                                  step_size = 250 )
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
# map the scores back to the original enhancer per gene system
aaa <- match(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence,
             ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_uniq_seq)
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece <- list()
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Seq <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece$Chopped_Seq[aaa]
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Score <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece$Chopped_Score[aaa]
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Seq) <- names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence)
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Score) <- names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence)
for(i in 1:length(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Seq)){
  names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Seq[[i]]) <- paste(names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Seq)[i], c(1:length(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Seq[[i]])), sep = "_")
  rownames(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Score[[i]]) <- paste(names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Seq)[i], c(1:length(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Seq[[i]])), sep = "_")
}

# for each gene show enhancers binding affinity for different TFs in a heatmap

aabind <- do.call(rbind, ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece$Chopped_Score)
aa <- unlist(lapply(strsplit(rownames(aabind), split="\\."), "[[", 1))
aarowcol <- integer(length(aa))
for(i in 1:length(unique(aa))){
  aarowcol[aa %in% unique(aa)[i]] <- i
}
length(unique(aarowcol))
aarowcol2 <- col_vector[aarowcol]


cor(apply(aabind, 2, median), unlist(lapply(TF.motifs.Shrinked.t, ncol)))
infContentCalc()
aa <- lapply(X = TF.motifs.Shrinked.t, infContentCalc)
aa <- unlist(lapply(X = aa, sum))
cor(aa, unlist(lapply(TF.motifs.Shrinked.t, ncol)))

aabind2 <-aabind * aa
boxplot.matrix(aabind, las= 2, outline = F)
boxplot.matrix(log(aabind), las = 2)

boxplot(aabind[, 3], outline = F)
sum(log(aabind[, 3]) > -15)
aach <- aarowcol[log(aabind[, 3]) > -12]
length(unique(aach))

sort(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`9687.4`$Forward$ESR1[, 4], decreasing = T)


length(unique(aarowcol))
png("Enhancer_affinity.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)   

#plot the heatmap using gplot heatmap.2

heatmap.2(log(aabind[,c(1, 3)]),
          #cellnote = format(round(SSEmat, 2), nsmall = 2),  # same data set for cell labels
          main = "affinity matrix", # heat map title
          notecol="black",      # change font color of cell labels to black
          # density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(9,7),     # widens margins around plot
         # col=my_palette,       # use on color palette defined earlier
        #  breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",     # only draw a .dendrogram dendrogram
          RowSideColors = aarowcol,# grouping row-variables(genes) into different categories
          #ColSideColors = .ColSideColors,# grouping column-variables(models) into different categories
          Rowv = F,
          Colv = F,
          scale = "none",
          na.rm=TRUE
) 

dev.off() 

names(TF.motifs.Shrinked.t)
seqLogo::seqLogo(TF.motifs.Shrinked.t[[10]])
seqLogo::seqLogo(TF.motifs.TimeSeries.loose.t[[2]])

################################################################################################################
################################################################################################################
# Filter for enhancers that pass my threshold of
# look at the score of my enhancers within the previously calculated scores of ER-ER dimer sites:

ER.associated.reg.elements_uniq_pair_Concat_sort_un

aa <- do.call(c, ER.associated.reg.elements$sequence)
aaa <- aa[!duplicated(aa)]

aaa2 <- match(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence,
             aaa)
#get the list of motifhits of ER monomers at most 3 bp apart, for each sequence of interest: this is LLR threshold 3
ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th3 <- ER.associated.reg.elements_uniq_pair_Concat_sort_un$LLR_Thresh_3$UnsortedOutput$ESR1_ESR1[aaa2]
names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th3) <- names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence)

sum(unlist(lapply(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th3, length)) == 0)

#get the list of motifhits of ER monomers at most 3 bp apart, for each sequence of interest: this is LLR threshold 2
ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th2 <- ER.associated.reg.elements_uniq_pair_Concat_sort_un$LLR_Thresh_2$UnsortedOutput$ESR1_ESR1[aaa2]
names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th2) <- names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence)
which(unlist(lapply(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th2, length)) == 0)




# DECIDED TO KEEP ALL original enhancers, but can filter the sliding windows which do not contain any sites.
ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th2_uniq <- ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th2[!duplicated(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence)]

aaToRange <- function(x){
  if(nrow(x) > 0){
    aa <- matrix(nrow = nrow(x), ncol = 2)
    aa[, 1] <- x[, 1]
    aa[, 2] <- x[, 2] + 5
    return(aa)
  }else{
    return(numeric(0))
  }
}
ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th2_uniq_range <- lapply(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th2_uniq, aaToRange)

ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece_filtered <- EnhancerChopper(Enhancer_seq = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_uniq_seq,
                                                                                                           EnhancerGR = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_uniq_GR,
                                                                                                           motifScore_output = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore,
                                                                                                           motifList = TF.motifs.Shrinked.t,
                                                                                                           no_thresh = T,
                                                                                                           LLR_thresh = 0,
                                                                                                           piece_length = 1000,
                                                                                                           step_size = 250,
                                                                                                           must_contain = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th2_uniq_range)
aaa <- match(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence,
             ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_uniq_seq)
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered <- list()
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered$Chopped_Seq <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece_filtered$Chopped_Seq[aaa]
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered$Chopped_Score <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece_filtered$Chopped_Score[aaa]
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered$Chopped_GRanges <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece_filtered$Chopped_GRanges[aaa]
# create a per gene list of enhancer affinity matrices using the names
aa <- names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered$Chopped_Score)
aaa <- unlist(lapply(strsplit(aa, split = "\\."), "[[", 1))
aaaa <- unique(aaa)
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene <- list()
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq <- list()
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score <- list()
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges <- list()
for(i in 1:length(aaaa)){
  aaw <- which(aaa %in% aaaa[i])
  aa1 <- do.call(rbind, ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered$Chopped_Score[aaw])
  ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score[[i]] <- aa1
  aa2 <- do.call(c, ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered$Chopped_Seq[aaw])
  ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq[[i]] <- aa2
  aa3 <- lapply(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered$Chopped_GRanges[aaw], as.data.frame)
  aa3 <- do.call(rbind, aa3)
  aa3 <- makeGRangesFromDataFrame(aa3)
  ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges[[i]] <- aa3
}
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score) <- aaaa
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq) <- aaaa
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges) <- aaaa

# computing the same thing:ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene but using the wrapper function and not normalizing affinities by maxLLR
aa_example <- EnhancerChopper_wrapper(Entrez_IDs = rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52),
                                      .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                      .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                      .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                      .promoter.sequence=Promoter.Sequence.Char,
                                      #.Enhancer_seq=character(0),
                                      #.EnhancerGR=character(0),
                                      #.motifScore_output,
                                      .motifList=TF.motifs.Shrinked.t[3],
                                      .no_thresh = T,
                                      .LLR_thresh = 0,
                                      .piece_length = 1000,
                                      .step_size = 250,
                                      ..diff_max_LLR = F,
                                      dimer_filter_ER=T,
                                      dimer_list=ER.associated.reg.elements_uniq_pair_Concat_sort_un,
                                      dimer_filter_ER_LLR_th=2, 
                                      dimer_filter_corresponding_seq_nonProcessed=ER.associated.reg.elements$sequence,
                                      GRange_overlap=F,
                                      GRange_overlap_object=ReMapChIP.GRanges.list[52],
                                      GRange_overlap_percent=50)
######################
# for 172 genes after adding E2F1
my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1 <- EnhancerChopper_wrapper(Entrez_IDs = rownames(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172),
                                                                                            .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                                                                            .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                                                                            .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                                                                            .promoter.sequence=Promoter.Sequence.Char,
                                                                                            unique_motifScore = ER.associated.reg.elements_microarray_gte5nzAll_172_uniq_MotifScore,
                                                                                            #.Enhancer_seq=character(0),
                                                                                            #.EnhancerGR=character(0),
                                                                                            #.motifScore_output,
                                                                                            .motifList=TF.motifs.Shrinked.hocomoco_E2F1_added_linear,
                                                                                            .no_thresh = T,
                                                                                            .LLR_thresh = 0,
                                                                                            .piece_length = 1000,
                                                                                            .step_size = 250,
                                                                                            ..diff_max_LLR = F,
                                                                                            dimer_filter_ER=F,
                                                                                            dimer_list=ER.associated.reg.elements_uniq_pair_Concat_sort_un,
                                                                                            dimer_filter_ER_LLR_th=2, 
                                                                                            dimer_filter_corresponding_seq_nonProcessed=ER.associated.reg.elements$sequence,
                                                                                            GRange_overlap=T,
                                                                                            GRange_overlap_object=ReMapChIP.GRanges.list_merged[52],
                                                                                            GRange_overlap_percent=50)
names(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Seq$`10057`)
names(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Chopped_Seq$`10057`)

#####################################################################################################################
# binding all together
aabind2 <- do.call(rbind, ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered$Chopped_Score)
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
plotExpression(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, filename = "genes_gte4nonzeroKD_atleast1p1n_heatmap_2.png",.dendrogram = "both", .Colv = T, .distfun = my_dist_fuc)

################################################################################################################
################################################################################################################
# Filter the enhancers for the ones that overlap with an ER chip peak.
aa <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges
#aa <- lapply(aa, as.data.frame)
#aa2 <- do.call(rbind, aa)
aaint <- list()
for(i in 1:length(aa)){
  aaint[[i]] <-  OverlapInteractionBatchExtract(subjectCoord1=aa[[i]],
                                                subjectCoord2List=ReMapChIP.GRanges.list[52],
                                                mode = "both" )
}
aa2 <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq
aa3 <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score
aa1_filtered <- list()
aa2_filtered <- list()
aa3_filtered <- list()
for(i in 1:length(aa)){
  aaind <- which((aaint[[i]]$OverlapMat + aaint[[i]]$IntMat) != 0)
  print(length(aaind))
  aa1_filtered[[i]] <- aa[[i]][aaind]
  aa2_filtered[[i]] <- aa2[[i]][aaind]
  aa3_filtered[[i]] <- aa3[[i]][aaind, ]
}

ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_ChIPfilt <- list(Chopped_Seq=aa2_filtered,
                                                                                                            Chopped_Score=aa3_filtered,
                                                                                                            Chopped_GRanges=aa1_filtered)

################################################################################################################
################################################################################################################
#create an enhancer list like the one above but with lower resolution, meaning each enhancer 1500 bp and steps 500 bp

ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece_filtered_1500_500 <- EnhancerChopper(Enhancer_seq = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_uniq_seq,
                                                                                                                    EnhancerGR = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_uniq_GR,
                                                                                                                    motifScore_output = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore,
                                                                                                                    motifList = TF.motifs.Shrinked.t,
                                                                                                                    no_thresh = T,
                                                                                                                    LLR_thresh = 0,
                                                                                                                    piece_length = 1500,
                                                                                                                    step_size = 500,
                                                                                                                    must_contain = ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th2_uniq_range)
aaa <- match(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted$sequence,
             ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_uniq_seq)
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500 <- list()
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500$Chopped_Seq <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece_filtered_1500_500$Chopped_Seq[aaa]
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500$Chopped_Score <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece_filtered_1500_500$Chopped_Score[aaa]
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500$Chopped_GRanges <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece_filtered_1500_500$Chopped_GRanges[aaa]
# create a per gene list of enhancer affinity matrices using the names
aa <- names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500$Chopped_Score)
aaa <- unlist(lapply(strsplit(aa, split = "\\."), "[[", 1))
aaaa <- unique(aaa)
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene <- list()
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Seq <- list()
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Score <- list()
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_GRanges <- list()
for(i in 1:length(aaaa)){
  aaw <- which(aaa %in% aaaa[i])
  aa1 <- do.call(rbind, ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500$Chopped_Score[aaw])
  ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Score[[i]] <- aa1
  aa2 <- do.call(c, ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500$Chopped_Seq[aaw])
  ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Seq[[i]] <- aa2
  aa3 <- lapply(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500$Chopped_GRanges[aaw], as.data.frame)
  aa3 <- do.call(rbind, aa3)
  aa3 <- makeGRangesFromDataFrame(aa3)
  ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_GRanges[[i]] <- aa3
}
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Score) <- aaaa
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Seq) <- aaaa
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_GRanges) <- aaaa

############################################################################################################################################
# Filter the enhancers for the ones that overlap with an ER chip peak.
aa <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_GRanges
#aa <- lapply(aa, as.data.frame)
#aa2 <- do.call(rbind, aa)
aaint <- list()
for(i in 1:length(aa)){
  aaint[[i]] <-  OverlapInteractionBatchExtract(subjectCoord1=aa[[i]],
                                                subjectCoord2List=ReMapChIP.GRanges.list[52],
                                                mode = "both"
                                                #,...minOverlapPercent = 50
                                                )
}
aa2 <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Seq
aa3 <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Score
aa1_filtered <- list()
aa2_filtered <- list()
aa3_filtered <- list()
for(i in 1:length(aa)){
  aaind <- which((aaint[[i]]$OverlapMat + aaint[[i]]$IntMat) != 0)
  print(length(aaind))
  aa1_filtered[[i]] <- aa[[i]][aaind]
  aa2_filtered[[i]] <- aa2[[i]][aaind]
  if(length(aaind) == 1){
    print(paste("1 out of ", nrow(aa3[[i]]), "remained"))
    aa3_filtered[[i]] <- t(as.matrix(aa3[[i]][aaind, ]))
  }else if(length(aaind) == 0){
    print(paste("0 left out of ", nrow(aa3[[i]])))
    aa3_filtered[[i]] <- t(as.matrix(aa3[[i]][aaind, ]))
  }else{
    print(paste(length(aaind), " out of ", nrow(aa3[[i]]), "remained"))
    aa3_filtered[[i]] <- aa3[[i]][aaind, ]
  }
}

ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt <- list(Chopped_Seq=aa2_filtered,
                                                                                                                     Chopped_Score=aa3_filtered,
                                                                                                                     Chopped_GRanges=aa1_filtered)
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt$Chopped_Score) <- names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Score)
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt$Chopped_Seq) <-names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Seq)
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt$Chopped_GRanges) <-names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_GRanges)
#######################################################################################################
# Filter enhancers for overlap with ChiP peaks, criteria: at least 50% of the length of the enhancer should be within the chip peak

aa <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_GRanges
#aa <- lapply(aa, as.data.frame)
#aa2 <- do.call(rbind, aa)
aaint <- list()
for(i in 1:length(aa)){
  aaint[[i]] <-  OverlapInteractionBatchExtract(subjectCoord1=aa[[i]],
                                                subjectCoord2List=ReMapChIP.GRanges.list[52],
                                                mode = "both"
                                                ,...minOverlapPercent = 50,
                                                .int_piece_max_length = 10000)
}
aa2 <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Seq
aa3 <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Score
aa1_filtered <- list()
aa2_filtered <- list()
aa3_filtered <- list()
for(i in 1:length(aa)){
  # aaind <- which((aaint[[i]]$OverlapMat + aaint[[i]]$IntMat) > 1)
  aaind <- which((aaint[[i]]$OverlapMat) > 0)
  print(length(aaind))
  aa1_filtered[[i]] <- aa[[i]][aaind]
  aa2_filtered[[i]] <- aa2[[i]][aaind]
  if(length(aaind) == 1){
    print(paste("1 out of ", nrow(aa3[[i]]), "remained"))
    aa3_filtered[[i]] <- t(as.matrix(aa3[[i]][aaind, ]))
  }else if(length(aaind) == 0){
    print(paste("0 left out of ", nrow(aa3[[i]])))
    aa3_filtered[[i]] <- t(as.matrix(aa3[[i]][1, ]))
  }else{
    print(paste(length(aaind), " out of ", nrow(aa3[[i]]), "remained"))
    aa3_filtered[[i]] <- aa3[[i]][aaind, ]
  }
}

ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt_50p <- list(Chopped_Seq=aa2_filtered,
                                                                                                                     Chopped_Score=aa3_filtered,
                                                                                                                     Chopped_GRanges=aa1_filtered)
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt_50p$Chopped_Score) <- names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Score)
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt_50p$Chopped_Seq) <-names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Seq)
names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt_50p$Chopped_GRanges) <-names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_GRanges)

length(unlist(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt_50p$Chopped_Seq))

table(nchar(unlist(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500$Chopped_Seq)))

aa_bpp <- (lapply(aaint, "[[", 4))
max(unlist(lapply(aa_bpp, nrow)))
aa_bp <- matrix(nrow = max(unlist(lapply(aa_bpp, nrow))), ncol = length(aaint))
for(i in 1:length(aaint)){
  aa_bp[1:nrow(aaint[[i]]$IntMat), i] <- aaint[[i]]$IntMat
}
boxplot.matrix(aa_bp, las =2, main = "nu.of.inteacting chip peaks")

aa_bpo <- (lapply(aaint, "[[", 2))
max(unlist(lapply(aa_bpo, nrow)))
aa_bo <- matrix(nrow = max(unlist(lapply(aa_bpo, nrow))), ncol = length(aaint))
for(i in 1:length(aaint)){
  aa_bo[1:nrow(aaint[[i]]$OverlapMat), i] <- aaint[[i]]$OverlapMat
}
boxplot.matrix(aa_bo, las =2, main = "nu.of.overlapping chip peaks")

aa_er_ch_w <- width(ReMapChIP.GRanges.list[[52]])






