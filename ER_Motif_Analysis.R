# ER project, motif analysis
#in this script I write functions to compute and visualize motif distributions for given enhancers and genes

###################################################################################################################
##############################################                    #################################################
##############################################      Libraries     #################################################
##############################################                    #################################################
###################################################################################################################


###################################################################################################################
##############################################                    #################################################
##############################################      Functions     #################################################
##############################################                    #################################################
###################################################################################################################
# Functions used from other scripts:
# MotifScore and CompLLRRandom used from SyntheticEnhancers.R 
# Addpsudo from DNAshape.R

# Functions written in this script:
# 1) MotifVisualizer
# 2) findNextOccuring
# 3) motifPairInvestigator
# 4) EnrichedPairs
# 5) 
###################################################################################################################
###################################################################################################################
MotifHitVisualizer <- function(MotifScore_OutputList = numeric(0),
                               .motifList = TF.motifs.TimeSeries.loose.t,
                               seqList,
                               ChosenMotifIndex = numeric(0),
                               ChosenMotifIndex.Vicinity = T,
                               VicinityThresh = numeric(0),
                               LLRthresh = rep(0, length(.motifList)),
                               col_vec = col_vector,
                               plot = T,
                               pairInvestigateMat = matrix(nrow = 0, ncol = 4)){
  # this function gets a list of sequences, and a list of motifs and visulaizes the motif hits in the sequence as barplots
  # MotifScore_OutputList is a list of outputs of MotifScore function (in case it is precomputed)
  # motifList is a named list where each entry is a PWM
  # seqList is a list of sequences
  # ChosenMotifIndex if only a specific motif has to be visualized.
  # ChosenMotifIndex.Vicinity if the motif hits in the ("VicinityThresh (bp)") distance of the "ChosenMotifIndex" have to be visulized
  # VicinityThresh is an integer vector where each entry is the threshold distance for that chosen tf (each entry corresponds to the same index in ChosenMotifIndex) 
  # LLRthresh : a numeric vector: each entry corresponding to one motif in the .motifList. For each: visualize hits with LLR score higher than their corresponding LLRthresh
  # pairInvestigateMat is a matrix, each row indicates a pair of motifs to investigate in relation to eac other. first column: index of first motif (in .motifList), second column: index of second motif, third column: 1 if the have to be on the same strand, -1 if they have to be on opposite strands, forth column: max distance between the motifs
  # it outputs a list containing the motifScore outputs
  
  # compute motifScores if they are not already given as input
  if (length(MotifScore_OutputList) > 0){
    MotifScore_OutList <- MotifScore_OutputList
  }else{
    seqList <- lapply(X=seqList, FUN=toupper )
    MotifScore_OutList <- lapply(X=seqList,
                                 FUN=MotifScore,
                                 bg = c(0.25, 0.25, 0.25, 0.25),
                                 motifList = .motifList)
  }
  prev.Par <- par() #save the previous par setting
  #start plotting
  MotifScore_OutList_Plot <- MotifScore_OutList
  # set weak signals to zero
  for (cur.seq in 1:length(seqList)){
    for (tf in 1:length(.motifList)){
      MotifScore_OutList_Plot[[cur.seq]][[1]][[tf]][MotifScore_OutList_Plot[[cur.seq]][[1]][[tf]] < LLRthresh[tf]] <- 0
      MotifScore_OutList_Plot[[cur.seq]][[2]][[tf]][MotifScore_OutList_Plot[[cur.seq]][[2]][[tf]] < LLRthresh[tf]] <- 0
    }
  }
  
  if(length(ChosenMotifIndex) == 0){#if you want all motifs to be visualized
    MotifScore_OutList_Plot.Mod <- MotifScore_OutList_Plot

  }else{#if you want only certain motifs and maybe the ones within their VicinityThresh proximity to be visualized.
    MotifScore_OutList_Plot.Mod <- MotifScore_OutList_Plot
    if(ChosenMotifIndex.Vicinity){#if you want the ones in vicinity
      unchosen.tfs <- setdiff(c(1:length(.motifList)), ChosenMotifIndex)
      for(cur.seq in 1:length(seqList)){
        #create a rangemat that you want to visualize:
        seq.hits.rangemat <- matrix(nrow = 0, ncol = 2)
        for(chosen.tf in 1:length(ChosenMotifIndex)){
          cur.seq.hits.forw <- which(MotifScore_OutList_Plot.Mod[[cur.seq]][[1]][[ChosenMotifIndex[chosen.tf]]][, 1] > 0)
          cur.seq.hits.rev  <- which(MotifScore_OutList_Plot.Mod[[cur.seq]][[2]][[ChosenMotifIndex[chosen.tf]]][, 1] > 0)
          cur.seq.hits <- union(cur.seq.hits.forw, cur.seq.hits.rev)
          cur.seq.hits.rangemat <- matrix(nrow = length(cur.seq.hits), ncol = 2)
          cur.seq.hits.rangemat[, 1] <- cur.seq.hits - (VicinityThresh[chosen.tf] + 1)
          cur.seq.hits.rangemat[, 2] <- cur.seq.hits + VicinityThresh[chosen.tf] + ncol(.motifList[[ChosenMotifIndex[chosen.tf]]])
          seq.hits.rangemat <- rbind(seq.hits.rangemat, cur.seq.hits.rangemat)
        }#got the range of indices that I care about: now go over other TFs and considering their motif length modify this rangemat
        for(curr.unchosen.tf in 1:length(unchosen.tfs)){#get the index of the positions that we care about in other tfs
          cur.un.seq.hits.rangemat <- seq.hits.rangemat
          cur.un.seq.hits.rangemat[, 1] <- cur.un.seq.hits.rangemat[, 1] - (ncol(.motifList[[unchosen.tfs[curr.unchosen.tf]]])-1)
          cur.un.seq.hits.ind <- integer(0)
          for(i in 1:nrow(cur.un.seq.hits.rangemat)){
            cur.un.seq.hits.ind <- c(cur.un.seq.hits.ind, c(cur.un.seq.hits.rangemat[i, 1] : cur.un.seq.hits.rangemat[i, 2]))
          }
          cur.un.seq.hits.ind <- sort(cur.un.seq.hits.ind[!duplicated(cur.un.seq.hits.ind)], decreasing = F)#sort and remove duplicates
          cur.un.seq.hits.ind <- cur.un.seq.hits.ind[cur.un.seq.hits.ind < nrow(MotifScore_OutList_Plot.Mod[[cur.seq]][[1]][[unchosen.tfs[curr.unchosen.tf]]])] #remove the indices greater than length of motifscore result
          cur.un.seq.hits.ind.notCare <- setdiff(c(1 : nrow(MotifScore_OutList_Plot.Mod[[cur.seq]][[1]][[unchosen.tfs[curr.unchosen.tf]]])), cur.un.seq.hits.ind)
          MotifScore_OutList_Plot.Mod[[cur.seq]][[1]][[unchosen.tfs[curr.unchosen.tf]]][cur.un.seq.hits.ind.notCare, ] <- 0
          MotifScore_OutList_Plot.Mod[[cur.seq]][[2]][[unchosen.tfs[curr.unchosen.tf]]][cur.un.seq.hits.ind.notCare, ] <- 0
        }# end of loop over unchosen.tfs
      }#end of loop over sequences
    }else{#if you don't want the ones in vicinity
      for(cur.seq in 1:length(seqList)){
        MotifScore_OutList_Plot.Mod[[cur.seq]][[1]] <- MotifScore_OutList_Plot.Mod[[cur.seq]][[1]][ChosenMotifIndex]
        MotifScore_OutList_Plot.Mod[[cur.seq]][[2]] <- MotifScore_OutList_Plot.Mod[[cur.seq]][[2]][ChosenMotifIndex]
      }
    }
  }
  #now I have the modified motifScore output: MotifScore_OutList_Plot.Mod and I can plot it
  par(mfrow = c(length(seqList), 1), mar = c(0.1, 4, 0.1, 0.1))
  #find maximum LLR to set the ylim
  ylim0 <- 0
  for(i in 1:length(MotifScore_OutList_Plot.Mod[[1]][[1]])){
    if(MotifScore_OutList_Plot.Mod[[1]][[1]][[i]][1, 2] > ylim0){
      ylim0 <- MotifScore_OutList_Plot.Mod[[1]][[1]][[i]][1, 2]
    }
  }
  
  for(cur.seq in 1:length(seqList)){
    barplot(MotifScore_OutList_Plot.Mod[[cur.seq]][[1]][[1]][, 1], xlab = "", border=NA, las=2, col = 2, ylim = c(-1*ylim0,ylim0))
    for(cur.tf in 1:length(MotifScore_OutList_Plot.Mod[[cur.seq]][[1]])){
      barplot(      MotifScore_OutList_Plot.Mod[[cur.seq]][[1]][[cur.tf]][, 1], xlab = "", border=NA, las=2, col = col_vec[cur.tf], ylim = c(-1*ylim0,ylim0), add = T)
      barplot( -1 * MotifScore_OutList_Plot.Mod[[cur.seq]][[2]][[cur.tf]][, 1], xlab = "", border=NA, las=2, col = col_vec[cur.tf], ylim = c(-1*ylim0,ylim0), add = T)
    }
  }
  #switch back to previous par setting
  par(mfrow = prev.Par$mfrow, mar = prev.Par$mar)
  return(MotifScore_OutList_Plot.Mod)
}
#what else do you want?
#the number of non-overlapping occurances of one motif in certain threshold of another motif of the same or different type, one the same of different strand


###################################################################################################################
###################################################################################################################
# Example
aaSeq <- MCFEnhancers[1872,"sequ"]
aa <- MotifHitVisualizer(.motifList = TF.motifs.TimeSeries.loose.t, seqList = as.list(aaSeq), LLRthresh = c(3,3), ChosenMotifIndex = 1, ChosenMotifIndex.Vicinity = F, VicinityThresh = 5)

###################################################################################################################
###################################################################################################################
#this function will be used in the next function to find pairs
findNextOccuring <- function(MotifScoreModifiedList,
                             motif1Index, motif2Index,
                             motif1Length, motif2Length,
                             orient,
                             .seqList,
                             distTresh){
  # this function gets the modified output of motifScore, and the two motifs, finds
  # pairs with motif1 first, motif2 second in same or opposite strands
  # MotifScoreModified modified output of motifscore
  # motif1Index, motif2Index : index of motifs in their MotifScoreModifiedList
  # orient : if 1 looks for occurances on the same strand, if -1 looks for occurances on other strand, if 0 looks for both
  # .seqList list of sequences
  # distTresh is the distance between the end of one motif with the start of the next
  OccurList <- list()
  if (! orient %in% c(-1, 0, 1)){
    print("Orient is wrong. enter -1, 0, 1")
    return(-1)
  }
  for(cur.seq in 1:length(.seqList)){
    print(paste("finding pairs in sequence number", cur.seq))
    cur.1st.seq.hits.forw <- which(MotifScoreModifiedList[[cur.seq]][[1]][[motif1Index]][, 1] > 0)
    cur.1st.seq.hits.rev  <- which(MotifScoreModifiedList[[cur.seq]][[2]][[motif1Index]][, 1] > 0)
    cur.2nd.seq.hits.forw <- which(MotifScoreModifiedList[[cur.seq]][[1]][[motif2Index]][, 1] > 0)
    cur.2nd.seq.hits.rev  <- which(MotifScoreModifiedList[[cur.seq]][[2]][[motif2Index]][, 1] > 0)
    
    cur.seq.hits.rangemat.forw <- matrix(nrow = length(cur.1st.seq.hits.forw), ncol = 2)
    cur.seq.hits.rangemat.rev <-  matrix(nrow = length(cur.1st.seq.hits.rev),  ncol = 2)
    cur.seq.hits.rangemat.forw[, 1] <- cur.1st.seq.hits.forw + motif1Length
    cur.seq.hits.rangemat.forw[, 2] <- cur.1st.seq.hits.forw + motif1Length + distTresh
    rownames(cur.seq.hits.rangemat.forw) <- cur.1st.seq.hits.forw
    cur.seq.hits.rangemat.rev[, 1] <- cur.1st.seq.hits.rev + motif1Length
    # if((cur.1st.seq.hits.rev + motif1Length + distTresh) > )
    cur.seq.hits.rangemat.rev[, 2] <- cur.1st.seq.hits.rev + motif1Length + distTresh
    rownames(cur.seq.hits.rangemat.rev) <- cur.1st.seq.hits.rev
    if (orient != -1){
      print(paste("computing", "forwforw pairs"))
      forwforw <- matrix(nrow = 0, ncol = 4)
      if(nrow(cur.seq.hits.rangemat.forw) > 0){
        for(crow in 1:nrow(cur.seq.hits.rangemat.forw)){
          curoc <-  which(cur.2nd.seq.hits.forw %in% c(cur.seq.hits.rangemat.forw[crow, 1] : cur.seq.hits.rangemat.forw[crow, 2]))
          if(length(curoc) > 0){
            curocmat <- matrix(nrow = length(curoc), ncol = 4)
            colnames(curocmat) <- c("position1", "position2", "LLR1", "LLR2")
            for(curoccnnt in 1:length(curoc)){
              curocmat[curoccnnt, 1] <- cur.1st.seq.hits.forw[crow]
              curocmat[curoccnnt, 2] <- cur.2nd.seq.hits.forw[curoc[curoccnnt]]
              curocmat[curoccnnt, 3] <- MotifScoreModifiedList[[cur.seq]][[1]][[motif1Index]][curocmat[curoccnnt, 1], 1]
              curocmat[curoccnnt, 4] <- MotifScoreModifiedList[[cur.seq]][[1]][[motif2Index]][curocmat[curoccnnt, 2], 1]
            }
            forwforw <- rbind(forwforw, curocmat)
          }
        }
      }

      print(paste("computing", "revrev pairs"))
      revrev <- matrix(nrow = 0, ncol = 4)
      if (nrow(cur.seq.hits.rangemat.rev) > 0){
        for(crow in 1:nrow(cur.seq.hits.rangemat.rev)){
          curoc <-  which(cur.2nd.seq.hits.rev %in% c(cur.seq.hits.rangemat.rev[crow, 1] : cur.seq.hits.rangemat.rev[crow, 2]))
          if(length(curoc) > 0){
            curocmat <- matrix(nrow = length(curoc), ncol = 4)
            colnames(curocmat) <- c("position1", "position2", "LLR1", "LLR2")
            for(curoccnnt in 1:length(curoc)){
              curocmat[curoccnnt, 1] <- cur.1st.seq.hits.rev[crow]
              curocmat[curoccnnt, 2] <- cur.2nd.seq.hits.rev[curoc[curoccnnt]]
              curocmat[curoccnnt, 3] <- MotifScoreModifiedList[[cur.seq]][[2]][[motif1Index]][curocmat[curoccnnt, 1], 1]
              curocmat[curoccnnt, 4] <- MotifScoreModifiedList[[cur.seq]][[2]][[motif2Index]][curocmat[curoccnnt, 2], 1]
            }
            revrev <- rbind(revrev, curocmat)
          }
        }
      }
    }
    if(orient != 1){
      print(paste("computing", "forwrev pairs"))
      forwrev <- matrix(nrow = 0, ncol = 4)
      if(nrow(cur.seq.hits.rangemat.forw) > 0){
        for(crow in 1:nrow(cur.seq.hits.rangemat.forw)){
          curoc <-  which(cur.2nd.seq.hits.rev %in% c(cur.seq.hits.rangemat.forw[crow, 1] : cur.seq.hits.rangemat.forw[crow, 2]))
          if(length(curoc) > 0){
            curocmat <- matrix(nrow = length(curoc), ncol = 4)
            colnames(curocmat) <- c("position1", "position2", "LLR1", "LLR2")
            for(curoccnnt in 1:length(curoc)){
              curocmat[curoccnnt, 1] <- cur.1st.seq.hits.forw[crow]
              curocmat[curoccnnt, 2] <- cur.2nd.seq.hits.rev[curoc[curoccnnt]]
              curocmat[curoccnnt, 3] <- MotifScoreModifiedList[[cur.seq]][[1]][[motif1Index]][curocmat[curoccnnt, 1], 1]
              curocmat[curoccnnt, 4] <- MotifScoreModifiedList[[cur.seq]][[2]][[motif2Index]][curocmat[curoccnnt, 2], 1]
            }
            forwrev <- rbind(forwrev, curocmat)
          }
        }
      }

      print(paste("computing", "revforw pairs"))
      revforw <- matrix(nrow = 0, ncol = 4)
      if(nrow(cur.seq.hits.rangemat.rev) > 0){
        for(crow in 1:nrow(cur.seq.hits.rangemat.rev)){
          curoc <-  which(cur.2nd.seq.hits.forw %in% c(cur.seq.hits.rangemat.rev[crow, 1] : cur.seq.hits.rangemat.rev[crow, 2]))
          if(length(curoc) > 0){
            curocmat <- matrix(nrow = length(curoc), ncol = 4)
            colnames(curocmat) <- c("position1", "position2", "LLR1", "LLR2")
            for(curoccnnt in 1:length(curoc)){
              curocmat[curoccnnt, 1] <- cur.1st.seq.hits.rev[crow]
              curocmat[curoccnnt, 2] <- cur.2nd.seq.hits.forw[curoc[curoccnnt]]
              curocmat[curoccnnt, 3] <- MotifScoreModifiedList[[cur.seq]][[2]][[motif1Index]][curocmat[curoccnnt, 1], 1]
              curocmat[curoccnnt, 4] <- MotifScoreModifiedList[[cur.seq]][[1]][[motif2Index]][curocmat[curoccnnt, 2], 1]
            }
            revforw <- rbind(revforw, curocmat)
          }
        }
      }
    }
    OccurList[[cur.seq]] <- list()
    if(orient == 1){
      OccurList[[cur.seq]][[1]] <- forwforw
      OccurList[[cur.seq]][[2]] <- revrev
      names(OccurList[[cur.seq]]) <- c("forwforw", "revrev")
    }else if(orient == -1){
      OccurList[[cur.seq]][[1]] <- forwrev
      OccurList[[cur.seq]][[2]] <- revforw
      names(OccurList[[cur.seq]]) <- c("forwrev", "revforw")
    }else{
      OccurList[[cur.seq]][[1]] <- forwforw
      OccurList[[cur.seq]][[2]] <- revrev
      OccurList[[cur.seq]][[3]] <- forwrev
      OccurList[[cur.seq]][[4]] <- revforw
      names(OccurList[[cur.seq]]) <- c("forwforw", "revrev", "forwrev", "revforw")
    }
  }#end of loop over sequences
  return(OccurList)
}
##########################################################################################################################
###################################################################################################################

###################################################################################################################
###################################################################################################################
motifPairInvestigator <- function(MotifScore_OutputList = numeric(0),
                                  .motifList = TF.motifs.TimeSeries.loose.t,
                                  seqList,
                                  LLRthresh = rep(0, length(.motifList)),
                                  plot = T,
                                  pairInvestigateMat = matrix(nrow = 0, ncol = 4)
                                  #direction=T
                                  ){
  # this function gets a list of sequences, and a list of motifs, and index of pairs of motifs and
  # outputs the number, strength and position of co-occurance of those motifs
  # MotifScore_OutputList is a list of outputs of MotifScore function (in case it is precomputed)
  # .motifList is a named list where each entry is a PWM
  # seqList is a list of sequences
  # LLRthresh : a numeric vector: each entry corresponding to one motif in the .motifList. consider hits only stronger than this threshold
  # pairInvestigateMat is a matrix, each row indicates a pair of motifs to investigate in relation to each other. first column: index of first motif (in .motifList), second column: index of second motif, third column: 1 if the have to be on the same strand, -1 if they have to be on opposite strands, and 0 if it doesn't matter, forth column: max distance between the motifs (end of first to start of second)
  # it outputs a list that for each sequence: has one entry for each row of the pairInvestigateMat. that entry is a matrix with number of co-occurances as number of rows and 4 columns. First and second columns: strength(LLR) of first and second motifs, third and forth position of first and second motif hits
  # REMOVED FOR NOW---direction  note that by default direction is considred: the first TF has to occur first, can use the option direction = F to compute both directions
  
  # compute motifScores if they are not already given as input
  if (length(MotifScore_OutputList) > 0){
    MotifScore_OutList <- MotifScore_OutputList
  }else{
    print("computing motif scores")
    seqList <- lapply(X=seqList, FUN=toupper )
    print("computing motif scores")
    MotifScore_OutList <- lapply(X=seqList,
                                 FUN=MotifScore,
                                 bg = c(0.25, 0.25, 0.25, 0.25),
                                 motifList = .motifList)
  }
  #set LLRs less than LLRthresh to zero
  for (cur.seq in 1:length(seqList)){
    for (tf in 1:length(.motifList)){
      MotifScore_OutList[[cur.seq]][[1]][[tf]][MotifScore_OutList[[cur.seq]][[1]][[tf]] < LLRthresh[tf]] <- 0
      MotifScore_OutList[[cur.seq]][[2]][[tf]][MotifScore_OutList[[cur.seq]][[2]][[tf]] < LLRthresh[tf]] <- 0
    }
  }
  pairInteractionList <- list()
  #finding pairs
  print("finding pairs..")
  for(cur_pair in 1:nrow(pairInvestigateMat)){
    print(paste("pair number", cur_pair))
    pairInteractionList[[cur_pair]] <-   findNextOccuring(MotifScoreModifiedList=MotifScore_OutList,
                                              motif1Index=pairInvestigateMat[cur_pair, 1], motif2Index=pairInvestigateMat[cur_pair, 2],
                                              motif1Length=ncol(.motifList[[pairInvestigateMat[cur_pair, 1]]]), motif2Length=ncol(.motifList[[pairInvestigateMat[cur_pair, 2]]]),
                                              orient=pairInvestigateMat[cur_pair, 3],
                                              .seqList=seqList,
                                              distTresh=pairInvestigateMat[cur_pair, 4])
  }
  names(pairInteractionList) <- character(length = nrow(pairInvestigateMat))
  for(cur_pair in 1:nrow(pairInvestigateMat)){
    names(pairInteractionList)[cur_pair] <- paste0(names(.motifList)[pairInvestigateMat[cur_pair, 1]],
                                                   "_",
                                                   names(.motifList)[pairInvestigateMat[cur_pair, 2]])
  }
  return(pairInteractionList)
}
###################################################################################################################
###################################################################################################################
# example
aaSeq <- MCFEnhancers[1872:1873,"sequ"]
aa2 <- motifPairInvestigator(MotifScore_OutputList = numeric(0),
                             .motifList = TF.motifs.TimeSeries.loose.t,
                             seqList=as.list(aaSeq),
                             LLRthresh = c(2, 2),
                             pairInvestigateMat = rbind(c(1, 2, 1 ,5),
                                                        c(2, 1, 1 ,5),
                                                        c(1, 1, -1, 3))
                             )
###################################################################################################################
###################################################################################################################
removeRedundantPairs <- function(motifPairInvestigatorOutput, index = c(1,2)){
  #skip for now
}
###################################################################################################################
###################################################################################################################

EnrichedPairs <- function(motifPairInvestigatorOutput){
  # Ranks the sequences based on the number of each motifpair in their regulatory elements
  
  # motifPairInvestigatorOutput : output of motifPairInvestigator
  # outputs two lists of each three entries, first list: each entry is the index of the 
  # sequences sorted by the number of occurance of each pair
  # and second list is the corresponding motifPairInvestigatorOutput entry (concatanated)
  sortedIndexList <- list()
  unsortedOutput <- list()
  sortedOutput <- list()
  unsortedOutputCount <- list()
  
  for (pair_cu in 1:length(motifPairInvestigatorOutput)){
    print(paste("Going through", names(motifPairInvestigatorOutput)[pair_cu]))
    unsortedOutput[[pair_cu]] <- list()
    unsortedOutputCount[[pair_cu]] <- integer(length = length(motifPairInvestigatorOutput[[pair_cu]]))
    for(seq_cu in 1:length(motifPairInvestigatorOutput[[pair_cu]])){
      #for each sequence concatanate all types of orientations and name the rows correspondingly
      cu_nrows <- unlist(lapply(motifPairInvestigatorOutput[[pair_cu]][[seq_cu]], nrow))
      myRowNames <- character(0)
      for(cur_orient in 1:length(cu_nrows)){
        if(cu_nrows[cur_orient] > 0){
          myRowNames <- c(myRowNames, paste(names(cu_nrows)[cur_orient], c(1:cu_nrows[cur_orient]), sep="_"))
        }
      }
      unsortedOutput[[pair_cu]][[seq_cu]] <- do.call(rbind, motifPairInvestigatorOutput[[pair_cu]][[seq_cu]])
      rownames(unsortedOutput[[pair_cu]][[seq_cu]] ) <- myRowNames
      unsortedOutputCount[[pair_cu]][seq_cu] <- nrow(unsortedOutput[[pair_cu]][[seq_cu]])
    }#end of loop over sequences
    sortedIndexList[[pair_cu]] <- sort(unsortedOutputCount[[pair_cu]], decreasing = T, index.return = T)$ix
    sortedOutput[[pair_cu]] <- unsortedOutput[[pair_cu]][sortedIndexList[[pair_cu]]]
  }#end of loop over pairs
  names(sortedIndexList) <- names(motifPairInvestigatorOutput)
  names(sortedOutput) <- names(motifPairInvestigatorOutput)
  names(unsortedOutput) <- names(motifPairInvestigatorOutput)
  return(list(SortedIndex = sortedIndexList, SortedOutput = sortedOutput, UnsortedOutput = unsortedOutput))
}
###################################################################################################################
###################################################################################################################
#Example
aa <- EnrichedPairs(ER.associated.reg.elements_uniq_pair[[1]])
###################################################################################################################
###################################################################################################################
# next function starts here

###################################################################################################################
##############################################                    #################################################
##############################################      analysis      #################################################
##############################################                    #################################################
###################################################################################################################
# extracting new motifs for ESR1 and RARA:

TF.motifs.TimeSeries.loose.t <- list()

aa <- query(MotifDb, "ESR1")
seqLogo::seqLogo(aa$`Hsapiens-jaspar2016-ESR1-MA0112.2`)
TF.motifs.TimeSeries.loose.t[[1]] <- aa$`Hsapiens-jaspar2016-ESR1-MA0112.2`[, 6:11]
seqLogo::seqLogo(TF.motifs.TimeSeries.loose.t[[1]])

aa <- query(MotifDb, "RARA")
seqLogo::seqLogo(aa$`Hsapiens-HOCOMOCOv10-RARA_HUMAN.H10MO.C`)
TF.motifs.TimeSeries.loose.t[[2]] <- aa$`Hsapiens-HOCOMOCOv10-RARA_HUMAN.H10MO.C`[, 1:6]
seqLogo::seqLogo(TF.motifs.TimeSeries.loose.t[[2]])
names(TF.motifs.TimeSeries.loose.t) <- c("ESR1", "RARA")
#add colnames
for(i in 1:length(TF.motifs.TimeSeries.loose.t)){
  colnames(TF.motifs.TimeSeries.loose.t[[i]]) <- c(1:ncol(TF.motifs.TimeSeries.loose.t[[i]]))
}
#add psuedocount of 0.001
TF.motifs.TimeSeries.loose.t <- lapply(X = TF.motifs.TimeSeries.loose.t, FUN = Addpsudo)
############################################################################################################
# Gather all ER associated regulatory elements
aa  <- which(rownames(GSE78167.RNAseq.Avg.Norm) %in% Genes.Associated.REMAP.ER.Entrez)
aaa <- match(rownames(GSE78167.RNAseq.Avg.Norm)[aa], names(Vicinity100kb.Enhancer.plusPromoter.By.gene))
ER.associated.reg.elements <- WriteFastaOfBag(enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene[aaa],
                Enhancer.Sequence=MCFEnhancers[,"sequ"],
                promoter.promoter.int=Promoter.gene.PromoterIntList[aaa],
                promoter.sequence=Promoter.Sequence.Char,
                return.Length = T,
                returnListNotwrite = T)
aa <- do.call(c, ER.associated.reg.elements$sequence)
aaa <- aa[!duplicated(aa)]
########################################################################################################################
# Run motif scan with ER and RAR motifs
# 5 following lines are commented out to prevent accidents. they have been used.########################
# ER.associated.reg.elements_uniq_MotifScore <- list()
# for(i in 1:length(aaa)){
#   print(i)
#   ER.associated.reg.elements_uniq_MotifScore[[i]] <- MotifScore(seq = aaa[i], bg = c(0.25, 0.25, 0.25, 0.25), motifList = TF.motifs.TimeSeries.loose.t)
# }
########################################################################################################################
# Compute the number of hits of ER-ER (opposite strand), ER-RAR (same strand), and RAR-ER (same strand) with
# differnt thresholds for each unique enhancer and sotre them in a list
## apply the function to find different types of binding site combination to the generated unique motifscore results:
aa <- do.call(c, ER.associated.reg.elements$sequence)
aaa <- aa[!duplicated(aa)]
ER.associated.reg.elements_uniq_pair <- list()
for(i in 1:5){
  ER.associated.reg.elements_uniq_pair[[i]] <- motifPairInvestigator(MotifScore_OutputList = ER.associated.reg.elements_uniq_MotifScore,
                                                                     .motifList = TF.motifs.TimeSeries.loose.t,
                                                                     seqList=aaa,
                                                                     LLRthresh = c(i, i),
                                                                     pairInvestigateMat = rbind(c(1, 2, 1 ,5),
                                                                                                c(2, 1, 1 ,5),
                                                                                                c(1, 1, -1, 3))
                                                                     )
}
names(ER.associated.reg.elements_uniq_pair) <- paste0("LLR_Thresh_", c(1:5))
########################################################################################################################
# rank the unique sequences based on the number of hits for each of the 3 types of motif occurances for each threshold
# create a list with 5 entries for the five LLR thresholds,
# each entry is a list of size 2 with index and concatanated output for each type of motif occurence:
#  ER-RAR, RAR-ER, ER-ER each entry is a sorted list of sequences based on the number of hits for that type
ER.associated.reg.elements_uniq_pair_Concat_sort_un <- list()
for(i in 1:length(ER.associated.reg.elements_uniq_pair)){
  ER.associated.reg.elements_uniq_pair_Concat_sort_un[[i]] <- EnrichedPairs(ER.associated.reg.elements_uniq_pair[[i]])
}
names(ER.associated.reg.elements_uniq_pair_Concat_sort_un) <- names(ER.associated.reg.elements_uniq_pair)
########################################################################################################################
# create a matrix with nrow = number of sequences, ncol=3 (for each pair) for each threshold
ER.associated.reg.elements_uniq_OccuranceMatrixList <- list()
for (i in 1:length(ER.associated.reg.elements_uniq_pair_Concat_sort_un)){
  ER.associated.reg.elements_uniq_OccuranceMatrixList[[i]] <- matrix(nrow = length(ER.associated.reg.elements_uniq_pair_Concat_sort_un$LLR_Thresh_1$SortedIndex$ESR1_RARA),
                                                                     ncol = length(ER.associated.reg.elements_uniq_pair_Concat_sort_un$LLR_Thresh_1$SortedIndex))
  colnames(ER.associated.reg.elements_uniq_OccuranceMatrixList[[i]]) <- names(ER.associated.reg.elements_uniq_pair_Concat_sort_un$LLR_Thresh_1$SortedIndex)
  for(j in 1:ncol(ER.associated.reg.elements_uniq_OccuranceMatrixList[[i]])){
    ER.associated.reg.elements_uniq_OccuranceMatrixList[[i]][, j] <- unlist(lapply(ER.associated.reg.elements_uniq_pair_Concat_sort_un[[i]]$UnsortedOutput[[j]], nrow))
  }
}
names(ER.associated.reg.elements_uniq_OccuranceMatrixList) <- names(ER.associated.reg.elements_uniq_pair_Concat_sort_un)
########################################################################################################################
#sort the sequences based on number of hits of different pairs in different thresholds
aaas <- list()
for(i in 1:5){
  aaas[[i]] <- list()
  for(j in 1:3){
    aaas[[i]][[j]] <- sort(ER.associated.reg.elements_uniq_OccuranceMatrixList[[i]][, j], decreasing = T, index.return = T)$ix
  }
}
########################################################################################################################
# select sequences with no [ER-RAR or RAR-ER] and more than 2 ER-ER
# select sequences with more than 3 [ER-RAR or RAR-ER] and more than 3 ER-ER

########################## select sequences with no [ER-RAR or RAR-ER] and more than 2 ER-ER
aa <- do.call(c, ER.associated.reg.elements$sequence)
aaa <- aa[!duplicated(aa)]

aaxx_d <- which(ER.associated.reg.elements_uniq_OccuranceMatrixList$LLR_Thresh_4[, 3] > 1 &
    ER.associated.reg.elements_uniq_OccuranceMatrixList$LLR_Thresh_4[, 1] < 1 &
    ER.associated.reg.elements_uniq_OccuranceMatrixList$LLR_Thresh_4[, 2] < 1   )

ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence <- aaa[aaxx_d]

########################## select sequences with more than 3 [ER-RAR or RAR-ER] and more than 3 ER-ER
aaxx_u <- which(ER.associated.reg.elements_uniq_OccuranceMatrixList$LLR_Thresh_4[, 3] > 2 &
    (ER.associated.reg.elements_uniq_OccuranceMatrixList$LLR_Thresh_4[, 1] > 2 |
       ER.associated.reg.elements_uniq_OccuranceMatrixList$LLR_Thresh_4[, 2] > 2))

ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence <- aaa[aaxx_u]
####################################find the index of the enhacners in MCFEnhancers[,"sequ"] and promoters in Promoter.Sequence.Char
######separate promoters and enhancers based on their length find their index 

aapr <- which(unlist(lapply(ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence, nchar)) == 1000)
aaenh <- which(unlist(lapply(ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence, nchar)) != 1000)

ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.enh <- ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence[aaenh]
ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.pro <- ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence[aapr]


aapr <- which(unlist(lapply(ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence, nchar)) == 1000)
aaenh <- which(unlist(lapply(ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence, nchar)) != 1000)

ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.enh <- ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence[aaenh]
ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.pro <- ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence[aapr]

########################Find enhancers in MCFEnhancers[,"sequ"]
ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.enh.index <- integer(length(ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.enh))
for(i in 1:length(ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.enh)){
  ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.enh.index[i] <- which(toupper(MCFEnhancers[,"sequ"]) %in% ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.enh[i])
}
ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.enh.index <- integer(length(ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.enh))
for(i in 1:length(ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.enh)){
  ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.enh.index[i] <- which(toupper(MCFEnhancers[,"sequ"]) %in% ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.enh[i])
}
#########################Find promoters in Promoter.Sequence.Char

ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.pro.index <- integer(length = length(ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.pro))
for(i in 1:length(ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.pro)){
  ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.pro.index[i] <- which(toupper(Promoter.Sequence.Char) %in% ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.pro[i])
}

ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.pro.index <- integer(length = length(ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.pro))



######################################### find out the enrichment of ER and RAR chip peaks in each group of sequences #########################
Promoter.gene.RAR.ER.Interaction.Overlap.byPromoter.Mat # matirx of promoter interaction or overlap with ER-RAR ReChIP
Enhancer.gene.RAR.ER.Interaction.Overlap.byEnhancer.Mat # matirx of enhancer interaction or overlap with ER-RAR ReChIP
Genes.Associated.RAR.ER.Entrez # gene Entrez names associated with ER-RAR ReChIP

##################### Enrichment test for ER-RAR associated enhancers in the ERERgt3.ERRARgt3 enhacner set ######################
sum(which(Enhancer.gene.RAR.ER.Interaction.Overlap.byEnhancer.Mat == 1) %in% ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.enh.index)
#hypergeometric test: 33 overlap, 224 total ER-RAR associated enhancers, 8189 total ER associated enhancers that are not ER-RAR associated, 135 size of the motif deerived set
aa <- phyper(33, 224, 8189, 343) #pval = 2.0139e-11
######################################################################################################### 

##################### Enrichment test for ER-RAR associated enhancers in the ERERgt2.ERRAReq0 enhacner set ######################
sum(which(Enhancer.gene.RAR.ER.Interaction.Overlap.byEnhancer.Mat == 1) %in% ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.enh.index)
#hypergeometric test: 1 overlap, 224 total ER-RAR associated enhancers, 8189 total ER associated enhancers that are not ER-RAR associated, 22 size of the motif deerived set
aa <- phyper(1, 224, 8189, 94) #pval = 0.7191352

######################################################################################################### 

##################### Enrichment test for ER-RAR associated genes in the ERERgt3.ERRARgt3 gene set ######################
sum(Genes.Associated.RAR.ER.Entrez %in% ER.associated.genes.ERERgt3.ERRARgt3)
#hypergeometric test: 85 overlap, 336 total ER-RAR associated genes in my starting set, 5483 total ER associated genes that are not ER-RAR associated, 319 size of the motif deerived set
aa <- phyper(146, 336, 5483, 724) # pval = 0
######################################################################################################### 
intersect(Genes.Associated.RAR.ER.Entrez,ER.associated.genes.ERERgt2.ERRAReq0)
#hypergeometric test: 6 overlap, 336 total ER-RAR associated genes in my starting set, 5483 total ER associated genes that are not ER-RAR associated, 48 size of the motif deerived set
aa <- phyper(10, 336, 5483, 305) # pval = 0.9704475


########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
########## ER-RAR associated genes (Genes.Associated.RAR.ER.Entrez) are enriched in the ERERgt3.ERRARgt3 gene set (ER.associated.genes.ERERgt3.ERRARgt3) ############################## 
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 


######################################### find genes associated with each group of enhancers ########################################

#################### ERERgt2.ERRAReq0
aaenh <- matrix(0L, nrow = length(MCFEnhancers[,"sequ"]), ncol = 1)
aaenh[ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.enh.index, 1] <- 1
colnames(aaenh) <- "ERERgt2ERRAReq0"
aapro <- matrix(0L, nrow = length(Promoter.Sequence.Char), ncol = 1)
aapro[ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence.pro.index, 1] <- 1
colnames(aapro) <- "ERERgt2ERRAReq0"

ER.associated.genes.ERERgt2.ERRAReq0 <- ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList,
                                                             gene.enhancer = Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                                             enhancerChipBinary = aaenh,
                                                             promoterChipBinary = aapro)
ER.associated.genes.ERERgt2.ERRAReq0 <- rownames(ER.associated.genes.ERERgt2.ERRAReq0$AssociatedPeaksBinary)[which(ER.associated.genes.ERERgt2.ERRAReq0$AssociatedPeaksBinary == 1)]
#################### ERERgt3.ERRARgt3
aaenh <- matrix(0L, nrow = length(MCFEnhancers[,"sequ"]), ncol = 1)
aaenh[ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.enh.index, 1] <- 1
colnames(aaenh) <- "ERERgt3ERRARgt3"
aapro <- matrix(0L, nrow = length(Promoter.Sequence.Char), ncol = 1)
aapro[ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence.pro.index, 1] <- 1
colnames(aapro) <- "ERERgt3ERRARgt3"

ER.associated.genes.ERERgt3.ERRARgt3 <- ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList,
                                                             gene.enhancer = Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                                             enhancerChipBinary = aaenh,
                                                             promoterChipBinary = aapro)
ER.associated.genes.ERERgt3.ERRARgt3 <- rownames(ER.associated.genes.ERERgt3.ERRARgt3$AssociatedPeaksBinary)[which(ER.associated.genes.ERERgt3.ERRARgt3$AssociatedPeaksBinary == 1)]

########################################################################################################################
#Choose genes for GEMSTAT modeling:
length(intersect(rownames(GSE78167.RNAseq.Diff), ER.associated.genes.ERERgt2.ERRAReq0))
length(intersect(rownames(GSE78167.RNAseq.Diff), ER.associated.genes.ERERgt3.ERRARgt3))
# 22 genes in : 
ER.associated.genes.ERERgt2.ERRAReq0.GEMSTAT <- setdiff(intersect(rownames(GSE78167.RNAseq.Diff), ER.associated.genes.ERERgt2.ERRAReq0), Genes.Associated.RAR.ER.Entrez)
ER.associated.genes.ERERgt2.ERRAReq0.GEMSTAT <- ER.associated.genes.ERERgt2.ERRAReq0.GEMSTAT[-c(13,17,19)]
# 22 genes in : 
ER.associated.genes.ERERgt3.ERRARgt3.GEMSTAT <- intersect(intersect(ER.associated.genes.ERERgt3.ERRARgt3, GSE78167.RNAseq.Diff) , Genes.Associated.RAR.ER.Entrez)

# look at their expression pattern in the data set : are these differentially expressed? how many enhancers do each have? 
GSE78167.RNAseq.Diff #differentially expressed genes exp mat
########################################################################################################################

# ###map back from the unique entries to the original gene-based list 
# 
# for(gene_nu in 1:length(ER.associated.reg.elements$sequence)){
#   ER.associated.reg.elements_pergene_pairSites[[gene_nu]] <- list()
#   for(reg_el in 1:length(ER.associated.reg.elements$sequence[[gene_nu]])){
#     unique_index <- which(aaa %in% ER.associated.reg.elements$sequence[[gene_nu]][reg_el])
#     ER.associated.reg.elements_pergene_pairSites[[gene_nu]][[reg_el]] <- list(ESR1_RARA=ER.associated.reg.elements_uniq_pair[[1]][[unique_index]],
#                                                                               RARA_ESR1=ER.associated.reg.elements_uniq_pair[[2]][[unique_index]],
#                                                                               ESR1_ESR1=ER.associated.reg.elements_uniq_pair[[3]][[unique_index]])
#   }
#   
# }
# names(ER.associated.reg.elements_pergene_pairSites) <- names(ER.associated.reg.elements$sequence)
# #######################

################################################# Scan the two groups for ER-RAR and ER-only for AP1 motifs
################################################ There are two variants of this motif to be tested
TF.motifs.TimeSeries.AP1.t <- list()
# AP1 motif first variant
TF.motifs.TimeSeries.AP1.t[[1]] <- t(RemoveDegenRow(TF.motifs$JUN_1, background=c(0.25, 0.25, 0.25, 0.25), threshold=0.3))

# AP1 motif second variant
TF.motifs.TimeSeries.AP1.t[[2]] <- t(RemoveDegenRow(TF.motifs$JUN_2, background=c(0.25, 0.25, 0.25, 0.25), threshold=0.3))

names(TF.motifs.TimeSeries.AP1.t) <- c("JUN_1","JUN_2")
#sequences to be explored:
ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence
ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence

ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore <- list()
ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore[[1]] <- list()
ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore[[2]] <- list()
names(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore) <- c("ERERgt2_ERRARgt2", "ERERgt1.ERRAReq0")

for(i in 1:length(ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence)){
 print(i)
  ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore[[1]][[i]] <- MotifScore(seq = ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence[i],
                                                                              bg = c(0.25, 0.25, 0.25, 0.25),
                                                                              motifList = TF.motifs.TimeSeries.AP1.t)
}

for(i in 1:length(ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence)){
  print(i)
  ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore[[2]][[i]] <- MotifScore(seq = ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence[i],
                                                                                   bg = c(0.25, 0.25, 0.25, 0.25),
                                                                                   motifList = TF.motifs.TimeSeries.AP1.t)
}
ER.associated.reg.elements.ERERRAR.Filtered.sequence <- list(ERgt1RARgt1 = ER.associated.reg.elements.ERERgt3.ERRARgt3.sequence, ERgt2RAReq0=ER.associated.reg.elements.ERERgt2.ERRAReq0.sequence)
#######extracting motif hits given a threshold
#default threshold: max(LLR) - 4
ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits <- list()
for(i in 1:length(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore)){
  ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits[[i]] <- list()
  for (j in 1:length(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore[[i]])){
    ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits[[i]][[j]] <- ExtractMotifHits(seq = ER.associated.reg.elements.ERERRAR.Filtered.sequence[[i]][j],
                                                                                          motifList = TF.motifs.TimeSeries.AP1.t,
                                                                                          FinalScore = ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore[[i]][[j]],
                                                                                          threshold=exp(-4))
  }
}
names(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits) <- c("ERgt1RARgt1", "ERgt2RAReq0")

#another threshold: max(LLR) - 5
ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.LLR_5 <- list()
for(i in 1:length(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore)){
  ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.LLR_5[[i]] <- list()
  for (j in 1:length(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore[[i]])){
    ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.LLR_5[[i]][[j]] <- ExtractMotifHits(seq = ER.associated.reg.elements.ERERRAR.Filtered.sequence[[i]][j],
                                                                                          motifList = TF.motifs.TimeSeries.AP1.t,
                                                                                          FinalScore = ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifScore[[i]][[j]],
                                                                                          threshold=exp(-5))
  }
}
names(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.LLR_5) <- c("ERgt1RARgt1", "ERgt2RAReq0")



####################################################################################################
# create a list with two entries, each a matrix with nrow: num sequences in that group, ncol = 2 for the two motifs, each entry is the number of motif hits for that sequence
ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList <- list()
ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList[[1]] <- matrix(nrow = length(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits[[1]]), ncol = 2)
ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList[[2]] <- matrix(nrow = length(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits[[2]]), ncol = 2)

aa <- lapply(X = ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits$ERgt1RARgt1,"[[", 3)
for(i in 1:length(aa)){
  for(j in 1:2)
  ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList[[1]][i,j] <- length(unique(aa[[i]][[j]]))
}
aa <- lapply(X = ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits$ERgt2RAReq0,"[[", 3)
for(i in 1:length(aa)){
  for(j in 1:2)
    ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList[[2]][i,j] <- length(unique(aa[[i]][[j]]))
}

names(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList) <- c("ERgt1RARgt1", "ERgt2RAReq0")
#### same for -5 threshold

ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5 <- list()
ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5[[1]] <- matrix(nrow = length(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.LLR_5[[1]]), ncol = 2)
ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5[[2]] <- matrix(nrow = length(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.LLR_5[[2]]), ncol = 2)

aa <- lapply(X = ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.LLR_5$ERgt1RARgt1,"[[", 3)
for(i in 1:length(aa)){
  for(j in 1:2)
    ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5[[1]][i,j] <- length(unique(aa[[i]][[j]]))
}
aa <- lapply(X = ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.LLR_5$ERgt2RAReq0,"[[", 3)
for(i in 1:length(aa)){
  for(j in 1:2)
    ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5[[2]][i,j] <- length(unique(aa[[i]][[j]]))
}

names(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5) <- c("ERgt1RARgt1", "ERgt2RAReq0")


############ Removing enhancers with length more than 10kb
which(unlist(lapply(ER.associated.reg.elements.ERERRAR.Filtered.sequence$ERgt1RARgt1, nchar)) < 10000)
sum(unlist(lapply(ER.associated.reg.elements.ERERRAR.Filtered.sequence$ERgt2RAReq0, nchar)) > 10000)

ER.associated.reg.elements.ERERRAR.Filtered.sequence.LT10KB <- ER.associated.reg.elements.ERERRAR.Filtered.sequence
ER.associated.reg.elements.ERERRAR.Filtered.sequence.LT10KB$ERgt1RARgt1 <- ER.associated.reg.elements.ERERRAR.Filtered.sequence.LT10KB$ERgt1RARgt1[which(unlist(lapply(ER.associated.reg.elements.ERERRAR.Filtered.sequence$ERgt1RARgt1, nchar)) < 10000)]

ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LT10KB <- ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList
ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LT10KB$ERgt1RARgt1 <- ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LT10KB$ERgt1RARgt1[which(unlist(lapply(ER.associated.reg.elements.ERERRAR.Filtered.sequence$ERgt1RARgt1, nchar)) < 10000), ]

ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5.LT10KB <- ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5
ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5.LT10KB$ERgt1RARgt1 <- ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5.LT10KB$ERgt1RARgt1[which(unlist(lapply(ER.associated.reg.elements.ERERRAR.Filtered.sequence$ERgt1RARgt1, nchar)) < 10000), ]


boxplot.matrix(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LT10KB$ERgt1RARgt1)
boxplot.matrix(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LT10KB$ERgt2RAReq0)

boxplot.matrix(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5.LT10KB$ERgt1RARgt1)
boxplot.matrix(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5.LT10KB$ERgt2RAReq0)
######################################################################################################

sum(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5.LT10KB$ERgt1RARgt1[,2] > 0)
sum(ER.associated.reg.elements_ERRAR_Filtered_AP1_MotifHits.MatList.LLR_5.LT10KB$ERgt2RAReq0[,2] > 0)

