# ER_eqtl analysis
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(IRanges)
library(ChIPseeker)
library(dplyr)
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
#########################################################################################################
#########################################################################################################
##########################################     functions      ###########################################
#########################################################################################################
#########################################################################################################

# write a function to convert interaction data between genome versions
interaction_liftover <- function(int_a_gr, int_b_gr, my_chain){
  # int_a_gr is a Granges object for the first side of interactions
  # int_b_gr is a Granges object for the second side of interactions
  # my_chain : is the chain file for liftover
  library(rtracklayer)
  stopifnot(length(int_a_gr) == length(int_b_gr))
  int_a_gr_df <- data.frame(int_a_gr)
  int_b_gr_df <- data.frame(int_b_gr)
  int_a_gr_df$unq_id <- c(1:nrow(int_a_gr_df))
  int_b_gr_df$unq_id <- c(1:nrow(int_b_gr_df))
  int_a_gr <- makeGRangesFromDataFrame(int_a_gr_df,keep.extra.columns = T)
  int_b_gr <- makeGRangesFromDataFrame(int_b_gr_df,keep.extra.columns = T)
  int_a_gr_new <- unlist(liftOver(x = int_a_gr, chain = my_chain))
  int_b_gr_new <- unlist(liftOver(x = int_b_gr, chain = my_chain))
  
  dupl_index_a <- unique(int_a_gr_new$unq_id[duplicated(int_a_gr_new$unq_id)])
  dupl_index_b <- unique(int_b_gr_new$unq_id[duplicated(int_b_gr_new$unq_id)])
  
  # just initializing, at the end will remove the first entry
  dup_new_grange_a <- int_a_gr_new[1]
  dup_new_grange_b <- int_b_gr_new[1]
  for(i in 1:length(dupl_index_a)){
    cur_a_ind <- which(int_a_gr_new$unq_id %in% dupl_index_a[i])
    if(dupl_index_a[i] %in% dupl_index_b){
      cur_b_ind <- which(int_b_gr_new$unq_id %in% dupl_index_a[i])
      for(j in 1:length(cur_b_ind)){
        cur_b_grange <- rep(int_b_gr_new[cur_b_ind[j]], length(cur_a_ind))
        cur_a_grange <- int_a_gr_new[cur_a_ind]
        dup_new_grange_a <- c(dup_new_grange_a, cur_a_grange)
        dup_new_grange_b <- c(dup_new_grange_b, cur_b_grange)
      }
    }else{
      cur_b_grange <- rep(int_b_gr_new[int_b_gr_new$unq_id %in% dupl_index_a[i]], length(cur_a_ind))
      cur_a_grange <- int_a_gr_new[cur_a_ind]
      dup_new_grange_a <- c(dup_new_grange_a, cur_a_grange)
      dup_new_grange_b <- c(dup_new_grange_b, cur_b_grange)
      
    }
  }
  
  dupl_index_b_not_a <- setdiff(dupl_index_b, dupl_index_a)
  # go over the ones that are in b but not a
  for(i in 1:length(dupl_index_b_not_a)){
    cur_b_ind <- which(int_b_gr_new$unq_id %in% dupl_index_b_not_a[i])
    cur_b_grange <- int_b_gr_new[cur_b_ind]
    cur_a_grange <- rep(int_a_gr_new[int_a_gr_new$unq_id %in% dupl_index_b_not_a[i]], length(cur_b_ind))
    dup_new_grange_a <- c(dup_new_grange_a, cur_a_grange)
    dup_new_grange_b <- c(dup_new_grange_b, cur_b_grange)
  }
  dup_new_grange_a <- dup_new_grange_a[2:length(dup_new_grange_a)]
  dup_new_grange_b <- dup_new_grange_b[2:length(dup_new_grange_b)]
  stopifnot(length(dup_new_grange_a) == length(dup_new_grange_b), 
            all(dup_new_grange_a$unq_id == dup_new_grange_b$unq_id))
  dupl_index_all <- union(dupl_index_a, dupl_index_b)
  intersect_new <- intersect(int_a_gr_new$unq_id, int_b_gr_new$unq_id)
  uniq_data_a <- int_a_gr_new[! int_a_gr_new$unq_id %in% dupl_index_all]
  uniq_data_a <- uniq_data_a[uniq_data_a$unq_id %in% intersect_new]
  uniq_data_a <- uniq_data_a[sort(uniq_data_a$unq_id, decreasing = F, index.return=T)$ix]
  uniq_data_b <- int_b_gr_new[! int_b_gr_new$unq_id %in%dupl_index_all]
  uniq_data_b <- uniq_data_b[uniq_data_b$unq_id %in% intersect_new]
  uniq_data_b <- uniq_data_b[sort(uniq_data_b$unq_id, decreasing = F, index.return=T)$ix]
  stopifnot(length(uniq_data_a) == length(uniq_data_b),
            all(uniq_data_b$unq_id == uniq_data_a$unq_id))
  
  all_new_A <- c(uniq_data_a, dup_new_grange_a)
  all_new_B <- c(uniq_data_b, dup_new_grange_b)
  return(list(new_A = all_new_A, new_B = all_new_B))
  
}
#########################################################################################################
# example
aa_1 <-MCF7_interaction_hg19[,1:3] 
colnames(aa_1) <- c("chr", "start", "end")
#aa_1$reg_id <- c(1:nrow(aa_1))
aa_2 <-MCF7_interaction_hg19[,4:6] 
colnames(aa_2) <- c("chr", "start", "end")
#aa_2$reg_id <- c(1:nrow(aa_2))
aa_1 <- makeGRangesFromDataFrame(aa_1,keep.extra.columns = T)
aa_2 <- makeGRangesFromDataFrame(aa_2,keep.extra.columns = T)
aa_19_39_chain <- import.chain("~/Documents/Shayan/BioInf/Liftover_chain_files/hg19ToHg38.over.chain")

aatest_gr <- interaction_liftover(int_a_gr = aa_1, int_b_gr = aa_2, my_chain = aa_19_39_chain)


#########################################################################################################
write_variant_seq <- function(enhancer_GR, 
                              eqtl_GR,
                              eqtl_names=character(0),
                              eqtl_ref_alt, 
                              all_combs=F, 
                              my_genome=BSgenome.Hsapiens.UCSC.hg38,
                              label = numeric(0)){
  # finds the overlap between ehnancer and eqtl sequence
  # check if the reference is correct
  # change the enhancer that overlaps with an eqtl and write the sequence 
  
  # enhancer_GR: GRanges object containing the coordinates of enhancers
  # eqtl_GR: GRanges object containing the coordinates of eqtls
  # eqtl_names: names for each eqtl or snp
  # eqtl_ref_alt: two column dataframe/matrix. number of rows should be equal to the length of eqtl_GR, 
  #  first column is the reference allele, second col is the alternative allele (if there are more than one alternative alleles they should be separated by ",")
  # all_combs: logical. if True all combinations of the eqtls on an enhancer will be created, if False 
  #  eqtls will be incorporated only one at a time.
  # my_genome: is the genome used for extracting the sequences
  # label : optional. numeric vector of length equal to length(enhancer_GR).
  if(!dir.exists("Variant_seq")){
    system("mkdir Variant_seq")
  }
  if(!dir.exists("Variant_Labels")){
    system("mkdir Variant_Labels")
  }
  library(GenomicRanges)
  library(Biostrings)
  library(ShortRead)
  
  stopifnot(nrow(eqtl_ref_alt) == length(eqtl_GR),
            (length(label) == 0 | length(label) == length(enhancer_GR)),
            length(names(enhancer_GR)) == length(enhancer_GR))
  
  alter_length <- unlist(lapply(strsplit(eqtl_ref_alt[, 2], split = ","), length))
  alter_length[alter_length == 0] <- 1
  if(length(eqtl_names) == 0){
    eqtl_names <- c(1:length(eqtl_GR))
  }
  overlap_holder <- findOverlaps(query = enhancer_GR, subject = eqtl_GR)
  ol_enhancers_index <- unique(overlap_holder@from)
  wt_seq_all <- character(length(ol_enhancers_index))
  names(wt_seq_all) <- names(enhancer_GR)[ol_enhancers_index]
  var_seq_list <- list()
  if(length(label) > 0){
    wt_label <- label[ol_enhancers_index]
  }
  
  for(cur_enh in 1:length(ol_enhancers_index)){
     print("enh ind")
     print(ol_enhancers_index[cur_enh])
    cur_enh_eqtl_ind <- overlap_holder@to[overlap_holder@from %in% ol_enhancers_index[cur_enh]]
    cur_enh_eqtl_name <- eqtl_names[cur_enh_eqtl_ind]
    wt_seq_all[cur_enh] <- getSeq(x = my_genome, names = enhancer_GR[ol_enhancers_index[cur_enh]])

    wt_seq_char <- as.character(wt_seq_all[cur_enh])
    stopifnot(nchar(wt_seq_char) == 1000)
    var_seq_list[[cur_enh]] <- character(length = sum(alter_length[cur_enh_eqtl_ind]))
    names(var_seq_list[[cur_enh]]) <- character(length = sum(alter_length[cur_enh_eqtl_ind]))
    cur_cnt <- 1
    for(cur_qtl in 1:length(cur_enh_eqtl_ind)){
       print("eqtl ind")
       print(cur_enh_eqtl_ind[cur_qtl])
      eqtl_pos <- eqtl_GR[cur_enh_eqtl_ind[cur_qtl]]@ranges@start #- 1
      start_pos <- enhancer_GR[ol_enhancers_index[cur_enh]]@ranges@start
      seq_pos <- eqtl_pos - start_pos + 1
      if(seq_pos <= 0 | (seq_pos + nchar(eqtl_ref_alt[cur_enh_eqtl_ind[cur_qtl], 1]) - 1) > 1000){
        print("position 0")
        next()
      }
      ref_check <- substr(x = wt_seq_char, 
                          start = seq_pos,
                          stop = (seq_pos + nchar(eqtl_ref_alt[cur_enh_eqtl_ind[cur_qtl], 1]) - 1))
      # check if the reference seq is correct
       print("seq_pos_start")
       print(seq_pos)
       print("seq_pos_end")
       print((seq_pos + nchar(eqtl_ref_alt[cur_enh_eqtl_ind[cur_qtl], 1]) - 1))
       print("obtained ref seq")
       print(ref_check)
       print("provided ref seq")
       print(eqtl_ref_alt[cur_enh_eqtl_ind[cur_qtl], 1])
     # if(eqtl_ref_alt[cur_enh_eqtl_ind[cur_qtl], 1] != "null"){
       stopifnot(ref_check == eqtl_ref_alt[cur_enh_eqtl_ind[cur_qtl], 1])
     # }

      #stopifnot(cur_enh <= 5)
      my_alt <- unlist(strsplit(eqtl_ref_alt[cur_enh_eqtl_ind[cur_qtl], 2], split = ","))
      if(length(my_alt) == 0){
       my_alt <- ""
      }
      print("my_alt")
      print(my_alt)
      print(alter_length[cur_enh_eqtl_ind[cur_qtl]])
      for(cur_Alt in 1:alter_length[cur_enh_eqtl_ind[cur_qtl]]){
        tmp_Seq <- wt_seq_char
        if(nchar(eqtl_ref_alt[cur_enh_eqtl_ind[cur_qtl], 1]) == nchar(my_alt[cur_Alt])){
          substr(x = tmp_Seq, start = seq_pos, stop = (seq_pos + nchar(eqtl_ref_alt[cur_enh_eqtl_ind[cur_qtl], 1]) - 1)) <- my_alt[cur_Alt]
          var_seq_list[[cur_enh]][cur_cnt] <- tmp_Seq
          names(var_seq_list[[cur_enh]])[cur_cnt] <- paste0("eqtl_",eqtl_names[cur_enh_eqtl_ind[cur_qtl]], "_", cur_Alt)
          cur_cnt <- cur_cnt + 1
        }else{ # if length of ref and alt are not equal
          peice_1 <- substr(x = tmp_Seq, 
                            start = 1, 
                            stop = seq_pos - 1)
          peice_2 <- substr(x = tmp_Seq,
                            start = (seq_pos + nchar(eqtl_ref_alt[cur_enh_eqtl_ind[cur_qtl], 1])),
                            stop = nchar(tmp_Seq))
          tmp_Seq <- paste0(peice_1, my_alt[cur_Alt], peice_2)
          var_seq_list[[cur_enh]][cur_cnt] <- tmp_Seq
          names(var_seq_list[[cur_enh]])[cur_cnt] <- paste0("eqtl_",eqtl_names[cur_enh_eqtl_ind[cur_qtl]], "_", cur_Alt)
          cur_cnt <- cur_cnt + 1
        }
        
        
      }
    }
    print("nchar")
    print(nchar(var_seq_list[[cur_enh]]))
    writeFasta(DNAStringSet(var_seq_list[[cur_enh]]), width = max(nchar(var_seq_list[[cur_enh]])),
               file = paste0("Variant_seq/seq_var_",names(enhancer_GR)[ol_enhancers_index[cur_enh]],".fa"))
    if(length(label) > 0){
      cat(c("Rows", paste0("1", "\n")), sep = "\t",
          append = F, file = paste0("Variant_Labels/Label_var_",
                                    names(enhancer_GR)[ol_enhancers_index[cur_enh]],".txt"))
      for(i in 1:length(var_seq_list[[cur_enh]])){
        cat(names(var_seq_list[[cur_enh]])[i], paste0(wt_label[cur_enh],
                                                            rep("\n", as.integer(i != length(var_seq_list[[cur_enh]])))), 
            append = T
            , file = paste0("Variant_Labels/Label_var_",
                            names(enhancer_GR)[ol_enhancers_index[cur_enh]],".txt"), sep = "\t" )
      }
    }
  }
  writeFasta(DNAStringSet(wt_seq_all), width = 1000,
             file = paste0("Variant_seq/seq_WT.fa"))
  if(length(label) > 0){
    cat(c("Rows", paste0("1", "\n")), sep = "\t",
        append = F, file = "Variant_Labels/Label_WT.txt")
    for(i in 1:length(wt_label)){
      cat(names(enhancer_GR)[ol_enhancers_index[i]], paste0(wt_label[i],
                                                          rep("\n", as.integer(i != length(wt_label)))), 
          append = T
          , file = "Variant_Labels/Label_WT.txt", sep = "\t" )
    }
  }
  return(list(wt_seq = wt_seq_all, 
              var_seq= var_seq_list))
}

###############################################################################################################
# example
ER_gtex_eqtl_sig <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/Breast_Mammary_Tissue.v8.signif_variant_gene_pairs.txt", 
                               stringsAsFactors = F, header = T)
aa_eqtl <-  strsplit(ER_gtex_eqtl_sig$variant_id, split= "_")
aa_eqtl_chr <- unlist(lapply(aa_eqtl, "[[", 1))
aa_eqtl_st <- unlist(lapply(aa_eqtl, "[[", 2))
aa_eqtl_ref <- unlist(lapply(aa_eqtl, "[[", 3))
aa_eqtl_alt <- unlist(lapply(aa_eqtl, "[[", 4))

aa_ch <- nchar(aa_eqtl_ref)
which(aa_ch > 1)[1:5]
aa_eqtl[which(aa_ch > 1)[1:5]]

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_eqtl_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/sig_eqtl.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_qtl <- readPeakFile(aa_eqtl_Add, as = "GRanges")
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_15")
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)
aa_pos_neg_vars <- write_variant_seq(enhancer_GR = aa_pos_neg, 
                                     eqtl_GR = aa_qtl,
                                     eqtl_ref_alt = cbind(aa_eqtl_ref, aa_eqtl_alt),
                                     all_combs = F,
                                     my_genome = BSgenome.Hsapiens.UCSC.hg38, 
                                     label = c(rep(1, length(aa_pos)), rep(0, length(aa_neg))))

##############
# test cases
aa_pos_neg[1]

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/")
aa_qtl2 <- aa_snp[1:6]
aa_eqtl_ref <- character(6)
aa_eqtl_alt <- character(6)

start(aa_qtl2)[1] <- start(aa_pos_neg)[1] 
end(aa_qtl2)[1] <- start(aa_pos_neg)[1] 
aa_eqtl_ref[1] <- "A"
aa_eqtl_alt[1] <- "G,ATT"

start(aa_qtl2)[2] <- start(aa_pos_neg)[1] 
end(aa_qtl2)[2] <- start(aa_pos_neg)[1] 
aa_eqtl_ref[2] <- "AGC"
aa_eqtl_alt[2] <- "C"

end(aa_qtl2)[3] <- end(aa_pos_neg)[1] 
start(aa_qtl2)[3] <- end(aa_pos_neg)[1] 
aa_eqtl_ref[3] <- "C"
aa_eqtl_alt[3] <- "G,CGG"


end(aa_qtl2)[4] <- end(aa_pos_neg)[1] - 2
start(aa_qtl2)[4] <- end(aa_pos_neg)[1] - 2
aa_eqtl_ref[4] <- "CAC"
aa_eqtl_alt[4] <- "C"

start(aa_qtl2)[5] <- start(aa_pos_neg)[1] + 8
end(aa_qtl2)[5] <- start(aa_pos_neg)[1] + 8
aa_eqtl_ref[5] <- "G"
aa_eqtl_alt[5] <- "C,GTT"

start(aa_qtl2)[6] <- start(aa_pos_neg)[1] + 8
end(aa_qtl2)[6] <- start(aa_pos_neg)[1] + 8
aa_eqtl_ref[6] <- "GCC"
aa_eqtl_alt[6] <- "G"

aatester <- write_variant_seq(enhancer_GR = aa_pos_neg[1], 
                              eqtl_GR = aa_qtl2,
                              eqtl_ref_alt = cbind(aa_eqtl_ref, aa_eqtl_alt),
                              all_combs = F,
                              my_genome = BSgenome.Hsapiens.UCSC.hg38, 
                              label = c(rep(1, length(aa_pos[1]))))



######################################################################################################
# write a function to find out the TF responsible for a snp effect
# find the snp on sequence
# for each TF evaluate LLR at each position for both the WT and variant seq
# mark sites above thresh on each seq
# store changes for all points that are marked as BS in either of the sequences
# store gain or loss of a binding site for a TF

snp_investigator <- function(my_genome,
                             my_motifs,
                             snp_chr,
                             snp_start,
                             snp_ref,
                             snp_alt,
                             min_LLR){
  # my_genome : is the genome being used
  # my_motifs : is a list, where each entry is a motifs to scan the sequences, each row is 
  #  a base, each column is a position
  # snp_chr : is a character indicating the chromosome name : eg "chr1"
  # snp_start : is a numeric indicating the start position of the snp
  # snp_ref : is a character indicating the original sequence affected by the snp
  # snp_alt : is a character replacing the original sequence
  # min_LLR : is a named numeric vector containing the minimum LLR for each TF, 
  #  any site with LLR above this threshold is considered a site.

  stopifnot(length(names(my_motifs)) == length(my_motifs),
            length(my_motifs) == length(min_LLR), 
            all(names(min_LLR) == names(my_motifs)),
            is.character(snp_chr),
            length(snp_ref) == 1,
            length(snp_alt)==1, 
            all(unlist(lapply(my_motifs, nrow)) %in% 4))
  max_len <- max(unlist(lapply(my_motifs, ncol)))
  my_start <- snp_start - max_len + 1
  my_end <- snp_start + nchar(snp_ref) - 1 + max_len #- 1
  my_GR <- makeGRangesFromDataFrame(data.frame(chr=snp_chr, 
                                               start = my_start,
                                               end =my_end ))
  ref_Seq <- getSeq(x = my_genome, names =my_GR, as.character=T)
  ref_check <- substr(x = ref_Seq, 
                      start = max_len,
                      stop = max_len + nchar(snp_ref) - 1)
  print("ref_check")
  print(ref_check)
  print("snp_ref")
  print(snp_ref)
  stopifnot(ref_check == snp_ref)
  
  alt_Seq <- ref_Seq
  if(nchar(snp_ref) == nchar(snp_alt)){
    substr(x = alt_Seq, start = max_len, stop = max_len + nchar(snp_ref) - 1) <- snp_alt
  }else{
    peice_1 <- substr(x = alt_Seq, start = 1, stop = max_len - 1)
    peice_2 <- substr(x = alt_Seq, start = max_len + nchar(snp_ref), stop = nchar(alt_Seq))
    alt_Seq <- paste0(peice_1, snp_alt, peice_2)
    # print("snp_ref")
    # print(snp_ref)
    # print("snp_alt")
    # print(snp_alt)
    # print("ref_Seq")
    # print(ref_Seq)
    # print("alt_Seq")
    # print(alt_Seq)
    # print("max_len")
    # print(max_len)
  }
  
  stopifnot(ref_Seq != alt_Seq)
  # check if RevComp is actually working
  # ref_revcomp <- as.character(reverse(complement(DNAString(ref_Seq))))
  # alt_revcomp <- as.character(reverse(complement(DNAString(alt_Seq))))
  motifScore_ref <- MotifScore(seq = ref_Seq,
                               bg = c(0.25, 0.25, 0.25, 0.25),
                               motifList = my_motifs)
  # motifScore_ref_rc <- MotifScore(seq = ref_revcomp,
  #                              bg = c(0.25, 0.25, 0.25, 0.25),
  #                              motifList = my_motifs)
  motifScore_alt <- MotifScore(seq = alt_Seq,
                               bg = c(0.25, 0.25, 0.25, 0.25),
                               motifList = my_motifs)
  # motifScore_alt_rc <- MotifScore(seq = alt_revcomp,
  #                              bg = c(0.25, 0.25, 0.25, 0.25),
  #                              motifList = my_motifs)
  #  my_results <- list()
  
  my_results <- matrix(nrow = 0, ncol = 10)
  colnames(my_results) <- c("TF","old","new","change", "strand", "loss_or_gain", "affected_pos_old", "affected_pos_new", "ref", "alt" )
  # data.frame(TF=character(0), 
  #                               change = numeric(0),
  #                               strand = character(0), 
  #                               loss_or_gain=numeric(0), 
  #                               stringsAsFactors = F)
  for(cur_tf in 1:length(my_motifs)){
    ref_sc_f  <- motifScore_ref$Forward[[cur_tf]][, 1]
    ref_sc_r  <- motifScore_ref$Reverse[[cur_tf]][, 1]
    alt_sc_f  <- motifScore_alt$Forward[[cur_tf]][, 1]
    alt_sc_r  <- motifScore_alt$Reverse[[cur_tf]][, 1]
    
    # ref_sc_f  <- motifScore_ref$Forward[[cur_tf]][, 1]
    # ref_sc_r  <- motifScore_ref_rc$Forward[[cur_tf]][, 1]
    # alt_sc_f  <- motifScore_alt$Forward[[cur_tf]][, 1]
    # alt_sc_r  <- motifScore_alt_rc$Forward[[cur_tf]][, 1]
    
    
    if(sum(c(ref_sc_f, ref_sc_r, alt_sc_f, alt_sc_r) >= min_LLR[cur_tf], na.rm =T) == 0){
      print(paste0("There are no sites for ",  names(my_motifs)[cur_tf]))
      next()
    }else{ # if there is a site on either strand of ref or alt sequence
      print(paste("There is a site somewhere for", names(my_motifs)[cur_tf]))
      for(cur_strand in c(1:2)){
        if(sum(c(motifScore_ref[[cur_strand]][[cur_tf]][, 1], motifScore_alt[[cur_strand]][[cur_tf]][, 1]) >= min_LLR[cur_tf]) > 0){
          # print(paste("There is a site for ",names(my_motifs)[cur_tf], " on ", names(motifScore_ref)[cur_strand]))
          my_old <- max(motifScore_ref[[cur_strand]][[cur_tf]][, 1])
          my_old_pos <- which.max(motifScore_ref[[cur_strand]][[cur_tf]][, 1])
          my_new <- max(motifScore_alt[[cur_strand]][[cur_tf]][, 1])
          my_new_pos <- which.max(motifScore_alt[[cur_strand]][[cur_tf]][, 1])
          my_percent <- (my_new - my_old)/my_old
                           
          # print(my_percent)
          nu_alt <- sum(motifScore_alt[[cur_strand]][[cur_tf]][, 1] >= min_LLR[cur_tf])
          # print("nu_alt")
          # print(nu_alt)
          nu_ref <- sum(motifScore_ref[[cur_strand]][[cur_tf]][, 1] >= min_LLR[cur_tf])
          # print("nu_ref")
          # print(nu_ref)
          my_loss_gain <- nu_alt - nu_ref
          
         # if(max_len - my_old_pos >= 0){
            cur_mot_pos_1 <- max_len - my_old_pos + 1
            cur_mot_pos_2 <- max_len - my_new_pos + 1
          # }else{
          #   cur_mot_pos <- my_old_pos - max_len + 1
          # }
          
          my_res <- c(names(my_motifs)[cur_tf],
                      my_old,
                      my_new,
                      my_percent,
                      names(motifScore_alt)[cur_strand],
                      my_loss_gain, 
                      cur_mot_pos_1,
                      cur_mot_pos_2,
                      snp_ref,
                      snp_alt)
          names(my_res) <- colnames(my_results)
          # print(my_res)
          my_results <- rbind(my_results,
                              my_res)
        }#end if there is any site on this strand
      }#end loop over strand
    }# end if there was a site
  }# end of loop over TFs
  my_resultsdf <- data.frame(my_results, stringsAsFactors = F)
  my_resultsdf$change <- as.numeric(my_resultsdf$change)
  my_resultsdf$old <- as.numeric(my_resultsdf$old)
  my_resultsdf$new <- as.numeric(my_resultsdf$new)
  my_resultsdf$loss_or_gain <- as.numeric(my_resultsdf$loss_or_gain)
  return(my_resultsdf)
}
############################################################################################################
# example
aamin_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
aamaxL_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
for(i in 1:length(TF.motifs.Expanded_LLR_to_pval_pseudo)){
  aamaxL_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][1,1] / 1000
  aa <- which(TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][,2] > 0.0001)[1]
  if(aa == 1){
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][aa, 1] / 1000
  }else{
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][(aa-1), 1] / 1000
  }
}
names(aamin_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
names(aamaxL_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)

aa_names <- c("ESR1_2","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","RARA", "RELA","RUNX1","SP1","TFAP2C","YBX1")
TF.motifs.Expanded_pseudo_exp12 <- TF.motifs.Expanded_pseudo[aa_names]
TF.motifs.Expanded_pseudo_exp12_t <- lapply(TF.motifs.Expanded_pseudo_exp12, t)

aa_pos_dif_filtered_all2 <- aa_pos_dif_filtered_all
aa_neg_dif_filtered_all2 <- aa_neg_dif_filtered_all
aa_pos_dif_filtered_all2 <- aa_pos_dif_filtered_all2[rowSums(aa_pos_dif_filtered_all2 >= 0.05, na.rm = T) > 0,]
aa_neg_dif_filtered_all2 <- aa_neg_dif_filtered_all2[rowSums(aa_neg_dif_filtered_all2 >= 0.05, na.rm = T) > 0,]
aa_all_names <- c(rownames(aa_pos_dif_filtered_all2 ), rownames(aa_neg_dif_filtered_all2))
aa_all_names <- unlist(lapply(strsplit(aa_all_names, split = "_"), "[[", 2))
aa_all_qtl <- ER_gtex_eqtl_sig[as.numeric(aa_all_names),]
aa_eqtl <-  strsplit(aa_all_qtl$variant_id, split= "_")
aa_eqtl_chr <- unlist(lapply(aa_eqtl, "[[", 1))
aa_eqtl_st <- as.numeric(unlist(lapply(aa_eqtl, "[[", 2)))
aa_eqtl_ref <- unlist(lapply(aa_eqtl, "[[", 3))
aa_eqtl_alt <- unlist(lapply(aa_eqtl, "[[", 4))

aa <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                       my_motifs = TF.motifs.Expanded_pseudo_exp12_t,
                       snp_chr = aa_eqtl_chr[i],
                       snp_start = aa_eqtl_st[i],
                       snp_ref = aa_eqtl_ref[i],
                       snp_alt = aa_eqtl_alt[i],
                       min_LLR = aamin_LLR_lowest[aa_names])
############################################################################################################
########################################################################################
# write a function to read the variant output results and compute percentile change
compute_percentile_change <- function(WT_var_from_hal, Var_from_hal_neg, Var_from_hal_pos){
  # WT_RData_var : is the name of the variable in the RData file for WT saved from hal
  # Var_from_hal_neg : is the name of the variable for neg entries in the RData file for variants saved from hal
  # Var_from_hal_pos : is the name of the variable for pos entries in the RData file for variants saved from hal
  
  aa_neg_dif_percentile_4 <- list()
  aa_pos_dif_percentile_4 <- list()
  
  
  aa_WT_GT <- rownames(WT_var_from_hal)
  aa_WT_GT <- unlist(lapply(strsplit(aa_WT_GT, "_"), "[[", 1))
  aa_WT_GT[aa_WT_GT  == "pos"] <- 1
  aa_WT_GT[aa_WT_GT  == "neg"] <- 0
  aa_WT_GT <- as.numeric(aa_WT_GT)
  names(aa_WT_GT) <- rownames(WT_var_from_hal)
  
  # computing percentile changes for 
  for(i in 1:length(Var_from_hal_neg)){
    print(paste0("neg_",i," out of ", length(Var_from_hal_neg)))
    aa_neg_dif_percentile_4[[i]] <- matrix(nrow = nrow(Var_from_hal_neg[[i]]), 
                                           ncol = ncol(Var_from_hal_neg[[i]]))
    colnames(aa_neg_dif_percentile_4[[i]]) <- colnames(Var_from_hal_neg[[i]])
    rownames(aa_neg_dif_percentile_4[[i]]) <- rownames(Var_from_hal_neg[[i]])
    for(j in 1:nrow(Var_from_hal_neg[[i]])){# computing percentile change for each variant 
      for(kk in 1:ncol(Var_from_hal_neg[[i]])){# computing percentile change for each model
        aa_all_before <- sum(WT_var_from_hal[,kk] <= WT_var_from_hal[names(Var_from_hal_neg)[i],kk], na.rm = T)/sum(!is.na(WT_var_from_hal[,kk]))
        aa_all_after <- sum(WT_var_from_hal[,kk] <= Var_from_hal_neg[[i]][j, kk], na.rm = T)/sum(!is.na(WT_var_from_hal[,kk]))
        aa_neg_dif_percentile_4[[i]][j,kk] <- aa_all_after - aa_all_before
      }
    }
  }
  names(aa_neg_dif_percentile_4) <- names(Var_from_hal_neg)
  aa_neg_dif_percentile_4_all <- do.call(what = rbind, aa_neg_dif_percentile_4)
  
  for(i in 1:length(Var_from_hal_pos)){
    print(paste0("pos_",i," out of ", length(Var_from_hal_pos)))
    aa_pos_dif_percentile_4[[i]] <- matrix(nrow = nrow(Var_from_hal_pos[[i]]),
                                           ncol = ncol(Var_from_hal_pos[[i]]))
    colnames(aa_pos_dif_percentile_4[[i]]) <- colnames(Var_from_hal_pos[[i]])
    rownames(aa_pos_dif_percentile_4[[i]]) <- rownames(Var_from_hal_pos[[i]])
    for(j in 1:nrow(Var_from_hal_pos[[i]])){
      for(kk in 1:ncol(Var_from_hal_pos[[i]])){
        aa_all_before <- sum(WT_var_from_hal[,kk] <= WT_var_from_hal[names(Var_from_hal_pos)[i],kk], na.rm = T)/sum(!is.na(WT_var_from_hal[,kk]))
        aa_all_after <- sum(WT_var_from_hal[,kk] <= Var_from_hal_pos[[i]][j, kk], na.rm = T)/sum(!is.na(WT_var_from_hal[,kk]))
        aa_pos_dif_percentile_4[[i]][j,kk] <- aa_all_after - aa_all_before
      }
    }
  }
  names(aa_pos_dif_percentile_4) <- names(Var_from_hal_pos)
  aa_pos_dif_percentile_4_all <- do.call(what = rbind, aa_pos_dif_percentile_4)
  aa_dif_percentile_4_all <- rbind(aa_pos_dif_percentile_4_all, aa_neg_dif_percentile_4_all)
  return(aa_dif_percentile_4_all)
}
########################################################################################
# example
load("~/Documents/Shayan/BioInf/EstrogenReceptor/Final_variant_Set/Merged_all_output_exp18_hal_WT.RData")
load("~/Documents/Shayan/BioInf/EstrogenReceptor/Final_variant_Set/Merged_all_output_exp18_hal.RData")
pos_mergedVar_output_hal
neg_mergedVar_output_hal
WT_models_merged_all_hal
aaMerged_percentile_change <- compute_percentile_change(WT_var_from_hal = WT_models_merged_all_hal,
                                                        Var_from_hal_neg = neg_mergedVar_output_hal, 
                                                        Var_from_hal_pos = pos_mergedVar_output_hal)
#############################################################################################################
# write a function to conduct prediction tests and comparison to random, given a boolean condition and an aggregation mode
compare_ranking_to_random <- function(#comp_mode,
  all_percentile,
  rank_among_index= c(1:nrow(all_percentile)),
  TP, 
  sampling_range=numeric(0),
  do_boxplot=T, 
  percent_in_sample=F){
  # comp_mode can be one of: per_model , mean_pres , median_preds , max_preds , vote_models
  # all_percentile is the matrix containing percentile changes for all variants
  # rank_among_index : are the indecis among which random sampling and ranking will be done: used to restrict the full set
  # TP : is a boolean with size equal to the number of rows on all_percentile, true entries mark the True Positives for the test being done
  # sampling_range : is a numeric vector containing the length of the samples being drawn (and the number of ranked elements).
  # percent_in_sample : boolean, if True: it will compute scores as: (#TP in our sample) / sample_size
  stopifnot(#is.character(comp_mode), 
    nrow(all_percentile) == length(TP),
    max(rank_among_index) <= nrow(all_percentile), 
    sum(duplicated(rank_among_index)) == 0,
    !is.unsorted(rank_among_index))
  
  cur_mat <-  all_percentile[rank_among_index,]
  cur_lab <- TP[rank_among_index]
  aa_model_index_Sorted_percentile <- matrix(nrow = nrow(cur_mat),
                                             ncol = ncol(cur_mat))
  
  
  for(i in 1:ncol(aa_model_index_Sorted_percentile)){
    aa_model_index_Sorted_percentile[, i] <- sort(abs(cur_mat[, i]),decreasing = T, index.return = T)$ix
  }
  aa_max <- apply(abs(cur_mat), MARGIN = 1, FUN = max)
  aa_mea <- apply(abs(cur_mat), MARGIN = 1, FUN = mean)
  aa_med <- apply(abs(cur_mat), MARGIN = 1, FUN = median)
  
  aa_model_sort_ind_max <- sort(aa_max, decreasing = T, index.return = T)$ix
  aa_model_sort_ind_mea <- sort(aa_mea, decreasing = T, index.return = T)$ix
  aa_model_sort_ind_med <- sort(aa_med, decreasing = T, index.return = T)$ix
  
  aa_percModels <- rowSums(abs(cur_mat) > 0.001)/ncol(cur_mat)
  aa_perc_names <- sort(aa_percModels,decreasing = T,index.return=T)$ix
  
  if(length(sampling_range) == 0){
    sampling_range <- seq(1, length(all_percentile), 100)  #c(seq(1,100, 5), seq(110, length(aa_clinvar_index), 50))
  }
  aa_rep_num <- ncol(cur_mat)
  
  
  aa_random_Score <- matrix(nrow = aa_rep_num, ncol = length(sampling_range))
  colnames(aa_random_Score) <- sampling_range
  
  aa_model_Score_mean <- numeric(length = length(sampling_range))
  aa_model_Score_median <- numeric(length = length(sampling_range))
  aa_model_Score_max <- numeric(length = length(sampling_range))
  aa_model_Score_perc <- numeric(length = length(sampling_range))
  
  aa_model_Score_percentile <- matrix(nrow = length(sampling_range), ncol = ncol(cur_mat))
  colnames(aa_model_Score_percentile) <- colnames(cur_mat)
  rownames(aa_model_Score_percentile) <- sampling_range
  
  
  aa_tot <- sum(cur_lab, na.rm=T)
  for(i in 1:length(sampling_range)){
    for(j in 1:aa_rep_num){
      aasampl <- sample(x = c(1:length(cur_lab)), size = sampling_range[i], replace = F)
      if(percent_in_sample){
        aa_random_Score[j, i] <- sum(cur_lab[aasampl], na.rm=T)/sampling_range[i]
      }else{
        aa_random_Score[j, i] <- sum(cur_lab[aasampl], na.rm=T)/aa_tot
      }
      
    }
    for(j in 1:ncol(cur_mat)){
      if(percent_in_sample){
        aa_model_Score_percentile[i, j] <- sum(cur_lab[aa_model_index_Sorted_percentile[1:sampling_range[i],j]], na.rm=T)/sampling_range[i]
      }else{
        aa_model_Score_percentile[i, j] <- sum(cur_lab[aa_model_index_Sorted_percentile[1:sampling_range[i],j]], na.rm=T)/aa_tot
      }
      
    }
    if(percent_in_sample){
      aa_model_Score_mean[i] <- sum(cur_lab[aa_model_sort_ind_mea[1:sampling_range[i]]], na.rm=T)/sampling_range[i]
      aa_model_Score_median[i] <- sum(cur_lab[aa_model_sort_ind_med[1:sampling_range[i]]] , na.rm=T)/sampling_range[i]
      aa_model_Score_max[i] <- sum(cur_lab[aa_model_sort_ind_max[1:sampling_range[i]]], na.rm=T)/sampling_range[i]
      aa_model_Score_perc[i] <-sum(cur_lab[aa_perc_names[1:sampling_range[i]]], na.rm=T)/sampling_range[i]
    }else{
      aa_model_Score_mean[i] <- sum(cur_lab[aa_model_sort_ind_mea[1:sampling_range[i]]], na.rm=T)/aa_tot
      aa_model_Score_median[i] <- sum(cur_lab[aa_model_sort_ind_med[1:sampling_range[i]]] , na.rm=T)/aa_tot
      aa_model_Score_max[i] <- sum(cur_lab[aa_model_sort_ind_max[1:sampling_range[i]]], na.rm=T)/aa_tot
      aa_model_Score_perc[i] <-sum(cur_lab[aa_perc_names[1:sampling_range[i]]], na.rm=T)/aa_tot
      
    }
  }
  
  if(do_boxplot){
    boxplot.matrix(aa_random_Score, las = 2, xlab = "#predictions",
                   ylab = "TPR")
    points(aa_model_Score_mean, col = 2, pch = 16, cex = 0.7)
    points(aa_model_Score_median, col = 3, pch = 17, cex = 0.7)
    points(aa_model_Score_max, col = 4, pch = 18, cex = 0.7)
    points(aa_model_Score_perc, col = 5, pch = 19, cex = 0.7)
    legend(x = "topleft",legend = c("mean", "median", "max", "vote"), 
           pch = c(16,17,18,19), fill = c(2,3,4,5), cex = 0.5)
    
  }
  aa_allmat <- matrix(nrow = ncol(cur_mat), ncol = ncol(aa_random_Score) * 2)
  aac <- 1
  for(i in seq(1, ncol(aa_allmat), 2)){
    aa_allmat[, i] <- aa_random_Score[, aac]
    aa_allmat[, i+1] <- t(aa_model_Score_percentile)[, aac]
    aac <- aac + 1
  }
  if(do_boxplot){
    boxplot.matrix(aa_allmat, col = rep(c(2,3), ncol(cur_mat)),
                   xaxt = "n", xlab = "#predictions", ylab = "TPR")
    axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
         labels = sampling_range, las = 2)
    abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)
  }
  
  return(list(Random_scores = aa_random_Score,
              Model_scores = aa_model_Score_percentile,
              mean_ranked_score = aa_model_Score_mean, 
              median_ranked_score = aa_model_Score_median, 
              max_ranked_score=aa_model_Score_max, 
              percent_ranked_score= aa_model_Score_perc))
  
}
########################################################################################
# example


aalab_orig <- ER_SNPs_merged_allbypos_modified$is_gtex
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
aaTP <- aalab_orig[aasap]

aa_ex <- compare_ranking_to_random(all_percentile = Merged_percentile_change_noDUP,
                                   rank_among_index= c(1:nrow(Merged_percentile_change_noDUP)),
                                   TP = aaTP, 
                                   sampling_range=c(seq(10,300,10), seq(400,2000, 100), seq(3000,33000,1000)),
                                   do_boxplot=T)



################################################################################################################################################################################
################################################################################################################################################################################
ER_gtex_eqtl_sig <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/Breast_Mammary_Tissue.v8.signif_variant_gene_pairs.txt", 
                               stringsAsFactors = F, header = T)
aa_eqtl <-  strsplit(ER_gtex_eqtl_sig$variant_id, split= "_")
aa_eqtl_chr <- unlist(lapply(aa_eqtl, "[[", 1))
aa_eqtl_st <- unlist(lapply(aa_eqtl, "[[", 2))
aa_eqtl_ref <- unlist(lapply(aa_eqtl, "[[", 3))
aa_eqtl_alt <- unlist(lapply(aa_eqtl, "[[", 4))

aa_ch <- nchar(aa_eqtl_ref)
which(aa_ch > 1)[1:5]
aa_eqtl[which(aa_ch > 1)[1:5]]

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_eqtl_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/sig_eqtl.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_qtl <- readPeakFile(aa_eqtl_Add, as = "GRanges")
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_15")
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)
aa_pos_neg_vars <- write_variant_seq(enhancer_GR = aa_pos_neg, 
                                 eqtl_GR = aa_qtl,
                                 eqtl_ref_alt = cbind(aa_eqtl_ref, aa_eqtl_alt),
                                 all_combs = F,
                                 my_genome = BSgenome.Hsapiens.UCSC.hg38, 
                                 label = c(rep(1, length(aa_pos)), rep(0, length(aa_neg))))

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")

# write GEMSTAT jobs to run
# write hal job for KDs
aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas1[1:25]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas2[1:25]]
aaadd3 <- union(aaadd1, aaadd2)

aa_par_name <- paste0(aaadd3, ".txt.Filter")
aa_seq_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_15/Variant_seq/",
                                   recursive = T)
aa_lab_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_15/Variant_Labels/",
                           recursive = T)
aa_nam_spl <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_seq_files, split = "\\."), "[[", 1)), split = "_"), "[", c(3,4)), paste, collapse = "_"))


setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_15")
for(i in 1:length(aa_par_name)){
  for(j in 1:length(aa_seq_files)){
    cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr ",
          paste0("-s Variant_seq/", aa_seq_files[j]), 
          paste0("-e Variant_Labels/", aa_lab_files[j]),
          "-m motifs.wtmx -f TF_exp.tab", 
          paste0("-fo Varaint_out/", aa_nam_spl[j], "_", aaadd3[i], ".out"), 
          "-o DIRECT -c Coop/coop.par ", 
          paste0("-p Trained_par/", aa_par_name[i]),
          "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1\n"
          ), 
        sep = " ", append = (i!=1), file = "variant_jobs.job")
  }
}

# read output and find the ones with some effect
# first read WT for all models

aa_wt_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out_eqtl/",
                              pattern = "NA_NA_par*", full.names = T)
aa_wt_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out_eqtl/",
                              pattern = "NA_NA_par*", full.names = F)
aa_nam_spl <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_wt_files, split = "\\."), "[[", 1)), split = "_"), "[", c(3,4)), paste, collapse = "_"))

aa_t <- read_output_train_test_GEMSTAT_indiv(output_file = aa_wt_files_ful[1])
aa_wt_models <- matrix(nrow=nrow(aa_t$GT),
                       ncol=length(aa_wt_files))

colnames(aa_wt_models) <- aa_nam_spl
rownames(aa_wt_models) <- aa_t$pred$Rows
for(i in 1:ncol(aa_wt_models)){
  aa_wt_models[, i] <- read_output_train_test_GEMSTAT_indiv(output_file = aa_wt_files_ful[i])$pred$X1
}
aa_neg_out <- list()
aa_pos_out <- list()

aa_neg_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out_eqtl//",
                              pattern = "neg_*", full.names = T)
aa_neg_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out_eqtl//",
                          pattern = "neg_*", full.names = F)
aa_neg_sp_1 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_neg_files, split = "\\."), "[[", 1)), split = "_"), "[", c(1,2)), paste, collapse = "_"))
aa_neg_sp_2 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_neg_files, split = "\\."), "[[", 1)), split = "_"), "[", c(3,4)), paste, collapse = "_"))
aa_neg_uniq <- unique(aa_neg_sp_1)

aa_pos_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out_eqtl/",
                               pattern = "pos_*", full.names = T)
aa_pos_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out_eqtl/",
                           pattern = "pos_*", full.names = F)
aa_pos_sp_1 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_pos_files, split = "\\."), "[[", 1)), split = "_"), "[", c(1,2)), paste, collapse = "_"))
aa_pos_sp_2 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_pos_files, split = "\\."), "[[", 1)), split = "_"), "[", c(3,4)), paste, collapse = "_"))
aa_pos_uniq <- unique(aa_pos_sp_1)

aa_neg_dif <- list()
aa_pos_dif <- list()

for(i in 1:length(aa_neg_uniq)){
  aa_cur_ind <- which(aa_neg_sp_1 %in% aa_neg_uniq[i])
  aa_t <- read.table(file = aa_neg_files_ful[aa_cur_ind[1]], header = T, stringsAsFactors = F)
  aa_ind <- seq(2, nrow(aa_t), 2)
  aa_t <- aa_t[aa_ind, ]
  aa_neg_out[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models))
  colnames(aa_neg_out[[i]]) <- colnames(aa_wt_models)
  rownames(aa_neg_out[[i]]) <- aa_t$Rows
  aa_neg_dif[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models))
  colnames(aa_neg_dif[[i]]) <- colnames(aa_wt_models)
  rownames(aa_neg_dif[[i]]) <- aa_t$Rows
  for(j in 1:length(aa_cur_ind)){
    aa_table <- read.table(file =aa_neg_files_ful[aa_cur_ind[j]], 
                           header = T, stringsAsFactors = F)
    aa_neg_out[[i]][, aa_neg_sp_2[aa_cur_ind[j]]] <- aa_table[seq(2, nrow(aa_table), 2),"X1"]
  }
  for(j in 1:nrow(aa_neg_out[[i]])){
    aa_neg_dif[[i]][j, ] <- (aa_neg_out[[i]][j, ] - aa_wt_models[aa_neg_uniq[i], ])/ aa_wt_models[aa_neg_uniq[i], ]
    
  }
}
names(aa_neg_out) <- aa_neg_uniq
names(aa_neg_dif) <- aa_neg_uniq

for(i in 1:length(aa_pos_uniq)){
  aa_cur_ind <- which(aa_pos_sp_1 %in% aa_pos_uniq[i])
  aa_t <- read.table(file = aa_pos_files_ful[aa_cur_ind[1]], header = T, stringsAsFactors = F)
  aa_ind <- seq(2, nrow(aa_t), 2)
  aa_t <- aa_t[aa_ind, ]
  aa_pos_out[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models))  
  colnames(aa_pos_out[[i]]) <- colnames(aa_wt_models)
  rownames(aa_pos_out[[i]]) <- aa_t$Rows
  aa_pos_dif[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models))
  colnames(aa_pos_dif[[i]]) <- colnames(aa_wt_models)
  rownames(aa_pos_dif[[i]]) <- aa_t$Rows
  for(j in 1:length(aa_cur_ind)){
    aa_table <- read.table(file =aa_pos_files_ful[aa_cur_ind[j]], 
                           header = T, stringsAsFactors = F)
    aa_pos_out[[i]][, aa_pos_sp_2[aa_cur_ind[j]]] <- aa_table[seq(2, nrow(aa_table), 2),"X1"]
  }
  for(j in 1:nrow(aa_pos_out[[i]])){
    aa_pos_dif[[i]][j, ] <- (aa_pos_out[[i]][j, ] - aa_wt_models[aa_pos_uniq[i], ])/ aa_wt_models[aa_pos_uniq[i], ]
  }
  
}
names(aa_pos_out) <- aa_pos_uniq
names(aa_pos_dif) <- aa_pos_uniq
aa_pos_dif$pos_1[2, 20:40]
(aa_pos_out$pos_1[2,20:40] - aa_wt_models["pos_1", 20:40])/(aa_wt_models["pos_1", 20:40])
max(aa_pos_dif$pos_1, na.rm = T)

aa_pos_dif_all <- do.call(what = rbind, aa_pos_dif)
boxplot.matrix(t(aa_pos_dif_all)*100)

aa_neg_dif_all <- do.call(what = rbind, aa_neg_dif)
boxplot.matrix(t(aa_neg_dif_all))

z <- zClust(x=t(aa_pos_dif_all[,-c(17)]), scale="row", zlim=c(-3,3), method="average")
# require(RColorBrewer)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(9)
heatmap.2(z$data, trace='none',
          col=rev(cols),
          #Rowv=F,Colv=F,
          main="predicted change"
          ,Rowv=z$Rowv, 
          Colv=z$Colv
          )
# cols <- colorRampPalette(brewer.pal(10, "RdBu"))(9)
# heatmap.2(x = t(aa_pos_dif_all[,-c(17)]), 
#           trace='none',
#           Rowv = T,
#           Colv = T,
#           breaks = c(-1,-0.5,-0.1, 0.1,0.5, 1, 1.5, 2, 2.5, 3), col = rev(cols))
z <- zClust(x=t(aa_neg_dif_all[,-c(17)]), scale="row", zlim=c(-3,3), method="average")
# require(RColorBrewer)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(9)
heatmap.2(z$data, trace='none',
          col=rev(cols),
          Rowv=F,Colv=F,
          main="predicted change neg"
          # ,Rowv=z$Rowv, 
          # Colv=z$Colv
)
# find enhancers containing multiple eqtls where some don't affect expression and some do


aa_pos_dif_filtered <- aa_pos_dif
for(i in length(aa_pos_dif_filtered):1){
  if(max(abs(aa_pos_dif_filtered[[i]]), na.rm = T) < 0.05){
    aa_pos_dif_filtered[[i]] <- NULL
  }
}
aanrwo <- unlist(lapply(aa_pos_dif_filtered, nrow))
par(mar = c(2,3,2,2), mfrow = c(5, 3))
for(i in 1:length(aa_pos_dif_filtered)){
  boxplot.matrix(t(aa_pos_dif_filtered[[i]]), 
                 las= 2, 
                 main = names(aa_pos_dif_filtered)[i],
                 xaxt = "n")
  abline(h= seq(-3, 3, 0.05), col = 2, lwd = 0.5, lty = 3)
  abline(h= 0, col = 2, lwd = 1, lty = 1)
}
# boxplot.matrix(t(aa_pos_dif_filtered$pos_312), 
#                las= 2, main = "pos_312 predicted change")
# View(aa_pos_dif_filtered$pos_171)

aa_pos_dif_filtered_all <- do.call(rbind, aa_pos_dif_filtered)
aa_neg_dif_filtered <- aa_neg_dif
for(i in length(aa_neg_dif_filtered):1){
  if(max(abs(aa_neg_dif_filtered[[i]]), na.rm = T) < 0.05){
    aa_neg_dif_filtered[[i]] <- NULL
  }
}
aa_neg_dif_filtered_all <-  do.call(rbind, aa_neg_dif_filtered)
aanrwo <- unlist(lapply(aa_neg_dif_filtered, nrow))
par(mar = c(1,1,1,1), mfrow = c(10, 7))
for(i in 1:length(aa_neg_dif_filtered)){
  if(aanrwo[i] > 1){
    boxplot.matrix(t(aa_neg_dif_filtered[[i]]), 
                   las= 2, 
                   main = names(aa_neg_dif_filtered)[i],
                   xaxt = "n")
    abline(h= seq(-3, 3, 0.05), col = 2, lwd = 0.5, lty = 3)
    abline(h= 0, col = 2, lwd = 1, lty = 1)
  }

}


aa_pos_dif_filtered$pos_171
ER_gtex_eqtl_sig[1372231,]

# get a set of eqtls that are make least 5% difference in some model
aa_pos_dif_filtered_all2 <- aa_pos_dif_filtered_all
aa_neg_dif_filtered_all2 <- aa_neg_dif_filtered_all
aa_pos_dif_filtered_all2 <- aa_pos_dif_filtered_all2[rowSums(aa_pos_dif_filtered_all2 >= 0.05, na.rm = T) > 0,]
aa_neg_dif_filtered_all2 <- aa_neg_dif_filtered_all2[rowSums(aa_neg_dif_filtered_all2 >= 0.05, na.rm = T) > 0,]

####################################################
# only using predictions that are somewhat accurate (positives greater than median, negatives less than median)
aa_neg_dif_acc <- aa_neg_dif 
aa_pos_dif_acc <- aa_pos_dif 

for(i in 1:length(aa_neg_uniq)){
  for(j in 1:nrow(aa_neg_out[[i]])){
    aa_neg_dif_acc[[i]][j, ] <- (aa_neg_out[[i]][j, ] - Exp12_WT_model_prediction_Filtered[aa_neg_uniq[i], ])/ Exp12_WT_model_prediction_Filtered[aa_neg_uniq[i], ]
  }
}
names(aa_neg_dif_acc) <- aa_neg_uniq

for(i in 1:length(aa_pos_uniq)){
  for(j in 1:nrow(aa_pos_out[[i]])){
    aa_pos_dif_acc[[i]][j, ] <- (aa_pos_out[[i]][j, ] - Exp12_WT_model_prediction_Filtered[aa_pos_uniq[i], ])/ Exp12_WT_model_prediction_Filtered[aa_pos_uniq[i], ]
  }
}
names(aa_pos_dif_acc) <- aa_pos_uniq


aa_neg_dif_filtered_acc <- aa_neg_dif_acc
for(i in length(aa_neg_dif_filtered_acc):1){
  if(max(abs(aa_neg_dif_filtered_acc[[i]]), na.rm = T) < 0.05){
    aa_neg_dif_filtered_acc[[i]] <- NULL
  }
}

aa_pos_dif_filtered_acc <- aa_pos_dif_acc
for(i in length(aa_pos_dif_filtered_acc):1){
  if(max(abs(aa_pos_dif_filtered_acc[[i]]), na.rm = T) < 0.05){
    aa_pos_dif_filtered_acc[[i]] <- NULL
  }
}

aa_pos_dif_filtered_all_acc <- do.call(rbind, aa_pos_dif_filtered_acc)
aa_neg_dif_filtered_all_acc <- do.call(rbind, aa_neg_dif_filtered_acc)
nrow(aa_pos_dif_filtered_all_acc)
nrow(aa_neg_dif_filtered_all_acc)

#########################################################################################################
# remove enhancers with less confident predictions from the results
aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas1[1:25]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas2[1:25]]
aaadd3 <- union(aaadd1, aaadd2)

par(mfrow = c(6,7), mar = c(2,2,2,2))
for(i in 1:length(aaadd3)){
  boxplot(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_results[[paste0("output_", aaadd3[i])]]$pred$X1~
            E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_results[[paste0("output_", aaadd3[i])]]$GT$X1, 
          main = aaadd3[i])
}
par(mfrow = c(6,7), mar = c(2,2,2,2))
for(i in 1:length(aaadd3)){
  boxplot(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$pred$X1~
            E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$GT$X1, 
          main = aaadd3[i])
}
par(mfrow = c(6,7), mar = c(2,2,2,2))
for(i in 1:length(aaadd3)){
  boxplot(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$pred$X1~
            E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$GT$X1, 
          main = aaadd3[i])
}

png("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/predicted_exp_all.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

par(mfrow = c(6,7), mar = c(2,2,2,2))
for(i in 1:length(aaadd3)){
  boxplot(c(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_results[[paste0("output_", aaadd3[i])]]$pred$X1, 
            E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$pred$X1,
            E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$pred$X1)~
            c(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_results[[paste0("output_", aaadd3[i])]]$GT$X1,
              E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$GT$X1,
              E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$GT$X1)
            , 
          main = aaadd3[i])
}
dev.off()

aa_all_pred <- matrix(nrow = 2118, ncol = 42)
colnames(aa_all_pred) <- aaadd3
rownames(aa_all_pred) <- c(rownames(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_results$output_par_1_1$pred),
                           rownames(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_results$output_par_1_1_test$pred),
                           rownames(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_results$output_par_1_1_valid$pred))

for(i in 1:length(aaadd3)){
  aa_all_pred[, i] <- c(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_results[[paste0("output_", aaadd3[i])]]$pred$X1, 
                        E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$pred$X1,
                        E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$pred$X1)
}

aa_all_gt <- c(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_results[[paste0("output_", aaadd3[i])]]$GT$X1,
               E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$GT$X1,
               E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$GT$X1)
aa_all_pred_filtered <- aa_all_pred
# in each column 
for(i in 1:ncol(aa_all_pred_filtered)){
  aa_quantile <- quantile(aa_all_pred_filtered[, i])
  aa_all_pred_filtered[(aa_all_pred_filtered[,i] > aa_quantile[3] & aa_all_gt == 0), i] <- NA
  aa_all_pred_filtered[(aa_all_pred_filtered[,i] < aa_quantile[3] & aa_all_gt == 1), i] <- NA
}
barplot(sort(colSums(is.na(aa_all_pred_filtered))))
barplot(sort(rowSums(is.na(aa_all_pred_filtered)), decreasing = T))
sum(rowSums(is.na(aa_all_pred_filtered)) == 42 & aa_all_gt == 0)



png("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/predicted_exp_all_filtered_median.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

par(mfrow = c(6,7), mar = c(2,2,2,2))
for(i in 1:length(aaadd3)){
  boxplot(aa_all_pred_filtered[, i]~aa_all_gt
          , 
          main = aaadd3[i])
}
dev.off()
# filter the eqtl results 
GEMSTAT_Ensemble_train_SeqList[[4]]
GEMSTAT_Ensemble_test_SeqList[[2]]
GEMSTAT_Ensemble_validation_SeqList[[1]]


aa_train_match_pos <- match(GEMSTAT_Ensemble_train_SeqList[[4]], 
                            Positive_set_seq_list_char_1000bp[[4]])
aa_train_match_neg <- match(GEMSTAT_Ensemble_train_SeqList[[4]], 
                            Negative_set_seq_list_char_1000bp[[4]])
aa_train_match_all <- character(length = length(aa_train_match_pos))
for(i in 1:length(aa_train_match_pos)){
  if(!is.na(aa_train_match_pos[i])){
    aa_train_match_all[i] <- paste0("pos_", aa_train_match_pos[i])
  }else if(!is.na(aa_train_match_neg[i])){
    aa_train_match_all[i] <- paste0("neg_", aa_train_match_neg[i])
  }
}

aa_test_match_pos <- match(GEMSTAT_Ensemble_test_SeqList[[2]], 
                           Positive_set_seq_list_char_1000bp[[4]])
aa_test_match_neg <- match(GEMSTAT_Ensemble_test_SeqList[[2]],
                           Negative_set_seq_list_char_1000bp[[4]])
aa_test_match_all <- character(length = length(aa_test_match_pos))
for(i in 1:length(aa_test_match_pos)){
  if(!is.na(aa_test_match_pos[i])){
    aa_test_match_all[i] <- paste0("pos_", aa_test_match_pos[i])
  }
  if(!is.na(aa_test_match_neg[i])){
    aa_test_match_all[i] <- paste0("neg_", aa_test_match_neg[i])
  }
}
aa_valid_match_pos <- match(GEMSTAT_Ensemble_validation_SeqList[[1]],
                            Positive_set_seq_list_char_1000bp[[4]])
aa_valid_match_neg <- match(GEMSTAT_Ensemble_validation_SeqList[[1]],
                            Negative_set_seq_list_char_1000bp[[4]])
aa_valid_match_all <- character(length = length(aa_valid_match_pos))
for(i in 1:length(aa_valid_match_pos)){
  if(!is.na(aa_valid_match_pos[i])){
    aa_valid_match_all[i] <- paste0("pos_", aa_valid_match_pos[i])
  }else if(!is.na(aa_valid_match_neg[i])){
    aa_valid_match_all[i] <- paste0("neg_", aa_valid_match_neg[i])
  }
}

aa_train_test_valid_nam <- c(aa_train_match_all,
                             aa_test_match_all,
                             aa_valid_match_all)

aa_old_names <-c(names(GEMSTAT_Ensemble_train_SeqList[[4]]), 
                 GEMSTAT_Ensemble_test_SeqList[[2]],
                 GEMSTAT_Ensemble_validation_SeqList[[1]])
# naming WT model predictions

aa_all_pred <- matrix(nrow = 2118, ncol = 42)
colnames(aa_all_pred) <- aaadd3
rownames(aa_all_pred) <- aa_train_test_valid_nam

for(i in 1:length(aaadd3)){
  aa_all_pred[, i] <- c(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_results[[paste0("output_", aaadd3[i])]]$pred$X1, 
                        E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$pred$X1,
                        E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$pred$X1)
}

aa_all_gt <- c(E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_results[[paste0("output_", aaadd3[i])]]$GT$X1,
               E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$GT$X1,
               E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$GT$X1)
names(aa_all_gt) <- aa_train_test_valid_nam

Positive_set_gr_list_1000bp[[4]]
Negative_set_gr_list_1000bp[[4]]

Exp12_WT_model_prediction <- aa_all_pred
aasr <- sort(colnames(Exp12_WT_model_prediction), index.return = T)$ix
Exp12_WT_model_prediction <- Exp12_WT_model_prediction[, aasr]
Exp12_WT_model_prediction_Filtered <- Exp12_WT_model_prediction
for(i in 1:ncol(Exp12_WT_model_prediction_Filtered)){
  aa_quantile <- quantile(Exp12_WT_model_prediction_Filtered[, i])
  Exp12_WT_model_prediction_Filtered[(Exp12_WT_model_prediction_Filtered[,i] > aa_quantile[3] & aa_all_gt == 0), i] <- NA
  Exp12_WT_model_prediction_Filtered[(Exp12_WT_model_prediction_Filtered[,i] < aa_quantile[3] & aa_all_gt == 1), i] <- NA
  
}

# filtering results for experiment 16 results

aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[16]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[16]]$annot)[aas1[1:25]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[16]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[16]]$annot)[aas2[1:25]]
aaadd3 <- union(aaadd1, aaadd2)


aa_all_pred <- matrix(nrow = 2118, ncol = length(aaadd3))
colnames(aa_all_pred) <- aaadd3
rownames(aa_all_pred) <- aa_train_test_valid_nam

for(i in 1:length(aaadd3)){
  aa_all_pred[, i] <- c(E_RNA_GEMSTAT_Ensemble_Outlist[[16]]$train_results[[paste0("output_", aaadd3[i])]]$pred$X1, 
                        E_RNA_GEMSTAT_Ensemble_Outlist[[16]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$pred$X1,
                        E_RNA_GEMSTAT_Ensemble_Outlist[[16]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$pred$X1)
}

aa_all_gt <- c(E_RNA_GEMSTAT_Ensemble_Outlist[[16]]$train_results[[paste0("output_", aaadd3[i])]]$GT$X1,
               E_RNA_GEMSTAT_Ensemble_Outlist[[16]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$GT$X1,
               E_RNA_GEMSTAT_Ensemble_Outlist[[16]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$GT$X1)
names(aa_all_gt) <- aa_train_test_valid_nam

#Positive_set_gr_list_1000bp[[4]]
#Negative_set_gr_list_1000bp[[4]]

Exp16_WT_model_prediction <- aa_all_pred
aasr <- sort(colnames(Exp16_WT_model_prediction), index.return = T)$ix
Exp16_WT_model_prediction <- Exp16_WT_model_prediction[, aasr]
Exp16_WT_model_prediction_Filtered <- Exp16_WT_model_prediction
for(i in 1:ncol(Exp16_WT_model_prediction_Filtered)){
  aa_quantile <- quantile(Exp16_WT_model_prediction_Filtered[, i])
  Exp16_WT_model_prediction_Filtered[(Exp16_WT_model_prediction_Filtered[,i] > aa_quantile[3] & aa_all_gt == 0), i] <- NA
  Exp16_WT_model_prediction_Filtered[(Exp16_WT_model_prediction_Filtered[,i] < aa_quantile[3] & aa_all_gt == 1), i] <- NA
  
}
par(mfrow = c(7,7), mar = c(1,1,1,1))
for(i in 1:ncol(Exp16_WT_model_prediction)){
  boxplot(Exp16_WT_model_prediction[, i]~aa_all_gt, las = 2)
}
par(mfrow = c(7,7), mar = c(1,1,1,1))
for(i in 1:ncol(Exp16_WT_model_prediction)){
  boxplot(Exp16_WT_model_prediction_Filtered[, i]~aa_all_gt, las = 2)
}
colSums(is.na(Exp16_WT_model_prediction_Filtered))
sum(rowSums(is.na(Exp16_WT_model_prediction_Filtered))  >= 20)
#########################################################################################################
# filtering results for experiment 18 results

aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas1[1:150]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas2[1:150]]
aaadd3 <- union(aaadd1, aaadd2)
aaadd33 <- unlist(lapply(lapply(strsplit(aaadd3, "_"), "[", c(2,3,4)), paste, collapse="_"))

aa_all_pred <- matrix(nrow = 2118, ncol = length(aaadd33))
colnames(aa_all_pred) <- aaadd33
rownames(aa_all_pred) <- aa_train_test_valid_nam

for(i in 1:length(aaadd3)){
  aa_all_pred[, i] <- c(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$train_results[[paste0("output_", aaadd33[i])]]$pred$X1, 
                        E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$pred$X1,
                        E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$pred$X1)
}

aa_all_gt <- c(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$train_results[[paste0("output_", aaadd33[i])]]$GT$X1,
               E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$test_results[[paste0("output_", aaadd3[i], "_test")]]$GT$X1,
               E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_results[[paste0("output_", aaadd3[i], "_valid")]]$GT$X1)
names(aa_all_gt) <- aa_train_test_valid_nam

#Positive_set_gr_list_1000bp[[4]]
#Negative_set_gr_list_1000bp[[4]]

Exp18_WT_model_prediction <- aa_all_pred
aasr <- sort(colnames(Exp18_WT_model_prediction), index.return = T)$ix
Exp18_WT_model_prediction <- Exp18_WT_model_prediction[, aasr]
Exp18_WT_model_prediction_Filtered <- Exp18_WT_model_prediction
for(i in 1:ncol(Exp18_WT_model_prediction_Filtered)){
  aa_quantile <- quantile(Exp18_WT_model_prediction_Filtered[, i])
  Exp18_WT_model_prediction_Filtered[(Exp18_WT_model_prediction_Filtered[,i] > aa_quantile[3] & aa_all_gt == 0), i] <- NA
  Exp18_WT_model_prediction_Filtered[(Exp18_WT_model_prediction_Filtered[,i] < aa_quantile[3] & aa_all_gt == 1), i] <- NA
  
}
par(mfrow = c(16,16), mar = c(0.1,0.1,0.1,0.1))
for(i in 1:ncol(Exp18_WT_model_prediction)){
  boxplot(Exp18_WT_model_prediction[, i]~aa_all_gt, las = 2, xlab = "", ylab="", yaxt = "n")
}
par(mfrow = c(7,7),  mar = c(0.1,0.1,0.1,0.1))
for(i in 1:ncol(Exp18_WT_model_prediction)){
  boxplot(Exp18_WT_model_prediction_Filtered[, i]~aa_all_gt, las = 2, xlab = "", ylab="", yaxt = "n")
}
colSums(is.na(Exp18_WT_model_prediction_Filtered))
sum(rowSums(is.na(Exp18_WT_model_prediction_Filtered))  >= 200)
#########################################################################################################
# looking at parameters of chosen models

aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas1[1:150]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas2[1:150]]
aaadd3 <- union(aaadd1, aaadd2)
aaadd33 <- unlist(lapply(lapply(strsplit(aaadd3, "_"), "[", c(2,3,4)), paste, collapse="_"))


rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[1:10]
aa_annot_All <- E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot[paste0("trained_", aaadd33),]
aa_bind_All <- E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$binding[paste0("trained_", aaadd33),]
aa_alpha_All <- E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$alpha[paste0("trained_", aaadd33),]
aa_coop_All <- E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$coop[paste0("trained_", aaadd33),]
aa_qbtm_All <- E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$qbtm[paste0("trained_", aaadd33)]
aa_beta_All <- E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$beta[paste0("trained_", aaadd33)]
aa_logi_All <- E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$logistic[paste0("trained_", aaadd33),]

heatmap.2(cbind(log10(aa_bind_All),
                log10(aa_alpha_All),
                log10(aa_coop_All),
                log10(aa_qbtm_All),
                log10(aa_beta_All),
                log10(aa_logi_All)),
          Rowv = T, Colv = F, trace = "none")

aasrt <- sort(aa_all_sig_hold_perc_2_uniq$`0.05`[6,], index.return=T)$ix

hist(aa_pos_neg_out_snp_all_2[,aasrt[1]])
par(mfrow = c(2,2))
#hist(Exp18_WT_model_prediction_Filtered[,aasrt[1] ])

hist(Exp18_WT_model_prediction[unlist(lapply(strsplit(rownames(Exp18_WT_model_prediction), "_"), "[[", 1)) == "neg",
                               aasrt[1] ], main = "UnFiltered neg", xlim = range(Exp18_WT_model_prediction[,aasrt[1]]), breaks = 20)
hist(Exp18_WT_model_prediction[unlist(lapply(strsplit(rownames(Exp18_WT_model_prediction), "_"), "[[", 1)) == "pos",
                               aasrt[1] ], main = "UnFiltered pos", xlim = range(Exp18_WT_model_prediction[,aasrt[1]]), breaks = 20)

hist(Exp18_WT_model_prediction_Filtered[unlist(lapply(strsplit(rownames(Exp18_WT_model_prediction), "_"), "[[", 1)) == "neg",
                               aasrt[1] ], main = "Filtered neg", xlim = range(Exp18_WT_model_prediction[,aasrt[1]]), breaks = 20)
hist(Exp18_WT_model_prediction_Filtered[unlist(lapply(strsplit(rownames(Exp18_WT_model_prediction), "_"), "[[", 1)) == "pos",
                               aasrt[1] ], main = "Filtered pos", xlim = range(Exp18_WT_model_prediction[,aasrt[1]]), breaks = 20)

hist(aa_pos_neg_out_snp_all_2[,aasrt[1]])
hist(aa_all_dif_snp_perc_2[,aasrt[1]], breaks = 200)
hist(aa_all_dif_snp_2[,aasrt[1]])
sum(aa_all_dif_snp_2 > 1 , na.rm = T)
sum(aa_all_dif_snp > 1 , na.rm = T)


#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
# reading in eqtls

ER_gtex_eqtl_sig <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/Breast_Mammary_Tissue.v8.signif_variant_gene_pairs.txt", stringsAsFactors = F, header = T)
ER_gtex_eqtl_sig$variant_id[1:10]

aa_eqtl <-  strsplit(ER_gtex_eqtl_sig$variant_id, split= "_")
aa_eqtl_chr <- unlist(lapply(aa_eqtl, "[[", 1))
aa_eqtl_st <- unlist(lapply(aa_eqtl, "[[", 2))
aa_eqtl_ref <- unlist(lapply(aa_eqtl, "[[", 3))
aa_eqtl_alt <- unlist(lapply(aa_eqtl, "[[", 4))

pie(table(unlist(strsplit(aa_eqtl_alt, ""))))
pie(table(unlist(strsplit(aa_eqtl_ref, ""))))


aa_eqtl_df <- data.frame(chr = aa_eqtl_chr,
                         start = aa_eqtl_st, 
                         end = aa_eqtl_st,
                         refrence = aa_eqtl_ref,
                         alt = aa_eqtl_alt,
                         gene = ER_gtex_eqtl_sig$gene_id, 
                         tss_distance = ER_gtex_eqtl_sig$tss_distance,
                         ma_samples = ER_gtex_eqtl_sig$ma_samples,
                         ma_count = ER_gtex_eqtl_sig$ma_count,
                         maf= ER_gtex_eqtl_sig$maf, 
                         pval_nominal = ER_gtex_eqtl_sig$pval_nominal,
                         slope = ER_gtex_eqtl_sig$slope,
                         slope_se = ER_gtex_eqtl_sig$slope_se,
                         pval_nominal_threshold = ER_gtex_eqtl_sig$pval_nominal_threshold,
                         min_pval_nominal = ER_gtex_eqtl_sig$min_pval_nominal,
                         pval_beta = ER_gtex_eqtl_sig$pval_beta,
                         stringsAsFactors = F)

aa_eqtl_gr <- makeGRangesFromDataFrame(aa_eqtl_df, keep.extra.columns = T)
options(scipen=999)
write.table(aa_eqtl_gr, 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/sig_eqtl.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
write.table(Positive_set_gr_list_1000bp[[4]], 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
write.table(Negative_set_gr_list_1000bp[[4]], 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
options(scipen=0)


aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_eqtl_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/sig_eqtl.bed"


aaa_eqtl_pos4 <- bedtools_intersect(bedfile_names = c(aa_eqtl_Add, aa_pos_Add), 
                                                wa=F, wb=F, loj=F, wo=F, wao=F,
                                                u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                                S=F, sorted=F, merge_after_distance=0, 
                                                output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/eqtl_pos_4.bed", 
                                                read_output = T)
aaa_eqtl_neg4 <- bedtools_intersect(bedfile_names = c(aa_eqtl_Add, aa_neg_Add),
                                                wa=F, wb=F, loj=F, wo=F, wao=F,
                                                u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                                S=F, sorted=F, merge_after_distance=0, 
                                                output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/eqtl_neg_4.bed",
                                    read_output = T)

GTEX_eqtl_pos4_neg4  <- c(aaa_eqtl_pos4, aaa_eqtl_neg4)
start(GTEX_eqtl_pos4_neg4) <- start(GTEX_eqtl_pos4_neg4) - 1
# Accepted_chromosomes <- paste0("chr", c(c(1:22), "X", "Y"))
# aadiff <- setdiff(levels(seqnames(Negative_set_gr_list[[4]])) , Accepted_chromosomes)
# Negative_set_gr_list[[4]] <- Negative_set_gr_list[[4]][! seqnames(Negative_set_gr_list[[4]]) %in% aadiff]
# 
# 
# 
# li_H3K27ac_gr <- readPeakFile(li_H3K27ac, as = "GRanges")
# li_ER_union_gr <- readPeakFile(li_ER_union, as = "GRanges")
# li_erna_only_gr <- readPeakFile(aa_tmp_erna_only_afterE2, as = "GRanges")
# vennplot(Sets = list(H3K27Ac = li_H3K27ac_gr,
#                      ER = li_ER_union_gr
#                      #eRNA_union = li_erna_union_gr
#                      ,eRNA = li_erna_only_gr
#                      #,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
# ),
# by = "Vennerable")
# 
################
# how many enhancers overlap 

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_eqtl_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/sig_eqtl.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_qtl <- readPeakFile(aa_eqtl_Add, as = "GRanges")

aa_ov_pos <- findOverlaps(query = aa_pos, subject = aa_qtl)
length(unique(aa_ov_pos@from))
hist(as.numeric(table(aa_ov_pos@from)), breaks = 100)

################################################################################################################
ER_pancan_eqtl_cis <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/BRCA_tumor.cis_eQTL.csv",
                                 stringsAsFactors = F, header = T, sep = ",")
ER_pancan_eqtl_trans <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/BRCA_tumor.trans_eQTL.csv",
                                   stringsAsFactors = F, header = T, sep = ",")

aaref <- unlist(lapply(strsplit(ER_pancan_eqtl_cis$SNP_alleles, split="/"), "[[", 1))
aaalt <- unlist(lapply(strsplit(ER_pancan_eqtl_cis$SNP_alleles, split="/"), "[[", 2))
ER_pancan_eqtl_cis_df19 <- data.frame(chr = ER_pancan_eqtl_cis$SNP_chr,
                                      start = ER_pancan_eqtl_cis$SNP_position,
                                      end = ER_pancan_eqtl_cis$SNP_position, 
                                      ref= aaref, alt = aaalt,
                                      snp_name = ER_pancan_eqtl_cis$eQTLs,
                                      gene = ER_pancan_eqtl_cis$egenes,
                                      stringsAsFactors = F)
ER_pancan_eqtl_cis_gr19 <- makeGRangesFromDataFrame(df = ER_pancan_eqtl_cis_df19,
                                                    keep.extra.columns = T)
aa_19_39_chain <- import.chain("~/Documents/Shayan/BioInf/Liftover_chain_files/hg19ToHg38.over.chain")
ER_pancan_eqtl_cis_gr38 <-  liftOver(ER_pancan_eqtl_cis_gr19, chain = aa_19_39_chain)
ER_pancan_eqtl_cis_gr38 <- unlist(ER_pancan_eqtl_cis_gr38)

aaref <- unlist(lapply(strsplit(ER_pancan_eqtl_trans$X5, split="/"), "[[", 1))
aaalt <- unlist(lapply(strsplit(ER_pancan_eqtl_trans$X5, split="/"), "[[", 2))
ER_pancan_eqtl_trans_df19 <- data.frame(chr = ER_pancan_eqtl_trans$X3,
                                      start = ER_pancan_eqtl_trans$X4,
                                      end = ER_pancan_eqtl_trans$X4, 
                                      ref= aaref,
                                      alt = aaalt,
                                      snp_name = ER_pancan_eqtl_trans$X2,
                                      gene = ER_pancan_eqtl_trans$X6,
                                      stringsAsFactors = F)
ER_pancan_eqtl_trans_gr19 <- makeGRangesFromDataFrame(df = ER_pancan_eqtl_trans_df19, 
                                                      keep.extra.columns = T)
aa_19_39_chain <- import.chain("~/Documents/Shayan/BioInf/Liftover_chain_files/hg19ToHg38.over.chain")
ER_pancan_eqtl_trans_gr38 <-  liftOver(ER_pancan_eqtl_trans_gr19, chain = aa_19_39_chain)
ER_pancan_eqtl_trans_gr38 <- unlist(ER_pancan_eqtl_trans_gr38)




# options(scipen=999)
# write.table(aa_eqtl_gr, 
#             file="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/sig_eqtl.bed",
#             quote=F, sep="\t", row.names=F, col.names=F)
# options(scipen=0)
################################################################################################################
################################################################################################################
# looking at all SNPs, are eqtls enriched among the snps that are predicted to have some effect?
aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"

#aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/00-common_all.vcf"
aa_snp_all <- read.table(file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/00-common_cut.txt",
                         header = F, stringsAsFactors = F)
aa_snp_all <- rbind(c(1, 10177, "rs367896724", "A", "AC"), aa_snp_all)
colnames(aa_snp_all) <- c("CHROM", "POS", "ID", "REF", "ALT")
aa_snp_all$POS <- as.numeric(aa_snp_all$POS)
aa_snp_df <- data.frame(chr = paste0("chr",aa_snp_all$CHROM),
                         start = aa_snp_all$POS, 
                         end = aa_snp_all$POS + 1,
                         ref = aa_snp_all$REF,
                         alt = aa_snp_all$ALT,
                        ID = aa_snp_all$ID,
                         stringsAsFactors = F)
#rownames(aa_snp_df) <- aa_snp_all$ID
SNP_common_all_df <- aa_snp_df
rm(aa_snp_df)
rm(aa_snp_all)
save(list = c("SNP_common_all_df"), 
     file = "~/Documents/Shayan/BioInf/EstrogenReceptor/SNP_common_all_df.RData")

aa_snp_gr <- makeGRangesFromDataFrame(SNP_common_all_df, keep.extra.columns = T)
options(scipen=999)
write.table(aa_snp_gr, 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
options(scipen=0)
aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps.bed"
aaa_snp_pos4 <- bedtools_intersect(bedfile_names = c(aa_snp_Add, aa_pos_Add), 
                                    wa=F, wb=F, loj=F, wo=F, wao=F,
                                    u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                    S=F, sorted=F, merge_after_distance=0, 
                                    output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/snps_pos_4.bed", 
                                    read_output = T)
aaa_snp_neg4 <- bedtools_intersect(bedfile_names = c(aa_snp_Add, aa_neg_Add), 
                                   wa=F, wb=F, loj=F, wo=F, wao=F,
                                   u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                   S=F, sorted=F, merge_after_distance=0, 
                                   output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/snps_neg_4.bed", 
                                   read_output = T)
aa_snp_pos_neg_4 <- c(aaa_snp_pos4, aaa_snp_neg4)
options(scipen=999)

write.table(aa_snp_pos_neg_4, 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
options(scipen=0)

head(aaa_snp_pos4)
aaa_snp_pos4df <- as.data.frame(aaa_snp_pos4)
aaa_snp_neg4df <- as.data.frame(aaa_snp_neg4)

colnames(aaa_snp_pos4df)[c(8,9,10)] <- c("ref", "alt", "id")
colnames(aaa_snp_neg4df)[c(8,9,10)] <- c("ref", "alt", "id")
#aaa_snp_pos4df$ref <- 
aaa_snp_pos4df$alt <-   as.character(levels(aaa_snp_pos4df$alt)[as.numeric(aaa_snp_pos4df$alt)])
aaa_snp_neg4df$alt <-   as.character(levels(aaa_snp_neg4df$alt)[as.numeric(aaa_snp_neg4df$alt)])

#aa_snp_pos4_all <- aa_snp_df[((aa_snp_df$chr == aaa_snp_pos4df$seqnames) & (aa_snp_df$start == aaa_snp_pos4df$start)), ]
# open up the oneswith more than one alternative
aa_eqtl_st <- unlist(lapply(aa_eqtl, "[[", 2))
aa_eqtl_ref <- unlist(lapply(aa_eqtl, "[[", 3))
aa_eqtl_alt <- unlist(lapply(aa_eqtl, "[[", 4))

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_snp <- readPeakFile(aa_snp_Add, as = "GRanges")
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_15")
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

aadf <- as.data.frame(aa_snp_pos_neg_4)
aa_names <- as.character(levels(aadf$V8)[as.numeric(aadf$V8)])
aaa <- c(aaa_snp_pos4df$alt, aaa_snp_neg4df$alt)
aatab <- table(aaa)
aaspl <- unlist(strsplit(names(aatab), split = ","))
aatab2 <- table(aaspl)
aanam <- names(aatab2)

aa_aanam <- unique(unlist(strsplit(aanam, split = "")))

aaa_snp_neg4df$ref <- as.character(levels(aaa_snp_neg4df$ref)[as.numeric(aaa_snp_neg4df$ref)])
aa_ref_Alt <- cbind(c(aaa_snp_pos4df$ref, aaa_snp_neg4df$ref), c(aaa_snp_pos4df$alt, aaa_snp_neg4df$alt))
aa_nam <- c(aaa_snp_pos4df$ref, aaa_snp_neg4df$ref)
aatab <- table(aa_nam)

unique(unlist(strsplit(c(aaa_snp_pos4df$ref, aaa_snp_neg4df$ref), split = "")))
aa_pos_neg_vars_snp <- write_variant_seq(enhancer_GR = aa_pos_neg, 
                                     eqtl_GR = aa_snp,
                                     eqtl_names = aa_names,
                                     eqtl_ref_alt = aa_ref_Alt,
                                     all_combs = F,
                                     my_genome = BSgenome.Hsapiens.UCSC.hg38, 
                                     label = c(rep(1, length(aa_pos)), rep(0, length(aa_neg))))

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")

# write GEMSTAT jobs to run
# write hal job for KDs
aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas1[1:25]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas2[1:25]]
aaadd3 <- union(aaadd1, aaadd2)

aa_par_name <- paste0(aaadd3, ".txt.Filter")
aa_seq_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_15/Variant_seq/",
                           recursive = T)
aa_lab_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_15/Variant_Labels/",
                           recursive = T)
aa_nam_spl <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_seq_files, split = "\\."), "[[", 1)), split = "_"), "[", c(3,4)), paste, collapse = "_"))


setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_15")
for(i in 1:length(aa_par_name)){
  for(j in 1:length(aa_seq_files)){
    cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr ",
          paste0("-s Variant_seq/", aa_seq_files[j]), 
          paste0("-e Variant_Labels/", aa_lab_files[j]),
          "-m motifs.wtmx -f TF_exp.tab", 
          paste0("-fo Varaint_out/", aa_nam_spl[j], "_", aaadd3[i], ".out"), 
          "-o DIRECT -c Coop/coop.par ", 
          paste0("-p Trained_par/", aa_par_name[i]),
          "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1\n"
    ), 
    sep = " ", append = (i!=1), file = "variant_jobs_snp.job")
  }
}

# read output and find the ones with some effect
# first read WT for all models

aa_wt_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out/",
                              pattern = "NA_NA_par*", full.names = T)
aa_wt_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out/",
                          pattern = "NA_NA_par*", full.names = F)
aa_nam_spl <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_wt_files, split = "\\."), "[[", 1)), split = "_"), "[", c(3,4)), paste, collapse = "_"))

aa_t <- read_output_train_test_GEMSTAT_indiv(output_file = aa_wt_files_ful[1])
aa_wt_models_snp <- matrix(nrow=nrow(aa_t$GT),
                       ncol=length(aa_wt_files))

colnames(aa_wt_models_snp) <- aa_nam_spl
rownames(aa_wt_models_snp) <- aa_t$pred$Rows
for(i in 1:ncol(aa_wt_models_snp)){
  aa_wt_models_snp[, i] <- read_output_train_test_GEMSTAT_indiv(output_file = aa_wt_files_ful[i])$pred$X1
}
aa_neg_out_snp <- list()
aa_pos_out_snp <- list()

aa_neg_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out/",
                               pattern = "neg_*", full.names = T)
aa_neg_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out/",
                           pattern = "neg_*", full.names = F)
aa_neg_sp_1 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_neg_files, split = "\\."), "[[", 1)), split = "_"), "[", c(1,2)), paste, collapse = "_"))
aa_neg_sp_2 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_neg_files, split = "\\."), "[[", 1)), split = "_"), "[", c(3,4)), paste, collapse = "_"))
aa_neg_uniq <- unique(aa_neg_sp_1)

aa_pos_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out/",
                               pattern = "pos_*", full.names = T)
aa_pos_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Varaint_out/",
                           pattern = "pos_*", full.names = F)
aa_pos_sp_1 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_pos_files, split = "\\."), "[[", 1)), split = "_"), "[", c(1,2)), paste, collapse = "_"))
aa_pos_sp_2 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_pos_files, split = "\\."), "[[", 1)), split = "_"), "[", c(3,4)), paste, collapse = "_"))
aa_pos_uniq <- unique(aa_pos_sp_1)

aa_neg_dif_snp <- list()
aa_pos_dif_snp <- list()

for(i in 1:length(aa_neg_uniq)){
  print(i)
  # aa_cur_ind <- which(aa_neg_sp_1 %in% aa_neg_uniq[i])
  # aa_t <- read.table(file = aa_neg_files_ful[aa_cur_ind[1]], header = T, stringsAsFactors = F)
  # aa_ind <- seq(2, nrow(aa_t), 2)
  # aa_t <- aa_t[aa_ind, ]
  # aa_neg_out_snp[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  # colnames(aa_neg_out_snp[[i]]) <- colnames(aa_wt_models_snp)
  # rownames(aa_neg_out_snp[[i]]) <- aa_t$Rows
  # aa_neg_dif_snp[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  # colnames(aa_neg_dif_snp[[i]]) <- colnames(aa_wt_models_snp)
  # rownames(aa_neg_dif_snp[[i]]) <- aa_t$Rows
  # for(j in 1:length(aa_cur_ind)){
  #   aa_table <- read.table(file =aa_neg_files_ful[aa_cur_ind[j]], 
  #                          header = T, stringsAsFactors = F)
  #   aa_neg_out_snp[[i]][, aa_neg_sp_2[aa_cur_ind[j]]] <- aa_table[seq(2, nrow(aa_table), 2),"X1"]
  # }
  for(j in 1:nrow(aa_neg_out_snp[[i]])){
    aa_neg_dif_snp[[i]][j, ] <- (aa_neg_out_snp[[i]][j, ] - Exp12_WT_model_prediction_Filtered[aa_neg_uniq[i], ])/ Exp12_WT_model_prediction_Filtered[aa_neg_uniq[i], ]
    
  }
}
names(aa_neg_out_snp) <- aa_neg_uniq
names(aa_neg_dif_snp) <- aa_neg_uniq

for(i in 1:length(aa_pos_uniq)){
  print(i)
  # aa_cur_ind <- which(aa_pos_sp_1 %in% aa_pos_uniq[i])
  # aa_t <- read.table(file = aa_pos_files_ful[aa_cur_ind[1]], header = T, stringsAsFactors = F)
  # aa_ind <- seq(2, nrow(aa_t), 2)
  # aa_t <- aa_t[aa_ind, ]
  # aa_pos_out_snp[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))  
  # colnames(aa_pos_out_snp[[i]]) <- colnames(aa_wt_models_snp)
  # rownames(aa_pos_out_snp[[i]]) <- aa_t$Rows
  # aa_pos_dif_snp[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  # colnames(aa_pos_dif_snp[[i]]) <- colnames(aa_wt_models_snp)
  # rownames(aa_pos_dif_snp[[i]]) <- aa_t$Rows
  # for(j in 1:length(aa_cur_ind)){
  #   aa_table <- read.table(file =aa_pos_files_ful[aa_cur_ind[j]], 
  #                          header = T, stringsAsFactors = F)
  #   aa_pos_out_snp[[i]][, aa_pos_sp_2[aa_cur_ind[j]]] <- aa_table[seq(2, nrow(aa_table), 2),"X1"]
  # }
  for(j in 1:nrow(aa_pos_out_snp[[i]])){
    aa_pos_dif_snp[[i]][j, ] <- (aa_pos_out_snp[[i]][j, ] - Exp12_WT_model_prediction_Filtered[aa_pos_uniq[i], ])/ Exp12_WT_model_prediction_Filtered[aa_pos_uniq[i], ]
  }
  
}
names(aa_pos_out_snp) <- aa_pos_uniq
names(aa_pos_dif_snp) <- aa_pos_uniq
aa_pos_dif_snp$pos_1[2, 20:40]
(aa_pos_out_snp$pos_1[2,20:40] - aa_wt_models_snp["pos_1", 20:40])/(aa_wt_models_snp["pos_1", 20:40])
max(aa_pos_dif_snp$pos_1, na.rm = T)

# aa_pos_dif_snp_old_noFilt <- aa_pos_dif_snp # before filtering the outputs of WT models
# aa_neg_dif_snp_old_noFilt <- aa_neg_dif_snp
aa_pos_dif_all_snp <- do.call(what = rbind, aa_pos_dif_snp)
boxplot.matrix(t(aa_pos_dif_all_snp[1:5,])*100)

aa_neg_dif_all_snp <- do.call(what = rbind, aa_neg_dif_snp)
#boxplot.matrix(t(aa_neg_dif_all_snp))




z <- zClust(x=t(aa_pos_dif_all_snp[,-c(17)]), scale="row", zlim=c(-3,3), method="average")
# require(RColorBrewer)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(9)
heatmap.2(z$data, trace='none',
          col=rev(cols),
          #Rowv=F,Colv=F,
          main="predicted change"
          ,Rowv=z$Rowv, 
          Colv=z$Colv
)
z <- zClust(x=t(aa_neg_dif_all[,-c(17)]), scale="row", zlim=c(-3,3), method="average")
# require(RColorBrewer)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(9)
heatmap.2(z$data, trace='none',
          col=rev(cols),
          Rowv=F,Colv=F,
          main="predicted change neg"
          # ,Rowv=z$Rowv, 
          # Colv=z$Colv
)
# find enhancers containing multiple eqtls where some don't affect expression and some do


aa_pos_dif_filtered_snp <- aa_pos_dif_snp
for(i in length(aa_pos_dif_filtered_snp):1){
  if(max(abs(aa_pos_dif_filtered_snp[[i]]), na.rm = T) < 0.05){
    aa_pos_dif_filtered_snp[[i]] <- NULL
  }
}
aanrwo <- unlist(lapply(aa_pos_dif_filtered_snp, nrow))
par(mar = c(1,1,1,1), mfrow = c(10, 10))
for(i in 1:length(aa_pos_dif_filtered_snp)){
  boxplot.matrix(t(aa_pos_dif_filtered_snp[[i]]), 
                 las= 2, 
                 main = names(aa_pos_dif_filtered_snp)[i],
                 xaxt = "n", outline = F)
  abline(h= seq(-3, 3, 0.05), col = 2, lwd = 0.5, lty = 3)
  abline(h= 0, col = 2, lwd = 1, lty = 1)
}
# boxplot.matrix(t(aa_pos_dif_filtered_snp$pos_312), 
#                las= 2, main = "pos_312 predicted change")
# View(aa_pos_dif_filtered_snp$pos_171)

aa_neg_dif_filtered_snp <- aa_neg_dif_snp
for(i in length(aa_neg_dif_filtered_snp):1){
  if(max(abs(aa_neg_dif_filtered_snp[[i]]), na.rm = T) < 0.05){
    aa_neg_dif_filtered_snp[[i]] <- NULL
  }
}
aanrwo <- unlist(lapply(aa_neg_dif_filtered_snp, nrow))
par(mar = c(1,1,1,1), mfrow = c(10, 7))
for(i in 1:length(aa_neg_dif_filtered_snp)){
  if(aanrwo[i] > 1){
    boxplot.matrix(t(aa_neg_dif_filtered_snp[[i]]), 
                   las= 2, 
                   main = names(aa_neg_dif_filtered_snp)[i],
                   xaxt = "n")
    abline(h= seq(-3, 3, 0.05), col = 2, lwd = 0.5, lty = 3)
    abline(h= 0, col = 2, lwd = 1, lty = 1)
  }
  
}
################################################################################################################
# read in GWAS SNPs, and mark snps that are GWAS
GWAS_breast_Cancer <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/gwas-association-downloaded_2020-01-21-EFO_0000305-withChildTraits.tsv",
                                 header = T, stringsAsFactors = F)
# GWAS_breast_Cancer$SNPS[1:15]
# aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
# aa_snp <- readPeakFile(aa_snp_Add, as = "GRanges")
# aa_snp_names <- aa_snp$V10
# aa_snp_names <- as.character(levels(aa_snp_names)[as.numeric(aa_snp_names)])
# which(aa_snp_names %in% GWAS_breast_Cancer$SNPS)
# aa_GWAS_SNPS <- aa_snp_names[which(aa_snp_names %in% GWAS_breast_Cancer$SNPS)]
# 
aa <- as.numeric(GWAS_breast_Cancer$CHR_POS)
GWAS_breast_Cancer <- GWAS_breast_Cancer[!is.na(aa),]

GWAS_breast_Cancer_df38 <- data.frame(chr = paste0("chr", GWAS_breast_Cancer$CHR_ID),
                                      start = GWAS_breast_Cancer$CHR_POS,
                                      end = GWAS_breast_Cancer$CHR_POS, 
                                      #ref= aaref, alt = aaalt,
                                      snp_name = GWAS_breast_Cancer$SNPS,
                                      #gene = ER_pancan_eqtl_cis$egenes,
                                      stringsAsFactors = F)

GWAS_breast_Cancer_gr38 <- makeGRangesFromDataFrame(df = GWAS_breast_Cancer_df38,
                                                    keep.extra.columns = T)
################################################################################################################
# mark which SNPs are eqtls
# get the intersection of SNPs and eqtls
head(ER_gtex_eqtl_sig)

aa_eqtl_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/sig_eqtl.bed"
aa_qtl <- readPeakFile(aa_eqtl_Add, as = "GRanges")
ranges(aa_qtl) <- ranges(aa_qtl) + 1
aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_snp <- readPeakFile(aa_snp_Add, as = "GRanges")
ranges(aa_snp) <- ranges(aa_snp) + 1
aa_overlap_holder <- findOverlaps(query = aa_snp, subject = aa_qtl)
#ol_enhancers_index <- unique(overlap_holder@from)
aa_snp_qtl_marker <- numeric(length(aa_snp))
aa_snp_qtl_marker[aa_overlap_holder@from] <- 1
#aa_snp$V10 <- as.character(levels(aa_snp$V10)[as.numeric(aa_snp$V10)])
names(aa_snp_qtl_marker) <- aa_snp$V10
aa_snp_gwas_marker <- numeric(length(aa_snp))
names(aa_snp_gwas_marker) <-  aa_snp$V10
aa_snp_gwas_marker[aa_GWAS_SNPS] <- 1


aa_neg_dif_all_snp_names <- unlist(lapply(strsplit(rownames(aa_neg_dif_all_snp), split = "_"), "[[", 2))
aa_pos_dif_all_snp_names <- unlist(lapply(strsplit(rownames(aa_pos_dif_all_snp), split = "_"), "[[", 2))

aa <- aa_snp_qtl_marker[match(aa_neg_dif_all_snp_names, names(aa_snp_qtl_marker))]
aa2 <- aa_snp_gwas_marker[match(aa_neg_dif_all_snp_names, names(aa_snp_gwas_marker))]
aa_neg_dif_all_snp <- cbind(aa_neg_dif_all_snp, aa)
aa_neg_dif_all_snp <- cbind(aa_neg_dif_all_snp, aa2)
aa <- aa_snp_qtl_marker[match(aa_pos_dif_all_snp_names, names(aa_snp_qtl_marker))]
aa2 <- aa_snp_gwas_marker[match(aa_pos_dif_all_snp_names, names(aa_snp_gwas_marker))]
aa_pos_dif_all_snp <- cbind(aa_pos_dif_all_snp, aa)
aa_pos_dif_all_snp <- cbind(aa_pos_dif_all_snp, aa2)

#aa_pos_dif_all_snp <- aa_pos_dif_all_snp[,-c(17)]
#aa_neg_dif_all_snp <- aa_neg_dif_all_snp[,-c(17)]
par(mfrow = c(6,7), mar = c(1,1,1,1))
for(i in 1:(ncol(aa_pos_dif_all_snp)-1)){
  boxplot(aa_pos_dif_all_snp[,i]~aa_pos_dif_all_snp[, ncol(aa_pos_dif_all_snp)], main = colnames(aa_pos_dif_all_snp)[i])
}

par(mfrow = c(6,7), mar = c(1,1,1,1))
for(i in 1:(ncol(aa_neg_dif_all_snp)-1)){
  boxplot(aa_neg_dif_all_snp[,i]~aa_neg_dif_all_snp[, ncol(aa_neg_dif_all_snp)], main = colnames(aa_neg_dif_all_snp)[i])
}
# for several thresholds count the number of preds above and below by eqtl
aa_pos_sig_nonsig <- matrix(nrow = 6, ncol = ncol(aa_pos_dif_all_snp)-2)
aa_neg_sig_nonsig <- matrix(nrow = 6, ncol = ncol(aa_neg_dif_all_snp)-2)
colnames(aa_neg_sig_nonsig) <- colnames(aa_neg_dif_all_snp)[1:(ncol(aa_neg_dif_all_snp)-2)]
colnames(aa_pos_sig_nonsig) <- colnames(aa_pos_dif_all_snp)[1:(ncol(aa_pos_dif_all_snp)-2)]

rownames(aa_pos_sig_nonsig) <- c("sig_all", "sig_qtl", "nonsig_all", "nonsig_qtl", "pval", "gwas_sig")
rownames(aa_pos_sig_nonsig) <- c("sig_all", "sig_qtl", "nonsig_all", "nonsig_qtl", "pval", "gwas_sig")

aa_thresh <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
aa_pos_sig_hold <- list()
aa_neg_sig_hold <- list()
for(i in 1:length(aa_thresh)){
  aa_pos_sig_hold[[i]] <- aa_pos_sig_nonsig
  aa_neg_sig_hold[[i]] <- aa_neg_sig_nonsig
  for(j in 1:(ncol(aa_pos_sig_nonsig)-1)){
    aa_pos_sig_hold[[i]][1,j] <- sum(abs(aa_pos_dif_all_snp[,j]) >= aa_thresh[i], na.rm = T)
    aa_pos_sig_hold[[i]][2,j] <- sum(aa_pos_dif_all_snp[abs(aa_pos_dif_all_snp[,j]) >= aa_thresh[i],( ncol(aa_pos_dif_all_snp) - 1)], na.rm = T)
    aa_pos_sig_hold[[i]][3,j] <- sum(abs(aa_pos_dif_all_snp[,j]) < aa_thresh[i], na.rm = T)
    aa_pos_sig_hold[[i]][4,j] <- sum(aa_pos_dif_all_snp[abs(aa_pos_dif_all_snp[,j]) < aa_thresh[i], (ncol(aa_pos_dif_all_snp)-1)], na.rm = T)
    aa_pos_sig_hold[[i]][5,j] <- phyper(q = aa_pos_sig_hold[[i]][2,j],
                                       m = sum(aa_pos_dif_all_snp[,(ncol(aa_pos_dif_all_snp) - 1)] == 1),
                                       n = sum(aa_pos_dif_all_snp[,(ncol(aa_pos_dif_all_snp) - 1)] == 0), 
                                       k = aa_pos_sig_hold[[i]][1,j],lower.tail = F)
    aa_pos_sig_hold[[i]][6,j] <- sum(aa_pos_dif_all_snp[abs(aa_pos_dif_all_snp[,j]) >= aa_thresh[i], ncol(aa_pos_dif_all_snp)], na.rm = T)
    
    aa_neg_sig_hold[[i]][1,j] <- sum(abs(aa_neg_dif_all_snp[,j]) >= aa_thresh[i], na.rm = T)
    aa_neg_sig_hold[[i]][2,j] <- sum(aa_neg_dif_all_snp[abs(aa_neg_dif_all_snp[,j]) >= aa_thresh[i], (ncol(aa_neg_dif_all_snp)-1)], na.rm = T)
    aa_neg_sig_hold[[i]][3,j] <- sum(abs(aa_neg_dif_all_snp[,j]) < aa_thresh[i], na.rm = T)
    aa_neg_sig_hold[[i]][4,j] <- sum(aa_neg_dif_all_snp[abs(aa_neg_dif_all_snp[,j]) < aa_thresh[i], (ncol(aa_neg_dif_all_snp)-1)], na.rm = T)
    aa_neg_sig_hold[[i]][5,j] <- phyper(q = aa_neg_sig_hold[[i]][2,j],
                                        m = sum(aa_neg_dif_all_snp[,(ncol(aa_neg_dif_all_snp)-1)] == 1),
                                        n = sum(aa_neg_dif_all_snp[,(ncol(aa_neg_dif_all_snp)-1)] == 0), 
                                        k = aa_neg_sig_hold[[i]][1,j],lower.tail = F)
    aa_neg_sig_hold[[i]][6,j] <- sum(aa_neg_dif_all_snp[abs(aa_neg_dif_all_snp[,j]) >= aa_thresh[i], ncol(aa_neg_dif_all_snp)], na.rm = T)
  }
}
names(aa_pos_sig_hold) <- aa_thresh
names(aa_neg_sig_hold) <- aa_thresh
View(aa_pos_sig_hold$`0.05`)
aa_pos_sig_hold[[2]]
phyper(q = aa_pos_sig_hold[[i]][2,j],
       m = sum(aa_pos_dif_all_snp[,ncol(aa_pos_dif_all_snp)] == 1),
       n = sum(aa_pos_dif_all_snp[,ncol(aa_pos_dif_all_snp)] == 0), 
       k = aa_pos_sig_hold[[i]][1,j])

aa_neg_sig_hold$`1`[5,]



for(i in 1:length(aa_pos_sig_hold)){
  print(names(aa_pos_sig_hold)[i])
  print("pos")
  print(sum(((aa_pos_sig_hold$`0.5`)[[i]][5,] < 0.05) & (aa_pos_sig_hold[[i]][2,] > 0), na.rm =T))
  print("gwas_pos")
  print(sum(aa_pos_sig_hold[[i]][6,] > 0, na.rm=T))
  print("neg")
  print(sum((aa_neg_sig_hold[[i]][5,]< 0.05) & (aa_neg_sig_hold[[i]][2,] > 0), na.rm =T) )
  print("gwas_neg")
  print(sum(aa_neg_sig_hold[[i]][6,] > 0, na.rm=T))
}
View(aa_pos_sig_hold$`0.01`)
View(aa_pos_sig_hold$`0.5`)
View(aa_pos_sig_hold$`0.9`)
View(aa_pos_sig_hold$`1`)
View(aa_neg_sig_hold$`0.01`)
View(aa_neg_sig_hold$`0.8`)
View(aa_neg_sig_hold$`1`)

# looking at pos and negs together and computing hypergeom pval for combined
aa_all_sig_nonsig <- matrix(nrow = 6, ncol = ncol(aa_pos_dif_all_snp)-2)
colnames(aa_all_sig_nonsig) <- colnames(aa_neg_dif_all_snp)[1:(ncol(aa_neg_dif_all_snp)-2)]

aa_all_dif_snp <- rbind(aa_pos_dif_all_snp, aa_neg_dif_all_snp)
rownames(aa_all_sig_nonsig) <- c("sig_all", "sig_qtl", "nonsig_all", "nonsig_qtl", "pval", "gwas_sig")

aa_thresh <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
aa_all_sig_hold <- list()

for(i in 1:length(aa_thresh)){
  aa_all_sig_hold[[i]] <- aa_all_sig_nonsig
  for(j in 1:(ncol(aa_pos_sig_nonsig)-1)){
    aa_all_sig_hold[[i]][1,j] <- sum(abs(aa_all_dif_snp[,j]) >= aa_thresh[i], na.rm = T)
    aa_all_sig_hold[[i]][2,j] <- sum(aa_all_dif_snp[abs(aa_all_dif_snp[,j]) >= aa_thresh[i],( ncol(aa_all_dif_snp) - 1)], na.rm = T)
    aa_all_sig_hold[[i]][3,j] <- sum(abs(aa_all_dif_snp[,j]) < aa_thresh[i], na.rm = T)
    aa_all_sig_hold[[i]][4,j] <- sum(aa_all_dif_snp[abs(aa_all_dif_snp[,j]) < aa_thresh[i], (ncol(aa_all_dif_snp)-1)], na.rm = T)
    aa_all_sig_hold[[i]][5,j] <- phyper(q = aa_all_sig_hold[[i]][2,j],
                                        m = sum(aa_all_dif_snp[,(ncol(aa_all_dif_snp) - 1)] == 1),
                                        n = sum(aa_all_dif_snp[,(ncol(aa_all_dif_snp) - 1)] == 0), 
                                        k = aa_all_sig_hold[[i]][1,j],lower.tail = F)
    aa_all_sig_hold[[i]][6,j] <- sum(aa_all_dif_snp[abs(aa_all_dif_snp[,j]) >= aa_thresh[i], ncol(aa_all_dif_snp)], na.rm = T)
    
  }
}
names(aa_all_sig_hold) <- aa_thresh

for(i in 1:length(aa_all_sig_hold)){
  print(names(aa_all_sig_hold)[i])
  print(sum((aa_all_sig_hold[[i]][5,] < 0.05) & (aa_all_sig_hold[[i]][2,] > 0), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold[[i]][6,] > 0, na.rm=T))
}
View(aa_all_sig_hold$`0.7`)
sort(aa_all_sig_hold$`0.8`[5,], decreasing = F)

aaw <- which(aa_all_dif_snp[, "par_2646"] > 0.8 & aa_all_dif_snp[, (ncol(aa_all_dif_snp) - 1)] == 1)
#aaw <- which(aa_all_dif_snp[, "par_2646"] > 0.8 & aa_all_dif_snp[, (ncol(aa_all_dif_snp))] == 1)
aaw2 <- unlist(lapply(strsplit(x = names(aaw), split="_"), "[[", 2))

aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_snp <- readPeakFile(aa_snp_Add, as = "GRanges")
ranges(aa_snp) <- ranges(aa_snp) + 1
aa_snp <- as.data.frame(aa_snp)
aa_snp$V10 <- as.character(levels(aa_snp$V10)[as.numeric(aa_snp$V10)])
aa_snp$V9 <- as.character(levels(aa_snp$V9)[as.numeric(aa_snp$V9)])
aa_snp$V8 <- as.character(levels(aa_snp$V8)[as.numeric(aa_snp$V8)])
aa_snp$seqnames <- as.character(levels(aa_snp$seqnames)[as.numeric(aa_snp$seqnames)])


aa_snp2 <- aa_snp[aa_snp$V10 %in% aaw2,]

# aa_all_names <- c(rownames(aa_pos_dif_filtered_all2 ), rownames(aa_neg_dif_filtered_all2))
# aa_all_names <- unlist(lapply(strsplit(aa_all_names, split = "_"), "[[", 2))
#aa_all_qtl <- ER_gtex_eqtl_sig[as.numeric(aa_all_names),]
#aa_all_qtl_uniq <- unique(aa_all_qtl$variant_id)
#aa_eqtl <-  aa_snp$V1
aa_eqtl_chr <- aa_snp2$seqnames
aa_eqtl_st <- aa_snp2$start
aa_eqtl_ref <- aa_snp2$V8
aa_eqtl_alt <- aa_snp2$V9
aan <-   aa_snp2$V10
aa_eqtl_alt2 <- strsplit(aa_eqtl_alt, split = ",")
aarank <- as.numeric(unlist(lapply(strsplit(names(aaw), split = "_"), "[[", 3)))
names(aarank) <- aaw2
for(i in 1:length(aa_eqtl_alt2)){
  aa_eqtl_alt2[[i]] <- aa_eqtl_alt2[[i]][aarank[aan[i]]]
}
aa_eqtl_alt2 <- unlist(aa_eqtl_alt2)


aamin_LLR_low <- aamin_LLR - 0.1*aamin_LLR
aamin_LLR_lowest <- numeric(length(aamin_LLR))
names(aamin_LLR_lowest) <- names(aamin_LLR_low)
aa_TF.motifs.Expanded_pseudo_exp12_t_ps <- lapply(TF.motifs.Expanded_pseudo_exp12_t, Addpsudo)
aa_tf_res_eqtl2646 <- list()
for(i in 1:length(aa_eqtl_alt2)){
  print(i)
  aa_tf_res_eqtl2646[[i]] <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                                     my_motifs = aa_TF.motifs.Expanded_pseudo_exp12_t_ps,
                                     snp_chr = aa_eqtl_chr[i],
                                     snp_start = aa_eqtl_st[i]-1,
                                     snp_ref = aa_eqtl_ref[i],
                                     snp_alt = aa_eqtl_alt2[i],
                                     min_LLR = aamin_LLR[aa_names])
  
}
names(aa_tf_res_eqtl2646)  <- aa_snp2$V10

aa_tf_res_eqtl2646$rs3809260

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
# read in interaction data and convert from hg19 t0 hg38
MCF7_interaction_hg19 <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/ConformationData/4DGenome_HomoSapiens_hg19.txt",
                                    header = T,
                                    stringsAsFactors = F)
MCF7_interaction_hg19 <- MCF7_interaction_hg19[MCF7_interaction_hg19$Cell.Tissue == "MCF7",]

aa_1 <-MCF7_interaction_hg19[,c(c(1:3), c(7), c(11:15))] 
colnames(aa_1)[1:3] <- c("chr", "start", "end")
aa_2 <-MCF7_interaction_hg19[,c(c(4:6), c(8), c(11:15))] 
colnames(aa_2)[1:3] <- c("chr", "start", "end")


aa_1 <- makeGRangesFromDataFrame(aa_1,keep.extra.columns = T)
aa_2 <- makeGRangesFromDataFrame(aa_2,keep.extra.columns = T)
aa_19_39_chain <- import.chain("~/Documents/Shayan/BioInf/Liftover_chain_files/hg19ToHg38.over.chain")

MCF7_interaction_hg38 <- interaction_liftover(int_a_gr = aa_1, int_b_gr = aa_2, my_chain = aa_19_39_chain)

ReMapChIP.GRanges.list

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
#aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
#aa_snp <- readPeakFile(aa_snp_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
#start(aa_snp) <- start(aa_snp) - 1
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

Enhancer.ReMapchip.Overlap.byEnhancer #a list where each entry is the overlap list of the enhancers with one chip dataset

Enhancer.ReMapchip.Overlap.byEnhancer <- OverlapInteractionBatchExtract(subjectCoord1=aa_pos_neg,
                                                                        subjectCoord2List=ReMapChIP.GRanges.list, 
                                                                        interactingA = MCF7_interaction_hg38$new_A, 
                                                                        interactingB = MCF7_interaction_hg38$new_B,
                                                                        mode = "both")
Enhancer.ReMapchip.Overlap.byEnhancer$OverlapMat[1:5,1:5]
# for(i in 1:length(ReMapChIP.GRanges.list)){
#   Enhancer.ReMapchip.Overlap.byEnhancer[[i]] <- overlapExtractor(subjectCoord1 = aa_pos_neg,
#                                                                  subjectCoord2 =ReMapChIP.GRanges.list[[i]])
# }
#names(Enhancer.ReMapchip.Overlap.byEnhancer) <- names(ReMapChIP.GRanges.list)
# removing IMPET infered interactions and re-evaluating results
table(MCF7_interaction_hg38$new_A$Detection_Method)

MCF7_interaction_hg38_expOnly <- list()
MCF7_interaction_hg38_expOnly$new_A <- MCF7_interaction_hg38$new_A[! MCF7_interaction_hg38$new_A$Detection_Method == "IM-PET"]
MCF7_interaction_hg38_expOnly$new_B <- MCF7_interaction_hg38$new_B[! MCF7_interaction_hg38$new_B$Detection_Method == "IM-PET"]


aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
#aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
#aa_snp <- readPeakFile(aa_snp_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
#start(aa_snp) <- start(aa_snp) - 1
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

Enhancer.ReMapchip.Overlap.byEnhancer_experimental_only <- OverlapInteractionBatchExtract(subjectCoord1=aa_pos_neg,
                                                                                          subjectCoord2List=ReMapChIP.GRanges.list, 
                                                                                          interactingA = MCF7_interaction_hg38_expOnly$new_A, 
                                                                                          interactingB = MCF7_interaction_hg38_expOnly$new_B,
                                                                                          mode = "both")

rownames(Enhancer.ReMapchip.Overlap.byEnhancer_experimental_only$IntMat) <- names(aa_pos_neg)
aa_chip_int <- list()
for(i in 1:nrow(Enhancer.ReMapchip.Overlap.byEnhancer_experimental_only$IntMat)){
  aaw <- which(Enhancer.ReMapchip.Overlap.byEnhancer_experimental_only$IntMat[i, ] > 0)
  aa_chip_int[[i]] <- colnames(Enhancer.ReMapchip.Overlap.byEnhancer_experimental_only$IntMat)[aaw]
}
names(aa_chip_int) <- names(aa_pos_neg)
Enhancer_TF_interact_name_expOnly <- aa_chip_int
hist(unlist(lapply(Enhancer_TF_interact_name_expOnly, length)))

aa_chip_ovl <- list()
for(i in 1:nrow(Enhancer.ReMapchip.Overlap.byEnhancer_experimental_only$OverlapMat)){
  aaw <- which(Enhancer.ReMapchip.Overlap.byEnhancer_experimental_only$OverlapMat[i, ] > 0)
  aa_chip_ovl[[i]] <- colnames(Enhancer.ReMapchip.Overlap.byEnhancer_experimental_only$OverlapMat)[aaw]
}
names(aa_chip_ovl) <- names(aa_pos_neg)
Enhancer_TF_overlap_name_expOnly <- aa_chip_ovl
hist(unlist(lapply(Enhancer_TF_overlap_name_expOnly, length)))

Enhancer_TF_overlap_interact_name_expOnly <- list()
for(i in 1:length(Enhancer_TF_overlap_name_expOnly)){
  Enhancer_TF_overlap_interact_name_expOnly[[i]] <- union(Enhancer_TF_overlap_name_expOnly[[i]],
                                                          Enhancer_TF_interact_name_expOnly[[i]])
}
names(Enhancer_TF_overlap_interact_name_expOnly) <- names(Enhancer_TF_overlap_name_expOnly)
hist(unlist(lapply(Enhancer_TF_overlap_interact_name_expOnly, length)))
aaTF_enh_chip <- data.frame(Enhancer = character(0), TF = character(0))
for(i in 1:length(Enhancer_TF_overlap_interact_name_expOnly)){
  if(length(Enhancer_TF_overlap_interact_name_expOnly[[i]]) > 0){
    aacur <- data.frame(Enhancer =rep(names(Enhancer_TF_overlap_interact_name_expOnly)[i],
                                      length(Enhancer_TF_overlap_interact_name_expOnly[[i]]) ),
                        TF = Enhancer_TF_overlap_interact_name_expOnly[[i]])
    aaTF_enh_chip <- rbind(aaTF_enh_chip, aacur)
    
  }
}
aaTF_enh_chip$Enhancer <- levels(aaTF_enh_chip$Enhancer)[as.numeric(aaTF_enh_chip$Enhancer)]
aaTF_enh_chip$TF <- levels(aaTF_enh_chip$TF)[as.numeric(aaTF_enh_chip$TF)]
TF_enhancer_chip_network <- aaTF_enh_chip

#####################################################################################################################
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
genesAll <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
Promoter1kbdf_gene <- as.data.frame(resize(genesAll, width=1, fix='start'))
#make Grange and DataFrame for 
#1: promoter (1 kb upstream of tss)
Promoter1kbdf_gene$start <- Promoter1kbdf_gene$start - 1000
Promoter1kbdf_gene <- StartEndCoordinateCheck(inputRangedf = Promoter1kbdf_gene)
Promoter1kbGR_gene <- makeGRangesFromDataFrame(Promoter1kbdf_gene)
aa_list <- list()
aa_list[[1]] <- Promoter1kbGR_gene

enhancer_gene_ovlap_interact <- OverlapInteractionBatchExtract(subjectCoord1=aa_pos_neg,
                                                               subjectCoord2List=aa_list, 
                                                               interactingA = MCF7_interaction_hg38$new_A, 
                                                               interactingB = MCF7_interaction_hg38$new_B,
                                                               mode = "both")
names(enhancer_gene_ovlap_interact$OverlapList[[1]]) <- paste0(names(aa_pos_neg), "_")
names(enhancer_gene_ovlap_interact$IntList[[1]]) <- paste0(names(aa_pos_neg), "_")

aa_enh_geneprom_ov <- unlist(enhancer_gene_ovlap_interact$OverlapList[[1]])
aa_enh_geneprom_in <- unlist(enhancer_gene_ovlap_interact$IntList[[1]])

Enhancer_genePromoter_overlap_mat <- data.frame(enh_name = names(aa_enh_geneprom_ov),
                                         gene_entID = Promoter1kbdf_gene$gene_id[aa_enh_geneprom_ov],
                                         stringsAsFactors = F)
aaaaa <- unlist(lapply(lapply(strsplit(Enhancer_genePromoter_overlap_mat$enh_name, "_"), "[", c(1,2)), paste, collapse = "_"))
Enhancer_genePromoter_overlap_mat$enh_name <- aaaaa

Enhancer_genePromoter_interac_mat <- data.frame(enh_name = names(aa_enh_geneprom_in),
                                                gene_entID = Promoter1kbdf_gene$gene_id[aa_enh_geneprom_in],
                                                stringsAsFactors = F)
aaaaa <- unlist(lapply(lapply(strsplit(Enhancer_genePromoter_interac_mat$enh_name, "_"), "[", c(1,2)), paste, collapse = "_"))
Enhancer_genePromoter_interac_mat$enh_name <- aaaaa

enhancer_gene_ovlap_interact_experimental_only <- OverlapInteractionBatchExtract(subjectCoord1=aa_pos_neg,
                                                               subjectCoord2List=aa_list, 
                                                               interactingA = MCF7_interaction_hg38_expOnly$new_A, 
                                                               interactingB = MCF7_interaction_hg38_expOnly$new_B,
                                                               mode = "both")

names(enhancer_gene_ovlap_interact_experimental_only$OverlapList[[1]]) <- paste0(names(aa_pos_neg), "_")
names(enhancer_gene_ovlap_interact_experimental_only$IntList[[1]]) <- paste0(names(aa_pos_neg), "_")

aa_enh_geneprom_ov <- unlist(enhancer_gene_ovlap_interact_experimental_only$OverlapList[[1]])
aa_enh_geneprom_in <- unlist(enhancer_gene_ovlap_interact_experimental_only$IntList[[1]])

# Enhancer_genePromoter_overlap_mat <- data.frame(enh_name = names(aa_enh_geneprom_ov),
#                                                 gene_entID = Promoter1kbdf_gene$gene_id[aa_enh_geneprom_ov],
#                                                 stringsAsFactors = F)
# aaaaa <- unlist(lapply(lapply(strsplit(Enhancer_genePromoter_overlap_mat$enh_name, "_"), "[", c(1,2)), paste, collapse = "_"))
# Enhancer_genePromoter_overlap_mat$enh_name <- aaaaa

Enhancer_genePromoter_interac_mat_experimental_only <- data.frame(enh_name = names(aa_enh_geneprom_in),
                                                gene_entID = Promoter1kbdf_gene$gene_id[aa_enh_geneprom_in],
                                                stringsAsFactors = F)
aaaaa <- unlist(lapply(lapply(strsplit(Enhancer_genePromoter_interac_mat_experimental_only$enh_name, "_"), "[", c(1,2)), paste, collapse = "_"))
Enhancer_genePromoter_interac_mat_experimental_only$enh_name <- aaaaa

#####################################################################################################################
Vicinity20kbdf_gene <- as.data.frame(resize(genesAll, width=1, fix='start'))
Vicinity20kbdf_gene$end <- Vicinity20kbdf_gene$start + 10000
Vicinity20kbdf_gene$start <- Vicinity20kbdf_gene$start - 10000
Vicinity20kbdf_gene <- StartEndCoordinateCheck(inputRangedf = Vicinity20kbdf_gene)
Vicinity20kbGR_gene <- makeGRangesFromDataFrame(Vicinity20kbdf_gene)

aa_list <- list()
aa_list[[1]] <- Vicinity20kbGR_gene

enhancer_gene_vicinity_ovlap <- OverlapInteractionBatchExtract(subjectCoord1=aa_pos_neg,
                                                               subjectCoord2List=aa_list, 
                                                               interactingA = MCF7_interaction_hg38$new_A, 
                                                               interactingB = MCF7_interaction_hg38$new_B,
                                                               mode = "overlap")

names(enhancer_gene_vicinity_ovlap$OverlapList[[1]]) <- paste0(names(aa_pos_neg), "_")
aa_enh_gene_vicinity <- unlist(enhancer_gene_vicinity_ovlap$OverlapList[[1]])
Enhancer_gene_20kb_vicinity_mat <- data.frame(enh_name = names(aa_enh_gene_vicinity),
                                         gene_entID = Vicinity20kbdf_gene$gene_id[aa_enh_gene_vicinity],
                                         stringsAsFactors = F)
aaaaa <- unlist(lapply(lapply(strsplit(Enhancer_gene_20kb_vicinity_mat$enh_name, "_"), "[", c(1,2)), paste, collapse = "_"))
Enhancer_gene_20kb_vicinity_mat$enh_name <- aaaaa

###########################
# getting the same dataframes for gene - gene interactions and vicinity
# use 50kb up and 50kb downstream for vicinity

Promoter1kbGR_gene

Vicinity100kbdf_gene <- as.data.frame(resize(genesAll, width=1, fix='start'))
Vicinity100kbdf_gene$end <- Vicinity100kbdf_gene$start + 50000
Vicinity100kbdf_gene$start <- Vicinity100kbdf_gene$start - 50000
Vicinity100kbdf_gene <- StartEndCoordinateCheck(inputRangedf = Vicinity100kbdf_gene)
Vicinity100kbGR_gene <- makeGRangesFromDataFrame(Vicinity100kbdf_gene)

aa_list <- list()
aa_list[[1]] <- Vicinity100kbdf_gene

gene_100kb_Vicinity <- OverlapInteractionBatchExtract(subjectCoord1=Vicinity100kbdf_gene,
                                                      subjectCoord2List=aa_list, 
                                                      interactingA = MCF7_interaction_hg38_expOnly$new_A, 
                                                      interactingB = MCF7_interaction_hg38_expOnly$new_B,
                                                      mode = "overlap")



names(gene_100kb_Vicinity$OverlapList[[1]]) <- paste0(Vicinity100kbdf_gene$gene_id, "__")
aa_gene_gene_vicinity <- unlist(gene_100kb_Vicinity$OverlapList[[1]])
aa_nu1 <- (unlist(lapply(strsplit(names(aa_gene_gene_vicinity), "__"), "[[", 1)))
gene_gene_vicinity_mat_100kb <- data.frame(gene1_entID = aa_nu1,
                                      gene2_entID = Vicinity100kbdf_gene$gene_id[aa_gene_gene_vicinity],
                                      stringsAsFactors = F)

aa_list <- list()
aa_list[[1]] <- Promoter1kbdf_gene

promoter_promoter_interact_experimental_only <- OverlapInteractionBatchExtract(subjectCoord1=Promoter1kbdf_gene,
                                                                                 subjectCoord2List=aa_list, 
                                                                                 interactingA = MCF7_interaction_hg38_expOnly$new_A, 
                                                                                 interactingB = MCF7_interaction_hg38_expOnly$new_B,
                                                                                 mode = "interaction")
names(promoter_promoter_interact_experimental_only$IntList[[1]]) <- paste0(Promoter1kbdf_gene$gene_id, "__")
aa_gene_gene_ovl <- unlist(promoter_promoter_interact_experimental_only$IntList[[1]])
aa_nu1 <- (unlist(lapply(strsplit(names(aa_gene_gene_ovl), "__"), "[[", 1)))
promoter_promoter_interaction_mat <- data.frame(gene1_entID = aa_nu1,
                                           gene2_entID = Promoter1kbdf_gene$gene_id[aa_gene_gene_ovl],
                                           stringsAsFactors = F)


#####################################################################################################################
  
aa_list <- list()
aa_list[[1]] <- aa_pos_neg

enhancer_enhancer_interact <- OverlapInteractionBatchExtract(subjectCoord1=aa_pos_neg,
                                                                   subjectCoord2List=aa_list, 
                                                                   interactingA = MCF7_interaction_hg38$new_A, 
                                                                   interactingB = MCF7_interaction_hg38$new_B,
                                                                   mode = "both")

names(enhancer_enhancer_interact$OverlapList[[1]]) <- paste0(names(aa_pos_neg), "_")
names(enhancer_enhancer_interact$IntList[[1]]) <- paste0(names(aa_pos_neg), "_")

aa_enh_enh_ov <- unlist(enhancer_enhancer_interact$OverlapList[[1]])
aa_enh_enh_in <- unlist(enhancer_enhancer_interact$IntList[[1]])

Enhancer_Enhancer_overlap_mat <- data.frame(enh_name_1 = names(aa_enh_enh_ov),
                                            enh_name_2 = names(aa_pos_neg)[aa_enh_enh_ov],
                                                stringsAsFactors = F)
aaaaa <- unlist(lapply(lapply(strsplit(Enhancer_Enhancer_overlap_mat$enh_name_1, "_"), "[", c(1,2)), paste, collapse = "_"))
Enhancer_Enhancer_overlap_mat$enh_name_1 <- aaaaa
Enhancer_Enhancer_overlap_mat <- Enhancer_Enhancer_overlap_mat[! Enhancer_Enhancer_overlap_mat$enh_name_1 == Enhancer_Enhancer_overlap_mat$enh_name_2,]

Enhancer_Enhancer_interac_mat <- data.frame(enh_name_1 = names(aa_enh_enh_in),
                                            enh_name_2 = names(aa_pos_neg)[aa_enh_enh_in],
                                                stringsAsFactors = F)
aaaaa <- unlist(lapply(lapply(strsplit(Enhancer_Enhancer_interac_mat$enh_name_1, "_"), "[", c(1,2)), paste, collapse = "_"))
Enhancer_Enhancer_interac_mat$enh_name_1 <- aaaaa
Enhancer_Enhancer_interac_mat <- Enhancer_Enhancer_interac_mat[! Enhancer_Enhancer_interac_mat$enh_name_1 == Enhancer_Enhancer_interac_mat$enh_name_2,]





##############################################################################################################################
#Get the overlap of all enhancers with ER chip peaks
Enhancer.ERchip.Overlap.byEnhancer <- list()#a list where each entry is the overlap list of the enhancers with one chip dataset
for(i in 1:length(ERpickCistromeTimePoint)){
  Enhancer.ERchip.Overlap.byEnhancer[[i]] <- overlapExtractor(subjectCoord1 = MCFEnhancersGR,
                                                              subjectCoord2 =ERpickCistromeTimePoint[[i]] )
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







##################################
###gene
AssociatedChip.gene.20kb <-  ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList , gene.enhancer = Vicinity20kb.Enhancer.plusPromoter.By.gene.Extended ,enhancerChipBinary = Enhancer.ERchip.interORover.byEnhancer.Binary,promoterChipBinary = Promoter.gene.ERchip.interORover.byPromoter.Binary)
AssociatedChip.gene.200kb <- ChipPeakInvestigator(promoter.promoter = Promoter.gene.PromoterIntList , gene.enhancer = Vicinity200kb.Enhancer.plusPromoter.By.gene.Extended ,enhancerChipBinary = Enhancer.ERchip.interORover.byEnhancer.Binary,promoterChipBinary = Promoter.gene.ERchip.interORover.byPromoter.Binary)
###transcript
AssociatedChip.transcript.20kb  <-  ChipPeakInvestigator(promoter.promoter = Promoter.transcrpt.PromoterIntList , gene.enhancer = Vicinity20kb.Enhancer.plusPromoter.By.Transcript.Extended ,enhancerChipBinary = Enhancer.ERchip.interORover.byEnhancer.Binary,promoterChipBinary = Promoter.transcript.ERchip.interORover.byPromoter.Binary)
AssociatedChip.transcript.200kb <-  ChipPeakInvestigator(promoter.promoter = Promoter.transcrpt.PromoterIntList , gene.enhancer = Vicinity200kb.Enhancer.plusPromoter.By.Transcript.Extended ,enhancerChipBinary = Enhancer.ERchip.interORover.byEnhancer.Binary,promoterChipBinary = Promoter.transcript.ERchip.interORover.byPromoter.Binary)

#Filter for genes or transcripts that their association with a chip peak has been repeated in experiments of the same conditions
names(ERpickCistromeTimePoint)
sum(rowSums(AssociatedChip.gene.20kb[[2]][,c(c(1:9),16)]) > 7)
sum(rowSums(AssociatedChip.gene.20kb[[2]][,c(10,11,17,18)]) > 3)
sum(rowSums(AssociatedChip.gene.20kb[[2]][,12:15]) > 3)
length(unique(which(rowSums(AssociatedChip.gene.20kb[[2]][,c(c(2:6),c(8,9),16)]) >5),which(rowSums(AssociatedChip.gene.20kb[[2]][,c(11,17,18)]) > 1 ),which(rowSums(AssociatedChip.gene.20kb[[2]][,c(12,15)]) > 0 )))
# sum(rowSums(GeneERchipVicinityPlusInteractionBinary[,c(10,11,17,18)]) > 1 )
# sum(rowSums(GeneERchipVicinityPlusInteractionBinary[,12:15]) > 2 )
-c(1,7,13,14)



############################################################################################################
aamin_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
aamaxL_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
for(i in 1:length(TF.motifs.Expanded_LLR_to_pval_pseudo)){
  aamaxL_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][1,1] / 1000
  aa <- which(TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][,2] > 0.0001)[1]
  if(aa == 1){
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][aa, 1] / 1000
  }else{
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][(aa-1), 1] / 1000
  }
}
names(aamin_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
names(aamaxL_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
aathr <- 1- (aamin_LLR/aamaxL_LLR)
aa_names <- c("ESR1_2","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","RARA", "RELA","RUNX1","SP1","TFAP2C","YBX1")
TF.motifs.Expanded_pseudo_exp12 <- TF.motifs.Expanded_pseudo[aa_names]
TF.motifs.Expanded_pseudo_exp12_t <- lapply(TF.motifs.Expanded_pseudo_exp12, t)

aa_pos_dif_filtered_all2 <- aa_pos_dif_filtered_all_acc
aa_neg_dif_filtered_all2 <- aa_neg_dif_filtered_all_acc
aa_pos_dif_filtered_all2 <- aa_pos_dif_filtered_all2[rowSums(abs(aa_pos_dif_filtered_all2) >= 0.1, na.rm = T) > 0,]
aa_neg_dif_filtered_all2 <- aa_neg_dif_filtered_all2[rowSums(abs(aa_neg_dif_filtered_all2) >= 0.1, na.rm = T) > 0,]
aa_all_names <- c(rownames(aa_pos_dif_filtered_all2 ), rownames(aa_neg_dif_filtered_all2))
aa_all_names <- unlist(lapply(strsplit(aa_all_names, split = "_"), "[[", 2))
aa_all_qtl <- ER_gtex_eqtl_sig[as.numeric(aa_all_names),]
aa_all_qtl_uniq <- unique(aa_all_qtl$variant_id)
aa_eqtl <-  strsplit(aa_all_qtl_uniq, split= "_")
aa_eqtl_chr <- unlist(lapply(aa_eqtl, "[[", 1))
aa_eqtl_st <- as.numeric(unlist(lapply(aa_eqtl, "[[", 2)))
aa_eqtl_ref <- unlist(lapply(aa_eqtl, "[[", 3))
aa_eqtl_alt <- unlist(lapply(aa_eqtl, "[[", 4))

aamin_LLR_low <- aamin_LLR - 0.1*aamin_LLR
aamin_LLR_lowest <- numeric(length(aamin_LLR))
names(aamin_LLR_lowest) <- names(aamin_LLR_low)
aa_TF.motifs.Expanded_pseudo_exp12_t_ps <- lapply(TF.motifs.Expanded_pseudo_exp12_t, Addpsudo)
aa_tf_res <- list()
for(i in 1:length(aa_eqtl_alt)){
  print(i)
  aa_tf_res[[i]] <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                                     my_motifs = aa_TF.motifs.Expanded_pseudo_exp12_t_ps,
                                     snp_chr = aa_eqtl_chr[i],
                                     snp_start = aa_eqtl_st[i],
                                     snp_ref = aa_eqtl_ref[i],
                                     snp_alt = aa_eqtl_alt[i],
                                     min_LLR = aamin_LLR_low[aa_names])
  
}
names(aa_tf_res) <- aa_all_qtl_uniq
aarw <- unlist(lapply(aa_tf_res, nrow))
table(aarw)
# aarw
# 0  1  2  3  4 
# 24 47  7  1  1 
nrow(aa_tf_res[[5]])

aa_tf_res[[2]]

aa_tf_res[[1]]
nrow(aa_tf_res[[4]])
aa_all_qtl_uniq[4]
which(ER_gtex_eqtl_sig$variant_id %in% aa_all_qtl_uniq[4])
c("eqtl_1387997", "eqtl_1388480" ,"eqtl_1388870", "eqtl_1388924") %in% rownames(aa_pos_dif_filtered_all2)


nrow(aa_pos_dif_filtered$pos_185)
# > table(aarw)
# aarw
# 0  1  2  3 
# 40 34  4  2 


# aarw 0.05
# 0  1  2  3 
# 36 28 14  2 

aa_tf_res_zero <- list()
for(i in 1:length(aa_eqtl_alt)){
  print(i)
  aa_tf_res_zero[[i]] <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                                     my_motifs = aa_TF.motifs.Expanded_pseudo_exp12_t_ps,
                                     snp_chr = aa_eqtl_chr[i],
                                     snp_start = aa_eqtl_st[i],
                                     snp_ref = aa_eqtl_ref[i],
                                     snp_alt = aa_eqtl_alt[i],
                                     min_LLR = aamin_LLR_lowest[aa_names])
  
}
names(aa_tf_res_zero) <- aa_all_qtl_uniq
aarw <- unlist(lapply(aa_tf_res_zero, nrow))
table(aarw)

aa_tf_res[[17]]
aa_tf_res_zero[[17]]
aa_tf_res_zero$chr17_75550062_G_A_b38
aa_tf_res_zero$chr17_75549684_T_C_b38


aa_tat1 <- "ACACGCCCGCCACGTGAGTGCCTGGGCTCCCACCAGTCAGGTCCGCGCGACCCGGTCCCCAT"
aa_tat2 <- "ACACGCCCGCCACGTGAGTGCCTGGGCTCCCGCCAGTCAGGTCCGCGCGACCCGGTCCCCAT"
aa_tat1_ms <- MotifScore(seq = aa_tat1, 
                         motifList = aa_TF.motifs.Expanded_pseudo_exp12_t_ps,
                         bg = c(0.25,0.25,0.25,0.25))
aa_tat2_ms <- MotifScore(seq = aa_tat2, 
                         motifList = aa_TF.motifs.Expanded_pseudo_exp12_t_ps,
                         bg = c(0.25,0.25,0.25,0.25))

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

# run GEMSTAT exp16 models on all common snps

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_snp <- readPeakFile(aa_snp_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
start(aa_snp) <- start(aa_snp) - 1

library(BSgenome.Hsapiens.UCSC.hg38)
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_16")
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

aa_snpdf <- as.data.frame(aa_snp)
aa_ref <- as.character(levels(aa_snpdf$V8)[as.numeric(aa_snpdf$V8)])
aa_alt <- as.character(levels(aa_snpdf$V9)[as.numeric(aa_snpdf$V9)])

aa_ref_Alt <- cbind(aa_ref, aa_alt)
aa_names <-  as.character(levels(aa_snpdf$V10)[as.numeric(aa_snpdf$V10)])
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_16")
aa_pos_neg_vars_snp <- write_variant_seq(enhancer_GR = aa_pos_neg, 
                                         eqtl_GR = aa_snp,
                                         eqtl_names = aa_names,
                                         eqtl_ref_alt = aa_ref_Alt,
                                         all_combs = F,
                                         my_genome = BSgenome.Hsapiens.UCSC.hg38, 
                                         label = c(rep(1, length(aa_pos)), rep(0, length(aa_neg))))

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")

# write GEMSTAT jobs to run
# write hal job for KDs
aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[16]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[16]]$annot)[aas1[1:25]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[16]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[16]]$annot)[aas2[1:25]]
aaadd3 <- union(aaadd1, aaadd2)

aa_par_name <- paste0(aaadd3, ".txt.Filter")
aa_seq_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_16/Variant_seq/",
                           recursive = T)
aa_lab_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_16/Variant_Labels/",
                           recursive = T)
aa_nam_spl <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_seq_files, 
                                                                   split = "\\."), 
                                                          "[[", 1)), split = "_"),
                                   "[", c(3,4)),
                            paste, collapse = "_"))


setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_16")
for(i in 1:length(aa_par_name)){
  for(j in 1:length(aa_seq_files)){
    cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr ",
          paste0("-s Variant_seq/", aa_seq_files[j]), 
          paste0("-e Variant_Labels/", aa_lab_files[j]),
          "-m motifs.wtmx -f TF_exp.tab", 
          paste0("-fo Varaint_out/", aa_nam_spl[j], "_", aaadd3[i], ".out"), 
          "-o DIRECT -c Coop/coop.par ", 
          paste0("-p Trained_par/", aa_par_name[i]),
          "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1\n"
    ), 
    sep = " ", append = (i!=1), file = "variant_jobs_snp_exp16.job")
  }
}

################# ################# ################# ################# #################
# read results of the snps
# read output and find the ones with some effect
# first read WT for all models

aa_wt_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_16/Varaint_out/",
                              pattern = "NA_NA_par*", full.names = T)
aa_wt_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_16/Varaint_out/",
                          pattern = "NA_NA_par*", full.names = F)
aa_nam_spl <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_wt_files, split = "\\."), "[[", 1)), split = "_"), "[", c(3,4)), paste, collapse = "_"))

aa_t <- read_output_train_test_GEMSTAT_indiv(output_file = aa_wt_files_ful[1])
aa_wt_models_snp <- matrix(nrow=nrow(aa_t$GT),
                           ncol=length(aa_wt_files))

colnames(aa_wt_models_snp) <- aa_nam_spl
rownames(aa_wt_models_snp) <- aa_t$pred$Rows
for(i in 1:ncol(aa_wt_models_snp)){
  aa_wt_models_snp[, i] <- read_output_train_test_GEMSTAT_indiv(output_file = aa_wt_files_ful[i])$pred$X1
}
aa_neg_out_snp <- list()
aa_pos_out_snp <- list()

aa_neg_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_16/Varaint_out/",
                               pattern = "neg_*", full.names = T)
aa_neg_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_16/Varaint_out/",
                           pattern = "neg_*", full.names = F)
aa_neg_sp_1 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_neg_files, split = "\\."), "[[", 1)), split = "_"), "[", c(1,2)), paste, collapse = "_"))
aa_neg_sp_2 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_neg_files, split = "\\."), "[[", 1)), split = "_"), "[", c(3,4)), paste, collapse = "_"))
aa_neg_uniq <- unique(aa_neg_sp_1)

aa_pos_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_16/Varaint_out/",
                               pattern = "pos_*", full.names = T)
aa_pos_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_16/Varaint_out/",
                           pattern = "pos_*", full.names = F)
aa_pos_sp_1 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_pos_files, split = "\\."), "[[", 1)), split = "_"), "[", c(1,2)), paste, collapse = "_"))
aa_pos_sp_2 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_pos_files, split = "\\."), "[[", 1)), split = "_"), "[", c(3,4)), paste, collapse = "_"))
aa_pos_uniq <- unique(aa_pos_sp_1)

aa_neg_dif_snp <- list()
aa_pos_dif_snp <- list()

aa_neg_dif_snp_percentile <- list()
aa_pos_dif_snp_percentile <- list()

aa_Exp16_WT_model_prediction_Filtered_GT <- rownames(Exp16_WT_model_prediction_Filtered)
aa_Exp16_WT_model_prediction_Filtered_GT <- unlist(lapply(strsplit(aa_Exp16_WT_model_prediction_Filtered_GT, split = "_"), "[[", 1))
aa_Exp16_WT_model_prediction_Filtered_GT[aa_Exp16_WT_model_prediction_Filtered_GT == "pos"] <- 1
aa_Exp16_WT_model_prediction_Filtered_GT[aa_Exp16_WT_model_prediction_Filtered_GT == "neg"] <- 0
aa_Exp16_WT_model_prediction_Filtered_GT <- as.numeric(aa_Exp16_WT_model_prediction_Filtered_GT)
names(aa_Exp16_WT_model_prediction_Filtered_GT) <- rownames(Exp16_WT_model_prediction_Filtered)


for(i in 1:length(aa_neg_uniq)){
  print(i)
  aa_cur_ind <- which(aa_neg_sp_1 %in% aa_neg_uniq[i])
  aa_t <- read.table(file = aa_neg_files_ful[aa_cur_ind[1]], header = T, stringsAsFactors = F)
  aa_ind <- seq(2, nrow(aa_t), 2)
  aa_t <- aa_t[aa_ind, ]
  # aa_neg_out_snp[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  # colnames(aa_neg_out_snp[[i]]) <- colnames(aa_wt_models_snp)
  # rownames(aa_neg_out_snp[[i]]) <- aa_t$Rows
  aa_neg_dif_snp[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  colnames(aa_neg_dif_snp[[i]]) <- colnames(aa_wt_models_snp)
  rownames(aa_neg_dif_snp[[i]]) <- aa_t$Rows
  aa_neg_dif_snp_percentile[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  colnames(aa_neg_dif_snp_percentile[[i]]) <- colnames(aa_wt_models_snp)
  rownames(aa_neg_dif_snp_percentile[[i]]) <- aa_t$Rows
  # for(j in 1:length(aa_cur_ind)){
  #   aa_table <- read.table(file =aa_neg_files_ful[aa_cur_ind[j]],
  #                          header = T, stringsAsFactors = F)
  #   aa_neg_out_snp[[i]][, aa_neg_sp_2[aa_cur_ind[j]]] <- aa_table[seq(2, nrow(aa_table), 2),"X1"]
  # }
  for(j in 1:nrow(aa_neg_out_snp[[i]])){
    aa_neg_dif_snp[[i]][j, ] <- (aa_neg_out_snp[[i]][j, ] - Exp16_WT_model_prediction_Filtered[aa_neg_uniq[i], ])/ Exp16_WT_model_prediction_Filtered[aa_neg_uniq[i], ]
    for(kk in 1:ncol(aa_neg_out_snp[[i]])){
      if(!is.na(Exp16_WT_model_prediction_Filtered[aa_neg_uniq[i],kk])){
        aa_pos_perc_before <- sum(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 1,kk] <= Exp16_WT_model_prediction_Filtered[aa_neg_uniq[i],kk], na.rm = T)/sum(!is.na(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 1,kk]))
        aa_neg_perc_before <- sum(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 0,kk] >  Exp16_WT_model_prediction_Filtered[aa_neg_uniq[i],kk], na.rm = T)/sum(!is.na(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 0,kk]))
        aa_pos_perc_after <- sum(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 1,kk] <=  aa_neg_out_snp[[i]][j, kk], na.rm = T)/sum(!is.na(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 1,kk]))
        aa_neg_perc_after <- sum(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 0,kk] >   aa_neg_out_snp[[i]][j, kk], na.rm = T)/sum(!is.na(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 0,kk]))
        aa_neg_dif_snp_percentile[[i]][j,kk] <- -(aa_neg_perc_after - aa_neg_perc_before) + (aa_pos_perc_after - aa_pos_perc_before)
        # print("aa_pos_perc_before")
        # print(aa_pos_perc_before)
        # print("aa_neg_perc_before")
        # print(aa_neg_perc_before)
        # print("aa_pos_perc_after")
        # print(aa_pos_perc_after)
        # print("aa_neg_perc_after")
        # print(aa_neg_perc_after)
        
        
        
        
      }
    }
  }
}
names(aa_neg_out_snp) <- aa_neg_uniq
names(aa_neg_dif_snp) <- aa_neg_uniq
names(aa_neg_dif_snp_percentile) <- aa_neg_uniq
aa_neg_out_snp_all <- do.call(what = rbind, aa_neg_out_snp)
aa_neg_dif_all_snp_percentile <- do.call(what = rbind, aa_neg_dif_snp_percentile)


for(i in 1:length(aa_pos_uniq)){
  print(i)
  aa_cur_ind <- which(aa_pos_sp_1 %in% aa_pos_uniq[i])
  aa_t <- read.table(file = aa_pos_files_ful[aa_cur_ind[1]], header = T, stringsAsFactors = F)
  aa_ind <- seq(2, nrow(aa_t), 2)
  aa_t <- aa_t[aa_ind, ]
  # aa_pos_out_snp[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  # colnames(aa_pos_out_snp[[i]]) <- colnames(aa_wt_models_snp)
  # rownames(aa_pos_out_snp[[i]]) <- aa_t$Rows
  aa_pos_dif_snp[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  colnames(aa_pos_dif_snp[[i]]) <- colnames(aa_wt_models_snp)
  rownames(aa_pos_dif_snp[[i]]) <- aa_t$Rows
  aa_pos_dif_snp_percentile[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  colnames(aa_pos_dif_snp_percentile[[i]]) <- colnames(aa_wt_models_snp)
  rownames(aa_pos_dif_snp_percentile[[i]]) <- aa_t$Rows
  
  # for(j in 1:length(aa_cur_ind)){
  #   aa_table <- read.table(file =aa_pos_files_ful[aa_cur_ind[j]],
  #                          header = T, stringsAsFactors = F)
  #   aa_pos_out_snp[[i]][, aa_pos_sp_2[aa_cur_ind[j]]] <- aa_table[seq(2, nrow(aa_table), 2),"X1"]
  # }
  for(j in 1:nrow(aa_pos_out_snp[[i]])){
    aa_pos_dif_snp[[i]][j, ] <- (aa_pos_out_snp[[i]][j, ] - Exp16_WT_model_prediction_Filtered[aa_pos_uniq[i], ])/ Exp16_WT_model_prediction_Filtered[aa_pos_uniq[i], ]
    for(kk in 1:ncol(aa_pos_out_snp[[i]])){
      if(!is.na(Exp16_WT_model_prediction_Filtered[aa_pos_uniq[i],kk])){
        aa_pos_perc_before <- sum(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 1,kk] <= Exp16_WT_model_prediction_Filtered[aa_neg_uniq[i],kk], na.rm = T)/sum(!is.na(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 1,kk]))
        aa_neg_perc_before <- sum(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 0,kk] >  Exp16_WT_model_prediction_Filtered[aa_neg_uniq[i],kk], na.rm = T)/sum(!is.na(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 0,kk]))
        aa_pos_perc_after <- sum(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 1,kk] <=  aa_pos_out_snp[[i]][j, kk], na.rm = T)/sum(!is.na(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 1,kk]))
        aa_neg_perc_after <- sum(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 0,kk] >   aa_pos_out_snp[[i]][j, kk], na.rm = T)/sum(!is.na(Exp16_WT_model_prediction_Filtered[aa_Exp16_WT_model_prediction_Filtered_GT == 0,kk]))
        aa_pos_dif_snp_percentile[[i]][j,kk] <- -(aa_neg_perc_after - aa_neg_perc_before) + (aa_pos_perc_after - aa_pos_perc_before)
        # print("aa_pos_perc_before")
        # print(aa_pos_perc_before)
        # print("aa_neg_perc_before")
        # print(aa_neg_perc_before)
        # print("aa_pos_perc_after")
        # print(aa_pos_perc_after)
        # print("aa_neg_perc_after")
        # print(aa_neg_perc_after)
        # 
        # 
        # 
        
      }
    }
  }
  
}
names(aa_pos_out_snp) <- aa_pos_uniq
names(aa_pos_dif_snp) <- aa_pos_uniq
names(aa_pos_dif_snp_percentile) <- aa_pos_uniq


aa_pos_out_snp_all <- do.call(what = rbind, aa_pos_out_snp)

# aa_pos_dif_snp_old_noFilt <- aa_pos_dif_snp # before filtering the outputs of WT models
# aa_neg_dif_snp_old_noFilt <- aa_neg_dif_snp
aa_pos_dif_all_snp <- do.call(what = rbind, aa_pos_dif_snp)
aa_pos_dif_all_snp_percentile <- do.call(what = rbind, aa_pos_dif_snp_percentile)

par(mfrow = c(1,1), mar = c(4,4,4,4))
boxplot.matrix(t(aa_pos_dif_all_snp[1:5,])*100)

aa_neg_dif_all_snp <- do.call(what = rbind, aa_neg_dif_snp)
#boxplot.matrix(t(aa_neg_dif_all_snp))
sum(colSums(is.na(aa_pos_dif_all_snp)) == nrow(aa_pos_dif_all_snp))

aa_pos_neg_out_snp_all <- rbind(aa_pos_out_snp_all, aa_neg_out_snp_all)

# z <- zClust(x=t(aa_pos_dif_all_snp[,-c(17)]), scale="row", zlim=c(-3,3), method="average")
# # require(RColorBrewer)
# cols <- colorRampPalette(brewer.pal(10, "RdBu"))(9)
# heatmap.2(z$data, trace='none',
#           col=rev(cols),
#           #Rowv=F,Colv=F,
#           main="predicted change"
#           ,Rowv=z$Rowv, 
#           Colv=z$Colv
# )
# cols <- colorRampPalette(brewer.pal(10, "RdBu"))(9)
# heatmap.2(x = t(aa_pos_dif_all[,-c(17)]), 
#           trace='none',
#           Rowv = T,
#           Colv = T,
#           breaks = c(-1,-0.5,-0.1, 0.1,0.5, 1, 1.5, 2, 2.5, 3), col = rev(cols))
# z <- zClust(x=t(aa_neg_dif_all[,-c(17)]), scale="row", zlim=c(-3,3), method="average")
# # require(RColorBrewer)
# cols <- colorRampPalette(brewer.pal(10, "RdBu"))(9)
# heatmap.2(z$data, trace='none',
#           col=rev(cols),
#           Rowv=F,Colv=F,
#           main="predicted change neg"
#           # ,Rowv=z$Rowv, 
#           # Colv=z$Colv
# )
# find enhancers containing multiple eqtls where some don't affect expression and some do


aa_pos_dif_filtered_snp <- aa_pos_dif_snp
for(i in length(aa_pos_dif_filtered_snp):1){
  if(max(abs(aa_pos_dif_filtered_snp[[i]]), na.rm = T) < 0.05){
    aa_pos_dif_filtered_snp[[i]] <- NULL
  }
}
aanrwo <- unlist(lapply(aa_pos_dif_filtered_snp, nrow))
par(mar = c(1,1,1,1), mfrow = c(10, 10))
for(i in 1:length(aa_pos_dif_filtered_snp)){
  boxplot.matrix(t(aa_pos_dif_filtered_snp[[i]]), 
                 las= 2, 
                 main = names(aa_pos_dif_filtered_snp)[i],
                 xaxt = "n", outline = F)
  abline(h= seq(-3, 3, 0.05), col = 2, lwd = 0.5, lty = 3)
  abline(h= 0, col = 2, lwd = 1, lty = 1)
}
# boxplot.matrix(t(aa_pos_dif_filtered_snp$pos_312), 
#                las= 2, main = "pos_312 predicted change")
# View(aa_pos_dif_filtered_snp$pos_171)

aa_neg_dif_filtered_snp <- aa_neg_dif_snp
for(i in length(aa_neg_dif_filtered_snp):1){
  if(max(abs(aa_neg_dif_filtered_snp[[i]]), na.rm = T) < 0.05){
    aa_neg_dif_filtered_snp[[i]] <- NULL
  }
}
aanrwo <- unlist(lapply(aa_neg_dif_filtered_snp, nrow))
par(mar = c(1,1,1,1), mfrow = c(10, 7))
for(i in 1:length(aa_neg_dif_filtered_snp)){
  if(aanrwo[i] > 1){
    boxplot.matrix(t(aa_neg_dif_filtered_snp[[i]]), 
                   las= 2, 
                   main = names(aa_neg_dif_filtered_snp)[i],
                   xaxt = "n")
    abline(h= seq(-3, 3, 0.05), col = 2, lwd = 0.5, lty = 3)
    abline(h= 0, col = 2, lwd = 1, lty = 1)
  }
  
}


################# ################# ################# ################# #################
# check if eqtls are enriched in the affected ones
head(ER_gtex_eqtl_sig)

aa_eqtl_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/sig_eqtl.bed"
aa_qtl_gtex <- readPeakFile(aa_eqtl_Add, as = "GRanges")
start(aa_qtl_gtex) <- start(aa_qtl_gtex) -1
end(aa_qtl_gtex) <- end(aa_qtl_gtex) + 1
aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_snp <- readPeakFile(aa_snp_Add, as = "GRanges")
start(aa_snp) <- start(aa_snp) - 1
end(aa_snp) <- end(aa_snp) + 1
######## overlap with gtex eqtl
aa_overlap_holder <- findOverlaps(query = aa_snp, subject = aa_qtl_gtex)
#ol_enhancers_index <- unique(overlap_holder@from)
aa_snp_qtl_marker_gtex <- numeric(length(aa_snp))
aa_snp_qtl_marker_gtex[aa_overlap_holder@from] <- 1
aa_snp$V10 <- as.character(levels(aa_snp$V10)[as.numeric(aa_snp$V10)])
names(aa_snp_qtl_marker_gtex) <- aa_snp$V10

######## overlap with pancan eqtl

ranges(ER_pancan_eqtl_cis_gr38) <- ranges(ER_pancan_eqtl_cis_gr38) + 1
aa_overlap_holder <- findOverlaps(query = aa_snp, subject = ER_pancan_eqtl_cis_gr38)
aa_snp_qtl_marker_pancan_cis <- numeric(length(aa_snp))
aa_snp_qtl_marker_pancan_cis[aa_overlap_holder@from] <- 1
#aa_snp$V10 <- as.character(levels(aa_snp$V10)[as.numeric(aa_snp$V10)])
names(aa_snp_qtl_marker_pancan_cis) <- aa_snp$V10

ranges(ER_pancan_eqtl_trans_gr38) <- ranges(ER_pancan_eqtl_trans_gr38) + 1
aa_overlap_holder <- findOverlaps(query = aa_snp,
                                  subject = ER_pancan_eqtl_trans_gr38)
aa_snp_qtl_marker_pancan_trans <- numeric(length(aa_snp))
aa_snp_qtl_marker_pancan_trans[aa_overlap_holder@from] <- 1
names(aa_snp_qtl_marker_pancan_trans) <- aa_snp$V10

aa_snp_qtl_marker_pancan <- aa_snp_qtl_marker_pancan_trans + aa_snp_qtl_marker_pancan_cis
aa_snp_qtl_marker_gtex_pancan <- aa_snp_qtl_marker_pancan + aa_snp_qtl_marker_gtex
aa_snp_qtl_marker_gtex_pancan[aa_snp_qtl_marker_gtex_pancan > 1] <- 1


ER_pancan_cis_trans_gr38 <- c(ER_pancan_eqtl_cis_gr38, ER_pancan_eqtl_trans_gr38)
ranges(ER_pancan_cis_trans_gr38) <- ranges(ER_pancan_cis_trans_gr38) - 1

options(scipen=999)
write.table(ER_pancan_cis_trans_gr38, 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/ER_pancan_cis_Trans.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
options(scipen=0)
ER_pancan_cis_trans_gr38_unique <- ER_pancan_cis_trans_gr38[!duplicated(ER_pancan_cis_trans_gr38$snp_name)]
write.table(ER_pancan_cis_trans_gr38_unique, 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/ER_pancan_cis_Trans_unique.bed",
            quote=F, sep="\t", row.names=F, col.names=F)



aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_eqtl_pancan_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/ER_pancan_cis_Trans.bed"


aaa_eqtl_pancan_pos4 <- bedtools_intersect(bedfile_names = c(aa_eqtl_pancan_Add, aa_pos_Add), 
                                    wa=F, wb=F, loj=F, wo=F, wao=F,
                                    u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                    S=F, sorted=F, merge_after_distance=0, 
                                    output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/eqtl_pancan_pos_4.bed", 
                                    read_output = T)
aaa_eqtl_pancan_neg4 <- bedtools_intersect(bedfile_names = c(aa_eqtl_pancan_Add, aa_neg_Add),
                                    wa=F, wb=F, loj=F, wo=F, wao=F,
                                    u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                    S=F, sorted=F, merge_after_distance=0, 
                                    output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/eqtl_pancan_neg_4.bed",
                                    read_output = T)
start(aaa_eqtl_pancan_pos4) <- start(aaa_eqtl_pancan_pos4) - 1
start(aaa_eqtl_pancan_neg4) <- start(aaa_eqtl_pancan_neg4) - 1
ER_pancan_cis_trans_gr38_pos_neg_4 <-  c(aaa_eqtl_pancan_pos4, aaa_eqtl_pancan_neg4)

######## overlap with GWAS
ranges(GWAS_breast_Cancer_gr38) <- ranges(GWAS_breast_Cancer_gr38) + 1
aa_overlap_holder <- findOverlaps(query = aa_snp,
                                  subject = GWAS_breast_Cancer_gr38)
aa_snp_gwas_marker <- numeric(length(aa_snp))
aa_snp_gwas_marker[aa_overlap_holder@from] <- 1
names(aa_snp_gwas_marker) <-  aa_snp$V10


aa_all_dif_snp <- rbind(aa_pos_dif_all_snp, aa_neg_dif_all_snp)
aa_all_dif_snp_names <- unlist(lapply(strsplit(rownames(aa_all_dif_snp), split="_"), "[[", 2))


aa_snp_gwas_marker_match <-  aa_snp_gwas_marker[match(aa_all_dif_snp_names, names(aa_snp_gwas_marker))]
aa_snp_qtl_marker_pancan_match <- aa_snp_qtl_marker_pancan[match(aa_all_dif_snp_names, names(aa_snp_qtl_marker_pancan))]
aa_snp_qtl_marker_gtex_match <- aa_snp_qtl_marker_gtex[match(aa_all_dif_snp_names, names(aa_snp_qtl_marker_gtex))]
aa_snp_qtl_marker_gtex_pancan_match <- aa_snp_qtl_marker_gtex_pancan[match(aa_all_dif_snp_names, names(aa_snp_qtl_marker_gtex_pancan))]


aa_all_sig_nonsig <- matrix(nrow = 10, ncol = ncol(aa_pos_dif_all_snp))
colnames(aa_all_sig_nonsig) <- colnames(aa_neg_dif_all_snp)[1:(ncol(aa_neg_dif_all_snp))]

rownames(aa_all_sig_nonsig) <- c("sig_all", "sig_qtl_gtex", "nonsig_all", "nonsig_qtl", "pval_pancan", "pval_gtex", "pval_both", "gwas_sig", "sig_qtl_pancan", "sig_qtl_either")



aa_thresh <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
aa_all_sig_hold <- list()

for(i in 1:length(aa_thresh)){
  aa_all_sig_hold[[i]] <- aa_all_sig_nonsig
  for(j in 1:(ncol(aa_all_dif_snp))){
    aa_all_sig_hold[[i]][1,j] <- sum(abs(aa_all_dif_snp[,j]) >= aa_thresh[i], na.rm = T)
    aa_all_sig_hold[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold[[i]][3,j] <- sum(abs(aa_all_dif_snp[,j]) < aa_thresh[i], na.rm = T)
    aa_all_sig_hold[[i]][4,j] <- sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp[,j]) < aa_thresh[i]], na.rm = T)
    aa_all_sig_hold[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp[,j]) >= aa_thresh[i]], na.rm = T)
    
    if(sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold[[i]][5,j] <- phyper(q = sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp[,j]) >= aa_thresh[i]], na.rm = T),
                                          m = sum(aa_snp_qtl_marker_pancan_match == 1),
                                          n = sum(aa_snp_qtl_marker_pancan_match == 0), 
                                          k = aa_all_sig_hold[[i]][1,j],
                                          lower.tail = F)
    }

    if(sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold[[i]][6,j] <- phyper(q = sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp[,j]) >= aa_thresh[i]], na.rm = T),
                                          m = sum(aa_snp_qtl_marker_gtex_match == 1),
                                          n = sum(aa_snp_qtl_marker_gtex_match == 0), 
                                          k = aa_all_sig_hold[[i]][1,j],
                                          lower.tail = F)
    }

    if(sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold[[i]][7,j] <- phyper(q = sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp[,j]) >= aa_thresh[i]], na.rm = T),
                                          m = sum(aa_snp_qtl_marker_gtex_pancan_match == 1),
                                          n = sum(aa_snp_qtl_marker_gtex_pancan_match == 0), 
                                          k = aa_all_sig_hold[[i]][1,j],
                                          lower.tail = F)
    }

    
    aa_all_sig_hold[[i]][8,j] <- sum(aa_snp_gwas_marker_match[abs(aa_all_dif_snp[,j]) >= aa_thresh[i]], na.rm = T)
    
  }
}
names(aa_all_sig_hold) <- aa_thresh

for(i in 1:length(aa_all_sig_hold)){
  print(names(aa_all_sig_hold)[i])
  print("gtex")
  print(sum((aa_all_sig_hold[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold[[i]][8,] > 0, na.rm=T))
}
View(aa_all_sig_hold$`0.9`)



par(mfrow = c(1,1), mar = c(4,4,4,4))

E_RNA_GEMSTAT_Ensemble_Parlist[[16]]$logistic[1:5,]
# get the param list of the ones passed the filteringss
aa_bind <- E_RNA_GEMSTAT_Ensemble_Parlist[[16]]$binding[paste0(colnames(aa_neg_dif_all_snp), "_1"),]
aa_aplh <- E_RNA_GEMSTAT_Ensemble_Parlist[[16]]$alpha[paste0(colnames(aa_neg_dif_all_snp), "_1"),]
aa_coop <- E_RNA_GEMSTAT_Ensemble_Parlist[[16]]$coop[paste0(colnames(aa_neg_dif_all_snp), "_1"),]
aa_qbtm <- E_RNA_GEMSTAT_Ensemble_Parlist[[16]]$qbtm[paste0(colnames(aa_neg_dif_all_snp), "_1")]
aa_beta <- E_RNA_GEMSTAT_Ensemble_Parlist[[16]]$beta[paste0(colnames(aa_neg_dif_all_snp), "_1")]
aa_logi <- E_RNA_GEMSTAT_Ensemble_Parlist[[16]]$logistic[paste0(colnames(aa_neg_dif_all_snp), "_1"),]

aa_all_param <- cbind(log10(aa_bind), log10(aa_aplh), log10(aa_coop), aa_qbtm, aa_beta, aa_logi)
colnames(aa_all_param)[1:17] <- paste0(colnames(aa_all_param)[1:17], "_bind")
colnames(aa_all_param)[18:34] <- paste0(colnames(aa_all_param)[18:34], "_alpha")


aa_pca <- prcomp(aa_all_param, center = TRUE, scale = TRUE)
summary(aa_pca)

par(mfrow = c(1,1), mar = c(4,4,4,4))
plot(aa_pca$x[,1],aa_pca$x[,2], xlab="PC1 (18%)",
     ylab = "PC2 (12%)", main = "PC1 / PC2 - plot"
     , col =aacol[round(aacol_mat["0.9",]) + 3] 
     )
library(RColorBrewer)
aacol <- colorRampPalette(brewer.pal(9, "YlOrRd"))(7)
aacol_mat <- matrix(nrow = length(aa_all_sig_hold), 
                    ncol = ncol(aa_all_sig_hold[[1]]))
rownames(aacol_mat) <- names(aa_all_sig_hold)
colnames(aacol_mat) <- colnames(aa_all_sig_hold[[1]])

for(i in 1:nrow(aacol_mat)){
#  for(j in 1:ncol(aacol_mat)){
    aacol_mat[i,] <- -log10(aa_all_sig_hold[[i]][6,])
#  }
}
par(mfrow = c(5,3), mar = c(2,2,2,2))
for(i in 1:nrow(aacol_mat)){
  plot(aa_pca$x[,1],aa_pca$x[,2], xlab="PC1 (18%)", pch = 16, cex = 1.5,
       ylab = "PC2 (12%)", main = rownames(aacol_mat)[i]
       , col =aacol[round(aacol_mat[i,]) + 4] 
  )
}


library(ggbiplot)
par(mfrow = c(1,1), mar = c(8,4,4,4))
ggbiplot(aa_pca, labels = rownames(aa_pca$x))

par(mfrow = c(3,1), mar = c(8,4,4,4))
barplot(aa_all_param["par_5966_1",], las = 2, main = "par_5966")
barplot(aa_all_param["par_3191_1",], las = 2, main = "par_3191")
barplot(aa_all_param["par_296_1",], las = 2, main = "par_296")


#####################

aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_snp <- readPeakFile(aa_snp_Add, as = "data.frame")
aa_snp$V10 <- as.character(levels(aa_snp$V10)[as.numeric(aa_snp$V10)])
aa_snp$V9 <- as.character(levels(aa_snp$V9)[as.numeric(aa_snp$V9)])
aa_snp$V8 <- as.character(levels(aa_snp$V8)[as.numeric(aa_snp$V8)])
aa_snp$V1 <- as.character(levels(aa_snp$V1)[as.numeric(aa_snp$V1)])
aa_snp$V2 <- aa_snp$V2 - 1


aaa <- rownames(aa_all_dif_snp)[which(aa_all_dif_snp[,"par_5966"] > 0.9)]

aa_1 <- unlist(lapply(strsplit(aaa, split="_"), "[[", 2))
aa_2 <- unlist(lapply(strsplit(aaa, split="_"), "[[", 3))

aa_eqtl_chr <- character(length(aaa))
aa_eqtl_st <- numeric(length(aaa))
aa_eqtl_ref <- numeric(length(aaa))
aa_eqtl_alt <- numeric(length(aaa))

for(i in 1:length(aaa)){
  aa_eqtl_chr[i] <- aa_snp$V1[aa_snp$V10 == aa_1[i]]
  aa_eqtl_st[i] <- aa_snp$V2[aa_snp$V10 == aa_1[i]]
  aa_eqtl_ref[i] <- aa_snp$V8[aa_snp$V10 == aa_1[i]]
  aaaaalt <- unlist(strsplit(aa_snp$V9[aa_snp$V10 == aa_1[i]], split = ","))
  aa_eqtl_alt[i] <- aaaaalt[as.numeric(aa_2[i])]
}



aamin_LLR_low <- aamin_LLR - 0.1*aamin_LLR
aamin_LLR_lowest <- numeric(length(aamin_LLR))
names(aamin_LLR_lowest) <- names(aamin_LLR_low)
aa_TF.motifs.Expanded_pseudo_exp12_t_ps <- lapply(TF.motifs.Expanded_pseudo_exp12_t, Addpsudo)
aa_snp_inves <- list()
for(i in 1:length(aa_eqtl_alt)){
  print(i)
  aa_snp_inves[[i]] <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                                              my_motifs = TF.motifs.Expanded_pseudo_exp12_t,
                                              snp_chr = aa_eqtl_chr[i],
                                              snp_start = aa_eqtl_st[i]-1,
                                              snp_ref = aa_eqtl_ref[i],
                                              snp_alt = aa_eqtl_alt[i],
                                              min_LLR = aamin_LLR_low[aa_names])
  
}



names(aa_snp_inves) <- aa_1
aarw <- unlist(lapply(aa_snp_inves, nrow))
table(aarw)


aanam <- names(aa_snp_qtl_marker_gtex_pancan_match)[aa_snp_qtl_marker_gtex_pancan_match == 1]

#aanam <- names(aa_snp_qtl_marker_gtex_pancan_match)[aa_snp_qtl_marker_gtex_match == 1]


aa11 <- intersect(aa_1, aanam)


aa11[5]
aa_snp_inves[[aa11[1]]]
aa_snp_inves[[aa11[2]]]
aa_snp_inves[[aa11[3]]]
aa_snp_inves[[aa11[4]]]
aa_snp_inves[[aa11[5]]]


aa_snp[aa_snp$V10 == aa11[5],]
aa_snp[aa_snp$V10 == aa11[4],]
# neg_793    chr17 [82022299, 82023298] 

aa_snp[aa_snp$V10 =="rs73593858",]
ER_gtex_eqtl_sig[ER_gtex_eqtl_sig$variant_id == "chr10_19246337_G_T_b38",]

aamin_LLR
aamaxL_LLR



boxplot(t(aa_pos_neg_out_snp_all[aaa,]), las = 2)


aaa <- rownames(aa_all_dif_snp)[which(aa_all_dif_snp[,"par_5966"] > 0.9)]
aa_1 <- unlist(lapply(strsplit(aaa, split="_"), "[[", 2))
aa_2 <- unlist(lapply(strsplit(aaa, split="_"), "[[", 3))
aanam <- names(aa_snp_qtl_marker_gtex_pancan_match)[aa_snp_qtl_marker_gtex_pancan_match == 1]
aamy <- aaa[aa_1%in%aanam]

for(j in 1:length(aamy)){
  for(i in 1:length(aa_pos_out_snp)){
    if(aamy[j] %in% rownames(aa_pos_out_snp[[i]])){
      print("########")
      print(aamy[j])
      print(names(aa_pos_out_snp)[i])
      print("########")
    }
  }
  for(i in 1:length(aa_neg_out_snp)){
    if(aamy[j] %in% rownames(aa_neg_out_snp[[i]])){
      print("########")
      print(aamy[j])
      print(names(aa_neg_out_snp)[i])
      print("########")
    }
  }
}
"eqtl_rs2246069_1"
"neg_1372"

"eqtl_rs73593858_1"
"neg_195"

"eqtl_rs8004904_1"
"neg_484"

"eqtl_rs4444402_2"
"neg_793"

"eqtl_rs4444402_3"
"neg_793"

"eqtl_rs4513161_2"
"neg_793"

"eqtl_rs4513161_3"
"neg_793"
aa_neg_out_snp$neg_793
aa_neg_dif_snp$neg_793

par(mfrow = c(4,2), mar = c(2,2,2,2))
par(mfrow = c(1,1), mar = c(4,4,4,4), oma = c(0,0,0,0))
boxplot.matrix(cbind(Exp16_WT_model_prediction_Filtered["neg_1372",],
                     aa_pos_neg_out_snp_all[aamy[1],]), outline = F,
               main = paste0("neg_1372","\n",aamy[1] ))

boxplot.matrix(cbind(Exp16_WT_model_prediction_Filtered["neg_195",],
                     aa_pos_neg_out_snp_all[aamy[2],]), outline = F, 
               main = paste0("neg_195","\n",aamy[2] ))

boxplot.matrix(cbind(Exp16_WT_model_prediction_Filtered["neg_484",],
                     aa_pos_neg_out_snp_all[aamy[3],]), outline = F, 
               main = paste0("neg_484","\n",aamy[3] ))

boxplot.matrix(cbind(Exp16_WT_model_prediction_Filtered["neg_793",],
                     aa_pos_neg_out_snp_all[aamy[4],]), outline = F, 
               main = paste0("neg_793","\n",aamy[4] ))

boxplot.matrix(cbind(Exp16_WT_model_prediction_Filtered["neg_793",],
                     aa_pos_neg_out_snp_all[aamy[5],]), outline = F, 
               main = paste0("neg_793","\n",aamy[5] ))

boxplot.matrix(cbind(Exp16_WT_model_prediction_Filtered["neg_793",],
                     aa_pos_neg_out_snp_all[aamy[6],]), outline = F,
               main = paste0("neg_793","\n",aamy[6] ))

boxplot.matrix(cbind(Exp16_WT_model_prediction_Filtered["neg_793",],
                     aa_pos_neg_out_snp_all[aamy[7],]), outline = F, 
               main = paste0("neg_793","\n",aamy[7] ))



############################################################################################################
# Check enrichments for percentile change


aa_all_dif_snp_perc <- rbind(aa_pos_dif_all_snp_percentile, aa_neg_dif_all_snp_percentile)
aa_all_dif_snp_perc_names <- unlist(lapply(strsplit(rownames(aa_all_dif_snp_perc), split="_"), "[[", 2))


aa_snp_gwas_marker_match <-  aa_snp_gwas_marker[match(aa_all_dif_snp_perc_names, names(aa_snp_gwas_marker))]
aa_snp_qtl_marker_pancan_match <- aa_snp_qtl_marker_pancan[match(aa_all_dif_snp_perc_names, names(aa_snp_qtl_marker_pancan))]
aa_snp_qtl_marker_gtex_match <- aa_snp_qtl_marker_gtex[match(aa_all_dif_snp_perc_names, names(aa_snp_qtl_marker_gtex))]
aa_snp_qtl_marker_gtex_pancan_match <- aa_snp_qtl_marker_gtex_pancan[match(aa_all_dif_snp_perc_names, names(aa_snp_qtl_marker_gtex_pancan))]


aa_all_sig_nonsig_perc <- matrix(nrow = 10, ncol = ncol(aa_pos_dif_all_snp))
colnames(aa_all_sig_nonsig_perc) <- colnames(aa_neg_dif_all_snp_percentile)[1:(ncol(aa_neg_dif_all_snp_percentile))]

rownames(aa_all_sig_nonsig_perc) <- c("sig_all", "sig_qtl_gtex", "nonsig_all", "nonsig_qtl", "pval_pancan", "pval_gtex", "pval_both", "gwas_sig", "sig_qtl_pancan", "sig_qtl_either")



aa_thresh <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9)
aa_all_sig_hold_perc <- list()

for(i in 1:length(aa_thresh)){
  aa_all_sig_hold_perc[[i]] <- aa_all_sig_nonsig_perc
  for(j in 1:(ncol(aa_all_dif_snp_perc))){
    aa_all_sig_hold_perc[[i]][1,j] <- sum(abs(aa_all_dif_snp_perc[,j]) >= aa_thresh[i], na.rm = T)
    aa_all_sig_hold_perc[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_perc[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc[[i]][3,j] <- sum(abs(aa_all_dif_snp_perc[,j]) < aa_thresh[i], na.rm = T)
    aa_all_sig_hold_perc[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_perc[,j]) < aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp_perc[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp_perc[,j]) >= aa_thresh[i]], na.rm = T)
    
    if(sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp_perc[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc[[i]][5,j] <- phyper(q = sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp_perc[,j]) >= aa_thresh[i]], na.rm = T),
                                          m = sum(aa_snp_qtl_marker_pancan_match == 1),
                                          n = sum(aa_snp_qtl_marker_pancan_match == 0), 
                                          k = aa_all_sig_hold_perc[[i]][1,j],
                                          lower.tail = F)
    }
    
    if(sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_perc[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc[[i]][6,j] <- phyper(q = sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_perc[,j]) >= aa_thresh[i]], na.rm = T),
                                          m = sum(aa_snp_qtl_marker_gtex_match == 1),
                                          n = sum(aa_snp_qtl_marker_gtex_match == 0), 
                                          k = aa_all_sig_hold_perc[[i]][1,j],
                                          lower.tail = F)
    }
    
    if(sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp_perc[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc[[i]][7,j] <- phyper(q = sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp_perc[,j]) >= aa_thresh[i]], na.rm = T),
                                          m = sum(aa_snp_qtl_marker_gtex_pancan_match == 1),
                                          n = sum(aa_snp_qtl_marker_gtex_pancan_match == 0), 
                                          k = aa_all_sig_hold_perc[[i]][1,j],
                                          lower.tail = F)
    }
    
    
    aa_all_sig_hold_perc[[i]][8,j] <- sum(aa_snp_gwas_marker_match[abs(aa_all_dif_snp_perc[,j]) >= aa_thresh[i]], na.rm = T)
    
  }
}
names(aa_all_sig_hold_perc) <- aa_thresh

for(i in 1:length(aa_all_sig_hold_perc)){
  print(names(aa_all_sig_hold_perc)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_perc[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_perc[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_perc[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_perc[[i]][8,] > 0, na.rm=T))
}
View(aa_all_sig_hold_perc$`0.5`)



##################################################################################################################################

# write GEMSTAT jobs to run
# write hal job for KDs
aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas1[1:150]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas2[1:150]]
aaadd3 <- union(aaadd1, aaadd2)

aa_par_name <- paste0(aaadd3, ".txt")
aa_seq_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/Variant_seq/",
                           recursive = T)
aa_lab_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/Variant_Labels/",
                           recursive = T)
aa_nam_spl <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_seq_files, 
                                                                   split = "\\."), 
                                                          "[[", 1)), split = "_"),
                                   "[", c(3,4)),
                            paste, collapse = "_"))


setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18")
for(i in 1:length(aa_par_name)){
  for(j in 1:length(aa_seq_files)){
    cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr ",
          paste0("-s Variant_seq/", aa_seq_files[j]), 
          paste0("-e Variant_Labels/", aa_lab_files[j]),
          "-m motifs.wtmx -f TF_exp.tab", 
          paste0("-fo Varaint_out/", aa_nam_spl[j], "_", aaadd3[i], ".out"), 
          "-o DIRECT -c Coop/coop.par ", 
          paste0("-p Trained_par/", aa_par_name[i]),
          "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1\n"
    ), 
    sep = " ", append = !(i==1 & j==1), file = "variant_jobs_snp_exp18.job")
  }
}

################# ################# ################# ################# #################
# read results of the snps
# read output and find the ones with some effect
# first read WT for all models

aa_wt_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/Varaint_out/",
                              pattern = "NA_NA_trained_par*", full.names = T)
aa_wt_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/Varaint_out/",
                          pattern = "NA_NA_trained_par*", full.names = F)
aa_nam_spl <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_wt_files, split = "\\."), "[[", 1)), split = "_"), "[", c(4,5)), paste, collapse = "_"))

aa_t <- read_output_train_test_GEMSTAT_indiv(output_file = aa_wt_files_ful[1])
aa_wt_models_snp <- matrix(nrow=nrow(aa_t$GT),
                           ncol=length(aa_wt_files))

colnames(aa_wt_models_snp) <- aa_nam_spl
rownames(aa_wt_models_snp) <- aa_t$pred$Rows
for(i in 1:ncol(aa_wt_models_snp)){
  aa_wt_models_snp[, i] <- read_output_train_test_GEMSTAT_indiv(output_file = aa_wt_files_ful[i])$pred$X1
}
aa_neg_out_snp_2 <- list()
aa_pos_out_snp_2 <- list()

aa_neg_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/Varaint_out/",
                               pattern = "neg_*", full.names = T)
aa_neg_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/Varaint_out/",
                           pattern = "neg_*", full.names = F)
aa_neg_sp_1 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_neg_files, split = "\\."), "[[", 1)), split = "_"), "[", c(1,2)), paste, collapse = "_"))
aa_neg_sp_2 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_neg_files, split = "\\."), "[[", 1)), split = "_"), "[", c(4,5)), paste, collapse = "_"))


# there was a problem: in 9648, had to fix
aa_neg_sp_1[9648] <- "neg_1033"
aa_neg_sp_111 <- c(aa_neg_sp_1[1:9648], aa_neg_sp_1[9650:length(aa_neg_sp_1)])
aa_neg_sp_1 <- aa_neg_sp_111
aa_neg_uniq <- unique(aa_neg_sp_1)

aa_pos_files_ful <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/Varaint_out/",
                               pattern = "pos_*", full.names = T)
aa_pos_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/Varaint_out/",
                           pattern = "pos_*", full.names = F)
aa_pos_sp_1 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_pos_files, split = "\\."), "[[", 1)), split = "_"), "[", c(1,2)), paste, collapse = "_"))
aa_pos_sp_2 <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_pos_files, split = "\\."), "[[", 1)), split = "_"), "[", c(4,5)), paste, collapse = "_"))
aa_pos_uniq <- unique(aa_pos_sp_1)

aa_neg_dif_snp_2 <- list()
aa_pos_dif_snp_2 <- list()

aa_neg_dif_snp_percentile_2 <- list()
aa_pos_dif_snp_percentile_2 <- list()

aa_Exp18_WT_model_prediction_Filtered_GT <- rownames(Exp18_WT_model_prediction_Filtered)
aa_Exp18_WT_model_prediction_Filtered_GT <- unlist(lapply(strsplit(aa_Exp18_WT_model_prediction_Filtered_GT, split = "_"), "[[", 1))
aa_Exp18_WT_model_prediction_Filtered_GT[aa_Exp18_WT_model_prediction_Filtered_GT == "pos"] <- 1
aa_Exp18_WT_model_prediction_Filtered_GT[aa_Exp18_WT_model_prediction_Filtered_GT == "neg"] <- 0
aa_Exp18_WT_model_prediction_Filtered_GT <- as.numeric(aa_Exp18_WT_model_prediction_Filtered_GT)
names(aa_Exp16_WT_model_prediction_Filtered_GT) <- rownames(Exp16_WT_model_prediction_Filtered)





for(i in 1:length(aa_neg_uniq)){
  print(i)
  aa_cur_ind <- which(aa_neg_sp_1 %in% aa_neg_uniq[i])
  aa_t <- read.table(file = aa_neg_files_ful[aa_cur_ind[1]], header = T, stringsAsFactors = F)
  aa_ind <- seq(2, nrow(aa_t), 2)
  aa_t <- aa_t[aa_ind, ]
  aa_neg_out_snp_2[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  colnames(aa_neg_out_snp_2[[i]]) <- colnames(aa_wt_models_snp)
  rownames(aa_neg_out_snp_2[[i]]) <- aa_t$Rows
  aa_neg_dif_snp_2[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  colnames(aa_neg_dif_snp_2[[i]]) <- colnames(aa_wt_models_snp)
  rownames(aa_neg_dif_snp_2[[i]]) <- aa_t$Rows
  aa_neg_dif_snp_percentile_2[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  colnames(aa_neg_dif_snp_percentile_2[[i]]) <- colnames(aa_wt_models_snp)
  rownames(aa_neg_dif_snp_percentile_2[[i]]) <- aa_t$Rows
  for(j in 1:length(aa_cur_ind)){
    aa_table <- read.table(file =aa_neg_files_ful[aa_cur_ind[j]],
                           header = T, stringsAsFactors = F)
    aa_neg_out_snp_2[[i]][, aa_neg_sp_2[aa_cur_ind[j]]] <- aa_table[seq(2, nrow(aa_table), 2),"X1"]
  }
  for(j in 1:nrow(aa_neg_out_snp_2[[i]])){
    aa_neg_dif_snp_2[[i]][j, ] <- (aa_neg_out_snp_2[[i]][j, ] - Exp18_WT_model_prediction_Filtered[aa_neg_uniq[i], ])/ Exp18_WT_model_prediction_Filtered[aa_neg_uniq[i], ]
    for(kk in 1:ncol(aa_neg_out_snp_2[[i]])){
      if(!is.na(Exp18_WT_model_prediction_Filtered[aa_neg_uniq[i],kk])){
        aa_pos_perc_before <- sum(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 1,kk] <= Exp18_WT_model_prediction_Filtered[aa_neg_uniq[i],kk], na.rm = T)/sum(!is.na(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 1,kk]))
        aa_neg_perc_before <- sum(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 0,kk] >  Exp18_WT_model_prediction_Filtered[aa_neg_uniq[i],kk], na.rm = T)/sum(!is.na(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 0,kk]))
        aa_pos_perc_after <- sum(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 1,kk] <=  aa_neg_out_snp_2[[i]][j, kk], na.rm = T)/sum(!is.na(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 1,kk]))
        aa_neg_perc_after <- sum(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 0,kk] >   aa_neg_out_snp_2[[i]][j, kk], na.rm = T)/sum(!is.na(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 0,kk]))
        aa_neg_dif_snp_percentile_2[[i]][j,kk] <- -(aa_neg_perc_after - aa_neg_perc_before) + (aa_pos_perc_after - aa_pos_perc_before)
        # print("aa_pos_perc_before")
        # print(aa_pos_perc_before)
        # print("aa_neg_perc_before")
        # print(aa_neg_perc_before)
        # print("aa_pos_perc_after")
        # print(aa_pos_perc_after)
        # print("aa_neg_perc_after")
        # print(aa_neg_perc_after)
        
        
        
        
      }
    }
  }
}
names(aa_neg_out_snp_2) <- aa_neg_uniq
names(aa_neg_dif_snp_2) <- aa_neg_uniq
names(aa_neg_dif_snp_percentile_2) <- aa_neg_uniq
aa_neg_out_snp_all_2 <- do.call(what = rbind, aa_neg_out_snp_2)
aa_neg_dif_snp_percentile_2 <- do.call(what = rbind, aa_neg_dif_snp_percentile_2)

aa_neg_dif_snp_2$neg_990[,10]
aa_neg_dif_snp_percentile_2[rownames(aa_neg_dif_snp_2$neg_990), 10]


for(i in 1:length(aa_pos_uniq)){
  print(i)
  aa_cur_ind <- which(aa_pos_sp_1 %in% aa_pos_uniq[i])
  aa_t <- read.table(file = aa_pos_files_ful[aa_cur_ind[1]], header = T, stringsAsFactors = F)
  aa_ind <- seq(2, nrow(aa_t), 2)
  aa_t <- aa_t[aa_ind, ]
 # aa_pos_out_snp_2[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
#  colnames(aa_pos_out_snp_2[[i]]) <- colnames(aa_wt_models_snp)
#  rownames(aa_pos_out_snp_2[[i]]) <- aa_t$Rows
  aa_pos_dif_snp_2[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  colnames(aa_pos_dif_snp_2[[i]]) <- colnames(aa_wt_models_snp)
  rownames(aa_pos_dif_snp_2[[i]]) <- aa_t$Rows
  aa_pos_dif_snp_percentile_2[[i]] <- matrix(nrow = nrow(aa_t), ncol = ncol(aa_wt_models_snp))
  colnames(aa_pos_dif_snp_percentile_2[[i]]) <- colnames(aa_wt_models_snp)
  rownames(aa_pos_dif_snp_percentile_2[[i]]) <- aa_t$Rows
  
  # for(j in 1:length(aa_cur_ind)){
  #   aa_table <- read.table(file =aa_pos_files_ful[aa_cur_ind[j]],
  #                          header = T, stringsAsFactors = F)
  #   aa_pos_out_snp_2[[i]][, aa_pos_sp_2[aa_cur_ind[j]]] <- aa_table[seq(2, nrow(aa_table), 2),"X1"]
  # }
  for(j in 1:nrow(aa_pos_out_snp_2[[i]])){
    aa_pos_dif_snp_2[[i]][j, ] <- (aa_pos_out_snp_2[[i]][j, ] - Exp18_WT_model_prediction_Filtered[aa_pos_uniq[i], ])/ Exp18_WT_model_prediction_Filtered[aa_pos_uniq[i], ]
    for(kk in 1:ncol(aa_pos_out_snp_2[[i]])){
      if(!is.na(Exp18_WT_model_prediction_Filtered[aa_pos_uniq[i],kk])){
        aa_pos_perc_before <- sum(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 1,kk] <= Exp18_WT_model_prediction_Filtered[aa_pos_uniq[i],kk], na.rm = T)/sum(!is.na(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 1,kk]))
        aa_neg_perc_before <- sum(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 0,kk] >  Exp18_WT_model_prediction_Filtered[aa_pos_uniq[i],kk], na.rm = T)/sum(!is.na(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 0,kk]))
        aa_pos_perc_after <-  sum(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 1,kk] <=  aa_pos_out_snp_2[[i]][j, kk], na.rm = T)/sum(!is.na(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 1,kk]))
        aa_neg_perc_after <-  sum(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 0,kk] >   aa_pos_out_snp_2[[i]][j, kk], na.rm = T)/sum(!is.na(Exp18_WT_model_prediction_Filtered[aa_Exp18_WT_model_prediction_Filtered_GT == 0,kk]))
        aa_pos_dif_snp_percentile_2[[i]][j,kk] <- -(aa_neg_perc_after - aa_neg_perc_before) + (aa_pos_perc_after - aa_pos_perc_before)
      }
    }
  }
  
}
names(aa_pos_out_snp_2) <- aa_pos_uniq
names(aa_pos_dif_snp_2) <- aa_pos_uniq
names(aa_pos_dif_snp_percentile_2) <- aa_pos_uniq

aa_pos_dif_snp_2$pos_124[,1]
aa_pos_dif_snp_percentile_2$pos_124[,1]


png(filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Plots/exp18_plot_pos_neg_filtered.png",   
    width = 8*300,        # 5 x 300 pixels
    height = 30*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
par(mfrow = c(31,8), mar = c(1,1,1,1))
for(i in 1:ncol(Exp18_WT_model_prediction_Filtered)){
  boxplot(Exp18_WT_model_prediction_Filtered[, i]~aa_Exp18_WT_model_prediction_Filtered_GT)
}
dev.off()

png(filename = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Plots/exp18_plot_pos_neg.png",   
    width = 8*300,        # 5 x 300 pixels
    height = 30*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
par(mfrow = c(31,8), mar = c(1,1,1,1))
for(i in 1:ncol(Exp18_WT_model_prediction_Filtered)){
  boxplot(Exp18_WT_model_prediction[, i]~aa_Exp18_WT_model_prediction_Filtered_GT)
}
dev.off()


aa_pos_out_snp_all_2 <- do.call(what = rbind, aa_pos_out_snp_2)

# aa_pos_dif_snp_old_noFilt <- aa_pos_dif_snp # before filtering the outputs of WT models
# aa_neg_dif_snp_old_noFilt <- aa_neg_dif_snp
aa_pos_dif_all_snp_2 <- do.call(what = rbind, aa_pos_dif_snp_2)
aa_pos_dif_all_snp_percentile_2 <- do.call(what = rbind, aa_pos_dif_snp_percentile_2)

par(mfrow = c(1,1), mar = c(4,4,4,4))
boxplot.matrix(t(aa_pos_dif_all_snp_2[1:5,])*100)

aa_neg_dif_all_snp_2 <- do.call(what = rbind, aa_neg_dif_snp_2)
#boxplot.matrix(t(aa_neg_dif_all_snp))
sum(colSums(is.na(aa_pos_dif_all_snp_2)) == nrow(aa_pos_dif_all_snp_2))

aa_pos_neg_out_snp_all_2 <- rbind(aa_pos_out_snp_all_2, aa_neg_out_snp_all_2)




aa_all_dif_snp_2 <- rbind(aa_pos_dif_all_snp_2, aa_neg_dif_all_snp_2)
aa_all_dif_snp_names <- unlist(lapply(strsplit(rownames(aa_all_dif_snp), split="_"), "[[", 2))


aa_snp_gwas_marker_match <-  aa_snp_gwas_marker[match(aa_all_dif_snp_names, names(aa_snp_gwas_marker))]
aa_snp_qtl_marker_pancan_match <- aa_snp_qtl_marker_pancan[match(aa_all_dif_snp_names, names(aa_snp_qtl_marker_pancan))]
aa_snp_qtl_marker_gtex_match <- aa_snp_qtl_marker_gtex[match(aa_all_dif_snp_names, names(aa_snp_qtl_marker_gtex))]
aa_snp_qtl_marker_gtex_pancan_match <- aa_snp_qtl_marker_gtex_pancan[match(aa_all_dif_snp_names, names(aa_snp_qtl_marker_gtex_pancan))]


aa_all_sig_nonsig_2 <- matrix(nrow = 10, ncol = ncol(aa_pos_dif_all_snp_2))
colnames(aa_all_sig_nonsig_2) <- colnames(aa_neg_dif_all_snp_2)[1:(ncol(aa_neg_dif_all_snp_2))]

rownames(aa_all_sig_nonsig_2) <- c("sig_all", "sig_qtl_gtex", "nonsig_all", "nonsig_qtl", "pval_pancan", "pval_gtex", "pval_both", "gwas_sig", "sig_qtl_pancan", "sig_qtl_either")



aa_thresh <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
aa_all_sig_hold_2 <- list()
aa_all_sig_hold_2_uniq <- list()

for(i in 1:length(aa_thresh)){
  print(i)
  aa_all_sig_hold_2[[i]] <- aa_all_sig_nonsig_2
  aa_all_sig_hold_2_uniq[[i]] <- aa_all_sig_nonsig_2
  for(j in 1:(ncol(aa_all_dif_snp_2))){
    aa_all_sig_hold_2[[i]][1,j] <- sum(abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i], na.rm = T)
    aa_m_uniq_up <- unique(aa_all_dif_snp_names[which(abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i])])
    aa_all_sig_hold_2_uniq[[i]][1,j] <- length(aa_m_uniq_up)
    aa_all_sig_hold_2[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_2_uniq[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex[match(aa_m_uniq_up, names(aa_snp_qtl_marker_gtex))])
    aa_all_sig_hold_2[[i]][3,j] <- sum(abs(aa_all_dif_snp_2[,j]) < aa_thresh[i], na.rm = T)
    aa_m_uniq_down <- unique(aa_all_dif_snp_names[which(abs(aa_all_dif_snp_2[,j]) < aa_thresh[i])])
    aa_all_sig_hold_2_uniq[[i]][3,j] <- length(aa_m_uniq_down)
    aa_all_sig_hold_2[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_2[,j]) < aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_2_uniq[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex[match(aa_m_uniq_down, names(aa_snp_qtl_marker_gtex))])
    aa_all_sig_hold_2[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_2_uniq[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan[match(aa_m_uniq_up, names(aa_snp_qtl_marker_pancan))])
    aa_all_sig_hold_2[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_2_uniq[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan[match(aa_m_uniq_up, names(aa_snp_qtl_marker_gtex_pancan))])
    if(sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_2[[i]][5,j] <- phyper(q = sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i]], na.rm = T),
                                          m = sum(aa_snp_qtl_marker_pancan_match == 1),
                                          n = sum(aa_snp_qtl_marker_pancan_match == 0), 
                                          k = aa_all_sig_hold_2[[i]][1,j],
                                          lower.tail = F)
      aa_all_sig_hold_2_uniq[[i]][5,j] <- phyper(q = aa_all_sig_hold_2_uniq[[i]][9,j],
                                                 m = sum(aa_snp_qtl_marker_pancan == 1),
                                                 n = sum(aa_snp_qtl_marker_pancan == 0), 
                                                 k = aa_all_sig_hold_2_uniq[[i]][1,j],
                                                 lower.tail = F)
    }
    
    if(sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_2[[i]][6,j] <- phyper(q = sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i]], na.rm = T),
                                          m = sum(aa_snp_qtl_marker_gtex_match == 1),
                                          n = sum(aa_snp_qtl_marker_gtex_match == 0), 
                                          k = aa_all_sig_hold_2[[i]][1,j],
                                          lower.tail = F)
      aa_all_sig_hold_2_uniq[[i]][6,j] <- phyper(q = aa_all_sig_hold_2_uniq[[i]][2,j],
                                            m = sum(aa_snp_qtl_marker_gtex == 1),
                                            n = sum(aa_snp_qtl_marker_gtex == 0), 
                                            k = aa_all_sig_hold_2_uniq[[i]][1,j],
                                            lower.tail = F)
    }
    
    if(sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_2[[i]][7,j] <- phyper(q = sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i]], na.rm = T),
                                          m = sum(aa_snp_qtl_marker_gtex_pancan_match == 1),
                                          n = sum(aa_snp_qtl_marker_gtex_pancan_match == 0), 
                                          k = aa_all_sig_hold_2[[i]][1,j],
                                          lower.tail = F)
      aa_all_sig_hold_2_uniq[[i]][7,j] <- phyper(q = aa_all_sig_hold_2_uniq[[i]][10,j],
                                            m = sum(aa_snp_qtl_marker_gtex_pancan == 1),
                                            n = sum(aa_snp_qtl_marker_gtex_pancan == 0), 
                                            k = aa_all_sig_hold_2_uniq[[i]][1,j],
                                            lower.tail = F)
    }
    
    
    aa_all_sig_hold_2[[i]][8,j] <- sum(aa_snp_gwas_marker_match[abs(aa_all_dif_snp_2[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_2_uniq[[i]][8,j] <- sum(aa_snp_gwas_marker[match(aa_m_uniq_up, names(aa_snp_gwas_marker))], na.rm = T)
    
  }
}
names(aa_all_sig_hold_2) <- aa_thresh
names(aa_all_sig_hold_2_uniq) <- aa_thresh

for(i in 1:length(aa_all_sig_hold_2)){
  print(names(aa_all_sig_hold_2)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_2[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_2[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_2[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_2[[i]][8,] > 0, na.rm=T))
}

for(i in 1:length(aa_all_sig_hold_2_uniq)){
  print(names(aa_all_sig_hold_2_uniq)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_2_uniq[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_2_uniq[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_2_uniq[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_2_uniq[[i]][8,] > 0, na.rm=T))
  print("################")
}
aaxxx <- sort(aa_all_sig_hold_2_uniq$`0.25`[5, ], decreasing = F, index.return=T)$ix
View(aa_all_sig_hold_2_uniq$`0.25`[, aaxxx[1:10]])
##########

aa_all_dif_snp_perc_2 <- rbind(aa_pos_dif_all_snp_percentile_2, aa_neg_dif_snp_percentile_2)
aa_all_dif_snp_perc_names <- unlist(lapply(strsplit(rownames(aa_all_dif_snp_perc_2), split="_"), "[[", 2))


aa_snp_gwas_marker_match <-  aa_snp_gwas_marker[match(aa_all_dif_snp_perc_names, names(aa_snp_gwas_marker))]
aa_snp_qtl_marker_pancan_match <- aa_snp_qtl_marker_pancan[match(aa_all_dif_snp_perc_names, names(aa_snp_qtl_marker_pancan))]
aa_snp_qtl_marker_gtex_match <- aa_snp_qtl_marker_gtex[match(aa_all_dif_snp_perc_names, names(aa_snp_qtl_marker_gtex))]
aa_snp_qtl_marker_gtex_pancan_match <- aa_snp_qtl_marker_gtex_pancan[match(aa_all_dif_snp_perc_names, names(aa_snp_qtl_marker_gtex_pancan))]


aa_all_sig_nonsig_perc_2 <- matrix(nrow = 10, ncol = ncol(aa_pos_dif_all_snp_2))
colnames(aa_all_sig_nonsig_perc_2) <- colnames(aa_neg_dif_snp_percentile_2)[1:(ncol(aa_neg_dif_snp_percentile_2))]

rownames(aa_all_sig_nonsig_perc_2) <- c("sig_all", "sig_qtl_gtex", "nonsig_all", "nonsig_qtl", "pval_pancan", "pval_gtex", "pval_both", "gwas_sig", "sig_qtl_pancan", "sig_qtl_either")



aa_thresh <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
               1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9)
aa_all_sig_hold_perc_2 <- list()
aa_all_sig_hold_perc_2_uniq <- list()
for(i in 1:length(aa_thresh)){
  print(aa_thresh[i])
  aa_all_sig_hold_perc_2[[i]] <- aa_all_sig_nonsig_perc_2
  aa_all_sig_hold_perc_2_uniq[[i]] <- aa_all_sig_nonsig_perc_2
  for(j in 1:(ncol(aa_all_dif_snp_perc_2))){
    aa_all_sig_hold_perc_2[[i]][1,j] <- sum(abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i], na.rm = T)
    aa_m_uniq_up <- unique(aa_all_dif_snp_perc_names[which(abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i])])
    aa_all_sig_hold_perc_2_uniq[[i]][1,j] <- length(aa_m_uniq_up)
    aa_all_sig_hold_perc_2[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_2_uniq[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex[match(aa_m_uniq_up, names(aa_snp_qtl_marker_gtex))])
    aa_all_sig_hold_perc_2[[i]][3,j] <- sum(abs(aa_all_dif_snp_perc_2[,j]) < aa_thresh[i], na.rm = T)
    aa_m_uniq_down <- unique(aa_all_dif_snp_perc_names[which(abs(aa_all_dif_snp_perc_2[,j]) < aa_thresh[i])])
    aa_all_sig_hold_perc_2_uniq[[i]][3,j] <- length(aa_m_uniq_down)
    aa_all_sig_hold_perc_2[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_perc_2[,j]) < aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_2_uniq[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex[match(aa_m_uniq_down, names(aa_snp_qtl_marker_gtex))])
    aa_all_sig_hold_perc_2[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_2_uniq[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan[match(aa_m_uniq_up, names(aa_snp_qtl_marker_pancan))])
    aa_all_sig_hold_perc_2[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_2_uniq[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan[match(aa_m_uniq_up, names(aa_snp_qtl_marker_gtex_pancan))])
    
    if(sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc_2[[i]][5,j] <- phyper(q = sum(aa_snp_qtl_marker_pancan_match[abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i]], na.rm = T),
                                               m = sum(aa_snp_qtl_marker_pancan_match == 1),
                                               n = sum(aa_snp_qtl_marker_pancan_match == 0), 
                                               k = aa_all_sig_hold_perc_2[[i]][1,j],
                                               lower.tail = F)
      aa_all_sig_hold_perc_2_uniq[[i]][5,j] <- phyper(q = aa_all_sig_hold_perc_2_uniq[[i]][9,j],
                                                      m = sum(aa_snp_qtl_marker_pancan == 1),
                                                      n = sum(aa_snp_qtl_marker_pancan == 0), 
                                                      k = aa_all_sig_hold_perc_2_uniq[[i]][1,j],
                                                      lower.tail = F)
    }
    
    if(sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc_2[[i]][6,j] <- phyper(q = sum(aa_snp_qtl_marker_gtex_match[abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i]], na.rm = T),
                                               m = sum(aa_snp_qtl_marker_gtex_match == 1),
                                               n = sum(aa_snp_qtl_marker_gtex_match == 0), 
                                               k = aa_all_sig_hold_perc_2[[i]][1,j],
                                               lower.tail = F)
      aa_all_sig_hold_perc_2_uniq[[i]][6,j] <- phyper(q = aa_all_sig_hold_perc_2_uniq[[i]][2,j],
                                                 m = sum(aa_snp_qtl_marker_gtex == 1),
                                                 n = sum(aa_snp_qtl_marker_gtex == 0), 
                                                 k = aa_all_sig_hold_perc_2_uniq[[i]][1,j],
                                                 lower.tail = F)
      
    }
    
    if(sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc_2[[i]][7,j] <- phyper(q = sum(aa_snp_qtl_marker_gtex_pancan_match[abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i]], na.rm = T),
                                               m = sum(aa_snp_qtl_marker_gtex_pancan_match == 1),
                                               n = sum(aa_snp_qtl_marker_gtex_pancan_match == 0), 
                                               k = aa_all_sig_hold_perc_2[[i]][1,j],
                                               lower.tail = F)
      aa_all_sig_hold_perc_2_uniq[[i]][7,j] <- phyper(q = aa_all_sig_hold_perc_2_uniq[[i]][10,j],
                                                 m = sum(aa_snp_qtl_marker_gtex_pancan == 1),
                                                 n = sum(aa_snp_qtl_marker_gtex_pancan == 0), 
                                                 k = aa_all_sig_hold_perc_2_uniq[[i]][1,j],
                                                 lower.tail = F)
    }
    
    
    aa_all_sig_hold_perc_2[[i]][8,j] <- sum(aa_snp_gwas_marker_match[abs(aa_all_dif_snp_perc_2[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_2_uniq[[i]][8,j] <- sum(aa_snp_gwas_marker[match(aa_m_uniq_up, names(aa_snp_gwas_marker))], na.rm = T)
    
  }
}
names(aa_all_sig_hold_perc_2) <- aa_thresh
names(aa_all_sig_hold_perc_2_uniq) <- aa_thresh

for(i in 1:length(aa_all_sig_hold_perc_2)){
  print(names(aa_all_sig_hold_perc_2)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_perc_2[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_perc_2[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_perc_2[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_perc_2[[i]][8,] > 0, na.rm=T))
  print("##############")
}

for(i in 1:length(aa_all_sig_hold_perc_2_uniq)){
  print(names(aa_all_sig_hold_perc_2_uniq)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_perc_2_uniq[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_perc_2_uniq[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_perc_2_uniq[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_perc_2_uniq[[i]][8,] > 0, na.rm=T))
  print("##############")
}
View(aa_all_sig_hold_perc_2_uniq$`0.05`[])

aax <- sort(aa_all_sig_hold_perc_2_uniq$`0.05`[6,], decreasing = F, index.return =T)$ix

View(aa_all_sig_hold_perc_2_uniq$`0.05`[, aax[1:16]])

which(aa_all_sig_hold_perc_2_uniq$`0.4`[5,] < 0.05)
summary(aa_all_sig_hold_perc_2_uniq$`0.05`[8,])
summary(aa_all_sig_hold_perc_2_uniq[[3]][8,])
which(aa_all_sig_hold_perc_2_uniq[[1]][8,] == 2)
aa_all_sig_hold_perc_2_uniq[[4]][,136]
View(aa_all_sig_hold_perc_2_uniq[[1]][, which(aa_all_sig_hold_perc_2_uniq[[1]][8,] == 2)])

####################################
# identify the TFs responsible for the eqtls found to be impactful bt models

aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_snp <- readPeakFile(aa_snp_Add, as = "data.frame")
aa_snp$V10 <- as.character(levels(aa_snp$V10)[as.numeric(aa_snp$V10)])
aa_snp$V9 <- as.character(levels(aa_snp$V9)[as.numeric(aa_snp$V9)])
aa_snp$V8 <- as.character(levels(aa_snp$V8)[as.numeric(aa_snp$V8)])
aa_snp$V1 <- as.character(levels(aa_snp$V1)[as.numeric(aa_snp$V1)])
aa_snp$V2 <- aa_snp$V2 - 1

aasrt <- sort(aa_all_sig_hold_perc_2_uniq$`0.05`[6,], index.return=T)$ix
aaa <- rownames(aa_all_dif_snp_perc_2)[which(abs(aa_all_dif_snp_perc_2[,aasrt[1]]) > 0.05)]

aa_1 <- unlist(lapply(strsplit(aaa, split="_"), "[[", 2))
aa_2 <- unlist(lapply(strsplit(aaa, split="_"), "[[", 3))


aa_q <- which(aa_1 %in% names(aa_snp_qtl_marker_gtex)[aa_snp_qtl_marker_gtex == 1])

aaa <- aaa[aa_q]
aa_1 <- aa_1[aa_q]
aa_2 <- aa_2[aa_q]

aa_eqtl_chr <- character(length(aaa))
aa_eqtl_st <- numeric(length(aaa))
aa_eqtl_ref <- numeric(length(aaa))
aa_eqtl_alt <- numeric(length(aaa))

for(i in 1:length(aaa)){
  aa_eqtl_chr[i] <- aa_snp$V1[aa_snp$V10 == aa_1[i]]
  aa_eqtl_st[i] <- aa_snp$V2[aa_snp$V10 == aa_1[i]]
  aa_eqtl_ref[i] <- aa_snp$V8[aa_snp$V10 == aa_1[i]]
  aaaaalt <- unlist(strsplit(aa_snp$V9[aa_snp$V10 == aa_1[i]], split = ","))
  aa_eqtl_alt[i] <- aaaaalt[as.numeric(aa_2[i])]
}



aamin_LLR_low <- aa_exp18_minLLR - 0.1*aa_exp18_minLLR
aamin_LLR_lowest <- numeric(length(aa_exp18_minLLR))
names(aamin_LLR_lowest) <- names(aamin_LLR_low)

aa_snp_inves_2 <- list()
for(i in 1:length(aa_eqtl_alt)){
  print(i)
  aa_snp_inves_2[[i]] <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                                        my_motifs = TF.motifs.Expanded_pseudo_exp12_t,
                                        snp_chr = aa_eqtl_chr[i],
                                        snp_start = aa_eqtl_st[i]-1,
                                        snp_ref = aa_eqtl_ref[i],
                                        snp_alt = aa_eqtl_alt[i],
                                        min_LLR = aa_exp18_minLLR[aa_names])
  
}

aaaaaa <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                                        my_motifs = TF.motifs.Expanded_pseudo_exp12_t,
                                        snp_chr = aa_eqtl_chr[i],
                                        snp_start = aa_eqtl_st[i]-1,
                                        snp_ref = aa_eqtl_ref[i],
                                        snp_alt = aa_eqtl_alt[i],
                                        min_LLR = aa_exp18_minLLR[aa_names])

names(aa_snp_inves_2) <- aaa
aarw <- unlist(lapply(aa_snp_inves_2, nrow))
table(aarw)

i<-9
aa_snp_inves_2[[i]]
names(aa_snp_inves_2)[i]
# find out how many are on positive seq and how many on negative


aa_dpn <- character(length(aaa))
names(aa_dpn) <- aaa
for(j in 1:length(aaa)){
  for(i in 1:length(aa_pos_out_snp_2)){
    if(aaa[j] %in% rownames(aa_pos_out_snp_2[[i]])){
      aa_dpn[j] <- names(aa_pos_out_snp_2)[i]
      print("########")
      print(aaa[j])
      print(names(aa_pos_out_snp_2)[i])
      print("########")
    }
  }
  for(i in 1:length(aa_neg_out_snp_2)){
    if(aaa[j] %in% rownames(aa_neg_out_snp_2[[i]])){
      aa_dpn[j] <- names(aa_neg_out_snp_2)[i]
      print("########")
      print(aaa[j])
      print(names(aa_neg_out_snp_2)[i])
      print("########")
    }
  }
}

####################






aa_pos_neg
Enhancer.ReMapchip.Overlap.byEnhancer
enhancer_gene_ovlap_interact
enhancer_enhancer_interact

rownames(Enhancer.ReMapchip.Overlap.byEnhancer$OverlapMat) <- names(aa_pos_neg)
rownames(Enhancer.ReMapchip.Overlap.byEnhancer$IntMat ) <- names(aa_pos_neg)


aa_qtl_chip_ov <- Enhancer.ReMapchip.Overlap.byEnhancer$OverlapMat[aa_dpn,]
aa_qtl_chip_in <- Enhancer.ReMapchip.Overlap.byEnhancer$IntMat[aa_dpn,]

aa_snp_inves_2[[1]]


aa_resp_TF <- matrix(nrow = length(aa_snp_inves_2), ncol = 4)
for(i in 1:length(aa_snp_inves_2)){
  aacnt <- 1
  for(j in 1:nrow(aa_snp_inves_2[[i]])){
    if((! aa_snp_inves_2[[i]][j,4] == 0) & (! duplicated(aa_snp_inves_2[[i]][,1])[j]) ){
      
      aa_resp_TF[i, aacnt] <- aa_snp_inves_2[[i]][j, 1]
      aacnt <- aacnt + 1
    }
  }
}
rownames(aa_resp_TF) <- names(aa_snp_inves_2)




aa_allTF <- unique(as.character(aa_resp_TF))
aa_allTF <- aa_allTF[! is.na(aa_allTF)]
aa_allTF[aa_allTF == "ESR1_2"] <- "ESR1"

aaint <- intersect(aa_allTF, colnames(aa_qtl_chip_ov))
intersect(colnames(aa_qtl_chip_ov), aa_names)
aa_resp_TF[aa_resp_TF == "ESR1_2"] <- "ESR1"
write.csv(x = aa_resp_TF, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/eqtl_TF.csv")

aa_Target_list <- list()
for(i in 1:length(aaint)){
  aa_Target_list[[i]] <- character(0)
  for(j in 1:nrow(aa_resp_TF)){
    if(aaint[i] %in% aa_resp_TF[j, ]){
      aa_Target_list[[i]] <- c(aa_Target_list[[i]], aa_dpn[rownames(aa_resp_TF)[j]]) 
    }
  }
}

names(aa_Target_list) <- aaint
aa_Target_list

for(i in 1:length(aa_Target_list)){
  print(names(aa_Target_list)[i])
  print("overlap")
  print(aa_qtl_chip_ov[aa_Target_list[[i]], names(aa_Target_list)[i]])
  print("interact")
  print(aa_qtl_chip_in[aa_Target_list[[i]], names(aa_Target_list)[i]])
  print("################")
}
i <- 2
aa_qtl_chip_ov[aa_Target_list[[i]], names(aa_Target_list)[i]]



names(aa_dpn)[aa_dpn == "pos_185"]
aa_resp_TF["eqtl_rs8076354_1",]


aa_snp_inves_2[["eqtl_rs7674850_1"]]
aa_dpn["eqtl_rs7674850_1"]

aa_snp_inves_2[["eqtl_rs4818018_2"]]
aa_dpn["eqtl_rs4818018_2"]

aa_snp_inves_2[["eqtl_rs10917374_1"]] # check this one with ER RAR binding chip
aa_dpn["eqtl_rs10917374_1"]

aa_snp_inves_2[["eqtl_rs574134_1"]] 
aa_dpn["eqtl_rs574134_1"]

aa_snp_inves_2[["eqtl_rs58757492_2"]] 
aa_dpn["eqtl_rs58757492_2"]


aasrt <- sort(aa_all_sig_hold_perc_2_uniq$`0.05`[6,], index.return=T)$ix
aaa <- rownames(aa_all_dif_snp_perc_2)[which(abs(aa_all_dif_snp_perc_2[,aasrt[1]]) > 0.05)]

i <- 54
aa_snp_inves_2[[i]]
names(aa_snp_inves_2)[i]
anam <- names(aa_snp_inves_2)[i]
aa_dpn[anam]

boxplot(cbind(Exp18_WT_model_prediction_Filtered[aa_dpn[anam],], 
              aa_pos_neg_out_snp_all_2[anam,]), outline=F)
Exp18_WT_model_prediction_Filtered[aa_dpn[anam],aasrt[1]]
aa_pos_neg_out_snp_all_2[anam,aasrt[1]]
aa_all_dif_snp_perc_2[anam,aasrt[1]]
aa_all_dif_snp_2[anam,aasrt[1]]

aa_qtl_chip_ov[aa_dpn[anam], aa_resp_TF[anam, 1]]
aa_qtl_chip_ov[aa_dpn[anam], aa_resp_TF[anam, 2]]
enhancer_gene_ovlap_interact$IntList[[1]][[aa_dpn[anam]]]
enhancer_gene_ovlap_interact$OverlapList[[1]][[aa_dpn[anam]]]
aa_dpn[anam]
# colnames(enhancer_gene_ovlap_interact$IntMat)[1]
# rownames(enhancer_gene_ovlap_interact$IntMat) <- names(aa_pos_neg)
# rownames(enhancer_gene_ovlap_interact$OverlapMat) <-  names(aa_pos_neg)
# names(enhancer_gene_ovlap_interact$IntList[[1]]) <- names(aa_pos_neg)
# names(enhancer_gene_ovlap_interact$OverlapList[[1]]) <-  names(aa_pos_neg)


aa_snp_inves_2[[5]]


################################################
# exploring what is wrong with aa_wt_models_snp not being the same as Exp18_WT_model_prediction_Filtered

aamch <- match(rownames(aa_wt_models_snp),rownames(Exp18_WT_model_prediction_Filtered))
aaww <- which(rowSums(abs(aa_wt_models_snp - Exp18_WT_model_prediction_Filtered[aamch, ]) > 1e-2, na.rm = T) > 0)

aanewnames <- unlist(lapply(lapply(strsplit(colnames(Exp18_WT_model_prediction_Filtered), "_"), "[", c(1,2)), paste, collapse = "_"))
sum(aanewnames == colnames(aa_wt_models_snp))

aww <- which(abs(Exp18_WT_model_prediction_Filtered["neg_1566",] - aa_wt_models_snp["neg_1566",]) > 0.01)

View(rbind(Exp18_WT_model_prediction_Filtered["neg_1566",aww],aa_wt_models_snp["neg_1566",aww] ))

# run GEMSTAT exp18 models on all common snps

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_snp <- readPeakFile(aa_snp_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
start(aa_snp) <- start(aa_snp) - 2
end(aa_snp) <- end(aa_snp) - 1

library(BSgenome.Hsapiens.UCSC.hg38)
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18")
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

# removing duplicates
aa_snpdf <- as.data.frame(aa_snp)
aa_is_dup1 <- duplicated(aa_snpdf[,1:(ncol(aa_snpdf) - 1)])
aaremove <- which(aa_is_dup1 > 0)
aa_snp <- aa_snp[-c(aaremove)]
aa_snpdf <- as.data.frame(aa_snp)

aa_ref <- as.character(levels(aa_snpdf$V8)[as.numeric(aa_snpdf$V8)])
aa_alt <- as.character(levels(aa_snpdf$V9)[as.numeric(aa_snpdf$V9)])

aa_ref_Alt <- cbind(aa_ref, aa_alt)
aa_names <-  as.character(levels(aa_snpdf$V10)[as.numeric(aa_snpdf$V10)])
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18")
aa_pos_neg_vars_snp <- write_variant_seq(enhancer_GR = aa_pos_neg, 
                                         eqtl_GR = aa_snp,
                                         eqtl_names = aa_names,
                                         eqtl_ref_alt = aa_ref_Alt,
                                         all_combs = F,
                                         my_genome = BSgenome.Hsapiens.UCSC.hg38, 
                                         label = c(rep(1, length(aa_pos)), rep(0, length(aa_neg))))

#setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
aa <- as.character(getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = aa_snp[2]))

load("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/snp_output_exp18_hal_WT.RData")
load("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/snp_output_exp18_hal.RData")
pos_snp_output_hal
neg_snp_output_hal
WT_models_snp_hal

aa_neg_dif_snp_3 <- list()
aa_pos_dif_snp_3 <- list()

aa_neg_dif_snp_percentile_3 <- list()
aa_pos_dif_snp_percentile_3 <- list()

aa_neg_dif_snp_percentile_4 <- list()
aa_pos_dif_snp_percentile_4 <- list()

# aa_Exp18_WT_model_prediction_Filtered_GT <- rownames(Exp18_WT_model_prediction_Filtered)
# aa_Exp18_WT_model_prediction_Filtered_GT <- unlist(lapply(strsplit(aa_Exp18_WT_model_prediction_Filtered_GT, split = "_"), "[[", 1))
# aa_Exp18_WT_model_prediction_Filtered_GT[aa_Exp18_WT_model_prediction_Filtered_GT == "pos"] <- 1
# aa_Exp18_WT_model_prediction_Filtered_GT[aa_Exp18_WT_model_prediction_Filtered_GT == "neg"] <- 0
# aa_Exp18_WT_model_prediction_Filtered_GT <- as.numeric(aa_Exp18_WT_model_prediction_Filtered_GT)
# names(aa_Exp16_WT_model_prediction_Filtered_GT) <- rownames(Exp16_WT_model_prediction_Filtered)

aa_WT_GT <- rownames(WT_models_snp_hal)
aa_WT_GT <- unlist(lapply(strsplit(aa_WT_GT, "_"), "[[", 1))
aa_WT_GT[aa_WT_GT  == "pos"] <- 1
aa_WT_GT[aa_WT_GT  == "neg"] <- 0
aa_WT_GT <- as.numeric(aa_WT_GT)
names(aa_WT_GT) <- rownames(WT_models_snp_hal)

for(i in 1:length(neg_snp_output_hal)){
  print(i)
  aa_neg_dif_snp_3[[i]] <- matrix(nrow = nrow(neg_snp_output_hal[[i]]), ncol = ncol(neg_snp_output_hal[[i]]))
  colnames(aa_neg_dif_snp_3[[i]]) <- colnames(neg_snp_output_hal[[i]])
  rownames(aa_neg_dif_snp_3[[i]]) <- rownames(neg_snp_output_hal[[i]])
  
  aa_neg_dif_snp_percentile_3[[i]] <- matrix(nrow = nrow(neg_snp_output_hal[[i]]), ncol = ncol(neg_snp_output_hal[[i]]))
  colnames(aa_neg_dif_snp_percentile_3[[i]]) <- colnames(neg_snp_output_hal[[i]])
  rownames(aa_neg_dif_snp_percentile_3[[i]]) <- rownames(neg_snp_output_hal[[i]])
  
  aa_neg_dif_snp_percentile_4[[i]] <- matrix(nrow = nrow(neg_snp_output_hal[[i]]), ncol = ncol(neg_snp_output_hal[[i]]))
  colnames(aa_neg_dif_snp_percentile_4[[i]]) <- colnames(neg_snp_output_hal[[i]])
  rownames(aa_neg_dif_snp_percentile_4[[i]]) <- rownames(neg_snp_output_hal[[i]])
  
  for(j in 1:nrow(neg_snp_output_hal[[i]])){
    aa_neg_dif_snp_3[[i]][j, ] <- (neg_snp_output_hal[[i]][j, ] - WT_models_snp_hal[names(neg_snp_output_hal)[i], ])/ WT_models_snp_hal[names(neg_snp_output_hal)[i], ]
    for(kk in 1:ncol(neg_snp_output_hal[[i]])){
      #if(!is.na(Exp18_WT_model_prediction_Filtered[names(neg_snp_output_hal)[i],kk])){
        aa_pos_perc_before <- sum(WT_models_snp_hal[aa_WT_GT == 1,kk] <= WT_models_snp_hal[names(neg_snp_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 1,kk]))
        aa_neg_perc_before <- sum(WT_models_snp_hal[aa_WT_GT == 0,kk] >  WT_models_snp_hal[names(neg_snp_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 0,kk]))
        aa_pos_perc_after <-  sum(WT_models_snp_hal[aa_WT_GT == 1,kk] <=  neg_snp_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 1,kk]))
        aa_neg_perc_after <-  sum(WT_models_snp_hal[aa_WT_GT == 0,kk] >   neg_snp_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 0,kk]))
        aa_neg_dif_snp_percentile_3[[i]][j,kk] <- -(aa_neg_perc_after - aa_neg_perc_before) + (aa_pos_perc_after - aa_pos_perc_before)
        aa_all_before <- sum(WT_models_snp_hal[,kk] <= WT_models_snp_hal[names(neg_snp_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[,kk]))
        aa_all_after <- sum(WT_models_snp_hal[,kk] <= neg_snp_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[,kk]))
        aa_neg_dif_snp_percentile_4[[i]][j,kk] <- aa_all_after - aa_all_before
      #}
    }
  }
}
names(aa_neg_dif_snp_3) <- names(neg_snp_output_hal)
names(aa_neg_dif_snp_percentile_3) <- names(neg_snp_output_hal)
names(aa_neg_dif_snp_percentile_4) <- names(neg_snp_output_hal)

neg_snp_output_hal_all <- do.call(what = rbind, neg_snp_output_hal)
aa_neg_dif_all_snp_3 <- do.call(what = rbind, aa_neg_dif_snp_3)
aa_neg_dif_snp_all_percentile_3 <- do.call(what = rbind, aa_neg_dif_snp_percentile_3)
aa_neg_dif_snp_all_percentile_4 <- do.call(what = rbind, aa_neg_dif_snp_percentile_4)



for(i in 1:length(pos_snp_output_hal)){
  print(i)
  aa_pos_dif_snp_3[[i]] <- matrix(nrow = nrow(pos_snp_output_hal[[i]]), ncol = ncol(pos_snp_output_hal[[i]]))
  colnames(aa_pos_dif_snp_3[[i]]) <- colnames(pos_snp_output_hal[[i]])
  rownames(aa_pos_dif_snp_3[[i]]) <- rownames(pos_snp_output_hal[[i]])
  
  aa_pos_dif_snp_percentile_3[[i]] <- matrix(nrow = nrow(pos_snp_output_hal[[i]]), ncol = ncol(pos_snp_output_hal[[i]]))
  colnames(aa_pos_dif_snp_percentile_3[[i]]) <- colnames(pos_snp_output_hal[[i]])
  rownames(aa_pos_dif_snp_percentile_3[[i]]) <- rownames(pos_snp_output_hal[[i]])
  
  aa_pos_dif_snp_percentile_4[[i]] <- matrix(nrow = nrow(pos_snp_output_hal[[i]]), ncol = ncol(pos_snp_output_hal[[i]]))
  colnames(aa_pos_dif_snp_percentile_4[[i]]) <- colnames(pos_snp_output_hal[[i]])
  rownames(aa_pos_dif_snp_percentile_4[[i]]) <- rownames(pos_snp_output_hal[[i]])
  
  for(j in 1:nrow(pos_snp_output_hal[[i]])){
    aa_pos_dif_snp_3[[i]][j, ] <- (pos_snp_output_hal[[i]][j, ] - WT_models_snp_hal[names(pos_snp_output_hal)[i], ])/ WT_models_snp_hal[names(pos_snp_output_hal)[i], ]
    for(kk in 1:ncol(pos_snp_output_hal[[i]])){
      #if(!is.na(Exp18_WT_model_prediction_Filtered[names(pos_snp_output_hal)[i],kk])){
      aa_pos_perc_before <- sum(WT_models_snp_hal[aa_WT_GT == 1,kk] <= WT_models_snp_hal[names(pos_snp_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 1,kk]))
      aa_neg_perc_before <- sum(WT_models_snp_hal[aa_WT_GT == 0,kk] >  WT_models_snp_hal[names(pos_snp_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 0,kk]))
      aa_pos_perc_after <-  sum(WT_models_snp_hal[aa_WT_GT == 1,kk] <=  pos_snp_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 1,kk]))
      aa_neg_perc_after <-  sum(WT_models_snp_hal[aa_WT_GT == 0,kk] >   pos_snp_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 0,kk]))
      aa_pos_dif_snp_percentile_3[[i]][j,kk] <- -(aa_neg_perc_after - aa_neg_perc_before) + (aa_pos_perc_after - aa_pos_perc_before)
      aa_all_before <- sum(WT_models_snp_hal[,kk] <= WT_models_snp_hal[names(pos_snp_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[,kk]))
      aa_all_after <- sum(WT_models_snp_hal[,kk] <= pos_snp_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[,kk]))
      aa_pos_dif_snp_percentile_4[[i]][j,kk] <- aa_all_after - aa_all_before
      #}
    }
  }
}
names(aa_pos_dif_snp_3) <- names(pos_snp_output_hal)
names(aa_pos_dif_snp_percentile_3) <- names(pos_snp_output_hal)
names(aa_pos_dif_snp_percentile_4) <- names(pos_snp_output_hal)

pos_snp_output_hal_all <- do.call(what = rbind, pos_snp_output_hal)
aa_pos_dif_all_snp_3 <- do.call(what = rbind, aa_pos_dif_snp_3)
aa_pos_dif_snp_all_percentile_3 <- do.call(what = rbind, aa_pos_dif_snp_percentile_3)
aa_pos_dif_snp_all_percentile_4 <- do.call(what = rbind, aa_pos_dif_snp_percentile_4)


aa_pos_neg_out_snp_all_3 <- rbind(pos_snp_output_hal_all, neg_snp_output_hal_all)
aa_all_dif_snp_3 <- rbind(aa_pos_dif_all_snp_3, aa_neg_dif_all_snp_3)
aa_all_dif_snp_perc_3 <- rbind(aa_pos_dif_snp_all_percentile_3, aa_neg_dif_snp_all_percentile_3)
aa_all_dif_snp_perc_4 <- rbind(aa_pos_dif_snp_all_percentile_4, aa_neg_dif_snp_all_percentile_4)


aa_all_dif_snp_names_3 <- unlist(lapply(strsplit(rownames(aa_all_dif_snp_3), split="_"), "[[", 2))

aa_snp_gwas_marker_match_3 <-  aa_snp_gwas_marker[match(aa_all_dif_snp_names_3, names(aa_snp_gwas_marker))]
aa_snp_qtl_marker_pancan_match_3 <- aa_snp_qtl_marker_pancan[match(aa_all_dif_snp_names_3, names(aa_snp_qtl_marker_pancan))]
aa_snp_qtl_marker_gtex_match_3 <- aa_snp_qtl_marker_gtex[match(aa_all_dif_snp_names_3, names(aa_snp_qtl_marker_gtex))]
aa_snp_qtl_marker_gtex_pancan_match_3 <- aa_snp_qtl_marker_gtex_pancan[match(aa_all_dif_snp_names_3, names(aa_snp_qtl_marker_gtex_pancan))]


aa_all_sig_nonsig_3 <- matrix(nrow = 10, ncol = ncol(aa_pos_dif_all_snp_3))
colnames(aa_all_sig_nonsig_3) <- colnames(aa_neg_dif_all_snp_3)

rownames(aa_all_sig_nonsig_3) <- c("sig_all", "sig_qtl_gtex", "nonsig_all", "nonsig_qtl", "pval_pancan", "pval_gtex", "pval_both", "gwas_sig", "sig_qtl_pancan", "sig_qtl_either")



aa_thresh <- c(0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
aa_all_sig_hold_3 <- list()
aa_all_sig_hold_3_uniq <- list()

for(i in 1:length(aa_thresh)){
  print(i)
  aa_all_sig_hold_3[[i]] <- aa_all_sig_nonsig_3
  aa_all_sig_hold_3_uniq[[i]] <- aa_all_sig_nonsig_3
  for(j in 1:(ncol(aa_all_dif_snp_3))){
    aa_all_sig_hold_3[[i]][1,j] <- sum(abs(aa_all_dif_snp_3[,j]) >= aa_thresh[i], na.rm = T)
    aa_m_uniq_up <- unique(aa_all_dif_snp_names_3[which(abs(aa_all_dif_snp_3[,j]) >= aa_thresh[i])])
    aa_all_sig_hold_3_uniq[[i]][1,j] <- length(aa_m_uniq_up)
    aa_all_sig_hold_3[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_3[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_3_uniq[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex[match(aa_m_uniq_up, names(aa_snp_qtl_marker_gtex))])
    aa_all_sig_hold_3[[i]][3,j] <- sum(abs(aa_all_dif_snp_3[,j]) < aa_thresh[i], na.rm = T)
    aa_m_uniq_down <- unique(aa_all_dif_snp_names_3[which(abs(aa_all_dif_snp_3[,j]) < aa_thresh[i])])
    aa_all_sig_hold_3_uniq[[i]][3,j] <- length(aa_m_uniq_down)
    aa_all_sig_hold_3[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_3[,j]) < aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_3_uniq[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex[match(aa_m_uniq_down, names(aa_snp_qtl_marker_gtex))])
    aa_all_sig_hold_3[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan_match_3[abs(aa_all_dif_snp_3[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_3_uniq[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan[match(aa_m_uniq_up, names(aa_snp_qtl_marker_pancan))])
    aa_all_sig_hold_3[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan_match_3[abs(aa_all_dif_snp_3[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_3_uniq[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan[match(aa_m_uniq_up, names(aa_snp_qtl_marker_gtex_pancan))])
    if(sum(aa_snp_qtl_marker_pancan_match_3[abs(aa_all_dif_snp_3[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_3[[i]][5,j] <- phyper(q = sum(aa_snp_qtl_marker_pancan_match_3[abs(aa_all_dif_snp_3[,j]) >= aa_thresh[i]], na.rm = T),
                                            m = sum(aa_snp_qtl_marker_pancan_match_3 == 1),
                                            n = sum(aa_snp_qtl_marker_pancan_match_3 == 0), 
                                            k = aa_all_sig_hold_3[[i]][1,j],
                                            lower.tail = F)
      aa_all_sig_hold_3_uniq[[i]][5,j] <- phyper(q = aa_all_sig_hold_3_uniq[[i]][9,j],
                                                 m = sum(aa_snp_qtl_marker_pancan[names(aa_snp_qtl_marker_pancan) %in% names(aa_snp_qtl_marker_pancan_match_3)] == 1),
                                                 n = sum(aa_snp_qtl_marker_pancan[names(aa_snp_qtl_marker_pancan) %in% names(aa_snp_qtl_marker_pancan_match_3)] == 0), 
                                                 k = aa_all_sig_hold_3_uniq[[i]][1,j],
                                                 lower.tail = F)
    }
    
    if(sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_3[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_3[[i]][6,j] <- phyper(q = sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_3[,j]) >= aa_thresh[i]], na.rm = T),
                                            m = sum(aa_snp_qtl_marker_gtex_match_3 == 1),
                                            n = sum(aa_snp_qtl_marker_gtex_match_3 == 0), 
                                            k = aa_all_sig_hold_3[[i]][1,j],
                                            lower.tail = F)
      aa_all_sig_hold_3_uniq[[i]][6,j] <- phyper(q = aa_all_sig_hold_3_uniq[[i]][2,j],
                                                 m = sum(aa_snp_qtl_marker_gtex[names(aa_snp_qtl_marker_gtex) %in% names(aa_snp_qtl_marker_gtex_match_3)] == 1),
                                                 n = sum(aa_snp_qtl_marker_gtex[names(aa_snp_qtl_marker_gtex) %in% names(aa_snp_qtl_marker_gtex_match_3)] == 0), 
                                                 k = aa_all_sig_hold_3_uniq[[i]][1,j],
                                                 lower.tail = F)
    }
    
    if(sum(aa_snp_qtl_marker_gtex_pancan_match_3[abs(aa_all_dif_snp_3[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_3[[i]][7,j] <- phyper(q = aa_all_sig_hold_3[[i]][10,j],
                                            m = sum(aa_snp_qtl_marker_gtex_pancan_match_3 == 1),
                                            n = sum(aa_snp_qtl_marker_gtex_pancan_match_3 == 0), 
                                            k = aa_all_sig_hold_3[[i]][1,j],
                                            lower.tail = F)
      aa_all_sig_hold_3_uniq[[i]][7,j] <- phyper(q = aa_all_sig_hold_3_uniq[[i]][10,j],
                                                 m = sum(aa_snp_qtl_marker_gtex_pancan[names(aa_snp_qtl_marker_gtex_pancan) %in% names(aa_snp_qtl_marker_gtex_pancan_match_3)] == 1),
                                                 n = sum(aa_snp_qtl_marker_gtex_pancan[names(aa_snp_qtl_marker_gtex_pancan) %in% names(aa_snp_qtl_marker_gtex_pancan_match_3)] == 0), 
                                                 k = aa_all_sig_hold_3_uniq[[i]][1,j],
                                                 lower.tail = F)
    }
    
    
    aa_all_sig_hold_3[[i]][8,j] <- sum(aa_snp_gwas_marker_match_3[abs(aa_all_dif_snp_3[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_3_uniq[[i]][8,j] <- sum(aa_snp_gwas_marker[match(aa_m_uniq_up, names(aa_snp_gwas_marker))], na.rm = T)
    
  }
}
names(aa_all_sig_hold_3) <- aa_thresh
names(aa_all_sig_hold_3_uniq) <- aa_thresh

for(i in 1:length(aa_all_sig_hold_3)){
  print(names(aa_all_sig_hold_3)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_3[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_3[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_3[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_3[[i]][8,] > 0, na.rm=T))
  print("################")
}
sort(aa_all_sig_hold_3$`0.01`[5,])[1:10]
sort(aa_all_sig_hold_3$`0.01`[6,])[1:10]
sort(aa_all_sig_hold_3$`0.01`[7,])[1:10]

for(i in 1:length(aa_all_sig_hold_3_uniq)){
  print(names(aa_all_sig_hold_3_uniq)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_3_uniq[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_3_uniq[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_3_uniq[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_3_uniq[[i]][8,] > 0, na.rm=T))
  print("################")
}

sort(aa_all_sig_hold_3_uniq$`0.01`[5,])[1:10]
sort(aa_all_sig_hold_3_uniq$`0.01`[6,])[1:10]
sort(aa_all_sig_hold_3_uniq$`0.01`[7,])[1:10]

aaxxx <- sort(aa_all_sig_hold_3_uniq$`0.25`[5, ], decreasing = F, index.return=T)$ix
View(aa_all_sig_hold_3_uniq$`0.25`[, aaxxx[1:10]])
##########

aa_all_dif_snp_perc_3
aa_all_dif_snp_perc_4
aa_all_dif_snp_perc_names_3 <- unlist(lapply(strsplit(rownames(aa_all_dif_snp_perc_3), split="_"), "[[", 2))




aa_all_sig_nonsig_perc_3 <- matrix(nrow = 10, ncol = ncol(aa_all_dif_snp_perc_3))
colnames(aa_all_sig_nonsig_perc_3) <- colnames(aa_all_dif_snp_perc_3)

rownames(aa_all_sig_nonsig_perc_3) <- c("sig_all", "sig_qtl_gtex", "nonsig_all",
                                        "nonsig_qtl", "pval_pancan", "pval_gtex",
                                        "pval_both", "gwas_sig", "sig_qtl_pancan", "sig_qtl_either")



aa_thresh <- c(0.001 ,0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
               1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9)
aa_all_sig_hold_perc_3 <- list()
aa_all_sig_hold_perc_3_uniq <- list()
for(i in 1:length(aa_thresh)){
  print(aa_thresh[i])
  aa_all_sig_hold_perc_3[[i]] <- aa_all_sig_nonsig_perc_3
  aa_all_sig_hold_perc_3_uniq[[i]] <- aa_all_sig_nonsig_perc_3
  for(j in 1:(ncol(aa_all_dif_snp_perc_3))){
    aa_all_sig_hold_perc_3[[i]][1,j] <- sum(abs(aa_all_dif_snp_perc_3[,j]) >= aa_thresh[i], na.rm = T)
    aa_m_uniq_up <- unique(aa_all_dif_snp_perc_names_3[which(abs(aa_all_dif_snp_perc_3[,j]) >= aa_thresh[i])])
    aa_all_sig_hold_perc_3_uniq[[i]][1,j] <- length(aa_m_uniq_up)
    aa_all_sig_hold_perc_3[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_perc_3[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_3_uniq[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex[match(aa_m_uniq_up, names(aa_snp_qtl_marker_gtex))])
    aa_all_sig_hold_perc_3[[i]][3,j] <- sum(abs(aa_all_dif_snp_perc_3[,j]) < aa_thresh[i], na.rm = T)
    aa_m_uniq_down <- unique(aa_all_dif_snp_perc_names_3[which(abs(aa_all_dif_snp_perc_3[,j]) < aa_thresh[i])])
    aa_all_sig_hold_perc_3_uniq[[i]][3,j] <- length(aa_m_uniq_down)
    aa_all_sig_hold_perc_3[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_perc_3[,j]) < aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_3_uniq[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex[match(aa_m_uniq_down, names(aa_snp_qtl_marker_gtex))])
    aa_all_sig_hold_perc_3[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan_match_3[abs(aa_all_dif_snp_perc_3[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_3_uniq[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan[match(aa_m_uniq_up, names(aa_snp_qtl_marker_pancan))])
    aa_all_sig_hold_perc_3[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan_match_3[abs(aa_all_dif_snp_perc_3[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_3_uniq[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan[match(aa_m_uniq_up, names(aa_snp_qtl_marker_gtex_pancan))])
    
    if(sum(aa_snp_qtl_marker_pancan_match_3[abs(aa_all_dif_snp_perc_3[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc_3[[i]][5,j] <- phyper(q = aa_all_sig_hold_perc_3[[i]][9,j],
                                                 m = sum(aa_snp_qtl_marker_pancan_match_3 == 1),
                                                 n = sum(aa_snp_qtl_marker_pancan_match_3 == 0), 
                                                 k = aa_all_sig_hold_perc_3[[i]][1,j],
                                                 lower.tail = F)
      aa_all_sig_hold_perc_3_uniq[[i]][5,j] <- phyper(q = aa_all_sig_hold_perc_3_uniq[[i]][9,j],
                                                      m = sum(aa_snp_qtl_marker_pancan[names(aa_snp_qtl_marker_pancan) %in% names(aa_snp_qtl_marker_pancan_match_3)] == 1),
                                                      n = sum(aa_snp_qtl_marker_pancan[names(aa_snp_qtl_marker_pancan) %in% names(aa_snp_qtl_marker_pancan_match_3)] == 0), 
                                                      k = aa_all_sig_hold_perc_3_uniq[[i]][1,j],
                                                      lower.tail = F)
    }
    
    if(sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_perc_3[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc_3[[i]][6,j] <- phyper(q = aa_all_sig_hold_perc_3[[i]][2,j],
                                                 m = sum(aa_snp_qtl_marker_gtex_match_3 == 1),
                                                 n = sum(aa_snp_qtl_marker_gtex_match_3 == 0), 
                                                 k = aa_all_sig_hold_perc_3[[i]][1,j],
                                                 lower.tail = F)
      aa_all_sig_hold_perc_3_uniq[[i]][6,j] <- phyper(q = aa_all_sig_hold_perc_3_uniq[[i]][2,j],
                                                      m = sum(aa_snp_qtl_marker_gtex[names(aa_snp_qtl_marker_gtex) %in% names(aa_snp_qtl_marker_gtex_match_3)] == 1),
                                                      n = sum(aa_snp_qtl_marker_gtex[names(aa_snp_qtl_marker_gtex) %in% names(aa_snp_qtl_marker_gtex_match_3)] == 0), 
                                                      k = aa_all_sig_hold_perc_3_uniq[[i]][1,j],
                                                      lower.tail = F)
      
    }
    
    if(sum(aa_snp_qtl_marker_gtex_pancan_match_3[abs(aa_all_dif_snp_perc_3[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc_3[[i]][7,j] <- phyper(q = aa_all_sig_hold_perc_3[[i]][10,j],
                                                 m = sum(aa_snp_qtl_marker_gtex_pancan_match_3 == 1),
                                                 n = sum(aa_snp_qtl_marker_gtex_pancan_match_3 == 0), 
                                                 k = aa_all_sig_hold_perc_3[[i]][1,j],
                                                 lower.tail = F)
      aa_all_sig_hold_perc_3_uniq[[i]][7,j] <- phyper(q = aa_all_sig_hold_perc_3_uniq[[i]][10,j],
                                                      m = sum(aa_snp_qtl_marker_gtex_pancan[names(aa_snp_qtl_marker_gtex_pancan) %in% names(aa_snp_qtl_marker_gtex_pancan_match_3)] == 1),
                                                      n = sum(aa_snp_qtl_marker_gtex_pancan[names(aa_snp_qtl_marker_gtex_pancan) %in% names(aa_snp_qtl_marker_gtex_pancan_match_3)] == 0), 
                                                      k = aa_all_sig_hold_perc_3_uniq[[i]][1,j],
                                                      lower.tail = F)
    }
    
    
    aa_all_sig_hold_perc_3[[i]][8,j] <- sum(aa_snp_gwas_marker_match_3[abs(aa_all_dif_snp_perc_3[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_3_uniq[[i]][8,j] <- sum(aa_snp_gwas_marker[match(aa_m_uniq_up, names(aa_snp_gwas_marker))], na.rm = T)
    
  }
}
names(aa_all_sig_hold_perc_3) <- aa_thresh
names(aa_all_sig_hold_perc_3_uniq) <- aa_thresh

for(i in 1:length(aa_all_sig_hold_perc_3)){
  print(names(aa_all_sig_hold_perc_3)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_perc_3[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_perc_3[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_perc_3[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_perc_3[[i]][8,] > 0, na.rm=T))
  print("##############")
}

sort(aa_all_sig_hold_perc_3$`0.01`[5,])[1:10]
sort(aa_all_sig_hold_perc_3$`0.01`[6,])[1:10]
sort(aa_all_sig_hold_perc_3$`0.01`[7,])[1:10]


for(i in 1:length(aa_all_sig_hold_perc_3_uniq)){
  print(names(aa_all_sig_hold_perc_3_uniq)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_perc_3_uniq[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_perc_3_uniq[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_perc_3_uniq[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_perc_3_uniq[[i]][8,] > 0, na.rm=T))
  print("##############")
}

sort(aa_all_sig_hold_perc_3_uniq$`0.05`[5,])[1:10]
sort(aa_all_sig_hold_perc_3_uniq$`0.05`[6,])[1:10]
sort(aa_all_sig_hold_perc_3_uniq$`0.05`[7,])[1:10]
#############
aa_all_sig_nonsig_perc_4 <- matrix(nrow = 10, ncol = ncol(aa_all_dif_snp_perc_4))
colnames(aa_all_sig_nonsig_perc_4) <- colnames(aa_all_dif_snp_perc_4)

rownames(aa_all_sig_nonsig_perc_4) <- c("sig_all", "sig_qtl_gtex", "nonsig_all",
                                        "nonsig_qtl", "pval_pancan", "pval_gtex",
                                        "pval_both", "gwas_sig", "sig_qtl_pancan", "sig_qtl_either")



aa_thresh <- c(0.001, 0.005 ,0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1)
aa_all_sig_hold_perc_4 <- list()
aa_all_sig_hold_perc_4_uniq <- list()
for(i in 1:length(aa_thresh)){
  print(aa_thresh[i])
  aa_all_sig_hold_perc_4[[i]] <- aa_all_sig_nonsig_perc_4
  aa_all_sig_hold_perc_4_uniq[[i]] <- aa_all_sig_nonsig_perc_4
  for(j in 1:(ncol(aa_all_dif_snp_perc_4))){
    aa_all_sig_hold_perc_4[[i]][1,j] <- sum(abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i], na.rm = T)
    aa_m_uniq_up <- unique(aa_all_dif_snp_perc_names_3[which(abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i])])
    aa_all_sig_hold_perc_4_uniq[[i]][1,j] <- length(aa_m_uniq_up)
    aa_all_sig_hold_perc_4[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_4_uniq[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex[match(aa_m_uniq_up, names(aa_snp_qtl_marker_gtex))])
    aa_all_sig_hold_perc_4[[i]][3,j] <- sum(abs(aa_all_dif_snp_perc_4[,j]) < aa_thresh[i], na.rm = T)
    aa_m_uniq_down <- unique(aa_all_dif_snp_perc_names_3[which(abs(aa_all_dif_snp_perc_4[,j]) < aa_thresh[i])])
    aa_all_sig_hold_perc_4_uniq[[i]][3,j] <- length(aa_m_uniq_down)
    aa_all_sig_hold_perc_4[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_perc_4[,j]) < aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_4_uniq[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex[match(aa_m_uniq_down, names(aa_snp_qtl_marker_gtex))])
    aa_all_sig_hold_perc_4[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_4_uniq[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan[match(aa_m_uniq_up, names(aa_snp_qtl_marker_pancan))])
    aa_all_sig_hold_perc_4[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_4_uniq[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan[match(aa_m_uniq_up, names(aa_snp_qtl_marker_gtex_pancan))])
    
    if(sum(aa_snp_qtl_marker_pancan_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc_4[[i]][5,j] <- phyper(q = aa_all_sig_hold_perc_4[[i]][9,j],
                                                 m = sum(aa_snp_qtl_marker_pancan_match_3 == 1),
                                                 n = sum(aa_snp_qtl_marker_pancan_match_3 == 0), 
                                                 k = aa_all_sig_hold_perc_4[[i]][1,j],
                                                 lower.tail = F)
      aa_all_sig_hold_perc_4_uniq[[i]][5,j] <- phyper(q = aa_all_sig_hold_perc_4_uniq[[i]][9,j],
                                                      m = sum(aa_snp_qtl_marker_pancan[names(aa_snp_qtl_marker_pancan) %in% names(aa_snp_qtl_marker_pancan_match_3)] == 1),
                                                      n = sum(aa_snp_qtl_marker_pancan[names(aa_snp_qtl_marker_pancan) %in% names(aa_snp_qtl_marker_pancan_match_3)] == 0), 
                                                      k = aa_all_sig_hold_perc_4_uniq[[i]][1,j],
                                                      lower.tail = F)
    }
    
    if(sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc_4[[i]][6,j] <- phyper(q = aa_all_sig_hold_perc_4[[i]][2,j],
                                                 m = sum(aa_snp_qtl_marker_gtex_match_3 == 1),
                                                 n = sum(aa_snp_qtl_marker_gtex_match_3 == 0), 
                                                 k = aa_all_sig_hold_perc_4[[i]][1,j],
                                                 lower.tail = F)
      aa_all_sig_hold_perc_4_uniq[[i]][6,j] <- phyper(q = aa_all_sig_hold_perc_4_uniq[[i]][2,j],
                                                      m = sum(aa_snp_qtl_marker_gtex[names(aa_snp_qtl_marker_gtex) %in% names(aa_snp_qtl_marker_gtex_match_3)] == 1),
                                                      n = sum(aa_snp_qtl_marker_gtex[names(aa_snp_qtl_marker_gtex) %in% names(aa_snp_qtl_marker_gtex_match_3)] == 0), 
                                                      k = aa_all_sig_hold_perc_4_uniq[[i]][1,j],
                                                      lower.tail = F)
      
    }
    
    if(sum(aa_snp_qtl_marker_gtex_pancan_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]], na.rm = T) > 0){
      aa_all_sig_hold_perc_4[[i]][7,j] <- phyper(q = aa_all_sig_hold_perc_4[[i]][10,j],
                                                 m = sum(aa_snp_qtl_marker_gtex_pancan_match_3 == 1),
                                                 n = sum(aa_snp_qtl_marker_gtex_pancan_match_3 == 0), 
                                                 k = aa_all_sig_hold_perc_4[[i]][1,j],
                                                 lower.tail = F)
      aa_all_sig_hold_perc_4_uniq[[i]][7,j] <- phyper(q = aa_all_sig_hold_perc_4_uniq[[i]][10,j],
                                                      m = sum(aa_snp_qtl_marker_gtex_pancan[names(aa_snp_qtl_marker_gtex_pancan) %in% names(aa_snp_qtl_marker_gtex_pancan_match_3)] == 1),
                                                      n = sum(aa_snp_qtl_marker_gtex_pancan[names(aa_snp_qtl_marker_gtex_pancan) %in% names(aa_snp_qtl_marker_gtex_pancan_match_3)] == 0), 
                                                      k = aa_all_sig_hold_perc_4_uniq[[i]][1,j],
                                                      lower.tail = F)
    }
    
    
    aa_all_sig_hold_perc_4[[i]][8,j] <- sum(aa_snp_gwas_marker_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_4_uniq[[i]][8,j] <- sum(aa_snp_gwas_marker[match(aa_m_uniq_up, names(aa_snp_gwas_marker))], na.rm = T)
    
  }
}
names(aa_all_sig_hold_perc_4) <- aa_thresh
names(aa_all_sig_hold_perc_4_uniq) <- aa_thresh

for(i in 1:length(aa_all_sig_hold_perc_4)){
  print(names(aa_all_sig_hold_perc_4)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_perc_4[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_perc_4[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_perc_4[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_perc_4[[i]][8,] > 0, na.rm=T))
  print("##############")
}
sort(aa_all_sig_hold_perc_4$`0.05`[7,])[1:10]

for(i in 1:length(aa_all_sig_hold_perc_4_uniq)){
  print(names(aa_all_sig_hold_perc_4_uniq)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_perc_4_uniq[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_perc_4_uniq[[i]][5,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_perc_4_uniq[[i]][7,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_perc_4_uniq[[i]][8,] > 0, na.rm=T))
  print("##############")
}

sort(aa_all_sig_hold_perc_4_uniq$`0.005`[7,])[1:10]

aa_all_sig_hold_perc_4_uniq$`0.005`[, "par_4493"]



# looking at the correlation of various scores

aacor_row <- numeric(length = nrow(aa_all_dif_snp_perc_4))
aacor_col <- numeric(length = ncol(aa_all_dif_snp_perc_4))

for(i in 1:nrow(aa_all_dif_snp_perc_4)){
  aacor_row[i] <- cor(aa_all_dif_snp_perc_4[i,], aa_all_dif_snp_perc_3[i,])
}
par(mfrow = c(1,1), mar = c(4,4,4,4))
hist(aacor_row)


for(i in 1:ncol(aa_all_dif_snp_perc_4)){
  aacor_col[i] <- cor(aa_all_dif_snp_perc_4[,i], aa_all_dif_snp_perc_3[,i])
}
hist(aacor_col)

######################
# evaluate a different hypergeom test: what percent of models predict any change (> 0.001)
aa_all_dif_snp_perc_4
aa_any_change <- abs(aa_all_dif_snp_perc_4) > 0.001
aa_any_change <- rowSums(aa_any_change)
aa_any_change <- aa_any_change/244
hist(aa_any_change, breaks = 150)

aa_thresh <- seq(0,0.99,0.05)
aa_thresh_pval_gtex <- numeric(length = length(aa_thresh))
aa_thresh_pval_pancan <- numeric(length = length(aa_thresh))
aa_thresh_pval_both <- numeric(length = length(aa_thresh))
names(aa_thresh_pval_gtex) <- aa_thresh
names(aa_thresh_pval_pancan) <- aa_thresh
names(aa_thresh_pval_both) <- aa_thresh

for(i in 1:length(aa_thresh)){
  aa_thresh_pval_gtex[i] <- phyper(q = sum(aa_snp_qtl_marker_gtex_match_3[aa_any_change > aa_thresh[i]]),
                                   m = sum(aa_snp_qtl_marker_gtex_match_3), 
                                   n = sum(aa_snp_qtl_marker_gtex_match_3 == 0),
                                   k = sum(aa_any_change > aa_thresh[i]), lower.tail = F)
  aa_thresh_pval_pancan[i] <- phyper(q = sum(aa_snp_qtl_marker_pancan_match_3[aa_any_change > aa_thresh[i]]),
                                   m = sum(aa_snp_qtl_marker_pancan_match_3), 
                                   n = sum(aa_snp_qtl_marker_pancan_match_3 == 0),
                                   k = sum(aa_any_change > aa_thresh[i]), lower.tail = F)
  aa_thresh_pval_both[i] <- phyper(q = sum(aa_snp_qtl_marker_gtex_pancan_match_3[aa_any_change > aa_thresh[i]]),
                                   m = sum(aa_snp_qtl_marker_gtex_pancan_match_3), 
                                   n = sum(aa_snp_qtl_marker_gtex_pancan_match_3 == 0),
                                   k = sum(aa_any_change > aa_thresh[i]), lower.tail = F)
}


aa_thresh_pval_gtex
aa_thresh_pval_pancan
aa_thresh_pval_both

sum(aa_any_change > 0.75 & aa_snp_qtl_marker_gtex_pancan_match_3 == 1)





aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_snp <- readPeakFile(aa_snp_Add, as = "data.frame")
aa_snp$V10 <- as.character(levels(aa_snp$V10)[as.numeric(aa_snp$V10)])
aa_snp$V9 <- as.character(levels(aa_snp$V9)[as.numeric(aa_snp$V9)])
aa_snp$V8 <- as.character(levels(aa_snp$V8)[as.numeric(aa_snp$V8)])
aa_snp$V1 <- as.character(levels(aa_snp$V1)[as.numeric(aa_snp$V1)])
aa_snp$V2 <- aa_snp$V2 - 2
aa_snp$V3 <- aa_snp$V3 - 1





aaa <- names(aa_any_change)[aa_any_change > 0.75 & aa_snp_qtl_marker_gtex_pancan_match_3 == 1]
aa_1 <- unlist(lapply(strsplit(aaa, split="_"), "[[", 2))
aa_2 <- unlist(lapply(strsplit(aaa, split="_"), "[[", 3))


aa_eqtl_chr <- character(length(aaa))
aa_eqtl_st <- numeric(length(aaa))
aa_eqtl_ref <- numeric(length(aaa))
aa_eqtl_alt <- numeric(length(aaa))

for(i in 1:length(aaa)){
  aa_eqtl_chr[i] <- aa_snp$V1[aa_snp$V10 == aa_1[i]]
  aa_eqtl_st[i] <- aa_snp$V2[aa_snp$V10 == aa_1[i]]
  aa_eqtl_ref[i] <- aa_snp$V8[aa_snp$V10 == aa_1[i]]
  aaaaalt <- unlist(strsplit(aa_snp$V9[aa_snp$V10 == aa_1[i]], split = ","))
  aa_eqtl_alt[i] <- aaaaalt[as.numeric(aa_2[i])]
}



aamin_LLR_low <- aa_exp18_minLLR - 0.1

aa_snp_inves_3 <- list()
for(i in 1:length(aa_eqtl_alt)){
  print(i)
  aa_snp_inves_3[[i]] <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                                          my_motifs = TF.motifs.Expanded_pseudo_exp12_t,
                                          snp_chr = aa_eqtl_chr[i],
                                          snp_start = aa_eqtl_st[i],
                                          snp_ref = aa_eqtl_ref[i],
                                          snp_alt = aa_eqtl_alt[i],
                                          min_LLR = aamin_LLR_low[aa_names])
  
}

names(aa_snp_inves_3) <- aaa

aa_snp_inves_3_mat <- do.call(rbind, aa_snp_inves_3)
aa_snp_inves_3_mat <- aa_snp_inves_3_mat[aa_snp_inves_3_mat[, 4] != 0,]
#View(aa_snp_inves_3_mat)

aan <- unlist(lapply(strsplit(rownames(aa_snp_inves_3_mat), "\\."), "[[", 1))
aa_snp_inves_3_mat <- aa_snp_inves_3_mat[aan %in%  rownames(aa_all_dif_snp_perc_4),]
aan <- unlist(lapply(strsplit(rownames(aa_snp_inves_3_mat), "\\."), "[[", 1))

aamacx <- numeric(length(aan) )
aamean <- numeric(length(aan) )
for(i in 1:length(aamacx)){
  aamacx[i] <- max(abs(aa_all_dif_snp_perc_4[aan[i], ])) * sign(aa_all_dif_snp_perc_4[aan[i], which.max(abs(aa_all_dif_snp_perc_4[aan[i], ]))])
  aamean[i] <- mean(abs(aa_all_dif_snp_perc_4[aan[i], abs(aa_all_dif_snp_perc_4[aan[i],]) > 0.001 ]))
}
 
aa_snp_inves_3_mat <- cbind(aa_snp_inves_3_mat, aamacx)
aa_snp_inves_3_mat <- cbind(aa_snp_inves_3_mat, aamean)
colnames(aa_snp_inves_3_mat)[c(11,12)] <- c("max_percentile_change", "mean_percentile_change")


aa_dpn <- character(length(aan))
names(aa_dpn) <- aan
for(j in 1:length(aan)){
  for(i in 1:length(pos_snp_output_hal)){
    if(aan[j] %in% rownames(pos_snp_output_hal[[i]])){
      aa_dpn[j] <- names(pos_snp_output_hal)[i]
      print("########")
      print(aan[j])
      print(names(pos_snp_output_hal)[i])
      print("########")
    }
  }
  for(i in 1:length(neg_snp_output_hal)){
    if(aan[j] %in% rownames(neg_snp_output_hal[[i]])){
      aa_dpn[j] <- names(neg_snp_output_hal)[i]
      print("########")
      print(aan[j])
      print(names(neg_snp_output_hal)[i])
      print("########")
    }
  }
}

aa_snp_inves_3_mat <- cbind(aa_snp_inves_3_mat, aa_dpn)
aa_snp_inves_3_mat$TF[aa_snp_inves_3_mat$TF == "ESR1_2"] <- "ESR1"


aa_qtl_chip_ov <- Enhancer.ReMapchip.Overlap.byEnhancer$OverlapMat[aa_dpn,]
aa_qtl_chip_in <- Enhancer.ReMapchip.Overlap.byEnhancer$IntMat[aa_dpn,]

aa_qtl_chip_ovl_evid <- array(dim = nrow(aa_qtl_chip_ov))
aa_qtl_chip_int_evid <- array(dim =nrow(aa_qtl_chip_in))
for(i in 1:length(aa_qtl_chip_ovl_evid)){
  if(aa_snp_inves_3_mat$TF[i] %in% colnames(aa_qtl_chip_ov)){
    aa_qtl_chip_ovl_evid[i] <- aa_qtl_chip_ov[i, aa_snp_inves_3_mat$TF[i]]
    aa_qtl_chip_int_evid[i] <- aa_qtl_chip_in[i, aa_snp_inves_3_mat$TF[i]]
  }
}

aa_snp_inves_3_mat <- cbind(aa_snp_inves_3_mat, cbind(aa_qtl_chip_ovl_evid, aa_qtl_chip_int_evid))
colnames(aa_snp_inves_3_mat)[c(13:15)] <- c("enhancer", "chip_overlap_evidence", "chip_interaction_evidence")



names(enhancer_gene_ovlap_interact$IntList[[1]]) <- names(aa_pos_neg)
names(enhancer_gene_ovlap_interact$OverlapList[[1]]) <- names(aa_pos_neg)
names(enhancer_enhancer_interact$IntList[[1]]) <- names(aa_pos_neg)
aa_ovl_gene <- character(length = nrow(aa_snp_inves_3_mat))
aa_int_gene <- character(length = nrow(aa_snp_inves_3_mat))
aa_int_enha <- character(length = nrow(aa_snp_inves_3_mat))
for(i in 1: nrow(aa_snp_inves_3_mat)){
  aa_ovl_gene[i] <- paste(enhancer_gene_ovlap_interact$OverlapList[[1]][[aa_snp_inves_3_mat$enhancer[i]]], collapse = "__")
  aa_int_gene[i] <- paste(enhancer_gene_ovlap_interact$IntList[[1]][[aa_snp_inves_3_mat$enhancer[i]]], collapse = "__")
  aa_int_enha[i] <- paste(enhancer_enhancer_interact$IntList[[1]][[aa_snp_inves_3_mat$enhancer[i]]], collapse = "__")
}
aa_snp_inves_3_mat <- cbind(aa_snp_inves_3_mat, cbind(cbind(aa_int_gene, aa_ovl_gene), aa_int_enha))
colnames(aa_snp_inves_3_mat)[16:18] <- c("interacting_Gene", "overlapping_Gene", "interacting_enhancer")

write.csv(x = aa_snp_inves_3_mat, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/eqtl_TF_75pModels.csv")
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# changing the universe of SNPs to all SNPs instead of the ones overlapping our enhancers
head(ER_gtex_eqtl_sig)
aa_allsnp_file="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps.bed"
aa_eqtl_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/sig_eqtl.bed"

aa_qtl <- readPeakFile(aa_eqtl_Add, as = "GRanges")
start(aa_qtl) <- start(aa_qtl) - 2
#aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_snp_all <- readPeakFile(aa_allsnp_file, as = "GRanges")
start(aa_snp_all) <- start(aa_snp_all) -1
aa_overlap_holder <- findOverlaps(query = aa_snp_all, subject = aa_qtl)
#ol_enhancers_index <- unique(overlap_holder@from)
aa_snp_qtl_marker_all <- numeric(length(aa_snp_all))
aa_snp_qtl_marker_all[aa_overlap_holder@from] <- 1
#aa_snp$V10 <- as.character(levels(aa_snp$V10)[as.numeric(aa_snp$V10)])
names(aa_snp_qtl_marker) <- aa_snp_all$V8
aa_snp_gwas_marker <- numeric(length(aa_snp))
names(aa_snp_gwas_marker) <-  aa_snp$V10
aa_snp_gwas_marker[aa_GWAS_SNPS] <- 1



ranges(ER_pancan_eqtl_cis_gr38) <- ranges(ER_pancan_eqtl_cis_gr38) + 1
aa_overlap_holder_paanc <- findOverlaps(query = aa_snp_all, subject = ER_pancan_eqtl_cis_gr38)
aa_snp_qtl_marker_pancan_cis <- numeric(length(aa_snp_all))
aa_snp_qtl_marker_pancan_cis[aa_overlap_holder_paanc@from] <- 1
#aa_snp$V10 <- as.character(levels(aa_snp$V10)[as.numeric(aa_snp$V10)])
names(aa_snp_qtl_marker_pancan_cis) <- aa_snp$V10

ranges(ER_pancan_eqtl_trans_gr38) <- ranges(ER_pancan_eqtl_trans_gr38) + 1
aa_overlap_holder_paanctrans <- findOverlaps(query = aa_snp_all,
                                  subject = ER_pancan_eqtl_trans_gr38)
aa_snp_qtl_marker_pancan_trans <- numeric(length(aa_snp_all))
aa_snp_qtl_marker_pancan_trans[aa_overlap_holder_paanctrans@from] <- 1
names(aa_snp_qtl_marker_pancan_trans) <- aa_snp$V10

aa_snp_qtl_marker_pancan_new <- aa_snp_qtl_marker_pancan_trans + aa_snp_qtl_marker_pancan_cis
aa_snp_qtl_marker_pancan_new[aa_snp_qtl_marker_pancan_new > 0] <- 1
aa_snp_qtl_marker_gtex_pancan_new <- aa_snp_qtl_marker_all + aa_snp_qtl_marker_pancan_new


######## overlap with GWAS
ranges(GWAS_breast_Cancer_gr38) <- ranges(GWAS_breast_Cancer_gr38) + 1
aa_overlap_holder_gwas_new <- findOverlaps(query = aa_snp_all,
                                  subject = GWAS_breast_Cancer_gr38)
aa_snp_gwas_marker <- numeric(length(aa_snp_all))
aa_snp_gwas_marker[aa_overlap_holder_gwas_new@from] <- 1
names(aa_snp_gwas_marker) <-  aa_snp$V10





aa_all_sig_nonsig_perc_4 <- matrix(nrow = 12, ncol = ncol(aa_all_dif_snp_perc_4))
colnames(aa_all_sig_nonsig_perc_4) <- colnames(aa_all_dif_snp_perc_4)

rownames(aa_all_sig_nonsig_perc_4) <- c("sig_all", "sig_qtl_gtex", "nonsig_all",
                                        "nonsig_qtl", "pval_pancan", "pval_gtex",
                                        "pval_either", "gwas_sig", "sig_qtl_pancan",
                                        "sig_qtl_either", "sig_qtl_both", "pval_both")



aa_thresh <- c(0.001, 0.005 ,0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1)
aa_all_sig_hold_perc_4_univ2 <- list()
for(i in 1:length(aa_thresh)){
  print(aa_thresh[i])
  aa_all_sig_hold_perc_4_univ2[[i]] <- aa_all_sig_nonsig_perc_4
  for(j in 1:(ncol(aa_all_dif_snp_perc_4))){
    aa_all_sig_hold_perc_4_univ2[[i]][1,j] <- sum(abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i], na.rm = T)
    aa_all_sig_hold_perc_4_univ2[[i]][2,j] <- sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_4_univ2[[i]][3,j] <- sum(abs(aa_all_dif_snp_perc_4[,j]) < aa_thresh[i], na.rm = T)
    aa_all_sig_hold_perc_4_univ2[[i]][4,j] <- sum(aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_perc_4[,j]) < aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_4_univ2[[i]][9,j] <- sum(aa_snp_qtl_marker_pancan_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_4_univ2[[i]][10,j] <- sum(aa_snp_qtl_marker_gtex_pancan_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]], na.rm = T)
    aa_all_sig_hold_perc_4_univ2[[i]][11,j] <- sum((aa_snp_qtl_marker_pancan_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]]) & (aa_snp_qtl_marker_gtex_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]]), na.rm = T)

    if(aa_all_sig_hold_perc_4_univ2[[i]][9,j] > 0){
      aa_all_sig_hold_perc_4_univ2[[i]][5,j] <- phyper(q = aa_all_sig_hold_perc_4_univ2[[i]][9,j],
                                                 m = 363422,
                                                 n = 36939556, 
                                                 k = aa_all_sig_hold_perc_4_univ2[[i]][1,j],
                                                 lower.tail = F)
    }
    
    if(aa_all_sig_hold_perc_4_univ2[[i]][2,j] > 0){
      aa_all_sig_hold_perc_4_univ2[[i]][6,j] <- phyper(q = aa_all_sig_hold_perc_4_univ2[[i]][2,j],
                                                 m = 935185,
                                                 n = 36367793, 
                                                 k = aa_all_sig_hold_perc_4_univ2[[i]][1,j],
                                                 lower.tail = F)
    }
    
    if(aa_all_sig_hold_perc_4_univ2[[i]][10,j] > 0){
      aa_all_sig_hold_perc_4_univ2[[i]][7,j] <- phyper(q = aa_all_sig_hold_perc_4_univ2[[i]][10,j],
                                                 m = 1113895,
                                                 n = 36189083, 
                                                 k = aa_all_sig_hold_perc_4_univ2[[i]][1,j],
                                                 lower.tail = F)
    }
    if(aa_all_sig_hold_perc_4_univ2[[i]][11,j] > 0){
      aa_all_sig_hold_perc_4_univ2[[i]][12,j] <- phyper(q = aa_all_sig_hold_perc_4_univ2[[i]][11,j],
                                                       m = 184712,
                                                       n = 37118266, 
                                                       k = aa_all_sig_hold_perc_4_univ2[[i]][1,j],
                                                       lower.tail = F)
    }
    
    
    aa_all_sig_hold_perc_4_univ2[[i]][8,j] <- sum(aa_snp_gwas_marker_match_3[abs(aa_all_dif_snp_perc_4[,j]) >= aa_thresh[i]], na.rm = T)

  }
}
names(aa_all_sig_hold_perc_4_univ2) <- aa_thresh

for(i in 1:length(aa_all_sig_hold_perc_4_univ2)){
  print(names(aa_all_sig_hold_perc_4_univ2)[i])
  print("gtex")
  print(sum((aa_all_sig_hold_perc_4_univ2[[i]][6,] < 0.05), na.rm =T))
  print("pancan")
  print(sum((aa_all_sig_hold_perc_4_univ2[[i]][5,] < 0.05), na.rm =T))
  print("either")
  print(sum((aa_all_sig_hold_perc_4_univ2[[i]][7,] < 0.05), na.rm =T))
  print("both")
  print(sum((aa_all_sig_hold_perc_4_univ2[[i]][12,] < 0.05), na.rm =T))
  print("gwas")
  print(sum(aa_all_sig_hold_perc_4_univ2[[i]][8,] > 0, na.rm=T))
  print("##############")
}
i<-1
names(aa_all_sig_hold_perc_4_univ2)[i]
sort(aa_all_sig_hold_perc_4_univ2[[i]][5,])[1:10]
sort(aa_all_sig_hold_perc_4_univ2[[i]][6,])[1:10]
sort(aa_all_sig_hold_perc_4_univ2[[i]][7,])[1:10]
sort(aa_all_sig_hold_perc_4_univ2[[i]][12,])[1:10]

#### evaluate pvalue for percent models that agree there is a change

aa_all_dif_snp_perc_4
aa_any_change <- abs(aa_all_dif_snp_perc_4) > 0.001
aa_any_change <- rowSums(aa_any_change)
aa_any_change <- aa_any_change/244
hist(aa_any_change, breaks = 150)

aa_thresh <- seq(0,0.99,0.05)
aa_thresh_pval_gtex <- numeric(length = length(aa_thresh))
aa_thresh_pval_pancan <- numeric(length = length(aa_thresh))
aa_thresh_pval_either <- numeric(length = length(aa_thresh))
aa_thresh_pval_both <- numeric(length = length(aa_thresh))

names(aa_thresh_pval_gtex) <- aa_thresh
names(aa_thresh_pval_pancan) <- aa_thresh
names(aa_thresh_pval_both) <- aa_thresh
names(aa_thresh_pval_either) <- aa_thresh

for(i in 1:length(aa_thresh)){
  aa_thresh_pval_gtex[i] <- phyper(q = sum(aa_snp_qtl_marker_gtex_match_3[aa_any_change > aa_thresh[i]]),
                                   m = 935185, 
                                   n = 36367793,
                                   k = sum(aa_any_change > aa_thresh[i]), lower.tail = F)
  aa_thresh_pval_pancan[i] <- phyper(q = sum(aa_snp_qtl_marker_pancan_match_3[aa_any_change > aa_thresh[i]]),
                                     m = 363422, 
                                     n = 36939556,
                                     k = sum(aa_any_change > aa_thresh[i]), lower.tail = F)
  aa_thresh_pval_either[i] <- phyper(q = sum(aa_snp_qtl_marker_gtex_pancan_match_3[aa_any_change > aa_thresh[i]]),
                                   m = 1113895, 
                                   n = 36189083,
                                   k = sum(aa_any_change > aa_thresh[i]), lower.tail = F)
  aa_thresh_pval_both[i] <- phyper(q = sum((aa_snp_qtl_marker_gtex_match_3[aa_any_change > aa_thresh[i]]) & (aa_snp_qtl_marker_pancan_match_3[aa_any_change > aa_thresh[i]])),
                                   m = 184712, 
                                   n = 37118266,
                                   k = sum(aa_any_change > aa_thresh[i]), lower.tail = F)
}



barplot(-log10(aa_thresh_pval_gtex), las = 2)
barplot(-log10(aa_thresh_pval_pancan), las = 2)
barplot(-log10(aa_thresh_pval_either), las = 2)
barplot(-log10(aa_thresh_pval_both), las = 2)

sum(aa_any_change > 0.75 & aa_snp_qtl_marker_gtex_pancan_match_3 == 1)
#########################
#write an excell file for all SNPs that are predicted to be effective --- add a column to show if it's a gtex, or pancan eqtl
sum(aa_any_change > 0.75)





aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_snp <- readPeakFile(aa_snp_Add, as = "data.frame")
aa_snp$V10 <- as.character(levels(aa_snp$V10)[as.numeric(aa_snp$V10)])
aa_snp$V9 <- as.character(levels(aa_snp$V9)[as.numeric(aa_snp$V9)])
aa_snp$V8 <- as.character(levels(aa_snp$V8)[as.numeric(aa_snp$V8)])
aa_snp$V1 <- as.character(levels(aa_snp$V1)[as.numeric(aa_snp$V1)])
aa_snp$V2 <- aa_snp$V2 - 2
aa_snp$V3 <- aa_snp$V3 - 1





aaa <- names(aa_any_change)[aa_any_change > 0.75]
aaa <- aaa[!duplicated(aaa)]
aa_1 <- unlist(lapply(strsplit(aaa, split="_"), "[[", 2))
aa_2 <- unlist(lapply(strsplit(aaa, split="_"), "[[", 3))


aa_eqtl_chr <- character(length(aaa))
aa_eqtl_st <- numeric(length(aaa))
aa_eqtl_ref <- numeric(length(aaa))
aa_eqtl_alt <- numeric(length(aaa))

for(i in 1:length(aaa)){
  aa_eqtl_chr[i] <- unique(aa_snp$V1[aa_snp$V10 == aa_1[i]])
  aa_eqtl_st[i] <- unique(aa_snp$V2[aa_snp$V10 == aa_1[i]])
  aa_eqtl_ref[i] <- unique(aa_snp$V8[aa_snp$V10 == aa_1[i]])
  aaaaalt <- unlist(strsplit(unique(aa_snp$V9[aa_snp$V10 == aa_1[i]]), split = ","))
  aa_eqtl_alt[i] <- aaaaalt[as.numeric(aa_2[i])]
}



aamin_LLR_low <- aa_exp18_minLLR - 0.1

aa_snp_inves_4 <- list()
for(i in 1:length(aa_eqtl_alt)){
  print(i)
  aa_snp_inves_4[[i]] <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                                          my_motifs = TF.motifs.Expanded_pseudo_exp12_t,
                                          snp_chr = aa_eqtl_chr[i],
                                          snp_start = aa_eqtl_st[i],
                                          snp_ref = aa_eqtl_ref[i],
                                          snp_alt = aa_eqtl_alt[i],
                                          min_LLR = aamin_LLR_low[aa_names])
  
}

names(aa_snp_inves_4) <- aaa

aa_snp_inves_4_mat <- do.call(rbind, aa_snp_inves_4)
aa_snp_inves_4_mat <- aa_snp_inves_4_mat[aa_snp_inves_4_mat[, 4] != 0,]
#View(aa_snp_inves_3_mat)

aan <- unlist(lapply(strsplit(rownames(aa_snp_inves_4_mat), "\\."), "[[", 1))
aa_snp_inves_4_mat <- aa_snp_inves_4_mat[aan %in%  rownames(aa_all_dif_snp_perc_4),]
aan <- unlist(lapply(strsplit(rownames(aa_snp_inves_4_mat), "\\."), "[[", 1))

aamacx <- numeric(length(aan) )
aamean <- numeric(length(aan) )
for(i in 1:length(aamacx)){
  aamacx[i] <- max(abs(aa_all_dif_snp_perc_4[aan[i], ])) * sign(aa_all_dif_snp_perc_4[aan[i], which.max(abs(aa_all_dif_snp_perc_4[aan[i], ]))])
  aamean[i] <- mean(abs(aa_all_dif_snp_perc_4[aan[i], abs(aa_all_dif_snp_perc_4[aan[i],]) > 0.001 ]))
}

aa_snp_inves_4_mat <- cbind(aa_snp_inves_4_mat, aamacx)
aa_snp_inves_4_mat <- cbind(aa_snp_inves_4_mat, aamean)
colnames(aa_snp_inves_4_mat)[c(11,12)] <- c("max_percentile_change", "mean_percentile_change")


aa_dpn <- character(length(aan))
names(aa_dpn) <- aan
for(j in 1:length(aan)){
  for(i in 1:length(pos_snp_output_hal)){
    if(aan[j] %in% rownames(pos_snp_output_hal[[i]])){
      aa_dpn[j] <- names(pos_snp_output_hal)[i]
      print("########")
      print(aan[j])
      print(names(pos_snp_output_hal)[i])
      print("########")
    }
  }
  for(i in 1:length(neg_snp_output_hal)){
    if(aan[j] %in% rownames(neg_snp_output_hal[[i]])){
      aa_dpn[j] <- names(neg_snp_output_hal)[i]
      print("########")
      print(aan[j])
      print(names(neg_snp_output_hal)[i])
      print("########")
    }
  }
}

aa_snp_inves_4_mat <- cbind(aa_snp_inves_4_mat, aa_dpn)
aa_snp_inves_4_mat$TF[aa_snp_inves_4_mat$TF == "ESR1_2"] <- "ESR1"


aa_qtl_chip_ov <- Enhancer.ReMapchip.Overlap.byEnhancer$OverlapMat[aa_dpn,]
aa_qtl_chip_in <- Enhancer.ReMapchip.Overlap.byEnhancer$IntMat[aa_dpn,]

aa_qtl_chip_ovl_evid <- array(dim = nrow(aa_qtl_chip_ov))
aa_qtl_chip_int_evid <- array(dim =nrow(aa_qtl_chip_in))
for(i in 1:length(aa_qtl_chip_ovl_evid)){
  if(aa_snp_inves_4_mat$TF[i] %in% colnames(aa_qtl_chip_ov)){
    aa_qtl_chip_ovl_evid[i] <- aa_qtl_chip_ov[i, aa_snp_inves_4_mat$TF[i]]
    aa_qtl_chip_int_evid[i] <- aa_qtl_chip_in[i, aa_snp_inves_4_mat$TF[i]]
  }
}

aa_snp_inves_4_mat <- cbind(aa_snp_inves_4_mat, cbind(aa_qtl_chip_ovl_evid, aa_qtl_chip_int_evid))
colnames(aa_snp_inves_4_mat)[c(13:15)] <- c("enhancer", "chip_overlap_evidence", "chip_interaction_evidence")



#names(enhancer_gene_ovlap_interact$IntList[[1]]) <- names(aa_pos_neg)
#names(enhancer_gene_ovlap_interact$OverlapList[[1]]) <- names(aa_pos_neg)
names(enhancer_enhancer_interact$IntList[[1]]) <- names(aa_pos_neg)
names(enhancer_enhancer_interact$OverlapList[[1]]) <- names(aa_pos_neg)

aa_ovl_gene <- character(length = nrow(aa_snp_inves_4_mat))
aa_int_gene <- character(length = nrow(aa_snp_inves_4_mat))
aa_int_enha <- character(length = nrow(aa_snp_inves_4_mat))
aa_ovl_enha <- character(length = nrow(aa_snp_inves_4_mat))
for(i in 1: nrow(aa_snp_inves_4_mat)){
  aa_ovl_gene[i] <- paste(enhancer_gene_ovlap_interact$OverlapList[[1]][[aa_snp_inves_4_mat$enhancer[i]]], collapse = "__")
  aa_int_gene[i] <- paste(enhancer_gene_ovlap_interact$IntList[[1]][[aa_snp_inves_4_mat$enhancer[i]]], collapse = "__")
  aa_int_enha[i] <- paste(enhancer_enhancer_interact$IntList[[1]][[aa_snp_inves_4_mat$enhancer[i]]], collapse = "__")
  aa_ovl_enha[i] <- paste(enhancer_enhancer_interact$OverlapList[[1]][[aa_snp_inves_4_mat$enhancer[i]]], collapse = "__")
}
aa_snp_inves_4_mat <- cbind(aa_snp_inves_4_mat, cbind(cbind(aa_int_gene, aa_ovl_gene), cbind(aa_int_enha, aa_ovl_enha)))
colnames(aa_snp_inves_4_mat)[16:19] <- c("interacting_Gene", "overlapping_Gene", "interacting_enhancer", "overlapping_enhancer")

aars <- unlist(lapply(strsplit(unlist(lapply(strsplit(rownames(aa_snp_inves_4_mat), "\\."), "[[", 1)), "_"), "[[", 2))
aars_gtex <- aa_snp_qtl_marker_gtex_match_3[match(aars, names(aa_snp_qtl_marker_gtex_match_3))]
aars_pancan <- aa_snp_qtl_marker_pancan_match_3[match(aars, names(aa_snp_qtl_marker_pancan_match_3))]
aars_GWAS <- numeric(length(aars))
names(aars_GWAS) <- aars
aars_GWAS[aars %in% aa_GWAS_SNPS] <- 1
aa_snp_inves_4_mat <- cbind(aa_snp_inves_4_mat, aars_gtex, aars_pancan, aars_GWAS)
colnames(aa_snp_inves_4_mat)[20:22] <- c("GTEX_eqtl", "PanCan_eqtl","BC_GWAS")

write.csv(x = aa_snp_inves_4_mat, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/eqtl_TF_75pModels_all.csv")
##############################################################################################################################
##############################################################################################################################
# do hypergeom only for enhancers --> only one test:

length(unique(unlist(lapply(strsplit(names(aa_any_change), "_"), "[[", 2))))

aa_pval_enh_gtex <- phyper(q = sum(aa_snp_qtl_marker_gtex_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_match_3))]),
                           m = 935185, 
                           n = 36367793,
                           k = 26382,
                           lower.tail = F)
aa_enhancer_gtex_FPR <- sum(aa_snp_qtl_marker_gtex_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_match_3))] == 0)/36367793
aa_enhancer_gtex_FNR <- (935185 - sum(aa_snp_qtl_marker_gtex_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_match_3))] == 1))/935185
aa_enhancer_gtex_TPR <- sum(aa_snp_qtl_marker_gtex_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_match_3))])/935185
aa_enhancer_gtex_TNR <- (36367793 - sum(aa_snp_qtl_marker_gtex_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_match_3))] == 0))/36367793

aa_pval_enh_pancan <- phyper(q = sum(aa_snp_qtl_marker_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_pancan_match_3))]),
                             m = 363422, 
                             n = 36939556,
                             k = 26382,
                             lower.tail = F)
aa_enhancer_pancan_FPR <- sum(aa_snp_qtl_marker_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_pancan_match_3))] == 0)/36939556
aa_enhancer_pancan_FNR <- (363422 - sum(aa_snp_qtl_marker_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_pancan_match_3))] == 1))/363422
aa_enhancer_pancan_TPR <- sum(aa_snp_qtl_marker_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_pancan_match_3))])/363422
aa_enhancer_pancan_TNR <- (36939556 - sum(aa_snp_qtl_marker_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_pancan_match_3))] == 0))/36939556

aa_pval_enh_either <- phyper(q = sum(aa_snp_qtl_marker_gtex_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_pancan_match_3))]),
                             m = 1113895, 
                             n = 36189083,
                             k = 26382,
                             lower.tail = F)
aa_enhancer_either_FPR <- sum(aa_snp_qtl_marker_gtex_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_pancan_match_3))] == 0)/36189083
aa_enhancer_either_FNR <- (1113895 - sum(aa_snp_qtl_marker_gtex_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_pancan_match_3))] == 1))/1113895
aa_enhancer_either_TPR <- sum(aa_snp_qtl_marker_gtex_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_pancan_match_3))])/1113895
aa_enhancer_either_TNR <- (36189083 - sum(aa_snp_qtl_marker_gtex_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_pancan_match_3))] == 0))/36189083

aaxx <- aa_snp_qtl_marker_gtex_match_3 + aa_snp_qtl_marker_pancan_match_3 

aa_pval_enh_both <- phyper(q = sum(aaxx[!duplicated(names(aaxx))] == 2),
                           m = 184712, 
                           n = 37118266,
                           k = 26382, 
                           lower.tail = F)
aa_enhancer_both_FPR
aa_enhancer_both_FNR
aa_enhancer_both_TPR
aa_enhancer_both_TNR


aa_pval_enh_gtex
aa_pval_enh_pancan
aa_pval_enh_either
aa_pval_enh_both

aa_enhancer_gtex_FPR
aa_enhancer_FNR 
aa_enhancer_TPR 
aa_enhancer_TNR 

######## what is gained using models: what is the Number of False positives, what is the number of False negatives 



aa_all_dif_snp_perc_4
aa_any_change <- abs(aa_all_dif_snp_perc_4) > 0.001
aa_any_change <- rowSums(aa_any_change)
aa_any_change <- aa_any_change/244
hist(aa_any_change, breaks = 150)

aa_thresh <- seq(0,0.99,0.05)
aa_thresh_pval_gtex_FPR <- numeric(length = length(aa_thresh))
aa_thresh_pval_pancan_FPR <- numeric(length = length(aa_thresh))
aa_thresh_pval_either_FPR <- numeric(length = length(aa_thresh))
aa_thresh_pval_both_FPR <- numeric(length = length(aa_thresh))

aa_thresh_pval_gtex_FNR <- numeric(length = length(aa_thresh))
aa_thresh_pval_pancan_FNR <- numeric(length = length(aa_thresh))
aa_thresh_pval_either_FNR <- numeric(length = length(aa_thresh))
aa_thresh_pval_both_FNR <- numeric(length = length(aa_thresh))

aa_thresh_pval_gtex_TNR <- numeric(length = length(aa_thresh))
aa_thresh_pval_pancan_TNR <- numeric(length = length(aa_thresh))
aa_thresh_pval_either_TNR <- numeric(length = length(aa_thresh))
aa_thresh_pval_both_TNR <- numeric(length = length(aa_thresh))

aa_thresh_pval_gtex_TPR <- numeric(length = length(aa_thresh))
aa_thresh_pval_pancan_TPR <- numeric(length = length(aa_thresh))
aa_thresh_pval_either_TPR <- numeric(length = length(aa_thresh))
aa_thresh_pval_both_TPR <- numeric(length = length(aa_thresh))



names(aa_thresh_pval_gtex_FPR) <- aa_thresh
names(aa_thresh_pval_pancan_FPR) <- aa_thresh
names(aa_thresh_pval_both_FPR) <- aa_thresh
names(aa_thresh_pval_either_FPR) <- aa_thresh

names(aa_thresh_pval_gtex_FNR) <- aa_thresh
names(aa_thresh_pval_pancan_FNR) <- aa_thresh
names(aa_thresh_pval_both_FNR) <- aa_thresh
names(aa_thresh_pval_either_FNR) <- aa_thresh

names(aa_thresh_pval_gtex_TPR) <- aa_thresh
names(aa_thresh_pval_pancan_TPR) <- aa_thresh
names(aa_thresh_pval_both_TPR) <- aa_thresh
names(aa_thresh_pval_either_TPR) <- aa_thresh

names(aa_thresh_pval_gtex_TNR) <- aa_thresh
names(aa_thresh_pval_pancan_TNR) <- aa_thresh
names(aa_thresh_pval_both_TNR) <- aa_thresh
names(aa_thresh_pval_either_TNR) <- aa_thresh

for(i in 1:length(aa_thresh)){
  aa_thresh_pval_gtex[i] <- phyper(q = sum(aa_snp_qtl_marker_gtex_match_3[aa_any_change > aa_thresh[i]]),
                                   m = 935185, 
                                   n = 36367793,
                                   k = sum(aa_any_change > aa_thresh[i]), lower.tail = F)
  aang <- sum(aa_snp_qtl_marker_gtex_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_match_3))] == 0)
  aaps <- sum(aa_snp_qtl_marker_gtex_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_match_3))] == 1)
  aa_thresh_pval_gtex_FPR[i] <- (sum(aa_snp_qtl_marker_gtex_match_3[aa_any_change > aa_thresh[i]] == 0))/36367793
  aa_thresh_pval_gtex_FNR[i] <- ((935185 - aaps) + sum(aa_snp_qtl_marker_gtex_match_3[aa_any_change < aa_thresh[i]] == 1))/935185
  aa_thresh_pval_gtex_TPR[i] <- (sum(aa_snp_qtl_marker_gtex_match_3[aa_any_change > aa_thresh[i]]))/935185
  aa_thresh_pval_gtex_TNR[i] <- ((36367793 - aang) + sum(aa_snp_qtl_marker_gtex_match_3[aa_any_change < aa_thresh[i]] == 0))/(36367793)
  
  aa_thresh_pval_pancan[i] <- phyper(q = sum(aa_snp_qtl_marker_pancan_match_3[aa_any_change > aa_thresh[i]]),
                                     m = 363422, 
                                     n = 36939556,
                                     k = sum(aa_any_change > aa_thresh[i]), lower.tail = F)
  aang <- sum(aa_snp_qtl_marker_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_pancan_match_3))] == 0)
  aaps <- sum(aa_snp_qtl_marker_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_pancan_match_3))] == 1)
  aa_thresh_pval_pancan_FPR[i] <- (sum(aa_snp_qtl_marker_pancan_match_3[aa_any_change > aa_thresh[i]] == 0))/36939556
  aa_thresh_pval_pancan_FNR[i] <- ((363422 - aaps) + sum(aa_snp_qtl_marker_pancan_match_3[aa_any_change < aa_thresh[i]] == 1))/363422
  aa_thresh_pval_pancan_TPR[i] <- (sum(aa_snp_qtl_marker_pancan_match_3[aa_any_change > aa_thresh[i]]))/363422
  aa_thresh_pval_pancan_TNR[i] <- ((36939556 - aang) + sum(aa_snp_qtl_marker_pancan_match_3[aa_any_change < aa_thresh[i]] == 0))/(36939556)
  
  aa_thresh_pval_either[i] <- phyper(q = sum(aa_snp_qtl_marker_gtex_pancan_match_3[aa_any_change > aa_thresh[i]]),
                                     m = 1113895, 
                                     n = 36189083,
                                     k = sum(aa_any_change > aa_thresh[i]), lower.tail = F)
  aang <- sum(aa_snp_qtl_marker_gtex_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_pancan_match_3))] == 0)
  aaps <- sum(aa_snp_qtl_marker_gtex_pancan_match_3[!duplicated(names(aa_snp_qtl_marker_gtex_pancan_match_3))] == 1)
  
  aa_thresh_pval_either_FPR[i] <- (sum(aa_snp_qtl_marker_gtex_pancan_match_3[aa_any_change > aa_thresh[i]] == 0))/36189083
  aa_thresh_pval_either_FNR[i] <- ((1113895 - aaps) +  sum(aa_snp_qtl_marker_gtex_pancan_match_3[aa_any_change < aa_thresh[i]] == 1))/1113895
  aa_thresh_pval_either_TPR[i] <- (sum(aa_snp_qtl_marker_gtex_pancan_match_3[aa_any_change > aa_thresh[i]]))/1113895
  aa_thresh_pval_either_TNR[i] <- ((36189083 - aang) + sum(aa_snp_qtl_marker_gtex_pancan_match_3[aa_any_change < aa_thresh[i]] == 0))/(36189083)
  
  aa_thresh_pval_both[i] <- phyper(q = sum((aa_snp_qtl_marker_gtex_match_3[aa_any_change > aa_thresh[i]]) & (aa_snp_qtl_marker_pancan_match_3[aa_any_change > aa_thresh[i]])),
                                   m = 184712, 
                                   n = 37118266,
                                   k = sum(aa_any_change > aa_thresh[i]), lower.tail = F)
  
}

par(mfrow = c(4,3), mar = c(2,2,2,2))
plot(aa_thresh_pval_gtex_FPR,
     ylim = c(min(aa_thresh_pval_gtex_FPR), aa_enhancer_gtex_FPR),
     main = "FPR gtex", ylab = "")
abline(h=aa_enhancer_gtex_FPR, col = 2)

plot(aa_thresh_pval_pancan_FPR,
     ylim = c(min(aa_thresh_pval_pancan_FPR), aa_enhancer_pancan_FPR),
     main = "FPR pancan", ylab = "")
abline(h=aa_enhancer_pancan_FPR, col = 2)

plot(aa_thresh_pval_either_FPR,
     ylim = c(8.906001e-05, 0.000699), 
     main = "FPR either", ylab = "")
abline(h=aa_enhancer_either_FPR, col = 2)


plot(aa_thresh_pval_gtex_FNR,
     #ylim = c(min(aa_thresh_pval_gtex_FNR), aa_enhancer_gtex_FNR),
     main = "FNR gtex", ylab = "")
abline(h=aa_enhancer_gtex_FNR, col = 2)

plot(aa_thresh_pval_pancan_FNR,
     #ylim = c(min(aa_thresh_pval_pancan_FPR), aa_enhancer_pancan_FPR),
     main = "FNR pancan", ylab = "")
abline(h=aa_enhancer_pancan_FNR, col = 2)

plot(aa_thresh_pval_either_FNR,
     #ylim = c(8.906001e-05, 0.000699), 
     main = "FNR either", ylab = "")
abline(h=aa_enhancer_either_FNR, col = 2)


plot(aa_thresh_pval_gtex_TNR,
     #ylim = c(min(aa_thresh_pval_gtex_FNR), aa_enhancer_gtex_FNR),
     main = "TNR gtex", ylab = "")
abline(h=aa_enhancer_gtex_TNR, col = 2)

plot(aa_thresh_pval_pancan_TNR,
     #ylim = c(min(aa_thresh_pval_pancan_FPR), aa_enhancer_pancan_FPR),
     main = "TNR pancan", ylab = "")
abline(h=aa_enhancer_pancan_TNR, col = 2)

plot(aa_thresh_pval_either_TNR,
     #ylim = c(8.906001e-05, 0.000699), 
     main = "TNR either", ylab = "")
abline(h=aa_enhancer_either_TNR, col = 2)

plot(aa_thresh_pval_gtex_TPR, 
     ylim = c(min(aa_thresh_pval_gtex_TPR), aa_enhancer_gtex_TPR),
     main = "TPR gtex", ylab = "")
abline(h=aa_enhancer_gtex_TPR, col = 2)

plot(aa_thresh_pval_pancan_TPR,
     ylim = c(min(aa_thresh_pval_pancan_TPR), aa_enhancer_pancan_TPR),
     main = "TPR pancan", ylab = "")
abline(h=aa_enhancer_pancan_TPR, col = 2)

plot(aa_thresh_pval_either_TPR,
     ylim = c(min(aa_thresh_pval_either_TPR), aa_enhancer_either_TPR), 
     main = "TPR either", ylab = "")
abline(h=aa_enhancer_either_TPR, col = 2)
############################################################################################################

# check TPR, FPR, FNR, TNR at the same number of predictions from each method:
#  sample 

# create two vectors: length equal to number of enhancers, for each enhancer note the name of snps/eqtls on that enhancer
# get the number of unique ones from the above two vectors

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_snp <- readPeakFile(aa_snp_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
start(aa_snp) <- start(aa_snp) - 2
end(aa_snp) <- end(aa_snp) - 1

names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

# removing duplicates
aa_snpdf <- as.data.frame(aa_snp)
aa_is_dup1 <- duplicated(aa_snpdf[,1:(ncol(aa_snpdf) - 1)])
aaremove <- which(aa_is_dup1 > 0)
aa_snp <- aa_snp[-c(aaremove)]
aa_snpdf <- as.data.frame(aa_snp)
aa_snp <- makeGRangesFromDataFrame(aa_snpdf, keep.extra.columns = T)

aa_names <-  as.character(levels(aa_snpdf$V10)[as.numeric(aa_snpdf$V10)])

aa_overlap_holder <- findOverlaps(query = aa_pos_neg, subject = aa_snp)
aa_ol_enhancers_index <- unique(aa_overlap_holder@from)

enhancer_snp_list <- list()
enhancer_eqtl_list <- list()
for(i in 1:length(aa_pos_neg)){
  if(i %in% aa_ol_enhancers_index){
    enhancer_snp_list[[i]] <- unique(aa_names[aa_overlap_holder@to[aa_overlap_holder@from %in% i]])
    enhancer_eqtl_list[[i]] <- enhancer_snp_list[[i]][aa_snp_qtl_marker_gtex_pancan[match(enhancer_snp_list[[i]], names(aa_snp_qtl_marker_gtex_pancan))] == 1]
  }else{
    enhancer_snp_list[[i]] <- character(0)
    enhancer_eqtl_list[[i]] <- character(0)
  }
  
}
names(enhancer_snp_list) <- names(aa_pos_neg)
names(enhancer_eqtl_list) <- names(aa_pos_neg)

enhancer_snp_size <- unlist(lapply(enhancer_snp_list, length))
enhancer_qtl_size <- unlist(lapply(enhancer_eqtl_list, length))

aasample_step <- 50
aasample_size <- seq(5,length(aa_pos_neg), aasample_step)

aa_iter_nu <- 10
aa_enh_TPR_either <- numeric(length(aasample_size)*aa_iter_nu)
aa_pred_num <- numeric(length(aasample_size)*aa_iter_nu)
aa_TP_num <- numeric(length(aasample_size)*aa_iter_nu)
aa_TPR_model_max <- numeric(length(aasample_size)*aa_iter_nu)
aa_TPR_model_med <- numeric(length(aasample_size)*aa_iter_nu)
aa_TPR_model_mean <- numeric(length(aasample_size)*aa_iter_nu)
aa_TPR_model_perc <- numeric(length(aasample_size)*aa_iter_nu)

aa_model_TPR_either <- matrix(nrow = length(aa_enh_TPR_either), 
                              ncol = ncol(aa_all_dif_snp_perc_4))

aa_max_row <- apply(X = abs(aa_all_dif_snp_perc_4), MARGIN = 1, FUN = max)
aa_med_row <- apply(X = abs(aa_all_dif_snp_perc_4), MARGIN = 1, FUN = median)
aa_mea_row <- apply(X = abs(aa_all_dif_snp_perc_4), MARGIN = 1, FUN = mean)
aa_rowname_sorted_max <- rownames(aa_all_dif_snp_perc_4)[sort(aa_max_row,
                                                              decreasing = T,index.return=T)$ix]
aa_rowname_sorted_mea <- rownames(aa_all_dif_snp_perc_4)[sort(aa_mea_row, 
                                                              decreasing = T,index.return=T)$ix]
aa_rowname_sorted_med <- rownames(aa_all_dif_snp_perc_4)[sort(aa_med_row, 
                                                              decreasing = T,index.return=T)$ix]

aa_rowname_sorted_max_sp <- unique(unlist(lapply(strsplit(aa_rowname_sorted_max, "_"), "[[", 2)))
aa_rowname_sorted_mea_sp <- unique(unlist(lapply(strsplit(aa_rowname_sorted_mea, "_"), "[[", 2)))
aa_rowname_sorted_med_sp <- unique(unlist(lapply(strsplit(aa_rowname_sorted_med, "_"), "[[", 2)))

aa_percModels <- rowSums(abs(aa_all_dif_snp_perc_4) > 0.001)/ncol(aa_all_dif_snp_perc_4)
aa_perc_names <- rownames(aa_all_dif_snp_perc_4)[sort(aa_percModels,
                                                      decreasing = T,index.return=T)$ix]

aa_perc_names_sp <- unique(unlist(lapply(strsplit(aa_perc_names, "_"), "[[", 2)))

aa_all_dif_snp_perc_4_sind <- matrix(nrow = nrow(aa_all_dif_snp_perc_4), ncol = ncol(aa_all_dif_snp_perc_4))
for(i in 1:ncol(aa_all_dif_snp_perc_4)){
  aa_all_dif_snp_perc_4_sind[, i] <- rownames(aa_all_dif_snp_perc_4)[sort(abs(aa_all_dif_snp_perc_4[, j]),decreasing = T,index.return=T)$ix]
  aa_all_dif_snp_perc_4_sind[, i] <- unlist(lapply(strsplit(aa_all_dif_snp_perc_4_sind[, i], "_"), "[[", 2))
       
}

for(i in 1:length(aa_enh_TPR_either)){
  print(i)
  aa_cur_size <- aasample_size[floor(i/aa_iter_nu) + 1 * (as.numeric(i != length(aa_enh_TPR_either)))]
  aa_Samp <- sample(x = c(1:length(aa_pos_neg)), size = aa_cur_size, replace = F)
  aa_TP_num[i] <- sum(enhancer_qtl_size[aa_Samp])
  aa_pred_num[i] <- sum(enhancer_snp_size[aa_Samp])
  aa_enh_TPR_either[i] <- aa_TP_num[i]/aa_pred_num[i]
  
  # computing for median, max, and mean of all models


  aa_rowname_sorted_max_c <- aa_rowname_sorted_max_sp[c(1:aa_pred_num[i])]
  aa_rowname_sorted_mea_c <- aa_rowname_sorted_mea_sp[c(1:aa_pred_num[i])]
  aa_rowname_sorted_med_c <- aa_rowname_sorted_med_sp[c(1:aa_pred_num[i])]
  
  aaTPModel_max <- sum(aa_snp_qtl_marker_gtex_pancan[match(aa_rowname_sorted_max_c, 
                                                           names(aa_snp_qtl_marker_gtex_pancan))], na.rm = T)
  aaTPModel_med <- sum(aa_snp_qtl_marker_gtex_pancan[match(aa_rowname_sorted_med_c,
                                                           names(aa_snp_qtl_marker_gtex_pancan))], na.rm = T)
  aaTPModel_mea <- sum(aa_snp_qtl_marker_gtex_pancan[match(aa_rowname_sorted_mea_c,
                                                           names(aa_snp_qtl_marker_gtex_pancan))], na.rm = T)
  aa_TPR_model_max[i] <- aaTPModel_max/aa_pred_num[i]
  aa_TPR_model_med[i] <- aaTPModel_med/aa_pred_num[i]
  aa_TPR_model_mean[i] <- aaTPModel_mea/aa_pred_num[i]
  # computing for percent models predicting some change greater than 0.001

  aa_perc_names_c <- aa_perc_names_sp[c(1:aa_pred_num[i])]
  aaTPModel_perc <- sum(aa_snp_qtl_marker_gtex_pancan[match(aa_perc_names_c,
                                                            names(aa_snp_qtl_marker_gtex_pancan))], na.rm = T)
  aa_TPR_model_perc[i] <- aaTPModel_perc/aa_pred_num[i]
  
 # computing rates per model
  for(j in 1:ncol(aa_model_TPR_either)){
    print(j)
    aaModelPred <- unique(aa_all_dif_snp_perc_4_sind[,j])[c(1:aa_pred_num[i])]
    
    aaTPModel <- sum(aa_snp_qtl_marker_gtex_pancan[match(aaModelPred, 
                                                         names(aa_snp_qtl_marker_gtex_pancan))], na.rm=T)
    aa_model_TPR_either[i, j] <- aaTPModel/aa_pred_num[i]
  }

  
}

par(mfrow = c(2,2), mar = c(4,4,4,4))
plot(aa_enh_TPR_either, aa_TPR_model_max, main="max", ylab = "Model", xlab = "enhancer")
abline(a = 0, b = 1, col=2, lty = 2)
plot(aa_enh_TPR_either, aa_TPR_model_med, main="median", ylab = "Model", xlab = "enhancer")
abline(a = 0, b = 1, col=2, lty = 2)
plot(aa_enh_TPR_either, aa_TPR_model_mean, main="mean", ylab = "Model", xlab = "enhancer")
abline(a = 0, b = 1, col=2, lty = 2)
plot(aa_enh_TPR_either, aa_TPR_model_perc, main="percentile", ylab = "Model", xlab = "enhancer")
abline(a = 0, b = 1, col=2, lty = 2)


par(mfrow = c(2,2), mar = c(4,4,4,4))
plot(aa_pred_num, aa_enh_TPR_either, main="max", ylab = "TPR", xlab = "#predictions" , type = "p"
     , xlim = c(0,7500)
     )
points(aa_pred_num, aa_TPR_model_max, col = 2)

plot(aa_pred_num, aa_enh_TPR_either, main="median", 
     ylab = "TPR", xlab = "#predictions", type = "p"
     , xlim = c(0,7500)
     )
points(aa_pred_num, aa_TPR_model_med, col = 2)

plot(aa_pred_num, aa_enh_TPR_either, main="mean",
     ylab = "TPR", xlab = "#predictions", type = "p"
     , xlim = c(0,7500)
     )
points(aa_pred_num, aa_TPR_model_mean, col = 2)

plot(aa_pred_num, aa_enh_TPR_either, main="perc", 
     ylab = "TPR", xlab = "#predictions", type = "p"
     , xlim = c(0,7500)
     )
points(aa_pred_num, aa_TPR_model_perc, col = 2)



par(mfrow = c(10, 10), mar =c(1, 1, 1, 1))
for(i in 1:ncol(aa_model_TPR_either)){
  plot(aa_enh_TPR_either, aa_model_TPR_either[, i], 
       main=colnames(aa_all_dif_snp_perc_4)[i],
       xaxt = "n", yaxt = "n")
  abline(a = 0, b = 1, col=2, lty = 2)
  
}

par(mfrow = c(10, 10), mar =c(1, 1, 1, 1))
for(i in 1:ncol(aa_model_TPR_either)){
  plot(aa_pred_num,aa_enh_TPR_either, 
       main=colnames(aa_all_dif_snp_perc_4)[i],
       xaxt = "n", yaxt = "n", type="l")
  lines(aa_pred_num, aa_model_TPR_either[, i], col = 2)
  abline(a = 0, b = 1, col=2, lty = 2)
  
}
aa_snp_qtl_marker_gtex
aa_snp_qtl_marker_pancan
aa_snp_qtl_marker_gtex_pancan

##################################################################################################
##################################################################################################
##################################################################################################
# create a set of sequences containing variants to test:
# somatic mutations from cosmic <- get interesect with my enhancers


#somatic from cosmic
ER_somatic_mutation_all <- read.csv("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/V90_38_NONCODINGVARIANTS.csv",
                                    sep = ",", header = T, stringsAsFactors = F)

aa <- ER_somatic_mutation_all$GENOME_POSITION
aa_chr <- unlist(lapply(strsplit(aa, ":"), "[[", 1))
aa_chr[aa_chr == "23"] <- "X"
aa_chr[aa_chr == "24"] <- "Y"
aa_chr <- paste0("chr", aa_chr)
aa_nch <- unlist(lapply(strsplit(aa, ":"), "[[", 2))
aa_nch_sp <- strsplit(aa_nch, "-")
aa_sta <- unlist(lapply(aa_nch_sp, "[[", 1))
aa_end <- unlist(lapply(aa_nch_sp, "[[", 2))

aa_cosmic_df <- data.frame(chr = aa_chr,
                         start = aa_sta, 
                         end = aa_end,
                         stringsAsFactors = F)

colnames(ER_somatic_mutation_all)
unique(ER_somatic_mutation_all$FATHMM_MKL_NON_CODING_GROUPS)
aa_cosmic_df <- cbind(aa_cosmic_df, ER_somatic_mutation_all[, c(1,2,3,4,5,8,9,10,11,12,13,14,15,17,18,19,20,21,23,25,26,27,28)])

aa_cosmic_gr <- makeGRangesFromDataFrame(aa_cosmic_df, keep.extra.columns = T)
options(scipen=999)
write.table(aa_cosmic_gr, 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/cosmic_breast_somatic.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
options(scipen=0)

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_somatic_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/cosmic_breast_somatic.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_somatic <- readPeakFile(aa_somatic_Add, as = "GRanges")

aa_ov_pos <- findOverlaps(query = aa_pos, subject = aa_somatic)
length(unique(aa_ov_pos@to))
aa_ov_neg <- findOverlaps(query = aa_neg, subject = aa_somatic)
length(unique(aa_ov_neg@to))
hist(as.numeric(table(aa_ov_pos@from)), breaks = 100)


aa_somatic_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/cosmic_breast_somatic.bed"
aaa_somatic_pos4 <- bedtools_intersect(bedfile_names = c(aa_somatic_Add, aa_pos_Add), 
                                   wa=F, wb=F, loj=F, wo=F, wao=F,
                                   u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                   S=F, sorted=F, merge_after_distance=0, 
                                   output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/somatic_pos_4.bed", 
                                   read_output = T)
aaa_somatic_neg4 <- bedtools_intersect(bedfile_names = c(aa_somatic_Add, aa_neg_Add), 
                                   wa=F, wb=F, loj=F, wo=F, wao=F,
                                   u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                   S=F, sorted=F, merge_after_distance=0, 
                                   output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/somatic_neg_4.bed", 
                                   read_output = T)
start(aaa_somatic_pos4) <- start(aaa_somatic_pos4) - 1
start(aaa_somatic_neg4) <- start(aaa_somatic_neg4) - 1

aa_somatic_pos_neg_4 <- c(aaa_somatic_pos4, aaa_somatic_neg4)


options(scipen=999)

write.table(aa_somatic_pos_neg_4, 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/somatic_pos_neg4.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
options(scipen=0)

# write somatic mutations

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_somatic_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/somatic_pos_neg4.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_somatic <- readPeakFile(aa_somatic_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
start(aa_somatic) <- start(aa_somatic) - 1
#end(aa_snp) <- end(aa_snp) - 1

library(BSgenome.Hsapiens.UCSC.hg38)
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC")
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

# removing duplicates
aa_somatic_df <- as.data.frame(aa_somatic)
aa_is_dup1 <- duplicated(aa_somatic_df)
aaremove <- which(aa_is_dup1 > 0)
aa_somatic <- aa_somatic[-c(aaremove)]
aa_somatic_df <- as.data.frame(aa_somatic)

aa_ref <- as.character(levels(aa_somatic_df$V22)[as.numeric(aa_somatic_df$V22)])
aa_alt <- as.character(levels(aa_somatic_df$V23)[as.numeric(aa_somatic_df$V23)])

aa_ref_Alt <- cbind(aa_ref, aa_alt)
unique(aa_ref)
unique(aa_alt)
aa_names <-  paste0("somatic_", c(1:length(aa_somatic)))
# fix nulls in reg and alt: in alt just set them to "", in ref find the ref sequence, replace ref eith ref sequence and update alt to include the ref seq
aa_ref_Alt2 <- aa_ref_Alt
for(i in 1:nrow(aa_ref_Alt)){
  if(aa_ref_Alt[i, 1] == "null"){
    aa_ref_Alt[i, 1] <- as.character(getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = aa_somatic[i]))
    aa_ref_Alt[i, 2] <- paste0(aa_ref_Alt[i, 1], aa_ref_Alt[i, 2])
  }else if(aa_ref_Alt[i, 2] == "null"){
    aa_ref_Alt[i, 2] <- ""
  }
}
aaaa <- unlist(lapply(strsplit(aa_ref_Alt[, 2], split = ","), length))
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC")
aa_pos_neg_vars_snp <- write_variant_seq(enhancer_GR = aa_pos_neg, 
                                         eqtl_GR = aa_somatic,
                                         eqtl_names = aa_names,
                                         eqtl_ref_alt = aa_ref_Alt,
                                         all_combs = F,
                                         my_genome = BSgenome.Hsapiens.UCSC.hg38, 
                                         label = c(rep(1, length(aa_pos)), rep(0, length(aa_neg))))

# write GEMSTAT jobs to run
# write hal job for KDs
aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas1[1:150]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas2[1:150]]
aaadd3 <- union(aaadd1, aaadd2)

aa_par_name <- paste0(aaadd3, ".txt")
aa_seq_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/Variant_seq/",
                           recursive = T)
aa_lab_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/Variant_Labels/",
                           recursive = T)
aa_nam_spl <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_seq_files, 
                                                                   split = "\\."), 
                                                          "[[", 1)), split = "_"),
                                   "[", c(3,4)),
                            paste, collapse = "_"))


setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC")
for(i in 1:length(aa_par_name)){
  for(j in 1:length(aa_seq_files)){
    cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr ",
          paste0("-s Variant_seq/", aa_seq_files[j]), 
          paste0("-e Variant_Labels/", aa_lab_files[j]),
          "-m motifs.wtmx -f TF_exp.tab", 
          paste0("-fo Varaint_out/", aa_nam_spl[j], "_", aaadd3[i], ".out"), 
          "-o DIRECT -c Coop/coop.par ", 
          paste0("-p Trained_par/", aa_par_name[i]),
          "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1\n"
    ), 
    sep = " ", append = !(i==1 & j==1), file = "variant_jobs_cosmic_exp18.job")
  }
}

################# ################# ################# ################# #################
# read cosmic results
load("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/hal_cosmic.RData")
pos_cosmic_output_hal
neg_cosmic_output_hal
WT_models_cosmic_hal

aa_neg_dif_cosmic_3 <- list()
aa_pos_dif_cosmic_3 <- list()

# aa_neg_dif_cosmic_percentile_3 <- list()
# aa_pos_dif_cosmic_percentile_3 <- list()

aa_neg_dif_cosmic_percentile_4 <- list()
aa_pos_dif_cosmic_percentile_4 <- list()

# aa_Exp18_WT_model_prediction_Filtered_GT <- rownames(Exp18_WT_model_prediction_Filtered)
# aa_Exp18_WT_model_prediction_Filtered_GT <- unlist(lapply(strsplit(aa_Exp18_WT_model_prediction_Filtered_GT, split = "_"), "[[", 1))
# aa_Exp18_WT_model_prediction_Filtered_GT[aa_Exp18_WT_model_prediction_Filtered_GT == "pos"] <- 1
# aa_Exp18_WT_model_prediction_Filtered_GT[aa_Exp18_WT_model_prediction_Filtered_GT == "neg"] <- 0
# aa_Exp18_WT_model_prediction_Filtered_GT <- as.numeric(aa_Exp18_WT_model_prediction_Filtered_GT)
# names(aa_Exp16_WT_model_prediction_Filtered_GT) <- rownames(Exp16_WT_model_prediction_Filtered)

aa_WT_GT <- rownames(WT_models_cosmic_hal)
aa_WT_GT <- unlist(lapply(strsplit(aa_WT_GT, "_"), "[[", 1))
aa_WT_GT[aa_WT_GT  == "pos"] <- 1
aa_WT_GT[aa_WT_GT  == "neg"] <- 0
aa_WT_GT <- as.numeric(aa_WT_GT)
names(aa_WT_GT) <- rownames(WT_models_cosmic_hal)

for(i in 1:length(neg_cosmic_output_hal)){
  print(i)
  aa_neg_dif_cosmic_3[[i]] <- matrix(nrow = nrow(neg_cosmic_output_hal[[i]]), ncol = ncol(neg_cosmic_output_hal[[i]]))
  colnames(aa_neg_dif_cosmic_3[[i]]) <- colnames(neg_cosmic_output_hal[[i]])
  rownames(aa_neg_dif_cosmic_3[[i]]) <- rownames(neg_cosmic_output_hal[[i]])
  
  # aa_neg_dif_snp_percentile_3[[i]] <- matrix(nrow = nrow(neg_snp_output_hal[[i]]), ncol = ncol(neg_snp_output_hal[[i]]))
  # colnames(aa_neg_dif_snp_percentile_3[[i]]) <- colnames(neg_snp_output_hal[[i]])
  # rownames(aa_neg_dif_snp_percentile_3[[i]]) <- rownames(neg_snp_output_hal[[i]])
  
  aa_neg_dif_cosmic_percentile_4[[i]] <- matrix(nrow = nrow(neg_cosmic_output_hal[[i]]), ncol = ncol(neg_cosmic_output_hal[[i]]))
  colnames(aa_neg_dif_cosmic_percentile_4[[i]]) <- colnames(neg_cosmic_output_hal[[i]])
  rownames(aa_neg_dif_cosmic_percentile_4[[i]]) <- rownames(neg_cosmic_output_hal[[i]])
  
  for(j in 1:nrow(neg_cosmic_output_hal[[i]])){
    aa_neg_dif_cosmic_3[[i]][j, ] <- (neg_cosmic_output_hal[[i]][j, ] - WT_models_cosmic_hal[names(neg_cosmic_output_hal)[i], ])/ WT_models_cosmic_hal[names(neg_cosmic_output_hal)[i], ]
    for(kk in 1:ncol(neg_cosmic_output_hal[[i]])){
      #if(!is.na(Exp18_WT_model_prediction_Filtered[names(neg_snp_output_hal)[i],kk])){
      # aa_pos_perc_before <- sum(WT_models_snp_hal[aa_WT_GT == 1,kk] <= WT_models_snp_hal[names(neg_cosmic_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 1,kk]))
      # aa_neg_perc_before <- sum(WT_models_snp_hal[aa_WT_GT == 0,kk] >  WT_models_snp_hal[names(neg_cosmic_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 0,kk]))
      # aa_pos_perc_after <-  sum(WT_models_snp_hal[aa_WT_GT == 1,kk] <=  neg_cosmic_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 1,kk]))
      # aa_neg_perc_after <-  sum(WT_models_snp_hal[aa_WT_GT == 0,kk] >   neg_cosmic_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 0,kk]))
      # aa_neg_dif_snp_percentile_3[[i]][j,kk] <- -(aa_neg_perc_after - aa_neg_perc_before) + (aa_pos_perc_after - aa_pos_perc_before)
      aa_all_before <- sum(WT_models_cosmic_hal[,kk] <= WT_models_cosmic_hal[names(neg_cosmic_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[,kk]))
      aa_all_after <- sum(WT_models_cosmic_hal[,kk] <= neg_cosmic_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[,kk]))
      aa_neg_dif_cosmic_percentile_4[[i]][j,kk] <- aa_all_after - aa_all_before
      #}
    }
  }
}
names(aa_neg_dif_cosmic_3) <- names(neg_cosmic_output_hal)
#names(aa_neg_dif_snp_percentile_3) <- names(neg_cosmic_output_hal)
names(aa_neg_dif_cosmic_percentile_4) <- names(neg_cosmic_output_hal)

neg_cosmic_output_hal_all <- do.call(what = rbind, neg_cosmic_output_hal)
aa_neg_dif_all_cosmic_3 <- do.call(what = rbind, aa_neg_dif_cosmic_3)
#aa_neg_dif_snp_all_percentile_3 <- do.call(what = rbind, aa_neg_dif_snp_percentile_3)
aa_neg_dif_cosmic_all_percentile_4 <- do.call(what = rbind, aa_neg_dif_cosmic_percentile_4)



for(i in 1:length(pos_cosmic_output_hal)){
  print(i)
  aa_pos_dif_cosmic_3[[i]] <- matrix(nrow = nrow(pos_cosmic_output_hal[[i]]), ncol = ncol(pos_cosmic_output_hal[[i]]))
  colnames(aa_pos_dif_cosmic_3[[i]]) <- colnames(pos_cosmic_output_hal[[i]])
  rownames(aa_pos_dif_cosmic_3[[i]]) <- rownames(pos_cosmic_output_hal[[i]])
  
  # aa_pos_dif_snp_percentile_3[[i]] <- matrix(nrow = nrow(pos_cosmic_output_hal[[i]]), ncol = ncol(pos_cosmic_output_hal[[i]]))
  # colnames(aa_pos_dif_snp_percentile_3[[i]]) <- colnames(pos_cosmic_output_hal[[i]])
  # rownames(aa_pos_dif_snp_percentile_3[[i]]) <- rownames(pos_cosmic_output_hal[[i]])
  
  aa_pos_dif_cosmic_percentile_4[[i]] <- matrix(nrow = nrow(pos_cosmic_output_hal[[i]]), ncol = ncol(pos_cosmic_output_hal[[i]]))
  colnames(aa_pos_dif_cosmic_percentile_4[[i]]) <- colnames(pos_cosmic_output_hal[[i]])
  rownames(aa_pos_dif_cosmic_percentile_4[[i]]) <- rownames(pos_cosmic_output_hal[[i]])
  
  for(j in 1:nrow(pos_cosmic_output_hal[[i]])){
    aa_pos_dif_cosmic_3[[i]][j, ] <- (pos_cosmic_output_hal[[i]][j, ] - WT_models_cosmic_hal[names(pos_cosmic_output_hal)[i], ])/ WT_models_cosmic_hal[names(pos_cosmic_output_hal)[i], ]
    for(kk in 1:ncol(pos_cosmic_output_hal[[i]])){
      #if(!is.na(Exp18_WT_model_prediction_Filtered[names(pos_snp_output_hal)[i],kk])){
      # aa_pos_perc_before <- sum(WT_models_cosmic_hal[aa_WT_GT == 1,kk] <= WT_models_cosmic_hal[names(pos_snp_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[aa_WT_GT == 1,kk]))
      # aa_neg_perc_before <- sum(WT_models_cosmic_hal[aa_WT_GT == 0,kk] >  WT_models_cosmic_hal[names(pos_snp_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[aa_WT_GT == 0,kk]))
      # aa_pos_perc_after <-  sum(WT_models_cosmic_hal[aa_WT_GT == 1,kk] <=  pos_snp_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[aa_WT_GT == 1,kk]))
      # aa_neg_perc_after <-  sum(WT_models_cosmic_hal[aa_WT_GT == 0,kk] >   pos_snp_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[aa_WT_GT == 0,kk]))
      # aa_pos_dif_snp_percentile_3[[i]][j,kk] <- -(aa_neg_perc_after - aa_neg_perc_before) + (aa_pos_perc_after - aa_pos_perc_before)
      aa_all_before <- sum(WT_models_cosmic_hal[,kk] <= WT_models_cosmic_hal[names(pos_cosmic_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[,kk]))
      aa_all_after <- sum(WT_models_cosmic_hal[,kk] <= pos_cosmic_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[,kk]))
      aa_pos_dif_cosmic_percentile_4[[i]][j,kk] <- aa_all_after - aa_all_before
      #}
    }
  }
}
names(aa_pos_dif_cosmic_3) <- names(pos_cosmic_output_hal)
#names(aa_pos_dif_snp_percentile_3) <- names(pos_cosmic_output_hal)
names(aa_pos_dif_cosmic_percentile_4) <- names(pos_cosmic_output_hal)

pos_cosmic_output_hal_all <- do.call(what = rbind, pos_cosmic_output_hal)
aa_pos_dif_all_cosmic_3 <- do.call(what = rbind, aa_pos_dif_cosmic_3)
#aa_pos_dif_snp_all_percentile_3 <- do.call(what = rbind, aa_pos_dif_snp_percentile_3)
aa_pos_dif_cosmic_all_percentile_4 <- do.call(what = rbind, aa_pos_dif_cosmic_percentile_4)


aa_pos_neg_out_cosmic_all_3 <- rbind(pos_cosmic_output_hal_all, neg_cosmic_output_hal_all)
aa_all_dif_cosmic_3 <- rbind(aa_pos_dif_all_cosmic_3, aa_neg_dif_all_cosmic_3)
#aa_all_dif_cosmic_perc_3 <- rbind(aa_pos_dif_cosmic_all_percentile_3, aa_neg_dif_cosmic_all_percentile_3)
aa_all_dif_cosmic_perc_4 <- rbind(aa_pos_dif_cosmic_all_percentile_4, aa_neg_dif_cosmic_all_percentile_4)


aa_all_dif_cosmic_names_3 <- unlist(lapply(strsplit(rownames(aa_all_dif_cosmic_3), split="_"), "[[", 3))

aa_snp_gwas_marker_match_3 <-  aa_snp_gwas_marker[match(aa_all_dif_snp_names_3, names(aa_snp_gwas_marker))]
aa_snp_qtl_marker_pancan_match_3 <- aa_snp_qtl_marker_pancan[match(aa_all_dif_snp_names_3, names(aa_snp_qtl_marker_pancan))]
aa_snp_qtl_marker_gtex_match_3 <- aa_snp_qtl_marker_gtex[match(aa_all_dif_snp_names_3, names(aa_snp_qtl_marker_gtex))]
aa_snp_qtl_marker_gtex_pancan_match_3 <- aa_snp_qtl_marker_gtex_pancan[match(aa_all_dif_cosmic_names_3, names(aa_snp_qtl_marker_gtex_pancan))]
#################
# find out whether high ranking snps of my model have higher scores:
# look at spearman correlation between the score for each model and the FATHMM score

nrow(aa_all_dif_cosmic_3)
nrow(aa_all_dif_cosmic_perc_4)
aa_somat_index <- as.numeric(unlist(lapply(strsplit(rownames(aa_all_dif_cosmic_3), "_"), "[[", 3)))
aa_somatic_scores_noncoding <- as.numeric(levels(aa_somatic_df$V25)[as.numeric(aa_somatic_df$V25)])
aa_somatic_scores_coding <- as.numeric(levels(aa_somatic_df$V26)[as.numeric(aa_somatic_df$V26)])

aa_speman_nc <- numeric(length = ncol(aa_all_dif_cosmic_perc_4))
aa_speman_c <- numeric(length = ncol(aa_all_dif_cosmic_perc_4))

for(i in 1:ncol(aa_all_dif_cosmic_perc_4)){
  print(i)
  aa_speman_nc[i] <- cor.test(aa_somatic_scores_noncoding[aa_somat_index],
                              abs(aa_all_dif_cosmic_perc_4)[,i],
                              method = "spearman")$estimate
  aa_speman_c[i] <- cor.test(aa_somatic_scores_coding[aa_somat_index],
                              abs(aa_all_dif_cosmic_perc_4)[,i],
                              method = "spearman")$estimate
}
hist(aa_speman_nc)
hist(aa_speman_c)
for(i in 1:ncol(aa_all_dif_cosmic_3)){
  print(i)
  aa_speman_nc[i] <- cor.test(aa_somatic_scores_noncoding[aa_somat_index],
                              abs(aa_all_dif_cosmic_3)[,i],
                              method = "spearman")$estimate
  aa_speman_c[i] <- cor.test(aa_somatic_scores_coding[aa_somat_index],
                             abs(aa_all_dif_cosmic_3)[,i],
                             method = "spearman")$estimate
}
hist(aa_speman_nc)
hist(aa_speman_c)

aa_model_index_Sorted_percentile <- matrix(nrow = nrow(aa_all_dif_cosmic_perc_4),
                                           ncol = ncol(aa_all_dif_cosmic_perc_4))
aa_model_index_Sorted_change <- matrix(nrow = nrow(aa_all_dif_cosmic_3),
                                           ncol = ncol(aa_all_dif_cosmic_3))
for(i in 1:ncol(aa_model_index_Sorted_percentile)){
  aa_model_index_Sorted_percentile[, i] <- aa_somat_index[sort(abs(aa_all_dif_cosmic_perc_4[, i]),
                                                decreasing = T, index.return = T)$ix]
  aa_model_index_Sorted_change[, i] <- aa_somat_index[sort(abs(aa_all_dif_cosmic_3[, i]),
                                            decreasing = T, index.return = T)$ix]
}
aa_max <- apply(abs(aa_all_dif_cosmic_perc_4), MARGIN = 1, FUN = max)
aa_mea <- apply(abs(aa_all_dif_cosmic_perc_4), MARGIN = 1, FUN = mean)
aa_med <- apply(abs(aa_all_dif_cosmic_perc_4), MARGIN = 1, FUN = median)

aa_model_sort_ind_max <- aa_somat_index[sort(aa_max, decreasing = T, index.return = T)$ix]
aa_model_sort_ind_mea <- aa_somat_index[sort(aa_mea, decreasing = T, index.return = T)$ix]
aa_model_sort_ind_med <- aa_somat_index[sort(aa_med, decreasing = T, index.return = T)$ix]

aa_percModels <- rowSums(abs(aa_all_dif_cosmic_perc_4) > 0.001)/ncol(aa_all_dif_cosmic_perc_4)
aa_perc_names <- aa_somat_index[sort(aa_percModels,decreasing = T,index.return=T)$ix]

aasample_size <- c(seq(1,100, 5), seq(110, length(aa_somat_index), 50))
aa_rep_num <- 244
aa_random_Score <- matrix(nrow = aa_rep_num, ncol = length(aasample_size))
colnames(aa_random_Score) <- aasample_size
aath <- 0.7

aa_model_Score_mean <- numeric(length = length(aasample_size))
aa_model_Score_median <- numeric(length = length(aasample_size))
aa_model_Score_max <- numeric(length = length(aasample_size))
aa_model_Score_perc <- numeric(length = length(aasample_size))

aa_model_Score_percentile <- matrix(nrow = length(aasample_size), ncol = ncol(aa_all_dif_cosmic_3))
colnames(aa_model_Score_percentile) <- colnames(aa_all_dif_cosmic_3)
rownames(aa_model_Score_percentile) <- aasample_size

aa_model_Score_change <- matrix(nrow = length(aasample_size), ncol = ncol(aa_all_dif_cosmic_3))
colnames(aa_model_Score_change) <- colnames(aa_all_dif_cosmic_3)
rownames(aa_model_Score_change) <- aasample_size

aa_tot <- sum(aa_somatic_scores_noncoding[aa_somat_index] >= aath, na.rm=T)
for(i in 1:length(aasample_size)){
  for(j in 1:aa_rep_num){
    aasampl <- sample(x = aa_somatic_scores_noncoding[aa_somat_index], size = aasample_size[i], replace = F)
    aa_random_Score[j, i] <- sum(aasampl >= aath, na.rm=T)/aa_tot
  }
  for(j in 1:ncol(aa_all_dif_cosmic_perc_4)){
    aa_model_Score_percentile[i, j] <- sum(aa_somatic_scores_noncoding[aa_model_index_Sorted_percentile[1:aasample_size[i],j]] >= aath, na.rm=T)/aa_tot
    aa_model_Score_change[i, j] <- sum(aa_somatic_scores_noncoding[aa_model_index_Sorted_change[1:aasample_size[i],j]] >= aath, na.rm=T)/aa_tot
    
  }
  aa_model_Score_mean[i] <- sum(aa_somatic_scores_noncoding[aa_model_sort_ind_mea[1:aasample_size[i]]] >= aath, na.rm=T)/aa_tot
  aa_model_Score_median[i] <- sum(aa_somatic_scores_noncoding[aa_model_sort_ind_med[1:aasample_size[i]]] >= aath, na.rm=T)/aa_tot
  aa_model_Score_max[i] <- sum(aa_somatic_scores_noncoding[aa_model_sort_ind_max[1:aasample_size[i]]] >= aath, na.rm=T)/aa_tot
  aa_model_Score_perc[i] <-sum(aa_somatic_scores_noncoding[aa_perc_names[1:aasample_size[i]]] >= aath, na.rm=T)/aa_tot
}

[,1:22]
boxplot.matrix(aa_random_Score[,1:22], las = 2, xlab = "#predictions",
               ylab = "TPR", main = "non_normal_snps")
points(aa_model_Score_mean[1:22], col = 2, pch = 16, cex = 0.7)
points(aa_model_Score_median[1:22], col = 3, pch = 17, cex = 0.7)
points(aa_model_Score_max[1:22], col = 4, pch = 18, cex = 0.7)
points(aa_model_Score_perc[1:22], col = 5, pch = 19, cex = 0.7)
legend(x = "topleft",legend = c("mean", "median", "max", "vote"), 
       pch = c(16,17,18,19), fill = c(2,3,4,5), cex = 0.5)

boxplot.matrix(t(aa_model_Score_percentile))

aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_random_Score) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, 104, 2)){
  aa_allmat[, i] <- aa_random_Score[, aac]
  aa_allmat[, i+1] <- t(aa_model_Score_percentile)[, aac]
  aac <- aac + 1
}
boxplot.matrix(aa_allmat, col = rep(c(2,3), 52),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "FATHHM_score_clsf_0_70")
axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
     labels = aasample_size, las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)

aa_allmat2 <- matrix(nrow = 244, ncol = ncol(aa_random_Score) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, 104, 2)){
  aa_allmat2[, i] <- aa_random_Score[, aac]
  aa_allmat2[, i+1] <- t(aa_model_Score_change)[, aac]
  aac <- aac + 1
}
boxplot.matrix(aa_allmat2[,1:44], col = rep(c(2,3), 52),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "FATHHM_score_clsf_0_70")
axis(side = 1, at = seq(1.5,ncol(aa_allmat2[,1:44]), 2), 
     labels = aasample_size[1:22], las = 2)
abline(v = seq(0.5, ncol(aa_allmat2[,1:44])+1, 2), col = 4)

### check for normal SNPs vs somatic mutations: does method find normal mutations to be more effective than somatic ones?
aasample_size <- c(seq(1,100, 5), seq(110, length(aa_somat_index), 50))
aa_rep_num <- 244
aa_random_Score <- matrix(nrow = aa_rep_num, ncol = length(aasample_size))
colnames(aa_random_Score) <- aasample_size
aath <- 0.7

aa_model_Score_mean <- numeric(length = length(aasample_size))
aa_model_Score_median <- numeric(length = length(aasample_size))
aa_model_Score_max <- numeric(length = length(aasample_size))
aa_model_Score_perc <- numeric(length = length(aasample_size))

aa_model_Score_percentile <- matrix(nrow = length(aasample_size), ncol = ncol(aa_all_dif_cosmic_3))
colnames(aa_model_Score_percentile) <- colnames(aa_all_dif_cosmic_3)
rownames(aa_model_Score_percentile) <- aasample_size

aa_model_Score_change <- matrix(nrow = length(aasample_size), ncol = ncol(aa_all_dif_cosmic_3))
colnames(aa_model_Score_change) <- colnames(aa_all_dif_cosmic_3)
rownames(aa_model_Score_change) <- aasample_size
aa_somatic_snp_vs_no <- levels(aa_somatic_df$V24)[as.numeric(aa_somatic_df$V24)]

aa_tot <- sum(aa_somatic_snp_vs_no[aa_somat_index] == "n", na.rm=T)
for(i in 1:length(aasample_size)){
  for(j in 1:aa_rep_num){
    aasampl <- sample(x = aa_somatic_snp_vs_no[aa_somat_index], size = aasample_size[i], replace = F)
    aa_random_Score[j, i] <- sum(aasampl == "n", na.rm=T)/aa_tot
  }
  for(j in 1:ncol(aa_all_dif_cosmic_perc_4)){
    aa_model_Score_percentile[i, j] <- sum(aa_somatic_snp_vs_no[aa_model_index_Sorted_percentile[1:aasample_size[i],j]] == "n", na.rm=T)/aa_tot
    aa_model_Score_change[i, j] <- sum(aa_somatic_snp_vs_no[aa_model_index_Sorted_change[1:aasample_size[i],j]] == "n", na.rm=T)/aa_tot
    
  }
  aa_model_Score_mean[i] <- sum(aa_somatic_snp_vs_no[aa_model_sort_ind_mea[1:aasample_size[i]]] == "n", na.rm=T)/aa_tot
  aa_model_Score_median[i] <- sum(aa_somatic_snp_vs_no[aa_model_sort_ind_med[1:aasample_size[i]]] == "n", na.rm=T)/aa_tot
  aa_model_Score_max[i] <- sum(aa_somatic_snp_vs_no[aa_model_sort_ind_max[1:aasample_size[i]]] == "n", na.rm=T)/aa_tot
  aa_model_Score_perc[i] <-sum(aa_somatic_snp_vs_no[aa_perc_names[1:aasample_size[i]]] == "n", na.rm=T)/aa_tot
}

### repeat the above analysis but with all common snps included
nrow(aa_all_dif_cosmic_3)
nrow(aa_all_dif_cosmic_perc_4)
aa_all_dif_cosmic_and_snp_perc_4
aa_somat_index <- as.numeric(unlist(lapply(strsplit(rownames(aa_all_dif_cosmic_3), "_"), "[[", 3)))
aa_somatic_snp_vs_no <- levels(aa_somatic_df$V24)[as.numeric(aa_somatic_df$V24)]

aa_all_dif_cosmic_and_snp_perc_4 <- rbind(aa_all_dif_cosmic_perc_4,
                                          aa_all_dif_snp_perc_4)
aa_all_dif_cosmic_and_snp_perc_4_shuffind <- sample(c(1:nrow(aa_all_dif_cosmic_and_snp_perc_4)),
                                                    size = nrow(aa_all_dif_cosmic_and_snp_perc_4), replace = F)

aa_all_dif_cosmic_and_snp_perc_4_shuffled <- aa_all_dif_cosmic_and_snp_perc_4[aa_all_dif_cosmic_and_snp_perc_4_shuffind, ]
#aa_all_dif_snp_perc_4

aa_model_index_Sorted_percentile <- matrix(nrow = nrow(aa_all_dif_cosmic_and_snp_perc_4),
                                           ncol = ncol(aa_all_dif_cosmic_and_snp_perc_4))
aa_model_index_Sorted_change <- matrix(nrow = nrow(aa_all_dif_cosmic_3),
                                       ncol = ncol(aa_all_dif_cosmic_3))
for(i in 1:ncol(aa_model_index_Sorted_percentile)){
  aa_model_index_Sorted_percentile[, i] <- sort(abs(aa_all_dif_cosmic_and_snp_perc_4_shuffled[, i]),
                                                               decreasing = T, index.return = T)$ix
  # aa_model_index_Sorted_change[, i] <- sort(abs(aa_all_dif_cosmic_3[, i]),
  #                                                          decreasing = T, index.return = T)$ix
}
aa_max <- apply(abs(aa_all_dif_cosmic_and_snp_perc_4_shuffled), MARGIN = 1, FUN = max)
aa_mea <- apply(abs(aa_all_dif_cosmic_and_snp_perc_4_shuffled), MARGIN = 1, FUN = mean)
aa_med <- apply(abs(aa_all_dif_cosmic_and_snp_perc_4_shuffled), MARGIN = 1, FUN = median)

aa_model_sort_ind_max <- sort(aa_max, decreasing = T, index.return = T)$ix
aa_model_sort_ind_mea <- sort(aa_mea, decreasing = T, index.return = T)$ix
aa_model_sort_ind_med <- sort(aa_med, decreasing = T, index.return = T)$ix

aa_percModels <- rowSums(abs(aa_all_dif_cosmic_and_snp_perc_4_shuffled) > 0.001)/ncol(aa_all_dif_cosmic_and_snp_perc_4_shuffled)
aa_perc_names <- sort(aa_percModels,decreasing = T,index.return=T)$ix


aa_non_snp <- aa_somatic_snp_vs_no[aa_somat_index]
aa_non_snp_all <- c(aa_non_snp, rep("y", nrow(aa_all_dif_snp_perc_4)))
aa_non_snp_all_shuffled <- aa_non_snp_all[aa_all_dif_cosmic_and_snp_perc_4_shuffind]


aasample_size <- c(seq(1,500, 10),
                   seq(550, 2000, 50),
                   seq(2000, 20000, 1000))
aa_rep_num <- 244
aa_random_Score <- matrix(nrow = aa_rep_num, ncol = length(aasample_size))
colnames(aa_random_Score) <- aasample_size
aath <- 0.7

aa_model_Score_mean <- numeric(length = length(aasample_size))
aa_model_Score_median <- numeric(length = length(aasample_size))
aa_model_Score_max <- numeric(length = length(aasample_size))
aa_model_Score_perc <- numeric(length = length(aasample_size))

aa_model_Score_percentile <- matrix(nrow = length(aasample_size), ncol = ncol(aa_all_dif_cosmic_3))
colnames(aa_model_Score_percentile) <- colnames(aa_all_dif_cosmic_3)
rownames(aa_model_Score_percentile) <- aasample_size

aa_model_Score_change <- matrix(nrow = length(aasample_size), ncol = ncol(aa_all_dif_cosmic_3))
colnames(aa_model_Score_change) <- colnames(aa_all_dif_cosmic_3)
rownames(aa_model_Score_change) <- aasample_size
#aa_somatic_snp_vs_no <- levels(aa_somatic_df$V24)[as.numeric(aa_somatic_df$V24)]

aa_tot <- sum(aa_non_snp_all == "n", na.rm=T)
for(i in 1:length(aasample_size)){
  print(i)
  for(j in 1:aa_rep_num){
    aasampl <- sample(x = aa_non_snp_all, size = aasample_size[i], replace = F)
    aa_random_Score[j, i] <- sum(aasampl == "n", na.rm=T)/aa_tot
  }
  for(j in 1:ncol(aa_all_dif_cosmic_and_snp_perc_4)){
    aa_model_Score_percentile[i, j] <- sum(aa_non_snp_all_shuffled[aa_model_index_Sorted_percentile[1:aasample_size[i],j]] == "n", na.rm=T)/aa_tot
    #aa_model_Score_change[i, j] <- sum(aa_non_snp_all[aa_model_index_Sorted_change[1:aasample_size[i],j]] == "n", na.rm=T)/aa_tot
    
  }
  aa_model_Score_mean[i] <- sum(aa_non_snp_all_shuffled[aa_model_sort_ind_mea[1:aasample_size[i]]] == "n", na.rm=T)/aa_tot
  aa_model_Score_median[i] <- sum(aa_non_snp_all_shuffled[aa_model_sort_ind_med[1:aasample_size[i]]] == "n", na.rm=T)/aa_tot
  aa_model_Score_max[i] <- sum(aa_non_snp_all_shuffled[aa_model_sort_ind_max[1:aasample_size[i]]] == "n", na.rm=T)/aa_tot
  aa_model_Score_perc[i] <-sum(aa_non_snp_all_shuffled[aa_perc_names[1:aasample_size[i]]] == "n", na.rm=T)/aa_tot
}

[,1:22]
boxplot.matrix(aa_random_Score, las = 2, xlab = "#predictions",
               ylab = "TPR", main = "BreastCancer_somatic_detection")
points(aa_model_Score_mean, col = 2, pch = 16, cex = 0.7)
points(aa_model_Score_median, col = 3, pch = 17, cex = 0.7)
points(aa_model_Score_max, col = 4, pch = 18, cex = 0.7)
points(aa_model_Score_perc, col = 5, pch = 19, cex = 0.7)
legend(x = "topleft",legend = c("mean", "median", "max", "vote"), 
       pch = c(16,17,18,19), fill = c(2,3,4,5), cex = 0.5)

boxplot.matrix(aa_random_Score[,1:72], las = 2, xlab = "#predictions",
               ylab = "TPR", main = "BreastCancer_somatic_detection", ylim = c(0, 0.08))
points(aa_model_Score_mean[1:72], col = 2, pch = 16, cex = 0.7)
points(aa_model_Score_median[1:72], col = 3, pch = 17, cex = 0.7)
points(aa_model_Score_max[1:72], col = 4, pch = 18, cex = 0.7)
points(aa_model_Score_perc[1:72], col = 5, pch = 19, cex = 0.7)
legend(x = "topleft",legend = c("mean", "median", "max", "vote"), 
       pch = c(16,17,18,19), fill = c(2,3,4,5), cex = 0.5)

boxplot.matrix(t(aa_model_Score_percentile))

aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_random_Score) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, ncol(aa_allmat), 2)){
  aa_allmat[, i] <- aa_random_Score[, aac]
  aa_allmat[, i+1] <- t(aa_model_Score_percentile)[, aac]
  aac <- aac + 1
}
boxplot.matrix(aa_allmat, col = rep(c(2,3), ncol(aa_random_Score)),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "BreastCancer_Somatic_prediction")
axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
     labels = aasample_size, las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)


sum(aa_any_change > 0.75)


######################################################################################################
# investigate the top somatic mutations

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_somatic_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC/somatic_pos_neg4.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_somatic <- readPeakFile(aa_somatic_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
start(aa_somatic) <- start(aa_somatic) - 1
#end(aa_snp) <- end(aa_snp) - 1

library(BSgenome.Hsapiens.UCSC.hg38)
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/COSMIC")
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

# removing duplicates
aa_somatic_df <- as.data.frame(aa_somatic)
aa_is_dup1 <- duplicated(aa_somatic_df)
aaremove <- which(aa_is_dup1 > 0)
aa_somatic <- aa_somatic[-c(aaremove)]
aa_somatic_df <- as.data.frame(aa_somatic)

aa_ref <- as.character(levels(aa_somatic_df$V22)[as.numeric(aa_somatic_df$V22)])
aa_alt <- as.character(levels(aa_somatic_df$V23)[as.numeric(aa_somatic_df$V23)])

aa_ref_Alt <- cbind(aa_ref, aa_alt)
#unique(aa_ref)
#unique(aa_alt)
aa_names <-  paste0("somatic_", c(1:length(aa_somatic)))
# fix nulls in reg and alt: in alt just set them to "", in ref find the ref sequence, replace ref eith ref sequence and update alt to include the ref seq
aa_ref_Alt2 <- aa_ref_Alt
for(i in 1:nrow(aa_ref_Alt)){
  if(aa_ref_Alt[i, 1] == "null"){
    aa_ref_Alt[i, 1] <- as.character(getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = aa_somatic[i]))
    aa_ref_Alt[i, 2] <- paste0(aa_ref_Alt[i, 1], aa_ref_Alt[i, 2])
  }else if(aa_ref_Alt[i, 2] == "null"){
    aa_ref_Alt[i, 2] <- ""
  }
}

aa_percModels <- rowSums(abs(aa_all_dif_cosmic_perc_4) > 0.001)/ncol(aa_all_dif_cosmic_perc_4)
aa_somat_index <- as.numeric(unlist(lapply(strsplit(rownames(aa_all_dif_cosmic_3), "_"), "[[", 3)))
aa_perc_names <- aa_somat_index[sort(aa_percModels,decreasing = T,index.return=T)$ix]


aa_eqtl_chr <- as.character(levels(aa_somatic_df$seqnames)[as.numeric(aa_somatic_df$seqnames)])[aa_perc_names[1:300]]
aa_eqtl_st <- aa_somatic_df$start[aa_perc_names[1:300]]
aa_eqtl_ref <- aa_ref_Alt[aa_perc_names[1:300], 1]
aa_eqtl_alt <- aa_ref_Alt[aa_perc_names[1:300], 2]
aa_names <- c("ESR1_2","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","RARA", "RELA","RUNX1","SP1","TFAP2C","YBX1")


aamin_LLR_low <- aa_exp18_minLLR - 0.1

aa_snp_inves_cosmic_300 <- list()
for(i in 1:length(aa_eqtl_alt)){
  print(i)
  aa_snp_inves_cosmic_300[[i]] <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                                          my_motifs = TF.motifs.Expanded_pseudo_exp12_t,
                                          snp_chr = aa_eqtl_chr[i],
                                          snp_start = aa_eqtl_st[i],
                                          snp_ref = aa_eqtl_ref[i],
                                          snp_alt = aa_eqtl_alt[i],
                                          min_LLR = aamin_LLR_low[aa_names])
  
}

names(aa_snp_inves_cosmic_300) <- aa_perc_names[1:300]
aarn <- unlist(lapply(aa_snp_inves_cosmic_300, nrow))
table(aarn)

aa_snp_inves_cosmic_300 <- aa_snp_inves_cosmic_300[!duplicated(names(aa_snp_inves_cosmic_300))]

aa_somatic_df
aa_snp_inves_cosmic_300_mat <- do.call(rbind, aa_snp_inves_cosmic_300)
aa_snp_inves_cosmic_300_mat <- aa_snp_inves_cosmic_300_mat[aa_snp_inves_cosmic_300_mat[, 4] != 0,]
#View(aa_snp_inves_3_mat)

rownames(aa_all_dif_cosmic_perc_4)[1:5]

aa_somatic_scores_noncoding <- as.numeric(levels(aa_somatic_df$V25)[as.numeric(aa_somatic_df$V25)])
#aa_somatic_scores_coding <- as.numeric(levels(aa_somatic_df$V26)[as.numeric(aa_somatic_df$V26)])
aa_somatic_snp_vs_no <- levels(aa_somatic_df$V24)[as.numeric(aa_somatic_df$V24)]

aan <- as.numeric(unlist(lapply(strsplit(rownames(aa_snp_inves_cosmic_300_mat), "\\."), "[[", 1)))
aaFATHMMscore <- aa_somatic_scores_noncoding[aan]
aa_somatic_snp_vs_no_curr <- aa_somatic_snp_vs_no[aan]
aa_somatic_snp_vs_no_curr[aa_somatic_snp_vs_no_curr == "n"] <- "somatic"
aa_somatic_snp_vs_no_curr[aa_somatic_snp_vs_no_curr == "y"] <- "normal"


aamacx <- numeric(length(aan) )
aamean <- numeric(length(aan) )
aa_all_dif_cosmic_perc_4_filtered <- aa_all_dif_cosmic_perc_4[!duplicated(rownames(aa_all_dif_cosmic_perc_4)),]
aaperc_num <- as.numeric(unlist(lapply(strsplit(rownames(aa_all_dif_cosmic_perc_4_filtered), "_"), "[[", 3)))

for(i in 1:length(aamacx)){
  aacur <- which(aaperc_num == aan[i])
  aamacx[i] <- max(abs(aa_all_dif_cosmic_perc_4_filtered[aacur, ])) * sign(aa_all_dif_cosmic_perc_4_filtered[aacur, which.max(abs(aa_all_dif_cosmic_perc_4_filtered[aacur, ]))])
  aamean[i] <- mean(abs(aa_all_dif_cosmic_perc_4_filtered[aacur, abs(aa_all_dif_cosmic_perc_4_filtered[aacur,]) > 0.001 ]))
}

aa_snp_inves_cosmic_300_mat <- cbind(aa_snp_inves_cosmic_300_mat, aamacx)
aa_snp_inves_cosmic_300_mat <- cbind(aa_snp_inves_cosmic_300_mat, aamean)
colnames(aa_snp_inves_cosmic_300_mat)[c(11,12)] <- c("max_percentile_change", "mean_percentile_change")


aa_dpn <- character(length(aan))
names(aa_dpn) <- aan
aan_somatic <- paste0("eqtl_somatic_", aan, "_1")
for(j in 1:length(aan_somatic)){
  for(i in 1:length(pos_cosmic_output_hal)){
    if(aan_somatic[j] %in% rownames(pos_cosmic_output_hal[[i]])){
      aa_dpn[j] <- names(pos_cosmic_output_hal)[i]
      print("########")
      print(aan_somatic[j])
      print(names(pos_cosmic_output_hal)[i])
      print("########")
    }
  }
  for(i in 1:length(neg_cosmic_output_hal)){
    if(aan_somatic[j] %in% rownames(neg_cosmic_output_hal[[i]])){
      aa_dpn[j] <- names(neg_cosmic_output_hal)[i]
      print("########")
      print(aan_somatic[j])
      print(names(neg_cosmic_output_hal)[i])
      print("########")
    }
  }
}

aa_snp_inves_cosmic_300_mat <- cbind(aa_snp_inves_cosmic_300_mat, aa_dpn)
aa_snp_inves_cosmic_300_mat$TF[aa_snp_inves_cosmic_300_mat$TF == "ESR1_2"] <- "ESR1"


aa_qtl_chip_ov <- Enhancer.ReMapchip.Overlap.byEnhancer$OverlapMat[aa_dpn,]
aa_qtl_chip_in <- Enhancer.ReMapchip.Overlap.byEnhancer$IntMat[aa_dpn,]

aa_qtl_chip_ovl_evid <- array(dim = nrow(aa_qtl_chip_ov))
aa_qtl_chip_int_evid <- array(dim =nrow(aa_qtl_chip_in))
for(i in 1:length(aa_qtl_chip_ovl_evid)){
  if(aa_snp_inves_cosmic_300_mat$TF[i] %in% colnames(aa_qtl_chip_ov)){
    aa_qtl_chip_ovl_evid[i] <- aa_qtl_chip_ov[i, aa_snp_inves_cosmic_300_mat$TF[i]]
    aa_qtl_chip_int_evid[i] <- aa_qtl_chip_in[i, aa_snp_inves_cosmic_300_mat$TF[i]]
  }
}

aa_snp_inves_cosmic_300_mat <- cbind(aa_snp_inves_cosmic_300_mat, cbind(aa_qtl_chip_ovl_evid, aa_qtl_chip_int_evid))
colnames(aa_snp_inves_cosmic_300_mat)[c(13:15)] <- c("enhancer", "chip_overlap_evidence", "chip_interaction_evidence")



#names(enhancer_gene_ovlap_interact$IntList[[1]]) <- names(aa_pos_neg)
#names(enhancer_gene_ovlap_interact$OverlapList[[1]]) <- names(aa_pos_neg)
names(enhancer_enhancer_interact$IntList[[1]]) <- names(aa_pos_neg)
#names(enhancer_enhancer_interact$OverlapList[[1]]) <- names(aa_pos_neg)

aa_ovl_gene <- character(length = nrow(aa_snp_inves_cosmic_300_mat))
aa_int_gene <- character(length = nrow(aa_snp_inves_cosmic_300_mat))
aa_int_enha <- character(length = nrow(aa_snp_inves_cosmic_300_mat))
aa_ovl_enha <- character(length = nrow(aa_snp_inves_cosmic_300_mat))
for(i in 1: nrow(aa_snp_inves_cosmic_300_mat)){
  aa_ovl_gene[i] <- paste(enhancer_gene_ovlap_interact$OverlapList[[1]][[aa_snp_inves_cosmic_300_mat$enhancer[i]]], collapse = "__")
  aa_int_gene[i] <- paste(enhancer_gene_ovlap_interact$IntList[[1]][[aa_snp_inves_cosmic_300_mat$enhancer[i]]], collapse = "__")
  aa_int_enha[i] <- paste(enhancer_enhancer_interact$IntList[[1]][[aa_snp_inves_cosmic_300_mat$enhancer[i]]], collapse = "__")
  #aa_ovl_enha[i] <- paste(enhancer_enhancer_interact$OverlapList[[1]][[aa_snp_inves_4_mat$enhancer[i]]], collapse = "__")
}
aa_snp_inves_cosmic_300_mat <- cbind(aa_snp_inves_cosmic_300_mat, cbind(cbind(aa_int_gene, aa_ovl_gene), cbind(aa_int_enha, aa_ovl_enha)))
colnames(aa_snp_inves_cosmic_300_mat)[16:19] <- c("interacting_Gene", "overlapping_Gene", "interacting_enhancer", "overlapping_enhancer")

# aars <- unlist(lapply(strsplit(unlist(lapply(strsplit(rownames(aa_snp_inves_cosmic_300_mat), "\\."), "[[", 1)), "_"), "[[", 2))
# aars_gtex <- aa_snp_qtl_marker_gtex_match_3[match(aars, names(aa_snp_qtl_marker_gtex_match_3))]
# aars_pancan <- aa_snp_qtl_marker_pancan_match_3[match(aars, names(aa_snp_qtl_marker_pancan_match_3))]
# aars_GWAS <- numeric(length(aars))
# names(aars_GWAS) <- aars
# aars_GWAS[aars %in% aa_GWAS_SNPS] <- 1
aa_somatic_snp_vs_no_curr 
aaFATHMMscore
aa_snp_inves_cosmic_300_mat <- cbind(aa_snp_inves_cosmic_300_mat, aaFATHMMscore, aa_somatic_snp_vs_no_curr)
colnames(aa_snp_inves_cosmic_300_mat)[20:21] <- c("FATHMM_noncoding", "Somatic_vs_normal")

write.csv(x = aa_snp_inves_cosmic_300_mat, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/eqtl_TF_somatic_cosmic_vote_400.csv")


table(aa_snp_inves_cosmic_300_mat$TF[aa_snp_inves_cosmic_300_mat$FATHMM_noncoding > 0.7])
table(aa_snp_inves_cosmic_300_mat$Somatic_vs_normal)
table(aa_snp_inves_cosmic_300_mat$Somatic_vs_normal[aa_snp_inves_cosmic_300_mat$FATHMM_noncoding > 0.7])

aawh <- which(aa_snp_inves_cosmic_300_mat$FATHMM_noncoding > 0.7 & (aa_snp_inves_cosmic_300_mat$chip_overlap_evidence > 0 | aa_snp_inves_cosmic_300_mat$chip_interaction_evidence > 0))
write.csv(x = aa_snp_inves_cosmic_300_mat[aawh,], file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/eqtl_TF_somatic_cosmic_vote_400_highFathmm_chip_evidence.csv")
View(aa_snp_inves_cosmic_300_mat[aawh,])


######################################################################################################
#######################################################################################################################
######################################################################################################
# clinvar mutants
aa_all_ClinVar <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/variant_summary.txt",
                             header = T, sep = "\t", stringsAsFactors = F)
aa_all_ClinVar <- aa_all_ClinVar[aa_all_ClinVar$Assembly == "GRCh38",]


aa_chr <- aa_all_ClinVar$Chromosome
aa_chr <- paste0("chr", aa_chr)
aa_sta <- aa_all_ClinVar$Start
aa_end <- aa_all_ClinVar$Stop

aa_clinvar_df <- data.frame(chr = aa_chr,
                           start = aa_sta, 
                           end = aa_end,
                           stringsAsFactors = F)

colnames(aa_all_ClinVar)

aa_clinvar_df <- cbind(aa_clinvar_df, aa_all_ClinVar[, c(c(1:18), c(22:31))])

aa_clinvar_gr <- makeGRangesFromDataFrame(aa_clinvar_df, keep.extra.columns = T)
options(scipen=999)
write.table(aa_clinvar_gr, 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/clinvar_hg38.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
options(scipen=0)

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_clinvar_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/clinvar_hg38.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_clinvar <- readPeakFile(aa_clinvar_Add, as = "GRanges")

start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
start(aa_clinvar) <- start(aa_clinvar) - 1

aa_ov_pos <- findOverlaps(query = aa_pos, subject = aa_clinvar)
length(unique(aa_ov_pos@to))
aa_ov_neg <- findOverlaps(query = aa_neg, subject = aa_clinvar)
length(unique(aa_ov_neg@from))
hist(as.numeric(table(aa_ov_pos@to)), breaks = 100)


aa_clinvar_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/clinvar_hg38.bed"
aaa_clinvar_pos4 <- bedtools_intersect(bedfile_names = c(aa_clinvar_Add, aa_pos_Add), 
                                       wa=F, wb=F, loj=F, wo=F, wao=F,
                                       u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                       S=F, sorted=F, merge_after_distance=0, 
                                       output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/clinvar_pos_4.bed", 
                                       read_output = T)
aaa_clinvar_neg4 <- bedtools_intersect(bedfile_names = c(aa_clinvar_Add, aa_neg_Add), 
                                       wa=F, wb=F, loj=F, wo=F, wao=F,
                                       u=T, c=F, v=F, f=1e-9, .F=1e-9, r=F, e=F, s=F,
                                       S=F, sorted=F, merge_after_distance=0, 
                                       output_filename="~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/clinvar_neg_4.bed", 
                                       read_output = T)
start(aaa_clinvar_pos4) <- start(aaa_clinvar_pos4) - 1
start(aaa_clinvar_neg4) <- start(aaa_clinvar_neg4) - 1

aa_clinvar_pos_neg_4 <- c(aaa_clinvar_pos4, aaa_clinvar_neg4)


options(scipen=999)

write.table(aa_clinvar_pos_neg_4, 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/clinvar_pos_neg4.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
options(scipen=0)


sum(paste0("rs", unique(aa_clinvar_pos_neg_4$V15)) %in% names(aa_snp_gwas_marker))
#######################################################################################
# write clinvar mutations

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_clinvar_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/clinvar_pos_neg4.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_clinvar <- readPeakFile(aa_clinvar_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
start(aa_clinvar) <- start(aa_clinvar) - 1
#end(aa_snp) <- end(aa_snp) - 1

library(BSgenome.Hsapiens.UCSC.hg38)
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/")
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

# removing duplicates
aa_clinvar_df <- as.data.frame(aa_clinvar)
aa_is_dup1 <- duplicated(aa_clinvar_df)
aaremove <- which(aa_is_dup1 > 0)
aa_clinvar <- aa_clinvar[-c(aaremove)]
aa_clinvar_df <- as.data.frame(aa_clinvar)
table(aa_clinvar_df$V9)
aaremove <- which(aa_clinvar_df$V26 == "na")
aa_clinvar <- aa_clinvar[-c(aaremove)]
aa_clinvar_df <- as.data.frame(aa_clinvar)
table(aa_clinvar_df$V9)
#summary(aa_clinvar_df$width[aa_clinvar_df$V9 == "duplication"])
options(scipen=999)
write.table(aa_clinvar, 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/clinvar_pos_neg4_onlyWith_RefAlt.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
options(scipen=0)



aa_ref <- as.character(levels(aa_clinvar_df$V26)[as.numeric(aa_clinvar_df$V26)])
aa_alt <- as.character(levels(aa_clinvar_df$V27)[as.numeric(aa_clinvar_df$V27)])
unique(unlist(strsplit(aa_ref, "")))
unique(unlist(strsplit(aa_alt, "")))

aa_ref_Alt <- cbind(aa_ref, aa_alt)
#unique(aa_ref)
#unique(aa_alt)
aa_names <-  paste0("clinvar_", c(1:length(aa_clinvar)))
# fixing for deletions: coordinates should be updated 
aaw <- which(aa_clinvar$V9 == "deletion")
start(aa_clinvar[aaw]) <- start(aa_clinvar[aaw]) - 1
aa_clinvar_df <- as.data.frame(aa_clinvar)


#aaaa <- unlist(lapply(strsplit(aa_ref_Alt[, 2], split = ","), length))

aa_clinvar[685]
getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = aa_clinvar[685])
# 
start(aa_clinvar)[131] <- start(aa_clinvar)[131] - 1
aa_clinvar_df <- as.data.frame(aa_clinvar)


start(aa_clinvar)[233] <- start(aa_clinvar)[233] - 1
aa_clinvar_df <- as.data.frame(aa_clinvar)

start(aa_clinvar)[685] <- start(aa_clinvar)[685] - 1
aa_clinvar_df <- as.data.frame(aa_clinvar)


aaw <- which(aa_clinvar$V9 == "short repeat")
aa_clinvar[aaw[7:length(aaw)]]

# i <- 20
# aa_clinvar[aaw[i]]$V26
# aa_clinvar[aaw[i]]$V27
# getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = aa_clinvar[aaw[i]])
# start(aa_clinvar)[aaw[i]] <- start(aa_clinvar)[aaw[i]] - 1

aa_clinvar_df <- as.data.frame(aa_clinvar)
sort(table(aa_clinvar_df$V22))

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/")
aa_pos_neg_vars_clinvar <- write_variant_seq(enhancer_GR = aa_pos_neg, 
                                         eqtl_GR = aa_clinvar,
                                         eqtl_names = aa_names,
                                         eqtl_ref_alt = aa_ref_Alt,
                                         all_combs = F,
                                         my_genome = BSgenome.Hsapiens.UCSC.hg38, 
                                         label = c(rep(1, length(aa_pos)), 
                                                   rep(0, length(aa_neg))))

# write GEMSTAT jobs to run
# write hal job for KDs
aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas1[1:150]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas2[1:150]]
aaadd3 <- union(aaadd1, aaadd2)

aa_par_name <- paste0(aaadd3, ".txt")
aa_seq_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/Variant_seq/",
                           recursive = T)
aa_lab_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/Variant_Labels/",
                           recursive = T)
aa_nam_spl <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_seq_files, 
                                                                   split = "\\."), 
                                                          "[[", 1)), split = "_"),
                                   "[", c(3,4)),
                            paste, collapse = "_"))


setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar")
for(i in 1:length(aa_par_name)){
  for(j in 1:length(aa_seq_files)){
    cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr ",
          paste0("-s Variant_seq/", aa_seq_files[j]), 
          paste0("-e Variant_Labels/", aa_lab_files[j]),
          "-m motifs.wtmx -f TF_exp.tab", 
          paste0("-fo Varaint_out/", aa_nam_spl[j], "_", aaadd3[i], ".out"), 
          "-o DIRECT -c Coop/coop.par ", 
          paste0("-p Trained_par/", aa_par_name[i]),
          "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1\n"
    ), 
    sep = " ", append = !(i==1 & j==1), file = "variant_jobs_clinvar_exp18.job")
  }
}
# reading results for Clinvar
load("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/clinvar_output_exp18_hal_WT.RData")
load("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/clinvar_output_exp18_hal.RData")
pos_clinvar_output_hal
neg_clinvar_output_hal
WT_models_clinvar_hal

aa_neg_dif_clinvar_3 <- list()
aa_pos_dif_clinvar_3 <- list()

# aa_neg_dif_cosmic_percentile_3 <- list()
# aa_pos_dif_cosmic_percentile_3 <- list()

aa_neg_dif_clinvar_percentile_4 <- list()
aa_pos_dif_clinvar_percentile_4 <- list()

# aa_Exp18_WT_model_prediction_Filtered_GT <- rownames(Exp18_WT_model_prediction_Filtered)
# aa_Exp18_WT_model_prediction_Filtered_GT <- unlist(lapply(strsplit(aa_Exp18_WT_model_prediction_Filtered_GT, split = "_"), "[[", 1))
# aa_Exp18_WT_model_prediction_Filtered_GT[aa_Exp18_WT_model_prediction_Filtered_GT == "pos"] <- 1
# aa_Exp18_WT_model_prediction_Filtered_GT[aa_Exp18_WT_model_prediction_Filtered_GT == "neg"] <- 0
# aa_Exp18_WT_model_prediction_Filtered_GT <- as.numeric(aa_Exp18_WT_model_prediction_Filtered_GT)
# names(aa_Exp16_WT_model_prediction_Filtered_GT) <- rownames(Exp16_WT_model_prediction_Filtered)

aa_WT_GT <- rownames(WT_models_clinvar_hal)
aa_WT_GT <- unlist(lapply(strsplit(aa_WT_GT, "_"), "[[", 1))
aa_WT_GT[aa_WT_GT  == "pos"] <- 1
aa_WT_GT[aa_WT_GT  == "neg"] <- 0
aa_WT_GT <- as.numeric(aa_WT_GT)
names(aa_WT_GT) <- rownames(WT_models_clinvar_hal)

for(i in 1:length(neg_clinvar_output_hal)){
  print(i)
  aa_neg_dif_clinvar_3[[i]] <- matrix(nrow = nrow(neg_clinvar_output_hal[[i]]), ncol = ncol(neg_clinvar_output_hal[[i]]))
  colnames(aa_neg_dif_clinvar_3[[i]]) <- colnames(neg_clinvar_output_hal[[i]])
  rownames(aa_neg_dif_clinvar_3[[i]]) <- rownames(neg_clinvar_output_hal[[i]])
  
  # aa_neg_dif_snp_percentile_3[[i]] <- matrix(nrow = nrow(neg_snp_output_hal[[i]]), ncol = ncol(neg_snp_output_hal[[i]]))
  # colnames(aa_neg_dif_snp_percentile_3[[i]]) <- colnames(neg_snp_output_hal[[i]])
  # rownames(aa_neg_dif_snp_percentile_3[[i]]) <- rownames(neg_snp_output_hal[[i]])
  
  aa_neg_dif_clinvar_percentile_4[[i]] <- matrix(nrow = nrow(neg_clinvar_output_hal[[i]]), ncol = ncol(neg_clinvar_output_hal[[i]]))
  colnames(aa_neg_dif_clinvar_percentile_4[[i]]) <- colnames(neg_clinvar_output_hal[[i]])
  rownames(aa_neg_dif_clinvar_percentile_4[[i]]) <- rownames(neg_clinvar_output_hal[[i]])
  
  for(j in 1:nrow(neg_clinvar_output_hal[[i]])){
    aa_neg_dif_clinvar_3[[i]][j, ] <- (neg_clinvar_output_hal[[i]][j, ] - WT_models_clinvar_hal[names(neg_clinvar_output_hal)[i], ])/ WT_models_clinvar_hal[names(neg_clinvar_output_hal)[i], ]
    for(kk in 1:ncol(neg_clinvar_output_hal[[i]])){
      #if(!is.na(Exp18_WT_model_prediction_Filtered[names(neg_snp_output_hal)[i],kk])){
      # aa_pos_perc_before <- sum(WT_models_snp_hal[aa_WT_GT == 1,kk] <= WT_models_snp_hal[names(neg_cosmic_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 1,kk]))
      # aa_neg_perc_before <- sum(WT_models_snp_hal[aa_WT_GT == 0,kk] >  WT_models_snp_hal[names(neg_cosmic_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 0,kk]))
      # aa_pos_perc_after <-  sum(WT_models_snp_hal[aa_WT_GT == 1,kk] <=  neg_cosmic_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 1,kk]))
      # aa_neg_perc_after <-  sum(WT_models_snp_hal[aa_WT_GT == 0,kk] >   neg_cosmic_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_snp_hal[aa_WT_GT == 0,kk]))
      # aa_neg_dif_snp_percentile_3[[i]][j,kk] <- -(aa_neg_perc_after - aa_neg_perc_before) + (aa_pos_perc_after - aa_pos_perc_before)
      aa_all_before <- sum(WT_models_clinvar_hal[,kk] <= WT_models_clinvar_hal[names(neg_clinvar_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_clinvar_hal[,kk]))
      aa_all_after <- sum(WT_models_clinvar_hal[,kk] <= neg_clinvar_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_clinvar_hal[,kk]))
      aa_neg_dif_clinvar_percentile_4[[i]][j,kk] <- aa_all_after - aa_all_before
      #}
    }
  }
}
names(aa_neg_dif_clinvar_3) <- names(neg_clinvar_output_hal)
#names(aa_neg_dif_snp_percentile_3) <- names(neg_cosmic_output_hal)
names(aa_neg_dif_clinvar_percentile_4) <- names(neg_clinvar_output_hal)

neg_clinvar_output_hal_all <- do.call(what = rbind, neg_clinvar_output_hal)
aa_neg_dif_all_clinvar_3 <- do.call(what = rbind, aa_neg_dif_clinvar_3)
#aa_neg_dif_snp_all_percentile_3 <- do.call(what = rbind, aa_neg_dif_snp_percentile_3)
aa_neg_dif_clinvar_all_percentile_4 <- do.call(what = rbind, aa_neg_dif_clinvar_percentile_4)



for(i in 1:length(pos_clinvar_output_hal)){
  print(i)
  aa_pos_dif_clinvar_3[[i]] <- matrix(nrow = nrow(pos_clinvar_output_hal[[i]]), ncol = ncol(pos_clinvar_output_hal[[i]]))
  colnames(aa_pos_dif_clinvar_3[[i]]) <- colnames(pos_clinvar_output_hal[[i]])
  rownames(aa_pos_dif_clinvar_3[[i]]) <- rownames(pos_clinvar_output_hal[[i]])
  
  # aa_pos_dif_snp_percentile_3[[i]] <- matrix(nrow = nrow(pos_cosmic_output_hal[[i]]), ncol = ncol(pos_cosmic_output_hal[[i]]))
  # colnames(aa_pos_dif_snp_percentile_3[[i]]) <- colnames(pos_cosmic_output_hal[[i]])
  # rownames(aa_pos_dif_snp_percentile_3[[i]]) <- rownames(pos_cosmic_output_hal[[i]])
  
  aa_pos_dif_clinvar_percentile_4[[i]] <- matrix(nrow = nrow(pos_clinvar_output_hal[[i]]), ncol = ncol(pos_clinvar_output_hal[[i]]))
  colnames(aa_pos_dif_clinvar_percentile_4[[i]]) <- colnames(pos_clinvar_output_hal[[i]])
  rownames(aa_pos_dif_clinvar_percentile_4[[i]]) <- rownames(pos_clinvar_output_hal[[i]])
  
  for(j in 1:nrow(pos_clinvar_output_hal[[i]])){
    aa_pos_dif_clinvar_3[[i]][j, ] <- (pos_clinvar_output_hal[[i]][j, ] - WT_models_clinvar_hal[names(pos_clinvar_output_hal)[i], ])/ WT_models_clinvar_hal[names(pos_clinvar_output_hal)[i], ]
    for(kk in 1:ncol(pos_clinvar_output_hal[[i]])){
      #if(!is.na(Exp18_WT_model_prediction_Filtered[names(pos_snp_output_hal)[i],kk])){
      # aa_pos_perc_before <- sum(WT_models_cosmic_hal[aa_WT_GT == 1,kk] <= WT_models_cosmic_hal[names(pos_snp_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[aa_WT_GT == 1,kk]))
      # aa_neg_perc_before <- sum(WT_models_cosmic_hal[aa_WT_GT == 0,kk] >  WT_models_cosmic_hal[names(pos_snp_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[aa_WT_GT == 0,kk]))
      # aa_pos_perc_after <-  sum(WT_models_cosmic_hal[aa_WT_GT == 1,kk] <=  pos_snp_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[aa_WT_GT == 1,kk]))
      # aa_neg_perc_after <-  sum(WT_models_cosmic_hal[aa_WT_GT == 0,kk] >   pos_snp_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_cosmic_hal[aa_WT_GT == 0,kk]))
      # aa_pos_dif_snp_percentile_3[[i]][j,kk] <- -(aa_neg_perc_after - aa_neg_perc_before) + (aa_pos_perc_after - aa_pos_perc_before)
      aa_all_before <- sum(WT_models_clinvar_hal[,kk] <= WT_models_clinvar_hal[names(pos_clinvar_output_hal)[i],kk], na.rm = T)/sum(!is.na(WT_models_clinvar_hal[,kk]))
      aa_all_after <- sum(WT_models_clinvar_hal[,kk] <= pos_clinvar_output_hal[[i]][j, kk], na.rm = T)/sum(!is.na(WT_models_clinvar_hal[,kk]))
      aa_pos_dif_clinvar_percentile_4[[i]][j,kk] <- aa_all_after - aa_all_before
      #}
    }
  }
}
names(aa_pos_dif_clinvar_3) <- names(pos_clinvar_output_hal)
#names(aa_pos_dif_snp_percentile_3) <- names(pos_cosmic_output_hal)
names(aa_pos_dif_clinvar_percentile_4) <- names(pos_clinvar_output_hal)

pos_clinvar_output_hal_all <- do.call(what = rbind, pos_clinvar_output_hal)
aa_pos_dif_all_clinvar_3 <- do.call(what = rbind, aa_pos_dif_clinvar_3)
#aa_pos_dif_snp_all_percentile_3 <- do.call(what = rbind, aa_pos_dif_snp_percentile_3)
aa_pos_dif_clinvar_all_percentile_4 <- do.call(what = rbind, aa_pos_dif_clinvar_percentile_4)


aa_pos_neg_out_clinvar_all_3 <- rbind(pos_clinvar_output_hal_all, neg_clinvar_output_hal_all)
aa_all_dif_clinvar_3 <- rbind(aa_pos_dif_all_clinvar_3, aa_neg_dif_all_clinvar_3)
#aa_all_dif_cosmic_perc_3 <- rbind(aa_pos_dif_cosmic_all_percentile_3, aa_neg_dif_cosmic_all_percentile_3)
aa_all_dif_clinvar_perc_4 <- rbind(aa_pos_dif_clinvar_all_percentile_4, aa_neg_dif_clinvar_all_percentile_4)


aa_all_dif_clinvar_names_3 <- unlist(lapply(strsplit(rownames(aa_all_dif_clinvar_3), split="_"), "[[", 3))

####### looking at proportion of pathogenic variants found by models against random
nrow(aa_all_dif_clinvar_3)
nrow(aa_all_dif_clinvar_perc_4)
aa_clinvar_index <- as.numeric(unlist(lapply(strsplit(rownames(aa_all_dif_clinvar_3), "_"), "[[", 3)))
aa_clinvar_patho <- aa_clinvar_df$V15


aa_model_index_Sorted_percentile <- matrix(nrow = nrow(aa_all_dif_clinvar_perc_4),
                                           ncol = ncol(aa_all_dif_clinvar_perc_4))
aa_model_index_Sorted_change <- matrix(nrow = nrow(aa_all_dif_clinvar_3),
                                       ncol = ncol(aa_all_dif_clinvar_3))
for(i in 1:ncol(aa_model_index_Sorted_percentile)){
  aa_model_index_Sorted_percentile[, i] <- aa_clinvar_index[sort(abs(aa_all_dif_clinvar_perc_4[, i]),
                                                               decreasing = T, index.return = T)$ix]
  aa_model_index_Sorted_change[, i] <- aa_clinvar_index[sort(abs(aa_all_dif_clinvar_3[, i]),
                                                           decreasing = T, index.return = T)$ix]
}
aa_max <- apply(abs(aa_all_dif_clinvar_perc_4), MARGIN = 1, FUN = max)
aa_mea <- apply(abs(aa_all_dif_clinvar_perc_4), MARGIN = 1, FUN = mean)
aa_med <- apply(abs(aa_all_dif_clinvar_perc_4), MARGIN = 1, FUN = median)

aa_model_sort_ind_max <- aa_clinvar_index[sort(aa_max, decreasing = T, index.return = T)$ix]
aa_model_sort_ind_mea <- aa_clinvar_index[sort(aa_mea, decreasing = T, index.return = T)$ix]
aa_model_sort_ind_med <- aa_clinvar_index[sort(aa_med, decreasing = T, index.return = T)$ix]

aa_percModels <- rowSums(abs(aa_all_dif_clinvar_perc_4) > 0.001)/ncol(aa_all_dif_clinvar_perc_4)
aa_perc_names <- aa_clinvar_index[sort(aa_percModels,decreasing = T,index.return=T)$ix]

aasample_size <- c(seq(1,100, 5), seq(110, length(aa_clinvar_index), 50))
aa_rep_num <- 244
aa_random_Score <- matrix(nrow = aa_rep_num, ncol = length(aasample_size))
colnames(aa_random_Score) <- aasample_size
aath <- 0.7

aa_model_Score_mean <- numeric(length = length(aasample_size))
aa_model_Score_median <- numeric(length = length(aasample_size))
aa_model_Score_max <- numeric(length = length(aasample_size))
aa_model_Score_perc <- numeric(length = length(aasample_size))

aa_model_Score_percentile <- matrix(nrow = length(aasample_size), ncol = ncol(aa_all_dif_clinvar_3))
colnames(aa_model_Score_percentile) <- colnames(aa_all_dif_clinvar_3)
rownames(aa_model_Score_percentile) <- aasample_size

aa_model_Score_change <- matrix(nrow = length(aasample_size), ncol = ncol(aa_all_dif_clinvar_3))
colnames(aa_model_Score_change) <- colnames(aa_all_dif_clinvar_3)
rownames(aa_model_Score_change) <- aasample_size

aa_tot <- sum(aa_clinvar_patho[aa_clinvar_index] == 1, na.rm=T)
for(i in 1:length(aasample_size)){
  for(j in 1:aa_rep_num){
    aasampl <- sample(x = aa_clinvar_patho[aa_clinvar_index], size = aasample_size[i], replace = F)
    aa_random_Score[j, i] <- sum(aasampl == 1, na.rm=T)/aa_tot
  }
  for(j in 1:ncol(aa_all_dif_clinvar_perc_4)){
    aa_model_Score_percentile[i, j] <- sum(aa_clinvar_patho[aa_model_index_Sorted_percentile[1:aasample_size[i],j]] == 1, na.rm=T)/aa_tot
    aa_model_Score_change[i, j] <- sum(aa_clinvar_patho[aa_model_index_Sorted_change[1:aasample_size[i],j]] == 1, na.rm=T)/aa_tot
    
  }
  aa_model_Score_mean[i] <- sum(aa_clinvar_patho[aa_model_sort_ind_mea[1:aasample_size[i]]] == 1, na.rm=T)/aa_tot
  aa_model_Score_median[i] <- sum(aa_clinvar_patho[aa_model_sort_ind_med[1:aasample_size[i]]] == 1, na.rm=T)/aa_tot
  aa_model_Score_max[i] <- sum(aa_clinvar_patho[aa_model_sort_ind_max[1:aasample_size[i]]] == 1, na.rm=T)/aa_tot
  aa_model_Score_perc[i] <-sum(aa_clinvar_patho[aa_perc_names[1:aasample_size[i]]] == 1, na.rm=T)/aa_tot
}

[,1:22]
boxplot.matrix(aa_random_Score, las = 2, xlab = "#predictions",
               ylab = "TPR", main = "clinvar_pathogenic")
points(aa_model_Score_mean, col = 2, pch = 16, cex = 0.7)
points(aa_model_Score_median, col = 3, pch = 17, cex = 0.7)
points(aa_model_Score_max, col = 4, pch = 18, cex = 0.7)
points(aa_model_Score_perc, col = 5, pch = 19, cex = 0.7)
legend(x = "topleft",legend = c("mean", "median", "max", "vote"), 
       pch = c(16,17,18,19), fill = c(2,3,4,5), cex = 0.5)

boxplot.matrix(aa_random_Score[,1:22], las = 2, xlab = "#predictions",
               ylab = "TPR", main = "clinvar_pathogenic")
points(aa_model_Score_mean[1:22], col = 2, pch = 16, cex = 0.7)
points(aa_model_Score_median[1:22], col = 3, pch = 17, cex = 0.7)
points(aa_model_Score_max[1:22], col = 4, pch = 18, cex = 0.7)
points(aa_model_Score_perc[1:22], col = 5, pch = 19, cex = 0.7)
legend(x = "topleft",legend = c("mean", "median", "max", "vote"), 
       pch = c(16,17,18,19), fill = c(2,3,4,5), cex = 0.5)

boxplot.matrix(t(aa_model_Score_percentile))

aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_random_Score) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, 94, 2)){
  aa_allmat[, i] <- aa_random_Score[, aac]
  aa_allmat[, i+1] <- t(aa_model_Score_percentile)[, aac]
  aac <- aac + 1
}
boxplot.matrix(aa_allmat, col = rep(c(2,3), 47),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "clinvar_pathogenic")
axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
     labels = aasample_size, las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)

boxplot.matrix(aa_allmat[,1:44], col = rep(c(2,3), 47),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "clinvar_pathogenic")
axis(side = 1, at = seq(1.5,ncol(aa_allmat[,1:44]), 2), 
     labels = aasample_size[1:22], las = 2)
abline(v = seq(0.5, ncol(aa_allmat[,1:44])+1, 2), col = 4)

aa_allmat2 <- matrix(nrow = 244, ncol = ncol(aa_random_Score) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, 94, 2)){
  aa_allmat2[, i] <- aa_random_Score[, aac]
  aa_allmat2[, i+1] <- t(aa_model_Score_change)[, aac]
  aac <- aac + 1
}
boxplot.matrix(aa_allmat2[,1:44], col = rep(c(2,3), 47),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "clinvar_pathogenic")
axis(side = 1, at = seq(1.5,ncol(aa_allmat2[,1:44]), 2), 
     labels = aasample_size[1:22], las = 2)
abline(v = seq(0.5, ncol(aa_allmat2[,1:44])+1, 2), col = 4)
###################################################
# repeat the above analysis with all common snps: Does it pick pathogenic ones more often?

aa_clinvar_index <- as.numeric(unlist(lapply(strsplit(rownames(aa_all_dif_clinvar_3), "_"), "[[", 3)))
aa_clinvar_patho <- aa_clinvar_df$V15[aa_clinvar_index]
aa_clinvar_patho_all <- c(aa_clinvar_patho, rep(0, nrow(aa_all_dif_snp_perc_4)))
aa_all_dif_clinvar_and_snp_perc_4 <- rbind(aa_all_dif_clinvar_perc_4,
                                           aa_all_dif_snp_perc_4)

aa_all_dif_clinvar_and_snp_perc_4_shuffind <- sample(c(1:nrow(aa_all_dif_clinvar_and_snp_perc_4)), 
                                                     size = nrow(aa_all_dif_clinvar_and_snp_perc_4),
                                                     replace = F)
aa_all_dif_clinvar_and_snp_perc_4_shuffled <- aa_all_dif_clinvar_and_snp_perc_4[aa_all_dif_clinvar_and_snp_perc_4_shuffind,]
aa_clinvar_patho_all_shuffled <- aa_clinvar_patho_all[aa_all_dif_clinvar_and_snp_perc_4_shuffind]

aa_model_index_Sorted_percentile <- matrix(nrow = nrow(aa_all_dif_clinvar_and_snp_perc_4),
                                           ncol = ncol(aa_all_dif_clinvar_and_snp_perc_4))
aa_model_index_Sorted_change <- matrix(nrow = nrow(aa_all_dif_clinvar_3),
                                       ncol = ncol(aa_all_dif_clinvar_3))
for(i in 1:ncol(aa_model_index_Sorted_percentile)){
  aa_model_index_Sorted_percentile[, i] <- sort(abs(aa_all_dif_clinvar_and_snp_perc_4_shuffled[, i]),
                                                                 decreasing = T, index.return = T)$ix
 # aa_model_index_Sorted_change[, i] <- sort(abs(aa_all_dif_clinvar_3[, i]),
  #                                                           decreasing = T, index.return = T)$ix
}
aa_max <- apply(abs(aa_all_dif_clinvar_and_snp_perc_4_shuffled), MARGIN = 1, FUN = max)
aa_mea <- apply(abs(aa_all_dif_clinvar_and_snp_perc_4_shuffled), MARGIN = 1, FUN = mean)
aa_med <- apply(abs(aa_all_dif_clinvar_and_snp_perc_4_shuffled), MARGIN = 1, FUN = median)

aa_model_sort_ind_max <- sort(aa_max, decreasing = T, index.return = T)$ix
aa_model_sort_ind_mea <- sort(aa_mea, decreasing = T, index.return = T)$ix
aa_model_sort_ind_med <- sort(aa_med, decreasing = T, index.return = T)$ix

aa_percModels <- rowSums(abs(aa_all_dif_clinvar_and_snp_perc_4_shuffled) > 0.001)/ncol(aa_all_dif_clinvar_and_snp_perc_4_shuffled)
aa_perc_names <- sort(aa_percModels,decreasing = T,index.return=T)$ix

#aasample_size <- c(seq(1,100, 5), seq(110, length(aa_clinvar_index), 50))
aasample_size <- c(seq(1,500, 10),
                   seq(550, 2000, 50),
                   seq(2000, 20000, 1000))
aa_rep_num <- 244
aa_random_Score <- matrix(nrow = aa_rep_num, ncol = length(aasample_size))
colnames(aa_random_Score) <- aasample_size
aath <- 0.7

aa_model_Score_mean <- numeric(length = length(aasample_size))
aa_model_Score_median <- numeric(length = length(aasample_size))
aa_model_Score_max <- numeric(length = length(aasample_size))
aa_model_Score_perc <- numeric(length = length(aasample_size))

aa_model_Score_percentile <- matrix(nrow = length(aasample_size), ncol = ncol(aa_all_dif_clinvar_3))
colnames(aa_model_Score_percentile) <- colnames(aa_all_dif_clinvar_3)
rownames(aa_model_Score_percentile) <- aasample_size

aa_model_Score_change <- matrix(nrow = length(aasample_size), ncol = ncol(aa_all_dif_clinvar_3))
colnames(aa_model_Score_change) <- colnames(aa_all_dif_clinvar_3)
rownames(aa_model_Score_change) <- aasample_size

aa_tot <- sum(aa_clinvar_patho_all_shuffled == 1, na.rm=T)
for(i in 1:length(aasample_size)){
  for(j in 1:aa_rep_num){
    aasampl <- sample(x = aa_clinvar_patho_all_shuffled, size = aasample_size[i], replace = F)
    aa_random_Score[j, i] <- sum(aasampl == 1, na.rm=T)/aa_tot
  }
  for(j in 1:ncol(aa_all_dif_clinvar_and_snp_perc_4)){
    aa_model_Score_percentile[i, j] <- sum(aa_clinvar_patho_all_shuffled[aa_model_index_Sorted_percentile[1:aasample_size[i],j]] == 1, na.rm=T)/aa_tot
    #aa_model_Score_change[i, j] <- sum(aa_clinvar_patho_all[aa_model_index_Sorted_change[1:aasample_size[i],j]] == 1, na.rm=T)/aa_tot
    
  }
  aa_model_Score_mean[i] <- sum(aa_clinvar_patho_all_shuffled[aa_model_sort_ind_mea[1:aasample_size[i]]] == 1, na.rm=T)/aa_tot
  aa_model_Score_median[i] <- sum(aa_clinvar_patho_all_shuffled[aa_model_sort_ind_med[1:aasample_size[i]]] == 1, na.rm=T)/aa_tot
  aa_model_Score_max[i] <- sum(aa_clinvar_patho_all_shuffled[aa_model_sort_ind_max[1:aasample_size[i]]] == 1, na.rm=T)/aa_tot
  aa_model_Score_perc[i] <-sum(aa_clinvar_patho_all_shuffled[aa_perc_names[1:aasample_size[i]]] == 1, na.rm=T)/aa_tot
}


[,1:22]
boxplot.matrix(aa_random_Score, las = 2, xlab = "#predictions",
               ylab = "TPR", main = "Clinvar_pathogenic_detection")
points(aa_model_Score_mean, col = 2, pch = 16, cex = 0.7)
points(aa_model_Score_median, col = 3, pch = 17, cex = 0.7)
points(aa_model_Score_max, col = 4, pch = 18, cex = 0.7)
points(aa_model_Score_perc, col = 5, pch = 19, cex = 0.7)
legend(x = "topleft",legend = c("mean", "median", "max", "vote"), 
       pch = c(16,17,18,19), fill = c(2,3,4,5), cex = 0.5)

boxplot.matrix(aa_random_Score[,1:72], las = 2, xlab = "#predictions",
               ylab = "TPR", main = "Clinvar_pathogenic_detection")
points(aa_model_Score_mean[1:72], col = 2, pch = 16, cex = 0.7)
points(aa_model_Score_median[1:72], col = 3, pch = 17, cex = 0.7)
points(aa_model_Score_max[1:72], col = 4, pch = 18, cex = 0.7)
points(aa_model_Score_perc[1:72], col = 5, pch = 19, cex = 0.7)
legend(x = "topleft",legend = c("mean", "median", "max", "vote"), 
       pch = c(16,17,18,19), fill = c(2,3,4,5), cex = 0.5)

boxplot.matrix(t(aa_model_Score_percentile))

aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_random_Score) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, ncol(aa_allmat), 2)){
  aa_allmat[, i] <- aa_random_Score[, aac]
  aa_allmat[, i+1] <- t(aa_model_Score_percentile)[, aac]
  aac <- aac + 1
}
boxplot.matrix(aa_allmat[,1:120], col = rep(c(2,3), ncol(aa_random_Score)),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "Clinvar_pathogenic_detection")
axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
     labels = aasample_size, las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)

boxplot.matrix(aa_allmat, col = rep(c(2,3), ncol(aa_random_Score)),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "Clinvar_pathogenic_detection")
axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
     labels = aasample_size, las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)

for(i in 1:20){
  print(i)
  print(sum(aa_clinvar_patho_all_shuffled[aa_perc_names[1:aasample_size[i]]] == 1, na.rm=T))
  print(sum(aa_clinvar_patho_all_shuffled[aa_model_sort_ind_mea[1:aasample_size[i]]] == 1, na.rm=T))
  print("#######")
}

##################################
# investigate the top clinvar mutations

aa_clinvar[1]

library(BSgenome.Hsapiens.UCSC.hg38)
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)


aa_ref <- as.character(levels(aa_clinvar_df$V26)[as.numeric(aa_clinvar_df$V26)])
aa_alt <- as.character(levels(aa_clinvar_df$V27)[as.numeric(aa_clinvar_df$V27)])
unique(unlist(strsplit(aa_ref, "")))
unique(unlist(strsplit(aa_alt, "")))

aa_ref_Alt <- cbind(aa_ref, aa_alt)

aa_percModels <- rowSums(abs(aa_all_dif_clinvar_perc_4) > 0.001)/ncol(aa_all_dif_clinvar_perc_4)
aa_clinvar_index <- as.numeric(unlist(lapply(strsplit(rownames(aa_all_dif_clinvar_3), "_"), "[[", 3)))
aa_perc_names <- aa_clinvar_index[sort(aa_percModels,decreasing = T,index.return=T)$ix]
table(aa_clinvar_df$V15[aa_perc_names[1:400]])

aa_eqtl_chr <- as.character(levels(aa_clinvar_df$seqnames)[as.numeric(aa_clinvar_df$seqnames)])[aa_perc_names[1:400]]
aa_eqtl_st <- aa_clinvar_df$start[aa_perc_names[1:400]]
aa_eqtl_ref <- aa_ref_Alt[aa_perc_names[1:400], 1]
aa_eqtl_alt <- aa_ref_Alt[aa_perc_names[1:400], 2]
aa_names <- c("ESR1_2","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","RARA", "RELA","RUNX1","SP1","TFAP2C","YBX1")


aamin_LLR_low <- aa_exp18_minLLR - 0.1

aa_snp_inves_clinvar_400 <- list()
for(i in 1:length(aa_eqtl_alt)){
  print(i)
  aa_snp_inves_clinvar_400[[i]] <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                                                   my_motifs = TF.motifs.Expanded_pseudo_exp12_t,
                                                   snp_chr = aa_eqtl_chr[i],
                                                   snp_start = aa_eqtl_st[i],
                                                   snp_ref = aa_eqtl_ref[i],
                                                   snp_alt = aa_eqtl_alt[i],
                                                   min_LLR = aamin_LLR_low[aa_names])
  
}

names(aa_snp_inves_clinvar_400) <- aa_perc_names[1:400]
aarn <- unlist(lapply(aa_snp_inves_clinvar_400, nrow))
table(aarn)




aa_snp_inves_clinvar_400 <- aa_snp_inves_clinvar_400[!duplicated(names(aa_snp_inves_clinvar_400))]

aa_clinvar_df
aa_clinvar_patho <- aa_clinvar_df$V15
aa_clinvar_patho_name <- as.character(levels(aa_clinvar_df$V14)[as.numeric(aa_clinvar_df$V14)])

aa_snp_inves_clinvar_400_mat <- do.call(rbind, aa_snp_inves_clinvar_400)
aa_snp_inves_clinvar_400_mat <- aa_snp_inves_clinvar_400_mat[aa_snp_inves_clinvar_400_mat[, 4] != 0,]
#View(aa_snp_inves_3_mat)

rownames(aa_all_dif_clinvar_perc_4)[1:5]


aan <- as.numeric(unlist(lapply(strsplit(rownames(aa_snp_inves_clinvar_400_mat), "\\."), "[[", 1)))
aa_clinvar_patho_curr <- aa_clinvar_patho[aan]
aa_clinvar_patho_name_curr <- aa_clinvar_patho_name[aan]


aamacx <- numeric(length(aan) )
aamean <- numeric(length(aan) )
aa_all_dif_clinvar_perc_4_filtered <- aa_all_dif_clinvar_perc_4[!duplicated(rownames(aa_all_dif_clinvar_perc_4)),]
aaperc_num <- as.numeric(unlist(lapply(strsplit(rownames(aa_all_dif_clinvar_perc_4_filtered), "_"), "[[", 3)))

for(i in 1:length(aamacx)){
  aacur <- which(aaperc_num == aan[i])
  aamacx[i] <- max(abs(aa_all_dif_clinvar_perc_4_filtered[aacur, ])) * sign(aa_all_dif_clinvar_perc_4_filtered[aacur, which.max(abs(aa_all_dif_clinvar_perc_4_filtered[aacur, ]))])
  aamean[i] <- mean(abs(aa_all_dif_clinvar_perc_4_filtered[aacur, abs(aa_all_dif_clinvar_perc_4_filtered[aacur,]) > 0.001 ]))
}

aa_snp_inves_clinvar_400_mat <- cbind(aa_snp_inves_clinvar_400_mat, aamacx)
aa_snp_inves_clinvar_400_mat <- cbind(aa_snp_inves_clinvar_400_mat, aamean)
colnames(aa_snp_inves_clinvar_400_mat)[c(11,12)] <- c("max_percentile_change", "mean_percentile_change")


aa_dpn <- character(length(aan))
names(aa_dpn) <- aan
aan_somatic <- paste0("eqtl_clinvar_", aan, "_1")
for(j in 1:length(aan_somatic)){
  for(i in 1:length(pos_clinvar_output_hal)){
    if(aan_somatic[j] %in% rownames(pos_clinvar_output_hal[[i]])){
      aa_dpn[j] <- names(pos_clinvar_output_hal)[i]
      print("########")
      print(aan_somatic[j])
      print(names(pos_clinvar_output_hal)[i])
      print("########")
    }
  }
  for(i in 1:length(neg_clinvar_output_hal)){
    if(aan_somatic[j] %in% rownames(neg_clinvar_output_hal[[i]])){
      aa_dpn[j] <- names(neg_clinvar_output_hal)[i]
      print("########")
      print(aan_somatic[j])
      print(names(neg_clinvar_output_hal)[i])
      print("########")
    }
  }
}

aa_snp_inves_clinvar_400_mat <- cbind(aa_snp_inves_clinvar_400_mat, aa_dpn)
aa_snp_inves_clinvar_400_mat$TF[aa_snp_inves_clinvar_400_mat$TF == "ESR1_2"] <- "ESR1"


aa_qtl_chip_ov <- Enhancer.ReMapchip.Overlap.byEnhancer$OverlapMat[aa_dpn,]
aa_qtl_chip_in <- Enhancer.ReMapchip.Overlap.byEnhancer$IntMat[aa_dpn,]

aa_qtl_chip_ovl_evid <- array(dim = nrow(aa_qtl_chip_ov))
aa_qtl_chip_int_evid <- array(dim =nrow(aa_qtl_chip_in))
for(i in 1:length(aa_qtl_chip_ovl_evid)){
  if(aa_snp_inves_clinvar_400_mat$TF[i] %in% colnames(aa_qtl_chip_ov)){
    aa_qtl_chip_ovl_evid[i] <- aa_qtl_chip_ov[i, aa_snp_inves_clinvar_400_mat$TF[i]]
    aa_qtl_chip_int_evid[i] <- aa_qtl_chip_in[i, aa_snp_inves_clinvar_400_mat$TF[i]]
  }
}

aa_snp_inves_clinvar_400_mat <- cbind(aa_snp_inves_clinvar_400_mat, cbind(aa_qtl_chip_ovl_evid, aa_qtl_chip_int_evid))
colnames(aa_snp_inves_clinvar_400_mat)[c(13:15)] <- c("enhancer", "chip_overlap_evidence", "chip_interaction_evidence")



#names(enhancer_gene_ovlap_interact$IntList[[1]]) <- names(aa_pos_neg)
#names(enhancer_gene_ovlap_interact$OverlapList[[1]]) <- names(aa_pos_neg)
names(enhancer_enhancer_interact$IntList[[1]]) <- names(aa_pos_neg)
#names(enhancer_enhancer_interact$OverlapList[[1]]) <- names(aa_pos_neg)

aa_ovl_gene <- character(length = nrow(aa_snp_inves_clinvar_400_mat))
aa_int_gene <- character(length = nrow(aa_snp_inves_clinvar_400_mat))
aa_int_enha <- character(length = nrow(aa_snp_inves_clinvar_400_mat))
aa_ovl_enha <- character(length = nrow(aa_snp_inves_clinvar_400_mat))
for(i in 1: nrow(aa_snp_inves_clinvar_400_mat)){
  aa_ovl_gene[i] <- paste(enhancer_gene_ovlap_interact$OverlapList[[1]][[aa_snp_inves_clinvar_400_mat$enhancer[i]]], collapse = "__")
  aa_int_gene[i] <- paste(enhancer_gene_ovlap_interact$IntList[[1]][[aa_snp_inves_clinvar_400_mat$enhancer[i]]], collapse = "__")
  aa_int_enha[i] <- paste(enhancer_enhancer_interact$IntList[[1]][[aa_snp_inves_clinvar_400_mat$enhancer[i]]], collapse = "__")
  #aa_ovl_enha[i] <- paste(enhancer_enhancer_interact$OverlapList[[1]][[aa_snp_inves_4_mat$enhancer[i]]], collapse = "__")
}
aa_snp_inves_clinvar_400_mat <- cbind(aa_snp_inves_clinvar_400_mat, cbind(cbind(aa_int_gene, aa_ovl_gene), cbind(aa_int_enha, aa_ovl_enha)))
colnames(aa_snp_inves_clinvar_400_mat)[16:19] <- c("interacting_Gene", "overlapping_Gene", "interacting_enhancer", "overlapping_enhancer")

aa_clinvar_patho_name_curr 
aa_clinvar_patho_curr
aa_snp_inves_clinvar_400_mat <- cbind(aa_snp_inves_clinvar_400_mat, aa_clinvar_patho_name_curr, aa_clinvar_patho_curr)
colnames(aa_snp_inves_clinvar_400_mat)[20:21] <- c("pathogen_status", "pathogen_simple")

write.csv(x = aa_snp_inves_clinvar_400_mat, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/eqtl_TF_somatic_clinvar_vote_400.csv")


table(aa_snp_inves_clinvar_400_mat$TF[aa_snp_inves_clinvar_400_mat$pathogen_simple == 1])
table(aa_snp_inves_clinvar_400_mat$Somatic_vs_normal)
table(aa_snp_inves_clinvar_400_mat$Somatic_vs_normal[aa_snp_inves_clinvar_400_mat$FATHMM_noncoding > 0.7])

aawh <- which(aa_snp_inves_clinvar_400_mat$pathogen_simple == 1 & (aa_snp_inves_clinvar_400_mat$chip_overlap_evidence > 0 | aa_snp_inves_clinvar_400_mat$chip_interaction_evidence > 0))

View(aa_snp_inves_clinvar_400_mat[aawh,])
write.csv(x = aa_snp_inves_clinvar_400_mat[aawh,], file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/eqtl_TF_somatic_clinvar_vote_400_patho_withChIP.csv")


#################
#PharmGKB
# var_drug_ann
# var_fa_ann
# var_pheno_ann
aa_pharm_drug <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/PharmGKB/annotations/var_drug_ann.tsv", stringsAsFactors = F)
which(aa_pharm_drug$Variant %in% names(aa_snp_qtl_marker_gtex))

aa_pharm_fa <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/PharmGKB/annotations/var_fa_ann.tsv", stringsAsFactors = F)
which(aa_pharm_fa$Variant %in% names(aa_snp_qtl_marker_gtex))

aa_pharm_pheno <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/PharmGKB/annotations/var_pheno_ann.tsv", stringsAsFactors = F)
which(aa_pharm_pheno$Variant %in% names(aa_snp_qtl_marker_gtex))
aa_pharm_pheno[which(aa_pharm_pheno$Variant %in% names(aa_snp_qtl_marker_gtex)),]
#################

# HGMD variants
# can't download the public version of the dataset
#################
# drug related

# get the ensemble id of genes that overlap or interact with some gene
(Promoter1kbdf_gene$gene_id[1:5]) 
aa_ovlgenes <- unique(unlist(enhancer_gene_ovlap_interact$OverlapList[[1]]))
aa_intgenes <- unlist(enhancer_gene_ovlap_interact$IntList[[1]])
aa_intgenes <- aa_intgenes[!duplicated(aa_intgenes)]
aa_ovlgenes <- Promoter1kbdf_gene$gene_id[aa_ovlgenes]
aa_intgenes <- Promoter1kbdf_gene$gene_id[aa_intgenes]
aa_all_genes <- union(aa_ovlgenes,aa_intgenes )

names(enhancer_gene_ovlap_interact$OverlapList[[1]]) <- names(aa_pos_neg)
names(enhancer_gene_ovlap_interact$IntList[[1]]) <- names(aa_pos_neg)

aa_ovlgenes <- unlist(enhancer_gene_ovlap_interact$OverlapList[[1]])
aa_intgenes <- unlist(enhancer_gene_ovlap_interact$IntList[[1]])
aa_ovlgenes_name <- Promoter1kbdf_gene$gene_id[aa_ovlgenes]
names(aa_ovlgenes_name) <- names(aa_ovlgenes)
aa_intgenes_name <- Promoter1kbdf_gene$gene_id[aa_intgenes]
names(aa_intgenes_name) <- names(aa_intgenes)
aa_all_genes_names <- c(aa_ovlgenes_name, aa_intgenes_name)
## 
# read drug genes
aa_files_full <- list.files("~/Documents/Shayan/BioInf/EstrogenReceptor/drug_response/", full.names = T)
aa_files <- list.files("~/Documents/Shayan/BioInf/EstrogenReceptor/drug_response/", full.names = F)
aa_drug_response <- list()
for(i in 1:length(aa_files_full)){
  print(i)
  aa_drug_response[[i]] <- read.delim(file = aa_files_full[i], stringsAsFactors = F)
  aa_drug_response[[i]] <- aa_drug_response[[i]][aa_drug_response[[i]]$PVAL <= 0.05,]
  aa_drug_response[[i]]$ENSEMBL_new <- as.character(as.numeric(unlist(lapply(strsplit(aa_drug_response[[i]]$ENSEMBL, "ENSG"), "[[", 2))))
}
names(aa_drug_response) <- aa_files

aa_drug_response_related <- list()
for(i in 1:length(aa_drug_response)){
  aa_drug_response_related[[i]] <- aa_drug_response[[i]][aa_drug_response[[i]]$ENSEMBL_new %in% aa_all_genes, ]
}
names(aa_drug_response_related) <- names(aa_drug_response)
unlist(lapply(aa_drug_response_related, nrow))
# find enhancers realted to each drug
aa_drug_response_related_all <- do.call(rbind, aa_drug_response_related)
aa_enh_per_gene <- list()
for(i in 1:nrow(aa_drug_response_related_all)){
  aa_enh_per_gene[[i]] <- names(aa_all_genes_names)[aa_all_genes_names %in% aa_drug_response_related_all$ENSEMBL_new[i]]
}
names(aa_enh_per_gene) <- paste0(rownames(aa_drug_response_related_all), "_", aa_drug_response_related_all$ENSEMBL_new)

# find snps that are related to each drug
aa_snp_per_gene <- list()
aa_eqt_per_gene <- list()
for(i in 1:length(aa_enh_per_gene)){
  aa_snp_per_gene[[i]] <- unlist(enhancer_snp_list[names(enhancer_snp_list) %in% aa_enh_per_gene[[i]]])
  aa_eqt_per_gene[[i]] <- unlist(enhancer_eqtl_list[names(enhancer_eqtl_list) %in% aa_enh_per_gene[[i]]])
}
names(aa_snp_per_gene) <- names(aa_enh_per_gene)
names(aa_eqt_per_gene) <- names(aa_enh_per_gene)

aa_dr <- unlist(lapply(strsplit(names(aa_enh_per_gene), "\\."), "[[", 1))
aa_dr_uniq <- unique(aa_dr)
aa_drug_snp <- list()
aa_drug_eqtl <- list()

for(i in 1:length(aa_dr_uniq)){
  aa_drug_snp[[i]] <- unlist(aa_snp_per_gene[aa_dr %in% aa_dr_uniq[i]])
  #aa_drug_snp[[i]] <- aa_drug_snp[[i]][!duplicated(aa_drug_snp[[i]])]
  aa_drug_eqtl[[i]] <- unlist(aa_eqt_per_gene[aa_dr %in% aa_dr_uniq[i]])
  #aa_drug_eqtl[[i]] <- aa_drug_eqtl[[i]][!duplicated(aa_drug_eqtl[[i]])]
}
names(aa_drug_snp) <- aa_dr_uniq
names(aa_drug_eqtl) <- aa_dr_uniq


# find the relative rank of the snps that are related to each drug
aa_max_row <- apply(X = abs(aa_all_dif_snp_perc_4), MARGIN = 1, FUN = max)
aa_med_row <- apply(X = abs(aa_all_dif_snp_perc_4), MARGIN = 1, FUN = median)
aa_mea_row <- apply(X = abs(aa_all_dif_snp_perc_4), MARGIN = 1, FUN = mean)
aa_rowname_sorted_max <- rownames(aa_all_dif_snp_perc_4)[sort(aa_max_row,
                                                              decreasing = T,index.return=T)$ix]
aa_rowname_sorted_mea <- rownames(aa_all_dif_snp_perc_4)[sort(aa_mea_row, 
                                                              decreasing = T,index.return=T)$ix]
aa_rowname_sorted_med <- rownames(aa_all_dif_snp_perc_4)[sort(aa_med_row, 
                                                              decreasing = T,index.return=T)$ix]

aa_rowname_sorted_max_sp <- unique(unlist(lapply(strsplit(aa_rowname_sorted_max, "_"), "[[", 2)))
aa_rowname_sorted_mea_sp <- unique(unlist(lapply(strsplit(aa_rowname_sorted_mea, "_"), "[[", 2)))
aa_rowname_sorted_med_sp <- unique(unlist(lapply(strsplit(aa_rowname_sorted_med, "_"), "[[", 2)))

aa_percModels <- rowSums(abs(aa_all_dif_snp_perc_4) > 0.001)/ncol(aa_all_dif_snp_perc_4)
aa_perc_names <- rownames(aa_all_dif_snp_perc_4)[sort(aa_percModels,
                                                      decreasing = T,index.return=T)$ix]

aa_perc_names_sp <- unique(unlist(lapply(strsplit(aa_perc_names, "_"), "[[", 2)))

aa_drug_snp_percentile_mean <- list()
aa_drug_snp_percentile_max <- list()
aa_drug_snp_percentile_med <- list()
aa_drug_snp_percentile_vote <- list()


for(i in 1:length(aa_drug_snp)){
  aa_drug_snp_percentile_mean[[i]] <- numeric(length = length(aa_drug_snp[[i]]))
  aa_drug_snp_percentile_max[[i]] <- numeric(length = length(aa_drug_snp[[i]]))
  aa_drug_snp_percentile_med[[i]] <- numeric(length = length(aa_drug_snp[[i]]))
  aa_drug_snp_percentile_vote[[i]] <- numeric(length = length(aa_drug_snp[[i]]))
  if(length(aa_drug_snp[[i]]) > 0){
    for(j in 1:length(aa_drug_snp[[i]])){
      aa_drug_snp_percentile_mean[[i]][j] <- 1 - which(aa_rowname_sorted_mea_sp == aa_drug_snp[[i]][j])/length(aa_rowname_sorted_mea_sp)
      aa_drug_snp_percentile_max[[i]][j] <- 1 -which(aa_rowname_sorted_max_sp == aa_drug_snp[[i]][j])/length(aa_rowname_sorted_mea_sp)
      aa_drug_snp_percentile_med[[i]][j] <- 1 - which(aa_rowname_sorted_med_sp == aa_drug_snp[[i]][j])/length(aa_rowname_sorted_mea_sp)
      aa_drug_snp_percentile_vote[[i]][j] <- 1 - which(aa_perc_names_sp == aa_drug_snp[[i]][j])/length(aa_rowname_sorted_mea_sp)
    }
  }
  
}

max(unlist(lapply(aa_drug_snp, length)))

aa_mat_max <- matrix(nrow = max(unlist(lapply(aa_drug_snp, length))),
                 ncol = length(aa_drug_snp))
aa_mat_mea <- matrix(nrow = max(unlist(lapply(aa_drug_snp, length))),
                     ncol = length(aa_drug_snp))
aa_mat_med <- matrix(nrow = max(unlist(lapply(aa_drug_snp, length))),
                     ncol = length(aa_drug_snp))
aa_mat_vot <- matrix(nrow = max(unlist(lapply(aa_drug_snp, length))),
                     ncol = length(aa_drug_snp))

colnames(aa_mat_max) <- paste0(names(aa_drug_snp), "__" ,unlist(lapply(aa_drug_snp, length)))
colnames(aa_mat_mea) <- paste0(names(aa_drug_snp), "__" ,unlist(lapply(aa_drug_snp, length)))
colnames(aa_mat_med) <- paste0(names(aa_drug_snp), "__" ,unlist(lapply(aa_drug_snp, length)))
colnames(aa_mat_vot) <- paste0(names(aa_drug_snp), "__" ,unlist(lapply(aa_drug_snp, length)))

for(i in 1:ncol(aa_mat_max)){
  if(length(aa_drug_snp[[i]]) > 0){
    aa_mat_max[1:length(aa_drug_snp[[i]]), i] <- aa_drug_snp_percentile_max[[i]]
    aa_mat_mea[1:length(aa_drug_snp[[i]]), i] <- aa_drug_snp_percentile_mean[[i]]
    aa_mat_med[1:length(aa_drug_snp[[i]]), i] <- aa_drug_snp_percentile_med[[i]]
    aa_mat_vot[1:length(aa_drug_snp[[i]]), i] <- aa_drug_snp_percentile_vote[[i]]
  }
  
}
par(mfrow = c(1,1), mar = c(10,6,4,4))
boxplot.matrix(aa_mat_max, las = 2, ylab = "percentile among common SNPs\nthat fall under ER enhancers", main = "max_over_ensemble")
boxplot.matrix(aa_mat_mea, las = 2, ylab = "percentile among common SNPs\nthat fall under ER enhancers", main = "mean_over_ensemble")
boxplot.matrix(aa_mat_med, las = 2, ylab = "percentile among common SNPs\nthat fall under ER enhancers", main = "median_over_ensemble")
boxplot.matrix(aa_mat_vot, las = 2, ylab = "percentile among common SNPs\nthat fall under ER enhancers", main = "vote_over_ensemble")

#names(aa_drug_snp)
aa_breast_use <- c("not_used", "not_used", "not_used", "not_used", "used", "used", "not_used", "used", "used", "used", "not_used", "used", "used", "used", "used", "used", "not_used", "not_used", "used", "used", "used", "not_used", "used")
aa_cancer_use <- c("cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "non_cancer", "cancer", "cancer", "non_cancer", "cancer", "cancer", "cancer", "cancer", "non_cancer", "cancer")
names(aa_breast_use) <- names(aa_drug_snp)
names(aa_cancer_use) <- names(aa_drug_snp)

aa_col <- aa_breast_use
aa_col[aa_col == "not_used"] <- 2
aa_col[aa_col == "used"] <- 3
par(mfrow = c(1,1), mar = c(10,6,4,4))
boxplot.matrix(aa_mat_max, las = 2, col = aa_col, ylab = "percentile among common SNPs\nthat fall under ER enhancers", main = "max_over_ensemble")
boxplot.matrix(aa_mat_mea, las = 2, col = aa_col, ylab = "percentile among common SNPs\nthat fall under ER enhancers", main = "mean_over_ensemble")
boxplot.matrix(aa_mat_med, las = 2, col = aa_col, ylab = "percentile among common SNPs\nthat fall under ER enhancers", main = "median_over_ensemble")
boxplot.matrix(aa_mat_vot, las = 2, col = aa_col, ylab = "percentile among common SNPs\nthat fall under ER enhancers", main = "vote_over_ensemble")



###################################################
# repeat drug response analysis at gene level\



###################################################
# create a dataframe with one row for each snp, enhancer, gene connection
#    add drug info, ranks by models , chip info, knockdown info,
#    repsonsible TF, change in motif score, snp position, enhancer position
#  first: create a dataframe, with one row for each snp-enhancer connection, add genome coordinates of each
# include snps and somatic mutations





enhancer_snp_list


aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_snp_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/all_snps_pos_neg4.bed"
aa_eqtl_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/sig_eqtl.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
aa_snp <- readPeakFile(aa_snp_Add, as = "GRanges")
aa_eqtl <- readPeakFile(aa_eqtl_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
start(aa_snp) <- start(aa_snp) - 2
end(aa_snp) <- end(aa_snp) - 1

names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

# removing duplicates
aa_snpdf <- as.data.frame(aa_snp)
aa_is_dup1 <- duplicated(aa_snpdf[,1:(ncol(aa_snpdf) - 1)])
aaremove <- which(aa_is_dup1 > 0)
aa_snp <- aa_snp[-c(aaremove)]
aa_snpdf <- as.data.frame(aa_snp)
aa_snp <- makeGRangesFromDataFrame(aa_snpdf, keep.extra.columns = T)
Somatic_pos_neg <- aa_somatic_df
Common_SNP_pos_neg <- aa_snpdf

Somatic_pos_neg_GR <- aa_somatic
Common_SNP_pos_neg_GR <- aa_snp

aa_overlap_holder <- findOverlaps(query = aa_pos_neg, subject = aa_snp)
aa_ol_enhancers_index <- unique(aa_overlap_holder@from)

enhancer_snp_list <- list()
enhancer_eqtl_list <- list()
for(i in 1:length(aa_pos_neg)){
  if(i %in% aa_ol_enhancers_index){
    enhancer_snp_list[[i]] <- unique(aa_names[aa_overlap_holder@to[aa_overlap_holder@from %in% i]])
    enhancer_eqtl_list[[i]] <- enhancer_snp_list[[i]][aa_snp_qtl_marker_gtex_pancan[match(enhancer_snp_list[[i]], names(aa_snp_qtl_marker_gtex_pancan))] == 1]
  }else{
    enhancer_snp_list[[i]] <- character(0)
    enhancer_eqtl_list[[i]] <- character(0)
  }
  
}
names(enhancer_snp_list) <- names(aa_pos_neg)
names(enhancer_eqtl_list) <- names(aa_pos_neg)

Common_SNP_pos_neg_GR
GTEX_eqtl_pos4_neg4
ER_pancan_cis_trans_gr38_pos_neg_4
Somatic_pos_neg_GR



aa_Q <- GTEX_eqtl_pos4_neg4
ranges(aa_Q) <- ranges(aa_Q) + 1
aa_S <- Common_SNP_pos_neg_GR
ranges(aa_S) <- ranges(aa_S) + 1
aa_overlap_holder <- findOverlaps(query = aa_Q, subject = aa_S)

length(unique(aa_overlap_holder@from))
length(unique(aa_overlap_holder@to))

# second: create a dataframe, with one row per enhancer-gene connection: note if the connection is through proximity or overlap, or interaction, also note the source of interaction
# next merge the above two data frames

length(unique(ER_pancan_cis_trans_gr38_pos_neg_4))

#################

# create a combined dataframe starting with all common
aa_initial <- as.data.frame(Common_SNP_pos_neg_GR)
aa_initial$V9 <- as.character(levels(aa_initial$V9)[as.numeric(aa_initial$V9)])
aa_initial$V8 <- as.character(levels(aa_initial$V8)[as.numeric(aa_initial$V8)])
aa_initial$V10 <- as.character(levels(aa_initial$V10)[as.numeric(aa_initial$V10)])

aaasp <- strsplit(aa_initial$V9, ",")
aalla <- unlist(lapply(aaasp, length))
# first expand to put each alternative allele in a separate row
aa_initial_expanded <- aa_initial[rep(row.names(aa_initial), aalla),]
for(i in 1:nrow(aa_initial)){
  if(aalla[i] > 1){
    aa_initial_expanded[row.names(aa_initial)[i], "V9"] <- aaasp[[i]][1]
    for(j in 2:aalla[i]){
      aa_initial_expanded[paste(row.names(aa_initial)[i], j-1, sep = "."), "V9"] <- aaasp[[i]][j]
    }
  }
}
aa_initial_expanded$seqnames <- as.character(levels(aa_initial_expanded$seqnames)[as.numeric(aa_initial_expanded$seqnames)])
aa_initial_expanded_2 <- aa_initial_expanded
# add GTEX_eqtl_pos4_neg4 to aa_initial_expanded
GTEX_eqtl_pos4_neg4_df <- as.data.frame(GTEX_eqtl_pos4_neg4)

GTEX_eqtl_pos4_neg4_df$seqnames <- as.character(levels(GTEX_eqtl_pos4_neg4_df$seqnames)[as.numeric(GTEX_eqtl_pos4_neg4_df$seqnames)])
GTEX_eqtl_pos4_neg4_df$V6 <- as.character(levels(GTEX_eqtl_pos4_neg4_df$V6)[as.numeric(GTEX_eqtl_pos4_neg4_df$V6)])
GTEX_eqtl_pos4_neg4_df$V7 <- as.character(levels(GTEX_eqtl_pos4_neg4_df$V7)[as.numeric(GTEX_eqtl_pos4_neg4_df$V7)])
GTEX_eqtl_pos4_neg4_df$V8 <- as.character(levels(GTEX_eqtl_pos4_neg4_df$V8)[as.numeric(GTEX_eqtl_pos4_neg4_df$V8)])
GTEX_eqtl_pos4_neg4_df$GTEX_index <-  c(1:nrow(GTEX_eqtl_pos4_neg4_df))


aa_initial_expanded <- aa_initial_expanded_2
aa_initial_expanded <- aa_initial_expanded[,c(1,2,3,10,11,12)]
# aa_initial_expanded$GTEX_index <- numeric(length = nrow(aa_initial_expanded))
# aa_initial_expanded$Gene_cor <- character(length = nrow(aa_initial_expanded))
names(aa_initial_expanded)[c(4,5,6)] <- c("ref", "alt", "SNP_names")



# GTEX_index is the index of the eqtl in  GTEX_eqtl_pos4_neg4_df
aa_GTEX_eqtl_pos4_neg4_df <- GTEX_eqtl_pos4_neg4_df
aa_GTEX_eqtl_pos4_neg4_df <- aa_GTEX_eqtl_pos4_neg4_df[,c(1,2,3,8,9,10, 21)]
names(aa_GTEX_eqtl_pos4_neg4_df)[c(4,5,6)] <- c("ref", "alt", "Gene_cor_Gtex")




# for(i in 1:200){
#   print(i)
#   aaw <- which(aa_initial_expanded[, c(1,2,3,6,7)] %in% aa_GTEX_eqtl_pos4_neg4_df[i,c(1,2,3,6,7)]) 
#   if(length(aaw) == 1){
#     print("found")
#     aa_initial_expanded$GTEX_index[aaw] <- c(aa_initial_expanded$GTEX_index[aaw], i)
#     aa_initial_expanded$Gene_cor[aaw] <- c(aa_initial_expanded$Gene_cor[aaw], aa_GTEX_eqtl_pos4_neg4_df$V8[i])
#   }else if(length(aaw) == 0){
#     aa_initial_expanded <- rbind(aa_initial_expanded, aa_GTEX_eqtl_pos4_neg4_df[i, c(1,2,3,4,5,6,7,9,10,8)])
#   }
#   
# }




aa_initial_expanded$common_initial_expanded_id <- c(1:nrow(aa_initial_expanded))
colnames(aa_initial_expanded)

aa_common_gtex <- full_join(x = aa_initial_expanded,
                            y = aa_GTEX_eqtl_pos4_neg4_df,
                            by = c("seqnames", "start", "end", "ref" , "alt"))
# next add pancan:
ER_pancan_cis_trans_gr38_pos_neg_4df <- as.data.frame(ER_pancan_cis_trans_gr38_pos_neg_4)
head(ER_pancan_cis_trans_gr38_pos_neg_4df$V9)
ER_pancan_cis_trans_gr38_pos_neg_4df$seqnames <- as.character(levels(ER_pancan_cis_trans_gr38_pos_neg_4df$seqnames)[as.numeric(ER_pancan_cis_trans_gr38_pos_neg_4df$seqnames)])
ER_pancan_cis_trans_gr38_pos_neg_4df$V6 <- as.character(levels(ER_pancan_cis_trans_gr38_pos_neg_4df$V6)[as.numeric(ER_pancan_cis_trans_gr38_pos_neg_4df$V6)])
ER_pancan_cis_trans_gr38_pos_neg_4df$V7 <- as.character(levels(ER_pancan_cis_trans_gr38_pos_neg_4df$V7)[as.numeric(ER_pancan_cis_trans_gr38_pos_neg_4df$V7)])
ER_pancan_cis_trans_gr38_pos_neg_4df$V8 <- as.character(levels(ER_pancan_cis_trans_gr38_pos_neg_4df$V8)[as.numeric(ER_pancan_cis_trans_gr38_pos_neg_4df$V8)])
ER_pancan_cis_trans_gr38_pos_neg_4df$V9 <- as.character(levels(ER_pancan_cis_trans_gr38_pos_neg_4df$V9)[as.numeric(ER_pancan_cis_trans_gr38_pos_neg_4df$V9)])
ER_pancan_cis_trans_gr38_pos_neg_4df$pancan_index <-  c(1:nrow(ER_pancan_cis_trans_gr38_pos_neg_4df))
names(ER_pancan_cis_trans_gr38_pos_neg_4df)[c(8,9,10,11)] <- c( "ref" , "alt", "snp_name_pancan", "Gene_cor_pancan")
aaER_pancan_cis_trans_gr38_pos_neg_4df <- ER_pancan_cis_trans_gr38_pos_neg_4df
aaER_pancan_cis_trans_gr38_pos_neg_4df <- ER_pancan_cis_trans_gr38_pos_neg_4df[,c(1,2,3,8,9,10,11,12)]
aa_common_gtex_pancan <- full_join(x = aa_common_gtex,
                            y = aaER_pancan_cis_trans_gr38_pos_neg_4df,
                            by = c("seqnames", "start", "end", "ref" , "alt"))
# pancan_index is index in this: ER_pancan_cis_trans_gr38_pos_neg_4df

# next add somatic:
head(Somatic_pos_neg)
colnames(Somatic_pos_neg)[10:32] <- colnames(ER_somatic_mutation_all)[c(1,2,3,4,5,8,9,10,11,12,13,14,15,17,18,19,20,21,23,25,26,27,28)]

aa_ref <- as.character(levels(aa_somatic_df$V22)[as.numeric(aa_somatic_df$V22)])
aa_alt <- as.character(levels(aa_somatic_df$V23)[as.numeric(aa_somatic_df$V23)])

aa_ref_Alt <- cbind(aa_ref, aa_alt)
#unique(aa_ref)
#unique(aa_alt)
#aa_names <-  paste0("somatic_", c(1:length(aa_somatic)))
# fix nulls in reg and alt: in alt just set them to "", in ref find the ref sequence, replace ref eith ref sequence and update alt to include the ref seq
aa_ref_Alt2 <- aa_ref_Alt
for(i in 1:nrow(aa_ref_Alt)){
  if(aa_ref_Alt[i, 1] == "null"){
    aa_ref_Alt[i, 1] <- as.character(getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = aa_somatic[i]))
    aa_ref_Alt[i, 2] <- paste0(aa_ref_Alt[i, 1], aa_ref_Alt[i, 2])
  }else if(aa_ref_Alt[i, 2] == "null"){
    aa_ref_Alt[i, 2] <- ""
  }
}

Somatic_pos_neg$ref <- aa_ref_Alt[, 1]
Somatic_pos_neg$alt <- aa_ref_Alt[, 2]
Somatic_pos_neg$somatic_index <- c(1:nrow(Somatic_pos_neg))
aaSomatic_pos_neg <- Somatic_pos_neg[, c(1,2,3,23,26,27,33,34,35)]

colnames(aaSomatic_pos_neg)[c(5)] <- c("SNP_somatic_y_n")
aaSomatic_pos_neg$seqnames <- as.character(levels(aaSomatic_pos_neg$seqnames)[as.numeric(aaSomatic_pos_neg$seqnames)])
aaSomatic_pos_neg$MUTATION_SOMATIC_STATUS <- as.character(levels(aaSomatic_pos_neg$MUTATION_SOMATIC_STATUS)[as.numeric(aaSomatic_pos_neg$MUTATION_SOMATIC_STATUS)])
aaSomatic_pos_neg$SNP_somatic_y_n <- as.character(levels(aaSomatic_pos_neg$SNP_somatic_y_n)[as.numeric(aaSomatic_pos_neg$SNP_somatic_y_n)])
aaSomatic_pos_neg$FATHMM_MKL_NON_CODING_SCORE <- as.numeric(levels(aaSomatic_pos_neg$FATHMM_MKL_NON_CODING_SCORE)[as.numeric(aaSomatic_pos_neg$FATHMM_MKL_NON_CODING_SCORE)])

aa_common_gtex_pancan_somatic <- full_join(x = aa_common_gtex_pancan,
                                           y = aaSomatic_pos_neg,
                                           by = c("seqnames", "start", "end", "ref" , "alt"))


# next add clinvar:

Clinvar_pos_neg_df <- aa_clinvar_df
names(Clinvar_pos_neg_df)[10:37] <- colnames(aa_all_ClinVar)[c(c(1:18), c(22:31))]
Clinvar_pos_neg_df$Clinvar_index <- c(1:nrow(Clinvar_pos_neg_df))
aaClinvar_pos_neg_df <- Clinvar_pos_neg_df[,c(1,2,3,11,13,16,17,19,28,29,38)]
names(aaClinvar_pos_neg_df)[c(4,5,8,9,10)] <- c("Type_clinvar", "GeneID_clinvar", "SNP_name_Clinvar", "ref", "alt")
aaClinvar_pos_neg_df$seqnames <- as.character(levels(aaClinvar_pos_neg_df$seqnames)[as.numeric(aaClinvar_pos_neg_df$seqnames)])
aaClinvar_pos_neg_df$Type_clinvar <- as.character(levels(aaClinvar_pos_neg_df$Type_clinvar)[as.numeric(aaClinvar_pos_neg_df$Type_clinvar)])
aaClinvar_pos_neg_df$ClinicalSignificance <- as.character(levels(aaClinvar_pos_neg_df$ClinicalSignificance)[as.numeric(aaClinvar_pos_neg_df$ClinicalSignificance)])
aaClinvar_pos_neg_df$ClinicalSignificance <- as.character(levels(aaClinvar_pos_neg_df$ClinicalSignificance)[as.numeric(aaClinvar_pos_neg_df$ClinicalSignificance)])
aaClinvar_pos_neg_df$ref <- as.character(levels(aaClinvar_pos_neg_df$ref)[as.numeric(aaClinvar_pos_neg_df$ref)])
aaClinvar_pos_neg_df$alt <- as.character(levels(aaClinvar_pos_neg_df$alt)[as.numeric(aaClinvar_pos_neg_df$alt)])

aaClinvar_pos_neg_df$Clinvar_index[1:5]

# aa_clinvar_df <- cbind(aa_clinvar_df, aa_all_ClinVar[, c(c(1:18), c(22:31))])

aa_common_gtex_pancan_somatic_clinical <- full_join(x = aa_common_gtex_pancan_somatic,
                                           y = aaClinvar_pos_neg_df,
                                           by = c("seqnames", "start", "end", "ref" , "alt"))


ER_SNPs_merged <- aa_common_gtex_pancan_somatic_clinical
nrow(ER_SNPs_merged)



aamg <- colnames(ER_SNPs_merged)[6:ncol(ER_SNPs_merged)]
ER_SNPs_merged_allbypos <- aggregate(x = ER_SNPs_merged[aamg],
                                            by = ER_SNPs_merged[c("seqnames", "start", "end", "ref", "alt")], 
                                            FUN = c)
ER_SNPs_merged_allbypos$pos_uniq_id <- c(1:nrow(ER_SNPs_merged_allbypos))
# fix some mistakes in the file:
ER_SNPs_merged_allbypos$ref[30208] <- "A"
ER_SNPs_merged_allbypos$alt[30208] <- "C"
ER_SNPs_merged_allbypos$ref[19610] <- "T"
ER_SNPs_merged_allbypos$alt[19610] <- "C"


# check if ref and alt are provided correctly
aacheck_ref <- numeric(length = nrow(ER_SNPs_merged_allbypos))

#aacheck_alt <- numeric(length = nrow(ER_SNPs_merged_allbypos))
5694
aaxind <- which(!is.na(ER_SNPs_merged_allbypos$snp_name_pancan))
aacheck_ref2 <- numeric(length = length(aaxind))
for(i in 1:length(aaxind)){
  print(i)
  aaseq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = ER_SNPs_merged_allbypos_GR[aaxind[i]], as.character = T)
  aaseq <- unlist(strsplit(aaseq, ""))[1]
  aa_tss <- unlist(strsplit(ER_SNPs_merged_allbypos$ref[aaxind[i]], ""))[1]
  if(aaseq == aa_tss){
    aacheck_ref2[i] <- 1
  }else{
    print("not equal")
  }
  
}
# fix the ones where ref is wrong:
aawrong <- aaxind[aacheck_ref2 == 0]
ER_SNPs_merged_allbypos$ref[aawrong[1]] <- "T"
ER_SNPs_merged_allbypos$alt[aawrong[1]] <- "G"

ER_SNPs_merged_allbypos[aawrong[2],]
ER_SNPs_merged_allbypos$ref[aawrong[2]] <- "G"
ER_SNPs_merged_allbypos$alt[aawrong[2]] <- "A"

ER_SNPs_merged_allbypos[aawrong[3],]
ER_SNPs_merged_allbypos$ref[aawrong[3]] <- "C"
ER_SNPs_merged_allbypos$alt[aawrong[3]] <- "A"

ER_SNPs_merged_allbypos[aawrong[4],]
ER_SNPs_merged_allbypos$ref[aawrong[4]] <- "C"
ER_SNPs_merged_allbypos$alt[aawrong[4]] <- "A"

ER_SNPs_merged_allbypos[aawrong[5],]
ER_SNPs_merged_allbypos$start[aawrong[5]] <- ER_SNPs_merged_allbypos$start[aawrong[5]] - 3
ER_SNPs_merged_allbypos$end[aawrong[5]] <- ER_SNPs_merged_allbypos$end[aawrong[5]] - 1
ER_SNPs_merged_allbypos$ref[aawrong[5]] <- "TCA"
ER_SNPs_merged_allbypos$alt[aawrong[5]] <- "T"

ER_SNPs_merged_allbypos[aawrong[6],]
ER_SNPs_merged_allbypos$ref[aawrong[6]] <- "T"
ER_SNPs_merged_allbypos$alt[aawrong[6]] <- "C"

ER_SNPs_merged_allbypos[aawrong[7],]
ER_SNPs_merged_allbypos$ref[aawrong[7]] <- "G"
ER_SNPs_merged_allbypos$alt[aawrong[7]] <- "A"

ER_SNPs_merged_allbypos[aawrong[8],]
ER_SNPs_merged_allbypos$ref[aawrong[8]] <- "C"
ER_SNPs_merged_allbypos$alt[aawrong[8]] <- "T"

ER_SNPs_merged_allbypos[aawrong[9],]
ER_SNPs_merged_allbypos$ref[aawrong[9]] <- "C"
ER_SNPs_merged_allbypos$alt[aawrong[9]] <- "A"

ER_SNPs_merged_allbypos[aawrong[10],]
ER_SNPs_merged_allbypos$ref[aawrong[10]] <- "G"
ER_SNPs_merged_allbypos$alt[aawrong[10]] <- "A"

ER_SNPs_merged_allbypos[aawrong[11],]
ER_SNPs_merged_allbypos$ref[aawrong[11]] <- "C"
ER_SNPs_merged_allbypos$alt[aawrong[11]] <- "T"

ER_SNPs_merged_allbypos[aawrong[12],]
ER_SNPs_merged_allbypos$ref[aawrong[12]] <- "G"
ER_SNPs_merged_allbypos$alt[aawrong[12]] <- "A"

ER_SNPs_merged_allbypos[aawrong[13],]
ER_SNPs_merged_allbypos$ref[aawrong[13]] <- "C"
ER_SNPs_merged_allbypos$alt[aawrong[13]] <- "T"

ER_SNPs_merged_allbypos[aawrong[14],]
ER_SNPs_merged_allbypos$ref[aawrong[14]] <- "G"
ER_SNPs_merged_allbypos$alt[aawrong[14]] <- "A"

ER_SNPs_merged_allbypos[aawrong[15],]
ER_SNPs_merged_allbypos$ref[aawrong[15]] <- "G"
ER_SNPs_merged_allbypos$alt[aawrong[15]] <- "C"

ER_SNPs_merged_allbypos[aawrong[16],]
ER_SNPs_merged_allbypos$ref[aawrong[16]] <- "G"
ER_SNPs_merged_allbypos$alt[aawrong[16]] <- "A"

ER_SNPs_merged_allbypos[aawrong[17],]
ER_SNPs_merged_allbypos$ref[aawrong[17]] <- "G"
ER_SNPs_merged_allbypos$alt[aawrong[17]] <- "A"

ER_SNPs_merged_allbypos[aawrong[18],]
ER_SNPs_merged_allbypos$ref[aawrong[18]] <- "A"
ER_SNPs_merged_allbypos$alt[aawrong[18]] <- "G"

ER_SNPs_merged_allbypos[aawrong[19],]
ER_SNPs_merged_allbypos$ref[aawrong[19]] <- "A"
ER_SNPs_merged_allbypos$alt[aawrong[19]] <- "G"

ER_SNPs_merged_allbypos[aawrong[20],]
ER_SNPs_merged_allbypos$ref[aawrong[20]] <- "A"
ER_SNPs_merged_allbypos$alt[aawrong[20]] <- "G"

ER_SNPs_merged_allbypos[aawrong[21],]
ER_SNPs_merged_allbypos$ref[aawrong[21]] <- "G"
ER_SNPs_merged_allbypos$alt[aawrong[21]] <- "A"

ER_SNPs_merged_allbypos[aawrong[22],]
ER_SNPs_merged_allbypos$ref[aawrong[22]] <- "A"
ER_SNPs_merged_allbypos$alt[aawrong[22]] <- "G"

ER_SNPs_merged_allbypos[aawrong[23],]
ER_SNPs_merged_allbypos$ref[aawrong[23]] <- "G"
ER_SNPs_merged_allbypos$alt[aawrong[23]] <- "A"

ER_SNPs_merged_allbypos[aawrong[24],]
ER_SNPs_merged_allbypos$ref[aawrong[24]] <- "A"
ER_SNPs_merged_allbypos$alt[aawrong[24]] <- "G"

ER_SNPs_merged_allbypos[aawrong[25],]
ER_SNPs_merged_allbypos$ref[aawrong[25]] <- "A"
ER_SNPs_merged_allbypos$alt[aawrong[25]] <- "G"

ER_SNPs_merged_allbypos[aawrong[26],]
ER_SNPs_merged_allbypos$ref[aawrong[26]] <- "T"
ER_SNPs_merged_allbypos$alt[aawrong[26]] <- "C"

ER_SNPs_merged_allbypos[aawrong[27],]
ER_SNPs_merged_allbypos$ref[aawrong[27]] <- "G"
ER_SNPs_merged_allbypos$alt[aawrong[27]] <- "A"

ER_SNPs_merged_allbypos[aawrong[28],]
ER_SNPs_merged_allbypos$ref[aawrong[28]] <- "T"
ER_SNPs_merged_allbypos$alt[aawrong[28]] <- "C"

ER_SNPs_merged_allbypos[aawrong[29],]
ER_SNPs_merged_allbypos$ref[aawrong[29]] <- "T"
ER_SNPs_merged_allbypos$alt[aawrong[29]] <- "C"

ER_SNPs_merged_allbypos[aawrong[30],]
ER_SNPs_merged_allbypos$ref[aawrong[30]] <- "A"
ER_SNPs_merged_allbypos$alt[aawrong[30]] <- "C"

aatsss <- ER_SNPs_merged_allbypos_GR[aawrong[5]]
start(aatsss) <- start(aatsss) - 3
end(aatsss) <- end(aatsss) - 1
getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = aatsss, as.character = T)
############### change entires with multiple entries to their unique values if possible
aatstt <- ER_SNPs_merged_allbypos

aatstt$SNP_names <- lapply(ER_SNPs_merged_allbypos$SNP_names, unique)
aatstt$common_initial_expanded_id <-  lapply(ER_SNPs_merged_allbypos$common_initial_expanded_id, unique)
aatstt$Gene_cor_Gtex <- lapply(ER_SNPs_merged_allbypos$Gene_cor_Gtex, unique)
aatstt$GTEX_index <- lapply(ER_SNPs_merged_allbypos$GTEX_index, unique)
aatstt$snp_name_pancan <- lapply(ER_SNPs_merged_allbypos$snp_name_pancan, unique)
aatstt$Gene_cor_pancan <- lapply(ER_SNPs_merged_allbypos$Gene_cor_pancan, unique)
aatstt$pancan_index <- lapply(ER_SNPs_merged_allbypos$pancan_index, unique)
aatstt$MUTATION_SOMATIC_STATUS <- lapply(ER_SNPs_merged_allbypos$MUTATION_SOMATIC_STATUS, unique)
aatstt$SNP_somatic_y_n <- lapply(ER_SNPs_merged_allbypos$SNP_somatic_y_n, unique)
#aatstt$FATHMM_MKL_NON_CODING_SCORE
aatstt$somatic_index <- lapply(ER_SNPs_merged_allbypos$somatic_index, unique)
aatstt$Type_clinvar <- lapply(ER_SNPs_merged_allbypos$Type_clinvar, unique)
aatstt$GeneID_clinvar <- lapply(ER_SNPs_merged_allbypos$GeneID_clinvar, unique)
aatstt$ClinicalSignificance <- lapply(ER_SNPs_merged_allbypos$ClinicalSignificance, unique)
aatstt$ClinSigSimple <- lapply(ER_SNPs_merged_allbypos$ClinSigSimple, unique)
aatstt$SNP_name_Clinvar <- lapply(ER_SNPs_merged_allbypos$SNP_name_Clinvar, unique)
aatstt$Clinvar_index <- lapply(ER_SNPs_merged_allbypos$Clinvar_index, unique)
aaxxx <- aatstt$FATHMM_MKL_NON_CODING_SCORE
aaxxx2 <- (lapply(aaxxx, as.character))
aaxxx2 <- lapply(aaxxx2, unique)
aaxxx2 <- lapply(aaxxx2, as.numeric)
aatstt$FATHMM_MKL_NON_CODING_SCORE <- aaxxx2
aatstt$is_gtex <- !is.na(aatstt$GTEX_index)
aatstt$is_pancan <- !is.na(aatstt$pancan_index)
aatstt$is_somatic <- !is.na(aatstt$somatic_index)
aatstt$is_somatic_only <- aatstt$SNP_somatic_y_n == "n"
aatstt$is_somatic_gt_FATHMM0_7 <- aatstt$FATHMM_MKL_NON_CODING_SCORE >= 0.7
aatstt$is_clinical_pathogenic <- aatstt$ClinSigSimple == 1

ER_SNPs_merged_allbypos_modified <- aatstt
ER_SNPs_merged_allbypos_modified$is_clinical <- !is.na(ER_SNPs_merged_allbypos_modified$ClinSigSimple)
ER_SNPs_merged_allbypos_modified$is_common <- !is.na(ER_SNPs_merged_allbypos_modified$common_initial_expanded_id)
ER_SNPs_merged_allbypos_modified$is_fathmm_scored <- !is.na(aatstt$is_somatic_gt_FATHMM0_7)
# write variants for all the snps and evaluate them all in one place
aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1

library(BSgenome.Hsapiens.UCSC.hg38)

names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

# removing duplicates

aa_ref <- ER_SNPs_merged_allbypos$ref
aa_alt <- ER_SNPs_merged_allbypos$alt
ER_SNPs_merged_allbypos_GR <- makeGRangesFromDataFrame(ER_SNPs_merged_allbypos[, 1:3])


aa_ref_Alt <- cbind(aa_ref, aa_alt)
aa_names <-  ER_SNPs_merged_allbypos$pos_uniq_id
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/Final_variant_Set")
aa_pos_neg_vars_merged <- write_variant_seq(enhancer_GR = aa_pos_neg, 
                                            eqtl_GR = ER_SNPs_merged_allbypos_GR,
                                            eqtl_names = ER_SNPs_merged_allbypos$pos_uniq_id,
                                            eqtl_ref_alt = aa_ref_Alt,
                                            all_combs = F,
                                            my_genome = BSgenome.Hsapiens.UCSC.hg38, 
                                            label = c(rep(1, length(aa_pos)), rep(0, length(aa_neg))))


# write hal jobs for all the variants:
aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas1[1:150]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas2[1:150]]
aaadd3 <- union(aaadd1, aaadd2)

aa_par_name <- paste0(aaadd3, ".txt")
aa_seq_files <- list.files(path = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Final_variant_Set/Variant_seq/",
                           recursive = T)
aa_lab_files <- list.files(path = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Final_variant_Set/Variant_Labels/",
                           recursive = T)
aa_nam_spl <- unlist(lapply(lapply(strsplit(unlist(lapply(strsplit(aa_seq_files, 
                                                                   split = "\\."), 
                                                          "[[", 1)), split = "_"),
                                   "[", c(3,4)),
                            paste, collapse = "_"))


setwd("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Final_variant_Set")
for(i in 1:length(aa_par_name)){
  for(j in 1:length(aa_seq_files)){
    cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr ",
          paste0("-s Variant_seq/", aa_seq_files[j]), 
          paste0("-e Variant_Labels/", aa_lab_files[j]),
          "-m motifs.wtmx -f TF_exp.tab", 
          paste0("-fo Varaint_out/", aa_nam_spl[j], "_", aaadd3[i], ".out"), 
          "-o DIRECT -c Coop/coop.par ", 
          paste0("-p Trained_par/", aa_par_name[i]),
          "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1\n"
    ), 
    sep = " ", append = !(i==1 & j==1), file = "variant_jobs_merged_all_exp18.job")
  }
}
######################
# read the results



########################################################################################
load("~/Documents/Shayan/BioInf/EstrogenReceptor/Final_variant_Set/Merged_all_output_exp18_hal_WT.RData")
load("~/Documents/Shayan/BioInf/EstrogenReceptor/Final_variant_Set/Merged_all_output_exp18_hal.RData")
pos_mergedVar_output_hal
neg_mergedVar_output_hal
WT_models_merged_all_hal
Merged_percentile_change <- compute_percentile_change(WT_var_from_hal = WT_models_merged_all_hal,
                                                      Var_from_hal_neg = neg_mergedVar_output_hal, 
                                                      Var_from_hal_pos = pos_mergedVar_output_hal)
###################################################################################################
######################################################################################## gtex
aalab_orig <- ER_SNPs_merged_allbypos_modified$is_gtex
aa_som_rem <- which(ER_SNPs_merged_allbypos_modified$is_somatic)
aa_cli_rem <- which(ER_SNPs_merged_allbypos_modified$is_clinical)

aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[! aasap %in% union(aa_som_rem,aa_cli_rem),]
aasap <- aasap[! aasap %in% union(aa_som_rem,aa_cli_rem)]
aaTP <- aalab_orig[aasap]




aa_gtex_ornot <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
                                   rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
                                   TP = aaTP, 
                                   sampling_range=seq(100,5000, 100),#c(seq(10,300,10), seq(400,2000, 100), seq(3000,33000,1000)),
                                   do_boxplot=F)


# boxplot.matrix(aa_gtex_ornot$Random_scores[,10:50], las = 2, xlab = "#predictions",
#                ylab = "TPR", main = "gtex ranking")
# points(aa_gtex_ornot$mean_ranked_score[10:60], col = 2, pch = 16, cex = 0.7)
# points(aa_gtex_ornot$median_ranked_score[10:60], col = 3, pch = 17, cex = 0.7)
# points(aa_gtex_ornot$max_ranked_score[10:60], col = 4, pch = 18, cex = 0.7)
# points(aa_gtex_ornot$percent_ranked_score[10:60], col = 5, pch = 19, cex = 0.7)
# legend(x = "topleft",legend = c("mean", "median", "max", "vote"), 
#        pch = c(16,17,18,19), fill = c(2,3,4,5), cex = 0.5)


aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_gtex_ornot$Random_scores) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, ncol(aa_allmat), 2)){
  aa_allmat[, i] <- aa_gtex_ornot$Random_scores[, aac]
  aa_allmat[, i+1] <- t(aa_gtex_ornot$Model_scores)[, aac]
  aac <- aac + 1
}


aas <- sort(c(seq(1,100, 6),seq(2,100, 6)))

boxplot.matrix(aa_allmat[,aas], col = rep(c(2,3), length(aas)/2),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "gtex ranking", outline=F)
axis(side = 1, at = seq(1.5,ncol(aa_allmat[,aas]), 2), 
     labels = seq(100,5000, 100)[seq(2,100, 6)/2],#c(seq(10,300,10), seq(400,2000, 100), seq(3000,33000,1000)),
     las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = "grey", lwd = 0.7, lty = 3)
legend(x = "topleft",legend = c("location-based", "model-based"), fill = c(2,3), cex = 0.6)

######################################################## pancan
aalab_orig <- ER_SNPs_merged_allbypos_modified$is_pancan
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
aaTP <- aalab_orig[aasap]

aalab_orig <-  ER_SNPs_merged_allbypos_modified$is_pancan
aa_som_rem <- which(ER_SNPs_merged_allbypos_modified$is_somatic)
aa_cli_rem <- which(ER_SNPs_merged_allbypos_modified$is_clinical)

aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[! aasap %in% union(aa_som_rem,aa_cli_rem),]
aasap <- aasap[! aasap %in% union(aa_som_rem,aa_cli_rem)]
aaTP <- aalab_orig[aasap]

aa_pancan_ornot <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
                                             rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
                                             TP = aaTP, 
                                             sampling_range=seq(10,10000, 100),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)),
                                             do_boxplot=T)


aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_pancan_ornot$Random_scores) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, ncol(aa_allmat), 2)){
  aa_allmat[, i] <- aa_pancan_ornot$Random_scores[, aac]
  aa_allmat[, i+1] <- t(aa_pancan_ornot$Model_scores)[, aac]
  aac <- aac + 1
}
par(mfrow = c(1,1), mar = c(6,6,4,4))
boxplot.matrix(aa_allmat[,1:100], col = rep(c(2,3), ncol(aa_pancan_ornot$Random_scores)),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "pancan qtl")
axis(side = 1,
     at = seq(1.5,ncol(aa_allmat), 2), 
     labels = seq(10,10000, 100),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)),
     las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)
legend(x = "topleft",legend = c("location-based", "model-based"), fill = c(2,3), cex = 0.6)

######################################################## gtex or pancan
aalab_orig <- ER_SNPs_merged_allbypos_modified$is_pancan | ER_SNPs_merged_allbypos_modified$is_gtex
aa_som_rem <- which(ER_SNPs_merged_allbypos_modified$is_somatic)
aa_cli_rem <- which(ER_SNPs_merged_allbypos_modified$is_clinical)


aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[! aasap %in% union(aa_som_rem,aa_cli_rem),]
aasap <- aasap[! aasap %in% union(aa_som_rem,aa_cli_rem)]
aaTP <- aalab_orig[aasap]
aa_pancan_gtex_ornot <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
                                             rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
                                             TP = aaTP, 
                                             sampling_range=seq(10,10000, 100), #c(seq(10,400,10), seq(500,3000, 100), seq(3000,31000,1000)),
                                             do_boxplot=T)


aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_pancan_gtex_ornot$Random_scores) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, ncol(aa_allmat), 2)){
  aa_allmat[, i] <- aa_pancan_gtex_ornot$Random_scores[, aac]
  aa_allmat[, i+1] <- t(aa_pancan_gtex_ornot$Model_scores)[, aac]
  aac <- aac + 1
}
par(mfrow = c(1,1), mar = c(6,6,4,4))
boxplot.matrix(aa_allmat[,1:100], col = rep(c(2,3), ncol(aa_pancan_gtex_ornot$Random_scores)),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "gtex/pancan eqtl")
axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
     labels = seq(10,10000, 100),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)), 
     las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)
legend(x = "topleft",legend = c("location-based", "model-based"), fill = c(2,3), cex = 0.6)

######################################################## somatic
aalab_orig <- ER_SNPs_merged_allbypos_modified$is_somatic
aa_som_keep <- which(ER_SNPs_merged_allbypos_modified$is_somatic)
#aa_cli_rem <- which(ER_SNPs_merged_allbypos_modified$is_clinical)
aa_common_keep <-  which(ER_SNPs_merged_allbypos_modified$is_common)


aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[ aasap %in% union(aa_som_keep,aa_common_keep),]
aasap <- aasap[aasap %in% union(aa_som_keep,aa_common_keep)]
aaTP <- aalab_orig[aasap]
aa_somatic_ornot <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
                                                  rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
                                                  TP = aaTP, 
                                                  sampling_range=seq(100,5000, 100), #c(seq(10,400,10), seq(500,3000, 100), seq(3000,31000,1000)),
                                                  do_boxplot=T)


aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_somatic_ornot$Random_scores) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, ncol(aa_allmat), 2)){
  aa_allmat[, i] <- aa_somatic_ornot$Random_scores[, aac]
  aa_allmat[, i+1] <- t(aa_somatic_ornot$Model_scores)[, aac]
  aac <- aac + 1
}

aas <- sort(c(seq(1,100, 6),seq(2,100, 6)))

boxplot.matrix(aa_allmat[,aas], col = rep(c(2,3), length(aas)/2),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "somatic ranking")
axis(side = 1, at = seq(1.5,ncol(aa_allmat[,aas]), 2), 
     labels = seq(100,5000, 100)[seq(2,100, 6)/2],#c(seq(10,300,10), seq(400,2000, 100), seq(3000,33000,1000)),
     las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = "grey", lwd = 0.7, lty = 3)



par(mfrow = c(1,1), mar = c(6,6,4,4))
boxplot.matrix(aa_allmat[,1:100], col = rep(c(2,3), ncol(aa_somatic_ornot$Random_scores)),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "breast cancer somatic from cosmic")
axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
     labels = seq(10,10000, 100),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)), 
     las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)
legend(x = "topleft",legend = c("location-based", "model-based"), fill = c(2,3), cex = 0.6)

######################################################## Somatic_only
aalab_orig <- ER_SNPs_merged_allbypos_modified$is_somatic_only
aalab_orig[is.na(aalab_orig)] <- F
aa_som_keep <- which(ER_SNPs_merged_allbypos_modified$is_somatic)
#aa_cli_rem <- which(ER_SNPs_merged_allbypos_modified$is_clinical)
aa_common_keep <-  which(ER_SNPs_merged_allbypos_modified$is_common)


aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[ aasap %in% union(aa_som_keep,aa_common_keep),]
aasap <- aasap[aasap %in% union(aa_som_keep,aa_common_keep)]
aaTP <- aalab_orig[aasap]
aa_somaticOnly_ornot <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
                                              rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
                                              TP = aaTP, 
                                              sampling_range=seq(10,10000, 100), #c(seq(10,400,10), seq(500,3000, 100), seq(3000,31000,1000)),
                                              do_boxplot=T)


aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_somaticOnly_ornot$Random_scores) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, ncol(aa_allmat), 2)){
  aa_allmat[, i] <- aa_somaticOnly_ornot$Random_scores[, aac]
  aa_allmat[, i+1] <- t(aa_somaticOnly_ornot$Model_scores)[, aac]
  aac <- aac + 1
}
par(mfrow = c(1,1), mar = c(6,6,4,4))
boxplot.matrix(aa_allmat[,1:100], col = rep(c(2,3), ncol(aa_somaticOnly_ornot$Random_scores)),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "breast cancer somatic from cosmic")
axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
     labels = seq(10,10000, 100),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)), 
     las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)
legend(x = "topleft",legend = c("location-based", "model-based"), fill = c(2,3), cex = 0.6)

######################################################## high FATHMM
aalab_orig <- ER_SNPs_merged_allbypos_modified$is_somatic_gt_FATHMM0_7
#aalab_orig[is.na(aalab_orig)] <- F
aa_som_keep <- which(ER_SNPs_merged_allbypos_modified$is_fathmm_scored)
#aa_cli_rem <- which(ER_SNPs_merged_allbypos_modified$is_clinical)
#aa_common_keep <-  which(ER_SNPs_merged_allbypos_modified$is_common)


aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[ aasap %in% aa_som_keep,]
aasap <- aasap[aasap %in% aa_som_keep]
aaTP <- aalab_orig[aasap]
aa_somaticFathmm_ornot <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
                                                  rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
                                                  TP = aaTP, 
                                                  sampling_range=seq(10,1315, 40), #c(seq(10,400,10), seq(500,3000, 100), seq(3000,31000,1000)),
                                                  do_boxplot=T)


aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_somaticFathmm_ornot$Random_scores) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, ncol(aa_allmat), 2)){
  aa_allmat[, i] <- aa_somaticFathmm_ornot$Random_scores[, aac]
  aa_allmat[, i+1] <- t(aa_somaticFathmm_ornot$Model_scores)[, aac]
  aac <- aac + 1
}
par(mfrow = c(1,1), mar = c(6,6,4,4))
boxplot.matrix(aa_allmat[,1:20], col = rep(c(2,3), ncol(aa_somaticFathmm_ornot$Random_scores)),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "FATHMM significant among somatic")
axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
     labels = seq(10,1315, 40),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)), 
     las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)
legend(x = "topleft",legend = c("location-based", "model-based"), fill = c(2,3), cex = 0.6)

######################################################## Clinical pathogenic among clinical + common

aalab_orig <- ER_SNPs_merged_allbypos_modified$is_clinical_pathogenic
#aa_som_keep <- which(ER_SNPs_merged_allbypos_modified$is_somatic)
aa_cli_rem <- which(ER_SNPs_merged_allbypos_modified$is_clinical)
aa_common_keep <-  which(ER_SNPs_merged_allbypos_modified$is_common)


aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[ aasap %in% union(aa_cli_rem,aa_common_keep),]
aasap <- aasap[aasap %in% union(aa_cli_rem,aa_common_keep)]
aaTP <- aalab_orig[aasap]
aaTP[is.na(aaTP)] <- F
aa_clinical_patho_amongAll_ornot <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
                                              rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
                                              TP = aaTP, 
                                              sampling_range=seq(10,5000, 100), #c(seq(10,400,10), seq(500,3000, 100), seq(3000,31000,1000)),
                                              do_boxplot=T)


aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_clinical_patho_amongAll_ornot$Random_scores) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, ncol(aa_allmat), 2)){
  aa_allmat[, i] <- aa_clinical_patho_amongAll_ornot$Random_scores[, aac]
  aa_allmat[, i+1] <- t(aa_clinical_patho_amongAll_ornot$Model_scores)[, aac]
  aac <- aac + 1
}
par(mfrow = c(1,1), mar = c(6,6,4,4))


aas <- sort(c(seq(1,100, 6),seq(2,100, 6)))

boxplot.matrix(aa_allmat[,aas], col = rep(c(2,3), length(aas)/2),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "clinvar ranking")
axis(side = 1, at = seq(1.5,ncol(aa_allmat[,aas]), 2), 
     labels = seq(10,5000, 100)[seq(2,100, 6)/2],#c(seq(10,300,10), seq(400,2000, 100), seq(3000,33000,1000)),
     las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = "grey", lwd = 0.7, lty = 3)


boxplot.matrix(aa_allmat[,1:50], col = rep(c(2,3), ncol(aa_clinical_patho_amongAll_ornot$Random_scores)),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "clinvar patho among all")
axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
     labels = seq(10,5000, 100),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)), 
     las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)
legend(x = "topleft",legend = c("location-based", "model-based"), fill = c(2,3), cex = 0.6)

sum(ER_SNPs_merged_allbypos_modified$is_gtex)/sum(ER_SNPs_merged_allbypos_modified$is_common | ER_SNPs_merged_allbypos_modified$is_gtex)

######################################################## Clinical pathogenic among clinical 

aalab_orig <- ER_SNPs_merged_allbypos_modified$is_clinical_pathogenic
#aa_som_keep <- which(ER_SNPs_merged_allbypos_modified$is_somatic)
aa_cli_rem <- which(ER_SNPs_merged_allbypos_modified$is_clinical)
#aa_common_keep <-  which(ER_SNPs_merged_allbypos_modified$is_common)


aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[ aasap %in% aa_cli_rem,]
aasap <- aasap[aasap %in% aa_cli_rem]
aaTP <- aalab_orig[aasap]
#aaTP[is.na(aaTP)] <- F
aa_clinical_patho_amongClin_ornot <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
                                                              rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
                                                              TP = aaTP, 
                                                              sampling_range=seq(10,1420, 20), #c(seq(10,400,10), seq(500,3000, 100), seq(3000,31000,1000)),
                                                              do_boxplot=T)


aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_clinical_patho_amongClin_ornot$Random_scores) * 2)
#colnames(aa_allmat) <- aasample_size
aac <- 1
for(i in seq(1, ncol(aa_allmat), 2)){
  aa_allmat[, i] <- aa_clinical_patho_amongClin_ornot$Random_scores[, aac]
  aa_allmat[, i+1] <- t(aa_clinical_patho_amongClin_ornot$Model_scores)[, aac]
  aac <- aac + 1
}
par(mfrow = c(1,1), mar = c(6,6,4,4))
boxplot.matrix(aa_allmat[,1:20], col = rep(c(2,3), ncol(aa_clinical_patho_amongClin_ornot$Random_scores)),
               xaxt = "n", xlab = "#predictions", ylab = "TPR", 
               main = "clinvar patho among clinical")
axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
     labels = seq(10,1420, 20),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)), 
     las = 2)
abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)
################################################################################################################ 
# find out what enhancers do the snps fall under
#Merged_percentile_change_noDUP
#ER_SNPs_merged_allbypos_modified

library(BSgenome.Hsapiens.UCSC.hg38)

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)




aaoverlap_holder <- findOverlaps(query = ER_SNPs_merged_allbypos_GR, subject = aa_pos_neg)
aa_snpind <- unique(aaoverlap_holder@from)
aasll <- cbind(aaoverlap_holder@from, aaoverlap_holder@to)

aasll[, 2] <- names(aa_pos_neg)[aasll[, 2]]
aasll <- as.data.frame(aasll,stringsAsFactors =F)
aasll_agg <- aggregate(aasll["V2"], by = aasll["V1"], FUN = c)
aasll_agg$V1 <- as.numeric(aasll_agg$V1)
aaset <- sort(aasll_agg$V1, index.return = T)$ix
aasll_agg <- aasll_agg[aaset,]
rownames(aasll_agg) <- as.character(aasll_agg$V1)
aaln <- unlist(lapply(aasll_agg$V2, length))


ER_SNPs_merged_allbypos_modified$overlapping_enh <- aasll_agg$V2
############################################################### 
## all  overlap or interaction between enhancers, genes and their 20kb vicinity:
# enhancer_gene_ovlap_interact
# enhancer_enhancer_interact

Enhancer_Enhancer_interac_mat
Enhancer_Enhancer_overlap_mat

Enhancer_genePromoter_interac_mat
Enhancer_genePromoter_overlap_mat

Enhancer_gene_20kb_vicinity_mat
###############################################################
# add genes in vicinity and interacting genes to the ER_SNPs_merged_allbypos_modified dataframe
ER_SNPs_merged_allbypos_modified$overlapping_enh[1:5]
aa_snp_enh_gene_ov <- list()
aa_snp_enh_gene_in <- list()
aa_snp_enh_enha_ov <- list()
aa_snp_enh_enha_in <- list()

for(i in 1:nrow(ER_SNPs_merged_allbypos_modified)){
  aa_snp_enh_gene_ov[[i]] <- unique(unlist(Enhancer_gene_20kb_vicinity_mat$gene_entID[Enhancer_gene_20kb_vicinity_mat$enh_name %in% ER_SNPs_merged_allbypos_modified$overlapping_enh[[i]]]))
  aa_snp_enh_gene_in[[i]] <- unique(unlist(Enhancer_genePromoter_interac_mat_experimental_only$gene_entID[Enhancer_genePromoter_interac_mat$enh_name %in% ER_SNPs_merged_allbypos_modified$overlapping_enh[[i]]]))
  aa_snp_enh_enha_ov[[i]] <- unique(unlist(Enhancer_Enhancer_overlap_mat$enh_name_2[Enhancer_Enhancer_overlap_mat$enh_name_1 %in% ER_SNPs_merged_allbypos_modified$overlapping_enh[[i]]]))
  aa_snp_enh_enha_in[[i]] <- unique(unlist(Enhancer_Enhancer_interac_mat$enh_name_2[Enhancer_Enhancer_interac_mat$enh_name_1 %in% ER_SNPs_merged_allbypos_modified$overlapping_enh[[i]]]))
}

ER_SNPs_merged_allbypos_modified$Enh_gene20kb_ovl <- aa_snp_enh_gene_ov
ER_SNPs_merged_allbypos_modified$Enh_geneProm_int <- aa_snp_enh_gene_in
ER_SNPs_merged_allbypos_modified$Enh_Enh_int <- aa_snp_enh_enha_in
ER_SNPs_merged_allbypos_modified$Enh_Enh_ovl <- aa_snp_enh_enha_ov



##############################################################################################################################
# rank genes based on SNP scores of all vaiants

aa_ovlgenes <- unique(unlist(Enhancer_gene_20kb_vicinity_mat$gene_entID))
aa_intgenes <- unique(unlist(Enhancer_genePromoter_interac_mat$gene_entID))
aa_all_genes <- union(aa_ovlgenes,aa_intgenes )


## 
# read drug genes
aa_files_full <- list.files("~/Documents/Shayan/BioInf/EstrogenReceptor/drug_response/", full.names = T)
aa_files <- list.files("~/Documents/Shayan/BioInf/EstrogenReceptor/drug_response/", full.names = F)
aa_drug_response <- list()
aa_drug_response_20 <- list()
for(i in 1:length(aa_files_full)){
  print(i)
  aa_drug_response[[i]] <- read.delim(file = aa_files_full[i], stringsAsFactors = F)
  aa_drug_response_20[[i]] <- aa_drug_response[[i]][aa_drug_response[[i]]$PVAL <= 0.2,]
  aa_drug_response[[i]] <- aa_drug_response[[i]][aa_drug_response[[i]]$PVAL <= 0.05,]
  
  aa_drug_response[[i]]$ENSEMBL_new <- as.character(as.numeric(unlist(lapply(strsplit(aa_drug_response[[i]]$ENSEMBL, "ENSG"), "[[", 2))))
  aa_drug_response_20[[i]]$ENSEMBL_new <- as.character(as.numeric(unlist(lapply(strsplit(aa_drug_response_20[[i]]$ENSEMBL, "ENSG"), "[[", 2))))
}
names(aa_drug_response) <- aa_files

aa_drug_response_related <- list()
aa_drug_response_related_20 <- list()
for(i in 1:length(aa_drug_response)){
  aa_drug_response_related[[i]] <- aa_drug_response[[i]][aa_drug_response[[i]]$ENSEMBL_new %in% aa_all_genes, ]
  aa_drug_response_related_20[[i]] <- aa_drug_response_20[[i]][aa_drug_response_20[[i]]$ENSEMBL_new %in% aa_all_genes, ]
}
names(aa_drug_response_related) <- names(aa_drug_response)
names(aa_drug_response_related_20) <- names(aa_drug_response)

unlist(lapply(aa_drug_response_related, nrow))
unlist(lapply(aa_drug_response_related_20, nrow))


# ER_SNPs_merged_allbypos_modified$Enh_gene20kb_ovl
# ER_SNPs_merged_allbypos_modified$Enh_geneProm_int
aaasst <- mapply(union, ER_SNPs_merged_allbypos_modified$Enh_gene20kb_ovl, ER_SNPs_merged_allbypos_modified$Enh_geneProm_int)
ER_SNPs_merged_allbypos_modified$Enh_gene_both <- aaasst



aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
#aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[ aasap %in% aa_cli_rem,]

aa_max <- apply(abs(Merged_percentile_change_noDUP), MARGIN = 1, FUN = max)
aa_mea <- apply(abs(Merged_percentile_change_noDUP), MARGIN = 1, FUN = mean)
aa_med <- apply(abs(Merged_percentile_change_noDUP), MARGIN = 1, FUN = median)

aa_model_sort_ind_max <- aasap[sort(aa_max, decreasing = T, index.return = T)$ix]
aa_model_sort_ind_mea <- aasap[sort(aa_mea, decreasing = T, index.return = T)$ix]
aa_model_sort_ind_med <- aasap[sort(aa_med, decreasing = T, index.return = T)$ix]

aa_percModels <- rowSums(abs(Merged_percentile_change_noDUP) > 0.001)/ncol(Merged_percentile_change_noDUP)
aa_perc_names <- aasap[sort(aa_percModels,decreasing = T,index.return=T)$ix]


aa_gene_ranked_max <- unlist(ER_SNPs_merged_allbypos_modified$Enh_gene_both[aa_model_sort_ind_max])
aa_gene_ranked_mea <- unlist(ER_SNPs_merged_allbypos_modified$Enh_gene_both[aa_model_sort_ind_mea])
aa_gene_ranked_med <- unlist(ER_SNPs_merged_allbypos_modified$Enh_gene_both[aa_model_sort_ind_med])
aa_gene_ranked_per <- unlist(ER_SNPs_merged_allbypos_modified$Enh_gene_both[aa_perc_names])

aa_gene_ranked_max <- aa_gene_ranked_max[!duplicated(aa_gene_ranked_max)]
aa_gene_ranked_mea <- aa_gene_ranked_mea[!duplicated(aa_gene_ranked_mea)]
aa_gene_ranked_med <- aa_gene_ranked_med[!duplicated(aa_gene_ranked_med)]
aa_gene_ranked_per <- aa_gene_ranked_per[!duplicated(aa_gene_ranked_per)]


aa_drug_gene_percentile_mean <- list()
aa_drug_gene_percentile_max <- list()
aa_drug_gene_percentile_med <- list()
aa_drug_gene_percentile_vote <- list()


for(i in 1:length(aa_drug_response_related_20)){
  aa_drug_gene_percentile_mean[[i]] <- numeric(length = nrow(aa_drug_response_related_20[[i]]))
  aa_drug_gene_percentile_max[[i]] <- numeric(length = nrow(aa_drug_response_related_20[[i]]))
  aa_drug_gene_percentile_med[[i]] <- numeric(length = nrow(aa_drug_response_related_20[[i]]))
  aa_drug_gene_percentile_vote[[i]] <- numeric(length = nrow(aa_drug_response_related_20[[i]]))
  if(nrow(aa_drug_response_related_20[[i]]) > 0){
    for(j in 1:nrow(aa_drug_response_related_20[[i]])){
      aax <- which(aa_gene_ranked_mea == aa_drug_response_related_20[[i]]$ENSEMBL_new[j])
      if(length(aax) == 0){
        aax <- length(aa_gene_ranked_mea)
      }
      aa_drug_gene_percentile_mean[[i]][j] <- 1 - aax/length(aa_gene_ranked_mea)
      
      aax <- which(aa_gene_ranked_max == aa_drug_response_related_20[[i]]$ENSEMBL_new[j])
      if(length(aax) == 0){
        aax <- length(aa_gene_ranked_max)
      }
      aa_drug_gene_percentile_max[[i]][j] <- 1 - aax/length(aa_gene_ranked_max)
      
      aax <- which(aa_gene_ranked_med == aa_drug_response_related_20[[i]]$ENSEMBL_new[j])
      if(length(aax) == 0){
        aax <- length(aa_gene_ranked_med)
      }
      aa_drug_gene_percentile_med[[i]][j] <- 1 - aax/length(aa_gene_ranked_med)
      
      aax <- which(aa_gene_ranked_per == aa_drug_response_related_20[[i]]$ENSEMBL_new[j])
      if(length(aax) == 0){
        aax <- length(aa_gene_ranked_per)
      }
      aa_drug_gene_percentile_vote[[i]][j] <- 1 - aax/length(aa_gene_ranked_per)
      }
  }
}

#names(aa_drug_gene_percentile_vote) <- names(aa_drug_response_related)
aa_mat_max <- matrix(nrow = max(unlist(lapply(aa_drug_response_related_20, nrow))),
                     ncol = length(aa_drug_response_related_20))
aa_mat_mea <- matrix(nrow = max(unlist(lapply(aa_drug_response_related_20, nrow))),
                     ncol = length(aa_drug_response_related_20))
aa_mat_med <- matrix(nrow = max(unlist(lapply(aa_drug_response_related_20, nrow))),
                     ncol = length(aa_drug_response_related_20))
aa_mat_vot <- matrix(nrow = max(unlist(lapply(aa_drug_response_related_20, nrow))),
                     ncol = length(aa_drug_response_related_20))

colnames(aa_mat_max) <- paste0(names(aa_drug_response_related_20), "__" ,unlist(lapply(aa_drug_response_related_20, nrow)))
colnames(aa_mat_mea) <- paste0(names(aa_drug_response_related_20), "__" ,unlist(lapply(aa_drug_response_related_20, nrow)))
colnames(aa_mat_med) <- paste0(names(aa_drug_response_related_20), "__" ,unlist(lapply(aa_drug_response_related_20, nrow)))
colnames(aa_mat_vot) <- paste0(names(aa_drug_response_related_20), "__" ,unlist(lapply(aa_drug_response_related_20, nrow)))

for(i in 1:ncol(aa_mat_max)){
  if(nrow(aa_drug_response_related_20[[i]]) > 0){
    aa_mat_max[1:nrow(aa_drug_response_related_20[[i]]), i] <- aa_drug_gene_percentile_max[[i]]
    aa_mat_mea[1:nrow(aa_drug_response_related_20[[i]]), i] <- aa_drug_gene_percentile_mean[[i]]
    aa_mat_med[1:nrow(aa_drug_response_related_20[[i]]), i] <- aa_drug_gene_percentile_med[[i]]
    aa_mat_vot[1:nrow(aa_drug_response_related_20[[i]]), i] <- aa_drug_gene_percentile_vote[[i]]
  }
  
}

#names(aa_drug_snp)
aa_breast_use <- c("not_used", "not_used", "not_used", "not_used", "used", "used", "not_used", "used", "used", "used", "not_used", "used", "used", "used", "used", "used", "not_used", "not_used", "used", "used", "used", "not_used", "used")
aa_cancer_use <- c("cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "non_cancer", "cancer", "cancer", "non_cancer", "cancer", "cancer", "cancer", "cancer", "non_cancer", "cancer")
names(aa_breast_use) <- names(aa_drug_snp)
names(aa_cancer_use) <- names(aa_drug_snp)

aa_col <- aa_breast_use
aa_col[aa_col == "not_used"] <- 2
aa_col[aa_col == "used"] <- 3

par(mfrow = c(1,1), mar = c(10,6,4,4))
boxplot.matrix(aa_mat_max, las = 2,
               ylab = "percentile among ranked genes",
               main = "max_over_ensemble", col = aa_col)
boxplot.matrix(aa_mat_mea, las = 2, 
               ylab = "percentile among ranked genes", 
               main = "mean_over_ensemble", col = aa_col)
boxplot.matrix(aa_mat_med, las = 2, 
               ylab = "percentile among ranked genes",
               main = "median_over_ensemble", col = aa_col)
boxplot.matrix(aa_mat_vot, las = 2, 
               ylab = "percentile among ranked genes",
               main = "vote_over_ensemble"
               #, col = aa_col
               )
##############################################################################################################################


# investigate all variants with some predicted change


#ER_SNPs_merged_allbypos_modified
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
Merged_percentile_change_noDUP_withInd <- cbind(aasap, Merged_percentile_change_noDUP)
rownames(Merged_percentile_change_noDUP_withInd) <- NULL
colnames(Merged_percentile_change_noDUP_withInd)[1] <- "pos_uniq_id"
Merged_percentile_change_noDUP_withInd <- as.data.frame(Merged_percentile_change_noDUP_withInd)
aars <- rowSums(Merged_percentile_change_noDUP_withInd == 0)
Merged_percentile_change_noDUP_withInd_someImpact <- Merged_percentile_change_noDUP_withInd[aars != 244,]

aatmpmat <- as.matrix(Merged_percentile_change_noDUP_withInd[, 2:ncol(Merged_percentile_change_noDUP_withInd)])
aa_max <- apply(abs(aatmpmat), MARGIN = 1, FUN = max)
for(i in 1:nrow(aatmpmat)){
  aa_max[i] <- aa_max[i] * sign(aatmpmat[i, which.max(abs(aatmpmat[i,]))])
}
aa_mea <- apply(aatmpmat, MARGIN = 1, FUN = mean)
aa_med <- apply(aatmpmat, MARGIN = 1, FUN = median)
aa_perc <- rowSums(abs(aatmpmat) > 0.001)/ncol(aatmpmat)
Merged_percentile_change_noDUP_withInd$median_change <- aa_med
Merged_percentile_change_noDUP_withInd$mean_change <- aa_mea
Merged_percentile_change_noDUP_withInd$max_change <- aa_max
Merged_percentile_change_noDUP_withInd$percent_changed <- aa_perc






ER_SNPs_merged_allbypos_modified_and_modelScores <- full_join(x = ER_SNPs_merged_allbypos_modified,
                                                              y = Merged_percentile_change_noDUP_withInd,
                                                              by = c("pos_uniq_id"))

#Index of variants to be evaluated: Merged_percentile_change_noDUP_withInd_someImpact$pos_uniq_id
# save necessary things and run this in parallel
library(BSgenome.Hsapiens.UCSC.hg38)
ER_SNPs_merged_allbypos_modified_and_modelScores_filtered <- ER_SNPs_merged_allbypos_modified_and_modelScores[ER_SNPs_merged_allbypos_modified_and_modelScores$pos_uniq_id %in% Merged_percentile_change_noDUP_withInd_someImpact$pos_uniq_id,]
save(list = c("snp_investigator", "MotifScore", "CompLLRRandom",
              "ER_SNPs_merged_allbypos_modified_and_modelScores_filtered",
              "TF.motifs.Expanded_pseudo_exp12_t", "aa_exp18_minLLR"),
     file = "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/snp_investigator_workplace.RData")


aa_ref <- ER_SNPs_merged_allbypos_modified_and_modelScores_filtered$ref
aa_alt <- ER_SNPs_merged_allbypos_modified_and_modelScores_filtered$alt
# unique(unlist(strsplit(aa_ref, "")))
# unique(unlist(strsplit(aa_alt, "")))

aa_ref_Alt <- cbind(aa_ref, aa_alt)


aa_eqtl_chr <- ER_SNPs_merged_allbypos_modified_and_modelScores_filtered$seqnames
aa_eqtl_st <- ER_SNPs_merged_allbypos_modified_and_modelScores_filtered$start
aa_eqtl_ref <- aa_ref_Alt[, 1]
aa_eqtl_alt <- aa_ref_Alt[, 2]
aa_names <- c("ESR1_2","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","RARA", "RELA","RUNX1","SP1","TFAP2C","YBX1")

aamin_LLR_low <- aa_exp18_minLLR - 0.1

ER_all_snp_investigate <- list()
for(i in 1:length(aa_eqtl_alt)){
  print(i)
  ER_all_snp_investigate[[i]] <- snp_investigator(my_genome = BSgenome.Hsapiens.UCSC.hg38,
                                                    my_motifs = TF.motifs.Expanded_pseudo_exp12_t,
                                                    snp_chr = aa_eqtl_chr[i],
                                                    snp_start = aa_eqtl_st[i],
                                                    snp_ref = aa_eqtl_ref[i],
                                                    snp_alt = aa_eqtl_alt[i],
                                                    min_LLR = aamin_LLR_low[aa_names])
  
}
names(ER_all_snp_investigate) <- ER_SNPs_merged_allbypos_modified_and_modelScores_filtered$pos_uniq_id
save(list = c("ER_all_snp_investigate"), 
     file = "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ER_all_variants_filtered_investigated.RData")
aarn <- unlist(lapply(ER_all_snp_investigate, nrow))
table(aarn)
load("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ER_all_variants_filtered_investigated.RData")


ER_all_snp_investigate_mat <- do.call(rbind, ER_all_snp_investigate)
ER_all_snp_investigate_mat <- ER_all_snp_investigate_mat[ER_all_snp_investigate_mat[, 4] != 0,]

aaxcv <- as.numeric(unlist(lapply(strsplit(rownames(ER_all_snp_investigate_mat), "\\."), "[[", 1)))
table(as.numeric(table(aaxcv)))
ER_all_snp_investigate_mat$pos_uniq_id <- aaxcv
aamaxL_LLR

aachange_max <- numeric(length = nrow(ER_all_snp_investigate_mat))
for(i in 1:length(aachange_max)){
  aachange_max[i] <- (ER_all_snp_investigate_mat$new[i] - ER_all_snp_investigate_mat$old[i])/aamaxL_LLR[ER_all_snp_investigate_mat$TF[i]]
}
ER_all_snp_investigate_mat$change_percent_max <- aachange_max
ER_all_snp_investigate_mat_withcoord <- left_join(ER_all_snp_investigate_mat,
                                                  ER_SNPs_merged_allbypos_modified[, c(c(1:3), 23)], 
                                                  by = c("pos_uniq_id"))
Change_in_binding_snp_byTF <- list()
aa_names <- c("ESR1_2","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","RARA", "RELA","RUNX1","SP1","TFAP2C","YBX1")

for(i in 1:length(aa_names)){
  Change_in_binding_snp_byTF[[i]] <- ER_all_snp_investigate_mat_withcoord[ER_all_snp_investigate_mat_withcoord$TF %in% aa_names[i],]
}
names(Change_in_binding_snp_byTF) <- aa_names


###################################################################### aggregate by position only
############################################################################################################################################
# checking how many snps have the same corrdinates but different ref/alts

aamg <- c("ref", "alt", "pos_uniq_id")
ER_SNPs_merged_allbypos_modified_byPOSonly <- aggregate(x = ER_SNPs_merged_allbypos_modified[aamg],
                                                        by = ER_SNPs_merged_allbypos_modified[c("seqnames", "start", "end")], 
                                                        FUN = c)


aatb <- unlist(lapply(ER_SNPs_merged_allbypos_modified_byPOSonly$pos_uniq_id , length))

ER_SNPs_merged_allbypos_modified_byPOSonly[aatb > 5,]

# in order to create expression prediction change plots for GViz: I'll create two types of data
# first: summaries of model predictions (mean, median, max, percent) can be added to the main dataframe: ER_SNPs_merged_allbypos_modified
ER_SNPs_merged_allbypos_modified_and_modelScores$median_change
ER_SNPs_merged_allbypos_modified_and_modelScores$mean_change
ER_SNPs_merged_allbypos_modified_and_modelScores$max_change
ER_SNPs_merged_allbypos_modified_and_modelScores$percent_changed



# second for full model prediction visualization (boxplot) we need to choose one alternative variant per position:
# for the ER_SNPs_merged_allbypos_modified_byPOSonly dataframe we will go through each row and choose one pos_uniq_id with the highest median change among models

#Merged_percentile_change_noDUP_withInd
aaln <- unlist(lapply(ER_SNPs_merged_allbypos_modified_byPOSonly$pos_uniq_id, length))
aa_visual_expr_change_mat_posUniq <- unlist(ER_SNPs_merged_allbypos_modified_byPOSonly$pos_uniq_id[aaln == 1])
for(i in 1:nrow(ER_SNPs_merged_allbypos_modified_byPOSonly)){
  if(aaln[i] > 1){
    aacurpos <- ER_SNPs_merged_allbypos_modified_byPOSonly$pos_uniq_id[[i]]
    aacur_ind <- match(aacurpos, ER_SNPs_merged_allbypos_modified_and_modelScores$pos_uniq_id)
    aamx <- which.max(abs(ER_SNPs_merged_allbypos_modified_and_modelScores$median_change[aacur_ind]))
    aa_visual_expr_change_mat_posUniq <- c(aa_visual_expr_change_mat_posUniq, aacurpos[aamx])
  }
}

sum(duplicated(aa_visual_expr_change_mat_posUniq))
aa_visual_expr_change_mat_posUniq <- sort(aa_visual_expr_change_mat_posUniq)
One_alternative_index <- aa_visual_expr_change_mat_posUniq
################################################################################
# add the overlapping TFs as an extra column to ER_SNPs_merged_allbypos_modified
#ER_SNPs_merged_allbypos_modified
#Enhancer_TF_overlap_name_expOnly

aa_all_ovl_TF <- list()
for(i in 1:nrow(ER_SNPs_merged_allbypos_modified)){
  aacneh <- unlist(ER_SNPs_merged_allbypos_modified$overlapping_enh[[i]])
  names(aacneh) <- NULL
  aa_cur <- unique(unlist(Enhancer_TF_overlap_name_expOnly[which(names(Enhancer_TF_overlap_name_expOnly) %in% aacneh)]))
  aa_all_ovl_TF[[i]] <- aa_cur
}

ER_SNPs_merged_allbypos_modified$overlapping_TFs <- aa_all_ovl_TF
rm(aa_all_ovl_TF)

################################################################################
# add a column for responsible TF to ER_SNPs_merged_allbypos_modified

# first aggregate by pos_uniq_id
aamg <- c("TF", "change", "change_percent_max")
ER_all_snp_investigate_mat$TF[ER_all_snp_investigate_mat$TF == "ESR1_2"] <- "ESR1"
ER_all_snp_investigate_mat_by_uniq_id <- aggregate(x = ER_all_snp_investigate_mat[aamg],
                                                        by = ER_all_snp_investigate_mat[c("pos_uniq_id")], 
                                                        FUN = c)

names(ER_all_snp_investigate_mat_by_uniq_id)[c(2,3,4)] <- c("Responsible_TF","LLR_percent_change_old", "LLR_percent_change_max")

aaER_SNPs_merged_allbypos_modified_resp <- full_join(x = ER_SNPs_merged_allbypos_modified,
                                                              y = ER_all_snp_investigate_mat_by_uniq_id,
                                                              by = c("pos_uniq_id"))
ER_SNPs_merged_allbypos_modified <- aaER_SNPs_merged_allbypos_modified_resp
rm(aaER_SNPs_merged_allbypos_modified_resp)

aarepsup <- list()
for(i in 1:nrow(ER_SNPs_merged_allbypos_modified)){
  aarepsup[[i]] <- ER_SNPs_merged_allbypos_modified$Responsible_TF[[i]] %in% ER_SNPs_merged_allbypos_modified$overlapping_TFs[[i]]
}
ER_SNPs_merged_allbypos_modified$Responsible_TF_support <- aarepsup


# filter some eqtls to get some which satisfy the following criteria
#is pancan breast cancer eqtl,
# some change predicted, 
# correlating gene is in 20kb or interacts with the enhancer
# responsible TF does have chip support

ER_SNPs_merged_allbypos_modified$is_pancan
ER_SNPs_merged_allbypos_modified$Gene_cor_pancan_entrez
ER_SNPs_merged_allbypos_modified$Enh_gene_both
ER_SNPs_merged_allbypos_modified$Enh_geneProm_int
ER_SNPs_merged_allbypos_modified$Responsible_TF_support


aawpancan <- which(ER_SNPs_merged_allbypos_modified$is_pancan)
aa_is_aroud <- list()
aa_is_int <- list()
aa_supp <- numeric(length = length(aawpancan))
for(i in 1:length(aawpancan)){
  aa_is_aroud[[i]] <- intersect(ER_SNPs_merged_allbypos_modified$Enh_gene20kb_ovl[[aawpancan[i]]],
                                ER_SNPs_merged_allbypos_modified$Gene_cor_pancan_entrez[[aawpancan[i]]])
  aa_is_int[[i]] <- intersect(ER_SNPs_merged_allbypos_modified$Enh_geneProm_int[[aawpancan[i]]],
                              ER_SNPs_merged_allbypos_modified$Gene_cor_pancan_entrez[[aawpancan[i]]])
  aa_supp[i] <- sum(ER_SNPs_merged_allbypos_modified$Responsible_TF_support[[aawpancan[i]]])
  
}

aaa <- unlist(lapply(aa_is_aroud, length))
aacand <- aawpancan[which(aaa > 0 & aa_supp > 0)]
aacand_on <- aawpancan[aa_supp > 0]


View(ER_SNPs_merged_allbypos_modified[aacand,])
View(ER_SNPs_merged_allbypos_modified[aacand_on,])

ER_SNPs_merged_allbypos_modified$overlapping_TFs[[32887]]

aawsom <- which(ER_SNPs_merged_allbypos_modified$is_somatic)
# aa_is_aroud <- list()
# aa_is_int <- list()
aa_supp <- numeric(length = length(aawsom))
for(i in 1:length(aawsom)){
  # aa_is_aroud[[i]] <- intersect(ER_SNPs_merged_allbypos_modified$Enh_gene20kb_ovl[[aawsom[i]]],
  #                               ER_SNPs_merged_allbypos_modified$Gene_cor_pancan_entrez[[aawsom[i]]])
  # aa_is_int[[i]] <- intersect(ER_SNPs_merged_allbypos_modified$Enh_geneProm_int[[aawsom[i]]],
  #                             ER_SNPs_merged_allbypos_modified$Gene_cor_pancan_entrez[[aawsom[i]]])
  aa_supp[i] <- sum(ER_SNPs_merged_allbypos_modified$Responsible_TF_support[[aawsom[i]]])
  
}

aacand2 <- aawsom[which(aa_supp > 0)]

View(ER_SNPs_merged_allbypos_modified[aacand2,])
View(ER_SNPs_merged_allbypos_modified[aacand,])


aacand_3_4 <- aawsom[which(aa_supp >= 3)]
View(ER_SNPs_merged_allbypos_modified[aacand_3_4,])

for(i in 1:length(aacand2)){
  print(i)
  if(abs(ER_SNPs_merged_allbypos_modified$median_change[aacand2[i]]) >= 0.04){
    aa_my_an_list <- list(aitrack, 
                          agtrack, 
                          aaot,
                          aasnp_track,
                          aaot_change,
                          aaexpression_change_median_track,
                          aa_tmp_erna_beforeE2_AnnotationTrack,
                          aa_tmp_erna_only_afterE2_AnnotationTrack,
                          aa_ER_AnnotationTrack,
                          aa_H3k27ac_AnnotationTrack)
    aacnt <- 1
    if(length(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacand2[i]]]) > 0){
      aactf <- unique(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacand2[i]]])
      for(j in 1:length(aactf)){
        if(aactf[j] %in% names(aa_remap_track_list)){
          aa_my_an_list[[length(aa_my_an_list) + 1]] <- aa_remap_track_list[[which(names(aa_remap_track_list) %in% aactf[j])]]
        }
      }
    }
    for(j in 1:length(ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand2[i]]])){
      aa_curenh <- which(names(aa_pos_neg) == ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand2[i]]][j])
      png(paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/Visual/3_",
                 ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand2[i]]][j] ,
                 ".png"),    # create PNG for the heat map        
          width = 8*300,        # 5 x 300 pixels
          height = 11*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8)        # smaller font size
      
      plotTracks(aa_my_an_list, 
                 from = start(aa_pos_neg)[aa_curenh], 
                 to = end(aa_pos_neg[aa_curenh]),
                 chromosome = as.character(seqnames(aa_pos_neg[aa_curenh])))
      dev.off()
    }
  }

}

View(ER_SNPs_merged_allbypos_modified[11642,])
View(ER_SNPs_merged_allbypos_modified[32887,])


aagt_gene <- unique(unlist(ER_SNPs_merged_allbypos_modified$Gene_cor_Gtex[which(ER_SNPs_merged_allbypos_modified$is_gtex)]))
aagt_gene2 <- unlist(lapply(strsplit(aagt_gene, "\\."), "[[", 1))
aagt_gene3 <- as.character(as.numeric(unlist(lapply(strsplit(aagt_gene2, "ENSG"), "[[", 2))))
aagtgene_convert <- cbind(aagt_gene, aagt_gene3)

aagt <- which(ER_SNPs_merged_allbypos_modified$is_gtex)
aagt_newlist <- list()
for(i in 1:nrow(ER_SNPs_merged_allbypos_modified)){
  if(i %in% aagt){
    
    aa_new <- aagtgene_convert[match(unlist(ER_SNPs_merged_allbypos_modified$Gene_cor_Gtex[[i]]), aagtgene_convert[, 1]) , 2]
    names(aa_new) <- NULL
    aagt_newlist[[i]] <- aa_new
  }else{
    aagt_newlist[[i]] <- NA
  }
}
ER_SNPs_merged_allbypos_modified$Gene_cor_Gtex_formatted <- aagt_newlist

aagt <- which(ER_SNPs_merged_allbypos_modified$is_gtex)
aa_is_aroud <- list()
aa_is_int <- list()
aa_supp <- numeric(length = length(aagt))
for(i in 1:length(aagt)){
  aa_is_aroud[[i]] <- intersect(ER_SNPs_merged_allbypos_modified$Enh_gene20kb_ovl[[aagt[i]]],
                                ER_SNPs_merged_allbypos_modified$Gene_cor_Gtex_formatted[[aagt[i]]])
  aa_is_int[[i]] <- intersect(ER_SNPs_merged_allbypos_modified$Enh_geneProm_int[[aagt[i]]],
                              ER_SNPs_merged_allbypos_modified$Gene_cor_Gtex_formatted[[aagt[i]]])
  aa_supp[i] <- sum(ER_SNPs_merged_allbypos_modified$Responsible_TF_support[[aagt[i]]])
  
}

aaa <- unlist(lapply(aa_is_aroud, length))
aaa2 <- unlist(lapply(aa_is_int, length))
aacand3 <- aagt[which(aa_supp > 0)]

View(ER_SNPs_merged_allbypos_modified[aacand3,])


for(i in 1:length(aacand)){
  print(i)
  if(abs(ER_SNPs_merged_allbypos_modified$median_change[aacand[i]]) >= 0.04){
    aa_my_an_list <- list(aitrack, 
                          agtrack, 
                          aaot,
                          aasnp_track,
                          aaot_change,
                          aaexpression_change_median_track,
                          #aa_tmp_erna_beforeE2_AnnotationTrack,
                          aa_tmp_erna_only_afterE2_AnnotationTrack,
                          aa_ER_AnnotationTrack,
                          aa_H3k27ac_AnnotationTrack)
    aacnt <- 1
    if(length(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacand[i]]]) > 0){
      aactf <- unique(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacand[i]]])
      for(j in 1:length(aactf)){
        if(aactf[j] %in% names(aa_remap_track_list)){
          aa_my_an_list[[length(aa_my_an_list) + 1]] <- aa_remap_track_list[[which(names(aa_remap_track_list) %in% aactf[j])]]
        }
      }
    }
    for(j in 1:length(ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand[i]]])){
      aa_curenh <- which(names(aa_pos_neg) == ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand[i]]][j])
      png(paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/Visual/1_",
                 ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand[i]]][j] ,
                 ".png"),    # create PNG for the heat map        
          width = 8*300,        # 5 x 300 pixels
          height = 11*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8)        # smaller font size
      
      plotTracks(aa_my_an_list, 
                 from = start(aa_pos_neg)[aa_curenh], 
                 to = end(aa_pos_neg[aa_curenh]),
                 chromosome = as.character(seqnames(aa_pos_neg[aa_curenh])))
      dev.off()
    }
  }
  
}

for(i in 1:length(aacand_on)){
  print(i)
  if(abs(ER_SNPs_merged_allbypos_modified$median_change[aacand_on[i]]) >= 0.04){
    aa_my_an_list <- list(aitrack, 
                          agtrack, 
                          aaot,
                          aasnp_track,
                          aaot_change,
                          aaexpression_change_median_track,
                          aa_tmp_erna_beforeE2_AnnotationTrack,
                          aa_tmp_erna_only_afterE2_AnnotationTrack,
                          aa_ER_AnnotationTrack,
                          aa_H3k27ac_AnnotationTrack)
    aacnt <- 1
    if(length(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacand_on[i]]]) > 0){
      aactf <- unique(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacand_on[i]]])
      for(j in 1:length(aactf)){
        if(aactf[j] %in% names(aa_remap_track_list)){
          aa_my_an_list[[length(aa_my_an_list) + 1]] <- aa_remap_track_list[[which(names(aa_remap_track_list) %in% aactf[j])]]
        }
      }
    }
    for(j in 1:length(ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand_on[i]]])){
      aa_curenh <- which(names(aa_pos_neg) == ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand_on[i]]][j])
      png(paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/Visual/1_",
                 ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand_on[i]]][j] ,
                 ".png"),    # create PNG for the heat map        
          width = 8*300,        # 5 x 300 pixels
          height = 11*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8)        # smaller font size
      
      plotTracks(aa_my_an_list, 
                 from = start(aa_pos_neg)[aa_curenh], 
                 to = end(aa_pos_neg[aa_curenh]),
                 chromosome = as.character(seqnames(aa_pos_neg[aa_curenh])))
      dev.off()
    }
  }
  
}


for(i in 1:length(aacand3)){
  print(i)
  if(abs(ER_SNPs_merged_allbypos_modified$median_change[aacand3[i]]) >= 0.04){
    aa_my_an_list <- list(aitrack, 
                          agtrack, 
                          aaot,
                          aasnp_track,
                          aaot_change,
                          aaexpression_change_median_track,
                          aa_tmp_erna_beforeE2_AnnotationTrack,
                          aa_tmp_erna_only_afterE2_AnnotationTrack,
                          aa_ER_AnnotationTrack,
                          aa_H3k27ac_AnnotationTrack)
    aacnt <- 1
    if(length(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacand3[i]]]) > 0){
      aactf <- unique(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacand3[i]]])
      for(j in 1:length(aactf)){
        if(aactf[j] %in% names(aa_remap_track_list)){
          aa_my_an_list[[length(aa_my_an_list) + 1]] <- aa_remap_track_list[[which(names(aa_remap_track_list) %in% aactf[j])]]
        }
      }
    }
    for(j in 1:length(ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand3[i]]])){
      aa_curenh <- which(names(aa_pos_neg) == ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand3[i]]][j])
      png(paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/Visual/2_",
                 ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacand3[i]]][j] ,
                 ".png"),    # create PNG for the heat map        
          width = 8*300,        # 5 x 300 pixels
          height = 11*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8)        # smaller font size
      
      plotTracks(aa_my_an_list, 
                 from = start(aa_pos_neg)[aa_curenh], 
                 to = end(aa_pos_neg[aa_curenh]),
                 chromosome = as.character(seqnames(aa_pos_neg[aa_curenh])))
      dev.off()
    }
  }
  
}

aascheme <- getScheme()



aawcli <- which(ER_SNPs_merged_allbypos_modified$is_clinical_pathogenic)
# aa_is_aroud <- list()
# aa_is_int <- list()
aa_supp <- numeric(length = length(aawcli))
for(i in 1:length(aawcli)){
  # aa_is_aroud[[i]] <- intersect(ER_SNPs_merged_allbypos_modified$Enh_gene20kb_ovl[[aawsom[i]]],
  #                               ER_SNPs_merged_allbypos_modified$Gene_cor_pancan_entrez[[aawsom[i]]])
  # aa_is_int[[i]] <- intersect(ER_SNPs_merged_allbypos_modified$Enh_geneProm_int[[aawsom[i]]],
  #                             ER_SNPs_merged_allbypos_modified$Gene_cor_pancan_entrez[[aawsom[i]]])
  aa_supp[i] <- sum(ER_SNPs_merged_allbypos_modified$Responsible_TF_support[[aawcli[i]]])
  
}
aacondcli <- aawcli[which(aa_supp > 0)]


View(ER_SNPs_merged_allbypos_modified[aacondcli,])


for(i in 1:length(aacondcli)){
  print(i)
  if(abs(ER_SNPs_merged_allbypos_modified$median_change[aacondcli[i]]) >= 0.04){
    aa_my_an_list <- list(aitrack, 
                          agtrack, 
                          aaot,
                          aasnp_track,
                          aaot_change,
                          aaexpression_change_median_track,
                          aa_tmp_erna_beforeE2_AnnotationTrack,
                          aa_tmp_erna_only_afterE2_AnnotationTrack,
                          aa_ER_AnnotationTrack,
                          aa_H3k27ac_AnnotationTrack)
    aacnt <- 1
    if(length(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacondcli[i]]]) > 0){
      aactf <- unique(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacondcli[i]]])
      for(j in 1:length(aactf)){
        if(aactf[j] %in% names(aa_remap_track_list)){
          aa_my_an_list[[length(aa_my_an_list) + 1]] <- aa_remap_track_list[[which(names(aa_remap_track_list) %in% aactf[j])]]
        }
      }
    }
    for(j in 1:length(ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacondcli[i]]])){
      aa_curenh <- which(names(aa_pos_neg) == ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacondcli[i]]][j])
      png(paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/Visual/4_",
                 ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacondcli[i]]][j] ,
                 ".png"),    # create PNG for the heat map        
          width = 8*300,        # 5 x 300 pixels
          height = 11*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8)        # smaller font size
      
      plotTracks(aa_my_an_list, 
                 from = start(aa_pos_neg)[aa_curenh], 
                 to = end(aa_pos_neg[aa_curenh]),
                 chromosome = as.character(seqnames(aa_pos_neg[aa_curenh])))
      dev.off()
    }
  }
  
}
################################################################################################################
# filter for variants that are worth exploring:

ER_SNPs_merged_allbypos_modified$is_clinical_pathogenic[is.na(ER_SNPs_merged_allbypos_modified$is_clinical_pathogenic)] <- F

aasum <- sum((ER_SNPs_merged_allbypos_modified$is_gtex + 
               ER_SNPs_merged_allbypos_modified$is_pancan + 
               ER_SNPs_merged_allbypos_modified$is_somatic + 
               ER_SNPs_merged_allbypos_modified$is_clinical_pathogenic) > 0)

aa_sta <- which((ER_SNPs_merged_allbypos_modified$is_gtex + 
                      ER_SNPs_merged_allbypos_modified$is_pancan + 
                      ER_SNPs_merged_allbypos_modified$is_somatic + 
                      ER_SNPs_merged_allbypos_modified$is_clinical_pathogenic) > 0)
aa_sta_chip <- numeric(length = length(aa_sta))
for(i in 1:length(aa_sta)){
  aa_sta_chip[i] <- sum(ER_SNPs_merged_allbypos_modified$Responsible_TF_support[[aa_sta[i]]])
}

aa_sta <- aa_sta[aa_sta_chip > 0]

sort(table(unlist(ER_SNPs_merged_allbypos_modified$Responsible_TF[aa_sta])))

aaposuniq <- ER_SNPs_merged_allbypos_modified$pos_uniq_id[aa_sta]
aaER_all_snp_investigate_mat <- ER_all_snp_investigate_mat
aaER_all_snp_investigate_mat$TF[aaER_all_snp_investigate_mat$TF == "ESR1"] <- "ESR1_2"

aaold_max <- numeric(length = nrow(ER_all_snp_investigate_mat))
aanew_max <- numeric(length = nrow(ER_all_snp_investigate_mat))
for(i in 1:length(aaold_max)){
  aaold_max[i] <- (aaER_all_snp_investigate_mat$old[i])/aamaxL_LLR[aaER_all_snp_investigate_mat$TF[i]]
  aanew_max[i] <- (aaER_all_snp_investigate_mat$new[i])/aamaxL_LLR[aaER_all_snp_investigate_mat$TF[i]]
}
ER_all_snp_investigate_mat$old_percent_max <- aaold_max
ER_all_snp_investigate_mat$new_percent_max <- aanew_max
hist(ER_all_snp_investigate_mat$old_percent_max[ER_all_snp_investigate_mat$TF == "ESR1"])

aa_int <- ER_all_snp_investigate_mat[ER_all_snp_investigate_mat$pos_uniq_id %in% aaposuniq,]
aa_int2 <- aa_int

aa_int2 <- aa_int2[(aa_int$old_percent_max >= 0.7 | aa_int$new_percent_max >= 0.7),]
aa_int2 <- aa_int2[abs(aa_int2$change_percent_max) >= 0.4,]
aa_int_ER_SNPs_merged_allbypos_modified <- ER_SNPs_merged_allbypos_modified[ER_SNPs_merged_allbypos_modified$pos_uniq_id %in% aa_int2$pos_uniq_id,]

View(aa_int_ER_SNPs_merged_allbypos_modified)
for(i in 1:nrow(aa_int_ER_SNPs_merged_allbypos_modified)){
  aa_c_p <- aa_int_ER_SNPs_merged_allbypos_modified$pos_uniq_id[i]
  aa_ctf <- aa_int2$TF[aa_int2$pos_uniq_id == aa_c_p]
  aawhich <- which(aa_int_ER_SNPs_merged_allbypos_modified$Responsible_TF[[i]] == aa_ctf)[1]
  aa_int_ER_SNPs_merged_allbypos_modified$Responsible_TF[[i]] <- aa_int_ER_SNPs_merged_allbypos_modified$Responsible_TF[[i]][aawhich]
  aa_int_ER_SNPs_merged_allbypos_modified$Responsible_TF_support[[i]] <- aa_int_ER_SNPs_merged_allbypos_modified$Responsible_TF_support[[i]][aawhich]
}
View(aa_int_ER_SNPs_merged_allbypos_modified)
aasup <- numeric(nrow(aa_int_ER_SNPs_merged_allbypos_modified))
for(i in 1:nrow(aa_int_ER_SNPs_merged_allbypos_modified)){
  aasup[i] <- sum(aa_int_ER_SNPs_merged_allbypos_modified$Responsible_TF_support[[i]])
}
aa_int_ER_SNPs_merged_allbypos_modified <- aa_int_ER_SNPs_merged_allbypos_modified[aasup > 0,]
View(aa_int_ER_SNPs_merged_allbypos_modified[,c("pos_uniq_id","is_gtex", "is_pancan","is_somatic", "is_common" ,"overlapping_enh", "median_change","Responsible_TF","LLR_percent_change_max","Enh_gene20kb_ovl", "Enh_geneProm_int")])
table(unlist(aa_int_ER_SNPs_merged_allbypos_modified$Responsible_TF))

aacond_int <- which(ER_SNPs_merged_allbypos_modified$pos_uniq_id %in% aa_int_ER_SNPs_merged_allbypos_modified$pos_uniq_id)

for(i in 1:length(aacond_int)){
  print(i)
  if(abs(ER_SNPs_merged_allbypos_modified$median_change[aacond_int[i]]) >= 0.15){
    aa_my_an_list <- list(aitrack, 
                          agtrack, 
                          aaot,
                          aasnp_track,
                          aaot_change,
                          aaexpression_change_median_track,
                          #aa_tmp_erna_beforeE2_AnnotationTrack,
                          aa_tmp_erna_only_afterE2_AnnotationTrack,
                          aa_ER_AnnotationTrack,
                          aa_H3k27ac_AnnotationTrack)
    aacnt <- 1
    if(length(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacond_int[i]]]) > 0){
      aactf <- unique(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacond_int[i]]])
      for(j in 1:length(aactf)){
        if(aactf[j] %in% names(aa_remap_track_list)){
          aa_my_an_list[[length(aa_my_an_list) + 1]] <- aa_remap_track_list[[which(names(aa_remap_track_list) %in% aactf[j])]]
        }
      }
    }
    for(j in 1:length(ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacond_int[i]]])){
      aa_curenh <- which(names(aa_pos_neg) == ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacond_int[i]]][j])
      png(paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/Visual/5_",
                 ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacond_int[i]]][j] ,
                 ".png"),    # create PNG for the heat map        
          width = 4*300,        # 5 x 300 pixels
          height = 9*300,
          res = 300,            # 300 pixels per inch
          pointsize = 9)        # smaller font size
      
      plotTracks(aa_my_an_list, 
                 from = start(aa_pos_neg)[aa_curenh], 
                 to = end(aa_pos_neg[aa_curenh]),
                 chromosome = as.character(seqnames(aa_pos_neg[aa_curenh])), 
                # background.panel = "#FFFEDB",
                 background.title = "darkblue")
      dev.off()
    }
  }
  
}

aa_int2_ref <- aa_int2[aa_int2$pos_uniq_id %in% aa_int_ER_SNPs_merged_allbypos_modified$pos_uniq_id,]
View(aa_int2_ref)
aa_int_ER_SNPs_merged_allbypos_modified$pos_uniq_id
aa_int_ER_SNPs_merged_allbypos_modified_posneg <- aa_int_ER_SNPs_merged_allbypos_modified
aa_int_ER_SNPs_merged_allbypos_modified_posneg <- aa_int_ER_SNPs_merged_allbypos_modified_posneg[ c(2,3,4,5,7, 9,10, 12, 13), ]
View(aa_int_ER_SNPs_merged_allbypos_modified_posneg[, c("pos_uniq_id","is_gtex", "is_pancan","is_somatic", "is_common" ,"overlapping_enh", "median_change","Responsible_TF","LLR_percent_change_max","Enh_gene20kb_ovl", "Enh_geneProm_int")])

View(ER_SNPs_merged_allbypos_modified[ER_SNPs_merged_allbypos_modified$seqnames == "chr12" &
                                        ER_SNPs_merged_allbypos_modified$start > 101697156 & 
                                        ER_SNPs_merged_allbypos_modified$start < 101697206,])

aatst_df <- data.frame(chr = "chr20", start=50709255-10,end=50709255+10)
aatst_gr <- makeGRangesFromDataFrame(aatst_df)
aa <- as.character(getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = aatst_gr))
seqLogo::seqLogo(reverseComplement(TF.motifs.Expanded_new_pseudo_t$FOXA1))

##################################################################################################################
# create the the list of finalists for variant impact testing:
aa_posuniqud <- c(32887, 6840, 5985, 7000, 11642, 31880, 27789)
View(ER_SNPs_merged_allbypos_modified[match(aa_posuniqud, ER_SNPs_merged_allbypos_modified$pos_uniq_id),])
View(ER_all_snp_investigate_mat[ER_all_snp_investigate_mat$pos_uniq_id %in% aa_posuniqud,])


################################################################################################################
################################################################################################################
# plot the extended version of neg_383, pos_118, and pos_197
THree_chosen_extended_seq[1]
nchar(THree_chosen_extended_seq[1])


aabinding_ylim <- rbind(c(-0.02, 1.02), c(-0.02, 1.02), c(-0.02, 0.85))
aabindingchange_ylim <- rbind(c(-0.02, 0.4), c(-0.5, 0.02), c(-0.7, 0.02))
aaactivity_ylim <- rbind(c(-0.02,0.15), c(-0.2, 0.02), c(-0.2, 0.02))
rownames(aabinding_ylim) <- names(THree_chosen_extended_seq)
rownames(aabindingchange_ylim) <-  names(THree_chosen_extended_seq)
rownames(aaactivity_ylim) <-  names(THree_chosen_extended_seq)

displayPars(aa_tmp_erna_only_afterE2_AnnotationTrack) <- list(rotation.title = 0, size = 0.3)
displayPars(aa_ER_AnnotationTrack) <- list(rotation.title = 0, size = 0.3)
displayPars(aa_H3k27ac_AnnotationTrack) <- list(rotation.title = 0, size = 0.3)

for(i in 1:length(aa_remap_track_list)){
  displayPars(aa_remap_track_list[[i]]) <- list(rotation.title = 0, size = 0.3)
}

displayPars(aasnp_track) <-  list(rotation.title = 0, size = 0.3)
displayPars(agtrack) <- list(scale = 0.2)
aammap <- c(2, 1, 3)
for(i in 1:length(aacond_int_chosen)){
  print(i)
  if(abs(ER_SNPs_merged_allbypos_modified$median_change[aacond_int_chosen[i]]) >= 0.15){
    
    for(aactrack in 1:length(aatracklist)){
      displayPars(aatracklist[[aactrack]]) <- list(alpha.title = 1, alpha = 1,
                                                   ylim = aabinding_ylim[aammap[i],],
                                                   size = 1)
    }
    aaot <- OverlayTrack(trackList = aatracklist)
    
    for(aactrack in 1:length(aabinding_change_track_list)){
      displayPars(aabinding_change_track_list[[aactrack]]) <- list(alpha.title = 1, alpha = 1, 
                                                                   lwd.grid = 0.65, lty.grid = 2, 
                                                                   col.grid = "lightgray",
                                                                   ylim = aabindingchange_ylim[aammap[i],], 
                                                                   size = 0.4, 
                                                                   rotation.title = 90)
    }
    aaot_change <- OverlayTrack(trackList = aabinding_change_track_list)
    
    displayPars(aaexpression_change_median_track) <- list(alpha.title = 1, alpha = 1, 
                                                          lwd.grid = 0.65, lty.grid = 2, 
                                                          col.grid = "lightgray",
                                                          ylim = aaactivity_ylim[aammap[i], ],
                                                          size = 0.4,
                                                          rotation.title = 90)
    aacond_int_chosen <- aacond_int[c(3,4,7)]
    
    
    aa_my_an_list <- list(aitrack, 
                          agtrack, 
                          aaot,
                          aasnp_track,
                          aaot_change,
                          aaexpression_change_median_track,
                          #aa_tmp_erna_beforeE2_AnnotationTrack,
                          aa_tmp_erna_only_afterE2_AnnotationTrack,
                          aa_ER_AnnotationTrack,
                          aa_H3k27ac_AnnotationTrack)
    aacnt <- 1
    if(length(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacond_int_chosen[i]]]) > 0){
      aactf <- unique(ER_SNPs_merged_allbypos_modified$Responsible_TF[[aacond_int_chosen[i]]])
      for(j in 1:length(aactf)){
        if(aactf[j] %in% names(aa_remap_track_list)){
          aa_my_an_list[[length(aa_my_an_list) + 1]] <- aa_remap_track_list[[which(names(aa_remap_track_list) %in% aactf[j])]]
        }
      }
    }
    for(j in 1:length(ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacond_int_chosen[i]]])){
      aa_curenh <- which(names(THree_chosen_extended) == ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacond_int_chosen[i]]][j])
      if(length(aa_curenh) > 0){
        png(paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/Visual/7_",
                   ER_SNPs_merged_allbypos_modified$overlapping_enh[[aacond_int_chosen[i]]][j] ,
                   ".png"),    # create PNG for the heat map        
            width = 4*300,        # 5 x 300 pixels
            height = 6*300,
            res = 300,            # 300 pixels per inch
            pointsize = 9)        # smaller font size
        
        plotTracks(aa_my_an_list, 
                   from = start(THree_chosen_extended)[aa_curenh], 
                   to = end(THree_chosen_extended[aa_curenh]),
                   chromosome = as.character(seqnames(THree_chosen_extended[aa_curenh])), 
                   # background.panel = "#FFFEDB",
                   background.title = "darkblue",
                   title.width = 1.2)
        dev.off()
      }

    }
  }
  
}


# compute the number of TPs if sampling randomly from common snps and the eqtls
ER_gtex_eqtl_sig
ER_pancan_eqtl_cis
ER_somatic_mutation_all


aa_clinvar_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/ClinVar/clinvar_hg38.bed"
aa_clinvar <- readPeakFile(aa_clinvar_Add, as = "GRanges")

aauniq <- unique(ER_gtex_eqtl_sig$variant_id)
unlist(lapply(strsplit(aauniq, "_"), "[[", 1))

aa_qtl_df <- data.frame(chr = unlist(lapply(strsplit(aauniq, "_"), "[[", 1)),
                        start = as.numeric(unlist(lapply(strsplit(aauniq, "_"), "[[", 2))), 
                        end = as.numeric(unlist(lapply(strsplit(aauniq, "_"), "[[", 2))) + 1, 
                        gtex_id = c(1:length(aauniq)))
aa_qtl_gr <- makeGRangesFromDataFrame(aa_qtl_df)
write.table(aa_qtl_gr, 
            file="~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/sig_eqtl_uniq.bed",
            quote=F, sep="\t", row.names=F, col.names=F)

nu_uniq_common_snp <- 37302978 - 705534
nu_uniq_gtex_breast_Cancer_eqtls <- length(aauniq)
nu_uniq_gtex_overlap_common <- 870778
nu_unique_pancan_eqtls <- length(ER_pancan_cis_trans_gr38_unique)
nu_uniq_pancan_overlap_common <- 327977
aa_gtex_random_ratio <- nu_uniq_gtex_breast_Cancer_eqtls/(nu_uniq_common_snp + (nu_uniq_gtex_breast_Cancer_eqtls - nu_uniq_gtex_overlap_common))
nu_somatic_breast_cancer <- nrow(ER_somatic_mutation_all)
nu_somatic_common_overlap <- 97529
nu_pancan_gtex_overlap <- 173402
nu_clinvar_vars <- 655934
nu_clinvar_vars_common_intersect <- 150187

# simulate random
aa <- numeric(length(nu_uniq_common_snp + (nu_uniq_gtex_breast_Cancer_eqtls - nu_uniq_gtex_overlap_common)))
aa[1:nu_uniq_gtex_breast_Cancer_eqtls] <- 1



ER_somatic_mutation_all$FATHMM_MKL_NON_CODING_SCORE <- as.numeric(ER_somatic_mutation_all$FATHMM_MKL_NON_CODING_SCORE)
nu_somatic_breast_cancer_fathmmScored <- sum(!is.na(ER_somatic_mutation_all$FATHMM_MKL_NON_CODING_SCORE))
nu_somatic_breast_cancer_fathmmScored_sig7 <- sum(ER_somatic_mutation_all$FATHMM_MKL_NON_CODING_SCORE > 0.7, na.rm = T)




aaa1 <- "AAGCAGAGCTTGCAGTGAGCCGAGATCGCACCACTGCATTCCAGCCTGGGTGACAGAACAAGACTCCATCTCAAATAAAACAAAAAAAAGATATTCTCATACTTCGGCCTCCTGAGTAGCTGGGATTACAGGCACCAGCCACCACGATGGGCTAATTTTTTGTATTTTTAGTAGAGATGGGATTTCACCATGTTGGCCAGGGTGGTCTCAACCTCCTGGCCTCAAGTGATCCCCCTCCTCGACCTCCCAAAATGCTGGGATTACAGGCATAGGCTGCTGTGCCCAGACGTGAGTTACGTTTCAAAAGGAGTGGGGGATCTCCAAGTGGCTCCCCACAAGTGGATCTGCAAGTTCTCACTCCCCTCTTTCTGGGACTCAACTCAAAGTCTCCTCTTCCCTGTTTTATTGTCATCACTGCGCTTGTCAGGAACTGAAATGACCTCCTTTGCTTGTCTGTTGATCGACTGTCTCTCGTCCCCGGCCCCAGCGCCTCAAATGGCAGGGAGACCTATGTAGGTCGGCGTCGGGCCAACCCTCTCCTGTGCCTCTCGCCTGGGGCCTGGCACTGCGCAGGGGTCAATAAACGGGTGGTTGGCTAACTGCCCTGCCACCGTGCCTGCCAGAGCGGCCTCCTCCACCCTCCCTGGGTGCCACTCGGCCCCCACTCTGGTGCCAAGGCCCTGGGTTTGCCTCTTGGTCTTTGGCCAATCGGTCAACAATGCCCCTTGGGCACCGGGCACTTCCCCCACCCTACAAAGCGGCAAGTGGGCGACTGCGCATGCCCCTCCTGAGCCTCCACCAGGCCCAGTGAGGCGAGGAAGTGACCCCACCGCGGCCTCTCGGCTGCCAACCCCAGAACGCAGATCGAGCCCCTCAACGCCTCGGAGCCCCCGGGCCTCCTTCTGCGCAGCGCCTCCTGGTTCTCTCGCTCGGTTTCCGCCGCTCCCGGGCCACCCCAAGAGTTGCGCCGTCCAACAAAAATATCCACCCAGGTCGCGCCCTCCTGGCCCTCCGAGTCCTGGCGGGGAACTCGGAGGCGCCCTGGGCCGTGGCGCCCGTTTTCCTGGTCTCCGGGCACGGCCGCAGGGGCGCTCGCGGCCCGCGATGCCCGCAGCCCGGGC"
aaa2 <- "AAGCAGAGCTTGCAGTGAGCCGAGATCGCACCACTGCATTCCAGCCTGGGTGACAGAACAAGACTCCATCTCAAATAAAACAAAAAAAAGATATTCTCATACTTCGGCCTCCTGAGTAGCTGGGATTACAGGCACCAGCCACCACGATGGGCTAATTTTTTGTATTTTTAGTAGAGATGGGATTTCACCATGTTGGCCAGGGTGGTCTCAACCTCCTGGCCTCAAGTGATCCCCCTCCTCGACCTCCCAAAATGCTGGGATTACAGGCATAGGCTGCTGTGCCCAGACGTGAGTTACGTTTCAAAAGGAGTGGGGGATCTCCAAGTGGCTCCCCACAAGTGGATCTGCAAGTTCTCACTCCCCTCTTTCTGGGACTCAACTCAAAGTCTCCTCTTCCCTGTTTTATTGTCATCACTGCGCTTGTCAGGAACTGAAATGACCTCCTTTGCTTGTCTGTTGATCGACTGTCTCTCGTCCCCGGCCCCAGCGCCTCAAATGGCAGGGAGACCTATGTAGGTCGGCGTCGGGCCAACCCTCTCCTGTGCCTCTCGCCTGGGGCCTGGCACTGCGCAGGGGTCAATAAACGGGTGGTTGGCTAACTGCCCTGCCACCGTGCCTGCCAGAGCGGCCTCCTCCACCCTCCCTGGGTGCCACTCGGCCCCCACTCTGGTGCTCGACCCCTGGGTTTGCCTCTTGGTCTTTGGCCAATCGGTCAACAATGCCCCTTGGGCACCGGGCACTTCCCCCACCCTACAAAGCGGCAAGTGGGCGACTGCGCATGCCCCTCCTGAGCCTCCACCAGGCCCAGTGAGGCGAGGAAGTGACCCCACCGCGGCCTCTCGGCTGCCAACCCCAGAACGCAGATCGAGCCCCTCAACGCCTCGGAGCCCCCGGGCCTCCTTCTGCGCAGCGCCTCCTGGTTCTCTCGCTCGGTTTCCGCCGCTCCCGGGCCACCCCAAGAGTTGCGCCGTCCAACAAAAATATCCACCCAGGTCGCGCCCTCCTGGCCCTCCGAGTCCTGGCGGGGAACTCGGAGGCGCCCTGGGCCGTGGCGCCCGTTTTCCTGGTCTCCGGGCACGGCCGCAGGGGCGCTCGCGGCCCGCGATGCCCGCAGCCCGGGC"
aaa3 <- "AAGCAGAGCTTGCAGTGAGCCGAGATCGCACCACTGCATTCCAGCCTGGGTGACAGAACAAGACTCCATCTCAAATAAAACAAAAAAAAGATATTCTCATACTTCGGCCTCCTGAGTAGCTGGGATTACAGGCACCAGCCACCACGATGGGCTAATTTTTTGTATTTTTAGTAGAGATGGGATTTCACCATGTTGGCCAGGGTGGTCTCAACCTCCTGGCCTCAAGTGATCCCCCTCCTCGACCTCCCAAAATGCTGGGATTACAGGCATAGGCTGCTGTGCCCAGACGTGAGTTACGTTTCAAAAGGAGTGGGGGATCTCCAAGTGGCTCCCCACAAGTGGATCTGCAAGTTCTCACTCCCCTCTTTCTGGGACTCAACTCAAAGTCTCCTCTTCCCTGTTTTATTGTCATCACTGCGCTTGTCAGGAACTGAAATGACCTCCTTTGCTTGTCTGTTGATCGACTGTCTCTCGTCCCCGGCCCCAGCGCCTCAAATGGCAGGGAGACCTATGTAGGTCGGCGTCGGGCCAACCCTCTCCTGTGCCTCTCGCCTGGGGCCTGGCACTGCGCAGGGGTCAATAAACGGGTGGTTGGCTAACTGCCCTGCCACCGTGCCTGCCAGAGCGGCCTCCTCCACCCTCCCTGGGTGCCACTCGGCCCCCACTCTGGTGCTCGACCCCTGGGTTTGCCTCTTGGTCTTTGGCCAATCGGTCAACAATGCCCCTTGGGCACCGGGCACTTCCCCCACCCTACAAAGCGGCAAGTGGGCGACTGCGCATGCCCCTCCTGAGCCTCCACCAGGCCCAGTGAGGCGAGGAAGTGACCCCACCGCGGCCTCTCGGCTGCCAACCCCAGAACGCAGATCGAGCCCCTCAACGCCTCGGAGCCCCCGGGCCTCCTTCTGCGCAGCGCCTCCTGGTTCTCTCGCTCGGTTTCCGCCGCTCCCGGGCCACCCCAAGAGTTGCGCCGTCCAACAAAAATATCCACCCAGGTCGCGCCCTCCTGGCCCTCCGAGTCCTGGCGGGGAACTCGGAGGCGCCCTGGGCCGTGGCGCCCGTTTTCCTGGTCTCCGGGCACGGCCGCAGGGGCGCTCGCGGCCCGCGATGCCCGCAGCCCGGGC"
as.character(THree_chosen_extended_seq[3]) == aaa1
aaa2 == aaa3



