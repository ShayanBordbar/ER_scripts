# visializing enhancers with snps on them
# What are the things that you want to be visualized?
# sites above threshold on the enhancer for all TFs --> can be used from GEMSTAT annotations
# change in binding 
# chip tracks where available
# e-RNA where avialable
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
library(Gviz)
library(RColorBrewer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# I will use GViz package for visualizing these datasets

#install the sushi 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Sushi")
# library(Sushi)
# Sushi_data = data(package = 'Sushi')
# data(list = Sushi_data$results[,3])
# Sushi_data$results[,3]
# head(Sushi_DNaseI.bedgraph)
#################################################################################################################
#################################################################################################################
#####################################            Functions             ##########################################
#################################################################################################################
#################################################################################################################
# write a function to take in annotations from GEMSTAT, threshold on LLR for calling sites, MAXLLR for each TF,
# the position of enhancers and 
# write a bed file with intensity being the LLR for each TF separately
write_bed_from_annotation <- function(annotation_file,
                                      TF_names,
                                       MaxLLR,
                                      # LLR_thresh,
                                      enhancer_gr){
  # annotation_file : is the address of output of GEMSTAT seqannot
  # TF_names : is a character vector with names of TFs that we want to write a bed file for
  # MaxLLR : is a named numeric vector the MAX LLR that can be achieved for the TFs motif
#  # LLR_thresh : is a named numeric vector the minimum accepter LLR for the TFs
  # enhancer_gr: is a Granges object containing the chr, start and end for all enhancers
  # this function assumes that the annotations have been trimmed for the right thresholds --> thresholds can be givn as input to seqannot
  
  stopifnot(#length(TF_names) == length(MaxLLR),
#            length(TF_names) == length(LLR_thresh),
            all(names(MaxLLR) == TF_names)#,
            #all(names(LLR_thresh) == TF_names)
            )
  
  my_annot <- annotation_reader(annot_file =annotation_file, 
                                        TF_names = TF_names)
  stopifnot(all(TF_names %in% names(my_annot[[1]])), 
            all(names(my_annot) %in% names(enhancer_gr)))
  TF_df_list <- list() # this will contain one entry for each TF: That entry is a dataframe with chr, start, end and 1 - (LLRscore/LLRmax)
  for(i in 1:length(TF_names)){
    TF_df_list[[i]] <- matrix(nrow = 0, ncol = 5)
    colnames(TF_df_list[[i]]) <- c("chr", "start", "end", "strand", "score")
  }
  for(cur_enh in 1:length(my_annot)){
    for(cur_tf in 1:length(my_annot[[cur_enh]])){
      if(nrow(my_annot[[cur_enh]][[cur_tf]]) > 0){
        cur_star <- as.numeric(unlist(lapply(strsplit(my_annot[[cur_enh]][[cur_tf]]$V1, "\\.."), "[[", 1)))
        cur_stop <- as.numeric(unlist(lapply(strsplit(my_annot[[cur_enh]][[cur_tf]]$V1, "\\.."), "[[", 2)))
        cur_score <- 1-  my_annot[[cur_enh]][[cur_tf]]$V4/MaxLLR[cur_tf]
        gr_index <- which(names(enhancer_gr) %in% names(my_annot)[cur_enh])
        enh_St <- start(enhancer_gr[gr_index])
        cur_star <- cur_star + (enh_St  - 1)
        cur_stop <- cur_stop + (enh_St  - 1)
        cur_strand <- my_annot[[cur_enh]][[cur_tf]]$V2
        cur_chr <- rep(as.character(seqnames(enhancer_gr[gr_index])), length(cur_star))
        cur_df <- cbind(cur_chr, cur_star, cur_stop, cur_strand,cur_score)
        TF_df_list[[cur_tf]] <- rbind(TF_df_list[[cur_tf]], cur_df)
      }
    } # end of loop over TFs
  }# end of loop over enhancers
  names(TF_df_list) <- TF_names
  TF_df_list_df <- list()
  for(i in 1:length(TF_df_list)){
    TF_df_list_df[[i]] <- data.frame(chr = TF_df_list[[i]][, 1],
                                     start = as.numeric(TF_df_list[[i]][, 2]), 
                                     end = as.numeric(TF_df_list[[i]][, 3]), 
                                     strand = TF_df_list[[i]][, 4],score = as.numeric(TF_df_list[[i]][, 5]))
  }
  names(TF_df_list_df) <- TF_names
  return(TF_df_list_df)
}
#################################################################################################################
# example
aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)

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
cat(aathr[aa_names],sep = "\n", file = "~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/ann_thr.txt")
system("~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/annotation_do.sh")

# aa_pos_neg_annot2 <- annotation_reader(annot_file ="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/pos_neg_2118_th.ann", 
#                                        TF_names = names(TF.motifs.Expanded_pseudo_exp12))

aa_names <- c("ESR1_2","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","RARA", "RELA","RUNX1","SP1","TFAP2C","YBX1")

aa_test <- write_bed_from_annotation(annotation_file = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/pos_neg_2118_th.ann",
                                     TF_names = aa_names,
                                     MaxLLR = aamaxL_LLR[aa_names],
                                     # LLR_thresh,
                                     enhancer_gr = aa_pos_neg)
aa_test$ESR1_2[1:10,]
#################################################################################################################
# get annotations for all sequences from GEMSTAT

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)
aa_pos_neg_Seq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = aa_pos_neg)
pos_neg_set4_2118_sequence <- aa_pos_neg_Seq


writeFasta(aa_pos_neg_Seq, width = 1000,
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/pos_neg_2118.fa")

aa_names <- c("ESR1_2","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","RARA", "RELA","RUNX1","SP1","TFAP2C","YBX1")



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
aathr[aa_names]
aamaxL_LLR[aa_names] -  aamin_LLR[aa_names]
cat(aathr[aa_names],sep = "\n", file = "~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/ann_thr.txt")


WT_binding_bed_df <- write_bed_from_annotation(annotation_file = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/pos_neg_2118_th.ann",
                                     TF_names = aa_names,
                                     MaxLLR = aamaxL_LLR[aa_names],
                                     # LLR_thresh,
                                     enhancer_gr = aa_pos_neg)
aax <- WT_binding_bed_df
aax$ESR1_2$strand <- as.character(levels(aax$ESR1_2$strand)[as.numeric(aax$ESR1_2$strand)])
aax$ESR1_2$strand[aax$ESR1_2$strand == "-"] <- -1
aax$ESR1_2$strand[aax$ESR1_2$strand == "+"] <- 1
aax$ESR1_2$strand <- as.numeric(aax$ESR1_2$strand)

aax$ESR1_2[aax$ESR1_2$chr == as.character(seqnames(aa_pos_neg[i])) & aax$ESR1_2$start > start(aa_pos_neg)[i] & aax$ESR1_2$start < end(aa_pos_neg)[i],]



########################################################################################################################
# get fasta sequence for the extended versions of the three chosen enhancers:
#neg_383
#chr12:101696596-101697634
#pos_118
#chr14:77399430-77400995
#pos_197
#chr19:10501558-10502678

library(Biostrings)
THree_chosen_extended <- makeGRangesFromDataFrame(df = data.frame(chr = c("chr12", "chr14", "chr19"),
                                                                  start = c(101696596, 77399430, 10501558),
                                                                  end = c(101697634, 77400995, 10502678)))
names(THree_chosen_extended) <- c("neg_383", "pos_118", "pos_197")
THree_chosen_extended_seq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = THree_chosen_extended)


library(ShortRead)
writeFasta(THree_chosen_extended_seq, width = 1000,
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/neg_383_pos_118_pos_197_extended.fa")
aa_names <- c("ESR1_2","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","RARA", "RELA","RUNX1","SP1","TFAP2C","YBX1")


WT_binding_bed_Extended_three <- write_bed_from_annotation(annotation_file = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/pos_neg_extended_3_th.ann",
                                               TF_names = aa_names,
                                               MaxLLR = aamaxL_LLR[aa_names],
                                               # LLR_thresh,
                                               enhancer_gr = THree_chosen_extended)
WT_binding_bed_df2 <- WT_binding_bed_df
for(i in 1:length(WT_binding_bed_df)){
  WT_binding_bed_df[[i]] <- WT_binding_bed_df[[i]][!duplicated(WT_binding_bed_df[[i]]),]
  if(nrow(WT_binding_bed_Extended_three[[i]]) > 0){
    for(j in 1:nrow(WT_binding_bed_Extended_three[[i]])){
      aatmp <- rbind(WT_binding_bed_df[[i]], WT_binding_bed_Extended_three[[i]][j,])
      if(sum(duplicated(aatmp)) == 0){
        print(paste0("added ", names(WT_binding_bed_df)[i], " ", j))
        WT_binding_bed_df[[i]] <- aatmp
      }
    }
  }

}

WT_binding_bed_Extended_three$ESR1_2 %in% WT_binding_bed_df$ESR1_2


# check if the variants have created any new sites:
aa_annot <- annotation_reader(annot_file ="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Enhancer_visual/pos_neg_extended_3_th_var_wt.ann", 
                              TF_names = aa_names)

aa_ident <- matrix(nrow = 3, ncol = 17)
colnames(aa_ident) <- names(aa_annot$neg_383_wt)
rownames(aa_ident) <- c("neg_383", "pos_118", "pos_197")
for(i in 1:3){
  for(j in 1:length(aa_annot[[2*i]])){
    aa_ident[i, j] <- identical(aa_annot[[(2*i) - 1]][[j]], aa_annot[[(2*i)]][[j]])
  }
}
rowSums(aa_ident)
identical(aa_annot$neg_383_wt$PGR, aa_annot$neg_383_variant$PGR)
identical(aa_annot$pos_118_wt$JUN_1, aa_annot$pos_118_variant$JUN_1)

aa_annot$pos_118_wt$NKX3_1
aa_annot$pos_118_variant$NKX3_1
aasq <- as.character(THree_chosen_extended_seq)
substr(aasq[2], 1149, 1161)


########################################################################################################################

# visualization one enhancer as an example and using it to create a function to visualize any enhancer
# ReMapChIP.GRanges.list is the one containing granges for all remap chips

library(RColorBrewer)
n <- 17
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
aaaaaacolor = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]


#display.brewer.all() 
#colorRampPalette(brewer.pal(11, "Set3"))(17)

#aamycoltest <- colorRampPalette(brewer.pal(11, "Set3"))(17)
#aamycoltest <- sample(x = col_vector, size = n, replace = F)
aamycoltest <- sample(x = aaaaaacolor, size = n, replace = F)
barplot(rep(5, 17), names.arg = c(1:17), col = aamycoltest, las = 2)


aatracklist <- list()
aamycol <- aamycoltest
for(i in 1:length(WT_binding_bed_df)){
  aatracklist[[i]] <- DataTrack(data = WT_binding_bed_df[[i]]$score,
                                start = WT_binding_bed_df[[i]]$start, 
                                end = WT_binding_bed_df[[i]]$end,
                                chromosome = WT_binding_bed_df[[i]]$chr,
                                genome = "hg38", 
                                name = "Binding strength",  
                                groups = factor(names(WT_binding_bed_df)[i], 
                                                                levels = names(WT_binding_bed_df)), 
                                col = aamycol,
                                legend = TRUE, 
                                ylim = c(-0.02, 1.02),
                                cex.legend = 0.7, 
                                type= "hist")
}

#for(i in 1:length(aatracklist)){
#  displayPars(aatracklist[[i]]) <- list(alpha.title = 1, alpha = 1)
#}


aaot <- OverlayTrack(trackList = aatracklist)
agtrack <- GenomeAxisTrack(cex = 0.55)
aachr <- as.character(unique(WT_binding_bed_df$ESR1_2$chr))
aastrack <- SequenceTrack(Hsapiens, chromosome = aachr)
aitrack <- IdeogramTrack(genome = "hg38", chromosome = aachr)
aaylims <- extendrange(c(0,1))

which(names(aa_pos_neg) == "neg_384")
which(names(aa_pos_neg) == "pos_118")
which(names(aa_pos_neg) == "pos_197")
i <- 197
plotTracks(list(aitrack, agtrack, aaot, aastrack),
           from = start(aa_pos_neg)[i], 
           to = end(aa_pos_neg[i]), 
           chromosome = as.character(seqnames(aa_pos_neg[i]))
          # ylim = aaylims
           #, type = c("hist")
           #, window = 30, 
          #,background.panel = "#FFFEDB",
          ,background.title = "darkblue"
           )
#  snp
aasnp_track <- AnnotationTrack(ER_SNPs_merged_allbypos_GR, name = "Variants", genome = "hg38", stacking = "dense")
i <- 1
plotTracks(list(aitrack, agtrack, aaot, aasnp_track), from = start(aa_pos_neg)[i], to = end(aa_pos_neg[i]), chromosome = as.character(seqnames(aa_pos_neg[i])),
           ylim = aaylims, type = c("hist")
           #, window = 30
)
#  change of binding due to snp
#ER_all_snp_investigate_mat_withcoord
#Change_in_binding_snp_byTF
aabinding_change_track_list <- list()
for(i in 1:length(WT_binding_bed_df)){
  aabinding_change_track_list[[i]] <- DataTrack(data = Change_in_binding_snp_byTF[[i]]$change_percent_max,
                                                start = Change_in_binding_snp_byTF[[i]]$start, 
                                end = Change_in_binding_snp_byTF[[i]]$end, 
                                chromosome = Change_in_binding_snp_byTF[[i]]$seqnames,
                                genome = "hg38", 
                                name = "Binding\n change", 
                                groups = factor(names(Change_in_binding_snp_byTF)[i], 
                                                                            levels = names(Change_in_binding_snp_byTF)), 
                                col = aamycol,
                                legend = FALSE,
                                ylim = c(-0.7, 0.4), 
                                grid = T,
                                type = "hist")
}

for(i in 1:length(aabinding_change_track_list)){
  displayPars(aabinding_change_track_list[[i]]) <- list(alpha.title = 1, alpha = 1, lwd.grid = 0.65, lty.grid = 2, col.grid = "lightgray")
}


aaot_change <- OverlayTrack(trackList = aabinding_change_track_list)

#displayPars(aaot_change) <- list(ylim = c(-1.5, 1.5))
aaylims <- extendrange(c(0,1))
i <- 197
plotTracks(list(aitrack, agtrack, aaot, aasnp_track, aaot_change), 
           from = start(aa_pos_neg)[i], 
           to = end(aa_pos_neg[i]),
           chromosome = as.character(seqnames(aa_pos_neg[i])),
           background.title = "darkblue"
          # ylim = aaylims,
           #type = c("hist")
           #, window = 30
)



# change of prediction due to snp: histplot (multiple alternatives are allowed) or boxplot(one alternative has to be chosen)

aaexpression_change_mean_track <- DataTrack(data = ER_SNPs_merged_allbypos_modified_and_modelScores$mean_change,
                                            start = ER_SNPs_merged_allbypos_modified_and_modelScores$start, 
                                            end = ER_SNPs_merged_allbypos_modified_and_modelScores$end,
                                            chromosome = ER_SNPs_merged_allbypos_modified_and_modelScores$seqnames, 
                                            genome = "hg38", 
                                            name = "mean",  
                                            type = "hist",grid = T,
                                            #ylim = extendrange(range(ER_SNPs_merged_allbypos_modified_and_modelScores$mean_change, na.rm = T))
                                            ylim = c(-0.1,0.1)
)
displayPars(aaexpression_change_mean_track) <- list(alpha.title = 1, alpha = 0.9, lwd.grid = 0.5, lty.grid = 2, col.grid = "lightgray")

aaexpression_change_median_track <- DataTrack(data = ER_SNPs_merged_allbypos_modified_and_modelScores$median_change,
                                              start = ER_SNPs_merged_allbypos_modified_and_modelScores$start, 
                                              end = ER_SNPs_merged_allbypos_modified_and_modelScores$end,
                                              chromosome = ER_SNPs_merged_allbypos_modified_and_modelScores$seqnames, 
                                              genome = "hg38", 
                                              name = "Activity\n change",  
                                              type = "hist",grid = T,
                                              ylim = c(-0.2,0.15)
                                              #ylim = extendrange(range(ER_SNPs_merged_allbypos_modified_and_modelScores$median_change, na.rm = T))
)
displayPars(aaexpression_change_median_track) <- list(alpha.title = 1, alpha = 1, lwd.grid = 0.65, lty.grid = 2, col.grid = "lightgray")

aaexpression_change_max_track <- DataTrack(data = ER_SNPs_merged_allbypos_modified_and_modelScores$max_change,
                                           start = ER_SNPs_merged_allbypos_modified_and_modelScores$start, 
                                           end = ER_SNPs_merged_allbypos_modified_and_modelScores$end,
                                           chromosome = ER_SNPs_merged_allbypos_modified_and_modelScores$seqnames, 
                                           genome = "hg38", 
                                           name = "max",  
                                           type = "hist",grid = T,
                                           ylim = c(-0.1,0.1)
                                           #ylim = extendrange(range(ER_SNPs_merged_allbypos_modified_and_modelScores$max_change, na.rm = T))
)
displayPars(aaexpression_change_max_track) <- list(alpha.title = 1, alpha = 0.9, lwd.grid = 0.5, lty.grid = 2, col.grid = "lightgray")

aaexpression_change_percent_track <- DataTrack(data = ER_SNPs_merged_allbypos_modified_and_modelScores$percent_changed,
                                               start = ER_SNPs_merged_allbypos_modified_and_modelScores$start, 
                                               end = ER_SNPs_merged_allbypos_modified_and_modelScores$end,
                                               chromosome = ER_SNPs_merged_allbypos_modified_and_modelScores$seqnames, 
                                               genome = "hg38", 
                                               name = "percent",  
                                               type = "hist",grid = T,
                                               #ylim = extendrange(range(ER_SNPs_merged_allbypos_modified_and_modelScores$percent_changed, na.rm = T))
                                               ylim = c(-0.2, 0.2)
)

displayPars(aaexpression_change_percent_track) <- list(alpha.title = 1, alpha = 0.9, lwd.grid = 0.5, lty.grid = 2, col.grid = "lightgray")


i <- 197
plotTracks(list(aitrack, 
                agtrack,
                aaot,
                aasnp_track,
                aaot_change, 
                aaexpression_change_mean_track,
                aaexpression_change_median_track,
                aaexpression_change_max_track, 
                aaexpression_change_percent_track), 
           from = start(aa_pos_neg)[i], 
           to = end(aa_pos_neg[i]),
           chromosome = as.character(seqnames(aa_pos_neg[i])),
           # ylim = aaylims,
           type = c("hist")
           #, window = 30
)

aapp <- unlist(lapply(strsplit(colnames(ER_SNPs_merged_allbypos_modified_and_modelScores), "_"), "[[", 1))
aacolmn <- which(aapp == "par")
aa_datmat <- as.matrix(ER_SNPs_merged_allbypos_modified_and_modelScores[One_alternative_index,aacolmn])
aaexpression_change_boxplot_track <- DataTrack(data = t(aa_datmat),
                                               start = ER_SNPs_merged_allbypos_modified_and_modelScores$start[One_alternative_index], 
                                               end = ER_SNPs_merged_allbypos_modified_and_modelScores$end[One_alternative_index],
                                               chromosome = ER_SNPs_merged_allbypos_modified_and_modelScores$seqnames[One_alternative_index], 
                                               genome = "hg38", 
                                               name = "change",  
                                               # groups = factor(names(Change_in_binding_snp_byTF)[i], 
                                               #                levels = names(Change_in_binding_snp_byTF)), 
                                               #col = aamycol,
                                               # legend = TRUE,
                                               type = "boxplot",grid = T,
                                               ylim = c(-0.1,0.1)
                                               #ylim = extendrange(range(aa_datmat, na.rm = T))
)
displayPars(aaexpression_change_boxplot_track) <- list(lwd.grid = 0.5, lty.grid = 2, col.grid = "lightgray")

i <- 118
plotTracks(list(aitrack, 
                agtrack, 
                aaot,
                aasnp_track,
                aaot_change, 
                aaexpression_change_mean_track,
                aaexpression_change_median_track,
                aaexpression_change_max_track,
#                aaexpression_change_percent_track, 
                aaexpression_change_boxplot_track
                ), 
           from = start(aa_pos_neg)[i], 
           to = end(aa_pos_neg[i]),
           chromosome = as.character(seqnames(aa_pos_neg[i]))
           # ylim = aaylims,
           #type = c("hist")
           #, window = 30
)

########
# remap chip
load("~/Documents/Shayan/BioInf/EstrogenReceptor/Validation/Remap_GRanges_merged.RData")
i <- 54
#length(ReMapChIP.GRanges.list_merged[[i]])
#length(ReMapChIP.GRanges.list[[i]])



names(ReMapChIP.GRanges.list_merged) <- names(ReMapChIP.GRanges.list)
aa_remap_track_list <- list()
aach <- names(ReMapChIP.GRanges.list)
aa_names <- c("ESR1","FOXA1","GATA3","JUN","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","RARA", "RELA","RUNX1","SP1","TFAP2C","YBX1")
intersect(aach,aa_names)
setdiff(aa_names,aach)
aaint <- intersect(aach,aa_names)
aaint <- c(aaint, "NCOA1", "NCOA2", "NCOA3", "AR",  "RAD21", "EP300")
aaint <- aaint[-c(3)]

# get the intersection with the enhancers and plot those only
# for(i in 1:length(aaint)){
#   print(i)
#   print(aaint[i])
#   print(length(ReMapChIP.GRanges.list_merged[[aaint[i]]]))
#   aa_remap_track_list[[i]] <- AnnotationTrack(range = ReMapChIP.GRanges.list_merged[[aaint[i]]],
#                                               name= aaint[i],
#                                               genome = "hg38")
# }
# names(aa_remap_track_list) <- aaint

MCF7_YBX1_encode <- readPeakFile("~/Documents/Shayan/BioInf/EstrogenReceptor/Chip-Data/MCF7_YBX1_Encode.bed", as = "GRanges")
aa_ybx1_track <- AnnotationTrack(MCF7_YBX1_encode,
                                 name= "YBX1", 
                                 genome = "hg38")
ReMapChIP.GRanges.list[[86]] <- MCF7_YBX1_encode
names(ReMapChIP.GRanges.list)[86] <- "YBX1"


MCF7_PBX1_encode <- readPeakFile("~/Documents/Shayan/BioInf/EstrogenReceptor/Chip-Data/MCF7_PBX1_noTreat_GSM692743_PBX1_MCF7_red_m16_10_5.bed_peaks.bed", 
                                 as = "GRanges")
MCF7_PBX1_encode <- MCF7_PBX1_encode[MCF7_PBX1_encode$V5 >= 186.94]

aa_18_38_chain <- import.chain("~/Documents/Shayan/BioInf/Liftover_chain_files/hg18ToHg38.over.chain")

MCF7_PBX1_encode_hg38 <- liftOver(MCF7_PBX1_encode, chain = aa_18_38_chain)
MCF7_PBX1_encode_hg38 <- unlist(MCF7_PBX1_encode_hg38)
ReMapChIP.GRanges.list[[87]] <- MCF7_PBX1_encode_hg38
names(ReMapChIP.GRanges.list)[87] <- "PBX1"

aa_pbx1_track <- AnnotationTrack(MCF7_PBX1_encode_hg38,
                                 name= "PBX1", 
                                 genome = "hg38")

MCF7_SP1_encode <- readPeakFile("~/Documents/Shayan/BioInf/EstrogenReceptor/Chip-Data/GSE92014_ENCFF577EMC_optimal_idr_thresholded_peaks_GRCh38_MCF7_SP1.bed", as = "GRanges")
ReMapChIP.GRanges.list[[88]] <- MCF7_SP1_encode
names(ReMapChIP.GRanges.list)[88] <- "SP1"


RARa_ERa_Binding_hg18 <- read.delim("~/Documents/Shayan/BioInf/EstrogenReceptor/Chip-Data/ERa_RARa_coBinding.csv",header = T, sep = ",")
RARa_ERa_Binding_hg18 <- makeGRangesFromDataFrame(RARa_ERa_Binding_hg18)
aa_18_38_chain <- import.chain("~/Documents/Shayan/BioInf/Liftover_chain_files/hg18ToHg38.over.chain")
MCF7_RARa_ERa_encode_hg38 <- liftOver(RARa_ERa_Binding_hg18, chain = aa_18_38_chain)
MCF7_RARa_ERa_encode_hg38 <- unlist(MCF7_RARa_ERa_encode_hg38)
ReMapChIP.GRanges.list[[89]] <- MCF7_RARa_ERa_encode_hg38
names(ReMapChIP.GRanges.list)[89] <- "RARA"
############# 
# for all remaps write the intersect with pos neg set 4 in a different bed file and use only those to create the tracks
ReMapChIP.GRanges.list

aa_pos_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/pos4_1000.bed"
aa_neg_Add <- "~/Documents/Shayan/BioInf/EstrogenReceptor/eqtl/neg4_1000.bed"
aa_pos <- readPeakFile(aa_pos_Add, as = "GRanges")
aa_neg <- readPeakFile(aa_neg_Add, as = "GRanges")
start(aa_pos) <- start(aa_pos) - 1
start(aa_neg) <- start(aa_neg) - 1
names(aa_pos) <- paste0("pos_", c(1:length(aa_pos)))
names(aa_neg) <- paste0("neg_", c(1:length(aa_neg)))
aa_pos_neg <- c(aa_pos, aa_neg)
aa_pos_neg$V7 <- names(aa_pos_neg)


ReMapChIP.GRanges.list_overlapping <- list()
for(i in 1:length(ReMapChIP.GRanges.list)){
  print(i)
  aa_cur_ov <- findOverlaps(query = ReMapChIP.GRanges.list[[i]], subject = aa_pos_neg)
  ReMapChIP.GRanges.list_overlapping[[i]] <-  ReMapChIP.GRanges.list[[i]][aa_cur_ov@from]
}
names(ReMapChIP.GRanges.list_overlapping)  <- names(ReMapChIP.GRanges.list)

ReMapChIP.GRanges.list_overlapping[[52]] <- li_ER_union_gr
aa_remap_track_list <- list()
for(i in 1:length(ReMapChIP.GRanges.list_overlapping)){
  print(i)
  aa_remap_track_list[[i]] <- AnnotationTrack(range = ReMapChIP.GRanges.list_overlapping[[i]],
                                              name= names(ReMapChIP.GRanges.list_overlapping)[i],
                                              genome = "hg38",
                                              stacking = "dense")
}
names(aa_remap_track_list) <-  names(ReMapChIP.GRanges.list_overlapping)







############# eRNA
li_H3K27ac <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/H3K27AC/Li2013/33099_H3K27AC_CHIPSEQ_E2_LI2013_filtered_5_100bpmerged.bed"
li_ER_union <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/ESR1/Li2013/merged_ER_e2_1h_union_rep1_rep2_filtered_5_100bpmerged.bed"  # 34057
li_erna_intersect <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/e_RNA_rep1_rep2_intersecting_reads.bed"
aa_tmp_erna_only_afterE2 <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/Positive_Negative_sets/4/tmp_instersect_eRNA_OnlyAfterE2.bed"
aa_non_E2_eRNA <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Li2013/intersect_ETOH1h_rep1_rep2.bed"

li_H3K27ac_gr <- readPeakFile(li_H3K27ac, as = "GRanges")
li_ER_union_gr <- readPeakFile(li_ER_union, as = "GRanges")
li_erna_only_gr <- readPeakFile(aa_tmp_erna_only_afterE2, as = "GRanges")
aa_li_none2_Erna_gr <- readPeakFile(aa_non_E2_eRNA, as = "GRanges")


aa_tmp_erna_only_afterE2_AnnotationTrack <- AnnotationTrack(li_erna_only_gr,
                                                            name= "eRNA after E2",
                                                            genome = "hg38", 
                                                            stacking = "dense")
aa_tmp_erna_beforeE2_AnnotationTrack <- AnnotationTrack(aa_li_none2_Erna_gr,
                                                            name= "eRNA before E2",
                                                            genome = "hg38", 
                                                            stacking = "dense")

i <- 1
plotTracks(list(aitrack, 
                agtrack, 
                aaot,
                aasnp_track,
                aaot_change, 
                #aaexpression_change_mean_track,
                aaexpression_change_median_track,
                #aaexpression_change_max_track,
                aaexpression_change_percent_track, 
                #aaexpression_change_boxplot_track
                aa_tmp_erna_only_afterE2_AnnotationTrack,
                aa_tmp_erna_beforeE2_AnnotationTrack
), 
from = start(aa_pos_neg)[i], 
to = end(aa_pos_neg[i]),
chromosome = as.character(seqnames(aa_pos_neg[i])),
# ylim = aaylims,
type = c("hist")
#, window = 30
)
# adding H3k27 and ER chip
Exp18_used_TF_names <- aa_names

aa_ER_AnnotationTrack <- AnnotationTrack(li_ER_union_gr,
                                                            name= "ER",
                                                            genome = "hg38", 
                                                            stacking = "dense")
aa_H3k27ac_AnnotationTrack <- AnnotationTrack(li_H3K27ac_gr,
                                                        name= "H3k27Ac",
                                                        genome = "hg38", 
                                                        stacking = "dense")


which(names(aa_pos_neg) == "neg_660")
which(names(aa_pos_neg) == "pos_61")
which(names(aa_pos_neg) == "neg_767")
which(names(aa_pos_neg) == "neg_400")
i <- 849
plotTracks(list(aitrack, 
                agtrack, 
                aaot,
                aasnp_track,
                aaot_change, 
                #aaexpression_change_mean_track,
                aaexpression_change_median_track,
                #aaexpression_change_max_track,
                aaexpression_change_percent_track, 
                #aaexpression_change_boxplot_track,
                aa_tmp_erna_beforeE2_AnnotationTrack,
                aa_tmp_erna_only_afterE2_AnnotationTrack,
                aa_ER_AnnotationTrack,
                aa_H3k27ac_AnnotationTrack,
                aa_ybx1_track,
                aa_pbx1_track
), 
from = start(aa_pos_neg)[i], 
to = end(aa_pos_neg[i]),
chromosome = as.character(seqnames(aa_pos_neg[i]))
# ylim = aaylims,
#type = c("hist")
#, window = 30
)
################################################################################################################
# list some interesting enhancers and plot them
# plot the chip peaks for ER, e2 before and after, and the responsible TF





################################################################################################################
################################################################################################################
####################################         GViz examples         ############################################
################################################################################################################
################################################################################################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Gviz")
library(Gviz)

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
data(cpgIslands)
class(cpgIslands)

aachr <- as.character(unique(seqnames(cpgIslands)))
aagen<- genome(cpgIslands)

atrack <- AnnotationTrack(cpgIslands, name = "CpG")
plotTracks(atrack)

agtrack <- GenomeAxisTrack()
plotTracks(list(agtrack, atrack))

aitrack <- IdeogramTrack(genome = aagen, chromosome = aachr)
plotTracks(list(aitrack, agtrack, atrack))



plotTracks(list(itrack, gtrack, atrack, grtrack, strack), 
           from = 26591822, to = 26591852, cex = 0.8)


set.seed(255)
lim <- c(26700000, 26750000)
coords <- sort(c(lim[1], 
                 sample(seq(from = lim[1], to = lim[2]), 99), 
                 lim[2]))
dat <- runif(100, min = -10, max = 10)
dtrack <- DataTrack(data = dat, start = coords[-length(coords)],
                    end = coords[-1], chromosome = chr, genome = gen, 
                    name = "Uniform")
plotTracks(list(itrack, gtrack, atrack, grtrack, dtrack), 
           from = lim[1], to = lim[2])


WT_binding_bed_df
aachr <- as.character(unique(WT_binding_bed_df$ESR1_2$chr))
aastrack <- SequenceTrack(Hsapiens, chromosome = aachr)

aadtrack <- DataTrack(data = WT_binding_bed_df$ESR1_2$score, start = WT_binding_bed_df$ESR1_2$start,
                      end = WT_binding_bed_df$ESR1_2$end, chromosome = WT_binding_bed_df$ESR1_2$chr, genome = "hg38", 
                      name = "ER binding affinity")
aitrack <- IdeogramTrack(genome = "hg38", chromosome = aachr)
plotTracks(list(aitrack, agtrack, atrack))

agtrack <- GenomeAxisTrack()




dp <- availableDisplayPars(grtrack)
tail(dp)

plotTracks(list(itrack, gtrack, atrack, grtrack), 
           background.panel = "#FFFEDB", background.title = "darkblue")

plotTracks(dTrack, groups = rep(c("control", "treated"), each = 3),
           type = c("a", "p", "confint"))
plotTracks(dTrack, groups = rep(c("control", "treated"), each = 3),
           type = c("a", "p"), legend = TRUE)


plotTracks(dTrack.big, type = "hist", window = -1, windowSize = 2500)



aTrack <- AnnotationTrack(start = c(10, 40, 120), width = 15, 
                          chromosome = "chrX", 
                          strand = c("+", "*", "-"),
                          id = c("Huey", "Dewey", "Louie"), 
                          genome = "hg19", name = "foo")
plotTracks(aTrack)


plotTracks(aTrack, shape = "ellipse", featureAnnotation = "id", 
           fontcolor.feature = "darkblue")

aTrack.groups <- AnnotationTrack(start = c(50, 180, 260, 460, 860, 1240), 
                                 width = c(15, 20, 40, 100, 200, 20),
                                 chromosome = "chrX",
                                 strand = rep(c("+", "*", "-"), 
                                              c(1, 3, 2)),
                                 group = rep(c("Huey", "Dewey", "Louie"), 
                                             c(1, 3, 2)),
                                 genome = "hg19", name = "foo")
plotTracks(aTrack.groups, groupAnnotation = "group")

plotTracks(aTrack.groups, groupAnnotation = "group", 
           just.group = "right")

#overlay
dat <- runif(100, min = -2, max = 22)
dtrack2 <- DataTrack(data = dat, start = coords[-length(coords)], 
                     end = coords[-1], chromosome = chr, genome = gen, 
                     name = "Uniform2",  groups = factor("sample 2", 
                                                         levels = c("sample 1", "sample 2")), legend = TRUE)
displayPars(dtrack) <- list(groups = factor("sample 1", 
                                            levels = c("sample 1", "sample 2")), legend = TRUE)
ot <- OverlayTrack(trackList=list(dtrack2, dtrack))
ylims <- extendrange(range(c(values(dtrack), values(dtrack2))))
plotTracks(list(itrack, gtrack, ot), from = lim[1], to = lim[2], 
           ylim = ylims, type = c("smooth", "p"))


grtrack <- GeneRegionTrack(geneModels, genome = gen, chromosome = chr, 
                           name = "Gene Model",
                           transcriptAnnotation = "symbol",
                           background.title = "brown")
head(displayPars(grtrack))
displayPars(grtrack) <- list(background.panel = "#FFFEDB", col = NULL)
head(displayPars(grtrack))


getOption("Gviz.scheme")
## [1] "default"
scheme <- getScheme()
scheme$GeneRegionTrack$fill <- "salmon"
scheme$GeneRegionTrack$col <- NULL
scheme$GeneRegionTrack$transcriptAnnotation <- "transcript"
addScheme(scheme, "myScheme")
options(Gviz.scheme = "myScheme")
grtrack <- GeneRegionTrack(geneModels, genome = gen,
                           chromosome = chr, name = "Gene Model")
plotTracks(grtrack)



################################################################################################################
################################################################################################################
####################################         SUSHI examples         ############################################
################################################################################################################
################################################################################################################

plotBedgraph(Sushi_DNaseI.bedgraph,chrom = "chr11",chromstart = 1955000,chromend = 1960000,colorbycol= SushiColors(5))
labelgenome(chrom = "chr11",chromstart = 1955000,chromend= 1960000 ,n=4,scale="Mb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)
plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph,chrom = "chr11",chromstart= 1955000,chromend= 1960000,
             transparency=.50,color=SushiColors(2)(2)[1])
plotBedgraph(Sushi_DNaseI.bedgraph,chrom = "chr11",chromstart= 1955000,chromend= 1960000,
             transparency=.50,color=SushiColors(2)(2)[2],overlay=TRUE,
             rescaleoverlay=TRUE)
labelgenome(chrom = "chr11",chromstart= 1955000,chromend= 1960000,n=3,scale="Kb")


i <- 4

plotBed(beddata =aax$ESR1_2,
        chrom = as.character(seqnames(aa_pos_neg[i])), 
        chromstart = start(aa_pos_neg)[i], 
        chromend = end(aa_pos_neg)[i],
        colorby = aax$ESR1_2$strand,
        colorbycol = SushiColors(2)
        ,row  = "auto",
        wiggle=0.001
        #       ,splitstrand=TRUE
)
labelgenome(chrom = as.character(seqnames(aa_pos_neg[i])), 
            chromstart = start(aa_pos_neg)[i],
            chromend = end(aa_pos_neg)[i] 
            ,n=20,scale="bp")



plotBed(
  beddata = Sushi_ChIPSeq_pol2.bed,chrom = "chr11",chromstart = 2281200,
  chromend =2282200,colorby    = Sushi_ChIPSeq_pol2.bed$strand,
  colorbycol = SushiColors(2),row  = "auto",wiggle=0.001,splitstrand=TRUE)
labelgenome(chrom= "chr11",chromstart = 2281200,chromend=2282200,n=2,scale="Kb")


i <-5
plotBedgraph(aax$ESR1_2,
             chrom = as.character(seqnames(aa_pos_neg[i])), 
             chromstart = start(aa_pos_neg)[i], 
             chromend = end(aa_pos_neg)[i]
             #,colorbycol= SushiColors(5)
)
labelgenome(chrom = as.character(seqnames(aa_pos_neg[i])), 
            chromstart = start(aa_pos_neg)[i],
            chromend = end(aa_pos_neg)[i] 
            ,n=100,scale="bp", las = 2)

mtext("site strength",side=2,line=1.75,cex=1,font=2)


################################################################################################################
################################################################################################################

# find interesting enhancers to visualize


ER_SNPs_merged_allbypos_modified <- cbind(ER_SNPs_merged_allbypos_modified, ER_SNPs_merged_allbypos_modified_and_modelScores[, c(283:286)])


ER_SNPs_merged_allbypos_modified$is_pancan
ER_SNPs_merged_allbypos_modified$Gene_cor_pancan
ER_SNPs_merged_allbypos_modified$Enh_gene20kb_ovl
ER_SNPs_merged_allbypos_modified$Enh_geneProm_int
ER_SNPs_merged_allbypos_modified$Enh_gene_both

aa_all_symbol <- unique(unlist(ER_SNPs_merged_allbypos_modified$Gene_cor_pancan))
aa_all_symbol <- aa_all_symbol[!is.na(aa_all_symbol)]


aa_all_symbol_ent <- geneNameToEntrez(aa_all_symbol)
sum(duplicated(aa_all_symbol_ent$ENTREZID))

aaw <- which(ER_SNPs_merged_allbypos_modified$is_pancan)

aaentr <- list()
for(i in 1:length(aaw)){
  aaaal <- unlist(ER_SNPs_merged_allbypos_modified$Gene_cor_pancan[aaw[i]])
  names(aaaal) <- NULL
  aaentr[[i]] <- aa_all_symbol_ent$ENTREZID[aa_all_symbol_ent$ALIAS %in% aaaal]
}

aa_all_pancan_cor_Gene_entrez <- list()
for(i in 1:nrow(ER_SNPs_merged_allbypos_modified)){
  if(i %in% aaw){
    aac <- which(aaw == i)
    aa_all_pancan_cor_Gene_entrez[[i]] <- aaentr[[aac]]
  }else{
    aa_all_pancan_cor_Gene_entrez[[i]] <- NA
  }
  
}
ER_SNPs_merged_allbypos_modified$Gene_cor_pancan_entrez <- aa_all_pancan_cor_Gene_entrez


aaw <- which(ER_SNPs_merged_allbypos_modified$is_pancan)
aa_is_instc <- list()
for(i in 1:length(aaw)){
  aa_cor <- unlist(ER_SNPs_merged_allbypos_modified$Gene_cor_pancan_entrez[[aaw[i]]])
  names(aa_cor) <- NULL
  aa_int <- unlist(ER_SNPs_merged_allbypos_modified$Enh_geneProm_int[[aaw[i]]])
  names(aa_int) <- NULL
  aa_is_instc[[i]] <- intersect(aa_cor,aa_int)
}
table(unlist(lapply(aa_is_instc, length)))

ER_SNPs_merged_allbypos_modified$Enh_Enh_int


aaw <- which(ER_SNPs_merged_allbypos_modified$is_gtex)
ER_SNPs_merged_allbypos_modified$Gene_cor_Gtex[aaw[1]]









