# creating ITF jobs
# working with ITF results
##########################################################################################
# write inputs to perform iTF on the positive set for all pairs of TFs
# write motifs:
MotifWriter(motif.List = TF.motifs.Expanded_pseudo_count,  pseudo = 0.001,
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/motifs_expanded_pseudo_cnt")
# write fasta of sequences positive

# write fasta of sequences negative
length(Negative_set_seq_list_char_1000bp[[4]])
length(Positive_set_seq_list_char_1000bp[[4]])
writeFasta(DNAStringSet(Negative_set_seq_list_char_1000bp[[4]]), width = 1000,
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/negative_4.fa")
writeFasta(DNAStringSet(Positive_set_seq_list_char_1000bp[[4]]), width = 1000, 
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/positive_4.fa")
# meta for positive
# sample_inputs/trl.fa	TRL	TRL	0-15
aa_al <- character(0)
for(i in c(seq(5, 50, 5), seq(60, 100, 10), seq(120, 200, 20), seq(250, 500, 50), seq(600, 1000, 200))){
  aa_al <- c(aa_al, paste0("0-", i))
}
#aa_al <- paste(aa_al, collapse = ", ")

for(i in 1:length(TF.motifs.Expanded_pseudo_count)){
  for(j in i:length(TF.motifs.Expanded_pseudo_count)){
    cat(c("Inputs/positive_4.fa",
          names(TF.motifs.Expanded_pseudo_count)[i], 
          names(TF.motifs.Expanded_pseudo_count)[j],
          paste0(aa_al, "\n")), 
        sep = "\t",
        file = "positive_meta_1.meta",
        append = T)
  }
}
# meta for negative

for(i in 1:length(TF.motifs.Expanded_pseudo_count)){
  for(j in i:length(TF.motifs.Expanded_pseudo_count)){
    cat(c("Inputs/negative_4.fa",
          names(TF.motifs.Expanded_pseudo_count)[i], 
          names(TF.motifs.Expanded_pseudo_count)[j],
          paste0(aa_al, "\n")), 
        sep = "\t",
        file = "negative_meta_1.meta",
        append = T)
  }
}

# jobs
# ./wrap_iTFs.pl Inputs/motifs_expanded_pseudo_cnt.wtmx Inputs/positive_meta_1.meta outdir_pos_1 -thresh 4,5,6

# ./wrap_iTFs.pl Inputs/motifs_expanded_pseudo_cnt.wtmx Inputs/negative_meta_1.meta outdir_neg_1 -thresh 4,5,6

###############
# reading ITF results
aa <- read.table("~/Documents/iTFs_Mac/outdir_pos_1/summary/AR_ESR1_3.d0-5.txt", 
                 stringsAsFactors = F, header = T)
# shit
# writing meta files to run on hal
aa_al <- character(0)
for(i in c(seq(5, 50, 5), seq(60, 100, 10), seq(120, 200, 20),
           seq(250, 500, 50), seq(600, 1000, 200))){
  aa_al <- c(aa_al, paste0("0-", i))
}
aa_al <- paste(aa_al, collapse = ",")
setwd("~/Documents/iTFs_Mac/Inputs_hal/")
for(i in 1:length(TF.motifs.Expanded_pseudo_count)){
  for(j in i:length(TF.motifs.Expanded_pseudo_count)){
    cat(c("Input/positive_4.fa",
          names(TF.motifs.Expanded_pseudo_count)[i], 
          names(TF.motifs.Expanded_pseudo_count)[j],
          aa_al), 
        sep = "\t",
        file = paste0("positive_meta_",i, "_", j,".meta" ),
        append = F)
  }
}
# meta for negative

for(i in 1:length(TF.motifs.Expanded_pseudo_count)){
  for(j in i:length(TF.motifs.Expanded_pseudo_count)){
    cat(c("Input/negative_4.fa",
          names(TF.motifs.Expanded_pseudo_count)[i], 
          names(TF.motifs.Expanded_pseudo_count)[j],
          aa_al), 
        sep = "\t",
        file = paste0("negative_meta_",i, "_", j,".meta" ),
        append = F)
  }
}


# write job files for running on hal
#./wrap_iTFs.pl Inputs/motifs_expanded_pseudo_cnt.wtmx Inputs/negative_meta_1.meta outdir_neg_1 -thresh 4,5,6
aaf_p <- list.files(path = ".", full.names = F, pattern = "positive_*")
aaf_n <- list.files(path = ".", full.names = F, pattern = "negative_*")

for(i in 1:length(aaf_p)){
  cat(c(paste0("./wrap_iTFs.pl Input/motifs_expanded_pseudo_cnt.wtmx Input/", aaf_p[i]),
        "outdir_pos_hal_1 -thresh 4,5,6\n"), sep = " ", append = T, file = "hal_iTF.job")
}
for(i in 1:length(aaf_n)){
  cat(c(paste0("./wrap_iTFs.pl Input/motifs_expanded_pseudo_cnt.wtmx Input/", aaf_n[i]),
        "outdir_neg_hal_1 -thresh 4,5,6\n"), sep = " ", append = T, file = "hal_iTF.job")
}
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/")

###############
# read TOMTOM results
aatomtom <- read.table("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/tomtom.csv", 
                       sep = ",", header = T, stringsAsFactors = F)
aatomtom[1,1] <- "Query_ID"
names(aatomtom) <- aatomtom[1,]
aatomtom <- aatomtom[-c(1),]
TF.motifs.Expanded_pseudo_count_TOMTOM_res <- aatomtom
aacolsnum <- colnames(TF.motifs.Expanded_pseudo_count_TOMTOM_res)[c(3,4,5,6,7)]
TF.motifs.Expanded_pseudo_count_TOMTOM_res[aacolsnum] <- sapply(TF.motifs.Expanded_pseudo_count_TOMTOM_res[aacolsnum], as.numeric)
TF.motifs.Expanded_pseudo_count_TOMTOM_res <- TF.motifs.Expanded_pseudo_count_TOMTOM_res[1:1473,]
aatomtom <- TF.motifs.Expanded_pseudo_count_TOMTOM_res
sum(aatomtom$`q-value` <= 0.2)


###########
## reading ITF results second round
aa_fnames_full <- list.files("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/ITF_Results/summary", full.names = T)
aa_fnames <- list.files("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/ITF_Results/summary", full.names = F)


aa1 <- read.table(file = aa_fnames_full[1], header = T, stringsAsFactors = F)
View(aa1)
aggregate(x = aa1[c("var_name", "varPos_end", "bcPos_st", "bc_seq20bp")],
          by = aa1[c("lib_name")], 
          FUN = c)
# read the results in this format:
#  create a list with one entry for each threshold
#  each entry is a list containing 561 entries for each pair
#  each entry of the previous list is a matrix (9*29): rows corresponding to, distance p-value, orientation-specific distance p-values, orientation p-values and 
#   each column correspondng to a distance. entries are negative log pvalues
Positive_set_ITF_results <- list()
for(i in 1:3){
  Positive_set_ITF_results[[i]] <- list()
}
names(Positive_set_ITF_results) <- c(4,5,6)

aaal <- c(seq(5, 50, 5), seq(60, 100, 10), seq(120, 200, 20),
          seq(250, 500, 50), seq(600, 1000, 200))
for(k in 1:length(Positive_set_ITF_results)){
  aacnt <- 1
  aa_name <- character(0)
  for(i in 1:length(TF.motifs.Expanded_pseudo_count)){
    for(j in i:length(TF.motifs.Expanded_pseudo_count)){
      aa_name <- c(aa_name, paste(names(TF.motifs.Expanded_pseudo_count)[i],
                                  names(TF.motifs.Expanded_pseudo_count)[j], sep = "_"))
      Positive_set_ITF_results[[k]][[aacnt]] <- matrix(nrow = 9, ncol = 29)
      colnames(Positive_set_ITF_results[[k]][[aacnt]]) <- aaal
      rownames(Positive_set_ITF_results[[k]][[aacnt]]) <- c("dist_all", paste0("dist_", c(1:4)), paste0("ori_", c(1:4)))
      aacnt <- aacnt + 1
    }
  }
  names(Positive_set_ITF_results[[k]]) <- aa_name
}

aa_pair_all <- unlist(lapply(strsplit(aa_fnames, split = "\\."), "[[", 1))
aa <- strsplit(unlist(lapply(strsplit(aa_fnames, split = "\\."), "[[", 2)), split = "-")
aa_col_all <- unlist(lapply(aa, "[[", 2))
for(aa_file in 1:length(aa_fnames_full)){
  aa_curfile <- read.table(file = aa_fnames_full[aa_file], 
                           header = T, stringsAsFactors = F)
  aa_pair_name <- aa_pair_all[aa_file]
  aa_col_name <- aa_col_all[aa_file]
  aa_thr_uniq <- unique(aa_curfile$fimo_thr)
  for(aa_thr in 1:length(aa_thr_uniq)){
    aa_cur_tab <- aa_curfile[aa_curfile$fimo_thr == aa_thr_uniq[aa_thr], ]
    for(a_row in 1:nrow(aa_cur_tab)){
      stopifnot(nrow(aa_cur_tab) == 5)
      if(a_row == 1){
        Positive_set_ITF_results[[as.character(aa_thr_uniq[aa_thr])]][[aa_pair_name]][a_row, aa_col_name] <- -log10(aa_cur_tab$cond_dist_pval[a_row])
      }else{
        Positive_set_ITF_results[[as.character(aa_thr_uniq[aa_thr])]][[aa_pair_name]][a_row, aa_col_name] <- -log10(aa_cur_tab$cond_dist_pval[a_row])
        Positive_set_ITF_results[[as.character(aa_thr_uniq[aa_thr])]][[aa_pair_name]][a_row + 4, aa_col_name] <- -log10(aa_cur_tab$ori_pval[a_row])
      }
    }
  }
}
aa_nu_sig_list <- list()
for(i in 1:3){
  aa_nu_sig_list[[i]] <- matrix(nrow = 9, ncol = 561)
  colnames(aa_nu_sig_list[[i]]) <- names(Positive_set_ITF_results[[1]])
  rownames(aa_nu_sig_list[[i]]) <- rownames(Positive_set_ITF_results[[1]][[1]])
  
}
for(i in 1:length(Positive_set_ITF_results)){
  for(j in 1:length(Positive_set_ITF_results[[i]])){
    aa_nu_sig_list[[i]][,j] <- rowSums(Positive_set_ITF_results[[i]][[j]] >= 2)
  }
}
summary(aa_nu_sig_list[[1]][1,])
sum(aa_nu_sig_list[[3]][1:5,] > 0)

which(colSums(aa_nu_sig_list[[3]]) > 0)
aa_nu_sig_list[[3]][,which(colSums(aa_nu_sig_list[[3]]) > 0)]
# how to show:
#  pick one threshold
#  for pairs with significance
#  for each pair : line plot: 5 lines for the "allTOgether + four" orientations (if the orientation reaches significance).
#  y axis: negative log10 pvalue, x axis: distance

# read negative results
load("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/ITF_Results/Negative_set_ITF_results.RData")

aa_nu_sig_list2 <- list()
for(i in 1:3){
  aa_nu_sig_list2[[i]] <- matrix(nrow = 9, ncol = 561)
  colnames(aa_nu_sig_list2[[i]]) <- names(Negative_set_ITF_results[[1]])
  rownames(aa_nu_sig_list2[[i]]) <- rownames(Negative_set_ITF_results[[1]][[1]])
  
}
for(i in 1:length(Negative_set_ITF_results)){
  for(j in 1:length(Negative_set_ITF_results[[i]])){
    aa_nu_sig_list2[[i]][,j] <- rowSums(Negative_set_ITF_results[[i]][[j]] >= 2)
  }
}
sum(is.na(aa_nu_sig_list2[[1]]))
sum(aa_nu_sig_list2[[1]] > 0, na.rm = T)
sum(aa_nu_sig_list2[[2]] > 0, na.rm = T)
sum(aa_nu_sig_list2[[3]] > 0, na.rm = T)
sum(aa_nu_sig_list2[[1]][1:5,] > 0, na.rm = T)
aa_nu_sig_list2[[3]][, which(colSums(aa_nu_sig_list2[[3]][1:5,]) > 0)]

Negative_set_ITF_results$`6`$CEBPB_ESR1_1



aaw <- which(colSums(aa_nu_sig_list2[[1]][1:5,]) > 0)

par(mfrow = c(3,3), mar= c(3,3,1,1))
for(i in 1:length(aaw)){
  plot(Negative_set_ITF_results[[1]][[aaw[i]]][1, ],
       ylim = range(Negative_set_ITF_results$`4`[[aaw[i]]][1:5, ], na.rm = T),
       col = 1, main = names(aaw)[i], ylab = "", xlab = "",
       type = "l", xaxt = "n", lwd = 1.5)
  axis(side = 1, at = c(1:ncol(Negative_set_ITF_results[[1]][[aaw[i]]]))[seq(1, 29, 3)],
       labels = colnames(Negative_set_ITF_results[[1]][[aaw[i]]])[seq(1, 29, 3)]
       , las = 2)
  for(j in 2:5){
    lines(Negative_set_ITF_results[[1]][[aaw[i]]][j,], col = j,lwd = 1.5)
  }
}

par(mfrow = c(3,3), mar= c(3,3,1,1))
for(i in 1:length(aaw)){
  plot(Positive_set_ITF_results[[1]][[aaw[i]]][1, ],
       ylim = range(Positive_set_ITF_results$`4`[[aaw[i]]][1:5, ], na.rm = T),
       col = 1, main = names(aaw)[i], ylab = "", xlab = "", type = "l", xaxt = "n", lwd = 1.5)
  axis(side = 1, at = c(1:ncol(Positive_set_ITF_results[[1]][[aaw[i]]]))[seq(1, 29, 3)],
       labels = colnames(Positive_set_ITF_results[[1]][[aaw[i]]])[seq(1, 29, 3)]
       , las = 2)
  for(j in 2:5){
    lines(Positive_set_ITF_results[[1]][[aaw[i]]][j,], col = j, lwd = 1.5)
  }
}



