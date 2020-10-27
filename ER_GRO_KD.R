library(ChIPseeker)
# knockdown GRO-seq data

# write the GRanges of pos and neg sequences with their names
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

options(scipen=999)
write.table(aa_pos_neg, 
              file="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Liu2014/pos_neg_4.bed",
              quote=F, sep="\t", row.names=F, col.names=F)

options(scipen=0)

###########
# Find out the pecentile changes after kd of GATA3, RARA, TFAP2C

aaseq <- ShortRead::readFasta("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/Variant_seq/seq_WT.fa")
aaseq2 <- as.character(aaseq@id) #  these are the rownames identifying each enhancer

aa_KD_gata3 <- matrix(nrow = length(aaseq2),
                      ncol = length(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list))
rownames(aa_KD_gata3) <- aaseq2
colnames(aa_KD_gata3) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list)


for(i in 1:ncol(aa_KD_gata3)){
  aa_cur_1 <- Exp18_KD_percentile_change_list[[i]][, 3]
  aa_cur_2 <- Exp18_KD_percentile_change_list[[i]][, 19]
  aass <- apply(X = cbind(aa_cur_1, aa_cur_2), MARGIN = 1, FUN = function(x) x[which.max(abs(x))])
  aa_KD_gata3[, i] <- aass
}


aa_KD_RARA <- matrix(nrow = length(aaseq2),
                      ncol = length(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list))
rownames(aa_KD_RARA) <- aaseq2
colnames(aa_KD_RARA) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list)


for(i in 1:ncol(aa_KD_RARA)){
  aa_cur_1 <- Exp18_KD_percentile_change_list[[i]][, 12]
  aa_cur_2 <- Exp18_KD_percentile_change_list[[i]][, 22]
  aass <- apply(X = cbind(aa_cur_1, aa_cur_2), MARGIN = 1, FUN = function(x) x[which.max(abs(x))])
  aa_KD_RARA[, i] <- aass
}

aa_KD_TFAP2C <- matrix(nrow = length(aaseq2),
                     ncol = length(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list))
rownames(aa_KD_TFAP2C) <- aaseq2
colnames(aa_KD_TFAP2C) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list)


for(i in 1:ncol(aa_KD_TFAP2C)){
  aa_cur_1 <- Exp18_KD_percentile_change_list[[i]][, 16]
  aa_cur_2 <- Exp18_KD_percentile_change_list[[i]][, 23]
  aass <- apply(X = cbind(aa_cur_1, aa_cur_2), MARGIN = 1, FUN = function(x) x[which.max(abs(x))])
  aa_KD_TFAP2C[, i] <- aass
}


aa_KD_gata3
aa_KD_RARA
aa_KD_TFAP2C

aa_KD_gata3_mean <- apply(aa_KD_gata3, 1, mean)
aa_KD_RARA_mean <- apply(aa_KD_RARA, 1, mean)
aa_KD_TFAP2C_mean <- apply(aa_KD_TFAP2C, 1, mean)
###############
# create a TF-enhancer network through KD data
aaseq <- ShortRead::readFasta("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_18/Variant_seq/seq_WT.fa")
aaseq2 <- as.character(aaseq@id) #  these are the rownames identifying each enhancer



TF_enhancer_effect_list <- list()
for(aa_cur_feat in 1:ncol(Exp18_KD_percentile_change_list[[1]])){
  print(aa_cur_feat)
  TF_enhancer_effect_list[[aa_cur_feat]] <- matrix(nrow = length(aaseq2),
                                                   ncol = length(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list))
  rownames(TF_enhancer_effect_list[[aa_cur_feat]] ) <- aaseq2
  colnames(TF_enhancer_effect_list[[aa_cur_feat]] ) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list)
  for(aa_cur_model in 1:ncol(TF_enhancer_effect_list[[aa_cur_feat]])){
    TF_enhancer_effect_list[[aa_cur_feat]][,aa_cur_model] <- Exp18_KD_percentile_change_list[[aa_cur_model]][, aa_cur_feat]
  }
  
}
names(TF_enhancer_effect_list) <- colnames(Exp18_KD_percentile_change_list[[1]])
TF_enhancer_effect_list_mean <- do.call(cbind, lapply(TF_enhancer_effect_list, function(x) apply(x, 1, mean)))
colSums(abs(TF_enhancer_effect_list_mean) > 0.05)
ER_KD_TF_enh_network <- data.frame(TF = character(0),
                                   Enhancer = character(0))
aa_thr <- 0.05
for(i in 1:17){
  aa_cur_enh <- which(abs(TF_enhancer_effect_list_mean[, i]) >= aa_thr)
  if(length(aa_cur_enh) > 0){
    aadf <- data.frame(TF = rep(colnames(TF_enhancer_effect_list_mean)[i], length(aa_cur_enh)),
                       Enhancer = rownames(TF_enhancer_effect_list_mean)[aa_cur_enh])
    ER_KD_TF_enh_network <- rbind(ER_KD_TF_enh_network, aadf)
  }

}
ER_KD_TF_enh_network$TF <- levels(ER_KD_TF_enh_network$TF)[as.numeric(ER_KD_TF_enh_network$TF)]
ER_KD_TF_enh_network$Enhancer <- levels(ER_KD_TF_enh_network$Enhancer)[as.numeric(ER_KD_TF_enh_network$Enhancer)]

head(ER_KD_TF_enh_network)
aacoopkd <- rownames(TF_enhancer_effect_list_mean)[which(TF_enhancer_effect_list_mean[,23] < -0.05)]
aapadp2ctkd
intersect(aacoopkd, aapadp2ctkd)


aa_tfap2c_interesting_Enhs <- intersect(aacoopkd, aapadp2ctkd)

TFAP2C_coop_network <-  ER_KD_TF_enh_network[ER_KD_TF_enh_network$Enhancer %in% aa_tfap2c_interesting_Enhs,]


###############
# read the intersect with ctrl and # read the intersect with GATA3 kd

aa_sictrl_e2_Add <- "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Liu2014/pos_neg_sictrl_E2_intersect.bed"
aa_sigata3_e2_Add <- "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Liu2014/pos_neg_siGATA3_E2_intersect.bed"
aa_sictrl_e2 <- readPeakFile(aa_sictrl_e2_Add, as = "GRanges")
aa_sigata3_e2 <- readPeakFile(aa_sigata3_e2_Add, as = "GRanges")

# read the intersect with ctrl and # read the intersect with RARs kd

aa_shctrl_1_e2_Add <- "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Liu2014/pos_neg_shctrl_E2_intersect.bed"
aa_shRARs_e2_Add <- "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Liu2014/pos_neg_shRARs_E2_intersect.bed"
aa_shctrl_1_e2 <- readPeakFile(aa_shctrl_1_e2_Add, as = "GRanges")
aa_shRARs_e2 <- readPeakFile(aa_shRARs_e2_Add, as = "GRanges")

# read the intersect with ctrl and # read the intersect with AP2g kd

aa_shctrl_2_e2_Add <- "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Liu2014/pos_neg_shctrl_2_E2_intersect.bed"
aa_shAP2g_e2_Add <- "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GROseq/Liu2014/pos_neg_shAP2g_E2_intersect.bed"
aa_shctrl_2_e2 <- readPeakFile(aa_shctrl_2_e2_Add, as = "GRanges")
aa_shAP2g_e2 <- readPeakFile(aa_shAP2g_e2_Add, as = "GRanges")

aa_sictrl_e2$V9 <- as.character(levels(aa_sictrl_e2$V9)[as.numeric(aa_sictrl_e2$V9)])
aa_sigata3_e2$V9 <- as.character(levels(aa_sigata3_e2$V9)[as.numeric(aa_sigata3_e2$V9)])
aa_shctrl_1_e2$V9 <- as.character(levels(aa_shctrl_1_e2$V9)[as.numeric(aa_shctrl_1_e2$V9)])
aa_shRARs_e2$V9 <- as.character(levels(aa_shRARs_e2$V9)[as.numeric(aa_shRARs_e2$V9)])
aa_shctrl_2_e2$V9 <- as.character(levels(aa_shctrl_2_e2$V9)[as.numeric(aa_shctrl_2_e2$V9)])
aa_shAP2g_e2$V9 <- as.character(levels(aa_shAP2g_e2$V9)[as.numeric(aa_shAP2g_e2$V9)])




aa_gata_down_enh <- setdiff(aa_sictrl_e2$V9, aa_sigata3_e2$V9)
aa_gata_up_enh <- setdiff(aa_sigata3_e2$V9, aa_sictrl_e2$V9)

aa_RAR_down_enh <- setdiff(aa_shctrl_1_e2$V9, aa_shRARs_e2$V9)
aa_RAR_up_enh <- setdiff(aa_shRARs_e2$V9, aa_shctrl_1_e2$V9)

aa_AP2g_down_enh <- setdiff(aa_shctrl_2_e2$V9, aa_shAP2g_e2$V9)
aa_AP2g_up_enh <- setdiff(aa_shAP2g_e2$V9, aa_shctrl_2_e2$V9)

###############


# look to see how many of the ac

aath <- -0.05

sum(aa_KD_gata3_mean < aath)
sum(aa_KD_gata3_mean[aa_gata_down_enh] < aath )
length(aa_gata_down_enh)/length(aa_KD_gata3_mean)
sum(aa_KD_gata3_mean[aa_gata_down_enh] < aath )/sum(aa_KD_gata3_mean < aath)

sum(aa_KD_RARA_mean < aath)
sum(aa_KD_RARA_mean[aa_RAR_down_enh] < aath)
length(aa_RAR_down_enh)/length(aa_KD_RARA_mean)
sum(aa_KD_RARA_mean[aa_RAR_down_enh] < aath )/sum(aa_KD_RARA_mean < aath)



sum(aa_KD_TFAP2C_mean < aath)
sum(aa_KD_TFAP2C_mean[aa_AP2g_down_enh] < aath ,na.rm = T)
length(aa_AP2g_down_enh)/length(aa_KD_TFAP2C_mean)
sum(aa_KD_TFAP2C_mean[aa_AP2g_down_enh] < aath ,na.rm = T)/sum(aa_KD_TFAP2C_mean < aath)

aa_KD_TFAP2C_mean_actdown <- aa_KD_TFAP2C_mean[aa_AP2g_down_enh]
aapadp2ctkd <- names(aa_KD_TFAP2C_mean_actdown)[which(aa_KD_TFAP2C_mean_actdown < -0.05)]
aath <- -0.05
phyper(q = sum(aa_KD_TFAP2C_mean[aa_AP2g_down_enh] < aath ,na.rm = T),
       m = sum(aa_KD_TFAP2C_mean < aath), 
       n = sum(aa_KD_TFAP2C_mean >= aath), 
       k = length(intersect(aa_AP2g_down_enh, names(aa_KD_TFAP2C_mean))),
       lower.tail = F)

aa_shctrl_2_e2$V9
phyper(q = sum(aa_KD_TFAP2C_mean[aa_AP2g_down_enh] < aath ,na.rm = T),
       m = sum(aa_KD_TFAP2C_mean[names(aa_KD_TFAP2C_mean) %in% aa_shctrl_2_e2$V9] < aath), 
       n = sum(aa_KD_TFAP2C_mean[names(aa_KD_TFAP2C_mean) %in% aa_shctrl_2_e2$V9] >= aath), 
       k = length(intersect(aa_AP2g_down_enh, names(aa_KD_TFAP2C_mean))),
       lower.tail = F)

aadfdf <- data.frame(row.names = names(aa_KD_TFAP2C_mean),
                     KD_predicted = aa_KD_TFAP2C_mean < aath,
                     KD_observed = names(aa_KD_TFAP2C_mean) %in% aa_AP2g_down_enh)


png(filename = "~/Desktop/ER_paper/fig5_d.png", 
    width = 4*300, 
    height = 4*300, 
    res = 300,
    pointsize = 8)
VennDiagram::draw.pairwise.venn(area1 = sum(aadfdf$KD_predicted), 
                                area2 = sum(aadfdf$KD_observed), 
                                cross.area = sum(rowSums(aadfdf) == 2),
                                category = c("Predicted", "Observed"),
                                fill = c("blue", "red"),
                                cex = 2,
                                cat.cex = 2,
                                cat.just = list(c(-1, -2), c(0.5, -0.5)),
                                #cat.col = c("blue", "red"),
                                fontfamily = c("Arial"),
                                cat.fontfamily = c("Arial", "Arial"),
                                scaled = T, 
                                margin = 0.1)

dev.off()
aath <- -0.05
phyper(q = sum(aa_KD_TFAP2C_mean[aa_AP2g_down_enh] < aath ,na.rm = T),
       m = sum(aa_KD_TFAP2C_mean < aath), 
       n = sum(aa_KD_TFAP2C_mean >= aath), 
       k = length(intersect(aa_AP2g_down_enh, names(aa_KD_TFAP2C_mean))),
       lower.tail = F)

aadfdf <- data.frame(row.names = names(aa_KD_TFAP2C_mean),
                     KD_predicted = aa_KD_TFAP2C_mean < aath,
                     KD_observed = names(aa_KD_TFAP2C_mean) %in% aa_AP2g_down_enh)

png(filename = "~/Desktop/ER_paper/fig5_d_2_a.png", 
    width = 4*300, 
    height = 4*300, 
    res = 300,
    pointsize = 8)
VennDiagram::draw.pairwise.venn(area1 = sum(aadfdf$KD_predicted), 
                                area2 = sum(aadfdf$KD_observed), 
                                cross.area = sum(rowSums(aadfdf) == 2),
                                category = c("Predicted", "Observed"),
                                fill = c("blue", "red"),
                                cex = 2,
                                cat.cex = 2,
                                cat.just = list(c(1, -1.5), c(0, -1.5)),
                                #cat.col = c("blue", "red"),
                                fontfamily = c("Arial"),
                                cat.fontfamily = c("Arial", "Arial"),
                                scaled = T, 
                                margin = 0.1)

dev.off()

aadfdf <- data.frame(row.names = names(aa_KD_TFAP2C_mean)[names(aa_KD_TFAP2C_mean) %in% aa_shctrl_2_e2$V9],
                     KD_predicted = aa_KD_TFAP2C_mean[names(aa_KD_TFAP2C_mean) %in% aa_shctrl_2_e2$V9] < aath,
                     KD_observed = names(aa_KD_TFAP2C_mean)[names(aa_KD_TFAP2C_mean) %in% aa_shctrl_2_e2$V9] %in% aa_AP2g_down_enh)

png(filename = "~/Desktop/ER_paper/fig5_d_2_a_fixed.png", 
    width = 4*300, 
    height = 4*300, 
    res = 300,
    pointsize = 8)
VennDiagram::draw.pairwise.venn(area1 = sum(aadfdf$KD_predicted), 
                                area2 = sum(aadfdf$KD_observed), 
                                cross.area = sum(rowSums(aadfdf) == 2),
                                category = c("Predicted", "Observed"),
                                fill = c("blue", "red"),
                                cex = 2,
                                cat.cex = 2,
                                cat.just = list(c(1, -1.5), c(0, -1.5)),
                                #cat.col = c("blue", "red"),
                                fontfamily = c("Arial"),
                                cat.fontfamily = c("Arial", "Arial"),
                                scaled = T, 
                                margin = 0.1)

dev.off()
#########################################
# patient eRNA data, read and see if you can define a binary vector for each patient, with one entry
#  per enhancer that shows the activity of the enhancer in that patient.
# Then define various states, for each TF KD, TF overexpression and then see if those TFs/(maybe some enhancers can stratify) and stratify the patients
# compare patient profile with TF KD, OE, using hamming distance.



aa_pos_neg["neg_383"]

aa_pos_neg["pos_118"]

aa_pos_neg["pos_197"]




