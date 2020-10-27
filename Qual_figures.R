#  Qual paper figures
install.packages("Rtsne")
install.packages("plotrix")
install.packages("gapplot")
library(Rtsne)
library(plotrix)
library(gapplot)
library(ggplot2)

# figure 2
# RF input features, performance and important output features

# 2.a)
### sumLR
aa_newmat <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]]),
                               nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]])),
                    ncol = 66)
aat_test <- list()
colnames(aa_newmat) <- character(66)
for(i in 1:ncol(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]])){
  aaln1 <- length(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]][, i])
  aaln2 <- length(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]][, i])
  aa_normal <- scale(c(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]][, i], 
                       Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]][, i]))
  aa_newmat[1:nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]]),(2*i - 1)] <- aa_normal[1:aaln1]
  colnames(aa_newmat)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]])[i], "Pos", sep = "_")
  aa_newmat[1:nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]]),(2*i)] <- aa_normal[(aaln1+1):length(aa_normal)]
  colnames(aa_newmat)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]])[i], "Neg", sep = "_")
  aat_test[[i]] <- t.test(aa_normal[1:aaln1], aa_normal[(aaln1+1):length(aa_normal)])
  print(colnames(aa_newmat)[2*i - 1])
  print(colnames(aa_newmat)[2*i])
  print(aat_test[[i]] )
}

par(mfrow = c(1, 1), mar = c(4,4,4,4))

boxplot.matrix(aa_newmat[, c(7,8 ,9,10, 31,32, 35,36)],
               #las = 2, 
               main = "Likelihood ratio",
               at= c(1,2, 3.5,4.5, 6,7, 8.5,9.5),
               xaxt = "n", 
               ylab = "normalized sum of LR",
               outline = F,
              col = rep(c("green", "red"), 4))
axis(side = 1, at = c(1.5, 4, 6.5, 9), labels = c("ER dimer", "ER monomer", 
                                                        "NR2F2", "NR5A2"))
legend("topleft", legend = c("positive", "negative"), 
       fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, 
       bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)


###2.b adjacency
aa_newmat2 <- matrix(nrow = max(nrow(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]]),
                                nrow(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]])),
                     ncol = sum(ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]]),
                                ncol(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]])))
colnames(aa_newmat2) <- character(sum(ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]]),
                                      ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]])))
aat_test2 <- list()
aatoplot <- numeric(0)
for(i in 1:ncol(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]])){
  aaln1 <- length(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]][, i])
  aaln2 <- length(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]][, i])
  aa_normal <- scale(c(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]][, i], 
                       Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]][, i]))
  aa_newmat2[1:nrow(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]]),(2*i - 1)] <- aa_normal[1:aaln1]
  colnames(aa_newmat2)[2*i - 1] <- paste(colnames(Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]])[i], "Pos", sep = "_")
  aa_newmat2[1:nrow(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]]),(2*i)] <- aa_normal[(aaln1+1):length(aa_normal)]
  colnames(aa_newmat2)[2*i] <- paste(colnames(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]])[i], "Neg", sep = "_")
  aat_test2[[i]] <- t.test(aa_normal[1:aaln1],
                           aa_normal[(aaln1+1):length(aa_normal)])
  if(aat_test2[[i]]$p.value <= 1e-4){
    aatoplot <- c(aatoplot, c((2*i - 1), (2*i)))
    print(2*i - 1)
    print(colnames(aa_newmat2)[2*i - 1])
    print(2*i)
    print(colnames(aa_newmat2)[2*i])
    print(aat_test2[[i]] )
  }

  
}
par(mfrow = c(1, 1), mar = c(7,4,4,4))
boxplot.matrix(aa_newmat2[, c(17,18 , 27,28, 37,38, 49,50, 217,218,
                              247,248, 269,270, 271,272,
                              467,468, 559,560, 583,584, 
                              587,588, 593,594, 1061,1062)],
               boxwex = 0.25,
               notch = F,
               #las = 2, 
               main = "Pairwise adjacency",
               ylab = "normalized adjacency score",
               at= as.numeric(rbind(seq(0, 14, by = 0.75), seq(0.25, 14, by = 0.75)))[1:28],
               xaxt = "n", 
               outline = F,
               col = rep(c("green", "red"), 16))
axis(side = 1, at = colMeans(rbind(seq(0, 14, by = 0.75), 
                                   seq(0.25, 14, by = 0.75)))[1:14], las =2,
     labels = c("ER:ER", "Half-ER:ER", "FOXA1:ER", "GATA3:ER", "NR2C1:ER", 
                "NR2F2:ER", "NR2F2:NR2C1", "NR2F2:NR2F2", "POU5F1:ER","PPARG:ER",
                "PPARG:NR2F2", "PPARG:NR5A2", "PPARG:PGR", "YBX1:ER"))
legend("bottomright", legend = c("positive", "negative"), 
       fill = c("green", "red"), 
       cex = 0.8, inset=.02,  horiz=F, 
       bty = "n", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)

# 2.c random forest performance

aa_gglist_4 <- list()

for(i in c(2)){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  
  aatrainind <- match(c(data_partition_mut_GEMSTAT[[1]], data_partition_mut_GEMSTAT[[2]]),
                      rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]))
  aagridLen <- c(5, 30, 50)
  aa_train_data <- aa_all_data[aatrainind, ]
  aa_test_data <-  aa_all_data[-aatrainind, ]
  
  aaindex_pos <- aa_test_data$label_set == "Pos"
  aaindex_neg <- aa_test_data$label_set == "Neg"
  
  aarfdown <- predict(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[i]], 
                      aa_test_data, type = "prob")
  
  aapreds_list <- list(aarfdown$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_test_data$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = c("RF test ROC\n", paste0("AUROC: " ,
                                                                    format(aauc@y.values[[1]], digits = 2))), lwd = 1.2)
  abline(a = 0, b = 1, col = "gray", lty = 3, lwd = 0.7)
  # legend(x = "bottomright", cex = 0.6,
  #        legend = c(paste0("AUROC: " ,
  #                          format(aauc@y.values[[1]], digits = 2))))
  
  aa_model_list <- list(RF_down = my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[i]])
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_test_data))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  for(the_pr in aamodel_list_pr){
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    num_mod <- num_mod + 1
  }
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  aacustom_col <- "black"
  
  aa_gglist_4[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(color = "black", size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    #scale_color_manual(values = "", labels = "") +
    geom_abline(intercept = sum(aa_test_data$label_set == "Pos")/nrow(aa_test_data),
                slope = 0, color = "gray", size = 1) +
    ggtitle(label = paste0("RF test PRC\n ", paste("AUPRC:" ,
                                                   format(aamodel_list_pr_auc, digits = 2))))+
    theme_light() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text( size = 10),
          axis.text.y = element_text( size = 10))
    
}
names(aa_gglist_4) <-  names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2)
aa_gglist_4[[2]]
print(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[2]])
aa <- varImp(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[2]],scale = F)
aax <- aa$importance
aaxn <- rownames(aax)
aax <- aax$Overall
names(aax) <- aaxn
aax[]
par(mfrow = c(1,1), mar = c(4,10,4,4))

barplot(aax[(sort(aax, decreasing = T, index.return=T)$ix[1:20])], horiz = F, las = 2,
        main = "", xaxt = "n")

# figure 3
E_RNA_GEMSTAT_Ensemble_Outlist[[7]]
E_RNA_GEMSTAT_Ensemble_Parlist[[7]]
par(mfrow = c(1, 2), mar = c(4,4,4,4))
aacol <- integer(length(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$train_ROC))
aacol[] <- 1
aacol[E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_ROC >= 0.60 & E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_PRC >= 0.32] <- 2
# plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_ROC , E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$train_ROC,
#      main = "AUROC", xlab = "test", ylab = "training", pch = 16, cex = 0.8, col = aacol)
# plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_PRC, E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$train_PRC, 
#      main = "AUPRC", xlab = "test", ylab = "training", pch = 16, cex = 0.8, col = aacol)

par(mfrow = c(1, 2), mar = c(4,4,4,4))
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_ROC , E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$train_ROC,
     main = "AUROC", xlab = "validation", ylab = "training", pch = 16, cex = 0.8, col = aacol)
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_PRC, E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$train_PRC, 
     main = "AUPRC", xlab = "validation", ylab = "training", pch = 16, cex = 0.8, col = aacol)

# plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_ROC[aacol == 2] , E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_ROC[aacol == 2],
#      main = "AUROC", xlab = "validation", ylab = "test", pch = 16, cex = 0.8)
# plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_PRC[aacol == 2], E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_PRC[aacol == 2],
#      main = "AUPRC", xlab = "validation", ylab = "test", pch = 16, cex = 0.8)

# plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_PRC[aacol == 2], E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_ROC[aacol == 2],
#      main = "", xlab = "validation AUPRC", ylab = "validation AUROC", pch = 16, cex = 0.8)

par(mfrow = c(1, 1), mar = c(4,4,4,4))
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_PRC[aacol == 2], 
     E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_ROC[aacol == 2], 
     main = "", xlab = "test AUPRC", ylab = "test AUROC", pch = 16, cex = 0.8)


plot(aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$valid_results$output_log_par_420_1_valid$ROC_curve,
     main = "valid ROC")
plot(aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$valid_results$output_log_par_420_1_valid$PRC_curve,
     main = "valid PRC")
plot(aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_results$output_log_par_420_1_test$ROC_curve,
     main = "test ROC")
plot(aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_results$output_log_par_420_1_test$PRC_curve,
     main = "test PRC")




aachosenModel <- which(aacol == 2)
aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered <- E_RNA_GEMSTAT_Ensemble_Outlist[[7]]
aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$train_results <- aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$train_results[aachosenModel]
aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_results <- aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_results[aachosenModel]
aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$valid_results <- aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$valid_results[aachosenModel]

aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_ROC <- aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_ROC[aachosenModel]
aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_PRC <- aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_PRC[aachosenModel]
aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$train_ROC <- aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$train_ROC[aachosenModel]
aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$train_PRC <- aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$train_PRC[aachosenModel]
aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$valid_ROC <- aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$valid_ROC[aachosenModel]
aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$valid_PRC <- aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$valid_PRC[aachosenModel]

aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered <- E_RNA_GEMSTAT_Ensemble_Parlist[[7]]
aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$annot <- aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$annot[aachosenModel, ]
aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$binding <- aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$binding[aachosenModel, ]
aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$alpha <- aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$alpha[aachosenModel, ]
aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$coop <- aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$coop[aachosenModel, ]
aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$qbtm <- aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$qbtm[aachosenModel]
aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$beta <- aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$beta[aachosenModel]
aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$logistic <- aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$logistic[aachosenModel,]

aa_all_param <- cbind(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$binding,
                      aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$alpha,
                      aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$coop,
                      aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$qbtm, 
                      aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$beta,
                      aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$logistic)
colnames(aa_all_param)[1:18] <- paste(colnames(aa_all_param)[1:18], "_binding")
colnames(aa_all_param)[19:36] <- paste(colnames(aa_all_param)[19:36], "_activation")
colnames(aa_all_param)[47] <- "qbtm"
colnames(aa_all_param)[48] <- "beta"
aaparam_pca <- prcomp(aa_all_param[,-c(20)], scale. = F, center = T)
summary(aaparam_pca)
# library(ggbiplot)
# ggbiplot(aaparam_pca,ellipse=F,circle=F, col=aamydf2$test_ROC
#          ,obs.scale = 1, var.scale = 1, var.axes=T
#         # , labels=rownames(mtcars), groups=mtcars.country
#          ) 
#+
  #geom_col(data = aamydf2, aes(color=test_ROC))
aamydf2 <- as.data.frame(cbind(aaparam_pca$x[, 1:2], aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_ROC))
names(aamydf2) <- c("pca1", "pca2", "Test_AUROC")
ggplot(aamydf2, aes(x=pca1, y=pca2, color=Test_AUROC)) +
  geom_point(size=2) + 
  ggtitle(label = "PCA analysis of GEMSTAT parameters") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

aa_perf_srt <-names(aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_PRC)[sort(aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_PRC,
                                                                               decreasing = T, index.return=T)$ix]
aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
aa_perf_srt22 <- paste0("log_", aa_perf_srt2)

# heatmap.2(log10(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$binding[aa_perf_srt22, ]), Rowv = T, Colv = T,
#           dendrogram = 'both', trace = 'none')
# heatmap.2(log10(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$alpha[aa_perf_srt22, ] + 1e-2), Rowv = F, Colv = T,
#           dendrogram = 'column', trace = 'none')
# heatmap.2(log10(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$coop[aa_perf_srt22, ]), Rowv = F, Colv = T,
#           dendrogram = 'column', trace = 'none', margins = c(10,10))

aa_alpar <- cbind(log10(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$binding[aa_perf_srt22, ]), 
                  log10(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$alpha[aa_perf_srt22, ] + 1e-2),
                  log10(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$coop[aa_perf_srt22, ]))
colnames(aa_alpar)[1:18] <- paste(colnames(aa_alpar)[1:18], "_binding")
colnames(aa_alpar)[19:36] <- paste(colnames(aa_alpar)[19:36], "_activation")

heatmap.2(aa_alpar,
          Rowv = T, Colv = T,dendrogram = 'both',
          trace = 'none', margins = c(8,5), scale = "none")

# plot(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$qbtm[aa_perf_srt22],
#      main = "qbtm_sorted_by_perf_exp5_hER", pch = 16, ylab = "")
# plot(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$beta[aa_perf_srt22],
#      main = "beta_sorted_by_perf_exp5_hER", pch = 16)
# 
# plot(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$logistic[aa_perf_srt22, 1], 
#      ylab="", main = "logistic_bias_sorted_by_perf_exp5_hER")
# plot(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$logistic[aa_perf_srt22, 2],
#      ylab="", main = "logistic_coeff_sorted_by_perf_exp5_hER")



# aatsne <- Rtsne(aa_all_param[,-c(20)])
# par(mfrow = c(1, 1), mar = c(4,4,4,4))
# aamydf <- as.data.frame(cbind(aatsne$Y, aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$test_ROC))
# names(aamydf) <- c("tsne1", "tsne2", "test_ROC")
# ggplot(aamydf, aes(x=tsne1, y=tsne2, color=test_ROC)) +
#  geom_point(size=3) + 
#    theme_bw()

# aa1 <- unlist(lapply(lapply( strsplit( names(aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$train_results), split= "_"), "[", c(2,3,4) ), paste, collapse="_"))
# aa2 <- unlist(lapply(lapply( strsplit( rownames(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$annot), split= "_"), "[", c(2,3,4) ), paste, collapse="_"))
# sum(aa1 == aa2)
# a)GEMSTAT input feautures,
# b) GEMSTAT ensemble results,
# c) best models test performance ?? which ones? just one or more?
# d) how many of false positives are actually positive in one replicate, how many of True negatives are present in at least one replicate.

# figure 4
# parameter sensitivity heatmap
# parameter PCA colored by performance
#  parameter values of best performing models?
# TF, coop knock down of best training models,
par(mfrow = c(3, 4), mar = c(2,2,2,2))
for(i in c(3, 5,6,8,9, 11, 12, 13, 18,19,20,21)){
  plot( E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list$output_par_420_1[, 2], 
        E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list$output_par_420_1[, i], 
        #xaxt = "s", yaxt = "s",
        ylab="KD", xlab = "WT", labels =F, tick = F,
        main = colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list$output_par_420_1)[i])
   mtext(text = "KD", side = 2, outer = F, line = 0.5, cex = 0.8)
   mtext(text = "WT", side = 1, outer = F, line = 0.5, cex = 0.8)
  #axis(side = 1, xlab = "WT",    )
}

# figure 5
aa_diff_mean_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                           ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2) * 2)
rownames(aa_diff_mean_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
aa <- (cbind(paste(colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30], "pos", sep = "_"),
             paste(colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30], "neg", sep = "_")))
colnames(aa_diff_mean_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_max_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                          ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2)* 2)
rownames(aa_diff_max_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_max_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_meangt1_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                              ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2)*2)
rownames(aa_diff_meangt1_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_meangt1_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_nugt1_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                            ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2)*2)
rownames(aa_diff_nugt1_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_nugt1_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_list_ex <- list()
for(i in 1:length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)){
  #print(sum(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2] == 0))
  aa_diff_list_ex[[i]] <- (E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]] - E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2])/E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2]
  
  for(j in 3:ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])){
    aaone <- E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 1] == 1
    aazero <- E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 1] == 0
    aa_diff_mean_sep[i, (2*(j-2) - 1)] <- mean((aa_diff_list_ex[[i]][aaone,j]))
    aavv <- aa_diff_list_ex[[i]][aaone,j]
    aa_diff_max_sep[i, (2*(j-2) - 1)] <- max(abs(aavv)) * sign(aavv[which.max(abs(aavv))])
    aa_diff_meangt1_sep[i, (2*(j-2) - 1)] <- mean(aa_diff_list_ex[[i]][(abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aaone ,j])
    aa_diff_nugt1_sep[i, (2*(j-2) - 1)] <- sum((abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aaone)/sum(aaone)
    aa_diff_mean_sep[i, (2*(j-2) )] <- mean((aa_diff_list_ex[[i]][aazero,j]))
    aavv <- aa_diff_list_ex[[i]][aazero,j]
    aa_diff_max_sep[i, (2*(j-2))] <- max(abs(aavv)) * sign(aavv[which.max(abs(aavv))])
    aa_diff_meangt1_sep[i, (2*(j-2))] <- mean(aa_diff_list_ex[[i]][(abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aazero ,j])
    aa_diff_nugt1_sep[i, (2*(j-2))] <- sum((abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aazero)/sum(aazero)
  }
  #names(aa_diff_list_ex[[i]]) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])[3:30]
}
names(aa_diff_list_ex) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)

par(mfrow = c(1,1), mar = c(7,6,4,2))
boxplot.matrix(aa_diff_mean_sep * 100, las = 2, outline= F, xaxt = 'n',
               main = "In silico perturbation", ylab = "mean % change in predicted expression",
               col = rep(c("green", "red"), 28))
# axis(side = 1, at = seq(1.5, 56.5, 2),las= 2,
#      labels = c("ER", "Half_ER", "FOXA1", "JUN", "LEF1", "NFIB", "NKX3-1",
#                 "NR2F2", "NR3C1", "NR5A2", "PBX1", "PGR", "PPARD", "PPARG",
#                 "RARA", "RUNX1", "SP1", "YBX1", "ER:ER", "ER:Half_ER",
#                 "ER:FOXA1", "ER:NR2F2", "ER:PGR", "ER:PPARD", "ER:PPARG",
#                 "ER:YBX1", "NR2F2:NR2F2","NR2F2:RUNX1"))
abline(h = seq(-100, 100, 5), col = "gray", lty = 3, lwd = 0.8)
abline(v = seq(0.5, 57.5, 2), col = "gray",lty = 2, lwd = 1 )
legend("topright", legend = c("positive", "negative"), 
       fill = c("green", "red"), 
       cex = 0.6, inset=.02,  horiz=T, 
       bty = "o", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)

boxplot.matrix(aa_diff_nugt1_sep * 100, las = 2, outline= F, xaxt = 'n', 
               ylab = "% enhancers affected by more than 5%", col = rep(c("green", "red"), 28))
abline(h = seq(0, 100, 5), col = "gray", lty = 3, lwd = 0.8)
abline(v = seq(0.5, 60.5, 2), col = "gray",lty = 2, lwd = 1 )
# axis(side = 1, at = seq(1.5, 56.5, 2),las= 2,
#      labels = c("ER", "Half_ER", "FOXA1", "JUN", "LEF1", "NFIB", "NKX3-1",
#                 "NR2F2", "NR3C1", "NR5A2", "PBX1", "PGR", "PPARD", "PPARG",
#                 "RARA", "RUNX1", "SP1", "YBX1", "ER:ER", "ER:Half_ER", 
#                 "ER:FOXA1", "ER:NR2F2", "ER:PGR", "ER:PPARD", "ER:PPARG", 
#                 "ER:YBX1", "NR2F2:NR2F2","NR2F2:RUNX1"))
# legend("topright", legend = c("positive", "negative"), 
#        fill = c("green", "red"), 
#        cex = 0.8, inset=.02,  horiz=F, 
#        bty = "o", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)


par(mfrow = c(1,1), mar = c(7,6,4,2))
#aa_diff_meangt1_sep_plot <- aa_diff_meangt1_sep
#aa_diff_meangt1_sep_plot[aa_diff_meangt1_sep_plot > 0.5]
boxplot.matrix(aa_diff_meangt1_sep * 100, las = 2, outline= F,xaxt = "n"
               , ylab = "mean % change in predicted expression\n of enhancers affected more than 5%",
               col = rep(c("green", "red"), 28)
               ,ylim = c(-50, 50)
)
abline(h = seq(-100, 300, 5), col = "gray", lty = 3, lwd = 0.8)
abline(v = seq(0.5, 60.5, 2), col = "gray",lty = 2, lwd = 1 )
abline(h = 0, col = "gray",lty = 2, lwd = 2 )
#axis.break(axis = 2,breakpos = 50, style = "slash")
axis(side = 1, at = seq(1.5, 56.5, 2),las= 2,
     labels = c("ER", "Half_ER", "FOXA1", "JUN", "LEF1", "NFIB", "NKX3-1",
                "NR2F2", "NR3C1", "NR5A2", "PBX1", "PGR", "PPARD", "PPARG",
                "RARA", "RUNX1", "SP1", "YBX1", "ER:ER", "ER:Half_ER",
                "ER:FOXA1", "ER:NR2F2", "ER:PGR", "ER:PPARD", "ER:PPARG",
                "ER:YBX1", "NR2F2:NR2F2","NR2F2:RUNX1"))
# legend("topright", legend = c("positive", "negative"),
#        fill = c("green", "red"),
#        cex = 0.8, inset=.02,  horiz=F,
#        bty = "o", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)



boxplot.matrix(aa_diff_max_sep*100, las = 2, outline= F,xaxt = "n",
               ylab= "max % change in predicted expression", col = rep(c("green", "red"), 28)
               #,ylim = c(-1, 2)
)
abline(h = seq(-500, 600, 10), col = "gray", lty = 3, lwd = 0.8)
abline(v = seq(0.5, 60.5, 2), col = "gray",lty = 3, lwd = 0.9 )
#axis.break(axis = 2,breakpos = 50, style = "slash")
axis(side = 1, at = seq(1.5, 56.5, 2),las= 2,
     labels = c("ER", "Half_ER", "FOXA1", "JUN", "LEF1", "NFIB", "NKX3-1",
                "NR2F2", "NR3C1", "NR5A2", "PBX1", "PGR", "PPARD", "PPARG",
                "RARA", "RUNX1", "SP1", "YBX1", "ER:ER", "ER:Half_ER", 
                "ER:FOXA1", "ER:NR2F2", "ER:PGR", "ER:PPARD", "ER:PPARG", 
                "ER:YBX1", "NR2F2:NR2F2","NR2F2:RUNX1"))
legend("topright", legend = c("positive", "negative"),
       fill = c("green", "red"),
       cex = 0.8, inset=.02,  horiz=F,
       bty = "o", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)


################################################################################################################
################################################################################################################
# ER paper figures

#fig1: schematic
#fig2: same as the qual fig2
# fig3:
aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas1[1:150]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas2[1:150]]
aaadd3 <- union(aaadd1, aaadd2)

aanew <- paste0("output_", unlist(lapply(lapply(strsplit(aaadd3, "_"), "[", c(2,3,4)),paste, collapse="_")))

E_RNA_GEMSTAT_Ensemble_Outlist[[18]]
E_RNA_GEMSTAT_Ensemble_Parlist[[18]]
par(mfrow = c(1, 2), mar = c(4,4,4,4))

aacol <- integer(length(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$train_ROC))
aacol[] <- 1
aacol[names(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$train_ROC) %in% aanew] <- 2
#aacol[E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_ROC >= 0.60 & E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_PRC >= 0.32] <- 2
# plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_ROC , E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$train_ROC,
#      main = "AUROC", xlab = "test", ylab = "training", pch = 16, cex = 0.8, col = aacol)
# plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_PRC, E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$train_PRC, 
#      main = "AUPRC", xlab = "test", ylab = "training", pch = 16, cex = 0.8, col = aacol)

par(mfrow = c(1, 2), mar = c(4,4,4,4))
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_ROC , E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$train_ROC,
     main = "AUROC", xlab = "validation", ylab = "training", pch = 16, cex = 0.8, col = aacol)
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_PRC, E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$train_PRC, 
     main = "AUPRC", xlab = "validation", ylab = "training", pch = 16, cex = 0.8, col = aacol)

# plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_ROC[aacol == 2] , E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_ROC[aacol == 2],
#      main = "AUROC", xlab = "validation", ylab = "test", pch = 16, cex = 0.8)
# plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_PRC[aacol == 2], E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_PRC[aacol == 2],
#      main = "AUPRC", xlab = "validation", ylab = "test", pch = 16, cex = 0.8)

# plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_PRC[aacol == 2], E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$valid_ROC[aacol == 2],
#      main = "", xlab = "validation AUPRC", ylab = "validation AUROC", pch = 16, cex = 0.8)

par(mfrow = c(1, 1), mar = c(4,4,4,4))
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$test_PRC[aacol == 2], 
     E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$test_ROC[aacol == 2], 
     main = "", xlab = "test AUPRC", ylab = "test AUROC", pch = 16, cex = 0.8)

aaw <- which(aacol == 2)
which.max(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$test_PRC)
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_results$output_trained_par_1606_1_valid$ROC_curve,
     main = "valid ROC")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_results$output_trained_par_1606_1_valid$PRC_curve,
     main = "valid PRC")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$test_results$output_trained_par_3278_1_test$ROC_curve,
     main = "test ROC")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$test_results$output_trained_par_3278_1_test$PRC_curve,
     main = "test PRC")


# fig4 


aa_average_percentile_change_bylabel <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list), 
                                               ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list[[i]])-2)*2)
aa_average_percentile_change_gt5_bylabel <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list), 
                                                   ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list[[i]])-2)*2)
aa_average_percentile_change_bylabel_noabs <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list), 
                                                     ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list[[i]])-2)*2)
aa_average_percentile_change_gt5_bylabel_noabs <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list), 
                                                         ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list[[i]])-2)*2)

aa_average_percentile_change_percent_bylabel <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list), 
                                                       ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list[[i]])-2)*2)
colnames(aa_average_percentile_change_bylabel) <- character(ncol(aa_average_percentile_change_bylabel))
colnames(aa_average_percentile_change_gt5_bylabel) <- character(ncol(aa_average_percentile_change_gt5_bylabel))
colnames(aa_average_percentile_change_bylabel_noabs) <- character(ncol(aa_average_percentile_change_bylabel_noabs))
colnames(aa_average_percentile_change_gt5_bylabel_noabs) <- character(ncol(aa_average_percentile_change_gt5_bylabel_noabs))
colnames(aa_average_percentile_change_percent_bylabel) <- character(ncol(aa_average_percentile_change_percent_bylabel))



aacnames <- c("ER", "FOXA1", "GATA3" ,"JUN", "LEF1", "NKX3-1",
              "NR3C1", "NR5A2","PAX2", "PBX1", "PGR",
              "RARA","RELA", "RUNX1", "SP1","TFAP2C", "YBX1",
              "ER:FOXA1", "ER:GATA3", "ER:PBX1", "ER:PGR", "ER:RARA","ER:TFAP2C",
              "ER:YBX1", "JUN_1:JUN_1","YBX1:RUNX1")

for(i in 1:length(Exp18_KD_percentile_change_list)){
  for(j in 1:ncol(Exp18_KD_percentile_change_list[[i]])){
    awpos <- which(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list[[i]][, 1] == 1)
    awneg <- which(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list[[i]][, 1] == 0)
    aa_average_percentile_change_bylabel[i, 2*j - 1] <- mean(abs(Exp18_KD_percentile_change_list[[i]][awpos, j]))
    aa_average_percentile_change_bylabel[i, 2*j    ] <- mean(abs(Exp18_KD_percentile_change_list[[i]][awneg, j]))
    
    aa_average_percentile_change_bylabel_noabs[i, 2*j - 1] <- mean(Exp18_KD_percentile_change_list[[i]][awpos, j])
    aa_average_percentile_change_bylabel_noabs[i, 2*j    ] <- mean(Exp18_KD_percentile_change_list[[i]][awneg, j])
    
    
    # aa_average_percentile_change[i, j] <- mean(abs(Exp18_KD_percentile_change_list[[i]][, j]))
    aasum <- which(abs(Exp18_KD_percentile_change_list[[i]][, j]) >= 0.05)
    aasump <- intersect(aasum, awpos)
    aasumn <- intersect(aasum, awneg)
    aa_average_percentile_change_percent_bylabel[i, 2*j - 1] <- length(aasump)/length(awpos)*100
    aa_average_percentile_change_percent_bylabel[i, 2*j ] <- length(aasumn)/length(awneg)*100
    if(length(aasump) >= 5){
      aa_average_percentile_change_gt5_bylabel[i, 2*j - 1] <- mean(abs(Exp18_KD_percentile_change_list[[i]][aasump, j]))
      aa_average_percentile_change_gt5_bylabel_noabs[i, 2*j - 1] <- mean(Exp18_KD_percentile_change_list[[i]][aasump, j])
    }
    if(length(aasumn) >= 5){
      aa_average_percentile_change_gt5_bylabel[i, 2*j ] <- mean(abs(Exp18_KD_percentile_change_list[[i]][aasumn, j]))
      aa_average_percentile_change_gt5_bylabel_noabs[i, 2*j ] <- mean(Exp18_KD_percentile_change_list[[i]][aasumn, j])
      
    }
    colnames(aa_average_percentile_change_bylabel)[c(2*j - 1, 2*j)] <- c(paste0(aacnames[j], "_pos"), paste0(aacnames[j], "_neg"))
    colnames(aa_average_percentile_change_bylabel_noabs)[c(2*j - 1, 2*j)] <- c(paste0(aacnames[j], "_pos"), paste0(aacnames[j], "_neg"))
    colnames(aa_average_percentile_change_gt5_bylabel)[c(2*j - 1, 2*j)] <- c(paste0(aacnames[j], "_pos"), paste0(aacnames[j], "_neg"))
    colnames(aa_average_percentile_change_gt5_bylabel_noabs)[c(2*j - 1, 2*j)] <- c(paste0(aacnames[j], "_pos"), paste0(aacnames[j], "_neg"))
    colnames(aa_average_percentile_change_percent_bylabel)[c(2*j - 1, 2*j)] <- c(paste0(aacnames[j], "_pos"), paste0(aacnames[j], "_neg"))
  }
}

par(mfrow = c(1,1), mar = c(8,5,2,2))

boxplot.matrix(aa_average_percentile_change_bylabel_noabs*100,
               las = 2, main = "In-silico TF purterbation", 
               col = rep(c("green", "red"), 26), xaxt = "n",  
               ylab = "mean percentile change in predicted expression", outline=F)
abline(h = seq(-100, 300, 5), col = "gray", lty = 3, lwd = 0.8)
abline(v = seq(0.5, 60.5, 2), col = "gray",lty = 2, lwd = 1 )
abline(h = 0, col = "gray",lty = 2, lwd = 2 )
#axis.break(axis = 2,breakpos = 50, style = "slash")
axis(side = 1, at = seq(1.5, 52.5, 2),las= 2,
     labels =  c("ER", "FOXA1", "GATA3" ,"JUN", "LEF1", "NKX3-1",
                 "NR3C1", "NR5A2","PAX2", "PBX1", "PGR",
                 "RARA","RELA", "RUNX1", "SP1","TFAP2C", "YBX1",
                 "ER:FOXA1", "ER:GATA3", "ER:PBX1", "ER:PGR", "ER:RARA","ER:TFAP2C",
                 "ER:YBX1", "JUN_1:JUN_1","YBX1:RUNX1"))
legend("bottomright", legend = c("positive", "negative"),
       fill = c("green", "red"),
       cex = 0.8, inset=.02,  horiz=F,
       bty = "o", xjust = 1, pt.cex = 0.5, x.intersp = 0.2)

boxplot.matrix(aa_average_percentile_change_percent_bylabel,
               las = 2, main = "", ylab = "% enhancers affected by more than 5 percentile",
               col = rep(c("green", "red"), 26), xaxt = "n", outline = F)
abline(h = seq(-100, 300, 5), col = "gray", lty = 3, lwd = 0.8)
abline(v = seq(0.5, 60.5, 2), col = "gray",lty = 2, lwd = 1 )
#abline(h = 0, col = "gray",lty = 2, lwd = 2 )
#axis.break(axis = 2,breakpos = 50, style = "slash")
axis(side = 1, at = seq(1.5, 52.5, 2),las= 2,
     labels =  c("ER", "FOXA1", "GATA3" ,"JUN", "LEF1", "NKX3-1",
                 "NR3C1", "NR5A2","PAX2", "PBX1", "PGR",
                 "RARA","RELA", "RUNX1", "SP1","TFAP2C", "YBX1",
                 "ER:FOXA1", "ER:GATA3", "ER:PBX1", "ER:PGR", "ER:RARA","ER:TFAP2C",
                 "ER:YBX1", "JUN_1:JUN_1","YBX1:RUNX1"))



boxplot.matrix(aa_average_percentile_change_gt5_bylabel_noabs *100,
               las = 2, main = "", 
               col = rep(c("green", "red"), 26), xaxt = "n"
               ,ylab = "mean percentile change in predicted expression\n of enhancers affected more than 5 percentile"
               , outline=F)
abline(h = seq(-100, 300, 5), col = "gray", lty = 3, lwd = 0.8)
abline(v = seq(0.5, 60.5, 2), col = "gray",lty = 2, lwd = 1 )
abline(h = 0, col = "gray",lty = 2, lwd = 2 )
#axis.break(axis = 2,breakpos = 50, style = "slash")
axis(side = 1, at = seq(1.5, 52.5, 2),las= 2,
     labels =  c("ER", "FOXA1", "GATA3" ,"JUN", "LEF1", "NKX3-1",
                 "NR3C1", "NR5A2","PAX2", "PBX1", "PGR",
                 "RARA","RELA", "RUNX1", "SP1","TFAP2C", "YBX1",
                 "ER:FOXA1", "ER:GATA3", "ER:PBX1", "ER:PGR", "ER:RARA","ER:TFAP2C",
                 "ER:YBX1", "JUN_1:JUN_1","YBX1:RUNX1"))






aats <- data.frame(mean_percentile=as.numeric(aa_average_percentile_change_bylabel_noabs) * 100,
                   mean_percentile_gt5 =as.numeric(aa_average_percentile_change_gt5_bylabel_noabs) * 100,
                   mean_percent_affected_gt5=as.numeric(aa_average_percentile_change_percent_bylabel),
                   experiment = as.factor(rep(colnames(aa_average_percentile_change_bylabel_noabs), 
                                              each = nrow(aa_average_percentile_change_bylabel_noabs))),
                   model_name=as.factor(rep(names(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list),
                                            ncol(aa_average_percentile_change_bylabel_noabs)))
                    ,
                    KD= as.factor(unlist(lapply(strsplit(rep(colnames(aa_average_percentile_change_bylabel_noabs), 
                                                      each = nrow(aa_average_percentile_change_bylabel_noabs)), "_"), "[[", 1)))
                   )


#compute t-test between different classes of enhancers on each variable
aa_sp_first <- unlist(mapply(lapply(strsplit(levels(aats$experiment)[as.numeric(aats$experiment)], "_"), 
                     function(x) x[1:(length(x) - 1)]), FUN = paste, collapse = "_"))
aa_sp_last <- unlist(lapply(strsplit(levels(aats$experiment)[as.numeric(aats$experiment)], "_"), 
                                    function(x) x[length(x)]))
aa_ttest_frac <- array(dim =  length(unique(aa_sp_first)))
names(aa_ttest_frac) <- unique(aa_sp_first)

aa_ttest_pcgt5 <- array(dim =  length(unique(aa_sp_first)))
names(aa_ttest_pcgt5) <- unique(aa_sp_first)

aa_wilc_frac <- array(dim =  length(unique(aa_sp_first)))
names(aa_wilc_frac) <- unique(aa_sp_first)

aa_wilc_pcgt5 <- array(dim =  length(unique(aa_sp_first)))
names(aa_wilc_pcgt5) <- unique(aa_sp_first)


sink("~/Documents/Shayan/BioInf/EstrogenReceptor/ttest_sink.txt")
for(i in 1:length(aa_ttest_frac)){
  print(names(aa_ttest_frac)[i])
  aa_x_ind <- which((aa_sp_first %in% names(aa_ttest_frac)[i]) & (aa_sp_last  %in% "pos"))
  aa_y_ind <-  which((aa_sp_first %in% names(aa_ttest_frac)[i]) & (aa_sp_last  %in% "neg"))
  aa_cur_f1 <- t.test(x = aats$mean_percent_affected_gt5[aa_x_ind], 
                     y = aats$mean_percent_affected_gt5[aa_y_ind],
                     alternative = "two.sided")
  print("frac ttest")
  print(aa_cur_f1)
  aa_ttest_frac[i] <- (aa_cur_f1$p.value)
  
  aa_cur_f3 <- wilcox.test(x = aats$mean_percent_affected_gt5[aa_x_ind], 
                      y = aats$mean_percent_affected_gt5[aa_y_ind],
                      alternative = "two.sided")
  print("frac wilc")
  print(aa_cur_f3)
  aa_wilc_frac[i] <- (aa_cur_f3$p.value)
  if((sum(!is.na(aats$mean_percentile_gt5[aa_x_ind]))>=5) &  (sum(!is.na(aats$mean_percentile_gt5[aa_y_ind]))>=5)){
    aa_cur_f2 <- t.test(x = aats$mean_percentile_gt5[aa_x_ind], 
                        y = aats$mean_percentile_gt5[aa_y_ind],
                        alternative = "two.sided")
    aa_ttest_pcgt5[i] <- (aa_cur_f2$p.value)
    print("percentile gt5")
    print(aa_cur_f2)

    
    aa_cur_f4 <- wilcox.test(x = aats$mean_percentile_gt5[aa_x_ind], 
                             y = aats$mean_percentile_gt5[aa_y_ind],
                             alternative = "two.sided")
    print("gt5 wilc")
    print(aa_cur_f4)
    aa_wilc_pcgt5[i] <- (aa_cur_f4$p.value)
  }
}
aa_ttest_frac_adj <- p.adjust(aa_ttest_frac, method = "BH")
aa_ttest_pcgt5_adj <- p.adjust(aa_ttest_pcgt5, method = "BH")
aa_wilc_frac_adj <- p.adjust(aa_wilc_frac, method = "BH")
aa_wilc_pcgt5_adj <- p.adjust(aa_wilc_pcgt5, method = "BH")

sink()


aatgc <- summarySE(aats, measurevar=c("mean_percentile"), 
                   groupvars=c("experiment"), na.rm=F)

aatgc_alt <- summarySE(aats, measurevar=c("mean_percentile_gt5"), 
                   groupvars=c("KD"), na.rm=T)
aatgc_alt2 <- summarySE(aats, measurevar=c("mean_percent_affected_gt5"), 
                       groupvars=c("KD"), na.rm=T)

aatgc$class <- as.factor(rep(c("Negative", "Positive"), nrow(aatgc)/2))
aatgc$KD <- unlist(lapply(strsplit(as.character(levels(aatgc$experiment)[as.numeric(aatgc$experiment)]), "_"), "[[", 1))
aatgc$KD[c(21,22)] <- c("JUN_1:JUN_1", "JUN_1:JUN_1")
# aatgc$KD <- as.factor(aatgc$KD)
aatgc2 <- aatgc[match(colnames(aa_average_percentile_change_bylabel_noabs),levels(aatgc$experiment)[as.numeric(aatgc$experiment)]),]
aatgc <- aatgc2
aatgc$KD <- factor(aatgc$KD, levels = unique(aatgc$KD))

aa_fig4_a <- ggplot(aatgc,
                    aes(x=KD,
                        y=mean_percentile,
                        fill = class)) + 
  geom_bar(#mapping = aes(group = KD),
    position=position_dodge(width = .8), 
    width = 0.9,
    stat="identity", show.legend = T ) +
  geom_errorbar(aes(ymin=mean_percentile-sd,
                    ymax=mean_percentile+sd),
                width=.3,                    # Width of the error bars
                #position="dodge"
                position=position_dodge(.8)
  ) + 
  ggtitle("In-silico Knockdown")+
  ylab("mean percentile change in predicted expression\n ")+
  theme_minimal()+
  scale_fill_brewer(name="Enhancer class", palette = "Set1") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,
                                   colour="black", size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_text(colour="black", size = 10, face = "bold"),
        legend.justification=c(1,1), 
        legend.position=c(1,1),
        legend.text = element_text(colour="black", size = 10
                                   #, face = "bold"
                                   ),
        legend.background = element_rect(fill="gray90", size=.5, linetype=NULL),
        axis.title.y = element_text(colour="black", size = 11, face = "bold"), 
        plot.title = element_text(colour="black", size = 12,
                                  face = "bold", hjust = .5, vjust = .5)

  ) +
  geom_vline(xintercept = seq(0.5,53,1), color = "grey", size=0.2) 
  #scale_fill_manual()
# + 
#   facet_wrap(~class)
aa_fig4_a

aatgc <- summarySE(aats, measurevar=c("mean_percent_affected_gt5"), 
                   groupvars=c("experiment"), na.rm=F)
aatgc$class <- as.factor(rep(c("Negative", "Positive"), nrow(aatgc)/2))
aatgc$KD <- unlist(lapply(strsplit(as.character(levels(aatgc$experiment)[as.numeric(aatgc$experiment)]), "_"), "[[", 1))
aatgc$KD[c(21,22)] <- c("JUN_1:JUN_1", "JUN_1:JUN_1")
# aatgc$KD <- as.factor(aatgc$KD)
aatgc2 <- aatgc[match(colnames(aa_average_percentile_change_bylabel_noabs),levels(aatgc$experiment)[as.numeric(aatgc$experiment)]),]
aatgc2 <- aatgc2[1:(nrow(aatgc2) - 4),]
aatgc <- aatgc2
aatgc$KD <- factor(aatgc$KD, levels = unique(aatgc$KD))

aa_fig4_b <- ggplot(aatgc,
                    aes(x=KD,
                        y=mean_percent_affected_gt5,
                        fill = class)) + 
  geom_hline(yintercept = seq(0,60,20), color = "grey", size=0.2)+
  geom_errorbar(aes(ymin= mean_percent_affected_gt5-sd,
                    ymax=mean_percent_affected_gt5+sd),
                width=.5,                    # Width of the error bars
                #position="dodge"
                position=position_dodge(.8)
  ) + 
  geom_bar(#mapping = aes(group = KD),
    position=position_dodge(width = .8), 
    width = 0.9,
    stat="identity", show.legend = T ) +


  ggtitle("In-silico Knockdown")+
  ylab("% enhancers affected by more than 5 percentile\n ")+
  theme_minimal()+
  scale_fill_brewer(name="Enhancer class", palette = "Set1") +
  theme(axis.title.x=element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.grid.major = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
         legend.title=element_text(colour="black", size = 11, face = "bold"),
         legend.justification=c(1,1), 
         legend.position=c(1,1),
         legend.text = element_text(colour="black", size = 11
                                   , face = "bold"
        ),
        legend.background = element_rect(fill="gray90", size=.5, linetype=NULL)
        ,axis.title.y = element_text(colour="black", size = 11),
        plot.title = element_text(colour="black", size = 13, face = "bold",
                                  hjust = .5, vjust = .5)
        
  ) +
  geom_vline(xintercept = seq(0.5,53,1), color = "grey", size=0.2)
# + 
#   facet_wrap(~class)

png("~/Desktop/ER_paper/fig_4_a.png",         
    width = 8*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 14)        #  font size
aa_fig4_b
dev.off()

aatgc <- summarySE(aats, measurevar=c("mean_percentile_gt5"), 
                   groupvars=c("experiment"), na.rm=T)
aatgc$class <- as.factor(rep(c("Negative", "Positive"), nrow(aatgc)/2))
aatgc$KD <- unlist(lapply(strsplit(as.character(levels(aatgc$experiment)[as.numeric(aatgc$experiment)]), "_"), "[[", 1))
aatgc$KD[c(21,22)] <- c("JUN_1:JUN_1", "JUN_1:JUN_1")
# aatgc$KD <- as.factor(aatgc$KD)
aatgc2 <- aatgc[match(colnames(aa_average_percentile_change_bylabel_noabs),levels(aatgc$experiment)[as.numeric(aatgc$experiment)]),]
aatgc2 <- aatgc2[1:(nrow(aatgc2) - 4),]
aatgc <- aatgc2
aatgc$KD <- factor(aatgc$KD, levels = unique(aatgc$KD))

aa_fig4_c <- ggplot(aatgc,
                    aes(x=KD,
                        y=mean_percentile_gt5,
                        fill = class)) + 
  geom_hline(yintercept = seq(-20,20,20), color = "grey", size=0.2)+
  geom_errorbar(aes(ymin=mean_percentile_gt5-  sd,
                    ymax=mean_percentile_gt5+  sd),
                width=.5,                    # Width of the error bars
                #position="dodge"
                position=position_dodge(.8)
  ) + 
  geom_bar(#mapping = aes(group = KD),
    position=position_dodge(width = .8), 
    width = 0.9,
    stat="identity", show.legend = F ) +

  ggtitle("")+
  ylab("Mean percentile change in  expression\n")+
  theme_minimal()+
  scale_fill_brewer(name="Enhancer class", palette = "Set1") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,colour="black", size = 11, face = "bold"),
        panel.grid.major = element_blank()#,
       # axis.text.x=element_blank(),
       # axis.ticks.x=element_blank()#,
        # legend.title=element_text(colour="black", size = 10, face = "bold"),
        # legend.justification=c(1,1), 
        # legend.position=c(1,1),
        # legend.text = element_text(colour="black", size = 10
        #, face = "bold"
        #),
        #legend.background = element_rect(fill="gray90", size=.5, linetype=NULL)
       ,axis.title.y = element_text(colour="black", size = 11, face = "bold")
        
  ) +
  geom_vline(xintercept = seq(0.5,53,1), color = "grey", size=0.2) 
# + 
#   facet_wrap(~class)
png("~/Desktop/ER_paper/fig_4_b_short.png",         
    width = 8*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 14)        #  font size
aa_fig4_c
dev.off()

# adding KDs dor three TFs KDs
# output_trained_par_3278_1_test
par(mfrow = c(1, 1), mar = c(4,4,4,4))
plot( E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_1026_1[, 2], 
      E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_1026_1[, "ESR1_2"], 
      xaxt = "n", yaxt = "n",ylim = c(0.05, 0.633), xlim = c(0.05, 0.633), 
      ylab="", xlab = "",  pch = 20, cex = 0.6
      #,main = colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_420_1)[3]
      )
mtext(text = "ER KD", side = 2, outer = F, line = 0.5, cex = 1)
mtext(text = "WT", side = 1, outer = F, line = 0.5, cex = 1)
plot( E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_1026_1[, 2], 
      E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_1026_1[, "NKX3_1"], 
      xaxt = "n", yaxt = "n",
      ylim = range(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_1026_1[, 2]),
      xlim = range(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_1026_1[, 2]),
      ylab="", xlab = "",  pch = 20, cex = 0.6
      #,main = colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_420_1)[3]
      )
mtext(text = "NKX3_1 KD", side = 2, outer = F, line = 0.5, cex = 1)
mtext(text = "WT", side = 1, outer = F, line = 0.5, cex = 1)

plot( E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_1026_1[, 2], 
      E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_1026_1[, "ESR1_2__FOXA1"], 
      xaxt = "n", yaxt = "n",
      ylim = range(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_1026_1[, 2]),
      xlim = range(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_1026_1[, 2]),
      ylab="", xlab = "",  pch = 20, cex = 0.6
      #,main = colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[7]]$output_mat_list$output_par_420_1)[3]
)
mtext(text = "ER:FOXA1 interaction KD", side = 2, outer = F, line = 0.5, cex = 1)
mtext(text = "WT", side = 1, outer = F, line = 0.5, cex = 1)

#mtext(text = "ER", side = 3, outer = F, line = 1, cex = 1)
#axis(side = 1, xlab = "WT",    )
#}
# panel d is made in ER_GRO_KD.R script starting line 143


# grid.arrange(aa_fig4_a, aa_fig4_b, aa_fig4_c, ncol=1)
# 
# aa_fig4_a
# aa_fig4_b
# aa_fig4_c



################################## fig5
# aalab_orig <- ER_SNPs_merged_allbypos_modified$is_gtex
# aa_som_rem <- which(ER_SNPs_merged_allbypos_modified$is_somatic)
# aa_cli_rem <- which(ER_SNPs_merged_allbypos_modified$is_clinical)
# aa_gt <- which(ER_SNPs_merged_allbypos_modified$is_gtex)
# aa_cm <- which(ER_SNPs_merged_allbypos_modified$is_common)
# aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
# Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
# aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
# aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[aasap %in% union(aa_gt,aa_cm),]
# aasap <- aasap[aasap %in% union(aa_gt,aa_cm)]
# aaTP <- aalab_orig[aasap]
# 
# 
# 
# 
# aa_gtex_ornot_pc <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
#                                            rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
#                                            TP = aaTP, 
#                                            sampling_range=seq(500,5000, 500),#c(seq(10,300,10), seq(400,2000, 100), seq(3000,33000,1000)),
#                                            do_boxplot=F, 
#                                            percent_in_sample = T)
# 
# 
# # aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_gtex_ornot$Random_scores) * 2)
# # #colnames(aa_allmat) <- aasample_size
# # aac <- 1
# # for(i in seq(1, ncol(aa_allmat), 2)){
# #   aa_allmat[, i] <- aa_gtex_ornot$Random_scores[, aac]
# #   aa_allmat[, i+1] <- t(aa_gtex_ornot$Model_scores)[, aac]
# #   aac <- aac + 1
# # }
# # aas <- sort(c(seq(1,100, 6),seq(2,100, 6)))
# # 
# # boxplot.matrix(aa_allmat[,aas], col = rep(c(2,3), length(aas)/2),
# #                xaxt = "n", xlab = "#predictions", ylab = "TPR", 
# #                main = "gtex ranking")
# # axis(side = 1, at = seq(1.5,ncol(aa_allmat[,aas]), 2), 
# #      labels = seq(100,5000, 100)[seq(2,100, 6)/2],#c(seq(10,300,10), seq(400,2000, 100), seq(3000,33000,1000)),
# #      las = 2)
# # abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = "grey", lwd = 0.7, lty = 3)
# # legend(x = "topleft",legend = c("location-based", "model-based"), fill = c(2,3), cex = 0.6)
# 
# # simulate random
# aa <- numeric(length = (nu_uniq_common_snp + (nu_uniq_gtex_breast_Cancer_eqtls - nu_uniq_gtex_overlap_common)))
# aa[1:nu_uniq_gtex_breast_Cancer_eqtls] <- 1
# aa_samp_size <- seq(500,5000, 500)
# aa_score <- numeric(length(aa_samp_size)*244)
# 
# 
# aa_cnt <- 1
# for(i in 1:length(aa_samp_size)){
#   
#   for(j in 1:244){
#     print("i")
#     print(i)
#     print("j")
#     print(j)
#     aasamp <- sample(x = aa, size = aa_samp_size[i], replace = F)
#     aa_score[aa_cnt] <- sum(aasamp)/aa_samp_size[i]
#     aa_cnt <- aa_cnt + 1
#   }
# }
# 
# 
# aa1 <- data.frame(Score=as.numeric(aa_gtex_ornot_pc$Random_scores), 
#         Samp_nu = factor(rep(colnames(aa_gtex_ornot_pc$Random_scores),
#             each = nrow(aa_gtex_ornot_pc$Random_scores)),
#             levels = colnames(aa_gtex_ornot_pc$Random_scores)), 
#         Model_nu = as.numeric(rep(c(1:nrow(aa_gtex_ornot_pc$Random_scores)),
#                                   ncol(aa_gtex_ornot_pc$Random_scores))),
#         Method= rep("Location-based", length(as.numeric(aa_gtex_ornot_pc$Random_scores))))
# 
# aa2 <- data.frame(Score=as.numeric(t(aa_gtex_ornot_pc$Model_scores)), 
#                   Samp_nu = factor(as.numeric(rep(colnames(t(aa_gtex_ornot_pc$Model_scores)),
#                                            each = nrow(t(aa_gtex_ornot_pc$Model_scores)))),
#                                    levels = colnames(t(aa_gtex_ornot_pc$Model_scores))), 
#                   Model_nu = as.numeric(rep(c(1:nrow(t(aa_gtex_ornot_pc$Model_scores))),
#                                             ncol(t(aa_gtex_ornot_pc$Model_scores)))),
#                   Method= rep("Model-based", length(as.numeric(t(aa_gtex_ornot_pc$Model_scores)))))
# 
# 
# aa3 <- data.frame(Score=aa_score, 
#                   Samp_nu = factor(rep(as.character(aa_samp_size),
#                                        each = 244),
#                                    levels = as.character(aa_samp_size)), 
#                   Model_nu = as.numeric(rep(c(1:244),
#                                             length(aa_samp_size))),
#                   Method= rep("Random", length(aa_score)))
# 
# aaal_gtex <- rbind(aa3, aa1, aa2)
# 
# ggplot(data = aaal_gtex, aes(x = Samp_nu, y = Score, fill = Method))+
#   geom_boxplot()+
#   ggtitle("Breast Cancer GTEX eqtl ranking")+
#   ylab("TPR")+
#   xlab("#predictions")+
#   theme_minimal()
#   
# ######################################################## pancanqtl
# aalab_orig <- ER_SNPs_merged_allbypos_modified$is_pancan
# aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
# Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
# aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
# aaTP <- aalab_orig[aasap]
# 
# aalab_orig <-  ER_SNPs_merged_allbypos_modified$is_pancan
# aalab_pancan <- which(ER_SNPs_merged_allbypos_modified$is_pancan)
# aalab_cm <- which(ER_SNPs_merged_allbypos_modified$is_common)
# aa_som_rem <- which(ER_SNPs_merged_allbypos_modified$is_somatic)
# aa_cli_rem <- which(ER_SNPs_merged_allbypos_modified$is_clinical)
# 
# aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
# Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
# aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
# aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[aasap %in% union(aalab_pancan,aalab_cm),]
# aasap <- aasap[aasap %in% union(aalab_pancan,aalab_cm)]
# aaTP <- aalab_orig[aasap]
# 
# aa_pancan_ornot_pc <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
#                                              rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
#                                              TP = aaTP, 
#                                              sampling_range=seq(500,5000, 500),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)),
#                                              do_boxplot=F,
#                                              percent_in_sample = T)
# 
# 
# # simulate random
# aa <- numeric(length = (nu_uniq_common_snp + (nu_unique_pancan_eqtls - nu_uniq_pancan_overlap_common)))
# aa[1:nu_unique_pancan_eqtls] <- 1
# aa_samp_size <- seq(500,5000, 500)
# aa_score <- numeric(length(aa_samp_size)*244)
# 
# 
# aa_cnt <- 1
# for(i in 1:length(aa_samp_size)){
#   
#   for(j in 1:244){
#     print("i")
#     print(i)
#     print("j")
#     print(j)
#     aasamp <- sample(x = aa, size = aa_samp_size[i], replace = F)
#     aa_score[aa_cnt] <- sum(aasamp)/aa_samp_size[i]
#     aa_cnt <- aa_cnt + 1
#   }
# }
# 
# 
# aa1 <- data.frame(Score=as.numeric(aa_pancan_ornot_pc$Random_scores), 
#                   Samp_nu = factor(rep(colnames(aa_pancan_ornot_pc$Random_scores),
#                                        each = nrow(aa_pancan_ornot_pc$Random_scores)),
#                                    levels = colnames(aa_pancan_ornot_pc$Random_scores)), 
#                   Model_nu = as.numeric(rep(c(1:nrow(aa_pancan_ornot_pc$Random_scores)),
#                                             ncol(aa_pancan_ornot_pc$Random_scores))),
#                   Method= rep("Location-based", length(as.numeric(aa_pancan_ornot_pc$Random_scores))))
# 
# aa2 <- data.frame(Score=as.numeric(t(aa_pancan_ornot_pc$Model_scores)), 
#                   Samp_nu = factor(as.numeric(rep(colnames(t(aa_pancan_ornot_pc$Model_scores)),
#                                                   each = nrow(t(aa_pancan_ornot_pc$Model_scores)))),
#                                    levels = colnames(t(aa_pancan_ornot_pc$Model_scores))), 
#                   Model_nu = as.numeric(rep(c(1:nrow(t(aa_pancan_ornot_pc$Model_scores))),
#                                             ncol(t(aa_pancan_ornot_pc$Model_scores)))),
#                   Method= rep("Model-based", length(as.numeric(t(aa_pancan_ornot_pc$Model_scores)))))
# 
# 
# aa3 <- data.frame(Score=aa_score, 
#                   Samp_nu = factor(rep(as.character(aa_samp_size),
#                                        each = 244),
#                                    levels = as.character(aa_samp_size)), 
#                   Model_nu = as.numeric(rep(c(1:244),
#                                             length(aa_samp_size))),
#                   Method= rep("Random", length(aa_score)))
# 
# aaal_pancan <- rbind(aa3, aa1, aa2)
# 
# ggplot(data = aaal_pancan, aes(x = Samp_nu, y = Score, fill = Method))+
#   geom_boxplot()+
#   ggtitle("Breast Cancer pancan eqtl ranking")+
#   ylab("TPR")+
#   xlab("#predictions")+
#   theme_minimal()

######################################################## gtex or pancan
aalab_orig <- ER_SNPs_merged_allbypos_modified$is_pancan | ER_SNPs_merged_allbypos_modified$is_gtex
aa_pan_keep <- which(ER_SNPs_merged_allbypos_modified$is_pancan)
aa_com_keep <- which(ER_SNPs_merged_allbypos_modified$is_common)
aa_gte_keep <- which(ER_SNPs_merged_allbypos_modified$is_gtex)

aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[aasap %in% union(aa_gte_keep,union(aa_pan_keep,aa_com_keep)),]
aasap <- aasap[aasap %in% union(aa_gte_keep,union(aa_pan_keep,aa_com_keep))]
aaTP <- aalab_orig[aasap]
aa_pancan_gtex_ornot_pc <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
                                                  rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
                                                  TP = aaTP, 
                                                  sampling_range=seq(500,5000, 500), #c(seq(10,400,10), seq(500,3000, 100), seq(3000,31000,1000)),
                                                  do_boxplot=F, 
                                                  percent_in_sample = T)


# simulate random
aa <- numeric(length = (nu_uniq_common_snp + (nu_unique_pancan_eqtls - nu_uniq_pancan_overlap_common) + (nu_uniq_gtex_breast_Cancer_eqtls - nu_uniq_gtex_overlap_common)))
aa[1:(nu_unique_pancan_eqtls + nu_uniq_gtex_breast_Cancer_eqtls - nu_pancan_gtex_overlap)] <- 1
aa_samp_size <- seq(500,5000, 500)
aa_score <- numeric(length(aa_samp_size)*244)


aa_cnt <- 1
for(i in 1:length(aa_samp_size)){
  
  for(j in 1:244){
    print("i")
    print(i)
    print("j")
    print(j)
    aasamp <- sample(x = aa, size = aa_samp_size[i], replace = F)
    aa_score[aa_cnt] <- sum(aasamp)/aa_samp_size[i]
    aa_cnt <- aa_cnt + 1
  }
}


aa1 <- data.frame(Score=as.numeric(aa_pancan_gtex_ornot_pc$Random_scores), 
                  Samp_nu = factor(rep(colnames(aa_pancan_gtex_ornot_pc$Random_scores),
                                       each = nrow(aa_pancan_gtex_ornot_pc$Random_scores)),
                                   levels = colnames(aa_pancan_gtex_ornot_pc$Random_scores)), 
                  Model_nu = as.numeric(rep(c(1:nrow(aa_pancan_gtex_ornot_pc$Random_scores)),
                                            ncol(aa_pancan_gtex_ornot_pc$Random_scores))),
                  Method= rep("Location-based", length(as.numeric(aa_pancan_gtex_ornot_pc$Random_scores))))

aa2 <- data.frame(Score=as.numeric(t(aa_pancan_gtex_ornot_pc$Model_scores)), 
                  Samp_nu = factor(as.numeric(rep(colnames(t(aa_pancan_gtex_ornot_pc$Model_scores)),
                                                  each = nrow(t(aa_pancan_gtex_ornot_pc$Model_scores)))),
                                   levels = colnames(t(aa_pancan_gtex_ornot_pc$Model_scores))), 
                  Model_nu = as.numeric(rep(c(1:nrow(t(aa_pancan_gtex_ornot_pc$Model_scores))),
                                            ncol(t(aa_pancan_gtex_ornot_pc$Model_scores)))),
                  Method= rep("Model-based", length(as.numeric(t(aa_pancan_gtex_ornot_pc$Model_scores)))))


aa3 <- data.frame(Score=aa_score, 
                  Samp_nu = factor(rep(as.character(aa_samp_size),
                                       each = 244),
                                   levels = as.character(aa_samp_size)), 
                  Model_nu = as.numeric(rep(c(1:244),
                                            length(aa_samp_size))),
                  Method= rep("Random", length(aa_score)))

aaal_pancan_gtex <- rbind(aa2, aa1, aa3)

t.test(aaal_pancan_gtex$Score[aaal_pancan_gtex$Method == "Model-based" & aaal_pancan_gtex$Samp_nu == 1500], 
       aaal_pancan_gtex$Score[aaal_pancan_gtex$Method == "Location-based" & aaal_pancan_gtex$Samp_nu == 1500], 
       paired = F)$p.value
t.test(aaal_pancan_gtex$Score[aaal_pancan_gtex$Method == "Random" & aaal_pancan_gtex$Samp_nu == 1500], 
       aaal_pancan_gtex$Score[aaal_pancan_gtex$Method == "Location-based" & aaal_pancan_gtex$Samp_nu == 1500], 
       paired = F)$p.value


table(aaal_pancan_gtex$Samp_nu)
ggplot(data = aaal_pancan_gtex, aes(x = Samp_nu, y = Score, fill = Method))+
  geom_boxplot(color = "grey30", outlier.shape = NA)+
  scale_fill_brewer(palette = "Pastel1") +
  ggtitle("Breast-related eQTL ranking")+
  ylab("Precision")+
  xlab("# Predictions")+
  theme_minimal() + 
  theme(axis.title.x=element_text(colour="black", size = 11, face = "bold"), 
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        #panel.grid.major = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.y=element_line(),
        axis.ticks.x=element_line(),
        #axis.ticks.length.y = el,
        legend.title=element_text(colour="black", size = 10, face = "bold"),
        legend.justification=c(1,1), 
        legend.position=c(1,1),
        axis.line = element_line( color = "grey", size=0.2),
        legend.text = element_text(colour="black", size = 10
                                   #, face = "bold"
        ),
        legend.background = element_rect(fill="gray90", size=.5, linetype=NULL),
        axis.title.y = element_text(colour="black", size = 11, face = "bold"), 
        plot.title = element_text(colour="black", size = 12, face = "bold", hjust = .5, vjust = .5)
        
  ) + 
  geom_vline(xintercept = seq(1.5,20,1), color = "grey", size=0.2)
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
aa_somatic_ornot_pc <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
                                              rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
                                              TP = aaTP, 
                                              sampling_range=seq(500,5000, 500), #c(seq(10,400,10), seq(500,3000, 100), seq(3000,31000,1000)),
                                              do_boxplot=F,
                                              percent_in_sample = T)


# aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_somatic_ornot_pc$Random_scores) * 2)
# #colnames(aa_allmat) <- aasample_size
# aac <- 1
# for(i in seq(1, ncol(aa_allmat), 2)){
#   aa_allmat[, i] <- aa_somatic_ornot_pc$Random_scores[, aac]
#   aa_allmat[, i+1] <- t(aa_somatic_ornot_pc$Model_scores)[, aac]
#   aac <- aac + 1
# }
# 
# aas <- sort(c(seq(1,100, 6),seq(2,100, 6)))
# 
# boxplot.matrix(aa_allmat[,aas], col = rep(c(2,3), length(aas)/2),
#                xaxt = "n", xlab = "#predictions", ylab = "TPR", 
#                main = "somatic ranking")
# axis(side = 1, at = seq(1.5,ncol(aa_allmat[,aas]), 2), 
#      labels = seq(100,5000, 100)[seq(2,100, 6)/2],#c(seq(10,300,10), seq(400,2000, 100), seq(3000,33000,1000)),
#      las = 2)
# abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = "grey", lwd = 0.7, lty = 3)
# 
# 
# 
# par(mfrow = c(1,1), mar = c(6,6,4,4))
# boxplot.matrix(aa_allmat[,1:100], col = rep(c(2,3), ncol(aa_somatic_ornot$Random_scores)),
#                xaxt = "n", xlab = "#predictions", ylab = "TPR", 
#                main = "breast cancer somatic from cosmic")
# axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
#      labels = seq(10,10000, 100),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)), 
#      las = 2)
# abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)
# legend(x = "topleft",legend = c("location-based", "model-based"), fill = c(2,3), cex = 0.6)
# simulate random
aa <- numeric(length = (nu_uniq_common_snp + (nu_somatic_breast_cancer - nu_somatic_common_overlap)))
aa[1:nu_somatic_breast_cancer] <- 1
aa_samp_size <- seq(500,5000, 500)
aa_score <- numeric(length(aa_samp_size)*244)


aa_cnt <- 1
for(i in 1:length(aa_samp_size)){
  
  for(j in 1:244){
    print("i")
    print(i)
    print("j")
    print(j)
    aasamp <- sample(x = aa, size = aa_samp_size[i], replace = F)
    aa_score[aa_cnt] <- sum(aasamp)/aa_samp_size[i]
    aa_cnt <- aa_cnt + 1
  }
}


aa1 <- data.frame(Score=as.numeric(aa_somatic_ornot_pc$Random_scores), 
                  Samp_nu = factor(rep(colnames(aa_somatic_ornot_pc$Random_scores),
                                       each = nrow(aa_somatic_ornot_pc$Random_scores)),
                                   levels = colnames(aa_somatic_ornot_pc$Random_scores)), 
                  Model_nu = as.numeric(rep(c(1:nrow(aa_somatic_ornot_pc$Random_scores)),
                                            ncol(aa_somatic_ornot_pc$Random_scores))),
                  Method= rep("Location-based", length(as.numeric(aa_somatic_ornot_pc$Random_scores))))

aa2 <- data.frame(Score=as.numeric(t(aa_somatic_ornot_pc$Model_scores)), 
                  Samp_nu = factor(as.numeric(rep(colnames(t(aa_somatic_ornot_pc$Model_scores)),
                                                  each = nrow(t(aa_somatic_ornot_pc$Model_scores)))),
                                   levels = colnames(t(aa_somatic_ornot_pc$Model_scores))), 
                  Model_nu = as.numeric(rep(c(1:nrow(t(aa_somatic_ornot_pc$Model_scores))),
                                            ncol(t(aa_somatic_ornot_pc$Model_scores)))),
                  Method= rep("Model-based", length(as.numeric(t(aa_somatic_ornot_pc$Model_scores)))))


aa3 <- data.frame(Score=aa_score, 
                  Samp_nu = factor(rep(as.character(aa_samp_size),
                                       each = 244),
                                   levels = as.character(aa_samp_size)), 
                  Model_nu = as.numeric(rep(c(1:244),
                                            length(aa_samp_size))),
                  Method= rep("Random", length(aa_score)))

aaal_somatic <- rbind(aa2, aa1, aa3)

t.test(aaal_somatic$Score[aaal_somatic$Method == "Model-based" & aaal_somatic$Samp_nu == 1500], 
       aaal_somatic$Score[aaal_somatic$Method == "Random" & aaal_somatic$Samp_nu == 1500], 
       paired = F)
t.test(aaal_somatic$Score[aaal_somatic$Method == "Random" & aaal_somatic$Samp_nu == 1500], 
       aaal_somatic$Score[aaal_somatic$Method == "Location-based" & aaal_somatic$Samp_nu == 1500], 
       paired = F)

ggplot(data = aaal_somatic, aes(x = Samp_nu, y = Score, fill = Method))+
  geom_boxplot(show.legend = T, color = "grey30", outlier.shape = NA)+
  scale_fill_brewer(palette = "Pastel1") +
  ggtitle("Breast Cancer somatic variant ranking")+
  ylab("Precision")+
  xlab("# Predictions")+
  theme_minimal() + 
  theme(axis.title.x=element_text(colour="black", size = 11, face = "bold"), 
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        #panel.grid.major = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.y=element_line(),
        axis.ticks.x=element_line(),
        #axis.ticks.length.y = el,
        legend.title=element_text(colour="black", size = 10, face = "bold"),
        legend.justification=c(1,1), 
        legend.position=c(1,1),
        axis.line = element_line( color = "grey", size=0.2),
        legend.text = element_text(colour="black", size = 10
                                   #, face = "bold"
        ),
        legend.background = element_rect(fill="gray90", size=.5, linetype=NULL),
        axis.title.y = element_text(colour="black", size = 11, face = "bold"), 
        plot.title = element_text(colour="black", size = 12, face = "bold", hjust = .5, vjust = .5)
        
  ) + 
  geom_vline(xintercept = seq(1.5,20,1), color = "grey", size=0.2)



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
aa_somaticFathmm_ornot_pc <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
                                                    rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
                                                    TP = aaTP, 
                                                    sampling_range=seq(50,500, 50), #c(seq(10,400,10), seq(500,3000, 100), seq(3000,31000,1000)),
                                                    do_boxplot=F, 
                                                    percent_in_sample = T)

# 
# aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_somaticFathmm_ornot$Random_scores) * 2)
# #colnames(aa_allmat) <- aasample_size
# aac <- 1
# for(i in seq(1, ncol(aa_allmat), 2)){
#   aa_allmat[, i] <- aa_somaticFathmm_ornot$Random_scores[, aac]
#   aa_allmat[, i+1] <- t(aa_somaticFathmm_ornot$Model_scores)[, aac]
#   aac <- aac + 1
# }
# par(mfrow = c(1,1), mar = c(6,6,4,4))
# boxplot.matrix(aa_allmat[,1:20], col = rep(c(2,3), ncol(aa_somaticFathmm_ornot$Random_scores)),
#                xaxt = "n", xlab = "#predictions", ylab = "TPR", 
#                main = "FATHMM significant among somatic")
# axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
#      labels = seq(10,1315, 40),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)), 
#      las = 2)
# abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)
# legend(x = "topleft",legend = c("location-based", "model-based"), fill = c(2,3), cex = 0.6)
aa <- numeric(length = nu_somatic_breast_cancer_fathmmScored)
aa[1:nu_somatic_breast_cancer_fathmmScored_sig7] <- 1
aa_samp_size <- seq(50,500, 50)
aa_score <- numeric(length(aa_samp_size)*244)
aa_cnt <- 1
for(i in 1:length(aa_samp_size)){
  
  for(j in 1:244){
    print("i")
    print(i)
    print("j")
    print(j)
    aasamp <- sample(x = aa, size = aa_samp_size[i], replace = F)
    aa_score[aa_cnt] <- sum(aasamp)/aa_samp_size[i]
    aa_cnt <- aa_cnt + 1
  }
}

aa1 <- data.frame(Score=as.numeric(aa_somaticFathmm_ornot_pc$Random_scores), 
                  Samp_nu = factor(rep(colnames(aa_somaticFathmm_ornot_pc$Random_scores),
                                       each = nrow(aa_somaticFathmm_ornot_pc$Random_scores)),
                                   levels = colnames(aa_somaticFathmm_ornot_pc$Random_scores)), 
                  Model_nu = as.numeric(rep(c(1:nrow(aa_somaticFathmm_ornot_pc$Random_scores)),
                                            ncol(aa_somaticFathmm_ornot_pc$Random_scores))),
                  Method= rep("Location-based", length(as.numeric(aa_somaticFathmm_ornot_pc$Random_scores))))

aa2 <- data.frame(Score=as.numeric(t(aa_somaticFathmm_ornot_pc$Model_scores)), 
                  Samp_nu = factor(as.numeric(rep(colnames(t(aa_somaticFathmm_ornot_pc$Model_scores)),
                                                  each = nrow(t(aa_somaticFathmm_ornot_pc$Model_scores)))),
                                   levels = colnames(t(aa_somaticFathmm_ornot_pc$Model_scores))), 
                  Model_nu = as.numeric(rep(c(1:nrow(t(aa_somaticFathmm_ornot_pc$Model_scores))),
                                            ncol(t(aa_somaticFathmm_ornot_pc$Model_scores)))),
                  Method= rep("Model-based", length(as.numeric(t(aa_somaticFathmm_ornot_pc$Model_scores)))))
aa3 <- data.frame(Score=aa_score, 
                  Samp_nu = factor(rep(as.character(aa_samp_size),
                                       each = 244),
                                   levels = as.character(aa_samp_size)), 
                  Model_nu = as.numeric(rep(c(1:244),
                                            length(aa_samp_size))),
                  Method= rep("Random", length(aa_score)))
aaal_fathmm <- rbind(aa2, aa1, aa3)

aacnnnt <- 300
t.test(aaal_fathmm$Score[aaal_fathmm$Method == "Model-based" & aaal_fathmm$Samp_nu == aacnnnt], 
       aaal_fathmm$Score[aaal_fathmm$Method == "Location-based" & aaal_fathmm$Samp_nu == aacnnnt], 
       paired = F)$p.value

t.test(aaal_fathmm$Score[aaal_fathmm$Method == "Random" & aaal_fathmm$Samp_nu == 150], 
       aaal_fathmm$Score[aaal_fathmm$Method == "Location-based" & aaal_fathmm$Samp_nu == 150], 
       paired = F)
library(ggsignif)

ggplot(data = aaal_fathmm, aes(x = Samp_nu, y = Score, fill = Method))+
  geom_boxplot(color = "grey30", outlier.shape = NA)+
 #scale_fill_manual(values=c("blue", "green"))+
  #scale_fill_grey()+
  scale_fill_brewer(palette = "Pastel1") +
  ggtitle("FATHMM significant among somatic")+
  ylab("Precision")+
  xlab("# Predictions")+
  theme_minimal() + 
  theme(axis.title.x=element_text(colour="black", size = 11, face = "bold"), 
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        #panel.grid.major = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.y=element_line(),
        axis.ticks.x=element_line(),
        #axis.ticks.length.y = el,
        legend.title=element_text(colour="black", size = 10, face = "bold"),
        legend.justification=c(1,1), 
        legend.position=c(1,1),
        axis.line = element_line( color = "grey", size=0.2),
        legend.text = element_text(colour="black", size = 10
                                   #, face = "bold"
        ),
        legend.background = element_rect(fill="gray90", size=.5, linetype=NULL),
        axis.title.y = element_text(colour="black", size = 11, face = "bold"), 
        plot.title = element_text(colour="black", size = 12, face = "bold", hjust = .5, vjust = .5)
        
  ) + 
  geom_vline(xintercept = seq(1.5,20,1), color = "grey", size=0.2)

######################################################## Clinical pathogenic among clinical + common

# aalab_orig <- ER_SNPs_merged_allbypos_modified$is_clinical_pathogenic
# #aa_som_keep <- which(ER_SNPs_merged_allbypos_modified$is_somatic)
# aa_cli_rem <- which(ER_SNPs_merged_allbypos_modified$is_clinical)
# aa_common_keep <-  which(ER_SNPs_merged_allbypos_modified$is_common)
# 
# 
# aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change), "_"), "[[", 2)))
# Merged_percentile_change_noDUP <- Merged_percentile_change[!duplicated(aasap),]
# aasap <- as.numeric(unlist(lapply(strsplit(rownames(Merged_percentile_change_noDUP), "_"), "[[", 2)))
# aaMerged_percentile_change_noDUP <- Merged_percentile_change_noDUP[ aasap %in% union(aa_cli_rem,aa_common_keep),]
# aasap <- aasap[aasap %in% union(aa_cli_rem,aa_common_keep)]
# aaTP <- aalab_orig[aasap]
# aaTP[is.na(aaTP)] <- F
# aa_clinical_patho_amongAll_ornot <- compare_ranking_to_random(all_percentile = aaMerged_percentile_change_noDUP,
#                                                               rank_among_index= c(1:nrow(aaMerged_percentile_change_noDUP)),
#                                                               TP = aaTP, 
#                                                               sampling_range=seq(10,5000, 100), #c(seq(10,400,10), seq(500,3000, 100), seq(3000,31000,1000)),
#                                                               do_boxplot=T)
# 
# 
# aa_allmat <- matrix(nrow = 244, ncol = ncol(aa_clinical_patho_amongAll_ornot$Random_scores) * 2)
# #colnames(aa_allmat) <- aasample_size
# aac <- 1
# for(i in seq(1, ncol(aa_allmat), 2)){
#   aa_allmat[, i] <- aa_clinical_patho_amongAll_ornot$Random_scores[, aac]
#   aa_allmat[, i+1] <- t(aa_clinical_patho_amongAll_ornot$Model_scores)[, aac]
#   aac <- aac + 1
# }
# par(mfrow = c(1,1), mar = c(6,6,4,4))
# 
# 
# aas <- sort(c(seq(1,100, 6),seq(2,100, 6)))
# 
# boxplot.matrix(aa_allmat[,aas], col = rep(c(2,3), length(aas)/2),
#                xaxt = "n", xlab = "#predictions", ylab = "TPR", 
#                main = "clinvar ranking")
# axis(side = 1, at = seq(1.5,ncol(aa_allmat[,aas]), 2), 
#      labels = seq(10,5000, 100)[seq(2,100, 6)/2],#c(seq(10,300,10), seq(400,2000, 100), seq(3000,33000,1000)),
#      las = 2)
# abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = "grey", lwd = 0.7, lty = 3)
# 
# 
# boxplot.matrix(aa_allmat[,1:50], col = rep(c(2,3), ncol(aa_clinical_patho_amongAll_ornot$Random_scores)),
#                xaxt = "n", xlab = "#predictions", ylab = "TPR", 
#                main = "clinvar patho among all")
# axis(side = 1, at = seq(1.5,ncol(aa_allmat), 2), 
#      labels = seq(10,5000, 100),#c(seq(10,400,10), seq(500,3000, 100), seq(3000,33000,1000)), 
#      las = 2)
# abline(v = seq(0.5, ncol(aa_allmat)+1, 2), col = 4)
# legend(x = "topleft",legend = c("location-based", "model-based"), fill = c(2,3), cex = 0.6)

######################################################## ######################################################## 
######################################################## ######################################################## 
# Visualizing an enhancer for figure 6
# done in ER_eQTL.R script: starting at line 8894


######################################################## ######################################################## 
######################################################## ######################################################## 
# supplementary figure for parameters

aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas1[1:150]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[18]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$annot)[aas2[1:150]]
aaadd3 <- union(aaadd1, aaadd2)

heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$binding[aaadd3, ]),
          Rowv = T, 
          Colv = T,
          dendrogram = 'both',
          trace = 'none',
          main = "Binding parameter values (log10)")

heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$alpha[aaadd3, ]),
          Rowv = T,
          Colv = T,
          dendrogram = 'both', 
          trace = 'none',
          main = "Activation parameter values (log10)")

heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$coop[aaadd3, ]), 
          Rowv = T, 
          Colv = T,
         dendrogram = 'both',
         trace = 'none',
         main = "Coop parameter values (log10)", 
         margins = c(10,10))

heatmap.2(log10(cbind(E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$alpha[aaadd3, ],
                      E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$binding[aaadd3, ],
                      E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$coop[aaadd3, ])), 
          Rowv = T, 
          Colv = T,
          dendrogram = 'none',
          trace = 'none',
          main = "all parameter values (log10)", 
          margins = c(10,10))

aaalpha <- E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$alpha[aaadd3, ]
colnames(aaalpha) <- paste0(colnames(aaalpha), "_alpha")
aabind <- E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$binding[aaadd3, ]
colnames(aabind) <- paste0(colnames(aabind), "_binding")
aaasam <- cbind(aaalpha,
                aabind,
                E_RNA_GEMSTAT_Ensemble_Parlist[[18]]$coop[aaadd3, ])
rownames(aaasam) <- NULL

png("~/Desktop/ER_paper/fig3f.png",         
    width = 9*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 14)        #  font size

aa_col <- colorspace::diverge_hcl(n = 30)
heatmap.2(log10(aaasam), 
          Rowv = T, 
          Colv = F,
          dendrogram = 'none',
          trace = 'none',
          main = "",
          margins = c(9,1), col = aa_col)
dev.off()





