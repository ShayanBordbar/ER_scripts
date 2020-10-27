# Linear model for finding the enhancers to work with among all the chopped enhancers of a gene.
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
########################################################################################################################
########################################                       #########################################################
########################################       LIBRARIES       #########################################################
########################################                       #########################################################
########################################################################################################################
library(sigmoid)
library(stringdist)
library(gplots)
library(nloptr)
########################################################################################################################
########################################                       #########################################################
########################################       Functions       #########################################################
########################################                       #########################################################
########################################################################################################################
# all functions in this script:

# 1. createOptimDataset
# 2. obj_func_linear_minEnh
# 3. objFuncEval
# 4. obj_func_chosen_enh
# 5. SimAnnAcceptProb
# 6. SimAnnObjOptim
# 7. SimAnnStepJobConstructor
# 8. DAGmanConstructor
# 9. GreedyObjOptim
# 10. GreedyJobConstructor
# 11. confusion_constructor
# 12. CreateRandomPrediction
# 13. tcol
# 14. ggcols
# 15. Enhancer_vote_barplot
# 16. PerformanceHeattMap_General
# 17. SilicoConcenOptim
# 18. insertOnes
# 19. Enhancer_Score_plot
# 20. TF_KD_evaluator
########################################################################################################################
createOptimDataset <-function(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                              TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun,
                              ChoppedEnhancerScores){
  # This function takes in the 
  # GeneExpMat : gene expression matrix which contains one value (in {-1, 0, 1, NA}) for each
  # (Treatment vs control) condition. one row per gene, one column per condition
  # TFExpMat : which contains one value (in {0, 0.33, 0.66, 1}) for each condition (one for control, one for treatment). one row per TF, one column per condition.
  # ChoppedEnhancerScores : is a named list with length equal to number of genes. each entry 
  # contains a matrix of affinity values of all enhancers corresponding to that gene for each TF.
  
  # outputs the dataset used for training the objective function, dataset is a list containg two entries, the first one 
  # is a matrix which contains one row per gene_(Treatment vs control)condition  pair (for non NA values of the gene). columns correspond to TF
  # expression value in that condition, and the genes observed expression value
  # second entry of the output list is ChoppedEnhancerScores sorted to be in the same order as the expression matrix
  
  output_mat <- matrix(nrow = (nrow(GeneExpMat) * ncol(GeneExpMat)), ncol = (2 * nrow(TFExpMat) + 1))
  rownames(output_mat) <- character(nrow(output_mat))
  colnames(output_mat) <- c(paste(rownames(TFExpMat), "C", sep = "_"),
                            paste(rownames(TFExpMat), "T", sep = "_"),
                            "Observed")
  
  for(cur_cond in 1:ncol(GeneExpMat)){
    cur_rows <- c(1:nrow(GeneExpMat)) + ((cur_cond - 1) * nrow(GeneExpMat))
   
    output_mat[cur_rows, 1:(2 * nrow(TFExpMat))] <-  matrix(rep(c(TFExpMat[, (2*cur_cond - 1)],
                                                                  TFExpMat[, (2*cur_cond)]),
                                                                length(cur_rows)),
                                                            nrow = length(cur_rows),
                                                            byrow = T)
    output_mat[cur_rows, (2 * nrow(TFExpMat) + 1)] <- GeneExpMat[, cur_cond]
    rownames(output_mat)[cur_rows] <- paste(rownames(GeneExpMat), colnames(GeneExpMat)[cur_cond], sep = ".")
  }
  #remove NAs
  output_mat <- output_mat[!is.na(output_mat[, ncol(output_mat)]), ]
  
  srt_ind <- match(names(ChoppedEnhancerScores), rownames(GeneExpMat))
  ChoppedEnhancerScores_sorted <- ChoppedEnhancerScores[srt_ind]
  return(list(gene_conditon_mat = output_mat,
              Affinity_scores = ChoppedEnhancerScores_sorted))
}
########################################################################################################################
########################################################################################################################
# example
aa <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                         TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun,
                         ChoppedEnhancerScores = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score) 

########################################################################################################################
########################################################################################################################
obj_func_linear_minEnh <- function(my_par, my_data){
  # this is the objective function which takes data processed in "createOptimDataset" and a set of parameters and outputs a scalar
  # my_data is the output of "createOptimDataset"
  # my_par is a numeric vector of length equal to "number of TFs plus one": one parameter per TF plus one intercept
  
  nu_TFs <- (ncol(my_data$gene_conditon_mat) - 1)/2
  nu_genes <- length(my_data$Affinity_scores)
  names(my_par) <- c(colnames(my_data$gene_conditon_mat)[1:nu_TFs], "intercept")
  # coefficient multiplied to the exponential in logSumExp
  exp_coeff <- 5
  pseudoExp <- 10e-6
  #names of the genes in each row of dataset:
  data_gene_name <- unlist(lapply(strsplit(rownames(my_data$gene_conditon_mat), split="\\."), "[[", 1))
  data_cond_name <- unlist(lapply(strsplit(rownames(my_data$gene_conditon_mat), split="\\."), "[[", 2))
  # Create a vector to hold gene losses
  gene_losses <- numeric(nu_genes)
  # loop over genes
  for(cur_gene in 1:nu_genes){
    nu_enhancers <- nrow(my_data$Affinity_scores[[cur_gene]])
    cur_gene_losses <- numeric(nu_enhancers)
    cur_gene_inds <- which(data_gene_name %in% names(my_data$Affinity_scores)[cur_gene])
    for(cur_enh in 1:nu_enhancers){
      cur_loss <- 0
      for(cur_cond in 1:length(cur_gene_inds)){
        # Current enhancer's TF affinity
        TF_aff <- my_data$Affinity_scores[[cur_gene]][cur_enh, ]
        # Current condition's control and treatment TF expressions times the affinity, concatenated with 1 for intercept 
        TF_C <- c((my_data$gene_conditon_mat[cur_gene_inds[cur_cond], 1:nu_TFs] * TF_aff), 1)
        TF_T <- c((my_data$gene_conditon_mat[cur_gene_inds[cur_cond], (nu_TFs+1):(2*nu_TFs)] * TF_aff), 1)
        # print("TF_T:")
        # print(TF_T)
        # print("TF_C :")
        # print(TF_C)
        # print("my Par:")
        # print(my_par)
         #current prediction
        #cur_pred <- 2 * sigmoid(log((relu(my_par %*% TF_T) + pseudoExp) / (relu(my_par %*% TF_C) + pseudoExp))) - 1
        
        F_treat <- relu(my_par %*% TF_T) + pseudoExp
        F_control <- relu(my_par %*% TF_C) + pseudoExp
        cur_pred <- (F_treat - F_control) / (F_treat + F_control)
        # print("#################################  ###############################")
        # print("my_par %*% TF_T :")
        # print(my_par %*% TF_T)
        # print("my_par %*% TF_C :")
        # print(my_par %*% TF_C)
        # print("#################################  ###############################")
        # print("#################################  ###############################")
        # print(paste("current prediction for", "condition", cur_cond, "/", length(cur_gene_inds), "enhancer",  cur_enh, "/", nu_enhancers, "gene", cur_gene, "/", nu_genes, ":" ))
        # print(cur_pred)
        # print("#################################  ###############################")
        #current observation
        cur_obs <- my_data$gene_conditon_mat[cur_gene_inds[cur_cond], (2*nu_TFs + 1)]
        cur_loss_cond <- (cur_pred - cur_obs)^2
        cur_loss <- cur_loss + cur_loss_cond
      }# end of loop over conditions
      cur_gene_losses[cur_enh] <- cur_loss
    }# end of loop over enhancers
    # print("current gene loss:")
    # print(cur_gene_losses)
    # max(-exp_coeff * cur_gene_losses)  # log sum exp trick
    gene_losses[cur_gene] <-  -log(sum(exp(-exp_coeff * cur_gene_losses)))
    # print("processed:")
    # print(gene_losses[cur_gene])
  }# end of loop over genes
  return(sum(gene_losses))
}


########################################################################################################################
########################################################################################################################
#example
aa <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                         TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun,
                         ChoppedEnhancerScores = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score) 
# sink("obj_func_linear_minEnh_my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52.txt")
(aaresult2 <- optim(par = runif(n = ((ncol(aa$gene_conditon_mat) - 1)/2 + 1), min = -10, max = 10), fn = obj_func_linear_minEnh, my_data = aa))
# sink()
aaaFirstOptim <- aaresult2
########################################################################################################################
########################################################################################################################
objFuncEval <- function(my_optimized_par,
                        my_data,
                        input_gene_expMat,
                        round_to_int = F,
                        chosen_enh = integer(0),
                        log_sum_exp=F,
                        fixed_threshold_value,
                        intercept_pergene=F,
                        constrained=F,
                        sigmoid_on_conditionLoss = F,
                        consider_label_weights = F){
  # This function takes the inputs to objective function and the 
  # optimized parameters and returns the prediction matrix (row per gene, column per condition),
  # and the index of enhancer used for each gene
  # round_to_int : boolean, if True rounds predictions to closest integer
  # chosen_enh : an integer vector containing one integer per gene which shows the index of the used enhancer within the enhancers associated with that gene.
  # log_sum_exp : if True uses the log sum exp approach to compute the loss, otherwise uses the error of the best enhacner (smallest error)
  # fixed_threshold_value : numeric vector of length (length(unique_labels) - 1), has the thresholds to be used for evaluation. USED in case of predictions
  # intercept_pergene : if True the model has one intercept parameter per gene.
  # constrained : if Ture Constrained optimization has been used and there is no need to apply relU or pseudoExp
  # sigmoid_on_conditionLoss : logical if True, each condition's loss will be subject to a sigmoid function to not
  #  care about losses greater than certain amount.
  # consider_label_weights : logical if True it considers the label weights in evaluating the errors
  nu_TFs <- (ncol(my_data$gene_conditon_mat) - 1)/2
  nu_genes <- length(my_data$Affinity_scores)
  nu_intercepts <- ifelse(intercept_pergene, nu_genes, 1)
  #names of the genes in each row of dataset:
  data_gene_name <- unlist(lapply(strsplit(rownames(my_data$gene_conditon_mat), split="\\."), "[[", 1))
  data_gene_name_unique <- unique(data_gene_name)
  #name of the conditions
  data_cond_name <- unlist(lapply(strsplit(rownames(my_data$gene_conditon_mat), split="\\."), "[[", 2))
  
  if(intercept_pergene){
    names(my_optimized_par) <- c(colnames(my_data$gene_conditon_mat)[1:nu_TFs], paste("intercept", data_gene_name_unique, sep = "_"))
  }else{
    names(my_optimized_par) <- c(colnames(my_data$gene_conditon_mat)[1:nu_TFs], "intercept")
  }
  
  nu_conds_pergene <- rowSums(!(t(apply(input_gene_expMat, MARGIN = 1, FUN = is.na))))
  #
  if (length(chosen_enh) > 0){
    aa <- my_data
    for(i in 1:length(aa$Affinity_scores)){
      aa$Affinity_scores[[i]] <- matrix(nrow = 1, ncol = ncol(my_data$Affinity_scores[[i]]))
      aa$Affinity_scores[[i]][1, ] <- my_data$Affinity_scores[[i]][chosen_enh[i], ]
    }
    my_data <- aa
  }
  # to store the chosen enhancer
  Chosen_Enh <- integer(nu_genes)
  
  names(Chosen_Enh) <- names(my_data$Affinity_scores)
  # to store the predicted expression
  Final_exp <- numeric(nrow(my_data$gene_conditon_mat))
  # coefficient multiplied to the exponential in sigmoids of loss
  exp_coeff1 <- 1.2
  exp_coeff2 <- 20
  # For avoiding zero in denom in LFC
  pseudoExp <- 10e-6
  # Create a vector to hold gene losses
  gene_losses <- numeric(nu_genes)
  gene_losses_perEnh <- list()
  gene_losses_perEnh_round <- list()
  gene_round_acc_perEnh <- list()
  names(gene_losses) <- names(my_data$Affinity_scores)
  
  loss_weight <- numeric(nrow(my_data$gene_conditon_mat))
  if(consider_label_weights){
    balance_table <- table(my_data$gene_conditon_mat[, ncol(my_data$gene_conditon_mat)])
    for(i in 1:length(balance_table)){
      loss_weight[my_data$gene_conditon_mat[, ncol(my_data$gene_conditon_mat)] %in% as.numeric(names(balance_table)[i])] <- max(balance_table)/balance_table[i]
    }
  }else{
    loss_weight[] <- 1
  }
  # loop over genes
  for(cur_gene in 1:nu_genes){
    nu_enhancers <- nrow(my_data$Affinity_scores[[cur_gene]])
    cur_gene_losses <- numeric(nu_enhancers)
    cur_gene_losses_round <- numeric(nu_enhancers)
    cur_gene_round_acc <- numeric(nu_enhancers)
    cur_gene_inds <- which(data_gene_name %in% names(my_data$Affinity_scores)[cur_gene])
    #to store the current predicted expressions
    cur_exp_mat <- matrix(nrow = length(cur_gene_inds), ncol = nu_enhancers)
    cur_exp_mat_round <- matrix(nrow = length(cur_gene_inds), ncol = nu_enhancers)
    
    parameter_extention <- rep(0, nu_intercepts)
    if(intercept_pergene){
      parameter_extention[cur_gene] <- 1
    }else{
      parameter_extention <- 1
    }
    
    for(cur_enh in 1:nu_enhancers){
      cur_loss <- 0
      cur_loss_round <- 0
      cur_round_acc <- 0
      # Current enhancer's TF affinity
      TF_aff <- my_data$Affinity_scores[[cur_gene]][cur_enh, ]
      for(cur_cond in 1:length(cur_gene_inds)){
        # Current condition's control and treatment TF expressions times the affinity, concatenated with 1 or (a vector of length nu_genes) for intercept 
        TF_C <- c((my_data$gene_conditon_mat[cur_gene_inds[cur_cond], 1:nu_TFs] * TF_aff), parameter_extention)
        TF_T <- c((my_data$gene_conditon_mat[cur_gene_inds[cur_cond], (nu_TFs+1):(2*nu_TFs)] * TF_aff), parameter_extention)
        #current prediction
        #cur_pred <- 2 * sigmoid(log((relu(my_optimized_par %*% TF_T) + pseudoExp) / (relu(my_optimized_par %*% TF_C) + pseudoExp))) - 1
        if(constrained){
          F_treat <- my_optimized_par %*% TF_T 
          F_control <- my_optimized_par %*% TF_C
        }else{
          F_treat <- relu(my_optimized_par %*% TF_T) + pseudoExp
          F_control <- relu(my_optimized_par %*% TF_C) + pseudoExp
        }
        pow_ratio <- (F_control/ F_treat)^(exp_coeff1/log(2))
        cur_pred <- (2 / (1 + pow_ratio)) - 1
        #cur_pred <- (F_treat - F_control) / (F_treat + F_control)
        
        cur_pred_round <- round(cur_pred)
        # cur_pred <- 2 * sigmoid(log((my_optimized_par %*% TF_T)/(my_optimized_par %*% TF_C))) - 1
        # print("#################################  ###############################")
        # print(paste("current prediction for", "condition", cur_cond, "/", length(cur_gene_inds), "enhancer",  cur_enh, "/", nu_enhancers, "gene", cur_gene, "/", nu_genes, ":" ))
        # print(cur_pred)
        # print("#################################  ###############################")
        cur_exp_mat[cur_cond, cur_enh] <- cur_pred
        cur_exp_mat_round[cur_cond, cur_enh] <- cur_pred_round
        #current observation
        cur_obs <- my_data$gene_conditon_mat[cur_gene_inds[cur_cond], (2*nu_TFs + 1)]
       
        if(sigmoid_on_conditionLoss){
          cur_loss_sig <- 0.5 * (1/(1 + exp(-1 * exp_coeff2 *(cur_pred - cur_obs)^2)) -0.5)
          cur_loss_cond <- 0.5 * cur_loss_sig * loss_weight[cur_gene_inds[cur_cond]]
        }else{
          cur_loss_cond <- 0.5 * (cur_pred - cur_obs)^2 * loss_weight[cur_gene_inds[cur_cond]]
        }
        cur_acc_cond_round <- (cur_pred_round == cur_obs)
        cur_loss_cond_round <- (cur_pred_round - cur_obs)^2
        cur_loss <- cur_loss + cur_loss_cond
        cur_loss_round <- cur_loss_round + cur_loss_cond_round
        cur_round_acc <- cur_round_acc + cur_acc_cond_round
      }# end of loop over conditions
      cur_gene_losses[cur_enh] <- cur_loss
      cur_gene_losses_round[cur_enh] <- cur_loss_round
      cur_gene_round_acc[cur_enh] <- cur_round_acc
    }# end of loop over enhancers
    Chosen_Enh[cur_gene] <- which.min(cur_gene_losses)
    gene_losses_perEnh[[cur_gene]] <- cur_gene_losses
    gene_losses_perEnh_round[[cur_gene]] <- cur_gene_losses_round
    gene_round_acc_perEnh[[cur_gene]] <- cur_gene_round_acc/nu_conds_pergene[cur_gene]
    # print("######  ##############   #########")
    # print(paste("current gene losses: ", cur_gene_losses))
    # print("#################################  ###############################")
    cur_gene_pred <- cur_exp_mat[, which.min(cur_gene_losses)]
    Final_exp[cur_gene_inds] <- cur_gene_pred
    if(log_sum_exp){
      # max(-exp_coeff * cur_gene_losses) is for preventing over flow and insuring the largest value will be e^0 = 1
      gene_losses[cur_gene] <- -log(sum(exp(-exp_coeff * cur_gene_losses)))
    }else{
      gene_losses[cur_gene] <- min(cur_gene_losses)
    }

    # if(round_to_int){ # round the predicted expression mat
    # 
    # }
    # print("processed:")
    # print(gene_losses[cur_gene])
  }# end of loop over genes
  names(gene_losses_perEnh) <- rownames(input_gene_expMat)
  names(gene_losses_perEnh_round) <- rownames(input_gene_expMat)
  names(gene_round_acc_perEnh) <- rownames(input_gene_expMat)
  
  #create the predicted gene expression mat
  PredictedGeneExpMat <- matrix(nrow = nrow(input_gene_expMat), ncol = ncol(input_gene_expMat))
  rownames(PredictedGeneExpMat) <- rownames(input_gene_expMat)
  colnames(PredictedGeneExpMat) <- colnames(input_gene_expMat)
  for(i in 1:nrow(PredictedGeneExpMat)){
    for(j in 1:ncol(PredictedGeneExpMat)){
      cur_gene_cond <- which(data_gene_name %in% rownames(PredictedGeneExpMat)[i] & (data_cond_name %in% colnames(PredictedGeneExpMat)[j]))
      if(length(cur_gene_cond) == 1){
        PredictedGeneExpMat[i, j] <- Final_exp[cur_gene_cond]
      }else if(length(cur_gene_cond) > 1){
        print("something is wrong, more than one matches to gene_cond")
      }
    }
  }
  if(round_to_int){ # compute a new error indicating the (number of true predictions) / (number of wrong predictions + number of True predictions) 
    ## Rounded_error_rate is not actually error but the opposite
    if(length(fixed_threshold_value) > 0){
      stopifnot(length(fixed_threshold_value) == (length(table(input_gene_expMat)) - 1))
      print(length(PredictedGeneExpMat))
      print(length(input_gene_expMat))
      Final_exp_discretized <- prediction_discretizer(prediction = PredictedGeneExpMat,
                                                      label=input_gene_expMat, 
                                                      fixed_thresh = fixed_threshold_value)
    }else{
      Final_exp_discretized <- prediction_discretizer(prediction = PredictedGeneExpMat,
                                                      label=input_gene_expMat)
    }

    
    PredictedGeneExpMat_discrete <- Final_exp_discretized$round_prediction
    Rounded_error_rate <- sum(PredictedGeneExpMat_discrete == input_gene_expMat, na.rm = T)/sum(!is.na(input_gene_expMat))
    return(list(Predicted_expression_mat=PredictedGeneExpMat_discrete,
                Enhancer_index=Chosen_Enh,
                losses = gene_losses,
                loss_perGene_perEnh = gene_losses_perEnh,
                Rounded_Error = Rounded_error_rate,
                accuracy_perGene_perEnh=gene_round_acc_perEnh,
                loss_perGene_perEnh_round = gene_losses_perEnh_round,
                Predicted_expression_mat_realValue=PredictedGeneExpMat,
                discretization_threshold=Final_exp_discretized$threshods))
  }
  
  return(list(Predicted_expression_mat=PredictedGeneExpMat,
              Enhancer_index=Chosen_Enh,
              losses = gene_losses,
              loss_perGene_perEnh = gene_losses_perEnh,
              accuracy_perGene_perEnh=gene_round_acc_perEnh,
              loss_perGene_perEnh_round = gene_losses_perEnh_round))
  
}

# simulation to find out best coeff
# aaF_control <- seq(0, 20, length.out = 400) + 1e-6
# aaF_treat <- seq(0, 20, length.out = 400) + 1e-6
# aaexp_coeff1 <- c(1, 1.2, 1.5, 1.8, 2, 3, 4, 5)
# for(k in 1:length(aaexp_coeff1)){
#   aapow_ratio <- matrix(nrow = length(aaF_control), ncol = length(aaF_treat))
#   aacur_pred <- matrix(nrow = length(aaF_control), ncol = length(aaF_treat))
#   rownames(aacur_pred) <- format(x = aaF_control, digits = 2, scientific = F) 
#   colnames(aacur_pred) <- format(x = aaF_treat, digits = 2, scientific = F) 
#   for(i in 1:length(aaF_control)){
#     for(j in 1:length(aaF_treat)){
#       aapow_ratio[i, j] <- (aaF_control[i]/ aaF_treat[j])^(aaexp_coeff1[k]/log(2))
#       aacur_pred[i, j] <- (2 / (1 + aapow_ratio[i, j])) - 1
#     }
#   }
#   heatmap.2(aacur_pred,
#             Rowv = F,
#             Colv = F, 
#             dendrogram = "none",
#             trace = "none",
#             breaks = seq(-1, 1, length.out = 200),
#             col = colorspace::diverge_hcl(199), 
#             main = paste0("aaexp_coeff1 : ", aaexp_coeff1[k]))
#   print("aaexp_coeff1")
#   print(aaexp_coeff1[k])
#   print(paste0("less than -0.5 = ", sum(aacur_pred <= -0.5) / length(aacur_pred), 
#                " between -0.5 and 0.5 = ", 
#                sum(aacur_pred > -0.5 & aacur_pred < 0.5)/length(aacur_pred), 
#                " greater than 0.5 = ", sum(aacur_pred >= 0.5)/length(aacur_pred))
#         )
#   print("###################################################")
# }



########################################################################################################################
########################################################################################################################
#example
aa <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                         TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun,
                         ChoppedEnhancerScores = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score) 

aa_start_par <- runif(n = ((ncol(aa$gene_conditon_mat) - 1)/2 + 1), min = -1000, max = 1000)
aa_st_enh <- 
aaresult2 <- optim(par = runif(n = ((ncol(aa$gene_conditon_mat) - 1)/2 + 1), min = -1000, max = 1000), fn = obj_func_chosen_enh, my_data = aa, chosen_enh_ind = )
aaresult3 <- optim(par = runif(n = ((ncol(aa$gene_conditon_mat) - 1)/2 + 1), min = -1000, max = 1000), fn = obj_func_chosen_enh, my_data = aa, chosen_enh_ind = )
aaaSecondOptim <- aaresult2
aaa <- objFuncEval(my_optimized_par = aaresult2$par, my_data = aa, input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)

plotExpression(expMat = aaa$Predicted_expression_mat, .dendrogram = "none", .Rowv = F, filename = "Prediction_minEnh_2.png")
plotExpression(expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, .dendrogram = "none", .Rowv = F, filename = "RealMat.png")


#####

########################################################################################################################
########################################################################################################################
obj_func_chosen_enh <- function(my_par, 
                                my_data, 
                                chosen_enh_ind, 
                                intercept_pergene=F,
                                constrained=F, 
                                return_grad = T, 
                                sig_over_loss=F){
  # my_data is the output of "createOptimDataset"
  # chosen_enh_ind is an integer vector of length = number of genes. contains the index of
  #    the chosen enhancer among the enhancers associated with this gene 
  # my_par is a numeric vector of length equal to "number of TFs + number of intercepts": 
  #    one parameter per TF plus one intercept (for all or per gene)
  
  nu_TFs <- (ncol(my_data$gene_conditon_mat) - 1)/2
  nu_genes <- length(my_data$Affinity_scores)
  nu_intercepts <- ifelse(intercept_pergene, nu_genes, 1)
  
  # names of the genes in each row of dataset:
  data_gene_name <- unlist(lapply(strsplit(rownames(my_data$gene_conditon_mat), split="\\."), "[[", 1))
  data_gene_name_unique <- unique(data_gene_name)
  stopifnot((intercept_pergene & length(my_par) == (nu_TFs + nu_genes)) | (!intercept_pergene &  length(my_par) == (nu_TFs + 1)))
  
  if(intercept_pergene){
    names(my_par) <- c(colnames(my_data$gene_conditon_mat)[1:nu_TFs], 
                       paste("intercept", data_gene_name_unique, sep = "_"))
  }else{
    names(my_par) <- c(colnames(my_data$gene_conditon_mat)[1:nu_TFs], "intercept")
  }
  pseudoExp <- 10e-6

  # Create a vector to hold gene losses
  gene_losses <- numeric(nu_genes)
  
  #create a weight vector for observation, such that at the loss of each condition is multiplied by the corresponding weight value
  loss_weight <- numeric(nrow(my_data$gene_conditon_mat))
  balance_table <- table(my_data$gene_conditon_mat[, ncol(my_data$gene_conditon_mat)])
  for(i in 1:length(balance_table)){
    loss_weight[my_data$gene_conditon_mat[, ncol(my_data$gene_conditon_mat)] %in% as.numeric(names(balance_table)[i])] <- max(balance_table)/balance_table[i]
  }
  # loop over genes
  for(cur_gene in 1:nu_genes){
    cur_chosen_enh <- chosen_enh_ind[cur_gene]
    cur_gene_cond_inds <- which(data_gene_name %in% names(my_data$Affinity_scores)[cur_gene])
    # Current enhancer's TF affinity
    TF_aff <- my_data$Affinity_scores[[cur_gene]][cur_chosen_enh, ]
    cur_loss <- 0
    for(cur_cond in 1:length(cur_gene_cond_inds)){
      
      parameter_extention <- rep(0, nu_intercepts)
      if(intercept_pergene){
        parameter_extention[cur_gene] <- 1
      }else{
        parameter_extention <- 1
      }
      # Current condition's control and treatment TF expressions times the affinity, concatenated with 1 for intercept 
      TF_C <- c((my_data$gene_conditon_mat[cur_gene_cond_inds[cur_cond], 1:nu_TFs] * TF_aff), parameter_extention)
      TF_T <- c((my_data$gene_conditon_mat[cur_gene_cond_inds[cur_cond], (nu_TFs+1):(2*nu_TFs)] * TF_aff), parameter_extention)
        # print("TF_T:")
        # print(TF_T)
        # print("TF_C :")
        # print(TF_C)
        # print("my Par:")
        # print(my_par)
        #current prediction
        
      #cur_pred <- 2 * sigmoid(log((relu(my_par %*% TF_T) + pseudoExp) / (relu(my_par %*% TF_C) + pseudoExp))) - 1
      F_treat <- relu(my_par %*% TF_T) + pseudoExp
      F_control <- relu(my_par %*% TF_C) + pseudoExp
      sig_coeff <- 1.2
      pow_ratio <- (F_control/ F_treat)^(sig_coeff/log(2))
      cur_pred <- (2 / (1 + pow_ratio)) - 1
      # cur_pred <- (F_treat - F_control) / (F_treat + F_control)
        # print("#################################  ###############################")
        # print("my_par %*% TF_T :")
        # print(my_par %*% TF_T)
        # print("my_par %*% TF_C :")
        # print(my_par %*% TF_C)
        # print("#################################  ###############################")
        # print("#################################  ###############################")
        # print(paste("current prediction for", "condition", cur_cond, "/", length(cur_gene_inds), "enhancer",  cur_enh, "/", nu_enhancers, "gene", cur_gene, "/", nu_genes, ":" ))
        # print(cur_pred)
        # print("#################################  ###############################")
        #current observation
      cur_obs <- my_data$gene_conditon_mat[cur_gene_cond_inds[cur_cond], (2*nu_TFs + 1)]
      if(sig_over_loss){
        cur_loss_sig <- 0.5 * (1/(1 + exp(-20*(cur_pred - cur_obs)^2)) -0.5)
      }else{
        cur_loss_sig <- (cur_pred - cur_obs)^2
      }
      cur_loss_cond <- 0.5 * cur_loss_sig * loss_weight[cur_gene_cond_inds[cur_cond]]
      cur_loss <- cur_loss + cur_loss_cond
    }# end of loop over conditions
    gene_losses[cur_gene] <- cur_loss
    # print("current gene loss:")
    # print(cur_gene_enh_loss)
    # print("processed:")
    # print(gene_losses[cur_gene])
  }# end of loop over genes
  return(sum(gene_losses))
}

# aaobs = -1
# aapred <- numeric(200)
# aax <- numeric(200)
# for(i in c(1:200)){
#   aapred[i] <- seq(-1, 1, length.out = 200)[i]
#   aax[i] <- 0.5 * (1/(1 + exp(-20*(aaobs - aapred[i])^2)) -0.5)
# }
# plot(aapred[aapred > 0.5], aax[aapred > 0.5])
# plot(aapred, aax)
# aax[aapred <= -0.407 & aapred >= -0.55]
########################################################################################################################
########################################################################################################################

obj_func_chosen_enh_gradient <- function(my_par, 
                                         my_data,
                                         chosen_enh_ind, 
                                         intercept_pergene=F, 
                                         constrained=T,
                                         return_grad=T){
  # my_data is the output of "createOptimDataset"
  # chosen_enh_ind is an integer vector of length = number of genes. contains the index of
  #    the chosen enhancer among the enhancers associated with this gene 
  # my_par is a numeric vector of length equal to "number of TFs + ifelse(intercept_pergene, nu_genes, 1)": 
  #    one parameter per TF plus one (or nu_genes) intercept
  # intercept_pergene if True the model uses one intercept per gene
  # constrained: if True it uses constrains to keep predicted expression values positive, otherwise it uses relU
  # return_grad : if True computes and returns gradients
  nu_TFs <- (ncol(my_data$gene_conditon_mat) - 1)/2
  nu_genes <- length(my_data$Affinity_scores)
  nu_intercepts <- ifelse(intercept_pergene, nu_genes, 1)


  # names of the genes in each row of dataset:
  data_gene_name <- unlist(lapply(strsplit(rownames(my_data$gene_conditon_mat), split="\\."), "[[", 1))
  data_gene_name_unique <- unique(data_gene_name)
  stopifnot((intercept_pergene & length(my_par) == (nu_TFs + nu_genes)) | (!intercept_pergene &  length(my_par) == (nu_TFs + 1)))
  
  # if(!intercept_pergene){
  #   # print("this is not yet implemented for one intercept par for all genes")
  #   # stop()
  #   
  # }
  if(intercept_pergene){
    names(my_par) <- c(colnames(my_data$gene_conditon_mat)[1:nu_TFs], 
                       paste("intercept", data_gene_name_unique, sep = "_"))
  }else{
    names(my_par) <- c(colnames(my_data$gene_conditon_mat)[1:nu_TFs], "intercept")
  }
  pseudoExp <- 10e-6
  
  # Create a vector to hold gene losses
  gene_losses <- numeric(nu_genes)
  # Create a vector to hold gradients
  grad_holder <- numeric(nu_TFs + nu_intercepts)
  
  
  #create a weight vector for observation, such that at the loss of each condition is multiplied by the corresponding weight value
  loss_weight <- numeric(nrow(my_data$gene_conditon_mat))
  balance_table <- table(my_data$gene_conditon_mat[, ncol(my_data$gene_conditon_mat)])
  for(i in 1:length(balance_table)){
    loss_weight[my_data$gene_conditon_mat[, ncol(my_data$gene_conditon_mat)] %in% as.numeric(names(balance_table)[i])] <- max(balance_table)/balance_table[i]
  }
  # loop over genes
  for(cur_gene in 1:nu_genes){
    cur_chosen_enh <- chosen_enh_ind[cur_gene]
    cur_gene_cond_inds <- which(data_gene_name %in% names(my_data$Affinity_scores)[cur_gene])
    # Current enhancer's TF affinity
    TF_aff <- my_data$Affinity_scores[[cur_gene]][cur_chosen_enh, ]
    cur_loss <- 0
    for(cur_cond in 1:length(cur_gene_cond_inds)){
      parameter_extention <- rep(0, nu_intercepts)
      if(intercept_pergene){
        parameter_extention[cur_gene] <- 1
      }else{
        parameter_extention <- 1
      }
      # Current condition's control and treatment TF expressions times the affinity, concatenated with 1 for intercept 
      TF_C <- c((my_data$gene_conditon_mat[cur_gene_cond_inds[cur_cond], 1:nu_TFs] * TF_aff), parameter_extention)
      TF_T <- c((my_data$gene_conditon_mat[cur_gene_cond_inds[cur_cond], (nu_TFs+1):(2*nu_TFs)] * TF_aff), parameter_extention)
      # print("TF_T:")
      # print(TF_T)
      # print("TF_C :")
      # print(TF_C)
      # print("my Par:")
      # print(my_par)
      #current prediction
      
      #cur_pred <- 2 * sigmoid(log((relu(my_par %*% TF_T) + pseudoExp) / (relu(my_par %*% TF_C) + pseudoExp))) - 1
      
      # F_treat <- relu(my_par %*% TF_T) + pseudoExp
      if(constrained){
        F_treat <-   my_par %*% TF_T
        F_control <- my_par %*% TF_C
      }else{
        F_treat_raw <- my_par %*% TF_T
        F_treat <- relu(F_treat_raw) + pseudoExp
        F_control_raw <- my_par %*% TF_C
        F_control <- relu(F_control_raw) + pseudoExp
      }

      
      #cur_pred <- (F_treat - F_control) / (F_treat + F_control) # This is wrong should be changed in all cases
      pow_ratio <- (F_control/ F_treat)^(1/log(2))
      cur_pred <- (2 / (1 + pow_ratio)) - 1
      # print("#################################  ###############################")
      # print("my_par %*% TF_T :")
      # print(my_par %*% TF_T)
      # print("my_par %*% TF_C :")
      # print(my_par %*% TF_C)
      # print("#################################  ###############################")
      # print("#################################  ###############################")
      # print(paste("current prediction for", "condition", cur_cond, "/", length(cur_gene_inds), "enhancer",  cur_enh, "/", nu_enhancers, "gene", cur_gene, "/", nu_genes, ":" ))
      # print(cur_pred)
      # print("#################################  ###############################")
      #current observation
      cur_obs <- my_data$gene_conditon_mat[cur_gene_cond_inds[cur_cond], (2*nu_TFs + 1)]
      cur_loss_cond <- 0.5 * (cur_pred - cur_obs)^2 * loss_weight[cur_gene_cond_inds[cur_cond]]
      cur_loss <- cur_loss + cur_loss_cond
      
      if(return_grad){
        # Compute gradients
        # update TF parmaeter gradients
        for(cur_tf in 1:nu_TFs){
          if(constrained){
            r_Fc_r_alpha_t <- TF_C[cur_tf]
            r_Ft_r_alpha_t <- TF_T[cur_tf]
            # partial_control <- TF_C[cur_tf] * F_treat
            # partial_treat   <- TF_T[cur_tf] * F_control
          }else{
            r_Fc_r_alpha_t <- ifelse(F_control_raw >= 0, TF_C[cur_tf], 0)
            r_Ft_r_alpha_t <- ifelse(F_control_raw >= 0, TF_T[cur_tf], 0)
            # partial_control <- ifelse(F_control_raw >= 0, 2 * TF_C[cur_tf] * F_treat, 0)
            # partial_treat   <- ifelse(F_treat_raw >= 0, 2 * TF_T[cur_tf] * F_control, 0)
          }
          r_Pred_r_powrat <- -2 / ((1 + pow_ratio)^2)
          r_powrat_r_alpha_t <- (1 / log(2)) * (F_control/ F_treat)^((1/log(2)) - 1) * (r_Fc_r_alpha_t * F_treat - r_Ft_r_alpha_t * F_control) / (F_treat^2)
          r_Pred_r_alpha_t <- r_Pred_r_powrat * r_powrat_r_alpha_t
#          partial_grad <- (partial_treat - partial_control)/(1 + F_control/F_treat)^2
          
          cur_grad <- (cur_pred - cur_obs) * loss_weight[cur_gene_cond_inds[cur_cond]] * r_Pred_r_alpha_t
          grad_holder[cur_tf] <- grad_holder[cur_tf] + cur_grad
        }
        # update intercept parmaeter gradients
        if(intercept_pergene){
          print("gradient computation for one intercept per gene has not been implemented yet")
          stop()
         #  if(constrained){
         #    int_p_treat <- 1
         #    int_p_control <- 1
         #  }else{
         #    int_p_treat <-   ifelse(F_treat_raw >= 0, 1, 0)
         #    int_p_control <- ifelse(F_control_raw >= 0, 1, 0)
         #  }
         #  
         # # partial_intercept <- 2* (F_control*int_p_treat - F_treat*int_p_control)/(1 + F_control/F_treat)^2
         #  intercept_grad <- (cur_pred - cur_obs) * loss_weight[cur_gene_cond_inds[cur_cond]] * partial_intercept
         #  grad_holder[nu_TFs + cur_gene] <- grad_holder[nu_TFs + cur_gene] + intercept_grad
        }else{ # one intercept for all
          if(constrained){
            int_p_treat <- 1
            int_p_control <- 1
          }else{
            int_p_treat <-   ifelse(F_treat_raw >= 0, 1, 0)
            int_p_control <- ifelse(F_control_raw >= 0, 1, 0)
          }
          r_powrat_r_beta <- (1 / log(2)) * (F_control/ F_treat)^((1/log(2)) - 1) * nu_TFs * ((F_treat * int_p_control) - (F_control * int_p_treat)) / F_treat^2
          r_Pred_r_beta <- r_Pred_r_powrat * r_powrat_r_beta
          intercept_grad <- (cur_pred - cur_obs) * loss_weight[cur_gene_cond_inds[cur_cond]] * r_Pred_r_beta
          grad_holder[nu_TFs + 1] <- grad_holder[nu_TFs + 1] + intercept_grad
        }
      }
    }# end of loop over conditions
    gene_losses[cur_gene] <- cur_loss
    # print("current gene loss:")
    # print(cur_gene_enh_loss)
    # print("processed:")
    # print(gene_losses[cur_gene])
  }# end of loop over genes
  if(return_grad){
    #my_Res <- sum(gene_losses)
    #attr(my_Res, "gradient") <- grad_holder
    #my_Res
    return(list("objective" = sum(gene_losses), "gradient" = grad_holder))
  }else{
    return(sum(gene_losses))
  }
}

########################################################################################################################
########################################################################################################################
obj_func_chosen_enh_gradient_only <- function(my_par, 
                                         my_data,
                                         chosen_enh_ind, 
                                         intercept_pergene=F, 
                                         constrained=T,
                                         return_grad=T){
  # my_data is the output of "createOptimDataset"
  # chosen_enh_ind is an integer vector of length = number of genes. contains the index of
  #    the chosen enhancer among the enhancers associated with this gene 
  # my_par is a numeric vector of length equal to "number of TFs + ifelse(intercept_pergene, nu_genes, 1)": 
  #    one parameter per TF plus one (or nu_genes) intercept
  # intercept_pergene if True the model uses one intercept per gene
  # constrained: if True it uses constrains to keep predicted expression values positive, otherwise it uses relU
  # return_grad : if True computes and returns gradients
  nu_TFs <- (ncol(my_data$gene_conditon_mat) - 1)/2
  nu_genes <- length(my_data$Affinity_scores)
  nu_intercepts <- ifelse(intercept_pergene, nu_genes, 1)
  
  
  # names of the genes in each row of dataset:
  data_gene_name <- unlist(lapply(strsplit(rownames(my_data$gene_conditon_mat), split="\\."), "[[", 1))
  data_gene_name_unique <- unique(data_gene_name)
  stopifnot((intercept_pergene & length(my_par) == (nu_TFs + nu_genes)) | (!intercept_pergene &  length(my_par) == (nu_TFs + 1)))
  
  # if(!intercept_pergene){
  #   # print("this is not yet implemented for one intercept par for all genes")
  #   # stop()
  #   
  # }
  if(intercept_pergene){
    names(my_par) <- c(colnames(my_data$gene_conditon_mat)[1:nu_TFs], 
                       paste("intercept", data_gene_name_unique, sep = "_"))
  }else{
    names(my_par) <- c(colnames(my_data$gene_conditon_mat)[1:nu_TFs], "intercept")
  }
  pseudoExp <- 10e-6
  
  # Create a vector to hold gene losses
  gene_losses <- numeric(nu_genes)
  # Create a vector to hold gradients
  grad_holder <- numeric(nu_TFs + nu_intercepts)
  
  
  #create a weight vector for observation, such that at the loss of each condition is multiplied by the corresponding weight value
  loss_weight <- numeric(nrow(my_data$gene_conditon_mat))
  balance_table <- table(my_data$gene_conditon_mat[, ncol(my_data$gene_conditon_mat)])
  for(i in 1:length(balance_table)){
    loss_weight[my_data$gene_conditon_mat[, ncol(my_data$gene_conditon_mat)] %in% as.numeric(names(balance_table)[i])] <- max(balance_table)/balance_table[i]
  }
  # loop over genes
  for(cur_gene in 1:nu_genes){
    cur_chosen_enh <- chosen_enh_ind[cur_gene]
    cur_gene_cond_inds <- which(data_gene_name %in% names(my_data$Affinity_scores)[cur_gene])
    # Current enhancer's TF affinity
    TF_aff <- my_data$Affinity_scores[[cur_gene]][cur_chosen_enh, ]
    cur_loss <- 0
    for(cur_cond in 1:length(cur_gene_cond_inds)){
      parameter_extention <- rep(0, nu_intercepts)
      if(intercept_pergene){
        parameter_extention[cur_gene] <- 1
      }else{
        parameter_extention <- 1
      }
      # Current condition's control and treatment TF expressions times the affinity, concatenated with 1 for intercept 
      TF_C <- c((my_data$gene_conditon_mat[cur_gene_cond_inds[cur_cond], 1:nu_TFs] * TF_aff), parameter_extention)
      TF_T <- c((my_data$gene_conditon_mat[cur_gene_cond_inds[cur_cond], (nu_TFs+1):(2*nu_TFs)] * TF_aff), parameter_extention)
      # print("TF_T:")
      # print(TF_T)
      # print("TF_C :")
      # print(TF_C)
      # print("my Par:")
      # print(my_par)
      #current prediction
      
      #cur_pred <- 2 * sigmoid(log((relu(my_par %*% TF_T) + pseudoExp) / (relu(my_par %*% TF_C) + pseudoExp))) - 1
      
      # F_treat <- relu(my_par %*% TF_T) + pseudoExp
      if(constrained){
        F_treat <-   my_par %*% TF_T
        F_control <- my_par %*% TF_C
      }else{
        F_treat_raw <- my_par %*% TF_T
        F_treat <- relu(F_treat_raw) + pseudoExp
        F_control_raw <- my_par %*% TF_C
        F_control <- relu(F_control_raw) + pseudoExp
      }
      
      
      #cur_pred <- (F_treat - F_control) / (F_treat + F_control) # This is wrong should be changed in all cases
      pow_ratio <- (F_control/ F_treat)^(1/log(2))
      cur_pred <- (2 / (1 + pow_ratio)) - 1
      # print("#################################  ###############################")
      # print("my_par %*% TF_T :")
      # print(my_par %*% TF_T)
      # print("my_par %*% TF_C :")
      # print(my_par %*% TF_C)
      # print("#################################  ###############################")
      # print("#################################  ###############################")
      # print(paste("current prediction for", "condition", cur_cond, "/", length(cur_gene_inds), "enhancer",  cur_enh, "/", nu_enhancers, "gene", cur_gene, "/", nu_genes, ":" ))
      # print(cur_pred)
      # print("#################################  ###############################")
      #current observation
      cur_obs <- my_data$gene_conditon_mat[cur_gene_cond_inds[cur_cond], (2*nu_TFs + 1)]
      # cur_loss_cond <- 0.5 * (cur_pred - cur_obs)^2 * loss_weight[cur_gene_cond_inds[cur_cond]]
      # cur_loss <- cur_loss + cur_loss_cond
      
      if(return_grad){
        # Compute gradients
        # update TF parmaeter gradients
        for(cur_tf in 1:nu_TFs){
          if(constrained){
            r_Fc_r_alpha_t <- TF_C[cur_tf]
            r_Ft_r_alpha_t <- TF_T[cur_tf]
            # partial_control <- TF_C[cur_tf] * F_treat
            # partial_treat   <- TF_T[cur_tf] * F_control
          }else{
            r_Fc_r_alpha_t <- ifelse(F_control_raw >= 0, TF_C[cur_tf], 0)
            r_Ft_r_alpha_t <- ifelse(F_control_raw >= 0, TF_T[cur_tf], 0)
            # partial_control <- ifelse(F_control_raw >= 0, 2 * TF_C[cur_tf] * F_treat, 0)
            # partial_treat   <- ifelse(F_treat_raw >= 0, 2 * TF_T[cur_tf] * F_control, 0)
          }
          r_Pred_r_powrat <- -2 / ((1 + pow_ratio)^2)
          r_powrat_r_alpha_t <- (1 / log(2)) * (F_control/ F_treat)^((1/log(2)) - 1) * (r_Fc_r_alpha_t * F_treat - r_Ft_r_alpha_t * F_control) / (F_treat^2)
          r_Pred_r_alpha_t <- r_Pred_r_powrat * r_powrat_r_alpha_t
          #          partial_grad <- (partial_treat - partial_control)/(1 + F_control/F_treat)^2
          
          cur_grad <- (cur_pred - cur_obs) * loss_weight[cur_gene_cond_inds[cur_cond]] * r_Pred_r_alpha_t
          grad_holder[cur_tf] <- grad_holder[cur_tf] + cur_grad
        }
        # update intercept parmaeter gradients
        if(intercept_pergene){
          print("gradient computation for one intercept per gene has not been implemented yet")
          stop()
          #  if(constrained){
          #    int_p_treat <- 1
          #    int_p_control <- 1
          #  }else{
          #    int_p_treat <-   ifelse(F_treat_raw >= 0, 1, 0)
          #    int_p_control <- ifelse(F_control_raw >= 0, 1, 0)
          #  }
          #  
          # # partial_intercept <- 2* (F_control*int_p_treat - F_treat*int_p_control)/(1 + F_control/F_treat)^2
          #  intercept_grad <- (cur_pred - cur_obs) * loss_weight[cur_gene_cond_inds[cur_cond]] * partial_intercept
          #  grad_holder[nu_TFs + cur_gene] <- grad_holder[nu_TFs + cur_gene] + intercept_grad
        }else{ # one intercept for all
          if(constrained){
            int_p_treat <- 1
            int_p_control <- 1
          }else{
            int_p_treat <-   ifelse(F_treat_raw >= 0, 1, 0)
            int_p_control <- ifelse(F_control_raw >= 0, 1, 0)
          }
          r_powrat_r_beta <- (1 / log(2)) * (F_control/ F_treat)^((1/log(2)) - 1) * nu_TFs * ((F_treat * int_p_control) - (F_control * int_p_treat)) / F_treat^2
          r_Pred_r_beta <- r_Pred_r_powrat * r_powrat_r_beta
          intercept_grad <- (cur_pred - cur_obs) * loss_weight[cur_gene_cond_inds[cur_cond]] * r_Pred_r_beta
          grad_holder[nu_TFs + 1] <- grad_holder[nu_TFs + 1] + intercept_grad
        }
      }
    }# end of loop over conditions
    gene_losses[cur_gene] <- cur_loss
    # print("current gene loss:")
    # print(cur_gene_enh_loss)
    # print("processed:")
    # print(gene_losses[cur_gene])
  }# end of loop over genes
  if(return_grad){
    return(grad_holder)
  }else{
    print("ONLY GRADIENT!!")
  }
}
########################################################################################################################
########################################################################################################################
# function for inequality constrains
obj_func_chosen_enh_ineq <- function(my_par, my_data, chosen_enh_ind, intercept_pergene=T){
  nu_TFs <- (ncol(my_data$gene_conditon_mat) - 1)/2
  nu_genes <- length(my_data$Affinity_scores)
  data_gene_name <- unlist(lapply(strsplit(rownames(my_data$gene_conditon_mat), split="\\."), "[[", 1))
  data_gene_name_unique <- unique(data_gene_name)
  inequality_neg <- numeric(nrow(my_data$gene_conditon_mat)*2)
  ineq_count = 1
  grad_holder <- numeric((nu_TFs + nu_genes) * nrow(my_data$gene_conditon_mat) * 2)
  # loop over genes
  for(cur_gene in 1:nu_genes){
    cur_chosen_enh <- chosen_enh_ind[cur_gene]
    cur_gene_cond_inds <- which(data_gene_name %in% names(my_data$Affinity_scores)[cur_gene])
    # Current enhancer's TF affinity
    TF_aff <- my_data$Affinity_scores[[cur_gene]][cur_chosen_enh, ]
    cur_loss <- 0
    for(cur_cond in 1:length(cur_gene_cond_inds)){
      # Current condition's control and treatment TF expressions times the affinity, concatenated with 1 for intercept 
      parameter_extention <- rep(0, nu_genes)
      if(intercept_pergene){
        parameter_extention[cur_gene] <- 1
      }else{
        parameter_extention <- 1
      }
      TF_C <- c((my_data$gene_conditon_mat[cur_gene_cond_inds[cur_cond], 1:nu_TFs] * TF_aff), parameter_extention)
      TF_T <- c((my_data$gene_conditon_mat[cur_gene_cond_inds[cur_cond], (nu_TFs+1):(2*nu_TFs)] * TF_aff), parameter_extention)
      # print("TF_T:")
      # print(TF_T)
      # print("TF_C :")
      # print(TF_C)
      # print("my Par:")
      # print(my_par)
      #current prediction
      
      #cur_pred <- 2 * sigmoid(log((relu(my_par %*% TF_T) + pseudoExp) / (relu(my_par %*% TF_C) + pseudoExp))) - 1
      
      # F_treat <- relu(my_par %*% TF_T) + pseudoExp : instead of this format I will just use the simple format, but I'll add constrains to the optimization
      F_treat <- my_par %*% TF_T
      F_control <- my_par %*% TF_C
      
      
      inequality_neg[ineq_count] <- -1 * F_treat
      inequality_neg[ineq_count + 1] <- 1e-9 - F_control
      
      # print("#################################  ###############################")
      # print("my_par %*% TF_T :")
      # print(my_par %*% TF_T)
      # print("my_par %*% TF_C :")
      # print(my_par %*% TF_C)
      # print("#################################  ###############################")
      # print("#################################  ###############################")
      # print(paste("current prediction for", "condition", cur_cond, "/", length(cur_gene_inds), "enhancer",  cur_enh, "/", nu_enhancers, "gene", cur_gene, "/", nu_genes, ":" ))
      # print(cur_pred)
      # print("#################################  ###############################")

      # Compute gradients
      if(intercept_pergene){
        # update TF parmaeter gradients
        for(cur_tf in 1:nu_TFs){
         grad_holder[(ineq_count-1 * length(my_par)) + cur_tf] <- grad_holder[cur_tf] - TF_T[cur_tf]
         grad_holder[(ineq_count * length(my_par)) + cur_tf] <- grad_holder[cur_tf] - TF_C[cur_tf]
        }
        # update intercept parmaeter gradients
        grad_holder[(ineq_count-1 * length(my_par)) + nu_TFs + cur_gene] <- grad_holder[nu_TFs + cur_gene] - 1
        grad_holder[(ineq_count* length(my_par)) + nu_TFs + cur_gene] <- grad_holder[nu_TFs + cur_gene] - 1
        
        ineq_count = ineq_count + 2
      }else{
        print("gradient computation for one intercept for all has not been implemented yet")
        grad_holder[] <- 1
      }
      
    }# end of loop over conditions
    #gene_losses[cur_gene] <- cur_loss
    # print("current gene loss:")
    # print(cur_gene_enh_loss)
    # print("processed:")
    # print(gene_losses[cur_gene])
  }# end of loop over genes
  #print(grad_holder)
  return( list( "constraints"=inequality_neg, "jacobian"=grad_holder ) )
}
########################################################################################################################
########################################################################################################################
#example

# initial values
aa_start_par <- runif(n = ((ncol(ER_52_opt_simAnn_2_input$gene_conditon_mat) - 1)/2  + 52), min = -10, max = 10)

# lower and upper bounds of control
#lb <- c( 1, 1, 1, 1 )
#ub <- c( 5, 5, 5, 5 )
local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel"  = 1.0e-7 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel"  = 1.0e-7,
              "maxeval"   = 1000,
              "local_opts" = local_opts )
aa_time_grad <- proc.time()
res <- nloptr( x0=aa_start_par,
               eval_f=obj_func_chosen_enh_gradient,
               my_data = ER_52_opt_simAnn_2_input,
               chosen_enh_ind = Sim_Ann_weighted_148_restart_enhancers[1, ], 
               intercept_pergene=T, 
         #      lb=lb,
          #     ub=ub,
               eval_g_ineq=obj_func_chosen_enh_ineq,
          #     eval_g_eq=eval_g_eq,
               opts=opts)
proc.time()
print( res )

aa_time_grad2 <- proc.time()
# initial values
aa_start_par <- runif(n = ((ncol(ER_52_opt_simAnn_2_input$gene_conditon_mat) - 1)/2  + 1), min = -10, max = 10)
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel"  = 1.0e-7,
              "maxeval"   = 1000,
              "local_opts" = local_opts )
res2 <- nloptr( x0=aa_start_par,
               eval_f=obj_func_chosen_enh_gradient,
               my_data = ER_52_opt_simAnn_2_input,
               chosen_enh_ind = Sim_Ann_weighted_148_restart_enhancers[1, ],
               intercept_pergene=F, 
               constrained=F,
               return_grad = T,
               #      lb=lb,
               #     ub=ub,
               #eval_g_ineq=obj_func_chosen_enh_ineq,
               #     eval_g_eq=eval_g_eq,
               opts=opts)
proc.time()
print( res2 )
aa_grad_eval <- objFuncEval(my_optimized_par = res2$solution,
                            my_data = ER_52_opt_simAnn_2_input,
                            input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                            round_to_int = T,
                            chosen_enh = integer(0),
                            log_sum_exp=F,
                            intercept_pergene=F,
                            constrained=F)
aa_grad_eval2 <- objFuncEval(my_optimized_par = res$solution,
                            my_data = ER_52_opt_simAnn_2_input,
                            input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                            round_to_int = F,
                            chosen_enh = integer(0),
                            log_sum_exp=F,
                            intercept_pergene=T,
                            constrained=T)

sum(aa_grad_eval$losses)
sum(aa_grad_eval2$losses)
plotExpression(expMat = aa_grad_eval2$Predicted_expression_mat, .dendrogram = "none", .Rowv = F, .Colv = F, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),filename = "grad_test.png")
plotExpression(expMat = aa_grad_eval$Predicted_expression_mat, .dendrogram = "none", .Rowv = F, .Colv = F, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),filename = "grad_test_1.png")

aa_time_optim <- proc.time()
aa_optim <- optim(par = aa_start_par, fn = obj_func_chosen_enh_gradient, my_data = ER_52_opt_simAnn_2_input, chosen_enh_ind = Sim_Ann_weighted_148_restart_enhancers[1, ],intercept_pergene=T, constrained=F, return_grad=F)
proc.time()
aa_optim_eval <- objFuncEval(my_optimized_par = aa_optim$par,
                             my_data = ER_52_opt_simAnn_2_input,
                             input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                             round_to_int = T,
                             chosen_enh = Sim_Ann_weighted_148_restart_enhancers[1, ],
                             log_sum_exp=F,
                             intercept_pergene=T,
                             constrained=F)
aa_time_sbplx <- proc.time()
aa_sbplx <- sbplx(x0 = aa_start_par,
                  fn = obj_func_chosen_enh_gradient,
                  my_data = ER_52_opt_simAnn_2_input, 
                  chosen_enh_ind = Sim_Ann_weighted_148_restart_enhancers[1, ],
                  intercept_pergene=T, 
                  constrained=F,
                  return_grad=F)
proc.time()
aa_sbplx_eval <- objFuncEval(my_optimized_par = aa_sbplx$par,
                             my_data = ER_52_opt_simAnn_2_input,
                             input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                             round_to_int = T,
                             chosen_enh = integer(0),
                             log_sum_exp=F,
                             intercept_pergene=T,
                             constrained=F)
########################################################################################################################
########################################################################################################################
SimAnnAcceptProb <- function(cur_sol, prev_sol, temp){
  return(exp((prev_sol - cur_sol)/temp))
}
########################################################################################################################
########################################################################################################################
SimAnnObjOptim <- function(.my_data,
                           initial_chosen_enh=integer(0),
                           initial_par=numeric(0),
                           .intercept_pergene,
                           starting_temp = 1,
                           par_range = c(-10,10),
                           min_temp,
                           iter_per_temp,
                           max_gene_per_iter=integer(0),
                           .sig_over_loss=F,
                           ALPHA,
                           my_algorithm = "NLOPT_LN_COBYLA"){
  # .my_data : is the output of "createOptimDataset"
  # initial_chosen_enh : is an initial set of chosen enhancers, one per gene
  # initial_par is a numeric vector of size (number_of_TFs + 1) which specifies the initial parameters
  # par_range is a matrix containing the range from which the starting parameters will be taken from, each row corresponds to one parameter and first column is lower bound and second column is upper bound
  # .intercept_pergene : logical if Ture there is one intercept per gene
  # min_temp is the minimum temperature of simulated annealing
  # iter_per_temp is the number of iterations in each temprature
  # max_gene_per_iter : is the maximum number of genes to go through during each iteration. defaults to all genes if not specified.
  # .sig_over_loss : WILL BE PASSED TO obj_func_chosen_enh. Logical, if True it uses a sigmoid on the errors of each gene-cond pair
  # ALPHA is the coefficient by which the temprature decreases
  library(nloptr)
  nu_genes <- length(.my_data$Affinity_scores)
  if(length(max_gene_per_iter) == 0){
    max_gene_per_iter <- nu_genes
  }
  nu_enh_per_gene <- unlist(lapply(.my_data$Affinity_scores, nrow))
  nu_TFs <- ncol(.my_data$Affinity_scores[[1]])
  # initialize initial_chosen_enh if it's already not initialized
  if(length(initial_chosen_enh) == 0){
    initial_chosen_enh <- integer(nu_genes)
    for(i in 1:nu_genes){
      initial_chosen_enh[i] <- sample(c(1:nu_enh_per_gene[i]), size = 1)
    }
  }
  # initialize parameters if it's already not initialized
  nu_pars <- integer(0)
  if(.intercept_pergene){
    nu_pars <- nu_genes + nu_TFs
  }else{
    nu_pars <- 1 + nu_TFs
  }
  
  if(length(par_range) == 2){
    par_range <- cbind(rep(par_range[1], nu_pars),
                       rep(par_range[2], nu_pars))
  }else{
    stopifnot(is.matrix(par_range),
              nrow(par_range) == nu_pars)
  }
  my_lb <- par_range[, 1]
  my_ub <- par_range[, 2]
  if(length(initial_par) == 0){
    my_par <- numeric(nu_pars)
    for(i in 1:nrow(par_range)){
      my_par[i] <- runif(1, min = my_lb[i], max = my_ub[i])
    }
  }else{
    my_par <- initial_par
  }
  #initialize the list to store the results.
  cnt <- 1
  Result_list <- list()
  Result_list$Enhancer_index <- list()
  Result_list$parameters <- list()
  Result_list$obj_value <- list()
  # Do the first optimization
  local_opts <- list(algorithm = "NLOPT_LD_MMA",
                     xtol_rel=1e-04)
  opts <- list( "algorithm" = my_algorithm
                ,"xtol_rel"  = 1.0e-04,
                "maxeval"   = 100
                ,"local_opts" = local_opts
                # ,      
                # check_derivatives = T, 
                # check_derivatives_tol = 1e-04,
                # check_derivatives_print='all'
  )
  First_optim <- nloptr( x0=my_par,
                         eval_f=obj_func_chosen_enh,
                         my_data = .my_data,
                         chosen_enh_ind = initial_chosen_enh,
                         intercept_pergene=.intercept_pergene, 
                         constrained=F,
                         return_grad = F, 
                         sig_over_loss = .sig_over_loss,
                         # lb=my_lb,
                         # ub=my_ub,
                         #eval_g_ineq=obj_func_chosen_enh_ineq,
                         #     eval_g_eq=eval_g_eq,
                         opts=opts)
  # First_optim <- optim(par = my_par,
  #                      fn = obj_func_chosen_enh,
  #                      my_data = .my_data,
  #                      chosen_enh_ind = initial_chosen_enh)
  Result_list$Enhancer_index[[cnt]] <- initial_chosen_enh
  Result_list$parameters[[cnt]] <- First_optim$solution
  Result_list$obj_value[[cnt]] <- First_optim$objective
  cnt <- cnt + 1
  cnt_all <- 1
  Temp <- starting_temp

  while(Temp > min_temp){
    for(cur_iter in 1:iter_per_temp){
      #print("cur_iter:")
      #print(cur_iter)
      shuffled_gene_index <- sample(c(1:nu_genes), size = nu_genes, replace = F)
      for(cur_gene in 1:max_gene_per_iter){
        #print("cur_gene")
        #print(cur_gene)
        cur_chosen_enh_temp <- Result_list$Enhancer_index[[cnt-1]]
        if(nu_enh_per_gene[shuffled_gene_index[cur_gene]] > 1){
          cur_new_enh <- sample(x = setdiff(c(1:nu_enh_per_gene[shuffled_gene_index[cur_gene]]),
                                            Result_list$Enhancer_index[[cnt-1]][shuffled_gene_index[cur_gene]]),
                                size = 1)
          cur_chosen_enh_temp[shuffled_gene_index[cur_gene]] <- cur_new_enh
          
          cur_optim <- nloptr( x0=Result_list$parameters[[cnt-1]],
                               eval_f=obj_func_chosen_enh,
                               my_data = .my_data,
                               chosen_enh_ind = cur_chosen_enh_temp,
                               intercept_pergene=.intercept_pergene, 
                               constrained=F,
                               return_grad = F, 
                               sig_over_loss = .sig_over_loss,
                               # lb=my_lb,
                               # ub=my_ub,
                               #eval_g_ineq=obj_func_chosen_enh_ineq,
                               #     eval_g_eq=eval_g_eq,
                               opts=opts)
          
          # cur_optim <- optim(par = Result_list$parameters[[cnt-1]],
          #                    fn = obj_func_chosen_enh,
          #                    my_data = .my_data,
          #                    chosen_enh_ind = cur_chosen_enh_temp)
          cnt_all <- cnt_all + 1
          print(paste("number of conducted optims:", cnt_all))
          accprob <- SimAnnAcceptProb(cur_sol = cur_optim$objective,
                                      prev_sol = Result_list$obj_value[[cnt-1]],
                                      temp = Temp)
          if(accprob > runif(1)){
            Result_list$Enhancer_index[[cnt]] <- cur_chosen_enh_temp
            Result_list$parameters[[cnt]] <- cur_optim$solution
            Result_list$obj_value[[cnt]] <- cur_optim$objective
            print(paste("number of accepted solutions:", cnt))
            cnt <- cnt + 1
          }
        }

      }
    }
    Temp <- Temp * ALPHA
    print("temp is")
    print(Temp)
  }
  return(Result_list)
}
########################################################################################################################
########################################################################################################################
#example
aa <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[c(1:4), ],
                         TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01,
                         ChoppedEnhancerScores = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score[c(1:4)]) 

aa_sim_res <- SimAnnObjOptim(.my_data = aa,
                             initial_chosen_enh=integer(0),
                             initial_par=numeric(0),
                             min_temp = 0.1,
                             iter_per_temp = 1,
                             ALPHA = 0.9) 

########################################################################################################################
########################################################################################################################
SimAnnStepJobConstructor <- function(script_name,
                                     nu_jobs,
                                     file_name,
                                     start_temp=1,
                                     min_temp,
                                     alPha,
                                     alpha_Step=1){
  # script_name : is the name of R script to be run
  # nu_jobs : is the number of jobs to be run in parallel
  # file_name : is the prefix name of the files to be written
  # start_temp : is the temperature to start with
  # min_temp : is the minimum temperature of the last job
  # alpha : is the coefficient by which to decrease temperature
  # alpha_Step : is the number of times the job reaches to the step to multiply the temperature by alpha
  
  # this function creates a set of files each containing job list, split
  # into many starting and end temperatures as a means for checkpointing. 
  # each job has three arguments, a number which will be used in random number generation,
  # second number is the starting temp of that job and 
  # third argument is minimum temp of that job
  
  nu_job_type <- ceiling(log((min_temp/start_temp), base = alPha))
  nu_job_type <- ceiling(nu_job_type/alpha_Step)
  for(jb_tp in 1:nu_job_type){
    cur_st_temp <- start_temp * ((alPha)^((jb_tp-1) * alpha_Step))
    cur_end_temp <- start_temp * ((alPha)^(jb_tp * alpha_Step))
    cur_end_temp <- (cur_end_temp + (cur_end_temp/alPha))/2
    for(jb_rep in 1:nu_jobs){
      cat(c("Rscript --vanilla", script_name, jb_rep, cur_st_temp, cur_end_temp, "\n"),
          file = paste(file_name, jb_tp, sep = "_"),
          sep = " ",
          append = T)
    }
  }
}

########################################################################################################################
########################################################################################################################
# example
SimAnnStepJobConstructor(script_name="ER_Optimizer_SimAnn_hal.R",
                         nu_jobs=10,
                         file_name="ER_Optimizer_SimAnn_hal_job",
                         start_temp=1,
                         min_temp=0.01,
                         alPha=0.8,
                         alpha_Step=1)
########################################################################################################################
########################################################################################################################
DAGmanConstructor <- function(submit_prefix,
                              start_ind=1,
                              end_ind,
                              my_order = c(1:(end_ind-start_ind + 1)),
                              filename="my_DagMan"){
  # this function creates inputs for DAGMAN in order to be fed to
  #   "condor_submit_dag" for sequential processing of jobs
  # submit_prefix (a character vector of length one) is the prefix of the submit job. for example
  #   if the submit file names is "job_1.submit" then the prefix would be "job"
  # start_ind (an integer) is the index of the first submit job e.g. if
  #   we start with job_1.submit then it would be 1
  # end_ind is (an integer) the index of the last submit job e.g. if
  #   the last submit file is job_10.submit then it would be 10
  # my_order is an integer vector of length equal to number of submit jobs, containing
  #   the order by which the jobs have to be executed. 1 being the first job. many jobs
  #   can have the same priority if needed
  # filename is the name of the file to write to
  # it's important to note that the submit job file names
  #  should be formatted like this: "submit_prefix"_"index"."submit" e.g. job_1.submit
  
  nu_jobs <- end_ind - start_ind + 1
  # job_ind holds the index of all jobs
  job_ind <- c(start_ind:end_ind)
  # writing the first part introducing the jobs and naming them
  for( cur_job_ind in 1:nu_jobs){
    cat(c("JOB",
          paste("job", job_ind[cur_job_ind], sep = "_"),
          paste0(submit_prefix, "_", job_ind[cur_job_ind],
                 ".submit",
                 "\n")), sep = " ", file = filename ,append = T)
  }
  # writing the second part for ordering
  uniqe_prio <- unique(my_order)
  for(line_nu in 1:(length(uniqe_prio) - 1)){
    cur_parents <- job_ind[which(my_order %in% uniqe_prio[line_nu])]
    cur_childs <- job_ind[which(my_order %in% uniqe_prio[line_nu+1])]
    cat(c("PARENT", paste("job", cur_parents, sep = "_"), " "), file = filename, sep = " ", append = T)
    cat(c("CHILD", paste("job", cur_childs, sep = "_"), "\n"), file = filename, sep = " ", append = T)
  }
}
########################################################################################################################
########################################################################################################################
#example
DAGmanConstructor(submit_prefix = "job",
                  start_ind=1,
                  end_ind=10,
                  my_order = c(rep(1,3), rep(2,3), rep(3,3), 4),
                  filename="my_DagMan")
########################################################################################################################
########################################################################################################################
GreedyObjOptim <- function(.my_data,
                           initial_chosen_enh=integer(0),
                           initial_par=numeric(0),
                           par_range = c(-10,10),
                           gene_number){
  # This function performs a greedy search over all enhancers of each gene to find the one that best matches
  # I will write it to take gene number as input in order to be able to break jobs in multiple parts for feeding into hal
  
  # .my_data : is the output of "createOptimDataset"
  # initial_chosen_enh : is an initial set of chosen enhancers, one per gene
  # initial_par is a numeric vector of size (number_of_TFs + 1) which specifies the initial parameters
  # par_range is the range from which the starting parameters will be taken from
  # gene_number is the number of the index of the gene within the input data
  nu_genes <- length(.my_data$Affinity_scores)
  #print(nu_genes)
  nu_enh_per_gene <- unlist(lapply(.my_data$Affinity_scores, nrow))
  #print(nu_enh_per_gene)
  nu_TFs <- ncol(.my_data$Affinity_scores[[1]])
  # initialize initial_chosen_enh if it's already not initialized
  if(length(initial_chosen_enh) == 0){
    initial_chosen_enh <- integer(nu_genes)
    for(i in 1:nu_genes){
      #print(nu_enh_per_gene[i])
      initial_chosen_enh[i] <- sample(x = c(1:nu_enh_per_gene[i]), size = 1, replace = F)
    }
  }
  # initialize parameters if it's already not initialized
  if(length(initial_par) == 0){
    my_par <- runif((nu_TFs + 1), min = par_range[1], max = par_range[2])
  }else{
    my_par <- initial_par
  }
  # create a vector to store results
  obj_res_holder <- numeric(length = nu_enh_per_gene[gene_number])
  obj_param_holder <- matrix(nrow = nu_enh_per_gene[gene_number], ncol = length(my_par))
  # Do the optimization for each enhancer of the specified gene
  chosen_enh_set <- initial_chosen_enh
  for(cur_enh in 1:nu_enh_per_gene[gene_number]){
    chosen_enh_set[gene_number] <- cur_enh
    Enh_obj_optim <- optim(par = my_par,
                           fn = obj_func_chosen_enh,
                           my_data = .my_data,
                           chosen_enh_ind = chosen_enh_set)
    obj_res_holder[cur_enh] <- Enh_obj_optim$value
    obj_param_holder[cur_enh, ] <- Enh_obj_optim$par
  }
  Best_enh <- which.min(obj_res_holder)
  chosen_enh_set[gene_number] <- Best_enh
  return(list(Enhancer_index=chosen_enh_set,
              gene_ind=gene_number,
              parameters=obj_param_holder[Best_enh, ],
              value=obj_res_holder[Best_enh]))
}
########################################################################################################################
########################################################################################################################
# example
aa <- GreedyObjOptim (.my_data = ER_52_opt_simAnn_input_ChipFilt_1500_500,
                      initial_chosen_enh=integer(0),
                      initial_par=numeric(0),
                      par_range = c(-10,10),
                      gene_number = 1)
########################################################################################################################
########################################################################################################################
GreedyJobConstructor <- function(script_name="ER_Optimizer_Greedy_hal.R",
                                 nu_jobs,
                                 file_name="ER_Optimizer_Greedy_hal_job",
                                 nu_genes,
                                 repeat_rounds,
                                 fixed_genes=integer(0)){
  # script_name : is the name of R script to be run
  # nu_jobs : is the number of jobs to be run in parallel
  # file_name : is the prefix name of the files to be written
  # nu_genes : is the number of genes in the dataset
  # repeat_rounds : is the number of times that you want the search to be done over all genes.
  # fixed_genes : is the index of the genes to be in the same position in all random starts, these are usually long jobs. fixed to be done all in one step instead of slowing down all steps
  # This function writes job files to implement the greedy search function above. each job takes
  # two arguments, one the random start number (which is from 1 to nu_jobs) and a 
  # gene number which is from 1 to nu_genes. jobs will be written in way that can be fed to DAGman
  gene_reading_index <- matrix(nrow = nu_jobs, ncol = (nu_genes * repeat_rounds))
  for(cur_job in 1:nu_jobs){
    for(cur_rep in 1:repeat_rounds){
      aa <- c(1:nu_genes) + (cur_rep - 1)*nu_genes
      gene_reading_index[cur_job, aa] <- sample(x = c(1:nu_genes), size = nu_genes, replace = F)
    }
  }
  
  # fix genes in certain positions for all jobs 
  if(length(fixed_genes) > 0){
    for(i in 1:length(fixed_genes)){
      chosenPos <- sample(x = c(1:nu_genes), size = repeat_rounds)
      for(cur_r in 1:length(chosenPos)){
        chosenPos[cur_r] <- chosenPos[cur_r] + nu_genes * (cur_r-1)
      }
      for(ii in 1:nrow(gene_reading_index)){
        curr_pos <- which(gene_reading_index[ii, ] == fixed_genes[i])
        chosen_pos_content <- gene_reading_index[ii, chosenPos]
        gene_reading_index[ii, chosenPos] <- fixed_genes[i]
        gene_reading_index[ii, curr_pos] <- chosen_pos_content
      }
    }
  }
  
  for(job_tp in 1:ncol(gene_reading_index)){
    for(jb_rep in 1:nu_jobs){
      cat(c("Rscript --vanilla", script_name, jb_rep, gene_reading_index[jb_rep, job_tp], "\n"),
          file = paste(file_name, job_tp, sep = "_"),
          sep = " ",
          append = T)
    }
  }
}
########################################################################################################################
########################################################################################################################
# example
aa <- getwd()
setwd("hal_job_greedy_148_5rep/")
GreedyJobConstructor(script_name="ER_Optimizer_Greedy_hal.R",
                     nu_jobs = 10,
                     file_name="ER_Optimizer_Greedy_hal_job",
                     nu_genes=3,
                     repeat_rounds=2)
setwd(aa)

GreedyJobConstructor(script_name="ER_Optimizer_Greedy_Conc_modif.R",
                     nu_jobs = 12,
                     file_name="ER_Optimizer_Greedy_hal_job",
                     nu_genes=12,
                     repeat_rounds=2,
                     fixed_genes = c(1,2,5,6))
########################################################################################################################
########################################################################################################################
confusion_constructor <- function(predictionMat_List, realMat){
  # predictionMat_List : is a list of prediction matrices
  # realMat : is the real expression matrix
  # This function creates a confusion matrix per prediction in prediction list, and returns them in a list
  
  # The list that'll contain confusion matrices
  confusion_list <- list()
  element_table_real <- table(as.numeric(realMat))
  #print("real table")
  #print(element_table_real)
  nu_elements <- length(element_table_real)
  # create a template for all confusion matrices
  template_mat <- matrix(nrow = nu_elements + 2, ncol = nu_elements + 2)
  rownames(template_mat) <-c(paste("Predicted", names(element_table_real)), "Actual_sum", "sensitivity") 
  colnames(template_mat) <- c(paste("Actual", names(element_table_real)), "Predicted_sum", "precision")
  real_element_position <- list()
  for(i in 1:nu_elements){
    real_element_position[[i]] <- which(as.numeric(realMat) %in% as.numeric(names(element_table_real)[i]))
  }
  names(real_element_position) <- names(element_table_real)
  
  for(cur_pred in 1:length(predictionMat_List)){
    confusion_list[[cur_pred]] <- template_mat
    cur_element_prisiton <- list()
    for(i in 1:nu_elements){
      cur_element_prisiton[[i]] <- which(as.numeric(predictionMat_List[[cur_pred]]) %in% as.integer(names(element_table_real)[i]))
    }
    names(cur_element_prisiton) <- names(element_table_real)
    for(cur_elem_p in 1:nu_elements){
      for(cur_elem_r in 1:nu_elements)
      confusion_list[[cur_pred]][cur_elem_p, cur_elem_r] <- length(intersect(cur_element_prisiton[[cur_elem_p]],
                                                                             real_element_position[[cur_elem_r]]))
    }
    #add sums and sensitivity and precision
    # actual sum
    confusion_list[[cur_pred]][nu_elements+1, 1:nu_elements] <- colSums(confusion_list[[cur_pred]][1:nu_elements, 1:nu_elements])
    # predicted sum 
    confusion_list[[cur_pred]][1:nu_elements, nu_elements+1] <- rowSums(confusion_list[[cur_pred]][1:nu_elements, 1:nu_elements])
    #sensitivity and precision
    for(i in 1:nu_elements){
      #sensitivity
      confusion_list[[cur_pred]][nu_elements+2, i] <- confusion_list[[cur_pred]][i, i]/confusion_list[[cur_pred]][nu_elements + 1, i]
      #precision
      confusion_list[[cur_pred]][i, nu_elements+2] <- confusion_list[[cur_pred]][i, i]/confusion_list[[cur_pred]][i, nu_elements + 1]
    }
  }
  return(confusion_list)
}

########################################################################################################################
########################################################################################################################
#example
aa <- confusion_constructor(predictionMat_List=Sim_Ann_148_restart_ExpMat[1:2], realMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
########################################################################################################################
########################################################################################################################
CreateRandomPrediction <- function(real_exp_mat, num){
  # real_exp_mat is the real expression matrix
  # num in the number of random predictions to produce
  # function returns the accuracy and confusion matrix for each of the random predictions
  
  Random_pred_list <- list()
  Random_pred_list_mat <- list()
  Random_confusion_list <- list()
  Random_precision <- numeric(num)
  real_exp_vec <- as.numeric(real_exp_mat)
  non_na_pos <- which(!is.na(real_exp_vec))
  element_table <- table(real_exp_vec)
  template <- numeric(length = length(real_exp_vec))
  template[is.na(real_exp_vec)] <- NA
  for (nu_pred in 1:num){
    my_samples <- list()
    my_non_na_pos <- non_na_pos
    Random_pred_list[[nu_pred]] <- template
    for(elem in 1:length(element_table)){
      my_samples[[elem]] <- sample(x = my_non_na_pos, size = element_table[elem], replace = F)
      my_non_na_pos <- setdiff(my_non_na_pos, my_samples[[elem]])
      Random_pred_list[[nu_pred]][my_samples[[elem]]] <- as.numeric(names(element_table)[elem])
    }
    assertthat::are_equal(element_table, table(Random_pred_list[[nu_pred]]))
    Random_precision[nu_pred] <- sum(Random_pred_list[[nu_pred]] == real_exp_vec, na.rm = T) / sum(!is.na(real_exp_vec))
    Random_pred_list_mat[[nu_pred]] <- matrix(Random_pred_list[[nu_pred]], byrow = F, ncol = ncol(real_exp_mat))
  }
  # calculate the confusion matrices
  Random_confusion_list <- confusion_constructor(predictionMat_List = Random_pred_list_mat,
                                                 realMat = real_exp_mat)
  return(list(Precision = Random_precision, Confusion_mat_list = Random_confusion_list))
}
########################################################################################################################
########################################################################################################################
# example
aa <- CreateRandomPrediction(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, num = 10)
aa$Precision
aa$Confusion_mat_list[[1]]
########################################################################################################################
########################################################################################################################
tcol <- function(color, trans = 255) {
  
  if (length(color) != length(trans) & 
      !any(c(length(color), length(trans)) == 1)) 
    stop('Vector lengths not correct')
  if (length(color) == 1 & length(trans) > 1) 
    color <- rep(color, length(trans))
  if (length(trans) == 1 & length(color) > 1) 
    trans <- rep(trans, length(color))
  
  res <- paste0('#', apply(apply(rbind(col2rgb(color)), 2, function(x) 
    format(as.hexmode(x), 2)), 2, paste, collapse = ''))
  res <- unlist(unname(Map(paste0, res, as.character(as.hexmode(trans)))))
  res[is.na(color)] <- NA
  return(res)
}
########################################################################################################################
########################################################################################################################
ggcols <- function (n, l = 65, c = 100) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = l, c = c)[1:n]
}
########################################################################################################################
########################################################################################################################
Enhancer_vote_barplot <- function(Enhancer_matrix, nu_enh_per_gene, .ylim = numeric(0)){
  # Enhancer_matrix is a matrix where each row represents a model, each
  #  column represents a gene. each entry is the enhancer index (within 
  #  that gene's enhancers) chosen for that gene in that model
  # nu_enh_per_gene is an integer vector containing the total number of enhancers per gene.
  # Plots a bar per gene, divides the bar into fragments representing enhancers of the gene
  input_enh <- matrix(0L, nrow = max(Enhancer_matrix), ncol = ncol(Enhancer_matrix))
  for(i in 1:ncol(Enhancer_matrix)){
    cur_table <- table(Enhancer_matrix[, i])
    input_enh[as.integer(names(cur_table)), i] <- cur_table
  }
  input_enh <- apply(X = input_enh, MARGIN = 2, FUN = sort, decreasing = T)
  if(length(.ylim )==0){
    .ylim = c(-nrow(input_enh), nrow(Enhancer_matrix))
  }
  barplot(input_enh,
          col = sample(ggcols(150), 150, replace = F),
          ylim = .ylim,
          space = 0.5,
          names.arg = colnames(Enhancer_matrix),
          las = 2)
  barplot(-nu_enh_per_gene, space = 0.5,las = 2,  add = T)
  #return()
  return(input_enh)
}
########################################################################################################################
########################################################################################################################
#example
aaa <- unlist(lapply(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score, nrow))
Enhancer_vote_barplot(Sim_Ann_148_restart_enhancers, aaa)
########################################################################################################################
########################################################################################################################
Enhancer_Affinity_heatmap <- function(Enhancer_matrix,
                                      affinity_score_list,
                                      .col_vec=col_vector,
                                      exportplot=T,
                                      filename="affinity_heatmap.png",
                                      .dendrogram="column",
                                      .ColSideColors=numeric(0),
                                      .Rowv=F,
                                      .Colv=T){
  # Enhancer_matrix is a matrix where each row represents a model, each
  #  column represents a gene. each entry is the enhancer index (within 
  #  that gene's enhancers) chosen for that gene in that model
  # affinity_score_list : a list where each entry corresponds to a gene and contains a matrix of 
  #  enhancer affinities in which each row is an enhancer, each column is a TF
  # filename : plot file name
  enh_ind_list <- list()
  affinity_mat <- matrix(nrow = 0, ncol = ncol(affinity_score_list[[1]]))
  side_col <- character(0)
  for(i in 1:ncol(Enhancer_matrix)){
    enh_ind_list[[i]] <- unique(Enhancer_matrix[, i])
    affinity_mat <- rbind(affinity_mat, affinity_score_list[[i]][enh_ind_list[[i]], ])
    side_col <- c(side_col, rep(.col_vec[i%%length(.col_vec)], length(enh_ind_list[[i]])))
  }
  
  # create the heatmap
  #colors for heatmap
  my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
  #export or not?
  if(exportplot){
    png(filename,    # create PNG for the heat map        
        width = 8*300,        # 5 x 300 pixels
        height = 8*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  #plot the heatmap using gplot heatmap.2
  
  # col_breaks = c(seq(0.0, 0.33, length=100),     # for red
  #                seq(0.34, 0.66, length=100),   # for yellow
  #                seq(0.67, 1, length=100))         # for green
  if(length(.ColSideColors)==0){
    .ColSideColors <- rep("white", ncol(affinity_mat))
  }
  heatmap.2(scale(affinity_mat),
            #cellnote = format(round(SSEmat, 2), nsmall = 2),  # same data set for cell labels
            main = "enh_affinity", # heat map title
            notecol="black",      # change font color of cell labels to black
            # density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(9,7),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            dendrogram=.dendrogram,     # only draw a .dendrogram dendrogram
            RowSideColors = side_col,# grouping row-variables(genes) into different categories
            ColSideColors = .ColSideColors,# grouping column-variables(models) into different categories
            Rowv = .Rowv,
            Colv = .Colv,
            scale = "none"
  ) 
  if(exportplot){
    dev.off() 
  }
  return(affinity_mat)
}
########################################################################################################################
########################################################################################################################
#example
aa_GR_loss <- rowSums(do.call(rbind, lapply(Optim_Greedy_conc_modif_Ensemble_Results_eval[[4]], "[[", 3)))
aa_GR_loss_sort_ind <- sort(aa_GR_loss, decreasing = F, index.return=T)$ix
aa_GR_enh <- do.call(rbind, lapply(Optim_Greedy_conc_modif_Ensemble_Results_eval[[4]], "[[", 2))
aaa <- unlist(lapply(ER_52_opt_simAnn_2_input$Affinity_scores, nrow))
par(mfrow=c(1,1))
Enhancer_vote_barplot(aa_GR_enh[aa_GR_loss_sort_ind[1:10], ], aaa/15, .ylim = c(-10, 10))

aa <- Enhancer_Affinity_heatmap(Enhancer_matrix = aa_GR_enh[aa_GR_loss_sort_ind[1:10], ],
                          affinity_score_list = ER_52_opt_simAnn_2_input$Affinity_scores,
                          .col_vec=col_vector,
                          exportplot=T,
                          filename="affinity_heatmap_row_col_norm.png",
                          .dendrogram="column",
                          .ColSideColors=numeric(0),
                          .Rowv=F,
                          .Colv=T)
nrow(aa)
########################################################################################################################
########################################################################################################################

PerformanceHeattMap_General <- function(real_exp_mat,
                                        prediction_mat_list,
                                        my_perf_mat=numeric(0),
                                        gene_index = c(1:nrow(real_exp_mat)),
                                        model_index= c(1:length(prediction_mat_list)),
                                        .Colv = T, .Rowv = F, .dendrogram = "col",exportplot = T,
                                        filename = "PerformanceHeatMap.png",
                                        .RowSideColors = character(0), .ColSideColors = character(0), nu_col=100){
  # real_exp_mat is the real expression matrix
  # prediction_mat_list : is the list of predicted expression matrices
  # gene_index : index of genes to be plotted
  # model_index : index of models to be plotted
  # .Colv : logical, if the columns be reordered by clustering
  # .Rowv : logical, if Rows be reorder by clustering
  # .dendrogram : can be "col", "row", "both", "none"
  # exportplot : logical, if to export plot 
  # filename : name of the exported plot
  # .RowSideColors is a charachter vector with length equal to number of genes indicating the color of the sidebar for each gene
  # .ColSideColors is a charachter vector with length equal to number of models indicating the color of the sidebar for each model
  # nu_col : number of colors making gradient of each color
  eq_func <- function(x, y=real_exp_mat){
    all_acc <- numeric(length = nrow(x))
    for(i in 1:nrow(x)){
      all_acc[i] <- sum(x[i, ] == y[i, ], na.rm = T)/sum(!is.na(x[i, ]))
    }
    return(all_acc)
  }
  
  # if(length(.RowSideColors) == 0){
  #   .RowSideColors = rep("white", length(gene_index))
  # }
  
  if(length(.ColSideColors) == 0){
    .ColSideColors = rep("white", length(model_index))
  }
  nu_models_total <- length(prediction_mat_list)
  nu_genes_total <- nrow(real_exp_mat)
  
  
  if(length(my_perf_mat)==0){
    my_acc_list <- lapply(X = prediction_mat_list, FUN = eq_func)
    Performance_mat_total <- do.call(cbind, my_acc_list)
    rownames(Performance_mat_total) <- rownames(real_exp_mat)
    my_Perf_mat <- Performance_mat_total[gene_index, model_index]
  }else{
    my_Perf_mat <- my_perf_mat[gene_index, model_index]
    rownames(Performance_mat_total) <- rownames(real_exp_mat)
  }
  ## get the number of non-NAs per gene and assign the siderowcolor based on that
  num_non_na <- rowSums(!is.na(real_exp_mat))
  if(length(.RowSideColors) ==0){
    .RowSideColors <- character(length(num_non_na))
    .RowSideColors[num_non_na <= 20] <- "red"
    .RowSideColors[num_non_na > 20  &  num_non_na <= 30] <- "yellow"
    .RowSideColors[num_non_na > 30] <- "green"
  }
  
  #colors for heatmap
  my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 3*nu_col - 1)
  #export or not?
  if(exportplot){
    png(filename,    # create PNG for the heat map        
        width = 8*300,        # 5 x 300 pixels
        height = 8*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  #plot the heatmap using gplot heatmap.2
  
  col_breaks = c(seq(0.0, 0.33, length=nu_col),     # for red
                 seq(0.34, 0.66, length=nu_col),   # for yellow
                 seq(0.67, 1, length=nu_col))         # for green
  heatmap.2(my_Perf_mat,
            #cellnote = format(round(SSEmat, 2), nsmall = 2),  # same data set for cell labels
            main = "performance", # heat map title
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
            scale = "none"
  ) 
  if(exportplot){
    dev.off() 
  }
  return(my_Perf_mat)
}
########################################################################################################################
########################################################################################################################
#example
aaa <- unlist(lapply(strsplit(colnames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), split ="_"), "[[", 1))
aaa_col_color <- character(length(Sim_Ann_148_restart_ExpMat))
for(i in 1:length(unique(aaa))){
  aaa_col_color[aaa %in% unique(aaa)[i]] <- col_vector[i]
}
aa <- PerformanceHeattMap_General(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                  prediction_mat_list = Sim_Ann_148_restart_ExpMat,
                                  .Colv = T, .Rowv = T, .dendrogram = "both",exportplot = T,
                                  filename = "PerformanceHeatMap.png",
                                  .RowSideColors = character(0), .ColSideColors = aaa_col_color)


########################################################################################################################
########################################################################################################################
SilicoConcenOptim <- function(model_param=numeric(0),
                              init_enh=integer(0),
                              Modif_vec,
                              ..my_data,
                              .gene_number){
  # model_param is a numeric vector containing the starting parameters of the next optimization
  # init_enh is an integer vector, specifying the index of each gene's enhancer.
  # Modif_vec is an integer vector with the same length as the number of TFs. where
  # if 0: no modification, if 1: change to 1 in treatment conditions and to 0 in control conditions, if -1: change to 0 in treatment conditions and 1 in control conditions.
  # ..my_data is the output of createOptimDataset() function. for the current dataset. the input TF conc is not modified.
  # .gene_number is the index of the gene whose enhancers are subject to change.
  # outputs the same as GreedyObjOptim() function. the only change is to use the modified TF concentration.
  
  nu_TFs <- ncol(..my_data$Affinity_scores[[1]])
  modified_data <- ..my_data
  # setting controls to zero if 1
  modified_data$gene_conditon_mat[, which(Modif_vec == 1)] <- 0
  # setting treatments to one if 1
  modified_data$gene_conditon_mat[, (which(Modif_vec == 1) + nu_TFs)] <- 1
  # setting controls to one if -1
  modified_data$gene_conditon_mat[, which(Modif_vec == -1)] <- 1
  # setting treatments to zero if -1
  modified_data$gene_conditon_mat[, (which(Modif_vec == -1) + nu_TFs)] <- 0
  
  # running greedy optimization with the modified data
  mod_Res <- GreedyObjOptim(.my_data = modified_data,
                            initial_chosen_enh = init_enh,
                            initial_par = model_param,
                            gene_number = .gene_number)
  return(mod_Res)
}

########################################################################################################################
########################################################################################################################
#example
aa <- SilicoConcenOptim(..my_data = ER_52_opt_simAnn_input,
                        init_enh=my_initial_chosen_enh,
                        model_param=my_initial_par,
                        .gene_number=1,
                        Modif_vec=Optim_Greedy_conc_modif_vec[1, ])
########################################################################################################################
########################################################################################################################
insertOnes <- function(template_vec, insert_pos, insert_vec=numeric(0), cols = c("grey", "red")){
  # this is a helper for Enhancer_Score_plot. to insert elements of "insert_vec" at positions specified by "insert_pos", into template_vec
  # insert_pos is an integer vector containg the positions of template_vec after whoch the corresponding "insert_vec" element should be inserted
  # cols is a vector containing two colors, first one for the normal entries, second one for spacing entries
  
  if(length(insert_vec) == 0){
    insert_vec <- integer(length(insert_pos))
    insert_vec[] <- 1
  }
  my_contigs <- list()
  for(i in 1:length(insert_pos)){
    if(i == 1){
      my_contigs[[i]] <- template_vec[1:insert_pos[i]]
    }else{
      my_contigs[[i]] <- template_vec[(insert_pos[i-1] + 1):insert_pos[i]]
    }
  }
  my_contigs[[i+1]] <- template_vec[(insert_pos[i] + 1):length(template_vec)]
  
  new_vec <- numeric(0)
  my_col_vec <- character(0)
  for(i in 1:length(my_contigs)){
    new_vec <- c(new_vec, my_contigs[[i]])
    my_col_vec <- c(my_col_vec, rep(cols[1], length(my_contigs[[i]])))
    if(i != length(my_contigs)){
      new_vec <- c(new_vec, insert_vec[i])
      my_col_vec <- c(my_col_vec, cols[2])
    }
  }
  results <- rbind(new_vec, my_col_vec)
  rownames(results) <- c("new_vec", "colors")
  return(results)
}
########################################################################################################################
########################################################################################################################
# example
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_accuracy_perGene_perEnh, enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges)
aaa <- insertOnes(template_vec = aa$score$`10057`, insert_pos = which(aa$overlap$`10057` == 0), insert_vec=numeric(0), cols = c("grey", "red"))
aaa <- insertOnes(template_vec = aa$score$`10125`, insert_pos = which(aa$overlap$`10125` == 0), insert_vec=numeric(0), cols = c("grey", "red"))

########################################################################################################################
########################################################################################################################

Enhancer_Score_plot <- function(accuracy_per_model_per_gene_per_enh,
                                enhancer_Granges,
                                model_index=c(1:length(accuracy_per_model_per_gene_per_enh)),
                                export_plot = T,
                                filename = "enhancer_score_barplot.png",
                                .cols = c("grey", "red"),
                                loss = F,
                                real_exp_mat=numeric(0),
                                .ylim = integer(0),
                                draw_plot=T){
  # loss_per_model_per_gene_per_enh is a list containing one entry per model. that entry is
  #  a list containing one entry for each gene. that entry is a numeric vector indicating the
  #  accuracy for that gene given each enhancer.
  # enhancer_Granges is a list containing one entry per gene. that entry is a GRanges object of the enhancers of that gene.
  # model_index is the index of the models that performance is averaged over.
  # .cols : is a character vector of colors, first entry used for bars of normal entries. second entry used for spacer bars
  # loss : logical. if True then the first entry loss_per_model_per_gene_per_enh will be treated as squared error, instead of accuracy
  # real_exp_mat is the real expression matrix
  # .ylim is the y limits of the plots. if not provided it default to c(0, 1) or c(0, 1.5*(expected_loss(gene)))
  
  # it calculates the average score of each enhancer over all models in model_index. creates
  #  a barplot per gene indicating that score. overlapping enhancers have bars close to each
  #  other while non-overlapping enhancers have some empty space in between their bars.
  if(loss){ 
    # get the square root if it's loss --> removed this step because I calculate the standard deviation of sum of squared errors and can't just take sqrt for that as well. in order to keep the expected value and its sd in the same scale i don't take the sqrt at this step.
    # for(i in 1:length(accuracy_per_model_per_gene_per_enh)){
    #   accuracy_per_model_per_gene_per_enh[[i]] <-  lapply(accuracy_per_model_per_gene_per_enh[[i]], sqrt)
    # }
    expected_loss <- numeric(nrow(real_exp_mat))
    names(expected_loss) <- names(enhancer_Granges)
    std_loss <- numeric(nrow(real_exp_mat))
    names(std_loss) <- names(enhancer_Granges)
  }

  
  nu_genes <- length(enhancer_Granges)
  nu_ones <- numeric(nu_genes)
  nu_minus_ones <- numeric(nu_genes)
  nu_zero <- numeric(nu_genes)
  nu_check <- numeric(nu_genes)
  nu_nonNA <- numeric(nu_genes)
  avg_score <- list()
  my_overlap_indicator <- list()
  accuracy_per_model_per_gene_per_enh_working <- accuracy_per_model_per_gene_per_enh[model_index]

  for(cur_gene in 1:length(accuracy_per_model_per_gene_per_enh_working[[1]])){
    #print("cur_gene")
    #print(cur_gene)
    # calculating the average
    cur_gene_list <- lapply(accuracy_per_model_per_gene_per_enh_working, "[[", cur_gene)
    cur_gene_mat <- do.call(rbind, cur_gene_list)
    cur_gene_avg_vec <- apply(cur_gene_mat, MARGIN = 2, FUN = mean)
    avg_score[[cur_gene]] <- cur_gene_avg_vec
    
    # overlap calc
    cur_gene_overlap <- findOverlaps(query=enhancer_Granges[[cur_gene]], subject=enhancer_Granges[[cur_gene]])
    
    # a binary vector indicating if the current enhancer overlaps with the next enhancer in the GRanges object
    cur_overlap_binary <- integer(length = length(enhancer_Granges[[cur_gene]]) - 1)
    if(length(cur_overlap_binary) > 0){
      for(cur_enh in 1:(length(enhancer_Granges[[cur_gene]]) - 1)){
        #print("cur_enh")
        #print(cur_enh)
        cur_enh_overlap <- cur_gene_overlap@to[cur_gene_overlap@from == cur_enh]
        if(length(cur_enh_overlap) > 0){
          #print("yes")
          cur_overlap_binary[cur_enh] <- as.integer((cur_enh + 1) %in% cur_enh_overlap)
        }
      }# end loop over enhancers
    }
    my_overlap_indicator[[cur_gene]] <- cur_overlap_binary
  }#end loop over genes
  names(avg_score) <- names(enhancer_Granges)
  # print(my_overlap_indicator)
  names(my_overlap_indicator) <- names(enhancer_Granges)
  
  # create the list for bar plot by basically
  # inserting a 1 number after the scores for enhancers that do not overlap with
  # the next enhancer and creating a color vector as well
  avg_score_plot_list <- list()
  color_vector_per_gene <- list()
  for(cur_gene in 1:length(avg_score)){
    
    # compute the expected score by shuffling for each gene
    nu_ones[cur_gene] <- sum(real_exp_mat[cur_gene, ] == 1, na.rm = T)
    nu_minus_ones[cur_gene] <- sum(real_exp_mat[cur_gene, ] == -1, na.rm = T)
    nu_nonNA[cur_gene] <- sum(!is.na(real_exp_mat[cur_gene, ]))
    nu_zero[cur_gene] <- nu_nonNA[cur_gene] - nu_ones[cur_gene] - nu_minus_ones[cur_gene]
    nu_check[cur_gene] <- sum(real_exp_mat[cur_gene, ] == 0, na.rm = T)
    if(nu_check[cur_gene] != nu_zero[cur_gene]){
      print("something is wrong in counting the number of 1, -1, 0")
      return(0)
    }
    if(loss){
      expected_loss[cur_gene] <- 2 * ((nu_ones[cur_gene] + nu_minus_ones[cur_gene]) - ((nu_ones[cur_gene] - nu_minus_ones[cur_gene])^2)/nu_nonNA[cur_gene])
      std_loss[cur_gene] <- nu_nonNA[cur_gene] * (-(expected_loss[cur_gene]/nu_nonNA[cur_gene])^2 + (expected_loss[cur_gene]/nu_nonNA[cur_gene]) + (24 * nu_ones[cur_gene] * nu_minus_ones[cur_gene] / (nu_nonNA[cur_gene])^2))
      # print("(expected_loss[cur_gene])^2")
      # print((expected_loss[cur_gene])^2)
      # print("(expected_loss[cur_gene])")
      # print((expected_loss[cur_gene]))
      # print("(24 * nu_ones[cur_gene] * nu_minus_ones[cur_gene] / (nu_nonNA[cur_gene])^2)")
      # print((24 * nu_ones[cur_gene] * nu_minus_ones[cur_gene] / (nu_nonNA[cur_gene])^2))
      # print("std_loss[cur_gene]")
      # print(std_loss[cur_gene])
      std_loss[cur_gene] <- sqrt(std_loss[cur_gene])
    }
    
    if(length(avg_score[[cur_gene]]) > 1){ #if gene has more than one enhancer
      all_non_overlapp <- which(my_overlap_indicator[[cur_gene]] == 0)
      if(length(all_non_overlapp) > 0){#if there is any non overlapping enhancer
        if(loss){
          my_insert_vec <- integer(length(all_non_overlapp))
          my_insert_vec[] <- expected_loss[cur_gene] 
        }else{
          my_insert_vec <- numeric(0)
        }

        cur_inserted <- insertOnes(template_vec = avg_score[[cur_gene]], insert_pos =all_non_overlapp, insert_vec=my_insert_vec, cols = .cols)
        avg_score_plot_list[[cur_gene]] <- as.numeric(cur_inserted[1, ])
        color_vector_per_gene[[cur_gene]] <- cur_inserted[2, ]
      }else{#if there is no non overlapping enhancer
        avg_score_plot_list[[cur_gene]] <- avg_score[[cur_gene]]
        color_vector_per_gene[[cur_gene]] <- rep(.cols[1], length(avg_score[[cur_gene]]))
      }

    }else{ #if gene has only one enhancer
      avg_score_plot_list[[cur_gene]] <- avg_score[[cur_gene]]
      color_vector_per_gene[[cur_gene]] <- .cols[1]
    }
    

  }#end of loop over genes
  
  if(draw_plot){
    if(export_plot){
      png(filename,    # create PNG for the heat map        
          width = 8*300,        # 5 x 300 pixels
          height = 1.5*nu_genes*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8)        # smaller font size
    }
    par(mfrow = c(nu_genes, 1), mar= c(3, 5, 2, 1))
    
    for(i in 1:nu_genes){
      #print(avg_score_plot_list[[i]])
      if(length(.ylim) == 0){
        if(loss){
          .ylim <- c(0, 1.5*expected_loss[i])
        }else{
          .ylim <- c(0, 1)
        }
      }
      
      barplot(height = avg_score_plot_list[[i]],
              col = color_vector_per_gene[[i]],
              ylab = names(avg_score)[i],
              #yaxt = "n",
              cex.axis = 2,
              ylim = .ylim
      )
      if(loss){
        abline(h=expected_loss[i], col=3, lty=5, lwd = 1.2)
        abline(h= expected_loss[i] + std_loss[i] , col=3, lty=3, lwd = 0.8)
        abline(h= expected_loss[i] - std_loss[i] , col=3, lty=3, lwd = 0.8)
        abline(h= expected_loss[i] + 2 * std_loss[i] , col=3, lty=3, lwd = 0.8)
        abline(h= expected_loss[i] - 2 * std_loss[i] , col=3, lty=3, lwd = 0.8)
        # abline(h=0.9*(expected_loss[i]), col=3, lty=8, lwd = 0.5)
        # abline(h=0.8*(expected_loss[i]), col=3, lty=8, lwd = 0.5)
        # abline(h=0.7*(expected_loss[i]), col=3, lty=8, lwd = 0.5)
        # abline(h=0.6*(expected_loss[i]), col=3, lty=8, lwd = 0.5)
        # abline(h=0.5*(expected_loss[i]), col=3, lty=8, lwd = 0.5)
        # abline(h=0.4*(expected_loss[i]), col=3, lty=8, lwd = 0.5)
        # abline(h=0.3*(expected_loss[i]), col=3, lty=8, lwd = 0.5)
        # abline(h=0.2*(expected_loss[i]), col=3, lty=8, lwd = 0.5)
        # abline(h=0.1*(expected_loss[i]), col=3, lty=8, lwd = 0.5)
        
        mtext(paste("0:", nu_zero[i], ", 1:", nu_ones[i], ", -1:", nu_minus_ones[i]), side = 3, line = -0.1)
      }else{
        abline(h=seq(0,1, 0.2), col=3, lty=8, lwd = 0.5)
        mtext(paste("0:", nu_zero[i], ", 1:", nu_ones[i], ", -1:", nu_minus_ones[i]), side = 3, line = -0.1)
        
      }
    }
    if(export_plot){
      dev.off()
    }
    
  }

  return(list(score=avg_score,
              overlap=my_overlap_indicator,
              plotted_bar_list=avg_score_plot_list,
              expected_score_shuffling=expected_loss,
              expected_std_shuffling =std_loss))
}

########################################################################################################################
########################################################################################################################
#example
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_accuracy_perGene_perEnh,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          real_exp_mat =my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52 )

aa_sort <- sort(unlist(Sim_Ann_148_restart_round_precision), index.return=T, decreasing = T)$ix
Sim_Ann_148_restart_loss_perGene_perEnh_rounded <- lapply(Sim_Ann_148_restart_all_models_eval, "[[", 7)
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_loss_perGene_perEnh_rounded,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = aa_sort[1],
                          filename = "enhancer_score_barplot_top1model_class_loss_calib_round.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
########################################################################################################################
########################################################################################################################
TF_KD_evaluator <- function(my_optimized_par_mat,
                            .my_data,
                            .input_gene_expMat,
                            .round_to_int = T,
                            TF_index,
                            WT_evaluation=numeric(0)
                            ){
  # my_optimized_par_mat : is a matrix of parameter values, each row
  #  corresponds to a model, each column to a parameter. this matrix has to be named
  # .my_data : input data to the optimization
  # .input_gene_expMat : original gene expression matrix
  # round_to_int : round predictions to closest integer. provide a
  #  performance estimate by rounding the predictions and counting the correct calls.
  # TF_index : index of the TFs to be knockDown, index refers to TF input's order in .my_data
  # WT_evaluation is a list with length equal to number of models (nrow(my_optimized_par_mat)) containing
  #  the results of objFuncEval on the input for WT conditions, be careful that the
  #  .round_to_int parameter should've been the same as the one currently used
  # it outputs a list, containing one entry per TF KD experiment. within that entry:
  #  average performance of each model before and after TF KD,
  #  average performance of each gene before and after TF KD
  
  nu_genes <- nrow(.input_gene_expMat)
  nu_cond <- ncol(.input_gene_expMat)
  nu_enh <- sum(unlist(lapply(.my_data$Affinity_scores, nrow)))
  nu_model <- nrow(my_optimized_par_mat)
  nu_TFs <- ncol(.my_data$Affinity_scores[[1]])
  
  if(length(WT_evaluation) == 0){ # create WT evaluations if not already available
    WT_evaluation <- list()
    print("Wild type evaluation")
    for(i in 1:nu_model){
      print(i)
      WT_evaluation[[i]] <- objFuncEval(my_optimized_par = my_optimized_par_mat[i,],
                                        my_data = .my_data,
                                        input_gene_expMat = .input_gene_expMat ,
                                        round_to_int = .round_to_int)
    }
    names(WT_evaluation) <- rownames(my_optimized_par_mat)
  }
  
  KD_evaluation <- list()
  print("KD  evaluation")
  for(cur_tf in 1:length(TF_index)){
    KD_evaluation[[cur_tf]] <- list()
    modified_data <- .my_data
    modified_data$gene_conditon_mat[, TF_index[cur_tf]] <- 0 # set the control condition expression to zero
    modified_data$gene_conditon_mat[, (TF_index[cur_tf] + nu_TFs)] <- 0 # set the tretment condition expression to zero
    for(cur_mod in 1:nu_model){
      print(paste("tf number", cur_tf, "model number", cur_mod))
      KD_evaluation[[cur_tf]][[cur_mod]] <- objFuncEval(my_optimized_par = my_optimized_par_mat[cur_mod, ],
                                                        my_data = modified_data,
                                                        input_gene_expMat = .input_gene_expMat ,
                                                        round_to_int = .round_to_int)
    }
    names(KD_evaluation[[cur_tf]]) <- rownames(my_optimized_par_mat)
  }
  names(KD_evaluation) <- colnames(.my_data$Affinity_scores[[1]])[TF_index]
  
  # for each TF compute the average change in performance
  average_loss_change_per_model <- list()
  average_Accuracy_change_per_model <- list()
  average_loss_change_per_gene <- list()
  average_Accuracy_change_per_gene <- list()
  average_loss_change_per_enhancer <- list()
  average_Accuracy_change_per_enhancer <- list()
  for(cur_tf in 1:length(TF_index)){
    average_loss_change_per_model[[cur_tf]] <- numeric(nu_model)
    
    if(!.round_to_int){
      print("there will be problems with accuracy measurements because of not rounding to int,
            either set .round_to_int=T or change the objFuncEval function such that it outputs
            the rounded accuracy even when round_to_int is set to False")
    }
    average_Accuracy_change_per_model[[cur_tf]] <- numeric(nu_model)
    
    average_loss_change_per_gene[[cur_tf]] <- numeric(nu_genes)
    average_Accuracy_change_per_gene[[cur_tf]] <- numeric(nu_genes)
    average_loss_change_per_enhancer[[cur_tf]] <- numeric(nu_enh)
    average_Accuracy_change_per_enhancer[[cur_tf]] <- numeric(nu_enh)
    for(cur_mod in 1:nu_model){
      average_loss_change_per_model[[cur_tf]][cur_mod] <- (sum(WT_evaluation[[cur_mod]]$losses) - sum(KD_evaluation[[cur_tf]][[cur_mod]]$losses))/length(WT_evaluation[[cur_mod]]$losses)
      #print(WT_evaluation[[cur_mod]]$Rounded_Error)
      #print(KD_evaluation[[cur_mod]]$Rounded_Error)
      average_Accuracy_change_per_model[[cur_tf]][cur_mod] <- WT_evaluation[[cur_mod]]$Rounded_Error - KD_evaluation[[cur_tf]][[cur_mod]]$Rounded_Error
      average_loss_change_per_gene[[cur_tf]] <- average_loss_change_per_gene[[cur_tf]] + (WT_evaluation[[cur_mod]]$losses - KD_evaluation[[cur_tf]][[cur_mod]]$losses)
      average_Accuracy_change_per_gene[[cur_tf]] <- average_Accuracy_change_per_gene[[cur_tf]] + (unlist(lapply(WT_evaluation[[cur_mod]]$accuracy_perGene_perEnh, max)) - unlist(lapply(KD_evaluation[[cur_tf]][[cur_mod]]$accuracy_perGene_perEnh, max)))
      average_loss_change_per_enhancer[[cur_tf]] <- average_loss_change_per_enhancer[[cur_tf]] + (unlist(WT_evaluation[[cur_mod]]$loss_perGene_perEnh) - unlist(KD_evaluation[[cur_tf]][[cur_mod]]$loss_perGene_perEnh))
      average_Accuracy_change_per_enhancer[[cur_tf]] <- average_Accuracy_change_per_enhancer[[cur_tf]] + (unlist(WT_evaluation[[cur_mod]]$accuracy_perGene_perEnh) - unlist(KD_evaluation[[cur_tf]][[cur_mod]]$accuracy_perGene_perEnh))
    }
    average_loss_change_per_gene[[cur_tf]] <- average_loss_change_per_gene[[cur_tf]]/nu_model
    average_Accuracy_change_per_gene[[cur_tf]] <- average_Accuracy_change_per_gene[[cur_tf]]/nu_model
    average_loss_change_per_enhancer[[cur_tf]] <- average_loss_change_per_enhancer[[cur_tf]]/nu_model
    average_Accuracy_change_per_enhancer[[cur_tf]] <- average_Accuracy_change_per_enhancer[[cur_tf]]/nu_model
    #name them
    names(average_loss_change_per_model[[cur_tf]]) <- rownames(my_optimized_par_mat)
    names(average_Accuracy_change_per_model[[cur_tf]]) <- rownames(my_optimized_par_mat)
    names(average_loss_change_per_gene[[cur_tf]]) <- rownames(.input_gene_expMat)
    names(average_Accuracy_change_per_gene[[cur_tf]]) <- rownames(.input_gene_expMat)
  }
  names(average_loss_change_per_model) <- names(KD_evaluation)
  names(average_Accuracy_change_per_model) <- names(KD_evaluation)
  names(average_loss_change_per_gene) <- names(KD_evaluation)
  names(average_Accuracy_change_per_gene) <- names(KD_evaluation)
  names(average_loss_change_per_enhancer) <- names(KD_evaluation)
  names(average_Accuracy_change_per_enhancer) <- names(KD_evaluation)
  
  
  return(list(WT_eval=WT_evaluation,
              KD_eval=KD_evaluation,
              loss_change_per_model=average_loss_change_per_model,
              loss_change_per_gene=average_loss_change_per_gene,
              loss_change_per_enhancer=average_loss_change_per_enhancer,
              accu_change_per_model=average_Accuracy_change_per_model,
              accu_change_per_gene=average_Accuracy_change_per_gene,
              accu_change_per_enhancer=average_Accuracy_change_per_enhancer))
}

aa <- objFuncEval(my_optimized_par = Sim_Ann_148_restart_parameters[1,],
            my_data = ER_52_opt_simAnn_input,
            input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52 ,
            round_to_int = F, log_sum_exp = F)
########################################################################################################################
########################################################################################################################
#example
aa <- TF_KD_evaluator(my_optimized_par_mat = Sim_Ann_148_restart_parameters[1:3,],
                      .my_data=ER_52_opt_simAnn_input,
                      .input_gene_expMat=my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                      .round_to_int = T,
                      TF_index=c(1, 5),
                      WT_evaluation=numeric(0))

########################################################################################################################
########################################################################################################################
objFuncEval_multiEnh <- function(){
  # I want this objective function to create the output expression matrix
  #  using the parameters learned in training. If enhancers are not 
  
}

########################################################################################################################
########################################################################################################################


########################################################################################################################
########################################################################################################################
# next function starts here

########################################################################################################################
########################################                       #########################################################
########################################        Analysis       #########################################################
########################################                       #########################################################
########################################################################################################################

# run the optimization for minimum enhancer
(aaresult2 <- optim(par = runif(n = ((ncol(aa$gene_conditon_mat) - 1)/2 + 1), min = -100, max = 100), fn = obj_func_linear_minEnh, my_data = aa))

aa <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                         TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun,
                         ChoppedEnhancerScores = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score) 

aaLinearOptMinEnh <- list()
aaEvaluatedOptMinEnh <- list()
# next for line commented out to prevent accidental runs
# for(i in 1:100){
#   print(i)
#   aaLinearOptMinEnh[[i]] <- optim(par = runif(n = ((ncol(aa$gene_conditon_mat) - 1)/2 + 1), min = -20, max = 20), fn = obj_func_linear_minEnh, my_data = aa)
# }
for(i in 1:27){
  print(i)
  aaEvaluatedOptMinEnh[[i]] <- objFuncEval(my_optimized_par = aaLinearOptMinEnh[[i]]$par, my_data = aa, input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
}
# look at optimized obj function value
plot(unlist(lapply(aaLinearOptMinEnh, "[[", 2)), type = "l")

# look at error per enhancer, create enhancer error variance matrix rows genes, columns optimization runs
aaEnhancerErrorVar <- matrix(nrow = 52, ncol = 27)
for(i in 1:27){
  for (j in 1:52){
    aaEnhancerErrorVar[j, i] <- var(aaEvaluatedOptMinEnh[[i]]$loss_perGene_perEnh[[j]])
  }
}
aaEvaluatedOptMinEnh[[13]]$loss_perGene_perEnh[[6]]

for(i in 1:27){
  plotExpression(expMat = aaEvaluatedOptMinEnh[[i]]$Predicted_expression_mat, .dendrogram = "none", .Rowv = F, filename = paste("Prediction_",i,".png", sep = ""))
  
}

aaaaaa_optim_dataset <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                           TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun,
                                           ChoppedEnhancerScores = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score) 

aa_linear_hal_Res <- list()
aa_linear_hal_Res_eval <- list()
aa <- list.files("linear_hal_res/")
for(i in 1:length(aa)){
  load(paste0("linear_hal_res/", aa[i]))
  aaa <- which(unlist(lapply(OptimizedExp, length)) == 5)
  aa_linear_hal_Res[[i]] <- OptimizedExp[[aaa]]
  aa_linear_hal_Res_eval[[i]] <- objFuncEval(my_optimized_par = aa_linear_hal_Res[[i]]$par, my_data = aaaaaa_optim_dataset, input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
  remove(OptimizedExp)
}

load("linear_hal_res/Optim_res_1.RData")
aa <- OptimizedExp
remove(OptimizedExp)
load("linear_hal_res/Optim_res_2.RData")
aa2 <- OptimizedExp

for(i in 1:length(aa_linear_hal_Res_eval)){
  plotExpression(expMat = aa_linear_hal_Res_eval[[i]]$Predicted_expression_mat, .dendrogram = "none", .Rowv = F, filename = paste("Prediction_hal_",i,".png", sep = ""))
  
}
####################################################################################################v####################
##################################################################################################################################
##############################  SIMULATED ANNEALING ANALYSIS ############################################################
########################################################################################################################
##############################################################################################################
#creating input for hal runs
ER_52_opt_simAnn_input <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                             TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01,
                                             ChoppedEnhancerScores = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score) 
# write the hal job file all at once
for(i in 1:148){
  cat(c("Rscript --vanilla ER_Optimizer_SimAnn_hal.R", i, "\n"),
      file = "ER_Optimizer_SimAnn_hal_job",
      sep = " ",
      append = T)
}
####################
# write hal jobs for different temperatures
SimAnnStepJobConstructor(script_name="ER_Optimizer_SimAnn_hal.R",
                         nu_jobs=148,
                         file_name="ER_Optimizer_SimAnn_hal_job",
                         start_temp=1,
                         min_temp=1e-5,
                         alPha=0.95,
                         alpha_Step=1)
##############################################################
step_counter <- 1
nu_jobs <- 148
##############################################################

cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_script", sep="")
for(i in 1:225){
  cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
        paste("ER_Optimizer_SimAnn_hal_job", i, sep = "_"),
        nu_jobs,
        paste("tmp_ER_Optimizer_SimAnn_hal_job", i, sep = "_"),
        ">",
        paste("tmp_ER_Optimizer_SimAnn_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_script",
      sep=" ",
      append = T)
}
cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_exec", sep="")
for(i in 1:225){
  cat(c("chmod +x",
        paste("tmp_ER_Optimizer_SimAnn_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_exec",
      sep=" ",
      append = T)
}
##############
DAGmanConstructor(submit_prefix = "tmp_ER_Optimizer_SimAnn_hal_job",
                  start_ind=7,
                  end_ind=225,
                  #order = c(rep(1,3), rep(2,3), rep(3,3), 4),
                  filename="my_DagMan_7_225")
#####################################################################################################################
#####################################################################################################################
# weighted sim ann run:
#creating input for hal runs
# normalize using max LR
aaa_max_consLLR <- numeric(length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward))
for(i in 1:length(aaa_max_consLLR)){
  aaa_max_consLLR[i] <- ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward[[i]][1, 2]
}
names(aaa_max_consLLR) <- names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward)

aaa_max_consLR <- unlist(lapply(aaa_max_consLLR, exp))

aa_norm_score <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score
for(i in 1:length(aa_norm_score)){
  for(j in 1:nrow(aa_norm_score[[i]])){
    aa_norm_score[[i]][j, ] <- aa_norm_score[[i]][j, ] * aaa_max_consLR
  }
}

ER_52_opt_simAnn_2_input <- ER_52_opt_simAnn_input
ER_52_opt_simAnn_2_input$Affinity_scores <- aa_norm_score

####################
# write hal jobs for different temperatures
aa_wd <- getwd()
setwd("Hal_simAnn_weighted_jobs/")
SimAnnStepJobConstructor(script_name="ER_Optimizer_SimAnn_weighted_hal.R",
                         nu_jobs=148,
                         file_name="ER_Optimizer_SimAnn_weighted_hal_job",
                         start_temp=1,
                         min_temp=1e-6,
                         alPha=0.95,
                         alpha_Step=1)

##############################################################
step_counter <- 1
nu_jobs <- 148
##############################################################

cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_script", sep="")
for(i in 1:270){
  cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
        paste("ER_Optimizer_SimAnn_weighted_hal_job", i, sep = "_"),
        nu_jobs,
        paste("tmp_ER_Optimizer_SimAnn_weighted_hal_job", i, sep = "_"),
        ">",
        paste("ER_Optimizer_SimAnn_weighted_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_script",
      sep=" ",
      append = T)
}
cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_exec", sep="")
for(i in 1:270){
  cat(c("chmod +x",
        paste("ER_Optimizer_SimAnn_weighted_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_exec",
      sep=" ",
      append = T)
}
##############
DAGmanConstructor(submit_prefix = "ER_Optimizer_SimAnn_weighted_hal_job",
                  start_ind=1,
                  end_ind=270,
                  #order = c(rep(1,3), rep(2,3), rep(3,3), 4),
                  filename="my_DagMan_1_270")
DAGmanConstructor(submit_prefix = "ER_Optimizer_SimAnn_weighted_hal_job",
                  start_ind=190,
                  end_ind=270,
                  #order = c(rep(1,3), rep(2,3), rep(3,3), 4),
                  filename="my_DagMan_190_270")
setwd(aa_wd)
#####################################################################################################################
#####################################################################################################################
# Enhancer seq, scores and GRANGES filtered for the ones that overlap or interact with one ER chip peak
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_ChIPfilt

ER_52_opt_simAnn_input_ChipFilt <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                                      TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01,
                                                      ChoppedEnhancerScores = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_ChIPfilt$Chopped_Score) 
# Enhancer seq, scores and GRANGES filtered for the ones that overlap or interact with one ER chip peak and pieces of length 1500 and step size of 500
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt

ER_52_opt_simAnn_input_ChipFilt_1500_500 <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                                               TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01,
                                                               ChoppedEnhancerScores = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt$Chopped_Score) 

sum(unlist(lapply(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene_ChIPfilt$Chopped_Seq, length)))
########################################################################################################################
########################################################################################################################
# Reading Simulated annealing results
Sim_Ann_148_restart_results <- list()
for(i in 1:148){
  load(paste0("SimAnn_Results/Optim_sim_Ann_", i, ".RData"))
  Sim_Ann_148_restart_results[[i]] <- SimAnnRes[[length(SimAnnRes)]]
  remove(SimAnnRes)
}

aa_enh <- lapply(Sim_Ann_148_restart_results, "[[", 1)
names(aa_enh) <- c(1:length(aa_enh))
for(i in 1:length(aa_enh)){
  aa_enh[[i]] <- do.call(rbind, aa_enh[[i]])
  rownames(aa_enh[[i]]) <- c(1:nrow(aa_enh[[i]]))
}

aa_par <- lapply(Sim_Ann_148_restart_results, "[[", 2)
names(aa_par) <- c(1:length(aa_par))
for(i in 1:length(aa_par)){
  aa_par[[i]] <- do.call(rbind, aa_par[[i]])
  rownames(aa_par[[i]]) <- c(1:nrow(aa_par[[i]]))
}
aa_obj <- lapply(Sim_Ann_148_restart_results, "[[", 3)
names(aa_obj) <- c(1:length(aa_obj))
for(i in 1:length(aa_obj)){
  names(aa_obj[[i]]) <- c(1:length(aa_obj[[i]]))
}

aa_obj <- unlist(aa_obj)
aa_enh <- do.call(rbind, aa_enh)
rownames(aa_enh) <- names(aa_obj)
aa_par <- do.call(rbind, aa_par)
rownames(aa_par) <- names(aa_obj)
########################### Grab the last step of each random start in sim Ann
aa_names_1 <- unlist(lapply(strsplit(rownames(aa_par), split = "\\."), "[[", 1))
aa_names_2 <- unlist(lapply(strsplit(rownames(aa_par), split = "\\."), "[[", 2))
aa_names_1_unique <- unique(aa_names_1)
aa_max_ind <- integer(length(aa_names_1_unique))

for(i in 1:length(aa_names_1_unique)){
  aa_ind <- which(aa_names_1 %in% aa_names_1_unique[i])
  aa_indm <- which.max(aa_names_2[aa_ind])
  aa_max_ind[i] <- aa_ind[aa_indm]
}

aa_obj_last <- aa_obj[aa_max_ind]
aa_par_last <- aa_par[aa_max_ind, ]
boxplot.matrix(aa_par_last )
aa_enh_last <- aa_enh[aa_max_ind, ]
boxplot.matrix(aa_enh_last)
aa_s_ind <- sort(aa_obj_last, decreasing = F, index.return = T)$ix
aa_obj_last_sort <- aa_obj_last[aa_s_ind]
aa_enh_last_sort <- aa_enh_last[aa_s_ind, ]
aa_par_last_sort <- aa_par_last[aa_s_ind, ]
#######
setwd("SimAnn_Results_plot/All_148_last_run")
aa_pred_expMat <- list()
for(i in 1:nrow(Sim_Ann_148_restart_parameters)){
  aa_pred_expMat[[i]] <- objFuncEval(my_optimized_par = Sim_Ann_148_restart_parameters[i, ],
                                     my_data = ER_52_opt_simAnn_input,
                                     input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                     round_to_int = T)
  
  plotExpression(expMat = aa_pred_expMat[[i]]$Predicted_expression_mat, .dendrogram = "none", .Rowv = F, .Colv = F, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),filename = paste0("sim_ann_",i, "_" ,rownames(Sim_Ann_148_restart_parameters)[i], ".png"))
}

plotExpression(expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, .dendrogram = "both", .Rowv = T, .Colv = T,filename = "real_exp_2.png", .distfun = my_dist_fuc, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5))

names(aa_pred_expMat) <- rownames(aa_par_last_sort)

Sim_Ann_148_restart_ExpMat <- lapply(Sim_Ann_148_restart_all_models_eval, "[[", 1)
Sim_Ann_148_restart_loss_accum <- lapply(Sim_Ann_148_restart_all_models_eval, "[[", 3)
Sim_Ann_148_restart_loss_perGene_perEnh <- lapply(Sim_Ann_148_restart_all_models_eval, "[[", 4)
Sim_Ann_148_restart_round_accuracy <- lapply(Sim_Ann_148_restart_all_models_eval, "[[", 5)
Sim_Ann_148_restart_parameters <- aa_par_last_sort
Sim_Ann_148_restart_enhancers <- aa_enh_last_sort
Sim_Ann_148_restart_accuracy_perGene_perEnh <- lapply(Sim_Ann_148_restart_all_models_eval, "[[", 6)
Sim_Ann_148_restart_loss_perGene_perEnh_rounded <- lapply(Sim_Ann_148_restart_all_models_eval, "[[", 7)
colnames(Sim_Ann_148_restart_enhancers) <- rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
Sim_Ann_148_restart_objVal <- aa_obj_last_sort

setwd("../..")
# plot difference from real
setwd("SimAnn_Results_plot/All_148_last_run_comp/")
for(i in 1:length(aa_pred_expMat)){
  aa_comp <- abs(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52 - aa_pred_expMat[[i]]$Predicted_expression_mat)
  plotExpression(expMat = aa_comp, .dendrogram = "none", .Rowv = F, .Colv = F,filename = paste0("sim_ann_com_",i, "_" ,rownames(aa_par_last_sort)[i], ".png"))
  
}
setwd("../..")
aa_rounded_error <- unlist(lapply(aa_pred_expMat, "[[", 5))
names(aa_rounded_error) <- names(aa_pred_expMat)
plot(aa_rounded_error, ylim = c(0, 1), ylab = "performance", xlab = "models")
abline(h = 0.55, col = 2)
abline(h = seq(0.30, 0.65, 0.05), col = 2)
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### Create confusion matrix for each model
# number of -1, 0, 1 s that have been predicted correctly versus the number present in the dataset
# columns are actual classes
# rows are predicted classes
Sim_Ann_148_restart_confusion <- confusion_constructor(predictionMat_List=Sim_Ann_148_restart_ExpMat,
                                                       realMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
names(Sim_Ann_148_restart_confusion) <- names(Sim_Ann_148_restart_ExpMat)
Sim_Ann_148_restart_confusion$`83.51`
which.max(Sim_Ann_148_restart_round_precision)
Sim_Ann_148_restart_confusion[[122]]

aa <- lapply(Sim_Ann_148_restart_confusion, as.numeric)
aaa <- do.call(rbind, aa)

aa_S <-sort(aaa[, 1], decreasing = T, index.return=T)$ix
Sim_Ann_148_restart_confusion[[aa_S[3]]]

##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### 
# Create random predictions to compare with sim ann predictions
Sim_Ann_148_restart_random_pred <- CreateRandomPrediction(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, num = 148)
par(mfrow = c(1,1), mar = c(6,6,6,6))
plot(c(unlist(Sim_Ann_148_restart_round_precision),
       Sim_Ann_148_restart_random_pred$Precision),
     ylim = c(0, 1),
     col = c(rep(3, 148), rep(2, 148)),
     ylab = "Accuracy",
     xlab = "model index",
     cex = 0.8,
     pch = 16)
legend(200, 1, legend=c("Sim_Ann", "Random"),
       col=c("green", "red"), pch=16, cex=1.3,
       bty = "n", x.intersp = 0.1,
       y.intersp = 0.2, pt.cex = 1.3, box.col = 1)
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# Draw heatmap of the ensemble
aaa <- unlist(lapply(strsplit(colnames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), split ="_"), "[[", 1))
aaa_col_color <- character(length(Sim_Ann_148_restart_ExpMat))
for(i in 1:length(unique(aaa))){
  aaa_col_color[aaa %in% unique(aaa)[i]] <- col_vector[i]
}
aa <- PerformanceHeattMap_General(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                  prediction_mat_list = Sim_Ann_148_restart_ExpMat,
                                  .Colv = T, .Rowv = T, .dendrogram = "both",exportplot = T,
                                  filename = "PerformanceHeatMap.png",
                                  .RowSideColors = character(0), .ColSideColors = aaa_col_color)
PerformanceHeattMap_General(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                            prediction_mat_list = Sim_Ann_148_restart_ExpMat,
                            .Colv = F, .Rowv = T, .dendrogram = "row",exportplot = T,
                            filename = "PerformanceHeatMap.png",
                            .RowSideColors = character(0), .ColSideColors = aaa_col_color)
boxplot.matrix(aa, las = 2)
aaa <- colSums(aa)
which.max(aaa)
aaa <- apply(aa, 2, median)
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# draw heatmap of the parameters

aa_abs_param <- Sim_Ann_148_restart_parameters
aa_added <- range(aa_abs_param)[1]
log10(-aa_added)
aa_abs_param <- aa_abs_param - range(aa_abs_param)[1] + 1
aa_abs_param <- log10(aa_abs_param)
range(aa_abs_param)

aa <- kmeans(aa_abs_param, centers = 3)

# MDS of parameters

aaad <- dist(aa_abs_param) # euclidean distances between the rows
aafit <- cmdscale(aaad,eig=TRUE, k=3) # k is the number of dim
aafit # view results

# plot solution 
aax <- aafit$points[,1]
aay <- aafit$points[,2]
plot(aax, aay)
text(aax, aay, labels = row.names(aa_abs_param), cex=.7)
aaa_breaks = c(seq(0, 2.21, length=2),     # for red
               seq(2.22, 5, length=2),           # for yellow
               seq(5.01, 10, length=2)) 



ParameterHeatMap(parameterMat = aa_abs_param[1:100,],
                 exportplot = T,
                 filename = "parameterHeatMap.png",
                 .Rowv = T, .Colv = T,
                 .dendrogram = "both",
                 logTransform = F,
                 #.cellnote = format(round(log10(parameterMat), 2), nsmall = 2),
                 .RowSideColors = character(0),
                 col_breaks = aaa_breaks,
                 nu_colors = 5
                 )
aaa_breaks = c(seq(range(Sim_Ann_148_restart_parameters)[1], 0, length=100),     # for red
               seq(0.1, 1000, length=100),           # for yellow
               seq(1000.1, range(Sim_Ann_148_restart_parameters)[2], length=100))

####### Normalize parameters using the LLR(max) of each motif: the affinity scores for each TF is sum(LR/LR(max)), hence here I divide the parameters by LR(max)

aaa_max_consLLR <- numeric(length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward))
for(i in 1:length(aaa_max_consLLR)){
  aaa_max_consLLR[i] <- ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward[[i]][1, 2]
}
names(aaa_max_consLLR) <- names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward)

aaa_max_consLR <- unlist(lapply(aaa_max_consLLR, exp))
aaa_max_consLR <- aaa_max_consLR/min(aaa_max_consLR)

aa_norm_par <- Sim_Ann_148_restart_parameters
for(i in 1:nrow(aa_norm_par)){
  aa_norm_par[i, 1:19] <- aa_norm_par[i, 1:19]/aaa_max_consLR
}
boxplot.matrix(aa_norm_par, outline = T, las = 2, ylim = c(-25, 2000))


aaa_breaks = c(
  seq(range(aa_norm_par[1:80,])[1], -1, length=5),    
               seq(-0.99, 0, length=3),           
               seq(0.01, 1, length=3),
               seq(1.01, 10, length = 5),
               seq(10.01, 100, length = 5),
               seq(100.01, 500, length = 5)
               #seq(500.01, 2000, length = 20)
               ) 
aa_norm_par2 <- aa_norm_par
aa_norm_par2[aa_norm_par2 > 2000] <- 2000

aa_norm_par2 <- cbind(aa_norm_par2, unlist(Sim_Ann_148_restart_round_accuracy) * 100)
colnames(aa_norm_par2)[ncol(aa_norm_par2)] <- "accuracy_percent"

png(filename = "parameters_and_accuracy.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # smaller font size
heatmap.2(x = aa_norm_par2,
          Rowv = T, Colv = T,
          dendrogram = "both",
          #rowsep = aa_rowsep-1, sepwidth = c(4,4), sepcolor = "green",
          #RowSideColors = aa_rowsidecol,
          trace="none", na.rm = T, symbreaks=F , symkey=F, symm = F, breaks = aaa_breaks,col = colorRampPalette(c("green", "black", "red"))(n = 25))
dev.off()

aa <- lapply(Sim_Ann_148_restart_confusion, as.numeric)
aaa <- do.call(rbind, aa)

aa_S <-sort(aaa[, 1], decreasing = T, index.return=T)$ix
boxplot.matrix(aa_norm_par2[aa_S[1:20], ], las = 2, ylim = c(-50, 200))
aa_S <-sort(aaa[, 5], decreasing = T, index.return=T)$ix
boxplot.matrix(aa_norm_par2[aa_S[1:20], ], las = 2, ylim = c(-50, 200))
aa_S <-sort(aaa[, 9], decreasing = T, index.return=T)$ix
boxplot.matrix(aa_norm_par2[aa_S[1:20], ], las = 2, ylim = c(-50, 200))

# ParameterHeatMap(parameterMat = aa_norm_par2[1:80,],
#                  exportplot = T,
#                  filename = "parameterHeatMap.png",
#                  .Rowv = T, .Colv = T,
#                  .dendrogram = "both",
#                  logTransform = F,
#                  #.cellnote = format(round(log10(parameterMat), 2), nsmall = 2),
#                  .RowSideColors = character(0),
#                  col_breaks = aaa_breaks,
#                  col_break_ch_nu = 6
#                 
# )

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### Find enhancers per gene
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

# going over all models with accuracy higher than 0.55
# within each model, for each gene select the top "max(3, no.of.enhancers)" for the gene
# aggregate results over all models.
# see how many enhancers have been discarded on average
aa_par_last_sort_filtered <-  aa_par_last_sort[aa_rounded_error > 0.6, ]
aa_enh_last_sort_filtered <-  aa_enh_last_sort[aa_rounded_error > 0.6, ]
aa_obj_last_sort_filtered <-  aa_obj_last_sort[aa_rounded_error > 0.6]
aa_pred_expMat_filtered <- aa_pred_expMat[aa_rounded_error > 0.6]
aa_enhancer_list <- list()
for(i in 1:52){
  aa_enhancer_list[[i]] <- integer(0)
}
names(aa_enhancer_list) <- names(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score)

for (ii_model in 1:nrow(aa_par_last_sort_filtered)){
  #print("model")
  #print(ii_model)
  for(ii_gene in 1:length(aa_pred_expMat_filtered[[ii_model]]$loss_perGene_perEnh)){
    #print("gene")
    #print(ii_gene)
    aa_enhancer_list[[ii_gene]] <- union(aa_enhancer_list[[ii_gene]],
                                         sort(aa_pred_expMat_filtered[[ii_model]]$loss_perGene_perEnh[[ii_gene]],
                                              index.return=T,
                                              decreasing = F)$ix[1:(min(1, length(aa_pred_expMat_filtered[[ii_model]]$loss_perGene_perEnh[[ii_gene]])))])
    #print("is.na")
    #print(sum(is.na(aa_enhancer_list[[ii_gene]])))
  }
}

aa_enh_ratio <- unlist(lapply(aa_enhancer_list, length)) / unlist(lapply(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score
                                                                         , nrow))
barplot(unlist(lapply(aa_enhancer_list, length)), las = 2)
barplot(aa_enh_ratio, las=2, ylab = "ratio", xlab="genes")
#####################################################################################################################
##### Find enhancers per gene##### Find enhancers per gene##### Find enhancers per gene
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# Using voting, for models better than 0.55 collect votes on enhancers per gene, see which enhancers have the highest votes for each gene
aa_table_list <- list()
for(i in 1:ncol(Sim_Ann_148_restart_enhancers)){
  aa_table_list[[i]] <- sort(table(Sim_Ann_148_restart_enhancers[, i]), decreasing = T)
}
aaa <- unlist(lapply(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score, nrow))
Enhancer_vote_barplot(Sim_Ann_148_restart_enhancers[1:50, ], aaa/4, .ylim=c(-50, 50))
#####################################################################################################################
######### plotting the average score of each enhancer
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_accuracy_perGene_perEnh,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          #model_index = c(1:75),
                          filename = "enhancer_score_barplot_allmodels.png")
aaa <- unlist(aa$score)
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_accuracy_perGene_perEnh,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = c(1:75),
                          filename = "enhancer_score_barplot_top75models.png")
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_accuracy_perGene_perEnh,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = c(1),
                          filename = "enhancer_score_barplot_top1models.png")
aa_sort <- sort(unlist(Sim_Ann_148_restart_round_precision), index.return=T, decreasing = T)$ix
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_accuracy_perGene_perEnh,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = aa_sort[1],
                          filename = "enhancer_score_barplot_top1model_class.png")
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_accuracy_perGene_perEnh,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = aa_sort[1:5],
                          filename = "enhancer_score_barplot_top5model_class.png")
#####################################################################################################################
#### look at the loss in comparison to expected loss
aa_sort <- sort(unlist(Sim_Ann_148_restart_round_precision), index.return=T, decreasing = T)$ix
Sim_Ann_148_restart_loss_perGene_perEnh_rounded <- lapply(Sim_Ann_148_restart_all_models_eval, "[[", 7)
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_loss_perGene_perEnh_rounded,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = aa_sort[1:5],
                          filename = "enhancer_score_barplot_top5model_class_loss_calib_round.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_loss_perGene_perEnh_rounded,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = aa_sort[1:10],
                          filename = "enhancer_score_barplot_top10model_class_loss_calib_round.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)

aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_loss_perGene_perEnh_rounded,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          #model_index = c(1:10),
                          filename = "enhancer_score_barplot_all_loss_calib_rounded.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_loss_perGene_perEnh_rounded,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = 1,
                          filename = "enhancer_score_barplot_top1_loss_calib_rounded.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_loss_perGene_perEnh_rounded,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = aa_sort[1],
                          filename = "enhancer_score_barplot_top1_class_loss_calib_rounded.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)

#### looking at the raw loss (not rounded)
Sim_Ann_148_restart_loss_perGene_perEnh <- lapply(Sim_Ann_148_restart_all_models_eval, "[[", 4)
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_loss_perGene_perEnh,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = c(1:10),
                          filename = "enhancer_score_barplot_top10_loss_calib.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)

#####################################################################################################################
########### looking at affinity of different enhancers for TFs
# normalize using max LR
aaa_max_consLLR <- numeric(length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward))
for(i in 1:length(aaa_max_consLLR)){
  aaa_max_consLLR[i] <- ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward[[i]][1, 2]
}
names(aaa_max_consLLR) <- names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward)

aaa_max_consLR <- unlist(lapply(aaa_max_consLLR, exp))

aa_norm_score <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score
for(i in 1:length(aa_norm_score)){
  for(j in 1:nrow(aa_norm_score[[i]])){
    aa_norm_score[[i]][j, ] <- aa_norm_score[[i]][j, ] * aaa_max_consLR
  }
}
#par(mfrow= c(5, 2), mar= c(1, 4, 1, 1))
png("enhancer_affinity_boxplot.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = length(aa_norm_score)*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)        # smaller font size
par(mfrow = c(length(aa_norm_score), 1) , mar= c(0.1, 4, 0.1, 16))

for(i in 1:(length(aa_norm_score)-1)){
  boxplot.matrix(log10(aa_norm_score[[i]]), las = 2, cex = 0.5, xaxt = "n")
  mtext(paste("gene:", names(aa_norm_score)[i], "nu_enh:", nrow(aa_norm_score[[i]])), side = 4, line = 1, adj = 0, las = 1)
}
par(mar= c(4, 4, 0.1, 16))
boxplot.matrix(log10(aa_norm_score[[i+1]]), las = 2, cex = 0.5)
mtext(paste("gene:", names(aa_norm_score)[i+1], "nu_enh:", nrow(aa_norm_score[[i+1]])), side = 4, line = 1, adj = 0, las = 1)
dev.off()

#### heat map of enhancer affinities
aa_norm_score_bind <- do.call(what = rbind, aa_norm_score)
quantile(log10(aa_norm_score_bind))
aaa_breaks = c(
  seq(range(log10(aa_norm_score_bind))[1], 0, length=3),    
  seq(0.00001, 1, length=2),           
  seq(1.00001, 2, length=2),
  seq(2.00001, 4, length = 3),
  seq(4.00001, 6, length = 3),
  seq(6.00001, range(log10(aa_norm_score_bind))[2], length = 3)
  #seq(500.01, 2000, length = 20)
) 
aa_sc <- unlist(lapply(aa_norm_score, nrow))
aa_rowsidecol <- character(0)
for(i in 1:length(aa_sc)){
  aa_rowsidecol <- c(aa_rowsidecol, rep(col_vector[i%%length(col_vector)], aa_sc[i]))
}

ParameterHeatMap(parameterMat = log10(aa_norm_score_bind),
                 exportplot = T,
                 filename = "enhancer_affinity_heatmap.png",
                 .Rowv = F, .Colv = F,
                 .dendrogram = "none",
                 logTransform = F,
                 #.cellnote = format(round(log10(parameterMat), 2), nsmall = 2),
                 .RowSideColors = aa_rowsidecol,
                 col_breaks = aaa_breaks,
                 col_break_ch_nu = 6
                 
)

#### Decide on a smaller set of TFs based on simulated KD results
aa_tf_filtered_ind <- c(3, 4, 5, 7, 11, 12, 13, 16, 18, 19)
colnames(aa_norm_score_bind)[aa_tf_filtered_ind]
aa_norm_score_bind_filtered <- aa_norm_score_bind[, aa_tf_filtered_ind]

#### normalizing each column
par(mfrow = c(4, 3), mar = c(3, 2, 2, 1))
for(i in 1:ncol(aa_norm_score_bind_filtered)){
  hist(log10(aa_norm_score_bind_filtered[, i]), breaks = 100, main = "")
  mtext(colnames(aa_norm_score_bind_filtered)[i], side = 3)
}
library(preprocessCore)
aa_norm_score_bind_filtered_quantile <- normalize.quantiles(log10(aa_norm_score_bind_filtered))
rownames(aa_norm_score_bind_filtered_quantile) <- rownames(aa_norm_score_bind_filtered)
colnames(aa_norm_score_bind_filtered_quantile) <- colnames(aa_norm_score_bind_filtered)
par(mfrow = c(4, 3), mar = c(3, 2, 2, 1))
for(i in 1:ncol(aa_norm_score_bind_filtered)){
  hist(aa_norm_score_bind_filtered_quantile[, i], breaks = 100, main = "")
  mtext(colnames(aa_norm_score_bind_filtered_quantile)[i], side = 3)
}
aa_norm_score_bind_filtered_zscore <- scale(log10(aa_norm_score_bind_filtered))
rownames(aa_norm_score_bind_filtered_zscore) <- rownames(aa_norm_score_bind_filtered)
colnames(aa_norm_score_bind_filtered_zscore) <- colnames(aa_norm_score_bind_filtered)
par(mfrow = c(4, 3), mar = c(3, 2, 2, 1))
for(i in 1:ncol(aa_norm_score_bind_filtered)){
  hist(aa_norm_score_bind_filtered_zscore[, i], breaks = 100, main = "")
  mtext(colnames(aa_norm_score_bind_filtered_zscore)[i], side = 3)
}

par(mfrow = c(3, 1), mar = c(3, 2, 2, 1))
boxplot.matrix(log10(aa_norm_score_bind_filtered))
boxplot.matrix(aa_norm_score_bind_filtered_quantile)
boxplot.matrix(aa_norm_score_bind_filtered_zscore)



##### having a bar next to each row indicating the avg score of the enhancer
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_loss_perGene_perEnh,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = c(1:10),
                          filename = "enhancer_score_barplot_top10_loss_calib.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
aa_score_sorted <- list()
aa_score_sorted_ind <- list()
for(i in 1:length(aa$score)){
  aa_score_sorted[[i]] <- sort(aa$score[[i]], decreasing = F)
  aa_score_sorted_ind[[i]] <- sort(aa$score[[i]], index.return=T, decreasing = F)$ix
}
aaa <- unlist(aa_score_sorted)
png("enhancer_score_bar.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 2*length(aa_norm_score)*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)        # smaller font size
par(mfrow = c(1,1), mar = c(1,1,1,1))
barplot(sqrt(aaa), border = F, horiz = T, las = 2, col =  aa_rowsidecol, width = 2)
abline(v = seq(1,10, 1), lty = 3, col = 2)
dev.off()
#### sort enhancers of each gene (rows of affinity matrix) by performance within that gene
aa_norm_score_sort <- aa_norm_score
for(i in 1:length(aa_norm_score)){
  aa_norm_score_sort[[i]] <- aa_norm_score[[i]][aa_score_sorted_ind[[i]], ]
}
aa_norm_score_sort_bind <- do.call(rbind, aa_norm_score_sort)
aa_norm_score_sort_bind_filtered <- aa_norm_score_sort_bind[, aa_tf_filtered_ind]
aa_norm_score_sort_bind_filtered_zscore <-  scale(log10(aa_norm_score_sort_bind_filtered))

aa_sc <- unlist(lapply(aa_norm_score, nrow))
aa_rowsidecol <- character(0)
for(i in 1:length(aa_sc)){
  aa_rowsidecol <- c(aa_rowsidecol, rep(col_vector[i%%length(col_vector)], aa_sc[i]))
}
aa_rowsep <-  which(!duplicated(aa_rowsidecol))

png(filename = "heatmap_binding.png",    # create PNG for the heat map        
    width = 8*500,        # 5 x 300 pixels
    height = 16*500,
    res = 500,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
heatmap.2(x = aa_norm_score_sort_bind_filtered_zscore,
          Rowv = F, Colv = T,
          dendrogram = "column",
          rowsep = aa_rowsep-1, sepwidth = c(4,4), sepcolor = "green",
          RowSideColors = aa_rowsidecol,
          trace="none")
dev.off()
#####################################################################################################################################################
# compute rank correlation for each gene and each TF: between average performance of enhacners, and affinity of enhancers for that TF

aaa_max_consLLR <- numeric(length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward))
for(i in 1:length(aaa_max_consLLR)){
  aaa_max_consLLR[i] <- ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward[[i]][1, 2]
}
names(aaa_max_consLLR) <- names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward)

aaa_max_consLR <- unlist(lapply(aaa_max_consLLR, exp))

aa_norm_score <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score
for(i in 1:length(aa_norm_score)){
  for(j in 1:nrow(aa_norm_score[[i]])){
    aa_norm_score[[i]][j, ] <- aa_norm_score[[i]][j, ] * aaa_max_consLR
  }
}
aa_norm_score_bind <- do.call(rbind, aa_norm_score)
aa_score <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_loss_perGene_perEnh,
                                enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                                model_index = c(1:10),
                                filename = "enhancer_score_barplot_top10_loss_calib.png",
                                loss = T,
                                real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)

aa_score_aff <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_accuracy_perGene_perEnh,
                                enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                                model_index = c(1:10),
                                filename = "enhancer_score_barplot_top10_accuracy_calib.png",
                                loss = F,
                                real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
aa_score_aff_max <- unlist(lapply(aa_score_aff$score, max))
aa_score_aff_max <- 2 * aa_score_aff_max -1
aa_tf_gene_correlation <- matrix(nrow = length(aa_norm_score), ncol = ncol(aa_norm_score[[1]]))
rownames(aa_tf_gene_correlation) <- names(aa_norm_score)
colnames(aa_tf_gene_correlation) <- colnames(aa_norm_score[[1]])
for(i in 1:nrow(aa_tf_gene_correlation)){
  for(j in 1:ncol(aa_tf_gene_correlation)){
    aa_x <- aa_score_aff$score[[i]]
    aa_y <- aa_norm_score[[i]][, j]
    aa_tf_gene_correlation[i, j] <- cor(x = aa_x, y = aa_y, method = "spearman")
  }
}
aa_tf_gene_correlation <- cbind(aa_tf_gene_correlation, aa_score_aff_max)
colnames(aa_tf_gene_correlation)[ncol(aa_tf_gene_correlation)] <- "accuracy"
aa_nu_enh<-  unlist(lapply(aa_score$score, length))
aa_xx <- mapply(FUN = paste, rownames(aa_tf_gene_correlation), as.character(aa_nu_enh), sep="_")
rownames(aa_tf_gene_correlation) <- aa_xx
png(filename = "heatmap_cor_per_Aff.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # smaller font size
aa_tf_gene_correlation[is.na(aa_tf_gene_correlation)] <- 0
heatmap.2(x = aa_tf_gene_correlation,
          Rowv = T, Colv = T,
          dendrogram = "both",
          #rowsep = aa_rowsep-1, sepwidth = c(4,4), sepcolor = "green",
          #RowSideColors = aa_rowsidecol,
          trace="none", na.rm = T)
dev.off()
# repeat but just for the filtered list of TFs:


png(filename = "heatmap_cor_accu_Aff_filtered.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # smaller font size
aa_tf_gene_correlation[is.na(aa_tf_gene_correlation)] <- 0
heatmap.2(x = aa_tf_gene_correlation[, c(aa_tf_filtered_ind, ncol(aa_tf_gene_correlation))],
          Rowv = T, Colv = T,
          dendrogram = "both",
          #rowsep = aa_rowsep-1, sepwidth = c(4,4), sepcolor = "green",
          #RowSideColors = aa_rowsidecol,
          trace="none", na.rm = T)
dev.off()

#####################################################################################################################
#######################################                                      ########################################
#######################################      in silico TF KD evaluation      ########################################
#######################################                                      ########################################
#####################################################################################################################
######## setting TFs concentrations to zero one by one
aa <- TF_KD_evaluator(my_optimized_par_mat = Sim_Ann_148_restart_parameters,
                      .my_data=ER_52_opt_simAnn_input,
                      .input_gene_expMat=my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                      .round_to_int = T,
                      TF_index=c(1:19),
                      WT_evaluation=numeric(0))
Sim_Ann_148_restart_TF_KD <- aa
Sim_Ann_148_restart_all_models_eval <- aa$WT_eval

aa <- Sim_Ann_148_restart_TF_KD$loss_change_per_model
aaa <- do.call(rbind, aa)
rowMeans(aaa)
boxplot.matrix(t(aaa), las = 2, main = "(sum_loss(WT) - sum_loss(KD))/nu_genes")
abline(h = seq(-6, 2, 0.5), col = 2, lty = 8)
aa <- Sim_Ann_148_restart_TF_KD$accu_change_per_model
aaa <- do.call(rbind, aa)
rowMeans(aaa)
boxplot.matrix(t(aaa), las = 2, main = "(sum_acc(WT) - sum_acc(KD))/nu_genes")
abline(h = seq(-0.1, 0.25, 0.05), col = 2, lty = 8)
aa <- Sim_Ann_148_restart_TF_KD$loss_change_per_gene
aaa <- do.call(rbind, aa)
rowMeans(aaa)
boxplot.matrix(t(aaa), las = 2, main = "(sum_loss(WT) - sum_loss(KD))/nu_models", outline = T)
abline(h = seq(-7, 2, 0.5), col = 2, lty = 8)

aa <- Sim_Ann_148_restart_TF_KD$accu_change_per_gene
aaa <- do.call(rbind, aa)
rowMeans(aaa)
boxplot.matrix(t(aaa), las = 2, main = "(sum_accu(WT) - sum_accu(KD))/nu_models")
abline(h = seq(-0.1, 0.4, 0.05), col = 2, lty = 8)


#####################################################################################################################
#######################################                                      ########################################
#######################################    CREATING INPUT FOR GEMSTAT RUN    ########################################
#######################################                                      ########################################
#####################################################################################################################
# exp 1
#write sequence
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_SimAnnFiltered <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene
for(i in 1:length(aa_enhancer_list)){
  ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_SimAnnFiltered$Chopped_Seq[[i]] <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_SimAnnFiltered$Chopped_Seq[[i]][sort(aa_enhancer_list[[i]])]
  ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_SimAnnFiltered$Chopped_Score[[i]] <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_SimAnnFiltered$Chopped_Score[[i]][sort(aa_enhancer_list[[i]]), ]
  ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_SimAnnFiltered$Chopped_GRanges[[i]] <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_SimAnnFiltered$Chopped_GRanges[[i]][sort(aa_enhancer_list[[i]])]
}
WriteFastaOfBag(write_from_list = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_SimAnnFiltered$Chopped_Seq,
                output.File.Name = "testListSeq")

# write expression
aa_map_c <- unlist(lapply(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene_SimAnnFiltered$Chopped_Seq, length))
aa_map = integer(0)
for(i in 0:(length(aa_map_c)-1)){
  aa_map <- c(aa_map, rep(i, aa_map_c[i+1]))
}
ExpressionWriter(expressionMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                 from_log_fold_change = T,
                 output.File.Name = "log_enh_test",
                 multi_enhancer_map = aa_map)
# write motifs
MotifWriter(motif.List = TF.motifs.Shrinked.count, output.File.Name = "motifs")
# write TF expression
ExpressionWriter(expressionMat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01,
                 output.File.Name = "TF_exp")
plotExpression(expMat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01, .dendrogram = "none", .Rowv = F, .Colv = F,filename = "TF_exp_dJun_ER01.png", colorPoints = c(-1, 0.2, 0.40, 0.70, 1.1))

#write mapping file
for(i in 1:length(aa_map)){
  cat(c(aa_map[i]), file = "my_grouping_real", sep = "\n", append = T)
}
#write control/exp file
aa_cnt_trt <- integer(0)
for(i in 0:(ncol(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)-1)){
  aa_cnt_trt <- c(aa_cnt_trt, rep(i, 2))
}
for(i in 1:length(aa_cnt_trt)){
  cat(c(aa_cnt_trt[i]), file = "my_ctrl_treat", sep = "\n", append = T)
}
#####################################################################################################################
# exp 2
# write sequence
# get the enhancers that do better than expected by chance, by at least one standard deviation
# if no such enhancers exists for that gene: take the ones that do better than expectation
# if no enhancer with better than random expectation performance was not found, remove that gene from data set
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene

aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_148_restart_loss_perGene_perEnh,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = c(1:10),
                          filename = "enhancer_score_barplot_top10_loss_calib.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)

aa_enhancer_list <- list()
for(i in 1:length(aa$score)){
  aa_wh <- which(aa$score[[i]] <= (aa$expected_score_shuffling[i] - aa$expected_std_shuffling[i]))
  if(length(aa_wh) > 10){
    aa_wh_sortind <- sort(aa$score[[i]][aa_wh], decreasing = F, index.return=T)$ix
    aa_wh <- aa_wh[aa_wh_sortind[1:10]]
  }
  if(length(aa_wh) == 0){
    aa_wh <- which(aa$score[[i]] <= aa$expected_score_shuffling[i])
    if(length(aa_wh) == 0){
      print(i)
      print("bad gene")
    }else{
      aa_wh_sortind <- sort(aa$score[[i]][aa_wh], decreasing = F, index.return=T)$ix
      if(length(aa_wh) > 3){
        aa_wh <- aa_wh[aa_wh_sortind[1:3]]
      }
    } 
  }
  aa_enhancer_list[[i]] <- unlist(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq[[i]])[aa_wh]
}
names(aa_enhancer_list) <- names(aa$expected_score_shuffling)

aa_size <- unlist(lapply(aa_enhancer_list, length))
aa_enhancer_list <- aa_enhancer_list[aa_size > 0]
WriteFastaOfBag(write_from_list = aa_enhancer_list,
                output.File.Name = "regElem_10max")

# write expression
my_exp_mat <- my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[aa_size > 0,]
aa_map_c <- unlist(lapply(aa_enhancer_list, length))
aa_map = integer(0)
for(i in 0:(length(aa_map_c)-1)){
  aa_map <- c(aa_map, rep(i, aa_map_c[i+1]))
}
ExpressionWriter(expressionMat = my_exp_mat,
                 from_log_fold_change = T,
                 output.File.Name = "log_enh_test_exp2",
                 multi_enhancer_map = aa_map)
# write motifs
MotifWriter(motif.List = TF.motifs.Shrinked.count, output.File.Name = "motifs")
# write TF expression
ExpressionWriter(expressionMat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01,
                 output.File.Name = "TF_exp")
plotExpression(expMat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01, .dendrogram = "none", .Rowv = F, .Colv = F,filename = "TF_exp_dJun_ER01.png", colorPoints = c(-1, 0.2, 0.40, 0.70, 1.1))

#write mapping file
for(i in 1:length(aa_map)){
  cat(c(aa_map[i]), file = "my_grouping_real_exp2", sep = "\n", append = T)
}
#write control/exp file
aa_cnt_trt <- integer(0)
for(i in 0:(ncol(my_exp_mat)-1)){
  aa_cnt_trt <- c(aa_cnt_trt, rep(i, 2))
}
for(i in 1:length(aa_cnt_trt)){
  cat(c(aa_cnt_trt[i]), file = "my_ctrl_treat", sep = "\n", append = T)
}
###########################################################################################################################
# exp3
# just 2 genes and 5 TFs, for experimenting with GEMSTAT
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_weighted_148_restart_loss_perGene_perEnh,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = c(1:10),
                          filename = "enhancer_score_barplot_top10_loss_calib.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)

aa_enhancer_list <- list()
for(i in 1:length(aa$score)){
  aa_wh <- which(aa$score[[i]] <= (aa$expected_score_shuffling[i] - aa$expected_std_shuffling[i]))
  if(length(aa_wh) > 3){
    aa_wh_sortind <- sort(aa$score[[i]][aa_wh], decreasing = F, index.return=T)$ix
    aa_wh <- aa_wh[aa_wh_sortind[1:3]]
  }
  if(length(aa_wh) == 0){
    aa_wh <- which(aa$score[[i]] <= aa$expected_score_shuffling[i])
    if(length(aa_wh) == 0){
      print(i)
      print("bad gene")
    }else{
      aa_wh_sortind <- sort(aa$score[[i]][aa_wh], decreasing = F, index.return=T)$ix
      if(length(aa_wh) > 3){
        aa_wh <- aa_wh[aa_wh_sortind[1:3]]
      }
    } 
  }
  aa_enhancer_list[[i]] <- unlist(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq[[i]])[aa_wh]
}
names(aa_enhancer_list) <- names(aa$expected_score_shuffling)

aa_size <- unlist(lapply(aa_enhancer_list, length))
aa_enhancer_list <- aa_enhancer_list[aa_size > 0]
WriteFastaOfBag(write_from_list = aa_enhancer_list[1:3],
                output.File.Name = "reg2Elem_3max")

# write expression
my_exp_mat <- my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[1:3,]
aa_map_c <- unlist(lapply(aa_enhancer_list[1:3], length))
aa_map = integer(0)
for(i in 0:(length(aa_map_c)-1)){
  aa_map <- c(aa_map, rep(i, aa_map_c[i+1]))
}
ExpressionWriter(expressionMat = my_exp_mat,
                 from_log_fold_change = T,
                 output.File.Name = "log_enh_test_exp3",
                 multi_enhancer_map = aa_map)
# write motifs
MotifWriter(motif.List = TF.motifs.Shrinked.count[1:5], output.File.Name = "motifs_5")
# write TF expression
ExpressionWriter(expressionMat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01[1:5,],
                 output.File.Name = "TF_exp_5")
plotExpression(expMat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01, .dendrogram = "none", .Rowv = F, .Colv = F,filename = "TF_exp_dJun_ER01.png", colorPoints = c(-1, 0.2, 0.40, 0.70, 1.1))

#write mapping file
for(i in 1:length(aa_map)){
  cat(c(aa_map[i]), file = "my_grouping_real_exp3", sep = "\n", append = T)
}
#write control/exp file
aa_cnt_trt <- integer(0)
for(i in 0:(ncol(my_exp_mat)-1)){
  aa_cnt_trt <- c(aa_cnt_trt, rep(i, 2))
}
for(i in 1:length(aa_cnt_trt)){
  cat(c(aa_cnt_trt[i]), file = "my_ctrl_treat_3", sep = "\n", append = T)
}
# write weights
GEMSTAT_weight_writer(expmat =my_exp_mat ,
                      "third_weights",
                      multi_enhancer_map = aa_map,
                      fold_change_to_real = T)
###########################################################################################################################
##############################                                       ######################################################
##############################         Greedy search ANALYSIS        ######################################################
##############################                                       ######################################################
###########################################################################################################################
aa <- getwd()
setwd("hal_job_greedy_148_5rep/")
GreedyJobConstructor(script_name="ER_Optimizer_Greedy_hal.R",
                     nu_jobs = 148,
                     file_name="ER_Optimizer_Greedy_hal_job",
                     nu_genes=52,
                     repeat_rounds=5)

#create the submit file creator
cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_script_greedy", sep="")
for(i in 1:260){
  cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
        paste("ER_Optimizer_Greedy_hal_job", i, sep = "_"),
        nu_jobs,
        paste("tmp_ER_Optimizer_Greedy_hal_job", i, sep = "_"),
        ">",
        paste("ER_Optimizer_Greedy_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_script_greedy",
      sep=" ",
      append = T)
}
#making them executable
cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_greedy_submit_exec", sep="")
for(i in 1:260){
  cat(c("chmod +x",
        paste("ER_Optimizer_Greedy_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_greedy_submit_exec",
      sep=" ",
      append = T)
}
#create the dagman file
DAGmanConstructor(submit_prefix = "ER_Optimizer_Greedy_hal_job",
                  start_ind=1,
                  end_ind=260,
                  filename="my_DagMan")
setwd(aa)

step_counter <- 1
nu_jobs <- 148
nu_genes <- 52
save.image(file = "Inputs_greedy.RData")
##############################################################################################################
####################################################                             #############################
####################################################   Reading Greedy results    #############################
####################################################                             #############################
##############################################################################################################
Greedy_148_restart_results <- list()
for(i in 1:148){
  load(paste0("greedy_results/Optim_Greedy_", i, ".RData"))
  Greedy_148_restart_results[[i]] <- GreedyRes[[length(GreedyRes)]]
  remove(GreedyRes)
}

Greedy_148_restart_enh <- lapply(Greedy_148_restart_results, "[[", 1)
names(Greedy_148_restart_enh) <- c(1:length(Greedy_148_restart_enh))

Greedy_148_restart_par <- lapply(Greedy_148_restart_results, "[[", 3)
names(Greedy_148_restart_par) <- c(1:length(Greedy_148_restart_par))


Greedy_148_restart_obj <- lapply(Greedy_148_restart_results, "[[", 4)
names(Greedy_148_restart_obj) <- c(1:length(Greedy_148_restart_obj))



Greedy_148_restart_obj <- unlist(Greedy_148_restart_obj)
Greedy_148_restart_enh <- do.call(rbind, Greedy_148_restart_enh)
rownames(Greedy_148_restart_enh) <- names(Greedy_148_restart_obj)
Greedy_148_restart_par <- do.call(rbind, Greedy_148_restart_par)
rownames(Greedy_148_restart_par) <- names(Greedy_148_restart_obj)

#look at chosen enhancers
aaa <- unlist(lapply(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_1500_500_perGene$Chopped_Score, nrow))
Enhancer_vote_barplot(Enhancer_matrix = Greedy_148_restart_enh, nu_enh_per_gene = aaa)

#create the prediction matrices
setwd("Greedy_Results_plot/All_148_predictions")
for(i in 1:nrow(aa_par_last_sort)){
  aa_pred_expMat[[i]] <- objFuncEval(my_optimized_par = aa_par_last_sort[i, ],
                                     my_data = ER_52_opt_simAnn_input,
                                     input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                     round_to_int = T)
  
  plotExpression(expMat = aa_pred_expMat[[i]]$Predicted_expression_mat, .dendrogram = "none", .Rowv = F, .Colv = F,filename = paste0("sim_ann_",i, "_" ,rownames(aa_par_last_sort)[i], ".png"))
}

########################################################################################################################
##############################################################################################################


a <- c(0:36)
for(i in 1:length(a)){
  cat(c(a[i]), file = "my_grouping_1", sep = "\n", append = T)
}


aa <- c(rep(0,4), rep(1, 4), rep(2, 9), rep(3,3), rep(4, 10), rep(5, 7))

for(i in 1:length(aa)){
  cat(c(aa[i]), file = "my_grouping_2", sep = "\n", append = T)
}
######################
aaabind <- do.call(rbind,ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score)

boxplot(aaabind)
aaabind
apply(aaabind, 2, median)
aatest <- aa_par_last_sort[1:50,1:19]
boxplot.matrix(aa_par_last_sort[1:50,1:19], las = 2, outline = F)
for(i in 1:nrow(aatest)){
  aatest[i, ] <- aa_par_last_sort[i,1:19] * apply(aaabind, 2, median)
}
boxplot.matrix(aatest, las = 2, outline= F)
boxplot.matrix(log(abs(aatest)), las = 2, outline= F)



boxplot.matrix(aaabind, las = 2)


boxplot.matrix(log(aaabind), las = 2, ylim = c(-30, 5))


aaa <- unlist(lapply(strsplit(colnames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), split ="_"), "[[", 1))
aaa_col_color <- character(length(aaa))
for(i in 1:length(unique(aaa))){
  aaa_col_color[aaa %in% unique(aaa)[i]] <- col_vector[i]
}
plotExpression(expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, .dendrogram = "both", .Rowv = T, .Colv = T,filename = "real_exp_2.png", .distfun = my_dist_fuc, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5), .ColSideColors = aaa_col_color)

########################################################################################################################
###################################################                                     ################################
###################################################    weighted simulated annealing     ################################
###################################################                                     ################################
########################################################################################################################
#####################################################################################################################
# weighted sim ann run:
#creating input for hal runs
# normalize using max LR
aaa_max_consLLR <- numeric(length(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward))
for(i in 1:length(aaa_max_consLLR)){
  aaa_max_consLLR[i] <- ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward[[i]][1, 2]
}
names(aaa_max_consLLR) <- names(ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore$`10057.1`$Forward)

aaa_max_consLR <- unlist(lapply(aaa_max_consLLR, exp))

aa_norm_score <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Score
for(i in 1:length(aa_norm_score)){
  for(j in 1:nrow(aa_norm_score[[i]])){
    aa_norm_score[[i]][j, ] <- aa_norm_score[[i]][j, ] * aaa_max_consLR
  }
}

ER_52_opt_simAnn_2_input <- ER_52_opt_simAnn_input
ER_52_opt_simAnn_2_input$Affinity_scores <- aa_norm_score
ER_52_opt_simAnn_input$Affinity_scores$`10057`
ER_52_opt_simAnn_2_input$Affinity_scores$`10057`
aa_w <-which(unlist(lapply(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Seq, length)) > 40)
my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Score$`11045`
aa <- my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Score[[6]]
cor.plot(cor(t(aa))) 
heatmap.2(cor(t(aa)), Rowv = F, 
          Colv = F, na.rm = T, dendrogram = "none", trace = "none", 
          breaks = seq(-1, 1, length.out = 21), margins = c(10,10),
          col = colorspace::diverge_hcl(20), main = "enhancer_aff_Cor")

heatmap.2(log10(aa), Rowv = T, 
          Colv = T, na.rm = T, dendrogram = "both", trace = "none", 
          breaks = seq(0, 5, length.out = 21), margins = c(10,10),
          col = colorspace::terrain_hcl(20), main = "enhancerLLR")

aadist <- as.matrix(dist(aa, method = "maximum"))
heatmap.2(aadist, Rowv = F, 
          Colv = F, na.rm = T, dendrogram = "none", trace = "none", 
          breaks = seq(0, 3.75, length.out = 21), margins = c(10,10),
          col = colorspace::terrain_hcl(20), main = "enhancerLLR")
sum(aadist > 0 & aadist < 0.5)
####################
# write hal jobs for different temperatures
aa_wd <- getwd()
setwd("Hal_simAnn_weighted_jobs/")
SimAnnStepJobConstructor(script_name="ER_Optimizer_SimAnn_weighted_hal.R",
                         nu_jobs=148,
                         file_name="ER_Optimizer_SimAnn_weighted_hal_job",
                         start_temp=1,
                         min_temp=1e-6,
                         alPha=0.95,
                         alpha_Step=1)

##############################################################
step_counter <- 1
nu_jobs <- 148
##############################################################

cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_script", sep="")
for(i in 1:270){
  cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
        paste("ER_Optimizer_SimAnn_weighted_hal_job", i, sep = "_"),
        nu_jobs,
        paste("tmp_ER_Optimizer_SimAnn_weighted_hal_job", i, sep = "_"),
        ">",
        paste("ER_Optimizer_SimAnn_weighted_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_script",
      sep=" ",
      append = T)
}
cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_exec", sep="")
for(i in 1:270){
  cat(c("chmod +x",
        paste("ER_Optimizer_SimAnn_weighted_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_exec",
      sep=" ",
      append = T)
}
##############
DAGmanConstructor(submit_prefix = "ER_Optimizer_SimAnn_weighted_hal_job",
                  start_ind=1,
                  end_ind=270,
                  #order = c(rep(1,3), rep(2,3), rep(3,3), 4),
                  filename="my_DagMan_1_270")
DAGmanConstructor(submit_prefix = "ER_Optimizer_SimAnn_weighted_hal_job",
                  start_ind=190,
                  end_ind=270,
                  #order = c(rep(1,3), rep(2,3), rep(3,3), 4),
                  filename="my_DagMan_190_270")
setwd(aa_wd)
############################################################################################################
# setting up weighted sim ann for 172 genes
ER_172_opt_simAnn_input <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                         TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01,
                         ChoppedEnhancerScores = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Score) 
# after changing the TF expression matrix to be base on raw fold changes
ER_172_opt_simAnn_input_BORAW <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                                                    TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times,
                                                    ChoppedEnhancerScores = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Score) 

####################
# write hal jobs for different temperatures
aa_wd <- getwd()
setwd("Hal_simAnn_weighted_jobs_2/")
SimAnnStepJobConstructor(script_name="ER_Optimizer_SimAnn_weighted_hal.R",
                         nu_jobs=148,
                         file_name="ER_Optimizer_SimAnn_weighted_hal_job",
                         start_temp=1,
                         min_temp=1e-6,
                         alPha=0.95,
                         alpha_Step=1)
# write hal jobs for different temperatures AFTER after changing the TF expression matrix to be base on raw fold changes
aa_wd <- getwd()
setwd("Hal_simAnn_weighted_jobs_3/")
SimAnnStepJobConstructor(script_name="ER_Optimizer_SimAnn_weighted_hal_MAX20.R",
                         nu_jobs=148,
                         file_name="ER_Optimizer_SimAnn_weighted_hal_job",
                         start_temp=1,
                         min_temp=1e-6,
                         alPha=0.95,
                         alpha_Step=1)

##############################################################
step_counter <- 1
nu_jobs <- 148
##############################################################

cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_script", sep="")
for(i in 1:270){
  cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
        paste("ER_Optimizer_SimAnn_weighted_hal_job", i, sep = "_"),
        nu_jobs,
        paste("tmp_ER_Optimizer_SimAnn_weighted_hal_job", i, sep = "_"),
        ">",
        paste("ER_Optimizer_SimAnn_weighted_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_script",
      sep=" ",
      append = T)
}
cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_exec", sep="")
for(i in 1:270){
  cat(c("chmod +x",
        paste("ER_Optimizer_SimAnn_weighted_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_exec",
      sep=" ",
      append = T)
}
#### after using raw fold changes for TF expression matrix
cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_script", sep="")
for(i in 1:270){
  cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
        paste("ER_Optimizer_SimAnn_weighted_hal_job", i, sep = "_"),
        nu_jobs,
        paste("tmp_ER_Optimizer_SimAnn_weighted_hal_job", i, sep = "_"),
        ">",
        paste("ER_Optimizer_SimAnn_weighted_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_script",
      sep=" ",
      append = T)
}
cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_exec", sep="")
for(i in 1:270){
  cat(c("chmod +x",
        paste("ER_Optimizer_SimAnn_weighted_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_exec",
      sep=" ",
      append = T)
}
##############
DAGmanConstructor(submit_prefix = "ER_Optimizer_SimAnn_weighted_hal_job",
                  start_ind=1,
                  end_ind=270,
                  #order = c(rep(1,3), rep(2,3), rep(3,3), 4),
                  filename="my_DagMan_1_270")
save(list = c("ER_172_opt_simAnn_input", 
       "nu_jobs", "obj_func_chosen_enh",
       "SimAnnAcceptProb", "SimAnnObjOptim",
       "SimAnnStepJobConstructor", "step_counter"), 
     file = "Inputs_simAnn.RData")
setwd(aa_wd)

save(list = c("ER_172_opt_simAnn_input_BORAW", 
              "nu_jobs", "obj_func_chosen_enh",
              "SimAnnAcceptProb", "SimAnnObjOptim",
              "SimAnnStepJobConstructor", "step_counter"), 
     file = "Inputs_simAnn.RData")

## 
# test run
aa_SimAnnRes <- list()
aa_time <- proc.time()
aa_SimAnnRes[[1]] <- SimAnnObjOptim(.my_data = ER_172_opt_simAnn_input,
                                    initial_chosen_enh=numeric(0),
                                    initial_par=numeric(0),
                                    starting_temp=1,
                                    min_temp=0.9,
                                    iter_per_temp=1,
                                    ALPHA=0.95) 
proc.time() - aa_time
plot(unlist(aa_SimAnnRes[[1]]$obj_value))
# after using Raw fold change for TF expression
aa_SimAnnRes <- list()
aa_time <- proc.time()
aa_SimAnnRes[[1]] <- SimAnnObjOptim(.my_data = ER_172_opt_simAnn_input_BORAW,
                                    initial_chosen_enh=numeric(0),
                                    initial_par=numeric(0),
                                    starting_temp=1,
                                    min_temp=0.95,
                                    iter_per_temp=1,
                                    ALPHA=0.95,
                                    .intercept_pergene = T, 
                                    max_gene_per_iter = 5, 
                                    .sig_over_loss=T,
                                    my_algorithm = "NLOPT_LN_COBYLA") 
proc.time() - aa_time
plot(unlist(aa_SimAnnRes[[1]]$obj_value))


############################################################################################################
# create another dagman for a subset of the genes becasue it takes really long for all genes


table(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172)
aa1 <- apply(X = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
             MARGIN = 1, FUN = (function(x) sum(x != 0 & !is.na(x))))
hist(aa1, breaks = 10)
sum(aa1 >= 10)
my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_top63 <- my_CommonDifExpMat_16_ERassoc_gte5nzAll_172[aa1 >= 10, ]
plotExpression(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_top63, filename = "genes_gte5nonzero_172_top63_heatmap.png")

ER_63_opt_simAnn_input <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_top63,
                                              TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01,
                                              ChoppedEnhancerScores = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Score[aa1 >= 10]) 
ER_63_opt_simAnn_input_test109 <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172[aa1 < 10, ],
                                                     TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01,
                                                     ChoppedEnhancerScores = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Score[aa1 < 10]) 
####################
# write hal jobs for different temperatures
aa_wd <- getwd()
setwd("Hal_simAnn_weighted_jobs_63_genes/")
SimAnnStepJobConstructor(script_name="ER_Optimizer_SimAnn_weighted_hal.R",
                         nu_jobs=148,
                         file_name="ER_Optimizer_SimAnn_weighted_hal_job",
                         start_temp=1,
                         min_temp=1e-6,
                         alPha=0.95,
                         alpha_Step=1)
##############################################################
step_counter <- 1
nu_jobs <- 148
##############################################################

cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_script", sep="")
for(i in 1:270){
  cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
        paste("ER_Optimizer_SimAnn_weighted_hal_job", i, sep = "_"),
        nu_jobs,
        paste("tmp_ER_Optimizer_SimAnn_weighted_hal_job", i, sep = "_"),
        ">",
        paste("ER_Optimizer_SimAnn_weighted_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_script",
      sep=" ",
      append = T)
}
cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_exec", sep="")
for(i in 1:270){
  cat(c("chmod +x",
        paste("ER_Optimizer_SimAnn_weighted_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_exec",
      sep=" ",
      append = T)
}
##############
DAGmanConstructor(submit_prefix = "ER_Optimizer_SimAnn_weighted_hal_job",
                  start_ind=1,
                  end_ind=270,
                  #order = c(rep(1,3), rep(2,3), rep(3,3), 4),
                  filename="my_DagMan_1_270")
save(list = c("ER_63_opt_simAnn_input", 
              "nu_jobs", "obj_func_chosen_enh",
              "SimAnnAcceptProb", "SimAnnObjOptim",
              "SimAnnStepJobConstructor", "step_counter"), 
     file = "Inputs_simAnn.RData")
setwd(aa_wd)
## 
# test run
aa_SimAnnRes_63 <- list()
aa_SimAnnRes_63[[1]] <- SimAnnObjOptim(.my_data = ER_63_opt_simAnn_input,
                                    initial_chosen_enh=numeric(0),
                                    initial_par=numeric(0),
                                    starting_temp=1,
                                    min_temp=0.9,
                                    iter_per_temp=1,
                                    ALPHA=0.95) 
############################################################################################################
# reading the results of weighted sim ann run for 63 genes
aa_cur_dir <- getwd()
setwd("Hal_simAnn_weighted_jobs_63_genes/Results")

# Reading Simulated annealing results
Sim_Ann_weighted_148_restart_63_genes_results <- list()
for(i in 1:148){
  load(paste0("Optim_sim_Ann_", i, ".RData"))
  Sim_Ann_weighted_148_restart_63_genes_results[[i]] <- SimAnnRes
  remove(SimAnnRes)
}
setwd(aa_cur_dir)

### what is the best model, plot the losses in each of the 148 runs
aa_res_hold <- list()
aa_res_hold_agg <- list()
for(i in 1:length(Sim_Ann_weighted_148_restart_63_genes_results)){
  aa_res_hold[[i]] <- lapply(Sim_Ann_weighted_148_restart_63_genes_results[[i]], "[[", 3)
  aa_res_hold[[i]] <- lapply(aa_res_hold[[i]], unlist)
  names(aa_res_hold[[i]]) <- c(1:length(aa_res_hold[[i]]))
  for(j in 1:length(aa_res_hold[[i]])){
    names(aa_res_hold[[i]][[j]]) <- c(1:length(aa_res_hold[[i]][[j]]))
  }
  aa_res_hold_agg[[i]] <- do.call(c, aa_res_hold[[i]])
}
plot(aa_res_hold_agg[[3]])
# find the maximum performing model per restart
aa_res_hold_agg_sort <- lapply(aa_res_hold_agg, sort, decreasing = F)
aa_res_hold_best <- unlist(lapply(aa_res_hold_agg_sort, "[", 1))
aa_res_hold_best_names <- strsplit(names(aa_res_hold_best), "\\.")
Sim_Ann_weighted_148_restart_63_genes_results_best <- list()
Sim_Ann_weighted_148_restart_63_genes_results_best$Enhancer_index <- matrix(nrow = length(aa_res_hold_best_names),
                                                                            ncol = length(Sim_Ann_weighted_148_restart_63_genes_results[[1]][[1]]$Enhancer_index[[1]]))
Sim_Ann_weighted_148_restart_63_genes_results_best$parameters <- matrix(nrow = length(aa_res_hold_best_names),
                                                                        ncol = length(Sim_Ann_weighted_148_restart_63_genes_results[[1]][[1]]$parameters[[1]]))
Sim_Ann_weighted_148_restart_63_genes_results_best$obj_value <- numeric(length = length(aa_res_hold_best_names))
for(i in 1:length(aa_res_hold_best_names)){
  Sim_Ann_weighted_148_restart_63_genes_results_best$Enhancer_index[i, ] <- Sim_Ann_weighted_148_restart_63_genes_results[[i]][[as.integer(aa_res_hold_best_names[[i]][1])]][[1]][[as.integer(aa_res_hold_best_names[[i]][2])]]
  Sim_Ann_weighted_148_restart_63_genes_results_best$parameters[i, ] <- Sim_Ann_weighted_148_restart_63_genes_results[[i]][[as.integer(aa_res_hold_best_names[[i]][1])]][[2]][[as.integer(aa_res_hold_best_names[[i]][2])]]
  Sim_Ann_weighted_148_restart_63_genes_results_best$obj_value[i] <- Sim_Ann_weighted_148_restart_63_genes_results[[i]][[as.integer(aa_res_hold_best_names[[i]][1])]][[3]][[as.integer(aa_res_hold_best_names[[i]][2])]]
}

#sort models based on performance
aa <- sort(Sim_Ann_weighted_148_restart_63_genes_results_best$obj_value, decreasing = F, index.return=T)$ix
Sim_Ann_weighted_148_restart_63_genes_results_best$Enhancer_index <- Sim_Ann_weighted_148_restart_63_genes_results_best$Enhancer_index[aa, ]
Sim_Ann_weighted_148_restart_63_genes_results_best$parameters <- Sim_Ann_weighted_148_restart_63_genes_results_best$parameters[aa, ]
Sim_Ann_weighted_148_restart_63_genes_results_best$obj_value <- Sim_Ann_weighted_148_restart_63_genes_results_best$obj_value[aa]
rownames(Sim_Ann_weighted_148_restart_63_genes_results_best$Enhancer_index) <- aa
rownames(Sim_Ann_weighted_148_restart_63_genes_results_best$parameters) <- aa
names(Sim_Ann_weighted_148_restart_63_genes_results_best$obj_value) <- aa

## evaulate the models
Sim_Ann_weighted_148_restart_63_genes_results_best_eval<- list()
for(i in 1:nrow(Sim_Ann_weighted_148_restart_63_genes_results_best$parameters)){
  Sim_Ann_weighted_148_restart_63_genes_results_best_eval[[i]] <- objFuncEval(my_optimized_par = Sim_Ann_weighted_148_restart_63_genes_results_best$parameters[i, ],
                                                                             my_data = ER_63_opt_simAnn_input,
                                                                             input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_top63,
                                                                             round_to_int = T)
}
table(Sim_Ann_weighted_148_restart_63_genes_results_best_eval[[1]]$Predicted_expression_mat)
### are the models better than random, create confusion mats
Sim_Ann_weighted_148_restart_63_genes_confusion <- confusion_constructor(predictionMat_List=lapply(Sim_Ann_weighted_148_restart_63_genes_results_best_eval, "[[", 1),
                                                                realMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_top63)
names(Sim_Ann_weighted_148_restart_63_genes_confusion) <- names(Sim_Ann_weighted_148_restart_63_genes_results_best$obj_value)

Sim_Ann_weighted_148_restart_63_genes_random_pred <- CreateRandomPrediction(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_top63, num = 148)
par(mfrow = c(1,1), mar = c(6,6,6,6))
plot(c(unlist(lapply(Sim_Ann_weighted_148_restart_63_genes_results_best_eval,"[[", 5)),
       Sim_Ann_weighted_148_restart_63_genes_random_pred$Precision),
     ylim = c(0, 1),
     col = c(rep(3, 148), rep(2, 148)),
     ylab = "Accuracy",
     xlab = "model index",
     cex = 0.8,
     pch = 16,
     main="Linear Model Train Accuracy")
legend(200, 1, legend=c("Sim_Ann_weighted", "shuffled"),
       col=c("green", "red"), pch=16, cex=1.1,
       bty = "n", x.intersp = 0.1,
       y.intersp = 0.6, pt.cex = 1.1, box.col = 1)

# Draw heatmap of the ensemble
aa <- PerformanceHeattMap_General(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_top63,
                                  prediction_mat_list = lapply(Sim_Ann_weighted_148_restart_63_genes_results_best_eval, "[[", 1),
                                  .Colv = F, .Rowv = T, .dendrogram = "row",exportplot = T,
                                  filename = "PerformanceHeatMap_weighted_63genes.png",
                                  .RowSideColors = character(0))
############## scoring each enhancer
aa_sort <- sort(unlist(lapply(Sim_Ann_weighted_148_restart_63_genes_results_best_eval,"[[", 5)), decreasing = T, index.return=T)$ix
aa1 <- apply(X = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
             MARGIN = 1, FUN = (function(x) sum(x != 0 & !is.na(x))))
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = lapply(Sim_Ann_weighted_148_restart_63_genes_results_best_eval,"[[", 7),
                          enhancer_Granges = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_GRanges[aa1 >= 10],
                          model_index = aa_sort[1],
                          filename = "weightedModel_enhancer_score_barplot_top1_class_loss_calib_rounded.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_top63)
############## scoring the other genes : the ones not included in the linear model
Sim_Ann_weighted_148_restart_63_genes_results_best_eval_test_109<- list()
for(i in 1:nrow(Sim_Ann_weighted_148_restart_63_genes_results_best$parameters)){
  Sim_Ann_weighted_148_restart_63_genes_results_best_eval_test_109[[i]] <- objFuncEval(my_optimized_par = Sim_Ann_weighted_148_restart_63_genes_results_best$parameters[i, ],
                                                                                       my_data = ER_63_opt_simAnn_input_test109,
                                                                                       input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172[aa1 < 10,],
                                                                                       round_to_int = T, 
                                                                                       fixed_threshold_value = Sim_Ann_weighted_148_restart_63_genes_results_best_eval[[i]]$discretization_threshold)
}
### are the models better than random, create confusion mats
Sim_Ann_weighted_148_restart_63_genes_confusion_test109 <- confusion_constructor(predictionMat_List=lapply(Sim_Ann_weighted_148_restart_63_genes_results_best_eval_test_109, "[[", 1),
                                                                         realMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172[aa1 < 10,])
names(Sim_Ann_weighted_148_restart_63_genes_confusion_test109) <- names(Sim_Ann_weighted_148_restart_63_genes_results_best$obj_value)

Sim_Ann_weighted_148_restart_63_genes_random_pred_test109 <- CreateRandomPrediction(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172[aa1 < 10,], num = 148)
par(mfrow = c(1,1), mar = c(6,6,6,6))

plot(c(unlist(lapply(Sim_Ann_weighted_148_restart_63_genes_results_best_eval_test_109,"[[", 5)),
       Sim_Ann_weighted_148_restart_63_genes_random_pred_test109$Precision),
     ylim = c(0, 1),
     col = c(rep(3, 148), rep(2, 148)),
     ylab = "Accuracy",
     xlab = "model index",
     cex = 0.8,
     pch = 16, main = "linear model Test Accuracy")
legend(200, 1, legend=c("Sim_Ann_weighted", "shuffled"),
       col=c("green", "red"), pch=16, cex=1.1,
       bty = "n", x.intersp = 0.1,
       y.intersp = 0.6, pt.cex = 1.1, box.col = 1)

# Gather the results for all genes together, only for models with decent performance on test data.
Sim_Ann_weighted_148_restart_63_genes_results_best_eval_172<- list()
for(i in 1:nrow(Sim_Ann_weighted_148_restart_63_genes_results_best$parameters)){
  Sim_Ann_weighted_148_restart_63_genes_results_best_eval_172[[i]] <- objFuncEval(my_optimized_par = Sim_Ann_weighted_148_restart_63_genes_results_best$parameters[i, ],
                                                                                       my_data = ER_172_opt_simAnn_input,
                                                                                       input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                                                                                       round_to_int = T, 
                                                                                       fixed_threshold_value = Sim_Ann_weighted_148_restart_63_genes_results_best_eval[[i]]$discretization_threshold)
}
max(Sim_Ann_weighted_148_restart_63_genes_random_pred_test109$Precision)
aa <- which(unlist(lapply(Sim_Ann_weighted_148_restart_63_genes_results_best_eval_test_109,"[[", 5)) > 0.5)
hist(unlist(lapply(Sim_Ann_weighted_148_restart_63_genes_results_best_eval_test_109,"[[", 5)))

aa_dat <- data.frame( x=as.numeric(unlist(lapply(Sim_Ann_weighted_148_restart_63_genes_results_best_eval_test_109,"[[", 5))), 
                      above_50= as.numeric(unlist(lapply(Sim_Ann_weighted_148_restart_63_genes_results_best_eval_test_109,"[[", 5))) >= 0.5 )
qplot(x,data=aa_dat,geom="histogram",
      fill=above_50, main = "Linear model test accuracy")


Sim_Ann_weighted_148_restart_63_genes_results_best_eval_172 <- Sim_Ann_weighted_148_restart_63_genes_results_best_eval_172[aa]
names(Sim_Ann_weighted_148_restart_63_genes_results_best_eval_172) <- rownames(Sim_Ann_weighted_148_restart_63_genes_results_best$parameters)[aa]
Sim_Ann_weighted_148_restart_63_genes_results_best_parameter_172 <- Sim_Ann_weighted_148_restart_63_genes_results_best$parameters[aa, ]
colnames(Sim_Ann_weighted_148_restart_63_genes_results_best_parameter_172) <- c(rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01), "intercept")
names(TF.motifs.Shrinked.hocomoco_E2F1_added.count) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01)
save(list = c("TF.motifs.Shrinked.hocomoco_E2F1_added.count",
              "TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01", 
              "Sim_Ann_weighted_148_restart_63_genes_results_best_eval_172",
              "Sim_Ann_weighted_148_restart_63_genes_results_best_parameter_172",
              "my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1",
              "my_CommonDifExpMat_16_ERassoc_gte5nzAll_172", 
              "TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added"),
     file = "to_be_added_to_ensemble_construction_workplace_172.RData")
############################################################################################################
# setting up sim ann for 172 genes only for ER as the TF

# getting TF expression mat for only ER
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times_ERonly <- matrix(nrow = 1, 
                                                                                           ncol = ncol(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times))
rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times_ERonly) <- "ESR1"
colnames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times_ERonly) <- colnames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times)
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times_ERonly[1, ]  <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times[4, ]

# getting the motif scores for only ER
my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1_FilteredScore_ERonly <- my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Score
for(i in 1:length(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1_FilteredScore_ERonly)){
  aatmp <- matrix(nrow = nrow(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1_FilteredScore_ERonly[[i]]), ncol = 1)
  rownames(aatmp) <- rownames(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1_FilteredScore_ERonly[[i]])
  colnames(aatmp) <- "ESR1"
  aatmp[, 1] <- my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1_FilteredScore_ERonly[[i]][, 4]
  my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1_FilteredScore_ERonly[[i]] <- aatmp
}
names(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1_FilteredScore_ERonly) <- names(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Score)

my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1_FilteredScore_ERonly$`100507584`

ER_172_opt_simAnn_input_BORAW_ERonly <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                                                           TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times_ERonly,
                                                           ChoppedEnhancerScores = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1_FilteredScore_ERonly) 
aaSimAnnRes <- list()
aa_parRange <- cbind(c(-20, rep(-10, 172)), c(20, rep(10, 172)))
aa_init_enh = rep(1, 172)
aa_init_par = runif(173, min = -10,max = 10)

aa_time2 <- proc.time()
aaSimAnnRes[[1]] <- SimAnnObjOptim(.my_data = ER_172_opt_simAnn_input_BORAW_ERonly,
                                    initial_chosen_enh=aa_init_enh,
                                   .intercept_pergene = T,
                                   par_range = aa_parRange,
                                    initial_par=aa_init_par,
                                    starting_temp=1,
                                    min_temp=0.95,
                                    iter_per_temp=1,
                                    ALPHA=0.95) 
proc.time() - aa_time2

local_opts <- list(algorithm = "NLOPT_LD_MMA",
                   xtol_rel=1e-04)
opts <- list( "algorithm" = "NLOPT_LN_BOBYQA"
              ,"xtol_rel"  = 1.0e-04,
              "maxeval"   = 100
              ,"local_opts" = local_opts
              # ,      
              # check_derivatives = T, 
              # check_derivatives_tol = 1e-04,
              # check_derivatives_print='all'
)
aa_time2 <- proc.time()
aa <- nloptr( x0=aa_init_par,
              eval_f=obj_func_chosen_enh,
              my_data = ER_172_opt_simAnn_input_BORAW_ERonly,
              chosen_enh_ind = aa_init_enh,
              intercept_pergene=T, 
              constrained=F,
              return_grad = F, 
              lb=aa_parRange[,1],
              ub=aa_parRange[,2],
              #eval_g_ineq=obj_func_chosen_enh_ineq,
              #     eval_g_eq=eval_g_eq,
              opts=opts)
proc.time() - aa_time2
aa$solution

aa <- nloptr.get.default.options()
aa_all_algs <- unlist(strsplit(aa$possible_values[1], split = ", "))
aaa <- strsplit(aa_all_algs, split = "_")
aa_ap <- which(unlist(lapply(aaa, "[[", 2)) == "LD" | unlist(lapply(aaa, "[[", 2)) == "GD")
aa_ap2 <- which(unlist(lapply(aaa, "[[", 2)) == "LN"| unlist(lapply(aaa, "[[", 2)) == "GN")
aa_time_list <- list()

aa_test_input <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                                                       TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_BORAW_ER10times,
                                                       ChoppedEnhancerScores = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Score)
# aaSimAnnRes <- list()
#aa_parRange <- cbind(c(-20, rep(-10, 172)), c(20, rep(10, 172)))
# aa_init_enh = rep(1, 172)

aa_start_par <- runif(n = ((ncol(aa_test_input$gene_conditon_mat) - 1)/2  + 172), 
                      min = -10, max = 10)


aa_test_enh <- integer(length(aa_test_input$Affinity_scores))
for(i in 1:length(aa_test_enh)){
  aa_test_enh[i] <- sample(x = c(1:unlist(lapply(aa_test_input$Affinity_scores, nrow))[i]), size = 1)
}

local_opts$algorithm <- "NLOPT_LD_MMA"
local_opts$xtol_rel <- 1.0e-04

aa_res_nd2 <- list()
aa_time_nd2 <- list()

for(i in 20:length(aa_ap2)){
  print(aa_all_algs[aa_ap2[i]])
  opts <- list( "algorithm" = aa_all_algs[aa_ap2[i]]
                ,"xtol_rel"  = 1.0e-04,
                "maxeval"   = 100
                ,"local_opts" = local_opts
                # ,      
                # check_derivatives = T, 
                # check_derivatives_tol = 1e-04,
                # check_derivatives_print='all'
  )
  
  aa_time <- proc.time()
  aa_res_nd2[[i]] <- nloptr( x0=aa_start_par,
                            eval_f=obj_func_chosen_enh,
                            my_data = aa_test_input,
                            chosen_enh_ind = aa_test_enh,
                            intercept_pergene=T, 
                            constrained=F,
                            return_grad = F,
                            sig_over_loss = F,
                            # lb=aa_parRange[,1],
                            # ub=aa_parRange[,2],
                            opts=opts)
  
  aa_time_nd2[[i]] <- proc.time() - aa_time
  print(aa_res_nd2[[i]])
  print(aa_time_nd2[[i]] )
}


aa_time_nd2[[20]]
aa_res_nd2[[20]]$objective

aa_obj_n <- unlist(lapply(aa_res_nd2, "[[", 17))[is.finite(unlist(lapply(aa_res_nd2, "[[", 17)))]
names(aa_res_nd2) <- aa_all_algs[aa_ap2]
unlist(lapply(aa_time_nd2, "[[", 3))[is.finite(unlist(lapply(aa_res_nd2, "[[", 17)))]
aa_start_par <- runif(n = ((ncol(aa_test_input$gene_conditon_mat) - 1)/2  + 172), 
                      min = -10, max = 10)
aa_erpar <- seq(-2, 1, length.out = 20)
for( i in 1:length(aa_erpar)){
  aa_opt_eval <- objFuncEval(my_optimized_par = c(aa_erpar[i],aa_start_par[2:173]),
                             # aa_res_nd2[[20]]$solution,
                             my_data = aa_test_input,
                             input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                             round_to_int = T,
                             chosen_enh = aa_test_enh, 
                             intercept_pergene = T,
                             sigmoid_on_conditionLoss = F,
                             consider_label_weights = T,
                             fixed_threshold_value = numeric(0))  
  print("########################################")
  print("loss")
  print(sum(aa_opt_eval$losses))
  print(paste0("ER par is : ", aa_erpar[i]) )
  print("less than -0.5")
  print(sum(aa_opt_eval$Predicted_expression_mat_realValue <= -0.5, na.rm = T))
  print("more than 0.5")
  print(sum(aa_opt_eval$Predicted_expression_mat_realValue >= 0.5, na.rm = T))
  print("between -0.5 and 0.5")
  print(sum(aa_opt_eval$Predicted_expression_mat_realValue > -0.5 & aa_opt_eval$Predicted_expression_mat_realValue < 0.5, na.rm = T))
  hist(aa_opt_eval$Predicted_expression_mat_realValue[!is.na(aa_opt_eval$Predicted_expression_mat_realValue)], breaks = 100)
  print("########################################")
  
}
aa_all_aff <- do.call(rbind,aa_test_input$Affinity_scores )

############################################################################################################
######### reading the weighted sim ann results
aa_cur_dir <- getwd()
setwd("SimAnn_Weighted_Results/")

# Reading Simulated annealing results
Sim_Ann_weighted_148_restart_results <- list()
for(i in 1:148){
  load(paste0("Optim_sim_Ann_", i, ".RData"))
  Sim_Ann_weighted_148_restart_results[[i]] <- SimAnnRes
  remove(SimAnnRes)
}
setwd(aa_cur_dir)


Sim_Ann_weighted_148_restart_results_last <- list()
for(i in 1:length(Sim_Ann_weighted_148_restart_results)){
  Sim_Ann_weighted_148_restart_results_last[[i]] <- Sim_Ann_weighted_148_restart_results[[i]][[length(Sim_Ann_weighted_148_restart_results[[i]])]]
}
##### Grab the last step of each random start in sim Ann
Sim_Ann_weighted_148_restart_results_last_last <- list()
for(i in 1:length(Sim_Ann_weighted_148_restart_results_last)){
  Sim_Ann_weighted_148_restart_results_last_last[[i]] <- list()
  Sim_Ann_weighted_148_restart_results_last_last[[i]]$Enhancer_index <- Sim_Ann_weighted_148_restart_results_last[[i]]$Enhancer_index[[length(Sim_Ann_weighted_148_restart_results_last[[i]]$Enhancer_index)]]
  Sim_Ann_weighted_148_restart_results_last_last[[i]]$parameters <- Sim_Ann_weighted_148_restart_results_last[[i]]$parameters[[length(Sim_Ann_weighted_148_restart_results_last[[i]]$parameters)]]
  Sim_Ann_weighted_148_restart_results_last_last[[i]]$obj_value <- Sim_Ann_weighted_148_restart_results_last[[i]]$obj_value[[length(Sim_Ann_weighted_148_restart_results_last[[i]]$obj_value)]]
}

Sim_Ann_weighted_148_restart_results_last_last_conc <- list()
aa <- lapply(Sim_Ann_weighted_148_restart_results_last_last, "[[", 1)
Sim_Ann_weighted_148_restart_results_last_last_conc$Enhancer_index <- do.call(rbind, aa)
aa <- lapply(Sim_Ann_weighted_148_restart_results_last_last, "[[", 2)
Sim_Ann_weighted_148_restart_results_last_last_conc$parameters  <- do.call(rbind, aa)
aa <- lapply(Sim_Ann_weighted_148_restart_results_last_last, "[[", 3)
Sim_Ann_weighted_148_restart_results_last_last_conc$obj_value  <- do.call(c, aa)

remove(Sim_Ann_weighted_148_restart_results_last)
remove(Sim_Ann_weighted_148_restart_results_last_last)
#sort models based on performance
aa <- sort(Sim_Ann_weighted_148_restart_results_last_last_conc$obj_value, decreasing = F, index.return=T)$ix
Sim_Ann_weighted_148_restart_results_last_last_conc$Enhancer_index <- Sim_Ann_weighted_148_restart_results_last_last_conc$Enhancer_index[aa, ]
Sim_Ann_weighted_148_restart_results_last_last_conc$parameters <- Sim_Ann_weighted_148_restart_results_last_last_conc$parameters[aa, ]
Sim_Ann_weighted_148_restart_results_last_last_conc$obj_value <- Sim_Ann_weighted_148_restart_results_last_last_conc$obj_value[aa]
rownames(Sim_Ann_weighted_148_restart_results_last_last_conc$Enhancer_index) <- aa
rownames(Sim_Ann_weighted_148_restart_results_last_last_conc$parameters) <- aa
names(Sim_Ann_weighted_148_restart_results_last_last_conc$obj_value) <- aa

#######
#setwd("SimAnn_Results_plot/Weighted_simAnn/")

names(Sim_Ann_weighted_148_restart_all_models_eval) <- rownames(Sim_Ann_weighted_148_restart_results_last_last_conc$parameters)

plotExpression(expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, .dendrogram = "both", .Rowv = T, .Colv = T,filename = "real_exp_2.png", .distfun = my_dist_fuc, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5))


Sim_Ann_weighted_148_restart_ExpMat <- lapply(Sim_Ann_weighted_148_restart_all_models_eval, "[[", 1)
Sim_Ann_weighted_148_restart_loss_accum <- lapply(Sim_Ann_weighted_148_restart_all_models_eval, "[[", 3)
Sim_Ann_weighted_148_restart_loss_perGene_perEnh <- lapply(Sim_Ann_weighted_148_restart_all_models_eval, "[[", 4)
Sim_Ann_weighted_148_restart_round_accuracy <- lapply(Sim_Ann_weighted_148_restart_all_models_eval, "[[", 5)
Sim_Ann_weighted_148_restart_parameters <- Sim_Ann_weighted_148_restart_results_last_last_conc$parameters
Sim_Ann_weighted_148_restart_enhancers <- lapply(Sim_Ann_weighted_148_restart_all_models_eval, "[[", 2)
Sim_Ann_weighted_148_restart_enhancers <- do.call(rbind, Sim_Ann_weighted_148_restart_enhancers)
Sim_Ann_weighted_148_restart_accuracy_perGene_perEnh <- lapply(Sim_Ann_weighted_148_restart_all_models_eval, "[[", 6)
Sim_Ann_weighted_148_restart_loss_perGene_perEnh_rounded <- lapply(Sim_Ann_weighted_148_restart_all_models_eval, "[[", 7)
Sim_Ann_weighted_148_restart_objVal <- Sim_Ann_weighted_148_restart_results_last_last_conc$obj_value

# compute accuracy using automatic thresholding
# needs re-evaluation since these predictions are already rounded
Sim_Ann_weighted_148_restart_all_models_eval_realValue <- list()
for(i in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  Sim_Ann_weighted_148_restart_all_models_eval_realValue[[i]] <- objFuncEval(my_optimized_par = Sim_Ann_weighted_148_restart_parameters[i, ],
                                                                             my_data = ER_52_opt_simAnn_2_input,
                                                                             input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                                                             round_to_int = F)
}


Sim_Ann_weighted_148_restart_round_accuracy_automatic <- numeric(length(Sim_Ann_weighted_148_restart_all_models_eval))
Sim_Ann_weighted_148_restart_round_automatic <- list()
for(i in 1:length(Sim_Ann_weighted_148_restart_round_accuracy_automatic)){
  Sim_Ann_weighted_148_restart_round_automatic[[i]] <- prediction_discretizer(prediction = Sim_Ann_weighted_148_restart_all_models_eval_realValue[[i]]$Predicted_expression_mat,
                                                                              label = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
}
for(i in 1:length(Sim_Ann_weighted_148_restart_round_accuracy_automatic)){
  Sim_Ann_weighted_148_restart_round_accuracy_automatic[i] <- sum(Sim_Ann_weighted_148_restart_round_automatic[[i]]$round_prediction == my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, na.rm = T)/sum(!is.na(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52))
}
hist(Sim_Ann_weighted_148_restart_round_accuracy_automatic)

# Check which enhancers have  been updated in the evaluation process
aa <- numeric(nrow(Sim_Ann_weighted_148_restart_enhancers))
aa2 <- matrix(nrow=nrow(Sim_Ann_weighted_148_restart_enhancers), ncol=length(Sim_Ann_weighted_148_restart_loss_perGene_perEnh[[1]]))
for(i in 1:nrow(Sim_Ann_weighted_148_restart_enhancers)){
  aa[i] <- sum(Sim_Ann_weighted_148_restart_enhancers[i,] == Sim_Ann_weighted_148_restart_results_last_last_conc$Enhancer_index[i,])
  for(j in 1:ncol(aa2)){
    aa2[i, j] <- Sim_Ann_weighted_148_restart_loss_perGene_perEnh[[i]][[j]][Sim_Ann_weighted_148_restart_enhancers[i,j]] - Sim_Ann_weighted_148_restart_loss_perGene_perEnh[[i]][[j]][Sim_Ann_weighted_148_restart_results_last_last_conc$Enhancer_index[i,j]]
  }
}
boxplot.matrix(t(aa2))
# This analysis showed that in many cases the best enhancer was not chosen and was refined in evaluations --> this is at 
#  least in part because the optimized  objective function is weighted while the one evaluated is not.

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### Create confusion matrix for each model
# number of -1, 0, 1 s that have been predicted correctly versus the number present in the dataset
# columns are actual classes
# rows are predicted classes
Sim_Ann_weighted_148_restart_confusion <- confusion_constructor(predictionMat_List=Sim_Ann_weighted_148_restart_ExpMat,
                                                                realMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
names(Sim_Ann_weighted_148_restart_confusion) <- names(Sim_Ann_weighted_148_restart_results_last_last_conc$obj_value)
Sim_Ann_weighted_148_restart_confusion[[1]]
which.max(Sim_Ann_weighted_148_restart_round_accuracy)
sort(unlist(Sim_Ann_weighted_148_restart_round_accuracy), decreasing = T, index.return=T)$ix
Sim_Ann_weighted_148_restart_confusion[[20]]
Sim_Ann_weighted_148_restart_confusion[[40]]
aa <- lapply(Sim_Ann_weighted_148_restart_confusion, as.numeric)
aaa <- do.call(rbind, aa)

aa_S <-sort(aaa[, 9], decreasing = T, index.return=T)$ix
Sim_Ann_weighted_148_restart_confusion[[aa_S[1]]]

##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### 
# Create random predictions to compare with sim ann predictions
Sim_Ann_weighted_148_restart_random_pred <- CreateRandomPrediction(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, num = 148)
par(mfrow = c(1,1), mar = c(6,6,6,6))
plot(c(unlist(Sim_Ann_weighted_148_restart_round_accuracy),
       Sim_Ann_weighted_148_restart_random_pred$Precision),
     ylim = c(0, 1),
     col = c(rep(3, 148), rep(2, 148)),
     ylab = "Accuracy",
     xlab = "model index",
     cex = 0.8,
     pch = 16)
legend(200, 1.4, legend=c("Sim_Ann_weighted", "shuffled"),
       col=c("green", "red"), pch=16, cex=1.1,
       bty = "n", x.intersp = 0.1,
       y.intersp = 0.2, pt.cex = 1.1, box.col = 1)
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# Draw heatmap of the ensemble

aa <- PerformanceHeattMap_General(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                  prediction_mat_list = Sim_Ann_weighted_148_restart_ExpMat,
                                  .Colv = T, .Rowv = T, .dendrogram = "both",exportplot = T,
                                  filename = "PerformanceHeatMap_weighted.png",
                                  .RowSideColors = character(0))

boxplot.matrix(aa, las = 2)
aaa <- colSums(aa)
which.max(aaa)
aaa <- apply(aa, 2, median)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# draw heatmap of the parameters
colnames(Sim_Ann_weighted_148_restart_parameters) <- colnames(Sim_Ann_148_restart_parameters)
aa_abs_param <- Sim_Ann_weighted_148_restart_parameters
# aa_added <- range(aa_abs_param)[1]
# log10(-aa_added)
# aa_abs_param <- aa_abs_param - range(aa_abs_param)[1] + 1
# aa_abs_param <- log10(aa_abs_param)
# range(aa_abs_param)


aaa_breaks = c(
  seq(range(aa_abs_param[1:80,])[1], -1, length=5),    
  seq(-0.99, 0, length=3),           
  seq(0.01, 1, length=3),
  seq(1.01, 10, length = 5),
  seq(10.01, 100, length = 5),
  seq(100.01, 500, length = 5)
  #seq(500.01, 2000, length = 20)
) 
# aa_norm_par2 <- aa_norm_par
# aa_norm_par2[aa_norm_par2 > 2000] <- 2000

aa_abs_param <- cbind(aa_abs_param, unlist(Sim_Ann_weighted_148_restart_round_accuracy) * 100)
colnames(aa_abs_param)[ncol(aa_abs_param)] <- "accuracy_percent"

png(filename = "parameters_and_accuracy_weighted.png",    # create PNG for the heat map        
    width = 8*300,        # 8 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # smaller font size
heatmap.2(x = aa_abs_param,
          Rowv = T, Colv = T,
          dendrogram = "both",
          #rowsep = aa_rowsep-1, sepwidth = c(4,4), sepcolor = "green",
          #RowSideColors = aa_rowsidecol,
          trace="none", na.rm = T, symbreaks=F , symkey=F, symm = F, breaks = aaa_breaks,col = colorRampPalette(c("green", "black", "red"))(n = 25))
dev.off()

aa <- lapply(Sim_Ann_148_restart_confusion, as.numeric)
aaa <- do.call(rbind, aa)

aa_S <-sort(aaa[, 1], decreasing = T, index.return=T)$ix
boxplot.matrix(aa_norm_par2[aa_S[1:20], ], las = 2, ylim = c(-50, 200))
aa_S <-sort(aaa[, 5], decreasing = T, index.return=T)$ix
boxplot.matrix(aa_norm_par2[aa_S[1:20], ], las = 2, ylim = c(-50, 200))
aa_S <-sort(aaa[, 9], decreasing = T, index.return=T)$ix
boxplot.matrix(aa_norm_par2[aa_S[1:20], ], las = 2, ylim = c(-50, 200))
##### Find enhancers per gene##### Find enhancers per gene##### Find enhancers per gene
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# Using voting, for models better than 0.55 collect votes on enhancers per gene, see which enhancers have the highest votes for each gene
aa_table_list <- list()
for(i in 1:ncol(Sim_Ann_weighted_148_restart_enhancers)){
  aa_table_list[[i]] <- sort(table(Sim_Ann_weighted_148_restart_enhancers[, i]), decreasing = T)
}
aaa <- unlist(lapply(ER_52_opt_simAnn_2_input$Affinity_scores, nrow))
Enhancer_vote_barplot(Sim_Ann_weighted_148_restart_enhancers[1:50, ], aaa/4, .ylim=c(-50, 50))
############## scoring each enhancer
aa_sort <- sort(unlist(Sim_Ann_weighted_148_restart_round_accuracy), decreasing = T, index.return=T)$ix
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_weighted_148_restart_loss_perGene_perEnh_rounded,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = aa_sort[1],
                          filename = "weightedModel_enhancer_score_barplot_top1_class_loss_calib_rounded.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_weighted_148_restart_loss_perGene_perEnh_rounded,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = aa_sort[1:10],
                          filename = "weightedModel_enhancer_score_barplot_top10_class_loss_calib_rounded.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)

aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = Sim_Ann_weighted_148_restart_loss_perGene_perEnh,
                          enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                          model_index = c(1:10),
                          filename = "weightedModel_enhancer_score_barplot_top10_loss_calib.png",
                          loss = T,
                          real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
##################
#####################################################################################################################
#######################################                                      ########################################
#######################################      in silico TF KD evaluation      ########################################
#######################################                                      ########################################
#####################################################################################################################
######## setting TFs concentrations to zero one by one
Sim_Ann_weighted_148_restart_TF_KD <- TF_KD_evaluator(my_optimized_par_mat = Sim_Ann_weighted_148_restart_parameters,
                                                      .my_data=ER_52_opt_simAnn_2_input,
                                                      .input_gene_expMat=my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                                      .round_to_int = T,
                                                      TF_index=c(1:19),
                                                      WT_evaluation=numeric(0))
par(mfrow=c(1,1), mar = rep(4,4))
aa <- Sim_Ann_weighted_148_restart_TF_KD$loss_change_per_model
aaa <- do.call(rbind, aa)
rowMeans(aaa)
boxplot.matrix(t(aaa), las = 2, main = "(sum_loss(WT) - sum_loss(KD))/nu_genes")
abline(h = seq(-6, 2, 0.5), col = 2, lty = 8)
aa <- Sim_Ann_weighted_148_restart_TF_KD$accu_change_per_model
aaa <- do.call(rbind, aa)
rowMeans(aaa)
boxplot.matrix(t(aaa), las = 2, main = "(sum_acc(WT) - sum_acc(KD))/nu_genes")
abline(h = seq(-0.1, 0.25, 0.05), col = 2, lty = 8)
aa <- Sim_Ann_weighted_148_restart_TF_KD$loss_change_per_gene
aaa <- do.call(rbind, aa)
rowMeans(aaa)
boxplot.matrix(t(aaa), las = 2, main = "(sum_loss(WT) - sum_loss(KD))/nu_models", outline = T)
abline(h = seq(-7, 2, 0.5), col = 2, lty = 8)

aa <- Sim_Ann_weighted_148_restart_TF_KD$accu_change_per_gene
aaa <- do.call(rbind, aa)
rowMeans(aaa)
boxplot.matrix(t(aaa), las = 2, main = "(sum_accu(WT) - sum_accu(KD))/nu_models")
abline(h = seq(-0.1, 0.4, 0.05), col = 2, lty = 8)

########################################################################################################################
########################################################################################################################
###################################################                                     ################################
###################################################    TF concentration modification    ################################
###################################################                                     ################################
########################################################################################################################
########################################################################################################################
# TF concentration modification experiments

# Set of initial enhancers: most voted in sim ann (75  high performing models)
aa_table_list <- list()
for(i in 1:ncol(Sim_Ann_148_restart_enhancers)){
  aa_table_list[[i]] <- sort(table(Sim_Ann_148_restart_enhancers[1:75, i]), decreasing = T)
}
aa_mvEnh_ind <- integer(length(aa_table_list))
for(i in 1:length(aa_mvEnh_ind)){
  aa_mvEnh_ind[i] <- as.integer(names(aa_table_list[[i]])[1])
}

Sim_Ann_148_restart_enhancers_mostVoted <- aa_mvEnh_ind

# Set of initial parameters: median of top 75 models
Sim_Ann_148_restart_parameters_median_top75 <- apply(X = Sim_Ann_148_restart_parameters[1:75, ], MARGIN = 2, FUN = median)

# data
ER_52_opt_simAnn_input

# modification vector
aa <- c(1, 10, 13, 14, 15, 17)
aaa <- matrix(nrow = 0, ncol = 19)
for(i in 1:length(aa)){
  aaaa <- integer(19)
  aaaa[aa[i]] <- 1
  aaa <- rbind(aaa, aaaa)
  aaaa[aa[i]] <- -1
  aaa <- rbind(aaa, aaaa)
}
Optim_Greedy_conc_modif_vec <- aaa
colnames(Optim_Greedy_conc_modif_vec) <- names(TF.motifs.Shrinked.t)
rownames(Optim_Greedy_conc_modif_vec) <- c(1:nrow(Optim_Greedy_conc_modif_vec))

aa <- getwd()
setwd("TF_conc_modif/")
aaa <- unlist(lapply(ER_52_opt_simAnn_input$Affinity_scores, nrow))

GreedyJobConstructor(script_name="ER_Optimizer_Greedy_Conc_modif.R",
                     nu_jobs = 12,
                     file_name="ER_Optimizer_Conc_modif_Greedy_hal_job",
                     nu_genes=52,
                     repeat_rounds=2,
                     fixed_genes = which(aaa > 20))

#create the submit file creator
cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_submit_script_Conc_modif_greedy", sep="")
for(i in 1:104){
  cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
        paste("ER_Optimizer_Conc_modif_Greedy_hal_job", i, sep = "_"),
        nu_jobs,
        paste("tmp_ER_Optimizer_Conc_modif_Greedy_hal_job", i, sep = "_"),
        ">",
        paste("ER_Optimizer_Conc_modif_Greedy_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_submit_script_Conc_modif_greedy",
      sep=" ",
      append = T)
}
#making them executable
cat(c("#!/bin/bash", "\n"),file="ER_Optim_hal_Conc_modif_greedy_submit_exec", sep="")
for(i in 1:104){
  cat(c("chmod +x",
        paste("ER_Optimizer_Conc_modif_Greedy_hal_job_", i, ".submit", sep = "")
        ,"\n"),
      file="ER_Optim_hal_Conc_modif_greedy_submit_exec",
      sep=" ",
      append = T)
}
#create the dagman file
DAGmanConstructor(submit_prefix = "ER_Optimizer_Conc_modif_Greedy_hal_job",
                  start_ind=1,
                  end_ind=104,
                  filename="my_DagMan_tf")
setwd(aa)

step_counter <- 1
nu_jobs <- 12
nu_genes <- 52
save.image(file = "Inputs_greedy_concModif.RData")
########################################################################################################################
# reading the TF concentration modif results
aaff <- list.files(path = "Hal_Greedy_conc_modif_res")
Optim_Greedy_conc_modif_Results <- list()
for(i in 1:length(aaff)){
  load(paste0("Hal_Greedy_conc_modif_res/", aaff[i]))
  Optim_Greedy_conc_modif_Results[[i]] <- GreedyRes_conc
  remove(GreedyRes_conc)
}
par(mfrow = c(4, 3), mar = c(1,4,1,1))
for(i in 1:12){
  plot(unlist(lapply(Optim_Greedy_conc_modif_Results[[i]], "[[",4)))
}

Optim_Greedy_conc_modif_Results_last <- list()
for(i in 1:length(Optim_Greedy_conc_modif_Results)){
  Optim_Greedy_conc_modif_Results_last[[i]] <- Optim_Greedy_conc_modif_Results[[i]][[length(Optim_Greedy_conc_modif_Results[[i]])]]
}
names(Optim_Greedy_conc_modif_Results_last) <- c("ARp", "ARn", "GRp", "GRn", "PGRp", "PGRn", "RARAp","RARAn","RARGp", "RARGn", "RXRAp","RXRAn")
Optim_Greedy_conc_modif_Results_last_param <- lapply(Optim_Greedy_conc_modif_Results_last, "[[", 3)
Optim_Greedy_conc_modif_Results_last_param <- do.call(rbind, Optim_Greedy_conc_modif_Results_last_param)
Optim_Greedy_conc_modif_Results_last_obj <- lapply(Optim_Greedy_conc_modif_Results_last, "[[", 4)
Optim_Greedy_conc_modif_Results_last_obj <- do.call(c, Optim_Greedy_conc_modif_Results_last_obj)

Optim_Greedy_conc_modif_Results_eval <- list()

for(i in 1:nrow(Optim_Greedy_conc_modif_Results_last_param)){
  print(i)
  aa_modif_vec <- Optim_Greedy_conc_modif_vec[i, ]
  aa <- ER_52_opt_simAnn_input
  # setting controls to zero if 1
  aa$gene_conditon_mat[, which(aa_modif_vec == 1)] <- 0
  # setting treatments to one if 1
  aa$gene_conditon_mat[, (which(aa_modif_vec == 1) + 19)] <- 1
  # setting controls to one if -1
  aa$gene_conditon_mat[, which(aa_modif_vec == -1)] <- 1
  # setting treatments to zero if -1
  aa$gene_conditon_mat[, (which(aa_modif_vec == -1) + 19)] <- 0
  
  Optim_Greedy_conc_modif_Results_eval[[i]] <- objFuncEval(my_optimized_par = Optim_Greedy_conc_modif_Results_last_param[i, ],
                                                                   my_data = aa,
                                                                   input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                                                   round_to_int = T)
  
  #plotExpression(expMat = Sim_Ann_weighted_148_restart_all_models_eval[[i]]$Predicted_expression_mat, .dendrogram = "none", .Rowv = F, .Colv = F, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),filename = paste0("sim_ann_",i, "_" ,rownames(Sim_Ann_weighted_148_restart_results_last_last_conc$parameters)[i], ".png"))
}
names(Optim_Greedy_conc_modif_Results_eval) <- names(Optim_Greedy_conc_modif_Results_last)

Optim_Greedy_conc_modif_Results_ExpMat <- lapply(Optim_Greedy_conc_modif_Results_eval, "[[", 1)
Optim_Greedy_conc_modif_Results_loss_accum <- lapply(Optim_Greedy_conc_modif_Results_eval, "[[", 3)
Optim_Greedy_conc_modif_Results_loss_perGene_perEnh <- lapply(Optim_Greedy_conc_modif_Results_eval, "[[", 4)
Optim_Greedy_conc_modif_Results_round_accuracy <- lapply(Optim_Greedy_conc_modif_Results_eval, "[[", 5)
Optim_Greedy_conc_modif_Results_parameters <- Optim_Greedy_conc_modif_Results_last_param
Optim_Greedy_conc_modif_Results_enhancers <- lapply(Optim_Greedy_conc_modif_Results_eval, "[[", 2)
Optim_Greedy_conc_modif_Results_enhancers <- do.call(rbind, Optim_Greedy_conc_modif_Results_enhancers)
Optim_Greedy_conc_modif_Results_accuracy_perGene_perEnh <- lapply(Optim_Greedy_conc_modif_Results_eval, "[[", 6)
Optim_Greedy_conc_modif_Results_loss_perGene_perEnh_rounded <- lapply(Optim_Greedy_conc_modif_Results_eval, "[[", 7)
Optim_Greedy_conc_modif_Results_objVal <- Optim_Greedy_conc_modif_Results_last_obj

Optim_Greedy_conc_modif_Results_confusion <- confusion_constructor(predictionMat_List=Optim_Greedy_conc_modif_Results_ExpMat,
                                                                realMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)

names(Optim_Greedy_conc_modif_Results_confusion) <- names(Optim_Greedy_conc_modif_Results_eval)

Optim_Greedy_conc_modif_Results_confusion$ARn
unlist(Optim_Greedy_conc_modif_Results_round_accuracy)

length(Optim_Greedy_conc_modif_Results[[1]][[1]])
colnames(Optim_Greedy_conc_modif_vec)
Optim_Greedy_conc_modif_vec[2,]
# Evaluating the initial point for the  optimizations, to compare with the obtained results.
# enhancer
Sim_Ann_148_restart_enhancers_mostVoted
# Set of initial parameters: median of top 75 models
Sim_Ann_148_restart_parameters_median_top75
# evaluate 
aa_eval2 <- objFuncEval(my_optimized_par = Sim_Ann_148_restart_parameters_median_top75,
                        my_data = ER_52_opt_simAnn_input,
                        input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                        round_to_int = T,
                        chosen_enh = Sim_Ann_148_restart_enhancers_mostVoted)
aa_eval2$Rounded_Error
########################################################################################################################
########################################################################################################################
# run 148 dagman each containing 12 jobs for TF concentraion modifications
# use the weighted models, and their enhancers as the starting points.

# Set of initial enhancers: enhancers best performing by each of the weighted models
Sim_Ann_weighted_148_restart_enhancers

# Set of initial parameters: parameters of the weighted models
Sim_Ann_weighted_148_restart_parameters
 
# data
ER_52_opt_simAnn_2_input

# modification vector
aa <- c(1, 10, 13, 14, 15, 17)
aaa <- matrix(nrow = 0, ncol = 19)
for(i in 1:length(aa)){
  aaaa <- integer(19)
  aaaa[aa[i]] <- 1
  aaa <- rbind(aaa, aaaa)
  aaaa[aa[i]] <- -1
  aaa <- rbind(aaa, aaaa)
}
Optim_Greedy_conc_modif_vec <- aaa
colnames(Optim_Greedy_conc_modif_vec) <- names(TF.motifs.Shrinked.t)
rownames(Optim_Greedy_conc_modif_vec) <- c(1:nrow(Optim_Greedy_conc_modif_vec))

step_counter <- 1
nu_jobs <- 12
nu_genes <- 52

aa_wd <- getwd()
setwd("TF_conc_modif/TF_conc_148")
aaa <- unlist(lapply(ER_52_opt_simAnn_2_input$Affinity_scores, nrow))

for(i in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  GreedyJobConstructor(script_name=paste0("ER_Optimizer_Greedy_Conc_modif_",i ,".R"),
                       nu_jobs = 12,
                       file_name=paste0("ER_Optimizer_Conc_modif_Greedy_hal_job_", i),
                       nu_genes=52,
                       repeat_rounds=1,
                       fixed_genes = which(aaa > 20))
}


#create the submit file creator

for(j in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  cat(c("#!/bin/bash", "\n"),file=paste0("ER_Optim_hal_submit_script_Conc_modif_greedy_", j), sep="")
  for(i in 1:52){
    cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
          paste("ER_Optimizer_Conc_modif_Greedy_hal_job", j, i, sep = "_"),
          nu_jobs,
          paste("tmp_ER_Optimizer_Conc_modif_Greedy_hal_job", j, i, sep = "_"),
          ">",
          paste("ER_Optimizer_Conc_modif_Greedy_hal_job_", j,"_", i, ".submit", sep = "")
          ,"\n"),
        file=paste0("ER_Optim_hal_submit_script_Conc_modif_greedy_", j),
        sep=" ",
        append = T)
  }
}

#putting submit commands in one file
cat(c("#!/bin/bash", "\n"),file="Run_create_submit_files", sep="")
for(j in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  cat(c("./",paste("ER_Optim_hal_submit_script_Conc_modif_greedy_",j, sep = "")
        ,"\n"),
      file="Run_create_submit_files",
      sep="",
      append = T)
}

#making them executable
for(i in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  cat(c("#!/bin/bash", "\n"),file=paste0("ER_Optim_hal_Conc_modif_greedy_submit_exec_",i), sep="")
  for(j in 1:52){
    cat(c("chmod +x",
          paste("ER_Optimizer_Conc_modif_Greedy_hal_job_", i,"_",j, ".submit", sep = "")
          ,"\n"),
        file=paste0("ER_Optim_hal_Conc_modif_greedy_submit_exec_",i),
        sep=" ",
        append = T)
  }
}

# put all commands in one file
cat(c("#!/bin/bash", "\n"),file="Run_make_executable_submit_files", sep="")
for(j in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  cat(c("./",paste("ER_Optim_hal_Conc_modif_greedy_submit_exec_",j, sep = "")
        ,"\n"),
      file="Run_make_executable_submit_files",
      sep="",
      append = T)
}


#create the dagman file
for(i in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  DAGmanConstructor(submit_prefix = paste0("ER_Optimizer_Conc_modif_Greedy_hal_job_", i),
                    start_ind=1,
                    end_ind=52,
                    filename=paste0("my_DagMan_tf_", i))
}


cat(c("#!/bin/bash", "\n"),file="my_DagMan_all_148", sep="")
for(i in 1:148){
  cat(c("condor_submit_dag  my_DagMan_tf_",
        i
        ,"\n"),
      file="my_DagMan_all_148",
      sep="",
      append = T)
}

setwd(aa_wd)

step_counter <- 1
nu_jobs <- 12
nu_genes <- 52
save.image(file = "Inputs_greedy_concModif.RData")
#### reading and rewriting the R scripts

aa <- list()
for(i in 1:148){
  aa[[i]] <- readLines("ER_Optimizer_Greedy_Conc_modif_1.R")
  aa[[i]][11] <- paste0('if(file.exists(paste0(\"Optim_Greedy_conc_modif\","',"_",i,'_" , args[1], \".RData\"))){')
  aa[[i]][12] <- paste0('\tload(paste0(\"Optim_Greedy_conc_modif\","',"_",i,'_" , args[1], \".RData\"))')
  aa[[i]][17] <- paste0("\tmy_initial_chosen_enh <- Sim_Ann_weighted_148_restart_enhancers[", i, " ,]")
  aa[[i]][18] <- paste0("\tmy_initial_par <- Sim_Ann_weighted_148_restart_parameters[", i, " ,]")
  aa[[i]][22] <- paste0("set.seed((as.integer(args[1]) + (step_counter - 1)*nu_jobs) + ",(i-1),"* nu_genes * nu_jobs)")
  aa[[i]][24] <- "GreedyRes_conc[[step_counter]] <- SilicoConcenOptim(..my_data = ER_52_opt_simAnn_2_input,"
  aa[[i]][31] <- paste0('save.image(paste0(\"Optim_Greedy_conc_modif\","',"_",i,'_" , args[1], \".RData\"))')
  writeLines(aa[[i]], paste0("Modified/ER_Optimizer_Greedy_Conc_modif_", i, ".R"))
}
########################################################################################################################
# Reading the concentration modification 148 results
# create a list of length 12 where each entry corresponds to one TF conc modification
# each entry of the above list is a list of 148 entries.
# each of the 148 entries will contain the results of the modification for the model.


aa <- list.files("Hal_Greedy_conc_modif_res/Ensemble_Res/", pattern = "Optim*")
aaa <- strsplit(aa, split="_")
aaa1 <- lapply(aaa, "[[", 5)
aaa1 <- unlist(aaa1)
aaa1 <- as.integer(aaa1)
aaa2 <- lapply(aaa, "[[", 6)
aaa2 <- unlist(aaa2)
aaa2 <- strsplit(aaa2, split="\\.")
aaa2 <- unlist(lapply(aaa2, "[[", 1))
aaa2 <- as.integer(aaa2)
aa_model_index <- list()
for(i in 1:length(unique(aaa2))){
  aa_cur_modif_ind <- which(aaa2 %in% i)
  aa_cur_model <- aaa1[aa_cur_modif_ind]
  #print(length(unique(aa_cur_model)))
  aa_cur_modif_ind_sort <- sort(aa_cur_model, decreasing = F, index.return=T)$ix
  aa_cur_modif_sort_ind <- aa_cur_modif_ind[aa_cur_modif_ind_sort]
  aa_model_index[[i]] <- aa_cur_modif_sort_ind
  #print(aaa1[aa_cur_modif_sort_ind])
}

Optim_Greedy_conc_modif_Ensemble_Results <- list()
for(i in 1:length(aa_model_index)){
  Optim_Greedy_conc_modif_Ensemble_Results[[i]] <- list()
  for(j in 1:length(aa_model_index[[i]])){
    load(paste0("Hal_Greedy_conc_modif_res/Ensemble_Res/",aa[aa_model_index[[i]][j]]))
    Optim_Greedy_conc_modif_Ensemble_Results[[i]][[j]] <- GreedyRes_conc[[52]]
    remove(GreedyRes_conc)
  }
}
#looking at performance over iterations
Optim_Greedy_conc_modif_Ensemble_Perf_iterationsList <- list()

for(i in 1:length(aa_model_index)){
  Optim_Greedy_conc_modif_Ensemble_Perf_iterationsList[[i]] <- matrix(nrow=148, ncol = 52)
  for(j in 1:length(aa_model_index[[i]])){
    load(paste0("Hal_Greedy_conc_modif_res/Ensemble_Res/",aa[aa_model_index[[i]][j]]))
    Optim_Greedy_conc_modif_Ensemble_Perf_iterationsList[[i]][j, ] <- unlist(lapply(GreedyRes_conc, "[[", 4))
    remove(GreedyRes_conc)
  }
}


par(mfrow = c(3, 4))
for(j in 1:12){
  plot(Optim_Greedy_conc_modif_Ensemble_Perf_iterationsList[[j]][1,], ylim=range(Optim_Greedy_conc_modif_Ensemble_Perf_iterationsList[[j]]), type="l")
  for(i in 2:148){
    lines(Optim_Greedy_conc_modif_Ensemble_Perf_iterationsList[[j]][i,], col=col_vector[i%%length(col_vector)])
  }
}
# compare to loss of the seed models

#evaluate each model and then compare losses with the evaluated seed model
Sim_Ann_weighted_148_restart_all_models_eval
Optim_Greedy_conc_modif_Ensemble_Results_eval <- list()
Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracyDif <- matrix(nrow = 148, ncol=12)
Optim_Greedy_conc_modif_Ensemble_Results_SumErrorDif <- matrix(nrow = 148, ncol=12)
for(i in 1:length(Optim_Greedy_conc_modif_Ensemble_Results)){
  Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]] <- list()
  aa_modif_vec <- Optim_Greedy_conc_modif_vec[i, ]
  aa <- ER_52_opt_simAnn_2_input
  aa$gene_conditon_mat[, which(aa_modif_vec == 1)] <- 0
  aa$gene_conditon_mat[, (which(aa_modif_vec == 1) + 19)] <- 1
  aa$gene_conditon_mat[, which(aa_modif_vec == -1)] <- 1
  aa$gene_conditon_mat[, (which(aa_modif_vec == -1) + 19)] <- 0
  for(j in 1:length(Optim_Greedy_conc_modif_Ensemble_Results[[i]])){
    print(paste("modif nu", i, "model nu", j))
    Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]][[j]] <- objFuncEval(my_optimized_par = Optim_Greedy_conc_modif_Ensemble_Results[[i]][[j]]$parameters,
                                                                           my_data = aa,
                                                                           input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                                                           round_to_int = T)
    Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracyDif[j, i] <- Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]][[j]]$Rounded_Error - Sim_Ann_weighted_148_restart_all_models_eval[[j]]$Rounded_Error
    Optim_Greedy_conc_modif_Ensemble_Results_SumErrorDif[j, i] <- sum(Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]][[j]]$losses) - sum(Sim_Ann_weighted_148_restart_all_models_eval[[j]]$losses)
  }
}
colnames(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracyDif) <-  c("ARp", "ARn", "GRp", "GRn", "PGRp", "PGRn", "RARAp","RARAn","RARGp", "RARGn", "RXRAp","RXRAn")
colnames(Optim_Greedy_conc_modif_Ensemble_Results_SumErrorDif) <- c("ARp", "ARn", "GRp", "GRn", "PGRp", "PGRn", "RARAp","RARAn","RARGp", "RARGn", "RXRAp","RXRAn")
boxplot.matrix(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracyDif, main="Acc(modif) - Acc(nonModif)", las=2 )
abline(h=seq(-0.15, 0.10, 0.01), col = 2, lty=3)
boxplot.matrix(Optim_Greedy_conc_modif_Ensemble_Results_SumErrorDif, main="sumLosses(modif) - sumLosses(nonModif)", las=2 )
abline(h=seq(-100, 200, 10), col = 2, lty=3)
ER_52_opt_simAnn_2_input
Optim_Greedy_conc_modif_vec

Optim_Greedy_conc_modif_Ensemble_Results_SumError <- matrix(nrow=148, ncol=12)
Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy <- matrix(nrow=148, ncol=12)
aa_roundError <- numeric(148)
aa_sumLoss <- numeric(148)
for(i in 1:12){
  for(j in 1:148){
    Optim_Greedy_conc_modif_Ensemble_Results_SumError[j, i] <- sum(Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]][[j]]$losses)
    Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[j, i] <- Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]][[j]]$Rounded_Error
    if(i ==1){
      aa_roundError[j] <- Sim_Ann_weighted_148_restart_all_models_eval[[j]]$Rounded_Error
      aa_sumLoss[j] <- sum(Sim_Ann_weighted_148_restart_all_models_eval[[j]]$losses)
    }
  }
}
par(mfrow = c(3, 4))
aa_labels <- colnames(Optim_Greedy_conc_modif_Ensemble_Results_SumErrorDif)
for(i in 1:12){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_SumError[, i], aa_sumLoss, xlab = aa_labels[i], xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError),ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError), ylab="WT_sum_Loss")
  abline(a = c(0, 1), col = 2)
}
for(i in 1:12){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, i], aa_roundError, xlab = aa_labels[i], xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy), ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy),ylab="WT_accuracy")
  abline(a = c(0, 1), col = 2)
  
}

########
# re-evaluating the concentration modifications with non-rounded values
Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value <- list()
for(i in 1:length(Optim_Greedy_conc_modif_Ensemble_Results)){
  Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value[[i]] <- list()
  aa_modif_vec <- Optim_Greedy_conc_modif_vec[i, ]
  aa <- ER_52_opt_simAnn_2_input
  aa$gene_conditon_mat[, which(aa_modif_vec == 1)] <- 0
  aa$gene_conditon_mat[, (which(aa_modif_vec == 1) + 19)] <- 1
  aa$gene_conditon_mat[, which(aa_modif_vec == -1)] <- 1
  aa$gene_conditon_mat[, (which(aa_modif_vec == -1) + 19)] <- 0
  for(j in 1:length(Optim_Greedy_conc_modif_Ensemble_Results[[i]])){
    print(paste("modif nu", i, "model nu", j))
    Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value[[i]][[j]] <- objFuncEval(my_optimized_par = Optim_Greedy_conc_modif_Ensemble_Results[[i]][[j]]$parameters,
                                                                           my_data = aa,
                                                                           input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                                                           round_to_int = F)
    # Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracyDif[j, i] <- Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]][[j]]$Rounded_Error - Sim_Ann_weighted_148_restart_all_models_eval[[j]]$Rounded_Error
    # Optim_Greedy_conc_modif_Ensemble_Results_SumErrorDif[j, i] <- sum(Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]][[j]]$losses) - sum(Sim_Ann_weighted_148_restart_all_models_eval[[j]]$losses)
  }
}
Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat <- list()
Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat_accuracy <- list()

for(i in 1:length(Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value)){
  Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat[[i]] <- list()
  Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat_accuracy[[i]] <- numeric(length(Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value[[i]]))
  for(j in 1:length(Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value[[i]])){
    Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat[[i]][[j]] <- prediction_discretizer(prediction = Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value[[i]][[j]]$Predicted_expression_mat,
                                                                                                    label = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
    Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat_accuracy[[i]][j] <- sum(Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat[[i]][[j]]$round_prediction == my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, na.rm = T)/sum(!is.na(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52))
  }
}

boxplot.matrix(do.call(cbind, Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat_accuracy))

########
### plot for only top 10 best performing models in training
aa_top10_ind <- sort(aa_sumLoss, decreasing = F, index.return=T)$ix[1:10]
par(mfrow = c(3, 4))
aa_labels <- colnames(Optim_Greedy_conc_modif_Ensemble_Results_SumErrorDif)
for(i in 1:12){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_SumError[aa_top10_ind, i], aa_sumLoss[aa_top10_ind], xlab = aa_labels[i], xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError[aa_top10_ind, ]),ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError[aa_top10_ind, ]), ylab="WT_sum_Loss")
  abline(a = c(0, 1), col = 2)
}
for(i in 1:12){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[aa_top10_ind, i], aa_roundError[aa_top10_ind], xlab = aa_labels[i], xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[aa_top10_ind, ]), ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[aa_top10_ind, ]),ylab="WT_accuracy")
  abline(a = c(0, 1), col = 2)
  
}
par(mfrow= c(3, 2))
for(i in 1:6){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_SumError[aa_top10_ind, 2*i - 1],
       Optim_Greedy_conc_modif_Ensemble_Results_SumError[aa_top10_ind, 2*i],
       xlab = aa_labels[2*i - 1],
       xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError[aa_top10_ind, ]),
       ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError[aa_top10_ind, ]),
       ylab=aa_labels[2*i])
  abline(a = c(0, 1), col = 2)
}
for(i in 1:6){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[aa_top10_ind, 2*i - 1],
       Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[aa_top10_ind, 2*i],
       xlab = aa_labels[2*i - 1],
       xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[aa_top10_ind, ]),
       ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[aa_top10_ind, ]),
       ylab=aa_labels[2*i])
  abline(a = c(0, 1), col = 2)
}
########################################################################################################################
########################################################################################################################
# compute wilcoxon rank test, one sided
aa_roundError <- numeric(148)
aa_sumLoss <- numeric(148)
for(j in 1:148){
      aa_roundError[j] <- Sim_Ann_weighted_148_restart_all_models_eval[[j]]$Rounded_Error
      aa_sumLoss[j] <- sum(Sim_Ann_weighted_148_restart_all_models_eval[[j]]$losses)
  }

aa_wilcox <- list()
for(i in 1:ncol(Optim_Greedy_conc_modif_Ensemble_Results_SumError)){
  aa_wilcox[[i]] <- wilcox.test(x=Optim_Greedy_conc_modif_Ensemble_Results_SumError[, i],
                                y=aa_sumLoss,
                                paired=F,
                               # mu = 1,
                                alternative="less"
                                #,exact=T
                               )
}
aa <- unlist(lapply(aa_wilcox, "[[", 3))
names(aa) <- colnames(Optim_Greedy_conc_modif_Ensemble_Results_SumErrorDif)
aa

for(i in 1:12){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_SumError[, i],
       aa_sumLoss, xlab = aa_labels[i],
       xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError),
       ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError),
       ylab="WT_sum_Loss",
       main=format(aa[i], nsmall = 3)
       )
  abline(a = c(0, 1), col = 2)
}

aa_wilcox2 <- list()
for(i in 1:ncol(Optim_Greedy_conc_modif_Ensemble_Results_SumError)){
  aa_wilcox2[[i]] <- wilcox.test(x=Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, i],
                                y=aa_roundError,
                                paired=F,
                                # mu = 1,
                                alternative="greater"
                                #,exact=T
  )
}
aa2 <- unlist(lapply(aa_wilcox2, "[[", 3))
names(aa2) <- colnames(Optim_Greedy_conc_modif_Ensemble_Results_SumErrorDif)
aa2

par(mfrow = c(3, 4))
aa_labels <- colnames(Optim_Greedy_conc_modif_Ensemble_Results_SumErrorDif)

for(i in 1:12){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, i],
       aa_roundError, xlab = aa_labels[i],
       xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy),
       ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy),
       ylab="WT_accuracy", 
       main=format(aa2[i], nsmall = 3))
  abline(a = c(0, 1), col = 2)
  
}

########################################################################################################################
########################################################################################################################
############ How well do the modified conc models do on the RNA seq test set?
aa <- match(rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52),
            rownames(my_CommonDifExpMat_RNAseq))
RNAseq_Test_1_geneExp <- my_CommonDifExpMat_RNAseq[aa, ]
aa <- integer(0)
for(i in 1:nrow(RNAseq_Test_1_geneExp)){
  if(all(is.na(RNAseq_Test_1_geneExp[i, ]))){
    aa <- c(aa, i)
  }
}
RNAseq_Test_1_geneExp <- RNAseq_Test_1_geneExp[-aa, ]
aa_removed <- aa

Sim_Ann_weighted_148_ConcModif_RNASEQ_test_1_chosen_enh <- list()
for(i in 1:length(Optim_Greedy_conc_modif_Ensemble_Results_eval)){
  Sim_Ann_weighted_148_ConcModif_RNASEQ_test_1_chosen_enh[[i]] <- list()
  aa_modif_vec <- Optim_Greedy_conc_modif_vec[i, ]
  aa <- RNAseq_Test_1_dataset
  aa$gene_conditon_mat[, which(aa_modif_vec == 1)] <- 0
  aa$gene_conditon_mat[, (which(aa_modif_vec == 1) + 19)] <- 1
  aa$gene_conditon_mat[, which(aa_modif_vec == -1)] <- 1
  aa$gene_conditon_mat[, (which(aa_modif_vec == -1) + 19)] <- 0
  aa_cur_enh <-  do.call(rbind, lapply(Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]], "[[", 2))
  for(j in 1:length(Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]])){
    print(paste("modif nu", i, "model nu", j))
    Sim_Ann_weighted_148_ConcModif_RNASEQ_test_1_chosen_enh[[i]][[j]] <- objFuncEval(my_optimized_par = Optim_Greedy_conc_modif_Ensemble_Results[[i]][[j]]$parameters,
                                                                                     my_data = aa,
                                                                                     input_gene_expMat = RNAseq_Test_1_geneExp,
                                                                                     round_to_int = T,
                                                                                     chosen_enh = aa_cur_enh[j, -aa_removed])
  }
}
Sim_Ann_weighted_148_ConcModif_RNASEQ_test_1_chosen_enh[[1]][[1]]$Rounded_Error
#aa_random_pred <- CreateRandomPrediction(real_exp_mat = RNAseq_Test_1_geneExp, num = 148)
#par(mfrow = c(1,1), mar = c(6,6,6,6))

Optim_Greedy_conc_modif_Ensemble_Results_SumError_RNASEQ_test_chosenEnh <- matrix(nrow=148, ncol=12)
Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy_RNASEQ_test_chosenEnh <- matrix(nrow=148, ncol=12)
aa_roundError <- numeric(148)
aa_sumLoss <- numeric(148)
for(i in 1:12){
  for(j in 1:148){
    Optim_Greedy_conc_modif_Ensemble_Results_SumError_RNASEQ_test_chosenEnh[j, i] <- sum(Sim_Ann_weighted_148_ConcModif_RNASEQ_test_1_chosen_enh[[i]][[j]]$losses)
    Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy_RNASEQ_test_chosenEnh[j, i] <- Sim_Ann_weighted_148_ConcModif_RNASEQ_test_1_chosen_enh[[i]][[j]]$Rounded_Error
    if(i ==1){
      aa_roundError[j] <- Sim_Ann_weighted_148_RNASEQ_test_1_chosen_enh[[j]]$Rounded_Error 
      aa_sumLoss[j] <- sum(Sim_Ann_weighted_148_RNASEQ_test_1_chosen_enh[[j]]$losses)
    }
  }
}
par(mfrow = c(3, 4))
aa_labels <- colnames(Optim_Greedy_conc_modif_Ensemble_Results_SumErrorDif)
for(i in 1:12){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_SumError_RNASEQ_test_chosenEnh[, i], aa_sumLoss, xlab = aa_labels[i], xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError_RNASEQ_test_chosenEnh),ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError_RNASEQ_test_chosenEnh), ylab="WT_sum_Loss_test")
  abline(a = c(0, 1), col = 2)
}
for(i in 1:12){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy_RNASEQ_test_chosenEnh[, i], aa_roundError, xlab = aa_labels[i], xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy_RNASEQ_test_chosenEnh), ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy_RNASEQ_test_chosenEnh),ylab="WT_accuracy_test")
  abline(a = c(0, 1), col = 2)
  
}
??optimx
########################################################################################################################
########################################################################################################################
# conducting knock down exps on evaluated models of concentration modification results
# for each modification exp just knock down the TF that its concentration has changed,fix the enhancer and evaluate the rsults
Optim_Greedy_conc_modif_Ensemble_Results_KD_SumError <- matrix(nrow = length(Optim_Greedy_conc_modif_Ensemble_Results_eval[[1]]),
                                                               ncol=length(Optim_Greedy_conc_modif_Ensemble_Results_eval))
Optim_Greedy_conc_modif_Ensemble_Results_KD_RoundAccuracy <- matrix(nrow = length(Optim_Greedy_conc_modif_Ensemble_Results_eval[[1]]),
                                                               ncol=length(Optim_Greedy_conc_modif_Ensemble_Results_eval))

Optim_Greedy_conc_modif_Ensemble_Results_KD_eval <- list()
for(i in 1:length(Optim_Greedy_conc_modif_Ensemble_Results_eval)){
  Optim_Greedy_conc_modif_Ensemble_Results_KD_eval[[i]] <- list()
  aa_modif_vec <- Optim_Greedy_conc_modif_vec[i, ]
  aa <- ER_52_opt_simAnn_2_input
  aa$gene_conditon_mat[, which(aa_modif_vec == 1)] <- 0
  aa$gene_conditon_mat[, (which(aa_modif_vec == 1) + 19)] <- 0
  aa$gene_conditon_mat[, which(aa_modif_vec == -1)] <- 0
  aa$gene_conditon_mat[, (which(aa_modif_vec == -1) + 19)] <- 0
  for(j in 1:length(Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]])){
    print(paste("modif nu", i, "model nu", j))
    Optim_Greedy_conc_modif_Ensemble_Results_KD_eval[[i]][[j]] <- objFuncEval(my_optimized_par = Optim_Greedy_conc_modif_Ensemble_Results[[i]][[j]]$parameters,
                                                                           my_data = aa,
                                                                           input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                                                           round_to_int = T,
                                                                           chosen_enh = Optim_Greedy_conc_modif_Ensemble_Results_eval[[i]][[j]]$Enhancer_index)
    Optim_Greedy_conc_modif_Ensemble_Results_KD_SumError[j, i] <- sum(Optim_Greedy_conc_modif_Ensemble_Results_KD_eval[[i]][[j]]$losses)
    Optim_Greedy_conc_modif_Ensemble_Results_KD_RoundAccuracy[j, i] <- Optim_Greedy_conc_modif_Ensemble_Results_KD_eval[[i]][[j]]$Rounded_Error
    
  }
}
par(mfrow=c(3, 4))
for(i in 1:12){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_KD_SumError[, i],
       Optim_Greedy_conc_modif_Ensemble_Results_SumError[, i],
       xlab =paste(aa_labels[i], "KD",sep = "_"),
       xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_KD_SumError),
       ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_KD_SumError),
       ylab=paste(aa_labels[i], "modif",sep = "_"))
  abline(a = c(0, 1), col = 2)
}

par(mfrow=c(3, 4))
for(i in 1:12){
  plot(Optim_Greedy_conc_modif_Ensemble_Results_KD_RoundAccuracy[, i],
       Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, i],
       xlab =paste(aa_labels[i], "KD",sep = "_"),
       xlim=c(0.34, 0.66),
       ylim=c(0.34, 0.66),
       ylab=paste(aa_labels[i], "modif",sep = "_"))
  abline(a = c(0, 1), col = 2)
}
########################################################################################################################
########################################################################################################################
# looking more closely at enhancers chosen by GR modif top models
# are the same enhancers chosen by the top models
# what is the composition of the enhacers chosen by the top models for each gene.


aa_GR_loss <- rowSums(do.call(rbind, lapply(Optim_Greedy_conc_modif_Ensemble_Results_eval[[4]], "[[", 3)))
aa_GR_loss_sort_ind <- sort(aa_GR_loss, decreasing = F, index.return=T)$ix
aa_GR_enh <- do.call(rbind, lapply(Optim_Greedy_conc_modif_Ensemble_Results_eval[[4]], "[[", 2))
aaa <- unlist(lapply(ER_52_opt_simAnn_2_input$Affinity_scores, nrow))
par(mfrow=c(1,1))
aa <- Enhancer_vote_barplot(aa_GR_enh[aa_GR_loss_sort_ind[1:10], ], aaa/15, .ylim = c(-10, 10))
# get the list of chosen enhacner affinities 
which(aa[1,] >= 3)
my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52_gene_gte3common_index <- which(aa[1,] >= 3)
plotExpression(expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[which(aa[1,] >= 3),], 
               .dendrogram = "none", 
               .Rowv = F, 
               .Colv = F,
               colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),
               filename = "real_shrinked_commonEnh_gte3.png")
sum(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[which(aa[1,] >= 3),] == 1, na.rm = T)
sum(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[which(aa[1,] >= 3),] == 0, na.rm = T)
sum(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[which(aa[1,] >= 3),] == -1, na.rm = T)
sum(is.na(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[which(aa[1,] >= 3),]), na.rm = T)


sum(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52 == 1, na.rm = T)
sum(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52 == 0, na.rm = T)
sum(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52 == -1, na.rm = T)
sum(is.na(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52))

my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3 <- my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[which(aa[1,] >= 3), ]
aa_random_pred <- CreateRandomPrediction(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3, num = 148)
plot(aa_random_pred$Precision)

# compute the accuracy results for this subset: my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52_gene_gte3common_index
# random results
Sim_Ann_weighted_148_restart_random_pred_39 <- CreateRandomPrediction(real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3, num = 148)

# weighted sim ann results
Sim_Ann_weighted_148_restart_round_accuracy_automatic_39 <- numeric(length(Sim_Ann_weighted_148_restart_round_accuracy_automatic))
for(i in 1:length(Sim_Ann_weighted_148_restart_round_accuracy_automatic)){
  Sim_Ann_weighted_148_restart_round_accuracy_automatic_39[i] <- sum(Sim_Ann_weighted_148_restart_round_automatic[[i]]$round_prediction[my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52_gene_gte3common_index, ] == my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3, na.rm = T)/sum(!is.na(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3))
}

# Concentration modification results
Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat_accuracy_39 <- list()

for(i in 1:length(Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value)){
  #Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat[[i]] <- list()
  Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat_accuracy_39[[i]] <- numeric(length(Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value[[i]]))
  for(j in 1:length(Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value[[i]])){
    # Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat[[i]][[j]] <- prediction_discretizer(prediction = Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value[[i]][[j]]$Predicted_expression_mat,
    #                                                                                                      label = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
    Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat_accuracy_39[[i]][j] <- sum(Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat[[i]][[j]]$round_prediction[my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52_gene_gte3common_index,] == my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3, na.rm = T)/sum(!is.na(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3))
  }
}

# double modif results
Optim_Greedy_conc_modif_Double_Results_automated_thresh_accuracy_39 <- numeric(length(Optim_Greedy_conc_modif_Double_Results_eval_real_value))
for(i in 1:length(Optim_Greedy_conc_modif_Double_Results_eval_real_value)){
  # Optim_Greedy_conc_modif_Double_Results_automated_thresh[[i]] <- prediction_discretizer(prediction = Optim_Greedy_conc_modif_Double_Results_eval_real_value[[i]]$Predicted_expression_mat,
  #                                                                                        label=my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
  Optim_Greedy_conc_modif_Double_Results_automated_thresh_accuracy_39[i] <- sum(Optim_Greedy_conc_modif_Double_Results_automated_thresh[[i]]$round_prediction[my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52_gene_gte3common_index, ] == my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3, na.rm = T)/sum(!is.na(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3))
}
########################################################################################################################
########################################################################################################################

# performing a double modification for GRn and RARp to see how does the model change
# run 148 dagman containing one job for the double TF concentraion modification
# use the weighted models, and their enhancers as the starting points.

# Set of initial enhancers: enhancers best performing by each of the weighted models
Sim_Ann_weighted_148_restart_enhancers

# Set of initial parameters: parameters of the weighted models
Sim_Ann_weighted_148_restart_parameters

# data
ER_52_opt_simAnn_2_input

# modification vector
aa <- c(1, 10, 13, 14, 15, 17)
aaa <- matrix(nrow = 0, ncol = 19)
for(i in 1:length(aa)){
  aaaa <- integer(19)
  aaaa[aa[i]] <- 1
  aaa <- rbind(aaa, aaaa)
  aaaa[aa[i]] <- -1
  aaa <- rbind(aaa, aaaa)
}
Optim_Greedy_conc_modif_vec_double <- matrix(aaa[4, ] + aaa[7, ], byrow = T, nrow=1)
colnames(Optim_Greedy_conc_modif_vec_double) <- names(TF.motifs.Shrinked.t)
rownames(Optim_Greedy_conc_modif_vec_double) <- c(1:nrow(Optim_Greedy_conc_modif_vec_double))

step_counter <- 1
nu_jobs <- 1
nu_genes <- 52

aa_wd <- getwd()
setwd("TF_conc_modif/Double_KD")
aaa <- unlist(lapply(ER_52_opt_simAnn_2_input$Affinity_scores, nrow))

for(i in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  GreedyJobConstructor(script_name=paste0("ER_Optimizer_Greedy_Conc_modif_",i ,".R"),
                       nu_jobs = 1,
                       file_name=paste0("ER_Optimizer_Conc_modif_double_Greedy_hal_job_", i),
                       nu_genes=52,
                       repeat_rounds=1,
                       fixed_genes = which(aaa > 20))
}


#create the submit file creator

for(j in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  cat(c("#!/bin/bash", "\n"),file=paste0("ER_Optim_hal_submit_script_Conc_modif_double_greedy_", j), sep="")
  for(i in 1:52){
    cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
          paste("ER_Optimizer_Conc_modif_double_Greedy_hal_job", j, i, sep = "_"),
          nu_jobs,
          paste("tmp_ER_Optimizer_Conc_modif_double_Greedy_hal_job", j, i, sep = "_"),
          ">",
          paste("ER_Optimizer_Conc_modif_double_Greedy_hal_job_", j,"_", i, ".submit", sep = "")
          ,"\n"),
        file=paste0("ER_Optim_hal_submit_script_Conc_modif_double_greedy_", j),
        sep=" ",
        append = T)
  }
}

#putting submit commands in one file
cat(c("#!/bin/bash", "\n"),file="Run_create_submit_files", sep="")
for(j in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  cat(c("./",paste("ER_Optim_hal_submit_script_Conc_modif_double_greedy_",j, sep = "")
        ,"\n"),
      file="Run_create_submit_files",
      sep="",
      append = T)
}

#making them executable
for(i in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  cat(c("#!/bin/bash", "\n"),file=paste0("ER_Optim_hal_Conc_modif_double_greedy_submit_exec_",i), sep="")
  for(j in 1:52){
    cat(c("chmod +x",
          paste("ER_Optimizer_Conc_modif_double_Greedy_hal_job_", i,"_",j, ".submit", sep = "")
          ,"\n"),
        file=paste0("ER_Optim_hal_Conc_modif_double_greedy_submit_exec_",i),
        sep=" ",
        append = T)
  }
}

# put all commands in one file
cat(c("#!/bin/bash", "\n"),file="Run_make_executable_submit_files", sep="")
for(j in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  cat(c("./",paste("ER_Optim_hal_Conc_modif_double_greedy_submit_exec_",j, sep = "")
        ,"\n"),
      file="Run_make_executable_submit_files",
      sep="",
      append = T)
}


#create the dagman file
for(i in 1:nrow(Sim_Ann_weighted_148_restart_parameters)){
  DAGmanConstructor(submit_prefix = paste0("ER_Optimizer_Conc_modif_double_Greedy_hal_job_", i),
                    start_ind=1,
                    end_ind=52,
                    filename=paste0("my_DagMan_tf_double_", i))
}


cat(c("#!/bin/bash", "\n"),file="my_DagMan_all_148", sep="")
for(i in 1:148){
  cat(c("condor_submit_dag  my_DagMan_tf_double_",
        i
        ,"\n"),
      file="my_DagMan_all_148",
      sep="",
      append = T)
}

setwd(aa_wd)

step_counter <- 1
nu_jobs <- 12
nu_genes <- 52
save.image(file = "Inputs_greedy_concModif.RData")
# use the already modified R scripts for the previous single modification runs
########################################################################################################################
########################################################################################################################
# read the double KD results

aa <- list.files("Hal_Greedy_conc_modif_res/Double_KD/", pattern = "Optim*")
aaa <- strsplit(aa, split="_")
aaa1 <- lapply(aaa, "[[", 5)
aaa1 <- unlist(aaa1)
aaa1 <- as.integer(aaa1)

aa_model_sorted_ind <- sort(aaa1, decreasing = F, index.return=T)$ix 



Optim_Greedy_conc_modif_DoubleKD_Results <- list()
for(i in 1:length(aa_model_sorted_ind)){
  load(paste0("Hal_Greedy_conc_modif_res/Double_KD/",aa[aa_model_sorted_ind[i]]))
  Optim_Greedy_conc_modif_DoubleKD_Results[[i]] <- GreedyRes_conc[[length(GreedyRes_conc)]]
  remove(GreedyRes_conc)
}

#looking at performance over iterations
Optim_Greedy_conc_modif_DoubleKD_Perf_iterationsList <- matrix(nrow=148, ncol = 52)
for(i in 1:length(aa_model_sorted_ind)){
    load(paste0("Hal_Greedy_conc_modif_res/Ensemble_Res/",aa[aa_model_sorted_ind[i]]))
  Optim_Greedy_conc_modif_DoubleKD_Perf_iterationsList[i, ] <- unlist(lapply(GreedyRes_conc, "[[", 4))
    remove(GreedyRes_conc)
}


par(mfrow = c(1, 1))
plot(Optim_Greedy_conc_modif_DoubleKD_Perf_iterationsList[1,], ylim=range(Optim_Greedy_conc_modif_DoubleKD_Perf_iterationsList), type="l")
for(i in 2:148){
    lines(Optim_Greedy_conc_modif_DoubleKD_Perf_iterationsList[i,], col=col_vector[i%%length(col_vector)])
}

#evaluate each model and then compare losses with the evaluated seed model
Sim_Ann_weighted_148_restart_all_models_eval
Optim_Greedy_conc_modif_Double_Results_eval <- list()
Optim_Greedy_conc_modif_Double_Results_RoundAccuracyDif <- numeric(148)
Optim_Greedy_conc_modif_Double_Results_SumErrorDif <- numeric(148)
Optim_Greedy_conc_modif_Double_Results_SumError <- numeric(148)
Optim_Greedy_conc_modif_Double_Results_RoundAccuracy <- numeric(148)

aa_modif_vec <- Optim_Greedy_conc_modif_vec_double
aa <- ER_52_opt_simAnn_2_input
aa$gene_conditon_mat[, which(aa_modif_vec == 1)] <- 0
aa$gene_conditon_mat[, (which(aa_modif_vec == 1) + 19)] <- 1
aa$gene_conditon_mat[, which(aa_modif_vec == -1)] <- 1
aa$gene_conditon_mat[, (which(aa_modif_vec == -1) + 19)] <- 0
for(i in 1:length(Optim_Greedy_conc_modif_DoubleKD_Results)){
    print(paste("model nu", i))
    Optim_Greedy_conc_modif_Double_Results_eval[[i]] <- objFuncEval(my_optimized_par = Optim_Greedy_conc_modif_DoubleKD_Results[[i]]$parameters,
                                                                           my_data = aa,
                                                                           input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                                                           round_to_int = T)
    Optim_Greedy_conc_modif_Double_Results_RoundAccuracyDif[i] <- Optim_Greedy_conc_modif_Double_Results_eval[[i]]$Rounded_Error - Sim_Ann_weighted_148_restart_all_models_eval[[i]]$Rounded_Error
    Optim_Greedy_conc_modif_Double_Results_SumErrorDif[i] <- sum(Optim_Greedy_conc_modif_Double_Results_eval[[i]]$losses) - sum(Sim_Ann_weighted_148_restart_all_models_eval[[i]]$losses)
    Optim_Greedy_conc_modif_Double_Results_SumError[i] <- sum(Optim_Greedy_conc_modif_Double_Results_eval[[i]]$losses)
    Optim_Greedy_conc_modif_Double_Results_RoundAccuracy[i] <- Optim_Greedy_conc_modif_Double_Results_eval[[i]]$Rounded_Error
}
#######
# re-evaluating to get real value predictions
Optim_Greedy_conc_modif_Double_Results_eval_real_value <- list()
aa_modif_vec <- Optim_Greedy_conc_modif_vec_double
aa <- ER_52_opt_simAnn_2_input
aa$gene_conditon_mat[, which(aa_modif_vec == 1)] <- 0
aa$gene_conditon_mat[, (which(aa_modif_vec == 1) + 19)] <- 1
aa$gene_conditon_mat[, which(aa_modif_vec == -1)] <- 1
aa$gene_conditon_mat[, (which(aa_modif_vec == -1) + 19)] <- 0
for(i in 1:length(Optim_Greedy_conc_modif_DoubleKD_Results)){
  print(paste("model nu", i))
  Optim_Greedy_conc_modif_Double_Results_eval_real_value[[i]] <- objFuncEval(my_optimized_par = Optim_Greedy_conc_modif_DoubleKD_Results[[i]]$parameters,
                                                                  my_data = aa,
                                                                  input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                                                                  round_to_int = F)
}

Optim_Greedy_conc_modif_Double_Results_automated_thresh <- list()
Optim_Greedy_conc_modif_Double_Results_automated_thresh_accuracy <- numeric(length(Optim_Greedy_conc_modif_Double_Results_eval_real_value))
for(i in 1:length(Optim_Greedy_conc_modif_Double_Results_eval_real_value)){
  Optim_Greedy_conc_modif_Double_Results_automated_thresh[[i]] <- prediction_discretizer(prediction = Optim_Greedy_conc_modif_Double_Results_eval_real_value[[i]]$Predicted_expression_mat,
                                                                                    label=my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)
  Optim_Greedy_conc_modif_Double_Results_automated_thresh_accuracy[i] <- sum(Optim_Greedy_conc_modif_Double_Results_automated_thresh[[i]]$round_prediction == my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, na.rm = T)/sum(!is.na(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52))
}
#######

aa_roundError <- numeric(148)
aa_sumLoss <- numeric(148)

for(j in 1:148){
      aa_roundError[j] <- Sim_Ann_weighted_148_restart_all_models_eval[[j]]$Rounded_Error
      aa_sumLoss[j] <- sum(Sim_Ann_weighted_148_restart_all_models_eval[[j]]$losses)
}

par(mfrow = c(1, 1))
plot(Optim_Greedy_conc_modif_Double_Results_SumError,
     aa_sumLoss,
     xlab = "RARAp_GRn",
     xlim=range(Optim_Greedy_conc_modif_Double_Results_SumError),
     ylim=range(Optim_Greedy_conc_modif_Double_Results_SumError),
     ylab="WT_sum_Loss")
abline(a = c(0, 1), col = 2)

Optim_Greedy_conc_modif_Double_Results_RoundAccuracy <- unlist(lapply(Optim_Greedy_conc_modif_Double_Results_eval,"[[", 5))
plot(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy,
     aa_roundError,
     xlab = "RARAp_GRn",
     xlim=range(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy),
     ylim=range(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy),
     ylab="WT_accuracy")
abline(a = c(0, 1), col = 2)

### plot for only top 10 best performing models in training
aa_top10_ind <- sort(aa_sumLoss, decreasing = F, index.return=T)$ix[1:10]
par(mfrow = c(1, 1))
plot(Optim_Greedy_conc_modif_Double_Results_SumError[aa_top10_ind],
     aa_sumLoss[aa_top10_ind],
     xlab = "RARAp_GRn",
     xlim=range(Optim_Greedy_conc_modif_Double_Results_SumError[aa_top10_ind]),
     ylim=range(Optim_Greedy_conc_modif_Double_Results_SumError[aa_top10_ind]),
     ylab="WT_sum_Loss")
abline(a = c(0, 1), col = 2)

plot(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy[aa_top10_ind],
     aa_roundError[aa_top10_ind],
     xlab = "RARAp_GRn",
     xlim=range(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy[aa_top10_ind]),
     ylim=range(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy[aa_top10_ind]),
     ylab="WT_accuracy")
abline(a = c(0, 1), col = 2)
  


aa_wilcox <-list()
aa_wilcox[[1]]<- wilcox.test(x=Optim_Greedy_conc_modif_Double_Results_SumError,
                        y=aa_sumLoss,
                        paired=F,
                        # mu = 1,
                        alternative="less"
                        #,exact=T
                        )
aa_wilcox[[2]]<- wilcox.test(x=Optim_Greedy_conc_modif_Double_Results_SumError,
                           y=Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 4],
                           paired=F,
                           # mu = 1,
                           alternative="less"
                           #,exact=T
)
aa_wilcox[[3]]<- wilcox.test(x=Optim_Greedy_conc_modif_Double_Results_SumError,
                           y=Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 7],
                           paired=F,
                           # mu = 1,
                           alternative="less"
                           #,exact=T
)

aa <- unlist(lapply(aa_wilcox, "[[", 3))

par(mfrow= c(2, 2))

plot(Optim_Greedy_conc_modif_Double_Results_SumError,
     aa_sumLoss,
      xlab = "RARAp_GRn_sum_Loss" ,
      xlim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError),
      ylim=range(Optim_Greedy_conc_modif_Ensemble_Results_SumError),
      ylab= "WT_sum_Loss",
     main=format(aa[1], nsmall = 3))
abline(a = c(0, 1), col = 2)
plot(Optim_Greedy_conc_modif_Double_Results_SumError,
     Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 4],
     xlab = "RARAp_GRn_sum_Loss" ,
     xlim=range(cbind(Optim_Greedy_conc_modif_Ensemble_Results_SumError,
                      Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 4])),
     ylim=range(cbind(Optim_Greedy_conc_modif_Ensemble_Results_SumError,
                      Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 4])),
     ylab= "GRn_sum_Loss", 
     main=format(aa[2], nsmall = 3))
abline(a = c(0, 1), col = 2)

plot(Optim_Greedy_conc_modif_Double_Results_SumError,
     Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 7],
     xlab = "RARAp_GRn_sum_Loss" ,
     xlim=range(cbind(Optim_Greedy_conc_modif_Ensemble_Results_SumError,
                      Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 7])),
     ylim=range(cbind(Optim_Greedy_conc_modif_Ensemble_Results_SumError,
                      Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 7])),
     ylab= "RARAp_sum_Loss",
     main=format(aa[3], nsmall = 3))
abline(a = c(0, 1), col = 2)

aa <- cbind(aa_sumLoss,
            Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 4],
            Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 7],
            Optim_Greedy_conc_modif_Double_Results_SumError)
colnames(aa) <- c("WT", "GR negative", "RARA positive", "Double modification")

par(mfrow = c(1,1), mar = c(9,4,4,4))
boxplot.matrix(aa,
               las=2, ylab = "Loss")

aa_wilcox2 <-list()
aa_wilcox2[[1]]<- wilcox.test(x=Optim_Greedy_conc_modif_Double_Results_RoundAccuracy,
                             y=aa_roundError,
                             paired=F,
                             # mu = 1,
                             alternative="greater"
                             #,exact=T
)
aa_wilcox2[[2]]<- wilcox.test(x=Optim_Greedy_conc_modif_Double_Results_RoundAccuracy,
                             y=Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, 4],
                             paired=F,
                             # mu = 1,
                             alternative="greater"
                             #,exact=T
)
aa_wilcox2[[3]]<- wilcox.test(x=Optim_Greedy_conc_modif_Double_Results_RoundAccuracy,
                             y=Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, 7],
                             paired=F,
                             # mu = 1,
                             alternative="greater"
                             #,exact=T
)

aa2 <- unlist(lapply(aa_wilcox2, "[[", 3))


par(mfrow= c(2, 2))
plot(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy,
     aa_roundError,
     xlab = "RARAp_GRn_accuracy" ,
       xlim=range(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy),
       ylim=range(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy),
       ylab="WT_accuracy", main=aa2[1])
abline(a = c(0, 1), col = 2)

plot(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy,
     Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, 4],
     xlab = "RARAp_GRn_accuracy" ,
     xlim=range(cbind(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy,
                      Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, 4])),
     ylim=range(cbind(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy,
                      Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, 4])),
     ylab= "GRn_accuracy", main=aa2[2])
abline(a = c(0, 1), col = 2)

plot(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy,
     Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, 7],
     xlab = "RARAp_GRn_accuracy" ,
     xlim=range(cbind(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy,
                      Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, 7])),
     ylim=range(cbind(Optim_Greedy_conc_modif_Double_Results_RoundAccuracy,
                      Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, 7])),
     ylab= "RARAp_accuracy", main=aa2[3])
abline(a = c(0, 1), col = 2)

########################################################################################################################
########################################################################################################################

# GR ChIP peak analysis
# in the concentration modification experiments we observed that setting GR concentration to 1 in control and
#  zero in E2 treatment conditions improved the training results. This observation can be supported by ChIP peak evidence.
#  If the GR ChIP peaks dissapear after E2 treatment on the set of enhancers chosen by the top 10 models. 
#  (or top ten with improvement after concentration modification)
# I don't have GR peaks after E2 treatment alone

GR_ChIP_list <- list()
GR_ChIP_list[[1]] <- ReMapChIP[ReMapChIP$V4 == "GSE72249.nr3c1.mcf7",]
GR_ChIP_list[[2]] <- ReMapChIP[ReMapChIP$V4 == "GSE72249.nr3c1.mcf7_dex",]
for(i in 1:length(GR_ChIP_list)){
  names(GR_ChIP_list[[i]]) <- c("Chr", "start", "end", "experiment","score","strand", "summit1","summit2","studyNumber" )
}
GR_ChIP_list_GRanges <- mapply(makeGRangesFromDataFrame, GR_ChIP_list)
names(GR_ChIP_list_GRanges) <- c("control", "dex")


aa_top10_ind_wt <- aa_top10_ind
aa_top10_ind_GR <-  sort(Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 4], decreasing = F, index.return=T)$ix[1:10] 
aa_top10_ind_RAR <- sort(Optim_Greedy_conc_modif_Ensemble_Results_SumError[, 7], decreasing = F, index.return=T)$ix[1:10] 

# getting the genes and enhancer indices that had better performance than their corresponding non modified model by (2.5 in loss)
aa_gene <- list()
aa_gene2 <- list()
aa <- matrix(nrow = 10, ncol = 52)
aa2 <- matrix(nrow = 10, ncol = 52)
aa_enh <- list()
aa_enh2 <- list()
for(i in 1:length(aa_top10_ind_GR)){
  aa[i, ] <- Optim_Greedy_conc_modif_Ensemble_Results_eval[[4]][[aa_top10_ind_GR[i]]]$losses - Sim_Ann_weighted_148_restart_all_models_eval[[aa_top10_ind_GR[i]]]$losses
  aa2[i, ] <- Optim_Greedy_conc_modif_Ensemble_Results_eval[[7]][[aa_top10_ind_RAR[i]]]$losses - Sim_Ann_weighted_148_restart_all_models_eval[[aa_top10_ind_RAR[i]]]$losses
  aa_gene[[i]] <- which(aa[i, ] < -5)
  aa_gene2[[i]] <-  which(aa2[i, ] < -5)
  aa_enh[[i]] <-  Optim_Greedy_conc_modif_Ensemble_Results_eval[[4]][[aa_top10_ind_GR[i]]]$Enhancer_index[aa_gene[[i]]]
  aa_enh2[[i]] <- Optim_Greedy_conc_modif_Ensemble_Results_eval[[7]][[aa_top10_ind_RAR[i]]]$Enhancer_index[aa_gene2[[i]]]
}
# create a Granges Object containing the coordinates of the above enhancers
aa_my_granges <- ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges
for(i in 1:length(aa_enh)){
  if(length(aa_enh[[i]]) > 0){
    for(j in 1:length(aa_enh[[i]])){
      if(i == 1 & j==1){
        aa_grange <- aa_my_granges[[aa_gene[[i]][j]]][aa_enh[[i]][j]]
      }else{
        aa_grange <- c(aa_grange, aa_my_granges[[aa_gene[[i]][j]]][aa_enh[[i]][j]])
      }
    }
  }

  if(length(aa_enh2[[i]]) > 0){
    for(k in 1:length(aa_enh2[[i]])){
      if(i == 1 & k==1){
        aa_grange2 <- aa_my_granges[[aa_gene2[[i]][k]]][aa_enh2[[i]][k]]
      }else{
        aa_grange2 <- c(aa_grange2, aa_my_granges[[aa_gene2[[i]][k]]][aa_enh2[[i]][k]])
      }
    }
  }
}
c(11, 12, 21, 32, 38, 47, 50)
c(6, 8, 18, 20, 35, 37)
aa_grange <- unique(aa_grange)
aa_grange2 <- unique(aa_grange2)

aa_overlap <- findOverlaps(aa_grange, GR_ChIP_list_GRanges$dex)
aa_grange[aa_overlap@from]

aa_overlap2 <- findOverlaps(aa_grange2, RARa.ERa.Binding.hg19.Granges)
aa_grange2[aa_overlap2@from]

colnames(aa) <- names(aa_my_granges)
colnames(aa2) <- names(aa_my_granges)

boxplot.matrix(aa, las=2, main="Loss_GRmod - Loss_wt top10")
boxplot.matrix(aa2, las=2, main="Loss_RARmod - Loss_wt top10")



########################################################################################################################
########################################################################################################################
###################################################                                     ################################
###################################################           RNA-seq test set          ################################
###################################################                                     ################################
########################################################################################################################
########################################################################################################################
# test the best models on the RNA-seq datasets, same genes different conditions
# create the dataset for test
aa <- match(rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), rownames(my_CommonDifExpMat_RNAseq))
RNAseq_Test_1_geneExp <- my_CommonDifExpMat_RNAseq[aa, ]
aa <- integer(0)
for(i in 1:nrow(RNAseq_Test_1_geneExp)){
  if(all(is.na(RNAseq_Test_1_geneExp[i, ]))){
    aa <- c(aa, i)
  }
}
RNAseq_Test_1_geneExp <- RNAseq_Test_1_geneExp[-aa, ]
plotExpression(expMat = RNAseq_Test_1_geneExp, .dendrogram = "none", .Rowv = F, .Colv = F, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),filename = "RNAseq_test_1.png")
RNAseq_Test_1_dataset <- createOptimDataset(GeneExpMat = RNAseq_Test_1_geneExp,
                                            TFExpMat = TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun_ER01,
                                            ChoppedEnhancerScores = ER_52_opt_simAnn_2_input$Affinity_scores[-aa])


# testing the weighted sim ANN MODELS
Sim_Ann_weighted_148_RNASEQ_test_1 <- list()
for(i in 1:nrow(Sim_Ann_weighted_148_restart_results_last_last_conc$parameters)){
  print(i)
  Sim_Ann_weighted_148_RNASEQ_test_1[[i]] <- objFuncEval(my_optimized_par = Sim_Ann_weighted_148_restart_results_last_last_conc$parameters[i, ],
                                                         my_data = RNAseq_Test_1_dataset,
                                                         input_gene_expMat = RNAseq_Test_1_geneExp,
                                                         round_to_int = T)
  
  #plotExpression(expMat = Sim_Ann_weighted_148_restart_all_models_eval[[i]]$Predicted_expression_mat, .dendrogram = "none", .Rowv = F, .Colv = F, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),filename = paste0("sim_ann_",i, "_" ,rownames(Sim_Ann_weighted_148_restart_results_last_last_conc$parameters)[i], ".png"))
}
Sim_Ann_weighted_148_RNASEQ_test_1[[1]]$accuracy_perGene_perEnh
table(RNAseq_Test_1_dataset$gene_conditon_mat[, 39])

aa_random_pred <- CreateRandomPrediction(real_exp_mat = RNAseq_Test_1_geneExp, num = 148)
par(mfrow = c(1,1), mar = c(6,6,6,6))
aa_acc <- lapply(Sim_Ann_weighted_148_RNASEQ_test_1,"[[", 5)
plot(c(unlist(aa_acc),
       aa_random_pred$Precision),
     ylim = c(0, 1),
     col = c(rep(3, 148), rep(2, 148)),
     ylab = "Accuracy",
     xlab = "model index",
     cex = 0.8,
     pch = 16)
legend(200, 1, legend=c("Sim_Ann", "Shuffled"),
       col=c("green", "red"), pch=16, cex=1.3,
       bty = "n", x.intersp = 0.1,
       y.intersp = 0.2, pt.cex = 1.3, box.col = 1)
aa_pred <- lapply(Sim_Ann_weighted_148_RNASEQ_test_1,"[[", 1)
aa_confusion <- confusion_constructor(predictionMat_List=aa_pred,
                                      realMat = RNAseq_Test_1_geneExp)
aa_num <- lapply(X = aa_confusion, as.numeric)
aa_num <- do.call(rbind, aa_num)

which.max(aa_num[, 15])
aa_confusion[[20]]
aa_confusion[[5]]
aa_sens <- aa_num[, c(5, 10, 15)]
aa_sens_sd <- apply(aa_sens, 1, sd)
aa_sens_sum <- rowSums(aa_sens)
sort(aa_sens_sum, decreasing = T, index.return=T)$ix
sort(aa_sens_sd, decreasing = F, index.return=T)$ix
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# Draw heatmap of the ensemble

aa <- PerformanceHeattMap_General(real_exp_mat = RNAseq_Test_1_geneExp,
                                  prediction_mat_list = lapply(Sim_Ann_weighted_148_RNASEQ_test_1,"[[", 1),
                                  .Colv = T, .Rowv = T, .dendrogram = "both",exportplot = T,
                                  filename = "PerformanceHeatMap_weighted_RNAseq_test_1.png",
                                  .RowSideColors = character(0))

boxplot.matrix(aa, las = 2)
aaa <- colSums(aa)
which.max(aaa)
aaa <- apply(aa, 2, median)
table(RNAseq_Test_1_geneExp)
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# testing the weighted sim ANN MODELS--- Forcing them to use the same enhancers they chose in training
aa <- match(rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52), rownames(my_CommonDifExpMat_RNAseq))
RNAseq_Test_1_geneExp <- my_CommonDifExpMat_RNAseq[aa, ]
aa <- integer(0)
for(i in 1:nrow(RNAseq_Test_1_geneExp)){
  if(all(is.na(RNAseq_Test_1_geneExp[i, ]))){
    aa <- c(aa, i)
  }
}
RNAseq_Test_1_geneExp <- RNAseq_Test_1_geneExp[-aa, ]

Sim_Ann_weighted_148_RNASEQ_test_1_chosen_enh <- list()
for(i in 1:nrow(Sim_Ann_weighted_148_restart_results_last_last_conc$parameters)){
  print(i)
  Sim_Ann_weighted_148_RNASEQ_test_1_chosen_enh[[i]] <- objFuncEval(my_optimized_par = Sim_Ann_weighted_148_restart_results_last_last_conc$parameters[i, ],
                                                                    my_data = RNAseq_Test_1_dataset,
                                                                    input_gene_expMat = RNAseq_Test_1_geneExp,
                                                                    round_to_int = T,
                                                                    chosen_enh = Sim_Ann_weighted_148_restart_enhancers[i, -aa])
  
  #plotExpression(expMat = Sim_Ann_weighted_148_restart_all_models_eval[[i]]$Predicted_expression_mat, .dendrogram = "none", .Rowv = F, .Colv = F, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),filename = paste0("sim_ann_",i, "_" ,rownames(Sim_Ann_weighted_148_restart_results_last_last_conc$parameters)[i], ".png"))
}
aa_random_pred <- CreateRandomPrediction(real_exp_mat = RNAseq_Test_1_geneExp, num = 148)
par(mfrow = c(1,1), mar = c(6,6,6,6))
aa_acc <- lapply(Sim_Ann_weighted_148_RNASEQ_test_1_chosen_enh,"[[", 5)
plot(c(unlist(aa_acc),
       aa_random_pred$Precision),
     ylim = c(0, 1),
     col = c(rep(3, 148), rep(2, 148)),
     ylab = "Accuracy",
     xlab = "model index",
     cex = 0.8,
     pch = 16)
legend(200, 1.4, legend=c("Sim_Ann", "Shuffled"),
       col=c("green", "red"), pch=16, cex=1.3,
       bty = "n", x.intersp = 0.1,
       y.intersp = 0.2, pt.cex = 1.3, box.col = 1)

aa_pred <- lapply(Sim_Ann_weighted_148_RNASEQ_test_1_chosen_enh,"[[", 1)
aa_confusion <- confusion_constructor(predictionMat_List=aa_pred,
                                      realMat = RNAseq_Test_1_geneExp)
aa_num <- lapply(X = aa_confusion, as.numeric)
aa_num <- do.call(rbind, aa_num)

which.max(aa_num[, 15])
aa_confusion[[5]]
aa_confusion[[20]]
aa_sens <- aa_num[, c(5, 10, 15)]
aa_sens_sd <- apply(aa_sens, 1, sd)
aa_sens_sum <- rowSums(aa_sens)
sort(aa_sens_sum, decreasing = T, index.return=T)$ix
sort(aa_sens_sd, decreasing = F, index.return=T)$ix
########################################################################################################################
###################################################        Testing on different set of genes        ######################################
########################################################################################################################
# can it predict for new set of genes?
# select a new set of genes, on which RNA-seq data is less sparse, build the chopped enhaners, affinity calculations, and evaluate the models
colnames(my_CommonDifExpMat_RNAseq)
aa_common <- my_CommonDifExpMat_RNAseq[ rownames(my_CommonDifExpMat_RNAseq) %in% Genes.Associated.REMAP.ER.Entrez, ]

aalak <- aa_common
aalak[is.na(aalak)] <- -100

aa <- apply(X = aalak, MARGIN = 1, FUN = (function(x) all(x == 0 | x == -100)))


#aa <- apply(X = my_CommonDifExpMat_16_ERassoc, MARGIN = 1, FUN = (function(x) all(is.na(x))))


aalak <- aalak[!aa, ]

aa4p <-  apply(X = aalak, MARGIN = 1, FUN = (function(x) sum(x == 1)))
aa4n <-  apply(X = aalak, MARGIN = 1, FUN = (function(x) sum(x == -1)))
aa2 <- apply(X = aalak, MARGIN = 1, FUN = (function(x) sum(x != 0 & x != -100)))
aa3 <- intersect(which(aa4p > 0 | aa4n > 0), which(aa2 >=4))

length(aa3[!rownames(aalak)[aa3] %in% rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)])

aa3 <- aa3[!rownames(aalak)[aa3] %in% rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52)]
RNAseq_Test_2_geneExp <- aalak[aa3, ]
RNAseq_Test_2_geneExp[RNAseq_Test_2_geneExp==-100] <- NA
plotExpression(expMat = RNAseq_Test_2_geneExp, .dendrogram = "none", .Rowv = F, .Colv = F, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),filename = "RNAseq_test_2.png")
######
RNAseq_Test_2_geneExp  ### is the test gene expression matrix for genes not included in training
######
###### Creating the chopped enhancers 
RNAseq_Test_2_enhancers <- EnhancerChopper_wrapper(Entrez_IDs = rownames(RNAseq_Test_2_geneExp),
                                                   .enhancer.per.gene=Vicinity100kb.Enhancer.plusPromoter.By.gene,
                                                   .Enhancer.Sequence=MCFEnhancers[,"sequ"],
                                                   .promoter.promoter.int=Promoter.gene.PromoterIntList,
                                                   .promoter.sequence=Promoter.Sequence.Char,
                                                   #.Enhancer_seq=character(0),
                                                   #.EnhancerGR=character(0),
                                                   #.motifScore_output,
                                                   .motifList=TF.motifs.Shrinked.t,
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

RNAseq_Test_2_dataset <- createOptimDataset(GeneExpMat = RNAseq_Test_2_geneExp,
                                            TFExpMat = TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun_ER01,
                                            ChoppedEnhancerScores = RNAseq_Test_2_enhancers$Chopped_Score)


# testing the weighted sim ANN MODELS
Sim_Ann_weighted_148_RNASEQ_test_2 <- list()
for(i in 1:nrow(Sim_Ann_weighted_148_restart_results_last_last_conc$parameters)){
  print(i)
  Sim_Ann_weighted_148_RNASEQ_test_2[[i]] <- objFuncEval(my_optimized_par = Sim_Ann_weighted_148_restart_results_last_last_conc$parameters[i, ],
                                                         my_data = RNAseq_Test_2_dataset,
                                                         input_gene_expMat = RNAseq_Test_2_geneExp,
                                                         round_to_int = T)
  
  #plotExpression(expMat = Sim_Ann_weighted_148_restart_all_models_eval[[i]]$Predicted_expression_mat, .dendrogram = "none", .Rowv = F, .Colv = F, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),filename = paste0("sim_ann_",i, "_" ,rownames(Sim_Ann_weighted_148_restart_results_last_last_conc$parameters)[i], ".png"))
}

aa <- unlist(lapply(Sim_Ann_weighted_148_RNASEQ_test_2, "[[", 5))
hist(aa)
aa_random_pred_2 <- CreateRandomPrediction(real_exp_mat = RNAseq_Test_2_geneExp, num = 148)
par(mfrow = c(1,1), mar = c(6,6,6,6))
aa_acc <- lapply(Sim_Ann_weighted_148_RNASEQ_test_2,"[[", 5)
plot(c(unlist(aa_acc),
       aa_random_pred_2$Precision),
     ylim = c(0, 1),
     col = c(rep(3, 148), rep(2, 148)),
     ylab = "Accuracy",
     xlab = "model index",
     cex = 0.8,
     pch = 16)
legend(200, 1.4, legend=c("Sim_Ann", "Shuffled"),
       col=c("green", "red"), pch=16, cex=1.1,
       bty = "n", x.intersp = 0.1,
       y.intersp = 0.2, pt.cex = 1.3, box.col = 1)
aa_pred_2 <- lapply(Sim_Ann_weighted_148_RNASEQ_test_2,"[[", 1)
aa_confusion <- confusion_constructor(predictionMat_List=aa_pred_2,
                                      realMat = RNAseq_Test_2_geneExp)
aa_num <- lapply(X = aa_confusion, as.numeric)
aa_num <- do.call(rbind, aa_num)

which.max(aa_num[, 1])
aa_confusion[[82]]
aa_confusion[[35]]
aa_sens <- aa_num[, c(5, 10, 15)]
aa_sens_sd <- apply(aa_sens, 1, sd)
aa_sens_sum <- rowSums(aa_sens)
sort(aa_sens_sum, decreasing = T, index.return=T)$ix
sort(aa_sens_sd, decreasing = F, index.return=T)$ix

aa_confusion[sort(aa_sens_sum, decreasing = T, index.return=T)$ix[1:10]]

## further analysis of test sets

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# Draw heatmap of the ensemble

aa <- PerformanceHeattMap_General(real_exp_mat = RNAseq_Test_2_geneExp,
                                  prediction_mat_list = lapply(Sim_Ann_weighted_148_RNASEQ_test_2, "[[", 1),
                                  .Colv = T, .Rowv = T, .dendrogram = "both",exportplot = T,
                                  filename = "PerformanceHeatMap_weighted_RNAseq_test_2.png",
                                  .RowSideColors = character(0))

boxplot.matrix(aa, las = 2)
aaa <- colSums(aa)
which.max(aaa)
aaa <- apply(aa, 2, median)
aa_ind <- sort(aaa, decreasing = T, index.return=T)$ix


aaa <- lapply(Sim_Ann_weighted_148_RNASEQ_test_2, "[[", 7)
aa <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = lapply(Sim_Ann_weighted_148_RNASEQ_test_2, "[[", 7),
                          enhancer_Granges = RNAseq_Test_2_enhancers$Chopped_GRanges,
                          model_index = aa_ind[1:30],
                          filename = "Enhancer_score_plot_weighted_RNAseq_test_2_top30_class.png",
                          loss = T,
                          real_exp_mat = RNAseq_Test_2_geneExp,
                          draw_plot = T)
######################################################################################################
# Checking if other optimization algs do better

aa_start_par <- runif(n = ((ncol(ER_52_opt_simAnn_2_input$gene_conditon_mat) - 1)/2 + 1), min = -1000, max = 1000)
aa_time_optim <- proc.time()
aa_optim <- optim(par = aa_start_par,
                  fn = obj_func_chosen_enh,
                  my_data = ER_52_opt_simAnn_2_input, 
                  chosen_enh_ind = Sim_Ann_weighted_148_restart_enhancers[1, ])
proc.time()
aa_time_sbplx <- proc.time()
aa_sbplx <- sbplx(x0 = aa_start_par,
                  fn = obj_func_chosen_enh,
                  my_data = ER_52_opt_simAnn_2_input,
                  chosen_enh_ind = Sim_Ann_weighted_148_restart_enhancers[1, ])
proc.time()

aa_optim$value
aa_sbplx$value


aa_test_input <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172[1:10,],
                                    TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01,
                                    ChoppedEnhancerScores = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Score[1:10]) 

aa_start_par <- runif(n = ((ncol(aa_test_input$gene_conditon_mat) - 1)/2  + 1), min = -10, max = 10)
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG"
              ,"xtol_rel"  = 1.0e-04,
              "maxeval"   = 100,
              "local_opts" = local_opts
              # ,      
              # check_derivatives = T, 
              # check_derivatives_tol = 1e-04,
              # check_derivatives_print='all'
              )

aa_test_enh <- integer(length(aa_test_input$Affinity_scores))
for(i in 1:length(aa_test_enh)){
  aa_test_enh[i] <- sample(x = c(1:unlist(lapply(aa_test_input$Affinity_scores, nrow))[i]), size = 1)
}

aa_time <- proc.time()
res2 <- nloptr( x0=aa_start_par,
                eval_f=obj_func_chosen_enh_gradient,
                my_data = aa_test_input,
                chosen_enh_ind = aa_test_enh,
                intercept_pergene=F, 
                constrained=F,
                return_grad = T, 
                #      lb=lb,
                #     ub=ub,
                #eval_g_ineq=obj_func_chosen_enh_ineq,
                #     eval_g_eq=eval_g_eq,
                opts=opts)

proc.time() - aa_time

aa_time2 <- proc.time()
aa_optim <- optim(par = aa_start_par,
                  fn = obj_func_chosen_enh,
                  my_data = aa_test_input, 
                  chosen_enh_ind = aa_test_enh)
proc.time() - aa_time2


proc.time()
print( res2 )
######################################################################################################
# test diffferent optimization algorithms in nloptr for the result and time consumed
aa <- nloptr.get.default.options()
aa_all_algs <- unlist(strsplit(aa$possible_values[1], split = ", "))
aaa <- strsplit(aa_all_algs, split = "_")
aa_ap <- which(unlist(lapply(aaa, "[[", 2)) == "LD" | unlist(lapply(aaa, "[[", 2)) == "GD")
aa_ap2 <- which(unlist(lapply(aaa, "[[", 2)) == "LN"| unlist(lapply(aaa, "[[", 2)) == "GN")
aa_time_list <- list()

aa_test_input <- createOptimDataset(GeneExpMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                                    TFExpMat=TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01,
                                    ChoppedEnhancerScores = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172_chopped_1000_250_E2F1$Filtered$Chopped_Score) 


aa_start_par <- runif(n = ((ncol(aa_test_input$gene_conditon_mat) - 1)/2  + 1), 
                      min = -10, max = 10)


aa_test_enh <- integer(length(aa_test_input$Affinity_scores))
for(i in 1:length(aa_test_enh)){
  aa_test_enh[i] <- sample(x = c(1:unlist(lapply(aa_test_input$Affinity_scores, nrow))[i]), size = 1)
}

aa_res <- list()
local_opts$algorithm <- "NLOPT_LD_MMA"
local_opts$xtol_rel <- 1.0e-04
for(i in 1:length(aa_ap)){
  print(aa_all_algs[aa_ap[i]])
  opts <- list( "algorithm" = aa_all_algs[aa_ap[i]]
                ,"xtol_rel"  = 1.0e-04,
                "maxeval"   = 100,
                "local_opts" = local_opts
                # ,      
                # check_derivatives = T, 
                # check_derivatives_tol = 1e-04,
                # check_derivatives_print='all'
  )
  
  aa_time <- proc.time()
  aa_res[[i]] <- nloptr( x0=aa_start_par,
                  eval_f=obj_func_chosen_enh_gradient,
                  my_data = aa_test_input,
                  chosen_enh_ind = aa_test_enh,
                  intercept_pergene=F, 
                  constrained=F,
                  return_grad = T, 
                  #      lb=lb,
                  #     ub=ub,
                  #eval_g_ineq=obj_func_chosen_enh_ineq,
                  #     eval_g_eq=eval_g_eq,
                  opts=opts)
  
  aa_time_list[[i]] <- proc.time() - aa_time
  print(aa_res[[i]])
  print(aa_time_list[[i]] )
}
names(aa_time_list) <- aa_all_algs[aa_ap]
names(aa_res) <-  aa_all_algs[aa_ap]

aa_res_nd <- list()
aa_time_nd <- list()
for(i in 1:length(aa_ap2)){
  print(aa_all_algs[aa_ap2[i]])
  opts <- list( "algorithm" = aa_all_algs[aa_ap2[i]]
                ,"xtol_rel"  = 1.0e-04,
                "maxeval"   = 100
                ,"local_opts" = local_opts
                # ,      
                # check_derivatives = T, 
                # check_derivatives_tol = 1e-04,
                # check_derivatives_print='all'
  )
  
  aa_time <- proc.time()
  aa_res_nd[[i]] <- nloptr( x0=aa_start_par,
                         eval_f=obj_func_chosen_enh,
                         my_data = aa_test_input,
                         chosen_enh_ind = aa_test_enh,
                         intercept_pergene=F, 
                         constrained=F,
                         return_grad = F, 
                         #      lb=lb,
                         #     ub=ub,
                         #eval_g_ineq=obj_func_chosen_enh_ineq,
                         #     eval_g_eq=eval_g_eq,
                         opts=opts)
  
  aa_time_nd[[i]] <- proc.time() - aa_time
  print(aa_res_nd[[i]])
  print(aa_time_nd[[i]] )
}

names(aa_time_nd) <- aa_all_algs[aa_ap2]
names(aa_res_nd) <-  aa_all_algs[aa_ap2]
aa_obj_d <- unlist(lapply(aa_res, "[[", 17))[is.finite(unlist(lapply(aa_res, "[[", 17)))]
aa_obj_n <- unlist(lapply(aa_res_nd, "[[", 17))[is.finite(unlist(lapply(aa_res_nd, "[[", 17)))]

aa_time_d <- unlist(lapply(aa_time_list, "[", 3))[is.finite(unlist(lapply(aa_res, "[[", 17)))]
aa_time_n <- unlist(lapply(aa_time_nd, "[", 3))[is.finite(unlist(lapply(aa_res_nd, "[[", 17)))]

aa_all_obj <- c(aa_obj_d, aa_obj_n)
aa_all_time <- c(aa_time_d, aa_time_n)
par(mar = c(14,4,1,4))
barplot(aa_all_obj/5, las = 2)
barplot(aa_all_time, las = 2)


# create barplot for comparing time and performance
aa_all_obj_time <- matrix(nrow = 2, ncol = length(aa_all_obj))
for(i in 1:ncol(aa_all_obj_time)){
  aa_all_obj_time[1, i] <- aa_all_time[i]
  aa_all_obj_time[2, i] <- aa_all_obj[i]/5 - aa_all_time[i]
}
colnames(aa_all_obj_time) <- names(aa_all_obj)
barplot(aa_all_obj_time, las = 2)

# evaluate the fast ones
aa_par_nd <- do.call(rbind, lapply(aa_res_nd, "[[", 18)[is.finite(unlist(lapply(aa_res_nd, "[[", 17)))]) 
aa_res_nd_eval <- list()
for(i in 1:nrow(aa_par_nd)){
  aa_res_nd_eval[[i]] <- objFuncEval(my_optimized_par = aa_par_nd[i, ],
                                     my_data = aa_test_input,
                                     input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                                     round_to_int = T,
                                     chosen_enh = aa_test_enh)
  
}
"NLOPT_LN_BOBYQA"
aa_par_d <- do.call(rbind, lapply(aa_res, "[[", 18)[is.finite(unlist(lapply(aa_res, "[[", 17)))]) 
aa_res_d_eval <- list()
for(i in 1:nrow(aa_par_d)){
  aa_res_d_eval[[i]] <- objFuncEval(my_optimized_par = aa_par_d[i, ],
                                     my_data = aa_test_input,
                                     input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                                     round_to_int = T, chosen_enh = aa_test_enh)
  
}
barplot(unlist(lapply(aa_res, "[[", 17))[is.finite(unlist(lapply(aa_res, "[[", 17)))], las = 2)
barplot(unlist(lapply(aa_res_nd, "[[", 17))[is.finite(unlist(lapply(aa_res_nd, "[[", 17)))], las = 2)


aa_time2 <- proc.time()
aa_optim <- optim(par = aa_start_par,
                  fn = obj_func_chosen_enh,
                  my_data = aa_test_input, 
                  chosen_enh_ind = aa_test_enh)
aa_opt_eval <- objFuncEval(my_optimized_par = aa_optim$par,
                           my_data = aa_test_input,
                           input_gene_expMat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                           round_to_int = T, chosen_enh = aa_test_enh)
proc.time() - aa_time2




aa_all_algs2 <-  c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
  "Brent")
aa_res2 <- list()
aa_time2_list <- list()

for(i in 1:length(aa_all_algs2[1:3])){
  print(aa_all_algs2[i])
  aa_time2 <- proc.time()
  aa_res2[[i]] <- optim(par = aa_start_par, method = aa_all_algs2[i],
                    fn = obj_func_chosen_enh, 
                    my_data = aa_test_input,
                    chosen_enh_ind = aa_test_enh,
                    gr = obj_func_chosen_enh_gradient_only,
                    intercept_pergene=F, 
                    constrained=F, 
                    return_grad = T)
  print(aa_res2[[i]])
  aa_time2_list[[i]] <- proc.time() - aa_time2
  print(aa_time2_list[[i]])
}

aa_time2 <- proc.time()
aa_r <- nlm(f = obj_func_chosen_enh_gradient,
            p = aa_start_par,
            print.level = 2, 
            my_data = aa_test_input,
            chosen_enh_ind = aa_test_enh,
            intercept_pergene=F, 
            constrained=F, 
            return_grad = T)
proc.time() - aa_time2

######
# Aparently "NLOPT_LN_BOBYQA" is faster and more accurate than other algs



