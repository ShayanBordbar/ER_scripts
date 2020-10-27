# GEMSTAT multi_enhacner weighted fold change input constructor and output parser
########################################################################################################################
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
########################################################################################################################
########################################                       #########################################################
########################################       LIBRARIES       #########################################################
########################################                       #########################################################
########################################################################################################################
library(sigmoid)
library(stringdist)
library(gplots)
########################################################################################################################
########################################                       #########################################################
########################################       Functions       #########################################################
########################################                       #########################################################
########################################################################################################################
# all functions in this script:
# 1. Enhancer_writer
# 2. RangeChecker
# 3. EnsembleParConstructor_Multi_enh
# 4. par_ff_lb_ub_coop_Writer_multi_enh
# 5. initialPar_Writer_multi_enh
# 6. freefixWriter_multi_enh
# 7. upper_lower_bound_Writer_multi_enh
# 8. coop_Writer_Multi_enh
# 9. Hal_job_writer_Multi_enh
# 10. GEMSTAT_init_BOlinear
# 11. GEMSTAT_output_Reader_multiEnh_ensemble
# 12. GEMSTAT_log_Reader_multiEnh
# 13. mean_RMSE_calculator
# 14. boxplot_grouped_by_label
# 15. prediction_discretizer
# 16. coop_parameter_creator_comb
# 17. KD_analyzer
# 18. bash_batch_run_zip_transfer_process
# 19. GenerateComb
# 20. count_site_from_annotation
# 21. count_site_from_annotation_combined
# 22. annotation_Writer
# 23. KD_Performer
# 24. read_KD_results
# 25. KD_Performer_by_params
# 26. GEMSTAT_input_copier_multienh
# 27. job_copier
# 28. round_custume_thresh
# 29. GEMSTAT_output_Reader_multiEnh
# 30. GeneFoldChange_by_TFFoldChange
# 31. set_row_col_nu
# 32. GeneFoldChange_by_TFFoldChange_model_RAP_wrapper
# 33. KD_by_affected_enhancers
# 34. KD_on_gene_evaluator
# 35. KD_on_gene_evaluator_multiTF_wrapper
# 36. annotation_overlap_finder
# 37. range_overlap
# 38. 
########################################################################################################################
Enhancer_writer <- function(.loss_per_model_per_gene_per_enh,
                            .enhancer_Granges,
                            .enhancer_sequence,
                            .model_index=c(1:length(.loss_per_model_per_gene_per_enh)),
                            .real_exp_mat,
                            max_per_gene=integer(0),
                            compare_to_expectation=T,
                            sd_from_expected=0,
                            seq_file_name="sequences.fa",
                            write_to_file=T){
  # gets the score of enhacners and chooses enhancers to feed to gemstat, chooses one among each 
  #  pair of overlapping enhancers
  # .loss_per_model_per_gene_per_enh : is a list where each entry is a list corresponding to a
  #  model. each model corresponding list contains
  #  one element per gene. That element contains a numeric vector indicating the loss value of
  #  each enhancer for that gene under the current model
  # .enhancer_Granges : is a list with one entry corresponding to each gene. that entry contains
  #  a GRanges object indicating th coordinates of each enhancer assigned to the gene
  # enhancer_sequence : is a list with one entry corresponding to each gene. that entry is a list
  #  where each of its entries is seauence of one of enhancers assigned to that gene. ordering must
  #  be the same as everything else
  # .model_index : index of models to be considred
  # .real_exp_mat : observed expression matrix
  # max_per_gene : maximum number of enhacners per gene
  # compare_to_expectation : if True, only uses enhancers that do better than expected score by
  #  shuffling
  # sd_from_expected : if compare_to_expectation is true, this option specifies how many
  #  standard deviations better that expected value shoud the error be, for the gene to be accepted.
  # name for the final sequence file
  # write_to_file : if True writes the sequences to a file
  # Uses Enhancer_Score_plot function to compute average scores for enhancers.
  enhancer_seq_list <- list()
  Chosen_GR_list <- list()
  enhancer_score_all <- Enhancer_Score_plot(accuracy_per_model_per_gene_per_enh = .loss_per_model_per_gene_per_enh,
                                            enhancer_Granges = .enhancer_Granges,
                                            model_index = .model_index,
                                            filename = "whatever.png",
                                            loss = T,
                                            real_exp_mat = .real_exp_mat,
                                            draw_plot = F)
  sorted_enhacner_index_list <- list() # for each gene holds the index of the enhancers sorted by performance
  for(cur_gene in 1:nrow(.real_exp_mat)){
    cur_sorted_ind <- sort(enhancer_score_all$score[[cur_gene]], decreasing = F, index.return=T)$ix
    if(compare_to_expectation){ # compared to expected value if specified 
      cur_sorted_filtered <- which(enhancer_score_all$score[[cur_gene]][cur_sorted_ind] <
                                     (enhancer_score_all$expected_score_shuffling[cur_gene] 
                                      - sd_from_expected * enhancer_score_all$expected_std_shuffling[cur_gene])) 
      cur_sorted_ind <- cur_sorted_ind[cur_sorted_filtered]
    }
    sorted_enhacner_index_list[[cur_gene]] <- cur_sorted_ind
    
  }# end of loop over genes
  names(sorted_enhacner_index_list) <- rownames(.real_exp_mat)
  sorted_enhacner_index_list_filtered <- sorted_enhacner_index_list
  # removing overlapping enhancers with lower scores
  for(cur_gene in 1:nrow(.real_exp_mat)){
    mark_for_removal = integer(length(enhancer_score_all$score[[cur_gene]])) #it will be removed if marked as 1
    if(length(sorted_enhacner_index_list[[cur_gene]]) > 0){
      cur_gene_overlap <- findOverlaps(.enhancer_Granges[[cur_gene]])
      for(cur_sort_enh in sorted_enhacner_index_list[[cur_gene]]){
        if(mark_for_removal[cur_sort_enh] == 0){
          cur_enh_overlap <- cur_gene_overlap@to[cur_gene_overlap@from == cur_sort_enh]
          mark_for_removal[cur_enh_overlap] <- 1
          mark_for_removal[cur_sort_enh] <- 0
        }
      } # end of loop over sorted enhancers 
      remaining_enh <- which(mark_for_removal == 0)
      sorted_remaining_enh <- sorted_enhacner_index_list_filtered[[cur_gene]] %in% remaining_enh
      sorted_enhacner_index_list_filtered[[cur_gene]] <- sorted_enhacner_index_list_filtered[[cur_gene]][sorted_remaining_enh]
      
      if(length(max_per_gene) > 0){# remove enhancers if there is cap on max number
        if(length(sorted_enhacner_index_list_filtered[[cur_gene]]) > max_per_gene){
          sorted_enhacner_index_list_filtered[[cur_gene]] <- sorted_enhacner_index_list_filtered[[cur_gene]][1:max_per_gene]
        }
      }
      # get the sequences in order
      enhancer_seq_list[[cur_gene]] <- unlist(.enhancer_sequence[[cur_gene]][sorted_enhacner_index_list_filtered[[cur_gene]]])
      Chosen_GR_list[[cur_gene]] <- .enhancer_Granges[[cur_gene]][sorted_enhacner_index_list_filtered[[cur_gene]]]
    }#end of if there is any enhancers left
  }#end of loop over genes
  names(enhancer_seq_list) <- names(sorted_enhacner_index_list_filtered)
  names(Chosen_GR_list) <- names(sorted_enhacner_index_list_filtered)
  # remove genes without good enhancers and write the enhacners
  nu_enh_per_gene <- unlist(lapply(enhancer_seq_list, length))
  removed_genes_index <- which(nu_enh_per_gene == 0)
  names(removed_genes_index) <- rownames(.real_exp_mat)[removed_genes_index]
  enhancer_seq_list <- enhancer_seq_list[nu_enh_per_gene > 0]
  if(write_to_file){
    WriteFastaOfBag(write_from_list = enhancer_seq_list,
                    output.File.Name = seq_file_name)
  }
  return(list(sequenceList = enhancer_seq_list,
              GRangesList = Chosen_GR_list,
              sorted_index_filtered_list=sorted_enhacner_index_list_filtered,
              enhancer_score_list=enhancer_score_all,
              removed_genes=removed_genes_index))
}
########################################################################################################################
########################################################################################################################
#example

aa <- Enhancer_writer(.loss_per_model_per_gene_per_enh = Sim_Ann_weighted_148_restart_loss_perGene_perEnh,
                      .enhancer_Granges = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges,
                      .enhancer_sequence = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq,
                      .model_index=c(1:50),
                      .real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52,
                      max_per_gene=integer(0),
                      compare_to_expectation=T,
                      sd_from_expected=1,
                      seq_file_name="sequences")
aa$sorted_index_filtered_list
########################################################################################################################
########################################################################################################################
RangeChecker <- function(vector_list, range_list, equal_allowed=T, return_detailed=F){
  # vector_list is a list where each entry is a vector of numeric values
  # range_list is a list where each entry is a range (two column matrix), corresponding to the same
  #  entry of vector_list
  # equal_allowed : if True it allows the entries to be equal to upper or lower bounds, otherwise
  #  it doesn't
  # return_detailed : if True returns all comparisions in a list
  # function checks if all vectors are in the range specified by the range_list, returns True if
  #  currect and False otherwise
  stopifnot(length(vector_list)==length(range_list),
            all(unlist(lapply(vector_list, length)) == unlist(lapply(range_list, nrow))))
  
  check_vec <- logical(length(vector_list))
  check_vec_list <- list()
  for(cur_ent in 1:length(vector_list)){
    check_vec_list[[cur_ent]] <- list()
    if(length(vector_list[[cur_ent]]) == 0){
      if(length(range_list[[cur_ent]]) == 0){
        check_vec[cur_ent] <- T
      }
    }else{
      if(equal_allowed){
        check_vec_list[[cur_ent]][[1]] <- (range_list[[cur_ent]][, 1] - vector_list[[cur_ent]]) <= 0
        check_vec_list[[cur_ent]][[2]] <- (range_list[[cur_ent]][, 2] - vector_list[[cur_ent]]) >= 0
        a <- all((range_list[[cur_ent]][, 1] - vector_list[[cur_ent]]) <= 0 )
        b <- all((range_list[[cur_ent]][, 2] - vector_list[[cur_ent]]) >= 0 )
      }else{
        check_vec_list[[cur_ent]][[1]] <- (range_list[[cur_ent]][, 1] - vector_list[[cur_ent]]) < 0
        check_vec_list[[cur_ent]][[2]] <- (range_list[[cur_ent]][, 2] - vector_list[[cur_ent]]) > 0
        a <- all((range_list[[cur_ent]][, 1] - vector_list[[cur_ent]]) < 0 )
        b <- all((range_list[[cur_ent]][, 2] - vector_list[[cur_ent]]) > 0 )
        
      }
      names(check_vec_list[[cur_ent]]) <- c("lower_bound", "upper_bound")
      check_vec[cur_ent] <- a&b
    }
  }
  if(return_detailed & !all(check_vec)){
    return(check_vec_list)
  }
  return(all(check_vec))
}
########################################################################################################################
########################################################################################################################
#example
aa <- list(c(1,2,3), c(6,0,1), numeric(0))
aa_r <- list(cbind(c(0,2,3), c(2,3,6)), cbind(c(0,0,0), c(6,6,6)), matrix(nrow = 0, ncol = 2))
RangeChecker(aa, aa_r)
########################################################################################################################
########################################################################################################################
EnsembleParConstructor_Multi_enh <- function(TFnames,
                                             annotation_rangeMat,
                                             annotation_logscale=F,
                                             binding_rangeMat,
                                             binding_logscale=F,
                                             activation_rangeMat,
                                             activation_logscale=F,
                                             coop_rangeMat,
                                             coop_logscale=F,
                                             qBTM_rangeMat,
                                             qBTM_logscale=F,
                                             log_Reg_bias_range,
                                             log_Reg_coeff_range,
                                             nu.sample,
                                             annotation_thresh_ff,
                                             initial_bind_w_ff,
                                             initial_alpha_ff,
                                             initial_coop_weight_ff,
                                             initial_qBTM_ff,
                                             log_Reg_bias_ff,
                                             log_Reg_coeff_ff,
                                             annotation_thresh,
                                             initial_bind_w,
                                             initial_alpha,
                                             initial_coop_weight,
                                             initial_qBTM,
                                             initial_logReg_bias,
                                             initial_logReg_coeff)
{
  # Create an ensemble of initial Par files for GEMSTAT
  # inputs to this function:
  # 1. TFnames : a character vector containing the name of TFs
  # 4. nu.sample : number of samples per compartment: we will have 2^number of parameters compartments,
  #  and "nu.sample" indicates how many samples shuld be retrieved from a compartment
  # 5. .coopertingTFIndex is a matrix where nrow = number of interacting pairs, ncol= 2 . Each row has
  #  the index of interacting TFs in one interaction (index in TF expression matrix)
  # log_Reg_bias_range : is the range of bias parameter of logisitc regression
  # log_Reg_coeff_range : is the range of coefficient parameter of logisitc regression
  number2binary = function(number, noBits) {
    # function creates the binary version of a number
    binary_vector = rev(as.numeric(intToBits(number)))
    if(missing(noBits)) {
      return(binary_vector)
    } else {
      binary_vector[-(1:(length(binary_vector) - noBits))]
    }
  }
  if(length(log_Reg_bias_range) > 0){
    log_Reg_bias_rangeMat <- matrix(log_Reg_bias_range, nrow = 1)
  }else{
    log_Reg_bias_rangeMat <- matrix(nrow = 0, ncol = 2)
  }
  if(length(log_Reg_coeff_range) > 0){
    log_Reg_coeff_rangeMat <- matrix(log_Reg_coeff_range, nrow = 1)
  }else{
    log_Reg_coeff_rangeMat <- matrix(nrow = 0, ncol = 2)
  }
  
  print("debug1")
  my_range_list <- list(annotation_rangeMat,
                        binding_rangeMat,
                        activation_rangeMat,
                        coop_rangeMat,
                        qBTM_rangeMat, 
                        log_Reg_bias_rangeMat,
                        log_Reg_coeff_rangeMat)
  all_free_fix <- c(annotation_thresh_ff,
                    initial_bind_w_ff,
                    initial_alpha_ff,
                    initial_coop_weight_ff,
                    initial_qBTM_ff, 
                    log_Reg_bias_ff,
                    log_Reg_coeff_ff)
  print("annotation_thresh_ff")
  print(annotation_thresh_ff)
  print("initial_bind_w_ff")
  print(initial_bind_w_ff)
  print("initial_alpha_ff")
  print(initial_alpha_ff)
  print("initial_coop_weight_ff")
  print(initial_coop_weight_ff)
  print("initial_qBTM_ff")
  print(initial_qBTM_ff)
  print("log_Reg_bias_ff")
  print(log_Reg_bias_ff)
  print("log_Reg_coeff_ff")
  print(log_Reg_coeff_ff)
  
  all_par.RangeMat <- rbind(annotation_rangeMat,
                            binding_rangeMat,
                            activation_rangeMat,
                            coop_rangeMat,
                            qBTM_rangeMat,
                            log_Reg_bias_rangeMat,
                            log_Reg_coeff_rangeMat)
  all_par_log_column <- c(rep(as.integer(annotation_logscale), nrow(annotation_rangeMat)), 
                          rep(as.integer(binding_logscale), nrow(binding_rangeMat)),
                          rep(as.integer(activation_logscale), nrow(activation_rangeMat)),
                          rep(as.integer(coop_logscale), nrow(coop_rangeMat)),
                          rep(as.integer(qBTM_logscale), nrow(qBTM_rangeMat)), 0, 0)
  all_par.RangeMat <- cbind(all_par.RangeMat, all_par_log_column)
  rownames(all_par.RangeMat) <- c(paste("annotation", c(1:nrow(annotation_rangeMat)), sep = "_"),
                              paste("binding", c(1:nrow(binding_rangeMat)), sep = "_"),
                              paste("activation", c(1:nrow(activation_rangeMat)), sep = "_"),
                              paste("coop", c(1:nrow(coop_rangeMat)), sep = "_"),
                              paste("qBTM", c(1:nrow(qBTM_rangeMat)), sep = "_"),
                              "logBias_1", "logCoeff_1"
                              )
  # print("rownames")
  # print(rownames(all_par.RangeMat))
  #print(all_free_fix)
  par.RangeMat <- all_par.RangeMat[all_free_fix==1, ]
  #print("rownames")
  #print(rownames(par.RangeMat))
  par.RangeMat
  Par.nu= nrow(par.RangeMat);
  Range_combination <- matrix(0, nrow=2^Par.nu, ncol=Par.nu);
  print(paste("Ensemble size is:", nrow(Range_combination) * nu.sample))
  for(i in 0:(2^Par.nu-1))
  {
    Range_combination[i+1,] = number2binary(i,Par.nu)
  }
  
  compartmentRanges <- list()
  for(i in 1:nrow(Range_combination)){
    rangeFile <- matrix(0, nrow=Par.nu, ncol=3);
    rangeFile[,3] <- par.RangeMat[, 3];
    for(j in 1:Par.nu){
      if(Range_combination[i, j]==0){
        rangeFile[j, 1] <- par.RangeMat[j, 1]
        rangeFile[j, 2] <- (par.RangeMat[j, 2]+par.RangeMat[j, 1])/2
      }
      else if(Range_combination[i,j]==1){
        rangeFile[j, 1] <- (par.RangeMat[j, 2]+par.RangeMat[j, 1])/2
        rangeFile[j, 2] <- par.RangeMat[j, 2]
      }
      else{
        rangeFile[j, 1] <- par.RangeMat[j, 1]
        rangeFile[j, 2] <- par.RangeMat[j, 2]
      }
    }
    compartmentRanges[[i]] <- rangeFile
  }
  
  
  #now sample from the ranges for each compartment
  compartmentSamples <- list()
  for(comp in 1:length(compartmentRanges)){
    ranges <- compartmentRanges[[comp]]
    Nsample <- nu.sample
    #n <- nrow(ranges)
    sample <- matrix(0, nrow=Nsample, ncol=Par.nu)
    for(i in 1:Par.nu){
      if(ranges[i, 3]==0){
        sample[, i] <- runif(Nsample, ranges[i, 1], ranges[i, 2])
      }else{
        sample[, i] <- 10^(runif(Nsample, ranges[i, 1], ranges[i, 2]))
      }
    }
    colnames(sample) <- rownames(par.RangeMat)
    compartmentSamples[[comp]] <- sample
  }
  # create a list of lists to hold the parameters in their initial format and
  #  transform the sampled parameter to their initial form using their assigned names
  #creating a template parameter list
  template_par_set <- list(annotation=annotation_thresh,
                           binding=initial_bind_w,
                           activation=initial_alpha,
                           coop=initial_coop_weight,
                           qBTM=initial_qBTM, 
                           logBias = initial_logReg_bias,
                           logCoeff = initial_logReg_coeff)
  # retreiving the index of the parameters which will change (were not fixed in the beginning) using
  #  their names
  changing_index_per_param <- list()
  #print("colnames")
  #print(colnames(compartmentSamples[[1]]))
  par_names <- unlist(lapply(strsplit(colnames(compartmentSamples[[1]]), split="_"), "[[", 1))
  par_numbers <-  as.integer(unlist(lapply(strsplit(colnames(compartmentSamples[[1]]),
                                                    split="_"), "[[", 2))) 
  # print("par_names")
  # print(par_names)
  # print("par_numbers")
  # print(par_numbers)
  for(cur_par in 1:length(template_par_set)){
    cur_par_ind <- which(par_names == names(template_par_set)[cur_par])
    changing_index_per_param[[cur_par]] <- cbind(par_numbers[cur_par_ind], cur_par_ind)
  }
  names(changing_index_per_param) <- names(template_par_set)
  #print(changing_index_per_param)
  
  # having the changing indices, now I go through the sampled parameters and create the list of
  #  parmaeters for the ensemble
  sample_parameter_list <- list()
  param_set_names <- character(0)
  counter <- 1
  for(cur_comp in 1:length(compartmentSamples)){
    for(cur_Samp in 1:nrow(compartmentSamples[[cur_comp]])){
      sample_parameter_list[[counter]] <- template_par_set
      for(cur_par in 1:length(sample_parameter_list[[counter]])){
        sample_parameter_list[[counter]][[cur_par]][changing_index_per_param[[cur_par]][, 1]] <- compartmentSamples[[cur_comp]][cur_Samp, changing_index_per_param[[cur_par]][, 2]]
      }
      param_set_names <- c(param_set_names, paste("par", cur_comp, cur_Samp, sep = "_"))
      counter <- counter + 1
    }
  }
  names(sample_parameter_list) <- param_set_names
  #check if all the parameters are within the specified ranges.
  check_vec <- logical(length(sample_parameter_list))
  for(i in 1:length(sample_parameter_list)){
    check_vec[i] <- RangeChecker(sample_parameter_list[[i]], my_range_list)
  }
  stopifnot(all(check_vec))
  return(sample_parameter_list)
}
#####################################
#####################################
#example
aa <- EnsembleParConstructor_Multi_enh(TFnames=names(TF.motifs.Expanded_new_pseudo_padded)[1:3],
                                       annotation_rangeMat=rbind(c(0.1,1), c(0.1,0.9), c(0.6,1)),
                                       annotation_logscale=F,
                                       binding_rangeMat=rbind(c(0.01,50), c(0.01,10), c(1,8)),
                                       binding_logscale=F,
                                       activation_rangeMat=rbind(c(0.01, 5), c(0.01, 10), c(0.1, 4)),
                                       activation_logscale=F,
                                       coop_rangeMat=rbind(c(0.1, 30), c(1, 8)),
                                       coop_logscale=F,
                                       qBTM_rangeMat=cbind(rep(0.001, 2), rep(1, 2)),
                                       qBTM_logscale=F,
                                       log_Reg_bias_range = c(0, 1),
                                       log_Reg_coeff_range = c(0.1, 50),
                                       log_Reg_bias_ff  = 1,
                                       log_Reg_coeff_ff = 1,
                                       nu.sample=2,
                                       annotation_thresh_ff=c(0, 0, 0),
                                       initial_bind_w_ff=c(1,1,1),
                                       initial_alpha_ff=c(1,0,1),
                                       initial_coop_weight_ff=c(0,0),
                                       initial_qBTM_ff=c(0,0),
                                       annotation_thresh=c(0.6, 0.7, 0.8),
                                       initial_bind_w=c(4,6,7),
                                       initial_alpha=c(1,2,3),
                                       initial_coop_weight=c(3,7),
                                       initial_qBTM=rep(0.1, 2), 
                                       initial_logReg_bias = 0.5,
                                       initial_logReg_coeff = 1)
aa$par_1_2
aa$par_53_1$qBTM

########################################################################################################################
########################################################################################################################
FactorInfoWriter <- function(.TF_names, role_vec, file_name ="factorInfo.txt"){
  # .TF_names is a character vector containing the names of the TFs
  # role_vec is the binary integer vector with length equal to the number of TFs, 1 shows the TF is
  #  an activator, 0 shows the TF  is a repressor
  stopifnot(length(.TF_names) == length(role_vec),
            all(role_vec %in% c(0,1)) )
  cat(c(.TF_names[1],
        as.integer(role_vec[1]==1),
        paste0(as.integer(role_vec[1]==0), "\n")),
      file = file_name, sep = "\t", append = F)
  if(length(.TF_names) > 1){
    for(tf_nu in 2:length(.TF_names)){
      cat(c(.TF_names[tf_nu],
            as.integer(role_vec[tf_nu]==1),
            paste0(as.integer(role_vec[tf_nu]==0), "\n")),
            file = file_name, sep = "\t", append = T)
      
    }
  }
}
########################################################################################################################
########################################################################################################################
# example
FactorInfoWriter(TF_names = c(names(TF.motifs.Shrinked.count)), role_vec = c(rep(1, 9), rep(0, 10)))
########################################################################################################################
########################################################################################################################

par_ff_lb_ub_coop_Writer_multi_enh <- function(.TFnames,
                                               .annotation_thresh=numeric(0),
                                               .annotation_range=numeric(0),
                                               .initial_bind_w=numeric(0),
                                               .bind_w_range=numeric(0),
                                               .initial_alpha=numeric(0),
                                               .alpha_range=numeric(0),
                                               .initial_log_reg_bias = numeric(0),
                                               .initial_log_reg_coeff = numeric(0),
                                               .initial_log_reg_bias_range = numeric(0),
                                               .initial_log_reg_coeff_range = numeric(0),
                                               TF_role_decided=F,
                                               factor_role_vector=integer(0),
                                               .coop_tf_mat=numeric(0),
                                               .initial_coop_weight=numeric(0),
                                               .coop_weight_range=numeric(0),
                                               .coop_type=character(0),
                                               .coop_dist=integer(0),
                                               .coop_orientation=integer(0),
                                               one_qbtm_per_enh=T,
                                               .initial_qBTM=numeric(0),
                                               .qBTMrange=numeric(0),
                                               one_beta_per_enh=F,
                                               .initial_pi_beta=numeric(0), 
                                               .pi_beta_upper=numeric(0),
                                               .pi_beta_lower=numeric(0),
                                               .annotation_thresh_ff=integer(0),
                                               .initial_bind_w_ff=integer(0),
                                               .initial_alpha_ff=integer(0),
                                               .initial_coop_weight_ff=integer(0),
                                               .initial_qBTM_ff=integer(0),
                                               .initial_pi_beta_ff=integer(0),
                                               .initial_log_reg_bias_ff = numeric(0),
                                               .initial_log_reg_coeff_ff = numeric(0),
                                               .filename_start=character(0),
                                               .filename_ff=character(0),
                                               .filename_upper=character(0),
                                               .filename_lower=character(0),
                                               .file_name_coop=character(0),
                                               .nu_enhacners,
                                               ensemble_mode=F,
                                               .nu_samples=0,
                                               .annotation_thresh_ff_ens=.annotation_thresh_ff,
                                               .initial_bind_w_ff_ens=.initial_bind_w_ff,
                                               .initial_alpha_ff_ens=.initial_alpha_ff,
                                               .initial_coop_weight_ff_ens=.initial_coop_weight_ff,
                                               .initial_qBTM_ff_ens=.initial_qBTM_ff,
                                               .initial_log_reg_bias_ff_ens = .initial_log_reg_bias_ff,
                                               .initial_log_reg_coeff_ff_ens = .initial_log_reg_coeff_ff,
                                               Write_to_file=T,
                                               create_folders=T,
                                               .create_bounds=T,
                                               .create_par=T,
                                               .create_ff=T,
                                               .create_coop=T,
                                               logistic_params = F
                                                 ){
  # this function will create the following files for GEMSTAT input:
  # start.par, upper.par, lower.par, freefix.txt and coop.txt
  # if the inputs are not provided it will fill in with default values
  
  # TFnames : a character vector containing the name of TFs
  # annotation_thresh : a nemeric vector (one entry per TF) containing the GEMSTAT threshold beyond 
  #  which sites are called. order corresponds to the TFnames
  # initial_bind_w : a nemeric vector (one entry per TF) containing the initial binding parameter for
  #  each TF
  # initial_alpha : a nemeric vector (one entry per TF) containing the initial activity (alpha) 
  #  parameter for each TF
  # coop_tf_mat : is a matrix where nrow = number of interacting pairs, ncol= 2 . Each row has the
  #  index of interacting TFs in one interaction (index in TFnames)
  # initial_coop_weight : a nemeric vector (one per interaction i.e. length(initial_coop_weight) ==
  #  nrow(coop_tf_mat)). each entry is the initial parameter for this interaction
  # initial_qBTM : a nemeric vector (one entry per enhancer), indiacting initialization of the the
  #  basal transcription of that enhacner
  # initial_pi_beta : a matrix with each row corresponding to one enhancer. two columns, first pi 
  #  value, second beta value for that enhacner. This can be just one row if the parameters are shared
  #  between enhacners.
  # nu_enhacners : is the number of enhancers given as input
  # annotation_thresh_ff : an integer vector (one entry per TF) containing the GEMSTAT threshold beyond
  #  which sites are called. order corresponds to the TFnames. . each value can be 0 or 1. if zero that
  #  parameter won't be optimized.
  # initial_bind_w_ff : an integer vector (one entry per TF) containing the initial binding parameter
  #  for each TF. . each value can be 0 or 1. if zero that parameter won't be optimized.
  # initial_alpha_ff : an integer vector (one entry per TF) containing the initial activity (alpha) 
  #  parameter for each TF. . each value can be 0 or 1. if zero that parameter won't be optimized.
  # initial_coop_weight_ff : an integer vector (one per interaction i.e. length(initial_coop_weight) 
  #  == nrow(coop_tf_mat)). each entry is the initial parameter for this interaction.. each value can
  #  be 0 or 1. if zero that parameter won't be optimized.
  # initial_qBTM_ff : a nemeric vector (one entry per enhancer), indiacting initialization of the the
  #  basal transcription of that enhacner. . each value can be 0 or 1. if zero that parameter won't be
  #  optimized.
  # initial_pi_beta_ff : a matrix with each row corresponding to one enhancer. two columns, first pi
  #  value, second beta value for that enhacner. This can be just one row if the parameters are shared
  #  between enhacners. each value can be 0 or 1. if zero that parameter won't be optimized.
  #### in all range matrices, first column is the lower bound, second column is the upper bound. except
  ####   the initial_pi_beta_upper, initial_pi_beta_lower
  # annotation_range : a matrix (one row per TF) containing the lower and upper bounds for each 
  #  annotation threshold. order corresponds to the TFnames
  # bind_w_range : a matrix (one row per TF) containing the lower and upper bounds for each binding
  #  weight
  # alpha_range : a matrix (one row per TF) containing the lower and upper bounds for each alpha weight
  # TF_role_decided : if True, the roles of the TFs are decided beforehand, a factorInfo file will be
  #  created and upper lower bounds, initial values and freefix will be adjusted accordingly. If False (default)
  #  no prior role is decided, factorInfo file will not be created. only the activator alpha range will
  #  be used for both repression and activation. 
  # factor_role_vector : is the binary integer vector with length equal to the number of TFs, 1 shows
  #  the TF is an activator, 0 shows the TF  is a repressor

  # coop_weight_range : a matrix (one row per interaction i.e. nrow(initial_coop_weight) == 
  #  nrow(coop_tf_mat)). each row is the lower and upper bound for the weight of that interaction
  # qBTM_range : a nemeric vector (one entry per enhancer), indiacting initialization of the the basal
  #  transcription of that enhacner
  # pi_beta_upper : a matrix with each row corresponding to one enhancer. two columns, first upper 
  #  limit for pi value, second upper limit for beta value for that enhacner. This can be just one row
  #  if the parameters are shared between enhacners.
  # pi_beta_lower : a matrix with each row corresponding to one enhancer. two columns, first lower
  #  limit for pi value, second lower limit for beta value for that enhacner. This can be just one row
  #  if the parameters are shared between enhacners.
  # coop_dist : is an integer vector, with length equal to  nrow(coop_tf_mat). it defines the distance
  #  threshold from end of first motif to the beginning of the second motif
  # coop_type : is a character vector, with length equal to  nrow(coop_tf_mat). each entry indicates
  #  the type of interaction and can take one of the following values:
  #  "SIMPLE", "DIMER", "HALF_DIRECTIONAL", "HELICAL", "HELICAL_DIRECTIONAL"
  # coop_orientation is an integer matrix, each row corresponding to an interaction, each entry can
  #  take one of three values: 0, 1, -1 which indicate the strand on which the motif appears 
  #  (only used for dimer interactions). if zero the interaction won't use this feature.
  # one_qbtm_per_enh : if True assigns one qBTM per enhancer. other wise only one for all
  # one_beta_per_enh : if True assigns one beta and one pi per enhancer. other wise only one for all
  # ensemble_mode : if True creates and ensemble (nu_samples) of parameters instead of just one
  # nu_samples : number of samples in the ensemble
  # .annotation_thresh_ff_ens : same format as .annotation_thresh_ff, if 1 that parameter will be
  #  changed in the ensemble, other wise fixed to its initial value
  # .initial_bind_w_ff_ens : same format as .initial_bind_w_ff, if 1 that parameter will be changed
  #  in the ensemble, other wise fixed to its initial value
  # .initial_alpha_ff_ens :  same format as .initial_alpha_ff, if 1 that parameter will be changed
  #  in the ensemble, other wise fixed to its initial value
  # .initial_coop_weight_ff_ens : same format as .initial_coop_weight_ff, if 1 that parameter will
  #  be changed in the ensemble, other wise fixed to its initial value
  # .initial_qBTM_ff_ens :  same format as .initial_qBTM_ff, if 1 that parameter will be changed 
  #  in the ensemble, other wise fixed to its initial value
  # Write_to_file if True writes outputs to file, otherwise returns the starting values and bounds
  # logistic_params : if True it initializes the bias and coeff parameters related to logisitic regression (if not already initialized), otherwise lets them be numeric(0)
  # check the format of inputs
  print("Checking format of inputs ...")
  stopifnot(is.character(.TFnames),
            (length(.annotation_thresh) == 0 | length(.annotation_thresh) == length(.TFnames)),
            (length(.annotation_thresh_ff) == 0 | length(.annotation_thresh_ff) == length(.TFnames)),
            (length(.annotation_range) == 0 | nrow(.annotation_range) == length(.TFnames)),
            (length(.initial_bind_w) == 0 | length(.initial_bind_w) == length(.TFnames)),
            (length(.bind_w_range) == 0 | nrow(.bind_w_range) == length(.TFnames)),
            (length(.initial_bind_w_ff) == 0 | length(.initial_bind_w_ff) == length(.TFnames)),
            (length(.initial_alpha) == 0 | length(.initial_alpha) == length(.TFnames)),
            (length(.initial_alpha_ff) == 0 | length(.initial_alpha_ff) == length(.TFnames)),
            (length(.alpha_range) == 0 | nrow(.alpha_range) == length(.TFnames)), 
            #(length(.initial_coop_weight == 0) | length(.initial_coop_weight) == nrow(.coop_tf_mat)),
            #(length(.initial_coop_weight_ff == 0) | length(.initial_coop_weight_ff) == nrow(.coop_tf_mat)),
            #(length(.coop_type == 0) | length(.coop_type) == nrow(.coop_tf_mat)),
            #(length(.coop_dist == 0) | length(.coop_dist) == nrow(.coop_tf_mat)),
            #(length(.coop_orientation== 0) | nrow(.coop_orientation) == nrow(.coop_tf_mat)),
            #(length(.coop_weight_range== 0) | nrow(.coop_weight_range) == nrow(.coop_tf_mat)),
            (length(.initial_qBTM) == 0 | length(.initial_qBTM) == .nu_enhacners | length(.initial_qBTM) == 1),
            (length(.initial_qBTM_ff) == 0 | length(.initial_qBTM_ff) == .nu_enhacners | length(.initial_qBTM_ff) == 1),
            (length(.qBTMrange) == 0 | nrow(.qBTMrange) == .nu_enhacners |  nrow(.qBTMrange) == 1),
            (length(.initial_pi_beta) == 0 | nrow(.initial_pi_beta) == .nu_enhacners | nrow(.initial_pi_beta)==1),
            (length(.initial_pi_beta_ff) == 0 | nrow(.initial_pi_beta_ff) == .nu_enhacners | nrow(.initial_pi_beta_ff)==1),
            (length(.pi_beta_upper) == 0 | nrow(.pi_beta_upper) == .nu_enhacners | nrow(.pi_beta_upper)==1),
            (length(.pi_beta_lower) == 0 | nrow(.pi_beta_lower) == .nu_enhacners | nrow(.pi_beta_lower)==1),
            length(.nu_enhacners)==1,
            length(.TFnames) > 0
           # ,((length(factor_role_vector) > 0 & TF_role_decided) | (length(factor_role_vector) == 0 & !TF_role_decided))
            )

  nu_TFs <- length(.TFnames)
  # fill in the inputs if not provided
  if(length(.annotation_thresh)==0){
    print("initializing annotation_thresh ...")
    .annotation_thresh <- rep(0.6, nu_TFs)
  }
  if(length(.annotation_range)==0){
    print("initializing annotation_range ...")
    .annotation_range <- cbind(rep(0.001, nu_TFs),
                               rep(1.0, nu_TFs))
    #print(.annotation_range)
  }
  if(length(.initial_bind_w)==0){
    print("initializing initial_bind_w ...")
    .initial_bind_w <- rep(1.0, nu_TFs)
  }
  if(length(.bind_w_range)==0){
    print("initializing bind_w_range ...")
    .bind_w_range <- cbind(rep(0.01, nu_TFs), rep(250.5, nu_TFs))
    #print(.bind_w_range)
  }
  if(length(.initial_alpha)==0){
    print("initializing initial_alpha ...")
    .initial_alpha <- rep(1.0, nu_TFs)
  }
  if(length(.alpha_range)==0){
    print("initializing alpha_range ...")
    .alpha_range <- cbind(rep(1e-6, nu_TFs), rep(100, nu_TFs))
    #print(.alpha_range)
  }
  if(length(.initial_log_reg_bias) == 0 & logistic_params){
    .initial_log_reg_bias <- 0.5
  }
  if(length(.initial_log_reg_coeff) == 0 & logistic_params){
    .initial_log_reg_coeff <- 10
  }
  #setting cooperativity parameters
  if(length(.coop_tf_mat)!=0){
    nu_ints <- nrow(.coop_tf_mat)
    if(length(.initial_coop_weight)==0){
      print("initializing .initial_coop_weight ...")
      .initial_coop_weight <- rep(1.0, nu_ints)
    }
    if(length(.coop_weight_range)==0){
      print("initializing .coop_weight_range ...")
      .coop_weight_range <- cbind(rep(0.01, nu_ints), rep(100.5, nu_ints))
    }
    if(length(.coop_type)==0){
      print("initializing .coop_type ...")
      .coop_type <- rep("SIMPLE", nu_ints)
    }
    if(length(.coop_dist)==0){
      print("initializing .coop_dist ...")
      .coop_dist <- rep(50, nu_ints)
    }
    if(length(.coop_orientation)==0){
      print("initializing .coop_orientation ...")
      .coop_orientation <- cbind(rep(0, nu_ints), rep(0, nu_ints))
    }
  }
  if(length(.initial_qBTM)==0){
    if(one_qbtm_per_enh){
      print("initializing .initial_qBTM (one for each enhancer)...")
      .initial_qBTM <- rep(0.01, .nu_enhacners)
    }else{
      print("initializing .initial_qBTM (one for all) ...")
      .initial_qBTM <- 0.01
    }
  }
  if(length(.qBTMrange)==0){
    if(one_qbtm_per_enh){
      print("initializing .qBTMrange (one for each enhancer)...")
      .qBTMrange <- cbind(rep(1e-3, .nu_enhacners), rep(1, .nu_enhacners))
    }else{
      print("initializing .qBTMrange (one for all) ...")
      .qBTMrange <- matrix(nrow=1, ncol = 2)
      .qBTMrange[1,1] <- 1e-3
      .qBTMrange[1,2] <- 1
    }
    
  }
  if(length(.initial_pi_beta)==0){
    if(one_beta_per_enh){
      print("initializing .initial_pi_beta (one for each enhancer)...")
      .initial_pi_beta <- cbind(rep(1e-10, .nu_enhacners), rep(1, .nu_enhacners))
    }else{
      print("initializing .initial_pi_beta (one for all) ...")
      .initial_pi_beta <- matrix(nrow=1, ncol = 2)
      .initial_pi_beta[1,1] <- 1e-10
      .initial_pi_beta[1,2] <- 1.0
    }
    
  }
  if(length(.pi_beta_upper) ==0){
    if(one_beta_per_enh){
      print("initializing .pi_beta_upper (one for each enhancer)...")
      .pi_beta_upper <- cbind(rep(1E10, .nu_enhacners), rep(100.0, .nu_enhacners))
    }else{
      print("initializing .pi_beta_upper (one for all) ...")
      .pi_beta_upper <- matrix(nrow=1, ncol = 2)
      .pi_beta_upper[1,1] <- 1E10
      .pi_beta_upper[1,2] <- 100.0
    }
  }
  if(length(.pi_beta_lower) ==0){
    if(one_beta_per_enh){
      print("initializing .pi_beta_lower (one for each enhancer)...")
      .pi_beta_lower <- cbind(rep(1E-75, .nu_enhacners), rep(1e-1, .nu_enhacners))
    }else{
      print("initializing .pi_beta_lower (one for all) ...")
      .pi_beta_lower <- matrix(nrow=1, ncol = 2)
      .pi_beta_lower[1,1] <- 1E-75
      .pi_beta_lower[1,2] <- 1e-1
    }
  }
  # intializing logistic reg range params
  if(length(.initial_log_reg_bias_range) == 0 & logistic_params){
    .initial_log_reg_bias_range <- c(0.01, 1)
  }
  if(length(.initial_log_reg_coeff_range) == 0 & logistic_params){
    .initial_log_reg_coeff_range <- c(1, 50)
  }
  # intializing logistic reg ff params
  if(length(.initial_log_reg_bias_ff) == 0 & logistic_params){
    .initial_log_reg_bias_ff <- 1
  }
  if(length(.initial_log_reg_coeff_ff) == 0 & logistic_params){
    .initial_log_reg_coeff_ff <- 1
  }
  if(length(.annotation_thresh_ff)==0){
    print("initializing .annotation_thresh_ff ...")
    .annotation_thresh_ff <- rep(0, nu_TFs)
    if(length(.annotation_thresh_ff_ens)==0){
      .annotation_thresh_ff_ens <- .annotation_thresh_ff
    }
  }
  if(length(.initial_bind_w_ff)==0){
    print("initializing .initial_bind_w_ff ...")
    .initial_bind_w_ff <- rep(1, nu_TFs)
    if(length(.initial_bind_w_ff_ens)==0){
      .initial_bind_w_ff_ens <- .initial_bind_w_ff
    }
  }
  if(length(.initial_alpha_ff)==0){
    print("initializing .initial_alpha_ff ...")
    .initial_alpha_ff <- rep(1, nu_TFs)
    if(length(.initial_alpha_ff_ens)==0){
      .initial_alpha_ff_ens <- .initial_alpha_ff
    }
  }
  if(length(.coop_tf_mat) != 0){
    if(length(.initial_coop_weight_ff)==0){
      print("initializing .initial_coop_weight_ff ...")
      .initial_coop_weight_ff <- rep(1, nu_ints)
      if(length(.initial_coop_weight_ff_ens)==0){
        .initial_coop_weight_ff_ens <- .initial_coop_weight_ff
      }
    }
  }
  if(length(.initial_qBTM_ff)==0){
    if(one_qbtm_per_enh){
      print("initializing .initial_qBTM_ff (one for each enhancer)...")
      .initial_qBTM_ff <- rep(1, .nu_enhacners)
      if(length(.initial_qBTM_ff_ens)==0){
        .initial_qBTM_ff_ens <- .initial_qBTM_ff
      }
    }else{
      print("initializing .initial_qBTM_ff (one for all) ...")
      .initial_qBTM_ff <- 1
      if(length(.initial_qBTM_ff_ens)==0){
        .initial_qBTM_ff_ens <- .initial_qBTM_ff
      }
    }
  }else if(length(.initial_qBTM_ff)==1){ #if only one number is provided for all qBTMs
    if(one_qbtm_per_enh){
      .initial_qBTM_ff <- rep(.initial_qBTM_ff, .nu_enhacners)
    }
  }
  
  if(length(.initial_pi_beta_ff)==0){
    if(one_beta_per_enh){
      print("initializing .initial_pi_beta_ff (one for each enhancer)...")
      .initial_pi_beta_ff <- cbind(rep(0, .nu_enhacners), rep(1, .nu_enhacners))
    }else{
      print("initializing .initial_pi_beta_ff (one for all) ...")
      .initial_pi_beta_ff <- matrix(nrow=1, ncol=2)
      .initial_pi_beta_ff[1,1] <- 0
      .initial_pi_beta_ff[1,2] <- 1
    }
  }else if(length(.initial_pi_beta_ff)==1){
    if(one_beta_per_enh){
      .initial_pi_beta_ff <- cbind(rep(0, .nu_enhacners), rep(.initial_pi_beta_ff, .nu_enhacners))
    }else{
      my_init <- .initial_pi_beta_ff
      .initial_pi_beta_ff <- matrix(nrow=1, ncol=2)
      .initial_pi_beta_ff[1,1] <- 0
      .initial_pi_beta_ff[1,2] <- my_init
    }
  }
  if(length(.filename_start)==0){
    print("initializing .filename_start ...")
    .filename_start <- "start.par"
  }
  if(length(.filename_ff)==0){
    print("initializing .filename_ff ]...")
    .filename_ff <- "ff.par"
  }
  if(length(.filename_upper)==0){
    print("initializing .filename_upper ]...")
    .filename_upper <- "upper.par"
  }
  if(length(.filename_lower)==0){
    print("initializing .filename_lower ]...")
    .filename_lower <- "lower.par"
  }
  if(length(.file_name_coop)==0){
    print("initializing .file_name_coop ]...")
    .file_name_coop <- "coop.par"
  }
  # check if all initialized parameters are within ranges (can't be equal to upper, lower bounds)
  
  if(length(.coop_weight_range) == 0){
    .coop_weight_range <- matrix(nrow =0, ncol = 2)
  }
  
  if(length(.initial_log_reg_bias_range) == 0){
    .initial_log_reg_bias_rangeMat <-  matrix(nrow =0, ncol = 2)
  }else{
    .initial_log_reg_bias_rangeMat <- matrix(.initial_log_reg_bias_range, nrow = 1)
  }
  
  if(length(.initial_log_reg_coeff_range) == 0){
    .initial_log_reg_coeff_rangeMat <-  matrix(nrow =0, ncol = 2)
  }else{
    .initial_log_reg_coeff_rangeMat <- matrix(.initial_log_reg_coeff_range, nrow = 1)
  }
  print(".qBTMrange")
  print(.qBTMrange)
  print(".initial_qBTM")
  print(.initial_qBTM)
  my_range_list <- list(.annotation_range, 
                        .bind_w_range,
                        .alpha_range,
                        .coop_weight_range,
                        .qBTMrange,
                        cbind(.pi_beta_lower[, 1],
                              .pi_beta_upper[, 1]), 
                        cbind(.pi_beta_lower[, 2],
                              .pi_beta_upper[, 2]), 
                        .initial_log_reg_bias_rangeMat,
                        .initial_log_reg_coeff_rangeMat)
  my_initial_list <- list(.annotation_thresh,
                          .initial_bind_w,
                          .initial_alpha,
                          .initial_coop_weight,
                          .initial_qBTM,
                          .initial_pi_beta[, 1],
                          .initial_pi_beta[, 2],
                          .initial_log_reg_bias, 
                          .initial_log_reg_coeff)
  range_checked <- RangeChecker(vector_list = my_initial_list,
                                range_list = my_range_list,
                                equal_allowed = F,
                                return_detailed = T)
  #return(range_checked)
  if(length(range_checked) == 1){
    if(range_checked){
      print("initial values were checked to be in the specified ranges")
    }
  }else{
    print("initial values are not in the specified ranges")
    print(range_checked)
    stop("initial values are not in the specified ranges")
  }
  if(Write_to_file){
    # write the start file
    prev_WD <- getwd()
    on.exit(setwd(prev_WD))
    if(create_folders){
      if(!dir.exists(paste0(prev_WD, "/bounds"))){
        dir.create(paste0(prev_WD, "/bounds"))
      }
      if(!dir.exists(paste0(prev_WD, "/free_fix"))){
        dir.create(paste0(prev_WD, "/free_fix"))
      }
      if(!dir.exists(paste0(prev_WD, "/Parameters"))){
        dir.create(paste0(prev_WD, "/Parameters"))
      }
    }
    if(ensemble_mode){
      print("writing the initial parameter files in ensemble mode ...")
      if(create_folders){
        setwd(paste0(prev_WD,"/Parameters"))
      }
      param_ensemble <- EnsembleParConstructor_Multi_enh(TFnames = .TFnames,
                                                         annotation_rangeMat = .annotation_range,
                                                         annotation_logscale=F,
                                                         binding_rangeMat = .bind_w_range,
                                                         binding_logscale=F,
                                                         activation_rangeMat = .alpha_range,
                                                         activation_logscale=F,
                                                         coop_rangeMat = .coop_weight_range,
                                                         coop_logscale=F,
                                                         qBTM_rangeMat = .qBTMrange,
                                                         qBTM_logscale=F,
                                                         nu.sample=.nu_samples,
                                                         log_Reg_bias_range = .initial_log_reg_bias_range,
                                                         log_Reg_coeff_range = .initial_log_reg_coeff_range,
                                                         log_Reg_bias_ff = .initial_log_reg_bias_ff_ens,
                                                         log_Reg_coeff_ff = .initial_log_reg_coeff_ff_ens,
                                                         initial_logReg_bias = .initial_log_reg_bias,
                                                         initial_logReg_coeff = .initial_log_reg_coeff,
                                                         annotation_thresh_ff = .annotation_thresh_ff_ens,
                                                         initial_bind_w_ff = .initial_bind_w_ff_ens,
                                                         initial_alpha_ff = .initial_alpha_ff_ens,
                                                         initial_coop_weight_ff = .initial_coop_weight_ff_ens,
                                                         initial_qBTM_ff = .initial_qBTM_ff_ens,
                                                         annotation_thresh=.annotation_thresh,
                                                         initial_bind_w=.initial_bind_w,
                                                         initial_alpha=.initial_alpha,
                                                         initial_coop_weight=.initial_coop_weight,
                                                         initial_qBTM=.initial_qBTM)
      if(.create_par){
        for(ens_elem in 1:length(param_ensemble)){
          initialPar_Writer_multi_enh(TFnames = .TFnames,
                                      annotation_thresh = param_ensemble[[ens_elem]][[1]],
                                      initial_bind_w = param_ensemble[[ens_elem]][[2]],
                                      initial_alpha = param_ensemble[[ens_elem]][[3]],
                                      coop_tf_mat = .coop_tf_mat,
                                      initial_coop_weight = param_ensemble[[ens_elem]][[4]],
                                      initial_qBTM=param_ensemble[[ens_elem]][[5]],
                                      initial_pi_beta = .initial_pi_beta,
                                      initial_log_reg_bias = param_ensemble[[ens_elem]]$logBias,
                                      initial_log_reg_coeff = param_ensemble[[ens_elem]]$logCoeff,
                                      nu_enhacners=.nu_enhacners,
                                      filename = paste0(names(param_ensemble)[ens_elem], ".par"))
        }
      }

    }else{
      print("writing the initial parameter file ...")
      if(create_folders){
        setwd(paste0(prev_WD,"/Parameters"))
      }
      if(.create_par){
        initialPar_Writer_multi_enh(TFnames = .TFnames,
                                    annotation_thresh = .annotation_thresh,
                                    initial_bind_w = .initial_bind_w,
                                    initial_alpha = .initial_alpha,
                                    coop_tf_mat = .coop_tf_mat,
                                    initial_coop_weight = .initial_coop_weight,
                                    initial_qBTM=.initial_qBTM,
                                    initial_pi_beta = .initial_pi_beta,
                                    initial_log_reg_bias = .initial_log_reg_bias,
                                    initial_log_reg_coeff = .initial_log_reg_coeff,
                                    nu_enhacners=.nu_enhacners,
                                    filename = .filename_start)
      }
    }
    
    if(.create_ff){
      # write the free fix file
      print("writing the free-fix file ...")
      if(create_folders){
        setwd(paste0(prev_WD, "/free_fix"))
      }
      
      freefixWriter_multi_enh(TFnames = .TFnames,
                              annotation_thresh_ff=.annotation_thresh_ff,
                              initial_bind_w_ff=.initial_bind_w_ff,
                              initial_alpha_ff=.initial_alpha_ff,
                              coop_tf_mat=.coop_tf_mat,
                              initial_coop_weight_ff=.initial_coop_weight_ff,
                              initial_qBTM_ff=.initial_qBTM_ff,
                              initial_pi_beta_ff=.initial_pi_beta_ff,
                              initial_log_reg_bias_ff = .initial_log_reg_bias_ff, 
                              initial_log_reg_coeff_ff = .initial_log_reg_coeff_ff,
                              nu_enhacners=.nu_enhacners,
                              filename = .filename_ff,
                              .TF_role_decided = TF_role_decided,
                              .factor_role_vector = factor_role_vector)
      
    }
    if(.create_bounds){
      # write the upper/lower bound file
      print("writing the upper/lower bound files ...")
      if(create_folders){
        setwd(paste0(prev_WD, "/bounds"))
      }
      upper_lower_bound_Writer_multi_enh(TFnames = .TFnames,
                                         annotation_range=.annotation_range,
                                         initial_bind_w_range=.bind_w_range,
                                         initial_alpha_range=.alpha_range,
                                         coop_tf_mat=.coop_tf_mat,
                                         initial_coop_weight_range=.coop_weight_range,
                                         initial_qBTM_range=.qBTMrange,
                                         initial_pi_beta_upper=.pi_beta_upper,
                                         initial_pi_beta_lower=.pi_beta_lower,
                                         initial_log_reg_bias_range = .initial_log_reg_bias_range,
                                         initial_log_reg_coeff_range = .initial_log_reg_coeff_range,
                                         nu_enhacners=.nu_enhacners,
                                         filename_upper=.filename_upper,
                                         filename_lower=.filename_lower,
                                         .TF_role_decided = TF_role_decided)
    }
    # write the cooperativities
    if(length(.coop_tf_mat) !=0 & .create_coop){
      print("writing the cooperativity file ...")
      if(create_folders){
        if(!dir.exists(paste0(prev_WD, "/Coop"))){
          dir.create(paste0(prev_WD, "/Coop"))
        }
        setwd(paste0(prev_WD, "/Coop"))
      }

      coop_Writer_Multi_enh(TFnames = .TFnames,
                            coop_tf_mat = .coop_tf_mat,
                            coop_type=.coop_type,
                            coop_dist=.coop_dist,
                            coop_orientation=.coop_orientation,
                            filename=.file_name_coop)
    }
    if(create_folders){
      setwd(prev_WD)
    }
  }else{ #if we don't want to write to file
    if(ensemble_mode){
      param_ensemble <- EnsembleParConstructor_Multi_enh(TFnames = .TFnames,
                                                         annotation_rangeMat = .annotation_range,
                                                         annotation_logscale=F,
                                                         binding_rangeMat = .bind_w_range,
                                                         binding_logscale=F,
                                                         activation_rangeMat = .alpha_range,
                                                         activation_logscale=F,
                                                         coop_rangeMat = .coop_weight_range,
                                                         coop_logscale=F,
                                                         qBTM_rangeMat = .qBTMrange,
                                                         qBTM_logscale=F,
                                                         nu.sample=.nu_samples,
                                                         log_Reg_bias_range = .initial_log_reg_bias_range,
                                                         log_Reg_coeff_range = .initial_log_reg_coeff_range,
                                                         log_Reg_bias_ff = .initial_log_reg_bias_ff_ens,
                                                         log_Reg_coeff_ff = .initial_log_reg_coeff_ff_ens,
                                                         initial_logReg_bias = .initial_log_reg_bias,
                                                         initial_logReg_coeff = .initial_log_reg_coeff,
                                                         annotation_thresh_ff = .annotation_thresh_ff_ens,
                                                         initial_bind_w_ff = .initial_bind_w_ff_ens,
                                                         initial_alpha_ff = .initial_alpha_ff_ens,
                                                         initial_coop_weight_ff = .initial_coop_weight_ff_ens,
                                                         initial_qBTM_ff = .initial_qBTM_ff_ens,
                                                         annotation_thresh=.annotation_thresh,
                                                         initial_bind_w=.initial_bind_w,
                                                         initial_alpha=.initial_alpha,
                                                         initial_coop_weight=.initial_coop_weight,
                                                         initial_qBTM=.initial_qBTM)
      return(list(parameter_ensemble=param_ensemble,
                  nu_enhacners=.nu_enhacners,
                  annotation_thresh_ff=.annotation_thresh_ff,
                  initial_bind_w_ff=.initial_bind_w_ff,
                  initial_alpha_ff=.initial_alpha_ff,
                  initial_coop_weight_ff=.initial_coop_weight_ff,
                  initial_qBTM_ff=.initial_qBTM_ff,
                  initial_pi_beta_ff=.initial_pi_beta_ff,
                  initial_log_reg_bias_ff = .initial_log_reg_bias_ff,
                  initial_log_reg_coeff_ff = .initial_log_reg_coeff_ff,
                  annotation_range=.annotation_range,
                  initial_bind_w_range=.bind_w_range,
                  initial_alpha_range=.alpha_range,
                  initial_coop_weight_range=.coop_weight_range,
                  initial_qBTM_range=.qBTMrange,
                  initial_pi_beta_upper=.pi_beta_upper,
                  initial_pi_beta_lower=.pi_beta_lower,
                  initial_log_reg_bias_range = .initial_log_reg_bias_range,
                  initial_log_reg_coeff_range = .initial_log_reg_coeff_range,
                  coop_type=.coop_type,
                  coop_dist=.coop_dist,
                  coop_orientation=.coop_orientation)
             )
    }else{
      return(list(TFnames=.TFnames,
                  annotation_thresh = .annotation_thresh,
                  initial_bind_w = .initial_bind_w,
                  initial_alpha = .initial_alpha,
                  initial_log_reg_bias = .initial_log_reg_bias,
                  initial_log_reg_coeff = .initial_log_reg_coeff,
                  coop_tf_mat = .coop_tf_mat,
                  initial_coop_weight = .initial_coop_weight,
                  initial_qBTM=.initial_qBTM,
                  initial_pi_beta = .initial_pi_beta,
                  nu_enhacners=.nu_enhacners,
                  annotation_thresh_ff=.annotation_thresh_ff,
                  initial_bind_w_ff=.initial_bind_w_ff,
                  initial_alpha_ff=.initial_alpha_ff,
                  initial_coop_weight_ff=.initial_coop_weight_ff,
                  initial_qBTM_ff=.initial_qBTM_ff,
                  initial_pi_beta_ff=.initial_pi_beta_ff,
                  initial_log_reg_bias_ff = .initial_log_reg_bias_ff,
                  initial_log_reg_coeff_ff = .initial_log_reg_coeff_ff,
                  annotation_range=.annotation_range,
                  initial_bind_w_range=.bind_w_range,
                  initial_alpha_range=.alpha_range,
                  initial_coop_weight_range=.coop_weight_range,
                  initial_qBTM_range=.qBTMrange,
                  initial_pi_beta_upper=.pi_beta_upper,
                  initial_pi_beta_lower=.pi_beta_lower,
                  initial_log_reg_bias_range = .initial_log_reg_bias_range,
                  initial_log_reg_coeff_range = .initial_log_reg_coeff_range,
                  coop_type=.coop_type,
                  coop_dist=.coop_dist,
                  coop_orientation=.coop_orientation))
    }

  }
}
#############################################################################
#############################################################################
#example
par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=names(TF.motifs.Expanded_new_pseudo_padded)[1:3],
                                   .annotation_thresh=numeric(0),
                                   .annotation_range=numeric(0),
                                   .initial_bind_w=numeric(0),
                                   .bind_w_range=numeric(0),
                                   .initial_alpha=numeric(0),
                                   .initial_log_reg_bias = numeric(0),
                                   .initial_log_reg_coeff = numeric(0),
                                   .initial_log_reg_bias_range = numeric(0),
                                   .initial_log_reg_coeff_range = numeric(0),
                                   .alpha_range=numeric(0),
                                   .coop_tf_mat=rbind(c(1,2), c(2, 3)),
                                   .initial_coop_weight=numeric(0),
                                   .coop_weight_range=numeric(0),
                                   .coop_type=character(0),
                                   .coop_dist=integer(0),
                                   .coop_orientation=integer(0),
                                   one_qbtm_per_enh=F,
                                   .initial_qBTM=numeric(0),
                                   .qBTMrange=numeric(0),
                                   one_beta_per_enh=F,
                                   .initial_pi_beta=numeric(0), 
                                   .pi_beta_upper=numeric(0),
                                   .pi_beta_lower=numeric(0),
                                   .annotation_thresh_ff=integer(0),
                                   .initial_bind_w_ff=integer(0),
                                   .initial_alpha_ff=integer(0),
                                   .initial_coop_weight_ff=integer(0),
                                   .initial_qBTM_ff=integer(0),
                                   .initial_pi_beta_ff=integer(0),
                                   .filename_start=character(0),
                                   .filename_ff=character(0),
                                   .filename_upper=character(0),
                                   .filename_lower=character(0),
                                   .file_name_coop=character(0),
                                   .initial_log_reg_bias_ff = numeric(0),
                                   .initial_log_reg_coeff_ff = numeric(0),
                                   .nu_enhacners=5,
                                   ensemble_mode = T,
                                   .nu_samples = 2,
                                   .annotation_thresh_ff_ens = c(0,0,0),
                                   .initial_coop_weight_ff_ens = c(0,0),
                                   .initial_qBTM_ff_ens = 0,
                                   .initial_log_reg_bias_ff_ens = 0,
                                   .initial_log_reg_coeff_ff_ens = 0,
                                   Write_to_file = T,
                                   create_folders = T,
                                   TF_role_decided = F,
                                   factor_role_vector=integer(0),  
                                   .create_bounds=T,
                                   .create_par=T,
                                   .create_ff=T,
                                   .create_coop=T,
                                   logistic_params = T)

########################################################################################################################
########################################################################################################################
initialPar_Writer_multi_enh <- function(TFnames,
                                        annotation_thresh,
                                        initial_bind_w,
                                        initial_alpha,
                                        coop_tf_mat,
                                        initial_coop_weight,
                                        initial_qBTM,
                                        initial_pi_beta,
                                        nu_enhacners,
                                        initial_log_reg_bias,
                                        initial_log_reg_coeff,
                                        filename
                                        # ,.TF_role_decided=F,
                                        # .factor_role_vector=numeric(0)
                                        ){
  # TFnames : a character vector containing the name of TFs
  # annotation_thresh : a nemeric vector (one entry per TF) containing the GEMSTAT threshold beyond
  #  which sites are called. order corresponds to the TFnames
  # initial_bind_w : a nemeric vector (one entry per TF) containing the initial binding parameter for
  #  each TF
  # initial_alpha : a nemeric vector (one entry per TF) containing the initial activity (alpha) 
  #  parameter for each TF
  # coop_tf_mat : is a matrix where nrow = number of interacting pairs, ncol= 2 . Each row has the
  #  index of interacting TFs in one interaction (index in TFnames)
  # initial_coop_weight : a nemeric vector (one per interaction i.e. length(initial_coop_weight) ==
  #  nrow(coop_tf_mat)). each entry is the initial parameter for this interaction
  # initial_qBTM : a nemeric vector (one entry per enhancer), indiacting initialization of the the
  #  basal transcription of that enhacner
  # initial_pi_beta : a matrix with each row corresponding to one enhancer. two columns, first pi
  #  value, second beta value for that enhacner. This can be just one row if the parameters are shared
  #  between enhacners.
  # nu_enhacners : is the number of enhancers given as input
  # filename : name of the initial parameter file.
  # .TF_role_decided : if True, the roles of the TFs are decided beforehand, a factorInfo file will be
  #  created and upper lower bounds, initial values and freefix will be adjusted accordingly. If False (default)
  #  no prior role is decided, factorInfo file will not be created. only the activator alpha range will
  #  be used for both repression and activation. 
  # .factor_role_vector : is the binary integer vector with length equal to the number of TFs, 1 shows
  #  the TF is an activator, 0 shows the TF  is a repressor
  # initial_log_reg_bias : is a positive float, indicating the initial bias of the logisitic classifier if available
  # initial_log_reg_coeff : is a positive float, indicating the initial coeff of the logisitic classifier if available
  # check the format of the inputs
  stopifnot(length(initial_coop_weight) == nrow(coop_tf_mat),
            is.character(TFnames),
            is.numeric(annotation_thresh),
            is.numeric(initial_bind_w),
            is.numeric(initial_alpha),
            length(annotation_thresh) == length(TFnames),
            length(initial_bind_w) == length(TFnames),
            length(initial_alpha) == length(TFnames),
            (length(initial_qBTM) == 1 | length(initial_qBTM) == nu_enhacners),
            (nrow(initial_pi_beta) == 1 | nrow(initial_pi_beta) == nu_enhacners),
            (length(initial_log_reg_bias) <= 1), 
            (length(initial_log_reg_coeff) <= 1),
            length(initial_log_reg_bias) == length(initial_log_reg_coeff)
            )
  
  # write parameters for the TFs
  cat('{"tfs":{', file = filename, sep = "")
  for(tf_nu in 1:length(TFnames)){
    cat(c('"', TFnames[tf_nu], '"',
          ':{"annot_thresh":', annotation_thresh[tf_nu],
          ',"maxbind":', initial_bind_w[tf_nu],
          ',"alpha_a":', ifelse(test = initial_alpha[tf_nu] >= 1, yes = initial_alpha[tf_nu], no = 1),
          ',"alpha_r":', ifelse(test = initial_alpha[tf_nu] <= 1, yes = initial_alpha[tf_nu], no = 1),
          '}', rep(',', as.integer(tf_nu != length(TFnames)))), sep = "", file = filename, append = T)
  }
  cat('}', file = filename, sep = "", append = T)
  
  # writing cooperativities
  if(length(initial_coop_weight) > 0){
    cat(',"inter":{', sep = "", file = filename, append = T)
    for(int_nu in 1:nrow(coop_tf_mat)){
      cat(c('"', TFnames[coop_tf_mat[int_nu, 1]],
            ":", TFnames[coop_tf_mat[int_nu, 2]], '":',
            initial_coop_weight[int_nu],
            rep(',', as.integer(int_nu != nrow(coop_tf_mat)))), sep = "", file = filename, append = T)
    }
    cat('}', file = filename, sep = "", append = T)
  }
  
  #writing qbtms
  cat(',"qbtm":[', sep = "", file = filename, append = T)
  for(cur_qbtm in 1:length(initial_qBTM)){
    cat(c(initial_qBTM[cur_qbtm],
          rep(',', as.integer(cur_qbtm != length(initial_qBTM)))), sep="", file = filename, append = T)
  }
  cat(c('],', '\n'), sep = "", file = filename, append = T)
  
  #writing pi and beta for enhancers
  cat(c('"enh":[', "\n"), sep="", file=filename, append = T)
  for(cur_enh in 1:nrow(initial_pi_beta)){
    cat(c('{"pi":', initial_pi_beta[cur_enh, 1],
          ',"beta":', initial_pi_beta[cur_enh, 2], '}',
          rep(',', as.integer(cur_enh != nrow(initial_pi_beta))),
          rep('\n', as.integer(cur_enh == nrow(initial_pi_beta)))),
        sep="", file = filename, append = T)
  }
  cat(c(paste0(']', rep(",", length(initial_log_reg_bias))), '\n'), sep="", file = filename, append = T)
  if(length(initial_log_reg_bias) > 0){
    cat(c('"log_Reg":[', "\n"), sep="", file=filename, append = T)
    cat(c('{"bias":', initial_log_reg_bias,
          ',"coeff":', initial_log_reg_coeff, '}\n'),
        sep="", file = filename, append = T)
    cat(c(']', '\n'), sep="", file = filename, append = T)
  }
  cat('}', file = filename, sep = "", append = T)
}
########################################################################################################################
# library(rjson)
# aa <- fromJSON(file="Seeded_GEMSTAT_ens/Experiment_1_to_50/Experiment_11/Inputs/Parameters/start_10.par")
# aa2 <- fromJSON(file="Seeded_GEMSTAT_ens/Experiment_2129/Inputs/Parameters/start_10.par")
########################################################################################################################
########################################################################################################################
#example
initialPar_Writer_multi_enh(TFnames=names(TF.motifs.Expanded_new_pseudo_padded)[1:3],
                            annotation_thresh=c(0.6, 0.7, 0.8),
                            initial_bind_w=c(4,6,7),
                            initial_alpha=c(1,2,3),
                            coop_tf_mat=rbind(c(1,2),c(2,3)),
                            initial_coop_weight=c(3,7),
                            initial_qBTM=rep(0.1, 5),
                            initial_pi_beta=cbind(rep(0.001, 5), rep(1, 5)),
                            nu_enhacners=5,
                            initial_log_reg_bias = 0.5,
                            initial_log_reg_coeff = 10,
                            filename="start_test.txt")
aa <- fromJSON(file="start_test.txt")

########################################################################################################################
########################################################################################################################

freefixWriter_multi_enh <- function(TFnames,
                                    annotation_thresh_ff,
                                    initial_bind_w_ff,
                                    initial_alpha_ff,
                                    coop_tf_mat,
                                    initial_coop_weight_ff,
                                    initial_qBTM_ff,
                                    initial_pi_beta_ff,
                                    initial_log_reg_bias_ff,
                                    initial_log_reg_coeff_ff,
                                    nu_enhacners,
                                    filename
                                    ,.TF_role_decided=F,
                                    .factor_role_vector=numeric(0)
                                    ){
  # TFnames : a character vector containing the name of TFs
  # annotation_thresh_ff : an integer vector (one entry per TF) containing the GEMSTAT threshold
  #  beyond which sites are called. order corresponds to the TFnames. . each value can be 0 or 1. if
  #  zero that parameter won't be optimized.
  # initial_bind_w_ff : an integer vector (one entry per TF) containing the initial binding parameter
  #  for each TF. . each value can be 0 or 1. if zero that parameter won't be optimized.
  # initial_alpha_ff : an integer vector (one entry per TF) containing the initial activity (alpha)
  #  parameter for each TF. . each value can be 0 or 1. if zero that parameter won't be optimized.
  # coop_tf_mat : is a matrix where nrow = number of interacting pairs, ncol= 2 . Each row has the
  #  index of interacting TFs in one interaction (index in TFnames)
  # initial_coop_weight_ff : an integer vector (one per interaction i.e. length(initial_coop_weight) ==
  #  nrow(coop_tf_mat)). each entry is the initial parameter for this interaction.. each value can be
  #  0 or 1. if zero that parameter won't be optimized.
  # initial_qBTM_ff : a nemeric vector (one entry per enhancer), indiacting initialization of 
  #  the basal transcription of that enhacner. . each value can be 0 or 1. if zero that parameter 
  #  won't be optimized.
  # initial_pi_beta_ff : a matrix with each row corresponding to one enhancer. two columns, first
  #  pi value, second beta value for that enhacner. This can be just one row if the parameters are
  #  shared between enhacners. each value can be 0 or 1. if zero that parameter won't be optimized.
  # initial_log_reg_bias_ff : 0 or 1. if zero  parameter won't be optimized.
  # initial_log_reg_coeff_ff : 0 or 1. if zero  parameter won't be optimized.
  # nu_enhacners : is the number of enhancers given as input
  # filename : name of the initial parameter file.
  # .TF_role_decided : if True, the roles of the TFs are decided beforehand, a factorInfo file will be
  #  created and upper lower bounds, initial values and freefix will be adjusted accordingly. If False (default)
  #  no prior role is decided, factorInfo file will not be created. only the activator alpha range will
  #  be used for both repression and activation. 
  # .factor_role_vector : is the binary integer vector with length equal to the number of TFs, 1 shows
  #  the TF is an activator, 0 shows the TF  is a repressor
  
  # check the format of the inputs
  stopifnot(length(initial_coop_weight_ff) == nrow(coop_tf_mat),
            is.character(TFnames),
            #is.integer(annotation_thresh_ff),
            #is.integer(initial_bind_w_ff),
            #is.integer(initial_alpha_ff),
            length(annotation_thresh_ff) == length(TFnames),
            length(initial_bind_w_ff) == length(TFnames),
            length(initial_alpha_ff) == length(TFnames),
            (length(initial_qBTM_ff) == 1 | length(initial_qBTM_ff) == nu_enhacners),
            (nrow(initial_pi_beta_ff) == 1 | nrow(initial_pi_beta_ff) == nu_enhacners)
            ,all(cbind(annotation_thresh_ff, initial_bind_w_ff, initial_alpha_ff) %in% c(0, 1)),
            (all(initial_coop_weight_ff %in% c(0, 1)) | length(coop_tf_mat) == 0),
            all(initial_qBTM_ff %in% c(0, 1)),
            all(initial_pi_beta_ff %in% c(0, 1))
  )
  
  # write parameters for the TFs
  cat('{"tfs":{', file = filename, sep = "")
  for(tf_nu in 1:length(TFnames)){
    if(length(.factor_role_vector) > 0){
      my_act <- as.integer(.factor_role_vector[tf_nu] == 1)
      my_rep <- as.integer(.factor_role_vector[tf_nu] == 0)
    }else if(.TF_role_decided){
      my_act <- initial_alpha_ff[tf_nu]
      my_rep <- initial_alpha_ff[tf_nu]
    }else{
      my_act <- initial_alpha_ff[tf_nu]
      my_rep <- 0
    }
    cat(c('"', TFnames[tf_nu], '"',
          ':{"annot_thresh":', annotation_thresh_ff[tf_nu],
          ',"maxbind":', initial_bind_w_ff[tf_nu],
          ',"alpha_a":', my_act,
          ',"alpha_r":', my_rep,
          '}', rep(',', as.integer(tf_nu != length(TFnames)))), sep = "", file = filename, append = T)
  }
  cat('}', file = filename, sep = "", append = T)
  
  # writing cooperativities
  if(length(coop_tf_mat) > 0){
    cat(',"inter":{', sep = "", file = filename, append = T)
    for(int_nu in 1:nrow(coop_tf_mat)){
      cat(c('"', TFnames[coop_tf_mat[int_nu, 1]],
            ":", TFnames[coop_tf_mat[int_nu, 2]], '":',
            initial_coop_weight_ff[int_nu],
            rep(',', as.integer(int_nu != nrow(coop_tf_mat)))), sep = "", file = filename, append = T)
    }
    cat('}', file = filename, sep = "", append = T)
  }
  
  #writing qbtms
  cat(',"qbtm":[', sep = "", file = filename, append = T)
  for(cur_qbtm in 1:length(initial_qBTM_ff)){
    cat(c(initial_qBTM_ff[cur_qbtm],
          rep(',', as.integer(cur_qbtm != length(initial_qBTM_ff)))), sep="", file = filename, append = T)
  }
  cat(c('],', '\n'), sep = "", file = filename, append = T)
  
  #writing pi and beta for enhancers
  cat(c('"enh":[', "\n"), sep="", file=filename, append = T)
  for(cur_enh in 1:nrow(initial_pi_beta_ff)){
    cat(c('{"pi":', initial_pi_beta_ff[cur_enh, 1],
          ',"beta":', initial_pi_beta_ff[cur_enh, 2], '}',
          rep(',', as.integer(cur_enh != nrow(initial_pi_beta_ff))),
          rep('\n', as.integer(cur_enh == nrow(initial_pi_beta_ff)))),
        sep="", file = filename, append = T)
  }
  cat(c(paste0(']', rep(",", length(initial_log_reg_bias_ff))), '\n'), 
      sep="", file = filename, append = T)
  if(length(initial_log_reg_bias_ff) > 0){
    cat(c('"log_Reg":[', "\n"), sep="", file=filename, append = T)
    cat(c('{"bias":', initial_log_reg_bias_ff,
          ',"coeff":', initial_log_reg_coeff_ff, '}\n'),
        sep="", file = filename, append = T)
    cat(c(']', '\n'), sep="", file = filename, append = T)
  }
  cat('}', file = filename, sep = "", append = T)
}
########################################################################################################################
########################################################################################################################
#example
freefixWriter_multi_enh(TFnames=names(TF.motifs.Expanded_new_pseudo_padded)[1:3],
                        annotation_thresh_ff=c(0, 0, 0),
                        initial_bind_w_ff=c(1,1,1),
                        initial_alpha_ff=c(1,1,1),
                        coop_tf_mat=numeric(0),
                        initial_coop_weight_ff=c(1, 1),
                        initial_qBTM_ff = rep(1, 5),
                        initial_pi_beta_ff = cbind(rep(0, 5), rep(1, 5)),
                        nu_enhacners=5,
                        initial_log_reg_bias_ff = 1, 
                        initial_log_reg_coeff_ff = 1,
                        filename= "ff_test.txt")
aa <- fromJSON(file = "ff_test.txt")
freefixWriter_multi_enh(TFnames=names(TF.motifs.Expanded_new_pseudo_padded)[1:3],
                        annotation_thresh_ff=c(0, 0, 0),
                        initial_bind_w_ff=c(1,1,1),
                        initial_alpha_ff=c(1,1,1),
                        coop_tf_mat=numeric(0),
                        initial_coop_weight_ff=c(1, 1),
                        initial_qBTM_ff = rep(1, 5),
                        initial_pi_beta_ff = matrix(c(0, 1), nrow = 1),
                        nu_enhacners=5,
                        initial_log_reg_bias_ff = 1, 
                        initial_log_reg_coeff_ff = 1,
                        filename= "ff_test.txt")
########################################################################################################################
########################################################################################################################
upper_lower_bound_Writer_multi_enh <- function(TFnames,
                                               annotation_range,
                                               initial_bind_w_range,
                                               initial_alpha_range,
                                               coop_tf_mat,
                                               initial_coop_weight_range,
                                               initial_qBTM_range,
                                               initial_pi_beta_upper,
                                               initial_pi_beta_lower,
                                               initial_log_reg_bias_range,
                                               initial_log_reg_coeff_range,
                                               nu_enhacners,
                                               filename_upper,
                                               filename_lower
                                               ,.TF_role_decided=F,
                                               .factor_role_vector=numeric(0)
                                               ){
  # TFnames : a character vector containing the name of TFs
  # in all range matrices, first column is the lower bound, second column is the upper bound. except
  #  the initial_pi_beta_upper, initial_pi_beta_lower
  # annotation_range : a matrix (one row per TF) containing the lower and upper bounds for each
  #  annotation threshold. order corresponds to the TFnames
  # initial_bind_w_range : a matrix (one row per TF) containing the lower and upper bounds for each
  #  binding weight
  # initial_alpha_range : a matrix (one row per TF) containing the lower and upper bounds for each
  #  alpha weight
  # coop_tf_mat : is a matrix where nrow = number of interacting pairs, ncol= 2 . Each row has the 
  #  index of interacting TFs in one interaction (index in TFnames)
  # initial_coop_weight_range : a matrix (one row per interaction i.e. nrow(initial_coop_weight) ==
  #  nrow(coop_tf_mat)). each row is the lower and upper bound for the weight of that interaction
  # initial_qBTM_range : a nemeric vector (one entry per enhancer), indiacting initialization of the
  #  the basal transcription of that enhacner
  # initial_pi_beta_upper : a matrix with each row corresponding to one enhancer. two columns, first
  #  upper limit for pi value, second upper limit for beta value for that enhacner. This can be just
  #  one row if the parameters are shared between enhacners.
  # initial_pi_beta_lower : a matrix with each row corresponding to one enhancer. two columns, first
  #  lower limit for pi value, second lower limit for beta value for that enhacner. This can be just
  #  one row if the parameters are shared between enhacners.
  # initial_log_reg_bias_range : a numeric vector of length 2. first and second entry are the lower and upper bound  for bias parameter, respectively. # related to logisitic regression
  # initial_log_reg_coeff_range : a numeric vector of length 2. first and second entry are the lower and upper bound  for coeff parameter, respectively.# related to logisitic regression
  # nu_enhacners : is the number of enhancers given as input
  # filename : name of the initial parameter file.
  # .TF_role_decided : if True, the roles of the TFs are decided beforehand, a factorInfo file will be
  #  created and upper lower bounds, initial values and freefix will be adjusted accordingly. If False (default)
  #  no prior role is decided, factorInfo file will not be created. only the activator alpha range will
  #  be used for both repression and activation. 
  # .factor_role_vector : is the binary integer vector with length equal to the number of TFs, 1 shows
  #  the TF is an activator, 0 shows the TF  is a repressor
  
  # check the format of the inputs
  stopifnot(nrow(initial_coop_weight_range) == nrow(coop_tf_mat),
            is.character(TFnames),
            nrow(annotation_range) == length(TFnames),
            nrow(initial_bind_w_range) == length(TFnames),
            nrow(initial_alpha_range) == length(TFnames),
            (nrow(initial_qBTM_range) == 1 | nrow(initial_qBTM_range) == nu_enhacners),
            (nrow(initial_pi_beta_upper) == 1 | nrow(initial_pi_beta_upper) == nu_enhacners),
            (nrow(initial_pi_beta_lower) == 1 | nrow(initial_pi_beta_lower) == nu_enhacners), 
            (length(initial_log_reg_bias_range) == 0 | length(initial_log_reg_bias_range) == 2),
            (length(initial_log_reg_coeff_range) == 0 | length(initial_log_reg_coeff_range) == 2),
            length(initial_log_reg_coeff_range) == length(initial_log_reg_bias_range)
  )
  
  
  # write bounds for the TFs parameters
  cat('{"tfs":{', file = filename_lower, sep = "")
  cat('{"tfs":{', file = filename_upper, sep = "")
  for(tf_nu in 1:length(TFnames)){
    cat(c('"', TFnames[tf_nu], '"',
          ':{"annot_thresh":', annotation_range[tf_nu, 1],
          ',"maxbind":', initial_bind_w_range[tf_nu, 1],
          ',"alpha_a":', ifelse(test = .TF_role_decided, yes = 0.999, no = initial_alpha_range[tf_nu, 1]),
          ',"alpha_r":', initial_alpha_range[tf_nu, 1],
          '}', rep(',', as.integer(tf_nu != length(TFnames)))),
        sep = "", file = filename_lower, append = T)
  }
  for(tf_nu in 1:length(TFnames)){
    cat(c('"', TFnames[tf_nu], '"',
          ':{"annot_thresh":', annotation_range[tf_nu, 2],
          ',"maxbind":', initial_bind_w_range[tf_nu, 2],
          ',"alpha_a":', initial_alpha_range[tf_nu, 2],
          ',"alpha_r":', ifelse(test = .TF_role_decided, yes = 1.001, no = initial_alpha_range[tf_nu, 2]),
          '}', rep(',', as.integer(tf_nu != length(TFnames)))),
        sep = "", file = filename_upper, append = T)
  }
  cat('}', file = filename_lower, sep = "", append = T)
  cat('}', file = filename_upper, sep = "", append = T)
  
  # writing bounds for cooperativities
  if(length(coop_tf_mat) > 0){
    cat(',"inter":{', sep = "", file = filename_lower, append = T)
    cat(',"inter":{', sep = "", file = filename_upper, append = T)
    for(int_nu in 1:nrow(coop_tf_mat)){
      cat(c('"', TFnames[coop_tf_mat[int_nu, 1]],
            ":", TFnames[coop_tf_mat[int_nu, 2]], '":',
            initial_coop_weight_range[int_nu, 1],
            rep(',', as.integer(int_nu != nrow(coop_tf_mat)))),
          sep = "", file = filename_lower, append = T)
      cat(c('"', TFnames[coop_tf_mat[int_nu, 1]],
            ":", TFnames[coop_tf_mat[int_nu, 2]], '":',
            initial_coop_weight_range[int_nu, 2],
            rep(',', as.integer(int_nu != nrow(coop_tf_mat)))),
          sep = "", file = filename_upper, append = T)
    }
    cat('}', file = filename_lower, sep = "", append = T)
    cat('}', file = filename_upper, sep = "", append = T)
  }

  #writing bounds on qbtms
  cat(',"qbtm":[', sep = "", file = filename_lower, append = T)
  cat(',"qbtm":[', sep = "", file = filename_upper, append = T)
  for(cur_qbtm in 1:nrow(initial_qBTM_range)){
    cat(c(initial_qBTM_range[cur_qbtm, 1],
          rep(',', as.integer(cur_qbtm != nrow(initial_qBTM_range)))),
        sep="", file = filename_lower, append = T)
    cat(c(initial_qBTM_range[cur_qbtm, 2],
          rep(',', as.integer(cur_qbtm != nrow(initial_qBTM_range)))),
        sep="", file = filename_upper, append = T)
  }
  cat(c('],', '\n'), sep = "", file = filename_lower, append = T)
  cat(c('],', '\n'), sep = "", file = filename_upper, append = T)
  #writing bounds on pi and beta for enhancers
  cat(c('"enh":[', "\n"), sep="", file=filename_lower, append = T)
  cat(c('"enh":[', "\n"), sep="", file=filename_upper, append = T)
  for(cur_enh in 1:nrow(initial_pi_beta_lower)){
    cat(c('{"pi":', initial_pi_beta_lower[cur_enh, 1],
          ',"beta":', initial_pi_beta_lower[cur_enh, 2], '}',
          rep(',', as.integer(cur_enh != nrow(initial_pi_beta_lower))),
          rep('\n', as.integer(cur_enh == nrow(initial_pi_beta_lower)))),
        sep="", file = filename_lower, append = T)
    cat(c('{"pi":', initial_pi_beta_upper[cur_enh, 1],
          ',"beta":', initial_pi_beta_upper[cur_enh, 2], '}',
          rep(',', as.integer(cur_enh != nrow(initial_pi_beta_upper))),
          rep('\n', as.integer(cur_enh == nrow(initial_pi_beta_upper)))),
        sep="", file = filename_upper, append = T)
  }
  # cat(c(']', '\n'), sep="", file = filename_lower, append = T)
  # cat(c(']', '\n'), sep="", file = filename_upper, append = T)
  
  cat(c(paste0(']', rep(",", length(initial_log_reg_bias_range)/2)), '\n'), 
      sep="", file = filename_lower, append = T)
  cat(c(paste0(']', rep(",", length(initial_log_reg_bias_range)/2)), '\n'), 
      sep="", file = filename_upper, append = T)
  if(length(initial_log_reg_bias_range) > 0){
    cat(c('"log_Reg":[', "\n"), sep="", file=filename_lower, append = T)
    cat(c('"log_Reg":[', "\n"), sep="", file=filename_upper, append = T)
    cat(c('{"bias":', initial_log_reg_bias_range[1],
          ',"coeff":', initial_log_reg_coeff_range[1], '}\n'),
        sep="", file = filename_lower, append = T)
    cat(c('{"bias":', initial_log_reg_bias_range[2],
          ',"coeff":', initial_log_reg_coeff_range[2], '}\n'),
        sep="", file = filename_upper, append = T)
    cat(c(']', '\n'), sep="", file = filename_lower, append = T)
    cat(c(']', '\n'), sep="", file = filename_upper, append = T)
  }
  
  cat('}', file = filename_lower, sep = "", append = T)
  cat('}', file = filename_upper, sep = "", append = T)
}
########################################################################################################################
########################################################################################################################
#example
upper_lower_bound_Writer_multi_enh(TFnames=names(TF.motifs.Expanded_new_pseudo_padded)[1:3],
                                   annotation_range=rbind(c(0.1,1), c(0.1,0.9), c(0.6,1)),
                                   initial_bind_w_range = rbind(c(0.01,50), c(0.01,10), c(1,3)),
                                   initial_alpha_range=rbind(c(0.01, 5), c(0.01, 10), c(0.1, 3)),
                                   coop_tf_mat=rbind(c(1,2),c(2,3)),
                                   initial_coop_weight_range=rbind(c(0.1, 30), c(1, 8)),
                                   initial_qBTM_range=cbind(rep(0.001, 5), rep(0.1, 5)),
                                   initial_pi_beta_upper=cbind(rep(100, 5), rep(10, 5)),
                                   initial_pi_beta_lower=cbind(rep(0.001, 5), rep(0.01, 5)),
                                   initial_log_reg_bias_range = c(0, 1), 
                                   initial_log_reg_coeff_range = c(0.01, 100),
                                   nu_enhacners=5,
                                   filename_upper="upper_test.par",
                                   filename_lower="lower_test.par")
aa <- fromJSON(file = "lower_test.par")
aa$log_Reg

########################################################################################################################
########################################################################################################################
coop_Writer_Multi_enh <- function(TFnames,
                                  coop_tf_mat,
                                  coop_type,
                                  coop_dist,
                                  coop_orientation,
                                  filename){
  # TFnames is a character vector containing the name of the TFs in the same order that their 
  #  expression matrix and motif file is constructed
  # coop_tf_mat : is a matrix where each row corresponds to an interaction, first column is the
  #  index of the first TF, second column is the index of the second TF
  # coop_dist : is an integer vector, with length equal to  nrow(coop_tf_mat). it defines the distance
  #  threshold from end of first motif to the beginning of the second motif
  # coop_type : is a character vector, with length equal to  nrow(coop_tf_mat). each entry indicates 
  #  the type of interaction and can take one of the following values:
  #  "SIMPLE", "DIMER", "HALF_DIRECTIONAL", "HELICAL", "HELICAL_DIRECTIONAL"
  # coop_orientation is an integer matrix, each row corresponding to an interaction, each entry can
  #  take one of three values: 0, 1, -1 which indicate the strand on which the motif appears (only 
  #  used for dimer interactions). if zero the interaction won't use this feature.
  
  stopifnot(is.character(TFnames),
            (length(TFnames) > 0),
            (nrow(coop_tf_mat) > 0),
            (ncol(coop_tf_mat)==2),
            (sum(coop_tf_mat > length(TFnames))==0),
            (length(coop_type) == nrow(coop_tf_mat)),
            (sum(coop_type %in% c("SIMPLE", "DIMER", 
                                  "HALF_DIRECTIONAL",
                                  "HELICAL",
                                  "HELICAL_DIRECTIONAL")) == length(coop_type)),
            (length(coop_dist) == nrow(coop_tf_mat)),
            (nrow(coop_orientation)==nrow(coop_tf_mat)),
            (ncol(coop_orientation)==2),
            (sum(coop_orientation %in% c(-1, 0, 1)) == nrow(coop_tf_mat) * 2),
            (sum(abs(coop_orientation[coop_type == "SIMPLE",])) == 0), 
            is.character(filename),
            (length(filename)==1)
            )
  
  
  for(cptf in 1:nrow(coop_tf_mat)){
    cat(paste(TFnames[coop_tf_mat[cptf, 1]], 
              "\t",
              TFnames[coop_tf_mat[cptf, 2]],
              "\t",
              coop_type[cptf],
              "\t",
              coop_dist[cptf],
              #if(coop_orientation[cptf, ] %in% c(1, -1))
              rep("\t+",as.integer(coop_orientation[cptf, 1] == 1)),
              rep("\t-",as.integer(coop_orientation[cptf, 1] == -1)),
              rep("\t+",as.integer(coop_orientation[cptf, 2] == 1)),
              rep("\t-",as.integer(coop_orientation[cptf, 2] == -1)),
              sep = ""),
        file = filename, sep = "\n", append = T)
  }
  
}
########################################################################################################################
########################################################################################################################
#example
coop_Writer_Multi_enh(TFnames=names(TF.motifs.Shrinked.t),
                      coop_tf_mat=rbind(c(1,3), c(4, 5), c(14, 8)),
                      coop_type = c('DIMER', 'DIMER', 'SIMPLE' ),
                      coop_dist=c(10,20,30),
                      coop_orientation=rbind(c(1,-1), c(1, 1), c(0, 0)),
                      filename="test_coop.txt")
########################################################################################################################
########################################################################################################################
Hal_job_writer_Multi_enh <- function(exp.nu,
                                     model_nu,
                                     same_gene_model_nu=model_nu,
                                     seqName,
                                     expressionName,
                                     #weightName=character(0),
                                     #Enh_gene_map_name=character(0),
                                     #treat_control_name=character(0),
                                     motif_file_names="motifs.wtmx",
                                     TF_expression_name="TF_exp",
                                     Shared_dir="/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/",
                                     home_dir = character(0),
                                     GEMSTAT_call = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_MultiEnh5/seq2expr",
                                     #par.exp.nu = exp.nu,
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
                                     .ff=character(0),
                                     .po=character(0),
                                     .lower_bound=character(0),
                                     .upper_bound=character(0),
                                     .no_gt_out=character(0),
                                     .softmin_groups=character(0),
                                     .train_weights=character(0),
                                     .control_treat_map=character(0),
                                     .one_beta=character(0),
                                     job_file_name){
  # Function to write hal Jobs for GEMSTAT runs for multi_enhancer version
  # exp.nu is the number if the experiment
  # seqName is the prefix for sequence file names
  # expressionName is the prefix for expression file names
  # weightName is the prefix for weight file names.
  # Enh_gene_map_name is the prefix for Enh_gene_map file names.
  # treat_control_name is the prefix for treat_control file name.
  # This function assumes a certain hirarchy of directories for Inputs:
  # This is shared for all: Shared_dir/
  # sequence:        Shared_dir/Experiment_exp.nu/Input/Sequence/seqName_i.fa
  # gene expression: Shared_dir/Experiment_exp.nu/Input/Gene_expression/expressionName_i.tab
  # motifs:          Shared_dir/Experiment_exp.nu/Input/Motifs/motifs.wtmx
  # TF expression:   Shared_dir/Experiment_exp.nu/Input/TF_Expression/TF_Expression_quantile.tab
  # Shared_dir is the root directory where the data and GEMSTAT program are stored
  # GEMSTAT program to call: e.g. Shared_dir/src_MultiEnh5/seq2expr
  # GEMSTAT_call : address of the gemstat version to be called
  # -a annFile
  # -o modelOption
  # -c coopFile 
  # -i factorInfoFile
  # -r repressionFile 
  # -oo objOption 
  # -mc maxContact 
  # -p (prefix for) parFile 
  # -rt repressionDistThr 
  # -na nAlternations 
  # -ct coopDistThr 
  # -oq if any argument provided for this, there will be a separate qBTM assigned to each CRM otherwise
  #  they will share a qBTM
  # -sigma factorIntSigma
  # -softmin_groups : prefix for the name of enhancer to gene mapping file
  # -control_treat_map: prefix for the name of treatment/control meta file
  # .one_beta : should be used when one beta is used for all enhancers
  # .no_gt_out if used it causes the output to not include ground truth exp values
  # .ff the free fix file
  # .po is a file to write the parameter ouputs to
  # .lower_bound : lower bound file
  # .upper_bound : upper bound file
  # .train_weights : prefix for the name of weights files
  # same_gene_model_nu : is a number only to be used in the case that different models
  #  have the same genes and only one ff, ub, lb, and coop and gene expression file, and 
  #  mapping file and weight file is being generated for all models
  
  
  sequence_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                        "/Inputs/Sequence/",
                        seqName, rep(paste0("_", model_nu),
                                     length(model_nu)),".fa", sep="")
  gene_expression_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                               "/Inputs/Gene_expression/",
                               expressionName, rep(paste0("_", same_gene_model_nu),
                                                   length(model_nu)),
                               ".tab", sep="")
  motif_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                     "/Inputs/Motifs/", motif_file_names, sep="")
  TF_expression_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                             "/Inputs/TF_Expression/",TF_expression_name,".tab", sep="")
  output_file_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                           "/Outputs/", "Experiment_", exp.nu, rep(paste0("_", model_nu),
                                                                   length(model_nu)), ".out", sep="")
  log_file_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                        "/Outputs/", "Experiment_", exp.nu, rep(paste0("_", model_nu),
                                                                length(model_nu)),".log", sep="")
  if(length(.train_weights) > 0){
    weight_file_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                             "/Inputs/Weights/",
                             .train_weights, rep(paste0("_", same_gene_model_nu), length(model_nu)),
                             ".weights", sep="")
  }else{
    weight_file_dir <- character(0)
  }
  if(length(.softmin_groups) > 0){
    Enh_gene_map_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                              "/Inputs/Enh_Gene_Map/",
                              .softmin_groups, rep(paste0("_", same_gene_model_nu),
                                                   length(model_nu)),
                              ".txt", sep="")
  }else{
    Enh_gene_map_dir  <- character(0)
  }
  if(length(.control_treat_map) > 0){
    control_treat_map_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                              "/Inputs/Treat_Control_Map/",
                              .control_treat_map, ".txt", sep="")
  }else{
    control_treat_map_dir <- character(0)
  }
  if(length(.p) > 0){
    parFile_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                         "/Inputs/Parameters/",
                         .p, rep(paste0("_", model_nu),
                                 length(model_nu)), ".par", sep="")
  }else{
    parFile_dir <- character(0)
  }
  if(length(.lower_bound) > 0){
    lb_file_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                         "/Inputs/bounds/",
                         .lower_bound, rep(paste0("_", same_gene_model_nu),
                                           length(model_nu)), ".par", sep="")
  }else{
    lb_file_dir <- character(0)
  }
  if(length(.upper_bound) > 0){
    ub_file_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                         "/Inputs/bounds/",
                         .upper_bound, rep(paste0("_", same_gene_model_nu),
                                           length(model_nu)), ".par", sep="")
  }else{
    ub_file_dir <- character(0)
  }
  if(length(.ff) > 0){
    ff_file_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                         "/Inputs/free_fix/",
                         .ff, rep(paste0("_", same_gene_model_nu),
                                  length(model_nu)), ".par", sep="")
  }else{
    ff_file_dir <- character(0)
  }
  if(length(.c) > 0){
    coop_file_dir <- paste(Shared_dir, "Experiment_", exp.nu,
                           "/Inputs/Coop/",
                           .c, rep(paste0("_", same_gene_model_nu),
                                   length(model_nu)), ".par", sep="")
  }else{
    coop_file_dir <- character(0)
  }
  
  #writing all
  cat(c(GEMSTAT_call, "-s", sequence_dir,
        "-e", gene_expression_dir,
        "-m", motif_dir,
        "-f", TF_expression_dir,
        "-fo", output_file_dir,
        rep("-a", length(.a)), .a,
        rep("-o", length(.o)), .o,
        rep("-c", (length(.c) > 0)), coop_file_dir,
        rep("-i", length(.i)), .i,
        rep("-r", length(.r)), .r,
        rep("-oo", length(.oo)), .oo,
        rep("-mc", length(.mc)), .mc,
        rep("-p", length(.p)), parFile_dir,
        rep("-rt", length(.rt)), .rt,
        rep("-na", length(.na)), .na,
        rep("-ct", length(.ct)), .ct,
        rep("-oq", length(.oq)), .oq,
        rep("-onebeta", length(.one_beta)), .one_beta,
        rep("-sigma", length(.sigma)), .sigma,
        rep("-softmin_groups", length(.softmin_groups)), Enh_gene_map_dir,
        rep("-control_treat_map", length(.control_treat_map)), control_treat_map_dir,
        rep("-train_weights", length(.train_weights)), weight_file_dir,
        rep("-lower_bound", length(.lower_bound)), lb_file_dir,
        rep("-upper_bound", length(.upper_bound)), ub_file_dir,
        rep("-ff", length(.ff)), ff_file_dir,
        " >> ", log_file_dir, "\n"),
      file=job_file_name, sep=" ", append=T)
}
########################################################################################################################
########################################################################################################################
# Example
Hal_job_writer_Multi_enh(exp.nu=1,
                         model_nu=integer(0),
                         seqName="alaki",
                         expressionName="alakitar",
                         motif_file_names="motifs.wtmx",
                         Shared_dir="/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/",
                         home_dir = character(0),
                         GEMSTAT_call = "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_MultiEnh5/seq2expr",
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
                         .ff=character(0),
                         .po=character(0),
                         .lower_bound=character(0),
                         .upper_bound=character(0),
                         .no_gt_out=character(0),
                         .softmin_groups=character(0),
                         .train_weights=character(0),
                         .control_treat_map=character(0),
                         job_file_name="jj")
########################################################################################################################
########################################################################################################################

GEMSTAT_init_BOlinear <- function(.dir=getwd(),
                                  .exp.nu,
                                  TF_names,
                                  model_evaluations,
                                  model_gene_ind=c(1:length(model_evaluations[[1]][[4]])),
                                  model_TF_ind = c(1:length(TF_names)),
                                  .compare_to_expectation=T,
                                  model_parameters,
                                  nu_enh_per_gene,
                                  TF_KD_evaluation,
                                  quantile_nu,
                                  enhancer_GRang,
                                  enhancer_Seq,
                                  real_exp_mat,
                                  TF_expression_mat,
                                  motif_list,
                                  annotation_Filtered_list=numeric(0),
                                  annotation_thresh=numeric(0),
                                  annotation_range=numeric(0),
                                  initial_bind_w=numeric(0),
                                  bind_w_range=numeric(0),
                                  initial_alpha=numeric(0),
                                  alpha_range=numeric(0),
                                  coop_tf_mat=numeric(0),
                                  initial_coop_weight=numeric(0),
                                  coop_weight_range=numeric(0),
                                  coop_type=character(0),
                                  coop_dist=integer(0),
                                  coop_orientation=integer(0),
                                  .one_qbtm_per_enh=T,
                                  initial_qBTM=numeric(0),
                                  qBTMrange=numeric(0),
                                  .one_beta_per_enh=F,
                                  initial_pi_beta=numeric(0), 
                                  pi_beta_upper=numeric(0),
                                  pi_beta_lower=numeric(0),
                                  annotation_thresh_ff=integer(0),
                                  initial_bind_w_ff=integer(0),
                                  initial_alpha_ff=integer(0),
                                  initial_coop_weight_ff=integer(0),
                                  initial_qBTM_ff=integer(0),
                                  initial_pi_beta_ff=integer(0),
                                  nu_enhacners,
                                  .ensemble_mode = F,
                                  nu_samples = 0,
                                  annotation_thresh_ff_ens = annotation_thresh_ff,
                                  initial_coop_weight_ff_ens = initial_coop_weight_ff,
                                  initial_qBTM_ff_ens = initial_qBTM_ff,
                                  na=5,
                                  .job_file_name="job",
                                  .GEMSTAT_call,
                                  job_shared_dir=.dir,
                                  TF_role_dec=F,
                                  create_bounds=T,
                                  create_coop=T,
                                  create_enh_gene_map=T,
                                  create_ff=T,
                                  create_expmat=T, 
                                  create_motifs=T,
                                  create_annotation=T,
                                  create_params=T,
                                  create_seq=T,
                                  create_TFexp=T,
                                  create_TFinfo=T,
                                  create_treatContMap=T,
                                  create_weights=T,
                                  create_SubmitFile=T,
                                  create_jobFile=T,
                                  OnlyTFexpModif=F, 
                                  label_to_kill=integer(0),
                                  all_models_same_genes = F){
  # takes parameters and evaluations of a linear model and initializes PARAMETER, enhacner sequence,
  #  and expression files.
  # creates one or several initial parameter sets seeded in each linear model
  # gets the number of enhancers per model
  # read the explanation for the parameters not ecplained here in the 
  #  "par_ff_lb_ub_coop_Writer_multi_enh" function.
  # nu_enh_per_gene : max number of enhancers assigned to each gene. can be one integer for all
  #  enhancers, or an integer vector with length equal to number of genes. [UPDATE: FOR NOW JUST A
  #  SINGLE INTEGER WHICH INDICATES THE MAX NUMBER OF ENHANCERS PER GENE]
  # model_evaluations : is a list where each entry is the ouput of "objFuncEval" function for one model
  # model_parameters : is a matrix with nrow equal to the number of models and ncol equal to number of
  #  parameters
  # TF_KD_evaluation : is the output of "TF_KD_evaluator" for all TFs
  # quantile_nu : is the number of quantiles used for mappings between linear and GEMSTAT
  #  parameters. negative and positive values will be treated separately (i.e. quantile_nu quantiles
  #  for positives and quantile_nu quantiles for negative values)
  # annotation_Filtered_list : is the third entry of the output of count_site_from_annotation_combined 
  #  function. is a list that contains one entry per model, each of those entries contain one entry pergene.
  # na : .na number of alteration in GEMSTAT call
  # job_shared_dir : if you want the directory used for writing jobs to be different from where the
  #  files are stored, put it here
  # .TF_role_decided : if True, the roles of the TFs are decided beforehand, a factorInfo file will be
  #   created and upper lower bounds, initial values and freefix will be adjusted accordingly. If False (default)
  #   no prior role is decided, factorInfo file will not be created. only the activator alpha range will
  #   be used for both repression and activation. 
  # .compare_to_expectation : if True only considers enhancers that have done better than random expectation. (if True some genes might be removed)
  # model_TF_ind : index of the TFs in the linear model, to be used in gemstat model. eg. if 
  #  we only want the first three TFs of linear model in the GEMSTAT model c(1, 2, 3)
  # create_* : if True and all other conditions are satisfied * will be produced. if False it won't be produced.
  # OnlyTFexpModif : if True it will only create a new TF experssion matrix and everything else has to be copied from somewhere else
  # all_models_same_genes : logical. If true means all models are using th same gene set. If False different
  #  models will have their own ff, expression matrix, ub and lb and so on ..
  #create the default ranges unless specified
  ############################
  # the following parameters:
  
  # initial_qBTM
  # qBTMrange
  # initial_pi_beta
  # pi_beta_upper
  # pi_beta_lower
  # initial_qBTM_ff
  # initial_pi_beta_ff
  # initial_qBTM_ff_ens
  
  # can not be specified as before since we are not aware of the number of enhancers beforehand. hence
  #  they take either numeric(0), or just one number for the parameter which will be used for all
  #  enhancers. (in case of ranges they take two numbers)
  ############################
  stopifnot(!(all_models_same_genes & .compare_to_expectation))
  
  if(OnlyTFexpModif){
    if(dir.exists(paste0(.dir,"/Experiment_", .exp.nu))){
      print("Directory for this experiment number already exists")
      return(0)
    }
    dir.create(paste0(.dir,"/Experiment_", .exp.nu))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Outputs"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Gene_expression"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Sequence"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Motifs"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/TF_Expression"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Weights"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Enh_Gene_Map"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Treat_Control_Map"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/bounds"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/free_fix"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Coop"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/TF_info"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Parameters"))
    dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Filtered_Annotations"))
    if(create_TFexp){
      setwd(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/TF_Expression"))
      ExpressionWriter(expressionMat = TF_expression_mat,
                       output.File.Name = "TF_exp")
    }
    return(0)
  }
  if(.one_qbtm_per_enh){
    my_qBTMrange1 <- matrix(rep(qBTMrange, nu_enhacners), nrow= nu_enhacners, byrow = T)
  }else{
    my_qBTMrange1 <- qBTMrange
  }
  my_default_par <- par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=TF_names,
                                                       .annotation_thresh=annotation_thresh,
                                                       .annotation_range=annotation_range,
                                                       .initial_bind_w=initial_bind_w,
                                                       .bind_w_range=bind_w_range,
                                                       .initial_alpha=initial_alpha,
                                                       .alpha_range=alpha_range,
                                                       .coop_tf_mat=coop_tf_mat,
                                                       .initial_coop_weight=initial_coop_weight,
                                                       .coop_weight_range=coop_weight_range,
                                                       .coop_type=coop_type,
                                                       .coop_dist=coop_dist,
                                                       .coop_orientation=coop_orientation,
                                                       one_qbtm_per_enh=.one_qbtm_per_enh,
                                                       .initial_qBTM=initial_qBTM,
                                                       .qBTMrange=my_qBTMrange1,
                                                       one_beta_per_enh=.one_beta_per_enh,
                                                       .initial_pi_beta=initial_pi_beta, 
                                                       .pi_beta_upper=pi_beta_upper,
                                                       .pi_beta_lower=pi_beta_lower,
                                                       .annotation_thresh_ff=annotation_thresh_ff,
                                                       .initial_bind_w_ff=initial_bind_w_ff,
                                                       .initial_alpha_ff=initial_alpha_ff,
                                                       .initial_coop_weight_ff=initial_coop_weight_ff,
                                                       .initial_qBTM_ff=initial_qBTM_ff,
                                                       .initial_pi_beta_ff=initial_pi_beta_ff,
                                                       .filename_start=character(0),
                                                       .filename_ff=character(0),
                                                       .filename_upper=character(0),
                                                       .filename_lower=character(0),
                                                       .file_name_coop=character(0),
                                                       .nu_enhacners=nu_enhacners,
                                                       ensemble_mode = .ensemble_mode,
                                                       .nu_samples = nu_samples,
                                                       .annotation_thresh_ff_ens = annotation_thresh_ff_ens,
                                                       .initial_coop_weight_ff_ens = initial_coop_weight_ff_ens,
                                                       .initial_qBTM_ff_ens = initial_qBTM_ff_ens,
                                                       Write_to_file = F,
                                                       TF_role_decided = TF_role_dec)
  
  # First converting the parameters to GEMSTAT scale: using quantiles, positive and negative separately
  # parameters to convert: all TF parameters to alpha (activator or repressor) parameters, intercept parameter to qBTM parameter
  # creating the quantiles for positive and negative values separately for all parameters together
  # .dir : is the shared directory where the input and output data is going to be stored
  Linear_model_quantiles <- matrix(nrow=2, ncol=quantile_nu+1)
  rownames(Linear_model_quantiles) <- c('pos', 'neg')
  stopifnot(length(model_TF_ind) == length(TF_names))
  my_model_parameters <- cbind(model_parameters[, model_TF_ind], model_parameters[, ncol(model_parameters)])
  TF_parameters <- my_model_parameters[ ,1:(ncol(my_model_parameters) - 1)]
  Linear_model_quantiles[1, ] <- quantile(TF_parameters[TF_parameters >= 0],
                                          seq(0, 1, length.out = quantile_nu+1))
  Linear_model_quantiles[2, ] <- quantile(TF_parameters[TF_parameters <  0],
                                          seq(0, 1, length.out = quantile_nu+1))
  
  full_transformed_mat <- matrix(nrow = nrow(my_model_parameters), ncol = ncol(my_model_parameters))
  #quantiling the TF parameters
  associated_value_pos <- seq(max(min(my_default_par$initial_alpha_range[, 1]), 1),
                              max(my_default_par$initial_alpha_range[, 2]),
                              length.out = quantile_nu + 2)
  # to avoid hitting the lower and upper bounds
  associated_value_pos <- associated_value_pos[2:(length(associated_value_pos) - 1)] 
  associated_value_neg <- seq(min(my_default_par$initial_alpha_range[, 1]),
                              min(1, max(my_default_par$initial_alpha_range[, 2])),
                              length.out = quantile_nu + 2)
  # to avoid hitting the lower and upper bounds
  associated_value_neg <- associated_value_neg[2: (length(associated_value_neg) - 1)] 
  # create the TF role matrix
  TF_role_mat <- matrix(nrow = nrow(my_model_parameters), ncol = length(TF_names))
  colnames(TF_role_mat) <- TF_names
  
  for(par_nu in 1:(ncol(my_model_parameters) - 1)){
    cur_pos_par <- (my_model_parameters[, par_nu] >= 0)
    cur_neg_par <- (my_model_parameters[, par_nu] <= 0)
    TF_role_mat[, par_nu] <- as.integer((my_model_parameters[, par_nu] >= 0))
    transformed_par <- numeric(nrow(my_model_parameters))

    for(cur_qu in 1:quantile_nu){
      transformed_par[(cur_pos_par == 1 &
                         (my_model_parameters[, par_nu] >= Linear_model_quantiles[1, cur_qu]) &
                         (my_model_parameters[, par_nu] < Linear_model_quantiles[1, cur_qu+1]))] <- min(associated_value_pos[cur_qu],
                                                                                                     (my_default_par$initial_alpha_range[par_nu, 2] - 10e-7))
      transformed_par[(cur_neg_par == 1 &
                         (my_model_parameters[, par_nu] >= Linear_model_quantiles[2, cur_qu]) &
                         (my_model_parameters[, par_nu] < Linear_model_quantiles[2, cur_qu+1]))] <- max(associated_value_neg[cur_qu],
                                                                                                     (my_default_par$initial_alpha_range[par_nu, 1]+ 10e-7))
      if(cur_qu == quantile_nu){
        transformed_par[(cur_pos_par == 1 &
                           (my_model_parameters[, par_nu] >= Linear_model_quantiles[1, cur_qu]) &
                           (my_model_parameters[, par_nu] <= Linear_model_quantiles[1, cur_qu+1]))] <- min(associated_value_pos[cur_qu],
                                                                                                       (my_default_par$initial_alpha_range[par_nu, 2] - 10e-7))
        transformed_par[(cur_neg_par == 1 &
                           (my_model_parameters[, par_nu] >= Linear_model_quantiles[2, cur_qu]) &
                           (my_model_parameters[, par_nu] <= Linear_model_quantiles[2, cur_qu+1]))] <- max(associated_value_neg[cur_qu],
                                                                                                       (my_default_par$initial_alpha_range[par_nu, 1]+ 10e-7))
      }
    }
    full_transformed_mat[, par_nu] <- transformed_par
  }
  #quantiling the qBTM parameter
  qBTM_parameters <- my_model_parameters[ , ncol(my_model_parameters)]
  qBTM_transformed <- numeric(length(qBTM_parameters))
  Linear_model_quantiles_qBTM <- quantile(qBTM_parameters, seq(0, 1, length.out = quantile_nu+1))
  associated_value_qBTM <- seq(max(0, min(my_default_par$initial_qBTM_range)),
                               min(1,  max(my_default_par$initial_qBTM_range)),
                               length.out = quantile_nu + 2)

  associated_value_qBTM <- associated_value_qBTM[2:(length(associated_value_qBTM) - 1)]
  for(cur_qu in 1:quantile_nu){
    qBTM_transformed[(qBTM_parameters >= Linear_model_quantiles_qBTM[cur_qu]) &
                       (qBTM_parameters < Linear_model_quantiles_qBTM[cur_qu+1])] <- associated_value_qBTM[cur_qu]
    if(cur_qu == quantile_nu){
      qBTM_transformed[(qBTM_parameters >= Linear_model_quantiles_qBTM[cur_qu]) &
                         (qBTM_parameters <= Linear_model_quantiles_qBTM[cur_qu+1])] <- associated_value_qBTM[cur_qu]
    }
  }
  full_transformed_mat[, ncol(full_transformed_mat)] <- qBTM_transformed
  # print("full_transformed_mat")
  # print(full_transformed_mat)
  # now I can start creating files for GEMSTAT models initiated by a line in full_transformed_mat
  # create the directories
  prev_dir <- getwd()
  on.exit(setwd(prev_dir))
  if(dir.exists(paste0(.dir,"/Experiment_", .exp.nu))){
    print("Directory for this experiment number already exists")
    return(0)
  }
  dir.create(paste0(.dir,"/Experiment_", .exp.nu))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Outputs"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Gene_expression"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Sequence"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Motifs"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/TF_Expression"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Weights"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Enh_Gene_Map"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Treat_Control_Map"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/bounds"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/free_fix"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Coop"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/TF_info"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Parameters"))
  dir.create(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Filtered_Annotations"))
  # DEPENDENT ON THE SEED:
  All_model_loss_per_gene_per_enh <- lapply(model_evaluations, "[[", 4)
  Chosen_model_loss_per_gene_per_enh <- list()
  for(i in 1:length(All_model_loss_per_gene_per_enh)){
    Chosen_model_loss_per_gene_per_enh[[i]] <- All_model_loss_per_gene_per_enh[[i]][model_gene_ind]
  }
  for(cur_model in 1:nrow(my_model_parameters)){
    # Write enhancers
    setwd(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Sequence"))
    Enhacner_Sequence_written <- Enhancer_writer(.loss_per_model_per_gene_per_enh = Chosen_model_loss_per_gene_per_enh,
                                                 .enhancer_Granges = enhancer_GRang,
                                                 .enhancer_sequence = enhancer_Seq,
                                                 .model_index=cur_model,
                                                 .real_exp_mat = real_exp_mat,
                                                 max_per_gene=nu_enh_per_gene,
                                                 compare_to_expectation=.compare_to_expectation,
                                                 sd_from_expected=0,
                                                 seq_file_name=paste0("sequences_", cur_model),
                                                 write_to_file = create_seq)
    # remove the genes that were not included in the sequence writing part
    new_nu_enhancers <- sum(unlist(lapply(Enhacner_Sequence_written$sequenceList, length)))
    print("old number of enhancers: ")
    print(nu_enhacners)
    print("new number of enhancers: ")
    print(new_nu_enhancers)
    print("Enhacner_Sequence_written$removed_genes")
    print(Enhacner_Sequence_written$removed_genes)
    if(length(Enhacner_Sequence_written$removed_genes) > 0){
      new_real_exp_mat <- real_exp_mat[-Enhacner_Sequence_written$removed_genes, ]
    }else{
      new_real_exp_mat <- real_exp_mat
    }
    # Write expression
    setwd(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Gene_expression"))
    aa_map_c <- unlist(lapply(Enhacner_Sequence_written$sequenceList, length))
    aa_map = integer(0)
    for(i in 0:(length(aa_map_c)-1)){
      aa_map <- c(aa_map, rep(i, aa_map_c[i+1]))
    }

    if(create_expmat){
      if(all_models_same_genes & nu_enh_per_gene == 1){
        if(cur_model == 1){
          ExpressionWriter(expressionMat = new_real_exp_mat,
                           from_log_fold_change = T,
                           output.File.Name = paste0("expression_", cur_model),
                           multi_enhancer_map = aa_map)
        }
      }else{
        ExpressionWriter(expressionMat = new_real_exp_mat,
                         from_log_fold_change = T,
                         output.File.Name = paste0("expression_", cur_model),
                         multi_enhancer_map = aa_map)
      }
    }

    # Write TF info file
    if(TF_role_dec & create_TFinfo){
      setwd(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/TF_info"))
      FactorInfoWriter(.TF_names = TF_names,
                       role_vec = TF_role_mat[cur_model, ],
                       file_name = paste0("factor_info_", cur_model))
    }
    # Write weights
    setwd(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Weights"))
    if(create_weights){
      if(all_models_same_genes & nu_enh_per_gene == 1){
        if(cur_model == 1){
          GEMSTAT_weight_writer(expmat =new_real_exp_mat ,
                                file_name =  paste0("weights_", cur_model),
                                multi_enhancer_map = aa_map,
                                fold_change_to_real = T,
                                to_kill = label_to_kill)
        }
      }else{
        GEMSTAT_weight_writer(expmat =new_real_exp_mat ,
                              file_name =  paste0("weights_", cur_model),
                              multi_enhancer_map = aa_map,
                              fold_change_to_real = T,
                              to_kill = label_to_kill)
      }
    }

    # Write enhancer_gene_map
    #write mapping file
    setwd(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Enh_Gene_Map"))
    if(create_enh_gene_map){
      if(all_models_same_genes & nu_enh_per_gene == 1){
        if(cur_model == 1){
          for(i in 1:length(aa_map)){
            cat(c(aa_map[i]), file = paste0("geneEnhMap_", cur_model,".txt"), sep = "\n", append = T)
          }
        }
      }else{
        for(i in 1:length(aa_map)){
          cat(c(aa_map[i]), file = paste0("geneEnhMap_", cur_model,".txt"), sep = "\n", append = T)
        }
      }
    }

    
    # Write the ub, lb, coop, initial values
    setwd(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs"))
    if(.one_qbtm_per_enh){
      my_init_qBTM <- rep(full_transformed_mat[cur_model, ncol(full_transformed_mat)],
                          new_nu_enhancers)
    }else{
      my_init_qBTM <- full_transformed_mat[cur_model, ncol(full_transformed_mat)]
    }
    
    # check size of initial_qBTM_ff
    if(length(initial_qBTM_ff) > 0){
      print("setting new initial_qBTM_ff_ens")
      if(length(initial_qBTM_ff) != 1){
        print(paste("WARNING initial_qBTM_ff has length", length(initial_qBTM_ff) ,
                    "Which is not expected"))
      }else{
        if(.one_qbtm_per_enh){
          new_initial_qBTM_ff <- rep(initial_qBTM_ff, new_nu_enhancers)
        }else{
          new_initial_qBTM_ff <- initial_qBTM_ff
        }
      }
    }else{
      new_initial_qBTM_ff <-numeric(0)
    }
    print("new_initial_qBTM_ff")
    print(new_initial_qBTM_ff)
    # check size of initial_qBTM_ff_ens
    if(length(initial_qBTM_ff_ens) > 0){
      print("setting new initial_qBTM_ff_ens")
      if(length(initial_qBTM_ff_ens) != 1){
        print(paste("WARNING initial_qBTM_ff_ens has length", length(initial_qBTM_ff_ens) ,
                    "Which is not expected"))
      }else{
        if(.one_qbtm_per_enh){
          new_initial_qBTM_ff_ens <- rep(initial_qBTM_ff_ens, new_nu_enhancers)
        }else{
          new_initial_qBTM_ff_ens <- initial_qBTM_ff_ens
        }
      }
    }else{
      new_initial_qBTM_ff_ens <-numeric(0)
    }
    # check size of qBTMrange
    if(length(qBTMrange) > 0){
      if(length(qBTMrange) != 2){
        print(paste("WARNING qBTMrange has length", length(qBTMrange) , "Which is not expected"))
      }else{
        if(.one_qbtm_per_enh){
          new_qBTMrange <- matrix(rep(matrix(qBTMrange, nrow = 1), new_nu_enhancers),
                                  nrow=new_nu_enhancers, byrow = T)
        }else{
          new_qBTMrange <- qBTMrange
        }
      }
    }else{
      new_qBTMrange <- numeric(0)
    }

    if(TF_role_dec){
      my_TF_role_mat <- TF_role_mat[cur_model, ]
    }else{
      my_TF_role_mat <- integer(0)
    }

    if(all_models_same_genes & nu_enh_per_gene == 1){ # To prevent writing redundant lb, ub, ff and coop files
      if(cur_model == 1){
        par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=TF_names,
                                           .annotation_thresh=annotation_thresh,
                                           .annotation_range=annotation_range,
                                           .initial_bind_w=initial_bind_w,
                                           .bind_w_range=bind_w_range,
                                           .initial_alpha=full_transformed_mat[cur_model, (1:(ncol(full_transformed_mat) - 1))],
                                           .alpha_range=alpha_range,
                                           .coop_tf_mat=coop_tf_mat,
                                           .initial_coop_weight=initial_coop_weight,
                                           .coop_weight_range=coop_weight_range,
                                           .coop_type=coop_type,
                                           .coop_dist=coop_dist,
                                           .coop_orientation=coop_orientation,
                                           one_qbtm_per_enh=.one_qbtm_per_enh,
                                           .initial_qBTM=my_init_qBTM,
                                           .qBTMrange=new_qBTMrange,
                                           one_beta_per_enh=.one_beta_per_enh,
                                           .initial_pi_beta=initial_pi_beta, 
                                           .pi_beta_upper=pi_beta_upper,
                                           .pi_beta_lower=pi_beta_lower,
                                           .annotation_thresh_ff=annotation_thresh_ff,
                                           .initial_bind_w_ff=initial_bind_w_ff,
                                           .initial_alpha_ff=initial_alpha_ff,
                                           .initial_coop_weight_ff=initial_coop_weight_ff,
                                           .initial_qBTM_ff=new_initial_qBTM_ff,
                                           .initial_pi_beta_ff=initial_pi_beta_ff,
                                           .filename_start=paste0("start_", cur_model,".par"),
                                           .filename_ff=paste0("ff_", cur_model,".par"),
                                           .filename_upper=paste0("ub_", cur_model,".par"),
                                           .filename_lower=paste0("lb_", cur_model,".par"),
                                           .file_name_coop=paste0("coop_", cur_model,".par"),
                                           .nu_enhacners=new_nu_enhancers,
                                           ensemble_mode = F,
                                           .nu_samples = nu_samples,
                                           .annotation_thresh_ff_ens = annotation_thresh_ff_ens,
                                           .initial_coop_weight_ff_ens = initial_coop_weight_ff_ens,
                                           .initial_qBTM_ff_ens = new_initial_qBTM_ff_ens,
                                           Write_to_file = T,
                                           create_folders = T,
                                           TF_role_decided = TF_role_dec,
                                           factor_role_vector = my_TF_role_mat, 
                                           .create_bounds = create_bounds ,
                                           .create_ff = create_ff,
                                           .create_par = create_params,
                                           .create_coop = create_coop)
      }else{
        par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=TF_names,
                                           .annotation_thresh=annotation_thresh,
                                           .annotation_range=annotation_range,
                                           .initial_bind_w=initial_bind_w,
                                           .bind_w_range=bind_w_range,
                                           .initial_alpha=full_transformed_mat[cur_model, (1:(ncol(full_transformed_mat) - 1))],
                                           .alpha_range=alpha_range,
                                           .coop_tf_mat=coop_tf_mat,
                                           .initial_coop_weight=initial_coop_weight,
                                           .coop_weight_range=coop_weight_range,
                                           .coop_type=coop_type,
                                           .coop_dist=coop_dist,
                                           .coop_orientation=coop_orientation,
                                           one_qbtm_per_enh=.one_qbtm_per_enh,
                                           .initial_qBTM=my_init_qBTM,
                                           .qBTMrange=new_qBTMrange,
                                           one_beta_per_enh=.one_beta_per_enh,
                                           .initial_pi_beta=initial_pi_beta, 
                                           .pi_beta_upper=pi_beta_upper,
                                           .pi_beta_lower=pi_beta_lower,
                                           .annotation_thresh_ff=annotation_thresh_ff,
                                           .initial_bind_w_ff=initial_bind_w_ff,
                                           .initial_alpha_ff=initial_alpha_ff,
                                           .initial_coop_weight_ff=initial_coop_weight_ff,
                                           .initial_qBTM_ff=new_initial_qBTM_ff,
                                           .initial_pi_beta_ff=initial_pi_beta_ff,
                                           .filename_start=paste0("start_", cur_model,".par"),
                                           .filename_ff=paste0("ff_", cur_model,".par"),
                                           .filename_upper=paste0("ub_", cur_model,".par"),
                                           .filename_lower=paste0("lb_", cur_model,".par"),
                                           .file_name_coop=paste0("coop_", cur_model,".par"),
                                           .nu_enhacners=new_nu_enhancers,
                                           ensemble_mode = F,
                                           .nu_samples = nu_samples,
                                           .annotation_thresh_ff_ens = annotation_thresh_ff_ens,
                                           .initial_coop_weight_ff_ens = initial_coop_weight_ff_ens,
                                           .initial_qBTM_ff_ens = new_initial_qBTM_ff_ens,
                                           Write_to_file = T,
                                           create_folders = T,
                                           TF_role_decided = TF_role_dec,
                                           factor_role_vector = my_TF_role_mat, 
                                           .create_bounds = F ,
                                           .create_ff = F,
                                           .create_par = create_params,
                                           .create_coop = F)
      }
    }else{
      par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=TF_names,
                                         .annotation_thresh=annotation_thresh,
                                         .annotation_range=annotation_range,
                                         .initial_bind_w=initial_bind_w,
                                         .bind_w_range=bind_w_range,
                                         .initial_alpha=full_transformed_mat[cur_model, (1:(ncol(full_transformed_mat) - 1))],
                                         .alpha_range=alpha_range,
                                         .coop_tf_mat=coop_tf_mat,
                                         .initial_coop_weight=initial_coop_weight,
                                         .coop_weight_range=coop_weight_range,
                                         .coop_type=coop_type,
                                         .coop_dist=coop_dist,
                                         .coop_orientation=coop_orientation,
                                         one_qbtm_per_enh=.one_qbtm_per_enh,
                                         .initial_qBTM=my_init_qBTM,
                                         .qBTMrange=new_qBTMrange,
                                         one_beta_per_enh=.one_beta_per_enh,
                                         .initial_pi_beta=initial_pi_beta, 
                                         .pi_beta_upper=pi_beta_upper,
                                         .pi_beta_lower=pi_beta_lower,
                                         .annotation_thresh_ff=annotation_thresh_ff,
                                         .initial_bind_w_ff=initial_bind_w_ff,
                                         .initial_alpha_ff=initial_alpha_ff,
                                         .initial_coop_weight_ff=initial_coop_weight_ff,
                                         .initial_qBTM_ff=new_initial_qBTM_ff,
                                         .initial_pi_beta_ff=initial_pi_beta_ff,
                                         .filename_start=paste0("start_", cur_model,".par"),
                                         .filename_ff=paste0("ff_", cur_model,".par"),
                                         .filename_upper=paste0("ub_", cur_model,".par"),
                                         .filename_lower=paste0("lb_", cur_model,".par"),
                                         .file_name_coop=paste0("coop_", cur_model,".par"),
                                         .nu_enhacners=new_nu_enhancers,
                                         ensemble_mode = F,
                                         .nu_samples = nu_samples,
                                         .annotation_thresh_ff_ens = annotation_thresh_ff_ens,
                                         .initial_coop_weight_ff_ens = initial_coop_weight_ff_ens,
                                         .initial_qBTM_ff_ens = new_initial_qBTM_ff_ens,
                                         Write_to_file = T,
                                         create_folders = T,
                                         TF_role_decided = TF_role_dec,
                                         factor_role_vector = my_TF_role_mat, 
                                         .create_bounds = create_bounds,
                                         .create_ff = create_ff,
                                         .create_par = create_params,
                                         .create_coop = create_coop)
    }
  }
  # # INDEPENDENT OF THE SEED:
  # Write TF expression mat
  if(create_TFexp){
    setwd(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/TF_Expression"))
    ExpressionWriter(expressionMat = TF_expression_mat,
                     output.File.Name = "TF_exp")
  }
  
  if(create_motifs){
    # Write motifs
    setwd(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Motifs"))
    MotifWriter(motif.List = motif_list, output.File.Name = "motifs")
  }
  
  if(length(annotation_Filtered_list) >0){
    if(create_annotation){
      setwd(paste0(.dir,"/Experiment_", .exp.nu,
                   "/Inputs", "/Filtered_Annotations"))
      annotation_Writer(annotation_Filtered_list)
    }
  }
  # Write control-treatment map
  setwd(paste0(.dir,"/Experiment_", .exp.nu, "/Inputs", "/Treat_Control_Map"))
  aa_cnt_trt <- integer(0)
  for(i in 0:(ncol(real_exp_mat)-1)){
    aa_cnt_trt <- c(aa_cnt_trt, rep(i, 2))
  }
  if(create_treatContMap){
    for(i in 1:length(aa_cnt_trt)){
      cat(c(aa_cnt_trt[i]), file = "ctrl_treat_map.txt", sep = "\n", append = T)
    }
  }

  # write jobs and submit files
  setwd(paste0(.dir,"/Experiment_", .exp.nu))
  my_coop <- character(0)
  if(length(coop_tf_mat) > 0){
    my_coop <- "coop"
  }
  if(length(.job_file_name) == 1){
    .job_file_name <- rep(.job_file_name, nrow(my_model_parameters))
  }
  # TODO: fix the jobs so that if refers to the first file if the rest is redundant
  for(cur_model in 1:nrow(my_model_parameters)){
    if(.one_beta_per_enh){
      my_one_beta <- character(0)
    }else{
      my_one_beta <- 1
    }
    
    if(.one_qbtm_per_enh){
      my_oq <- 1
    }else{
      my_oq <- character(0)
    }
    if(TF_role_dec){
      my_facInfo <- paste0(job_shared_dir,"/Experiment_",
                           .exp.nu, "/Inputs", "/TF_info",
                           "/","factor_info_", cur_model)
    }else{
      my_facInfo <- character(0)
    }
    if(length(annotation_Filtered_list) > 0){
      my_annot <- paste0(job_shared_dir, "/Experiment_",
                         .exp.nu, "/Inputs/", "Filtered_Annotations/annot_",
                         cur_model, ".ann")
    }else{
      my_annot <- character(0)
    }
    if(create_jobFile){
      if(all_models_same_genes & nu_enh_per_gene == 1){
        Hal_job_writer_Multi_enh(exp.nu=.exp.nu,
                                 model_nu=cur_model,
                                 same_gene_model_nu = 1,
                                 seqName="sequences",
                                 expressionName="expression",
                                 motif_file_names="motifs.wtmx",
                                 Shared_dir=paste0(job_shared_dir, "/"),
                                 home_dir = character(0),
                                 GEMSTAT_call = .GEMSTAT_call,
                                 .a=my_annot,
                                 .o="DIRECT",
                                 .c=my_coop,
                                 .i= my_facInfo,
                                 .r=character(0),
                                 .oo=character(0),
                                 .mc=character(0),
                                 .p="start",
                                 .rt=character(0),
                                 .na=na,
                                 .ct=character(0),
                                 .oq=my_oq,
                                 .one_beta = my_one_beta,
                                 .sigma=character(0),
                                 .ff="ff",
                                 .po=character(0),
                                 .lower_bound="lb",
                                 .upper_bound="ub",
                                 .no_gt_out=character(0),
                                 .softmin_groups="geneEnhMap",
                                 .train_weights="weights",
                                 .control_treat_map="ctrl_treat_map",
                                 job_file_name=.job_file_name[cur_model])
      }else{
        Hal_job_writer_Multi_enh(exp.nu=.exp.nu,
                                 model_nu=cur_model,
                                 seqName="sequences",
                                 expressionName="expression",
                                 motif_file_names="motifs.wtmx",
                                 Shared_dir=paste0(job_shared_dir, "/"),
                                 home_dir = character(0),
                                 GEMSTAT_call = .GEMSTAT_call,
                                 .a=my_annot,
                                 .o="DIRECT",
                                 .c=my_coop,
                                 .i= my_facInfo,
                                 .r=character(0),
                                 .oo=character(0),
                                 .mc=character(0),
                                 .p="start",
                                 .rt=character(0),
                                 .na=na,
                                 .ct=character(0),
                                 .oq=my_oq,
                                 .one_beta = my_one_beta,
                                 .sigma=character(0),
                                 .ff="ff",
                                 .po=character(0),
                                 .lower_bound="lb",
                                 .upper_bound="ub",
                                 .no_gt_out=character(0),
                                 .softmin_groups="geneEnhMap",
                                 .train_weights="weights",
                                 .control_treat_map="ctrl_treat_map",
                                 job_file_name=.job_file_name[cur_model])
      }
    }
  } # end of loop over models
  # write script to create submit file
  unique_jobs <- unique(.job_file_name)
  if(create_SubmitFile){
    for(ujfn in 1:length(unique_jobs)){
      line_number <- length(readLines(unique_jobs[ujfn]))
      sys_com <- paste0("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
                        " ", unique_jobs[ujfn], " ",
                        line_number, " tmpjob_", unique_jobs[ujfn]," >",
                        unique_jobs[ujfn], ".submit")
      cat(c("#!/bin/bash", "\n"),
          file=paste0("hal_sub_creator_", unique_jobs[ujfn]), sep="")
      cat(sys_com,
          file=paste0("hal_sub_creator_", unique_jobs[ujfn]), sep="", append = T)
      sys_com <- paste0("chmod +x ", paste0("hal_sub_creator_", unique_jobs[ujfn]))
      system(sys_com)
    }
  }
}
########################################################################################################################
save(list = c("GEMSTAT_init_BOlinear", "Hal_job_writer_Multi_enh"),file = "Input_Constructor.RData")
########################################################################################################################
########################################################################################################################
#example
setwd("Seeded_GEMSTAT_ens/")
aa_gd <- getwd()
aa <- GEMSTAT_init_BOlinear(.exp.nu = 23000,
                      .dir=aa_gd,
                      .GEMSTAT_call = "/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr",
                      TF_names = names(TF.motifs.Shrinked.count)[c(1, 2, 3, 4, 4, 5)],
                      model_evaluations = Sim_Ann_weighted_148_restart_all_models_eval[1:2],
                      model_gene_ind = c(1,2),
                      model_TF_ind = c(1, 2, 3, 4, 4, 5),
                      #model_parameters = cbind(Sim_Ann_weighted_148_restart_parameters[1:2,1:3], Sim_Ann_weighted_148_restart_parameters[1:2,20]),
                      model_parameters = Sim_Ann_weighted_148_restart_parameters[1:2,],
                      nu_enh_per_gene = 1,
                      TF_KD_evaluation = Sim_Ann_weighted_148_restart_TF_KD,
                      quantile_nu = 4,
                      enhancer_GRang = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges[1:2],
                      enhancer_Seq = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq[1:2],
                      real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[1:2,],
                      TF_expression_mat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01[c(1, 2, 3, 4, 4, 5), ],
                      motif_list = TF.motifs.Shrinked.count[1:3],
                      annotation_thresh=numeric(0),
                      annotation_range=numeric(0),
                      initial_bind_w=numeric(0),
                      bind_w_range=numeric(0),
                      initial_alpha=numeric(0),
                      alpha_range=numeric(0),
                      coop_tf_mat=numeric(0),
                      initial_coop_weight=numeric(0),
                      coop_weight_range=numeric(0),
                      coop_type=character(0),
                      coop_dist=integer(0),
                      coop_orientation=integer(0),
                      .one_qbtm_per_enh=T,
                      initial_qBTM=numeric(0),
                      qBTMrange=numeric(0),
                      .one_beta_per_enh=T,
                      initial_pi_beta=numeric(0), 
                      pi_beta_upper=numeric(0),
                      pi_beta_lower=numeric(0),
                      annotation_thresh_ff=integer(0),
                      initial_bind_w_ff=integer(0),
                      initial_alpha_ff=integer(0),
                      initial_coop_weight_ff=integer(0),
                      initial_qBTM_ff=1,
                      initial_pi_beta_ff=integer(0),
                      nu_enhacners = length(unlist(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq[1:3])),
                      .ensemble_mode = F,
                      nu_samples = 0,
                      na = 2,
                      .job_file_name = "newjob",
                      .compare_to_expectation = T,
                      TF_role_dec = T
                      , create_ff = T, create_bounds = F, create_coop = F,
                      create_expmat = F, create_motifs = F, create_params = F,
                      create_seq = F, create_TFexp = F, create_TFinfo = F,
                      create_enh_gene_map = F, create_treatContMap = F, create_weights = F
                      #,annotation_thresh_ff_ens = annotation_thresh_ff,
                      #initial_coop_weight_ff_ens = initial_coop_weight_ff,
                     # initial_qBTM_ff_ens = initial_qBTM_ff
                      )
ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene
my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52
"/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/src_MultiEnh5/seq2expr"
aa <- lapply(Sim_Ann_weighted_148_restart_all_models_eval, "[[", 4)
aabb <- list()
for(i in 1:length(aa)){
  print(i)
  aabb[[i]] <- aa[[i]][c(1:2)]
}
########################################################################################################################
########################################################################################################################
GEMSTAT_input_copier_multienh <- function(lower_bounds=integer(0),
                                          upper_bounds=integer(0),
                                          coop=integer(0), 
                                          annotations=integer(0),
                                          Enh_gene_map=integer(0),
                                          free_fix=integer(0),
                                          Gene_Exp=integer(0), 
                                          Motifs=integer(0),
                                          Params=integer(0), 
                                          Params_from_output=F,
                                          seq=integer(0),
                                          TFexp=integer(0),
                                          TF_info=integer(0),
                                          Treat_cont_map=integer(0), 
                                          wights=integer(0),
                                          jobfile = integer(0),
                                          .prev_na=5, 
                                          .optim=T,
                                          cur_exp, 
                                          shared_dir){
  # each input is the the number of experiment from which the corresponding input should be copied
  # Params_from_output is a logical. if True it takes the optimized parameters of the mentioned experiment instead
  #  of initial parameters.
  # cur_exp is the number of current experiment
  # shared_dir is the directory where all experiments reside.
  # .prev_na the number of alternations in the job file that we want to copy from
  # .optim : True if we want to have optimization in the new job (sets na to 5), False if we want to just evaluate (sets na to 0)
  prev_dir <- getwd()
  on.exit(setwd(prev_dir))
  
  # copy lower bounds
  if(length(lower_bounds)==1){
    print(paste0("copying lower_bounds from inputs of experiment ", lower_bounds))
    cur_dir <- paste0(shared_dir, "/", "Experiment_", lower_bounds, "/Inputs/bounds/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/bounds/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    setwd(cur_dir)
    list_of_files <- list.files(paste0(shared_dir, "/", "Experiment_", lower_bounds, "/Inputs/bounds/" ), pattern = "lb_*")
    file.copy(from = list_of_files, to = dest_dir)
  }
  # copy upper bounds
  if(length(upper_bounds)==1){
    print(paste0("copying upper_bounds from inputs of experiment ", upper_bounds))
    cur_dir <- paste0(shared_dir, "/", "Experiment_", upper_bounds, "/Inputs/bounds/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/bounds/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    setwd(cur_dir)
    list_of_files <- list.files(paste0(shared_dir, "/", "Experiment_", upper_bounds, "/Inputs/bounds/" ), pattern = "ub_*")
    file.copy(from = list_of_files, to = dest_dir)
  }
  # copy coops
  if(length(coop)==1){
    cur_dir <- paste0(shared_dir, "/", "Experiment_", coop, "/Inputs/Coop/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/Coop/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    
    if(length(list.files(cur_dir)) > 0){
      print(paste0("copying coops from inputs of experiment ", coop))
      setwd(cur_dir)
      file.copy(from = list.files(cur_dir) , to = dest_dir)
    }
  }
  # copy annotations
  if(length(annotations) == 1){
    cur_dir <- paste0(shared_dir, "/", "Experiment_", annotations, "/Inputs/Filtered_Annotations/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/Filtered_Annotations/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    if(length(list.files(cur_dir)) > 0){
      print(paste0("copying annotations from inputs of experiment ", annotations))
      setwd(cur_dir)
      file.copy(from = list.files(cur_dir) , to = dest_dir)
    }
  }
  # copy Enh_gene_map
  if(length(Enh_gene_map)==1){
    print(paste0("copying Enh_gene_map from inputs of experiment ", Enh_gene_map))
    cur_dir <- paste0(shared_dir, "/", "Experiment_", Enh_gene_map, "/Inputs/Enh_Gene_Map/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/Enh_Gene_Map/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    setwd(cur_dir)
    file.copy(from = list.files(cur_dir) , to = dest_dir)
  }
  # copy free_fix
  if(length(free_fix)==1){
    print(paste0("copying free_fix from inputs of experiment ", free_fix))
    cur_dir <- paste0(shared_dir, "/", "Experiment_", free_fix, "/Inputs/free_fix/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/free_fix/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    setwd(cur_dir)
    file.copy(from = list.files(cur_dir) , to = dest_dir)
  }
  # copy gene_Exp
  if(length(Gene_Exp)==1){
    print(paste0("copying gene expression from inputs of experiment ", Gene_Exp))
    cur_dir <- paste0(shared_dir, "/", "Experiment_", Gene_Exp, "/Inputs/Gene_expression/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/Gene_expression/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    setwd(cur_dir)
    file.copy(from = list.files(cur_dir) , to = dest_dir)
  }
  # copy Motifs
  if(length(Motifs)==1){
    print(paste0("copying motifs from inputs of experiment ", Motifs))
    cur_dir <- paste0(shared_dir, "/", "Experiment_", Motifs, "/Inputs/Motifs/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/Motifs/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    setwd(cur_dir)
    file.copy(from = list.files(cur_dir) , to = dest_dir)
  }
  # copy Params
  if(length(Params)==1){
    if(Params_from_output){
      print(paste0("copying params from outputs of experiment ", Params))
      cur_dir <- paste0(shared_dir, "/", "Experiment_", Params, "/Outputs/" )
      setwd(cur_dir)
      list_of_files <- list.files(cur_dir, pattern = "*.Filter")
      dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/Parameters/" )
      stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
      file.copy(from = list_of_files, to = dest_dir)
      setwd(dest_dir)
      # rename files
      aaa <- strsplit(list_of_files, split = "_")
      aaaa <- lapply(aaa, unlist)
      aaaaa <- lapply(aaaa, strsplit, split="\\.")
      aaaaaa <- lapply(aaaaa, unlist)
      aaaaaaa <- unlist(lapply(aaaaaa, "[[", 3))
      file.rename(list_of_files, paste0("start_", aaaaaaa, ".par"))
    }else{
      print(paste0("copying params from inputs of experiment ", Params))
      cur_dir <- paste0(shared_dir, "/", "Experiment_", Params, "/Inputs/Parameters/" )
      dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/Parameters/" )
      stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
      setwd(cur_dir)
      file.copy(from = list.files(cur_dir) , to = dest_dir)
    }
  }
  # copy seq
  if(length(seq)==1){
    print(paste0("copying Sequence from inputs of experiment ", seq))
    cur_dir <- paste0(shared_dir, "/", "Experiment_", seq, "/Inputs/Sequence/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/Sequence/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    setwd(cur_dir)
    file.copy(from = list.files(cur_dir) , to = dest_dir)
  }
  
  # copy TFexp
  if(length(TFexp)==1){
    print(paste0("copying TF_expression from inputs of experiment ", TFexp))
    cur_dir <- paste0(shared_dir, "/", "Experiment_", TFexp, "/Inputs/TF_Expression/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/TF_Expression/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    setwd(cur_dir)
    file.copy(from = list.files(cur_dir) , to = dest_dir)
  }
  
  # copy TFinfo
  if(length(TF_info)==1){
    print(paste0("copying TF_info from inputs of experiment ", TF_info))
    cur_dir <- paste0(shared_dir, "/", "Experiment_", TF_info, "/Inputs/TF_info/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/TF_info/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    setwd(cur_dir)
    file.copy(from = list.files(cur_dir) , to = dest_dir)
  }
  
  # copy Treat_cont_map
  if(length(Treat_cont_map)==1){
    print(paste0("copying Treat_cont_map from inputs of experiment ", Treat_cont_map))
    cur_dir <- paste0(shared_dir, "/", "Experiment_", Treat_cont_map, "/Inputs/Treat_Control_Map/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/Treat_Control_Map/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    setwd(cur_dir)
    file.copy(from = list.files(cur_dir) , to = dest_dir)
  }
  
  # copy weights
  if(length(wights)==1){
    print(paste0("copying weigths from inputs of experiment ", wights))
    cur_dir <- paste0(shared_dir, "/", "Experiment_", wights, "/Inputs/Weights/" )
    dest_dir <- paste0(shared_dir, "/", "Experiment_", cur_exp, "/Inputs/Weights/" )
    stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
    setwd(cur_dir)
    file.copy(from = list.files(cur_dir) , to = dest_dir)
  }
  # copy jobfile
  if(length(jobfile)==1){
    print(paste0("copying jobfile from inputs of experiment ", jobfile))
    job_copier(prev_exp = jobfile,
               .cur_exp=cur_exp,
               optimize=.optim,
               prev_na=.prev_na,
               .shared_dir=shared_dir)
    
  }
} 
########################################################################################################################
########################################################################################################################
# example
aa <- GEMSTAT_init_BOlinear(.exp.nu = 230,
                            .dir=aa_gd,
                            .GEMSTAT_call = "/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr",
                            TF_names = names(TF.motifs.Shrinked.count)[c(1, 2, 3, 4, 4, 5)],
                            model_evaluations = Sim_Ann_weighted_148_restart_all_models_eval[1:2],
                            model_gene_ind = c(1,2),
                            model_TF_ind = c(1, 2, 3, 4, 4, 5),
                            #model_parameters = cbind(Sim_Ann_weighted_148_restart_parameters[1:2,1:3], Sim_Ann_weighted_148_restart_parameters[1:2,20]),
                            model_parameters = Sim_Ann_weighted_148_restart_parameters[1:2,],
                            nu_enh_per_gene = 1,
                            TF_KD_evaluation = Sim_Ann_weighted_148_restart_TF_KD,
                            quantile_nu = 4,
                            enhancer_GRang = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_GRanges[1:2],
                            enhancer_Seq = ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq[1:2],
                            real_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[1:2,],
                            TF_expression_mat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01[c(1, 2, 3, 4, 4, 5), ],
                            motif_list = TF.motifs.Shrinked.count[1:3],
                            annotation_thresh=numeric(0),
                            annotation_range=numeric(0),
                            initial_bind_w=numeric(0),
                            bind_w_range=numeric(0),
                            initial_alpha=numeric(0),
                            alpha_range=numeric(0),
                            coop_tf_mat=numeric(0),
                            initial_coop_weight=numeric(0),
                            coop_weight_range=numeric(0),
                            coop_type=character(0),
                            coop_dist=integer(0),
                            coop_orientation=integer(0),
                            .one_qbtm_per_enh=T,
                            initial_qBTM=numeric(0),
                            qBTMrange=numeric(0),
                            .one_beta_per_enh=T,
                            initial_pi_beta=numeric(0), 
                            pi_beta_upper=numeric(0),
                            pi_beta_lower=numeric(0),
                            annotation_thresh_ff=integer(0),
                            initial_bind_w_ff=integer(0),
                            initial_alpha_ff=integer(0),
                            initial_coop_weight_ff=integer(0),
                            initial_qBTM_ff=1,
                            initial_pi_beta_ff=integer(0),
                            nu_enhacners = length(unlist(ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene$Chopped_Seq[1:3])),
                            .ensemble_mode = F,
                            nu_samples = 0,
                            na = 2,
                            .job_file_name = "newjob",
                            .compare_to_expectation = T,
                            TF_role_dec = T
                            , create_ff = T, create_bounds = F, create_coop = F, create_expmat = F, create_motifs = F, create_params = F, create_seq = F, create_TFexp = F, create_TFinfo = F, create_enh_gene_map = F, create_treatContMap = F, create_weights = F
                            #,annotation_thresh_ff_ens = annotation_thresh_ff,
                            #initial_coop_weight_ff_ens = initial_coop_weight_ff,
                            # initial_qBTM_ff_ens = initial_qBTM_ff
)

GEMSTAT_input_copier_multienh(lower_bounds=37,
                              upper_bounds=38,
                              coop=28, 
                              Enh_gene_map=33,
                              free_fix=integer(0),
                              Gene_Exp=33, 
                              Motifs=33,
                              Params=33, 
                              Params_from_output=T,
                              seq=33,
                              TFexp=33,
                              TF_info=33,
                              Treat_cont_map=33, 
                              wights=33, 
                              cur_exp=230, 
                              shared_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens")
########################################################################################################################
########################################################################################################################
# job and submit job copier
job_copier <- function(prev_exp, .cur_exp, optimize=T, prev_na=10, .shared_dir){
  # prev_exp is an integer indicating the index of the experiment which we want to copy from
  # .cur_exp is the index of the experiment that we want to create the job for
  # optimize is a boolean, if True is sets na to 5 if false it sets na to 0
  # .shared_dir : is the shared directory between all experiments
  prev_wd <- getwd()
  on.exit(setwd(prev_wd))
  
  print("copying job file ")
  cur_dir <- paste0(.shared_dir, "/", "Experiment_", prev_exp)
  dest_dir <- paste0(.shared_dir, "/", "Experiment_", .cur_exp)
  stopifnot(dir.exists(cur_dir), dir.exists(dest_dir))
  setwd(cur_dir)
  file.copy(from = list.files(cur_dir, pattern = "*.job") , to = dest_dir)
  setwd(dest_dir)
  
  #aa <- readLines("based_on_linear_1enh_perGene.job")
  replace_command <- paste0("sed -i ",'"" ' ,"'s/", "Experiment_", prev_exp, "/", "Experiment_", .cur_exp, "/g' ", "based_on_linear_1enh_perGene.job")
  system(replace_command)
  if(optimize){
    replace_command <- paste0("sed -i ",'"" ' ,"'s/", "-na ", prev_na, "/", "-na ", 5, "/g' ", "based_on_linear_1enh_perGene.job")
  }else{
    replace_command <- paste0("sed -i ",'"" ' ,"'s/", "-na ", prev_na, "/", "-na ", 0, "/g' ", "based_on_linear_1enh_perGene.job")
  }
  system(replace_command)
}
########################################################################################################################
########################################################################################################################
# expample
job_copier(prev_exp = 282, .cur_exp=294, optimize=F, prev_na=5, .shared_dir="/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens")
  
########################################################################################################################
########################################################################################################################
round_custume_thresh <- function(input_mat, thresh=0.5, max_neg=-1, max_pos=1){
  # input_mat: matrix to be rounded
  ouput_mat <- input_mat
  ouput_mat[abs(ouput_mat) %% 1 > thresh & ouput_mat > 0 & !is.na(ouput_mat)] <- ceiling(ouput_mat[ouput_mat %% 1 > thresh & ouput_mat > 0 & !is.na(ouput_mat)])
  ouput_mat[abs(ouput_mat) %% 1 <= thresh & ouput_mat > 0 & !is.na(ouput_mat)] <- floor(ouput_mat[ouput_mat %% 1 <= thresh & ouput_mat > 0 & !is.na(ouput_mat)])
  
  ouput_mat[abs(ouput_mat) %% 1 > thresh & ouput_mat < 0 & !is.na(ouput_mat)] <-  floor(ouput_mat[abs(ouput_mat) %% 1 > thresh & ouput_mat < 0 & !is.na(ouput_mat)])
  ouput_mat[abs(ouput_mat) %% 1 <= thresh & ouput_mat < 0 & !is.na(ouput_mat)] <- ceiling(ouput_mat[abs(ouput_mat) %% 1 <= thresh & ouput_mat < 0 & !is.na(ouput_mat)])
  ouput_mat[ouput_mat < max_neg] <- max_neg
  ouput_mat[ouput_mat > max_pos] <- max_pos
  return(ouput_mat)
}
########################################################################################################################
########################################################################################################################
#example
aa <- round_custume_thresh(cbind(c(0.1, 0.4, -0.6, 0,  -0.3, 1),
                                 c(0.2, 8.2, -3.4, 9.5, 4.7, 11.9)),
                           thresh = 0.3)
########################################################################################################################
########################################################################################################################
GEMSTAT_output_Reader_multiEnh <- function(expression_file,
                                           fold_change=T,
                                           .fixed_thresh,
                                           ignore_lab=integer(0)){
  library(sigmoid)
  # to read expression Output of GEMSTAT
  # expression_file is the .out output of GEMSTAT, where each row is gene and each column is a condition
  # ignore_lab : the labels that you want to be ignored (will be set to NA in both ground truth and prediction)
  stopifnot(file.exists(expression_file), 
            is.logical(fold_change),
            is.numeric(.fixed_thresh), 
            is.integer(ignore_lab))
  read_exp_file <- read.table(file = expression_file,
                              header = T,
                              sep = "\t", 
                              stringsAsFactors = F)
  exp_mat_GT_Pr <- as.matrix(read_exp_file[, 2:ncol(read_exp_file)])
  rownames(exp_mat_GT_Pr) <- read_exp_file[, 1]
  
  # separate ground truth from predicted matrix
  rowname_disect <- strsplit(rownames(exp_mat_GT_Pr), split="_")
  rowname_len <- unlist(lapply(rowname_disect, length))
  stopifnot(all(rowname_len %in% c(2, 3))) # make sure all the names are either two or three parts
  GT_ind <- which(rowname_len %in% 3)
  PR_ind <- which(rowname_len %in% 2)
  exp_mat_GT <- exp_mat_GT_Pr[GT_ind, ]
  exp_mat_PR <- exp_mat_GT_Pr[PR_ind, ]
  my_col_num <- ncol(exp_mat_PR)
  # compute the loss for each enhancer to choose the best performing one
  if(fold_change){
    my_col_num <- ncol(exp_mat_PR)/2
    GT_FC_mat <- matrix(nrow = nrow(exp_mat_GT), ncol = (ncol(exp_mat_GT)/2))
    PR_FC_mat <- matrix(nrow = nrow(exp_mat_PR), ncol = (ncol(exp_mat_PR)/2))
    RMSE_PerEnh <- numeric(nrow(exp_mat_PR))
    nu_nonNA_perEnh <- numeric(nrow(exp_mat_PR))
    for(cur_gene in 1 : nrow(GT_FC_mat)){# loop over enhancers
      for(cur_cond in 1:ncol(GT_FC_mat)){# loop over conditions
        if(exp_mat_GT[cur_gene, 2* cur_cond] > 1){ # find the NAs
          if(exp_mat_GT[cur_gene, 2* cur_cond] < 10){
            print("something is wrong. The NA values in the expression matrix are not 617 (are actually less than 10)")
          }
          GT_FC_mat[cur_gene, cur_cond] <- NA
          PR_FC_mat[cur_gene, cur_cond] <- NA
          next()
        }
        GT_FC_mat[cur_gene, cur_cond] <- log2(exp_mat_GT[cur_gene, 2*cur_cond]/exp_mat_GT[cur_gene, (2*cur_cond - 1)])
        PR_FC_mat[cur_gene, cur_cond] <-  2 * sigmoid(log2(exp_mat_PR[cur_gene, 2*cur_cond]/exp_mat_PR[cur_gene, (2*cur_cond - 1)])) - 1

        
      }# end of loop over conditions
      if(length(ignore_lab) > 0){
        for(clab in ignore_lab){
          cl_ind <- which(GT_FC_mat[cur_gene,] %in% clab)
          GT_FC_mat[cur_gene, cl_ind] <- NA
          PR_FC_mat[cur_gene, cl_ind] <- NA
        }
      }
      #compute root mean squared error per enhancer
      nu_nonNA_perEnh[cur_gene] <- sum(!is.na(GT_FC_mat[cur_gene, ]))
      RMSE_PerEnh[cur_gene] <- sqrt(sum((GT_FC_mat[cur_gene, ] - PR_FC_mat[cur_gene, ])^2,
                                        na.rm=T)/nu_nonNA_perEnh[cur_gene])
    }# end of  loop over enhancers
    GT_Exp_mat <- GT_FC_mat
    PR_Exp_mat <- PR_FC_mat
  }else{ # for non-fold change setting
    print("computing losses per enhancer is not implemented yet for non fold_change setting")
    GT_Exp_mat <- exp_mat_GT
    PR_Exp_mat <- exp_mat_PR
    if(length(ignore_lab) > 0){
      for(clab in ignore_lab){
        GT_Exp_mat[GT_Exp_mat %in% clab] <- NA
        PR_Exp_mat[GT_Exp_mat %in% clab] <- NA
      }
    }
    return(list(ground_truth = GT_Exp_mat,
                prediction = PR_Exp_mat))
  }
  
  # choose one enhancer per group of enhancers corresponding to a gene (choose the one with minimum error)
  PR_rowname_disect <- rowname_disect[PR_ind]
  all_gene_name <- unlist(lapply(PR_rowname_disect, "[[", 1))
  all_enh_index <- as.integer(unlist(lapply(PR_rowname_disect, "[[", 2)))
  unique_gene_name <- unique(all_gene_name)
  GT_Exp_mat_unique <- matrix(nrow = length(unique_gene_name), ncol= my_col_num)
  PR_Exp_mat_unique <- matrix(nrow = length(unique_gene_name), ncol= my_col_num)
  rownames(GT_Exp_mat_unique) <- unique_gene_name
  rownames(PR_Exp_mat_unique) <- unique_gene_name
  RMSE_perGene <- numeric(length(unique_gene_name))
  nu_nonNA_perGene <- numeric(length(unique_gene_name))
  chosen_enh_perGene<- integer(length(unique_gene_name))
  names(RMSE_perGene) <- unique_gene_name
  names(nu_nonNA_perGene) <- unique_gene_name
  names(chosen_enh_perGene) <- unique_gene_name
  for(cur_gene in 1:length(unique_gene_name)){
    cur_gene_ind <- which(all_gene_name %in% unique_gene_name[cur_gene])
    cur_MSE_min_ind <- which.min(RMSE_PerEnh[cur_gene_ind])
    chosen_ind <- cur_gene_ind[cur_MSE_min_ind]
    GT_Exp_mat_unique[cur_gene, ] <- GT_Exp_mat[chosen_ind, ]
    PR_Exp_mat_unique[cur_gene, ] <- PR_Exp_mat[chosen_ind, ]
    RMSE_perGene[cur_gene] <- RMSE_PerEnh[chosen_ind]
    nu_nonNA_perGene[cur_gene] <- nu_nonNA_perEnh[chosen_ind]
    chosen_enh_perGene[cur_gene] <- all_enh_index[chosen_ind]
  }
  discretized_exp_thresh <- prediction_discretizer(prediction = PR_Exp_mat_unique,
                                                   label = GT_Exp_mat_unique,
                                                   fixed_thresh = .fixed_thresh)
  # GT_Exp_mat_unique_round <- round_custume_thresh(GT_Exp_mat_unique, thresh = round_thresh)
  # PR_Exp_mat_unique_round <- round_custume_thresh(PR_Exp_mat_unique, thresh = round_thresh)
  # PR_Exp_mat_unique_round[PR_Exp_mat_unique_round > 1] <- 1
  # PR_Exp_mat_unique_round[PR_Exp_mat_unique_round < -1] <- -1
  PR_Exp_mat_unique_round <- discretized_exp_thresh$round_prediction
  eq_func <- function(x, y){
    all_acc <- numeric(length = nrow(x))
    for(i in 1:nrow(x)){
      all_acc[i] <- sum(x[i, ] == y[i, ], na.rm = T)/sum(!is.na(x[i, ]))
    }
    return(all_acc)
  }
  Acc_round_all <- sum(PR_Exp_mat_unique_round == GT_Exp_mat_unique, na.rm = T)/sum(!is.na(PR_Exp_mat_unique_round))
  Acc_round_Pergene <- eq_func(x = PR_Exp_mat_unique_round, y = GT_Exp_mat_unique)
  return(list(Ground_truth=GT_Exp_mat_unique,
              Prediction=PR_Exp_mat_unique,
              Ground_truth_round=GT_Exp_mat_unique,
              Prediction_round=PR_Exp_mat_unique_round,
              RootMean_sq_error_perGene = RMSE_perGene,
              Accuracy_round_perGene = Acc_round_Pergene,
              Accuracy_round_all=Acc_round_all,
              Enhancer_index =chosen_enh_perGene,
              label_thresh=discretized_exp_thresh$threshods))
}
###################################################################################################################
###################################################################################################################
#example
aa_Read_example <- GEMSTAT_output_Reader_multiEnh(expression_file = "Seeded_GEMSTAT_ens/Experiment_20/Outputs/Experiment_20_1.out",fold_change = T)
cbind(aa_Read_example$Enhancer_index, aa_Read_example$Enhancer_index)
cbind(aa_Read_example$RootMean_sq_error_perGene, aa_Read_example$RootMean_sq_error_perGene)
aa_Read_example$
###################################################################################################################
###################################################################################################################
GEMSTAT_output_Reader_multiEnh_ensemble <- function(output_dir,
                                                    .fold_change=T,
                                                    ..fixed_thresh=numeric(0),
                                                    .ignore_lab=integer(0)){
  # fixed thresh is a list of thresholds, one for each model
  # .ignore_lab : integer vector of labels that you want to be ignored
  stopifnot(dir.exists(output_dir), 
            (length(..fixed_thresh) == 0 | is.list(..fixed_thresh)))
  my_out_files <- list.files(output_dir, pattern="*.out")
  stopifnot(length(my_out_files) > 0)
  my_out_files_split_1 <- strsplit(my_out_files, split="\\.")
  my_out_files_split_2 <- lapply(my_out_files_split_1, "[[", 1)
  my_out_files_split <- strsplit(unlist(my_out_files_split_2), split="_")
  model_nu <- as.integer(unlist(lapply(my_out_files_split, "[[", 3)))
  my_out_files_sorted <- my_out_files[sort(model_nu, decreasing = F, index.return=T)$ix]
  sorted_model_nu <- sort(model_nu, decreasing = F)
  
  
  prediction_raw_list <- list()
  prediction_round_list <- list()
  groundTruth_List <- list()
  RMSE_perGene_perModel <- list()
  Accuracy_perGene_perModel <- list()
  Enhancer_index_perModel <- list()
  modeled_genes <- list()
  threshold_list <- list()
  if(!.fold_change){
    non_fold_change <- list()
  }
  accuracy_all <- numeric(length(my_out_files_sorted))
  for(i in 1:length(my_out_files_sorted)){
    print(paste("model_nu ", sorted_model_nu[i]))
    if(! .fold_change){
      cur_output_parsed <- GEMSTAT_output_Reader_multiEnh(expression_file=paste0(output_dir, "/",my_out_files_sorted[i]),
                                                          fold_change = .fold_change,
                                                          .fixed_thresh = numeric(0),
                                                          ignore_lab = .ignore_lab)
      non_fold_change[[i]] <- cur_output_parsed
    }else{ # if fold change
      if(length(..fixed_thresh) > 0){
        cur_output_parsed <- GEMSTAT_output_Reader_multiEnh(expression_file=paste0(output_dir, "/",my_out_files_sorted[i]),
                                                            fold_change = .fold_change, 
                                                            .fixed_thresh = ..fixed_thresh[[sorted_model_nu[i]]], 
                                                            ignore_lab = .ignore_lab)
      }else{
        cur_output_parsed <- GEMSTAT_output_Reader_multiEnh(expression_file=paste0(output_dir, "/",my_out_files_sorted[i]),
                                                            fold_change = .fold_change,
                                                            .fixed_thresh = numeric(0),
                                                            ignore_lab = .ignore_lab)
      }
      
      RMSE_perGene_perModel[[i]] <- cur_output_parsed$RootMean_sq_error_perGene
      Accuracy_perGene_perModel[[i]] <- cur_output_parsed$Accuracy_round_perGene
      Enhancer_index_perModel[[i]] <- cur_output_parsed$Enhancer_index
      prediction_raw_list[[i]] <- cur_output_parsed$Prediction
      prediction_round_list[[i]] <- cur_output_parsed$Prediction_round
      groundTruth_List[[i]] <- cur_output_parsed$Ground_truth
      modeled_genes[[i]] <- rownames(cur_output_parsed$Ground_truth)
      accuracy_all[i] <- cur_output_parsed$Accuracy_round_all
      threshold_list[[i]] <- cur_output_parsed$label_thresh
    }

  }
  if(! .fold_change){
    names(non_fold_change) <- sort(model_nu, decreasing = F)
    return(non_fold_change)
  }
  
  #get the union of all genes
  all_genes_names <- Reduce(union, modeled_genes)
  RMSE_perGene_perModel_mat <- matrix(nrow = length(all_genes_names), ncol = length(my_out_files_sorted))
  Accuracy_perGene_perModel_mat <- matrix(nrow = length(all_genes_names), ncol = length(my_out_files_sorted))
  Enhancer_index_perModel_mat <- matrix(nrow = length(all_genes_names), ncol = length(my_out_files_sorted))
  rownames(RMSE_perGene_perModel_mat) <- all_genes_names
  rownames(Accuracy_perGene_perModel_mat) <- all_genes_names
  rownames(Enhancer_index_perModel_mat) <- all_genes_names
  colnames(RMSE_perGene_perModel_mat) <- sort(model_nu, decreasing = F)
  colnames(Accuracy_perGene_perModel_mat) <- sort(model_nu, decreasing = F)
  colnames(Enhancer_index_perModel_mat) <- sort(model_nu, decreasing = F)
  
  # Populate the error, acc, and enh matrices using the names
  for(i in 1:length(my_out_files_sorted)){
    cur_ind <- match(names(RMSE_perGene_perModel[[i]]), rownames(RMSE_perGene_perModel_mat))
    RMSE_perGene_perModel_mat[cur_ind, i] <- RMSE_perGene_perModel[[i]]
    Accuracy_perGene_perModel_mat[cur_ind, i] <- Accuracy_perGene_perModel[[i]]
    Enhancer_index_perModel_mat[cur_ind, i] <- Enhancer_index_perModel[[i]]
  }

  names(prediction_raw_list) <- sort(model_nu, decreasing = F)
  names(prediction_round_list) <- sort(model_nu, decreasing = F)
  names(groundTruth_List) <- sort(model_nu, decreasing = F)
  names(threshold_list) <-  sort(model_nu, decreasing = F)
  # create a random prediction per model, to see the shuffled performance in cases that the gene number is different
  Random_prediction_list <- list()
  for(cur_mod in 1:length(groundTruth_List)){
    if(length(table(groundTruth_List[[cur_mod]])) < 2){
      next()
    }
    Random_prediction_list[[cur_mod]] <- CreateRandomPrediction(real_exp_mat = groundTruth_List[[cur_mod]], num = 10)
  }
  return(list(Prediction_raw_list= prediction_raw_list,
              Prediction_round_list=prediction_round_list,
              GroundTruth_List = groundTruth_List,
              RMSE=RMSE_perGene_perModel_mat,
              Accuracy=Accuracy_perGene_perModel_mat,
              Enhancer_index=Enhancer_index_perModel_mat,
              Accuracy_All=accuracy_all, 
              Random_Prediction_list=Random_prediction_list,
              threshold=threshold_list))
}
########################################################################################################################
########################################################################################################################
# example
GEMSTAT_based_on_linear_exp1 <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = "Seeded_GEMSTAT_ens/Experiment_1/Outputs", .fold_change = T)
##### ##### ##### ##### ##### ##### ##### ##### ####
########################################################################################################################
########################################################################################################################
GEMSTAT_log_Reader_multiEnh <- function(log_files_dir,
                                        role_dec = F,
                                        TFinfo_dir = character(0)){
  library(rjson)
  # to read log Output of GEMSTAT
  # log_files_dir : the directory where (processed) log files live (.Filter files)
  # role_dec : if roles where decided, meaning a TFinfo file was given as input to GEMSTAT run
  # TFinfo_dir : the directory containing the TFinfo file corresponding to the lof file
  
  stopifnot(((role_dec & length(TFinfo_dir) > 0)|(!role_dec & length(TFinfo_dir) == 0)))
  
  my_log_files <- list.files(log_files_dir, pattern = "*.Filter")
  
  my_log_files_split_1 <- strsplit(my_log_files, split="\\.")
  my_log_files_split_2 <- lapply(my_log_files_split_1, "[[", 1)
  my_log_files_split <- strsplit(unlist(my_log_files_split_2), split="_")
  model_nu <- as.integer(unlist(lapply(my_log_files_split, "[[", 3)))
  my_log_files_sorted <- my_log_files[sort(model_nu, decreasing = F, index.return=T)$ix]
  nu_models <- length(my_log_files_sorted)
  annotation_pars <- list()
  binding_pars <- list()
  alpha_activ_pars <- list()
  alpha_repre_pars <- list()
  coop_pars <- list()
  qBTM_pars <- list()
  beta_pars <- list()
  
  for(cur_model in 1:nu_models){
    cur_jason <- fromJSON(file=paste0(log_files_dir, "/", my_log_files_sorted[cur_model]))
    cur_tf_pars <- cur_jason$tfs
    cur_coop_pars <-  cur_jason$inter
    cur_qBTM_pars <- cur_jason$qbtm
    cur_beta_pars <- cur_jason$enh
    annotation_pars[[cur_model]] <- unlist(lapply(cur_tf_pars, "[[", 1) )
    binding_pars[[cur_model]] <- unlist(lapply(cur_tf_pars, "[[", 2) )
    alpha_activ_pars[[cur_model]] <- unlist(lapply(cur_tf_pars, "[[", 3) )
    alpha_repre_pars[[cur_model]] <- unlist(lapply(cur_tf_pars, "[[", 4) )
    coop_pars[[cur_model]] <- unlist(cur_coop_pars)
    qBTM_pars[[cur_model]] <- cur_qBTM_pars
    beta_pars[[cur_model]] <- unlist(lapply(cur_beta_pars, "[[", 2) )
  }
  
  annotation_par_mat <- do.call(rbind, annotation_pars)
  rownames(annotation_par_mat) <- c(1:nrow(annotation_par_mat))
  colnames(annotation_par_mat) <- paste("annotation", names(annotation_pars[[1]]), sep = "_")
  
  binding_pars_mat <- do.call(rbind, binding_pars)
  rownames(binding_pars_mat) <- c(1:nrow(binding_pars_mat))
  colnames(binding_pars_mat) <- paste("Binding", names(binding_pars[[1]]), sep = "_")
  
  alpha_activ_pars_mat <- do.call(rbind, alpha_activ_pars)
  rownames(alpha_activ_pars_mat) <- c(1:nrow(binding_pars_mat))
  colnames(alpha_activ_pars_mat) <- paste("activation", names(binding_pars[[1]]), sep = "_")
  
  alpha_repre_pars_mat <- do.call(rbind, alpha_repre_pars)
  rownames(alpha_repre_pars_mat) <- c(1:nrow(binding_pars_mat))
  colnames(alpha_repre_pars_mat) <- paste("repression", names(binding_pars[[1]]), sep = "_")
  
  
  # parsing the alpha parameters into one per TF
  my_alpha_pars <- matrix(nrow = nrow(alpha_activ_pars_mat), ncol = ncol(alpha_activ_pars_mat))
  rownames(my_alpha_pars) <- c(1:nrow(my_alpha_pars))
  colnames(my_alpha_pars) <- paste("alpha", names(binding_pars[[1]]), sep = "_")
  if(role_dec){
    for(cur_model in 1:nrow(my_alpha_pars)){
      cur_info <- read.table(paste0(TFinfo_dir, "/", "factor_info_", cur_model) 
                             , header = F, sep = "\t", stringsAsFactors = F)
      for(cur_tf in 1:ncol(my_alpha_pars)){
        stopifnot(cur_info[cur_tf, 1] == names(binding_pars[[1]])[cur_tf])
        cur_role <- ifelse(cur_info[cur_tf, 2] == 1, yes = 1, no = 0)
        if(cur_role == 1){
          my_alpha_pars[cur_model, cur_tf] <- alpha_activ_pars_mat[cur_model, cur_tf]
        }else if(cur_role == 0){
          my_alpha_pars[cur_model, cur_tf] <- alpha_repre_pars_mat[cur_model, cur_tf]
        }
      }
    }
  }else{
    my_alpha_pars[] <- alpha_activ_pars_mat[]
  }
  
  
  return(list(annotation=annotation_par_mat,
              binding=binding_pars_mat,
              activation=alpha_activ_pars_mat,
              repression=alpha_repre_pars_mat,
              alpha_effective = my_alpha_pars,
              coop = coop_pars,
              qBTM=qBTM_pars,
              Beta= beta_pars))
}
########################################################################################################################
########################################################################################################################
#example
aa <- GEMSTAT_log_Reader_multiEnh(log_files_dir = "Seeded_GEMSTAT_ens/Experiment_1/Outputs")
aa2 <- GEMSTAT_log_Reader_multiEnh(log_files_dir = "Seeded_GEMSTAT_ens/Experiment_5/Outputs",
                                   role_dec = T,
                                   TFinfo_dir= "Seeded_GEMSTAT_ens/Experiment_5/Inputs/TF_info")

aa <- read.table("Seeded_GEMSTAT_ens/Experiment_8/Inputs/TF_info/factor_info_2", header = F, sep = "\t", stringsAsFactors = F)
########################################################################################################################
########################################################################################################################
mean_RMSE_calculator <- function(GEMSTAT_output_Reader_multiEnh_ensemble_output){
  # to calculate the mean RMSE per experiment
  aa2 <- lapply(GEMSTAT_output_Reader_multiEnh_ensemble_output, "[[", 4)
  aa3 <- lapply(aa2, colSums, na.rm=T)
  aa_isna <- function(x){
    xx <- !is.na(x)
    return(colSums(xx))
  }
  aa4 <- lapply(aa2, aa_isna)
  aa_div <- function(a, b){
    return(a/b)
  }
  aa5 <- mapply(FUN = aa_div, aa3, aa4)
  boxplot.matrix(aa5, las=2, ylab="mean RMSE")
  abline(h=seq(floor(range(aa5)[1]) , ceiling(range(aa5)[2]), 0.05), col =2, lty = 4, lwd=0.5)
  return(aa5)
}
########################################################################################################################
########################################################################################################################
# example
aa <- mean_RMSE_calculator(GEMSTAT_based_on_linear_exp_results)
########################################################################################################################
########################################################################################################################

boxplot_grouped_by_label <- function(GEMSTAT_output_Reader_multiEnh_ensemble_output, 
                                     export_plot=T, 
                                     filename="boxplot_by_label.png",
                                     return_mat=F,
                                     .outline=T,
                                     model_index=c(1:length(GEMSTAT_output_Reader_multiEnh_ensemble_output$GroundTruth_List))){
  # This function gets the output of "GEMSTAT_output_Reader_multiEnh_ensemble" function and produces a boxplot where
  #  each bar represents the values associated with one true label (so there are number_of_labels * number_of_models bars in total)
  # return_mat : if True it returns the matrix created for the boxplot
  # .outline : if True it shows the outliers
  # model_index : an integer vector containing the index of the models to be plotted, the returned 
  #  matrix will have all the models not just the oned specified here
  
  my_Output <- GEMSTAT_output_Reader_multiEnh_ensemble_output
  my_table_list <- aa_table_list <- lapply(my_Output$GroundTruth_List, table)
  aa_max_lab_nu <- max(unlist(my_table_list))
  my_label_count <- length(my_table_list[[1]])
  my_model_cnt <- length(my_Output$GroundTruth_List)
  my_output_bygroupMat <- matrix(nrow = aa_max_lab_nu,
                                 ncol= (my_label_count * my_model_cnt))
  my_colnames <- character(0)
  for(model_nu in 1:my_model_cnt){
    my_colnames <- c(my_colnames, paste(rep(paste("model", model_nu, sep="_"), my_label_count),
                                        names(my_table_list[[model_nu]]), sep = "|"))
    for(lab_nu in 1:my_label_count){ # loop over label numbers to populate the my_output_bygroupMat matrix
      cur_pred <- my_Output$Prediction_raw_list[[model_nu]][my_Output$GroundTruth_List[[model_nu]]
                                                            == as.integer(names(my_table_list[[model_nu]])[lab_nu])]
      cur_pred <- cur_pred[!is.na(cur_pred)]
      my_output_bygroupMat[1:my_table_list[[model_nu]][lab_nu],
                           my_label_count*(model_nu-1) + lab_nu] <-  cur_pred
    }# end of loop over labels
    
  }# end of loop over models
  colnames(my_output_bygroupMat) <- my_colnames
  if(export_plot){
    png(filename,    # create PNG for the heat map        
        width = (length(model_index)/3) *300,        # 5 x 300 pixels
        height = 8*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
  }
  prev_par <- par()
  on.exit(par(mfrow = prev_par$mfrow, mar = prev_par$mar))
  
  par(mfrow = c(1, 1), mar=c(6,4,2,2))
  print(length(model_index))
  plot_index = integer(0)
  for(i in 1:length(model_index)){
    plot_index <- c(plot_index, 
                    c((3*(model_index[i] - 1) + 1):(3*model_index[i])))
  }
  boxplot.matrix(my_output_bygroupMat[, plot_index],
                 col=rep(c(2:(my_label_count + 1)),
                         my_model_cnt),
                 outline=.outline,
                 las=2)
  abline(v= seq((my_label_count+0.5),
                (my_label_count*my_model_cnt + 0.5),
                my_label_count),
         col = (my_label_count + 3),
         lty = 1, lwd=1)
  abline(h = seq(-1, 1, 0.05),
         lty=4, lwd = 0.2, 
         col = (my_label_count + 3))
  if(export_plot){
    dev.off()
  }
  if(return_mat){
    return(my_output_bygroupMat)
  }
}
####################################################################################################
####################################################################################################
# example
boxplot_grouped_by_label(GEMSTAT_output_Reader_multiEnh_ensemble_output = GEMSTAT_based_on_linear_exp_results_sigmoid[[31]], 
                         export_plot=T, 
                         filename="boxplot_by_label.png",
                         return_mat=F)
####################################################################################################
####################################################################################################
prediction_discretizer <- function(prediction,
                                   label,
                                   fixed_thresh=numeric(0),
                                   ignore_label = integer(0),
                                   weighted = F){
  # prediction is a numeric vector containing the predicted expression for each sample
  # label is the actual label of that sample (-1, 0 or 1)
  # the function tries to find thresholds to descritize the predictions such that it achieves
  #  the highest training accuracy
  # fixed_thresh : numeric vector of length (length(unique_labels) - 1), has the thresholds to be used for evaluation
  # ignore_label : the labels that you want to ignore
  # weighted : is a logical, if True, the discretizer finds the thresholds such that the weighted sum of errors is minimal, if
  #  False it maximizes the accuracy regardless of the number of each label present in the ground truth dataset
  # it outputs the thresholds as well as the descrized output
  
  stopifnot(length(prediction) == length(label))
  if(length(ignore_label) > 0){
    for(i in 1:length(ignore_label)){
      label[label %in% ignore_label[i]] <- NA
      prediction[label %in% ignore_label[i]] <- NA
    }
  }
  
  my_table <- table(label)
  my_labels <- sort(as.integer(names(my_table)),
                    decreasing = F)
  my_weight_mat <- label
  if(weighted){
    for(i in 1:length(my_labels)){
      my_weight_mat[label %in% my_labels[i]] <- sum(my_table)/my_table[as.character(my_labels[i])]
    }
  }else{
    my_weight_mat[!is.na(label)] <- 1
  }

  
  if(length(my_labels) == 1){ # if there is only one label (this is just an easy fix, not the best solution)
    print("for just one label, simple using 1e-3 and -1e-3 as thresholds. this is temporary and only works for -1, 0, 1 system")
    return(list(threshods=c(-0.001, 0.001),
                round_prediction=round_custume_thresh(input_mat = prediction,
                                                      thresh = 0.001,
                                                      max_neg = -1,
                                                      max_pos = 1)))
  }
  predicted_lab <- integer(length = length(prediction))
  if(length(fixed_thresh) > 0){
    stopifnot(length(fixed_thresh) == (length(my_labels) - 1))
    best_thresh <- fixed_thresh
  }else{ # if a fix thresh is not provided: find the best threshold
    label_average <- numeric(0)
    for (cur_lab in 1:length(my_labels)){
      cur_Avg <- mean(prediction[label == my_labels[cur_lab]], na.rm=T)
      label_average <- c(label_average, cur_Avg)
    }
    names(label_average) <- my_labels
    # print(label_average)
    
    # checking if the average of labels are in order
    for(i in 1:(length(my_labels)-1)){
      if(label_average[i+1] <= label_average[i]){
        print(paste0("BAD predictions. average of predictions with label ",
                     my_labels[i+1],
                     " is less than the average of predictions with label  ",
                     my_labels[i], ". Setting threshold to 0.1."))
        return(list(threshods=c(-0.1, 0.1),
                    round_prediction=round_custume_thresh(input_mat = prediction,
                                                          thresh = 0.1,
                                                          max_neg = -1,
                                                          max_pos = 1)))
      }
    }
    
    potential_thresh <- list()
    best_thresh <- numeric(length = (length(my_labels)-1))
    aa_names <- character(0)
    for(i in 1:(length(my_labels) - 1)){
      aa_names <- c(aa_names, paste(my_labels[i], my_labels[i+1], sep = "_"))
    }
    names(best_thresh) <- aa_names
    for(cur_thresh in 1:(length(my_labels)-1)){
      # print("cur_thresh")
      # print(cur_thresh)
      my_cur_pot <- as.numeric(prediction[prediction > label_average[cur_thresh] &
                                            prediction < label_average[cur_thresh + 1]])
      potential_thresh[[cur_thresh]] <- my_cur_pot[!is.na(my_cur_pot)]
      # print(paste0("there are ", 
      #              length(potential_thresh[[cur_thresh]]),
      #              " candidates for the threshold between ",
      #              my_labels[cur_thresh],
      #              " and ",
      #              my_labels[cur_thresh + 1]))
      if(length(potential_thresh[[cur_thresh]]) == 0){
        print(paste0("BAD predictions. there are no candidates for threshold between ",
                     as.character(my_labels[cur_thresh+1]),
                     " and ",
                     as.character(my_labels[cur_thresh]),
                     ". Setting threshold to 0.1."))
        return(list(threshods=c(-0.1, 0.1),
                    round_prediction=round_custume_thresh(input_mat = prediction,
                                                          thresh = 0.1,
                                                          max_neg = -1,
                                                          max_pos = 1)))
      }
      
      # try candidate thresholds one by one and keep the ones with highest accuracy
      prev_loss <- length(label)
      for(cur_cand in 1:length(potential_thresh[[cur_thresh]])){
        trying_pred <- prediction
        trying_pred[trying_pred <= potential_thresh[[cur_thresh]][cur_cand]] <- my_labels[cur_thresh]
        trying_pred[trying_pred >  potential_thresh[[cur_thresh]][cur_cand]] <- my_labels[cur_thresh + 1]
        cur_pred <- trying_pred[label == my_labels[cur_thresh] | label == my_labels[cur_thresh + 1]]
        cur_label <- label[label == my_labels[cur_thresh] | label == my_labels[cur_thresh + 1]]
        cur_weight <- my_weight_mat[label == my_labels[cur_thresh] | label == my_labels[cur_thresh + 1]]
        
        # check if sizes and the NAs of label and prediction match
        stopifnot(length(cur_label) == length(cur_pred)
                  , sum(is.na(cur_label) == is.na(cur_pred)) == length(cur_label))
        cur_loss <- sum((cur_label != cur_pred) * cur_weight , na.rm = T)
        
        # /(sum(!is.na(cur_pred)))
        # set a new best threshold if the accuracy is higher
        if(cur_loss < prev_loss){
          prev_loss <- cur_loss
          best_thresh[cur_thresh] <- potential_thresh[[cur_thresh]][cur_cand]
          #print(paste0("best current acc, in distinguishing between ", my_labels[cur_thresh], " and ", my_labels[cur_thresh + 1]))
          #print(prev_loss)
        }
      } # end of loop over candidates
    }# end of loop over thresholds 
  } # end of computing best thresh if its not given
  final_prediction <- prediction
  if(length(best_thresh) == 1){
    final_prediction[prediction <= best_thresh] <- my_labels[1]
    final_prediction[prediction > best_thresh] <- my_labels[2]
  }else{
    for(i in 1:(length(best_thresh))){
      if(i ==1){
        final_prediction[prediction <= best_thresh[i]] <- my_labels[i]
        #final_prediction[prediction > best_thresh[i] & prediction <= best_thresh[i+1]] <- my_labels[i+1]
      }else{
        final_prediction[prediction <= best_thresh[i] & prediction >  best_thresh[i-1]] <- my_labels[i]
        #final_prediction[prediction >  best_thresh[i] & prediction <= best_thresh[i+1]] <- my_labels[i+1]
      }
    } # end of loop over thresholds
    final_prediction[prediction > best_thresh[length(best_thresh)]] <- my_labels[length(my_labels)]
  } # if there is more than one threshold
  
  final_Accuracy <- sum(final_prediction == label, na.rm=T)/sum(!is.na(label))
  print("final_Accuracy")
  print(final_Accuracy)
  return(list(threshods=best_thresh,
              round_prediction=final_prediction))
}

save(list = c("GEMSTAT_output_Reader_multiEnh", 
              "GEMSTAT_output_Reader_multiEnh_ensemble",
              "prediction_discretizer", 
              "read_KD_results"), 
     file = "Updated_read_expression_functions.RData")
####################################################################################################
####################################################################################################
# example
aa <- prediction_discretizer(prediction = GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[53]]$Prediction_raw_list[[25]],
                             label = GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[53]]$GroundTruth_List[[25]])
####################################################################################################
####################################################################################################
coop_parameter_creator_comb <- function(TF_names,
                                        TF_index, 
                                        excluded_pairs=numeric(0),
                                        given_coop_mat=numeric(0),
                                        new_weight_range=c(0.001, 100),
                                        new_initial_weight = 1,
                                        new_coop_type="SIMPLE",
                                        new_coop_dist=50,
                                        new_coop_orientation=c(0, 0),
                                        base_mat=matrix(nrow=0, ncol=2),
                                        base_weight_range=matrix(nrow=0, ncol=2),
                                        base_initial_coop_weight=numeric(0),
                                        base_coop_type=character(0),
                                        base_coop_dist=integer(0),
                                        base_coop_orientation=matrix(nrow=0, ncol=2),
                                        only_pair_index=F
){
  # TF_names is character vector containing the names of the TFs
  # TF_index is an integer vector index of TFs which I want to consider within the TF_names
  # excluded_pairs is a list containing tuples that I don't want to consider.
  # base_mat is a matrix to be used as the base coop mat. with two columns. (eg for the case
  #  that I want some coop to be present in all matrices)
  # given_coop_mat :is a matrix the same format as base_mat, if provided the function won't compute the combinations
  # only_pair_index : if True only computes and returns the pair index (used for updating the experiment data files)
  # a function to create all possible coop_tf_mat, my_coop_weight_range, my_initial_coop_weight, my_coop_type, my_coop_dist, my_coop_orientation
  
  stopifnot(all(TF_index <= length(TF_names)),
            (length(excluded_pairs) ==0 | is.list(excluded_pairs)))
  
  if(length(given_coop_mat) > 0){
    pair_ind <- given_coop_mat
  }else{
    pair_ind <- matrix(nrow=0, ncol = 2)
    TF_index <- sort(TF_index, decreasing = F)
    for(i in 1:(length(TF_index) - 1)){
      for(j in (i+1):length(TF_index)){
        rm_flag<-0
        if(length(excluded_pairs) > 0){ # remove the excluded pairs
          
          for(k in 1:length(excluded_pairs)){
            if(TF_index[i] == excluded_pairs[[k]][1] & TF_index[j] == excluded_pairs[[k]][2]){
              rm_flag <- 1
              break()
            }
          }
        }
        if(rm_flag == 0){
          pair_ind <- rbind(pair_ind, c(TF_index[i], TF_index[j]))
        }
      }
    }
  }
  if(only_pair_index){
    return(pair_ind)
  }
  
  coop_mat_list <- list()
  new_index <- integer(nrow(pair_ind))
  
  for(cur_row in 1:nrow(pair_ind)){
    cur_coop_mat <- rbind(base_mat, pair_ind[cur_row, ])
    cur_coop_order <- order(cur_coop_mat[,1], cur_coop_mat[,2])
    new_index[cur_row] <- which(cur_coop_order == nrow(cur_coop_mat))
    cur_coop_mat <- cur_coop_mat[cur_coop_order, ]
    coop_mat_list[[cur_row]] <- cur_coop_mat
    if(is.null(nrow(coop_mat_list[[cur_row]]))){
      coop_mat_list[[cur_row]] <- matrix(coop_mat_list[[cur_row]],
                                         byrow = T,
                                         ncol = 2)
    }
  }
  
  # updating the other coop associated things
  weight_range_list <- list()
  initial_coop_weight_list <- list()
  coop_type_list <- list()
  coop_dist_list <- list()
  coop_orientation_list <- list()
  for(i in 1:length(new_index)){
    if(new_index[i] == 1){
      weight_range_list[[i]] <- rbind(new_weight_range, 
                                      base_weight_range)
      initial_coop_weight_list[[i]] <- c(new_initial_weight,
                                         base_initial_coop_weight)
      coop_type_list[[i]] <- c(new_coop_type, base_coop_type)
      coop_dist_list[[i]] <- c(new_coop_dist, base_coop_dist)
      coop_orientation_list[[i]] <- rbind(new_coop_orientation, base_coop_orientation)
    }else if(new_index[i] == nrow(coop_mat_list[[i]])){
      weight_range_list[[i]] <- rbind(base_weight_range,
                                      new_weight_range)
      initial_coop_weight_list[[i]] <- c(base_initial_coop_weight, new_initial_weight)
      coop_type_list[[i]] <- c(base_coop_type, new_coop_type)
      coop_dist_list[[i]] <- c(base_coop_dist, new_coop_dist)
      coop_orientation_list[[i]] <- rbind(base_coop_orientation, new_coop_orientation)
    }else{
      weight_range_list[[i]] <- rbind(base_weight_range[1:(new_index[i]-1), ],
                                      new_weight_range,
                                      base_weight_range[new_index[i]:nrow(base_weight_range), ])
      initial_coop_weight_list[[i]] <- c(base_initial_coop_weight[1:(new_index[i]-1)],
                                         new_initial_weight,
                                         base_initial_coop_weight[new_index[i]:length(base_initial_coop_weight)])
      coop_type_list[[i]] <- c(base_coop_type[1:(new_index[i]-1)], 
                               new_coop_type,
                               base_coop_type[new_index[i]:length(base_coop_type)])
      coop_dist_list[[i]] <- c(base_coop_dist[1:(new_index[i]-1)],
                               new_coop_dist,
                               base_coop_dist[new_index[i]:length(base_coop_dist)])
      coop_orientation_list[[i]] <- rbind(base_coop_orientation[1:(new_index[i]-1), ],
                                          new_coop_orientation,
                                          base_coop_orientation[new_index[i]:nrow(base_coop_orientation),])
    }
  }
  return(list(Coop_mat_list = coop_mat_list,
              Weight_range_list=weight_range_list,
              Initial_weight_list=initial_coop_weight_list,
              Coop_type_list=coop_type_list,
              Coop_dist_list=coop_dist_list,
              Coop_orientation_list=coop_orientation_list))
}
########################################################################################################################
########################################################################################################################
#example
aa_cc <- coop_parameter_creator_comb(TF_names=names(TF.motifs.Shrinked.halfsites.count),
                                     TF_index=c(1:20), 
                                     excluded_pairs=list(c(13, 14)),
                                     new_weight_range=c(0.001, 100),
                                     new_initial_weight = 1,
                                     new_coop_type="SIMPLE",
                                     new_coop_dist=50,
                                     new_coop_orientation=c(0, 0),
                                     base_mat=rbind(c(1, 1), c(3, 3), c(10, 10), c(13, 14), c(15, 15), c(16, 16)),
                                     base_weight_range=cbind(c(100, 100, 100, 100, 100, 100),
                                                             c(500, 500, 500, 500, 500, 500)),
                                     base_initial_coop_weight=c(300, 300, 300, 300, 300, 300),
                                     base_coop_type=rep("DIMER", 6),
                                     base_coop_dist=c(4, 4, 4, 4, 6, 7),
                                     base_coop_orientation=rbind(c(1, -1), c(1, -1), c(1, -1), c(1, 1), c(1, 1), c(1, 1))
)
####################################################################################################
####################################################################################################
annotation_reader <- function(annot_file, TF_names){
  # annot_file is the address of the annotation file
  # TF_names names of the TF which we want to include in the output
  annot_raw <- read.delim(file =annot_file, header = F, sep = "\t", stringsAsFactors = F)
  annot_raw_na <- is.na(annot_raw)
  enh_name_index <- which(rowSums(annot_raw_na) > 0)
  enh_names <- annot_raw[enh_name_index, 1]
  enh_names <- gsub(">","", enh_names)
  nu_enh <- length(enh_names)
  annot_per_enh <- list()
  
  #separeting the annotations by enh
  # TODO: UPDATE so separetes by genes and enhancers, for now treating the enhancers as genes
  for(cur_enh in 1:(length(enh_name_index) - 1)){
    annot_per_enh[[cur_enh]] <- annot_raw[(enh_name_index[cur_enh] + 1):(enh_name_index[cur_enh+1] - 1),]
    rownames(annot_per_enh[[cur_enh]]) <- c(1:nrow(annot_per_enh[[cur_enh]]))
  }
  annot_per_enh[[length(enh_name_index)]] <- annot_raw[(enh_name_index[length(enh_name_index)] + 1):(nrow(annot_raw)), ]
  names(annot_per_enh) <- enh_names
  #separating within each gene by TFs
  annot_per_enh_perTF <- list()
  for(cur_enh in 1:nu_enh){
    annot_per_enh_perTF[[cur_enh]] <- list()
    for(cur_TF in 1:length(TF_names)){
      cur_ind <- which(annot_per_enh[[cur_enh]][, 3] == TF_names[cur_TF])
      annot_per_enh_perTF[[cur_enh]][[cur_TF]] <- annot_per_enh[[cur_enh]][cur_ind, ]
    }
    names(annot_per_enh_perTF[[cur_enh]]) <- TF_names
  }
  names(annot_per_enh_perTF) <- enh_names
  return(annot_per_enh_perTF)
}
####################################################################################################
####################################################################################################
#example

GEMSTAT_based_on_linear_exp_annotation_54[[i]] <- annotation_reader(annot_file = paste0("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_54/Annotation/Experiment_54_", i, ".annot"), 
                                                                    TF_names = names(TF.motifs.Shrinked.halfsites.count))
####################################################################################################
####################################################################################################
KD_analyzer <- function(base_experiment, KD_experiments, KD_names, top_nu=15, model_index=integer(0), plot_res=T, fixed_ylim=T){
  # base_experiment is a list containing the results of the base experiment (output of GEMSTAT_output_Reader_multiEnh_ensemble for that experiment)
  # KD_experiments is a list of lists, each contatining the output of GEMSTAT_output_Reader_multiEnh_ensemble for one KD experiment
  # KD_names : names of KD experiments
  # top_nu : the top 'top_nu'  models will be considered in terms of accuracy. if model_index is not provided. otherwise modelindex will be used
  # model_index : the index of models to be used. if not provided top_nu best performing models will be chosen
  # plot_res : if True plots the results
  # fixed_ylim : if True it fixes the ylim of plots to c(0, 1), if False they will be set automatically
  stopifnot(length(KD_experiments) == length(KD_names))
  
  if(length(model_index) == 0){
    aa_goodind <- sort(base_experiment$Accuracy_All,
                       decreasing = T, index.return=T)$ix[1:top_nu]
  }else{
    aa_goodind <- model_index
  }
  
  # defining list to store number of correctly classified and incorrectly classified datapoints for each model
  aa_miscalssified_list <- list()
  aa_miscalssified_list_sep <- list()
  aa_corcalssified_list <- list()
  aa_ratio_list <- list()
  aa_ratio_list_sep <- list()
  
  for(aa_cmodif in 1:length(KD_experiments)){
    aa_miscalssified_list[[aa_cmodif]] <- matrix(nrow = length(aa_goodind), ncol = 3)
    colnames(aa_miscalssified_list[[aa_cmodif]]) <- c(-1, 0, 1)
    aa_corcalssified_list[[aa_cmodif]] <- matrix(nrow = length(aa_goodind), ncol = 3)
    colnames(aa_corcalssified_list[[aa_cmodif]]) <- c(-1, 0, 1)
    aa_ratio_list[[aa_cmodif]]         <- matrix(nrow = length(aa_goodind), ncol = 3)
    colnames(aa_ratio_list[[aa_cmodif]]) <- c(-1, 0, 1)
    aa_miscalssified_list_sep[[aa_cmodif]] <- matrix(nrow = length(aa_goodind), ncol = 6)
    aa_ratio_list_sep[[aa_cmodif]]         <- matrix(nrow = length(aa_goodind), ncol = 6)
    colnames(aa_miscalssified_list_sep[[aa_cmodif]]) <- c("-1_0", "-1_1", "0_-1", "0_1", "1_0", "1_-1")
    colnames(aa_ratio_list_sep[[aa_cmodif]]) <- colnames(aa_miscalssified_list_sep[[aa_cmodif]])
    
    rownames(aa_miscalssified_list[[aa_cmodif]]) <- aa_goodind
    rownames(aa_miscalssified_list_sep[[aa_cmodif]]) <- aa_goodind
    rownames(aa_corcalssified_list[[aa_cmodif]]) <- aa_goodind
    rownames(aa_ratio_list[[aa_cmodif]])<- aa_goodind
    rownames(aa_ratio_list_sep[[aa_cmodif]])<- aa_goodind
    # looping over models
    for(aa_cmodel in 1:length(aa_goodind)){
      aa_gt <- base_experiment$GroundTruth_List[[aa_goodind[aa_cmodel]]]
      aa_afterKD <- KD_experiments[[aa_cmodif]]$Prediction_round_list[[aa_goodind[aa_cmodel]]]
      aa_beforeKD <- base_experiment$Prediction_round_list[[aa_goodind[aa_cmodel]]]
      aa_missed <- aa_gt[aa_beforeKD == aa_gt & aa_afterKD != aa_gt]
      aa_allcor <- aa_gt[aa_beforeKD == aa_gt]
      
      aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 1] <- sum(aa_allcor == -1, na.rm = T)
      aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 1] <- sum(aa_missed == -1, na.rm = T)
      aa_ratio_list[[aa_cmodif]][aa_cmodel, 1] <- aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 1]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 1]
      
      aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 1] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == -1 & aa_afterKD == 0, na.rm = T)
      aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 2] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == -1 & aa_afterKD == 1, na.rm = T)
      aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 1] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 1]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 1]
      aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 2] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 2]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 1]
      
      aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 2] <- sum(aa_allcor ==  0, na.rm = T)
      aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 2] <- sum(aa_missed ==  0, na.rm = T)
      aa_ratio_list[[aa_cmodif]][aa_cmodel, 2] <- aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 2]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 2]
      
      
      aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 3] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == 0 & aa_afterKD == -1, na.rm = T)
      aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 4] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == 0 & aa_afterKD == 1, na.rm = T)
      aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 3] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 3]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 2]
      aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 4] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 4]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 2]
      
      aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 3] <- sum(aa_allcor ==  1, na.rm = T)
      aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 3] <- sum(aa_missed ==  1, na.rm = T)
      aa_ratio_list[[aa_cmodif]][aa_cmodel, 3] <- aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 3]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 3]
      
      aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 5] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == 1 & aa_afterKD == 0, na.rm = T)
      aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 6] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == 1 & aa_afterKD == -1, na.rm = T)
      aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 5] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 5]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 3]
      aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 6] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 6]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 3]
    } # end of loop over models
  } # end of loop over KD experiments
  names(aa_miscalssified_list) <- KD_names
  names(aa_corcalssified_list) <- KD_names
  names(aa_ratio_list) <- KD_names
  names(aa_ratio_list_sep) <- KD_names
  names(aa_miscalssified_list_sep) <-  KD_names
  
  if(plot_res){
    prev_par <- par()
    on.exit(par(mfrow = prev_par$mfrow, mar = prev_par$mar))
    
    par(mfrow=c(4, 6), mar = c(3, 3, 4, 1))
    for(i in 1:length(KD_experiments)){
      if(!fixed_ylim){
        my_ylim = range(aa_ratio_list_sep[[i]])
      }else{
        my_ylim <- c(0, 1)
      }
      boxplot.matrix(aa_ratio_list_sep[[i]],
                     main = names(aa_ratio_list)[i],
                     ylim = my_ylim, 
                     col = c(2,2, 3,3, 4,4), 
                     las = 2)
    }
    
    for(i in 1:length(KD_experiments)){
      if(!fixed_ylim){
        my_ylim = range(aa_ratio_list[[i]])
      }else{
        my_ylim <- c(0, 1)
      }
      boxplot.matrix(aa_ratio_list[[i]],
                     main = names(aa_ratio_list)[i],
                     ylim = my_ylim, 
                     col = c(2, 3, 4), 
                     las = 2)
    }
  }
  return(list(ratio_list_separated=aa_ratio_list_sep, ratio_list_combined=aa_ratio_list))
}
####################################################################################################
####################################################################################################
# example
aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[56:77],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=15, model_index=integer(0), plot_res=T, fixed_ylim = F)


####################################################################################################
####################################################################################################
bash_batch_run_zip_transfer_process <- function(exp_index, cr_run=T, cr_zip=T, cr_transfer=T, cr_process=T, cr_Rscript=T){
  # creates bash scripts for batch running, transfering, zippting and processing of experiments
  stopifnot(length(exp_index) > 1)
  if(cr_run){
    print("run file")
    cat("#!/bin/bash\n", file = paste0("run_", exp_index[1], 
                                       "_", exp_index[length(exp_index)], ".sh")
        , append = F)
    for(i in exp_index){
      cat(c(paste0("cd Experiment_", i),
            "./hal_sub_creator_based_on_linear_1enh_perGene.job",
            "chmod +x based_on_linear_1enh_perGene.job.submit", 
            "./based_on_linear_1enh_perGene.job.submit",
            "cd .."), sep = "\n",
          file = paste0("run_", exp_index[1], "_", 
                        exp_index[length(exp_index)], ".sh"),
          append = T)
    }
  }
  
  if(cr_zip){
    print("zip file")
    cat("#!/bin/bash\n", file = paste0("zip_", exp_index[1], 
                                       "_", exp_index[length(exp_index)], ".sh"),
        append = F)
    for(i in exp_index){
      cat(c(paste0("cd Experiment_", i),
            paste0("zip exp", i, "_out.zip", " Outputs/*"),
            "cd .."), sep = "\n",
          file =  paste0("zip_", exp_index[1], "_", exp_index[length(exp_index)], ".sh"),
          append = T)
    }
  }
  
  if(cr_transfer){
    print("transfer file")
    cat("#!/bin/bash\n", file = paste0("transfer_", exp_index[1], "_", exp_index[length(exp_index)], ".sh")
        , append = F)
    for(i in  exp_index){
      cat(c(paste0("rsync Experiment_", i,
                   "/exp", i,
                   "_out.zip ~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_",
                   i,"/")), sep = "\n",
          file = paste0("transfer_", exp_index[1], "_", exp_index[length(exp_index)], ".sh"),
          append = T)
    }
  }
  
  if(cr_process){
    print("process file")
    cat("#!/bin/bash\n", file = paste0("process_", exp_index[1],
                                       "_", exp_index[length(exp_index)],
                                       ".sh"),
        append = F)
    for(i in  exp_index){
      cat(c(paste0("cd ~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_", i, "/"),
            "rm -r Outputs", 
            paste0("unzip exp",i,"_out.zip"),
            "cp ../log_filter.sh .",
            "./log_filter.sh"), sep = "\n",
          file = paste0("process_", exp_index[1],
                        "_", exp_index[length(exp_index)],
                        ".sh"),
          append = T)
    }
  }
  
  if(cr_Rscript){
    print("run Rscripts")
    cat("#!/bin/bash\n", file = paste0("Seeded_GEMSTAT_ens/Rscript_", exp_index[1],
                                       "_", exp_index[length(exp_index)],
                                       ".sh"),
        append = F)
    for(i in exp_index){
      cat(paste0("Rscript Ensemble_based_on_linear_creator_noRole_exp", i, ".R\n"),
          file = paste0("Seeded_GEMSTAT_ens/Rscript_", exp_index[1],
                        "_", exp_index[length(exp_index)],
                        ".sh"),
          append = T)
    }
  }
  
}
####################################################################################################
####################################################################################################
#example
bash_batch_run_zip_transfer_process(exp_index = c(599:620))
####################################################################################################
####################################################################################################
GenerateComb <- function(Set,TotalNumberOfTFs=12, my_alphabet=c("0", "1")){
  # This Function Generates a matrix where each row is a MotBinary for all combinations
  #  of the motifs in a given set
  # Set is the set of indecis of Motifs present in this combination eg : c(1,6,11)
  # TotalNumberOfTFs is the total number of TFs in the pwmList that is going to be used
  
  library(tcR)
  k <- length(Set)
  MotBinary <- matrix(0L,nrow = (length(my_alphabet)^k),ncol = TotalNumberOfTFs)
  AllCombs <- generate.kmers(.k = k,.alphabet = my_alphabet)
  AllCombMat <- matrix(0L,nrow = length(AllCombs),ncol = k)
  for(i in 1:length(AllCombs)){
    AllCombMat[i,] <- as.numeric(substring(AllCombs[i], seq(1,k,1), seq(1,k,1)))
  }
  # for(i in 1:nrow(AllCombMat)){
  #   MotBinary[i,c(Set * AllCombMat[i,])] <- 1
  # }
  return(AllCombMat)
}


####################################################################################################
####################################################################################################
count_site_from_annotation <- function(annotation_list, 
                                       TF_names,
                                       TF_index=c(1:length(TF_names)),
                                       homoDimer, 
                                       dimer_orientation,
                                       homoDimer_distance,
                                       heterodimer_pair,
                                       heterodimer_distance,
                                       MAXLLR,
                                       annotation_thresh){
  # This function takes in a list where each entry is an output of annotation_reader function for one
  #  model and counts the number of occurances of each motif for each enhancer, taking homodimers and heterodimers into account
  # annotation_list : is a list where each entry is an output of annotation_reader function
  # TF_names : character vector containing names of all the TFs
  # TF_index : integer vecrtor containing the index of the TFs that we want to count
  # homoDimer : logical vector with length equal to number of TFs, True if
  #  the TF is a homodimer, F if otherwise (monomer or heterodimer)
  # homoDimer_distance : integer vector with length equal to number of TFs. contains 0 for non-homodimer
  #  TFs and a number indicating the number of bases between the two monomers for homodimers
  # heterodimer_pair : is a list where each entry is a tuple contaning the index of two motifs which form
  #  a heterodimer. the first number in each tuple has to be smaller than the second item
  # heterodimer_distance <- is an integer vector with length equal to length of heterodimer_pair
  #  list. showing the number of bases between the two motifs 
  # dimer_orientation is an integer vector, each entry can take one of the four values 0, 1, 2, 3
  #  0 is used for entries which are not homodimers, 1 means both instances on the
  #  same strand ((+,+) or (-,-)), 2 means (+, -), 3 means (-, +)
  # MAXLLR is a numeric vector containing the MAXLLR for each motif. length equal to length of TF_names
  # annotation_thresh is a numeric vector containing the threshold for  (LLR_max - LLR_current)/LLR_max . same as threshold fed to GEMSTAT
  
  # checking the sizes of inputs
  stopifnot(length(annotation_list) > 0,
            length(annotation_list[[1]][[1]]) == length(TF_names), 
            all(TF_index <= length(TF_names)) ,
            length(homoDimer) == length(TF_names),
            length(homoDimer_distance) == length(TF_names),
            length(heterodimer_pair) == length(heterodimer_distance),
            length(MAXLLR) == length(TF_names), 
            length(annotation_thresh) == length(TF_names))
  
  nu_models <- length(annotation_list)
  nu_genes <- length(annotation_list[[1]])
  nu_TFs <- length(TF_names)
  homo_index <- which(homoDimer)
  hetero_second <- integer(0)
  for(i in heterodimer_pair){
    hetero_second <- c(hetero_second, i[2])
  }
  # getting the occurances of each motif above its threshold for all genes and all models
  all_occur_ab_th <- list()
  tf_gene_site_mat <- list()
  for(i in 1:nu_models){
    tf_gene_site_mat[[i]] <- matrix(0L, nrow = nu_genes, ncol = nu_TFs)
    colnames(tf_gene_site_mat[[i]]) <- TF_names
  }
  for(cur_tf in TF_index){
    #print(paste0("TF number ", cur_tf))
    tmp_occur <- list()
    all_occur_ab_th[[cur_tf]] <- list()
    for(cur_model in 1:nu_models){
      #print(paste0("model number ", cur_model))
      tmp_occur[[cur_model]] <- lapply(annotation_list[[cur_model]], "[[", cur_tf)
      all_occur_ab_th[[cur_tf]][[cur_model]] <- list()
      for(cur_gene in 1:nu_genes){
        #print(paste0("gene number ", cur_gene))
        all_occur_ab_th[[cur_tf]][[cur_model]][[cur_gene]] <- tmp_occur[[cur_model]][[cur_gene]][tmp_occur[[cur_model]][[cur_gene]]$V4 <= annotation_thresh[cur_tf]*MAXLLR[cur_tf], ]
        if (cur_tf %in% homo_index){ # processing homodimers
          print(paste0("processing homodimer", TF_names[cur_tf]))
          # creating the dimer oreientation based on input
          if(dimer_orientation[cur_tf] == 1){
            cur_orien <- list(c("+","+"),c("-","-"))
          }else if(dimer_orientation[cur_tf] == 2){
            cur_orien <- list(c("+","-"))
          }else if(dimer_orientation[cur_tf] == 3){
            cur_orien <- list(c("-","+"))
          }
          
          if(nrow(all_occur_ab_th[[cur_tf]][[cur_model]][[cur_gene]]) > 1){ # if there are any above threshold sites
            aa_pos <- do.call(rbind ,strsplit(all_occur_ab_th[[cur_tf]][[cur_model]][[cur_gene]][,1],
                                              split = "\\.."))
            for(cur_site in 1:(nrow(aa_pos) - 1)){ # loop over sites
              if (((as.integer(aa_pos[cur_site+1, 1]) - as.integer(aa_pos[cur_site, 2])) <= (homoDimer_distance[cur_tf] + 1)) &
                  ((as.integer(aa_pos[cur_site+1, 1]) - as.integer(aa_pos[cur_site, 2])) > 0) &
                  list(c(all_occur_ab_th[[cur_tf]][[cur_model]][[cur_gene]][cur_site, 2],
                         all_occur_ab_th[[cur_tf]][[cur_model]][[cur_gene]][cur_site + 1, 2])) %in%  cur_orien){ # if the sites are in the right distance and orientation
                tf_gene_site_mat[[cur_model]][cur_gene, cur_tf] <- tf_gene_site_mat[[cur_model]][cur_gene, cur_tf] + 1
              } # if the sites are in the right distance and orientation
            }# end of loop over sites
          } # if there are any above threshold sites
          # end of processing homodimers
        }else if(cur_tf %in% hetero_second){
          
          # determine the index of the pair 
          cur_tf_het_ind <- which(hetero_second == cur_tf)
          other_tf <- heterodimer_pair[[cur_tf_het_ind]][1]
          print(paste0("processing heterodimer", TF_names[cur_tf], " and ", TF_names[other_tf]))
          tf_gene_site_mat[[cur_model]][cur_gene, other_tf] <- 0
          if(nrow(all_occur_ab_th[[cur_tf  ]][[cur_model]][[cur_gene]]) > 0 &
             nrow(all_occur_ab_th[[other_tf]][[cur_model]][[cur_gene]]) > 0){ # if both motifs of heterodimer have sites
            aa_pos_1 <- do.call(rbind ,strsplit(all_occur_ab_th[[other_tf]][[cur_model]][[cur_gene]][,1],
                                                split = "\\.."))
            aa_pos_2 <- do.call(rbind ,strsplit(all_occur_ab_th[[cur_tf]][[cur_model]][[cur_gene]][,1],
                                                split = "\\.."))
            for(cur_site_1 in 1:nrow(aa_pos_1)){
              for(cur_site_2 in 1:nrow(aa_pos_2)){
                if( (as.integer(aa_pos_2[cur_site_2, 1]) - as.integer(aa_pos_1[cur_site_1, 2])) <= (heterodimer_distance[cur_tf_het_ind] + 1) &
                    (as.integer(aa_pos_2[cur_site_2, 1]) - as.integer(aa_pos_1[cur_site_1, 2])) > 0 & 
                    (all_occur_ab_th[[cur_tf]][[cur_model]][[cur_gene]][cur_site_2, 2] == 
                     all_occur_ab_th[[other_tf]][[cur_model]][[cur_gene]][cur_site_1, 2]) ){
                  tf_gene_site_mat[[cur_model]][cur_gene, cur_tf] <- tf_gene_site_mat[[cur_model]][cur_gene, cur_tf] + 1
                  tf_gene_site_mat[[cur_model]][cur_gene, other_tf] <- tf_gene_site_mat[[cur_model]][cur_gene, other_tf] + 1
                } # if they are in the right distance and orientation
              } # end  of loop over sites of second half
            } # end  of loop over sites of first half
          } # if both motifs of heterodimer have sites
          # end of processing heterodimers
        }else{
          print("nrow(all_occur_ab_th[[cur_tf]][[cur_model]][[cur_gene]])")
          print(nrow(all_occur_ab_th[[cur_tf]][[cur_model]][[cur_gene]]))
          tf_gene_site_mat[[cur_model]][cur_gene, cur_tf] <- nrow(all_occur_ab_th[[cur_tf]][[cur_model]][[cur_gene]])
        }
        
      }#end of loop over genes
      
    } # end of loop over models
  } # end of loop over TFs
  
  
  return(list(Sites_Per_TF = all_occur_ab_th,
              Site_count_mat = tf_gene_site_mat))
}
#########################################################################################################################
#########################################################################################################################
# example
aa_homodimer <- c(T, F, T, F, F, F, F, F, F, T, F, F, F, F, T, T, F, F, F, F)
names(aa_homodimer) <- names(TF.motifs.Shrinked.halfsites.count_2)
aa_dimer_orien <- c(2, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0)
names(aa_dimer_orien) <- names(TF.motifs.Shrinked.halfsites.count_2)
aa_homoDimer_distance <- c(3, 0, 3, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 5, 6, 0, 0, 0, 0)
names(aa_homoDimer_distance) <- names(TF.motifs.Shrinked.halfsites.count_2)
TF.motifs.Shrinked.halfsites.count_2_maxLLR <- c(6.90008, 11.0531, 6.2838, 10.705, 7.12241, 12.2036,
                                                 10.3949, 16.4493, 12.0708, 6.94399, 8.91828, 11.544,
                                                 3.76517, 6.73403, 7.49103, 7.95261, 7.64841, 7.49922,
                                                 10.0583, 10.7687)
aa_newthr3 <- c(0.17, 0.71, 0.6, 0.65, 0.22,0.48, 0.6, 0.70, 0.65, 0.19, 0.52, 0.81,
                0.14, 0.16, 0.44, 0.48, 0.38, 0.34, 0.57, 0.57)

names(aa_newthr3) <- names(TF.motifs.Shrinked.halfsites.count)
aa <- count_site_from_annotation(annotation_list = aa_GEMSTAT_based_on_linear_exp_annotation_305, 
                                 TF_names = names(TF.motifs.Shrinked.halfsites.count_2),
                                 TF_index=c(1:20),
                                 homoDimer = aa_homodimer, 
                                 dimer_orientation = aa_dimer_orien,
                                 homoDimer_distance=aa_homoDimer_distance,
                                 heterodimer_pair = list(c(13, 14)),
                                 heterodimer_distance = 3,
                                 MAXLLR = TF.motifs.Shrinked.halfsites.count_2_maxLLR,
                                 annotation_thresh = aa_newthr3)
boxplot.matrix(aa$Site_count_mat[[13]], las =2)

####################################################################################################
####################################################################################################
count_site_from_annotation_combined <- function(annotation_list, 
                                                TF_names,
                                                TF_index=c(1:length(TF_names)),
                                                homoDimer, 
                                                dimer_orientation,
                                                homoDimer_distance,
                                                heterodimer_pair,
                                                heterodimer_distance,
                                                MAXLLR,
                                                annotation_thresh,
                                                LLR_to_pVal_list=numeric(0),
                                                pVal_thresh=numeric(0)){
  # This function takes in a list where each entry is an output of annotation_reader function for one
  #  model and counts the number of occurances of each motif for each enhancer, taking homodimers and heterodimers into account
  # annotation_list : is a list where each entry is an output of annotation_reader function
  # TF_names : character vector containing names of all the TFs
  # TF_index : integer vecrtor containing the index of the TFs that we want to count
  # homoDimer : logical vector with length equal to number of TFs, True if
  #  the TF is a homodimer, F if otherwise (monomer or heterodimer)
  # homoDimer_distance : integer vector with length equal to number of TFs. contains 0 for non-homodimer
  #  TFs and a number indicating the number of bases between the two monomers for homodimers
  # heterodimer_pair : is a list where each entry is a tuple contaning the index of two motifs which form
  #  a heterodimer. the first number in each tuple has to be smaller than the second item
  # heterodimer_distance <- is an integer vector with length equal to length of heterodimer_pair
  #  list. showing the number of bases between the two motifs 
  # dimer_orientation is an integer vector, each entry can take one of the four values 0, 1, 2, 3
  #  0 is used for entries which are not homodimers, 1 means both instances on the
  #  same strand ((+,+) or (-,-)), 2 means (+, -), 3 means (-, +)
  # MAXLLR is a numeric vector containing the MAXLLR for each motif. length equal to length of TF_names
  # annotation_thresh is a numeric vector containing the threshold for  (LLR_max - LLR_current)/LLR_max . same as threshold fed to GEMSTAT
  # LLR_to_pVal_list : is a list where each entry is output of LLR2PVAL: a matrix with two columns where,
  #  first column is  LLR*1000 and second column is thr associated p-value. This has to be a named list and the names should correspond to TF names
  # pVal_thresh : is a threshold on p-values, if provided this will be used to compute the
  #  annotation_thresh for each TF and the annotation_thresh input will be overwritten
  
  # checking the sizes of inputs
  stopifnot(length(annotation_list) > 0,
            length(annotation_list[[1]][[1]]) == length(TF_names), 
            all(TF_index <= length(TF_names)) ,
            length(homoDimer) == length(TF_names),
            length(homoDimer_distance) == length(TF_names),
            length(heterodimer_pair) == length(heterodimer_distance),
            length(MAXLLR) == length(TF_names), 
            length(annotation_thresh) == length(TF_names),
            (length(LLR_to_pVal_list) == 0 | length(LLR_to_pVal_list) == length(TF_names)),
            (length(pVal_thresh) == 0 | (length(pVal_thresh) > 0 & length(LLR_to_pVal_list) > 0)))
  
  nu_models <- length(annotation_list)
  nu_genes <- length(annotation_list[[1]])
  gene_names <- names(annotation_list[[1]])
  nu_TFs <- length(TF_names)
  homo_index <- which(homoDimer)
  hetero_first <- integer(0)
  hetero_second <- integer(0)
  for(i in heterodimer_pair){
    stopifnot(i[2] > i[1])
    hetero_first  <- c(hetero_first, i[1])
    hetero_second <- c(hetero_second, i[2])
  }
  if(length(pVal_thresh) > 0){
    stopifnot(length(names(LLR_to_pVal_list)) > 0,
              all(names(LLR_to_pVal_list) == TF_names))
    my_LLR_thresh <- numeric(nu_TFs)
    for(i in 1:nu_TFs){
      aa <- which(LLR_to_pVal_list[[i]][,2] > pVal_thresh)[1]
      my_LLR_thresh[i] <- LLR_to_pVal_list[[i]][(aa-1), 1] / 1000
    }
    names(my_LLR_thresh) <- TF_names
    annotation_thresh <- (MAXLLR - my_LLR_thresh)/MAXLLR
    names(annotation_thresh) <- names(my_LLR_thresh) 
    print("LLR threshold")
    print(my_LLR_thresh)
    print("GEMSTAT annotation_thresh")
    print(annotation_thresh)
  }
  
  
  # create a list of matracies to store the count of TF occurences for each gene
  tf_gene_site_mat <- list()
  for(i in 1:nu_models){
    tf_gene_site_mat[[i]] <- matrix(0L, nrow = nu_genes, ncol = nu_TFs)
    colnames(tf_gene_site_mat[[i]]) <- TF_names
    rownames(tf_gene_site_mat[[i]]) <- gene_names
  }
  
  # getting the occurances of each motif for all genes and all models
  all_occur <- list()
  for(cur_tf in TF_index){
    all_occur[[cur_tf]] <- list()
    for(cur_model in 1:nu_models){
      all_occur[[cur_tf]][[cur_model]] <- lapply(annotation_list[[cur_model]], "[[", cur_tf)
    }
  }
  all_occur_right_pos <- list()
  all_occur_ab_th <- list()
  for(cur_tf in TF_index){
    all_occur_right_pos[[cur_tf]] <- list()
    #print(paste0("TF number ", cur_tf))
    # If homodimer or heterodimer first filter for position, then for LLR
    if(cur_tf %in% c(homo_index, unlist(heterodimer_pair))){ # if TF is homo or heterodimer
      if(cur_tf %in% hetero_first){ # skip if its the first part of a heterodimer
        next()
      }
      if(cur_tf %in% homo_index){ # procesing homodimers
        print(paste0("processing homodimer ", TF_names[cur_tf]))
        
        # creating the dimer oreientation based on input
        if(dimer_orientation[cur_tf] == 1){
          cur_orien <- list(c("+","+"),c("-","-"))
        }else if(dimer_orientation[cur_tf] == 2){
          cur_orien <- list(c("+","-"))
        }else if(dimer_orientation[cur_tf] == 3){
          cur_orien <- list(c("-","+"))
        }
        
        for(cur_model in 1:nu_models){
          all_occur_right_pos[[cur_tf]][[cur_model]] <- list()
          for(cur_gene in 1:nu_genes){
            if(nrow(all_occur[[cur_tf]][[cur_model]][[cur_gene]]) > 1){ # if there is more than one hit for the motif on this gene
              aa_pos <- do.call(rbind ,strsplit(all_occur[[cur_tf]][[cur_model]][[cur_gene]][,1],
                                                split = "\\.."))
              aa_keep <- logical(nrow(aa_pos))
              for(cur_site in 1:(nrow(aa_pos) - 1)){ # loop over sites
                if (((as.integer(aa_pos[cur_site+1, 1]) - as.integer(aa_pos[cur_site, 2])) <= (homoDimer_distance[cur_tf] + 1)) &
                    ((as.integer(aa_pos[cur_site+1, 1]) - as.integer(aa_pos[cur_site, 2])) > 0) &
                    list(c(all_occur[[cur_tf]][[cur_model]][[cur_gene]][cur_site, 2],
                           all_occur[[cur_tf]][[cur_model]][[cur_gene]][cur_site + 1, 2])) %in%  cur_orien &
                    ((all_occur[[cur_tf]][[cur_model]][[cur_gene]]$V4[cur_site] + 
                      all_occur[[cur_tf]][[cur_model]][[cur_gene]]$V4[cur_site + 1]) <= 
                     (2*MAXLLR[cur_tf] - my_LLR_thresh[cur_tf]))){ # if the sites are in the right distance and orientation
                  tf_gene_site_mat[[cur_model]][cur_gene, cur_tf] <- tf_gene_site_mat[[cur_model]][cur_gene, cur_tf] + 1
                  aa_keep[c(cur_site, (cur_site+1))] <- T
                } # if the sites are in the right distance and orientation
              }# end of loop over sites
              all_occur_right_pos[[cur_tf]][[cur_model]][[cur_gene]] <- all_occur[[cur_tf]][[cur_model]][[cur_gene]][aa_keep, ]
            } else if(nrow(all_occur[[cur_tf]][[cur_model]][[cur_gene]]) == 1){ # if there is only one site for a monomer
              all_occur_right_pos[[cur_tf]][[cur_model]][[cur_gene]] <- matrix(nrow = 0, 
                                                                               ncol = ncol(all_occur[[cur_tf]][[cur_model]][[cur_gene]]))
            }
          }# end of loop over genes
          names(all_occur_right_pos[[cur_tf]][[cur_model]]) <- gene_names
        } # end of loop over models
        # end of procesing a homodimer 
      } else if(cur_tf %in% hetero_second){
        # determine the index of the pair 
        cur_tf_het_ind <- which(hetero_second == cur_tf)
        other_tf <- hetero_first[cur_tf_het_ind]
        print(paste0("processing heterodimer", TF_names[cur_tf], " and ", TF_names[other_tf]))
        tf_gene_site_mat[[cur_model]][cur_gene, other_tf] <- 0
        for(cur_model in 1:nu_models){
          all_occur_right_pos[[cur_tf]][[cur_model]] <- list()
          for(cur_gene in 1:nu_genes){
            if(nrow(all_occur[[cur_tf  ]][[cur_model]][[cur_gene]]) > 0 &
               nrow(all_occur[[other_tf]][[cur_model]][[cur_gene]]) > 0){# if both motifs of heterodimer have sites
              aa_pos_1 <- do.call(rbind ,strsplit(all_occur[[other_tf]][[cur_model]][[cur_gene]][,1],
                                                  split = "\\.."))
              aa_pos_2 <- do.call(rbind ,strsplit(all_occur[[cur_tf]][[cur_model]][[cur_gene]][,1],
                                                  split = "\\.."))
              aa_keep_1 <- logical(nrow(aa_pos_1))
              aa_keep_2 <- logical(nrow(aa_pos_2))
              for(cur_site_1 in 1:nrow(aa_pos_1)){
                for(cur_site_2 in 1:nrow(aa_pos_2)){
                  if( (as.integer(aa_pos_2[cur_site_2, 1]) - as.integer(aa_pos_1[cur_site_1, 2])) <= (heterodimer_distance[cur_tf_het_ind] + 1) &
                      (as.integer(aa_pos_2[cur_site_2, 1]) - as.integer(aa_pos_1[cur_site_1, 2])) > 0 & 
                      (all_occur[[cur_tf]][[cur_model]][[cur_gene]][cur_site_2, 2] == 
                       all_occur[[other_tf]][[cur_model]][[cur_gene]][cur_site_1, 2]) & 
                      ((all_occur[[cur_tf]][[cur_model]][[cur_gene]]$V4[cur_site_2] + 
                        all_occur[[other_tf]][[cur_model]][[cur_gene]]$V4[cur_site_1]) <= 
                       (MAXLLR[cur_tf] + MAXLLR[other_tf] - my_LLR_thresh[cur_tf]))){
                    tf_gene_site_mat[[cur_model]][cur_gene, cur_tf] <- tf_gene_site_mat[[cur_model]][cur_gene, cur_tf] + 1
                    tf_gene_site_mat[[cur_model]][cur_gene, other_tf] <- tf_gene_site_mat[[cur_model]][cur_gene, other_tf] + 1
                    aa_keep_1[cur_site_1] <- T
                    aa_keep_2[cur_site_2] <- T
                  } # if they are in the right distance and orientation
                } # end  of loop over sites of second half
              } # end  of loop over sites of first half
              all_occur_right_pos[[cur_tf  ]][[cur_model]][[cur_gene]] <- all_occur[[cur_tf  ]][[cur_model]][[cur_gene]][aa_keep_2,]
              all_occur_right_pos[[other_tf]][[cur_model]][[cur_gene]] <- all_occur[[other_tf]][[cur_model]][[cur_gene]][aa_keep_1,]
            }else{ # if both motifs of heterodimer do not have sites
              all_occur_right_pos[[cur_tf  ]][[cur_model]][[cur_gene]] <- matrix(nrow = 0, ncol=  ncol(all_occur[[cur_tf  ]][[cur_model]][[cur_gene]]))
              all_occur_right_pos[[other_tf]][[cur_model]][[cur_gene]] <- matrix(nrow = 0, ncol=  ncol(all_occur[[other_tf]][[cur_model]][[cur_gene]]))
            }
          } # end of loop over genes
          names(all_occur_right_pos[[cur_tf]][[cur_model]]) <- gene_names
          names(all_occur_right_pos[[other_tf]][[cur_model]]) <- gene_names
        } # end of loop over models
        # end of processing heterodimers
      }
      
    } else{ # if TF is not a homo or heterodimer
      for(cur_model in 1:nu_models){
        all_occur_right_pos[[cur_tf]][[cur_model]] <- list()
        for(cur_gene in 1:nu_genes){
          all_occur_right_pos[[cur_tf]][[cur_model]][[cur_gene]] <- all_occur[[cur_tf]][[cur_model]][[cur_gene]][all_occur[[cur_tf]][[cur_model]][[cur_gene]]$V4 <= annotation_thresh[cur_tf]*MAXLLR[cur_tf], ]
          tf_gene_site_mat[[cur_model]][cur_gene, cur_tf] <- nrow(all_occur_right_pos[[cur_tf]][[cur_model]][[cur_gene]])
        } # end of loop over genes
        names(all_occur_right_pos[[cur_tf]][[cur_model]]) <- gene_names
      } # end of loop over models
    } # end of processing one monomer
  }# end of loop over TFs
  
  # TODO : turn the new annotations into an annotation file readable by GEMSTAT
  all_annot_per_model_per_gene <- list()
  for(cur_model in 1:nu_models){
    all_annot_per_model_per_gene[[cur_model]] <- list()
    for(cur_gene in 1:nu_genes){
      all_annot_per_model_per_gene[[cur_model]][[cur_gene]] <- matrix(nrow = 0, 
                                                                      ncol = ncol(all_occur_right_pos[[1]][[1]][[1]]))
      for(cur_tf in 1:nu_TFs){
        all_annot_per_model_per_gene[[cur_model]][[cur_gene]] <- rbind(all_annot_per_model_per_gene[[cur_model]][[cur_gene]],
                                                                       all_occur_right_pos[[cur_tf]][[cur_model]][[cur_gene]])
      }
      my_ind <- sort(as.integer(rownames(all_annot_per_model_per_gene[[cur_model]][[cur_gene]])), decreasing = F, index.return=T)$ix
      all_annot_per_model_per_gene[[cur_model]][[cur_gene]] <- all_annot_per_model_per_gene[[cur_model]][[cur_gene]][my_ind, ]
    } # end of loop over genes
    names(all_annot_per_model_per_gene[[cur_model]]) <- gene_names
  } # end of loop over models

  
  
  return(list(Sites_Per_TF = all_occur_right_pos,
              Site_count_mat = tf_gene_site_mat,
              Annot_per_model_perGene = all_annot_per_model_per_gene))
  
}


####################################################################################################
####################################################################################################
 # TODO : test an example
# example
aa_homodimer <-   c(T, F, T, F, F, F, F, F, T, F, F, F, T, T, F, T, F, F)
names(aa_homodimer) <- names(TF.motifs.Shrinked.hocomoco.count)
aa_dimer_orien <-        c(2, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 0, 1, 0, 0)
names(aa_dimer_orien) <- names(TF.motifs.Shrinked.hocomoco.count)
aa_homoDimer_distance <- c(3, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 5, 5, 0, 1, 0, 0)
names(aa_homoDimer_distance) <- names(TF.motifs.Shrinked.hocomoco.count)
TF.motifs.Shrinked.hocomoco.count_maxLLR <- c(6.13309, 10.0693, 6.28522, 11.3454, 7.91607, 9.98706,
                                                 8.40569, 12.5907, 7.2447, 10.2185, 13.5565, 11.9268,
                                                 6.69131, 7.37613, 9.94809, 6.93567, 12.5285, 9.73563)
names(TF.motifs.Shrinked.hocomoco.count_maxLLR) <- names(TF.motifs.Shrinked.hocomoco.count)
aa_newthr3 <- c(0.17, 0.71, 0.6, 0.65, 0.22,0.48, 0.6, 0.70, 0.65, 0.19, 0.52, 0.81,
                0.14, 0.16, 0.44, 0.48, 0.38, 0.34)

names(aa_newthr3) <- names(TF.motifs.Shrinked.hocomoco.count)
aa <- count_site_from_annotation_combined(annotation_list = GEMSTAT_based_on_linear_exp_annotation_1397, 
                                          TF_names = names(TF.motifs.Shrinked.hocomoco.count),
                                          TF_index=c(1:18),
                                          homoDimer = aa_homodimer, 
                                          dimer_orientation = aa_dimer_orien,
                                          homoDimer_distance=aa_homoDimer_distance,
                                          heterodimer_pair = numeric(0),
                                          heterodimer_distance = numeric(0),
                                          MAXLLR = TF.motifs.Shrinked.hocomoco.count_maxLLR,
                                          annotation_thresh = aa_newthr3,
                                          LLR_to_pVal_list = TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval,
                                          pVal_thresh = 0.0001)
GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m4 <- aa
#boxplot.matrix(aa$Site_count_mat[[13]], las = 2)
aa <- count_site_from_annotation_combined(annotation_list = GEMSTAT_based_on_linear_exp_annotation_1397, 
                                          TF_names = names(TF.motifs.Shrinked.hocomoco.count),
                                          TF_index=c(1:18),
                                          homoDimer = aa_homodimer, 
                                          dimer_orientation = aa_dimer_orien,
                                          homoDimer_distance=aa_homoDimer_distance,
                                          heterodimer_pair = numeric(0),
                                          heterodimer_distance = numeric(0),
                                          MAXLLR = TF.motifs.Shrinked.hocomoco.count_maxLLR,
                                          annotation_thresh = aa_newthr3,
                                          LLR_to_pVal_list = TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval,
                                          pVal_thresh = 0.001)
GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3 <- aa


for(i in 1:length(TF.motifs.Shrinked.hocomoco)){
  seqLogo::seqLogo(t(TF.motifs.Shrinked.hocomoco[[i]]))
}
boxplot.matrix(aa$Site_count_mat[[3]], las =2, outline = T)
sum(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat[[13]][, 3] > 0)
sum(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m4$Site_count_mat[[13]][, 3] > 0)

j=0
for(i in 1:148){
  print(i)
  j = j + sum(colSums(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m4$Site_count_mat[[i]]) == 0)
  print(sum(colSums(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m4$Site_count_mat[[i]]) == 0))
}

####################################################################################################
####################################################################################################
annotation_Writer <- function(annotation_list, f_name=character(0)){
  # annotation_list is the third entry of the output of "count_site_from_annotation_combined" function.
  # f_name: name of the filtered annot file
  nu_models <- length(annotation_list)
  for(cur_model in 1:nu_models){
    nu_genes <- length(annotation_list[[cur_model]])
    for(cur_gene in 1:nu_genes){
      if(length(f_name) > 0){
        aaf_name <- f_name
      }else{
        aaf_name <-  paste0("annot_", cur_model, ".ann")
      }
      cat(paste0(">", names(annotation_list[[cur_model]])[cur_gene], "\n"), 
          file = aaf_name, append = T)
      write.table(as.matrix(annotation_list[[cur_model]][[cur_gene]]),
                  file = aaf_name, 
                  quote = F, row.names = F, col.names = F, sep = "\t", append = T)
      
    }
  }
}
####################################################################################################
####################################################################################################
# example
setwd("Seeded_GEMSTAT_ens/Experiment_1397/Annotation/Filtered/")
annotation_Writer(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Annot_per_model_perGene)
setwd("../../../..")

####################################################################################################
####################################################################################################
KD_Performer <- function(parent_Exp, 
                         TF_index,
                         model_index=integer(0),
                         shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                         prev_na,
                         run_GEMSTAT=F){
  # parent_Exp : is an integer indicating the name of the Experiment where KDs are going to be performed on
  # TF_index: is a list where each entry is an integer vector of the indecis of TFs that
  #  you want to Knockdown in that experiment
  # shared_dir : is the directori where the experiment leaves
  # number of na in parent experiment
  # run_GEMSTAT : logical, if True runs GEMSTAT
  # model_index : is an integer vector, containing the number of models to perform the KD on. if nothing is provided all models are considered
  outside_directory <- getwd()
  on.exit(setwd(outside_directory))
  # go to the directory of the parent Experiment
  Experiment_directory <- paste0(shared_dir, "/Experiment_", parent_Exp)
  setwd(Experiment_directory)
  # create a folder for KD Exps
  if (! dir.exists("KnockDown")){
    system("mkdir KnockDown")
  }
  KD_directory <- paste0(shared_dir, "/Experiment_", parent_Exp,"/KnockDown")
  # create a folder for KD TF expression
  setwd(KD_directory)
  #############################################################
  # Copy parameters from output of the parent experiment
  if (! dir.exists("Parameters")){
    system("mkdir Parameters")
  }
  KD_par_dir <-  paste0(KD_directory, "/Parameters")
  WT_Output_dir <- paste0(Experiment_directory, "/Outputs")
  print(paste0("copying params from outputs of experiment ", parent_Exp))
  setwd(WT_Output_dir)
  list_of_files <- list.files(WT_Output_dir, pattern = "*.Filter")
  stopifnot(dir.exists(KD_par_dir), 
            dir.exists(WT_Output_dir))
  file.copy(from = list_of_files, to = KD_par_dir)
  setwd(KD_par_dir)
  # rename files
  aaa <- strsplit(list_of_files, split = "_")
  aaaa <- lapply(aaa, unlist)
  aaaaa <- lapply(aaaa, strsplit, split="\\.")
  aaaaaa <- lapply(aaaaa, unlist)
  aaaaaaa <- unlist(lapply(aaaaaa, "[[", 3))
  file.rename(list_of_files, paste0("start_", aaaaaaa, ".par"))
  #############################################################
  # create TF exp matrix and job file for each KD experiment
  
  for (cur_kd in 1:length(TF_index)){
    setwd(KD_directory)
    if(dir.exists(paste(TF_index[[cur_kd]], collapse = "_"))){
      next()
    }
    system(paste0("mkdir ", paste(TF_index[[cur_kd]], collapse = "_")))
    cur_KD_directory <- paste0(KD_directory, "/", paste(TF_index[[cur_kd]], collapse = "_"))
    setwd(cur_KD_directory)
    system("mkdir TF_Expression")
    system("mkdir Outputs")
    # create KD TF expression matrix for each TF in TF_index
    setwd(paste0(cur_KD_directory, "/TF_Expression"))
    WT_TF_expression <- read.delim(file = paste0(Experiment_directory, "/Inputs/TF_Expression/TF_exp.tab"),
                                   sep = "\t", header = F, stringsAsFactors = F )
    all_na <- which(colSums(is.na(WT_TF_expression)) == nrow(WT_TF_expression))
    if(length(all_na) > 0){
      WT_TF_expression <- WT_TF_expression[, -all_na]
    }
   
    colnames(WT_TF_expression) <- WT_TF_expression[1,]
    rownames(WT_TF_expression) <- WT_TF_expression[, 1]
    WT_TF_expression <- WT_TF_expression[2:nrow(WT_TF_expression), 2:ncol(WT_TF_expression)]
    
    my_TF_expMat <- WT_TF_expression
    my_TF_expMat[TF_index[[cur_kd]], ] <- 0
    ExpressionWriter(expressionMat = my_TF_expMat,
                     output.File.Name = "TF_exp")
    # Copying the job file and fixing it
    setwd(Experiment_directory)
    file.copy(from = list.files(Experiment_directory, pattern = "*.job") , to = cur_KD_directory)
    setwd(cur_KD_directory)
    # writing the replace commands
    #gemstat dir
    replace_command <- paste0("sed -i ",'"" ' ,"'s#",
                              "/shared-mounts/sinhas/tabebor2/GEMSTAT_git/my_fork/GEMSTAT/src/seq2expr",
                              "#",
                              "/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr",
                              "#g' ",
                              "based_on_linear_1enh_perGene.job")
    system(replace_command)
    # shared dir
    replace_command <- paste0("sed -i ",'"" ' ,"'s#",
                              "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/based_on_linear/Experiment_", parent_Exp,
                              "#",
                              Experiment_directory,
                              "#g' ",
                              "based_on_linear_1enh_perGene.job")
    system(replace_command)
    # change TF_expression file
    replace_command <- paste0("sed -i ",'"" ' ,"'s#",
                              "/Inputs/TF_Expression/TF_exp.tab",
                              "#",
                              "/KnockDown/", paste(TF_index[[cur_kd]], collapse = "_") ,"/TF_Expression", "/TF_exp.tab",
                              "#g' ",
                              "based_on_linear_1enh_perGene.job")
    system(replace_command)
    # use the output parameters
    replace_command <- paste0("sed -i ",'"" ' ,"'s#",
                              "Inputs/Parameters",
                              "#",
                              "KnockDown/Parameters",
                              "#g' ",
                              "based_on_linear_1enh_perGene.job")
    system(replace_command)
    # change the na to zero 
    replace_command <- paste0("sed -i ",'"" ' ,"'s#",
                              "-na ", prev_na, "#",
                              "-na ", 0,
                              "#g' ",
                              "based_on_linear_1enh_perGene.job")
    system(replace_command)
    
    # outputs
    replace_command <- paste0("sed -i ",'"" ' ,"'s#",
                              "Outputs",
                              "#",
                              "KnockDown/", paste(TF_index[[cur_kd]], collapse = "_"), "/Outputs",
                              "#g' ",
                              "based_on_linear_1enh_perGene.job")
    system(replace_command)
    add_bash_command <- paste0("sed -i ", '"" ', "'1s:^:#!/bin/bash\\'$'\\n:' based_on_linear_1enh_perGene.job")
    system(add_bash_command)
    # remove lines corresponding to models that are not wanted
    if(length(model_index) > 0){
      to_be_Removed <- sort(setdiff(c(1:length(list_of_files)), 
                                    model_index),
                            decreasing = T)
      for(cur_rm in 1:length(to_be_Removed)){
        rm_command <- paste0('sed -i "" ', "'", to_be_Removed[cur_rm] + 1, "d' based_on_linear_1enh_perGene.job")
        system(rm_command)
      }
    }
    executable_command <- "chmod +x based_on_linear_1enh_perGene.job"
    system(executable_command)
    if(run_GEMSTAT){
      run_command <- "./based_on_linear_1enh_perGene.job"
      system(run_command)
    }
  }
}
########################################################################################################################
########################################################################################################################
# example
KD_Performer(parent_Exp=758, 
             TF_index=list(1, 2, c(4, 5)),
             shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
             prev_na=5, run_GEMSTAT = T)
########################################################################################################################
########################################################################################################################
read_KD_results <- function(experiment_nu, 
                            experiment_output, 
                            shared_dir, 
                            KD_name = character(0),
                            coop_KD=F){
  # uses these functions: 
  # GEMSTAT_output_Reader_multiEnh_ensemble, 
  # experiment_nu : is an integer showing the number of the experiment
  # experiment_output: output of GEMSTAT_output_Reader_multiEnh_ensemble for the parent experiment
  # shared_dir : is the directory where all experiments live
  # KD_name : the name of the KD to be read. if provided it only reads the results of this KD and 
  #  not the rest, if not provided it reads all KDs in the folder
  # coop_KD : logical, read coop KD results or normal KD results
  add_slash <- ifelse(substr(shared_dir, nchar(shared_dir), 
                             nchar(shared_dir)) == "/", yes = "", no = "/")
  parent_dir <- paste0(shared_dir, add_slash, "Experiment_", experiment_nu)
  parent_content <- list.dirs(paste0(shared_dir, add_slash, "Experiment_", experiment_nu),
                              recursive = F, full.names = F)
  stopifnot("KnockDown" %in% parent_content)
  KD_dir <- paste0(parent_dir, "/", "KnockDown")
  KD_content <- list.dirs(path = KD_dir,
                          full.names = F,
                          recursive = F)
  if(coop_KD){ # reading coop KD results
    stopifnot("Coop_KD" %in% KD_content)
    coop_KD_dir <- paste0(KD_dir, "/", "Coop_KD")
    coop_KD_content <- list.dirs(path = coop_KD_dir,
                                 full.names = F,
                                 recursive = F)
    if(length(KD_name) > 0){
      stopifnot(all(KD_name %in% coop_KD_content))
      coop_KD_content <- KD_name
    }
    KF_output <- list()
    for(cur_KD in 1:length(coop_KD_content)){
      cur_KD_dir <- paste0(coop_KD_dir, "/", coop_KD_content[cur_KD])
      KF_output[[cur_KD]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0(cur_KD_dir,
                                                                                         "/Outputs"),
                                                                     .fold_change = T,
                                                                     ..fixed_thresh = experiment_output$threshold)
    }
    names(KF_output) <- coop_KD_content
    return(KF_output)
  }else{ # if not reading the coop KD results (reading TF KD results)
    KD_content <- setdiff(KD_content, "Parameters")
    
    un_desired <- unlist(lapply(lapply(lapply(X = strsplit(KD_content, split="_"),
                                              as.numeric), is.na), sum)) > 0
    if(sum(un_desired) > 0){
      print(paste0("not considering the following entries as KD items: ",
                   paste0(KD_content[un_desired], collapse = " ")))
    }
    KD_content <- KD_content[!un_desired]
    if(length(KD_name) > 0){
      stopifnot(all(KD_name %in% KD_content))
      KD_content <- KD_name
    }
    
    # stop if there is any irrelevant file or directory
    stopifnot(! any(is.na(as.numeric(unlist(strsplit(KD_content, split="_"))))))
    stopifnot(length(experiment_output$threshold) > 0) # make sure a threshold is set
    KF_output <- list()
    for(cur_KD in 1:length(KD_content)){
      cur_KD_dir <- paste0(KD_dir, "/", KD_content[cur_KD])
      KF_output[[cur_KD]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0(cur_KD_dir,
                                                                                         "/Outputs"),
                                                                     .fold_change = T, 
                                                                     ..fixed_thresh = experiment_output$threshold)
      
    }
    names(KF_output) <- KD_content
    return(KF_output)
  }
}
########################################################################################################################
########################################################################################################################
#example
aa <- read_KD_results(experiment_nu = 758, 
                      shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens")
########################################################################################################################
########################################################################################################################
KD_Performer_by_params <- function(parent_Exp,
                                   TF_names,
                                   Coop_mat_list, 
                                   shared_dir, 
                                   prev_na,
                                   run_GEMSTAT){
  # This function Performs KDs on cooperativity (for now, later other changes can be added as well)
  # parent_Exp : is an integer indicating the name of the Experiment where KDs are going to be performed on
  # TF_names : is the name of the TFs
  # Coop_mat_list: is a list where each entry is a matrix, containing the indices of the
  #  cooperating TFs. each column of the matrix indicates a TF and each row is an interaction to be removed 
  # shared_dir : is the directory where the experiment leaves
  # prev_na : nu_nanumber of na in parent experiment
  # run_GEMSTAT : logical, if True runs GEMSTAT
  
  library(rjson)
  outside_directory <- getwd()
  on.exit(setwd(outside_directory))
  # go to the directory of the parent Experiment
  Experiment_directory <- paste0(shared_dir, "/Experiment_", parent_Exp)
  setwd(Experiment_directory)
  
  # create a folder for KD Exps
  if (! dir.exists("KnockDown")){
    system("mkdir KnockDown")
  }
  KD_directory <- paste0(shared_dir, "/Experiment_", parent_Exp,"/KnockDown")
  setwd(KD_directory)
  
  # create a folder for Coop KDs
  if (! dir.exists("Coop_KD")){
    system("mkdir Coop_KD")
  }
  Coop_KD_directory <- paste0(KD_directory, "/", "Coop_KD")
  setwd(Coop_KD_directory)
  # Go through the list of provided coops
  for(cur_coop_mat in 1:length(Coop_mat_list)){
    # if there is only one row and its not formatted as a matrix, reformat it
    if(is.null(nrow(Coop_mat_list[[cur_coop_mat]]))){
      Coop_mat_list[[cur_coop_mat]] <- matrix(Coop_mat_list[[cur_coop_mat]], nrow=1, byrow = T)
    }
    # create a folder for the current coop KD
    cur_folder_name <- character(0)
    for(i in 1:nrow(Coop_mat_list[[cur_coop_mat]])){
      cur_folder_name <- paste0(cur_folder_name,
                                TF_names[Coop_mat_list[[cur_coop_mat]][i, 1]],
                                "_",
                                TF_names[Coop_mat_list[[cur_coop_mat]][i, 2]],
                                "__")
    }
    cur_folder_name <- substr(cur_folder_name, 1, (nchar(cur_folder_name)-2))
    cur_coop_KD_dir <- paste0(Coop_KD_directory, "/", cur_folder_name)
    if(! dir.exists(cur_coop_KD_dir)){
      system(paste0("mkdir ", cur_coop_KD_dir))
    }
    setwd(cur_coop_KD_dir)
    # create a parameter and an output directory
    if(! dir.exists("Parameters")){
      system("mkdir Parameters")
    }
    cur_coop_KD_par_dir <- paste0(cur_coop_KD_dir, "/", "Parameters")
    if(! dir.exists("Outputs")){
      system("mkdir Outputs")
    }
    cur_coop_KD_output_dir <- paste0(cur_coop_KD_dir, "/", "Outputs")
    
    
    # read the output parameter files for each model
    WT_Output_dir <- paste0(Experiment_directory, "/Outputs")
    print(paste0("reading params from outputs of experiment ", parent_Exp))
    setwd(WT_Output_dir)
    list_of_files <- list.files(WT_Output_dir, pattern = "*.Filter")
    
    ### renaming the parameters
    aaa <- strsplit(list_of_files, split = "_")
    aaaa <- lapply(aaa, unlist)
    aaaaa <- lapply(aaaa, strsplit, split="\\.")
    aaaaaa <- lapply(aaaaa, unlist)
    aaaaaaa <- unlist(lapply(aaaaaa, "[[", 3))
    new_par_names <- paste0("start_", aaaaaaa, ".par")
    ###
    for(cur_model in 1:length(list_of_files)){
      print(list_of_files[cur_model])
      cur_par_file <- fromJSON(file = list_of_files[cur_model])
      if(typeof(cur_par_file$qbtm) != "list"){
        cur_par_file$qbtm <- list(cur_par_file$qbtm)
      }
      all_coops <- names(cur_par_file$inter)
      # go through each pair
      for(cur_row in 1:nrow(Coop_mat_list[[cur_coop_mat]])){
        ref_name <- paste(TF_names[Coop_mat_list[[cur_coop_mat]][cur_row, ]], collapse = ":")
        alt_name <- paste(TF_names[rev(Coop_mat_list[[cur_coop_mat]][cur_row, ])], collapse = ":")
        cur_existance <- c(ref_name, alt_name) %in% all_coops
        stopifnot(sum(cur_existance) == 1 | (sum(cur_existance) == 2 & ref_name==alt_name))
        cur_par_file$inter[[which(all_coops %in% c(ref_name, alt_name))]] <- 1
      } # end of loop over rows of one KD list entry
      cat(toJSON(cur_par_file), append = F, 
          file = paste0(cur_coop_KD_par_dir, "/", new_par_names[cur_model]))
    } # end of loop over models
    # Copying the job file and fixing it
    setwd(Experiment_directory)
    file.copy(from = list.files(Experiment_directory, pattern = "*.job") , to = cur_coop_KD_dir)
    setwd(cur_coop_KD_dir)
    # writing the replace commands
    #gemstat dir
    replace_command <- paste0("sed -i ",'"" ' ,"'s#",
                              "/shared-mounts/sinhas/tabebor2/GEMSTAT_git/my_fork/GEMSTAT/src/seq2expr",
                              "#",
                              "/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr",
                              "#g' ",
                              "based_on_linear_1enh_perGene.job")
    system(replace_command)
    # shared dir
    replace_command <- paste0("sed -i ",'"" ' ,"'s#",
                              "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/based_on_linear/Experiment_", parent_Exp,
                              "#",
                              Experiment_directory,
                              "#g' ",
                              "based_on_linear_1enh_perGene.job")
    system(replace_command)
    # use the output parameters
    replace_command <- paste0("sed -i ",'"" ' ,"'s#",
                              "Inputs/Parameters",
                              "#",
                              "KnockDown/Coop_KD/",cur_folder_name,"/Parameters",
                              "#g' ",
                              "based_on_linear_1enh_perGene.job")
    system(replace_command)
    # change the na to zero 
    replace_command <- paste0("sed -i ",'"" ' ,"'s#",
                              "-na ", prev_na, "#",
                              "-na ", 0,
                              "#g' ",
                              "based_on_linear_1enh_perGene.job")
    system(replace_command)
    # outputs
    replace_command <- paste0("sed -i ",'"" ' ,"'s#",
                              "Outputs",
                              "#",
                              "KnockDown/Coop_KD/", cur_folder_name, "/Outputs",
                              "#g' ",
                              "based_on_linear_1enh_perGene.job")
    system(replace_command)
    add_bash_command <- paste0("sed -i ", '"" ', "'1s:^:#!/bin/bash\\'$'\\n:' based_on_linear_1enh_perGene.job")
    system(add_bash_command)
    executable_command <- "chmod +x based_on_linear_1enh_perGene.job"
    system(executable_command)
    if(run_GEMSTAT){
      run_command <- "./based_on_linear_1enh_perGene.job"
      system(run_command)
    }
  } # end of loop over KD list entries
}
########################################################################################################################
########################################################################################################################
# example
KD_Performer_by_params(parent_Exp=777,
                       TF_names = names(TF.motifs.Shrinked.halfsites.count_2),
                       Coop_mat_list=list(c(3, 3), c(10, 10)), 
                       shared_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens", 
                       prev_na=5,
                       run_GEMSTAT=T)
aa_coop_kd <- read_KD_results(experiment_nu=777, 
                              shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens", 
                              coop_KD=T)

########################################################################################################################
########################################################################################################################
GeneFoldChange_by_TFFoldChange <- function(Gene_exp_mat,
                                           TF_exp_mat,
                                           TF_gene_site_mat,
                                           site_nu_thresh=1,
                                           TF_index=c(1:nrow(TF_exp_mat))){
  # Gene_exp_mat : is the gene expression matrix, has one row per gene and one column per pair of
  #  control and treatment conditions
  # TF_Exp_mat : is the TF expression matrix, has one row per TF, but it can have one column per pair of
  #  Treatment and Control conditions, or one column for each condition in which case the fold changes
  #  will be calculated. if one column per conditions, then the treatment of each expeiment should follow the control
  # TF_gene_site_mat : is a matrix which contains the number of sites of each TF in each gene's enhancer. each
  #  Row indicates a genes and each column inidicates a TF. Rows and columns have to be named
  # site_nu_thresh : minimum number of sites of a TF within an enhancer, in order for the enhancer to be considred a
  #  target for that TF
  # TF_index : is an integer vector containing the index of the TFs to be plotted.
  
  # This function plots a boxplot for each TF in TF_index. each plot will have a separate box for each
  #  fold change value of the TF. This box contains the Fold change value of the genes with at least
  #  site_nu_thresh number of sites for the TF, in the conditions that the TF had an specific fold change value.
  
  # check the size of inputs
  stopifnot(nrow(Gene_exp_mat) == nrow(TF_gene_site_mat),
            (ncol(Gene_exp_mat) == ncol(TF_exp_mat) | (2*ncol(Gene_exp_mat)) == ncol(TF_exp_mat)),
            nrow(TF_exp_mat) == ncol(TF_gene_site_mat),
            site_nu_thresh > 0,
            length(TF_index) > 0,
            length(rownames(TF_gene_site_mat)) == nrow(TF_gene_site_mat), 
            length(colnames(TF_gene_site_mat)) == ncol(TF_gene_site_mat))
  
  nu_genes <- nrow(Gene_exp_mat)
  gene_names <- rownames(Gene_exp_mat)
  nu_TFs <- nrow(TF_exp_mat)
  TF_names <- rownames(TF_exp_mat)
  
  # compute the TF log fold change if not given
  if((2*ncol(Gene_exp_mat)) == ncol(TF_exp_mat)){ 
    TF_Exp_FC <- matrix(nrow = nrow(TF_exp_mat),
                        ncol = ncol(Gene_exp_mat))
    rownames(TF_Exp_FC) <- rownames(TF_exp_mat)
    colnames(TF_Exp_FC) <- colnames(Gene_exp_mat)
    for(i in 1:ncol(Gene_exp_mat)){
      TF_Exp_FC[,i] <- log2((TF_exp_mat[,2*i] + 0.001)/(TF_exp_mat[,2*i - 1]+ 0.001))
    }
    TF_Exp_FC <- TF_Exp_FC[TF_index, ]
  } else{
    TF_Exp_FC <- TF_exp_mat[TF_index, ]
  }
  
  # trim the TF_gene_site_mat for the chosen TFs
  TF_gene_site_mat_chosen <- TF_gene_site_mat[, TF_index]
  
  # compute the index of genes with at least site_nu_thresh number of sites for the TFs specidied in TF_index
  TF_target_gene_index <- list()
  for(i in 1:length(TF_index)){
    TF_target_gene_index[[i]] <- which(TF_gene_site_mat_chosen[, i] >= site_nu_thresh)
  }
  names(TF_target_gene_index) <- TF_names[TF_index]
  
  # plot 
  prev_par <- par()
  on.exit(par(mfrow = prev_par$mfrow,
              mar = prev_par$mar))
  par(mfrow = set_row_col_nu(length(TF_index)),
      mar = c(3, 3, 3, 3))
  for(cur_tf in 1:length(TF_index)){
    cur_tf_table <- table(TF_Exp_FC[cur_tf, ])
    cur_tf_unique <- sort(unique(TF_Exp_FC[cur_tf, ]))
    stopifnot(all.equal(cur_tf_unique,
                        as.numeric(names(cur_tf_table))))
    cond_all <- integer(0)
    exp_all <- numeric(0)
    name_all <- character(0)
    for(cur_tab in 1:length(cur_tf_table)){
      cur_conds <- which(TF_Exp_FC[cur_tf, ] %in% cur_tf_unique[cur_tab])
      stopifnot(length(cur_conds) == cur_tf_table[cur_tab])
      cond_all <- c(cond_all, rep(cur_tab,
                                  (length(TF_target_gene_index[[cur_tf]]) * cur_tf_table[cur_tab])))
      name_all <- c(name_all, rep(format(as.numeric(names(cur_tf_table)[cur_tab]), digits = 3), 
                                  (length(TF_target_gene_index[[cur_tf]]) * cur_tf_table[cur_tab])))

      exp_all <- c(exp_all, as.numeric(Gene_exp_mat[TF_target_gene_index[[cur_tf]],
                                                    cur_conds]))
    }
    boxplot(exp_all~name_all, 
            main = paste(TF_names[TF_index[cur_tf]], 
                         length(TF_target_gene_index[[cur_tf]]),
                         paste(cur_tf_table, collapse = ","),
                         sep = "__"),
            las = 2)
  }
}
########################################################################################################################
########################################################################################################################
# example
my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_48_gene_gte3WTGRcommon
# TF expression mat is RAP specific
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_ER01_monomer_hocomoco
# TF_gene_site_mat is model speicific
GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat[[1]]


GeneFoldChange_by_TFFoldChange(Gene_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_48_gene_gte3WTGRcommon,
                                           TF_exp_mat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_ER01_monomer_hocomoco,
                                           TF_gene_site_mat = GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat[[1]],
                                           site_nu_thresh=1)
########################################################################################################################
########################################################################################################################
set_row_col_nu <- function(nu_all, more_rows = T){
  # returns the number of rows and columns for a plot given total nu_all number of subplots
  # nu_all : an integer, indicating the number of subplots
  # more_rows : is a logical if True the number of rows will be greater than the number of columns
  aa_sq <- sqrt(nu_all)
  if((floor(aa_sq) * ceiling(aa_sq)) >= nu_all){
    if(more_rows){
      aa_nrow <- ceiling(aa_sq)
      aa_ncol <- floor(aa_sq)
    }else{
      aa_nrow <- floor(aa_sq)
      aa_ncol <- ceiling(aa_sq)
    }
  }else{
    aa_nrow <- ceiling(aa_sq)
    aa_ncol <- ceiling(aa_sq)
  }
  return(c(aa_nrow, aa_ncol))
}
########################################################################################################################
########################################################################################################################
# example
set_row_col_nu(88)
########################################################################################################################
########################################################################################################################
GeneFoldChange_by_TFFoldChange_model_RAP_wrapper <- function(.Gene_exp_mat,
                                                             .TF_exp_mat,
                                                             TF_gene_site_mat_list, 
                                                             .site_nu_thresh,
                                                             .TF_index = c(1:nrow(.TF_exp_mat)),
                                                             RAP_raw_mat=numeric(0), 
                                                             RAP_list = list(),
                                                             RAP_nu=integer(0),
                                                             model_nu){
  # .Gene_exp_mat : is the gene expression matrix, has one row per gene and one column per pair of
  #  control and treatment conditions
  # .TF_Exp_mat : is the TF expression matrix, has one row per TF, but it can have one column per pair of
  #  Treatment and Control conditions, or one column for each condition in which case the fold changes
  #  will be calculated. if one column per conditions, then the treatment of each expeiment should follow the control
  # TF_gene_site_mat_list : is a list where each entry corresponds to a model. each entry is a matrix which contains the number of sites of each TF in each gene's enhancer. each
  #  Row indicates a genes and each column inidicates a TF. Rows and columns have to be named
  # .site_nu_thresh : minimum number of sites of a TF within an enhancer, in order for the enhancer to be considred a
  #  target for that TF
  # .TF_index : is an integer vector containing the index of the TFs to be plotted.
  # RAP_raw_mat: is the matrix where each row corresponds to a RAP and each column corresponds to a nuclear receptor.
  #  columns are named after TFs. each enry indicates the receptor activation profile of a certain TF and it can take values of (0, 1, 2)
  # RAP_list:  is a list where each entry is a list and corresponds to a TF. the entries of the list
  #  for each TF are the options that that TF can take in the combinations: here is an example:
  # aa_modif_all_options_p20_1398 <- list(list(1, -1, 21) , # this is for TF number 1. its saying that TF can take "0,1" , "1,0", and "1,1" configurations in "control, treatment" pairs
  #                                       list(9, -9, 29),
  #                                       list(12, -12, 32),
  #                                       list(13,-13, 33),
  #                                       list(14,-14, 34),
  #                                       list(16,-16, 36))
  # RAP_nu : is an integer indicating the row of RAP_raw_mat to be considered.
  # model_nu : is an integer indicating the index of the model to be used.
  
  # checking the size of inputs
  stopifnot(model_nu <= length(TF_gene_site_mat_list))
  
  # if RAPs are used, then update the TF expression matrix
  if(length(RAP_raw_mat) > 0){
    RAP_list_dif <- max(unlist(RAP_list[[1]])) - abs(min(unlist(RAP_list[[1]])))
    stopifnot(all(colnames(RAP_raw_mat) %in% rownames(.TF_exp_mat)), 
              RAP_nu < nrow(RAP_raw_mat),
              all(RAP_raw_mat %in% c(1, 2, 3)))
    RAP_all_List <- list()
    for(i in 1:nrow(RAP_raw_mat)){
      RAP_all_List[[i]] <- integer(0)
      for(j in 1:ncol(RAP_raw_mat)){
        RAP_all_List[[i]] <- c(RAP_all_List[[i]], RAP_list[[j]][RAP_raw_mat[i, j]])
      }
      RAP_all_List[[i]] <- unlist(RAP_all_List[[i]])
    }
    new_TF_exp <- .TF_exp_mat
    
    for(j in 1:length(RAP_all_List[[RAP_nu]])){
      if(RAP_all_List[[RAP_nu]][j] < 0){
        new_TF_exp[-RAP_all_List[[RAP_nu]][j], ] <- rep(c(1,0), (ncol(.TF_exp_mat) / 2))
      }else if(RAP_all_List[[RAP_nu]][j] <= nrow(.TF_exp_mat)){
        new_TF_exp[ RAP_all_List[[RAP_nu]][j], ] <- rep(c(0,1), (ncol(.TF_exp_mat) / 2))
      }else{
        new_TF_exp[RAP_all_List[[RAP_nu]][j] - RAP_list_dif, ] <- rep(1, ncol(.TF_exp_mat))
      }
    }
  }else{
    new_TF_exp <- .TF_exp_mat
  }
  
  GeneFoldChange_by_TFFoldChange(Gene_exp_mat = .Gene_exp_mat,
                                 TF_exp_mat = new_TF_exp,
                                 TF_gene_site_mat = TF_gene_site_mat_list[[model_nu]],
                                 site_nu_thresh=.site_nu_thresh,
                                 TF_index = .TF_index)
}
########################################################################################################################
########################################################################################################################
# example
aa_exp_names_2 <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 1))) - 1399
aa_model_names <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 2)))
aa_model_names_org <- aa_model_index_1398_all[aa_model_names]

colnames(aa_modif_all_comb_p20_raw_1398)[2] <- "NR3C1"

GeneFoldChange_by_TFFoldChange_model_RAP_wrapper(.Gene_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_48_gene_gte3WTGRcommon,
                                                 .TF_exp_mat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_ER01_monomer_hocomoco,
                                                 TF_gene_site_mat_list = GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat, 
                                                 .site_nu_thresh = 1,
                                                 # .TF_index = c(1:nrow(.TF_exp_mat)),
                                                 RAP_raw_mat=aa_modif_all_comb_p20_raw_1398, 
                                                 RAP_list = aa_modif_all_options_p20_1398,
                                                 RAP_nu=aa_exp_names_2[1],
                                                 model_nu = aa_model_names_org[1])

########################################################################################################################
########################################################################################################################
KD_by_affected_enhancers <- function(WT_Pred_mat_list, 
                                     KD_Pred_mat_list,
                                     real_mat,
                                     TF_gene_site_mat_list,
                                     site_nu_thresh=1,
                                     model_index, 
                                     TF_index,
                                     plot_sc = F){
  # WT_Pred_mat_list : is a list where each entry corresponds to a model and is a matrix with same dims as real_mat, contains the predictions of the WT model
  # KD_Pred_mat_list : is a list where each entry corresponds to a model and is a matrix with same dims as real_mat, contains the predictions of the KD model
  # real_mat : is the gene expression matrix: each row a gene, each column a conditions. entries are 0, 1, -1 or NA
  # TF_gene_site_mat_list : is a list where each entry corresponds to a model. each entry is a matrix which contains the
  #  number of sites of each TF in each gene's enhancer. each
  #  Row indicates a genes and each column inidicates a TF. Rows and columns have to be named
  # site_nu_thresh : minimum number of sites of a TF within an enhancer, in order for the enhancer to be
  #  considred a target for that TF
  # model_index: is and integer vector containing the index of the models in WT_Pred_mat_list and KD_Pred_mat_list in the
  #  original models (what TF_gene_site_mat_list is sorted by)
  # TF_index is an integer indicating the TF that's knockdown. (index corresponds to TF_gene_site_mat_list entry column order)
  # plot_sc : if T it plots the results
  stopifnot(all(dim(real_mat) == dim(WT_Pred_mat_list[[1]])), 
            all(dim(KD_Pred_mat_list[[1]]) == dim(WT_Pred_mat_list[[1]])),
            length(WT_Pred_mat_list) == length(KD_Pred_mat_list),
            length(WT_Pred_mat_list) == length(model_index), 
            all(model_index <= length(TF_gene_site_mat_list)), 
            length(TF_index) == 1, 
            TF_index > 0,
            TF_index <= ncol(TF_gene_site_mat_list[[1]]))
  
  nu_models <- length(WT_Pred_mat_list)
  nu_genes <- nrow(real_mat)
  Acc_by_model_WT <- numeric(nu_models)
  Acc_by_model_KD <- numeric(nu_models)
  model_agreement <- numeric(nu_models)
  num_gene_per_model <- integer(nu_models)
  for(cur_model in 1:nu_models){
    cur_target_genes <- which(TF_gene_site_mat_list[[cur_model]][, TF_index] >= site_nu_thresh)
    num_gene_per_model[cur_model] <- length(cur_target_genes)
    Acc_by_model_WT[cur_model] <- sum(WT_Pred_mat_list[[cur_model]][cur_target_genes, ] == real_mat[cur_target_genes, ] , na.rm = T)/sum(!is.na(real_mat[cur_target_genes, ]))
    Acc_by_model_KD[cur_model] <- sum(KD_Pred_mat_list[[cur_model]][cur_target_genes, ] == real_mat[cur_target_genes, ] , na.rm = T)/sum(!is.na(real_mat[cur_target_genes, ]))
    model_agreement[cur_model] <- sum(KD_Pred_mat_list[[cur_model]][cur_target_genes, ] == WT_Pred_mat_list[[cur_model]][cur_target_genes, ] , na.rm = T)/sum(!is.na(real_mat[cur_target_genes, ]))
  }
  Acc_by_Model_diff <- Acc_by_model_WT - Acc_by_model_KD
  if(plot_sc){
    plot(num_gene_per_model,
         Acc_by_Model_diff,
         xlab = "Nu_targeted_genes", 
         ylab = "KD effect")
  }
  return(list(Nu_target_genes=num_gene_per_model, 
              KD_effect=Acc_by_Model_diff, 
              WT_KD_agreement = model_agreement))
}


########################################################################################################################
########################################################################################################################
# example
aa_WT_Pred_mat_list <- list()
aa_KD_Pred_mat_list <- list()
aa_exp_names <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 1)))
aa_model_names <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 2)))
aa_model_names_org <- aa_model_index_1398_all[aa_model_names]
aa_TF_index <- 3
for(i in 1:length(aa_exp_names)){
  aa_WT_Pred_mat_list[[i]] <- GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[aa_exp_names[i]]]$Prediction_round_list[[aa_model_names[i]]]
  aa_KD_Pred_mat_list[[i]] <- GEMSTAT_based_on_linear_exp_KD_Results[[aa_exp_names[i]]][[as.character(aa_TF_index)]]$Prediction_round_list[[aa_model_names[i]]]
}
aa <- KD_by_affected_enhancers(WT_Pred_mat_list = aa_WT_Pred_mat_list, 
                               KD_Pred_mat_list = aa_KD_Pred_mat_list,
                               real_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_48_gene_gte3WTGRcommon,
                               TF_gene_site_mat_list = GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat,
                               site_nu_thresh=1,
                               model_index = aa_model_names_org, 
                               TF_index=aa_TF_index,plot_sc = T)
########################################################################################################################
########################################################################################################################
########################################################################################################################
# write a function that goes through KD of all TFs in different RAPs and averages the 
#  effect of a TF KD on the gene's predicted expression over all models and RAPs, creates a heatmap where each row is a gene model combination:
# eg rows: gene_1_model_1, gene_1_model_2, gene_1_model_3, ...
# eg cols : TF_1, TF_2
KD_on_gene_evaluator <- function(WT_Prediction_list,
                                 KD_Prediction_list, 
                                 Avg_over_models, 
                                 .ct_treat_map = integer(0)){
  # WT_Prediction_list is a named list where each entry is a matrix of rounded predictions from the WT model, with 
  #  rows representing genes and column representing conditions
  # KD_Prediction_list is a named list of entries that each contain a Matrix with the same size as WT_Prediction_list, which
  #  are the predictions of a KD model. It's important that all entries are from the same KD experiment since
  #  the output might be averaged over all models
  # Avg_over_models : is a logical. if True it averages the change to gene's expression over multiple outputs of
  # .ct_treat_map : if provided it averages over over treatment and control conditions separately. this is an integer vector containing zero for contorl and 1 for treatmnet conditions.

  if(Avg_over_models){
    print(" It's important that all entries of KD_Prediction are from the same KD experiment since the output is averaged over all models")
  }
  
  gene_cond <- dim(WT_Prediction_list[[1]])
  size_mat_wt <- do.call(what = rbind,
                         args = lapply(WT_Prediction_list, dim))
  size_mat_kd <- do.call(what = rbind,
                         args = lapply(KD_Prediction_list, dim))
  # print("gene_cond")
  # print(gene_cond)
  # print("size_mat_wt")
  # print(size_mat_wt)
  # print("size_mat_kd")
  # print(size_mat_kd)
  
  nu_models <- length(KD_Prediction_list)
  gene_names <- rownames(WT_Prediction_list[[1]])
  ### DEBUG  ###   ###   ###   ### 
  stopifnot(nrow(size_mat_kd) == nu_models,
            nrow(size_mat_wt) == nu_models)
  ###   ###   ###   ###   ###   ### 
  stopifnot((all(size_mat_wt[, 1] %in%  gene_cond[1]) & all(size_mat_wt[, 2] %in%  gene_cond[2])),
            (all(size_mat_kd[, 1] %in%  gene_cond[1]) & all(size_mat_kd[, 2] %in%  gene_cond[2])),
            length(KD_Prediction_list) == length(WT_Prediction_list),
            is.logical(Avg_over_models),
            all(unlist(lapply(WT_Prediction_list, is.matrix))), 
            all(unlist(lapply(KD_Prediction_list, is.matrix))), 
            length(names(WT_Prediction_list)) == nu_models,
            length(names(KD_Prediction_list)) == nu_models, 
            all(names(KD_Prediction_list) == names(WT_Prediction_list)))
  if(length(.ct_treat_map) > 0){
    Diff_List_couple <- list()
    Diff_List_couple[[1]] <- list()
    Diff_List_couple[[2]] <- list()
    names(Diff_List_couple) <- c("control", "treatment")
    for(cur_model in 1:nu_models){
      Diff_List_couple[[1]][[cur_model]] <- KD_Prediction_list[[cur_model]][, .ct_treat_map %in% 0] - WT_Prediction_list[[cur_model]][, .ct_treat_map %in% 0]
      Diff_List_couple[[2]][[cur_model]] <- KD_Prediction_list[[cur_model]][, .ct_treat_map %in% 1] - WT_Prediction_list[[cur_model]][, .ct_treat_map %in% 1]
    }
    Diff_List_couple[[1]] <- lapply(Diff_List_couple[[1]], rowMeans, na.rm = T)
    Diff_List_couple[[2]] <- lapply(Diff_List_couple[[2]], rowMeans, na.rm = T)
    if(Avg_over_models){
      final_vector_control <- rowMeans(do.call(cbind, Diff_List_couple[[1]]))
      final_vector_treat <- rowMeans(do.call(cbind, Diff_List_couple[[2]]))
      gene_index <- c(1:length(final_vector_control))
      names(final_vector_control) <- gene_names[gene_index]
      names(final_vector_treat) <- gene_names[gene_index]
      stopifnot(length(final_vector_control) == gene_cond[1])
      return(list(Dif_vector = c(final_vector_control, final_vector_treat),
                  Gene_index = c(rep(1, length(final_vector_control)), rep(2, length(final_vector_treat))), # using this to distinguish control and treatment instead of distinguishing genes
                  Averaged = Avg_over_models))
    }else{
      print("not supported to have Avg_over_models== F and provide a .ct_treat_map")
      return(0)
    }
  }else{
    Diff_List <- list()
    for(cur_model in 1:nu_models){
      Diff_List[[cur_model]] <- KD_Prediction_list[[cur_model]] - WT_Prediction_list[[cur_model]]
    }
    
    Diff_List <- lapply(Diff_List, rowMeans, na.rm = T)
    if(Avg_over_models){
      final_vector <- rowMeans(do.call(cbind, Diff_List))
      gene_index <- c(1:length(final_vector))
      names(final_vector) <- gene_names[gene_index]
      stopifnot(length(final_vector) == gene_cond[1])
    }else{
      # Format the final vector such that entries for the same gene on different models are next to eachother
      final_vector <- as.numeric(do.call(rbind, Diff_List))
      gene_index <- as.numeric(matrix(data = rep(c(1:gene_cond[1]), nu_models),
                                      byrow = T,
                                      nrow = nu_models))
      model_index <- rep(c(1:nu_models), gene_cond[1])
      names(final_vector) <- paste(gene_index, model_index, sep = "_")
    }
    
    return(list(Dif_vector = final_vector,
                Gene_index = gene_index,
                Averaged = Avg_over_models))
  } # if not provided a ct map

}
########################################################################################################################
########################################################################################################################
# Example
aa_pr <- aa_WT2[[2131]]$Prediction_round_list
aa_prkd <- aa_KD2[[2131]]$`4`$Prediction_round_list
aa <- KD_on_gene_evaluator(WT_Prediction_list = aa_pr ,
                           KD_Prediction_list = aa_prkd,
                           Avg_over_models=F)

########################################################################################################################
########################################################################################################################
KD_on_gene_evaluator_multiTF_wrapper <- function(.WT_Prediction_list,
                                                 KD_Prediction_list_list,
                                                 .Avg_over_models=F,
                                                 KD_names,
                                                 plt=T,
                                                 plot_gene_index=integer(0), 
                                                 .breaks = numeric(0),
                                                 .colors = numeric(0),
                                                 ct_treat_map = integer(0)){
  # .WT_Prediction_list is a list where each entry is a matrix of rounded predictions from the WT model, with 
  #  rows representing genes and column representing conditions
  # KD_Prediction_list_list is a list where each entry is a list of entries that each contain a Matrix with the same size as WT_Prediction_list, which
  #  are the predictions of a KD model. It's important that all entries are from the same KD experiment since
  #  the output might be averaged over all models
  # .Avg_over_models : is a logical. if True it averages the change to gene's expression over multiple outputs of
  # plt : is a logical if True it plots a heatmap of the output
  # plot_gene_index : is the index of the genes to be plotted
  # ct_treat_map : if provided it averages over over treatment and control conditions separately. this is an integer vector containing zero for contorl and 1 for treatmnet conditions.
  # This function is a wrapper for KD_on_gene_evaluator. it outputs a matrix where each row corresponds to either
  #  a gene or a gene_model tuple. each column corresponds to a knockdown in KD_Prediction_list_list. entries are
  #  average change in predicted expression after that KD
  
  nu_KDs <- length(KD_Prediction_list_list)
  stopifnot(length(KD_names) == nu_KDs)
  Res_hold <- list()
  for(cur_kd in 1:nu_KDs){
    Res_hold[[cur_kd]] <- KD_on_gene_evaluator(WT_Prediction_list = .WT_Prediction_list ,
                                               KD_Prediction_list = KD_Prediction_list_list[[cur_kd]],
                                               Avg_over_models=.Avg_over_models,
                                               .ct_treat_map = ct_treat_map)
  }
  diff_mat <- do.call(cbind, lapply(Res_hold, "[[", 1))
  g_index <- Res_hold[[1]]$Gene_index
  colnames(diff_mat) <- KD_names
  rownames(diff_mat) <- names(Res_hold[[1]]$Dif_vector)
  nu_genes <- length(unique(g_index))
  if(plt){
    all_colors <- rgb(runif(nu_genes), 
                      runif(nu_genes), 
                      runif(nu_genes))
    my_rowside_color <- all_colors[g_index]
    if(ncol(diff_mat) == 1){
      diff_mat <- cbind(diff_mat,
                        rep(0, nrow(diff_mat)))
      colnames(diff_mat)[2] <- "dummy"
    }
    if(length(plot_gene_index) > 0){
      diff_mat_plot <- diff_mat[g_index %in% plot_gene_index, ]
      my_rowside_color <- my_rowside_color[g_index %in% plot_gene_index]
    }else{
      diff_mat_plot <- diff_mat
    }
    if(length(.breaks) == 0){
      my_range <- c(-abs(max(range(diff_mat_plot))),
                    abs(max(range(diff_mat_plot))))
      my_breaks <- seq(my_range[1],
                       my_range[2],
                       length.out = 21)
    }else{
      my_breaks <- .breaks
    }
    if(length(.colors)==0){
      #my_cols <- rainbow(n=20, alpha=1)
      my_cols <- colorspace::diverge_hsv(20)
    }else{
      my_cols <- .colors
    }
    heatmap.2(diff_mat_plot, 
              Rowv = T,
              Colv = F, 
              dendrogram = "row",
              trace = "none",
              breaks = my_breaks,
              # col = rgb(r = seq(0, 1, length.out = 100),
              #           g = seq(0.1, 0, length.out = 100), 
              #           blue = seq(1, 0, length.out = 100)),
              RowSideColors = my_rowside_color,
              col = my_cols,
              margins = c(8,8), main = "KD effect per gene")
  }
  return(diff_mat)
}
########################################################################################################################
########################################################################################################################
# Example
library(gplots)
aa_pr <- aa_WT2[[2139]]$Prediction_round_list
aa_prkd <- lapply(aa_KD2[[2139]], "[[", 2)
aa <- KD_on_gene_evaluator_multiTF_wrapper(.WT_Prediction_list = aa_pr ,
                                           KD_Prediction_list_list = aa_prkd,
                                           .Avg_over_models=T, 
                                           KD_names = c("ER_KD"))
########################################################################################################################
########################################################################################################################
# write a function that takes the annotation per model per gene, takes the index of one TF and outputs:
# the number of overlaps of the sites of that TF with other TFs per model per gene
# overall number of overlaps with each TF per model
annotation_overlap_finder <- function(annoatation_witer_input, TF_names, TF_index){
  # annoatation_witer_input : is the input to annotation_Writer function. is a list that
  #  contains one entry per model and that entry contains one entry per gene which contains the
  #  annotations 
  # TF_names : is a character vector containing the names of all TFs used in the model
  # TF_index : is an integer containing the index of the TF to be examined
  
  stopifnot(TF_index <= length(TF_names), 
            length(TF_index) == 1)
  nu_models <- length(annoatation_witer_input)
  nu_genes <- length(annoatation_witer_input[[1]])
  des_tf <- TF_names[TF_index]
  res_holder <- list()
  for(cur_model in 1:nu_models){
    res_holder[[cur_model]] <- list()
    for(cur_gene in 1:nu_genes){
      res_holder[[cur_model]][[cur_gene]] <- character(0)
      cur_annot <- annoatation_witer_input[[cur_model]][[cur_gene]]
      if(! des_tf %in% cur_annot[, 3]){
        next()
      }
      # get the positions of the current hits of the desired TF
      des_pos <- cur_annot[cur_annot[, 3] %in% des_tf ,1]
      
      des_pos_split <- strsplit(des_pos, split = "\\..")
      des_pos1 <- as.integer(unlist(lapply(des_pos_split, "[[", 1)))
      des_pos2 <- as.integer(unlist(lapply(des_pos_split, "[[", 2)))
      des_pos_mat <- cbind(des_pos1, 
                           des_pos2)
      nondes_pos_split <-  strsplit(cur_annot[! (cur_annot[, 3] %in% des_tf) ,1],
                                    split = "\\..")
      nondes_pos_split_1 <- as.integer(unlist(lapply(nondes_pos_split, "[[", 1)))
      nondes_pos_split_2 <- as.integer(unlist(lapply(nondes_pos_split, "[[", 2)))
      nondes_pos_mat <- cbind(nondes_pos_split_1,
                              nondes_pos_split_2)
      nondes_pos_mat_names <- cur_annot[! (cur_annot[, 3] %in% des_tf) ,3]
      which_overlap <- nondes_pos_mat_names[range_overlap(r_a = des_pos_mat, 
                                                          r_b = nondes_pos_mat)]
      res_holder[[cur_model]][[cur_gene]] <- c(res_holder[[cur_model]][[cur_gene]],
                                               which_overlap)
    } #end of loop over genes
  } # end of loop over models
  return(res_holder)
}
########################################################################################################################
########################################################################################################################
# Example
aa <- annotation_overlap_finder(annoatation_witer_input = GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Annot_per_model_perGene,
                                TF_names = names(TF.motifs.Shrinked.hocomoco_E2F1_added.count),
                                TF_index = 4)
########################################################################################################################
########################################################################################################################
range_overlap <- function(r_a, r_b){
  # r_a is a matrix with 2 columns and arbitrary number of rows, each row indicates a
  #  range with first entry marking the beginning of the range and the second entry marking the end. 
  # r_b has the same format as r_a
  # it returns the index of the rows in r_b which overlap with rows in r_a
  
  # checking the format of inputs
  stopifnot(length(r_a) >= 2,
            length(r_b) >= 2)
  if(length(r_a) == 2){
    r_a <- matrix(r_a, byrow = T, ncol = 2)
  }
  if(length(r_b) == 2){
    r_b <- matrix(r_b, byrow = T, ncol = 2)
  }
  stopifnot(is.matrix(r_a), 
            is.matrix(r_b), 
            all(r_a[, 1] < r_a[, 2]),
            all(r_b[, 1] < r_b[, 2]))
  overlap_res <- logical(nrow(r_b))
  for(i in 1:nrow(r_b)){
    for(j in 1:nrow(r_a)){
      cur_ol <- (r_a[j, 1] <= r_b[i, 2] & r_a[j, 2] >= r_b[i, 1])
      if(cur_ol){
        overlap_res[i] <- T
        break()
      }
    }
  }
  return(which(overlap_res))
}
########################################################################################################################
########################################################################################################################
# example

range_overlap(r_a = cbind(c(1, 3), c(4, 6)), 
              r_b = cbind(c(1,4,5,7,8), c(2,5,9,8,10)))
########################################################################################################################
########################################################################################################################
# next function starts here





####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

# workplace to be loaded in hal

TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10 <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10[10, seq(1, 90, 2)] <- rep(1, ncol(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10)/2)
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10[10, seq(2, 90, 2)] <- rep(0, ncol(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10)/2)

plotExpression(expMat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10,
               exportplot = T,
               filename = "TF_Exp_GRn.png",
               .dendrogram = "none",
               .Rowv = F, .Colv = F, colorPoints = c(-1, 0, 0.33, 0.66, 1))


TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01 <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01[14, seq(1, 90, 2)] <- rep(0, ncol(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01)/2)
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01[14, seq(2, 90, 2)] <- rep(1, ncol(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01)/2)

plotExpression(expMat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01,
               exportplot = T,
               filename = "TF_Exp_RARAp.png",
               .dendrogram = "none",
               .Rowv = F, .Colv = F, colorPoints = c(-1, 0, 0.33, 0.66, 1))

TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10_RARA01 <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10_RARA01[14, seq(1, 90, 2)] <- rep(0, ncol(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01)/2)
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10_RARA01[14, seq(2, 90, 2)] <- rep(1, ncol(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01)/2)

plotExpression(expMat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10_RARA01,
               exportplot = T,
               filename = "TF_Exp_RARAp_GRn.png",
               .dendrogram = "none",
               .Rowv = F, .Colv = F, colorPoints = c(-1, 0, 0.33, 0.66, 1))

Sim_Ann_weighted_148_restart_parameters_GR10 <- do.call(rbind, lapply(Optim_Greedy_conc_modif_Ensemble_Results[[4]], "[[", 3))
colnames(Sim_Ann_weighted_148_restart_parameters_GR10) <- colnames(Sim_Ann_weighted_148_restart_parameters)

Sim_Ann_weighted_148_restart_parameters_RAR10 <- do.call(rbind, lapply(Optim_Greedy_conc_modif_Ensemble_Results[[7]], "[[", 3))
colnames(Sim_Ann_weighted_148_restart_parameters_RAR10) <- colnames(Sim_Ann_weighted_148_restart_parameters)

Sim_Ann_weighted_148_restart_parameters_GR10_RAR10 <- do.call(rbind, lapply(Optim_Greedy_conc_modif_DoubleKD_Results, "[[", 3))
colnames(Sim_Ann_weighted_148_restart_parameters_GR10_RAR10) <- colnames(Sim_Ann_weighted_148_restart_parameters)

hal_workplace_exp1 <- ls()
#rm(list=setdiff(ls(), c(hal_workplace_exp1, "MotifWriter")))
########################################################################################################################
########################################################################################################################
####################################################################################################
# Reading exp 1 - 4 results
GEMSTAT_based_on_linear_exp_results <- list()
GEMSTAT_based_on_linear_exp_pars <- list()
GEMSTAT_based_on_linear_exp1 <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = "Seeded_GEMSTAT_ens/Experiment_1/Outputs", .fold_change = T)
GEMSTAT_based_on_linear_exp2 <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = "Seeded_GEMSTAT_ens/Experiment_2/Outputs", .fold_change = T)
GEMSTAT_based_on_linear_exp3 <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = "Seeded_GEMSTAT_ens/Experiment_3/Outputs", .fold_change = T)
GEMSTAT_based_on_linear_exp4 <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = "Seeded_GEMSTAT_ens/Experiment_4/Outputs", .fold_change = T)
GEMSTAT_based_on_linear_exp_results[[1]] <- GEMSTAT_based_on_linear_exp1
GEMSTAT_based_on_linear_exp_results[[2]] <- GEMSTAT_based_on_linear_exp2
GEMSTAT_based_on_linear_exp_results[[3]] <- GEMSTAT_based_on_linear_exp3
GEMSTAT_based_on_linear_exp_results[[4]] <- GEMSTAT_based_on_linear_exp4

for(i in 1:4){
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
}
names(GEMSTAT_based_on_linear_exp_pars) <- paste0("experiment_", c(1:4))
names(GEMSTAT_based_on_linear_exp_results) <- paste0("experiment_", c(1:4))

par(mfrow = c(2, 2))

aa_ind <- sort(apply(GEMSTAT_based_on_linear_exp1$Accuracy, MARGIN = 2, median, na.rm=T), decreasing = T, index.return=T)$ix
boxplot.matrix(GEMSTAT_based_on_linear_exp1$Accuracy[, aa_ind], ylim =c(0, 1), main= "WT")

aa_ind <- sort(apply(GEMSTAT_based_on_linear_exp2$Accuracy, MARGIN = 2, median, na.rm=T), decreasing = T, index.return=T)$ix
boxplot.matrix(GEMSTAT_based_on_linear_exp2$Accuracy[, aa_ind], ylim =c(0, 1), main= "GR")

aa_ind <- sort(apply(GEMSTAT_based_on_linear_exp3$Accuracy, MARGIN = 2, median, na.rm=T), decreasing = T, index.return=T)$ix
boxplot.matrix(GEMSTAT_based_on_linear_exp3$Accuracy[, aa_ind], ylim =c(0, 1), main= "RARA")

aa_ind <- sort(apply(GEMSTAT_based_on_linear_exp4$Accuracy, MARGIN = 2, median, na.rm=T), decreasing = T, index.return=T)$ix
boxplot.matrix(GEMSTAT_based_on_linear_exp4$Accuracy[, aa_ind], ylim =c(0, 1), main= "GR_RARA")

par(mfrow = c(2, 2))

aa_ind <- sort(apply(GEMSTAT_based_on_linear_exp1$RMSE, MARGIN = 2, median, na.rm=T), decreasing = F, index.return=T)$ix
boxplot.matrix(GEMSTAT_based_on_linear_exp1$RMSE[, aa_ind], ylim =c(0, 2), main= "WT")

aa_ind <- sort(apply(GEMSTAT_based_on_linear_exp2$RMSE, MARGIN = 2, median, na.rm=T), decreasing = F, index.return=T)$ix
boxplot.matrix(GEMSTAT_based_on_linear_exp2$RMSE[, aa_ind], ylim =c(0, 2), main= "GR")

aa_ind <- sort(apply(GEMSTAT_based_on_linear_exp3$RMSE, MARGIN = 2, median, na.rm=T), decreasing = F, index.return=T)$ix
boxplot.matrix(GEMSTAT_based_on_linear_exp3$RMSE[, aa_ind], ylim =c(0, 2), main= "RARA")

aa_ind <- sort(apply(GEMSTAT_based_on_linear_exp4$RMSE, MARGIN = 2, median, na.rm=T), decreasing = F, index.return=T)$ix
boxplot.matrix(GEMSTAT_based_on_linear_exp4$RMSE[, aa_ind], ylim =c(0, 2), main= "GR_RARA")


par(mfrow = c(1, 1))
plot(c(GEMSTAT_based_on_linear_exp1$Accuracy_All, 
       GEMSTAT_based_on_linear_exp2$Accuracy_All,
       GEMSTAT_based_on_linear_exp3$Accuracy_All, 
       GEMSTAT_based_on_linear_exp4$Accuracy_All))

plot(c(sort(unlist(Sim_Ann_weighted_148_restart_round_accuracy)),
       sort(Sim_Ann_weighted_148_restart_random_pred$Precision), c(sort(GEMSTAT_based_on_linear_exp1$Accuracy_All), 
                                                                   sort(GEMSTAT_based_on_linear_exp2$Accuracy_All),
                                                                   sort(GEMSTAT_based_on_linear_exp3$Accuracy_All), 
                                                                   sort(GEMSTAT_based_on_linear_exp4$Accuracy_All))),
     ylim = c(0, 1),
     col = c(rep(3, 148), rep(2, 148), rep(4, 148), rep(5, 148), rep(6, 148), rep(7, 148)),
     ylab = "Accuracy",
     xlab = "model index",
     cex = 0.8,
     pch = 16)
legend(500, 1, legend=c("Sim_Ann_weighted", "shuffled", "WT_GEM", "GR_GEM", "RARA_GEM", "GR_RARA_GEM"),
       col=c(3, 2, 4, 5, 6, 7), pch=16, cex=1.1,
       bty = "n", x.intersp = 0.1,
       y.intersp = 0.2, pt.cex = 1.1, box.col = 1)

####################################################################################################
# reading results experiments 5-12

for(i in 5:8){
  
  GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                       role_dec = T,
                                                                       TFinfo_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Inputs/TF_info"))
  
}
names(GEMSTAT_based_on_linear_exp_results)[5:8] <- paste0("experiment_", c(5:8))
names(GEMSTAT_based_on_linear_exp_pars)[5:8] <- paste0("experiment_", c(5:8))

plot(unlist(lapply(GEMSTAT_based_on_linear_exp_results, "[[", 7)))
plotExpression(expMat = GEMSTAT_based_on_linear_exp_results$experiment_5$Prediction_round_list$`2`, .dendrogram = "none", .Rowv = F, .Colv = F, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),filename = "pred.png")
plot(c(sort(unlist(Sim_Ann_weighted_148_restart_round_accuracy)),
       sort(Sim_Ann_weighted_148_restart_random_pred$Precision), c(sort(GEMSTAT_based_on_linear_exp1$Accuracy_All), 
                                                                   sort(GEMSTAT_based_on_linear_exp2$Accuracy_All),
                                                                   sort(GEMSTAT_based_on_linear_exp3$Accuracy_All), 
                                                                   sort(GEMSTAT_based_on_linear_exp4$Accuracy_All)) ,
       unlist(lapply(GEMSTAT_based_on_linear_exp_results[5:8], "[[", 7)) ),
     ylim = c(0, 1),
     col = c(rep(3, 148), rep(2, 148), rep(4, 148), rep(5, 148), rep(6, 148), rep(7, 148), rep(8, 148), rep(9, 148), rep(10, 148), rep(11, 148)),
     ylab = "Accuracy",
     xlab = "model index",
     cex = 0.8,
     pch = 16)


for(i in 9:12){
  GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
}
names(GEMSTAT_based_on_linear_exp_results)[9:12] <- paste0("experiment_", c(9:12))
names(GEMSTAT_based_on_linear_exp_pars)[9:12] <- paste0("experiment_", c(9:12))

plot(unlist(lapply(GEMSTAT_based_on_linear_exp_results, "[[", 7)))
####################################################################################################
####################################################################################################
#gathering all accuracy results in one matrix
Linear_and_GEMSTAT_Accuracy_all <- matrix(nrow = 148, ncol = 1)
Linear_and_GEMSTAT_Accuracy_all[,1] <- Sim_Ann_weighted_148_restart_random_pred$Precision
colnames(Linear_and_GEMSTAT_Accuracy_all)[1] <- "Random"
Linear_and_GEMSTAT_Accuracy_all <- cbind(Linear_and_GEMSTAT_Accuracy_all,
                                         unlist(Sim_Ann_weighted_148_restart_round_accuracy),
                                         Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, 4],
                                         Optim_Greedy_conc_modif_Ensemble_Results_RoundAccuracy[, 7],
                                         Optim_Greedy_conc_modif_Double_Results_RoundAccuracy)
colnames(Linear_and_GEMSTAT_Accuracy_all)[2:5] <- c("linear_wt", "linear_GR", "linear_RAR", "linear_double")

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results[1:28], "[[", 7))
Linear_and_GEMSTAT_Accuracy_all <- cbind(Linear_and_GEMSTAT_Accuracy_all, aa)
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all, las=2)
####################################################################################################
# reading results experiment 13
# this is the same as exp 5 only has one beta for all genes.


GEMSTAT_based_on_linear_exp_results[[13]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 13, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[13]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 13, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[13] <- paste0("experiment_", 13)

plot(unlist(lapply(GEMSTAT_based_on_linear_exp_results, "[[", 7)))

####################################################################################################
# reading results experiment 14
# this is the same as exp 13 but only has one qBTM for all genes.

GEMSTAT_based_on_linear_exp_results[[14]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 14, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[14]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 14, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[14] <- paste0("experiment_", 14)

plot(unlist(lapply(GEMSTAT_based_on_linear_exp_results, "[[", 7)))

####################################################################################################
# reading results experiment 15 (GR-RAR modified no coop)
# this is the same as exp 14 but only uses genes that were predicted better than random in the linear model.
GEMSTAT_based_on_linear_exp_results[[15]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 15, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[15]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 15, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[15] <- paste0("experiment_", 15)

plot(unlist(lapply(GEMSTAT_based_on_linear_exp_results, "[[", 7)))

####################################################################################################
# reading results experiment 16 (GR-RAR modified)
# this is the same as exp 14 but only uses factorInfo (predetermined factor roles)
GEMSTAT_based_on_linear_exp_results[[16]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_",
                                                                                                         16, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[16]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 16, "/Outputs"),
                                                                      role_dec = T, TFinfo_dir = "Seeded_GEMSTAT_ens/Experiment_16/Inputs/TF_info")
names(GEMSTAT_based_on_linear_exp_results)[16] <- paste0("experiment_", 16)

plot(unlist(lapply(GEMSTAT_based_on_linear_exp_results, "[[", 7)))

####################################################################################################
# reading results experiment 17 (GR-RAR modified)
# this is the same as exp 15 but only uses factorInfo (predetermined factor roles)
GEMSTAT_based_on_linear_exp_results[[17]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_",
                                                                                                         17, "/Outputs"),
                                                                                     .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[17]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 17, "/Outputs"),
                                                                      role_dec = T, TFinfo_dir = "Seeded_GEMSTAT_ens/Experiment_17/Inputs/TF_info")
names(GEMSTAT_based_on_linear_exp_results)[17] <- paste0("experiment_", 17)

plot(unlist(lapply(GEMSTAT_based_on_linear_exp_results[c(15,17)], "[[", 7)))
plot(sort(unlist(lapply(GEMSTAT_based_on_linear_exp_results[15], "[[", 7))))
points(sort(unlist(lapply(GEMSTAT_based_on_linear_exp_results[17], "[[", 7))), col=2)



####################################################################################################
# reading results experiment 18  (non modified ER-GR)
# this is the same as exp 15 but only uses coop between ER and GR

GEMSTAT_based_on_linear_exp_results[[18]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 18, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[18]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 18, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[18] <- paste0("experiment_", 18)
names(GEMSTAT_based_on_linear_exp_pars)[18] <- paste0("experiment_", 18)


####################################################################################################
# reading results experiment 19  (GR modified ER-GR)
# this is the same as exp 15 but uses the inputs of GR modified models and coop between ER and GR
GEMSTAT_based_on_linear_exp_results[[19]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 19, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[19]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 19, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[19] <- paste0("experiment_", 19)
names(GEMSTAT_based_on_linear_exp_pars)[19] <- paste0("experiment_", 19)

####################################################################################################
# reading results experiment 20 (GR modified no coop)
# this is the same as exp 19 but no coop
GEMSTAT_based_on_linear_exp_results[[20]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 20, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[20]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 20, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[20] <- paste0("experiment_", 20)
names(GEMSTAT_based_on_linear_exp_pars)[20] <- paste0("experiment_", 20)

plot(unlist(lapply(GEMSTAT_based_on_linear_exp_results, "[[", 7)))
plot(sort(unlist(lapply(GEMSTAT_based_on_linear_exp_results[15], "[[", 7))))
points(sort(unlist(lapply(GEMSTAT_based_on_linear_exp_results[20], "[[", 7))), col=2)


####################################################################################################
# reading results experiment 21 (non-modified ER-RARA)
# this is the same as exp 15 but uses coop between ER and RARA
GEMSTAT_based_on_linear_exp_results[[21]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 21, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[21]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 21, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[21] <- paste0("experiment_", 21)
names(GEMSTAT_based_on_linear_exp_pars)[21] <- paste0("experiment_", 21)


####################################################################################################
# reading results experiment 22 (RARA modified ER-RARA)
# this is the same as exp 15 but uses  the inputs of RARA modified models and coop between ER and RARA
GEMSTAT_based_on_linear_exp_results[[22]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 22, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[22]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 22, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[22] <- paste0("experiment_", 22)
names(GEMSTAT_based_on_linear_exp_pars)[22] <- paste0("experiment_", 22)


####################################################################################################
# reading results experiment 23 (RARA modified no coop)
# this is the same as exp 22 but uses no coop

GEMSTAT_based_on_linear_exp_results[[23]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 23, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[23]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 23, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[23] <- paste0("experiment_", 23)
names(GEMSTAT_based_on_linear_exp_pars)[23] <- paste0("experiment_", 23)

plot(unlist(lapply(GEMSTAT_based_on_linear_exp_results, "[[", 7)))
plot(sort(unlist(lapply(GEMSTAT_based_on_linear_exp_results[15], "[[", 7))))
points(sort(unlist(lapply(GEMSTAT_based_on_linear_exp_results[23], "[[", 7))), col=2)
boxplot.matrix(cbind(unlist(lapply(GEMSTAT_based_on_linear_exp_results[15], "[[", 7)),
                     unlist(lapply(GEMSTAT_based_on_linear_exp_results[23], "[[", 7))))
####################################################################################################
# reading results experiment 24 (RARA-GR modified ER-GR, ER-RARA)
# this is the same as exp 15 but uses  the inputs of RARA-GR modified models and coop between ER and RARA , ER and GR
GEMSTAT_based_on_linear_exp_results[[24]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 24,
                                                                                                         "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[24]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 24, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[24] <- paste0("experiment_", 24)
names(GEMSTAT_based_on_linear_exp_pars)[24] <- paste0("experiment_", 24)

####################################################################################################
# reading results experiment 25 (RARA-GR modified ER-RARA)
# this is the same as exp 24 but uses  the inputs of RARA-GR modified models and coop between ER and RARA

####################################################################################################
# reading results experiment 26 (RARA-GR modified ER-GR)
# this is the same as exp 24 RARA-GR modified but uses coop between  ER and GR

####################################################################################################
# reading results experiment 27 (no modif no coop)
# this is the same as exp 24 but uses  no modif and no coop

####################################################################################################
# reading results experiment 28 (no modif, ER-RARA, ER-GR)
# this is the same as exp 15 but uses coop between ER and RARA , ER and GR

for(i in 25:28){
  GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 
                                                                                                          i, "/Outputs"), .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
ncol(Linear_and_GEMSTAT_Accuracy_all)
setwd("Seeded_GEMSTAT_ens/Plot_best_per_exp/")
for(i in 1:length(GEMSTAT_based_on_linear_exp_results)){
  print(i)
  aa_max_ind <- which.max(GEMSTAT_based_on_linear_exp_results[[i]][[7]])
  plotExpression(expMat = GEMSTAT_based_on_linear_exp_results[[i]]$Prediction_round_list[[aa_max_ind]], 
                 .dendrogram = "none", 
                 .Rowv = F, 
                 .Colv = F,
                 colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),
                 filename = paste0("pred_", i, ".png"))
}
setwd("../..")


plotExpression(expMat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52, 
               .dendrogram = "none", 
               .Rowv = F, 
               .Colv = F,
               colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),
               filename = "real.png")


####################################################################################################
#create a dataframe containing the information about each experiment

GEMSTAT_experiment_info_df <- matrix(nrow = 28, ncol = 8)
colnames(GEMSTAT_experiment_info_df) <- c("name", "linear_Base", "qBTM", "Beta", "TFinfo", "Cooperativity", "nu_free_par", "comments")
GEMSTAT_experiment_info_df[, 1] <- paste("experiment", c(1:28), sep = "_")

GEMSTAT_experiment_info_df[, 2] <- c("WT", "GR", "RAR", "Double",
                                     "WT", "GR", "RAR","Double",
                                     "WT", "GR", "RAR", "Double",
                                     "Double","Double", "Double", "Double",
                                     "Double", "WT", "GR","GR",
                                     "WT", "RAR", "RAR", "Double",
                                     "Double","Double", "WT", "WT")

GEMSTAT_experiment_info_df[, 3] <- c("per_gene", "per_gene", "per_gene", "per_gene",
                                     "per_gene", "per_gene", "per_gene", "per_gene",
                                     "per_gene", "per_gene", "per_gene", "per_gene",
                                     "per_gene", "one", "one", "one",
                                     "one", "one", "one", "one",
                                     "one", "one", "one", "one",
                                     "one", "one", "one", "one")
GEMSTAT_experiment_info_df[, 4] <- c("one", "one","one","one",
                                     "per_gene", "per_gene", "per_gene", "per_gene",
                                     "per_gene", "per_gene", "per_gene", "per_gene",
                                     "one", "one", "one", "one",
                                     "one", "one", "one", "one", 
                                     "one", "one", "one", "one", 
                                     "one", "one", "one", "one")
GEMSTAT_experiment_info_df[, 5] <- c("no", "no", "no", "no",
                                     "yes", "yes", "yes", "yes",
                                     "no", "no", "no", "no",
                                     "no", "no", "no", "yes",
                                     "yes", "no", "no", "no", 
                                     "no", "no", "no", "no",
                                     "no", "no", "no", "no")

GEMSTAT_experiment_info_df[, 6] <- c("none", "none", "none", "none",
                                     "none", "none", "none", "none",
                                     "none", "none", "none", "none",
                                     "none", "none", "none", "none",
                                     "none", "ER-GR_SIMPLE_dist_50", "ER-GR_SIMPLE_dist_50", "none",
                                     "ER-RARA_SIMPLE_dist_50", "ER-RARA_SIMPLE_dist_50", "none", "ER-GR_SIMPLE_dist_50|ER-RARA_SIMPLE_dist_50",
                                     "ER-RARA_SIMPLE_dist_50", "ER-GR_SIMPLE_dist_50", "none", "ER-GR_SIMPLE_dist_50|ER-RARA_SIMPLE_dist_50")
aa_remove_bad_genes <- c("yes", "yes", "yes", "yes",
                         "no", "no", "no", "no", 
                         "no", "no", "no", "no",
                         "no", "no", "yes", "no", 
                         "yes", "yes", "yes", "yes",
                         "yes", "yes", "yes", "yes",
                         "yes", "yes", "yes", "yes") 
GEMSTAT_experiment_info_df <- cbind(GEMSTAT_experiment_info_df, aa_remove_bad_genes)
colnames(GEMSTAT_experiment_info_df)[9] <- "remove_linear_bad_genes"

####################################################################################################
####################################################################################################
Linear_and_GEMSTAT_Accuracy_all # contains all accuracies so far
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all, las = 2, ylab="accuracy")
colnames(Linear_and_GEMSTAT_Accuracy_all)[2:5] <- c("linear_wt", "linear_GR", "linear_RAR", "linear_double")
colnames(Linear_and_GEMSTAT_Accuracy_all)[6:33] <- paste("experiment", c(1:28), sep = "_")


####################################################################################################
####################################################################################################

####################################################################################################
# reading results experiment 29 (no modif, no coop, filtred set of genes(based on common enhancers))
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[28,])
GEMSTAT_experiment_info_df[29, 1] <- "experiment_29"
GEMSTAT_experiment_info_df[29, 8] <- "Filtered set of genes based on (at least 3) common enhancers among top10 linear models"

####################################################################################################
# reading results experiment 30 (no modif, ER-GR, filtred set of genes(based on common enhancers))
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[29,])
GEMSTAT_experiment_info_df[30, 1] <- "experiment_30"
GEMSTAT_experiment_info_df[30, 6]  <- "ER-GR_SIMPLE_dist_50"

####################################################################################################
# reading results experiment 31 (no modif, ER-RARA, filtred set of genes(based on common enhancers))
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[30,])
GEMSTAT_experiment_info_df[31, 1] <- "experiment_31"
GEMSTAT_experiment_info_df[31, 6]  <- "ER-RARA_SIMPLE_dist_50"
####################################################################################################
# reading results experiment 32 (no modif, ER-GR, ER-RARA filtred set of genes(based on common enhancers))
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[31,])
GEMSTAT_experiment_info_df[32, 1] <- "experiment_32"
GEMSTAT_experiment_info_df[32, 6]  <- "ER-GR_SIMPLE_dist_50|ER-RARA_SIMPLE_dist_50"

for(i in 29:32){
  GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results[29:32], "[[", 7))
Linear_and_GEMSTAT_Accuracy_all <- cbind(Linear_and_GEMSTAT_Accuracy_all, aa)
colnames(Linear_and_GEMSTAT_Accuracy_all)[34:37] <- paste("experiment", c(29:32), sep = "_")
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all, las=2)

####################################################################################################
####################################################################################################
# Creating half motifs for nuclear receptors
# monomer
library(seqLogo)

TF.motifs.Shrinked.halfsites.count <- list()
# AR
# seqLogo::seqLogo(TF.motifs.Shrinked.t[[1]])  
TF.motifs.Shrinked.halfsites.count[[1]] <- TF.motifs.Shrinked.count[[1]][1:6, ]
# CEBPB
TF.motifs.Shrinked.halfsites.count[[2]] <- TF.motifs.Shrinked.count[[2]]
# ESR1
# seqLogo::seqLogo(TF.motifs.Shrinked.t[[3]])  
TF.motifs.Shrinked.halfsites.count[[3]] <- TF.motifs.Shrinked.count[[3]][1:6, ]
# FOXA1
TF.motifs.Shrinked.halfsites.count[[4]] <- TF.motifs.Shrinked.count[[4]]
# GATA3
TF.motifs.Shrinked.halfsites.count[[5]] <- TF.motifs.Shrinked.count[[5]]
# JUN_1
TF.motifs.Shrinked.halfsites.count[[6]] <- TF.motifs.Shrinked.count[[6]]
# JUN_2
TF.motifs.Shrinked.halfsites.count[[7]] <- TF.motifs.Shrinked.count[[7]]
# LEF1
TF.motifs.Shrinked.halfsites.count[[8]] <- TF.motifs.Shrinked.count[[8]]
# NKX3-1
TF.motifs.Shrinked.halfsites.count[[9]] <- TF.motifs.Shrinked.count[[9]]
# NR3C1
# seqLogo::seqLogo(TF.motifs.Shrinked.t[[10]])  
TF.motifs.Shrinked.halfsites.count[[10]] <- TF.motifs.Shrinked.count[[10]][1:6, ]
# NR5A2
TF.motifs.Shrinked.halfsites.count[[11]] <- TF.motifs.Shrinked.count[[11]]
# PBX1
TF.motifs.Shrinked.halfsites.count[[12]] <- TF.motifs.Shrinked.count[[12]]
# PGR_s1
# seqLogo::seqLogo(TF.motifs.Shrinked.t[[13]])  
TF.motifs.Shrinked.halfsites.count[[13]] <- TF.motifs.Shrinked.count[[13]][1:4, ]
# PGR_s2
TF.motifs.Shrinked.halfsites.count[[14]] <- TF.motifs.Shrinked.count[[13]][8:13, ]
# RARA
TF.motifs.Shrinked.halfsites.count[[15]] <- TF.motifs.Shrinked.count[[14]][1:6, ]
# RARG
TF.motifs.Shrinked.halfsites.count[[16]] <- TF.motifs.Shrinked.count[[15]][1:6, ]
# RUNX1
TF.motifs.Shrinked.halfsites.count[[17]] <- TF.motifs.Shrinked.count[[16]]
# RXRA
TF.motifs.Shrinked.halfsites.count[[18]] <- TF.motifs.Shrinked.count[[17]][8:13,]
# SP1
TF.motifs.Shrinked.halfsites.count[[19]] <- TF.motifs.Shrinked.count[[18]]
# TFAP2C
TF.motifs.Shrinked.halfsites.count[[20]] <- TF.motifs.Shrinked.count[[19]]

names(TF.motifs.Shrinked.halfsites.count)[1:12] <- names(TF.motifs.Shrinked.count)[1:12]
names(TF.motifs.Shrinked.halfsites.count)[13] <- "PGR_s1"
names(TF.motifs.Shrinked.halfsites.count)[14] <- "PGR_s2"
names(TF.motifs.Shrinked.halfsites.count)[15:20] <- names(TF.motifs.Shrinked.count)[14:19]

# Modifying TF expression profile
# WT
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer <-rbind(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01[1:13, ],
                                                                                  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01[13, ],
                                                                                  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01[14:19, ])
rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer)[13:14] <- c("PGR_s1", "PGR_s2")
# GR
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10_monomer <-rbind(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10[1:13, ],
                                                                                  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10[13, ],
                                                                                  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10[14:19, ])
rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10_monomer)[13:14] <- c("PGR_s1", "PGR_s2")

# RARA
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01_monomer <-rbind(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01[1:13, ],
                                                                                  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01[13, ],
                                                                                  TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01[14:19, ])
rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_RARA01_monomer)[13:14] <- c("PGR_s1", "PGR_s2")

# Double
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10_RARA01_monomer <-rbind(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10_RARA01[1:13, ],
                                                                                       TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10_RARA01[13, ],
                                                                                       TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10_RARA01[14:19, ])
rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_GR10_RARA01_monomer)[13:14] <- c("PGR_s1", "PGR_s2")
####################################################################################################
####################################################################################################
# experiment 33
# same as experiment 32 but with no coop between ER-RAR and no coop between ER and GR. using monomer motifs (coop between pairs of monomers)
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[32,])
GEMSTAT_experiment_info_df[33, 1] <- "experiment_33"
GEMSTAT_experiment_info_df[33, 6]  <- "AR-AR_DIMER_dist_3|ER-ER_DIMER_dist_3|NR3C1-NR3C1_DIMER_dist_3|PGRs1-PGRs2_DIMER_dist_3|RARA-RARA_DIMER_dist_5|RARG-RARG_DIMER_dist_6"
GEMSTAT_experiment_info_df[33, 8] <- "Using monomer motifs for nuclear receptors and bounding for high interaction between self-pairs and low binding for each monomer"
i  <- 33
GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[i] <- paste0("experiment_", i)
names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)

####################################################################################################
####################################################################################################
# experiment 34
# simple model (same as exp 27) with no coop. changing the annotation thresholds for TFs (giving them as free parameters). using dimer TFs.
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[27,])
GEMSTAT_experiment_info_df[34, 1] <- "experiment_34"
GEMSTAT_experiment_info_df[34, 8] <- "setting annotation thresholds as free parameters to see what it does."

i  <- 34
GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[i] <- paste0("experiment_", i)
names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)


aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results[33:34], "[[", 7))

Linear_and_GEMSTAT_Accuracy_all <- cbind(Linear_and_GEMSTAT_Accuracy_all, aa)
colnames(Linear_and_GEMSTAT_Accuracy_all)[38:39] <- paste("experiment", c(33:34), sep = "_")
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all, las=2)

####################################################################################################
####################################################################################################
# experiment 35
# simple monomer model (same as exp 33) with no coop. changing na to 20
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[33,])
GEMSTAT_experiment_info_df[35, 1] <- "experiment_35"
GEMSTAT_experiment_info_df[35, 8] <- "monomer with na=20."
i  <- 35
GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[i] <- paste0("experiment_", i)
names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)

####################################################################################################
####################################################################################################
# experiment 36
# simple model (same as exp 34) with no coop. changing the annotation thresholds for TFs (giving them as free parameters). using dimer TFs.
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[34,])
GEMSTAT_experiment_info_df[36, 1] <- "experiment_36"
GEMSTAT_experiment_info_df[36, 8] <- "setting annotation thresholds as free parameters to see what it does, and na=20."

i  <- 36
GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[i] <- paste0("experiment_", i)
names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results[35:36], "[[", 7))

Linear_and_GEMSTAT_Accuracy_all <- cbind(Linear_and_GEMSTAT_Accuracy_all, aa)
colnames(Linear_and_GEMSTAT_Accuracy_all)[40:41] <- paste("experiment", c(35:36), sep = "_")
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all, las=2)
####################################################################################################
####################################################################################################
# experiment 37
# simple monomer model (same as exp 33) with no coop. widening the range of binding, activation and coop parameters
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[35,])
GEMSTAT_experiment_info_df[37, 1] <- "experiment_37"
GEMSTAT_experiment_info_df[37, 8] <- "widening the range of binding, activation and coop parameters."

####################################################################################################
####################################################################################################
# experiment 38
# simple monomer model (same as exp 37) with no coop. wide range. GR modified workplace
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[37,])
GEMSTAT_experiment_info_df[38, 1] <- "experiment_38"
GEMSTAT_experiment_info_df[38, 2] <- "GR"

####################################################################################################
####################################################################################################
# experiment 39
# simple monomer model (same as exp 37) with no coop. wide range. RARA modified workplace
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[37,])
GEMSTAT_experiment_info_df[39, 1] <- "experiment_39"
GEMSTAT_experiment_info_df[39, 2] <- "RAR"

####################################################################################################
####################################################################################################
# experiment 40
# simple monomer model (same as exp 37) with no coop. wide range. GR-RARA modified workplace
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[37,])
GEMSTAT_experiment_info_df[40, 1] <- "experiment_40"
GEMSTAT_experiment_info_df[40, 2] <- "Double"

for(i in 37:40){
  GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results[37:40], "[[", 7))
Linear_and_GEMSTAT_Accuracy_all <- cbind(Linear_and_GEMSTAT_Accuracy_all, aa)
colnames(Linear_and_GEMSTAT_Accuracy_all)[42:45] <- paste("experiment", c(37:40), sep = "_")
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all, las=2)

####################################################################################################
####################################################################################################
# experiment 41
# simple monomer model (same as exp 37) with no coop. wide range. alpha act upper bound 500 instead of 10e+6 in exp 37
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[37,])
GEMSTAT_experiment_info_df[41, 1] <- "experiment_41"
GEMSTAT_experiment_info_df[41, 8] <- "widening the range of binding, activation and coop parameters. alpha act upper is 500"

i  <- 41
GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[i] <- paste0("experiment_", i)
names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results[41], "[[", 7))
Linear_and_GEMSTAT_Accuracy_all <- cbind(Linear_and_GEMSTAT_Accuracy_all, aa)
colnames(Linear_and_GEMSTAT_Accuracy_all)[46] <- paste("experiment", 41, sep = "_")
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all, las=2)
aa <- mean_RMSE_calculator(GEMSTAT_based_on_linear_exp_results)

####################################################################################################
####################################################################################################
# experiment 42
# using the final parameters of experiment 29, fixing them. only changing annotation threshold (everything else is fixed)
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[29,])
GEMSTAT_experiment_info_df[42, 1] <- "experiment_42"
GEMSTAT_experiment_info_df[42, 8] <- "used the optimized parameters of exp 29 and fix all parameters except annotation thresholds"



####################################################################################################
####################################################################################################
# experiment 43
# using the final parameters of experiment 33, fixing them. only changing annotation threshold (everything else is fixed)
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[33,])
GEMSTAT_experiment_info_df[43, 1] <- "experiment_43"
GEMSTAT_experiment_info_df[43, 8] <- "used the optimized parameters of exp 33 and fix all parameters except annotation thresholds"

####################################################################################################
####################################################################################################
# experiment 44
# using the final parameters of experiment 41, fixing them. only changing annotation threshold (everything else is fixed)
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[41,])
GEMSTAT_experiment_info_df[44, 1] <- "experiment_44"
GEMSTAT_experiment_info_df[44, 8] <- "used the optimized parameters of exp 41 and fix all parameters except annotation thresholds"

for(i in 42:44){
  GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results[42:44], "[[", 7))
Linear_and_GEMSTAT_Accuracy_all <- cbind(Linear_and_GEMSTAT_Accuracy_all, aa)
colnames(Linear_and_GEMSTAT_Accuracy_all)[47:49] <- paste("experiment", c(42:44), sep = "_")
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all, las=2)
abline(h=seq(0, 1, 0.01), col=2, lty=2, lwd=0.5)
aa <- mean_RMSE_calculator(GEMSTAT_based_on_linear_exp_results)

ParameterHeatMap(GEMSTAT_based_on_linear_exp_pars[[1]]$binding)
heatmap
####################################################################################################
####################################################################################################
# experiment 45
# same as experiment 33
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[33,])
GEMSTAT_experiment_info_df[45, 1] <- "experiment_45"
GEMSTAT_experiment_info_df[45, 8] <- "same as experiment 33 but using one qBTM per gene and widening the range of qBTM to c(1e-5, 5)"

i  <- 45
GEMSTAT_based_on_linear_exp_results[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results)[i] <- paste0("experiment_", i)
names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results[i], "[[", 7))
Linear_and_GEMSTAT_Accuracy_all <- cbind(Linear_and_GEMSTAT_Accuracy_all, aa)
colnames(Linear_and_GEMSTAT_Accuracy_all)[i + 5] <- paste("experiment", i, sep = "_")
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all, las=2)
abline(h=seq(0, 1, 0.01), col = 2, lwd = 0.5, lty = 4)
aa <- mean_RMSE_calculator(GEMSTAT_based_on_linear_exp_results)

aa <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)

View(GEMSTAT_based_on_linear_exp_results[[45]]$Prediction_raw_list$`62`)
#####
GEMSTAT_based_on_linear_exp_results_sigmoid <- list()
for(i in 9:45){
  GEMSTAT_based_on_linear_exp_results_sigmoid[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
  
}

which.max(GEMSTAT_based_on_linear_exp_results_sigmoid[[41]]$Accuracy_All)
View(GEMSTAT_based_on_linear_exp_results_sigmoid[[41]]$Prediction_round_list[[145]])

aa <- mean_RMSE_calculator(GEMSTAT_based_on_linear_exp_results_sigmoid)
aa <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 29, "/Outputs"), .fold_change = T, .round_thresh = 0.4)

aa2 <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 29, "/Outputs"), .fold_change = T, .round_thresh = 0.3)

aa3 <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 29, "/Outputs"), .fold_change = T, .round_thresh = 0.2)
aa4 <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 29, "/Outputs"), .fold_change = T, .round_thresh = 0.1)


boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid[[29]]$Accuracy_All, aa$Accuracy_All, aa2$Accuracy_All, aa3$Accuracy_All, aa4$Accuracy_All))
GEMSTAT_based_on_linear_exp_results_sigmoid_thr1 <- list()
for(i in 1:45){
  GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T, .round_thresh = 0.1)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_thr1, "[[", 7))
boxplot.matrix(cbind(Linear_and_GEMSTAT_Accuracy_all[,1:5], aa), las=2)
abline(h=seq(0,1, 0.01), col=2, lty=4, lwd=0.5)

aa <- GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[31]]$Prediction_raw_list[[129]]
sum(aa > -0.05 & aa < 0.05, na.rm = T)
sum(aa >= 0.05, na.rm = T)

GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[31]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 31, "/Outputs"), .fold_change = T, .round_thresh = 0.1)
aa <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", 31, "/Outputs"), .fold_change = T, .round_thresh = 0.05)
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[31]]$Accuracy_All, aa$Accuracy_All))


plotExpression(expMat = aa$Prediction_round_list[[129]], .dendrogram = "none", .Rowv = F, .Colv = F, colorPoints = c(-2, -0.8, -0.2, 0.2 ,1.5),filename = "pred.png")
####################################################################################################
# for each of the original experiments plot raw values predicted grouped by their real value

GEMSTAT_based_on_linear_exp_results_sigmoid_prediction_bygroup <- list()

for(i in 1:length(GEMSTAT_based_on_linear_exp_results_sigmoid)){
  aa_table_list <- lapply(GEMSTAT_based_on_linear_exp_results_sigmoid[[i]]$GroundTruth_List, table)
  aa <- max(unlist(aa_table_list))
  GEMSTAT_based_on_linear_exp_results_sigmoid_prediction_bygroup[[i]] <- matrix(nrow = aa,
                                                                                ncol= (3 *length(GEMSTAT_based_on_linear_exp_results_sigmoid[[i]]$GroundTruth_List)))
  for(j in 1:length(GEMSTAT_based_on_linear_exp_results_sigmoid[[i]]$GroundTruth_List)){
    aagg <- GEMSTAT_based_on_linear_exp_results_sigmoid[[i]]$Prediction_raw_list[[j]][GEMSTAT_based_on_linear_exp_results_sigmoid[[i]]$GroundTruth_List[[j]] == as.integer(names(aa_table_list[[j]])[1])]
    GEMSTAT_based_on_linear_exp_results_sigmoid_prediction_bygroup[[i]][1:aa_table_list[[j]][1], 3*(j-1) + 1] <-  aagg[!is.na(aagg)]
    aagg <- GEMSTAT_based_on_linear_exp_results_sigmoid[[i]]$Prediction_raw_list[[j]][GEMSTAT_based_on_linear_exp_results_sigmoid[[i]]$GroundTruth_List[[j]] == as.integer(names(aa_table_list[[j]])[2])]
    GEMSTAT_based_on_linear_exp_results_sigmoid_prediction_bygroup[[i]][1:aa_table_list[[j]][2], 3*(j-1) + 2] <-  aagg[!is.na(aagg)]
    aagg <- GEMSTAT_based_on_linear_exp_results_sigmoid[[i]]$Prediction_raw_list[[j]][GEMSTAT_based_on_linear_exp_results_sigmoid[[i]]$GroundTruth_List[[j]] == as.integer(names(aa_table_list[[j]])[3])]
    GEMSTAT_based_on_linear_exp_results_sigmoid_prediction_bygroup[[i]][1:aa_table_list[[j]][3], 3*j] <-        aagg[!is.na(aagg)]
    
  }
}
boxplot.matrix(GEMSTAT_based_on_linear_exp_results_sigmoid_prediction_bygroup[[31]][,1:90], col=rep(c(2, 3, 4), 40), outline=T)

####################################################################################################
####################################################################################################
# experiment 46
# same as experiment 33 but with ER-GR cooperativity
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[33,])
GEMSTAT_experiment_info_df[46, 1] <- "experiment_46"
GEMSTAT_experiment_info_df[46, 6] <- "AR-AR_DIMER_dist_4|ER-ER_DIMER_dist_4|NR3C1-NR3C1_DIMER_dist_4|PGRs1-PGRs2_DIMER_dist_4|RARA-RARA_DIMER_dist_6|RARG-RARG_DIMER_dist_7|ER-NR3C1_SIMPLE_dist_50"

####################################################################################################
####################################################################################################
# experiment 47
# same as experiment 33 but with ER-RARA cooperativity
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[33,])
GEMSTAT_experiment_info_df[47, 1] <- "experiment_47"
GEMSTAT_experiment_info_df[47, 6] <- "AR-AR_DIMER_dist_4|ER-ER_DIMER_dist_4|NR3C1-NR3C1_DIMER_dist_4|PGRs1-PGRs2_DIMER_dist_4|RARA-RARA_DIMER_dist_6|RARG-RARG_DIMER_dist_7|ER-RARA_SIMPLE_dist_50"

####################################################################################################
####################################################################################################
# experiment 48
# same as experiment 33 but with ER-RARA and ER-GR cooperativity
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[33,])
GEMSTAT_experiment_info_df[48, 1] <- "experiment_48"
GEMSTAT_experiment_info_df[48, 6] <- "AR-AR_DIMER_dist_4|ER-ER_DIMER_dist_4|NR3C1-NR3C1_DIMER_dist_4|PGRs1-PGRs2_DIMER_dist_4|RARA-RARA_DIMER_dist_6|RARG-RARG_DIMER_dist_7|ER-NR3C1_SIMPLE_dist_50|ER-RARA_SIMPLE_dist_50"

####################################################################################################
####################################################################################################
# experiment 49
# same as experiment 33 but with coop distances increased by one (since in GEMSTAT code it uses < not <= for distance threshold)
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[33,])
GEMSTAT_experiment_info_df[49, 1] <- "experiment_49"
GEMSTAT_experiment_info_df[49, 6] <- "AR-AR_DIMER_dist_4|ER-ER_DIMER_dist_4|NR3C1-NR3C1_DIMER_dist_4|PGRs1-PGRs2_DIMER_dist_4|RARA-RARA_DIMER_dist_6|RARG-RARG_DIMER_dist_7"


for(i in 46:49){
  GEMSTAT_based_on_linear_exp_results_sigmoid[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
  GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T, .round_thresh = 0.1)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_thr1)  <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}




aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_thr1, "[[", 7))
colnames(aa) <- c(paste("Exp", c(1:ncol(aa)), sep="_"))
boxplot.matrix(cbind(Linear_and_GEMSTAT_Accuracy_all[,1:5], aa), las=2, main="threshold 0.1")
abline(h=seq(0,1, 0.01), col=2, lty=4, lwd=0.5)

boxplot_grouped_by_label(GEMSTAT_output_Reader_multiEnh_ensemble_output = GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[49]], 
                         export_plot=T, 
                         filename="boxplot_by_label_thr1_49.png",
                         return_mat=F )



aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid, "[[", 7))
colnames(aa) <- c(paste("Exp", c(1:ncol(aa)), sep="_"))
boxplot.matrix(cbind(Linear_and_GEMSTAT_Accuracy_all[,1:5], aa), las=2, main="threshold 0.5" )
abline(h=seq(0,1, 0.01), col=2, lty=4, lwd=0.5)

GEMSTAT_based_on_linear_exp_results_sigmoid_thr05 <- list()
for(i in 1:49){
  GEMSTAT_based_on_linear_exp_results_sigmoid_thr05[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                    .fold_change = T, .round_thresh = 0.05)
}
aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_thr05, "[[", 7))
boxplot.matrix(cbind(Linear_and_GEMSTAT_Accuracy_all[,1:5], aa), las=2)
abline(h=seq(0,1, 0.01), col=2, lty=4, lwd=0.5)


####################################################################################################
####################################################################################################
# preparing RNA-seq test set for GEMSTAT
aa <- match(rownames(RNAseq_Test_1_geneExp), rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3))
RNAseq_Test_1_geneExp_common3Enh <- RNAseq_Test_1_geneExp[!is.na(aa),]
RNAseq_Test_1_geneExp_common3Enh_indmatch <- match(rownames(RNAseq_Test_1_geneExp_common3Enh), rownames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52))

TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer <-rbind(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun_ER01[1:13, ],
                                                                              TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun_ER01[13, ],
                                                                              TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun_ER01[14:19, ])
rownames(TF.Exp.Shrinked.RNAseq_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer)[13:14] <- c("PGR_s1", "PGR_s2")

####################################################################################################
####################################################################################################
# experiment 50
# using optimized parameters of experiment 49 on RNA seq test dataset. no training
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[49,])
GEMSTAT_experiment_info_df[50, 1] <- "experiment_50"
GEMSTAT_experiment_info_df[50, 8] <- "This is a test on RNA-seq data using optimized parameters of experiment 49"
GEMSTAT_experiment_info_df[50, 9] <- "no"
i <- 50
GEMSTAT_based_on_linear_exp_results_sigmoid[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T, .round_thresh = 0.1)
GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
names(GEMSTAT_based_on_linear_exp_results_sigmoid)[i] <- paste0("experiment_", i)
names(GEMSTAT_based_on_linear_exp_results_sigmoid_thr1)[i]  <- paste0("experiment_", i)
names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)

boxplot_grouped_by_label(GEMSTAT_output_Reader_multiEnh_ensemble_output = GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[50]], 
                         export_plot=T, 
                         filename="boxplot_by_label_thr1_50.png",
                         return_mat=F )
####################################################################################################
GEMSTAT_experiment_info_df <- cbind(GEMSTAT_experiment_info_df,
                                    c(rep(52, 28), rep(39, 21), 31))
colnames(GEMSTAT_experiment_info_df)[10] <- "starting number of genes"
####################################################################################################
####################################################################################################
# experiment 51
# same as exp 27 but not allowing removal of bad prediced genes in linear
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[27,])
GEMSTAT_experiment_info_df[51, 1] <- "experiment_51"
GEMSTAT_experiment_info_df[51, 9] <- "no"

####################################################################################################
####################################################################################################
# experiment 52
# same as exp 29 but not allowing removal of bad prediced genes in linear
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[29,])
GEMSTAT_experiment_info_df[52, 1] <- "experiment_52"
GEMSTAT_experiment_info_df[52, 9] <- "no"
####################################################################################################
####################################################################################################
# experiment 53
# same as exp 49 but woth 52 genes and not allowing removal of bad prediced genes in linear
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[49,])
GEMSTAT_experiment_info_df[53, 1] <- "experiment_53"
GEMSTAT_experiment_info_df[53, 9] <- "no"
GEMSTAT_experiment_info_df[53, 10] <- "52"
####################################################################################################
####################################################################################################
# experiment 54
# same as exp 49 but not allowing removal of bad prediced genes in linear
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[49,])
GEMSTAT_experiment_info_df[54, 1] <- "experiment_54"
GEMSTAT_experiment_info_df[54, 9] <- "no"
####################################################################################################
####################################################################################################
# experiment 55
# same as 53 but not allowing removal of bad prediced genes in linear
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[53,])
GEMSTAT_experiment_info_df[55, 1] <- "experiment_55"
GEMSTAT_experiment_info_df[55, 9] <- "yes"

# reading exp51-55
for(i in 51:55){
  GEMSTAT_based_on_linear_exp_results_sigmoid[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T)
  GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"), .fold_change = T, .round_thresh = 0.1)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_thr1)[i]  <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_thr1, "[[", 7))
colnames(aa) <- c(paste("Exp", c(1:ncol(aa)), sep="_"))
boxplot.matrix(cbind(Linear_and_GEMSTAT_Accuracy_all[,1:5], aa), las=2, main="threshold 0.1")
abline(h=seq(0,1, 0.01), col=2, lty=4, lwd=0.5)

boxplot_grouped_by_label(GEMSTAT_output_Reader_multiEnh_ensemble_output = GEMSTAT_based_on_linear_exp_results_sigmoid_thr1[[53]], 
                         export_plot=T, 
                         filename="boxplot_by_label_thr1_53.png",
                         return_mat=F )

########## Computing the results with automatic threshold setting
GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh <- list()
for(i in 1:55){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
}

aa1 <- cbind(Sim_Ann_weighted_148_restart_random_pred$Precision,
             Sim_Ann_weighted_148_restart_round_accuracy_automatic,
             do.call(cbind, Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat_accuracy),
             Optim_Greedy_conc_modif_Double_Results_automated_thresh_accuracy)
colnames(aa1) <- c("Random","WT","ARp", "ARn", "GRp", "GRn", "PGRp",
                   "PGRn", "RARAp","RARAn","RARGp", "RARGn", "RXRAp","RXRAn","GRn_RARAp")
GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[50]]$Accuracy_All <- c(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[50]]$Accuracy_All, c(0.33, 0.33))
aa2 <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh, "[[", 7))
colnames(aa2) <- c(paste("Exp", c(1:ncol(aa2)), sep="_"))

Linear_and_GEMSTAT_Accuracy_all_automatized <- cbind(aa1, aa2)
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all_automatized, las=2, main="threshold automatic" )
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)

# View(GEMSTAT_experiment_info_df[46:55,c(1,6, 8, 9,10)])
# aa_srt <- sort(aa2[,46], decreasing = T, index.return=T)$ix
# View(GEMSTAT_based_on_linear_exp_pars[[46]]$binding[aa_srt[1:10], ])
# View(GEMSTAT_based_on_linear_exp_pars[[46]]$alpha_effective[aa_srt[1:10], ])
# GEMSTAT_based_on_linear_exp_pars[[46]]$coop[aa_srt[1:10]]
##############################################################################################################
##############################################################################################################
aa <- PerformanceHeattMap_General(real_exp_mat =GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$GroundTruth_List[[1]] ,
                                  prediction_mat_list = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Prediction_round_list,
                                  .Colv = T, .Rowv = T, .dendrogram = "both",exportplot = T,
                                  filename = "PerformanceHeatMap_model_54.png",
                                  .RowSideColors = character(0))
##############################################################################################################
# parameters of experiment number 54
aa_sr <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Accuracy_All, decreasing = T, index.return=T)$ix

aa_all_coop <- do.call(rbind, GEMSTAT_based_on_linear_exp_pars[[54]]$coop)

aaa_breaks = c(
  seq(1e-3, 1, length=5),    
  seq(1.01, 5, length=5),           
  seq(5.01, 10, length=5)
  # seq(1.01, 10, length = 5),
  # seq(10.01, 100, length = 5),
  # seq(100.01, 500, length = 5)
  #seq(500.01, 2000, length = 20)
)
png(filename = "parameters.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # smaller font size
heatmap.2(x = cbind(GEMSTAT_based_on_linear_exp_pars[[54]]$binding[aa_sr[1:15],]/25,
                    GEMSTAT_based_on_linear_exp_pars[[54]]$alpha_effective[aa_sr[1:15],],
                    (aa_all_coop[aa_sr[1:15], ] - 100)/40 ),
          Rowv = T, Colv = F,
          dendrogram = "row",
          #rowsep = aa_rowsep-1, sepwidth = c(4,4), sepcolor = "green",
          #RowSideColors = aa_rowsidecol,
          trace="none", na.rm = T, symbreaks=F , symkey=F, symm = F,
          breaks = aaa_breaks,
          col = colorRampPalette(c("green", "black", "red"))(n = 14), margins = c(7, 7))
dev.off()
##############################################################################################################
##############################################################################################################
# starting a new set of experiments.

# 20 experiments setting the concentration of each TF to 0 one by one. no optimization. to see
#  what parameters are the models actually using. 
# Creating the R script to create the inputs and job files.
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 22:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp54.R")
  aa[[i]][14] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  if(i == 21){
    aa[[i]][15] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][15] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][15] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }

  aa[[i]][18] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 55 + i, ",")
  aa[[i]][37] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][67] <-  "                      na = 0,"
  aa[[i]][75] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=T, create_jobFile=T"
  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=54,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=54,")
  aa[[i]] <- c(aa[[i]], "                              coop=54, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=54,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=54,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=54, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=54,")
  aa[[i]] <- c(aa[[i]], "                              Params=54, ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=T,")
  aa[[i]] <- c(aa[[i]], "                              seq=54,")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=54,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=54, ")
  aa[[i]] <- c(aa[[i]], "                              wights=54, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 55+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  55 + i, ".R"))
}
setwd("..")

###### Writing information for the experiments
for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[54,])
  GEMSTAT_experiment_info_df[55+i, 8] <- paste("in silico knockdown of TF", names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[55+i, 1] <- paste0("experiment_", 55+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[54,])
GEMSTAT_experiment_info_df[55+21, 8] <- "in silico knockdown of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[55+21, 1] <- paste0("experiment_", 55+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[54,])
GEMSTAT_experiment_info_df[55+22, 8] <- "in silico knockdown of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[55+22, 1] <- paste0("experiment_", 55+22)

###### ###### ###### reading the results for insilico knockdowns (no training) ###### ###### ###### 

for(i in 56:77){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa1 <- cbind(Sim_Ann_weighted_148_restart_random_pred_39$Precision,
             Sim_Ann_weighted_148_restart_round_accuracy_automatic_39,
             do.call(cbind, Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat_accuracy_39),
             Optim_Greedy_conc_modif_Double_Results_automated_thresh_accuracy_39)
colnames(aa1) <- c("Random","WT","ARp", "ARn", "GRp", "GRn", "PGRp",
                   "PGRn", "RARAp","RARAn","RARGp", "RARGn", "RXRAp","RXRAn","GRn_RARAp")
aa2 <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh, "[[", 7))
colnames(aa2) <- c(paste("Exp", c(1:ncol(aa2)), sep="_"))
par(mar = c(8,4,4,4))
aa2_my8 <- aa2[,c(51, 52, 27, 29, 53, 54, 55, 49)]
colnames(aa2_my8) <- c("dimer_52_no", "dimer_39_no", "dimer_52_yes", "dimer_39_yes", "monomer_52_no", "monomer_39_no", "monomer_52_yes", "monomer_39_yes")
boxplot.matrix(aa2_my8, las=2)


Linear_and_GEMSTAT_Accuracy_all_automatized_linear39 <- cbind(aa1, aa2)
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all_automatized_linear39, las=2, main="threshold automatic" )
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
colnames(Linear_and_GEMSTAT_Accuracy_all_automatized_linear39)[71:90] <- names(TF.motifs.Shrinked.halfsites.count)
colnames(Linear_and_GEMSTAT_Accuracy_all_automatized_linear39)[91] <- "PGR_s1 + PGR_s2"
colnames(Linear_and_GEMSTAT_Accuracy_all_automatized_linear39)[92] <- "JUN_1 + JUN_2"
boxplot.matrix(Linear_and_GEMSTAT_Accuracy_all_automatized_linear39[, c(c(1:15), 69, c(71:92))], las=2, main="threshold automatic" )
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)

### fix the threshold before reading
GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh_fixed <- list()
for(i in 56:77){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh_fixed[[i-55]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$threshold)
  # GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  # names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  # names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh_fixed, "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:15]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
### find out the composition of datapoints (-1, 0, 1) within the points that were predicted correctly in original model but
###  are now misclassified in each of the KD experiments.



aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:15]
aa_miscalssified_list <- list()
aa_miscalssified_list_sep <- list()
aa_corcalssified_list <- list()
aa_ratio_list <- list()
aa_ratio_list_sep <- list()
for(aa_cmodif in 1:22){
  aa_miscalssified_list[[aa_cmodif]] <- matrix(nrow = length(aa_goodind), ncol = 3)
  colnames(aa_miscalssified_list[[aa_cmodif]]) <- c(-1, 0, 1)
  aa_corcalssified_list[[aa_cmodif]] <- matrix(nrow = length(aa_goodind), ncol = 3)
  colnames(aa_corcalssified_list[[aa_cmodif]]) <- c(-1, 0, 1)
  aa_ratio_list[[aa_cmodif]]         <- matrix(nrow = length(aa_goodind), ncol = 3)
  colnames(aa_ratio_list[[aa_cmodif]]) <- c(-1, 0, 1)
  
  aa_miscalssified_list_sep[[aa_cmodif]] <- matrix(nrow = length(aa_goodind), ncol = 6)
  aa_ratio_list_sep[[aa_cmodif]]         <- matrix(nrow = length(aa_goodind), ncol = 6)
  colnames(aa_miscalssified_list_sep[[aa_cmodif]]) <- c("-1_0", "-1_1", "0_-1", "0_1", "1_0", "1_-1")
  colnames(aa_ratio_list_sep[[aa_cmodif]]) <- colnames(aa_miscalssified_list_sep[[aa_cmodif]])
  
  rownames(aa_miscalssified_list[[aa_cmodif]]) <- aa_goodind
  rownames(aa_miscalssified_list_sep[[aa_cmodif]]) <- aa_goodind
  rownames(aa_corcalssified_list[[aa_cmodif]]) <- aa_goodind
  rownames(aa_ratio_list[[aa_cmodif]])<- aa_goodind
  rownames(aa_ratio_list_sep[[aa_cmodif]])<- aa_goodind
  
  for(aa_cmodel in 1:length(aa_goodind)){
    aa_gt <- GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$GroundTruth_List[[aa_goodind[aa_cmodel]]]
    aa_afterKD <- GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[55 + aa_cmodif]]$Prediction_round_list[[aa_goodind[aa_cmodel]]]
    aa_beforeKD <- GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Prediction_round_list[[aa_goodind[aa_cmodel]]]
    aa_missed <- aa_gt[aa_beforeKD == aa_gt & aa_afterKD != aa_gt]
    aa_allcor <- aa_gt[aa_beforeKD == aa_gt]
    
    aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 1] <- sum(aa_allcor == -1, na.rm = T)
    aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 1] <- sum(aa_missed == -1, na.rm = T)
    aa_ratio_list[[aa_cmodif]][aa_cmodel, 1] <- aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 1]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 1]
   
    aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 1] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == -1 & aa_afterKD == 0, na.rm = T)
    aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 2] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == -1 & aa_afterKD == 1, na.rm = T)
    aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 1] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 1]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 1]
    aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 1] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 2]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 1]
    
    aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 2] <- sum(aa_allcor ==  0, na.rm = T)
    aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 2] <- sum(aa_missed ==  0, na.rm = T)
    aa_ratio_list[[aa_cmodif]][aa_cmodel, 2] <- aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 2]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 2]
    
        
    aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 3] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == 0 & aa_afterKD == -1, na.rm = T)
    aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 4] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == 0 & aa_afterKD == 1, na.rm = T)
    aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 3] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 3]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 2]
    aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 4] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 4]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 2]
    
    aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 3] <- sum(aa_allcor ==  1, na.rm = T)
    aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 3] <- sum(aa_missed ==  1, na.rm = T)
    aa_ratio_list[[aa_cmodif]][aa_cmodel, 3] <- aa_miscalssified_list[[aa_cmodif]][aa_cmodel, 3]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 3]

    aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 5] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == 1 & aa_afterKD == 0, na.rm = T)
    aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 6] <- sum(aa_beforeKD == aa_gt & aa_afterKD != aa_gt & aa_gt == 1 & aa_afterKD == -1, na.rm = T)
    aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 5] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 5]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 3]
    aa_ratio_list_sep[[aa_cmodif]][aa_cmodel, 6] <- aa_miscalssified_list_sep[[aa_cmodif]][aa_cmodel, 6]/aa_corcalssified_list[[aa_cmodif]][aa_cmodel, 3]
    
  }
}
names(aa_miscalssified_list) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")
names(aa_corcalssified_list) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")
names(aa_ratio_list) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")
names(aa_ratio_list_sep) <- names(aa_ratio_list) 
names(aa_miscalssified_list_sep) <-  names(aa_ratio_list) 

par(mfrow=c(4, 6), mar = c(3, 3, 4, 1))
for(i in 1:22){
  boxplot.matrix(aa_ratio_list_sep[[i]],
                 main = names(aa_ratio_list)[i],
                 ylim = c(0, 1), 
                 col = c(2,2, 3,3, 4,4), 
                 las = 2)
}
par(mfrow=c(4, 6), mar = c(3, 3, 4, 1))
for(i in 1:22){
  boxplot.matrix(aa_ratio_list[[i]],
                 main = names(aa_ratio_list)[i],
                 ylim = c(0, 1), 
                 col = c(2, 3, 4), 
                 las = 2)
}

par(mfrow = c(1, 1))
aa_tab <- table(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$GroundTruth_List[[aa_goodind[aa_cmodel]]])
aa_cor <- aa_corcalssified_list[[1]]
for(i in 1:nrow(aa_cor)){
  aa_cor[i, ] <- aa_cor[i, ]/aa_tab
}

boxplot.matrix(aa_cor)
plot(aa_cor[1,], type = "l", ylim = c(0.4,0.9))
for(i in 2:nrow(aa_cor)){
  lines(aa_cor[i,], col=col_vector[i])
}

########################################################################################################################
########################################################################################################################
# 12 experiments corresponding to the TF conc modification
aa_names <- c("ARp", "ARn", "GRp", "GRn", "PGRp", "PGRn",
             "RARAp","RARAn","RARGp", "RARGn", "RXRAp","RXRAn")
for(i in 1:12){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[54,])
  GEMSTAT_experiment_info_df[77+i, 8] <- "concentration modification experiment"
  GEMSTAT_experiment_info_df[77+i, 2] <- aa_names[i]
  GEMSTAT_experiment_info_df[77+i, 1] <- paste0("experiment_", 77+i)
}


# Creating the R script to create the inputs and job files.
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
aa_modif <- list(1, -1 ,10, -10, c(13, 14), c(-13, -14), 15,-15, 16,-16, 18,-18)

for(i in 1:12){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp54.R")
  aa[[i]][14] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  for(j in 1:length(aa_modif[[i]])){
    if(aa_modif[[i]][j] > 0){
      aa[[i]][15 + j-1] <- paste0('aa_myTFexpmat[',aa_modif[[i]][j],', ] <- rep(c(0,1), 45)')
    }else{
      aa[[i]][15 + j-1] <- paste0('aa_myTFexpmat[',-aa_modif[[i]][j],', ] <- rep(c(1,0), 45)')
    }
  }
  
  
  aa[[i]][18] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 77 + i, ",")
  aa[[i]][27] <- paste0("                      model_evaluations = Optim_Greedy_conc_modif_Ensemble_Results_eval[[", i, "]],")
  aa[[i]][30] <- paste0("                      model_parameters = do.call(rbind, lapply(Optim_Greedy_conc_modif_Ensemble_Results[[", i, ']], "[[", 3)),')
  aa[[i]][37] <- "                      TF_expression_mat = aa_myTFexpmat,"
  #aa[[i]][67] <-  "                      na = 0,"
  aa[[i]][75] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=T,create_seq=T, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=T, create_jobFile=T"
  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=54,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=54,")
  aa[[i]] <- c(aa[[i]], "                              coop=54, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=54,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=54,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=54, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=54,")
  aa[[i]] <- c(aa[[i]], "                              Params=numeric(0), ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=F,")
  aa[[i]] <- c(aa[[i]], "                              seq=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=54,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=54, ")
  aa[[i]] <- c(aa[[i]], "                              wights=54, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 77+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  77 + i, ".R"))
}
setwd("..")
######
# reading the concentration modification results

for(i in c(79, 89)){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa1 <- cbind(Sim_Ann_weighted_148_restart_random_pred_39$Precision,
             Sim_Ann_weighted_148_restart_round_accuracy_automatic_39)
colnames(aa1) <- c("Random","WT_linear")
aa2 <- do.call(cbind, Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat_accuracy_39)
colnames(aa2) <- paste(c("ARp", "ARn", "GRp", "GRn", "PGRp", "PGRn",
                         "RARAp","RARAn","RARGp", "RARGn", "RXRAp","RXRAn"),
                       "linear", sep = "_")
aa3 <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh, "[[", 7))
colnames(aa3) <- c(paste("Exp", c(1:ncol(aa3)), sep="_"))
aa4 <- aa3[, c(54, c(78:89))]
colnames(aa4) <- c("WT_GEMSTAT", paste(c("ARp", "ARn", "GRp", "GRn", "PGRp",
                                         "PGRn", "RARAp","RARAn","RARGp", "RARGn",
                                         "RXRAp","RXRAn"), "GEMSTAT", sep="_"))
par(mar = c(10,4,4,4))
boxplot.matrix(cbind(aa1, aa2, aa4), las= 2, main="concentration_modif accuracy", cex=0.5)
abline(h=seq(0, 1, 0.05), col=2, lty=4, lwd=0.5)
########################################################################################################################
########################################################################################################################
# cooperativity
# running 190 experiments for all cooperativity pairs
# only try for the top 15 models.

# write coop parameters, save it to the input construction workplace

Exp_90_278_coop_parameter_combination <- coop_parameter_creator_comb(TF_names=names(TF.motifs.Shrinked.halfsites.count),
                                                                     TF_index=c(1:20), 
                                                                     excluded_pairs=list(c(13, 14)),
                                                                     new_weight_range=c(0.001, 100),
                                                                     new_initial_weight = 1,
                                                                     new_coop_type="SIMPLE",
                                                                     new_coop_dist=50,
                                                                     new_coop_orientation=c(0, 0),
                                                                     base_mat=rbind(c(1, 1), c(3, 3), c(10, 10), c(13, 14), c(15, 15), c(16, 16)),
                                                                     base_weight_range=cbind(c(100, 100, 100, 100, 100, 100),
                                                                                             c(500, 500, 500, 500, 500, 500)),
                                                                     base_initial_coop_weight=c(300, 300, 300, 300, 300, 300),
                                                                     base_coop_type=rep("DIMER", 6),
                                                                     base_coop_dist=c(4, 4, 4, 4, 6, 7),
                                                                     base_coop_orientation=rbind(c(1, -1), c(1, -1), c(1, -1), c(1, 1), c(1, 1), c(1, 1))
)

# Creating the R script to create the inputs and job files.
setwd("Seeded_GEMSTAT_ens/")

aa <- list()
for(i in 2:length(Exp_90_278_coop_parameter_combination$Coop_mat_list)){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp54.R")
  aa[[i]][8] <-  paste0("my_coop_tf_mat <-        Exp_90_278_coop_parameter_combination$Coop_mat_list[[", i,"]]")
  aa[[i]][9] <-  paste0("my_coop_weight_range <-  Exp_90_278_coop_parameter_combination$Weight_range_list[[", i,"]]")
  aa[[i]][10] <- paste0("my_initial_coop_weight<- Exp_90_278_coop_parameter_combination$Initial_weight_list[[", i,"]]")
  aa[[i]][11] <- paste0("my_coop_type <-          Exp_90_278_coop_parameter_combination$Coop_type_list[[", i,"]]")
  aa[[i]][12] <- paste0("my_coop_dist <-          Exp_90_278_coop_parameter_combination$Coop_dist_list[[", i,"]]")
  aa[[i]][13] <- paste0("my_coop_orientation <-   Exp_90_278_coop_parameter_combination$Coop_orientation_list[[", i,"]]")
  aa[[i]][18] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 89 + i, ",")
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  89 + i, ".R"))
}
setwd("..")
########################################################################################################################
# update the info dataframe
aa_pair_index <- coop_parameter_creator_comb(TF_names=names(TF.motifs.Shrinked.halfsites.count),
                                             TF_index=c(1:20), 
                                             excluded_pairs=list(c(13, 14)),
                                             new_weight_range=c(0.001, 100),
                                             new_initial_weight = 1,
                                             new_coop_type="SIMPLE",
                                             new_coop_dist=50,
                                             new_coop_orientation=c(0, 0),
                                             base_mat=rbind(c(1, 1), c(3, 3), c(10, 10), c(13, 14), c(15, 15), c(16, 16)),
                                             base_weight_range=cbind(c(100, 100, 100, 100, 100, 100),
                                                                     c(500, 500, 500, 500, 500, 500)),
                                             base_initial_coop_weight=c(300, 300, 300, 300, 300, 300),
                                             base_coop_type=rep("DIMER", 6),
                                             base_coop_dist=c(4, 4, 4, 4, 6, 7),
                                             base_coop_orientation=rbind(c(1, -1), c(1, -1), c(1, -1), c(1, 1), c(1, 1), c(1, 1)),
                                             only_pair_index = T)

for(i in 1:189){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[54,])
  GEMSTAT_experiment_info_df[89+i, 8] <- paste0("coop added: ",
                                                names(TF.motifs.Shrinked.halfsites.count)[aa_pair_index[i, 1]], "_",
                                                names(TF.motifs.Shrinked.halfsites.count)[aa_pair_index[i, 2]])
  GEMSTAT_experiment_info_df[89+i, 6] <- paste0(GEMSTAT_experiment_info_df[89+i, 6], 
                                                "|",
                                                names(TF.motifs.Shrinked.halfsites.count)[aa_pair_index[i, 1]],
                                                "-",
                                                names(TF.motifs.Shrinked.halfsites.count)[aa_pair_index[i, 2]], 
                                                "_SIMPLE_dist_50")
  GEMSTAT_experiment_info_df[89+i, 1] <- paste0("experiment_", 89+i)
}
########################################################################################################################
# reading the coop results

for(i in 90:112){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
for(i in 130:152){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
for(i in 153:169){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

for(i in 170:182){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
for(i in 183:193){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
for(i in 194:203){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
for(i in 210:229){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

for(i in 250:274){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa1 <- cbind(Sim_Ann_weighted_148_restart_random_pred_39$Precision,
             Sim_Ann_weighted_148_restart_round_accuracy_automatic_39)
colnames(aa1) <- c("Random","WT_linear")
aa2 <- do.call(cbind, Optim_Greedy_conc_modif_Ensemble_Results_eval_real_value_automat_accuracy_39)
colnames(aa2) <- paste(c("ARp", "ARn", "GRp", "GRn", "PGRp", "PGRn",
                         "RARAp","RARAn","RARGp", "RARGn", "RXRAp","RXRAn"),
                       "linear", sep = "_")
aa3 <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(54, c(c(90:112), c(130:203), c(210:229), c(250:274)))], "[[", 7))
#colnames(aa3) <- c(paste("Exp", c(1:ncol(aa3)), sep="_"))
#aa4 <- aa3[, c(54, c(c(90:112), c(130:152), c(170:182)))]
aa_names <- character(0)
aa_ind <- c((c(90:112) - 89), (c(130:203) - 89), (c(210:229) - 89), (c(250:274) - 89))
for(i in 1:length(aa_ind)){
  aa_names <- c(aa_names, paste(names(TF.motifs.Shrinked.halfsites.count)[aa_pair_index[aa_ind[i], 1]],
                                names(TF.motifs.Shrinked.halfsites.count)[aa_pair_index[aa_ind[i], 2]], sep = "--"))
}
colnames(aa3) <- c("WT_GEMSTAT", aa_names)
par(mar = c(10,4,4,4))
boxplot.matrix(cbind(aa1, aa2, aa3), las= 2, main="coop_one_by_one accuracy", cex=0.5)
abline(h=seq(0, 1, 0.05), col=2, lty=4, lwd=0.5)

aa_ind[which(aa_names == "ESR1--PGR_s2")] + 89
sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[137]]$Accuracy_All, decreasing = T, index.return=T)$ix
GEMSTAT_based_on_linear_exp_pars[[137]]$coop[[25]]
GEMSTAT_based_on_linear_exp_pars[[137]]$binding[25,]
GEMSTAT_based_on_linear_exp_pars[[137]]$alpha_effective[25,]


aa_ind[which(aa_names == "CEBPB--FOXA1")] + 89
sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[110]]$Accuracy_All, decreasing = T, index.return=T)$ix
GEMSTAT_based_on_linear_exp_pars[[110]]$coop[[25]]
GEMSTAT_based_on_linear_exp_pars[[110]]$alpha_effective[25,]

aabbn <- aa_ind[which(aa_names == "GATA3--PGR_s1")] + 89
sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[aabbn]]$Accuracy_All, decreasing = T, index.return=T)$ix
GEMSTAT_based_on_linear_exp_pars[[aabbn]]$coop[[25]]
GEMSTAT_based_on_linear_exp_pars[[aabbn]]$alpha_effective[25,]


aa_sr1 <- integer(0)
for(i in c(c(54:112), c(130:152), c(170:193), c(210:229), c(250:274))){
  aa <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]]$Accuracy_All, decreasing = T, index.return=T)$ix[1:5]
  aa_sr1 <- c(aa_sr1, aa)
}
sort(table(aa_sr1))
########################################################################################################################
########################################################################################################################
# exp 279
# knock down all TFs together (sanity check)
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[54,])
GEMSTAT_experiment_info_df[279, 8] <- "KD of All TFs at the same time"
GEMSTAT_experiment_info_df[279, 1] <- "experiment_279"

#### reading the results
for(i in 279:279){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
########################################################################################################################
########################################################################################################################
########################################################################################################################

# Finding out why is PGR so crucial for the models?
########################################################################################################################
######### its affinity in the chosen enhancers
aa_aff <- ER_52_opt_simAnn_2_input$Affinity_scores[my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52_gene_gte3common_index]
aa_enh_ind <- matrix(nrow = 148, ncol = 39)
for(i in 1:nrow(aa_enh_ind)){
  aa_enh_ind[i, ] <- Sim_Ann_weighted_148_restart_all_models_eval[[i]]$Enhancer_index[my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52_gene_gte3common_index]
}
aa_aff_enh <- list()
for(i in 1:nrow(aa_enh_ind)){
  aa_aff_enh[[i]] <- matrix(nrow = 39, ncol = 19)
  for(j in 1:nrow(aa_aff_enh[[i]])){
    aa_aff_enh[[i]][j, ] <- aa_aff[[j]][aa_enh_ind[i, j],]
  }
}
aa_aff_enh_all <- do.call(rbind, aa_aff_enh)
colnames(aa_aff_enh_all) <- names(TF.motifs.Shrinked.count)
boxplot.matrix(aa_aff_enh_all, las = 2, outline = F)
####### Read annotations of experiment number 54

MotifWriter(motif.List = TF.motifs.Shrinked.halfsites.count, output.File.Name = "Seeded_GEMSTAT_ens/Experiment_54/NewMotif/motifs")
GEMSTAT_based_on_linear_exp_annotation_54 <- list()
for(i in 1:148){
  GEMSTAT_based_on_linear_exp_annotation_54[[i]] <- annotation_reader(annot_file = paste0("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_54/Annotation_newMotif/Experiment_54_", i, ".annot"), 
                                                                      TF_names = names(TF.motifs.Shrinked.halfsites.count))
}

aa_sr <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Accuracy_All, decreasing = T, index.return = T)$ix
aa_by_model_TF <- list()
for(aa_model in 1:148){
  aa_by_model_TF[[aa_model]] <- list()
  for(i in 1:length(TF.motifs.Shrinked.halfsites.count)){
    aa_by_model_TF[[aa_model]][[i]] <- do.call(rbind, lapply(GEMSTAT_based_on_linear_exp_annotation_54[[aa_sr[aa_model]]], "[[", i)) 
  }
  names(aa_by_model_TF[[aa_model]]) <- names(TF.motifs.Shrinked.halfsites.count)
}
names(aa_by_model_TF) <- aa_sr
aa_nu_thr_6 <-  matrix( nrow = 148, ncol = length(TF.motifs.Shrinked.halfsites.count))
for(aa_model in 1:148){
  for(i in 1:length(TF.motifs.Shrinked.halfsites.count)){
    aa_nu_thr_6[aa_model, i] <- sum(aa_by_model_TF[[aa_model]][[i]][, 5] > 0.6)
  }
}
colnames(aa_nu_thr_6) <- names(TF.motifs.Shrinked.halfsites.count)
rownames(aa_nu_thr_6) <- names(aa_by_model_TF) 

aa_newthr <- c(0.83, 0.29, 0.55, 0.35, 0.78, 0.52, 0.72, 0.30, 0.35,
               0.81, 0.48, 0.54, 0.86, 0.84, 0.56, 0.52, 0.62, 0.66, 0.43, 0.43)
aa_newthr2 <- c(0.83, 0.29, 0.55, 0.35, 0.78, 0.52, 0.27, 0.30, 0.35,
               0.81, 0.48, 0.20, 0.86, 0.84, 0.56, 0.52, 0.62, 0.66, 0.43, 0.43)
aa_newthr3 <- c(0.17, 0.71, 0.6, 0.65, 0.22,0.48, 0.6, 0.70, 0.65, 0.19, 0.52, 0.81,
                0.14, 0.16, 0.44, 0.48, 0.38, 0.34, 0.57, 0.57)

names(aa_newthr3) <- names(TF.motifs.Shrinked.halfsites.count)

aa_nu_thr_cu <-  matrix( nrow = 148, ncol = length(TF.motifs.Shrinked.halfsites.count))
for(aa_model in 1:148){
  for(i in 1:length(TF.motifs.Shrinked.halfsites.count)){
    aa_nu_thr_cu[aa_model, i] <- sum(aa_by_model_TF[[aa_model]][[i]][, 4] <= TF.motifs.Shrinked.halfsites.count_2_maxLLR[i](1- aa_newthr3[i]))
  }
}
colnames(aa_nu_thr_cu) <- names(TF.motifs.Shrinked.halfsites.count)
rownames(aa_nu_thr_cu) <- names(aa_by_model_TF) 


aa_nu_thr_cu_max <- matrix( nrow = 148, ncol = length(TF.motifs.Shrinked.halfsites.count))
for(aa_model in 1:148){
  for(i in 1:length(TF.motifs.Shrinked.halfsites.count)){
    aa_nu_thr_cu_max[aa_model, i] <- max(aa_by_model_TF[[aa_model]][[i]][, 5])
  }
}

colnames(aa_nu_thr_cu_max) <- names(TF.motifs.Shrinked.halfsites.count)
rownames(aa_nu_thr_cu_max) <- names(aa_by_model_TF) 
boxplot.matrix(aa_nu_thr_cu_max, las = 2, main = "site_count_thr_cu_all_models")

boxplot.matrix(aa_nu_thr_cu, las = 2, main = "site_count_thr_cu_all_models")
boxplot.matrix(aa_nu_thr_cu[, -c(13,14)], las = 2, main = "site_count_thr_cu_all_models")


boxplot.matrix(aa_nu_thr_6, las = 2, main = "site_count_thr6_all_models")
boxplot.matrix(aa_nu_thr_6[, -c(13)], las = 2, main = "site_count_thr6_all_models_noPGRS1")
boxplot.matrix(aa_nu_thr_6[1:15, ], las = 2, main = "site_count_thr6_top15_models")
boxplot.matrix(aa_nu_thr_6[1:15, -c(13)], las = 2, main = "site_count_thr6_top15_models_noPGRS1")


# create one matrix per Model : one row for each gene, one column for each TF, entries are the number of sites above threshold
aa_gene_tf_matlist <- list()
for(i in 1:length())
GEMSTAT_based_on_linear_exp_annotation_54
aa <- lapply(GEMSTAT_based_on_linear_exp_annotation_54[[25]], "[[", 3)

GEMSTAT_based_on_linear_exp_annotation_52 <- list()
for(i in 1:148){
  GEMSTAT_based_on_linear_exp_annotation_52[[i]] <- annotation_reader(annot_file = paste0("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_52/Annotation/Experiment_52_", i, ".annot"), 
                                                                      TF_names = names(TF.motifs.Shrinked.count))
}
aa <- lapply(GEMSTAT_based_on_linear_exp_annotation_52[[25]], "[[", 13)
### looking if a different ER motif (hocomoco improves the situation)

aa_ER_hoco <- read.table("Motifs/Hocomoco/ER/ESR1_M6217_1.02.txt", header=T, sep ="\t")
aa_ER_hoco <-  aa_ER_hoco[,2:5]
seqLogo::seqLogo(t(aa_ER_hoco))

aa_ER_hoco_cnt <- PWMtoCount(aa_ER_hoco)
aa_ER_hoco_cnt2 <- aa_ER_hoco_cnt[10:15,]


MotifWriter(list(aa_ER_hoco_cnt2),output.File.Name = "Seeded_GEMSTAT_ens/Experiment_52/ER_new_motif_half")
aa_GEMSTAT_based_on_linear_exp_annotation_52_newER <- list()
for(i in 1:148){
  aa_GEMSTAT_based_on_linear_exp_annotation_52_newER[[i]] <- annotation_reader(annot_file = paste0("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_52/Annotation_newER/Experiment_52_", i, ".annot"), 
                                                                      TF_names = "ER")
}


aa_ER_hoco_rev <- PWMEnrich::reverseComplement(t(aa_ER_hoco))
seqLogo::seqLogo(aa_ER_hoco_rev)
aa_ER_hoco_cnt_rev <- PWMtoCount(t(aa_ER_hoco_rev))
aa_ER_hoco_cnt_rev2 <- aa_ER_hoco_cnt[1:6,]
TF.motifs.Shrinked.halfsites.count_2 <- TF.motifs.Shrinked.halfsites.count
TF.motifs.Shrinked.halfsites.count_2$ESR1 <- aa_ER_hoco_cnt_rev2
save(TF.motifs.Shrinked.halfsites.count_2, file = "Seeded_GEMSTAT_ens/new_loose_motif.RData")
#########################################################################################################################
# reading annotation experiment 305
aa_GEMSTAT_based_on_linear_exp_annotation_305 <- list()
for(i in 1:148){
  aa_GEMSTAT_based_on_linear_exp_annotation_305[[i]] <- annotation_reader(annot_file = paste0("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_305/Annotation/Experiment_52_", i, ".annot"), 
                                                                               TF_names = names(TF.motifs.Shrinked.halfsites.count_2))
}
TF.motifs.Shrinked.halfsites.count_2_maxLLR <- c(6.90008, 11.0531, 6.2838, 10.705, 7.12241, 12.2036,
                                                 10.3949, 16.4493, 12.0708, 6.94399, 8.91828, 11.544,
                                                 3.76517, 6.73403, 7.49103, 7.95261, 7.64841, 7.49922,
                                                 10.0583, 10.7687)
aa <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[305]]$Accuracy_All, decreasing = T, index.return=T)$ix
aa_all <- list()
for(aa_tf in c(18)){
  aa_ER <- list()
  aa_ER_ab <- list()
  aa_er_site_nu <- matrix(0L, nrow = length(aa_GEMSTAT_based_on_linear_exp_annotation_305), ncol = 39)
  for(i in 1:length(aa_GEMSTAT_based_on_linear_exp_annotation_305)){
    aa_ER[[i]] <- lapply(aa_GEMSTAT_based_on_linear_exp_annotation_305[[i]], "[[", aa_tf)
    aa_ER_ab[[i]] <- list()
    
    for(j in 1:length(aa_ER[[i]])){
      if(nrow(aa_ER[[i]][[j]]) > 1){
        aa_ER_ab[[i]][[j]] <- aa_ER[[i]][[j]][aa_ER[[i]][[j]]$V4 <= aa_newthr3[aa_tf]*TF.motifs.Shrinked.halfsites.count_2_maxLLR[aa_tf],]
        if(nrow(aa_ER_ab[[i]][[j]]) > 1){
          aa_pos <- do.call(rbind ,strsplit(aa_ER_ab[[i]][[j]][,1], split = "\\.."))
          for(k in 1:(nrow(aa_pos) - 1)){
            if(as.integer(aa_pos[k+1, 1]) - as.integer(aa_pos[k, 2]) <= 2 &
               as.integer(aa_pos[k+1, 1]) - as.integer(aa_pos[k, 2]) > 0 &      
               (
                 (aa_ER_ab[[i]][[j]][k,   2] == "+" &
                  aa_ER_ab[[i]][[j]][k+1, 2] == "+") 
                  |
                  (aa_ER_ab[[i]][[j]][k,   2] == "-" &
                   aa_ER_ab[[i]][[j]][k+1, 2] == "-")
               ) 
               ){
              aa_er_site_nu[i, j] <- aa_er_site_nu[i, j] + 1
            }
          } # end of loop over matches
        } # if there is any match better than 0.55
      } # if there is any match
      
    } # end of loop over genes
  } # end of loop over models
  aa_all[[aa_tf]] <- aa_er_site_nu
}
names(aa_all) <- names(TF.motifs.Shrinked.halfsites.count_2)

boxplot.matrix(aa_all[[18]])
rownames(aa_er_site_nu) <- c(1:148)
colnames(aa_er_site_nu) <- names(aa_ER[[i]])
par(mfrow = c(1,1), mar= c(6,4,4,4))
boxplot.matrix(aa_er_site_nu, las = 2)
which.max(rowSums(aa_er_site_nu))

aa_sum_nz <- integer(148)
for(i in 1:148){
  aa_sum_nz[i] <- sum(aa_er_site_nu[i,] != 0)
}
hist(aa_sum_nz, main= "number of genes that have ER sites")
#########################################################################################################################
#########################################################################################################################

#########################################################################################################################
########## its concentration profile
plotExpression(expMat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer,
               exportplot = T,
               filename = "TF_Exp.png",
               .dendrogram = "none",
               .Rowv = F, .Colv = F, colorPoints = c(-1, 0, 0.33, 0.66, 1))
aa_fc <- matrix(nrow = 20, ncol = 45)
for(i in 1:45){
  aa_fc[,i] <- log2((TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer[,2*i] + 0.001)/(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer[,2*i - 1]+ 0.001))
}
colnames(aa_fc) <- colnames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3)
rownames(aa_fc) <- rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer)
plotExpression(expMat = aa_fc,
               exportplot = T,
               filename = "TF_Exp_FC.png",
               .dendrogram = "none",
               .Rowv = F, .Colv = F, colorPoints = c(-1, 0, 0.33, 0.66, 1))
aa_cond_table <- matrix(nrow = 3, ncol = 45)
rownames(aa_cond_table) <- c(1, 0 , -1)
colnames(aa_cond_table) <- colnames(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3)
for(i in 1:ncol(aa_cond_table)){
  for(j in 1:nrow(aa_cond_table)){
    aa_cond_table[j, i] <- sum(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3[, i] == as.integer(rownames(aa_cond_table)[j]), na.rm=T)/sum(!is.na(my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_39_enh_Common_gte3[, i]))
  }
}

plotExpression(expMat = rbind(aa_cond_table,aa_fc),
               exportplot = T,
               filename = "TF_Exp_FC_table.png",
               .dendrogram = "none",
               .Rowv = F, .Colv = F, colorPoints = c(-1, 0, 0.33, 0.66, 1))
########################################################################################################################
########################################################################################################################
##### Converting LLRs to Pvalues for each TF and choosing an LLR corresponding to a certain p-val for each TF

## Writing Each TF in a different file
for(i in 1:length(TF.motifs.Shrinked.halfsites.count)){
  MotifWriter(motif.List = TF.motifs.Shrinked.halfsites.count[i], 
              output.File.Name = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/MotifLLR2Pvalue/ER_monomers/",names(TF.motifs.Shrinked.halfsites.count)[i]))
}



########################################################################################################################
########################################################################################################################
# exp 280
# changing the annotation thresholds, fizing parameters for TFs that don't have binding sites with the new thresholds
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[54,])
GEMSTAT_experiment_info_df[280, 8] <- "setting annotation threshold to c(0.83, 0.29, 0.55, 0.35, 0.78,0.52,
0.72, 0.30, 0.35,0.81, 0.48, 0.54, 0.86, 0.84, 0.56, 0.52, 0.62, 0.66, 0.43, 0.43) and fixing parameters
for CEBPB, JUN_1, LEF1, NKX3-1, PBX1, TFAP2C; and changed na to 5"
GEMSTAT_experiment_info_df[280, 6] <- paste0(GEMSTAT_experiment_info_df[280, 6], "|RXRA-RXRA_DIMER_dist_2")
GEMSTAT_experiment_info_df[280, 1] <- "experiment_280"

#### reading the results
for(i in 280:280){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa3 <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(54, 280)], "[[", 7))
boxplot.matrix(aa3)
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
########################################################################################################################
########################################################################################################################
# exp 281
# same as 280 but setting na to 10
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[280,])
GEMSTAT_experiment_info_df[281, 8] <- "same as experiment 280 but changed na to 10"
GEMSTAT_experiment_info_df[281, 1] <- "experiment_281"

#### reading the results
for(i in 281:281){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa3 <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(54, 280, 281)], "[[", 7))
boxplot.matrix(aa3)
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
########################################################################################################################
########################################################################################################################
# exp 282
# same as 280 but changing annotation thresh for JUN_2 and PBX1, and setting PBX1 parameters free
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[280,])
GEMSTAT_experiment_info_df[282, 8] <- "setting annotation threshold to c(0.83, 0.29, 0.55, 0.35, 0.78,0.52, 0.27, 0.30, 0.35,0.81, 0.48, 0.20, 0.86, 0.84, 0.56, 0.52, 0.62, 0.66, 0.43, 0.43) and fixing parameters for CEBPB, JUN_1, LEF1, NKX3-1, TFAP2C; and changed na to 5"
GEMSTAT_experiment_info_df[282, 1] <- "experiment_282"

#### reading the results
for(i in 282:282){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa3 <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(54, 282, 305)], "[[", 7))
boxplot.matrix(aa3)
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)

sort(aa3[, 4], decreasing = T, index.return=T)$ix[1:10]
########################################################################################################################
# KD of each TF one by one
for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[282,])
  GEMSTAT_experiment_info_df[282+i, 8] <- paste("in silico knockdown On Experiment 282 of TFs ", names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[282+i, 1] <- paste0("experiment_", 282+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[282,])
GEMSTAT_experiment_info_df[282+21, 8] <- "in silico knockdown On Experiment 282 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[282+21, 1] <- paste0("experiment_", 282+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[282,])
GEMSTAT_experiment_info_df[282+22, 8] <- "in silico knockdown On Experiment 282 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[282+22, 1] <- paste0("experiment_", 282+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp282.R")
  aa[[i]][18] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  if(i == 21){
    aa[[i]][19] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][19] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][19] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 282 + i, ",")
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][82] <-  "                      na = 0,"
  aa[[i]][90] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F"
  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=282,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=282,")
  aa[[i]] <- c(aa[[i]], "                              coop=282, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=282,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=282,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=282, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=282,")
  aa[[i]] <- c(aa[[i]], "                              Params=282, ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=T,")
  aa[[i]] <- c(aa[[i]], "                              seq=282,")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=282,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=282, ")
  aa[[i]] <- c(aa[[i]], "                              wights=282, ")
  aa[[i]] <- c(aa[[i]], "                              jobfile=282, ")
  aa[[i]] <- c(aa[[i]], "                              .prev_na=5, ")
  aa[[i]] <- c(aa[[i]], "                              .optim=F, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 282+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  282 + i, ".R"))
}
setwd("..")

for(i in 283:304){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                                      .fold_change = T,
                                                                                                                      .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[282]]$threshold)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[283:304], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[282]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[282]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[282]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:15]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[282]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[282]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[282]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[282]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[283:304],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=15, model_index=integer(0), plot_res=T, fixed_ylim = F)
########################################################################################################################
########################################################################################################################
# Exp 305
# trying a new motif for ER: TF.motifs.Shrinked.halfsites.count_2
# TF.motifs.Shrinked.halfsites.count_2
# annotation thresholds: c(0.83, 0.29, 0.05, 0.35, 0.78,0.52, 0.27, 0.30, 0.35, 0.81, 0.48, 0.20, 0.86, 0.84, 0.56, 0.52, 0.62, 0.66, 0.43, 0.43)
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[282,])
GEMSTAT_experiment_info_df[305, 8] <- "setting annotation threshold to c(0.83, 0.29, 0.05, 0.35, 0.78,0.52, 0.27, 0.30, 0.35, 0.81, 0.48, 0.20, 0.86, 0.84, 0.56, 0.52, 0.62, 0.66, 0.43, 0.43) and fixing parameters for CEBPB, JUN_1, LEF1, NKX3-1, TFAP2C; and changed na to 5; trying a new motif for ER: trying a new (loose) motif for ER: TF.motifs.Shrinked.halfsites.count_2"
GEMSTAT_experiment_info_df[305, 1] <- paste0("experiment_", 305)

for(i in 305:305){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
par(mfrow = c(1, 1), mar = c(7,4,2,2))
aa <- cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Accuracy_All,
            GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[305]]$Accuracy_All)
colnames(aa) <- c("original", "loose_motif")
boxplot.matrix(aa,
               las=2)
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
########################################################################################################################
########################################################################################################################
# KD of each TF one by one
for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[305,])
  GEMSTAT_experiment_info_df[305+i, 8] <- paste("in silico knockdown On Experiment 305 of TFs ", names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[305+i, 1] <- paste0("experiment_", 305+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[305,])
GEMSTAT_experiment_info_df[305+21, 8] <- "in silico knockdown On Experiment 305 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[305+21, 1] <- paste0("experiment_", 305+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[305,])
GEMSTAT_experiment_info_df[305+22, 8] <- "in silico knockdown On Experiment 305 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[305+22, 1] <- paste0("experiment_", 305+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp305.R")
  aa[[i]][18] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  if(i == 21){
    aa[[i]][19] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][19] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][19] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 305 + i, ",")
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][82] <-  "                      na = 0,"
  aa[[i]][90] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F"
  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=305,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=305,")
  aa[[i]] <- c(aa[[i]], "                              coop=305, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=305,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=305,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=305, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=305,")
  aa[[i]] <- c(aa[[i]], "                              Params=305, ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=T,")
  aa[[i]] <- c(aa[[i]], "                              seq=305,")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=305,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=305, ")
  aa[[i]] <- c(aa[[i]], "                              wights=305, ")
  aa[[i]] <- c(aa[[i]], "                              jobfile=305, ")
  aa[[i]] <- c(aa[[i]], "                              .prev_na=5, ")
  aa[[i]] <- c(aa[[i]], "                              .optim=F, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 305+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  305 + i, ".R"))
}
setwd("..")

for(i in 306:327){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[305]]$threshold)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[306:327], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[305]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[305]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[305]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:15]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[305]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[305]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[305]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[305]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[306:327],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=15, model_index=integer(0), plot_res=T, fixed_ylim = F)
########################################################################################################################
########################################################################################################################
# experiment 328
# same as experiment 305 but assigning zero weights to all datapoints with 0 and -1 labels.
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[305,])
GEMSTAT_experiment_info_df[328, 1] <- "experiment_328" 
GEMSTAT_experiment_info_df[328, 8] <- "same as experiment 305 but assigning zero weights to all datapoints with 0 and -1 labels. " 

for(i in 328:328){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .ignore_lab = c(0, -1))
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}


aa <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[328]]$Accuracy_All, decreasing = T, index.return=T)$ix
GEMSTAT_based_on_linear_exp_pars[[328]]$alpha_effective[25,]
hist(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[328]]$Accuracy_All, main="performance on 1s only")
########################################################################################################################
########################################################################################################################
# experiment 329
# same as experiment 305 but assigning zero weights to all datapoints with -1 label.
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[305,])
GEMSTAT_experiment_info_df[329, 1] <- "experiment_329" 
GEMSTAT_experiment_info_df[329, 8] <- "same as experiment 305 but assigning zero weights to all datapoints with -1 label." 
for(i in 329:329){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .ignore_lab = c(-1))
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
hist(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[329]]$Accuracy_All, main="performance on 1s and 0s")
aa <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[329]]$Accuracy_All, decreasing = T, index.return=T)$ix
GEMSTAT_based_on_linear_exp_pars[[329]]$binding[31,]


########################################################################################################################
########################################################################################################################
# KD of experiment 328
for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[328,])
  GEMSTAT_experiment_info_df[329+i, 8] <- paste("in silico knockdown On Experiment 328 of TFs ", names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[329+i, 1] <- paste0("experiment_", 329+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[328,])
GEMSTAT_experiment_info_df[329+21, 8] <- "in silico knockdown On Experiment 328 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[329+21, 1] <- paste0("experiment_", 329+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[328,])
GEMSTAT_experiment_info_df[329+22, 8] <- "in silico knockdown On Experiment 328 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[329+22, 1] <- paste0("experiment_", 329+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp328.R")
  aa[[i]][18] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  if(i == 21){
    aa[[i]][19] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][19] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][19] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 329 + i, ",")
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][82] <-  "                      na = 0,"
  aa[[i]][91] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F"
  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=328,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=328,")
  aa[[i]] <- c(aa[[i]], "                              coop=328, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=328,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=328,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=328, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=328,")
  aa[[i]] <- c(aa[[i]], "                              Params=328, ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=T,")
  aa[[i]] <- c(aa[[i]], "                              seq=328,")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=328,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=328, ")
  aa[[i]] <- c(aa[[i]], "                              wights=328, ")
  aa[[i]] <- c(aa[[i]], "                              jobfile=328, ")
  aa[[i]] <- c(aa[[i]], "                              .prev_na=5, ")
  aa[[i]] <- c(aa[[i]], "                              .optim=F, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 329+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  329 + i, ".R"))
}
setwd("..")

# reading the results
for(i in 330:351){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[328]]$threshold,
                                                                                                             .ignore_lab = c(-1, 0))
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[330:351], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[328]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[328]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[328]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:15]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[328]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[328]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.7, 1), ylim = c(0.7, 1), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[328]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[328]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[330:351],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=15, model_index=integer(0), plot_res=T, fixed_ylim = F)
########################################################################################################################
########################################################################################################################
# KD of experiment 329
for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[329,])
  GEMSTAT_experiment_info_df[351+i, 8] <- paste("in silico knockdown On Experiment 329 of TFs ", names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[351+i, 1] <- paste0("experiment_", 351+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[329,])
GEMSTAT_experiment_info_df[351+21, 8] <- "in silico knockdown On Experiment 329 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[351+21, 1] <- paste0("experiment_", 329+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[329,])
GEMSTAT_experiment_info_df[351+22, 8] <- "in silico knockdown On Experiment 329 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[351+22, 1] <- paste0("experiment_", 351+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp329.R")
  aa[[i]][18] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  if(i == 21){
    aa[[i]][19] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][19] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][19] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 351 + i, ",")
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][82] <-  "                      na = 0,"
  aa[[i]][91] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F"
  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=329,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=329,")
  aa[[i]] <- c(aa[[i]], "                              coop=329, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=329,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=329,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=329, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=329,")
  aa[[i]] <- c(aa[[i]], "                              Params=329, ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=T,")
  aa[[i]] <- c(aa[[i]], "                              seq=329,")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=329,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=329, ")
  aa[[i]] <- c(aa[[i]], "                              wights=329, ")
  aa[[i]] <- c(aa[[i]], "                              jobfile=329, ")
  aa[[i]] <- c(aa[[i]], "                              .prev_na=5, ")
  aa[[i]] <- c(aa[[i]], "                              .optim=F, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 351+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  351 + i, ".R"))
}
setwd("..")

for(i in 352:373){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[329]]$threshold,
                                                                                                             .ignore_lab = c(-1))
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[352:373], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")
par(mfrow=c(1, 1), mar = c(5,4, 4, 1))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[329]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[329]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[329]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:15]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[329]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[329]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.9), ylim = c(0.4, 0.9), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[329]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[329]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[352:373],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=15, model_index=integer(0), plot_res=T, fixed_ylim = F)

########################################################################################################################
########################################################################################################################
# experiment 374
# same as experiment 305 but assigning zero weights to all datapoints with 0 label.
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[305,])
GEMSTAT_experiment_info_df[374, 1] <- "experiment_374" 
GEMSTAT_experiment_info_df[374, 8] <- "same as experiment 305 but assigning zero weights to all datapoints with 0 label." 
for(i in 374:374){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .ignore_lab = c(0))
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
par(mfrow = c(1, 1), mar = c(4,4,4,4))
hist(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[374]]$Accuracy_All, main="performance_1s_and_-1s")
########################################################################################################################
########################################################################################################################
# KD of exp 374
for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[374,])
  GEMSTAT_experiment_info_df[374+i, 8] <- paste("in silico knockdown On Experiment 374 of TFs ", names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[374+i, 1] <- paste0("experiment_", 374+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[374,])
GEMSTAT_experiment_info_df[374+21, 8] <- "in silico knockdown On Experiment 374 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[374+21, 1] <- paste0("experiment_", 374+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[374,])
GEMSTAT_experiment_info_df[374+22, 8] <- "in silico knockdown On Experiment 374 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[374+22, 1] <- paste0("experiment_", 374+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp374.R")
  aa[[i]][18] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  if(i == 21){
    aa[[i]][19] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][19] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][19] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 374 + i, ",")
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][82] <-  "                      na = 0,"
  aa[[i]][91] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F"
  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=374,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=374,")
  aa[[i]] <- c(aa[[i]], "                              coop=374, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=374,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=374,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=374, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=374,")
  aa[[i]] <- c(aa[[i]], "                              Params=374, ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=T,")
  aa[[i]] <- c(aa[[i]], "                              seq=374,")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=374,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=374, ")
  aa[[i]] <- c(aa[[i]], "                              wights=374, ")
  aa[[i]] <- c(aa[[i]], "                              jobfile=374, ")
  aa[[i]] <- c(aa[[i]], "                              .prev_na=5, ")
  aa[[i]] <- c(aa[[i]], "                              .optim=F, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 374+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  374 + i, ".R"))
}
setwd("..")

for(i in 375:396){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[374]]$threshold,
                                                                                                             .ignore_lab = c(0))
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[375:396], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")
par(mfrow=c(1, 1), mar = c(5,4, 4, 1))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[374]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[374]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[374]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:15]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[374]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[374]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.9), ylim = c(0.4, 0.9), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[374]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[374]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[375:396],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=15, model_index=integer(0), plot_res=T, fixed_ylim = F)

########################################################################################################################
########################################################################################################################
# exp 397
# just realized that the annotation threshold is (LLR_site - LLR_max)/LLR_max  --> which is 1 - (what I though it is before).

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[305,])
GEMSTAT_experiment_info_df[397, 1] <- "experiment_397"
GEMSTAT_experiment_info_df[397, 8] <- "setting annotation threshold to c(0.17, 0.71, 0.6, 0.65, 0.22,0.48, 0.6, 0.70, 0.65, 0.19, 0.52, 0.81, 0.14, 0.16, 0.44, 0.48, 0.38, 0.34, 0.57, 0.57). Not fixing any of TF parameters. trying a new motif for ER: trying a new (loose) motif for ER: TF.motifs.Shrinked.halfsites.count_2"
c(0.17, 0.71, 0.6, 0.65, 0.22,0.48, 0.6, 0.70, 0.65, 0.19, 0.52, 0.81, 0.14, 0.16, 0.44, 0.48, 0.38, 0.34, 0.57, 0.57)

for(i in 397:397){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
########################################################################################################################
########################################################################################################################
# KD of Experiment 397
for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[397,])
  GEMSTAT_experiment_info_df[397+i, 8] <- paste("in silico knockdown On Experiment 397 of TFs ", names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[397+i, 1] <- paste0("experiment_", 397+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[397,])
GEMSTAT_experiment_info_df[397+21, 8] <- "in silico knockdown On Experiment 397 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[397+21, 1] <- paste0("experiment_", 397+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[397,])
GEMSTAT_experiment_info_df[397+22, 8] <- "in silico knockdown On Experiment 397 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[397+22, 1] <- paste0("experiment_", 397+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp397.R")
  aa[[i]][18] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  if(i == 21){
    aa[[i]][19] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][19] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][19] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 397 + i, ",")
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][82] <-  "                      na = 0,"
  aa[[i]][90] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F"
  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=397,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=397,")
  aa[[i]] <- c(aa[[i]], "                              coop=397, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=397,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=397,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=397, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=397,")
  aa[[i]] <- c(aa[[i]], "                              Params=397, ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=T,")
  aa[[i]] <- c(aa[[i]], "                              seq=397,")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=397,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=397, ")
  aa[[i]] <- c(aa[[i]], "                              wights=397, ")
  aa[[i]] <- c(aa[[i]], "                              jobfile=397, ")
  aa[[i]] <- c(aa[[i]], "                              .prev_na=5, ")
  aa[[i]] <- c(aa[[i]], "                              .optim=F, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 397+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  397 + i, ".R"))
}
setwd("..")

for(i in 398:419){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[397]]$threshold)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[398:419], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[397]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[397]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[397]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:15]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[397]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[397]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[397]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[397]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[398:419],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=15, model_index=integer(0), plot_res=T, fixed_ylim = F)
########################################################################################################################
########################################################################################################################
# Experiment 420 
# same as experiment 397, only changing the ub on cooperativity to 25000
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[397,])
GEMSTAT_experiment_info_df[420, 8] <- "same as experiment 397, only changing the ub on cooperativity to 25000"


for(i in 420:420){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(54,397, 420)], "[[", 7))
par(mfrow= c(1,1), mar = c(5,4,4,4))
boxplot.matrix(aa)

########################################################################################################################
########################################################################################################################

# creating all combinations of concentration modifications.
aa_modif_all_options <- list(list(1, -1) ,
                          list(10, -10),
                          list(c(13, 14), c(-13, -14)),
                          list(15,-15),
                          list(16,-16),
                          list(18,-18))
aa_my_comb <- rbind(rep(0, 6), GenerateMotBinary(c(1:6), 6))
aa_my_comb <- aa_my_comb + 1 # in this matrix one means concentration will be the same as ER and 2 means it will be opposite to ER
aa_modif_all_comb <- list()
for(i in 1:nrow(aa_my_comb)){
  aa_modif_all_comb[[i]] <- integer(0)
  for(j in 1:ncol(aa_my_comb)){
    aa_modif_all_comb[[i]] <- c(aa_modif_all_comb[[i]], aa_modif_all_options[[j]][aa_my_comb[i, j]])
  }
  aa_modif_all_comb[[i]] <- unlist(aa_modif_all_comb[[i]])
}

# choosing the models to work on

aa_sr1 <- integer(0)
for(i in c(c(54, 397), c(78:278))){
  if(length(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]]) > 0){
    aa <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]]$Accuracy_All, decreasing = T, index.return=T)$ix[1:5]
    aa_sr1 <- c(aa_sr1, aa)
  }
}
sort(table(aa_sr1), decreasing = T)[1:10]

# Following is the index of the chosen models
aa_model_index <- c(2,6,13,15,25,25,33,44,56,67,77,121,141)

# updating the data file
aa_name_conv <- names(TF.motifs.Shrinked.halfsites.count_2[abs(sort(unlist(aa_modif_all_options)))])
aa_name_conv[1:7] <- paste(aa_name_conv[1:7],"n",sep = "_")
aa_name_conv[8:14] <- paste(aa_name_conv[8:14],"p",sep = "_")
aa_name_conv <- cbind(sort(unlist(aa_modif_all_options)), aa_name_conv)
for(i in 1:length(aa_modif_all_comb)){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[420,])
  GEMSTAT_experiment_info_df[420+i, 1] <- paste0("experiment_", 420+i)
  GEMSTAT_experiment_info_df[420+i, 8] <- paste("Concentration Modification. Base experiment 420.  my_model_index <- c(2,6,13,15,25,25,33,44,56,67,77,121,141). my_qBTMrange <- t(matrix(c(1e-5, 1))) .", "conc modif: ",
                                                paste(aa_name_conv[aa_name_conv[, 1] %in% as.character(aa_modif_all_comb[[i]]), 2], collapse = " "))
  
}



# Creating the R script to create the inputs and job files.
setwd("Seeded_GEMSTAT_ens/")
aa <- list()

for(i in 1:length(aa_modif_all_comb)){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp420.R")
  aa[[i]][17] <- "my_qBTMrange <- t(matrix(c(1e-5, 1)))"
  aa[[i]][18] <- "my_model_index <- c(2,6,13,15,25,25,33,44,56,67,77,121,141)"
  aa[[i]][19] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  for(j in 1:length(aa_modif_all_comb[[i]])){
    if(aa_modif_all_comb[[i]][j] > 0){
      aa[[i]][20 + j-1] <- paste0('aa_myTFexpmat[',aa_modif_all_comb[[i]][j],', ] <- rep(c(0,1), 45)')
    }else{
      aa[[i]][20 + j-1] <- paste0('aa_myTFexpmat[',-aa_modif_all_comb[[i]][j],', ] <- rep(c(1,0), 45)')
    }
  }
  
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 420 + i, ",")
  aa[[i]][42] <- "                      model_evaluations = Sim_Ann_weighted_148_restart_all_models_eval[my_model_index],"
  aa[[i]][45] <- "                      model_parameters = Sim_Ann_weighted_148_restart_parameters[my_model_index, ],"
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][68] <- "                      qBTMrange=my_qBTMrange,"
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  420 + i, ".R"))
}
setwd("..")

# reading the results
for(i in 421:484){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[421:484], "[[", 7))
boxplot.matrix(aa, las= 2)
abline(h=seq(0, 1, 0.05), col = 2, lty = 4)

which.max(aa[, 28])

aa2 <- matrix(nrow = 148, ncol = (ncol(aa) + 2))
aa2[, 1] <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[420], "[[", 7))
for(i in 1:ncol(aa)){
  aa2[1:13, i+1] <- aa[, i] 
}
aa2[, ncol(aa2)] <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[485], "[[", 7))
boxplot.matrix(cbind(do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(54, 397)], "[[", 7)), aa2))
abline(h=seq(0, 1, 0.05), col = 2, lty = 4)


########################################################################################################################
########################################################################################################################
# Experiment 485 
# same as experiment 420, only changing the ub on qBTM to 5
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[420,])
GEMSTAT_experiment_info_df[485, 8] <- paste("same as experiment 420, only changing the ub on qBTM to 5")
GEMSTAT_experiment_info_df[485, 1] <- paste0("experiment_", 485)

for(i in 485:485){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(54,397, 420, 485)], "[[", 7))
par(mfrow= c(1,1), mar = c(5,4,4,4))
boxplot.matrix(aa)
########################################################################################################################
########################################################################################################################
# Experiment 486
# same as experiment 397, only changing the ub on binding to 25  and upper bind on coop to 2500

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[397,])
GEMSTAT_experiment_info_df[486, 8] <- paste("same as experiment 397, only changing the ub on binding to 25  and upper bind on coop to 2500")
GEMSTAT_experiment_info_df[486, 1] <- paste0("experiment_", 486)

for(i in 486:486){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(54,397, 420, 485, 486)], "[[", 7))
boxplot.matrix(aa)
abline(h=seq(0, 1, 0.05), col = 2, lty = 4)

########################################################################################################################
########################################################################################################################
# Experiment 487
# same as experiment 397, only changing the ub on qBTM to 1

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[397,])
GEMSTAT_experiment_info_df[487, 8] <- paste("same as experiment 397, only changing the ub on qBTM to 1")
GEMSTAT_experiment_info_df[487, 1] <- paste0("experiment_", 487)

for(i in 487:487){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(54,397, 420, 485, 486, 487)], "[[", 7))
boxplot.matrix(aa)
abline(h=seq(0, 1, 0.05), col = 2, lty = 4)

########################################################################################################################
########################################################################################################################
# Experiment 488
# same as experiment 54, only adding RXRA coop and na is 5 instead of 10 for 13 good models and using the new ER motif
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[54,])
GEMSTAT_experiment_info_df[488, 8] <- paste("same as experiment 54, only adding RXRA coop and na is 5 instead of 10 for 13 good models (model list the same as Experiment 421) and using the new ER motif")
GEMSTAT_experiment_info_df[488, 1] <- paste0("experiment_", 488)

for(i in 488:488){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
par(mfrow=c(1, 1), mar = c(4,4,4,4))
hist(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[54]]$Accuracy_All, breaks = 13)

########################################################################################################################
########################################################################################################################
# Experiment 489
# same as experiment 397, only using the new ER motif and only for top 13 models

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[397,])
GEMSTAT_experiment_info_df[489, 8] <- paste("same as experiment 397, only using the new ER motif and only for top 13 models (model list the same as Experiment 421)")
GEMSTAT_experiment_info_df[489, 1] <- paste0("experiment_", 489)

for(i in 489:489){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
hist(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[489]]$Accuracy_All, breaks = 13)
########################################################################################################################
########################################################################################################################
# Experiment 490
# same as experiment 397, only using the new ER motif and only for top 13 models and changing coop limit to 2500
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[397,])
GEMSTAT_experiment_info_df[490, 8] <- paste("same as experiment 397, only using the new ER motif and only for top 13 models and changing coop limit to 2500 (model list the same as Experiment 421)")
GEMSTAT_experiment_info_df[490, 1] <- paste0("experiment_", 490)

for(i in 490:490){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
hist(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[490]]$Accuracy_All, breaks = 13)
########################################################################################################################
########################################################################################################################

# doing conc modifications on exp 489
# creating all combinations of concentration modifications/ this time with the new ER motif and higher 
aa_modif_all_options <- list(list(1, -1) ,
                             list(10, -10),
                             list(c(13, 14), c(-13, -14)),
                             list(15,-15),
                             list(16,-16),
                             list(18,-18))
aa_my_comb <- rbind(rep(0, 6), GenerateMotBinary(c(1:6), 6))
aa_my_comb <- aa_my_comb + 1 # in this matrix one means concentration will be the same as ER and 2 means it will be opposite to ER
aa_modif_all_comb <- list()
for(i in 1:nrow(aa_my_comb)){
  aa_modif_all_comb[[i]] <- integer(0)
  for(j in 1:ncol(aa_my_comb)){
    aa_modif_all_comb[[i]] <- c(aa_modif_all_comb[[i]], aa_modif_all_options[[j]][aa_my_comb[i, j]])
  }
  aa_modif_all_comb[[i]] <- unlist(aa_modif_all_comb[[i]])
}

# Following is the index of the chosen models
aa_model_index <- c(2,6,13,15,25,25,33,44,56,67,77,121,141)

# updating the data file
aa_name_conv <- names(TF.motifs.Shrinked.halfsites.count_2[abs(sort(unlist(aa_modif_all_options)))])
aa_name_conv[1:7] <- paste(aa_name_conv[1:7],"n",sep = "_")
aa_name_conv[8:14] <- paste(aa_name_conv[8:14],"p",sep = "_")
aa_name_conv <- cbind(sort(unlist(aa_modif_all_options)), aa_name_conv)
for(i in 1:length(aa_modif_all_comb)){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[489,])
  GEMSTAT_experiment_info_df[490+i, 1] <- paste0("experiment_", 490+i)
  GEMSTAT_experiment_info_df[490+i, 8] <- paste("Concentration Modification. Base experiment 489.  my_model_index <- c(2,6,13,15,25,25,33,44,56,67,77,121,141).", "conc modif: ",
                                                paste(aa_name_conv[aa_name_conv[, 1] %in% as.character(aa_modif_all_comb[[i]]), 2], collapse = " "))
  
}

# Creating the R script to create the inputs and job files.
setwd("Seeded_GEMSTAT_ens/")
aa <- list()

for(i in 1:length(aa_modif_all_comb)){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp489.R")
 # aa[[i]][17] <- "my_qBTMrange <- t(matrix(c(1e-5, 1)))"
 # aa[[i]][18] <- "my_model_index <- c(2,6,13,15,25,25,33,44,56,67,77,121,141)"
  aa[[i]][19] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  for(j in 1:length(aa_modif_all_comb[[i]])){
    if(aa_modif_all_comb[[i]][j] > 0){
      aa[[i]][20 + j-1] <- paste0('aa_myTFexpmat[',aa_modif_all_comb[[i]][j],', ] <- rep(c(0,1), 45)')
    }else{
      aa[[i]][20 + j-1] <- paste0('aa_myTFexpmat[',-aa_modif_all_comb[[i]][j],', ] <- rep(c(1,0), 45)')
    }
  }
  
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 490 + i, ",")
#  aa[[i]][42] <- "                      model_evaluations = Sim_Ann_weighted_148_restart_all_models_eval[my_model_index],"
 # aa[[i]][45] <- "                      model_parameters = Sim_Ann_weighted_148_restart_parameters[my_model_index, ],"
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
#  aa[[i]][68] <- "                      qBTMrange=my_qBTMrange,"
  aa[[i]][90] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F"
   aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
   aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=489,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=489,")
  aa[[i]] <- c(aa[[i]], "                              coop=489, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=489,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=489,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=489, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=489,")
  aa[[i]] <- c(aa[[i]], "                              Params=489, ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=F,")
  aa[[i]] <- c(aa[[i]], "                              seq=489,")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=489,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=489, ")
  aa[[i]] <- c(aa[[i]], "                              wights=489, ")
  aa[[i]] <- c(aa[[i]], "                              jobfile=489, ")
  aa[[i]] <- c(aa[[i]], "                              .prev_na=5, ")
  aa[[i]] <- c(aa[[i]], "                              .optim=T, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 490+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  # 
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  490 + i, ".R"))
}
setwd("..")

# reading the results
for(i in 491:554){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
par(mfrow = c(1, 1), mar = c(6,4,4,4))
aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[491:554], "[[", 7))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[489]]$Accuracy_All, aa), las= 2)
abline(h=seq(0, 1, 0.01), col = 2, lty = 4)

# get the max of each column:
aa_max <- apply(X = aa, MARGIN = 2, max)
aa_max_ind <- sort(aa_max, decreasing = T, index.return=T)$ix
aa_comb <- do.call(rbind, aa_modif_all_comb)
aa_comb[aa_max_ind, ]
aa_comb_bin <- matrix(as.numeric(aa_comb[aa_max_ind, ] > 0), nrow = 64)
aa_comb_bin <- aa_comb_bin[, -3]
colnames(aa_comb_bin) <- c("AR", "GR", "PGR", "RARA", "RARG", "RXRA")

heatmap.2(aa_comb_bin, Rowv = F, Colv = F, dendrogram = "none")


aa2 <- matrix(nrow = 148, ncol = (ncol(aa) + 2))
aa2[, 1] <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[420], "[[", 7))
for(i in 1:ncol(aa)){
  aa2[1:13, i+1] <- aa[, i] 
}
aa2[, ncol(aa2)] <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[485], "[[", 7))
boxplot.matrix(cbind(do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(54, 397)], "[[", 7)), aa2))
abline(h=seq(0, 1, 0.05), col = 2, lty = 4)


########################################################################################################################
########################################################################################################################
# doing KD for experiment 521


for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[521,])
  GEMSTAT_experiment_info_df[554+i, 8] <- paste("in silico knockdown On Experiment 521 of TFs ", 
                                                names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[554+i, 1] <- paste0("experiment_", 554+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[521,])
GEMSTAT_experiment_info_df[554+21, 8] <- "in silico knockdown On Experiment 521 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[554+21, 1] <- paste0("experiment_", 554+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[521,])
GEMSTAT_experiment_info_df[554+22, 8] <- "in silico knockdown On Experiment 521 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[554+22, 1] <- paste0("experiment_", 554+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp521.R")
  #aa[[i]][18] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  if(i == 21){
    aa[[i]][27] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][27] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][27] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 554 + i, ",")
 # aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][82] <-  "                      na = 0,"
#  aa[[i]][90] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F"
#  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
#  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]][93] <- "GEMSTAT_input_copier_multienh(lower_bounds=521,"
  aa[[i]][94] <-  "                              upper_bounds=521,"
  aa[[i]][95] <-   "                              coop=521, "
  aa[[i]][96] <-   "                              Enh_gene_map=521,"
  aa[[i]][97] <-   "                              free_fix=521,"
  aa[[i]][98] <-   "                              Gene_Exp=521, "
  aa[[i]][99] <-   "                              Motifs=521,"
  aa[[i]][100] <-   "                              Params=521, "
  aa[[i]][101] <-   "                              Params_from_output=T,"
  aa[[i]][102] <-   "                              seq=521,"
  aa[[i]][103] <-   "                              TFexp=numeric(0),"
  aa[[i]][104] <-   "                              TF_info=521,"
  aa[[i]][105] <-   "                              Treat_cont_map=521, "
  aa[[i]][106] <-   "                              wights=521, "
  aa[[i]][107] <-   "                              jobfile=521, "
  aa[[i]][108] <-   "                              .prev_na=5, "
  aa[[i]][109] <-   "                              .optim=F, "
  aa[[i]][110] <-   paste0("                              cur_exp=", 554+i, ",")
  aa[[i]][111] <-   "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")"
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  554 + i, ".R"))
}
setwd("..")

for(i in 555:576){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[521]]$threshold)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[555:576], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[521]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[521]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[521]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:13]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[521]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[521]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[521]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[397]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[398:419],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=15, model_index=integer(0), plot_res=T, fixed_ylim = F)
########################################################################################################################
########################################################################################################################
# doing KD for experiment 551
for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[551,])
  GEMSTAT_experiment_info_df[576+i, 8] <- paste("in silico knockdown On Experiment 551 of TFs ", 
                                                names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[576+i, 1] <- paste0("experiment_", 576+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[551,])
GEMSTAT_experiment_info_df[576+21, 8] <- "in silico knockdown On Experiment 551 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[576+21, 1] <- paste0("experiment_", 576+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[551,])
GEMSTAT_experiment_info_df[576+22, 8] <- "in silico knockdown On Experiment 551 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[576+22, 1] <- paste0("experiment_", 576+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp551.R")
  #aa[[i]][18] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  if(i == 21){
    aa[[i]][27] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][27] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][27] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 576 + i, ",")
  # aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][82] <-  "                      na = 0,"
  #  aa[[i]][90] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F"
  #  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
  #  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]][93] <- "GEMSTAT_input_copier_multienh(lower_bounds=551,"
  aa[[i]][94] <-  "                              upper_bounds=551,"
  aa[[i]][95] <-   "                              coop=551, "
  aa[[i]][96] <-   "                              Enh_gene_map=551,"
  aa[[i]][97] <-   "                              free_fix=551,"
  aa[[i]][98] <-   "                              Gene_Exp=551, "
  aa[[i]][99] <-   "                              Motifs=551,"
  aa[[i]][100] <-   "                              Params=551, "
  aa[[i]][101] <-   "                              Params_from_output=T,"
  aa[[i]][102] <-   "                              seq=551,"
  aa[[i]][103] <-   "                              TFexp=numeric(0),"
  aa[[i]][104] <-   "                              TF_info=551,"
  aa[[i]][105] <-   "                              Treat_cont_map=551, "
  aa[[i]][106] <-   "                              wights=551, "
  aa[[i]][107] <-   "                              jobfile=551, "
  aa[[i]][108] <-   "                              .prev_na=5, "
  aa[[i]][109] <-   "                              .optim=F, "
  aa[[i]][110] <-   paste0("                              cur_exp=", 576+i, ",")
  aa[[i]][111] <-   "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")"
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  576 + i, ".R"))
}
setwd("..")

for(i in 577:598){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[551]]$threshold)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[577:598], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[551]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[551]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[551]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:13]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[551]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[551]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[521]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[551]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[577:598],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=13, model_index=integer(0), plot_res=T, fixed_ylim = F)

########################################################################################################################
########################################################################################################################
# doing KD for experiment 519

for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[519,])
  GEMSTAT_experiment_info_df[598+i, 8] <- paste("in silico knockdown On Experiment 519 of TFs ", 
                                                names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[598+i, 1] <- paste0("experiment_", 598+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[519,])
GEMSTAT_experiment_info_df[598+21, 8] <- "in silico knockdown On Experiment 519 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[598+21, 1] <- paste0("experiment_", 598+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[519,])
GEMSTAT_experiment_info_df[598+22, 8] <- "in silico knockdown On Experiment 519 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[598+22, 1] <- paste0("experiment_", 598+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp519.R")
  #aa[[i]][18] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  if(i == 21){
    aa[[i]][27] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][27] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][27] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 598 + i, ",")
  # aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][82] <-  "                      na = 0,"
  #  aa[[i]][90] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F"
  #  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
  #  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]][93] <- "GEMSTAT_input_copier_multienh(lower_bounds=519,"
  aa[[i]][94] <-  "                              upper_bounds=519,"
  aa[[i]][95] <-   "                              coop=519, "
  aa[[i]][96] <-   "                              Enh_gene_map=519,"
  aa[[i]][97] <-   "                              free_fix=519,"
  aa[[i]][98] <-   "                              Gene_Exp=519, "
  aa[[i]][99] <-   "                              Motifs=519,"
  aa[[i]][100] <-   "                              Params=519, "
  aa[[i]][101] <-   "                              Params_from_output=T,"
  aa[[i]][102] <-   "                              seq=519,"
  aa[[i]][103] <-   "                              TFexp=numeric(0),"
  aa[[i]][104] <-   "                              TF_info=519,"
  aa[[i]][105] <-   "                              Treat_cont_map=519, "
  aa[[i]][106] <-   "                              wights=519, "
  aa[[i]][107] <-   "                              jobfile=519, "
  aa[[i]][108] <-   "                              .prev_na=5, "
  aa[[i]][109] <-   "                              .optim=F, "
  aa[[i]][110] <-   paste0("                              cur_exp=", 598+i, ",")
  aa[[i]][111] <-   "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")"
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  598 + i, ".R"))
}
setwd("..")

for(i in 599:620){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[519]]$threshold)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[599:620], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

par(mfrow = c(1, 1), mar = c(4,4,4,4))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[519]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[519]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[519]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:13]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[519]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[519]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[519]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[519]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[599:620],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=13, model_index=integer(0), plot_res=T, fixed_ylim = F)


########################################################################################################################
########################################################################################################################
# doing KD for experiment 501
for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[501,])
  GEMSTAT_experiment_info_df[620+i, 8] <- paste("in silico knockdown On Experiment 501 of TFs ", 
                                                names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[620+i, 1] <- paste0("experiment_", 620+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[501,])
GEMSTAT_experiment_info_df[620+21, 8] <- "in silico knockdown On Experiment 501 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[620+21, 1] <- paste0("experiment_", 620+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[501,])
GEMSTAT_experiment_info_df[620+22, 8] <- "in silico knockdown On Experiment 501 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[620+22, 1] <- paste0("experiment_", 620+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp501.R")
  if(i == 21){
    aa[[i]][27] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][27] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][27] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 620 + i, ",")
  aa[[i]][82] <-  "                      na = 0,"
  aa[[i]][93] <- "GEMSTAT_input_copier_multienh(lower_bounds=501,"
  aa[[i]][94] <-  "                              upper_bounds=501,"
  aa[[i]][95] <-   "                              coop=501, "
  aa[[i]][96] <-   "                              Enh_gene_map=501,"
  aa[[i]][97] <-   "                              free_fix=501,"
  aa[[i]][98] <-   "                              Gene_Exp=501, "
  aa[[i]][99] <-   "                              Motifs=501,"
  aa[[i]][100] <-   "                              Params=501, "
  aa[[i]][101] <-   "                              Params_from_output=T,"
  aa[[i]][102] <-   "                              seq=501,"
  aa[[i]][103] <-   "                              TFexp=numeric(0),"
  aa[[i]][104] <-   "                              TF_info=501,"
  aa[[i]][105] <-   "                              Treat_cont_map=501, "
  aa[[i]][106] <-   "                              wights=501, "
  aa[[i]][107] <-   "                              jobfile=501, "
  aa[[i]][108] <-   "                              .prev_na=5, "
  aa[[i]][109] <-   "                              .optim=F, "
  aa[[i]][110] <-   paste0("                              cur_exp=", 620+i, ",")
  aa[[i]][111] <-   "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")"
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  620 + i, ".R"))
}
setwd("..")
bash_batch_run_zip_transfer_process(exp_index = c(621:642))

for(i in 621:642){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[501]]$threshold)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[621:642], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

par(mfrow = c(1, 1), mar = c(4,4,4,4))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[501]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[501]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[501]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:13]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[501]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[501]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[501]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[501]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[621:642],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=13, model_index=integer(0), plot_res=T, fixed_ylim = F)

########################################################################################################################
########################################################################################################################
# doing KD for experiment 516

for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[516,])
  GEMSTAT_experiment_info_df[642+i, 8] <- paste("in silico knockdown On Experiment 516 of TFs ", 
                                                names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[642+i, 1] <- paste0("experiment_", 642+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[516,])
GEMSTAT_experiment_info_df[642+21, 8] <- "in silico knockdown On Experiment 516 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[642+21, 1] <- paste0("experiment_", 642+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[516,])
GEMSTAT_experiment_info_df[642+22, 8] <- "in silico knockdown On Experiment 516 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[642+22, 1] <- paste0("experiment_", 642+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp516.R")
  if(i == 21){
    aa[[i]][27] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][27] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][27] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 642 + i, ",")
  aa[[i]][82] <-  "                      na = 0,"
  aa[[i]][93] <- "GEMSTAT_input_copier_multienh(lower_bounds=516,"
  aa[[i]][94] <-  "                              upper_bounds=516,"
  aa[[i]][95] <-   "                              coop=516, "
  aa[[i]][96] <-   "                              Enh_gene_map=516,"
  aa[[i]][97] <-   "                              free_fix=516,"
  aa[[i]][98] <-   "                              Gene_Exp=516, "
  aa[[i]][99] <-   "                              Motifs=516,"
  aa[[i]][100] <-   "                              Params=516, "
  aa[[i]][101] <-   "                              Params_from_output=T,"
  aa[[i]][102] <-   "                              seq=516,"
  aa[[i]][103] <-   "                              TFexp=numeric(0),"
  aa[[i]][104] <-   "                              TF_info=516,"
  aa[[i]][105] <-   "                              Treat_cont_map=516, "
  aa[[i]][106] <-   "                              wights=516, "
  aa[[i]][107] <-   "                              jobfile=516, "
  aa[[i]][108] <-   "                              .prev_na=5, "
  aa[[i]][109] <-   "                              .optim=F, "
  aa[[i]][110] <-   paste0("                              cur_exp=", 642+i, ",")
  aa[[i]][111] <-   "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")"
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  642 + i, ".R"))
}
setwd("..")
bash_batch_run_zip_transfer_process(exp_index = c(643:664))

for(i in 643:664){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[516]]$threshold)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[643:664], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

par(mfrow = c(1, 1), mar = c(4,4,4,4))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[516]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[516]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[516]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:13]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[516]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[516]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[516]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[516]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[643:664],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=13, model_index=integer(0), plot_res=T, fixed_ylim = F)

########################################################################################################################
########################################################################################################################
# doing conc modifications on exp 489
# creating all combinations of concentration modifications/ this time with the new ER motif
aa_modif_all_options_p20 <- list(list(1, -1, 21) ,
                                 list(10, -10, 31),
                                 list(c(13, 14), c(-13, -14), c(33, 34)),
                                 list(15,-15, 35),
                                 list(16,-16, 36),
                                 list(18,-18, 38))
aa <- GenerateComb(Set = c(1:6), TotalNumberOfTFs = 6, my_alphabet = c("1", "2", "3"))
aa_keep <- logical(nrow(aa))
for(i in 1:nrow(aa)){
  aa_keep[i] <- 3 %in% aa[i,]
}
aa <- aa[aa_keep, ]
aa_modif_all_comb_p20_raw <- aa
aa_modif_all_comb_p20 <- list()
for(i in 1:nrow(aa)){
  aa_modif_all_comb_p20[[i]] <- integer(0)
  for(j in 1:ncol(aa)){
    aa_modif_all_comb_p20[[i]] <- c(aa_modif_all_comb_p20[[i]], aa_modif_all_options_p20[[j]][aa[i, j]])
  }
  aa_modif_all_comb_p20[[i]] <- unlist(aa_modif_all_comb_p20[[i]])
}

# Following is the index of the chosen models
aa_model_index <- c(2,6,13,15,25,25,33,44,56,67,77,121,141)

# updating the data file
aa_name_conv <- names(TF.motifs.Shrinked.halfsites.count_2[abs(sort(unlist(aa_modif_all_options_p20)))[1:14]])
aa_name_conv <- c(aa_name_conv, names(TF.motifs.Shrinked.halfsites.count_2[abs(sort(unlist(aa_modif_all_options_p20)))[15:21] - 20]))
aa_name_conv[1:7] <- paste(aa_name_conv[1:7],"n",sep = "_")
aa_name_conv[8:14] <- paste(aa_name_conv[8:14],"p",sep = "_")
aa_name_conv[15:21] <- paste(aa_name_conv[15:21],"all",sep = "_")
aa_name_conv <- cbind(sort(unlist(aa_modif_all_options_p20)), aa_name_conv)
for(i in 1:length(aa_modif_all_comb_p20)){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[489,])
  GEMSTAT_experiment_info_df[664+i, 1] <- paste0("experiment_", 664+i)
  GEMSTAT_experiment_info_df[664+i, 8] <- paste("Concentration Modification. Base experiment 489.  my_model_index <- c(2,6,13,15,25,25,33,44,56,67,77,121,141).", "conc modif: ",
                                                paste(aa_name_conv[aa_name_conv[, 1] %in% as.character(aa_modif_all_comb_p20[[i]]), 2], collapse = " "))
  
}

boxplot.matrix(aa$Site_count_mat[[2]], las=2)
# Creating the R script to create the inputs and job files.
setwd("Seeded_GEMSTAT_ens/")
aa <- list()

for(i in 1:length(aa_modif_all_comb_p20)){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp489.R")
  aa[[i]][19] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer"
  for(j in 1:length(aa_modif_all_comb_p20[[i]])){
    if(aa_modif_all_comb_p20[[i]][j] < 0){
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[',-aa_modif_all_comb_p20[[i]][j],', ] <- rep(c(1,0), 45)')
    }else if(aa_modif_all_comb_p20[[i]][j] < 21){
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[',aa_modif_all_comb_p20[[i]][j],', ] <- rep(c(0,1), 45)')
    }else{
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[',aa_modif_all_comb_p20[[i]][j] - 20,', ] <- rep(1, 90)')
    }
  }
  
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 664 + i, ",")
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][90] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F"
  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=489,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=489,")
  aa[[i]] <- c(aa[[i]], "                              coop=489, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=489,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=489,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=489, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=489,")
  aa[[i]] <- c(aa[[i]], "                              Params=489, ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=F,")
  aa[[i]] <- c(aa[[i]], "                              seq=489,")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=489,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=489, ")
  aa[[i]] <- c(aa[[i]], "                              wights=489, ")
  aa[[i]] <- c(aa[[i]], "                              jobfile=489, ")
  aa[[i]] <- c(aa[[i]], "                              .prev_na=5, ")
  aa[[i]] <- c(aa[[i]], "                              .optim=T, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 664+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  # 
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  664 + i, ".R"))
}
setwd("..")
bash_batch_run_zip_transfer_process(exp_index = c(665:1329))

for(i in 665:1329){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
par(mfrow = c(1, 1), mar = c(6,4,4,4))
aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[665:1329], "[[", 7))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[489]]$Accuracy_All, aa), las= 2)
abline(h=seq(0, 1, 0.01), col = 2, lty = 4)
aaaa <- apply(X = aa, MARGIN = 2,which.max)
aaaa_max <- apply(X = aa, MARGIN = 2,max)
aaaa_max_ind <- sort(aaaa_max, decreasing = T, index.return = T)$ix
#aa_modif_all_comb_p20_mat <- do.call(what = rbind, aa_modif_all_comb_p20)
aa_modif_all_comb_p20_raw_sort <- aa_modif_all_comb_p20_raw[aaaa_max_ind,]
#aa_modif_all_comb_p20_raw_sort <- aa_modif_all_comb_p20_raw_sort[, -3]
colnames(aa_modif_all_comb_p20_raw_sort) <- c("AR", "GR", "PGR", "RARA", "RARG", "RXRA")

heatmap.2(aa_modif_all_comb_p20_raw_sort, Rowv = F, Colv = F, dendrogram = "none", trace = 'none')

aa_top <- aa_modif_all_comb_p20_raw_sort[1:60,]
aa_top_table <- numeric(6*3)
for(i in 1:ncol(aa_top)){
  aa_tab <- table(aa_top[, i])
  for(j in 1:length(aa_tab)){
    aa_top_table[3 * (i - 1) + as.numeric(names(aa_tab)[j])] <- aa_tab[j]
  }
}

barplot(aa_top_table, 
        col=c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3), rep(6, 3)), 
        names.arg = c("ARp", "ARn", "ARall", "GRp", "GRn","GRall","PGRp", "PGRn",
                      "PGRall","RARAp", "RARAn","RARAall","RARGp", "RARGn","RARGall",
                      "RXRAp", "RXRAn", "RXRAall"), las=2, main = "top60_experiments")
########################################################################################################################
########################################################################################################################
aaaa_max_ind[1:10] + 664
#  1042  973  759  691 1039  897  970  943  765  911
# KD of experiment 1042

for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[1042,])
  GEMSTAT_experiment_info_df[1329+i, 8] <- paste("in silico knockdown On Experiment 1042 of TFs ", 
                                                names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[1329+i, 1] <- paste0("experiment_", 1329+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[1042,])
GEMSTAT_experiment_info_df[1329+21, 8] <- "in silico knockdown On Experiment 1042 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[1329+21, 1] <- paste0("experiment_", 1329+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[1042,])
GEMSTAT_experiment_info_df[1329+22, 8] <- "in silico knockdown On Experiment 1042 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[1329+22, 1] <- paste0("experiment_", 1329+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp1042.R")
  if(i == 21){
    aa[[i]][27] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][27] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][27] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 1329 + i, ",")
  aa[[i]][82] <-  "                      na = 0,"
  aa[[i]][93] <- "GEMSTAT_input_copier_multienh(lower_bounds=1042,"
  aa[[i]][94] <-  "                              upper_bounds=1042,"
  aa[[i]][95] <-   "                              coop=1042, "
  aa[[i]][96] <-   "                              Enh_gene_map=1042,"
  aa[[i]][97] <-   "                              free_fix=1042,"
  aa[[i]][98] <-   "                              Gene_Exp=1042, "
  aa[[i]][99] <-   "                              Motifs=1042,"
  aa[[i]][100] <-   "                              Params=1042, "
  aa[[i]][101] <-   "                              Params_from_output=T,"
  aa[[i]][102] <-   "                              seq=1042,"
  aa[[i]][103] <-   "                              TFexp=numeric(0),"
  aa[[i]][104] <-   "                              TF_info=1042,"
  aa[[i]][105] <-   "                              Treat_cont_map=1042, "
  aa[[i]][106] <-   "                              wights=1042, "
  aa[[i]][107] <-   "                              jobfile=1042, "
  aa[[i]][108] <-   "                              .prev_na=5, "
  aa[[i]][109] <-   "                              .optim=F, "
  aa[[i]][110] <-   paste0("                              cur_exp=", 1329+i, ",")
  aa[[i]][111] <-   "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")"
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  1329 + i, ".R"))
}
setwd("..")
bash_batch_run_zip_transfer_process(exp_index = c(1330:1351))

for(i in 1330:1351){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1042]]$threshold)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[1330:1351], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

par(mfrow = c(1, 1), mar = c(4,4,4,4))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1042]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1042]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1042]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:13]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1042]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1042]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1042]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1042]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[1330:1351],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=13, model_index=integer(0), plot_res=T, fixed_ylim = F)
########################################################################################################################

# KD of experiment 973

for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[973,])
  GEMSTAT_experiment_info_df[1351+i, 8] <- paste("in silico knockdown On Experiment 973 of TFs ", 
                                                 names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[1351+i, 1] <- paste0("experiment_", 1351+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[973,])
GEMSTAT_experiment_info_df[1351+21, 8] <- "in silico knockdown On Experiment 973 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[1351+21, 1] <- paste0("experiment_", 1351+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[973,])
GEMSTAT_experiment_info_df[1351+22, 8] <- "in silico knockdown On Experiment 973 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[1351+22, 1] <- paste0("experiment_", 1351+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp973.R")
  if(i == 21){
    aa[[i]][27] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][27] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][27] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 1351 + i, ",")
  aa[[i]][82] <-  "                      na = 0,"
  aa[[i]][93] <- "GEMSTAT_input_copier_multienh(lower_bounds=973,"
  aa[[i]][94] <-  "                              upper_bounds=973,"
  aa[[i]][95] <-   "                              coop=973, "
  aa[[i]][96] <-   "                              Enh_gene_map=973,"
  aa[[i]][97] <-   "                              free_fix=973,"
  aa[[i]][98] <-   "                              Gene_Exp=973, "
  aa[[i]][99] <-   "                              Motifs=973,"
  aa[[i]][100] <-   "                              Params=973, "
  aa[[i]][101] <-   "                              Params_from_output=T,"
  aa[[i]][102] <-   "                              seq=973,"
  aa[[i]][103] <-   "                              TFexp=numeric(0),"
  aa[[i]][104] <-   "                              TF_info=973,"
  aa[[i]][105] <-   "                              Treat_cont_map=973, "
  aa[[i]][106] <-   "                              wights=973, "
  aa[[i]][107] <-   "                              jobfile=973, "
  aa[[i]][108] <-   "                              .prev_na=5, "
  aa[[i]][109] <-   "                              .optim=F, "
  aa[[i]][110] <-   paste0("                              cur_exp=", 1351+i, ",")
  aa[[i]][111] <-   "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")"
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  1351 + i, ".R"))
}
setwd("..")
bash_batch_run_zip_transfer_process(exp_index = c(1352:1373))

for(i in 1352:1373){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[973]]$threshold)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[1352:1373], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

par(mfrow = c(1, 1), mar = c(4,4,4,4))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[973]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[973]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[973]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:13]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[973]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[973]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[973]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[973]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[1352:1373],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=13, model_index=integer(0), plot_res=T, fixed_ylim = F)

########################################################################################################################
########################################################################################################################
# KD of exp 759

for(i in 1:20){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[759,])
  GEMSTAT_experiment_info_df[1373+i, 8] <- paste("in silico knockdown On Experiment 759 of TFs ", 
                                                 names(TF.motifs.Shrinked.halfsites.count)[i], ". No training was done (na=0)")
  GEMSTAT_experiment_info_df[1373+i, 1] <- paste0("experiment_", 1373+i)
}
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[759,])
GEMSTAT_experiment_info_df[1373+21, 8] <- "in silico knockdown On Experiment 759 of TFs PGR_s1, and PGR_s2. No training was done (na=0)"
GEMSTAT_experiment_info_df[1373+21, 1] <- paste0("experiment_", 1373+21)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[759,])
GEMSTAT_experiment_info_df[1373+22, 8] <- "in silico knockdown On Experiment 759 of TFs JUN_1, and JUN_2 No training was done (na=0)"
GEMSTAT_experiment_info_df[1373+22, 1] <- paste0("experiment_", 1373+22)

# create the Rscripts
setwd("Seeded_GEMSTAT_ens/")
aa <- list()
for(i in 1:22){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp759.R")
  if(i == 21){
    aa[[i]][27] <- 'aa_myTFexpmat[c(13, 14), ] <- 0'
  }else if(i == 22){
    aa[[i]][27] <- 'aa_myTFexpmat[c(6, 7), ] <- 0'
  }else{
    aa[[i]][27] <- paste0('aa_myTFexpmat[',i,', ] <- 0')
  }
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 1373 + i, ",")
  aa[[i]][82] <-  "                      na = 0,"
  aa[[i]][93] <- "GEMSTAT_input_copier_multienh(lower_bounds=759,"
  aa[[i]][94] <-  "                              upper_bounds=759,"
  aa[[i]][95] <-   "                              coop=759, "
  aa[[i]][96] <-   "                              Enh_gene_map=759,"
  aa[[i]][97] <-   "                              free_fix=759,"
  aa[[i]][98] <-   "                              Gene_Exp=759, "
  aa[[i]][99] <-   "                              Motifs=759,"
  aa[[i]][100] <-   "                              Params=759, "
  aa[[i]][101] <-   "                              Params_from_output=T,"
  aa[[i]][102] <-   "                              seq=759,"
  aa[[i]][103] <-   "                              TFexp=numeric(0),"
  aa[[i]][104] <-   "                              TF_info=759,"
  aa[[i]][105] <-   "                              Treat_cont_map=759, "
  aa[[i]][106] <-   "                              wights=759, "
  aa[[i]][107] <-   "                              jobfile=759, "
  aa[[i]][108] <-   "                              .prev_na=5, "
  aa[[i]][109] <-   "                              .optim=F, "
  aa[[i]][110] <-   paste0("                              cur_exp=", 1373+i, ",")
  aa[[i]][111] <-   "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")"
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  1373 + i, ".R"))
}
setwd("..")
bash_batch_run_zip_transfer_process(exp_index = c(1374:1395))

for(i in 1374:1395){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T,
                                                                                                             .round_thresh = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[759]]$threshold)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[1374:1395], "[[", 7))
colnames(aa) <- c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2")

par(mfrow = c(1, 1), mar = c(4,4,4,4))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[759]]$Accuracy_All, aa),
               las=2, main="allmodels_KD")
abline(h=seq(0,1, 0.05), col=2, lty=4, lwd=0.5)
aa2 <- aa - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[759]]$Accuracy_All
boxplot.matrix(aa2, las=2, main="all_tfKD_accuracy_diff")
abline(h=seq(-0.5,0.5, 0.05), col=2, lty=4, lwd=0.5)
aa_goodind <- sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[759]]$Accuracy_All,
                   decreasing = T, index.return=T)$ix[1:13]
aa_gmod <- aa[aa_goodind, ]
aa2_gmod <- aa_gmod - GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[759]]$Accuracy_All[aa_goodind]
boxplot.matrix(aa2_gmod, las=2, main="top15_tfKD_accuracy_diff")
abline(h=seq(-1,1, 0.05), col=2, lty=4, lwd=0.5)

par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[759]]$Accuracy_All[aa_goodind],
       aa_gmod[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}
par(mfrow=c(4, 6), mar = c(1, 1, 4, 1))
for(i in 1:22){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[759]]$Accuracy_All,
       aa[, i], xlim = c(0.4, 0.7), ylim = c(0.4, 0.7), main = colnames(aa_gmod)[i])
  abline(c(0, 1), col=2)
}

aa <- KD_analyzer(base_experiment = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[759]],
                  KD_experiments = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[1374:1395],
                  KD_names = c(names(TF.motifs.Shrinked.halfsites.count), "PGR_s1 + PGR_s2", "JUN_1 + JUN_2"),
                  top_nu=13, model_index=integer(0), plot_res=T, fixed_ylim = F)
########################################################################################################################
########################################################################################################################
# perform KD experiment


GEMSTAT_based_on_linear_exp_KD_Results <- list()
for(i in 759:1329){
  print(paste0("Experiment ", i))
  KD_Performer(parent_Exp=i, 
               TF_index=list(3, 13, 14, c(13, 14)),
               shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
               prev_na=5,
               run_GEMSTAT = T)
  GEMSTAT_based_on_linear_exp_KD_Results[[i]] <- read_KD_results(experiment_nu = i, 
                                        shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens")
}

########################################################################################################################
########################################################################################################################
# Perform KD for the 64 0, 1 experiments
for(i in 491:554){
  print(paste0("Experiment ", i))
  KD_Performer(parent_Exp=i, 
               TF_index=list(3, 13, 14, c(13, 14)),
               shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
               prev_na=5,
               run_GEMSTAT = T)
  GEMSTAT_based_on_linear_exp_KD_Results[[i]] <- read_KD_results(experiment_nu = i, 
                                                                 shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens")
}
########################################################################################################################
# Find Experiments where ER KD causes significant change in performance and the accuracy is high in WT
aa_WT_acc <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(c(491:554), c(665:1329))], "[[", 7))
aa_ER_KD_Acc <- do.call(cbind, lapply(lapply(X = GEMSTAT_based_on_linear_exp_KD_Results[c(c(491:554), c(665:1329))], "[[", 4), "[[", 7))
aa_acc_diff <- aa_WT_acc - aa_ER_KD_Acc
aa_acc_diff_6 <- aa_acc_diff
aa_acc_diff_6[aa_WT_acc < 0.6] <- NA
boxplot.matrix(aa_acc_diff_6, las = 2)
sum(aa_acc_diff_6>= 0.06, na.rm =T) 
colnames(aa_acc_diff_6) <- c(c(491:554), c(665:1329))
colnames(aa_WT_acc) <- c(c(491:554), c(665:1329))
colnames(aa_ER_KD_Acc) <- c(c(491:554), c(665:1329))

aa_acc_diff_6_10 <- aa_acc_diff_6
aa_acc_diff_6_6 <- aa_acc_diff_6

aa_acc_diff_6_10[aa_acc_diff_6_10 < 0.1] <- NA
aa_acc_diff_6_6[aa_acc_diff_6_6 < 0.06] <- NA

aa_acc_diff_6_10_f <- aa_acc_diff_6_10[,! colSums(is.na(aa_acc_diff_6_10)) == nrow(aa_acc_diff_6_10)]
aa_acc_diff_6_6_f  <- aa_acc_diff_6_6[,! colSums(is.na(aa_acc_diff_6_6)) == nrow(aa_acc_diff_6_6)]

GEMSTAT_experiment_info_df[as.numeric(colnames(aa_acc_diff_6_10_f)) , 8]
colnames(aa_modif_all_comb_p20_raw)<-c("AR", "GR", "PGR", "RARA", "RARG", "RXRA")
aa_modif_all_comb_mat <- do.call(rbind, aa_modif_all_comb)
aa_modif_all_comb_mat[aa_modif_all_comb_mat > 0] <- 1
aa_modif_all_comb_mat[aa_modif_all_comb_mat < 0] <- 2
aa_modif_all_comb_mat <- aa_modif_all_comb_mat[, -c(3)]
colnames(aa_modif_all_comb_mat)<-c("AR", "GR", "PGR", "RARA", "RARG", "RXRA")
aa_modif_all_comb_mat_all <- rbind(aa_modif_all_comb_mat, aa_modif_all_comb_p20_raw)
rownames(aa_modif_all_comb_mat_all) <- c(c(491:554), c(665:1329))


heatmap.2(aa_modif_all_comb_mat_all[colnames(aa_acc_diff_6_10_f),]
          , Rowv = T, Colv = T, dendrogram = "both",
          trace = 'none', distfun = my_dist_fuc)
heatmap.2(aa_modif_all_comb_mat_all[colnames(aa_acc_diff_6_6_f),]
          , Rowv = T, Colv = T, dendrogram = "both",
          trace = 'none', distfun = my_dist_fuc)

aa_top <- aa_modif_all_comb_mat_all[colnames(aa_acc_diff_6_6_f),]
aa_top_table <- numeric(6*3)
for(i in 1:ncol(aa_top)){
  aa_tab <- table(aa_top[, i])
  for(j in 1:length(aa_tab)){
    aa_top_table[3 * (i - 1) + as.numeric(names(aa_tab)[j])] <- aa_tab[j]
  }
}
par(mfrow = c(1, 1), mar = rep(6, 4))
barplot(aa_top_table, 
        col=c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3), rep(6, 3)), 
        names.arg = c("ARp", "ARn", "ARall", "GRp", "GRn","GRall","PGRp", "PGRn",
                      "PGRall","RARAp", "RARAn","RARAall","RARGp", "RARGn","RARGall",
                      "RXRAp", "RXRAn", "RXRAall"), las=2, main = "WT acc > 0.6, ER KD drop > 0.06")

par(mfrow = c(4, 4), mar = c(2, 3, 1, 1))
for(i in 1:ncol(aa_acc_diff_6_10_f)){
  plot(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]], 
       aa_ER_KD_Acc[, colnames(aa_acc_diff_6_10_f)[i]],
       xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD")
  abline(c(0, 1), col=2)
}
# what are the parameters for the models 
aa_mod_ind <- list()
for(i in 1:ncol(aa_acc_diff_6_10_f)){
  aa_mod_ind[[i]] <- which(! is.na(aa_acc_diff_6_10_f[, i]))
}
names(aa_mod_ind) <- colnames(aa_acc_diff_6_10_f)

# look at the parameters
aa_bindings <- matrix(nrow = length(unlist(aa_mod_ind)), ncol =20)
colnames(aa_bindings) <- names(GEMSTAT_based_on_linear_exp_pars[[691]]$binding[8,])
rownames(aa_bindings) <- character(nrow(aa_bindings))

aa_alphas <- matrix(nrow = length(unlist(aa_mod_ind)), ncol =20)
colnames(aa_alphas) <- names(GEMSTAT_based_on_linear_exp_pars[[691]]$alpha_effective[8,])
rownames(aa_alphas) <- character(nrow(aa_alphas))

aa_coops <- matrix(nrow = length(unlist(aa_mod_ind)), ncol =7)
colnames(aa_coops) <- names(GEMSTAT_based_on_linear_exp_pars[[691]]$coop[[8]])
rownames(aa_coops) <- character(nrow(aa_coops))

aa_cnt = 1
for(i in 1:length(aa_mod_ind)){
  for(j in 1:length(aa_mod_ind[[i]])){
    aa_bindings[aa_cnt, ] <- GEMSTAT_based_on_linear_exp_pars[[as.numeric(names(aa_mod_ind))[i]]]$binding[aa_mod_ind[[i]][j], ]
    aa_alphas[aa_cnt, ] <- GEMSTAT_based_on_linear_exp_pars[[as.numeric(names(aa_mod_ind))[i]]]$alpha_effective[aa_mod_ind[[i]][j], ]
    aa_coops[aa_cnt, ] <- GEMSTAT_based_on_linear_exp_pars[[as.numeric(names(aa_mod_ind))[i]]]$coop[[aa_mod_ind[[i]][j]]]
    rownames(aa_bindings)[aa_cnt] <- paste0(as.numeric(names(aa_mod_ind))[i], "_", aa_mod_ind[[i]][j])
    rownames(aa_alphas)[aa_cnt] <- paste0(as.numeric(names(aa_mod_ind))[i], "_", aa_mod_ind[[i]][j])
    rownames(aa_coops)[aa_cnt] <- paste0(as.numeric(names(aa_mod_ind))[i], "_", aa_mod_ind[[i]][j])
    aa_cnt <- aa_cnt + 1
  }
  
}
par(mfrow = c(1, 1), mar = c(6,4,4,4))
boxplot.matrix(aa_bindings, las = 2)
boxplot.matrix(aa_alphas, las = 2)
boxplot.matrix(aa_coops, las = 2)

aa_exp_names <- unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 1))
aa_col <- sample(x = col_vector, size = 5, replace = F)
heatmap.2(cbind(aa_alphas, aa_modif_all_comb_mat_all[aa_exp_names,])
          , Rowv = T, Colv = F, dendrogram = "none",
          trace = 'none', breaks = c(0, 1 , 2, 3, 7, 10),
          margins = c(8,8), col = col_vector[11:15], main = "alphaERsensitive+concComb")

heatmap.2(aa_coops
          , Rowv = T, Colv = T, dendrogram = "none",
          trace = 'none', breaks = c(100, 300, 500, 700, 1000),
          margins = c(8,8), col = col_vector[13:16], main="CoopERsensitive")
heatmap.2(cbind(aa_bindings, aa_modif_all_comb_mat_all[aa_exp_names,])
          , Rowv = T, Colv = F, dendrogram = "none",
          trace = 'none', breaks = c(0, 1, 2, 3, 10, 50),
          col = col_vector[11:15], margins = c(8,8), 
          main="BindingERsensitive")

# pair the parameters and concentration profiles to see how does that effect the parameters

# are these models also sensitive to PGR KD?
aa_PGR_KD_Acc <- do.call(cbind, lapply(lapply(X = GEMSTAT_based_on_linear_exp_KD_Results[c(c(491:554), c(665:1329))], "[[", 2), "[[", 7))
aa_PGR1_KD_Acc <- do.call(cbind, lapply(lapply(X = GEMSTAT_based_on_linear_exp_KD_Results[c(c(491:554), c(665:1329))], "[[", 1), "[[", 7))
aa_PGR2_KD_Acc <- do.call(cbind, lapply(lapply(X = GEMSTAT_based_on_linear_exp_KD_Results[c(c(491:554), c(665:1329))], "[[", 3), "[[", 7))

colnames(aa_PGR_KD_Acc) <- c(c(491:554), c(665:1329))
colnames(aa_PGR1_KD_Acc) <- c(c(491:554), c(665:1329))
colnames(aa_PGR2_KD_Acc) <- c(c(491:554), c(665:1329))


par(mfrow = c(4, 4), mar = c(2, 3, 4, 1))
plot(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[1]], 
     aa_ER_KD_Acc[, colnames(aa_acc_diff_6_10_f)[1]],
     xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", pch=1, col=2)
points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[1]],
       aa_PGR_KD_Acc[, colnames(aa_acc_diff_6_10_f)[1]],
       col = 4 , pch = 2)
abline(c(0, 1), col=1, lwd = 0.6)
par(xpd=T)
legend(legend = c("ER_KD", "PGR_KD"),
       fill = c(2, 4), x =0.4, y=1
       #, pch = c(19, 16)
       ,bty = "n", y.intersp = 0.6, cex = 0.7)
par(xpd=F)
for(i in 2:ncol(aa_acc_diff_6_10_f)){
  plot(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]], 
       aa_ER_KD_Acc[, colnames(aa_acc_diff_6_10_f)[i]],
       xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", pch=1, col=2)
  points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]],
         aa_PGR_KD_Acc[, colnames(aa_acc_diff_6_10_f)[i]],
         col = 4 , pch = 2)
  # points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        aa_PGR1_KD_Acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        col = 4 , pch = 19)
  # points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        aa_PGR2_KD_Acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        col = 6, pch = 17)
  
  abline(c(0, 1), col=1, lwd = 0.6)
  }

# Perform All other KDs on the 15 experiments identified as ER sensitive


# Perform Cooperativity Knock down to see if its actually important
aa_coop_kd <- list()
aa_cnt = 0
for(i in as.numeric(colnames(aa_acc_diff_6_10_f))){
  aa_cnt =aa_cnt + 1
  # KD_Performer_by_params(parent_Exp=i,
  #                        TF_names = names(TF.motifs.Shrinked.halfsites.count_2),
  #                        Coop_mat_list=list(c(3, 3), c(13, 14)),
  #                        shared_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
  #                        prev_na=5,
  #                        run_GEMSTAT=T)
  aa_coop_kd[[aa_cnt]] <- read_KD_results(experiment_nu=i, 
                                shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens", 
                                coop_KD=T)
}
names(aa_coop_kd) <- colnames(aa_acc_diff_6_10_f)
GEMSTAT_based_on_linear_exp_KD_COOP_Results <- list()
for(i in 1:length(aa_coop_kd)){
  GEMSTAT_based_on_linear_exp_KD_COOP_Results[[as.numeric(names(aa_coop_kd))[i]]] <- aa_coop_kd[[i]]
}
aa_coop_ER_kd_acc <- do.call(cbind, lapply(lapply(X = aa_coop_kd, "[[", 1), "[[", 7))
aa_coop_PGR_kd_acc <- do.call(cbind, lapply(lapply(X = aa_coop_kd, "[[", 2), "[[", 7))

par(mfrow = c(4, 4), mar = c(2, 3, 5, 1))
# plot(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[1]],
#      aa_ER_KD_Acc[, colnames(aa_acc_diff_6_10_f)[1]],
#      xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD",
#      col=2, pch=1)
plot(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[1]],
     aa_PGR_KD_Acc[, colnames(aa_acc_diff_6_10_f)[1]],
     xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", col = 4, pch=1)

# points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[1]],
#        aa_coop_ER_kd_acc[, 1],
#        col = 3, pch =4)
points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[1]],
       aa_coop_PGR_kd_acc[, 1],
       col = 6, pch = 4)
abline(c(0, 1), col=1, lwd = 0.6)
par(xpd=T)
legend(legend = c("PGR_KD", "PGR_coop_KD"),
       fill = c(4, 6), x =0.4, y=1.25
       #, pch = c(19, 16)
       ,bty = "n", y.intersp = 0.6, cex = 0.9)
par(xpd=F)
 for(i in 2:ncol(aa_acc_diff_6_10_f)){
  # plot(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #      aa_ER_KD_Acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #      xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", col=2, pch=1)
  plot(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]],
       aa_PGR_KD_Acc[, colnames(aa_acc_diff_6_10_f)[i]],
       xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", col = 4, pch=1)
  # points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        aa_PGR_KD_Acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        col = 6 , pch = 19, cex = 0.6)
  # points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        aa_PGR1_KD_Acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        col = 4 , pch = 19)
  # points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        aa_PGR2_KD_Acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        col = 6, pch = 17)
  # points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        aa_coop_ER_kd_acc[, i],
  #        col = 3, pch = 4)
  # plot(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]],
  #        aa_coop_ER_kd_acc[, i],
  #        col = 3, pch = 17, xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD")
  points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[i]],
         aa_coop_PGR_kd_acc[, i],
         col = 6, pch = 4)
  abline(c(0, 1), col=1, lwd = 0.6)
}

# compare the results of TF KD with coop KD
# Ideally the results should be the same as the TF Knockdown
########################################################################################################################
########################################################################################################################
# new experiment with lower upper bounds for binding and higher lower bound for cooperativity
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[397,])
GEMSTAT_experiment_info_df[1396, 8] <- paste("same as Experiment 397 but changing the bounds on binding and cooperativity of dimers. binding lb=0.00001, binding ub=0.001, coop lb = 1000000, coop ub=1000250")
GEMSTAT_experiment_info_df[1396, 1] <- paste0("experiment_", 1396)


for(i in 1396:1396){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
                                                                                                             
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
# Perform TF KD
KD_Performer(parent_Exp=1396, 
             TF_index=list(3, 13, 14, c(13, 14)),
             shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
             prev_na=5,
             run_GEMSTAT = F)
GEMSTAT_based_on_linear_exp_KD_Results[[1396]] <- read_KD_results(experiment_nu = 1396, 
                                                               shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens")

# Perform coop KD
KD_Performer_by_params(parent_Exp=1396,
                       TF_names = names(TF.motifs.Shrinked.halfsites.count_2),
                       Coop_mat_list=list(c(3, 3), c(13, 14)),
                       shared_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                       prev_na=5,
                       run_GEMSTAT=T)
GEMSTAT_based_on_linear_exp_KD_COOP_Results[[1396]] <-  read_KD_results(experiment_nu=1396, 
                                                                shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens", 
                                                                coop_KD=T)


par(mfrow = c(1, 1), mar = c(4,4,4,4))
plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1396]]$Accuracy_All, 
     GEMSTAT_based_on_linear_exp_KD_COOP_Results[[1396]]$ESR1_ESR1$Accuracy_All,
     xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", xlab = "WT", pch = 4, col = 2)
points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1396]]$Accuracy_All,
       GEMSTAT_based_on_linear_exp_KD_Results[[1396]]$`3`$Accuracy_All, col = 2, pch = 1)

points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1396]]$Accuracy_All, 
     GEMSTAT_based_on_linear_exp_KD_COOP_Results[[1396]]$PGR_s1_PGR_s2$Accuracy_All,
     col = 3, pch = 4
#     ,xlim = c(0.2, 0.7), ylim = c(0.2, 0.7)
     )
points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1396]]$Accuracy_All,
       GEMSTAT_based_on_linear_exp_KD_Results[[1396]]$`13_14`$Accuracy_All, col = 3, pch = 1)
abline(c(0, 1), col=1, lwd = 0.6)
par(xpd=T)
legend(legend = c("ER_KD", "ER_coop_KD","PGR_KD", "PGR_coop_KD"),
       col = c(2, 2, 3, 3), pch = c(1, 4, 1, 4),x="topleft"
       #, pch = c(19, 16)
       ,bty = "n", y.intersp = 0.6, cex = 0.9)
par(xpd=F)



########################################################################################################################
########################################################################################################################
# create a TF to gene map --> for each of the 13 models: create a matrix whereeach row is a TF, each column is a gene and
#  the entry shows number of binding sites above specified threshold
aa_homodimer <- c(T, F, T, F, F, F, F, F, F, T, F, F, F, F, T, T, F, F, F, F)
names(aa_homodimer) <- names(TF.motifs.Shrinked.halfsites.count_2)
aa_dimer_orien <- c(2, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0)
names(aa_dimer_orien) <- names(TF.motifs.Shrinked.halfsites.count_2)
aa_homoDimer_distance <- c(3, 0, 3, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 5, 6, 0, 0, 0, 0)
names(aa_homoDimer_distance) <- names(TF.motifs.Shrinked.halfsites.count_2)
TF.motifs.Shrinked.halfsites.count_2_maxLLR <- c(6.90008, 11.0531, 6.2838, 10.705, 7.12241, 12.2036,
                                                 10.3949, 16.4493, 12.0708, 6.94399, 8.91828, 11.544,
                                                 3.76517, 6.73403, 7.49103, 7.95261, 7.64841, 7.49922,
                                                 10.0583, 10.7687)
aa_newthr3 <- c(0.17, 0.71, 0.6, 0.65, 0.22,0.48, 0.6, 0.70, 0.65, 0.19, 0.52, 0.81,
                0.14, 0.16, 0.44, 0.48, 0.38, 0.34, 0.57, 0.57)

names(aa_newthr3) <- names(TF.motifs.Shrinked.halfsites.count)
aa <- count_site_from_annotation(annotation_list = aa_GEMSTAT_based_on_linear_exp_annotation_305, 
                                 TF_names = names(TF.motifs.Shrinked.halfsites.count_2),
                                 TF_index=c(13, 14),
                                 homoDimer = aa_homodimer, 
                                 dimer_orientation = aa_dimer_orien,
                                 homoDimer_distance=aa_homoDimer_distance,
                                 heterodimer_pair = list(c(13, 14)),
                                 heterodimer_distance = 3,
                                 MAXLLR = TF.motifs.Shrinked.halfsites.count_2_maxLLR,
                                 annotation_thresh = aa_newthr3)
aa_site_Count <- aa
par(mfrow = c(1, 1), mar= c(4, 4,4 ,4))
boxplot.matrix(aa_site_Count$Site_count_mat[[13]], las =2)

par(mfrow = c(5, 4), mar= c(2, 2,2 ,1))
for(i in 1:length(aa_site_Count$Sites_Per_TF)){
  hist(unlist(lapply(aa_site_Count$Sites_Per_TF[[i]][[25]], nrow)), main = names(TF.motifs.Shrinked.halfsites.count_2)[i])
}
par(mfrow = c(5, 4), mar= c(2, 2,2 ,1))
for(i in 1:length(aa_site_Count$Sites_Per_TF)){
  hist(aa_site_Count$Site_count_mat[[25]][, i], main = names(TF.motifs.Shrinked.halfsites.count_2)[i])
}


c(2,6,13,15,25,25,33,44,56,67,77,121,141)
par(mfrow = c(1, 1), mar= c(4, 4,4 ,4))
heatmap.2(aa_site_Count$Site_count_mat[[13]],
          Rowv = T, Colv = T, dendrogram = "both", trace = "none",
          breaks = c(0, 1, 2, 3, 5, 10, 20), col = col_vector[11:16])
ncol(aa_site_Count$Site_count_mat[[1]])
########################################################################################################################
########################################################################################################################
# Experiment 1397
# Experiment with New set of Motifs: All from HocoMoco
# PGR is not dimer anymore
# only one TF for JUN
# Gene expression matrix is changed to union of the ones who did well on the WT linear experiment and 
#  the ones who did well in GR modification experiment

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[397,])
GEMSTAT_experiment_info_df[1397, 8] <- paste("Experiment with New set of Motifs: All from HocoMoco, PGR is not dimer anymore, only one TF for JUN, Gene expression matrix is changed to union of the ones who did well on the WT linear experiment and the ones who did well in GR modification experiment, bounds are different, didn't actually perform this yet")
GEMSTAT_experiment_info_df[1397, 6] <- "AR-AR_DIMER_dist_4|ER-ER_DIMER_dist_4|NR3C1-NR3C1_DIMER_dist_4|RARA-RARA_DIMER_dist_6|RARG-RARG_DIMER_dist_6|RXRA-RXRA_DIMER_dist_2"
GEMSTAT_experiment_info_df[1397, 1] <- paste0("experiment_", 1397)
GEMSTAT_experiment_info_df[1397,10] <- 48
my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52_gene_gte3WTGRcommon_index <- sort(union(CS598_gene_index, 
                           my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52_gene_gte3common_index))

my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_48_gene_gte3WTGRcommon <- my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52[my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_52_gene_gte3WTGRcommon_index,]
TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_ER01_monomer_hocomoco <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_dJun_ER01_monomer[c(c(1:6), c(8:13), c(15:20)),]
rownames(TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_ER01_monomer_hocomoco)[12] <- "PGR"
GEMSTAT_based_on_linear_exp_annotation_1397 <- list()
for(i in 1:148){
  print(i)
  GEMSTAT_based_on_linear_exp_annotation_1397[[i]] <- annotation_reader(annot_file = paste0("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_1397_annotation/Annotation/Experiment_1397_", i, ".annot"), 
                                                                        TF_names = names(TF.motifs.Shrinked.hocomoco.count))
  
}

aa <- count_site_from_annotation(annotation_list = GEMSTAT_based_on_linear_exp_annotation_1397, 
                                 TF_names = names(TF.motifs.Shrinked.hocomoco.count),
                                 TF_index=c(1:18),
                                 homoDimer = aa_homodimer, 
                                 dimer_orientation = aa_dimer_orien,
                                 homoDimer_distance=aa_homoDimer_distance,
                                 heterodimer_pair = list(c(13, 14)),
                                 heterodimer_distance = 3,
                                 MAXLLR = TF.motifs.Shrinked.halfsites.count_2_maxLLR,
                                 annotation_thresh = aa_newthr3)

for(i in 1397:1397){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}



########################################################################################################################
# experiment 1398
# same as exp 1397 only with manual annotations fed in
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[1397,])
GEMSTAT_experiment_info_df[1398, 8] <- paste("same as exp 1397 only with manual annotations fed in")
GEMSTAT_experiment_info_df[1398, 1] <- paste0("experiment_", 1398)

aa_homodimer <-   c(T, F, T, F, F, F, F, F, T, F, F, F, T, T, F, T, F, F)
names(aa_homodimer) <- names(TF.motifs.Shrinked.hocomoco.count)
aa_dimer_orien <-        c(2, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 0, 1, 0, 0)
names(aa_dimer_orien) <- names(TF.motifs.Shrinked.hocomoco.count)
aa_homoDimer_distance <- c(3, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 5, 5, 0, 1, 0, 0)
names(aa_homoDimer_distance) <- names(TF.motifs.Shrinked.hocomoco.count)
TF.motifs.Shrinked.hocomoco.count_maxLLR <- c(6.13309, 10.0693, 6.28522, 11.3454, 7.91607, 9.98706,
                                              8.40569, 12.5907, 7.2447, 10.2185, 13.5565, 11.9268,
                                              6.69131, 7.37613, 9.94809, 6.93567, 12.5285, 9.73563)
names(TF.motifs.Shrinked.hocomoco.count_maxLLR) <- names(TF.motifs.Shrinked.hocomoco.count)
aa_newthr3 <- c(0.17, 0.71, 0.6, 0.65, 0.22,0.48, 0.6, 0.70, 0.65, 0.19, 0.52, 0.81,
                0.14, 0.16, 0.44, 0.48, 0.38, 0.34)

names(aa_newthr3) <- names(TF.motifs.Shrinked.hocomoco.count)
aa <- count_site_from_annotation_combined(annotation_list = GEMSTAT_based_on_linear_exp_annotation_1397, 
                                          TF_names = names(TF.motifs.Shrinked.hocomoco.count),
                                          TF_index=c(1:18),
                                          homoDimer = aa_homodimer, 
                                          dimer_orientation = aa_dimer_orien,
                                          homoDimer_distance=aa_homoDimer_distance,
                                          heterodimer_pair = numeric(0),
                                          heterodimer_distance = numeric(0),
                                          MAXLLR = TF.motifs.Shrinked.hocomoco.count_maxLLR,
                                          annotation_thresh = aa_newthr3,
                                          LLR_to_pVal_list = TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval,
                                          pVal_thresh = 0.0001)
GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m4 <- aa
#boxplot.matrix(aa$Site_count_mat[[13]], las = 2)
aa <- count_site_from_annotation_combined(annotation_list = GEMSTAT_based_on_linear_exp_annotation_1397, 
                                          TF_names = names(TF.motifs.Shrinked.hocomoco.count),
                                          TF_index=c(1:18),
                                          homoDimer = aa_homodimer, 
                                          dimer_orientation = aa_dimer_orien,
                                          homoDimer_distance=aa_homoDimer_distance,
                                          heterodimer_pair = numeric(0),
                                          heterodimer_distance = numeric(0),
                                          MAXLLR = TF.motifs.Shrinked.hocomoco.count_maxLLR,
                                          annotation_thresh = aa_newthr3,
                                          LLR_to_pVal_list = TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval,
                                          pVal_thresh = 0.001)
GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3 <- aa
save(list = c("GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3",
              "GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m4",
              "GEMSTAT_init_BOlinear",
              "annotation_Writer", 
              "GEMSTAT_input_copier_multienh"), file = "Annotation_manual.RData")


for(i in 1398:1398){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}
# Perform TF KD
KD_Performer(parent_Exp=1398, 
             TF_index=as.list(c(1:18)),
             shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
             prev_na=5,
             run_GEMSTAT = T)
GEMSTAT_based_on_linear_exp_KD_Results[[1398]] <- read_KD_results(experiment_nu = 1398, 
                                                                  shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens")

# Perform coop KD
KD_Performer_by_params(parent_Exp=1398,
                       TF_names = names(TF.motifs.Shrinked.hocomoco.count),
                       Coop_mat_list=list(c(1, 1), c(3,3),c(9, 9), c(13, 13), c(14, 14), c(16, 16)),
                       shared_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                       prev_na=5,
                       run_GEMSTAT=T)
GEMSTAT_based_on_linear_exp_KD_COOP_Results[[1398]] <-  read_KD_results(experiment_nu=1398, 
                                                                        shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens", 
                                                                        coop_KD=T)

par(mfrow = c(1, 1), mar = c(4,4,4,4))
plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All, 
     GEMSTAT_based_on_linear_exp_KD_COOP_Results[[1398]]$ESR1_ESR1$Accuracy_All,
     xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", xlab = "WT", pch = 4, col = 2)
points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All,
       GEMSTAT_based_on_linear_exp_KD_Results[[1398]]$`3`$Accuracy_All, col = 2, pch = 1)

points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All, 
       GEMSTAT_based_on_linear_exp_KD_COOP_Results[[1398]]$AR_AR$Accuracy_All,
       col = 3, pch = 4
       #     ,xlim = c(0.2, 0.7), ylim = c(0.2, 0.7)
)
points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All,
       GEMSTAT_based_on_linear_exp_KD_Results[[1398]]$`1`$Accuracy_All, col = 3, pch = 1)
abline(c(0, 1), col=1, lwd = 0.6)
par(xpd=T)
legend(legend = c("ER_KD", "ER_coop_KD","PGR_KD", "PGR_coop_KD"),
       col = c(2, 2, 3, 3), pch = c(1, 4, 1, 4),x="topleft"
       #, pch = c(19, 16)
       ,bty = "n", y.intersp = 0.6, cex = 0.9)
par(xpd=F)



boxplot.matrix(do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_KD_COOP_Results[[1398]], "[[", 7)), las = 2)
aa_x <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_KD_Results[[1398]], "[[", 7))
colnames(aa_x) <- names(TF.motifs.Shrinked.hocomoco.count)
boxplot.matrix(aa_x, las = 2)

par(mfrow = c(4, 5), mar = c(1,1,3,1))
for(i in 1:length(GEMSTAT_based_on_linear_exp_KD_Results[[1398]])){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All, 
       GEMSTAT_based_on_linear_exp_KD_Results[[1398]][[i]]$Accuracy_All,
       xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", xlab = "WT", pch = 1, col = 2, main=names(TF.motifs.Shrinked.hocomoco.count)[i])
  if(i %in% c(1, 3, 9, 13, 14, 16)){
    points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All,
           GEMSTAT_based_on_linear_exp_KD_COOP_Results[[1398]][[which(c(1, 3, 9, 13, 14, 16) %in% i)]]$Accuracy_All, col = 2, pch = 4)
  }

  abline(c(0, 1), col=1, lwd = 0.6)
}


# binding site heatmaps of top 25 models
GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3
c(2,6,13,15,25,25,33,44,56,67,77,121,141)
aa_model_index_1398
par(mfrow = c(1, 1), mar= c(4, 4,4 ,4))
heatmap.2(aa_site_Count$Site_count_mat[[aa_model_index_1398[3]]],
          Rowv = F, Colv = F, dendrogram = "none", trace = "none",
          breaks = c(0, 0.5, 1.5, 3, 5, 10, 20), col = col_vector[11:16])
heatmap.2(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat[[aa_model_index_1398[3]]],
          Rowv = F, Colv = F, dendrogram = "none", trace = "none",
          breaks = c(0, 0.5, 1.5, 3, 5, 10, 20), col = col_vector[11:16])
aa_nz <- numeric(length(aa_model_index_1398))
for(i in 1:length(aa_model_index_1398)){
  aa_nz[i] <- sum(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat[[aa_model_index_1398[i]]][, 3] > 0)/
    nrow(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat[[aa_model_index_1398[i]]])
}
hist(aa_nz, breaks = 10)

aa_nz_all <- numeric(148)
for(i in 1:148){
  aa_nz_all[i] <- sum(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat[[i]][, 3] > 0)/
    nrow(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat[[i]])
}

aa_dat <- data.frame( x=aa_nz_all, above_60=aa_nz_all>0.6 )
library(ggplot2)
qplot(x,data=aa_dat,geom="histogram",fill=above_60, cex = 3)

sum(aa_nz_all > 0.60)
aa_model_index_1398_ERBS <- which(aa_nz_all > 0.6)

par(mfrow = c(4, 5), mar = c(1,1,3,1))
for(i in 1:length(GEMSTAT_based_on_linear_exp_KD_Results[[1398]])){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All[aa_nz_all > 0.6], 
       GEMSTAT_based_on_linear_exp_KD_Results[[1398]][[i]]$Accuracy_All[aa_nz_all > 0.6],
       xlim = c(0.3, 0.65), ylim = c(0.3, 0.65), ylab = "KD", xlab = "WT", pch = 1, col = 4,
       main=names(TF.motifs.Shrinked.hocomoco.count)[i])
  if(i %in% c(1, 3, 9, 13, 14, 16)){
    points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All[aa_nz_all > 0.6],
           GEMSTAT_based_on_linear_exp_KD_COOP_Results[[1398]][[which(c(1, 3, 9, 13, 14, 16) %in% i)]]$Accuracy_All[aa_nz_all > 0.6], col = 2, pch = 4)
  }
  
  abline(c(0, 1), col=1, lwd = 0.6)
}

########################################################################################################################
########################################################################################################################
# experiment 1399: just like 1398 but with a subset of models : aa_model_index_1398_all
aaaa <-sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All, decreasing = T, index.return=T)$ix[1:25]
aa_model_index_1398 <- aaaa
aa_model_index_1398_ERBS
aa_model_index_1398_all <- union(aa_model_index_1398, aa_model_index_1398_ERBS)

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[1398,])
GEMSTAT_experiment_info_df[1399, 1] <- paste0("experiment_", 1399)
GEMSTAT_experiment_info_df[1399, 8] <- "Same as exp 1398, only the model inedx is changed. my_model_index <- c(21,56,61,142,49,57,115,9,134,117,127,94,55,118,72,1,30,51,77,113,13,27,34,84,74,3,11,15,16,20,33,42,43,44,47,64,80,92,93,101,108,120,121,130,132,137)."
                                               


########################################################################################################################
########################################################################################################################
# receptor activation profile (729) 01, 10, 11 for the top 25 models of experiment 1399
# doing conc modifications on exp 1398
# creating all combinations of concentration modifications/ this time with the new ER motif
aa_modif_all_options_p20_1398 <- list(list(1, -1, 21),
                                 list(9, -9, 29),
                                 list(12, -12, 32),
                                 list(13,-13, 33),
                                 list(14,-14, 34),
                                 list(16,-16, 36))
aa <- GenerateComb(Set = c(1:6), TotalNumberOfTFs = 6, my_alphabet = c("1", "2", "3"))
aa_modif_all_comb_p20_raw_1398 <- aa

heatmap.2(aa_modif_all_comb_p20_raw_1398)
rownames(aa_modif_all_comb_p20_raw_1398) <- c(1400:2128)
aa_modif_all_comb_p20_1398 <- list()
for(i in 1:nrow(aa)){
  aa_modif_all_comb_p20_1398[[i]] <- integer(0)
  for(j in 1:ncol(aa)){
    aa_modif_all_comb_p20_1398[[i]] <- c(aa_modif_all_comb_p20_1398[[i]], aa_modif_all_options_p20_1398[[j]][aa[i, j]])
  }
  aa_modif_all_comb_p20_1398[[i]] <- unlist(aa_modif_all_comb_p20_1398[[i]])
}

# Following is the index of the chosen models
aaaa <-sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All, decreasing = T, index.return=T)$ix[1:25]
aa_model_index_1398 <- aaaa
aa_model_index_1398_ERBS
aa_model_index_1398_all <- union(aa_model_index_1398, aa_model_index_1398_ERBS)

# updating the data file
aa_name_conv <- names(TF.motifs.Shrinked.hocomoco.count[abs(sort(unlist(aa_modif_all_options_p20_1398)))[1:12]])
aa_name_conv <- c(aa_name_conv, names(TF.motifs.Shrinked.hocomoco.count[abs(sort(unlist(aa_modif_all_options_p20_1398)))[13:18] - 20]))
aa_name_conv[1:6] <- paste(aa_name_conv[1:6],"n",sep = "_")
aa_name_conv[7:12] <- paste(aa_name_conv[7:12],"p",sep = "_")
aa_name_conv[13:18] <- paste(aa_name_conv[13:18],"all",sep = "_")
aa_name_conv <- cbind(sort(unlist(aa_modif_all_options_p20_1398)), aa_name_conv)
for(i in 1:length(aa_modif_all_comb_p20_1398)){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[1399,])
  GEMSTAT_experiment_info_df[1399+i, 1] <- paste0("experiment_", 1399+i)
  GEMSTAT_experiment_info_df[1399+i, 8] <- paste("Concentration Modification. Base experiment 1398.  my_model_index <- c(21,56,61,142,49,57,115,9,134,117,127,94,55,118,72,1,30,51,77,113,13,27,34,84,74,3,11,15,16,20,33,42,43,44,47,64,80,92,93,101,108,120,121,130,132,137).", "conc modif: ",
                                                paste(aa_name_conv[aa_name_conv[, 1] %in% as.character(aa_modif_all_comb_p20_1398[[i]]), 2], collapse = " "))
  
}

# Creating the R script to create the inputs and job files.
setwd("Seeded_GEMSTAT_ens/")
aa <- list()

for(i in 1:length(aa_modif_all_comb_p20_1398)){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp1399.R")
#  aa[[i]][16] <- "my_model_index <-  c(21,56,61,142,49,57,115,9,134,117,127,94,55,118,72,1,30,51,77,113,13,27,34,84,74,3,11,15,16,20,33,42,43,44,47,64,80,92,93,101,108,120,121,130,132,137)"
  aa[[i]][19] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_ER01_monomer_hocomoco"
  
  for(j in 1:length(aa_modif_all_comb_p20_1398[[i]])){
    if(aa_modif_all_comb_p20_1398[[i]][j] < 0){
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[',-aa_modif_all_comb_p20_1398[[i]][j],', ] <- rep(c(1,0), 45)')
    }else if(aa_modif_all_comb_p20_1398[[i]][j] < 21){
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[',aa_modif_all_comb_p20_1398[[i]][j],', ] <- rep(c(0,1), 45)')
    }else{
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[',aa_modif_all_comb_p20_1398[[i]][j] - 20,', ] <- rep(1, 90)')
    }
  }
  
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 1399 + i, ",")
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][91] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F, create_annotation = F"
  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=1399,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=1399,")
  aa[[i]] <- c(aa[[i]], "                              coop=1399, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=1399,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=1399,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=1399, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=1399,")
  aa[[i]] <- c(aa[[i]], "                              annotations=1399,")
  aa[[i]] <- c(aa[[i]], "                              Params=1399, ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=F,")
  aa[[i]] <- c(aa[[i]], "                              seq=1399,")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=1399,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=1399, ")
  aa[[i]] <- c(aa[[i]], "                              wights=1399, ")
  aa[[i]] <- c(aa[[i]], "                              jobfile=1399, ")
  aa[[i]] <- c(aa[[i]], "                              .prev_na=5, ")
  aa[[i]] <- c(aa[[i]], "                              .optim=T, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 1399+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  # 
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  1399 + i, ".R"))
}
setwd("..")
bash_batch_run_zip_transfer_process(exp_index = c(1400:1599))
bash_batch_run_zip_transfer_process(exp_index = c(1600:1799))
bash_batch_run_zip_transfer_process(exp_index = c(1800:1999))
bash_batch_run_zip_transfer_process(exp_index = c(2000:2128))




for(i in 1400:2128){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
  # 
  print("####################################################################################################")
  print("####################################################################################################")
  print("####################################################################################################")
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print("####################################################################################################")
  print("####################################################################################################")
  print("####################################################################################################")
  # # Perform TF KD
  # KD_Performer(parent_Exp=i, 
  #              TF_index=list(3),
  #              shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
  #              prev_na=5,
  #              run_GEMSTAT = F)
  GEMSTAT_based_on_linear_exp_KD_Results[[i]] <- read_KD_results(experiment_nu = i,
                                                               shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens")
  
  # Perform coop KD
  # KD_Performer_by_params(parent_Exp=i,
  #                        TF_names = names(TF.motifs.Shrinked.hocomoco.count),
  #                        Coop_mat_list=list(c(3,3)),
  #                        shared_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
  #                        prev_na=5,
  #                        run_GEMSTAT=F)
  GEMSTAT_based_on_linear_exp_KD_COOP_Results[[i]] <-  read_KD_results(experiment_nu=i,
                                                                       shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                                                                       coop_KD=T)
  
  
}
par(mfrow = c(1, 1), mar = c(6,4,4,4))
#c(c(1400:1900), c(2000:2128))
aa <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(1400:2128)], "[[", 7))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All[aa_model_index_1398_all], aa), las= 2, main="Receptor activation profile performance")
abline(h=seq(0, 1, 0.01), col = 2, lty = 4)
heatmap.2(aa_modif_all_comb_p20_raw_1398,
          Rowv = F, Colv = F, dendrogram = "none", trace = 'none')
aa_dat <- data.frame( x=as.numeric(aa), above_60= as.numeric(aa) >= 0.6 )
#library(ggplot2)
qplot(x,data=aa_dat,geom="histogram",fill=above_60, cex = 3)

aaaa <- apply(X = aa, MARGIN = 2,which.max)
aaaa_max <- apply(X = aa, MARGIN = 2,max)
aaaa_max_ind <- sort(aaaa_max, decreasing = T, index.return = T)$ix
#aa_modif_all_comb_p20_mat <- do.call(what = rbind, aa_modif_all_comb_p20)
#aa_modif_all_comb_p20_raw_1398 <- aa_modif_all_comb_p20_raw_1398[c(1400:2128),]
aa_modif_all_comb_p20_raw_sort_1398 <- aa_modif_all_comb_p20_raw_1398[aaaa_max_ind,]
#aa_modif_all_comb_p20_raw_sort <- aa_modif_all_comb_p20_raw_sort[, -3]
colnames(aa_modif_all_comb_p20_raw_sort_1398) <- c("AR", "GR", "PGR", "RARA", "RARG", "RXRA")

heatmap.2(aa_modif_all_comb_p20_raw_sort_1398, 
          Rowv = F, Colv = F, dendrogram = "none",
          trace = 'none', main = "RAP sorted by Performance")

aa_top <- aa_modif_all_comb_p20_raw_sort_1398[1:60,]
aa_top_table <- numeric(6*3)
for(i in 1:ncol(aa_top)){
  aa_tab <- table(aa_top[, i])
  for(j in 1:length(aa_tab)){
    aa_top_table[3 * (i - 1) + as.numeric(names(aa_tab)[j])] <- aa_tab[j]
  }
}

barplot(aa_top_table, 
        col=c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3), rep(6, 3)), 
        names.arg = c("ARp", "ARn", "ARall", "GRp", "GRn","GRall","PGRp", "PGRn",
                      "PGRall","RARAp", "RARAn","RARAall","RARGp", "RARGn","RARGall",
                      "RXRAp", "RXRAn", "RXRAall"), las=2, main = "top60_experiments")
### Filter for experiments that are ER sensitive

aa_WT_acc_1398 <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[c(1400:2128)], "[[", 7))
aa_ER_KD_Acc_1398 <- do.call(cbind, lapply(lapply(X = GEMSTAT_based_on_linear_exp_KD_Results[c(1400:2128)], "[[", 1), "[[", 7))
aa_coop_ER_kd_acc_1398 <- do.call(cbind, lapply(lapply(X = GEMSTAT_based_on_linear_exp_KD_COOP_Results[c(1400:2128)], "[[", 1), "[[", 7))
colnames(aa_coop_ER_kd_acc_1398) <- colnames(aa_ER_KD_Acc_1398)
aa_acc_diff_1398 <- aa_WT_acc_1398 - aa_ER_KD_Acc_1398

aa_dat <- data.frame( x=as.numeric(aa_acc_diff_1398), above_10= as.numeric(aa_acc_diff_1398) >= 0.1 )
#library(ggplot2)
qplot(x,data=aa_dat,geom="histogram",
      fill=above_10, cex = 3, main = "WT_Accuracy - ER_KD_Accuracy")


aa_acc_diff_6_1398 <- aa_acc_diff_1398
aa_acc_diff_6_1398[aa_WT_acc_1398 < 0.6] <- NA
boxplot.matrix(aa_acc_diff_6_1398, las = 2)
sum(aa_acc_diff_6_1398>= 0.06, na.rm =T) 
colnames(aa_acc_diff_6_1398) <- c(1400:2128)
colnames(aa_WT_acc_1398) <- c(1400:2128)
colnames(aa_ER_KD_Acc_1398) <- c(1400:2128)

aa_acc_diff_6_10_1398 <- aa_acc_diff_6_1398
aa_acc_diff_6_6_1398 <- aa_acc_diff_6_1398

aa_acc_diff_6_10_1398[aa_acc_diff_6_10_1398 < 0.1] <- NA
aa_acc_diff_6_6_1398[aa_acc_diff_6_6_1398 < 0.06] <- NA

aaKD_coop_dif <-aa_ER_KD_Acc_1398[!is.na(aa_acc_diff_6_10_1398)] -aa_coop_ER_kd_acc_1398[!is.na(aa_acc_diff_6_10_1398)]
par(mfrow = c(1, 1), mar = rep(6, 4))
hist(aaKD_coop_dif, main = "ERKD_acc - ERCOOPKD_acc", breaks = 10)
aa_acc_diff_6_10_f_1398 <- aa_acc_diff_6_10_1398[,! colSums(is.na(aa_acc_diff_6_10_1398)) == nrow(aa_acc_diff_6_10_1398)]
aa_acc_diff_6_6_f_1398  <- aa_acc_diff_6_6_1398[,! colSums(is.na(aa_acc_diff_6_6_1398)) == nrow(aa_acc_diff_6_6_1398)]

GEMSTAT_experiment_info_df[as.numeric(colnames(aa_acc_diff_6_10_f_1398)) , 8]
colnames(aa_modif_all_comb_p20_raw_1398)<-c("AR", "GR", "PGR", 
                                            "RARA", "RARG", "RXRA")


heatmap.2(aa_modif_all_comb_p20_raw_1398[colnames(aa_acc_diff_6_10_f_1398),]
          , Rowv = F, Colv = F, dendrogram = "none",
          trace = 'none', distfun = my_dist_fuc,
          main = "RAP acc > 0.6 and ER-sensitive", col = c("green", "Red", "White"))

heatmap.2(aa_modif_all_comb_p20_raw_1398[colnames(aa_acc_diff_6_6_f_1398),]
          , Rowv = T, Colv = T, dendrogram = "both",
          trace = 'none', distfun = my_dist_fuc, col = c("green", "Red", "White"))

aa_top <- aa_modif_all_comb_p20_raw_1398[colnames(aa_acc_diff_6_10_f_1398),]
aa_top_table <- numeric(6*3)
for(i in 1:ncol(aa_top)){
  aa_tab <- table(aa_top[, i])
  for(j in 1:length(aa_tab)){
    aa_top_table[3 * (i - 1) + as.numeric(names(aa_tab)[j])] <- aa_tab[j]
  }
}
par(mfrow = c(1, 1), mar = rep(6, 4))
barplot(aa_top_table, 
        col=c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3), rep(6, 3)), 
        names.arg = c("ARp", "ARn", "ARall", "GRp", "GRn","GRall","PGRp", "PGRn",
                      "PGRall","RARAp", "RARAn","RARAall","RARGp", "RARGn","RARGall",
                      "RXRAp", "RXRAn", "RXRAall"), las=2, main = "WT acc > 0.6, ER KD drop > 0.1")

aaaaaa <- aa_WT_acc_1398[, colnames(aa_acc_diff_6_10_f_1398)]
par(mfrow = c(6, 4), mar = c(2, 3, 1, 1))
for(i in 1:ncol(aa_acc_diff_6_10_f_1398)){
  plot(aa_WT_acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[i]], 
       aa_ER_KD_Acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[i]],
       xlim = c(0.3, 0.7), ylim = c(0.3, 0.7), ylab = "KD", col = 2)
  abline(c(0, 1), col=1)
  points(aa_WT_acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[i]],
         aa_coop_ER_kd_acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[i]], pch = 4, col = 3)
}

# perform Other KnockDowns for these models
aa_GEMSTAT_based_on_linear_exp_KD_Results <- list()
aa_GEMSTAT_based_on_linear_exp_KD_COOP_Results <- list()
for(i in as.integer(colnames(aa_acc_diff_6_10_f_1398))){
  print("####################################################################################################")
  print("####################################################################################################")
  print("####################################################################################################")
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print("####################################################################################################")
  print("####################################################################################################")
  print("####################################################################################################")
  # Perform TF KD
  KD_Performer(parent_Exp=i,
               TF_index=as.list(c(c(1:2), c(4:18))),
               shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
               prev_na=5,
               run_GEMSTAT = T)
  aa_GEMSTAT_based_on_linear_exp_KD_Results[[i]] <- read_KD_results(experiment_nu = i,
                                                                 shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens")
  
  # Perform coop KD
  KD_Performer_by_params(parent_Exp=i,
                         TF_names = names(TF.motifs.Shrinked.hocomoco.count),
                         Coop_mat_list=list(c(1, 1), c(9,9), c(13, 13), c(14, 14), c(16,16)),
                         shared_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                         prev_na=5,
                         run_GEMSTAT=T)
  aa_GEMSTAT_based_on_linear_exp_KD_COOP_Results[[i]] <-  read_KD_results(experiment_nu=i,
                                                                       shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                                                                       coop_KD=T)
  
  
}
GEMSTAT_based_on_linear_exp_KD_Results[as.integer(colnames(aa_acc_diff_6_10_f_1398))] <- aa_GEMSTAT_based_on_linear_exp_KD_Results[as.integer(colnames(aa_acc_diff_6_10_f_1398))]
GEMSTAT_based_on_linear_exp_KD_COOP_Results[as.integer(colnames(aa_acc_diff_6_10_f_1398))] <- aa_GEMSTAT_based_on_linear_exp_KD_COOP_Results[as.integer(colnames(aa_acc_diff_6_10_f_1398))]




####################################################################################################
# plot the kncokdown data
# boxplot: one box for each TF or COOP KD. each point is a model. only models with performance higher than 0.6 and 
aa_KD_holder <- list()
aa_Coop_holder <- list()
aa_names <- as.integer(names(aa_GEMSTAT_based_on_linear_exp_KD_Results[[2084]]))
for(i in 1:length(aa_GEMSTAT_based_on_linear_exp_KD_Results[[2084]])){
  
  aa_KD_holder[[i]] <- do.call(cbind, lapply(lapply(aa_GEMSTAT_based_on_linear_exp_KD_Results[as.integer(colnames(aa_acc_diff_6_10_f_1398))], "[[", which(aa_names %in% i)), "[[", 7))
}
names(aa_KD_holder) <- names(TF.motifs.Shrinked.hocomoco.count)

for(i in 1:length(aa_GEMSTAT_based_on_linear_exp_KD_COOP_Results[[2084]])){
  aa_Coop_holder[[i]] <- do.call(cbind, lapply(lapply(aa_GEMSTAT_based_on_linear_exp_KD_COOP_Results[as.integer(colnames(aa_acc_diff_6_10_f_1398))], "[[", i), "[[", 7))
}
names(aa_Coop_holder) <- names(aa_GEMSTAT_based_on_linear_exp_KD_COOP_Results[[2084]])


aa_KD_eff_mat <- matrix(nrow = sum(!is.na(aa_acc_diff_6_10_f_1398)),
                        ncol = (length(aa_KD_holder)+length(aa_Coop_holder)))

colnames(aa_KD_eff_mat) <- c(names(aa_KD_holder), names(aa_Coop_holder))
for(i in 1:length(aa_KD_holder)){
  aa_KD_eff_mat[, i] <- aa_KD_holder[[i]][!is.na(aa_acc_diff_6_10_f_1398)]
}
for(i in 1:length(aa_Coop_holder)){
  aa_KD_eff_mat[, length(aa_KD_holder) + i] <- aa_Coop_holder[[i]][!is.na(aa_acc_diff_6_10_f_1398)]
}
boxplot.matrix(aa_KD_eff_mat, 
               las = 2,
               main = "KD Accuracy of Filtered models")

# difference with WT

aa_KD_diff_mat <- matrix(nrow = sum(!is.na(aa_acc_diff_6_10_f_1398)),
                        ncol = (length(aa_KD_holder)+length(aa_Coop_holder)))

colnames(aa_KD_diff_mat) <- c(names(aa_KD_holder), names(aa_Coop_holder))
for(i in 1:length(aa_KD_holder)){
  aa_KD_diff_mat[, i] <- aaaaaa[!is.na(aa_acc_diff_6_10_f_1398)] -
    aa_KD_holder[[i]][!is.na(aa_acc_diff_6_10_f_1398)]
}


aaaaaa <- aa_WT_acc_1398[, colnames(aa_acc_diff_6_10_f_1398)]
for(i in 1:length(aa_Coop_holder)){
  aa_KD_diff_mat[, length(aa_KD_holder) + i] <- aaaaaa[!is.na(aa_acc_diff_6_10_f_1398)] -  aa_Coop_holder[[i]][!is.na(aa_acc_diff_6_10_f_1398)]
}
boxplot.matrix(aa_KD_diff_mat, 
               las = 2,
               main = "Acc_wt - Acc_KD  of Filtered models")
abline(h=seq(-0.05, 0.2, 0.01), col = 2, lwd = 0.4, lty = 4)



aa_KD_eff_mat_rank <- aa_KD_eff_mat
for(i in 1:nrow(aa_KD_eff_mat)){
  aa_KD_eff_mat_rank[i, ] <- rank(aa_KD_eff_mat[i, ])
}

summary
boxplot.matrix(aa_KD_eff_mat_rank,
               las = 2,
               main = "KD effect Rank of Filtered models")
abline(h=seq(0, 24, 1), col = 2, lwd = 0.5, lty = 4)



aaKD_coop_dif <-aa_ER_KD_Acc_1398[!is.na(aa_acc_diff_6_10_1398)] -aa_coop_ER_kd_acc_1398[!is.na(aa_acc_diff_6_10_1398)]


# for each receptor activation profile
####################################################################################################
# what are the parameters for the models 
aa_mod_ind <- list()
for(i in 1:ncol(aa_acc_diff_6_10_f_1398)){
  aa_mod_ind[[i]] <- which(! is.na(aa_acc_diff_6_10_f_1398[, i]))
}
names(aa_mod_ind) <- colnames(aa_acc_diff_6_10_f_1398)
aa_base_linear_model <- aa_model_index_1398_all[unlist(aa_mod_ind)]

# look at the parameters
aa_bindings <- matrix(nrow = length(unlist(aa_mod_ind)), ncol = 18)
colnames(aa_bindings) <- names(GEMSTAT_based_on_linear_exp_pars[[1400]]$binding[8,])
rownames(aa_bindings) <- character(nrow(aa_bindings))

aa_alphas <- matrix(nrow = length(unlist(aa_mod_ind)), ncol =18)
colnames(aa_alphas) <- names(GEMSTAT_based_on_linear_exp_pars[[1400]]$alpha_effective[8,])
rownames(aa_alphas) <- character(nrow(aa_alphas))

aa_coops <- matrix(nrow = length(unlist(aa_mod_ind)), ncol =6)
colnames(aa_coops) <- names(GEMSTAT_based_on_linear_exp_pars[[1400]]$coop[[8]])
rownames(aa_coops) <- character(nrow(aa_coops))

aa_cnt = 1
for(i in 1:length(aa_mod_ind)){
  for(j in 1:length(aa_mod_ind[[i]])){
    aa_bindings[aa_cnt, ] <- GEMSTAT_based_on_linear_exp_pars[[as.numeric(names(aa_mod_ind))[i]]]$binding[aa_mod_ind[[i]][j], ]
    aa_alphas[aa_cnt, ] <- GEMSTAT_based_on_linear_exp_pars[[as.numeric(names(aa_mod_ind))[i]]]$alpha_effective[aa_mod_ind[[i]][j], ]
    aa_coops[aa_cnt, ] <- GEMSTAT_based_on_linear_exp_pars[[as.numeric(names(aa_mod_ind))[i]]]$coop[[aa_mod_ind[[i]][j]]]
    rownames(aa_bindings)[aa_cnt] <- paste0(as.numeric(names(aa_mod_ind))[i], "_", aa_mod_ind[[i]][j])
    rownames(aa_alphas)[aa_cnt] <- paste0(as.numeric(names(aa_mod_ind))[i], "_", aa_mod_ind[[i]][j])
    rownames(aa_coops)[aa_cnt] <- paste0(as.numeric(names(aa_mod_ind))[i], "_", aa_mod_ind[[i]][j])
    aa_cnt <- aa_cnt + 1
  }
  
}
par(mfrow = c(1, 1), mar = c(6,4,4,4))
boxplot.matrix(aa_bindings, las = 2)
boxplot.matrix(aa_alphas, las = 2)
boxplot.matrix(aa_coops, las = 2)

aa_exp_names <- unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 1))
aa_col <- sample(x = col_vector, size = 5, replace = F)
aa_raw_RAP_comb_chosen <- aa_modif_all_comb_p20_raw_1398[aa_exp_names,]

heatmap.2(cbind(aa_alphas, aa_modif_all_comb_p20_raw_1398[aa_exp_names,])
          , Rowv = T, Colv = F, dendrogram = "row",
          trace = 'none', breaks = c(0, 0.5, 1, 2.5, 3, 7, 10),
          margins = c(8,8), col = col_vector[11:16], main = "alphaERsensitive+concComb_1398", 
          RowSideColors = col_vector[(aa_base_linear_model + 10)%%length(col_vector)])

heatmap.2(aa_coops
          , Rowv = T, Colv = T, dendrogram = "none",
          trace = 'none', breaks = c(100000, 300000, 600000, 900000, 1000000),
          margins = c(8,8), col = col_vector[13:16], main="CoopERsensitive", 
          RowSideColors = col_vector[(aa_base_linear_model + 10)%%length(col_vector)])

heatmap.2(cbind(log10(aa_coops), aa_modif_all_comb_p20_raw_1398[aa_exp_names,])
          , Rowv = F, Colv = F, dendrogram = "none",
          trace = 'none', breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5,  5.2, 5.4, 5.6, 5.8, 6),
          margins = c(8,8), main="CoopERsensitive", col = col_vector[51:60],
          RowSideColors = col_vector[(aa_base_linear_model + 10)%%length(col_vector)])


heatmap.2(cbind(aa_bindings, aa_modif_all_comb_p20_raw_1398[aa_exp_names,])
          , Rowv = T, Colv = F, dendrogram = "none",
          trace = 'none', breaks = c(0, 1, 2, 3, 10, 50),
          col = col_vector[11:15], margins = c(8,8), 
          main="BindingERsensitive")

# pair the parameters and concentration profiles to see how does that effect the parameters


par(mfrow = c(4, 4), mar = c(2, 3, 4, 1))
plot(aa_WT_acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[1]], 
     aa_ER_KD_Acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[1]],
     xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", pch=1, col=2)
# points(aa_WT_acc[, colnames(aa_acc_diff_6_10_f)[1]],
#        aa_PGR_KD_Acc[, colnames(aa_acc_diff_6_10_f)[1]],
#        col = 4 , pch = 2)
abline(c(0, 1), col=1, lwd = 0.6)
# par(xpd=T)
# legend(legend = c("ER_KD", "PGR_KD"),
#        fill = c(2, 4), x =0.4, y=1
#        #, pch = c(19, 16)
#        ,bty = "n", y.intersp = 0.6, cex = 0.7)
# par(xpd=F)



for(i in 2:ncol(aa_acc_diff_6_10_f_1398)){
  plot(aa_WT_acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[i]], 
       aa_ER_KD_Acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[i]],
       xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", pch=1, col=2)

  
  abline(c(0, 1), col=1, lwd = 0.6)
}


par(mfrow = c(4, 4), mar = c(2, 3, 5, 1))
plot(aa_WT_acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[1]],
     aa_ER_KD_Acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[1]],
     xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD",
     col=2, pch=1)

par(xpd=F)
for(i in 2:ncol(aa_acc_diff_6_10_f_1398)){
  plot(aa_WT_acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[i]],
       aa_ER_KD_Acc_1398[, colnames(aa_acc_diff_6_10_f_1398)[i]],
       xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", col=2, pch=1)

  abline(c(0, 1), col=1, lwd = 0.6)
}

##########################################
# create model summary:
# 1. linear base model number
# 2. % genes with >= 1 ERBS
# 3. a Table with one row per TF and the following columns:
#  KD rank, RAP, Alpha, Binding, selfCoop
aa_base_linear_model <- aa_model_index_1398_all[unlist(aa_mod_ind)]
aa_nz_all_mat <- matrix(nrow = 148, ncol = 18)
for(i in 1:148){
  for(j in 1:18){
    aa_nz_all_mat[i,j] <- sum(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat[[i]][, j] > 0)/
      nrow(GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat[[i]])
    
  }
}
colnames(aa_nz_all_mat) <- names(TF.motifs.Shrinked.hocomoco.count)
par(mfrow = c(1,1), mar = c(6,6,6,6))
boxplot.matrix(aa_nz_all_mat, las =2, main = "%Genes containing at least one BS")
#aa_nz_all[aa_base_linear_model]
aa_GEMSTAT_filtered_summary <- list()
for(i in 1:length(aa_base_linear_model)){
  aa_GEMSTAT_filtered_summary[[i]] <- list()
  aa_GEMSTAT_filtered_summary[[i]][[1]] <- aa_base_linear_model[i]
  aa_GEMSTAT_filtered_summary[[i]][[2]] <- aa_nz_all_mat[aa_GEMSTAT_filtered_summary[[i]][[1]], ]
  names(aa_GEMSTAT_filtered_summary[[i]][[2]]) <- names(TF.motifs.Shrinked.hocomoco.count)
  aa_GEMSTAT_filtered_summary[[i]][[3]] <- matrix(nrow = ncol(aa_alphas), ncol =5)
  rownames(aa_GEMSTAT_filtered_summary[[i]][[3]]) <- names(TF.motifs.Shrinked.hocomoco.count)
  colnames(aa_GEMSTAT_filtered_summary[[i]][[3]]) <- c("KD_rank", "RAP", "Alpha", "Binding", "SelfCoop")
  names(aa_GEMSTAT_filtered_summary[[i]]) <- c("Base_Model", "Percent_ERBS", "TF_mat")
  
  for(j in 1:nrow(aa_GEMSTAT_filtered_summary[[i]][[3]])){
    # setting rank
    aa_GEMSTAT_filtered_summary[[i]][[3]][j, 1] <- aa_KD_eff_mat_rank[i, j]
    # setting RAP
    if(names(TF.motifs.Shrinked.hocomoco.count)[j] %in% colnames(aa_raw_RAP_comb_chosen)){
      aa_GEMSTAT_filtered_summary[[i]][[3]][j, 2] <- aa_raw_RAP_comb_chosen[i, names(TF.motifs.Shrinked.hocomoco.count)[j]]
    }
    # set alpha
    aa_GEMSTAT_filtered_summary[[i]][[3]][j, 3] <- aa_alphas[i, j]
    
    # set binding 
    aa_GEMSTAT_filtered_summary[[i]][[3]][j, 4] <- aa_bindings[i, j]
    
    # set coop 
    if( names(TF.motifs.Shrinked.hocomoco.count)[j] %in% unlist(lapply(strsplit(colnames(aa_coops), split = ":"), "[[", 1))){
      aa_GEMSTAT_filtered_summary[[i]][[3]][j, 5] <- aa_coops[i, which(unlist(lapply(strsplit(colnames(aa_coops), split = ":"), "[[", 1)) %in% names(TF.motifs.Shrinked.hocomoco.count)[j] )]
    }
  }
}
names(aa_GEMSTAT_filtered_summary) <- paste(colnames(aa_acc_diff_6_10_f_1398),
                                            unlist(aa_mod_ind), sep="_")
aa_GEMSTAT_filtered_summary$`1434_34`$TF_mat

aa_GEMSTAT_filtered_summary_per_TF <- list()
for(i in 1:nrow(aa_GEMSTAT_filtered_summary[[1]]$TF_mat)){
  aa_GEMSTAT_filtered_summary_per_TF[[i]] <- matrix(nrow = length(aa_GEMSTAT_filtered_summary),
                                                    ncol = ncol(aa_GEMSTAT_filtered_summary[[1]]$TF_mat) + 1)
  rownames(aa_GEMSTAT_filtered_summary_per_TF[[i]]) <- names(aa_GEMSTAT_filtered_summary)
  colnames(aa_GEMSTAT_filtered_summary_per_TF[[i]]) <- c(colnames(aa_GEMSTAT_filtered_summary[[1]]$TF_mat), "Percent_gene_with_site")
  for(j in 1:length(aa_GEMSTAT_filtered_summary)){
    aa_GEMSTAT_filtered_summary_per_TF[[i]][j, 1:ncol(aa_GEMSTAT_filtered_summary[[j]]$TF_mat)] <- aa_GEMSTAT_filtered_summary[[j]]$TF_mat[i, ]
    aa_GEMSTAT_filtered_summary_per_TF[[i]][j, ncol(aa_GEMSTAT_filtered_summary_per_TF[[i]])] <- aa_GEMSTAT_filtered_summary[[j]][[2]][i]
  }
  
}
names(aa_GEMSTAT_filtered_summary_per_TF) <- names(TF.motifs.Shrinked.hocomoco.count)
aa_GEMSTAT_filtered_summary_per_TF$PGR

# plot per RAP
for(i in c(1, 9, 12, 13, 14, 16)){
  par(mfrow = c(6,3), mar = c(3,2,1,1))
  for(j in 1:18){
    boxplot(aa_GEMSTAT_filtered_summary_per_TF[[j]][, 6]~aa_GEMSTAT_filtered_summary_per_TF[[i]][, 2],
            main=paste(names(TF.motifs.Shrinked.hocomoco.count)[j], names(TF.motifs.Shrinked.hocomoco.count)[i], sep = "_"))
    abline(h = 0, col = 2, lwd= 1, lty = 3)
  }
}

for(i in c(1, 9, 12, 13, 14, 16)){
  par(mfrow = c(2,3), mar = c(3,2,1,1))
  for(j in c(1, 3, 9, 13, 14, 16)){
    boxplot(log10(aa_GEMSTAT_filtered_summary_per_TF[[j]][, 5])~aa_GEMSTAT_filtered_summary_per_TF[[i]][, 2],
            main=paste(names(TF.motifs.Shrinked.hocomoco.count)[j], names(TF.motifs.Shrinked.hocomoco.count)[i], sep = "_"))
    abline(h = 0, col = 2, lwd= 1, lty = 3)
  }
}

for(i in c(1, 9, 12, 13, 14, 16)){
  par(mfrow = c(2,3), mar = c(3,2,1,1))
  for(j in c(1, 9, 12, 13, 14, 16)){
    boxplot(aa_GEMSTAT_filtered_summary_per_TF[[j]][, 2]~aa_GEMSTAT_filtered_summary_per_TF[[i]][, 2],
            main=paste(names(TF.motifs.Shrinked.hocomoco.count)[j], names(TF.motifs.Shrinked.hocomoco.count)[i], sep = "_"))
   # abline(h = 0, col = 2, lwd= 1, lty = 3)
  }
}


aa_cov <- cor(aa_raw_RAP_comb_chosen)
heatmap.2(aa_cov, Rowv = F, Colv = F, dendrogram = "none", trace = "none", 
          breaks = c(-1, -0.5, -0.2, 0.2, 0.5, 1),
          col = col_vector[11:15])


library(ggm)
aa_pcor <- pcor.test(aa_raw_RAP_comb_chosen[, 1], aa_raw_RAP_comb_chosen[, 2], aa_raw_RAP_comb_chosen[, 3:6])

aa_pcor_all <- matrix(nrow = ncol(aa_raw_RAP_comb_chosen), 
                      ncol = ncol(aa_raw_RAP_comb_chosen))
for(i in 1:(ncol(aa_raw_RAP_comb_chosen) - 1)){
  for(j in (i+1):ncol(aa_raw_RAP_comb_chosen)){
    aa_cnd <- colnames(aa_raw_RAP_comb_chosen)[c(c(i , j), setdiff(c(1:ncol(aa_raw_RAP_comb_chosen)), c(i, j)))]
    aa_pcor_all[i, j] <- pcor(aa_cnd, var(aa_raw_RAP_comb_chosen))
    aa_pcor_all[j, i] <- aa_pcor_all[i, j]
  }
}
colnames(aa_pcor_all) <- colnames(aa_raw_RAP_comb_chosen)
rownames(aa_pcor_all) <- colnames(aa_raw_RAP_comb_chosen)
heatmap.2(aa_pcor_all, Rowv = F, Colv = F, dendrogram = "none", trace = "none"
#          ,  breaks = c(-1, -0.5, -0.2, 0.2, 0.5, 1),
 #         col = col_vector[11:15]
)

########################################################################################################################
aa_exp_names_2 <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 1))) - 1399
aa_model_names <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 2)))
aa_model_names_org <- aa_model_index_1398_all[aa_model_names]

colnames(aa_modif_all_comb_p20_raw_1398)[2] <- "NR3C1"

GeneFoldChange_by_TFFoldChange_model_RAP_wrapper(.Gene_exp_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_48_gene_gte3WTGRcommon,
                                                 .TF_exp_mat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_ER01_monomer_hocomoco,
                                                 TF_gene_site_mat_list = GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat, 
                                                 .site_nu_thresh = 1,
                                                 # .TF_index = c(1:nrow(.TF_exp_mat)),
                                                 RAP_raw_mat=aa_modif_all_comb_p20_raw_1398, 
                                                 RAP_list = aa_modif_all_options_p20_1398,
                                                 RAP_nu=aa_exp_names_2[13],
                                                 model_nu = aa_model_names_org[13])
########################################################################################################################
aa_WT_Pred_mat_list <- list()
aa_KD_Pred_mat_list <- list()
aa_exp_names <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 1)))
aa_model_names <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 2)))
aa_model_names_org <- aa_model_index_1398_all[aa_model_names]
par(mfrow = set_row_col_nu(length(TF.motifs.Shrinked.hocomoco.count)),
    mar= c(4, 3, 3, 3))
aa <- list()
for(j in 1:length(TF.motifs.Shrinked.hocomoco.count)){
  aa_TF_index <- j
  for(i in 1:length(aa_exp_names)){
    aa_WT_Pred_mat_list[[i]] <- GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[aa_exp_names[i]]]$Prediction_round_list[[aa_model_names[i]]]
    aa_KD_Pred_mat_list[[i]] <- GEMSTAT_based_on_linear_exp_KD_Results[[aa_exp_names[i]]][[as.character(aa_TF_index)]]$Prediction_round_list[[aa_model_names[i]]]
  }
  aa[[j]] <- KD_by_affected_enhancers(WT_Pred_mat_list = aa_WT_Pred_mat_list, 
                                 KD_Pred_mat_list = aa_KD_Pred_mat_list,
                                 real_mat = my_CommonDifExpMat_16_ERassoc_gte4nzKD_atl1p1n_48_gene_gte3WTGRcommon,
                                 TF_gene_site_mat_list = GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3$Site_count_mat,
                                 site_nu_thresh=1,
                                 model_index = aa_model_names_org, 
                                 TF_index=aa_TF_index,plot_sc = F)
  plot(aa[[j]][[1]], aa[[j]][[3]], 
       main = names(TF.motifs.Shrinked.hocomoco.count)[j], 
       ylab = "KD effect", xlab = "#targetted genes")
  abline(h = seq(-0.1, 0.5, 0.01), lty = 3, lwd = 0.6, col = 3)
}

########################################################################################################################################################################################
##########################################################################################################################################
# Experiment 2129
# new set of genes (172), based on linear models trained only on 62 of the 172 genes. filtered for test performance above 0.5 --> 78 models survived
# added a E2F1 TF to the model. currently uses 19 TFs

GEMSTAT_based_on_linear_exp_annotation_172_78models <- list()
for(i in 1:78){
  print(i)
  GEMSTAT_based_on_linear_exp_annotation_172_78models[[i]] <- annotation_reader(annot_file = paste0("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_2/Annotation/Experiment_2_", i, ".annot"), 
                                                                        TF_names = names(TF.motifs.Shrinked.hocomoco_E2F1_added.count))
  
}
TF.motifs.Shrinked.hocomoco_E2F1_added.count_MaxLLR <- c(6.13309, 10.0693, 11.8223, 6.28522, 11.3454, 7.91607, 9.98706,
                                              8.40569, 12.5907, 7.2447, 10.2185, 13.5565, 11.9268,
                                              6.69131, 7.37613, 9.94809, 6.93567, 12.5285, 9.73563)
names(TF.motifs.Shrinked.hocomoco_E2F1_added.count_MaxLLR) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)

aa_homodimer <-   c(T, F, F, T, F, F, F, F, F, T, F, F, F, T, T, F, T, F, F)
names(aa_homodimer) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
aa_dimer_orien <-        c(2, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 0, 1, 0, 0)
names(aa_dimer_orien) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
aa_homoDimer_distance <- c(3, 0, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 5, 5, 0, 1, 0, 0)
names(aa_homoDimer_distance) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)

aa <- count_site_from_annotation_combined(annotation_list = GEMSTAT_based_on_linear_exp_annotation_172_78models, 
                                          TF_names = names(TF.motifs.Shrinked.hocomoco_E2F1_added.count),
                                          TF_index=c(1:19),
                                          homoDimer = aa_homodimer, 
                                          dimer_orientation = aa_dimer_orien,
                                          homoDimer_distance=aa_homoDimer_distance,
                                          heterodimer_pair = numeric(0),
                                          heterodimer_distance = numeric(0),
                                          MAXLLR = TF.motifs.Shrinked.hocomoco_E2F1_added.count_MaxLLR,
                                          annotation_thresh = rep(1, 19),
                                          LLR_to_pVal_list = TF.motifs.Shrinked.hocomoco.count_dimer_concat_motif_pval_E2F1_added,
                                          pVal_thresh = 0.001)
GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3 <- aa
aa_ER_S <- matrix(nrow = 78, ncol = 172)
aa_ER_S_nz <- numeric(78)
for(i in 1:nrow(aa_ER_S)){
  aa_ER_S[i, ] <- GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Site_count_mat[[i]][, 4]
  aa_ER_S_nz[i] <- sum(aa_ER_S[i, ] > 0)
}
sum(colSums(aa_ER_S) == 0)
hist(aa_ER_S_nz/172, main = "fraction of genes with at least one ER site")
heatmap.2(aa_ER_S,
          Rowv = F, Colv = F, dendrogram = "none", trace = "none",
          breaks = c(0, 0.5, 1.5, 3, 5, 10, 20), col = col_vector[11:16])

heatmap.2(aa$Site_count_mat[[1]],
          Rowv = F, Colv = F, dendrogram = "none", trace = "none",
          breaks = c(0, 0.5, 1.5, 3, 5, 10, 20), col = col_vector[11:16])
save(list = c("GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3"),
     file = "annotation_manual_172genes.RData")

GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[1397,])
GEMSTAT_experiment_info_df[2129, 1] <- paste0("experiment_", 2129)
GEMSTAT_experiment_info_df[2129, 8] <- "Experiment with New set of Motifs: All from HocoMoco, PGR is not dimer anymore, only one TF for JUN, added E2F1 TF. Gene expression matrix is changed to contain 172 genes. there are 78 linear models"
GEMSTAT_experiment_info_df[2129, 10] <- 172

#GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh <-list()
#GEMSTAT_based_on_linear_exp_pars <- list()
for(i in 2129:2129){
  GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)
  GEMSTAT_based_on_linear_exp_pars[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh)[i] <- paste0("experiment_", i)
  names(GEMSTAT_based_on_linear_exp_pars)[i] <- paste0("experiment_", i)
}

GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Random_Prediction_list[[1]]$Precision
hist(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All, breaks = 20)
max(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All)
sum(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All > 0.5)
aa_dat <- data.frame( x= c(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All,
                           GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Random_Prediction_list[[1]]$Precision), 
                      Shuffled= as.logical(c(rep(0, 78), rep(1, 10))) )
qplot(x,data=aa_dat,geom="histogram",
      fill=Shuffled, main = "GEMSTAT model training accuracy")


KD_Performer(parent_Exp=2129, 
             TF_index=as.list(c(1:19)),
             shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
             prev_na=5,
             run_GEMSTAT = T)
aa_KD_2129 <- list()
aa_KD_2129[[1]] <- read_KD_results(experiment_nu = 2129, 
                              shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                              experiment_output = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]])

aa_KD_coop_2129 <- list()
# Perform coop KD
KD_Performer_by_params(parent_Exp=2129,
                       TF_names = names(TF.motifs.Shrinked.hocomoco_E2F1_added.count),
                       Coop_mat_list=list(c(1, 1), c(4,4),c(10, 10), c(14, 14), c(15, 15), c(17, 17)),
                       shared_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                       prev_na=5,
                       run_GEMSTAT=T)
aa_KD_coop_2129[[1]] <-  read_KD_results(experiment_nu=2129, 
                                    shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens", 
                                    coop_KD=T,
                                    experiment_output = GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]])
# GEMSTAT_based_on_linear_exp_KD_Results <- list()
# GEMSTAT_based_on_linear_exp_KD_COOP_Results <- list()
GEMSTAT_based_on_linear_exp_KD_Results[[2129]] <- aa_KD_2129[[1]]
GEMSTAT_based_on_linear_exp_KD_COOP_Results[[2129]] <- aa_KD_coop_2129[[1]]

par(mfrow = c(1, 1), mar = c(4,4,4,4))
plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All, 
     GEMSTAT_based_on_linear_exp_KD_COOP_Results[[2129]]$ESR1_ESR1$Accuracy_All,
     xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", xlab = "WT", pch = 4, col = 2, cex = 2)
points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All,
       GEMSTAT_based_on_linear_exp_KD_Results[[2129]]$`4`$Accuracy_All, col = 3, pch = 1,cex = 2)

# points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All, 
#        GEMSTAT_based_on_linear_exp_KD_COOP_Results[[2129]]$AR_AR$Accuracy_All,
#        col = 3, pch = 4
       #     ,xlim = c(0.2, 0.7), ylim = c(0.2, 0.7)
# )
# points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[1398]]$Accuracy_All,
#        GEMSTAT_based_on_linear_exp_KD_Results[[1398]]$`1`$Accuracy_All, col = 3, pch = 1)
abline(c(0, 1), col=1, lwd = 0.6)
par(xpd=T)
legend(legend = c("ER_KD", "ER_coop_KD"),
       col = c(3, 2), pch = c(1, 4),x="topleft"
       #, pch = c(19, 16)
       ,bty = "n", y.intersp = 0.6, cex =1)
par(xpd=F)



boxplot.matrix(do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_KD_COOP_Results[[2129]], "[[", 7)), las = 2)
aa_x <- do.call(cbind, lapply(GEMSTAT_based_on_linear_exp_KD_Results[[2129]], "[[", 7))
colnames(aa_x) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[as.integer(names(GEMSTAT_based_on_linear_exp_KD_Results[[2129]]))]
boxplot.matrix(aa_x, las = 2)

par(mfrow = c(4, 5), mar = c(1,1,3,1))
for(i in 1:length(GEMSTAT_based_on_linear_exp_KD_Results[[2129]])){
  plot(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All, 
       GEMSTAT_based_on_linear_exp_KD_Results[[2129]][[as.character(i)]]$Accuracy_All,
       xlim = c(0.2, 0.7), ylim = c(0.2, 0.7), ylab = "KD", 
       xlab = "WT", pch = 1, col = 2,
       main=names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[i], cex = 1.5)
  if(i %in% c(1, 4, 10, 14, 15, 17)){
    points(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All,
           GEMSTAT_based_on_linear_exp_KD_COOP_Results[[2129]][[which(c(1, 4, 10, 14, 15, 17) %in% i)]]$Accuracy_All, col = 3, pch = 4)
  }
  
  abline(c(0, 1), col=1, lwd = 0.6)
}


# plot the difference of predictions of KD models with the WT model. as opossed to comparing their accuracy
aaaa <-sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All, decreasing = T, index.return=T)$ix[1:78]

aa_WT_KD_agreement <- matrix(nrow = length(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All[aaaa]),
                             ncol = length(GEMSTAT_based_on_linear_exp_KD_Results[[2129]]))
for(i in 1:nrow(aa_WT_KD_agreement)){
  aa_wt <- GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Prediction_round_list[[aaaa[i]]]
  for(j in 1:ncol(aa_WT_KD_agreement)){
    aa_k <- GEMSTAT_based_on_linear_exp_KD_Results[[2129]][[as.character(j)]]$Prediction_round_list[[aaaa[i]]]
    aa_WT_KD_agreement[i, j] <- sum(aa_wt == aa_k, na.rm = T)/sum(!is.na(aa_wt))
  }
}
rownames(aa_WT_KD_agreement) <- aaaa
colnames(aa_WT_KD_agreement) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)

par(mfrow = c(1, 1), mar = c(4,4,4,4))
boxplot.matrix(aa_WT_KD_agreement, las = 2, main = "WT_KD_prediction_agreement_all")

########################################################################################################################
#exp 2130
#using these models only: c(52,62,61,22,57,48,19,60,5,37,10,25,29,23,28,70,12,21,53,42,16,9,36,50,4)
# not running, just create inputs so the modification experiments can copy from this
GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[2129,])
GEMSTAT_experiment_info_df[2130, 1] <- paste0("experiment_", 2130)
GEMSTAT_experiment_info_df[2130, 8] <- "same as 2129 only using a subset of the models: c(52,62,61,22,57,48,19,60,5,37,10,25,29,23,28,70,12,21,53,42,16,9,36,50,4)"
GEMSTAT_experiment_info_df[2130, 10] <- 172


########################################################################################################################
# Doing RAP 729 experiments on outputs of best models in experiment 2129
aa_modif_all_options_p20_2129 <- list(list(1, -1, 21),
                                      list(10, -10, 30),
                                      list(13,-13, 33),
                                      list(14,-14, 34),
                                      list(15, -15, 35),
                                      list(17,-17, 37))
aa <- GenerateComb(Set = c(1:6), TotalNumberOfTFs = 6, my_alphabet = c("1", "2", "3"))
aa_modif_all_comb_p20_raw_2129 <- aa

heatmap.2(aa_modif_all_comb_p20_raw_2129)
rownames(aa_modif_all_comb_p20_raw_2129) <- c(2131:2859)
aa_modif_all_comb_p20_2129 <- list()
for(i in 1:nrow(aa)){
  aa_modif_all_comb_p20_2129[[i]] <- integer(0)
  for(j in 1:ncol(aa)){
    aa_modif_all_comb_p20_2129[[i]] <- c(aa_modif_all_comb_p20_2129[[i]], aa_modif_all_options_p20_2129[[j]][aa[i, j]])
  }
  aa_modif_all_comb_p20_2129[[i]] <- unlist(aa_modif_all_comb_p20_2129[[i]])
}

# Following is the index of the chosen models
aaaa <-sort(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All, decreasing = T, index.return=T)$ix[1:25]
#aa_model_index_2129 <- aaaa
aa_model_index_2129_2 <- aaaa
# updating the data file
aa_name_conv <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count[abs(sort(unlist(aa_modif_all_options_p20_2129)))[1:12]])
aa_name_conv <- c(aa_name_conv, names(TF.motifs.Shrinked.hocomoco_E2F1_added.count[abs(sort(unlist(aa_modif_all_options_p20_2129)))[13:18] - 20]))
aa_name_conv[1:6] <- paste(aa_name_conv[1:6],"n",sep = "_")
aa_name_conv[7:12] <- paste(aa_name_conv[7:12],"p",sep = "_")
aa_name_conv[13:18] <- paste(aa_name_conv[13:18],"all",sep = "_")
aa_name_conv <- cbind(sort(unlist(aa_modif_all_options_p20_2129)), aa_name_conv)
for(i in 1:length(aa_modif_all_comb_p20_2129)){
  GEMSTAT_experiment_info_df <- rbind(GEMSTAT_experiment_info_df, GEMSTAT_experiment_info_df[2129,])
  GEMSTAT_experiment_info_df[2130+i, 1] <- paste0("experiment_", 2130+i)
  GEMSTAT_experiment_info_df[2130+i, 8] <- paste("Concentration Modification. Base experiment 2129.  my_model_index <- c(4,5,9,10,12,16,19,21,22,23,25,28,29,36,37,42,48,50,52,53,57,60,61,62,70).", "conc modif: ",
                                                 paste(aa_name_conv[aa_name_conv[, 1] %in% as.character(aa_modif_all_comb_p20_2129[[i]]), 2], collapse = " "))
  
}

# Creating the R script to create the inputs and job files.
setwd("Seeded_GEMSTAT_ens/")
aa <- list()

for(i in 1:length(aa_modif_all_comb_p20_2129)){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp2129.R")
  aa[[i]][15] <- "my_model_index <-  c(52,62,61,22,57,48,19,60,5,37,10,25,29,23,28,70,12,21,53,42,16,9,36,50,4)"
  aa[[i]][19] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01"
  
  for(j in 1:length(aa_modif_all_comb_p20_2129[[i]])){
    if(aa_modif_all_comb_p20_2129[[i]][j] < 0){
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[',-aa_modif_all_comb_p20_2129[[i]][j],', ] <- rep(c(1,0), 45)')
    }else if(aa_modif_all_comb_p20_2129[[i]][j] < 21){
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[', aa_modif_all_comb_p20_2129[[i]][j],', ] <- rep(c(0,1), 45)')
    }else{
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[',aa_modif_all_comb_p20_2129[[i]][j] - 20,', ] <- rep(1, 90)')
    }
  }
  
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 2130 + i, ",")
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][94] <-  "                      ,create_bounds=F, create_coop=F, create_enh_gene_map=F,create_ff=F,create_expmat=F, create_motifs=F, create_params=F,create_seq=F, create_TFexp=T, create_TFinfo=F, create_treatContMap=F,create_weights=F, create_SubmitFile=F, create_jobFile=F, create_annotation = F"
  aa[[i]] <- c(aa[[i]], "                      ,OnlyTFexpModif=T")
  aa[[i]] <- c(aa[[i]], "                      )")
  aa[[i]] <- c(aa[[i]], "GEMSTAT_input_copier_multienh(lower_bounds=2130,")
  aa[[i]] <- c(aa[[i]], "                              upper_bounds=2130,")
  aa[[i]] <- c(aa[[i]], "                              coop=2130, ")
  aa[[i]] <- c(aa[[i]], "                              Enh_gene_map=2130,")
  aa[[i]] <- c(aa[[i]], "                              free_fix=2130,")
  aa[[i]] <- c(aa[[i]], "                              Gene_Exp=2130, ")
  aa[[i]] <- c(aa[[i]], "                              Motifs=2130,")
  aa[[i]] <- c(aa[[i]], "                              annotations=2130,")
  aa[[i]] <- c(aa[[i]], "                              Params=2130, ")
  aa[[i]] <- c(aa[[i]], "                              Params_from_output=F,")
  aa[[i]] <- c(aa[[i]], "                              seq=2130,")
  aa[[i]] <- c(aa[[i]], "                              TFexp=numeric(0),")
  aa[[i]] <- c(aa[[i]], "                              TF_info=2130,")
  aa[[i]] <- c(aa[[i]], "                              Treat_cont_map=2130, ")
  aa[[i]] <- c(aa[[i]], "                              wights=2130, ")
  aa[[i]] <- c(aa[[i]], "                              jobfile=2130, ")
  aa[[i]] <- c(aa[[i]], "                              .prev_na=5, ")
  aa[[i]] <- c(aa[[i]], "                              .optim=T, ")
  aa[[i]] <- c(aa[[i]], paste0("                              cur_exp=", 2130+i, ","))
  aa[[i]] <- c(aa[[i]], "                              shared_dir=\"/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens\")")
  # 
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  2130 + i, ".R"))
}
setwd("..")
bash_batch_run_zip_transfer_process(exp_index = c(2131:2859))
bash_batch_run_zip_transfer_process(exp_index = c(2131:2650))
bash_batch_run_zip_transfer_process(exp_index = c(2651:2859))
# rename the annotation files in all folders
aa_model_index_2129
aa_rm <- setdiff(c(1:78), aa_model_index_2129)
aa_rmch <- paste0("annot_", aa_rm, ".ann")
for(i in 1:729){
  setwd(paste0("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/",
               "Experiment_", (2130 + i), "/Inputs/Filtered_Annotations"))
  file.remove(aa_rmch)
  for(j in 1:length(aa_model_index_2129)){
    file.rename(from = paste0("annot_", aa_model_index_2129[j], ".ann"),
                to   = paste0("annot_", j, ".ann"))
  }
}
setwd("../../..")

#####################
# GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh
# GEMSTAT_based_on_linear_exp_pars
aa_WT2 <- list()
aa_WT2_par <- list()
aa_KD2 <- list()
for(i in 2131:2859){
  aa_WT2[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                                                                             .fold_change = T)

  aa_WT2_par[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(aa_WT2)[i] <- paste0("experiment_", i)
  names(aa_WT2_par)[i] <- paste0("experiment_", i)
  #
  print("####################################################################################################")
  print("####################################################################################################")
  print("####################################################################################################")
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print("####################################################################################################")
  print("####################################################################################################")
  print("####################################################################################################")
  # Perform TF KD
  # KD_Performer(parent_Exp=i,
  #              TF_index=list(4),
  #              shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
  #              prev_na=5,
  #              run_GEMSTAT = T)
  aa_KD2[[i]] <- read_KD_results(experiment_nu = i,
                                 shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                                 experiment_output = aa_WT2[[i]])

  # Perform coop KD
  # KD_Performer_by_params(parent_Exp=i,
  #                        TF_names = names(TF.motifs.Shrinked.hocomoco.count),
  #                        Coop_mat_list=list(c(4,4)),
  #                        shared_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
  #                        prev_na=5,
  #                        run_GEMSTAT=F)
  # GEMSTAT_based_on_linear_exp_KD_COOP_Results[[i]] <-  read_KD_results(experiment_nu=i,
  #                                                                      shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
  #                                                                      coop_KD=T)
}

par(mfrow = c(1, 1), mar = c(6,4,4,4))
aa <- do.call(cbind, lapply(aa_WT2[c(2131:2859)], "[[", 7))
boxplot.matrix(cbind(GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh[[2129]]$Accuracy_All[aa_model_index_2129], aa),
               las= 2, main="Receptor activation profile performance")
abline(h=seq(0, 1, 0.01), col = 2, lty = 4)
heatmap.2(aa_modif_all_comb_p20_raw_2129,
          Rowv = F, Colv = F, dendrogram = "none", trace = 'none')
aa_q <- quantile(as.numeric(aa),probs = seq(0, 1, length.out = 11))
aa_dat <- data.frame( x=as.numeric(aa), above_90percentile= as.numeric(aa) >= aa_q[10] )
#library(ggplot2)
qplot(x,data=aa_dat,geom="histogram",fill=above_90percentile)

aaaa <- apply(X = aa, MARGIN = 2,which.max)
aaaa_max <- apply(X = aa, MARGIN = 2,max)
aaaa_max_ind <- sort(aaaa_max, decreasing = T, index.return = T)$ix
aa_modif_all_comb_p20_raw_sort_2129 <- aa_modif_all_comb_p20_raw_2129[aaaa_max_ind,]
#aa_modif_all_comb_p20_raw_sort <- aa_modif_all_comb_p20_raw_sort[, -3]
colnames(aa_modif_all_comb_p20_raw_sort_2129) <- c("AR", "GR", "PGR", "RARA", "RARG", "RXRA")

heatmap.2(aa_modif_all_comb_p20_raw_sort_2129, 
          Rowv = F, Colv = F, dendrogram = "none",
          trace = 'none', main = "RAP sorted by Performance")

aa_top <- aa_modif_all_comb_p20_raw_sort_2129[1:100,]
aa_top_table <- numeric(6*3)
for(i in 1:ncol(aa_top)){
  aa_tab <- table(aa_top[, i])
  for(j in 1:length(aa_tab)){
    aa_top_table[3 * (i - 1) + as.numeric(names(aa_tab)[j])] <- aa_tab[j]
  }
}

barplot(aa_top_table, 
        col=c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3), rep(6, 3)), 
        names.arg = c("ARp", "ARn", "ARall", "GRp", "GRn","GRall","PGRp", "PGRn",
                      "PGRall","RARAp", "RARAn","RARAall","RARGp", "RARGn","RARGall",
                      "RXRAp", "RXRAn", "RXRAall"), las=2, main = "top100_experiments")
### Filter for experiments that are ER sensitive

aa_WT_acc_2129 <- do.call(cbind, lapply(aa_WT2[c(2131:2859)], "[[", 7))
aa_ER_KD_Acc_2129 <- do.call(cbind, lapply(lapply(X = aa_KD2[c(2131:2859)], "[[", "4"), "[[", 7))
#aa_coop_ER_kd_acc_1398 <- do.call(cbind, lapply(lapply(X = GEMSTAT_based_on_linear_exp_KD_COOP_Results[c(1400:2128)], "[[", 1), "[[", 7))
#colnames(aa_coop_ER_kd_acc_1398) <- colnames(aa_ER_KD_Acc_1398)
aa_acc_diff_2129 <- aa_WT_acc_2129 - aa_ER_KD_Acc_2129

aa_dat <- data.frame( x=as.numeric(aa_acc_diff_2129),
                      above_10= as.numeric(aa_acc_diff_2129) >= 0.1 )
#library(ggplot2)
qplot(x,data=aa_dat,geom="histogram",
      fill=above_10, main = "WT_Accuracy - ER_KD_Accuracy")

# find the ER sensitive models based on the model prediction agreement between WT and KD instead of comparing their accuracy
aa_WT_KD_agreement <- matrix(nrow = length(aa_WT2[[2131]]$Accuracy_All),
                             ncol = 729)
for(i in 1:ncol(aa_WT_KD_agreement)){ # over RAPs
  for(j in 1:nrow(aa_WT_KD_agreement)){ # over models
    aa_wt <- aa_WT2[[2130 + i]]$Prediction_round_list[[j]]
    aa_k  <- aa_KD2[[2130 + i]][['4']]$Prediction_round_list[[j]]
    aa_WT_KD_agreement[j, i] <- sum(aa_wt == aa_k, na.rm = T)/sum(!is.na(aa_wt))
  }
}
rownames(aa_WT_KD_agreement) <- c(1:25)
colnames(aa_WT_KD_agreement) <- c(1:729)


par(mfrow = c(1, 1), mar = c(4,4,4,4))
boxplot.matrix(aa_WT_KD_agreement, las = 2, main = "WT_KD_prediction_agreement_all")

aa_q <- quantile(as.numeric(aa_WT_KD_agreement),probs = seq(0, 1, length.out = 11))
aa_dat <- data.frame( x=as.numeric(aa_WT_KD_agreement), below_10percentile= as.numeric(aa_WT_KD_agreement) <= aa_q[2] )
#library(ggplot2)
qplot(x,data=aa_dat,geom="histogram",fill=below_10percentile, main = "%agreement between KD and WT models")


aa_WT_KD_agreement_2129_5241 <- aa_WT_KD_agreement
aa_WT_KD_agreement_2129_5241[aa_WT_acc_2129 < 0.5241461] <- NA

aa_WT_KD_agreement_2129_5468 <- aa_WT_KD_agreement
aa_WT_KD_agreement_2129_5468[aa_WT_acc_2129 < 0.5468198] <- NA
hist(aa_WT_acc_2129[aa_WT_KD_agreement_2129_5468 < 0.6])
aa_WT_KD_agreement_2129_5468_6 <- aa_WT_KD_agreement_2129_5468
aa_WT_KD_agreement_2129_5468_7 <- aa_WT_KD_agreement_2129_5468

# aa_WT_KD_agreement_2129_5468_6 contains only the percentage agreement of models with accuracy higher than 0.5468 and agreement less than 0.6
aa_WT_KD_agreement_2129_5468_6[aa_WT_KD_agreement_2129_5468_6 > 0.6] <- NA
aa_WT_KD_agreement_2129_5468_7[aa_WT_KD_agreement_2129_5468_7 > 0.7] <- NA

colnames(aa_WT_KD_agreement_2129_5468_6) <- c(2131:2859)
colnames(aa_WT_KD_agreement_2129_5468) <- c(2131:2859)
colnames(aa_WT_KD_agreement_2129_5468_7) <- c(2131:2859)



boxplot.matrix(aa_WT_KD_agreement_2129_5468, las = 2)
sum(aa_WT_KD_agreement_2129_5241 < 0.50, na.rm =T) 
colnames(aa_WT_KD_agreement_2129_5241) <- c(2131:2859)
colnames(aa_WT_acc_2129) <- c(2131:2859)
colnames(aa_ER_KD_Acc_2129) <- c(2131:2859)

aa_WT_KD_agreement_2129_5241_10p <- aa_WT_KD_agreement_2129_5241
aa_WT_KD_agreement_2129_5241_20p <- aa_WT_KD_agreement_2129_5241
aa_WT_KD_agreement_2129_5241_40p <- aa_WT_KD_agreement_2129_5241
aa_WT_KD_agreement_2129_5241_lt80Perc <- aa_WT_KD_agreement_2129_5241
quantile(aa_WT_KD_agreement_2129_5241_10p, na.rm=T, probs = seq(0, 1, length.out = 11))
aa_WT_KD_agreement_2129_5241_lt80Perc[aa_WT_KD_agreement_2129_5241 > 0.8] <- NA
aa_WT_KD_agreement_2129_5241_10p[aa_WT_KD_agreement_2129_5241_10p > 0.8539458] <- NA
aa_WT_KD_agreement_2129_5241_20p[aa_WT_KD_agreement_2129_5241_20p > 0.8763251] <- NA
aa_WT_KD_agreement_2129_5241_40p[aa_WT_KD_agreement_2129_5241_40p > 0.8924028] <- NA


# aaKD_coop_dif <-aa_ER_KD_Acc_1398[!is.na(aa_acc_diff_6_10_1398)] -aa_coop_ER_kd_acc_1398[!is.na(aa_acc_diff_6_10_1398)]
# par(mfrow = c(1, 1), mar = rep(6, 4))
# hist(aaKD_coop_dif, main = "ERKD_acc - ERCOOPKD_acc", breaks = 10)
aa_WT_KD_agreement_2129_5241_10p <- aa_WT_KD_agreement_2129_5241_10p[,! colSums(is.na(aa_WT_KD_agreement_2129_5241_10p)) == nrow(aa_WT_KD_agreement_2129_5241_10p)]
aa_WT_KD_agreement_2129_5241_20p  <- aa_WT_KD_agreement_2129_5241_20p[,! colSums(is.na(aa_WT_KD_agreement_2129_5241_20p)) == nrow(aa_WT_KD_agreement_2129_5241_20p)]
aa_WT_KD_agreement_2129_5241_40p <- aa_WT_KD_agreement_2129_5241_40p[, ! colSums(is.na(aa_WT_KD_agreement_2129_5241_40p)) == nrow(aa_WT_KD_agreement_2129_5241_40p)]
aa_WT_KD_agreement_2129_5241_lt80Perc <- aa_WT_KD_agreement_2129_5241_lt80Perc[, ! colSums(is.na(aa_WT_KD_agreement_2129_5241_lt80Perc)) == nrow(aa_WT_KD_agreement_2129_5241_lt80Perc)]
aa_WT_KD_agreement_2129_5468_6 <- aa_WT_KD_agreement_2129_5468_6[, ! colSums(is.na(aa_WT_KD_agreement_2129_5468_6)) == nrow(aa_WT_KD_agreement_2129_5468_6)]
aa_WT_KD_agreement_2129_5468_7 <- aa_WT_KD_agreement_2129_5468_7[, ! colSums(is.na(aa_WT_KD_agreement_2129_5468_7)) == nrow(aa_WT_KD_agreement_2129_5468_7)]

colnames(aa_modif_all_comb_p20_raw_2129)<-c("AR", "GR", "PGR", 
                                            "RARA", "RARG", "RXRA")


heatmap.2(aa_modif_all_comb_p20_raw_2129[colnames(aa_WT_KD_agreement_2129_5241_lt80Perc),]
          , Rowv = T, Colv = F, dendrogram = "none",
          trace = 'none', distfun = my_dist_fuc,
          main = "RAP acc > 0.5241 and ER-sensitive", col = c("green", "Red", "White"))

heatmap.2(aa_modif_all_comb_p20_raw_2129[colnames(aa_WT_KD_agreement_2129_5241_10p),]
          , Rowv = T, Colv = T, dendrogram = "both",
          trace = 'none', distfun = my_dist_fuc, col = c("green", "Red", "White"))

heatmap.2(aa_modif_all_comb_p20_raw_2129[colnames(aa_WT_KD_agreement_2129_5241_40p),]
          , Rowv = T, Colv = T, dendrogram = "both",
          trace = 'none', distfun = my_dist_fuc, col = c("green", "Red", "White"))

heatmap.2(aa_modif_all_comb_p20_raw_2129[colnames(aa_WT_KD_agreement_2129_5468_6),]
          , Rowv = T, Colv = T, dendrogram = "both",
          trace = 'none', distfun = my_dist_fuc, col = c("green", "Red", "White"))

heatmap.2(aa_modif_all_comb_p20_raw_2129[colnames(aa_WT_KD_agreement_2129_5468_7),]
          , Rowv = T, Colv = T, dendrogram = "both",
          trace = 'none', distfun = my_dist_fuc, col = c("green", "Red", "White"))

aa_top <- aa_modif_all_comb_p20_raw_2129[colnames(aa_WT_KD_agreement_2129_5241_lt80Perc),]
aa_top_table <- numeric(6*3)
for(i in 1:ncol(aa_top)){
  aa_tab <- table(aa_top[, i])
  for(j in 1:length(aa_tab)){
    aa_top_table[3 * (i - 1) + as.numeric(names(aa_tab)[j])] <- aa_tab[j]
  }
}
par(mfrow = c(1, 1), mar = rep(6, 4))
barplot(aa_top_table, 
        col=c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3), rep(6, 3)), 
        names.arg = c("ARp", "ARn", "ARall", "GRp", "GRn","GRall","PGRp", "PGRn",
                      "PGRall","RARAp", "RARAn","RARAall","RARGp", "RARGn","RARGall",
                      "RXRAp", "RXRAn", "RXRAall"), las=2, 
        main = "WT acc > 0.52, ER KD agreement < 0.8")

aaaaaa <- aa_WT_acc_2129[, colnames(aa_WT_KD_agreement_2129_5468_6)]
par(mfrow = c(6, 4), mar = c(2, 3, 1, 1))
for(i in 1:ncol(aa_WT_KD_agreement_2129_5468_6)){
  plot(aa_WT_acc_2129[, colnames(aa_WT_KD_agreement_2129_5468_6)[i]], 
       aa_ER_KD_Acc_2129[, colnames(aa_WT_KD_agreement_2129_5468_6)[i]],
       xlim = c(0.4, 0.6), ylim = c(0.4, 0.6), ylab = "KD", col = 2)
  abline(c(0, 1), col=1)
  # points(aa_WT_acc_1398[, colnames(aa_WT_KD_agreement_2129_5241_40p)[i]],
  #        aa_coop_ER_kd_acc_1398[, colnames(aa_WT_KD_agreement_2129_5241_40p)[i]], pch = 4, col = 3)
}

aa_mod_ind <- list()
for(i in 1:ncol(aa_WT_KD_agreement_2129_5241_lt80Perc)){
  aa_mod_ind[[i]] <- which(! is.na(aa_WT_KD_agreement_2129_5241_lt80Perc[, i]))
}
names(aa_mod_ind) <- colnames(aa_WT_KD_agreement_2129_5241_lt80Perc)
#aa_base_linear_model <- aa_model_index_2129[unlist(aa_mod_ind)]

# look at the parameters
aa_bindings <- matrix(nrow = length(unlist(aa_mod_ind)), ncol = 19)
colnames(aa_bindings) <- colnames(aa_WT2_par[[2131]]$binding)
rownames(aa_bindings) <- character(nrow(aa_bindings))

aa_alphas <- matrix(nrow = length(unlist(aa_mod_ind)), ncol =19)
colnames(aa_alphas) <- colnames(aa_WT2_par[[2131]]$alpha_effective)
rownames(aa_alphas) <- character(nrow(aa_alphas))

aa_coops <- matrix(nrow = length(unlist(aa_mod_ind)), ncol =6)
colnames(aa_coops) <- names(aa_WT2_par[[2131]]$coop[[8]])
rownames(aa_coops) <- character(nrow(aa_coops))

aa_cnt = 1
for(i in 1:length(aa_mod_ind)){
  for(j in 1:length(aa_mod_ind[[i]])){
    aa_bindings[aa_cnt, ] <- aa_WT2_par[[as.numeric(names(aa_mod_ind))[i]]]$binding[aa_mod_ind[[i]][j], ]
    aa_alphas[aa_cnt, ] <- aa_WT2_par[[as.numeric(names(aa_mod_ind))[i]]]$alpha_effective[aa_mod_ind[[i]][j], ]
    aa_coops[aa_cnt, ] <- aa_WT2_par[[as.numeric(names(aa_mod_ind))[i]]]$coop[[aa_mod_ind[[i]][j]]]
    rownames(aa_bindings)[aa_cnt] <- paste0(as.numeric(names(aa_mod_ind))[i], "_", aa_mod_ind[[i]][j])
    rownames(aa_alphas)[aa_cnt] <- paste0(as.numeric(names(aa_mod_ind))[i], "_", aa_mod_ind[[i]][j])
    rownames(aa_coops)[aa_cnt] <- paste0(as.numeric(names(aa_mod_ind))[i], "_", aa_mod_ind[[i]][j])
    aa_cnt <- aa_cnt + 1
  }
  
}
par(mfrow = c(1, 1), mar = c(6,4,4,4))
boxplot.matrix(aa_bindings, las = 2)
boxplot.matrix(aa_alphas, las = 2)
boxplot.matrix(aa_coops, las = 2)

#aa_exp_names <- unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 1))
#aa_col <- sample(x = col_vector, size = 5, replace = F)

# accept only the models that have ER as an activator
# aa_bindings <- aa_bindings[aa_alphas[, 4] > 1, ]
# aa_coops <- aa_coops[aa_alphas[, 4] > 1, ]
# aa_alphas <- aa_alphas[aa_alphas[, 4] > 1, ]

aa_exp_names<- unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 1))
aa_model_names <- unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 2))
aa_raw_RAP_comb_chosen <- aa_modif_all_comb_p20_raw_2129[aa_exp_names,]

aa_exp_model <- list()
for(i in 1:length(unique(aa_exp_names))){
  aa_exp_model[[i]] <- aa_model_names[aa_exp_names == unique(aa_exp_names)[i]]
}
names(aa_exp_model) <- unique(aa_exp_names)

#aa_exp_model <- aa_exp_model[12:length(unique(aa_exp_names))]

heatmap.2(cbind(aa_alphas, aa_modif_all_comb_p20_raw_2129[aa_exp_names,])
          , Rowv = T, Colv = F, dendrogram = "row",
          trace = 'none', breaks = c(0, 0.5, 1, 2.5, 3, 7, 10),
          margins = c(8,8), col = col_vector[11:16], main = "alphaERsensitive+concComb_2129"
          #,RowSideColors = col_vector[(aa_base_linear_model + 10)%%length(col_vector)]
          )

heatmap.2(aa_coops
          , Rowv = T, Colv = T, dendrogram = "none",
          trace = 'none', breaks = c(100000, 300000, 600000, 900000, 1000000),
          margins = c(8,8), col = col_vector[13:16], main="CoopERsensitive"
          #,RowSideColors = col_vector[(aa_base_linear_model + 10)%%length(col_vector)]
          )

heatmap.2(cbind(log10(aa_coops), aa_modif_all_comb_p20_raw_2129[aa_exp_names,])
          , Rowv = F, Colv = F, dendrogram = "none",
          trace = 'none', breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5,  5.2, 5.4, 5.6, 5.8, 6),
          margins = c(8,8), main="CoopERsensitive", col = col_vector[51:60]
          #,RowSideColors = col_vector[(aa_base_linear_model + 10)%%length(col_vector)]
          )


heatmap.2(cbind(aa_bindings, aa_modif_all_comb_p20_raw_2129[aa_exp_names,])
          , Rowv = T, Colv = F, dendrogram = "none",
          trace = 'none', breaks = c(0, 1, 2, 3, 10, 50),
          col = col_vector[11:15], margins = c(8,8), 
          main="BindingERsensitive")

# perform Other KnockDowns for these models
# save(list = c("aa_exp_model", "KD_Performer"), file = "aa_exp_KDperf.RData")

aadone <- unlist(lapply(aa_KD2, length))
aadone <- which(aadone == 19)
aa_done_kds <- list()
for(i in 1:length(aadone)){
  aa_done_kds[[i]] <- names(aa_KD2[[aadone[1]]]$`19`$threshold)
}
names(aa_done_kds) <- aadone

aa_exp_model_exist <- aa_exp_model

for(i in 1:length(aa_exp_model)){
  if(names(aa_exp_model)[i] %in% names(aa_done_kds)){
    for(j in 1:length(aa_exp_model[[i]])){
      aa_exp_model_exist[[i]][j] <- aa_exp_model[[i]][j] %in% aa_done_kds[[names(aa_exp_model)[i]]]
    }
  }else{
    aa_exp_model_exist[[i]][] <- F
  }

}


aa_exp_model
aa_GEMSTAT_based_on_linear_exp_KD_Results <- list()
aa_GEMSTAT_based_on_linear_exp_KD_COOP_Results <- list()
for(i in 1:length(aa_exp_model)){
  if(! all(aa_exp_model_exist[[i]])){
    print("####################################################################################################")
    print("####################################################################################################")
    print("####################################################################################################")
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print("####################################################################################################")
    print("####################################################################################################")
    print("####################################################################################################")
    # Perform TF KD
    #  if(i %in% unique(aa_exp_names)[12:length(unique(aa_exp_names))]){
    KD_Performer(parent_Exp= as.integer(names(aa_exp_model))[i],
                 TF_index=as.list(c(c(1:3), c(5:19))),
                 model_index = aa_exp_model[[i]],
                 shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                 prev_na=5,
                 run_GEMSTAT = T)
    #  }
    
    aa_GEMSTAT_based_on_linear_exp_KD_Results[[i]] <- read_KD_results(experiment_nu =  as.integer(names(aa_exp_model))[i],
                                                                      shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                                                                      experiment_output = aa_WT2[[as.integer(names(aa_exp_model))[i]]])
  }
}


for(i in 1:length(aa_exp_model)){
  if(! all(aa_exp_model_exist[[i]])){
    aa_KD2[[as.integer(names(aa_exp_model))[i]]] <- aa_GEMSTAT_based_on_linear_exp_KD_Results[[i]]
  }
}

save(list = c("aa_KD2", "aa_WT2", "aa_WT2_par"), 
     file = "WT_KD_after_debug_thresh.RData")



aa <- lsf.str()
save(list = aa, file = "All_functions_02_11_2019.RData")
# length(aa_KD2)
# create a list of WT and a list of list of KD predictions with the same names
aa_wt_list <- list()
aa_cnt <- 1
aa_mynames <- character(0)
for(i in 1:length(aa_exp_model)){
  for(j in 1:length(aa_exp_model[[i]])){
    aa_wt_list[[aa_cnt]] <- aa_WT2[[as.integer(names(aa_exp_model))[i]]]$Prediction_round_list[[aa_exp_model[[i]][j]]]
    aa_mynames <- c(aa_mynames, paste0(names(aa_exp_model)[i], "_", aa_exp_model[[i]][j]))
    aa_cnt <- aa_cnt + 1
  }
}
names(aa_wt_list) <- aa_mynames

aa_kd_list_list <- list()
for(i in 1:length(TF.motifs.Shrinked.hocomoco_E2F1_added.count)){
  aa_kd_list_list[[i]] <- list()
  aa_cnt <- 1
  aa_mynames <- character(0)
  for(j in 1:length(aa_exp_model)){
    for(k in 1:length(aa_exp_model[[j]])){
        aa_kd_list_list[[i]][[aa_cnt]] <- aa_KD2[[as.integer(names(aa_exp_model))[j]]][[as.character(i)]]$Prediction_round_list[[as.character(aa_exp_model[[j]][k])]]
        aa_mynames <- c(aa_mynames, paste0(names(aa_exp_model)[j], "_", aa_exp_model[[j]][k]))
        aa_cnt <- aa_cnt + 1
        #stopifnot(length(aa_KD2[[as.integer(names(aa_exp_model))[j]]][[as.character(i)]]$Prediction_round_list) == length(aa_exp_model[[j]]))
    }# end of loop over mods
    
  } # end of loop over exps
  names(aa_kd_list_list[[i]]) <- aa_mynames
}
names(aa_kd_list_list) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)

sum(do.call(rbind, aa_KD2[[2139]][['4']]$threshold) == do.call(rbind, aa_WT2[[2139]]$threshold))

library(gplots)
#aa_pr <- aa_WT2[[2139]]$Prediction_round_list
#aa_prkd <- lapply(aa_KD2[[2139]], "[[", 2)
aa <- KD_on_gene_evaluator_multiTF_wrapper(.WT_Prediction_list = aa_wt_list ,
                                           KD_Prediction_list_list = aa_kd_list_list,
                                           .Avg_over_models=T, 
                                         KD_names = names(aa_kd_list_list),
                                           plt = T)
aa1 <- KD_on_gene_evaluator_multiTF_wrapper(.WT_Prediction_list = aa_wt_list ,
                                           KD_Prediction_list_list = aa_kd_list_list,
                                           .Avg_over_models=F, 
                                           KD_names = names(aa_kd_list_list),
                                           plt = T, plot_gene_index = c(1, 2, 3, 4))

aa_bind_comb <- cbind(aa_bindings,
                      aa_modif_all_comb_p20_raw_2129[aa_exp_names,])


summary(aa_bind_comb[aa_bind_comb[, 22] == 1, 13])
summary(aa_bind_comb[aa_bind_comb[, 22] == 2, 13])
par(mfrow = c(1,1), mar = c(4,4,4,4))
boxplot(aa_bind_comb[,13] ~ aa_bind_comb[, 22],
        xlab = "PGR_RAP", ylab = "PGR binding par",
        main = "PGR binding grouped by PGR RAP")

aa_alpha_comb <- cbind(aa_alphas,
                      aa_modif_all_comb_p20_raw_2129[aa_exp_names,])


summary(aa_alpha_comb[aa_alpha_comb[, 22] == 1, 13])
summary(aa_alpha_comb[aa_alpha_comb[, 22] == 2, 13])
par(mfrow = c(1,1), mar = c(4,4,4,4))
boxplot(aa_alpha_comb[,13] ~ aa_alpha_comb[, 22],
        xlab = "PGR_RAP", ylab = "PGR alpha par",
        main = "PGR alpha grouped by PGR RAP")

sort(as.integer(unique(unlist(aa_exp_model))))
GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Annot_per_model_perGene[[1]][[3]]
GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Annot_per_model_perGene[[12]][[19]]
1:length(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
for(iii in c(13)){
  aa <- annotation_overlap_finder(annoatation_witer_input = GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Annot_per_model_perGene,
                                  TF_names = names(TF.motifs.Shrinked.hocomoco_E2F1_added.count),
                                  TF_index = iii)
  
  aa_t_models <- sort(as.integer(unique(unlist(aa_exp_model))))
  par(mfrow = c(4, 2), mar = c(1, 1,1,1))
  for(i in aa_t_models){
    aa_gal <- sum(GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Site_count_mat[[i]][, iii] > 0)
    aa_ov <- sum(unlist(lapply(aa[[i]], length)) > 0)
    pie(table(unlist(aa[[i]])),
        main = paste0(aa_ov, "/", aa_gal))
    
  }
}


par(mfrow = c(4, 2), mar = c(1, 1,1,1))
for(i in aa_t_models){
  aa_gal <- sum(GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Site_count_mat[[i]][, 4] > 0)
  aa_ov <- sum(unlist(lapply(aa[[i]], length)) == 1)
  pie(table(unlist(aa[[i]][unlist(lapply(aa[[i]], length))==1])),
      main = paste0(aa_ov, "/", aa_gal))
  
}
sum(GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Site_count_mat[[1]][, 4] == 0)
pie(table(unlist(aa[[1]])))
table(unlist(aa[[i]][unlist(lapply(aa[[i]], length))==1]))




##########################################
# create model summary:
# 1. linear base model number
# 2. % genes with >= 1 ERBS
# 3. a Table with one row per TF and the following columns:
#  KD rank, RAP, Alpha, Binding, selfCoop
aa_base_linear_model <- aa_model_index_2129[as.integer(aa_model_names)]

aa_nz_all_mat <- matrix(nrow = 78, ncol = 19)
for(i in 1:78){
  for(j in 1:19){
    aa_nz_all_mat[i,j] <- sum(GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Site_count_mat[[i]][, j] > 0)/
      nrow(GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Site_count_mat[[i]])
    
  }
}
colnames(aa_nz_all_mat) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
colnames(aa_raw_RAP_comb_chosen_fil)[2] <- "NR3C1"
par(mfrow = c(1,1), mar = c(6,6,6,6))
boxplot.matrix(aa_nz_all_mat, las =2, main = "%Genes containing at least one BS")
#aa_nz_all[aa_base_linear_model]
aa_GEMSTAT_filtered_summary <- list()
for(i in 1:length(aa_base_linear_model)){
  aa_GEMSTAT_filtered_summary[[i]] <- list()
  aa_GEMSTAT_filtered_summary[[i]][[1]] <- aa_base_linear_model[i]
  aa_GEMSTAT_filtered_summary[[i]][[2]] <- aa_nz_all_mat[aa_GEMSTAT_filtered_summary[[i]][[1]], ]
  names(aa_GEMSTAT_filtered_summary[[i]][[2]]) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
  aa_GEMSTAT_filtered_summary[[i]][[3]] <- matrix(nrow = ncol(aa_alphas), ncol =5)
  rownames(aa_GEMSTAT_filtered_summary[[i]][[3]]) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
  colnames(aa_GEMSTAT_filtered_summary[[i]][[3]]) <- c("KD_rank", "RAP", "Alpha", "Binding", "SelfCoop")
  names(aa_GEMSTAT_filtered_summary[[i]]) <- c("Base_Model", "Percent_ERBS", "TF_mat")
  
  for(j in 1:nrow(aa_GEMSTAT_filtered_summary[[i]][[3]])){
    # setting rank
    #aa_GEMSTAT_filtered_summary[[i]][[3]][j, 1] <- aa_KD_eff_mat_rank[i, j]
    aa_GEMSTAT_filtered_summary[[i]][[3]][j, 1] <- NA
    # setting RAP
    if(names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[j] %in% colnames(aa_raw_RAP_comb_chosen)){
      aa_GEMSTAT_filtered_summary[[i]][[3]][j, 2] <- aa_raw_RAP_comb_chosen[i, names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[j]]
    }
    # set alpha
    aa_GEMSTAT_filtered_summary[[i]][[3]][j, 3] <- aa_alphas[i, j]
    
    # set binding 
    aa_GEMSTAT_filtered_summary[[i]][[3]][j, 4] <- aa_bindings[i, j]
    
    # set coop 
    if( names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[j] %in% unlist(lapply(strsplit(colnames(aa_coops), split = ":"), "[[", 1))){
      aa_GEMSTAT_filtered_summary[[i]][[3]][j, 5] <- aa_coops[i, which(unlist(lapply(strsplit(colnames(aa_coops), split = ":"), "[[", 1)) %in% names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[j] )]
    }
  }
}
names(aa_GEMSTAT_filtered_summary) <- rownames(aa_alphas)

aa_GEMSTAT_filtered_summary$`1434_34`$TF_mat

aa_GEMSTAT_filtered_summary_per_TF <- list()
for(i in 1:nrow(aa_GEMSTAT_filtered_summary[[1]]$TF_mat)){
  aa_GEMSTAT_filtered_summary_per_TF[[i]] <- matrix(nrow = length(aa_GEMSTAT_filtered_summary),
                                                    ncol = ncol(aa_GEMSTAT_filtered_summary[[1]]$TF_mat) + 1)
  rownames(aa_GEMSTAT_filtered_summary_per_TF[[i]]) <- names(aa_GEMSTAT_filtered_summary)
  colnames(aa_GEMSTAT_filtered_summary_per_TF[[i]]) <- c(colnames(aa_GEMSTAT_filtered_summary[[1]]$TF_mat), "Percent_gene_with_site")
  for(j in 1:length(aa_GEMSTAT_filtered_summary)){
    aa_GEMSTAT_filtered_summary_per_TF[[i]][j, 1:ncol(aa_GEMSTAT_filtered_summary[[j]]$TF_mat)] <- aa_GEMSTAT_filtered_summary[[j]]$TF_mat[i, ]
    aa_GEMSTAT_filtered_summary_per_TF[[i]][j, ncol(aa_GEMSTAT_filtered_summary_per_TF[[i]])] <- aa_GEMSTAT_filtered_summary[[j]][[2]][i]
  }
  
}
names(aa_GEMSTAT_filtered_summary_per_TF) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)

# plot per RAP
for(i in c(1, 10, 13, 14, 15, 17)){
  par(mfrow = c(5,4), mar = c(3,2,1,1))
  for(j in 1:19){
    boxplot(aa_GEMSTAT_filtered_summary_per_TF[[j]][, 3]~aa_GEMSTAT_filtered_summary_per_TF[[i]][, 2],
            main=paste(names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[j],
                       names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[i], sep = "_"))
    abline(h = 0, col = 2, lwd= 1, lty = 3)
  }
}

for(i in c(1, 10, 13, 14, 15, 17)){
  par(mfrow = c(2,3), mar = c(3,2,1,1))
  for(j in c(1, 10, 14, 15, 17)){
    boxplot(log10(aa_GEMSTAT_filtered_summary_per_TF[[j]][, 5])~aa_GEMSTAT_filtered_summary_per_TF[[i]][, 2],
            main=paste(names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[j],
                       names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[i], sep = "_"))
    abline(h = 0, col = 2, lwd= 1, lty = 3)
  }
}

for(i in c(1, 10, 13, 14, 15, 17)){
  par(mfrow = c(2,3), mar = c(3,2,1,1))
  for(j in c(1, 10, 13, 14, 15, 17)){
    boxplot(aa_GEMSTAT_filtered_summary_per_TF[[j]][, 2]~aa_GEMSTAT_filtered_summary_per_TF[[i]][, 2],
            main=paste(names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[j],
                       names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[i], sep = "_"))
    # abline(h = 0, col = 2, lwd= 1, lty = 3)
  }
}


########################################################################################################################
aa_exp_names <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 1))) - 2130
aa_model_names <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 2)))
aa_model_names_org <- aa_model_index_2129[aa_model_names]

colnames(aa_modif_all_comb_p20_raw_2129)[2] <- "NR3C1"

GeneFoldChange_by_TFFoldChange_model_RAP_wrapper(.Gene_exp_mat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                                                 .TF_exp_mat = TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01,
                                                 TF_gene_site_mat_list = GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Site_count_mat, 
                                                 .site_nu_thresh = 1,
                                                 # .TF_index = c(1:nrow(.TF_exp_mat)),
                                                 RAP_raw_mat=aa_modif_all_comb_p20_raw_2129, 
                                                 RAP_list = aa_modif_all_options_p20_2129,
                                                 RAP_nu=aa_exp_names[13],
                                                 model_nu = aa_model_names_org[13])
########################################################################################################################
aa_WT_Pred_mat_list <- list()
aa_KD_Pred_mat_list <- list()
aa_exp_names <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphasl), split = "_"), "[[", 1)))
aa_model_names <- as.numeric( unlist(lapply(strsplit(rownames(aa_alphas), split = "_"), "[[", 2)))
aa_model_names_org <- aa_model_index_2129[aa_model_names]
par(mfrow = set_row_col_nu(length(TF.motifs.Shrinked.hocomoco_E2F1_added.count)),
    mar= c(4, 4, 4, 3))
aa <- list()
for(j in 1:length(TF.motifs.Shrinked.hocomoco_E2F1_added.count)){
  aa_TF_index <- j
  for(i in 1:length(aa_exp_names)){
    aa_WT_Pred_mat_list[[i]] <- aa_WT2[[aa_exp_names[i]]]$Prediction_round_list[[aa_model_names[i]]]
    aa_KD_Pred_mat_list[[i]] <- aa_KD2[[aa_exp_names[i]]][[as.character(aa_TF_index)]]$Prediction_round_list[[as.character(aa_model_names[i])]]
  }
  aa[[j]] <- KD_by_affected_enhancers(WT_Pred_mat_list = aa_WT_Pred_mat_list, 
                                      KD_Pred_mat_list = aa_KD_Pred_mat_list,
                                      real_mat = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172,
                                      TF_gene_site_mat_list = GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Site_count_mat,
                                      site_nu_thresh=1,
                                      model_index = aa_model_names_org, 
                                      TF_index=aa_TF_index,
                                      plot_sc = F)
  plot(aa[[j]][[1]], aa[[j]][[3]], 
       main = names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)[j], 
       ylab = "agree_target_genes", xlab = "#targetted genes", ylim = c(0.4, 1))
  abline(h = seq(0.4, 1, 0.1), lty = 3, lwd = 0.6, col = 3)
}
########################################################################################################################
# read the raw KD results for the 59 models
aa_raw_kd <- list()
for(i in 1:length(rownames(aa_alphas))){
  aa_e <- unlist(strsplit(rownames(aa_alphas)[i], split = "_"))[1]
  aa_m <- unlist(strsplit(rownames(aa_alphas)[i], split = "_"))[2]
  aa_raw_kd[[i]] <- list()
  for(j in 1:length(TF.motifs.Shrinked.hocomoco_E2F1_added.count)){
    aa_raw_kd[[i]][[j]] <- GEMSTAT_output_Reader_multiEnh(expression_file = paste0("Seeded_GEMSTAT_ens/Experiment_", aa_e, "/KnockDown/", j, "/Outputs/Experiment_", aa_e, "_", aa_m ,".out"), 
                                                          .fixed_thresh = numeric(0), 
                                                          fold_change = F)
  }
  names(aa_raw_kd[[i]]) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)
}
names(aa_raw_kd) <- rownames(aa_alphas)

aa_raw_wt <- list()
for(i in 1:length(rownames(aa_alphas))){
  aa_e <- unlist(strsplit(rownames(aa_alphas)[i], split = "_"))[1]
  aa_m <- unlist(strsplit(rownames(aa_alphas)[i], split = "_"))[2]
  aa_raw_wt[[i]] <- list()
  aa_raw_wt[[i]] <- GEMSTAT_output_Reader_multiEnh(expression_file = paste0("Seeded_GEMSTAT_ens/Experiment_", aa_e, "/Outputs/Experiment_", aa_e, "_", aa_m ,".out"), 
                                                          .fixed_thresh = numeric(0), 
                                                          fold_change = F)$prediction
  
}
names(aa_raw_wt) <- rownames(aa_alphas)

aa_kd_list_list_raw <- list()
for(i in 1:length(TF.motifs.Shrinked.hocomoco_E2F1_added.count)){
  aa_kd_list_list_raw[[i]] <- list()
  for(j in 1:length(rownames(aa_alphas))){
    aa_kd_list_list_raw[[i]][[j]] <- aa_raw_kd[[j]][[i]]$prediction
  } # end of loop over exp_models
  names(aa_kd_list_list_raw[[i]]) <- rownames(aa_alphas)
}
names(aa_kd_list_list_raw) <- names(TF.motifs.Shrinked.hocomoco_E2F1_added.count)

library(gplots)
aa <- KD_on_gene_evaluator_multiTF_wrapper(.WT_Prediction_list = aa_raw_wt ,
                                           KD_Prediction_list_list = aa_kd_list_list_raw,
                                           .Avg_over_models=T, 
                                           KD_names = names(aa_kd_list_list_raw),
                                           plt = T, 
                                           #.breaks = seq(-0.4, 0.2, length.out = 11), 
                                          # .colors = rainbow( n = 10),
                                           ct_treat_map = rep(c(0, 1), 45))
par(mfrow = c(1, 1), mar = c(4,4,4,4))
boxplot.matrix(aa, las = 2)


######## Looking at thresholds trying to debug the fold change conversion since it seems to be buggy
aacnt <- 0
aa_which <- integer(729)
for(i in 2131:2859){
  aa <- (do.call(rbind, aa_WT2[[i]]$threshold))
  aacnt = aacnt + sum(aa[, 2] < 0)
  aa_which[i - 2130 ] <- sum(aa[, 2] < 0)
}
sum(aa_which > 0)

aa_WT_newdis <- list()
aa_KD_newdis <- list()

for(i in as.integer(names(aa_exp_model[1:2]))){
  aa_WT_newdis[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                               .fold_change = T, ..fixed_thresh = numeric(0))
  aa_KD_newdis[[i]] <- read_KD_results(experiment_nu =  i,
                                       shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                                       experiment_output = aa_WT_newdis[[i]])
  
}
# check whether the KD of ER goes down

sum(aa_KD_newdis[[2305]]$`4`$Prediction_raw_list$`12` == aa_KD2[[2305]]$`4`$Prediction_raw_list$`12`, na.rm = T)
sum(aa_KD_newdis[[2305]]$`4`$Prediction_round_list$`12` == aa_KD2[[2305]]$`4`$Prediction_round_list$`12`, na.rm = T)

aa_KD_newdis[[2305]]$`4`$threshold$`12`
aa_KD2[[2305]]$`4`$threshold$`12`
aa_raw_wt1 <- GEMSTAT_output_Reader_multiEnh(expression_file = paste0("Seeded_GEMSTAT_ens/Experiment_", 2305, "/Outputs/Experiment_", 2305, "_", 12 ,".out"), 
                                                 .fixed_thresh = numeric(0), 
                                                 fold_change = F)
aat <- aa_raw_wt1$ground_truth
aat2 <- matrix(nrow= 172, ncol = 45)
for(i in 1:45){
  aat2[, i] <- log2(aat[, 2*i] / aat[, (2*i - 1)])
}
sum(aat2 == my_CommonDifExpMat_16_ERassoc_gte5nzAll_172, na.rm = T) / sum(!is.na(my_CommonDifExpMat_16_ERassoc_gte5nzAll_172))
########################################################################################################################
#Perform double KD
aa_exp_model
aa_GEMSTAT_based_on_linear_exp_KD_Results <- list()
aa_GEMSTAT_based_on_linear_exp_KD_COOP_Results <- list()
for(i in 1:length(aa_exp_model)){
    print("####################################################################################################")
    print("####################################################################################################")
    print("####################################################################################################")
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print(i)
    print("####################################################################################################")
    print("####################################################################################################")
    print("####################################################################################################")
    # Perform TF KD
    #  if(i %in% unique(aa_exp_names)[12:length(unique(aa_exp_names))]){
    KD_Performer(parent_Exp= as.integer(names(aa_exp_model))[i],
                 TF_index=list(c(4, 6),
                               c(4, 7),
                               c(4, 11),
                               c(4, 13),
                               c(4, 17), 
                               c(4, 18), 
                               c(4, 19)),
                 model_index = aa_exp_model[[i]],
                 shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                 prev_na=5,
                 run_GEMSTAT = T)
    #  }
    
    aa_GEMSTAT_based_on_linear_exp_KD_Results[[i]] <- read_KD_results(experiment_nu =  as.integer(names(aa_exp_model))[i],
                                                                      shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                                                                      experiment_output = aa_WT2[[as.integer(names(aa_exp_model))[i]]],
                                                                      KD_name = c("4_6", "4_7", "4_11", "4_13", "4_17", "4_18", "4_19"))
}

load("Double_KD_res.RData")

for(i in 1:length(aa_exp_model)){
  aa_KD2[[as.integer(names(aa_exp_model)[i])]] <- c(aa_KD2[[as.integer(names(aa_exp_model)[i])]], 
                                                    aa_GEMSTAT_based_on_linear_exp_KD_Results[[i]])
}

aa_wt_list <- list()
aa_cnt <- 1
aa_mynames <- character(0)
for(i in 1:length(aa_exp_model)){
  for(j in 1:length(aa_exp_model[[i]])){
    aa_wt_list[[aa_cnt]] <- aa_WT2[[as.integer(names(aa_exp_model))[i]]]$Prediction_round_list[[aa_exp_model[[i]][j]]]
    aa_mynames <- c(aa_mynames, paste0(names(aa_exp_model)[i], "_", aa_exp_model[[i]][j]))
    aa_cnt <- aa_cnt + 1
  }
}
names(aa_wt_list) <- aa_mynames

aa_kd_list_list <- list()
aanames <- c(c(1:19), 
             c("4_6", "4_7", "4_11", "4_13", "4_17", "4_18", "4_19"))
for(i in 1:length(aanames)){
  aa_kd_list_list[[i]] <- list()
  aa_cnt <- 1
  aa_mynames <- character(0)
  for(j in 1:length(aa_exp_model)){
    for(k in 1:length(aa_exp_model[[j]])){
      aa_kd_list_list[[i]][[aa_cnt]] <- aa_KD2[[as.integer(names(aa_exp_model))[j]]][[aanames[i]]]$Prediction_round_list[[as.character(aa_exp_model[[j]][k])]]
      aa_mynames <- c(aa_mynames, paste0(names(aa_exp_model)[j], "_", aa_exp_model[[j]][k]))
      aa_cnt <- aa_cnt + 1
      #stopifnot(length(aa_KD2[[as.integer(names(aa_exp_model))[j]]][[as.character(i)]]$Prediction_round_list) == length(aa_exp_model[[j]]))
    }# end of loop over mods
    
  } # end of loop over exps
  names(aa_kd_list_list[[i]]) <- aa_mynames
}
names(aa_kd_list_list) <- c(names(TF.motifs.Shrinked.hocomoco_E2F1_added.count), 
                            c("ER_GATA3", "ER_JUN", "ER_NR5A2", "ER_PGR", "ER_RXRA", "ER_SP1", "ER_TFAP2C"))


library(gplots)
#aa_pr <- aa_WT2[[2139]]$Prediction_round_list
#aa_prkd <- lapply(aa_KD2[[2139]], "[[", 2)
aa <- KD_on_gene_evaluator_multiTF_wrapper(.WT_Prediction_list = aa_wt_list ,
                                           KD_Prediction_list_list = aa_kd_list_list,
                                           .Avg_over_models=T, 
                                           KD_names = names(aa_kd_list_list),
                                           plt = T)

aamymat <- aa[, c("ESR1", "PGR", "ER_PGR")]
aamybreaks <- seq(-2, 2, length.out = 21)
aamycol <- colorspace::diverge_hsv(20)
heatmap.2(aamymat,
          Rowv = T, 
          Colv = F, 
          dendrogram = 'row',
          breaks = aamybreaks, 
          col = aamycol,
          margins = c(8, 8), trace = F)
########################################################################################################################
# add cooperativity between ER and PGR
# Do it for the RAPs that have ER and PGR in the same phase
# Do it for all 78 models.

aa_modif_all_comb_p20_raw_2129_CHOSEN <- aa_modif_all_comb_p20_raw_2129[as.integer(names(aa_exp_model)) - 2130, ]
aa_modif_all_comb_p20_2129_CHOSEN <- aa_modif_all_comb_p20_2129[as.integer(names(aa_exp_model)) - 2130]
aa_modif_all_comb_p20_2129_CHOSEN <- aa_modif_all_comb_p20_2129_CHOSEN[aa_modif_all_comb_p20_raw_2129_CHOSEN[, 3] %in% c(1, 3)]
aa_modif_all_comb_p20_raw_2129_CHOSEN <- aa_modif_all_comb_p20_raw_2129_CHOSEN[aa_modif_all_comb_p20_raw_2129_CHOSEN[, 3] %in% c(1, 3), ]

# create 37 new experiments each for 25 models
# create the scripts for creating the input files.
# Creating the R script to create the inputs and job files.
setwd("Seeded_GEMSTAT_ens/")
#

aa <- list()
for(i in 1:length(aa_modif_all_comb_p20_2129_CHOSEN)){
  aa[[i]] <- readLines("Ensemble_based_on_linear_creator_noRole_exp2129.R")
  aa[[i]][9] <- "my_coop_tf_mat <-       rbind(c(1, 1), c(4, 4), c(4, 13), c(10, 10), c(14, 14), c(15, 15), c(17, 17))"
  aa[[i]][10] <- "my_coop_weight_range <- cbind(c(rep(100000, 2), 0.01, rep(100000, 4)), c(rep(1000000, 2), 100, rep(1000000, 4)))"
  aa[[i]][11] <- "my_initial_coop_weight<- c(rep(500000, 2), 1 ,rep(500000, 4))"
  aa[[i]][12] <- 'my_coop_type <- c(rep("DIMER", 2), "SIMPLE", rep("DIMER", 4))'
  aa[[i]][13] <- "my_coop_dist <- c(4, 4, 50 ,4, 6, 6, 2)"
  aa[[i]][14] <- "my_coop_orientation <- rbind(c(1, -1), c(1, -1), c(0, 0), c(1, -1), c(1, 1), c(1, 1), c(1, 1))"
  aa[[i]][15] <- "my_model_index <-  c(52,62,61,22,57,48,19,60,5,37,10,25,29,23,28,70,12,21,53,42,16,9,36,50,4)"
  aa[[i]][19] <- "aa_myTFexpmat <- TF.Exp.Shrinked.microarray_ENTREZ_PerDataset_unique_mat_2_ER01"
  
  for(j in 1:length(aa_modif_all_comb_p20_2129_CHOSEN[[i]])){
    if(aa_modif_all_comb_p20_2129_CHOSEN[[i]][j] < 0){
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[',-aa_modif_all_comb_p20_2129_CHOSEN[[i]][j],', ] <- rep(c(1,0), 45)')
    }else if(aa_modif_all_comb_p20_2129_CHOSEN[[i]][j] < 21){
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[', aa_modif_all_comb_p20_2129_CHOSEN[[i]][j],', ] <- rep(c(0,1), 45)')
    }else{
      aa[[i]][19 + j] <- paste0('aa_myTFexpmat[',aa_modif_all_comb_p20_2129_CHOSEN[[i]][j] - 20,', ] <- rep(1, 90)')
    }
  }
  
  
  aa[[i]][33] <- paste0("GEMSTAT_init_BOlinear(.exp.nu = ", 2859 + i, ",")
  aa[[i]][52] <- "                      TF_expression_mat = aa_myTFexpmat,"
  aa[[i]][54] <- "                       annotation_Filtered_list = GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3$Annot_per_model_perGene[my_model_index],"
  writeLines(aa[[i]], paste0("Ensemble_based_on_linear_creator_noRole_exp",  2859 + i, ".R"))
}
setwd("..")
bash_batch_run_zip_transfer_process(exp_index = c(2860:2896))

aa_WT3 <- list()
aa_WT3_par <- list()
aa_KD3 <- list()
aaKDCOOP3 <- list()
for(i in 2860:2896){
  aa_WT3[[i]] <- GEMSTAT_output_Reader_multiEnh_ensemble(output_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"),
                                                         .fold_change = T)
  
  aa_WT3_par[[i]] <- GEMSTAT_log_Reader_multiEnh(log_files_dir = paste0("Seeded_GEMSTAT_ens/Experiment_", i, "/Outputs"))
  names(aa_WT3)[i] <- paste0("experiment_", i)
  names(aa_WT3_par)[i] <- paste0("experiment_", i)
  #
  print("####################################################################################################")
  print("####################################################################################################")
  print("####################################################################################################")
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print(i)
  print("####################################################################################################")
  print("####################################################################################################")
  print("####################################################################################################")
  # Perform TF KD
  KD_Performer(parent_Exp=i,
               TF_index=list(4),
               shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
               prev_na=5,
               run_GEMSTAT = T)
  aa_KD3[[i]] <- read_KD_results(experiment_nu = i,
                                 shared_dir = "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                                 experiment_output = aa_WT3[[i]])
  
  #Perform coop KD
  KD_Performer_by_params(parent_Exp=i,
                         TF_names = names(TF.motifs.Shrinked.hocomoco_E2F1_added.count),
                         Coop_mat_list=list(c(4,4), c(4, 13)),
                         shared_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                         prev_na=5,
                         run_GEMSTAT=T)
  aaKDCOOP3[[i]] <-  read_KD_results(experiment_nu=i,
                                     shared_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens",
                                     coop_KD=T,
                                     experiment_output = aa_WT3[[i]])
}


load("ERPGRCOOPKD.RData")

aa <- do.call(cbind, lapply(aa_WT3[c(2860:2896)], "[[", 7))
boxplot.matrix(cbind(aa_WT2[[2129]]$Accuracy_All[aa_model_index_2129], aa),
               las= 2, main="Coop_RAP performance")
abline(h=seq(0, 1, 0.01), col = 2, lty = 4)
heatmap.2(aa_modif_all_comb_p20_raw_2129,
          Rowv = F, Colv = F, dendrogram = "none", trace = 'none')
aa_q <- quantile(as.numeric(aa),probs = seq(0, 1, length.out = 11))
aa_dat <- data.frame( x=as.numeric(aa), above_90percentile= as.numeric(aa) >= aa_q[10] )
#library(ggplot2)
qplot(x,data=aa_dat,geom="histogram",fill=above_90percentile)

aaaa <- apply(X = aa, MARGIN = 2,which.max)
aaaa_max <- apply(X = aa, MARGIN = 2,max)
aaaa_max_ind <- sort(aaaa_max, decreasing = T, index.return = T)$ix
aa_modif_all_comb_p20_raw_sort_2860 <- aa_modif_all_comb_p20_raw_2129_CHOSEN[aaaa_max_ind,]
#aa_modif_all_comb_p20_raw_sort <- aa_modif_all_comb_p20_raw_sort[, -3]
colnames(aa_modif_all_comb_p20_raw_sort_2860) <- c("AR", "GR", "PGR", "RARA", "RARG", "RXRA")

heatmap.2(aa_modif_all_comb_p20_raw_sort_2860, 
          Rowv = F, Colv = F, dendrogram = "none",
          trace = 'none', main = "RAP sorted by Performance")

aa_top <- aa_modif_all_comb_p20_raw_sort_2860[1:10,]
aa_top_table <- numeric(6*3)
for(i in 1:ncol(aa_top)){
  aa_tab <- table(aa_top[, i])
  for(j in 1:length(aa_tab)){
    aa_top_table[3 * (i - 1) + as.numeric(names(aa_tab)[j])] <- aa_tab[j]
  }
}

barplot(aa_top_table, 
        col=c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3), rep(6, 3)), 
        names.arg = c("ARp", "ARn", "ARall", "GRp", "GRn","GRall","PGRp", "PGRn",
                      "PGRall","RARAp", "RARAn","RARAall","RARGp", "RARGn","RARGall",
                      "RXRAp", "RXRAn", "RXRAall"), las=2, main = "top100_experiments")
### Filter for experiments that are ER sensitive

aa_WT_acc_2860 <- do.call(cbind, lapply(aa_WT3[c(2860:2896)], "[[", 7))
aa_ER_KD_Acc_2860 <- do.call(cbind, lapply(lapply(X = aa_KD3[c(2860:2896)], "[[", "4"), "[[", 7))
aa_coop_ER_kd_acc_2860 <- do.call(cbind, lapply(lapply(X = aaKDCOOP3[c(2860:2896)], "[[", "ESR1_ESR1"), "[[", 7))
aa_coop_ERPGR_kd_acc_2860 <- do.call(cbind, lapply(lapply(X = aaKDCOOP3[c(2860:2896)], "[[", "ESR1_PGR"), "[[", 7))
#colnames(aa_coop_ER_kd_acc_1398) <- colnames(aa_ER_KD_Acc_1398)
aa_acc_diff_2860 <- aa_WT_acc_2860 - aa_coop_ERPGR_kd_acc_2860

aa_dat <- data.frame( x=as.numeric(aa_acc_diff_2860),
                      above_10= as.numeric(aa_acc_diff_2860) >= 0.1 )
#library(ggplot2)
qplot(x,data=aa_dat,geom="histogram",
      fill=above_10, main = "WT_Accuracy - ER_KD_Accuracy")

# find the ER sensitive models based on the model prediction agreement between WT and KD instead of comparing their accuracy
aa_WT_KD_agreement <- matrix(nrow = length(aa_WT3[[2860]]$Accuracy_All),
                             ncol = 37)
for(i in 1:ncol(aa_WT_KD_agreement)){ # over RAPs
  for(j in 1:nrow(aa_WT_KD_agreement)){ # over models
    aa_wt <- aa_WT3[[2859 + i]]$Prediction_round_list[[j]]
    aa_k  <- aaKDCOOP3[[2859 + i]][[2]]$Prediction_round_list[[j]]
    aa_WT_KD_agreement[j, i] <- sum(aa_wt == aa_k, na.rm = T)/sum(!is.na(aa_wt))
  }
}
rownames(aa_WT_KD_agreement) <- c(1:25)
colnames(aa_WT_KD_agreement) <- c(1:37)


par(mfrow = c(1, 1), mar = c(4,4,4,4))
boxplot.matrix(aa_WT_KD_agreement, las = 2, main = "WT_KD_prediction_agreement_all")



aa_wt_list <- list()
aa_cnt <- 1
aa_mynames <- character(0)
for(i in 1:length(aa_exp_model)){
  for(j in 1:length(aa_exp_model[[i]])){
    aa_wt_list[[aa_cnt]] <- aa_WT2[[as.integer(names(aa_exp_model))[i]]]$Prediction_round_list[[aa_exp_model[[i]][j]]]
    aa_mynames <- c(aa_mynames, paste0(names(aa_exp_model)[i], "_", aa_exp_model[[i]][j]))
    aa_cnt <- aa_cnt + 1
  }
}
names(aa_wt_list) <- aa_mynames

aa_kd_list_list <- list()
aanames <- c(c(1:19), 
             c("4_6", "4_7", "4_11", "4_13", "4_17", "4_18", "4_19"))
for(i in 1:length(aanames)){
  aa_kd_list_list[[i]] <- list()
  aa_cnt <- 1
  aa_mynames <- character(0)
  for(j in 1:length(aa_exp_model)){
    for(k in 1:length(aa_exp_model[[j]])){
      aa_kd_list_list[[i]][[aa_cnt]] <- aa_KD2[[as.integer(names(aa_exp_model))[j]]][[aanames[i]]]$Prediction_round_list[[as.character(aa_exp_model[[j]][k])]]
      aa_mynames <- c(aa_mynames, paste0(names(aa_exp_model)[j], "_", aa_exp_model[[j]][k]))
      aa_cnt <- aa_cnt + 1
      #stopifnot(length(aa_KD2[[as.integer(names(aa_exp_model))[j]]][[as.character(i)]]$Prediction_round_list) == length(aa_exp_model[[j]]))
    }# end of loop over mods
    
  } # end of loop over exps
  names(aa_kd_list_list[[i]]) <- aa_mynames
}
names(aa_kd_list_list) <- c(names(TF.motifs.Shrinked.hocomoco_E2F1_added.count), 
                            c("ER_GATA3", "ER_JUN", "ER_NR5A2", "ER_PGR", "ER_RXRA", "ER_SP1", "ER_TFAP2C"))


library(gplots)
#aa_pr <- aa_WT2[[2139]]$Prediction_round_list
#aa_prkd <- lapply(aa_KD2[[2139]], "[[", 2)
aa <- KD_on_gene_evaluator_multiTF_wrapper(.WT_Prediction_list = aa_wt_list ,
                                           KD_Prediction_list_list = aa_kd_list_list,
                                           .Avg_over_models=T, 
                                           KD_names = names(aa_kd_list_list),
                                           plt = T)



###########################################################################################


########################################################################################################################
######################################                                       ###########################################
######################################          Organize Workplace           ###########################################
######################################                                       ###########################################
########################################################################################################################
# saving and removing unnecessary variables 

save(list = c(c("AssociatedChip.gene.100kb.Interaction",
              "AssociatedChip.gene.100kb.Overlap",
              "AssociatedChip.gene.100kb.RAR.ER", 
              "ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece",
              "ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered",
              "ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene",
              "ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece",
              "ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece_filtered",
              "ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore",
              "ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted",
              "ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th2",
              "ER.associated.reg.elements_uniq_OccuranceMatrixList", 
              "ER.associated.reg.elements_uniq_pair",
              "ER.associated.reg.elements_uniq_pair_Concat_sort_un",
              "ERpickCistrome",
              "ERpickCistromeTimePoint",
              "ERpickCistromeTimePointdf", 
              "Exp_90_278_coop_parameter_combination",
              "GEMSTAT_based_on_linear_exp_annotation_1397",
              "GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3",
              "Gene.ER.Enhancer.annotation.19TFs"
              ,"Gene.ER.Enhancer.annotation.33TFs",
              "Gene.ER.Enhancer.Length", 
              "Genes.Clustered.Based.On.Expression", 
              "GR_ChIP_list", 
              "GR_ChIP_list_GRanges",
              "Greedy_148_restart_enh",
              "Greedy_148_restart_obj", 
              "Greedy_148_restart_par",
              "Greedy_148_restart_results",
              "motif.map", 
              "motif.map.expanded",
              "motifHits.ByGene.RegElement",
              "motifHits.ByGene.RegElement.expanded",
              "Promoter.gene.RAR.ER.Interaction.Overlap.byPromoter",
              "Promoter.gene.RAR.ER.Interaction.Overlap.byPromoter.Mat",
              "ReMapChIP.GRanges.list", 
              "ReMapChIP_expName", 
              "ReMapChIP_FactorName_unique",
              "TimeSerDataMat",
              "TimeSerDataMatNorm"),
              ls(pattern = "DataSet.Exp*"),
              ls(pattern = "GSE*"), 
              ls(pattern = "Optim_Greedy*"),
              ls(pattern = "Sim_Ann_148*"), 
              ls(pattern = "TF.Expression*"),
              ls(pattern = "TF.Index*"),
              ls(pattern = "TF.motifs.TimeSeries*"),
              ls(pattern = "Vicinity200kb*"), 
              ls(pattern = "Vicinity20kb*")
              )
     , file = "Dropped_from_workplace_01_16_2019.RData")

remove_variables_01_16_2019 <- c(c("AssociatedChip.gene.100kb.Interaction",
                                   "AssociatedChip.gene.100kb.Overlap",
                                   "AssociatedChip.gene.100kb.RAR.ER", 
                                   "ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece",
                                   "ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered",
                                   "ER.associated.reg.elements_gte4nzKD_atl1p1n_52_seq_chopped_scoredPerPiece_filtered_perGene",
                                   "ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece",
                                   "ER.associated.reg.elements_gte4nzKD_atl1p1n_52_uniq_seq_chopped_scoredPerPiece_filtered",
                                   "ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_uniq_MotifScore",
                                   "ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted",
                                   "ER.associated.reg.elements_microarray_gte4nzKD_atl1p1n_52_unlisted_monomer_Score_Th2",
                                   "ER.associated.reg.elements_uniq_OccuranceMatrixList", 
                                   "ER.associated.reg.elements_uniq_pair",
                                   "ER.associated.reg.elements_uniq_pair_Concat_sort_un",
                                   "ERpickCistrome",
                                   "ERpickCistromeTimePoint",
                                   "ERpickCistromeTimePointdf", 
                                   "Exp_90_278_coop_parameter_combination",
                                   "GEMSTAT_based_on_linear_exp_annotation_1397",
                                   "GEMSTAT_based_on_linear_exp_annotation_1397_filtered_m3",
                                   "Gene.ER.Enhancer.annotation.19TFs"
                                   , "Gene.ER.Enhancer.annotation.33TFs",
                                   "Gene.ER.Enhancer.Length", 
                                   "Genes.Clustered.Based.On.Expression", 
                                   "GR_ChIP_list", 
                                   "GR_ChIP_list_GRanges",
                                   "Greedy_148_restart_enh",
                                   "Greedy_148_restart_obj", 
                                   "Greedy_148_restart_par",
                                   "Greedy_148_restart_results",
                                   "motif.map", 
                                   "motif.map.expanded",
                                   "motifHits.ByGene.RegElement",
                                   "motifHits.ByGene.RegElement.expanded",
                                   "Promoter.gene.RAR.ER.Interaction.Overlap.byPromoter",
                                   "Promoter.gene.RAR.ER.Interaction.Overlap.byPromoter.Mat",
                                   "ReMapChIP.GRanges.list", 
                                   "ReMapChIP_expName", 
                                   "ReMapChIP_FactorName_unique",
                                   "TimeSerDataMat",
                                   "TimeSerDataMatNorm"),
                                 ls(pattern = "DataSet.Exp*"),
                                 ls(pattern = "GSE*"), 
                                 ls(pattern = "Optim_Greedy*"),
                                 ls(pattern = "Sim_Ann_148*"), 
                                 ls(pattern = "TF.Expression*"),
                                 ls(pattern = "TF.Index*"),
                                 ls(pattern = "TF.motifs.TimeSeries*"),
                                 ls(pattern = "Vicinity200kb*"), 
                                 ls(pattern = "Vicinity20kb*"))

remove(list = remove_variables_01_16_2019)
# variables that were saved and removed: remove_variables_01_16_2019

save(ER.associated.reg.elements_microarray_gte5nzAll_172_uniq_MotifScore, 
     file = "ER.associated.reg.elements_microarray_gte5nzAll_172_uniq_MotifScore.RData")
rm(ER.associated.reg.elements_microarray_gte5nzAll_172_uniq_MotifScore)

save(list = c("GEMSTAT_based_on_linear_exp_KD_Results", 
              "GEMSTAT_based_on_linear_exp_KD_COOP_Results", 
              "GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh"), 
     file = "GEMSTAT_based_on_linear_WT_KD_COOPKD_results.RData")

rm(list = c("GEMSTAT_based_on_linear_exp_KD_Results", 
            "GEMSTAT_based_on_linear_exp_KD_COOP_Results", 
            "GEMSTAT_based_on_linear_exp_results_sigmoid_automat_thresh"))

save.image("ER_Project_workplace_no_bigfile.RData")

save(list = c("my_Descretized_List", "my_Descretized_List_RNAseq"),
     file = "my_Descretized_List_microarray_RNAseq.RData")
rm(list = c("my_Descretized_List", "my_Descretized_List_RNAseq"))

save(list = c("GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3"),
     file = "GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3.RData")
rm(GEMSTAT_based_on_linear_exp_annotation_172_78models_filtered_m3)

remove_Sim_ann_results <- ls(pattern = "Sim_Ann*")
rm(list = ls(pattern = "Sim_Ann*"))

remove_Expkd <- ls(pattern = "Exp.kd*")
rm(list = ls(pattern = "Exp.kd*"))

rm(list = c("AssociatedChip.transcript.200kb",
            "AssociatedChip.transcript.20kb",
            "AssociatedChip.gene.20kb"))

save.image("ER_Project_workplace_no_bigfile_after_BORAW.RData")

########################################################################################################################
########################################################################################################################
#####################################                        ###########################################################
#####################################      bash scripts      ###########################################################
#####################################                        ###########################################################
########################################################################################################################
########################################################################################################################

cat(c("#!/bin/bash\n"), file = "remove_copy_annot_2131_2859.sh", append = F)
for(i in 2131:2859){
  cat(c(paste0("cd ", "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/based_on_linear/Experiment_",
             i, "/Inputs/Filtered_Annotations"),
        "rm *",
        paste0("rm ", "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/based_on_linear/Experiment_",
                       i, "/Outputs/*") 
        
        ),
      file = "remove_copy_annot_2131_2859.sh", append = T, sep = "\n")
  for(j in 1:length(aa_model_index_2129)){
    cat(paste0("cp /shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/based_on_linear/Experiment_2129/Inputs/Filtered_Annotations/annot_", 
               aa_model_index_2129[j],".ann ", "annot_", j,".ann", "\n"), 
        file = "remove_copy_annot_2131_2859.sh", append = T)
  }
}

cat(c("#!/bin/bash\n"), file = "remove_copy_annot_2131_2859_mine.sh", append = F)
for(i in 2131:2859){
  cat(c(paste0("cd ", "/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_",
               i, "/Inputs/Filtered_Annotations"),
        "rm *"
        # paste0("rm ", "/shared-mounts/sinhas/tabebor2/ER_Project/Model/GEMSTAT/based_on_linear/Experiment_",
        #        i, "/Outputs/*") 
        
  ),
  file = "remove_copy_annot_2131_2859_mine.sh", append = T, sep = "\n")
  for(j in 1:length(aa_model_index_2129)){
    cat(paste0("cp /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_2129/Inputs/Filtered_Annotations/annot_", 
               aa_model_index_2129[j],".ann ", "annot_", j,".ann", "\n"), 
        file = "remove_copy_annot_2131_2859_mine.sh", append = T)
  }
}



cat(c("#!/bin/bash\n"), file = "remove_2301_2390.sh", append = F)
for(i in 2301:2390){
  cat(paste0("rm -r ", "Experiment_",
             i,
             "\n"),
      file = "remove_2301_2390.sh", append = T)
}

cat(c("#!/bin/bash\n"), file = "run_script_2_40.run", append = F)
for(i in 2:40){
  cat(paste0("Rscript Ensemble_based_on_linear_creator_noRole_exp",
             89+i,
             ".R\n"),
      file = "run_script_2_40.run", append = T)
}
cat(c("#!/bin/bash\n"), file = "run_script_41_80.run", append = F)
for(i in 41:80){
  cat(paste0("Rscript Ensemble_based_on_linear_creator_noRole_exp",
             89+i,
             ".R\n"),
      file = "run_script_41_80.run", append = T)
}
cat(c("#!/bin/bash\n"), file = "run_script_81_120.run", append = F)
for(i in 81:120){
  cat(paste0("Rscript Ensemble_based_on_linear_creator_noRole_exp",
             89+i,
             ".R\n"),
      file = "run_script_81_120.run", append = T)
}
cat(c("#!/bin/bash\n"), file = "run_script_121_160.run", append = F)
for(i in 121:160){
  cat(paste0("Rscript Ensemble_based_on_linear_creator_noRole_exp",
             89+i,
             ".R\n"),
      file = "run_script_121_160.run", append = T)
}
cat(c("#!/bin/bash\n"), file = "run_script_161_189.run", append = F)
for(i in 161:189){
  cat(paste0("Rscript Ensemble_based_on_linear_creator_noRole_exp",
             89+i,
             ".R\n"),
      file = "run_script_161_189.run", append = T)
}



cat("#!/bin/bash\n", file = "run_283_304.sh", append = F)
for(i in 283:304){
  cat(c(paste0("cd Experiment_", i),
        "./hal_sub_creator_based_on_linear_1enh_perGene.job",
        "chmod +x based_on_linear_1enh_perGene.job.submit", 
        "./based_on_linear_1enh_perGene.job.submit",
        "cd .."), sep = "\n",
      file = "run_283_304.sh",
      append = T)
}
cat("#!/bin/bash\n", file = "run_577_598.sh", append = F)
for(i in 577:598){
  cat(c(paste0("cd Experiment_", i),
        "./hal_sub_creator_based_on_linear_1enh_perGene.job",
        "chmod +x based_on_linear_1enh_perGene.job.submit", 
        "./based_on_linear_1enh_perGene.job.submit",
        "cd .."), sep = "\n",
      file = "run_577_598.sh",
      append = T)
}
cat("#!/bin/bash\n", file = "zip_577_598.sh", append = F)
for(i in 577:598){
  cat(c(paste0("cd Experiment_", i),
        paste0("zip exp", i, "_out.zip", " Outputs/*"),
        "cd .."), sep = "\n",
      file = "zip_577_598.sh",
      append = T)
}
cat("#!/bin/bash\n", file = "transfer_577_598.sh", append = F)
for(i in  577:598){
  cat(c(paste0("rsync Experiment_", i,
               "/exp", i,
               "_out.zip ~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_",
               i,"/")), sep = "\n",
      file = "transfer_577_598.sh",
      append = T)
}
cat("#!/bin/bash\n", file = "process_577_598.sh", append = F)
for(i in  577:598){
  cat(c(paste0("cd ~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_", i, "/"),
        "rm -r Outputs", 
        paste0("unzip exp",i,"_out.zip"),
        "cp ../Experiment_89/log_filter.sh .",
        "./log_filter.sh"), sep = "\n",
      file = "process_577_598.sh",
      append = T)
}
cat("#!/bin/bash\n", file = "zip_194_203.sh", append = F)
for(i in 194:203){
  cat(c(paste0("cd Experiment_", i),
        paste0("zip exp", i, "_out.zip", " Outputs/*"),
        "cd .."), sep = "\n",
      file = "zip_194_203.sh",
      append = T)
}

cat("#!/bin/bash\n", file = "transfer_194_203.sh", append = F)
for(i in 194:203){
  cat(c(paste0("rsync Experiment_", i,
               "/exp", i,
               "_out.zip ~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_",
               i,"/")), sep = "\n",
      file = "transfer_194_203.sh",
      append = T)
}

cat("#!/bin/bash\n", file = "process_194_203.sh", append = F)
for(i in 194:203){
  cat(c(paste0("cd ~/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_", i, "/"),
        "rm -r Outputs", 
        paste0("unzip exp",i,"_out.zip"),
        "cp ../Experiment_89/log_filter.sh .",
        "./log_filter.sh"), sep = "\n",
      file = "process_194_203.sh",
      append = T)
}

# get annotation files
cat("#!/bin/bash\n", file = "Seeded_GEMSTAT_ens/Experiment_54/annotation_job_exp54.sh", append = F)
for(i in 1:148){
  cat(paste0("/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seqannot -s /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_54/Inputs/Sequence/sequences_", 
             i, 
             ".fa  -m /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_54/Inputs/Motifs/motifs.wtmx  >>  /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_54/Annotation/Experiment_54_"
             , i,".annot\n"), file = "Seeded_GEMSTAT_ens/Experiment_54/annotation_job_exp54.sh", append = T)
}

cat("#!/bin/bash\n", file = "Seeded_GEMSTAT_ens/Experiment_52/annotation_job_exp52.sh", append = F)
for(i in 1:148){
  cat(paste0("/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seqannot -s /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_52/Inputs/Sequence/sequences_", 
             i, 
             ".fa  -m /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_52/Inputs/Motifs/motifs.wtmx  >>  /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_52/Annotation/Experiment_52_"
             , i,".annot\n"), file = "Seeded_GEMSTAT_ens/Experiment_52/annotation_job_exp52.sh", append = T)
}
cat("#!/bin/bash\n", file = "Seeded_GEMSTAT_ens/Experiment_52/annotation_job_exp52_newER.sh", append = F)
for(i in 1:148){
  cat(paste0("/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seqannot -s /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_52/Inputs/Sequence/sequences_", 
             i, 
             ".fa  -m /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_52/ER_new_motif.wtmx  >>  /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_52/Annotation_newER/Experiment_52_"
             , i,".annot\n"), file = "Seeded_GEMSTAT_ens/Experiment_52/annotation_job_exp52_newER.sh", append = T)
}
cat("#!/bin/bash\n", file = "Seeded_GEMSTAT_ens/Experiment_52/annotation_job_exp52_newER_half.sh", append = F)
for(i in 1:148){
  cat(paste0("/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seqannot -s /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_52/Inputs/Sequence/sequences_", 
             i, 
             ".fa  -m /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_52/ER_new_motif_half.wtmx  >>  /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_52/Annotation_newER_half/Experiment_52_"
             , i,".annot\n"), file = "Seeded_GEMSTAT_ens/Experiment_52/annotation_job_exp52_newER_half.sh", append = T)
}

cat("#!/bin/bash\n", file = "Seeded_GEMSTAT_ens/Experiment_305/annotation_job_exp305.sh", append = F)
for(i in 1:148){
  cat(paste0("/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seqannot -s /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_305/Inputs/Sequence/sequences_", 
             i, 
             ".fa  -m /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_305/Inputs/Motifs/motifs.wtmx  >>  /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_305/Annotation/Experiment_52_"
             , i,".annot\n"), file = "Seeded_GEMSTAT_ens/Experiment_305/annotation_job_exp305.sh", append = T)
}

cat("#!/bin/bash\n", file = "Seeded_GEMSTAT_ens/Experiment_1397/annotation_job.sh", append = F)
for(i in 1:148){
  cat(paste0("/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seqannot -s /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_1397/Inputs/Sequence/sequences_", 
             i, 
             ".fa  -m /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_1397/Inputs/Motifs/motifs.wtmx  >>  /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_1397/Annotation/Experiment_1397_"
             , i,".annot\n"), file = "Seeded_GEMSTAT_ens/Experiment_1397/annotation_job.sh", append = T)
}

cat("#!/bin/bash\n", file = "Seeded_GEMSTAT_ens/Experiment_2/annotation_job.sh", append = F)
for(i in 1:78){
  cat(paste0("/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seqannot -s /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_2/Inputs/Sequence/sequences_", 
             i, 
             ".fa  -m /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_2/Inputs/Motifs/motifs.wtmx  >>  /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/Seeded_GEMSTAT_ens/Experiment_2/Annotation/Experiment_2_"
             , i,".annot\n"), file = "Seeded_GEMSTAT_ens/Experiment_2/annotation_job.sh", append = T)
}

# create bash to run LLR2Pval for all TFs
cat("#!/bin/bash\n", file = "MotifLLR2Pvalue/LLR2Pval_ERmonomer.sh", append = F)
for(i in 1:length(TF.motifs.Shrinked.halfsites.count)){
  cat(paste0("/shared-mounts/sinhas/lib/myperlfunc/MotifLLR2Pvalue//motif_pvalue2.pl ER_monomers/",
             names(TF.motifs.Shrinked.halfsites.count)[i], ".wtmx ER_monomers\n"),
      file = "MotifLLR2Pvalue/LLR2Pval_ERmonomer.sh", append = T)
}

cat("#!/bin/bash\n", file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/expanded_2pval.sh", append = F)
for(i in 1:length(TF.motifs.Expanded.count)){
  cat(paste0("./motif_pvalue2.pl Motifs_to_pval/",
             names(TF.motifs.Expanded.count)[i], ".wtmx Motifs_to_pval\n"),
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/expanded_2pval.sh",
      append = T)
}

cat("#!/bin/bash\n", file = "MotifLLR2Pvalue/LLR2Pval_Hocomoco.sh", append = F)
for(i in 1:length(TF.motifs.Shrinked.hocomoco.count)){
  cat(paste0("./motif_pvalue2.pl HocoMoco_Shrinked/",
             names(TF.motifs.Shrinked.hocomoco.count)[i], ".wtmx HocoMoco_Shrinked\n"),
      file = "MotifLLR2Pvalue/LLR2Pval_Hocomoco.sh", append = T)
}
cat("#!/bin/bash\n", file = "Seeded_GEMSTAT_ens/Rscript283_304.sh", append = F)
for(i in 283:304){
  cat(paste0("Rscript Ensemble_based_on_linear_creator_noRole_exp", i, ".R\n"),
      file = "Seeded_GEMSTAT_ens/Rscript283_304.sh", append = T)
}
cat("#!/bin/bash\n", file = "Seeded_GEMSTAT_ens/Rscript306_327.sh", append = F)
for(i in 306:327){
  cat(paste0("Rscript Ensemble_based_on_linear_creator_noRole_exp", i, ".R\n"),
      file = "Seeded_GEMSTAT_ens/Rscript306_327.sh", append = T)
}
cat("#!/bin/bash\n", file = "Seeded_GEMSTAT_ens/Rscript577_598.sh", append = F)
for(i in 577:598){
  cat(paste0("Rscript Ensemble_based_on_linear_creator_noRole_exp", i, ".R\n"),
      file = "Seeded_GEMSTAT_ens/Rscript577_598.sh", append = T)
}


for(i in 2101:2128){
  cat(paste0("Experiment_", i), file = "Seeded_GEMSTAT_ens/exp_zip_2101_2128", append = T, sep = "\n")
}

########################################################################################################################
########################################################################################################################
########################################################################################################################

aap1 <- apply(X = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172, MARGIN = 1, FUN = (function(x) sum(x==1, na.rm = T)))
aam1 <- apply(X = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172, MARGIN = 1, FUN = (function(x) sum(x==-1, na.rm = T)))
aa00 <- apply(X = my_CommonDifExpMat_16_ERassoc_gte5nzAll_172, MARGIN = 1, FUN = (function(x) sum(x==0, na.rm = T)))

boxplot.matrix(cbind(aap1, aam1, aa00))
hist(aap1, breaks = 20)
boxplot(aap1[aap1 > 0 & aam1 > 0], breaks = 20)
boxplot(aap1)

boxplot(aam1[aap1 > 0 & aam1 > 0], breaks = 20)
boxplot(aam1)

sum(aam1 == 0 & aap1 > 0)


hist(log(aap1/aam1)[aap1 > 0 & aam1 > 0], breaks = 25)
