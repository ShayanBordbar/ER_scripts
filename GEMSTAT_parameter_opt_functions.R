# GEMSTAT motif optimizer
# Optimizing ER motif in a gradient based way
################################################################################################################################################
###################################################                     ########################################################################
###################################################      Functions      #######################################################################
###################################################                     ########################################################################
################################################################################################################################################


read_output_train_test_GEMSTAT_indiv_onlyPRC <- function(output_file){
  # takes GEMSTAT ouput and returns AUPRC
  output_all <- read.table(file = output_file ,
                           header = T, stringsAsFactors = F)
  output_gt <- output_all[seq(1, nrow(output_all), 2),]
  output_model <- output_all[seq(2, nrow(output_all), 2),]
  my_PRC <- pr.curve(scores.class0 = output_model$X1[output_gt$X1 == 1],
                     scores.class1 = output_model$X1[output_gt$X1 == 0],
                     curve = TRUE)
  return(my_PRC[[2]])
}
################################################################################################################################################
################################################################################################################################################
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
################################################################################################################################################
################################################################################################################################################
pwm_gradient_calc <- function(my_motif_list,
                              motif_index,
                              grad_calc_step, 
                              constrained_rows,
                              inputs_dir, 
                              .GEMSTAT_call, 
                              trained_par_address){
  # my_motif_list: is a list where each entry is a probability based PWM, where each row is a position and columns are: A, C, G, T: each row sums to one
  # motif_index : is the index of the motif to be modified
  # grad_calc_step: is the step size for the change used to compute gradient
  # constrained_rows: is a logical vector, True for rows that optimization should be constrained and F otherwise.
  #  constrain is to keep the consensus for that position unchanged.
  #  # consens: is a character vector containing the consensus sequence for the motif
  # inputs_dir: is the directory where GEMSTAT input files are located (it should end with /)
  # GEMSTAT_call is the command used to run GEMSTAT for evaluation: it shouldn't have the motif part and output file address part as well as the parameter part
  
  stopifnot(ncol(my_motif_list[[motif_index]]) == 4)
  new_motif_list <- list()
  pos_cnt <- 0
  
  aan<- unlist(strsplit(trained_par_address, split = "\\/"))
  aan <- aan[length(aan)]
  model_name <- unlist(strsplit(aan, split = "\\."))[1]
  if(!dir.exists(paste0("Motif_opt_",model_name))){
    dir.create(paste0("Motif_opt_",model_name))
  }
  
  my_gradient <- matrix(nrow = nrow(my_motif_list[[motif_index]]), 
                        ncol = ncol(my_motif_list[[motif_index]]))
  my_motif <- my_motif_list[[motif_index]]
  MotifWriter(motif.List = my_motif_list, pseudo = 0.001,
              output.File.Name = paste0(inputs_dir, paste0("Motif_opt_",model_name, "/"), "base_motif"))
  base_perf <- updated_motif_eval(GEMSTAT_call = .GEMSTAT_call,
                                  updated_motif_address = paste0(inputs_dir, paste0("Motif_opt_",model_name, "/"),  "base_motif.wtmx"), 
                                  trained_par_add =trained_par_address,
                                  output_file_address = paste0(inputs_dir, paste0("Motif_opt_",model_name, "/"), "out_base"))
  print(paste("base performance is :", base_perf))
  for(cur_pos in 1:nrow(my_motif)){
    cur_cns <- which.max(my_motif[cur_pos,])
    for(cur_base in 1:4){
      pos_cnt <- pos_cnt + 1
      new_motif_list[[pos_cnt]] <- my_motif
      new_motif_list[[pos_cnt]][cur_pos, cur_base] <- ifelse(((cur_cns == cur_base) | (!constrained_rows[cur_pos])), 
                                                             max(0.001, my_motif[cur_pos, cur_base] + grad_calc_step), 
                                                             max(0.001, min(my_motif[cur_pos, cur_cns],
                                                                            my_motif[cur_pos, cur_base] + grad_calc_step)))
      new_motif_list[[pos_cnt]][cur_pos, ] <- new_motif_list[[pos_cnt]][cur_pos, ]/sum(new_motif_list[[pos_cnt]][cur_pos, ])
      #print(paste("new motif ", cur_pos, "_", cur_base))
      #print(new_motif_list[[pos_cnt]])
      new_mot_list <- my_motif_list
      new_mot_list[[motif_index]] <- new_motif_list[[pos_cnt]]
      MotifWriter(motif.List = new_mot_list, pseudo = 0.001,
                  output.File.Name = paste0(inputs_dir, paste0("Motif_opt_",model_name, "/"), "grad_", cur_pos, "_", cur_base))
      my_gradient[cur_pos, cur_base] <- updated_motif_eval(GEMSTAT_call = .GEMSTAT_call,
                                                           trained_par_add =trained_par_address,
                                                           updated_motif_address = paste0(inputs_dir, paste0("Motif_opt_",model_name, "/"), "grad_", cur_pos, "_", cur_base, ".wtmx"), 
                                                           output_file_address = paste0(inputs_dir, paste0("Motif_opt_",model_name, "/"), "out_", cur_pos, "_", cur_base))
      # print(paste0("gradient for ",  cur_pos, "_", cur_base))
      # print(my_gradient[cur_pos, cur_base] - base_perf) 
    }
  }
  my_gradient <- my_gradient - base_perf
  # print("full gradient mat")
  # print(my_gradient)
  return(my_gradient)
}
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
CountToProbPWM <- function(motif){
  stopifnot(ncol(motif) == 4)
  newMotif <- motif/rowSums(motif)
  return(newMotif)
}
################################################################################################################################################
################################################################################################################################################

#aa_newmot <- lapply(TF.motifs.Expanded_pseudo_count, CountToProbPWM)
################################################################################################################################################
pwm_update <- function(my_motif,
                       gradient_mat,
                       grad_move_step, 
                       constrained_rows){
  # my_motif: is a probability based PWM, where each row is a position and columns are: A, C, G, T: each row sums to one
  # gradient_mat : is a matrix with the same dimensions as pwm, indicates the gradient of change
  # grad_move_step: is the step size for the change in gradient direction
  # constrained_rows: is a logical vector, True for rows that optimization should be constrained and F otherwise.
  #  constrain is to keep the consensus for that position unchanged.
  my_motif_updated <- my_motif
  for(cur_pos in 1:nrow(my_motif)){
    cur_cns <- which.max(my_motif[cur_pos,])
    cns_new <- max(0.0001, min(1, my_motif[cur_pos, cur_cns] + 0.01 * grad_move_step * gradient_mat[cur_pos, cur_cns]))
    my_motif_updated[cur_pos, cur_cns] <- cns_new
    for(cur_base in setdiff(c(1:4), cur_cns)){
      if(constrained_rows[cur_pos]){
        base_new <- max(0.0001, min((my_motif_updated[cur_pos, cur_cns] - 0.01),(my_motif[cur_pos, cur_base] + 0.01 * grad_move_step * gradient_mat[cur_pos, cur_base])))
      }else{
        base_new <- max(0.0001, min(1, (my_motif[cur_pos, cur_base] + 0.01 * grad_move_step * gradient_mat[cur_pos, cur_base])))
      }
      my_motif_updated[cur_pos, cur_base] <- base_new
    }
    my_motif_updated[cur_pos,] <- my_motif_updated[cur_pos,]/sum(my_motif_updated[cur_pos,])
  }
  return(my_motif_updated)
}

################################################################################################################################################
################################################################################################################################################

updated_motif_eval <- function(GEMSTAT_call, 
                               trained_par_add,
                               updated_motif_address, 
                               output_file_address){
  # GEMSTAT_call is the command used to run GEMSTAT for evaluation: it shouldn't have the motif part and output file address part
  # updated_motif_address is a character containing the address of the updated motif file
  # trained_par_add: is the address of the trained parameter to be used
  # it runs GEMSTAT using the new motif and evaluates AUPRC
  
  system(paste0(GEMSTAT_call," -p ", trained_par_add ," -m ", updated_motif_address, " -fo ", output_file_address))
  my_PRC <- read_output_train_test_GEMSTAT_indiv_onlyPRC(output_file = output_file_address)
  print("current PRC")
  print(my_PRC)
  return(my_PRC)
}
################################################################################################################################################
################################################################################################################################################
# 
# random_batch_constructor <- function(train_size, batch_size=130, nu_batches=100){
#   # given the size of positive (pos_size) and negative sets (neg_size), and the
#   #  number of positive (pos_nu) and negatives(neg_nu) in each batch, creates a list of
#   #  length (nu_batches) batches of training indexes. each entry of the list is a list
#   #  containing two integer vectors. first and second entries are indices of examples in
#   #  the positive and negative sets, respectively.
#   shuff_num <- ceiling(train_size/nu_batches)
#   
#   for(cur_epoch in 1:nu_epochs){
#     train_ind <- sample(x = c(1:train_size), size = train_size, replace = F)
#     
#     
#   }
# }
################################################################################################################################################
################################################################################################################################################

GEMSTAT_motif_optimizer <- function(initial_motif_list,
                                    .motif_index, 
                                    nu_steps,
                                    ..GEMSTAT_call,
                                    .grad_calc_step,
                                    .grad_move_step,
                                    .constrained_rows,
                                    .inputs_dir, 
                                    .trained_par_address){
  library(PRROC)
  working_mot_list <- initial_motif_list
  motif_holder_list <- list()
  for(cur_step in 1:nu_steps){
    print(paste("step", cur_step))
    print("computing gradient")
    cur_grad <- pwm_gradient_calc(my_motif_list = working_mot_list,
                                  motif_index = .motif_index,
                                  grad_calc_step = .grad_calc_step,
                                  constrained_rows = .constrained_rows,
                                  inputs_dir = .inputs_dir, 
                                  .GEMSTAT_call = ..GEMSTAT_call, 
                                  trained_par_address = .trained_par_address)
    working_mot_list[[.motif_index]] <- pwm_update(my_motif = working_mot_list[[.motif_index]], 
                                                   gradient_mat = cur_grad, 
                                                   grad_move_step = .grad_move_step, 
                                                   constrained_rows = .constrained_rows)
    motif_holder_list[[cur_step]] <- working_mot_list[[.motif_index]]
    
  }
  return(list(last_set=working_mot_list, all_steps = motif_holder_list))
}
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
motif_optimizer_seq_label_weight_creator <- function(nu_steps,
                                                     batch_size_motif,
                                                     batch_size_GEMSTAT,
                                                     training_seq_pool,
                                                     training_seq_label){
  # nu_steps : is an integer indicating the number of steps for Optimization process
  # batch_size_motif: is an integer indicating the size of the training set in motif optimization
  # batch_size_GEMSTAT : is an integer indicating the size of the training set in GEMSTAT optimization
  # training_seq_pool: character vector containing the sequences in training set
  # training_seq_label: integer vector indicating the labels: each entry is either 0 or 1
  library(ShortRead)
  library(Biostrings)
  stopifnot(length(names(training_seq_pool)) == length(training_seq_pool), 
            length(training_seq_pool) == length(training_seq_label),
            batch_size_motif <= length(training_seq_label), 
            batch_size_GEMSTAT <= length(training_seq_label), 
            all(training_seq_label %in% c(0, 1)))
  
  if(!dir.exists("Motif_opt_input")){
    dir.create("Motif_opt_input")
  }
  if(!dir.exists("Motif_opt_output")){
    dir.create("Motif_opt_output")
  }
  total_seq_nu <- batch_size_motif * nu_steps
  total_seq_nu_GEM <- batch_size_GEMSTAT * nu_steps
  training_seq_nu <- length(training_seq_pool)
  nu_shuff <- ceiling(total_seq_nu/training_seq_nu)
  nu_shuff_GEMS <- ceiling(total_seq_nu_GEM/training_seq_nu)
  all_seq <- character(0)
  all_label <- integer(0)
  for(c_s in 1:nu_shuff){ # creating the seq set for motif optimization
    cur_samp <- sample(x = c(1:training_seq_nu), 
                       size = training_seq_nu, 
                       replace = F)
    all_seq <- c(all_seq, training_seq_pool[cur_samp])
    all_label <- c(all_label, training_seq_label[cur_samp])
  }
  
  all_seq_GEMS <- character(0)
  all_label_GEMS <- integer(0)
  for(c_s in 1:nu_shuff_GEMS){ # creating the seq set for motif optimization
    cur_samp <- sample(x = c(1:training_seq_nu), 
                       size = training_seq_nu, 
                       replace = F)
    all_seq_GEMS <- c(all_seq_GEMS, training_seq_pool[cur_samp])
    all_label_GEMS <- c(all_label_GEMS, training_seq_label[cur_samp])
  }
  
  ind_cnt_motif <- 1
  ind_cnt_GEMS <- 1
  seq_list <- list()
  label_list <- list()
  seq_list_GEMS <- list()
  label_list_GEMS <- list()
  for(cur_step in 1:nu_steps){
    seq_list[[cur_step]] <- all_seq[(ind_cnt_motif):(ind_cnt_motif + batch_size_motif - 1)]
    label_list[[cur_step]] <- all_label[(ind_cnt_motif):(ind_cnt_motif + batch_size_motif - 1)]
    
    seq_list_GEMS[[cur_step]] <- all_seq_GEMS[(ind_cnt_GEMS):(ind_cnt_GEMS + batch_size_GEMSTAT - 1)]
    label_list_GEMS[[cur_step]] <- all_label_GEMS[(ind_cnt_GEMS):(ind_cnt_GEMS + batch_size_GEMSTAT - 1)]
    
    # write sequence
    writeFasta(DNAStringSet(seq_list[[cur_step]]), width = 1000,
               file = paste0("Motif_opt_input/seq_step_",cur_step,".fa"))
    writeFasta(DNAStringSet(seq_list_GEMS[[cur_step]]), width = 1000,
               file = paste0("Motif_opt_input/seq_GEMS_step_",cur_step,".fa"))
    # write label
    aalab1 <- as.numeric(label_list[[cur_step]])
    aalab1 <- matrix(aalab1, nrow = length(aalab1))
    rownames(aalab1) <- names(seq_list[[cur_step]])
    cat(c("Rows", paste0("1", "\n")), sep = "\t",
        append = F, file = paste0("Motif_opt_input/Label_step_", cur_step, ".tab"))
    for(i in 1:length(aalab1)){
      cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                      rep("\n", as.integer(i != length(aalab1)))), 
          append = T,
          file = paste0("Motif_opt_input/Label_step_", cur_step, ".tab"), sep = "\t")
    }
    
    aalab1 <- as.numeric(label_list_GEMS[[cur_step]])
    aalab1 <- matrix(aalab1, nrow = length(aalab1))
    rownames(aalab1) <- names(seq_list_GEMS[[cur_step]])
    cat(c("Rows", paste0("1", "\n")), sep = "\t",
        append = F, file = paste0("Motif_opt_input/Label_GEMS_step_", cur_step, ".tab"))
    for(i in 1:length(aalab1)){
      cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                      rep("\n", as.integer(i != length(aalab1)))), 
          append = T,
          file = paste0("Motif_opt_input/Label_GEMS_step_", cur_step, ".tab"), sep = "\t" )
    }
    
    # write weights --> only for GEMSTAT training
    aalab1 <- as.numeric(label_list_GEMS[[cur_step]])
    aalab1 <- matrix(aalab1, nrow = length(aalab1))
    rownames(aalab1) <- names(seq_list_GEMS[[cur_step]])
    zero_cnt <- sum(aalab1 == 0)
    one_cnt <-  sum(aalab1 == 1)
    aalab1[aalab1 == 0] <-  one_cnt/(zero_cnt + one_cnt)
    aalab1[aalab1 == 1] <- zero_cnt/(zero_cnt + one_cnt)
    cat(c("Rows", paste0("1", "\n")), sep = "\t",
        append = F, file = paste0("Motif_opt_input/Weight_GEMS_step_", cur_step, ".tab"))
    for(i in 1:length(aalab1)){
      cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                      rep("\n", as.integer(i != length(aalab1)))), 
          append = T,
          file = paste0("Motif_opt_input/Weight_GEMS_step_", cur_step, ".tab"), sep = "\t" )
    }
    
    ind_cnt_motif <- ind_cnt_motif + batch_size_motif
    ind_cnt_GEMS <- ind_cnt_GEMS + batch_size_GEMSTAT
  }
}

################################################################################################################################################
################################################################################################################################################
# example
# setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8")
# aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[4]]), split = "_")), "[[", 1))
# aalab1[aalab1 == "pos"] <- 1
# aalab1[aalab1 == "neg"] <- 0
# aalab1 <- as.numeric(aalab1)
# motif_optimizer_seq_label_weight_creator(nu_steps = 3,
#                                          batch_size_motif = 200,
#                                          batch_size_GEMSTAT = 400,
#                                          training_seq_pool = GEMSTAT_Ensemble_train_SeqList[[4]],
#                                          training_seq_label = aalab1)
# 
# aa_f1 <- list.files( path = "Motif_opt_input/", pattern = "seq_GEMS_step_*", full.names = T)
# aa_f2 <- list.files( path = "Motif_opt_input/", pattern = "seq_step_*", full.names = T)
# 
# aa_seq1 <- list()
# aa_seq2 <- list()
# 
# for(i in 1:3){
#   aa_seq1[[i]] <- readLines(aa_f1[i])
#   aa_seq2[[i]] <- readLines(aa_f2[i])
# }
# for(i in 1:3){
#   print("##########################")
#   aaa <- aa_seq1[[i]][seq(1, length(aa_seq1[[i]]), 2)]
#   aaa <- unlist(lapply(strsplit(aaa, split = "_"), "[[", 1))
#   print(table(aaa))
#   print("#############")
#   aaa <- aa_seq2[[i]][seq(1, length(aa_seq2[[i]]), 2)]
#   aaa <- unlist(lapply(strsplit(aaa, split = "_"), "[[", 1))
#   print(table(aaa))
#   print("##########################")
# }



################################################################################################################################################
################################################################################################################################################
pwm_gradient_motif_writer <- function(my_motif_list_RData_address,
                                      motif_index,
                                      grad_calc_step, 
                                      constrained_rows=logical(0),
                                      model_number_name,
                                      step_nu,
                                      main_directory ="/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/"){
  # my_motif_list_RData_address: is the address of RData file containing an object called:
  #  original_motif_list which is a list where each entry is a probability based PWM, where each row is a position and columns are: A, C, G, T: each row sums to one
  # motif_index : is the index of the motif to be modified
  # grad_calc_step: is the step size for the change used to compute gradient
  # constrained_rows: is a logical vector, True for rows that optimization should be constrained and F otherwise.
  #  constrain is to keep the consensus for that position unchanged.
  #  # consens: is a character vector containing the consensus sequence for the motif
  # model_number_name: is a character or integer identifying the model, it will be used in creating file names
  # inputs_dir: is the directory where GEMSTAT input files are located (it should end with /)
  # trained_par_address : is the address of the trained par to be used
  # .GEMSTAT_call: is the command to call GEMSTAT, without containing the -p -m -fo -po flags: these will be added afterwards
  # step_nu : integer indicating the step of the optimization: to include in job names
  load(my_motif_list_RData_address)
  my_motif_list <- original_motif_list
  stopifnot(ncol(my_motif_list[[motif_index]]) == 4)
  new_motif_list <- list()
  
  if( step_nu == 1){
    if(!dir.exists(paste0(main_directory, "Motif_opt_grad"))){
      dir.create(paste0(main_directory, "Motif_opt_grad"))
    }
  }
  
  if(length(constrained_rows) == 0){
    constrained_rows <- rep(T, nrow(my_motif_list[[motif_index]]))
  }
  
  # aan <- unlist(strsplit(trained_par_address, split = "\\/"))
  # aan <- aan[length(aan)]
  # model_name <- unlist(strsplit(aan, split = "\\."))[1]
  model_name <- model_number_name
  
  my_motif <- my_motif_list[[motif_index]]
  
  if(step_nu == 1){ # otherwise this file is written by the function that updates motifs
    MotifWriter(motif.List = my_motif_list, pseudo = 0.001,
                output.File.Name = paste0(main_directory, "Motif_opt_grad/", "base_motif_", model_name))
  }
  
  
  pos_cnt <- 0
  for(cur_pos in 1:nrow(my_motif)){
    cur_cns <- which.max(my_motif[cur_pos,])
    for(cur_base in 1:4){
      pos_cnt <- pos_cnt + 1
      new_motif_list[[pos_cnt]] <- my_motif
      
      if(cur_cns == cur_base){
        if(constrained_rows[cur_pos]){
          new_base <- max(0.0001, min((my_motif[cur_pos, cur_cns] - 0.01), (my_motif[cur_pos, cur_base] + grad_calc_step)))
        }else{
          new_base <-  max(0.0001, (my_motif[cur_pos, cur_base] + grad_calc_step))
        }
        
      }else{
        new_base <-  max(0.0001, (my_motif[cur_pos, cur_base] + grad_calc_step))
      }
      
      new_motif_list[[pos_cnt]][cur_pos, cur_base] <- new_base 
      new_motif_list[[pos_cnt]][cur_pos, ] <- new_motif_list[[pos_cnt]][cur_pos, ]/sum(new_motif_list[[pos_cnt]][cur_pos, ])
      new_mot_list <- my_motif_list
      new_mot_list[[motif_index]] <- new_motif_list[[pos_cnt]]
      MotifWriter(motif.List = new_mot_list, pseudo = 0.001,
                  output.File.Name = paste0("/scratch/", 
                                            "grad_", model_name, "_",
                                            cur_pos, "_", cur_base))
      
    }
  }
  system(paste0("cp /scratch/grad_", model_name, "_* ", main_directory ,"Motif_opt_grad/"))
  
}
################################################################################################################################################
norm_vec <- function(x) sqrt(sum(x^2))
# create a script that calls this function at every step, with the right arguments
################################################################################################################################################
pwm_gradient_calculator_motif_update <- function(my_motif_list_RData_address,
                                                 my_par_address,
                                                 my_weight_address,
                                                 motif_index,
                                                 grad_move_step,
                                                 grad_calc_step,
                                                 constrained_rows=logical(0), 
                                                 model_number_name, 
                                                 step_nu,
                                                 main_directory){
  # step_nu : is an integer indicating the step number in order to identify the right motif list , as well as right output files to read in
  # my_motif_list_RData_address: is the address of RData file containing the motifs to work on: the motif is named: original_motif_list
  # model_number_name: is a character or integer identifying the model, it will be used in creating file names
  # my_par_address: address of the parameter file to read parameters of the logistic regression
  # my_weight_address : address of the weight file for training. using the weights to compute the obj func
  # main_directory: is the directory where GEMSTAT input files are located (it should end with /) eg "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/"
  
  library(PRROC)
  library(rjson)
  load(my_motif_list_RData_address)
  my_motif_list <- original_motif_list
  # aan <- unlist(strsplit(trained_par_address, split = "\\/"))
  # aan <- aan[length(aan)]
  # model_name <- unlist(strsplit(aan, split = "\\."))[1]
  model_name <- model_number_name
  
  if(length(constrained_rows) == 0){
    constrained_rows <- rep(T, nrow(my_motif_list[[motif_index]]))
  }
  
  system(paste0("cp ", main_directory, "Motif_opt_output/base_out_", model_name, " /scratch/"))
  system(paste0("cp ", main_directory, "Motif_opt_output/grad_out_", model_name, "*", " /scratch/"))
  #print("copied to scratch.")
  
  cur_par_file <- fromJSON(file = my_par_address)
  my_bias <- cur_par_file$log_Reg[[1]]$bias
  my_coeff <- cur_par_file$log_Reg[[1]]$coeff
  my_beta <- cur_par_file$enh[[1]]$beta
  my_weights <- read.table(my_weight_address , header = T, sep = "\t", stringsAsFactors = F)
  
  base_perf <- read_output_train_test_GEMSTAT_indiv_logistic(output_file = paste0("/scratch/base_out_", model_name),
                                                             bias_par = my_bias,
                                                             coeff_par = my_coeff,
                                                             weights = my_weights,
                                                             beta = my_beta)
  system(paste0("rm /scratch/base_out_", model_name))
  #print("read base output from scratch.")
  out_files <- list.files(path = "/scratch/",
                          pattern = paste0("grad_out_", model_name, "*"))
  # print("out_files")
  # print(out_files)
  my_out_split <- strsplit(out_files, split = "_")
  my_ln <- length(my_out_split[[1]])
  out_row <- as.integer(unlist(lapply(strsplit(out_files, split = "_"), "[[", (my_ln-1))))
  out_bas <- as.integer(unlist(lapply(strsplit(out_files, split = "_"), "[[", my_ln)))
  PRC_holder <- numeric(length = length(out_files))
  gradient_holder <- matrix(nrow = nrow(my_motif_list[[motif_index]]), ncol = 4)
  for(cur_out in 1:length(out_files)){
    cur_perf <- read_output_train_test_GEMSTAT_indiv_logistic(output_file = paste0("/scratch/", out_files[cur_out]),
                                                              bias_par = my_bias, 
                                                              coeff_par = my_coeff,
                                                              weights = my_weights,
                                                              beta = my_beta)
    gradient_holder[out_row[cur_out], out_bas[cur_out]] <- -(cur_perf - base_perf)/grad_calc_step
    system(paste0("rm /scratch/", out_files[cur_out]))
  }
  print("gradient_holder before normalization")
  print(gradient_holder)
  gradient_holder <- gradient_holder/norm_vec(c(gradient_holder))
  print("gradient_holder after normalization")
  print(gradient_holder)
  stopifnot(sum(is.na(gradient_holder)) == 0)
  #print("read other outputs from scratch.")
  print("motif_before_update")
  print(my_motif_list[[motif_index]])
  updated_motif <- pwm_update(my_motif = my_motif_list[[motif_index]],
                              gradient_mat = gradient_holder,
                              grad_move_step = grad_move_step, 
                              constrained_rows = constrained_rows)
  print("motif_after_update")
  print(updated_motif)
  working_mot_list <- my_motif_list
  working_mot_list[[motif_index]] <- updated_motif
  MotifWriter(motif.List = working_mot_list, pseudo = 0.001,
              output.File.Name = paste0(main_directory,"Motif_opt_grad/", "base_motif_", model_name))
  original_motif_list  <- working_mot_list
  save(list = c("original_motif_list"), 
       file = paste0(main_directory,"Motif_opt_grad/base_motif_", model_name, "_", step_nu, ".RData"))
}
################################################################################################################################################
################################################################################################################################################
read_output_train_test_GEMSTAT_indiv_logistic <- function(output_file, bias_par, coeff_par, weights, beta=numeric(0)){
  # weights : a matrix with two columns, first column is the sequence name, second column is the weight
  # output_file : address of the output file
  # bias_par: a numeric used in logistic regression
  # coeff_par: a numeric used in logistic regression
  # normalize_output_01 : if True first brings the output between 0 and 1 and then calculates the loss
  # beta: if provided it first divides predictions by beta and then calculates the loss
  output_all <- read.table(file = output_file ,
                           header = T, stringsAsFactors = F)
  output_gt <- output_all[seq(1, nrow(output_all), 2),]
  output_model <- output_all[seq(2, nrow(output_all), 2),]
  output_gt[,2] <- 2*output_gt[,2] - 1 
  # if(normalize_output_01){
  #   output_model[, 2] <- output_model[, 2]/max(abs(output_model[, 2]))
  # }
  if(length(beta) == 1){
    output_model[, 2] <- output_model[, 2]/beta
  }
  # print(nrow(output_model))
  # print(nrow(weights))
  # stopifnot(nrow(output_model) == nrow(weights),
  #           all(output_model[, 1] == weights[, 1]))
  # quick fix just for now
  #
  aa_lab <- unlist(lapply(strsplit(output_model$Rows, split = "_"), "[[", 1))
  aa_lab_table <- table(aa_lab)
  my_weights <- matrix(nrow = nrow(output_model), ncol = 2)
  my_weights[, 1] <- output_model[, 1]
  for(clab in 1:length(aa_lab_table)){
    my_weights[aa_lab %in% names(aa_lab_table)[clab], 2] <- as.numeric(aa_lab_table[setdiff(names(aa_lab_table), names(aa_lab_table)[clab])])/sum(aa_lab_table)
  }
  #print(my_weights)
  weights <- my_weights
  
  aa_error <- 0
  for(cur_enh in 1:nrow(output_model)){
    cur_LL = as.numeric(weights[cur_enh, 2]) * log(1 + exp( -1.0 * output_gt[cur_enh, 2] * coeff_par * (output_model[cur_enh, 2] - bias_par) ))
    aa_error <- aa_error + cur_LL
  }
  return(aa_error)
}

#as.numeric(0.2025 * log(1 + exp( 18.1284 * (76.161 - 0.98) )))

################################################################################################################################################
# example
# aa <- fromJSON(file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Trained_par/par_3203_1.txt.Filter")
# aa2 <- read.table("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Out/output_par_1_1.txt" , header = T, sep = "\t", stringsAsFactors = F)
# read_output_train_test_GEMSTAT_indiv_logistic(output_file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Out/output_par_3203_1.txt",
#                                               bias_par = aa$log_Reg[[1]]$bias, 
#                                               coeff_par = aa$log_Reg[[1]]$coeff,
#                                               weights = aa2)
# aa22 <- aa2[seq(2, nrow(aa2), 2),]
# aax <- (table(unlist(lapply(strsplit(aa22$Rows, split = "_"), "[[", 1))))
################################################################################################################################################
parameter_hash <- function(binding_mat, alpha_mat, coop_mat,
                           bind_ens, alpha_ens, coop_ens, 
                           binding_range, alpha_range, coop_range,
                           binding_log, alpha_log, coop_log){
  # binding_mat : is a matrix where each row is a model, each column is binding parameter, order is the same as bind_ens
  # alpha_mat: same as binding_mat but for alphas
  # coop_mat: same as binding_mat but for coops
  # bind_ens: is an integer vector, same length ans ncol of binding_mat, showing whether the space was divided based on this parameter?
  # alpha_ens: same as bind_ens but for alphas
  # coop_ens: same as bind_ens but for coops
  # binding_range: is a matrix with nrow=2 and ncol==ncol(binding_mat), 
  #  first row indicates the lower bound for this param, second row shows lower bound on binding params
  # alpha_range: same as binding_range but for alphas
  # coop_range: same as binding_range but for coops
  # binding_log: logical, if True it indicates that the space division was done in log10 scale for binding params
  # alpha_log: logical, if True it indicates that the space division was done in log10 scale for alpha params
  # coop_log: logical, if True it indicates that the space division was done in log10 scale for coop params
  
  
  stopifnot(ncol(binding_mat) == length(bind_ens), 
            ncol(alpha_mat) == length(alpha_ens),
            ncol(coop_mat) == length(coop_ens))
  
  if(binding_log){
    my_binding_mat <- log10(binding_mat)
    break_points_binding <- colMeans(log10(binding_range))
  }else{
    my_binding_mat <- binding_mat
    break_points_binding <- colMeans(binding_range)
  }
  if(alpha_log){
    my_alpha_mat <- log10(alpha_mat)
    break_points_alpha <- colMeans(log10(alpha_range))
  }else{
    my_alpha_mat <- alpha_mat
    break_points_alpha <- colMeans(alpha_range)
  }
  if(coop_log){
    my_coop_mat <- log10(coop_mat)
    break_points_coop <- colMeans(log10(coop_range))
  }else{
    my_coop_mat <- coop_mat
    break_points_coop <- colMeans(coop_range)
  }
  all_pars <- cbind(cbind(my_binding_mat, my_alpha_mat), my_coop_mat)
  all_points<- c(bind_ens, alpha_ens, coop_ens)
  
  
  print(break_points_binding)
  print(break_points_alpha)
  print(break_points_coop)
  
  break_points <- c(break_points_binding,
                    break_points_alpha, 
                    break_points_coop)
  
  filtered_par_mat <- all_pars[, all_points==1]
  filtered_break_poitns <- break_points[all_points==1]
  bit_matrix <- matrix(nrow = nrow(filtered_par_mat), ncol = ncol(filtered_par_mat))
  for(c_mod in 1:nrow(filtered_par_mat)){
    bit_matrix[c_mod, ] <- as.integer(filtered_par_mat[c_mod, ] < filtered_break_poitns)
  }
  return(bit_matrix)
}
################################################################################################################################################
# # example
# aa_bind_ff_ens <- c(0,0,1,0,0,0,1,0,1,0,0,0,0,1,1,0,1,0)
# aa_alpha_ff_ens<- c(0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,0)
# 
# # choose three distinct models
# aa_perf_srt <- names(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_PRC)[sort(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_PRC,
#                                                                          decreasing = T, index.return=T)$ix]
# aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
# aa_perf_srt22 <- paste0("log_", aa_perf_srt2)
# #E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$binding[aa_perf_srt22[1:100],]
# 
# aa_p_hash <- parameter_hash(binding_mat = E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$binding[aa_perf_srt22[1:100],],
#                             alpha_mat = E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$alpha[aa_perf_srt22[1:100],],
#                             coop_mat = E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$coop[aa_perf_srt22[1:100],],
#                             bind_ens = aa_bind_ff_ens,
#                             alpha_ens = aa_alpha_ff_ens, 
#                             coop_ens = rep(0, 10), 
#                             binding_range = rbind(rep(0.01, 18), rep(250.5, 18)),
#                             alpha_range = rbind(rep(1e-6, 18), rep(100, 18)),
#                             coop_range = rbind(rep(0.01, 10), rep(100.5, 10)),
#                             binding_log = F, 
#                             alpha_log = F, 
#                             coop_log = F)
# 
# aa_dup_hash <- duplicated(aa_p_hash)
# aa_chosen_model <- character(0)
# aacnt <- 1
# aa_cur <- 1
# while(aacnt <= 5){
#   if(!aa_dup_hash[aa_cur]){
#     aa_chosen_model <- c(aa_chosen_model, rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$binding[aa_perf_srt22[1:100],])[aa_cur])
#     aacnt <- aacnt + 1
#   }
#   aa_cur <- aa_cur + 1
#   
# }
# aa_chosen_model

################################################################################################################################################
motif_optimizer_job_writer <- function(nu_steps,
                                       trained_par_address,
                                       model_number_name,
                                       motif_length,
                                       grad_move_step,
                                       grad_calc_step,
                                       GEMSTAT_call,
                                       motif_index,
                                       annotation_fixed = F,
                                       main_directory ="/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/"){
  # .GEMSTAT_call: is the command to call GEMSTAT, without containing the -s -e  -train_weights  -na -p -m -fo -po flags: these will be added afterwards
  # motif_length : is the length of the motif being optimized
  # trained_par_address : is the address from the main directory to the trained parameters of the model
  # model_number_name: is a character or integer identifying the model, it will be used in creating file names
  # grad_calc_step: is the step size for the change used to compute gradient
  # grad_move_step: is the step size for the change in gradient direction
  # nu_steps : integer, number of steps to run the optimization
  # motif_index : is the index of the motif to be opimized within the motif files
  # annotation_fixed: if True the annotation for the training sequences should be available precomputed for the GEM runs to use: it's in annotation folder and is named the same as the sequences
  # main_directory: character, indicating the address of the experiment.
  
  
  # aan <- unlist(strsplit(trained_par_address, split = "\\/"))
  # aan <- aan[length(aan)]
  # model_name <- unlist(strsplit(aan, split = "\\."))[1]
  model_name <- model_number_name
  
  #grad_motif_names <- sort(paste0(paste0("grad_", model_name,"_", c(1:motif_length), "_"), rep(c(1:4), motif_length)))
  #grad_motif_address <- paste0(main_directory, "Motif_opt_grad/", grad_motif_names, ".wtmx" )
  
  base_output_address <- paste0(main_directory, "Motif_opt_output/base_out_", model_name)
  #grad_output_names <- sort(paste0(paste0("grad_out_", model_name ,"_", c(1:motif_length), "_"), rep(c(1:4), motif_length)))
  #grad_output_address <- paste0(main_directory, "Motif_opt_output/", grad_output_names )
  
  
  for(cur_step in 1:nu_steps){
    # base model run
    
    job_file_name_1 <- paste0("grad_pwm_cons_", model_name, "_", cur_step, ".job")
    job_file_name_2 <- paste0("grad_GEM_run_", model_name, "_", cur_step, ".job")
    job_file_name_3 <- paste0("grad_calc_pwm_update_", model_name, "_", cur_step, ".job")
    job_file_name_4 <- paste0("GEMSTAT_train_", model_name, "_", cur_step, ".job")
    
    updated_motif_address_RData_current <- paste0(main_directory,
                                                  "Motif_opt_grad/base_motif_",
                                                  model_name, "_", cur_step, ".RData")
    cur_motif_opt_annotation_address <-  paste0(main_directory,
                                                "Annotation_steps/base_motif_",
                                                model_name, "_", cur_step, ".RData")
    if(cur_step == 1){
      updated_motif_address_RData_previous <- paste0(main_directory,
                                                     "base_motif.RData")
      updated_par_address_previous <- paste0(main_directory, trained_par_address)
      base_motif_address <- paste0(main_directory,  "motifs.wtmx")
    }else{
      updated_motif_address_RData_previous <- paste0(main_directory,
                                                     "Motif_opt_grad/base_motif_",
                                                     model_name, "_", (cur_step - 1), ".RData")
      updated_par_address_previous <- paste0(main_directory, "Motif_opt_output/", "updated_par_",
                                             model_name, "_", (cur_step - 1), ".txt")
      base_motif_address <- paste0(main_directory, "Motif_opt_grad/", "base_motif_", model_name, ".wtmx")
    }

    
    updated_par_address_current  <- paste0(main_directory, "Motif_opt_output/",
                                           "updated_par_",model_name, "_", cur_step, ".txt")
    
    train_seq_file_address_grad <- paste0(main_directory, "Motif_opt_input/seq_step_", cur_step, ".fa")
    train_seq_file_address_GS <- paste0(main_directory, "Motif_opt_input/seq_GEMS_step_", cur_step, ".fa")
    train_lab_file_address_grad <- paste0(main_directory, "Motif_opt_input/Label_step_", cur_step, ".tab")
    train_lab_file_address_GS <- paste0(main_directory, "Motif_opt_input/Label_GEMS_step_", cur_step, ".tab")
    train_wei_file_address_GS <- paste0(main_directory, "Motif_opt_input/Weight_GEMS_step_", cur_step, ".tab")
    
    train_seq_file_address_grad_annot <- paste0(main_directory, "Annotation_steps/seq_step_", cur_step, ".annot")
    train_seq_file_address_GS_annot <- paste0(main_directory, "Annotation_steps/seq_GEMS_step_", cur_step, ".annot")
    na_test <- 0
    na_train <- 2
    
    #first job (construct motifs)
    cat(c("Rscript", "--vanilla", "pwm_gradient_motif_writer_script.R",
          updated_motif_address_RData_previous,
          motif_index,
          grad_calc_step, 
          model_number_name, 
          paste0(cur_step, "\n")), 
        sep = " ", file = job_file_name_1, append = F)
    
    #second job (run GEMSTAT for gradients)
    cat(c(GEMSTAT_call,
          rep(paste0("-a ", train_seq_file_address_grad_annot), as.integer(annotation_fixed)),
          "-s", train_seq_file_address_grad,
          "-e", train_lab_file_address_grad,
          "-p", updated_par_address_previous, 
          "-m", base_motif_address, 
          "-na", na_test,
          "-fo", paste0(base_output_address, "\n")), 
        sep = " ", file = job_file_name_2, append = F)
    for(cur_pos in 1:motif_length){
      for(cur_base in 1:4){
        cat(c(GEMSTAT_call, 
              rep(paste0("-a ", train_seq_file_address_grad_annot), as.integer(annotation_fixed)),
              "-s", train_seq_file_address_grad,
              "-e", train_lab_file_address_grad,
              "-p", updated_par_address_previous, 
              "-m",  paste0(main_directory, "Motif_opt_grad/grad_",model_name, "_", cur_pos, "_", cur_base, ".wtmx" ), 
              "-na", na_test,
              "-fo", paste0(main_directory, "Motif_opt_output/grad_out_",model_name, "_", cur_pos, "_", cur_base, "\n")), 
            sep = " ", file = job_file_name_2, append = T)
      }
    }
    
    #Third job (gradient calculation, motif update)
    cat(c("Rscript", "--vanilla", "pwm_gradient_calculator_motif_update_script.R",
          updated_motif_address_RData_previous,
          motif_index,
          grad_calc_step, 
          grad_move_step,
          model_number_name, 
          cur_step,
          updated_par_address_previous,
          paste0(train_wei_file_address_GS, "\n")), 
        sep = " ", file = job_file_name_3, append = F)
    
    #Fourth job (GEMSTAT bound check -- GEMSTAT training)
    lb_address <- paste0(main_directory, "bounds/lower.par")
    ub_address <- paste0(main_directory, "bounds/upper.par")
    cat(c("Rscript", "--vanilla", "GEMSTAT_param_bound_checker_script.R", 
          updated_par_address_previous, lb_address, paste0(ub_address, "\n")),
        sep = " ", file = job_file_name_4, append = F)
    cat(c(GEMSTAT_call, 
          rep(paste0("-a ", train_seq_file_address_GS_annot), as.integer(annotation_fixed)),
          "-s", train_seq_file_address_GS,
          "-e", train_lab_file_address_GS,
          "-train_weights", train_wei_file_address_GS,
          "-p", updated_par_address_previous, 
          "-m", paste0(main_directory, "Motif_opt_grad/", "base_motif_", model_name, ".wtmx"), 
          "-na", na_train,
          "-po", updated_par_address_current,
          "-fo", paste0(base_output_address, "\n")), 
        sep = " ", file = job_file_name_4, append = T)
  }
}
################################################################################################################################################
# example
# aa_call <-  "/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1 -a Annotation/train.ann "
# 
# setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8")
# motif_optimizer_job_writer(nu_steps = 3,
#                            trained_par_address= "Trained_par/log_par_3118_1.txt.Filter",
#                            model_number_name = 3118,
#                            motif_length = 18,
#                            grad_move_step = 1,
#                            grad_calc_step = 0.01,
#                            GEMSTAT_call = aa_call,
#                            motif_index = 1,
#                            main_directory ="/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/")
# 
################################################################################################################################################

motif_optimizer_job_2_submit <- function(model_name_number,
                                         motif_length,
                                         nu_steps,
                                         output_file_name = paste0("job_to_submit_motif_optimizer_",model_name_number,".sh")){
  my_job_types <- c("grad_pwm_cons_", 
                    "grad_GEM_run_",
                    "grad_calc_pwm_update_", 
                    "GEMSTAT_train_")
  my_job_names <- paste0(my_job_types, model_name_number)
  nu_process <- c(1, (motif_length*4 + 1), 1, 1)
  sub_cnt <- 1
  submit_exe_file <- paste0("submit_to_Exec_", model_name_number, ".sh")
  cat("#!/bin/bash\n", file = output_file_name, sep = "", append = F)
  cat("#!/bin/bash\n", file = submit_exe_file, sep = "", append = F)
  for(cur_step in 1:nu_steps){
    for(cur_job in 1:length(my_job_names)){
      cur_job_pre <- paste0(my_job_names[cur_job], "_", cur_step)
      cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
            paste0(cur_job_pre, ".job"),
            nu_process[cur_job],
            paste0("tmp_", cur_job_pre),
            ">", 
            paste0("job_", model_name_number, "_", sub_cnt, ".submit", "\n")),
          sep = " ", file = output_file_name, append = T)
      cat(c("chmod +x", paste0("job_", model_name_number, "_", sub_cnt, ".submit", "\n")), 
          sep = " ", append = T, file = submit_exe_file)
      sub_cnt <- sub_cnt + 1
    }
  }
}
################################################################################################################################################
# example
# motif_optimizer_job_2_submit(model_name_number = 3118,
#                              motif_length = 18,
#                              nu_steps = 10,
#                              output_file_name = "job_to_submit_motif_optimizer.sh")
# 

################################################################################################################################################

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
                 "\n")), sep = " ", file = filename ,append = (cur_job_ind != 1))
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
# DAGmanConstructor(submit_prefix = "job_3118",
#                   start_ind=1,
#                   end_ind=40,
#                   my_order = c(1:40),
#                   filename="my_DagMan")
########################################################################################################################
########################################################################################################################

########################################################################################################################
########################################################################################################################

GEMSTAT_param_bound_checker <- function(par_file, ub_file, lb_file){
  library(rjson)
  # read par
  cur_par_file <- fromJSON(file = par_file)
  if(typeof(cur_par_file$qbtm) != "list"){
    cur_par_file$qbtm <- list(cur_par_file$qbtm)
  }
  # read ub
  cur_ub_file <- fromJSON(file = ub_file)
  if(typeof(cur_ub_file$qbtm) != "list"){
    cur_ub_file$qbtm <- list(cur_ub_file$qbtm)
  }
  # read lb
  cur_lb_file <- fromJSON(file = lb_file)
  if(typeof(cur_lb_file$qbtm) != "list"){
    cur_lb_file$qbtm <- list(cur_lb_file$qbtm)
  }
  
  updated_cur_par_file <- cur_par_file
  my_par <- cur_par_file
  
  nu_TFs <- length(my_par$tfs)
  TF_annot <- numeric(nu_TFs)
  TF_annot_ub <- numeric(nu_TFs)
  TF_annot_lb <- numeric(nu_TFs)
  TF_binding <- numeric(nu_TFs)
  TF_binding_ub <- numeric(nu_TFs)
  TF_binding_lb <- numeric(nu_TFs)
  TF_alpha <- numeric(nu_TFs)
  TF_alpha_ub <- numeric(nu_TFs)
  TF_alpha_lb <- numeric(nu_TFs)
  names(TF_annot) <- names(my_par$tfs)
  names(TF_annot_ub) <- names(my_par$tfs)
  names(TF_annot_lb) <- names(my_par$tfs)
  names(TF_binding) <- names(my_par$tfs)
  names(TF_binding_lb) <- names(my_par$tfs)
  names(TF_binding_ub) <- names(my_par$tfs)
  names(TF_alpha) <- names(my_par$tfs)
  names(TF_alpha_lb) <- names(my_par$tfs)
  names(TF_alpha_ub) <- names(my_par$tfs)
  for(cur_tf in 1:nu_TFs){
    TF_annot[cur_tf] <- my_par$tfs[[cur_tf]]$annot_thresh
    TF_annot_lb[cur_tf] <- cur_lb_file$tfs[[cur_tf]]$annot_thresh
    TF_annot_ub[cur_tf] <- cur_ub_file$tfs[[cur_tf]]$annot_thresh
    if(TF_annot[cur_tf] >= TF_annot_ub[cur_tf]){
      updated_cur_par_file$tfs[[cur_tf]]$annot_thresh <- cur_ub_file$tfs[[cur_tf]]$annot_thresh - 1e-3
    }
    if(TF_annot[cur_tf] <= TF_annot_lb[cur_tf]){
      updated_cur_par_file$tfs[[cur_tf]]$annot_thresh <- cur_lb_file$tfs[[cur_tf]]$annot_thresh + 1e-3
    }
    
    
    TF_binding[cur_tf] <- my_par$tfs[[cur_tf]]$maxbind
    TF_binding_lb[cur_tf] <- cur_lb_file$tfs[[cur_tf]]$maxbind
    TF_binding_ub[cur_tf] <- cur_ub_file$tfs[[cur_tf]]$maxbind
    if(TF_binding[cur_tf] >= TF_binding_ub[cur_tf]){
      updated_cur_par_file$tfs[[cur_tf]]$maxbind <- cur_ub_file$tfs[[cur_tf]]$maxbind - 1e-3
    }
    if(TF_binding[cur_tf] <= TF_binding_lb[cur_tf]){
      updated_cur_par_file$tfs[[cur_tf]]$maxbind <- cur_lb_file$tfs[[cur_tf]]$maxbind + 1e-3
    }
    TF_alpha[cur_tf] <- my_par$tfs[[cur_tf]]$alpha_a
    TF_alpha_lb[cur_tf] <- cur_lb_file$tfs[[cur_tf]]$alpha_a
    TF_alpha_ub[cur_tf] <- cur_ub_file$tfs[[cur_tf]]$alpha_a
    if(TF_alpha[cur_tf] >= TF_alpha_ub[cur_tf]){
      updated_cur_par_file$tfs[[cur_tf]]$alpha_a <- cur_ub_file$tfs[[cur_tf]]$alpha_a - 1e-3
    }
    if(TF_alpha[cur_tf] <= TF_alpha_lb[cur_tf]){
      updated_cur_par_file$tfs[[cur_tf]]$alpha_a <- cur_lb_file$tfs[[cur_tf]]$alpha_a + 1e-3
    }
  }
  
  # intercations
  Coop_par <- unlist(my_par$inter)
  Coop_par_lb <- unlist(cur_lb_file$inter)
  Coop_par_ub <- unlist(cur_ub_file$inter)
  if(length(Coop_par) > 0){
    for(c_coop in 1:length(Coop_par)){
      if(Coop_par[c_coop] >= Coop_par_ub[c_coop]){
        updated_cur_par_file$inter[[c_coop]] <- cur_ub_file$inter[[c_coop]] - 1e-3
      }
      if(Coop_par[c_coop] <= Coop_par_lb[c_coop]){
        updated_cur_par_file$inter[[c_coop]] <- cur_lb_file$inter[[c_coop]] + 1e-3
      }
    }
  }
  
  # qBTM
  qbtm_par <- unlist(my_par$qbtm)
  qbtm_par_lb <- unlist(cur_lb_file$qbtm)
  qbtm_par_ub <- unlist(cur_ub_file$qbtm)
  for(c_q in 1:length(qbtm_par)){
    if(qbtm_par[c_q] >= qbtm_par_ub[c_q]){
      updated_cur_par_file$qbtm[[c_q]] <- cur_ub_file$qbtm[[c_q]] - 1e-3
    }
    if(qbtm_par[c_q] <= qbtm_par_lb[c_q]){
      updated_cur_par_file$qbtm[[c_q]] <- cur_lb_file$qbtm[[c_q]] + 1e-3
    }
  }
  
  # beta
  beta_par <- my_par$enh[[1]]$beta
  beta_par_lb <- cur_lb_file$enh[[1]]$beta
  beta_par_ub <- cur_ub_file$enh[[1]]$beta
  if(beta_par >= beta_par_ub){
    updated_cur_par_file$enh[[1]]$beta <- cur_ub_file$enh[[1]]$beta - 1e-3
  }
  if(beta_par <= beta_par_lb){
    updated_cur_par_file$enh[[1]]$beta <- cur_ub_file$enh[[1]]$beta - 1e-3
  }
  
  # logistic
  logistic_par <- unlist(my_par$log_Reg[[1]])
  logistic_par_lb <- unlist(cur_lb_file$log_Reg[[1]])
  logistic_par_ub <- unlist(cur_ub_file$log_Reg[[1]])
  for(c_l in 1:length(logistic_par)){
    if(logistic_par[c_l] >= logistic_par_ub[c_l]){
      updated_cur_par_file$log_Reg[[1]][[c_l]] <- cur_ub_file$log_Reg[[1]][[c_l]] - 1e-3
    }
    if(logistic_par[c_l] <= logistic_par_lb[c_l]){
      updated_cur_par_file$log_Reg[[1]][[c_l]] <- cur_lb_file$log_Reg[[1]][[c_l]] + 1e-3
    }
  }
  cat(toJSON(updated_cur_par_file), append = F, 
      file = par_file)
}
########################################################################################################################
# example
# GEMSTAT_param_bound_checker(par_file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/log_par_3118_1.txt copy.Filter",
#                             ub_file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/bounds/upper.par",
#                             lb_file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/bounds/lower.par")
########################################################################################################################




