# GEMSTAT motif optimizer
# Optimizing ER motif in a gradient based way
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
################################################################################################################################################
###################################################                     ########################################################################
###################################################      Functions      #######################################################################
###################################################                     ########################################################################
################################################################################################################################################


read_output_train_test_GEMSTAT_indiv_onlyPRC <- function(output_file){
  
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
    cns_new <- max(0.001, min(1, my_motif[cur_pos, cur_cns] + grad_move_step*gradient_mat[cur_pos, cur_cns]))
    # print(paste0("cur_pos is ", cur_pos, " cur cons is ", cur_cns,
    #              " old value is ", my_motif[cur_pos, cur_cns]," new value is ", cns_new))
    # print(paste0("gradient value is: ", gradient_mat[cur_pos, cur_cns]) )
    # print(paste0("grad_move_step is  ", grad_move_step))
    my_motif_updated[cur_pos, cur_cns] <- cns_new
    for(cur_base in setdiff(c(1:4), cur_cns)){
      base_new <- ifelse((!constrained_rows[cur_pos]),
                         max(0.001, min(1, my_motif[cur_pos, cur_base] + grad_move_step*gradient_mat[cur_pos, cur_base])),
                         max(0.001, min(my_motif[cur_pos, cur_cns],
                                        my_motif[cur_pos, cur_base] + grad_move_step*gradient_mat[cur_pos, cur_base])))
      
      # print(paste0("cur_pos is ", cur_pos, " cur base is ", cur_base, 
      #              " old value was ", my_motif[cur_pos, cur_base], 
      #              " new value is ", base_new))
      # print(paste0("gradient value is: ", gradient_mat[cur_pos, cur_base]) )
      # print(paste0("grad_move_step is  ", grad_move_step))
      my_motif_updated[cur_pos, cur_base] <- base_new
    }
    my_motif_updated[cur_pos,] <- my_motif_updated[cur_pos,]/sum(my_motif_updated[cur_pos,])
  }
  # print("updated motif")
  # print(my_motif_updated)
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

# example
aa_call <-  "/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr -s Sequence/train_seq.fa -e Labels/Label_train.txt -f TF_exp.tab  -o DIRECT -c Coop/coop.par  -p Trained_par/log_par_3118_1.txt.Filter -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1 -train_weights weight_train.txt -a Annotation/train.ann "
aa_names <- c("ESR1_2","ESR1_3","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","PPARD","RARA", "RELA","RUNX1","SP1","YBX1")
names(TF.motifs.Expanded_pseudo)[names(TF.motifs.Expanded_pseudo) == "NKX3-1"] <- "NKX3_1"

c(rep(T, 6), rep(F, 3), rep(T, 6), rep(F, 3))

save(list = c("TF.motifs.Expanded_pseudo"),
     file = "~/Documents/Shayan/BioInf/EstrogenReceptor/TFmotifsExpanded_pseudo.RData")

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8")
sink("optimizer_out_9.txt")

aa <- GEMSTAT_motif_optimizer(initial_motif_list=TF.motifs.Expanded_pseudo[aa_names],
                              .motif_index=1, 
                              nu_steps=100,
                              ..GEMSTAT_call=aa_call,
                              .grad_calc_step=0.01,
                              .grad_move_step=100,
                              .constrained_rows=c(rep(T, 6), rep(F, 3), rep(T, 6), rep(F, 3)),
                              .inputs_dir="~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/")
sink()
aa_first_optimized2 <- aa

MotifWriter(motif.List = aa_first_optimized2, pseudo = 0.001,
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/first_opt2")

aai <- c(1143, 1384, 2917, 3118, 4076)
i <- aai[5]
aaa <- read_output_train_test_GEMSTAT_indiv(output_file = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/Motif_opt_log_par_",
                                                                 i, "_1/optimized_motif_2iter_train_1_0.out"), .plot = T)

plot(aaa$ROC_curve)
plot(aaa$PRC_curve)
aaa2 <- read_output_train_test_GEMSTAT_indiv(output_file = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/Motif_opt_log_par_",
                                                                  i, "_1/optimized_motif_2iter_test_1_0.out"), .plot = T)
plot(aaa2$ROC_curve)
plot(aaa2$PRC_curve)

aaa3 <- read_output_train_test_GEMSTAT_indiv(output_file = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/Motif_opt_log_par_",
                                                                  i, "_1/optimized_motif_2iter_valid_1_0.out"), .plot = T)
plot(aaa3$ROC_curve)
plot(aaa3$PRC_curve)


aaa_af <- read_output_train_test_GEMSTAT_indiv(output_file = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/Motif_opt_log_par_",
                                                                 i, "_1/optimized_motif_2iter_train_1_1.out"), .plot = T)

plot(aaa_af$ROC_curve)
plot(aaa_af$PRC_curve)
aaa2_af <- read_output_train_test_GEMSTAT_indiv(output_file = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/Motif_opt_log_par_",
                                                                  i, "_1/optimized_motif_2iter_test_1_1.out"), .plot = T)
plot(aaa2_af$ROC_curve)
plot(aaa2_af$PRC_curve)

aaa3_af <- read_output_train_test_GEMSTAT_indiv(output_file = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/Motif_opt_log_par_",
                                                                  i, "_1/optimized_motif_2iter_valid_1_1.out"), .plot = T)
plot(aaa3_af$ROC_curve)
plot(aaa3_af$PRC_curve)



E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$train_PRC["output_par_4076_1"]
E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$test_PRC["output_log_par_4076_1_test"]
E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_PRC["output_log_par_4076_1_valid"]
E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$train_ROC["output_par_4076_1"]
E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$test_ROC["output_log_par_4076_1_test"]
E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_ROC["output_log_par_4076_1_valid"]

plot(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$train_results[["output_par_3118_1"]]$ROC_curve)
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$train_results[["output_par_3118_1"]]$PRC_curve)

plot(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_results[["output_log_par_3118_1_valid"]]$ROC_curve)
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_results[["output_log_par_3118_1_valid"]]$PRC_curve)

plot(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$test_results[["output_log_par_3118_1_test"]]$ROC_curve)
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$test_results[["output_log_par_3118_1_test"]]$PRC_curve)

# plot for each model, two line plots each showing AUC under PR or ROC in each step

aa_Acc_PRC <- matrix(nrow = 15, ncol = 3)
aa_Acc_ROC <- matrix(nrow = 15, ncol = 3)
rownames(aa_Acc_PRC) <- sort(paste0(paste0("model_", 
                               c(1143, 1384, 2917, 3118, 4076)), rep(c("_train", "_test", "_valid"), 5)))
rownames(aa_Acc_ROC) <- rownames(aa_Acc_PRC)
colnames(aa_Acc_PRC) <- c("0_0", "1_0", "1_1")
colnames(aa_Acc_ROC) <- c("0_0", "1_0", "1_1")

aa_Acc_ROC[1, ] <- c(0.626, 0.615, 0.643)
aa_Acc_PRC[1, ] <- c(0.351,0.347, 0.362)
aa_Acc_ROC[2, ] <- c(0.709,0.708, 0.738)
aa_Acc_PRC[2, ] <- c(0.38, 0.441, 0.462)
aa_Acc_ROC[3, ] <- c(0.630,0.623, 0.626)
aa_Acc_PRC[3, ] <- c(0.37, 0.327, 0.33)

aa_Acc_ROC[4, ] <- c(0.606, 0.608, 0.650)
aa_Acc_PRC[4, ] <- c(0.349,0.356, 0.365)
aa_Acc_ROC[5, ] <- c(0.703,0.725, 0.753)
aa_Acc_PRC[5, ] <- c(0.416, 0.465, 0.487)
aa_Acc_ROC[6, ] <- c(0.616,0.631, 0.633)
aa_Acc_PRC[6, ] <- c(0.369, 0.348, 0.328)

aa_Acc_ROC[7, ] <- c(0.615, 0.628, 0.655)
aa_Acc_PRC[7, ] <- c(0.357,0.357, 0.350)
aa_Acc_ROC[8, ] <- c(0.710,0.709, 0.745)
aa_Acc_PRC[8, ] <- c(0.397, 0.453, 0.476)
aa_Acc_ROC[9, ] <- c(0.629,0.636, 0.622)
aa_Acc_PRC[9, ] <- c(0.370, 0.340, 0.329)

aa_Acc_ROC[10, ] <- c(0.625, 0.631, 0.635)
aa_Acc_PRC[10, ] <- c(0.358,0.372, 0.344)
aa_Acc_ROC[11, ] <- c(0.709,0.705, 0.743)
aa_Acc_PRC[11, ] <- c(0.389, 0.440, 0.459)
aa_Acc_ROC[12, ] <- c(0.636,0.647, 0.627)
aa_Acc_PRC[12, ] <- c(0.376, 0.364, 0.328)

aa_Acc_ROC[13, ] <- c(0.627, 0.625, 0.660)
aa_Acc_PRC[13, ] <- c(0.356,0.372, 0.369)
aa_Acc_ROC[14, ] <- c(0.709,0.707, 0.747)
aa_Acc_PRC[14, ] <- c(0.388, 0.468, 0.470)
aa_Acc_ROC[15, ] <- c(0.636,0.637, 0.635)
aa_Acc_PRC[15, ] <- c(0.370, 0.356, 0.329)



par(mfrow = c(5,2), mar = c(2,2,1,1))
for(i in 1:5){
  plot(aa_Acc_ROC[3*(i-1) + 1,], type = "l", ylim = range(aa_Acc_ROC))
  lines(aa_Acc_ROC[3*(i-1) + 2,], col = 2)
  lines(aa_Acc_ROC[3*(i-1) + 3,], col = 3)
  
  plot(aa_Acc_PRC[3*(i-1) + 1,], type = "l", ylim = range(aa_Acc_PRC))
  lines(aa_Acc_PRC[3*(i-1) + 2,], col = 2)
  lines(aa_Acc_PRC[3*(i-1) + 3,], col = 3)
}


aa_fldr <- c("Motif_opt_log_par_1143_1",
             "Motif_opt_log_par_1384_1",
             "Motif_opt_log_par_2917_1", 
             "Motif_opt_log_par_3118_1",
             "Motif_opt_log_par_4076_1")

aa_all_res <- list()
aa_cnt <- 1
aa_cnt_3 <- 0
aa_Acc_ROC_new <- matrix(nrow = 15, ncol = 10)
colnames(aa_Acc_ROC_new) <- character(length = 10)
rownames(aa_Acc_ROC_new) <- character(length = 15)
aa_Acc_PRC_new <- matrix(nrow = 15, ncol = 10)
colnames(aa_Acc_PRC_new) <- character(length = 10)
rownames(aa_Acc_PRC_new) <- character(length = 15)
for(i in 1:5){
  for(j in (i-1):(i)){
    aa_cnt2 <- 1
    aa_cnt_3 <- aa_cnt_3 + 1
    colnames(aa_Acc_ROC_new)[aa_cnt_3] <- paste0(i, "_", j)
      
    
    for(k in aa_fldr){
      for(l in c("test", "train", "valid")){
        aa_all_res[[aa_cnt]] <- read_output_train_test_GEMSTAT_indiv(output_file = paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/",k, "/",
                                                                                          "optimized_motif_2iter_",l,"_",i,"_", j,".out"), .plot = T)
        aa_Acc_ROC_new[aa_cnt2, aa_cnt_3] <- aa_all_res[[aa_cnt]]$ROC_curve$auc
        aa_Acc_PRC_new[aa_cnt2, aa_cnt_3] <- aa_all_res[[aa_cnt]]$PRC_curve$auc.integral
        if(i == 1 & j == 0){
          rownames(aa_Acc_ROC_new)[aa_cnt2] <- paste0(k, "_", l)
        }else{
          stopifnot(rownames(aa_Acc_ROC_new)[aa_cnt2] == paste0(k, "_", l))
        }
        aa_cnt <- aa_cnt + 1
        aa_cnt2 <- aa_cnt2 + 1
        
      }
    }
  }
}

# aa_Acc_ROC_new <- cbind(aa_Acc_ROC[, 1],aa_Acc_ROC_new)
# aa_Acc_PRC_new <- cbind(aa_Acc_PRC[, 1],aa_Acc_PRC_new)
# colnames(aa_Acc_ROC_new)[1] <- "0_0"
# colnames(aa_Acc_PRC_new)[1] <- "0_0"
par(mfrow = c(5,2), mar = c(3,3,1,1))
for(i in 1:5){
  plot(aa_Acc_ROC_new[3*(i-1) + 1,], type = "l", ylim = range(aa_Acc_ROC_new), lwd = 2, xaxt = "n")
  axis(side = 1, at = c(1:11), labels = colnames(aa_Acc_ROC_new), las =2)
  lines(aa_Acc_ROC_new[3*(i-1) + 2,], col = 2, lwd = 2)
  lines(aa_Acc_ROC_new[3*(i-1) + 3,], col = 3, lwd = 2)
  abline(h=seq(0.60, 0.80, 0.02), col = 4, lwd = 0.7, lty = 2)
  
  plot(aa_Acc_PRC_new[3*(i-1) + 1,], type = "l", ylim = range(aa_Acc_PRC_new), las = 2, lwd = 2, xaxt = "n")
  axis(side = 1, at = c(1:11), labels = colnames(aa_Acc_ROC_new), las =2)
  lines(aa_Acc_PRC_new[3*(i-1) + 2,], col = 2, lwd = 2)
  lines(aa_Acc_PRC_new[3*(i-1) + 3,], col = 3, lwd = 2)
  abline(h=seq(0.30, 0.55, 0.02), col = 4, lwd = 0.7, lty = 2)
}




################################################################################################################################################
################################################################################################################################################
# re-write the pipeline, such that it can be used in dagman setting on hal: write jobs for 
#  each step, write a dagman file to prioratize ruuning the jobs
# Given the number of steps and size of the batches create one set of sequence, label,
#  weights for each step, and write a job that calls this script
# Starting with step one, write a script that given a seed motif (probability), and gradient_step_size creates all
#  motifs for computing the gradient, then writes a job to run GEMSTAT given the above sequences and motifs
# job to compute gradient and create a new motif for the next step
# write jobs to do these for the given number of steps
# write dagman to prioratize the jobs
################################################################################################################################################
motif_optimizer_seq_label_weight_creator <- function(nu_steps,
                                                     batch_size, 
                                                     training_seq_pool,
                                                     training_seq_label){
  library(ShortRead)
  library(Biostrings)
  stopifnot(length(names(training_seq_pool)) == length(training_seq_pool), 
            length(training_seq_pool) == length(training_seq_label),
            batch_size < length(training_seq_label), 
            all(training_seq_label %in% c(0, 1)))
  
  if(!dir.exists("Motif_opt_input")){
    dir.create("Motif_opt_input")
  }
  if(!dir.exists("Motif_opt_output")){
    dir.create("Motif_opt_output")
  }
  total_seq_nu <- batch_size * nu_steps
  training_seq_nu <- length(training_seq_pool)
  nu_shuff <- ceiling(total_seq_nu/training_seq_nu)
  all_seq <- character(0)
  all_label <- integer(0)
  for(c_s in 1:nu_shuff){
    cur_samp <- sample(x = c(1:training_seq_nu), 
                       size = training_seq_nu, 
                       replace = F)
    all_seq <- c(all_seq, training_seq_pool[cur_samp])
    all_label <- c(all_label, training_seq_label[cur_samp])
  }
  ind_cnt <- 1
  seq_list <- list()
  label_list <- list()
  for(cur_step in 1:nu_steps){
    seq_list[[cur_step]] <- all_seq[(ind_cnt):(ind_cnt + batch_size - 1)]
    label_list[[cur_step]] <- all_label[(ind_cnt):(ind_cnt + batch_size - 1)]
    # write sequence
    writeFasta(DNAStringSet(seq_list[[cur_step]]), width = 1000,
               file = paste0("Motif_opt_input/seq_step_",cur_step,".fa"))
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
          file = paste0("Motif_opt_input/Label_step_", cur_step, ".tab", sep = "\t" ))
    }
    # write weights --> skip weights for now, see if there are errors
    ind_cnt <- ind_cnt + batch_size
  }
}
################################################################################################################################################
# write inputs beforehand
# create motifs for gradient --> a job to do this given the input motif and gradient step
# create jobs to run GEMSTAT with created motifs
# run jobs
################################################################################################################################################
pwm_gradient_motif_job_writer <- function(my_motif_list,
                                          motif_index,
                                          grad_calc_step, 
                                          constrained_rows,
                                          trained_par_address,
                                          #inputs_dir="./", 
                                          .GEMSTAT_call){
  # my_motif_list: is a list where each entry is a probability based PWM, where each row is a position and columns are: A, C, G, T: each row sums to one
  # motif_index : is the index of the motif to be modified
  # grad_calc_step: is the step size for the change used to compute gradient
  # constrained_rows: is a logical vector, True for rows that optimization should be constrained and F otherwise.
  #  constrain is to keep the consensus for that position unchanged.
  #  # consens: is a character vector containing the consensus sequence for the motif
  # inputs_dir: is the directory where GEMSTAT input files are located (it should end with /)
  # trained_par_address : is the address of the trained par to be used
  # .GEMSTAT_call: is the command to call GEMSTAT, without containing the -p -m -fo flags: these will be added afterwards
  stopifnot(ncol(my_motif_list[[motif_index]]) == 4)
  new_motif_list <- list()
  pos_cnt <- 0
  if(!dir.exists("Motif_opt_grad")){
    dir.create("Motif_opt_grad")
  }
  
  aan<- unlist(strsplit(trained_par_address, split = "\\/"))
  aan <- aan[length(aan)]
  model_name <- unlist(strsplit(aan, split = "\\."))[1]
  
  my_motif <- my_motif_list[[motif_index]]
  
  MotifWriter(motif.List = my_motif_list, pseudo = 0.001,
              output.File.Name = paste0("Motif_opt_grad/", "base_motif_", model_name))
  job_file_name <- paste0("gradient_GEM_run_", model_name, ".job")
  
  cat(paste0(GEMSTAT_call," -p ", trained_par_address," -m ", paste0("Motif_opt_grad/", "base_motif_", model_name, ".wtmx"),
             " -fo ", paste0("Motif_opt_output/base_out_", model_name), "\n"),
      file = job_file_name, append = F)

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
      new_mot_list <- my_motif_list
      new_mot_list[[motif_index]] <- new_motif_list[[pos_cnt]]
      MotifWriter(motif.List = new_mot_list, pseudo = 0.001,
                  output.File.Name = paste0("Motif_opt_grad/", "grad_", model_name, "_", cur_pos, "_", cur_base))
      cat(paste0(GEMSTAT_call," -p ", trained_par_address, " -m ", paste0("Motif_opt_grad/", "grad_", model_name, "_", cur_pos, "_", cur_base, ".wtmx"),
                 " -fo ", paste0("Motif_opt_output/out_", model_name, "_",cur_pos, "_",cur_base), "\n"),
          file = job_file_name, append = T)
    }
  }

}
################################################################################################################################################

# create a script that calls this function at every step, with the right arguments
################################################################################################################################################
pwm_gradient_calculator_motif_update <- function(my_motif_list,
                                                 motif_index,
                                                 grad_move_step, 
                                                 constrained_rows, 
                                                 trained_par_address){
  aan<- unlist(strsplit(trained_par_address, split = "\\/"))
  aan <- aan[length(aan)]
  model_name <- unlist(strsplit(aan, split = "\\."))[1]
  
  base_perf <- read_output_train_test_GEMSTAT_indiv_onlyPRC(output_file = paste0("Motif_opt_output/base_out_", model_name))
  
  out_files <- list.files(path = "Motif_opt_output/",
                          pattern = paste0("out_", model_name, "*"))
  out_row <- as.integer(unlist(lapply(strsplit(out_files, split = "_"), "[[", 3)))
  out_bas <- as.integer(unlist(lapply(strsplit(out_files, split = "_"), "[[", 4)))
  PRC_holder <- numeric(length = length(out_files))
  gradient_holder <- matrix(nrow = nrow(my_motif_list[[motif_index]]), ncol = 4)
  for(cur_out in 1:length(out_files)){
    gradient_holder[out_row[cur_out], out_bas[cur_out]] <- read_output_train_test_GEMSTAT_indiv_onlyPRC(output_file = paste0("Motif_opt_output/", out_files[cur_out])) - base_perf
  }
  updated_motif <- pwm_update(my_motif = my_motif_list[[motif_index]],
                              gradient_mat = gradient_holder,
                              grad_move_step = grad_move_step, 
                              constrained_rows = constrained_rows)
  working_mot_list <- my_motif_list
  working_mot_list[[motif_index]] <- updated_motif
  MotifWriter(motif.List = working_mot_list, pseudo = 0.001,
              output.File.Name = paste0("Motif_opt_grad/", "base_motif_", model_name))
}
# if(exists("myFirstFun", mode = "function"))
#   source("MyUtils.R")

## !/usr/bin/env Rscript
#  args = commandArgs(trailingOnly=TRUE)



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
  
  break_points <- c(break_points_binding, break_points_alpha, break_points_coop)
  
  filtered_par_mat <- all_pars[, all_points==1]
  filtered_break_poitns <- break_points[all_points==1]
  bit_matrix <- matrix(nrow = nrow(filtered_par_mat), ncol = ncol(filtered_par_mat))
  for(c_mod in 1:nrow(filtered_par_mat)){
    bit_matrix[c_mod, ] <- as.integer(filtered_par_mat[c_mod, ] < filtered_break_poitns)
  }
  return(bit_matrix)
}
################################################################################################################################################
################################################################################################################################################
# example
aa_bind_ff_ens <- c(0,0,1,0,0,0,1,0,1,0,0,0,0,1,1,0,1,0)
aa_alpha_ff_ens<- c(0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,0)

# choose three distinct models
aa_perf_srt <- names(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_PRC)[sort(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_PRC,
                                                                         decreasing = T, index.return=T)$ix]
aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
aa_perf_srt22 <- paste0("log_", aa_perf_srt2)
#E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$binding[aa_perf_srt22[1:100],]

aa_p_hash <- parameter_hash(binding_mat = E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$binding[aa_perf_srt22[1:100],],
                            alpha_mat = E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$alpha[aa_perf_srt22[1:100],],
                            coop_mat = E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$coop[aa_perf_srt22[1:100],],
                            bind_ens = aa_bind_ff_ens,
                            alpha_ens = aa_alpha_ff_ens, 
                            coop_ens = rep(0, 10), 
                            binding_range = rbind(rep(0.01, 18), rep(250.5, 18)),
                            alpha_range = rbind(rep(1e-6, 18), rep(100, 18)),
                            coop_range = rbind(rep(0.01, 10), rep(100.5, 10)),
                            binding_log = F, 
                            alpha_log = F, 
                            coop_log = F)

################################################################################################################################################
################################################################################################################################################
names(aa_newmot) <- names(TF.motifs.Expanded_pseudo_count)
aa_newmot[[1]]


MotifWriter(motif.List = aa_newmot[aa_names], pseudo = 0.001,
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/base_motif_prob")
MotifWriter(motif.List = aa_newmot[aa_names], pseudo = 0.001,
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/base_motif_prob")

# pick top 50 models from experimennt 8
# E_RNA_GEMSTAT_Ensemble_Parlist[[8]]
# E_RNA_GEMSTAT_Ensemble_Outlist[[8]]

# write a fucntion that takes a motif list, index of the TFs to be modified, 
#  constrains for each position, a gradient direction
#  and creates variants of it by changing every cell 

# coop 
# "ESR1_2", "ESR1_3" 10
# "ESR1_2", "FOXA1" 45
# "ESR1_2", "GATA3" 35
# "ESR1_2", "PBX1" 70
# "ESR1_2", "PGR" 50
# "ESR1_2", "PPARD" 20
# "ESR1_2", "RARA" 20
# "ESR1_2", "YBX1" 100
# "JUN_1", "JUN_1" 15
# "RUNX1", "YBX1"  20

aa_perf_srt <-names(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_PRC)[sort(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_PRC,
                                                                        decreasing = T, index.return=T)$ix]
aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
aa_perf_srt22 <- paste0("log_", aa_perf_srt2)

# ## using DAMO
# # write ESR1_2 motif in their style
# ncol(TF.motifs.Expanded_pseudo_t[[4]])
# 
# seqLogo::seqLogo(TF.motifs.Expanded_pseudo_t[[4]][,1:15])
# aa_mymotif <- TF.motifs.Expanded_pseudo_t[[4]][,1:15]
# aa_a <- c("A", "C", "G", "T")
# cat(c("> ESR1_2\n"), file = "ESR1_2_Damo.pwm", append = "F")
# for(i in 1:nrow(aa_mymotif)){
#   cat(c(aa_a[i], "|" , c(aa_mymotif[i,], "\n")), sep = " ",
#       file = "ESR1_2_Damo.pwm", append = "T" )
# }
# 
# aalab <- unlist(lapply(strsplit(data_partition_mut_GEMSTAT_trimmed[[1]], split="_"), "[[", 1))
# aa_train_ind_pos <- data_partition_mut_GEMSTAT_trimmed[[1]][sample(x = which(aalab %in% "pos"), size = 270, replace = F)]
# aa_train_ind_neg <- data_partition_mut_GEMSTAT_trimmed[[1]][sample(x = which(aalab %in% "neg"), size = 1002, replace = F)]
# 
# aamatchp <- match(aa_train_ind_pos, names(Positive_set_seq_list_char_1000bp[[4]]))
# aamatchn <- match(aa_train_ind_neg, names(Negative_set_seq_list_char_1000bp[[4]]))
# 
# aatpos <- Positive_set_seq_list_char_1000bp[[4]][aamatchp]
# aatneg <- Negative_set_seq_list_char_1000bp[[4]][aamatchn]
# aatneg2 <- sample(aatneg, size = length(aatpos), replace = F)
# 
# 
# writeFasta(DNAStringSet(aatpos), width = 1000,
#            file = "~/Documents/DAMO/DAMO/Test/pos_seq.fa")
# writeFasta(DNAStringSet(aatneg2), width = 1000,
#            file = "~/Documents/DAMO/DAMO/Test/neg_seq_1.fa")

################################################################################################################################################
# do a small experiment ---> take three models (good performing-distinct) --->
#  perform two steps of motif optimization on each ---> re-train thermodynamic params --> 
#  look at performance --> re-optimize motifs --> retrain GEMSTAT --> repeat 2 times.





aa_bind_ff_ens <- c(0,0,1,0,0,0,1,0,1,0,0,0,0,1,1,0,1,0)
aa_alpha_ff_ens<- c(0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,0)

# choose three distinct models
aa_perf_srt <- names(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_PRC)[sort(E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_PRC,
                                                                         decreasing = T, index.return=T)$ix]
aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
aa_perf_srt22 <- paste0("log_", aa_perf_srt2)
#E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$binding[aa_perf_srt22[1:100],]

aa_p_hash <- parameter_hash(binding_mat = E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$binding[aa_perf_srt22[1:500],],
                            alpha_mat = E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$alpha[aa_perf_srt22[1:500],],
                            coop_mat = E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$coop[aa_perf_srt22[1:500],],
                            bind_ens = aa_bind_ff_ens,
                            alpha_ens = aa_alpha_ff_ens, 
                            coop_ens = rep(0, 10), 
                            binding_range = rbind(rep(0.01, 18), rep(250.5, 18)),
                            alpha_range = rbind(rep(1e-6, 18), rep(100, 18)),
                            coop_range = rbind(rep(0.01, 10), rep(100.5, 10)),
                            binding_log = F, 
                            alpha_log = F, 
                            coop_log = F)

aa_dup_hash <- duplicated(aa_p_hash)
aa_chosen_model <- character(0)
aacnt <- 1
aa_cur <- 1
while(aacnt <= 25){
  if(!aa_dup_hash[aa_cur]){
    aa_chosen_model <- c(aa_chosen_model, rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[8]]$binding[aa_perf_srt22[1:500],])[aa_cur])
    aacnt <- aacnt + 1
  }
  aa_cur <- aa_cur + 1
  
}
aa_chosen_model
motif_optim_chosen_models_1 <- aa_chosen_model
# now for each of these models, do two rounds of motif optimization and then a round of GEMSTAT opt

############################################################################################################
# creating output for hal optimization
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8")

aa_names <- c("ESR1_2","ESR1_3","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","PPARD","RARA", "RELA","RUNX1","SP1","YBX1")
names(TF.motifs.Expanded_pseudo)[names(TF.motifs.Expanded_pseudo) == "NKX3-1"] <- "NKX3_1"

original_motif_list <- TF.motifs.Expanded_pseudo[aa_names]
save(list = c("original_motif_list"), file = "base_motif.RData")

aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[4]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
motif_optimizer_seq_label_weight_creator(nu_steps = 100,
                                         batch_size_motif = 200,
                                         batch_size_GEMSTAT = 400,
                                         training_seq_pool = GEMSTAT_Ensemble_train_SeqList[[4]],
                                         training_seq_label = aalab1)



aa_call_noannot <-  "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1"
aa_call_annot <-  "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1 -a Annotation/train.ann "

# -a Annotation/train.ann 


motif_optim_chosen_models_1
aa_sp_name <-  unlist(lapply(strsplit(motif_optim_chosen_models_1, split = "_"), "[[", 3))

for(aa_model in 1:length(motif_optim_chosen_models_1)){
  aa_adres <- paste0("Trained_par/", motif_optim_chosen_models_1[aa_model], ".txt.Filter")
  aa_nu <- aa_sp_name[aa_model]
  aa_dag <- paste0("job_", aa_nu)
  aa_dagf <- paste0("my_DagMan_", aa_nu)
  
  motif_optimizer_job_writer(nu_steps = 100,
                             trained_par_address= aa_adres,
                             model_number_name = aa_nu,
                             motif_length = 18,
                             grad_move_step = 1,
                             grad_calc_step = 0.01,
                             GEMSTAT_call = aa_call_noannot,
                             motif_index = 1,
                             main_directory ="/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/")
  motif_optimizer_job_2_submit(model_name_number = aa_nu,
                               motif_length = 18,
                               nu_steps = 100)
  DAGmanConstructor(submit_prefix = aa_dag,
                    start_ind=1,
                    end_ind=400,
                    my_order = c(1:400),
                    filename=aa_dagf)
}

# create script to make these exec and run dags
cat(c("#!/bin/bash\n"), file = "job_2_sub_executable_dagrun.sh", append = F)
for(aa_model in 1:length(motif_optim_chosen_models_1)){
  aa_nu <- aa_sp_name[aa_model]
  cat(c("chmod", "+x", paste0("job_to_submit_motif_optimizer_",aa_nu, ".sh\n" )),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("./job_to_submit_motif_optimizer_",aa_nu, ".sh\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c("chmod", "+x", paste0("submit_to_Exec_",aa_nu, ".sh\n" )),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("./submit_to_Exec_",aa_nu, ".sh\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("condor_submit_dag my_DagMan_",aa_nu, "\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
}

############################################################################################################
# write jobs to evaluate performance of the models in each step on training/test/validation sets

aa_sp_name <-  unlist(lapply(strsplit(motif_optim_chosen_models_1, split = "_"), "[[", 3))

# write motifs from RData

for(aa_model in 1:length(motif_optim_chosen_models_1)){
  aa_nu <- aa_sp_name[aa_model]
  for(aa_step in 1:100){
    aa_read_adres  <- paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Motif_opt_grad/base_motif_", 
                             aa_nu, "_", aa_step, ".RData")
    aa_write_adres <- paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Optimized_motifs/motif_", 
                            aa_nu, "_", aa_step)
    # job to write motifs
    cat(c("Rscript --vanilla motif_writer_script.R", aa_read_adres, paste0(aa_write_adres, "\n")), 
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)), file = "write_motifs.job")
    # job to run GEMSTAT on training set
    aa_GEM_call <- "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1"
    aa_train_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Motif_opt_input/seq_GEMS_step_", aa_step, ".fa")
    aa_train_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Motif_opt_input/Label_GEMS_step_", aa_step, ".tab")
    
    aa_test_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Sequence/test_seq.fa")
    aa_test_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Labels/Label_test.txt")
    
    aa_valid_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Sequence/valid_seq.fa")
    aa_valid_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Labels/Label_validation.txt")
    
    
    aa_cur_params <- paste0("-p /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Motif_opt_output/updated_par_",aa_nu,"_",aa_step,".txt")
    aa_cur_motifs <- paste0("-m ", aa_write_adres, ".wtmx")
    aa_na <- "-na 0"
    aa_cur_out_train <- paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Model_eval_tr_ts_val/train_out_", aa_nu, "_", aa_step, ".txt")
    aa_cur_out_test <-  paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Model_eval_tr_ts_val/test_out_", aa_nu, "_", aa_step, ".txt")
    aa_cur_out_valid <- paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/Model_eval_tr_ts_val/valid_out_", aa_nu, "_", aa_step, ".txt")
    
    # job to run GEMSTAT on test set
    # job to run GEMSTAT on valid set
    cat(c(aa_GEM_call, aa_train_seq, aa_train_lab, aa_cur_params,
          aa_cur_motifs, aa_na, paste0(aa_cur_out_train, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_train_eval.job")
    cat(c(aa_GEM_call, aa_test_seq, aa_test_lab, aa_cur_params,
          aa_cur_motifs, aa_na,  paste0(aa_cur_out_test, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_test_eval.job")
    cat(c(aa_GEM_call, aa_valid_seq, aa_valid_lab, aa_cur_params,
          aa_cur_motifs, aa_na,  paste0(aa_cur_out_valid, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_valid_eval.job")
  }
}
# write job2submit scripts

cat("#!/bin/bash\n", file = "write_motifs_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "write_motifs.job",
      "2500",
      paste0("tmp_", "write_motifs"),
      ">", 
      paste0("write_motifs.submit", "\n")),
    sep = " ", file = "write_motifs_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_train_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_train_eval.job",
      "2500",
      paste0("tmp_", "GEM_train_eval"),
      ">", 
      paste0("GEM_train_eval.submit", "\n")),
    sep = " ", file = "GEM_train_eval_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_test_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_test_eval.job",
      "2500",
      paste0("tmp_", "GEM_test_eval"),
      ">", 
      paste0("GEM_test_eval.submit", "\n")),
    sep = " ", file = "GEM_test_eval_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_valid_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_valid_eval.job",
      "2500",
      paste0("tmp_", "GEM_valid_eval"),
      ">", 
      paste0("GEM_valid_eval.submit", "\n")),
    sep = " ", file = "GEM_valid_eval_job2submit.sh", append = T)
############################################################################################################
# read performance in terms of AUROC and AUPRC for all models

aa_sp_name <-  unlist(lapply(strsplit(motif_optim_chosen_models_1, split = "_"), "[[", 3))
motif_optim_perf_auroc_train <- matrix(nrow = length(aa_sp_name), ncol = 100)
rownames(motif_optim_perf_auroc_train) <- aa_sp_name

motif_optim_perf_auroc_test <- matrix(nrow = length(aa_sp_name), ncol = 100)
rownames(motif_optim_perf_auroc_test) <- aa_sp_name

motif_optim_perf_auroc_valid <- matrix(nrow = length(aa_sp_name), ncol = 100)
rownames(motif_optim_perf_auroc_valid) <- aa_sp_name

motif_optim_perf_auprc_train <- matrix(nrow = length(aa_sp_name), ncol = 100)
rownames(motif_optim_perf_auprc_train) <- aa_sp_name

motif_optim_perf_auprc_test <- matrix(nrow = length(aa_sp_name), ncol = 100)
rownames(motif_optim_perf_auprc_test) <- aa_sp_name

motif_optim_perf_auprc_valid <- matrix(nrow = length(aa_sp_name), ncol = 100)
rownames(motif_optim_perf_auprc_valid) <- aa_sp_name

aa_Train_perf <- list()
aa_Test_perf <- list()
aa_valid_perf <- list()
aa_cnt <- 1
for(aa_model in 1:length(motif_optim_chosen_models_1)){
  print(aa_model)
  aa_nu <- aa_sp_name[aa_model]
  for(aa_step in 1:100){
    
    aa_tr_file <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/Model_eval_tr_ts_val/train_out_", aa_nu, "_", aa_step, ".txt")
    aa_ts_file <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/Model_eval_tr_ts_val/test_out_", aa_nu, "_", aa_step, ".txt")
    aa_va_file <-  paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/Model_eval_tr_ts_val/valid_out_", aa_nu, "_", aa_step, ".txt")
    if (file.exists(aa_tr_file)){
      aa_Train_perf[[aa_cnt]] <- read_output_train_test_GEMSTAT_indiv(output_file = aa_tr_file)
      motif_optim_perf_auroc_train[aa_model, aa_step] <- aa_Train_perf[[aa_cnt]]$ROC_curve[[2]]
      motif_optim_perf_auprc_train[aa_model, aa_step] <- aa_Train_perf[[aa_cnt]]$PRC_curve[[2]]
      
    }
    if (file.exists(aa_ts_file)){
      aa_Test_perf[[aa_cnt]] <- read_output_train_test_GEMSTAT_indiv(output_file = aa_ts_file)
      motif_optim_perf_auroc_test[aa_model, aa_step] <- aa_Test_perf[[aa_cnt]]$ROC_curve[[2]]
      motif_optim_perf_auprc_test[aa_model, aa_step] <- aa_Test_perf[[aa_cnt]]$PRC_curve[[2]]
      
    }
    if (file.exists(aa_va_file)){
      aa_valid_perf[[aa_cnt]] <- read_output_train_test_GEMSTAT_indiv(output_file =aa_va_file)
      motif_optim_perf_auroc_valid[aa_model, aa_step] <- aa_valid_perf[[aa_cnt]]$ROC_curve[[2]]
      motif_optim_perf_auprc_valid[aa_model, aa_step] <- aa_valid_perf[[aa_cnt]]$PRC_curve[[2]]
      
    }
    
    aa_cnt <- aa_cnt + 1
    
    
  }
}
############################################################################################################
# plot the results
boxplot.matrix(motif_optim_perf_auprc_train)
boxplot.matrix(motif_optim_perf_auprc_test)
boxplot.matrix(motif_optim_perf_auprc_valid)

aa_range_prc <- range(cbind(motif_optim_perf_auprc_train, motif_optim_perf_auprc_test, motif_optim_perf_auprc_valid), na.rm = T)
aa_range_roc <- range(cbind(motif_optim_perf_auroc_train, motif_optim_perf_auroc_test, motif_optim_perf_auroc_valid), na.rm = T)

aa_sp_name <-  unlist(lapply(strsplit(motif_optim_chosen_models_1, split = "_"), "[[", 3))


motif_optim_perf_auprc_compare <- matrix(0L,nrow = length(aa_sp_name), 
                                         ncol = 100 )
motif_optim_perf_auroc_compare <- matrix(0L,nrow = length(aa_sp_name), 
                                         ncol = 100 )

par(mfrow = c(5, 5), mar = c(2,2,2,2))
for(aa_model in 1:length(motif_optim_chosen_models_1)){
  plot(motif_optim_perf_auprc_train[aa_model,], type = "l", col = 1, ylim = aa_range_prc, main = aa_sp_name[aa_model])
  abline(h = E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$train_PRC[paste0("output_par_",aa_sp_name[aa_model], "_1")], col = 1, lty = 2)
  
  
  lines(motif_optim_perf_auprc_test[aa_model,], col = 2)
  aa_base_test <- E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$test_PRC[paste0("output_",motif_optim_chosen_models_1[aa_model], "_test")]
  abline(h = aa_base_test, col = 2, lty = 2)
  motif_optim_perf_auprc_compare[aa_model, ] <- motif_optim_perf_auprc_compare[aa_model, ] + as.integer(motif_optim_perf_auprc_test[aa_model,] > aa_base_test)
  
  
  aa_base_valid <- E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_PRC[paste0("output_",motif_optim_chosen_models_1[aa_model], "_valid")]
  lines(motif_optim_perf_auprc_valid[aa_model,], col = 3)
  abline(h = aa_base_valid, col = 3, lty = 2)
  motif_optim_perf_auprc_compare[aa_model, ] <- motif_optim_perf_auprc_compare[aa_model, ] + as.integer(motif_optim_perf_auprc_valid[aa_model,] > aa_base_valid)
  
}

par(mfrow = c(5, 5), mar = c(2,2,2,2))
for(aa_model in 1:length(motif_optim_chosen_models_1)){
  plot(motif_optim_perf_auroc_train[aa_model,], type = "l", col = 1, ylim = aa_range_roc, main = aa_sp_name[aa_model])
  abline(h = E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$train_ROC[paste0("output_par_",aa_sp_name[aa_model], "_1")], col = 1, lty = 2)
  
  lines(motif_optim_perf_auroc_test[aa_model,], col = 2)
  aa_base_test <- E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$test_ROC[paste0("output_",motif_optim_chosen_models_1[aa_model], "_test")]
  abline(h = aa_base_test, col = 2, lty = 2)
  motif_optim_perf_auroc_compare[aa_model, ] <- motif_optim_perf_auroc_compare[aa_model, ] + as.integer(motif_optim_perf_auroc_test[aa_model,] > aa_base_test)
  
  
  lines(motif_optim_perf_auroc_valid[aa_model,], col = 3)
  aa_base_valid <- E_RNA_GEMSTAT_Ensemble_Outlist[[8]]$valid_ROC[paste0("output_",motif_optim_chosen_models_1[aa_model], "_valid")]
  abline(h = aa_base_valid, col = 3, lty = 2)
  motif_optim_perf_auroc_compare[aa_model, ] <- motif_optim_perf_auroc_compare[aa_model, ] + as.integer(motif_optim_perf_auroc_valid[aa_model,] > aa_base_valid)
  
}


aaa <- motif_optim_perf_auprc_compare
which(rowSums(aaa == 2, na.rm = T) > 0)
which(motif_optim_perf_auprc_compare[8, ] == 2)

############################################################################################################
# write Jpegs of motif logos

aa_cur_model <- 457
for(i in 1:100){
  aapngfname <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/Optimized_motifs_100_logo/logo_",aa_cur_model,"_", i, ".png")
  load(paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_8/Motif_opt_grad/base_motif_", aa_cur_model,"_", i, ".RData" ))
  png(filename = aapngfname, width = 600, height = 300, units = "px" )
  seqLogo::seqLogo(pwm = t(original_motif_list[[1]]))
  dev.off()
}

############################################################################################################
# found a bug, rerunning dagmans

aa_sp_name <-  unlist(lapply(strsplit(motif_optim_chosen_models_1, split = "_"), "[[", 3))

# create script to make these exec and run dags
cat(c("#!/bin/bash\n"), file = "job_2_sub_executable_dagrun_2.sh", append = F)
for(aa_model in 1:length(motif_optim_chosen_models_1)){
  aa_nu <- aa_sp_name[aa_model]
  cat(c(paste0("condor_submit_dag -force my_DagMan_",aa_nu, "\n")),
      sep = " ", file = "job_2_sub_executable_dagrun_2.sh", append = T)
}



load("~/Desktop/base_motif_1384_5.RData")
seqLogo::seqLogo(pwm = t(original_motif_list[[2]]))
############################################################################################################
# need to change the objective function from AUPRC to logistic regression, there is a reason that people usually don't use it.
# changed
# running
# creating output for hal optimization
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12")

aa_names <- c("ESR1_2","FOXA1","GATA3","JUN_1","LEF1","NKX3_1","NR3C1",
              "NR5A2", "PAX2","PBX1","PGR","RARA", "RELA","RUNX1","SP1","TFAP2C","YBX1")
names(TF.motifs.Expanded_pseudo)[names(TF.motifs.Expanded_pseudo) == "NKX3-1"] <- "NKX3_1"

original_motif_list <- TF.motifs.Expanded_pseudo[aa_names]
save(list = c("original_motif_list"), file = "base_motif.RData")

aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[4]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
motif_optimizer_seq_label_weight_creator(nu_steps = 20,
                                         batch_size_motif = 200,
                                         batch_size_GEMSTAT = 400,
                                         training_seq_pool = GEMSTAT_Ensemble_train_SeqList[[4]],
                                         training_seq_label = aalab1)



aa_call_noannot <-  "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1"
aa_call_annot <-  "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1 -a Annotation/train.ann "

# -a Annotation/train.ann 


aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas1[1:25]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas2[1:25]]
aaadd3 <- union(aaadd1, aaadd2)

aa_sp_name <-  unlist(lapply(strsplit(aaadd3, split = "_"), "[[", 2))

for(aa_model in 1:length(aaadd3)){
  aa_adres <- paste0("Trained_par/", aaadd3[aa_model], ".txt.Filter")
  aa_nu <- aa_sp_name[aa_model]
  aa_dag <- paste0("job_", aa_nu)
  aa_dagf <- paste0("my_DagMan_", aa_nu)
  
  motif_optimizer_job_writer(nu_steps = 20,
                             trained_par_address= aa_adres,
                             model_number_name = aa_nu,
                             motif_length = nrow(TF.motifs.Expanded_pseudo$ESR1_2),
                             grad_move_step = 1,
                             grad_calc_step = 0.01,
                             GEMSTAT_call = aa_call_noannot,
                             motif_index = 1,
                             main_directory ="/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/")
  motif_optimizer_job_2_submit(model_name_number = aa_nu,
                               motif_length = nrow(TF.motifs.Expanded_pseudo$ESR1_2),
                               nu_steps = 20)
  DAGmanConstructor(submit_prefix = aa_dag,
                    start_ind=1,
                    end_ind=80,
                    my_order = c(1:80),
                    filename=aa_dagf)
}

# create script to make these exec and run dags
cat(c("#!/bin/bash\n"), file = "job_2_sub_executable_dagrun.sh", append = F)
for(aa_model in 1:length(aaadd3)){
  aa_nu <- aa_sp_name[aa_model]
  cat(c("chmod", "+x", paste0("job_to_submit_motif_optimizer_",aa_nu, ".sh\n" )),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("./job_to_submit_motif_optimizer_",aa_nu, ".sh\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c("chmod", "+x", paste0("submit_to_Exec_",aa_nu, ".sh\n" )),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("./submit_to_Exec_",aa_nu, ".sh\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("condor_submit_dag my_DagMan_",aa_nu, "\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
}
# rerunning because of a bug
cat(c("#!/bin/bash\n"), file = "job_2_sub_executable_dagrun_rerun.sh", append = F)
for(aa_model in 1:length(aaadd3)){
  aa_nu <- aa_sp_name[aa_model]
  # cat(c("chmod", "+x", paste0("job_to_submit_motif_optimizer_",aa_nu, ".sh\n" )),
  #     sep = " ", file = "job_2_sub_executable_dagrun_rerun.sh", append = T)
  # cat(c(paste0("./job_to_submit_motif_optimizer_",aa_nu, ".sh\n")),
  #     sep = " ", file = "job_2_sub_executable_dagrun_rerun.sh", append = T)
  # cat(c("chmod", "+x", paste0("submit_to_Exec_",aa_nu, ".sh\n" )),
  #     sep = " ", file = "job_2_sub_executable_dagrun_rerun.sh", append = T)
  # cat(c(paste0("./submit_to_Exec_",aa_nu, ".sh\n")),
  #     sep = " ", file = "job_2_sub_executable_dagrun_rerun.sh", append = T)
  cat(c(paste0("condor_submit_dag -force my_DagMan_",aa_nu, "\n")),
      sep = " ", file = "job_2_sub_executable_dagrun_rerun.sh", append = T)
}

# write jobs to evaluate performance of the models in each step on training/test/validation sets

aa_sp_name <-  unlist(lapply(strsplit(aaadd3, split = "_"), "[[", 2))

# write motifs from RData

for(aa_model in 1:length(aaadd3)){
  aa_nu <- aa_sp_name[aa_model]
  for(aa_step in 1:20){
    aa_read_adres  <- paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_grad/base_motif_", 
                             aa_nu, "_", aa_step, ".RData")
    aa_write_adres <- paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Optimized_motifs/motif_", 
                             aa_nu, "_", aa_step)
    # job to write motifs
    cat(c("Rscript --vanilla motif_writer_script.R", aa_read_adres, paste0(aa_write_adres, "\n")), 
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)), file = "write_motifs.job")
    # job to run GEMSTAT on training set
    aa_GEM_call <- "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1"
    aa_train_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_input/seq_GEMS_step_", aa_step, ".fa")
    aa_train_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_input/Label_GEMS_step_", aa_step, ".tab")
    
    aa_test_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Sequence/test_seq.fa")
    aa_test_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Labels/Label_test.txt")
    
    aa_valid_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Sequence/valid_seq.fa")
    aa_valid_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Labels/Label_validation.txt")
    
    
    aa_cur_params <- paste0("-p /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_output/updated_par_",aa_nu,"_",aa_step,".txt")
    aa_cur_motifs <- paste0("-m ", aa_write_adres, ".wtmx")
    aa_na <- "-na 0"
    aa_cur_out_train <- paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/train_out_", aa_nu, "_", aa_step, ".txt")
    aa_cur_out_test <-  paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/test_out_", aa_nu, "_", aa_step, ".txt")
    aa_cur_out_valid <- paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/valid_out_", aa_nu, "_", aa_step, ".txt")
    
    # job to run GEMSTAT on test set
    # job to run GEMSTAT on valid set
    cat(c(aa_GEM_call, aa_train_seq, aa_train_lab, aa_cur_params,
          aa_cur_motifs, aa_na, paste0(aa_cur_out_train, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_train_eval.job")
    cat(c(aa_GEM_call, aa_test_seq, aa_test_lab, aa_cur_params,
          aa_cur_motifs, aa_na,  paste0(aa_cur_out_test, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_test_eval.job")
    cat(c(aa_GEM_call, aa_valid_seq, aa_valid_lab, aa_cur_params,
          aa_cur_motifs, aa_na,  paste0(aa_cur_out_valid, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_valid_eval.job")
  }
}
# write job2submit scripts

cat("#!/bin/bash\n", file = "write_motifs_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "write_motifs.job",
      "840",
      paste0("tmp_", "write_motifs"),
      ">", 
      paste0("write_motifs.submit", "\n")),
    sep = " ", file = "write_motifs_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_train_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_train_eval.job",
      "840",
      paste0("tmp_", "GEM_train_eval"),
      ">", 
      paste0("GEM_train_eval.submit", "\n")),
    sep = " ", file = "GEM_train_eval_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_test_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_test_eval.job",
      "840",
      paste0("tmp_", "GEM_test_eval"),
      ">", 
      paste0("GEM_test_eval.submit", "\n")),
    sep = " ", file = "GEM_test_eval_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_valid_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_valid_eval.job",
      "840",
      paste0("tmp_", "GEM_valid_eval"),
      ">", 
      paste0("GEM_valid_eval.submit", "\n")),
    sep = " ", file = "GEM_valid_eval_job2submit.sh", append = T)
############################################################################################################
# read performance in terms of AUROC and AUPRC for all models

aa_sp_name <-  unlist(lapply(strsplit(aaadd3, split = "_"), "[[", 2))
motif_optim_perf_auroc_train_12 <- matrix(nrow = length(aa_sp_name), ncol = 50)
rownames(motif_optim_perf_auroc_train_12) <- aa_sp_name

motif_optim_perf_auroc_test_12 <- matrix(nrow = length(aa_sp_name), ncol = 50)
rownames(motif_optim_perf_auroc_test_12) <- aa_sp_name

motif_optim_perf_auroc_valid_12 <- matrix(nrow = length(aa_sp_name), ncol = 50)
rownames(motif_optim_perf_auroc_valid_12) <- aa_sp_name

motif_optim_perf_auprc_train_12 <- matrix(nrow = length(aa_sp_name), ncol = 50)
rownames(motif_optim_perf_auprc_train_12) <- aa_sp_name

motif_optim_perf_auprc_train_12_logistloss <- matrix(nrow = length(aa_sp_name), ncol = 50)
rownames(motif_optim_perf_auprc_train_12_logistloss) <- aa_sp_name

motif_optim_perf_auprc_test_12 <- matrix(nrow = length(aa_sp_name), ncol = 50)
rownames(motif_optim_perf_auprc_test_12) <- aa_sp_name

motif_optim_perf_auprc_test_12_logistloss <- matrix(nrow = length(aa_sp_name), ncol = 50)
rownames(motif_optim_perf_auprc_test_12_logistloss) <- aa_sp_name

motif_optim_perf_auprc_valid_12 <- matrix(nrow = length(aa_sp_name), ncol = 50)
rownames(motif_optim_perf_auprc_valid_12) <- aa_sp_name

motif_optim_perf_auprc_valid_12_logistloss <- matrix(nrow = length(aa_sp_name), ncol = 50)
rownames(motif_optim_perf_auprc_valid_12_logistloss) <- aa_sp_name

aa_Train_perf_12 <- list()
aa_Test_perf_12 <- list()
aa_valid_perf_12 <- list()
aa_cnt <- 1
for(aa_model in 1:length(aa_sp_name)){
  print(aa_model)
  aa_nu <- aa_sp_name[aa_model]
  for(aa_step in 1:50){
    aa_tr_file <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/train_out_", aa_nu, "_", aa_step, ".txt")
    aa_ts_file <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/test_out_", aa_nu, "_", aa_step, ".txt")
    aa_va_file <-  paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/valid_out_", aa_nu, "_", aa_step, ".txt")
    aa_par_file <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Motif_opt_output/updated_par_", aa_nu, "_", aa_step, ".txt") 
    if (file.exists(aa_tr_file)){# & !(aa_step == 3 & aa_model == 28) & !(aa_step == 9 & aa_model == 40)
      aa_Train_perf_12[[aa_cnt]] <- read_output_train_test_GEMSTAT_indiv(output_file = aa_tr_file)
      motif_optim_perf_auroc_train_12[aa_model, aa_step] <- aa_Train_perf_12[[aa_cnt]]$ROC_curve[[2]]
      motif_optim_perf_auprc_train_12[aa_model, aa_step] <- aa_Train_perf_12[[aa_cnt]]$PRC_curve[[2]]
      aapar <- fromJSON(file = aa_par_file)
      aa_bias <- aapar$log_Reg[[1]]$bias
      aa_coeff <- aapar$log_Reg[[1]]$coeff
      aa_beta <- aapar$enh[[1]]$beta
      motif_optim_perf_auprc_train_12_logistloss[aa_model, aa_step] <- read_output_train_test_GEMSTAT_indiv_logistic(output_file = aa_tr_file, 
                                                                                                                     bias_par =aa_bias,
                                                                                                                     coeff_par =aa_coeff,
                                                                                                                     weights = aa_tr_file,
                                                                                                                     #normalize_output_01 = T
                                                                                                                     beta = aa_beta)
    }
    if (file.exists(aa_ts_file)){#& !(aa_step == 3 & aa_model == 28) & !(aa_step == 9 & aa_model == 40)
      aa_Test_perf_12[[aa_cnt]] <- read_output_train_test_GEMSTAT_indiv(output_file = aa_ts_file)
      motif_optim_perf_auroc_test_12[aa_model, aa_step] <- aa_Test_perf_12[[aa_cnt]]$ROC_curve[[2]]
      motif_optim_perf_auprc_test_12[aa_model, aa_step] <- aa_Test_perf_12[[aa_cnt]]$PRC_curve[[2]]
      aapar <- fromJSON(file = aa_par_file)
      aa_bias <- aapar$log_Reg[[1]]$bias
      aa_coeff <- aapar$log_Reg[[1]]$coeff
      motif_optim_perf_auprc_test_12_logistloss[aa_model, aa_step] <- read_output_train_test_GEMSTAT_indiv_logistic(output_file = aa_ts_file, 
                                                                                                                     bias_par = aa_bias,
                                                                                                                     coeff_par = aa_coeff,
                                                                                                                     weights = aa_ts_file,
                                                                                                                    #normalize_output_01 = T,
                                                                                                                    beta = aa_beta)
      
    }
    if (file.exists(aa_va_file)){#& !(aa_step == 3 & aa_model == 28) & !(aa_step == 9 & aa_model == 40)
      aa_valid_perf_12[[aa_cnt]] <- read_output_train_test_GEMSTAT_indiv(output_file =aa_va_file)
      motif_optim_perf_auroc_valid_12[aa_model, aa_step] <- aa_valid_perf_12[[aa_cnt]]$ROC_curve[[2]]
      motif_optim_perf_auprc_valid_12[aa_model, aa_step] <- aa_valid_perf_12[[aa_cnt]]$PRC_curve[[2]]
      aapar <- fromJSON(file = aa_par_file)
      aa_bias <- aapar$log_Reg[[1]]$bias
      aa_coeff <- aapar$log_Reg[[1]]$coeff
      motif_optim_perf_auprc_valid_12_logistloss[aa_model, aa_step] <- read_output_train_test_GEMSTAT_indiv_logistic(output_file = aa_va_file, 
                                                                                                                     bias_par = aa_bias,
                                                                                                                     coeff_par = aa_coeff,
                                                                                                                     weights = aa_va_file,
                                                                                                                     #normalize_output_01 = T,
                                                                                                                     beta = aa_beta)
      
    }
    aa_cnt <- aa_cnt + 1
  }
}
############################################################################################################
# plot the results
par(mfrow = c(3, 3), mar = c(2,2,2,2))
boxplot.matrix(motif_optim_perf_auprc_train_12, main= " AUPRC train")
boxplot.matrix(motif_optim_perf_auroc_train_12, main= " AUROC train")
boxplot.matrix(motif_optim_perf_auprc_test_12, main= " AUPRC test")
boxplot.matrix(motif_optim_perf_auroc_test_12, main= " AUROC test")
boxplot.matrix(motif_optim_perf_auprc_valid_12, main= " AUPRC valid")
boxplot.matrix(motif_optim_perf_auroc_valid_12, main= " AUROC valid")
boxplot.matrix(motif_optim_perf_auprc_train_12_logistloss, main = "logist loss train")
boxplot.matrix(motif_optim_perf_auprc_valid_12_logistloss, main = "logist loss valid")
boxplot.matrix(motif_optim_perf_auprc_test_12_logistloss, main = "logist loss test")

aa_range_prc <- range(cbind(motif_optim_perf_auprc_train_12, motif_optim_perf_auprc_test_12, motif_optim_perf_auprc_valid_12), na.rm = T)
aa_range_roc <- range(cbind(motif_optim_perf_auroc_train_12, motif_optim_perf_auroc_test_12, motif_optim_perf_auroc_valid_12), na.rm = T)
aa_range_logit <- range(cbind(motif_optim_perf_auprc_train_12_logistloss, motif_optim_perf_auprc_test_12_logistloss, motif_optim_perf_auprc_valid_12_logistloss), na.rm = T)
aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas1[1:25]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas2[1:25]]
aaadd3 <- union(aaadd1, aaadd2)

aa_sp_name <-  unlist(lapply(strsplit(aaadd3, split = "_"), "[[", 2))


motif_optim_perf_auprc_compare_12 <- matrix(0L,nrow = length(aa_sp_name), 
                                         ncol = 50 )
motif_optim_perf_auroc_compare_12 <- matrix(0L,nrow = length(aa_sp_name), 
                                         ncol = 50 )
motif_optim_perf_auprc_compare_12_valid <- matrix(0L,nrow = length(aa_sp_name), 
                                            ncol = 50 )
motif_optim_perf_auroc_compare_12_valid <- matrix(0L,nrow = length(aa_sp_name), 
                                                  ncol = 50 )


par(mfrow = c(6, 7), mar = c(2,2,2,2))
for(aa_model in 1:length(aaadd3)){
  plot(motif_optim_perf_auprc_train_12[aa_model,], type = "l", col = 1, ylim = aa_range_prc, main = aa_sp_name[aa_model])
  abline(h = E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_PRC[paste0("output_par_",aa_sp_name[aa_model], "_1")], col = 1, lty = 2)
  
  
  lines(motif_optim_perf_auprc_test_12[aa_model,], col = 2)
  aa_base_test <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_PRC[paste0("output_",aaadd3[aa_model], "_test")]
  abline(h = aa_base_test, col = 2, lty = 2)
  motif_optim_perf_auprc_compare_12[aa_model, ] <- motif_optim_perf_auprc_compare_12[aa_model, ] + as.integer(motif_optim_perf_auprc_test_12[aa_model,] > aa_base_test)
  
  
  aa_base_valid <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_PRC[paste0("output_",aaadd3[aa_model], "_valid")]
  lines(motif_optim_perf_auprc_valid_12[aa_model,], col = 3)
  abline(h = aa_base_valid, col = 3, lty = 2)
  motif_optim_perf_auprc_compare_12[aa_model, ] <- motif_optim_perf_auprc_compare_12[aa_model, ] + as.integer(motif_optim_perf_auprc_valid_12[aa_model,] > aa_base_valid)
  motif_optim_perf_auprc_compare_12_valid[aa_model, ] <- motif_optim_perf_auprc_valid_12[aa_model,] - aa_base_valid
}

which(rowSums(motif_optim_perf_auprc_compare_12_valid>0, na.rm=T) > 0)
rownames(motif_optim_perf_auprc_test_12)[which(rowSums(motif_optim_perf_auprc_compare_12_valid>0, na.rm=T) > 0)]
E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_PRC[paste0("output_par_",7084, "_1_valid")]

par(mfrow = c(6, 7), mar = c(2,2,2,2))
for(aa_model in 1:length(aaadd3)){
  plot(motif_optim_perf_auroc_train_12[aa_model,], type = "l", col = 1, ylim = aa_range_roc, main = aa_sp_name[aa_model])
  abline(h = E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_ROC[paste0("output_par_",aa_sp_name[aa_model], "_1")], col = 1, lty = 2)
  
  lines(motif_optim_perf_auroc_test_12[aa_model,], col = 2)
  aa_base_test <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_ROC[paste0("output_",aaadd3[aa_model], "_test")]
  abline(h = aa_base_test, col = 2, lty = 2)
  motif_optim_perf_auroc_compare_12[aa_model, ] <- motif_optim_perf_auroc_compare_12[aa_model, ] + as.integer(motif_optim_perf_auroc_test_12[aa_model,] > aa_base_test)
  
  
  lines(motif_optim_perf_auroc_valid_12[aa_model,], col = 3)
  aa_base_valid <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_ROC[paste0("output_",aaadd3[aa_model], "_valid")]
  abline(h = aa_base_valid, col = 3, lty = 2)
  motif_optim_perf_auroc_compare_12[aa_model, ] <- motif_optim_perf_auroc_compare_12[aa_model, ] + as.integer(motif_optim_perf_auroc_valid_12[aa_model,] > aa_base_valid)
  motif_optim_perf_auroc_compare_12_valid[aa_model, ] <- motif_optim_perf_auroc_valid_12[aa_model,] - aa_base_valid
  
}



aaa1 <- motif_optim_perf_auroc_compare_12_valid
aaa2 <- motif_optim_perf_auprc_compare_12_valid
aaa1[motif_optim_perf_auroc_compare_12_valid > 0 &
             motif_optim_perf_auprc_compare_12_valid > 0] <- 1000
aaa2[motif_optim_perf_auroc_compare_12_valid > 0 &
             motif_optim_perf_auprc_compare_12_valid > 0] <- 1000
which(rowSums(aaa1 == 1000 , na.rm = T) > 0)
which(rowSums(aaa2 == 1000, na.rm = T) > 0)
# 
rownames(motif_optim_perf_auroc_valid_12)[which(rowSums(aaa2 == 1000 , na.rm = T) > 0)]

motif_optim_perf_auroc_compare_12_valid[which(motif_optim_perf_auroc_compare_12_valid > 0 &
                                                motif_optim_perf_auprc_compare_12_valid > 0)]
motif_optim_perf_auprc_compare_12_valid[which(motif_optim_perf_auroc_compare_12_valid > 0 &
                                                motif_optim_perf_auprc_compare_12_valid > 0)]



aaa1 <- motif_optim_perf_auprc_compare_12
which(rowSums(aaa1 == 2, na.rm = T) > 0)
which(motif_optim_perf_auprc_compare_12[8, ] == 2)
rownames(motif_optim_perf_auprc_valid_12)[which(rowSums(aaa1 == 2, na.rm = T) > 0)]
which(motif_optim_perf_auprc_compare_12[31, ] == 2)

aaa2 <- motif_optim_perf_auroc_compare_12
which(rowSums(aaa2 == 2, na.rm = T) > 0)
which(motif_optim_perf_auroc_compare_12[1, ] == 2)
rownames(motif_optim_perf_auroc_valid_12)[which(rowSums(aaa2 == 2, na.rm = T) > 0)]
which(motif_optim_perf_auroc_compare_12[31, ] == 2)


par(mfrow = c(6, 7), mar = c(2,2,2,2))
for(aa_model in 1:length(aaadd3)){
  plot(motif_optim_perf_auprc_train_12_logistloss[aa_model,], type = "l", col = 1, ylim = aa_range_logit, main = aa_sp_name[aa_model])
#  abline(h = E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$train_PRC[paste0("output_par_",aa_sp_name[aa_model], "_1")], col = 1, lty = 2)
  
  
  lines(motif_optim_perf_auprc_test_12_logistloss[aa_model,], col = 2)
#  aa_base_test <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$test_PRC[paste0("output_",aaadd3[aa_model], "_test")]
#  abline(h = aa_base_test, col = 2, lty = 2)
#  motif_optim_perf_auprc_compare_12[aa_model, ] <- motif_optim_perf_auprc_compare_12[aa_model, ] + as.integer(motif_optim_perf_auprc_test_12[aa_model,] > aa_base_test)
  
  
#  aa_base_valid <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_PRC[paste0("output_",aaadd3[aa_model], "_valid")]
  lines(motif_optim_perf_auprc_valid_12_logistloss[aa_model,], col = 3)
#  abline(h = aa_base_valid, col = 3, lty = 2)
#  motif_optim_perf_auprc_compare_12[aa_model, ] <- motif_optim_perf_auprc_compare_12[aa_model, ] + as.integer(motif_optim_perf_auprc_valid_12[aa_model,] > aa_base_valid)
#  motif_optim_perf_auprc_compare_12_valid[aa_model, ] <- motif_optim_perf_auprc_valid_12[aa_model,] - aa_base_valid
}

############################################################################################################
# write Jpegs of motif logos

aa_cur_model <- 5057
aa <- list()
for(i in 1:20){
  aapngfname <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Optimized_motifs_100_logo/logo_",aa_cur_model,"_", i, ".png")
  aacf <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Motif_opt_grad/base_motif_", aa_cur_model,"_", i, ".RData" )
  if(file.exists(aacf)){
    load(aacf)
    png(filename = aapngfname, width = 600, height = 300, units = "px" )
    aa[[i]] <- original_motif_list[[1]]
    seqLogo::seqLogo(pwm = t(original_motif_list[[1]]))
    dev.off()
  }else{
    print(paste0("model ", aa_cur_model, " step ", i, " not present."))
  }
}
max(abs(TF.motifs.Expanded_pseudo[[aa_names[1]]] - aa[[1]]))
max(abs(aa[[2]] - aa[[1]]))
max(abs(aa[[3]] - aa[[2]]))
max(abs(TF.motifs.Expanded_pseudo[[aa_names[1]]] - aa[[1]]))
############################################################################################################
# run the same thing but with annotation
aa_call_noannot <-  "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1"
aa_call_annot <-  "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1 -a Annotation/train.ann "

# -a Annotation/train.ann 

aa_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Motif_opt_input", pattern = "seq*", full.names = T)
aa_files_s  <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Motif_opt_input", pattern = "seq*", full.names = F)
aa_files_s <- unlist(lapply(strsplit(aa_files_s, split = "\\."), "[[", 1))
aa_files_s <- paste0(aa_files_s, ".annot")
system("mkdir Annotation_steps")
cat(c("#!/bin/bash\n"), file = "annotation_steps.sh", append = F)
for(i in 1:length(aa_files)){
  cat(c("/Users/Shayan/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seqannot", 
        "-s", aa_files[i],
        "-m /Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/motifs.wtmx",
        ">>",  paste0("/Users/Shayan/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Annotation_steps/", aa_files_s[i], "\n")),
      sep = " ", append = T, file = "annotation_steps.sh")
}



setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Jobs_with_annot/")

aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas1[1:25]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas2[1:25]]
aaadd3 <- union(aaadd1, aaadd2)

aa_sp_name <-  unlist(lapply(strsplit(aaadd3, split = "_"), "[[", 2))

for(aa_model in 1:length(aaadd3)){
  aa_adres <- paste0("Trained_par/", aaadd3[aa_model], ".txt.Filter")
  aa_nu <- aa_sp_name[aa_model]
  aa_dag <- paste0("job_", aa_nu)
  aa_dagf <- paste0("my_DagMan_", aa_nu)
  
  motif_optimizer_job_writer(nu_steps = 20,
                             trained_par_address= aa_adres,
                             model_number_name = aa_nu,
                             motif_length = nrow(TF.motifs.Expanded_pseudo$ESR1_2),
                             grad_move_step = 1,
                             grad_calc_step = 0.01,
                             GEMSTAT_call = aa_call_noannot,
                             annotation_fixed = T,
                             motif_index = 1,
                             main_directory ="/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/")
  motif_optimizer_job_2_submit(model_name_number = aa_nu,
                               motif_length = nrow(TF.motifs.Expanded_pseudo$ESR1_2),
                               nu_steps = 20)
  DAGmanConstructor(submit_prefix = aa_dag,
                    start_ind=1,
                    end_ind=80,
                    my_order = c(1:80),
                    filename=aa_dagf)
}

# create script to make these exec and run dags
cat(c("#!/bin/bash\n"), file = "job_2_sub_executable_dagrun.sh", append = F)
for(aa_model in 1:length(aaadd3)){
  aa_nu <- aa_sp_name[aa_model]
  cat(c("chmod", "+x", paste0("job_to_submit_motif_optimizer_",aa_nu, ".sh\n" )),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("./job_to_submit_motif_optimizer_",aa_nu, ".sh\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c("chmod", "+x", paste0("submit_to_Exec_",aa_nu, ".sh\n" )),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("./submit_to_Exec_",aa_nu, ".sh\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("condor_submit_dag my_DagMan_",aa_nu, "\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
}
###################
# write evaluation jobs with annotation
aa_sp_name <-  unlist(lapply(strsplit(aaadd3, split = "_"), "[[", 2))

# write motifs from RData

for(aa_model in 1:length(aaadd3)){
  aa_nu <- aa_sp_name[aa_model]
  for(aa_step in 1:20){
    aa_read_adres  <- paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_grad/base_motif_", 
                             aa_nu, "_", aa_step, ".RData")
    aa_write_adres <- paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Optimized_motifs/motif_", 
                             aa_nu, "_", aa_step)
    # job to write motifs
    cat(c("Rscript --vanilla motif_writer_script.R", aa_read_adres, paste0(aa_write_adres, "\n")), 
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)), file = "write_motifs.job")
    # job to run GEMSTAT on training set
    aa_GEM_call <- "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1"
    aa_train_annot <- paste0("-a /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Annotation_steps/seq_GEMS_step_", aa_step, ".annot")
    aa_train_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_input/seq_GEMS_step_", aa_step, ".fa")
    aa_train_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_input/Label_GEMS_step_", aa_step, ".tab")
    
    aa_test_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Sequence/test_seq.fa")
    aa_test_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Labels/Label_test.txt")
    aa_test_annot <- paste0("-a /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Annotation/test.ann")
    
    aa_valid_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Sequence/valid_seq.fa")
    aa_valid_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Labels/Label_validation.txt")
    aa_valid_annot <- paste0("-a /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Annotation/valid.ann")
    
    
    aa_cur_params <- paste0("-p /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_output/updated_par_",aa_nu,"_",aa_step,".txt")
    aa_cur_motifs <- paste0("-m ", aa_write_adres, ".wtmx")
    aa_na <- "-na 0"
    aa_cur_out_train <- paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/train_out_", aa_nu, "_", aa_step, ".txt")
    aa_cur_out_test <-  paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/test_out_", aa_nu, "_", aa_step, ".txt")
    aa_cur_out_valid <- paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/valid_out_", aa_nu, "_", aa_step, ".txt")
    
    # job to run GEMSTAT on test set
    # job to run GEMSTAT on valid set
    cat(c(aa_GEM_call, aa_train_seq, aa_train_lab, aa_train_annot, aa_cur_params,
          aa_cur_motifs, aa_na, paste0(aa_cur_out_train, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_train_eval.job")
    cat(c(aa_GEM_call, aa_test_seq, aa_test_lab, aa_test_annot, aa_cur_params,
          aa_cur_motifs, aa_na,  paste0(aa_cur_out_test, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_test_eval.job")
    cat(c(aa_GEM_call, aa_valid_seq, aa_valid_lab, aa_valid_annot, aa_cur_params,
          aa_cur_motifs, aa_na,  paste0(aa_cur_out_valid, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_valid_eval.job")
  }
}
# write job2submit scripts

cat("#!/bin/bash\n", file = "write_motifs_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "write_motifs.job",
      "840",
      paste0("tmp_", "write_motifs"),
      ">", 
      paste0("write_motifs.submit", "\n")),
    sep = " ", file = "write_motifs_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_train_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_train_eval.job",
      "840",
      paste0("tmp_", "GEM_train_eval"),
      ">", 
      paste0("GEM_train_eval.submit", "\n")),
    sep = " ", file = "GEM_train_eval_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_test_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_test_eval.job",
      "840",
      paste0("tmp_", "GEM_test_eval"),
      ">", 
      paste0("GEM_test_eval.submit", "\n")),
    sep = " ", file = "GEM_test_eval_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_valid_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_valid_eval.job",
      "840",
      paste0("tmp_", "GEM_valid_eval"),
      ">", 
      paste0("GEM_valid_eval.submit", "\n")),
    sep = " ", file = "GEM_valid_eval_job2submit.sh", append = T)

############################################################################################################
#visulaze results in terms of gemstat obj func

############################################################################################################
# run a toy example to find out what is wrong with motif optimization
# print gradient matrices
# fix annotations


############################
# running without stoch grad dec
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[4]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12")

motif_optimizer_seq_label_weight_creator(nu_steps = 20,
                                         batch_size_motif = 1272,
                                         batch_size_GEMSTAT = 1272,
                                         training_seq_pool = GEMSTAT_Ensemble_train_SeqList[[4]],
                                         training_seq_label = aalab1)




cat(c("#!/bin/bash\n"), file = "job_2_sub_executable_dagrun_rerun.sh", append = F)
for(aa_model in 1:length(aaadd3)){
  aa_nu <- aa_sp_name[aa_model]
  # cat(c("chmod", "+x", paste0("job_to_submit_motif_optimizer_",aa_nu, ".sh\n" )),
  #     sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  # cat(c(paste0("./job_to_submit_motif_optimizer_",aa_nu, ".sh\n")),
  #     sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  # cat(c("chmod", "+x", paste0("submit_to_Exec_",aa_nu, ".sh\n" )),
  #     sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  # cat(c(paste0("./submit_to_Exec_",aa_nu, ".sh\n")),
  #     sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("condor_submit_dag -force my_DagMan_",aa_nu, "\n")),
      sep = " ", file = "job_2_sub_executable_dagrun_rerun.sh", append = T)
}



######
norm_vec <- function(x) sqrt(sum(x^2))
norm_vec(aaa)
aaa <- matrix(runif(20), nrow = 4)
c(aaa)


aaa <- original_motif_list[[1]]
norm(aaa, type= "O")
norm_vec(aaa)
#########################################################################################################
# problem with training is now solved


# run a 50 step motif optimization of the full data set with annotation
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[4]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12")

motif_optimizer_seq_label_weight_creator(nu_steps = 50,
                                         batch_size_motif = 1272,
                                         batch_size_GEMSTAT = 1272,
                                         training_seq_pool = GEMSTAT_Ensemble_train_SeqList[[4]],
                                         training_seq_label = aalab1)



aa_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Motif_opt_input", pattern = "seq*", full.names = T)
aa_files_s  <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Motif_opt_input", pattern = "seq*", full.names = F)
aa_files_gem <- paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_input/", aa_files_s)
aa_files_sann <- unlist(lapply(strsplit(aa_files_s, split = "\\."), "[[", 1))
aa_files_sann <- paste0(aa_files_sann, ".annot")
system("mkdir Annotation_steps")
#cat(c("#!/bin/bash\n"), file = "annotation_steps.job", append = F)
for(i in 1:length(aa_files_gem)){
  cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seqannot", 
        "-s", aa_files_gem[i],
        "-m /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/motifs.wtmx",
        ">>",  paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Annotation_steps/", aa_files_sann[i], "\n")),
      sep = " ", append = (i != 1), file = "annotation_steps.job")
}

cat("#!/bin/bash\n", file = "annotation_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "annotation_steps.job",
      "50",
      paste0("tmp_", "annotation_steps"),
      ">", 
      paste0("annotation_steps.submit", "\n")),
    sep = " ", file = "annotation_job2submit.sh", append = T)


setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Jobs_with_annot/")

aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas1[1:25]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas2[1:25]]
aaadd3 <- union(aaadd1, aaadd2)

aa_sp_name <-  unlist(lapply(strsplit(aaadd3, split = "_"), "[[", 2))

for(aa_model in 1:length(aaadd3)){
  aa_adres <- paste0("Trained_par/", aaadd3[aa_model], ".txt.Filter")
  aa_nu <- aa_sp_name[aa_model]
  aa_dag <- paste0("job_", aa_nu)
  aa_dagf <- paste0("my_DagMan_", aa_nu)
  
  motif_optimizer_job_writer(nu_steps = 50,
                             trained_par_address= aa_adres,
                             model_number_name = aa_nu,
                             motif_length = nrow(TF.motifs.Expanded_pseudo$ESR1_2),
                             grad_move_step = 1,
                             grad_calc_step = 0.01,
                             GEMSTAT_call = aa_call_noannot,
                             annotation_fixed = T,
                             motif_index = 1,
                             main_directory ="/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/")
  motif_optimizer_job_2_submit(model_name_number = aa_nu,
                               motif_length = nrow(TF.motifs.Expanded_pseudo$ESR1_2),
                               nu_steps = 50)
  DAGmanConstructor(submit_prefix = aa_dag,
                    start_ind=1,
                    end_ind=200,
                    my_order = c(1:200),
                    filename=aa_dagf)
}

# create script to make these exec and run dags
cat(c("#!/bin/bash\n"), file = "job_2_sub_executable_dagrun.sh", append = F)
for(aa_model in 1:length(aaadd3)){
  aa_nu <- aa_sp_name[aa_model]
  cat(c("chmod", "+x", paste0("job_to_submit_motif_optimizer_",aa_nu, ".sh\n" )),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("./job_to_submit_motif_optimizer_",aa_nu, ".sh\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c("chmod", "+x", paste0("submit_to_Exec_",aa_nu, ".sh\n" )),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("./submit_to_Exec_",aa_nu, ".sh\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("condor_submit_dag my_DagMan_",aa_nu, "\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
}

aa_sp_name <-  unlist(lapply(strsplit(aaadd3, split = "_"), "[[", 2))

# write motifs from RData

for(aa_model in 1:length(aaadd3)){
  aa_nu <- aa_sp_name[aa_model]
  for(aa_step in 1:50){
    aa_read_adres  <- paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_grad/base_motif_", 
                             aa_nu, "_", aa_step, ".RData")
    aa_write_adres <- paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Optimized_motifs/motif_", 
                             aa_nu, "_", aa_step)
    # job to write motifs
    cat(c("Rscript --vanilla motif_writer_script.R", aa_read_adres, paste0(aa_write_adres, "\n")), 
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)), file = "write_motifs.job")
    # job to run GEMSTAT on training set
    aa_GEM_call <- "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1"
    aa_train_annot <- paste0("-a /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Annotation_steps/seq_GEMS_step_", aa_step, ".annot")
    aa_train_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_input/seq_GEMS_step_", aa_step, ".fa")
    aa_train_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_input/Label_GEMS_step_", aa_step, ".tab")
    
    aa_test_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Sequence/test_seq.fa")
    aa_test_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Labels/Label_test.txt")
    aa_test_annot <- paste0("-a /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Annotation/test.ann")
    
    aa_valid_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Sequence/valid_seq.fa")
    aa_valid_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Labels/Label_validation.txt")
    aa_valid_annot <- paste0("-a /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Annotation/valid.ann")
    
    
    aa_cur_params <- paste0("-p /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_output/updated_par_",aa_nu,"_",aa_step,".txt")
    aa_cur_motifs <- paste0("-m ", aa_write_adres, ".wtmx")
    aa_na <- "-na 0"
    aa_cur_out_train <- paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/train_out_", aa_nu, "_", aa_step, ".txt")
    aa_cur_out_test <-  paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/test_out_", aa_nu, "_", aa_step, ".txt")
    aa_cur_out_valid <- paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/valid_out_", aa_nu, "_", aa_step, ".txt")
    
    # job to run GEMSTAT on test set
    # job to run GEMSTAT on valid set
    cat(c(aa_GEM_call, aa_train_seq, aa_train_lab, aa_train_annot, aa_cur_params,
          aa_cur_motifs, aa_na, paste0(aa_cur_out_train, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_train_eval.job")
    cat(c(aa_GEM_call, aa_test_seq, aa_test_lab, aa_test_annot, aa_cur_params,
          aa_cur_motifs, aa_na,  paste0(aa_cur_out_test, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_test_eval.job")
    cat(c(aa_GEM_call, aa_valid_seq, aa_valid_lab, aa_valid_annot, aa_cur_params,
          aa_cur_motifs, aa_na,  paste0(aa_cur_out_valid, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_valid_eval.job")
  }
}
# write job2submit scripts

cat("#!/bin/bash\n", file = "write_motifs_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "write_motifs.job",
      "2100",
      paste0("tmp_", "write_motifs"),
      ">", 
      paste0("write_motifs.submit", "\n")),
    sep = " ", file = "write_motifs_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_train_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_train_eval.job",
      "2100",
      paste0("tmp_", "GEM_train_eval"),
      ">", 
      paste0("GEM_train_eval.submit", "\n")),
    sep = " ", file = "GEM_train_eval_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_test_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_test_eval.job",
      "2100",
      paste0("tmp_", "GEM_test_eval"),
      ">", 
      paste0("GEM_test_eval.submit", "\n")),
    sep = " ", file = "GEM_test_eval_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_valid_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_valid_eval.job",
      "2100",
      paste0("tmp_", "GEM_valid_eval"),
      ">", 
      paste0("GEM_valid_eval.submit", "\n")),
    sep = " ", file = "GEM_valid_eval_job2submit.sh", append = T)







##############################################################################################################


# run a 50 step motif optimization of the full data set without annotation


setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Jobs_with_no_annot/")

aasum1 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_PRC
aas1 <- sort(aasum1, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas1[1:25]]
aasum2 <- E_RNA_GEMSTAT_Ensemble_Outlist[[12]]$valid_ROC
aas2 <- sort(aasum2, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(E_RNA_GEMSTAT_Ensemble_Parlist[[12]]$annot)[aas2[1:25]]
aaadd3 <- union(aaadd1, aaadd2)

aa_sp_name <-  unlist(lapply(strsplit(aaadd3, split = "_"), "[[", 2))

for(aa_model in 1:length(aaadd3)){
  aa_adres <- paste0("Trained_par/", aaadd3[aa_model], ".txt.Filter")
  aa_nu <- aa_sp_name[aa_model]
  aa_dag <- paste0("job_", aa_nu)
  aa_dagf <- paste0("my_DagMan_", aa_nu)
  
  motif_optimizer_job_writer(nu_steps = 50,
                             trained_par_address= aa_adres,
                             model_number_name = aa_nu,
                             motif_length = nrow(TF.motifs.Expanded_pseudo$ESR1_2),
                             grad_move_step = 1,
                             grad_calc_step = 0.01,
                             GEMSTAT_call = aa_call_noannot,
                             annotation_fixed = F,
                             motif_index = 1,
                             main_directory ="/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/")
  motif_optimizer_job_2_submit(model_name_number = aa_nu,
                               motif_length = nrow(TF.motifs.Expanded_pseudo$ESR1_2),
                               nu_steps = 50)
  DAGmanConstructor(submit_prefix = aa_dag,
                    start_ind=1,
                    end_ind=200,
                    my_order = c(1:200),
                    filename=aa_dagf)
}

# create script to make these exec and run dags
cat(c("#!/bin/bash\n"), file = "job_2_sub_executable_dagrun.sh", append = F)
for(aa_model in 1:length(aaadd3)){
  aa_nu <- aa_sp_name[aa_model]
  cat(c("chmod", "+x", paste0("job_to_submit_motif_optimizer_",aa_nu, ".sh\n" )),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("./job_to_submit_motif_optimizer_",aa_nu, ".sh\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c("chmod", "+x", paste0("submit_to_Exec_",aa_nu, ".sh\n" )),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("./submit_to_Exec_",aa_nu, ".sh\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
  cat(c(paste0("condor_submit_dag my_DagMan_",aa_nu, "\n")),
      sep = " ", file = "job_2_sub_executable_dagrun.sh", append = T)
}

aa_sp_name <-  unlist(lapply(strsplit(aaadd3, split = "_"), "[[", 2))


# write motifs from RData

for(aa_model in 1:length(aaadd3)){
  aa_nu <- aa_sp_name[aa_model]
  for(aa_step in 1:50){
    aa_read_adres  <- paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_grad/base_motif_", 
                             aa_nu, "_", aa_step, ".RData")
    aa_write_adres <- paste0("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Optimized_motifs/motif_", 
                             aa_nu, "_", aa_step)
    # job to write motifs
    cat(c("Rscript --vanilla motif_writer_script.R", aa_read_adres, paste0(aa_write_adres, "\n")), 
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)), file = "write_motifs.job")
    # job to run GEMSTAT on training set
    aa_GEM_call <- "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -f TF_exp.tab  -o DIRECT -c Coop/coop.par -lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -onebeta 1"
 #   aa_train_annot <- paste0("-a /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Annotation_steps/seq_GEMS_step_", aa_step, ".annot")
    aa_train_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_input/seq_GEMS_step_", aa_step, ".fa")
    aa_train_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_input/Label_GEMS_step_", aa_step, ".tab")
    
    aa_test_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Sequence/test_seq.fa")
    aa_test_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Labels/Label_test.txt")
  #  aa_test_annot <- paste0("-a /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Annotation/test.ann")
    
    aa_valid_seq <- paste0("-s /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Sequence/valid_seq.fa")
    aa_valid_lab <- paste0("-e /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Labels/Label_validation.txt")
  #  aa_valid_annot <- paste0("-a /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Annotation/valid.ann")
    
    
    aa_cur_params <- paste0("-p /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Motif_opt_output/updated_par_",aa_nu,"_",aa_step,".txt")
    aa_cur_motifs <- paste0("-m ", aa_write_adres, ".wtmx")
    aa_na <- "-na 0"
    aa_cur_out_train <- paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/train_out_", aa_nu, "_", aa_step, ".txt")
    aa_cur_out_test <-  paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/test_out_", aa_nu, "_", aa_step, ".txt")
    aa_cur_out_valid <- paste0("-fo /shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_12/Model_eval_tr_ts_val/valid_out_", aa_nu, "_", aa_step, ".txt")
    
    # job to run GEMSTAT on test set
    # job to run GEMSTAT on valid set
    cat(c(aa_GEM_call, 
          aa_train_seq, 
          aa_train_lab,
         # aa_train_annot,
          aa_cur_params,
          aa_cur_motifs,
          aa_na, 
          paste0(aa_cur_out_train, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_train_eval.job")
    cat(c(aa_GEM_call, 
          aa_test_seq, 
          aa_test_lab,
          #aa_test_annot,
          aa_cur_params,
          aa_cur_motifs, 
          aa_na, 
          paste0(aa_cur_out_test, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_test_eval.job")
    cat(c(aa_GEM_call, 
          aa_valid_seq,
          aa_valid_lab, 
         # aa_valid_annot,
          aa_cur_params,
          aa_cur_motifs, 
          aa_na, 
          paste0(aa_cur_out_valid, "\n")),
        sep = " ", append = (! (aa_model == 1 & aa_step == 1)) , file = "GEM_valid_eval.job")
  }
}
# write job2submit scripts

cat("#!/bin/bash\n", file = "write_motifs_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "write_motifs.job",
      "2100",
      paste0("tmp_", "write_motifs"),
      ">", 
      paste0("write_motifs.submit", "\n")),
    sep = " ", file = "write_motifs_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_train_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_train_eval.job",
      "2100",
      paste0("tmp_", "GEM_train_eval"),
      ">", 
      paste0("GEM_train_eval.submit", "\n")),
    sep = " ", file = "GEM_train_eval_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_test_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_test_eval.job",
      "2100",
      paste0("tmp_", "GEM_test_eval"),
      ">", 
      paste0("GEM_test_eval.submit", "\n")),
    sep = " ", file = "GEM_test_eval_job2submit.sh", append = T)

cat("#!/bin/bash\n", file = "GEM_valid_eval_job2submit.sh", sep = "", append = F)
cat(c("/shared-mounts/sinhas/tabebor2/create_jobs_submit.pl",
      "GEM_valid_eval.job",
      "2100",
      paste0("tmp_", "GEM_valid_eval"),
      ">", 
      paste0("GEM_valid_eval.submit", "\n")),
    sep = " ", file = "GEM_valid_eval_job2submit.sh", append = T)



##############################################################################################################
# Chose model with maximum validation AUPRC: 5057 step 1



aa_cur_model <- 5057
aa <- list()
for(i in 1:1){
  aapngfname <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Optimized_motifs_100_logo/logo_",aa_cur_model,"_", i, ".png")
  aacf <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_12/Motif_opt_grad/base_motif_", aa_cur_model,"_", i, ".RData" )
  if(file.exists(aacf)){
    load(aacf)
#    png(filename = aapngfname, width = 600, height = 300, units = "px" )
    aa[[i]] <- original_motif_list[[1]]
#    seqLogo::seqLogo(pwm = t(original_motif_list[[1]]))
#    dev.off()
  }else{
#    print(paste0("model ", aa_cur_model, " step ", i, " not present."))
  }
}
max(abs(TF.motifs.Expanded_pseudo[[aa_names[1]]] - aa[[1]]))
max(abs(aa[[2]] - aa[[1]]))
max(abs(aa[[3]] - aa[[2]]))
max(abs(TF.motifs.Expanded_pseudo[[aa_names[1]]] - aa[[1]]))


TF.motifs.Expanded_pseudo_optim <- TF.motifs.Expanded_pseudo
TF.motifs.Expanded_pseudo_optim[["ESR1_2"]] <- aa[[1]]
##############################################################################################################






