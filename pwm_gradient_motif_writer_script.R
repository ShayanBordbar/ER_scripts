# !/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
source("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/GEMSTAT_parameter_opt_functions.R")
my_my_motif_list_RData_address <- args[1]
my_motif_index <- as.integer(args[2])
my_grad_calc_step <- as.numeric(args[3])
my_model_number_name <- args[4]
my_step_nu <- as.integer(args[5])






pwm_gradient_motif_writer(my_motif_list_RData_address = my_my_motif_list_RData_address,
                          motif_index = my_motif_index, 
                          grad_calc_step = my_grad_calc_step,
                          constrained_rows = c(rep(T, 6), rep(F, 3), rep(T, 6), rep(F, 3)),
                          model_number_name = my_model_number_name,
                          step_nu = my_step_nu,
                          main_directory = "/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/")


