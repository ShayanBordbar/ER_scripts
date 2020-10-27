# !/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
source("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/GEMSTAT_parameter_opt_functions.R")
my_par <- args[1]
my_lb <- args[2]
my_ub <- args[3]


GEMSTAT_param_bound_checker(par_file =my_par,
                            ub_file = my_ub,
                            lb_file = my_lb)
