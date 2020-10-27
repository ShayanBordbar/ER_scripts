# !/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
source("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/E_RNA_GEMSTAT_Logistic/Ensemble/Experiment_8/GEMSTAT_parameter_opt_functions.R")

my_motif_RData <- args[1] 
my_motif_address <- args[2] 
load(my_motif_RData) # contains original_motif_list

MotifWriter(motif.List = original_motif_list,
            pseudo = 0.001,
            output.File.Name = my_motif_address)