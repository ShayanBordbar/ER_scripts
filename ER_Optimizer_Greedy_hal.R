#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(sigmoid)
# step_counter has been set to one in the initial workplace and is being incremented by one at every step
# nu_jobs has been saved in the initial workplace
# nu_genes is the number of genes, should be set in initial workplace

if(file.exists(paste0("Optim_Greedy_", args[1], ".RData"))){
	load(paste0("Optim_Greedy_", args[1], ".RData"))
	my_initial_chosen_enh <- GreedyRes[[step_counter-1]]$Enhancer_index
	my_initial_par <- GreedyRes[[step_counter-1]]$parameters
}else{
	load("Inputs_greedy.RData")
	my_initial_chosen_enh <- integer(0)
	my_initial_par <- numeric(0)

}
set.seed((as.integer(args[1]) + (step_counter - 1)*nu_jobs))
GreedyRes <- list()
GreedyRes[[step_counter]] <- GreedyObjOptim(.my_data = ER_52_opt_simAnn_input_ChipFilt_1500_500,
                                             initial_chosen_enh=my_initial_chosen_enh,
                                             initial_par=my_initial_par,
                                             gene_number=as.integer(args[2])) 
step_counter <- step_counter + 1

save.image(paste0("Optim_Greedy_", args[1], ".RData"))