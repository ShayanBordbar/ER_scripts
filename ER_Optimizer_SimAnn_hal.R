#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(sigmoid)
# step_counter has been set to one in the initial workplace and is being incremented by one at every step
# nu_jobs has been saved in the initial workplace

if(file.exists(paste0("Optim_sim_Ann_", args[1], ".RData"))){
	load(paste0("Optim_sim_Ann_", args[1], ".RData"))
	last_cnt <- length(SimAnnRes[[step_counter-1]]$Enhancer_index)
	my_initial_chosen_enh <- SimAnnRes[[step_counter-1]]$Enhancer_index[[last_cnt]]
	my_initial_par <- SimAnnRes[[step_counter-1]]$parameters[[last_cnt]]
}else{
	load("Inputs_simAnn.RData")
	my_initial_chosen_enh <- integer(0)
	my_initial_par <- numeric(0)

}
set.seed((as.integer(args[1]) + (step_counter - 1)*nu_jobs))
SimAnnRes <- list()
SimAnnRes[[step_counter]] <- SimAnnObjOptim(.my_data = ER_52_opt_simAnn_input,
                                                   initial_chosen_enh=my_initial_chosen_enh,
                                                   initial_par=my_initial_par,
                                                   starting_temp=as.numeric(args[2]),
                                                   min_temp=as.numeric(args[3]),
                                                   iter_per_temp=1,
                                                   ALPHA=0.95) 
step_counter <- step_counter + 1

save.image(paste0("Optim_sim_Ann_", args[1], ".RData"))