#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# arg[1] spcifies the the row of Optim_Greedy_conc_modif_vec to be read. TF modification index
# arg[2] specifies the gene number 

library(sigmoid)
# step_counter has been set to one in the initial workplace and is being incremented by one at every step
# nu_jobs has been saved in the initial workplace
# nu_genes is the number of genes, should be set in initial workplace

if(file.exists(paste0("Optim_Greedy_conc_modif", args[1], ".RData"))){
	load(paste0("Optim_Greedy_conc_modif", args[1], ".RData"))
	my_initial_chosen_enh <- GreedyRes_conc[[step_counter-1]]$Enhancer_index
	my_initial_par <- GreedyRes_conc[[step_counter-1]]$parameters
}else{
	load("Inputs_greedy_concModif.RData")
	my_initial_chosen_enh <- Sim_Ann_148_restart_enhancers_mostVoted
	my_initial_par <- Sim_Ann_148_restart_parameters_median_top75
	GreedyRes_conc <- list()

}
set.seed((as.integer(args[1]) + (step_counter - 1)*nu_jobs))

GreedyRes_conc[[step_counter]] <- SilicoConcenOptim(..my_data = ER_52_opt_simAnn_input,
                                                    init_enh=my_initial_chosen_enh,
                                                    model_param=my_initial_par,
                                                    .gene_number=as.integer(args[2]),
                                                    Modif_vec=Optim_Greedy_conc_modif_vec[as.integer(args[1]), ])
step_counter <- step_counter + 1

save.image(paste0("Optim_Greedy_conc_modif", args[1], ".RData"))
